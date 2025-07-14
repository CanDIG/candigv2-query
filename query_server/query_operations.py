from flask import Flask
import config
import copy
import re
import requests
import connexion
import secrets
import urllib
from flask import request, Response
from authx.auth import get_user_id, get_auth_token, is_user_candig_authorized, verify_service_token
from candigv2_logging.logging import CanDIGLogger


logger = CanDIGLogger(__file__)


PAGE_SIZE = 10000000

app = Flask(__name__)
app.config['SECRET_KEY'] = secrets.token_bytes(32)

def get_service_info():
    return {
        "id": "org.candig.query",
        "name": "CanDIG query service",
        "type": {
            "group": "org.candig",
            "artifact": "query",
            "version": "v0.1.0"
        },
        "description": "A query microservice for operating with HTSGet & Katsu",
        "organization": {
            "name": "CanDIG",
            "url": "https://www.distributedgenomics.ca"
        },
        "version": "0.1.0"
    }


def safe_get_response_json(response, name):
    if not response.ok:
        raise Exception(f"Could not get {name} response: {response.status_code} {response.text}")
    return response.json()

def get_headers():
    # Add a service token to the headers so that other services will know this is from the query service:
    headers = {}
    for k in connexion.request.headers.keys():
        headers[k] = connexion.request.headers[k]
    headers["X-Service-Token"] = config.SERVICE_TOKEN
    return headers


# Grab a list of donors matching a given filter from the given URL
def get_donors_from_katsu(url, param_name, parameter_list, headers, therapy_type=None, keep_all=False):
    permissible_donors = set()
    all_results = []
    for parameter in parameter_list:
        # TODO: Fix the page_size call here -- use a consume_all() query like in the frontend
        parameters = {
            param_name: parameter,
            'page_size': PAGE_SIZE
        }
        if therapy_type != None:
            parameters['systemic_therapy_type'] = therapy_type
        treatments = requests.get(f"{url}?{urllib.parse.urlencode(parameters)}", headers=headers)
        results = safe_get_response_json(treatments, f'Katsu {param_name}')['items']
        permissible_donors |= set([result['submitter_donor_id'] for result in results])
        if keep_all:
            all_results.extend(results)

    # If we are required to return all results, query at least once
    if not parameter_list and keep_all:
        parameters = {
            'page_size': PAGE_SIZE
        }
        if therapy_type != None:
            parameters['systemic_therapy_type'] = therapy_type
        treatments = requests.get(f"{url}?{urllib.parse.urlencode(parameters)}", headers=headers)
        results = safe_get_response_json(treatments, f'Katsu {param_name}')['items']
        permissible_donors |= set([result['submitter_donor_id'] for result in results])
        all_results.extend(results)
    return permissible_donors, all_results

def add_or_increment(dict, key):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1

def get_summary_stats(donors, primary_sites, treatments):
    # Perform (and cache) summary statistics
    age_at_diagnosis = {}
    donors_by_id = {}
    primary_site_count = {}
    patients_per_program = {}
    treatment_type_count = {}
    for donor in donors:
        # A donor's date of birth is defined as the (negative) interval between actual DOB and the date of first diagnosis
        # So we just use that info
        donors_by_id[donor["submitter_donor_id"]] = donor
        if donor['date_of_birth'] and donor['date_of_birth']['month_interval']:
            age = abs(donor['date_of_birth']['month_interval']) // 12
            age = age // 10 * 10
            if age < 20:
                add_or_increment(age_at_diagnosis, '0-19 Years')
            elif age > 79:
                add_or_increment(age_at_diagnosis, '80+ Years')
            else:
                add_or_increment(age_at_diagnosis, f'{age}-{age+9} Years')

        program_id = donor['program_id']
        add_or_increment(patients_per_program, program_id)

        # primary sites
        if donor['submitter_donor_id'] in primary_sites:
            if primary_sites[donor['submitter_donor_id']] is not None:
                for primary_site in primary_sites[donor['submitter_donor_id']]:
                    add_or_increment(primary_site_count, str(primary_site))
            else:
                add_or_increment(primary_site_count, 'None')

        if donor['submitter_donor_id'] in treatments and treatments[donor['submitter_donor_id']] is not None:
            for treatment_type in treatments[donor['submitter_donor_id']]:
                add_or_increment(treatment_type_count, treatment_type)

    return {
        'age_at_diagnosis': age_at_diagnosis,
        'treatment_type_count': treatment_type_count,
        'primary_site_count': primary_site_count,
        'patients_per_program': patients_per_program
    }

def query_htsget_gene(headers, gene_array):
    for g in gene_array:
        payload = {
            'query': {
                'requestParameters': {
                    'gene_id': g
                }
            },
            'meta': {
                'apiVersion': 'v2'
            }
        }

        return safe_get_response_json(requests.post(
            f"{config.HTSGET_URL}/beacon/v2/g_variants",
            headers=headers,
            json=payload), 'HTSGet Gene')

def query_htsget_pos(headers, assembly, chrom, start=0, end=10000000):
    payload = {
        'query': {
            'requestParameters': {
                'assemblyId': assembly,
                'referenceName': chrom,
                'start': [start],
                'end': [end]
            }
        },
        'meta': {
            'apiVersion': 'v2'
        }
    }

    return safe_get_response_json(requests.post(
        f"{config.HTSGET_URL}/beacon/v2/g_variants",
        headers=headers,
        json=payload), 'HTSGet position')

# Figure out whether to use gene search or position search
def query_htsget(headers, gene, assembly, chrom):
    if gene != "":
        return query_htsget_gene(headers, gene)
    else:
        search = re.search(r'(chr)*([XY0-9]{1,2}):(\d+)-(\d+)', chrom)
        return query_htsget_pos(headers, assembly, search.group(2), int(search.group(3)), int(search.group(4)))

# Recursively deep search an object or list for any values less than the aggregate threshold
# NB: does not handle tuples
def censor_response(object):
    if type(object) is list:
        return [censor_response(value) for value in object]
    elif isinstance(object, dict):
        new_dict = {}
        for key, val in object.items():
            new_dict[key] = censor_response(val)
        return new_dict
    elif isinstance(object, int):
        return f"<{config.AGGREGATE_COUNT_THRESHOLD}" if object < config.AGGREGATE_COUNT_THRESHOLD else object

    # Unknown -- leave as-is
    return object

# The return value does not like None being used as a key, so this helper function recursively
# goes through the dictionary provided, and changes all keys to strings
# NB: This overwrites any keys that were previously not strings, and can cause data deletion
# if there was two keys e.g. 12 and "12"
def fix_dicts(to_fix):
    if isinstance(to_fix, dict):
        new_dict = {}
        for key, value in to_fix.items():
            new_dict[str(key)] = fix_dicts(value)
        return new_dict
    elif isinstance(to_fix, list):
        new_list = []
        for value in to_fix:
            new_list.append(fix_dicts(value))
        return new_list
    else:
        return to_fix

# Helper function to format the response from a query into a format the data portal understands
def format_query_response(donors, genomic_query, summary_stats, page, page_size):
    # Determine which part of the filtered donors to send back
    full_data = {
        'results': donors[(page*page_size):((page+1)*page_size)],
        'genomic': genomic_query,
        'count': len(donors),
        'summary': summary_stats
    }
    # full_data['genomic_query_info'] = genomic_query_info

    # Add prev and next parameters to the repsonse, appending a session ID.
    # Essentially we want to go session ID -> list of donors
    # and then paginate the list of donors, calling donors_with_clinical_data on each before returning
    return fix_dicts(full_data), 200

@app.route('/query')
def query(treatment="", primary_site="", drug_name="", systemic_therapy_type="", chrom="", gene="", page=0, page_size=10, assembly="hg38", exclude_programs=[], session_id=""):
    headers = get_headers()

    # NB: We're still doing table joins here, which is probably not where we want to do them
    # We're grabbing (and storing in memory) all the donor data in Katsu with the below request

    # Query the appropriate Katsu endpoint
    url = f"{config.KATSU_URL}/v3/authorized/query/"

    param_mapping = [
        (treatment, "treatment_type"),
        (primary_site, "primary_site"),
        (drug_name, "systemic_therapy_drug_name"),
        (systemic_therapy_type, "systemic_therapy_type"),
        (exclude_programs, "exclude_programs")
    ]
    params = {
        'page_size': PAGE_SIZE
    }
    for param in param_mapping:
        if param[0] == "" or param[0] == []:
            continue
        params[param[1]] = param[0]

    full_url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    donors_req = requests.get(full_url, headers=headers)
    if not donors_req.ok:
        if donors_req.status_code == 401:
            # 401 Unauthorized: the token is invalid
            ret_val = format_query_response([], [], get_summary_stats([], {}, {}), page, page_size)
            return ret_val[0], 401
        else:
            err_msg = f"Could not got Katsu donors response: {donors_req.status_code} {donors_req.text}"
            logger.error(err_msg)
            # Do not forward the response from Katsu in case of compromising information (due to X-Service-Token)
            raise Exception(err_msg)
    donors = donors_req.json()['items']

    # Filter on excluded programs
    donors = [donor for donor in donors if donor['program_id'] not in exclude_programs]

    # Note: We get three extra things from /authorized/query that aren't part of the Donors object:
    # 1) submitter_sample_ids
    # 2) primary_site
    # 3) treatment_type
    # These are used to build the summary information without needing to re-query Katsu
    # For the purposes of the return value, let's remove all three of these into their own variables
    summary_info = {}
    summary_headers = ['submitter_sample_ids', 'primary_site', 'treatment_type']
    for header in summary_headers:
        summary_info[header] = {}
        for donor in donors:
            summary_info[header][donor['submitter_donor_id']] = donor[header]
            del donor[header]

    # Now we combine this with HTSGet, if any
    genomic_query = []
    # genomic_query_info = None
    if gene != "" or chrom != "":
        try:
            htsget = query_htsget(headers, gene, assembly, chrom)

            # We need to be able to map sample registrations, so we'll grab it from Katsu
            samplereg_query_req = requests.get(f"{config.KATSU_URL}/v3/authorized/sample_registrations/?page_size=10000000", headers=headers)
            samplereg_query = safe_get_response_json(samplereg_query_req, 'Katsu sample registrations')
            samplereg_mapping = {}
            for samplereg in samplereg_query['items']:
                samplereg_mapping[samplereg['submitter_sample_id']] = (samplereg['submitter_donor_id'], samplereg['tumour_normal_designation'])

            # genomic_query_info contains ALL matches from every dataset
            # This is meant to be used to fill out the summary stats ONLY
            # However, that part isn't covered in this PR (it's in DIG-1372 (https://candig.atlassian.net/browse/DIG-1372))
            # and does not yet function
            # genomic_query_info = htsget['query_info']
            # for program in genomic_query_info:
            #    sample_ids = genomic_query_info[program]

            htsget_found_donors = {}
            response = htsget['estimatedResults'] if 'estimatedResults' in htsget else {}
            caseLevelData = []
            for program in response.keys():
                if type(response[program]) is not list:
                    continue
                for item in response[program]:
                    submitter_sample_id = item["submitter_sample_id"]
                    case_data = {
                        "program_id": program,
                        "submitter_sample_id": submitter_sample_id,
                        "variant_count": item["variant_count"]
                    }
                    if submitter_sample_id in samplereg_mapping:
                        case_data['donor_id'] = samplereg_mapping[submitter_sample_id][0]
                        case_data['tumour_normal_designation'] = samplereg_mapping[submitter_sample_id][1]
                    else:
                        logger.error(f"Could not find donor mapping for {case_data}")
                        case_data['donor_id'] = submitter_sample_id
                        case_data['tumour_normal_designation'] = 'Tumour'
                    htsget_found_donors[case_data['donor_id']] = 1
                    caseLevelData.append(case_data)

            # Filter clinical results based on genomic results
            donors = [donor for donor in donors if donor['submitter_donor_id'] in htsget_found_donors]
            katsu_allowed_donors = {}
            for donor in donors:
                katsu_allowed_donors[f"{donor['program_id']}~{donor['submitter_donor_id']}"] = 1
            for case_data in caseLevelData:
                if ('donor_id' in case_data and 'program_id' in case_data and
                    f"{case_data['program_id']}~{case_data['donor_id']}" in katsu_allowed_donors):
                    genomic_query.append(case_data)

        except Exception as ex:
            logger.error(f"Error while reading HTSGet response: {ex}")

    # TODO: Cache the above list of donor IDs and summary statistics
    summary_stats = get_summary_stats(donors, summary_info['primary_site'], summary_info['treatment_type'])

    return format_query_response(donors, genomic_query, summary_stats, page, page_size)


def is_discovery_allowed():
    if "X-Service-Token" in connexion.request.headers:
        tokens = connexion.request.headers["X-Service-Token"].split(",")
        for token in tokens:
            if verify_service_token(service="federation", token=token):
                return True, 200
        else:
            return {"error": "Request claims to be from federation but it's not"}, 403
    if not is_user_candig_authorized(connexion.request):
        return {"error": "User is not CanDIG authorized"}, 403
    return True, 200


@app.route('/genomic_completeness')
def genomic_completeness():
    is_allowed, status_code = is_discovery_allowed()
    if status_code != 200:
        return is_allowed, status_code

    headers = get_headers()

    programs = safe_get_response_json(requests.get(f"{config.HTSGET_URL}/ga4gh/drs/v1/programs",
            # Reuse their bearer token
            headers=headers), 'HTSGet programs')
    retVal = {}
    for program_id in programs:
        program = safe_get_response_json(requests.get(f"{config.HTSGET_URL}/ga4gh/drs/v1/programs/{program_id}",
        # Reuse their bearer token
        headers=headers), 'HTSGet program statistics')
        if program_id not in retVal:
            retVal[program_id] = program["statistics"]

    return retVal, 200

@app.route('/discovery/programs')
def discovery_programs():
    is_allowed, status_code = is_discovery_allowed()
    if status_code != 200:
        return is_allowed, status_code

    headers = get_headers()
    # Grab all programs from Katsu
    url = f"{config.KATSU_URL}/v3/discovery/programs/"
    r = safe_get_response_json(requests.get(url, headers=headers), 'Katsu sample registrations')

    # Aggregate all of the programs' return values into one value for the entire site
    site_summary_stats = {
        'required_but_missing': {},
        'summary_cases': {
            'total_cases': 0,
            'complete_cases': 0
        }
    }
    for program in r:
        if 'metadata' not in program:
            logger.error(f"Strange result from Katsu: no metadata in {program}")
            continue
        metadata = program['metadata']

        # There's five metadata categories we care about:
        # schemas_used is a set, schemas_not_used is the inverse of that set
        if 'summary_cases' in metadata:
            try:
                site_summary_stats['summary_cases']['complete_cases'] += metadata['summary_cases']['complete_cases']
                site_summary_stats['summary_cases']['total_cases'] += metadata['summary_cases']['total_cases']
            except:
                logger.error(f"Strange result from Katsu: unreadable summary_cases in {program} metadata")

        if 'required_but_missing' not in metadata:
            # Unreadable result; we cannot continue
            continue
        required_but_missing = metadata['required_but_missing']
        try:
            for field in required_but_missing:
                # Assuming these are of the form 'treatment_setting': {'total': 1, 'missing': 0}
                if field in site_summary_stats['required_but_missing']:
                    for category in required_but_missing[field]:
                        if category in site_summary_stats['required_but_missing'][field]:
                            for instance in required_but_missing[field][category]:
                                site_summary_stats['required_but_missing'][field][category][instance] += required_but_missing[field][category][instance]
                        else:
                            site_summary_stats['required_but_missing'][field][category] = copy.deepcopy(required_but_missing[field][category])
                else:
                    site_summary_stats['required_but_missing'][field] = copy.deepcopy(required_but_missing[field])
        except Exception as ex:
            logger.error(f"Unable to parse required fields result from Katsu: {ex}")

    # Return both the site's aggregated return value and each individual programs'
    ret_val = {
        'site': site_summary_stats,
        'programs': r
    }

    return fix_dicts(ret_val), 200

@app.route('/discovery')
def discovery():
    is_allowed, status_code = is_discovery_allowed()
    if status_code != 200:
        return is_allowed, status_code

    headers = get_headers()
    headers.pop("Authorization", None)

    # Extract from query parameters
    target_service = request.args.get("targetService", "katsu")
    target_path = request.args.get("targetPath")

    if target_service == "katsu":
        url = f"{config.KATSU_URL}/{target_path}"
        response = requests.get(url, headers=headers)

        outheaders = {"Content-Type": "application/json"}

        if response.ok:
            return Response(response=response.text, status=200, headers=outheaders)
        else:
            return {"error": "Failed to fetch data from Katsu"}, response.status_code

    return {"error": "Invalid target service"}, 400

@app.route('/discovery/query')
def discovery_query(treatment="", primary_site="", drug_name="", chrom="", gene="", assembly="hg38", exclude_programs=[]):
    is_allowed, status_code = is_discovery_allowed()
    if status_code != 200:
        return is_allowed, status_code

    url = f"{config.KATSU_URL}/v3/explorer/donors/"
    headers = get_headers()

    param_mapping = [
        (treatment, "treatment_type"),
        (primary_site, "primary_site"),
        (drug_name, "systemic_therapy_drug_name"),
        (exclude_programs, "exclude_programs")
    ]
    params = {
        "page_size": PAGE_SIZE
    }
    for param in param_mapping:
        if param[0] == "" or param[0] == []:
            continue
        params[param[1]] = param[0]

    full_url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    donors = safe_get_response_json(requests.get(full_url, headers=headers), 'Katsu explorer donors')

    # Cross reference with HTSGet, if necessary
    if gene != "" or chrom != "":
        # First, we need to map all Katsu-identified specimens
        samplereg_mapping = {}
        for donor in donors:
            if 'submitter_sample_ids' in donor and type(donor['submitter_sample_ids']) is list:
                for sample_id in donor['submitter_sample_ids']:
                    samplereg_mapping[sample_id] = donor

        try:
            htsget = query_htsget(headers, gene, assembly, chrom)

            htsget_found_donors = {}
            response = htsget['estimatedResults'] if 'estimatedResults' in htsget else {}
            for program in response.keys():
                if type(response[program]) is not list:
                    continue
                for item in response[program]:
                    submitter_sample_id = item["submitter_sample_id"]
                    if submitter_sample_id in samplereg_mapping:
                        donor_id = samplereg_mapping[submitter_sample_id][0]
                        htsget_found_donors[donor_id] = 1
                    else:
                        logger.error(f"Could not find donor mapping for {item}")
                        donor_id = submitter_sample_id
                        htsget_found_donors[donor_id] = 1
            # Filter clinical results based on genomic results
            donors = [donor for donor in donors if f"{donor['program_id']}~{donor['submitter_donor_id']}" in htsget_found_donors]

        except Exception as ex:
            logger.error(f"Error while querying HTSGet: {ex}")

    # Assemble summary statistics
    # NB: Do we need this split up into site-vs-program as well?
    summary_stats = {
        'age_at_diagnosis': {},
        'treatment_type_count': {},
        'primary_site_count': {},
        'patients_per_program': {}
    }
    summary_stat_mapping = [
        ('age_at_diagnosis', 'age_at_diagnosis'),
        ('treatment_type_count', 'treatment_type'),
        ('patients_per_program', 'program_id'),
        ('primary_site_count', 'primary_site')
    ]
    for donor in donors:
        for mapping in summary_stat_mapping:
            if type(donor[mapping[1]]) is list:
                for item in donor[mapping[1]]:
                    add_or_increment(summary_stats[mapping[0]], item)
            else:
                add_or_increment(summary_stats[mapping[0]], donor[mapping[1]])

    # Censor if necessary
    summary_stats = censor_response(summary_stats)

    return fix_dicts(summary_stats), 200

@app.route('/whoami')
def whoami():
    # Grab information about the currently logged-in user
    return { 'key': get_user_id(connexion.request, opa_url = config.OPA_URL) }
