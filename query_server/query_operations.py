from flask import request, Flask
import config
import copy
import re
import requests
import secrets
import urllib
from authx.auth import get_user_id, get_auth_token
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

def safe_get_request_json(request, name):
    if not request.ok:
        raise Exception(f"Could not get {name} response: {request.status_code} {request.text}")
    return request.json()

# Grab a list of donors matching a given filter from the given URL
def get_donors_from_katsu(url, param_name, parameter_list, headers):
    permissible_donors = set()
    for parameter in parameter_list:
        # TODO: Fix the page_size call here -- use a consume_all() query like in the frontend
        parameters = {
            param_name: parameter,
            'page_size': PAGE_SIZE
        }
        treatments = requests.get(f"{url}?{urllib.parse.urlencode(parameters)}", headers=headers)
        results = safe_get_request_json(treatments, f'Katsu {param_name}')['items']
        permissible_donors |= set([result['submitter_donor_id'] for result in results])
    return permissible_donors

def add_or_increment(dict, key):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1

def get_summary_stats(donors, headers):
    # Perform (and cache) summary statistics
    age_at_diagnosis = {}
    donors_by_id = {}
    primary_site_count = {}
    patients_per_cohort = {}
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

        # Cancer types
        if donor['primary_site']:
            for cancer_type in donor['primary_site']:
                if cancer_type in primary_site_count:
                    primary_site_count[cancer_type] += 1
                else:
                    primary_site_count[cancer_type] = 1
        program_id = donor['program_id']
        if program_id in patients_per_cohort:
            patients_per_cohort[program_id] += 1
        elif program_id is not None:
            patients_per_cohort[program_id] = 1

    # Treatment types
    # http://candig.docker.internal:8008/v2/authorized/treatments/
    treatments = requests.get(f"{config.KATSU_URL}/v2/authorized/treatments/?page_size={PAGE_SIZE}",
        headers=headers)
    treatments = safe_get_request_json(treatments, 'Katsu treatments')['items']
    treatment_type_count = {}
    for treatment in treatments:
        if (treatment["submitter_donor_id"] in donors_by_id and
            "treatment_type" in treatment and
            treatment["treatment_type"] is not None):
            try:
                for treatment_type in treatment["treatment_type"]:
                    add_or_increment(treatment_type_count, treatment_type)
            except TypeError as e:
                logger.log_message("ERROR", f"Could not grab summary treatment stats: {e}")
                pass

    return {
        'age_at_diagnosis': age_at_diagnosis,
        'treatment_type_count': treatment_type_count,
        'primary_site_count': primary_site_count,
        'patients_per_cohort': patients_per_cohort
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

        return safe_get_request_json(requests.post(
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

    return safe_get_request_json(requests.post(
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

@app.route('/query')
def query(treatment="", primary_site="", chemotherapy="", immunotherapy="", hormone_therapy="", chrom="", gene="", page=0, page_size=10, assembly="hg38", exclude_cohorts=[], session_id=""):
    # Add a service token to the headers so that other services will know this is from the query service:
    headers = {}
    for k in request.headers.keys():
        headers[k] = request.headers[k]
    headers["X-Service-Token"] = config.SERVICE_TOKEN

    # NB: We're still doing table joins here, which is probably not where we want to do them
    # We're grabbing (and storing in memory) all the donor data in Katsu with the below request

    # Query the appropriate Katsu endpoint
    params = { 'page_size': PAGE_SIZE }
    url = f"{config.KATSU_URL}/v2/authorized/donors/"
    if primary_site != "":
        params['primary_site'] = primary_site
    r = safe_get_request_json(requests.get(f"{url}?{urllib.parse.urlencode(params, True)}",
        # Reuse their bearer token
        headers=headers), 'Katsu Donors')
    donors = r['items']

    # Filter on excluded cohorts
    donors = [donor for donor in donors if donor['program_id'] not in exclude_cohorts]

    # Will need to look into how to go about this -- ideally we implement this into the SQL in Katsu's side
    filters = [
        (treatment, f"{config.KATSU_URL}/v2/authorized/treatments/", 'treatment_type'),
        (chemotherapy, f"{config.KATSU_URL}/v2/authorized/chemotherapies/", 'drug_name'),
        (immunotherapy, f"{config.KATSU_URL}/v2/authorized/immunotherapies/", 'drug_name'),
        (hormone_therapy, f"{config.KATSU_URL}/v2/authorized/hormone_therapies/", 'drug_name')
    ]
    for (this_filter, url, param_name) in filters:
        if this_filter != "":
            permissible_donors = get_donors_from_katsu(
                url,
                param_name,
                this_filter,
                headers
            )
            donors = [donor for donor in donors if donor['submitter_donor_id'] in permissible_donors]

    # Now we combine this with HTSGet, if any
    genomic_query = []
    # genomic_query_info = None
    if gene != "" or chrom != "":
        try:
            htsget = query_htsget(headers, gene, assembly, chrom)

            # We need to be able to map specimens, so we'll grab it from Katsu
            specimen_query_req = requests.get(f"{config.KATSU_URL}/v2/authorized/sample_registrations/?page_size=10000000", headers=headers)
            specimen_query = safe_get_request_json(specimen_query_req, 'Katsu sample registrations')
            specimen_mapping = {}
            for specimen in specimen_query['items']:
                specimen_mapping[specimen['submitter_sample_id']] = (specimen['submitter_donor_id'], specimen['tumour_normal_designation'])

            # genomic_query_info contains ALL matches from every dataset
            # This is meant to be used to fill out the summary stats ONLY
            # However, that part isn't covered in this PR (it's in DIG-1372 (https://candig.atlassian.net/browse/DIG-1372))
            # and does not yet function
            # genomic_query_info = htsget['query_info']
            # for cohort in genomic_query_info:
            #    sample_ids = genomic_query_info[cohort]

            htsget_found_donors = {}
            for response in htsget['response']:
                for case_data in response['caseLevelData']:
                    if 'biosampleId' not in case_data:
                        logger.log_message("ERROR", f"Could not parse htsget response for {case_data}")
                        continue
                    id = case_data['biosampleId'].split('~')
                    if len(id) > 1:
                        case_data['program_id'] = id[0]
                        submitter_specimen_id = id[1]
                        case_data['submitter_specimen_id'] = submitter_specimen_id
                        if submitter_specimen_id in specimen_mapping:
                            case_data['donor_id'] = specimen_mapping[submitter_specimen_id][0]
                            case_data['tumour_normal_designation'] = specimen_mapping[submitter_specimen_id][1]
                        else:
                            logger.log_message("ERROR", f"Could not find donor mapping for {case_data}")
                            case_data['donor_id'] = submitter_specimen_id
                            case_data['tumour_normal_designation'] = 'Tumour'
                        htsget_found_donors[case_data['donor_id']] = 1
                    else:
                        logger.log_message("ERROR", f"Could not parse biosampleId for {case_data}")
                        case_data['program_id'] = None
                        case_data['donor_id'] = None
                        case_data['submitter_specimen_id'] = case_data['biosampleId']
                        case_data['tumour_normal_designation'] = 'Tumour'
                    case_data['position'] = response['variation']['location']['interval']['start']['value']
            # Filter clinical results based on genomic results
            donors = [donor for donor in donors if donor['submitter_donor_id'] in htsget_found_donors]
            katsu_allowed_donors = {}
            for donor in donors:
                katsu_allowed_donors[f"{donor['program_id']}~{donor['submitter_donor_id']}"] = 1
            for response in htsget['response']:
                for case_data in response['caseLevelData']:
                    if ('donor_id' in case_data and 'program_id' in case_data and
                        f"{case_data['program_id']}~{case_data['donor_id']}" in katsu_allowed_donors):
                        genomic_query.append(case_data)

        except Exception as ex:
            logger.log_message("ERROR", f"Error while reading HTSGet response: {ex}")

    # TODO: Cache the above list of donor IDs and summary statistics
    summary_stats = get_summary_stats(donors, headers)

    # Determine which part of the filtered donors to send back
    full_data = {
        'results': [donor for donor in donors[(page*page_size):((page+1)*page_size)]],
        'genomic': genomic_query,
        'count': len(donors),
        'summary': summary_stats
    }
    # full_data['genomic_query_info'] = genomic_query_info

    # Add prev and next parameters to the repsonse, appending a session ID.
    # Essentially we want to go session ID -> list of donors
    # and then paginate the list of donors, calling donors_with_clinical_data on each before returning
    return fix_dicts(full_data), 200

@app.route('/genomic_completeness')
def genomic_completeness():
    # Add a service token to the headers so that Katsu will know this is from the query service:
    headers = {}
    for k in request.headers.keys():
        headers[k] = request.headers[k]
    headers["X-Service-Token"] = config.SERVICE_TOKEN

    samples = safe_get_request_json(requests.get(f"{config.HTSGET_URL}/htsget/v1/samples",
            # Reuse their bearer token
            headers=headers), 'HTSGet cohort statistics')

    retVal = {}
    for sample in samples:
        program_id = sample['cohort']
        if program_id not in retVal:
            retVal[program_id] = { 'genomes': 0, 'transcriptomes': 0, 'all': 0 }
        if len(sample['genomes']) > 0 and len(sample['transcriptomes']) > 0:
            retVal[program_id]['all'] += 1
        if len(sample['genomes']) > 0:
            retVal[program_id]['genomes'] += 1
        if len(sample['transcriptomes']) > 0:
            retVal[program_id]['transcriptomes'] += 1

    return retVal, 200

@app.route('/discovery/programs')
def discovery_programs():
    # Grab all programs from Katsu
    url = f"{config.KATSU_URL}/v2/discovery/programs/"
    r = safe_get_request_json(requests.get(url), 'Katsu sample registrations')

    # Aggregate all of the programs' return values into one value for the entire site
    site_summary_stats = {
        'schemas_used': set(),
        'schemas_not_used': set(),
        'required_but_missing': {},
        'cases_missing_data': set(),
        'summary_cases': {
            'total_cases': 0,
            'complete_cases': 0
        }
    }
    unused_schemas = set()
    unused_initialized = False
    for program in r:
        if 'metadata' not in program:
            logger.log_message("ERROR", f"Strange result from Katsu: no metadata in {program}")
            continue
        metadata = program['metadata']

        # There's five metadata categories we care about:
        # schemas_used is a set, schemas_not_used is the inverse of that set
        if not unused_initialized:
            unused_initialized = True
            unused_schemas = set(metadata['schemas_not_used'])
        if 'schemas_used' in metadata:
            site_summary_stats['schemas_used'] |= set(metadata['schemas_used'])
        if 'cases_missing_data' in metadata:
            site_summary_stats['cases_missing_data'] |= set(metadata['cases_missing_data'])
        if 'summary_cases' in metadata:
            try:
                site_summary_stats['summary_cases']['complete_cases'] += metadata['summary_cases']['complete_cases']
                site_summary_stats['summary_cases']['total_cases'] += metadata['summary_cases']['total_cases']
            except:
                logger.log_message("ERROR", f"Strange result from Katsu: unreadable summary_cases in {program} metadata")

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
            logger.log_message("ERROR", f"Unable to parse required fields result from Katsu: {ex}")

    for schema in site_summary_stats['schemas_used']:
        unused_schemas.discard(schema)
    site_summary_stats['schemas_not_used'] = list(unused_schemas)
    site_summary_stats['schemas_used'] = list(site_summary_stats['schemas_used'])
    site_summary_stats['cases_missing_data'] = list(site_summary_stats['cases_missing_data'])

    # Return both the site's aggregated return value and each individual programs'
    ret_val = {
        'site': site_summary_stats,
        'programs': r
    }

    return fix_dicts(ret_val), 200

@app.route('/discovery/query')
def discovery_query(treatment="", primary_site="", chemotherapy="", immunotherapy="", hormone_therapy="", chrom="", gene="", assembly="hg38", exclude_cohorts=[]):
    url = f"{config.KATSU_URL}/v2/explorer/donors/"
    headers = {}
    for k in request.headers.keys():
        headers[k] = request.headers[k]
    headers["X-Service-Token"] = config.SERVICE_TOKEN

    param_mapping = [
        (treatment, "treatment_type"),
        (primary_site, "primary_site"),
        (chemotherapy, "chemotherapy_drug_name"),
        (immunotherapy, "immunotherapy_drug_name"),
        (hormone_therapy, "hormone_therapy_drug_name"),
        (exclude_cohorts, "exclude_cohorts")
    ]
    params = {}
    for param in param_mapping:
        if param[0] == "" or param[0] == []:
            continue
        params[param[1]] = param[0]

    full_url = f"{url}?{urllib.parse.urlencode(params, doseq=True)}"
    donors = safe_get_request_json(requests.get(full_url, headers=headers), 'Katsu explorer donors')

    # Cross reference with HTSGet, if necessary
    if gene != "" or chrom != "":
        # First, we need to map all Katsu-identified specimens
        specimen_mapping = {}
        for donor in donors:
            if 'submitter_sample_ids' in donor and type(donor['submitter_sample_ids']) is list:
                for sample_id in donor['submitter_sample_ids']:
                    specimen_mapping[f"{donor['program_id']}~{sample_id}"] = donor

        try:
            htsget = query_htsget(headers, gene, assembly, chrom)

            htsget_found_donors = {}
            for program_id in htsget['query_info']:
                for sample_id in htsget['query_info'][program_id]:
                    # NB: We're allowing the entire donor as long as any specimen matches -- is that what we want?
                    merged_id = f"{program_id}~{sample_id}"
                    if merged_id in specimen_mapping:
                        found_donor = specimen_mapping[merged_id]
                        htsget_found_donors[f"{found_donor['program_id']}~{found_donor['submitter_donor_id']}"] = 1
                    else:
                        logger.log_message("ERROR", f"Could not find specimen identified in HTSGet: {merged_id}")
            # Filter clinical results based on genomic results
            donors = [donor for donor in donors if f"{donor['program_id']}~{donor['submitter_donor_id']}" in htsget_found_donors]

        except Exception as ex:
            logger.log_message("ERROR", f"Error while querying HTSGet: {ex}")

    # Assemble summary statistics
    # NB: Do we need this split up into site-vs-program as well?
    summary_stats = {
        'age_at_diagnosis': {},
        'treatment_type_count': {},
        'primary_site_count': {},
        'patients_per_cohort': {}
    }
    summary_stat_mapping = [
        ('age_at_diagnosis', 'age_at_diagnosis'),
        ('treatment_type_count', 'treatment_type'),
        ('patients_per_cohort', 'program_id'),
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
    logger.log_message("DEBUG", config.OPA_URL)
    logger.log_message("DEBUG", config.AUTHZ)
    token = get_auth_token(request)
    headers = {
        "Authorization": f"Bearer {token}"
    }
    response = requests.post(
        config.OPA_URL + f"/v1/data/idp/user_key",
        headers=headers,
        json={
            "input": {
                    "token": token
                }
            }
        )
    logger.log_message("DEBUG", response)
    return { 'key': get_user_id(request, opa_url = config.OPA_URL) }
