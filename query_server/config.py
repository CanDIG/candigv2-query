import configparser
import os
from authx.auth import create_service_token
from candigv2_logging.logging import CanDIGLogger


logger = CanDIGLogger(__file__)

config = configparser.ConfigParser(interpolation=None)
config.read(os.path.abspath(f"{os.path.dirname(os.path.realpath(__file__))}/../config.ini"))

AUTHZ = config['authz']
QUERY_URL = os.getenv("QUERY_URL", f"http://localhost:{config['DEFAULT']['Port']}")
AGGREGATE_COUNT_THRESHOLD = int(os.getenv("AGGREGATE_COUNT_THRESHOLD", "5"))

PORT = config['DEFAULT']['Port']

HTSGET_URL = config['DEFAULT']['CANDIG_HTSGET_URL']
KATSU_URL = config['DEFAULT']['CANDIG_KATSU_URL']
OPA_URL = config['authz']['CANDIG_OPA_URL']

DEBUG_MODE = False
if os.getenv("DEBUG_MODE", "1") == "1":
    DEBUG_MODE = True

service_token_path = os.path.abspath(f"{os.path.dirname(os.path.realpath(__file__))}/../service_token")
try:
    # We will need to create a long-lived service token to ensure that every request is coming from us
    # To prevent concurrency issues, we'll generate one during startup and use it for every request (nb: insecure?)
    if not os.path.isfile(service_token_path):
        with open(service_token_path, "w") as service_token_file:
            SERVICE_TOKEN = create_service_token()
            service_token_file.write(SERVICE_TOKEN)
    else:
        with open(service_token_path, "r") as service_token_file:
            SERVICE_TOKEN = service_token_file.read()
    if DEBUG_MODE:
        print(f"SERVICE_TOKEN: {SERVICE_TOKEN}")

except:
    logger.error("Could not obtain a service token")
    SERVICE_TOKEN = ""
