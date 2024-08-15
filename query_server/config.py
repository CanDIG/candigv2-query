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

try:
    SERVICE_TOKEN = create_service_token()
except:
    logger.log_message("ERROR", "Could not obtain a service token")
    SERVICE_TOKEN = ""
