import configparser
import os
from authx.auth import create_service_token

config = configparser.ConfigParser(interpolation=None)
config.read(os.path.abspath(f"{os.path.dirname(os.path.realpath(__file__))}/../config.ini"))

AUTHZ = config['authz']
QUERY_URL = os.getenv("QUERY_URL", f"http://localhost:{config['DEFAULT']['Port']}")

PORT = config['DEFAULT']['Port']

HTSGET_URL = config['DEFAULT']['CANDIG_HTSGET_URL']
KATSU_URL = config['DEFAULT']['CANDIG_KATSU_URL']

DEBUG_MODE = False
if os.getenv("DEBUG_MODE", "1") == "1":
    DEBUG_MODE = True

try:
    SERVICE_TOKEN = create_service_token()
except:
    print("Could not obtain a service token")
    SERVICE_TOKEN = ""
