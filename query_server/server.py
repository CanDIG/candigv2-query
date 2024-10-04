from flask_cors import CORS
import connexion
import secrets
import candigv2_logging.logging
from config import PORT, DEBUG_MODE


candigv2_logging.logging.initialize()

# Create the application instance
app = connexion.FlaskApp(__name__, specification_dir='./', options={"swagger_url": "/api"})
app.app.config['SECRET_KEY'] = secrets.token_bytes(32)
CORS(app.app)

app.add_api('openapi.yaml', pythonic_params=True, strict_validation=True)

if __name__ == '__main__':
    app.run(port = PORT, debug=DEBUG_MODE)
