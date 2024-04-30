from flask_cors import CORS
import connexion
import secrets
import logging
from config import PORT, DEBUG_MODE

# Create the application instance
app = connexion.FlaskApp(__name__, specification_dir='./')
app.app.config['SECRET_KEY'] = secrets.token_bytes(32)
CORS(app.app)

app.add_api('openapi.yaml', pythonic_params=True, options={"swagger_url": "/api"})

if __name__ == '__main__':
    logging.basicConfig(filename='record.log', level=logging.DEBUG, format=f'%(asctime)s %(levelname)s %(name)s %(threadName)s : %(message)s')
    app.run(port = PORT, debug=DEBUG_MODE)
