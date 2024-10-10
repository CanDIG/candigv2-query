# Query Microservice

Query API for Katsu and HTSGet that allows for caching and pagination of calls.

[openapi spec on Swagger](https://editor.swagger.io/?url=https://raw.githubusercontent.com/CanDIG/query/main/query_server/openapi.yaml)

## Stack

- [Flask](http://flask.pocoo.org/)

## Installation

The server software can be installed in a virtual environment:
```
python setup.py install
```

## Running

This application can be configured by way of the config.ini file in the root of the project.
The server can be run with: 

```
python query_server/server.py
```

This application can also be set up in a docker container. A docker-compose file and Dockerfile are provided.
