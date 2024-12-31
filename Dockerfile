ARG venv_python
FROM python:${venv_python}

LABEL Maintainer="CanDIG Project"
LABEL "candigv2"="query_app"

USER root

RUN groupadd -r candig && useradd -rm candig -g candig

RUN apt-get update

COPY requirements.txt /app/query_server/requirements.txt

RUN pip install --no-cache-dir -r /app/query_server/requirements.txt

COPY . /app/query_server

WORKDIR /app/query_server

RUN chown -R candig:candig /app/query_server

USER candig

ENTRYPOINT ["bash", "entrypoint.sh"]
