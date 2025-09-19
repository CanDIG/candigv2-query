#!/usr/bin/env bash

sed -i s@\<HTSGET_URL\>@$CANDIG_HTSGET_URL@ config.ini
sed -i s@\<DRS_URL\>@$CANDIG_DRS_URL@ config.ini
sed -i s@\<KATSU_URL\>@$CANDIG_KATSU_URL@ config.ini
sed -i s@\<OPA_URL\>@$OPA_URL@ config.ini

cd query_server
gunicorn -k uvicorn.workers.UvicornWorker server:app
