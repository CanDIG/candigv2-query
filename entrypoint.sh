#!/usr/bin/env bash

if [[ -f "initial_setup" ]]; then
    sed -i s@\<HTSGET_URL\>@$CANDIG_HTSGET_URL@ config.ini
    sed -i s@\<KATSU_URL\>@$CANDIG_KATSU_URL@ config.ini
    sed -i s@\<OPA_URL\>@$OPA_URL@ config.ini

    rm initial_setup
fi

cd query_server
gunicorn server:app