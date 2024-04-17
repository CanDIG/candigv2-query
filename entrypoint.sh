#!/usr/bin/env bash
OPA_SECRET=$(cat /run/secrets/opa-service-token)

if [[ -f "initial_setup" ]]; then
    sed -i s@\<HTSGET_URL\>@$CANDIG_HTSGET_URL@ config.ini
    sed -i s@\<KATSU_URL\>@$CANDIG_KATSU_URL@ config.ini

    rm initial_setup
fi

uwsgi uwsgi.ini
