#!/usr/bin/env bash

RETCODE=0;
python3 -m pytest -m "not slow" --runrust --codecheck "$@"
RETCODE=$RETCODE+$?;

#python -m pytest tests/unit/* "$@"
#python -m pytest tests/integration/* "$@"
exit $(($RETCODE))
