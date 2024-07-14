#!/usr/bin/env bash

RETCODE=0;
python -m pytest -m "not slow" --runrust --codecheck --durations=10 "$@"
RETCODE=$RETCODE+$?;

#python -m pytest tests/unit/* "$@"
#python -m pytest tests/integration/* "$@"
exit $(($RETCODE))
