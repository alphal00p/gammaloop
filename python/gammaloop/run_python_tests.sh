#!/usr/bin/env bash
python -m pytest -m "not slow" --runrust --codecheck "$@"
#python -m pytest tests/unit/* "$@"
#python -m pytest tests/integration/* "$@"
