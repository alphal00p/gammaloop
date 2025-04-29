#/usr/bin/env python3
import os
import sys
import importlib.metadata
from packaging.requirements import Requirement

requirements_file_path = os.path.abspath( os.path.join( os.path.dirname(__file__), os.path.pardir, 'requirements.txt') )

class MissingDependency(Exception):
    pass

def check_requirements(file_path):
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            # Skip comments and blank lines
            if not line or line.startswith('#'):
                continue
            req = Requirement(line)
            try:
                installed_version = importlib.metadata.version(req.name)
                if req.specifier and not req.specifier.contains(installed_version, prereleases=True):
                    raise MissingDependency(f"{req.name} version {installed_version} does not meet {req.specifier}")
                else:
                    print(f"{req.name} version {installed_version} is OK")
            except importlib.metadata.PackageNotFoundError:
                raise MissingDependency(f"{req.name} is not installed.")

check_requirements(requirements_file_path)
