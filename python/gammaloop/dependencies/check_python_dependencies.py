import importlib.metadata
import re
import sys

def parse_semver(version_str):
    """Extracts the major, minor, and patch versions from a semantic versioning string, ignoring pre-release and build metadata."""
    if version_str:
        match = re.match(r'(\d+)\.(\d+)\.(\d+)', version_str)
        if match:
            return tuple(int(part) for part in match.groups())
    return None

def compare_versions(v1, v2):
    """Compares two semantic version tuples."""
    for a, b in zip(v1, v2):
        if a < b:
            return -1
        elif a > b:
            return 1
    return 0

def check_requirement(pkg_name, specifier, version_str):
    """Checks if the installed package meets the version requirement, or if it's installed at all if no version is specified."""
    try:
        installed_version_str = importlib.metadata.version(pkg_name)
        installed_version = parse_semver(installed_version_str)
        required_version = parse_semver(version_str) if version_str else None

        # If no version is specified, the presence of the package is enough.
        if required_version is None:
            return True

        if required_version:
            comparison = compare_versions(installed_version, required_version)
            if specifier in ['==', '==='] and comparison != 0:
                return False
            elif specifier == '!=' and comparison == 0:
                return False
            elif specifier == '>' and comparison <= 0:
                return False
            elif specifier == '<' and comparison >= 0:
                return False
            elif specifier == '>=' and comparison < 0:
                return False
            elif specifier == '<=' and comparison > 0:
                return False
        return True
    except importlib.metadata.PackageNotFoundError:
        return False

def check_requirements(requirements_path='../requirements.txt'):
    with open(requirements_path, 'r') as f:
        lines = f.readlines()

    problems = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            # Adjusted regex to make version specifier and number optional
            match = re.match(r'([a-zA-Z0-9_-]+)([><=!]=?)?(\d+(?:\.\d+){0,2})?', line)
            if match:
                pkg_name, specifier, version_str = match.groups()
                if not check_requirement(pkg_name, specifier or '', version_str):
                    problem = f'{pkg_name}'
                    if specifier and version_str:
                        problem += f'{specifier}{version_str}'
                    problems.append(problem)

    if problems:
        print("Packages that do not meet the version requirements or are not installed:")
        for problem in problems:
            print(problem)
        sys.exit(1)
    else:
        print("All packages meet the version requirements or are installed.")
        sys.exit(0)

# Example usage
check_requirements()

