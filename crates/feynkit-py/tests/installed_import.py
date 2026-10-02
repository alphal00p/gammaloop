"""Import must leave license activation to the caller's first Symbolica operation."""

import os
import subprocess
import sys

environment = os.environ.copy()
environment.pop("SYMBOLICA_LICENSE", None)
environment.pop("SYMBOLICA_HIDE_BANNER", None)

# A fresh process is essential: another test may already have initialized the
# global license manager. Importing E also registers the native community classes.
result = subprocess.run(
    [
        sys.executable,
        "-u",
        "-c",
        (
            "from symbolica import E, set_license_key; "
            "from symbolica.community import hepkit; "
            "print('imported', flush=True); "
            "E('x')"
        ),
    ],
    env=environment,
    capture_output=True,
    text=True,
    check=True,
    timeout=60,
)
assert result.stdout.startswith("imported\n"), result.stdout
assert "restricted Symbolica instance" in result.stdout, result.stdout
assert not result.stderr, result.stderr
print("Import defers the license check until the first symbolic operation")
