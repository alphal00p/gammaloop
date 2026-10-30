"""Inject the private license only immediately before exec of a direct runtime ELF."""

import os
import sys
from pathlib import Path

env = dict(os.environ)
assert "SYMBOLICA_LICENSE" not in env and "SYMBOLICA_LICENSE_SIGNED" not in env
config = Path(env.get("XDG_CONFIG_HOME", str(Path.home() / ".config")))
license_value = (config / "symbolica/license").read_text().strip()
assert license_value
env["SYMBOLICA_LICENSE"] = license_value
os.execve(sys.argv[1], sys.argv[1:], env)
