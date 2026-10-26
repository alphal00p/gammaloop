import os as _os

# Select the installed, authenticated fonts before native render caches initialize.
_os.environ.setdefault(
    "SYMBOLICA_DATA_DIRECTORY",
    _os.path.join(_os.path.dirname(__file__), "_data", "fonts"),
)

from symbolica.core import *  # noqa: E402, F403
