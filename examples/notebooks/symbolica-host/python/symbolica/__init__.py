import os as _os

# Select the installed, authenticated fonts before native render caches initialize.
_os.environ.setdefault(
    "SYMBOLICA_DATA_DIRECTORY",
    _os.path.join(_os.path.dirname(__file__), "_data", "fonts"),
)

from symbolica.core import *  # noqa: E402, F403


def __getattr__(name):
    # `import symbolica.community.hepkit.vakint as v` after `import symbolica`
    # finds the native module already registered in sys.modules and skips its
    # parent packages; import a parent on first access instead.
    if name == "community":
        import importlib

        return importlib.import_module("symbolica.community")
    raise AttributeError(f"module 'symbolica' has no attribute {name!r}")
