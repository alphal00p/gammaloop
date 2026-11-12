"""Community modules of the shared notebook host (graph, HEPKit, tensor)."""


def __getattr__(name):
    # See symbolica.__getattr__: parents of a native submodule that is already
    # registered in sys.modules are imported on first access.
    if name in ("graph", "hepkit", "tensor"):
        import importlib

        return importlib.import_module(f"symbolica.community.{name}")
    raise AttributeError(f"module 'symbolica.community' has no attribute {name!r}")
