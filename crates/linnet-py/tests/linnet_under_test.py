"""The installed Linnet bindings this suite tests.

The standalone ``linnet-py`` wheel provides ``linnet_py``; the Symbolica
Community wheel provides ``symbolica.community.graph``. The suite tests
whichever is installed; set ``LINNET_PY_MODULE`` to one of these names when
both are.
"""

import importlib
import importlib.util
import os
import unittest

MODULES = ("linnet_py", "symbolica.community.graph")


def _installed(name):
    try:
        return importlib.util.find_spec(name) is not None
    except ModuleNotFoundError:
        return False


def _select():
    requested = os.environ.get("LINNET_PY_MODULE")
    if requested in MODULES:
        return requested
    installed = [name for name in MODULES if _installed(name)]
    if requested or len(installed) != 1:
        raise ImportError(
            f"set LINNET_PY_MODULE to one of {MODULES}; installed: {installed}"
        )
    return installed[0]


#: The module every class reports, imported by the tests as ``graph``.
NAME = _select()
graph = importlib.import_module(NAME)
#: Only the Community build accepts Symbolica expressions as keys of their own kind.
COMMUNITY = NAME == "symbolica.community.graph"
community_only = unittest.skipUnless(COMMUNITY, "needs the Symbolica Community wheel")

#: Builds a field label: a Symbolica symbol where the build accepts expressions.
if COMMUNITY:
    from symbolica import S as field
else:
    field = str
