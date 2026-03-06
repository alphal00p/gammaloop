# IMPORTANT: keep this ordering because check_gammaloop_dependencies can add the gammaloop dependencies to sys.path
from gammaloop import check_gammaloop_dependencies  # nopep8
check_gammaloop_dependencies()  # nopep8
from gammaloop.misc.common import *  # nopep8 # type: ignore

success = register_symbolica()
if not success:
    raise GammaLoopError(
        "Could not register Symbolica license key; this is necessary for running tests.")
