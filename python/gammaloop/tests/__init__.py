from gammaloop.misc.common import *

success = register_symbolica()
if not success:
    raise GammaLoopError("Could not register Symbolica license key; this is necessary for running tests.")