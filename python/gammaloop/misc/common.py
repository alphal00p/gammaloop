import os
import sys
import logging
import symbolica
from enum import Enum
from gammaloop.misc.utils import setup_logging

GL_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.path.pardir))
DATA_PATH = os.path.abspath(os.path.join(GL_PATH,'data'))
GL_CONSOLE_HANDLER = setup_logging()

global GL_IS_SYMBOLICA_REGISTERED
global GL_DEBUG
GL_DEBUG = False
GL_IS_SYMBOLICA_REGISTERED = None

logger = logging.getLogger('GammaLoop')

try:
    from symbolica import Expression, Transformer, set_license_key
except ImportError:
    logger.critical("Could not import Symbolica, please install symbolica with 'python -m pip install symbolica'.")
    sys.exit(1)

def register_symbolica() -> bool:
    global GL_IS_SYMBOLICA_REGISTERED
    # Check if is_licensed exist for backward compatibility. Later on, it can be removed.
    if GL_IS_SYMBOLICA_REGISTERED is not None or (hasattr(symbolica,'is_licensed') and symbolica.is_licensed()):
        GL_IS_SYMBOLICA_REGISTERED = True
        return GL_IS_SYMBOLICA_REGISTERED

    if 'SYMBOLICA_LICENSE' not in os.environ:
        license_file_path = os.path.join(DATA_PATH,'config','symbolica_license.txt')
        if os.path.isfile(license_file_path):
            with open(license_file_path,'r') as f:
                symbolica_license = f.read().split('\n')[0].strip()
                if symbolica_license != "<PASE_YOUR_SYMBOLICA_LICENSE_HERE>":
                    try:
                        set_license_key(symbolica_license)
                        GL_IS_SYMBOLICA_REGISTERED = True
                        return True
                    except:
                        logger.critical(f"Could not set Symbolica license key from file '{license_file_path}'.")
                        GL_IS_SYMBOLICA_REGISTERED = True
                        return False
        else:
            GL_IS_SYMBOLICA_REGISTERED = False
            return False
    else:
        try:
            set_license_key(os.environ['SYMBOLICA_LICENSE'])
            GL_IS_SYMBOLICA_REGISTERED = True
            return True
        except:
            logger.critical("Could not set Symbolica license key from environment variable 'SYMBOLICA_LICENSE'.")
            GL_IS_SYMBOLICA_REGISTERED = False
            return False

class Side(Enum):
    LEFT = 0
    RIGHT = 1

class GammaLoopError(Exception):
    pass