import os
import sys
import logging
from enum import Enum
import yaml
import symbolica  # pylint: disable=import-error
from gammaloop.misc.utils import setup_logging

pjoin = os.path.join

GL_PATH = os.path.abspath(os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.path.pardir))
DATA_PATH = os.path.abspath(os.path.join(GL_PATH, 'data'))
GL_CONSOLE_HANDLER = setup_logging()

GL_DEBUG = False
GL_IS_SYMBOLICA_REGISTERED = None

GAMMALOOP_CONFIG_PATHS = [
    pjoin(DATA_PATH, 'config', 'gammaloop_config.yaml'),
    pjoin(os.path.expanduser('~'), '.gammaloop_config.yaml'),
]

logger = logging.getLogger('GammaLoop')

try:
    # pylint: disable=unused-import
    from symbolica import Expression, Transformer, set_license_key
except ImportError:
    logger.critical(
        "Could not import Symbolica, please install symbolica with 'python -m pip install symbolica'.")
    sys.exit(1)


def load_configuration(config_path, quiet: bool = False) -> dict:
    if not quiet:
        logger.info(
            "Updating gammaLoop configuration from '%s'.", config_path)
    try:
        gammaloop_config = yaml.load(open(config_path, 'r', encoding='utf-8'),
                                     Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        logger.critical(
            "Could not load gammaLoop configuration from '%s'. Error: \n%s", config_path, str(exc))
        sys.exit(1)
    return gammaloop_config


def register_symbolica() -> bool:
    # pylint: disable=global-statement
    global GL_IS_SYMBOLICA_REGISTERED
    # Check if is_licensed exist for backward compatibility. Later on, it can be removed.
    if GL_IS_SYMBOLICA_REGISTERED is not None or (hasattr(symbolica, 'is_licensed') and symbolica.is_licensed()):
        GL_IS_SYMBOLICA_REGISTERED = True
        return GL_IS_SYMBOLICA_REGISTERED

    if 'SYMBOLICA_LICENSE' not in os.environ:
        symbolica_license = None
        for path in GAMMALOOP_CONFIG_PATHS:
            if os.path.exists(path):
                gammaloop_config = load_configuration(path, quiet=True)
                if 'symbolica_license' in gammaloop_config and gammaloop_config['symbolica_license'] != "<PASTE_YOUR_SYMBOLICA_LICENSE_HERE>":
                    symbolica_license = gammaloop_config['symbolica_license']
                    break
        if symbolica_license is not None:
            try:
                set_license_key(symbolica_license)
                GL_IS_SYMBOLICA_REGISTERED = True
                return True
            except:  # pylint: disable=bare-except
                logger.critical(
                    "Could not set Symbolica license key from gammaloop configuration file specifying license '%s'.", symbolica_license)
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
        except:  # pylint: disable=bare-except
            logger.critical(
                "Could not set Symbolica license key from environment variable 'SYMBOLICA_LICENSE' specifying license '%s'.", os.environ['SYMBOLICA_LICENSE'])
            GL_IS_SYMBOLICA_REGISTERED = False
            return False


class Side(Enum):
    LEFT = 0
    RIGHT = 1

    def __str__(self):
        if self == Side.LEFT:
            return "LEFT"
        elif self == Side.RIGHT:
            return "RIGHT"
        else:
            raise GammaLoopError(f"Unknown side: {self}")


class GammaLoopError(Exception):
    pass
