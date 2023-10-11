import os
import sys
import logging
from enum import Enum
import yaml
from typing import Any
import symbolica  # pylint: disable=import-error # type: ignore
from gammaloop.misc.utils import setup_logging

pjoin = os.path.join

GL_PATH = os.path.abspath(os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.path.pardir))
DATA_PATH = os.path.abspath(os.path.join(GL_PATH, 'data'))
GL_CONSOLE_HANDLER = setup_logging()

GL_DEBUG = False
gl_is_symbolica_registered = None

GAMMALOOP_CONFIG_PATHS = [
    pjoin(DATA_PATH, 'config', 'gammaloop_config.yaml'),
    pjoin(os.path.expanduser('~'), '.gammaloop_config.yaml'),
]

# Useful for pyLance to work fine when doing something like sum([[1,2],[4,5]],[]) to flatten the nested list
EMPTY_LIST: list[Any] = []

logger = logging.getLogger('GammaLoop')

try:
    # pylint: disable=unused-import
    from symbolica import Expression, Transformer, set_license_key  # type: ignore
except ImportError:
    logger.critical(
        "Could not import Symbolica, please install symbolica with 'python -m pip install symbolica'.")
    sys.exit(1)


def load_configuration(config_path: str, quiet: bool = False) -> dict[str, Any]:
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
    global gl_is_symbolica_registered
    # Check if is_licensed exist for backward compatibility. Later on, it can be removed.
    if gl_is_symbolica_registered is not None or (hasattr(symbolica, 'is_licensed') and symbolica.is_licensed()):
        gl_is_symbolica_registered = True
        return gl_is_symbolica_registered

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
                gl_is_symbolica_registered = True
                return True
            except:  # pylint: disable=bare-except
                logger.critical(
                    "Could not set Symbolica license key from gammaloop configuration file specifying license '%s'.", symbolica_license)
                gl_is_symbolica_registered = True
                return False
        else:
            gl_is_symbolica_registered = False
            return False
    else:
        try:
            set_license_key(os.environ['SYMBOLICA_LICENSE'])
            gl_is_symbolica_registered = True
            return True
        except:  # pylint: disable=bare-except
            logger.critical(
                "Could not set Symbolica license key from environment variable 'SYMBOLICA_LICENSE' specifying license '%s'.", os.environ['SYMBOLICA_LICENSE'])
            gl_is_symbolica_registered = False
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
