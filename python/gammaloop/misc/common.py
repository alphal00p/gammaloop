import os
import sys
import logging
import subprocess
from enum import Enum
import yaml
from typing import Any
from gammaloop.misc.utils import setup_logging, Colour  # type: ignore
from gammaloop import check_gammaloop_dependencies
pjoin = os.path.join

GL_PATH = os.path.abspath(os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.path.pardir))
DATA_PATH = os.path.abspath(os.path.join(GL_PATH, 'data'))
GL_CONSOLE_HANDLER = setup_logging()

global GL_DEBUG
GL_DEBUG = False

gl_is_symbolica_registered = None

GAMMALOOP_CONFIG_PATHS = [
    pjoin(DATA_PATH, 'config', 'gammaloop_config.yaml'),
    pjoin(os.path.expanduser('~'), '.gammaloop_config.yaml'),
]

# Useful for pyLance to work fine when doing something like sum([[1,2],[4,5]],[]) to flatten the nested list
EMPTY_LIST: list[Any] = []


def get_git_revision() -> str | None:
    try:
        rev = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
        if int(subprocess.check_output('git status --porcelain -uno | wc -l', shell=True).decode('ascii').strip()) > 0:
            dirty_flag = '-modified'
        else:
            dirty_flag = ''
        return f"{rev}{dirty_flag}"
    except:
        return None


GIT_REVISION = get_git_revision()


class GammaLoopWarning(Enum):
    FloatInExpression = 100
    DroppingEpsilonTerms = 101

    def __str__(self):
        if self == GammaLoopWarning.FloatInExpression:
            return "FloatInExpression"
        elif self == GammaLoopWarning.DroppingEpsilonTerms:
            return "DroppingEpsilonTerms"
        else:
            raise GammaLoopError(f"Unknown side: {self}")


GL_WARNINGS_ISSUED: set[GammaLoopWarning] = set()

logger = logging.getLogger('GammaLoop')

try:
    import symbolica  # type: ignore # pylint: disable=import-error
    # pylint: disable=unused-import
    from symbolica import Expression, Transformer, set_license_key  # type: ignore
    if not str(os.path.abspath(symbolica.__file__)).startswith(str(os.path.join(GL_PATH, 'dependencies'))):
        logger.warning(
            """\n\n/*
|
| Symbolica is not being imported from the virtual environment of gammaloop.
| Instead it was loaded from %s%s%s.
| This may lead to unexpected behaviour. You can avoid this by running the following command:
|
| %sgammaloop --build_dependencies with_venv && source %s%s
|
\\*
""",
            Colour.RED, os.path.abspath(symbolica.__file__), Colour.END,
            Colour.GREEN, os.path.abspath(os.path.join(GL_PATH, 'dependencies', 'venv', 'bin', 'activate')), Colour.END)
except ImportError:
    logger.critical(
        "Could not import Symbolica, please run\n\n%sgammaloop --build_dependencies with_venv && source %s%s\n", Colour.GREEN,
        os.path.abspath(os.path.join(
            GL_PATH, 'dependencies', 'venv', 'bin', 'activate')),
        Colour.END)
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
    if gl_is_symbolica_registered is not None or (hasattr(symbolica, 'is_licensed') and symbolica.is_licensed()):  # type: ignore # nopep8
        gl_is_symbolica_registered = True
        return True

    if 'SYMBOLICA_LICENSE' not in os.environ or os.environ['SYMBOLICA_LICENSE'] == '':
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
            except Exception:  # pylint: disable=bare-except
                logger.critical(
                    "Could not set Symbolica license key from gammaloop configuration file specifying license '%s'.", symbolica_license)
                gl_is_symbolica_registered = False
                return False
        else:
            try:
                os.environ['SYMBOLICA_LICENSE'] = 'GAMMALOOP_USER'
                set_license_key('GAMMALOOP_USER')
                gl_is_symbolica_registered = True
                return True
            except:  # pylint: disable=bare-except
                logger.critical(
                    "The version of symbolica loaded requires a license. Specify via the environment variable 'SYMBOLICA_LICENSE' or through the config files of gammaloop.")
                gl_is_symbolica_registered = False
                return False
    else:
        try:
            set_license_key(os.environ['SYMBOLICA_LICENSE'])
            gl_is_symbolica_registered = True
            return True
        except:  # pylint: disable=bare-except
            if hasattr(symbolica, 'is_licensed') and symbolica.is_licensed():
                gl_is_symbolica_registered = True
                return True
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
