import os
import sys
from subprocess import Popen
from enum import StrEnum

__version__ = "0.0.1"

GL_PATH = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

DEPENDENCIES_CHECKED = False


class CLIColour(StrEnum):
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def check_gammaloop_dependencies(clean_dependencies=False, build_dependencies=False, no_gammaloop_python_venv=False):

    global DEPENDENCIES_CHECKED
    if DEPENDENCIES_CHECKED and all(opt is False for opt in [clean_dependencies, build_dependencies, no_gammaloop_python_venv]):
        return

    gammaloop_root_path = os.path.abspath(GL_PATH)

    if clean_dependencies:
        print("Cleaning gammaloop dependencies in '%s' using %s./bin/build_dependencies.sh clean%s ..." % (
            os.path.abspath(os.path.join(gammaloop_root_path, 'dependencies')), CLIColour.GREEN, CLIColour.END))
        Popen([f"bash {os.path.join('.', 'bin','build_dependencies.sh')} clean",],
              shell=True, cwd=gammaloop_root_path).wait()
        print("%sDependencies folder cleaned successfully.%s" %
              (CLIColour.GREEN, CLIColour.END))
        if not build_dependencies:
            sys.exit(0)

    if build_dependencies:
        print("Building gammaloop dependencies in '%s' using the %s./bin/build_dependencies.sh%s script ..." % (
            os.path.abspath(os.path.join(gammaloop_root_path, 'dependencies')), CLIColour.GREEN, CLIColour.END))
        Popen([f"bash {os.path.join('.', 'bin', 'build_dependencies.sh')}",],
              shell=True, cwd=gammaloop_root_path).wait()
        if not os.path.isfile(os.path.join(gammaloop_root_path, 'dependencies', 'INSTALLED')):
            print("%sCould not build the dependencies. Find more information in '%s'.%s" % (
                CLIColour.RED, os.path.join(gammaloop_root_path, 'dependencies', 'dependency_build.log.log'), CLIColour.END))
            sys.exit(1)
        else:
            print("%sDependencies built successfully.%s" %
                  (CLIColour.GREEN, CLIColour.END))
            sys.exit(0)

    if not os.path.isfile(os.path.join(gammaloop_root_path, 'dependencies', 'INSTALLED')):
        print("\nGammaloop dependencies are %snot installed%s. Run '%sgammaloop --build_dependencies%s' to install them first. Exiting.\n" % (
            CLIColour.RED, CLIColour.END, CLIColour.GREEN, CLIColour.END))
        sys.exit(1)

    venv_path = os.path.abspath(os.path.join(
        gammaloop_root_path, 'dependencies', 'venv'))

    if venv_path != os.path.abspath(os.path.join(os.path.dirname(sys.executable), os.path.pardir)):
        # Do not warn about this for now as it is not strictly necessary as of now
        # print("%sWARNING:%s It is recommended to run gammaloop within its Python virtual environment, by issuing the command:\n\n%ssource %s/bin/activate%s\n" % (
        #    CLIColour.YELLOW, CLIColour.END, CLIColour.GREEN, venv_path, CLIColour.END))
        pass

    if not no_gammaloop_python_venv:
        if not os.path.isfile(os.path.join(venv_path, 'site_paths.txt')):
            print("%sWARNING:%s Could not find Python virtual environment list of sites in file '%s'.\nConsider running '%sgammaloop --clean_dependencies --build_dependencies%s' to re-install gammaloop dependencies.\n" % (
                CLIColour.YELLOW, CLIColour.END, os.path.join(venv_path, 'venv_site_paths.txt'), CLIColour.GREEN, CLIColour.END))
        else:
            try:
                with open(os.path.join(venv_path, 'site_paths.txt'), 'r') as f:
                    site_paths = eval(f.read())
            except Exception as e:
                site_paths = None
                print("%sWARNING:%s Could not extract list of site paths of the Python virtual environment of gammaloop from file '%s' (error: %s%s%s).\nConsider running '%sgammaloop --clean_dependencies --build_dependencies%s' to re-install gammaloop dependencies.\n" % (
                    CLIColour.YELLOW, CLIColour.END, os.path.join(venv_path, 'site_paths.txt'), CLIColour.RED, str(e), CLIColour.END, CLIColour.GREEN, CLIColour.END))
            if site_paths is not None:
                site_paths = [sp for sp in site_paths if sp.startswith(
                    venv_path) and (len(sys.path) == 0 or sp != sys.path[0])]
                if len(site_paths) > 0:
                    print("%sINFO:%s The following paths have been automatically and temporarily added by gammaloop to your PYTHONPATH:\n%s" % (
                        CLIColour.GREEN, CLIColour.END, '\n'.join('%s%s%s' % (CLIColour.BLUE, site_path, CLIColour.END) for site_path in site_paths)))
                    for site_path in site_paths:
                        sys.path.insert(0, site_path)

    DEPENDENCIES_CHECKED = True


def cli():
    from argparse import ArgumentParser  # pylint: disable=import-outside-toplevel
    import logging  # pylint: disable=import-outside-toplevel

    argparser = ArgumentParser(
        prog='γLoop',
        description='Compute differential cross-section with Local Unitarity')
    argparser.add_argument('command_file', type=str,
                           help='Input file listing commands to run', nargs='?')
    argparser.add_argument('--debug', '-d', action='store_true',
                           default=False, help='Enable debug mode')
    argparser.add_argument('--build_dependencies', '-bd', action='store_true',
                           default=False, help='Build dependencies')
    argparser.add_argument('--clean_dependencies', '-cd', action='store_true',
                           default=False, help='Clean dependencies build')
    argparser.add_argument('--no_gammaloop_python_venv', '-no_venv', action='store_true',
                           default=False, help='Do not add python virtual environment installed by gammaloop to the PYTHONPATH')
    argparser.add_argument('--quiet', '-q', action='store_true',
                           default=False, help='Enable quiet mode')
    argparser.add_argument('--no_file_handler', '-nfh', action='store_true',
                           default=False, help='Disable logging to files')
    argparser.add_argument('--show_venv_activate_path', '-venv', action='store_true',
                           default=False, help='Show gammaloop Python venv activate path')
    argparser.add_argument('--command', '-c', dest='run_command', action='store_true',
                           default=False, help='Directly run command specified')
    argparser.add_argument('--logging_format', '-lf', dest='logging_format', default='long',
                           choices=['long', 'short', 'min', 'none'], help='Type of prefix for the logging format')
    args = argparser.parse_args()

    from . import misc
    misc.LOGGING_PREFIX_FORMAT = args.logging_format

    if args.show_venv_activate_path:
        venv_path = os.path.abspath(os.path.join(
            GL_PATH, 'dependencies', 'venv', 'bin', 'activate'))
        if not os.path.isfile(venv_path):
            print("%sCould not find the gammaloop Python virtual environment activate script at '%s'.%s" % (
                CLIColour.RED, venv_path, CLIColour.END))
            print("Make sur to first run %sgammaloop --build_dependencies%s%s" %
                  (CLIColour.GREEN, venv_path, CLIColour.END))
            sys.exit(1)
        else:
            print(venv_path)
            sys.exit(0)

    # Before importing anything of gammaloop, check the dependencies with proper arguments
    check_gammaloop_dependencies(args.clean_dependencies, args.build_dependencies,
                                 args.no_gammaloop_python_venv)

    from .interface.gammaloop_interface import CommandList
    from .interface.gammaloop_interface import GammaLoop
    import gammaloop.misc.common as common
    from .misc.common import GL_CONSOLE_HANDLER, register_symbolica, logger, DATA_PATH, GammaLoopError

    if args.command_file is None:
        logger.critical(
            "%sCommand file or string not specified.%s", CLIColour.RED, CLIColour.END)
        argparser.print_help()
        sys.exit(1)

    if args.no_file_handler:
        os.environ['GL_ENABLE_FILE_HANDLERS'] = 'FALSE'
    else:
        os.environ['GL_ENABLE_FILE_HANDLERS'] = 'TRUE'

    success = register_symbolica()
    if not success:
        logger.warning("Could not register Symbolica license key; you will be facing a prompting banner from this dependency.\n" +
                       "If you want to avoid this, please set the environment variable 'SYMBOLICA_LICENSE' to your Symbolica license key or paste it in the file ")

    gamma_loop = GammaLoop()

    if args.debug:
        common.GL_DEBUG = True
        GL_CONSOLE_HANDLER.setLevel(logging.DEBUG)

    if args.quiet:
        common.GL_DEBUG = False
        GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)

    if args.run_command:
        commands = CommandList.from_string(args.command_file)
    else:
        if not os.path.isfile(args.command_file):
            command_file = None
            paths_to_try = [
                os.path.abspath(os.path.join(
                    GL_PATH, args.command_file)),
                os.path.abspath(os.path.join(
                    GL_PATH, 'examples', 'cards', args.command_file)),
                os.path.abspath(os.path.join(
                    DATA_PATH, 'run_cards', args.command_file)),
            ]
            for path in paths_to_try:
                if os.path.isfile(path):
                    command_file = path
                    break
            if command_file is None:
                raise GammaLoopError(
                    f"File '{args.command_file}' does not exist.")
        elif args.command_file is not None:
            command_file = args.command_file
        else:
            raise GammaLoopError(
                "No command file or string specified.")
        commands = CommandList.from_file(command_file)

    gamma_loop.run(commands)
