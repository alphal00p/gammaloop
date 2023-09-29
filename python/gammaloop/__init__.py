import os
from .interface.gammaloop_interface import CommandList
from .interface.gammaloop_interface import GammaLoop
from .misc.common import GL_DEBUG, GL_CONSOLE_HANDLER, register_symbolica, logger, DATA_PATH, GammaLoopError

__version__ = "0.0.1"


def cli():
    from argparse import ArgumentParser  # pylint: disable=import-outside-toplevel
    import logging  # pylint: disable=import-outside-toplevel

    argparser = ArgumentParser(
        prog='γLoop',
        description='Compute differential cross-section with Local Unitarity')
    argparser.add_argument('command_file', type=str,
                           help='Input file listing commands to run')
    argparser.add_argument('--debug', '-d', action='store_true',
                           default=False, help='Enable debug mode')
    argparser.add_argument('--quiet', '-q', action='store_true',
                           default=False, help='Enable quiet mode')
    argparser.add_argument('--no_file_handler', '-nfh', action='store_true',
                           default=False, help='Disable logging to files')
    argparser.add_argument('--command', '-c', dest='run_command', action='store_true',
                           default=False, help='Directly run command specified')
    args = argparser.parse_args()

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
        GL_DEBUG = True  # pylint: disable=redefined-outer-name, unused-variable, invalid-name
        GL_CONSOLE_HANDLER.setLevel(logging.DEBUG)

    if args.quiet:
        GL_DEBUG = False  # pylint: disable=invalid-name
        GL_CONSOLE_HANDLER.setLevel(logging.CRITICAL)

    if args.run_command:
        commands = CommandList.from_string(args.command_file)
    else:
        if not os.path.isfile(args.command_file):
            command_file = os.path.join(
                DATA_PATH, 'run_cards', args.command_file)
            if not os.path.isfile(command_file):
                raise GammaLoopError(
                    f"File '{args.command_file}' does not exist.")
        else:
            command_file = args.command_file
        commands = CommandList.from_file(command_file)

    gamma_loop.run(commands)
