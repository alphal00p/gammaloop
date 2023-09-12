#!/usr/bin/env python3
import argparse
import logging
import sys
from gammaloop.misc.common import *
from gammaloop.misc.utils import Colour
from argparse import ArgumentParser
import gammaloop.base_objects.model as model
import gammaloop.base_objects.graph as graph
from pathlib import Path
import importlib
from typing import Sequence
# This is the pyo3 binding of the gammaloop rust engine
import gammaloop._gammaloop as gl_rust

AVAILABLE_COMMANDS = [
    'import_model',
    'export_model',
    'import_qgraph',
    'help'
]

class CommandList(list):

    @staticmethod
    def from_file(filename: str):
        cf = CommandList()
        with open(filename, 'r') as f:
            cf.parse(f.read())
        return cf

    @staticmethod
    def from_string(cmd_str: str):
        cf = CommandList()
        cf.parse(cmd_str)
        return cf

    def add_command(self, cmd: str) -> None:
        cmd_split: list[str] = cmd.split(' ', 1)
        if len(cmd_split[0])>1 and cmd_split[0].startswith('#'):
            return
        if cmd_split[0] not in AVAILABLE_COMMANDS:
            raise GammaLoopError("Unknown command: {}".format(cmd_split[0]))
        if len(cmd_split) >= 2:
            self.append((cmd_split[0], cmd_split[1]))
        else:
            self.append((cmd_split[0], ''))

    def parse(self, cmd_str: str) -> None:
        multi_line = ''
        for line in sum([ t.split('\n') for t in cmd_str.strip().split(';') ],[]):
            line = line.strip()
            if line == '':
                continue
            if line.endswith('\\\n'):
                multi_line += line[:-2]
                continue
            elif line.endswith('\\'):
                multi_line += line[:-1]
                continue

            if line.endswith('\n'):
                multi_line += line[:-1]
            else:
                multi_line += line
            if multi_line != '':
                self.add_command(multi_line)
                multi_line = ''
        if multi_line != '':
            self.add_command(multi_line)


class GammaLoop(object):

    def __init__(self):
        self.model : model.Model = model.Model('NotLoaded')

        # Initialize a gammaloop rust engine worker which will be used throughout the session
        self.rust_worker : gl_rust.Worker = gl_rust.Worker.new()
        logger.debug(f"Successfully initialized GammaLoop rust worker from library {gl_rust.__file__}.")

    def run(self, cmds: CommandList) -> None:

        for cmd, args in cmds:
            logger.debug(f"Running command '{cmd}' with arguments '{args}'")
            match cmd:
                case 'import_model':
                    self.do_import_model(args)
                case 'export_model':
                    self.do_export_model(args)
                case 'import_qgraph':
                    self.do_import_qgraph(args)
                case 'help':
                    self.do_help(args)
                case _:
                    raise GammaLoopError(f"Invalid command '{cmd}'")

    # import_model command
    import_model_parser = ArgumentParser(prog='import_model')
    import_model_parser.add_argument('model', metavar='model', type=str,
                    help='UFO model name or or YAML file path to load')
    import_model_parser.add_argument('--format', '-f', type=str, default='ufo',
                    choices=['yaml', 'ufo'], help='Format of the model for load.')
    def do_import_model(self, str_args: str) -> None:
        if str_args=='help':
            self.import_model_parser.print_help()
            return
        args = self.import_model_parser.parse_args(str_args.split(' '))

        match args.format:
            case 'ufo':
                self.model = model.Model.from_UFOModel(args.model)
            case 'yaml':
                with open(args.model,'r') as f:
                    self.model = model.Model.from_yaml(f.read())

        # Assign the UFO model to our rust worker
        self.rust_worker.load_model_from_yaml_str(self.model.to_yaml())

        logger.info("Successfully loaded model '{}'.".format(self.model.name))

    # export_model command
    export_model_parser = ArgumentParser(prog='export_model')
    export_model_parser.add_argument('model_file_path', metavar='model_file_path', type=str,
                    help='Model file path to export to')
    export_model_parser.add_argument('--format', '-f', type=str, default='yaml',
                    choices=['yaml', 'ufo'], help='Format to export the model in.')
    def do_export_model(self, str_args: str) -> None:
        if str_args=='help':
            self.export_model_parser.print_help()
            return
        args = self.export_model_parser.parse_args(str_args.split(' '))
        
        if self.model.name == 'NotLoaded':
            raise GammaLoopError("No model loaded. Please load a model first with 'import_model' command.")

        match args.format:
            case 'ufo':
                raise GammaLoopError("Exporting to UFO format is not yet supported.")
            case 'yaml':
                with open(args.model_file_path,'w') as f:
                    f.write(self.model.to_yaml())

        logger.info(f"Successfully exported model '{self.model.name}' to '{args.model_file_path}'.")

    # import_qgraph command
    import_qgraf_parser = ArgumentParser(prog='import_qgraph')
    import_qgraf_parser.add_argument('qgraph_file', metavar='qgraph_file', type=str,
                    help='Path to the qgraph python output to load')
    import_qgraf_parser.add_argument('--no_compile', '-nc', action='store_true',
                    default=False, help='Prevent compilation of qgraph python output.')
    def do_import_qgraph(self, str_args: str) -> None:
        if str_args=='help':
            self.import_qgraf_parser.print_help()
            return
        args = self.import_qgraf_parser.parse_args(str_args.split(' '))

        if self.model.name == 'NotLoaded':
            raise GammaLoopError("No model loaded. Please load a model first with 'import_model' command.")

        if not os.path.isfile(args.qgraph_file):
            raise GammaLoopError(f"File '{args.qgraph_file}' does not exist")

        qgraph_path = Path(args.qgraph_file)
        sys.path.insert(0, str(qgraph_path.parent))

        if not args.no_compile:
            # compile the file first before importing it with optimization flag and same executable as used for gammaLoop.
            # This avoids that memory isn't freed after compiling when using the __import__ directly
            logger.info("Compiling imported supergraphs.")
            subprocess.run([sys.executable, '-O', '-m', p.stem], cwd=qgraph_path.parent)

        qgraph_loaded_module = __import__(qgraph_path.stem)
        # Reload to avoid border effects if this is the second time qgraph output is loaded within this same Python session.
        importlib.reload(qgraph_loaded_module)
        logger.info("Imported {} graphs from qgraph output '{}'.".format(len(qgraph_loaded_module.graphs),args.qgraph_file))
        del sys.path[0]

        graphs = []
        for qgraph_object in qgraph_loaded_module.graphs:
            graphs.append(graph.Graph.from_qgraph(self.model, qgraph_object))
        logger.info("Successfully loaded {} supergraphs.".format(len(qgraph_loaded_module.graphs)))
        
    # help command
    help_parser = ArgumentParser(prog='help')
    help_parser.add_argument('cmd', metavar='cmd', type=str,
                    help='Specify command to print help for')
    def do_help(self, str_args: str) -> None:
        if str_args == '':
            print("Available commands are: %s" % ', '.join(f'{Colour.GREEN}{a}{Colour.END}' for a in AVAILABLE_COMMANDS))
            return
        
        if str_args=='help':
            self.help_parser.print_help()
            return

        args = self.help_parser.parse_args(str_args.split(' '))
        self.run(CommandList.from_string(f"{args.cmd} help"))