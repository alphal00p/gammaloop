#!/usr/bin/env python3
import sys
from pathlib import Path
import importlib
from argparse import ArgumentParser
import subprocess
import os
import yaml
from pprint import pformat
from gammaloop.misc.common import GammaLoopError, logger, Side, pjoin
from gammaloop.misc.utils import Colour
import gammaloop.base_objects.model as model
import gammaloop.base_objects.graph as graph
import gammaloop.cross_section.cross_section as cross_section
import gammaloop.cross_section.supergraph as supergraph
from gammaloop.exporters.exporters import AmplitudesExporter, CrossSectionsExporter, OutputMetaData
# This is the pyo3 binding of the gammaloop rust engine
import gammaloop._gammaloop as gl_rust  # pylint: disable=import-error, no-name-in-module

# pylint: disable=unused-variable

AVAILABLE_COMMANDS = [
    'import_model',
    'export_model',
    'import_graphs',
    'output',
    'help',
    'launch',
    'info',
    'integrate',
    'inspect',
    'test_ir_limits',
    'test_uv_limits',
]


def split_str_args(str_args) -> list[str]:
    return str_args.split(' ') if str_args != '' else []


class CommandList(list):

    @staticmethod
    def from_file(filename: str):
        command_file = CommandList()
        with open(filename, 'r', encoding='utf-8') as file:
            command_file.parse(file.read())
        return command_file

    @staticmethod
    def from_string(cmd_str: str):
        command_file = CommandList()
        command_file.parse(cmd_str)
        return command_file

    def add_command(self, cmd: str) -> None:
        cmd = cmd.strip()
        if len(cmd) == 0:
            return
        if cmd[0] == '!':
            self.append(('shell_run', cmd[1:]))
            return
        cmd_split: list[str] = cmd.split(' ', 1)
        if len(cmd_split[0]) > 1 and cmd_split[0].startswith('#'):
            return
        if cmd_split[0] not in AVAILABLE_COMMANDS:
            raise GammaLoopError(f"Unknown command: {cmd_split[0]}")
        if len(cmd_split) >= 2:
            self.append((cmd_split[0], cmd_split[1]))
        else:
            self.append((cmd_split[0], ''))

    def parse(self, cmd_str: str) -> None:
        multi_line = ''
        for line in sum([t.split('\n') for t in cmd_str.strip().split(';')], []):
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

    def nice_string(self) -> str:
        return '\n'.join([
            (f'{cmd} {args}' if cmd != 'shell_run' else f'!{args}').strip() for cmd, args in self
        ])


class GammaLoop(object):

    def __init__(self):
        self.model: model.Model = model.Model('NotLoaded')

        # Initialize a gammaloop rust engine worker which will be used throughout the session
        self.rust_worker: gl_rust.Worker = gl_rust.Worker.new()
        logger.debug(
            "Successfully initialized GammaLoop rust worker from library %s.", gl_rust.__file__)

        self.cross_sections: cross_section.CrossSectionList = cross_section.CrossSectionList()
        self.amplitudes: cross_section.AmplitudeList = cross_section.AmplitudeList()

        self.launched_output: Path | None = None
        self.command_history: CommandList = CommandList()

    def get_model_from_rust_worker(self) -> model.Model:
        return model.Model.from_yaml(self.rust_worker.get_model())

    def get_cross_sections_from_rust_worker(self) -> cross_section.CrossSectionList:
        if self.model.is_empty():
            raise GammaLoopError(
                "No model loaded. Please load a model first with 'import_model' command.")
        cross_sections = cross_section.CrossSectionList()
        yaml_cross_sections = yaml.safe_load(
            self.rust_worker.get_cross_sections())
        for yaml_cross_section in yaml_cross_sections:
            cross_sections.add_cross_section(cross_section.CrossSection.from_serializable_dict(
                self.model, yaml_cross_section))
        return cross_sections

    def get_amplitudes_from_rust_worker(self) -> cross_section.AmplitudeList:
        if self.model.is_empty():
            raise GammaLoopError(
                "No model loaded. Please load a model first with 'import_model' command.")
        amplitudes = cross_section.AmplitudeList()
        yaml_amplitudes = yaml.safe_load(
            self.rust_worker.get_amplitudes())
        for yaml_amplitude in yaml_amplitudes:
            amplitudes.add_amplitude(cross_section.Amplitude.from_serializable_dict(
                self.model, yaml_amplitude))
        return amplitudes

    def run(self, cmds: CommandList) -> None:

        for cmd, args in cmds:
            logger.debug("Running command '%s' with arguments '%s'", cmd, args)
            self.command_history.append((cmd, args))
            match cmd:
                case 'shell_run':  # Triggered by special character '!' at the beginning of the command
                    subprocess.run(args, shell=True,
                                   cwd=os.getcwd(), check=True)
                case _:
                    if cmd in AVAILABLE_COMMANDS:
                        if not hasattr(self, f'do_{cmd}'):
                            raise GammaLoopError(
                                f"Command '{cmd}' not yet implemented.")
                        getattr(self, f'do_{cmd}')(args)
                    else:
                        raise GammaLoopError(f"Invalid command '{cmd}'")

    # import_model command
    import_model_parser = ArgumentParser(prog='import_model')
    import_model_parser.add_argument('model', metavar='model', type=str,
                                     help='UFO model name or or YAML file path to load')
    import_model_parser.add_argument('--format', '-f', type=str, default='ufo',
                                     choices=['yaml', 'ufo'], help='Format of the model for load.')

    def do_import_model(self, str_args: str) -> None:
        if str_args == 'help':
            self.import_model_parser.print_help()
            return
        args = self.import_model_parser.parse_args(split_str_args(str_args))

        match args.format:
            case 'ufo':
                self.model = model.Model.from_ufo_model(args.model)
            case 'yaml':
                with open(args.model, 'r', encoding='utf-8') as file:
                    self.model = model.Model.from_yaml(file.read())

        # Assign the UFO model to our rust worker
        self.rust_worker.load_model_from_yaml_str(self.model.to_yaml())

        logger.info("Successfully loaded model '%s'.", self.model.name)

    # export_model command
    export_model_parser = ArgumentParser(prog='export_model')
    export_model_parser.add_argument('model_file_path', metavar='model_file_path', type=str,
                                     help='Model file path to export to')
    export_model_parser.add_argument('--format', '-f', type=str, default='yaml',
                                     choices=['yaml', 'ufo'], help='Format to export the model in.')

    def do_export_model(self, str_args: str) -> None:
        if str_args == 'help':
            self.export_model_parser.print_help()
            return
        args = self.export_model_parser.parse_args(split_str_args(str_args))

        if self.model.is_empty():
            raise GammaLoopError(
                "No model loaded. Please load a model first with 'import_model' command.")

        match args.format:
            case 'ufo':
                raise GammaLoopError(
                    "Exporting to UFO format is not yet supported.")
            case 'yaml':
                with open(args.model_file_path, 'w', encoding='utf-8') as file:
                    file.write(self.model.to_yaml())

        logger.info("Successfully exported model '%s' to '%s'.",
                    self.model.name, args.model_file_path)

    # import_graphs command
    import_graphs_parser = ArgumentParser(prog='import_graphs')
    import_graphs_parser.add_argument('file_path', metavar='file_path', type=str,
                                      help='Path to the qgraph python output to load')
    import_graphs_parser.add_argument('--no_compile', '-nc', action='store_true',
                                      default=False, help='Prevent compilation of qgraph python output.')
    import_graphs_parser.add_argument('--format', '-f', type=str, default='yaml',
                                      choices=['yaml', 'qgraph'], help='Format to import the graphs in.')

    def do_import_graphs(self, str_args: str) -> None:
        if str_args == 'help':
            self.import_graphs_parser.print_help()
            return
        args = self.import_graphs_parser.parse_args(split_str_args(str_args))

        if self.model.is_empty():
            raise GammaLoopError(
                "No model loaded. Please load a model first with 'import_model' command.")

        if not os.path.isfile(args.file_path):
            raise GammaLoopError(f"File '{args.file_path}' does not exist")

        file_path = Path(os.path.abspath(args.file_path))
        sys.path.insert(0, str(file_path.parent))

        if not args.no_compile:
            # compile the file first before importing it with optimization flag and same executable as used for gammaLoop.
            # This avoids that memory isn't freed after compiling when using the __import__ directly
            logger.info("Compiling imported supergraphs.")
            subprocess.run([sys.executable, '-O', '-m',
                           file_path.stem], cwd=file_path.parent, check=True)

        qgraph_loaded_module = __import__(file_path.stem)
        # Reload to avoid border effects if this is the second time qgraph output is loaded within this same Python session.
        importlib.reload(qgraph_loaded_module)
        logger.info("Imported %s graphs from qgraph output '%s'.",
                    len(qgraph_loaded_module.graphs), args.file_path)
        del sys.path[0]

        graphs = []
        for i_qg, qgraph_object in enumerate(qgraph_loaded_module.graphs):
            graphs.append(graph.Graph.from_qgraph(
                self.model, qgraph_object, name=f"{file_path.stem}_{i_qg}"))
        logger.info("Successfully loaded %s graphs.",
                    len(qgraph_loaded_module.graphs))

        # Now determine if it is a supergraph or an amplitude graph
        graph_type = None
        for a_graph in graphs:
            if len(a_graph.get_incoming_edges()) == 0 and len(a_graph.get_outgoing_edges()) == 0:
                if graph_type is None:
                    graph_type = 'supergraph'
                elif graph_type != 'supergraph':
                    raise GammaLoopError(
                        "Mixed type of graphs (supergraph, forward scattering graph and/or amplitude graph) found in qgraph output loaded.")
            elif len(a_graph.get_incoming_edges()) == len(a_graph.get_outgoing_edges()) == len(a_graph.external_connections) and not any(None in c for c in a_graph.external_connections):
                if graph_type is None:
                    graph_type = 'forward_scattering_graph'
                elif graph_type != 'forward_scattering_graph':
                    raise GammaLoopError(
                        "Mixed type of graphs (supergraph, forward scattering graph and/or amplitude graph) found in qgraph output loaded.")
            else:
                if graph_type is None:
                    graph_type = 'amplitude'
                elif graph_type != 'amplitude':
                    raise GammaLoopError(
                        "Mixed type of graphs (supergraph, forward scattering graph and/or amplitude graph) found in qgraph output loaded.")

        if graph_type == 'supergraph':
            raise GammaLoopError(
                "Processing of supergraphs not yet supported.")

        if graph_type == 'forward_scattering_graph':
            logger.warning(
                "Processing of forward scattering graphs not fully implemented yet.")
            self.cross_sections = cross_section.CrossSectionList([cross_section.CrossSection(
                f'{file_path.stem}',
                # Wrap the forward scattering graphs within a dummy supergraph
                [cross_section.supergraph.SuperGraph(
                    sg_id=i, graph=graph.Graph.empty_graph('DUMMY'), multiplicity=1.0,
                    topology_class=[],
                    cuts=[
                        cross_section.supergraph.SuperGraphCut(
                            cut_edges=[],
                            forward_scattering_graph=cross_section.supergraph.ForwardScatteringGraph(
                                sg_id=i,
                                sg_cut_id=0,
                                multiplicity=1.0,
                                graph=g,
                                cuts=[]  # Will be filled in later
                            )
                        )]
                ) for i, g in enumerate(graphs)]
            )])

        if graph_type == 'amplitude':
            self.amplitudes = cross_section.AmplitudeList([cross_section.Amplitude(
                file_path.stem,
                [
                    supergraph.AmplitudeGraph(
                        sg_id=0, sg_cut_id=0, fs_cut_id=i, amplitude_side=Side.LEFT,
                        graph=g
                    )
                ]
            ) for i, g in enumerate(graphs)])

    # output command
    output_parser = ArgumentParser(prog='output')
    output_parser.add_argument(
        'output_path', type=str, help='Path to output the cross section to')

    def do_output(self, str_args: str) -> None:
        if str_args == 'help':
            self.output_parser.print_help()
            return
        args = self.output_parser.parse_args(split_str_args(str_args))

        if Path(args.output_path).exists():
            raise GammaLoopError(
                f"Output path '{args.output_path}' already exists, clean it up first.")

        if self.model.is_empty():
            raise GammaLoopError(
                "No model loaded. Please load a model first with 'import_model' command.")

        if len(self.cross_sections) == 0 and len(self.amplitudes) == 0:
            raise GammaLoopError("No process generated yet.")

        if len(self.cross_sections) > 0:
            cross_section_exporter = CrossSectionsExporter(self, args)
            cross_section_exporter.export(
                args.output_path, self.cross_sections)

        if len(self.amplitudes) > 0:
            amplitude_exporter = AmplitudesExporter(self, args)
            amplitude_exporter.export(args.output_path, self.amplitudes)

    #
    # Run interface type of commands below (those bound to a particular output already generated)
    #

    # launch command
    launch_parser = ArgumentParser(prog='launch')
    launch_parser.add_argument(
        'path_to_launch', type=str, help='Path to launch a given run command to')

    def do_launch(self, str_args: str) -> None:
        if str_args == 'help':
            self.launch_parser.print_help()
            return
        args = self.launch_parser.parse_args(split_str_args(str_args))

        if not Path(pjoin(args.path_to_launch, 'output_metadata.yaml')).exists():
            raise GammaLoopError(
                f"Cross-section or amplitude path '{args.path_to_launch}' does not appear valid.")

        with open(pjoin(args.path_to_launch, 'output_metadata.yaml'), 'r', encoding='utf-8') as file:
            output_metadata = OutputMetaData.from_yaml_str(file.read())

        # Sync the model with the one used in the output
        with open(pjoin(args.path_to_launch, 'sources', 'model', f"{output_metadata['model_name']}.yaml"), 'r', encoding='utf-8') as file:
            yaml_model = file.read()
            self.model = model.Model.from_yaml(yaml_model)
            self.rust_worker.load_model_from_yaml_str(yaml_model)

        # Depending on the type of output, sync cross-section or amplitude
        if output_metadata['output_type'] == 'amplitudes':
            self.amplitudes = cross_section.AmplitudeList()
            self.rust_worker.reset_amplitudes()
            for amplitude_name in output_metadata['contents']:
                with open(pjoin(args.path_to_launch, 'sources', 'amplitudes', f'{amplitude_name}', 'amplitude.yaml'), 'r', encoding='utf-8') as file:
                    amplitude_yaml = file.read()
                    self.amplitudes.add_amplitude(
                        cross_section.Amplitude.from_yaml_str(self.model, amplitude_yaml))
                    self.rust_worker.add_amplitude_from_yaml_str(
                        amplitude_yaml)

        if output_metadata['output_type'] == 'cross_sections':
            self.cross_sections = cross_section.CrossSectionList()
            self.rust_worker.reset_cross_sections()
            for cross_section_name in output_metadata['contents']:
                with open(pjoin(args.path_to_launch, 'sources', 'cross_sections', f'{cross_section_name}', 'cross_section.yaml'), 'r', encoding='utf-8') as file:
                    cross_section_yaml = file.read()
                    self.cross_sections.add_cross_section(
                        cross_section.CrossSection.from_yaml_str(self.model, cross_section_yaml))
                    self.rust_worker.add_cross_section_from_yaml_str(
                        cross_section_yaml)

        self.launched_output = Path(os.path.abspath(args.path_to_launch))

    # info command
    info_parser = ArgumentParser(prog='info')
    info_parser.add_argument('--full', '-f', action='store_true',
                             default=False, help='Show detailed output.')

    def do_info(self, str_args: str) -> None:
        if str_args == 'help':
            self.info_parser.print_help()
            return
        args = self.info_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        # To improve and completement with process definition when available, test results etc...
        if len(self.amplitudes) > 0:
            n_tot_amplitude_graphs = sum(
                [len(amp.amplitude_graphs) for amp in self.amplitudes])
            logger.info("Output '%s' contains %s amplitudes with a total of %s forward scattering graphs.",
                        self.launched_output.stem, len(self.amplitudes), n_tot_amplitude_graphs)
            if args.full:
                logger.info("Amplitudes:\n%s",
                            pformat(self.amplitudes.to_serializable()))
        if len(self.cross_sections) > 0:
            n_tot_supergraphs = sum([len(cs.supergraphs)
                                    for cs in self.cross_sections])
            n_tot_forward_scattering_graphs = sum(
                [len(sg.cuts) for cs in self.cross_sections for sg in cs.supergraphs])
            logger.info("Output '%s' contains %s cross-sections with a total of %s supergraphs and %s forward scattering graphs.",
                        self.launched_output.stem, len(self.cross_sections), n_tot_supergraphs, n_tot_forward_scattering_graphs)
            if args.full:
                logger.info("Cross-sections:\n%s",
                            pformat(self.cross_sections.to_serializable()))

    # inspect command
    inspect_parser = ArgumentParser(prog='inspect')

    def do_inspect(self, str_args: str) -> None:
        if str_args == 'help':
            self.inspect_parser.print_help()
            return
        args = self.inspect_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        raise GammaLoopError("Command not implemented yet")

    # integrate command
    integrate_parser = ArgumentParser(prog='integrate')

    def do_integrate(self, str_args: str) -> None:
        if str_args == 'help':
            self.integrate_parser.print_help()
            return
        args = self.integrate_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        raise GammaLoopError("Command not implemented yet")

    # test_ir_limits
    test_ir_limits_parser = ArgumentParser(prog='test_ir_limits')

    def do_test_ir_limits(self, str_args: str) -> None:
        if str_args == 'help':
            self.test_ir_limits_parser.print_help()
            return
        args = self.test_ir_limits_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        raise GammaLoopError("Command not implemented yet")

    # test_uv_limits
    test_uv_limits_parser = ArgumentParser(prog='test_ir_limits')

    def do_test_uv_limits(self, str_args: str) -> None:
        if str_args == 'help':
            self.test_uv_limits_parser.print_help()
            return
        args = self.test_uv_limits_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        raise GammaLoopError("Command not implemented yet")

    # help command
    help_parser = ArgumentParser(prog='help')
    help_parser.add_argument('cmd', metavar='cmd', type=str,
                             help='Specify command to print help for')

    def do_help(self, str_args: str) -> None:
        if str_args == '':
            print(
                f"Available commands are: {', '.join(f'{Colour.GREEN}{a}{Colour.END}' for a in AVAILABLE_COMMANDS)}")
            return

        if str_args == 'help':
            self.help_parser.print_help()
            return

        args = self.help_parser.parse_args(split_str_args(str_args))
        self.run(CommandList.from_string(f"{args.cmd} help"))
