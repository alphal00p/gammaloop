#!/usr/bin/env python3
import sys
from pathlib import Path
import importlib
from argparse import ArgumentParser
import subprocess
import os
from typing import Any
from pprint import pformat
import yaml
import shutil
import copy
from gammaloop.misc.common import GammaLoopError, logger, Side, pjoin, load_configuration, GAMMALOOP_CONFIG_PATHS, gl_is_symbolica_registered, GL_PATH, GL_WARNINGS_ISSUED, GammaLoopWarning
from gammaloop.misc.utils import Colour, verbose_yaml_dump
from gammaloop.base_objects.model import Model, InputParamCard
from gammaloop.base_objects.param_card import ParamCard, ParamCardWriter
from gammaloop.base_objects.graph import Graph
import gammaloop.cross_section.cross_section as cross_section
import gammaloop.cross_section.supergraph as supergraph
from gammaloop.exporters.exporters import AmplitudesExporter, CrossSectionsExporter, OutputMetaData, update_run_card_in_output
# This is the pyo3 binding of the gammaloop rust engine
import gammaloop._gammaloop as gl_rust  # pylint: disable=import-error, no-name-in-module # type: ignore

# pylint: disable=unused-variable

AVAILABLE_COMMANDS = [
    'import_model',
    'export_model',
    'import_graphs',
    'show_settings',
    'output',
    'help',
    'launch',
    'info',
    'integrate',
    'inspect',
    'hpc_run',
    'test_ir_limits',
    'test_uv_limits',
    'set',
    'reset'
]


class GammaLoopConfiguration(object):

    def __init__(self, path: str | None = None, quiet=False):
        self._shorthands: dict[str, str] = {}
        self._config: dict[str, Any] = {
            'symbolica': {
                'license': "GAMMALOOP_USER"
            },
            'drawing': {
                'mode': 'feynmp',
                'combined_graphs_pdf_grid_shape': [3, 2],
                'feynmp': {
                    'show_edge_labels': True,
                    'show_particle_names': True,
                    'show_edge_names': True,
                    'show_edge_momenta': True,
                    'show_edge_composite_momenta': True,
                    'label_size': "20pt",
                    'label_distance': 13.0,
                    'draw_arrow_for_all_edges': True,
                    'draw_lmb': True,
                    'arc_max_distance': 1.0,
                    'external_legs_tension': 3.0,
                    'default_tension': 1.0,
                    'show_vertex_labels': True,
                    'use_vertex_names': True,
                    'vertex_size': 5.0,
                    'vertex_shape': "circle",
                    'line_width': 1.75,
                    'arrow_size_for_single_line': 1.5,
                    'arrow_size_for_double_line': 2.0,
                    'line_color': "black",
                    'label_color': "black",
                    'non_lmb_color': "blue",
                    'lmb_color': "red",
                    'anchor_tension': 2.0,
                    'caption_size': "40pt",
                }
            },
            'export_settings': {
                'write_default_settings': False,
            },
            'run_settings': {
                'General': {
                    'debug': 0,
                    'use_ltd': False
                },
                'Integrand': {
                    'type': 'gamma_loop'
                },
                'Kinematics': {
                    'e_cm': 3.0,
                    'externals': {
                        'type': 'constant',
                        'momenta': None  # Will be set automatically
                    }
                },
                'Parameterization': {
                    'mode': 'spherical',
                    'mapping': 'linear',
                    'b': 10.0
                },
                'Integrator': {
                    'n_bins': 16,
                    'bin_number_evolution': None,
                    'min_samples_for_update': 100,
                    'n_start': 1000000,
                    'n_increase': 0,
                    'n_max': 1000000000,
                    'integrated_phase': 'real',
                    'learning_rate': 1.5,
                    'train_on_avg': False,
                    'show_max_wgt_info': False,
                    'max_prob_ratio': 0.01,
                    'seed': 0
                },
                'Observables': [],
                'Selectors': [],
                'Stability': {
                    'rotation_axis': 'x',
                    'levels': [
                        {
                            'precision': 'Double',
                            'required_precision_for_re': 1.e-5,
                            'required_precision_for_im': 1.e-5,
                            'escalate_for_large_weight_threshold': 0.9
                        },
                        {
                            'precision': 'Quad',
                            'required_precision_for_re': 1.e-5,
                            'required_precision_for_im': 1.e-5,
                            'escalate_for_large_weight_threshold': -1.0
                        }
                    ]
                },
                'sampling': {
                    'type': 'default'
                }
            }
        }

        if path is None:
            for config_path in GAMMALOOP_CONFIG_PATHS:
                if os.path.exists(config_path):
                    self.update(load_configuration(config_path, quiet), '')
        # Allowing for '' to be provided as a path in order to force the default configuration to be kept.
        elif path != '':
            self.update(load_configuration(path, quiet))

        self.update_shorthands()

    def update_shorthands(self, cur_path=None, root_dict: dict[str, Any] | None = None) -> None:
        if cur_path is None:
            cur_path = []
        if root_dict is None:
            root_dict = self._config
        for key, value in root_dict.items():
            if not isinstance(key, str):
                continue
            # We cannot create short-hands for conflicting option names that would yield ambiguities
            for depth in range(1, len(cur_path)+1):
                shorthand = '.'.join(cur_path[depth:] + [key,])
                if shorthand not in self._shorthands:
                    self._shorthands[shorthand] = '.'.join(cur_path + [key,])
            if isinstance(value, dict) and not key.endswith('_dict'):
                self.update_shorthands(
                    cur_path=cur_path + [key,], root_dict=value)

    def _update_config_chunk(self, root_path: str, config_chunk: dict[str, Any], updater: Any) -> None:
        for key, value in updater.items():
            if root_path == '':
                setting_path = key
            else:
                setting_path = '.'.join([root_path, key])
            if not isinstance(key, str):
                raise GammaLoopError(
                    f"Invalid path for setting {setting_path}")
            if key not in config_chunk:
                # Allow to create new keys for the rust run settings configuration
                if 'run_settings' in root_path.split('.'):
                    config_chunk[key] = value
                    self.update_shorthands(cur_path=root_path.split(
                        '.')+[key,], root_dict=config_chunk)
                    continue
                raise GammaLoopError(
                    f"No settings {setting_path} in gammaloop configuration.")
            if isinstance(value, dict):
                if key.endswith('_dict'):
                    config_chunk[key] = updater
                    continue
                if any(not isinstance(k, str) for k in updater):
                    if all(isinstance(k, str) for k in config_chunk[key]):
                        raise GammaLoopError(
                            f"Invalid value for setting {setting_path}:\n{pformat(updater)}")
                    else:
                        config_chunk[key] = updater
                        continue
                else:
                    self._update_config_chunk(
                        setting_path, config_chunk[key], value)
            else:
                if value is not None and config_chunk[key] is not None and type(value) is not type(config_chunk[key]):
                    raise GammaLoopError(
                        f"Invalid value for setting {setting_path}. Default value of type '{type(config_chunk[key]).__name__}' is:\n{pformat(config_chunk[key])}\nand you supplied this value of type '{type(value).__name__}':\n{pformat(value)}")
                config_chunk[key] = value
                continue

    def update(self, new_setting: Any, path: str = '', allow_shorthands: bool = True) -> None:
        context = self._config
        for key in path.split('.'):
            if key == "":
                break
            if key not in context:
                raise GammaLoopError(
                    f"No settings '{path}' in gammaloop configuration.")
            context = context[key]
        self._update_config_chunk(path, context, new_setting)

    def set_setting(self, path: str, new_setting: Any, allow_shorthands: bool = True) -> str:
        if allow_shorthands and path in self._shorthands:
            return self.set_setting(self._shorthands[path], new_setting, allow_shorthands=False)

        p = path.split('.')
        self.update({p[-1]: new_setting}, path=('.'.join(p[:-1])
                    if len(p) > 1 else ''), allow_shorthands=allow_shorthands)
        return path

    def get_setting(self, path: str, allow_shorthands: bool = True) -> Any:
        if allow_shorthands and path in self._shorthands:
            return self.get_setting(self._shorthands[path], allow_shorthands=False)
        context = self._config
        for key in path.split('.'):
            if key == "":
                break
            if key not in context:
                raise GammaLoopError(
                    f"No settings '{path}' in gammaloop configuration.")
            context = context[key]
        return context

    def __getitem__(self, key: str) -> Any:
        if key in self._config:
            return self._config[key]
        else:
            raise GammaLoopError(f"Unknown gammloop setting '{key}'")

    def __str__(self) -> str:
        return verbose_yaml_dump(self._config)


def split_str_args(str_args: str) -> list[str]:
    return str_args.split(' ') if str_args != '' else []


class CommandList(list[tuple[str, str]]):

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
        if cmd.startswith('#'):
            return
        if cmd_split[0] not in AVAILABLE_COMMANDS:
            raise GammaLoopError(f"Unknown command: {cmd_split[0]}")
        if len(cmd_split) >= 2:
            self.append((cmd_split[0], cmd_split[1]))
        else:
            self.append((cmd_split[0], ''))

    def parse(self, cmd_str: str) -> None:
        multi_line = ''
        for line in sum([t.split('\n') for t in cmd_str.strip().split(';')], []):  # type: ignore
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
        self.model: Model = Model('NotLoaded')
        self.model_directory: Path = Path('NotLoaded')

        # Initialize a gammaloop rust engine worker which will be used throughout the session
        self.rust_worker: gl_rust.Worker = gl_rust.Worker.new()
        logger.debug(
            "Successfully initialized GammaLoop rust worker from library %s.", gl_rust.__file__)  # type: ignore

        self.cross_sections: cross_section.CrossSectionList = cross_section.CrossSectionList()
        self.amplitudes: cross_section.AmplitudeList = cross_section.AmplitudeList()

        self.config: GammaLoopConfiguration = GammaLoopConfiguration()
        self.default_config: GammaLoopConfiguration = copy.deepcopy(
            self.config)
        self.launched_output: Path | None = None
        self.command_history: CommandList = CommandList()

        if gl_is_symbolica_registered is False:
            raise GammaLoopError("Symbolica is not registered and since gammaLoop uses it both within Python and Rust, multiple instances of Symbolica are necessary.\n"
                                 + "Please register Symbolica by setting the environment variable 'SYMBOLICA_LICENSE' or by adding it to the gammaloop configuration file.\n"
                                 + "Also make sure you have a working internet connection. Alternatively, revert the license setting to the default value 'GAMMALOOP_USER'.")

    def get_model_from_rust_worker(self) -> Model:
        return Model.from_yaml(self.rust_worker.get_model())

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
            logger.debug("Running command '%s%s%s' with arguments '%s%s%s'",
                         Colour.GREEN, cmd, Colour.END, Colour.BLUE, args, Colour.END)
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

    # reset command
    reset_parser = ArgumentParser(prog='set')
    reset_parser.add_argument('path', metavar='path', type=str,
                              help='Setting path to reset default value for')

    def do_reset(self, str_args: str) -> None:
        if str_args == 'help':
            self.reset_parser.print_help()
            return
        args = self.reset_parser.parse_args(split_str_args(str_args))
        default_value = self.default_config.get_setting(args.path)
        full_path = self.config.set_setting(args.path, default_value)
        logger.info("Setting '%s%s%s' successfully reset its default value: %s%s%s",
                    Colour.GREEN, full_path, Colour.END, Colour.BLUE, repr(default_value), Colour.END)

    # set command
    set_parser = ArgumentParser(prog='set')
    set_parser.add_argument('path', metavar='path', type=str,
                            help='Setting path to set value to')
    set_parser.add_argument(
        'value', metavar='value', type=str, help='value to set as valid python syntax')

    def do_set(self, str_args: str) -> None:
        if str_args == 'help':
            self.set_parser.print_help()
            return
        args = self.set_parser.parse_args(split_str_args(str_args))

        try:
            config_value = eval(args.value)  # pylint: disable=eval-used
        except Exception as exc:
            raise GammaLoopError(
                f"Invalid value '{args.value}' for setting '{args.path}'. Error:\n{exc}") from exc

        full_path = self.config.set_setting(args.path, config_value)
        str_setting = pformat(config_value)
        logger.info("Setting '%s%s%s' to:%s%s%s%s", Colour.GREEN,
                    full_path, Colour.END, Colour.BLUE, '\n' if len(str_setting) > 80 else ' ', str_setting, Colour.END)

    # show_settings command
    show_settings = ArgumentParser(prog='show_settings')
    show_settings.add_argument('path', metavar='path', type=str,
                               help='Setting path to show')

    def do_show_settings(self, str_args: str) -> None:
        if str_args == 'help':
            self.show_settings.print_help()
            return
        args = self.show_settings.parse_args(split_str_args(str_args))

        setting = self.config.get_setting(args.path)
        str_setting: str = pformat(setting)
        logger.info("Current value of setting %s%s%s:%s%s%s%s",
                    Colour.GREEN, args.path, Colour.END, '\n' if len(str_setting) > 80 else ' ', Colour.BLUE, str_setting, Colour.END)

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

        model_directory: str = os.path.dirname(args.model)
        if model_directory == '':
            model_directory = pjoin(GL_PATH, 'data', 'models')
        else:
            model_directory = os.path.abspath(model_directory)

        split_model_specification = os.path.basename(args.model).split('-')
        model_name: str = ''
        model_restriction: str | None = None
        if len(split_model_specification) > 3:
            raise GammaLoopError(
                "Model names cannot contain hyphens and must be specified with format <model_name>[-<restriction_name>].")
        elif len(split_model_specification) == 2:
            model_name, model_restriction = split_model_specification
        elif len(split_model_specification) == 1:
            model_name = split_model_specification[0]
            if args.format == 'ufo':
                if os.path.isfile(pjoin(model_directory, model_name, 'restrict_default.dat')):
                    model_restriction = 'default'
            else:
                if os.path.isfile(pjoin(model_directory, 'restrict_default.dat')):
                    model_restriction = 'default'

        match args.format:
            case 'ufo':
                self.model = Model.from_ufo_model(
                    pjoin(model_directory, model_name))
                model_restriction_dir = pjoin(model_directory, model_name)
            case 'yaml':
                with open(pjoin(model_directory, model_name), 'r', encoding='utf-8') as file:
                    self.model = Model.from_yaml(file.read())
                model_restriction_dir = model_directory
            case _:
                raise GammaLoopError(f"Invalid model format: '{args.format}'")

        # Make sure to issue warnings again if they were issued before
        if GammaLoopWarning.FloatInExpression in GL_WARNINGS_ISSUED:
            GL_WARNINGS_ISSUED.remove(GammaLoopWarning.FloatInExpression)

        self.model.restriction = model_restriction

        if model_restriction not in [None, 'full']:
            if not os.path.isfile(pjoin(model_restriction_dir, f'restrict_{model_restriction}.dat')):
                raise GammaLoopError(
                    f"Restriction file 'restrict_{model_restriction}.dat' not found for model '{model_name}' in directory '{pjoin(model_directory,model_name)}'.")
            else:
                param_card = ParamCard(
                    pjoin(model_restriction_dir, f'restrict_{model_restriction}.dat'))
        else:
            model_restriction = 'full'
            if not os.path.isfile(pjoin(model_restriction_dir, 'restrict_full.dat')):
                self.model.apply_input_param_card(
                    InputParamCard.default_from_model(self.model), simplify=False)
                ParamCardWriter.write(
                    Path(pjoin(model_restriction_dir, 'restrict_full.dat')), self.model, generic=True)
            param_card = ParamCard(
                pjoin(model_restriction_dir, 'restrict_full.dat'))

        if model_restriction != 'full':
            self.model.apply_input_param_card(InputParamCard.from_param_card(
                param_card, self.model), simplify=True)
        else:
            self.model.apply_input_param_card(InputParamCard.from_param_card(
                param_card, self.model), simplify=False)

        self.model_directory = Path(model_directory)

        # Assign the UFO model to our rust worker
        self.rust_worker.load_model_from_yaml_str(self.model.to_yaml())
        logger.info("Successfully loaded model '%s'.",
                    self.model.get_full_name())

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
            case _:
                raise GammaLoopError(
                    "Invalid model format: '%s' for exporting model.", args.format)
        logger.info("Successfully exported model '%s' to '%s'.",
                    self.model.get_full_name(), args.model_file_path)

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

        match args.format:
            case 'yaml':
                try:
                    all_raw_graphs: list[Any] = yaml.safe_load(
                        open(file_path, 'r', encoding='utf-8'))['graphs']
                except Exception as exc:
                    raise GammaLoopError(
                        f"Error while loading graphs from YAML file '{args.file_path}'. Error:\n{exc}") from exc

            case 'qgraph':
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
                all_raw_graphs: list[Any] = qgraph_loaded_module.graphs

            case _:
                raise GammaLoopError(
                    "Invalid graph format: '%s' for importing graphs.", args.format)

        graphs: list[Graph] = []
        for i_qg, qgraph_object in enumerate(all_raw_graphs):
            graphs.append(Graph.from_qgraph(
                self.model, qgraph_object, name=f"{file_path.stem}_{i_qg}"))
        logger.info("Successfully loaded %s graphs.",
                    len(all_raw_graphs))

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
                "%sProcessing of forward scattering graphs not fully implemented yet.%s", Colour.RED, Colour.END)
            self.cross_sections = cross_section.CrossSectionList([cross_section.CrossSection(
                f'{file_path.stem}',
                # Wrap the forward scattering graphs within a dummy supergraph
                [cross_section.supergraph.SuperGraph(
                    sg_id=i, graph=Graph.empty_graph('DUMMY'), multiplicity=1.0,
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
                    for i, g in enumerate(graphs)]
            )])

    # output command
    output_parser = ArgumentParser(prog='output')
    output_parser.add_argument(
        'output_path', type=str, help='Path to output the cross section to')
    output_parser.add_argument('-exp', '--expression', default=False, action='store_true',
                               help='Generate expression associated to the graph and output it to a text file')
    output_parser.add_argument('-mr', '--model_replacements', default=False, action='store_true',
                               help='Generate coupling replacements and output it to a text file')
    output_parser.add_argument('-ef', '--expression_format', type=str, default='file',
                               choices=['file', 'mathematica', 'latex'], help='Format to export symbolica objects in the numerator output.')

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

        if args.model_replacements:
            self.rust_worker.export_coupling_replacement_rules(
                args.output_path, args.expression_format)
            logger.info(
                "Model replacement rules exported to model directory.")

        if len(self.cross_sections) == 0 and len(self.amplitudes) == 0:
            raise GammaLoopError("No process generated yet.")

        if len(self.cross_sections) > 0:
            cross_section_exporter = CrossSectionsExporter(self, args)
            cross_section_exporter.export(
                args.output_path, self.cross_sections)
            logger.info("Cross-sections exported to '%s'.", args.output_path)

        if len(self.amplitudes) > 0:
            amplitude_exporter = AmplitudesExporter(self, args)
            amplitude_exporter.export(args.output_path, self.amplitudes)
            if args.expression:
                amplitude_exporter.export_expression(
                    args.output_path, self.amplitudes, args.expression_format)

            logger.info("Amplitudes exported to '%s'.", args.output_path)
    #
    # Run interface type of commands below (those bound to a particular output already generated)
    #

    # launch command
    launch_parser = ArgumentParser(prog='launch')
    launch_parser.add_argument(
        'path_to_launch', type=str, help='Path to launch a given run command to')
    launch_parser.add_argument('--no_overwrite_model', '-nom', action='store_true', default=False,
                               help='Do not overwrite model with the new parameter and coupling values derived from the param card.')
    launch_parser.add_argument('--no_overwrite_run_settings', '-nors', action='store_true', default=False,
                               help='Do not overwrite the run settings with the new settings derived from the run card.')

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
            self.model = Model.from_yaml(file.read())
        self.model.apply_input_param_card(InputParamCard.from_param_card(
            ParamCard(pjoin(args.path_to_launch, 'cards', 'param_card.dat')), self.model), simplify=False)
        processed_yaml_model = self.model.to_yaml()
        if not args.no_overwrite_model:
            with open(pjoin(args.path_to_launch, 'sources', 'model', f"{output_metadata['model_name']}.yaml"), 'w', encoding='utf-8') as file:
                file.write(processed_yaml_model)
        self.rust_worker.load_model_from_yaml_str(processed_yaml_model)

        # Sync the run_settings from the interface to the output or vice-versa
        if not args.no_overwrite_run_settings:
            run_settings = load_configuration(
                pjoin(args.path_to_launch, 'cards', 'run_card.yaml'), True)
            # Do not overwrite the external momenta setting in gammaLoop if they are set constant in the process dir
            if run_settings['Kinematics']['externals']['type'] == 'constant':
                del run_settings['Kinematics']
            self.config.update({'run_settings': run_settings})

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
            self.rust_worker.load_amplitudes_derived_data(
                pjoin(args.path_to_launch, 'sources', 'amplitudes'))

            self.rust_worker.load_amplitude_integrands(
                pjoin(args.path_to_launch, 'cards', 'run_card.yaml'))

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
    inspect_parser.add_argument(
        "integrand", type=str, help="Integrand to inspect.")
    inspect_parser.add_argument('--use-f128', '-f128', action='store_true',
                                default=False, help='Use f128 precision for the inspection.')
    inspect_parser.add_argument('--point', '-p', nargs="+", type=float,
                                default=None, help='Point to inspect.')
    inspect_parser.add_argument(
        '--term', '-t', nargs="+", type=int, default=None, help="term to inspect.")
    inspect_parser.add_argument('--is-momentum-space', '-ms', action='store_true',
                                default=False, help='Inspect in momentum space.')
    inspect_parser.add_argument(
        '--force-radius', '-fr', action='store_true', default=False, help='Force radius to be used.')

    def do_inspect(self, str_args: str) -> None:
        if str_args == 'help':
            self.inspect_parser.print_help()
            return
        args = self.inspect_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        self.rust_worker.inspect_integrand(
            args.integrand, args.point, args.term, args.force_radius, args.is_momentum_space, args.use_f128)

        # raise GammaLoopError("Command not implemented yet")

    # integrate command
    integrate_parser = ArgumentParser(prog='integrate')
    integrate_parser.add_argument(
        "integrand", type=str, help="Integrand to integrate.")
    integrate_parser.add_argument('--cores', '-c', type=int, default=1,)
    integrate_parser.add_argument(
        '--target', '-t', type=str, default=None)
    integrate_parser.add_argument(
        '--restart', '-r', action='store_true',)

    def do_integrate(self, str_args: str) -> None:
        if str_args == 'help':
            self.integrate_parser.print_help()
            return
        args = self.integrate_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        # Update the run card in the output with the current configuration
        update_run_card_in_output(self.launched_output, self.config)
        self.rust_worker.load_amplitude_integrands(
            pjoin(self.launched_output, 'cards', 'run_card.yaml'))

        target: tuple[float, float] | None = None
        if args.target is not None:
            try:
                tmp: Any = eval(args.target)
                target = (float(tmp[0]), float(tmp[1]))
            except Exception as exc:
                raise GammaLoopError(
                    f"Invalid target '{args.target}'. It should be a 2-tuple identifying a complex numbers. Error:\n{exc}") from exc

        result_output_path = self.launched_output.joinpath(
            "runs").joinpath("run.yaml")

        workspace_path = self.launched_output.joinpath("workspace")
        if args.restart and os.path.exists(workspace_path):
            shutil.rmtree(workspace_path)

        if not os.path.exists(workspace_path):
            os.mkdir(workspace_path)

        self.rust_worker.integrate_integrand(
            args.integrand, args.cores, str(result_output_path), str(workspace_path), target)

        # nuke the workspace if integration finishes
        # For now leave the possibility of restarting where integration left off.
        # Maybe in the future add an option to automatically clean the workspace after running is completed or
        # specify a "run_tag" that allows to have mutliple workspace concurrently active
        # shutil.rmtree(workspace_path)

    # test_ir_limits
    test_ir_limits_parser = ArgumentParser(prog='test_ir_limits')

    def do_test_ir_limits(self, str_args: str) -> None:
        if str_args == 'help':
            self.test_ir_limits_parser.print_help()
            return
        _args = self.test_ir_limits_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        # Update the run card in the output with the current configuration
        update_run_card_in_output(self.launched_output, self.config)
        self.rust_worker.load_amplitude_integrands(
            pjoin(self.launched_output, 'cards', 'run_card.yaml'))

        raise GammaLoopError("Command not implemented yet")

    # test_uv_limits
    test_uv_limits_parser = ArgumentParser(prog='test_ir_limits')

    def do_test_uv_limits(self, str_args: str) -> None:
        if str_args == 'help':
            self.test_uv_limits_parser.print_help()
            return
        _args = self.test_uv_limits_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        # Update the run card in the output with the current configuration
        update_run_card_in_output(self.launched_output, self.config)
        self.rust_worker.load_amplitude_integrands(
            pjoin(self.launched_output, 'cards', 'run_card.yaml'))

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

    hpc_parser = ArgumentParser(prog='hpc_run')
    hpc_parser.add_argument('integrand', type=str,
                            help="Integrand to integrate", default=None)

    def do_hpc_run(self, str_args: str) -> None:
        args = self.hpc_parser.parse_args(split_str_args(str_args))
        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        # Update the run card in the output with the current configuration
        update_run_card_in_output(self.launched_output, self.config)
        self.rust_worker.load_amplitude_integrands(
            pjoin(self.launched_output, 'cards', 'run_card.yaml'))

        # create a workspace for batch input/output files
        workspace_path = self.launched_output.joinpath("workspace")

        if not workspace_path.exists():
            workspace_path.mkdir()

        from hyperqueue import LocalCluster, Job  # type: ignore
        from hyperqueue.cluster import WorkerConfig  # type: ignore

        print("starting hpc test run")
        self.rust_worker.load_master_node(args.integrand)

        # this will be loaded from settings, eventually with dynamic points that can increase with iterations
        n_iterations = 10
        n_tasks = 4
        n_points_per_task = 1000000
        n_cores = 1

        export_grid = False
        output_accumulator = False

        # this local cluster is just for testing
        with LocalCluster() as cluster:
            cluster.start_worker()

            client = cluster.client()  # type: ignore

            for _ in range(n_iterations):
                task_ids = [i*(n_iterations + 1) for i in range(n_tasks)]
                job = Job()

                for id in task_ids:
                    input_file = workspace_path.joinpath(
                        "job_{}".format(str(id)))

                    output_file = workspace_path.joinpath(
                        "job_{}_out".format(str(id)))

                    self.rust_worker.write_batch_input(
                        n_cores, n_points_per_task, export_grid, output_accumulator, str(workspace_path), id)

                    job.program(["bin/gammaloop_rust_cli", "batch", "--batch_input_file={}".format(str(input_file)),
                                "--name=massless_triangle", "--process_file=triangle", "--output_name={}".format(str(output_file))])

                submitted = client.submit(job)
                client.wait_for_jobs([submitted])

                for id in task_ids:
                    self.rust_worker.process_batch_output(
                        str(workspace_path), id)

                self.rust_worker.update_iter()
                self.rust_worker.display_master_node_status()

        shutil.rmtree(workspace_path)
