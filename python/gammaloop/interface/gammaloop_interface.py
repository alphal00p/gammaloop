#!/usr/bin/env python3
import platform
import sys
from pathlib import Path
import importlib
from argparse import ArgumentParser, BooleanOptionalAction, Namespace
import subprocess
import os
import time
from typing import Any, Dict
from pprint import pformat
import yaml  # type: ignore
import shutil
import copy
import logging
import pydot
from gammaloop import __version__
from gammaloop.misc.common import GammaLoopError, logger, Side, pjoin, load_configuration, GL_CONSOLE_HANDLER, GAMMALOOP_CONFIG_PATHS, gl_is_symbolica_registered, GL_PATH, GL_WARNINGS_ISSUED, GIT_REVISION, GammaLoopWarning
from gammaloop.misc.utils import Colour, verbose_yaml_dump, format_elapsed
from gammaloop.base_objects.model import Model, InputParamCard
from gammaloop.base_objects.param_card import ParamCard, ParamCardWriter
from gammaloop.base_objects.graph import Graph
from gammaloop.base_objects.process import Process
import gammaloop.cross_section.cross_section as cross_section
import gammaloop.cross_section.supergraph as supergraph
from gammaloop.exporters.exporters import AmplitudesExporter, CrossSectionsExporter, OutputMetaData, update_run_card_in_output
# This is the pyo3 binding of the gammaloop rust engine
import gammaloop._gammaloop as gl_rust
import gammaloop.interface.debug_display as debug_display
# pylint: disable=unused-variable


class TaggedScalar(yaml.ScalarNode):
    def __init__(self, tag, value):
        super().__init__(tag, value)


def tagged_scalar_representer(dumper, data):
    return data


yaml.add_representer(TaggedScalar, tagged_scalar_representer)

AVAILABLE_COMMANDS = [
    'import_model',
    'export_model',
    'import_graphs',
    'export_graphs',
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
    'set_model_param',
    'reset',
    'define_graph',
    'generate',
    'generate_amplitude',
    'display_debug_log'
]


class GammaLoopConfiguration(object):

    def __init__(self, path: str | None = None, quiet=False):
        self._shorthands: dict[str, str] = {}
        self._config: dict[str, Any] = {
            'symbolica': {
                'license': "GAMMALOOP_USER"
            },
            'drawing': {
                'modes': ['feynmp', "dot"],
                'combined_graphs_pdf_grid_shape': [3, 2],
                'dot': {
                    'layout': 'neato',
                    'show_overall_factor': True,
                    'graph_options': {
                        'fontsize': 10,
                        'ratio': 1.5
                    },
                    'node_options': {
                        'fontsize': 7,
                        'shape': 'circle',
                        'margin': 0,
                        'height': 0.01,
                        'penwidth': 0.6
                    },
                    'edge_options': {
                        'fontsize': 7,
                        'arrowsize': 0.3,
                        'penwidth': 0.6,
                    },
                },
                'feynmp': {
                    'reverse_outgoing_edges_order': True,
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
                'compile_cff': True,
                'numerator_settings': {
                    'dump_expression': 'Mathematica',
                    'eval_settings': {
                        'type': 'Joint',
                        'cpe_rounds': 1,
                        'compile_options': {
                            'subtype': 'Compiled',
                        }
                    },
                    'global_numerator': None,
                    'global_prefactor': {'color': "1", 'colorless': "1"},
                    'dump_expression': 'Mathematica',
                    'gamma_algebra': 'Concrete',
                    'parse_mode': 'Polynomial'
                },
                'cpe_rounds_cff': 1,
                'compile_separate_orientations': False,
                'gammaloop_compile_options': {
                    'inline_asm': True,
                    'optimization_level': 3,
                    'fast_math': True,
                    'unsafe_math': True,
                    'compiler': 'g++',
                    'custom': []
                },
                'tropical_subgraph_table_settings': {
                    'panic_on_fail': False,
                    'target_omega': 1.0
                }
            },
            'run_settings': {
                'General': {
                    'debug': 0,
                    'use_ltd': False,
                    'force_orientations': None,
                    'load_compiled_cff': True,
                    'load_compiled_numerator': True,
                    'joint_numerator_eval': True,
                    'amplitude_prefactor': {
                        're': 0.0,
                        'im': 1.0
                    },
                    'load_compiled_separate_orientations': False
                },
                'Integrand': {
                    'type': 'gamma_loop'
                },
                'Kinematics': {
                    'e_cm': 3.0,
                    'externals': {
                        'type': 'constant',
                        'data': {
                            'momenta': None,  # Will be set automatically
                            'helicities': None,
                        }
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
                    'discrete_dim_learning_rate': 1.5,
                    'continuous_dim_learning_rate': 0.0,
                    'train_on_avg': False,
                    'show_max_wgt_info': False,
                    'max_prob_ratio': 0.01,
                    'seed': 0
                },
                'Observables': [],
                'Selectors': [],
                'Stability': {
                    'rotation_axis': [{'type': 'x'}],
                    'levels': [
                        {
                            'precision': 'Double',
                            'required_precision_for_re': 1.e-7,
                            'required_precision_for_im': 1.e-7,
                            'escalate_for_large_weight_threshold': 0.9
                        },
                        {
                            'precision': 'Quad',
                            'required_precision_for_re': 1.e-10,
                            'required_precision_for_im': 1.e-10,
                            'escalate_for_large_weight_threshold': -1.0
                        }
                    ],
                    'rotate_numerator': False,
                },
                'sampling': {
                    'type': 'discrete_graph_sampling',
                    'subtype': 'tropical',
                    'upcast_on_failure': False
                },
                'subtraction': {
                    'local_ct_settings': {
                        'uv_localisation': {
                            'sliver_width': 10.0,
                            'dynamic_width': False,
                            'gaussian_width': 1.0
                        },
                        'dampen_integrable_singularity': {
                            'type': 'exponential'
                        }
                    },
                    'integrated_ct_settings': {
                        'range': {
                            'type': 'infinite',
                            'h_function_settings': {
                                'function': 'poly_exponential',
                                'sigma': 1.0,
                                'enabled_dampening': True,
                                'power': None,
                            }

                        }
                    },
                    'overlap_settings': {
                        'force_global_center': None,
                        'check_global_center': True,
                        'try_origin': False,
                        'try_origin_all_lmbs': False,
                    }
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

    def _update_config_chunk(self, root_path: str, config_chunk: dict[str, Any] | None, updater: Any) -> None:

        for key, value in updater.items():
            if root_path == '':
                setting_path = key
            else:
                setting_path = '.'.join([root_path, key])
            if not isinstance(key, str):
                raise GammaLoopError(
                    f"Invalid path for setting {setting_path}")
            if config_chunk is None:
                raise GammaLoopError(
                    f"Settting key '{key}' not found in parameters. You can force the setting of that key by wrapping the value to be set within quotes.")
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

                    if config_chunk[key] is None:
                        config_chunk[key] = value
                    else:
                        self._update_config_chunk(
                            setting_path, config_chunk[key], value)

            else:
                if value is not None and config_chunk[key] is not None and type(value) is not type(config_chunk[key]):
                    if isinstance(value, str) and isinstance(config_chunk[key], dict):
                        try:
                            value = eval(value)
                        except:
                            raise GammaLoopError(f"Invalid value for setting {setting_path}. It is a string that needs to evaluate to a python dictionary:\n{pformat(updater)}")  # nopep8
                        if not isinstance(value, dict):
                            raise GammaLoopError(f"Invalid value for setting {setting_path}. It is a string that needs to evaluate to a python dictionary:\n{pformat(updater)}")  # nopep8
                    else:
                        raise GammaLoopError(
                            f"Invalid value for setting {setting_path}. Default value of type '{type(config_chunk[key]).__name__}' is:\n{pformat(config_chunk[key])}\nand you supplied this value of type '{type(value).__name__}':\n{pformat(value)}")
                config_chunk[key] = value
                continue

    def update(self, new_setting: Any, path: str = '', allow_shorthands: bool = True) -> None:
        context = self._config
        p_split = path.split('.')
        for key in p_split[:-1]:
            if key == "":
                break
            if key not in context:
                raise GammaLoopError(
                    f"No settings '{path}' in gammaloop configuration.")
            context = context[key]
        if path != '' and isinstance(new_setting, dict):
            context[p_split[-1]].update(new_setting)
        else:
            if path != '':
                context = context[p_split[-1]]
            self._update_config_chunk(path, context, new_setting)

    def set_setting(self, path: str, new_setting: Any, allow_shorthands: bool = True) -> str:
        p = path.split('.')
        if allow_shorthands and p[0] in self._shorthands:
            return self.set_setting('.'.join([self._shorthands[p[0]],]+p[1:]), new_setting, allow_shorthands=False)

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

    @ staticmethod
    def from_file(filename: str):
        command_file = CommandList()
        with open(filename, 'r', encoding='utf-8') as file:
            command_file.parse(file.read())
        return command_file

    @ staticmethod
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
        if cmd.startswith('#'):
            return
        # Allow inline comments
        cmd_split: list[str] = cmd.split('#', 1)[0].strip().split(' ', 1)

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
        self.rust_worker: gl_rust.Worker = gl_rust.Worker()

        logger.debug('Starting interface of GammaLoop v%s%s',
                     __version__,
                     (f' (git rev.: {GIT_REVISION})' if GIT_REVISION is not None else ''))
        logger.debug(
            "Successfully initialized GammaLoop rust worker from library %s.", gl_rust.__file__)  # type: ignore

        self.cross_sections: cross_section.CrossSectionList = cross_section.CrossSectionList()
        self.amplitudes: cross_section.AmplitudeList = cross_section.AmplitudeList()

        self.config: GammaLoopConfiguration = GammaLoopConfiguration()
        self.default_config: GammaLoopConfiguration = copy.deepcopy(
            self.config)
        self.launched_output: Path | None = None
        self.process: Process | None = None
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
        'value', metavar='value', type=str, nargs="*", help='value to set as valid python syntax')
    set_parser.add_argument(
        '--no_update_run_card', '-n', default=False, action='store_true',
        help='Do not update the run card of the current launched output. This is useful for batching changes.')

    def do_set(self, str_args: str) -> None:
        if str_args == 'help':
            self.set_parser.print_help()
            return
        args = self.set_parser.parse_args(split_str_args(str_args))

        try:
            config_value = eval(''.join(args.value))  # nopep8 # pylint: disable=eval-used
        except Exception as exc:
            raise GammaLoopError(
                f"Invalid value '{args.value}' for setting '{args.path}'. Error:\n{exc}") from exc

        full_path = self.config.set_setting(args.path, config_value)
        str_setting = pformat(config_value)
        logger.info("Setting '%s%s%s' to:%s%s%s%s", Colour.GREEN,
                    full_path, Colour.END, Colour.BLUE, '\n' if len(str_setting) > 80 else ' ', str_setting, Colour.END)
        if not args.no_update_run_card and self.launched_output is not None:
            # Update the run card in the output with the current configuration
            update_run_card_in_output(self.launched_output, self.config)

    # show_settings command
    set_model_param_settings = ArgumentParser(prog='set_model_param')
    set_model_param_settings.add_argument('param', metavar='param', type=str,
                                          help='Model external parameter name to modify')
    set_model_param_settings.add_argument('value', metavar='value', type=float,
                                          help='Value to assign to that model parameter')
    set_model_param_settings.add_argument('--no_update', '-nu',
                                          default=False, action='store_true', help='Do not update dependent model parameters yet')
    set_model_param_settings.add_argument('--no_overwrite', '-no',
                                          default=False, action='store_true', help='Do not overwrite the param card and YAML model on disk with the new value')

    def do_set_model_param(self, str_args: str) -> None:
        if str_args == 'help':
            self.set_model_param_settings.print_help()
            return
        args = self.set_model_param_settings.parse_args(
            split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        if self.model.is_empty():
            raise GammaLoopError(
                "Command set_model_param can only be called once a model has been loaded.")

        input_card = InputParamCard.from_model(self.model)

        if args.param != 'update_only':
            if args.param not in input_card:
                raise GammaLoopError(
                    f"Model parameter '{args.param}' not found in model '{self.model.get_full_name()}'. Available external parameters are:\n\
                        {', '.join(p for p in input_card.keys())}")

            input_card[args.param] = complex(args.value)

        self.model.apply_input_param_card(
            input_card, simplify=False, update=not args.no_update)

        if args.param != 'update_only':
            logger.info("Setting model parameter '%s%s%s' to %s%s%s", Colour.GREEN,
                        args.param, Colour.END, Colour.BLUE, args.value, Colour.END)
        else:
            logger.info("Updating all dependent model parameters")

        if not args.no_update:
            processed_yaml_model = self.model.to_yaml()
            self.rust_worker.load_model_from_yaml_str(processed_yaml_model)

            if not args.no_overwrite:
                ParamCardWriter.write(
                    self.launched_output.joinpath('cards', 'param_card.dat'), self.model, generic=True)
                logger.debug("Successfully updated param card '%s' with new value of parameter '%s%s%s' to %s%f%s",
                             self.launched_output.parent.joinpath('cards', 'param_card.dat'), Colour.GREEN, args.param, Colour.END, Colour.BLUE, args.value, Colour.END)
                with open(self.launched_output.joinpath('output_metadata.yaml'), 'r', encoding='utf-8') as file:
                    output_metadata = OutputMetaData.from_yaml_str(file.read())
                self.launched_output.joinpath('cards', 'param_card.dat')
                with open(self.launched_output.joinpath('sources', 'model', f"{output_metadata['model_name']}.yaml"), 'w', encoding='utf-8') as file:
                    file.write(processed_yaml_model)
                logger.debug("Successfully updated YAML model sources '%s'.", self.launched_output.joinpath(
                    'sources', 'model', f"{output_metadata['model_name']}.yaml"))

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
                    f"Restriction file 'restrict_{model_restriction}.dat' not found for model '{model_name}' in directory '{pjoin(model_directory, model_name)}'.")
            else:
                param_card = ParamCard(
                    pjoin(model_restriction_dir, f'restrict_{model_restriction}.dat'))
        else:
            model_restriction = 'full'
            if not os.path.isfile(pjoin(model_restriction_dir, 'restrict_full.dat')):
                self.model.apply_input_param_card(
                    InputParamCard.from_model(self.model), simplify=False)
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

    # generate command
    generate_parser = ArgumentParser(prog='generate')
    generate_parser.add_argument(
        'process', metavar='process', type=str, nargs="+", help='Process to generate.')
    generate_parser.add_argument('--clear_existing_processes', '-clear', default=False, action='store_true',
                                 help='Clear existing processes stored before adding this one..')
    generate_parser.add_argument('--num_threads', '-nt', default=None, type=int,
                                 help='Number of threads to parallelize the generation on. (default: all available cores)')
    generate_parser.add_argument('--max_multiplicity_for_fast_cut_filter', '-mmfcf', default=6, type=int,
                                 help='Specify the maximum cut multiplicity before switching to an alternative Cutkosky cut filter that is faster for large multiplicities. (default: %(default)s)')
    generate_parser.add_argument('--graph_prefix', '-gp', type=str, default="GL",
                                 help='Graph name prefix. default: "GL"')
    generate_parser.add_argument('--max_n_bridges', '-mnb', type=int, default=None,
                                 help='Specify the maximum number of bridges for the graphs to generate. Set negative to disable. (default: 0)')
    generate_parser.add_argument('--filter_self_loop', default=False, action=BooleanOptionalAction,
                                 help='Filter all self-loops directly during generation.')
    generate_parser.add_argument('--numerator_aware_isomorphism_grouping', '-num_grouping', default="group_identical_graphs_up_to_scalar_rescaling", type=str,
                                 choices=["no_grouping", "only_detect_zeroes", "group_identical_graphs_up_to_sign",
                                          "group_identical_graphs_up_to_scalar_rescaling"],
                                 help='Group identical diagrams after generation and including numerator (default: group_identical_graphs_up_to_scalar_rescaling)')
    generate_parser.add_argument('--number_of_fermion_loops', '-nfl', default=None, type=int,
                                 help='Number of fermion loops to consider in the amplitude generation. (default: any)')
    generate_parser.add_argument('--number_of_factorized_loop_subtopologies', '-nfactl', default=None, type=int,
                                 help='Number of factorizable loops (seoarated graph ears) to consider in the amplitude generation. (default: any)')
    # Tadpole filter

    generate_parser.add_argument('--filter_cross_section_tadpoles', default=None, action=BooleanOptionalAction,
                                 help='Filter tadpoles when sewing cross-section.')
    generate_parser.add_argument('--filter_tadpoles', default=None, action=BooleanOptionalAction,
                                 help='Filter tadpole diagrams.')
    generate_parser.add_argument('--veto_tadpole_attached_to_massive', dest='veto_tadpoles_attached_to_massive_lines', default=None, action=BooleanOptionalAction,
                                 help='Filter tadpole diagrams attached to massive lines.')
    generate_parser.add_argument('--veto_tadpole_attached_to_massless', dest='veto_tadpoles_attached_to_massless_lines', default=None, action=BooleanOptionalAction,
                                 help='Filter tadpole diagrams attached to massive lines.')
    generate_parser.add_argument('--veto_scaleless_tadpole', dest='veto_only_scaleless_tadpoles', default=None, action=BooleanOptionalAction,
                                 help='Filter scaleless! tadpole diagrams.')
    # Snail filter
    generate_parser.add_argument('--filter_snails', default=None, action=BooleanOptionalAction,
                                 help='Filter snail diagrams.')
    generate_parser.add_argument('--veto_snail_attached_to_massive', dest='veto_snails_attached_to_massive_lines', default=None, action=BooleanOptionalAction,
                                 help='Filter snail diagrams attached to massive lines.')
    generate_parser.add_argument('--veto_snail_attached_to_massless', dest='veto_snails_attached_to_massless_lines', default=None, action=BooleanOptionalAction,
                                 help='Filter snail diagrams attached to massive lines.')
    generate_parser.add_argument('--veto_scaleless_snail', dest='veto_only_scaleless_snails', default=None, action=BooleanOptionalAction,
                                 help='Filter scaleless! snail diagrams.')
    # Selfenergy filter
    generate_parser.add_argument('--filter_selfenergies', default=None, action=BooleanOptionalAction,
                                 help='Filter out (external) self energy contributions.')
    generate_parser.add_argument('--veto_selfenergy_of_massive_lines', dest='veto_self_energy_of_massive_lines', default=None, action=BooleanOptionalAction,
                                 help='Filter snail diagrams attached to massive lines.')
    generate_parser.add_argument('--veto_selfenergy_of_massless_lines', dest='veto_self_energy_of_massless_lines', default=None, action=BooleanOptionalAction,
                                 help='Filter snail diagrams attached to massive lines.')
    generate_parser.add_argument('--veto_scaleless_selfenergy', dest='veto_only_scaleless_self_energy', default=None, action=BooleanOptionalAction,
                                 help='Filter scaleless! tadpole diagrams.')
    # Symmetrization options
    generate_parser.add_argument('--symmetrize_initial_states', default=None, action=BooleanOptionalAction,
                                 help='Symmetrize initial states in diagram generation. (default: Automatic)')
    generate_parser.add_argument('--symmetrize_final_states', default=None, action=BooleanOptionalAction,
                                 help='Symmetrize final states in diagram generation. (default: Automatic)')
    generate_parser.add_argument('--symmetrize_left_right_states', '-slrs', default=None, action=BooleanOptionalAction,
                                 help='Symmetrize left and right forward scattering states in cross-section diagram generation. (default: Automatic)')
    generate_parser.add_argument('--allow_symmetrization_of_external_fermions_in_amplitudes', '-symferm', default=False, action=BooleanOptionalAction,
                                 help='Allow symmetrization of external fermions for amplitude generation. Disabled by default as it requires care from the user because of the fermion swap negative signs. (default: False)')

    # Cross-section cut options
    generate_parser.add_argument('--n_cut_blobs', '-ncb', type=int, nargs=2, default=[1, 1],
                                 help='Set the range of cut blobs on either side of the cut for cross-section generation')
    generate_parser.add_argument('--n_cut_spectators', '-ncs', type=int, nargs=2, default=[0, 0],
                                 help='Set the range of cut spectators on either side of the cut for cross-section generation')

    # Selection options
    generate_parser.add_argument('--loop_momentum_bases', '-lmbs', default=None, type=str,
                                 help='Specify LMB to use with a string corresponding to a python dictionary with format "{graph_name: list[edge_names]}" (default: Automatic)')
    generate_parser.add_argument('--select_graphs', '-selected_graphs', default=None, type=str, nargs="+",
                                 help='Select only the graphs with the specified names (default: all)')
    generate_parser.add_argument('--veto_graphs', '-veto_graphs', default=None, type=str, nargs="+",
                                 help='Veto graphs with the specified names (default: no veto)')

    # grouping options
    generate_parser.add_argument('--compare_canonized_numerator', '-can_num', default=None, action=BooleanOptionalAction,
                                 help='Enable comparison of canonized numerators when grouping diagrams. (default: Automatic)')
    generate_parser.add_argument('--number_of_samples_for_numerator_comparisons', '-n_num_samples', default=None, type=int,
                                 help='Number of numerical samples to consider when performing comparison of numerators when grouping graphs (default: Automatic)')
    generate_parser.add_argument('--consider_internal_masses_only_in_numerator_isomorphisms', '-masses_only', default=None, action=BooleanOptionalAction,
                                 help='Consider internal masses only (and not PDGS) when grouping isomorphic diagrams. (default: Automatic)')
    generate_parser.add_argument('--numerical_samples_seed', '-sseed', type=int, default=None,
                                 help='Seed for generation of numerical prime samples when grouping diagrams. (default: Automatic)')
    generate_parser.add_argument('--fully_numerical_substitution_when_comparing_numerators', '-full_num', default=None, action=BooleanOptionalAction,
                                 help='Substitute all parameters and square roots to num numerical values when sampling numerators for grouping diagrams. (default: Automatic)')

    # generate_amplitude command, without the --amplitude option
    generate_amplitude_parser = copy.deepcopy(generate_parser)
    generate_amplitude_parser.prog = 'generate_amplitude'

    generate_parser.add_argument('--amplitude', '-a', default=False,
                                 action='store_true', help='Generate an amplitude to this contribution')

    def do_generate(self, input_args: str | Namespace) -> None:

        if isinstance(input_args, str):
            if input_args == 'help':
                self.generate_parser.print_help()
                return
            split_args = split_str_args(input_args)
            args = self.generate_parser.parse_args(split_args)
        else:
            split_args = None
            args = input_args

        if args.clear_existing_processes:
            if args.amplitude:
                self.amplitudes.clear()
            else:
                self.cross_sections.clear()

        parsed_process = Process.from_input_args(self.model, args)

        if args.amplitude and parsed_process.cross_section_orders is not None:
            logger.warning(
                "Amplitude generation requested but cross section orders are also specified. Ignoring cross section orders.")
            parsed_process.cross_section_orders = None

        if args.amplitude and parsed_process.cross_section_loop_count is not None:
            logger.warning(
                "Amplitude generation requested but cross section loop count is also specified. Ignoring cross section loop count.")
            parsed_process.cross_section_loop_count = None

        # if parsed_process.perturbative_orders is not None:
        #     logger.warning(
        #         f"The nature of the specific perturbative orders specified ({' '.join('%s=%s' % (k, v) for k, v in parsed_process.perturbative_orders.items())}) is not yet supported. Only the loop count will be derived from them.")

        self.process = parsed_process

        generation_args = []
        if split_args is not None:
            aggregate = False
            for arg in split_args:
                if aggregate:
                    generation_args.append(arg)
                else:
                    if arg.startswith('-'):
                        generation_args.append(arg)
                        aggregate = True
        generation_args_str = ' '.join(generation_args)
        logger.info("Generating diagrams for process: %s%s%s%s",
                    Colour.GREEN,
                    self.process,
                    Colour.END,
                    "" if generation_args_str == '' else f" {Colour.BLUE}{generation_args_str}{Colour.END}"  # nopep8
                    )  # nopep8
        t_start = time.time()
        all_graphs: list[Graph] = self.process.generate_diagrams(
            self.rust_worker, self.model, args,
            self.config.get_setting("global_prefactor.color"),
            self.config.get_setting("global_prefactor.colorless"),
        )
        if len(all_graphs) > 0:
            logger.info("A total of %s%s%s graphs%s have been generated in %s%s%s.",
                        Colour.GREEN, Colour.BOLD, len(all_graphs), Colour.END,
                        Colour.GREEN, format_elapsed(time.time()-t_start), Colour.END)
        else:
            raise GammaLoopError(
                f"No graphs were generated for process:\n{(self.process)}.")

        if args.amplitude:
            self.amplitudes.add_amplitude(cross_section.Amplitude(
                f"{args.graph_prefix}_{self.process.process_shell_name()}",
                [
                    supergraph.AmplitudeGraph(
                        sg_id=0, sg_cut_id=0, fs_cut_id=i, amplitude_side=Side.LEFT,
                        # This is *not* the symmetry factor, but instead what would come out of grouping similar diagrams (like flavour multiplicity)
                        multiplicity="1",
                        graph=g
                    )
                    for i, g in enumerate(all_graphs)
                ]
            ))
        else:
            logger.warning(
                "%sProcessing of forward scattering graphs not fully implemented yet.%s", Colour.RED, Colour.END)
            self.cross_sections.add_cross_section(cross_section.CrossSection(
                f"{args.graph_prefix}_{self.process.process_shell_name()}",
                # Wrap the forward scattering graphs within a dummy supergraph
                [cross_section.supergraph.SuperGraph(
                    sg_id=i, graph=Graph.empty_graph('DUMMY'), multiplicity="1",
                    topology_class=[],
                    cuts=[
                        cross_section.supergraph.SuperGraphCut(
                            cut_edges=[],
                            forward_scattering_graph=cross_section.supergraph.ForwardScatteringGraph(
                                sg_id=i,
                                sg_cut_id=0,
                                multiplicity="1",
                                graph=g,
                                cuts=[]  # Will be filled in later
                            )
                        )]
                ) for i, g in enumerate(all_graphs)]
            ))

    def do_generate_amplitude(self, str_args: str) -> None:
        if str_args == 'help':
            self.generate_amplitude_parser.print_help()
            return
        args = self.generate_parser.parse_args(
            split_str_args(str_args))
        args.amplitude = True
        self.do_generate(args)

    # export_graph command
    export_graphs_parser = ArgumentParser(prog='export_graph')
    export_graphs_parser.add_argument(
        'output_path', metavar='output_path', type=str, help='Filename to write exported graphs to.')
    export_graphs_parser.add_argument('--graph_names', '-gn', type=str, nargs='*',
                                      help='Graph names to export [default: all].')
    export_graphs_parser.add_argument('--graph_container_name', '-gct', type=str, default=None,
                                      help='Graph container name to export [default: automatic].')

    def do_export_graphs(self, str_args: str) -> None:
        if str_args == 'help':
            self.export_graphs_parser.print_help()
            return
        args = self.export_graphs_parser.parse_args(split_str_args(str_args))

        graphs: list[tuple[Graph, dict[str, Any]]] = []
        if len(self.amplitudes) == 0 and len(self.cross_sections) == 0:
            raise GammaLoopError(
                "No graphs loaded. Please load graphs first with 'import_graphs' command.")

        for a in self.amplitudes:
            if args.graph_container_name is None or a.name == args.graph_container_name:
                graphs.extend((g.graph, {"multiplicity_factor": g.multiplicity})
                              for g in a.amplitude_graphs)
                break

        if len(graphs) == 0:
            for xs in self.cross_sections:
                if args.graph_container_name is None or xs.name == args.graph_container_name:
                    for xs_g in xs.supergraphs:
                        if xs_g.graph.is_empty():
                            for fwd_scattering_cut in xs_g.cuts:
                                fwd_scattering_cut.forward_scattering_graph
                                graphs.extend([(isr_cut.forward_scattering_graph.graph, {
                                              "multiplicity_factor": isr_cut.forward_scattering_graph.multiplicity}) for isr_cut in xs_g.cuts])
                        else:
                            graphs.append(
                                (xs_g.graph, {"multiplicity_factor": xs_g.multiplicity}))
                    break

        with open(args.output_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(g.to_pydot(attr).to_string()
                    for g, attr in graphs))

        if not args.graph_names is None:
            graphs = [(g, attr)
                      for g, attr in graphs if g.name in args.graph_names]

        if len(graphs) == 0:
            raise GammaLoopError(
                "No graphs found with the specified container name or graph name(s).")

        logger.info("A total of %s Graphs successfully exported to file '%s'.",
                    len(graphs), args.output_path)

    # import_graphs command
    import_graphs_parser = ArgumentParser(prog='import_graphs')
    import_graphs_parser.add_argument('file_path', metavar='file_path', type=str,
                                      help='Path to the qgraph python output to load')
    import_graphs_parser.add_argument('--no_compile', '-nc', action='store_true',
                                      default=False, help='Prevent compilation of qgraph python output.')
    import_graphs_parser.add_argument('--format', '-f', type=str, default='dot',
                                      choices=['dot', 'yaml', 'qgraph'], help='Format to import the graphs in.')
    import_graphs_parser.add_argument('--ids_to_load', '-ids', type=int, default=None, nargs='+',
                                      help='Graph indices to import (default: all).')

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

        graphs: list[tuple[Graph, dict[str, Any]]] = []
        match args.format:

            case 'dot':
                pydot_graphs = pydot.graph_from_dot_file(
                    file_path, encoding='utf-8')
                if pydot_graphs is None:
                    raise GammaLoopError(
                        "Failed to load graphs from dot file '%s'.", args.file_path)
                if args.ids_to_load is not None:
                    pydot_graphs = [
                        qg for i_g, qg in enumerate(pydot_graphs) if i_g in args.ids_to_load]
                graphs = [
                    (
                        Graph.from_pydot(self.model, pydot_g),
                        pydot_g.get_attributes()
                    ) for pydot_g in pydot_graphs
                ]

            case 'yaml':
                try:
                    all_raw_graphs: list[Any] = yaml.safe_load(
                        open(file_path, 'r', encoding='utf-8'))['graphs']
                except Exception as exc:
                    raise GammaLoopError(
                        f"Error while loading graphs from YAML file '{args.file_path}'. Error:\n{exc}") from exc
                if args.ids_to_load is not None:
                    all_raw_graphs = [
                        qg for i_g, qg in enumerate(all_raw_graphs) if i_g in args.ids_to_load]
                for i_qg, qgraph_object in enumerate(all_raw_graphs):
                    new_graph = Graph.from_qgraph(
                        self.model, qgraph_object, name=f"{file_path.stem}_{i_qg}")
                    attributes = qgraph_object
                    graphs.append((new_graph, attributes))

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
                if args.ids_to_load is not None:
                    all_raw_graphs = [
                        qg for i_g, qg in enumerate(all_raw_graphs) if i_g in args.ids_to_load]
                for i_qg, qgraph_object in enumerate(all_raw_graphs):
                    new_graph = Graph.from_qgraph(
                        self.model, qgraph_object, name=f"{file_path.stem}_{i_qg}")
                    attributes = qgraph_object
                    graphs.append((new_graph, attributes))

            case _:
                raise GammaLoopError(
                    "Invalid graph format: '%s' for importing graphs.", args.format)

        logger.info("Successfully loaded %s graphs.", len(graphs))

        # Now determine if it is a supergraph or an amplitude graph
        graph_type = None
        for a_graph, _attributes in graphs:
            if len(a_graph.get_incoming_edges()) == 0 and len(a_graph.get_outgoing_edges()) == 0:
                if graph_type is None:
                    graph_type = 'amplitude'
                elif graph_type != 'amplitude':
                    raise GammaLoopError(
                        "Mixed type of graphs (amplitude, forward scattering graph and/or amplitude graph) found in {}.".format(args.file_path))
            elif len(a_graph.get_incoming_edges()) == len(a_graph.get_outgoing_edges()) == len(a_graph.external_connections) and not any(None in c for c in a_graph.external_connections):
                if graph_type is None:
                    graph_type = 'forward_scattering_graph'
                elif graph_type != 'forward_scattering_graph':
                    raise GammaLoopError(
                        "Mixed type of graphs (supergraph, forward scattering graph and/or amplitude graph) found in {}.".format(args.file_path))
            else:
                if graph_type is None:
                    graph_type = 'amplitude'
                elif graph_type != 'amplitude':
                    raise GammaLoopError(
                        "Mixed type of graphs (supergraph, forward scattering graph and/or amplitude graph) found in {}.".format(args.file_path))

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
                    sg_id=i, graph=Graph.empty_graph('DUMMY'), multiplicity="1",
                    topology_class=[],
                    cuts=[
                        cross_section.supergraph.SuperGraphCut(
                            cut_edges=[],
                            forward_scattering_graph=cross_section.supergraph.ForwardScatteringGraph(
                                sg_id=i,
                                sg_cut_id=0,
                                multiplicity=attributes.get(
                                    "multiplicity_factor", "1"),
                                graph=g,
                                cuts=[]  # Will be filled in later
                            )
                        )]
                ) for i, (g, attributes) in enumerate(graphs)]
            )])

        if graph_type == 'amplitude':
            self.amplitudes = cross_section.AmplitudeList([cross_section.Amplitude(
                file_path.stem,
                [
                    supergraph.AmplitudeGraph(
                        sg_id=0, sg_cut_id=0, fs_cut_id=i, amplitude_side=Side.LEFT,
                        multiplicity=attributes.get(
                            "multiplicity_factor", "1"),
                        graph=g
                    )
                    for i, (g, attributes) in enumerate(graphs)]
            )])

    define_graph_parser = ArgumentParser(prog='define_graph')
    define_graph_parser.add_argument(
        '--name', '-n', type=str, default='graph', help='Name of the graph')
    define_graph_parser.add_argument(
        '--virtual-edges', '-ve', type=str, help='List of virtual edges to generate the graph from')
    define_graph_parser.add_argument(
        '--external-edges', '-ee', type=str, help='List of external edges to generate the graph from')
    define_graph_parser.add_argument(
        '--multiplicity_factor', '-mf', type=str, default="1", help="Multiplicity factor of the graph (default: '%(default)s')")
    define_graph_parser.add_argument(
        '--overall_factor', '-of', type=str, default="1", help="Overall factor of the graph (default: '%(default)s')")

    def do_define_graph(self, str_args: str) -> None:
        if str_args == 'help':
            self.define_graph_parser.print_help()
            return

        args = self.define_graph_parser.parse_args(split_str_args(str_args))

        try:
            virtual_edges = eval(args.virtual_edges)
        except Exception as exc:
            raise GammaLoopError(
                f"Invalid value '{args.virtual_edges}' for virtual edges. Error:\n{exc}") from exc

        try:
            external_edges = eval(args.external_edges)
        except Exception as exc:
            raise GammaLoopError(
                f"Invalid value '{args.external_edges}' for external edges. Error:\n{exc}") from exc

        graph: dict[str, Any] = {
        }

        graph["edges"] = {}
        graph["nodes"] = {}
        graph["overall_factor"] = args.overall_factor
        graph["multiplicity_factor"] = args.multiplicity_factor

        for external_edge in external_edges:
            type = external_edge[0]
            internal_node = external_edge[1]

            external_node = int("10" + str(internal_node))
            momentum_name = "p" + str(internal_node)

            graph["nodes"][external_node] = {
                "PDGs": (1000,),
                "momenta": (momentum_name,),
                "indices": (),
                "vertex_id": -1,
                "edge_ids": (external_node,)
            }

            graph["edges"][external_node] = {}
            if type == "in":
                graph["edges"][external_node]["type"] = "in"
                graph["edges"][external_node]["vertices"] = (
                    external_node, internal_node)
            elif type == "out":
                graph["edges"][external_node]["type"] = "out"
                graph["edges"][external_node]["vertices"] = (
                    internal_node, external_node)
            else:
                raise ValueError("Unknown external edge type")

            graph["edges"][external_node]["name"] = momentum_name
            graph["edges"][external_node]["PDG"] = 1000
            graph["edges"][external_node]["momentum"] = ""
            graph["edges"][external_node]["indices"] = ()

            if internal_node not in graph["nodes"]:
                graph["nodes"][internal_node] = {
                    "PDGs": (1000,),
                    "momenta": (),
                    "indices": (),
                    "vertex_id": 0,
                    "edge_ids": (external_node,)
                }
            else:
                graph["nodes"][internal_node]["edge_ids"] += (
                    int(1000),)
                graph["nodes"][internal_node]["PDGs"] += (1000,)

        for virtual_edge_id, virtual_edge in enumerate(virtual_edges):
            pdg = virtual_edge[0]
            internal_node1 = virtual_edge[1]
            internal_node2 = virtual_edge[2]

            edge_id = virtual_edge_id + 1
            graph["edges"][edge_id] = {
                "name": "q" + str(edge_id),
                "PDG": pdg,
                "type": "virtual",
                "momentum": "",
                "indices": (),
                "vertices": (internal_node1, internal_node2)
            }

            for node in [internal_node1, internal_node2]:
                if node not in graph["nodes"]:
                    graph["nodes"][node] = {
                        "PDGs": (pdg,),
                        "momenta": (),
                        "indices": (),
                        "vertex_id": 0,
                        "edge_ids": (edge_id,)
                    }
                else:
                    graph["nodes"][node]["edge_ids"] += (edge_id,)
                    graph["nodes"][node]["PDGs"] += (pdg,)

        gammaloop_graph = Graph.from_qgraph(self.model, graph, args.name)
        self.amplitudes = cross_section.AmplitudeList([cross_section.Amplitude(
            args.name,
            [
                supergraph.AmplitudeGraph(
                    sg_id=0, sg_cut_id=0, fs_cut_id=0, amplitude_side=Side.LEFT,
                    multiplicity=args.multiplicity_factor,
                    graph=gammaloop_graph
                )
            ]
        )])

    # output command
    output_parser = ArgumentParser(prog='output')
    output_parser.add_argument(
        'output_path', type=str, help='Path to output the cross section to')
    output_parser.add_argument('-nev', '--no_evaluators', default=False, action='store_true',
                               help='Do not generate evaluators')
    output_parser.add_argument('-exp', '--expression', default=False, action='store_true',
                               help='Generate expression associated to the graph and output it to a text file')
    output_parser.add_argument('-mr', '--model_replacements', default=False, action='store_true',
                               help='Generate coupling replacements and output it to a text file')
    output_parser.add_argument('-ef', '--expression_format', type=str, default='file',
                               choices=['file', 'mathematica', 'latex'], help='Format to export symbolica objects in the numerator output.')
    output_parser.add_argument('-ow', '--overwrite_output', default=False, action='store_true',
                               help='Overwrite output if already existing.')

    output_parser.add_argument('-yo', '--yaml_only', default=False, action='store_true',
                               help='Only output yaml.')

    def do_output(self, str_args: str) -> None:
        if str_args == 'help':
            self.output_parser.print_help()
            return
        args = self.output_parser.parse_args(split_str_args(str_args))

        if Path(args.output_path).exists():
            if args.overwrite_output:
                shutil.rmtree(args.output_path)
            else:
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

            if args.expression:
                amplitude_exporter.export_expression(
                    args.output_path, self.amplitudes, args.expression_format)
            amplitude_exporter.export(
                args.output_path, self.amplitudes, args.no_evaluators)

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
    launch_parser.add_argument('--load_run_settings', '-lrs', action='store_true', default=False,
                               help='Load the run settings from the run_card.yaml in the process output.')

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
        if args.load_run_settings:
            run_settings = load_configuration(
                pjoin(args.path_to_launch, 'cards', 'run_card.yaml'), True)
            # Do not overwrite the external momenta setting in gammaLoop if they are set constant in the process dir
            if run_settings['Kinematics']['externals']['type'] == 'constant':
                del run_settings['Kinematics']
            self.config.update({'run_settings': run_settings})
        else:
            update_run_card_in_output(args.path_to_launch, self.config)

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
                args.path_to_launch)

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

    def sync_worker_with_output(self, no_sync=False) -> None:
        if no_sync or self.launched_output is None:
            return

        # Read metadata
        with open(pjoin(self.launched_output, 'output_metadata.yaml'), 'r', encoding='utf-8') as file:
            output_metadata = OutputMetaData.from_yaml_str(file.read())

        # Sync the model with the one used in the output
        with open(pjoin(self.launched_output, 'sources', 'model', f"{output_metadata['model_name']}.yaml"), 'r', encoding='utf-8') as file:
            self.model = Model.from_yaml(file.read())
        self.model.apply_input_param_card(InputParamCard.from_param_card(
            ParamCard(pjoin(self.launched_output, 'cards', 'param_card.dat')), self.model), simplify=False)
        processed_yaml_model = self.model.to_yaml()
        self.rust_worker.load_model_from_yaml_str(processed_yaml_model)

        run_settings = load_configuration(
            pjoin(self.launched_output, 'cards', 'run_card.yaml'), True)
        # Do not overwrite the external momenta setting in gammaLoop if they are set constant in the process dir
        if run_settings['Kinematics']['externals']['type'] == 'constant':
            del run_settings['Kinematics']
        self.config.update({'run_settings': run_settings})

        # Depending on the type of output, sync cross-section or amplitude
        if output_metadata['output_type'] == 'amplitudes':
            self.amplitudes = cross_section.AmplitudeList()
            self.rust_worker.reset_amplitudes()
            for amplitude_name in output_metadata['contents']:
                with open(pjoin(self.launched_output, 'sources', 'amplitudes', f'{amplitude_name}', 'amplitude.yaml'), 'r', encoding='utf-8') as file:
                    amplitude_yaml = file.read()
                    self.amplitudes.add_amplitude(
                        cross_section.Amplitude.from_yaml_str(self.model, amplitude_yaml))
                    self.rust_worker.add_amplitude_from_yaml_str(
                        amplitude_yaml)
            self.rust_worker.load_amplitudes_derived_data(
                pjoin(self.launched_output))

            self.rust_worker.load_amplitude_integrands(
                pjoin(self.launched_output, 'cards', 'run_card.yaml'))

        elif output_metadata['output_type'] == 'cross_sections':
            self.cross_sections = cross_section.CrossSectionList()
            self.rust_worker.reset_cross_sections()
            for cross_section_name in output_metadata['contents']:
                with open(pjoin(self.launched_output, 'sources', 'cross_sections', f'{cross_section_name}', 'cross_section.yaml'), 'r', encoding='utf-8') as file:
                    cross_section_yaml = file.read()
                    self.cross_sections.add_cross_section(
                        cross_section.CrossSection.from_yaml_str(self.model, cross_section_yaml))
                    self.rust_worker.add_cross_section_from_yaml_str(
                        cross_section_yaml)

        self.rust_worker.sync()

    # inspect command
    inspect_parser = ArgumentParser(prog='inspect')
    inspect_parser.add_argument(
        "integrand", type=str, help="Integrand to inspect.")
    inspect_parser.add_argument('--use-f128', '-f128', action='store_true',
                                default=False, help='Use f128 precision for the inspection.')
    inspect_parser.add_argument('--point', '-p', nargs="+", type=float,
                                default=[], help='Point to inspect.')
    inspect_parser.add_argument(
        '--term', '-t', nargs="+", type=int, default=tuple([0,]), help="term to inspect.")
    inspect_parser.add_argument('--is-momentum-space', '-ms', action='store_true',
                                default=False, help='Inspect in momentum space.')
    inspect_parser.add_argument(
        '--force-radius', '-fr', action='store_true', default=False, help='Force radius to be used.')
    inspect_parser.add_argument(
        '--no_sync', '-ns', action='store_true', default=False,
        help='Do not sync rust worker with the process output (safe to do if not config change was issued since launch).')
    inspect_parser.add_argument('--last_max_weight', '-lmw', action='store_true',
                                default=False, help='Inspect the max weight point of the previous run')

    def do_inspect(self, str_args: str) -> Dict[str, Any]:
        if str_args == 'help':
            self.inspect_parser.print_help()
            return {}

        args = self.inspect_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        self.sync_worker_with_output(args.no_sync)

        if args.last_max_weight:
            workspace_path = self.launched_output.joinpath("workspace")
            res: tuple[float, float] = self.rust_worker.inspect_lmw_integrand(
                args.integrand, str(workspace_path), args.use_f128)
        else:
            res: tuple[float, float] = self.rust_worker.inspect_integrand(
                args.integrand, args.point, args.term, args.force_radius, args.is_momentum_space, args.use_f128)

        log_res: dict[str, Any] = {}
        log_res['final_result'] = complex(res[0], res[1])

        for file in os.listdir('log.glog'):
            file_location = pjoin('log.glog', file)
            file_dict = debug_display.parse_log_impl(file_location)
            log_res[file] = file_dict

        return log_res

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
    integrate_parser.add_argument(
        '--no_sync', '-ns', action='store_true', default=False,
        help='Do not sync rust worker with the process output (safe to do if not config change was issued since launch).')

    def do_integrate(self, str_args: str) -> list[complex]:
        if str_args == 'help':
            self.integrate_parser.print_help()
            return []
        args = self.integrate_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        self.sync_worker_with_output(args.no_sync)

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

        res = self.rust_worker.integrate_integrand(
            args.integrand, args.cores, str(result_output_path), str(workspace_path), target)

        return [complex(r[0], r[1]) for r in res]

        # nuke the workspace if integration finishes
        # For now leave the possibility of restarting where integration left off.
        # Maybe in the future add an option to automatically clean the workspace after running is completed or
        # specify a "run_tag" that allows to have mutliple workspace concurrently active
        # shutil.rmtree(workspace_path)

    # test_ir_limits
    test_ir_limits_parser = ArgumentParser(prog='test_ir_limits')
    test_ir_limits_parser.add_argument(
        '--no_sync', '-ns', action='store_true', default=False,
        help='Do not sync rust worker with the process output (safe to do if not config change was issued since launch).')

    def do_test_ir_limits(self, str_args: str) -> None:
        if str_args == 'help':
            self.test_ir_limits_parser.print_help()
            return
        args = self.test_ir_limits_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        self.sync_worker_with_output(args.no_sync)

        raise GammaLoopError("Command not implemented yet")

    # test_uv_limits
    test_uv_limits_parser = ArgumentParser(prog='test_ir_limits')
    test_uv_limits_parser.add_argument(
        '--no_sync', '-ns', action='store_true', default=False,
        help='Do not sync rust worker with the process output (safe to do if not config change was issued since launch).')

    def do_test_uv_limits(self, str_args: str) -> None:
        if str_args == 'help':
            self.test_uv_limits_parser.print_help()
            return
        args = self.test_uv_limits_parser.parse_args(split_str_args(str_args))

        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        self.sync_worker_with_output(args.no_sync)

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
    hpc_parser.add_argument(
        '--no_sync', '-ns', action='store_true', default=False,
        help='Do not sync rust worker with the process output (safe to do if not config change was issued since launch).')

    def do_hpc_run(self, str_args: str) -> None:
        args = self.hpc_parser.parse_args(split_str_args(str_args))
        if self.launched_output is None:
            raise GammaLoopError(
                "No output launched. Please launch an output first with 'launch' command.")

        self.sync_worker_with_output(args.no_sync)

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

    log_parser = ArgumentParser(prog='display_debug_log')
    log_parser.add_argument('--log_file', '-lf', type=str,
                            help='Log file to display', default="log.glog")
    log_parser.add_argument(
        '-eval', '-e', type=str, default=None)
    log_parser.add_argument('--subtraction', '-s',  action='store_true',
                            default=False, help='Show subtraction debug info')

    def do_display_debug_log(self, str_args: str) -> None:
        args = self.log_parser.parse_args(split_str_args(str_args))

        if args.log_file is None:
            raise GammaLoopError(
                "No log file to display, please provide a file using -lf")

        general_debug_dict = debug_display.build_general_debug_dict(
            args.log_file)
        debug_display.display_general(general_debug_dict)

        if args.eval is not None:
            tmp = eval(args.eval)
            file_list = ["{}_{}.jsonl".format(
                tmp_elem[0], tmp_elem[1]) for tmp_elem in tmp]
        else:
            file_list = os.listdir(args.log_file)

        for file in file_list:
            if file == "general.jsonl":
                continue

            eval_dict = debug_display.build_eval_debug_dict(
                args.log_file, file)

            logger.info("Debug info for for rotation '%s%s%s'",
                        Colour.BLUE, file, Colour.END, )

            debug_display.display_eval_default(eval_dict)

            if args.subtraction:
                logger.info("")
                logger.info("subtraction: ")
                logger.info("")
                debug_display.display_subtraction_data(eval_dict)
