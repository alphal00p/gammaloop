from __future__ import annotations
import argparse
import os
from pathlib import Path
import shutil
import yaml
import copy
import math
from typing import Any, TYPE_CHECKING

from gammaloop.cross_section.cross_section import AmplitudeList, CrossSectionList
from gammaloop.base_objects.graph import Vertex, EdgeType
from gammaloop.misc.common import DATA_PATH, pjoin, GammaLoopError, logger, load_configuration
import gammaloop.misc.utils as utils
import gammaloop.interface.gammaloop_interface as gammaloop_interface

if TYPE_CHECKING:
    from gammaloop.interface.gammaloop_interface import GammaLoop


class OutputMetaData(dict[str, Any]):
    def __init__(self, *args: Any, **opts: Any):
        super(OutputMetaData, self).__init__(*args, **opts)

    def to_yaml_str(self):
        return utils.verbose_yaml_dump(dict(self))

    @staticmethod
    def from_yaml_str(yaml_str: str):
        return OutputMetaData(yaml.safe_load(yaml_str))


def update_run_card_in_output(process_dir: Path, settings: gammaloop_interface.GammaLoopConfiguration):
    process_config = gammaloop_interface.GammaLoopConfiguration(path='')
    process_config.update(
        {'run_settings': load_configuration(pjoin(process_dir, 'cards', 'run_card.yaml'), True)})

    run_config = copy.deepcopy(settings.get_setting('run_settings'))
    # Do not update settings meant to be automatically adjusted if not explicitly set
    if run_config['Kinematics']['externals']['type'] == 'constant' and run_config['Kinematics']['externals']['momenta'] is None:
        del run_config['Kinematics']
    process_config.update({'run_settings': run_config})

    with open(pjoin(process_dir, 'cards', 'run_card.yaml'), 'w', encoding='utf-8') as file:
        file.write(utils.verbose_yaml_dump(
            process_config.get_setting('run_settings')))


def split_str_args(str_args: str) -> list[str]:
    return str_args.split(' ') if str_args != '' else []


class GammaLoopExporter(object):

    def __init__(self, gammaloop_interface: GammaLoop, _args: argparse.Namespace):
        self.gammaloop: GammaLoop = gammaloop_interface
        self.configuration_for_process = copy.deepcopy(self.gammaloop.config)
        # Process args argument further here if needed

    def generic_export(self, export_root: Path):

        os.makedirs(export_root, exist_ok=True)

        # Build the output structure
        os.makedirs(pjoin(export_root, 'cards'))
        shutil.copy(pjoin(self.gammaloop.model_directory, self.gammaloop.model.name,
                          f"restrict_{'full' if self.gammaloop.model.restriction is None else self.gammaloop.model.restriction}.dat"),
                    pjoin(export_root, 'cards', 'param_card.dat'))
        with open(pjoin(export_root, 'cards', 'run_card.yaml'), 'w', encoding='utf-8') as file:
            file.write(utils.verbose_yaml_dump(
                self.configuration_for_process.get_setting('run_settings')))
        with open(pjoin(export_root, 'cards', 'proc_card.gL'), 'w', encoding='utf-8') as file:
            file.write(self.gammaloop.command_history.nice_string())

        os.makedirs(pjoin(export_root, 'sources'), exist_ok=True)
        os.makedirs(pjoin(export_root, 'sources', 'model'), exist_ok=True)
        with open(pjoin(export_root, 'sources', 'model', f'{self.gammaloop.model.name}.yaml'), 'w', encoding='utf-8') as file:
            file.write(self.gammaloop.model.to_yaml())

        os.makedirs(pjoin(export_root, 'runs'))

        return OutputMetaData({
            'model_name': self.gammaloop.model.name,
        })

    def finalize_drawing(self, drawings_path: Path, drawing_file_paths: list[Path]):

        shutil.copy(
            pjoin(DATA_PATH, 'templates', 'drawing', 'combine_pages.py'),
            pjoin(drawings_path, 'combine_pages.py'))

        with open(pjoin(DATA_PATH, 'templates', 'drawing', 'makefile'), 'r', encoding='utf-8') as makefile_in:
            with open(pjoin(drawings_path, 'makefile'), 'w', encoding='utf-8') as makefile_out:
                all_targets: list[str] = []
                for drawing_file_path in drawing_file_paths:
                    if drawing_file_path.suffix != '.tex':
                        raise GammaLoopError(
                            "Finalization of diagram drawings only supports latex format.")
                    drawing_name = pjoin(os.path.relpath(
                        drawing_file_path.parent, drawings_path), drawing_file_path.stem)
                    all_targets.append(f'{drawing_name}.pdf')
                makefile_out.write(makefile_in.read().format(
                    all_targets=' '.join(all_targets),
                    output_name='feynman_diagrams.pdf',
                    n_rows=self.gammaloop.config['drawing']['combined_graphs_pdf_grid_shape'][0],
                    n_columns=self.gammaloop.config['drawing']['combined_graphs_pdf_grid_shape'][1]
                ))

    def build_external_momenta_from_connections(self, external_connections: list[tuple[Vertex | None, Vertex | None]], e_cm: float) -> tuple[float, list[list[float]]]:

        # (orientation_sign, mass, direction [in/out] )
        conf: list[tuple[int, float, int]] = []
        for conn in external_connections:
            match conn:
                case (None, v2) if v2 is not None:
                    sign = 1
                    direction = 1
                    match v2.edges[0].edge_type:
                        case EdgeType.INCOMING:
                            sign = -1
                            direction = 1
                        case EdgeType.OUTGOING:
                            sign = 1
                            direction = -1
                        case _:
                            raise GammaLoopError(
                                "Invalid external connection.")
                    mass = self.gammaloop.model.get_parameter(
                        v2.vertex_info.get_particles()[0].mass.name).value
                    if mass is None:
                        raise GammaLoopError(
                            "Explicit default value of the mass of external particle not defined.")
                    conf.append((sign, mass.real, direction))
                case (v1, _) if v1 is not None:
                    sign = 1
                    direction = 1
                    match v1.edges[0].edge_type:
                        case EdgeType.INCOMING:
                            sign = 1
                            direction = 1
                        case EdgeType.OUTGOING:
                            sign = -1
                            direction = -1
                        case _:
                            raise GammaLoopError(
                                "Invalid external connection.")
                    mass = self.gammaloop.model.get_parameter(
                        v1.vertex_info.get_particles()[0].mass.name).value
                    if mass is None:
                        raise GammaLoopError(
                            "Explicit default value of the mass of external particle not defined.")
                    conf.append((sign, mass.real, direction))
                case _:
                    raise GammaLoopError("Invalid external connection.")
        if len(conf) == 0:
            raise GammaLoopError("Invalid external connection.")

        # ensure some phase-space available
        if e_cm <= sum(m for (_es, m, _dir) in conf):
            e_cm = 2.*sum(m for (_es, m, _dir) in conf)

        externals = []
        if len(conf) == 1:
            if conf[0][1] == 0.:
                externals.append([conf[0][0]*conf[0][2]*e_cm, 0., 0., 0.])
            else:
                externals.append(
                    [conf[0][0]*conf[0][2]*conf[0][1], 0., 0., 0.])

        # We should do something smarter here (use an actual phase-space generator), but these defaults are not too important for now
        elif len(conf) == 2:
            if conf[0][1] == 0. and conf[1][1] == 0.:
                externals.extend([
                    [conf[0][0]*conf[0][2]*e_cm/2., 0.,
                        0., conf[0][0]*conf[0][2]*e_cm/2.],
                    [conf[1][0]*conf[1][2]*e_cm/2., 0.,
                        0., -conf[1][0]*conf[1][2]*e_cm/2.],
                ])
            else:
                externals.extend([
                    [
                        conf[0][0]*conf[0][2] *
                        math.sqrt((e_cm/2.)**2+conf[0][1]**2),
                        0., 0., conf[0][0]*conf[0][2]*e_cm/2.
                    ],
                    [
                        conf[1][0]*conf[0][2] *
                        math.sqrt((e_cm/2.)**2+conf[1][1]**2),
                        0., 0., -conf[1][0]*conf[0][2]*e_cm/2.
                    ],
                ])
        else:
            target_kin_energy = (
                e_cm-sum(m for (_es, m, _dir) in conf))/float(len(conf))
            first_leg_momentum = [0., 0., 0., 0.]
            for i_ext, (orientation_sign, mass, direction_sign) in enumerate(conf[1:]):
                spatial_part = [
                    1.+3*i_ext,
                    2.+3*i_ext,
                    3.+3*i_ext,
                ]
                spatial_part_sq = sum(x**2 for x in spatial_part)
                rescaling = math.sqrt(
                    (2.*target_kin_energy*mass + target_kin_energy**2)/spatial_part_sq)
                leg_momentum = [
                    mass+target_kin_energy,
                    rescaling*spatial_part[0],
                    rescaling*spatial_part[1],
                    rescaling*spatial_part[2],
                ]
                for j in range(4):
                    first_leg_momentum[j] -= direction_sign*leg_momentum[j]

                externals.append([orientation_sign*pi for pi in leg_momentum])

            # note: this last leg will not be onshell!
            externals.insert(0,
                             [conf[0][0]*conf[0][2] * pi for pi in first_leg_momentum])

        return (e_cm, externals)


class AmplitudesExporter(GammaLoopExporter):

    def __init__(self, gammaloop_interface: gammaloop_interface.GammaLoop, args: argparse.Namespace):
        GammaLoopExporter.__init__(self, gammaloop_interface, args)
        # Further processing here for additional info if need be

    def adjust_run_settings(self, amplitudes: AmplitudeList):

        if self.configuration_for_process.get_setting('run_settings.Kinematics.externals.type') == 'constant' and \
                self.configuration_for_process.get_setting('run_settings.Kinematics.externals.momenta') is None:
            if len(amplitudes) == 0 or len(amplitudes[0].amplitude_graphs) == 0:
                logger.warning(
                    "Could not identify external momenta structure.")
                return
            e_cm, external_momenta = self.build_external_momenta_from_connections(
                amplitudes[0].amplitude_graphs[0].graph.external_connections,
                self.configuration_for_process.get_setting(
                    'run_settings.Kinematics.e_cm')
            )
            self.configuration_for_process.set_setting(
                'run_settings.Kinematics.externals.momenta', external_momenta)
            self.configuration_for_process.set_setting(
                'run_settings.Kinematics.e_cm', e_cm)

    def export_expression(self, export_root: Path, amplitudes: AmplitudeList, format: str):
        for amplitude in amplitudes:
            os.makedirs(pjoin(export_root, 'sources',
                        'amplitudes', f'{amplitude.name}', 'expressions'))

        self.gammaloop.rust_worker.export_expressions(str(export_root), format)
    


    def export(self, export_root: Path, amplitudes: AmplitudeList):

        # Tweak the run configuration for the particular process exported before attending to the generic export
        self.adjust_run_settings(amplitudes)

        output_data = super(AmplitudesExporter,
                            self).generic_export(export_root)
        output_data['output_type'] = 'amplitudes'
        output_data['contents'] = [amplitude.name for amplitude in amplitudes]

        with open(pjoin(export_root, 'output_metadata.yaml'), 'w', encoding='utf-8') as file:
            file.write(output_data.to_yaml_str())

        os.makedirs(pjoin(export_root, 'sources', 'amplitudes'))
        self.gammaloop.rust_worker.reset_cross_sections()
        for amplitude in amplitudes:
            os.makedirs(pjoin(export_root, 'sources',
                        'amplitudes', f'{amplitude.name}'))
            amplitude_yaml = amplitude.to_yaml_str()
            # Already writing the file below is not necessary as it will be overwritten by the rust export, but it is useful for debugging
            with open(pjoin(export_root, 'sources', 'amplitudes', f'{amplitude.name}', 'amplitude.yaml'), 'w', encoding='utf-8') as file:
                file.write(amplitude_yaml)
            drawings_path = pjoin(
                export_root, 'sources', 'amplitudes', f'{amplitude.name}', 'drawings')
            os.makedirs(drawings_path)
            drawing_file_paths = amplitude.draw(
                self.gammaloop.model, drawings_path, **self.gammaloop.config['drawing'])
            self.finalize_drawing(Path(drawings_path), drawing_file_paths)

            self.gammaloop.rust_worker.add_amplitude_from_yaml_str(
                amplitude_yaml)

        # Now address the rust export aspect
        self.gammaloop.rust_worker.export_amplitudes(
            str(export_root), [amp.name for amp in amplitudes])

        if self.gammaloop.config['export_settings']['write_default_settings']:
            self.gammaloop.rust_worker.write_default_settings(str(export_root))


class CrossSectionsExporter(GammaLoopExporter):

    def __init__(self, gammaloop_interface: gammaloop_interface.GammaLoop, args: argparse.Namespace):
        GammaLoopExporter.__init__(self, gammaloop_interface, args)
        # Further processing here for additional info if need be

    def export(self, export_root: Path, cross_sections: CrossSectionList):

        output_data = super(CrossSectionsExporter,
                            self).generic_export(export_root)
        output_data['output_type'] = 'cross_sections'
        output_data['contents'] = [
            cross_section.name for cross_section in cross_sections]

        with open(pjoin(export_root, 'output_metadata.yaml'), 'w', encoding='utf-8') as file:
            file.write(output_data.to_yaml_str())

        os.makedirs(pjoin(export_root, 'sources', 'cross_sections'))
        self.gammaloop.rust_worker.reset_cross_sections()
        for one_cross_section in cross_sections:
            os.makedirs(pjoin(export_root, 'sources',
                        'cross_sections', f'{one_cross_section.name}'))
            yaml_xs = one_cross_section.to_yaml_str()
            # Already writing the file below is not necessary as it will be overwritten by the rust export, but it is useful for debugging
            with open(pjoin(export_root, 'sources', 'cross_sections', f'{one_cross_section.name}', 'cross_section.yaml'), 'w', encoding='utf-8') as file:
                file.write(yaml_xs)
            drawings_path = pjoin(
                export_root, 'sources', 'cross_sections', f'{one_cross_section.name}', 'drawings')
            os.makedirs(drawings_path)
            drawing_file_paths = one_cross_section.draw(
                self.gammaloop.model, drawings_path, **self.gammaloop.config['drawing'])
            self.finalize_drawing(Path(drawings_path), drawing_file_paths)

            self.gammaloop.rust_worker.add_cross_section_from_yaml_str(yaml_xs)

        # Now address the rust export aspect
        self.gammaloop.rust_worker.export_cross_sections(
            str(export_root), [cs.name for cs in cross_sections])
