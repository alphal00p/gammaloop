from __future__ import annotations
import argparse
import os
from pathlib import Path
import shutil
import yaml

import gammaloop
import gammaloop.cross_section as cross_section
from gammaloop.misc.common import DATA_PATH, pjoin
import gammaloop.misc.utils as utils


class OutputMetaData(dict):
    def __init__(self, *args, **opts):
        super(OutputMetaData, self).__init__(*args, **opts)

    def to_yaml_str(self):
        return utils.verbose_yaml_dump(dict(self))

    @staticmethod
    def from_yaml_str(yaml_str: str):
        return OutputMetaData(yaml.safe_load(yaml_str))


class GammaLoopExporter(object):

    def __init__(self, gammaloop_interface: gammaloop.gamma_loop_interface.GammaLoop, _args: argparse.Namespace):
        self.gammaloop: gammaloop.gamma_loop_interface.GammaLoop = gammaloop_interface
        # Process args argument further here if needed

    def generic_export(self, export_root: Path):

        os.makedirs(export_root)

        # Build the output structure
        os.makedirs(pjoin(export_root, 'Cards'))
        with open(pjoin(export_root, 'Cards', 'param_card.yaml'), 'w', encoding='utf-8') as file:
            file.write("TODO")
        shutil.copy(pjoin(DATA_PATH, 'run_cards', 'rust_run_config.yaml'), pjoin(
            export_root, 'Cards', 'run_card.yaml'))
        with open(pjoin(export_root, 'Cards', 'proc_card.gL'), 'w', encoding='utf-8') as file:
            file.write(self.gammaloop.command_history.nice_string())

        os.makedirs(pjoin(export_root, 'sources'))
        os.makedirs(pjoin(export_root, 'sources', 'model'))
        with open(pjoin(export_root, 'sources', 'model', f'{self.gammaloop.model.name}.yaml'), 'w', encoding='utf-8') as file:
            file.write(self.gammaloop.model.to_yaml())

        os.makedirs(pjoin(export_root, 'runs'))


class AmplitudesExporter(GammaLoopExporter):

    def __init__(self, gammaloop_interface: gammaloop.gamma_loop_interface.GammaLoop, args: argparse.Namespace):
        GammaLoopExporter.__init__(self, gammaloop_interface, args)
        # Further processing here for additional info if need be

    def export(self, export_root: Path, amplitudes: cross_section.AmplitudeList):

        super(AmplitudesExporter, self).generic_export(export_root)

        output_data = OutputMetaData({
            'model_name': self.gammaloop.model.name,
            'output_type': 'amplitudes',
            'amplitudes': [amplitude.name for amplitude in amplitudes]
        })
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
                export_root, 'sources', 'amplitudes', f'{amplitude.name}', 'drawings.pdf')
            amplitude.draw(self.gammaloop.model, drawings_path)
            self.gammaloop.rust_worker.add_amplitude_from_yaml_str(
                amplitude_yaml)

        # Now address the rust export aspect
        self.gammaloop.rust_worker.export_amplitudes(
            str(export_root), [amp.name for amp in amplitudes])


class CrossSectionsExporter(GammaLoopExporter):

    def __init__(self, gammaloop_interface: gammaloop.gamma_loop_interface.GammaLoop, args: argparse.Namespace):
        GammaLoopExporter.__init__(self, gammaloop_interface, args)
        # Further processing here for additional info if need be

    def export(self, export_root: Path, cross_sections: cross_section.CrossSectionList):

        super(CrossSectionsExporter, self).generic_export(export_root)

        output_data = OutputMetaData({
            'model_name': self.gammaloop.model.name,
            'output_type': 'cross_sections',
            'cross_sections': [cross_section.name for cross_section in cross_sections]
        })
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
                export_root, 'sources', 'cross_sections', f'{one_cross_section.name}', 'drawings.pdf')
            one_cross_section.draw(self.gammaloop.model, drawings_path)
            self.gammaloop.rust_worker.add_cross_section_from_yaml_str(yaml_xs)

        # Now address the rust export aspect
        self.gammaloop.rust_worker.export_cross_sections(
            str(export_root), [cs.name for cs in cross_sections])
