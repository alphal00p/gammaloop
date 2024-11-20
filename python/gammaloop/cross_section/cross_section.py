from __future__ import annotations

from pathlib import Path
from typing import Any
import yaml
import os
import gammaloop.cross_section.supergraph as supergraph
from gammaloop.base_objects.model import Model
import gammaloop.misc.utils as utils


class CrossSection(object):

    # This class is just a stub for now, but it will be used to store the cross section information, incl. process definition.
    def __init__(self, name: str, supergraphs: list[supergraph.SuperGraph]):
        self.name: str = name
        self.supergraphs: list[supergraph.SuperGraph] = supergraphs

    def draw(self, model: Model, drawings_path: str, **drawing_options: Any) -> dict[str, list[Path | None]]:

        if len(self.supergraphs) == 0:
            return {mode: [] for mode in drawing_options["modes"]}

        drawing_file_paths_per_mode: dict[str, list[Path | None]] = {}
        for mode in drawing_options["modes"]:
            drawings_path_for_mode = os.path.join(drawings_path, mode)
            os.makedirs(drawings_path_for_mode, exist_ok=True)
            drawing_file_paths_per_mode[mode] = []
            for super_graph in self.supergraphs:
                drawing_file_paths_per_mode[mode].append(super_graph.draw(
                    mode, model, drawings_path_for_mode, None, **drawing_options))

        return drawing_file_paths_per_mode

    @staticmethod
    def from_serializable_dict(model: Model, cross_section_dict: dict[str, Any]) -> CrossSection:

        return CrossSection(
            cross_section_dict['name'],
            [supergraph.SuperGraph.from_serializable_dict(
                model, sg_dict) for sg_dict in cross_section_dict['supergraphs']]
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'name': self.name,
            'supergraphs': [sg.to_serializable_dict() for sg in self.supergraphs],
        }

    @staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> CrossSection:
        cross_section_dict = yaml.safe_load(yaml_str)
        return CrossSection.from_serializable_dict(model, cross_section_dict)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump(self.to_serializable_dict())


class Amplitude(object):

    def __init__(self, name: str, amplitude_graphs: list[supergraph.AmplitudeGraph]):
        self.name: str = name
        self.amplitude_graphs: list[supergraph.AmplitudeGraph] = amplitude_graphs

    def draw(self, model: Model, drawings_path: str, **drawing_options: Any) -> dict[str, list[Path | None]]:

        if len(self.amplitude_graphs) == 0:
            return {mode: [] for mode in drawing_options["modes"]}

        drawing_file_paths_per_mode: dict[str, list[Path | None]] = {}
        for mode in drawing_options["modes"]:
            drawings_path_for_mode = os.path.join(drawings_path, mode)
            os.makedirs(drawings_path_for_mode, exist_ok=True)
            drawing_file_paths_per_mode[mode] = []
            for amplitude_graph in self.amplitude_graphs:
                drawing_file_paths_per_mode[mode].append(
                    amplitude_graph.draw(mode, model, drawings_path_for_mode, f'{amplitude_graph.fs_cut_id}_{amplitude_graph.graph.name}', **drawing_options))

        return drawing_file_paths_per_mode

    @staticmethod
    def from_serializable_dict(model: Model, amplitude_dict: dict[str, Any]) -> Amplitude:
        return Amplitude(
            amplitude_dict['name'],
            [supergraph.AmplitudeGraph.from_serializable_dict(
                model, amp_graph_dict) for amp_graph_dict in amplitude_dict['amplitude_graphs']]
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'name': self.name,
            'amplitude_graphs': [amp_graph.to_serializable_dict() for amp_graph in self.amplitude_graphs],
        }

    @staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> Amplitude:
        amplitude_dict = yaml.safe_load(yaml_str)
        return Amplitude.from_serializable_dict(model, amplitude_dict)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump(self.to_serializable_dict())


class CrossSectionList(list[CrossSection]):

    def __init__(self, *args: Any):
        super(CrossSectionList, self).__init__(*args)

    @staticmethod
    def from_serializable(model: Model, cross_section_list: list[dict[str, Any]]) -> CrossSectionList:
        return CrossSectionList([CrossSection.from_serializable_dict(model, cross_section_dict) for cross_section_dict in cross_section_list])

    def to_serializable(self) -> list[dict[str, Any]]:
        return [cross_section.to_serializable_dict() for cross_section in self]

    @staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> CrossSectionList:
        cross_section_list = yaml.safe_load(yaml_str)
        return CrossSectionList.from_serializable(model, cross_section_list)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump([cs.to_serializable_dict() for cs in self])

    def add_cross_section(self, cross_section: CrossSection) -> None:
        self.append(cross_section)


class AmplitudeList(list[Amplitude]):

    def __init__(self, *args: list[Any]):
        super(AmplitudeList, self).__init__(*args)

    @staticmethod
    def from_serializable(model: Model, amplitude_list: list[dict[str, Any]]) -> AmplitudeList:
        return AmplitudeList([Amplitude.from_serializable_dict(model, amplitude_dict) for amplitude_dict in amplitude_list])

    def to_serializable(self) -> list[dict[str, Any]]:
        return [amplitude.to_serializable_dict() for amplitude in self]

    @staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> AmplitudeList:
        amplitude_list = yaml.safe_load(yaml_str)
        return AmplitudeList.from_serializable(model, amplitude_list)

    def add_amplitude(self, amplitude: Amplitude) -> None:
        self.append(amplitude)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump([a.to_serializable_dict() for a in self])
