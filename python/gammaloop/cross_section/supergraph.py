from __future__ import annotations

from typing import Any
from pathlib import Path
from gammaloop.base_objects.model import Model
from gammaloop.base_objects.graph import Graph, Edge
from gammaloop.misc.common import Side, pjoin
import gammaloop.cross_section.cross_section as cross_section


class SuperGraphCut(object):
    def __init__(self, cut_edges: list[Edge], forward_scattering_graph: ForwardScatteringGraph):
        self.cut_edges: list[Edge] = cut_edges
        self.forward_scattering_graph: ForwardScatteringGraph = forward_scattering_graph

    @staticmethod
    def from_serializable_dict(model: Model, graph: Graph, sg_dict: dict[str, Any]) -> SuperGraphCut:

        return SuperGraphCut(
            [graph.get_edge(edge_name) for edge_name in sg_dict['cut_edges']],
            ForwardScatteringGraph.from_serializable_dict(
                model, sg_dict['forward_scattering_graph'])
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'cut_edges': [edge.name for edge in self.cut_edges],
            'forward_scattering_graph': self.forward_scattering_graph.to_serializable_dict(),
        }


class SuperGraph(object):
    def __init__(self, sg_id: int, graph: Graph, multiplicity: str, topology_class: list[int], cuts: list[SuperGraphCut]):
        self.sg_id: int = sg_id
        self.graph: Graph = graph
        self.multiplicity: str = multiplicity
        self.topology_class: list[int] = topology_class
        self.cuts: list[SuperGraphCut] = cuts

    def draw(self, drawing_mode: str | None, model: Model, drawings_path: str, file_name: str | None = None, **drawing_options: dict[str, Any]) -> Path | None:
        if len(self.graph.edges) == 0 and len(self.cuts) > 0:
            return self.cuts[0].forward_scattering_graph.draw(drawing_mode, model, drawings_path, file_name, specify_cut=False, **drawing_options)

        if file_name is None:
            file_name = f'{self.sg_id}_{self.graph.name}'
        return self.graph.draw(drawing_mode, model, pjoin(drawings_path, file_name), caption=f"Supergraph {self.graph.name.replace('_', ' ')}", diagram_id='#%d' % self.sg_id, **drawing_options)

    @ staticmethod
    def from_serializable_dict(model: Model, sg_dict: dict[str, Any]) -> SuperGraph:

        graph = Graph.from_serializable_dict(model, sg_dict['graph'])
        return SuperGraph(
            sg_dict['sg_id'],
            graph,
            sg_dict['multiplicity'],
            sg_dict['topology_class'],
            [SuperGraphCut.from_serializable_dict(
                model, graph, cut_dict) for cut_dict in sg_dict['cuts']]
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'sg_id': self.sg_id,
            'graph': self.graph.to_serializable_dict(),
            'multiplicity': self.multiplicity,
            'topology_class': self.topology_class,
            'cuts': [cut.to_serializable_dict() for cut in self.cuts],
        }


class ForwardScatteringGraphCut(object):
    def __init__(self, cut_edges: list[Edge], amplitudes: tuple[cross_section.Amplitude, cross_section.Amplitude]):
        self.cut_edges: list[Edge] = cut_edges
        self.amplitudes: tuple[cross_section.Amplitude,
                               cross_section.Amplitude] = amplitudes

    @ staticmethod
    def from_serializable_dict(graph: Graph, model: Model, fsgc_dict: dict[str, Any]) -> ForwardScatteringGraphCut:

        return ForwardScatteringGraphCut(
            [graph.get_edge(edge_name)
             for edge_name in fsgc_dict['cut_edges']],
            (cross_section.Amplitude.from_serializable_dict(model, fsgc_dict['amplitudes'][0]), cross_section.Amplitude.from_serializable_dict(
                model, fsgc_dict['amplitudes'][1]))
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'cut_edges': [edge.name for edge in self.cut_edges],
            'amplitudes': [amplitude.to_serializable_dict() for amplitude in self.amplitudes],
        }


class ForwardScatteringGraph(object):
    def __init__(self, sg_id: int, sg_cut_id: int, graph: Graph, multiplicity: str, cuts: list[ForwardScatteringGraphCut]):
        self.sg_id: int = sg_id
        self.sg_cut_id: int = sg_cut_id
        self.graph: Graph = graph
        self.multiplicity: str = multiplicity
        self.cuts: list[ForwardScatteringGraphCut] = cuts

    def draw(self, drawing_mode: str | None, model: Model, drawings_path: str, file_name: str | None, specify_cut: bool = True, **drawing_options: dict[str, Any]) -> Path | None:

        if file_name is None:
            file_name = f'{self.sg_id}_{self.sg_cut_id}_{self.graph.name}'

        if specify_cut:
            g_id = f'#(sg={self.sg_id},sg_cut={self.sg_cut_id})'
        else:
            g_id = f'#{self.sg_id}'
        return self.graph.draw(drawing_mode, model, pjoin(drawings_path, file_name),
                               caption=f"FSG {self.graph.name.replace(f'_{self.sg_id}', '').replace('_', ' ')}", diagram_id=g_id, **drawing_options)

    @staticmethod
    def from_serializable_dict(model: Model, amplitude_graph_dict: dict[str, Any]) -> ForwardScatteringGraph:
        graph = Graph.from_serializable_dict(
            model, amplitude_graph_dict['graph'])
        return ForwardScatteringGraph(
            amplitude_graph_dict['sg_id'],
            amplitude_graph_dict['sg_cut_id'],
            graph,
            amplitude_graph_dict['multiplicity'],
            [ForwardScatteringGraphCut.from_serializable_dict(
                graph, model, cut_dict) for cut_dict in amplitude_graph_dict['cuts']]
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'sg_id': self.sg_id,
            'sg_cut_id': self.sg_cut_id,
            'graph': self.graph.to_serializable_dict(),
            'multiplicity': self.multiplicity,
            'cuts': [cut.to_serializable_dict() for cut in self.cuts],
        }


class AmplitudeGraph(object):
    def __init__(self, sg_id: int, sg_cut_id: int, fs_cut_id: int, amplitude_side: Side, graph: Graph, multiplicity: str, multi_channeling_structure: list[int] = []):
        self.sg_id: int = sg_id
        self.sg_cut_id: int = sg_cut_id
        self.fs_cut_id: int = fs_cut_id
        self.amplitude_side: Side = amplitude_side
        self.graph: Graph = graph
        self.multiplicity: str = multiplicity
        self.multi_channeling_channels: list[int] = []

    def draw(self, drawing_mode: str | None, model: Model, drawings_path: str, file_name: str | None, **drawing_options: dict[str, Any]) -> Path | None:

        if file_name is None:
            file_name = f'{self.sg_id}_{self.sg_cut_id}_{self.fs_cut_id}_{str(self.amplitude_side).lower()}_{self.graph.name}'
            g_id = f'#(sg={self.fs_cut_id},sg_cut={self.sg_cut_id},fs_cut={self.fs_cut_id},side={str(self.amplitude_side).lower()})'
        else:
            g_id = f'#{self.fs_cut_id}'
        return self.graph.draw(drawing_mode, model, pjoin(drawings_path, file_name),
                               caption=f"Amplitude {self.graph.name.replace(f'_{self.fs_cut_id}', '').replace('_', ' ')}", diagram_id=g_id, **drawing_options)

    @staticmethod
    def from_serializable_dict(model: Model, amplitude_graph_dict: dict[str, Any]) -> AmplitudeGraph:

        return AmplitudeGraph(
            amplitude_graph_dict['sg_id'],
            amplitude_graph_dict['sg_cut_id'],
            amplitude_graph_dict['fs_cut_id'],
            Side[amplitude_graph_dict['amplitude_side']],
            Graph.from_serializable_dict(model, amplitude_graph_dict['graph']),
            amplitude_graph_dict['multiplicity'],
            amplitude_graph_dict['multi_channeling_channels']
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'sg_id': self.sg_id,
            'sg_cut_id': self.sg_cut_id,
            'fs_cut_id': self.fs_cut_id,
            'amplitude_side': str(self.amplitude_side),
            'graph': self.graph.to_serializable_dict(),
            'multiplicity': self.multiplicity,
            'multi_channeling_channels': self.multi_channeling_channels,
        }
