from __future__ import annotations

import os
import itertools
from pathlib import Path
from enum import StrEnum
import re
from typing import Any
import yaml
import pydot

from gammaloop.misc.common import GammaLoopError, DATA_PATH, pjoin, logger, EMPTY_LIST  # pylint: disable=unused-import # type: ignore
import gammaloop.misc.utils as utils
# type: ignore # pylint: disable=unused-import
from gammaloop.base_objects.model import Model, Propagator, VertexRule, Particle, Parameter


class EdgeType(StrEnum):
    INCOMING = 'in'
    OUTGOING = 'out'
    VIRTUAL = 'virtual'

# Abstract class for all vertex info


class VertexInfo(object):

    @staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict[str, Any]) -> VertexInfo:
        if serializable_dict['type'] == 'interacton_vertex_info':
            return InteractonVertexInfo.from_serializable_dict(model, serializable_dict)
        elif serializable_dict['type'] == 'external_vertex_info':
            return ExternalVertexInfo.from_serializable_dict(model, serializable_dict)
        elif serializable_dict['type'] == 'unspecified':
            return VertexInfo()
        else:
            raise GammaLoopError(
                f"Unknown vertex info type: {serializable_dict['type']}")

    def get_particles(self) -> list[Particle]:
        raise NotImplementedError()

    def get_vertex_rule_name(self) -> str | None:
        raise NotImplementedError()

    def is_external(self) -> bool:
        return False

    def get_type(self) -> str:
        return 'unspecified'

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'type': self.get_type()
        }


class InteractonVertexInfo(VertexInfo):
    def __init__(self, vertex_rule: VertexRule):
        self.vertex_rule: VertexRule = vertex_rule

    def get_type(self) -> str:
        return 'interacton_vertex_info'

    def get_vertex_rule_name(self) -> str | None:
        return self.vertex_rule.name

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'type': self.get_type(),
            'vertex_rule': self.vertex_rule.name
        }

    @staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict[str, Any]) -> InteractonVertexInfo:
        return InteractonVertexInfo(model.get_vertex_rule(serializable_dict['vertex_rule']))

    def get_particles(self) -> list[Particle]:
        return self.vertex_rule.particles


class ExternalVertexInfo(VertexInfo):
    def __init__(self, direction: EdgeType, particle: Particle):

        self.direction: EdgeType = direction
        self.particle: Particle = particle

    def is_external(self) -> bool:
        return True

    def get_type(self) -> str:
        return 'external_vertex_info'

    def get_vertex_rule_name(self) -> str | None:
        return None

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'type': self.get_type(),
            'particle': self.particle.name,
            'direction': self.direction.value,
        }

    @staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict[str, Any]) -> ExternalVertexInfo:
        return ExternalVertexInfo(
            EdgeType(serializable_dict['direction']),
            model.get_particle(serializable_dict['particle'])
        )

    def get_particles(self) -> list[Particle]:
        return [self.particle,]


class Vertex(object):
    def __init__(self, name: str, vertex_info: VertexInfo | None = None, edges: list[Edge] | None = None):
        self.name: str = name
        if vertex_info is None:
            self.vertex_info: VertexInfo = VertexInfo()
        else:
            self.vertex_info: VertexInfo = vertex_info
        if edges is None:
            self.edges: list[Edge] = []
        else:
            self.edges: list[Edge] = edges

    @staticmethod
    def default() -> Vertex:
        return Vertex('dummy_vertex')

    def set_vertex_rule_from_model(self, model: Model) -> None:
        assert self.edges is not None

        if len(self.edges) == 1:
            self.vertex_info = ExternalVertexInfo(
                EdgeType.OUTGOING if self.edges[0].edge_type == EdgeType.OUTGOING else EdgeType.INCOMING, self.edges[0].particle)
            return

        if self.vertex_info.get_type() == 'unspecified':
            vertex_particles = [e.particle if e.vertices[1] is self else e.particle.get_anti_particle(
                model) for e in self.edges]
            interactions = model.get_vertices_from_particles(vertex_particles)
            if len(interactions) == 0:
                raise GammaLoopError(
                    f"No interaction found for vertex {self.name} with particles [{','.join(p.name for p in vertex_particles)}]")
            elif len(interactions) > 1:
                raise GammaLoopError(
                    f"Found more than one valid interaction for vertex {self.name} with particles [{','.join(p.name for p in vertex_particles)}]. Importing from qgraph is not supported in this case.")
            self.vertex_info = InteractonVertexInfo(interactions[0])
        elif not isinstance(self.vertex_info, InteractonVertexInfo):
            raise GammaLoopError(
                "Vertex rule already set to a different type.")

    def sort_edges_according_to_vertex_info(self, model: Model) -> None:
        vertex_particles: list[Particle] = self.vertex_info.get_particles()
        if len(vertex_particles) <= 1:
            return
        sorted_edges: list[Edge] = []
        vertex_particles_per_edge: list[tuple[Particle, Edge]] = []
        self_loops = []
        for e in self.edges:
            if e.vertices[0] == e.vertices[1]:
                # self-loop treatment
                if e in self_loops:
                    continue
                self_loops.append(e)
                vertex_particles_per_edge.append((e.particle, e))
                vertex_particles_per_edge.append(
                    (e.particle.get_anti_particle(model), e))
            else:
                if e.vertices[1] is self:
                    vertex_particles_per_edge.append((e.particle, e))
                elif e.vertices[0] is self:
                    vertex_particles_per_edge.append(
                        (e.particle.get_anti_particle(model), e))
                else:
                    raise GammaLoopError(
                        f"Vertex {self.name} is not connected to edge {e.name}.")
        for part in vertex_particles:
            for i, (e_p, _edge) in enumerate(vertex_particles_per_edge):
                if e_p == part:
                    sorted_edges.append(vertex_particles_per_edge.pop(i)[1])
                    break
            else:
                print(sorted_edges, len(sorted_edges))
                raise GammaLoopError(
                    f"Particle {part.name} not found in vertex {self.name} with particles [{','.join(p.name for p in vertex_particles)}].")
        if len(vertex_particles_per_edge) > 0:
            raise GammaLoopError(
                f"Not all particles were found in vertex {self.name} with particles [{','.join(p.name for p in vertex_particles)}].")
        self.edges = sorted_edges

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'name': self.name,
            'vertex_info': self.vertex_info.to_serializable_dict(),
            'edges': [e.name for e in self.edges]
        }

    def to_pydot(self) -> pydot.Node:
        node_attr = {}
        if self.vertex_info.get_vertex_rule_name() is not None:
            node_attr['int_id'] = self.vertex_info.get_vertex_rule_name()
        return pydot.Node(
            self.name,
            **node_attr
        )

    def draw(self, graph: Graph, _model: Model, constant_definitions: dict[str, str], show_vertex_labels: bool = False,
             use_vertex_names: bool = False, vertex_size: float = 2.5, vertex_shape: str = 'circle', **_opts: dict[str, Any]) -> list[str]:
        # Possible shape choices are 'circle', 'square', 'triangle', 'diamond', 'pentagon', 'hexagon', 'triagram',
        # 'tetragram', 'pentagram', 'hexagram', 'triacross', 'cross', 'pentacross' and 'hexacross'.

        constant_definitions['vertexShape'] = vertex_shape
        constant_definitions['vertexSize'] = f'{vertex_size:.2f}thick'
        constant_definitions['blobSize'] = f'{vertex_size*3:.2f}'

        drawing: list[str] = []
        if self.vertex_info.is_external():
            return drawing

        vertex_pos: int = graph.get_vertex_position(self.name)
        prefix = r'\showVertexLabel{'
        if not show_vertex_labels:
            prefix = r'\hideVertexLabel{'
        drawing.append(r'%s\fmflabel{%s}{v%d}}' % (  # pylint: disable=consider-using-f-string
            prefix, self.name if use_vertex_names else r'$v_{%s}$' % vertex_pos, vertex_pos))  # pylint: disable=consider-using-f-string

        # One could tag the form factors vertices with blobs here
        is_blob = False
        if not is_blob:
            drawing.append(
                r'\efmfv{decor.shape=\vertexShape,decor.filled=full,label.dist=\labelDistance,decor.size=\vertexSize}{v%d}' % (  # pylint: disable=consider-using-f-string
                    vertex_pos))
        else:
            drawing.append(
                r'\efmfv{decor.shape=\vertexShape,decor.filled=shaded,label.dist=\labelDistance,decor.size=\blobSize}{v%d}' % (  # pylint: disable=consider-using-f-string
                    vertex_pos))

        return drawing


class Edge(object):
    def __init__(self, name: str, edge_type: EdgeType, particle: Particle, propagator: Propagator | None = None, vertices: tuple[Vertex, Vertex] | None = None):
        self.name: str = name
        self.edge_type: EdgeType = edge_type
        self.particle: Particle = particle
        if propagator is None:
            self.propagator: Propagator = Propagator.from_particle(
                particle, 'Feynman')
        else:
            self.propagator: Propagator = propagator
        if vertices is None:
            self.vertices: tuple[Vertex, Vertex] = (
                Vertex.default(), Vertex.default())
        else:
            self.vertices: tuple[Vertex, Vertex] = vertices

    @staticmethod
    def default() -> Edge:
        return Edge('dummy_edge', EdgeType.VIRTUAL, Particle.default())

    def get_vertex_positions(self, graph: Graph) -> tuple[int, int]:
        return (
            graph.get_vertex_position(self.vertices[0].name),
            graph.get_vertex_position(self.vertices[1].name),
        )

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'name': self.name,
            'edge_type': self.edge_type.value,
            'particle': self.particle.name,
            'propagator': self.propagator.name,
            'vertices': [self.vertices[0].name, self.vertices[1].name]
        }

    def to_pydot(self, **attributes) -> pydot.Edge:
        return pydot.Edge(
            self.vertices[0].name, self.vertices[1].name, pdg=self.particle.pdg_code, **attributes)

    def draw(self, graph: Graph, _model: Model, constant_definitions: dict[str, str], show_edge_labels: bool = True, show_particle_names: bool = True, show_edge_names: bool = True, show_edge_composite_momenta: bool = False,
             label_size: str = '15pt', label_distance: float = 13.0, show_edge_momenta: bool = True, draw_arrow_for_all_edges: bool = True, draw_lmb: bool = True,
             arc_max_distance: float = 1., external_legs_tension: float = 3., default_tension: float = 1.0, line_width: float = 1.0,
             arrow_size_for_single_line: float = 2.5, arrow_size_for_double_line: float = 2.5, line_color: str = 'black', label_color: str = 'black', lmb_color: str = 'red', non_lmb_color: str = 'blue', **_opts: Any) -> list[str]:
        constant_definitions['arcMaxDistance'] = f'{arc_max_distance:.2f}'
        constant_definitions['externalLegTension'] = f'{external_legs_tension:.2f}'  # nopep8
        constant_definitions['defaultTension'] = f'{default_tension:.2f}'
        constant_definitions['lineWidth'] = f'{line_width:.2f}'
        constant_definitions['lineColor'] = line_color
        constant_definitions['lmbColor'] = lmb_color
        constant_definitions['labelColor'] = label_color
        constant_definitions['labelSize'] = label_size
        constant_definitions['nonLmbColor'] = non_lmb_color
        constant_definitions['arrowSize'] = f'{arrow_size_for_single_line:.2f}'
        constant_definitions['doubleArrowSize'] = f'{arrow_size_for_double_line:.2f}'  # nopep8
        constant_definitions['labelDistance'] = f'{label_distance:.2f}'

        template_line = r'\efmf{%(line_type)s%(comma)s%(options)s}{%(left_vertex)s,%(right_vertex)s}'

        if graph.loop_momentum_basis is None:
            raise GammaLoopError(
                "Specify a loop momentum basis drawing a graph.")

        lmb_position = graph.loop_momentum_basis.index(
            self) if self in graph.loop_momentum_basis else None
        line_options: dict[str, str] = {}

        if self.particle.is_ghost():
            line_type = 'dots_arrow'
        else:
            match self.particle.spin:
                case 1:
                    if draw_arrow_for_all_edges:
                        line_type = 'dashes_arrow'
                    else:
                        line_type = 'dashes'
                case 2:
                    line_type = 'plain_arrow'
                case 3:
                    if self.particle.color == 1:
                        line_type = 'wiggly'
                    else:
                        line_type = 'curly'
                case 4:
                    line_type = 'zigzag'
                case 5:
                    line_type = 'dbl_wiggly'
                case _:
                    line_type = 'dbl_wiggly'

        if self.particle.is_massive() and abs(self.particle.spin) != 3:
            line_type = f'dbl_{line_type}'

        line_options['width'] = r'\lineWidth'
        line_options['fore'] = r'\lineColor'
        line_options['tension'] = r'\defaultTension'

        edge_name_label = ''
        if show_edge_names:
            edge_name_label = self.name  # pylint: disable=consider-using-f-string
        # Sanitize edge_name_label for latex
        edge_name_label = edge_name_label.replace('"', '')
        edge_name_label = edge_name_label.replace('_', '')

        label_components: list[str] = []
        if show_particle_names:
            # \overline{x} sadly does not work with the \fmfX{} command.
            # We must therefore change by appending an x
            label_components.append(
                re.sub(r"\\overline\{([^}]*)\}", r"\1x", self.particle.texname))

        if show_edge_momenta:
            if self.edge_type != EdgeType.VIRTUAL:
                for i_p, (e_1, e_2) in enumerate(graph.external_connections):
                    if self.edge_type == EdgeType.INCOMING and e_1 == self.vertices[0]:
                        label_components.append(r'p_{%d}' % (  # pylint: disable=consider-using-f-string
                            i_p+1))
                        break
                    elif self.edge_type == EdgeType.OUTGOING and e_2 == self.vertices[1]:
                        label_components.append(r'p_{%d}' % (  # pylint: disable=consider-using-f-string
                            i_p+1))
                        break
            else:
                if lmb_position is not None:
                    label_components.append(r'k_{%d}' % (  # pylint: disable=consider-using-f-string
                        lmb_position))
        if show_edge_composite_momenta:
            mom_sig = graph.get_edge_signature(self.name)
            if (mom_sig[0].count(0) + mom_sig[1].count(0)) < len(mom_sig[0])+len(mom_sig[1])-1:
                str_mom = ''.join(
                    (('+' if s > 0 else '-' if abs(s)
                     == 1 else f'{s:+d}')+f'k_{{{i_lm}}}')
                    for i_lm, s in enumerate(mom_sig[0]) if s != 0)+''.join(
                    (('+' if s > 0 else '-' if abs(s)
                     == 1 else f'{s:+d}')+f'p_{{{i_ext+1}}}')
                    for i_ext, s in enumerate(mom_sig[1]) if s != 0)
                if str_mom.startswith('+'):
                    str_mom = str_mom[1:]
                if any(char in str_mom for char in ['+', '-']):
                    str_mom = f'({str_mom})'
                label_components.append(str_mom)

        line_options['label'] = '%s%s$%s%s$' % (  # pylint: disable=consider-using-f-string
            r"\setLabelColor{\labelColor}" if lmb_position is None else r"\setLabelColor{\lmbColor}",
            edge_name_label,
            '|' if edge_name_label != '' and len(label_components) > 0 else '',
            '|'.join(label_components)
        )
        line_options['label.dist'] = r'\labelDistance'
        if show_edge_labels:
            line_options['label'] = r'\showEdgeLabel{%s}' % line_options[  # pylint: disable=consider-using-f-string
                'label']
        else:
            line_options['label'] = r'\hideEdgeLabel{%s}' % line_options[  # pylint: disable=consider-using-f-string
                'label']

        multi_edges = [e for e in graph.edges if set(
            e.vertices) == set(self.vertices)]
        self_position = multi_edges.index(self)
        base_side = 'left' if self.vertices[0] == multi_edges[0].vertices[0] else 'right'
        other_side = 'right' if self.vertices[0] == multi_edges[0].vertices[0] else 'left'
        if len(multi_edges) == 2:
            match self_position:
                case 0:
                    line_options[base_side] = str(0.3)
                case 1:
                    line_options[other_side] = str(0.3)
                case _:
                    raise GammaLoopError("Unreachable code.")
        elif len(multi_edges) == 3:
            match self_position:
                case 0:
                    line_options[base_side] = str(0.4)
                case 1:
                    pass
                case 2:
                    line_options[other_side] = str(0.4)
                case _:
                    raise GammaLoopError("Unreachable code.")
        elif len(multi_edges) > 3:
            if len(multi_edges) % 2 == 1:
                max_range = float(len(multi_edges)//2)
                self_relative_position = (
                    self_position-len(multi_edges)//2) / float(len(multi_edges)//2)
            else:
                max_range = (len(multi_edges)-1)/2.0
            self_relative_position = (self_position-max_range)/max_range
            if self_relative_position < 0.:
                line_options[base_side] = str(
                    abs(self_relative_position))
            elif self_relative_position == 0.:
                pass
            elif self_relative_position > 0.:
                line_options[other_side] = str(
                    abs(self_relative_position))
        if 'left' in line_options:
            line_options['left'] = line_options["left"]+r'*\arcMaxDistance'
        if 'right' in line_options:
            line_options['right'] = line_options["right"]+r'*\arcMaxDistance'

        replace_dict = {}
        replace_dict['line_type'] = line_type.replace('_arrow', '')
        replace_dict['left_vertex'] = f'v{graph.get_vertex_position(self.vertices[0].name)}'  # nopep8
        replace_dict['right_vertex'] = f'v{graph.get_vertex_position(self.vertices[1].name)}'  # nopep8
        replace_dict['options'] = ','.join(
            f'{k}={v}' for k, v in line_options.items())
        replace_dict['comma'] = ',' if len(replace_dict['options']) > 0 else ''
        drawing: list[str] = []
        drawing.append(
            f"% == Drawing edge '{self.name}' (e{graph.get_edge_position(self.name)}) with particle {self.particle.name}")
        drawing.append(template_line % replace_dict)

        if ('arrow' not in line_type and draw_arrow_for_all_edges) or 'arrow' in line_type:
            if 'label' in line_options:
                del line_options['label']
            if any(t in line_type for t in ['dbl', 'zigzag', 'wiggly', 'curly']):
                replace_dict['line_type'] = 'my_double_phantom_arrow'
            else:
                replace_dict['line_type'] = 'my_phantom_arrow'
            line_options['tension'] = '0'
            line_options['fore'] = r'\nonLmbColor'
            replace_dict['options'] = ','.join(
                f'{k}={v}' for k, v in line_options.items())
            replace_dict['comma'] = ',' if len(
                replace_dict['options']) > 0 else ''
            drawing.append(
                "% Draw an arrow to indicate orientations of vectors as well")
            drawing.append(template_line % replace_dict)

        if draw_lmb and lmb_position is not None:
            if 'label' in line_options:
                del line_options['label']
            if any(t in line_type for t in ['dbl', 'zigzag', 'wiggly', 'curly']):
                replace_dict['line_type'] = 'my_double_phantom_arrow'
            else:
                replace_dict['line_type'] = 'my_phantom_arrow'
            line_options['tension'] = '0'
            line_options['fore'] = r'\lmbColor'
            replace_dict['options'] = ','.join(
                f'{k}={v}' for k, v in line_options.items())
            replace_dict['comma'] = ',' if len(
                replace_dict['options']) > 0 else ''
            drawing.append("% Draw a colored arrow for this LMB edge")
            drawing.append(template_line % replace_dict)

        if self.edge_type != EdgeType.VIRTUAL or all(s == 0 for s in graph.get_edge_signature(self.name)[0]):
            # Add extra tension to externals
            if 'label' in line_options:
                del line_options['label']
            replace_dict['line_type'] = 'phantom'
            line_options['tension'] = r'\externalLegTension'
            line_options['fore'] = 'black'
            replace_dict['options'] = ','.join(
                f'{k}={v}' for k, v in line_options.items())
            replace_dict['comma'] = ',' if len(
                replace_dict['options']) > 0 else ''
            drawing.append(
                "% Draw a phantom additional tension for external leg")
            drawing.append(template_line % replace_dict)

        return drawing


class Graph(object):
    def __init__(self, name: str, vertices: list[Vertex], edges: list[Edge], external_connections: list[tuple[Vertex | None, Vertex | None]],
                 loop_momentum_basis: list[Edge] | None = None, overall_factor: str | int = "1",
                 edge_signatures: dict[str, tuple[list[int], list[int]]] | None = None):
        self.name: str = name
        self.vertices: list[Vertex] = vertices
        self.edges: list[Edge] = edges
        if edge_signatures is None:
            self.edge_signatures: dict[str,
                                       tuple[list[int], list[int]]] = {}
        else:
            self.edge_signatures: dict[str,
                                       tuple[list[int], list[int]]] = edge_signatures

        if isinstance(overall_factor, str):
            overall_factor_str = overall_factor.replace('"', '')
        else:
            overall_factor_str = str(overall_factor)
        self.overall_factor: str = overall_factor_str

        # For forward scattering graphs, keep track of the bipartite map, i.e. which in and out externals will carry identical momenta.
        self.external_connections: list[tuple[
            Vertex | None, Vertex | None]] = external_connections
        # if loop_momentum_basis is None:
        #     self.loop_momentum_basis: list[Edge] = []
        # else:
        self.loop_momentum_basis: list[Edge] | None = loop_momentum_basis
        self.name_to_position: dict[str, dict[str, int]] = {}

    def are_edges_and_vertices_list_consistent(self) -> None:
        """ Tests that the list of edges assigned to each vertex is consistent with the list of vertices assigned to each edge."""
        edge_ids_per_vertex_id: dict[str, set[str]] = {
            v.name: set() for v in self.vertices}
        for edge in self.edges:
            for vertex in edge.vertices:
                edge_ids_per_vertex_id[vertex.name].add(edge.name)
        for vertex in self.vertices:
            if set(e.name for e in vertex.edges) != edge_ids_per_vertex_id[vertex.name]:
                raise GammaLoopError(
                    f"Graph '{self.name}' has inconsistent list of edges and vertices. Vertex '{vertex.name}' is defined with edges {set(e.name for e in vertex.edges)} but has the following edges pointing to it {edge_ids_per_vertex_id[vertex.name]}.")

    def get_sorted_incoming_edges(self) -> list[Edge]:
        sorted_incoming_edges: list[Edge] = []
        for (u, _v) in self.external_connections:  # pylint: disable=invalid-name
            if u is None:
                continue
            for e in self.get_incoming_edges():  # pylint: disable=invalid-name
                if e.vertices[0] == u:
                    sorted_incoming_edges.append(e)
                    break
        return sorted_incoming_edges

    def get_sorted_outgoing_edges(self) -> list[Edge]:
        sorted_outgoing_edges: list[Edge] = []
        for (_u, v) in self.external_connections:  # pylint: disable=invalid-name
            if v is None:
                continue
            for e in self.get_outgoing_edges():  # pylint: disable=invalid-name
                if e.vertices[1] == v:
                    sorted_outgoing_edges.append(e)
                    break
        return sorted_outgoing_edges

    def get_edge_signature(self, edge_name: str) -> tuple[list[int], list[int]]:
        return self.edge_signatures[edge_name]

    @ staticmethod
    def empty_graph(name: str) -> Graph:
        return Graph(name, [], [], [], loop_momentum_basis=[], overall_factor="1", edge_signatures={})

    def is_empty(self) -> bool:
        return len(self.vertices) == 0 and len(self.edges) == 0

    def get_edge(self, edge_name: str) -> Edge:
        return self.edges[self.name_to_position['edges'][edge_name]]

    def get_vertex(self, vertex_name: str) -> Vertex:
        return self.vertices[self.name_to_position['vertices'][vertex_name]]

    def get_edge_position(self, edge_name: str) -> int:
        return self.name_to_position['edges'][edge_name]

    def get_vertex_position(self, vertex_name: str) -> int:
        return self.name_to_position['vertices'][vertex_name]

    def get_position(self, item: Edge | Vertex) -> int:
        if isinstance(item, Edge):
            return self.get_edge_position(item.name)
        else:
            return self.get_vertex_position(item.name)

    def get_incoming_edges(self) -> list[Edge]:
        return [e for e in self.edges if e.edge_type == EdgeType.INCOMING]

    def get_outgoing_edges(self) -> list[Edge]:
        return [e for e in self.edges if e.edge_type == EdgeType.OUTGOING]

    def get_external_edges(self) -> list[Edge]:
        return self.get_incoming_edges() + self.get_outgoing_edges()

    def synchronize_name_map(self) -> None:
        self.name_to_position['edges'] = {}
        for i_e, edge in enumerate(self.edges):
            self.name_to_position['edges'][edge.name] = i_e
        self.name_to_position['vertices'] = {}
        for i_v, vertex in enumerate(self.vertices):
            self.name_to_position['vertices'][vertex.name] = i_v

    def to_serializable_dict(self) -> dict[str, Any]:
        return {
            'name': self.name,
            'vertices': [v.to_serializable_dict() for v in self.vertices],
            'edges': [e.to_serializable_dict() for e in self.edges],
            'external_connections': [[e[0].name if e[0] is not None else None,
                                      e[1].name if e[1] is not None else None] for e in self.external_connections],
            'overall_factor': self.overall_factor,
            'loop_momentum_basis': [e.name for e in self.loop_momentum_basis] if self.loop_momentum_basis is not None else None,
            'edge_signatures': sorted([(k, v) for k, v in self.edge_signatures.items()], key=lambda x: x[0])
        }

    @ staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict[str, Any]) -> Graph:

        graph_vertices: list[Vertex] = []
        vertex_name_to_vertex_map: dict[str, Vertex] = {}
        for v_graph in serializable_dict['vertices']:
            v_name = v_graph['name']
            vert = Vertex(v_name)
            vertex_name_to_vertex_map[v_name] = vert
            graph_vertices.append(vert)

        graph_edges: list[Edge] = []
        edge_name_to_edge_map = {}
        for e_qgraph in serializable_dict['edges']:
            e_name = e_qgraph['name']
            try:
                props = model.get_propagator(e_qgraph['propagator'])
            except KeyError:
                props = None
            edge = Edge(e_name, EdgeType(e_qgraph['edge_type']), model.get_particle(e_qgraph['particle']), props,
                        (vertex_name_to_vertex_map[e_qgraph['vertices'][0]],
                        vertex_name_to_vertex_map[e_qgraph['vertices'][1]])
                        )
            edge_name_to_edge_map[e_name] = edge
            graph_edges.append(edge)

        for (vert, v_graph) in zip(graph_vertices, serializable_dict['vertices']):
            vert.edges = [edge_name_to_edge_map[e_name]
                          for e_name in v_graph['edges']]
            vert.vertex_info = VertexInfo.from_serializable_dict(
                model, v_graph['vertex_info'])
            vert.sort_edges_according_to_vertex_info(model)

        external_connections = [(vertex_name_to_vertex_map[e[0]] if e[0] is not None else None,
                                 vertex_name_to_vertex_map[e[1]] if e[1] is not None else None) for e in serializable_dict['external_connections']]

        overall_factor = serializable_dict['overall_factor']
        loop_momentum_basis: list[Edge] = [edge_name_to_edge_map[e_name]
                                           for e_name in serializable_dict['loop_momentum_basis']]
        graph = Graph(serializable_dict['name'], graph_vertices, graph_edges,
                      external_connections, loop_momentum_basis, overall_factor,
                      dict(serializable_dict['edge_signatures']))
        graph.synchronize_name_map()

        return graph

    @ staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> Graph:
        graph_dict = yaml.safe_load(yaml_str)
        return Graph.from_serializable_dict(model, graph_dict)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump(self.to_serializable_dict())

    def to_pydot(self, graph_attributes: dict[str, Any]) -> pydot.Graph:
        graph = pydot.Graph(graph_type='digraph',
                            graph_name=self.name,
                            overall_factor=self.overall_factor,
                            **graph_attributes)
        for vertex in self.vertices:
            graph.add_node(vertex.to_pydot())

        sorted_incoming_edges: list[tuple[Edge, dict[str, Any]]] = []
        sorted_outgoing_edges: list[tuple[Edge, dict[str, Any]]] = []
        for i_ext, (u, v) in enumerate(self.external_connections):
            if u is not None:
                for e in self.get_incoming_edges():
                    if e.vertices[0] == u:
                        sorted_incoming_edges.append(
                            (e, {"name": f"p{i_ext+1}_IN" if v is not None else f"p{i_ext+1}", "mom": f"p{i_ext+1}"}))
                        break
            if v is not None:
                for e in self.get_outgoing_edges():
                    if e.vertices[1] == v:
                        sorted_outgoing_edges.append(
                            (e, {"name": f"p{i_ext+1}_OUT" if u is not None else f"p{i_ext+1}", "mom": f"p{i_ext+1}"}))
                        break
        sorted_virtual_edges: list[tuple[Edge, dict[str, Any]]] = []
        for i_e_virt, e in enumerate([e for e in self.edges if e.edge_type == EdgeType.VIRTUAL]):
            edge_attrs: dict[str, Any] = {"name": f"q{i_e_virt+1}"}
            if self.loop_momentum_basis is not None and e in self.loop_momentum_basis:
                edge_attrs["lmb_index"] = self.loop_momentum_basis.index(e)
            sorted_virtual_edges.append((e, edge_attrs))
        sorted_virtual_edges = sorted(
            sorted_virtual_edges, key=lambda e: e[0].name)
        for edge, edge_attr in sorted_incoming_edges + sorted_outgoing_edges + sorted_virtual_edges:
            graph.add_edge(edge.to_pydot(**edge_attr))

        return graph

    @staticmethod
    def from_pydot(model: Model, g: pydot.Graph) -> Graph:
        """ Imports graph from a pydot graph object. If some information is missing in the graph supplied, it will be modified in-place to make it consistent."""

        # Sanity check
        for e in g.get_edges():
            if 'pdg' not in e.get_attributes():
                raise GammaLoopError(
                    f"Dot Edge {e} does not have a PDG assigned.")

        pydot_edges: list[pydot.Edge] = g.get_edges()

        # Add nodes that may not be specified
        existing_node_names = {n.get_name() for n in g.get_nodes()}
        for e in pydot_edges:
            if e.get_source() not in existing_node_names:
                existing_node_names.add(e.get_source())
                g.add_node(pydot.Node(e.get_source()))
            if e.get_destination() not in existing_node_names:
                existing_node_names.add(e.get_destination())
                g.add_node(pydot.Node(e.get_destination()))
        pydot_nodes: list[pydot.Node] = g.get_nodes()

        adjacency_map: dict[str, list[tuple[pydot.Edge, EdgeType]]] = {
            n.get_name(): [
                (e, EdgeType.INCOMING if e.get_destination()
                 == n.get_name() else EdgeType.OUTGOING)
                for e in pydot_edges if n.get_name() in [e.get_source(), e.get_destination()]
            ]
            for n in pydot_nodes
        }

        # The logic below is just to assign names if they are absent
        incoming_momenta = []
        outgoing_momenta = []
        # First pass to collect momenta names
        for n in pydot_nodes:
            if len(adjacency_map[n.get_name()]) == 1:
                half_edge, direction = adjacency_map[n.get_name()][0]
                if "mom" not in half_edge.get_attributes():
                    raise GammaLoopError(
                        f"External dot half-edge {half_edge} does not have a momentum assigned.")
                if direction == EdgeType.INCOMING:
                    incoming_momenta.append(half_edge.get_attributes()["mom"])
                else:
                    outgoing_momenta.append(half_edge.get_attributes()["mom"])
        # Second pass to assign edge names if missing
        i_internal_edge = 0
        for e in pydot_edges:
            e_attr = e.get_attributes()
            if len(adjacency_map[e.get_source()]) == 1:
                if "name" not in e_attr:
                    if e_attr["mom"] in incoming_momenta:
                        e_attr["name"] = f"{e_attr['mom']}_IN"
                    else:
                        e_attr["name"] = e_attr['mom']
            elif len(adjacency_map[e.get_destination()]) == 1:
                if "name" not in e_attr:
                    if e_attr["mom"] in outgoing_momenta:
                        e_attr["name"] = f"{e_attr['mom']}_OUT"
                    else:
                        e_attr["name"] = e_attr['mom']
            else:
                i_internal_edge += 1
                if "name" not in e_attr:
                    e_attr["name"] = f"q{i_internal_edge}"

        # Sort edge by name, with first incoming externals, then outgoing externals, then virtual edges.
        sorted_pydot_incoming_edges = sorted(
            [e for e in pydot_edges if len(adjacency_map[e.get_source()]) == 1], key=lambda e: e.get_attributes()["name"])
        sorted_pydot_outgoing_edges = sorted(
            [e for e in pydot_edges if len(adjacency_map[e.get_destination()]) == 1], key=lambda e: e.get_attributes()["name"])
        for e in sorted_pydot_incoming_edges+sorted_pydot_outgoing_edges:
            if "mom" not in e.get_attributes():
                raise GammaLoopError(
                    f"External dot half-edge {e} does not have a momentum assigned.")
        sorted_pydot_virtual_edges = sorted(
            [e for e in pydot_edges if len(adjacency_map[e.get_source()]) > 1 and len(
                adjacency_map[e.get_destination()]) > 1],
            key=lambda e: e.get_attributes()["name"])
        pydot_edges = sorted_pydot_incoming_edges + \
            sorted_pydot_outgoing_edges + sorted_pydot_virtual_edges

        graph_vertices: list[Vertex] = []
        pydot_vertex_name_to_vertex_map: dict[str, Vertex] = {}

        for n in pydot_nodes:
            vert = Vertex(n.get_name())
            pydot_vertex_name_to_vertex_map[n.get_name()] = vert
            graph_vertices.append(vert)

        graph_edges: list[Edge] = []
        pydot_edge_name_to_edge_map: dict[str, Edge] = {}
        for e in pydot_edges:
            e_attr = e.get_attributes()
            if 'propagator' in e_attr:
                props = model.get_propagator(e_attr['propagator'])
            else:
                props = None
            if len(adjacency_map[e.get_source()]) == 1:
                edge_type = EdgeType.INCOMING
            elif len(adjacency_map[e.get_destination()]) == 1:
                edge_type = EdgeType.OUTGOING
            else:
                edge_type = EdgeType.VIRTUAL
            edge = Edge(e_attr['name'], edge_type, model.get_particle_from_pdg(int(e_attr['pdg'])), props,
                        (pydot_vertex_name_to_vertex_map[e.get_source()],
                        pydot_vertex_name_to_vertex_map[e.get_destination()])
                        )
            pydot_edge_name_to_edge_map[e_attr['name']] = edge
            graph_edges.append(edge)

        for (vert, node) in zip(graph_vertices, pydot_nodes):
            vert.edges = [pydot_edge_name_to_edge_map[e.get_attributes()["name"]] for (
                e, e_type) in adjacency_map[node.get_name()]]
            node_attr = node.get_attributes()
            if "int_id" in node_attr:
                # Sadly the pydot parser adds quotes around string attributes for nodes (but not for edges)
                processed_int_id = node_attr["int_id"].strip('"')
                vert.vertex_info = InteractonVertexInfo(
                    model.get_vertex_rule(processed_int_id))
            else:
                vert.set_vertex_rule_from_model(model)
            vert.sort_edges_according_to_vertex_info(model)

        external_connections: dict[str, list[Vertex | None]] = {}
        for n in sorted_pydot_incoming_edges:
            momentum = n.get_attributes()["mom"]
            if momentum in external_connections:
                external_connections[momentum][0] = pydot_vertex_name_to_vertex_map[n.get_source(
                )]
            else:
                external_connections[momentum] = [
                    pydot_vertex_name_to_vertex_map[n.get_source()], None]
        for n in sorted_pydot_outgoing_edges:
            momentum = n.get_attributes()["mom"]
            if momentum in external_connections:
                external_connections[momentum][1] = pydot_vertex_name_to_vertex_map[n.get_destination(
                )]
            else:
                external_connections[momentum] = [
                    None, pydot_vertex_name_to_vertex_map[n.get_destination()]]

        external_connections_for_graph: list[tuple[Vertex | None, Vertex | None]] = [(v[0], v[1]) for _, v in sorted(
            list(external_connections.items()), key=lambda el: el[0])]

        overall_factor = g.get_attributes().get("overall_factor", "1")
        graph = Graph(g.get_name(), graph_vertices, graph_edges,
                      external_connections_for_graph, None, overall_factor, None)
        graph.synchronize_name_map()

        # Enforce specified LMB if available
        # NO CHECKS PERFORMED!
        specified_lmb = [(int(e.get_attributes()["lmb_index"]), pydot_edge_name_to_edge_map[e.get_attributes()[
                          "name"]]) for e in sorted_pydot_virtual_edges if 'lmb_index' in e.get_attributes()]
        if len(specified_lmb) > 0:
            graph.loop_momentum_basis = [e for _, e in sorted(
                specified_lmb, key=lambda el: el[0])]
        else:
            all_lmbs = graph.generate_loop_momentum_bases(
                target_n_lmbs_to_generate=1)
            if len(all_lmbs) == 0:
                raise GammaLoopError(
                    f"Could not find any valid loop momentum basis for graph '{g.get_name()}'.")
            else:
                # Arbitrarily pick the first valid momentum basis as the chosen one.
                graph.loop_momentum_basis = all_lmbs[0]

        graph.edge_signatures = {graph.edges[i].name: sig for i, sig in enumerate(
            graph.generate_momentum_flow())}

        # Make sure vertices are defined with a consistent set of edge names
        graph.are_edges_and_vertices_list_consistent()

        return graph

    @staticmethod
    def from_qgraph(model: Model, qgraph_object: dict[str, Any], name: str = 'default') -> Graph:
        """ Imports graph from a stylicized qgraph python output file. Will be deprecated when using in-house gammaloop graph generation. """

        graph_vertices: list[Vertex] = []
        vertex_qgraph_index_to_vertex_map: dict[str, Vertex] = {}
        qgraph_vertices = list(qgraph_object['nodes'].items())
        for i_vertex, (v_qgraph_index, v_qgraph) in enumerate(qgraph_vertices):
            vert = Vertex(f"v{i_vertex+1}")
            vertex_qgraph_index_to_vertex_map[v_qgraph_index] = vert
            graph_vertices.append(vert)

        graph_edges: list[Edge] = []
        edge_qgraph_index_to_edge_map = {}
        qgraph_edges = list(qgraph_object['edges'].items())
        for e_qgraph_index, e_qgraph in qgraph_edges:
            try:
                props = model.get_propagator(e_qgraph['propagator'])
            except KeyError:
                props = None
            edge = Edge(e_qgraph['name'], EdgeType(e_qgraph['type']), model.get_particle_from_pdg(e_qgraph['PDG']), props,
                        (vertex_qgraph_index_to_vertex_map[e_qgraph['vertices'][0]],
                        vertex_qgraph_index_to_vertex_map[e_qgraph['vertices'][1]])
                        )
            edge_qgraph_index_to_edge_map[e_qgraph_index] = edge
            graph_edges.append(edge)

        for (vert, (v_qgraph_index, v_qgraph)) in zip(graph_vertices, qgraph_vertices):
            vert.edges = [edge_qgraph_index_to_edge_map[e_qgraph_index]
                          for e_qgraph_index in v_qgraph['edge_ids']]
            vert.set_vertex_rule_from_model(model)
            vert.sort_edges_according_to_vertex_info(model)

        external_connections: dict[str,
                                   list[Vertex | None]] = {}
        for node_id, node in qgraph_vertices:
            # Hack necessary because of inconsistent qgraph output style :'(
            if isinstance(node['momenta'], str):
                node['momenta'] = [node['momenta'],]
            if len(node['momenta']) != 1:
                continue
            if qgraph_object['edges'][node['edge_ids'][0]]['type'] == 'in':
                if node['momenta'][0] in external_connections:
                    if external_connections[node['momenta'][0]][0] is not None:
                        raise GammaLoopError(
                            f"Multiple incoming external nodes with identical mometum label {node['momenta'][0]} assigned ")
                    external_connections[node['momenta'][0]
                                         ][0] = vertex_qgraph_index_to_vertex_map[node_id]
                else:
                    external_connections[node['momenta'][0]] = [
                        vertex_qgraph_index_to_vertex_map[node_id], None]
            else:
                if node['momenta'][0] in external_connections:
                    if external_connections[node['momenta'][0]][1] is not None:
                        raise GammaLoopError(
                            f"Multiple incoming external nodes with identical mometum label {node['momenta'][0]} assigned ")
                    external_connections[node['momenta'][0]
                                         ][1] = vertex_qgraph_index_to_vertex_map[node_id]
                else:
                    external_connections[node['momenta'][0]] = [
                        None, vertex_qgraph_index_to_vertex_map[node_id]]

        external_connections_for_graph: list[tuple[Vertex | None, Vertex | None]] = [(v[0], v[1]) for _, v in sorted(
            list(external_connections.items()), key=lambda el: el[0])]

        graph = Graph(name, graph_vertices, graph_edges,
                      external_connections_for_graph, None, qgraph_object['overall_factor'], None)
        graph.synchronize_name_map()

        all_lmbs = graph.generate_loop_momentum_bases(
            target_n_lmbs_to_generate=1)
        if len(all_lmbs) == 0:
            raise GammaLoopError(
                f"Could not find any valid loop momentum basis for graph '{name}'.")
        else:
            # Arbitrarily pick the first valid momentum basis as the chosen one.
            graph.loop_momentum_basis = all_lmbs[0]

        graph.edge_signatures = {graph.edges[i].name: sig for i, sig in enumerate(
            graph.generate_momentum_flow())}

        # Make sure vertices are defined with a consistent set of edge names
        graph.are_edges_and_vertices_list_consistent()

        return graph

    def get_edge_map(self) -> list[tuple[int, int]]:
        result: list[tuple[int, int]] = []
        for edge in self.edges:
            result.append(
                (
                    self.name_to_position['vertices'][edge.vertices[0].name],
                    self.name_to_position['vertices'][edge.vertices[1].name]
                )
            )
        return result

    def get_adjacency_map(self, edge_map: list[tuple[int, int]]) -> dict[int, list[tuple[tuple[int, tuple[int, int]], int]]]:
        result: dict[int, list[tuple[tuple[int, tuple[int, int]], int]]] = {}
        for i_e, edge in enumerate(edge_map):
            result.setdefault(edge[0], []).append(((i_e, edge), edge[1]))
            result.setdefault(edge[1], []).append(((i_e, edge), edge[0]))
        return result

    def generate_loop_momentum_bases(self, target_n_lmbs_to_generate: int | None = None) -> list[list[Edge]]:
        """ Return all valid loop momentum bases for the graph"""

        # logger.debug("Generating loop momentum bases for graph %s", self.name)

        # This function need to be called for each disconnected subgraph of self.
        vertices_not_encountered: list[int] = list(range(len(self.vertices)))
        all_spanning_tree_lists_from_disconnected_components: list[list[tuple[int, ...]]] = [
        ]
        while len(vertices_not_encountered) > 0:
            seed_vertex_position = vertices_not_encountered[0]
            accumulated_edge_sequence: list[int] = []
            excluded_edges: list[int] = []
            seen_edge_sequences: set[tuple[int, ...]] | None = None
            spanning_trees_for_this_subgraph: list[tuple[int, ...]] = []
            edge_map = self.get_edge_map()
            adjacency_map = self.get_adjacency_map(edge_map)
            # Use the code below to profile this computationally intensive function
            # import cProfile as profile
            # profile.runctx(
            #     'utils.generate_spanning_trees(spanning_trees_for_this_subgraph, linear_edge_map, adjacency_map, {seed_vertex_position},accumulated_edge_sequence, excluded_edges, seen_edge_sequences)', globals(), locals())
            # import sys
            # sys.exit(0)
            utils.generate_spanning_trees(spanning_trees_for_this_subgraph, edge_map, adjacency_map, {seed_vertex_position},
                                          accumulated_edge_sequence, excluded_edges, seen_edge_sequences,
                                          target_number_of_st_to_find=target_n_lmbs_to_generate, follow=False)

            for spanning_tree in spanning_trees_for_this_subgraph:
                for e_id in spanning_tree:
                    for side in [0, 1]:
                        try:
                            vertices_not_encountered.remove(
                                self.get_vertex_position(self.edges[e_id].vertices[side].name))
                        except ValueError:
                            pass
            all_spanning_tree_lists_from_disconnected_components.append(
                spanning_trees_for_this_subgraph)

        # Now build the cartesian product of all lists of spanning tree components from disconnected subgraphs
        spanning_trees = [tuple(sum((list(st) for st in st_combination), EMPTY_LIST)) for st_combination in
                          itertools.product(*all_spanning_tree_lists_from_disconnected_components)]

        # logger.debug("At total of %d loop momentum bases were found for graph %s", len(
        #    spanning_trees), self.name)

        return [[e for e_i, e in enumerate(self.edges) if e_i not in spanning_tree] for spanning_tree in spanning_trees]

    def generate_momentum_flow(self) -> list[tuple[list[int], list[int]]]:
        if self.loop_momentum_basis is None:
            raise GammaLoopError(
                "Specify a loop momentum basis before generating a momentum flow.")
        if len(self.external_connections) == 0:
            # In that case simply route to the first edge as sink
            sink_edge = 0
        else:
            if self.external_connections[-1][1] is not None:
                sink_edge = self.get_edge_position(
                    self.get_sorted_outgoing_edges()[-1].name)
            else:
                sink_edge = self.get_edge_position(
                    self.get_sorted_incoming_edges()[-1].name)
        external_edges: list[int] = []
        incoming_edges = iter(self.get_sorted_incoming_edges())
        outgoing_edges = iter(self.get_sorted_outgoing_edges())
        for (in_v, out_v) in self.external_connections:
            if in_v is not None:
                external_edges.append(
                    self.get_edge_position(next(incoming_edges).name))
            if out_v is not None:
                external_edges.append(
                    self.get_edge_position(next(outgoing_edges).name))
        signatures = utils.generate_momentum_flow(self.get_edge_map(),
                                                  [self.get_edge_position(lmb_e.name) for lmb_e in self.loop_momentum_basis], sink_edge, external_edges)

        # Merge identical external momenta
        merged_signatures: list[tuple[list[int], list[int]]] = []
        for sig in signatures:
            external_sig = iter(sig[1])
            merged_external_sig = [0 for _ in range(
                len(self.external_connections))]
            for i_ext, (in_edge, out_edge) in enumerate(self.external_connections):
                if in_edge is not None:
                    merged_external_sig[i_ext] += next(external_sig)
                if out_edge is not None:
                    merged_external_sig[i_ext] += next(external_sig)
            merged_signatures.append((sig[0], merged_external_sig))

        return merged_signatures

    def draw(self, drawing_mode: str | None, model: Model, file_path_without_extension: str, caption: str | None = None, diagram_id: str | None = None, **drawing_options: Any) -> Path | None:
        match drawing_mode:
            case 'feynmp':
                return self.draw_feynmp(model, file_path_without_extension, caption, diagram_id, **drawing_options['feynmp'])
            case 'dot':
                return self.draw_dot(model, file_path_without_extension, caption, diagram_id, **drawing_options['dot'])
            case None:
                return None
            case _:
                raise GammaLoopError(
                    f"Feynman drawing mode '{drawing_mode}' is not supported. Currently only 'feynmp' is supported.")

    def draw_dot(self, model: Model, file_path_without_extension: str, caption: str | None = None, diagram_id: str | None = None, **drawing_options: Any) -> Path:

        replace_dict: dict[str, Any] = {}
        replace_dict['graph_name'] = os.path.basename(
            file_path_without_extension)
        if caption:
            replace_dict['label'] = caption
        else:
            replace_dict['label'] = replace_dict['graph_name']
        replace_dict['label'] += f" {diagram_id}" if diagram_id else ''
        replace_dict['layout'] = drawing_options['layout']
        overall_factor_evaluated_str = utils.expression_to_string(
            utils.evaluate_graph_overall_factor(self.overall_factor), canonical=True)
        if drawing_options.get('show_overall_factor', True):
            replace_dict['label'] += f" x ( {overall_factor_evaluated_str} )"
        replace_dict['graph_options'] = f'\noverall_factor="{self.overall_factor}",\noverall_factor_evaluated="{overall_factor_evaluated_str}"'  # nopep8
        if len(drawing_options['graph_options']) == 0:
            replace_dict['graph_options'] += '\n'
        else:
            replace_dict['graph_options'] += ',\n' + ','.join(
                [f"{o}={str(v)}" for o, v in drawing_options['graph_options'].items()]) + '\n'
        replace_dict['node_options'] = ','.join(
            f"{o}={str(v)}" for o, v in drawing_options['node_options'].items())
        edge_options = dict(drawing_options['edge_options'])
        edge_options['penwidth'] = edge_options.get('penwidth', 0.6)
        replace_dict['edge_options'] = ','.join(
            f"{o}={str(v)}" for o, v in edge_options.items())
        edges_str = []
        for edge in self.edges:
            edge_repl_dict = {
                "start": edge.vertices[0].name,
                "end": edge.vertices[1].name
            }
            edge_repl_dict["label"] = f"{edge.name} | {edge.particle.name}"
            # Paint gluons in red
            if edge.particle.pdg_code == 21:
                edge_repl_dict["color"] = "red"
            elif edge.particle.is_ghost():
                if edge.particle.color != 1:
                    edge_repl_dict["color"] = "magenta"
                else:
                    edge_repl_dict["color"] = "cyan"
            elif edge.particle.spin % 2 == 0:
                edge_repl_dict["color"] = "black"
            else:
                edge_repl_dict["color"] = "blue"
            if edge.particle.is_massive():
                edge_repl_dict["penwidth"] = edge_options['penwidth']*2
            else:
                edge_repl_dict["penwidth"] = edge_options['penwidth']
            edge_repl_dict["style"] = "solid"
            if self.loop_momentum_basis:
                if edge in self.loop_momentum_basis:
                    edge_repl_dict["style"] = "dashed"
            edges_str.append(
                '"{start:s}" -> "{end:s}" [label="{label:s}",color="{color:s}",penwidth="{penwidth}",style="{style:s}"];'.format(**edge_repl_dict))
        replace_dict['edges'] = '\n'.join(edges_str)

        with open(pjoin(DATA_PATH, 'templates', 'drawing', 'drawing.dot.template'), 'r', encoding='utf-8') as f:
            template = f.read()
        output_path = file_path_without_extension+'.dot'
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(template.format(**replace_dict))

        return Path(output_path)

    def draw_feynmp(self, model: Model, file_path_without_extension: str, caption: str | None = None, diagram_id: str | None = None, caption_size: str = '20pt', **drawing_options: Any) -> Path:

        with open(pjoin(DATA_PATH, 'templates', 'drawing', 'drawing.tex.template'), 'r', encoding='utf-8') as file:
            template = file.read()
        source_file_path = Path(file_path_without_extension+'.tex')

        constant_definitions: dict[str, str] = {}

        with open(source_file_path, 'w', encoding='utf-8') as file:
            replace_dict = {}
            replace_dict['drawing_name'] = os.path.basename(
                file_path_without_extension)
            if len(self.get_incoming_edges()) > 0 and len(self.get_outgoing_edges()) > 0:
                replace_dict['incoming_vertices'] = r"\fmfleft{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[0].name)}' for e in self.get_incoming_edges())
                outgoing_edges = self.get_outgoing_edges()
                if drawing_options.get('reverse_outgoing_edges_order', False):
                    outgoing_edges = outgoing_edges[::-1]
                replace_dict['outgoing_vertices'] = r"\fmfright{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[1].name)}' for e in outgoing_edges)
            elif len(self.get_incoming_edges()) > 0:
                replace_dict['incoming_vertices'] = r"\fmfsurround{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[0].name)}' for e in self.get_incoming_edges())
                replace_dict['outgoing_vertices'] = "% None."
            elif len(self.get_outgoing_edges()) > 0:
                replace_dict['outgoing_vertices'] = r"\fmfsurround{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[0].name)}' for e in self.get_incoming_edges())
                replace_dict['incoming_vertices'] = "% None."
            else:
                replace_dict['incoming_vertices'] = "% None."
                replace_dict['outgoing_vertices'] = "% None."
                # The anchor tension will do the job of spreading the graph open
                pass

            edge_map = self.get_edge_map()
            longest_cycle_edges = utils.find_longest_cycle(edge_map)
            if len(longest_cycle_edges) == 0:
                # Tree-level graphs
                longest_cycle_edges = [self.get_edge_position(
                    self.get_sorted_incoming_edges()[0].name),]
            longest_cycle = utils.edges_cycle_to_vertex_cycle([
                self.edges[i_e].get_vertex_positions(self) for i_e in longest_cycle_edges])
            incoming_edges_paths: list[tuple[Edge, list[int]]] = []
            for edge in self.get_sorted_incoming_edges():
                incoming_edges_paths.append((edge, utils.find_shortest_path(
                    edge_map, self.get_vertex_position(
                        edge.vertices[0].name),
                    longest_cycle)))

            outgoing_edges_paths: list[tuple[Edge, list[int]]] = []
            for edge in self.get_sorted_outgoing_edges():
                outgoing_edges_paths.append((edge, utils.find_shortest_path(
                    edge_map, self.get_vertex_position(
                        edge.vertices[1].name),
                    longest_cycle)))

            replace_dict['incoming_edge_definitions'] = '\n'.join(
                sum((edge.draw(self, model, constant_definitions, **drawing_options) for edge, _ in incoming_edges_paths), EMPTY_LIST))
            replace_dict['outgoing_edge_definitions'] = '\n'.join(
                sum((edge.draw(self, model, constant_definitions, **drawing_options) for edge, _ in outgoing_edges_paths), EMPTY_LIST))
            internal_edges_drawing: list[str] = []
            for edge in self.edges:
                if edge.edge_type != EdgeType.VIRTUAL:
                    continue
                internal_edges_drawing.extend(
                    edge.draw(self, model, constant_definitions, **drawing_options))
            replace_dict['internal_edge_definitions'] = '\n'.join(
                internal_edges_drawing)

            # Now add external tension
            replace_dict['external_tension_definition'] = '\n'.join(
                self.draw_anchor_tension(incoming_edges_paths, outgoing_edges_paths, constant_definitions, longest_cycle, **drawing_options))

            replace_dict['vertex_definitions'] = '\n'.join(
                sum((vertex.draw(self, model, constant_definitions, **drawing_options) for vertex in self.vertices), EMPTY_LIST))

            constant_definitions['captionSize'] = caption_size
            if caption:
                replace_dict['caption'] = caption
            else:
                replace_dict['caption'] = replace_dict['drawing_name']
            replace_dict['caption'] = replace_dict['caption'].replace(
                '_', ' ').replace('#', r'\#')
            if diagram_id is not None:
                replace_dict['diagram_id'] = diagram_id.replace('#', r'\#')
            else:
                replace_dict['diagram_id'] = ''
            replace_dict['constant_definitions'] = '\n'.join(
                r'\newcommand{\%s}{%s}' % (k, v) for k, v in constant_definitions.items()  # pylint: disable=consider-using-f-string
            )

            file.write(template % replace_dict)

        return source_file_path

    def draw_anchor_tension(self, incoming_edges_paths: list[tuple[Edge, list[int]]], outgoing_edges_paths: list[tuple[Edge, list[int]]], constant_definitions: dict[str, str],
                            longest_cycle: list[int] | None, anchor_tension: float = 2.0, anchor_lines: str = 'phantom',
                            show_anchor_vertices: bool = False, anchor_color: str = 'red', **_opts: Any) -> list[str]:

        constant_definitions['anchorTension'] = f'{anchor_tension:.2f}'
        constant_definitions['anchorLine'] = anchor_lines
        constant_definitions['anchorColor'] = anchor_color

        drawing: list[str] = []

        # Collect how many attachment point there is from the left and from the right
        left_attachment_points = utils.remove_duplicates(
            [p[-1] for _e_in, p in incoming_edges_paths])
        # Remove duplicates without changing order
        right_attachment_points = utils.remove_duplicates(
            [p[-1] for _e_out, p in outgoing_edges_paths])

        external_route: list[int] | None = None
        if len(left_attachment_points) == 1 and len(right_attachment_points) == 1:
            external_route = [left_attachment_points[0],
                              right_attachment_points[0]]
        elif len(left_attachment_points) == 0 and len(right_attachment_points) == 2:
            external_route = [right_attachment_points[0],
                              right_attachment_points[1]]
        elif len(left_attachment_points) == 2 and len(right_attachment_points) == 0:
            external_route = [left_attachment_points[0],
                              left_attachment_points[1]]
        elif len(left_attachment_points) + len(right_attachment_points) > 2:
            # No additional tension needed in this case
            return drawing

        if longest_cycle is None or len(longest_cycle) == 2:
            return drawing

        connections = {}
        if external_route is None:
            # Use surrounding anchors, but start the cycle with the first left attachment point if there is one right if there is not.
            if len(left_attachment_points) > 0:
                start = longest_cycle.index(left_attachment_points[0])
            elif len(right_attachment_points) > 0:
                start = longest_cycle.index(right_attachment_points[0])
            else:
                start = 0
            longest_cycle = longest_cycle[start:]+longest_cycle[:start]
            drawing.append(r'\fmfsurroundn{s}{%s}' % len(  # pylint: disable=consider-using-f-string'
                longest_cycle))
            for i_v in range(1, len(longest_cycle)+1):
                drawing.append(r'\%sVertexLabel{\fmflabel{$s_%d$}{s%d}}' % (  # pylint: disable=consider-using-f-string'
                    'show' if show_anchor_vertices else 'hide', i_v, i_v))
            for i_s, i_v in enumerate(longest_cycle):
                connections[self.vertices[i_v]] = f's{i_s+1}'
        else:
            index_left = longest_cycle.index(external_route[0])
            index_right = longest_cycle.index(external_route[1])
            # print('longest_cycle', longest_cycle)
            # print('shortest_path_left', shortest_path_left)
            # print('shortest_path_right', shortest_path_right)

            flipped = False
            if index_left > index_right:
                flipped = True
                index_left, index_right = index_right, index_left

            top_vertices = longest_cycle[index_left+1: index_right]
            bot_vertices = (longest_cycle[index_right +
                                          1:]+longest_cycle[:index_left])
            if flipped:
                top_vertices.reverse()
            else:
                bot_vertices.reverse()

            # print('top_vertices', top_vertices)
            # print('bot_vertices', bot_vertices)
            connections: dict[Vertex, str] = {}
            if len(top_vertices) > 0:
                for i_v, t_v in enumerate(top_vertices):
                    connections[self.vertices[t_v]] = f't{i_v+1}'
                drawing.append(r'\fmftopn{t}{%s}' % len(  # pylint: disable=consider-using-f-string'
                    top_vertices))
                for i_v in range(1, len(top_vertices)+1):
                    drawing.append(r'\%sAnchorLabel{\fmflabel{$t_%d$}{t%d}}' % (  # pylint: disable=consider-using-f-string'
                        'show' if show_anchor_vertices else 'hide', i_v, i_v))
            if len(bot_vertices) > 0:
                for i_v, b_v in enumerate(bot_vertices):
                    connections[self.vertices[b_v]] = f'b{i_v+1}'
                drawing.append(r'\fmfbottomn{b}{%s}' % len(  # pylint: disable=consider-using-f-string'
                    bot_vertices))
                for i_v in range(1, len(bot_vertices)+1):
                    drawing.append(r'\%sAnchorLabel{\fmflabel{$b_%d$}{b%d}}' % (  # pylint: disable=consider-using-f-string'
                        'show' if show_anchor_vertices else 'hide', i_v, i_v))
        for vertex, anchor in connections.items():
            drawing.append(r'\efmf{\anchorLine,tension=\anchorTension,fore=\anchorColor}{%s,%s}' % (  # pylint: disable=consider-using-f-string'
                anchor, f'v{self.get_vertex_position(vertex.name)}'))

        return drawing
