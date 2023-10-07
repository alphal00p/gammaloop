from __future__ import annotations

import os
import itertools
from pathlib import Path
from enum import StrEnum
import yaml

import gammaloop.base_objects as base_objects
from gammaloop.misc.common import GammaLoopError, DATA_PATH, pjoin, logger  # pylint: disable=unused-import
import gammaloop.misc.utils as utils
from gammaloop.base_objects.model import Model, VertexRule, Particle, Parameter


class EdgeType(StrEnum):
    INCOMING = 'in'
    OUTGOING = 'out'
    VIRTUAL = 'virtual'

# Abstract class for all vertex info


class VertexInfo(object):

    @staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict) -> VertexInfo:
        if serializable_dict['type'] == 'interacton_vertex_info':
            return InteractonVertexInfo.from_serializable_dict(model, serializable_dict)
        elif serializable_dict['type'] == 'external_vertex_info':
            return ExternalVertexInfo.from_serializable_dict(model, serializable_dict)
        else:
            raise GammaLoopError(
                f"Unknown vertex info type: {serializable_dict['type']}")

    def get_particles(self):
        raise NotImplementedError()

    def is_external(self):
        return False


class InteractonVertexInfo(VertexInfo):
    def __init__(self, vertex_rule: VertexRule):
        self.vertex_rule: VertexRule = vertex_rule

    def get_type(self) -> str:
        return 'interacton_vertex_info'

    def to_serializable_dict(self) -> dict:
        return {
            'type': self.get_type(),
            'vertex_rule': self.vertex_rule.name
        }

    @staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict) -> InteractonVertexInfo:
        return InteractonVertexInfo(model.get_vertex_rule(serializable_dict['vertex_rule']))

    def get_particles(self) -> list[Particle]:
        return self.vertex_rule.particles


class ExternalVertexInfo(VertexInfo):
    def __init__(self, direction: EdgeType, particle: base_objects.model.Particle):

        self.direction: EdgeType = direction
        self.particle: base_objects.model.Particle = particle

    def is_external(self):
        return True

    def get_type(self) -> str:
        return 'external_vertex_info'

    def to_serializable_dict(self) -> dict:
        return {
            'type': self.get_type(),
            'particle': self.particle.name,
            'direction': self.direction.value,
        }

    @staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict) -> ExternalVertexInfo:
        return ExternalVertexInfo(
            EdgeType(serializable_dict['direction']),
            model.get_particle(serializable_dict['particle'])
        )

    def get_particles(self) -> list[Particle]:
        return [self.particle,]


class Vertex(object):
    def __init__(self, name: str, vertex_info: VertexInfo | None = None, edges: list[Edge] | None = None):
        self.name: str = name
        self.vertex_info: VertexInfo | None = vertex_info
        self.edges: list[Edge] | None = edges

    def set_vertex_rule_from_model(self, model: Model) -> None:
        if len(self.edges) == 1:
            self.vertex_info = ExternalVertexInfo(
                EdgeType.OUTGOING if self.edges[0].edge_type == EdgeType.OUTGOING else EdgeType.INCOMING, self.edges[0].particle)
            return

        if self.vertex_info is None:
            vertex_particles = [e.particle if e.vertices[1] is self else e.particle.get_anti_particle(
                model) for e in self.edges]
            orig_vertex_particles = list(vertex_particles)

            interactions = model.get_vertices_from_particles(vertex_particles)
            if len(interactions) == 0:
                raise GammaLoopError(
                    f"No interaction found for vertex {self.name} with particles [{','.join(p.name for p in orig_vertex_particles)}]")
            elif len(interactions) > 1:
                raise GammaLoopError(
                    f"Found more than one valid interaction for vertex {self.name} with particles [{','.join(p.name for p in orig_vertex_particles)}]. Importing from qgraph is not supported in this case.")
            self.vertex_info = InteractonVertexInfo(interactions[0])
        elif not isinstance(self.vertex_info, InteractonVertexInfo):
            raise GammaLoopError(
                "Vertex rule already set to a different type.")

        self.sort_edges_according_to_vertex_info(model)

    def sort_edges_according_to_vertex_info(self, model: Model) -> None:
        vertex_particles = self.vertex_info.get_particles()
        if len(vertex_particles) <= 1:
            return
        sorted_edges = []
        vertex_particles = [(e.particle, e) if e.vertices[1] is self else (
            e.particle.get_anti_particle(model), e) for e in self.edges]
        orig_vertex_particles = [v[0] for v in vertex_particles]
        for part in orig_vertex_particles:
            for i, (e_p, _edge) in enumerate(vertex_particles):
                if e_p == part:
                    sorted_edges.append(vertex_particles.pop(i)[1])
                    break
            else:
                raise GammaLoopError(
                    f"Particle {part[0].name} not found in vertex {self.name} with particles [{','.join(p.name for p in orig_vertex_particles)}].")
        if len(vertex_particles) > 0:
            raise GammaLoopError(
                f"Not all particles were found in vertex {self.name} with particles [{','.join(p.name for p in orig_vertex_particles)}].")
        self.edges = sorted_edges

    def to_serializable_dict(self) -> dict:
        return {
            'name': self.name,
            'vertex_info': self.vertex_info.to_serializable_dict(),
            'edges': [e.name for e in self.edges]
        }

    def draw(self, graph, _model: Model, constant_definitions: dict[str, str], show_vertex_labels=False,
             use_vertex_names=False, vertex_size=2.5, vertex_shape='circle', **_opts) -> list[str]:
        # Possible shape choices are 'circle', 'square', 'triangle', 'diamond', 'pentagon', 'hexagon', 'triagram',
        # 'tetragram', 'pentagram', 'hexagram', 'triacross', 'cross', 'pentacross' and 'hexacross'.

        constant_definitions['vertexShape'] = vertex_shape
        constant_definitions['vertexSize'] = f'{vertex_size:.2f}thick'
        constant_definitions['blobSize'] = f'{vertex_size*3:.2f}'

        drawing = []
        if self.vertex_info.is_external():
            return drawing

        vertex_pos = graph.get_vertex_position(self.name)
        prefix = r'\showVertexLabel{'
        if not show_vertex_labels:
            prefix = r'\hideVertexLabel{'
        drawing.append(r'%s\fmflabel{%s}{v%d}}' % (  # pylint: disable=consider-using-f-string
            prefix, self.name if use_vertex_names else r'$v_{%s}$' % vertex_pos, vertex_pos))  # pylint: disable=consider-using-f-string

        # One could tag the form factors vertices with blobs here
        is_blob = False
        if not is_blob:
            drawing.append(
                r'\efmfv{decor.shape=\vertexShape,decor.filled=full,decor.size=\vertexSize}{v%d}' % (  # pylint: disable=consider-using-f-string
                    vertex_pos))
        else:
            drawing.append(
                r'\efmfv{decor.shape=\vertexShape,decor.filled=shaded,decor.size=\blobSize}{v%d}' % (  # pylint: disable=consider-using-f-string
                    vertex_pos))

        return drawing


class Edge(object):
    def __init__(self, name: str, edge_type: EdgeType, particle: Particle, vertices: (Vertex, Vertex) | None = None):
        self.name: str = name
        self.edge_type: EdgeType = edge_type
        self.particle: Particle = particle
        self.vertices: (Vertex, Vertex) | None = vertices

    def get_vertex_positions(self, graph: Graph) -> (int, int):
        return (
            graph.get_vertex_position(self.vertices[0].name),
            graph.get_vertex_position(self.vertices[1].name),
        )

    def to_serializable_dict(self) -> dict:
        return {
            'name': self.name,
            'edge_type': self.edge_type.value,
            'particle': self.particle.name,
            'vertices': [self.vertices[0].name, self.vertices[1].name]
        }

    def draw(self, graph: Graph, _model: Model, constant_definitions: dict[str, str], show_edge_labels=True, show_particle_names=True, show_edge_names=True,
             show_edge_momenta=True, draw_arrow_for_all_edges=True, draw_lmb=True, arc_max_distance=1., external_legs_tension=3., default_tension=1.0, line_width=1.0,
             line_color='black', label_color='blue', lmb_color='red', **_opts) -> list[str]:

        constant_definitions['arcMaxDistance'] = f'{arc_max_distance:.2f}'
        constant_definitions['externalLegTension'] = f'{external_legs_tension:.2f}'
        constant_definitions['defaultTension'] = f'{default_tension:.2f}'
        constant_definitions['lineWidth'] = f'{line_width:.2f}'
        constant_definitions['lineColor'] = line_color
        constant_definitions['lmbColor'] = lmb_color
        constant_definitions['labelColor'] = label_color

        template_line = r'\efmf{%(line_type)s%(comma)s%(options)s}{%(left_vertex)s,%(right_vertex)s}'
        lmb_position = graph.loop_momentum_basis.index(
            self) if self in graph.loop_momentum_basis else None
        line_options = {}

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

        if self.particle.is_massive() and abs(self.particle.spin) != 3:
            line_type = f'dbl_{line_type}'

        line_options['width'] = r'\lineWidth'
        line_options['fore'] = r'\lineColor'
        line_options['tension'] = r'\defaultTension'

        edge_name_label = ''
        if show_edge_names:
            edge_name_label = self.name  # pylint: disable=consider-using-f-string
        label_components = []
        if show_particle_names:
            label_components.append(self.particle.texname)

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

        line_options['label'] = '%s%s$%s%s$' % (  # pylint: disable=consider-using-f-string
            r"\setLabelColor{\labelColor}" if lmb_position is None else r"\setLabelColor{\lmbColor}",
            edge_name_label,
            '|' if edge_name_label != '' and len(label_components) > 0 else '',
            '|'.join(label_components)
        )

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
        elif len(multi_edges) == 3:
            match self_position:
                case 0:
                    line_options[base_side] = str(0.4)
                case 1:
                    pass
                case 2:
                    line_options[other_side] = str(0.4)
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
        replace_dict['line_type'] = line_type
        replace_dict['left_vertex'] = f'v{graph.get_vertex_position(self.vertices[0].name)}'
        replace_dict['right_vertex'] = f'v{graph.get_vertex_position(self.vertices[1].name)}'
        replace_dict['options'] = ','.join(
            f'{k}={v}' for k, v in line_options.items())
        replace_dict['comma'] = ',' if len(replace_dict['options']) > 0 else ''
        drawing = []
        drawing.append(
            f"% == Drawing edge '{self.name}' (e{graph.get_edge_position(self.name)}) with particle {self.particle.name}")
        drawing.append(template_line % replace_dict)

        if 'arrow' not in line_type and draw_arrow_for_all_edges:
            if 'label' in line_options:
                del line_options['label']
            replace_dict['line_type'] = 'phantom_arrow'
            line_options['tension'] = 0
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
            replace_dict['line_type'] = 'phantom_arrow'
            line_options['tension'] = 0
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
    def __init__(self, name: str, vertices: list[Vertex], edges: list[Edge], external_connections: list[(Vertex | None, Vertex | None)],
                 loop_momentum_basis: list[Edge] | None = None, overall_factor: float = 1.0, edge_signatures: list[(list[int], list[int])] | None = None):
        self.name: str = name
        self.vertices: list[Vertex] = vertices
        self.edges: list[Edge] = edges
        self.edge_signatures: list[(list[int], list[int])] = edge_signatures
        self.overall_factor = overall_factor
        # For forward scattering graphs, keep track of the bipartite map, i.e. which in and out externals will carry identical momenta.
        self.external_connections: list[(
            Vertex | None, Vertex | None)] = external_connections
        self.loop_momentum_basis: list[Edge] = loop_momentum_basis
        self.name_to_position: dict[str, dict[str, int]] = {}

    def get_sorted_incoming_edges(self):
        sorted_incoming_edges = []
        for (u, _v) in self.external_connections:  # pylint: disable=invalid-name
            if u is None:
                continue
            for e in self.get_incoming_edges():  # pylint: disable=invalid-name
                if e.vertices[0] == u:
                    sorted_incoming_edges.append(e)
                    break
        return sorted_incoming_edges

    def get_sorted_outgoing_edges(self):
        sorted_outgoing_edges = []
        for (_u, v) in self.external_connections:  # pylint: disable=invalid-name
            if v is None:
                continue
            for e in self.get_outgoing_edges():  # pylint: disable=invalid-name
                if e.vertices[1] == v:
                    sorted_outgoing_edges.append(e)
                    break
        return sorted_outgoing_edges

    def get_edge_signature(self, edge_name):
        return self.edge_signatures[self.get_edge_position(edge_name)]

    @ staticmethod
    def empty_graph(name: str) -> Graph:
        return Graph(name, [], [], [], loop_momentum_basis=[], overall_factor=1.0, edge_signatures=[])

    def get_edge(self, edge_name: str) -> Edge:
        return self.edges[self.name_to_position['edges'][edge_name]]

    def get_vertex(self, vertex_name: str) -> Parameter:
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

    def to_serializable_dict(self) -> dict:
        return {
            'name': self.name,
            'vertices': [v.to_serializable_dict() for v in self.vertices],
            'edges': [e.to_serializable_dict() for e in self.edges],
            'external_connections': [[e[0].name if e[0] is not None else None,
                                      e[1].name if e[1] is not None else None] for e in self.external_connections],
            'overall_factor': self.overall_factor,
            'loop_momentum_basis': [e.name for e in self.loop_momentum_basis],
            'edge_signatures': self.edge_signatures
        }

    @ staticmethod
    def from_serializable_dict(model: Model, serializable_dict: dict) -> Graph:

        graph_vertices = []
        vertex_name_to_vertex_map = {}
        for v_graph in serializable_dict['vertices']:
            v_name = v_graph['name']
            vert = Vertex(v_name)
            vertex_name_to_vertex_map[v_name] = vert
            graph_vertices.append(vert)

        graph_edges = []
        edge_name_to_edge_map = {}
        for e_qgraph in serializable_dict['edges']:
            e_name = e_qgraph['name']
            edge = Edge(e_name, EdgeType(e_qgraph['edge_type']), model.get_particle(e_qgraph['particle']),
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

        loop_momentum_basis = [edge_name_to_edge_map[e_name]
                               for e_name in serializable_dict['loop_momentum_basis']]
        graph = Graph(serializable_dict['name'], graph_vertices, graph_edges,
                      external_connections, serializable_dict['overall_factor'],
                      loop_momentum_basis, serializable_dict['edge_signatures'])
        graph.synchronize_name_map()

        return graph

    @ staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> Graph:
        graph_dict = yaml.safe_load(yaml_str)
        return Graph.from_serializable_dict(model, graph_dict)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump(self.to_serializable_dict())

    @ staticmethod
    def from_qgraph(model: Model, qgraph_object, name: str = 'default') -> Graph:
        """ Imports graph form a stylicized qgraph python output file. Will be deprecated when using in-house gammaloop graph generation. """

        graph_vertices = []
        vertex_qgraph_index_to_vertex_map = {}
        qgraph_vertices = list(qgraph_object['nodes'].items())
        for i_vertex, (v_qgraph_index, v_qgraph) in enumerate(qgraph_vertices):
            vert = Vertex(f"v{i_vertex+1}")
            vertex_qgraph_index_to_vertex_map[v_qgraph_index] = vert
            graph_vertices.append(vert)

        graph_edges = []
        edge_qgraph_index_to_edge_map = {}
        qgraph_edges = list(qgraph_object['edges'].items())
        for e_qgraph_index, e_qgraph in qgraph_edges:
            edge = Edge(e_qgraph['name'], EdgeType(e_qgraph['type']), model.get_particle_from_pdg(e_qgraph['PDG']),
                        (vertex_qgraph_index_to_vertex_map[e_qgraph['vertices'][0]],
                        vertex_qgraph_index_to_vertex_map[e_qgraph['vertices'][1]])
                        )
            edge_qgraph_index_to_edge_map[e_qgraph_index] = edge
            graph_edges.append(edge)

        for (vert, (v_qgraph_index, v_qgraph)) in zip(graph_vertices, qgraph_vertices):
            vert.edges = [edge_qgraph_index_to_edge_map[e_qgraph_index]
                          for e_qgraph_index in v_qgraph['edge_ids']]
            vert.set_vertex_rule_from_model(model)

        external_connections = {}
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

        external_connections = [v for k, v in sorted(
            list(external_connections.items()), key=lambda el: el[0])]

        graph = Graph(name, graph_vertices, graph_edges,
                      external_connections, float(qgraph_object['overall_factor']))
        graph.synchronize_name_map()

        all_lmbs = graph.generate_loop_momentum_bases(
            target_n_lmbs_to_generate=1)
        if len(all_lmbs) == 0:
            raise GammaLoopError(
                f"Could not find any valid loop momentum basis for graph '{name}'.")
        else:
            # Arbitrarily pick the first valid momentum basis as the chosen one.
            graph.loop_momentum_basis = all_lmbs[0]

        graph.edge_signatures = graph.generate_momentum_flow()

        return graph

    def get_edge_map(self) -> list[(int, int)]:
        result = []
        for edge in self.edges:
            result.append(
                (
                    self.name_to_position['vertices'][edge.vertices[0].name],
                    self.name_to_position['vertices'][edge.vertices[1].name]
                )
            )
        return result

    def get_adjacency_map(self, edge_map: list[(int, int)]) -> dict[int, list[((int, (int, int)), int)]]:
        result = {}
        for i_e, edge in enumerate(edge_map):
            result.setdefault(edge[0], []).append(((i_e, edge), edge[1]))
            result.setdefault(edge[1], []).append(((i_e, edge), edge[0]))
        return result

    def generate_loop_momentum_bases(self, target_n_lmbs_to_generate: int | None = None) -> list[list[Edge]]:
        """ Return all valid loop momentum bases for the graph"""

        # logger.debug("Generating loop momentum bases for graph %s", self.name)

        # This function need to be called for each disconnected subgraph of self.
        vertices_not_encountered: list[list[tuple[int]]] = list(
            range(len(self.vertices)))
        all_spanning_tree_lists_from_disconnected_components = []
        while len(vertices_not_encountered) > 0:
            seed_vertex_position = vertices_not_encountered[0]
            seen_edge_sequences: set[int] = set()
            accumulated_edge_sequence: list[int] = []
            excluded_edges: list[int] = []
            seen_edge_sequences = None
            spanning_trees_for_this_subgraph: list[tuple[int]] = []
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
        spanning_trees = [tuple(sum((list(st) for st in st_combination), [])) for st_combination in
                          itertools.product(*all_spanning_tree_lists_from_disconnected_components)]

        # logger.debug("At total of %d loop momentum bases were found for graph %s", len(
        #    spanning_trees), self.name)

        return [[e for e_i, e in enumerate(self.edges) if e_i not in spanning_tree] for spanning_tree in spanning_trees]

    def generate_momentum_flow(self) -> list[(list[int], list[int])]:
        if self.loop_momentum_basis is None:
            raise GammaLoopError(
                "Specify a loop momentum basis before generating a momentum flow.")
        if len(self.external_connections) == 0:
            sink_edge = None
        else:
            if self.external_connections[-1][1] is not None:
                sink_edge = self.get_edge_position(
                    self.get_sorted_outgoing_edges()[-1].name)
            else:
                sink_edge = self.get_edge_position(
                    self.get_sorted_incoming_edges()[-1].name)
        external_edges = []
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
        merged_signatures = []
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

    def draw(self, model: Model, file_path_without_extension=str | None, caption=None, **drawing_options) -> Path:
        match drawing_options['mode']:
            case 'feynmp':
                return self.draw_feynmp(model, file_path_without_extension, caption, **drawing_options['feynmp'])
            case _:
                raise GammaLoopError(
                    "Feynman drawing mode '%d' is not supported. Currently only 'feynmp' is supported." % drawing_options['mode'])

    def draw_feynmp(self, model: Model, file_path_without_extension=str | None, caption=None, caption_size='20pt', **drawing_options) -> Path:

        with open(pjoin(DATA_PATH, 'templates', 'drawing', 'drawing.tex'), 'r', encoding='utf-8') as file:
            template = file.read()
        source_file_path = Path(file_path_without_extension+'.tex')

        constant_definitions = {}

        with open(source_file_path, 'w', encoding='utf-8') as file:
            replace_dict = {}
            replace_dict['drawing_name'] = os.path.basename(
                file_path_without_extension)
            if len(self.get_incoming_edges()) > 0 and len(self.get_outgoing_edges()) > 0:
                replace_dict['incoming_vertices'] = r"\fmfleft{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[0].name)}' for e in self.get_incoming_edges())
                replace_dict['outgoing_vertices'] = r"\fmfright{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[1].name)}' for e in self.get_outgoing_edges())
            elif len(self.get_incoming_edges()) > 0:
                replace_dict['incoming_vertices'] = r"\fmfsurround{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[0].name)}' for e in self.get_incoming_edges())
                replace_dict['outgoing_vertices'] = "% None."
            elif len(self.get_outgoing_edges()) > 0:
                replace_dict['outgoing_vertices'] = r"\fmfsurround{%s}" % ','.join(  # pylint: disable=consider-using-f-string
                    f'v{self.get_vertex_position(e.vertices[0].name)}' for e in self.get_incoming_edges())
                replace_dict['incoming_vertices'] = "% None."
            else:
                # The anchor tension will do the job of spreading the graph open
                pass

            edge_map = self.get_edge_map()
            longest_cycle_edges = utils.find_longest_cycle(edge_map)
            longest_cycle = utils.edges_cycle_to_vertex_cycle([
                self.edges[i_e].get_vertex_positions(self) for i_e in longest_cycle_edges])
            incoming_edges_paths: list[(Edge, tuple[int])] = []
            for edge in self.get_sorted_incoming_edges():
                incoming_edges_paths.append((edge, utils.find_shortest_path(
                    edge_map, self.get_vertex_position(
                        edge.vertices[0].name),
                    longest_cycle)))

            outgoing_edges_paths: list[(Edge, tuple[int])] = []
            for edge in self.get_sorted_outgoing_edges():
                outgoing_edges_paths.append((edge, utils.find_shortest_path(
                    edge_map, self.get_vertex_position(
                        edge.vertices[1].name),
                    longest_cycle)))

            replace_dict['incoming_edge_definitions'] = '\n'.join(
                sum((edge.draw(self, model, constant_definitions, **drawing_options) for edge, _ in incoming_edges_paths), []))
            replace_dict['outgoing_edge_definitions'] = '\n'.join(
                sum((edge.draw(self, model, constant_definitions, **drawing_options) for edge, _ in outgoing_edges_paths), []))
            internal_edges_drawing = []
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
                sum((vertex.draw(self, model, constant_definitions, **drawing_options) for vertex in self.vertices), []))

            constant_definitions['captionSize'] = caption_size
            if caption:
                replace_dict['caption'] = caption
            else:
                replace_dict['caption'] = replace_dict['drawing_name']
            replace_dict['caption'] = replace_dict['caption'].replace(
                '_', ' ').replace('#', r'\#')

            replace_dict['constant_definitions'] = '\n'.join(
                r'\newcommand{\%s}{%s}' % (k, v) for k, v in constant_definitions.items()  # pylint: disable=consider-using-f-string
            )

            file.write(template % replace_dict)

        return source_file_path

    def draw_anchor_tension(self, incoming_edges_paths: list[Edge, tuple[int]], outgoing_edges_paths: list[Edge, tuple[int]], constant_definitions: dict[str, str],
                            longest_cycle: list[int], anchor_tension: float = 2.0, anchor_lines='phantom',
                            show_anchor_vertices=False, anchor_color='red', **_opts) -> list[str]:

        constant_definitions['anchorTension'] = f'{anchor_tension:.2f}'
        constant_definitions['anchorLine'] = anchor_lines
        constant_definitions['anchorColor'] = anchor_color

        drawing = []

        # Collect how many attachment point there is from the left and from the right
        left_attachment_points = utils.remove_duplicates(
            [p[-1] for e_in, p in incoming_edges_paths])
        # Remove duplicates without changing order
        right_attachment_points = utils.remove_duplicates(
            [p[-1] for e_out, p in outgoing_edges_paths])

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
