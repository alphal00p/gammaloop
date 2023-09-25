from __future__ import annotations

from enum import StrEnum
import yaml

import gammaloop.base_objects as base_objects
from gammaloop.misc.common import GammaLoopError
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


class InteractonVertexInfo(VertexInfo):
    def __init__(self, vertex_rule: VertexRule):
        self.vertex_rule: VertexRule = vertex_rule

    def to_serializable_dict(self) -> dict:
        return {
            'type': 'interacton_vertex_info',
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

    def to_serializable_dict(self) -> dict:
        return {
            'type': 'external_vertex_info',
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


class Edge(object):
    def __init__(self, name: str, edge_type: EdgeType, particle: Particle, vertices: (Vertex, Vertex) | None = None):
        self.name: str = name
        self.edge_type: EdgeType = edge_type
        self.particle: Particle = particle
        self.vertices: (Vertex, Vertex) | None = vertices

    def to_serializable_dict(self) -> dict:
        return {
            'name': self.name,
            'edge_type': self.edge_type.value,
            'particle': self.particle.name,
            'vertices': [self.vertices[0].name, self.vertices[1].name]
        }


class Graph(object):
    def __init__(self, name: str, vertices: list[Vertex], edges: list[Edge], external_connections: list[(Vertex | None, Vertex | None)], overall_factor: float = 1.0):
        self.name: str = name
        self.vertices: list[Vertex] = vertices
        self.edges: list[Edge] = edges
        self.overall_factor = overall_factor
        # For forward scattering graphs, keep track of the bipartite map, i.e. which in and out externals will carry identical momenta.
        self.external_connections: list[(
            Vertex | None, Vertex | None)] = external_connections
        self.name_to_position: dict[str, dict[str, int]] = {}

    @staticmethod
    def empty_graph(name: str) -> Graph:
        return Graph(name, [], [], [], 1.0)

    def get_edge(self, edge_name: str) -> Edge:
        return self.edges[self.name_to_position['edges'][edge_name]]

    def get_vertex(self, vertex_name: str) -> Parameter:
        return self.vertices[self.name_to_position['vertices'][vertex_name]]

    def get_incoming_edges(self) -> list[Edge]:
        return [e for e in self.edges if e.edge_type == EdgeType.INCOMING]

    def get_outgoing_edges(self) -> list[Edge]:
        return [e for e in self.edges if e.edge_type == EdgeType.OUTGOING]

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
            'overall_factor': self.overall_factor
        }

    @staticmethod
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

        graph = Graph(serializable_dict['name'], graph_vertices, graph_edges,
                      external_connections, serializable_dict['overall_factor'])
        graph.synchronize_name_map()

        return graph

    @staticmethod
    def from_yaml_str(model: Model, yaml_str: str) -> Graph:
        graph_dict = yaml.safe_load(yaml_str)
        return Graph.from_serializable_dict(model, graph_dict)

    def to_yaml_str(self) -> str:
        return utils.verbose_yaml_dump(self.to_serializable_dict())

    @staticmethod
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

        return graph
