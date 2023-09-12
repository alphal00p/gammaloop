from __future__ import annotations

import gammaloop.base_objects.model as model
from gammaloop.misc.common import *

# Abstract class for all vertex info
class VertexInfo(object):
    pass

class InteractonVertexInfo(VertexInfo):
    def __init__(self, vertex_rule: model.VertexRule):
        self.vertex_rule: model.VertexRule = vertex_rule

class Vertex(object):
    def __init__(self, name: str, vertex_info: VertexInfo, edges: list[Edge] | None):
        self.name: name = name
        self.vertex_rule: VertexInfo = vertex_info
        self.edges: list[Edge] = edges

class Edge(object):
    def __init__(self, name: str, particle: model.Particle, vertices: (Vertex, Vertex) | None):
        self.id: name = name
        self.particle: model.Particle = particle
        self.vertices: (Vertex, Vertex) = vertices

# TODO create serializable counterpart and export to yaml, similarly as done for the model.

class Graph(object):
    def __init__(self, name: str, vertices: list[Vertex], edges: list[Edge]):
        self.name: str = name
        self.vertices: list[Vertex] = vertices
        self.edges: list[Edge] = edges
        self.name_to_position: dict[str,dict[str,int]] = {}
    
    def get_edge(self, edge_name: str) -> Edge:
        return self.edges[self.name_to_position['edges'][edge_name]]
    
    def get_vertex(self, vertex_name: str) -> Parameter:
        return self.vertices[self.name_to_position['vertices'][vertex_name]]

    def synchronize_name_map(self) -> None:
        for i_e, edge in enumerate(self.edges):
            self.name_to_position['edges'][edge.name] = i_e
        for i_v, vertex in enumerate(self.vertices):
            self.name_to_position['vertices'][vertex.name] = i_v


    @staticmethod
    def from_qgraph(model: model.Model, qgraph_object) -> Graph:
        """ Imports graph form a stylicized qgraph python output file. Will be deprecated when using in-house gammaloop graph generation. """
        logger.info("Importing graph from qgraph python output")
        # TODO
        return Graph('default',[],[])