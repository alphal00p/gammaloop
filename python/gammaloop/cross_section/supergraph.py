from __future__ import annotations
from gammaloop.base_objects.graph import Graph, Edge, Vertex

class SuperGraphCut(object):
    def __init__(self, cut_edges: list[Edges], forward_scattering_graph: ForwardScatteringGraph):
        self.cut_edges: list[Edges] = cut_edges
        self.forward_scattering_graph: ForwardScatteringGraph = forward_scattering_graph

class SuperGraph(object):
    def __init__(self, id: int, graph: Graph, cuts: list[SuperGraphCut]):
        self.id: int = id
        self.graph: Graph = graph
        self.cuts: list[SuperGraphCut] = cuts

class ForwardScatteringGraphCut(object):
    def __init__(self, cut_edges: list[Edges], amplitudes: (Amplitude, Amplitude)):
        self.cut_edges: list[Edges] = cut_edges
        self.amplitudes: (Amplitude, Amplitude) = amplitudes

class ForwardScatteringGraph(object):
    def __init__(self, id: (int, int), graph: Graph, cuts: list[ForwardScatteringGraphCut]):
        # ID format (SG_id, SG_cut_id)
        self.id: (int, int) = id
        self.graph: Graph = graph
        self.cuts: list[ForwardScatteringGraphCut] = cuts

class Amplitude(object):
    def __init__(self, id: (int, int, int, int), graph: Graph):
        # ID format (SG_id, SG_cut_id, FS_cut_id, cut_side)
        self.id: (int, int, int, int) = id
        self.graph: Graph = graph