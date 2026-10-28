"""Compare edge symmetries for one directed or undirected self-loop.

Run with the Symbolica build under investigation. No FeynKit import is needed.
"""

from symbolica import Graph

for attached_external in (False, True):
    for directed in (False, True):
        graph = Graph()
        vertex = graph.add_node(1)
        if attached_external:
            graph.add_edge(vertex, graph.add_node(2), data=3)
        graph.add_edge(vertex, vertex, directed=directed, data=4)
        _, _, symmetry, _ = graph.canonize()
        print(f"external={attached_external}, directed={directed}: {symmetry}")
