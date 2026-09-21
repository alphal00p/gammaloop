"""Exercise an installed FeynKit host with a separately built Linnet extension.

Pass the community module name as the first argument (``hep`` for that host).
"""

import gc
import importlib
from pathlib import Path
import sys
import weakref

import linnet_py

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = model.generate_diagrams(
    ["scalar_0"], ["scalar_0", "scalar_0"], loops=1
).diagrams
diagram = next(
    candidate
    for candidate in diagrams
    if len(
        {
            edge.source if edge.source is not None else edge.target
            for edge in candidate.external_edges
        }
    )
    > 1
)
graph = diagram.to_linnet()
assert type(graph) is linnet_py.Graph
assert graph is diagram.to_linnet()
assert graph.n_nodes == len(diagram.vertices)
assert graph.n_edges == len(diagram.edges)
assert all(vertex.interaction for vertex in diagram.vertices)
assert all(edge.is_dangling for edge in diagram.external_edges)
for edge in graph.edges():
    if edge.data.is_external:
        half = edge.source if edge.source is not None else edge.sink
        assert half.flow == (
            linnet_py.Flow.Sink
            if edge.data.external_state == "incoming"
            else linnet_py.Flow.Source
        )

assert all(not edge.is_external for edge in diagram.internal_edges)
assert sorted(half.data for half in graph.half_edges()) == list(
    range(graph.n_half_edges)
)

full = graph.full_subgraph()
empty = graph.empty_subgraph()
internal = graph.filter(edge=lambda edge: not edge.data.is_external)
assert diagram.numerator_expression(subgraph=full) == diagram.numerator_expression()
assert diagram.denominator_expression(subgraph=full) == diagram.denominator_expression()
assert diagram.denominator_expression(subgraph=empty) == 1
assert diagram.numerator_expression(subgraph=empty) == 1
assert diagram.loop_count_of(subgraph=full) == diagram.loop_count
assert diagram.loop_count_of(subgraph=internal) == diagram.loop_count
external = graph.filter(edge=lambda edge: edge.data.is_external)
external_vertices = {
    edge.source if edge.source is not None else edge.target
    for edge in diagram.external_edges
}
assert len(diagram.connected_components(external)) == len(external_vertices)
assert diagram.loop_count_of(subgraph=external) == 0
assert diagram.denominator_expression(subgraph=external) == 1
paired = next(edge for edge in graph.edges() if not edge.data.is_external)
boundary_half = graph.subgraph(half_edges=[paired.source.index])
assert diagram.denominator_expression(subgraph=boundary_half) == 1
assert diagram.numerator_expression(subgraph=boundary_half) == 1
assert diagram.is_connected(full)
assert len(diagram.connected_components(full)) == 1
assert diagram.momentum_basis(subgraph=full).loop_edges
assert diagram.loop_momentum_bases(limit=1, subgraph=full)
assert len(diagram.loop_momentum_basis.momentum_replacements()) == len(diagram.edges)
assert (
    diagram.build_cff(subgraph=full).to_expression()
    == diagram.build_cff().to_expression()
)

restored = fk.FeynmanDiagram.from_json(model, diagram.to_json())
assert restored.to_json() == diagram.to_json()
from_dot = fk.FeynmanDiagram.from_dot(model, diagram.to_dot())
assert from_dot.to_json() == diagram.to_json()
try:
    diagram.compatible_momentum_basis(restored.loop_momentum_basis)
except fk.DiagramError:
    pass
else:
    raise AssertionError("parent bases must retain diagram instance ownership")
assert diagram.compatible_momentum_basis(diagram.loop_momentum_basis)

for invalid in [restored.to_linnet().full_subgraph(), object()]:
    try:
        diagram.numerator_expression(subgraph=invalid)
    except (ValueError, TypeError):
        pass
    else:
        raise AssertionError("foreign selections must be rejected")

# Even an in-place edit that keeps the number of half-edges invalidates selections.
graph.reverse_edge(0)
try:
    diagram.numerator_expression(subgraph=full)
except ValueError:
    pass
else:
    raise AssertionError("stale selections must be rejected")
fresh = diagram.to_linnet()
assert fresh is not graph
assert (
    diagram.numerator_expression(subgraph=fresh.full_subgraph())
    == diagram.numerator_expression()
)

# A mutable analysis payload may point back at its owner; that cycle is collectable.
cycle_diagram = fk.FeynmanDiagram.from_json(model, diagram.to_json())
cycle_graph = cycle_diagram.to_linnet()
class Payload:
    pass


payload = Payload()
payload.diagram = cycle_diagram
cycle_graph.edge(0).data = payload
reference = weakref.ref(payload)
del cycle_diagram, cycle_graph, payload
gc.collect()
assert reference() is None
print("installed graph interoperability checks passed")
