"""Exercise an installed FeynKit host with a separately built Linnet extension.

Pass the community module name as the first argument (``hep`` for that host).
"""

import gc
import importlib
import sys
import weakref
from pathlib import Path

import linnet

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    loops=1,
    max_vertices=3,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=True,
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
assert type(graph) is linnet.Graph
assert type(graph).__module__ == "linnet"
assert graph is diagram.to_linnet()
assert graph.n_nodes == len(diagram.vertices)
assert graph.n_edges == len(diagram.edges)
assert all(vertex.interaction for vertex in diagram.vertices)
assert all(edge.is_dangling for edge in diagram.external_edges)
for edge in graph.edges():
    if edge.data.is_external:
        half = edge.source if edge.source is not None else edge.sink
        assert half.flow == (
            linnet.Flow.Sink
            if edge.data.external_state == "incoming"
            else linnet.Flow.Source
        )

assert all(not edge.is_external for edge in diagram.internal_edges)
assert sorted(half.data for half in graph.half_edges()) == list(
    range(graph.n_half_edges)
)

full = graph.full_subgraph()
assert type(full).__module__ == "linnet"
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

cycles, covered = diagram.cycle_basis(internal)
assert len(cycles) == diagram.loop_count
assert isinstance(covered, linnet.Subgraph)
assert isinstance(diagram.bridges(full), linnet.Subgraph)
assert isinstance(diagram.boundary(internal), linnet.Subgraph)
assert diagram.all_spanning_forests(internal)
assert diagram.all_bonds(subgraph=full)
assert diagram.all_cuts([0], [1])
assert isinstance(diagram.depth_first_traverse(0, subgraph=full), linnet.TraversalTree)
assert isinstance(
    diagram.breadth_first_traverse(0, subgraph=full), linnet.TraversalTree
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
except (ValueError, ReferenceError):
    pass
else:
    raise AssertionError("stale selections must be rejected")
fresh = diagram.to_linnet()
assert fresh is not graph
assert (
    diagram.numerator_expression(subgraph=fresh.full_subgraph())
    == diagram.numerator_expression()
)

cross_section = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    kind="cross_section",
    loops=1,
    max_vertices=2,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=True,
).diagrams[0]
assert len(cross_section.vertices) == 2
assert all(not edge.is_dangling for edge in cross_section.external_edges)
assert len(cross_section.external_edges) == 1
cross_graph = cross_section.to_linnet()
assert cross_graph.n_nodes == 2 and cross_graph.n_edges == 3
for cut in cross_section.cuts:
    assert cut.left.loop_count == cut.right.loop_count == 0
    assert len(cut.particles) == len(cut.edges) == 2
    assert all(particle.name == "scalar_0" for particle in cut.particles)
    assert isinstance(cut.left.subgraph, linnet.Subgraph)
    assert isinstance(cut.right.subgraph, linnet.Subgraph)
    assert cross_section.numerator_expression(subgraph=cut.left.subgraph) != 1
    assert len(cut.propagators()) == 2
    assert len(cut.propagators(edge_powers={cut.edges[0].id: 2})) == 2
for candidate in cross_section.topology_threshold_candidates:
    assert isinstance(candidate.left, linnet.Subgraph)
    assert isinstance(candidate.right, linnet.Subgraph)
assert "is_cut:" in cross_section.to_linnest()
assert (
    fk.FeynmanDiagram.from_json(model, cross_section.to_json()).to_json()
    == cross_section.to_json()
)
assert (
    fk.FeynmanDiagram.from_dot(model, cross_section.to_dot()).to_json()
    == cross_section.to_json()
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
