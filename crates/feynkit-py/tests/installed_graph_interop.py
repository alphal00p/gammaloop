"""Exercise HepKit and graph bindings in one installed Community extension.

The public API is provided by ``symbolica.community.hepkit``.
"""

import gc
import json
import weakref
import xml.etree.ElementTree as ET
from pathlib import Path

from symbolica.community import graph as linnet
from symbolica.community import hepkit as fk
from symbolica.community.tensor import TensorExpression

model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagrams = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_diagrams(loops=1, max_vertices=3, allow_self_loops=True)
    .diagrams
)
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
graph = diagram.to_graph()
assert type(graph) is linnet.Graph
assert type(graph).__module__ == "symbolica.community.graph"
assert graph is diagram.to_graph()
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

raw_full = graph.full_subgraph()
assert type(raw_full).__module__ == "symbolica.community.graph"
full = diagram.subgraph(raw_full)
empty = diagram.subgraph()
internal = diagram.filter(edge=lambda edge: not edge.is_external)
assert isinstance(full, fk.Subgraph)
assert isinstance(full, fk.FeynmanDiagram)
assert full.numerator_expression() == diagram.numerator_expression()
assert full.denominator_expression() == diagram.denominator_expression()
assert empty.denominator_expression() == TensorExpression(1)
assert empty.numerator_expression() == TensorExpression(1)
assert full.loop_count == diagram.loop_count
assert internal.loop_count == diagram.loop_count
external = diagram.filter(edge=lambda edge: edge.is_external)
external_vertices = {
    edge.source if edge.source is not None else edge.target
    for edge in diagram.external_edges
}
assert len(external.connected_components()) == len(external_vertices)
assert external.loop_count == 0
assert external.denominator_expression() == TensorExpression(1)
paired = next(edge for edge in graph.edges() if not edge.data.is_external)
boundary_half = diagram.subgraph(half_edges=[paired.source.data])
assert boundary_half.denominator_expression() == TensorExpression(1)
assert boundary_half.numerator_expression() == TensorExpression(1)
assert full.is_connected()
assert len(full.connected_components()) == 1
assert full.momentum_basis().loop_edges
assert full.loop_momentum_bases(limit=1)
assert len(diagram.loop_momentum_basis.momentum_replacements()) == len(diagram.edges)
assert (
    full.cross_free_family().to_expression()
    == diagram.cross_free_family().to_expression()
)

cycles, covered = internal.cycle_basis()
assert len(cycles) == diagram.loop_count
assert isinstance(covered, fk.Subgraph)
assert isinstance(full.bridges(), fk.Subgraph)
assert isinstance(internal.boundary(), fk.Subgraph)
assert internal.all_spanning_forests()
assert full.all_bonds()
assert diagram.all_cuts([0], [1])
assert isinstance(full.depth_first_traverse(0), fk.TraversalTree)
assert isinstance(full.breadth_first_traverse(0), fk.TraversalTree)


def native_ids(selection):
    return frozenset(half.data for half in selection.to_half_edges())


# Direct Rust analysis agrees with optional Python interoperability even when
# the exported half-edge numbering differs from the native diagram numbering.
assert full.bridges().half_edge_indices() == sorted(native_ids(graph.bridges(raw_full)))
assert {
    frozenset(component.half_edge_indices())
    for component in full.connected_components()
} == {native_ids(component) for component in graph.connected_components(raw_full)}
assert {
    frozenset(forest.half_edge_indices()) for forest in full.all_spanning_forests()
} == {native_ids(forest) for forest in graph.all_spanning_forests(raw_full)}
assert {frozenset(bond.half_edge_indices()) for bond in full.all_bonds()} == {
    native_ids(bond) for bond in graph.all_bonds(subgraph=raw_full)
}
assert internal.boundary().half_edge_indices() == sorted(
    native_ids(graph.boundary(internal.linnet_selection))
)
for half in graph.half_edges():
    assert diagram.subgraph(half_edges=[half.data]).linnet_selection == graph.subgraph(
        half_edges=[half.index]
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

for invalid in [restored.to_graph().full_subgraph(), object()]:
    try:
        diagram.subgraph(invalid)
    except (ValueError, TypeError):
        pass
    else:
        raise AssertionError("foreign selections must be rejected")

# Even an in-place edit that keeps the number of half-edges invalidates selections.
graph.reverse_edge(0)
try:
    diagram.subgraph(raw_full)
except (ValueError, ReferenceError):
    pass
else:
    raise AssertionError("stale selections must be rejected")
# Physics views retain immutable diagram ownership independently of analysis edits.
assert full.numerator_expression() == diagram.numerator_expression()
fresh = diagram.to_graph()
assert fresh is not graph
assert (
    diagram.subgraph(fresh.full_subgraph()).numerator_expression()
    == diagram.numerator_expression()
)

cross_section = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_cross_section(loops=1, max_vertices=2, allow_self_loops=True)
    .diagrams[0]
)
assert len(cross_section.vertices) == 2
assert all(not edge.is_dangling for edge in cross_section.external_edges)
assert len(cross_section.external_edges) == 1
cross_graph = cross_section.to_graph()
assert cross_graph.n_nodes == 2 and cross_graph.n_edges == 3
for cut in cross_section.cuts:
    assert cut.left.loop_count == cut.right.loop_count == 0
    assert len(cut.particles) == len(cut.edges) == 2
    assert all(particle.name == "scalar_0" for particle in cut.particles)
    assert isinstance(cut.left.subgraph, fk.Subgraph)
    assert isinstance(cut.right.subgraph, fk.Subgraph)
    assert cut.left.subgraph.numerator_expression() != TensorExpression(1)
    assert len(cut.propagators()) == 2
    assert len(cut.propagators(edge_powers={cut.edges[0].id: 2})) == 2
for candidate in cross_section.topology_threshold_candidates:
    assert isinstance(candidate.left, fk.Subgraph)
    assert isinstance(candidate.right, fk.Subgraph)
assert any(
    "is_cut" in json.loads(element.get("data-linnet-detail", "{}"))
    for element in ET.fromstring(cross_section.render().to_svg()).iter()
)
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
cycle_graph = cycle_diagram.to_graph()


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
