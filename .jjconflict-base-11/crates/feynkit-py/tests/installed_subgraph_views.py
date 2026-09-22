"""Check physics subgraph inheritance, ownership, and explicit excision."""

import gc
import importlib
import inspect
import sys
from pathlib import Path

import linnet

fk = importlib.import_module(
    f"symbolica.community.{sys.argv[1] if len(sys.argv) > 1 else 'feynkit'}"
)
model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagram = next(
    candidate
    for candidate in model.generate_diagrams(
        ["scalar_0"],
        ["scalar_0", "scalar_0"],
        loops=1,
        max_vertices=3,
        vertex_allow=["V_3_SCALAR_000"],
        allow_self_loops=False,
    ).diagrams
    if len(
        {
            edge.source if edge.source is not None else edge.target
            for edge in candidate.external_edges
        }
    )
    == 3
)
snapshot = diagram.to_json()
graph = diagram.to_linnet()
raw_full = graph.full_subgraph()
full = diagram.subgraph(raw_full)
internal = diagram.filter(edge=lambda edge: not edge.data.is_external)
external = diagram.filter(edge=lambda edge: edge.data.is_external)
empty = diagram.subgraph()

assert issubclass(fk.Subgraph, fk.FeynmanDiagram)
assert type(internal) is fk.Subgraph
assert isinstance(internal, fk.FeynmanDiagram)
assert internal.to_linnet() is graph
assert type(internal.linnet_selection) is linnet.Subgraph
assert internal.original.to_json() == snapshot
assert internal.loop_count == diagram.loop_count
assert empty.loop_count == 0
assert empty.numerator_expression() == empty.denominator_expression() == 1
assert full.numerator_expression() == diagram.numerator_expression()
assert full.denominator_expression() == diagram.denominator_expression()

# Nested selection and selection algebra retain the immutable physics owner.
assert internal.subgraph(raw_full).linnet_selection == internal.linnet_selection
assert internal.filter(edge=lambda edge: edge.data.is_external).loop_count == 0
assert (
    internal.filter(edge=lambda edge: edge.data.is_external).denominator_expression()
    == 1
)
assert (internal | external).linnet_selection == full.linnet_selection
assert (internal & external).linnet_selection == empty.linnet_selection
assert (full - external).linnet_selection == internal.linnet_selection
assert (full ^ external).linnet_selection == internal.linnet_selection
assert (~internal).linnet_selection == external.linnet_selection
assert all(isinstance(item, fk.Subgraph) for item in internal.connected_components())
assert isinstance(internal.boundary(), fk.Subgraph)
assert isinstance(internal.bridges(), fk.Subgraph)
assert "<svg" in internal.to_svg()
assert "<figure" in internal._repr_html_()
assert diagram.to_json() == snapshot

# Full excision preserves all physics but starts an independent owner.
copied = full.excise()
assert type(copied) is fk.FeynmanDiagram
assert copied.to_json() == snapshot
assert copied.to_linnet() is not graph
try:
    copied.subgraph(internal)
except (ValueError, TypeError):
    pass
else:
    raise AssertionError("an excised diagram reused the parent selection owner")
try:
    copied.numerator_expression(lmb=diagram.loop_momentum_basis)
except fk.DiagramError:
    pass
else:
    raise AssertionError("an excised diagram reused the parent momentum-basis owner")

# A single interaction retains all of its particle slots as boundary legs.
region = diagram.subgraph(nodes=[0])
excised = region.excise()
assert type(excised) is fk.FeynmanDiagram
assert len(excised.vertices) == 1
assert len(excised.edges) == 3
assert all(edge.is_dangling for edge in excised.edges)
assert excised.loop_count == 0
assert excised.symmetry_factor == 1
assert excised.overall_factor_expression() == 1
assert excised.numerator_prefactor_expression() == 1
assert excised.projector_expression() == 1
assert excised.cuts == []
excised.validate()
assert (
    fk.FeynmanDiagram.from_json(model, excised.to_json()).to_json() == excised.to_json()
)
assert (
    fk.FeynmanDiagram.from_dot(model, excised.to_dot()).to_json() == excised.to_json()
)
assert diagram.to_json() == snapshot

# Partial views require explicit excision before complete-topology operations.
for operation in (region.to_json, region.to_dot):
    try:
        operation()
    except fk.DiagramError as error:
        assert "excise" in str(error)
    else:
        raise AssertionError("a partial view silently serialized its parent diagram")

# Legacy keyword filters are absent from the physics API, including subclasses.
for name in (
    "numerator_expression",
    "denominator_expression",
    "momentum_basis",
    "loop_momentum_bases",
    "compatible_momentum_basis",
    "contracted_momentum_basis",
    "superficial_degree_of_divergence",
    "uv_expansion",
    "uv_counterterm",
    "build_cff",
    "tensor_reduce",
    "reduce_tensor_numerator",
    "connected_components",
    "cycle_basis",
    "all_spanning_forests",
    "all_bonds",
    "depth_first_traverse",
    "breadth_first_traverse",
):
    assert "subgraph" not in inspect.signature(getattr(diagram, name)).parameters
    assert "subgraph" not in inspect.signature(getattr(internal, name)).parameters
assert not hasattr(diagram, "loop_count_of")

# Editing an exported graph invalidates canonical imports but not physics views.
before = internal.numerator_expression()
graph.reverse_edge(0)
assert internal.numerator_expression() == before
assert internal.original.to_json() == snapshot
try:
    diagram.subgraph(raw_full)
except (ValueError, ReferenceError):
    pass
else:
    raise AssertionError("a stale canonical selection was accepted")
assert type(internal.linnet_selection) is linnet.Subgraph

# The native Arc keeps the source available after its Python wrapper is dropped.
source = fk.FeynmanDiagram.from_json(model, snapshot)
retained = source.filter(edge=lambda edge: not edge.data.is_external)
expected = retained.numerator_expression()
del source
gc.collect()
assert retained.numerator_expression() == expected
assert retained.original.to_json() == snapshot
assert type(retained.excise()) is fk.FeynmanDiagram
print("installed physics subgraph checks passed")
