"""Check shared canonical exports and Python GC across physics views."""

import gc
import weakref
from pathlib import Path

from symbolica.community import graph as linnet

from symbolica.community import hepkit as fk

model = fk.Model(Path(__file__).parent / "fixtures/scalars_2p_3p.json")
diagram = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_diagrams(loops=1, max_vertices=3, allow_self_loops=True)
    .diagrams[0]
)
snapshot = diagram.to_json()

for first_caller in ("parent", "view"):
    parent = fk.FeynmanDiagram.from_json(model, snapshot)
    view = parent.filter(edge=lambda edge: not edge.is_external)
    old_graph = parent.to_graph()
    stale = old_graph.full_subgraph()
    old_graph.reverse_edge(0)
    first = parent if first_caller == "parent" else view
    refreshed = first.to_graph()
    assert refreshed is not old_graph
    assert refreshed is parent.to_graph()
    assert refreshed is view.to_graph()
    assert refreshed is view.original.to_graph()
    imported = parent.subgraph(view.linnet_selection)
    assert imported.linnet_selection == view.linnet_selection
    assert imported.numerator_expression() == view.numerator_expression()
    try:
        parent.subgraph(stale)
    except (ValueError, ReferenceError):
        pass
    else:
        raise AssertionError("a stale canonical selection was accepted after refresh")

# Native cut factories can create views before any graph has been exported.
cross_section = (
    model.process(
        ["scalar_0"], ["scalar_0", "scalar_0"], vertex_allow=["V_3_SCALAR_000"]
    )
    .generate_cross_section(loops=1, max_vertices=2, allow_self_loops=True)
    .diagrams[0]
)
side = cross_section.cuts[0].left.subgraph
original = side.original
assert side.to_graph() is cross_section.to_graph() is original.to_graph()


class BackReference:
    def __init__(self, view):
        self.view = view

    def __call__(self, edge):
        return linnet.EdgeDrawing(label=edge.data.particle_name)


def cyclic_view(kind):
    parent = fk.FeynmanDiagram.from_json(model, snapshot)
    view = parent.filter(edge=lambda edge: not edge.is_external)
    graph = parent.to_graph()
    payload = BackReference(view)
    if kind == "payload":
        graph.edge(0).data = payload
    else:
        graph.render_config = graph.render_config.overlay(
            linnet.RenderSettings(selectors=linnet.DrawingSelectors(edge=payload))
        )
    return weakref.ref(payload), view.original


for kind in ("payload", "callback"):
    reference, survivor = cyclic_view(kind)
    gc.collect()
    assert reference() is not None, "a live physics owner lost its shared export"
    assert reference().view.to_graph() is survivor.to_graph()
    del survivor
    gc.collect()
    assert reference() is None, f"the shared export retained its {kind} cycle"

# Native result wrappers can also retain the optional export through their
# immutable owner; GC must see each counted holder reference.
for result_kind in ("tree", "partition"):
    parent = fk.FeynmanDiagram.from_json(model, cross_section.to_json())
    result = (
        parent.depth_first_traverse(0)
        if result_kind == "tree"
        else parent.all_cuts([0], [1])[0]
    )
    graph = parent.to_graph()
    payload = BackReference(result)
    graph.edge(0).data = payload
    reference = weakref.ref(payload)
    del parent, result, graph, payload
    gc.collect()
    assert reference() is None, f"the native {result_kind} retained its export cycle"

print("installed shared subgraph cache checks passed")
