"""Check shared canonical exports and Python GC across physics views."""

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
diagram = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    loops=1,
    max_vertices=3,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=True,
).diagrams[0]
snapshot = diagram.to_json()

for first_caller in ("parent", "view"):
    parent = fk.FeynmanDiagram.from_json(model, snapshot)
    view = parent.filter(edge=lambda edge: not edge.data.is_external)
    old_graph = parent.to_linnet()
    stale = old_graph.full_subgraph()
    old_graph.reverse_edge(0)
    first = parent if first_caller == "parent" else view
    refreshed = first.to_linnet()
    assert refreshed is not old_graph
    assert refreshed is parent.to_linnet()
    assert refreshed is view.to_linnet()
    assert refreshed is view.original.to_linnet()
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
cross_section = model.generate_diagrams(
    ["scalar_0"],
    ["scalar_0", "scalar_0"],
    kind="cross_section",
    loops=1,
    max_vertices=2,
    vertex_allow=["V_3_SCALAR_000"],
    allow_self_loops=True,
).diagrams[0]
side = cross_section.cuts[0].left.subgraph
original = side.original
assert side.to_linnet() is cross_section.to_linnet() is original.to_linnet()


class BackReference:
    def __init__(self, view):
        self.view = view

    def __call__(self, edge):
        return linnet.EdgeDrawing(label=edge.data.particle_name)


def cyclic_view(kind):
    parent = fk.FeynmanDiagram.from_json(model, snapshot)
    view = parent.filter(edge=lambda edge: not edge.data.is_external)
    graph = parent.to_linnet()
    payload = BackReference(view)
    if kind == "payload":
        graph.edge(0).data = payload
    else:
        graph.render_config = graph.render_config.overlay(
            linnet.RenderConfig(selectors=linnet.DrawingSelectors(edge=payload))
        )
    return weakref.ref(payload), view.original


for kind in ("payload", "callback"):
    reference, survivor = cyclic_view(kind)
    gc.collect()
    assert reference() is not None, "a live physics owner lost its shared export"
    assert reference().view.to_linnet() is survivor.to_linnet()
    del survivor
    gc.collect()
    assert reference() is None, f"the shared export retained its {kind} cycle"

print("installed shared subgraph cache checks passed")
