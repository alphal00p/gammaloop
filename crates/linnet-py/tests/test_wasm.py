"""Runtime smoke test shared by native and Pyodide wheel jobs."""

import unittest
from dataclasses import dataclass

import linnet_py as lp

DOT = r"""digraph browser {
  0;
  1;

  exte0 [style=invis];
  exte0 -> 0:0 [id=0, particle="a", pos="x:@-left!"];
  0:1 -> 1:2 [id=1, particle="g"];
  exte1 [style=invis];
  1:3 -> exte1 [id=2, particle="a", pos="x:@+right!"];
}"""


@dataclass(eq=False)
class NodeData:
    name: str | None


@dataclass(eq=False)
class EdgeData:
    attributes: dict[str, str]


def grouped_placement(attributes: dict[str, str]) -> dict[str, object] | None:
    raw = attributes.get("pos") or attributes.get("pin")
    if raw is None:
        return None

    placement: dict[str, object] = {"mode": lp.Placement.Pin}
    for component in raw.strip().strip('"()').split(","):
        axis, separator, coordinate = component.strip().partition(":")
        if separator == "" or axis not in {"x", "y"}:
            continue
        group = coordinate.strip().removesuffix("!")
        if group.startswith("@"):
            group = group[1:]
        elif group.startswith(("+@", "-@")):
            group = group[0] + group[2:]
        else:
            continue
        side = group[0] if group.startswith(("+", "-")) else None
        placement[axis] = {
            "kind": "group",
            "name": group[1:] if side is not None else group,
            "side": side,
        }
    return placement if len(placement) > 1 else None


def codec() -> lp.DotCodec:
    def encode_node(value: lp.NodeValue) -> lp.DotVertexData:
        return lp.DotVertexData(name=value.data.name)

    def decode_node(value: lp.DotVertexData) -> lp.NodeValue:
        return lp.NodeValue(data=NodeData(value.name))

    def encode_edge(value: lp.EdgeValue) -> lp.DotEdgeData:
        return lp.DotEdgeData(statements=value.data.attributes)

    def decode_edge(value: lp.DotEdgeData) -> lp.EdgeValue:
        attributes = dict(value.statements)
        placement = grouped_placement(attributes)
        return lp.EdgeValue(
            data=EdgeData(attributes),
            drawing=(
                lp.EdgeDrawing(placement=placement) if placement is not None else None
            ),
        )

    def encode_half_edge(_value: lp.HalfEdgeValue) -> lp.DotHalfEdgeData:
        return lp.DotHalfEdgeData()

    def decode_half_edge(_value: lp.DotHalfEdgeData) -> lp.HalfEdgeValue:
        return lp.HalfEdgeValue()

    return lp.DotCodec(
        encode_node=encode_node,
        decode_node=decode_node,
        encode_edge=encode_edge,
        decode_edge=decode_edge,
        encode_half_edge=encode_half_edge,
        decode_half_edge=decode_half_edge,
    )


def wasm_runtime_smoke() -> None:
    graph = lp.Graph.from_dot(DOT, codec())
    assert (graph.n_nodes, graph.n_edges, graph.n_half_edges) == (2, 3, 4)

    incoming = graph.edge(0).drawing.placement
    outgoing = graph.edge(2).drawing.placement
    assert incoming["x"] == {"kind": "group", "name": "left", "side": "-"}
    assert outgoing["x"] == {"kind": "group", "name": "right", "side": "+"}

    selected_nodes: list[int] = []
    selected_edges: list[int] = []
    selected_half_edges: list[int] = []

    def node_drawing(node: lp.Node) -> lp.NodeDrawing:
        selected_nodes.append(node.index)
        return lp.NodeDrawing(label=lp.MathSymbol("n", subscript=node.index))

    def edge_drawing(edge: lp.Edge) -> lp.EdgeDrawing:
        selected_edges.append(edge.index)
        return lp.EdgeDrawing(label=lp.MathSymbol("p", subscript=edge.index))

    def half_edge_drawing(half_edge: lp.HalfEdge) -> lp.HalfEdgeDrawing:
        selected_half_edges.append(half_edge.index)
        return lp.HalfEdgeDrawing(
            style={
                "stroke": lp.Stroke(
                    paint=lp.Color("black"),
                    thickness=lp.Length.pt(0.5),
                )
            }
        )

    graph.render_config = lp.RenderConfig(
        layouts=lp.LayoutOptions(
            algorithm=lp.LayoutAlgorithm.Force,
            direction=lp.LayoutDirection.Right,
            seed=7,
            steps=40,
            directional_force=0.4,
            label_steps=8,
        ),
        selectors=lp.DrawingSelectors(
            node=node_drawing,
            edge=edge_drawing,
            source=half_edge_drawing,
            sink=half_edge_drawing,
        ),
        drawing=lp.DrawOptions(show_half_edge_ids=True),
    )

    prepared = graph.prepare_render()
    assert selected_nodes == list(range(graph.n_nodes))
    assert selected_edges == list(range(graph.n_edges))
    assert selected_half_edges == list(range(graph.n_half_edges))
    assert "_linnet_template.render" in prepared.typst_source
    assert "graph-spec-path" in prepared.typst_source

    svg = prepared.to_svg()
    assert svg.lstrip().startswith("<svg")
    assert svg.rstrip().endswith("</svg>")


class WasmRuntimeTests(unittest.TestCase):
    def test_dot_selectors_and_typst_render(self) -> None:
        wasm_runtime_smoke()
