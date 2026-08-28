"""A Python-only physics theme built from Linnet's generic rendering API.

This file is an example, not part of ``linnet_py``.  The graph payloads,
particle vocabulary, colors, line decorations, labels, and momentum arrows all
belong to this application-level render configuration.
"""

from dataclasses import dataclass

import linnet_py as lp


@dataclass(eq=False)
class Propagator:
    """Arbitrary application data inspected by the drawing selectors."""

    particle: str


def physics_render_settings() -> lp.RenderConfig:
    """Return one custom physics-flavored ``RenderConfig``."""

    black = lp.Color("black")
    blue = lp.Color("blue")
    orange = lp.Color("orange")
    momentum_stroke = lp.Stroke(
        paint=lp.Color("red"),
        thickness=lp.Length.pt(0.45),
    )

    def node_drawing(node: lp.Node) -> lp.NodeDrawing:
        return lp.NodeDrawing(label=lp.MathSymbol("n", subscript=node.index))

    def edge_drawing(edge: lp.Edge) -> lp.EdgeDrawing:
        particle = edge.data.particle
        paint = {"photon": blue, "scalar": orange}.get(particle, black)
        stroke = lp.Stroke(
            paint=paint,
            thickness=lp.Length.pt(0.65),
            dash=(
                lp.Dash(lp.DashPattern.Dashed) if particle == "scalar" else lp.INHERIT
            ),
        )
        return lp.EdgeDrawing(
            label=lp.MathSymbol("p", subscript=edge.index),
            label_style={"fill": paint},
            style={"stroke": stroke},
            decoration=lp.Pattern.Wave if particle == "photon" else lp.INHERIT,
        )

    def half_edge_drawing(half_edge: lp.HalfEdge) -> lp.HalfEdgeDrawing:
        particle_layer = {}
        if half_edge.edge.data.particle == "fermion":
            particle_layer = {
                "mark": lp.Mark.barbed(),
                "mark-position": lp.MarkPosition.CenterIfDangling,
                "mark-orientation": lp.MarkOrientation.Edge,
            }

        edge = half_edge.edge
        arrow_half = edge.sink if edge.sink is not None else edge.source
        momentum_layer = {
            "offset": 0.46,
            "length": 5.0,
            "ratio": 0.5,
            "resolve-length": lp.EdgeLengthResolution.Min,
            "offset-side": "label",
            "stroke": momentum_stroke,
        }
        if arrow_half is not None and half_edge.index == arrow_half.index:
            momentum_layer["mark"] = lp.Mark(
                end=lp.MarkSymbol.Straight,
                stroke=momentum_stroke,
                scale=0.75,
            )
        return lp.HalfEdgeDrawing(style=(particle_layer, momentum_layer))

    return lp.RenderConfig(
        layouts=lp.LayoutOptions(
            algorithm=lp.LayoutAlgorithm.StableLayered,
            direction=lp.LayoutDirection.Right,
            label_steps=60,
        ),
        drawing=lp.DrawOptions(
            node_fill=lp.Color("white"),
            node_stroke=lp.Stroke(paint=black, thickness=lp.Length.pt(0.6)),
        ),
        selectors=lp.DrawingSelectors(
            node=node_drawing,
            edge=edge_drawing,
            source=half_edge_drawing,
            sink=half_edge_drawing,
        ),
    )


def example_graph() -> lp.Graph:
    incoming = lp.node("incoming")
    interaction = lp.node("interaction")
    outgoing = lp.node("outgoing")
    graph = lp.build(
        incoming,
        interaction,
        outgoing,
        lp.edge(lp.sink(incoming), "incoming fermion", data=Propagator("fermion")),
        lp.edge(
            lp.source(incoming),
            "fermion",
            lp.sink(interaction),
            data=Propagator("fermion"),
        ),
        lp.edge(
            lp.source(interaction),
            "photon",
            lp.sink(outgoing),
            data=Propagator("photon"),
        ),
        lp.edge(lp.source(outgoing), "outgoing scalar", data=Propagator("scalar")),
        name="python-only-physics-theme",
    )
    graph.render_config = physics_render_settings()
    return graph


if __name__ == "__main__":
    example_graph().render("physics-render-settings.svg")
