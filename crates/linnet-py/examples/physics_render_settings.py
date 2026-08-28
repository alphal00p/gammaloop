# ruff: noqa: B018, PLR1711  # Cell outputs and empty returns are Marimo syntax.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    from dataclasses import dataclass

    import linnet_py as lp
    import marimo as mo

    @dataclass(eq=False)
    class DotInteraction:
        """Arbitrary application data reconstructed from one DOT vertex."""

        name: str | None
        index: int | None
        payload: list[int] | None
        attributes: dict[str, str]

    @dataclass(eq=False)
    class Propagator:
        """Arbitrary application data inspected by the drawing selectors."""

        edge_id: int | None
        particle: str
        momentum: str | None
        payload: list[int] | None
        attributes: dict[str, str]
        local_attributes: dict[str, str]

    @dataclass(eq=False)
    class DotPort:
        """Arbitrary application data reconstructed from one DOT endpoint."""

        index: int | None
        statement: str | None
        payload: list[int] | None
        port_label: str | None
        compass: str | None

    def physics_dot_codec() -> lp.DotCodec:
        """Map ordinary physics DOT attributes to arbitrary Python payloads."""

        def encode_node(value: lp.NodeValue) -> lp.DotVertexData:
            record = value.data
            return lp.DotVertexData(
                name=record.name,
                index=record.index,
                payload=record.payload,
                statements=record.attributes,
            )

        def decode_node(value: lp.DotVertexData) -> lp.NodeValue:
            return lp.NodeValue(
                data=DotInteraction(
                    value.name,
                    value.index,
                    value.payload,
                    dict(value.statements),
                )
            )

        def encode_edge(value: lp.EdgeValue) -> lp.DotEdgeData:
            record = value.data
            attributes = dict(record.attributes)
            attributes["particle"] = record.particle
            if record.momentum is None:
                attributes.pop("lmb_rep", None)
            else:
                attributes["lmb_rep"] = record.momentum
            return lp.DotEdgeData(
                edge_id=record.edge_id,
                payload=record.payload,
                statements=attributes,
                local_statements=record.local_attributes,
            )

        def decode_edge(value: lp.DotEdgeData) -> lp.EdgeValue:
            attributes = dict(value.statements)
            return lp.EdgeValue(
                data=Propagator(
                    value.edge_id,
                    attributes.get("particle", "fermion"),
                    attributes.get("lmb_rep"),
                    value.payload,
                    attributes,
                    dict(value.local_statements),
                )
            )

        def encode_half_edge(value: lp.HalfEdgeValue) -> lp.DotHalfEdgeData:
            record = value.data
            return lp.DotHalfEdgeData(
                statement=record.statement,
                index=record.index,
                payload=record.payload,
                port_label=record.port_label,
                compass=record.compass,
            )

        def decode_half_edge(value: lp.DotHalfEdgeData) -> lp.HalfEdgeValue:
            return lp.HalfEdgeValue(
                data=DotPort(
                    value.index,
                    value.statement,
                    value.payload,
                    value.port_label,
                    value.compass,
                )
            )

        return lp.DotCodec(
            encode_node=encode_node,
            decode_node=decode_node,
            encode_edge=encode_edge,
            decode_edge=decode_edge,
            encode_half_edge=encode_half_edge,
            decode_half_edge=decode_half_edge,
        )

    def physics_render_settings(
        layout_algorithm: lp.LayoutAlgorithm,
        *,
        show_momenta: bool,
        show_indices: bool,
    ) -> lp.RenderConfig:
        """Return one custom physics-flavored ``RenderConfig``."""

        black = lp.Color("black")
        blue = lp.Color("blue")
        # GammaLoop lightens sink halves by 45%; Color stores the resulting value.
        light_black = lp.Color.rgb(115, 115, 115)
        light_blue = lp.Color.rgb(115, 179, 234)
        particle_kinds = {
            "a": "photon",
            "g": "gluon",
            "H": "scalar",
            "t": "fermion",
        }
        source_paints = {"fermion": blue}
        sink_paints = {"fermion": light_blue}
        particle_patterns = {
            "photon": {
                "pattern": lp.Pattern.Wave,
                "pattern-amplitude": 0.14,
                "pattern-wavelength": 0.55,
            },
            "gluon": {
                "pattern": lp.Pattern.Coil,
                "pattern-amplitude": 0.14,
                "pattern-wavelength": 0.55,
                "pattern-coil-longitudinal-scale": 1.6,
            },
        }
        scalar_dash = lp.Dash.pattern((lp.Length.em(0.1), lp.Length.em(0.45)))
        fermion_mark = lp.Mark(
            end=lp.MarkSymbol.Barbed,
            fill=black,
            stroke=lp.Stroke(paint=black, thickness=lp.Length.pt(0.2)),
            scale=0.75,
            anchor=lp.Anchor.Center,
            shorten_to=lp.AUTO,
        )
        momentum_stroke = lp.Stroke(
            paint=black,
            thickness=lp.Length.pt(0.55),
            cap=lp.StrokeCap.Round,
        )

        def node_drawing(node: lp.Node) -> lp.NodeDrawing:
            label = lp.MathSymbol("n", subscript=node.index) if show_indices else None
            return lp.NodeDrawing(label=label)

        def edge_drawing(edge: lp.Edge) -> lp.EdgeDrawing:
            label = lp.MathSymbol("p", subscript=edge.index) if show_indices else None
            return lp.EdgeDrawing(
                label=label,
                label_style={"fill": black},
            )

        def half_edge_drawing(
            half_edge: lp.HalfEdge,
        ) -> lp.HalfEdgeDrawing:
            edge = half_edge.edge
            particle = edge.data.particle
            kind = particle_kinds.get(particle, particle)
            is_sink = edge.sink is not None and half_edge.index == edge.sink.index
            paints = sink_paints if is_sink else source_paints
            paint = paints.get(
                kind,
                light_black if is_sink else black,
            )
            thickness = lp.Length.pt(1.0 if kind in {"fermion", "scalar"} else 0.55)
            stroke_options = {
                "paint": paint,
                "thickness": thickness,
                "cap": lp.StrokeCap.Round,
            }
            if kind == "scalar":
                stroke_options["dash"] = scalar_dash
            particle_layer = {
                "stroke": lp.Stroke(**stroke_options),
                **particle_patterns.get(kind, {}),
            }
            if kind == "fermion":
                particle_layer.update(
                    {
                        "mark": fermion_mark,
                        "mark-position": lp.MarkPosition.CenterIfDangling,
                        "mark-orientation": lp.MarkOrientation.Edge,
                    }
                )

            layers = [particle_layer]

            if show_momenta:
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
                layers.append(momentum_layer)

            return lp.HalfEdgeDrawing(style=tuple(layers))

        if layout_algorithm == lp.LayoutAlgorithm.Force:
            layouts = lp.LayoutOptions(
                algorithm=lp.LayoutAlgorithm.Force,
                direction=lp.LayoutDirection.Right,
                seed=19,
                steps=320,
                directional_force=0.45,
                label_steps=80,
            )
        else:
            layouts = lp.LayoutOptions(
                algorithm=lp.LayoutAlgorithm.StableLayered,
                direction=lp.LayoutDirection.Right,
                label_steps=80,
            )

        return lp.RenderConfig(
            layouts=layouts,
            drawing=lp.DrawOptions(
                node_fill=lp.Color("white"),
                node_stroke=lp.Stroke(
                    paint=black,
                    thickness=lp.Length.pt(0.6),
                ),
            ),
            selectors=lp.DrawingSelectors(
                node=node_drawing,
                edge=edge_drawing,
                source=half_edge_drawing,
                sink=half_edge_drawing,
            ),
        )

    default_dot = r"""digraph GL05 {
      num = "1";
      overall_factor = "(AutG(1))^(-1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(4)";
      overall_factor_evaluated = "-4";
      projector = "epsilon(0)*epsilon(1)*epsilon-bar(2)*epsilon-bar(3)";

      0 [int_id="V_137"];
      1 [int_id="V_137"];
      2 [int_id="V_134"];
      3 [int_id="V_134"];
      4 [int_id="V_134"];
      5 [int_id="V_134"];

      exte0 [style=invis];
      exte0 -> 3:0 [id=0, dir=none, lmb_rep="P(0,a___)", particle="a", pin="x:@-left"];
      exte1 [style=invis];
      exte1 -> 2:1 [id=1, dir=none, lmb_rep="P(1,a___)", particle="a", pin="x:@-left"];
      exte2 [style=invis];
      5:2 -> exte2 [id=2, dir=none, lmb_rep="P(2,a___)", particle="a", pin="x:@+right"];
      exte3 [style=invis];
      4:3 -> exte3 [id=3, dir=none, lmb_rep="-P(2,a___)+P(0,a___)+P(1,a___)", particle="a", pin="x:@+right"];

      0:4 -> 1:5 [id=4, lmb_id="0", lmb_rep="K(0,a___)", particle="t"];
      0:6 -> 1:7 [id=5, dir=none, lmb_rep="-K(0,a___)+K(1,a___)", particle="g"];
      5:8 -> 0:9 [id=6, lmb_rep="K(1,a___)", particle="t"];
      1:10 -> 4:11 [id=7, lmb_id="1", lmb_rep="K(1,a___)", particle="t"];
      3:12 -> 2:13 [id=8, lmb_rep="-P(1,a___)+K(1,a___)+P(2,a___)", particle="t"];
      2:14 -> 5:15 [id=9, lmb_rep="K(1,a___)+P(2,a___)", particle="t"];
      4:16 -> 3:17 [id=10, lmb_rep="-P(0,a___)-P(1,a___)+K(1,a___)+P(2,a___)", particle="t"];
    }"""
    return default_dot, lp, mo, physics_dot_codec, physics_render_settings


@app.cell
def _(mo):
    mo.md(r"""
    # Physics rendering from editable DOT

    This notebook parses ordinary DOT into a native Linnet `Graph`, then applies
    a physics-flavored rendering configuration assembled entirely in Python.
    Its particle geometry and independent momentum-arrow layer mirror
    GammaLoop's current Linnest template, but remain an example rather than a
    built-in Linnet mode.

    The notebook-local `DotCodec` maps `particle`, momentum, vertex, and port
    records into arbitrary Python dataclass instances. Only the selectors'
    typed drawing results cross into Typst. Invisible DOT vertices represent
    external legs, and edge direction retains the underlying source/sink flow
    used by the momentum arrows.

    The read-only panel below the diagram shows the exact generated Typst
    entrypoint compiled for the current view.
    """)
    return


@app.cell
def _(default_dot, lp, mo):
    dot_source = mo.ui.code_editor(
        value=default_dot,
        language="text",
        min_height=440,
        max_height=700,
        debounce=400,
        label="Editable physics DOT",
    )
    layout_algorithm = mo.ui.dropdown(
        options={
            "Force": lp.LayoutAlgorithm.Force,
            "Stable layered": lp.LayoutAlgorithm.StableLayered,
        },
        value="Force",
        label="Layout",
    )
    show_momenta = mo.ui.checkbox(value=True, label="Momentum arrows")
    show_indices = mo.ui.checkbox(value=True, label="pᵢ / nᵢ labels")
    mo.vstack(
        [
            mo.hstack(
                [layout_algorithm, show_momenta, show_indices],
                justify="start",
                gap=1.5,
            ),
            dot_source,
        ]
    )
    return dot_source, layout_algorithm, show_indices, show_momenta


@app.cell
def _(
    dot_source,
    layout_algorithm,
    lp,
    physics_dot_codec,
    physics_render_settings,
    show_indices,
    show_momenta,
):
    try:
        graph = lp.Graph.from_dot(dot_source.value, physics_dot_codec())
        graph.render_config = physics_render_settings(
            layout_algorithm.value,
            show_momenta=show_momenta.value,
            show_indices=show_indices.value,
        )
        graph.render_config.title = graph.name or "Parsed physics graph"
        parse_error = None
    except (RuntimeError, TypeError, ValueError) as error:
        graph = None
        parse_error = f"{type(error).__name__}: {error}"
    return graph, parse_error


@app.cell
def _(graph):
    if graph is None:
        prepared_render = None
        render_error = None
        rendered_svg = None
        typst_source = None
    else:
        try:
            prepared_render = graph.prepare_render()
            typst_source = prepared_render.typst_source
            rendered_svg = prepared_render.to_svg()
            render_error = None
        except (OSError, RuntimeError, TypeError, ValueError) as error:
            prepared_render = None
            render_error = f"{type(error).__name__}: {error}"
            rendered_svg = None
            typst_source = None
    return prepared_render, render_error, rendered_svg, typst_source


@app.cell
def _(mo, parse_error, render_error, rendered_svg):
    if parse_error is not None:
        _output = mo.callout(
            mo.md(f"`{parse_error}`"),
            kind="danger",
            title="DOT parsing failed",
        )
    elif render_error is not None:
        _output = mo.callout(
            mo.md(f"`{render_error}`"),
            kind="danger",
            title="Typst rendering failed",
        )
    else:
        _output = mo.Html(rendered_svg)
    _output
    return


@app.cell
def _(mo, parse_error, render_error, typst_source):
    if typst_source is None:
        _source_panel = mo.md(
            "Fix the DOT or rendering error to inspect the generated Typst."
        )
    else:
        _source_viewer = mo.ui.code_editor(
            value=typst_source,
            language="text",
            disabled=True,
            min_height=260,
            max_height=560,
            show_copy_button=True,
            label="Generated Typst entrypoint",
        )
        _source_panel = mo.vstack(
            [
                mo.md(r"""
                ## Live generated Typst

                This is the exact staged entrypoint used for the SVG above.
                It updates with the DOT and rendering controls.
                """),
                _source_viewer,
            ]
        )
    _source_panel
    return


@app.cell
def _(graph, mo):
    if graph is None:
        _details = mo.md("Fix the DOT input to inspect its parsed records.")
    else:
        _rows = []
        for _edge in graph.edges():
            _source = (
                (_edge.source.node.name or f"n{_edge.source.node.index}")
                if _edge.source is not None
                else "external"
            )
            _sink = (
                (_edge.sink.node.name or f"n{_edge.sink.node.index}")
                if _edge.sink is not None
                else "external"
            )
            _rows.append(
                {
                    "pᵢ": _edge.index,
                    "source": _source,
                    "sink": _sink,
                    "particle": _edge.data.particle,
                    "momentum": _edge.data.momentum or "",
                    "orientation": str(_edge.orientation)
                    .removeprefix("Orientation.")
                    .lower(),
                }
            )
        _details = mo.vstack(
            [
                mo.md(f"""
                ## Parsed native graph

                **{graph.n_nodes} nodes**, **{graph.n_edges} edges**,
                **{graph.n_half_edges} half-edges**, and
                **{graph.external_half_edges().n_half_edges} external legs**.
                The table reads the arbitrary `Propagator` objects reconstructed
                by the codec; these objects are never serialized to Typst.
                """),
                mo.ui.table(
                    _rows,
                    pagination=False,
                    selection=None,
                    show_download=False,
                    wrapped_columns=["momentum"],
                ),
            ]
        )
    _details
    return


if __name__ == "__main__":
    app.run()
