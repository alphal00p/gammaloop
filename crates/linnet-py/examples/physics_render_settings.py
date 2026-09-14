# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "linnet-py==0.1.0",
#     "marimo==0.24.0",
#     "typst==0.15.0",
# ]
# ///

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

    @dataclass(frozen=True)
    class ForceSimulation:
        """Typed controls for the notebook's force-layout pass."""

        steps: int
        seed: int
        directional_force: float
        spring_strength: float
        beta: float
        dangling_repulsion: float
        dangling_centroid_repulsion: float
        edge_edge_repulsion: float
        label_steps: int

    def grouped_placement(attributes: dict[str, str]) -> dict[str, object] | None:
        """Translate GammaLoop's grouped DOT coordinates to typed drawing data."""

        raw = attributes.get("pos") or attributes.get("pin")
        if raw is None:
            return None

        placement: dict[str, object] = {"mode": lp.Placement.Pin}
        for component in raw.strip().strip('"()').split(","):
            axis, separator, coordinate = component.strip().partition(":")
            if separator == "" or axis not in {"x", "y"}:
                continue
            coordinate = coordinate.strip().removesuffix("!")
            if coordinate.startswith("@"):
                group = coordinate[1:]
            elif coordinate.startswith(("+@", "-@")):
                group = coordinate[0] + coordinate[2:]
            else:
                continue
            side = group[0] if group.startswith(("+", "-")) else None
            name = group[1:] if side is not None else group
            placement[axis] = {
                "kind": "group",
                "name": name,
                "side": side,
            }

        return placement if len(placement) > 1 else None

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
            return lp.DotEdgeData(
                edge_id=record.edge_id,
                payload=record.payload,
                statements=attributes,
                local_statements=record.local_attributes,
            )

        def decode_edge(value: lp.DotEdgeData) -> lp.EdgeValue:
            attributes = dict(value.statements)
            placement = grouped_placement(attributes)
            return lp.EdgeValue(
                data=Propagator(
                    value.edge_id,
                    attributes.get("particle", "fermion"),
                    value.payload,
                    attributes,
                    dict(value.local_statements),
                ),
                drawing=(
                    lp.EdgeDrawing(placement=placement)
                    if placement is not None
                    else None
                ),
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

    def diagram_render_settings(
        layout_algorithm: lp.LayoutAlgorithm,
        *,
        force_simulation: ForceSimulation,
        feynman_styling: bool,
        show_half_edge_ids: bool,
        show_momenta: bool,
        show_indices: bool,
    ) -> lp.RenderConfig:
        """Return either a bare or Python-defined Feynman ``RenderConfig``."""

        if layout_algorithm == lp.LayoutAlgorithm.Force:
            layouts = lp.LayoutOptions(
                algorithm=lp.LayoutAlgorithm.Force,
                direction=lp.LayoutDirection.Right,
                seed=force_simulation.seed,
                steps=force_simulation.steps,
                directional_force=force_simulation.directional_force,
                spring_strength=force_simulation.spring_strength,
                beta=force_simulation.beta,
                dangling_repulsion=force_simulation.dangling_repulsion,
                dangling_centroid_repulsion=(
                    force_simulation.dangling_centroid_repulsion
                ),
                edge_edge_repulsion=force_simulation.edge_edge_repulsion,
                label_steps=force_simulation.label_steps,
            )
        else:
            layouts = lp.LayoutOptions(
                algorithm=lp.LayoutAlgorithm.StableLayered,
                direction=lp.LayoutDirection.Right,
                label_steps=force_simulation.label_steps,
            )

        if not feynman_styling:
            return lp.RenderConfig(
                layouts=layouts,
                drawing=lp.DrawOptions(show_half_edge_ids=show_half_edge_ids),
            )

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

        return lp.RenderConfig(
            layouts=layouts,
            drawing=lp.DrawOptions(
                show_half_edge_ids=show_half_edge_ids,
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
      ext [style=invis];
      ext -> 3 [dir=none, particle="a", pin="x:@-left"];
      ext -> 2 [dir=none, particle="a", pin="x:@-left"];
      5 -> ext [dir=none, particle="a", pin="x:@+right"];
      4 -> ext [dir=none, particle="a", pin="x:@+right"];

      0 -> 1 [particle="t"];
      0 -> 1 [dir=none, particle="g"];
      5 -> 0 [particle="t"];
      1 -> 4 [particle="t"];
      3 -> 2 [particle="t"];
      2 -> 5 [particle="t"];
      4 -> 3 [particle="t"];
    }"""
    return (
        ForceSimulation,
        default_dot,
        diagram_render_settings,
        lp,
        mo,
        physics_dot_codec,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # DOT rendering with optional Feynman styling

    This notebook parses ordinary DOT into a native Linnet `Graph`. Its default
    view applies a physics-flavored rendering configuration assembled entirely
    in Python, but the Feynman-styling toggle can leave the same graph as a bare
    generic Linnet rendering instead.

    The notebook-local `DotCodec` maps particle, vertex, and port records into
    arbitrary Python dataclass instances. Only the selectors' typed drawing
    results cross into Typst. Invisible DOT vertices represent
    external legs, and edge direction retains the underlying source/sink flow
    used by the momentum arrows.

    Bare mode retains the selected layout and grouped DOT coordinates while
    omitting the particle patterns, colors, direction marks, momentum arrows,
    and physics `pᵢ` / `nᵢ` labels. It uses generic `eᵢ` / `nᵢ` structural IDs
    instead. The optional `hᵢ` half-edge IDs work in either mode; momentum and
    physics-index controls apply only to Feynman mode.

    The force-simulation panel updates after a slider is released. Its label
    relaxation setting applies after either layout; the other settings apply
    only when the Force layout is selected.

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
        label="Editable DOT",
    )
    layout_algorithm = mo.ui.dropdown(
        options={
            "Force": lp.LayoutAlgorithm.Force,
            "Stable layered": lp.LayoutAlgorithm.StableLayered,
        },
        value="Force",
        label="Layout",
    )
    feynman_styling = mo.ui.checkbox(
        value=True,
        label="Feynman diagram styling",
    )
    show_momenta = mo.ui.checkbox(value=True, label="Momentum arrows")
    show_indices = mo.ui.checkbox(value=True, label="pᵢ / nᵢ labels")
    show_half_edge_ids = mo.ui.checkbox(value=False, label="hᵢ half-edge IDs")
    mo.vstack(
        [
            mo.hstack(
                [
                    feynman_styling,
                    layout_algorithm,
                    show_momenta,
                    show_indices,
                    show_half_edge_ids,
                ],
                justify="start",
                wrap=True,
                gap=1.5,
            ),
            dot_source,
        ]
    )
    return (
        dot_source,
        feynman_styling,
        layout_algorithm,
        show_half_edge_ids,
        show_indices,
        show_momenta,
    )


@app.cell
def _(mo):
    force_steps = mo.ui.slider(
        40,
        640,
        20,
        320,
        debounce=True,
        show_value=True,
        label="Force iterations",
    )
    force_seed = mo.ui.slider(
        0,
        99,
        1,
        19,
        debounce=True,
        show_value=True,
        label="Seed",
    )
    directional_force = mo.ui.slider(
        0,
        5,
        0.05,
        0.45,
        debounce=True,
        show_value=True,
        label="Directional force",
    )
    spring_strength = mo.ui.slider(
        1,
        100,
        1,
        11,
        debounce=True,
        show_value=True,
        label="Spring strength",
    )
    beta = mo.ui.slider(
        0,
        250,
        5,
        50,
        debounce=True,
        show_value=True,
        label="Vertex repulsion (β)",
    )
    dangling_repulsion = mo.ui.slider(
        0,
        10,
        0.1,
        5,
        debounce=True,
        show_value=True,
        label="External-leg repulsion",
    )
    dangling_centroid_repulsion = mo.ui.slider(
        0,
        5,
        0.05,
        1.25,
        debounce=True,
        show_value=True,
        label="External legs from node centroid",
    )
    edge_edge_repulsion = mo.ui.slider(
        0,
        0.5,
        0.01,
        0.1,
        debounce=True,
        show_value=True,
        label="Edge–edge repulsion",
    )
    label_steps = mo.ui.slider(
        0,
        200,
        10,
        80,
        debounce=True,
        show_value=True,
        label="Label relaxation steps",
    )
    mo.vstack(
        [
            mo.md(
                "### Force simulation\n\n"
                "Changes render after releasing a slider. The first eight controls "
                "apply only to Force; label relaxation runs after either layout."
            ),
            force_steps,
            force_seed,
            directional_force,
            spring_strength,
            beta,
            dangling_repulsion,
            dangling_centroid_repulsion,
            edge_edge_repulsion,
            label_steps,
        ],
        gap=0.75,
    )
    return (
        beta,
        dangling_repulsion,
        dangling_centroid_repulsion,
        directional_force,
        edge_edge_repulsion,
        force_seed,
        force_steps,
        label_steps,
        spring_strength,
    )


@app.cell
def _(
    ForceSimulation,
    beta,
    dangling_repulsion,
    dangling_centroid_repulsion,
    directional_force,
    edge_edge_repulsion,
    force_seed,
    force_steps,
    label_steps,
    spring_strength,
):
    force_simulation = ForceSimulation(
        steps=force_steps.value,
        seed=force_seed.value,
        directional_force=directional_force.value,
        spring_strength=spring_strength.value,
        beta=beta.value,
        dangling_repulsion=dangling_repulsion.value,
        dangling_centroid_repulsion=dangling_centroid_repulsion.value,
        edge_edge_repulsion=edge_edge_repulsion.value,
        label_steps=label_steps.value,
    )
    return force_simulation


@app.cell
def _(
    diagram_render_settings,
    dot_source,
    feynman_styling,
    force_simulation,
    layout_algorithm,
    lp,
    physics_dot_codec,
    show_half_edge_ids,
    show_indices,
    show_momenta,
):
    try:
        graph = lp.Graph.from_dot(dot_source.value, physics_dot_codec())
        graph.render_config = diagram_render_settings(
            layout_algorithm.value,
            force_simulation=force_simulation,
            feynman_styling=feynman_styling.value,
            show_half_edge_ids=show_half_edge_ids.value,
            show_momenta=show_momenta.value,
            show_indices=show_indices.value,
        )
        graph.render_config.title = graph.name or "Parsed DOT graph"
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
    return render_error, rendered_svg, typst_source


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
        _output = mo.Html(
            f'<div data-linnet-render-ready="physics">{rendered_svg}</div>'
        )
    _output
    return


@app.cell
def _(mo, typst_source):
    if typst_source is None:
        _source_panel = mo.md(
            "Fix the DOT or rendering error to inspect the generated Typst."
        )
    else:
        _source_panel = mo.vstack(
            [
                mo.md(r"""
                ## Live generated Typst

                This is the exact staged entrypoint used for the SVG above.
                It updates with the DOT and rendering controls.
                """),
                mo.ui.code_editor(
                    value=typst_source,
                    language="text",
                    disabled=True,
                    min_height=320,
                    max_height=700,
                    label="Generated Typst (read-only)",
                ),
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
                ),
            ]
        )
    _details
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
