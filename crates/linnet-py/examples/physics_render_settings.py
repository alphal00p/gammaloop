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
    import base64
    import io
    import json
    import os
    import tempfile
    import zipfile
    from dataclasses import dataclass
    from pathlib import Path

    import linnet_py as lp
    import marimo as mo
    import typst

    @dataclass(eq=False)
    class DotInteraction:
        """Arbitrary application data reconstructed from one DOT vertex."""

        name: str | None
        index: int | None
        payload: list[int] | None
        attributes: dict[str, str]

    @dataclass(eq=False)
    class Propagator:
        """Particle data inspected in the parsed-edge table."""

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
            return lp.EdgeValue(
                data=Propagator(
                    value.edge_id,
                    attributes.get("particle", "fermion"),
                    value.payload,
                    attributes,
                    dict(value.local_statements),
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
        layout_algorithm: str,
        *,
        force_simulation: ForceSimulation,
        custom_forces: bool,
        mode: str,
        show_half_edge_ids: bool,
        show_momenta: bool,
        show_momentum_labels: bool,
    ) -> dict[str, object]:
        """Supply native options to GammaLoop's existing Typst figure template."""

        layout = {
            "steps": force_simulation.steps,
            "seed": force_simulation.seed,
            "layout-algo": layout_algorithm,
        }
        if custom_forces:
            layout.update(
                {
                    "directional-force": force_simulation.directional_force,
                    "k-spring": force_simulation.spring_strength,
                    "beta": force_simulation.beta,
                    "gamma-dangling": force_simulation.dangling_repulsion,
                    "gamma-dangling-centroid": force_simulation.dangling_centroid_repulsion,
                    "gamma-ee": force_simulation.edge_edge_repulsion,
                    "label-steps": force_simulation.label_steps,
                }
            )
        # GammaLoop lightens sink halves by 45%; the shared Typst callbacks own
        # this and all other particle, label, arrow, and placement conventions.
        return {
            "layouts": [layout],
            "options": {
                "mode": mode,
                "momentum-arrows": show_momenta,
                "show-momentum": show_momentum_labels,
                "debug": show_half_edge_ids,
            },
        }

    default_dot = r"""digraph GL05 {
      ext [style=invis];
      ext -> 3 [dir=none, particle="a"];
      ext -> 2 [dir=none, particle="a"];
      // Order the outgoing legs along the outer cycle to avoid a crossing.
      5:12 -> ext [dir=none, particle="a"];
      4:15 -> ext [dir=none, particle="a"];

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
        Path,
        base64,
        default_dot,
        diagram_render_settings,
        io,
        json,
        lp,
        mo,
        os,
        physics_dot_codec,
        tempfile,
        typst,
        zipfile,
    )


@app.cell
def _(Path, base64, io, os, tempfile, zipfile):
    # The exporter embeds a deterministic archive from GammaLoop's save-dot
    # command. Native runs can reuse its public ZIP through this environment path.
    drawing_bundle = None
    if drawing_bundle is None:
        _bundle_path = os.environ.get("GAMMALOOP_DRAWING_BUNDLE")
        if _bundle_path is None:
            raise RuntimeError(
                "Set GAMMALOOP_DRAWING_BUNDLE to the exported "
                "public/gammaloop-drawing.zip when running this notebook locally."
            )
        _archive = Path(_bundle_path).read_bytes()
    else:
        _archive = base64.b64decode(drawing_bundle)
    drawing_workspace = tempfile.TemporaryDirectory(prefix="gammaloop-notebook-")
    drawing_root = Path(drawing_workspace.name)
    with zipfile.ZipFile(io.BytesIO(_archive)) as _zip:
        _zip.extractall(drawing_root)
    return drawing_root, drawing_workspace


@app.cell
def _(mo):
    mo.md(r"""
    # GammaLoop DOT drawing

    This notebook uses the same Typst templates and Standard Model particle map
    as `save dot` followed by `just draw`. Python supplies the edited DOT and
    control values; GammaLoop styles the particles and orders external legs,
    and Linnest measures, lays out, and draws the graph.

    **Automatic** placement recognizes amplitude legs and cross-section cut
    pairs. Existing X/Y positions are preserved, and external depth is pinned
    to zero. Particle labels face outward; optional momentum labels use `qₑ`
    with the original edge ID. Debug mode adds node and half-edge IDs.

    Open **Layout settings** for the sliders. GammaLoop's mode-specific presets
    stay active unless you enable custom force parameters. Edits render after
    typing pauses; slider changes render on release.
    """)
    return


@app.cell
def _(default_dot, mo):
    example = mo.ui.dropdown(
        options={
            "Amplitude": default_dot,
            "Cross-section": r"""digraph Cut {
  ext [style=invis];
  ext -> a [particle="e-", is_cut=0];
  ext -> b [particle="e+", is_cut=1];
  c -> ext [particle="e+", is_cut=1];
  d -> ext [particle="e-", is_cut=0];
  a -> b [particle="t"];
  b -> c [particle="t"];
  c -> d [particle="t"];
  d -> a [particle="t"];
  a -> c [particle="g", dir=none];
}""",
        },
        value="Amplitude",
        label="Example",
    )
    mode = mo.ui.dropdown(
        options={
            "Automatic": "auto",
            "Amplitude": "amplitude",
            "Cross-section": "cross-section",
            "No external preparation": "generic",
        },
        value="Automatic",
        label="External placement",
    )
    show_momenta = mo.ui.checkbox(value=False, label="Momentum arrows")
    show_momentum_labels = mo.ui.checkbox(value=False, label="Momentum labels qₑ")
    show_half_edge_ids = mo.ui.checkbox(value=False, label="Debug IDs")
    mo.hstack(
        [example, mode, show_momenta, show_momentum_labels, show_half_edge_ids],
        justify="start",
        wrap=True,
        gap=1.5,
    )
    return example, mode, show_half_edge_ids, show_momentum_labels, show_momenta


@app.cell
def _(example, mo):
    dot_source = mo.ui.code_editor(
        value=example.value,
        language="text",
        min_height=280,
        max_height=520,
        debounce=400,
        label="Editable DOT",
    )
    dot_source
    return (dot_source,)


@app.cell
def _(mo):
    layout_algorithm = mo.ui.dropdown(
        options={"Force": "force", "Stable layered": "stable-layered"},
        value="Force",
        label="Layout algorithm",
    )
    custom_forces = mo.ui.checkbox(
        value=False, label="Override GammaLoop force presets"
    )
    force_steps = mo.ui.slider(
        0,
        2400,
        100,
        1200,
        debounce=True,
        show_value=True,
        label="Force iterations",
    )
    force_seed = mo.ui.slider(
        0,
        99,
        1,
        42,
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
    mo.accordion(
        {
            "Layout settings": mo.vstack(
                [
                    layout_algorithm,
                    force_steps,
                    force_seed,
                    custom_forces,
                    mo.md(
                        "The remaining sliders apply only when force overrides are enabled. "
                        "Otherwise the renderer chooses GammaLoop's amplitude or cross-section presets."
                    ),
                    directional_force,
                    spring_strength,
                    beta,
                    dangling_repulsion,
                    dangling_centroid_repulsion,
                    edge_edge_repulsion,
                    label_steps,
                ],
                gap=0.75,
            ),
        }
    )
    return (
        custom_forces,
        layout_algorithm,
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
def _(dot_source, lp, physics_dot_codec):
    try:
        graph = lp.Graph.from_dot(dot_source.value, physics_dot_codec())
        parse_error = None
    except (RuntimeError, TypeError, ValueError) as error:
        graph = None
        parse_error = f"{type(error).__name__}: {error}"
    return graph, parse_error


@app.cell
def _(
    custom_forces,
    diagram_render_settings,
    dot_source,
    drawing_root,
    force_simulation,
    json,
    layout_algorithm,
    mode,
    show_half_edge_ids,
    show_momentum_labels,
    show_momenta,
    typst,
):
    try:
        _config = diagram_render_settings(
            layout_algorithm.value,
            force_simulation=force_simulation,
            custom_forces=custom_forces.value,
            mode=mode.value,
            show_half_edge_ids=show_half_edge_ids.value,
            show_momentum_labels=show_momentum_labels.value,
            show_momenta=show_momenta.value,
        )
        _config_json = json.dumps(
            json.dumps(_config, ensure_ascii=False), ensure_ascii=False
        )
        typst_source = (
            '#import "drawings/templates/figure.typ": render\n\n'
            f'#render(json(bytes({_config_json})) + (data-path: "/graph.dot",))\n'
        )
        (drawing_root / "graph.dot").write_text(dot_source.value, encoding="utf-8")
        (drawing_root / "main.typ").write_text(typst_source, encoding="utf-8")
        _svg = typst.compile(
            str(drawing_root / "main.typ"),
            root=str(drawing_root),
            format="svg",
            package_path=str(drawing_root / "typst-packages"),
            package_cache_path=str(drawing_root / "typst-packages"),
        )
        if isinstance(_svg, list):
            rendered_svg = "".join(_page.decode("utf-8") for _page in _svg)
        else:
            rendered_svg = _svg.decode("utf-8")
        render_error = None
    except (OSError, RuntimeError, TypeError, ValueError) as error:
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
                    "Edge ID": _edge.index,
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
                This table uses Linnet's Python DOT parser for inspection.
                The drawing receives the original DOT directly, so its styling
                and external placement come entirely from GammaLoop's Typst renderer.
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


if __name__ == "__main__":
    app.run()
