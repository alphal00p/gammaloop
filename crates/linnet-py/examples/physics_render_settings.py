# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "linnet==0.1.0",
#     "symbolica==3.0.0",
#     "marimo==0.24.0",
#     "typst==0.15.0",
# ]
# ///

# ruff: noqa: B018, PLR1711  # Cell outputs and empty returns are Marimo syntax.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="FeynKit DOT rendering")


@app.cell
def _():
    import html
    from pathlib import Path

    import linnet as lp
    import marimo as mo
    import symbolica.community.feynkit as fk

    # The documentation exporter bundles this Standard Model fixture.
    # Native sessions resolve the same model from this checkout.
    model = fk.Model(
        Path(__file__).resolve().parents[3]
        / "crates/feynkit-model/tests/fixtures/sm.json"
    )
    return fk, html, lp, mo, model


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # FeynKit DOT rendering

    Edit compact physics DOT. FeynKit resolves `particle` names or `pdg` codes
    against the Standard Model, infers interaction slots, and renders the parsed
    diagram with Linnet/Linnest. The preview and table inspect the same object.
    Annotated DOT from `FeynmanDiagram.to_dot()` is also accepted.

    The amplitude is a photon-scattering diagram with a top loop and gluon
    exchange. The cross-section describes photon-mediated electron–positron
    annihilation into muons. Its `is_cut` tags pair initial-state legs;
    `final_state="mu-,mu+"` selects physical cuts with those particles.
    Use `int_id` when particles alone do not identify an interaction uniquely.
    Numerators default to one; compact import does not generate Feynman rules.

    **Automatic** placement puts incoming legs on the left and outgoing legs
    on the right, with a centroid bias that separates them from the diagram.
    Cross-sections open initial-state connections into left/right legs by default;
    matching pairs share movable rows. Disable **Split initial state** to draw
    the sewn graph. Internal labels keep a uniform gap to their text bounds,
    sliding along their curves to avoid overlaps. Momentum labels use the diagram's
    stored loop momentum basis. Debug mode adds node and half-edge IDs.

    Hover for physics information; click to inspect an edge or vertex.
    Shift-click toggles selection. The outer quarters of an internal edge
    select its half-edges; the middle selects the whole edge. The transparent
    SVG follows the documentation's light or dark theme.

    Open **Layout settings** for the sliders. FeynKit's mode-specific presets
    stay active unless you enable custom force parameters. Edits render after
    typing pauses; slider changes render on release.
    """)
    return


@app.cell(hide_code=True)
def _():
    examples = {
        "Amplitude": r"""digraph GL05 {
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
}""",
        "Cross-section": r"""digraph MuonPair {
  graph [final_state="mu-,mu+"];
  ext [style=invis];
  ext -> eL [particle="e-", is_cut=0];
  ext -> eL [particle="e+", is_cut=1];
  eR -> ext [particle="e-", is_cut=0];
  eR -> ext [particle="e+", is_cut=1];

  eL -> muL [particle="a"];
  muR -> eR [particle="a"];
  muL -> muR [particle="mu-", lmb_id=0];
  muR -> muL [particle="mu-"];
}""",
    }
    return (examples,)


@app.cell(hide_code=True)
def _(examples, mo):
    example = mo.ui.dropdown(options=examples, value="Amplitude", label="Example")
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
    split_initial_state = mo.ui.checkbox(value=True, label="Split initial state")
    show_momenta = mo.ui.checkbox(value=False, label="Momentum arrows")
    show_momentum_labels = mo.ui.checkbox(value=False, label="Momentum labels")
    show_half_edge_ids = mo.ui.checkbox(value=False, label="Debug IDs")
    mo.hstack(
        [
            example,
            mode,
            split_initial_state,
            show_momenta,
            show_momentum_labels,
            show_half_edge_ids,
        ],
        justify="start",
        wrap=True,
        gap=1.5,
    )
    return (
        example,
        mode,
        show_half_edge_ids,
        show_momentum_labels,
        show_momenta,
        split_initial_state,
    )


@app.cell(hide_code=True)
def _(example, mo):
    dot_source = mo.ui.code_editor(
        value=example.value,
        language="text",
        min_height=280,
        max_height=520,
        debounce=400,
        label="Editable FeynKit DOT",
    )
    dot_source
    return (dot_source,)


@app.cell(hide_code=True)
def _(lp, mo):
    layout_algorithm = mo.ui.dropdown(
        options={
            "Force": lp.LayoutAlgorithm.Force,
            "Stable layered": lp.LayoutAlgorithm.StableLayered,
        },
        value="Force",
        label="Layout algorithm",
    )
    custom_forces = mo.ui.checkbox(value=False, label="Override FeynKit force presets")
    force_steps = mo.ui.slider(
        0,
        2400,
        100,
        100,
        debounce=True,
        show_value=True,
        label="Steps per epoch (30 epochs)",
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
                        "Otherwise the renderer chooses FeynKit's amplitude or cross-section presets."
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
def _(dot_source, fk, model):
    try:
        diagram = fk.FeynmanDiagram.from_dot(model, dot_source.value)
        parse_error = None
    except (fk.DiagramError, RuntimeError, TypeError, ValueError) as error:
        diagram = None
        parse_error = f"{type(error).__name__}: {error}"
    return diagram, parse_error


@app.cell
def _(
    beta,
    custom_forces,
    dangling_centroid_repulsion,
    dangling_repulsion,
    directional_force,
    edge_edge_repulsion,
    force_seed,
    force_steps,
    label_steps,
    layout_algorithm,
    lp,
    mode,
    show_half_edge_ids,
    show_momentum_labels,
    show_momenta,
    spring_strength,
    split_initial_state,
):
    _forces = {}
    if custom_forces.value:
        _forces = {
            "directional_force": directional_force.value,
            "spring_strength": spring_strength.value,
            "beta": beta.value,
            "dangling_repulsion": dangling_repulsion.value,
            "dangling_centroid_repulsion": dangling_centroid_repulsion.value,
            "edge_edge_repulsion": edge_edge_repulsion.value,
            "label_steps": label_steps.value,
        }
    # The shared physics template lightens sink halves by 45% and owns the
    # particle, label, arrow, and placement conventions.
    render_config = lp.RenderConfig(
        layouts=lp.LayoutOptions(
            algorithm=layout_algorithm.value,
            steps=force_steps.value,
            seed=force_seed.value,
            **_forces,
        ),
        template_options={
            "mode": mode.value,
            "split-initial-state": split_initial_state.value,
            "momentum-arrows": show_momenta.value,
            "show-momentum": show_momentum_labels.value,
            "debug": show_half_edge_ids.value,
        },
    )
    route_momenta = show_momenta.value or show_momentum_labels.value
    return render_config, route_momenta


@app.cell
def _(diagram, fk, render_config, route_momenta):
    rendered_html = None
    typst_source = None
    render_error = None
    if diagram is not None:
        try:
            rendered_html = diagram.to_html(config=render_config, momenta=route_momenta)
            typst_source = diagram.to_linnest(
                config=render_config, momenta=route_momenta
            )
        except (fk.DiagramError, OSError, RuntimeError, TypeError, ValueError) as error:
            render_error = f"{type(error).__name__}: {error}"
    return render_error, rendered_html, typst_source


@app.cell(hide_code=True)
def _(html, mo, parse_error, render_error, rendered_html):
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
            title="FeynKit rendering failed",
        )
    else:
        # The SVG manages iframe height; Marimo islands do not provide
        # mo.iframe's global resize callback.
        _output = mo.Html(
            '<div data-notebook-ready="physics_render_settings">'
            '<iframe title="Interactive FeynKit diagram" '
            'sandbox="allow-scripts allow-same-origin" '
            'style="width:100%;height:500px;border:0" '
            f'srcdoc="{html.escape(rendered_html, quote=True)}"></iframe></div>'
        )
    _output
    return


@app.cell(hide_code=True)
def _(mo, typst_source):
    if typst_source is None:
        _source_panel = mo.md(
            "Fix the DOT or rendering error to inspect the generated Typst."
        )
    else:
        _source_panel = mo.accordion(
            {
                "Live generated Typst": mo.vstack(
                    [
                        mo.md(
                            "FeynKit's Typst source for the diagram and rendering "
                            "controls above. The SVG's hover and selection behavior "
                            "is added after Typst rendering."
                        ),
                        mo.ui.code_editor(
                            value=typst_source,
                            language="text",
                            disabled=True,
                            min_height=320,
                            max_height=700,
                            label="Generated Typst (read-only)",
                        ),
                    ]
                ),
            }
        )
    _source_panel
    return


@app.cell(hide_code=True)
def _(diagram, mo):
    if diagram is None:
        _details = mo.md("Fix the DOT input to inspect the parsed diagram.")
    else:
        _details = mo.vstack(
            [
                mo.md(f"""
                ## Parsed FeynKit diagram

                **{len(diagram.vertices)} vertices**, **{len(diagram.edges)} edges**,
                **{diagram.loop_count} loops**, and **{len(diagram.cuts)} physical cuts**.
                This table reads the same `FeynmanDiagram` that produced the preview.
                """),
                mo.ui.table(
                    [
                        {
                            "Edge ID": edge.id,
                            "source": edge.source
                            if edge.source is not None
                            else "external",
                            "target": edge.target
                            if edge.target is not None
                            else "external",
                            "particle": edge.particle_name,
                            "PDG": edge.particle_pdg,
                            "external state": edge.external_state
                            if edge.is_external
                            else "internal",
                        }
                        for edge in diagram.edges
                    ],
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
