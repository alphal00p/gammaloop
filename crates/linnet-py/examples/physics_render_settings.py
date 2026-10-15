# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "symbolica==3.0.1",
#     "marimo==0.24.0",
# ]
# ///

# ruff: noqa: B018, PLR1711  # Cell outputs and empty returns are Marimo syntax.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="FeynKit DOT rendering")


@app.cell
def _():
    import html

    from symbolica.community import graph as lp
    import marimo as mo
    import symbolica.community.hepkit as fk

    # The Standard Model is embedded, so native and browser runs need no model file.
    # Native sessions and browser exports load identical model definitions.
    model = fk.Model.standard_model()
    return fk, html, lp, mo, model


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # FeynKit DOT rendering

    Edit physics DOT to inspect a Standard Model diagram. FeynKit draws the graph
    directly to SVG and typesets its labels with Typst. Incoming and outgoing legs,
    particles, momentum labels, and interaction details come from the same diagram.

    ImPrEd layout runs the requested passes in order. A second pass continues from
    the first drawing. Layered layout is also available. Portable Typst exports
    embed the captured SVG, preserving the preview's geometry.
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
    split_initial_state = mo.ui.checkbox(value=True, label="Split initial state")
    show_momenta = mo.ui.checkbox(value=False, label="Momentum arrows")
    show_momentum_labels = mo.ui.checkbox(value=False, label="Momentum labels")
    show_indices = mo.ui.checkbox(value=False, label="Node and edge indices")
    mo.hstack(
        [
            example,
            split_initial_state,
            show_momenta,
            show_momentum_labels,
            show_indices,
        ],
        wrap=True,
    )
    return (
        example,
        split_initial_state,
        show_momenta,
        show_momentum_labels,
        show_indices,
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
def _(mo):
    layout_algorithm = mo.ui.dropdown(
        options=["ImPrEd", "Layered"], value="ImPrEd", label="Layout"
    )
    first_steps = mo.ui.slider(
        0, 500, 10, 100, debounce=True, show_value=True, label="First pass steps"
    )
    second_steps = mo.ui.slider(
        0, 500, 10, 50, debounce=True, show_value=True, label="Second pass steps"
    )
    spacing = mo.ui.slider(
        1.0, 5.0, 0.1, 2.4, debounce=True, show_value=True, label="Spacing"
    )
    mo.accordion(
        {
            "Layout settings": mo.vstack(
                [layout_algorithm, first_steps, second_steps, spacing]
            )
        }
    )
    return layout_algorithm, first_steps, second_steps, spacing


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
    fk,
    lp,
    layout_algorithm,
    first_steps,
    second_steps,
    spacing,
    show_momenta,
    show_momentum_labels,
    show_indices,
    split_initial_state,
):
    if layout_algorithm.value == "ImPrEd":
        layouts = lp.LayoutSettings(
            impred_steps=first_steps.value, impred_spacing=spacing.value
        ).then(impred_steps=second_steps.value, impred_spacing=spacing.value)
    else:
        layouts = lp.LayoutSettings(algorithm=lp.LayoutAlgorithm.Dot)
    render_config = lp.RenderSettings(layouts=layouts)
    diagram_style = fk.DiagramStyle(
        split_initial_state=split_initial_state.value,
        momentum_arrows=show_momenta.value,
        show_momentum=show_momentum_labels.value,
        show_node_index=show_indices.value,
        show_edge_index=show_indices.value,
    )
    return render_config, diagram_style


@app.cell
def _(diagram, diagram_style, fk, render_config):
    rendered_html = None
    typst_source = None
    render_error = None
    if diagram is not None:
        try:
            drawing = diagram.render(config=render_config, style=diagram_style)
            rendered_html = drawing.to_html()
            typst_source = drawing.to_linnest()
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
                "Portable Typst export": mo.vstack(
                    [
                        mo.md(
                            "Portable Typst source embeds the captured SVG. "
                            "Graph layout and drawing have already run natively."
                        ),
                        mo.ui.code_editor(
                            value=typst_source,
                            language="text",
                            disabled=True,
                            min_height=320,
                            max_height=700,
                            label="Portable Typst (read-only)",
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
