import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="Spenso + Idenso live showcase")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Spenso + Idenso, live

    Build typed tensor expressions with **Spenso**, transform their ordinary
    Symbolica atoms with **Idenso**, and render the same semantics through the
    configurable Typst display. Change any control or edit any cell: Marimo
    recomputes only the affected section.

    The notebook uses one combined Symbolica Community extension, so Spenso,
    Idenso, and Symbolica all share the same expression type and symbol state.
    """)


@app.cell
def _():
    import pydot
    from symbolica.community import idenso, spenso
    from symbolica.community.spenso import (
        AUTO,
        DisplaySettings,
        Representation,
        Tensor,
        TensorName,
        as_tensor,
        dot,
        trace,
    )

    return (
        AUTO,
        DisplaySettings,
        Representation,
        Tensor,
        TensorName,
        as_tensor,
        dot,
        idenso,
        pydot,
        spenso,
        trace,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.callout(
        mo.md(
            "The native backend and the optional Typst renderer loaded "
            "successfully. The controls below affect every mathematical view."
        ),
        kind="success",
        title="One kernel, live display",
    )


@app.cell(hide_code=True)
def _(mo):
    layout = mo.ui.radio(
        options={
            "Ports": "ports",
            "Schoonschip": "schoonschip",
            "Function call": "call",
        },
        value="Ports",
        inline=True,
        label="Tensor layout",
    )
    show_dimensions = mo.ui.switch(value=False, label="Show dimensions")
    parentheses = mo.ui.switch(value=True, label="Parentheses")
    commas = mo.ui.switch(value=False, label="Commas")
    symbol_scripts = mo.ui.switch(value=True, label="Symbol scripts")
    index_gap = mo.ui.dropdown(
        options={"Tight": "0em", "Default": "0.08em", "Open": "0.18em"},
        value="Default",
        label="Index spacing",
    )
    factor_gap = mo.ui.dropdown(
        options={"Tight": "0.04em", "Default": "0.12em", "Open": "0.24em"},
        value="Default",
        label="Factor spacing",
    )
    mo.vstack(
        [
            layout,
            mo.hstack(
                [show_dimensions, parentheses, commas, symbol_scripts],
                justify="start",
                wrap=True,
            ),
            mo.hstack([index_gap, factor_gap], justify="start", wrap=True),
        ],
        gap=1,
    )
    return (
        commas,
        factor_gap,
        index_gap,
        layout,
        parentheses,
        show_dimensions,
        symbol_scripts,
    )


@app.cell
def _(
    DisplaySettings,
    commas,
    factor_gap,
    index_gap,
    layout,
    parentheses,
    show_dimensions,
    symbol_scripts,
):
    display_settings = DisplaySettings(
        tensor_layout=layout.value,
        show_dimensions=show_dimensions.value,
        parentheses=parentheses.value,
        commas=commas.value,
        symbol_scripts=symbol_scripts.value,
        index_gap=index_gap.value,
        factor_gap=factor_gap.value,
    )
    return (display_settings,)


@app.cell
def _(AUTO, Representation, TensorName, dot, trace):
    mink = Representation.mink(4)
    bis = Representation.bis(4)
    mu = mink("mu")
    nu = mink("nu")

    p = TensorName.vector("p", is_linear=True, tags=["kinematics"])
    q = TensorName.vector("q", is_linear=True, tags=["kinematics"])
    gamma = TensorName.gamma()(mink, bis, bis)
    mass = TensorName("m")()

    kinematic_factor = dot(p(1, mink), q(2, mink)) + mass * mass
    dirac_trace = trace(
        bis,
        gamma(mu, AUTO, AUTO),
        gamma(nu, AUTO, AUTO),
    )
    amplitude = kinematic_factor * dirac_trace
    return amplitude, mu, nu, p, q


@app.cell
def _(mo):
    def source_block(source, language=""):
        return mo.md(f"```{language}\n{source}\n```")

    return (source_block,)


@app.cell(hide_code=True)
def _(amplitude, display_settings, mo, show_dimensions, source_block):
    amplitude_html = amplitude.to_html(settings=display_settings)
    amplitude_typst = amplitude.to_typst(show_dimensions.value)
    amplitude_atom = amplitude.to_expression().to_canonical_string()
    mo.vstack(
        [
            mo.md(
                "## 1. Typed construction and three layouts\n\n"
                "The scalar product, rank-zero mass, and closed gamma chain remain "
                "one `TensorExpression`. Only the presentation changes."
            ),
            mo.ui.tabs(
                {
                    "Rendered": mo.Html(amplitude_html),
                    "Typst source (ports)": source_block(amplitude_typst, "typst"),
                    "Exact Symbolica atom": source_block(amplitude_atom, "text"),
                }
            ),
        ],
        gap=1,
    )


@app.cell
def _(TensorName, as_tensor, idenso, mu, nu, p):
    metric = TensorName.g()

    # Crossing to ordinary Symbolica explicitly leaves the repeated index visible
    # to Idenso instead of asking Spenso to choose a tensor-aware contraction.
    metric_product = metric(mu, nu).to_expression() * p(1, mu).to_expression()
    simplified_atom = idenso.simplify_metrics(metric_product)
    simplified_tensor = as_tensor(simplified_atom)
    return metric_product, simplified_atom, simplified_tensor


@app.cell(hide_code=True)
def _(
    display_settings,
    metric_product,
    mo,
    simplified_atom,
    simplified_tensor,
    source_block,
    spenso,
):
    before_metric = spenso.to_html(metric_product, settings=display_settings)
    after_metric = simplified_tensor.to_html(settings=display_settings)
    mo.vstack(
        [
            mo.md(
                "## 2. Idenso transformation, then Spenso reinference\n\n"
                "`simplify_metrics` accepts and returns an ordinary Symbolica "
                "`Expression`. `as_tensor` validates the result and restores its "
                "typed external interface."
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md("**Before**"), mo.Html(before_metric)]),
                    mo.vstack([mo.md("**After**"), mo.Html(after_metric)]),
                ],
                widths="equal",
                align="start",
                gap=2,
            ),
            mo.accordion(
                {
                    "Exact atoms": mo.hstack(
                        [
                            source_block(metric_product.to_canonical_string(), "text"),
                            source_block(simplified_atom.to_canonical_string(), "text"),
                        ],
                        widths="equal",
                        align="start",
                    )
                }
            ),
        ],
        gap=1,
    )


@app.cell
def _(as_tensor, idenso, mu, p, q):
    indexed_product = p(1, mu).to_expression() * q(2, mu).to_expression()
    dotted_atom = idenso.to_dots(indexed_product)
    dotted_tensor = as_tensor(dotted_atom)
    return dotted_atom, dotted_tensor, indexed_product


@app.cell(hide_code=True)
def _(
    display_settings,
    dotted_atom,
    dotted_tensor,
    indexed_product,
    mo,
    source_block,
):
    mo.accordion(
        {
            "A second Idenso rewrite: contracted indices → dot notation": mo.vstack(
                [
                    mo.Html(dotted_tensor.to_html(settings=display_settings)),
                    mo.hstack(
                        [
                            source_block(indexed_product.to_canonical_string(), "text"),
                            source_block(dotted_atom.to_canonical_string(), "text"),
                        ],
                        widths="equal",
                        align="start",
                    ),
                ]
            )
        }
    )


@app.cell
def _(Representation, Tensor, TensorName, display_settings, pydot):
    euc = Representation.euc(2)
    u = Tensor.dense(TensorName.vector("u")(euc), [1.0, 2.0])
    v = Tensor.dense(TensorName.vector("v")(euc), [3.0, 4.0])

    contraction = u * v
    contraction_math = contraction.to_html(settings=display_settings)
    contraction_dot = contraction.to_dot()
    graph = pydot.graph_from_dot_data(contraction_dot)[0]
    contraction_graph = graph.create_svg().decode()
    contraction.execute()
    contraction_result = contraction.result_scalar()
    return (
        contraction_dot,
        contraction_graph,
        contraction_math,
        contraction_result,
    )


@app.cell(hide_code=True)
def _(
    contraction_dot,
    contraction_graph,
    contraction_math,
    contraction_result,
    mo,
    source_block,
):
    mo.vstack(
        [
            mo.md(
                "## 3. Concrete tensors and an executable network\n\n"
                "The same multiplication operator promotes concrete tensors to a "
                "`TensorNetwork`. Its mathematical structure, execution graph, and "
                "computed scalar are all available without conflating their APIs."
            ),
            mo.ui.tabs(
                {
                    "Tensor view": mo.Html(contraction_math),
                    "Network graph": mo.Html(contraction_graph),
                    "DOT source": source_block(contraction_dot, "dot"),
                }
            ),
            mo.callout(
                mo.md(f"The contraction evaluates to **{contraction_result}**."),
                kind="success",
                title="Execution result",
            ),
        ],
        gap=1,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---

    Edit the constructors, add an Idenso pass, or replace the dense vectors
    above. The display controls remain ordinary `DisplaySettings`, so the same
    code works in scripts and notebooks outside Marimo.
    """)


if __name__ == "__main__":
    app.run()
