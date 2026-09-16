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
    return


@app.cell
def _():
    import pydot
    import symbolica as sy
    from symbolica.community import idenso, spenso
    from symbolica.community.spenso import (
        AUTO,
        DisplaySettings,
        Representation,
        Tensor,
        TensorExpression,
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
        TensorExpression,
        TensorName,
        as_tensor,
        dot,
        idenso,
        pydot,
        spenso,
        sy,
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
    return


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
def _(AUTO, Representation, TensorExpression, TensorName, dot, trace):
    mink = Representation.mink(4)
    bis = Representation.bis(4)
    mu = mink("mu")
    nu = mink("nu")

    p = TensorName.vector("p", is_linear=True, tags=["kinematics"])
    q = TensorName.vector("q", is_linear=True, tags=["kinematics"])
    gamma = TensorExpression.gamma(4)
    mass = TensorName("m")()

    kinematic_factor = dot(p(1, mink), q(2, mink)) + mass * mass
    dirac_trace = trace(
        bis,
        gamma(AUTO, AUTO, mu),
        gamma(AUTO, AUTO, nu),
    )
    amplitude = kinematic_factor * dirac_trace
    return amplitude, mu, nu, p, q


@app.cell
def _(mo):
    def source_block(source, language=""):
        return mo.md(f"```{language}\n{source}\n```")

    return (source_block,)


@app.cell
def _(Representation, TensorExpression, TensorName, sy):
    # `raw` keeps compact rank-one vectors in their contextual form. That is
    # what lets a vector appear inside a tensor port or a gamma factor, just as
    # it does in the Typst notation examples.
    def raw(head, *arguments):
        def atom(value):
            return value.to_expression() if hasattr(value, "to_expression") else value

        return atom(head)(*(atom(argument) for argument in arguments))

    atlas_mink = Representation.mink(4)
    atlas_euc = Representation.euc(3)
    atlas_bis = Representation.bis(4)
    atlas_lor = Representation("Lor", 4, is_self_dual=False)
    atlas_cof = Representation.cof(3)
    atlas_coad = Representation.coad(8)
    atlas_cos = Representation.cos(6)

    atlas_mu = atlas_mink(sy.S("mu"))
    atlas_nu = atlas_mink(sy.S("nu"))
    atlas_a = atlas_bis("a")
    atlas_b = atlas_bis("b")

    atlas_T = TensorName("T")
    atlas_A = TensorName("A")
    atlas_H = TensorName("H")
    atlas_D = TensorName("D")
    atlas_Palette = TensorName("Palette")
    atlas_p = TensorName.vector("p", is_linear=True, tags=["kinematics"])
    atlas_q = TensorName.vector("q", is_linear=True, tags=["kinematics"])
    atlas_chi = TensorName.vector("chi", is_linear=True)

    atlas_gamma = TensorName.gamma().to_expression()
    atlas_gamma0 = sy.S("spenso::gamma0")
    atlas_gamma5 = TensorName.gamma5().to_expression()
    atlas_projp = TensorName.projp().to_expression()
    atlas_in = sy.S("spenso::in")
    atlas_out = sy.S("spenso::out")

    atlas_p1 = raw(atlas_p, 1, atlas_mink)
    atlas_q2 = raw(atlas_q, 2, atlas_mink)
    atlas_dot = raw(sy.S("spenso::dot"), atlas_p1, atlas_q2)
    atlas_interleaved = atlas_T(atlas_mu, atlas_a, atlas_nu, atlas_b)
    atlas_layout_expression = raw(
        atlas_A,
        atlas_lor(1),
        raw(atlas_p, 1, atlas_lor),
        atlas_lor(2),
        atlas_lor(3).dual(),
    )
    atlas_two_mink = raw(
        atlas_A,
        atlas_mu,
        atlas_p1,
        atlas_nu,
        atlas_q2,
        atlas_mink(3),
    )
    atlas_heterogeneous = raw(
        atlas_H,
        atlas_mu,
        atlas_p1,
        atlas_a,
        raw(atlas_chi, 2, atlas_bis),
    )
    atlas_mixed_polarity = raw(
        atlas_D,
        atlas_mink("i"),
        raw(atlas_p, 1, atlas_lor),
        atlas_mink("j"),
        raw(atlas_q, 2, atlas_lor.dual()),
    )

    atlas_chain_factors = (
        raw(atlas_gamma, atlas_in, atlas_out, atlas_mu),
        raw(atlas_gamma, atlas_in, atlas_out, atlas_p1),
        raw(atlas_gamma, atlas_in, atlas_out, atlas_nu),
    )
    atlas_open_chain = raw(
        sy.S("spenso::chain"),
        atlas_bis("u"),
        atlas_bis("v"),
        *atlas_chain_factors,
    )
    atlas_explicit_chain = raw(
        sy.S("spenso::chain"),
        atlas_bis("a"),
        atlas_bis("b"),
        *atlas_chain_factors,
    )
    atlas_trace = raw(
        sy.S("spenso::trace"),
        atlas_bis,
        raw(sy.S("spenso::cyclic"), *atlas_chain_factors),
    )

    atlas_colour = TensorExpression.t(8, 3)(
        atlas_coad("A"),
        atlas_cof("i"),
        atlas_cof("j").dual(),
    )
    atlas_palettes = atlas_Palette(
        atlas_mink(1),
        atlas_mink(5),
        atlas_euc(1),
        atlas_bis(1),
        atlas_cof(1),
        atlas_coad(1),
        atlas_cos(1),
        atlas_mink("zeta_7"),
    )
    atlas_label_kinds = atlas_T(atlas_mink("mu"), atlas_mink(sy.S("mu")))

    atlas_groups = {
        "Ports and compact vectors": (
            ("Compact vector dot product", atlas_dot),
            ("Interleaved Lorentz and spinor ports", atlas_interleaved),
            ("Two Lorentz ports with compact vectors", atlas_two_mink),
            ("Heterogeneous bra / ket ports", atlas_heterogeneous),
            ("Mixed base and dual rows", atlas_mixed_polarity),
        ),
        "Chains and traces": (
            ("Explicit gamma tensor", raw(atlas_gamma, atlas_a, atlas_b, atlas_mu)),
            ("Open gamma chain", atlas_open_chain),
            ("End-labelled gamma chain", atlas_explicit_chain),
            ("Closed gamma trace", atlas_trace),
        ),
        "Dirac and colour heads": (
            ("Gamma zero", raw(atlas_gamma0, atlas_in, atlas_out)),
            ("Gamma five", raw(atlas_gamma5, atlas_in, atlas_out)),
            ("Positive chiral projector", raw(atlas_projp, atlas_in, atlas_out)),
            ("Colour generator", atlas_colour),
        ),
        "Representations and labels": (
            ("Automatic index palettes and wrapping", atlas_palettes),
            ("Identifier versus literal math label", atlas_label_kinds),
        ),
    }
    return atlas_groups, atlas_layout_expression


@app.cell(hide_code=True)
def _(
    DisplaySettings,
    atlas_groups,
    atlas_layout_expression,
    commas,
    display_settings,
    factor_gap,
    index_gap,
    mo,
    parentheses,
    show_dimensions,
    spenso,
    symbol_scripts,
):
    def card(label, expression, settings=display_settings):
        return mo.vstack(
            [mo.md(f"**{label}**"), mo.Html(spenso.to_html(expression, settings=settings))],
            gap=0.5,
        )

    def gallery(cases):
        rows = []
        for start in range(0, len(cases), 2):
            row = [card(*case) for case in cases[start : start + 2]]
            if len(row) == 1:
                row.append(mo.md(""))
            rows.append(mo.hstack(row, widths="equal", align="start", gap=2))
        return mo.vstack(rows, gap=2)

    def layout_settings(tensor_layout):
        return DisplaySettings(
            tensor_layout=tensor_layout,
            show_dimensions=show_dimensions.value,
            parentheses=parentheses.value,
            commas=commas.value,
            symbol_scripts=symbol_scripts.value,
            index_gap=index_gap.value,
            factor_gap=factor_gap.value,
        )

    comparison = mo.hstack(
        [
            card("Ports", atlas_layout_expression, layout_settings("ports")),
            card("Schoonschip", atlas_layout_expression, layout_settings("schoonschip")),
            card("Function call", atlas_layout_expression, layout_settings("call")),
        ],
        widths="equal",
        align="start",
        gap=2,
    )
    mo.vstack(
        [
            mo.md(
                "## 1. Notation atlas\n\n"
                "The layout selector now has a genuinely mixed tensor to work on. "
                "The comparison holds the atom fixed; only its notation changes."
            ),
            mo.ui.tabs(
                {
                    "Layout comparison": comparison,
                    "Ports and vectors": gallery(atlas_groups["Ports and compact vectors"]),
                    "Chains and traces": gallery(atlas_groups["Chains and traces"]),
                    "Dirac and colour": gallery(atlas_groups["Dirac and colour heads"]),
                    "Representations": gallery(atlas_groups["Representations and labels"]),
                }
            ),
            mo.accordion(
                {
                    "Why the compact examples are raw atoms": mo.md(
                        "A compact vector is contextual notation: it lives inside a "
                        "tensor or gamma factor rather than exposing its own port. "
                        "The cell above constructs those exact Symbolica atoms and "
                        "sends them straight to the same renderer."
                    )
                }
            ),
        ],
        gap=1,
    )
    return


@app.cell(hide_code=True)
def _(amplitude, display_settings, mo, show_dimensions, source_block):
    amplitude_html = amplitude.to_html(settings=display_settings)
    amplitude_typst = amplitude.to_typst(show_dimensions.value)
    amplitude_atom = amplitude.to_expression().to_canonical_string()
    mo.vstack(
        [
            mo.md(
                "## 2. Typed construction and three layouts\n\n"
                "The scalar product, rank-zero mass, and closed gamma chain remain "
                "one `TensorExpression`. Only the presentation changes.\n\n"
                "The exact atom intentionally exposes `in` and `out`: `AUTO` "
                "materializes the two chain endpoints. Gamma follows its storage "
                "order everywhere: `(spinor-in, spinor-out, Lorentz)`."
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
    return


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
                "## 3. Idenso transformation, then Spenso reinference\n\n"
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
    return


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
    return


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
                "## 4. Concrete tensors and an executable network\n\n"
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
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ---

    Edit the constructors, add an Idenso pass, or replace the dense vectors
    above. The display controls remain ordinary `DisplaySettings`, so the same
    code works in scripts and notebooks outside Marimo.
    """)
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


@app.cell
def _():
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


@app.cell
def _():
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


@app.cell
def _():
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


@app.cell
def _():
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


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
