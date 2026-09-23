# Run with a combined Symbolica host exposing symbolica.community.spenso:
#   python -m marimo edit examples/notebooks/spenso_current_kernel.py
# Additional packages: marimo==0.24.0, typst==0.15.0.
# The host build lives in examples/notebooks/symbolica-host; a stock Symbolica
# wheel need not expose this checkout's typed Spenso API.

import marimo

__generated_with = "0.24.0"
app = marimo.App(
    width="medium",
    app_title="Spenso: from a vertex to a current",
)


@app.cell
def _():
    import marimo as mo
    from symbolica import Expression
    from symbolica.community.spenso import AUTO, trace

    return AUTO, Expression, mo, trace


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # From a vertex to a current with Spenso

    Derive the vector current of a fermion–antifermion–vector vertex, with
    the coupling omitted:

    \[
      K^\mu = \sum_{a,b=0}^3 \bar J_a\,(\gamma^\mu)_{ab}\,J_b.
    \]

    This is the example in the **Deriving current kernels with Spenso**
    section of the [pyAmpliCol paper](https://github.com/ValentinHirschi/pyAmpliColPaper/blob/074fa4f78de9a52ef773ee11e452675502474a62/sections/implementation-design.tex#L187).
    We supply all eight symbolic input components, contract the spinor
    indices, and keep the Lorentz index open. Here `Jbar` is an already
    barred current supplied independently of `J`; no conjugation is implicit.

    The current API generates unknown components automatically and provides
    simplification methods on
    `TensorExpression`. Gamma's slots are **(row, column, Lorentz)**:
    `gamma(a, b, mu)` implements the equation above. Reversing `a, b`, as
    in the paper's listing, transposes the matrix in this contraction.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Suggested paper listing — copy the next cell

    The next cell is self-contained Python. Copy its body into the paper;
    no other notebook cell or Marimo helper is needed. Unknown currents acquire
    symbolic components automatically. The exact HEP library supplies the Weyl
    matrices; `components()` exposes the same generated input symbols afterward.

    **Suggested caption:** Deriving the four symbolic components of
    $K^\mu=\bar J\gamma^\mu J$ with Spenso. Repeated spinor indices are
    contracted, while the Lorentz index remains open.
    """)
    return


@app.cell
def _():
    from symbolica.community.spenso import (
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
    )

    spinor = Representation.bis(4)
    Jbar = TensorName("Jbar")(spinor)
    J = TensorName("J")(spinor)
    gamma = TensorExpression.gamma(4)

    # Gamma slots: row, column, Lorentz. Only mu remains open.
    current = Jbar("a") * gamma("a", "b", "mu") * J("b")
    library = TensorLibrary.hep_lib_atom()
    network = current.to_network(library=library)
    network.execute(library=library)
    kernel = network.result_tensor(library=library)
    bar_components, components = Jbar.components(), J.components()
    print(kernel)  # Components in the order mu = 0, 1, 2, 3.
    return (
        Representation,
        TensorExpression,
        bar_components,
        components,
        kernel,
        network,
    )


@app.cell(hide_code=True)
def _(kernel, mo, network):
    mo.vstack(
        [
            mo.md("## The indexed rule and its four component expressions"),
            mo.Html(network.structure().to_html()),
            mo.Html(kernel.to_html()),
            mo.md(
                "`kernel[:]` reads the components in logical order, μ = 0, 1, 2, 3. "
                "`to_dense()` only changes storage in place; it does not return a list."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Check every component exactly

    The built-in HEP library uses Weyl matrices with metric
    $g=\mathrm{diag}(1,-1,-1,-1)$:
    $\gamma^0=\left(\begin{smallmatrix}0&I\\I&0\end{smallmatrix}\right)$ and
    $\gamma^i=\left(\begin{smallmatrix}0&\sigma^i\\-\sigma^i&0\end{smallmatrix}\right)$.
    An independent matrix calculation checks all four symbolic polynomials,
    including their signs and the imaginary entries of $\gamma^2$.
    """)
    return


@app.cell
def _(Expression, Representation, bar_components, components, kernel):
    assert kernel.structure().interface == (Representation.mink(4)("mu"),)
    assert len(kernel) == 4
    _i = Expression.I
    _weyl = (
        ((0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 0), (0, 1, 0, 0)),
        ((0, 0, 0, 1), (0, 0, 1, 0), (0, -1, 0, 0), (-1, 0, 0, 0)),
        ((0, 0, 0, -_i), (0, 0, _i, 0), (0, _i, 0, 0), (-_i, 0, 0, 0)),
        ((0, 0, 1, 0), (0, 0, 0, -1), (-1, 0, 0, 0), (0, 1, 0, 0)),
    )
    for _mu, _matrix in enumerate(_weyl):
        _expected = sum(
            bar_components[_a] * _matrix[_a][_b] * components[_b]
            for _a in range(4)
            for _b in range(4)
        )
        assert (kernel[_mu] - _expected).expand() == 0
    symbolic_check = "All four components agree exactly with the Weyl matrices."
    return (symbolic_check,)


@app.cell(hide_code=True)
def _(mo):
    _intro = mo.md("""
    ## Reuse the numerical kernel

    Build the evaluator once, then vary the input without repeating the tensor
    contraction. For `Jbar = [1, 2, 3, 4]` and `J = [5, 6, 7, 8]`, the result
    is `[62, -16, 4j, 0]`. The slider rescales `J`.
    """)
    scale = mo.ui.slider(-2, 2, step=0.25, value=1, label="Scale of J", show_value=True)
    mo.vstack([_intro, scale])
    return (scale,)


@app.cell
def _(bar_components, components, kernel):
    evaluator = kernel.evaluator(
        constants={},
        funs={},
        params=[*bar_components, *components],
        n_cores=1,
    )
    return (evaluator,)


@app.cell
def _(evaluator, mo, scale, symbolic_check):
    _inputs = [1, 2, 3, 4, *[scale.value * _x for _x in (5, 6, 7, 8)]]
    numeric_current = evaluator.evaluate_complex([_inputs])[0]
    assert numeric_current[:] == [scale.value * _x for _x in (62, -16, 4j, 0)]
    mo.vstack(
        [
            mo.Html(numeric_current.to_html()),
            mo.callout(
                mo.md(symbolic_check + " Numerical evaluation also checked."),
                kind="success",
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## An Idenso identity, before component evaluation

    A single gamma matrix needs no Dirac simplification. To show a useful
    symbolic rewrite, close a two-gamma chain instead:
    $\mathrm{tr}(\gamma^\mu\gamma^\nu)=4g^{\mu\nu}$.
    `AUTO` supplies the spinor endpoints that `trace` closes.
    """)
    return


@app.cell
def _(AUTO, Representation, TensorExpression, mo, trace):
    lorentz = Representation.mink(4)
    _gamma = TensorExpression.gamma(4)
    gamma_trace = trace(
        Representation.bis(4),
        _gamma(AUTO, AUTO, lorentz("mu")),
        _gamma(AUTO, AUTO, lorentz("nu")),
    )
    reduced_trace = gamma_trace.simplify_gamma().simplify_metrics()
    assert (
        reduced_trace.to_expression() - 4 * lorentz.g("mu", "nu").to_expression()
    ).expand() == 0
    mo.vstack([mo.Html(gamma_trace.to_html()), mo.Html(reduced_trace.to_html())])
    return


if __name__ == "__main__":
    app.run()
