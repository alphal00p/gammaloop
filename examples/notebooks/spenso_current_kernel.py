# Run with a combined Symbolica host exposing symbolica.community.spenso:
#   python -m marimo edit examples/notebooks/spenso_current_kernel.py
# Additional packages: marimo==0.24.0, typst==0.15.0.
# The host build lives in examples/notebooks/symbolica-host; a stock Symbolica
# wheel need not expose this checkout's typed Spenso API.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium", app_title="Spenso: from a vertex to a current")


@app.cell
def _():
    import marimo as mo
    from symbolica import Expression, S
    from symbolica.community.spenso import (
        AUTO,
        Representation,
        Tensor,
        TensorExpression,
        TensorLibrary,
        TensorName,
        TensorNetwork,
        trace,
    )

    return (
        AUTO,
        Expression,
        Representation,
        S,
        Tensor,
        TensorExpression,
        TensorLibrary,
        TensorName,
        TensorNetwork,
        mo,
        trace,
    )


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

    The current API uses `Tensor.dense` and simplification methods on
    `TensorExpression`. Gamma's slots are **(row, column, Lorentz)**:
    `gamma(a, b, mu)` implements the equation above. Reversing `a, b`, as
    in the paper's listing, transposes the matrix in this contraction.
    """)


@app.cell
def _(Representation, S, Tensor, TensorLibrary, TensorName, TensorNetwork):
    spinor = Representation.bis(4)
    lorentz = Representation.mink(4)
    a, b, mu = spinor("a"), spinor("b"), lorentz("mu")
    Jbar, J = TensorName("Jbar"), TensorName("J")
    bar_components = S("bar0", "bar1", "bar2", "bar3")
    components = S("j0", "j1", "j2", "j3")

    # Unresolved representation slots define reusable component tensors.
    library = TensorLibrary.hep_lib_atom()
    library.register(Tensor.dense(Jbar(spinor), bar_components))
    library.register(Tensor.dense(J(spinor), components))
    gamma = library["spenso::gamma"]

    # Repeated spinor labels contract; mu is the only remaining port.
    current = Jbar(a) * gamma(a, b, mu) * J(b)
    network = TensorNetwork(current, library=library)
    network.execute(library=library)
    kernel = network.result_tensor(library=library)
    assert kernel.structure().interface == (mu,)
    assert len(kernel) == 4
    return bar_components, components, current, kernel, lorentz


@app.cell(hide_code=True)
def _(current, kernel, mo):
    mo.vstack(
        [
            mo.md("## The indexed rule and its four component expressions"),
            mo.Html(current.to_html()),
            mo.Html(kernel.to_html()),
            mo.md(
                "`kernel[:]` reads the components in logical order, μ = 0, 1, 2, 3. "
                "`to_dense()` only changes storage in place; it does not return a list."
            ),
        ]
    )


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


@app.cell
def _(Expression, bar_components, components, kernel):
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
    return (numeric_current,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## An Idenso identity, before component evaluation

    A single gamma matrix needs no Dirac simplification. To show a useful
    symbolic rewrite, close a two-gamma chain instead:
    $\mathrm{tr}(\gamma^\mu\gamma^\nu)=4g^{\mu\nu}$.
    `AUTO` supplies the spinor endpoints that `trace` closes.
    """)


@app.cell
def _(AUTO, Representation, TensorExpression, lorentz, mo, trace):
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
    return (reduced_trace,)


if __name__ == "__main__":
    app.run()
