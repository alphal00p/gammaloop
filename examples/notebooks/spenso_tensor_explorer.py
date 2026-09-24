"""Explore a concrete rank-three spinor current and its matrix slices."""

# ruff: noqa: B018, PLR1711 -- marimo uses final expressions and explicit cell returns.

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    from symbolica import S
    from symbolica.community.spenso import (
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
    )

    return Representation, S, TensorExpression, TensorLibrary, TensorName, mo


@app.cell
def _(Representation, TensorExpression, TensorName):
    spinor = Representation.bis(4)
    Jbar = TensorName("Jbar", print={"typst": "macron(J)", "latex": r"\bar{J}"})(spinor)
    J = TensorName("J")(spinor)
    gamma = TensorExpression.gamma(4)
    current = Jbar(1) * gamma(1, 2, 1) * gamma(2, 3, 2) * gamma(3, 4, 3) * J(4)
    current
    return (current,)


@app.cell
def _(TensorLibrary, current):
    library = TensorLibrary.hep_lib_atom()
    network = current.to_network(library=library)
    network.execute(library=library)
    kernel = network.result_tensor(library=library)
    kernel
    return (kernel,)


@app.cell
def _(Representation, mo):
    _lorentz = Representation.mink(4)
    mo.hstack(
        [mo.as_html(_lorentz.name), mo.as_html(_lorentz), mo.as_html(_lorentz(1))],
        wrap=True,
    )
    return


@app.cell
def _(TensorExpression):
    TensorExpression.gamma(4)(1, 2, 1).structure
    return


@app.cell
def _(Representation, S, TensorName):
    TensorName("explorer::A")(
        S("x"), 7, Representation.mink(4)(1), Representation.bis(4)
    ).structure
    return


@app.cell
def _(kernel):
    kernel.structure
    return


if __name__ == "__main__":
    app.run()
