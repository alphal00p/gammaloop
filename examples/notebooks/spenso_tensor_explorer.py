"""Explore a concrete rank-three spinor current and its matrix slices."""

# ruff: noqa: B018, PLR1711 -- marimo uses final expressions and explicit cell returns.

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    from symbolica.community.spenso import (
        DisplaySettings,
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
    )

    return (
        DisplaySettings,
        Representation,
        TensorExpression,
        TensorLibrary,
        TensorName,
        mo,
    )


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
def _(mo):
    show_static_matrix = mo.ui.checkbox(
        label="Also show the static mathematical matrix"
    )
    show_static_matrix
    return (show_static_matrix,)


@app.cell
def _(DisplaySettings, kernel, mo, show_static_matrix):
    mo.stop(not show_static_matrix.value)
    kernel.formatted(settings=DisplaySettings(tensor_view="matrix"))
    return


if __name__ == "__main__":
    app.run()
