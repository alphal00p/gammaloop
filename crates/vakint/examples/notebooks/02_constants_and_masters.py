"""Explore the engine's active registry through native inspection records."""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · Constants & masters")


@app.cell
def _():
    import marimo as mo
    from shared_gallery import Gallery
    from symbolica import E
    from symbolica.community.hepkit.vakint import ConstantSymbol, MasterSymbol

    return ConstantSymbol, E, Gallery, MasterSymbol, mo


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "inspect",
        "Meet the ingredients",
        "Inspect constants, exact Laurent coefficients, achieved precision and "
        "faithful master graphs. Every record belongs to the engine's active data snapshot.",
    )
    return


@app.cell
def _(Gallery, mo):
    load_catalogue = mo.ui.run_button(
        label="Load the official 1L–5L catalogue", kind="success"
    )
    mo.vstack(
        [
            Gallery.note(
                mo,
                "Offline by default",
                "The built-in 1L–3L registry is ready immediately. Loading 5L is explicit: a configured checkout is validated locally; otherwise the official latest packs download once and are reused on disk.",
            ),
            load_catalogue,
        ]
    )
    return (load_catalogue,)


@app.cell
def _(Gallery, load_catalogue):
    engine = Gallery.engine()
    load_report = None
    if load_catalogue.value:
        load_report = Gallery.load_data(engine, 5)
    constant_symbols = engine.constant_symbols()
    master_symbols = engine.master_symbols()
    return constant_symbols, engine, load_report, master_symbols


@app.cell(hide_code=True)
def _(Gallery, constant_symbols, load_report, master_symbols, mo):
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    (len(constant_symbols), "registered constant identities"),
                    (len(master_symbols), "loaded master identities"),
                    ("1L–5L" if load_report else "1L–3L", "active catalogue"),
                ],
            ),
            mo.accordion({"Validated load report": mo.json(load_report)})
            if load_report
            else mo.md(""),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "01 / constants",
        "Exact definition or finite source",
        "Availability describes the source; achieved precision describes the returned number.",
    )
    return


@app.cell
def _(constant_symbols, mo):
    _constant_options = [str(symbol.expression) for symbol in constant_symbols]
    _default = next(
        (label for label in _constant_options if label.endswith("z3")),
        next(iter(_constant_options)),
    )
    constant_request = (
        mo.md("{symbol}\n\n{precision}")
        .batch(
            symbol=mo.ui.dropdown(
                _constant_options,
                value=_default,
                searchable=True,
                full_width=True,
                label="Constant identity",
            ),
            precision=mo.ui.number(
                start=2, stop=21000, value=32, step=1, label="Requested decimal digits"
            ),
        )
        .form(submit_button_label="Inspect constant", show_clear_button=False)
    )
    constant_request
    return (constant_request,)


@app.cell
def _(ConstantSymbol, constant_request, engine):
    _request = constant_request.value or {
        "symbol": "z3",
        "precision": 32,
    }
    constant = engine.get_constant(
        ConstantSymbol(_request["symbol"]), prec=int(_request["precision"])
    )
    constant
    return (constant,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "02 / masters",
        "The whole Laurent series",
        "Lookup returns every available epsilon order, independently of the engine's expansion-depth setting.",
    )
    return


@app.cell
def _(master_symbols, mo):
    _master_options = [str(symbol.expression) for symbol in master_symbols]
    master_request = (
        mo.md("{symbol}\n\n{precision}")
        .batch(
            symbol=mo.ui.dropdown(
                _master_options,
                value=next(iter(_master_options)),
                searchable=True,
                full_width=True,
                label="Qualified master identity",
            ),
            precision=mo.ui.number(
                start=2, stop=21000, value=32, step=1, label="Requested decimal digits"
            ),
        )
        .form(submit_button_label="Inspect master", show_clear_button=False)
    )
    master_request
    return (master_request,)


@app.cell
def _(MasterSymbol, engine, master_request, master_symbols):
    _request = master_request.value or {
        "symbol": str(master_symbols[0].expression),
        "precision": 32,
    }
    master = engine.get_master(
        MasterSymbol(_request["symbol"]), prec=int(_request["precision"])
    )
    return (master,)


@app.cell(hide_code=True)
def _(Gallery, master, mo):
    master_view = mo.ui.tabs(
        {"Master & coefficients": "", "Native graph": "", "Canonical DOT": ""}
    )
    _status = (
        Gallery.note(
            mo,
            "Computed master",
            f"{len(master.coefficients)} recorded Laurent orders, in the stored unit-mass normalization.",
        )
        if master.is_computed
        else Gallery.note(
            mo,
            "Declared, not computed",
            "This catalogue master is a valid reduction terminal. Its empty coefficient map means no numerical Laurent data is available; it does not mean the coefficients are zero.",
            pending=True,
        )
    )
    mo.vstack([_status, master_view])
    return (master_view,)


@app.cell(hide_code=True)
def _(Gallery, master, master_view, mo):
    if master_view.value == "Native graph":
        master_panel = mo.as_html(master.graph.render(momenta=True))
    elif master_view.value == "Canonical DOT":
        master_panel = mo.ui.code_editor(
            master.graph.to_dot(), language="dot", disabled=True
        )
    else:
        master_panel = mo.Html(Gallery.math_font_css() + master.to_html())
    master_panel
    return


@app.cell(hide_code=True)
def _(engine, master, mo):
    mo.stop(master.is_computed)
    try:
        engine.numerical_evaluation(master.symbol.expression, {})
    except ValueError as _error:
        uncomputed_error = str(_error)
    else:
        raise AssertionError("A declared uncomputed master unexpectedly evaluated")
    mo.vstack(
        [
            mo.md("**Attempting numerical substitution of this declared terminal**"),
            mo.callout(uncomputed_error, kind="warn"),
        ]
    )
    return (uncomputed_error,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "03 / error contract",
        "Unknown identities and honest precision",
        "Failures carry context. An uncertain zero remains inspectable without being exported as an exact zero.",
    )
    return


@app.cell
def _(ConstantSymbol, MasterSymbol, engine, mo):
    _errors = []
    for _kind, _call in [
        (
            "Unknown constant",
            lambda: engine.get_constant(ConstantSymbol("gallery_missing_constant")),
        ),
        (
            "Unknown master",
            lambda: engine.get_master(MasterSymbol("gallery_missing_master")),
        ),
    ]:
        try:
            _call()
        except ValueError as _error:
            _errors.append({"request": _kind, "graceful error": str(_error)})
        else:
            raise AssertionError(f"{_kind} unexpectedly resolved")
    mo.ui.table(_errors, selection=None, show_download=False)
    return


@app.cell
def _(E, engine, mo):
    cancellation, _statistical_error = engine.numerical_evaluation(
        E("x-y", default_namespace="vakint"), {"x": 1.0, "y": 1.0}
    )
    try:
        cancellation.to_expression()
    except ValueError as _error:
        conversion_message = str(_error)
    else:
        raise AssertionError(
            "An uncertain zero unexpectedly became an exact expression"
        )
    mo.vstack([mo.as_html(cancellation), mo.callout(conversion_message, kind="warn")])
    return cancellation, conversion_message


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            mo.md(
                "A higher request never invents source digits. Vakint evaluates once, reports componentwise achieved accuracy, and preserves that accuracy when conversion is possible."
            ),
            Gallery.footer(mo),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
