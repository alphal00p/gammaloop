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
    _loops = (Gallery.catalogue_summary() or {}).get("loops", {})

    def _option(loop, cost):
        _size = Gallery.size_text(_loops.get(loop, {}).get("bytes"))
        _masters = _loops.get(loop, {}).get("masters")
        _facts = ", ".join(
            fact
            for fact in (
                f"{_masters} masters" if _masters else None,
                f"{_size} unpacked" if _size else None,
                cost,
            )
            if fact
        )
        return f"1L–{loop}L · adds the {loop}L pack ({_facts})"

    _four = _option(4, "about 10 s")
    _options = {_four: 4}
    if Gallery.five_loop_load_enabled():
        _options[_option(5, "about a minute and over 10 GB of memory")] = 5
    catalogue_choice = mo.ui.radio(
        _options, value=_four, label="Official catalogue to load"
    )
    load_catalogue = mo.ui.run_button(label="Load this catalogue", kind="success")
    _interrupt = (
        " In a marimo edit session the Interrupt button (Ctrl-C in a terminal) "
        "stops a load; the previous registry then stays active."
        if Gallery.interruptible()
        else ""
    )
    # Directory names only: the page never shows the server's paths.
    _source = (
        "the operator's external-data checkout "
        f"({Gallery.data_directory().name}; validated locally, never downloaded)"
        if Gallery.configured_checkout()
        else "the official release packs, installed once into the download "
        f"cache ({Gallery.data_directory().name}), never into the server's "
        "working directory"
    )
    mo.vstack(
        [
            Gallery.note(
                mo,
                "Offline by default",
                "The built-in 1L–3L registry is ready immediately. Loading is "
                f"explicit and runs in this page. It reads {_source}. Loading reads "
                "data only; it reduces nothing." + _interrupt,
            ),
            catalogue_choice,
            load_catalogue,
            mo.md("")
            if Gallery.five_loop_load_enabled()
            else Gallery.note(
                mo,
                "No five-loop catalogue on this server",
                "Loading the 5L pack runs inside this page for about a minute and "
                "holds over 10 GB for the session, so the operator enables it with "
                "VAKINT_GALLERY_ALLOW_5L=1. The five-loop page explores the 5L plan "
                "and a recorded run without it.",
                pending=True,
            ),
        ]
    )
    return catalogue_choice, load_catalogue


@app.cell
def _(Gallery, catalogue_choice, load_catalogue, mo):
    engine = Gallery.engine()
    load_report = load_error = None
    if load_catalogue.value:
        # A live status line shows while the native load runs (the host's
        # default progress="auto" uses marimo's own status spinner). A failed or
        # interrupted load leaves the engine's previous registry active.
        try:
            load_report = Gallery.load_data(engine, int(catalogue_choice.value))
        except (OSError, ValueError, RuntimeError) as _error:
            load_error = str(_error)
    constant_symbols = engine.constant_symbols()
    master_symbols = engine.master_symbols()
    master_options = Gallery.master_options(master_symbols)
    return (
        constant_symbols,
        engine,
        load_error,
        load_report,
        master_options,
        master_symbols,
    )


@app.cell(hide_code=True)
def _(Gallery, constant_symbols, load_error, load_report, master_symbols, mo):
    _loops = sorted(item["loop_count"] for item in (load_report or {}).get("loops", []))
    mo.vstack(
        [
            mo.callout(
                f"Loading failed: {load_error}. The built-in 1L–3L registry stays "
                "active. An operator selects another release with "
                "VAKINT_GALLERY_DATA_VERSION or a local checkout with "
                "VAKINT_GALLERY_ARTIFACTS.",
                kind="warn",
            )
            if load_error
            else mo.md(""),
            Gallery.metrics(
                mo,
                [
                    (len(constant_symbols), "registered constant identities"),
                    (len(master_symbols), "loaded master identities"),
                    (
                        f"1L–{max(_loops)}L" if _loops else "1L–3L",
                        "active catalogue",
                    ),
                    (
                        load_report["version"]
                        if load_report
                        else Gallery.builtin_data_version() or "built-in",
                        "external-data release",
                    ),
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
    constant_symbol = mo.ui.dropdown(
        _constant_options,
        value=_default,
        searchable=True,
        full_width=True,
        label="Constant identity",
    )
    constant_precision = mo.ui.number(
        start=2, stop=21000, value=32, step=1, label="Requested decimal digits"
    )
    mo.hstack([constant_symbol, constant_precision], widths=[3, 1], align="end")
    return constant_precision, constant_symbol


@app.cell
def _(ConstantSymbol, constant_precision, constant_symbol, engine):
    constant = engine.get_constant(
        ConstantSymbol(constant_symbol.value), prec=int(constant_precision.value)
    )
    constant
    return (constant,)


@app.cell(hide_code=True)
def _(Gallery, constant, mo):
    # The exact definition in the live pager, as this cell's own output.
    mo.stop(constant.definition is None)
    mo.vstack(
        [
            mo.md("**Exact definition**"),
            Gallery.formula(mo, constant.definition),
        ]
    )
    return


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
def _(master_options, mo):
    master_symbol = mo.ui.dropdown(
        master_options,
        value=next(iter(master_options)),
        searchable=True,
        full_width=True,
        label="Master identity (loops · family · propagator powers)",
    )
    master_precision = mo.ui.number(
        start=2, stop=21000, value=32, step=1, label="Requested decimal digits"
    )
    mo.hstack([master_symbol, master_precision], widths=[3, 1], align="end")
    return master_precision, master_symbol


@app.cell
def _(MasterSymbol, engine, master_precision, master_symbol):
    master = engine.get_master(
        MasterSymbol(master_symbol.value), prec=int(master_precision.value)
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
        master_panel = Gallery.frame(
            mo, master.graph.render(momenta=True).to_html(), height=420
        )
    elif master_view.value == "Canonical DOT":
        master_panel = mo.ui.code_editor(
            master.graph.to_dot(), language="dot", disabled=True
        )
    else:
        master_panel = Gallery.frame(mo, master.to_html(), height=560)
    master_panel
    return


@app.cell
def _(master, mo):
    _orders = {
        f"ε^{power}": coefficient
        for power, coefficient in sorted(master.coefficients.items())
        if coefficient.symbolic_expression.expression is not None
    }
    mo.stop(not _orders)
    master_order = mo.ui.radio(
        list(_orders),
        # The finite part when recorded, else the most singular order.
        value="ε^0" if "ε^0" in _orders else next(iter(_orders)),
        inline=True,
        label="Exact coefficient of",
    )
    master_orders = _orders
    master_order
    return master_order, master_orders


@app.cell(hide_code=True)
def _(Gallery, master_order, master_orders, mo):
    # Exact Laurent coefficients in the live pager (this cell's own output).
    Gallery.formula(
        mo, master_orders[master_order.value].symbolic_expression.expression
    )
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
def _(ConstantSymbol, Gallery, MasterSymbol, engine, mo):
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
    Gallery.table(mo, _errors)
    return


@app.cell
def _(E, Gallery, engine, mo):
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
    mo.vstack(
        [
            mo.as_html(cancellation),
            mo.callout(Gallery.short_numbers(conversion_message), kind="warn"),
        ]
    )
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
