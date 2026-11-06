"""Exact numerator cancellation, checked with the Ordered engine, and its step cache."""

# Marimo generated return markers are intentional.
# ruff: noqa: PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · Pinches & closure")


@app.cell
def _():
    from time import perf_counter

    import marimo as mo
    from shared_gallery import Gallery
    from symbolica import E
    from symbolica.community.hepkit import (
        FeynmanDiagram,
        IntegralFamily,
        Kinematics,
        Model,
    )
    from symbolica.community.hepkit.vakint import integral_from_diagram

    return (
        E,
        FeynmanDiagram,
        Gallery,
        IntegralFamily,
        Kinematics,
        Model,
        integral_from_diagram,
        mo,
        perf_counter,
    )


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "closure",
        "Cancel a denominator. Check the reduction.",
        "Independent pinch constructions test exact reduction identities.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            Gallery.section(mo, "01 / identity", "A fast, native closure check"),
            mo.md(r"""
        $$\left(\prod_i D_i^{b_i}\right) I(a_1,\ldots,a_n)
        =I(a_1-b_1,\ldots,a_n-b_n).$$

        The numerator contains the **full massive quadratic** $D_i$, not just
        $q_i^2$. This two-loop figure-eight control uses two separately imported
        native graphs, with a physical quartic scalar vertex and an explicit unit-mass binding.
        The campaign suite applies the same idea to connected five-loop pinches,
        several cancellations, raised powers, and alternative routings.
        """),
        ]
    )
    return


@app.cell
def _(mo):
    cancellation_count = mo.ui.slider(
        start=1,
        stop=3,
        step=1,
        value=1,
        label="Factors of the propagator in the numerator",
        show_value=True,
    )
    remaining_power = mo.ui.slider(
        start=1,
        stop=3,
        step=1,
        value=1,
        label="Propagator power after cancellation",
        show_value=True,
    )
    expand_numerator = mo.ui.checkbox(label="Expand the numerator before ingress")
    closure_button = mo.ui.run_button(label="Check exact closure", kind="success")
    mo.hstack(
        [
            mo.vstack([cancellation_count, remaining_power, expand_numerator]),
            closure_button,
        ],
        justify="space-between",
    )
    return cancellation_count, closure_button, expand_numerator, remaining_power


@app.cell
def _(
    E,
    FeynmanDiagram,
    Gallery,
    IntegralFamily,
    Kinematics,
    Model,
    cancellation_count,
    closure_button,
    expand_numerator,
    integral_from_diagram,
    mo,
    perf_counter,
    remaining_power,
):
    mo.stop(not closure_button.value)
    _model = Model.phi4()
    _topology = 'digraph figure_eight { a -> a [particle="phi", lmb_id=0]; a -> a [particle="phi", lmb_id=1]; }'
    source_graph = FeynmanDiagram.from_dot(_model, _topology)
    target_graph = FeynmanDiagram.from_dot(_model, _topology)
    _source_family = source_graph.propagator_family(kinematics=Kinematics())
    _target_family = target_graph.propagator_family(kinematics=Kinematics())
    _loop = _source_family.loop_momenta[0]
    _mass_squared = (
        _source_family.kinematics.scalar_product(_loop, _loop)
        - _source_family.denominators[0]
    ).expand()
    _mass_symbols = _mass_squared.get_all_symbols()
    if len(_mass_symbols) != 1 or _mass_squared != _mass_symbols[0] ** 2:
        raise ValueError("Expected the scalar model's single mass parameter")
    _mass = _mass_symbols[0]
    _source_family = IntegralFamily(
        _source_family.loop_momenta,
        [],
        [
            denominator.replace(_mass, E("1"))
            for denominator in _source_family.denominators
        ],
        kinematics=_source_family.kinematics,
    )
    _target_family = IntegralFamily(
        _target_family.loop_momenta,
        [],
        [
            denominator.replace(_mass, E("1"))
            for denominator in _target_family.denominators
        ],
        kinematics=_target_family.kinematics,
    )
    _power = int(cancellation_count.value)
    _remaining = int(remaining_power.value)
    cancelling_numerator = _source_family.denominators[0] ** _power
    if expand_numerator.value:
        cancelling_numerator = cancelling_numerator.expand()
    source_integral = integral_from_diagram(
        source_graph,
        _source_family,
        cancelling_numerator,
        powers=[_remaining + _power, 1],
        parameter_substitutions={_mass: E("1")},
    ).to_expression()
    target_integral = integral_from_diagram(
        target_graph,
        _target_family,
        E("1"),
        powers=[_remaining, 1],
        parameter_substitutions={_mass: E("1")},
    ).to_expression()
    # Live progress and Interrupt come from the host's default progress="auto".
    _engine = Gallery.engine()
    closure_execution = Gallery.execution(_engine)
    _start = perf_counter()
    source_result = _engine.evaluate_integral(source_integral)
    target_result = _engine.evaluate_integral(target_integral)
    closure_seconds = perf_counter() - _start
    difference = (source_result - target_result).expand().together()
    if difference != E("0"):
        raise AssertionError(f"Numerator/power closure failed: {difference}")
    return (
        cancelling_numerator,
        closure_execution,
        closure_seconds,
        difference,
        source_graph,
        source_integral,
        source_result,
        target_graph,
        target_integral,
        target_result,
    )


@app.cell(hide_code=True)
def _(
    Gallery,
    cancelling_numerator,
    closure_execution,
    closure_seconds,
    difference,
    mo,
    source_graph,
    source_integral,
    source_result,
    target_graph,
    target_integral,
    target_result,
):
    # Native static MathML survives tab remounts without widget comm models.
    # These values are short (a few items), so each preview is complete.
    closure_view = mo.ui.tabs(
        {
            "The identity": mo.vstack(
                [
                    mo.md("**Cancelling numerator**"),
                    Gallery.native_formula(mo, cancelling_numerator),
                    Gallery.legend(
                        mo,
                        "dot(K(i,mink(4)),K(j,mink(4))) is the native Minkowski "
                        "product k_i·k_j of loop momenta; the mass is bound to 1.",
                    ),
                    mo.md("**Source and target integrands**"),
                    Gallery.native_formula(mo, source_integral),
                    Gallery.native_formula(mo, target_integral),
                ]
            ),
            "Retained-master results": mo.vstack(
                [
                    Gallery.native_formula(mo, source_result),
                    Gallery.native_formula(mo, target_result),
                    mo.md("**Difference**"),
                    Gallery.native_formula(mo, difference),
                ]
            ),
            "Independent native graphs": mo.hstack(
                [mo.as_html(source_graph), mo.as_html(target_graph)]
            ),
        }
    )
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    ("Exact", "difference vanishes symbolically"),
                    (f"{closure_seconds:.3f} s", "this interactive paired reduction"),
                    *Gallery.execution_metrics(closure_execution),
                ],
            ),
            closure_view,
        ]
    )
    return (closure_view,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            Gallery.section(
                mo,
                "02 / five-loop tests",
                "Strong checks within a bounded runtime",
                "Each selected paired test has a 240-second target and a 300-second deadline, with no retries.",
            ),
            Gallery.evidence_card(
                mo, "closure", "Five-loop exact closure suite", optional=True
            ),
            Gallery.note(
                mo,
                "What the suite checks",
                "Connected pinches, multiple cancelling numerators, partial and complete cancellation of dotted lines, the two dotted rung attachments, a gluon rung, alternate routings, nonunit mass, and expanded/factored numerators. Recorded guards remain part of the result.",
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            Gallery.section(
                mo,
                "03 / caching",
                "Ordered step caches: replay one-step records",
                "Vakint reduces with the Ordered scheduler only. Its persistent cache "
                "stores weight-independent one-step records, never completed results.",
            ),
            mo.md(r"""
        A replayed step is identical to a computed one, so results never depend
        on the cache; grouping and cancellations below a key are kept. The
        engine options are `ordered_cache` (`"read-write"`, `"read-only"` or
        `"off"`), `ordered_cache_directory` (None keeps records in memory) and
        `ordered_cache_ram_bytes`; `load_ordered_cache()`,
        `save_ordered_cache()`, `ordered_cache_statistics()` and
        `warm_ordered_cache([...])` manage it explicitly. Records are kept per
        reduction program, rule tables and cache version, so a cache written by
        another version is never replayed.

        Page 05 keeps one such cache for the 3L/4L ladders: rerun a case to see
        its hits replace the misses of the first run.
        """),
            Gallery.evidence_card(
                mo,
                "cache-build",
                "Accumulated depth-first cache and source provenance",
                optional=True,
            ),
            Gallery.evidence_card(
                mo,
                "held-out",
                "Held-out timings, root hits and intermediate hits",
                optional=True,
            ),
            Gallery.note(
                mo,
                "History",
                "Earlier versions also reduced depth-first with completed-result "
                "caches. That strategy and its caches were removed from Vakint: a "
                "completed terminal image replaces a frontier key and skips every "
                "grouping and cancellation below it, which the Ordered scheduler "
                "never does.",
            ),
            Gallery.footer(mo),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
