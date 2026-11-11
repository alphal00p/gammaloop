"""State of the single-scale vacuum basis, 1L-5L: what Vakint returns exactly and what numerically.

The status tables are the benchmark generator's own output
(crates/vakint/docs/benchmark/tools/overview.py, the code behind the precision
tables of benchmark.pdf), read from its JSON or run on the configured data tree;
see basis_status.py for the order of sources. The live examples evaluate
shipped master graphs with Vakint on the loaded data: the four-loop one always
(on the button, about ten seconds), the five-loop ones only where the operator
enabled the 5L catalogue (VAKINT_GALLERY_ALLOW_5L=1: about a minute and over
10 GB). VAKINT_GALLERY_RUN_EXAMPLES=1 evaluates them without a button press
(headless runs and static exports).
"""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · State of the basis")


@app.cell(hide_code=True)
def _():
    import marimo as mo

    import basis_status as bs
    from shared_gallery import Gallery

    return Gallery, bs, mo


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "basis",
        "State of the vacuum-integral basis, one to five loops",
        "Which Laurent coefficients Vakint returns exactly, which in named "
        "constants known to a finite number of digits, and which only "
        "numerically, with how many digits. The tables are those of the "
        "benchmark; the examples evaluate live on the loaded data.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, bs, mo):
    status = bs.load_status(Gallery.configured_checkout())
    overview = status["overview"]
    _tree = (overview or {}).get("tree") or {}
    data_root = Gallery.data_root()
    data_basis, data_version = bs.data_basis(data_root)
    _metrics = [
        (_tree.get("version") or "—", "data release in the tables"),
        (bs.five_loop_basis(overview) or "—", "five-loop basis in the tables"),
        (
            {
                "file": "generator output",
                "generated": "generated now",
                "committed": "committed",
            }.get(status["origin"], "—"),
            "status tables",
        ),
        (
            f"{data_version or Gallery.builtin_data_version() or 'built-in'}"
            + (f" · {data_basis}" if data_basis else ""),
            "data for the live examples",
        ),
    ]
    mo.vstack(
        [
            Gallery.metrics(mo, _metrics),
            *[mo.callout(note, kind="warn") for note in status["notes"]],
            Gallery.note(
                mo,
                "Where the tables come from",
                f"From {status['source'] or 'nowhere'}: "
                + (
                    f"{bs.tree_label(overview)}"
                    + (
                        f", plus {overview['overlay']['title']} (marked preview / pending)"
                        if overview and overview.get("overlay")
                        else ""
                    )
                    + "."
                    if overview
                    else "no generator output is available."
                )
                + " The same generator writes the precision tables of benchmark.pdf, "
                "so the two cannot disagree.",
            ),
        ]
    )
    return data_basis, data_root, data_version, overview, status


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.legend(
        mo,
        "S = exact, every constant evaluates to any precision · "
        "S*(X) = exact in named constants, X = digits of the least precise one · "
        "N(X) = numerical only, X digits (published digits are the measured claim "
        "minus one) · 0 = exact zero · — = not available · "
        "'not every input' = some inputs reach coefficients the data lacks; "
        "Vakint then fails closed at that order instead of guessing.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, bs, mo, overview, status):
    mo.stop(overview is None)
    _rows, _columns = bs.result_rows(overview)
    _notes = bs.result_notes(status["tables"], _rows)
    mo.vstack(
        [
            Gallery.section(
                mo,
                "01 / per loop and order",
                "What Vakint returns for any vacuum integral at L loops",
                "Each cell is the worst class over every input at that loop order and "
                "ε order: an exact cell is exact for every input.",
            ),
            Gallery.table(mo, _rows, _columns),
            mo.accordion(
                {
                    "Notes to the table": mo.md(
                        "\n".join(f"- **[{mark}]** {text}" for mark, text in _notes)
                        or "No notes."
                    ),
                    "Coefficient counts per loop and order (S / S* / N / 0)": Gallery.table(
                        mo, bs.summary_rows(overview)
                    ),
                }
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, bs, mo, overview):
    _basis = bs.five_loop_basis(overview)
    status_filter = mo.ui.dropdown(
        {
            "every class": "all",
            "exact (S)": "exact",
            "exact in named constants (S*)": "named",
            "numerical (N)": "numerical",
            "exact zero": "zero",
            "not available": "missing",
        },
        value="every class",
        label="Class",
    )
    order_filter = mo.ui.dropdown(
        {"any order": "any", **{bs.eps(o): str(o) for o in range(-5, 1)}},
        value="any order",
        label="at ε order",
    )
    mo.vstack(
        [
            Gallery.section(
                mo,
                "02 / five-loop masters",
                f"The {_basis or 'five-loop'} masters, order by order",
                "Filter by class and ε order. 'pending' marks a result derived and "
                "reviewed but not yet in a released data tree.",
            ),
            mo.hstack([status_filter, order_filter], justify="start"),
        ]
    )
    return order_filter, status_filter


@app.cell(hide_code=True)
def _(Gallery, bs, mo, order_filter, overview, status_filter):
    mo.stop(overview is None)
    _rows, _columns = bs.master_rows(
        overview,
        5,
        status=None if status_filter.value == "all" else status_filter.value,
        order=None if order_filter.value == "any" else int(order_filter.value),
    )
    _total = sum(1 for r in overview["table1"] if r["loop"] == 5)
    mo.vstack(
        [
            mo.md(f"**{len(_rows)} of {_total}** five-loop masters match."),
            Gallery.table(mo, _rows, _columns),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, bs, data_basis, mo):
    run_examples = mo.ui.run_button(label="Evaluate the examples", kind="success")
    example_digits = mo.ui.number(
        start=10, stop=2000, value=60, step=10, label="Requested decimal digits"
    )
    _five = Gallery.five_loop_load_enabled()
    mo.vstack(
        [
            Gallery.section(
                mo,
                "03 / live examples",
                "Evaluate it yourself on the loaded data",
                "Each input is a shipped master graph (SHA-256 checked), built with "
                "HEPKit and integral_from_diagram, then evaluated by Vakint with "
                "symbolic masters in the data's own normalization (FMFTandMATAD, "
                "unit mass). Every class below is read off Vakint's output.",
            ),
            mo.hstack([example_digits, run_examples], justify="start", align="end"),
            mo.md("")
            if _five
            else Gallery.note(
                mo,
                "Five-loop examples need the 5L catalogue",
                "Loading the 5L pack takes about a minute and over 10 GB, so the "
                "operator enables it with VAKINT_GALLERY_ALLOW_5L=1. The four-loop "
                "example runs everywhere.",
                pending=True,
            ),
            mo.md("")
            if data_basis in (None, "BASIS-v3")
            else Gallery.note(
                mo,
                f"The loaded data are {data_basis}",
                "The five-loop examples show what this release returns. The exact "
                "poles described here need BASIS-v3 with the 84 adopted pole forms "
                "(vakint-external-data v0.3.0): run the page on such a tree "
                "(VAKINT_GALLERY_ARTIFACTS) to see them.",
                pending=True,
            ),
        ]
    )
    return example_digits, run_examples


@app.cell(hide_code=True)
def _(Gallery, bs, data_root, example_digits, mo, run_examples):
    # No mo.stop here: the cells below read four_loop (None = not evaluated).
    four_loop = None
    _key = "4L-spurious"
    if not (run_examples.value or bs.autorun()):
        _view = mo.md("*Press “Evaluate the examples” (four-loop: about ten seconds).*")
    else:
        try:
            four_loop = bs.evaluate_example(
                Gallery, _key, precision=int(example_digits.value)
            )
            _poles, _count = bs.retain_poles(Gallery, _key)
        except (OSError, ValueError, RuntimeError) as _error:
            _view = mo.callout(f"The four-loop example failed: {_error}", kind="warn")
        else:
            _names = bs.master_names(data_root, 4)
            _pole_rows = [
                {
                    "master": _names.get(master, master),
                    "most singular term of its reduction coefficient": bs.eps(power),
                }
                for master, power in sorted(_poles.items())
            ]
            _view = mo.vstack(
                [
                    mo.md(f"### {bs.EXAMPLES[_key]['title']}"),
                    mo.md(bs.EXAMPLES[_key]["point"]),
                    mo.md(
                        f"Master P15 with one dot on edge 0, data {four_loop['data']}. "
                        f"Its reduction (Retain mode) has {_count} masters, "
                        f"{len(_poles)} of them with a 1/ε reduction coefficient:"
                    ),
                    Gallery.table(mo, _pole_rows),
                    Gallery.table(mo, bs.example_table_rows(four_loop)),
                ]
            )
    mo.vstack([_view])
    return (four_loop,)


@app.cell(hide_code=True)
def _(Gallery, four_loop, mo):
    mo.stop(four_loop is None or four_loop["missing"] is not None)
    _pole = next(o for o in four_loop["orders"] if o["order"] == -1)
    _finite = next((o for o in four_loop["orders"] if o["order"] == 0), None)
    mo.vstack(
        [
            mo.md(
                "**ε⁻¹, exact** — no four-loop constant survives "
                f"({', '.join(_finite['named']) if _finite and _finite['named'] else 'none'} "
                "appear only at ε⁰):"
            ),
            Gallery.formula(mo, _pole["expression"]),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, bs, mo):
    five_loop_examples = mo.ui.multiselect(
        {
            bs.EXAMPLES[k]["title"]
            + (
                " (heavy: " + bs.EXAMPLES[k]["heavy"] + ")"
                if "heavy" in bs.EXAMPLES[k]
                else ""
            ): k
            for k in bs.EXAMPLES
            if bs.EXAMPLES[k]["loops"] == 5
        },
        value=[
            bs.EXAMPLES[k]["title"]
            + (
                " (heavy: " + bs.EXAMPLES[k]["heavy"] + ")"
                if "heavy" in bs.EXAMPLES[k]
                else ""
            )
            for k in bs.default_examples()
        ],
        label="Five-loop examples to evaluate (the heavy ones run in this page's kernel)",
        full_width=True,
    )
    five_loop_examples if Gallery.five_loop_load_enabled() else mo.md("")
    return (five_loop_examples,)


@app.cell(hide_code=True)
def _(Gallery, bs, data_basis, example_digits, five_loop_examples, mo, run_examples):
    # No mo.stop here: the cell below reads five_loop ({} = not evaluated).
    five_loop = {}
    _blocks = []
    if Gallery.five_loop_load_enabled() and (run_examples.value or bs.autorun()):
        for _key in five_loop_examples.value:
            _needs = bs.EXAMPLES[_key].get("needs")
            if _needs and data_basis != _needs:
                _blocks.append(mo.md(f"### {bs.EXAMPLES[_key]['title']}"))
                _blocks.append(
                    Gallery.note(
                        mo,
                        f"Needs {_needs} data",
                        f"The loaded data are {data_basis or 'not a 5L tree'}. "
                        "This example shows a BASIS-v3 property, so the page does not "
                        "evaluate it on these data: on BASIS-v2, P108 is itself a "
                        "master with numerical poles, and T1-8175 is no master and "
                        "reduces only slowly (a static export on v0.2.2 ran out of "
                        "45 GB after 54 minutes). " + bs.EXAMPLES[_key]["point"],
                        pending=True,
                    )
                )
                continue
            try:
                five_loop[_key] = bs.evaluate_example(
                    Gallery, _key, precision=int(example_digits.value)
                )
            except (ValueError, RuntimeError, OSError) as _error:
                five_loop[_key] = {"error": str(_error)}
            _result = five_loop[_key]
            _blocks.append(mo.md(f"### {bs.EXAMPLES[_key]['title']}"))
            _blocks.append(mo.md(bs.EXAMPLES[_key]["point"]))
            if "error" in _result:
                _blocks.append(
                    mo.callout(f"Vakint error: {_result['error']}", kind="warn")
                )
            elif _result["missing"]:
                _blocks.append(
                    mo.callout(
                        "Vakint fails closed: "
                        + Gallery.missing_orders_text(_result["missing"])
                        + ".",
                        kind="warn",
                    )
                )
            else:
                _blocks.append(
                    mo.md(
                        f"Data {_result['data']} · evaluation "
                        f"{Gallery.duration(_result['seconds']['evaluate'])} "
                        f"(+ {Gallery.duration(_result['seconds']['load'])} "
                        "engine and data)"
                    )
                )
                _blocks.append(Gallery.table(mo, bs.example_table_rows(_result)))
    mo.vstack(_blocks)
    return (five_loop,)


@app.cell(hide_code=True)
def _(Gallery, five_loop, mo):
    _shown = []
    for _key, _label in (("5L-P108", "P108"), ("5L-T1-8175", "T1-8175")):
        _result = five_loop.get(_key) or {}
        _order = next(
            (o for o in _result.get("orders") or [] if o["order"] == -2), None
        )
        if _order is not None and _order["expression"] is not None:
            _shown.append(
                mo.md(
                    f"**{_label}, ε⁻² ({_order['class']} = {_order['exact'] or 'numerical'}):**"
                )
            )
            _shown.append(Gallery.formula(mo, _order["expression"]))
    mo.stop(not _shown)
    mo.vstack(_shown)
    return


@app.cell(hide_code=True)
def _(Gallery, bs, mo, overview):
    mo.stop(overview is None)
    mo.vstack(
        [
            Gallery.section(
                mo,
                "04 / named constants",
                "The constants in the exact results",
                "Built-in constants evaluate to any precision; the finite-digit ones "
                "are stored with the digits shown. The traditional forms are "
                "verified to the full stored digits but not yet adopted in the "
                "data, which still uses the names in the first column.",
            ),
            Gallery.table(mo, bs.constants_rows(overview)),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, bs, mo, overview):
    _method = bs.method_paragraphs(overview) if overview else []
    mo.vstack(
        [
            Gallery.section(
                mo,
                "05 / method",
                "How the exact poles were derived and checked",
            ),
            mo.md(
                r"""
**The pole of a vacuum integral is simpler than its finite part.** The
R-operation of renormalization theory removes the ultraviolet divergences of a
Feynman integral subgraph by subgraph. Applied to an L-loop vacuum integral it
writes every pole as products of lower-loop integrals plus one overall
counterterm. So the poles at L loops need only lower-loop finite parts, which
are known exactly, plus a single number per pole order: the overall
counterterm.

**Part A from the R-operation, the counterterm from R\*.** For each five-loop
master, the subgraph part of the poles ("part A") is assembled exactly from
the one- to four-loop results. The overall counterterm was first identified
from the high-precision numerical values (an integer-relation fit with a
safety margin of tens of digits), then derived independently with R\*, a
variant of the R-operation that also handles infrared divergences and gives
the counterterm directly from simpler propagator integrals. For the masters
with a negative superficial degree of divergence the counterterm is zero and
the poles are exact outright. Two coefficients of master P6 follow exactly
from the requirement that the former masters reduce consistently.

**A numerical gate.** Every adopted form was evaluated to high precision and
compared with the stored numerical value of the same coefficient; it enters
the data only if they agree to the digits that value claims (published digits
are the claim minus one).

**Two independent reviews.** One reviewer re-evaluated all 84 forms and
re-identified the fitted ones blind; the other reviewed the derivation and
re-ran R\* on 85 masters. No coefficient was found wrong, and no pole rests on
a fit alone any more.

**Four loops.** The same argument makes every four-loop pole exact: the
Bessel-moment constants (K5Moment, PR7D2, ...) that appear in four-loop finite
parts must cancel at ε⁻¹, and they do, on every reduction record checked
(more than 340 000 cancelling products). The first example above shows it live.
"""
            ),
            mo.accordion(
                {
                    f"Benchmark method notes: {title}": mo.md(text)
                    for title, text in _method
                }
            )
            if _method
            else mo.md(""),
            Gallery.footer(mo),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
