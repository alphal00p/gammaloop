"""Explore the five-loop zero-momentum QCD fixture: input, plan and saved result.

Explore-only: opening or using this page never starts a five-loop reduction.
A server operator can enable one bounded, cancellable zero-momentum run with
VAKINT_GALLERY_ALLOW_5L=1 and a configured runner (see gallery.typ).
"""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · Five-loop: explore")


@app.cell
def _():
    import marimo as mo
    from shared_gallery import Gallery
    from symbolica.community.hepkit import FeynmanDiagram, Model

    return FeynmanDiagram, Gallery, Model, mo


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "five-loop",
        "Explore a five-loop zero-momentum reduction",
        "One closed massless quark trace and four gluon rungs. Inspect the input, "
        "the Ordered reduction plan and the saved exact result. This page never "
        "starts a five-loop reduction by itself.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    ("5", "vacuum loops"),
                    ("14", "total denominator powers"),
                    ("5", "original momentum rank"),
                    ("9", "zero-momentum input: A − r_mom"),
                    ("640", "declared 5L campaign masters; uncomputed"),
                ],
            ),
            Gallery.section(
                mo, "01 / scientific input", "Keep the complete zero-momentum input"
            ),
        ]
    )
    return


@app.cell
def _(FeynmanDiagram, Gallery, Model):
    qcd_model = Model.qcd()
    ladder = FeynmanDiagram.from_dot(qcd_model, Gallery.ladder_dot(5))
    ladder
    return ladder, qcd_model


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            mo.md(r"""
        $$g_{\mu\nu}\,\delta_{ab}\,
        \Pi^{ab,\mu\nu}(q)\big|_{q=0}.$$

        Set the external momentum to zero and retain **both current vertices**
        and all contributions. The metric/color contraction is unaveraged.
        Insert $M_{\mathrm{UV}}$ only in the vacuum denominators; the quark
        slash numerators stay massless. Native diagram factors are included once.
        """),
            mo.ui.table(
                [
                    {
                        "contribution": "Complete zero-momentum contraction",
                        "A": 14,
                        "r_mom": 5,
                        "A − r_mom": 9,
                    },
                ],
                selection=None,
                pagination=False,
                show_download=False,
            ),
            Gallery.note(
                mo,
                "One graph, precisely defined",
                "This is one graph's zero-momentum vacuum contribution, not a gauge-complete sum, beta function, or fully R-subtracted observable. The displayed graph supplies its geometry; the saved scalar input preserves the specialized numerator, routing and native factors. Twelve vacuum incidence edges coalesce to nine distinct quadratic denominator forms.",
            ),
            mo.accordion(
                {
                    "Momentum rank, saved indices and UV degree": mo.md(r"""
                The table uses the original momentum rank $r_{\mathrm{mom}}$.
                Saved predicates use signed denominator indices, their index rank
                $r_{\mathrm{idx}}$, and $D_{\mathrm{key}}=\sum_i a_i$.

                At fixed $M_{\mathrm{UV}}$, the five-loop UV degrees are
                $5d+2r_{\mathrm{mom}}-2A$ before lowering and
                $5d-2D_{\mathrm{key}}$ for an individual scalar term.
                A coefficient $(M_{\mathrm{UV}}^2)^j$ balances engineering
                mass dimension through $D_{\mathrm{key}}-j$; it does not
                change that scalar term's loop-momentum UV degree.
                """),
                    "Deferred quadratic fixture · no Run action": mo.md(r"""
                The complete quadratic Taylor coefficient and its prepared input
                are retained for a later decision. Its transverse contraction is
                $\delta_{ab}(g_{\mu\nu}-q_\mu q_\nu/q^2)/[8(d-1)q^2]$;
                the original groups have $(A,r_{\mathrm{mom}})=(14,4),(15,5),(16,6)$.
                Quadratic reduction and retesting are **deferred**. No completed
                quadratic reduction, cache comparison or Run action is offered here.
                """),
                }
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
                "02 / campaign",
                "Recorded predicates are the authority",
                "The recorded regions have D = 9 or 10 and native index rank up to six. Membership follows each saved signed-key predicate.",
            ),
            Gallery.evidence_card(
                mo,
                "campaign",
                "Campaign identity and mapped-input census",
                optional=True,
            ),
            Gallery.note(
                mo,
                "Reduction outside certification",
                "Vakint still attempts individually justified reductions outside recorded coverage. After exact combination it warns about surviving non-catalogue integrals. Retain mode preserves those terms; master substitution reports unavailable data explicitly.",
            ),
            Gallery.evidence_card(
                mo,
                "warning-controls",
                "Recorded coverage advisory and completed reduction",
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
                "03 / measurements",
                "Ordered reduction of the saved input",
                "Verified scalar-input loading, external data and native ordered reduction are measured separately. Fresh scientific preparation is independent validation.",
            ),
            Gallery.evidence_card(
                mo,
                "zero-momentum",
                "Complete zero-momentum physical reduction",
                optional=True,
            ),
            mo.md(
                "Final comparisons use three fresh processes per configuration with the same verified scalar input, medians and ranges, fixed workers/affinity, and exact result/guard checks. Original graph and numerator preparation are reported separately. Warm in-process reuse is reported separately."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    five_loop_plan = Gallery.five_loop_plan()
    if five_loop_plan is None:
        _body = Gallery.note(
            mo,
            "Reduction plan not available here",
            "Set VAKINT_GALLERY_ARTIFACTS to a vakint-external-data checkout to "
            "list the plan. Only its JSON index, 5L/ibps/selection.json (a few "
            "MB), is read; no IBP payload is decoded and nothing is downloaded. "
            "A 5L pack already installed in the gallery's data cache is used too.",
            pending=True,
        )
    else:
        _body = mo.vstack(
            [
                Gallery.metrics(
                    mo,
                    [
                        (five_loop_plan["owners"], "Apply owners (one phase each)"),
                        (
                            f"{five_loop_plan['route_records']:,}",
                            "saved Route records",
                        ),
                        (
                            f"{five_loop_plan['bundle_bytes'] / 1024**2:,.0f} MiB",
                            "owner rule bundles",
                        ),
                        (five_loop_plan["version"], "external-data version"),
                    ],
                ),
                mo.ui.table(
                    five_loop_plan["supports"],
                    selection=None,
                    pagination=False,
                    show_download=False,
                ),
            ]
        )
    mo.vstack(
        [
            Gallery.section(
                mo,
                "04 / plan",
                "The Ordered reduction plan, without reducing",
                "Owners are sectors whose IBP rules eliminate integrals; the Ordered "
                "scheduler reduces each owner's support exactly once, from the "
                "largest support down. Read from 5L/ibps/selection.json.",
            ),
            _body,
        ]
    )
    return (five_loop_plan,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "05 / reproduce",
        "Reduce the verified saved input (server opt-in)",
        "Opening this page never starts the physical reduction or builds a cache.",
    )
    return


@app.cell
def _(Gallery, mo):
    _enabled = Gallery.five_loop_runs_enabled()
    run_benchmark = mo.ui.run_button(
        label="Run the five-loop zero-momentum reduction",
        kind="warn",
        disabled=not _enabled,
    )
    mo.vstack(
        [
            Gallery.note(
                mo,
                "Explore-only unless the operator enables it",
                "A five-loop reduction takes minutes on many cores and over 10 GB "
                "of memory. This server does not allow it: explore the plan above "
                "and the saved result below instead. Operators enable one bounded, "
                "cancellable zero-momentum run with VAKINT_GALLERY_ALLOW_5L=1 and a "
                "configured runner.",
                pending=True,
            )
            if not _enabled
            else Gallery.note(
                mo,
                "Enabled on this server · substantial cost",
                f"The run reduces the verified zero-momentum input with the Ordered "
                f"engine on at most {Gallery.max_cores()} cores, at low priority, "
                f"with a {Gallery.max_seconds():g} s deadline. Cancel stops every "
                "process of the run. The quadratic fixture has no Run action.",
                pending=True,
            ),
            run_benchmark,
        ]
    )
    return (run_benchmark,)


@app.cell
def _(Gallery, mo, run_benchmark):
    benchmark_job = None
    if run_benchmark.value:
        try:
            benchmark_job = Gallery.start_benchmark("zero-momentum")
        except (OSError, ValueError, RuntimeError) as _error:
            mo.stop(True, mo.callout(str(_error), kind="warn"))
    # A reloaded page re-attaches to a five-loop run that is still going.
    benchmark_job = benchmark_job or Gallery.active_job(kind="five-loop")
    mo.stop(benchmark_job is None)
    return (benchmark_job,)


@app.cell
def _(Gallery, benchmark_job, mo):
    cancel_benchmark = mo.ui.button(
        label="Cancel the five-loop run",
        kind="danger",
        on_click=lambda _: Gallery.cancel_job(benchmark_job),
    )
    job_refresh = mo.ui.refresh(
        options=[2, 5, 15], default_interval=2, label="Background job status"
    )
    mo.hstack([cancel_benchmark, job_refresh], justify="start")
    return cancel_benchmark, job_refresh


@app.cell(hide_code=True)
def _(Gallery, benchmark_job, cancel_benchmark, job_refresh, mo):
    job_refresh.value
    interactive_status = Gallery.benchmark_status(benchmark_job)
    _report = cancel_benchmark.value or Gallery.cancel_report(benchmark_job)
    mo.vstack(
        [
            Gallery.note(
                mo,
                f"Five-loop run {str(interactive_status['state']).replace('_', ' ')}",
                "The run is a separate process group; only one gallery job runs at a time."
                + (
                    f" Cancel signalled {len(_report['signalled'])} process(es); "
                    + (
                        "none survives."
                        if not _report["survivors"]
                        else f"still alive: {_report['survivors']}"
                    )
                    if _report
                    else ""
                ),
                pending=interactive_status["state"] != "complete",
            ),
            mo.ui.table(
                [
                    {"": label, "value": value}
                    for label, value in Gallery.progress_rows(interactive_status)
                ],
                selection=None,
                pagination=False,
                show_download=False,
            ),
            mo.accordion(
                {"Raw job status": mo.json(Gallery.display_receipt(interactive_status))}
            ),
        ]
    )
    return (interactive_status,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    mo.vstack(
        [
            Gallery.section(
                mo,
                "06 / saved coefficients",
                "Load the result without repeating its reduction",
                "The exact retained Atom and 300-digit scalar coefficients are available offline in the companion example.",
            ),
            Gallery.metrics(
                mo,
                [
                    ("55", "individual master identities"),
                    ("58", "formal master products"),
                    ("449", "scalar Laurent coefficients"),
                    ("300", "validated exported decimal digits"),
                ],
            ),
            Gallery.note(
                mo,
                "Masters remain formal",
                "The scalar coefficients use unit mass and renormalization scale, retain GC_11^10, and include epsilon^5 with a coefficientwise O(epsilon^6) remainder. Five-loop master values are not evaluated; their own poles prevent interpreting this as a complete amplitude expansion through epsilon^5.",
            ),
            mo.md("""
            From the repository root, verify all saved values without IBPs:

            ```sh
            python crates/vakint/examples/five_loop_q0_coefficients.py verify
            ```

            Native Atom and exact-ratio JSON preserve component precision.
            Numerical `.expr` is a human display.

            [Retained native Atom](https://github.com/alphal00p/gammaloop/raw/refs/heads/vakint_rustred/crates/vakint/examples/data/five_loop_q0/q0-retained.atom)
            · [Coefficient records](https://github.com/alphal00p/gammaloop/blob/vakint_rustred/crates/vakint/examples/data/five_loop_q0/coefficient-table.json)
            · [Conventions and commands](https://github.com/alphal00p/gammaloop/blob/vakint_rustred/crates/vakint/docs/five_loop_q0_coefficients.pdf)
            """),
        ]
    )
    return


@app.cell
def _(mo):
    inspect_saved = mo.ui.run_button(
        label="Inspect the saved result here (offline, no reduction)", kind="success"
    )
    inspect_saved
    return (inspect_saved,)


@app.cell(hide_code=True)
def _(Gallery, inspect_saved, mo):
    mo.stop(not inspect_saved.value)
    saved_result = Gallery.saved_five_loop_result()
    _summary = saved_result.summary()
    _products = saved_result.coefficient_table["products"]
    _rows = []
    for _row in _products[:12]:
        _leading = min(_row["coefficients"], key=lambda number: number["epsilon_power"])
        _rows.append(
            {
                "formal master product": _row["master_product"][:90]
                + ("…" if len(_row["master_product"]) > 90 else ""),
                "leading ε power": _leading["epsilon_power"],
                "leading coefficient (real, 30 digits)": _leading["real"][
                    "display_300_significant_digits"
                ][:32],
            }
        )
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    (_summary["master_symbols"], "individual master identities"),
                    (_summary["master_products"], "formal master products"),
                    (_summary["scalar_coefficients"], "scalar Laurent coefficients"),
                    (
                        f"ε^{_summary['epsilon_minimum']} … ε^{_summary['epsilon_maximum']}",
                        "recorded orders per product",
                    ),
                ],
            ),
            mo.ui.table(_rows, selection=None, pagination=False, show_download=False),
            mo.md(
                f"First 12 of {len(_products)} products; every file was checked "
                "against its recorded sha256. The masters themselves stay formal."
            ),
        ]
    )
    return (saved_result,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.footer(mo)
    return


if __name__ == "__main__":
    app.run()
