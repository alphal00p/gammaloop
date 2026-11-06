"""Explore the five-loop zero-momentum QCD fixture: input, plan and a recorded run.

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
        "the Ordered reduction plan and a recorded run with its exact retained "
        "result. This page never starts a five-loop reduction by itself.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    # Catalogue counts come from the configured (or installed) data, never literals.
    _catalogue = Gallery.catalogue_summary() or {}
    _five = _catalogue.get("loops", {}).get(5)
    _masters = (
        [
            (
                _five["masters"],
                (
                    f"5L catalogue masters ({_five['computed']} with numerical "
                    f"values) · data {_catalogue.get('version')}"
                ),
            )
        ]
        if _five
        else []
    )
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    ("5", "vacuum loops"),
                    ("14", "total denominator powers"),
                    ("5", "original momentum rank"),
                    ("9", "zero-momentum input: A − r_mom"),
                    *_masters,
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
            Gallery.table(
                mo,
                [
                    {
                        "contribution": "Complete zero-momentum contraction",
                        "A": 14,
                        "r_mom": 5,
                        "A − r_mom": 9,
                    },
                ],
                numeric=("A", "r_mom", "A − r_mom"),
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
            "A 5L pack already installed in the download cache is used too.",
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
                Gallery.table(
                    mo,
                    five_loop_plan["supports"],
                    numeric=(
                        "owner support (lines)",
                        "owners",
                        "rule bundles (MiB)",
                    ),
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
    if Gallery.five_loop_runs_enabled():
        run_benchmark = mo.ui.run_button(
            label="Run the five-loop zero-momentum reduction", kind="warn"
        )
        _note = Gallery.note(
            mo,
            "Enabled on this server · substantial cost",
            f"The run reduces the verified zero-momentum input with the Ordered "
            f"engine on at most {Gallery.max_cores()} cores, at low priority, "
            f"with a {Gallery.max_seconds():g} s deadline. Cancel stops every "
            "process of the run. The quadratic fixture has no Run action.",
            pending=True,
        )
    else:
        # No disabled button: the note says why and what to look at instead.
        run_benchmark = None
        _note = Gallery.note(
            mo,
            "Explore-only on this server",
            "A five-loop reduction takes minutes on several cores and about 15 GB "
            "of memory, so this server does not run it. The recorded run below "
            "shows its monitor, timings and exact retained result. Operators "
            "enable one bounded, cancellable zero-momentum run with "
            "VAKINT_GALLERY_ALLOW_5L=1 and a configured runner.",
            pending=True,
        )
    mo.vstack([_note, run_benchmark] if run_benchmark is not None else [_note])
    return (run_benchmark,)


@app.cell
def _(Gallery, mo, run_benchmark):
    benchmark_job = None
    if run_benchmark is not None and run_benchmark.value:
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
            Gallery.job_monitor_card(
                mo, interactive_status, benchmark_job, title="Five-loop run"
            ),
            Gallery.note(
                mo,
                "A separate process group",
                "Only one gallery job runs at a time; Cancel stops all of its processes."
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
            mo.accordion(
                {"Raw job status": mo.json(Gallery.display_receipt(interactive_status))}
            ),
        ]
    )
    return (interactive_status,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    five_loop_fixture = Gallery.five_loop_fixture()
    mo.stop(
        five_loop_fixture is None,
        Gallery.note(
            mo,
            "No recorded five-loop run in this evidence directory",
            "data/five-loop-q0 is listed in data/manifest.json of the gallery "
            "sources; a VAKINT_GALLERY_EVIDENCE directory without it shows nothing here.",
            pending=True,
        ),
    )
    _report = five_loop_fixture["report"]
    _closure = _report.get("catalogue_closure") or {}
    _memory = _report.get("memory_after_retained_reduction") or {}
    _data = _report.get("data") or {}
    _execution = _report.get("execution") or {}
    _cache = _report.get("ordered_step_cache") or {}
    _provenance = _report.get("gallery_provenance") or {}
    mo.vstack(
        [
            Gallery.section(
                mo,
                "06 / recorded run",
                "The zero-momentum reduction, as it ran",
                "One complete release-build run of the saved input in the current "
                "131-master basis: its native monitor at completion, timings and "
                "exact retained result. Nothing is recomputed here.",
            ),
            Gallery.metrics(
                mo,
                [
                    (
                        f"{_closure.get('master_count')} · "
                        + ("closed" if _closure.get("closed") else "NOT closed"),
                        "catalogue masters in the result",
                    ),
                    (
                        f"{_report.get('retained_reduction_seconds', 0):,.0f} s",
                        "retained reduction (whole public call)",
                    ),
                    (
                        f"{_report.get('loading_seconds', 0):,.0f} s",
                        "external data and saved input loading",
                    ),
                    (
                        Gallery.gib(_memory.get("process_peak_rss_bytes")),
                        "peak memory (whole process)",
                    ),
                    *Gallery.execution_metrics(
                        f"{_execution.get('ordered_weight_arithmetic')} · "
                        f"{_execution.get('n_cores')} cores"
                    ),
                ],
            ),
            Gallery.monitor_card(
                mo,
                five_loop_fixture["model"],
                title="Recorded monitor at completion",
            ),
            Gallery.note(
                mo,
                "What this run establishes",
                f"Data {_data.get('version')} (manifest sha256 "
                f"{str(_data.get('manifest_sha256'))[:16]}…), {_report.get('guards', {}).get('count', 0):,} "
                "exported nonzero guards, "
                f"{_cache.get('misses', 0):,} Ordered step-cache misses and "
                f"{_cache.get('hits', 0):,} hits (a cold in-memory cache). The "
                "masters stay formal symbols; the result is the exact retained "
                "linear combination, identical for every core count and Ordered "
                "weight arithmetic. Recorded "
                f"{_provenance.get('recorded', 'by the integration gate')}.",
            ),
            Gallery.evidence_links(mo, "five-loop-q0"),
        ]
    )
    return (five_loop_fixture,)


@app.cell
def _(five_loop_fixture, mo):
    inspect_saved = mo.ui.run_button(
        label="Show the exact retained result (offline, no reduction)",
        kind="success",
    )
    mo.stop(five_loop_fixture["retained"] is None)
    inspect_saved
    return (inspect_saved,)


@app.cell(hide_code=True)
def _(Gallery, five_loop_fixture, inspect_saved, mo):
    mo.stop(not inspect_saved.value)
    from symbolica import E as _E

    _path, _digest = five_loop_fixture["retained"]
    saved_result = _E(Gallery.checked_bytes(_path, _digest).decode())
    # The live pager is this cell's own output: page through all terms and
    # open the numbered omission markers.
    mo.vstack(
        [
            Gallery.formula(mo, saved_result),
            Gallery.legend(
                mo,
                "RustRedMaster(rustred_family_…, powers) are the formal catalogue "
                "masters; GC_11 is the UFO strong coupling, mursq = μ_R², mUV the "
                "auxiliary UV mass and ε the dimensional regulator (d = 4 − 2ε). "
                f"sha256 {_digest} verified before parsing.",
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
