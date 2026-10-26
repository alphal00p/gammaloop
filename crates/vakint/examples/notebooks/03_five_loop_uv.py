"""The zero-momentum QCD fixture, its provenance, and cache comparisons."""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · Five-loop UV reduction")


@app.cell
def _():
    import os

    import marimo as mo
    from shared_gallery import Gallery
    from symbolica.community.hepkit import FeynmanDiagram, Model

    return FeynmanDiagram, Gallery, Model, mo, os


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "five-loop",
        "A five-loop zero-momentum reduction",
        "One closed massless quark trace and four gluon rungs. Follow the "
        "unaveraged metric/color contraction through the authenticated campaign "
        "and compare its uncached, light-cache and heavy-cache reductions.",
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
    ladder = FeynmanDiagram.from_dot(
        qcd_model,
        (
            Gallery.directory.parents[1] / "tests/fixtures/qcd_ladder/ladder.dot"
        ).read_text(),
    )
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
                mo, "campaign", "Campaign identity and mapped-input census"
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
                "Uncached versus cumulative-cache reduction",
                "Verified scalar-input loading, external data, native reduction and cache attachment are measured separately. Fresh scientific preparation is independent validation.",
            ),
            Gallery.evidence_card(
                mo, "zero-momentum", "Complete zero-momentum physical reduction"
            ),
            Gallery.evidence_card(
                mo,
                "zero-momentum-running-cache",
                "Optimized zero-momentum run · cumulative cache progress",
            ),
            Gallery.evidence_card(
                mo, "zero-momentum-light-cache", "Official light cache · at most 95 MB"
            ),
            Gallery.evidence_card(
                mo, "zero-momentum-heavy-cache", "Full heavy cache · hosted separately"
            ),
            Gallery.note(
                mo,
                "Two cache tiers, separately measured",
                "The full accumulated heavy cache is preserved without an aggregate size cap and supplied as a compressed local artifact for separate hosting. Only the trimmed light cache, selected by measured rule frequency and usefulness and limited to 95 MB, is published in the official repository and releases. A cached run uses the server's explicitly configured immutable base; light and heavy results remain distinct.",
            ),
            mo.md(
                "Final comparisons use three fresh processes per configuration with the same verified scalar input, medians and ranges, fixed workers/affinity, and exact result/guard checks. Original graph and numerator preparation are reported separately. Warm in-process reuse is reported separately."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.section(
        mo,
        "04 / reproduce",
        "Reduce the verified saved input",
        "Opening this page never starts the physical benchmark or builds a cache.",
    )
    return


@app.cell
def _(Gallery, mo, os):
    use_cache = mo.ui.checkbox(
        label="Attach the configured light or heavy cache",
        disabled=not bool(os.environ.get("VAKINT_GALLERY_CACHE")),
    )
    run_benchmark = mo.ui.run_button(
        label="Run reduction of the saved scalar input",
        kind="warn",
        disabled=not Gallery.benchmark_available(),
    )
    mo.vstack(
        [
            Gallery.note(
                mo,
                "Zero-momentum only",
                "This action reduces the verified zero-momentum input. The quadratic fixture is deferred and has no gallery Run action.",
            ),
            use_cache,
            Gallery.note(
                mo,
                "May take substantial time",
                "This invokes the compiled Rust reduction with the complete verified zero-momentum scalar input, preserving all contributions and contraction conditions. It reuses graph and trace preparation and saves completed native work in a fresh per-run cache delta, including durable prefixes of failed reductions. Reduction can still take minutes or hours; there is no default time cutoff. A single interactive run is not a controlled three-run benchmark.",
            ),
            run_benchmark,
            mo.md(
                "The server has no fully configured runner and prepared scalar-input collection."
                if not Gallery.benchmark_available()
                else "The runner uses the configured verified scalar-input collection. Run logs are retained in its development workspace; fresh scientific generation is validated separately."
            ),
        ]
    )
    return run_benchmark, use_cache


@app.cell
def _(Gallery, mo, run_benchmark, use_cache):
    mo.stop(not run_benchmark.value)
    try:
        benchmark_job = Gallery.start_benchmark("zero-momentum", cached=use_cache.value)
    except (OSError, ValueError, RuntimeError) as _error:
        mo.stop(True, mo.callout(str(_error), kind="warn"))
    return (benchmark_job,)


@app.cell
def _(benchmark_job, mo):
    job_refresh = mo.ui.refresh(
        options=[5, 15, 30], default_interval=5, label="Background job status"
    )
    job_refresh
    return (job_refresh,)


@app.cell(hide_code=True)
def _(Gallery, benchmark_job, job_refresh, mo):
    job_refresh.value
    interactive_status = Gallery.benchmark_status(benchmark_job)
    mo.vstack(
        [
            Gallery.note(
                mo,
                f"Benchmark {interactive_status['state']}",
                "The worker runs independently of the notebook UI. Only one gallery benchmark runs at a time.",
                pending=interactive_status["state"] != "complete",
            ),
            mo.json(interactive_status),
        ]
    )
    return (interactive_status,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.footer(mo)
    return


if __name__ == "__main__":
    app.run()
