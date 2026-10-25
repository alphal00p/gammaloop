"""The complete quadratic QCD fixture, its provenance, and measured reductions."""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · Five-loop UV reduction")


@app.cell
def _():
    import os

    import marimo as mo
    from symbolica.community.hepkit import FeynmanDiagram, Model

    from shared_gallery import Gallery

    return FeynmanDiagram, Gallery, Model, mo, os


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "five-loop",
        "A five-loop UV coefficient",
        "One closed massless quark trace, four gluon rungs, and the complete "
        "quadratic transverse Taylor coefficient. Follow the scientific input "
        "through the authenticated reduction campaign.",
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
                    ("14 → 16", "total denominator powers after Taylor expansion"),
                    ("10", "unlowered quadratic coefficient: A − r_mom"),
                    ("640", "declared 5L campaign masters; uncomputed"),
                ],
            ),
            Gallery.section(
                mo, "01 / scientific input", "Keep the full quadratic fixture"
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
        $$\frac{\delta_{ab}}{8(d-1)q^2}
        \left(g_{\mu\nu}-\frac{q_\mu q_\nu}{q^2}\right).$$

        Expand **numerator and propagator product together** through order
        $q^2$, then remove the auxiliary external direction by exact tensor
        reduction. Insert $M_{\mathrm{UV}}$ only in the vacuum denominators;
        the quark slash numerators stay massless. Native graph factors are
        included once.
        """),
            mo.ui.table(
                [
                    {
                        "contribution": "Quadratic · numerator order 2",
                        "A": 14,
                        "r_mom": 4,
                        "A − r_mom": 10,
                    },
                    {
                        "contribution": "Quadratic · mixed expansion",
                        "A": 15,
                        "r_mom": 5,
                        "A − r_mom": 10,
                    },
                    {
                        "contribution": "Quadratic · propagator order 2",
                        "A": 16,
                        "r_mom": 6,
                        "A − r_mom": 10,
                    },
                    {
                        "contribution": "Separate zero-momentum diagnostic",
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
                "This is a single graph's transverse quadratic coefficient, not a gauge-complete sum, beta function, or fully R-subtracted observable. Twelve vacuum incidence edges coalesce to nine distinct quadratic denominator forms.",
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
                """)
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
                "Uncached versus compact-cache reduction",
                "Verified scalar-input loading, external data, native reduction and cache attachment are measured separately. Fresh scientific preparation is independent validation.",
            ),
            Gallery.evidence_card(mo, "quadratic", "Complete quadratic fixture"),
            Gallery.evidence_card(
                mo, "zero-momentum", "Zero-momentum metric/color diagnostic"
            ),
            Gallery.evidence_card(
                mo,
                "zero-momentum-running-cache",
                "Optimized zero-momentum run · cumulative cache progress",
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
    benchmark_choice = mo.ui.dropdown(
        {
            "Complete quadratic coefficient": "quadratic",
            "Zero-momentum diagnostic": "zero-momentum",
        },
        value="Complete quadratic coefficient",
        label="Physical projection",
        full_width=True,
    )
    use_cache = mo.ui.checkbox(
        label="Attach the configured prepared cache",
        disabled=not bool(os.environ.get("VAKINT_GALLERY_CACHE")),
    )
    run_benchmark = mo.ui.run_button(
        label="Run reduction of the saved scalar input",
        kind="warn",
        disabled=not Gallery.benchmark_available(),
    )
    mo.vstack(
        [
            benchmark_choice,
            use_cache,
            Gallery.note(
                mo,
                "May take substantial time",
                "This invokes the compiled Rust reduction with the complete verified scalar input, preserving all contributions and projector conditions. It does not repeat graph generation, gamma traces or Taylor preparation. Reduction can still take minutes or hours; there is no default time cutoff. A single interactive run is not a controlled three-run benchmark.",
            ),
            run_benchmark,
            mo.md(
                "The server has no fully configured runner and prepared scalar-input collection."
                if not Gallery.benchmark_available()
                else "The runner uses the configured verified scalar-input collection. Run logs are retained in its development workspace; fresh scientific generation is validated separately."
            ),
        ]
    )
    return benchmark_choice, run_benchmark, use_cache


@app.cell
def _(Gallery, benchmark_choice, mo, run_benchmark, use_cache):
    mo.stop(not run_benchmark.value)
    try:
        benchmark_job = Gallery.start_benchmark(
            benchmark_choice.value, cached=use_cache.value
        )
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
