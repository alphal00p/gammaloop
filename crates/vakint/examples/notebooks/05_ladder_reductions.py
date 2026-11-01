"""Reduce the three- and four-loop QCD ladders in a cancellable background job.

Interactive: choose a case and press "Reduce in the background"; Cancel stops it.
The page keeps showing its last run (``?job=<id>`` in the URL survives a
reload); "Recent jobs" reopens any finished run.

Script mode runs one case headlessly and waits for it, for example
``python 05_ladder_reductions.py --case 3L-zero-momentum --n-cores 2``
(optional ``--arithmetic``, ``--max-seconds``, ``--cancel-after SECONDS``).
It prints one JSON line (job, state, retained sha256 and whether it matches the
recorded fixture, closure, seconds) and exits 0 only if the run reproduced the
fixture, or was cancelled as ``--cancel-after`` asked; 1 otherwise, 2 if the
run could not start (it never watches or cancels another session's job), 130
after Ctrl-C (which cancels the job).
"""

# Marimo cell outputs and generated return markers are intentional.
# ruff: noqa: B018, PLR1711
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Vakint · Reduce 3L & 4L")


@app.cell
def _():
    import marimo as mo
    from shared_gallery import Gallery

    return Gallery, mo


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "reduce",
        "Reduce three- and four-loop ladders",
        "The massless-quark QCD ladder self-energy with two or three gluon rungs, "
        "reduced end to end to retained masters with the Ordered engine. Each run "
        "is a separate, bounded process that you can cancel at any time.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    ladder_cases = Gallery.ladder_cases()
    _rows = [
        {
            "case": case["label"],
            "prepared input": f"{case['input']['bytes'] / 1024:,.0f} KB",
            "external data": case["external_data"],
            "recorded data loading (s)": case["recorded"]["data_loading_seconds"],
            "recorded reduction (s)": case["recorded"]["retained_reduction_seconds"],
            "catalogue masters": case["catalogue_masters"],
            "here": "ready"
            if Gallery.ladder_case_available(key)
            else "input not provided",
        }
        for key, case in ladder_cases.items()
    ]
    mo.vstack(
        [
            Gallery.section(
                mo,
                "01 / cases",
                "Prepared scalar inputs, verified before every run",
                "Each input is the rung fixture's prepared numerator and denominators "
                "(graph, Dirac and color algebra, Taylor projection and vacuum "
                "routing already done). Its sha256 is checked before reduction.",
            ),
            mo.ui.table(_rows, selection=None, pagination=False, show_download=False),
            Gallery.note(
                mo,
                "What a run computes",
                "The retained-master result in MSbar with dot-product notation and "
                "L+1 epsilon orders, as in the rung test. Every Ordered weight "
                "arithmetic and every core count return the identical expression, so "
                "its sha256 is compared with the recorded fixture. Recorded timings are "
                "one release-build run at 4 cores; the 4L cases first load the 4L "
                "external-data pack (about 140 MB, installed once).",
            ),
        ]
    )
    return (ladder_cases,)


@app.cell
def _(Gallery, ladder_cases, mo):
    _labels = {
        case["label"]: key
        for key, case in ladder_cases.items()
        if Gallery.ladder_case_available(key)
    }
    _default = next(
        case["label"]
        for key, case in ladder_cases.items()
        if case["default"] and key in _labels.values()
    )
    _cores = Gallery.max_cores()
    _deadline = Gallery.max_seconds()
    run_form = mo.md(
        "{case}\n\n{arithmetic}\n\n{n_cores}\n\n{max_seconds}\n\n{ordered_cache}"
    ).batch(
        case=mo.ui.dropdown(_labels, value=_default, label="Case", full_width=True),
        arithmetic=mo.ui.dropdown(
            list(Gallery.arithmetics),
            value="modular_polynomial",
            label="Ordered weight arithmetic",
        ),
        n_cores=mo.ui.slider(
            start=1,
            stop=_cores,
            value=min(2, _cores),
            label=f"n_cores (at most {_cores} on this server)",
            show_value=True,
        ),
        max_seconds=mo.ui.number(
            start=10,
            stop=_deadline,
            value=min(600, _deadline),
            step=10,
            label="Deadline (seconds)",
        ),
        ordered_cache=mo.ui.checkbox(
            value=True, label="Persistent Ordered step cache (warm reruns)"
        ),
    )
    run_button = mo.ui.run_button(label="Reduce in the background", kind="success")
    mo.vstack(
        [
            Gallery.section(mo, "02 / run", "Start a bounded background reduction"),
            run_form,
            run_button,
            Gallery.note(
                mo,
                "Bounded and cancellable",
                f"A run is its own process group on at most {_cores} CPUs, at low "
                "priority and with a wall-clock deadline. One run at a time per "
                "server. Cancel stops every process of the run. The page keeps "
                "showing your last run, also after a reload; a fresh page "
                "re-attaches to a run that is still going, and Recent jobs below "
                "reopens any finished one.",
            ),
        ]
    )
    return run_button, run_form


@app.cell
def _(Gallery, mo):
    get_finished_job, set_finished_job = mo.state(None)
    # The run this session shows: its last one, or ?job=<id> from the URL.
    get_shown_job, set_shown_job = mo.state(
        Gallery.ladder_job(mo.query_params().get("job"))
    )
    return get_finished_job, get_shown_job, set_finished_job, set_shown_job


@app.cell
def _(Gallery, mo, run_button, run_form):
    # Script mode runs exactly the case named on the command line; the report
    # cell below prints its outcome and sets the exit status.
    headless = mo.app_meta().mode == "script"
    run_request = started_job = start_error = None
    try:
        if headless:
            run_request = Gallery.headless_request(mo.cli_args(), run_form.value)
        elif run_button.value:
            run_request = {**run_form.value, "cancel_after": None}
        if run_request is not None:
            started_job = Gallery.start_reduction(
                run_request["case"],
                n_cores=int(run_request["n_cores"]),
                arithmetic=run_request["arithmetic"],
                max_seconds=float(run_request["max_seconds"]),
                ordered_cache=bool(run_request["ordered_cache"]),
            )
    except (OSError, ValueError, RuntimeError) as _error:
        start_error = str(_error)
    if headless and started_job is None:
        # Refused (unknown case, missing input, another run going): report it
        # and exit; never watch or cancel another session's job.
        Gallery.headless_exit(run_request, start_error)
    mo.callout(start_error, kind="warn") if start_error else mo.md("")
    return headless, run_request, started_job


@app.cell
def _(Gallery, get_shown_job, headless, mo, set_shown_job, started_job):
    if headless:
        # A headless run watches, reports and cancels only its own job.
        watched_job = started_job
    else:
        # Keep showing this session's last run when a control changes, else
        # re-attach to a run still going (a reload, another session).
        watched_job = (
            started_job or get_shown_job() or Gallery.active_job(kind="ladder")
        )
        if watched_job is not None and watched_job != get_shown_job():
            set_shown_job(watched_job)
            mo.query_params()["job"] = watched_job.rsplit("/", 1)[-1]
    return (watched_job,)


@app.cell
def _(Gallery, get_finished_job, mo, watched_job):
    mo.stop(
        watched_job is None,
        mo.md("_No reduction is running. Choose a case and start one above._"),
    )
    if (
        get_finished_job() == watched_job
        or Gallery.job_status(watched_job).get("state") in Gallery.final_states
    ):
        # The run is over: no Cancel, and no more status polling.
        cancel_button = job_refresh = None
    else:
        cancel_button = mo.ui.button(
            label="Cancel this run",
            kind="danger",
            on_click=lambda _: Gallery.cancel_job(watched_job),
        )
        job_refresh = mo.ui.refresh(
            options=[1, 2, 5], default_interval=1, label="Status refresh"
        )
    mo.hstack(
        [element for element in (cancel_button, job_refresh) if element is not None],
        justify="start",
    )
    return cancel_button, job_refresh


@app.cell(hide_code=True)
def _(
    Gallery,
    cancel_button,
    get_finished_job,
    headless,
    job_refresh,
    mo,
    run_request,
    set_finished_job,
    watched_job,
):
    job_refresh is not None and job_refresh.value
    if headless:
        # Ctrl-C cancels the job rather than leaving it to its deadline.
        job_status = Gallery.wait_job(
            watched_job,
            cancel_after=run_request["cancel_after"],
            cancel_on_interrupt=True,
        )
    else:
        job_status = Gallery.job_status(watched_job)
    _state = job_status.get("state")
    _final = _state in Gallery.final_states
    if _final and get_finished_job() != watched_job:
        set_finished_job(watched_job)
    _report = (
        (cancel_button.value if cancel_button is not None else None)
        or job_status.get("cancel_report")
        or Gallery.cancel_report(watched_job)
    )
    _parts = [
        Gallery.section(mo, "03 / progress", "The background run"),
        mo.ui.table(
            [
                {"": label, "value": value}
                for label, value in Gallery.progress_rows(job_status)
            ],
            selection=None,
            pagination=False,
            show_download=False,
        ),
    ]
    if _report:
        _survivors = _report.get("survivors") or []
        _parts.append(
            Gallery.note(
                mo,
                "Cancellation",
                f"Signalled {len(_report.get('signalled') or [])} process(es); "
                + (
                    f"still alive: {_survivors}"
                    if _survivors
                    else "no process of the run survives."
                ),
                pending=bool(_survivors),
            )
        )
    _parts.append(
        mo.accordion({"Raw job status": mo.json(Gallery.display_receipt(job_status))})
    )
    mo.vstack(_parts)
    return (job_status,)


@app.cell(hide_code=True)
def _(Gallery, headless, job_status, run_request, watched_job):
    # Script mode: one JSON line on stdout; a non-zero exit status unless the
    # run reproduced its recorded fixture (or was cancelled as asked).
    if headless:
        Gallery.headless_report(watched_job, job_status, run_request)
    return


@app.cell(hide_code=True)
def _(Gallery, job_status, ladder_cases, mo, watched_job):
    mo.stop(job_status.get("state") not in Gallery.final_states)
    ladder_result = Gallery.job_result(watched_job)
    _status = job_status
    if (
        _status.get("state") != "complete"
        or not ladder_result
        or ladder_result.get("schema") != "vakint.gallery.ladder-run.v1"
        or "error" in ladder_result
    ):
        _message = (ladder_result or {}).get("error") or _status.get("message") or ""
        mo.stop(
            True,
            mo.vstack(
                [
                    Gallery.note(
                        mo,
                        f"Run {str(_status.get('state')).replace('_', ' ')}",
                        str(_message)[-3000:],
                        pending=True,
                    )
                ]
            ),
        )
    _case = ladder_cases[ladder_result["case"]]
    _retained = ladder_result["retained"]
    _closure = ladder_result["catalogue_closure"]
    _seconds = ladder_result["seconds"]
    _metrics = [
        (
            "identical" if _retained["matches_recorded"] else "DIFFERENT",
            "retained sha256 vs recorded fixture",
        ),
        (
            f"{_closure['master_count']} · closed"
            if _closure["closed"]
            else f"{_closure['master_count']} · NOT closed",
            "catalogue masters in the result",
        ),
        (f"{_seconds.get('Ordered reduction', 0):.3f} s", "Ordered reduction"),
        (f"{_seconds['total']:.1f} s", "whole job (import, data, reduction)"),
        (ladder_result["execution"], "arithmetic and cores"),
    ]
    _cpu = ladder_result.get("cpu_seconds", {})
    _timings = [
        {"stage": stage, "wall seconds": value, "CPU seconds": _cpu.get(stage)}
        for stage, value in _seconds.items()
    ]
    _reference = [
        {
            "ε power": row["epsilon_power"],
            "real": row["real"],
            "imaginary": row["imaginary"],
        }
        for row in _case["numerical_reference"]["coefficients"]
    ]
    _text = (Gallery.checked_job(watched_job) / _retained["file"]).read_text()
    try:
        from symbolica import E

        _formula = Gallery.native_formula(mo, E(_text))
    except (ImportError, ValueError, TypeError):
        _formula = mo.ui.code_editor(_text[:20000], language="python", disabled=True)
    _cache = ladder_result.get("ordered_cache")
    mo.vstack(
        [
            Gallery.section(
                mo,
                "04 / result",
                _case["label"],
                f"sha256 {_retained['sha256']} · data: "
                + Gallery.compact_json(ladder_result.get("data")),
            ),
            Gallery.metrics(mo, _metrics),
            mo.accordion(
                {
                    "Retained-master result": _formula,
                    "Stage timings (wall and CPU seconds)": mo.ui.table(
                        _timings, selection=None, pagination=False, show_download=False
                    ),
                    "Recorded numerical reference (not recomputed here)": mo.vstack(
                        [
                            mo.md(
                                "High-precision Laurent coefficients of this input at "
                                "mUV = μ = G = 1 (MSbar), from the rung fixture."
                            ),
                            mo.ui.table(
                                _reference,
                                selection=None,
                                pagination=False,
                                show_download=False,
                            ),
                        ]
                    ),
                    "Ordered step cache": mo.json(_cache)
                    if _cache
                    else mo.md("Not used for this run."),
                    "Conventions": mo.md(_case["conventions"]),
                }
            ),
        ]
    )
    return (ladder_result,)


@app.cell(hide_code=True)
def _(Gallery, get_finished_job, mo, watched_job):
    # Listed again when a run starts or ends; a link reopens a run's result.
    get_finished_job()
    _rows = []
    for _work in Gallery.jobs()[:8]:
        try:
            _status = Gallery.job_status(_work)
        except (OSError, ValueError):
            continue
        _name = _work.rsplit("/", 1)[-1]
        _rows.append(
            {
                "job": _name,
                "state": _status.get("state"),
                "elapsed (s)": round(_status.get("elapsed_seconds") or 0, 1),
                "result": "shown above"
                if _work == watched_job
                else Gallery.job_link(mo, _name)
                if Gallery.ladder_job(_name)
                else "",
            }
        )
    mo.vstack(
        [
            mo.accordion(
                {
                    "Recent jobs on this server": mo.ui.table(
                        _rows, selection=None, pagination=False, show_download=False
                    )
                    if _rows
                    else mo.md("None yet.")
                }
            ),
            Gallery.footer(mo),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
