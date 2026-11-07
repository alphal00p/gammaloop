"""Reduce the three- and four-loop QCD ladders in a cancellable background job.

Interactive: choose a case (its input graph is drawn natively) and press
"Reduce in the background"; Cancel stops it. The page keeps showing its last
run (``?job=<id>`` in the URL survives a reload); "Recent jobs" reopens any
finished run. In marimo's edit mode a case can also run in the notebook's own
kernel, with the live monitor card; marimo's Interrupt button stops it there.
The last panel, "Input graph -> Laurent series", sets every parameter of the
input, the requested digits and the epsilon orders, and evaluates the whole
input in the kernel to its numerical Laurent coefficients (symbolic result
cached per input; marimo's Interrupt stops it in edit mode).

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


@app.cell(hide_code=True)
def _():
    from pathlib import Path

    import marimo as mo
    from shared_gallery import Gallery

    return Gallery, Path, mo


@app.cell(hide_code=True)
def _(Gallery, mo):
    Gallery.header(
        mo,
        "reduce",
        "Reduce three- and four-loop ladders",
        "The massless-quark QCD ladder self-energy with two or three gluon rungs, "
        "reduced end to end to retained masters with the Ordered engine. Each run "
        "is a separate, bounded process that you can cancel at any time; the last "
        "panel turns the input graph into its numerical Laurent series.",
    )
    return


@app.cell(hide_code=True)
def _(Gallery, mo):
    ladder_cases = Gallery.ladder_cases()
    _rows = [
        {
            "case": case["label"]
            + (
                ""
                if Gallery.ladder_case_available(key)
                else " (input not provided here)"
            ),
            "prepared input": f"{case['input']['bytes'] / 1024:,.0f} KB",
            "external data": case["external_data"],
            "recorded data loading (s)": case["recorded"]["data_loading_seconds"],
            "recorded reduction (s)": case["recorded"]["retained_reduction_seconds"],
            "catalogue masters": case["catalogue_masters"],
        }
        for key, case in ladder_cases.items()
    ]
    _four = (Gallery.catalogue_summary() or {}).get("loops", {}).get(4, {})
    _size = Gallery.size_text(_four.get("bytes"))
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
            Gallery.table(
                mo,
                _rows,
                numeric=(
                    "prepared input",
                    "recorded data loading (s)",
                    "recorded reduction (s)",
                    "catalogue masters",
                ),
            ),
            Gallery.note(
                mo,
                "What a run computes",
                "The retained-master result in MSbar with dot-product notation and "
                "L+1 epsilon orders, as in the rung test. Every Ordered weight "
                "arithmetic and every core count return the identical expression, so "
                "its sha256 is compared with the recorded fixture. Recorded timings are "
                "one release-build run at 4 cores; the 4L cases first load the 4L "
                "external-data pack"
                + (
                    f" ({_size} unpacked, installed once)."
                    if _size
                    else " (installed once)."
                ),
            ),
        ]
    )
    return (ladder_cases,)


@app.cell(hide_code=True)
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
    # One selection drives the graph preview, the runs and the Laurent panel.
    case_choice = mo.ui.dropdown(_labels, value=_default, label="Case", full_width=True)
    mo.vstack(
        [
            Gallery.section(
                mo,
                "02 / input",
                "The input graph",
                "The ladder self-energy of the selected case, drawn natively by "
                "FeynKit from its DOT source with its loop momenta. Both "
                "projections of an L-loop ladder start from this graph.",
            ),
            case_choice,
        ]
    )
    return (case_choice,)


@app.cell(hide_code=True)
def _(Gallery, case_choice, ladder_cases, mo):
    # Redrawn when the selection changes (page 01's FeynKit view).
    _case = ladder_cases[case_choice.value]
    input_diagram = Gallery.ladder_diagram(_case["loop_count"])
    _projection = (
        "the quadratic Taylor coefficient in the external momentum q, "
        "transverse projection and angular average"
        if _case["projection"] == "quadratic"
        else "zero external momentum, unaveraged metric/color contraction"
    )
    mo.vstack(
        [
            Gallery.metrics(
                mo,
                [
                    (input_diagram.loop_count, "independent loops"),
                    (len(input_diagram.internal_edges), "internal particle lines"),
                    (len(input_diagram.external_edges), "external gluons"),
                    (_case["catalogue_masters"], "catalogue masters in its result"),
                ],
            ),
            Gallery.diagram_frame(mo, input_diagram),
            Gallery.legend(
                mo,
                f"Prepared input of this case: {_projection}; massless quark "
                "lines, a common UV mass in the vacuum denominators.",
            ),
            mo.accordion(
                {
                    "Compact source DOT": Gallery.dot_source(
                        mo, Gallery.ladder_dot(_case["loop_count"])
                    )
                }
            ),
        ]
    )
    return (input_diagram,)


@app.cell(hide_code=True)
def _(Gallery, mo):
    _cores = Gallery.max_cores()
    _deadline = Gallery.max_seconds()
    run_form = mo.md(
        "{arithmetic}\n\n{n_cores}\n\n{max_seconds}\n\n{ordered_cache}"
    ).batch(
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
    # Edit mode only: marimo's Interrupt reaches this kernel there (the served
    # app mode has no interrupt, so it keeps to cancellable background jobs).
    kernel_button = (
        mo.ui.run_button(label="Reduce in this kernel (Interrupt stops it)")
        if Gallery.mode(mo) == "edit"
        else None
    )
    mo.vstack(
        [
            Gallery.section(
                mo, "03 / run", "Reduce the selected case to retained masters"
            ),
            run_form,
            mo.hstack(
                [run_button] + ([kernel_button] if kernel_button is not None else []),
                justify="start",
            ),
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
    return kernel_button, run_button, run_form


@app.cell(hide_code=True)
def _(Gallery, mo):
    get_finished_job, set_finished_job = mo.state(None)
    # The run this session shows: its last one, or ?job=<id> from the URL.
    get_shown_job, set_shown_job = mo.state(
        Gallery.ladder_job(mo.query_params().get("job"))
    )
    return get_finished_job, get_shown_job, set_finished_job, set_shown_job


@app.cell(hide_code=True)
def _(Gallery, case_choice, mo, run_button, run_form, set_shown_job):
    # Script mode runs exactly the case named on the command line; the report
    # cell below prints its outcome and sets the exit status.
    headless = mo.app_meta().mode == "script"
    run_request = started_job = start_error = None
    try:
        if headless:
            run_request = Gallery.headless_request(mo.cli_args(), run_form.value)
        elif run_button.value:
            run_request = {
                **run_form.value,
                "case": case_choice.value,
                "cancel_after": None,
            }
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
    if started_job is not None and not headless:
        set_shown_job(started_job)
    # No output (no empty gap) unless the start was refused.
    mo.callout(start_error, kind="warn") if start_error else None
    return headless, run_request, started_job


@app.cell(hide_code=True)
def _(Gallery, case_choice, kernel_button, mo, run_form, set_shown_job):
    # In-kernel run (edit mode): the monitor card redraws this cell's output
    # about once a second; Interrupt stops the native call within a moment.
    mo.stop(kernel_button is None or not kernel_button.value)
    _request = {**run_form.value, "case": case_choice.value}
    if not Gallery.interruptible():
        mo.output.replace(
            Gallery.note(
                mo,
                "This host has no live progress API",
                "The reduction runs, but neither the monitor nor Interrupt reaches "
                "the native call before it returns. Use the background job instead.",
                pending=True,
            )
        )
    try:
        kernel_job = Gallery.kernel_work(
            _request["case"],
            n_cores=int(_request["n_cores"]),
            arithmetic=_request["arithmetic"],
        )
    except (OSError, ValueError) as _error:
        mo.stop(True, mo.callout(str(_error), kind="warn"))
    _progress = Gallery.kernel_progress(mo, title="In this kernel")
    _interrupted = False
    try:
        Gallery.reduce_in_kernel(kernel_job, progress=_progress)
    except KeyboardInterrupt:
        # Recorded as cancelled; the kernel and its engine stay usable.
        _interrupted = True
    set_shown_job(kernel_job)
    # After an Interrupt marimo re-runs no other cell (this state change
    # included): the run's final card and outcome stay in this cell.
    Gallery.kernel_outcome(mo, kernel_job, progress=_progress, interrupted=_interrupted)
    return (kernel_job,)


@app.cell(hide_code=True)
def _(Gallery, get_shown_job, headless, mo, set_shown_job, started_job):
    if headless:
        # A headless run watches, reports and cancels only its own job.
        watched_job = started_job
    else:
        # Keep showing this session's last run when a control changes, else
        # re-attach to a run still going (a reload, another session).
        watched_job = (
            get_shown_job() or started_job or Gallery.active_job(kind="ladder")
        )
        if watched_job is not None and watched_job != get_shown_job():
            set_shown_job(watched_job)
        _name = (watched_job or "").rsplit("/", 1)[-1]
        if Gallery.ladder_job(_name):
            # Background jobs survive a reload; in-kernel runs are per session.
            mo.query_params()["job"] = _name
    return (watched_job,)


@app.cell(hide_code=True)
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
    _kernel = job_status.get("kind") == "kernel"
    _parts = [
        Gallery.section(
            mo,
            "04 / progress",
            "The in-kernel run" if _kernel else "The background run",
        ),
        Gallery.job_monitor_card(mo, job_status, watched_job),
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
def _(Gallery, Path, job_status, ladder_cases, mo, watched_job):
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
        (f"{_seconds['total']:.1f} s", "whole run (import, data, reduction)"),
        *Gallery.execution_metrics(ladder_result["execution"]),
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
    _cache = ladder_result.get("ordered_cache")
    _data = ladder_result.get("data") or {}
    _data_text = (
        f"data {_data['version']} · {Path(str(_data.get('path', ''))).name}"
        if _data.get("version")
        else str(_data.get("source", "built-in 1L–3L registry"))
    )
    mo.vstack(
        [
            Gallery.section(
                mo,
                "05 / result",
                _case["label"],
                f"sha256 {_retained['sha256']} · {_data_text}",
            ),
            Gallery.metrics(mo, _metrics),
            mo.accordion(
                {
                    "Stage timings (wall and CPU seconds)": Gallery.table(
                        mo, _timings, numeric=("wall seconds", "CPU seconds")
                    ),
                    "Recorded numerical reference (not recomputed here)": mo.vstack(
                        [
                            mo.md(
                                "High-precision Laurent coefficients of this input at "
                                "mUV = μ = G = 1 (MSbar), from the rung fixture."
                            ),
                            Gallery.table(mo, _reference, numeric=("ε power",)),
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
def _(Gallery, ladder_result, mo, watched_job):
    # The retained-master result in the live pager, as this cell's own output
    # (an accordion would unmount and dispose it).
    _text = (
        Gallery.checked_job(watched_job) / ladder_result["retained"]["file"]
    ).read_text()
    retained_result = None
    try:
        from symbolica import E as _E

        retained_result = _E(_text)
        _formula = Gallery.formula(mo, retained_result)
    except (ImportError, ValueError, TypeError):
        _formula = mo.ui.code_editor(_text[:20000], language="python", disabled=True)
    mo.vstack(
        [
            mo.md("**Retained-master result**"),
            _formula,
            Gallery.legend(
                mo,
                "RustRedMaster(rustred_family_…, powers) are the catalogue masters; "
                "GC_11 is the UFO strong coupling, mursq = μ_R², mUV the auxiliary UV "
                "mass and ε the dimensional regulator (d = 4 − 2ε).",
            ),
        ]
    )
    return (retained_result,)


@app.cell(hide_code=True)
def _(Gallery, mo, retained_result):
    mo.stop(retained_result is None)
    retained_term = Gallery.term_selector(mo, retained_result)
    mo.stop(retained_term is None)
    retained_term
    return (retained_term,)


@app.cell(hide_code=True)
def _(Gallery, mo, retained_result, retained_term):
    # One term at a time, each part in its own view (this cell's output).
    Gallery.term_view(mo, retained_result, retained_term.value)
    return


@app.cell(hide_code=True)
def _(Gallery, case_choice, ladder_cases, mo):
    # The panel's controls follow the selected case (a new case starts with
    # fresh values and no result).
    _case = ladder_cases[case_choice.value]
    _loops = _case["loop_count"]
    laurent_parameters = Gallery.laurent_parameters(case_choice.value)
    _fields = {
        f"p{_index}": mo.ui.text(
            value="1",
            label=f"${_parameter['latex']}$ · {_parameter['meaning']}",
            full_width=True,
        )
        for _index, _parameter in enumerate(laurent_parameters)
    }
    _in_page = Gallery.laurent_in_page(mo, case_choice.value)
    laurent_controls = mo.md(
        "**Parameters**, exact decimals or fractions (e.g. `91.1876`, `1/3`)\n\n"
        + "\n\n".join(f"{{{key}}}" for key in _fields)
        + "\n\n**Accuracy and ε orders**\n\n{digits}\n\n{orders}"
    ).batch(
        digits=mo.ui.number(
            start=1,
            stop=Gallery.laurent_max_digits,
            value=32,
            step=1,
            label="Requested decimal digits",
        ),
        orders=mo.ui.range_slider(
            start=-_loops,
            stop=Gallery.laurent_highest_order,
            step=1,
            value=[-_loops, 0],
            show_value=True,
            label="ε orders (lowest, highest)",
        ),
        **_fields,
    )
    laurent_button = mo.ui.run_button(
        label="Evaluate the Laurent series",
        kind="success",
        disabled=not _in_page,
    )
    _orders_note = (
        f"The series starts at {Gallery.epsilon_text(-_loops)}. Vakint returns "
        "every order its master data supports; an order that needs master "
        "coefficients beyond the data fails closed and is shown with Vakint's "
        "MissingMasterOrders reason instead of a number. Exact parameter values "
        "never limit the precision; this input has no external-momentum or other "
        "numerator parameter left (q = 0, or its normalized q² coefficient), and "
        "the constants of the result (π, Clausen values, master coefficients) "
        "are evaluated by Vakint itself."
    )
    mo.vstack(
        [
            Gallery.section(
                mo,
                "06 / numbers",
                "Input graph → Laurent series",
                "What everything above boils down to: once every parameter has a "
                "value, the input graph is a Laurent series in ε with numerical "
                "coefficients. Vakint reduces the whole input with symbolic "
                "masters (kept for this session per input and order count), "
                "then evaluates it numerically at the requested precision in "
                "this kernel.",
            ),
            laurent_controls,
            laurent_button,
            Gallery.note(mo, "Which orders", _orders_note),
            *(
                []
                if _in_page
                else [
                    Gallery.note(
                        mo,
                        "Not evaluated in this application",
                        "Without marimo's Interrupt (application mode) this page "
                        "evaluates only the quick cases in its kernel. Open the "
                        "notebook in marimo's edit mode for this case.",
                        pending=True,
                    )
                ]
            ),
        ]
    )
    return laurent_button, laurent_controls, laurent_parameters


@app.cell(hide_code=True)
def _(case_choice, mo):
    # The last series of this case; a new case starts without one.
    case_choice.value
    get_laurent, set_laurent = mo.state(None)
    return get_laurent, set_laurent


@app.cell(hide_code=True)
def _(
    Gallery,
    case_choice,
    laurent_button,
    laurent_controls,
    laurent_parameters,
    mo,
    set_laurent,
):
    # Runs on the button only. The monitor card redraws this cell's output
    # while Vakint loads data and reduces; marimo's Interrupt stops it and this
    # cell then shows the interrupted card itself (marimo re-runs no other
    # cell after an Interrupt). A finished series is shown by the next cell.
    mo.stop(not laurent_button.value)
    _request = laurent_controls.value
    _progress = Gallery.kernel_progress(mo, title="Input graph → Laurent series")
    laurent = None
    _outcome = None
    try:
        laurent = Gallery.laurent_series(
            case_choice.value,
            {
                _parameter["name"]: _request[f"p{_index}"]
                for _index, _parameter in enumerate(laurent_parameters)
            },
            digits=_request["digits"],
            orders=_request["orders"],
            progress=_progress,
        )
    except KeyboardInterrupt:
        _outcome = mo.vstack(
            [
                Gallery.monitor_card(
                    mo,
                    _progress.model,
                    state="interrupted",
                    title="Input graph → Laurent series",
                ),
                Gallery.note(
                    mo,
                    "Interrupted in this kernel",
                    "Nothing of this evaluation is kept; the kernel and its "
                    "engines stay usable. Press the button again to restart it.",
                ),
            ]
        )
    except (OSError, ValueError, RuntimeError) as _error:
        _outcome = Gallery.note(
            mo, "The evaluation did not complete", str(_error)[-3000:], pending=True
        )
    if laurent is not None:
        # The next cell shows the series; this cell keeps no progress card.
        set_laurent(laurent)
        mo.output.clear()
    _outcome
    return


@app.cell(hide_code=True)
def _(Gallery, get_laurent, input_diagram, mo):
    # The last series stays while values change, until the next evaluation;
    # its caption names the values it was computed for.
    _laurent = get_laurent()
    mo.stop(
        _laurent is None,
        Gallery.legend(
            mo,
            "Set the values and press “Evaluate the Laurent series”; nothing "
            "runs before.",
        ),
    )
    Gallery.laurent_view(mo, _laurent, input_diagram)
    return


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
                    "Recent jobs on this server": Gallery.table(
                        mo, _rows, numeric=("elapsed (s)",)
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
