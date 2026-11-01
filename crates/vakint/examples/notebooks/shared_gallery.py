"""Presentation, evidence and background jobs shared by the Vakint apps.

Native physics objects stay in their owning HEPKit/Symbolica classes. This module
only handles presentation, evidence files, the ladder cases and the explicit,
cancellable background jobs (``gallery_jobs``). Opening a page never starts a
reduction; every heavy run is an explicit, bounded and cancellable job.
"""

from __future__ import annotations

import hashlib
import html
import json
import os
import sys
import tempfile
import time
from itertools import pairwise
from pathlib import Path
from urllib.parse import quote

import gallery_jobs


class Gallery:
    directory = Path(__file__).resolve().parent
    ladder_directory = directory / "data" / "ladder"
    display_receipt_bytes = 8 * 1024**2
    pages = (
        ("graphs", "01", "Graphs & numerators"),
        ("inspect", "02", "Constants & masters"),
        ("five-loop", "03", "Five-loop: explore"),
        ("closure", "04", "Pinches & closure"),
        ("reduce", "05", "Reduce 3L & 4L"),
    )
    arithmetics = ("modular_polynomial", "modular_sampled", "exact")
    default_max_cores = 4
    default_max_seconds = 1800

    @classmethod
    def header(cls, mo, page, title, subtitle):
        prefix = os.environ.get("VAKINT_GALLERY_URL_PREFIX", "").rstrip("/")
        navigation = "".join(
            f'<a class="{"active" if key == page else ""}" '
            f'href="{html.escape(prefix)}/{key}/">'
            f"<span>{number}</span>{html.escape(label)}</a>"
            for key, number, label in cls.pages
        )
        return mo.Html(
            cls.styles()
            + '<header class="vg-hero"><div class="vg-eyebrow">'
            + "HEPKIT × VAKINT <span>INTERACTIVE FIELD NOTES</span></div>"
            + f"<h1>{html.escape(title)}</h1><p>{html.escape(subtitle)}</p>"
            + f'</header><nav class="vg-nav">{navigation}</nav>'
        )

    @staticmethod
    def math_font_css():
        prefix = os.environ.get("VAKINT_GALLERY_URL_PREFIX", "").rstrip("/")
        return (
            '<style>@font-face {font-family:"STIX Two Math";'
            f'src:url("{html.escape(prefix)}/assets/vakint-math.woff2") format("woff2");'
            "font-display:swap;}</style>"
        )

    @staticmethod
    def native_formula(mo, expression):
        """Keep native HTML and omission notes without a view-owned widget model."""
        payload = expression._repr_html_()
        if not isinstance(payload, str) or not payload.strip():
            raise ValueError("This expression has no native HTML presentation")
        return mo.Html('<div class="vg-native-formula">' + payload + "</div>")

    @staticmethod
    def styles():
        return (
            Gallery.math_font_css()
            + """<style>
        :root { --vg-navy:#16324c; --vg-teal:#087e86; --vg-muted:#61758b;
          --vg-line:#dbe5eb; --vg-paper:#f4f8fa; --vg-amber:#9f5707; }
        .vg-hero { background:linear-gradient(125deg,#102e46,#17566a);
          color:#fff; padding:2.5rem 2.7rem; border-radius:22px; margin:1rem 0; }
        .vg-eyebrow { font:600 11px/1.8 ui-sans-serif,system-ui;
          letter-spacing:.15em; color:#a8e6e1; }
        .vg-eyebrow span { margin-left:1.5rem; color:#cedce5; }
        .vg-hero h1 { font:600 clamp(29px,4vw,43px)/1.15 ui-sans-serif,system-ui;
          letter-spacing:-.035em; margin:1rem 0 .7rem; color:#fff; }
        .vg-hero p { font-size:16px; line-height:1.65; max-width:780px;
          color:#dce9ef; margin:0; }
        .vg-nav { display:flex; gap:7px; flex-wrap:wrap; margin:1.2rem 0 2.4rem; }
        .vg-nav a { color:var(--vg-muted); font-size:13px; text-decoration:none;
          padding:10px 13px; border:1px solid var(--vg-line); border-radius:9px; }
        .vg-nav span { font-size:10px; opacity:.65; margin-right:7px; }
        .vg-nav a:hover,.vg-nav a.active { background:#e7f5f3; color:#086c71;
          border-color:#acd9d3; }
        .vg-metrics { display:grid; grid-template-columns:repeat(auto-fit,minmax(155px,1fr));
          gap:12px; margin:1.2rem 0; }
        .vg-metric { background:var(--vg-paper); border:1px solid var(--vg-line);
          border-radius:15px; padding:20px 22px; }
        .vg-metric strong { display:block; color:var(--vg-navy); font-size:28px;
          font-weight:600; letter-spacing:-.03em; line-height:1.3; }
        .vg-metric span { display:block; font-size:12px; color:var(--vg-muted);
          margin-top:7px; line-height:1.5; }
        .vg-section { border-top:1px solid var(--vg-line); padding-top:1.7rem;
          margin-top:2rem; color:var(--vg-navy); }
        .vg-section small { color:var(--vg-teal); font-size:10px;
          letter-spacing:.15em; text-transform:uppercase; }
        .vg-section h2 { font-size:24px; letter-spacing:-.025em; margin:.35rem 0; }
        .vg-native-formula { max-width:100%; overflow-x:auto; }
        .vg-section p { color:var(--vg-muted); font-size:14px; line-height:1.6; }
        .vg-note { border:1px solid #cfe5e2; border-left:4px solid #45a69a;
          background:#f0f9f7; border-radius:10px; padding:17px 20px;
          margin:1rem 0; color:#234f54; font-size:14px; line-height:1.65;
          white-space:pre-line; overflow-wrap:anywhere; }
        .vg-note.pending { border-color:#efd7b8; border-left-color:#d99740;
          background:#fff8ec; color:#795124; }
        .vg-note strong { display:block; margin-bottom:4px; }
        .vg-foot { color:var(--vg-muted); border-top:1px solid var(--vg-line);
          font-size:11px; padding-top:1rem; margin-top:2rem; line-height:1.7; }
        .vakint-card { max-width:100%; overflow-x:auto; }
        .vakint-card pre,.vakint-card .value { max-height:18rem; overflow:auto; }
        @media(max-width:600px) { .vg-hero {padding:1.6rem;} .vg-nav a {flex:1 1 40%;}
          .vg-eyebrow span {display:block;margin:0;} }
        @media(prefers-color-scheme:dark) { :root { --vg-navy:#cbddea;
          --vg-muted:#a9bdcc; --vg-paper:#172936; --vg-line:#345064; }
          .vg-note {background:#173834;color:#c5e7df;} }
        </style>"""
        )

    @staticmethod
    def section(mo, number, title, description=""):
        return mo.Html(
            '<section class="vg-section">'
            f"<small>{html.escape(number)}</small><h2>{html.escape(title)}</h2>"
            f"<p>{html.escape(description)}</p></section>"
        )

    @staticmethod
    def metrics(mo, values):
        return mo.Html(
            '<div class="vg-metrics">'
            + "".join(
                '<div class="vg-metric">'
                f"<strong>{html.escape(str(value))}</strong>"
                f"<span>{html.escape(label)}</span></div>"
                for value, label in values
            )
            + "</div>"
        )

    @staticmethod
    def note(mo, title, message, *, pending=False):
        return mo.Html(
            f'<div class="vg-note{" pending" if pending else ""}">'
            f"<strong>{html.escape(title)}</strong>{html.escape(message)}</div>"
        )

    @classmethod
    def footer(cls, mo):
        runtime = cls.receipt("runtime")
        provenance = ""
        if isinstance(runtime, dict):
            provenance = "<br>" + " · ".join(
                f"{html.escape(str(label))}: {html.escape(str(value))}"
                for label, value in runtime.get("display", [])
            )
        return mo.Html(
            '<footer class="vg-foot">Native Symbolica expressions · shared HEPKit '
            "diagrams · precision-aware Vakint records<br>"
            '<a href="https://www.hepkit.org/gallery/notebooks/four_loop_numerator.py">'
            "HEPKit numerator tutorial</a> · "
            '<a href="https://github.com/alphal00p/vakint-external-data">'
            "External data and provenance</a>" + provenance + "</footer>"
        )

    @classmethod
    def max_cores(cls):
        """Ceiling on any gallery engine's or job's ``n_cores``.

        VAKINT_GALLERY_MAX_CORES (default 4), never more than the CPUs this
        server may use. A shared server stays responsive while a run works."""
        try:
            allowed = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            allowed = os.cpu_count() or 1
        try:
            limit = int(
                os.environ.get("VAKINT_GALLERY_MAX_CORES", cls.default_max_cores)
            )
        except ValueError:
            limit = cls.default_max_cores
        return max(1, min(limit, allowed))

    @classmethod
    def bounded_cores(cls, n_cores=None):
        """Requested cores (else VAKINT_GALLERY_N_CORES, else the ceiling),
        clamped to ``1..max_cores()``."""
        if n_cores is None and os.environ.get("VAKINT_GALLERY_N_CORES"):
            n_cores = int(os.environ["VAKINT_GALLERY_N_CORES"])
        ceiling = cls.max_cores()
        return ceiling if n_cores is None else max(1, min(int(n_cores), ceiling))

    @classmethod
    def engine(
        cls,
        precision=32,
        mode="retain",
        n_cores=None,
        ordered_weight_arithmetic=None,
        normalization="FMFTandMATAD",
        **options,
    ):
        """RustRed-only Ordered engine on at most ``max_cores()`` workers.

        Cores and the Ordered weight arithmetic (exact | modular_sampled |
        modular_polynomial) come from the arguments, else VAKINT_GALLERY_N_CORES
        and VAKINT_GALLERY_ARITHMETIC, else the core ceiling and Vakint's default
        arithmetic (modular_polynomial). Other Vakint options pass through."""
        from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

        if ordered_weight_arithmetic is None:
            ordered_weight_arithmetic = os.environ.get("VAKINT_GALLERY_ARITHMETIC")
        return Vakint(
            evaluation_order=[
                VakintEvaluationMethod.new_rustred_method(master_mode=mode)
            ],
            run_time_decimal_precision=precision,
            integral_normalization_factor=normalization,
            form_exe_path="/gallery-uses-rustred-only",
            n_cores=cls.bounded_cores(n_cores),
            ordered_weight_arithmetic=ordered_weight_arithmetic,
            **options,
        )

    @staticmethod
    def execution(engine):
        """The engine's Ordered weight arithmetic and worker cores, for display."""
        cores = engine.n_cores
        return f"{engine.ordered_weight_arithmetic} · {cores} core{'' if cores == 1 else 's'}"

    @staticmethod
    def data_directory():
        """External data: a configured checkout (VAKINT_GALLERY_ARTIFACTS), else
        a download cache (VAKINT_GALLERY_DATA_CACHE, else the user cache
        directory). Never the server's working directory."""
        configured = os.environ.get("VAKINT_GALLERY_ARTIFACTS")
        if configured:
            return Path(configured).resolve()
        cache = os.environ.get("VAKINT_GALLERY_DATA_CACHE")
        if not cache:
            base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
            cache = Path(base) / "vakint-gallery"
        directory = Path(cache).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        return directory

    @classmethod
    def load_data(cls, engine, max_loop_count):
        # A configured checkout is recognized and validated before any network
        # request. Otherwise the official packs install once into the cache.
        return engine.load_official_external_data(
            max_loop_count,
            version="latest",
            local_artifacts_path=str(cls.data_directory()),
        )

    @classmethod
    def evidence_file(cls, name, artifact_index=None):
        """Resolve only manifest-listed files beneath the evidence directory."""
        root = Path(os.environ.get("VAKINT_GALLERY_EVIDENCE", cls.directory / "data"))
        manifest_path = root / "manifest.json"
        if manifest_path.stat().st_size > 1024**2:
            raise ValueError("Gallery evidence manifest exceeds its size limit")
        manifest = json.loads(manifest_path.read_text())
        if manifest["schema"] != 1:
            raise ValueError("Unsupported gallery evidence schema")
        entry = manifest["receipts"].get(name)
        if entry is None:
            return None
        if artifact_index is not None:
            artifacts = entry.get("artifacts", [])
            if not 0 <= artifact_index < len(artifacts):
                return None
            entry = artifacts[artifact_index]
        path = (root / entry["path"]).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError("Evidence path escapes its dataset")
        return path, entry["sha256"]

    @classmethod
    def receipt(cls, name):
        """Read a compact checked receipt; retain large evidence as downloads."""
        resolved = cls.evidence_file(name)
        if resolved is None:
            return None
        path, digest = resolved
        if path.stat().st_size > cls.display_receipt_bytes:
            raise ValueError(
                "The full evidence exceeds the display limit; download it or select a compact measurement summary"
            )
        payload = path.read_bytes()
        if hashlib.sha256(payload).hexdigest() != digest:
            raise ValueError(f"Evidence checksum mismatch: {name}")
        return json.loads(payload)

    @staticmethod
    def display_receipt(value):
        """Bound browser output without changing the stored scientific record."""
        budget = [400]

        def preview(item, depth=0):
            budget[0] -= 1
            if budget[0] < 0 or depth > 6:
                return "… preview limit; complete data in the downloadable record"
            if isinstance(item, str):
                return (
                    item
                    if len(item) <= 2048
                    else item[:2048] + f"… [{len(item):,} characters total]"
                )
            if isinstance(item, list):
                values = [preview(v, depth + 1) for v in item[:8]]
                return (
                    values
                    if len(item) <= 8
                    else {"count": len(item), "preview": values}
                )
            if isinstance(item, dict):
                values = {k: preview(v, depth + 1) for k, v in list(item.items())[:32]}
                if len(item) > 32:
                    values["additional_fields"] = len(item) - 32
                return values
            return item

        return preview(value)

    @classmethod
    def evidence_links(cls, mo, name):
        prefix = os.environ.get("VAKINT_GALLERY_URL_PREFIX", "").rstrip("/")
        url = f"{prefix}/evidence/{quote(name, safe='')}"
        links = [
            f'<a href="{html.escape(url)}" download>Download the complete measurement record</a>'
        ]
        try:
            if cls.evidence_file(name) is None:
                return mo.md("")
            for index in range(16):
                resolved = cls.evidence_file(name, index)
                if resolved is None:
                    break
                path, digest = resolved
                links.append(
                    f'<a href="{html.escape(url)}/{index}" download title="SHA-256 {html.escape(digest)}">{html.escape(path.name)}</a>'
                )
        except (OSError, ValueError):
            return mo.md("")
        return mo.Html('<p class="vg-foot">' + " · ".join(links) + "</p>")

    @classmethod
    def evidence_card(cls, mo, name, title, *, optional=False):
        """A checked receipt's card; an absent ``optional`` receipt shows nothing."""
        try:
            receipt = cls.receipt(name)
        except (OSError, ValueError) as error:
            return mo.vstack(
                [
                    cls.note(
                        mo,
                        f"{title} · evidence not displayed",
                        str(error),
                        pending=True,
                    ),
                    cls.evidence_links(mo, name),
                ]
            )
        if receipt is None and optional:
            return mo.md("")
        if receipt is None:
            return cls.note(
                mo,
                f"{title} — awaiting validated measurements",
                "Measurements will appear here once their exact reduction results have been checked.",
                pending=True,
            )
        details = mo.vstack(
            [
                mo.accordion(
                    {
                        title + " · bounded measurement preview": mo.json(
                            cls.display_receipt(receipt)
                        )
                    }
                ),
                cls.evidence_links(mo, name),
            ]
        )
        schema = receipt.get("schema") if isinstance(receipt, dict) else None
        if schema == "vakint.gallery.coverage-advisory.v1":
            outcome = receipt["outcome"]
            return mo.vstack(
                [
                    cls.note(
                        mo,
                        "Recorded coverage advisory · reduction continued",
                        receipt["recorded_warning"],
                        pending=True,
                    ),
                    cls.metrics(
                        mo,
                        [
                            (
                                outcome["catalogue_master_count"],
                                "declared master identities",
                            ),
                            (outcome["unknown_count"], "surviving unknown integrals"),
                        ],
                    ),
                    cls.note(
                        mo, "What this example establishes", receipt["scope_note"]
                    ),
                    details,
                ]
            )
        if schema == "vakint.gallery.physical-status.v1":
            native = receipt.get("running_cache", {})
            fields = [
                (receipt.get("state", "Not recorded").replace("_", " "), title),
            ]
            for key, label in (
                ("accepted_records", "Cache records enqueued"),
                ("committed_records", "Cache records acknowledged durable"),
            ):
                if native.get(key) is not None:
                    fields.append((f"{native[key]:,}", label))
            return mo.vstack(
                [
                    cls.metrics(mo, fields),
                    cls.note(
                        mo,
                        "Recorded progress; final comparison pending",
                        f"Observed {receipt.get('observed_utc', 'at an unrecorded time')}. This immutable snapshot does not refresh automatically. Durable intermediate records are not a completed reduction or a measured speedup. Fresh-process reuse of the same input is a seen-input check, not a held-out benchmark.",
                        pending=True,
                    ),
                    details,
                ]
            )
        if schema == "vakint.qcd-ladder.benchmark.v1":
            prepared_input = receipt.get("input_mode") == "prepared"
            completed = receipt.get("reduction_completed") is True
            state = (
                "Retained symbolic result" if completed else "Reduction not completed"
            )
            if not completed and receipt.get("retained_result_available"):
                state = "Partial retained result"
            if "reduction_seconds" not in receipt:
                state = (
                    "Saved scalar input verified"
                    if prepared_input
                    else "Preparation recorded"
                )
            fields = [(state, title)]
            closure = receipt.get("catalogue_closure", {})
            if closure.get("unknown_count"):
                fields.append(
                    (
                        closure["unknown_count"],
                        "surviving non-catalogue integrals; full list downloadable",
                    )
                )
            for key, label in (
                ("reduction_seconds", "Native reduction · seconds"),
                ("preparation_seconds", "Fresh native preparation · seconds"),
                ("input_loading_seconds", "Verified input loading · seconds"),
                (
                    "end_to_end_seconds",
                    "Saved input to result · seconds"
                    if prepared_input
                    else "Graph to result · seconds",
                ),
            ):
                if receipt.get(key) is not None:
                    fields.append((f"{receipt[key]:,.3f}", label))
            timing_rows = []
            preparation = receipt.get("preparation") or {}
            for key, label in (
                ("graph_generation_seconds", "Graph and Feynman rules"),
                ("scalar_contraction_seconds", "Native trace contraction"),
                ("taylor_seconds", "Taylor stage (when applicable)"),
                ("angular_seconds", "Tensor reduction"),
                ("ingress_seconds", "Vacuum integral construction"),
            ):
                if key in preparation:
                    timing_rows.append({"Stage": label, "Seconds": preparation[key]})
            for key, label in (
                ("input_loading_seconds", "Verified saved scalar input"),
                ("data_loading_seconds", "External data"),
                ("cache_loading_seconds", "Persistent cache"),
                ("reduction_seconds", "Native reduction"),
            ):
                if receipt.get(key) is not None:
                    timing_rows.append({"Stage": label, "Seconds": receipt[key]})
            parts = [cls.metrics(mo, fields)]
            if prepared_input:
                parts.append(
                    cls.note(
                        mo,
                        "Verified saved scalar input",
                        "The complete prepared numerator and projector conditions were loaded and checked. Original graph, trace and Taylor preparation are historical provenance, excluded from this run's timings.",
                    )
                )
            if timing_rows:
                parts.append(mo.ui.table(timing_rows, selection=None, pagination=False))
            if receipt.get("reduction_error"):
                parts.append(
                    cls.note(
                        mo,
                        "Reduction stopped",
                        cls.display_receipt(receipt["reduction_error"]),
                        pending=True,
                    )
                )
            warnings = receipt.get("warnings", [])
            for warning in warnings[:4]:
                parts.append(
                    cls.note(
                        mo,
                        "Reduction warning · preview",
                        cls.display_receipt(warning),
                        pending=True,
                    )
                )
            if len(warnings) > 4:
                parts.append(
                    cls.note(
                        mo,
                        f"{len(warnings):,} warnings recorded",
                        "The complete warnings and unknown-integral identities are preserved in the downloadable measurement record.",
                        pending=True,
                    )
                )
            return mo.vstack([*parts, details])
        if schema == "vakint.five-loop-cache-build.v1":
            fields = [
                (
                    "Complete" if receipt.get("complete") else "Incomplete",
                    "Selected cache workload",
                ),
                (f"{receipt['elapsed_seconds'] / 3600:.2f} h", "Preparation time"),
                (str(receipt["selected_training_roots"]), "Selected roots"),
            ]
            if receipt.get("compressed_bytes") is not None:
                fields.append(
                    (
                        f"{receipt['compressed_bytes'] / 1048576:.2f} MiB",
                        "Compressed cache",
                    )
                )
            return mo.vstack([cls.metrics(mo, fields), details])
        return details

    # -- Background jobs ---------------------------------------------------
    @staticmethod
    def work_root():
        """Job directories: VAKINT_GALLERY_WORK, else a per-user temporary root."""
        configured = os.environ.get("VAKINT_GALLERY_WORK")
        if configured:
            return Path(configured).resolve()
        return Path(tempfile.gettempdir()).resolve() / f"vakint-gallery-{os.getuid()}"

    @classmethod
    def job_cpus(cls, n_cores):
        """The first ``n_cores`` CPUs this server may use: a job's affinity."""
        try:
            allowed = sorted(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            return None
        return allowed[: max(1, n_cores)]

    @classmethod
    def max_seconds(cls, requested=None):
        """A job's wall-clock deadline: requested, else VAKINT_GALLERY_MAX_SECONDS,
        else 30 minutes; at most the server's own ceiling."""
        try:
            ceiling = float(
                os.environ.get("VAKINT_GALLERY_MAX_SECONDS", cls.default_max_seconds)
            )
        except ValueError:
            ceiling = float(cls.default_max_seconds)
        if requested is None:
            return ceiling
        requested = float(requested)
        if not requested > 0:
            raise ValueError("A job deadline must be a positive number of seconds")
        return min(requested, ceiling)

    @classmethod
    def job_request(cls, kind, command, n_cores, max_seconds, **fields):
        """Command, bounded threads, CPUs, niceness and deadline of one job."""
        return {
            "kind": kind,
            "command": command,
            "n_cores": n_cores,
            "environment": {
                "RAYON_NUM_THREADS": str(n_cores),
                "VAKINT_N_CORES": str(n_cores),
            },
            "cpus": cls.job_cpus(n_cores),
            "nice": 10,
            "max_seconds": max_seconds,
            **fields,
        }

    @classmethod
    def jobs(cls):
        return gallery_jobs.jobs(cls.work_root())

    @classmethod
    def active_job(cls, kind=None):
        """The running job (of one kind: "ladder" or "five-loop"), also after a
        reload or from another session."""
        return gallery_jobs.active_job(cls.work_root(), kind=kind)

    @classmethod
    def checked_job(cls, path):
        work = Path(path).resolve()
        if not work.is_relative_to(cls.work_root()):
            raise ValueError("Job path escapes its workspace")
        return work

    @classmethod
    def job_status(cls, path):
        """The supervisor's status, with live stage, native progress and memory.

        The supervisor writes its status every five seconds; while the job runs,
        the page overlays the files the command updates continuously."""
        import benchmark_job

        work = cls.checked_job(path)
        current = gallery_jobs.status(work)
        if current.get("state") in gallery_jobs.FINAL_STATES:
            return current
        live = dict(current)
        stage = benchmark_job.read_stage(work)
        if stage:
            live["stage"] = stage
        live["native_progress"] = benchmark_job.read_native_progress(work)
        if isinstance(current.get("worker_pid"), int):
            live["resources"] = benchmark_job.worker_resources(current["worker_pid"])
        return live

    @classmethod
    def cancel_job(cls, path):
        """Stop every process of a job (SIGTERM, then SIGKILL) and confirm it."""
        return gallery_jobs.cancel(cls.checked_job(path))

    @classmethod
    def wait_job(cls, path, *, timeout=None, cancel_after=None, poll_seconds=0.5):
        """Block until a job ends (script mode). ``cancel_after`` seconds cancels
        it exactly as the Cancel button does; the report joins the status."""
        started = time.monotonic()
        report = None
        while True:
            current = cls.job_status(path)
            if current.get("state") in gallery_jobs.FINAL_STATES:
                return (
                    current if report is None else {**current, "cancel_report": report}
                )
            elapsed = time.monotonic() - started
            if cancel_after is not None and report is None and elapsed >= cancel_after:
                report = cls.cancel_job(path)
                continue
            if timeout is not None and elapsed >= timeout:
                return current
            time.sleep(poll_seconds)

    @classmethod
    def cancel_report(cls, path):
        """What the last Cancel of this job signalled, or None."""
        report = cls.checked_job(path) / "cancel-report.json"
        try:
            return json.loads(report.read_text())
        except (OSError, ValueError):
            return None

    @staticmethod
    def compact_json(value):
        return json.dumps(value, ensure_ascii=False) if value else "none"

    @classmethod
    def job_result(cls, path):
        """The job's own receipt (``benchmark.json``), or None."""
        receipt = cls.checked_job(path) / "benchmark.json"
        if not receipt.is_file() or receipt.stat().st_size > cls.display_receipt_bytes:
            return None
        return json.loads(receipt.read_text())

    @staticmethod
    def progress_rows(status):
        """A compact progress summary of a job status for the browser."""
        rows = [("State", str(status.get("state", "unknown")).replace("_", " "))]
        if status.get("elapsed_seconds") is not None:
            rows.append(("Elapsed", f"{status['elapsed_seconds']:,.1f} s"))
        stage = status.get("stage") or {}
        if stage.get("stage"):
            rows.append(("Stage", str(stage["stage"])))
        native = (status.get("native_progress") or {}).get("snapshot") or {}
        monitor = native.get("monitor") if isinstance(native, dict) else None
        monitor = monitor if isinstance(monitor, dict) else {}
        phase = monitor.get("phase") or native.get("caller_phase")
        if phase:
            rows.append(("Native phase", str(phase).replace("_", " ")))
        ordered = native.get("ordered_reduction") or monitor.get("ordered_reduction")
        if isinstance(ordered, dict):
            if ordered.get("total_phases"):
                rows.append(
                    (
                        "Ordered phases",
                        f"{ordered.get('completed_phases', 0)} / {ordered['total_phases']}"
                        + (
                            f" ({ordered['phase_kind']})"
                            if ordered.get("phase_kind")
                            else ""
                        ),
                    )
                )
            for key, label in (
                ("frontier_keys", "Frontier keys"),
                ("processed_keys", "Processed keys"),
                ("terminal_keys", "Terminal keys"),
            ):
                if isinstance(ordered.get(key), int):
                    rows.append((label, f"{ordered[key]:,}"))
        resources = status.get("resources") or {}
        if resources.get("worker_rss_bytes"):
            rows.append(
                ("Worker RSS", f"{resources['worker_rss_bytes'] / 1024**3:.2f} GiB")
            )
        return rows

    # -- Ladder reductions (3L, 4L) -----------------------------------------
    @staticmethod
    def ladder_dot(loop_count):
        """DOT of the massless-quark QCD ladder self-energy with L-1 gluon rungs.

        The rung fixture's graph (``ladder_target`` in tests/fixtures/qcd_ladder);
        L=5 is the four-rung ``ladder.dot`` with the same edges."""
        if loop_count < 1:
            raise ValueError("A quark-loop ladder needs at least one loop")
        rungs = range(1, loop_count)
        rail = ["left", *(f"t{rung}" for rung in rungs), "right"]
        rail += [*(f"b{rung}" for rung in reversed(rungs)), "left"]
        lines = [
            "digraph massless_quark_ladder_self_energy {",
            "    incoming [style=invis];",
            "    outgoing [style=invis];",
            '    incoming -> left [particle="g"];',
            '    right -> outgoing [particle="g"];',
            *(f'    {a} -> {b} [particle="d"];' for a, b in pairwise(rail)),
            *(f'    t{rung} -> b{rung} [particle="g"];' for rung in rungs),
            "}",
        ]
        return "\n".join(lines) + "\n"

    @classmethod
    def ladder_cases(cls):
        """The prepared ladder inputs and their recorded results (data/ladder)."""
        manifest = json.loads((cls.ladder_directory / "cases.json").read_text())
        if manifest.get("schema") != "vakint.gallery.ladder-cases.v1":
            raise ValueError("Unsupported ladder case manifest")
        return manifest["cases"]

    @classmethod
    def ladder_input(cls, key):
        """Path of a case's prepared scalar input; the worker verifies its sha256.

        Small inputs ship in data/ladder. Larger ones (4L quadratic, 1.7 MB) are
        read from VAKINT_GALLERY_LADDER_INPUTS when the server provides them."""
        case = cls.ladder_cases().get(key)
        if case is None:
            raise ValueError(f"Unknown ladder case: {key}")
        name = case["input"]["file"]
        if case["input"]["shipped"]:
            path = cls.ladder_directory / name
        else:
            configured = os.environ.get("VAKINT_GALLERY_LADDER_INPUTS")
            if not configured:
                raise ValueError(
                    f"{case['label']}: its prepared input is not shipped; set "
                    f"VAKINT_GALLERY_LADDER_INPUTS to a directory holding {name}"
                )
            path = Path(configured).resolve() / name
        if not path.is_file():
            raise ValueError(f"Prepared input is missing: {path}")
        return path

    @classmethod
    def ladder_case_available(cls, key):
        try:
            cls.ladder_input(key)
        except (OSError, ValueError):
            return False
        return True

    @classmethod
    def start_reduction(
        cls,
        case,
        *,
        n_cores=None,
        arithmetic=None,
        max_seconds=None,
        ordered_cache=True,
    ):
        """Reduce one ladder case in a background job; return its directory.

        The job runs reduce_worker.py on at most ``max_cores()`` CPUs with a
        wall-clock deadline, and stops on ``cancel_job``."""
        cases = cls.ladder_cases()
        if case not in cases:
            raise ValueError(f"Unknown ladder case: {case}")
        cls.ladder_input(case)
        arithmetic = arithmetic or "modular_polynomial"
        if arithmetic not in cls.arithmetics:
            raise ValueError(f"Unknown Ordered weight arithmetic: {arithmetic}")
        cores = cls.bounded_cores(n_cores)
        deadline = cls.max_seconds(max_seconds)
        root = cls.work_root()
        spec = {
            "case": case,
            "n_cores": cores,
            "arithmetic": arithmetic,
            "max_seconds": deadline,
            "data_directory": str(cls.data_directory())
            if cases[case]["loop_count"] >= 4
            or os.environ.get("VAKINT_GALLERY_ARTIFACTS")
            else None,
            "ordered_cache_directory": str(root / "ordered-cache")
            if ordered_cache
            else None,
        }
        return gallery_jobs.launch(
            root,
            case,
            lambda work: cls.job_request(
                "ladder",
                [sys.executable, str(cls.directory / "reduce_worker.py"), str(work)],
                cores,
                deadline,
                spec=spec,
            ),
        )

    # -- Five-loop exploration (no reduction) -------------------------------
    @classmethod
    def five_loop_plan(cls, root=None):
        """Summarize the 5L Ordered reduction plan from its JSON index, or None.

        Reads only ``5L/ibps/selection.json`` (a few MB) of a configured checkout
        or of an installed official pack; never downloads or decodes rules."""
        if root is None:
            configured = os.environ.get("VAKINT_GALLERY_ARTIFACTS")
            cache = os.environ.get("VAKINT_GALLERY_DATA_CACHE")
            base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
            root = Path(configured or cache or Path(base) / "vakint-gallery")
        root = Path(root)
        candidates = [root / "5L/ibps/selection.json"]
        candidates += sorted(
            root.glob(".vakint-external-data/*/5L/ibps/selection.json"), reverse=True
        )
        path = next(
            (candidate for candidate in candidates if candidate.is_file()), None
        )
        if path is None or path.stat().st_size > 64 * 1024**2:
            return None
        selection = json.loads(path.read_text())
        manifest = path.parents[1] / "manifest.json"
        version = (
            json.loads(manifest.read_text()).get("version", "unknown")
            if manifest.is_file()
            else "unknown"
        )
        supports = {}
        for owner in selection["owners"]:
            support = str(owner["mask"]).count("1")
            row = supports.setdefault(support, {"owners": 0, "bytes": 0})
            row["owners"] += 1
            row["bytes"] += int(owner["bytes"])
        return {
            "version": version,
            "path": str(path),
            "owners": len(selection["owners"]),
            "route_records": int(selection.get("route_record_count", 0)),
            "bundle_bytes": int(selection.get("total_bundle_bytes", 0)),
            "supports": [
                {
                    "owner support (lines)": support,
                    "owners": row["owners"],
                    "rule bundles (MiB)": round(row["bytes"] / 1024**2, 1),
                }
                for support, row in sorted(supports.items(), reverse=True)
            ],
        }

    @classmethod
    def saved_five_loop_result(cls):
        """The saved, checksummed 5L zero-momentum result (no reduction)."""
        examples = str(cls.directory.parent)
        if examples not in sys.path:
            sys.path.append(examples)
        import symbolica
        from q0_coefficients import SavedQ0Result

        return SavedQ0Result.load(api=symbolica)

    # -- Five-loop benchmark (explore-only unless the server opts in) --------
    @staticmethod
    def benchmark_available():
        return bool(
            os.environ.get("VAKINT_GALLERY_BENCHMARK")
            and os.environ.get("VAKINT_GALLERY_WORK")
            and os.environ.get("VAKINT_GALLERY_ARTIFACTS")
            and os.environ.get("VAKINT_GALLERY_PREPARED_INPUT_ROOT")
        )

    @classmethod
    def five_loop_runs_enabled(cls):
        """5L reductions are opt-in per server: VAKINT_GALLERY_ALLOW_5L=1 plus
        a configured runner. Otherwise the five-loop page is explore-only."""
        return (
            os.environ.get("VAKINT_GALLERY_ALLOW_5L") == "1"
            and cls.benchmark_available()
        )

    @staticmethod
    def prepared_input_directory(projection):
        """Select server-owned scalar input; the Rust loader validates its contents."""
        if projection not in ("quadratic", "zero-momentum"):
            raise ValueError("Unknown physical projection")
        configured = os.environ.get("VAKINT_GALLERY_PREPARED_INPUT_ROOT")
        if not configured:
            raise ValueError(
                "This gallery server has no configured prepared scalar inputs"
            )
        directory = Path(configured).resolve() / projection
        if not directory.is_dir():
            raise ValueError(f"Prepared scalar-input directory is missing: {directory}")
        manifest = directory / "prepared-input.json"
        if not manifest.is_file():
            raise ValueError(f"Prepared scalar-input manifest is missing: {manifest}")
        return directory

    @classmethod
    def start_benchmark(cls, projection, *, n_cores=None, max_seconds=None):
        """Start the opt-in five-loop job; it is bounded and cancellable."""
        if projection != "zero-momentum":
            raise ValueError(
                "Only the zero-momentum fixture is enabled; quadratic reduction is deferred"
            )
        if not cls.benchmark_available():
            raise ValueError(
                "This gallery server needs a benchmark runner, workspace, external data and prepared scalar inputs"
            )
        if os.environ.get("VAKINT_GALLERY_ALLOW_5L") != "1":
            raise ValueError(
                "Five-loop reductions are disabled on this server (explore-only); "
                "the operator enables them with VAKINT_GALLERY_ALLOW_5L=1"
            )
        prepared_input = cls.prepared_input_directory(projection)
        cores = cls.bounded_cores(n_cores)
        deadline = cls.max_seconds(max_seconds)

        def request(work):
            command = [
                os.environ["VAKINT_GALLERY_BENCHMARK"],
                "--projection",
                projection,
                "--artifacts",
                os.environ["VAKINT_GALLERY_ARTIFACTS"],
                "--work-dir",
                str(work),
                "--progress-snapshot",
                str(work / "monitor.json"),
                "--prepared-input",
                str(prepared_input),
            ]
            return cls.job_request("five-loop", command, cores, deadline)

        return gallery_jobs.launch(cls.work_root(), projection, request)

    @classmethod
    def benchmark_status(cls, path):
        return cls.job_status(path)
