"""Presentation and authenticated evidence shared by the four Vakint apps.

Native physics objects stay in their owning HEPKit/Symbolica classes. This module
only handles presentation, evidence files, and the explicit benchmark launcher.
"""

from __future__ import annotations

import hashlib
import html
import json
import os
import subprocess
import sys
import uuid
from pathlib import Path
from urllib.parse import quote


class Gallery:
    directory = Path(__file__).resolve().parent
    display_receipt_bytes = 8 * 1024**2
    pages = (
        ("graphs", "01", "Graphs & numerators"),
        ("inspect", "02", "Constants & masters"),
        ("five-loop", "03", "Five-loop UV reduction"),
        ("closure", "04", "Pinches & caching"),
    )

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

    @staticmethod
    def engine(precision=32, mode="retain"):
        from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

        return Vakint(
            evaluation_order=[
                VakintEvaluationMethod.new_rustred_method(master_mode=mode)
            ],
            run_time_decimal_precision=precision,
            integral_normalization_factor="FMFTandMATAD",
            form_exe_path="/gallery-uses-rustred-only",
        )

    @staticmethod
    def load_data(engine, max_loop_count):
        # A configured checkout is recognized and validated before any network
        # request. Otherwise official loading uses the caller's current directory.
        return engine.load_official_external_data(
            max_loop_count,
            version="latest",
            local_artifacts_path=os.environ.get("VAKINT_GALLERY_ARTIFACTS"),
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
    def evidence_card(cls, mo, name, title):
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
                ("taylor_seconds", "Quadratic expansion"),
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

    @staticmethod
    def benchmark_available():
        return bool(
            os.environ.get("VAKINT_GALLERY_BENCHMARK")
            and os.environ.get("VAKINT_GALLERY_WORK")
            and os.environ.get("VAKINT_GALLERY_ARTIFACTS")
            and os.environ.get("VAKINT_GALLERY_PREPARED_INPUT_ROOT")
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

    @staticmethod
    def start_benchmark(projection, *, cached=False):
        """Start one serialized background job, leaving the app kernel responsive."""
        import fcntl

        if projection not in ("quadratic", "zero-momentum"):
            raise ValueError("Unknown physical projection")
        if not Gallery.benchmark_available():
            raise ValueError(
                "This gallery server needs a benchmark runner, workspace, external data and prepared scalar inputs"
            )
        prepared_input = Gallery.prepared_input_directory(projection)
        if cached and not os.environ.get("VAKINT_GALLERY_CACHE"):
            raise ValueError("This gallery server has no configured prepared cache")
        root = Path(os.environ["VAKINT_GALLERY_WORK"]).resolve()
        root.mkdir(parents=True, exist_ok=True)
        lock = (root / "benchmark.lock").open("a+")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close()
            raise RuntimeError(
                "A gallery benchmark is already running; wait for it to complete"
            ) from None
        work = root / f"{projection}-{uuid.uuid4().hex[:12]}"
        work.mkdir(parents=True)
        command = [
            os.environ["VAKINT_GALLERY_BENCHMARK"],
            "--projection",
            projection,
            "--artifacts",
            os.environ["VAKINT_GALLERY_ARTIFACTS"],
            "--work-dir",
            str(work),
            "--prepared-input",
            str(prepared_input),
        ]
        if cached:
            command.extend(["--cache", os.environ["VAKINT_GALLERY_CACHE"]])
        (work / "request.json").write_text(json.dumps({"command": command}))
        (work / "status.json").write_text(json.dumps({"state": "starting"}))
        try:
            with (work / "worker.log").open("w") as log:
                subprocess.Popen(
                    [
                        sys.executable,
                        str(Gallery.directory / "benchmark_job.py"),
                        str(work),
                        str(lock.fileno()),
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    pass_fds=(lock.fileno(),),
                )
        finally:
            # The child owns the inherited descriptor until its benchmark ends.
            lock.close()
        return str(work)

    @staticmethod
    def benchmark_status(path):
        work = Path(path).resolve()
        root = Path(os.environ["VAKINT_GALLERY_WORK"]).resolve()
        if not work.is_relative_to(root):
            raise ValueError("Benchmark path escapes its workspace")
        return json.loads((work / "status.json").read_text())
