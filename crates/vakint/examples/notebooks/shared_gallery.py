"""Presentation, evidence and background jobs shared by the Vakint apps.

Native physics objects stay in their owning HEPKit/Symbolica classes. This module
only handles presentation, evidence files, the ladder cases and the explicit,
cancellable background jobs (``gallery_jobs``). Opening a page never starts a
reduction; every heavy run is an explicit, bounded and cancellable job.
"""

from __future__ import annotations

import decimal
import hashlib
import html
import json
import os
import re
import sys
import tempfile
import time
from fractions import Fraction
from itertools import pairwise
from pathlib import Path
from typing import ClassVar
from urllib.parse import quote

import gallery_jobs


class Gallery:
    directory = Path(__file__).resolve().parent
    ladder_directory = directory / "data" / "ladder"
    display_receipt_bytes = 8 * 1024**2
    pages = (
        ("graphs", "01", "Graphs & numerators", "01_graphs_and_numerators.py"),
        ("inspect", "02", "Constants & masters", "02_constants_and_masters.py"),
        ("five-loop", "03", "Five-loop: explore", "03_five_loop_uv.py"),
        ("closure", "04", "Pinches & closure", "04_pinches_and_cache.py"),
        ("reduce", "05", "Reduce 3L & 4L", "05_ladder_reductions.py"),
    )
    five_loop_fixture_directory = directory / "data" / "five-loop-q0"
    arithmetics = ("modular_polynomial", "modular_sampled", "exact")
    default_max_cores = 4
    default_max_seconds = 1800
    final_states = gallery_jobs.FINAL_STATES

    @staticmethod
    def mode(mo):
        """marimo's mode: "edit", "run" (the served gallery) or "script"."""
        try:
            return mo.app_meta().mode
        except (AttributeError, RuntimeError):
            return None

    @staticmethod
    def edit_link(source, **query):
        """marimo edit server URL of a notebook file, opened in marimo's app
        view (outputs only; its toolbar switches to the code view)."""
        extra = "".join(f"&{key}={quote(str(value))}" for key, value in query.items())
        return f"/?file={quote(source)}&view-as=present{extra}"

    @classmethod
    def header(cls, mo, page, title, subtitle):
        prefix = os.environ.get("VAKINT_GALLERY_URL_PREFIX", "").rstrip("/")
        # The gallery server mounts one route per app; marimo's edit server
        # opens each notebook file of its directory instead.
        edit = cls.mode(mo) == "edit" and not prefix

        def link(key, source):
            return cls.edit_link(source) if edit else f"{prefix}/{key}/"

        navigation = "".join(
            f'<a class="{"active" if key == page else ""}" '
            f'href="{html.escape(link(key, source))}" target="_self">'
            f"<span>{number}</span>{html.escape(label)}</a>"
            for key, number, label, source in cls.pages
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

    # -- Expressions: static previews and the live Symbolica/HEPKit pager ----
    pager_install_hint = 'pip install "symbolica[notebook-display]"'
    # Symbolica's own preview footers name the live-notebook route; the
    # gallery replaces only that sentence and keeps every omission note.
    _pager_sentences = re.compile(
        r"Use \.formatted\([^)]*\) in a live notebook with symbolica\[notebook\] "
        r"to navigate\.|Navigation requires a live notebook with "
        r"symbolica\[notebook-display\]\."
    )

    @staticmethod
    def is_tensor(value):
        return type(value).__name__ == "TensorExpression"

    @staticmethod
    def live_pager_available(tensor=False):
        """The live pager needs anywidget (and Typst for tensor pages)."""
        import importlib.util

        names = ("anywidget", "typst") if tensor else ("anywidget",)
        return all(importlib.util.find_spec(name) is not None for name in names)

    @classmethod
    def static_formula(cls, mo, expression, *, note=None):
        """Native static MathML preview with its omission notes; no widget model.

        For small values and for containers that unmount their content (tabs,
        accordions, tables). ``note`` replaces Symbolica's sentence about the
        live notebook route; the preview's item count and markers stay."""
        payload = expression._repr_html_()
        if not isinstance(payload, str) or not payload.strip():
            raise ValueError("This expression has no native HTML presentation")
        if "Page rendering failed" in payload and cls.is_tensor(expression):
            # Tensor pages need Typst; the Symbolica MathML of the same value does not.
            return cls.static_formula(mo, expression.to_expression(), note=note)
        if note is not None:
            payload = cls._pager_sentences.sub(html.escape(note), payload)
        return mo.Html('<div class="vg-native-formula">' + payload + "</div>")

    native_formula = static_formula

    @classmethod
    def formula(cls, mo, value, *, page_size=25):
        """The live Symbolica/HEPKit pager of a value: page through every item,
        open numbered omission markers, zoom with + / - / 0.

        Place the result directly in a cell's output (alone or in a vstack),
        never inside mo.ui.tabs, mo.accordion or a table: those unmount their
        panels, which disposes the pager's widget model. Without the display
        packages this falls back to the static preview and says what to install."""
        tensor = cls.is_tensor(value)
        if tensor and not cls.live_pager_available(tensor=True):
            value, tensor = value.to_expression(), False
        if not cls.live_pager_available():
            return cls.static_formula(
                mo,
                value,
                note=f"Live pager unavailable here: {cls.pager_install_hint}.",
            )
        if tensor:
            return value.paged(page_size=page_size)
        return value.formatted(page_size=page_size)

    @classmethod
    def product_factors(cls, value, *, minimum=4):
        """The factors of a value that is one long product, else None.

        A product never wraps in the pagers (they page and wrap sums): its
        factors listed one per line stay readable at any width."""
        try:
            from symbolica import AtomType

            tensor = cls.is_tensor(value)
            expression = value.to_expression() if tensor else value
            if expression.get_type() != AtomType.Mul or len(expression) < minimum:
                return None
            factors = list(expression)
        except (AttributeError, ImportError, TypeError, ValueError):
            return None
        if tensor:
            # Back in tensor notation (named indices, slots) for display.
            try:
                return [type(value)(factor) for factor in factors]
            except (TypeError, ValueError):
                pass
        return factors

    @classmethod
    def factor_list(cls, mo, factors):
        """Static formulas, one factor per line (safe inside an accordion)."""
        return mo.vstack([cls.static_formula(mo, factor) for factor in factors])

    # -- One term of a large retained result ----------------------------------
    _master_power = re.compile(r"\)\^(\d+)$")

    @classmethod
    def term_parts(cls, term):
        """(masters, other factors, coefficient) of one product term of a
        retained result: the formal master symbols, the small prefactors and
        the largest sum factor (the master's coefficient)."""
        from symbolica import AtomType

        factors = list(term) if term.get_type() == AtomType.Mul else [term]
        masters, others, coefficient = [], [], None
        for factor in factors:
            text = factor.to_canonical_string()
            if "RustRedMaster(" in text and len(text) < 4096:
                masters.append(factor)
            elif factor.get_type() == AtomType.Add and (
                coefficient is None or len(factor) > len(coefficient)
            ):
                if coefficient is not None:
                    others.append(coefficient)
                coefficient = factor
            else:
                others.append(factor)
        return masters, others, coefficient

    @classmethod
    def term_label(cls, master):
        text = master.to_canonical_string()
        label = cls.master_label(text)
        power = cls._master_power.search(text)
        return f"({label})^{power.group(1)}" if power else label

    @classmethod
    def term_selector(cls, mo, value, *, limit=500):
        """A dropdown over the top-level terms of a large sum (term -> index),
        or None for a value that is not a sum of several terms.

        The pinned Symbolica pager shows an oversized product as one numbered
        marker; opening one term at a time (``term_view``) gives its masters,
        prefactors and coefficient each their own view."""
        from symbolica import AtomType

        if value.get_type() != AtomType.Add or len(value) < 2:
            return None
        options = {}
        for index, term in enumerate(value):
            if index >= limit:
                break
            masters, _, coefficient = cls.term_parts(term)
            label = f"Term {index + 1}: " + (
                " × ".join(cls.term_label(master) for master in masters) or "no master"
            )
            if coefficient is not None:
                label += f" · coefficient of {len(coefficient):,} terms"
            options[label] = index
        shown = f"the first {limit:,} of " if len(value) > limit else ""
        return mo.ui.dropdown(
            options,
            value=next(iter(options)),
            label=f"Open one term ({shown}{len(value):,})",
            full_width=True,
        )

    @classmethod
    def term_view(cls, mo, value, index, *, page_size=25):
        """One top-level term of ``value``: its masters as labels, its small
        prefactors as a static formula and its coefficient in the live pager
        (a direct cell output, like ``formula``)."""
        import functools
        import operator

        term = value[index]
        masters, others, coefficient = cls.term_parts(term)
        parts = [
            cls.table(
                mo,
                [
                    {
                        "master": cls.term_label(master),
                        "symbol": master.to_canonical_string(),
                    }
                    for master in masters
                ],
                caption=f"Term {index + 1} of {len(value):,}: its formal masters",
            )
            if masters
            else mo.md(f"Term {index + 1} of {len(value):,} has no formal master.")
        ]
        if others:
            parts += [
                mo.md("**Prefactor**"),
                cls.static_formula(mo, functools.reduce(operator.mul, others)),
            ]
        if coefficient is not None:
            parts += [
                mo.md(f"**Coefficient** ({len(coefficient):,} terms)"),
                cls.formula(mo, coefficient, page_size=page_size),
            ]
        return mo.vstack(parts)

    @staticmethod
    def frame(mo, payload, *, height=400):
        """Isolated HTML (scripts, own styles) in an iframe that follows the
        notebook's light/dark theme instead of the operating system's."""
        mirror = (
            "<script>(function(){function sync(){try{var b=window.parent.document.body;"
            "var t=b.dataset.theme||(b.classList.contains('dark')?'dark':'light');"
            "var r=document.documentElement;r.dataset.theme=t;r.style.colorScheme=t;"
            "r.classList.remove('dark','light');r.classList.add(t)}catch(e){}}sync();"
            "try{new MutationObserver(sync).observe(window.parent.document.body,"
            "{attributes:true,attributeFilter:['class','data-theme']})}catch(e){}})();"
            "</script>"
        )
        return mo.iframe(
            mirror + Gallery.styles(frame=True) + payload, height=f"{int(height)}px"
        )

    @staticmethod
    def styles(frame=False):
        dark = (
            ':is(.dark,.dark-theme,[data-theme="dark"],[data-jp-theme-light="false"])'
        )
        light = (
            ':is(.light,.light-theme,[data-theme="light"],[data-jp-theme-light="true"])'
        )
        common = f"""
        .vg-native-formula {{ max-width:100%; overflow-x:auto; }}
        .vakint-card table {{ width:100%; }}
        .vakint-card td,.vakint-card dd {{ overflow-wrap:anywhere; }}
        /* Vakint cards follow the notebook theme, not the OS setting. */
        {light} .vakint-card {{ --vk-ink:#18364d; --vk-line:#d5e1e9;
          --vk-muted:#526777; background:#fff; }}
        {light} .vakint-card th {{ background:#edf4f8; }}
        {light} .vakint-card tbody tr:nth-child(even) {{ background:#f8fbfd; }}
        {dark} .vakint-card {{ --vk-ink:#e1edf4; --vk-muted:#aac0ce;
          --vk-line:#395361; background:#172b36; }}
        {dark} .vakint-card th {{ background:#263e4b; }}
        {dark} .vakint-card tbody tr:nth-child(even) {{ background:#1d3340; }}
        {dark} .vakint-card .vk-warning {{ background:#41341c; color:#f4d295; }}
        {dark} .vakint-card .vk-badge {{ background:#234137; color:#caefdf;
          border-color:#487563; }}
        """
        if frame:
            return (
                Gallery.math_font_css()
                + "<style>body{margin:0;font-family:system-ui,sans-serif}"
                + f"{dark} body{{color:#e1edf4}}"
                + common
                + "</style>"
            )
        return (
            Gallery.math_font_css()
            + """<style>
        :root { --vg-navy:#16324c; --vg-teal:#087e86; --vg-muted:#61758b;
          --vg-line:#dbe5eb; --vg-paper:#f4f8fa; --vg-amber:#9f5707;
          --vg-note-bg:#f0f9f7; --vg-note-ink:#234f54; --vg-note-line:#cfe5e2;
          --vg-pending-bg:#fff8ec; --vg-pending-ink:#795124;
          --vg-pending-line:#efd7b8; --vg-active-bg:#e7f5f3;
          --vg-active-ink:#086c71; --vg-track:#d3ebe8; --vg-fill:#087e86;
          --vg-good:#0b7d0b; --vg-critical:#c03535; --vg-serious:#b5532c;
          --vg-warning:#8a5b00; }
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
        .vg-nav a:hover,.vg-nav a.active { background:var(--vg-active-bg);
          color:var(--vg-active-ink); border-color:#acd9d3; }
        .vg-metrics { display:grid; grid-template-columns:repeat(auto-fit,minmax(155px,1fr));
          gap:12px; margin:1.2rem 0; }
        .vg-metric { background:var(--vg-paper); border:1px solid var(--vg-line);
          border-radius:15px; padding:20px 22px; min-width:0; }
        .vg-metric strong { display:block; color:var(--vg-navy);
          font-size:clamp(18px,2.4vw,28px); font-weight:600; letter-spacing:-.03em;
          line-height:1.3; overflow-wrap:break-word; }
        .vg-metric strong.long { font-size:clamp(15px,1.7vw,21px); }
        .vg-metric span { display:block; font-size:12px; color:var(--vg-muted);
          margin-top:7px; line-height:1.5; }
        .vg-section { border-top:1px solid var(--vg-line); padding-top:1.7rem;
          margin-top:2rem; color:var(--vg-navy); }
        .vg-section small { color:var(--vg-teal); font-size:10px;
          letter-spacing:.15em; text-transform:uppercase; }
        .vg-section h2 { font-size:24px; letter-spacing:-.025em; margin:.35rem 0; }
        .vg-section p { color:var(--vg-muted); font-size:14px; line-height:1.6; }
        .vg-note { border:1px solid var(--vg-note-line); border-left:4px solid #45a69a;
          background:var(--vg-note-bg); border-radius:10px; padding:17px 20px;
          margin:1rem 0; color:var(--vg-note-ink); font-size:14px; line-height:1.65;
          white-space:pre-line; overflow-wrap:anywhere; }
        .vg-note.pending { border-color:var(--vg-pending-line);
          border-left-color:#d99740; background:var(--vg-pending-bg);
          color:var(--vg-pending-ink); }
        .vg-note strong { display:block; margin-bottom:4px; }
        .vg-note code, .vg-legend code { font-size:12px; }
        .vg-legend { color:var(--vg-muted); font-size:12px; line-height:1.6;
          margin:.4rem 0 1rem; }
        .vg-foot { color:var(--vg-muted); border-top:1px solid var(--vg-line);
          font-size:11px; padding-top:1rem; margin-top:2rem; line-height:1.7; }
        .vg-table-wrap { max-width:100%; overflow-x:auto; margin:.8rem 0 1.2rem;
          border:1px solid var(--vg-line); border-radius:10px; }
        .vg-table { border-collapse:collapse; width:100%; font-size:13px;
          color:var(--vg-navy); }
        .vg-table th { text-align:left; font-weight:600; color:var(--vg-muted);
          background:var(--vg-paper); padding:9px 12px; font-size:12px;
          border-bottom:1px solid var(--vg-line); }
        .vg-table td { padding:9px 12px; vertical-align:top;
          border-top:1px solid var(--vg-line); overflow-wrap:anywhere;
          min-width:4.5rem; }
        .vg-table td.num { font-variant-numeric:tabular-nums; text-align:right;
          white-space:nowrap; }
        .vg-table th.num { text-align:right; }
        .vg-table caption { caption-side:bottom; text-align:left; padding:8px 12px;
          color:var(--vg-muted); font-size:12px; }
        .vg-monitor { border:1px solid var(--vg-line); border-radius:15px;
          background:var(--vg-paper); padding:18px 20px; margin:1rem 0;
          color:var(--vg-navy); }
        .vg-monitor-head { display:flex; flex-wrap:wrap; align-items:baseline;
          gap:6px 14px; }
        .vg-monitor-head strong { font-size:16px; font-weight:600; }
        .vg-monitor-head .vg-elapsed { margin-left:auto; color:var(--vg-muted);
          font-size:13px; font-variant-numeric:tabular-nums; }
        .vg-badge { display:inline-flex; gap:5px; align-items:center;
          font-size:12px; font-weight:600; padding:3px 9px; border-radius:999px;
          border:1px solid currentColor; }
        .vg-badge.running { color:var(--vg-teal); }
        .vg-badge.succeeded { color:var(--vg-good); }
        .vg-badge.failed { color:var(--vg-critical); }
        .vg-badge.interrupted,.vg-badge.stopped { color:var(--vg-serious); }
        .vg-badge.stale { color:var(--vg-warning); }
        .vg-meter { margin:14px 0 4px; }
        .vg-meter-track { height:8px; border-radius:4px; background:var(--vg-track);
          overflow:hidden; }
        .vg-meter-fill { height:100%; border-radius:4px; background:var(--vg-fill); }
        .vg-meter-label { display:flex; justify-content:space-between; gap:12px;
          font-size:12px; color:var(--vg-muted); margin-top:6px; }
        .vg-kpis { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
          gap:10px; margin-top:14px; }
        .vg-kpi { border-top:1px solid var(--vg-line); padding-top:8px; min-width:0; }
        .vg-kpi span { display:block; font-size:12px; color:var(--vg-muted); }
        .vg-kpi strong { display:block; font-size:clamp(16px,1.6vw,20px);
          font-weight:600; letter-spacing:-.02em; overflow-wrap:break-word;
          font-variant-numeric:tabular-nums; }
        .vg-kpi small { display:block; font-size:12px; color:var(--vg-muted); }
        .vg-monitor-foot { margin-top:12px; font-size:12px; color:var(--vg-muted);
          line-height:1.6; overflow-wrap:anywhere; }
        .vakint-card { max-width:100%; overflow-x:auto; }
        .vakint-card pre,.vakint-card .value { max-height:18rem; overflow:auto; }
        /* marimo's code view shows a cell whose code is hidden as a faded
           one-line strip: fold it to a thin rule so that outputs follow each
           other (hover opens the strip; a double-click or the cell menu's
           "Show code" reveals the code). The app view shows no code at all. */
        .marimo-cell [data-testid="cell-editor"].opacity-20,
        .marimo-cell .cm.opacity-20.h-8 { height:8px !important;
          min-height:0 !important; overflow:hidden; transition:height .12s; }
        .marimo-cell [data-testid="cell-editor"].opacity-20:hover,
        .marimo-cell .cm.opacity-20.h-8:hover { height:2rem !important; }
        /* Input graph -> Laurent series: the end of the pipeline at a glance. */
        .vg-pipeline { display:grid; align-items:center; gap:14px 20px;
          grid-template-columns:minmax(220px,4fr) auto minmax(380px,9fr);
          border:1px solid var(--vg-line); border-radius:15px;
          background:var(--vg-paper); padding:16px 20px; margin:1rem 0;
          color:var(--vg-navy); }
        .vg-pipe-graph { min-width:0; }
        .vg-pipe-graph iframe { width:100%; border:0; }
        .vg-pipe-arrow { display:flex; flex-direction:column; align-items:center;
          gap:8px; color:var(--vg-teal); }
        .vg-pipe-arrow b { font-size:46px; line-height:1; font-weight:400; }
        .vg-pipe-arrow small { font-size:11px; line-height:1.45; max-width:9.5rem;
          text-align:center; color:var(--vg-muted); }
        .vg-pipe-series { min-width:0; overflow-x:auto; font-size:15px; }
        .vg-pipe-series p { margin:.3rem 0; }
        .vg-pipe-label { font:600 10px/1.6 ui-sans-serif,system-ui;
          letter-spacing:.14em; text-transform:uppercase; color:var(--vg-teal); }
        .vg-pipe-caption { font-size:12px; line-height:1.6; color:var(--vg-muted);
          overflow-wrap:anywhere; }
        .vg-digits { border:1px solid var(--vg-line); border-radius:15px;
          background:var(--vg-paper); padding:14px 20px; margin:0 0 1rem;
          color:var(--vg-navy); }
        .vg-digit-row { display:grid; grid-template-columns:4.2rem minmax(0,1fr);
          gap:12px; padding:7px 0; border-top:1px solid var(--vg-line); }
        .vg-digits .vg-pipe-label + .vg-digit-row { border-top:0; }
        .vg-digit-k { font-style:italic; font-size:14px; white-space:nowrap;
          text-align:right; }
        .vg-digit-row code { background:none; padding:0; font-size:12.5px;
          line-height:1.65; white-space:normal; overflow-wrap:anywhere;
          font-variant-numeric:tabular-nums; color:inherit; }
        .vg-table td code { font-size:12px; overflow-wrap:anywhere; }
        /* One column when the card itself is narrow (sidebar open, small
           window): the series keeps the width its headline digits need. */
        .vg-pipe-wrap { container-type:inline-size; }
        @container (max-width:900px) { .vg-pipeline { grid-template-columns:1fr; }
          .vg-pipe-arrow b { transform:rotate(90deg); } }
        @media(max-width:600px) { .vg-hero {padding:1.6rem;} .vg-nav a {flex:1 1 40%;}
          .vg-eyebrow span {display:block;margin:0;} }
        /* Dark rules follow the notebook's theme (marimo body class / Jupyter
           attribute), never the operating system's preference. */
        DARK { --vg-navy:#cbddea; --vg-muted:#a9bdcc; --vg-paper:#172936;
          --vg-line:#345064; --vg-teal:#5cc6c0; --vg-note-bg:#173834;
          --vg-note-ink:#c5e7df; --vg-note-line:#2e5b55; --vg-pending-bg:#3a2e1a;
          --vg-pending-ink:#f1d3a2; --vg-pending-line:#6b5531;
          --vg-active-bg:#173834; --vg-active-ink:#9fe0da; --vg-track:#1f4547;
          --vg-fill:#4fb8b5; --vg-good:#5fd35f; --vg-critical:#ff8a8a;
          --vg-serious:#f2a07c; --vg-warning:#f6c453; }
        """.replace("DARK", dark)
            + common
            + "</style>"
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
        # Long identifiers (modular_polynomial) get a smaller face and break
        # after "_" rather than mid-word.
        long_value = ' class="long"'
        return mo.Html(
            '<div class="vg-metrics">'
            + "".join(
                '<div class="vg-metric">'
                f"<strong{long_value if len(str(value)) > 12 else ''}>"
                f"{html.escape(str(value)).replace('_', '_<wbr>')}</strong>"
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

    @staticmethod
    def legend(mo, text):
        """One small line explaining notation next to a formula."""
        return mo.Html(f'<p class="vg-legend">{html.escape(text)}</p>')

    @staticmethod
    def cell_html(value):
        """Table cell content: marimo Html (links, formulas) stays HTML,
        everything else is escaped text."""
        if hasattr(value, "text") and callable(getattr(value, "_mime_", None)):
            return value.text
        if value is None:
            return "—"
        if isinstance(value, float):
            return html.escape(f"{value:,.3f}".rstrip("0").rstrip("."))
        if isinstance(value, int) and not isinstance(value, bool):
            return f"{value:,}"
        return html.escape(str(value))

    @classmethod
    def table(cls, mo, rows, columns=None, *, caption=None, numeric=()):
        """A small read-only table that wraps long cells instead of clipping them.

        ``rows`` are dicts; ``columns`` (default: the first row's keys) fixes
        the order; ``numeric`` columns align right in tabular figures."""
        rows = list(rows)
        if columns is None:
            columns = list(rows[0]) if rows else []
        number_cell = ' class="num"'
        head = "".join(
            f'<th scope="col"{number_cell if c in numeric else ""}>'
            f"{html.escape(str(c))}</th>"
            for c in columns
        )
        body = "".join(
            "<tr>"
            + "".join(
                f"<td{number_cell if column in numeric else ''}>"
                f"{cls.cell_html(row.get(column))}</td>"
                for column in columns
            )
            + "</tr>"
            for row in rows
        )
        note = f"<caption>{html.escape(caption)}</caption>" if caption else ""
        return mo.Html(
            '<div class="vg-table-wrap"><table class="vg-table">'
            f"{note}<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"
        )

    # -- Number formatting ---------------------------------------------------
    @staticmethod
    def compact(value):
        """1,284 / 12.9K / 3.4M: auto-compact counts for metric tiles."""
        if value is None:
            return "—"
        value = float(value)
        for limit, suffix in ((1e9, "G"), (1e6, "M"), (1e4, "K")):
            if abs(value) >= limit:
                scaled = value / (limit if suffix != "K" else 1e3)
                return f"{scaled:,.1f}{suffix}"
        return (
            f"{value:,.0f}"
            if value.is_integer() or abs(value) >= 100
            else f"{value:,.1f}"
        )

    @staticmethod
    def gib(value):
        return "—" if value is None else f"{value / 1024**3:,.1f} GiB"

    @staticmethod
    def duration(seconds):
        """12.3 s below a minute, else 2:34 or 1:02:03."""
        if seconds is None:
            return "—"
        seconds = max(0.0, float(seconds))
        if seconds < 60:
            return f"{seconds:.1f} s"
        minutes, rest = divmod(round(seconds), 60)
        hours, minutes = divmod(minutes, 60)
        return f"{hours}:{minutes:02d}:{rest:02d}" if hours else f"{minutes}:{rest:02d}"

    _long_float = re.compile(r"(\d\.\d{2})\d{6,}(?:\.\.\.|…)?(e[+-]?\d+)?")

    @classmethod
    def short_numbers(cls, text):
        """Shorten printed binary floats (4.4408920985006261616945…e-16) to
        three significant digits; exact rationals and integers are untouched."""
        return cls._long_float.sub(lambda m: m.group(1) + (m.group(2) or ""), text)

    @classmethod
    def footer(cls, mo):
        runtime = cls.receipt("runtime")
        provenance = ""
        if isinstance(runtime, dict):
            provenance = "<br>" + " · ".join(
                f"{html.escape(str(label))}: {html.escape(str(value))}"
                for label, value in runtime.get("display", [])
            )
        pager = (
            "live expression pagers on"
            if cls.live_pager_available()
            else f"live expression pagers need {html.escape(cls.pager_install_hint)}"
        )
        return mo.Html(
            '<footer class="vg-foot">Native Symbolica expressions · shared HEPKit '
            f"diagrams · precision-aware Vakint records · {pager}<br>"
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
    def execution_metrics(execution):
        """Two metric tiles (arithmetic, cores) from an ``execution`` string, so
        a long arithmetic name never overflows one tile."""
        arithmetic, _, cores = str(execution).partition(" · ")
        tiles = [(arithmetic, "Ordered weight arithmetic")]
        if cores:
            count, _, unit = cores.partition(" ")
            tiles.append((count, f"worker {unit or 'cores'}"))
        return tiles

    @staticmethod
    def configured_checkout():
        """The operator's external-data checkout (VAKINT_GALLERY_ARTIFACTS), or None."""
        configured = os.environ.get("VAKINT_GALLERY_ARTIFACTS")
        return Path(configured).resolve() if configured else None

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

    @staticmethod
    def five_loop_load_enabled():
        """The in-page 1L-5L catalogue load is a server opt-in
        (VAKINT_GALLERY_ALLOW_5L=1): it runs in the page kernel for about a
        minute and holds over 10 GB for the session."""
        return os.environ.get("VAKINT_GALLERY_ALLOW_5L") == "1"

    # -- Live progress and interruption of in-kernel native calls ------------
    @staticmethod
    def supports_progress(method):
        """True when this host's Vakint call takes ``progress=`` (live progress
        and Ctrl-C / marimo Interrupt during the native call)."""
        signature = getattr(method, "__text_signature__", None)
        if not isinstance(signature, str):
            return False
        return re.search(r"\bprogress\s*=", signature) is not None

    @staticmethod
    def interrupt_options(interrupt):
        """``{"interrupt": interrupt}`` when this host's engine has the setting
        (``Vakint(interrupt="fine" | "coarse" | "off")``), else ``{}``."""
        if interrupt is None:
            return {}
        try:
            from symbolica.community.hepkit.vakint import Vakint
        except ImportError:
            return {}
        return {"interrupt": interrupt} if hasattr(Vakint, "interrupt") else {}

    @classmethod
    def progress_options(cls, method, progress):
        """``{"progress": progress}`` when the host supports it, else ``{}``.
        None keeps the host's default ("auto": marimo's own status spinner)."""
        if progress is None or not cls.supports_progress(method):
            return {}
        return {"progress": progress}

    @classmethod
    def interruptible(cls):
        """Whether in-kernel loads and reductions stop on Ctrl-C / Interrupt."""
        try:
            from symbolica.community.hepkit.vakint import Vakint
        except ImportError:
            return False
        return cls.supports_progress(Vakint.evaluate)

    # -- External data: version, catalogue counts and sizes ------------------
    _builtin_version = None

    @classmethod
    def builtin_data_version(cls, engine=None):
        """Release of the wheel's built-in 1L-3L registry (from a master's
        provenance), or None."""
        if cls._builtin_version is None:
            try:
                engine = engine or cls.engine(n_cores=1)
                master = engine.get_master(engine.master_symbols()[0], prec=2)
                version = master.provenance.version
            except (ImportError, IndexError, ValueError, AttributeError, TypeError):
                return None
            if not isinstance(version, str) or not version:
                return None
            Gallery._builtin_version = version
        return cls._builtin_version

    @classmethod
    def data_version(cls, engine=None):
        """The release a load requests. VAKINT_GALLERY_DATA_VERSION wins. A
        configured checkout is used as it is ("latest" accepts its manifest,
        offline). Downloads ask for the built-in registry's release, so the
        gallery never silently installs another published release."""
        configured = os.environ.get("VAKINT_GALLERY_DATA_VERSION")
        if configured:
            return configured
        if os.environ.get("VAKINT_GALLERY_ARTIFACTS"):
            return "latest"
        return cls.builtin_data_version(engine) or "latest"

    @classmethod
    def load_data(cls, engine, max_loop_count, *, progress=None):
        # A configured checkout is recognized and validated before any network
        # request. Otherwise the official packs install once into the cache.
        if max_loop_count >= 5 and not cls.five_loop_load_enabled():
            raise ValueError(
                "Loading the 5L catalogue is disabled on this server; the operator "
                "enables it with VAKINT_GALLERY_ALLOW_5L=1"
            )
        return engine.load_official_external_data(
            max_loop_count,
            version=cls.data_version(engine),
            local_artifacts_path=str(cls.data_directory()),
            **cls.progress_options(engine.load_official_external_data, progress),
        )

    @staticmethod
    def data_root(root=None):
        """The configured checkout, else the newest installed release in the
        download cache, else None. Never downloads."""
        if root is None:
            configured = os.environ.get("VAKINT_GALLERY_ARTIFACTS")
            cache = os.environ.get("VAKINT_GALLERY_DATA_CACHE")
            base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
            root = Path(configured or cache or Path(base) / "vakint-gallery")
        root = Path(root)
        if (root / "manifest.json").is_file():
            return root
        installed = sorted(
            (p.parent for p in root.glob(".vakint-external-data/*/manifest.json")),
            reverse=True,
        )
        return installed[0] if installed else None

    @classmethod
    def catalogue_summary(cls, root=None):
        """Release, and per loop count: masters, computed masters and unpacked
        bytes, read from the external-data manifests (JSON only); None when no
        checkout or installed release is available."""
        root = cls.data_root(root)
        if root is None:
            return None
        try:
            manifest = json.loads((root / "manifest.json").read_text())
            loops = {}
            for pack in manifest.get("loop_packs", []):
                path = root / pack["manifest"]
                if not path.is_file() or path.stat().st_size > 64 * 1024**2:
                    continue
                loop = json.loads(path.read_text())
                masters = loop.get("masters") or []
                loops[int(pack["loop_count"])] = {
                    "masters": len(masters),
                    "computed": sum(
                        1
                        for master in masters
                        if any(
                            isinstance(c, dict) and c.get("numerical")
                            for c in master.get("coefficients") or []
                        )
                    ),
                    "bytes": path.stat().st_size
                    + sum(
                        int(entry.get("bytes") or 0)
                        for entry in (loop.get("files") or {}).values()
                        if isinstance(entry, dict)
                    ),
                }
        except (OSError, ValueError, KeyError, TypeError):
            return None
        return {"version": manifest.get("version"), "path": str(root), "loops": loops}

    @staticmethod
    def size_text(size):
        """ "46 MB" (decimal, as archives are quoted), or None."""
        if not size:
            return None
        return f"{size / 1e9:.1f} GB" if size >= 1e9 else f"{size / 1e6:,.0f} MB"

    # -- Master identities as readable labels --------------------------------
    _master_symbol = re.compile(r"rustred_family_([0-9a-f]{8})[0-9a-f]*,([-\d,\s]+)\)")

    @classmethod
    def master_label(cls, symbol):
        """ "3L · family 56d4ac62 · powers (1,1,1,0,1,0)" from a master symbol.

        The arguments after the family are its propagator powers, one per
        loop-momentum scalar product, so L(L+1)/2 of them give the loop count."""
        text = str(getattr(symbol, "expression", symbol))
        match = cls._master_symbol.search(text)
        if match is None:
            return text
        powers = [p.strip() for p in match.group(2).split(",") if p.strip()]
        loops = (int((8 * len(powers) + 1) ** 0.5) - 1) // 2
        loops = f"{loops}L · " if loops * (loops + 1) // 2 == len(powers) else ""
        return f"{loops}family {match.group(1)} · powers ({','.join(powers)})"

    @classmethod
    def master_options(cls, symbols):
        """Dropdown options {readable label: symbol text}, by loops then family."""
        options = {}
        for symbol in symbols:
            text = str(symbol.expression)
            label = cls.master_label(text)
            while label in options:
                label += " ·"
            options[label] = text
        return dict(sorted(options.items()))

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
            recorded_with = (receipt.get("provenance") or {}).get("data_version")
            return mo.vstack(
                [
                    cls.note(
                        mo,
                        "Recorded coverage advisory · reduction continued"
                        + (
                            f" · recorded with data {recorded_with}"
                            if recorded_with
                            else ""
                        ),
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
                parts.append(cls.table(mo, timing_rows, numeric=("Seconds",)))
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
    def ladder_job(cls, name):
        """The ladder job named ``name`` (its directory name, e.g. from the page
        URL's ``?job=``) under this server's work root, running or finished;
        None for anything else."""
        if not isinstance(name, str) or not re.fullmatch(r"\w[\w.-]*", name):
            return None
        work = cls.work_root() / name
        if (work / "status.json").is_file() and gallery_jobs.job_kind(work) == "ladder":
            return str(work)
        return None

    @classmethod
    def job_link(cls, mo, name, label="show"):
        """A link that reopens the ladder page on job ``name`` (``?job=``).
        marimo's edit server needs the notebook file in the URL as well."""
        if cls.mode(mo) == "edit" and not os.environ.get("VAKINT_GALLERY_URL_PREFIX"):
            target = cls.edit_link("05_ladder_reductions.py", job=name)
        else:
            target = f"?job={quote(name)}"
        return mo.Html(
            f'<a href="{html.escape(target)}" target="_self">{html.escape(label)}</a>'
        )

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
        live["monitor"] = cls.job_monitor(work)
        if isinstance(current.get("worker_pid"), int):
            live["resources"] = benchmark_job.worker_resources(current["worker_pid"])
        return live

    @classmethod
    def cancel_job(cls, path):
        """Stop every process of a job (SIGTERM, then SIGKILL) and confirm it."""
        return gallery_jobs.cancel(cls.checked_job(path))

    @classmethod
    def wait_job(
        cls,
        path,
        *,
        timeout=None,
        cancel_after=None,
        poll_seconds=0.5,
        cancel_on_interrupt=False,
    ):
        """Block until a job ends (script mode). ``cancel_after`` seconds cancels
        it exactly as the Cancel button does; the report joins the status. With
        ``cancel_on_interrupt``, Ctrl-C cancels the job too (status
        ``interrupted``) instead of leaving it running until its deadline."""
        started = time.monotonic()
        report = None
        try:
            while True:
                current = cls.job_status(path)
                if current.get("state") in gallery_jobs.FINAL_STATES:
                    return (
                        current
                        if report is None
                        else {**current, "cancel_report": report}
                    )
                elapsed = time.monotonic() - started
                if (
                    cancel_after is not None
                    and report is None
                    and elapsed >= cancel_after
                ):
                    report = cls.cancel_job(path)
                    continue
                if timeout is not None and elapsed >= timeout:
                    return current
                time.sleep(poll_seconds)
        except KeyboardInterrupt:
            if not cancel_on_interrupt:
                raise
            report = cls.cancel_job(path)
            return {
                **cls.job_status(path),
                "cancel_report": report,
                "interrupted": True,
            }

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

    # -- Progress monitor card (background jobs and in-kernel calls) ---------
    monitor_schema = "vakint.monitor.v2"
    stale_seconds = 10.0
    stage_labels: ClassVar[dict[str, str]] = {
        "setup": "Setup",
        "data_load": "Data load",
        "ingress": "Ingress",
        "lowering": "Scalar lowering",
        "native": "Native reduction",
        "post_processing": "Post-processing",
        "output": "Output",
        "other": "Other",
    }
    # Job states and native statuses -> (badge class, icon, label).
    badges: ClassVar[dict[str, tuple[str, str, str]]] = {
        "running": ("running", "●", "Running"),
        "stopping": ("stopped", "■", "Stopping"),
        "complete": ("succeeded", "✓", "Complete"),
        "succeeded": ("succeeded", "✓", "Succeeded"),
        "cancelled": ("interrupted", "■", "Cancelled"),
        "interrupted": ("interrupted", "■", "Interrupted"),
        "timed_out": ("stopped", "■", "Timed out"),
        "failed": ("failed", "✕", "Failed"),
        "partial": ("failed", "✕", "Partial result"),
        "lost": ("failed", "✕", "Lost"),
    }

    @classmethod
    def job_monitor(cls, path):
        """A job's native monitor model (``monitor.json``, vakint.monitor.v2),
        read whole (the native publisher bounds it), with the file's age as
        ``snapshot_age_seconds``; None before the first snapshot or for any
        other document."""
        target = cls.checked_job(path) / "monitor.json"
        try:
            with target.open("rb") as stream:
                payload = stream.read(cls.display_receipt_bytes + 1)
                modified = os.fstat(stream.fileno()).st_mtime
            if len(payload) > cls.display_receipt_bytes:
                return None
            document = json.loads(payload)
        except (OSError, ValueError, RecursionError):
            return None
        monitor = document.get("monitor") if isinstance(document, dict) else None
        if not isinstance(monitor, dict) or monitor.get("schema") != cls.monitor_schema:
            return None
        model = monitor.get("model")
        if not isinstance(model, dict):
            return None
        return {**model, "snapshot_age_seconds": max(0.0, time.time() - modified)}

    @staticmethod
    def progress_model(progress):
        """A monitor model dict from a ``VakintProgress`` snapshot (or a dict).

        Native evaluations carry the full vakint.monitor.v2 model
        (``to_dict()``); data loads only report state, stage and elapsed time."""
        if progress is None or isinstance(progress, dict):
            return progress or {}
        to_dict = getattr(progress, "to_dict", None)
        model = to_dict() if callable(to_dict) else {}
        status = getattr(progress, "status", None)
        if model:
            return {**model, "status": status or model.get("status")}

        def get(name):
            return getattr(progress, name, None)

        return {
            "status": status,
            "elapsed_seconds": get("elapsed_seconds"),
            "stage": {"kind": get("stage"), "activity": get("activity")},
            "memory": {
                "rss_bytes": get("rss_bytes"),
                "peak_rss_bytes": get("peak_rss_bytes"),
            },
        }

    @staticmethod
    def _dict(value, *keys):
        for key in keys:
            value = value.get(key) if isinstance(value, dict) else None
        return value if isinstance(value, dict) else {}

    @staticmethod
    def _number(value):
        return (
            value
            if isinstance(value, (int, float)) and not isinstance(value, bool)
            else None
        )

    @classmethod
    def monitor_card(
        cls,
        mo,
        model,
        *,
        state=None,
        title=None,
        elapsed=None,
        activity=None,
        rss=None,
        peak_rss=None,
    ):
        """The simplified Vakint reduction monitor: the dashboard's most
        important metrics in one card, for background jobs and in-kernel calls.

        ``model`` is a vakint.monitor.v2 model (``monitor.json``'s
        ``monitor.model`` or ``VakintProgress.to_dict()``), possibly partial or
        None. ``state`` (a job state) overrides the native status; ``elapsed``,
        ``activity`` and the RSS values fill in what the model lacks (a job's
        coarse stage before the native monitor starts, the worker's memory)."""
        model = model if isinstance(model, dict) else {}
        number = cls._number
        status = str(state or model.get("status") or "running").lower()
        age = number(model.get("snapshot_age_seconds"))
        badge, icon, label = cls.badges.get(status, ("running", "●", status.title()))
        if badge == "running" and age is not None and age > cls.stale_seconds:
            badge, icon, label = "stale", "⚠", f"No update for {cls.duration(age)}"
        stage = cls._dict(model, "stage")
        ordered = cls._dict(model, "ordered")
        phase = cls._dict(ordered, "phase")
        kind = stage.get("kind") if isinstance(stage.get("kind"), str) else None
        headline = cls.stage_labels.get(kind, "") if kind else ""
        detail = stage.get("activity") or (
            stage.get("id").replace("_", " ")
            if isinstance(stage.get("id"), str) and stage.get("id") != kind
            else None
        )
        total = number(phase.get("total"))
        index = number(phase.get("index"))
        phase_kind = phase.get("kind") if isinstance(phase.get("kind"), str) else None
        if kind == "native" and total:
            if phase_kind == "done" or (index is not None and index >= total):
                detail = "Ordered schedule complete"
            else:
                parts = [f"{(phase_kind or 'reduce').title()} phase"]
                if number(phase.get("support")) is not None:
                    parts.append(f"support {phase['support']}")
                if number(phase.get("owner")) is not None:
                    parts.append(f"owner {phase['owner']}")
                detail = " · ".join(parts)
        if activity and (not headline or kind == "setup"):
            headline, detail = activity, None
        headline = " · ".join(str(p) for p in (headline, detail) if p) or "Starting"
        elapsed = number(model.get("elapsed_seconds")) if elapsed is None else elapsed

        parts = [
            '<div class="vg-monitor">',
            '<div class="vg-monitor-head">',
            (
                f'<span class="vg-badge {badge}"><span aria-hidden="true">{icon}</span>'
                f"{html.escape(label)}</span>"
            ),
            (
                f"<strong>{html.escape(title + ' · ' if title else '')}"
                f"{html.escape(headline)}</strong>"
            ),
            f'<span class="vg-elapsed">{html.escape(cls.duration(elapsed))} elapsed</span>'
            if elapsed is not None
            else "",
            "</div>",
        ]
        if total:
            done = min(int(index or 0), int(total))
            parts.append(
                '<div class="vg-meter">'
                f'<div class="vg-meter-track" role="progressbar" aria-valuemin="0" '
                f'aria-valuemax="{int(total)}" aria-valuenow="{done}" '
                f'aria-label="Ordered schedule position">'
                f'<div class="vg-meter-fill" style="width:{100 * done / total:.1f}%">'
                "</div></div>"
                f'<div class="vg-meter-label"><span>Ordered schedule: {done} / '
                f"{int(total)} phases done</span><span>a position in the schedule, "
                "not a share of the time</span></div></div>"
            )
        tiles = []
        rate = cls._dict(ordered, "recent_rate") or cls._dict(ordered, "request_rate")
        if number(rate.get("keys")) is not None and number(rate.get("seconds")):
            tiles.append(
                (
                    "Keys per second",
                    cls.compact(rate["keys"] / rate["seconds"]),
                    "recent"
                    if cls._dict(ordered, "recent_rate")
                    else "request average",
                )
            )
        if number(ordered.get("processed_keys")) is not None:
            inputs = number(ordered.get("input_keys"))
            tiles.append(
                (
                    "Processed keys",
                    cls.compact(ordered["processed_keys"]),
                    f"{inputs:,} input keys" if inputs is not None else "",
                )
            )
        frontier = cls._dict(ordered, "frontier")
        if number(frontier.get("keys")) is not None:
            peak = number(frontier.get("peak_keys"))
            tiles.append(
                (
                    "Frontier keys",
                    cls.compact(frontier["keys"]),
                    f"peak {cls.compact(peak)}" if peak is not None else "",
                )
            )
        hits = number(ordered.get("step_cache_hits"))
        misses = number(ordered.get("step_cache_misses"))
        if hits is not None and misses is not None:
            looked_up = hits + misses
            tiles.append(
                (
                    "Step cache",
                    f"{cls.compact(hits)} / {cls.compact(misses)}",
                    "hits / misses in the Ordered census"
                    + (f" · {hits / looked_up:.0%} hits" if looked_up else ""),
                )
            )
        memory = cls._dict(model, "memory")
        rss = number(memory.get("rss_bytes")) if rss is None else rss
        peak_rss = (
            number(memory.get("peak_rss_bytes")) if peak_rss is None else peak_rss
        )
        if rss is not None:
            cpu = number(memory.get("cpu_percent"))
            tiles.append(
                (
                    "Memory (RSS)",
                    cls.gib(rss),
                    " · ".join(
                        part
                        for part in (
                            f"peak {cls.gib(peak_rss)}" if peak_rss else "",
                            f"{cpu:.0f}% CPU" if cpu is not None else "",
                        )
                        if part
                    ),
                )
            )
        batch = cls._dict(model, "batch")
        if number(batch.get("total")):
            tiles.append(
                (
                    "Batch",
                    f"{batch.get('completed', 0)} / {batch['total']}",
                    str(batch.get("label") or ""),
                )
            )
        if tiles:
            parts.append(
                '<div class="vg-kpis">'
                + "".join(
                    f'<div class="vg-kpi"><span>{html.escape(name)}</span>'
                    f"<strong>{html.escape(value)}</strong>"
                    f"<small>{html.escape(note)}</small></div>"
                    for name, value, note in tiles
                )
                + "</div>"
            )
        notes = []
        census = number(ordered.get("census_age_seconds"))
        if badge in ("running", "stale") and kind == "native":
            if total and index is not None and index >= total:
                notes.append(
                    "The schedule is done; the modular lift and finalization that "
                    "follow report no live census."
                )
            elif census is not None and census > cls.stale_seconds:
                notes.append(f"Last Ordered census {cls.duration(census)} ago.")
        run = cls._dict(model, "run")
        if run:
            cores = number(run.get("n_cores"))
            notes.append(
                " · ".join(
                    str(part)
                    for part in (
                        run.get("arithmetic"),
                        f"{cores} core{'' if cores == 1 else 's'}" if cores else None,
                        f"step cache {run['step_cache']}"
                        if run.get("step_cache")
                        else None,
                    )
                    if part
                )
            )
        ledger = [
            f"{cls.stage_labels.get(item.get('kind'), item.get('kind'))} "
            f"{cls.duration(item['seconds'])}"
            for item in model.get("stages") or []
            if isinstance(item, dict)
            and number(item.get("seconds"))
            and item["seconds"] >= 0.05
        ]
        if ledger:
            notes.append("Stages: " + " · ".join(ledger))
        if notes:
            parts.append(
                '<div class="vg-monitor-foot">'
                + "<br>".join(html.escape(n) for n in notes if n)
                + "</div>"
            )
        parts.append("</div>")
        return mo.Html("".join(parts))

    @classmethod
    def job_monitor_card(cls, mo, status, path=None, *, title=None):
        """The monitor card of a background job: its state, coarse stage and
        worker memory from the supervisor, the rest from the native model."""
        model = status.get("monitor")
        if model is None and path is not None:
            model = cls.job_monitor(path)
        stage = status.get("stage") or {}
        native_started = bool(cls._dict(model, "stage").get("kind")) and (
            cls._dict(model, "stage").get("kind") != "setup"
        )
        resources = status.get("resources") or {}
        state = str(status.get("state", "running"))
        elapsed = cls._number(status.get("elapsed_seconds"))
        updated = cls._number(status.get("updated_unix_seconds"))
        if elapsed is not None and updated is not None and state == "running":
            # The supervisor rewrites its status every five seconds.
            elapsed += max(0.0, time.time() - updated)
        return cls.monitor_card(
            mo,
            model,
            state=state,
            title=title,
            elapsed=elapsed,
            activity=None
            if native_started
            else (stage.get("stage") if isinstance(stage.get("stage"), str) else None),
            rss=cls._number(resources.get("worker_rss_bytes")),
            peak_rss=cls._number(resources.get("worker_peak_rss_bytes")),
        )

    @classmethod
    def kernel_progress(cls, mo, *, title=None):
        """A ``progress=`` callback that redraws this cell's output as the
        monitor card, for in-kernel calls (edit mode). Its ``model`` attribute
        keeps the last snapshot (the final one, after the call)."""

        def show(progress):
            show.model = cls.progress_model(progress)
            mo.output.replace(cls.monitor_card(mo, show.model, title=title))

        show.model = None
        return show

    @classmethod
    def kernel_outcome(
        cls, mo, work, *, progress=None, interrupted=False, title="In this kernel"
    ):
        """The final card and a note for an in-kernel run, as the output of
        the cell that ran it.

        marimo runs no other cell after an Interrupt (not even the cells a
        state change would re-run), so an interrupted run must show itself
        here; ``progress`` is the run's ``kernel_progress`` callback."""
        status = cls.job_status(work)
        state = str(status.get("state") or ("cancelled" if interrupted else "failed"))
        elapsed = cls._number(status.get("elapsed_seconds"))
        model = getattr(progress, "model", None)
        if model:
            card = cls.monitor_card(
                mo,
                model,
                state="interrupted" if interrupted else state,
                title=title,
                elapsed=elapsed,
            )
        else:
            card = cls.job_monitor_card(mo, status, work, title=title)
        model = model or {}
        phase = cls._dict(model, "ordered", "phase")
        kind = cls._dict(model, "stage").get("kind")
        where = cls.stage_labels.get(kind) if isinstance(kind, str) else None
        if kind == "native" and cls._number(phase.get("total")):
            where = (
                f"{where} · {str(phase.get('kind') or 'reduce').title()} "
                f"{phase.get('index', 0)}/{phase['total']}"
            )
        taken = f" after {cls.duration(elapsed)}" if elapsed is not None else ""
        if interrupted:
            note = cls.note(
                mo,
                "Interrupted in this kernel",
                f"Stopped{taken}{f' during {where}' if where else ''}. The kernel "
                "and its engine stay usable: start the run again or pick another "
                "case. marimo re-runs no other cell after an Interrupt, so "
                "04 / progress and 05 / result below still show the run before "
                "this one until the next run or page change.",
            )
        elif state == "complete":
            note = cls.note(
                mo,
                "Finished in this kernel",
                f"Complete{taken}. 04 / progress and 05 / result below show this run.",
            )
        else:
            note = cls.note(
                mo,
                "The in-kernel run failed",
                str(status.get("message") or "See the kernel console for the error."),
            )
        return mo.vstack([card, note])

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

    # -- Graph drawings and DOT source ----------------------------------------
    @staticmethod
    def dot_source(mo, dot):
        """Read-only DOT source as plain text. marimo's code editor has no DOT
        grammar: declaring language="dot" shows "Language dot not supported"."""
        return mo.ui.code_editor(dot, language="text", disabled=True)

    @staticmethod
    def diagram_height(loop_count):
        return {1: 260, 2: 380, 3: 380}.get(int(loop_count), 480)

    @classmethod
    def diagram_frame(cls, mo, diagram, *, height=None):
        """The interactive FeynKit drawing of a native diagram, with its loop
        momenta, in a theme-following iframe (page 01's view). Its script
        never touches the notebook's DOM, so a new selection simply redraws."""
        return cls.frame(
            mo,
            diagram.render(momenta=True).to_html(),
            height=height or cls.diagram_height(diagram.loop_count),
        )

    _ladder_diagrams: ClassVar[dict] = {}

    @classmethod
    def ladder_diagram(cls, loop_count):
        """The native QCD ``FeynmanDiagram`` of ``ladder_dot(loop_count)``."""
        loop_count = int(loop_count)
        if loop_count not in cls._ladder_diagrams:
            from symbolica.community.hepkit import FeynmanDiagram, Model

            cls._ladder_diagrams[loop_count] = FeynmanDiagram.from_dot(
                Model.qcd(), cls.ladder_dot(loop_count)
            )
        return cls._ladder_diagrams[loop_count]

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

    @classmethod
    def kernel_work(cls, case, *, n_cores=None, arithmetic=None):
        """A fresh run directory for an in-kernel reduction of ``case``."""
        cases = cls.ladder_cases()
        if case not in cases:
            raise ValueError(f"Unknown ladder case: {case}")
        cls.ladder_input(case)
        arithmetic = arithmetic or "modular_polynomial"
        if arithmetic not in cls.arithmetics:
            raise ValueError(f"Unknown Ordered weight arithmetic: {arithmetic}")
        stamp = time.strftime("%Y%m%dT%H%M%S", time.gmtime())
        work = cls.work_root() / "kernel" / f"{case}-{stamp}-{os.getpid()}"
        work.mkdir(parents=True, exist_ok=False)
        spec = {
            "case": case,
            "n_cores": cls.bounded_cores(n_cores),
            "arithmetic": arithmetic,
            "data_directory": str(cls.data_directory())
            if cases[case]["loop_count"] >= 4 or cls.configured_checkout()
            else None,
            # In memory: an in-kernel run never writes the jobs' persistent cache.
            "ordered_cache_directory": None,
        }
        (work / "request.json").write_text(json.dumps({"kind": "kernel", "spec": spec}))
        return str(work)

    @classmethod
    def reduce_in_kernel(cls, work, *, progress=None):
        """Reduce a ``kernel_work`` run in this kernel (marimo edit mode).

        Interrupt (or Ctrl-C in a terminal) stops the native call: the run is
        recorded as cancelled and the KeyboardInterrupt propagates. The final
        ``status.json`` makes the run look like a finished job to the page."""
        import reduce_worker

        work = cls.checked_job(work)
        spec = json.loads((work / "request.json").read_text())["spec"]
        started = time.monotonic()
        state, message = "failed", None
        try:
            state = (
                "complete"
                if reduce_worker.run(work, spec, progress=progress) == 0
                else "failed"
            )
        except KeyboardInterrupt:
            state, message = "cancelled", "Interrupted in this kernel"
            raise
        finally:
            (work / "status.json").write_text(
                json.dumps(
                    {
                        "state": state,
                        "kind": "kernel",
                        "elapsed_seconds": time.monotonic() - started,
                        "updated_unix_seconds": time.time(),
                        "interrupted": state == "cancelled",
                        **({"message": message} if message else {}),
                    }
                )
            )
        return str(work)

    # -- Input graph -> Laurent series (this kernel) ---------------------------
    laurent_max_digits = 500
    laurent_highest_order = 2
    # Digits per coefficient in the headline series (it stays one readable
    # line); more requested digits are listed in full below it.
    laurent_headline_digits = 32
    # Exact parameter values: at most this many digits and |exponent|.
    laurent_exact_limit = 1000
    # Without marimo's Interrupt (application mode) the page evaluates in its
    # kernel only cases whose recorded data load and reduction take at most this.
    laurent_app_mode_seconds = 15.0
    # Engines kept per session; a 4L one holds the 4L pack (about 2.5 GB).
    laurent_engine_slots = 2
    # Known parameters: (plain label, LaTeX, meaning); others show their name.
    laurent_labels: ClassVar[dict[str, tuple[str, str, str]]] = {
        "vakint::mUV": (
            "m_UV",
            r"m_{\mathrm{UV}}",
            "common UV mass of every vacuum propagator",
        ),
        "vakint::mursq": ("μ_R²", r"\mu_R^2", "MSbar renormalization scale squared"),
        "UFO::G": (
            "G",
            "G",
            "strong coupling of the UFO model; the QCD vertices carry GC_11 = i·G",
        ),
    }
    _ladder_sources: ClassVar[dict] = {}
    _laurent_symbolic: ClassVar[dict] = {}
    _laurent_engines: ClassVar[list] = []
    _number_words: ClassVar[dict[int, str]] = {
        1: "one",
        2: "two",
        3: "three",
        4: "four",
        5: "five",
    }

    @classmethod
    def ladder_source(cls, case):
        """``case``'s prepared input with d = 4 − 2ε bound, verified like a run
        (file, canonical text and bound input sha256); parsed once a session."""
        if case not in cls._ladder_sources:
            from symbolica import E

            spec = cls.ladder_cases().get(case)
            if spec is None:
                raise ValueError(f"Unknown ladder case: {case}")
            path = cls.ladder_input(case)
            payload = path.read_bytes()
            if hashlib.sha256(payload).hexdigest() != spec["input"]["sha256"]:
                raise ValueError(f"{path.name}: sha256 differs from the recorded input")
            source = E(payload.decode(), default_namespace="qcd_ladder")
            canonical = source.to_canonical_string().encode()
            if hashlib.sha256(canonical).hexdigest() != spec["input"]["sha256"]:
                raise ValueError(
                    "The parsed input does not reproduce its canonical text"
                )
            binding = spec["dimension_binding"]
            bound = source.replace(E(binding["source"]), E(binding["target"]))
            digest = hashlib.sha256(bound.to_canonical_string().encode()).hexdigest()
            if digest != spec["evaluated_input_sha256"]:
                raise ValueError(
                    "The dimension-bound input differs from the recorded one"
                )
            cls._ladder_sources[case] = bound
        return cls._ladder_sources[case]

    @classmethod
    def laurent_parameters(cls, case):
        """Every parameter a numerical Laurent series of ``case`` needs, in
        page order: the input's own symbols other than ε (a UFO coupling
        replaced by the model parameters it is made of) and μ_R², which Vakint's
        MSbar normalization adds. Dicts: name, key (for numerical_evaluation),
        label, latex, meaning."""
        from symbolica import E
        from symbolica.community.hepkit import Model

        model = Model.qcd()
        names = set()
        for symbol in cls.ladder_source(case).get_all_symbols(False):
            name = symbol.get_name()
            if name == "vakint::ε" or name.startswith("symbolica::"):
                continue
            if name.startswith("UFO::GC_"):
                names.update(
                    part.get_name()
                    for part in model.expand_couplings(E(name)).get_all_symbols(False)
                    if not part.get_name().startswith("symbolica::")
                )
            else:
                names.add(name)
        names.add("vakint::mursq")
        known = list(cls.laurent_labels)
        rows = []
        for name in sorted(
            names, key=lambda n: (known.index(n) if n in known else len(known), n)
        ):
            short = name.rsplit("::", 1)[-1]
            label, latex, meaning = cls.laurent_labels.get(
                name,
                (
                    short,
                    r"\mathtt{" + short.replace("_", r"\_") + "}",
                    "parameter of the input",
                ),
            )
            rows.append(
                {
                    "name": name,
                    "key": name.removeprefix("vakint::"),
                    "label": label,
                    "latex": latex,
                    "meaning": meaning,
                }
            )
        return rows

    _exponent = re.compile(r"[eE]\s*([+-]?[\d_]+)\s*$")

    @classmethod
    def exact_number(cls, text):
        """An exact rational from a decimal or a fraction ("91.1876", "1/3",
        "2.5e-3"), so that a parameter never limits the precision. At most
        ``laurent_exact_limit`` digits and |exponent| (checked before the
        conversion, which would otherwise build the full power of ten)."""
        text = str(text).strip()
        limit = cls.laurent_exact_limit
        too_large = ValueError(
            f"“{text[:40]}{'…' if len(text) > 40 else ''}” is too large: use at "
            f"most {limit} digits and an exponent between −{limit} and {limit}"
        )
        if len(text) > limit + 40:
            raise too_large
        exponent = cls._exponent.search(text)
        try:
            power = abs(int(exponent.group(1))) if exponent else 0
        except ValueError:
            power = 0  # malformed: Fraction reports it below
        if power > limit:
            raise too_large
        try:
            value = Fraction(text)
        except (ValueError, ZeroDivisionError):
            raise ValueError(
                f"“{text}” is not an exact decimal or fraction (e.g. 91.1876 or 1/3)"
            ) from None
        # About 3.33 bits per decimal digit (a margin on top of the limit).
        bits = int(limit * 3.33) + 8
        if max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > (
            2 * bits
        ):
            raise too_large
        return value

    @staticmethod
    def exact_expression(value):
        """A Fraction as an exact Symbolica rational (substituted exactly by
        ``numerical_evaluation``; a Python float would cap it at 53 bits)."""
        from symbolica import E

        return E(f"{value.numerator}/{value.denominator}")

    @classmethod
    def laurent_values(cls, parameters, values):
        """``{name: Fraction}`` of every parameter from ``values`` ({name:
        text}); ValueError names the first missing, malformed or non-positive
        (masses, scales) one."""
        exact = {}
        for parameter in parameters:
            text = (values or {}).get(parameter["name"])
            if text is None or not str(text).strip():
                raise ValueError(f"Set a value for {parameter['label']}")
            try:
                value = cls.exact_number(text)
            except ValueError as error:
                raise ValueError(f"{parameter['label']}: {error}") from None
            if parameter["name"] in ("vakint::mUV", "vakint::mursq") and value <= 0:
                raise ValueError(f"{parameter['label']} must be positive")
            exact[parameter["name"]] = value
        return exact

    @classmethod
    def laurent_in_page(cls, mo, case):
        """Whether this page evaluates ``case`` in its kernel: always in edit
        mode (marimo's Interrupt stops it there), otherwise only cases whose
        recorded data load and reduction take at most laurent_app_mode_seconds."""
        if cls.mode(mo) == "edit":
            return True
        recorded = cls.ladder_cases()[case]["recorded"]
        seconds = (recorded.get("data_loading_seconds") or 0) + (
            recorded.get("retained_reduction_seconds") or 0
        )
        return seconds <= cls.laurent_app_mode_seconds

    @classmethod
    def laurent_form(cls, mo, case):
        """The Laurent panel's controls for ``case`` (a marimo batch): one exact
        text field per parameter (``p0``, ``p1``, … in ``laurent_parameters``
        order), the requested digits and the ε range [−L, highest order].

        The page builds one per case, once, in a cell that never re-runs:
        marimo replays the last value of each element id to a page that
        reconnects, and an element re-created by a re-run takes the same id
        with its fresh value, so the page would show values the kernel does
        not use."""
        loops = cls.ladder_cases()[case]["loop_count"]
        fields = {
            f"p{index}": mo.ui.text(
                value="1",
                label=f"${parameter['latex']}$ · {parameter['meaning']}",
                full_width=True,
            )
            for index, parameter in enumerate(cls.laurent_parameters(case))
        }
        return mo.md(
            "**Parameters**, exact decimals or fractions (e.g. `91.1876`, `1/3`)\n\n"
            + "\n\n".join(f"{{{key}}}" for key in fields)
            + "\n\n**Accuracy and ε orders**\n\n{digits}\n\n{orders}"
        ).batch(
            digits=mo.ui.number(
                start=1,
                stop=cls.laurent_max_digits,
                value=32,
                step=1,
                label="Requested decimal digits",
            ),
            orders=mo.ui.range_slider(
                start=-loops,
                stop=cls.laurent_highest_order,
                step=1,
                value=[-loops, 0],
                show_value=True,
                label="ε orders (lowest, highest)",
            ),
            **fields,
        )

    @classmethod
    def laurent_request(cls, case, value):
        """``laurent_series`` keyword arguments (values, digits, orders) from
        the value of ``case``'s ``laurent_form``."""
        value = value or {}
        return {
            "values": {
                parameter["name"]: value.get(f"p{index}")
                for index, parameter in enumerate(cls.laurent_parameters(case))
            },
            "digits": value.get("digits"),
            "orders": value.get("orders"),
        }

    @classmethod
    def laurent_digits(cls, digits):
        """The requested decimal digits as an int from 1 to laurent_max_digits;
        otherwise a ValueError that says what to set (an empty field is None)."""
        message = (
            "Set the requested decimal digits: a whole number from 1 to "
            f"{cls.laurent_max_digits}"
        )
        if digits is None or isinstance(digits, bool):
            raise ValueError(message)
        try:
            number = Fraction(str(digits).strip())
        except (ValueError, ZeroDivisionError):
            raise ValueError(message) from None
        if number.denominator != 1 or not 1 <= number <= cls.laurent_max_digits:
            raise ValueError(message)
        return int(number)

    @classmethod
    def laurent_data_loops(cls, loops):
        """Loops of the external data a case's engine loads, 0 for the
        built-in registry: as for a run (4L always, 3L from a configured checkout)."""
        return loops if loops >= 4 or cls.configured_checkout() else 0

    @classmethod
    def laurent_engine(
        cls, loops, *, precision, terms=None, progress=None, timings=None
    ):
        """A symbolic-master MSbar engine slot for ``loops``-loop inputs (a
        dict: data, digits, terms, engine, version): the most recent one
        holding the same data with ``precision`` digits and ``terms`` ε orders
        (None: any), else a new one (its data loaded with ``progress``; the
        seconds of the construction and data load are added to
        ``timings["load"]``). At most ``laurent_engine_slots`` engines stay
        alive."""
        data = cls.laurent_data_loops(loops)
        for index in range(len(cls._laurent_engines) - 1, -1, -1):
            slot = cls._laurent_engines[index]
            if (
                slot["data"] == data
                and slot["digits"] == int(precision)
                and terms in (None, slot["terms"])
            ):
                cls._laurent_engines.append(cls._laurent_engines.pop(index))
                return slot
        # Free a slot first: a 4L pack is large.
        surplus = len(cls._laurent_engines) - cls.laurent_engine_slots + 1
        del cls._laurent_engines[: max(0, surplus)]
        count = int(terms or loops + 1)
        started = time.perf_counter()
        engine = cls.engine(
            precision=int(precision),
            mode="symbolic",
            normalization="MSbar",
            use_dot_product_notation=True,
            number_of_terms_in_epsilon_expansion=count,
        )
        version = "built-in 1L–3L registry"
        if data:
            report = cls.load_data(engine, data, progress=progress)
            version = f"{data}L external data {report.get('version') or ''}".strip()
        if timings is not None:
            timings["load"] = timings.get("load", 0.0) + time.perf_counter() - started
        slot = {
            "data": data,
            "digits": int(precision),
            "terms": count,
            "engine": engine,
            "version": version,
        }
        cls._laurent_engines.append(slot)
        return slot

    @staticmethod
    def missing_orders(error):
        """The content of a MissingMasterOrdersError as plain data."""
        masters = [
            {
                "master": str(item[0]),
                "first_missing": item[1],
                "last_needed": item[2],
                "first_affected": item[3],
            }
            for item in getattr(error, "masters", None) or []
            if isinstance(item, (tuple, list)) and len(item) == 4
        ]
        affected = [
            m["first_affected"] for m in masters if isinstance(m["first_affected"], int)
        ]
        return {
            "message": str(error).removeprefix("Vakint error | "),
            "backend": getattr(error, "backend", None),
            "requested_order": getattr(error, "requested_order", None),
            "masters": masters,
            "first_affected": min(affected) if affected else None,
        }

    @classmethod
    def missing_orders_text(cls, missing):
        """What a fail-closed order lacks, e.g. "it needs ε² of 4 four-loop
        masters, which RustRed, data v0.2.1 does not provide"."""
        groups = {}
        for master in missing.get("masters") or []:
            loops = re.match(r"L(\d+)_", master["master"])
            word = cls._number_words.get(int(loops.group(1))) if loops else None
            key = (master["first_missing"], f"{word}-loop" if word else "")
            groups[key] = groups.get(key, 0) + 1
        if not groups:
            return "it needs master Laurent coefficients beyond the available data"
        parts = [
            f"{cls.epsilon_text(order)} of {count} {kind + ' ' if kind else ''}"
            f"master{'s' if count != 1 else ''}"
            for (order, kind), count in sorted(groups.items())
        ]
        backend = missing.get("backend") or "the loaded data"
        return f"it needs {', '.join(parts)}, which {backend} does not provide"

    @classmethod
    def laurent_symbolic(cls, case, terms, *, precision, progress=None, timings=None):
        """``((expression, missing), cached)``: ``case``'s whole input evaluated
        by Vakint with symbolic masters (Ordered reduction; masters as exact
        constants and data coefficients) through ``terms`` ε orders, couplings
        expanded by the QCD model; or ``missing`` (MissingMasterOrders as
        plain data) when Vakint fails closed. Kept per input and order count;
        an interrupted evaluation keeps nothing. ``timings`` as for
        ``laurent_engine``."""
        key = (case, int(terms))
        if key in cls._laurent_symbolic:
            return cls._laurent_symbolic[key], True
        from symbolica.community.hepkit import Model
        from symbolica.community.hepkit.vakint import MissingMasterOrdersError

        loops = cls.ladder_cases()[case]["loop_count"]
        source = cls.ladder_source(case)
        engine = cls.laurent_engine(
            loops,
            precision=precision,
            terms=int(terms),
            progress=progress,
            timings=timings,
        )["engine"]
        try:
            value = engine.evaluate(
                source, **cls.progress_options(engine.evaluate, progress)
            )
            entry = (Model.qcd().expand_couplings(value), None)
        except MissingMasterOrdersError as error:
            entry = (None, cls.missing_orders(error))
        cls._laurent_symbolic[key] = entry
        return entry, False

    @staticmethod
    def laurent_row(power, coefficient):
        """One Laurent coefficient as plain data: per component its decimal
        text, the significant digits Vakint vouches for, exact zero, error
        radius and error quality."""
        value = coefficient.numerical_value
        precision = coefficient.precision
        row = {"power": int(power)}
        for part, number, component in (
            ("real", value.real, precision.real),
            ("imaginary", value.imag, precision.imaginary),
        ):
            row[part] = {
                "text": str(number),
                "digits": component.significant_digits,
                "exact_zero": bool(component.exact_zero),
                "error": str(component.absolute_error),
                "quality": component.quality.name,
            }
        return row

    @staticmethod
    def laurent_warnings(warnings):
        """Vakint's precision warnings, one entry per message and sources with
        the ε orders it concerns."""
        grouped = {}
        for warning in warnings:
            sources = tuple(
                (source.identity, source.available_decimal_precision)
                for source in warning.sources
            )
            entry = grouped.setdefault(
                (warning.message, sources),
                {"message": warning.message, "sources": list(sources), "orders": []},
            )
            if warning.epsilon_power is not None:
                entry["orders"].append(int(warning.epsilon_power))
        return list(grouped.values())

    @classmethod
    def laurent_series(cls, case, values, *, digits, orders, progress=None):
        """The numerical Laurent series of ``case``'s whole input with every
        parameter set exactly (``values``: {name: text}), ``digits`` requested
        decimal digits and the ε orders ``orders`` = (lowest, highest).

        Vakint's own path: the symbolic-master evaluation (cached per input and
        order count, see ``laurent_symbolic``), then ``numerical_evaluation``
        at the requested precision with exact parameter values. Orders that
        fail closed keep Vakint's MissingMasterOrders reason and the series
        stops before them. Returns plain data for ``laurent_view``."""
        spec = cls.ladder_cases().get(case)
        if spec is None:
            raise ValueError(f"Unknown ladder case: {case}")
        loops = spec["loop_count"]
        digits = cls.laurent_digits(digits)
        try:
            lowest, highest = (int(order) for order in orders)
        except (TypeError, ValueError):
            raise ValueError("Choose the ε orders (lowest, highest)") from None
        lowest = max(lowest, -loops)
        if highest < lowest:
            raise ValueError(
                f"The series starts at {cls.epsilon_text(-loops)}; choose a higher order"
            )
        parameters = cls.laurent_parameters(case)
        exact = cls.laurent_values(parameters, values)
        # Wall time of the whole evaluation; "load" collects the engine
        # constructions and data loads inside it.
        timings = {"load": 0.0}
        began = time.perf_counter()
        (expression, missing), cached = cls.laurent_symbolic(
            case,
            highest + loops + 1,
            precision=digits,
            progress=progress,
            timings=timings,
        )
        through = highest
        if missing is not None:
            first = missing["first_affected"]
            through = first - 1 if isinstance(first, int) else -loops - 1
            if through >= -loops:
                (expression, again), reused = cls.laurent_symbolic(
                    case,
                    through + loops + 1,
                    precision=digits,
                    progress=progress,
                    timings=timings,
                )
                cached = cached and reused
                if again is not None:
                    expression, missing = None, again
        rows, warnings, symbols = [], [], []
        numerical_seconds = 0.0
        version = None
        if expression is not None:
            slot = cls.laurent_engine(
                loops, precision=digits, progress=progress, timings=timings
            )
            version = slot["version"]
            started = time.perf_counter()
            result, _ = slot["engine"].numerical_evaluation(
                expression,
                {
                    parameter["key"]: cls.exact_expression(exact[parameter["name"]])
                    for parameter in parameters
                },
            )
            numerical_seconds = time.perf_counter() - started
            rows = [
                cls.laurent_row(power, coefficient)
                for power, coefficient in sorted(result.coefficients.items())
            ]
            warnings = cls.laurent_warnings(result.warnings)
            names = {parameter["name"] for parameter in parameters} | {"vakint::ε"}
            symbols = sorted(
                {
                    symbol.get_name()
                    for symbol in expression.get_all_symbols(False)
                    if symbol.get_name() not in names
                }
            )
        total_seconds = time.perf_counter() - began
        return {
            "case": case,
            "label": spec["label"],
            "loop_count": loops,
            "digits": digits,
            "orders": [lowest, highest],
            "through": through,
            "coefficients": rows,
            "missing": missing,
            "warnings": warnings,
            "parameters": [
                {**parameter, "value": str(values[parameter["name"]]).strip()}
                for parameter in parameters
            ],
            "constants": symbols,
            "data": version
            or (missing or {}).get("backend")
            or "the loaded master data",
            "cached": cached,
            "seconds": {
                "total": total_seconds,
                "load": timings["load"],
                "symbolic": max(
                    0.0, total_seconds - timings["load"] - numerical_seconds
                ),
                "numerical": numerical_seconds,
            },
        }

    @staticmethod
    def rounded_decimal(text, digits):
        """``text`` (a decimal number) rounded half-even to ``digits``
        significant digits: (sign, mantissa "d.ddd…", decimal exponent)."""
        digits = max(1, int(digits))
        value = decimal.Context(prec=digits).create_decimal(str(text).strip())
        if not value.is_finite():
            raise ValueError(f"Not a finite number: {text}")
        if value.is_zero():
            return "", "0", 0
        sign, figures, exponent = value.as_tuple()
        scientific = exponent + len(figures) - 1
        figures = "".join(map(str, figures)).ljust(digits, "0")
        mantissa = figures[0] + ("." + figures[1:] if digits > 1 else "")
        return ("-" if sign else ""), mantissa, scientific

    @classmethod
    def shown_digits(cls, component, requested):
        """Digits shown for a component: the requested ones, or fewer when
        Vakint vouches for fewer."""
        achieved = component.get("digits")
        return (
            int(requested)
            if achieved is None
            else max(1, min(int(requested), achieved))
        )

    @classmethod
    def component_text(cls, component, requested):
        """Plain text of a component: "0" for an exact zero, "0 ± r" for an
        uncertain one, else the vouched-for digits in e-notation."""
        if component["exact_zero"]:
            return "0"
        if component.get("digits") is None:
            sign, mantissa, exponent = cls.rounded_decimal(component["error"], 3)
            return f"0 ± {mantissa}e{exponent}"
        sign, mantissa, exponent = cls.rounded_decimal(
            component["text"], cls.shown_digits(component, requested)
        )
        return f"{sign}{mantissa}" + (f"e{exponent}" if exponent else "")

    @staticmethod
    def _grouped(mantissa, separator=r"\mkern3mu ", more=""):
        """``mantissa`` ("d.ddd…") with its decimals in groups of five. The
        LaTeX separator is a primitive kern, not a macro such as a thin space: KaTeX
        counts macro expansions (at most 1000 per formula)."""
        head, _, tail = mantissa.partition(".")
        if not tail:
            return head + more
        groups = [tail[i : i + 5] for i in range(0, len(tail), 5)]
        return head + "." + separator.join(groups) + more

    @classmethod
    def component_latex(cls, component, requested, limit=None):
        """(sign, LaTeX body) of a non-zero component, None for an exact zero.
        At most ``limit`` digits (the headline): the leading digits of the
        full value, then an ellipsis when more are requested and vouched for."""
        if component["exact_zero"]:
            return None
        if component.get("digits") is None:
            _, mantissa, exponent = cls.rounded_decimal(component["error"], 3)
            return "", r"0\pm " + mantissa + rf"\times10^{{{exponent}}}"
        wanted = cls.shown_digits(component, requested)
        sign, mantissa, exponent = cls.rounded_decimal(component["text"], wanted)
        more = ""
        if limit is not None and wanted > max(2, int(limit)):
            mantissa, more = mantissa[: max(2, int(limit)) + 1], r"\ldots"
        power = rf"\times10^{{{exponent}}}" if exponent else ""
        return sign, cls._grouped(mantissa, more=more) + power

    @staticmethod
    def component_limit(row, limit):
        """Headline digits per component of ``row``: ``limit``, halved when
        both parts are non-zero (the line then carries two numbers)."""
        if limit is None:
            return None
        both = not (row["real"]["exact_zero"] or row["imaginary"]["exact_zero"])
        return max(2, int(limit) // 2) if both else int(limit)

    @classmethod
    def coefficient_latex(cls, row, requested, limit=None):
        """A complex coefficient in LaTeX with an aligned sign column; at most
        ``limit`` digits per number (see ``component_limit``)."""
        limit = cls.component_limit(row, limit)
        real = cls.component_latex(row["real"], requested, limit)
        imaginary = cls.component_latex(row["imaginary"], requested, limit)
        if real is None and imaginary is None:
            return r"\phantom{-}0"

        def lead(sign):
            return "-" if sign else r"\phantom{-}"

        if imaginary is None:
            return lead(real[0]) + real[1]
        if real is None:
            return lead(imaginary[0]) + imaginary[1] + r"\,i"
        joint = " - " if imaginary[0] else " + "
        return lead(real[0]) + real[1] + joint + imaginary[1] + r"\,i"

    @classmethod
    def headline_truncated(cls, summary):
        """Whether the headline series shows fewer digits than the table of
        every requested digit (some coefficient longer than its headline)."""
        requested = summary["digits"]
        for row in cls.shown_rows(summary):
            limit = cls.component_limit(row, cls.laurent_headline_digits)
            for component in (row["real"], row["imaginary"]):
                if (
                    not component["exact_zero"]
                    and component.get("digits") is not None
                    and cls.shown_digits(component, requested) > limit
                ):
                    return True
        return False

    @classmethod
    def shown_rows(cls, summary):
        """The coefficients of ``summary`` inside its ε range, below the
        first order that fails closed."""
        lowest, highest = summary["orders"]
        top = min(highest, summary["through"])
        return [row for row in summary["coefficients"] if lowest <= row["power"] <= top]

    @classmethod
    def component_digits(cls, component, requested):
        """Plain text of a component with every shown digit, decimals in
        groups of five and the power of ten in superscript (wraps between
        groups): "−6.01945 29704 … × 10⁻⁷"; "0" for an exact zero."""
        if component["exact_zero"]:
            return "0"
        if component.get("digits") is None:
            _, mantissa, exponent = cls.rounded_decimal(component["error"], 3)
            return f"0 ± {mantissa} × 10{str(exponent).translate(cls._superscripts)}"
        sign, mantissa, exponent = cls.rounded_decimal(
            component["text"], cls.shown_digits(component, requested)
        )
        power = f" × 10{str(exponent).translate(cls._superscripts)}" if exponent else ""
        return ("−" if sign else "") + cls._grouped(mantissa, separator=" ") + power

    @classmethod
    def coefficient_digits(cls, row, requested):
        """A complex coefficient for the list of every requested digit."""
        real, imaginary = row["real"], row["imaginary"]
        if imaginary["exact_zero"]:
            return cls.component_digits(real, requested)
        part = cls.component_digits(imaginary, requested)
        if real["exact_zero"]:
            return part + " i"
        joint = " − " if part.startswith("−") else " + "
        return (
            cls.component_digits(real, requested)
            + joint
            + "("
            + part.removeprefix("−")
            + ") i"
        )

    @classmethod
    def laurent_plain(cls, summary):
        """The shown coefficients as copyable plain text, one per line, every
        shown digit in e-notation: "c[-3] = 0 + 6.0194…e-7*I"."""
        requested = summary["digits"]
        return "\n".join(
            f"c[{row['power']}] = {cls.component_text(row['real'], requested)}"
            f" + ({cls.component_text(row['imaginary'], requested)})*I"
            for row in cls.shown_rows(summary)
        )

    _superscripts = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

    @classmethod
    def epsilon_text(cls, power):
        """ε with a superscript power for plain text: ε⁻³, ε⁰, ε¹."""
        return "ε" + str(int(power)).translate(cls._superscripts)

    @staticmethod
    def epsilon_term(power):
        if power == 0:
            return "c_{0}"
        if power == 1:
            return r"c_{1}\,\varepsilon"
        return rf"c_{{{power}}}\,\varepsilon^{{{power}}}"

    @classmethod
    def laurent_latex(cls, summary):
        """The series of a ``laurent_series`` summary as KaTeX markdown: the
        sum with its 𝒪(ε) tail, then every coefficient (aligned digits, at
        most ``laurent_headline_digits`` each, so that a line stays readable
        and KaTeX stays far below its expansion limit); an order that fails
        closed shows that instead of a number."""
        lowest, highest = summary["orders"]
        through = summary["through"]
        requested = summary["digits"]
        limit = cls.laurent_headline_digits
        shown = cls.shown_rows(summary)
        failed = list(range(max(lowest, through + 1), highest + 1))
        terms = [cls.epsilon_term(row["power"]) for row in shown]
        if lowest > -summary["loop_count"]:
            terms.insert(0, r"\cdots")
        last = shown[-1]["power"] if shown else lowest - 1
        tail = last + 1
        terms.append(
            r"\mathcal{O}(\varepsilon)"
            if tail == 1
            else rf"\mathcal{{O}}(\varepsilon^{{{tail}}})"
        )
        lines = [
            rf"c_{{{row['power']}}} &= {cls.coefficient_latex(row, requested, limit)}"
            for row in shown
        ] + [
            rf"c_{{{power}}} &\quad \text{{fails closed: MissingMasterOrders}}"
            for power in failed
        ]
        text = "$$\n\\mathcal{I}(\\varepsilon) = " + " + ".join(terms) + "\n$$\n"
        if lines:
            text += (
                "\n$$\n\\begin{aligned}\n"
                + " \\\\\n".join(lines)
                + "\n\\end{aligned}\n$$\n"
            )
        return text

    @classmethod
    def laurent_digits_view(cls, mo, summary):
        """Every requested digit of the shown coefficients (as many as Vakint
        vouches for), one coefficient per row, wrapping between groups of five
        decimals; the power of ten and the i stay at the end of the value."""
        requested = summary["digits"]
        rows = "".join(
            '<div class="vg-digit-row"><span class="vg-digit-k">c<sub>'
            + html.escape(str(row["power"]).replace("-", "−"))
            + "</sub> =</span><code>"
            + html.escape(cls.coefficient_digits(row, requested))
            + "</code></div>"
            for row in cls.shown_rows(summary)
        )
        return mo.Html(
            '<div class="vg-digits"><div class="vg-pipe-label">Every requested '
            f"digit · {requested} digits, as many as Vakint vouches for</div>"
            + rows
            + "</div>"
        )

    @classmethod
    def laurent_view(cls, mo, summary, diagram):
        """The end of the pipeline at a glance: the input graph, an arrow,
        the Laurent series with numerical coefficients (headline digits);
        every requested digit below it when the headline is shorter; then the
        coefficient table with Vakint's precision, the fail-closed orders and
        warnings."""
        requested = summary["digits"]
        headline = cls.laurent_headline_digits
        truncated = cls.headline_truncated(summary)
        values = ", ".join(
            f"{parameter['label']} = {parameter['value']}"
            for parameter in summary["parameters"]
        )
        series = mo.md(cls.laurent_latex(summary))
        frame = cls.diagram_frame(mo, diagram, height=340)
        caption = (
            "Feynman rules · Dirac and color algebra · Ordered IBP reduction · "
            f"masters from the {summary['data']} · {requested}-digit evaluation"
        )
        pipeline = mo.Html(
            '<div class="vg-pipe-wrap"><div class="vg-pipeline">'
            '<div class="vg-pipe-graph">'
            '<div class="vg-pipe-label">Input graph</div>'
            + frame.text
            + f'<div class="vg-pipe-caption">{html.escape(summary["label"])}</div>'
            + '</div><div class="vg-pipe-arrow"><b>→</b>'
            + f"<small>{html.escape(caption)}</small></div>"
            + '<div class="vg-pipe-series"><div class="vg-pipe-label">'
            + "Laurent series of the whole input</div>"
            + series.text
            + '<div class="vg-pipe-caption">MSbar, d = 4 − 2ε · at '
            + html.escape(values)
            + (
                f" · leading {headline} digits of each coefficient; all "
                f"{requested} requested digits below"
                if truncated
                else ""
            )
            + "</div></div></div></div>"
        )
        shown = cls.shown_rows(summary)
        loops = summary["loop_count"]
        through = summary["through"]
        achieved = [
            component["digits"]
            for row in shown
            for component in (row["real"], row["imaginary"])
            if component["digits"] is not None
        ]

        def code(text):
            return mo.Html(f"<code>{html.escape(text)}</code>")

        def vouched(component):
            if component["exact_zero"]:
                return "exact 0"
            return "—" if component["digits"] is None else str(component["digits"])

        # The table rounds long values to the headline digits; the list
        # above it then carries every requested digit.
        table_digits = min(requested, headline)
        table_rows = [
            {
                "ε power": row["power"],
                "real part": code(cls.component_text(row["real"], table_digits)),
                "imaginary part": code(
                    cls.component_text(row["imaginary"], table_digits)
                ),
                "digits Vakint vouches for (re / im)": f"{vouched(row['real'])} / "
                f"{vouched(row['imaginary'])}",
                "error radius (re / im)": " / ".join(
                    "0"
                    if component["exact_zero"]
                    else cls.component_text(
                        {**component, "text": component["error"], "digits": 3},
                        3,
                    )
                    for component in (row["real"], row["imaginary"])
                ),
                "error quality": " / ".join(
                    sorted({row["real"]["quality"], row["imaginary"]["quality"]})
                ),
            }
            for row in shown
        ]
        seconds = summary["seconds"]
        total = seconds.get(
            "total", seconds.get("symbolic", 0.0) + seconds.get("numerical", 0.0)
        )
        notes = []
        if seconds.get("load", 0.0) >= 0.05:
            notes.append(f"engine and data load {cls.duration(seconds['load'])}")
        if summary["cached"]:
            notes.append("symbolic result reused")
        parts = [
            pipeline,
            *([cls.laurent_digits_view(mo, summary)] if truncated else []),
            cls.metrics(
                mo,
                [
                    (requested, "requested decimal digits"),
                    (
                        min(achieved) if achieved else "—",
                        "fewest digits Vakint vouches for (any shown coefficient)",
                    ),
                    (
                        f"{cls.epsilon_text(-loops)} … {cls.epsilon_text(through)}"
                        if through >= -loops
                        else "none",
                        "orders Vakint computes for this input",
                    ),
                    (
                        cls.duration(total),
                        "this evaluation, all of it"
                        + (f" ({'; '.join(notes)})" if notes else ""),
                    ),
                ],
            ),
            cls.table(
                mo,
                table_rows,
                numeric=("ε power",),
                caption="Each value shows the digits Vakint vouches for, at most the "
                + (
                    f"requested ones (rounded to {headline} here; every requested "
                    "digit is listed above). "
                    if requested > headline
                    else "requested ones. "
                )
                + "Bounded: a guaranteed error radius from tracked arithmetic; "
                "Estimated: a propagated estimate.",
            )
            if table_rows
            else cls.note(mo, "No order in the chosen range", "Widen the ε range."),
        ]
        missing = summary.get("missing")
        if missing is not None:
            first = through + 1
            parts.append(
                cls.note(
                    mo,
                    f"{cls.epsilon_text(first)} and higher fail closed",
                    f"Vakint computes this input through {cls.epsilon_text(through)}. "
                    f"For {cls.epsilon_text(first)} "
                    f"{cls.missing_orders_text(missing)}. Vakint raises "
                    "MissingMasterOrders instead of returning a number that would "
                    "silently lack those terms.",
                    pending=True,
                )
            )
        warnings = summary.get("warnings") or []
        if warnings:
            parts.append(
                cls.note(
                    mo,
                    f"Precision warnings from Vakint ({len(warnings)})",
                    "\n".join(
                        "• "
                        + warning["message"]
                        + (
                            " ("
                            + ", ".join(
                                map(cls.epsilon_text, sorted(set(warning["orders"])))
                            )
                            + ")"
                            if warning["orders"]
                            else ""
                        )
                        + (
                            "; sources: "
                            + ", ".join(
                                f"{identity} ({available} digits)"
                                if available is not None
                                else str(identity)
                                for identity, available in warning["sources"]
                            )
                            if warning["sources"]
                            else ""
                        )
                        for warning in warnings
                    ),
                    pending=True,
                )
            )
        else:
            parts.append(
                cls.legend(
                    mo,
                    "Vakint reports no precision warning: exact parameter values "
                    "and the data's master coefficients support every shown digit.",
                )
            )
        details = {
            "Parameters and constants of the evaluated series": cls.table(
                mo,
                [
                    {
                        "symbol": parameter["name"],
                        "role": parameter["meaning"],
                        "value": code(parameter["value"]),
                    }
                    for parameter in summary["parameters"]
                ]
                + [
                    {
                        "symbol": name if len(name) <= 40 else name[:36] + "…",
                        "role": "constant, evaluated by Vakint from its data at the "
                        "requested precision",
                        "value": "—",
                    }
                    for name in summary["constants"]
                ],
            ),
        }
        if shown:
            details["Coefficients as plain text (every shown digit)"] = mo.plain_text(
                cls.laurent_plain(summary)
            )
        if missing is not None:
            details["Vakint's MissingMasterOrders message"] = mo.plain_text(
                missing["message"]
            )
        parts.append(mo.accordion(details))
        return mo.vstack(parts)

    # -- Script mode of the ladder page ---------------------------------------
    headless_options = ("case", "n-cores", "arithmetic", "max-seconds", "cancel-after")
    headless_usage = (
        "usage: python 05_ladder_reductions.py --case CASE [--n-cores N] "
        "[--arithmetic exact|modular_sampled|modular_polynomial] "
        "[--max-seconds SECONDS] [--cancel-after SECONDS]"
    )

    @classmethod
    def headless_request(cls, arguments, defaults):
        """The one run a script-mode ladder page performs.

        ``arguments`` are ``mo.cli_args()``, ``defaults`` the page form's values.
        Raises ValueError, with the usage, for a missing case, an unknown option
        or a malformed number; nothing is started then."""
        arguments = dict(arguments)
        unknown = sorted(set(arguments) - set(cls.headless_options))
        if unknown:
            raise ValueError(
                "Unknown option "
                + ", ".join(f"--{key}" for key in unknown)
                + "; "
                + cls.headless_usage
            )
        case = arguments.get("case")
        if case is None or case == "" or isinstance(case, (bool, list)):
            raise ValueError(
                "Name one case with --case ("
                + ", ".join(cls.ladder_cases())
                + "); "
                + cls.headless_usage
            )

        def number(key, default, integral=False):
            value = arguments.get(key, default)
            if value is None:
                return None
            try:
                if isinstance(value, (bool, list)) or value == "":
                    raise ValueError
                value = float(value)
                if integral and not value.is_integer():
                    raise ValueError
            except (TypeError, ValueError):
                raise ValueError(
                    f"--{key} needs one {'whole ' if integral else ''}number; "
                    + cls.headless_usage
                ) from None
            return int(value) if integral else value

        cancel_after = number("cancel-after", None)
        if cancel_after is not None and not cancel_after >= 0:
            raise ValueError("--cancel-after needs a non-negative number of seconds")
        return {
            "case": str(case),
            "n_cores": number("n-cores", defaults["n_cores"], integral=True),
            "arithmetic": str(arguments.get("arithmetic") or defaults["arithmetic"]),
            "max_seconds": number("max-seconds", defaults["max_seconds"]),
            "ordered_cache": bool(defaults["ordered_cache"]),
            "cancel_after": cancel_after,
        }

    @staticmethod
    def headless_exit(request, error):
        """A refused script-mode start: one JSON line, the reason on stderr,
        exit status 2. Nothing else (no other session's job) is touched."""
        summary = {
            "ok": False,
            "state": "not started",
            "case": (request or {}).get("case"),
            "error": error,
        }
        print(json.dumps(summary, ensure_ascii=False), flush=True)
        print(f"05_ladder_reductions.py: {error}", file=sys.stderr, flush=True)
        raise SystemExit(2)

    @classmethod
    def headless_summary(cls, work, status, request):
        """The outcome of a script-mode run and its exit status.

        0 only if the run completed, its retained sha256 equals the recorded
        fixture and the result closes on catalogue masters, or if it was
        cancelled as ``--cancel-after`` asked with no surviving process. 130
        after Ctrl-C (the job is cancelled), 1 otherwise."""
        state = status.get("state")
        result = cls.job_result(work) or {}
        retained = result.get("retained") or {}
        closure = result.get("catalogue_closure") or {}
        seconds = result.get("seconds") or {}
        report = status.get("cancel_report") or cls.cancel_report(work)
        summary = {
            "job": str(work),
            "case": request.get("case"),
            "state": state,
            "execution": result.get("execution"),
            "data": result.get("data"),
            "retained_sha256": retained.get("sha256"),
            "recorded_sha256": retained.get("recorded_sha256"),
            "matches_recorded": retained.get("matches_recorded"),
            "catalogue_closure": {
                "master_count": closure.get("master_count"),
                "closed": closure.get("closed"),
            }
            if closure
            else None,
            "seconds": {
                key: seconds[key]
                for key in ("Ordered reduction", "total")
                if key in seconds
            },
        }
        if report:
            summary["cancel_report"] = {
                "signalled": len(report.get("signalled") or []),
                "survivors": report.get("survivors") or [],
            }
        if state != "complete" or "error" in result:
            summary["message"] = str(
                result.get("error") or status.get("message") or ""
            )[-1000:]
        reproduced = (
            state == "complete"
            and result.get("schema") == "vakint.gallery.ladder-run.v1"
            and "error" not in result
            and retained.get("matches_recorded") is True
            and closure.get("closed") is True
        )
        cancelled_as_asked = (
            state == "cancelled"
            and request.get("cancel_after") is not None
            and not status.get("interrupted")
            and bool(report)
            and not report.get("survivors")
        )
        summary["ok"] = reproduced or cancelled_as_asked
        if status.get("interrupted"):
            summary["interrupted"] = True
            return summary, 130
        return summary, 0 if summary["ok"] else 1

    @classmethod
    def headless_report(cls, work, status, request):
        """Print a script-mode run's one-line JSON summary; exit non-zero
        (SystemExit) unless it reproduced its fixture or was cancelled as asked."""
        summary, code = cls.headless_summary(work, status, request)
        print(json.dumps(summary, ensure_ascii=False), flush=True)
        if code:
            if summary.get("interrupted"):
                reason = "interrupted; the job was cancelled"
            elif summary["state"] == "complete":
                reason = (
                    "the retained result differs from the recorded fixture"
                    if summary["matches_recorded"] is not True
                    else "the run did not produce a closed, recorded result"
                )
            else:
                reason = f"the run ended {summary['state']}"
            message = summary.get("message")
            print(
                f"05_ladder_reductions.py: {reason}"
                + (f": {message}" if message else "")
                + f" ({work})",
                file=sys.stderr,
                flush=True,
            )
            raise SystemExit(code)
        return summary

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

    @staticmethod
    def checked_bytes(path, digest):
        """A file's bytes, only if they match the recorded sha256."""
        payload = Path(path).read_bytes()
        if hashlib.sha256(payload).hexdigest() != digest:
            raise ValueError(f"Evidence checksum mismatch: {Path(path).name}")
        return payload

    @classmethod
    def five_loop_fixture(cls):
        """The recorded current-basis five-loop zero-momentum run (no reduction):
        its report, the native monitor's final model and the retained result's
        path, each checked against data/manifest.json; None if not shipped."""
        report = cls.receipt("five-loop-q0")
        if report is None:
            return None
        fixture = {"report": report, "model": None, "retained": None}
        for index in range(8):
            resolved = cls.evidence_file("five-loop-q0", index)
            if resolved is None:
                break
            path, digest = resolved
            if path.name == "monitor-model.json":
                fixture["model"] = json.loads(cls.checked_bytes(path, digest))
            elif path.suffix == ".symbolica":
                fixture["retained"] = (path, digest)
        return fixture

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
