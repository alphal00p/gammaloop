#!/usr/bin/env python
"""Rebuild every data file of the gamma5 documentation site (DATA_FORMAT.md) in one command.

    SYMBOLICA_LICENSE=... python tools/build_site.py [options]

Run it with a Python that provides Symbolica with the HepKIT community modules. The licence is read from the
environment only; it is never written anywhere. The package ``python/hepkit_renorm`` of this repository and the
package ``gdocs`` next to this file are put on ``sys.path`` by the script itself.

Steps, in order (each can be skipped with ``--skip NAME``):

``viewer``     Symbolica pager module and math font in ``site/assets/viewer`` (``mathview.write_assets``).
``graphs``     ``data/classes.json`` and ``data/graphs/*.json`` from the automation's enumeration, drawn with HepKIT
               (``gdocs.graphs.build_graphs``). ``--reuse-graphs`` keeps the existing files instead; ``--graph-cache``
               reuses rendered drawings; ``--contributions`` first re-evaluates every graph (about 7 min, 32 cores).
``unbroken``   ``data/unbroken.json`` and the class labels, results and summaries (``gdocs.catalog``).
``groups``     group titles, report references and summaries in ``data/classes.json`` (``gdocs.theory``).
``theory``     ``data/theory.json`` (``gdocs.theory.build_theory``); the sources are checked first.
``broken``     ``data/broken.json`` (``gdocs.broken``), when that module exists.
``resources``  ``site/resources`` (copies refreshed from ``--sources`` when given, scrubbed of local paths, user and
               host names), ``site/resources/README.md`` and ``data/resources.json``.
``manifest``   ``data/manifest.json`` (texts from ``gdocs.theory.manifest_texts``; date from ``--date``, else the
               environment variable ``SOURCE_DATE``, else today).
``validate``   ``node tools/tests/check_site.mjs SITE --no-mock``.

``--galleries`` also rewrites the two check pages (``data/graphs_gallery.html``, ``assets/viewer/demo.html``).
Files that still carry the mock marker (``"mock": true``, written by ``tools/mock_data.py``) are removed after the
build; a mock ``broken.json`` without ``gdocs.broken`` is removed and the broken phase is left out of the manifest.
Every step that fails stops the build with exit status 1. A full build takes about one minute on 32 cores;
``--reuse-graphs`` saves about 15 s.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import getpass
import importlib
import importlib.util
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent            # gamma5_docs/tools
DOCS = HERE.parent                                # gamma5_docs
REPO = DOCS.parent                                # repository root
PKG = REPO / "python" / "hepkit_renorm"           # the automation package (read only)
DEFAULT_SITE = DOCS / "site"
STEPS = ["viewer", "graphs", "unbroken", "groups", "theory", "broken", "resources", "manifest", "validate"]

for p in (str(HERE), str(PKG)):
    if p not in sys.path:
        sys.path.insert(0, p)


# ================================================================================================ helpers
class BuildError(RuntimeError):
    pass


def log(msg: str) -> None:
    print(f"  {msg}", flush=True)


def _is_mock(path: Path) -> bool:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if isinstance(data, dict):
        return data.get("mock") is True
    if isinstance(data, list):
        return bool(data) and all(isinstance(x, dict) and x.get("mock") is True for x in data)
    return False


def _write_json(path: Path, data, *, indent=None) -> None:
    text = json.dumps(data, ensure_ascii=False, indent=indent, separators=None if indent else (",", ":"))
    path.write_text(text + "\n", encoding="utf-8")


def _size(n: int) -> str:
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


# ================================================================================================ steps
def step_viewer(site: Path, args) -> list:
    from gdocs import mathview as mv

    out = mv.write_assets(site)
    rows = [(p, "pager module" if p.suffix == ".js" else "math font") for p in out]
    rows.append((site / "assets" / "viewer" / "mathview.js", "static adapter (hand-written, kept)"))
    if args.galleries:
        from gdocs import demo_gallery

        demo_gallery.build(site)
        rows.append((site / "assets" / "viewer" / "demo.html", "viewer check page"))
    return rows


def step_graphs(site: Path, args) -> list:
    from gdocs import graphs

    cls = site / "data" / "classes.json"
    if args.reuse_graphs:
        if not cls.exists() or _is_mock(cls):
            raise BuildError("--reuse-graphs: data/classes.json is missing or mock; build the graphs once")
        doc = json.loads(cls.read_text(encoding="utf-8"))
        missing = [c["graphs"] for c in doc["classes"] if not (site / c["graphs"]).exists()]
        if missing:
            raise BuildError(f"--reuse-graphs: {len(missing)} graph files missing, e.g. {missing[:3]}")
        log(f"reusing classes.json and {len(doc['classes'])} graph files")
        s = {"classes": len(doc["classes"]), "graphs": sum(c["n_graphs"] for c in doc["classes"])}
    else:
        if args.contributions:
            graphs.compute_contributions(site, processes=args.processes or 16)
        s = graphs.build_graphs(site, processes=args.processes, cache=args.graph_cache, log=log)
    gdir = site / "data" / "graphs"
    rows = [(cls, f"{s['classes']} classes"),
            (gdir, f"{s['graphs']} graphs in {len(list(gdir.glob('*.json')))} files")]
    if args.galleries:
        rows.append((graphs.write_gallery(site), "graph check page"))
    return rows


def step_unbroken(site: Path, args) -> list:
    from gdocs import catalog

    counts = catalog.build_unbroken(site, classes=True, log=log)
    by_type = ", ".join(f"{k} {v}" for k, v in counts["by_type"].items())
    return [(site / "data" / "unbroken.json", f"{counts['terms']} terms ({by_type})"),
            (site / "data" / "classes.json", "labels, results, summaries")]


def step_groups(site: Path, args) -> list:
    from gdocs import theory

    doc = theory.merge_class_groups(site)
    n = sum(1 for g in doc["groups"] if g.get("summary"))
    return [(site / "data" / "classes.json", f"{n} of {len(doc['groups'])} groups with a summary")]


def step_theory(site: Path, args) -> list:
    from gdocs import theory

    problems = theory.check(site)
    if problems:
        raise BuildError("theory sources: " + "; ".join(problems))
    data = theory.build_theory(site)
    return [(site / "data" / "theory.json", f"{len(data['cards'])} cards")]


def step_broken(site: Path, args) -> list:
    path = site / "data" / "broken.json"
    if importlib.util.find_spec("gdocs.broken") is None:
        if path.exists() and _is_mock(path):
            path.unlink()
            log("gdocs.broken does not exist: removed the mock broken.json; the broken phase is left out")
            return [(path, "removed (mock; gdocs.broken missing)")]
        if path.exists():
            log("gdocs.broken does not exist: kept the existing broken.json")
            return [(path, "kept (gdocs.broken missing)")]
        log("gdocs.broken does not exist: no broken phase")
        return [(path, "not built (gdocs.broken missing)")]
    mod = importlib.import_module("gdocs.broken")
    if hasattr(mod, "build_broken"):
        mod.build_broken(site, log=log)
    elif hasattr(mod, "main"):
        rc = mod.main(["--site", str(site)])
        if rc:
            raise BuildError(f"gdocs.broken.main returned {rc}")
    else:
        raise BuildError("gdocs.broken has neither build_broken(site_dir, log=...) nor main(argv)")
    if not path.exists() or _is_mock(path):
        raise BuildError("gdocs.broken did not write a real data/broken.json")
    doc = json.loads(path.read_text(encoding="utf-8"))
    return [(path, f"{len(doc.get('terms', []))} terms")]


# ------------------------------------------------------------------------------------------------ resources
# Copies refreshed from the project directory (``--sources``: the directory that holds ``report/``,
# ``SCRATCHPAD/`` and this repository). (source relative to it, destination relative to site/resources).
def _source_plan(src: Path) -> list[tuple[Path, str]]:
    plan = []
    rep = src / "report"
    for name in ("one_loop_gamma5.typ", "macros.typ", "diagrams.typ", "refs.bib"):
        plan.append((rep / name, f"report_src/{name}"))
    for p in sorted((rep / "sections").glob("*.typ")):
        plan.append((p, f"report_src/sections/{p.name}"))
    for p in sorted((rep / "sections" / "classes").glob("*.typ")):
        if not p.name.startswith("_"):
            plan.append((p, f"report_src/sections/classes/{p.name}"))
    sp = src / "SCRATCHPAD"
    plan.append((sp / "scripts" / "build_report.py", "report_src/build_report.py"))
    for p in sorted((sp / "classes").glob("*/final_result.json")):
        plan.append((p, f"report_classes/{p.parent.name}/final_result.json"))
    for name in ("ct_assembly.json", "pole_consistency.json", "x_solve_engine.json", "x_solve_hand.json"):
        plan.append((sp / "phase4" / name, f"report_checks/{name}"))
    plan.append((sp / "phase4" / "final_counterterms.json", "final_counterterms.json"))
    for name in ("sm_unbroken_bmhv.json", "sm_unbroken_bmhv_native.json", "ENUMERATION.md", "README.md"):
        plan.append((sp / "model" / name, f"report_model/{name}"))
    plan.append((sp / "toolkit" / "validation" / "VALIDATION.md", "toolkit/VALIDATION.md"))
    plan.append((sp / "toolkit" / "README.md", "toolkit/README.md"))
    return plan


# Always copied from the repository.
REPO_PLAN = [(PKG / "hepkit_renorm" / "models" / "sm_one_family.toml", "automation/sm_one_family.toml")]

TEXT_SUFFIXES = {".json", ".txt", ".md", ".typ", ".bib", ".py", ".toml"}
MAX_FILE = 2_000_000
MAX_TOTAL = 12_000_000
_KEY = re.compile(r"S-\d{3}-\d{4}\.\d{2}\.\d{2}-[A-Z0-9]{16,}")
_ABS = re.compile(r"(?<![\w.<>-])/(?:shared|home|common|nix|tmp|Users|var/folders)/[^\s\"'<>)\]`,;]*")


def _rel_path(m: re.Match) -> str:
    """A local absolute path, rewritten relative to the directory it belongs to."""
    p = m.group(0)
    for anchor in ("/run_sm_one_family",):
        if anchor in p:
            rest = p.split(anchor, 1)[1].lstrip("/")
            return rest or "."
    if "/one_loop_explicit_calculation/" in p:
        return p.split("/one_loop_explicit_calculation/", 1)[1]
    if "/gamma5/" in p:
        return "<scratch>/" + p.split("/gamma5/", 1)[1]
    return "<local path>"


def scrub_text(text: str, names: list[str]) -> str:
    """Remove absolute local paths, the user name and the host name (``names``) from ``text``."""
    # the report's build script locates the report next to itself once copied
    text = re.sub(r'^REPORT = "[^"\n]*"$', 'REPORT = os.path.dirname(os.path.abspath(__file__))', text, flags=re.M)
    text = _ABS.sub(_rel_path, text)
    text = re.sub(r'("host"\s*:\s*)"[^"]*"', r'\1"(removed)"', text)
    for n in names:
        if n and len(n) >= 3:
            text = re.sub(rf"\b{re.escape(n)}\b", "<removed>", text)
    return text


# One title and one sentence per file; group order is display order.
RES_GROUPS = [
    ("report", "Report"),
    ("automation", "Automation: model specification and run outputs"),
    ("results", "Final numbers of the report"),
    ("classes", "Per-class results of the report"),
    ("checks", "Assembly and checks of the report"),
    ("model", "Hand-written HepKIT model of the report"),
    ("toolkit", "Stand-alone toolkit"),
    ("report_src", "Typst sources of the report"),
]
RES_INFO = {
    "report/one_loop_gamma5.pdf": ("report", "Report (PDF)",
                                   "The report “One-loop renormalisation of the Standard Model in the unbroken phase with a non-anticommuting gamma5” (PDF, 186 pages)."),
    "resources/report_headings.txt": ("report", "Report headings and pages",
                                      "First page of every labelled heading of the report, used for the page references of this site."),
    "resources/automation/sm_one_family.toml": ("automation", "Model specification (TOML)",
                                                "Input of the automation: gauge group, chiral multiplets, scalar doublet, classical four-dimensional Lagrangian, gauge fixing and regularisation of the one-family Standard Model."),
    "resources/run_sm_one_family/sm_one_family_tree.json": ("automation", "Tree-level model (HepKIT JSON)",
                                                            "Regularised classical action with ghosts, sources and the breaking insertion, in the Feynman gauge, as a HepKIT model."),
    "resources/run_sm_one_family/sm_one_family_renormalised.json": ("automation", "Renormalised model (HepKIT JSON)",
                                                                    "The renormalised model: classical action, breaking, singular and finite restoring counterterms and their breaking insertions, as a HepKIT model."),
    "resources/run_sm_one_family/sm_one_family_store.json": ("automation", "Order store",
                                                             "Exact order store of the run (all functionals of order one), the input of the next loop order."),
    "resources/run_sm_one_family/sm_one_family_order_data.json": ("automation", "Order data",
                                                                  "Matched loop data of order one as exact local functionals: the pole coefficients and the finite breaking."),
    "resources/run_sm_one_family/run_report.json": ("automation", "Run report (JSON)",
                                                    "Every step, check and result of the run: graph counts per class, renormalisation constants, beta functions, breaking and restoring coefficients."),
    "resources/run_sm_one_family/run_report.txt": ("automation", "Run report (text)",
                                                   "Readable summary of the run report: steps with timings, the fifteen checks and the results of order one."),
    "resources/run_sm_one_family/comparison.json": ("automation", "Comparison with the report",
                                                    "Comparison of the run with the report's final numbers, quantity by quantity."),
    "resources/final_counterterms.json": ("results", "Final counterterms",
                                          "The report's final numbers: singular counterterm, finite breaking, restoring coefficients and literature comparison."),
    "resources/report_checks/ct_assembly.json": ("checks", "Assembly of the singular counterterm",
                                                 "Over-determined assembly of the singular counterterm from the eighteen classes: monomial coefficients, renormalisation constants, beta functions and evanescent part."),
    "resources/report_checks/pole_consistency.json": ("checks", "Pole identity",
                                                      "Check of the pole identity for the Green's functions with one breaking insertion, leg assignment by leg assignment."),
    "resources/report_checks/x_solve_engine.json": ("checks", "Restoring counterterm (engine solution)",
                                                    "Solution of the restoring condition with the plane-wave engine: the finite breaking by class and the restoring coefficients."),
    "resources/report_checks/x_solve_hand.json": ("checks", "Restoring counterterm (hand solution)",
                                                  "Independent hand solution of the restoring condition under the normalisation conditions."),
    "resources/report_model/sm_unbroken_bmhv.json": ("model", "HepKIT model of the report",
                                                     "Hand-written HepKIT JSON model of the regularised unbroken Standard Model with sterile partners, used for the report's graph enumeration."),
    "resources/report_model/sm_unbroken_bmhv_native.json": ("model", "HepKIT model with true ghost numbers",
                                                            "Variant of the hand-written model with anticommuting (native) ghosts."),
    "resources/report_model/ENUMERATION.md": ("model", "Graph enumeration summary",
                                              "Counts of the one-loop one-particle-irreducible graphs per class and their comparison with hand counts."),
    "resources/report_model/README.md": ("model", "Notes on the hand-written model",
                                         "Encoding of the hand-written model and description of the enumeration scripts (which are not copied here)."),
    "resources/toolkit/VALIDATION.md": ("toolkit", "Validation notes of the toolkit",
                                        "Independent validation of the stand-alone BMHV toolkit bmhvkit: what was checked and how."),
    "resources/toolkit/README.md": ("toolkit", "Toolkit interface",
                                    "Interface of the stand-alone BMHV and one-loop toolkit bmhvkit (the code itself is not copied here)."),
    "resources/report_src/build_report.py": ("report_src", "Report build script",
                                             "Assembles the class section and compiles the report with Typst."),
    "resources/report_src/one_loop_gamma5.typ": ("report_src", "Main Typst file",
                                                 "Title, abstract and the list of included sections."),
    "resources/report_src/macros.typ": ("report_src", "Notation macros", "Typst macros of the report's notation."),
    "resources/report_src/diagrams.typ": ("report_src", "Diagram macros", "Typst (CeTZ) macros that draw the report's Feynman graphs."),
    "resources/report_src/refs.bib": ("report_src", "Bibliography", "Bibliography of the report (BibTeX)."),
}


def _class_titles(site: Path) -> dict:
    try:
        doc = json.loads((site / "data" / "classes.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return {g["key"]: (g.get("title") or g["key"], (g.get("report") or [{}])[0]) for g in doc.get("groups", [])}


def _describe(rel: str, site: Path, ctitles: dict) -> tuple[str, str, str]:
    if rel in RES_INFO:
        return RES_INFO[rel]
    m = re.match(r"resources/report_classes/([^/]+)/final_result\.json$", rel)
    if m:
        key = m.group(1)
        title, ref = ctitles.get(key, (key, {}))
        where = f" ({ref['label']})" if ref.get("label") else ""
        return ("classes", f"{key}: final result",
                f"Final numbers of the class {key}, {title}{where}: graphs, loop functions, local operators, checks and caveats of the two independent calculations.")
    m = re.match(r"resources/report_src/sections/(?:classes/)?([^/]+)\.typ$", rel)
    if m:
        name = m.group(1)
        head = _typ_heading(site / rel)
        return ("report_src", f"{name}.typ", f"Typst source of the report: {head}." if head else "Typst source of the report.")
    raise BuildError(f"resources: no description for {rel}; add it to RES_INFO in build_site.py")


def _typ_heading(path: Path) -> str:
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            m = re.match(r"^(=+)\s+(.*?)\s*(<[^>]+>)?\s*$", line)
            if m:
                t = re.sub(r"\$[^$]*\$", "", m.group(2))
                t = re.sub(r"#\w+(\([^)]*\))?|[_*]", "", t)
                t = re.sub(r"\s+", " ", t).strip(" :,")
                return t
    except OSError:
        pass
    return ""


def step_resources(site: Path, args) -> list:
    res = site / "resources"
    res.mkdir(exist_ok=True)
    plan = list(REPO_PLAN)
    if args.sources:
        src = Path(args.sources).resolve()
        if not (src / "report").is_dir() or not (src / "SCRATCHPAD").is_dir():
            raise BuildError(f"--sources {src}: expected the directories report/ and SCRATCHPAD/ in it")
        plan += _source_plan(src)
    copied = 0
    for s, rel in plan:
        if not s.is_file():
            raise BuildError(f"resources: source file missing: {s.name} (for {rel})")
        if s.stat().st_size > MAX_FILE:
            raise BuildError(f"resources: {s.name} is {_size(s.stat().st_size)}, above the 2 MB limit")
        dst = res / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(s, dst)
        copied += 1
    # scrub every text file in place (idempotent)
    names = [getpass.getuser(), socket.gethostname().split(".")[0]]
    scrubbed = 0
    for p in sorted(res.rglob("*")):
        if not p.is_file() or p.suffix not in TEXT_SUFFIXES or p == res / "README.md":
            continue
        t = p.read_text(encoding="utf-8")
        u = scrub_text(t, names)
        if u != t:
            p.write_text(u, encoding="utf-8")
            scrubbed += 1
        if _KEY.search(u):
            raise BuildError(f"resources: {p.relative_to(site)} contains a licence key")
        if p.suffix == ".json":
            json.loads(u)  # still valid after scrubbing
    # list
    ctitles = _class_titles(site)
    entries = []
    files = [site / "report" / "one_loop_gamma5.pdf"] + sorted(
        p for p in res.rglob("*") if p.is_file() and p != res / "README.md")
    total = 0
    for p in files:
        rel = p.relative_to(site).as_posix()
        if any(part.startswith(".") or part == "__pycache__" for part in p.relative_to(site).parts):
            raise BuildError(f"resources: hidden or cache file {rel}")
        if rel.startswith("resources/"):
            total += p.stat().st_size
            if p.stat().st_size > MAX_FILE:
                raise BuildError(f"resources: {rel} is {_size(p.stat().st_size)}, above the 2 MB limit")
        group, title, descr = _describe(rel, site, ctitles)
        entries.append({"path": rel, "title": title, "description": descr, "bytes": p.stat().st_size,
                        "group": dict(RES_GROUPS)[group]})
    if total > MAX_TOTAL:
        raise BuildError(f"resources: {_size(total)} in total, above the 12 MB limit")
    order = {g: i for i, (_, g) in enumerate(RES_GROUPS)}
    entries.sort(key=lambda e: (order[e["group"]], 0 if e["path"] in RES_INFO else 1, e["path"]))
    _write_json(site / "data" / "resources.json", entries, indent=1)
    readme = _resources_readme(entries, total)
    (res / "README.md").write_text(readme, encoding="utf-8")
    return [(site / "data" / "resources.json", f"{len(entries)} entries"),
            (res, f"{len(entries) - 1} files, {_size(total)}; {copied} copied, {scrubbed} scrubbed"),
            (res / "README.md", "description of the tree")]


def _resources_readme(entries: list, total: int) -> str:
    lines = ["# Resources of the gamma5 documentation site", "",
             "Light inputs and results behind the site, for download. Every file is at most 2 MB; the directory holds "
             f"{_size(total)} in total. This file and `data/resources.json` are written by `tools/build_site.py`.", "",
             "Absolute local paths, user and host names were removed from the copies: paths are relative to the "
             "directory the file describes (the run directory for the run outputs, the project directory "
             "`one_loop_explicit_calculation/` for the report's working material), and `<scratch>` stands for a "
             "temporary working directory. Nothing else was changed, except that the report's build script now "
             "locates the report next to itself.", "",
             "| directory | content | origin |", "|---|---|---|",
             "| `automation/` | model specification of the automation | `python/hepkit_renorm/hepkit_renorm/models/` of this repository |",
             "| `run_sm_one_family/` | outputs of the end-to-end run of `hepkit_renorm` on that specification: tree-level and renormalised HepKIT models, order store, order data, run report, comparison with the report | the run directory of the automation |",
             "| `final_counterterms.json` | the report's final numbers | `SCRATCHPAD/phase4/` of the project |",
             "| `report_headings.txt` | first page of every labelled heading of the report | the report build |",
             "| `report_classes/<key>/` | final result of each class of Green's functions, `<key>` as in the report's class subsections | `SCRATCHPAD/classes/<key>/` |",
             "| `report_checks/` | assembly of the singular counterterm, pole identity, the two solutions of the restoring condition | `SCRATCHPAD/phase4/` |",
             "| `report_model/` | hand-written HepKIT JSON model of the report and its enumeration summary | `SCRATCHPAD/model/` |",
             "| `toolkit/` | interface and validation notes of the stand-alone toolkit `bmhvkit` | `SCRATCHPAD/toolkit/` |",
             "| `report_src/` | Typst sources of the report, its bibliography and its build script | `report/`, `SCRATCHPAD/scripts/build_report.py` |",
             "", "The report itself is `../report/one_loop_gamma5.pdf`. Third-party documents (papers, preprints) are "
             "not included.", "", "## Files"]
    group = None
    for e in entries:
        if e["group"] != group:
            group = e["group"]
            lines += ["", f"### {group}", ""]
        path = e["path"][len("resources/"):] if e["path"].startswith("resources/") else "../" + e["path"]
        lines.append(f"- `{path}` ({_size(e['bytes'])}): {e['description']}")
    return "\n".join(lines) + "\n"


def step_manifest(site: Path, args) -> list:
    from gdocs import theory

    t = theory.manifest_texts()
    phases = []
    for pid, title in (("unbroken", "Unbroken phase"), ("broken", "Broken phase")):
        f = site / "data" / f"{pid}.json"
        if f.exists() and not _is_mock(f):
            phases.append({"id": pid, "title": title, "file": f"data/{pid}.json"})
        elif pid == "unbroken":
            raise BuildError("manifest: data/unbroken.json is missing or mock")
    pdf = site / "report" / "one_loop_gamma5.pdf"
    if not pdf.exists():
        raise BuildError("manifest: report/one_loop_gamma5.pdf is missing")
    doc = {"title": t["title"], "subtitle": t["subtitle"], "model": "sm_one_family", "loop_order": 1,
           "generated": args.date, "report_pdf": "report/one_loop_gamma5.pdf", "report_pages": _pdf_pages(pdf),
           "phases": phases, "classes": "data/classes.json", "theory": "data/theory.json",
           "resources": "data/resources.json", "status": t["status"], "types": t["types"],
           "conventions": t["conventions"]}
    _write_json(site / "data" / "manifest.json", doc, indent=1)
    return [(site / "data" / "manifest.json", f"{len(phases)} phases, generated {args.date}")]


def _pdf_pages(pdf: Path) -> int:
    """Page count of the PDF from its page tree (no dependency)."""
    data = pdf.read_bytes()
    counts = [int(m) for m in re.findall(rb"/Type\s*/Pages\b[^>]*?/Count\s+(\d+)", data)]
    counts += [int(m) for m in re.findall(rb"/Count\s+(\d+)[^>]*?/Type\s*/Pages\b", data)]
    if not counts:
        return 186
    return max(counts)


def step_validate(site: Path, args) -> list:
    node = shutil.which("node")
    if not node:
        raise BuildError("validate: node not found on PATH (use --skip validate to build without the validator)")
    r = subprocess.run([node, str(HERE / "tests" / "check_site.mjs"), str(site), "--no-mock", "--max", "20"],
                       capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip().splitlines()
    for line in out[-25:]:
        log(line)
    if r.returncode != 0:
        raise BuildError("the contract validator failed")
    summary = next((ln for ln in out if ln.startswith("checked ")), "passed")
    return [(site / "data", summary.split(": ", 1)[-1] if ": " in summary else summary)]


def remove_mock(site: Path) -> list:
    rows = []
    for p in sorted((site / "data").glob("*.json")) + sorted((site / "data" / "graphs").glob("*.json")):
        if _is_mock(p):
            p.unlink()
            rows.append((p, "removed (mock)"))
    return rows


# ================================================================================================ main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="build_site.py", description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", default=str(DEFAULT_SITE), help="site directory (default: gamma5_docs/site)")
    ap.add_argument("--date", default=os.environ.get("SOURCE_DATE") or _dt.date.today().isoformat(),
                    help="'generated' date of the manifest, YYYY-MM-DD (default: $SOURCE_DATE, else today)")
    ap.add_argument("--skip", action="append", default=[], choices=STEPS, metavar="STEP",
                    help=f"skip a step (repeatable): {', '.join(STEPS)}")
    ap.add_argument("--only", action="append", default=[], choices=STEPS, metavar="STEP",
                    help="run only these steps (repeatable)")
    ap.add_argument("--reuse-graphs", action="store_true", help="keep data/classes.json and data/graphs as they are")
    ap.add_argument("--graph-cache", default=os.environ.get("GDOCS_GRAPH_CACHE"),
                    help="directory of rendered graph drawings (default: $GDOCS_GRAPH_CACHE, else none)")
    ap.add_argument("-j", "--processes", type=int, default=None, help="worker processes of the graph build")
    ap.add_argument("--contributions", action="store_true",
                    help="re-evaluate every graph for its contribution flags first (about 7 min on 32 cores)")
    ap.add_argument("--sources", default=None, metavar="PROJECT_DIR",
                    help="refresh the resource copies from this directory (it holds report/ and SCRATCHPAD/)")
    ap.add_argument("--galleries", action="store_true", help="also rewrite the graph and viewer check pages")
    args = ap.parse_args(argv)
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args.date):
        ap.error(f"--date must be YYYY-MM-DD, got {args.date!r}")
    site = Path(args.site).resolve()
    if not (site / "index.html").exists():
        ap.error(f"{site} is not the site directory (no index.html)")
    needs_symbolica = {"viewer", "graphs", "unbroken", "groups", "theory", "broken", "manifest"}
    todo = [s for s in STEPS if s not in args.skip and (not args.only or s in args.only)]
    if needs_symbolica & set(todo) and not os.environ.get("SYMBOLICA_LICENSE"):
        print("build_site: SYMBOLICA_LICENSE is not set in the environment", file=sys.stderr)
        return 1

    funcs = {"viewer": step_viewer, "graphs": step_graphs, "unbroken": step_unbroken, "groups": step_groups,
             "theory": step_theory, "broken": step_broken, "resources": step_resources, "manifest": step_manifest,
             "validate": step_validate}
    table = []
    t_all = time.time()
    for name in todo:
        if name == "validate":
            table += [("mock", p, n) for p, n in remove_mock(site)]
        print(f"[{name}]", flush=True)
        t0 = time.time()
        try:
            rows = funcs[name](site, args)
        except Exception as exc:  # noqa: BLE001  (every failure stops the build)
            traceback.print_exc()
            print(f"\nbuild_site: step '{name}' FAILED: {exc}", file=sys.stderr)
            return 1
        dt = time.time() - t0
        log(f"done in {dt:.1f} s")
        table += [(name, p, n) for p, n in rows]
    if "validate" not in todo:
        table += [("mock", p, n) for p, n in remove_mock(site)]

    print(f"\nwritten ({time.time() - t_all:.0f} s):")
    w = max((len(_rel(site, p)) for _, p, _ in table), default=10)
    for step, p, note in table:
        size = ""
        if p.is_file():
            size = _size(p.stat().st_size)
        elif p.is_dir():
            size = _size(sum(q.stat().st_size for q in p.rglob("*") if q.is_file()))
        print(f"  {step:<10} {_rel(site, p):<{w}}  {size:>10}  {note}")
    return 0


def _rel(site: Path, p: Path) -> str:
    try:
        return p.relative_to(site).as_posix() or "."
    except ValueError:
        return p.name


if __name__ == "__main__":
    sys.exit(main())
