#!/usr/bin/env python3
"""Install gate and publication check of docs/benchmark.pdf.

  receipt.py gate    --out OUT --data-root ROOT [--accept-pending FILE]
      refuses (exit 1, reasons on stdout) to install OUT's results unless they are publishable:
      not a draft; no pending cell (or every pending cell listed in FILE, a user-approved exception:
      {"approved_by": "user", "decision": "<reports/user-decisions.md entry>", "cells": ["<case> <cell>", ...]});
      the measured binary carries its build commit, was built clean (no modified path outside
      crates/vakint/docs), without debug assertions, and the repository HEAD differs from that commit only
      under crates/vakint/docs; the tools and cases that produced the results are those of the tree;
      the data root's manifests are those that were measured; OUT's precision tables (tables.json) are 'final'
      and equal a regeneration from the data root (tools/precision_tables.py --verify).
  receipt.py write   --out OUT --data-root ROOT --pdf PDF --receipt FILE [--accept-pending FILE]
      writes the receipt (identities of data, binary, build, tools, results, graphs, precision tables and PDF).
  receipt.py check   --data-root DATA_TREE [--repo REPO]
      for publish.sh: the committed receipt, suite.json, graphs.json, tables.json, PRECISION-TABLES.md and
      benchmark.pdf agree with each other, with the data tree being published (manifest SHA-256 of the root and of
      every loop pack; the precision tables regenerated from it are identical) and with the repository (build
      commit an ancestor of HEAD, no later change outside crates/vakint/docs, clean work tree). Exit 0 with
      'benchmark receipt OK', else exit 1 with the reasons.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent              # crates/vakint/docs/benchmark
DOCS = BENCH.parent              # crates/vakint/docs
DOCS_PREFIX = "crates/vakint/docs/"


def sha256(path: Path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git(repo, *a):
    r = subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def manifests(root: Path):
    out = {"root": sha256(root / "manifest.json") if (root / "manifest.json").exists() else None}
    for level in range(1, 6):
        m = root / f"{level}L" / "manifest.json"
        if m.exists():
            out[f"{level}L"] = sha256(m)
    return out


def cell_counts(results):
    cells = [(c["id"], k, v) for c in results["cases"] for k, v in c["cells"].items()]
    return {s: [f"{i} {k}" for i, k, v in cells if v["status"] == s] for s in ("ok", "pending", "failed", "n/a")}


def tables_reasons(tables: Path, root: Path, tool_dir: Path = HERE):
    """precision tables: final, and identical to a regeneration from the data root (never stale)"""
    if not tables.exists():
        return [f"no precision tables {tables} (make-benchmark.sh writes them)"]
    reasons = []
    try:
        t = json.loads(tables.read_text())
    except ValueError as e:
        return [f"{tables} is not JSON ({e})"]
    if t.get("status") != "final":
        reasons.append(f"the precision tables are '{t.get('status')}', not 'final' (install from a clean git checkout of the data tree)")
    r = subprocess.run([sys.executable, "-I", str(tool_dir / "precision_tables.py"), "--data-root", str(root), "--verify", str(tables)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        reasons.append("precision tables: " + (r.stdout + r.stderr).strip().splitlines()[-1][:400])
    return reasons


def outside_docs(paths):
    return [p for p in paths if p and not p.startswith(DOCS_PREFIX)]


def gate(out: Path, root: Path, accept: Path | None):
    reasons = []
    results = json.loads((out / "results.json").read_text())
    ident = results["identity"]
    repo = Path(git(BENCH, "rev-parse", "--show-toplevel") or BENCH)
    if results.get("draft"):
        reasons.append("results are a DRAFT (assembled with --draft)")
    counts = cell_counts(results)
    exception = None
    if counts["pending"]:
        if not accept:
            reasons.append(f"{len(counts['pending'])} pending cells (only a user-approved exception may waive them): "
                           + ", ".join(counts["pending"]))
        else:
            exception = json.loads(accept.read_text())
            if exception.get("approved_by") != "user" or not exception.get("decision"):
                reasons.append(f"{accept}: an exception needs approved_by = 'user' and the decision reference")
            missing = sorted(set(counts["pending"]) - set(exception.get("cells", [])))
            if missing:
                reasons.append("pending cells not covered by the user-approved exception: " + ", ".join(missing))
    build = ident.get("build") or {}
    if not build.get("commit"):
        reasons.append("the binary does not carry its build commit (build it with docs/benchmark/build-benchmark-binary.sh)")
    if build.get("dirty_paths") != 0:
        reasons.append(f"the binary was built from a modified tree (dirty paths outside docs: {build.get('dirty_paths')})")
    if build.get("debug_assertions") or ident.get("debug_assertions"):
        reasons.append("the binary has debug assertions (timings need the release profile)")
    if build.get("commit"):
        head = git(repo, "rev-parse", "HEAD")
        if git(repo, "merge-base", "--is-ancestor", build["commit"], "HEAD") is None:
            reasons.append(f"build commit {build['commit'][:12]} is not an ancestor of HEAD {str(head)[:12]}")
        else:
            changed = outside_docs((git(repo, "diff", "--name-only", build["commit"], "HEAD") or "").splitlines())
            if changed:
                reasons.append(f"HEAD changed {len(changed)} paths outside {DOCS_PREFIX} since the measured build, e.g. "
                               + ", ".join(changed[:5]))
        dirty = outside_docs([l[3:] for l in (git(repo, "status", "--porcelain", "--untracked-files=normal") or "").splitlines()])
        if dirty:
            reasons.append(f"the work tree has {len(dirty)} modified paths outside {DOCS_PREFIX}, e.g. " + ", ".join(dirty[:5]))
    tools_now = {p.name: sha256(p) for p in sorted(HERE.glob("*.py"))}
    for name, digest in (ident.get("tools_sha256") or {}).items():
        if tools_now.get(name) != digest:
            reasons.append(f"tools/{name} differs from the version that produced the results")
    if ident.get("cases_sha256") != sha256(BENCH / "cases.json"):
        reasons.append("cases.json differs from the version that produced the results")
    now = manifests(root)
    if ident.get("root_manifest_sha256") != now["root"]:
        reasons.append("the data root's manifest.json differs from the measured one")
    for level, digest in (ident.get("loop_manifest_sha256") or {}).items():
        if now.get(level) != digest:
            reasons.append(f"the data root's {level}/manifest.json differs from the measured one")
    reasons += tables_reasons(out / "tables.json", root)
    return reasons, results, counts, exception


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["gate", "write", "check"])
    ap.add_argument("--out", type=Path)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--pdf", type=Path)
    ap.add_argument("--receipt", type=Path)
    ap.add_argument("--accept-pending", type=Path)
    ap.add_argument("--repo", type=Path)
    a = ap.parse_args(argv)
    root = a.data_root.resolve()
    if a.command in ("gate", "write"):
        reasons, results, counts, exception = gate(a.out.resolve(), root, a.accept_pending)
        if a.command == "gate":
            if reasons:
                print("benchmark install REFUSED:\n- " + "\n- ".join(reasons))
                return 1
            print("benchmark install gate OK")
            return 0
        if reasons:
            print("benchmark receipt REFUSED:\n- " + "\n- ".join(reasons))
            return 1
        ident = results["identity"]
        receipt = {
            "schema": "vakint.benchmark-receipt.v3",
            "installed_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "data_root": str(root), "data_git": ident.get("data_git"),
            "root_manifest_sha256": ident.get("root_manifest_sha256"),
            "loop_manifest_sha256": ident.get("loop_manifest_sha256"),
            "executable_sha256": ident.get("executable_sha256"), "build": ident.get("build"),
            "harness_source_sha256": ident.get("harness_source_sha256"), "vakint_version": ident.get("vakint_version"),
            "tools_sha256": ident.get("tools_sha256"), "cases_sha256": ident.get("cases_sha256"),
            "draft": results.get("draft"), "protocol": results.get("protocol"),
            "cells": {k.replace("/", ""): len(v) for k, v in counts.items()},
            "pending_cells": counts["pending"], "failed_cells": counts["failed"], "pending_exception": exception,
            "suite_sha256": sha256(BENCH / "data" / "suite.json"),
            "graphs_sha256": sha256(BENCH / "data" / "graphs.json"),
            "tables_sha256": sha256(BENCH / "data" / "tables.json"),
            "tables_content_sha256": json.loads((BENCH / "data" / "tables.json").read_text()).get("content_sha256"),
            "tables_md_sha256": sha256(BENCH / "data" / "PRECISION-TABLES.md"),
            "pdf_sha256": sha256(a.pdf),
        }
        a.receipt.write_text(json.dumps(receipt, indent=1) + "\n")
        print(f"wrote {a.receipt}")
        return 0
    # check (publish.sh)
    repo = (a.repo or Path(git(BENCH, "rev-parse", "--show-toplevel") or BENCH)).resolve()
    docs = repo / DOCS_PREFIX
    reasons = []
    rpath = docs / "benchmark" / "data" / "receipt.json"
    if not rpath.exists():
        print(f"benchmark receipt MISSING: {rpath}")
        return 1
    r = json.loads(rpath.read_text())
    if r.get("draft") is not False:
        reasons.append("the receipt is not a final (non-draft) run")
    if r.get("pending_cells") and not (r.get("pending_exception") or {}).get("approved_by") == "user":
        reasons.append("pending cells without a user-approved exception")
    for key, path in (("suite_sha256", docs / "benchmark" / "data" / "suite.json"),
                      ("graphs_sha256", docs / "benchmark" / "data" / "graphs.json"),
                      ("tables_sha256", docs / "benchmark" / "data" / "tables.json"),
                      ("tables_md_sha256", docs / "benchmark" / "data" / "PRECISION-TABLES.md"),
                      ("pdf_sha256", docs / "benchmark.pdf")):
        if not path.exists() or sha256(path) != r.get(key):
            reasons.append(f"{path.relative_to(repo)} does not match the receipt")
    reasons += tables_reasons(docs / "benchmark" / "data" / "tables.json", root, docs / "benchmark" / "tools")
    now = manifests(root)
    if now["root"] != r.get("root_manifest_sha256"):
        reasons.append("the data tree's manifest.json is not the one the benchmark measured")
    for level, digest in (r.get("loop_manifest_sha256") or {}).items():
        if now.get(level) != digest:
            reasons.append(f"the data tree's {level}/manifest.json is not the one the benchmark measured")
    build = (r.get("build") or {}).get("commit")
    if not build:
        reasons.append("the receipt has no build commit")
    elif git(repo, "merge-base", "--is-ancestor", build, "HEAD") is None:
        reasons.append(f"build commit {build[:12]} is not an ancestor of HEAD")
    else:
        changed = outside_docs((git(repo, "diff", "--name-only", build, "HEAD") or "").splitlines())
        if changed:
            reasons.append(f"{len(changed)} paths outside {DOCS_PREFIX} changed after the measured build, e.g. "
                           + ", ".join(changed[:5]))
    if git(repo, "status", "--porcelain", "--untracked-files=no"):
        reasons.append("the repository work tree is not clean")
    if reasons:
        print("benchmark receipt FAILED:\n- " + "\n- ".join(reasons))
        return 1
    print("benchmark receipt OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
