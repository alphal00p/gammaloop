#!/usr/bin/env python3
"""Install gate and publication check of docs/benchmark.pdf.

  receipt.py gate    --out OUT --data-root ROOT [--accept-pending FILE]
      refuses (exit 1, reasons on stdout) to install OUT's results unless they are publishable:
      not a draft; no pending cell (or every pending cell listed in FILE, a user-approved exception:
      {"approved_by": "user", "decision": "<reports/user-decisions.md entry>", "cells": ["<case> <cell>", ...]});
      the measured binary carries its build commit, was built clean (no modified path outside
      crates/vakint/docs), without debug assertions, and the repository HEAD differs from that commit only
      in neutral paths (rule 1); the tools and cases that produced the results are those of the tree;
      the data root's manifests are those that were measured; OUT's precision tables (tables.json) are 'final'
      and equal a regeneration from the data root (tools/precision_tables.py --verify).
  receipt.py write   --out OUT --data-root ROOT --pdf PDF --receipt FILE [--accept-pending FILE]
      writes the receipt (identities of data, binary, build, tools, results, graphs, precision tables and PDF).
  receipt.py check   --data-root DATA_TREE [--repo REPO]
      for publish.sh: the committed receipt, suite.json, graphs.json, tables.json, PRECISION-TABLES.md and
      benchmark.pdf agree with each other, with the data tree being published and with the repository (rules
      below). Exit 0 with 'benchmark receipt OK', else exit 1 with the reasons.
  receipt.py refresh --out OUT --data-root DATA_TREE (--pdf PDF --receipt FILE | --dry-run)
      make-benchmark.sh --refresh-tables --install (a data release after the timing run): the committed timing
      receipt (docs/benchmark/data/receipt.json) stays; refuses unless DATA_TREE's reduction payloads equal the measured tree's (rule 2), OUT's
      regenerated reduction graphs equal the committed ones and OUT's tables are final and regenerate from
      DATA_TREE; then records the tables' data tree, the new table / PDF hashes and the refresh in the receipt.

Rules (receipt schema v4: timings measured once on the release line, later data releases refresh the tables):
  1. Source. The timings stay valid while the repository differs from the measured build only in paths that do
     not enter the timed calls (NEUTRAL_PATHS): crates/vakint/docs/**, crates/vakint/external_data.typ,
     crates/vakint/AGENTS.md, nix/rust-workspace.nix and the vendored snapshot crates/vakint/data/** except its reduction payloads
     crates/vakint/data/<L>L/ibps/** (every timed call loads all loop packs from --data-root, so the vendored
     snapshot is not read; a pin commit that changes the vendored IBP payloads still refuses). Checked by the
     source digest (SHA-256 of `git ls-tree -r` without the neutral paths) of the build commit, recorded in the
     receipt, against HEAD, plus no modified non-neutral path in the work tree. The build commit need not be an
     ancestor of HEAD (a rebase onto a later pin commit keeps the digest); the check reports whether it is.
  2. Data. The data tree is the measured one (manifest SHA-256 of the root and of every loop pack), or differs
     from it only outside the reduction payloads: per loop pack, the payload digest covers the manifest's
     structure (families, master ids / powers / names, IBP / cache / route-preparation records) and every
     listed file except master values and documentation (VALUE_FILES), each re-hashed on disk against the
     manifest. A new release that changes only master values therefore needs only --refresh-tables.
  3. Tables. tables.json is 'final' and regenerates identically from the data tree (precision_tables.py).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent              # crates/vakint/docs/benchmark
DOCS = BENCH.parent              # crates/vakint/docs
DOCS_PREFIX = "crates/vakint/docs/"
RECEIPT_SCHEMA = "vakint.benchmark-receipt.v4"
# Rule 1: repository paths that do not enter a timed call (see the module docs).
NEUTRAL_PATHS = re.compile(r"^(crates/vakint/docs/.+|crates/vakint/external_data\.typ|crates/vakint/AGENTS\.md|nix/rust-workspace\.nix"
                           r"|crates/vakint/data/(?![0-9]+L/ibps/).+)$")
# Rule 2: files of a loop pack that are master values or documentation, not reduction payloads.
VALUE_FILES = re.compile(r"^(masters/[^/]+/(numerical|symbolic)\.json|constants\.json|coefficient-supplement\.json"
                         r"|numerical-precision-supplement\.json|coverage\.json|producer-update\.json|README\.md"
                         r"|ibps/README\.md|provenance/(?!refinement-layer\.json$).+)$")
PAYLOAD_MANIFEST_KEYS = ("schema", "loop_count", "dimension", "normalization", "unit_mass_squared", "families", "ibps",
                         "cache", "ordered_cache", "route_preparation")
PAYLOAD_MASTER_KEYS = ("master_id", "family_sha256", "loop_count", "powers", "name", "alias", "kind", "parent_graph")


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


def payload_digest(root: Path):
    """Rule 2: per loop pack, SHA-256 of its reduction payload record; (digests, problems). Every payload file is
    re-hashed on disk against the manifest (a tree whose files differ from its manifest is a problem)."""
    out, problems = {}, []
    for level in range(1, 6):
        m = root / f"{level}L" / "manifest.json"
        if not m.exists():
            continue
        d = json.loads(m.read_text())
        rec = {k: d.get(k) for k in PAYLOAD_MANIFEST_KEYS}
        rec["masters"] = [{k: x.get(k) for k in PAYLOAD_MASTER_KEYS} for x in d.get("masters", [])]
        files = {}
        for path, meta in sorted((d.get("files") or {}).items()):
            if VALUE_FILES.match(path):
                continue
            digest = (meta or {}).get("sha256")
            f = root / f"{level}L" / path
            if not f.is_file():
                problems.append(f"{level}L/{path} listed in the manifest is missing")
            elif sha256(f) != digest:
                problems.append(f"{level}L/{path} differs from its manifest sha256")
            files[path] = digest
        rec["files"] = files
        out[f"{level}L"] = hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest()
    return out, problems


def source_digest(repo, commit):
    """Rule 1: SHA-256 of `git ls-tree -r COMMIT` without the neutral paths (None if COMMIT is not in REPO)."""
    tree = git(repo, "ls-tree", "-r", "--full-tree", commit)
    if tree is None:
        return None
    keep = [l for l in tree.splitlines() if not NEUTRAL_PATHS.match(l.split("\t", 1)[1])]
    return hashlib.sha256("\n".join(keep).encode()).hexdigest()


def non_neutral(paths):
    return [p for p in paths if p and not NEUTRAL_PATHS.match(p)]


def dirty_paths(repo):
    """Modified, staged or untracked non-neutral paths of the work tree (porcelain v1 -z: 'XY path', renames add the
    source as the next entry)."""
    r = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "-z", "--untracked-files=normal"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return [f"(git status failed in {repo})"]
    entries, paths, skip = r.stdout.split("\0"), [], False
    for e in entries:
        if skip:
            skip = False
            continue
        if len(e) > 3:
            paths.append(e[3:])
            skip = e[0] in "RC"
    return non_neutral(paths)


def source_reasons(repo, build_commit, recorded_digest=None):
    """Rule 1 for HEAD of REPO: (reasons, note)."""
    reasons = []
    head = git(repo, "rev-parse", "HEAD")
    digest = recorded_digest or source_digest(repo, build_commit)
    if digest is None:
        return [f"build commit {build_commit[:12]} is not in {repo} and the receipt records no source digest"], ""
    if source_digest(repo, "HEAD") != digest:
        changed = (non_neutral((git(repo, "diff", "--name-only", build_commit, "HEAD") or "").splitlines())
                   if git(repo, "cat-file", "-e", f"{build_commit}^{{commit}}") is not None else [])
        reasons.append(f"HEAD {str(head)[:12]} differs from the measured build {build_commit[:12]} outside the neutral paths"
                       + (" (" + ", ".join(changed[:5]) + ")" if changed else ""))
    dirty = dirty_paths(repo)
    if dirty:
        reasons.append(f"the work tree has {len(dirty)} modified non-neutral paths, e.g. " + ", ".join(dirty[:5]))
    in_history = git(repo, "merge-base", "--is-ancestor", build_commit, "HEAD") is not None
    note = ("build commit in the history of HEAD" if in_history else
            "build commit not in the history of HEAD (rebased line); source tree equal outside the neutral paths")
    return reasons, note


def data_reasons(root: Path, measured_root_sha, measured_loops, measured_payload):
    """Rule 2: (reasons, note)."""
    now = manifests(root)
    same = now["root"] == measured_root_sha and all(now.get(k) == v for k, v in (measured_loops or {}).items())
    if same:
        return [], "data tree = the measured tree"
    if not measured_payload:
        return ["the data tree's manifests are not the measured ones and the receipt records no payload digest"], ""
    digests, problems = payload_digest(root)
    reasons = list(problems[:5])
    for level, digest in measured_payload.items():
        if digests.get(level) != digest:
            reasons.append(f"the data tree's {level} reduction payloads differ from the measured ones (re-run the benchmark)")
    return reasons, "data tree differs from the measured one only in master values / documentation (rule 2)"


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
        reasons += source_reasons(repo, build["commit"])[0]
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


def git_info(path: Path):
    if git(path, "rev-parse", "--is-inside-work-tree") != "true":
        return None
    return {"commit": git(path, "rev-parse", "HEAD"), "describe": git(path, "describe", "--tags", "--always", "--dirty"),
            "remote": git(path, "config", "--get", "remote.origin.url"),
            "dirty_files": len([l for l in (git(path, "status", "--porcelain") or "").splitlines() if l.strip()])}


def refresh(out: Path, root: Path):
    """make-benchmark.sh --refresh-tables --install: tables of a new data release, timings kept."""
    reasons = []
    rpath = BENCH / "data" / "receipt.json"
    if not rpath.exists():
        return [f"no committed receipt {rpath} (measure first: make-benchmark.sh ... --install)"], None
    r = json.loads(rpath.read_text())
    if r.get("draft") is not False:
        reasons.append("the committed receipt is not a final timing run")
    if not r.get("data_payload_sha256"):
        reasons.append(f"the committed receipt ({r.get('schema')}) records no payload digest: re-run the benchmark")
    else:
        reasons += data_reasons(root, r.get("root_manifest_sha256"), r.get("loop_manifest_sha256"),
                                r.get("data_payload_sha256"))[0]
    committed = json.loads((BENCH / "data" / "graphs.json").read_text())
    regen = json.loads((out / "graphs.json").read_text())
    if committed.get("loops") != regen.get("loops"):
        reasons.append("the reduction graphs of the data tree differ from the committed ones (re-run the benchmark)")
    if sha256(BENCH / "data" / "suite.json") != r.get("suite_sha256"):
        reasons.append("docs/benchmark/data/suite.json does not match the committed receipt")
    reasons += tables_reasons(out / "tables.json", root)
    return reasons, r


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["gate", "write", "check", "refresh"])
    ap.add_argument("--out", type=Path)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--pdf", type=Path)
    ap.add_argument("--receipt", type=Path)
    ap.add_argument("--accept-pending", type=Path)
    ap.add_argument("--repo", type=Path)
    ap.add_argument("--dry-run", action="store_true", help="refresh: the checks only, nothing written")
    a = ap.parse_args(argv)
    root = a.data_root.resolve()
    if a.command == "refresh":
        reasons, r = refresh(a.out.resolve(), root)
        if reasons:
            print("benchmark refresh REFUSED:\n- " + "\n- ".join(reasons))
            return 1
        if a.dry_run:
            print("benchmark refresh gate OK")
            return 0
        tables = json.loads((BENCH / "data" / "tables.json").read_text())
        stamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        r.update({
            "tables_data_root": str(root), "tables_data_git": git_info(root),
            "tables_root_manifest_sha256": manifests(root)["root"],
            "tables_sha256": sha256(BENCH / "data" / "tables.json"),
            "tables_content_sha256": tables.get("content_sha256"),
            "tables_md_sha256": sha256(BENCH / "data" / "PRECISION-TABLES.md"),
            "pdf_sha256": sha256(a.pdf),
        })
        r.setdefault("refreshes", []).append({"utc": stamp, "data_git": r["tables_data_git"],
                                              "root_manifest_sha256": r["tables_root_manifest_sha256"],
                                              "tables_content_sha256": r["tables_content_sha256"]})
        a.receipt.write_text(json.dumps(r, indent=1) + "\n")
        print(f"refreshed {a.receipt} (tables from {root})")
        return 0
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
        repo = Path(git(BENCH, "rev-parse", "--show-toplevel") or BENCH)
        payload, _ = payload_digest(root)
        receipt = {
            "schema": RECEIPT_SCHEMA,
            "installed_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "data_root": str(root), "data_git": ident.get("data_git"),
            "root_manifest_sha256": ident.get("root_manifest_sha256"),
            "loop_manifest_sha256": ident.get("loop_manifest_sha256"),
            "data_payload_sha256": payload,
            "executable_sha256": ident.get("executable_sha256"), "build": ident.get("build"),
            "build_source_digest": source_digest(repo, (ident.get("build") or {}).get("commit") or "HEAD"),
            "neutral_paths": NEUTRAL_PATHS.pattern, "value_files": VALUE_FILES.pattern,
            "host": {k: (ident.get("host") or {}).get(k) for k in ("hostname", "cpu_model", "logical_cpus", "cpus_used")},
            "contention": results.get("contention"),
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
    notes = []
    if r.get("schema") != RECEIPT_SCHEMA:
        reasons.append(f"receipt schema {r.get('schema')} (this check needs {RECEIPT_SCHEMA}: re-install the benchmark)")
    more, note = data_reasons(root, r.get("root_manifest_sha256"), r.get("loop_manifest_sha256"), r.get("data_payload_sha256"))
    reasons += more
    notes.append(note)
    build = (r.get("build") or {}).get("commit")
    if not build:
        reasons.append("the receipt has no build commit")
    elif not r.get("build_source_digest"):
        reasons.append("the receipt records no build source digest")
    else:
        more, note = source_reasons(repo, build, r["build_source_digest"])
        reasons += more
        notes.append(note)
    if reasons:
        print("benchmark receipt FAILED:\n- " + "\n- ".join(reasons))
        return 1
    print("benchmark receipt OK (" + "; ".join(n for n in notes if n) + ")")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
