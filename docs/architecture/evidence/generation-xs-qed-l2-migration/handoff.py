"""Derived migration receipt, scoped canonical docs, and source/evidence seal."""

import argparse
import datetime
import hashlib
import os
import re
import sys
import tomllib
from pathlib import Path

import common as M
import migrate as X

sys.path.append(str(M.FROZEN))
import docs as DOC

NOTE = M.WORKSPACE / "docs/architecture/generation-xs-qed-l2-migration.typ"
REGISTRY = M.WORKSPACE / "docs/developers.toml"


def receipt():
    runtime = M.load(M.OUT / "runtime/results.json")
    M.require(
        len(runtime) == 6 and len(M.load(M.OUT / "binding.json")["domains"]) == 12,
        "receipt",
        "complete actual matrix",
    )
    old_controls = {
        name: {
            "path": str(OLD / name / "results.json"),
            "sha256": M.sha(OLD / name / "results.json"),
            "rows": len(M.load(OLD / name / "results.json")),
        }
        for name in ("controls", "supplements", "native-parser-controls")
    }
    result = {
        "schema": "xs-qed-l2-migration/1",
        "status": "MIGRATION_PREPARATION_CANDIDATE_REQUIRES_NEW_INDEPENDENT_REVIEW",
        "scope": {
            "case": X.CASE,
            "literal_only": True,
            "asserted_reference": "Classes::None returned by classes::check; Auto checked separately",
            "grouping": "None",
            "runtime_pools": [1, 8, 32],
            "historical_mathematical_pool": 1,
            "no_MAIN_integration": True,
        },
        "approved_certificate": M.APPROVED,
        "approvals": {k: M.load(M.OUT / "approvals" / (k + ".json")) for k in X.REVIEW},
        "code": M.load(M.OUT / "code-audit.json"),
        "derivation": M.load(M.OUT / "derivation.json"),
        "runtime": runtime,
        "binding": {
            "path": str(M.OUT / "binding.json"),
            "sha256": M.sha(M.OUT / "binding.json"),
            "domains": M.load(M.OUT / "binding.json")["domains"],
        },
        "dependency_closure": M.load(M.OUT / "build/closure-comparison.json"),
        "builds": {
            name: {
                "qualification": str(M.OUT / "build" / name / "qualification-v2.json"),
                "sha256": M.sha(M.OUT / "build" / name / "qualification-v2.json"),
                "elf": M.load(M.OUT / "build" / name / "qualification-v2.json")[
                    "elf_sha256"
                ],
            }
            for name in ("generation_oracle", "oracle_xs_qed_l2_migration_capture")
        },
        "replay": {
            mode: {
                "summary": M.load(M.OUT / "replay" / mode / "summary.json"),
                "verdict": M.load(M.OUT / "replay" / mode / "verdict.json"),
                "cross_pairs_sha256": M.sha(
                    M.OUT / "replay" / mode / "cross-pairs.json"
                ),
                "execution_sha256": M.sha(M.OUT / "replay" / mode / "execution.json"),
            }
            for mode in ("None", "Auto")
        },
        "controls": {
            "reused_frozen_suites": old_controls,
            "reviewer_B_heldouts": {
                "path": str(
                    M.WORKSPACE / X.REVIEW["B"][0] / "heldouts/results-final.json"
                ),
                "sha256": "79791f15956a48a10de8821075ede8d6b638fdc4476191746a0b5b3ca8941ca9",
                "scope": "104 reviewer runs; not permanent candidate suite coverage",
            },
            "fresh_runtime_negatives": M.load(M.OUT / "negatives/results.json"),
        },
        "limitations": [
            "Named repeat-free axis multiset transport only; no positional consumer, free-index typing or factor locality theorem.",
            "Interface presence and typed-port order are frozen-anchor integrity; baseline anchors self-compare. New runtime compares every emitted mode to its own approved full payload and ordered owner capture.",
            "Historical raw reversible-cook association remains capture-bound; reviewer B independently corroborated a graph-local injective renaming, not a universal decoder proof.",
            "Assertion pool delivery inferred from actual env, source and dump filenames. Observed Run/request/model/options are separate observer-process facts joined by complete payload and source companions, not assertion in-memory instrumentation.",
            "Only Q and edge maps exercised; K/P/vertex absent. Exact post-map definitions and source-audited spelling changes are load-bearing; no general pre-map consistency theorem or native tagged-equivalence theorem.",
            "Native result hash identifies count, not pairs; actual pair-file hashes identify inputs. Both native kernels may tolerate definition subsets without Python gates; real license passed privately, validity not used as equality evidence.",
            "Committed 82+18 controls and 4 parser probes are reused by exact source/raw bindings; B held-outs supply additional predicate coverage. Only applicable baselines and new runtime comparisons rerun here.",
            "New pools8/32 are current runtime observations only, not historical mathematical certification. No other row, grouping, physics, dependency, MAIN or expectation migration beyond one literal.",
            "Older helper/build command gaps remain as attributed by A and B. A independently bound exact336/341 observer/assertion units and actual kernels; no retroactive command/hash claims. Initial migration setup and first build helper versions are retained distinctly.",
            "Full frozen source-inventory comparison rejects only the authorized registry insertion with zero credit. This is not the full outer validator. Portable finite validate consumes exact executable/source/capture pins and required manifests without changing frozen guards.",
            "Full docs not green: exact inherited first-gate list retained; unknown sunset owners, then unclassified source46, then unexercised catalog exporter remain later gates.",
            "Archived .rs sources remain outside Cargo auto targets but affect broad commonCargoSources/build/fmt/test/WASM/docs filters; no new derivation/NAR source-invariance claim.",
        ],
        "chronology": {
            "preparatory_commands": "Initial guidance, new-child, two failed imports, derive/audit/prepare and first build launcher have tool transcripts but no new top-level command.json. Nested source/archive/build commands are actually recorded. Later stages use the self-recording launcher; no retroactive transcript invention.",
            "preservation_before": "Started before one-literal edit, completed concurrently with build; immutable prior seals and per-command guards additionally bind prior evidence.",
            "license_Nix": "All new Nix commands use explicit validated nonempty immutable Git URI; both license variables unset. No metadata/derivation-show, private-source flake, default installable or legacy licensed result invocation. No incident copies inspected or removed.",
        },
    }
    for mode in result["replay"].values():
        s = mode["summary"]
        M.require(
            s["owners"] == 540
            and s["nontrivial_cross_permutations"] == 56
            and s["distinct_native_texts"] == 168,
            "receipt",
            "finite independent domains",
        )
    M.require(
        M.sha(result["controls"]["reviewer_B_heldouts"]["path"])
        == result["controls"]["reviewer_B_heldouts"]["sha256"],
        "heldouts",
        "required exact result",
    )
    M.write(M.HERE / "receipt.json", result)
    print("receipt", M.sha(M.HERE / "receipt.json"), flush=True)


OLD = X.OLD


def checks():
    root = M.OUT / "source"
    for label, argv in (
        ("fmt", ["rustfmt", "--edition", "2024", "--check", root / M.MAIN]),
        (
            "clippy",
            [
                "cargo",
                "clippy",
                "--locked",
                "--profile",
                "dev-optim",
                "-p",
                "feynkit-generator",
                "--test",
                "generation_oracle",
                "--",
                "--deny",
                "warnings",
            ],
        ),
    ):
        p, _ = M.run(
            ["nix", "develop", M.URI, "-c", *argv],
            M.OUT / "checks" / label,
            cwd=root,
            inputs=[root / M.MAIN, M.WORKSPACE / M.MAIN],
        )
        M.require(p.returncode == 0, "checks", label)


def guards():
    expected = {
        "generation-xs-qed-l2-candidate": M.APPROVED,
        "generation-oracle-family-migration": "3432ac3f36890ae2664d8571085f05e1dd41f27f",
        "generation-oracle-family-approved": "87fa13a93b1c0b3e5792fd56937486ad962b2b55",
        "generation-oracle-family-round2": "78c0754c8fb732878123c02999ced3fc506c5a1d",
        "generation-oracle-family-round1": "e814b6b4c3463d9aa1c018939365b624619af24b",
        "generation-oracle-family-historical-diagnostic": "5e61ef89f4df4e818884e40b8d0c46108eab75a8",
    }
    for name, target in expected.items():
        p, _ = M.run(
            [
                "jj",
                "log",
                "--ignore-working-copy",
                "--no-graph",
                "-r",
                name,
                "-T",
                "commit_id",
            ],
            M.OUT / "guards" / name,
        )
        M.require(p.returncode == 0 and p.stdout.strip() == target, "frozen-ref", name)
    p, _ = M.run(
        [
            "git",
            "--git-dir=/shared/lcnbr/dev/gammaloop/main/.git",
            "ls-tree",
            "-rz",
            "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
        ],
        M.OUT / "guards/historical-tree",
    )
    M.require(p.returncode == 0, "historical", "immutable tree inventory")
    root = M.WORKSPACE.parent / "generation-oracle-family-historical"
    files, gitlinks = {}, []
    for line in p.stdout.split("\0"):
        if not line:
            continue
        meta, path = line.split("\t", 1)
        mode, kind, blob = meta.split()
        if kind == "commit":
            gitlinks.append(path)
            continue
        file = root / path
        value = os.readlink(file).encode() if mode == "120000" else file.read_bytes()
        actual = hashlib.sha1(
            b"blob " + str(len(value)).encode() + b"\0" + value, usedforsecurity=False
        ).hexdigest()
        M.require(actual == blob, "historical-pristine", path)
        files[path] = {"blob": blob, "mode": mode}
    M.write(
        M.OUT / "guards/result.json",
        {
            "frozen_refs": expected,
            "historical_source": "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
            "historical_actual_blobs": files,
            "gitlinks_not_regular_blobs": gitlinks,
            "scope": "post-run direct content guard; no retroactive before claim; MAIN files/index never mutated",
        },
    )


def docs():
    table = tomllib.loads(REGISTRY.read_text())
    notes = [n for s in table["section"] for n in s["note"]]
    M.require(table["schema"] == 3, "docs", "schema3")
    for key in ("id", "source", "title"):
        M.require(len({n[key] for n in notes}) == len(notes), "docs", key)
    entry = next(n for n in notes if n["id"] == "generation-xs-qed-l2-migration")
    M.require(
        re.findall(r"^= (.+)$", NOTE.read_text(), re.MULTILINE) == [entry["title"]],
        "docs",
        "canonical title",
    )
    original = X.immutable(M.APPROVED, "docs/developers.toml")
    inserted = REGISTRY.read_text().replace(original, "")
    # The whole original registry is split only once by the inserted block.
    block = re.search(
        r'\[\[section.note\]\]\nid = "generation-xs-qed-l2-migration".*?(?=\[\[section.note\]\])',
        REGISTRY.read_text(),
        re.DOTALL,
    )[0]
    M.require(
        REGISTRY.read_text().replace(block, "", 1) == original,
        "registry-preservation",
        "all prior entry bytes",
    )
    del inserted
    missing = {
        str(p.relative_to(M.WORKSPACE))
        for p in (M.WORKSPACE / "docs/architecture").glob("*.typ")
    } - {n["source"] for n in notes}
    M.require(
        missing == {"docs/architecture/generation-integration-source46.typ"},
        "docs",
        "only inherited unregistered note",
    )
    root = M.OUT / "docs"
    root.mkdir(parents=True, exist_ok=True)
    wrapper = root / "main.typ"
    wrapper.write_text(
        '#document("note.html", title: ['
        + entry["title"]
        + '])[#include "/docs/architecture/generation-xs-qed-l2-migration.typ"]\n'
    )
    inputs = [NOTE, REGISTRY, M.HERE / "receipt.json", wrapper, M.HERE / "handoff.py"]
    commands = []
    for label, argv in (
        (
            "pdf",
            [
                DOC.TYPST,
                "compile",
                "--root",
                ".",
                "--creation-timestamp",
                "0",
                NOTE,
                root / "report.pdf",
            ],
        ),
        (
            "html",
            [
                DOC.TYPST,
                "compile",
                "--root",
                ".",
                "--features",
                "html,bundle",
                "--format",
                "bundle",
                "--creation-timestamp",
                "0",
                wrapper,
                root / "bundle",
            ],
        ),
    ):
        p, cmd = M.run(argv, root / label, inputs=inputs)
        M.require(p.returncode == 0, "docs-render", label)
        commands.append(cmd)
    html = root / "bundle/note.html"
    parser = DOC.SemanticHTML()
    parser.feed(html.read_text())
    M.require(
        parser.title[0] == entry["title"]
        and parser.title.count(entry["title"]) == 1
        and parser.headers == 2
        and not parser.active,
        "semantic-HTML",
        "title and semantic tables",
    )
    for text in ("17909042738907132227", "17540448598867080140", "714", "157"):
        M.require(
            text in html.read_text(), "semantic-HTML", "actual Row integer " + text
        )
    altered = root / "altered.typ"
    altered.write_text(
        NOTE.read_text().replace("17909042738907132227", "17909042738907132228")
    )
    # Preserve relative receipt import by binding the private variant only over
    # the note path in a no-network read-only-root sandbox.
    p, altered_cmd = M.run(
        [
            "bwrap",
            "--ro-bind",
            "/",
            "/",
            "--dev",
            "/dev",
            "--proc",
            "/proc",
            "--bind",
            root,
            root,
            "--ro-bind",
            altered,
            NOTE,
            "--unshare-net",
            "--die-with-parent",
            "--",
            DOC.TYPST,
            "compile",
            "--root",
            ".",
            "--creation-timestamp",
            "0",
            NOTE,
            root / "altered.pdf",
        ],
        root / "altered-command",
        inputs=[*inputs, altered],
    )
    M.require(
        p.returncode == 0 and M.sha(root / "altered.pdf") != M.sha(root / "report.pdf"),
        "altered-text",
        "actual rendered text changes PDF",
    )
    builder = M.WORKSPACE / "target/debug/alphal00p-docs-builder"
    p, builder_cmd = M.run(
        [
            "bwrap",
            "--ro-bind",
            "/",
            "/",
            "--dev",
            "/dev",
            "--proc",
            "/proc",
            "--tmpfs",
            M.WORKSPACE / "target",
            "--ro-bind",
            builder,
            builder,
            "--bind",
            root,
            root,
            "--unshare-net",
            "--die-with-parent",
            "--",
            builder,
            "check",
        ],
        root / "builder",
        inputs=[*inputs, builder],
    )
    expected = (OLD / "docs/builder/stderr.txt").read_text()
    M.require(
        p.returncode == 1 and p.stderr == expected,
        "full-docs-first-gate",
        "exact nine inherited paths, no inclusion-only acceptance",
    )
    M.write(
        root / "validation.json",
        {
            "entry": entry,
            "commands": commands,
            "altered_text_command": altered_cmd,
            "builder_command": builder_cmd,
            "full_docs_exact_stderr": p.stderr,
            "full_docs": "not green; first gate exact9; masked unknown sunset owners then source46 then unexercised catalogs, independently reviewed A",
            "pdf_sha256": M.sha(root / "report.pdf"),
            "html_sha256": M.sha(html),
            "altered_pdf_sha256": M.sha(root / "altered.pdf"),
            "semantic_headers": parser.headers,
        },
    )


def seal():
    X.audit()
    p, _ = M.run(
        ["jj", "diff", "--summary", "-r", "@"], M.OUT / "inventory/changed-path-command"
    )
    M.require(p.returncode == 0, "diff", "actual changed paths")
    changed = [line[2:] for line in p.stdout.splitlines()]
    allowed = [
        M.MAIN,
        "docs/developers.toml",
        str(NOTE.relative_to(M.WORKSPACE)),
        *[str(p.relative_to(M.WORKSPACE)) for p in M.HERE.iterdir() if p.is_file()],
    ]
    M.require(set(changed) == set(allowed), "changed-paths", str(changed))
    M.write(
        M.OUT / "inventory/changed-paths.json",
        {
            "base": M.APPROVED,
            "paths": changed,
            "code_paths": [M.MAIN],
            "other_paths": "new migration evidence/note and selective registry insertion only",
        },
    )
    own = [p for p in M.HERE.rglob("*") if p.is_file()] + [
        NOTE,
        REGISTRY,
        M.WORKSPACE / M.MAIN,
    ]
    sources = X.manifest(own)
    M.write(M.OUT / "inventory/sources.json", sources)
    paths = [
        p
        for p in M.OUT.rglob("*")
        if p.is_file()
        and not any(x in p.relative_to(M.OUT).parts for x in ("cache", "tmp", "state"))
        and not p.is_relative_to(M.OUT / "inventory")
        and p.name != "freeze.json"
        and not p.is_relative_to(M.OUT / "commands")
    ]
    artifacts = X.manifest(paths)
    M.write(
        M.OUT / "inventory/artifacts.json",
        {
            "files": artifacts,
            "excluded": "private cache/tmp/state; self inventories/freeze; command records separately closed",
        },
    )
    commands = [
        p
        for p in (M.OUT / "commands").rglob("*")
        if p.is_file()
        and p.name != "command.json"
        or p.is_file()
        and p.name == "command.json"
        and "exit" in M.load(p)
    ]
    M.write(M.OUT / "inventory/completed-commands.json", X.manifest(commands))
    roots = [
        M.OUT / "inventory/sources.json",
        M.OUT / "inventory/artifacts.json",
        M.OUT / "inventory/completed-commands.json",
        M.OUT / "preservation/after.json",
        M.OUT / "docs/validation.json",
    ]
    M.write(
        M.OUT / "freeze.json",
        {
            "status": "candidate; stop for NEW independent migration review",
            "approved": M.APPROVED,
            "roots": X.manifest(roots),
            "historical_holder": "read-only pristine, no edits by this task",
            "no_MAIN": True,
            "counts": {"source_paths": len(sources), "artifact_files": len(artifacts)},
        },
    )
    print("freeze", M.sha(M.OUT / "freeze.json"), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("receipt", "checks", "docs", "guards", "seal"))
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()
    if not args.execute:
        stamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%dT%H%M%S%f")
        p, _ = M.run(
            [sys.executable, "-B", Path(__file__), args.stage, "--execute"],
            M.OUT / "commands" / (stamp + "-" + args.stage),
        )
        print(p.stdout, end="", flush=True)
        print(p.stderr, end="", file=sys.stderr, flush=True)
        sys.exit(p.returncode)
    globals()[args.stage]()


if __name__ == "__main__":
    main()
