"""Freeze an auditable summary after all fresh semantic processes have exited."""

import collections
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tomllib

SPEC = importlib.util.spec_from_file_location(
    "certificate", Path(__file__).with_name("validate.py")
)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
HERE = Path(__file__).parent
OUT = Path(
    os.environ.get(
        "FAMILY_ROUND2_OUT", str(V.WORKSPACE / "target/oracle-numerator-family-round2")
    )
)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(data)
    return digest.hexdigest()


def main():
    controls = V.load(OUT / "controls-certified/results.json")
    V.require(
        all(row["matched"] for row in controls), "final-controls", "unmatched intent"
    )
    per_domain = collections.Counter((row["case"], row["mode"]) for row in controls)
    V.require(
        set(per_domain) == {(case, mode) for case in V.CASES for mode in V.MODES},
        "final-controls",
        "missing domain",
    )
    expected_counts = {
        "qed_eemm_L1": 40,
        "ctorder_qcd_uugg_L1": 42,
        "sm_eemm_L0": 39,
        "ym_g6_L0": 41,
    }
    V.require(
        all(count == expected_counts[case] for (case, _), count in per_domain.items()),
        "final-controls",
        "incomplete per-domain mutation matrix",
    )
    V.require(
        len({(r["case"], r["mode"], r["mutation"]) for r in controls}) == len(controls),
        "final-controls",
        "duplicate mutation result",
    )
    source = sha(HERE / "validate.py")
    for row in controls:
        dest = (
            OUT
            / "controls-certified"
            / row["case"]
            / row["mode"]
            / row["mutation"]
            / "result"
        )
        V.require(
            V.load(dest / "execution.json")["validator_sha256"] == source,
            "source-identity",
            str(dest),
        )
    baseline = V.load(OUT / "replay-final/baseline.json")
    V.require(
        len(baseline) == 8
        and all(
            row["exit"] == 0 and row["verdict"]["disposition"] == "accept"
            for row in baseline
        ),
        "baseline",
        "incomplete",
    )
    families = []
    for case in V.CASES:
        modes = {}
        for mode in V.MODES:
            root = OUT / "replay-final/baseline" / case / mode
            V.require(
                V.load(root / "execution.json")["validator_sha256"] == source,
                "source-identity",
                str(root),
            )
            structural = V.load(root / "structural.json")
            cross = [
                p
                for p in V.load(root / "pairs.json")
                if p["label"].startswith("cross/") and p["historical"] != p["current"]
            ]
            V.require(all(p["cyclic_f_literal"] for p in cross), "cyclic-f", case)
            modes[mode] = {
                "graphs": structural["graphs"],
                "bound_pairs": structural["pairs"],
                "symbol_definitions": len(structural["definitions"]),
                "execution": V.load(root / "execution.json"),
                "controls": per_domain[(case, mode)],
                "nonliteral_cross_owner_instances": len(cross),
                "nonliteral_cross_instances_literal_after_even_f_rotation": len(cross),
            }
        families.append(
            {
                "case": case,
                "status": "CORRECTED_CANDIDATE_PENDING_NEW_REVIEW",
                "modes": modes,
            }
        )
    b1 = V.load(OUT / "replay-final/B1.json")
    for row in b1:
        V.require(
            row["old_validator"]["exit"] == 0
            and all(x["disposition"] == "accept" for x in row["old_exact"]),
            "B1-regression",
            row["mutation"],
        )
        V.require(
            all(
                x["exit"] == 2 and x["verdict"]["disposition"] == "reject"
                for x in row["new"]
            ),
            "B1-regression",
            row["mutation"],
        )
    probes = V.load(OUT / "replay-final/parser.json")
    aggregate = V.load(OUT / "replay-final/aggregate-predicate.json")
    V.require(
        all(r["matched"] for r in probes + aggregate), "parser-aggregate", "unmatched"
    )

    registry = tomllib.loads((V.WORKSPACE / "docs/developers.toml").read_text())
    notes = [
        n
        for section in registry["section"]
        for n in section.get("note", [])
        if n["id"] == "generation-oracle-family-diagnosis"
    ]
    V.require(len(notes) == 1, "canonical-doc", "registry identity")
    note = notes[0]
    text = (V.WORKSPACE / note["source"]).read_text()
    V.require(
        [line for line in text.splitlines() if line.startswith("= ")]
        == ["= " + note["title"]],
        "canonical-doc",
        "title",
    )
    V.require("table.header" in text, "canonical-doc", "semantic table header")

    summary = {
        "schema": 2,
        "verdict": "CANDIDATE_ONLY_ROOT_TO_COMMISSION_REVIEW",
        "date": "2026-10-30",
        "round1": {
            "bookmark": "generation-oracle-family-round1",
            "revision": "e814b6b4c3463d9aa1c018939365b624619af24b",
            "review_sha256": "0ec5a10aeae24bea52792f9d48700d81fb05a42f233bfe8e5e101eec300bf008",
        },
        "historical": {
            "diagnostic_ref": "generation-oracle-family-historical-diagnostic",
            "diagnostic_revision": "5e61ef89f4df4e818884e40b8d0c46108eab75a8",
            "pristine_source": "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
            "workspace_left_pristine": True,
        },
        "candidate_ref": "generation-oracle-family-round2",
        "runtime_current": "b9f788922341e7bde2e801a054423ae70ed7e9ae",
        "scope": {
            "cases": list(V.CASES),
            "modes": list(V.MODES),
            "pool": 1,
            "grouping": "None",
            "kernels": [
                "c78cffc20f318df1fce50539a587cbe8ffc717c0",
                "1a1829b1f3d0a88b53b44085b3b9f34ef77945b3",
            ],
        },
        "families": families,
        "validator_sha256": source,
        "semantic_controls": {
            "total": len(controls),
            "matched": len(controls),
            "by_stage": dict(
                collections.Counter(r["observed"].get("stage") for r in controls)
            ),
            "parser_probes": len(probes),
            "aggregate_predicate_controls": len(aggregate),
            "B1_actual_mutants": len(b1),
        },
        "definition_boundary": "Complete parsed whole-symbol tokens after only the five audited mappings; per-expression, per-route, cross-route and emitted-file consistency. Symbolica-normalized equality under identical attribute-bearing definitions, fresh processes per order/kernel. CT/YM nonliteral equalities depend on even cyclic f rotations.",
        "derived": [
            "complete payload/typed-owner domains",
            "fixed owner/path and positional tuple schema",
            "live decoded/payload links",
            "current raw/live link",
            "explicit named-index port bijections",
            "singleton None group chains",
            "factor-preserving intra-route aggregate product under both kernels",
            "CT incident rule identities/couplings/degree one",
            "CT None/Auto occurrence and group-index permutation",
        ],
        "capture_bound": [
            "historical raw cooked expression and producer-decoding association anchored to actual frozen capture; decoded value additionally semantically anchored",
            "graph incidence/routing/physical/model/request/factor values qualified by complete frozen records, not a general physics derivation",
            "no positional downstream consumer theorem",
        ],
        "N1_N2": "No new diagnostic or pristine producer ELF compiled. Unchanged preserved round-one ELFs are replayed. Pinned development-shell preparation may realize ordinary documentation helpers with both license variables absent. Original final-build-log gaps remain explicitly attributed to reviewer reruns and source/compiler/kernel corroboration. attempt1 was a successful pre-decoding adapter; the later license-aborted enhanced attempt lacks a retained log.",
        "N3": "Executable owner/payload domains, raw/live/payload linkage, singleton groups and aggregate product checks added; mirrored aggregate corruption rejects through the real checker, with direct aggregate-predicate controls retained.",
        "N4": "Supplemental closure manifests cover complete relied-on current/historical evidence, review and preflight qualifications. Old manifests remain untouched.",
        "N5": "Both routes derive the same None->Auto position permutation: 109->113,110->109,111->114,112->110,113->115,114->111,115->116,116->112; full keys/IDs and corresponding singleton group/member records are retained.",
        "N6": "Both Rust diagnostics archived byte-for-byte outside ordinary Cargo targets. Reproducible overlay uses the frozen round-one/historical source paths; no production lint suppression or census edits.",
        "N7": "New Python commands use -I -B. Own unnecessary tracked-evidence bytecode removed; frozen target bytecode covered by prior manifests retained as historical evidence.",
        "N8": "CT occurrence records derived independently for both routes/modes: 59 V_36_CT and 209 V_135_CT, each graph exactly one CT-rule vertex, actual declared/used couplings with CT degree exactly one.",
        "N9": "Port equivalence is for named-index contraction, not positional downstream consumers.",
        "N10": "Pinned Typst report compiled, active canonical title/registry checked; direct existing docs-builder replay retains the same eight inherited mark-geometry HTML failures.",
        "no_new_producer_run": True,
        "no_new_diagnostic_build": True,
        "no_new_independent_review_launched": True,
        "syntax_and_controls": "Full-path controls use fresh staged inputs; setup failures cannot be counted as rejection. String-level parser mutations and actual kernel parses are distinct from the unchanged ELF's post-parse Atom corruptions. Literal parsed-tree identities do not require the ELF's power control (1^2=1).",
        "development_attempts": [
            "initial all-owner ELF run rejected a unit owner at its post-parse power control; literal AST identities now remain literal and only nonliteral pairs enter that checker",
            "initial strict lexer omitted Unicode epsilon names and was extended to consume actual Unicode identifiers",
            "first controls command had a syntax error before any mutant ran; preserved transcript and later attempt directory are not counted",
            "pre-final-source baseline and controls are retained as development outputs; credited baseline and controls record the exact final validator SHA",
        ],
        "validation": {
            "canonical_registry": note,
            "pdf_sha256": sha(OUT / "report.pdf"),
            "frozen_evidence_files": len(
                V.load(OUT / "relied-on-frozen-evidence.json")
            ),
            "source_invariance_records": len(V.load(OUT / "source-invariance.json")),
        },
    }
    V.write(HERE / "receipt.json", summary)
    V.write(OUT / "summary.json", summary)
    # Evidence identity manifests are supplementary. Exclude their own entries.
    entries = []
    for root in (
        HERE,
        OUT,
        V.WORKSPACE / "target/oracle-numerator-family-round2",
        V.WORKSPACE / "target/oracle-numerator-family-round2-final",
    ):
        for path in sorted(root.rglob("*")):
            if path.is_file() and path.name not in {
                "round2-manifest.json",
                "closure.json",
            }:
                entries.append(
                    {
                        "path": str(path.absolute()),
                        "sha256": sha(path),
                        "bytes": path.stat().st_size,
                        "symlink": str(path.readlink()) if path.is_symlink() else None,
                    }
                )
    for path in (V.WORKSPACE / note["source"], V.WORKSPACE / "docs/developers.toml"):
        entries.append(
            {
                "path": str(path),
                "sha256": sha(path),
                "bytes": path.stat().st_size,
                "symlink": None,
            }
        )
    V.write(OUT / "round2-manifest.json", entries)
    V.write(
        HERE / "closure.json",
        {
            "frozen_manifest": str(OUT / "relied-on-frozen-evidence.json"),
            "frozen_manifest_sha256": sha(OUT / "relied-on-frozen-evidence.json"),
            "round2_manifest": str(OUT / "round2-manifest.json"),
            "round2_manifest_sha256": sha(OUT / "round2-manifest.json"),
            "manifest_entries": len(entries),
        },
    )
    print(
        json.dumps(
            {
                "families": len(families),
                "controls": len(controls),
                "parser_probes": len(probes),
                "closure_entries": len(entries),
                "receipt_sha256": sha(HERE / "receipt.json"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
