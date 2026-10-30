"""Derive the receipt from final executions, never from expected totals."""

import collections

import common as C


def main():
    layers = {}
    for name in (
        "baseline",
        "legacy",
        "reviewer",
        "domain",
        "staged",
        "anchor",
        "infrastructure",
    ):
        path = C.OUT / "candidate" / name / "results.json"
        rows = C.load(path)
        C.require(rows and all(r["matched"] for r in rows), "finalize", name)
        for row in rows:
            execution = C.OUT / "candidate" / name
            C.require(
                row["command"]["source_hashes"] == C.source_hashes(),
                "finalize",
                "source drift " + str(execution),
            )
        layers[name] = {
            "path": str(path),
            "sha256": C.sha(path),
            "runs": len(rows),
            "dispositions": dict(
                collections.Counter(r["observed"]["disposition"] for r in rows)
            ),
            "stages": dict(collections.Counter(r["observed"]["stage"] for r in rows)),
        }
    domains = []
    for case in C.B.CASES:
        modes = []
        for mode in C.B.MODES:
            folder = C.OUT / "candidate/baseline" / case / mode / "result"
            summary = C.load(folder / "summary.json")
            pairs = C.load(folder / "cross-pairs.json")
            nonliteral = [p for p in pairs if not p["token_tree_equal"]]
            C.require(
                all(p["even_f_literal"] for p in nonliteral),
                "finalize",
                "nonliteral explanation",
            )
            execution = C.load(folder / "execution.json")
            C.require(
                execution["inputs_before"] == execution["inputs_after"]
                and execution["elf_hashes_before"] == execution["elf_hashes_after"],
                "finalize",
                "binding drift",
            )
            modes.append(
                {
                    **summary,
                    "status": "candidate; final baseline passes; independent review pending",
                    "nonliteral_all_even_f_cyclic": True,
                    "elapsed_seconds": execution["elapsed_seconds"],
                    "execution_sha256": C.sha(folder / "execution.json"),
                }
            )
        domains.append({"case": case, "modes": modes})
    dispositions = C.load(C.HERE / "dispositions.json")
    receipt = {
        "schema": 3,
        "title": "Four ungrouped oracle numerator certificates: cross-route alignment and native faithfulness",
        "scope": {
            "cases": list(C.B.CASES),
            "classes": list(C.B.MODES),
            "grouping": "None",
            "pool": 1,
        },
        "theorem": "Symbolica-normalized local/global equality under identical attribute-bearing definitions, after whole-symbol audited maps and exact owner/physical alignment, for these frozen captured domains only",
        "layers": layers,
        "domains": domains,
        "finding_dispositions": dispositions,
        "source_hashes": C.source_hashes(),
        "elf_bindings": {
            route: {
                "path": str(elf),
                "sha256": C.sha(elf),
                "kernel": "c78cffc20f318df1fce50539a587cbe8ffc717c0"
                if route == "historical"
                else "1a1829b1f3d0a88b53b44085b3b9f34ef77945b3",
                "archived_source": str(
                    C.HERE.parent / "round2/archived-sources/oracle_family_exact.rs"
                ),
                "archived_source_sha256": C.sha(
                    C.HERE.parent / "round2/archived-sources/oracle_family_exact.rs"
                ),
                "qualification": "unchanged retained ELF; A independent source/compiler/fingerprint qualification, not new build logs",
            }
            for route, elf in C.B.ELFS.items()
        },
        "frozen_refs": {
            "round1": "e814b6b4c3463d9aa1c018939365b624619af24b",
            "round2": "78c0754c8fb732878123c02999ced3fc506c5a1d",
            "historical_diagnostic": "5e61ef89f4df4e818884e40b8d0c46108eab75a8",
            "historical_pristine_source": "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
            "current_runtime": "b9f788922341e7bde2e801a054423ae70ed7e9ae",
        },
        "limits": [
            "historical raw cook/decoded association capture-bound, not universally derived",
            "named-axis transport not positional downstream-consumer correctness",
            "baseline anchors are self-comparisons; no extra algebra credit",
            "no source/physics universality, other rows, grouped cases, pool8/32 or migration",
            "original final ELF build logs absent; independent A provenance evidence attributed separately",
            "full docs remains blocked at inherited first gate; inherited source46 classification remains masked",
        ],
        "execution_history": {
            "round3_superseded_development_roots": [
                str(C.OUT / "attempts"),
                str(C.OUT / "final"),
            ],
            "only_credited_round3_root": str(C.OUT / "candidate"),
            "development_attempts": str(C.OUT / "attempts"),
            "old_unrecorded_commands": "No retroactive transcript, license-policy or pre-run hash credit. A independent post-run verification is separate evidence.",
            "superseded_development_roots": [
                "target/oracle-numerator-family-round2",
                "target/oracle-numerator-family-round2-final",
            ],
            "credited_round2": "target/oracle-numerator-family-round2-candidate",
            "build_final_logs": "Original preserved build-final.log files describe pre-decoding builds, not credited final binaries; attempt1 runs exited0; later license-aborted attempt has no retained log.",
        },
    }
    C.write(C.OUT / "receipt.json", receipt)
    C.write(C.HERE / "receipt.json", receipt)
    print({name: layer["runs"] for name, layer in layers.items()})


if __name__ == "__main__":
    main()
