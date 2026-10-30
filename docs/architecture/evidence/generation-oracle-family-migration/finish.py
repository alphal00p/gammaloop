"""Derive a migration receipt from retained executions, never assumed totals."""

import collections
import shutil

import common as M


def main():
    validation = M.load(M.OUT / "validation.json")
    audit = M.load(M.OUT / "audit.json")
    matrix = {}
    for route in ("old", "new"):
        path = M.OUT / "runtime-final" / route / "results.json"
        rows = M.load(path)
        M.require(
            {(r["case"], r["pool"]) for r in rows}
            == {(c, p) for c in M.CASES for p in M.POOLS},
            "receipt",
            "exact case/pool matrix",
        )
        for r in rows:
            M.require(r["selected"] == 1 and r["ignored"] == 0, "selection", r["case"])
            policy = r["command"]["environment"]
            M.require(
                all(policy["update_force_flags_absent"].values())
                and set(policy["oracle_flags"])
                == {"FEYNKIT_ORACLE_THREADS", "FEYNKIT_ORACLE_DUMP"},
                "flags",
                r["case"],
            )
            M.require(
                policy["oracle_flags"]["FEYNKIT_ORACLE_THREADS"] == str(r["pool"]),
                "pool",
                r["case"],
            )
        matrix[route] = {
            "path": str(path),
            "sha256": M.sha(path),
            "runs": len(rows),
            "selected": sum(r["selected"] for r in rows),
            "passed": sum(r["passed"] for r in rows),
            "failed": sum(r["failed"] for r in rows),
            "ignored": sum(r["ignored"] for r in rows),
            "per_case_pool": [
                {
                    k: r[k]
                    for k in (
                        "case",
                        "pool",
                        "selected",
                        "passed",
                        "failed",
                        "ignored",
                        "filtered",
                    )
                }
                for r in rows
            ],
        }
    replay = M.load(M.OUT / "replay/results.json")
    proof = []
    for r in replay:
        path = M.OUT / "replay" / r["kind"] / "results/results.json"
        rows = M.load(path)
        M.require(
            all(row["matched"] for row in rows), "receipt", "replay disposition/stage"
        )
        for row in rows:
            M.require(
                row["command"]["source_hashes"] == M.C.source_hashes(),
                "receipt",
                "unmodified approved checker closure",
            )
        proof.append(
            {
                "kind": r["kind"],
                "count": len(rows),
                "dispositions": dict(
                    collections.Counter(x["observed"]["disposition"] for x in rows)
                ),
                "sha256": M.sha(path),
            }
        )
    approvals = []
    for name, expected in (
        ("a", "45482a7440a8dd5bb808d41176b40b9950865c344fb03f07510411e025783bc8"),
        ("b", "6ceec9c3de1180b4901f9bca6548589adb16b4688e97e162177a61d864f0792d"),
    ):
        path = M.HERE / f"approval-{name}.json"
        M.require(M.sha(path) == expected, "approval", name)
        approvals.append(
            {"path": str(path), "sha256": expected, "verdict": M.load(path)["verdict"]}
        )
    context = (
        M.WORKSPACE.parent
        / "generation-integration-plan/target/wp3-family-migration-preflight-30f-87fa-20261030/preflight.json"
    )
    M.require(
        M.sha(context)
        == "3dc3a7f7b3bff69249feec35e5c7f54439b0b03ae7db711777dbb8f6820a22e9",
        "context",
        "root preflight",
    )
    shutil.copy2(context, M.OUT / "main-preflight-context.json")
    old, new, observer = (
        M.load(M.OUT / "provenance" / (p + ".json")) for p in ("old", "new", "observer")
    )
    project = lambda data: sorted(
        (
            r["unit"].rsplit("-", 1)[0],
            r["kind"],
            r["features"],
            r["profile"],
            r["rustc"],
            str(r["rustflags"]),
        )
        for r in data["closure"]
    )
    observer_projection_equal = project(old) == project(observer)
    M.require(observer["root"]["rustc"] == old["root"]["rustc"], "observer", "compiler")
    order = (
        M.WORKSPACE
        / "target/oracle-numerator-family-round3/qualification/mode-order.json"
    )
    qualifications = {
        "interfaces_ports_locality": "R3B-N1/N3: interface presence and typed-port ORDER capture-bound; theorem is named multiset transport, not free-index inference, locality or positional-consumer correctness. Migration strictly compares actual complete observer bytes/fields with approved captures and joins complete assertion payloads. Observer metadata is not in-memory instrumentation of the assertion ELF.",
        "Run_pool": "R3B-N2/A3-N7: capture.json Run debug text is provenance, unconsumed by semantic validator. Migration explicitly checks actual Some(1/8/32), class, Grouping None, prepared/built counts. Only pool1 is historically mathematically certified.",
        "mapping": "R3B-N6: Q/edge/vertex only of five declared audited maps occur; K/P absent. Historical tagged spellings source-audited, not natively proved equivalent.",
        "license": "R3B-N4: invalid/missing license can produce correct restricted-mode output; no license-validity theorem or inequality credit. Real private config used only at direct runtime exec, no network, private cache/tmp/state. Missing-license policy controls are infrastructure, not kernel inequalities.",
        "predicate_coverage": "R3B-N5: baseline anchors are self-comparisons. Fresh committed 531 suite does not permanently cover all predicates; B's independent additional held-outs are preserved/cited, not claimed freshly rerun here.",
        "mode_order_attribution": "A3-N3: old qualification artifact used final qualify.py with pre-final common.py a8610c37, 20261030T022846018000-final-qualify; A independently verified. Migration independently derives full correspondence from all fresh observations, then checks this certified mapping.",
        "inventory_chronology": "A3-N1/N2: development-before-v1 is strictly pre-run with same map; frozen-before overlapped first runs, per-run guards close consumed-input boundary; closure-development-v1 superseded. New global preservation snapshot began after initial runtime/replay starts, not retroactively pre-run. Strict actual producer build and per-run inputs/ELFs guards retained.",
        "original_history": "A3-N5/N6/N10: unrecorded source-edit/format/rename chronology not reconstructed; original B1 mutants content-identical reserializations, not byte-identical. Some migration helper versions evolved during preparation; actual versions/raw commands and pre-format copies retained; not all interactive edits are a transcript.",
        "build_gaps": "Original final producer build logs absent; independent A compiler/source/fingerprint corroboration retained with limits. New qualified build has its actual source-before/after, Cargo/compiler log and retained ELF. Mutable normal-cache old artifact was replaced by ordinary compilation; retained old ELF and old fingerprint snapshot preserved.",
        "inherited_docs": "A3-N4/N9: full builder first stops at eight inherited HTML plus preserved R2B commands.md; masked unknown sunset owners precede unregistered source46. No unrelated fixes. Original registry review_ref singleton and historical prose remain unchanged.",
        "seal_limits": "A3-N8: original seal self-reference exclusions retained. New detached seal/manifests explicitly exclude their own active record, not mathematical evidence.",
        "source_sets": "Root preflight static findings only: archives outside Cargo auto-targets and scoped generator/Idenso/core compile source, but affect broad commonCargoSources and broad build/fmt/test/WASM/docs sets; docs runtime entire docs tree. Archive paths/bytes unchanged; no filter repair, no general NAR/derivation identity claim. Actual combined sources require later qualification.",
    }
    receipt = {
        "schema": 1,
        "title": "Four certified numerator oracle migrations",
        "status": "SEPARATE_PREPARATION_CANDIDATE_REQUIRES_NEW_INDEPENDENT_MIGRATION_REVIEW",
        "certificate": {
            "revision": M.APPROVED,
            "checkpoint": "generation-oracle-family-approved",
            "approvals": approvals,
            "scope": "four separate cases x None/Auto, Grouping None, pool1, both historical/current locked kernels",
        },
        "migration_scope": {
            "cases": list(M.CASES),
            "runtime_pools": list(M.POOLS),
            "runtime_kernel": "1a1829b1f3d0a88b53b44085b3b9f34ef77945b3",
            "other_rows": "all other74 of remaining78 MAIN rows unchanged and uncertified",
            "MAIN": "not imported or executed",
        },
        "rows": M.load(M.OUT / "derivation.json")["rows"],
        "audit": audit,
        "matrix": matrix,
        "runtime_validation": {
            "path": str(M.OUT / "validation.json"),
            "sha256": M.sha(M.OUT / "validation.json"),
            "pairs": len(validation["rows"]),
            "observer_mode_bindings": len(validation["capture_bindings"]),
        },
        "provenance": {
            p: {
                "path": str(M.OUT / "provenance" / (p + ".json")),
                "sha256": M.sha(M.OUT / "provenance" / (p + ".json")),
                "elf_sha256": data["artifact_sha256"],
                "dependency_units": data["closure_units"],
            }
            for p, data in (("old", old), ("new", new), ("observer", observer))
        },
        "feature_qualification": {
            "assertion_old_new": "complete dependency units, root/compiler/profile/features/flags/config compare exactly; not merely root default/native",
            "observer_dependency_projection_equal": observer_projection_equal,
            "observer_root_profile": observer["root"]["profile"],
            "assertion_root_profile": old["root"]["profile"],
            "observer_note": "build vs test root profiles differ; raw complete per-unit records retained; byte-identical emitted/cache ELF binding",
        },
        "proof_replay": proof,
        "mode_order_reference": {"path": str(order), "sha256": M.sha(order)},
        "negative_controls": M.load(M.OUT / "negatives/results.json"),
        "qualifications": qualifications,
        "MAIN_context": {
            "path": str(context),
            "retained_copy": str(M.OUT / "main-preflight-context.json"),
            "sha256": M.sha(context),
            "authority": "read-only30f context, not import approval; fresh future four-line patch and combined review required",
        },
        "preparation_history": {
            "credited_runtime": "runtime-final/",
            "supporting_previous_runtime": "runtime/; final flags hardened with explicit update/force unsets and repeated24 selections",
            "failed": "three stalled Nix clients interrupted before Cargo/runtime; qualified inventory first failed on gitlink; first observer fingerprint lookup failed because normal cache forwards ELF path; no pass credit",
            "retention": "all prior reviews, large pair/control evidence, raw outputs and failed attempts retained; no cache/store deletion",
        },
    }
    M.write(M.HERE / "receipt.json", receipt)
    M.write(M.OUT / "receipt.json", receipt)
    shutil.copy2(M.OUT / "four-literal.patch", M.HERE / "four-literal.patch")
    print(
        "matrix",
        {k: {f: v[f] for f in ("runs", "passed", "failed")} for k, v in matrix.items()},
        "replay",
        sum(r["count"] for r in proof),
    )


if __name__ == "__main__":
    main()
