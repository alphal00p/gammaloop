"""Derived receipt for this case only; prior approvals are design references."""

import collections

import boundary as B
import common as M


def main():
    binding = M.load(M.OUT / "binding/runtime.json")
    qualification = M.load(M.OUT / "qualification/source-domain.json")
    baselines = []
    for mode in ("None", "Auto"):
        root = M.OUT / "proof" / mode
        verdict = M.load(root / "verdict.json")
        M.require(
            verdict["disposition"] == "accept" and verdict["stage"] == "complete",
            "receipt",
            mode,
        )
        pairs = M.load(root / "cross-pairs.json")
        M.require(
            all(p["token_tree_equal"] for p in pairs),
            "nonliteral-obstruction",
            "new unexplained owner relation",
        )
        baseline = M.load(root / "summary.json")
        baseline.update(
            verdict=verdict,
            execution_sha256=M.sha(root / "execution.json"),
            pair_sha256=M.sha(root / "cross-pairs.json"),
            native_faithfulness_commands_sha256=M.sha(
                root / "native-faithfulness/commands.json"
            ),
            native_aggregate_commands_sha256=M.sha(
                root / "native-aggregate/commands.json"
            ),
        )
        baselines.append(baseline)
    controls = []
    for label in ("controls", "supplements"):
        path = M.OUT / label / "results.json"
        rows = M.load(path)
        M.require(
            rows and all(r["matched"] for r in rows),
            "receipt",
            "all intended predicates " + label,
        )
        controls.append(
            {
                "suite": label,
                "runs": len(rows),
                "result_sha256": M.sha(path),
                "dispositions": dict(
                    collections.Counter(r["observed"]["disposition"] for r in rows)
                ),
                "stages": dict(
                    collections.Counter(r["observed"]["stage"] for r in rows)
                ),
            }
        )
    modes = []
    for route in B.ROOTS:
        path = M.OUT / "binding" / (route + "-modes.json")
        m = M.load(path)
        modes.append(
            {
                "route": route,
                "path": str(path),
                "sha256": M.sha(path),
                "byte_equal": m["byte_equal"],
                "moved": [p for p in m["mapping"] if p["none"] != p["auto"]],
            }
        )
    reused = M.C.source_hashes()
    expected = {
        str(
            B.ROUND3 / "validate.py"
        ): "25805ad1b7b24ae7b68e551dbe97c6cb962f2a9cf733c5ba9d723231a05e9a67",
        str(
            B.ROUND3 / "common.py"
        ): "c4d3c0b1094b21dae30519b3b19881ea053149332c9002d2a296e0116a8694c9",
        str(
            B.ROUND3.parent / "round2/validate.py"
        ): "b468b3a4e394424c8b8add875d3e682ac80f9be9895c82bcc232fee7ac5e9fcb",
    }
    M.require(
        all(M.sha(p) == h for p, h in expected.items()),
        "approved-mechanisms",
        "actual source binding",
    )
    current = M.WORKSPACE / M.MAIN
    donor = M.OUT / "source/current" / M.MAIN
    M.require(M.sha(current) == M.sha(donor), "expectations", "no changes after3432")
    registry = (M.WORKSPACE / "docs/developers.toml").read_text()
    old_registry = (M.OUT / "source/current/docs/developers.toml").read_text()
    start = registry.index('[[section.note]]\nid = "generation-xs-qed-l2"')
    end = registry.index("[[section.note]]", start + len("[[section.note]]"))
    M.require(
        registry[:start] + registry[end:] == old_registry,
        "registry",
        "one insertion; all prior classifications/scopes byte preserved",
    )
    schema = M.load(M.HERE / "schema.json")
    runtime_rows = []
    for route in M.REVISIONS:
        records = M.load(M.OUT / "runtime" / route / "results.json")
        runtime_rows.append(
            {
                "route": route,
                "commands": records,
                "results_sha256": M.sha(M.OUT / "runtime" / route / "results.json"),
            }
        )
    native = M.load(M.OUT / "native-parser-controls/results.json")
    mapping_counts = collections.Counter()
    raw_differences = []
    for mode in ("None", "Auto"):
        ht = M.load(
            B.ROOTS["historical"] / M.CASES[0] / f"classes-{mode}-t1/typed-owners.json"
        )
        ct = M.load(
            B.ROOTS["current"] / M.CASES[0] / f"classes-{mode}-t1/typed-owners.json"
        )
        cm = {g["id"]: {o["owner"]: o for o in g["owners"]} for g in ct}
        for graph in ht:
            for owner in graph["owners"]:
                other = cm[graph["id"]][owner["owner"]]
                if owner["expression"] != other["expression"]:
                    raw_differences.append(
                        {"mode": mode, "id": graph["id"], "owner": owner["owner"]}
                    )
                for text in [
                    owner["producer_decoded_expression"],
                    *owner["ordered_ports"],
                ]:
                    for token in M.C.B.TOKEN.finditer(text):
                        if token.group() in M.C.B.SPELLINGS:
                            mapping_counts[token.group()] += 1
    report = {
        "schema": 1,
        "title": "Ungrouped two-loop QED cross-section numerator certificate",
        "status": "CERTIFICATE_CANDIDATE_REQUIRES_NEW_INDEPENDENT_REVIEW",
        "scope": {
            "case": M.CASES[0],
            "classes": ["None", "Auto"],
            "grouping": "None",
            "pool": 1,
            "excluded": "all other rows, grouped modes, pools8/32, migration/expectation/physics changes, MAIN integration",
        },
        "sources": M.REVISIONS,
        "baselines": baselines,
        "runtime_binding": binding,
        "source_domain_qualification": qualification,
        "schema_contract": {
            "sha256": M.sha(M.HERE / "schema.json"),
            "type_paths": len(schema["types"]),
            "object_paths": len(schema["object_keys"]),
            "derived_files": len(schema["derivation_inputs"]),
            "limit": "fixed actual supported shape; independent of mutable anchors, not universal JSON schema",
        },
        "controls": controls,
        "native_layer_controls": {
            "runs": len(native),
            "sha256": M.sha(M.OUT / "native-parser-controls/results.json"),
            "limit": "direct native parser; full validator rejects trailing syntax earlier, no double pass credit",
        },
        "modes": modes,
        "runtime_execution": runtime_rows,
        "reused_design_sources": reused,
        "required_design_pins": expected,
        "exact_ELFs": {
            r: {"path": str(p), "sha256": M.sha(p)} for r, p in M.C.B.ELFS.items()
        },
        "exact_ELF_source": {
            "path": str(
                M.HERE.parent
                / "generation-oracle-family-diagnosis/round2/archived-sources/oracle_family_exact.rs"
            ),
            "sha256": M.sha(
                M.HERE.parent
                / "generation-oracle-family-diagnosis/round2/archived-sources/oracle_family_exact.rs"
            ),
            "qualification": "retained approved exact grammar/equality engines, source/compiler/kernel corroboration in approved reviews; original helper final build-log gap retained, not inherited for NEW pristine/observer builds which have actual logs",
        },
        "theorem": "For each of this case's two modes independently: complete cross-route physical/resolved-input and symmetric owner alignment, then all540 mapped decoded owners token-identical under identical whole-symbol definitions;168 distinct mapped texts/ports natively faithful in both kernels; each producer aggregate equals prefactor*all local factors in both kernels with sums whole. Named repeat-free axis multiset transport only.",
        "mapping_qualification": "only the actually observed graph symbol spellings are mapped by five declared whole-token source-audited definitions; exercised names/counts listed separately. No global string stripping, no generalized schema removal; tagged notation mapping audited, not a native equivalence theorem.",
        "exercised_mapping_counts_both_modes": dict(mapping_counts),
        "raw_tensor_boundary": {
            "different_owner_instances_both_modes": len(raw_differences),
            "first": raw_differences[0] if raw_differences else None,
            "last_equal": "complete resolved inputs and graph/cut/rule/particle/routing/source-factor physical inventory",
            "first_observed_difference": "raw historical reversible-cooked tensor strings versus current native strings; mapped producer-decoded owners literal-equal, no novel algebraic relation",
            "scope": "bounded observed cause of this numerator encoding delta, no general physical correctness theorem or inference for other rows",
        },
        "current_observer_warning": "one unused NumeratorOptions import in the private observer example build; retained actual source/log, no production suppression or ordinary tracked Cargo target introduced",
        "limits": [
            "historical raw reversible cook/decoded association is actual producer capture, not universal decoding theorem",
            "actual interface presence/order and typed ports are capture-bound; stricter within-route anchor guards added; no port/free-index typing, locality or positional downstream-consumer theorem",
            "Run/pool is explicit inspected capture metadata, not folder-name inference or active-worker/performance theorem",
            "physical graph/rule/cut/factor values observed and cross-route checked, not general physics theorem",
            "observer in-memory typed owners joined to assertion producers by complete exact payload/source companions, not assertion in-memory instrumentation",
            "baseline anchor pairs self-comparisons, no extra independent algebra credit",
            "real private license direct runtime children/no network; restricted-mode correctness does not establish license validity; infrastructure never an inequality",
            "the source-level asserted-reference contract shared but new full None/Auto correspondence independently derived; modes not globally byte-equal",
            "no whole-graph expansion/distribution/traces/numerical or parameter-at-one substitute; exact numeric bookkeeping heads are source-defined rational weights only",
        ],
        "failures": [
            "historical observer attempt1 could not compile Case::numerators (absent at2b1); retained exact failed source/build; retry uses exact historical default options; current serialized options independently equal",
            "pre-capture estimate20owners/graph corrected by actual6vertices9edges+3globals=18; actual540; no memory budget inferred",
            "initial apply_patch context mismatch aborted before affected edits; invalid jq exploratory queries had no proof credit",
            "first report validation expected4semanticheaders but the actual note has3tables; observed exacttitle/noactivecontent/3headers, corrected this harness expectation; entire first render attempt retained with zero final report-validation credit",
            "helper formatting and auxiliary additions during preparation retained by versioned source hashes; not an exhaustive interactive edit transcript",
        ],
        "prior_preservation": {
            "before": str(M.OUT / "preservation/before.json"),
            "after": str(M.OUT / "preservation/after.json"),
            "required": True,
            "no_pruning": True,
        },
        "pristine_current_expectations_sha256": M.sha(current),
        "registry_insert_only": True,
        "docs_limits": "inherited8HTML+retained R2B commands.md first gate, masked unknown sunset owners then unregistered source46; no unrelated repair or full docs-green claim",
        "review_context": {
            "migration_approval_path": ".local/review-family-migration-opus-20261030T050603Z/review.json",
            "sha256": "78d9b395b72a7c070daabc217bddcabb613c001c14f362d059bda1306bdc902a",
            "scope": "prior migration only; not this new-case approval",
        },
    }
    M.write(M.HERE / "receipt.json", report)
    M.write(M.OUT / "receipt.json", report)
    print("baselines", len(baselines), "controls", controls)


if __name__ == "__main__":
    main()
