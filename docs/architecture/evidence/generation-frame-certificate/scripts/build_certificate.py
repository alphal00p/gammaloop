import hashlib
import json
from pathlib import Path

root = Path("target/frame-certificate")
evidence = Path("docs/architecture/evidence/generation-frame-certificate")
diagnosis = json.loads((evidence / "diagnosis.json").read_text())
read = lambda name: json.loads((root / name).read_text())
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()

proofs = {
    "effective_model_and_graph_fingerprint_binding": read("model-binding-summary.json"),
    "full_final_graph_key_and_group_join": read("full-key-join-summary.json"),
    "axis_and_frame": read("axis-and-frame-summary.json"),
    "factorized_stage_expressions": read("pair-expression-summary.json"),
    "generated_signed_maps": read("generated-map-summary.json"),
    "signed_affine_routing_and_named_lmb": read("routing-summary.json"),
    "two_term_prefactor_transport": read("prefactor-summary.json"),
    "source_members_and_ratios": read("member-summary.json"),
    "merged_member_denominators": read("member-denominator-summary.json"),
    "cut_sign_and_source_provenance": read("cut-sign-provenance-summary.json"),
    "corruption_controls": read("corruption-summary.json"),
}
artifacts = [
    root / "frame-model.json",
    root / "all-pair-axis-maps.jsonl",
    root / "pair-graphproof-sha256.txt",
    root / "model-hooks-sha256.txt",
    root / "cut-sign-provenance-full.json",
    root / "current-staged-observation.patch",
    root / "historical-staged-observation.patch",
    root / "current-model-hook-source.patch",
    root / "historical-model-hook-source.patch",
    root / "historical-frozen/comparison-qed.json",
    root / "historical-frozen/representative-qed.json",
    root / "historical-frozen/representative-lr.json",
    root / "current-frozen/comparison-qed.json",
    root / "current-frozen/representative-qed.json",
    root / "current-frozen/representative-lr.json",
]
artifacts.extend(sorted((root / "binaries").iterdir()))
artifacts.extend(
    sorted(path for path in (root / "model-hooks-current").rglob("*") if path.is_file())
)
artifacts.extend(
    sorted(
        path for path in (root / "model-hooks-historical").rglob("*") if path.is_file()
    )
)
artifacts.extend(sorted(root.glob("*-summary.json")))
artifacts.extend(sorted(root.glob("*-run.log")))
artifacts.extend(sorted((evidence / "scripts").glob("*.py")))
artifact_hashes = {str(path): sha(path) for path in artifacts}
pipeline = dict(diagnosis["input_pipeline"])
pipeline["last_confirmed_shared_boundary"] = (
    "Complete source, comparison and representative physical graph records, signed frame maps, named final LMBs, local factors, live typed ports, group membership, cuts and source factors."
)
pipeline["first_observed_difference"] = (
    "Historical constructor presentation tags and current explicit logical-interface fields; both are certified semantic representation changes in this bounded case."
)
pipeline["stage_matrix"] = [
    {
        "stage": "model transformation and request",
        "status": "shared",
        "evidence": "Six live full effective-model receipts are literally equal. Each producer kernel recomputes the authoritative BLAKE3 ModelDefinition fingerprint; all 1,572 pair-stage and 512 final graph 32-byte owner fields match it.",
    },
    {
        "stage": "enumeration and assignment",
        "status": "shared",
        "evidence": "Complete 66/66/130 source pair captures and 60/66/130 final full graph-key joins.",
    },
    {
        "stage": "initial named build for per-diagram LMBs",
        "status": "shared",
        "evidence": "Live first-build and second-build named basis receipts retained for both producers; final internal-edge bijection replays each assignment.",
    },
    {
        "stage": "prefactor and named-LMB build",
        "status": "shared value in different constructor schema",
        "evidence": "Complete typed two-term Q grammar, 1,048 provisional signed transports, 1,572 effective graph tensor bindings, and exact within-route affine basis transport.",
    },
    {
        "stage": "comparison and representative re-preparation",
        "status": "shared value",
        "evidence": "Every actual pair captured before grouping, with exact old/current physical records, factors, maps and typed axis bijections.",
    },
    {
        "stage": "grouping and finalization",
        "status": "shared value",
        "evidence": "All 262 source members, 256 groups, six merged-member ratio identities, cuts, signs and denominator mappings checked.",
    },
    {
        "stage": "census serialization",
        "status": "different native-schema numerator digest",
        "evidence": "Only old/current constructor presentation and explicit interface serialization differ; original historical assertions pass and current hashes fail.",
    },
]

certificate = {
    "schema": 2,
    "verdict": "exact_bounded_factor_preserving_equivalence_established_pending_new_independent_Astra_review",
    "scope": "Exactly three original grouping-frame numerator-digest rows. No expectation or production algorithm change; no recertification.",
    "intermediate_diagnosis": {
        "commit": "67f36295d41cdf872d978692c8aff93d3f3b855e",
        "durable_ref": "refs/tags/generation-frame-intermediate-diagnosis",
        "file": str(evidence / "diagnosis.json"),
    },
    "round1_review": {
        "commit": "c5a19031d8a527ca61903728b6465187f3e87928",
        "durable_ref": "refs/tags/generation-frame-certificate-round1",
        "scope_a": "target/astra-certificate-provenance-20261028T190253Z/scope-a-review.json",
        "scope_b": "target/astra-frame-scope-b-review-20261028T2010Z/review.typ",
        "verdicts": ["REQUEST_CHANGES", "REQUEST_CHANGES"],
        "old_verifier_failures": {
            "scope_a": [
                "changed detached model accepted by denominator proof",
                "one-sided final edge numerator accepted by final join",
            ],
            "scope_b": [
                "trailing prefactor coefficient or scalar accepted by old parser",
                "changed contracted dimension accepted by old parser",
                "both-producer effective prefactor mutation accepted by old replay",
                "both-producer D65 routing coefficient mutation accepted by old replay",
                "final aggregate and swapped master changes accepted by old replay",
            ],
        },
        "new_semantic_replay": "19 injections rejected by their owning proof scripts without SHA guards; exact per-case outcomes are in proofs.corruption_controls.",
    },
    "sources": {
        "historical": "8d49e6ecf0d8fb1169e7433f2c0b41048aa92e60",
        "current_frozen_parent": "59b57cead76e576fa9f23e3996f5614de939fe35",
        "historical_kernel": "c78cffc20f318df1fce50539a587cbe8ffc717c0",
        "current_kernel": "1a1829b1f3d0a88b53b44085b3b9f34ef77945b3",
        "historical_original_capture_commit": "3718714c791d40fb746dd9abc6eeb61fe86548c4",
        "historical_original_capture_ref": "refs/tags/generation-frame-historical-original-capture",
        "constructor_source_sha256": {
            "historical_graph_symbols": sha(
                "../generation-frame-historical-diagnostic/crates/feynkit-graph/src/symbols.rs"
            ),
            "current_graph_symbols": sha("crates/feynkit-graph/src/symbols.rs"),
            "historical_tensor_printer": sha(
                "../generation-frame-historical-diagnostic/crates/spenso/src/network/tags.rs"
            ),
            "current_tensor_printer": sha("crates/spenso/src/network/tags.rs"),
        },
        "qualified_five_constructor_bridge": {
            "path": "../generation-integration-plan/target/wp3-drawing-transplant/frame-digest-boundary.json",
            "sha256": sha(
                "../generation-integration-plan/target/wp3-drawing-transplant/frame-digest-boundary.json"
            ),
            "semantic_argument": "The old q/e/v tags are consumed only as presentation labels by the old tensor printer. The current declarations retain the same names, rank-one/tensor/index/component tags and arguments and register q/e/v as presentation notation. K and P have the same bounded contract. No representation, dimension, variance, dummy scope or coefficient is removed.",
        },
        "restored_production_source_controls": {
            "current_staged_sha256": sha(
                "crates/feynkit-generator/src/generation/staged.rs"
            ),
            "historical_staged_sha256": sha(
                "../generation-frame-historical-diagnostic/crates/feynkit-generator/src/generation/staged.rs"
            ),
            "current_staged_matches_frozen_parent_sha256": "6c8ebb4abd42e875742bb1fdfdfa5c6fc0ecff4e9c0c4adff8f49306b4110b21",
            "historical_staged_matches_frozen_source_sha256": "aeac316e293542fa45ee665696eef1d5e5a7826c2dc64763a9714223495ff6c5",
            "current_nixci_enabled": False,
            "historical_nixci_enabled": False,
            "lockfiles_and_uv_changed": False,
        },
    },
    "producers": {
        "historical_unmodified": {
            "binary": diagnosis["binaries"]["historical_original"]["path"],
            "sha256": diagnosis["binaries"]["historical_original"]["sha256"],
            "original_assertions": "3 passed, 0 failed, record mode unset, one thread",
            "log": diagnosis["binaries"]["historical_original"]["log"],
            "profile_and_features": "ci-optim, default + native, original locked kernel",
        },
        "current_unmodified": {
            "binary": diagnosis["binaries"]["current"]["path"],
            "sha256": diagnosis["binaries"]["current"]["sha256"],
            "original_assertions": "0 passed, 3 expected numerator-digest failures, every other census field equal, record mode unset, one thread",
            "log": diagnosis["binaries"]["current"]["log"],
            "profile_and_features": "ci-optim, default + native, original locked kernel",
        },
        "observational_current": {
            "reused_cache_path_at_capture": "/common/build-cache/lcnbr/00/0198baa63d28f1/ci-optim/deps/generation_oracle-513b63cbc76b9164",
            "graphproof_sha256_at_capture": "0bdb82d0dea6cc8f456dbeada912aad96535b238df0e412abf5afb83d4ea3a0f",
            "temporary_source_patch": str(root / "current-staged-observation.patch"),
            "source_restored_to_parent": True,
            "frozen_final_capture_bytes_equal_original": True,
            "binary_retention": "The ordinary shared cache path was reused by final test recompilation. The exact patch, build identity, source control, logs and immutable capture bytes are retained.",
        },
        "observational_historical": {
            "reused_cache_path_at_capture": "/common/build-cache/lcnbr/42/739ef72d4e7b05/ci-optim/deps/generation_oracle-de23d835aff56301",
            "graphproof_sha256_at_capture": "019c405cb24b299027667b58a5b6df956252a9bd9414f05750eab684cb83d516",
            "temporary_source_patch": str(root / "historical-staged-observation.patch"),
            "source_restored_to_parent": True,
            "frozen_final_capture_bytes_equal_original": True,
            "binary_retention": "The ordinary shared cache path was reused by final test recompilation. The exact patch, build identity, source control, logs and immutable capture bytes are retained.",
        },
    },
    "immutable_binaries": {
        path.name: {"path": str(path), "sha256": sha(path)}
        for path in sorted((root / "binaries").iterdir())
    },
    "observational_source_patches": {
        "current_test_only": {
            "path": str(root / "current-model-hook-source.patch"),
            "sha256": sha(root / "current-model-hook-source.patch"),
        },
        "historical_test_only": {
            "path": str(root / "historical-model-hook-source.patch"),
            "sha256": sha(root / "historical-model-hook-source.patch"),
        },
        "production_generation_staged_restored": True,
    },
    "inputs_and_pipeline": pipeline,
    "rows": diagnosis["rows"],
    "complete_pair_capture": {
        "historical": str(root / "pair-captures-historical-graphproof"),
        "current": str(root / "pair-captures-current-graphproof"),
        "files_per_side": 1051,
        "sha256_manifest": str(root / "pair-graphproof-sha256.txt"),
        "counts_original_sources": {
            "comparison_qed": 66,
            "representative_qed": 66,
            "representative_lr": 130,
        },
        "all_old_current_source_comparison_representative_physical_records_equal": True,
    },
    "proofs": proofs,
    "derived_totals": {
        "full_final_graph_key_joins": 256,
        "final_groups": 256,
        "original_source_members_per_side": 262,
        "extra_merged_members_per_side": 6,
        "factorized_expression_fields": 12960,
        "live_local_tensor_owners": 12174,
        "live_typed_axes": 27576,
        "uniquely_permuted_owners": 2852,
        "current_explicit_interfaces_matching_live_logical_order": 2852,
        "provisional_signed_two_term_prefactor_transports_both_producers": 1048,
        "complete_effective_prefactor_tensor_bindings_both_producers": 1572,
        "within_route_exact_affine_edge_signature_checks_both_producers": 8904,
        "final_graph_master_representative_payload_bindings_both_producers": 512,
        "full_effective_model_fingerprint_kernel_probe_executions": 6,
        "left_right_external_swaps_per_side": 40,
        "cross_section_cut_records_per_stage": 144,
        "cross_section_cut_transports_per_producer": 288,
        "semantic_adversarial_rejections_without_sha_guards": 19,
        "zero_local_or_aggregate_factors_in_captured_pairs": 0,
        "removed_after_built_source_boundary": 0,
    },
    "model_and_merged_denominators": {
        "effective_model_definition_per_producer": [
            str(root / f"model-hooks-{side}/comparison-qed/effective-model.json")
            for side in ("historical", "current")
        ],
        "effective_model_sha256": sha(
            root / "model-hooks-current/comparison-qed/effective-model.json"
        ),
        "graph_model_field_contract": "Only a 32-byte fingerprint is serialized per graph; six independently captured full effective models and producer-kernel probes bind those owner bytes to the literal model definition.",
        "d_particle_ids": [7, 8],
        "s_particle_ids": [9, 10],
        "d_and_s_and_antiparticles_mass": "ZERO",
        "d_and_s_propagator_denominator": "(P(idx(1,1)))^2",
        "d_and_s_propagator_numerators": [
            "(𝑖)*PSlash(idx(1,2),idx(1,1))",
            "(-𝑖)*PSlash(idx(1,1),idx(1,2))",
        ],
        "six_merged_pairs_denominator_and_signed_basis_equality": True,
    },
    "typed_constructor_crosscheck": {
        "current_kernel_D65_prefactor": "passed after audited Q/vertex constructor binding; direct cross-tag parse rejected by Symbolica as expected",
        "historical_kernel_D65_prefactor": "passed after inverse audited constructor binding",
        "current_log": str(root / "d65-typed-bridge.log"),
        "historical_log": "/shared/lcnbr/dev/gammaloop/generation-frame-historical-diagnostic/target/frame-certificate-capture/d65-typed-old-kernel.log",
        "coverage_limit": "Only the D65 typed admission helper executed in both kernels. The all-graph factor and transport proofs are exact structural replay of qualified producer captures, not script executions inside both kernels.",
    },
    "final_validation": {
        "current_final_test_binary_sha256": "9a2297a7794aedb674a6480716035f443ddde72c4fac29be9ed883193c269d27",
        "historical_final_test_binary_sha256": "1eaba7417f8a53a12fd1d2db12f9e8416dbb62c591de12d3a66bf40220e99fc7",
        "current_original_rows": "0 passed, 3 expected numerator-digest failures, 1 diagnostic ignored",
        "historical_original_rows": "3 passed, 0 failed, 1 diagnostic ignored",
        "current_D65_typed": "1 passed",
        "historical_D65_typed": "1 passed",
        "logs": [
            str(root / "final-current-rows.log"),
            str(root / "final-current-d65.log"),
            "/shared/lcnbr/dev/gammaloop/generation-frame-historical-diagnostic/target/frame-certificate-capture/final-historical-rows.log",
            "/shared/lcnbr/dev/gammaloop/generation-frame-historical-diagnostic/target/frame-certificate-capture/final-historical-d65.log",
        ],
        "reused_source39_Astra_review_sha256": "dc4093c671d7fd08d3427da8d34fc8e178304e731cd5f58fd362ce81c4363672",
        "frozen_pair_manifest": "2102/2102 SHA-256 entries verified after cleanup",
        "frozen_pair_manifest_check_log": str(root / "pair-manifest-check.log"),
        "new_certificate_Astra_review": "required and not yet performed",
        "documentation_registry": "The builder accepted the investigation entry and continued to catalog export. End-to-end docs-check is unavailable here because the unrelated gammalooprs catalog exporter executable is absent; no full docs pass is claimed.",
    },
    "factor_preservation_audit": {
        "proof_actions": [
            "read frozen JSON and live metadata",
            "bind full effective models to 32-byte graph fingerprints using each producer kernel",
            "decode historical reversible index labels with historical CookSettings::uncook",
            "bind only audited constructor presentation names",
            "compare exact local and aggregate factorized strings through the actual final master",
            "check complete finite two-term Q prefactor maps and actual effective tensors",
            "replay exact signed affine routing under captured basis coordinates",
            "check graph maps and member ratios",
        ],
        "excluded_actions": [
            "whole-graph numerator expansion",
            "Dirac trace or algebra preparation",
            "distribution of original mass-plus-slash factors",
            "interface omission",
            "arbitrary coefficient normalization",
            "numerical-point inference",
        ],
        "aggregate_export_owner": "FeynmanDiagram::numerator derives a product with SymbolicTensor::multiply_many; named_product_expression uses Atom::mul_many and compact rank-one dot contraction. NumeratorInterface::export may re-admit that aggregate to decide optional interface serialization; admission-only DotNormalizer::nested_vectors, tensor-power lowering and normalize_products were inspected and do not distribute an original sum under graph multiplication or call Dirac algebra. The proof does not additionally prepare or expand the aggregate.",
        "audited_source_owners": [
            "crates/feynkit-graph/src/lib.rs:718",
            "crates/feynkit-graph/src/lib.rs:727",
            "crates/feynkit-graph/src/numerator_interface.rs:157",
            "crates/idenso/src/tensor/inference.rs:799",
            "crates/idenso/src/tensor/composition.rs:190",
            "crates/idenso/src/tensor/composition/product.rs:311",
            "crates/idenso/src/shorthands/schoonschip/normalize_dots.rs:102",
        ],
        "typed_zero_history": "No original local or aggregate factor in the 262 source/comparison/representative pairs is literal zero; no hidden-zero equivalence is inferred.",
    },
    "artifact_sha256": artifact_hashes,
    "limits": [
        "This is exact equivalence for the three frozen producer rows and their original members, not a general grouping or physics theorem.",
        "It certifies labelled tensor value under explicit axis bijections, not identity of private proof caches, archived intermediate histories, or future collision-check state.",
        "No performance or distribution result is claimed.",
        "Historical expected hashes and all production algorithms remain unchanged.",
        "New independent Astra review is required before any targeted native-schema digest recertification; combined-tree validation follows that review.",
    ],
}
(evidence / "certificate.json").write_text(
    json.dumps(certificate, ensure_ascii=False, indent=2) + "\n"
)
print("wrote", evidence / "certificate.json")
