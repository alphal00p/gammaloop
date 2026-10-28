"""Append physical capture/reuse or engineering evidence without qualifying performance.

This file-only adapter never opens cache chunks or executes a scientific kernel.
The native producer and supervisor own cache authentication and acceptance.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
from pathlib import Path

from merge_five_loop_results import (
    PRIMARY,
    digest,
    evidence,
    hashed_file,
    load_physical_scope,
    physical_completion,
    physical_scope_note,
    read,
    require,
    validate_physical_scope,
)


def seconds(value):
    require(
        type(value) in (int, float) and math.isfinite(value) and value >= 0,
        "missing or invalid observed duration",
    )
    return value


def reference(item, base):
    path = hashed_file(item, base)
    if "bytes" in item:
        require(
            type(item["bytes"]) is int and path.stat().st_size == item["bytes"],
            "evidence length differs",
        )
    return path


def json_record(item, base):
    """Read authenticated bounded metadata without opening its referenced payloads."""
    path = Path(item["path"])
    path = path if path.is_absolute() else base / path
    require(
        path.suffix == ".json"
        and type(item["bytes"]) is int
        and 0 < item["bytes"] == path.stat().st_size <= 1_000_000,
        "evidence must be bounded JSON metadata",
    )
    return {"evidence": item, "raw_receipt": read(reference(item, base))}


def strict_closed_capture(cache, *, physical_execution=None):
    """Validate capture close; physical collectors attest reduction success separately."""
    native = cache.get("report", {})
    counters = (
        "generation",
        "accepted_records",
        "committed_records",
        "pending_records",
        "pending_bytes",
        "capture_skipped_records",
    )
    errors = ("flush_error", "snapshot_error")
    physical_success = (
        isinstance(physical_execution, dict)
        and physical_execution.get("schema")
        == "vakint.final-instrumented-q0-replay-execution.v1"
        and physical_execution.get("status") == "PASS"
        and physical_execution.get("collector_exit_code") == 0
        and physical_execution.get("collector_wait_is_authoritative") is True
        and physical_execution.get("errors")
        == physical_execution.get("signals_sent")
        == []
    )
    close = cache.get("close", {})
    require(
        isinstance(native, dict)
        and isinstance(close, dict)
        and isinstance(cache.get("snapshot"), dict)
        and cache.get("terminal") is True
        and cache.get("flush_complete") is True
        and all(key in cache and cache[key] is None for key in errors)
        and (
            ("reduction_error" in cache and cache["reduction_error"] is None)
            or ("reduction_error" not in cache and physical_success)
        )
        and {"error", "report"} <= close.keys()
        and close["error"] is None
        and close["report"] == native
        and native.get("schema") == "rustred.running-application-cache.v1"
        and native.get("closed") is True
        and "writer_error" in native
        and native["writer_error"] is None
        and "last_error" in native
        and native["last_error"] is None
        and native.get("capture_metadata_exhausted") is False
        and all(type(native.get(key)) is int and native[key] >= 0 for key in counters)
        and native["generation"] > 0
        and native["accepted_records"] == native["committed_records"]
        and native["pending_records"]
        == native["pending_bytes"]
        == native["capture_skipped_records"]
        == 0
        and all(
            cache.get("flush_report", {}).get(key) == native[key] for key in counters
        )
        and cache.get("snapshot", {}).get("schema")
        == "rustred.running-application-cache-snapshot.v1",
        "strict capture lacks a complete clean flush/snapshot/close",
    )
    return native


def completed_run(config_path, *, expected_job_id="q0"):
    projections = {
        "q0": ("zero-momentum", "zero-momentum-metric"),
        "quadratic": ("quadratic", "quadratic-transverse"),
    }
    require(expected_job_id in projections, "unknown physical job identity")
    config_path = Path(config_path).resolve()
    config = read(config_path)
    require(
        config["schema"] == "vakint.inactive-capture-physical-user-unit.v1",
        "unknown physical campaign configuration schema",
    )
    root = Path(config["output_root"])
    require(root.is_absolute(), "output root must be absolute")
    completion_path = root / "completion.json"
    completion = read(completion_path)
    require(
        config["job_id"] == completion["job_id"] == expected_job_id
        and config["performance_eligible"] is False
        and completion["performance_eligible"] is False
        and completion["functional_acceptance"] is True
        and completion["process_zero_acceptance"] is True
        and completion["collector_exit_code"] == 0
        and completion["config_sha256"] == digest(config_path)
        and completion["pure_ibp_seconds"] is None
        and completion["native_child_restarted_or_signalled"] is False,
        "run lacks accepted immutable functional completion",
    )
    require(
        (config["projection"], config["prepared_input"]["manifest_projection"])
        == projections[expected_job_id],
        "physical job/projection identity differs",
    )
    require(
        not completion.get("acceptance_error")
        and not completion.get("authoritative_wait_error"),
        "completion has errors",
    )
    seconds(completion["collector_wall_seconds_through_actual_wait"])
    process = completion["collector"]
    require(
        all(
            type(process.get(key)) is int and process[key] > 0
            for key in ("pid", "start_ticks")
        ),
        "missing actual process identity",
    )
    launch_path = root / "launch.json"
    launch = read(launch_path)
    argv = launch["argv"]
    require(
        launch["job_id"] == expected_job_id
        and launch["config_sha256"] == digest(config_path)
        and all(
            launch["collector"].get(key) == process[key]
            for key in ("pid", "start_ticks")
        )
        and isinstance(argv, list)
        and all(isinstance(arg, str) for arg in argv)
        and argv.count("--projection") == 1
        and not any(arg.startswith("--projection=") for arg in argv),
        "physical launch identity/projection differs",
    )
    projection_index = argv.index("--projection") + 1
    require(
        projection_index < len(argv)
        and argv[projection_index] == projections[expected_job_id][0],
        "physical launch projection differs",
    )

    cohort = config["validated_cohort"]
    binary = reference(cohort["executable"], config_path.parent)
    require(argv[0] == str(binary), "physical launch executable differs")
    artifact_path = reference(cohort["artifact_receipt"], config_path.parent)
    artifact = read(artifact_path)
    for key, value, length in (
        ("source_head", artifact.get("source_head"), 40),
        ("native_tree", artifact.get("elf", {}).get("native_tree"), 40),
        (
            "source_freeze_sha256",
            artifact.get("elf", {}).get("source_freeze_sha256"),
            64,
        ),
    ):
        require(
            isinstance(value, str) and re.fullmatch(f"[0-9a-f]{{{length}}}", value),
            f"missing or malformed producer provenance: {key}",
        )
    freeze_reference = artifact.get("source_freeze")
    require(
        isinstance(freeze_reference, dict)
        and freeze_reference.get("sha256") == artifact["elf"]["source_freeze_sha256"],
        "source-freeze provenance reference differs",
    )
    freeze = read(reference(freeze_reference, artifact_path.parent))
    require(
        freeze.get("native_tree") == artifact["elf"]["native_tree"]
        and freeze.get("head") == artifact["source_head"],
        "source-freeze provenance identity differs",
    )
    require(
        artifact["state"] == "PASS_OPTIMIZED_CAPTURE_PRODUCER_BUILD_AND_TINY_GATE"
        and artifact["elf"]["sha256"] == digest(binary)
        and artifact["source_head"] == cohort["vakint_revision"]
        and artifact["elf"]["profile"]["opt_level"] == "3"
        and artifact["elf"]["profile"]["test"] is False
        and artifact["elf"]["profile"]["debug_assertions"] is False
        and {"opt-level=3", "lto=fat"}.issubset(
            artifact["actual_rustc_codegen"]["codegen_flags"]
        ),
        "producer does not bind the optimized executable",
    )
    require(
        reference(completion["post_executable"], root) == binary
        and reference(completion["post_artifact_receipt"], root) == artifact_path,
        "post-run executable/artifact differs",
    )
    data_path = reference(cohort["data_manifest"], config_path.parent)
    selected = completion["post_selected_data"]
    require(
        reference(selected["manifest"], root) == data_path
        and selected["revision"] == cohort["data_identity"]
        and selected["version"] == cohort["data_version"],
        "data identity differs",
    )
    for item in cohort["dependency_receipts"]:
        reference(item, config_path.parent)

    run = root / "run"
    payload_names = (
        "benchmark.json",
        "result.symbolica",
        "nonzero-conditions.symbolica",
        "unknown-masters.json",
        "running-cache.json",
    )
    for name in payload_names:
        item = completion["terminal_artifacts"][name]
        require(item["exists"] is True, "missing complete output")
        reference({**item, "path": str(run / name)}, root)
    report = read(run / "benchmark.json")
    outputs = completion["outputs"]
    require(
        report["schema"] == "vakint.qcd-ladder.benchmark.v1"
        and report["retained_result_available"] is True
        and report["reduction_completed"] is True
        and not report.get("reduction_error")
        and report["source_commit"] == cohort["vakint_revision"]
        and report["executable_sha256"] == digest(binary)
        and report["data_manifest_sha256"] == digest(data_path)
        and outputs["benchmark_sha256"] == digest(run / "benchmark.json"),
        "complete report/source/data binding differs",
    )
    require(
        report["input_mode"] == "prepared"
        and report["fresh_native_preparation"] is False
        and report["preparation_seconds"] is None
        and report["coefficient_representation"] == "sparse"
        and report["projection"] == config["prepared_input"]["manifest_projection"]
        and report["contributions"] == config["prepared_input"]["contributions"],
        "physical prepared-input identity differs",
    )
    prepared = config["prepared_input"]
    manifest_path = reference(prepared["manifest"], config_path.parent)
    manifest = read(manifest_path)
    require(
        report["prepared_input"]["manifest_sha256"] == digest(manifest_path)
        and report["prepared_input"]["manifest"] == manifest
        and manifest["schema"] == "vakint.qcd-ladder.prepared-input.v1"
        and manifest["projection"] == prepared["manifest_projection"]
        and manifest["input"] == prepared["input"]
        and manifest["nonzero_conditions_file"] == prepared["guards"]
        and manifest["contributions"] == prepared["contributions"],
        "prepared manifest differs",
    )
    for name, key in (
        ("vakint-input.symbolica", "input"),
        ("projector-conditions.symbolica", "guards"),
    ):
        reference({**prepared[key], "path": str(manifest_path.parent / name)}, root)
    require(
        report["input_sha256"] == prepared["input"]["sha256"]
        and report["input_bytes"] == prepared["input"]["bytes"],
        "scalar input differs",
    )
    for key in ("scientific_acceptance", "provenance_addendum"):
        reference(prepared[key], config_path.parent)
    require(
        report["result_sha256"]
        == outputs["result_sha256"]
        == digest(run / "result.symbolica")
        and report["result_bytes"]
        == outputs["result_bytes"]
        == (run / "result.symbolica").stat().st_size,
        "complete result differs",
    )
    guards = report["guards"]
    guard_path = run / "nonzero-conditions.symbolica"
    require(
        guards["file"] == guard_path.name
        and guards["includes_preparation_conditions"] is True
        and guards["sha256"] == outputs["guards_sha256"] == digest(guard_path)
        and guards["count"]
        == outputs["guards_count"]
        == len(guard_path.read_text().splitlines()),
        "complete ordered guards differ",
    )
    closure = report["catalogue_closure"]
    require(
        read(run / "unknown-masters.json") == []
        and closure["closed"] is True
        and closure["unknown_count"] == outputs["unknown_count"] == 0
        and closure["unknown_sha256"]
        == outputs["unknown_sha256"]
        == digest(run / "unknown-masters.json"),
        "declared catalogue is not closed",
    )
    require(
        report["persistent_cache"] is None
        and report["running_cache_attached"] is True
        and report["persistent_cache_attachments_before"]
        == {str(n): False for n in range(1, 6)},
        "unexpected initial static cache",
    )
    loaded = report["data"]
    require(
        loaded["load_prepared_caches"] is False
        and loaded["version"] == cohort["data_version"]
        and {row["loop_count"] for row in loaded["loops"]} == set(range(1, 6))
        and all(
            row["cache_loaded"] is False and row["source"] == "local"
            for row in loaded["loops"]
        ),
        "unexpected official cache loading",
    )
    cache = read(run / "running-cache.json")
    bound = outputs["running_cache"]
    require(
        reference(bound["report"], root) == run / "running-cache.json",
        "cache report differs",
    )
    descriptor = reference(bound["snapshot_descriptor"], root)
    index = reference(bound["snapshot_index"], root)
    snapshot = cache["snapshot"]
    require(
        read(descriptor) == snapshot
        and snapshot["schema"] == "rustred.running-application-cache-snapshot.v1"
        and snapshot["directory"] == config["running_cache_directory"]
        and Path(snapshot["delta"]["manifest_path"]).resolve() == index
        and snapshot["delta"]["manifest_bytes"] == index.stat().st_size,
        "snapshot identity differs",
    )
    fields = bound["native_fields"]
    expected = {
        "flush_complete": cache["flush_complete"],
        "writer_fault": cache["report"]["writer_error"],
        **{
            key: cache["report"][key]
            for key in (
                "pending_records",
                "pending_bytes",
                "accepted_records",
                "committed_records",
            )
        },
        **{key: snapshot["delta"][key] for key in ("generation", "manifest_digest")},
    }
    require(
        fields == expected
        and fields["flush_complete"] is True
        and fields["writer_fault"] is None
        and fields["pending_records"] == fields["pending_bytes"] == 0
        and fields["accepted_records"] == fields["committed_records"]
        and cache["report"]["closed"] is True
        and cache["close"]["error"] is None
        and cache["report"]["last_error"] is None
        and cache["report"]["capture_metadata_exhausted"] is False
        and all(
            cache.get(key) is None
            for key in ("flush_error", "snapshot_error", "receipt_error")
        ),
        "cache did not finish a clean durable checkpoint",
    )
    seconds(report["reduction_seconds"])
    return {
        "config": config,
        "completion": completion,
        "report": report,
        "cache_report": cache,
        "artifact": artifact,
        "evidence": {
            "config": evidence(config_path),
            "completion": evidence(completion_path),
            "launch": evidence(launch_path),
            "artifact": evidence(artifact_path),
        },
    }


def physical_capture(campaign, physical_scope=None):
    if physical_scope is not None:
        physical_scope = validate_physical_scope(physical_scope)
    campaign = Path(campaign).resolve()
    comparison_path = campaign / "seen-reuse-comparison.json"
    comparison = read(comparison_path)
    require(
        comparison["state"] == "PASS_ONE_FRESH_SEEN_REUSE_EXACT_PARITY"
        and all(
            comparison[key] is True
            for key in (
                "family_map_equal",
                "seen_not_heldout",
                "fresh_process",
                "native_snapshot_unchanged",
            )
        )
        and comparison["newly_captured_records"] == 0
        and comparison["pure_ibp_seconds"] is None,
        "seen reuse is not accepted",
    )
    producer = completed_run(campaign / "q0.final.json")
    reuse = completed_run(campaign / "q0.seen-reuse.final.json")
    a, b = producer["report"], reuse["report"]
    completed_examples = [
        case for case, projection in PRIMARY.items() if a["projection"] == projection
    ]
    require(
        len(completed_examples) == 1
        and (
            physical_scope is None
            or completed_examples[0] in physical_scope["required"]
        ),
        "completed physical projection is outside the required scope",
    )
    pc, rc = producer["config"], reuse["config"]
    allowed = {
        "capture_only_producer",
        "required_predecessor",
        "output_root",
        "running_cache_arguments",
        "seen_reuse_contract",
    }
    require(
        {k: v for k, v in pc.items() if k not in allowed}
        == {k: v for k, v in rc.items() if k not in allowed},
        "producer/reuse scientific configuration differs",
    )
    require(
        pc["capture_only_producer"] is True
        and rc["capture_only_producer"] is False
        and pc["running_cache_arguments"] == ["--running-cache-capture-only"]
        and rc["running_cache_arguments"] == []
        and pc["output_root"] != rc["output_root"],
        "producer/reuse mode differs",
    )
    require(
        reference(rc["required_predecessor"], campaign)
        == Path(producer["evidence"]["completion"]["path"]),
        "reuse predecessor differs",
    )
    contract = rc["seen_reuse_contract"]
    require(
        contract["producer_config_sha256"] == producer["evidence"]["config"]["sha256"]
        and contract["producer_completion_sha256"]
        == producer["evidence"]["completion"]["sha256"]
        and contract["fresh_process"] is True
        and contract["same_elf_input_settings_data_affinity"] is True,
        "seen reuse contract differs",
    )
    for label, run in (("producer", producer), ("reuse", reuse)):
        require(
            comparison[label + "_completion_sha256"]
            == run["evidence"]["completion"]["sha256"],
            "comparison completion differs",
        )
    require(
        tuple(
            producer["completion"]["collector"][key] for key in ("pid", "start_ticks")
        )
        != tuple(
            reuse["completion"]["collector"][key] for key in ("pid", "start_ticks")
        ),
        "reuse did not identify a fresh process",
    )
    require(
        a["conventions"] == b["conventions"] and a["families"] == b["families"],
        "settings/family map differs",
    )
    names = ("result.symbolica", "nonzero-conditions.symbolica", "unknown-masters.json")
    require(
        comparison["exact_complete_payload_sha256_equal"] == dict.fromkeys(names, True),
        "comparison lacks complete parity",
    )
    for name in names:
        require(
            producer["completion"]["terminal_artifacts"][name]["sha256"]
            == reuse["completion"]["terminal_artifacts"][name]["sha256"],
            "producer/reuse complete payload differs",
        )
    require(
        a["running_cache_capture_only"] is True
        and a["running_cache_read_backing"] is False
        and a["saved_application_cache"]["read_backing"] is False
        and a["saved_application_cache"]["persistent"]["lookups"] == 0
        and not any(a["persistent_cache_attachments_during_evaluation"].values())
        and b["running_cache_capture_only"] is False
        and b["running_cache_read_backing"] is True
        and b["saved_application_cache"]["read_backing"] is True
        and b["persistent_cache_attachments_during_evaluation"]
        == {str(n): n == 5 for n in range(1, 6)},
        "capture-only baseline or persistent reuse policy differs",
    )
    old, new = producer["cache_report"], reuse["cache_report"]
    require(
        old["report"]["accepted_records"] > 0
        and new["report"]["accepted_records"] == new["report"]["committed_records"] == 0
        and all(
            old["snapshot"]["delta"][key] == new["snapshot"]["delta"][key]
            for key in ("generation", "manifest_digest")
        ),
        "reuse changed cache snapshot",
    )
    persistent = b["saved_application_cache"]["persistent"]
    require(
        persistent == comparison["cached_persistent_statistics"]
        and all(
            type(persistent.get(key)) is int and persistent[key] >= 0
            for key in ("lookups", "route_hits", "apply_hits")
        )
        and persistent["lookups"] > 0
        and persistent["route_hits"] + persistent["apply_hits"] > 0,
        "reuse lacks matching successful persistent cache hits",
    )
    observed = {
        "producer_reduction": a["reduction_seconds"],
        "cached_reduction": b["reduction_seconds"],
        "producer_wall": producer["completion"][
            "collector_wall_seconds_through_actual_wait"
        ],
        "reuse_wall": reuse["completion"]["collector_wall_seconds_through_actual_wait"],
    }
    require(
        all(
            seconds(value) == comparison[key + "_seconds"]
            for key, value in observed.items()
        ),
        "raw comparison duration differs",
    )
    require(
        comparison["cache_preserved"] == pc["running_cache_directory"],
        "comparison cache identity differs",
    )
    result = {
        "schema": "vakint.physical-capture-report-run.v1",
        "status": "functional-seen-physical-reuse",
        "performance_eligible": False,
        "seen_not_heldout": True,
        "sample_count": 1,
        "observed_seconds": observed,
        "pure_ibp_seconds": None,
        "raw_producer": producer,
        "raw_reuse": reuse,
        "raw_comparison": comparison,
        "evidence": evidence(comparison_path),
        "cache": {
            "producer_accepted_records": old["report"]["accepted_records"],
            "newly_captured_records": 0,
            "snapshot": new["snapshot"],
            "cached_persistent_statistics": persistent,
        },
        "scope": "One fresh seen-input functional reuse; capture writes are included in the producer. No held-out or serialized performance aggregate.",
    }
    if physical_scope is not None:
        result.update(physical_completion(completed_examples, physical_scope))
    return result


def exact_scalar_parity(report, verified, proof_record):
    proof = proof_record["raw_receipt"]
    require(
        verified["full58_proof"]["sha256"] == proof_record["evidence"]["sha256"]
        and verified["full58_scalar_coefficient_parity"] is True
        and verified["original_native_baseline_authenticated"] is True
        and verified["Float_numeric_atoms_rejected"] is True
        and proof["schema"] == "vakint.complete-master-product-parity.v1"
        and proof["status"] == "PASS"
        and type(proof["master_count"]) is int
        and proof["master_count"] > 0
        and proof["master_count"]
        == verified["master_identity_union_count"]
        == report["catalogue_closure"]["master_count"]
        and proof["product_union_count"]
        == verified["full_scalar_coefficient_count"]
        == len(proof["full_scalar_proofs"])
        and proof["product_union_count"] > 0
        and len({row["product"] for row in proof["full_scalar_proofs"]})
        == len(proof["full_scalar_proofs"])
        and all(
            row["exact_difference_zero"] is True for row in proof["full_scalar_proofs"]
        )
        and proof["global_together_or_full_amplitude_difference"] is False
        and proof["unit_mass_or_coupling_replacement"] is False
        and proof["numerical_sampling"] is False,
        "coefficient proof differs from the full exact scalar-product comparison",
    )


def engineering_comparison(index_path):
    """Project accepted JSON receipts only; never reopen symbolic/cache payloads."""
    index_path = Path(index_path).resolve()
    index = read(index_path)
    require(
        index["schema"] == "vakint.q0-engineering-evidence-index.v1"
        and set(index["observations"]) == {"before", "after_off", "after_on"},
        "unknown engineering evidence index",
    )
    observations = {}
    for name, references in index["observations"].items():
        records = {
            kind: {
                "evidence": item,
                "raw_receipt": read(reference(item, index_path.parent)),
            }
            for kind, item in references.items()
        }
        report = records["benchmark"]["raw_receipt"]
        execution = records["execution"]["raw_receipt"]
        verified = records["verification"]["raw_receipt"]
        require(
            report["schema"] == "vakint.qcd-ladder.benchmark.v1"
            and execution["schema"]
            == "vakint.final-instrumented-q0-replay-execution.v1"
            and verified["schema"]
            in {
                "vakint.final-instrumented-q0-verification.v1",
                "vakint.coefficient-local-q0-verification.v1",
            }
            and execution["status"] == verified["status"] == "PASS"
            and execution["collector_exit_code"] == 0
            and execution["collector_wait_is_authoritative"] is True
            and not execution["errors"]
            and not execution["signals_sent"]
            and report["reduction_completed"] is True
            and report["retained_result_available"] is True
            and report["performance_claim"] is False
            and report["projection"] == "zero-momentum-metric"
            and report["input_mode"] == "prepared"
            and report["fresh_native_preparation"] is False
            and report["preparation_seconds"] is None
            and report["source_commit"]
            == report["compiled_cohort"]["main_revision"]
            == execution["source_head"]
            == verified["collector_source_head"],
            "engineering observation lacks accepted completion/source binding",
        )
        require(
            report["end_to_end_scope"]
            == "verified prepared-input loading through cold result serialization; original preparation excluded",
            "engineering end-to-end timer scope differs",
        )
        require(
            verified["whole_ordered_guard_Atoms_equal"] is True
            and verified["families_equal"] is True
            and verified["cache_identity_unchanged"] is True
            and verified["guard_count"] == report["guards"]["count"]
            and verified["whole_canonical_input_sha256"] == report["input_sha256"]
            and verified["whole_canonical_result_sha256"] == report["result_sha256"]
            and verified["cache_identity"]
            == report["cache_identity"]
            == report["cache_identity_after"]
            and verified["measured_public_reduction_seconds"]
            == report["reduction_seconds"]
            and verified["pure_ibp_seconds"] is None
            and verified["HEAVY_file_performance_claim"] is False
            and all(
                verified[key] == 0
                for key in (
                    "numerical_evaluations",
                    "master_values_evaluated",
                    "IBP_runs_in_verifier",
                )
            )
            and report["catalogue_closure"]["closed"] is True
            and report["catalogue_closure"]["unknown_count"] == 0,
            "engineering observation lacks complete exact result/guard checks",
        )
        displayed_phases = {
            "saved_rule_application",
            "retained_catalogue_exact_combination",
        }
        if name == "before":
            displayed_phases.add("retained_outer_exact_combination")
        require(
            displayed_phases <= report["stages"].keys()
            and report["stages"].keys() == verified["stages"].keys(),
            "engineering observation lacks complete verified phase inventory",
        )
        for phase, samples in report["stages"].items():
            require(
                isinstance(samples, list)
                and samples
                and verified["stages"][phase]["samples"] == len(samples)
                and verified["stages"][phase]["total_seconds"]
                == sum(seconds(value) for value in samples),
                "engineering phase differs from verified raw samples",
            )
        if name == "before":
            require(
                verified["whole_retained_Atom_equal"] is True
                and verified["whole_native_baseline_equal"] is True
                and verified["whole_raw_result_equal"] is True,
                "historical engineering baseline is not exact",
            )
        else:
            exact_scalar_parity(report, verified, records["coefficient_proof"])
            require(
                report["cache_hit_diagnostic_enabled"] == (name == "after_on"),
                "engineering coefficient proof or diagnostic mode differs",
            )
        observations[name] = {
            "public_reduction_seconds": seconds(report["reduction_seconds"]),
            "prepared_input_to_serialization_seconds": seconds(
                report["end_to_end_seconds"]
            ),
            "supervised_wall_seconds": seconds(execution["elapsed_seconds"]),
            "pure_ibp_seconds": None,
            "records": records,
        }
    before, off, on = (
        observations[name]["records"]["benchmark"]["raw_receipt"]
        for name in ("before", "after_off", "after_on")
    )
    require(
        before["input_sha256"] == off["input_sha256"] == on["input_sha256"]
        and before["guards"] == off["guards"] == on["guards"]
        and before["families"] == off["families"] == on["families"]
        and off["compiled_cohort"] == on["compiled_cohort"]
        and off["cache_identity"] == on["cache_identity"]
        and off["result_sha256"] == on["result_sha256"]
        and off["result_bytes"] == on["result_bytes"],
        "engineering observations use different input/guards or OFF/ON results",
    )
    return {
        "schema": "vakint.q0-engineering-comparison.v1",
        "status": "verified-engineering-observations",
        "performance_eligible": False,
        "seen_not_heldout": True,
        "matched_serialized_repetitions": False,
        "cache_only_speedup_claim": False,
        "observations": observations,
        "observed_before_over_after_off": {
            "public_reduction": before["reduction_seconds"] / off["reduction_seconds"],
            "prepared_input_to_serialization": before["end_to_end_seconds"]
            / off["end_to_end_seconds"],
        },
        "evidence": evidence(index_path),
        "scope": "Single seen-input engineering observations with different source/cache cohorts and concurrent background work; not final three-sample medians or a cache-only comparison. Prepared-input loading through serialization excludes original graph preparation; supervised wall also includes process startup/exit. Pure IBP and terminal-only times remain unmeasured.",
    }


def cache_file_engineering(index_path, base):
    """Append completed FILE/utility metadata; never hash a referenced cache payload."""
    index_path = Path(index_path).resolve()
    index = read(index_path)
    require(
        index["schema"] == "vakint.cache-file-engineering-evidence-index.v1"
        and set(index["records"])
        == {"heavy", "light", "route_budget", "additions_budget"},
        "unknown cache-file engineering index",
    )

    records = {
        name: json_record(item, index_path.parent)
        for name, item in index["records"].items()
    }
    prior = base["engineering_optimization_runs"][index["baseline_engineering_run"]]
    require(
        prior["schema"] == "vakint.q0-engineering-comparison.v1"
        and prior["performance_eligible"] is False,
        "unknown engineering baseline",
    )
    old = prior["observations"]["after_off"]["records"]["benchmark"]["raw_receipt"]
    heavy = records["heavy"]["raw_receipt"]
    report = heavy["full_actual_benchmark_metadata"]
    verified = heavy["full_actual_scientific_verification_metadata"]
    execution = heavy["full_actual_physical_completion_metadata"]
    require(
        heavy["schema"] == "vakint.full-q0-shared-file-engineering-final-handoff.v1"
        and heavy["status"]
        == "PASS_SINGLE_COLD_FULL_Q0_ENGINEERING_AND_FULL_SCIENTIFIC_VERIFICATION"
        and heavy["performance_claim"] is False
        and execution["status"] == verified["status"] == "PASS"
        and execution["collector_exit_code"] == 0
        and execution["collector_wait_is_authoritative"] is True
        and not execution["errors"]
        and not execution["signals_sent"]
        and report["schema"] == "vakint.qcd-ladder.benchmark.v1"
        and report["reduction_completed"] is True
        and report["performance_claim"] is False
        and report["warm_requested"] is False
        and all(
            report[key] == old[key]
            for key in (
                "compiled_cohort",
                "executable_sha256",
                "input_sha256",
                "projection",
                "input_mode",
                "fresh_native_preparation",
                "end_to_end_scope",
                "guards",
                "result_sha256",
                "result_bytes",
                "families",
                "catalogue_closure",
            )
        )
        and verified["collector_source_head"]
        == execution["source_head"]
        == report["source_commit"]
        and verified["full58_scalar_coefficient_parity"] is True
        and verified["whole_ordered_guard_Atoms_equal"] is True
        and verified["families_equal"] is True
        and verified["whole_canonical_result_sha256"] == report["result_sha256"]
        and verified["guard_count"] == report["guards"]["count"]
        and verified["master_identity_union_count"]
        == old["catalogue_closure"]["master_count"]
        and verified["full_scalar_coefficient_count"]
        == prior["observations"]["after_off"]["records"]["coefficient_proof"][
            "raw_receipt"
        ]["product_union_count"]
        and verified["running_capture_durably_closed"] is True
        and report["running_cache_read_delta"] is False
        and heavy["Native_work"]["native"]
        == report["saved_application_cache"]["native"]
        and heavy["Native_work"]["persistent"]
        == report["saved_application_cache"]["persistent"]
        and verified["full58_proof"]["sha256"]
        == prior["observations"]["after_off"]["records"]["coefficient_proof"][
            "evidence"
        ]["sha256"]
        and verified["pure_ibp_seconds"] is None,
        "FILE observation lacks completed full-parity baseline binding",
    )
    require(
        report["stages"].keys() == verified["stages"].keys()
        and {"saved_rule_application", "retained_catalogue_exact_combination"}
        <= report["stages"].keys(),
        "FILE observation lacks verified phase inventory",
    )
    for phase, samples in report["stages"].items():
        require(
            samples
            and verified["stages"][phase]["samples"] == len(samples)
            and sum(seconds(value) for value in samples)
            == verified["stages"][phase]["total_seconds"],
            "FILE phase differs from verified samples",
        )
    clocks = heavy["clock_summary_seconds"]
    require(
        all(
            clocks[summary] == seconds(source[field])
            for summary, source, field in (
                ("public_call", report, "reduction_seconds"),
                ("collector_E2E", report, "end_to_end_seconds"),
                ("data_loading", report, "data_loading_seconds"),
                (
                    "Native_FILE_loading_auth_decode_index",
                    report,
                    "cache_loading_seconds",
                ),
                (
                    "inclusive_controller_elapsed_after_preflight",
                    execution,
                    "elapsed_seconds",
                ),
            )
        )
        and clocks["decoder_only"] is None
        and clocks["pure_IBP"] is None,
        "FILE clock summaries differ from measured timer scopes",
    )
    legacy = heavy["legacy_FILE_contract"]
    boundaries = heavy["clock_boundaries"]
    require(
        legacy["externally_supplied_typed_FileIdentity"] is False
        and legacy["Native_framing_and_semantic_binding_checked_on_load"] is True
        and legacy["actual_before"]
        == legacy["actual_after"]
        == report["cache_identity"]
        == report["cache_identity_after"]
        and boundaries["wholeFILE_SHA_scan_count"] == 2
        and boundaries["SHA_durations_independently_measured"] is False
        and boundaries["controller_minus_collector_E2E_is_not_pure_SHA"] is True,
        "FILE identity or SHA timer contract differs",
    )

    light = records["light"]["raw_receipt"]
    budget = records["route_budget"]["raw_receipt"]
    additions = records["additions_budget"]["raw_receipt"]
    require(
        light["schema"] == "vakint.light-utility-engineering-pilots.v1"
        and light["status"] == "PASS"
        and light["single_math_worker"] is True
        and light["A_B_full_canonical_result_and_all_guards_equal"] is True
        and [pilot["label"] for pilot in light["pilots"]] == ["A", "B"]
        and budget["schema"] == "vakint.physical-route-budget-pilot-completion.v1"
        and additions["schema"]
        == "vakint.physical-route-budget-additions-pilot-completion.v1"
        and additions["final_performance_selection"] is False
        and all(item["final_held_out_used"] is False for item in (budget, additions)),
        "LIGHT evidence is not completed engineering-only selection",
    )
    candidate_a = json_record(budget["accepted_result"], index_path.parent)
    records["light_baseline"] = json_record(
        light["prior_measurement"], index_path.parent
    )
    light_baseline = records["light_baseline"]["raw_receipt"]
    require(
        light_baseline["success"] is True
        and light_baseline["case_id"] == light["case_id"]
        and light_baseline["identity"]["executable_sha256"]
        == light["compiled_executable_sha256"]
        and light_baseline["cache_statistics"]["native"]["rule_applications"]
        == light["baseline"]["new_native_rules"]
        and all(
            light_baseline[key] == light["baseline"][key]
            for key in (
                "evaluation_seconds",
                "saved_application_including_owner_seconds",
                "stages",
            )
        ),
        "LIGHT historical baseline binding differs",
    )
    candidate_b = next(
        item
        for item in additions["accepted_candidates"]
        if item["cache"] == light["pilots"][1]["cache_file"]
    )
    records["candidate_a"] = candidate_a
    records["candidate_b"] = json_record(candidate_b["result"], index_path.parent)
    candidates = [candidate_a["raw_receipt"], records["candidate_b"]["raw_receipt"]]
    for pilot, candidate in zip(light["pilots"], candidates):
        schema, additional = {
            "A": (
                "vakint.physical-route-budget-candidate-result.v1",
                "additional_keys",
            ),
            "B": (
                "vakint.physical-route-budget-additions-candidate-result.v1",
                "additional_apply_keys",
            ),
        }[pilot["label"]]
        count = candidate["physical_route_limit"] + candidate[additional]
        require(
            candidate["schema"] == schema
            and candidate["native_export_accepted"] is True
            and candidate["exact_record_and_ordered_guard_equality"] is True
            and candidate["final_held_out_used"] is False
            and candidate["cache"] == pilot["cache_file"]
            and 0
            < candidate["actual_compressed_bytes"]
            == pilot["cache_file"]["bytes"]
            <= 95_000_000
            and all(
                candidate["native_counts"][key] == count
                for key in (
                    "selected_records",
                    "reimported_records",
                    "complete_records",
                )
            )
            and candidate["native_counts"]["missing_keys"]
            == candidate["native_counts"]["partial_records"]
            == 0,
            "LIGHT native count/byte/reimport contract differs",
        )
        parity = pilot["parity"]
        capture = pilot["closed_running_capture"]
        require(
            all(
                parity[key] is True
                for key in (
                    "full_canonical_value_bytes_equal",
                    "all_guards_bytes_equal",
                    "canonical_master_families_equal",
                    "catalogue_closed",
                    "frozen_input_owner_roots_equal",
                    "fixed_reduction_settings_equal",
                )
            )
            and parity["unknown_count"] == 0
            and capture["report"]["closed"] is True
            and capture["report"]["pending_records"]
            == capture["report"]["pending_bytes"]
            == 0
            and capture["report"]["writer_error"] is None
            and capture["report"]["committed_records"]
            == capture["report"]["accepted_records"]
            == pilot["actual_cache_statistics"]["persistent"]["captured_records"]
            and capture["all_new_completed_nodes_preserved"] is True
            and capture["read_delta"] is False,
            "LIGHT result/guard parity or durable capture differs",
        )
        timing = pilot["timing"]
        require(
            timing["saved_application_including_owner_seconds"]
            == sum(
                seconds(value) for value in timing["stages"]["saved_rule_application"]
            )
            and timing["native_timing_scope"]["pure_ibp_seconds"] is None,
            "LIGHT saved-call samples or timing scope differ",
        )
        for metric in (
            "evaluation_seconds",
            "data_loading_seconds",
            "cache_loading_seconds",
        ):
            seconds(timing[metric])
        for value in candidate["phase_seconds"].values():
            seconds(value)
        require(
            {
                "load_saved_application",
                "encode_selected_completed_cache",
                "write_selected_cache",
                "reimport_selected_cache",
                "compare_selected_native_visitor",
            }
            <= candidate["phase_seconds"].keys(),
            "LIGHT export lacks measured phase inventory",
        )
    return {
        "schema": "vakint.cache-file-engineering-comparison.v1",
        "status": "verified-engineering-observations",
        "performance_eligible": False,
        "matched_serialized_repetitions": False,
        "final_light_selection": False,
        "baseline_engineering_run": index["baseline_engineering_run"],
        "heavy_observed_change_percent": {
            "public_call": 100
            * (report["reduction_seconds"] / old["reduction_seconds"] - 1),
            "collector_E2E": 100
            * (report["end_to_end_seconds"] / old["end_to_end_seconds"] - 1),
        },
        "records": records,
        "evidence": evidence(index_path),
        "scope": "Single completed cold FILE and tuning utility observations, with concurrent physical work and different cache/capture policies. No causal format advantage, serialized median, held-out selection or final utility winner. Existing timer boundaries, raw stages, null pure-IBP/decoder phases and original cohorts remain authoritative.",
    }


def final_heavy_cache(records):
    """Join accepted metadata only; the native helper owns exhaustive file validation."""
    sources = records["heavy_sources"]["raw_receipt"]
    census = records["heavy_census"]["raw_receipt"]
    exported = records["heavy_export"]["raw_receipt"]
    imported = records["heavy_reimport"]["raw_receipt"]
    controller = records["heavy_reimport_controller"]["raw_receipt"]
    terminals = records["heavy_reimport_terminals"]["raw_receipt"]
    artifact = exported["artifact"]
    details = imported["details"]
    expected = {row["id"]: row for row in sources["sources"]}
    visits = {row["source_id"]: row for row in details["source_visits"]}
    seals = {row["source_id"]: row for row in imported["source_seals"]}
    require(
        census["schema"] == "vakint.final-thirteen-source-census-postflight.v1"
        and census["status"] == "PASS_ACTUAL_FULL13_CENSUS_ORIGINAL_PIDFDS_TERMINAL"
        and exported["schema"] == "vakint.final-thirteen-source-export-postflight.v1"
        and exported["status"]
        == "PASS_ACTUAL_ONE_FINAL13_EXPORT_ORIGINAL_PIDFDS_TERMINAL"
        and exported["source_manifest"]["sha256"]
        == records["heavy_sources"]["evidence"]["sha256"]
        and exported["actual_accepted_census"]["sha256"]
        == records["heavy_census"]["evidence"]["sha256"]
        and sources["source_count"]
        == census["source_count"]
        == exported["source_count"]
        and sources["total_stored_record_references"]
        == census["input_all_version_references"]
        == exported["all13_input_references"]
        == details["input_records"]
        and census["unique_records"]
        == exported["actual_unique_records"]
        == exported["actual_writer_statistics"]["records"]
        and census["identical_duplicates"] == exported["identical_duplicates"]
        and census["unique_records"] + census["identical_duplicates"]
        == details["input_records"]
        and census["full_row_conflicts"] == 0
        and census["original12_source_rows_and_seals_unchanged_superset"] is True
        and exported["original12_source_rows_and_index_seals_unchanged_superset"]
        is True
        and exported["all13_source_seals_unchanged_and_equal_actual_census"] is True
        and exported[
            "all_census_record_duplicate_coefficient_payload_and_row_spool_metrics_match"
        ]
        is True,
        "final HEAVY source/census/export bindings differ",
    )
    require(
        imported["schema"] == "vakint.heavy-native-phase.v1"
        and imported["phase"] == "reimport"
        and imported["status"] == "PASS"
        and imported["performance_claim"] is False
        and imported["request_sha256"] == exported["request"]["sha256"]
        and imported["request_sha256"] == controller["request_sha256"]
        and details["all_exact_record_values_and_ordered_guards_equal"] is True
        and details["state_imports_shared"] is True
        and details["cache_path"] == artifact["path"]
        and details["cache_sha256"] == artifact["sha256"]
        and details["cache_identity"] == artifact["typed_Native_identity"]
        and all(
            artifact[key] == artifact["typed_Native_identity"][key]
            for key in ("compressed_bytes", "decoded_bytes")
        )
        and len(visits)
        == len(details["source_visits"])
        == len(seals)
        == len(imported["source_seals"])
        == sources["source_count"]
        and visits.keys() == seals.keys() == expected.keys()
        and all(
            visits[source_id]["emitted_references"] == row["stored_record_references"]
            and all(
                seals[source_id][key] == row[key]
                for key in (
                    "directory",
                    "generation",
                    "native_index_identity",
                    "stored_record_references",
                )
            )
            for source_id, row in expected.items()
        ),
        "final HEAVY lacks exhaustive source-reference/ordered-guard round-trip equality",
    )
    require(
        controller["schema"] == "vakint.heavy-finalization-controller-phase.v1"
        and controller["phase"] == "reimport"
        and controller["status"] == "PASS"
        and controller["child_returncode"] == 0
        and controller["signals_sent"] == []
        and controller["driver_receipt_sha256"]
        == records["heavy_reimport"]["evidence"]["sha256"]
        and controller["helper_sha256"] == imported["helper_executable"]["sha256"]
        and controller["child"]["pid"] == imported["pid"]
        and terminals["schema"]
        == "vakint.final-heavy-phase-original-pidfds-terminal.v1"
        and terminals["phase"] == "reimport"
        and terminals["all_original_PIDFDs_terminal"] is True
        and terminals["signals_sent"] == []
        and all(
            terminals["terminal"][name]["original_PIDFD_terminal"] is True
            and all(
                terminals["terminal"][name]["original"][key] == controller[role][key]
                for key in ("pid", "start_ticks", "uid")
            )
            for name, role in (("helper", "child"), ("controller", "controller"))
        ),
        "final HEAVY reimport lacks an actual successful wait and original terminals",
    )
    return {
        "artifact": artifact,
        "source_count": sources["source_count"],
        "record_references": details["input_records"],
        "unique_records": census["unique_records"],
        "identical_duplicates": census["identical_duplicates"],
        "unique_coefficients": exported["actual_writer_statistics"][
            "unique_coefficients"
        ],
        "native_phase_seconds": {
            "census": seconds(census["Native_seconds"]),
            "export": seconds(exported["Native_seconds"]),
            "reimport": seconds(imported["seconds"]),
        },
        "controller_phase_seconds": {
            "census": seconds(census["controller_elapsed_seconds"]),
            "export": seconds(exported["controller_elapsed_seconds"]),
            "reimport": seconds(controller["elapsed_seconds"]),
        },
        "cohort": {
            "native_revision": imported["native_revision"],
            "symbolica_revision": imported["symbolica_revision"],
            "request_main_revision_label": imported["main_revision"],
            "helper_executable": imported["helper_executable"],
            "qualification": "Native helper does not link Vakint; its request Main revision is provenance, not the final public consumer cohort.",
        },
        "comparison_scope": "Every original source reference, exact coefficient and original ordered guard vector including multiplicity; final file reimport and shared state proved by the existing Native implementation. No cache payload is read by this metadata adapter.",
    }


def final_cached_q0(records, heavy, reference_report):
    """Validate one final-file cached observation without deriving matched timings."""
    report = records["final_q0_benchmark"]["raw_receipt"]
    execution = records["final_q0_execution"]["raw_receipt"]
    verified = records["final_q0_verification"]["raw_receipt"]
    monitor = records["final_q0_monitor"]["raw_receipt"]
    terminals = records["final_q0_terminals"]["raw_receipt"]
    capture = records["final_q0_capture"]["raw_receipt"]
    require(
        report["schema"] == "vakint.qcd-ladder.benchmark.v1"
        and execution["schema"] == "vakint.final-instrumented-q0-replay-execution.v1"
        and verified["schema"]
        == "vakint.shared-file-physical-q0-engineering-verification.v1"
        and execution["status"] == verified["status"] == "PASS"
        and execution["collector_exit_code"] == 0
        and execution["collector_wait_is_authoritative"] is True
        and execution["errors"] == execution["signals_sent"] == []
        and report["reduction_completed"] is True
        and report["retained_result_available"] is True
        and report["performance_claim"] is False
        and report["warm_requested"] is False
        and report["projection"] == "zero-momentum-metric"
        and report["input_mode"] == "prepared"
        and report["fresh_native_preparation"] is False
        and report["preparation_seconds"] is None
        and report["input_sha256"] == reference_report["input_sha256"]
        and report["source_commit"]
        == report["compiled_cohort"]["main_revision"]
        == execution["source_head"]
        == verified["collector_source_head"]
        and report["compiled_cohort"]["native_revision"]
        == heavy["cohort"]["native_revision"]
        and report["compiled_cohort"]["kernel_revision"]
        == heavy["cohort"]["symbolica_revision"]
        and report["end_to_end_scope"]
        == "verified prepared-input loading through cold result serialization; original preparation excluded",
        "final cached q0 lacks a completed single saved-input/source observation",
    )
    exact_scalar_parity(report, verified, records["final_q0_coefficient_proof"])
    require(
        verified["whole_ordered_guard_Atoms_equal"] is True
        and verified["families_equal"] is True
        and verified["cache_identity_unchanged"] is True
        and verified["guard_count"]
        == report["guards"]["count"]
        == reference_report["guards"]["count"]
        and verified["whole_canonical_input_sha256"] == report["input_sha256"]
        and verified["whole_canonical_result_sha256"] == report["result_sha256"]
        and verified["cache_identity"]
        == report["cache_identity"]
        == report["cache_identity_after"]
        and report["cache_identity"]["cache_sha256"] == heavy["artifact"]["sha256"]
        and verified["file_attachment"]["path"] == heavy["artifact"]["path"]
        and verified["file_attachment"]["known_sha256"] == heavy["artifact"]["sha256"]
        and verified["measured_public_reduction_seconds"] == report["reduction_seconds"]
        and verified["pure_ibp_seconds"] is None
        and verified["HEAVY_file_performance_claim"] is False
        and all(
            verified[key] == 0
            for key in (
                "numerical_evaluations",
                "master_values_evaluated",
                "IBP_runs_in_verifier",
            )
        )
        and report["catalogue_closure"]["closed"] is True
        and report["catalogue_closure"]["unknown_count"] == 0
        and verified["running_capture_durably_closed"] is True
        and verified["running_cache_read_delta"] is False
        and report["running_cache_read_delta"] is False
        and capture["read_delta"] is False
        and report["saved_application_cache"]["native"]
        == verified["actual_native_work"]
        and report["saved_application_cache"]["persistent"]
        == verified["actual_persistent_work_and_captures"]
        and strict_closed_capture(capture, physical_execution=execution)[
            "committed_records"
        ]
        == 0
        and capture["path"] == verified["running_cache_path"],
        "final cached q0 lacks full exact science/file identity or a closed preserved capture",
    )
    require(
        report["stages"].keys() == verified["stages"].keys()
        and all(
            isinstance(samples, list)
            and samples
            and verified["stages"][phase]["samples"] == len(samples)
            and verified["stages"][phase]["total_seconds"]
            == sum(seconds(value) for value in samples)
            for phase, samples in report["stages"].items()
        ),
        "final cached q0 stage samples differ from the scientific verifier",
    )
    require(
        monitor["schema"] == "vakint.q0-monitor-responsiveness.v1"
        and monitor["status"] == "PASS"
        and monitor["fault_count"] == 0
        and monitor["performance_claim"] is False
        and monitor["actual_worker_terminal"] is True
        and seconds(monitor["max_conservative_freshness_bound_seconds"]) <= 30
        and all(
            monitor["worker"][key] == execution["worker"][key]
            for key in ("pid", "start_ticks")
        )
        and terminals["schema"] == "vakint.q0-monitor-original-three-pidfds-terminal.v1"
        and terminals["all_three_original_PIDFDs_terminal"] is True
        and terminals["signals_sent"] == []
        and all(
            terminals["terminal"][name]["original_PIDFD_terminal"] is True
            and terminals["terminal"][name]["original"]
            == terminals["registration"]["originals"][name]
            and all(
                terminals["terminal"][name]["original"][key] == original[key]
                for key in ("pid", "start_ticks")
            )
            and (
                "uid" not in original
                or terminals["terminal"][name]["original"]["uid"] == original["uid"]
            )
            for name, original in (
                ("worker", execution["worker"]),
                ("controller", execution["controller"]),
                ("guard", execution["guard_registration"]["guard"]),
            )
        ),
        "final cached q0 lacks a fresh full-job monitor and actual original terminals",
    )
    return {
        "sample_count": 1,
        "seen_input": True,
        "performance_eligible": False,
        "cohort": report["compiled_cohort"],
        "executable_sha256": report["executable_sha256"],
        "file_sha256": heavy["artifact"]["sha256"],
        "master_count": verified["master_identity_union_count"],
        "scalar_coefficient_count": verified["full_scalar_coefficient_count"],
        "ordered_public_guard_count": verified["guard_count"],
        "observed_seconds": {
            "input_loading": seconds(report["input_loading_seconds"]),
            "data_loading": seconds(report["data_loading_seconds"]),
            "cache_loading_authentication_decode_index": seconds(
                report["cache_loading_seconds"]
            ),
            "public_reduction": seconds(report["reduction_seconds"]),
            "prepared_input_to_serialization": seconds(report["end_to_end_seconds"]),
            "controller_inclusive": seconds(
                execution["controller_inclusive_elapsed_seconds"]
            ),
            "saved_native_call_inclusive": sum(
                report["stages"]["saved_rule_application"]
            ),
            "final_coefficient_combination": sum(
                report["stages"]["retained_catalogue_exact_combination"]
            ),
        },
        "pure_ibp_seconds": None,
        "isolated_sha_seconds": None,
        "actual_native_work": verified["actual_native_work"],
        "actual_persistent_work": verified["actual_persistent_work_and_captures"],
        "monitor": monitor,
        "identity_scope": verified["file_attachment"],
        "timing_scope": "One n=1 monitored final-file cached saved-input process. Retain original phase boundaries and both inherent FILE SHA scans; loading/decode/hash/CAS are not isolated IBP. No matched median, causal cache-only speedup or unseen-input result.",
    }


def delivery_checkpoint(index_path, base):
    """Record delivery decisions and functional evidence, never performance medians."""
    index_path = Path(index_path).resolve()
    index = read(index_path)
    required = {"light_freeze", "monitor", "authorization", "heavy_sources"}
    publication_records = {"publication", "release_assets"}
    heavy_records = {
        "heavy_census",
        "heavy_export",
        "heavy_reimport",
        "heavy_reimport_controller",
        "heavy_reimport_terminals",
    }
    q0_records = {
        "final_q0_benchmark",
        "final_q0_execution",
        "final_q0_verification",
        "final_q0_coefficient_proof",
        "final_q0_monitor",
        "final_q0_terminals",
        "final_q0_capture",
    }
    optional = (
        publication_records
        | heavy_records
        | q0_records
        | {"heldout_correctness", "published_consumers"}
    )
    keys = set(index["records"])
    require(
        index["schema"] == "vakint.q0-delivery-evidence-index.v1"
        and required <= keys <= required | optional
        and (not keys & publication_records or publication_records <= keys)
        and (not keys & heavy_records or heavy_records <= keys)
        and (not keys & q0_records or q0_records | heavy_records <= keys),
        "unknown delivery evidence index",
    )
    records = {
        name: json_record(item, index_path.parent)
        for name, item in index["records"].items()
    }
    freeze = records["light_freeze"]["raw_receipt"]
    monitor = records["monitor"]["raw_receipt"]
    authorization = records["authorization"]["raw_receipt"]
    sources = records["heavy_sources"]["raw_receipt"]
    require(
        freeze["schema"] == "vakint.light-cache.user-directed-freeze.v1"
        and freeze["performance_claim"] is False
        and freeze["future_tuning_allowed"] is False
        and freeze["all_heldout_outcomes_used_for_selection"] is False
        and 0 < freeze["compressed_bytes"] <= 95_000_000
        and freeze["selected_records"] > 0,
        "LIGHT lacks an accepted non-performance freeze",
    )
    study = base["engineering_optimization_runs"][index["cache_file_engineering_run"]]
    candidate = study["records"]["candidate_b"]["raw_receipt"]
    require(
        candidate["actual_compressed_bytes"] == freeze["compressed_bytes"]
        and candidate["actual_decoded_bytes"] == freeze["decoded_bytes"]
        and candidate["cache"]["sha256"] == freeze["sha256"]
        and candidate["native_counts"]["complete_records"]
        == freeze["selected_records"],
        "frozen LIGHT differs from the complete native export",
    )
    require(
        authorization["schema"] == "vakint.q0-delivery-authorization.v1"
        and authorization["new_uncached_repetitions_waived"] is True
        and authorization["required_physical_examples"] == ["zero-momentum"]
        and authorization["user_prompt_verbatim"].strip(),
        "historical-baseline waiver or physical scope is absent",
    )
    historical = base["physical_capture_runs"][index["historical_physical_run"]]
    require(
        historical["performance_eligible"] is False,
        "historical reference must retain its original non-performance cohort",
    )
    require(
        monitor["schema"] == "vakint.q0-monitor-responsiveness-actual-final-handoff.v1"
        and monitor["status"]
        == "PASS_ACTUAL_ONE_CACHED_Q0_MONITOR_AND_FULL55_58_146678_SCIENCE"
        and monitor["science"]["all_original_mass_coupling_normalization_preserved"]
        is True
        and monitor["science"]["input_result_guard_and_family_identity_authenticated"]
        is True
        and monitor["monitor"]["fault_count"] == 0
        and monitor["monitor"]["all_actual_worker_births_match"] is True
        and monitor["actual_OS_resources"]["signals_sent"] == [],
        "monitor lacks complete functional acceptance",
    )
    freshness = seconds(
        monitor["monitor"]["maximum_conservative_freshness_bound_seconds"]
    )
    require(
        freshness <= monitor["monitor"]["required_limit_seconds"]
        and sources["schema"] == "vakint.heavy-final-union-source-manifest.v1"
        and sources["status"] == "complete"
        and sources["pending_sources"] == []
        and len(sources["sources"]) == sources["source_count"]
        and sum(row["stored_record_references"] for row in sources["sources"])
        == sources["total_stored_record_references"],
        "monitor freshness or closed HEAVY source inventory differs",
    )
    checkpoint = {
        "schema": "vakint.q0-delivery-checkpoint.v1",
        "performance_eligible": False,
        "cache_performance_winner": False,
        "new_uncached_repetitions_waived": True,
        "historical_physical_run": index["historical_physical_run"],
        "historical_uncached_seconds": seconds(
            historical["observed_seconds"]["producer_reduction"]
        ),
        "cache_file_engineering_run": index["cache_file_engineering_run"],
        "records": records,
        "pending": {
            "final_heavy_unique_records": None,
            "final_heavy_compressed_bytes": None,
            "final_heavy_decoded_bytes": None,
            "final_heavy_sha256": None,
            "final_heavy_export_reimport_seconds": None,
            "published_release_url": None,
            "actual_download_and_offline_reuse": None,
        },
        "evidence": evidence(index_path),
    }
    if publication_records <= keys:
        publication = records["publication"]["raw_receipt"]
        assets = records["release_assets"]["raw_receipt"]
        version = publication["version"]
        require(
            publication["schema"] == "vakint.external-data.serial-publication.v1"
            and publication["status"] == "PASS"
            and [stage["stage"] for stage in publication["stages"]]
            == ["source-push", "immutable-then-latest"]
            and all(stage["exit_code"] == 0 for stage in publication["stages"])
            and publication["immutable_release"]
            == f"https://github.com/alphal00p/vakint-external-data/releases/tag/{version}"
            and publication["latest_alias"]
            == "https://github.com/alphal00p/vakint-external-data/releases/tag/data-latest"
            and assets["schema"] == "vakint.frozen-light-full-release-pack-preflight.v1"
            and assets["status"]
            == "PASS_TWO_INDEPENDENT_FULL_PACKS_AND_LOCAL_PUBLISHER_PREFLIGHT"
            and assets["source_commit"] == publication["source"]
            and assets["checksums_byte_identical"] is True
            and assets["local_publisher_checked_assets_both_passed"] is True
            and assets["all_assets_below_decimal2GB"] is True,
            "publication lacks successful ordered stages or source-bound release assets",
        )
        checkpoint["publication"] = {
            "version": version,
            "source_commit": publication["source"],
            "immutable_release": publication["immutable_release"],
            "latest_alias": publication["latest_alias"],
            "release_manifest": assets["assets"]["manifest.json"],
            "five_loop_pack": assets["assets"][
                f"vakint-external-data-{version}-5L.tar.gz"
            ],
        }
        del checkpoint["pending"]["published_release_url"]
    if "published_consumers" in records:
        require(
            publication_records <= keys, "consumer proof lacks publication bindings"
        )
        consumers = records["published_consumers"]["raw_receipt"]
        published = checkpoint["publication"]
        require(
            consumers["schema"] == "vakint.final-published-B-consumer-handoff.v1"
            and consumers["status"]
            == "PASS_ACTUAL_FRESH_IMMUTABLE_LATEST_OFFLINE_AND_CACHE_OPTOUT"
            and consumers["publication"]["sha256"]
            == records["publication"]["evidence"]["sha256"]
            and consumers["data_source_revision"] == published["source_commit"]
            and consumers["version"] == published["version"]
            and consumers["release_manifest_sha256"]
            == published["release_manifest"]["sha256"]
            and consumers["five_loop_archive"]["sha256"]
            == published["five_loop_pack"]["sha256"]
            and consumers["five_loop_archive"]["bytes"]
            == published["five_loop_pack"]["bytes"]
            and all(
                consumers["selected_B"][key] == freeze[freeze_key]
                for key, freeze_key in (
                    ("records", "selected_records"),
                    ("compressed_bytes", "compressed_bytes"),
                    ("decoded_bytes", "decoded_bytes"),
                    ("sha256", "sha256"),
                )
            )
            and consumers["performance_claim"] is False
            and consumers["five_loop_reduction_evaluate_calls"] == 0
            and consumers["signals_sent"] == []
            and consumers["catalogue"]["actual_get_master_alias_checks"]
            == consumers["catalogue"]["declared_uncomputed_5L"]
            == 640,
            "published consumer lacks authenticated release/cache/getter bindings",
        )
        stages = {row["label"]: row for row in consumers["stages"]}
        expected = {
            "fresh-immutable",
            "fresh-latest",
            "offline-immutable",
            "offline-latest",
            "offline-optout-control",
        }
        require(
            len(consumers["stages"]) == len(expected) and stages.keys() == expected,
            "published consumer stage inventory differs",
        )
        for label, row in stages.items():
            reused = label.startswith("offline-")
            prepared = label != "offline-optout-control"
            origin = "reused" if reused else "downloaded"
            require(
                row["status"] == "PASS"
                and row["original_worker_PIDFD_terminal"] is True
                and row["load_source"] == origin
                and row["loop_sources"] == {str(loop): origin for loop in range(1, 6)}
                and row["load_prepared_caches"] is prepared
                and row["loop_cache_attachments"]
                == {str(loop): prepared for loop in range(1, 6)}
                and (not reused or row["official_API_TCP_connects"] == 0)
                and (
                    row["file_statistics"]["records"] == freeze["selected_records"]
                    if prepared
                    else row["file_statistics"] is None
                ),
                "published download/reuse/cache-optout stage is incomplete",
            )
        checkpoint["published_consumers"] = {
            "fresh_downloads": ["immutable", "latest"],
            "fresh_process_offline_reuse": ["immutable", "latest"],
            "cache_optout_control": True,
            "all_640_get_master_checks": True,
            "network_scope": consumers["network_scope"],
            "qualified_community": consumers["qualified_community"],
        }
        del checkpoint["pending"]["actual_download_and_offline_reuse"]
    if "heldout_correctness" in records:
        heldout = records["heldout_correctness"]["raw_receipt"]
        require(
            heldout["schema"]
            == "vakint.single-heldout-correctness-and-replay-handoff.v1"
            and heldout["status"]
            == "PASS_FULL_VALUE_PUBLIC_GUARD_SET_USED_MASTERS_FAMILIES_CATALOGUE"
            and heldout["performance_claim"] is False
            and heldout["eligibility_or_selection_authority"] is False
            and heldout["full_unexported_result_families_map_claim"] is False
            and heldout["literal_public_vector_order_claim"] is False
            and heldout["LIGHT_B_sha256"] == freeze["sha256"]
            and heldout["all_original_worker_controller_PIDFD_terminal"] is True
            and heldout["strict_noReadDelta_both_closed"] is True
            and all(
                heldout[key] is True
                for key in (
                    "full_value_bytes_equal",
                    "public_canonical_guard_SET_equal",
                    "used_master_and_family_inventory_equal",
                    "catalogue_closed_equal",
                )
            )
            and [row["stage"] for row in heldout["rows"]] == ["reduction", "replay"],
            "held-out case lacks qualified functional completion",
        )
        first, replay = heldout["rows"]
        captures = {}
        for row in heldout["rows"]:
            capture_refs = [
                item
                for item in row["raw_outputs"]
                if Path(item["path"]).name == "running-cache.json"
            ]
            require(
                len(capture_refs) == 1, "held-out stage lacks one saved capture receipt"
            )
            capture = json_record(capture_refs[0], index_path.parent)
            raw_capture = capture["raw_receipt"]
            native = strict_closed_capture(raw_capture)
            require(
                native == row["closed_capture_report"]
                and raw_capture["path"] == row["closed_capture_path"]
                and raw_capture["read_delta"] is False,
                "held-out saved capture differs from its closed stage receipt",
            )
            captures[row["stage"]] = capture
        require(
            all(
                first[key] == replay[key]
                for key in (
                    "input_sha256",
                    "full_value_sha256",
                    "public_guard_set_sha256",
                    "public_guard_set_count",
                    "used_masters",
                    "used_family_sha256",
                    "raw_reduction_settings",
                    "raw_reduction_settings_source",
                    "conventions",
                )
            )
            and all(
                row["catalogue_closure"]["closed"] is True
                and row["catalogue_closure"]["unknown_count"] == 0
                and not row["signals_sent"]
                for row in heldout["rows"]
            )
            and replay["cache_statistics"]["native"]["rule_applications"] == 0
            and replay["cache_statistics"]["persistent"]["requested_root_hits"] > 0
            and replay["closed_capture_report"]["committed_records"] == 0
            and {
                source["native_index_identity"]["sha256"]
                for source in heldout["new_nonempty_capture_sources"]
            }
            <= {
                source["native_index_identity"]["sha256"]
                for source in sources["sources"]
            },
            "held-out exact output, replay or HEAVY capture membership differs",
        )
        checkpoint["heldout_correctness"] = {
            "case_id": heldout["case_id"],
            "source_commit": heldout["actual_source"],
            "public_guard_set_count": first["public_guard_set_count"],
            "master_count": len(first["used_masters"]),
            "used_family_count": len(first["used_family_sha256"]),
            "captured_records": first["closed_capture_report"]["committed_records"],
            "replay_requested_root_hits": replay["cache_statistics"]["persistent"][
                "requested_root_hits"
            ],
            "replay_rule_applications": replay["cache_statistics"]["native"][
                "rule_applications"
            ],
            "strict_capture_receipts": captures,
        }
    if heavy_records <= keys:
        checkpoint["final_heavy"] = final_heavy_cache(records)
        for key in tuple(checkpoint["pending"]):
            if key.startswith("final_heavy_"):
                del checkpoint["pending"][key]
    if q0_records <= keys:
        reference_report = study["records"]["heavy"]["raw_receipt"][
            "full_actual_benchmark_metadata"
        ]
        checkpoint["final_cached_q0"] = final_cached_q0(
            records, checkpoint["final_heavy"], reference_report
        )
    if "acceptance" in index:
        require(
            q0_records | heavy_records <= keys
            and set(index["acceptance"])
            == {"heavy", "q0", "heavy_handoff", "q0_handoff"},
            "final acceptance requires complete HEAVY and q0 evidence",
        )
        authority = {
            name: json_record(item, index_path.parent)
            for name, item in index["acceptance"].items()
        }
        native = authority["heavy"]["raw_receipt"]
        science = authority["q0"]["raw_receipt"]
        require(
            native["schema"] == "vakint.root-final13-native-cache-acceptance.v1"
            and native["status"]
            == "PASS_ACTUAL_FINAL13_CENSUS_EXPORT_EXHAUSTIVE_REIMPORT"
            and native["artifact"]["sha256"]
            == checkpoint["final_heavy"]["artifact"]["sha256"]
            and native["consolidated_handoff"]["sha256"]
            == authority["heavy_handoff"]["evidence"]["sha256"]
            and science["schema"]
            == "vakint.root-final13-cached-q0-full-functional-acceptance.v1"
            and science["status"]
            == "PASS_ACTUAL_FULL55_58_146678_Q0_SCIENCE_MONITOR_AND_CACHE"
            and science["handoff"]["sha256"]
            == authority["q0_handoff"]["evidence"]["sha256"]
            and science["actual_small_proofs"]["full_science_proof"]["sha256"]
            == records["final_q0_verification"]["evidence"]["sha256"],
            "final Root acceptance differs from the actual cache/science proof",
        )
        checkpoint["accepted_completion_authorities"] = authority
    return checkpoint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--campaign", type=Path)
    source.add_argument("--engineering", type=Path)
    source.add_argument("--file-engineering", type=Path)
    source.add_argument("--delivery-checkpoint", type=Path)
    parser.add_argument("--physical-scope", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--select-current", action="store_true")
    args = parser.parse_args()
    require(not args.output.exists(), "output must be a new path")
    require(
        not (args.engineering or args.file_engineering or args.delivery_checkpoint)
        or not args.physical_scope,
        "engineering evidence cannot change physical completion scope",
    )
    physical_scope, scope_authorization = load_physical_scope(args.physical_scope)
    merged = copy.deepcopy(read(args.base))
    if args.delivery_checkpoint:
        run = delivery_checkpoint(args.delivery_checkpoint, merged)
    elif args.file_engineering:
        run = cache_file_engineering(args.file_engineering, merged)
    elif args.engineering:
        run = engineering_comparison(args.engineering)
    else:
        run = physical_capture(args.campaign, physical_scope)
    if scope_authorization is not None:
        run.update(
            physical_scope_authorization=scope_authorization,
            notes=[physical_scope_note(scope_authorization)],
        )
    identity = run["evidence"]
    if scope_authorization is not None:
        identity = [identity, scope_authorization["evidence"]]
    run_id = (
        "q0-delivery-"
        if args.delivery_checkpoint
        else "cache-file-engineering-"
        if args.file_engineering
        else "q0-engineering-"
        if args.engineering
        else "physical-capture-"
    ) + hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:24]
    namespace = (
        "delivery_checkpoints"
        if args.delivery_checkpoint
        else "engineering_optimization_runs"
        if args.engineering or args.file_engineering
        else "physical_capture_runs"
    )
    runs = merged.setdefault(namespace, {})
    require(run_id not in runs, "this evidence is already appended")
    runs[run_id] = run
    if args.select_current:
        selector = (
            "selected_delivery_checkpoint"
            if args.delivery_checkpoint
            else "selected_cache_file_engineering_run"
            if args.file_engineering
            else "selected_engineering_optimization_run"
            if args.engineering
            else "selected_physical_capture_run"
        )
        merged[selector] = run_id
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(merged, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(
        json.dumps(
            {
                "run_id": run_id,
                "output": str(args.output),
                "sha256": digest(args.output),
            }
        )
    )


if __name__ == "__main__":
    main()
