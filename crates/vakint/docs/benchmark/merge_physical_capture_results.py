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
            proof = records["coefficient_proof"]["raw_receipt"]
            require(
                verified["full58_proof"]["sha256"]
                == records["coefficient_proof"]["evidence"]["sha256"]
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
                    row["exact_difference_zero"] is True
                    for row in proof["full_scalar_proofs"]
                )
                and proof["global_together_or_full_amplitude_difference"] is False
                and proof["unit_mass_or_coupling_replacement"] is False
                and proof["numerical_sampling"] is False
                and report["cache_hit_diagnostic_enabled"] == (name == "after_on"),
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


def delivery_checkpoint(index_path, base):
    """Record delivery decisions and functional evidence, never performance medians."""
    index_path = Path(index_path).resolve()
    index = read(index_path)
    require(
        index["schema"] == "vakint.q0-delivery-evidence-index.v1"
        and set(index["records"])
        == {"light_freeze", "monitor", "authorization", "heavy_sources"},
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
    return {
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
