"""Append completed physical capture/reuse evidence without qualifying performance.

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

from merge_five_loop_results import digest, evidence, hashed_file, read, require


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


def physical_capture(campaign):
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
    return {
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--select-current", action="store_true")
    args = parser.parse_args()
    require(not args.output.exists(), "output must be a new path")
    run = physical_capture(args.campaign)
    merged = copy.deepcopy(read(args.base))
    run_id = (
        "physical-capture-"
        + hashlib.sha256(
            json.dumps(run["evidence"], sort_keys=True).encode()
        ).hexdigest()[:24]
    )
    runs = merged.setdefault("physical_capture_runs", {})
    require(run_id not in runs, "this evidence is already appended")
    runs[run_id] = run
    if args.select_current:
        merged["selected_physical_capture_run"] = run_id
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
