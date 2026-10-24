"""Append existing, completed 5L cache/suite evidence to a NEW workspace report.

This is a report adapter, not an authenticator or a reduction runner. The native
drivers own proof validation. Functional records never establish performance.
Nothing is written unless every supplied completed receipt passes the checks.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import statistics
from datetime import datetime, timedelta
from pathlib import Path

PRIMARY = {
    "quadratic-four-rung-uv": "quadratic-transverse",
    "zero-momentum-diagnostic": "zero-momentum-metric",
}
METRICS = (
    "data_loading_seconds",
    "cache_loading_seconds",
    "owner_proof_validation_seconds",
    "owner_preparation_seconds",
    "saved_application_including_owner_seconds",
    "ibp_excluding_owner_seconds",
    "evaluation_seconds",
    "end_to_end_seconds",
)
STAGES = (
    "frontend_canonicalization",
    "frontend_tensor_prepass",
    "frontend_admitted_witness_search",
    "parent_transport_and_scalar_lowering",
    "saved_rule_application",
    "terminal_geometry_and_lower_factors",
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def evidence(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": digest(path)}


def hashed_file(reference, base):
    require(isinstance(reference, dict), "missing hashed evidence reference")
    path = Path(reference["path"])
    path = path if path.is_absolute() else base / path
    require(digest(path) == reference["sha256"], f"changed evidence: {path}")
    return path.resolve()


def hashed_record(reference, base):
    path = hashed_file(reference, base)
    return path, read(path)


def artifact_binaries(artifact):
    """Use the compiler receipt's actual binaries, never a runtime source label."""
    binary = artifact.get("binary")
    if isinstance(binary, dict):
        return [binary]
    binaries = artifact.get("artifacts")
    require(isinstance(binaries, list) and binaries, "missing compiler artifacts")
    return binaries


COHORT_STAMPS = (
    ("main_revision", "VAKINT_CACHE_BUILD_GIT_REVISION"),
    ("native_revision", "VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
    ("kernel_revision", "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
    ("main_source_snapshot_sha256", "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
    (
        "native_source_snapshot_sha256",
        "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256",
    ),
)


def selected_inputs(workload):
    return {
        "source_census_sha256": workload["source_census_sha256"],
        **{
            kind: [
                {
                    "id": row["id"],
                    "powers": row["powers"],
                    "public_input_source_sha256": hashlib.sha256(
                        row["public_input"].encode()
                    ).hexdigest()
                    if row.get("public_input") is not None
                    else None,
                }
                for row in workload[kind]
            ]
            for kind in ("training", "tuning", "held_out")
        },
    }


class ReportEvidence:
    """Bind raw reports to existing launch/build and completed window receipts.

    A window receipt is an external acceptance decision, never manufactured by
    this adapter. Its explicit observations must agree with the referenced
    launchers and artifacts. Merely hashing a functional run does not qualify it.
    """

    def __init__(self, path):
        self.path = Path(path).resolve()
        self.index = read(path)
        require(
            self.index.get("schema") == "vakint.five-loop-report-evidence.v1",
            "unknown report evidence schema",
        )
        self.records = {}
        self.used = set()
        for item in self.index["records"]:
            record_path, record = hashed_record(item["record"], self.path.parent)
            key = str(record_path)
            require(key not in self.records, "duplicate report evidence")
            launch_path, launcher = hashed_record(item["launcher"], self.path.parent)
            artifact_path, artifact = hashed_record(item["artifact"], self.path.parent)
            self.records[key] = {
                "record": record,
                "launcher": launcher,
                "artifact": artifact,
                "references": {
                    "record": evidence(record_path),
                    "launcher": evidence(launch_path),
                    "artifact": evidence(artifact_path),
                },
                "window": item.get("window"),
                "case_id": item["case_id"],
            }

    def bind(self, path):
        record = self.records.get(str(Path(path).resolve()))
        require(record is not None, "report has no launcher/artifact evidence")
        self.used.add(str(Path(path).resolve()))
        return copy.deepcopy(record)

    def cache_qualification(self, identity, producer):
        """Check an accepted provenance bridge; never redo native authentication."""
        build = producer["raw_report"]
        require(
            all(
                identity[key] == build[key]
                for key in (
                    "cache_sha256",
                    "data_manifest_sha256",
                    "source_census_sha256",
                )
            ),
            "build/suite cache, data or census mismatch",
        )
        if (
            build.get("phase") != "pilot"
            and identity["workload_sha256"] == build["workload_sha256"]
            and identity["owner_proof"] == build["owner_proof"]
            and identity["build_git_revision"]
            == build["build_identity"]["build"]["git_revision"]
        ):
            return None
        path, accepted = hashed_record(
            self.index.get("cache_requalification"), self.path.parent
        )
        require(
            accepted.get("schema") == "vakint.cache-consumer-requalification.v1"
            and accepted.get("accepted") is True
            and all(
                accepted.get("checks", {}).get(check) is True
                for check in (
                    "canonical_consumer_proof",
                    "unchanged_selection",
                    "twelve_distinct_held_out_sets",
                    "training_raw_and_nonzero_owner_disjoint",
                    "cache_import_validated",
                )
            ),
            "cross-cohort cache reuse has no accepted requalification",
        )
        documents = {}
        for role, names in (
            (
                "producer",
                ("cache", "build_report", "workload", "owner_proof", "compiled_cohort"),
            ),
            (
                "consumer",
                (
                    "collector_artifacts",
                    "compiled_cohort",
                    "qualification_workload",
                    "final_workload",
                    "owner_proof",
                    "protected_union",
                    "training_route_qualification",
                    "alias_exclusion",
                ),
            ),
        ):
            documents[role] = {}
            for name in names:
                source = hashed_file(accepted[role][name], path.parent)
                if name != "cache":  # Native cache bytes are not a JSON document.
                    documents[role][name] = read(source)
        old, new = documents["producer"], documents["consumer"]
        require(
            accepted["producer"]["build_report"]["sha256"]
            == producer["evidence"]["sha256"]
            and old["build_report"] == build
            and accepted["producer"]["cache"]["sha256"]
            == producer["cache_file"]["sha256"]
            == identity["cache_sha256"],
            "requalification names another producer/cache",
        )
        binding = build["preparation"]["options"]["single_file"]["binding"]
        require(
            isinstance(binding, dict)
            and binding
            and accepted["semantic_binding"] == binding,
            "requalification native binding differs",
        )
        for role, doc, workload_name, expected in (
            ("producer", old, "workload", build),
            ("consumer", new, "final_workload", identity),
        ):
            workload = doc[workload_name]
            proof = doc["owner_proof"]
            revision = (
                expected["build_identity"]["build"]["git_revision"]
                if role == "producer"
                else expected["build_git_revision"]
            )
            require(
                doc["compiled_cohort"]["main_revision"] == revision,
                "requalification recorded revision differs from its own cohort",
            )
            require(
                accepted[role][workload_name]["sha256"] == expected["workload_sha256"]
                and accepted[role]["owner_proof"]["sha256"]
                == expected["owner_proof"]["sha256"]
                == workload["owner_proof"]["sha256"],
                "requalification workload/proof identity differs",
            )
            require(
                proof["compiled_cohort"]
                == doc["compiled_cohort"]
                == expected["owner_proof"]["compiled_cohort"]
                and proof["semantic_binding"]
                == doc["compiled_cohort"]["semantic_binding"]
                == binding,
                "requalification proof/cohort/binding differs",
            )
            require(
                workload["data_manifest_sha256"]
                == proof["data_manifest_sha256"]
                == expected["data_manifest_sha256"]
                and workload["source_census_sha256"]
                == expected["source_census_sha256"],
                "requalification changed frozen data/census",
            )
        selection = selected_inputs(old["workload"])
        require(
            selection
            == accepted["immutable_selection"]
            == selected_inputs(new["final_workload"])
            == selected_inputs(new["qualification_workload"]),
            "requalification changed frozen selection",
        )
        artifact = new["collector_artifacts"]
        require(
            all(
                artifact["stamps"].get(stamp) == new["compiled_cohort"][field]
                for field, stamp in COHORT_STAMPS
            )
            and new["compiled_cohort"]["main_revision"]
            == identity["build_git_revision"]
            and identity["executable_sha256"]
            == new["owner_proof"]["measurement_executable_sha256"]
            and identity["executable_sha256"]
            in {row["sha256"] for row in artifact_binaries(artifact)},
            "requalification consumer artifact/stamps differ",
        )
        return {"evidence": evidence(path), "raw_receipt": accepted}

    def final_selection(self, producer):
        """Keep a selected pilot's original build separate from the later decision."""
        path, selection = hashed_record(
            self.index.get("final_selection"), self.path.parent
        )
        require(
            selection.get("schema") == "vakint.cache-final-selection.v1"
            and selection.get("accepted") is True
            and selection.get("decision") == "retain_complete_candidate"
            and selection.get("producer_phase") == "pilot"
            and isinstance(selection.get("reason"), str)
            and selection["reason"].strip(),
            "pilot cache has no accepted final selection",
        )
        require(
            isinstance(selection.get("decided_utc"), str)
            and datetime.fromisoformat(selection["decided_utc"]).utcoffset()
            == timedelta(0),
            "final selection needs its actual UTC decision time",
        )
        refs = {
            key: evidence(hashed_file(selection.get(key), path.parent))
            for key in (
                "cache",
                "build_report",
                "producer_workload",
                "consumer_workload",
                "candidate_acceptance",
                "cache_requalification",
            )
        }
        bridge_path, bridge = hashed_record(
            self.index.get("cache_requalification"), self.path.parent
        )
        build = producer["raw_report"]
        require(
            refs["cache"] == producer["cache_file"]
            and refs["build_report"] == producer["evidence"]
            and refs["producer_workload"]["sha256"] == build["workload_sha256"]
            and refs["consumer_workload"]["sha256"]
            == bridge["consumer"]["final_workload"]["sha256"]
            and refs["cache_requalification"] == evidence(bridge_path),
            "final selection names another cache, build, workload or bridge",
        )
        candidate = read(refs["candidate_acceptance"]["path"])
        require(
            candidate.get("schema")
            == "vakint.strict-complete-small-candidate-acceptance.v1"
            and candidate.get("state")
            == "STRICT_COMPLETE_CANDIDATE_NOT_FINAL_SELECTION"
            and candidate["reports"]["build"]["sha256"]
            == producer["evidence"]["sha256"]
            and candidate["cache_file"]["sha256"] == producer["cache_file"]["sha256"]
            and candidate["workload_sha256"] == build["workload_sha256"]
            and candidate["selected_training_roots"]
            == candidate["completed_strict_roots"]
            == build["preparation"]["total_roots"],
            "final selection names another completed candidate",
        )
        _, workload = hashed_record(selection["consumer_workload"], path.parent)
        _, cohort = hashed_record(
            bridge["consumer"]["compiled_cohort"], bridge_path.parent
        )
        _, proof = hashed_record(bridge["consumer"]["owner_proof"], bridge_path.parent)
        identity = {
            "cache_sha256": refs["cache"]["sha256"],
            "data_manifest_sha256": workload["data_manifest_sha256"],
            "source_census_sha256": workload["source_census_sha256"],
            "workload_sha256": refs["consumer_workload"]["sha256"],
            "build_git_revision": cohort["main_revision"],
            "executable_sha256": proof["measurement_executable_sha256"],
            "owner_proof": {
                **bridge["consumer"]["owner_proof"],
                "compiled_cohort": cohort,
            },
        }
        # Reuse the full provenance/semantic checks, even for an unchanged cohort.
        qualification = self.cache_qualification(identity, producer)
        return {
            "evidence": evidence(path),
            "raw_receipt": selection,
            "cache_requalification": qualification,
        }

    def qualify(self, bundle, binary, case_id, variant, sample, supervision=None):
        launcher, artifact = bundle["launcher"], bundle["artifact"]
        require(bundle["case_id"] == case_id, "evidence belongs to another case")
        require(
            launcher.get("schema")
            not in (
                "vakint.primary-full-functional-run.v1",
                "vakint.q0-full-weighted-launch.v1",
            )
            and not launcher.get(
                "functional_only_not_isolated_three_run_performance", False
            )
            and launcher.get("final_performance_eligible") is True,
            "functional launcher is not final performance evidence",
        )
        require(
            launcher.get("fresh_process") is True
            and launcher.get("binary_sha256") == binary,
            "launcher does not identify this fresh executable",
        )
        require(
            launcher.get("case_id") == case_id
            and launcher.get("variant") == variant
            and launcher.get("sample") == sample,
            "launcher sample identity differs",
        )
        # Support the existing single-primary and multi-driver artifact receipts.
        binary_records = artifact_binaries(artifact)
        if isinstance(artifact.get("binary"), dict):
            binaries = [artifact["binary"]["sha256"]]
            codegen = artifact["actual_codegen"]
            require(
                str(codegen.get("optimization", codegen.get("opt_level"))) == "3"
                and codegen.get("fat_lto") is True
                and codegen.get("test") is False,
                "artifact is not the qualified release executable",
            )
        else:
            candidates = binary_records
            binaries = [
                item["sha256"]
                for item in candidates
                if item.get("normal_cli") is True
                and item["profile"]["opt_level"] == "3"
                and item["profile"]["test"] is False
                and item["actual_rustc"]["lto"] == "fat"
            ]
            binary_records = candidates
        for record in binary_records:
            path = record.get("immutable_path")
            if path and Path(path).exists():
                require(
                    digest(path) == record["sha256"], "available executable changed"
                )
        require(
            binary in binaries
            and artifact.get("source_drift") == []
            and artifact.get("source_unchanged", True) is True,
            "artifact binary/source mismatch",
        )
        stamps = artifact["stamps"]
        required = {
            "VAKINT_CACHE_BUILD_GIT_REVISION": 40,
            "VAKINT_CACHE_BUILD_RUSTRED_REVISION": 40,
            "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION": 40,
            "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256": 64,
            "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256": 64,
        }
        for field, length in required.items():
            value = stamps.get(field)
            require(
                isinstance(value, str)
                and len(value) == length
                and all(c in "0123456789abcdefABCDEF" for c in value),
                f"missing actual artifact provenance: {field}",
            )
        require(
            stamps.get("VAKINT_CACHE_BUILD_PROFILE") == "release",
            "artifact profile is not release",
        )
        window_path, window = hashed_record(bundle["window"], self.path.parent)
        require(
            window.get("schema") == "vakint.serialized-performance-window.v1"
            and window.get("final_performance_eligible") is True
            and window.get("accepted") is True
            and window.get("serialized") is True
            and window.get("all_processes_drained") is True
            and window.get("source_unchanged") is True,
            "no completed serialized performance window",
        )
        intervals = sorted(
            window["samples"], key=lambda row: row["started_unix_seconds"]
        )
        require(intervals, "empty performance window")
        require(
            len({row["report_sha256"] for row in intervals}) == len(intervals),
            "window duplicates a report",
        )
        window_start = number(window["started_unix_seconds"])
        window_end = number(window["finished_unix_seconds"])
        last_end = None
        matches = []
        for row in intervals:
            start, end = (
                number(row["started_unix_seconds"]),
                number(row["finished_unix_seconds"]),
            )
            require(
                window_start <= start < end <= window_end
                and (last_end is None or start >= last_end),
                "performance samples overlap or have invalid intervals",
            )
            require(
                row.get("exit_code") == 0 and row.get("fresh_process") is True,
                "window includes an unsuccessful or reused process",
            )
            last_end = end
            if row["report_sha256"] == bundle["references"]["record"]["sha256"]:
                matches.append(row)
        require(len(matches) == 1, "window omits or duplicates this report")
        row = matches[0]
        if supervision is not None:
            require(
                bundle["record"]["supervision"] == supervision,
                "sample supervision differs from its hashed report",
            )
            require(
                supervision.get("interval_scope")
                == "supervisor pre-spawn setup through supervisor-observed worker completion",
                "missing actual supervisor observation scope",
            )
            for field in ("child_pid", "supervisor_pid"):
                require(
                    type(supervision.get(field)) is int and supervision[field] > 0,
                    "missing actual supervisor process identity",
                )
            for field in ("started_unix_seconds", "finished_unix_seconds"):
                number(supervision.get(field))
            for field in (
                "child_pid",
                "supervisor_pid",
                "started_unix_seconds",
                "finished_unix_seconds",
                "interval_scope",
            ):
                require(
                    row.get(field) == launcher.get(field) == supervision[field],
                    "window/launcher differs from actual supervisor observations",
                )
        require(
            row["launcher_sha256"] == bundle["references"]["launcher"]["sha256"]
            and row["artifact_sha256"] == bundle["references"]["artifact"]["sha256"]
            and row["case_id"] == case_id
            and row["variant"] == variant
            and row["sample"] == sample,
            "window evidence binding differs",
        )
        affinity = window.get("cpu_affinity")
        threads = window.get("thread_environment")
        require(
            isinstance(affinity, str)
            and affinity
            and launcher.get("cpu_affinity") == affinity,
            "affinity not fixed by launcher/window",
        )
        require(
            isinstance(threads, dict)
            and all(
                threads.get(k) == "1"
                for k in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            )
            and launcher.get("thread_environment") == threads,
            "worker environment not fixed by launcher/window",
        )
        return {
            "window": evidence(window_path),
            "stamps": stamps,
            "binary_sha256": binary,
            "cpu_affinity": affinity,
            "thread_environment": threads,
        }


def number(value):
    require(
        type(value) in (int, float) and math.isfinite(value) and value >= 0,
        "missing, negative or nonfinite timing",
    )
    return value


def summarize(values):
    values = [number(value) for value in values]
    require(len(values) >= 3, "at least three actual samples required")
    return {
        "samples_seconds": values,
        "median_seconds": statistics.median(values),
        "min_seconds": min(values),
        "max_seconds": max(values),
    }


def build_record(path, cache, provenance=None):
    build = read(path)
    require(
        build.get("schema") == "vakint.five-loop-cache-build.v2", "unknown build schema"
    )
    prep = build["preparation"]
    require(
        build["phase"] in ("final", "pilot")
        and build["complete"] is True
        and build["publishable"] is True
        and prep["status"] == "complete",
        "only completed cache preparation may enter a measured report",
    )
    require(
        build["supervision"]["success"] is True
        and build["supervision"]["hard_deadline_exhausted"] is False,
        "cache preparation lacks successful supervision",
    )
    require(
        prep["total_roots"] > 0 and prep["completed_roots"] == prep["total_roots"],
        "selected training roots are incomplete",
    )
    require(
        all(
            prep.get(field) is None
            for field in ("failure", "resource_failure", "checkpoint_failure")
        ),
        "cache preparation carries a failure",
    )
    require(
        prep["native_profile"]["strict"] is True
        and prep["native_profile"]["unresolved_policy"] == "reject",
        "cache records must be strict complete reductions",
    )
    require(
        build["cache_sha256"] == digest(cache)
        and build["compressed_bytes"] == Path(cache).stat().st_size,
        "cache file differs from completed build",
    )
    require(
        build["within_50_mib_target"] is True, "size exception needs separate review"
    )
    producer = {
        "evidence": evidence(path),
        "cache_file": evidence(cache),
        "raw_report": build,
    }
    if build["phase"] == "pilot":
        require(provenance is not None, "pilot cache requires final-selection evidence")
        producer["final_selection"] = provenance.final_selection(producer)
    return producer


def suite_record(path, kind, producer, provenance):
    build = producer["raw_report"]
    path = Path(path).resolve()
    suite = read(path)
    require(
        suite.get("schema") == "vakint.five-loop-cache-suite.v1", "unknown suite schema"
    )
    identity = suite["identity"]
    require(
        identity["suite"] == kind
        and suite["complete"] is True
        and suite["provisional"] is False
        and suite["native_keys_authenticated"] is True
        and suite["execution_provenance_complete"] is True,
        "suite is not qualified",
    )
    require(
        suite[
            "final_held_out_comparison" if kind == "held-out" else "selection_allowed"
        ]
        is True,
        "driver did not qualify this suite",
    )
    count = identity["repetitions"]
    require(type(count) is int and count >= 3, "three repetitions required")
    requalification = provenance.cache_qualification(identity, producer)
    require(
        identity["build_profile"] == "release"
        and identity["rayon_num_threads"] == "1"
        and identity["cpu_affinity"]
        and identity["load_prepared_caches"] is False,
        "final measurement profile/worker/cache policy differs",
    )
    raw = suite["measurements"]
    ids = [case["id"] for case in suite["cases"]]
    require(
        len(ids) == len(set(ids)) and len(raw) == len(ids) * count * 2,
        "missing or duplicate case/variant samples",
    )
    require(
        kind != "held-out" or len(ids) >= 12, "at least twelve held-out cases required"
    )
    rows, seen, root_sets, contexts = [], set(), set(), set()
    for index, case in enumerate(suite["cases"]):
        case_id = case["id"]
        require(
            case["successful"] is True
            and case["exact_output_and_guard_parity"] is True
            and case["native_keys_authenticated"] is True,
            "case failed exact validation",
        )
        samples = [row for row in raw if row["case_id"] == case_id]
        require(len(samples) == 2 * count, "case sample count differs")
        require(
            len({row["input_sha256"] for row in samples}) == 1
            and len({row["value_sha256"] for row in samples}) == 1
            and len({row["guards_sha256"] for row in samples}) == 1,
            "case input/output/guards differ across samples or modes",
        )
        protected = samples[0]["frozen_protected_roots"]
        require(
            protected and all(len(key) == 15 for key in protected), "empty owner set"
        )
        root_sets.add(tuple(sorted(tuple(key) for key in protected)))
        require(
            isinstance(case.get("catalogue_closure"), dict),
            "final collector must report catalogue closure separately from Retain success",
        )
        output = {
            "id": case_id,
            "status": "measured",
            "catalogue_closure": case["catalogue_closure"],
            "requested_root_hits": case["requested_root_hits"],
            "intermediate_hits": case["intermediate_hits"],
            "backing_reuse": case["backing_reuse"],
            "variants": {},
        }
        for sample in samples:
            key = (case_id, sample["variant"], sample["sample"])
            require(
                key not in seen
                and sample["variant"] in ("uncached", "cached")
                and sample["sample"] in range(count),
                "duplicate or invalid sample",
            )
            seen.add(key)
            require(
                sample["identity"] == identity
                and sample["success"] is True
                and sample["worker_success"] is True
                and sample["fresh_process"] is True
                and sample["native_keys_authenticated"] is True
                and sample["supervision"]["success"] is True
                and sample["supervision"]["exit_code"] == 0
                and sample["supervision"]["deadline_exhausted"] is False,
                "unsuccessful or mismatched sample",
            )
            require(
                sample["owner_proof"] == identity["owner_proof"], "owner proof mismatch"
            )
            require(
                sample["frozen_protected_roots"] == protected
                and sample["actual_lowered_roots"] == sample["frozen_lowered_roots"],
                "actual root set differs",
            )
            before = sample["persistent_cache_attachments_before"]
            during = sample["persistent_cache_attachments_during_evaluation"]
            require(
                before == {str(loop): False for loop in range(1, 6)},
                "baseline cache present",
            )
            require(
                during
                == {
                    str(loop): sample["variant"] == "cached" and loop == 5
                    for loop in range(1, 6)
                },
                "wrong persistent backing",
            )
            require(
                sample["cache_statistics"]["binding"]
                == build["preparation"]["options"]["single_file"]["binding"],
                "cache semantic binding mismatch",
            )
            directory = (
                path.parent
                / f"case-{index:03}-{sample['sample']:02}-{sample['variant']}"
            )
            require(
                read(directory / "measurement.json") == sample, "saved sample differs"
            )
            bundle = provenance.bind(directory / "measurement.json")
            if requalification is not None:
                require(
                    bundle["references"]["artifact"]["sha256"]
                    == requalification["raw_receipt"]["consumer"][
                        "collector_artifacts"
                    ]["sha256"],
                    "sample artifact differs from accepted consumer",
                )
            qualified = provenance.qualify(
                bundle,
                identity["executable_sha256"],
                case_id,
                sample["variant"],
                sample["sample"],
                supervision=sample["supervision"],
            )
            require(
                qualified["stamps"]["VAKINT_CACHE_BUILD_GIT_REVISION"]
                == identity["build_git_revision"]
                and qualified["cpu_affinity"] == identity["cpu_affinity"],
                "suite/artifact/launcher context differs",
            )
            cohort = identity["owner_proof"]["compiled_cohort"]
            for field, stamp in COHORT_STAMPS:
                require(
                    cohort[field] == qualified["stamps"][stamp],
                    "suite compiled cohort differs",
                )
            contexts.add(json.dumps(qualified, sort_keys=True))
            output.setdefault("performance_evidence", []).append(qualified)
            for filename, field in (
                ("value.symbolica", "value_sha256"),
                ("guards.json", "guards_sha256"),
            ):
                require(
                    digest(directory / filename) == sample[field],
                    "result/guard sidecar changed",
                )
            closure = sample["catalogue_closure"]
            unknown = closure["unknown_masters"]
            require(
                unknown == sorted(set(unknown))
                and closure["unknown_count"] == len(unknown)
                and closure["closed"] is (not unknown)
                and read(directory / "unknown-masters.json") == unknown
                and digest(directory / "unknown-masters.json")
                == closure["unknown_sha256"],
                "unknown-master classification or sidecar changed",
            )
        unknown_union = sorted(
            {
                master
                for sample in samples
                for master in sample["catalogue_closure"]["unknown_masters"]
            }
        )
        require(
            output["catalogue_closure"]["unknown_masters"] == unknown_union
            and output["catalogue_closure"]["unknown_count"] == len(unknown_union)
            and output["catalogue_closure"]["closed"] is (not unknown_union),
            "case catalogue closure differs from samples",
        )
        for variant in ("uncached", "cached"):
            selected = sorted(
                (row for row in samples if row["variant"] == variant),
                key=lambda row: row["sample"],
            )
            metrics = {
                field: summarize(row[field] for row in selected) for field in METRICS
            }
            stages = {
                field: summarize(
                    sum(number(v) for v in row["stages"][field]) for row in selected
                )
                if all(field in row["stages"] for row in selected)
                else None
                for field in STAGES
            }
            expected = metrics["ibp_excluding_owner_seconds"]
            reported = case[variant]
            require(
                all(
                    expected[key] == reported[key]
                    for key in ("median_seconds", "min_seconds", "max_seconds")
                ),
                "driver median/range mismatch",
            )
            metrics.update(
                graph_preparation_seconds=None,
                pure_ibp_seconds=None,
                terminal_substitution_seconds=None,
            )
            output["variants"][variant] = {"metrics": metrics, "stages": stages}
        cached_samples = [row for row in samples if row["variant"] == "cached"]
        for field in ("requested_root_hits", "intermediate_hits"):
            require(
                sum(
                    row["cache_statistics"]["persistent"][field]
                    for row in cached_samples
                )
                == case[field],
                "driver cache-hit aggregate differs",
            )
        u = output["variants"]["uncached"]["metrics"]["ibp_excluding_owner_seconds"][
            "median_seconds"
        ]
        c = output["variants"]["cached"]["metrics"]["ibp_excluding_owner_seconds"][
            "median_seconds"
        ]
        output["uncached_over_cached_ratio"] = u / c if c > 0 else None
        require(
            output["uncached_over_cached_ratio"]
            == case.get("uncached_over_cached_ratio"),
            "driver ratio mismatch",
        )
        rows.append(output)
    require(
        len(contexts) == 1,
        "suite samples require one cohort and adopted performance window",
    )
    require(
        kind != "held-out" or len(root_sets) >= 12,
        "aliased owner sets cannot count as twelve",
    )
    return {
        "evidence": evidence(path),
        "identity": identity,
        "rows": rows,
        "raw_suite": suite,
        "cache_requalification": requalification,
        "metric": "saved application minus owner preparation",
        "not_separately_timed": [
            "pure IBP",
            "terminal substitution",
            "graph preparation",
        ],
    }


def primary_input(report, directory):
    """Bind persisted bytes/provenance; scientific validation belongs to the native loader."""
    mode = report.get("input_mode", "fresh")
    require(
        mode in ("fresh", "prepared")
        and report["fresh_native_preparation"] is (mode == "fresh"),
        "primary input/preparation mode differs",
    )
    persisted = report.get("prepared_input")
    if mode == "prepared":
        require(isinstance(persisted, dict), "prepared input requires its manifest")
        require(
            report.get("preparation") is None
            and report.get("preparation_seconds") is None,
            "historical preparation cannot become current-run work",
        )
        number(report.get("input_loading_seconds"))
    else:
        number(report["preparation_seconds"])
        require(
            report.get("input_loading_seconds") is None,
            "fresh preparation cannot report prepared-input loading",
        )
    input_path = Path(directory) / "vakint-input.symbolica"
    if persisted is not None:
        path, manifest = hashed_record(
            {
                "path": persisted["manifest_path"],
                "sha256": persisted["manifest_sha256"],
            },
            Path(directory),
        )
        require(
            path.name == "prepared-input.json"
            and manifest == persisted["manifest"]
            and manifest["schema"] == "vakint.qcd-ladder.prepared-input.v1"
            and manifest["definition"] == "massless-qcd-four-rung-native-uv-v1"
            and manifest["projection"] == report["projection"],
            "prepared manifest identity differs",
        )
        assembly = manifest["scientific_assembly_sha256"]
        require(
            isinstance(assembly, str)
            and len(assembly) == 64
            and all(c in "0123456789abcdefABCDEF" for c in assembly),
            "missing scientific assembly identity",
        )
        for name, entry in (
            ("vakint-input.symbolica", manifest["input"]),
            ("projector-conditions.symbolica", manifest["nonzero_conditions_file"]),
        ):
            payload = path.parent / name
            require(
                payload.stat().st_size == entry["bytes"]
                and digest(payload) == entry["sha256"],
                "prepared input/guard payload changed",
            )
        require(
            (path.parent / "projector-conditions.symbolica").read_text()
            == "\n".join(manifest["nonzero_conditions"])
            and report["prepared_nonzero_conditions"] == manifest["nonzero_conditions"],
            "prepared projector guards differ",
        )
        producer, consumer = (
            manifest["producer"]["fixture_sources_sha256"],
            report["consumer_fixture_sources_sha256"],
        )
        require(
            report["fixture_sources_sha256"] == producer
            and set(manifest["scientific_sources_sha256"])
            == {
                "physical.rs",
                "projector.rs",
                "vacuum.rs",
                "ladder.dot",
                "qcd.json",
                "phi3.json",
            }
            and all(
                producer.get(name) == consumer.get(name) == value
                for name, value in manifest["scientific_sources_sha256"].items()
            )
            and report["contributions"] == manifest["contributions"]
            and report["scalar_terms"] == manifest["scalar_terms"],
            "prepared producer/consumer definition differs",
        )
        number(manifest["preparation_seconds"])
        for value in manifest["timings"].values():
            number(value)
        if mode == "fresh":
            require(
                report["preparation_seconds"] == manifest["preparation_seconds"]
                and report["preparation"] == manifest["timings"],
                "fresh preparation receipt differs",
            )
        input_path = path.parent / "vakint-input.symbolica"
        require(
            report["input_sha256"] == manifest["input"]["sha256"],
            "prepared scalar identity differs",
        )
    require(digest(input_path) == report["input_sha256"], "primary input changed")
    return {"mode": mode, "prepared": persisted if mode == "prepared" else None}


def accepted_generator(reference, projection, base):
    """Read existing scientific acceptance; this does not issue a new acceptance."""
    path, accepted = hashed_record(reference, Path(base))
    expected = {
        "quadratic-transverse": (
            "vakint.current-quadratic-preparation.scientific-parity.v1",
            "PASS_ALL21_FRESH_CURRENT_DEPENDENCY_RAW_INPUT_PARITY",
            21,
        ),
        "zero-momentum-metric": (
            "vakint.current-q0-preparation.scientific-parity.v1",
            "PASS_Q0_FRESH_CURRENT_DEPENDENCY_RAW_INPUT_PARITY",
            1,
        ),
    }[projection]
    require(
        (accepted.get("schema"), accepted.get("status"), accepted.get("contributions"))
        == expected
        and accepted.get("collector_exit_code") == 0
        and accepted.get("source_drift") == []
        and accepted["input"]["exact_canonical_bytes_equal_accepted_input"] is True
        and accepted["projector_guards"]["bytes_equal_accepted"] is True,
        "primary generator acceptance is incomplete",
    )
    verified = evidence(path)
    if reference.get("provenance_addendum") is not None:
        extra_path, extra = hashed_record(reference["provenance_addendum"], Path(base))
        require(
            extra["schema"] == "vakint.current-q0-preparation.provenance-addendum.v1"
            and extra["status"]
            == "CORRECTS_FEYNKIT_OWNER_REVISION_ONLY_RAW_PARITY_UNCHANGED"
            and extra["original_scientific_parity"]["sha256"] == verified["sha256"]
            and extra["actual_compiled_workspace_revision"] == accepted["source_commit"]
            and extra["immutable_elf"]["sha256"] == accepted["immutable_elf_sha256"]
            and extra["input_sha256"] == accepted["input"]["sha256"]
            and extra["scientific_definition_and_raw_input_changes"] is False,
            "generator provenance addendum changed the scientific identity",
        )
        verified["provenance_addendum"] = evidence(extra_path)
    return verified, accepted


def primary_generator(report, launcher, base):
    """Bind the accepted physical input, not merely mutually equal samples."""
    reference, accepted = accepted_generator(
        launcher.get("generator_acceptance"), report["projection"], base
    )
    require(
        accepted["projection"] == report["projection"]
        and accepted["input"]["sha256"] == report["input_sha256"]
        and accepted["input"]["bytes"] == report["input_bytes"]
        and accepted["contributions"] == report["contributions"]
        and accepted["projector_guards"]["values"]
        == report["prepared_nonzero_conditions"]
        and accepted["projector_guards"]["sha256"]
        == hashlib.sha256(
            "\n".join(report["prepared_nonzero_conditions"]).encode()
        ).hexdigest(),
        "primary differs from the accepted complete physical fixture",
    )
    consumer = report.get(
        "consumer_fixture_sources_sha256", report["fixture_sources_sha256"]
    )
    require(
        all(
            consumer.get(name) == accepted["compiled_fixture_source_hashes"].get(name)
            for name in (
                "physical.rs",
                "projector.rs",
                "vacuum.rs",
                "ladder.dot",
                "qcd.json",
                "phi3.json",
            )
        ),
        "primary scientific sources differ from accepted generator",
    )
    if report.get("input_mode", "fresh") == "prepared":
        manifest = report["prepared_input"]["manifest"]
        require(
            manifest["producer"]["fixture_sources_sha256"]
            == accepted["compiled_fixture_source_hashes"]
            and manifest["producer"]["source_commit"] == accepted["source_commit"]
            and manifest["producer"]["executable_sha256"]
            == accepted["immutable_elf_sha256"],
            "prepared producer provenance differs from accepted generator",
        )
    return reference


def primary_record(bundle):
    """Retain original functional provenance; verify full outputs when present."""
    report = bundle["record"]
    case_id = bundle["case_id"]
    require(case_id in PRIMARY, "unknown primary case")
    schema = report.get("schema")
    require(
        schema
        in (
            "vakint.qcd-ladder.benchmark.v1",
            "vakint.q0-full-weighted-native-diagnostic.v2",
        ),
        "unknown primary report schema",
    )
    result = {
        "evidence": bundle["references"],
        "raw_report": report,
        "launcher": bundle["launcher"],
        "artifact": bundle["artifact"],
        "performance_eligible": False,
        "catalogue_closure": None,
        "status": "functional-incomplete",
    }
    if schema == "vakint.q0-full-weighted-native-diagnostic.v2":
        require(
            case_id == "zero-momentum-diagnostic" and bundle["window"] is None,
            "exploratory native q0 cannot become final primary performance",
        )
        artifact = bundle["artifact"]
        cohort = report["compiled_cohort"]
        require(
            all(
                cohort.get(key) == artifact.get(key) and cohort.get(key) is not None
                for key in (
                    "native",
                    "symbolica_revision",
                    "driver_source_sha256",
                    "input_sha256",
                    "data",
                    "native_profile",
                )
            )
            and report["input_sha256"] == artifact["input_sha256"],
            "exploratory q0 cohort differs",
        )
        bindings = bundle["launcher"]["bindings"]
        require(
            any(
                row["sha256"] == bundle["references"]["artifact"]["sha256"]
                for row in bindings
            )
            and any(row["sha256"] == artifact["binary_sha256"] for row in bindings),
            "exploratory launcher does not bind this artifact and executable",
        )
        result["status"] = "exploratory-functional"
        return result
    require(report["projection"] == PRIMARY[case_id], "primary projection differs")
    require(
        report["executable_sha256"]
        in {row["sha256"] for row in artifact_binaries(bundle["artifact"])}
        and report["executable_sha256"] == bundle["launcher"]["binary_sha256"],
        "primary report/launcher/artifact executable mismatch",
    )
    if report.get("retained_result_available") is not True:
        require(
            bundle["window"] is None,
            "incomplete primary cannot enter a performance window",
        )
        return result
    directory = Path(bundle["references"]["record"]["path"]).parent
    input_identity = primary_input(report, directory)
    result["input_mode"] = input_identity["mode"]
    result["prepared_input"] = input_identity["prepared"]
    require(
        digest(directory / "result.symbolica") == report["result_sha256"],
        "primary exact result changed",
    )
    closure = report["catalogue_closure"]
    unknown_path = directory / closure["unknown_file"]
    unknown = read(unknown_path)
    require(
        isinstance(unknown, list)
        and all(isinstance(item, str) for item in unknown)
        and unknown == sorted(set(unknown))
        and digest(unknown_path) == closure["unknown_sha256"]
        and closure["unknown_count"] == len(unknown)
        and closure["master_count"] >= len(unknown)
        and closure["closed"] is (not unknown),
        "invalid primary unknown inventory",
    )
    guards = report["guards"]
    guard_path = directory / guards["file"]
    require(
        guards["includes_preparation_conditions"] is True
        and digest(guard_path) == guards["sha256"]
        and len(guard_path.read_text().splitlines()) == guards["count"],
        "missing or changed complete primary guards",
    )
    result["catalogue_closure"] = {**closure, "unknown_masters": unknown}
    result["status"] = "functional-closed" if not unknown else "retained-incomplete"
    require(
        report["reduction_completed"] is (not unknown),
        "primary completion and declared catalogue closure disagree",
    )
    return result


def primary_records(provenance):
    output = {}
    for case_id in PRIMARY:
        bundles = [
            bundle
            for bundle in provenance.records.values()
            if bundle["case_id"] == case_id
        ]
        records = [primary_record(bundle) for bundle in bundles]
        provenance.used.update(
            bundle["references"]["record"]["path"] for bundle in bundles
        )
        item = {
            "status": "pending",
            "functional_records": records,
            "uncached": None,
            "cached": None,
            "catalogue_closure": None,
            "exact_output_and_guard_parity": None,
        }
        if records:
            item["status"] = "functional-only"
        candidates = [
            (bundle, record)
            for bundle, record in zip(bundles, records)
            if bundle["window"] is not None
        ]
        if not candidates:
            output[case_id] = item
            continue
        samples, signatures, contexts = {"uncached": {}, "cached": {}}, set(), set()
        for bundle, record in candidates:
            report, launcher = bundle["record"], bundle["launcher"]
            require(
                record["status"] == "functional-closed",
                "primary catalogue is incomplete",
            )
            variant, index = launcher["variant"], launcher["sample"]
            require(
                variant in samples
                and type(index) is int
                and index >= 0
                and index not in samples[variant],
                "duplicate or invalid primary sample",
            )
            qualified = provenance.qualify(
                bundle, report["executable_sha256"], case_id, variant, index
            )
            generator = primary_generator(
                report, launcher, Path(bundle["references"]["launcher"]["path"]).parent
            )
            record["generator_acceptance"] = generator
            require(
                report["source_commit"]
                == qualified["stamps"]["VAKINT_CACHE_BUILD_GIT_REVISION"]
                and report["source_dirty"] is False,
                "primary source/preparation context differs",
            )
            require(
                report["persistent_cache_attachments_before"]
                == {str(loop): False for loop in range(1, 6)}
                and report["persistent_cache_attachments_during_evaluation"]
                == {
                    str(loop): variant == "cached" and loop == 5 for loop in range(1, 6)
                },
                "primary persistent cache policy differs",
            )
            cache = launcher.get("cache_sha256")
            require(
                (variant == "uncached" and cache is None)
                or (
                    variant == "cached" and isinstance(cache, str) and len(cache) == 64
                ),
                "missing actual manual-cache identity",
            )
            if variant == "cached":
                require(
                    digest(report["persistent_cache"]) == cache, "primary cache changed"
                )
            mode = record["input_mode"]
            require(
                launcher.get("input_mode", "fresh") == mode,
                "primary launcher input mode differs",
            )
            manifest_sha = (
                record["prepared_input"]["manifest_sha256"]
                if mode == "prepared"
                else None
            )
            if mode == "prepared":
                require(
                    launcher.get("prepared_input_manifest_sha256") == manifest_sha,
                    "launcher did not bind the prepared input",
                )
            representation = report.get("coefficient_representation", "sparse")
            require(
                representation in ("sparse", "factorized"),
                "unknown coefficient representation",
            )
            signatures.add(
                json.dumps(
                    {
                        key: report[key]
                        for key in (
                            "projection",
                            "input_sha256",
                            "fixture_sources_sha256",
                            "data_manifest_sha256",
                            "result_sha256",
                            "guards",
                        )
                    },
                    sort_keys=True,
                )
            )
            contexts.add(
                json.dumps(
                    {
                        "stamps": qualified["stamps"],
                        "affinity": qualified["cpu_affinity"],
                        "threads": qualified["thread_environment"],
                        "representation": representation,
                        "input_mode": mode,
                        "generator_acceptance": generator,
                        "prepared_input_manifest_sha256": manifest_sha,
                        "consumer_fixture_sources_sha256": report.get(
                            "consumer_fixture_sources_sha256"
                        ),
                        "artifact": bundle["references"]["artifact"]["sha256"],
                        "window": qualified["window"]["sha256"],
                    },
                    sort_keys=True,
                )
            )
            record["performance_eligible"] = True
            record["performance_evidence"] = qualified
            samples[variant][index] = (report, cache)
        require(
            len(signatures) == 1 and len(contexts) == 1,
            "primary input/guards/result/representation/cohort parity differs",
        )
        require(
            set(samples["uncached"]) == set(samples["cached"])
            and len(samples["uncached"]) >= 3,
            "three matched fresh samples per mode required",
        )
        require(
            len({cache for _, cache in samples["cached"].values()}) == 1,
            "primary samples used different caches",
        )
        first_report = next(iter(samples["uncached"].values()))[0]
        item["generator_acceptance"] = generator
        item["input_mode"] = first_report.get("input_mode", "fresh")
        item["prepared_input"] = (
            first_report["prepared_input"] if item["input_mode"] == "prepared" else None
        )
        for variant, indexed in samples.items():
            reports = [indexed[index][0] for index in sorted(indexed)]
            for row in reports:
                for field in (
                    "data_loading_seconds",
                    "reduction_seconds",
                    "end_to_end_seconds",
                ):
                    number(row[field])
            metric_fields = (
                "preparation_seconds",
                "input_loading_seconds",
                "data_loading_seconds",
                "cache_loading_seconds",
                "reduction_seconds",
                "end_to_end_seconds",
            )
            metrics = {
                key: summarize(row[key] for row in reports)
                if all(row.get(key) is not None for row in reports)
                else None
                for key in metric_fields
            }
            metrics["owner_preparation_seconds"] = (
                summarize(
                    row["owner_preparation"]["preparation_during_reduction_seconds"]
                    for row in reports
                )
                if all(
                    row.get("owner_preparation", {}).get(
                        "preparation_during_reduction_seconds"
                    )
                    is not None
                    for row in reports
                )
                else None
            )
            metrics.update(
                pure_ibp_seconds=None,
                terminal_substitution_seconds=None,
                ibp_excluding_owner_seconds=None,
            )
            stages = {
                field: summarize(
                    sum(number(v) for v in row["stages"][field]) for row in reports
                )
                if all(field in row.get("stages", {}) for row in reports)
                else None
                for field in STAGES
            }
            preparation = {
                field: summarize(row["preparation"][field] for row in reports)
                if all(
                    (row.get("preparation") or {}).get(field) is not None
                    for row in reports
                )
                else None
                for field in (
                    "graph_generation_seconds",
                    "ingress_seconds",
                    "scalar_contraction_seconds",
                    "taylor_seconds",
                    "angular_seconds",
                )
            }
            item[variant] = {
                "metrics": metrics,
                "stages": stages,
                "preparation": preparation,
            }
        item["cache_sha256"] = next(iter(samples["cached"].values()))[1]
        item["data_manifest_sha256"] = reports[0]["data_manifest_sha256"]
        item.update(
            status="measured",
            catalogue_closure=True,
            exact_output_and_guard_parity=True,
        )
        output[case_id] = item
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--build", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--held-out", type=Path)
    parser.add_argument("--tuning", type=Path)
    parser.add_argument(
        "--evidence",
        type=Path,
        required=True,
        help="Hashed report/launcher/artifact/window index; functional records need no window",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--select-current", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    require(
        not output.exists(),
        "output must be a new path; never overwrite canonical/history files",
    )
    base = read(args.base)
    merged = copy.deepcopy(base)
    require(
        bool(args.build) == bool(args.cache),
        "build and cache must be supplied together",
    )
    require(
        not (args.held_out or args.tuning) or args.build,
        "suite measurements require their strict cache build",
    )
    provenance = ReportEvidence(args.evidence)
    require(args.build or provenance.records, "no report evidence supplied")
    build = build_record(args.build, args.cache, provenance) if args.build else None
    run = {
        "schema": "vakint.five-loop-report-run.v1",
        "cache": build,
        "held_out": suite_record(args.held_out, "held-out", build, provenance)
        if args.held_out
        else None,
        "tuning": suite_record(args.tuning, "tuning", build, provenance)
        if args.tuning
        else None,
        "primary": primary_records(provenance),
        "evidence_index": evidence(args.evidence),
        "scope": "Functional evidence is retained separately from eligible matched performance samples.",
    }
    require(
        provenance.used == set(provenance.records), "unused or unknown report evidence"
    )
    for primary in run["primary"].values():
        if primary["status"] == "measured":
            require(
                build is not None,
                "primary performance requires its completed strict cache build",
            )
            require(
                primary["cache_sha256"] == build["raw_report"]["cache_sha256"]
                and primary["data_manifest_sha256"]
                == build["raw_report"]["data_manifest_sha256"],
                "primary performance cache/data differs from completed build",
            )
    identity = json.dumps(
        [
            run["evidence_index"],
            build["evidence"] if build else None,
            run["held_out"]["evidence"] if run["held_out"] else None,
            run["tuning"]["evidence"] if run["tuning"] else None,
        ],
        sort_keys=True,
    )
    run_id = "five-loop-" + hashlib.sha256(identity.encode()).hexdigest()[:24]
    runs = merged.setdefault("five_loop_runs", {})
    require(run_id not in runs, "this evidence is already appended")
    runs[run_id] = run
    if args.select_current:
        merged["selected_five_loop_run"] = run_id
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        json.dump(merged, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(
        json.dumps({"run_id": run_id, "output": str(output), "sha256": digest(output)})
    )


if __name__ == "__main__":
    main()
