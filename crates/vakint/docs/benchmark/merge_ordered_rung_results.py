"""Append checked ordered/DFS physical-fixture records without rewriting history.

This adapter neither runs reductions nor proves symbolic equality. The fixture
or an independent comparator owns equality; its receipt must bind both reports.
Timings remain diagnostic unless an independently supplied protocol qualifies
matched fresh-process release measurements. Missing stages stay null.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import statistics
from datetime import datetime
from itertools import pairwise
from pathlib import Path

MATH_SCOPE = "vakint.rustred.reduction-math.v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def evidence(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": digest(path)}


def referenced(reference, base):
    require(isinstance(reference, dict), "missing hashed file reference")
    path = Path(reference["path"])
    path = path if path.is_absolute() else base / path
    require(path.is_file(), f"missing evidence: {path}")
    require(digest(path) == reference["sha256"], f"changed evidence: {path}")
    return path.resolve()


def seconds(value):
    require(
        type(value) in (int, float) and math.isfinite(value) and value >= 0,
        "invalid measured duration",
    )
    return value


def summary(values):
    return {
        "samples": values,
        "median_seconds": statistics.median(values),
        "min_seconds": min(values),
        "max_seconds": max(values),
    }


def reduction_math(report):
    """Return checked net math time, never infer it from public-call duration."""
    value = report.get("reduction_math_seconds")
    timing = report.get("reduction_math_timing")
    if timing is None:
        require(value is None, "math duration lacks its timing scope")
        return None
    require(
        timing["scope_id"] == MATH_SCOPE
        and type(timing["complete"]) is bool
        and type(timing["interval_count"]) is int
        and timing["interval_count"] >= 0,
        "invalid reduction-math scope/count",
    )
    gross = timing["gross_seconds"]
    preparation = timing["measured_owner_preparation_seconds"]
    for observed in (gross, preparation):
        if observed is not None:
            seconds(observed)
    if not timing["complete"]:
        require(
            value is None
            and isinstance(timing.get("incompleteness_reason"), str)
            and timing["incompleteness_reason"],
            "incomplete math telemetry must retain a reason and no net duration",
        )
        return None
    require(
        gross is not None
        and preparation is not None
        and preparation <= gross
        and (timing["interval_count"] > 0 or gross == 0)
        and type(report.get("worker_policy", {}).get("owner_preparation_workers"))
        is int
        and report.get("worker_policy", {}).get("owner_preparation_workers") == 1
        and timing.get("incompleteness_reason") is None,
        "complete math timing has invalid owner-preparation deduction or worker policy",
    )
    seconds(value)
    # Preserve Duration's nanosecond resolution and f64 serialization rounding.
    rounding = (
        2
        * max(
            math.ulp(x)
            for x in (value, gross, preparation, report["retained_reduction_seconds"])
        )
        + 1e-9
    )
    require(
        abs(value - (gross - preparation)) <= rounding
        and gross <= report["retained_reduction_seconds"] + rounding,
        "net math duration differs from its checked non-overlapping scope",
    )
    return value


def validate_report(path, expected_strategy):
    raw = read(path)
    require(
        raw["schema"] == "vakint.ordered-rung-comparison.v1"
        and raw["strategy"] == expected_strategy,
        "unexpected comparison report or strategy",
    )
    require(
        type(raw["loop_count"]) is int
        and raw["loop_count"] > 0
        and raw["gluon_rungs"] == raw["loop_count"] - 1,
        "invalid physical loop/rung count",
    )
    require(
        raw["projection"] in ("zero-momentum-metric", "quadratic-transverse"),
        "unknown physical projection",
    )
    require(raw["persistent_cache"] is False, "ordered comparison attaches a cache")
    require(
        type(raw["monitoring_enabled"]) is bool, "missing effective monitoring policy"
    )
    require(raw["reduction_completed"] is True, "physical reduction is incomplete")
    closure = raw["catalogue_closure"]
    require(
        closure["closed"] is True
        and closure["unknown"] == []
        and closure.get("unknown_count", 0) == 0
        and type(closure["master_count"]) is int
        and closure["master_count"] >= 0,
        "physical reduction has no accepted catalogue closure",
    )
    for key in ("input_sha256", "evaluated_input_sha256"):
        require(
            re.fullmatch(r"[0-9a-f]{64}", raw[key]),
            f"malformed {key} identity",
        )
    binding = raw["dimension_binding"]
    require(
        all(
            isinstance(binding[key], str) and binding[key]
            for key in ("scope", "source", "target")
        ),
        "missing exact dimension binding",
    )
    require(
        raw["producer"]["build_profile"] in ("release", "dev", "ci-optim"),
        "missing explicit build profile",
    )
    payload = raw["retained"]
    retained = referenced(
        {"path": payload["file"], "sha256": payload["sha256"]}, path.parent
    )
    require(retained.stat().st_size == payload["bytes"], "retained size differs")
    require(
        raw["normalization"] == "MSbar",
        "new physical reports require MSbar normalization",
    )
    for key in ("loading_seconds", "retained_reduction_seconds", "completed_seconds"):
        seconds(raw[key])
    reduction_math(raw)
    if raw["numerical_reference"] is not None:
        require(raw["loop_count"] <= 4, "uncomputed 5L numerical masters claimed")
        seconds(raw["numerical_reference"]["seconds"])
        numerical = raw["numerical_reference"]
        require(numerical["coefficients"], "empty numerical reference")
        orders = [item["epsilon_power"] for item in numerical["coefficients"]]
        require(
            all(type(order) is int for order in orders)
            and orders == sorted(set(orders)),
            "numerical reference orders are invalid or repeated",
        )
        require(
            numerical["point"]
            and numerical["epsilon_range"]
            and type(numerical["selection"]["requested_decimal_precision"]) is int
            and numerical["selection"]["requested_decimal_precision"] > 0,
            "missing numerical reference point/order/precision contract",
        )
    return raw


def cpu_set(value):
    cpus = set()
    for part in value.split(","):
        limits = [int(item) for item in part.split("-")]
        require(
            1 <= len(limits) <= 2 and 0 <= limits[0] <= limits[-1],
            "invalid CPU affinity",
        )
        cpus.update(range(limits[0], limits[-1] + 1))
    require(cpus and min(cpus) >= 0, "empty or negative CPU affinity")
    return cpus


def checked_option(command, name):
    require(command.count(name) == 1, f"missing or repeated launch option {name}")
    index = command.index(name) + 1
    require(index < len(command), f"missing launch value {name}")
    return command[index]


def source_cap(numerical):
    """Compare capability, allowing equivalent exact constant presentations."""
    limits = []
    for source in numerical["selection"]["surviving_dependencies"]:
        availability = source["available_precision"]
        if availability == "ArbitraryPrecision":
            continue
        fixed = re.fullmatch(r"FixedPrecision\((\d+)\)", availability)
        require(
            fixed is not None and int(fixed[1]) > 0,
            "unknown source precision capability",
        )
        limits.append(int(fixed[1]))
    return min(limits, default=None)


def validate_protocol(path, pairs):
    """Qualify timing scope, not infer a release profile from debug assertions."""
    raw = read(path)
    require(raw["schema"] == "vakint.ordered-rung-protocol.v1", "unknown protocol")
    require(
        raw["fresh_processes"] is True
        and raw["serialized_execution"] is True
        and raw["build_profile"] == "release"
        and raw["profile_overrides"] is False,
        "performance requires matched fresh serialized normal-release samples",
    )
    require(raw["cpu_affinity"], "missing fixed worker affinity")
    require(
        all(
            isinstance(raw["thread_environment"].get(key), str)
            and raw["thread_environment"][key]
            for key in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        ),
        "missing worker thread settings",
    )
    referenced(raw["executable"], path.parent)
    sources = raw["sources"]
    for name in ("main_revision", "native_revision", "kernel_revision"):
        require(re.fullmatch(r"[0-9a-f]{40}", sources[name]), f"invalid {name}")
    for name in (
        "main_source_tree_sha256",
        "native_source_tree_sha256",
        "cargo_lock_sha256",
    ):
        require(re.fullmatch(r"[0-9a-f]{64}", sources[name]), f"invalid {name}")
    build_path = referenced(raw["build_receipt"], path.parent)
    build = read(build_path)
    require(
        build["schema"] == "vakint.ordered-rung-build.v1"
        and build["success"] is True
        and build["build_profile"] == "release"
        and build["profile_overrides"] is False
        and build["executable"] == raw["executable"]
        and build["sources"] == sources
        and isinstance(build["command"], list)
        and all(isinstance(word, str) for word in build["command"])
        and "--release" in build["command"]
        and "--locked" in build["command"],
        "missing bound successful normal-release build receipt",
    )
    referenced(build["build_log"], build_path.parent)
    report_map = {
        record["evidence"]["sha256"]: record
        for pair in pairs
        for record in pair["reports"].values()
    }
    identities = set()
    accepted = set()
    windows = []
    for item in raw["samples"]:
        accepted.add(tuple(item[name]["report_sha256"] for name in ("dfs", "ordered")))
        for name in ("dfs", "ordered"):
            process = item[name]
            require(
                type(process["pid"]) is int
                and process["pid"] > 0
                and type(process["start_ticks"]) is int
                and process["start_ticks"] > 0
                and type(process["exit_code"]) is int
                and process["exit_code"] == 0
                and process["executable_sha256"] == raw["executable"]["sha256"],
                "missing successful original worker identity or executable binding",
            )
            identity = (process["pid"], process["start_ticks"])
            require(identity not in identities, "performance reuses a worker process")
            identities.add(identity)
            launch_path = referenced(process["launch"], path.parent)
            completion_path = referenced(process["completion"], path.parent)
            launch, completion = read(launch_path), read(completion_path)
            require(
                launch["schema"] == "vakint.ordered-rung-process.v1"
                and launch["guard_armed_before_native_exec"] is True
                and launch_path.parent == completion_path.parent
                and (launch["pid"], launch["start_ticks"]) == identity
                and launch["executable"]["sha256"] == raw["executable"]["sha256"]
                and set(launch["cpu_affinity"]) == cpu_set(raw["cpu_affinity"])
                and launch["thread_environment"] == raw["thread_environment"]
                and completion["success"] is True
                and type(completion["exit_code"]) is int
                and completion["exit_code"] == 0
                and completion["failure"] is None
                and completion["signals_sent_to_original_pidfd"] == [],
                "worker guard launch/completion binding differs",
            )
            require(process["report_sha256"] in report_map, "unbound protocol sample")
            report = report_map[process["report_sha256"]]
            command = launch["command"]
            projection = (
                "zero-momentum"
                if report["raw"]["projection"] == "zero-momentum-metric"
                else "quadratic"
            )
            require(
                Path(command[0]).resolve() == Path(raw["executable"]["path"]).resolve()
                and checked_option(command, "--engine") == name
                and checked_option(command, "--loops")
                == str(report["raw"]["loop_count"])
                and checked_option(command, "--projection") == projection
                and checked_option(command, "--build-profile") == "release"
                and command.count("--no-monitor") == 1
                and "--monitor" not in command
                and "--progress-snapshot" not in command
                and Path(checked_option(command, "--work-dir")).resolve()
                == Path(report["evidence"]["path"]).parent.parent,
                "guard command does not bind the reported case/engine/profile",
            )
            start, end = (
                datetime.fromisoformat(value)
                for value in (launch["utc"], completion["utc"])
            )
            require(
                start.tzinfo is not None and end.tzinfo is not None and end >= start,
                "worker has invalid serialized window timestamps",
            )
            windows.append((start, end))
    windows.sort()
    require(
        all(previous[1] <= following[0] for previous, following in pairwise(windows)),
        "qualified performance worker windows overlap",
    )
    require(len(accepted) == len(raw["samples"]), "duplicate qualified sample")
    for pair in pairs:
        require(
            (
                pair["reports"]["dfs"]["evidence"]["sha256"],
                pair["reports"]["ordered"]["evidence"]["sha256"],
            )
            in accepted,
            "protocol does not bind a supplied sample",
        )
        for item in pair["reports"].values():
            require(
                item["raw"]["monitoring_enabled"] is False,
                "qualified performance requires monitoring disabled",
            )
            require(
                item["raw"]["producer"]["debug_assertions"] is False
                and item["raw"]["producer"]["build_profile"] == "release"
                and item["raw"]["producer"]["commit"] == sources["main_revision"]
                and item["raw"]["producer"]["cargo_lock_sha256"]
                == sources["cargo_lock_sha256"]
                and item["raw"]["producer"]["executable_sha256"]
                == raw["executable"]["sha256"]
                and item["raw"]["producer"]["allowed_cpus"] == raw["cpu_affinity"]
                and item["raw"]["producer"]["comparison_process"]
                == "one engine in this process",
                "producer is not the matched fresh-process release cohort",
            )
            require(
                reduction_math(item["raw"]) is not None,
                "qualified performance requires complete reduction-math telemetry",
            )
            cache = item["raw"].get("cache_policy", {})
            require(
                cache.get("prepared_cache_loaded") is False
                and cache.get("previous_run_cache_restored") is False
                and isinstance(cache.get("depth_first_intrinsic_memoization"), str)
                and cache["depth_first_intrinsic_memoization"]
                and cache.get("ordered_result_cache") == "not supported",
                "qualified performance requires an explicit unchanged cache policy",
            )
    return {
        "evidence": evidence(path),
        "raw": raw,
        "build_receipt": {"evidence": evidence(build_path), "raw": build},
    }


def load_pair(path):
    raw = read(path)
    require(
        raw["schema"] == "vakint.ordered-rung-parity.v1"
        and raw["exact_retained_equality"] is True,
        "missing successful exact parity receipt",
    )
    reports = {}
    for name, strategy in (("dfs", "DepthFirst"), ("ordered", "Ordered")):
        report_path = referenced(raw["reports"][name], path.parent)
        report = validate_report(report_path, strategy)
        require(
            report["input_sha256"] == raw["input_sha256"]
            and report["retained"]["sha256"] == raw["retained"][name],
            "parity receipt does not bind the input/output",
        )
        reports[name] = {"evidence": evidence(report_path), "raw": report}
    dfs, ordered = (reports[name]["raw"] for name in ("dfs", "ordered"))
    for key in ("worker_policy", "cache_policy"):
        require(dfs.get(key) == ordered.get(key), f"paired {key} differs")
    for key in (
        "loop_count",
        "projection",
        "input_sha256",
        "evaluated_input_sha256",
        "dimension_binding",
        "monitoring_enabled",
        "data",
        "normalization",
        "producer",
    ):
        # Same-process/fresh-process mode labels can differ, but the complete
        # producer document must agree within each actual pair.
        require(dfs[key] == ordered[key], f"paired {key} differs")
    has_numerical = dfs["numerical_reference"] is not None
    require(
        has_numerical == (ordered["numerical_reference"] is not None),
        "only one engine produced a numerical reference",
    )
    require(
        raw["numerical_enclosures_agree"] is (True if has_numerical else None),
        "numerical parity scope differs from actual references",
    )
    if has_numerical:
        for key in ("point", "epsilon_range", "coupling_definitions"):
            require(
                dfs["numerical_reference"][key] == ordered["numerical_reference"][key],
                f"paired numerical {key} differs",
            )
        require(
            dfs["numerical_reference"]["selection"]["requested_decimal_precision"]
            == ordered["numerical_reference"]["selection"][
                "requested_decimal_precision"
            ],
            "paired numerical requested precision differs",
        )
        require(
            source_cap(dfs["numerical_reference"])
            == source_cap(ordered["numerical_reference"]),
            "paired numerical effective source cap differs",
        )
        require(
            [row["epsilon_power"] for row in dfs["numerical_reference"]["coefficients"]]
            == [
                row["epsilon_power"]
                for row in ordered["numerical_reference"]["coefficients"]
            ],
            "paired numerical Laurent orders differ",
        )
    return {"parity": {"evidence": evidence(path), "raw": raw}, "reports": reports}


def load_historical(path):
    """Keep the accepted old q0 cohort distinct from new matched measurements."""
    raw = read(path)
    require(
        raw["schema"] == "vakint.ordered-rung-historical-parity.v1"
        and raw["exact_retained_equality"] is True,
        "missing exact historical-reference parity",
    )
    acceptance_path = referenced(raw["historical"]["acceptance"], path.parent)
    acceptance = read(acceptance_path)
    require(
        acceptance["schema"]
        == "vakint.physical-q0.cold-completion-parent-acceptance.v1"
        and acceptance["status"] == "ACCEPTED_FUNCTIONAL_COLD_COMPLETE_CACHE_CLOSED",
        "historical baseline was not accepted complete",
    )
    baseline_path = referenced(raw["historical"]["report"], path.parent)
    baseline = read(baseline_path)
    require(
        baseline["schema"] == "vakint.qcd-ladder.benchmark.v1"
        and baseline["reduction_completed"] is True
        and baseline["retained_result_available"] is True
        and baseline["catalogue_closure"]["closed"] is True
        and baseline["catalogue_closure"]["unknown_count"] == 0
        and baseline["catalogue_closure"] == acceptance["catalogue_closure"],
        "historical reduction is incomplete",
    )
    saved = acceptance["output_artifacts_verified"]["benchmark.json"]
    require(
        digest(baseline_path) == saved["sha256"]
        and baseline_path.stat().st_size == saved["bytes"],
        "historical report differs from its acceptance",
    )
    result_path = referenced(raw["historical"]["retained"], path.parent)
    saved = acceptance["output_artifacts_verified"]["result.symbolica"]
    require(
        digest(result_path) == saved["sha256"] == baseline["result_sha256"]
        and result_path.stat().st_size == saved["bytes"] == baseline["result_bytes"],
        "historical retained result differs from its acceptance",
    )
    observed = seconds(baseline["reduction_seconds"])
    require(
        observed == seconds(acceptance["collector_reduction_seconds"]),
        "historical accepted timing differs",
    )
    ordered_path = referenced(raw["ordered_report"], path.parent)
    ordered = validate_report(ordered_path, "Ordered")
    require(
        ordered["loop_count"] == 5
        and ordered["projection"] == baseline["projection"] == "zero-momentum-metric"
        and ordered["input_sha256"] == baseline["input_sha256"] == raw["input_sha256"],
        "historical physical input/projection/normalization differs",
    )
    definition_path = referenced(raw["historical"]["definition"], path.parent)
    definition = read(definition_path)
    require(
        definition["schema"] == "vakint.ordered-rungs.legacy-definition-identity.v1"
        and definition["status"] == "PASS"
        and definition["normalization"] == "pySecDec"
        and definition["historical_scalar_input"]["sha256"] == raw["input_sha256"]
        and definition["historical_result_sha256"] == saved["sha256"],
        "historical normalization definition is not bound to the old result/input",
    )
    require(
        "normalization_conversion" in raw,
        "missing exact normalization conversion proof",
    )
    conversion = raw["normalization_conversion"]
    require(
        conversion["schema"] == "vakint.ordered-rung-normalization-conversion.v1"
        and conversion["source_normalization"] == "pySecDec"
        and conversion["target_normalization"] == "MSbar"
        and conversion["api"]
        == "VakintSettings::get_integral_normalization_factor_atom"
        and conversion["ratio_identity"] == "N_target/N_source"
        and conversion["loop_count"] == ordered["loop_count"]
        and conversion["dimension_binding"] == ordered["dimension_binding"]
        and all(
            isinstance(conversion[key], str) and conversion[key]
            for key in ("epsilon_symbol", "mu_r_sq_symbol")
        )
        and conversion["source_retained_sha256"] == saved["sha256"]
        and conversion["target_report_sha256"] == digest(ordered_path)
        and conversion["exact_conversion_verified"] is True
        and conversion["applied_before_laurent_truncation"] is True,
        "missing bound exact pySecDec-to-MSbar conversion proof",
    )
    factor = conversion["factor"]
    require(
        isinstance(factor["canonical_expression"], str)
        and factor["canonical_expression"]
        and hashlib.sha256(factor["canonical_expression"].encode()).hexdigest()
        == factor["sha256"],
        "normalization factor identity differs",
    )
    converted_path = referenced(conversion["converted_retained"], path.parent)
    require(
        converted_path != result_path
        and converted_path.stat().st_size == conversion["converted_retained"]["bytes"],
        "converted historical payload was not separately preserved",
    )
    require(
        raw["retained"]["historical"] == saved["sha256"]
        and raw["retained"]["historical_converted"]
        == conversion["converted_retained"]["sha256"]
        and raw["retained"]["ordered"] == ordered["retained"]["sha256"],
        "historical parity does not bind both retained payloads",
    )
    return {
        "parity": {"evidence": evidence(path), "raw": raw},
        "historical": {
            "acceptance": {"evidence": evidence(acceptance_path), "raw": acceptance},
            "report": {"evidence": evidence(baseline_path), "raw": baseline},
            "retained": evidence(result_path),
            "definition": {"evidence": evidence(definition_path), "raw": definition},
            "normalization_conversion": conversion,
            "public_retained_seconds": observed,
            "reduction_math_seconds": None,
            "master_count": acceptance["catalogue_closure"]["master_count"],
            "scope": acceptance["measurement_scope"],
        },
        "ordered": {
            "report": {"evidence": evidence(ordered_path), "raw": ordered},
            "public_retained_seconds": ordered["retained_reduction_seconds"],
            "reduction_math_seconds": reduction_math(ordered),
        },
        "exact_retained_equality": True,
        "performance_eligible": False,
        "ordered_over_dfs": None,
        "scope": "one historical capture-enabled pySecDec DFS observation and one new MSbar ordered observation; exact converted historical result compared; different cohorts; no causal or matched speed ratio",
    }


def append(base, run_id, receipts, protocol=None, historical_receipts=()):
    require(run_id and run_id not in base.get("ordered_rung_runs", {}), "run ID exists")
    pairs = [load_pair(Path(path).resolve()) for path in receipts]
    historical = [load_historical(Path(path).resolve()) for path in historical_receipts]
    require(pairs or historical, "at least one checked comparison is required")
    require(
        not protocol or pairs, "a fresh-process protocol cannot qualify old cohorts"
    )
    hashes = [item["parity"]["evidence"]["sha256"] for item in pairs + historical]
    require(len(hashes) == len(set(hashes)), "duplicate parity receipt")
    qualified = validate_protocol(Path(protocol).resolve(), pairs) if protocol else None
    grouped = {}
    for pair in pairs:
        raw = pair["reports"]["dfs"]["raw"]
        key = (raw["loop_count"], raw["projection"], raw["input_sha256"])
        grouped.setdefault(key, []).append(pair)
    rows = []
    for (loops, projection, identity), samples in sorted(grouped.items()):
        first = samples[0]["reports"]["dfs"]["raw"]
        for sample in samples[1:]:
            current = sample["reports"]["dfs"]["raw"]
            for name in ("producer", "data"):
                require(current[name] == first[name], f"sample cohort {name} differs")
            for name in (
                "worker_policy",
                "cache_policy",
                "dfs_request_memoization",
                "monitoring_enabled",
            ):
                require(
                    current.get(name) == first.get(name),
                    f"sample cohort {name} differs",
                )
            require(
                sample["parity"]["raw"]["numerical_enclosures_agree"]
                is samples[0]["parity"]["raw"]["numerical_enclosures_agree"],
                "mixed numerical-validation scope in one cohort",
            )
        metrics = {}
        for name in ("dfs", "ordered"):
            values = [sample["reports"][name]["raw"] for sample in samples]
            math_values = [reduction_math(row) for row in values]
            math_complete = all(value is not None for value in math_values)
            metrics[name] = {
                "reduction_math_seconds": summary(math_values)
                if math_complete
                else None,
                "math_gross_seconds": summary(
                    [row["reduction_math_timing"]["gross_seconds"] for row in values]
                )
                if math_complete
                else None,
                "measured_owner_preparation_seconds": summary(
                    [
                        row["reduction_math_timing"][
                            "measured_owner_preparation_seconds"
                        ]
                        for row in values
                    ]
                )
                if math_complete
                else None,
                "math_timing_samples": [
                    row.get("reduction_math_timing") for row in values
                ],
                "public_retained_seconds": summary(
                    [row["retained_reduction_seconds"] for row in values]
                ),
                "loading_seconds": summary([row["loading_seconds"] for row in values]),
                "collector_seconds": summary(
                    [row["completed_seconds"] for row in values]
                ),
                "pure_ibp_seconds": None,
                "routing_lowering_seconds": None,
                "terminal_substitution_seconds": None,
                "peak_rss_bytes": [
                    row["memory_after_retained_reduction"].get("process_peak_rss_bytes")
                    for row in values
                ],
                "ordered_progress": [row.get("ordered_progress") for row in values],
            }
        baseline = metrics["dfs"]["reduction_math_seconds"]
        ordered = metrics["ordered"]["reduction_math_seconds"]
        rows.append(
            {
                "loop_count": loops,
                "projection": projection,
                "input_sha256": identity,
                "sample_count": len(samples),
                "data": first["data"],
                "producer": first["producer"],
                "dfs_request_memoization": first["dfs_request_memoization"],
                "cache_policy": first.get("cache_policy"),
                "worker_policy": first.get("worker_policy"),
                "monitoring_enabled": first["monitoring_enabled"],
                "headline_timing_metric": "reduction_math_seconds",
                "math_scope_id": MATH_SCOPE,
                "metrics": metrics,
                "performance_eligible": qualified is not None,
                "ordered_over_dfs": (
                    ordered["median_seconds"] / baseline["median_seconds"]
                    if qualified is not None
                    and baseline is not None
                    and ordered is not None
                    and baseline["median_seconds"] > 0
                    else None
                ),
                "exact_retained_equality": True,
                "numerical_enclosures_agree": samples[0]["parity"]["raw"][
                    "numerical_enclosures_agree"
                ],
                "samples": samples,
            }
        )
    result = copy.deepcopy(base)
    result.setdefault("ordered_rung_runs", {})[run_id] = {
        "schema": "vakint.ordered-rung-report.v1",
        "persistent_caches": False,
        "dfs_request_memoization": "recorded per row from the actual DFS report; persistent caches remain off",
        "headline_timing_scope": "net checked non-overlapping reduction-math spans: request-local reducer/dispatch state, IBP/routing and exact weighted terminal-coefficient assembly/grouping; immutable artifact/program import and measured lazy owner preparation excluded",
        "timing_scope": "public retained evaluation includes routing/lowering/finalization; loading and later numerical reference excluded",
        "collector_scope": "prepared-input generation excluded; numerical reference included when present",
        "protocol": qualified,
        "summary": rows,
        "historical_comparisons": historical,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--parity", type=Path, action="append", default=[])
    parser.add_argument("--historical-parity", type=Path, action="append", default=[])
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--select-current", action="store_true")
    args = parser.parse_args()
    require(not args.output.exists(), "write a new output; never overwrite history")
    result = append(
        read(args.base), args.id, args.parity, args.protocol, args.historical_parity
    )
    if args.select_current:
        result["selected_ordered_rung_run"] = args.id
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
