"""Audit archived reports and native logs only; never load or evaluate GammaLoop."""

import gzip
import json
import re
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 650
root = Path(__file__).resolve().parent


def read(name):
    path = root / name
    if not path.exists():
        return gzip.decompress(
            path.with_suffix(path.suffix + ".gz").read_bytes()
        ).decode()
    return path.read_text()


cuts, joint = [
    json.loads(read(f"minimal_point/{mode}.json")) for mode in ("cuts", "cuts_joint")
]
left, right = cuts["reference_point"], joint["reference_point"]
assert left["halton_index"] == right["halton_index"] == 3363
assert left["cube_binary64_bits"] == right["cube_binary64_bits"]
assert left["cube_arb"] == right["cube_arb"]
ordinary, extended = [
    value["generator_foreign_inverse_rows"] for value in (left, right)
]
assert len(ordinary) == 6 and len(extended) == 7
for a, b in zip(ordinary, extended):
    assert a == b
    assert a["foreign_inverse"]["status"] == "ok"
failed = extended[-1]
assert failed["generator_name"] == "direct_joint_HZ"
assert failed["foreign_name"] == "lu_cut_1"
assert failed["foreign_inverse"]["status"] == "error"
assert (
    "joint: normalized ordinary fallback; disk policy declined"
    in failed["forward"]["diagnostics"]
)
assert len(failed["forward"]["point"]) == 12
error = failed["foreign_inverse"]["error"]
assert left["summed_reference"]["status"] == "ok"
for result in (right["summed_reference"], right["selected_reference"]["result"]):
    assert result["status"] == "error" and result["error"].endswith(error)
for field in ("cuts", "production_orientation_keys", "threshold_counterterm_metadata"):
    assert cuts["inventory"][field] == joint["inventory"][field]
assert len(joint["inventory"]["cuts"]) == 6
assert len(joint["inventory"]["production_orientation_keys"]) == 936
states = []
for folder in ("failed_preflight", "minimal_point"):
    before, after = [
        json.loads(read(f"{folder}/state_{time}.json")) for time in ("before", "after")
    ]
    assert before == after and len(before) == 35
    states.append(before)
assert states[0] == states[1]

log = read("fixture/comparison.log")
results = {}
for name in ("old_strict", "old_checked", "ray_strict", "ray_checked"):
    line = re.search(rf"^\s*{name}\s+(.*)$", log, re.MULTILINE).group(1).strip()
    results[name] = {
        "accepted": line.startswith("Ok("),
        "radius": re.search(r"solution: F\(VarFloat \{ float: ([^ ]+)", line).group(1),
        "derivative": re.search(
            r"derivative_at_solution: F\(VarFloat \{ float: ([^ ]+)", line
        ).group(1),
        "residual": re.search(
            r"error_of_function: F\(VarFloat \{ float: ([^ ]+)", line
        ).group(1),
        "iterations": int(re.search(r"num_iterations_used: (\d+)", line).group(1)),
    }
assert results["old_strict"] == results["old_checked"]
assert results["ray_strict"] == results["ray_checked"]
assert (
    not results["old_strict"]["accepted"] and results["old_strict"]["iterations"] == 17
)
assert results["ray_strict"]["accepted"] and results["ray_strict"]["iterations"] == 8
low, high = [
    Decimal(
        re.search(rf"^\s*original_residual_{side}\s+(\S+)", log, re.MULTILINE).group(1)
    )
    for side in ("lower", "upper")
]
budget = Decimal(64) * Decimal(2) ** -999 * Decimal(529)
assert 0 < low <= high <= budget
assert abs(Decimal(results["old_strict"]["residual"])) > budget
assert abs(Decimal(results["ray_strict"]["residual"])) <= budget
summary = {
    "scope": "Old build4 failure, exact callback comparison and repaired build5 point replay; full seven-mode acceptance is separate.",
    "halton_index": 3363,
    "forward_points_retained": 13,
    "foreign_cut1_inverse_successes": 12,
    "foreign_cut1_inverse_failures": 1,
    "ordinary_rows_identical_between_catalogues": 6,
    "failing_generator_policy": "normalized ordinary fallback; disk declined",
    "original_selected_and_summed_failure_reproduced": True,
    "state_hashes_unchanged": 35,
    "native_precision_bits": 1000,
    "root_results": results,
    "original_equation_certificate": {
        "lower": str(low),
        "upper": str(high),
        "width": str(high - low),
        "unchanged_budget_64_epsilon_529": str(budget),
        "upper_over_budget": str(high / budget),
        "interval_contains_zero": False,
        "meaning": "Approximate root satisfies the original routed equation budget; no exact-zero claim.",
    },
}

repaired = [
    json.loads(read(f"replay_build5/{mode}.json")) for mode in ("cuts", "cuts_joint")
]
for old, new in zip((cuts, joint), repaired):
    assert new["error"] is None
    assert old["correctness_settings"] == new["correctness_settings"]
    for field in (
        "cuts",
        "production_orientation_keys",
        "threshold_counterterm_metadata",
    ):
        assert old["inventory"][field] == new["inventory"][field]
    before, after = [report["reference_point"] for report in (old, new)]
    for field in (
        "halton_index",
        "cube",
        "cube_binary64_bits",
        "cube_arb",
        "summed_sample",
    ):
        assert before[field] == after[field]
    for a, b in zip(
        before["generator_foreign_inverse_rows"],
        after["generator_foreign_inverse_rows"],
        strict=True,
    ):
        assert a["forward"] == b["forward"]
        assert b["foreign_inverse"]["status"] == "ok"
    assert after["summed_reference"]["status"] == "ok"
selected = repaired[1]["reference_point"]["selected_reference"]
assert selected["sample"] == right["selected_reference"]["sample"]
assert selected["result"]["status"] == "ok"
assert selected["result"]["report"]["finite_points"] == 1
for time in ("before", "after"):
    assert json.loads(read(f"replay_build5/state_{time}.json")) == states[0]
build = json.loads(read("replay_build5/driver_build.json"))
old_build = json.loads(read("minimal_point/driver_build.json"))
assert build["source_sha256"] == old_build["source_sha256"]
assert build["returncode"] == 0
run = json.loads(read("replay_build5/run.json"))
assert run["returncode"] == 0
summary["build5_replay"] = {
    "forward_rows_exactly_unchanged": 13,
    "foreign_cut1_inverse_successes": 13,
    "summed_reference_successes": 2,
    "selected_joint_reference_successes": 1,
    "selected_joint_reported_normalization": selected["result"]["report"][
        "normalization"
    ],
    "same_driver_source_sha256": build["source_sha256"],
    "state_hashes_unchanged": 35,
    "elapsed_seconds": run["elapsed_seconds"],
    "scope": "Pointwise pipeline repair at the identical far-tail ordinary-fallback point; no one-point normalization, physical or full-preflight acceptance claim.",
}
print(json.dumps(summary, indent=2))
