#!/usr/bin/env python3
"""Audit current native GL638 reports; never load state or evaluate physics.

Usage: python audit_physical.py MANIFEST OUTPUT_DIRECTORY
Keeps the archived audit's Decimal350 total/cut checks, without comparing to
historical values or assuming bare and partition-weighted results are equal.
"""
from decimal import Decimal, localcontext
from pathlib import Path
import json
import math
import sys
import tomllib

D = Decimal
PHASES = ("re", "im")
BASELINE = "optimized_lmb"
CANDIDATE = "joint_hz_plus_lmb"


def read(path):
    return json.loads(path.read_text())


def number(value):
    if not isinstance(value, str):
        raise ValueError("native values must remain decimal strings")
    result = D(value)
    if not result.is_finite():
        raise ValueError(f"nonfinite native value: {value}")
    return result


def vector(value):
    return tuple(number(value[p]) for p in PHASES)


def norm(value):
    return sum(x * x for x in value).sqrt()


def metric(actual, expected, tolerance, minimum_scale=D(0)):
    delta = tuple(a - b for a, b in zip(actual, expected))
    scale = max(norm(actual), norm(expected), minimum_scale)
    error = norm(delta)
    limit = tolerance * scale
    components = {}
    for p, a, b, d in zip(PHASES, actual, expected, delta):
        component_scale = max(abs(a), abs(b))
        own_pass = abs(d) <= tolerance * component_scale
        components[p] = {
            "actual": str(a), "expected": str(b), "absolute_error": str(abs(d)),
            "relative_error": str(abs(d) / component_scale) if component_scale else "0",
            "relative_scale": str(component_scale), "own_relative_budget_passed": own_pass,
            "absolute_norm_budget": str(limit),
        }
    return {
        "passed": error <= limit, "tolerance": str(tolerance),
        "complex_error": str(error), "scale": str(scale), "absolute_budget": str(limit),
        "relative_complex_error": str(error / scale) if scale else "0",
        "components": components,
    }


def point_key(case):
    tokens = case["point_tokens"]
    values = [float(x) for x in tokens]
    if len(values) != 12 or not all(math.isfinite(x) for x in values):
        raise ValueError("expected twelve finite original binary64 coordinates")
    return tuple(x.hex() for x in values)


def events(evaluation):
    result = {}
    for event in evaluation["events"]:
        cut_id = event["cut_info"]["cut_id"]
        if cut_id in result:
            raise ValueError(f"duplicate physical cut {cut_id}")
        result[cut_id] = event
    if sorted(result) != list(range(6)):
        raise ValueError(f"expected six physical cuts, got {sorted(result)}")
    return result


def event_identity(event):
    # A selected raw row legitimately adds channel provenance. Everything else,
    # including graph/orientation/cut identity, must still match the bare row.
    info = dict(event["cut_info"])
    info.pop("sampling_channel_id", None)
    info.pop("sampling_channel_edge_ids", None)
    return event["group"], info


def compare(left, right, tolerance, weight=D(1)):
    actual = vector(left["evaluation"]["integrand_result"])
    expected = tuple(x * weight for x in vector(right["evaluation"]["integrand_result"]))
    total = metric(actual, expected, tolerance)
    left_events, right_events = events(left["evaluation"]), events(right["evaluation"])
    cuts = []
    for cut_id in range(6):
        a, b = left_events[cut_id], right_events[cut_id]
        same = event_identity(a) == event_identity(b)
        cut = metric(vector(a["weight"]), tuple(x * weight for x in vector(b["weight"])),
                     tolerance, max(norm(actual), norm(expected)))
        cuts.append({"cut_id": cut_id, "identity_matches": same, **cut})
    return {"passed": total["passed"] and all(c["passed"] and c["identity_matches"] for c in cuts),
            "weight_applied_to_expected": str(weight), "total": total, "cuts": cuts}


def audit(manifest_path, directory):
    manifest = read(manifest_path)
    result = {
        "passed": False, "decimal_precision": 350, "source_manifest": str(manifest_path),
        "criteria": {
            "physical_scope": "current candidate only; all936 original orientations, six cuts, CT on",
            "stability_scope": "production cards require check_on_norm=true (omitted saved field uses the Rust serde default true)",
            "tolerance": "sum of both final-precision requirements, min(Re,Im) per row",
            "total": "complex norm(error) <= tolerance * max(norm(actual),norm(expected))",
            "small_component": "same absolute complex-norm budget; own relative misses retained as diagnostics",
            "whole_zero": "exact equality required; no physical-unit absolute floor",
            "cut": "same complex criterion with scale at least either complete total norm",
            "event_sum": "one row precision budget; scale is max(sum norm,total norm)",
            "selected_raw": "forced-Arb selected = forced-Arb bare * independent canonical Arb w; no J or grid factor",
            "source": "exact original binary64 coordinate identity, not decimal redefinition",
            "limitations": "acceptance tolerances are not rigorous error bounds; no CT-star or historical-value claim",
        },
        "settings": {}, "retained_rows": [], "checks": [], "component_relative_misses": [],
    }
    checks = result["checks"]

    def check(kind, label, run):
        try:
            detail = run()
            record = {"kind": kind, "label": label, **detail}
        except (KeyError, TypeError, ValueError, ArithmeticError) as error:
            record = {"kind": kind, "label": label, "passed": False, "error": str(error)}
        checks.append(record)
        return record

    def require(condition, message):
        if not condition:
            raise ValueError(message)
        return {"passed": True}

    indexed = {}
    rows_with_budgets = []
    reports = {}
    for mode in manifest["proposal_order"]:
        report_path = directory / f"{mode}.json"
        if not report_path.exists():
            if mode in (BASELINE, CANDIDATE):
                check("inventory", mode, lambda: require(False, "required current mode report missing"))
            continue
        report = reports[mode] = read(report_path)
        with (directory / f"{mode}.settings.toml").open("rb") as stream:
            settings = tomllib.load(stream)
        levels = settings["stability"]["levels"]
        budgets = {}
        for level in levels:
            budget = min(D(str(level["required_precision_for_re"])), D(str(level["required_precision_for_im"])))
            if not budget.is_finite() or budget <= 0:
                raise ValueError("physics precision requirements must be positive and finite")
            precision = level["precision"]
            budgets[precision] = min(budgets.get(precision, budget), budget)
        result["settings"][mode] = {"check_on_norm": settings["stability"].get("check_on_norm", True),
                                   "precision_budgets": {p: str(v) for p, v in budgets.items()}}
        inventory = report["inventory"]
        check("inventory", mode, lambda: require(
            settings["stability"].get("check_on_norm", True) is True
            and len(inventory["production_orientation_keys"]) == 936
            and len(set(inventory["production_orientation_keys"])) == 936
            and inventory["exposed_orientation_selectors"] == 1
            and sorted(x["id"] for x in inventory["cuts"]) == list(range(6)),
            "wrong norm/orientation/cut scope"))
        check("coverage", mode, lambda: require(bool(report.get("replays")), "no retained replay rows"))
        for row in report.get("replays", []):
            result["retained_rows"].append({"mode": mode, "row": row})
            case, evaluation = row["case"], row["evaluation"]
            label = f"{mode}:{case['name']}"

            def validate_row():
                require(evaluation.get("valid") is True, evaluation.get("error", "invalid native evaluation"))
                metadata = evaluation["evaluation_metadata"]
                require(metadata["is_nan"] is False, "native metadata marks invalid value")
                stability = metadata["stability_results"]
                require(bool(stability), "missing final stability status")
                status = stability[-1]["status"]
                require(not (isinstance(status, dict) and "Unstable" in status)
                        and not (isinstance(status, str) and status.startswith("Unstable")),
                        "final stability status is unstable")
                require(case["mode"] == mode and case["threshold_ct_enabled"] is True
                        and case["momentum_space"] is True, "not the declared mode's original CT-on raw point")
                require(number(evaluation["integrator_weight"]) == 1
                        and evaluation["parameterization_jacobian"] is None, "raw input has unexpected outer factor")
                if case["forced_arb"]:
                    require(evaluation["precision"] == "Arb", "forced-Arb row did not return Arb")
                precision = evaluation["precision"]
                key = (mode, point_key(case), case.get("channel_name"), case["forced_arb"])
                require(key not in indexed, "duplicate mode/point/channel/precision-request row")
                total = vector(evaluation["integrand_result"])
                cut_events = events(evaluation)
                channel = case.get("channel_name")
                matches = [x for x in inventory["channels"] if x["label"] == channel]
                require(channel is None or len(matches) == 1, "selected canonical channel label missing/ambiguous")
                expected_channel = matches[0]["id"] if channel is not None else None
                require(all(e["cut_info"]["graph_id"] == 0
                            and e["cut_info"].get("sampling_channel_id") == expected_channel
                            for e in cut_events.values()), "physical event graph/channel provenance differs from request")
                event_sum = tuple(sum(vector(e["weight"])[i] for e in cut_events.values()) for i in range(2))
                summed = metric(event_sum, total, budgets[precision])
                indexed[key] = (row, budgets[precision])
                rows_with_budgets.append((key, row, budgets[precision]))
                return {"passed": summed["passed"], "event_sum": summed, "precision": precision}

            check("row_validity_and_event_sum", label, validate_row)

    for key, row, budget in rows_with_budgets:
        mode, point, channel, forced = key
        label = f"{mode}:{row['case']['name']}"
        if forced:
            check("coverage", label + ":ordinary_partner", lambda: require(
                (mode, point, channel, False) in indexed, "forced-Arb row lacks ordinary-stack partner"))
        if mode == BASELINE and channel is None:
            check("coverage", label + ":candidate_partner", lambda: require(
                (CANDIDATE, point, None, forced) in indexed, "baseline point lacks current joint bare partner"))
        if not forced:
            def precision_pair():
                other, other_budget = indexed[(mode, point, channel, True)]
                return compare(row, other, budget + other_budget)
            check("ordinary_vs_forced_arb", label, precision_pair)
        if channel is None and mode != BASELINE:
            def mode_pair():
                other, other_budget = indexed[(BASELINE, point, None, forced)]
                return compare(row, other, budget + other_budget)
            check("current_bare_vs_baseline", label, mode_pair)
        if channel is not None and forced:
            def partition_pair():
                bare, bare_budget = indexed[(mode, point, None, True)]
                oracle = row["canonical_raw_partition"]
                require(oracle["status"] == "inside_support" and oracle["precision"] == "Arb1000",
                        "selected raw point has no successful canonical Arb oracle")
                matches = [x for x in reports[mode]["inventory"]["channels"] if x["label"] == channel]
                require(len(matches) == 1 and oracle["channel_id"] == matches[0]["id"],
                        "oracle channel does not match canonical label")
                weight = number(oracle["weight"])
                require(0 < weight <= 1, "canonical selected partition must be finite in (0,1]")
                compared = compare(row, bare, budget + bare_budget, weight)
                return {**compared, "canonical_raw_partition": oracle}
            check("forced_arb_selected_equals_bare_times_w", label, partition_pair)

    required = {"ordinary_vs_forced_arb", "current_bare_vs_baseline", "forced_arb_selected_equals_bare_times_w"}
    for kind in sorted(required):
        check("coverage", kind, lambda kind=kind: require(any(c["kind"] == kind for c in checks), "required comparison absent"))
    for comparison in checks:
        values = [("total", comparison.get("total")), ("event_sum", comparison.get("event_sum"))]
        values.extend((f"cut{cut['cut_id']}", cut) for cut in comparison.get("cuts", []))
        for scope, value in values:
            if value:
                for phase, component in value["components"].items():
                    if not component["own_relative_budget_passed"]:
                        result["component_relative_misses"].append({"kind": comparison["kind"],
                            "label": comparison["label"], "scope": scope, "phase": phase, **component})
    result["passed"] = bool(checks) and all(c["passed"] for c in checks)
    result["counts"] = {"retained_rows": len(result["retained_rows"]), "checks": len(checks),
                        "failed_checks": sum(not c["passed"] for c in checks),
                        "component_relative_misses": len(result["component_relative_misses"])}
    return result


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: audit_physical.py MANIFEST OUTPUT_DIRECTORY")
    manifest_path, directory = Path(sys.argv[1]), Path(sys.argv[2])
    try:
        with localcontext() as context:
            context.prec = 350
            result = audit(manifest_path, directory)
    except (OSError, KeyError, TypeError, ValueError, ArithmeticError) as error:
        result = {"passed": False, "error": str(error), "scope": "report/schema failure; no physical run performed"}
    (directory / "physics_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed": result["passed"], "counts": result.get("counts"), "error": result.get("error")}))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
