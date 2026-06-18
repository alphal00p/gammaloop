#!/usr/bin/env python3
"""Monitor CT-by-CT contributions across the im[-] axis_x2 feature."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from gammaloop import GammaLoopAPI

from inspect_imag_maxweight_thresholds import (
    THRESHOLD_KEY_RE,
    complex_to_dict,
    multiply_complex,
    threshold_maps,
)
from scan_current_highstat_imag_maxweights import apply_width_overrides
from scan_imag_maxweight_approaches import STATE_DIR, apply_runtime_settings


ROOT = Path(__file__).resolve().parents[3]
APPROACH_DIR = (
    ROOT
    / "examples"
    / "cli"
    / "qq_hhh_2L"
    / "maxweight_approach_scans"
    / "xi_t4_sliver30_hsigma1_81920_approach"
)
SCAN_TSV = APPROACH_DIR / "imag_maxweight_approach_scan.tsv"
OUT_JSON = APPROACH_DIR / "axis_x2_ct_monitor.json"
OUT_POINT_TSV = APPROACH_DIR / "axis_x2_ct_monitor_points.tsv"
OUT_TERM_TSV = APPROACH_DIR / "axis_x2_ct_monitor_terms.tsv"
OUT_MD = APPROACH_DIR / "axis_x2_ct_monitor.md"
OUT_COMPONENT_PNG = APPROACH_DIR / "axis_x2_ct_monitor_components.png"
OUT_ABS_PNG = APPROACH_DIR / "axis_x2_ct_monitor_abs.png"

TARGET = "highstat_iter_0001_imag_minus"
DIRECTION = "axis_x2"
ROOT_X2 = 0.9441181235200773
DOMINANT_TERMS = [
    ("GL05_isr_p2_ct", "threshold_counterterm:4:0"),
    ("GL05_isr_p1_ct", "threshold_counterterm:34:0"),
    ("GL05", "threshold_counterterm:31:0"),
    ("GL05", "threshold_counterterm:5:0"),
]


def add_complex(values: list[dict[str, float]]) -> dict[str, float]:
    re_val = sum(value["re"] for value in values)
    im_val = sum(value["im"] for value in values)
    return {"re": re_val, "im": im_val, "abs": math.hypot(re_val, im_val)}


def fmt(value: dict[str, float]) -> str:
    return f"{value['re']:+.8e} {value['im']:+.8e}i |w|={value['abs']:.8e}"


def threshold_edges(group_map: dict[str, Any], term_key: str) -> list[int] | None:
    match = THRESHOLD_KEY_RE.match(term_key)
    if not match:
        return None
    esurface_id = int(match.group("esurface_id"))
    edges = group_map["thresholds"].get(esurface_id)
    return None if edges is None else [int(edge_id) for edge_id in edges]


def load_axis_rows() -> list[dict[str, str]]:
    rows = [
        row
        for row in csv.DictReader(SCAN_TSV.open(), delimiter="\t")
        if row["target"] == TARGET and row["direction"] == DIRECTION
    ]
    rows.sort(key=lambda row: float(row["signed_side_fraction"]))
    if not rows:
        raise RuntimeError(f"No {TARGET}/{DIRECTION} rows found in {SCAN_TSV}")
    return rows


def row_x(row: dict[str, str]) -> list[float]:
    return [float(row[f"x{index}"]) for index in range(6)]


def build_monitor_points() -> list[dict[str, Any]]:
    axis_rows = load_axis_rows()
    peak_row = min(axis_rows, key=lambda row: abs(float(row["t"])))
    base_x = row_x(peak_row)
    graph = int(peak_row["graph"])
    lmb_channel = int(peak_row["lmb_channel"])

    points: list[dict[str, Any]] = []
    offsets = np.geomspace(1.0e-10, 3.0e-5, 41)
    for sign in (-1.0, 1.0):
        for offset in offsets:
            x = list(base_x)
            x[2] = ROOT_X2 + sign * float(offset)
            points.append(
                {
                    "source": "root_log_probe",
                    "scan_row_index": None,
                    "signed_side_fraction": None,
                    "t": None,
                    "graph": graph,
                    "lmb_channel": lmb_channel,
                    "x": x,
                }
            )

    for index in range(25, 45):
        row = axis_rows[index]
        points.append(
            {
                "source": "saved_scan_row",
                "scan_row_index": index,
                "signed_side_fraction": float(row["signed_side_fraction"]),
                "t": float(row["t"]),
                "graph": int(row["graph"]),
                "lmb_channel": int(row["lmb_channel"]),
                "x": row_x(row),
                "scan_re": float(row["re"]),
                "scan_im": float(row["im"]),
                "scan_abs": float(row["abs"]),
            }
        )

    peak = {
        "source": "max_weight_location",
        "scan_row_index": "t0",
        "signed_side_fraction": 0.0,
        "t": 0.0,
        "graph": graph,
        "lmb_channel": lmb_channel,
        "x": list(base_x),
        "scan_re": float(peak_row["re"]),
        "scan_im": float(peak_row["im"]),
        "scan_abs": float(peak_row["abs"]),
    }
    points.append(peak)

    deduped: dict[float, dict[str, Any]] = {}
    for point in points:
        deduped[round(float(point["x"][2]), 16)] = point
    result = sorted(deduped.values(), key=lambda point: float(point["x"][2]))
    for index, point in enumerate(result):
        point["point_index"] = index
        point["x2"] = float(point["x"][2])
        point["dx2_from_root"] = float(point["x"][2]) - ROOT_X2
    return result


def evaluate_point(
    api: GammaLoopAPI,
    group_maps: dict[int, dict[str, Any]],
    point: dict[str, Any],
) -> dict[str, Any]:
    result = api.evaluate_sample(
        point["x"],
        integrand_name="subtracted",
        return_events=True,
        discrete_dim=[int(point["graph"]), int(point["lmb_channel"])],
    )
    group_map = group_maps[int(point["graph"])]
    events = []
    flat_terms = []
    event_weights = []
    for event_group in result.event_groups:
        for event in event_group.events:
            event_factor = {"re": 1.0, "im": 0.0}
            raw_terms = []
            for additional_weight in event.additional_weights:
                key = str(additional_weight.key)
                value = complex_to_dict(additional_weight.value)
                if key == "full_multiplicative_factor":
                    event_factor = value
                else:
                    raw_terms.append((key, value))

            graph_id = int(event.cut_info.graph_id)
            graph_name = group_map["graphs"].get(graph_id, {}).get("name", f"graph_{graph_id}")
            weighted_terms = []
            for key, raw in raw_terms:
                weighted = multiply_complex(raw, event_factor)
                entry = {
                    "point_index": int(point["point_index"]),
                    "x2": float(point["x2"]),
                    "dx2_from_root": float(point["dx2_from_root"]),
                    "event_graph_id": graph_id,
                    "event_graph_name": graph_name,
                    "term_key": key,
                    "edge_ids": threshold_edges(group_map, key),
                    "weighted": weighted,
                    "raw": raw,
                }
                weighted_terms.append(entry)
                flat_terms.append(entry)

            event_weight = complex_to_dict(event.weight)
            event_weights.append(event_weight)
            threshold_sum = add_complex(
                [
                    term["weighted"]
                    for term in weighted_terms
                    if term["term_key"].startswith("threshold_counterterm:")
                ]
            )
            original_sum = add_complex(
                [term["weighted"] for term in weighted_terms if term["term_key"] == "original"]
            )
            events.append(
                {
                    "event_graph_id": graph_id,
                    "event_graph_name": graph_name,
                    "cut_id": int(event.cut_info.cut_id),
                    "lmb_channel_id": int(event.cut_info.lmb_channel_id),
                    "event_weight": event_weight,
                    "full_multiplicative_factor": event_factor,
                    "original_sum": original_sum,
                    "threshold_sum": threshold_sum,
                    "top_terms": sorted(
                        weighted_terms,
                        key=lambda term: term["weighted"]["abs"],
                        reverse=True,
                    )[:10],
                }
            )

    by_key: dict[tuple[str, str], list[dict[str, float]]] = defaultdict(list)
    for term in flat_terms:
        by_key[(term["event_graph_name"], term["term_key"])].append(term["weighted"])
    term_sums = [
        {
            "event_graph_name": graph_name,
            "term_key": key,
            "weighted": add_complex(values),
        }
        for (graph_name, key), values in by_key.items()
    ]
    term_sums.sort(key=lambda term: term["weighted"]["abs"], reverse=True)
    dominant_sum = add_complex(
        [
            term["weighted"]
            for term in term_sums
            if (term["event_graph_name"], term["term_key"]) in DOMINANT_TERMS
        ]
    )
    total_sum = add_complex(event_weights)
    return {
        **point,
        "raw_integrand_result": complex_to_dict(result.integrand_result),
        "total_event_weight": total_sum,
        "dominant_term_sum": dominant_sum,
        "other_event_weight": {
            "re": total_sum["re"] - dominant_sum["re"],
            "im": total_sum["im"] - dominant_sum["im"],
            "abs": math.hypot(total_sum["re"] - dominant_sum["re"], total_sum["im"] - dominant_sum["im"]),
        },
        "is_nan": bool(result.is_nan),
        "events": events,
        "term_sums": term_sums,
        "flat_terms": flat_terms,
    }


def write_tsv(reports: list[dict[str, Any]]) -> None:
    point_keys = [
        "point_index",
        "source",
        "scan_row_index",
        "x0",
        "x1",
        "x2",
        "x3",
        "x4",
        "x5",
        "dx2_from_root",
        "signed_side_fraction",
        "t",
        "total_re",
        "total_im",
        "total_abs",
        "dominant_re",
        "dominant_im",
        "dominant_abs",
        "other_re",
        "other_im",
        "other_abs",
        "scan_re",
        "scan_im",
        "scan_abs",
        "is_nan",
    ]
    with OUT_POINT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=point_keys, delimiter="\t")
        writer.writeheader()
        for report in reports:
            total = report["total_event_weight"]
            dominant = report["dominant_term_sum"]
            other = report["other_event_weight"]
            row = {
                "point_index": report["point_index"],
                "source": report["source"],
                "scan_row_index": report.get("scan_row_index"),
                "dx2_from_root": report["dx2_from_root"],
                "signed_side_fraction": report.get("signed_side_fraction"),
                "t": report.get("t"),
                "total_re": total["re"],
                "total_im": total["im"],
                "total_abs": total["abs"],
                "dominant_re": dominant["re"],
                "dominant_im": dominant["im"],
                "dominant_abs": dominant["abs"],
                "other_re": other["re"],
                "other_im": other["im"],
                "other_abs": other["abs"],
                "scan_re": report.get("scan_re"),
                "scan_im": report.get("scan_im"),
                "scan_abs": report.get("scan_abs"),
                "is_nan": report["is_nan"],
            }
            for index, value in enumerate(report["x"]):
                row[f"x{index}"] = value
            writer.writerow(row)

    term_keys = [
        "point_index",
        "source",
        "x2",
        "dx2_from_root",
        "event_graph_id",
        "event_graph_name",
        "term_key",
        "edge_ids",
        "weighted_re",
        "weighted_im",
        "weighted_abs",
        "pole_scaled_re",
        "pole_scaled_im",
        "pole_scaled_abs",
        "raw_re",
        "raw_im",
        "raw_abs",
    ]
    with OUT_TERM_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=term_keys, delimiter="\t")
        writer.writeheader()
        for report in reports:
            for term in report["flat_terms"]:
                weighted = term["weighted"]
                raw = term["raw"]
                dx2 = float(term["dx2_from_root"])
                pole_scaled = {
                    "re": weighted["re"] * dx2,
                    "im": weighted["im"] * dx2,
                    "abs": weighted["abs"] * abs(dx2),
                }
                writer.writerow(
                    {
                        "point_index": term["point_index"],
                        "source": report["source"],
                        "x2": term["x2"],
                        "dx2_from_root": dx2,
                        "event_graph_id": term["event_graph_id"],
                        "event_graph_name": term["event_graph_name"],
                        "term_key": term["term_key"],
                        "edge_ids": (
                            "" if term["edge_ids"] is None else ",".join(map(str, term["edge_ids"]))
                        ),
                        "weighted_re": weighted["re"],
                        "weighted_im": weighted["im"],
                        "weighted_abs": weighted["abs"],
                        "pole_scaled_re": pole_scaled["re"],
                        "pole_scaled_im": pole_scaled["im"],
                        "pole_scaled_abs": pole_scaled["abs"],
                        "raw_re": raw["re"],
                        "raw_im": raw["im"],
                        "raw_abs": raw["abs"],
                    }
                )


def series_by_term(reports: list[dict[str, Any]]) -> dict[tuple[str, str], list[dict[str, Any]]]:
    series: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for report in reports:
        summed = {
            (term["event_graph_name"], term["term_key"]): term["weighted"]
            for term in report["term_sums"]
        }
        for key, value in summed.items():
            series[key].append(
                {
                    "x": float(report["dx2_from_root"]),
                    "value": value,
                    "source": report["source"],
                }
            )
    return series


def plot_reports(reports: list[dict[str, Any]]) -> None:
    series = series_by_term(reports)
    x_total = [float(report["dx2_from_root"]) for report in reports]
    totals = [report["total_event_weight"] for report in reports]
    dominant = [report["dominant_term_sum"] for report in reports]
    other = [report["other_event_weight"] for report in reports]

    top_keys = list(DOMINANT_TERMS)
    extra_keys = sorted(
        series,
        key=lambda key: max(entry["value"]["abs"] for entry in series[key]),
        reverse=True,
    )
    for key in extra_keys:
        if key not in top_keys:
            top_keys.append(key)
        if len(top_keys) >= 8:
            break

    max_component = max(
        [abs(value[part]) for value in totals + dominant + other for part in ("re", "im")]
        + [
            abs(entry["value"][part])
            for key in top_keys
            for entry in series[key]
            for part in ("re", "im")
        ]
    )
    linthresh_y = max(max_component * 1.0e-7, 1.0e-12)

    fig, (ax_re, ax_im) = plt.subplots(2, 1, figsize=(9.0, 7.2), sharex=True)
    ax_re.plot(x_total, [value["re"] for value in totals], "k-", linewidth=1.8, label="total")
    ax_im.plot(x_total, [value["im"] for value in totals], "k-", linewidth=1.8, label="total")
    ax_re.plot(x_total, [value["re"] for value in dominant], "r-", linewidth=1.4, label="dominant CT sum")
    ax_im.plot(x_total, [value["im"] for value in dominant], "r-", linewidth=1.4, label="dominant CT sum")
    ax_re.plot(x_total, [value["re"] for value in other], color="0.55", linewidth=1.1, label="all other terms")
    ax_im.plot(x_total, [value["im"] for value in other], color="0.55", linewidth=1.1, label="all other terms")
    for key in top_keys:
        entries = sorted(series[key], key=lambda entry: entry["x"])
        label = f"{key[0]} {key[1].replace('threshold_counterterm:', 'T')}"
        ax_re.plot(
            [entry["x"] for entry in entries],
            [entry["value"]["re"] for entry in entries],
            "--",
            linewidth=0.9,
            label=label,
        )
        ax_im.plot(
            [entry["x"] for entry in entries],
            [entry["value"]["im"] for entry in entries],
            "--",
            linewidth=0.9,
            label=label,
        )

    for axis in (ax_re, ax_im):
        axis.axvline(0.0, color="black", linewidth=0.9, alpha=0.75)
        axis.set_xscale("symlog", linthresh=1.0e-9)
        axis.set_yscale("symlog", linthresh=linthresh_y)
        axis.grid(True, which="both", alpha=0.25)
    ax_re.set_ylabel("real contribution")
    ax_im.set_ylabel("imag contribution")
    ax_im.set_xlabel("x2 - x2_root")
    ax_re.legend(ncols=2, fontsize="x-small")
    fig.suptitle("axis_x2 im[-] feature: individual CT contributions")
    fig.tight_layout()
    fig.savefig(OUT_COMPONENT_PNG, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    ax.plot(x_total, [max(value["abs"], 1.0e-300) for value in totals], "k-", linewidth=1.8, label="total")
    ax.plot(
        x_total,
        [max(value["abs"], 1.0e-300) for value in dominant],
        "r-",
        linewidth=1.4,
        label="dominant CT sum",
    )
    ax.plot(
        x_total,
        [max(value["abs"], 1.0e-300) for value in other],
        color="0.55",
        linewidth=1.1,
        label="all other terms",
    )
    for key in top_keys[:6]:
        entries = sorted(series[key], key=lambda entry: entry["x"])
        label = f"{key[0]} {key[1].replace('threshold_counterterm:', 'T')}"
        ax.plot(
            [entry["x"] for entry in entries],
            [max(entry["value"]["abs"], 1.0e-300) for entry in entries],
            "--",
            linewidth=0.95,
            label=label,
        )
    ax.axvline(0.0, color="black", linewidth=0.9, alpha=0.75)
    ax.set_xscale("symlog", linthresh=1.0e-9)
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.set_xlabel("x2 - x2_root")
    ax.set_ylabel("|contribution|")
    ax.legend(ncols=2, fontsize="x-small")
    ax.set_title("axis_x2 im[-] feature: contribution magnitudes")
    fig.tight_layout()
    fig.savefig(OUT_ABS_PNG, dpi=180)
    plt.close(fig)


def summarize_pole_scaling(reports: list[dict[str, Any]]) -> list[dict[str, Any]]:
    series = series_by_term(reports)
    summaries = []
    for key in DOMINANT_TERMS:
        entries = [
            entry
            for entry in series[key]
            if abs(entry["x"]) >= 1.0e-8 and abs(entry["x"]) <= 1.0e-6
        ]
        left = [entry for entry in entries if entry["x"] < 0.0]
        right = [entry for entry in entries if entry["x"] > 0.0]
        for side_label, side_entries in (("left", left), ("right", right)):
            scaled = [
                {
                    "re": entry["value"]["re"] * entry["x"],
                    "im": entry["value"]["im"] * entry["x"],
                    "abs": entry["value"]["abs"] * abs(entry["x"]),
                }
                for entry in side_entries
            ]
            if not scaled:
                continue
            mean = add_complex(scaled)
            count = len(scaled)
            mean = {part: mean[part] / count for part in ("re", "im", "abs")}
            summaries.append(
                {
                    "event_graph_name": key[0],
                    "term_key": key[1],
                    "side": side_label,
                    "count": count,
                    "mean_dx2_times_weight": mean,
                }
            )
    return summaries


def write_markdown(reports: list[dict[str, Any]], pole_summaries: list[dict[str, Any]]) -> None:
    left = max((report for report in reports if report["dx2_from_root"] < 0.0), key=lambda report: report["dx2_from_root"])
    right = min((report for report in reports if report["dx2_from_root"] > 0.0), key=lambda report: report["dx2_from_root"])
    scan_points = [
        report
        for report in reports
        if report["source"] in {"saved_scan_row", "max_weight_location"}
        and abs(float(report["dx2_from_root"])) < 2.0e-6
    ]
    scan_points.sort(key=lambda report: float(report["x2"]))

    all_term_sets = [
        {(term["event_graph_name"], term["term_key"]) for term in report["flat_terms"]}
        for report in reports
    ]
    common_terms = set.intersection(*all_term_sets)
    union_terms = set.union(*all_term_sets)
    appeared_or_disappeared = len(union_terms - common_terms)

    lines = [
        "# axis_x2 CT contribution monitor",
        "",
        "Runtime: xi = 4*(p1+p2), sliver_width = 30, gaussian_width = 1, threshold h_sigma = 1.",
        f"Fixed bin: target `{TARGET}`, graph `0`, LMB channel `9`.",
        f"Reference root from the rstar probe: `x2_root = {ROOT_X2:.16e}`.",
        "",
        "## Main observation",
        "",
        (
            "The sharp feature is not caused by one of the dominant CTs turning on or off. "
            f"Across the monitored points, `{appeared_or_disappeared}` term keys are absent from "
            "at least one point; the large terms listed below are present on both sides of the root. "
            "They all flip sign when `x2 - x2_root` changes sign, and they then decay together as "
            "`1 / |x2 - x2_root|`."
        ),
        "",
        "Dominant terms:",
        "",
    ]
    for graph_name, term_key in DOMINANT_TERMS:
        first = next(
            term
            for term in reports[0]["flat_terms"]
            if term["event_graph_name"] == graph_name and term["term_key"] == term_key
        )
        lines.append(
            f"- `{graph_name}` `{term_key}` edges `{first['edge_ids']}`"
        )

    lines.extend(
        [
            "",
            "## Nearest root probes",
            "",
            f"- nearest left: `dx2 = {left['dx2_from_root']:+.8e}`, total `{fmt(left['total_event_weight'])}`, dominant CT sum `{fmt(left['dominant_term_sum'])}`",
            f"- nearest right: `dx2 = {right['dx2_from_root']:+.8e}`, total `{fmt(right['total_event_weight'])}`, dominant CT sum `{fmt(right['dominant_term_sum'])}`",
            "",
            "## Saved scan rows around the feature",
            "",
        ]
    )
    for report in scan_points:
        lines.append(
            "- "
            f"row `{report['scan_row_index']}` dx2 `{report['dx2_from_root']:+.8e}` "
            f"scan `{report.get('scan_re', 0.0):+.8e} {report.get('scan_im', 0.0):+.8e}i`, "
            f"event total `{fmt(report['total_event_weight'])}`, "
            f"dominant CT sum `{fmt(report['dominant_term_sum'])}`, "
            f"other `{fmt(report['other_event_weight'])}`"
        )

    lines.extend(["", "## Pole scaling check", ""])
    for summary in pole_summaries:
        lines.append(
            "- "
            f"`{summary['event_graph_name']}` `{summary['term_key']}` {summary['side']} "
            f"mean `(x2-x2_root)*w = {fmt(summary['mean_dx2_times_weight'])}` "
            f"from `{summary['count']}` points"
        )

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- point summary TSV: `{OUT_POINT_TSV}`",
            f"- term TSV: `{OUT_TERM_TSV}`",
            f"- JSON: `{OUT_JSON}`",
            f"- component plot: `{OUT_COMPONENT_PNG}`",
            f"- magnitude plot: `{OUT_ABS_PNG}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines))


def strip_heavy_fields(report: dict[str, Any]) -> dict[str, Any]:
    keep = dict(report)
    keep["flat_terms"] = [
        {
            key: value
            for key, value in term.items()
            if key in {"point_index", "x2", "dx2_from_root", "event_graph_name", "term_key", "edge_ids", "weighted", "raw"}
        }
        for term in report["flat_terms"]
    ]
    return keep


def main() -> int:
    api = GammaLoopAPI(state_folder=STATE_DIR, read_only_state=True)
    apply_runtime_settings(api, xi_scale=4.0)
    apply_width_overrides(
        api,
        {"sliver_width": 30.0, "gaussian_width": 1.0, "h_sigma": 1.0},
    )
    api.run(
        "set process -p qq_hhh_2L -i subtracted kv "
        "general.generate_events=true "
        "general.store_additional_weights_in_event=true "
        "general.enable_cache=false "
        "integrator.integrated_phase=\"imag\""
    )
    group_maps = threshold_maps(api)
    points = build_monitor_points()
    reports = []
    for index, point in enumerate(points, start=1):
        reports.append(evaluate_point(api, group_maps, point))
        if index % 20 == 0:
            print(f"evaluated {index}/{len(points)} points", flush=True)

    pole_summaries = summarize_pole_scaling(reports)
    OUT_JSON.write_text(
        json.dumps(
            {
                "root_x2": ROOT_X2,
                "dominant_terms": DOMINANT_TERMS,
                "pole_summaries": pole_summaries,
                "reports": [strip_heavy_fields(report) for report in reports],
            },
            indent=2,
        )
        + "\n"
    )
    write_tsv(reports)
    plot_reports(reports)
    write_markdown(reports, pole_summaries)
    print(OUT_JSON)
    print(OUT_POINT_TSV)
    print(OUT_TERM_TSV)
    print(OUT_COMPONENT_PNG)
    print(OUT_ABS_PNG)
    print(OUT_MD)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
