#!/usr/bin/env python3
"""Compare event-level decomposition across the im[-] axis_x2 drop."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gammaloop import GammaLoopAPI

from inspect_imag_maxweight_thresholds import complex_to_dict, multiply_complex, threshold_maps
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
OUT_JSON = APPROACH_DIR / "axis_x2_drop_decomposition.json"
OUT_TSV = APPROACH_DIR / "axis_x2_drop_weight_terms.tsv"
OUT_MD = APPROACH_DIR / "axis_x2_drop_diagnostic.md"


@dataclass(frozen=True)
class Probe:
    label: str
    row_index: int
    row: dict[str, str]


def cadd(values: list[dict[str, float]]) -> dict[str, float]:
    re_val = sum(value["re"] for value in values)
    im_val = sum(value["im"] for value in values)
    return {"re": re_val, "im": im_val, "abs": math.hypot(re_val, im_val)}


def fmt(value: dict[str, float]) -> str:
    return f"{value['re']:+.8e} {value['im']:+.8e}i |w|={value['abs']:.8e}"


def load_probes() -> list[Probe]:
    rows = [
        row
        for row in csv.DictReader(SCAN_TSV.open(), delimiter="\t")
        if row["target"] == "highstat_iter_0001_imag_minus"
        and row["direction"] == "axis_x2"
    ]
    rows.sort(key=lambda row: float(row["signed_side_fraction"]))
    # The sharp local peak/drop is rows 32 -> 33 -> 34 in the saved 101-point scan.
    labels = {
        32: "left_rising",
        33: "local_peak_before_drop",
        34: "right_after_drop",
        35: "right_falling",
    }
    return [Probe(labels[index], index, rows[index]) for index in sorted(labels)]


def row_x(row: dict[str, str]) -> list[float]:
    return [float(row[f"x{i}"]) for i in range(6)]


def evaluate_probe(
    api: GammaLoopAPI,
    group_maps: dict[int, dict[str, Any]],
    probe: Probe,
) -> dict[str, Any]:
    row = probe.row
    graph = int(row["graph"])
    lmb_channel = int(row["lmb_channel"])
    result = api.evaluate_sample(
        row_x(row),
        integrand_name="subtracted",
        return_events=True,
        discrete_dim=[graph, lmb_channel],
    )
    group_map = group_maps[graph]
    events = []
    flat_terms = []
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
                weighted_terms.append({"key": key, "raw": raw, "weighted": weighted})
                flat_terms.append(
                    {
                        "probe": probe.label,
                        "row_index": probe.row_index,
                        "signed_side_fraction": float(row["signed_side_fraction"]),
                        "t": float(row["t"]),
                        "event_graph_id": graph_id,
                        "event_graph_name": graph_name,
                        "term_key": key,
                        "weighted": weighted,
                        "raw": raw,
                    }
                )

            threshold_terms = [
                term for term in weighted_terms if term["key"].startswith("threshold_counterterm:")
            ]
            original_terms = [term for term in weighted_terms if term["key"] == "original"]
            events.append(
                {
                    "event_graph_id": graph_id,
                    "event_graph_name": graph_name,
                    "cut_id": int(event.cut_info.cut_id),
                    "lmb_channel_id": int(event.cut_info.lmb_channel_id),
                    "event_weight": complex_to_dict(event.weight),
                    "full_multiplicative_factor": event_factor,
                    "original_sum": cadd([term["weighted"] for term in original_terms]),
                    "threshold_sum": cadd([term["weighted"] for term in threshold_terms]),
                    "all_terms_sum": cadd([term["weighted"] for term in weighted_terms]),
                    "top_terms": sorted(
                        weighted_terms,
                        key=lambda term: term["weighted"]["abs"],
                        reverse=True,
                    )[:12],
                }
            )
    events.sort(key=lambda event: event["event_weight"]["abs"], reverse=True)
    return {
        "label": probe.label,
        "row_index": probe.row_index,
        "signed_side_fraction": float(row["signed_side_fraction"]),
        "t": float(row["t"]),
        "x": row_x(row),
        "scan_re": float(row["re"]),
        "scan_im": float(row["im"]),
        "scan_abs": float(row["abs"]),
        "integrand_result": complex_to_dict(result.integrand_result),
        "is_nan": bool(result.is_nan),
        "events": events,
        "flat_terms": flat_terms,
    }


def write_tsv(reports: list[dict[str, Any]]) -> None:
    keys = [
        "probe",
        "row_index",
        "signed_side_fraction",
        "t",
        "event_graph_id",
        "event_graph_name",
        "term_key",
        "weighted_re",
        "weighted_im",
        "weighted_abs",
        "raw_re",
        "raw_im",
        "raw_abs",
    ]
    with OUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        for report in reports:
            for term in report["flat_terms"]:
                weighted = term["weighted"]
                raw = term["raw"]
                writer.writerow(
                    {
                        "probe": term["probe"],
                        "row_index": term["row_index"],
                        "signed_side_fraction": term["signed_side_fraction"],
                        "t": term["t"],
                        "event_graph_id": term["event_graph_id"],
                        "event_graph_name": term["event_graph_name"],
                        "term_key": term["term_key"],
                        "weighted_re": weighted["re"],
                        "weighted_im": weighted["im"],
                        "weighted_abs": weighted["abs"],
                        "raw_re": raw["re"],
                        "raw_im": raw["im"],
                        "raw_abs": raw["abs"],
                    }
                )


def summarize_changes(reports: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_probe = {report["label"]: report for report in reports}
    before = by_probe["local_peak_before_drop"]
    after = by_probe["right_after_drop"]
    before_terms = {
        (term["event_graph_name"], term["term_key"]): term["weighted"]
        for term in before["flat_terms"]
    }
    after_terms = {
        (term["event_graph_name"], term["term_key"]): term["weighted"]
        for term in after["flat_terms"]
    }
    changes = []
    for key in sorted(set(before_terms) | set(after_terms)):
        b = before_terms.get(key, {"re": 0.0, "im": 0.0, "abs": 0.0})
        a = after_terms.get(key, {"re": 0.0, "im": 0.0, "abs": 0.0})
        delta = {"re": a["re"] - b["re"], "im": a["im"] - b["im"]}
        delta["abs"] = math.hypot(delta["re"], delta["im"])
        changes.append(
            {
                "event_graph_name": key[0],
                "term_key": key[1],
                "before": b,
                "after": a,
                "delta_after_minus_before": delta,
                "ratio_abs_after_over_before": (
                    a["abs"] / b["abs"] if b["abs"] != 0.0 else math.inf
                ),
            }
        )
    changes.sort(key=lambda item: item["delta_after_minus_before"]["abs"], reverse=True)
    return changes


def write_markdown(reports: list[dict[str, Any]], changes: list[dict[str, Any]]) -> None:
    lines = [
        "# axis_x2 im[-] drop diagnostic",
        "",
        "Runtime: xi = 4*(p1+p2), sliver_width = 30, gaussian_width = 1, threshold h_sigma = 1.",
        "",
        "The compared rows bracket the sharp local drop on the im[-] axis_x2 approach line.",
        "",
        "## Probe points",
        "",
    ]
    for report in reports:
        lines.extend(
            [
                f"### {report['label']}",
                "",
                f"- row index: `{report['row_index']}`",
                f"- signed side fraction: `{report['signed_side_fraction']:+.16e}`",
                f"- t: `{report['t']:+.16e}`",
                f"- scan weight: `{report['scan_re']:+.8e} {report['scan_im']:+.8e}i`, |w| = `{report['scan_abs']:.8e}`",
                f"- API result: `{fmt(report['integrand_result'])}`",
                "",
                "Top event terms:",
                "",
            ]
        )
        for event in report["events"][:6]:
            lines.append(
                f"- graph `{event['event_graph_name']}` event `{fmt(event['event_weight'])}`, "
                f"original `{fmt(event['original_sum'])}`, threshold sum `{fmt(event['threshold_sum'])}`"
            )
            for term in event["top_terms"][:8]:
                lines.append(f"  - `{term['key']}`: `{fmt(term['weighted'])}`")
        lines.append("")

    lines.extend(
        [
            "## Largest term changes across local_peak_before_drop -> right_after_drop",
            "",
        ]
    )
    for change in changes[:30]:
        lines.append(
            "- "
            f"`{change['event_graph_name']}` `{change['term_key']}`: "
            f"before `{fmt(change['before'])}`, after `{fmt(change['after'])}`, "
            f"delta `{fmt(change['delta_after_minus_before'])}`, "
            f"|after|/|before| = `{change['ratio_abs_after_over_before']:.6e}`"
        )
    lines.append("")

    before_keys = {
        (term["event_graph_name"], term["term_key"])
        for report in reports
        if report["label"] == "local_peak_before_drop"
        for term in report["flat_terms"]
    }
    after_keys = {
        (term["event_graph_name"], term["term_key"])
        for report in reports
        if report["label"] == "right_after_drop"
        for term in report["flat_terms"]
    }
    disappeared = sorted(before_keys - after_keys)
    appeared = sorted(after_keys - before_keys)
    lines.extend(["## Presence changes", ""])
    lines.append(f"- terms disappeared: `{len(disappeared)}`")
    for graph_name, key in disappeared[:40]:
        lines.append(f"  - `{graph_name}` `{key}`")
    lines.append(f"- terms appeared: `{len(appeared)}`")
    for graph_name, key in appeared[:40]:
        lines.append(f"  - `{graph_name}` `{key}`")
    lines.append("")

    OUT_MD.write_text("\n".join(lines) + "\n")


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
    reports = [evaluate_probe(api, group_maps, probe) for probe in load_probes()]
    changes = summarize_changes(reports)
    OUT_JSON.write_text(json.dumps({"reports": reports, "changes": changes}, indent=2) + "\n")
    write_tsv(reports)
    write_markdown(reports, changes)
    print(OUT_JSON)
    print(OUT_TSV)
    print(OUT_MD)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
