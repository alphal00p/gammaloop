#!/usr/bin/env python3
"""Decompose qq_hhh_2L imaginary max-weight samples by threshold counterterm."""

from __future__ import annotations

import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gammaloop import GammaLoopAPI

from scan_imag_maxweight_approaches import (
    OUTPUT_DIR,
    STATE_DIR,
    TARGETS,
    apply_runtime_settings,
)


SCAN_JSON = OUTPUT_DIR / "imag_maxweight_approach_scan.json"
REPORT_DIR = OUTPUT_DIR / "threshold_decomposition"
THRESHOLD_KEY_RE = re.compile(r"^threshold_counterterm:(?P<esurface_id>\d+):(?P<overlap>\d+)$")


@dataclass(frozen=True)
class ProbePoint:
    label: str
    target: str
    graph_group: int
    lmb_channel: int
    x: tuple[float, float, float, float, float, float]
    source: str


def complex_to_dict(value: Any) -> dict[str, float]:
    re_val = float(value.re if hasattr(value, "re") else value.real)
    im_val = float(value.im if hasattr(value, "im") else value.imag)
    return {"re": re_val, "im": im_val, "abs": math.hypot(re_val, im_val)}


def multiply_complex(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    re_val = a["re"] * b["re"] - a["im"] * b["im"]
    im_val = a["re"] * b["im"] + a["im"] * b["re"]
    return {"re": re_val, "im": im_val, "abs": math.hypot(re_val, im_val)}


def add_complex(values: list[dict[str, float]]) -> dict[str, float]:
    re_val = sum(value["re"] for value in values)
    im_val = sum(value["im"] for value in values)
    return {"re": re_val, "im": im_val, "abs": math.hypot(re_val, im_val)}


def fmt_complex(value: dict[str, float]) -> str:
    return f"{value['re']:+.8e} {value['im']:+.8e}i"


def load_scan_maxima() -> list[ProbePoint]:
    if not SCAN_JSON.exists():
        return []
    rows = json.loads(SCAN_JSON.read_text())
    probes: list[ProbePoint] = []
    for target_label in sorted({row["target"] for row in rows}):
        for direction in sorted({row["direction"] for row in rows if row["target"] == target_label}):
            subset = [
                row
                for row in rows
                if row["target"] == target_label and row["direction"] == direction
            ]
            row = max(subset, key=lambda item: float(item["abs"]))
            probes.append(
                ProbePoint(
                    label=f"{target_label}_{direction}_scan_max",
                    target=target_label,
                    graph_group=int(row["graph"]),
                    lmb_channel=int(row["lmb_channel"]),
                    x=tuple(float(value) for value in row["x"]),  # type: ignore[arg-type]
                    source=(
                        "line_scan_max "
                        f"direction={direction} signed_side_fraction={row['signed_side_fraction']}"
                    ),
                )
            )
    return probes


def build_probe_points() -> list[ProbePoint]:
    probes = [
        ProbePoint(
            label=f"{target.label}_exact_cluster",
            target=target.label,
            graph_group=target.graph,
            lmb_channel=target.lmb_channel,
            x=target.x_peak,
            source=f"cluster {target.reported_component} {target.reported_value:+.8e}",
        )
        for target in TARGETS.values()
    ]
    probes.extend(load_scan_maxima())
    seen: set[tuple[int, int, tuple[float, ...]]] = set()
    unique: list[ProbePoint] = []
    for probe in probes:
        key = (probe.graph_group, probe.lmb_channel, tuple(round(value, 18) for value in probe.x))
        if key not in seen:
            seen.add(key)
            unique.append(probe)
    return unique


def threshold_maps(api: GammaLoopAPI) -> dict[int, dict[str, Any]]:
    info = api.get_integrand_info(integrand_name="subtracted")
    maps: dict[int, dict[str, Any]] = {}
    for group in info.graph_groups:
        maps[int(group.group_id)] = {
            "thresholds": {
                int(esurface.esurface_id): [int(edge_id) for edge_id in esurface.edge_ids]
                for esurface in group.threshold_esurfaces
            },
            "graphs": {
                int(graph.graph_id): {
                    "name": graph.name,
                    "is_master": bool(graph.is_master),
                }
                for graph in group.graphs
            },
            "lmbs": {
                int(lmb.channel_id): {
                    "basis_id": int(lmb.basis_id),
                    "edge_ids": [int(edge_id) for edge_id in lmb.edge_ids],
                    "matches_generation_basis": bool(lmb.matches_generation_basis),
                }
                for lmb in group.loop_momentum_bases
            },
        }
    return maps


def enrich_weight_key(key: str, group_map: dict[str, Any]) -> dict[str, Any]:
    match = THRESHOLD_KEY_RE.match(key)
    if not match:
        return {"kind": key}
    esurface_id = int(match.group("esurface_id"))
    return {
        "kind": "threshold_counterterm",
        "esurface_id": esurface_id,
        "overlap_group": int(match.group("overlap")),
        "edge_ids": group_map["thresholds"].get(esurface_id),
    }


def evaluate_probe(api: GammaLoopAPI, group_maps: dict[int, dict[str, Any]], probe: ProbePoint) -> dict[str, Any]:
    result = api.evaluate_sample(
        list(probe.x),
        integrand_name="subtracted",
        return_events=True,
        discrete_dim=[probe.graph_group, probe.lmb_channel],
    )
    raw_result = complex_to_dict(result.integrand_result)
    full_factor = {
        "re": float((1.0 if result.parameterization_jacobian is None else result.parameterization_jacobian))
        * float(result.integrator_weight),
        "im": 0.0,
    }
    weighted_result = multiply_complex(raw_result, full_factor)
    group_map = group_maps[probe.graph_group]
    events = []
    for event_group in result.event_groups:
        for event in event_group.events:
            cut_info = event.cut_info
            event_weight = complex_to_dict(event.weight)
            event_full_factor = full_factor
            raw_additional = []
            for additional_weight in event.additional_weights:
                key = str(additional_weight.key)
                value = complex_to_dict(additional_weight.value)
                if key == "full_multiplicative_factor":
                    event_full_factor = value
                    continue
                raw_additional.append(
                    {
                        "key": key,
                        "raw": value,
                        "weighted": multiply_complex(value, event_full_factor),
                        **enrich_weight_key(key, group_map),
                    }
                )
            raw_additional.sort(key=lambda entry: entry["weighted"]["abs"], reverse=True)
            ct_entries = [
                entry
                for entry in raw_additional
                if entry.get("kind") == "threshold_counterterm"
            ]
            original = next((entry for entry in raw_additional if entry["key"] == "original"), None)
            ct_weighted_sum = add_complex([entry["weighted"] for entry in ct_entries])
            raw_weighted_sum = add_complex([entry["weighted"] for entry in raw_additional])
            graph_info = group_map["graphs"].get(int(cut_info.graph_id), {})
            events.append(
                {
                    "graph_id": int(cut_info.graph_id),
                    "graph_name": graph_info.get("name"),
                    "graph_is_master": graph_info.get("is_master"),
                    "graph_group_id": int(cut_info.graph_group_id),
                    "cut_id": int(cut_info.cut_id),
                    "lmb_channel_id": int(cut_info.lmb_channel_id),
                    "lmb_channel_edge_ids": [int(edge_id) for edge_id in cut_info.lmb_channel_edge_ids],
                    "event_weight": event_weight,
                    "full_factor": event_full_factor,
                    "weighted_additional_sum": raw_weighted_sum,
                    "original_weighted": original["weighted"] if original else None,
                    "ct_weighted_sum": ct_weighted_sum,
                    "top_weights": raw_additional[:8],
                }
            )
    events.sort(key=lambda entry: entry["event_weight"]["abs"], reverse=True)
    return {
        "label": probe.label,
        "target": probe.target,
        "source": probe.source,
        "graph_group": probe.graph_group,
        "lmb_channel": probe.lmb_channel,
        "lmb": group_map["lmbs"].get(probe.lmb_channel),
        "x": list(probe.x),
        "raw_integrand_result": raw_result,
        "full_factor": full_factor,
        "weighted_integrand_result": weighted_result,
        "parameterization_jacobian": result.parameterization_jacobian,
        "integrator_weight": result.integrator_weight,
        "is_nan": bool(result.is_nan),
        "events": events,
    }


def write_tsv(reports: list[dict[str, Any]], path: Path) -> None:
    keys = [
        "probe",
        "source",
        "graph_group",
        "lmb_channel",
        "event_graph_id",
        "event_graph_name",
        "event_weight_re",
        "event_weight_im",
        "event_weight_abs",
        "rank",
        "weight_key",
        "kind",
        "esurface_id",
        "edge_ids",
        "raw_re",
        "raw_im",
        "raw_abs",
        "weighted_re",
        "weighted_im",
        "weighted_abs",
    ]
    with path.open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for report in reports:
            for event in report["events"]:
                for rank, entry in enumerate(event["top_weights"], start=1):
                    values = {
                        "probe": report["label"],
                        "source": report["source"],
                        "graph_group": report["graph_group"],
                        "lmb_channel": report["lmb_channel"],
                        "event_graph_id": event["graph_id"],
                        "event_graph_name": event["graph_name"],
                        "event_weight_re": event["event_weight"]["re"],
                        "event_weight_im": event["event_weight"]["im"],
                        "event_weight_abs": event["event_weight"]["abs"],
                        "rank": rank,
                        "weight_key": entry["key"],
                        "kind": entry.get("kind"),
                        "esurface_id": entry.get("esurface_id"),
                        "edge_ids": ",".join(str(edge_id) for edge_id in (entry.get("edge_ids") or [])),
                        "raw_re": entry["raw"]["re"],
                        "raw_im": entry["raw"]["im"],
                        "raw_abs": entry["raw"]["abs"],
                        "weighted_re": entry["weighted"]["re"],
                        "weighted_im": entry["weighted"]["im"],
                        "weighted_abs": entry["weighted"]["abs"],
                    }
                    handle.write("\t".join(str(values.get(key, "")) for key in keys) + "\n")


def write_markdown(reports: list[dict[str, Any]], path: Path) -> None:
    lines = [
        "# qq_hhh_2L Imaginary Max-Weight Threshold Decomposition",
        "",
        "Weights below are fixed graph-group/LMB-channel weights. The raw `additional_weights` "
        "entries are multiplied by the listed full multiplicative factor before comparison.",
        "",
    ]
    for report in reports:
        dominant_event = report["events"][0]
        dominant_weight = dominant_event["event_weight"]
        lines.extend(
            [
                f"## {report['label']}",
                "",
                f"- source: {report['source']}",
                f"- fixed bin: graph group {report['graph_group']}, LMB channel {report['lmb_channel']} "
                f"edges {report.get('lmb', {}).get('edge_ids')}",
                f"- weighted sample result: {fmt_complex(report['weighted_integrand_result'])} "
                f"(|w|={report['weighted_integrand_result']['abs']:.8e})",
                f"- dominant graph term: {dominant_event['graph_id']} "
                f"({dominant_event['graph_name']}), event weight {fmt_complex(dominant_weight)} "
                f"(|w|={dominant_weight['abs']:.8e})",
                "",
                "| rank | key | threshold edges | weighted value | |weighted| |",
                "| ---: | --- | --- | ---: | ---: |",
            ]
        )
        for rank, entry in enumerate(dominant_event["top_weights"][:6], start=1):
            edge_ids = entry.get("edge_ids")
            lines.append(
                f"| {rank} | `{entry['key']}` | {edge_ids if edge_ids is not None else ''} | "
                f"{fmt_complex(entry['weighted'])} | {entry['weighted']['abs']:.8e} |"
            )
        lines.append("")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    if not STATE_DIR.exists():
        raise SystemExit(f"State folder missing: {STATE_DIR}")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    api = GammaLoopAPI(state_folder=STATE_DIR, read_only_state=True)
    apply_runtime_settings(api)
    api.run(
        "set process -p qq_hhh_2L -i subtracted kv "
        "general.generate_events=true "
        "general.store_additional_weights_in_event=true "
        "general.enable_cache=false "
        "integrator.integrated_phase=\"imag\""
    )
    group_maps = threshold_maps(api)
    probes = build_probe_points()
    reports = []
    for index, probe in enumerate(probes, start=1):
        print(f"[{index}/{len(probes)}] evaluating {probe.label}", flush=True)
        reports.append(evaluate_probe(api, group_maps, probe))
    (REPORT_DIR / "threshold_decomposition.json").write_text(json.dumps(reports, indent=2) + "\n")
    write_tsv(reports, REPORT_DIR / "threshold_decomposition.tsv")
    write_markdown(reports, REPORT_DIR / "threshold_decomposition.md")
    print(f"Wrote threshold decomposition to {REPORT_DIR}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
