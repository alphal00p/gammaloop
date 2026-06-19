#!/usr/bin/env python3
"""Inspect fixed-bin low-stat qq_hhh_2L max-weight samples after threshold fixes."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gammaloop import GammaLoopAPI

from inspect_imag_maxweight_thresholds import complex_to_dict, enrich_weight_key, threshold_maps
from scan_imag_maxweight_approaches import apply_runtime_settings


ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_DIR = ROOT / "examples" / "cli" / "qq_hhh_2L"
STATE_DIR = EXAMPLE_DIR / "state_clean_probe"
OUTPUT_DIR = EXAMPLE_DIR / "maxweight_approach_scans" / "clean_lowstat_maxweight_checks"


@dataclass(frozen=True)
class ProbePoint:
    label: str
    graph_group: int
    lmb_channel: int
    x: tuple[float, float, float, float, float, float]
    source: str


PROBES = [
    ProbePoint(
        label="global_im_plus_4096",
        graph_group=1,
        lmb_channel=11,
        x=(
            0.71906345171460317,
            0.26040351173928078,
            0.80739917696325181,
            0.93515792609772008,
            0.57309382535847941,
            0.76652808410758222,
        ),
        source="integrate_debug_imag_4096 global im[+] +1.2437608010978304e1",
    ),
    ProbePoint(
        label="global_im_minus_4096",
        graph_group=1,
        lmb_channel=2,
        x=(
            0.72753172263837806,
            0.10085692334716734,
            0.79245886170786151,
            0.91491158631284397,
            0.90224993530646047,
            0.28701934851916089,
        ),
        source="integrate_debug_imag_4096 global im[-] -1.4665980292772911e1",
    ),
    ProbePoint(
        label="group0_im_plus_4096",
        graph_group=0,
        lmb_channel=6,
        x=(
            0.57160672727299811,
            0.027307468645824340,
            0.64649615051063458,
            0.59460577529663694,
            0.52756832492115968,
            0.97474224299388457,
        ),
        source="integrate_debug_imag_4096 group #0 im[+] +2.5115073073372836e0",
    ),
    ProbePoint(
        label="group0_im_minus_4096",
        graph_group=0,
        lmb_channel=9,
        x=(
            0.65050432546395509,
            0.74870642141957333,
            0.87489560807666433,
            0.51021409050547883,
            0.40840746778517723,
            0.78021812723544137,
        ),
        source="integrate_debug_imag_4096 group #0 im[-] -3.6578245054118312e0",
    ),
]


def multiply_by_real(value: dict[str, float], factor: float) -> dict[str, float]:
    re_val = value["re"] * factor
    im_val = value["im"] * factor
    return {"re": re_val, "im": im_val, "abs": math.hypot(re_val, im_val)}


def evaluate_probe(
    api: GammaLoopAPI, group_maps: dict[int, dict[str, Any]], probe: ProbePoint
) -> dict[str, Any]:
    result = api.evaluate_sample(
        list(probe.x),
        integrand_name="subtracted",
        return_events=True,
        discrete_dim=[probe.graph_group, probe.lmb_channel],
    )
    raw_result = complex_to_dict(result.integrand_result)
    full_factor = (1.0 if result.parameterization_jacobian is None else result.parameterization_jacobian) * float(
        result.integrator_weight
    )
    group_map = group_maps[probe.graph_group]
    events: list[dict[str, Any]] = []
    for event_group in result.event_groups:
        for event in event_group.events:
            weighted_terms: list[dict[str, Any]] = []
            event_full_factor = {"re": full_factor, "im": 0.0, "abs": abs(full_factor)}
            for additional_weight in event.additional_weights:
                key = str(additional_weight.key)
                value = complex_to_dict(additional_weight.value)
                if key == "full_multiplicative_factor":
                    event_full_factor = value
                    continue
                factor = float(event_full_factor["re"])
                weighted_terms.append(
                    {
                        "key": key,
                        "raw": value,
                        "weighted": multiply_by_real(value, factor),
                        **enrich_weight_key(key, group_map),
                    }
                )
            weighted_terms.sort(key=lambda item: item["weighted"]["abs"], reverse=True)
            cut_info = event.cut_info
            graph_info = group_map["graphs"].get(int(cut_info.graph_id), {})
            events.append(
                {
                    "graph_id": int(cut_info.graph_id),
                    "graph_name": graph_info.get("name"),
                    "lmb_channel_id": int(cut_info.lmb_channel_id),
                    "event_weight": complex_to_dict(event.weight),
                    "full_multiplicative_factor": event_full_factor,
                    "top_terms": weighted_terms[:12],
                }
            )
    events.sort(key=lambda item: item["event_weight"]["abs"], reverse=True)
    return {
        "label": probe.label,
        "source": probe.source,
        "graph_group": probe.graph_group,
        "lmb_channel": probe.lmb_channel,
        "x": probe.x,
        "raw_integrand_result": raw_result,
        "weighted_integrand_result": multiply_by_real(raw_result, full_factor),
        "parameterization_jacobian": result.parameterization_jacobian,
        "integrator_weight": result.integrator_weight,
        "is_nan": bool(result.is_nan),
        "events": events,
    }


def write_tsv(reports: list[dict[str, Any]], path: Path) -> None:
    keys = [
        "probe",
        "graph_group",
        "lmb_channel",
        "event_graph",
        "rank",
        "key",
        "esurface_id",
        "edge_ids",
        "weighted_re",
        "weighted_im",
        "weighted_abs",
    ]
    with path.open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for report in reports:
            for event in report["events"]:
                for rank, term in enumerate(event["top_terms"], start=1):
                    values = {
                        "probe": report["label"],
                        "graph_group": report["graph_group"],
                        "lmb_channel": report["lmb_channel"],
                        "event_graph": event["graph_name"],
                        "rank": rank,
                        "key": term["key"],
                        "esurface_id": term.get("esurface_id", ""),
                        "edge_ids": ",".join(str(edge_id) for edge_id in (term.get("edge_ids") or [])),
                        "weighted_re": term["weighted"]["re"],
                        "weighted_im": term["weighted"]["im"],
                        "weighted_abs": term["weighted"]["abs"],
                    }
                    handle.write("\t".join(str(values[key]) for key in keys) + "\n")


def main() -> None:
    if not STATE_DIR.exists():
        raise SystemExit(f"State folder missing: {STATE_DIR}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
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
    reports = [evaluate_probe(api, group_maps, probe) for probe in PROBES]
    (OUTPUT_DIR / "clean_lowstat_maxweight_checks.json").write_text(json.dumps(reports, indent=2) + "\n")
    write_tsv(reports, OUTPUT_DIR / "clean_lowstat_maxweight_checks.tsv")
    print(f"Wrote {len(reports)} clean low-stat max-weight checks to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
