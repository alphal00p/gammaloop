#!/usr/bin/env python3
"""Overlay threshold-CT existence probes on saved approach plots."""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from gammaloop import GammaLoopAPI

from inspect_imag_maxweight_thresholds import complex_to_dict
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
OUTPUT_TSV = APPROACH_DIR / "esurface_existence_overlay.tsv"
OUTPUT_JSON = APPROACH_DIR / "esurface_existence_overlay.json"


@dataclass(frozen=True)
class SurfaceProbe:
    color: str
    label: str
    graph_id: int
    graph_name: str
    key: str
    esurface_id: int
    edges: tuple[int, ...]
    equation: str
    signature_note: str


SURFACE_PROBES: dict[str, tuple[SurfaceProbe, SurfaceProbe]] = {
    "highstat_iter_0001_imag_plus": (
        SurfaceProbe(
            color="black",
            label="GL06_isr_p2_ct:T36",
            graph_id=5,
            graph_name="GL06_isr_p2_ct",
            key="threshold_counterterm:36:0",
            esurface_id=36,
            edges=(9, 14, 16),
            equation="OSE(9) + OSE(14) + OSE(16) + external_shift(T36) = 0",
            signature_note="LMB channel 14 edges [9,12]; threshold edges [9,14,16]",
        ),
        SurfaceProbe(
            color="red",
            label="GL06_isr_p2_ct:T13",
            graph_id=5,
            graph_name="GL06_isr_p2_ct",
            key="threshold_counterterm:13:0",
            esurface_id=13,
            edges=(11, 15),
            equation="OSE(11) + OSE(15) + external_shift(T13) = 0",
            signature_note="LMB channel 14 edges [9,12]; threshold edges [11,15]",
        ),
    ),
    "highstat_iter_0001_imag_minus": (
        SurfaceProbe(
            color="black",
            label="GL05_isr_p2_ct:T4",
            graph_id=2,
            graph_name="GL05_isr_p2_ct",
            key="threshold_counterterm:4:0",
            esurface_id=4,
            edges=(9, 12, 15),
            equation="OSE(9) + OSE(12) + OSE(15) + external_shift(T4) = 0",
            signature_note="LMB channel 9 edges [11,14]; threshold edges [9,12,15]",
        ),
        SurfaceProbe(
            color="red",
            label="GL05_isr_p1_ct:T34",
            graph_id=1,
            graph_name="GL05_isr_p1_ct",
            key="threshold_counterterm:34:0",
            esurface_id=34,
            edges=(15, 16),
            equation="OSE(15) + OSE(16) + external_shift(T34) = 0",
            signature_note="LMB channel 9 edges [11,14]; threshold edges [15,16]",
        ),
    ),
}


def load_rows() -> list[dict[str, Any]]:
    with SCAN_TSV.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def row_x(row: dict[str, Any]) -> list[float]:
    return [float(row[f"x{i}"]) for i in range(6)]


def event_has_surface(event: Any, probe: SurfaceProbe) -> tuple[bool, float]:
    if int(event.cut_info.graph_id) != probe.graph_id:
        return False, 0.0
    event_factor = {"re": 1.0, "im": 0.0}
    raw_value = None
    for additional_weight in event.additional_weights:
        key = str(additional_weight.key)
        value = complex_to_dict(additional_weight.value)
        if key == "full_multiplicative_factor":
            event_factor = value
        elif key == probe.key:
            raw_value = value
    if raw_value is None:
        return False, 0.0
    weighted_abs = raw_value["abs"] * math.hypot(event_factor["re"], event_factor["im"])
    return weighted_abs > 0.0, weighted_abs


def evaluate_existence(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
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

    out_rows: list[dict[str, Any]] = []
    for index, row in enumerate(rows, start=1):
        result = api.evaluate_sample(
            row_x(row),
            integrand_name="subtracted",
            return_events=True,
            discrete_dim=[int(row["graph"]), int(row["lmb_channel"])],
        )
        events = [event for group in result.event_groups for event in group.events]
        for probe in SURFACE_PROBES[row["target"]]:
            exists = False
            weighted_abs = 0.0
            for event in events:
                event_exists, event_weighted_abs = event_has_surface(event, probe)
                exists = exists or event_exists
                weighted_abs = max(weighted_abs, event_weighted_abs)
            out_rows.append(
                {
                    "target": row["target"],
                    "direction": row["direction"],
                    "signed_side_fraction": row["signed_side_fraction"],
                    "t": row["t"],
                    "surface_label": probe.label,
                    "color": probe.color,
                    "event_graph_id": probe.graph_id,
                    "event_graph_name": probe.graph_name,
                    "threshold_key": probe.key,
                    "esurface_id": probe.esurface_id,
                    "edge_ids": ",".join(str(edge_id) for edge_id in probe.edges),
                    "exists": int(exists),
                    "plot_value": 1.0 if exists else 1.0e-6,
                    "weighted_abs": weighted_abs,
                }
            )
        if index % 50 == 0:
            print(f"evaluated existence for {index}/{len(rows)} scan points", flush=True)
    return out_rows


def write_outputs(existence_rows: list[dict[str, Any]]) -> None:
    keys = [
        "target",
        "direction",
        "signed_side_fraction",
        "t",
        "surface_label",
        "color",
        "event_graph_id",
        "event_graph_name",
        "threshold_key",
        "esurface_id",
        "edge_ids",
        "exists",
        "plot_value",
        "weighted_abs",
    ]
    with OUTPUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(existence_rows)
    OUTPUT_JSON.write_text(json.dumps(existence_rows, indent=2) + "\n")


def signed(values: list[float]) -> list[float]:
    return values


def plot_target(scan_rows: list[dict[str, Any]], existence_rows: list[dict[str, Any]], target: str) -> Path:
    target_rows = [row for row in scan_rows if row["target"] == target]
    directions = sorted({row["direction"] for row in target_rows})
    reported = target_rows[0]["reported_component"]
    reported_value = float(target_rows[0]["reported_value"])
    graph = int(target_rows[0]["graph"])
    lmb_channel = int(target_rows[0]["lmb_channel"])

    fig, (ax_comp, ax_abs) = plt.subplots(2, 1, figsize=(13.0, 10.0), sharex=True)
    for direction in directions:
        subset = sorted(
            [row for row in target_rows if row["direction"] == direction],
            key=lambda row: float(row["signed_side_fraction"]),
        )
        x_axis = [float(row["signed_side_fraction"]) for row in subset]
        ax_comp.plot(x_axis, [float(row["re"]) for row in subset], linewidth=1.1, label=f"{direction} re")
        ax_comp.plot(
            x_axis,
            [float(row["im"]) for row in subset],
            "--",
            linewidth=1.1,
            label=f"{direction} im",
        )
        ax_abs.plot(x_axis, [float(row["abs"]) for row in subset], linewidth=1.2, label=direction)

    for probe in SURFACE_PROBES[target]:
        for direction in directions:
            subset = sorted(
                [
                    row
                    for row in existence_rows
                    if row["target"] == target
                    and row["direction"] == direction
                    and row["surface_label"] == probe.label
                ],
                key=lambda row: float(row["signed_side_fraction"]),
            )
            ax_abs.step(
                [float(row["signed_side_fraction"]) for row in subset],
                [float(row["plot_value"]) for row in subset],
                where="mid",
                color=probe.color,
                linewidth=1.5,
                alpha=0.85,
                label=f"{probe.label} exists" if direction == directions[0] else None,
            )

    for axis in (ax_comp, ax_abs):
        axis.axvline(0.0, color="black", linewidth=0.8, alpha=0.55)
        axis.set_xscale("symlog", linthresh=1.0e-6)
        axis.grid(True, which="both", alpha=0.25)
    ax_comp.set_yscale("symlog", linthresh=1.0e-7)
    ax_abs.set_yscale("log")
    ax_abs.set_ylim(1.0e-7, None)
    ax_comp.set_ylabel("fixed-bin weight component")
    ax_abs.set_ylabel("fixed-bin |weight| and existence")
    ax_abs.set_xlabel("signed fraction of available x-space line segment")
    ax_comp.legend(loc="best", fontsize=9)
    ax_abs.legend(loc="best", fontsize=9)
    fig.suptitle(
        f"{target}: graph {graph}, LMB channel {lmb_channel}, "
        f"reported {reported} {reported_value:.6e}\n"
        "black/red: selected threshold CT emitted (1) or absent (1e-6)",
        fontsize=15,
    )
    fig.tight_layout()
    out_path = APPROACH_DIR / f"{target}_approaches_with_existence.png"
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def main() -> int:
    rows = load_rows()
    existence_rows = evaluate_existence(rows)
    write_outputs(existence_rows)
    metadata = {
        "existence_convention": (
            "1 if the selected event-level threshold_counterterm key is emitted "
            "with nonzero weighted absolute value at that x point; 1e-6 otherwise"
        ),
        "surfaces": {
            target: [probe.__dict__ for probe in probes]
            for target, probes in SURFACE_PROBES.items()
        },
    }
    (APPROACH_DIR / "esurface_existence_overlay_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    for target in sorted(SURFACE_PROBES):
        print(plot_target(rows, existence_rows, target), flush=True)
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
