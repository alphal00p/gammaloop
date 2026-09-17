"""Summarize the existing map and Havana probe outputs, without model calls."""

import json
from decimal import Decimal
from pathlib import Path

base = Path(__file__).resolve().parent
report = json.loads((base / "inverse_replay.json").read_text())
trace = json.loads((base / "raw_native_trace.json").read_text())
raw = [Decimal(x) for x in report["common_raw_coordinates_decimal"]]
recorded = [Decimal(x) for row in trace["raw"] for x in row]
summary = {
    "physical_point_trace_max_abs_difference": str(
        max(abs(a - b) for a, b in zip(raw, recorded))
    ),
    "original_spike_estimator": [-310.7628500264454, -3397.274568766552],
    "physical_integrand": report[
        "physical_integrand_inferred_from_exact_original_draw"
    ],
    "modes": [],
    "estimator_scope": "Each channel row is the actual production conditional estimator including frozen learned channel and continuous weights. The full adaptive mixture estimator is a separate diagnostic, not the production estimator.",
}
for mode in report["modes"]:
    rows = []
    max_error = 0.0
    for channel in mode["channels"]:
        row = {
            k: channel[k]
            for k in ["id", "name", "support", "inverse_error"]
            if k in channel
        }
        if channel.get("support"):
            predicted = channel["predicted_conditional_estimator_at_exact_common_raw"]
            actual = channel["actual_production_evaluation_at_rounded_inverse_cube"]
            actual_value = actual.get("complete_estimator")
            if actual_value:
                error = max(abs((a - b) / b) for a, b in zip(actual_value, predicted))
                max_error = max(max_error, error)
                row["relative_prediction_replay_difference"] = error
            row.update(
                {
                    k: channel[k]
                    for k in [
                        "discrete_probability",
                        "continuous_cube_density",
                        "full_outer_inverse_adaptive_density",
                        "raw_map_density",
                        "raw_partition_alpha",
                        "J_times_alpha",
                        "rounded_cube_raw_roundtrip_max_abs_GeV",
                    ]
                }
            )
            row["conditional_estimator"] = predicted
            row["actual_estimator"] = actual_value
            row["actual_precision"] = actual.get("precision")
            row["actual_is_nan"] = actual.get("metadata", {}).get("is_nan")
            row["posterior_channel_probability_at_common_raw_point"] = (
                channel["adaptive_raw_density_contribution"]
                / mode["adaptive_mixture_raw_density"]
            )
        rows.append(row)
    item = {
        k: mode[k]
        for k in [
            "mode",
            "completed_iteration",
            "completed_points",
            "source_policy",
            "raw_map_density_sum",
            "adaptive_mixture_raw_density",
            "hypothetical_full_adaptive_mixture_estimator",
        ]
    }
    item["physical_integrand_binary64_input_control"] = mode[
        "physical_integrand_control_binary64_raw_input"
    ].get("complete_estimator")
    item["max_relative_conditional_prediction_replay_difference"] = max_error
    item["channels"] = rows
    if mode["mode"] == "optimized":
        item["original_saved_sample"] = mode["retained_negative_max_sample"]
    summary["modes"].append(item)
if len(summary["modes"]) == 2:
    optimized, augmented = summary["modes"]
    summary["augmented_over_optimized_raw_map_density"] = (
        augmented["raw_map_density_sum"] / optimized["raw_map_density_sum"]
    )
    summary["augmented_over_optimized_adaptive_raw_density"] = (
        augmented["adaptive_mixture_raw_density"]
        / optimized["adaptive_mixture_raw_density"]
    )
if len(summary["modes"]) == 2:
    augmented = summary["modes"][1]
    frozen = json.loads(
        (base / "frozen_augmented/integration_result.json").read_text()
    )["slots"][0]["integral"]
    channel3 = next(row for row in augmented["channels"] if row["id"] == 3)
    summary["augmented_channel3_single_sample_scale_at_completed_N"] = {
        component: {
            "single_sample_mean_contribution": channel3["conditional_estimator"][index]
            / augmented["completed_points"],
            "fraction_of_current_absolute_mean": abs(
                channel3["conditional_estimator"][index]
                / augmented["completed_points"]
                / frozen["result"][component]
            ),
            "ratio_to_current_error": abs(
                channel3["conditional_estimator"][index]
                / augmented["completed_points"]
                / frozen["error"][component]
            ),
        }
        for index, component in enumerate(["re", "im"])
    }
(base / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
