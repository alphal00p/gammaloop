"""Summarize independent seeds and plot convergence, cost, and efficiency."""

import argparse
import csv
import itertools
import json
import math
import statistics
from pathlib import Path

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results", type=Path)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
data = json.loads(args.results.read_text())
by_key = {(run["beta"], run["mode"], run["seed"]): run for run in data["runs"]}
args.output.mkdir(parents=True, exist_ok=True)
rows = []
for (beta, mode), group in itertools.groupby(
    sorted(data["runs"], key=lambda run: (run["beta"], run["mode"])),
    key=lambda run: (run["beta"], run["mode"]),
):
    runs = list(group)
    estimates = [run["value"] for run in runs]
    mean = statistics.mean(estimates)
    # Equal-weight replication avoids weighting estimates by their own noisy errors.
    mean_error = math.sqrt(sum(run["error"] ** 2 for run in runs)) / len(runs)
    rms_error = math.sqrt(statistics.mean(run["error"] ** 2 for run in runs))
    target = runs[0]["target"]
    (samples,) = {run["result"]["slots"][0]["integral"]["neval"] for run in runs}
    stats = [run["result"]["slots"][0]["integration_statistics"] for run in runs]
    rows.append(
        {
            "beta": beta,
            "mode": mode,
            "seeds": len(runs),
            "samples_per_run": samples,
            "target": runs[0]["target"],
            "mean": mean,
            "mean_error": mean_error,
            "mean_pull": (mean - target) / mean_error if target is not None else None,
            "seed_stddev": statistics.stdev(estimates) if len(runs) > 1 else None,
            "rms_error": rms_error,
            "relative_error_percent": 100
            * rms_error
            / abs(target if target is not None else mean),
            "relative_error_denominator": "target"
            if target is not None
            else "sample mean",
            "wall_seconds": statistics.mean(run["wall_seconds"] for run in runs),
            "cpu_seconds": statistics.mean(run["cpu_seconds"] for run in runs),
            "variance_times_wall": statistics.mean(
                run["variance_times_wall"] for run in runs
            ),
            "evaluation_us": 1e6
            * statistics.mean(stat["average_total_time_seconds"] for stat in stats),
            "mapping_us": 1e6
            * statistics.mean(
                stat["average_parameterization_time_seconds"] for stat in stats
            ),
            "physical_us": 1e6
            * statistics.mean(stat["average_integrand_time_seconds"] for stat in stats),
            "max_rejected_percent": max(
                stat["nan_or_unstable_percentage"] for stat in stats
            ),
            "mean_f64_percent": statistics.mean(
                stat["f64_percentage"] for stat in stats
            ),
            "mean_quad_percent": statistics.mean(
                stat["f128_percentage"] for stat in stats
            ),
            "mean_arb_percent": statistics.mean(
                stat["arb_percentage"] for stat in stats
            ),
            "fermi_sample_percent": statistics.mean(
                100
                * sum(
                    entry["processed_samples"]
                    for entry in (
                        run["result"]["slots"][0]["grid_breakdown"]["im"] or {}
                    ).get("entries", [])
                    if entry["bin_label"].startswith("fermi_")
                )
                / run["result"]["slots"][0]["integral"]["neval"]
                if mode != "fermi"
                else 100
                for run in runs
            ),
            "max_abs_run_pull": max(abs(run["pull"]) for run in runs)
            if target is not None
            else None,
            "rmse_to_reference": math.sqrt(
                statistics.mean((run["value"] - run["target"]) ** 2 for run in runs)
            )
            if target is not None
            else None,
            "runs_within_2sigma": sum(abs(run["pull"]) <= 2 for run in runs)
            if target is not None
            else None,
        }
    )
for row in rows:
    for baseline_mode in ("ordinary", "optimized"):
        matching = [
            candidate
            for candidate in rows
            if candidate["beta"] == row["beta"] and candidate["mode"] == baseline_mode
        ]
        if not matching:
            continue
        (baseline,) = matching
        assert baseline["samples_per_run"] == row["samples_per_run"]
        row[f"efficiency_vs_{baseline_mode}"] = (
            baseline["variance_times_wall"] / row["variance_times_wall"]
        )
        seed_ratios = [
            by_key[(run["beta"], baseline_mode, run["seed"])]["variance_times_wall"]
            / run["variance_times_wall"]
            for run in data["runs"]
            if (run["beta"], run["mode"]) == (row["beta"], row["mode"])
        ]
        # Observed replication range, not a statistical confidence interval.
        row[f"efficiency_seed_min_vs_{baseline_mode}"] = min(seed_ratios)
        row[f"efficiency_seed_max_vs_{baseline_mode}"] = max(seed_ratios)
        gain = (baseline["rms_error"] / row["rms_error"]) ** 2
        row[f"variance_gain_vs_{baseline_mode}"] = gain
        # Forecast only: assume every draw incurs the same additional physical cost.
        row[f"extra_cost_break_even_us_vs_{baseline_mode}"] = (
            max(
                0,
                1e6
                * (row["wall_seconds"] - gain * baseline["wall_seconds"])
                / (row["samples_per_run"] * (gain - 1)),
            )
            if gain > 1
            else None
        )

with (args.output / "summary.csv").open("w") as destination:
    writer = csv.DictWriter(destination, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
(args.output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")

fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")
zero_temperature = data["metadata"].get("zero_temperature", False)
mode_labels = {
    "optimized": "Ordinary LMBs",
    "optimized_fermi": "Ordinary LMBs + Fermi shells",
    "optimized_fermi_repeated": "Ordinary LMBs + repeated-momentum shell",
}
for mode in ("optimized", "optimized_fermi", "optimized_fermi_repeated"):
    selected = [row for row in rows if row["mode"] == mode]
    if not selected:
        continue
    betas = [row["beta"] for row in selected]
    label = mode_labels[mode]
    if zero_temperature:
        (row,) = selected
        position = list(mode_labels).index(mode)
        for axis, metric in zip(
            axes,
            ("relative_error_percent", "wall_seconds", "efficiency_vs_optimized"),
            strict=True,
        ):
            axis.bar(position, row[metric], label=label)
    else:
        axes[0].loglog(
            betas,
            [row["relative_error_percent"] for row in selected],
            "o-",
            label=label,
        )
        axes[1].semilogx(
            betas, [row["wall_seconds"] for row in selected], "o-", label=label
        )
        axes[2].loglog(
            betas,
            [row["efficiency_vs_optimized"] for row in selected],
            "o-",
            label=label,
        )
axes[0].set_ylabel("RMS reported relative error (%)")
axes[1].set_ylabel("Mean wall time (seconds)")
axes[2].set_ylabel("Efficiency vs ordinary LMBs (higher is better)")
axes[2].axhline(1, color="gray", linestyle="--", linewidth=0.8)
if not zero_temperature:
    axes[0].legend(fontsize=8)
for axis in axes:
    if zero_temperature:
        present = [
            mode for mode in mode_labels if any(row["mode"] == mode for row in rows)
        ]
        axis.set_xticks(
            [list(mode_labels).index(mode) for mode in present],
            [
                {
                    "optimized": "Ordinary",
                    "optimized_fermi": "+ Fermi shells",
                    "optimized_fermi_repeated": "+ Repeated shell",
                }[mode]
                for mode in present
            ],
            rotation=15,
        )
        axis.set_axisbelow(True)
    else:
        axis.set_xlabel(r"$\beta\mu$ ($\mu=1$)")
    axis.grid(alpha=0.2, which="both")
sample_counts = sorted({row["samples_per_run"] for row in rows})
sample_label = (
    f"{sample_counts[0]:,}"
    if len(sample_counts) == 1
    else f"{sample_counts[0]:,}–{sample_counts[-1]:,}"
)
precision_label = (
    "Arb only"
    if data["metadata"].get("arb_only")
    else "Double → Arb"
    if data["metadata"].get("skip_quad")
    else "Double → Quad → Arb"
)
fig.suptitle(
    f"{'Exact T=0' if zero_temperature else 'Thermal'} {data['metadata'].get('example', 'sunrise').replace('_', ' ')} · {sample_label} samples/run · "
    f"{rows[0]['seeds']} seeds · one worker"
    + (f"\n{precision_label}" if zero_temperature else "")
    + (
        "\nPRECISION FAILURE: efficiency is not interpretable"
        if data["metadata"].get("validation_note")
        else ""
    )
)
fig.savefig(args.output / "comparison.png", dpi=180)
fig.savefig(args.output / "comparison.pdf")
for row in rows:
    temperature_label = "T=0" if zero_temperature else f"{row['beta']:6g}"
    print(
        f"{temperature_label} {row['mode']:16s} error={row['relative_error_percent']:.3f}% "
        f"wall={row['wall_seconds']:.2f}s efficiency={row['efficiency_vs_optimized']:.2f}x "
        f"pull={row['mean_pull']} rejected={row['max_rejected_percent']:.3g}%"
    )
