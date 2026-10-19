"""Render saved scan results in Figure 15 conventions; never start integrations."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run_root", type=Path)
parser.add_argument(
    "--published-curves",
    type=Path,
    default=Path(__file__).resolve().parents[1] / "reference/published-curves.json",
)
args = parser.parse_args()
root = args.run_root.resolve()
os.environ.setdefault("MPLCONFIGDIR", str(root / "matplotlib-cache"))
for thread_variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[thread_variable] = "1"
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

manifest = json.loads((root / "manifest.json").read_text())
rows = json.loads((root / "summary.json").read_text())
references = json.loads((root / "reference/targets.json").read_text())
published = json.loads(args.published_curves.read_text())
lookup = {(round(r["E"], 9), round(r["k"], 9), r["component"]): r for r in references}
table = []
for raw in rows:
    row = dict(raw)
    reference = lookup.get(
        (round(row["k0_over_T"], 9), round(row["k_over_T"], 9), row["component"])
    )
    if reference and row.get("value") is not None and row.get("error") is not None:
        lo, nlo = reference["rho_lo"], reference["rho_nlo"]
        row.update(
            analytic_lo=lo,
            reference_nlo=nlo,
            reference_g3=lo + 9 * nlo,
            analytic_lo_plus_native_nlo_g3=lo + 9 * row["value"],
            error_g3=9 * row["error"],
            delta_nlo_relative=(row["value"] - nlo) / abs(nlo) if nlo else None,
            reference_tolerance=reference.get("reference_tolerance"),
        )
    table.append(row)
columns = [
    "name",
    "component",
    "k0_over_T",
    "k_over_T",
    "samples",
    "value",
    "error",
    "reference_nlo",
    "delta_nlo_relative",
    "pull",
    "relative_error_to_reference",
    "precision_reached",
    "final_stability_clean",
    "chi_sq_per_dof",
    "max_weight_impact",
    "analytic_lo",
    "reference_g3",
    "analytic_lo_plus_native_nlo_g3",
    "error_g3",
    "status",
]
with (root / "figure15-comparison.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(table)
summary = {
    "status": manifest.get("status", "running"),
    "source_commit": manifest["source_commit"],
    "started_utc": manifest["started_utc"],
    "hard_deadline_utc": manifest["hard_deadline_utc"],
    "available_components": len(table),
    "components_with_reference": sum("reference_nlo" in r for r in table),
    "precision_reached": sum(bool(r["precision_reached"]) for r in table),
    "within_half_percent_of_target": sum(
        r.get("delta_nlo_relative") is not None
        and abs(r["delta_nlo_relative"]) <= 0.005
        for r in table
    ),
    "discrepancy_watch": [
        {
            "name": r["name"],
            "component": r["component"],
            "pull": r["pull"],
            "samples": r["samples"],
            "chi_sq_per_dof": r.get("chi_sq_per_dof"),
            "max_weight_impact": r.get("max_weight_impact"),
        }
        for r in table
        if r.get("pull") is not None and abs(r["pull"]) >= 3
    ],
    "numerical_failures": [
        {"name": r["name"], "component": r["component"]}
        for r in table
        if not r["final_stability_clean"]
    ],
    "interpretation": "Native isolated NLO is compared independently. Figure overlay uses author analytic LO + 9*native NLO; the LO curve is not a new native integration. A three-sigma residual is a watch flag, not proof of a systematic discrepancy.",
}
(root / "comparison-report.json").write_text(json.dumps(summary, indent=2) + "\n")
colors = {0.0: "#c57915", 0.5: "#c12636", 1.0: "#8246a5", 1.5: "#3155ac"}
for zoom in (False, True):
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, layout="constrained")
    for ax, component in zip(axes, ("T", "L")):
        for curve in published["curves"]:
            if curve["component"] != component:
                continue
            is_lo = curve["order"] == "lo"
            k = curve["k_over_T"]
            for segment_index, segment in enumerate(curve["segments"]):
                ax.plot(
                    [vertex[0] for vertex in segment],
                    [vertex[1] for vertex in segment],
                    color="#999999" if is_lo else colors[k],
                    linewidth=0.8 if is_lo else 1.3,
                    linestyle="--" if is_lo else ":" if k == 0 else "-",
                    label=f"Published k/T={k:g}"
                    if not is_lo and segment_index == 0
                    else None,
                )
        for k, color in colors.items():
            native = sorted(
                (
                    r
                    for r in table
                    if r["component"] == component
                    and r["k_over_T"] == k
                    and "reference_g3" in r
                ),
                key=lambda r: r["k0_over_T"],
            )
            if native:
                ax.errorbar(
                    [r["k0_over_T"] for r in native],
                    [r["analytic_lo_plus_native_nlo_g3"] for r in native],
                    yerr=[r["error_g3"] for r in native],
                    fmt="o",
                    markersize=3,
                    capsize=2,
                    color=color,
                )
        ax.set_ylabel(r"$\rho_" + component + r"/(eT)^2$")
        ax.axhline(0, color="gray", linewidth=0.6)
        ax.grid(alpha=0.15)
        ax.set_xlim(0, 2.03)
    if zoom:
        axes[0].set_ylim(-0.4, 1.5)
        axes[1].set_ylim(-1.2, 0.7)
        fig.supxlabel("E/T; fixed display range, see CSV for complete error bars")
    else:
        axes[1].set_xlabel("E/T")
    axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle(
        "Figure 15 comparison, g=3\nLines: published curves; gray dashed: LO\nMarkers: analytic LO + 9 × native NLO, with MC error bars"
    )
    stem = "figure15-overlay-zoom" if zoom else "figure15-overlay"
    fig.savefig(root / (stem + ".png"), dpi=170)
    fig.savefig(root / (stem + ".pdf"))
    plt.close(fig)
fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
for component, marker in [("T", "o"), ("L", "s")]:
    for k, color in colors.items():
        selected = [
            r
            for r in table
            if r["component"] == component
            and r["k_over_T"] == k
            and r.get("pull") is not None
        ]
        if selected:
            ax.scatter(
                [r["k0_over_T"] for r in selected],
                [r["pull"] for r in selected],
                marker=marker,
                color=color,
                label=f"{component}, k/T={k:g}",
            )
ax.axhspan(-3, 3, color="gray", alpha=0.1)
ax.axhline(0, color="black", linewidth=0.7)
ax.set(
    xlabel="E/T",
    ylabel="(native NLO − author NLO) / MC error",
    title="Isolated NLO residuals; adaptive-error and rare-tail diagnostics remain in CSV",
)
ax.legend(fontsize=7, ncol=4)
ax.grid(alpha=0.15)
fig.savefig(root / "nlo-residuals.png", dpi=170)
plt.close(fig)
(root / "figure15-render.json").write_text(
    json.dumps(
        {
            "curve_source": published["source"],
            "curve_file_sha256": hashlib.file_digest(
                args.published_curves.open("rb"), "sha256"
            ).hexdigest(),
            "curve_method": published["method"],
            "curve_limitation": published["limitation"],
            "native_markers": [
                {
                    "name": row["name"],
                    "component": row["component"],
                    "E_over_T": row["k0_over_T"],
                    "k_over_T": row["k_over_T"],
                    "value_g3": row["analytic_lo_plus_native_nlo_g3"],
                    "error_g3": row["error_g3"],
                }
                for row in table
                if "reference_g3" in row
            ],
            "zoom_ylim": {"T": [-0.4, 1.5], "L": [-1.2, 0.7]},
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(summary))
