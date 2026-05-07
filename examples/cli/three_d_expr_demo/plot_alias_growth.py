#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def read_rows(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in [
            "j",
            "atom_ops",
            "atom_bytes",
            "aliased_ops",
            "aliased_bytes",
            "ops_ratio",
            "bytes_ratio",
            "alias_count",
        ]:
            row[key] = float(row[key]) if row[key] else float("nan")
    return rows


def case_title(path, rows):
    first = rows[0]
    graph = first["graph_name"]
    rep = first["representation"].upper()
    numerator = "with numerator" if first["with_numerator"] == "true" else "without numerator"
    return f"{graph} {rep} {numerator} ({Path(path).name})"


def set_dynamic_log_ylim(axis, *series):
    values = [
        value
        for data in series
        for value in data
        if math.isfinite(value) and value > 0.0
    ]
    if not values:
        return

    low = min(values)
    high = max(values)
    if low == high:
        axis.set_ylim(low / 1.5, high * 1.5)
        return

    pad = (math.log10(high) - math.log10(low)) * 0.08
    axis.set_ylim(10 ** (math.log10(low) - pad), 10 ** (math.log10(high) + pad))


def plot_case(pdf, path):
    rows = read_rows(path)
    if not rows:
        return

    x = [row["j"] for row in rows]
    atom_ops = [row["atom_ops"] for row in rows]
    aliased_ops = [row["aliased_ops"] for row in rows]
    atom_bytes = [row["atom_bytes"] for row in rows]
    aliased_bytes = [row["aliased_bytes"] for row in rows]
    ops_ratio = [row["ops_ratio"] for row in rows]
    bytes_ratio = [row["bytes_ratio"] for row in rows]

    fig, (top, bottom) = plt.subplots(
        2,
        1,
        figsize=(11.0, 8.2),
        gridspec_kw={"height_ratios": [2.2, 1.0]},
        constrained_layout=True,
    )
    fig.suptitle(case_title(path, rows), fontsize=14, fontweight="bold")

    top_bytes = top.twinx()
    top.plot(x, atom_ops, color="#8f2d2d", linewidth=2.0, label="Atom ops")
    top.plot(
        x,
        aliased_ops,
        color="#8f2d2d",
        linewidth=2.0,
        linestyle="--",
        label="Aliased ops",
    )
    top_bytes.plot(x, atom_bytes, color="#245c9c", linewidth=2.0, label="Atom bytes")
    top_bytes.plot(
        x,
        aliased_bytes,
        color="#245c9c",
        linewidth=2.0,
        linestyle="--",
        label="Aliased bytes",
    )
    top.set_ylabel("operations")
    top_bytes.set_ylabel("bytes")
    top.set_yscale("log")
    top_bytes.set_yscale("log")
    top.grid(True, which="both", alpha=0.25)

    handles, labels = top.get_legend_handles_labels()
    handles_b, labels_b = top_bytes.get_legend_handles_labels()
    top.legend(handles + handles_b, labels + labels_b, loc="upper left", ncols=2)

    bottom.plot(x, ops_ratio, color="#8f2d2d", linewidth=2.0, label="ops ratio")
    bottom.plot(x, bytes_ratio, color="#245c9c", linewidth=2.0, label="bytes ratio")
    bottom.axhline(1.0, color="black", linewidth=1.0, alpha=0.4)
    bottom.set_xlabel("partial sum size j")
    bottom.set_ylabel("Atom / AliasedAtom")
    bottom.set_yscale("log")
    set_dynamic_log_ylim(bottom, ops_ratio, bytes_ratio)
    bottom.grid(True, which="both", alpha=0.25)
    bottom.legend(loc="upper left")

    footer = rows[-1]
    fig.text(
        0.01,
        0.01,
        "final: "
        f"j={int(footer['j'])}, "
        f"aliases={int(footer['alias_count'])}, "
        f"ops ratio={footer['ops_ratio']:.3g}, "
        f"bytes ratio={footer['bytes_ratio']:.3g}",
        fontsize=9,
    )
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="+", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output) as pdf:
        for path in args.csv:
            plot_case(pdf, path)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
