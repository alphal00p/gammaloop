"""Analyze the frozen-grid convergence campaign of campaign.py.

Cumulative per-iteration archives give the sums S1 = n*mean and
S2 = n*((n-1)*error^2 + mean^2) of all samples so far. Differences with the last training
iteration isolate the frozen samples, whose mean and error do not depend on adaptation.
Disjoint blocks of frozen iterations give independent estimates at block sizes N = 2^k
iterations. For every configuration the script compares the empirical spread of these
block estimates with their internal errors, fits both against N^b, takes the slope
uncertainty from a bootstrap over seeds, and plots them with the cumulative estimate of
every seed.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

REFERENCE_FILE = Path(__file__).resolve().with_name("reference.json")
# Independent references (value, uncertainty) for Gamma_rho at T=0.5 (benchmark
# normalization), used when REFERENCE_FILE is absent.
REFERENCES = {
    "dt_rest1": (3.057651940618655e-4, 0.0),
    "dt_rest16": (5.282087736915636e-6, 0.0),
}


def frozen_increments(workspace, training):
    """Sample count and sums (n, S1, S2) of every frozen iteration."""
    archives = sorted((workspace / "results").glob("integration_result_iter_*.json"))
    sums = [[0.0, 0.0, 0.0]]
    for archive in archives:
        integral = json.loads(archive.read_text())["slots"][0]["integral"]
        n, mean = integral["neval"], integral["result"]["re"]
        error = integral["error"]["re"]
        sums.append([n, n * mean, n * ((n - 1) * error * error + mean * mean)])
    return np.diff(sums[training:], axis=0)


def estimate(increments):
    """Mean and error of the pooled samples of some iterations."""
    n, s1, s2 = np.sum(increments, axis=0)
    mean = s1 / n
    return mean, math.sqrt(max(s2 / n - mean * mean, 0.0) / (n - 1))


def blocks(seeds, size):
    """Estimates (mean, error) of disjoint blocks of `size` frozen iterations, per seed."""
    return [
        [
            estimate(increments[start : start + size])
            for start in range(0, len(increments) - size + 1, size)
        ]
        for increments in seeds
    ]


def fit_power(n, sigma, sigma_error=None):
    """Weighted least squares of log(sigma) = a + b log(N); returns b, its error and a."""
    x, y = np.log(n), np.log(sigma)
    w = (
        np.ones_like(x)
        if sigma_error is None
        else (sigma / np.maximum(sigma_error, 1e-300)) ** 2
    )
    design = np.vstack([np.ones_like(x), x]).T
    covariance = np.linalg.inv(design.T @ (design * w[:, None]))
    a, b = covariance @ design.T @ (w * y)
    residual = y - (a + b * x)
    dof = max(len(x) - 2, 1)
    chi2 = float(residual @ (w * residual)) / dof
    scale = chi2 if sigma_error is None else max(chi2, 1.0)
    return float(b), float(math.sqrt(covariance[1, 1] * scale)), float(a)


def block_scaling(seeds, batch):
    """Spread and median internal error of disjoint block estimates per block size."""
    length = min(len(increments) for increments in seeds)
    rows = []
    size = 1
    while size <= length // 2:
        estimates = [value for per_seed in blocks(seeds, size) for value in per_seed]
        means = np.array([mean for mean, _ in estimates])
        errors = np.array([error for _, error in estimates])
        spread = means.std(ddof=1)
        rows.append(
            {
                "n": size * batch,
                "blocks": len(estimates),
                "spread": float(spread),
                "spread_error": float(spread / math.sqrt(2 * (len(estimates) - 1))),
                "reported": float(np.median(errors)),
            }
        )
        size *= 2
    return rows


def columns(scaling, *keys):
    return [np.array([row[key] for row in scaling], dtype=float) for key in keys]


def spread_slope(scaling):
    """Weighted power-law fit of the block spread against the block size."""
    return fit_power(*columns(scaling, "n", "spread", "spread_error"))


def analyze(directory, training, references, bootstrap, rng, label=None):
    seeds = [
        frozen_increments(seed / "workspace", training)
        for seed in sorted(directory.glob("seed*"))
    ]
    seeds = [increments for increments in seeds if len(increments) >= 4]
    if len(seeds) < 2:
        return None
    length = min(len(increments) for increments in seeds)
    seeds = [increments[:length] for increments in seeds]
    batch = int(seeds[0][0][0])
    scaling = block_scaling(seeds, batch)
    n, spread, reported = columns(scaling, "n", "spread", "reported")
    slope_spread, slope_spread_error, _ = spread_slope(scaling)
    slope_reported, slope_reported_error, _ = fit_power(n, reported)
    resampled = []
    for _ in range(bootstrap):
        sample = [seeds[index] for index in rng.integers(0, len(seeds), len(seeds))]
        rows = block_scaling(sample, batch)
        if all(row["spread"] > 0 for row in rows):
            resampled.append(spread_slope(rows)[0])
    finals = [estimate(increments) for increments in seeds]
    final_mean = np.array([mean for mean, _ in finals])
    final_error = np.array([error for _, error in finals])
    combined_mean = float(np.average(final_mean, weights=final_error**-2))
    combined_error = float(np.sum(final_error**-2) ** -0.5)
    unweighted_mean = float(final_mean.mean())
    unweighted_error = float(final_mean.std(ddof=1) / math.sqrt(len(seeds)))
    summary = {
        "configuration": label or directory.name,
        "replicas": len(seeds),
        "frozen_samples_per_replica": int(np.sum(seeds[0], axis=0)[0]),
        "combined_estimate": combined_mean,
        "combined_error": combined_error,
        "unweighted_mean": unweighted_mean,
        "unweighted_mean_error": unweighted_error,
        "replica_chi2_per_dof": float(
            np.sum(((final_mean - combined_mean) / final_error) ** 2) / (len(seeds) - 1)
        ),
        "block_spread_over_reported": [float(value) for value in spread / reported],
        "block_spread_slope": [slope_spread, slope_spread_error],
        "block_spread_slope_bootstrap_68": (
            [float(value) for value in np.percentile(resampled, [16, 84])]
            if resampled
            else None
        ),
        "reported_error_slope": [slope_reported, slope_reported_error],
        "block_scaling": scaling,
        "replicas_final": [
            [float(m), float(e)] for m, e in zip(final_mean, final_error)
        ],
    }
    reference = references.get(directory.name.rsplit("_", 1)[0])
    if reference is not None:
        value, uncertainty = reference
        summary["reference"] = value
        summary["reference_error"] = uncertainty
        summary["combined_pull"] = (combined_mean - value) / math.hypot(
            combined_error, uncertainty
        )
        summary["unweighted_pull"] = (unweighted_mean - value) / math.hypot(
            unweighted_error, uncertainty
        )
    ladders = [
        [estimate(increments[: m + 1]) for m in range(length)] for increments in seeds
    ]
    return summary, ladders


def plot(output, result):
    summary, ladders = result
    figure = Figure(figsize=(11, 4.2), layout="tight")
    left, right = figure.subplots(1, 2)
    scaling = summary["block_scaling"]
    n, spread, spread_error, reported = columns(
        scaling, "n", "spread", "spread_error", "reported"
    )
    left.loglog(n, reported, "s-", color="C0", label="median internal error of blocks")
    left.errorbar(
        n,
        spread,
        yerr=spread_error,
        fmt="o",
        color="C1",
        label="spread of disjoint block estimates",
    )
    left.loglog(
        n, reported[0] * np.sqrt(n[0] / n), "k--", lw=1, label=r"$\propto N^{-1/2}$"
    )
    b, db = summary["block_spread_slope"]
    bootstrap = summary["block_spread_slope_bootstrap_68"]
    interval = f" [{bootstrap[0]:.2f}, {bootstrap[1]:.2f}]" if bootstrap else ""
    ratio = summary["block_spread_over_reported"]
    left.set_title(
        f"{summary['configuration']}: slope {b:.3f} ± {db:.3f}{interval}\n"
        f"spread / internal error {ratio[0]:.2f} … {ratio[-1]:.2f}",
        fontsize=10,
    )
    left.set_xlim(n[0] / 1.5, n[-1] * 1.5)
    left.set_xlabel("samples per block N")
    left.set_ylabel(r"uncertainty of $\Gamma_\rho$")
    left.legend(fontsize=8)
    batch = n[0]
    for ladder in ladders:
        counts = batch * np.arange(1, len(ladder) + 1)
        means = np.array([mean for mean, _ in ladder])
        errors = np.array([error for _, error in ladder])
        right.semilogx(counts, means, color="0.5", lw=0.8)
        right.fill_between(
            counts, means - errors, means + errors, color="C0", alpha=0.06
        )
    if "reference" in summary:
        value, uncertainty = summary["reference"], summary["reference_error"]
        right.axhline(value, color="C3", label="independent reference")
        right.axhspan(value - uncertainty, value + uncertainty, color="C3", alpha=0.2)
        right.legend(fontsize=8)
    right.set_xlabel("cumulative frozen samples per seed")
    right.set_ylabel(r"$\Gamma_\rho$ estimate ± reported error")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"{summary['configuration']}.{suffix}", dpi=150)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--training", type=int, default=10)
    parser.add_argument(
        "--references",
        type=Path,
        help="JSON object mapping <process>_<kinematics> to [value, uncertainty] "
        f"(default: {REFERENCE_FILE.name} next to this script, else built-in values)",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bootstrap", type=int, default=400)
    parser.add_argument(
        "--label", help="configuration name replacing the directory name"
    )
    args = parser.parse_args()
    directories = sorted(
        path for path in args.campaign.iterdir() if (path / "seed0").exists()
    )
    if args.label and len(directories) > 1:
        parser.error("--label requires a campaign with a single configuration")
    path = args.references or (REFERENCE_FILE if REFERENCE_FILE.exists() else None)
    references = (
        REFERENCES
        if path is None
        else {key: tuple(pair) for key, pair in json.loads(path.read_text()).items()}
    )
    output = args.output or args.campaign / "analysis"
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261102)
    summaries = []
    for directory in directories:
        result = analyze(
            directory, args.training, references, args.bootstrap, rng, args.label
        )
        if result is None:
            continue
        summaries.append(result[0])
        plot(output, result)
    (output / "summary.json").write_text(json.dumps(summaries, indent=1))
    for summary in summaries:
        pull = summary.get("combined_pull")
        ratio = summary["block_spread_over_reported"]
        slope, slope_error = summary["block_spread_slope"]
        print(
            f"{summary['configuration']:24s} {summary['combined_estimate']:+.6e}"
            f" ± {summary['combined_error']:.2e}"
            f"  block spread/reported {ratio[0]:.2f}..{ratio[-1]:.2f}"
            f"  spread slope {slope:+.3f}±{slope_error:.3f}"
            f" {summary['block_spread_slope_bootstrap_68']}"
            + (f"  pull {pull:+.2f}" if pull is not None else "")
        )


if __name__ == "__main__":
    main()
