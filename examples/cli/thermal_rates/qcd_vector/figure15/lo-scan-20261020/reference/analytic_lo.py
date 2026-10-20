#!/usr/bin/env python3
"""Independent massless one-loop thermal vector spectral functions.

Evaluate only smooth one-dimensional Fermi phase-space moments.  No NLO
master, reference binary, GSL, or third-party Python dependency is used.
Nc=3, T=1, one unit-charge flavor, with the electromagnetic charge factored out.
"""

from __future__ import annotations

import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD_TARGETS = Path("/tmp/qcd-fig15-scan-20261018T171925Z/reference/targets.json")
AUTHOR = Path("/tmp/qcd-fig15-reference-20261018/portable-build/gw3g-spectral-dca510d")
ARCHIVE = AUTHOR.parent / "spectral-v0.9-spf.zip"
NC = 3.0
ARCHIVE_SHA = "207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb"


@lru_cache(maxsize=None)
def gauss_legendre(n: int) -> tuple[tuple[float, float], ...]:
    """Compute the standard n-node Gauss-Legendre rule on [-1,1]."""
    nodes = []
    for i in range(1, (n + 1) // 2 + 1):
        z = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(30):
            previous, polynomial = 1.0, z
            for j in range(2, n + 1):
                previous, polynomial = polynomial, (
                    (2 * j - 1) * z * polynomial - (j - 1) * previous
                ) / j
            derivative = n * (z * polynomial - previous) / (z * z - 1)
            updated = z - polynomial / derivative
            if abs(updated - z) < 2e-16:
                z = updated
                break
            z = updated
        # Reevaluate the derivative at the final root.
        previous, polynomial = 1.0, z
        for j in range(2, n + 1):
            previous, polynomial = polynomial, (
                (2 * j - 1) * z * polynomial - (j - 1) * previous
            ) / j
        derivative = n * (z * polynomial - previous) / (z * z - 1)
        weight = 2 / ((1 - z * z) * derivative * derivative)
        nodes.append((-z, weight))
        if z != 0:
            nodes.append((z, weight))
    return tuple(sorted(nodes))


def spectral_lo(E: float, k: float, nodes: int = 64) -> dict[str, float]:
    """Return physical rho_T, rho_L, and rho_V=2rho_T+rho_L.

    Timelike: r=(E+k*x)/2 and W=1-nF(r)-nF(E-r).  Symmetry of
    W under x->-x gives kernels (1+x²) and (1-x²).
    Spacelike: p>=(k-E)/2, D=nF(p)-nF(p+E), t=2p+E.  The two
    crossed cuts give kernels (t²+k²) and (t²-k²), respectively.
    These are algebraically the author's psi0/1/2 LO moments, integrated
    directly to avoid the pinned source's truncated Li3 implementation.
    """
    if E < 0 or k < 0:
        raise ValueError("This reference covers nonnegative E and k only")
    if E == 0 or E == k:
        return {"T": 0.0, "L": 0.0, "V": 0.0, "tail_bound": 0.0}
    if k == 0:
        value = NC * E * E * math.tanh(E / 4) / (12 * math.pi)
        return {"T": value, "L": value, "V": 3 * value, "tail_bound": 0.0}
    K2 = (E - k) * (E + k)
    if E > k:
        tterms, lterms = [], []
        for x, weight in gauss_legendre(nodes):
            r = (E + k * x) / 2
            thermal = -math.expm1(-E) / (
                (1 + math.exp(-r)) * (1 + math.exp(-(E - r)))
            )
            tterms.append(weight * (1 + x * x) * thermal)
            lterms.append(weight * (1 - x * x) * thermal)
        rho_t = NC * K2 / (32 * math.pi) * math.fsum(tterms)
        rho_l = NC * K2 / (16 * math.pi) * math.fsum(lterms)
        tail = 0.0
    else:
        a = (k - E) / 2
        tterms, lterms = [], []
        edges = (0, 1, 2, 4, 8, 16, 32, 64)
        for left, right in zip(edges, edges[1:]):
            half = (right - left) / 2
            middle = (right + left) / 2
            for x, weight in gauss_legendre(nodes):
                p = a + middle + half * x
                ep = math.exp(-p)
                difference = ep * -math.expm1(-E) / (
                    (1 + ep) * (1 + ep * math.exp(-E))
                )
                t = 2 * p + E
                tterms.append(half * weight * (t * t + k * k) * difference)
                lterms.append(half * weight * (t * t - k * k) * difference)
        rho_t = -NC * K2 / (8 * math.pi * k**3) * math.fsum(tterms)
        rho_l = NC * K2 / (4 * math.pi * k**3) * math.fsum(lterms)
        # D <= exp(-p)*(1-exp(-E)); bound both omitted polynomial tails.
        pmax = a + edges[-1]
        integral_bound = math.exp(-pmax) * -math.expm1(-E) * (
            4 * (pmax * pmax + 2 * pmax + 2)
            + 4 * E * (pmax + 1) + E * E + k * k
        )
        tail = abs(NC * K2 / (4 * math.pi * k**3)) * integral_bound
    return {"T": rho_t, "L": rho_l, "V": 2 * rho_t + rho_l, "tail_bound": tail}


def trace_closed(E: float, k: float) -> float:
    if k == 0:
        return NC * E * E * math.tanh(E / 4) / (4 * math.pi)
    kp, km = (E + k) / 2, (E - k) / 2
    psi0 = float(km > 0) + 2 / k * (
        math.log1p(math.exp(-kp)) - math.log1p(math.exp(-abs(km)))
    )
    return NC * (E - k) * (E + k) * psi0 / (4 * math.pi)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name: str, value: object) -> None:
    (HERE / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main() -> None:
    assert sha256(ARCHIVE) == ARCHIVE_SHA
    grid = [(i / 4, k) for k in (0.0, 0.5, 1.0, 1.5)
            for i in range(1, 9) if i / 4 != k]
    grid += [(E, k) for E, k in (
        (0.45, 0.5), (0.55, 0.5), (0.95, 1.0),
        (1.05, 1.0), (1.45, 1.5), (1.55, 1.5))]
    grid = sorted(grid, key=lambda point: (point[1], point[0]))
    assert len(grid) == 35 and len(set(grid)) == 35
    rows, refinements, traces = [], [], []
    for E, k in grid:
        coarse, fine = spectral_lo(E, k, 64), spectral_lo(E, k, 96)
        for component in ("T", "L"):
            value = fine[component]
            numerical_bound = max(
                abs(value - coarse[component]),
                16 * math.ulp(value),
            ) + fine["tail_bound"]
            rows.append({
                "E": E, "k": k, "component": component,
                "rho_lo": value, "value": value,
                "reference_abs_error_estimate": numerical_bound,
                "source": "independent analytic LO Fermi phase-space moments",
                "source_archive_sha256": ARCHIVE_SHA,
            })
            refinements.append({
                "E": E, "k": k, "component": component,
                "rho_lo_64": coarse[component], "rho_lo_96": value,
                "absolute_difference": abs(value - coarse[component]),
                "relative_difference": abs(value - coarse[component]) / abs(value),
                "tail_bound": fine["tail_bound"],
            })
        closed = trace_closed(E, k)
        traces.append({"E": E, "k": k, "from_components": fine["V"],
                       "closed": closed, "absolute_difference": abs(fine["V"] - closed),
                       "relative_difference": abs(fine["V"] / closed - 1)})
    by_key = {(row["E"], row["k"], row["component"]): row for row in rows}
    prior = json.loads(OLD_TARGETS.read_text())
    comparisons = []
    for old in prior:
        new = by_key[(old["E"], old["k"], old["component"])]["rho_lo"]
        comparisons.append({
            "E": old["E"], "k": old["k"], "component": old["component"],
            "old_author_rho_lo": old["rho_lo"], "analytic_rho_lo": new,
            "absolute_difference": abs(new - old["rho_lo"]),
            "relative_difference": abs(new / old["rho_lo"] - 1),
        })
    assert len(prior) == len(comparisons) == 64
    # Legacy author special-function approximation is not the accuracy standard.
    assert max(row["relative_difference"] for row in comparisons) < 1e-5
    assert max(row["relative_difference"] for row in refinements) < 1e-12
    assert max(row["relative_difference"] for row in traces) < 1e-12
    assert all(row["rho_lo"] > 0 for row in rows if row["component"] == "T")
    assert all(row["rho_lo"] < 0 for row in rows
               if row["component"] == "L" and row["E"] < row["k"])
    assert all(row["rho_lo"] > 0 for row in rows if row["E"] > row["k"])
    rest_checks = []
    for E, k in grid:
        if k != 0:
            continue
        coarse, fine = spectral_lo(E, 0.0002, 96), spectral_lo(E, 0.0001, 96)
        exact = spectral_lo(E, 0.0)["T"]
        for component in ("T", "L"):
            extrapolated = (4 * fine[component] - coarse[component]) / 3
            rest_checks.append({"E": E, "component": component,
                                "exact_rest": exact,
                                "small_k_squared_extrapolation": extrapolated,
                                "relative_difference": abs(extrapolated / exact - 1)})
    assert max(row["relative_difference"] for row in rest_checks) < 1e-12
    provenance = {
        "normalization": {"T": 1, "Nc": 3, "flavors": 1,
                          "electric_charge_squared_factored_out": True,
                          "flavor_charge_squared_factored_out": True,
                          "rhoV_relation": "rhoV=2*rhoT+rhoL",
                          "order": "LO; no g^2 multiplier"},
        "source": {"paper": "https://arxiv.org/abs/1910.07552v2",
                   "source_archive_url": "https://zenodo.org/api/records/3478144/files/gw3g/spectral-v0.9-spf.zip/content",
                   "source_doi": "https://doi.org/10.5281/zenodo.3478144",
                   "source_tag": "v0.9-spf",
                   "source_commit": "dca510d2d2c8bee2288345e972ee3a54b37cdbd5",
                   "source_archive_sha256": ARCHIVE_SHA,
                   "archive_hash_verified_against_local_bytes": True,
                   "author_lo_expression_locations": [
                       "examples/dileptons.cpp:76", "examples/dileptons.cpp:164",
                       "src/therm.cpp:61-132"],
                   "source_file_sha256": {
                       str(path.relative_to(AUTHOR)): sha256(path) for path in (
                           AUTHOR / "examples/dileptons.cpp", AUTHOR / "src/therm.cpp")}},
        "method": "Direct smooth one-dimensional massless Fermi phase-space moments; Gauss-Legendre 64 vs96 nodes; exact analytic rest and zero/lightcone limits",
        "rest_formula": "rhoT=rhoL=Nc*E^2*tanh(E/4)/(12*pi)",
        "spacelike_convention": "rhoT>0,rhoL<0; ImPi00=-author_Rho00.lo; rhoL=-(E^2-k^2)/k^2*author_Rho00.lo",
        "precision": "binary64; refinement difference, 16 ulps, and an explicit exponential tail bound give per-target numerical estimates, not rigorous total roundoff bounds",
        "upstream_difference": "Author Li3 approximation is truncated; tiny differences from previous author LO values are expected and retained in validation.json",
        "old_targets_path": str(OLD_TARGETS), "old_targets_sha256": sha256(OLD_TARGETS),
        "generator_sha256": sha256(Path(__file__)),
    }
    write_json("targets.json", rows)
    write_json("provenance.json", provenance)
    validation = {
        "points": 35, "targets": len(rows), "old_targets_compared": len(comparisons),
        "max_old_relative_difference": max(r["relative_difference"] for r in comparisons),
        "max_old_absolute_difference": max(r["absolute_difference"] for r in comparisons),
        "max_refinement_relative_difference": max(r["relative_difference"] for r in refinements),
        "max_trace_relative_difference": max(r["relative_difference"] for r in traces),
        "sign_checks_passed": True, "exact_rest_used": True,
        "max_rest_limit_relative_difference": max(r["relative_difference"] for r in rest_checks),
        "previously_missing_points": [
            {"E": E, "k": k, **{c: spectral_lo(E, k, 96)[c] for c in ("T", "L", "V")}}
            for E, k in ((0.25, 1.0), (0.25, 1.5), (2.0, 0.0))],
        "old_target_comparisons": comparisons,
        "refinements": refinements, "trace_checks": traces, "rest_limit_checks": rest_checks,
    }
    write_json("validation.json", validation)
    print(json.dumps({k: v for k, v in validation.items()
                      if k not in ("old_target_comparisons", "refinements", "trace_checks", "rest_limit_checks")}, indent=2))


if __name__ == "__main__":
    main()
