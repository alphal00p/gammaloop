"""Independent finite-temperature reference for Gamma_rho(P) = -Im Pi_R(P)/E.

Massive scalar double triangle (dt) and triangle-box-triangle (tbt) self-energies with
the interaction -lam/2 (X H^2 + L H^2 + X L^2), lam=1, mH=0.6, mL=0.25, mu=0, Bose
statistics.  Pi_R = -kappa I_R with I_R the retarded continuation of the positive
Euclidean master sum-integral (unit cubic couplings), so Gamma = kappa Im I_R(E+i0)/E.

Method (see derivation.typ): exact Matsubara sums by iterated residues (matsubara.py),
external frequency kept bosonic until the end, continuation z -> E + i eps, radial
quadrature sinh-mapped around every eps-broadened threshold root (engine.py), vegas over
the remaining angles, polynomial extrapolation eps -> 0 with the full MC covariance.
No GammaLoop code, expression or evaluator is used.

Usage:  python thermal_reference.py --graph dt --energy 1 --momentum 0 --temperature 0.5 \
            --output out.json            (temperature 0 gives the vacuum)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import gvar
import numpy as np
import vegas

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from engine import Problem

# sha256 of the engine sources, hashed at import and recorded in every output JSON
CODE_SHA256 = {
    name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
    for name in ("thermal_reference.py", "engine.py", "matsubara.py")
}

_CACHE: dict = {}


class Integrand:
    """Picklable vegas lbatch integrand; each worker builds its Problem once."""

    fcntype = "lbatch"

    def __init__(self, kwargs, chunk=32):
        self.kwargs = kwargs
        self.chunk = chunk

    def problem(self):
        key = json.dumps(self.kwargs, sort_keys=True)
        if key not in _CACHE:
            _CACHE.clear()
            _CACHE[key] = Problem(**self.kwargs)
        return _CACHE[key]

    def __call__(self, u):
        import warnings

        warnings.simplefilter("ignore")
        p = self.problem()
        out = [p(u[i : i + self.chunk]) for i in range(0, u.shape[0], self.chunk)]
        return np.concatenate(out, axis=0)


def poly_weights(eps, degree):
    """Rows of the least-squares map from J(eps_i) to polynomial coefficients."""
    V = np.vander(np.asarray(eps), degree + 1, increasing=True)
    return np.linalg.pinv(V)


def gls_fit(eps, mean, cov, degree):
    """Correlated least squares J(eps) = sum_k c_k eps^k; returns c0, sigma(c0), chi2, dof."""
    eps = np.asarray(eps)
    V = np.vander(eps, degree + 1, increasing=True)
    # regularise the (nearly singular, strongly correlated) covariance mildly
    covr = cov + np.eye(len(eps)) * 1e-12 * np.max(np.diag(cov))
    ci = np.linalg.inv(covr)
    A = V.T @ ci @ V
    Ai = np.linalg.inv(A)
    c = Ai @ V.T @ ci @ mean
    r = mean - V @ c
    chi2 = float(r @ ci @ r)
    return (
        float(c[0]),
        float(np.sqrt(Ai[0, 0])),
        chi2,
        len(eps) - degree - 1,
        c.tolist(),
    )


def ols_fit(eps, mean, cov, degree):
    """Unweighted least squares; uncertainty propagated with the full covariance."""
    W = poly_weights(eps, degree)
    c = W @ mean
    var = W[0] @ cov @ W[0]
    return float(c[0]), float(np.sqrt(max(var, 0.0))), c.tolist()


def extrapolate(eps, mean, cov, primary_degree=2):
    """Polynomial eps->0 fits. Primary: OLS of degree primary_degree on all eps.

    The extrapolation uncertainty is the largest |shift| of J0 among the cubic fit and
    the two subset fits; each shift also carries its own correlated statistical sigma
    (shift_sigma), so a shift compatible with zero is visible as such."""
    eps = np.asarray(eps)
    mean = np.asarray(mean)
    cov = np.asarray(cov)
    n = len(eps)
    fits, wvec = {}, {}

    def add(name, idx, deg):
        W = poly_weights(eps[idx], deg)
        w = np.zeros(n)
        w[idx] = W[0]
        c = W @ mean[idx]
        fits[name] = {
            "value": float(c[0]),
            "stat": float(np.sqrt(max(w @ cov @ w, 0.0))),
            "coeffs": c.tolist(),
        }
        wvec[name] = w

    for deg in range(1, min(4, n)):
        add(f"ols_deg{deg}_all", np.arange(n), deg)
        if n - 1 > deg:
            add(f"ols_deg{deg}_drop_largest", np.arange(n - 1), deg)
        if n - 2 > deg:
            add(f"ols_deg{deg}_drop_two_smallest", np.arange(2, n), deg)
    deg = min(primary_degree, n - 1)
    pname = f"ols_deg{deg}_all"
    prim = fits[pname]
    shifts = []
    for key in (
        f"ols_deg{deg + 1}_all",
        f"ols_deg{deg}_drop_largest",
        f"ols_deg{deg}_drop_two_smallest",
    ):
        if key in fits:
            dw = wvec[key] - wvec[pname]
            fits[key]["shift"] = fits[key]["value"] - prim["value"]
            fits[key]["shift_sigma"] = float(np.sqrt(max(dw @ cov @ dw, 0.0)))
            shifts.append(abs(fits[key]["shift"]))
    syst = max(shifts, default=float("nan"))
    try:
        g = gls_fit(eps, mean, cov, deg)
        fits[f"gls_deg{deg}"] = {
            "value": g[0],
            "stat": g[1],
            "chi2": g[2],
            "dof": g[3],
            "coeffs": g[4],
        }
    except np.linalg.LinAlgError:
        pass
    return {
        "value": prim["value"],
        "stat": prim["stat"],
        "extrapolation": syst,
        "primary": pname,
        "fits": fits,
    }


def run(args):
    eps = sorted(args.eps)
    deg = min(args.degree, len(eps) - 1)
    w0 = poly_weights(eps, deg)[0]
    kwargs = {
        "graph": args.graph,
        "energy": args.energy,
        "momentum": args.momentum,
        "temperature": args.temperature,
        "eps_list": eps,
        "ng": args.ng,
        "nsub": args.nsub,
        "width_factor": args.width_factor,
        "extremum_window": args.extremum_window,
        "fit_weights": w0.tolist(),
    }
    f = Integrand(kwargs)
    p = Problem(**kwargs)
    dim = p.dim
    log = args.output.with_suffix(".log") if args.output else None

    def say(*a):
        print(*a, flush=True)
        if log is not None:
            with log.open("a") as fh:
                print(*a, file=fh)

    say(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} args={vars(args)}")
    say(f"# dim={dim} eps={eps} fit weights (deg {deg})={w0.tolist()}")
    np.random.seed(args.seed)
    gvar.ranseed(args.seed)
    integ = vegas.Integrator(
        [[0.0, 1.0]] * dim, nproc=args.nproc, min_neval_batch=args.batch
    )
    t0 = time.time()
    train = integ(f, nitn=args.nitn_train, neval=args.neval_train, alpha=0.5)
    say(f"# training ({time.time() - t0:.1f}s): {train[0]}  Q={train.Q:.2f}")
    say(train.summary())
    t1 = time.time()
    res = integ(f, nitn=args.nitn, neval=args.neval, adapt=False)
    t2 = time.time()
    say(f"# final ({t2 - t1:.1f}s): first component {res[0]}  Q={res.Q:.2f}")
    say(res.summary())
    vals = res[1:]
    mean = np.array([v.mean for v in vals])
    cov = gvar.evalcov(vals)
    extr = extrapolate(eps, mean, cov, primary_degree=deg)
    out = {
        "graph": args.graph,
        "energy": args.energy,
        "momentum": args.momentum,
        "temperature": args.temperature,
        "invariant_mass_squared": args.energy**2 - args.momentum**2,
        "Gamma": extr["value"],
        "stat_uncertainty": extr["stat"],
        "extrapolation_uncertainty": extr["extrapolation"],
        "first_component_direct": [res[0].mean, res[0].sdev],
        "eps": eps,
        "Gamma_eps_mean": mean.tolist(),
        "Gamma_eps_sdev": np.sqrt(np.diag(cov)).tolist(),
        "Gamma_eps_cov": cov.tolist(),
        "extrapolation": extr,
        "vegas": {
            "dim": dim,
            "nitn_train": args.nitn_train,
            "neval_train": args.neval_train,
            "nitn": args.nitn,
            "neval_per_itn": args.neval,
            "Q_final": res.Q,
            "chi2_dof_final": res.chi2 / max(res.dof, 1),
            "outer_samples_final": args.nitn * args.neval,
            "outer_samples_total": args.nitn * args.neval
            + args.nitn_train * args.neval_train,
            "seed": args.seed,
            "summary_final": res.summary(),
        },
        "quadrature": {
            "ng": args.ng,
            "nsub": args.nsub,
            "width_factor": args.width_factor,
            "extremum_window": args.extremum_window,
        },
        "runtime_s": {"train": t1 - t0, "final": t2 - t1, "total": t2 - t0},
        "nproc": args.nproc,
        "code_sha256": CODE_SHA256,
    }
    say(
        f"# RESULT {args.graph} E={args.energy} p={args.momentum} T={args.temperature}: "
        f"Gamma={extr['value']:.10e} stat={extr['stat']:.3e} extrap={extr['extrapolation']:.3e}"
    )
    for k, v in extr["fits"].items():
        say(
            f"#   {k}: {v['value']:.10e} +- {v['stat']:.3e}"
            + (f" chi2/dof={v['chi2']:.2f}/{v['dof']}" if "chi2" in v else "")
        )
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--graph", required=True, help="dt, tbt, or bubble:m1:m2 (one-loop control)"
    )
    ap.add_argument("--energy", type=float, required=True)
    ap.add_argument("--momentum", type=float, default=0.0)
    ap.add_argument("--temperature", type=float, default=0.5, help="0 = vacuum")
    ap.add_argument(
        "--eps", type=float, nargs="+", default=[0.002, 0.003, 0.004, 0.006, 0.008]
    )
    ap.add_argument(
        "--degree", type=int, default=2, help="polynomial degree of the eps->0 fit"
    )
    ap.add_argument("--neval", type=int, default=200000)
    ap.add_argument("--nitn", type=int, default=5)
    ap.add_argument("--neval-train", type=int, default=50000)
    ap.add_argument("--nitn-train", type=int, default=5)
    ap.add_argument("--nproc", type=int, default=48)
    ap.add_argument("--batch", type=int, default=20000)
    ap.add_argument("--ng", type=int, default=10)
    ap.add_argument("--nsub", type=int, default=6)
    ap.add_argument("--width-factor", type=float, default=0.5)
    ap.add_argument("--extremum-window", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args(argv)
    if args.nproc > 48:
        ap.error("at most 48 processes")
    run(args)


if __name__ == "__main__":
    main()
