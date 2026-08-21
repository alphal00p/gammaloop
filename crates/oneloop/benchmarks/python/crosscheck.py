#!/usr/bin/env python3
"""Cross-engine validation of the oneloop reducer (gitignored, local-only).

Reads the reductions emitted by `examples/emit_reductions.rs` and, for each family,
checks the oneloop IBP reduction against two independent engines:

  1. OneLOop (avh_olo, van Hameren) for the scalar master values A0/B0/C0/D0, and
  2. scipy Feynman-parameter Monte-Carlo for the *original* integral (direct).

For each family we form the full Laurent series

    V_reduced(eps) = sum_i c_i(d=4-2eps) * M_i(eps),   M_i from OneLOop,

and require (a) the 1/eps and 1/eps^2 poles cancel (the original is finite here) and
(b) the finite part equals the direct scipy value.  A running OneLOop-vs-feynalg tally
on the masters is a third, independent cross-check of the master values themselves.

    cargo run --release --example emit_reductions -p oneloop > /tmp/oneloop_reductions.txt
    python crosscheck.py /tmp/oneloop_reductions.txt
"""
import sys
from fractions import Fraction
from math import gamma, prod
import numpy as np
import sympy as sp
from scipy import integrate
from scipy.special import gammaln
import oneloop_bridge as ob

try:
    from feynalg.curated.analytic_integrals import (
        analytic_A0, analytic_B0, analytic_C0, analytic_D0,
    )
    HAVE_FEYNALG = True
except Exception:
    HAVE_FEYNALG = False

ob.set_renormalization_scale(1.0)      # mu^2 = 1, to match the bare scipy F
np.random.seed(20260718)
NS = 3_000_000
d = sp.Symbol("d")
eps = sp.Symbol("eps")

import re as _re
SGN_LQ = -1                            # dot(k,q)_Mink = -(k_E . q_E) for spacelike q (pinned)


def builder_for(num):
    """Tensor-oracle numerator builder for a mixed dot(k,q_j) label: an optional 'll'
    (= dot(k,k)) prefix followed by one or more q<i> factors -> the product of the
    matching LQ(q_i) moments (times LL() when 'll'-prefixed). Handles ARBITRARY rank
    (q1, q1q2, q1q2q3, q1q2q3q4, llq1, llq1q2, llq1q2q3, ...); else None. The moment
    tensor oracle (`tensor_direct_general`) integrates whatever product this builds."""
    m = _re.fullmatch(r"(ll)?((?:q\d)+)", num)
    if not m:
        return None
    has_ll = m.group(1) is not None
    idxs = [int(c) - 1 for c in _re.findall(r"q(\d)", m.group(2))]

    def build(qv, LL, LQ):
        expr = LL() if has_ll else 1
        for i in idxs:
            expr = expr * LQ(qv[i], SGN_LQ)
        return expr

    return build


# oneloop's MasterIntegral stores the pairwise Cayley invariants in LEXicographic order,
# not OneLOop's physical (adjacent-leg, diagonal) slot order.  So we identify each master's
# propagator lines by its MASSES and rebuild the args from the geometry (unambiguous), which
# tests oneloop's coefficients + topology selection without depending on the field ordering.
NMASS = {"A0": 1, "B0": 2, "C0": 3, "D0": 4}


def master_masses(kind, args):
    """The trailing entries of an emitted master's args are the propagator masses^2."""
    k = NMASS[kind]
    return [float(Fraction(x)) for x in args[-k:]]


def sij(offsets, i, j):
    d = np.array(offsets[i]) - np.array(offsets[j])
    return -float(d @ d)                      # spacelike invariant s_ij = -(r_i-r_j)^2


def lines_of(mmasses, cfg_masses):
    """Match each master mass^2 to a distinct config line (masses are distinct)."""
    idx, used = [], set()
    for mm in mmasses:
        best = min((k for k in range(len(cfg_masses)) if k not in used),
                   key=lambda k: abs(cfg_masses[k] - mm))
        idx.append(best)
        used.add(best)
    return sorted(idx)                        # sorted = cyclic order for a one-loop subtopo


def sinv_dict(pr):
    """Full pairwise {(i,j): s_ij} table from the emitted lex SINV line (i<j)."""
    n = len(pr["exps"].split(","))
    vals, d, t = pr["sinv"], {}, 0
    for i in range(n):
        for j in range(i + 1, n):
            d[(i, j)] = vals[t]
            t += 1
    return d


def sfun_for(pr):
    """(i,j) -> s_ij, from injected SINV (timelike) if present, else the offset geometry."""
    if "sinv" in pr:
        d = sinv_dict(pr)
        return lambda i, j: d[(min(i, j), max(i, j))]
    return lambda i, j: sij(pr["offsets"], i, j)


def geom_args(kind, lines, sfun, cfg_masses):
    m = [cfg_masses[i] for i in lines]
    if kind == "A0":
        return (m[0],)
    if kind == "B0":
        i, j = lines
        return (sfun(i, j), m[0], m[1])
    if kind == "C0":
        i, j, k = lines
        return (sfun(i, j), sfun(j, k), sfun(i, k), m[0], m[1], m[2])
    i, j, k, l = lines                        # box: legs (ij)(jk)(kl)(il), diagonals (ik)(jl)
    return (sfun(i, j), sfun(j, k), sfun(k, l), sfun(i, l),
            sfun(i, k), sfun(j, l), m[0], m[1], m[2], m[3])


def olo_master(kind, a):
    if kind == "A0":
        r = ob.one_point(a[0])
    elif kind == "B0":
        r = ob.two_point(*a)
    elif kind == "C0":
        r = ob.three_point(*a)
    else:
        r = ob.four_point(*a)
    return (complex(r.epsilon_minus_2), complex(r.epsilon_minus_1), complex(r.epsilon_0))


def feynalg_master(kind, a):
    if not HAVE_FEYNALG:
        return None
    try:
        v = {"A0": analytic_A0, "B0": analytic_B0,
             "C0": analytic_C0, "D0": analytic_D0}[kind](*a)
        return complex(v) if v is not None else None
    except Exception:
        return None


def scalar_laurent(lines, offsets, cfg_masses):
    """Laurent (em2, em1, e0) of the scalar integral (all exps=1) on `lines`.
    <=4 lines -> OneLOop master (exact, with poles); 5+ lines -> finite, scipy MC."""
    lines = sorted(lines)
    k = len(lines)
    if k <= 4:
        kind = {1: "A0", 2: "B0", 3: "C0", 4: "D0"}[k]
        return olo_master(kind, geom_args(kind, lines, lambda i, j: sij(offsets, i, j), cfg_masses))
    off_s = [offsets[i] for i in lines]
    m_s = [cfg_masses[i] for i in lines]
    v, _ = scipy_direct(off_s, m_s, [1] * k)         # finite N>=5 scalar
    return (0j, 0j, complex(v))


def true_numll(pr):
    """True Laurent of the k^2=dot(k,k) numerator integral: r0=0 so k^2 = D0 + m0^2,
    hence I[k^2] = I(drop line 0) + m0^2 * I(scalar).  Convention-free."""
    off, m = pr["offsets"], pr["masses"]
    n = len(m)
    m0 = m[0]
    drop = scalar_laurent(list(range(1, n)), off, m)
    full = scalar_laurent(list(range(n)), off, m)
    return tuple(drop[i] + m0 * full[i] for i in range(3))


# ---- direct scipy value of the original integral (Feynman MC, d=4) -----------
def scipy_direct(offsets, masses, exps):
    a = np.array(exps, dtype=float)
    A = int(sum(exps))
    n = len(exps)
    off = np.array(offsets)
    msq = np.array(masses)

    # Sampling density: flat Dir(a), UNLESS a line is (near-)massless. Then F -> 0 at the
    # x_l -> 1 corner (all other x -> 0) and F^(2-A) is near-singular there; flat Dir under-
    # samples it -> biased, heavy-tailed estimate. Concentrate on that corner via Dir(beta)
    # (small alpha on the massive lines) and reweight by Dir(a)/Dir(beta) to stay unbiased.
    massless = [i for i in range(n) if abs(msq[i]) < 1e-9]
    if massless:
        beta = a.copy()
        for i in range(n):
            if i not in massless:
                beta[i] = 0.3                    # pull massive-line x_i -> 0 (= corner)
        x = np.random.dirichlet(beta, size=NS)
        logB = lambda al: np.sum(gammaln(al)) - gammaln(np.sum(al))
        logw = (logB(beta) - logB(a)) + np.sum(
            (a - beta) * np.log(np.clip(x, 1e-300, None)), axis=1)   # Dir(a)/Dir(beta)
    else:
        x = np.random.dirichlet(a, size=NS)      # density prop to prod x_i^(a_i-1)
        logw = np.zeros(NS)

    F = x @ msq
    for i in range(n):
        for j in range(i + 1, n):
            dij = off[i] - off[j]
            F = F + x[:, i] * x[:, j] * float(dij @ dij)   # spacelike: + dist^2
    # I = (-1)^A Gamma(A - d/2)/Gamma(A) * < F^(d/2 - A) >_{Dir(a)}
    g = F ** (2.0 - A) * np.exp(logw)
    pref = (-1) ** A * gamma(A - 2.0) / gamma(A)
    val = pref * np.mean(g)
    err = abs(pref) * np.std(g) / np.sqrt(len(g))     # 1-sigma MC error of the mean
    return val, err


def scipy_quad(offsets, masses, exps, sfun=None):
    """Deterministic Feynman-parameter quadrature (n<=4), near-exact reference.
    I = (-1)^A Gamma(A-2)/prod Gamma(a_i) * INT_simplex prod x_i^(a_i-1) F^(2-A).
    F = sum x_i m_i^2 - sum_{i<j} x_i x_j s_ij; s_ij from sfun (Minkowski, may be timelike) or
    the spacelike offset geometry (s_ij = -dist^2, so -s_ij = +dist^2)."""
    a = exps
    A = sum(a)
    n = len(a)
    off = np.array(offsets)
    msq = np.array(masses)
    if sfun is None:
        sfun = lambda i, j: -float((off[i] - off[j]) @ (off[i] - off[j]))

    def integrand(xs):
        x = list(xs) + [1.0 - sum(xs)]
        F = sum(x[i] * msq[i] for i in range(n))
        for i in range(n):
            for j in range(i + 1, n):
                F += -x[i] * x[j] * sfun(i, j)
        w = prod(x[i] ** (a[i] - 1) for i in range(n))
        return w * F ** (2 - A)

    pref = (-1) ** A * gamma(A - 2.0) / prod(gamma(ai) for ai in a)
    if n == 2:
        v, _ = integrate.quad(lambda x0: integrand([x0]), 0, 1)
    elif n == 3:
        v, _ = integrate.dblquad(lambda x1, x0: integrand([x0, x1]),
                                 0, 1, 0, lambda x0: 1 - x0)
    elif n == 4:
        v, _ = integrate.tplquad(lambda x2, x1, x0: integrand([x0, x1, x2]),
                                 0, 1, 0, lambda x0: 1 - x0, 0, lambda x0, x1: 1 - x0 - x1)
    else:
        return None
    return pref * v


# ---- parse the emitted reduction file ---------------------------------------
def parse(path):
    procs, cur = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("PROCESS "):
            cur = {"terms": []}
            for kv in line[len("PROCESS "):].split():
                k, v = kv.split("=")
                cur[k] = v
        elif line.startswith("OFFSETS "):
            cur["offsets"] = [[float(c) for c in p.split(",")]
                              for p in line[len("OFFSETS "):].split(";")]
        elif line.startswith("MASSES "):
            cur["masses"] = [float(m) for m in line[len("MASSES "):].split(",")]
        elif line.startswith("SINV "):                 # explicit signed lex invariants (timelike)
            cur["sinv"] = [float(s) for s in line[len("SINV "):].split(",")]
        elif line.startswith("TERM "):
            body = line[len("TERM "):]
            assert body.startswith("coeff=( ")
            close = body.index(" ) ")            # outer delimiter (coeffs may contain "1/(d-3)")
            coeff = body[len("coeff=( "):close].strip()
            rest = body[close + 3:].split()
            cur["terms"].append((coeff, rest[0], rest[1:]))
        elif line.startswith("ENDPROCESS"):
            procs.append(cur)
            cur = None
    return procs


# eps-series of a coeff string in d, as (a0, a1, a2) with d = 4 - 2 eps
def coeff_series(coeff_str):
    expr = sp.sympify(coeff_str, locals={"d": d}).subs(d, 4 - 2 * eps)
    s = sp.series(expr, eps, 0, 3).removeO()
    p = sp.Poly(s, eps)
    return (complex(p.coeff_monomial(1)),
            complex(p.coeff_monomial(eps)),
            complex(p.coeff_monomial(eps ** 2)))


def two_oracle_agree(num, offsets, masses, exps, tv):
    """Independent second check: the covariant Minkowski-Gram oracle (a different computational
    path from the moment oracle) must agree with the primary value tv.  Returns (ok, tag)."""
    if not (num == "ll2" or _re.fullmatch(r"q\d|q\dq\d|llq\d", num)):
        return True, ""                                  # no covariant counterpart for this num
    from _tensor_oracle import tensor_covariant
    cv, _ = tensor_covariant(offsets, masses, exps, num)
    ok = abs(cv - tv) < 3e-3 * (1 + abs(tv))
    return ok, (" cov✓" if ok else " covX")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/oneloop_reductions.txt"
    procs = parse(path)
    print(f"{len(procs)} families from {path}\n")
    header = f"{'family':<16}{'cfg':>4}{'terms':>6}   {'|1/eps| resid':>13}   {'reduced (finite)':>20}   {'scipy direct':>20}   {'pull(sigma)':>11}   verdict"
    print(header)
    print("-" * len(header))

    fey_ok = fey_tot = 0
    n_pass = n_fail = 0
    for pr in procs:
        # V_reduced Laurent
        cm2 = cm1 = c0 = 0j
        for coeff, kind, args in pr["terms"]:
            a0, a1, a2 = coeff_series(coeff)
            if kind in ("A0", "B0"):
                # 1 invariant, symmetric in the masses -> emitted args are unambiguous
                a = [float(Fraction(x)) for x in args]
            else:
                # C0/D0: oneloop stores invariants lexicographically; rebuild the physical
                # (adjacent-leg, diagonal) OneLOop ordering from the geometry via the masses
                lines = lines_of(master_masses(kind, args), pr["masses"])
                a = geom_args(kind, lines, sfun_for(pr), pr["masses"])
            mm2, mm1, m0 = olo_master(kind, a)
            cm2 += a0 * mm2
            cm1 += a0 * mm1 + a1 * mm2
            c0 += a0 * m0 + a1 * mm1 + a2 * mm2
            # third check: OneLOop vs feynalg on this master's finite part
            fv = feynalg_master(kind, a)
            if fv is not None:
                fey_tot += 1
                if abs(fv - m0) <= 1e-6 * (1 + abs(m0)):
                    fey_ok += 1

        exps = [int(e) for e in pr["exps"].split(",")]
        num = pr.get("num", "1")
        reduced = c0.real
        n = len(exps)
        if num == "ll":
            # k^2 = D0 + m0^2 decomposition (may be divergent -> compare the full Laurent).
            tm2, tm1, t0 = true_numll(pr)
            direct = t0.real
            pole = abs(cm1 - tm1) + abs(cm2 - tm2)
            resid = pole + abs(c0 - t0)
            scale = 1 + abs(t0) + abs(tm1) + abs(tm2)
            tol = 5e-3 if n >= 5 else 1e-5          # n>=5 true value is scipy-MC limited
            metric, ok = f"{resid / scale:>9.2e} L", resid / scale < tol
        elif num in ("ll2", "ll3", "lld"):
            # (k^2)^p finite irreducible-tensor integral via the moment tensor oracle.
            # "lld" = k^2 on a finite dotted family (p=1 with raised propagator powers).
            from _tensor_oracle import tensor_direct
            p_tensor = 1 if num == "lld" else int(num[2])
            tv, terr = tensor_direct(pr["offsets"], pr["masses"], exps, p_tensor)
            direct = tv
            pole = abs(cm1) + abs(cm2)
            pull = abs(c0.real - tv) / (terr + 1e-300)
            agree, tag = two_oracle_agree(num, pr["offsets"], pr["masses"], exps, tv)
            metric = f"{pull:>7.1f}σ{tag}"
            ok = pole < 1e-6 * (1 + abs(c0.real)) and pull < 5.0 and agree
        elif num in ("ll2d", "ll3d"):
            # UV-divergent (k^2)^p: compare the full Laurent (1/eps pole + finite).
            from _tensor_oracle import tensor_laurent_llp
            tpole, tfin = tensor_laurent_llp(pr["offsets"], pr["masses"], exps, int(num[2]))
            direct = tfin
            pole = abs(cm1.real - tpole)                # 1/eps coefficient residual
            scale = 1 + abs(tfin) + abs(tpole)
            ok_p = pole < 5e-3 * scale
            ok_f = abs(c0.real - tfin) < 5e-3 * scale
            metric = f"P{cm1.real:+.3e}"
            ok = ok_p and ok_f
            print(f"{pr['name']:<16}{pr['cfg']:>4}  pole: oneloop {cm1.real:+.5e} vs "
                  f"oracle {tpole:+.5e} [{'OK' if ok_p else 'X'}]   finite: oneloop "
                  f"{c0.real:+.5e} vs oracle {tfin:+.5e} [{'OK' if ok_f else 'X'}]")
        elif builder_for(num):
            # dot(k,q_j) numerators (finite) via the general moment tensor oracle.
            from _tensor_oracle import tensor_direct_general
            tv, terr = tensor_direct_general(pr["offsets"], pr["masses"], exps, builder_for(num))
            direct = tv
            pole = abs(cm1) + abs(cm2)
            pull = abs(c0.real - tv) / (terr + 1e-300)
            agree, tag = two_oracle_agree(num, pr["offsets"], pr["masses"], exps, tv)
            metric = f"{pull:>7.1f}σ{tag}"
            ok = pole < 1e-6 * (1 + abs(c0.real)) and pull < 5.0 and agree
        elif pr["name"].startswith("tlx_"):
            # ABOVE threshold (complex, absorptive): scipy can't cross F=0.  Independent check via
            # the IBP mass-derivative I[a_d=2] = +/- d/dm_d^2 I[scalar], RHS finite-differenced from
            # OneLOop's COMPLEX master -> validates the reduction coefficients in the complex plane.
            d_line = next(i for i, e in enumerate(exps) if e > 1)
            kind = {3: "C0", 4: "D0"}[n]
            sf, ms, h = sfun_for(pr), list(pr["masses"]), 1e-4

            def scal0(shift):
                m2 = ms.copy(); m2[d_line] += shift
                return olo_master(kind, geom_args(kind, list(range(n)), sf, m2))[2]

            deriv = (scal0(h) - scal0(-h)) / (2 * h)
            resid = min(abs(c0 - deriv), abs(c0 + deriv))       # accept either IBP sign convention
            pole = abs(cm1) + abs(cm2)
            ok = pole < 1e-6 * (1 + abs(c0)) and resid < 1e-3 * (1 + abs(deriv))
            print(f"{pr['name']:<16}{pr['cfg']:>4}  reduced {c0.real:+.5e}{c0.imag:+.5e}j  vs "
                  f"d/dm^2 {deriv.real:+.5e}{deriv.imag:+.5e}j  |Im|={abs(c0.imag):.2e}  "
                  f"resid {resid:.1e} [{'OK' if ok else 'X'}]")
            reduced, direct, metric = c0.real, deriv.real, f"{resid:>9.2e} d"
        elif n <= 4:                                 # scalar/dotted, deterministic quadrature
            direct = float(scipy_quad(pr.get("offsets"), pr["masses"], exps, sfun_for(pr)))
            pole = abs(cm1) + abs(cm2)
            rel = abs(reduced - direct) / (abs(direct) + 1e-300)
            metric, ok = f"{rel:>9.2e} q", pole < 1e-6 * (1 + abs(reduced)) and rel < 2e-3
        else:                                        # scalar/dotted N>=5, MC (pull in sigma)
            dv, mcerr = scipy_direct(pr["offsets"], pr["masses"], exps)
            direct = dv.real
            pole = abs(cm1) + abs(cm2)
            pull = abs(reduced - direct) / (mcerr + 1e-300)
            metric, ok = f"{pull:>8.1f}σ ", pole < 1e-6 * (1 + abs(reduced)) and pull < 5.0
        n_pass += ok
        n_fail += (not ok)
        print(f"{pr['name']:<16}{pr['cfg']:>4}{len(pr['terms']):>6}   "
              f"{pole:>13.2e}   {reduced:>+20.10e}   {direct:>+20.10e}   "
              f"{metric:>11}   {'PASS' if ok else 'FAIL'}")

    print("-" * len(header))
    print(f"\n{n_pass}/{n_pass + n_fail} families PASS "
          f"(poles cancel to <1e-6; reduced finite part within 5σ of scipy MC)")
    if fey_tot:
        print(f"master finite-part agreement OneLOop vs feynalg: {fey_ok}/{fey_tot}")


if __name__ == "__main__":
    main()
