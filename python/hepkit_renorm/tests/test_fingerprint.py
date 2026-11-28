"""Exact plane-wave fingerprints (hepkit_renorm.fingerprint).

Golden data: ``tests/data/st_engine_golden.json`` (values of the numpy engines ``SCRATCHPAD/st`` in four dimensions
and ``SCRATCHPAD/phase4/pole_consistency`` in D = 6 at unit couplings, see tests/data/README.md). The functionals
of the numpy engines are rebuilt here from the core modules: the classical action of :mod:`hepkit_renorm.action`
plus the ghost and source terms of ``st/model.py::dens_S0`` and some basis monomials of ``st/bases.py``, written
term by term in the same order of the Grassmann factors.
"""

from __future__ import annotations

import itertools
import json
from fractions import Fraction

import numpy as np
import pytest
from conftest import DATA
from symbolica import Expression, S

from hepkit_renorm import linalg as la
from hepkit_renorm import tensors as tn
from hepkit_renorm.action import classical_action, classical_sectors, covariant_derivative
from hepkit_renorm.exact import second_field
from hepkit_renorm.fingerprint import Fingerprinter, Representation
from hepkit_renorm.functional import Density, LocalFunctional, Occurrence
from hepkit_renorm.groups import generator_expr

UNIT = {"g_s": 1, "g": 1, "gp": 1, "y_u": 1, "y_d": 1, "y_e": 1, "lambda": 1, "m2": 1}
HALF = Expression.num(1) / 2


# ---------------------------------------------------------------------------------------------- builders
def _mono(model, coeff, *occs):
    return Density.monomial(model.fields, coeff, *occs)


def _gen_terms(model, field, labels, slot_of_rep=None):
    """``[(coefficient, inner labels)]`` of ``g_r (T^a_R)_{label, inner}`` over the group slots and U(1) charges of
    ``field`` (``a`` the returned adjoint label, ``None`` for U(1)); ``slot_of_rep``: representation source."""
    info = model.fields[field]
    src = model.fields[slot_of_rep or field]
    out = []
    for k, slot in enumerate(info.slots):
        if slot.factor is None:
            continue
        fac = model.factors[slot.factor]
        a, j = tn.fresh(fac.adjoint_space()), tn.fresh(slot.space)
        inner = list(labels)
        inner[k] = j
        out.append((slot.factor, a, model.coupling(slot.factor) * generator_expr(fac, src.slots[k].rep, a, labels[k], j),
                    tuple(inner)))
    for fac_name, q in src.charges.items():
        out.append((fac_name, None, model.coupling(fac_name) * Expression.num(q.numerator) / q.denominator, tuple(labels)))
    return out


def _gh(factor, a):
    return Occurrence(f"gh_{factor}", () if a is None else (a,))


def source_terms(model) -> LocalFunctional:
    """The ghost and source terms of ``st/model.py::dens_S0`` (tilde S_0), with st's order of the factors."""
    f = model.fields
    total = Density(f)
    for name, fac in model.factors.items():
        mu = tn.fresh(tn.LORENTZ)
        a = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        total = total + _mono(model, 1, Occurrence(f"src_{name}", (mu,) + a)) * covariant_derivative(
            model, f"gh_{name}", mu, a)
        if not fac.abelian:
            a, b, c = (tn.fresh(fac.adjoint_space()) for _ in range(3))
            total = total + _mono(model, -HALF * model.coupling(name) * tn.f(a, b, c), Occurrence(f"src_gh_{name}", (a,)),
                                  _gh(name, b), _gh(name, c))
    for name, fm in model.fermions.items():
        pcb = tn.PR if fm.chirality == "L" else tn.PL
        s, t = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)
        labels = (s,) + model.slot_labels(name)[1:]
        for fac, a, gen, inner in _gen_terms(model, name, labels):
            # i Rbar wm P_c Psi  (order Rbar, omega, Psi)
            total = total + rbar_term(model, name, (fac, a, gen, inner), labels)
            # i omega Psibar P_cbar T R  (order omega, Psibar, R)
            total = total + _mono(model, tn.I * gen * pcb(s, t), _gh(fac, a), Occurrence(f"{name}bar", labels),
                                  Occurrence(f"src_{name}bar", (t,) + inner[1:]))
    for name in model.scalars:
        labels = model.slot_labels(name)
        for fac, a, gen, inner in _gen_terms(model, name, labels):
            total = total + _mono(model, tn.I * gen, Occurrence(f"src_{name}", labels), _gh(fac, a),
                                  Occurrence(name, inner))
            total = total + _mono(model, -tn.I * gen, _gh(fac, a), Occurrence(f"{name}dag", labels),
                                  Occurrence(f"src_{name}dag", inner))
    return LocalFunctional.integrate(total)


def rbar_term(model, name, gen_term, labels, swap=False) -> Density:
    """``i Rbar_f (g T^a omega^a) P_c Psi_f`` for one generator term, occurrences in st's order (``Rbar, omega,
    Psi``) or, with ``swap``, deliberately misordered as ``Rbar, Psi, omega`` without the Grassmann sign."""
    fac, a, gen, inner = gen_term
    pc = tn.PL if model.fermions[name].chirality == "L" else tn.PR
    t = tn.fresh(tn.SPINOR)
    occs = [Occurrence(f"src_{name}", labels), _gh(fac, a), Occurrence(name, (t,) + inner[1:])]
    if swap:
        occs = [occs[0], occs[2], occs[1]]
    return _mono(model, tn.I * gen * pc(labels[0], t), *occs)


def delta1(model) -> LocalFunctional:
    """eq:Delta1: ``omega^a [Psibar T^a ghat^mu P_cbar d_mu Psi + (d_mu Psibar) ghat^mu P_c T^a Psi]``."""
    total = Density(model.fields)
    for name, fm in model.fermions.items():
        pc, pcb = (tn.PL, tn.PR) if fm.chirality == "L" else (tn.PR, tn.PL)
        s, t, u = (tn.fresh(tn.SPINOR) for _ in range(3))
        mu = tn.fresh(tn.LORENTZ)
        labels = (s,) + model.slot_labels(name)[1:]
        for fac, a, gen, inner in _gen_terms(model, name, labels):
            total = total + _mono(model, gen * tn.gammahat(mu, s, t) * pcb(t, u), _gh(fac, a),
                                  Occurrence(f"{name}bar", labels), Occurrence(name, (u,) + inner[1:], (mu,)))
            total = total + _mono(model, gen * tn.gammahat(mu, s, t) * pc(t, u), _gh(fac, a),
                                  Occurrence(f"{name}bar", labels, (mu,)), Occurrence(name, (u,) + inner[1:]))
    return LocalFunctional.integrate(total)


def basis_monomials(model) -> dict[str, LocalFunctional]:
    """Some monomials of ``st/bases.py`` (unit couplings there; the coupling symbols are kept here)."""
    g, gp, gs = model.coupling("W"), model.coupling("B"), model.coupling("G")
    mu, nu, rho, sig = (tn.fresh(tn.LORENTZ) for _ in range(4))
    a, b, c = (tn.fresh(tn.adjoint_space("W", 3)) for _ in range(3))
    A = tn.fresh(tn.adjoint_space("G", 8))
    out = {}
    divg = _mono(model, 1, Occurrence("G", (mu, A), (mu,)))
    out["L_G"] = LocalFunctional.integrate((divg * divg).scale(-HALF))
    out["WBeps"] = LocalFunctional.integrate(_mono(model, g**2 * gp * tn.eps(mu, nu, rho, sig), Occurrence("W", (mu, a)),
                                                   Occurrence("W", (rho, a), (nu,)), Occurrence("B", (sig,))))
    out["WWWBeps"] = LocalFunctional.integrate(_mono(model, g**3 * gp * tn.eps(mu, nu, rho, sig) * tn.f(a, b, c),
                                                     Occurrence("W", (mu, a)), Occurrence("W", (nu, b)),
                                                     Occurrence("W", (rho, c)), Occurrence("B", (sig,))))
    i, j = tn.fresh(tn.fundamental_space("W", 2)), tn.fresh(tn.fundamental_space("W", 2))
    t_ = g * tn.T(a, i, j)
    out["J_W"] = LocalFunctional.integrate(
        _mono(model, tn.I * t_, Occurrence("Hdag", (i,), (mu,)), Occurrence("W", (mu, a)), Occurrence("H", (j,)))
        - _mono(model, tn.I * t_, Occurrence("Hdag", (i,)), Occurrence("W", (mu, a)), Occurrence("H", (j,), (mu,))))
    s, t, u = (tn.fresh(tn.SPINOR) for _ in range(3))
    ci, cj = tn.fresh(tn.fundamental_space("G", 3)), tn.fresh(tn.fundamental_space("G", 3))
    w = tn.fresh(tn.fundamental_space("W", 2))
    first = _mono(model, tn.I * gs * tn.T(A, ci, cj) * tn.PL(s, t), Occurrence("src_Q", (s, ci, w)),
                  Occurrence("gh_G", (A,)), Occurrence("Q", (t, cj, w)))
    second = _mono(model, tn.I * gs * tn.PR(u, s) * tn.T(A, ci, cj), Occurrence("gh_G", (A,)),
                   Occurrence("Qbar", (u, ci, w)), Occurrence("src_Qbar", (s, cj, w)))
    out["R_QG+"] = LocalFunctional.integrate(first + second)
    out["R_QG-"] = LocalFunctional.integrate(first - second)
    return out


# ---------------------------------------------------------------------------------------------- st bridge
_GAUGE = {"A": "", "rt": "src_", "w": "gh_", "z": "src_gh_"}
_RANGE = {"G": slice(0, 8), "W": slice(8, 11), "B": slice(11, 12)}


def _st_options(model, st_name):
    """``[(field, converter)]``: the fields of the package that the st field ``st_name`` stands for."""
    if st_name in _GAUGE:
        pre = _GAUGE[st_name]
        out = []
        for r, sl in _RANGE.items():
            lorentz = st_name in ("A", "rt")
            def conv(e, sl=sl, lorentz=lorentz, abelian=model.factors[r].abelian):
                x = e[sl]                                  # [a, (mu)]
                if lorentz:
                    x = (x * np.array([1, -1, -1, -1] + [-1] * (x.shape[-1] - 4))).T   # [mu, a] upper
                    return x[:, 0] if abelian else x
                return x[0] if abelian else x
            out.append((f"{pre}{r}", conv))
        return out
    simple = {"H": "H", "Hd": "Hdag", "sd": "src_H", "s": "src_Hdag"}
    if st_name in simple:
        return [(simple[st_name], lambda e: e)]
    kind, mult = st_name.split("_")
    field = {"Psi": mult, "Psib": f"{mult}bar", "Rb": f"src_{mult}", "R": f"src_{mult}bar"}[kind]
    dims = [d for d in (3 if "G" in model.fermions[mult].reps else 1, 2 if "W" in model.fermions[mult].reps else 1)]

    def conv(e):                                           # [component (colour*2 + isospin), spinor] -> [spinor, ...]
        x = e.reshape(dims + [e.shape[-1]])
        x = np.moveaxis(x, -1, 0)
        return x.reshape((e.shape[-1],) + tuple(d for d in dims if d > 1))
    return [(field, conv)]


def st_value(model, fp: Fingerprinter, functional: LocalFunctional, case: dict) -> int:
    """The package's value for a golden case: sum over the package fields each st leg stands for."""
    rep = fp.rep
    D = fp.dimension
    metric = np.array([1] + [-1] * (D - 1))
    moms = np.array(case["momenta_lower"], dtype=object) * metric          # upper components
    pols = [np.vectorize(lambda r, i: rep.complex(r, i), otypes=[object])(np.array(pr), np.array(pi))
            for pr, pi in case["polarisations"]]
    total = 0
    for choice in itertools.product(*(_st_options(model, n) for n in case["legs"])):
        cls = [c[0] for c in choice]
        conv = [c[1](e) for c, e in zip(choice, pols)]
        smp = fp.explicit_samples(cls, (moms % rep.p).astype(np.int64)[None], [np.asarray(x, dtype=object)[None] for x in conv])
        total = (total + int(fp.functional(functional, smp)[0])) % rep.p
    return total


@pytest.fixture(scope="module")
def golden():
    with open(DATA / "st_engine_golden.json") as fh:
        return json.load(fh)["cases"]


@pytest.fixture(scope="module")
def st_functionals(sm):
    s0 = LocalFunctional.integrate(classical_action(sm) + source_terms(sm))
    funcs = {"S0": s0, "Delta1": delta1(sm)}
    funcs.update(basis_monomials(sm))
    return funcs


def _expected(rep: Representation, case: dict) -> int:
    re, im = (Fraction(x) for x in case["value"])
    assert abs(float(re) - case["float"][0]) + abs(float(im) - case["float"][1]) < 1e-8 * max(1, abs(case["float"][0]))
    return rep.complex(re, im)


# ---------------------------------------------------------------------------------------------- tests
@pytest.mark.parametrize("D", [4, 6, 8])
def test_explicit_clifford_algebra(field, D):
    rep = Representation(field, D)
    p, S = rep.p, rep.spinor_dim
    g = rep.gamma.astype(object)
    one = np.eye(S, dtype=object)
    g5 = rep.gamma5.astype(object)
    for m, n in itertools.product(range(D), repeat=2):
        assert np.array_equal((g[m].dot(g[n]) + g[n].dot(g[m])) % p, one * (2 * rep.metric[m] * (m == n)) % p)
    for m in range(D):
        comm = (g[m].dot(g5) + (1 if m < 4 else -1) * g5.dot(g[m])) % p
        assert not comm.any()                                # {gammabar, g5} = 0, [gammahat, g5] = 0
    assert np.array_equal(g5.dot(g5) % p, one)
    pl, pr = rep.projectors
    assert np.array_equal((pl + pr) % p, rep.identity) and not (pl.astype(object).dot(pr.astype(object)) % p).any()
    assert rep.eps[0, 1, 2, 3] == 1 and rep.eps[1, 0, 2, 3] == p - 1
    assert int(np.trace(one)) * rep.trace_norm % p == 4


def test_samples_are_deterministic(sm, field):
    fp = Fingerprinter(sm, field)
    a, b = fp.samples(("Qbar", "Q", "G"), 5), fp.samples(("Qbar", "Q", "G"), 3)
    assert np.array_equal(a.polarisations[2][:3], b.polarisations[2]) and a.polarisations[0].shape == (5, 4, 3, 2)
    assert not np.any(sum(a.momenta.values()) % field.p)
    assert not np.array_equal(a.polarisations[1], fp.samples(("Qbar", "Q", "G"), 5, salt=1).polarisations[1])
    assert np.array_equal(fp.replace(i_sign=-1).samples(("Qbar", "Q", "G"), 3).momenta[1], b.momenta[1])
    with pytest.raises(ValueError):
        fp.explicit_samples(("H", "Hdag"), np.ones((1, 2, 4), dtype=np.int64), [np.ones((1, 2)), np.ones((1, 2))])


def test_agreement_with_numpy_engines(sm, field, golden, st_functionals):
    """V1-type check against the stand-alone engines (classical action, ghost and source terms, basis monomials,
    Delta_1 and the d-dimensional action in D = 6), including leg orders that differ by Grassmann signs."""
    for sign in (1, -1):         # i_sign = -1: the Galois-conjugate computation (inputs and value conjugated)
        fps = {D: Fingerprinter(sm, field, D, values=UNIT, i_sign=sign) for D in (4, 6)}
        for case in golden:
            fp = fps[case["dimension"]]
            got = st_value(sm, fp, st_functionals[case["functional"]], case)
            assert got == _expected(fp.rep, case), (sign, case["functional"], case["legs"])


def test_brst_module_against_fingerprints(sm, field):
    """Cross-check of the parallel ``brst`` module (skipped while its interface is not available): L_ext equals the
    source terms of the numpy engine, Delta_1 = s S_0 equals eq:Delta1 in D = 6, and S(S_0) = Delta_1."""
    brst = pytest.importorskip("hepkit_renorm.brst")
    try:
        br = brst.BRST(sm)
        ext, breaking, action = br.external, br.breaking, br.action
        st = br.slavnov_taylor(action)
    except (AttributeError, TypeError, NotImplementedError) as err:
        pytest.skip(f"brst interface not available: {err}")
    fp4, fp6 = Fingerprinter(sm, field, 4), Fingerprinter(sm, field, 6)
    assert fp4.equal(ext, source_terms(sm))
    assert fp6.equal(breaking, delta1(sm)) and fp4.vanishes(breaking)
    assert fp6.equal(st, breaking) and fp4.vanishes(st)


def test_wrong_signs_are_detected(sm, field, golden, st_functionals):
    fp = Fingerprinter(sm, field, values=UNIT)
    p = field.p
    cases = {(c["functional"], tuple(c["legs"])): c for c in golden if c["dimension"] == 4}
    # 1. Grassmann order: Rbar Psi omega instead of Rbar omega Psi flips the sign of every term
    labels = (tn.fresh(tn.SPINOR),) + sm.slot_labels("Q")[1:]
    gens = _gen_terms(sm, "Q", labels)
    right, wrong = (LocalFunctional.integrate(sum((rbar_term(sm, "Q", g, labels, swap) for g in gens[1:]),
                                                  rbar_term(sm, "Q", gens[0], labels, swap))) for swap in (False, True))
    case = cases[("S0", ("Rb_Q", "w", "Psi_Q"))]
    good = _expected(fp.rep, case)
    assert st_value(sm, fp, right, case) == good
    assert st_value(sm, fp, wrong, case) == (-good) % p
    # 2. relative sign of the hermitian-conjugate Yukawa term
    sec = classical_sectors(sm)
    yuk = sec["yukawa"]
    flipped = LocalFunctional.integrate(Density.from_monomials(
        sm.fields, [type(m)(-m.coefficient if "Ubar" in m.fields() else m.coefficient, m.occurrences) for m in yuk]))
    case = cases[("S0", ("Psib_U", "Psi_Q", "H"))]
    assert st_value(sm, fp, yuk, case) == _expected(fp.rep, case)
    assert st_value(sm, fp, flipped, case) != _expected(fp.rep, case)
    # 3. orientation of the Levi-Civita tensor
    basis = basis_monomials(sm)
    case = cases[("WBeps", ("A", "A", "A"))]
    assert st_value(sm, fp, basis["WBeps"], case) == _expected(fp.rep, case)
    w = [S(f"{x}_") for x in "abcd"]
    mirrored = LocalFunctional.integrate(Density.from_monomials(
        sm.fields, [type(m)(m.coefficient.replace(tn.eps(*w), tn.eps(w[1], w[0], w[2], w[3])), m.occurrences)
                    for m in basis["WBeps"]]))
    assert st_value(sm, fp, mirrored, case) == (-_expected(fp.rep, case)) % p
    # 4. a wrong sign of the U(1) coupling is visible in the U-U-(G, B) vertex
    wrong = Fingerprinter(sm, field, values={**UNIT, "gp": -1})
    case = cases[("S0", ("Psib_U", "Psi_U", "A"))]
    assert st_value(sm, wrong, sec["fermion"], case) != _expected(fp.rep, case)


@pytest.mark.parametrize("D", [4, 6])
def test_simplifier_is_sound(sm, field, s0_classical, D):
    """The symbolic BMHV simplification does not change any vertex function of S_0 (both modes)."""
    fp = Fingerprinter(sm, field, D)
    for mode in (True, False):
        simp = sm.simplifier(expand_projectors=mode)
        for cls in s0_classical.leg_classes():
            smp = fp.samples(cls, 2)
            raw = fp.functional(s0_classical, smp)
            v = s0_classical.vertex_function(list(smp.legs), simp)   # with momentum conservation
            assert np.array_equal(raw, fp.evaluate(v, smp)), (mode, cls)
            assert raw.any()


def test_integration_by_parts_and_total_derivatives(sm, field, s0_classical):
    fp = Fingerprinter(sm, field)
    assert fp.equal(s0_classical.integrate_by_parts("H"), s0_classical)
    mu = tn.fresh(tn.LORENTZ)
    labels = sm.slot_labels("H")
    dens = _mono(sm, 1, Occurrence("Hdag", labels), Occurrence("H", labels), Occurrence("B", (mu,)))
    assert fp.vanishes(LocalFunctional.integrate(dens.total_derivative(mu)))
    not_total = _mono(sm, 1, Occurrence("Hdag", labels), Occurrence("H", labels), Occurrence("B", (mu,), (mu,)))
    assert not fp.vanishes(LocalFunctional.integrate(not_total))


def test_evanescent_structures_need_d_above_four(sm, field):
    simp = sm.simplifier()
    fermion = classical_sectors(sm)["fermion"]
    ev = LocalFunctional.integrate(fermion.evanescent(simp))
    ph = LocalFunctional.integrate(fermion.physical(simp))
    fp4, fp6, fp8 = (Fingerprinter(sm, field, D) for D in (4, 6, 8))
    assert fp4.vanishes(ev) and not fp6.vanishes(ev)
    assert fp6.equal(LocalFunctional.integrate(ev + ph), fermion)
    assert fp4.vanishes(delta1(sm)) and not fp6.vanishes(delta1(sm))
    # (1/2) Abar dhat^2 Abar: V(G, G) = -gbar delta sphat(1, 1) after momentum conservation (sp path)
    mu, nu, rho, sig = (tn.fresh(tn.LORENTZ) for _ in range(4))
    a = tn.fresh(tn.adjoint_space("G", 8))
    ea = LocalFunctional.integrate(_mono(sm, HALF * tn.gbar(mu, nu) * tn.ghat(rho, sig), Occurrence("G", (mu, a)),
                                         Occurrence("G", (nu, a), (rho, sig))))
    smp = fp6.samples(("G", "G"), 3)
    (l1, a1), (l2, a2) = smp.legs[0].labels, smp.legs[1].labels
    expected = -tn.gbar(l1, l2) * tn.delta(a1, a2) * tn.sphat(1, 1)
    assert np.array_equal(fp6.functional(ea, smp), fp6.evaluate(expected, smp)) and fp6.functional(ea, smp).any()
    assert fp4.vanishes(ea)
    # three antisymmetrised hatted gammas need at least three hatted dimensions: zero in D = 6, not in D = 8
    s, t, u, v = (tn.fresh(tn.SPINOR) for _ in range(4))
    ls = (mu, nu, rho)
    anti = Expression.num(0)
    for perm in itertools.permutations(range(3)):
        sign = 1 if sum(1 for x, y in itertools.combinations(perm, 2) if x > y) % 2 == 0 else -1
        x, y, z = (ls[k] for k in perm)
        anti = anti + sign * tn.gammahat(x, s, t) * tn.gammahat(y, t, u) * tn.gammahat(z, u, v)
    ci = tn.fresh(tn.fundamental_space("G", 3))
    wi = tn.fresh(tn.fundamental_space("W", 2))
    three = LocalFunctional.integrate(_mono(sm, anti, Occurrence("Qbar", (s, ci, wi)), Occurrence("Q", (v, ci, wi)),
                                            Occurrence("B", (mu,)), Occurrence("B", (nu,), (rho,))))
    assert fp6.vanishes(three) and not fp8.vanishes(three)


def test_reconstruction_of_rational_coefficients(sm):
    """Known rational (and Gaussian-rational) coefficients of a combination are recovered exactly."""
    mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
    a = tn.fresh(tn.adjoint_space("G", 8))
    f_ = _mono(sm, 1, Occurrence("G", (nu, a), (mu,))) - _mono(sm, 1, Occurrence("G", (mu, a), (nu,)))
    T = LocalFunctional.integrate((f_ * f_).scale(-Expression.num(1) / 4))
    div = _mono(sm, 1, Occurrence("G", (mu, a), (mu,)))
    L = LocalFunctional.integrate((div * div).scale(-HALF))
    M = LocalFunctional.integrate(_mono(sm, HALF, Occurrence("G", (mu, a)), Occurrence("G", (mu, a))))
    sec = classical_sectors(sm)
    kin_q = LocalFunctional.integrate(Density.from_monomials(sm.fields, [m for m in sec["fermion"]
                                                                         if m.fields() == ("Q", "Qbar")]))
    basis = [T, L, M, kin_q, sec["potential"]]
    coeffs = [Fraction(3, 7), Fraction(-11, 6), Fraction(5), Fraction(2, 9), Fraction(-13, 4)]
    target = LocalFunctional(sm.fields)
    for c, b in zip(coeffs, basis):
        target = LocalFunctional.integrate(target + b.scale(Expression.num(c.numerator) / c.denominator))

    def solve_in(field):
        fp = Fingerprinter(sm, field)
        m = fp.matrix(basis)
        rhs = fp.matrix([target], layout=m.layout).values[:, 0]
        assert la.rank(m.values, field) == len(basis)
        return la.solve(m.values, rhs, field, require_unique=True)

    assert list(la.reconstruct_rationals(solve_in)) == coeffs
    # a dependent set: T + 2 L is in the span; its kernel vector is reconstructed
    fp = Fingerprinter(sm)
    m = fp.matrix([T, L, M, LocalFunctional.integrate(T + L.scale(2))])
    assert la.rank(m.values, fp.field) == 3
    ker = la.nullspace(m.values, fp.field)
    assert [fp.field.reconstruct(int(x)) for x in ker[0]] == [-1, -2, 0, 1]
    # Gaussian-rational coefficient: solve in both embeddings of i and split
    c = tn.I * Expression.num(3) / 5 + HALF
    z = LocalFunctional.integrate(kin_q.scale(c))
    sol = []
    for sign in (1, -1):
        fps = fp.replace(i_sign=sign)
        mm = fps.matrix([kin_q])
        sol.append(la.solve(mm.values, fps.matrix([z], layout=mm.layout).values[:, 0], fps.field)[0])
    re, im = la.split_gaussian(sol[0], sol[1], fp.field)
    assert (fp.field.reconstruct(int(re)), fp.field.reconstruct(int(im))) == (Fraction(1, 2), Fraction(3, 5))


def test_reconstruction_of_coupling_dependence(sm):
    """Coefficients of the scalar kinetic term on unit-coupling monomials, as polynomials in (g, g', g_s)."""
    scalar = classical_sectors(sm)["scalar"]
    mu = tn.fresh(tn.LORENTZ)
    a = tn.fresh(tn.adjoint_space("W", 3))
    i, j = tn.fresh(tn.fundamental_space("W", 2)), tn.fresh(tn.fundamental_space("W", 2))
    hwb = LocalFunctional.integrate(_mono(sm, tn.T(a, i, j), Occurrence("Hdag", (i,)), Occurrence("H", (j,)),
                                          Occurrence("W", (mu, a)), Occurrence("B", (mu,))))
    hbb = LocalFunctional.integrate(_mono(sm, 1, Occurrence("Hdag", (i,)), Occurrence("H", (i,)),
                                          Occurrence("B", (mu,)), Occurrence("B", (mu,))))

    def evaluate(field, point):
        fp = Fingerprinter(sm, field, values=dict(zip(("g", "gp", "g_s"), point)))
        out = []
        for basis in (hwb, hbb):
            cls = basis.leg_classes()
            m = fp.matrix([basis], cls)
            rhs = fp.matrix([scalar], layout=m.layout).values[:, 0]
            out.append(int(la.solve(m.values, rhs, field)[0]))
        return out

    polys = la.reconstruct_polynomials(evaluate, 3, 2)
    assert polys == [{(1, 1, 0): Fraction(1)}, {(0, 2, 0): Fraction(1, 4)}]


def test_second_prime_and_matrix_layout(sm, s0_classical):
    """Ranks agree in two independent primes; ``layout`` reproduces the rows."""
    sec = list(classical_sectors(sm).values())
    ranks = []
    for field in (None, second_field()):
        fp = Fingerprinter(sm, field)
        m = fp.matrix(sec + [s0_classical])
        assert np.array_equal(fp.matrix(sec + [s0_classical], layout=m.layout).values, m.values)
        ranks.append(la.rank(m.values, fp.field))
        assert m.block(("H", "Hdag")).shape[1] == len(sec) + 1
    assert ranks == [len(sec), len(sec)]
