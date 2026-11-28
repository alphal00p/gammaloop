"""Checks of the broken-phase catalogue (gdocs.broken).

(a) Exactness of the shift, by an independent route. For every role of the order store (tree ``S0``, singular
    ``16 pi^2 eps deltaS1``, finite ``16 pi^2 DeltaS1``, breaking ``Delta1``) and every broken-phase field content,
    the fingerprint (exact plane-wave vertex function over GF(p), D = 6) of the sum of the broken-phase terms
    (coefficient of that part times operator) equals the vertex function of the *unbroken* store functional
    evaluated with the polarisations of the broken legs written in the unbroken fields (eq:bl-rot, eq:bp-H,
    eq:bl-fermions, typed here from the report, not taken from gdocs.broken_fields) and with ``a`` legs ``H`` and
    ``b`` legs ``H^dagger`` at zero momentum with polarisation ``H_v``, weighted by ``1/(a! b!)`` (Taylor expansion of
    ``S[H_v + phi]``). The unbroken side never uses the substitution code of the catalogue.
(b) Statements of report Section 5.4: tree-level masses and mixing, the photon / Z-photon / gluon zero-momentum
    check (finite counterterm against the fermion-loop value quoted in the report), the tadpole and Goldstone mass
    counterterms, the neutral field-renormalisation matrix, the finite mixing and longitudinal terms.
(c) ``v = 0`` (before the vacuum condition) gives back the unbroken store functionals in the rotated basis
    (same independent route without constant legs).
Plus: the written ``data/broken.json`` carries the coefficients of the catalogue; operators of one content are
independent.

Run: ``SYMBOLICA_LICENSE=... PYTHONPATH=../../python/hepkit_renorm:. python -m pytest tests/test_broken.py`` from
gamma5_docs/tools (``GDOCS_SITE=<dir>`` checks another site directory). About two minutes.
"""

from __future__ import annotations

import json
import math
import os
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("symbolica")

from symbolica import Expression  # noqa: E402

from hepkit_renorm import Density, LocalFunctional, Occurrence  # noqa: E402
from hepkit_renorm import tensors as tn  # noqa: E402
from hepkit_renorm.exact import default_field  # noqa: E402
from hepkit_renorm.fingerprint import Fingerprinter  # noqa: E402

from gdocs import broken as br  # noqa: E402
from gdocs import broken_fields as bf  # noqa: E402
from gdocs import mathview as mv  # noqa: E402

SITE = Path(os.environ.get("GDOCS_SITE", br.SITE))
D = 6
NSAMP = 2


# ---------------------------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def ctx():
    return br.Ctx(SITE, log=lambda s: None)


@pytest.fixture(scope="module")
def terms(ctx):
    return br.catalogue(ctx, log=lambda s: None)


@pytest.fixture(scope="module")
def by_id(terms):
    return {t.id: t for t in terms}


@pytest.fixture(scope="module")
def field_():
    return default_field()


@pytest.fixture(scope="module")
def values(field_):
    return br.broken_values(field_, seed=3)


def _fps(field_, values, v_zero=False):
    vals = dict(values)
    if v_zero:
        vals["v"] = 0
        vals["m2"] = (values["m2"] + 12345) % field_.p     # m^2 independent of v
    fpb = Fingerprinter(bf.BrokenModel(), field_, D, values=vals)
    model = bf.unbroken_model()
    uvals = {k: vals[k] for k in ("g_s", "g", "gp", "y_u", "y_d", "y_e", "lambda", "m2", "xi")}
    fpu = Fingerprinter(model, field_, D, values=uvals)
    return fpb, fpu, vals


# ---------------------------------------------------------------------------------------------- the independent route
def _component_map(F, vals):
    """Unbroken-field polarisations of a broken leg (report eq:bl-rot, eq:bp-H, eq:bl-fermions):
    ``{broken field: [(unbroken field, fn(e) -> polarisation array)]}`` over GF(p) (object arrays)."""
    p = F.p
    i, r2 = F.i, F.sqrt(2)
    inv2 = pow(r2, p - 2, p)
    c, s = vals["c_W"], vals["s_W"]
    out = {}
    adj = {"Wp": [inv2, i * inv2 % p, 0], "Wm": [inv2, (-i * inv2) % p, 0], "Z": [0, 0, c], "A": [0, 0, s]}
    bco = {"Z": (-s) % p, "A": c}

    def w_pol(coef):
        return lambda e: (e.astype(object)[..., None] * np.array(coef, dtype=object)) % p

    def scale(k):
        return lambda e: (e.astype(object) * k) % p

    for pre in ("", "gh_", "agh_", "b_", "src_", "src_gh_"):
        for x in ("Wp", "Wm", "Z", "A"):
            nm = pre + x
            lst = [(pre + "W", w_pol(adj[x]))]
            if x in bco:
                lst.append((pre + "B", scale(bco[x])))
            out[nm] = lst
        out[pre + "G"] = [(pre + "G", lambda e: e.astype(object) % p)]

    def fund(vec):
        return lambda e: (e.astype(object)[..., None] * np.array(vec, dtype=object)) % p

    out["h"] = [("H", fund([0, inv2])), ("Hdag", fund([0, inv2]))]
    out["G0"] = [("H", fund([0, i * inv2 % p])), ("Hdag", fund([0, (-i * inv2) % p]))]
    out["Gp"] = [("H", fund([1, 0]))]
    out["Gm"] = [("Hdag", fund([1, 0]))]
    for k in (1, 2):
        out[f"src_H_{k}"] = [("src_H", fund([1, 0] if k == 1 else [0, 1]))]
        out[f"src_Hdag_{k}"] = [("src_Hdag", fund([1, 0] if k == 1 else [0, 1]))]
    return out


def _fermion_map(fpu):
    """Fermion legs: column fields ``Psi = P psi``, row fields ``Psibar = psibar P'`` (eq:bl-fermions)."""
    p = fpu.field.p
    PL, PR = (np.asarray(x, dtype=object) for x in fpu.rep.projectors)
    proj = {"PL": PL, "PR": PR}
    other = {"PL": "PR", "PR": "PL"}
    # (broken field) -> [(unbroken field, projector, SU(2) component or None)]
    table = {"u": [("Q", "PL", 0), ("U", "PR", None)], "up": [("Q", "PR", 0), ("U", "PL", None)],
             "d": [("Q", "PL", 1), ("D", "PR", None)], "dp": [("Q", "PR", 1), ("D", "PL", None)],
             "e": [("L", "PL", 1), ("E", "PR", None)], "ep": [("L", "PR", 1), ("E", "PL", None)],
             "nu": [("L", None, 0)]}
    out = {}
    for f, lst in table.items():
        col, row = [], []
        for X, P, comp in lst:
            def mk(P=P, comp=comp, bar=False):
                def fn(e):
                    e = e.astype(object)
                    if P is not None:
                        if bar:   # (psibar P')_s = psibar_t P'_{t s}
                            e = np.moveaxis(np.tensordot(e, proj[other[P]], axes=([1], [0])), -1, 1) % p
                        else:     # (P psi)_s = P_{s t} psi_t
                            e = np.moveaxis(np.tensordot(proj[P], e, axes=([1], [1])), 0, 1) % p
                    if comp is not None:
                        v = np.zeros(2, dtype=object)
                        v[comp] = 1
                        e = (e[..., None] * v) % p
                    return e
                return fn
            col.append((X, mk()))
            row.append((X + "bar", mk(bar=True)))
        out[f] = col
        out[f + "bar"] = row
    # fermion sources: SU(2) components (src_Q_1 is the upper component of Rbar_Q)
    for f in "QUDLE":
        for suff in ("", "bar"):
            for k in ((1, 2) if f in "QL" else (None,)):
                nm = f"src_{f}{suff}" + (f"_{k}" if k else "")

                def mk2(k=k):
                    def fn(e):
                        e = e.astype(object)
                        if k is None:
                            return e % p
                        v = np.zeros(2, dtype=object)
                        v[k - 1] = 1
                        return (e[..., None] * v) % p
                    return fn
                out[nm] = [(f"src_{f}{suff}", mk2())]
    return out


def _unbroken_value(fpu, F_unbroken, content, smp, vals, maps, constant=True):
    """The broken vertex function of ``content`` from the unbroken functional (module docstring, (a))."""
    p = fpu.field.p
    n = smp.n
    legs = list(content)
    options = [maps[x] for x in legs]
    total = np.zeros(n, dtype=object)
    hv = np.zeros(2, dtype=object)
    hv[1] = vals["v"] * pow(fpu.field.sqrt(2), p - 2, p) % p
    nmax = 4 - len(legs) if constant else 0
    import itertools
    for choice in itertools.product(*options):
        names = [c[0] for c in choice]
        pols = [fn(smp.polarisations[k]) for k, (_, fn) in enumerate(choice)]
        moms = np.stack([smp.momenta[k + 1] for k in range(len(legs))], axis=1).astype(object)
        for a in range(nmax + 1):
            for b in range(nmax + 1 - a):
                cls = names + ["H"] * a + ["Hdag"] * b
                if not cls:
                    continue
                m = np.concatenate([moms, np.zeros((n, a + b, D), dtype=object)], axis=1) if a + b else moms
                ps = pols + [np.tile(hv, (n, 1))] * (a + b)
                s = fpu.explicit_samples(cls, np.asarray(m, dtype=np.int64),
                                         [np.asarray(x, dtype=np.int64) for x in ps])
                val = fpu.functional(F_unbroken, s).astype(object)
                w = pow(math.factorial(a) * math.factorial(b), p - 2, p)
                total = (total + val * w) % p
    return total


def _broken_functional(terms, part, raw=False, v_zero=False):
    """Sum over the terms of (coefficient of ``part``) times operator (``finite``: the four-dimensional form)."""
    E_ = mv.parse("e")
    g, sw = bf.unbroken_model().param("g"), bf.SW
    total = Expression.num(0)
    for t in terms:
        parts = ("pole", "evanescent") if part == "pole" else (part,)
        for q in parts:
            src = t.raw if raw else t.coeffs
            if q not in src:
                continue
            c = src[q].replace(E_, g * sw)
            if v_zero:
                c = c.replace(bf.V, Expression.num(0))
            op = t.op4 if (q == "finite" and t.op4 is not None) else t.op
            total = total + c * op.expr
    return LocalFunctional(bf.BT, total)


ROLE = {"tree": "S0", "pole": "dS", "finite": "DS", "breaking": "Delta1"}


def _check_role(ctx, terms, field_, values, part, *, v_zero=False):
    fpb, fpu, vals = _fps(field_, values, v_zero=v_zero)
    store = getattr(ctx.c, ROLE[part])
    Fb = _broken_functional(terms, part, raw=v_zero, v_zero=v_zero)
    maps = _component_map(fpu.field, vals)
    maps.update(_fermion_map(fpu))
    # every content of the broken side and of the shifted store functional (the latter only to enumerate classes)
    contents = {tuple(c) for c in Fb.leg_classes()}
    contents |= {tuple(c) for c in bf.shift(store, constant=not v_zero).leg_classes() if c}
    bad = []
    for c in sorted(contents):
        smp = fpb.samples(list(c), NSAMP, salt=5)
        lhs = fpb.functional(Fb, smp).astype(object) % fpb.field.p
        rhs = _unbroken_value(fpu, store, c, smp, vals, maps, constant=not v_zero)
        if not np.array_equal(lhs % fpb.field.p, rhs % fpb.field.p):
            bad.append(c)
    return len(contents), bad


@pytest.mark.parametrize("part", ["tree", "pole", "finite", "breaking"])
def test_shift_exact(ctx, terms, field_, values, part):
    n, bad = _check_role(ctx, terms, field_, values, part)
    assert n > 0
    assert not bad, f"{part}: broken terms differ from the shifted store functional on {bad[:8]}"


@pytest.mark.parametrize("part", ["tree", "pole", "finite", "breaking"])
def test_v_zero_gives_unbroken(ctx, terms, field_, values, part):
    """(c): at v = 0 (before the vacuum condition) the catalogue is the unbroken store functional, rotated."""
    n, bad = _check_role(ctx, terms, field_, values, part, v_zero=True)
    assert n > 0
    assert not bad, f"{part} at v = 0: differs on {bad[:8]}"


# ---------------------------------------------------------------------------------------------- (b) the report's statements
P = bf.unbroken_model().param


def _ref(monos):
    """A broken-phase functional from ``[(coefficient, [(field, slots, derivatives)])]``."""
    dens = Density(bf.BT)
    for coeff, occs in monos:
        dens = dens + Density.monomial(bf.BT, coeff, *[Occurrence(f, tuple(sl), tuple(dv)) for f, sl, dv in occs])
    return LocalFunctional.integrate(dens)


MU, NU = tn.lor(1), tn.lor(2)
S1, S2 = tn.spin(1), tn.spin(2)


def _same(fp, A, B, content):
    smp = fp.samples(list(content), 3, salt=9)
    return np.array_equal(fp.functional(A, smp) % fp.field.p, fp.functional(B, smp) % fp.field.p)


def _part_functional(t, part):
    E_ = mv.parse("e")
    c = t.coeffs[part].replace(E_, P("g") * bf.SW)
    op = t.op4 if (part == "finite" and t.op4 is not None) else t.op
    return LocalFunctional(bf.BT, (c * op.expr))


def _terms_on(terms, content, part, no_derivative=False):
    out = Expression.num(0)
    for t in terms:
        if no_derivative and br._has_derivative(t.op):
            continue
        if tuple(t.content) == tuple(content) and part in t.coeffs:
            out = out + _part_functional(t, part).expr
    return LocalFunctional(bf.BT, out)


@pytest.fixture(scope="module")
def fp4(field_, values):
    return Fingerprinter(bf.BrokenModel(), field_, 4, values=values)


@pytest.fixture(scope="module")
def fp6(field_, values):
    return Fingerprinter(bf.BrokenModel(), field_, 6, values=values)


g_, gp_, gs_, lam_ = P("g"), P("gp"), P("g_s"), P("lambda")
yu, yd, ye = P("y_u"), P("y_d"), P("y_e")
v_, c_, s_, r2 = bf.V, bf.CW, bf.SW, bf.RT2
MW2 = g_ ** 2 * v_ ** 2 / 4
MZ2 = MW2 / c_ ** 2
MW, MZ = g_ * v_ / 2, g_ * v_ / (2 * c_)
Mh2 = 2 * lam_ * v_ ** 2
m_u, m_d, m_e = yu * v_ / r2, yd * v_ / r2, ye * v_ / r2
Ytil = 3 * yu ** 2 + 3 * yd ** 2 + ye ** 2
e_ = g_ * s_
half = Expression.num(1) / 2


def test_tree_masses(terms, fp6):
    """eq:bl-masses and eq:bl-L2V/L2S/L2F: M_W^2 W+.W- + 1/2 M_Z^2 Z.Z - 1/2 M_h^2 h^2 - m_f fbar f; no tree mass
    for the photon, the gluon, the Goldstone fields; no tree tadpole (eq:bl-vev)."""
    checks = {
        ("Wp", "Wm"): [(MW2, [("Wp", [MU], []), ("Wm", [MU], [])])],
        ("Z", "Z"): [(MZ2 / 2, [("Z", [MU], []), ("Z", [MU], [])])],
        ("h", "h"): [(-Mh2 / 2, [("h", [], []), ("h", [], [])])],
    }
    for content, mono in checks.items():
        assert _same(fp6, _tree_mass_part(terms, content), _ref(mono), content), content
    for f, m in (("u", m_u), ("d", m_d), ("e", m_e)):
        cf = tuple(sorted((f, f + "bar"), key=bf.BT.order))
        col = (tn.fund(tn.factor_symbol("G"), 1),) if f in "ud" else ()
        ref = _ref([(-m * tn.delta(S1, S2), [(f + "bar", (S1,) + col, []), (f, (S2,) + col, [])])])
        assert _same(fp6, _tree_mass_part(terms, cf), ref, cf), f
    for content in (("A", "A"), ("G", "G"), ("G0", "G0"), ("Gp", "Gm"), ("h",)):
        assert _tree_mass_part(terms, content).is_zero(), content


def _tree_mass_part(terms, content):
    """Tree-level part of the terms of ``content`` without derivatives (mass terms)."""
    out = Expression.num(0)
    for t in terms:
        if tuple(t.content) == tuple(content) and "tree" in t.coeffs and not br._has_derivative(t.op):
            out = out + _part_functional(t, "tree").expr
    return LocalFunctional(bf.BT, out)


def test_tree_mixing(terms, fp6):
    """eq:bl-L2mix: M_Z Z.dG0 + i M_W (W-.dG+ - W+.dG-)."""
    for content, mono in {
        ("Z", "G0"): [(MZ, [("Z", [MU], []), ("G0", [], [MU])])],
        ("Wm", "Gp"): [(tn.I * MW, [("Wm", [MU], []), ("Gp", [], [MU])])],
        ("Wp", "Gm"): [(-tn.I * MW, [("Wp", [MU], []), ("Gm", [], [MU])])],
    }.items():
        c = tuple(sorted(content, key=bf.BT.order))
        assert _same(fp6, _terms_on(terms, c, "tree"), _ref(mono), c), content
    assert _terms_on(terms, ("A", "G0"), "tree").is_zero()


def _coeff_of(fp, terms, content, part, mono):
    """The number k with (terms of ``content``, part) = k * mono on the samples (both non-zero)."""
    A = _terms_on(terms, content, part, no_derivative=True)
    B = _ref(mono)
    smp = fp.samples(list(content), 2, salt=13)
    a = fp.functional(A, smp) % fp.field.p
    b = fp.functional(B, smp) % fp.field.p
    pp = fp.field.p
    k = int(a[0]) * pow(int(b[0]), pp - 2, pp) % pp
    assert np.array_equal(a, (b * k) % pp), f"{content} {part} is not proportional to the reference monomial"
    return k


def _val(fp, expr):
    """Value of a scalar expression at the fingerprinter's parameter values (mod p)."""
    pp = fp.field.p
    vals = fp._values
    e = expr.replace(mv.parse("e"), g_ * s_)
    for name in ("g_s", "g", "gp", "y_u", "y_d", "y_e", "lambda", "m2", "xi", "v", "c_W", "s_W", "sqrt2"):
        sym = bf.BrokenModel().param(name)
        e = e.replace(sym, Expression.num(int(vals[sym.get_name()])))
    e = e.replace(Expression.I, Expression.num(int(fp.field.i))).expand()
    f = Fraction(e.to_canonical_string())
    return f.numerator * pow(f.denominator, pp - 2, pp) % pp


def test_zero_momentum_massless(terms, fp4):
    """Section 5.4.6 (eq:bl-M, eq:bl-dMVV, eq:bl-massless): the finite mass-like terms are 1/2 DeltaM^2_VV' A^V.A^V'
    with DeltaM^2_VV' = -2 sum_f N_c m_f^2 (v v' + 3 a a'); with the fermion-loop value quoted by the report,
    c = 2 v v' + 26/3 a a' per fermion and unit N_c, the photon, Z-photon and gluon two-point functions at zero
    momentum vanish."""
    pp = fp4.field.p
    T3 = {"u": half, "d": -half, "e": -half}
    q = {"u": Expression.num(2) / 3, "d": -Expression.num(1) / 3, "e": -Expression.num(1)}
    Nc = {"u": 3, "d": 3, "e": 1}
    m2 = {"u": m_u ** 2, "d": m_d ** 2, "e": m_e ** 2}
    vV = {"A": lambda f: e_ * q[f], "Z": lambda f: g_ / (2 * c_) * (T3[f] - 2 * s_ ** 2 * q[f])}
    aV = {"A": lambda f: Expression.num(0), "Z": lambda f: g_ / (2 * c_) * T3[f]}
    for V1, V2 in (("A", "A"), ("Z", "A"), ("Z", "Z")):
        sym = Expression.num(1) / 2 if V1 == V2 else Expression.num(1)
        k = _coeff_of(fp4, terms, (V1, V2), "finite", [(sym, [(V1, [MU], []), (V2, [MU], [])])])
        dM = sum((-2 * Nc[f] * m2[f] * (vV[V1](f) * vV[V2](f) + 3 * aV[V1](f) * aV[V2](f)) for f in "ude"),
                 Expression.num(0))
        assert k == _val(fp4, dM), (V1, V2)
        loop = sum((Nc[f] * m2[f] * (2 * vV[V1](f) * vV[V2](f) + Expression.num(26) / 3 * aV[V1](f) * aV[V2](f))
                    for f in "ude"), Expression.num(0))
        if V1 == "A":   # photon-photon and Z-photon: loop plus counterterm vanish
            assert (_val(fp4, loop) + k) % pp == 0, (V1, V2)
    # gluon: Tr(T^A T^B) = delta/2 replaces N_c; c = 2 v v' with v = g_s, counterterm -g_s^2 (m_u^2 + m_d^2) G.G / 2
    k = _coeff_of(fp4, terms, ("G", "G"), "finite",
                  [(tn.delta(tn.adj(tn.factor_symbol("G"), 1), tn.adj(tn.factor_symbol("G"), 2)) / 2,
                    [("G", [MU, tn.adj(tn.factor_symbol("G"), 1)], []), ("G", [MU, tn.adj(tn.factor_symbol("G"), 2)],
                                                                          [])])])
    assert k == _val(fp4, -gs_ ** 2 * (m_u ** 2 + m_d ** 2))
    loop = half * 2 * gs_ ** 2 * (m_u ** 2 + m_d ** 2)     # sum over u, d of (1/2) c m_q^2, c = 2 v v'
    assert (_val(fp4, loop) + k) % pp == 0      # 1/2 Delta M^2 G.G contributes Delta M^2 = k to Gamma


def test_tadpole_and_goldstone_masses(terms, by_id, fp6):
    """eq:bl-dS-scalar, eq:bl-dT: -v deltaT h with deltaT = v^2 (6 lambda^2 - lambda C_H - Y_4 + 3/16 G_4); the
    Goldstone mass counterterm is -deltaT (G+G- + 1/2 G0^2); h^2: -1/2 (M_h^2 dZ_H + dM_h^2 + deltaT)."""
    CH = Expression.num(3) / 4 * g_ ** 2 + Expression.num(1) / 4 * gp_ ** 2
    Y4 = 3 * yu ** 4 + 3 * yd ** 4 + ye ** 4
    G4 = 3 * g_ ** 4 + 2 * g_ ** 2 * gp_ ** 2 + gp_ ** 4
    dT = v_ ** 2 * (6 * lam_ ** 2 - lam_ * CH - Y4 + Expression.num(3) / 16 * G4)
    assert not any("tree" in t.coeffs for t in terms if tuple(t.content) == ("h",))
    assert _same(fp6, _terms_on(terms, ("h",), "pole"), _ref([(-v_ * dT, [("h", [], [])])]), ("h",))
    assert _same(fp6, _mass_pole(terms, ("Gp", "Gm")), _ref([(-dT, [("Gp", [], []), ("Gm", [], [])])]),
                 ("Gp", "Gm"))
    assert _same(fp6, _mass_pole(terms, ("G0", "G0")), _ref([(-dT / 2, [("G0", [], []), ("G0", [], [])])]),
                 ("G0", "G0"))
    # h^2: dM_h^2 = M_h^2 (dlambda/lambda + dZ_H), dlambda + 2 lambda dZ_H = 12 lambda^2 - 2 lambda C_H - Y4 + 3/16 G4
    dZH = 2 * CH - Ytil
    dlam = 12 * lam_ ** 2 - 2 * lam_ * CH - Y4 + Expression.num(3) / 16 * G4 - 2 * lam_ * dZH
    dMh2 = Mh2 * (dlam / lam_ + dZH)
    assert _same(fp6, _mass_pole(terms, ("h", "h")),
                 _ref([(-(Mh2 * dZH + dMh2 + dT) / 2, [("h", [], []), ("h", [], [])])]), ("h", "h"))


def _mass_pole(terms, content):
    out = Expression.num(0)
    for t in terms:
        if tuple(t.content) == tuple(content) and "pole" in t.coeffs and not br._has_derivative(t.op):
            out = out + _part_functional(t, "pole").expr
    return LocalFunctional(bf.BT, out)


def test_neutral_field_renormalisation(by_id):
    """eq:bl-dZZA: 16 pi^2 eps dZ_ZZ = (33 g^2 c^2 - 43 g'^2 s^2)/18, dZ_gamma gamma = -5/9 e^2,
    dZ_Z gamma = s c (33 g^2 + 43 g'^2)/18 (kinetic terms -1/4 [dZ_ZZ Z.Z + 2 dZ_Zg Z.A + dZ_gg A.A])."""
    k = {c: next(t for t in by_id.values() if tuple(t.content) == c and t.parent_ids[0] == "kin-W")
         for c in (("Z", "Z"), ("Z", "A"), ("A", "A"))}
    assert bf.equal(k[("Z", "Z")].coeffs["pole"], (33 * g_ ** 2 * c_ ** 2 - 43 * gp_ ** 2 * s_ ** 2) / 18)
    assert bf.equal(k[("A", "A")].coeffs["pole"], -Expression.num(5) / 9 * e_ ** 2)
    assert bf.equal(k[("Z", "A")].coeffs["pole"], s_ * c_ * (33 * g_ ** 2 + 43 * gp_ ** 2) / 18)
    for t in k.values():   # the operators are the transverse kinetic forms with tree coefficient 1 (Z, A)
        assert "tree" not in t.coeffs or bf.equal(t.coeffs["tree"], 1)


def test_mass_and_mixing_counterterms(terms, fp6, fp4):
    """eq:bl-dS-mass and eq:bl-dS-mix, and the finite mixing eq:bl-mix = -Ytil/3 times the tree mixing."""
    CH2 = Expression.num(3) / 2 * g_ ** 2 + half * gp_ ** 2     # 2 C_H
    checks = [
        (("Wp", "Wm"), "pole", MW2 * (CH2 - Ytil - 4 * g_ ** 2), [("Wp", [MU], []), ("Wm", [MU], [])]),
        (("Z", "Z"), "pole", MZ2 / 2 * (CH2 - Ytil - 4 * g_ ** 2 * c_ ** 2), [("Z", [MU], []), ("Z", [MU], [])]),
        (("Z", "A"), "pole", -2 * g_ ** 2 * s_ * c_ * MZ2, [("Z", [MU], []), ("A", [MU], [])]),
        (("Z", "G0"), "pole", MZ * (CH2 - Ytil - 2 * g_ ** 2 * c_ ** 2), [("Z", [MU], []), ("G0", [], [MU])]),
        (("A", "G0"), "pole", -2 * g_ ** 2 * s_ * c_ * MZ, [("A", [MU], []), ("G0", [], [MU])]),
        (("Wm", "Gp"), "pole", tn.I * MW * (CH2 - Ytil - 2 * g_ ** 2), [("Wm", [MU], []), ("Gp", [], [MU])]),
        (("Z", "G0"), "finite", -Ytil / 3 * MZ, [("Z", [MU], []), ("G0", [], [MU])]),
        (("Wm", "Gp"), "finite", -Ytil / 3 * tn.I * MW, [("Wm", [MU], []), ("Gp", [], [MU])]),
        (("Wp", "Wm"), "finite", -g_ ** 2 / 2 * (3 * m_u ** 2 + 3 * m_d ** 2 + m_e ** 2),
         [("Wp", [MU], []), ("Wm", [MU], [])]),
    ]
    for content, part, coeff, mono in checks:
        c = tuple(sorted(content, key=bf.BT.order))
        got = Expression.num(0)
        for t in terms:
            if tuple(t.content) == c and part in t.coeffs and not t.elementary and not br._is_ev(t) and \
                    t.kind != "unchanged":
                got = got + _part_functional(t, part).expr
        fp = fp4 if part == "finite" else fp6       # the finite counterterm is four-dimensional (bars omitted)
        assert _same(fp, LocalFunctional(bf.BT, got), _ref([(coeff, mono)]), c), (content, part)


def test_longitudinal_terms(by_id):
    """eq:bl-L in the mass basis: (dW+)(dW-) 2/3 g^2, (dZ)^2 (3 g^2 c^2 + 5 g'^2 s^2)/9, (dA)^2 8/9 e^2 (units lf)."""
    lw = {tuple(t.content): t for t in by_id.values() if t.parent_ids[0] == "long-W"}
    # the operators are written (dW+)(dW-) and 1/2 (dZ)^2, 1/2 (dA)^2, (dZ)(dA)
    assert bf.equal(lw[("Wp", "Wm")].coeffs["finite"], Expression.num(2) / 3 * g_ ** 2)
    assert bf.equal(lw[("Z", "Z")].coeffs["finite"], 2 * (3 * g_ ** 2 * c_ ** 2 + 5 * gp_ ** 2 * s_ ** 2) / 9)
    assert bf.equal(lw[("A", "A")].coeffs["finite"], 2 * Expression.num(8) / 9 * e_ ** 2)
    assert bf.equal(lw[("Z", "A")].coeffs["finite"], 2 * e_ * (3 * g_ * c_ - 5 * gp_ * s_) / 9)
    for c in (("Wp", "Wm"), ("Z", "Z"), ("A", "A"), ("Z", "A")):
        num, _ = br._num_of(br.mv.operator_parts(lw[c].op, table=bf.BT)[0][0])
        assert num == (Expression.num(1) / 2 if c[0] == c[1] else 1)


def test_fermion_currents(by_id):
    """eq:bl-Am and eq:bl-dS-ferm: tree couplings e q_f to the photon and g/c_W (T3 - s_W^2 q) on P_L,
    -g/c_W s_W^2 q on P_R to the Z; g/sqrt2 for W."""
    def term(parent, content):
        return next(t for t in by_id.values() if t.parent_ids[0] == parent and tuple(t.content) == content)

    assert bf.equal(term("vertex-Q-W", ("A", "u", "ubar")).coeffs["tree"], Expression.num(2) / 3 * e_)
    assert bf.equal(term("vertex-U-B", ("A", "u", "ubar")).coeffs["tree"], Expression.num(2) / 3 * e_)
    assert bf.equal(term("vertex-Q-W", ("Z", "u", "ubar")).coeffs["tree"],
                    g_ / c_ * (half - Expression.num(2) / 3 * s_ ** 2))
    assert bf.equal(term("vertex-U-B", ("Z", "u", "ubar")).coeffs["tree"],
                    -g_ / c_ * Expression.num(2) / 3 * s_ ** 2)
    assert bf.equal(term("vertex-L-W", ("Wp", "e", "nubar")).coeffs["tree"], g_ / r2)
    assert "tree" not in term("vertex-L-W", ("A", "nu", "nubar")).coeffs


def test_written_file(by_id):
    data = json.loads((SITE / "data" / "broken.json").read_text())
    assert "mock" not in data
    ids = [t["id"] for t in data["terms"]]
    assert len(ids) == len(set(ids)) == len(by_id)
    unb = {t["id"] for t in json.loads((SITE / "data" / "unbroken.json").read_text())["terms"]}
    for t in data["terms"]:
        bt = by_id[t["id"]]
        assert set(t["types"]) == set(bt.coeffs)
        for part, s in t["data"]["values"].items():
            assert bf.equal(br.parse_values(s), bt.coeffs[part]), (t["id"], part)
        assert t["parents"] and set(t["parents"]) <= unb
        if t["kind"] == "unchanged" and bt.identical:
            assert t["parents"] == [t["id"]]
    kinds = {k for k, _ in br.KINDS}
    assert {t["kind"] for t in data["terms"]} <= kinds


def test_independent_operators(terms):
    assert br.check_independence(terms) > 0


def test_parent_sums_and_chirality_groups(ctx, terms):
    """Physics review: the coefficient of every term is the sum over its parents of parent coefficient times
    factor; for a term that collects the two chirality projections of one bilinear this holds for each projection
    separately (and the remark says so), not for the sum over all parents."""
    merged = 0
    for t in terms:
        groups = [pl for _, pl in t.groups] if t.merged else [t.parents]
        if t.merged:
            merged += 1
            assert len(groups) >= 2 and sum(len(g) for g in groups) == len(t.parents), t.id
        for part, c in t.coeffs.items():
            for g in groups:
                s = Expression.num(0)
                for pid, lam, _ in g:
                    s = s + lam * ctx.coeffs[pid].get(part, Expression.num(0))
                assert bf.is_zero(bf.vev(s) - c), (t.id, part)
    assert merged >= 20
    data = {x["id"]: x for x in json.loads((SITE / "data" / "broken.json").read_text())["terms"]}
    for t in terms:
        rem = data[t.id]["remarks"] or ""
        if t.merged:
            assert "two chirality projections" in rem and "sum of the parents" not in rem, t.id
