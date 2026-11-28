"""Independent review, gate V2: Slavnov-Taylor identities of the derived tree-level action.

Statements (report 03_setup sec:brst, 03b sec:st-signs; ``SCRATCHPAD/st/test_signs.py``, ``test_b2.py``,
``test_b2_control.py``, ``test_blin.py``), each with a negative control that must fail:

* ``S(S_0) = s S_0 = Delta_1`` (eq:ST-S0, eq:Delta1) in ``D = 6`` and ``S(S_0) = 0`` in ``D = 4``, on every
  factor assignment of the 17 ghost-number-1 leg classes of ``SCRATCHPAD/st/test_signs.py``; ``Delta_1`` is the
  report's eq:V-delta, typed here. Controls: flipped sign of one pair of eq:ST; flipped BRST transformation.
* ``S(S_0 + T) - S(S_0) = s T`` for a source-free, non-invariant ``T`` (test_signs (b)).
* ``b^2 = 0`` in four dimensions on every field and source (test_b2); control: one flipped sign in eq:b-sources
  gives ``b^2 != 0`` on that source (test_b2_control).
* ``b_d G = S(S_0 + G) - S(S_0) - S(G)`` and ``S(S_0 + G) = Delta_1 + b_d G + S_q(G)`` for an even ``G`` containing
  every kind of source (test_blin, eq:m-S-exact).

All zero tests are exact evaluations with the independent evaluator (``tests/independent_support.py``) at random
kinematics and polarisations, class by class.
"""

from __future__ import annotations

import dataclasses
import itertools

import pytest
from independent_support import Evaluator, Sampler
from symbolica import Expression
from test_independent_v1_rules import (
    CHIRALITY,
    FACTOR_DIMS,
    HAND_PARAM,
    NONABELIAN,
    _adj,
    _generator,
    _rv,
    _spin,
    leg_pols,
    make_leg,
)

from hepkit_renorm import builtin_model
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import BRST
from hepkit_renorm.functional import Density, Derivation, FieldImage, LocalFunctional, Occurrence

GH1 = [['w', 'A'], ['w', 'A', 'A'], ['w', 'A', 'A', 'A'], ['w', 'Psib_Q', 'Psi_Q'], ['w', 'Psib_U', 'Psi_U', 'A'],
       ['w', 'Hd', 'H'], ['w', 'Hd', 'H', 'A'], ['w', 'Hd', 'H', 'Hd', 'H'], ['w', 'Psib_Q', 'Psi_U', 'Hd'],
       ['w', 'Psib_E', 'Psi_L', 'Hd'], ['w', 'w', 'Rb_Q', 'Psi_Q'], ['w', 'w', 'Psib_L', 'R_L'], ['w', 'w', 'rt', 'A'],
       ['w', 'w', 'w', 'z'], ['w', 'w', 'sd', 'H'], ['w', 'w', 'Hd', 's'], ['w', 'Psib_D', 'Psi_D', 'A', 'A']]
"""The 17 leg classes of SCRATCHPAD/st/test_signs.py (engine names; w, A, rt, z run over all gauge factors)."""


def engine_field(name: str, factor: str) -> str:
    """Package field of an engine leg name (``Rb_f`` = Rbar_f = Psi_f^*, ``R_f`` = Psibar_f^*, ``sd`` = sigma^dag =
    H^*, ``s`` = sigma = (H^dag)^*, ``rt`` = rho, ``z`` = zeta)."""
    fixed = {"Hd": "Hdag", "H": "H", "sd": "src_H", "s": "src_Hdag"}
    if name in fixed:
        return fixed[name]
    if name.startswith("Psib_"):
        return name[5:] + "bar"
    if name.startswith("Psi_"):
        return name[4:]
    if name.startswith("Rb_"):
        return "src_" + name[3:]
    if name.startswith("R_"):
        return "src_" + name[2:] + "bar"
    return {"w": "gh_", "A": "", "rt": "src_", "z": "src_gh_"}[name] + factor


def package_classes() -> list[tuple[str, ...]]:
    out = set()
    for legs in GH1:
        slots = [i for i, n in enumerate(legs) if n in ("w", "A", "rt", "z")]
        for facs in itertools.product("GWB", repeat=len(slots)):
            fmap = dict(zip(slots, facs))
            out.add(tuple(sorted(engine_field(n, fmap.get(i, "")) for i, n in enumerate(legs))))
    return sorted(out)


# ---------------------------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def model():
    return builtin_model("sm_one_family")


@pytest.fixture(scope="module")
def brst(model):
    return BRST(model)


@pytest.fixture(scope="module")
def params():
    smp = Sampler(202)
    vals = {n: smp.element() for n in ("g_s", "g", "gp", "y_u", "y_d", "y_e", "lambda", "m2", "xi")}
    vals.update({h: vals[p] for h, p in HAND_PARAM.items()})
    return vals


@pytest.fixture(scope="module")
def ev6(params):
    return Evaluator(6, params, FACTOR_DIMS)


@pytest.fixture(scope="module")
def ev4(params):
    return Evaluator(4, params, FACTOR_DIMS)


@pytest.fixture(scope="module")
def st_s0(brst):
    return brst.slavnov_taylor(brst.action)


# ---------------------------------------------------------------------------------------------- helpers
def class_value(model, ev, functional, cls, smp, order=None) -> int:
    """V of ``functional`` on the legs ``cls`` (in that order) at one random point."""
    legs = [make_leg(model, f, k) for k, f in enumerate(cls, start=1)]
    momenta = smp.momenta(len(cls), ev.D)
    pol = {}
    for leg in legs:
        pol.update(leg_pols(model, ev, leg, smp))
    return ev.value(functional.vertex_function(legs), pol, momenta)


def nonzero_classes(model, ev, functional, smp, classes=None) -> list[tuple[str, ...]]:
    """Leg classes (canonical order) on which ``functional`` does not vanish (``classes``: restrict to these)."""
    present = functional.leg_classes()
    if classes is not None:
        want = {tuple(sorted(c)) for c in classes}
        present = [c for c in present if tuple(sorted(c)) in want]
    return [c for c in present if class_value(model, ev, functional, c, smp) != 0]


def minus(a: LocalFunctional, b: LocalFunctional) -> LocalFunctional:
    return LocalFunctional.integrate(a - b)


def plus(*fs: LocalFunctional) -> LocalFunctional:
    total = fs[0]
    for f in fs[1:]:
        total = total + f
    return LocalFunctional.integrate(total)


def delta1_report(model, cls, legs):
    """eq:V-delta on (gh_r, fbar, f) with the ghost first; ``None`` if Delta_1 has no such vertex."""
    gh, fb, ff = cls
    r, f = gh[3:], ff
    if not (fb == f + "bar" and f in CHIRALITY and (r == "B" or r in NONABELIAN[f])):
        return None
    pc, pcb = (_rv("PL"), _rv("PR")) if CHIRALITY[f] == "L" else (_rv("PR"), _rv("PL"))
    lg, lb, lf = legs
    x, s = _rv("lor")(98_001), _rv("spin")(98_002)
    tg = _generator(model, list(zip(lb.labels, model.fields[lb.field].slots)), lf.labels, r, _adj(model, lg), f)
    sb, sf = _spin(model, lb), _spin(model, lf)
    return -Expression.I * tg * (_rv("gammahat")(x, sb, s) * _rv("p")(3, x) * pcb(s, sf)
                                 + _rv("gammahat")(x, sb, s) * _rv("p")(2, x) * pc(s, sf))


# ---------------------------------------------------------------------------------------------- S(S0)
def test_class_list_is_complete():
    classes = package_classes()
    assert len(GH1) == 17
    assert len(classes) > 60                    # every factor assignment of the 17 engine classes


def test_st_of_s0_is_delta1_in_six_and_zero_in_four_dimensions(model, brst, ev6, ev4, st_s0):
    smp = Sampler(31)
    classes = package_classes()
    assert nonzero_classes(model, ev4, st_s0, smp, classes) == []
    assert nonzero_classes(model, ev4, brst.breaking, smp, classes) == []
    seen_delta = 0
    for cls in classes:
        for functional in (st_s0, brst.breaking):
            gh = [x for x in cls if x.startswith("gh_")]
            fb = [x for x in cls if x.endswith("bar") and not x.startswith("src_")]
            ff = [x for x in cls if x not in gh + fb]
            ordered = gh + fb + ff
            if len(cls) == 3 and len(gh) == 1 and len(fb) == 1:
                legs = [make_leg(model, f, k) for k, f in enumerate(ordered, start=1)]
                exp = delta1_report(model, ordered, legs)
                if exp is not None:
                    seen_delta += 1
                    for _ in range(2):
                        momenta = smp.momenta(3, 6)
                        pol = {}
                        for leg in legs:
                            pol.update(leg_pols(model, ev6, leg, smp))
                        assert ev6.value(functional.vertex_function(legs), pol, momenta) == \
                            ev6.value(exp, pol, momenta) != 0, cls
                    continue
            if tuple(sorted(cls)) in {tuple(sorted(c)) for c in functional.leg_classes()}:
                assert class_value(model, ev6, functional, cls, smp) == 0, cls
    assert seen_delta == 2 * 3                  # (gh_r, Qbar, Q) for r = G, W, B, for S(S_0) and for s S_0


def test_st_sign_controls(model, brst, ev4):
    """Flipping the sign of any one pair of eq:ST (incl. the draft fermion and H^dag terms of 03b, which are the
    flips of the Psi and H^dag pairs) makes S(S_0) non-zero in four dimensions on the classes with that field."""
    smp = Sampler(32)
    classes = package_classes()
    for target in ("Q", "Qbar", "Hdag", "H", "G", "gh_W", "L"):
        pairs = [dataclasses.replace(pr, sign=-pr.sign) if pr.field == target else pr for pr in brst.pairs]
        assert pairs != list(brst.pairs)
        bad = brst.slavnov_taylor(brst.action, pairs)
        assert nonzero_classes(model, ev4, bad, smp, classes), target


@pytest.mark.parametrize("flip", ["Q", "Hdag", "gh_W"])
def test_brst_transformation_control(model, ev4, flip):
    """A wrong BRST transformation (one image negated) propagates consistently into L_gf + L_ext and the derived
    pairs, but S(S_0) = s S_0 no longer vanishes in four dimensions."""
    smp = Sampler(33)
    bad = BRST(model)
    good_s = BRST(model).s
    images = dict(good_s.images)
    img = images[flip]
    images[flip] = FieldImage(img.labels, -img.density)
    bad.__dict__["s"] = Derivation(1, images)
    assert nonzero_classes(model, ev4, bad.slavnov_taylor(bad.action), smp, package_classes())


def test_st_of_s0_plus_noninvariant_t(model, brst, ev6):
    """S(S_0 + T) - S(S_0) = s T for T = 3/10 Qbar gamma.B P_L Q + 7/10 Hdag H W.W (source free, not invariant)."""
    f = model.fields
    s, t, u = (tn.fresh(tn.SPINOR) for _ in range(3))
    i, w, j = tn.fresh(model.factors["G"].rep_space("F")), tn.fresh(model.factors["W"].rep_space("F")), \
        tn.fresh(model.factors["W"].rep_space("F"))
    mu, nu, a = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ), tn.fresh(model.factors["W"].adjoint_space())
    T = Density.monomial(f, Expression.num(3) / 10 * tn.gamma(mu, s, t) * tn.PL(t, u),
                         Occurrence("Qbar", (s, i, w)), Occurrence("Q", (u, i, w)), Occurrence("B", (mu,)))
    T = T + Density.monomial(f, Expression.num(7) / 10, Occurrence("Hdag", (j,)), Occurrence("H", (j,)),
                             Occurrence("W", (nu, a)), Occurrence("W", (nu, a)))
    T = LocalFunctional.integrate(T)
    assert not T.vanishes(model.simplifier())
    lhs = minus(brst.slavnov_taylor(plus(brst.action, T)), brst.slavnov_taylor(brst.action))
    diff = minus(lhs, brst.s(T))
    smp = Sampler(34)
    assert nonzero_classes(model, ev6, brst.s(T), smp)                 # s T is not zero ...
    assert nonzero_classes(model, ev6, diff, smp) == []                 # ... and equals the left side


# ---------------------------------------------------------------------------------------------- b^2 = 0
def density_nonzero(model, ev, density: Density, smp) -> list[tuple[str, ...]]:
    """Field classes on which a density with free labels does not vanish *as a density* (no integration by
    parts): its plane-wave coefficients at independent random momenta (no momentum conservation), with random
    polarisations on its free labels as well."""
    probe = object.__new__(LocalFunctional)        # the vertex-function code path, without the no-free-label check
    probe.fields, probe.expr = density.fields, density.expr
    free = {ev.label_key(lab): lab for lab in density.free_labels()}
    out = []
    for cls in sorted(set(density.field_content())):
        legs = [make_leg(model, f, k) for k, f in enumerate(cls, start=1)]
        momenta = {k: smp.vector(ev.D) for k in range(1, len(cls) + 1)}
        pol = {key: smp.vector(ev.dim(ev.space(lab))) for key, lab in free.items()}
        for leg in legs:
            pol.update(leg_pols(model, ev, leg, smp))
        if ev.value(probe.vertex_function(legs, momentum_conservation=False), pol, momenta) != 0:
            out.append(cls)
    return out


def b_squared(b: Derivation, name: str) -> tuple[Density, Density] | None:
    img = b.images.get(name)
    if img is None:
        return None
    return img.density, b(img.density)


def test_b_squared_vanishes_in_four_dimensions(model, brst, ev4):
    """b^2 X = 0 for every field and source X (b = the 4D limit of b_d = ST_{S_0})."""
    smp = Sampler(41)
    b = brst.b
    without = set()
    nontrivial = 0
    for name in model.fields:
        res = b_squared(b, name)
        if res is None:
            without.add(name)
            continue
        bx, bbx = res
        if density_nonzero(model, ev4, bx, smp):
            nontrivial += 1
        assert density_nonzero(model, ev4, bbx, smp) == [], name
    assert without == {"b_G", "b_W", "b_B", "gh_B"}       # s b = 0 and s omega^B = 0
    assert nontrivial == len(model.fields) - len(without)


@pytest.mark.parametrize("flip", ["src_H", "src_Qbar", "src_L", "src_Q", "src_Hdag", "src_W", "src_gh_G"])
def test_b_squared_control(model, brst, ev4, flip):
    """A flipped sign of b on one source (one equation of eq:b-sources) gives b^2 != 0 on that source."""
    smp = Sampler(42)
    images = dict(brst.b.images)
    img = images[flip]
    images[flip] = FieldImage(img.labels, -img.density)
    bad = Derivation(1, images)
    _, bbx = b_squared(bad, flip)
    assert density_nonzero(model, ev4, bbx, smp)
    _, good = b_squared(brst.b, flip)
    assert density_nonzero(model, ev4, good, smp) == []


# ---------------------------------------------------------------------------------------------- linearisation
def _g_with_every_source(model) -> LocalFunctional:
    """An even G containing every kind of source (after SCRATCHPAD/st/test_blin.py)."""
    f = model.fields
    G, W = model.factors["G"], model.factors["W"]
    a, b, c = (tn.fresh(G.adjoint_space()) for _ in range(3))
    mu = tn.fresh(tn.LORENTZ)
    s, t = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)
    i, j, k = (tn.fresh(G.rep_space("F")) for _ in range(3))
    w = tn.fresh(W.rep_space("F"))
    wi, wj = tn.fresh(W.rep_space("F")), tn.fresh(W.rep_space("F"))
    aw = tn.fresh(W.adjoint_space())
    nu, x = tn.fresh(tn.LORENTZ), tn.fresh(W.adjoint_space())
    q = lambda n, d: Expression.num(n) / d
    terms = [
        Density.monomial(f, q(2, 5) * tn.f(a, b, c), Occurrence("src_gh_G", (a,)), Occurrence("gh_G", (b,)),
                         Occurrence("gh_G", (c,))),
        Density.monomial(f, q(9, 10), Occurrence("src_G", (mu, a)), Occurrence("gh_G", (a,), (mu,))),
        Density.monomial(f, q(3, 5) * tn.I * tn.T(a, i, j) * tn.PL(s, t), Occurrence("src_Q", (s, i, w)),
                         Occurrence("gh_G", (a,)), Occurrence("Q", (t, j, w))),
        Density.monomial(f, q(1, 2) * tn.I * tn.PL(s, t), Occurrence("gh_B"), Occurrence("Ubar", (s, k)),
                         Occurrence("src_Ubar", (t, k))),
        Density.monomial(f, q(4, 5) * tn.I * tn.T(aw, wi, wj), Occurrence("src_H", (wi,)), Occurrence("gh_W", (aw,)),
                         Occurrence("H", (wj,))),
        Density.monomial(f, q(3, 10) * tn.I, Occurrence("gh_B"), Occurrence("Hdag", (wi,)),
                         Occurrence("src_Hdag", (wi,))),
        Density.monomial(f, q(7, 10), Occurrence("Hdag", (wj,)), Occurrence("H", (wj,)), Occurrence("W", (nu, x)),
                         Occurrence("W", (nu, x))),
    ]
    total = terms[0]
    for term in terms[1:]:
        total = total + term
    out = LocalFunctional.integrate(total)
    assert out.parity() == 0
    kinds = {f[n].kind for cls in out.field_content() for n in cls}
    assert "source" in kinds
    return out


def test_linearisation_property(model, brst, ev6):
    """b_d G = S(S0+G) - S(S0) - S(G) and S(S0+G) = Delta_1 + b_d G + S_q(G) (eq:m-S-exact), D = 6."""
    G = _g_with_every_source(model)
    sources = {n for cls in G.field_content() for n in cls if model.fields[n].kind == "source"}
    assert sources == {"src_gh_G", "src_G", "src_Q", "src_Ubar", "src_H", "src_Hdag"}
    st = brst.slavnov_taylor
    s0 = brst.action
    st_full = st(plus(s0, G))
    bdg = brst.b_d(G)
    lin = minus(minus(st_full, st(s0)), st(G))
    smp = Sampler(51)
    assert nonzero_classes(model, ev6, bdg, smp)
    assert nonzero_classes(model, ev6, minus(lin, bdg), smp) == []
    assert nonzero_classes(model, ev6, st(G), smp)                       # S_q(G) is not trivially zero
    exact = plus(brst.breaking, bdg, brst.quadratic_part(G))
    assert nonzero_classes(model, ev6, minus(st_full, exact), smp) == []
    # control: the linearised operator with one flipped source equation is not the linearisation
    images = dict(brst.b_d.images)
    img = images["src_Qbar"]
    images["src_Qbar"] = FieldImage(img.labels, -img.density)
    labels = model.slot_labels("Qbar")
    G2 = plus(G, LocalFunctional.integrate(Density.monomial(
        model.fields, Expression.I / 3, Occurrence("gh_B"), Occurrence("Qbar", labels), Occurrence("src_Qbar", labels))))
    bad = Derivation(1, images)(G2)
    lin2 = minus(minus(st(plus(s0, G2)), st(s0)), st(G2))
    assert nonzero_classes(model, ev6, minus(lin2, bad), smp)
