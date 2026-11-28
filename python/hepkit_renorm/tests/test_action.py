"""The regularised classical action of the one-family Standard Model (report 03_setup sec:action and
sec:feynman-rules): content, Feynman rules, hermiticity, rigid invariance, evanescent part."""

import pytest
from symbolica import Expression, S

from hepkit_renorm import Density, LocalFunctional, Occurrence, classical_sectors, standard_legs
from hepkit_renorm import tensors as tn
from hepkit_renorm.fingerprint import is_zero

I = tn.I
_w = S("w_")

FERMIONS = ("Q", "U", "D", "L", "E")
EXPECTED_CLASSES = {
    # gauge kinetic terms
    ("G", "G"), ("G", "G", "G"), ("G", "G", "G", "G"), ("W", "W"), ("W", "W", "W"), ("W", "W", "W", "W"), ("B", "B"),
    # fermion kinetic terms and gauge vertices (report eq:V-ffv)
    *{(f, f + "bar") for f in FERMIONS},
    *{(r, f, f + "bar") for f, rs in {"Q": "GWB", "U": "GB", "D": "GB", "L": "WB", "E": "B"}.items() for r in rs},
    # Higgs sector (eq:V-ssv, eq:V-ssss)
    ("H", "Hdag"), ("W", "H", "Hdag"), ("B", "H", "Hdag"), ("W", "W", "H", "Hdag"), ("W", "B", "H", "Hdag"),
    ("B", "B", "H", "Hdag"), ("H", "H", "Hdag", "Hdag"),
    # Yukawa terms and their conjugates (eq:V-yuk)
    ("Qbar", "U", "Hdag"), ("Q", "Ubar", "H"), ("Qbar", "D", "H"), ("Q", "Dbar", "Hdag"), ("Lbar", "E", "H"),
    ("L", "Ebar", "Hdag"),
}


def canon(classes):
    return {tuple(sorted(c)) for c in classes}


@pytest.fixture(scope="module")
def sectors(sm):
    return classical_sectors(sm)


def test_field_content_and_monomial_counts(sm, s0_classical, sectors):
    assert canon(s0_classical.leg_classes()) == canon(EXPECTED_CLASSES)
    assert {k: len(v) for k, v in sectors.items()} == {"gauge": 10, "fermion": 15, "scalar": 8, "potential": 2,
                                                       "yukawa": 6}
    assert len(s0_classical) == 41


def _digest_multiplet(name: str) -> str:
    bar = name.endswith("~")
    base = name.rstrip("~")
    base = {"W1": "W", "W2": "W", "W3": "W", "Q1": "Q", "Q2": "Q", "L1": "L", "L2": "L", "H1": "H", "H2": "H"}.get(
        base, base)
    if base in FERMIONS:
        return base + ("bar" if bar else "")
    if base == "H":
        return "Hdag" if bar else "H"
    return base


def test_vertex_families_agree_with_the_hand_written_model(s0_classical, golden_rules):
    """The multiplet-level families of the 174 hand-written rules without ghosts, sources and DELTA."""
    families = set()
    for rule in golden_rules["vertex_rules"]:
        if any(p.startswith(("gh", "rho", "zeta", "R", "sigma")) for p in rule["particles"]):
            continue
        families.add(tuple(sorted(_digest_multiplet(p) for p in rule["particles"])))
    ours = {c for c in canon(s0_classical.leg_classes()) if len(c) >= 3}
    assert ours == families


@pytest.fixture(scope="module")
def frule(sm, s0_classical):
    simp = sm.simplifier()

    def rule(names):
        legs = standard_legs(sm.fields, names)
        return [l.labels for l in legs], simp(I * s0_classical.vertex_function(legs, simp))

    return rule


def check(sm, value, reference, last_momentum=None):
    simp = sm.simplifier()
    if last_momentum is not None:
        others = sum((tn.p(k, _w) for k in range(1, last_momentum)), Expression.num(0))
        reference = reference.replace(tn.p(last_momentum, _w), -others)
    assert simp(value - reference) == 0


def test_gauge_self_couplings(sm, frule):
    p, g = tn.p, tn.g
    for r in ("G", "W"):
        c = sm.coupling(r)
        (m1, a1), (m2, a2), (m3, a3) = frule([r, r, r])[0]
        ref = c * tn.f(a1, a2, a3) * (g(m1, m2) * (p(1, m3) - p(2, m3)) + g(m2, m3) * (p(2, m1) - p(3, m1))
                                      + g(m3, m1) * (p(3, m2) - p(1, m2)))
        check(sm, frule([r, r, r])[1], ref, last_momentum=3)  # eq:V-vvv
        (m1, a), (m2, b), (m3, c_), (m4, e) = frule([r] * 4)[0]
        x = tn.fresh(sm.factors[r].adjoint_space())
        ff = (c**2 * tn.f(a, b, x) * tn.f(c_, e, x), c**2 * tn.f(a, c_, x) * tn.f(b, e, x),
              c**2 * tn.f(a, e, x) * tn.f(b, c_, x))
        ref = -I * (ff[0] * (g(m1, m3) * g(m2, m4) - g(m1, m4) * g(m2, m3))
                    + ff[1] * (g(m1, m2) * g(m3, m4) - g(m1, m4) * g(m2, m3))
                    + ff[2] * (g(m1, m2) * g(m3, m4) - g(m1, m3) * g(m2, m4)))
        check(sm, frule([r] * 4)[1], ref)  # eq:V-vvvv


def test_two_point_vertex_functions(sm, frule):
    """V_2(p) = pslash (p the incoming Psi momentum), Gamma_{H^dag H} = p^2 - m^2, Gamma_AA transverse."""
    (s1, *r1), (s2, *r2) = frule(["Qbar", "Q"])[0]
    deltas = tn.delta(r1[0], r2[0]) * tn.delta(r1[1], r2[1])
    check(sm, frule(["Qbar", "Q"])[1], I * tn.gamma(tn.lor(500), s1, s2) * tn.p(2, tn.lor(500)) * deltas, 2)
    (i,), (j,) = frule(["Hdag", "H"])[0]
    check(sm, frule(["Hdag", "H"])[1], I * (tn.sp(1, 1) - sm.param("m2")) * tn.delta(i, j))
    (m1, a1), (m2, a2) = frule(["G", "G"])[0]
    check(sm, frule(["G", "G"])[1],
          -I * (tn.g(m1, m2) * tn.sp(1, 1) - tn.p(1, m1) * tn.p(1, m2)) * tn.delta(a1, a2), 2)


@pytest.mark.parametrize("f, r", [(f, r) for f, rs in {"Q": "GWB", "U": "GB", "D": "GB", "L": "WB", "E": "B"}.items()
                                  for r in rs])
def test_fermion_gauge_vertices(sm, frule, f, r):
    """eq:V-ffv: i gbar^mu P_{c_f} Tg^a_f, the physical chirality only."""
    labs, value = frule([f + "bar", f, r])
    (s1, *g1), (s2, *g2), legA = labs
    proj = tn.PL if sm.fermions[f].chirality == "L" else tn.PR
    u = tn.fresh(tn.SPINOR)
    ref = I * tn.gammabar(legA[0], s1, u) * proj(u, s2) * sm.coupling(r)
    slots = sm.fields[f].slots[1:]
    for k, slot in enumerate(slots):
        if slot.factor == r:
            ref = ref * tn.T(legA[1], g1[k], g2[k])
        else:
            ref = ref * tn.delta(g1[k], g2[k])
    if sm.factors[r].abelian:
        q = sm.fields[f].charges["B"]
        ref = ref * Expression.num(q.numerator) / q.denominator
    check(sm, value, ref)


def test_higgs_vertices(sm, frule):
    g, gp, lam = sm.coupling("W"), sm.coupling("B"), sm.param("lambda")
    half = Expression.num(1) / 2
    # eq:V-ssv: H_j incoming p (leg 2), H^dag_i outgoing p' (leg 1 has incoming momentum -p')
    (i,), (j,), (mu, a) = frule(["Hdag", "H", "W"])[0]
    check(sm, frule(["Hdag", "H", "W"])[1], I * (tn.p(2, mu) - tn.p(1, mu)) * g * tn.T(a, i, j), 3)
    (i,), (j,), (mu,) = frule(["Hdag", "H", "B"])[0]
    check(sm, frule(["Hdag", "H", "B"])[1], I * (tn.p(2, mu) - tn.p(1, mu)) * gp * half * tn.delta(i, j), 3)
    # i g^{mu nu} {Tg^a, Tg^b}_{ij}
    (i,), (j,), (m1, a), (m2, b) = frule(["Hdag", "H", "W", "W"])[0]
    k = tn.fresh(sm.fields["H"].slots[0].space)
    anti = tn.T(a, i, k) * tn.T(b, k, j) + tn.T(b, i, k) * tn.T(a, k, j)
    check(sm, frule(["Hdag", "H", "W", "W"])[1], I * tn.g(m1, m2) * g**2 * anti)
    (i,), (j,), (m1, a), (m2,) = frule(["Hdag", "H", "W", "B"])[0]
    check(sm, frule(["Hdag", "H", "W", "B"])[1], I * tn.g(m1, m2) * g * gp * tn.T(a, i, j))
    (i,), (j,), (m1,), (m2,) = frule(["Hdag", "H", "B", "B"])[0]
    check(sm, frule(["Hdag", "H", "B", "B"])[1], I * tn.g(m1, m2) * gp**2 * half * tn.delta(i, j))
    # eq:V-ssss
    (i,), (j,), (k,), (l_,) = frule(["Hdag", "H", "Hdag", "H"])[0]
    check(sm, frule(["Hdag", "H", "Hdag", "H"])[1],
          -2 * I * lam * (tn.delta(i, j) * tn.delta(k, l_) + tn.delta(i, l_) * tn.delta(k, j)))


def test_yukawa_vertices(sm, frule):
    """eq:V-yuk (scalar leg labelled by the field operator that stands in the term)."""
    yu, yd, ye = (sm.param(n) for n in ("y_u", "y_d", "y_e"))
    (s1, a1, j), (s2, b2), (k,) = frule(["Qbar", "U", "Hdag"])[0]
    check(sm, frule(["Qbar", "U", "Hdag"])[1], -I * yu * tn.epsf(j, k) * tn.delta(a1, b2) * tn.PR(s1, s2))
    (s1, a1), (s2, b2, j), (k,) = frule(["Ubar", "Q", "H"])[0]
    check(sm, frule(["Ubar", "Q", "H"])[1], -I * yu * tn.epsf(j, k) * tn.delta(a1, b2) * tn.PL(s1, s2))
    (s1, a1, j), (s2, b2), (k,) = frule(["Qbar", "D", "H"])[0]
    check(sm, frule(["Qbar", "D", "H"])[1], -I * yd * tn.delta(j, k) * tn.delta(a1, b2) * tn.PR(s1, s2))
    (s1, a1), (s2, b2, j), (k,) = frule(["Dbar", "Q", "Hdag"])[0]
    check(sm, frule(["Dbar", "Q", "Hdag"])[1], -I * yd * tn.delta(j, k) * tn.delta(a1, b2) * tn.PL(s1, s2))
    (s1, j), (s2,), (k,) = frule(["Lbar", "E", "H"])[0]
    check(sm, frule(["Lbar", "E", "H"])[1], -I * ye * tn.delta(j, k) * tn.PR(s1, s2))
    (s1,), (s2, j), (k,) = frule(["Ebar", "L", "Hdag"])[0]
    check(sm, frule(["Ebar", "L", "Hdag"])[1], -I * ye * tn.delta(j, k) * tn.PL(s1, s2))


def test_hermiticity(sm, s0_classical, sectors):
    simp = sm.simplifier()
    assert LocalFunctional.integrate(s0_classical - s0_classical.conjugate()).vanishes(simp)
    for name, part in sectors.items():
        assert LocalFunctional.integrate(part - part.conjugate()).vanishes(simp), name
    assert sectors["yukawa"].conjugate() == sectors["yukawa"]  # h.c. terms generated: exactly hermitian
    anti = LocalFunctional.integrate(sectors["fermion"].scale(I))
    assert not is_zero(sm, LocalFunctional.integrate(anti - anti.conjugate()))  # certified non-zero


def test_rigid_invariance(sm, s0_classical, field):
    assert sm.rigid_defects(s0_classical, field) == []
    W = sm.fields["H"].slots[0].space
    i, j = tn.fresh(W), tn.fresh(W)
    eps_bilinear = Density.monomial(sm.fields, tn.epsf(i, j), Occurrence("Hdag", (i,)), Occurrence("H", (j,)))
    assert len(sm.rigid_defects(eps_bilinear, field)) == 1  # eps couples Fbar x F: not invariant
    s, c, w = sm.slot_labels("Q")
    t = tn.fresh(tn.SPINOR)
    charged = Density.monomial(sm.fields, tn.delta(s, t), Occurrence("Qbar", (s, c, w)), Occurrence("U", (t, c)),
                               Occurrence("H", (w,)))
    assert len(sm.rigid_defects(charged, field)) == 1  # hypercharges -1/6 + 2/3 + 1/2 != 0 (Qbar U Hdag is invariant)


def test_evanescent_part(sm, s0_classical, sectors):
    """Only the fermion kinetic term splits into a physical and an evanescent part among the fermion terms
    (report eq:Lferm-weyl); the gauge vertices, Yukawa terms and the potential are four-dimensional."""
    simp = sm.simplifier()
    ev = s0_classical.evanescent(simp)
    ferm_classes = {c for c in canon(ev.leg_classes()) if any(x.rstrip("bar") in FERMIONS for x in c if x != "B")}
    assert ferm_classes == {tuple(sorted((f, f + "bar"))) for f in FERMIONS}
    for name in ("potential", "yukawa"):
        assert sectors[name].evanescent(simp).is_zero(), name
    phys_fermion = sectors["fermion"].physical(simp)
    vertex_part = [m for m in sectors["fermion"].monomials() if len(m.occurrences) == 3]
    assert all(len(m.occurrences) in (2, 3) for m in phys_fermion.monomials())
    assert len(vertex_part) == 10
    # P + E reproduces the split action, and the 4D limit is the physical part at d = 4
    assert s0_classical.physical(simp) + ev == s0_classical.split_lorentz(simp)
