"""Gate V3: the classes of external legs (hepkit_renorm.classes) against the report.

Golden data: ``tests/data/report_classes.json`` (the 18 classes of report tab:classes with their multiplet-level
assignments from ``SCRATCHPAD/model/classes.py`` / ``enumeration.json``; see tests/data/README.md). The counts of
tab:Qbasis and the class lists of the plane-wave engine (``SCRATCHPAD/st/fingerprints.py``: 39 ghost-number-1
classes) and of report 04z (34 ghost-number-2 classes of the check ``b X_1 = 0``) are written out below.

Every difference between the derived classes and the golden assignments is asserted explicitly:

* ghost number 1, derived but not in ``classes.py``: ``omega omega rhotilde`` (seven gauge-factor assignments; the
  report lists ``gh gh tilde(srcA)`` among the five source functions of sec:g1_ccRF, the hand enumeration has no
  pattern for it) and the zero-momentum one-point function ``omega^B`` (``int omega^B`` is rigidly and Lorentz
  invariant; its only one-loop graph is a scaleless tadpole);
* in ``classes.py`` but not derived: ``omega_W W W G G``, ``omega_W W W W W``, ``omega_W W W B B`` and
  ``omega_W omega_W omega_W zeta_W``: every monomial vanishes identically (Bose / Grassmann symmetry against the
  SU(2) tensors), which the hand rules did not see; they have no diagrams (enumeration.json).
"""

from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction

import numpy as np
import pytest
from conftest import DATA

from hepkit_renorm import classes as C
from hepkit_renorm import linalg, model_from_dict
from hepkit_renorm.exact import default_field
from hepkit_renorm.groups import GaugeFactor, invariant_tensors

TOKEN = {"G": "G", "W": "W", "B": "B", "gh_G": "cG", "gh_W": "cW", "gh_B": "cB", "H": "H", "Hdag": "Hdag",
         "src_G": "rhoG", "src_W": "rhoW", "src_B": "rhoB", "src_gh_G": "zetaG", "src_gh_W": "zetaW",
         "src_gh_B": "zetaB", "src_H": "sigmadag", "src_Hdag": "sigma"}
for _f in "QUDLE":
    TOKEN.update({_f: _f, _f + "bar": _f + "bar", "src_" + _f: "R" + _f + "bar", "src_" + _f + "bar": "R" + _f})

# report tab:classes rows as unions of derived sectors
REPORT_SECTORS = {
    "g0_VV": ["g0_VV"], "g0_VVV": ["g0_VVV"], "g0_VVVV": ["g0_VVVV"], "g0_FF": ["g0_FF"], "g0_FFV": ["g0_FFV"],
    "g0_SS": ["g0_SS"], "g0_SSV_SSVV": ["g0_SSV", "g0_SSVV"], "g0_SSSS": ["g0_SSSS"], "g0_FFS": ["g0_FFS"],
    "g0_ghost_sources": ["g0_sources"], "g1_cV": ["g1_cV"], "g1_cVV": ["g1_cVV"],
    "g1_cVVV_cVVVV": ["g1_cVVV", "g1_cVVVV"], "g1_cFF": ["g1_cFF"], "g1_cFFV": ["g1_cFFV"],
    "g1_cSS_cSSV": ["g1_cSS", "g1_cSSV", "g1_cSSVV", "g1_cSSSS"], "g1_cFFS": ["g1_cFFS"], "g1_ccRF": ["g1_sources"]}

ONLY_DERIVED = {"cB", "cG cG rhoG", "cG cG rhoB", "cB cG rhoG", "cW cW rhoW", "cW cW rhoB", "cB cW rhoW",
                "cB cB rhoB"}
ONLY_HAND = {"G G W W cW", "W W W W cW", "B B W W cW", "cW cW cW zetaW"}


def _key(tokens) -> str:
    return " ".join(sorted(tokens))


@pytest.fixture(scope="module")
def golden() -> dict:
    with open(DATA / "report_classes.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def derived(sm) -> dict[int, list[C.LegClass]]:
    return {g: C.leg_classes(sm, g) for g in (0, 1, 2)}


def _golden_assignments(golden, g) -> dict[str, str]:
    return {_key(a): name for name, v in golden["classes"].items() if v["ghost_number"] == g
            for a in v["assignments"]}


# ---------------------------------------------------------------------------------------------- power counting
def test_power_counting_and_external_legs(sm):
    assert C.brst_dimension(sm) == 1
    assert [C.dimension_bound(sm, g) for g in (0, 1, 2)] == [4, 5, 6]
    ext = C.external_fields(sm)
    assert not any(n.startswith(("agh_", "b_")) for n in ext)
    assert len(ext) == len(sm.fields) - 6 and "src_G" in ext and "src_Hdag" in ext


# ---------------------------------------------------------------------------------------------- group traces
@pytest.mark.parametrize("group, reps", [
    ("SU(2)", ["F", "Fbar"]), ("SU(2)", ["F", "F"]), ("SU(2)", ["A", "A", "A"]), ("SU(2)", ["F", "F", "A"]),
    ("SU(2)", ["A", "A", "A", "A"]), ("SU(2)", ["F", "Fbar", "F", "Fbar"]), ("SU(3)", ["F", "Fbar", "A"]),
    ("SU(3)", ["F", "F", "F"]), ("SU(3)", ["A", "A", "A"]), ("SU(3)", ["F", "F", "Fbar", "Fbar"]),
    ("U(1)", [Fraction(1, 2), Fraction(-1, 3), Fraction(-1, 6)]), ("U(1)", [Fraction(1), Fraction(1)])])
def test_rigid_invariants_are_the_null_space_of_the_generators(group, reps):
    """sigma = 1: the character integral is the dimension of the exact null space of sum_k T^a_(k)."""
    factor = GaugeFactor("X", group, "gx")
    assert C.invariant_trace(factor, [(r, 1) for r in reps]) == len(invariant_tensors(factor, reps, default_field()))


@pytest.mark.parametrize("group, reps, perm", [
    ("SU(2)", ["A", "A", "A"], (1, 0, 2)), ("SU(2)", ["A", "A", "A"], (1, 2, 0)), ("SU(3)", ["A", "A", "A"], (1, 0, 2)),
    ("SU(3)", ["A", "A", "A"], (1, 2, 0)), ("SU(2)", ["A", "A", "A", "A"], (1, 0, 3, 2)),
    ("SU(2)", ["A", "A", "A", "A"], (1, 2, 3, 0)), ("SU(2)", ["F", "F", "A"], (1, 0, 2))])
def test_permutation_traces_on_the_null_space(group, reps, perm):
    """The trace of a permutation of identical legs on the explicit null space equals the character form."""
    field = default_field()
    factor = GaugeFactor("X", group, "gx")
    basis = invariant_tensors(factor, reps, field)
    b = linalg.as_field_matrix(np.stack([t.reshape(-1) for t in basis], axis=1), field)
    moved = linalg.as_field_matrix(np.stack([np.transpose(t, perm).reshape(-1) for t in basis], axis=1), field)
    trace = int(np.trace(linalg.solve(b, moved, field))) % field.p
    seen, cycles = set(), []
    for start in range(len(perm)):
        length, k = 0, start
        while k not in seen:
            seen.add(k)
            k, length = perm[k], length + 1
        if length:
            cycles.append((reps[start], length))
    assert trace == C.invariant_trace(factor, cycles) % field.p


# ---------------------------------------------------------------------------------------------- Lorentz x group x statistics
def test_structure_counts_of_report_completeness_paragraph(sm):
    """Report 03b, completeness of tab:Qbasis: momentum-space counting with Bose symmetry and epsilon."""
    sc = C.structure_counts
    assert sc(sm, ("G", "G"), 2) == (1, 0, 2)            # Q_M; Q_T, Q_L (epsilon p p = 0)
    assert sc(sm, ("G", "G", "G"), 1) == (0, 2)          # f and d: one each
    assert sc(sm, ("W", "W", "W"), 1) == (0, 1) and sc(sm, ("B", "B", "B"), 1) == (0, 1)
    assert sc(sm, ("W", "W", "B"), 1) == (0, 4)          # Q_B1, Q_B2, Q_B3 and Q_Beps
    assert sc(sm, ("G", "G", "W"), 1) == (0, 0) and sc(sm, ("G", "B", "B"), 1) == (0, 0)
    assert sc(sm, ("G", "G", "G", "G"), 0) == (3,) and sc(sm, ("G", "G", "G", "B"), 0) == (2,)  # d and f epsilon
    assert sc(sm, ("W", "W", "W", "B"), 0) == (1,) and sc(sm, ("W", "W", "W", "W"), 0) == (2,)
    assert sc(sm, ("Q", "Qbar"), 1) == (0, 1)             # only psibar gamma.p psi for one chirality
    assert sc(sm, ("Q", "Qbar"), 1, physical=False) == (2, 2)   # Dirac fields: {1, gamma5} x {1, pslash}
    assert sc(sm, ("H", "Hdag"), 2) == (1, 0, 1) and sc(sm, ("H", "H", "Hdag", "Hdag"), 0) == (1,)
    assert sc(sm, ("gh_G", "gh_G"), 4) == (0,) * 5      # omega^a P(p^2) omega^a: Grassmann antisymmetry
    assert sc(sm, ("gh_B", "gh_B", "src_gh_B"), 0) == (0,) and sc(sm, ("gh_G", "gh_G", "src_gh_G"), 0) == (1,)


# ---------------------------------------------------------------------------------------------- gate V3, ghost number 0
def test_ghost_number_0_classes_are_those_of_the_report(sm, derived, golden):
    mine = {_key(TOKEN[n] for n in c.fields): c for c in derived[0]}
    ref = _golden_assignments(golden, 0)
    assert set(mine) == set(ref) and len(mine) == 85
    for k, c in mine.items():
        assert c.key in REPORT_SECTORS[ref[k]]


def test_ghost_number_0_structure_counts_match_the_q_basis(derived):
    """Per sector, the number of independent monomials equals the rows of tab:Qbasis (97 in total), except for
    the matter-source monomials: the derived count has both hermitian combinations of each pair (src_f gh f,
    gh fbar src_fbar) and (src_H gh H, gh Hdag src_Hdag), the report's Q^{fr}_R, Q^{(r)}_sigma one (12 more)."""
    per = Counter()
    for c in derived[0]:
        per[c.key] += c.n_structures
    assert dict(per) == {"g0_VV": 9, "g0_VVV": 12, "g0_VVVV": 15, "g0_FF": 5, "g0_FFV": 10, "g0_SS": 2, "g0_SSV": 4,
                         "g0_SSVV": 4, "g0_SSSS": 1, "g0_FFS": 6, "g0_sources": 41}
    src = Counter()
    for c in derived[0]:
        if c.sector == "sources":
            src[c.signature] += c.n_structures
    # Q_rho_partial 3, Q_rho[t] 10, Q_zeta[t] 4; Q_R: 10 pairs, Q_sigma: 2 pairs
    assert dict(src) == {"c*V": 3, "cV*V": 10, "cc*c": 4, "cF*F": 20, "cS*S": 4}
    assert sum(per.values()) == 97 + 12


# ---------------------------------------------------------------------------------------------- gate V3, ghost number 1
def test_ghost_number_1_classes_are_those_of_the_report(sm, derived, golden):
    mine = {_key(TOKEN[n] for n in c.fields): c for c in derived[1]}
    ref = _golden_assignments(golden, 1)
    assert set(mine) - set(ref) == ONLY_DERIVED
    assert set(ref) - set(mine) == ONLY_HAND
    for k, c in mine.items():
        if k in ref:
            assert c.key in REPORT_SECTORS[ref[k]]
    inv = {v: k for k, v in TOKEN.items()}
    for k in ONLY_HAND:   # identically zero, for the physical components and for full Dirac fields alike
        legs = [inv[t] for t in k.split()]
        bound = C.dimension_bound(sm, 1) - sum(sm.fields[n].dimension for n in legs)
        assert not any(C.structure_counts(sm, legs, int(bound))), k
    assert {c.fields for c in derived[1] if c.signature == "c"} == {("gh_B",)}
    rho = {c.fields: c.structures for c in derived[1] if c.signature == "cc*V"}   # rhotilde^a t^{abc} omega^b d omega^c
    assert rho[("gh_G", "gh_G", "src_G")] == (0, 2) and rho[("gh_B", "gh_B", "src_B")] == (0, 1) and len(rho) == 7
    assert len(derived[1]) == 164 - len(ONLY_HAND) + len(ONLY_DERIVED)


def test_component_tuples_of_the_golden_assignments(derived, golden):
    """The 885 component tuples of enumeration.json belong to the 249 hand assignments; the four identically
    vanishing assignments account for the difference to the derived classes."""
    tuples = {_key(a.split()): n for v in golden["classes"].values() for a, n in v["component_tuples"].items()}
    assert sum(tuples.values()) == golden["total_component_tuples"] == 885
    mine = {_key(TOKEN[n] for n in c.fields) for g in (0, 1) for c in derived[g]}
    assert sum(n for k, n in tuples.items() if k in mine) == 885 - sum(tuples[k] for k in ONLY_HAND)


def test_superficial_degrees(derived):
    """tab:classes column dg(G) = D_g - sum dim."""
    degree = {"g0_VV": 2, "g0_VVV": 1, "g0_VVVV": 0, "g0_FF": 1, "g0_FFV": 0, "g0_SS": 2, "g0_SSV": 1, "g0_SSVV": 0,
              "g0_SSSS": 0, "g0_FFS": 0, "g1_cV": 3, "g1_cVV": 2, "g1_cVVV": 1, "g1_cVVVV": 0, "g1_cFF": 1,
              "g1_cFFV": 0, "g1_cSS": 2, "g1_cSSV": 1, "g1_cSSVV": 0, "g1_cSSSS": 0, "g1_cFFS": 0}
    for g in (0, 1):
        for c in derived[g]:
            if c.key in degree:
                assert c.degree == degree[c.key], c
            elif c.signature in ("c*V", "cF*F", "cS*S", "cV*V", "cc*c"):
                assert c.degree == (1 if c.signature == "c*V" else 0), c
            elif c.ghost_number == 1 and c.sector == "sources":
                assert c.degree == (1 if c.signature == "cc*V" else 0), c


# ---------------------------------------------------------------------------------------------- engine-level lists
def _combined(sm, c: C.LegClass) -> str:
    """The class with all gauge factors merged (the level of the plane-wave engine: w, A, rt, z)."""
    out = []
    for n in c.fields:
        fi = sm.fields[n]
        base = sm.fields[fi.source_of] if fi.is_source else fi
        if base.kind in ("gauge", "ghost"):
            out.append({("gauge", False): "A", ("ghost", False): "w", ("gauge", True): "rt",
                        ("ghost", True): "z"}[(base.kind, fi.is_source)])
        else:
            out.append(n)
    return " ".join(sorted(out))


ST_G1 = {"A w", "A A w", "A A A w", "A A A A w", "H Hdag w", "A H Hdag w", "A A H Hdag w", "H H Hdag Hdag w",
         "H src_H w w", "Hdag src_Hdag w w", "rt w w", "A rt w w", "w w w z"}
for _f in "QUDLE":
    ST_G1 |= {f"{_f} {_f}bar w", f"A {_f} {_f}bar w", f"{_f} src_{_f} w w", f"{_f}bar src_{_f}bar w w"}
ST_G1 |= {"H Q Ubar w", "Hdag Qbar U w", "D H Qbar w", "Dbar Hdag Q w", "E H Lbar w", "Ebar Hdag L w"}


def test_engine_level_ghost_number_1_classes(sm, derived):
    """st/fingerprints.py evaluates 39 classes; two of them (omega omega Rbar_E E and its conjugate) admit no
    non-zero monomial (omega^B omega^B = 0), and the derived list adds the one-point function omega^B."""
    assert len(ST_G1) == 39
    mine = {_combined(sm, c) for c in derived[1]}
    assert mine == (ST_G1 - {"E src_E w w", "Ebar src_Ebar w w"}) | {"w"}


def test_engine_level_ghost_number_2_classes(sm, derived):
    """Report 04z: b X_1 = 0 was checked on 34 classes (gh gh with gauge, Higgs, fermion and Yukawa legs, gh gh gh
    Rbar psi and conjugate). The derivation adds the five bosonic-source classes, which b maps into only from the
    bosonic-source components of X_1 (zero at one loop)."""
    report = {"A w w", "A A w w", "A A A w w", "A A A A w w", "H Hdag w w", "A H Hdag w w", "A A H Hdag w w",
              "H H Hdag Hdag w w", "H Q Ubar w w", "Hdag Qbar U w w", "D H Qbar w w", "Dbar Hdag Q w w",
              "E H Lbar w w", "Ebar Hdag L w w"}
    for f in "QUDLE":
        report |= {f"{f} {f}bar w w", f"A {f} {f}bar w w", f"{f} src_{f} w w w", f"{f}bar src_{f}bar w w w"}
    assert len(report) == 34
    mine = {_combined(sm, c) for c in derived[2]}
    assert mine == report | {"rt w w w", "A rt w w w", "H src_H w w w", "Hdag src_Hdag w w w", "w w w w z"}


# ---------------------------------------------------------------------------------------------- generic properties
def test_classes_are_closed_under_conjugation(sm, derived):
    for cls in derived.values():
        index = {c.fields: c for c in cls}
        for c in cls:
            partner = index[C.conjugate_fields(sm, c.fields)]
            assert partner.structures == c.structures and partner.signature == c.signature


def test_signatures_and_sectors(sm):
    assert C.kind_signature(sm, ["gh_W", "Qbar", "Q", "W"]) == "cFFV"
    assert C.kind_signature(sm, ["gh_G", "gh_G", "src_Q", "Q"]) == "ccF*F"
    assert C.sector(sm, ["src_H", "gh_W", "H"]) == "sources" and C.sector(sm, ["H", "Hdag", "B"]) == "SSV"
    with pytest.raises(ValueError):
        C.kind_signature(sm, ["agh_G"])


def test_classes_are_deterministic(sm, derived):
    again = C.leg_classes(sm, 1)
    assert [c.fields for c in again] == [c.fields for c in derived[1]]
    assert [c.structures for c in again] == [c.structures for c in derived[1]]


TOY = {"name": "toy", "gauge": {"X": {"group": "SU(2)", "coupling": "gx"}, "Y": {"group": "U(1)", "coupling": "gy"}},
       "fermions": {"P": {"chirality": "L", "X": "F", "Y": "1/3"}, "N": {"chirality": "R", "Y": "1/3"}},
       "scalars": {"phi": {"complex": False, "X": "A"}},
       "parameters": {"mphi": {}},
       "potential": [{"coupling": "mphi", "bilinears": [["phi", "phi"]]}],
       "gauge_fixing": {"kind": "linear_covariant", "xi": "1"}}


def test_a_model_with_a_real_adjoint_scalar():
    """Nothing is Standard-Model specific: SU(2) x U(1), a chiral doublet, a singlet, a real adjoint scalar."""
    m = model_from_dict(TOY)
    cls = {c.fields: c for c in C.leg_classes(m, 0)}
    assert cls[("phi", "phi")].structures == (1, 0, 1)
    assert ("phi", "phi", "phi") not in cls                       # eps^{abc}: no Bose-symmetric polynomial
    assert cls[("phi", "phi", "phi", "phi")].structures == (1,)
    assert ("P", "Pbar", "phi") not in cls                        # same-chirality bilinear is a vector
    assert cls[("X", "phi", "phi")].structures == (0, 1)           # eps^{abc} X^a phi^b d phi^c
    assert ("N", "Pbar", "phi") not in cls and ("N", "Nbar") in cls
    assert {c.key for c in C.leg_classes(m, 1)} >= {"g1_c", "g1_cV", "g1_cFF", "g1_sources"}
