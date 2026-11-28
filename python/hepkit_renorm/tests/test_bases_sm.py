"""Gate V4 for the one-family Standard Model: generated bases Q and N, the action of b, its kernel and image,
against report 03b (tab:Qbasis, tab:Nbasis, eq:closed, eq:btable-gauge, eq:btable-higgs) and the plane-wave engine's
b-table ``tests/data/st_btable_unit.json`` (copied from ``SCRATCHPAD/st/data/st_btable_unit.json``, see
``tests/data/README.md``).

Slow (about 1.5 min; the longest parts are the derived kernel description and b of the 353 monomials of N on the
280 ghost-number-2 classes): run with ``pytest --runslow tests/test_bases_sm.py``.

The report's named monomials are built in ``report_bases.py`` (test data only). Findings asserted here:

* Q: the exact count of hermitian real monomials is 109 = 97 + 12. The 97 of tab:Qbasis span a subspace; the 12
  others are ``i (M - M^dagger)`` for the matter-source pairs, i.e. ``i`` times the engine's anti-hermitian
  ``R_<f><r>-`` and ``sig_<r>-``: hermitian with real coefficients, so the hermitian basis contains them. Their
  b-images are independent of the others (rank 84 = 72 + 12), so they do not change the kernel (25).
* N: complete in every class (353); in the classes where the report claims completeness (cV; cVV at dimension 5)
  the counts agree; elsewhere tab:Nbasis lists the parts that b Q and X_1 need.
"""

import json
from fractions import Fraction

import numpy as np
import pytest
from conftest import DATA
from report_bases import Report

from hepkit_renorm import LocalFunctional, linalg
from hepkit_renorm.baction import b_closure, b_table, invariant_directions
from hepkit_renorm.bases import breaking_basis, counterterm_basis
from hepkit_renorm.brst import BRST
from hepkit_renorm.classes import sector
from hepkit_renorm.fingerprint import Fingerprinter, is_zero

pytestmark = pytest.mark.slow

UNIT = {"g_s": 1, "g": 1, "gp": 1, "m2": 1, "lambda": 1, "y_u": 1, "y_d": 1, "y_e": 1, "xi": 1}


@pytest.fixture(scope="module")
def brst(sm):
    return BRST(sm)


@pytest.fixture(scope="module")
def Q(sm):
    return counterterm_basis(sm)


@pytest.fixture(scope="module")
def N(sm):
    return breaking_basis(sm)


@pytest.fixture(scope="module")
def table(Q, N, brst):
    return b_table(Q, N, brst=brst)


@pytest.fixture(scope="module")
def report(sm):
    return Report(sm)


@pytest.fixture(scope="module")
def named(report):
    return report.q_basis()


@pytest.fixture(scope="module")
def named_images(named, brst):
    return {k: LocalFunctional.integrate(brst.b(f)) for k, f in named.items()}


def by_sector(model, basis, min_degree=0):
    out: dict[str, int] = {}
    for e in basis:
        if e.degree >= min_degree:
            key = sector(model, e.leg_class)
            out[key] = out.get(key, 0) + 1
    return out


# ---------------------------------------------------------------------------------------------- Q
def test_q_counts(sm, Q):
    """tab:Qbasis rows per sector; 29 + 12 source monomials (module docstring)."""
    assert len(Q) == 109
    assert by_sector(sm, Q) == {"VV": 9, "VVV": 12, "VVVV": 15, "FF": 5, "FFV": 10, "SS": 2, "SSV": 4, "SSVV": 4,
                                "SSSS": 1, "FFS": 6, "sources": 41}


def test_q_is_hermitian(sm, Q):
    for e in Q:
        assert is_zero(sm, LocalFunctional.integrate(e.functional - e.functional.conjugate())), e.name


def test_named_q_span(sm, Q, named):
    """The 109 engine monomials span Q; the 97 of tab:Qbasis are independent and hermitian; the engine's 12
    relative-sign variants are anti-hermitian, i times them hermitian and in Q."""
    fp = Fingerprinter(sm)
    names = list(named)
    coords = Q.coordinates([named[n] for n in names], fp)
    assert linalg.rank(coords, fp.field) == 109
    report97 = Report.report_names(names)
    assert len(report97) == 97
    assert linalg.rank(Q.coordinates([named[n] for n in report97], fp), fp.field) == 97
    anti = [n for n in names if n not in report97]
    for n in names:
        herm = is_zero(sm, LocalFunctional.integrate(named[n] - named[n].conjugate()))
        assert herm == (n in report97), n
    from hepkit_renorm import tensors as tn
    for n in anti:
        f = LocalFunctional.integrate(named[n].scale(tn.I))
        assert is_zero(sm, LocalFunctional.integrate(f - f.conjugate())), n


# ---------------------------------------------------------------------------------------------- N
def test_n_counts(sm, N):
    """Complete N: 353; the report's complete classes agree (cV: 6; cVV at dimension 5: 44), the others contain
    tab:Nbasis (counts there: cVVV 34, cVVVV 5, cSSV 17, cSSVV 9, sources 48)."""
    assert len(N) == 353
    full = by_sector(sm, N)
    assert full == {"c": 1, "cV": 6, "cVV": 50, "sources": 59, "cFF": 20, "cSS": 8, "cVVV": 84, "cFFV": 35,
                    "cSSV": 18, "cFFS": 16, "cVVVV": 41, "cSSVV": 13, "cSSSS": 2}
    assert by_sector(sm, N, min_degree=2)["cVV"] == 44
    report = {"cVVV": 34, "cVVVV": 5, "cSSV": 17, "cSSVV": 9, "sources": 48, "cFF": 20, "cFFV": 35, "cFFS": 16,
              "cSS": 8, "cSSSS": 2}
    assert all(full[k] >= v for k, v in report.items())


def test_named_n_in_span(sm, N, report):
    fp = Fingerprinter(sm)
    named_n = report.n_basis()
    coords = N.coordinates(list(named_n.values()), fp)
    assert linalg.rank(coords, fp.field) == len(named_n)


# ---------------------------------------------------------------------------------------------- b, kernel, image
def test_rank_and_kernel(sm, Q, table):
    assert len(Q) == 109 and table.rank == 84 and table.kernel.shape[0] == 25


def test_report_rank_and_kernel(sm, named, named_images):
    """On the 97 monomials of tab:Qbasis: rank 72, kernel 25 (report sec:Qbasis, sec:closed)."""
    names = list(named)
    report97 = Report.report_names(names)
    fp = Fingerprinter(sm)
    from hepkit_renorm.classes import leg_classes
    rows = [c.fields for c in leg_classes(sm, 1)]
    values = fp.matrix([named_images[n] for n in names], rows, 8).values
    assert linalg.rank(values, fp.field) == 84
    idx = [names.index(n) for n in report97]
    assert linalg.rank(values[:, idx], fp.field) == 72


def test_closed_directions_span_the_kernel(sm, Q, table, named, report):
    """eq:closed: the 25 operators are closed, independent, lie in the span of the 97 and span the kernel."""
    closed = report.closed(named)
    assert len(closed) == 25
    fp = table.fingerprinter
    coords = Q.coordinates(list(closed.values()), fp)
    for name, c in zip(closed, coords):
        assert table.is_closed(c), name
    assert linalg.rank(coords, fp.field) == 25
    assert linalg.rank(np.concatenate([coords, table.kernel]), fp.field) == 25
    report97 = Report.report_names(list(named))
    span97 = Q.coordinates([named[n] for n in report97], fp)
    assert linalg.rank(np.concatenate([span97, coords]), fp.field) == 97


def test_derived_description_of_the_kernel(sm, table, brst):
    """The kernel is spanned by derived generators: 15 b-exact directions (field renormalisations, from 18
    ghost-number -1 monomials), 8 coupling derivatives (three gauge couplings, m2, lambda, three Yukawas; the
    Yukawa phases are already b-exact combined with fermion phase rotations), 2 charge deformations
    (baryon and lepton number gauged with B; hypercharge is the coupling direction of g')."""
    k = invariant_directions(table, brst)
    assert k.dimension == 25 and k.remainder.shape[0] == 0
    assert k.rank_by_kind == {"exact": 15, "coupling": 23, "phase": 23, "charge": 25, "invariant": 25}
    exact = [d for d in k.directions if d.kind == "exact"]
    assert len(exact) == 18 and all(d.closed for d in exact)
    phases = {d.name: d.closed for d in k.directions if d.kind == "phase"}
    assert phases == {"phase g_s": False, "phase g": False, "phase gp": False, "phase y_u": True,
                      "phase y_d": True, "phase y_e": True}


def _st_table(field, names):
    doc = json.loads((DATA / "st_btable_unit.json").read_text())

    def value(x: str) -> int:
        if x.endswith("i"):
            return field(Fraction(x[:-1])) * field.i % field.p
        return field(Fraction(x))

    rows = [[value(d["B"].get(q, {}).get(n, "0")) for q in names] for d in doc.values() for n in d["N"]]
    return np.array(rows, dtype=np.int64)


def test_b_table_equals_the_engine_table_up_to_a_change_of_basis(sm, named_images, named):
    """st_btable_unit.json (unit couplings): the generated b on the 109 engine monomials has the same rank (84) and
    the same kernel, i.e. B_engine = M B_generated with M injective on the image."""
    names = list(named)
    fp = Fingerprinter(sm, values=UNIT)
    from hepkit_renorm.classes import leg_classes
    rows = [c.fields for c in leg_classes(sm, 1)]
    mine = fp.matrix([named_images[n] for n in names], rows, 8).values
    st = _st_table(fp.field, names)
    F = fp.field
    assert linalg.rank(mine, F) == linalg.rank(st, F) == 84
    k_mine, k_st = linalg.nullspace(mine, F), linalg.nullspace(st, F)
    assert k_mine.shape[0] == k_st.shape[0] == 25
    assert linalg.rank(np.concatenate([k_mine, k_st]), F) == 25


@pytest.mark.parametrize("lhs,rhs", [
    ("M_G", [(-1, "N1_G")]), ("M_W", [(-1, "N1_W")]), ("M_B", [(-1, "N1_B")]),
    ("L_G", [(-1, "N3_G"), (-1, "Nf2_G")]), ("L_B", [(-1, "N3_B")]),
    ("T_G", [(-1, "Nf1_G"), (1, "Nf2_G")]), ("T_W", [(-1, "Nf1_W"), (1, "Nf2_W")]),
    ("3f_G", [(1, "Nf1_G"), (-1, "Nf2_G"), (1, "Nff_G")]), ("3f_W", [(1, "Nf1_W"), (-1, "Nf2_W"), (1, "Nff_W")]),
    ("G4a", [(4, "NG4a")]), ("G4b", [(4, "NG4b")]),
    ("G4c", [("8/3 gs2", "Nff_G"), (Fraction(-4, 3), "NG4a"), (Fraction(8, 3), "NG4b")]),
    ("W4a", [(4, "NW4a")]), ("W4b", [("-4 g2", "Nff_W"), (4, "NW4a")]), ("B4", [(4, "NB4")]),
    ("GWa", [(2, "NGWa"), (2, "NGWb")]), ("GWb", [(2, "NGWc"), (2, "NGWd")]),
    ("WBa", [(2, "NWBa"), (2, "NWBb")]), ("GBb", [(2, "NGBc"), (2, "NGBd")]),
    ("SS", [(1, "NS_W"), (1, "NS_B")]), ("J_W", [(1, "NS_W"), (1, "NSV_WW"), (1, "NSV_BW")]),
    ("J_B", [(1, "NS_B"), (1, "NSV_WB"), (1, "NSV_BB")]), ("D_W", [(1, "NS2_W"), (1, "NSC")]),
    ("D_B", [(1, "NS2_B")]), ("HWW", [(1, "NSV_WW")]), ("HBB", [(1, "NSV_BB")]),
    ("HWB", [(1, "NSV_WB"), (1, "NSV_BW")]), ("HGG", [(2, "NSG")]),
    ("m", []), ("lam", []), ("yu+", []), ("yd-", []), ("ye+", []),
])
def test_named_b_table_rows(sm, report, named, named_images, lhs, rhs):
    """eq:btable-gauge and eq:btable-higgs: b Q_alpha = sum B N_beta on the report's named monomials (generic
    couplings; N^(rs)_{a..d} with the factor order of the report)."""
    from symbolica import Expression

    from hepkit_renorm import Density
    named_n = _named_n(report)
    total = Density(sm.fields) + named_images[lhs]
    for c, n in rhs:
        if isinstance(c, str):
            num, sym = c.split()
            factor = sm.coupling({"gs2": "G", "g2": "W"}[sym]) ** 2
            coeff = Expression.num(Fraction(num).numerator) / Fraction(num).denominator * factor
        else:
            coeff = Expression.num(Fraction(c).numerator) / Fraction(c).denominator
        total = total - named_n[n].scale(coeff)
    assert is_zero(sm, LocalFunctional.integrate(total), dimensions=(4,)), lhs


_N_CACHE: dict = {}


def _named_n(report):
    if "n" not in _N_CACHE:
        _N_CACHE["n"] = report.n_basis()
    return _N_CACHE["n"]


def test_b_squared_on_the_standard_model(sm, N, table, brst):
    """b X = 0 as fingerprint rows of b N on the 280 ghost-number-2 classes: every row of the b-table is
    annihilated (b^2 = 0 in four dimensions). Two samples per class (the longest test of the file)."""
    _, rows = b_closure(N, brst=brst, fingerprinter=table.fingerprinter, samples_per_class=2)
    p = table.field.p
    assert len(rows.layout) == 280
    assert not np.any((rows.values.astype(object) @ table.matrix.T.astype(object)) % p)
