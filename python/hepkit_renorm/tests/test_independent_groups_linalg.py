"""Independent review: group invariants and anomaly traces (report 03_setup tab:invariants, eq:SLR, eq:Dabc), and
the exact linear algebra of :mod:`hepkit_renorm.linalg` / :mod:`hepkit_renorm.exact` on adversarial input.

The expected group-theory numbers are typed from the report and recomputed here from independent generators
(``tests/independent_support.py``); the anomaly traces of the lepton-free content are worked out by hand in the
docstrings.
"""

from __future__ import annotations

import copy
import itertools
import random
import tomllib
from fractions import Fraction
from math import isqrt
from pathlib import Path

import numpy as np
import pytest
from independent_support import Evaluator, P, Sampler, fe, structure_constants, su_generators
from symbolica import Expression

import hepkit_renorm
from hepkit_renorm import linalg
from hepkit_renorm.exact import PrimeField, default_field, second_field
from hepkit_renorm.spec import model_from_dict

SPEC = Path(hepkit_renorm.__file__).parent / "models" / "sm_one_family.toml"

# tab:fields and tab:invariants (S_f^(r) at unit couplings; C_f as (G, W, B) coefficients of g_s^2, g^2, g'^2)
Y = {"Q": Fraction(1, 6), "U": Fraction(2, 3), "D": Fraction(-1, 3), "L": Fraction(-1, 2), "E": Fraction(-1),
     "H": Fraction(1, 2)}
REPS = {"Q": ("G", "W"), "U": ("G",), "D": ("G",), "L": ("W",), "E": (), "H": ("W",)}
LEFT = {"Q", "L"}
S_TABLE = {("Q", "G"): Fraction(1), ("U", "G"): Fraction(1, 2), ("D", "G"): Fraction(1, 2),
           ("Q", "W"): Fraction(3, 2), ("L", "W"): Fraction(1, 2), ("H", "W"): Fraction(1, 2),
           ("Q", "B"): Fraction(1, 6), ("U", "B"): Fraction(4, 3), ("D", "B"): Fraction(1, 3),
           ("L", "B"): Fraction(1, 2), ("E", "B"): Fraction(1), ("H", "B"): Fraction(1, 2)}
C_TABLE = {"Q": (Fraction(4, 3), Fraction(3, 4), Fraction(1, 36)), "U": (Fraction(4, 3), 0, Fraction(4, 9)),
           "D": (Fraction(4, 3), 0, Fraction(1, 9)), "L": (0, Fraction(3, 4), Fraction(1, 4)), "E": (0, 0, Fraction(1)),
           "H": (0, Fraction(3, 4), Fraction(1, 4))}
C_AD = {"G": Fraction(3), "W": Fraction(2), "B": Fraction(0)}


def multiplet_generators(m: str) -> np.ndarray:
    """The 12 unit-coupling generators (8 gluons, 3 W, 1 B) on the colour x SU(2) components of a multiplet."""
    tg, tw = su_generators(3), su_generators(2)
    nc, nw = (3 if "G" in REPS[m] else 1), (2 if "W" in REPS[m] else 1)
    one_c, one_w = np.eye(nc, dtype=object), np.eye(nw, dtype=object)
    gens = []
    for a in range(8):
        gens.append(np.kron(tg[a], one_w) % P if nc == 3 else np.zeros((nc * nw,) * 2, dtype=object))
    for i in range(3):
        gens.append(np.kron(one_c, tw[i]) % P if nw == 2 else np.zeros((nc * nw,) * 2, dtype=object))
    gens.append(np.eye(nc * nw, dtype=object) * fe(Y[m]) % P)
    return np.array(gens, dtype=object)


BLOCKS = {"G": range(8), "W": range(8, 11), "B": range(11, 12)}


def anomaly(content) -> np.ndarray:
    total = np.zeros((12, 12, 12), dtype=object)
    for m in content:
        t = multiplet_generators(m)
        tt = np.einsum("bij,cjk->bcik", t, t) % P
        anti = (tt + np.transpose(tt, (1, 0, 2, 3))) % P
        d = np.einsum("aki,bcik->abc", t, anti) % P
        total = (total + (d if m in LEFT else -d)) % P
    return total


def test_report_invariants_from_independent_generators():
    """tab:invariants, the adjoint Casimirs and eq:SLR, recomputed from explicit generators."""
    for m in Y:
        t = multiplet_generators(m)
        for r, blk in BLOCKS.items():
            for a, b in itertools.product(blk, blk):
                tr = int(np.trace(t[a].dot(t[b]) % P)) % P
                assert tr == (fe(S_TABLE.get((m, r), Fraction(0))) if a == b else 0), (m, r, a, b)
        cas = sum(t[a].dot(t[a]) for a in range(12)) % P
        expected = sum(c for c in C_TABLE[m])  # unit couplings
        assert (cas == np.eye(t.shape[1], dtype=object) * fe(expected) % P).all(), m
    for r, n in (("G", 3), ("W", 2)):
        f, _ = structure_constants(su_generators(n))
        c = np.einsum("ace,bce->ab", f, f) % P
        assert (c == np.eye(n * n - 1, dtype=object) * fe(C_AD[r]) % P).all()
    lh = {r: sum(S_TABLE.get((m, r), 0) for m in ("Q", "L")) for r in "GWB"}
    rh = {r: sum(S_TABLE.get((m, r), 0) for m in ("U", "D", "E")) for r in "GWB"}
    assert (lh, rh) == ({"G": 1, "W": 2, "B": Fraction(2, 3)}, {"G": 1, "W": 0, "B": Fraction(8, 3)})


@pytest.fixture(scope="module")
def sm():
    return hepkit_renorm.builtin_model("sm_one_family")


def _eval_invariant(expr, couplings):
    ev = Evaluator(4, couplings, {"G": 3, "W": 2})
    return ev._scalar(expr.expand(), {})


def test_package_invariants_against_report(sm):
    """Model.group_invariants against tab:invariants at random couplings (as polynomials in g_s, g, g')."""
    smp = Sampler(61)
    c = {"g_s": smp.element(), "g": smp.element(), "gp": smp.element()}
    sq = {"G": c["g_s"] ** 2 % P, "W": c["g"] ** 2 % P, "B": c["gp"] ** 2 % P}
    inv = sm.group_invariants(default_field())
    for (m, r), val in inv["S"].items():
        assert _eval_invariant(val, c) == fe(S_TABLE.get((m, r), Fraction(0))) * sq[r] % P, (m, r)
    assert {m for m, _ in inv["S"]} == set(Y)
    for m, val in inv["C"].items():
        exp = sum(fe(x) * sq[r] for x, r in zip(C_TABLE[m], "GWB")) % P
        assert _eval_invariant(val, c) == exp, m
    for r, val in inv["C_ad"].items():
        assert _eval_invariant(val, c) == fe(C_AD[r]) * sq[r] % P, r


def _without_leptons():
    with open(SPEC, "rb") as fh:
        doc = tomllib.load(fh)
    doc = copy.deepcopy(doc)
    doc["name"] = "sm_no_leptons"
    for m in ("L", "E"):
        del doc["fermions"][m]
    doc["yukawa"] = [y for y in doc["yukawa"] if y["coupling"] != "y_e"]
    del doc["parameters"]["y_e"]
    return model_from_dict(doc)


def test_anomaly_traces(sm):
    """eq:Dabc vanishes for the content of tab:fields (independently and in the package). Without leptons, by hand:
    D^{BBB} = 2 [6 (1/6)^3] - 2 [3 (2/3)^3 + 3 (-1/3)^3] = 1/18 - 14/9 = -3/2, D^{IJB} = 2 (1/6) (1/2) 3 delta = 1/2
    delta (and permutations), D^{ABC}, D^{ABB'} and the SU(2)^3 part vanish; sum_abc (D^{abc})^2 = 9/4 + 9/4 = 9/2."""
    assert not anomaly("QUDLE").any()
    d = anomaly("QUD")
    assert d[11, 11, 11] == fe(Fraction(-3, 2))
    for i, j in itertools.product(range(8, 11), repeat=2):
        for idx in ((i, j, 11), (i, 11, j), (11, i, j)):
            assert d[idx] == (fe(Fraction(1, 2)) if i == j else 0)
    assert not d[:8, :8, :].any() and not d[8:11, 8:11, 8:11].any()
    assert int(np.sum(d * d) % P) == fe(Fraction(9, 2))
    field = default_field()
    pkg = sm.anomaly_tensor(field)
    assert not np.any(np.asarray(pkg, dtype=object) % field.p)
    pkg2 = np.asarray(_without_leptons().anomaly_tensor(field), dtype=object)
    assert field.reconstruct(int(np.sum(pkg2 * pkg2) % field.p)) == Fraction(9, 2)
    assert field.reconstruct(int(pkg2[11, 11, 11])) == Fraction(-3, 2)
    assert field.reconstruct(int(pkg2[8, 8, 11])) == Fraction(1, 2)


# ---------------------------------------------------------------------------------------------- exact linear algebra
F = default_field()


def _rand(rng, shape, p=F.p):
    return rng.integers(0, p, size=shape, dtype=np.int64)


def test_rank_deficient_systems():
    rng = np.random.default_rng(71)
    u, v = _rand(rng, (20, 3)), _rand(rng, (3, 15))
    a = linalg.as_field_matrix(np.asarray(u, dtype=object).dot(np.asarray(v, dtype=object)), F)
    assert linalg.rank(a, F) == 3
    ker = linalg.nullspace(a, F)
    assert ker.shape == (12, 15)
    assert not (np.asarray(a, dtype=object).dot(np.asarray(ker.T, dtype=object)) % F.p).any()
    assert linalg.rank(ker, F) == 12
    x0 = _rand(rng, (15,))
    b = np.asarray(a, dtype=object).dot(np.asarray(x0, dtype=object)) % F.p
    x = linalg.solve(a, b, F)
    assert ((np.asarray(a, dtype=object).dot(np.asarray(x, dtype=object)) - b) % F.p == 0).all()
    with pytest.raises(ValueError, match="not unique"):
        linalg.solve(a, b, F, require_unique=True)
    # duplicated columns, zero rows and a zero matrix
    dup = np.concatenate([a, a[:, :2]], axis=1)
    assert linalg.rank(dup, F) == 3 and linalg.independent_columns(dup, F) == linalg.independent_columns(a, F)
    assert linalg.rank(np.zeros((4, 5), dtype=np.int64), F) == 0
    assert linalg.nullspace(np.zeros((0, 4), dtype=np.int64), F).shape == (4, 4)
    # entries far beyond 2^63 are reduced exactly (Python ints)
    big = np.array([[2**80 + 1, 2**70], [3 * 2**80 + 3, 3 * 2**70]], dtype=object)
    assert linalg.rank(big, F) == 1


def test_system_without_solution():
    rng = np.random.default_rng(72)
    a = linalg.as_field_matrix(np.asarray(_rand(rng, (6, 2)), dtype=object).dot(
        np.asarray(_rand(rng, (2, 4)), dtype=object)), F)
    good = np.asarray(a, dtype=object).dot(np.asarray(_rand(rng, (4,)), dtype=object)) % F.p
    bad = (good + np.eye(6, dtype=object)[0]) % F.p                  # leaves the 2-dimensional column space
    with pytest.raises(linalg.InconsistentSystemError) as err:
        linalg.solve(a, bad, F)
    assert err.value.columns == [0]
    both = np.stack([good, bad, good], axis=1)
    with pytest.raises(linalg.InconsistentSystemError) as err:
        linalg.solve(a, both, F)
    assert err.value.columns == [1]
    assert linalg.column_space_contains(a, good, F) and not linalg.column_space_contains(a, bad, F)


def test_large_denominators_need_several_primes():
    """x = A^-1 b with numerators and denominators of about 10^25: reconstructed exactly by CRT over several primes,
    and a single 31-bit prime is not enough."""
    rng = random.Random(73)
    n = 4
    a = [[rng.randint(-50, 50) for _ in range(n)] for _ in range(n)]
    x = [Fraction(rng.randint(10**24, 10**25) * (-1) ** k, rng.randint(10**24, 10**25)) for k in range(n)]
    b = [sum(Fraction(a[i][j]) * x[j] for j in range(n)) for i in range(n)]

    def compute(field: PrimeField):
        am = np.array([[field(v) for v in row] for row in a], dtype=object)
        bm = np.array([field(v) for v in b], dtype=object)
        return linalg.solve(am, bm, field, require_unique=True)

    got = linalg.reconstruct_rationals(compute)
    assert list(got) == x
    with pytest.raises(linalg.ReconstructionError):
        linalg.reconstruct_rationals(compute, max_primes=2)


def test_rational_reconstruction_edges():
    m = F.p * second_field().p
    for q in (Fraction(-7, 3), Fraction(0), Fraction(1, 10**9 + 7), Fraction(-(10**9), 3)):
        x = q.numerator * pow(q.denominator, -1, m) % m
        assert linalg.rational_reconstruction(x, m) == q
    bound = isqrt(F.p // 2)
    for q in (Fraction(-7, 3), Fraction(12345, 678), Fraction(-bound, bound - 1), Fraction(bound - 1, bound)):
        assert F.reconstruct(F(q)) == q
        assert linalg.rational_reconstruction(F(q), F.p) == q
    # a single prime can only certify rationals with |n|, d <= sqrt(p/2): a larger one comes back as some other
    # small rational, which is why reconstruct_rationals asks for a confirming prime (see the previous test)
    big = Fraction(2**30 + 1, 2**30 + 3)
    try:
        assert linalg.rational_reconstruction(F(big), F.p) != big
    except linalg.ReconstructionError:
        pass


def test_interpolation_degree_bound_is_checked():
    def cubic(field, pt):
        x, y = pt
        return (3 * x**3 + 5 * x * y + 7) % field.p
    with pytest.raises(linalg.InconsistentSystemError):
        linalg.interpolate(cubic, 2, 2, F)
    def values(field, pt):
        x, y = pt
        return np.array([(field(Fraction(1, 3)) * x**3 + 5 * x * y + 7) % field.p], dtype=object)
    polys = linalg.reconstruct_polynomials(values, 2, 3)
    assert polys == [{(3, 0): Fraction(1, 3), (1, 1): Fraction(5), (0, 0): Fraction(7)}]


def test_independent_evaluator_reads_numbers():
    """The independent evaluator reads Symbolica's complex rationals exactly."""
    from independent_support import IU, number
    for e, val in ((Expression.num(3) / 2 + Expression.I * 5 / 7, fe(Fraction(3, 2)) + fe(Fraction(5, 7)) * IU),
                   (Expression.I * 2 / 3, fe(Fraction(2, 3)) * IU), (-Expression.I, -IU),
                   (Expression.num(-7) / 3 - Expression.I * 2 / 9, fe(Fraction(-7, 3)) - fe(Fraction(2, 9)) * IU)):
        assert number(e.to_canonical_string()) == val % P
