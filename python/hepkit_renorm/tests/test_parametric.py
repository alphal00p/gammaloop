"""Parameter gradings, parametric fingerprints and graded reconstruction (hepkit_renorm.parametric)."""

from fractions import Fraction

import numpy as np
import pytest
from symbolica import Expression

from hepkit_renorm import LocalFunctional, linalg
from hepkit_renorm.bases import counterterm_basis
from hepkit_renorm.brst import BRST
from hepkit_renorm.classes import leg_classes
from hepkit_renorm.exact import default_field
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.parametric import GradingError, ParameterSpace, ParametricMatrix, parametric_rows, reconstruct


@pytest.fixture(scope="module")
def brst(sm):
    return BRST(sm)


@pytest.fixture(scope="module")
def space(sm, brst):
    return ParameterSpace(sm, brst.action, {"xi": 1})


def test_gradings_are_derived_from_the_action(space):
    """Loop weight n - 2 and mass dimension of every S_0 monomial fix the gradings of the parameters."""
    w = dict(zip(space.names, space.weights))
    d = dict(zip(space.names, space.dimensions))
    assert space.names == ("g_s", "g", "gp", "m2", "lambda", "y_u", "y_d", "y_e")
    assert w == {"g_s": 1, "g": 1, "gp": 1, "m2": 0, "lambda": 2, "y_u": 1, "y_d": 1, "y_e": 1}
    assert d == {"g_s": 0, "g": 0, "gp": 0, "m2": 2, "lambda": 0, "y_u": 0, "y_d": 0, "y_e": 0}
    # supports: weight 2 / dimension 0 are the 21 quadratic monomials of the six weight-1 couplings plus lambda
    assert len(space.support(2, 0)) == 22 and len(space.support(2, 2)) == 22 and len(space.support(4, 0)) == 148
    assert len(space.support(1, 2)) == 6 and space.support(0, 1) == []
    assert ParameterSpace(space.model, BRST(space.model).action, {n: 1 for n in space.names + ("xi",)}
                          ).support(4, 0) == [()]


def test_grade_of_functionals(sm, brst, space):
    """Generated basis monomials carry no parameters (grade 2 - n); dS_0/dg has grade -1."""
    Q = counterterm_basis(sm)
    for e in Q.elements[:20]:
        grade, _ = space.grade(e.functional)
        n = len(e.leg_class)
        assert grade == 2 - n
    deriv = LocalFunctional.integrate(brst.action.map_coefficients(
        lambda c: c.derivative(sm.param("g")).expand()))
    assert space.grade(deriv) == (-1, 0)
    with pytest.raises(GradingError):
        space.grade(LocalFunctional.integrate(Q.elements[0].functional + Q.elements[0].functional.scale(
            sm.param("g"))))


def test_parametric_rows_equal_direct_fingerprints(sm, brst, space):
    Q = counterterm_basis(sm)
    images = [LocalFunctional.integrate(brst.b(e.functional)) for e in Q.elements[:30]]
    fp = Fingerprinter(sm, values={"xi": 1})
    blocks = [(c.fields, fp.samples(list(c.fields), 3)) for c in leg_classes(sm, 1)[:25]]
    pm = parametric_rows(fp, space, images, blocks)
    for pt in space.random_points(fp.field, 2, seed=5):
        direct = Fingerprinter(sm, values=space.values(pt))
        expect = np.concatenate([np.stack([direct.functional(f, s) for f in images], axis=1) for _, s in blocks])
        assert np.array_equal(pm.at(pt), expect % fp.field.p)


def test_parametric_matrix_algebra():
    f = default_field()
    a = ParametricMatrix(f, (2, 2), {(1, 0): np.array([[1, 2], [3, 4]]), (0, 2): np.eye(2, dtype=np.int64)})
    b = ParametricMatrix.from_values(f, np.array([[5], [6]])).padded(2)
    st = ParametricMatrix.hstack([a, b])
    pt = (3, 7)
    assert np.array_equal(st.at(pt), np.concatenate([a.at(pt), b.at(pt)], axis=1))
    assert np.array_equal(a.at(pt), np.array([[3 + 49, 6], [9, 12 + 49]]))
    assert np.array_equal(a.columns([1]).at(pt)[:, 0], a.at(pt)[:, 1])
    assert np.array_equal(a.scale(-1).at(pt), (-a.at(pt)) % f.p)


def test_reconstruct_known_polynomials(space):
    """Two outputs on graded supports: a real polynomial and a complex one (both embeddings of i)."""
    g_s, g, gp = 0, 1, 2
    sup = space.support(2, 0)

    def prepare(field, i_sign):
        i = field.i if i_sign == 1 else field.p - field.i

        def evaluate(pt):
            x = [int(v) for v in pt]
            first = (field(Fraction(4, 3)) * x[g_s] ** 2 - field(Fraction(1, 36)) * x[gp] ** 2) % field.p
            second = (x[g] * x[gp] * (3 + 5 * i) - x[4]) % field.p      # 3 g g' + 5 i g g' - lambda
            return np.array([first, second])
        return evaluate

    real = reconstruct(prepare, [sup], space)[0]
    e = {k: v for k, v in real.items()}
    assert e == {(2, 0, 0, 0, 0, 0, 0, 0): Fraction(4, 3), (0, 0, 2, 0, 0, 0, 0, 0): Fraction(-1, 36)}
    both = reconstruct(prepare, [sup, sup], space, complex_values=True)
    assert both[1] == {(0, 1, 1, 0, 0, 0, 0, 0): (Fraction(3), Fraction(5)), (0, 0, 0, 0, 1, 0, 0, 0): (-1, 0)}
    assert space.expression(both[1]) == (space.symbols[1] * space.symbols[2] * (3 + 5 * Expression.I)
                                         - space.symbols[4]).expand()
    with pytest.raises(linalg.InconsistentSystemError):
        reconstruct(prepare, [space.support(1, 0), sup], space)       # g_s^2 is not of weight 1
