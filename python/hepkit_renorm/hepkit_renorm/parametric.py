"""Exact dependence on the parameters of the model (couplings, masses): parametric fingerprints and graded
reconstruction.

Every linear decision of the Slavnov-Taylor layer is taken over ``GF(p)`` at a point of parameter space
(:mod:`hepkit_renorm.fingerprint`, :mod:`hepkit_renorm.linalg`). The counterterms written into a model need their
coefficients as exact polynomials in the parameters. This module provides the two pieces that make that cheap and
exact for any model:

* **Parametric fingerprints** (:func:`parametric_rows`). A local functional ``F`` is a polynomial in the parameter
  symbols, ``F = sum_e p^e F_e`` with parameter-free ``F_e`` (:meth:`ParameterSpace.expand`). The fingerprints of
  the ``F_e`` are taken once, with the ordinary :class:`~hepkit_renorm.fingerprint.Fingerprinter`; the fingerprint
  of ``F`` at any parameter point is then the polynomial ``sum_e point^e fp(F_e)`` (:class:`ParametricMatrix`).
  The samples of a fingerprinter do not depend on the parameter values, so all points share the same rows.
* **Graded reconstruction** (:func:`reconstruct`). The outputs of an exact computation (solutions of linear systems
  assembled from parametric matrices) are evaluated at random points of ``GF(p)^n`` and interpolated on an exact
  monomial support; the rational coefficients are recovered over several primes
  (:func:`hepkit_renorm.linalg.reconstruct_rationals`: CRT, confirmed by one more prime). Complex coefficients
  are split into real and imaginary parts with the two embeddings of ``i``
  (:func:`hepkit_renorm.linalg.split_gaussian`).

Gradings (derived, :class:`ParameterSpace`)
-------------------------------------------
Two gradings of the parameters bound the supports. Both are read off the tree-level action ``S_0``:

* the *loop weight* ``w_p``: every monomial of ``S_0`` with ``n`` legs has parameter weight ``n - 2`` (couplings of
  three-point vertices weight 1, of four-point vertices weight 2, masses weight 0); a diagram of order ``l`` with
  ``n`` external legs then has weight ``n - 2 + 2 l``;
* the *mass dimension* ``dim_p``: every monomial of ``S_0`` has dimension 4 (``4 + g dim s`` at ghost number
  ``g``), shared between fields, derivatives and parameters.

The ``grade`` of a monomial is its parameter weight minus ``n - 2``; the ``excess`` its required parameter
dimension, ``4 + g dim s - sum dim(fields) - #derivatives``, minus the dimension of its parameters. At order ``l``
every counterterm and breaking functional has grade ``2 l`` and excess 0, so the coefficient of a basis element of
grade ``g_E`` and excess ``x_E`` is a polynomial of weight ``2 l - g_E`` and dimension ``x_E``
(:meth:`ParameterSpace.support`). Both systems of equations for ``w`` and ``dim`` are solved exactly; a model
without such gradings raises :class:`GradingError` (pass ``degree`` bounds instead).

Parameters may be *fixed* to rational values (the gauge-fixing parameter in Feynman gauge, or all of them for a
computation at one rational point, e.g. with loop data computed at that point): fixed parameters are substituted
before the expansion and are not variables of the reconstruction.
"""

from __future__ import annotations

import weakref
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from symbolica import AtomType, Expression

from . import linalg
from .classes import brst_dimension
from .exact import PrimeField
from .fingerprint import Fingerprinter, Samples
from .functional import LocalFunctional
from .spec import Model

__all__ = [
    "GradingError",
    "ParameterSpace",
    "ParametricMatrix",
    "Polynomial",
    "interpolate_outputs",
    "parametric_rows",
    "reconstruct",
]

Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Fraction | tuple[Fraction, Fraction]]
"""A polynomial in the free parameters: ``{exponents: coefficient}``; coefficients are :class:`Fraction` (real
reconstruction) or ``(re, im)`` pairs (complex reconstruction). Zero coefficients are omitted."""


class GradingError(ValueError):
    """The tree-level action admits no loop-weight or dimension grading of the parameters."""


def _num(q: Fraction | int) -> Expression:
    q = Fraction(q)
    return Expression.num(q.numerator) / q.denominator


def _monomial_exponents(key: Expression, symbols: Sequence[Expression]) -> Exponent:
    """Exponents of a product of powers of ``symbols`` (a key of ``Expression.coefficient_list``)."""
    names = {s.get_name(): k for k, s in enumerate(symbols)}
    out = [0] * len(symbols)
    if key.get_type() == AtomType.Num:
        return tuple(out)
    factors = list(key) if key.get_type() == AtomType.Mul else [key]
    for f in factors:
        if f.get_type() == AtomType.Pow:
            base, exp = list(f)
            out[names[base.get_name()]] += int(exp.to_canonical_string())
        else:
            out[names[f.get_name()]] += 1
    return tuple(out)


# ---------------------------------------------------------------------------------------------- parameter space
class ParameterSpace:
    """The free parameters of a computation, the fixed ones and the derived gradings (module docstring).

    ``fixed``: ``{name: rational value}`` (specification names, e.g. ``{"xi": 1}``); every other coupling and real
    parameter of the model is free. ``action``: the tree-level action used for the gradings (``S_0``).
    """

    def __init__(self, model: Model, action: LocalFunctional, fixed: Mapping[str, Fraction | int] | None = None
                 ) -> None:
        self.model = model
        self.fixed = {k: Fraction(v) for k, v in (fixed or {}).items()}
        every = [f.coupling for f in model.factors.values()] + list(model.parameters)
        unknown = set(self.fixed) - set(every)
        if unknown:
            raise KeyError(f"unknown parameters {sorted(unknown)}")
        for name, par in model.parameters.items():
            if not par.real and name not in self.fixed:
                raise NotImplementedError(f"complex parameter {name!r}: conj(p) as a variable is not implemented")
        self.names: tuple[str, ...] = tuple(n for n in every if n not in self.fixed)
        self.symbols: tuple[Expression, ...] = tuple(model.param(n) for n in self.names)
        self._fixed_symbols = [(model.param(n), _num(v)) for n, v in self.fixed.items()]
        self._cache: dict[int, tuple[weakref.ref, dict[Exponent, LocalFunctional]]] = {}
        self.weights, self.dimensions = self._gradings(action)

    def __len__(self) -> int:
        return len(self.names)

    # ------------------------------------------------------------------ expansion
    def substitute(self, expr: Expression) -> Expression:
        """``expr`` with the fixed parameters replaced by their values."""
        for sym, val in self._fixed_symbols:
            expr = expr.replace(sym, val)
        return expr

    def expand(self, functional: LocalFunctional) -> dict[Exponent, LocalFunctional]:
        """``F = sum_e p^e F_e`` with parameter-free ``F_e`` (fixed parameters substituted); cached per object."""
        hit = self._cache.get(id(functional))
        if hit is not None and hit[0]() is functional:
            return hit[1]
        expr = self.substitute(functional.expr).expand()
        out: dict[Exponent, LocalFunctional] = {}
        if not functional.is_zero():
            pairs = expr.coefficient_list(*self.symbols) if self.symbols else [(Expression.num(1), expr)]
            for key, coeff in pairs:
                e = _monomial_exponents(key, self.symbols)
                part = LocalFunctional(functional.fields, coeff)
                if not part.is_zero():
                    out[e] = part
        self._cache[id(functional)] = (weakref.ref(functional), out)
        return out

    # ------------------------------------------------------------------ gradings
    def _legs_and_dimension(self, functional: LocalFunctional):
        fields = self.model.fields
        ds = brst_dimension(self.model)
        for e, part in self.expand(functional).items():
            for m in part.monomials():
                legs = m.fields()
                g = sum(fields[n].ghost_number for n in legs)
                dim = sum((fields[n].dimension for n in legs), Fraction(0))
                ders = sum(len(o.derivatives) for o in m.occurrences)
                yield e, len(legs), 4 + g * ds - dim - ders

    def _gradings(self, action: LocalFunctional) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
        rows_w, rhs_w, rows_d, rhs_d = [], [], [], []
        seen = set()
        for e, n, required in self._legs_and_dimension(action):
            key = (e, n, required)
            if key in seen:
                continue
            seen.add(key)
            rows_w.append(e), rhs_w.append(n - 2)
            rows_d.append(e), rhs_d.append(required)
        return self._solve_grading(rows_w, rhs_w, "loop weight"), self._solve_grading(rows_d, rhs_d, "dimension")

    def _solve_grading(self, rows, rhs, what: str) -> tuple[Fraction, ...]:
        k = len(self.names)
        if k == 0:
            return ()
        a = np.array(rows, dtype=object).reshape(len(rows), k)
        b = np.array([Fraction(x) for x in rhs], dtype=object)

        def compute(field: PrimeField) -> np.ndarray:
            am = linalg.as_field_matrix(np.vectorize(field, otypes=[object])(a), field)
            bm = np.array([field(x) for x in b], dtype=object)
            try:
                sol = linalg.solve(am, bm, field, require_unique=True)
            except linalg.InconsistentSystemError as exc:
                raise GradingError(f"S_0 has no {what} grading of the parameters") from exc
            except ValueError as exc:
                raise GradingError(f"the {what} of some parameter is not fixed by S_0") from exc
            return sol

        return tuple(Fraction(x) for x in linalg.reconstruct_rationals(compute))

    def grade(self, functional: LocalFunctional) -> tuple[Fraction, Fraction]:
        """``(grade, excess)`` of a functional (module docstring); raises :class:`GradingError` if its monomials
        do not share them."""
        found = set()
        for e, n, required in self._legs_and_dimension(functional):
            w = sum(Fraction(x) * y for x, y in zip(e, self.weights))
            d = sum(Fraction(x) * y for x, y in zip(e, self.dimensions))
            found.add((w - (n - 2), required - d))
        if len(found) != 1:
            raise GradingError(f"functional with several (grade, excess) values {sorted(found)}")
        return found.pop()

    def support(self, weight: Fraction | int, dimension: Fraction | int) -> list[Exponent]:
        """All exponent vectors (non-negative) of the free parameters with loop weight ``weight`` and dimension
        ``dimension``, in a deterministic order. Raises :class:`GradingError` if the set is infinite (a parameter of
        weight and dimension zero)."""
        k = len(self.names)
        if k == 0:
            return [()]
        weight, dimension = Fraction(weight), Fraction(dimension)
        for w, d in zip(self.weights, self.dimensions):
            if w < 0 or d < 0:
                raise GradingError("negative parameter weights or dimensions are not supported")
            if w == 0 and d == 0:
                raise GradingError("a free parameter of weight and dimension zero: the support is not finite")
        out: list[Exponent] = []

        def rec(i: int, w: Fraction, d: Fraction, cur: tuple[int, ...]) -> None:
            if i == k:
                if w == 0 and d == 0:
                    out.append(cur)
                return
            wi, di = self.weights[i], self.dimensions[i]
            n = 0
            while w - n * wi >= 0 and d - n * di >= 0:
                rec(i + 1, w - n * wi, d - n * di, cur + (n,))
                n += 1

        if weight >= 0 and dimension >= 0:
            rec(0, weight, dimension, ())
        return sorted(out, key=lambda e: (sum(e), [-x for x in e]))

    # ------------------------------------------------------------------ points and values
    def values(self, point: Sequence[int]) -> dict[str, object]:
        """Fingerprinter values (``Fingerprinter(values=...)``) for a point of the free parameters (field
        elements) together with the fixed values."""
        out: dict[str, object] = dict(self.fixed)
        out.update({n: int(x) for n, x in zip(self.names, point)})
        return out

    def random_points(self, field: PrimeField, n: int, seed: int = 0) -> np.ndarray:
        """``n`` deterministic random points with non-zero coordinates."""
        rng = np.random.default_rng([seed, field.p, len(self.names), 7])
        return rng.integers(1, field.p, size=(n, len(self.names)), dtype=np.int64)

    def expression(self, poly: Mapping[Exponent, Fraction | tuple[Fraction, Fraction]]) -> Expression:
        """The polynomial as a Symbolica expression in the parameter symbols (``i`` as
        :data:`hepkit_renorm.tensors.I` for complex coefficients)."""
        from .tensors import I

        total = Expression.num(0)
        for e, c in poly.items():
            term = Expression.num(1)
            for sym, k in zip(self.symbols, e):
                term = term * sym ** k
            if isinstance(c, tuple):
                coeff = _num(c[0]) + I * _num(c[1])
            else:
                coeff = _num(c)
            total = total + coeff * term
        return total.expand()


# ---------------------------------------------------------------------------------------------- parametric matrices
@dataclass
class ParametricMatrix:
    """A matrix over ``GF(p)`` that is a polynomial in the free parameters: ``M(point) = sum_e point^e terms[e]``
    (exponents may be negative: the parameters are non-zero at every point)."""

    field: PrimeField
    shape: tuple[int, int]
    terms: dict[Exponent, np.ndarray]

    def at(self, point: Sequence[int]) -> np.ndarray:
        p = self.field.p
        out = np.zeros(self.shape, dtype=np.int64)
        for e, m in self.terms.items():
            c = 1
            for x, k in zip(point, e):
                c = c * pow(int(x), int(k), p) % p
            out = (out + m * c) % p
        return out

    @property
    def constant(self) -> bool:
        return all(not any(e) for e in self.terms)

    def columns(self, cols: Sequence[int]) -> ParametricMatrix:
        cols = list(cols)
        return ParametricMatrix(self.field, (self.shape[0], len(cols)), {e: m[:, cols] for e, m in self.terms.items()})

    def rows(self, rows: Sequence[int]) -> ParametricMatrix:
        rows = list(rows)
        return ParametricMatrix(self.field, (len(rows), self.shape[1]), {e: m[rows] for e, m in self.terms.items()})

    def scale(self, c: int) -> ParametricMatrix:
        p = self.field.p
        return ParametricMatrix(self.field, self.shape, {e: m * (c % p) % p for e, m in self.terms.items()})

    @staticmethod
    def hstack(parts: Sequence[ParametricMatrix]) -> ParametricMatrix:
        return ParametricMatrix._stack(parts, axis=1)

    @staticmethod
    def vstack(parts: Sequence[ParametricMatrix]) -> ParametricMatrix:
        return ParametricMatrix._stack(parts, axis=0)

    @staticmethod
    def _stack(parts: Sequence[ParametricMatrix], axis: int) -> ParametricMatrix:
        field = parts[0].field
        keys = sorted({e for m in parts for e in m.terms})
        terms = {}
        for e in keys:
            blocks = [m.terms.get(e, np.zeros(m.shape, dtype=np.int64)) for m in parts]
            terms[e] = np.concatenate(blocks, axis=axis)
        rows = sum(m.shape[0] for m in parts) if axis == 0 else parts[0].shape[0]
        cols = sum(m.shape[1] for m in parts) if axis == 1 else parts[0].shape[1]
        return ParametricMatrix(field, (rows, cols), terms)

    @staticmethod
    def from_values(field: PrimeField, values: np.ndarray) -> ParametricMatrix:
        """A constant matrix (no parameter dependence, or values at the computation's only point)."""
        values = np.asarray(values, dtype=np.int64) % field.p
        if values.ndim == 1:
            values = values[:, None]
        return ParametricMatrix(field, values.shape, {(): values} if values.size else {})

    def padded(self, n_vars: int) -> ParametricMatrix:
        """Exponent keys padded to ``n_vars`` entries (constant keys ``()`` become all zeros)."""
        terms: dict[Exponent, np.ndarray] = {}
        for e, m in self.terms.items():
            key = tuple(e) + (0,) * (n_vars - len(e))
            terms[key] = (terms[key] + m) % self.field.p if key in terms else m
        return ParametricMatrix(self.field, self.shape, terms)


def parametric_rows(fp: Fingerprinter, space: ParameterSpace, functionals: Sequence[LocalFunctional],
                    blocks: Sequence[tuple[tuple[str, ...], Samples]]) -> ParametricMatrix:
    """The fingerprints of ``functionals`` (columns) on the sample blocks (rows, in order) as a polynomial in the
    free parameters of ``space``."""
    n_rows = sum(s.n for _, s in blocks)
    terms: dict[Exponent, np.ndarray] = {}
    k = len(space)
    for j, f in enumerate(functionals):
        for e, part in space.expand(f).items():
            col = np.concatenate([fp.functional(part, s) for _, s in blocks]) if blocks else np.zeros(0, np.int64)
            if not np.any(col):
                continue
            key = tuple(e) + (0,) * (k - len(e))
            if key not in terms:
                terms[key] = np.zeros((n_rows, len(functionals)), dtype=np.int64)
            terms[key][:, j] = (terms[key][:, j] + col) % fp.field.p
    return ParametricMatrix(fp.field, (n_rows, len(functionals)), terms)


# ---------------------------------------------------------------------------------------------- reconstruction
def reconstruct(make_evaluator: Callable[[PrimeField, int], Callable[[np.ndarray], np.ndarray]],
                supports: Sequence[Sequence[Exponent]], space: ParameterSpace, *, complex_values: bool = False,
                max_primes: int = 6, seed: int = 0, extra_points: int = 3) -> list[Polynomial]:
    """Polynomials with rational (or complex rational) coefficients from exact evaluations.

    ``make_evaluator(field, i_sign)`` prepares the computation over a prime (and an embedding of ``i``) and returns
    ``evaluate(point) -> outputs`` (a vector over the field, one entry per output). Output ``j`` is a polynomial on
    ``supports[j]``. The same ``max(len(support)) + extra_points`` random points serve all outputs; the extra
    points test the supports (an inconsistent interpolation raises
    :class:`~hepkit_renorm.linalg.InconsistentSystemError`: the outputs are not polynomials on their supports).
    Coefficients are reconstructed over primes until one more prime confirms them.
    """
    supports = [tuple(tuple(e) for e in s) for s in supports]
    offsets = np.cumsum([0] + [len(s) for s in supports])

    def coefficients(field: PrimeField, i_sign: int) -> np.ndarray:
        return interpolate_outputs(make_evaluator(field, i_sign), supports, space, field, seed=seed,
                                   extra_points=extra_points)

    def compute(field: PrimeField) -> np.ndarray:
        if not complex_values:
            return coefficients(field, 1)
        plus, minus = coefficients(field, 1), coefficients(field, -1)
        re, im = linalg.split_gaussian(plus.astype(np.int64), minus.astype(np.int64), field)
        return np.concatenate([re, im])

    flat = linalg.reconstruct_rationals(compute, max_primes=max_primes) if offsets[-1] else np.zeros(0, object)
    n = int(offsets[-1])
    out: list[Polynomial] = []
    for j, sup in enumerate(supports):
        poly: Polynomial = {}
        for k, e in enumerate(sup):
            pos = int(offsets[j]) + k
            if complex_values:
                c = (Fraction(flat[pos]), Fraction(flat[n + pos]))
                if c != (0, 0):
                    poly[e] = c
            elif flat[pos] != 0:
                poly[e] = Fraction(flat[pos])
        out.append(poly)
    return out


def interpolate_outputs(evaluate: Callable[[np.ndarray], np.ndarray], supports: Sequence[Sequence[Exponent]],
                        space: ParameterSpace, field: PrimeField, *, seed: int = 0, extra_points: int = 3
                        ) -> np.ndarray:
    """The images over ``field`` of the coefficients of outputs that are polynomials on ``supports``: ``evaluate``
    at ``max(len(support)) + extra_points`` deterministic random points and solve the interpolation systems
    exactly (one prime of :func:`reconstruct`; the flat result lists the coefficients of output ``j`` on
    ``supports[j]`` in order). Raises :class:`~hepkit_renorm.linalg.InconsistentSystemError` if an output is not a
    polynomial on its support."""
    supports = [tuple(tuple(e) for e in s) for s in supports]
    distinct = sorted(set(supports), key=lambda s: (len(s), s))
    npts = max([len(s) for s in supports] + [0]) + extra_points
    offsets = np.cumsum([0] + [len(s) for s in supports])
    points = space.random_points(field, npts, seed) if len(space) else np.zeros((1, 0), dtype=np.int64)
    values = np.stack([np.asarray(evaluate(pt), dtype=np.int64) % field.p for pt in points])
    out = np.zeros(int(offsets[-1]), dtype=object)
    p = field.p
    for sup in distinct:
        idx = [j for j, s in enumerate(supports) if s == sup]
        if not sup:
            bad = [j for j in idx if np.any(values[:, j])]
            if bad:
                raise linalg.InconsistentSystemError(bad)
            continue
        a = np.array([[_power(pt, e, p) for e in sup] for pt in points], dtype=np.int64)
        sol = linalg.solve(a, values[:, idx], field, require_unique=True)
        for col, j in enumerate(idx):
            out[offsets[j]:offsets[j + 1]] = sol[:, col]
    return out


def _power(point: Sequence[int], e: Exponent, p: int) -> int:
    c = 1
    for x, k in zip(point, e):
        c = c * pow(int(x), int(k), p) % p
    return c


def evaluate_polynomial(poly: Mapping[Exponent, Fraction], point: Sequence[int], field: PrimeField) -> int:
    """The image of a (real) polynomial at a point of ``GF(p)``."""
    total = 0
    for e, c in poly.items():
        total = (total + field(c) * _power(point, e, field.p)) % field.p
    return total


def iter_exponents(polys: Iterable[Mapping[Exponent, object]]) -> list[Exponent]:
    """The union of the supports of several polynomials (sorted)."""
    return sorted({e for p in polys for e in p}, key=lambda e: (sum(e), [-x for x in e]))

