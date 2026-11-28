"""Gauge groups: explicit generators, structure constants, invariants and invariant tensors.

A gauge group is a product of factors SU(N) and U(1) (:class:`GaugeFactor`). Representations of a factor are
given by a string: ``"1"`` (singlet), ``"F"`` (fundamental), ``"Fbar"`` (antifundamental), ``"A"`` (adjoint) for
SU(N), and a rational charge for U(1).

Conventions (report, 03_setup eq:Tg..eq:invariants):

* SU(N) fundamental generators are the generalised Gell-Mann matrices divided by two, ordered as the Gell-Mann
  matrices for N = 3 and the Pauli matrices for N = 2: for k = 2..N the symmetric and antisymmetric off-diagonal
  pairs (j, k), j < k, then the diagonal ``diag(1, .., 1, -(k-1), 0, ..) / sqrt(2 k (k-1))``. Hence
  ``Tr(T^a T^b) = delta^{ab}/2`` and ``[T^a, T^b] = i f^{abc} T^c``, ``f^{123} = 1``.
* ``Fbar`` carries ``-(T^a)^* = -(T^a)^T``; the adjoint carries ``(T^a)_{bc} = -i f^{abc}``; a U(1) field of
  charge ``q`` carries ``q``.
* Coupling-weighted generators: ``cal(T)^a = g_r T^a`` for ``a`` in factor ``r`` (``g_r`` the factor's coupling);
  ``fg^{abc} = g_r f^{abc}``. The combined adjoint index of a product group runs over the factors in the order
  given, each contributing ``N^2 - 1`` values (one for U(1)).

Explicit tensors are computed over a prime field (:mod:`hepkit_renorm.exact`) in which ``i`` and the square roots
of the diagonal normalisations exist; every invariance and zero test is exact in that field. Rational invariants
are recovered by rational reconstruction (:meth:`PrimeField.reconstruct`).

Symbolic builders (:func:`generator_expr`, :func:`structure_expr`) return coefficients in the vocabulary of
:mod:`hepkit_renorm.tensors`; :class:`GroupEvaluator` turns such coefficients back into explicit arrays.
"""

from __future__ import annotations

import functools
import itertools
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from symbolica import AtomType, Expression

from . import tensors as tn
from .exact import PrimeField

__all__ = [
    "GaugeFactor",
    "GroupEvaluator",
    "Rep",
    "casimir_value",
    "complex_rational",
    "conjugate_rep",
    "dynkin_index",
    "generator_expr",
    "invariant_tensors",
    "product_generators",
    "rep_generators",
    "structure_expr",
    "su_d_tensor",
    "su_generators",
    "su_structure_constants",
]

Rep = str | Fraction
"""A representation of one factor: ``"1"``, ``"F"``, ``"Fbar"``, ``"A"`` (SU(N)) or a charge (U(1))."""


@dataclass(frozen=True)
class GaugeFactor:
    """One factor of the gauge group.

    ``name`` labels the factor and its gauge field (``G``, ``W``, ``B``), ``group`` is ``"SU(N)"`` or ``"U(1)"``,
    ``coupling`` the name of its coupling parameter.
    """

    name: str
    group: str
    coupling: str

    def __post_init__(self) -> None:
        if self.group != "U(1)" and not re.fullmatch(r"SU\((\d+)\)", self.group):
            raise ValueError(f"unsupported gauge group {self.group!r} (SU(N) or U(1))")
        if not self.abelian and self.n < 2:
            raise ValueError("SU(N) needs N >= 2")

    @property
    def abelian(self) -> bool:
        return self.group == "U(1)"

    @property
    def n(self) -> int:
        """N of SU(N); 1 for U(1)."""
        return 1 if self.abelian else int(self.group[3:-1])

    @property
    def dim_adjoint(self) -> int:
        return 1 if self.abelian else self.n**2 - 1

    def rep_dim(self, rep: Rep) -> int:
        """Dimension of a representation of this factor."""
        if self.abelian:
            return 1
        return {"1": 1, "F": self.n, "Fbar": self.n, "A": self.dim_adjoint}[_check_rep(self, rep)]

    def adjoint_space(self) -> tn.Space | None:
        """Index space of the adjoint (``None`` for U(1), whose adjoint index is dropped)."""
        return None if self.abelian else tn.adjoint_space(self.name, self.dim_adjoint)

    def rep_space(self, rep: Rep) -> tn.Space | None:
        """Index space of a matter slot in ``rep`` (``None`` for singlets and U(1))."""
        if self.abelian or rep == "1":
            return None
        if rep == "A":
            return self.adjoint_space()
        return tn.fundamental_space(self.name, self.n)


def _check_rep(factor: GaugeFactor, rep: Rep) -> Rep:
    if factor.abelian:
        return Fraction(rep)
    if rep not in ("1", "F", "Fbar", "A"):
        raise ValueError(f"representation {rep!r} of {factor.group} not supported")
    return rep


def conjugate_rep(factor: GaugeFactor, rep: Rep) -> Rep:
    """Representation of the conjugate field."""
    if factor.abelian:
        return -Fraction(rep)
    return {"1": "1", "F": "Fbar", "Fbar": "F", "A": "A"}[rep]


# ---------------------------------------------------------------------------------------------- explicit SU(N)
def _su_basis(n: int) -> list[tuple[str, int, int]]:
    out = []
    for k in range(1, n):
        for j in range(k):
            out += [("S", j, k), ("A", j, k)]
        out.append(("H", k, 0))
    return out


@functools.cache
def _su_generators_cached(n: int, p: int) -> np.ndarray:
    F = PrimeField(p)
    i = F.i
    half = F(Fraction(1, 2))
    gens = []
    for kind, j, k in _su_basis(n):
        m = [[0] * n for _ in range(n)]
        if kind == "S":
            m[j][k] = m[k][j] = half
        elif kind == "A":
            m[j][k] = (-i * half) % p
            m[k][j] = (i * half) % p
        else:  # diag(1, .., 1, -j, 0, ..) / sqrt(2 j (j + 1)) with j ones
            norm = F.inv(F.sqrt(2 * j * (j + 1)))
            for r in range(j):
                m[r][r] = norm
            m[j][j] = (-j * norm) % p
        gens.append(m)
    return np.array(gens, dtype=object)


def su_generators(n: int, field: PrimeField) -> np.ndarray:
    """Fundamental generators ``T[a, i, j]`` of SU(n), shape ``(n^2-1, n, n)``, over ``field``."""
    return _su_generators_cached(n, field.p).copy()


@functools.cache
def _su_fd_cached(n: int, p: int) -> tuple[np.ndarray, np.ndarray]:
    F = PrimeField(p)
    t = _su_generators_cached(n, p)
    prod = F.einsum("aij,bjk->abik", t, t)
    comm = F.reduce(prod - np.transpose(prod, (1, 0, 2, 3)))
    anti = F.reduce(prod + np.transpose(prod, (1, 0, 2, 3)))
    tr_c = F.einsum("abik,cki->abc", comm, t)
    tr_a = F.einsum("abik,cki->abc", anti, t)
    f = F.reduce(tr_c * ((-2 * F.i) % p))
    dd = F.reduce(tr_a * 2)
    return f, dd


def su_structure_constants(n: int, field: PrimeField) -> np.ndarray:
    """``f^{abc} = -2 i Tr([T^a, T^b] T^c)``, totally antisymmetric."""
    return _su_fd_cached(n, field.p)[0].copy()


def su_d_tensor(n: int, field: PrimeField) -> np.ndarray:
    """``d^{abc} = 2 Tr({T^a, T^b} T^c)``, totally symmetric (zero for SU(2))."""
    return _su_fd_cached(n, field.p)[1].copy()


def rep_generators(factor: GaugeFactor, rep: Rep, field: PrimeField) -> np.ndarray:
    """Unit-coupling generators of ``factor`` in ``rep``: array ``[a, i, j]`` of shape
    ``(dim_adjoint, dim_rep, dim_rep)``."""
    rep = _check_rep(factor, rep)
    if factor.abelian:
        return field.array([[[rep]]])
    n = factor.n
    if rep == "1":
        return field.zeros((factor.dim_adjoint, 1, 1))
    t = su_generators(n, field)
    if rep == "F":
        return t
    if rep == "Fbar":
        return field.reduce(-np.transpose(t, (0, 2, 1)))
    return field.reduce(su_structure_constants(n, field) * ((-field.i) % field.p))


def dynkin_index(factor: GaugeFactor, rep: Rep, field: PrimeField) -> Fraction:
    """``I(R)`` with ``Tr_R(T^a T^b) = I(R) delta^{ab}`` (unit coupling; ``q^2`` for U(1))."""
    t = rep_generators(factor, rep, field)
    tr = field.einsum("aij,bji->ab", t, t)
    val = tr[0, 0]
    expect = field.reduce(np.eye(factor.dim_adjoint, dtype=object) * val)
    if not field.is_zero(field.reduce(tr - expect)):
        raise ArithmeticError("trace form not proportional to delta")
    return field.reconstruct(int(val))


def casimir_value(factor: GaugeFactor, rep: Rep, field: PrimeField) -> Fraction:
    """``C_2(R)`` with ``sum_a T^a T^a = C_2(R) 1`` (unit coupling; ``q^2`` for U(1))."""
    t = rep_generators(factor, rep, field)
    c = field.einsum("aij,ajk->ik", t, t)
    val = c[0, 0]
    if not field.is_zero(field.reduce(c - np.eye(c.shape[0], dtype=object) * val)):
        raise ArithmeticError("Casimir not proportional to the identity")
    return field.reconstruct(int(val))


def invariant_tensors(factor: GaugeFactor, reps: Sequence[Rep], field: PrimeField) -> list[np.ndarray]:
    """Basis of the invariant tensors of ``rep_1 x ... x rep_n`` of one factor, as the exact null space of
    ``sum_k T^a_(k)`` (all ``a``) over ``field``. Tensor ``t[i_1, .., i_n]`` is invariant if
    ``sum_k (T^a_(k))_{i_k j} t[.., j, ..] = 0``. Basis vector ``m`` is normalised to ``1`` at the ``m``-th free
    component in row-major order (deterministic). Dense elimination: intended for products of moderate size
    (a few thousand components)."""
    gens = [rep_generators(factor, r, field) for r in reps]
    dims = [gg.shape[1] for gg in gens]
    size = int(np.prod(dims))
    rows = []
    for a in range(factor.dim_adjoint):
        op = np.zeros((size, size), dtype=object)
        for k, gk in enumerate(gens):
            mats = [np.eye(dd, dtype=object) for dd in dims]
            mats[k] = gk[a]
            m = mats[0]
            for mm in mats[1:]:
                m = np.kron(m, mm)
            op = op + m
        rows.append(field.reduce(op))
    big = np.concatenate(rows, axis=0)
    return [v.reshape(dims) for v in field.nullspace(big)]


def product_generators(factors: Sequence[GaugeFactor], reps: Sequence[Rep], field: PrimeField,
                       couplings: Sequence[int] | None = None) -> np.ndarray:
    """Generators of a multiplet of the product group on the tensor product of its slots (factor order,
    singlets and U(1) contributing dimension one), with the combined adjoint index over all factors (in order,
    ``dim_adjoint`` values each). ``couplings[r]`` multiplies the generators of factor ``r`` (default 1).
    Shape ``(sum dim_adjoint, D, D)``."""
    gens = [rep_generators(fa, r, field) for fa, r in zip(factors, reps)]
    dims = [gg.shape[1] for gg in gens]
    size = int(np.prod(dims))
    out = []
    for r, gg in enumerate(gens):
        c = 1 if couplings is None else couplings[r]
        for a in range(gg.shape[0]):
            mats = [np.eye(dd, dtype=object) for dd in dims]
            mats[r] = gg[a]
            m = mats[0]
            for mm in mats[1:]:
                m = np.kron(m, mm)
            out.append(field.reduce(m * c))
    return np.array(out, dtype=object).reshape(len(out), size, size)


# ---------------------------------------------------------------------------------------------- symbolic builders
def generator_expr(factor: GaugeFactor, rep: Rep, a: Expression | None, i: Expression | None,
                   j: Expression | None) -> Expression:
    """Unit-coupling generator ``(T^a_R)_{ij}`` as a coefficient: ``T(a,i,j)`` (F), ``-T(a,j,i)`` (Fbar),
    ``-i f(a,i,j)`` (A), the charge ``q`` (U(1), all labels ``None``), ``0`` for a singlet."""
    rep = _check_rep(factor, rep)
    if factor.abelian:
        return Expression.num(rep.numerator) / rep.denominator
    if rep == "1":
        return Expression.num(0)
    if rep == "F":
        return tn.T(a, i, j)
    if rep == "Fbar":
        return -tn.T(a, j, i)
    return -tn.I * tn.f(a, i, j)


def structure_expr(factor: GaugeFactor, a: Expression | None, b: Expression | None,
                   c: Expression | None) -> Expression:
    """Unit-coupling structure constants ``f^{abc}`` of the factor (``0`` for U(1))."""
    return Expression.num(0) if factor.abelian else tn.f(a, b, c)


# ---------------------------------------------------------------------------------------------- evaluation
_GROUP_HEADS = {tn.T.get_name(), tn.f.get_name(), tn.dsym.get_name(), tn.epsf.get_name()}


class GroupEvaluator:
    """Explicit evaluation of group heads (``T``, ``f``, ``dsym``, ``epsf``, ``delta`` on group labels) over a
    prime field.

    ``factors`` maps factor names to :class:`GaugeFactor`. Label values are 0-based component numbers.
    """

    def __init__(self, factors: Mapping[str, GaugeFactor], field: PrimeField, n_flavours: int = 1) -> None:
        self.factors = dict(factors)
        self.field = field
        self.n_flavours = n_flavours

    def label_dim(self, label: Expression) -> int:
        key = tn.label_space(label)
        if key == "flavour":
            return self.n_flavours
        if key in ("lorentz", "spinor"):
            raise ValueError("group evaluator does not evaluate Lorentz or spinor labels")
        kind, fac = key.split(":")
        fa = self.factors[fac]
        return fa.dim_adjoint if kind == "adjoint" else fa.n

    def _head_array(self, fn: Expression) -> tuple[np.ndarray, list[Expression]]:
        name = fn.get_name()
        args = list(fn)
        F = self.field
        if name == tn.delta.get_name():
            n = self.label_dim(args[0])
            return F.array(np.eye(n, dtype=int).tolist()), args
        fac = self.factors[tn.label_space(args[0]).split(":")[1]]
        if name == tn.T.get_name():
            return su_generators(fac.n, F), args
        if name == tn.f.get_name():
            return su_structure_constants(fac.n, F), args
        if name == tn.dsym.get_name():
            return su_d_tensor(fac.n, F), args
        if name == tn.epsf.get_name():
            n = fac.n
            e = np.zeros((n,) * n, dtype=object)
            for perm in itertools.permutations(range(n)):
                inv = sum(1 for x, y in itertools.combinations(perm, 2) if x > y)
                e[perm] = 1 if inv % 2 == 0 else F.p - 1
            return e, args
        raise ValueError(f"not a group head: {name}")

    def number(self, num: Expression) -> int:
        """Image of a Symbolica complex rational in the field."""
        z = complex_rational(num)
        return (self.field(z[0]) + self.field(z[1]) * self.field.i) % self.field.p

    def evaluate(self, term: Expression, open_labels: Sequence[Expression],
                 values: Mapping[str, int] | None = None) -> np.ndarray:
        """Explicit array of a product term of group heads, numbers and symbols, indexed by ``open_labels``.

        ``values`` maps symbol names (as printed by ``get_name()``) to field elements; every symbol of the term
        must have a value. Sums are not accepted: evaluate term by term.
        """
        F = self.field
        values = values or {}
        scalar = 1
        arrays: list[np.ndarray] = []
        labs: list[list[str]] = []
        factors = list(term) if term.get_type() == AtomType.Mul else [term]
        for fa in factors:
            t = fa.get_type()
            if t == AtomType.Num:
                scalar = scalar * self.number(fa) % F.p
            elif t == AtomType.Var:
                scalar = scalar * values[fa.get_name()] % F.p
            elif t == AtomType.Fn:
                arr, args = self._head_array(fa)
                arrays.append(arr)
                labs.append([a.to_canonical_string() for a in args])
            elif t == AtomType.Pow:
                base, ex = list(fa)
                k = int(str(ex))
                if base.get_type() == AtomType.Var:
                    scalar = scalar * pow(values[base.get_name()], k, F.p) % F.p
                elif base.get_type() == AtomType.Fn:
                    arr, args = self._head_array(base)
                    for _ in range(k):
                        arrays.append(arr)
                        labs.append([a.to_canonical_string() for a in args])
                else:
                    raise ValueError(f"cannot evaluate {fa}")
            else:
                raise ValueError(f"cannot evaluate {fa} (expand first)")
        letters: dict[str, str] = {}
        pool = iter("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")
        for ls in labs:
            for lab in ls:
                if lab not in letters:
                    letters[lab] = next(pool)
        out_keys = [o.to_canonical_string() for o in open_labels]
        for o in out_keys:
            if o not in letters:
                letters[o] = next(pool)
        if not arrays:
            if open_labels:
                raise ValueError("open labels absent from the term")
            return np.array(scalar, dtype=object)
        spec = ",".join("".join(letters[x] for x in ls) for ls in labs) + "->" + "".join(letters[o] for o in out_keys)
        res = F.einsum(spec, *arrays)
        return F.reduce(res * scalar)


def complex_rational(num: Expression) -> tuple[Fraction, Fraction]:
    """Real and imaginary part of a Symbolica complex rational number."""
    if num.get_type() != AtomType.Num:
        raise ValueError(f"not a number: {num}")
    re_part = ((num + num.conj()) / 2).expand()
    im_part = ((num - num.conj()) / (2 * tn.I)).expand()
    return Fraction(re_part.to_canonical_string()), Fraction(im_part.to_canonical_string())
