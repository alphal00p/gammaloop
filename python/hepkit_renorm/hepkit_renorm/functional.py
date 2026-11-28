"""Local functionals of fields and sources with exact Grassmann bookkeeping.

Representation
--------------
A :class:`Monomial` is ``coefficient * X_1 X_2 ... X_n``: a Grassmann-even coefficient (a Symbolica expression in
the vocabulary of :mod:`hepkit_renorm.tensors`) times an *ordered* product of field occurrences
``X_k = d_{mu_1} ... d_{mu_r} Phi_k`` (:class:`Occurrence`: field name, one label per index slot of the field, the
Lorentz labels of the partial derivatives). The order of the occurrences is the order of the Grassmann factors.
Every label is explicit; a label that appears twice in a monomial is contracted (with the metric for Lorentz
labels), a label that appears once is free.

A :class:`Density` is a sum of monomials (possibly with free labels, e.g. a field variation ``s Phi`` or an
equation of motion). A :class:`LocalFunctional` is the integral ``int d^d x`` of a density without free labels;
only for it are integration by parts, functional derivatives and vertex functions defined.

Internally a density is one canonical Symbolica expression, ``sum coefficient * occ(0, Phi_1, slots..) *
dv(0, mu) * occ(1, ...)``: occurrences are sorted by the field order of the :class:`FieldTable` (with the
Grassmann sign of the reordering), dummy labels are renamed canonically and equal monomials are merged. Two
densities are equal (``==``) iff their canonical expressions are identical; this is exact but does not use
integration by parts or tensor identities (use :meth:`LocalFunctional.vanishes` for that).

Sign conventions (report 02_method sec:m-objects, 03b sec:st-signs)
--------------------------------------------------------------------
* Left derivative: ``F[X + dX] = F + int dX dL F/dX``; right derivative: ``F + int (dR F/dX) dX``.
* A graded derivation of parity ``e`` acts from the left: ``delta(X_1 ... X_n) = sum_k (-1)^{e(|X_1|+..+|X_{k-1}|)}
  X_1 .. (delta X_k) .. X_n``, and commutes with partial derivatives.
* Vertex function: ``V_Y(Phi_1(p_1), ..., Phi_n(p_n))`` is the coefficient of ``theta_1 ... theta_n`` in
  ``Y[sum_k theta_k e_k exp(-i p_k.x)]`` with the polarisations stripped (free labels of the legs) and
  ``(2 pi)^d delta(sum p)`` removed: ``d_mu -> -i p_mu`` on the leg of incoming momentum ``p``. The Feynman rule is
  ``i V``. Exchanging two odd legs changes the sign of ``V``.
* Hermitian conjugation reverses the order of the occurrences (:mod:`hepkit_renorm.fields`).
"""

from __future__ import annotations

import itertools
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from dataclasses import field as dc_field

from symbolica import AtomType, Expression, S

from . import tensors as tn
from .fields import FieldTable

__all__ = [
    "Density",
    "Derivation",
    "FieldImage",
    "Leg",
    "LocalFunctional",
    "Monomial",
    "Occurrence",
    "field_symbol",
    "standard_legs",
]

_occ = S(f"{tn.NS}::occ")
_dv = S(f"{tn.NS}::dv")
_FIELD_NS = f"{tn.NS}_field::"
_FIELD_SYMS: dict[str, Expression] = {}


def field_symbol(name: str) -> Expression:
    """The Symbolica symbol that names a field inside the canonical expression."""
    if name not in _FIELD_SYMS:
        _FIELD_SYMS[name] = S(_FIELD_NS + name)
    return _FIELD_SYMS[name]


@dataclass(frozen=True)
class Occurrence:
    """``d_{derivatives} field[slots]``: one Grassmann factor of a monomial.

    ``slots`` holds one label per index slot of the field (order of :attr:`FieldInfo.slots`); ``derivatives``
    the Lorentz labels of the partial derivatives acting on the field (their order is irrelevant).
    """

    field: str
    slots: tuple[Expression, ...] = ()
    derivatives: tuple[Expression, ...] = ()


@dataclass(frozen=True)
class Monomial:
    """``coefficient * occurrences[0] * occurrences[1] * ...`` (ordered Grassmann product)."""

    coefficient: Expression
    occurrences: tuple[Occurrence, ...] = ()

    def fields(self) -> tuple[str, ...]:
        return tuple(o.field for o in self.occurrences)

    def separated(self) -> Monomial:
        """The same monomial with every slot and derivative label of every occurrence replaced by a fresh label
        tied to the old one by ``g`` (Lorentz) or ``delta`` in the coefficient: afterwards no two occurrences
        share a label directly, which is the form needed to strip polarisations."""
        return _separate(self)


@dataclass(frozen=True)
class Leg:
    """An external leg of a vertex function: field, the (free) labels of its slots, momentum number ``k``
    (its incoming momentum is ``p(k, mu)``)."""

    field: str
    labels: tuple[Expression, ...]
    momentum: int


def standard_legs(fields: FieldTable, names: Sequence[str]) -> list[Leg]:
    """Legs ``names[0], names[1], ...`` with momenta ``1, 2, ...`` and the deterministic labels
    ``space.label(10000 + 10 k + j)`` (leg ``k``, slot ``j``), which never collide with canonical or fresh ones."""
    legs = []
    for k, name in enumerate(names, start=1):
        labs = tuple(s.space.label(10_000 + 10 * k + j) for j, s in enumerate(fields[name].slots))
        legs.append(Leg(name, labs, k))
    return legs


# ---------------------------------------------------------------------------------------------- encoding
def _encode(m: Monomial) -> Expression:
    e = m.coefficient
    for pos, o in enumerate(m.occurrences):
        e = e * _occ(pos, field_symbol(o.field), *o.slots)
        for mu in o.derivatives:
            e = e * _dv(pos, mu)
    return e


def _terms(expr: Expression) -> list[Expression]:
    if expr.get_type() == AtomType.Add:
        return list(expr)
    if expr.get_type() == AtomType.Num and expr.to_canonical_string() == "0":
        return []
    return [expr]


def _head_name(power: Expression) -> str | None:
    base = next(iter(power))
    return base.get_name() if base.get_type() == AtomType.Fn else None


def _decode(term: Expression) -> Monomial:
    factors = list(term) if term.get_type() == AtomType.Mul else [term]
    coeff = Expression.num(1)
    occs: dict[int, list] = {}
    ders: dict[int, list] = {}
    for fa in factors:
        if fa.get_type() == AtomType.Fn and fa.get_name() == _occ.get_name():
            args = list(fa)
            pos = int(args[0].to_canonical_string())
            occs[pos] = [args[1].get_name()[len(_FIELD_NS):], tuple(args[2:])]
        elif fa.get_type() == AtomType.Fn and fa.get_name() == _dv.get_name():
            args = list(fa)
            ders.setdefault(int(args[0].to_canonical_string()), []).append(args[1])
        elif fa.get_type() == AtomType.Pow and _head_name(fa) == _dv.get_name():
            base, ex = list(fa)
            args = list(base)
            ders.setdefault(int(args[0].to_canonical_string()), []).extend([args[1]] * int(ex.to_canonical_string()))
        elif fa.get_type() == AtomType.Pow and _head_name(fa) == _occ.get_name():
            raise ValueError("repeated occurrence in a monomial")
        else:
            coeff = coeff * fa
    n = len(occs)
    if sorted(occs) != list(range(n)):
        raise ValueError("occurrence positions are not 0..n-1")
    return Monomial(coeff, tuple(Occurrence(occs[k][0], occs[k][1], tuple(ders.get(k, ()))) for k in range(n)))


def _perm_sign(seq: Sequence[int]) -> int:
    inv = sum(1 for i, j in itertools.combinations(range(len(seq)), 2) if seq[i] > seq[j])
    return -1 if inv % 2 else 1


def _odd_sign(fields: FieldTable, occs: Sequence[Occurrence], order: Sequence[int]) -> int:
    """Grassmann sign of reordering ``occs`` into ``[occs[i] for i in order]``."""
    odd_positions = [i for i in order if fields.odd(occs[i].field)]
    return _perm_sign(odd_positions)


def _canonical_term(fields: FieldTable, m: Monomial) -> Expression:
    """Canonical encoded form of one single-term monomial (0 if it vanishes by Grassmann antisymmetry)."""
    occs = m.occurrences
    base = sorted(range(len(occs)), key=lambda i: fields.order(occs[i].field))
    groups: list[list[int]] = []
    for i in base:
        if groups and occs[groups[-1][0]].field == occs[i].field:
            groups[-1].append(i)
        else:
            groups.append([i])
    best: tuple[str, Expression, set[int]] | None = None
    for choice in itertools.product(*(itertools.permutations(gr) for gr in groups)):
        order = [i for gr in choice for i in gr]
        sign = _odd_sign(fields, occs, order)
        mono = Monomial(m.coefficient, tuple(occs[i] for i in order))
        canon = tn.canonical_relabel(_encode(mono))
        key = canon.to_canonical_string()
        if best is None or key < best[0]:
            best = (key, canon, {sign})
        elif key == best[0]:
            best[2].add(sign)
    assert best is not None
    if len(best[2]) == 2:
        return Expression.num(0)
    return best[1] * best[2].pop()


def _normalize(fields: FieldTable, expr: Expression) -> Expression:
    total = Expression.num(0)
    for term in _terms(expr.expand()):
        total = total + _canonical_term(fields, _decode(term))
    return total


def _fresh_dummies(expr: Expression, keep: Iterable[Expression] = ()) -> Expression:
    """Rename every non-free label of each term (except ``keep``) to fresh labels."""
    keep_s = {k.to_canonical_string() for k in keep}
    out = Expression.num(0)
    for term in _terms(expr):
        counts: dict[str, list] = {}
        for lab in tn.labels_in(term):
            counts.setdefault(lab.to_canonical_string(), [lab, 0])[1] += 1
        mapping = []
        for k, (lab, c) in counts.items():
            if c >= 2 and k not in keep_s:
                mapping.append((lab, tn.fresh(_space_of(lab))))
        out = out + tn.rename(term, mapping)
    return out


def _space_of(label: Expression) -> tn.Space:
    key = tn.label_space(label)
    if key == "lorentz":
        return tn.LORENTZ
    if key == "spinor":
        return tn.SPINOR
    if key == "flavour":
        return tn.Space("flavour")
    kind, fac = key.split(":")
    return tn.Space("adjoint" if kind == "adjoint" else "fund", fac)


def _identity(a: Expression, b: Expression) -> Expression:
    return tn.g(a, b) if tn.label_space(a) == "lorentz" else tn.delta(a, b)


def _separate(m: Monomial) -> Monomial:
    """Equivalent monomial in which every slot and derivative label of every occurrence is fresh and linked to
    the rest by an explicit ``g`` / ``delta`` in the coefficient (no direct occurrence-occurrence contraction)."""
    coeff = m.coefficient
    occs = []
    for o in m.occurrences:
        slots = []
        for lab in o.slots:
            new = tn.fresh(_space_of(lab))
            coeff = coeff * _identity(lab, new)
            slots.append(new)
        ders = []
        for lab in o.derivatives:
            new = tn.fresh(tn.LORENTZ)
            coeff = coeff * tn.g(lab, new)
            ders.append(new)
        occs.append(Occurrence(o.field, tuple(slots), tuple(ders)))
    return Monomial(coeff, tuple(occs))


# ---------------------------------------------------------------------------------------------- Density
class Density:
    """A sum of monomials (a local field polynomial, possibly with free labels). Immutable."""

    def __init__(self, fields: FieldTable, expr: Expression | int = 0, *, _canonical: bool = False) -> None:
        self.fields = fields
        e = Expression.num(expr) if isinstance(expr, int) else expr
        self.expr: Expression = e if _canonical else _normalize(fields, e)

    # ------------------------------------------------------------------ construction and views
    @classmethod
    def from_monomials(cls, fields: FieldTable, monomials: Iterable[Monomial]):
        """Sum of monomials; checks that every occurrence has the right number of slots."""
        total = Expression.num(0)
        for m in monomials:
            for o in m.occurrences:
                if len(o.slots) != len(fields[o.field].slots):
                    raise ValueError(f"{o.field} needs {len(fields[o.field].slots)} slot labels, got {len(o.slots)}")
            total = total + _encode(m)
        return cls(fields, total)

    @classmethod
    def monomial(cls, fields: FieldTable, coefficient: Expression | int, *occurrences: Occurrence):
        """A single monomial."""
        c = Expression.num(coefficient) if isinstance(coefficient, int) else coefficient
        return cls.from_monomials(fields, [Monomial(c, tuple(occurrences))])

    def _new(self, expr: Expression, canonical: bool = False):
        return type(self)(self.fields, expr, _canonical=canonical)

    def monomials(self) -> tuple[Monomial, ...]:
        """The canonical monomials (one term each, canonical occurrence order and labels)."""
        return tuple(_decode(t) for t in _terms(self.expr))

    def __iter__(self):
        return iter(self.monomials())

    def __len__(self) -> int:
        return len(_terms(self.expr))

    def is_zero(self) -> bool:
        return len(self) == 0

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Density):
            return NotImplemented
        return (self.expr - other.expr).expand() == Expression.num(0)

    def __hash__(self) -> int:
        return hash(self.expr.to_canonical_string())

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.expr})"

    def free_labels(self) -> list[Expression]:
        """Free labels (each term must have the same set)."""
        out: dict[str, Expression] = {}
        for t in _terms(self.expr):
            for lab in tn.free_labels(t):
                out[lab.to_canonical_string()] = lab
        return list(out.values())

    def parity(self) -> int:
        """Grassmann parity (0 or 1); raises if the terms have different parities."""
        pars = {sum(self.fields.odd(o.field) for o in m.occurrences) % 2 for m in self.monomials()}
        if len(pars) > 1:
            raise ValueError("density of mixed parity")
        return pars.pop() if pars else 0

    def field_content(self) -> Counter:
        """``Counter`` of the field multisets (tuples in canonical order) of the monomials."""
        return Counter(m.fields() for m in self.monomials())

    # ------------------------------------------------------------------ algebra
    def __add__(self, other: Density) -> Density:
        self._check(other)
        return self._new(self.expr + other.expr, canonical=True)._renormalized()

    def __sub__(self, other: Density) -> Density:
        self._check(other)
        return self._new(self.expr - other.expr, canonical=True)._renormalized()

    def __neg__(self) -> Density:
        return self._new(-self.expr, canonical=True)._renormalized()

    def scale(self, c: Expression | int) -> Density:
        """Multiplication by a Grassmann-even scalar or tensor coefficient."""
        c = Expression.num(c) if isinstance(c, int) else c
        return self._new(self.expr * c)

    def _renormalized(self) -> Density:
        return self._new(self.expr.expand(), canonical=True)

    def _check(self, other: Density) -> None:
        if self.fields is not other.fields:
            raise ValueError("densities over different field tables")

    def __mul__(self, other: Density) -> Density:
        """Graded product ``self * other`` (occurrences of ``self`` to the left). Dummy labels of ``other`` are
        renamed apart; free labels are shared (and hence contracted if equal)."""
        self._check(other)
        right = _fresh_dummies(other.expr)
        total = Expression.num(0)
        left_monos = [_decode(t) for t in _terms(_fresh_dummies(self.expr))]
        right_monos = [_decode(t) for t in _terms(right)]
        for a in left_monos:
            for b in right_monos:
                total = total + _encode(Monomial(a.coefficient * b.coefficient, a.occurrences + b.occurrences))
        return Density(self.fields, total)

    def total_derivative(self, mu: Expression) -> Density:
        """``d_mu`` of the density (Leibniz rule over the occurrences)."""
        total = Expression.num(0)
        for m in self.monomials():
            for k, o in enumerate(m.occurrences):
                occs = list(m.occurrences)
                occs[k] = Occurrence(o.field, o.slots, o.derivatives + (mu,))
                total = total + _encode(Monomial(m.coefficient, tuple(occs)))
        return self._new(total)

    def rename(self, mapping: Sequence[tuple[Expression, Expression]]) -> Density:
        """Simultaneous renaming of (free) labels."""
        return self._new(tn.rename(self.expr, mapping))

    def apply(self, derivation: Derivation) -> Density:
        """Apply a graded derivation (acting from the left)."""
        return derivation(self)

    def map_coefficients(self, fn: Callable[[Expression], Expression]) -> Density:
        """Apply ``fn`` to every coefficient (e.g. a substitution of parameters)."""
        total = Expression.num(0)
        for m in self.monomials():
            total = total + _encode(Monomial(fn(m.coefficient), m.occurrences))
        return self._new(total)

    def conjugate(self) -> Density:
        """Hermitian conjugate: occurrences reversed and conjugated, coefficient by :func:`tensors.conjugate`."""
        total = Expression.num(0)
        for m in self.monomials():
            sign = 1
            occs = []
            for o in reversed(m.occurrences):
                info = self.fields[o.field]
                sign *= info.conj_sign
                occs.append(Occurrence(info.conjugate, o.slots, o.derivatives))
            total = total + sign * _encode(Monomial(tn.conjugate(m.coefficient), tuple(occs)))
        return self._new(total)

    def simplify(self, simplifier: tn.Simplifier) -> Density:
        """Apply the tensor simplifier to every monomial (metrics and deltas are absorbed into field slots and
        derivatives as well) and re-canonicalise."""
        return self._new(simplifier(self.expr))

    # ------------------------------------------------------------------ BMHV projections
    def split_lorentz(self, simplifier: tn.Simplifier) -> Density:
        """Write every Lorentz index of a field slot or derivative as the sum of its barred and hatted parts,
        ``X_mu = (gbar + ghat)_mu^nu X_nu``, and simplify. Afterwards every Lorentz label of a field or derivative
        is attached to a ``gbar`` or a ``ghat``."""
        total = Expression.num(0)
        for m in self.monomials():
            coeff = m.coefficient
            occs = []
            for o in m.occurrences:
                new_slots = []
                for lab in o.slots:
                    if tn.label_space(lab) == "lorentz":
                        nu = tn.fresh(tn.LORENTZ)
                        coeff = coeff * (tn.gbar(lab, nu) + tn.ghat(lab, nu))
                        new_slots.append(nu)
                    else:
                        new_slots.append(lab)
                new_ders = []
                for lab in o.derivatives:
                    nu = tn.fresh(tn.LORENTZ)
                    coeff = coeff * (tn.gbar(lab, nu) + tn.ghat(lab, nu))
                    new_ders.append(nu)
                occs.append(Occurrence(o.field, tuple(new_slots), tuple(new_ders)))
            total = total + _encode(Monomial(coeff, tuple(occs)))
        return self._new(simplifier(total))

    def physical(self, simplifier: tn.Simplifier) -> Density:
        """Physical projection ``P`` (report eq:Pproj): after :meth:`split_lorentz`, the terms without a hatted
        object (``ghat``, ``gammahat``, ``sphat``).
        Closed hatted traces have been evaluated (``ghat(mu,mu) = d - 4``) before the projection."""
        return self._project(simplifier, keep_hat=False)

    def evanescent(self, simplifier: tn.Simplifier) -> Density:
        """Evanescent part ``E = 1 - P``."""
        return self._project(simplifier, keep_hat=True)

    def four_dimensional_limit(self, simplifier: tn.Simplifier) -> Density:
        """``P`` followed by ``d -> 4``."""
        return self.physical(simplifier).map_coefficients(lambda c: c.replace(tn.d, Expression.num(4)))

    def _project(self, simplifier: tn.Simplifier, keep_hat: bool) -> Density:
        split = self.split_lorentz(simplifier)
        total = Expression.num(0)
        for t in _terms(split.expr):
            if tn.is_evanescent(t) == keep_hat:
                total = total + t
        return self._new(total, canonical=True)


# ---------------------------------------------------------------------------------------------- LocalFunctional
class LocalFunctional(Density):
    """``int d^d x`` of a density without free labels."""

    def __init__(self, fields: FieldTable, expr: Expression | int = 0, *, _canonical: bool = False) -> None:
        super().__init__(fields, expr, _canonical=_canonical)
        if self.free_labels():
            raise ValueError(f"a local functional has no free labels: {self.free_labels()}")

    @classmethod
    def integrate(cls, density: Density) -> LocalFunctional:
        """``int d^d x density``."""
        return cls(density.fields, density.expr, _canonical=True)

    def __mul__(self, other):  # the product of two integrals is not local
        raise TypeError("product of local functionals is not local; multiply densities instead")

    def density(self) -> Density:
        return Density(self.fields, self.expr, _canonical=True)

    # ------------------------------------------------------------------ derivatives
    def functional_derivative(self, field: str, labels: Sequence[Expression], side: str = "left") -> Density:
        """``dL F / d field(x)`` (``side="left"``) or ``dR F / d field(x)`` (``side="right"``); the result carries the
        free ``labels`` in the slots of ``field``. Euler-Lagrange: ``(-1)^r d_{mu_1..mu_r}`` acts on the rest."""
        if side not in ("left", "right"):
            raise ValueError("side must be 'left' or 'right'")
        info = self.fields[field]
        if len(labels) != len(info.slots):
            raise ValueError(f"{field} has {len(info.slots)} slots")
        result = Density(self.fields)
        expr = _fresh_dummies(self.expr)
        for term in _terms(expr):
            m = _separate(_decode(term))
            occs = m.occurrences
            for k, o in enumerate(occs):
                if o.field != field:
                    continue
                others = [self.fields.odd(x.field) for x in (occs[:k] if side == "left" else occs[k + 1:])]
                sign = -1 if (info.odd and sum(others) % 2) else 1
                rest = Monomial(m.coefficient, occs[:k] + occs[k + 1:])
                mapping = list(zip(o.slots, labels))
                piece = Density(self.fields, tn.absorb(tn.rename(_encode(rest), mapping)))
                for mu in o.derivatives:
                    piece = piece.total_derivative(mu)
                if len(o.derivatives) % 2:
                    sign = -sign
                result = result + piece.scale(sign)
        return result

    def integrate_by_parts(self, field: str) -> LocalFunctional:
        """Move all derivatives off the first occurrence of ``field`` in every monomial:
        ``int (d_mu X) R = - int X d_mu R``."""
        current = self
        while True:
            total = Expression.num(0)
            changed = False
            for m in current.monomials():
                k = next((i for i, o in enumerate(m.occurrences) if o.field == field), None)
                if k is None or not m.occurrences[k].derivatives:
                    total = total + _encode(m)
                    continue
                changed = True
                o = m.occurrences[k]
                mu, keep = o.derivatives[0], o.derivatives[1:]
                for j, x in enumerate(m.occurrences):
                    if j == k:
                        continue
                    occs = list(m.occurrences)
                    occs[k] = Occurrence(o.field, o.slots, keep)
                    occs[j] = Occurrence(x.field, x.slots, x.derivatives + (mu,))
                    total = total - _encode(Monomial(m.coefficient, tuple(occs)))
            current = LocalFunctional(self.fields, total)
            if not changed:
                return current

    # ------------------------------------------------------------------ vertex functions
    def vertex_function(self, legs: Sequence[Leg], simplifier: tn.Simplifier | None = None,
                        momentum_conservation: bool = True) -> Expression:
        """Vertex function ``V(legs)`` (module docstring; Feynman rule ``i V``).

        ``legs`` carry free labels and momentum numbers. With ``momentum_conservation`` the momentum of the last
        leg is eliminated, ``p_last = -sum(other p)``. With a ``simplifier`` the result is simplified and its
        dummies relabelled canonically, so that equal vertex functions compare equal.
        """
        want = sorted(l.field for l in legs)
        total = Expression.num(0)
        expr = _fresh_dummies(self.expr)
        for term in _terms(expr):
            m = _separate(_decode(term))
            if sorted(m.fields()) != want:
                continue
            occs = m.occurrences
            groups: dict[str, list[int]] = {}
            for j, l in enumerate(legs):
                groups.setdefault(l.field, []).append(j)
            occ_groups = {fn: [i for i, o in enumerate(occs) if o.field == fn] for fn in groups}
            names = list(groups)
            for choice in itertools.product(*(itertools.permutations(groups[n]) for n in names)):
                assign: dict[int, int] = {}  # occurrence -> leg
                for n, perm in zip(names, choice):
                    for i, j in zip(occ_groups[n], perm):
                        assign[i] = j
                odd_legs = [assign[i] for i in range(len(occs)) if self.fields.odd(occs[i].field)]
                sign = _perm_sign(odd_legs)
                mapping = []
                factor = Expression.num(sign)
                for i, o in enumerate(occs):
                    leg = legs[assign[i]]
                    mapping += list(zip(o.slots, leg.labels))
                    for mu in o.derivatives:
                        factor = factor * (-tn.I) * tn.p(leg.momentum, mu)
                total = total + tn.rename(m.coefficient, mapping) * factor
        if momentum_conservation and legs:
            last = legs[-1].momentum
            a = S("a_")
            others = Expression.num(0)
            for l in legs[:-1]:
                others = others + tn.p(l.momentum, a)
            total = total.replace(tn.p(last, a), -others)
        total = tn.absorb(total)
        if simplifier is None:
            return total
        return simplifier(total)

    def leg_classes(self) -> list[tuple[str, ...]]:
        """The field multisets of the monomials, each in canonical field order, sorted."""
        return sorted({m.fields() for m in self.monomials()}, key=lambda c: [self.fields.order(x) for x in c])

    def vertex_functions(self, simplifier: tn.Simplifier) -> dict[tuple[str, ...], Expression]:
        """``{leg class: V(standard_legs(class))}`` for every leg class of the functional (simplified)."""
        return {c: self.vertex_function(standard_legs(self.fields, c), simplifier) for c in self.leg_classes()}

    def vanishes(self, simplifier: tn.Simplifier) -> bool:
        """True if every vertex function vanishes after momentum conservation and simplification: the
        functional is zero up to integration by parts (and the identities the simplifier knows). A sufficient
        symbolic test only (``False`` is no proof of non-vanishing); the certified test, which needs the model, is
        :func:`hepkit_renorm.fingerprint.is_zero`."""
        return all(v.expand() == Expression.num(0) for v in self.vertex_functions(simplifier).values())


# ---------------------------------------------------------------------------------------------- derivations
@dataclass(frozen=True)
class FieldImage:
    """The image ``delta Phi`` of a field under a derivation: a density whose free labels are ``labels``, one per
    slot of ``Phi`` (they are replaced by the slot labels of each occurrence)."""

    labels: tuple[Expression, ...]
    density: Density


@dataclass(frozen=True)
class Derivation:
    """A graded derivation of parity ``parity`` acting from the left (module docstring), defined by its images on
    fields; fields without an image are annihilated. It commutes with partial derivatives."""

    parity: int
    images: Mapping[str, FieldImage] = dc_field(default_factory=dict)

    def __call__(self, d: Density) -> Density:
        fields = d.fields
        for name, img in self.images.items():
            if not img.density.is_zero() and img.density.parity() != (fields.odd(name) + self.parity) % 2:
                raise ValueError(f"image of {name} has the wrong parity")
        total = Density(fields)
        for term in _terms(_fresh_dummies(d.expr)):  # fresh dummies cannot collide with canonical image labels
            m = _decode(term)
            parity_left = 0
            for k, o in enumerate(m.occurrences):
                img = self.images.get(o.field)
                if img is not None and not img.density.is_zero():
                    sign = -1 if (self.parity and parity_left % 2) else 1
                    piece = Density(fields, _fresh_dummies(img.density.expr, keep=img.labels))
                    piece = piece.rename(list(zip(img.labels, o.slots)))
                    for mu in o.derivatives:
                        piece = piece.total_derivative(mu)
                    left = Monomial(m.coefficient, m.occurrences[:k])
                    right = Monomial(Expression.num(sign), m.occurrences[k + 1:])
                    total = total + _splice(fields, left, piece, right)
                parity_left += fields.odd(o.field)
        return type(d)(fields, total.expr, _canonical=True) if isinstance(d, LocalFunctional) else total


def _splice(fields: FieldTable, left: Monomial, middle: Density, right: Monomial) -> Density:
    total = Expression.num(0)
    for t in _terms(middle.expr):
        mm = _decode(t)
        total = total + _encode(Monomial(left.coefficient * mm.coefficient * right.coefficient,
                                         left.occurrences + mm.occurrences + right.occurrences))
    return Density(fields, total)
