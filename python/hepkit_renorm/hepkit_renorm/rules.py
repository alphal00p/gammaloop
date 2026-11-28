"""Feynman rules of a local functional: vertex rules, propagators, rule sets.

Everything here goes through :meth:`LocalFunctional.vertex_function`, the single code path for vertex functions
(report 03b sec:dictionary, 02_method sec:m-objects). This module adds the conventions a user of the rules needs.

Vertex rules (:class:`VertexRule`, :func:`vertex_rule`, :func:`tree_rules`)
---------------------------------------------------------------------------
* The vertex function ``V(Phi_1(p_1), ..., Phi_n(p_n))`` of a functional ``Y`` is the coefficient of ``theta_1 ...
  theta_n`` in ``Y[sum_k theta_k e_k exp(-i p_k.x)]`` with the polarisations ``e_k`` stripped: the free labels of
  leg ``k`` stand where its polarisation was. ``theta_k`` is Grassmann odd for an odd field, so ``V`` depends on the
  *leg order*: exchanging two odd legs flips its sign (:meth:`VertexRule.reordered`). The Feynman rule is ``i V``.
  A matrix-valued rule (legs ``Psibar`` and ``Psi``, or a source of either) carries the spinor and group labels of
  both legs; its row index is the label of the barred leg (``Psibar``, ``R_f``), as in ``Psibar M Psi``.
* Momenta: leg ``k`` has the momentum ``p(k, mu)``, incoming unless ``k`` is listed in ``outgoing``. For an outgoing
  leg ``p(k, mu)`` is minus the incoming momentum; this is the report's ``p'`` of an outgoing ``Psibar`` or
  ``H^dagger`` (eq:V-ssv, eq:V-delta). Momentum conservation eliminates the momentum of one leg (``eliminate``,
  the last by default, ``None`` for no elimination): ``sum_k sigma_k p(k) = 0`` with ``sigma = -1`` for outgoing legs.
* Canonical leg order (:func:`tree_rules`): the field order of the :class:`~hepkit_renorm.fields.FieldTable`, so
  that every rule of a rule set has a fixed Grassmann sign. The report writes some rules in a different order (for
  instance ``omegabar A omega`` for the ghost vertex, ``omega Psibar R_f`` for a source vertex): use
  :meth:`VertexRule.reordered` to obtain them; the sign is the sign of the permutation of the odd legs.
* Labels of the legs are the deterministic labels of :func:`~hepkit_renorm.functional.standard_legs`: rules of
  the same class built twice have identical labels, so two rules compare by the difference of their values.

Which leg classes are vertices: every class with at least three legs and every class with a BRST source (sources
never propagate, so ``rho omega`` is a two-point vertex). The other two-point classes form the quadratic form.

Propagators (:class:`Propagator`, :func:`propagators`)
------------------------------------------------------
For fields ``Phi_i`` with two-point vertex functions let ``K_ij(k) = V(Phi_i(-k), Phi_j(k))`` (leg ``i`` incoming
momentum ``-k``, leg ``j`` incoming momentum ``k``; ``k = p(0, mu)``). The propagator is ``D = i K^{-1}``:

* bosons (real fields, a complex scalar with its conjugate, gauge fields mixed with ``b``): ``sum_j K_ij D_jl = i
  delta_il``; ``D_jl(k)`` joins a ``Phi_j`` leg of one vertex, into which ``k`` flows, to a ``Phi_l`` leg of another.
  For a complex scalar ``D_{H Hdag} = i / K_{Hdag H}`` with ``k`` along the arrow (report eq:propagators).
* Grassmann-odd pairs ``(X, Xbar)`` (fermion and its Dirac conjugate, ghost and antighost):
  ``D_{X Xbar}(k) = i K_{Xbar X}(k)^{-1}``, the inverse of the matrix ``V(Xbar(-k), X(k))``, with ``k`` along the
  arrow (from the ``Xbar`` end to the ``X`` end); its first spinor index is the one of ``X``. With the antighost
  listed first this gives ``+i/k^2`` for ``omega omegabar`` (eq:propagators); the opposite order would flip it.

The inversion is exact. Bosonic quadratic forms must be covariant (``a g^{mu nu} + b k^mu k^nu`` between vector
fields, ``c k^mu`` between a vector and a scalar); spinor forms may contain ``gammabar.k``, ``gammahat.k`` and
``gamma5`` (the algebra they generate is inverted exactly, so a Dirac field with its sterile partner, whose kinetic
operator is the full ``kslash``, and its four-dimensional part ``kslashbar`` are handled alike). The internal
(group, flavour) structure must be a product of Kronecker deltas. Anything else raises ``NotImplementedError``; a
singular form (e.g. a gauge field without gauge fixing) raises ``ValueError``.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

from symbolica import AtomType, Expression, S

from . import tensors as tn
from .exact import PrimeField, default_field
from .fields import FieldTable
from .functional import Leg, LocalFunctional, standard_legs
from .spec import Model

__all__ = [
    "PROPAGATOR_MOMENTUM",
    "Propagator",
    "VertexRule",
    "is_vertex_class",
    "propagators",
    "quadratic_classes",
    "rule_table",
    "tree_rules",
    "vertex_rule",
]

PROPAGATOR_MOMENTUM = 0
"""Momentum number of propagators: ``k = p(0, mu)``."""

_a = S("a_")


def _terms(expr: Expression) -> list[Expression]:
    if expr.get_type() == AtomType.Add:
        return list(expr)
    if expr.get_type() == AtomType.Num and expr.to_canonical_string() == "0":
        return []
    return [expr]


def _factors(term: Expression) -> list[Expression]:
    return list(term) if term.get_type() == AtomType.Mul else [term]


def _is_zero(expr: Expression) -> bool:
    return expr.expand() == Expression.num(0)


def _perm_sign(seq: Sequence[int]) -> int:
    inv = sum(1 for i, j in itertools.combinations(range(len(seq)), 2) if seq[i] > seq[j])
    return -1 if inv % 2 else 1


# ---------------------------------------------------------------------------------------------- vertex rules
@dataclass(frozen=True)
class VertexRule:
    """The vertex function ``value = V(legs)`` of a functional (module docstring); Feynman rule ``i V``.

    ``parities`` are the Grassmann parities of the legs, ``outgoing`` the leg numbers whose ``p(k, mu)`` is the
    outgoing momentum, ``eliminated`` the leg whose momentum was removed by momentum conservation (``None``: none).
    """

    legs: tuple[Leg, ...]
    value: Expression
    parities: tuple[int, ...]
    outgoing: frozenset[int] = field(default_factory=frozenset)
    eliminated: int | None = None

    @property
    def fields(self) -> tuple[str, ...]:
        return tuple(leg.field for leg in self.legs)

    @property
    def feynman_rule(self) -> Expression:
        """``i V``."""
        return (tn.I * self.value).expand()

    def labels(self) -> list[tuple[Expression, ...]]:
        """The free labels of the legs, in leg order."""
        return [leg.labels for leg in self.legs]

    def reordered(self, order: Sequence[int]) -> VertexRule:
        """The same vertex with the legs in the order ``[legs[i] for i in order]``: labels and momenta stay
        attached to their legs, the value changes by the sign of the permutation of the odd legs."""
        if sorted(order) != list(range(len(self.legs))):
            raise ValueError("order must be a permutation of the leg positions")
        sign = _perm_sign([i for i in order if self.parities[i]])
        return VertexRule(tuple(self.legs[i] for i in order), (self.value * sign).expand(),
                          tuple(self.parities[i] for i in order), self.outgoing, self.eliminated)

    def to_dict(self) -> dict:
        """JSON-ready canonical form: legs (field, labels, momentum, outgoing), the value ``V`` as a canonical
        Symbolica string, and the conventions."""
        return {
            "legs": [{"field": leg.field, "labels": [lab.to_canonical_string() for lab in leg.labels],
                      "momentum": leg.momentum, "outgoing": leg.momentum in self.outgoing} for leg in self.legs],
            "vertex_function": self.value.to_canonical_string(),
            "convention": "V (Feynman rule i V); Grassmann order = leg order; p(k, mu) incoming unless outgoing",
            "eliminated_momentum": self.eliminated,
        }


def vertex_rule(functional: LocalFunctional, fields: Sequence[str], simplifier: tn.Simplifier | None = None, *,
                outgoing: Iterable[int] = (), eliminate: int | None = -1) -> VertexRule:
    """The vertex function of ``functional`` for the legs ``fields`` (in this order, momenta ``1..n``, labels of
    :func:`standard_legs`).

    ``outgoing``: leg numbers (1-based) whose ``p(k, mu)`` denotes the outgoing momentum. ``eliminate``: leg number
    whose momentum is eliminated by momentum conservation (``-1``: the last leg, ``None``: no elimination).
    The value is simplified with ``simplifier`` (canonical mode for comparisons, readable mode for display).
    """
    table = functional.fields
    legs = standard_legs(table, fields)
    n = len(legs)
    out = frozenset(outgoing)
    if not out <= set(range(1, n + 1)):
        raise ValueError("outgoing legs must be leg numbers 1..n")
    raw = functional.vertex_function(legs, None, momentum_conservation=False)
    for k in out:
        raw = raw.replace(tn.p(k, _a), -tn.p(k, _a))
    elim = None
    if eliminate is not None and n:
        elim = n if eliminate == -1 else eliminate
        if not 1 <= elim <= n:
            raise ValueError("eliminate must be a leg number")
        sigma = {k: (-1 if k in out else 1) for k in range(1, n + 1)}
        rest = Expression.num(0)
        for k in range(1, n + 1):
            if k != elim:
                rest = rest + sigma[k] * tn.p(k, _a)
        raw = raw.replace(tn.p(elim, _a), -sigma[elim] * rest)
    value = tn.absorb(raw)
    if simplifier is not None:
        value = simplifier(value)
    return VertexRule(tuple(legs), value, tuple(int(table.odd(f)) for f in fields), out, elim)


def is_vertex_class(fields: FieldTable, cls: Sequence[str]) -> bool:
    """A leg class is a vertex if it has at least three legs or contains a BRST source."""
    return len(cls) >= 3 or any(fields[f].is_source for f in cls)


def quadratic_classes(functional: LocalFunctional) -> list[tuple[str, ...]]:
    """The two-point classes without sources (the quadratic form that is inverted into propagators)."""
    return [c for c in functional.leg_classes() if not is_vertex_class(functional.fields, c)]


def _group_open_labels(model: Model, legs: Sequence[Leg]) -> list[Expression]:
    return [lab for leg in legs for lab, sl in zip(leg.labels, model.fields[leg.field].slots) if sl.factor is not None]


def tree_rules(functional: LocalFunctional, model: Model, simplifier: tn.Simplifier | None = None,
               field_: PrimeField | None = None) -> dict[tuple[str, ...], VertexRule]:
    """All non-vanishing vertex rules of ``functional``, one per vertex class (:func:`is_vertex_class`), legs in
    canonical field order, last momentum eliminated, values in the canonical form of ``simplifier`` (default
    ``model.simplifier()``).

    A class whose vertex function vanishes (after the simplifier, or after explicit evaluation of its group tensors
    over ``field_``, :meth:`Model.vanishes_modulo_group`) is dropped: e.g. the four-point classes of ``s S_0`` that
    cancel by ``[T^a, T^b] = i f^{abc} T^c``.
    """
    simp = simplifier or model.simplifier()
    fp = field_ or default_field()
    out: dict[tuple[str, ...], VertexRule] = {}
    for cls in functional.leg_classes():
        if not is_vertex_class(functional.fields, cls):
            continue
        rule = vertex_rule(functional, cls, simp)
        if _is_zero(rule.value) or model.vanishes_modulo_group(rule.value, _group_open_labels(model, rule.legs), fp):
            continue
        out[cls] = rule
    return out


# ---------------------------------------------------------------------------------------------- exact inversion
def _cancel(e: Expression) -> Expression:
    return e.together()


def _invert(matrix: list[list[Expression]]) -> list[list[Expression]]:
    """Exact inverse of a square matrix of rational functions (Gauss-Jordan; pivots tested with ``together``)."""
    n = len(matrix)
    a = [[_cancel(x) for x in row] + [Expression.num(int(i == j)) for j in range(n)] for i, row in enumerate(matrix)]
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c] != Expression.num(0)), None)
        if piv is None:
            raise ValueError("the quadratic form is singular (a gauge field without gauge fixing?)")
        a[c], a[piv] = a[piv], a[c]
        inv = _cancel(1 / a[c][c])
        a[c] = [_cancel(x * inv) for x in a[c]]
        for r in range(n):
            if r != c and a[r][c] != Expression.num(0):
                f = a[r][c]
                a[r] = [_cancel(x - f * y) for x, y in zip(a[r], a[c])]
    return [row[n:] for row in a]


# ---------------------------------------------------------------------------------------------- spinor algebra
# Elements of the algebra generated by bbar = gammabar.k, bhat = gammahat.k and gamma5: dicts {(x, y, z): coeff}
# for the word bbar^x bhat^y gamma5^z (x, y, z in {0, 1}); bbar^2 = kbar^2, bhat^2 = khat^2, bhat bbar = -bbar bhat,
# gamma5 bbar = -bbar gamma5, gamma5 bhat = bhat gamma5 (report eq:clifford, eq:g5comm).
_Word = tuple[int, int, int]
_WORDS: tuple[_Word, ...] = tuple(itertools.product((0, 1), repeat=3))


class _SpinorAlgebra:
    def __init__(self, kbar2: Expression, khat2: Expression) -> None:
        self.kbar2, self.khat2 = kbar2, khat2

    def mul_words(self, u: _Word, v: _Word) -> tuple[Expression, _Word]:
        x1, y1, z1 = u
        x2, y2, z2 = v
        sign = (-1) ** (z1 * x2 + y1 * x2)  # gamma5 past bbar, bhat past bbar
        c = Expression.num(sign)
        if x1 and x2:
            c = c * self.kbar2
        if y1 and y2:
            c = c * self.khat2
        return c, ((x1 + x2) % 2, (y1 + y2) % 2, (z1 + z2) % 2)

    def mul(self, a: Mapping[_Word, Expression], b: Mapping[_Word, Expression]) -> dict[_Word, Expression]:
        out: dict[_Word, Expression] = {}
        for u, cu in a.items():
            for v, cv in b.items():
                c, w = self.mul_words(u, v)
                out[w] = out.get(w, Expression.num(0)) + cu * cv * c
        return {w: c.expand() for w, c in out.items()}

    def inverse(self, a: Mapping[_Word, Expression]) -> dict[_Word, Expression]:
        """Right multiplication by ``a`` as an 8x8 matrix, inverted exactly; returns ``a^{-1}``."""
        cols = []
        for u in _WORDS:  # column u: image of the basis word u under x -> x a
            img = self.mul({u: Expression.num(1)}, a)
            cols.append([img.get(w, Expression.num(0)) for w in _WORDS])
        matrix = [[cols[j][i] for j in range(8)] for i in range(8)]
        inv = _invert(matrix)
        unit = _WORDS.index((0, 0, 0))
        return {w: inv[i][unit] for i, w in enumerate(_WORDS) if inv[i][unit] != Expression.num(0)}


# ---------------------------------------------------------------------------------------------- propagators
@dataclass(frozen=True)
class Propagator:
    """``value = D_{XY}(k)`` (module docstring) for ``fields = (X, Y)``; ``labels`` are the free labels of the two
    ends (``X`` first) and ``k = p(PROPAGATOR_MOMENTUM, mu)``."""

    fields: tuple[str, str]
    labels: tuple[tuple[Expression, ...], tuple[Expression, ...]]
    value: Expression


def _two_point(functional: LocalFunctional, a: str, b: str) -> tuple[list[Leg], Expression]:
    """``K_ab(k) = V(a(-k), b(k))`` (unsimplified) with ``k = p(0, mu)``."""
    legs = standard_legs(functional.fields, (a, b))
    raw = functional.vertex_function(legs, None, momentum_conservation=False)
    k = PROPAGATOR_MOMENTUM
    raw = raw.replace(tn.p(1, _a), -tn.p(k, _a)).replace(tn.p(2, _a), tn.p(k, _a))
    return legs, tn.absorb(raw)


def _internal(table: FieldTable, name: str) -> list[int]:
    return [i for i, s in enumerate(table[name].slots) if s.space.kind not in ("lorentz", "spinor")]


def _strip_internal(table: FieldTable, legs: Sequence[Leg], value: Expression, readable: tn.Simplifier,
                    ) -> tuple[Expression, Expression]:
    """Write ``value = R * prod delta(internal labels of leg 1, of leg 2)``; returns ``(R, prod delta)``."""
    i1, i2 = _internal(table, legs[0].field), _internal(table, legs[1].field)
    if [table[legs[0].field].slots[i].space for i in i1] != [table[legs[1].field].slots[i].space for i in i2]:
        raise NotImplementedError(f"two-point form {legs[0].field} {legs[1].field}: different internal spaces")
    deltas, norm = Expression.num(1), 1
    for j1, j2 in zip(i1, i2):
        deltas = deltas * tn.delta(legs[0].labels[j1], legs[1].labels[j2])
        norm *= table[legs[0].field].slots[j1].space.dim
    reduced = readable(value * deltas / norm)
    if not _is_zero(readable(reduced * deltas - value)):
        raise NotImplementedError(f"two-point form {legs[0].field} {legs[1].field} is not diagonal in its internal "
                                  "indices")
    return reduced, deltas


def _components(functional: LocalFunctional) -> list[list[str]]:
    """Connected components of the graph of the quadratic form, each in field order."""
    table = functional.fields
    adj: dict[str, set[str]] = {}
    for c in quadratic_classes(functional):
        x, y = c
        adj.setdefault(x, set()).add(y)
        adj.setdefault(y, set()).add(x)
    seen: set[str] = set()
    comps = []
    for start in sorted(adj, key=table.order):
        if start in seen:
            continue
        stack, comp = [start], set()
        while stack:
            v = stack.pop()
            if v not in comp:
                comp.add(v)
                stack.extend(adj[v] - comp)
        seen |= comp
        comps.append(sorted(comp, key=table.order))
    return comps


_ODD_PRIMARY = ("fermion", "ghost")
_ODD_PARTNER = ("fermion_bar", "antighost")


def propagators(functional: LocalFunctional, model: Model, simplifier: tn.Simplifier | None = None,
                ) -> list[Propagator]:
    """The propagators ``D = i K^{-1}`` of the quadratic part of ``functional`` (module docstring), every
    non-vanishing entry once, simplified with ``simplifier`` (default ``model.simplifier()``)."""
    table = functional.fields
    simp = simplifier or model.simplifier()
    readable = model.simplifier(expand_projectors=False)
    out: list[Propagator] = []
    for comp in _components(functional):
        if any(table.odd(f) for f in comp):
            out.extend(_odd_propagator(functional, comp, readable, simp))
        else:
            out.extend(_boson_propagators(functional, comp, readable, simp))
    return out


def _odd_propagator(functional, comp, readable, simp) -> list[Propagator]:
    table = functional.fields
    prim = [f for f in comp if table[f].kind in _ODD_PRIMARY]
    part = [f for f in comp if table[f].kind in _ODD_PARTNER]
    if len(comp) != 2 or len(prim) != 1 or len(part) != 1:
        raise NotImplementedError(f"odd quadratic form among {comp}: only (field, conjugate) pairs are supported")
    x, xbar = prim[0], part[0]
    legs, value = _two_point(functional, xbar, x)
    reduced, deltas = _strip_internal(table, legs, value, readable)
    k = PROPAGATOR_MOMENTUM
    spin_bar = [lab for lab, s in zip(legs[0].labels, table[xbar].slots) if s.space.kind == "spinor"]
    spin = [lab for lab, s in zip(legs[1].labels, table[x].slots) if s.space.kind == "spinor"]
    if spin:
        algebra = _SpinorAlgebra(tn.spbar(k, k), tn.sphat(k, k))
        element = _spinor_element(readable(reduced), spin_bar[0], spin[0], k)
        inv = algebra.inverse(element)
        d = tn.I * _spinor_expression(inv, spin[0], spin_bar[0], k)
    else:
        if tn.free_labels(reduced):
            raise NotImplementedError(f"two-point form {xbar} {x} has free Lorentz labels")
        d = tn.I / reduced
    return [Propagator((x, xbar), (legs[1].labels, legs[0].labels), simp(_cancel(d) * deltas))]


def _chain_factor(fa: Expression, used: dict) -> tuple[Expression, Expression, dict[_Word, Expression]]:
    """One spinor matrix factor -> (row label, column label, algebra element)."""
    one = Expression.num(1)
    half = one / 2
    name = fa.get_name()
    args = list(fa)
    if name == tn.delta.get_name():
        return args[0], args[1], {(0, 0, 0): one}
    if name == tn.gamma5.get_name():
        return args[0], args[1], {(0, 0, 1): one}
    if name == tn.PL.get_name():
        return args[0], args[1], {(0, 0, 0): half, (0, 0, 1): -half}
    if name == tn.PR.get_name():
        return args[0], args[1], {(0, 0, 0): half, (0, 0, 1): half}
    mu = args[0].to_canonical_string()
    if mu not in used:
        raise NotImplementedError(f"gamma matrix {fa} not contracted with the momentum")
    used[mu] = True
    if name == tn.gamma.get_name():
        return args[1], args[2], {(1, 0, 0): one, (0, 1, 0): one}
    if name == tn.gammabar.get_name():
        return args[1], args[2], {(1, 0, 0): one}
    if name == tn.gammahat.get_name():
        return args[1], args[2], {(0, 1, 0): one}
    raise NotImplementedError(f"unsupported spinor factor {fa}")


def _spinor_element(expr: Expression, row: Expression, col: Expression, k: int) -> dict[_Word, Expression]:
    """Decompose a spinor matrix ``M_{row col}`` (readable form, momentum ``p(k)``) into the algebra of
    :class:`_SpinorAlgebra`; scalar products are kept as coefficients."""
    algebra = _SpinorAlgebra(tn.spbar(k, k), tn.sphat(k, k))
    spinor_heads = {h.get_name() for h in (tn.delta, tn.gamma5, tn.PL, tn.PR, tn.gamma, tn.gammabar, tn.gammahat)}
    total: dict[_Word, Expression] = {}
    for term in _terms(expr):
        scalar, chain, moms = Expression.num(1), [], {}
        for fa in _factors(term):
            if fa.get_type() == AtomType.Fn and fa.get_name() in spinor_heads and any(
                    tn.label_space(x) == "spinor" for x in list(fa) if x.get_type() == AtomType.Fn):
                chain.append(fa)
            elif fa.get_type() == AtomType.Fn and fa.get_name() == tn.p.get_name():
                moms[list(fa)[1].to_canonical_string()] = False
            else:
                scalar = scalar * fa
        used = dict(moms)
        pieces = {}
        for fa in chain:
            r, c, el = _chain_factor(fa, used)
            pieces[r.to_canonical_string()] = (c, el)
        if not all(used.values()):
            raise NotImplementedError(f"momentum not contracted with a gamma matrix in {term}")
        element: dict[_Word, Expression] = {(0, 0, 0): Expression.num(1)}
        cur = row.to_canonical_string()
        for _ in range(len(pieces)):
            if cur not in pieces:
                raise NotImplementedError(f"broken spinor chain in {term}")
            c, el = pieces.pop(cur)
            element = algebra.mul(element, el)
            cur = c.to_canonical_string()
        if pieces or cur != col.to_canonical_string():
            raise NotImplementedError(f"spinor chain does not run from the barred to the unbarred leg in {term}")
        for w, cw in element.items():
            total[w] = (total.get(w, Expression.num(0)) + scalar * cw).expand()
    return {w: c for w, c in total.items() if not _is_zero(c)}


def _spinor_expression(element: Mapping[_Word, Expression], row: Expression, col: Expression, k: int) -> Expression:
    """Inverse of :func:`_spinor_element`: the matrix ``sum c_w w`` with labels ``(row, col)``."""
    out = Expression.num(0)
    for (x, y, z), c in element.items():
        factors = []
        if x:
            mu = tn.fresh(tn.LORENTZ)
            factors.append(lambda s, t, mu=mu: tn.gammabar(mu, s, t) * tn.p(k, mu))
        if y:
            mu = tn.fresh(tn.LORENTZ)
            factors.append(lambda s, t, mu=mu: tn.gammahat(mu, s, t) * tn.p(k, mu))
        if z:
            factors.append(lambda s, t: tn.gamma5(s, t))
        if not factors:
            out = out + c * tn.delta(row, col)
            continue
        labels = [row] + [tn.fresh(tn.SPINOR) for _ in factors[:-1]] + [col]
        term = c
        for f, (s, t) in zip(factors, itertools.pairwise(labels)):
            term = term * f(s, t)
        out = out + term
    return out


def _lorentz_label(table: FieldTable, leg: Leg) -> Expression | None:
    labs = [lab for lab, s in zip(leg.labels, table[leg.field].slots) if s.space.kind == "lorentz"]
    if len(labs) > 1:
        raise NotImplementedError(f"{leg.field}: bosonic fields with more than one Lorentz index")
    return labs[0] if labs else None


def _boson_propagators(functional, comp, readable, simp) -> list[Propagator]:
    """Real bosons, complex scalars with their conjugates and vector-scalar mixing (module docstring)."""
    table = functional.fields
    k = PROPAGATOR_MOMENTUM
    k2 = tn.sp(k, k)
    vec = [f for f in comp if any(s.space.kind == "lorentz" for s in table[f].slots)]
    sca = [f for f in comp if f not in vec]
    n_v = len(vec)
    order = vec + sca
    size = len(order)
    a_t = [[Expression.num(0)] * n_v for _ in range(n_v)]   # transverse block (vectors)
    long = [[Expression.num(0)] * size for _ in range(size)]  # longitudinal / scalar block
    present = {tuple(sorted(c, key=table.order)) for c in quadratic_classes(functional)}
    gs, pm, pn = S("hepkit_renorm_rules::G", "hepkit_renorm_rules::Pm", "hepkit_renorm_rules::Pn")
    for i, fi in enumerate(order):
        for j, fj in enumerate(order):
            if tuple(sorted((fi, fj), key=table.order)) not in present:
                continue
            legs, value = _two_point(functional, fi, fj)
            reduced, _ = _strip_internal(table, legs, value, readable)
            reduced = readable(reduced)
            mu, nu = _lorentz_label(table, legs[0]), _lorentz_label(table, legs[1])
            if mu is not None and nu is not None:
                r = reduced.replace(tn.g(mu, nu), gs).replace(tn.p(k, mu), pm).replace(tn.p(k, nu), pn)
                a = r.derivative(gs).replace(gs, 0).replace(pm, 0).replace(pn, 0)
                b = r.derivative(pm).derivative(pn).replace(gs, 0).replace(pm, 0).replace(pn, 0)
                rest = (r - a * gs - b * pm * pn).expand()
                a_t[i][j] = a
                long[i][j] = (a + b * k2).expand()
            elif mu is not None or nu is not None:
                lab = mu if mu is not None else nu
                r = reduced.replace(tn.p(k, lab), pm)
                c = r.derivative(pm).replace(pm, 0)
                rest = (r - c * pm).expand()
                # vector-scalar entries: K_vs = c k^mu enters M as c, K_sv = c' k^nu as c' k^2 (comment below)
                long[i][j] = (c * k2).expand() if mu is None else c
            else:
                rest = Expression.num(0)
                long[i][j] = reduced
            if not _is_zero(rest) or (mu is None and nu is None and tn.free_labels(reduced)):
                raise NotImplementedError(f"two-point form {fi} {fj} is not of the covariant form: {reduced}")
    # sum_j K_ij D_jl = i delta_il with K = a T + (a + b k^2) L on vectors, D = x T + y L + z k^mu + z' k^nu + w:
    # x = i a^{-1}; the matrix M = [[a + b k^2, c], [c' k^2, e]] gives U = i M^{-1} = [[y, z], [z' k^2, w]].
    x = _invert(a_t) if n_v else []
    u = _invert(long)
    out = []
    for j, fj in enumerate(order):
        for l_, fl in enumerate(order):
            legs = standard_legs(table, (fj, fl))
            mu, nu = _lorentz_label(table, legs[0]), _lorentz_label(table, legs[1])
            ujl = u[j][l_]
            if mu is not None and nu is not None:
                transverse = tn.g(mu, nu) - tn.p(k, mu) * tn.p(k, nu) / k2
                longitudinal = tn.p(k, mu) * tn.p(k, nu) / k2
                d = tn.I * (x[j][l_] * transverse + ujl * longitudinal)
            elif mu is not None:
                d = tn.I * ujl * tn.p(k, mu)
            elif nu is not None:
                d = tn.I * ujl / k2 * tn.p(k, nu)
            else:
                d = tn.I * ujl
            d = _cancel(d.expand())
            if _is_zero(d):
                continue
            deltas = _internal_deltas(table, legs)
            out.append(Propagator((fj, fl), (legs[0].labels, legs[1].labels), simp(d * deltas)))
    return out


def _internal_deltas(table: FieldTable, legs: Sequence[Leg]) -> Expression:
    out = Expression.num(1)
    for j1, j2 in zip(_internal(table, legs[0].field), _internal(table, legs[1].field)):
        out = out * tn.delta(legs[0].labels[j1], legs[1].labels[j2])
    return out


def rule_table(rules: Mapping[tuple[str, ...], VertexRule]) -> list[dict]:
    """The rules as a list of :meth:`VertexRule.to_dict`, in the order of the mapping."""
    return [r.to_dict() for r in rules.values()]

