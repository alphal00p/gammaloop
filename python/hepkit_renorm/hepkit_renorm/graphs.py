"""Graph records: the input of the loop back-ends (:mod:`hepkit_renorm.loops`).

A :class:`GraphRecord` is one Feynman graph of an amputated (1PI) vertex function ``Gamma_{Phi_1 ... Phi_n}``, in the
conventions of the report (02_method sec:m-objects, eq:m-Gamma): ``i Gamma = sum_G G`` and the value of ``G``
includes the Feynman rules ``i V`` of its vertices, the propagators ``D``, the symmetry factor and the Grassmann
sign (closed fermion and ghost loops included). The record holds data only; the loop integration is done by a
back-end. Nothing in this module depends on the loop order.

Contents of a record
--------------------
* ``legs``: the external legs of the vertex function, in the order of the Grassmann parameters ``theta_1 ...
  theta_n`` (the leg order of the class). Their labels are the free labels of the result and their momenta
  ``p(k, mu)`` are incoming; with :func:`~hepkit_renorm.functional.standard_legs` they are exactly the legs of the
  :class:`~hepkit_renorm.fingerprint.Fingerprinter` samples, so the value of the graph is evaluated on the same
  samples as the vertex functions of local operators.
* ``vertices``: each :class:`GraphVertex` carries a vertex function ``V`` (any expression in the vocabulary of
  :mod:`hepkit_renorm.tensors`, e.g. :class:`~hepkit_renorm.rules.VertexRule` values of ``S_0``, of ``Delta_1`` or of
  counterterms) with its own legs (field, labels, momentum number), and for every leg the :data:`End` it is joined
  to: an external leg (:class:`ExternalEnd`) or one end of a propagator (:class:`PropagatorEnd`). A leg number in
  ``outgoing`` denotes the outgoing momentum (as in :class:`~hepkit_renorm.rules.VertexRule`). ``dimension`` is the
  mass dimension of the local operator the vertex comes from (4 for ``S_0``; ``Delta_1 = s S_0`` has dimension 5,
  since ``s`` raises the dimension by one), ``order`` its loop order (0 for ``S_0``, ``L`` for a counterterm of
  order ``L``); both enter power counting and the order bookkeeping of higher loops.
* ``propagators``: each :class:`GraphPropagator` joins two vertex legs. ``fields = (X, Y)``: the ``X`` end is
  joined to a leg of field ``X``, the ``Y`` end to a leg of field ``Y`` (for a Dirac field ``(Psi, Psibar)``, for a
  ghost ``(omega, omegabar)``, for a complex scalar ``(H, H^dagger)``), as in :class:`~hepkit_renorm.rules.Propagator`.
  The momentum ``k`` of the propagator flows into the vertex at the ``X`` end (for an odd field: along the arrow,
  from the ``Y = Xbar`` end to the ``X`` end). The value is ``numerator / prod_j (k^2 - m_j^2)^{n_j}`` with the
  numerator in ``k = p(0, mu)`` (:data:`~hepkit_renorm.rules.PROPAGATOR_MOMENTUM`) and the full ``d``-dimensional
  ``k^2`` (never split) in the denominators.

Derived, not typed in
---------------------
* :meth:`GraphRecord.routing`: the momentum of every propagator as a combination of loop momenta
  ``p(LOOP_MOMENTUM + l, mu)`` and external momenta, from momentum conservation at the vertices (spanning tree; the
  chords carry the loop momenta, in propagator order).
* :meth:`GraphRecord.grassmann_sign`: the Wick sign. Let ``S`` be the concatenation, over the vertices in record
  order, of the Grassmann-odd legs of each vertex in its leg order (the order of the ``theta`` of its vertex
  function); every external odd leg is the parameter ``theta_j`` of its class leg ``j``, every internal odd leg the
  token of its propagator end. The sign is the parity of the permutation that brings ``S`` to ``(X_1 Y_1) (X_2 Y_2)
  ... theta_1 theta_2 ...``: each odd propagator as its pair (``X`` end, ``Y`` end), then the external parameters in
  class order. A closed loop of odd propagators gets its ``-1`` from this rule. The product of the vertex
  functionals is taken in record order, which matters only when more than one vertex is Grassmann odd.
* :meth:`GraphRecord.superficial_degree`: ``omega(G) = 4 L + sum_v (dim_v - sum_legs dim) + sum_e (dim X_e + dim Y_e
  - 4)`` (mass dimension of the graph at ``d = 4``; for vertices of dimension 4 it is ``4 - sum_ext dim``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction

from symbolica import AtomType, Expression, S

from . import tensors as tn
from .fields import FieldTable
from .functional import Leg
from .rules import PROPAGATOR_MOMENTUM, Propagator, VertexRule

__all__ = [
    "LOOP_MOMENTUM",
    "End",
    "ExternalEnd",
    "GraphPropagator",
    "GraphRecord",
    "GraphVertex",
    "Momentum",
    "PropagatorEnd",
]

LOOP_MOMENTUM = 900
"""Momentum number of the loop momenta: loop ``l`` is ``p(LOOP_MOMENTUM + l, mu)``."""

Momentum = Mapping[int, Fraction]
"""A momentum as ``{momentum number: coefficient}`` (external legs ``leg.momentum``, loops ``LOOP_MOMENTUM + l``)."""


@dataclass(frozen=True)
class ExternalEnd:
    """A vertex leg joined to the external leg ``leg`` (index into :attr:`GraphRecord.legs`)."""

    leg: int


@dataclass(frozen=True)
class PropagatorEnd:
    """A vertex leg joined to end ``side`` (0: the ``X`` end, into which the momentum flows; 1: the ``Y`` end) of
    propagator ``propagator`` (index into :attr:`GraphRecord.propagators`)."""

    propagator: int
    side: int


End = ExternalEnd | PropagatorEnd


@dataclass(frozen=True)
class GraphVertex:
    """A vertex of a graph: vertex function ``value`` (Feynman rule ``i value``) on ``legs``, joined to ``ends``."""

    legs: tuple[Leg, ...]
    value: Expression
    ends: tuple[End, ...]
    outgoing: frozenset[int] = frozenset()
    dimension: Fraction = Fraction(4)
    order: int = 0
    origin: str = ""

    def __post_init__(self) -> None:
        if len(self.legs) != len(self.ends):
            raise ValueError("one end per leg is required")
        object.__setattr__(self, "dimension", Fraction(self.dimension))

    @classmethod
    def from_rule(cls, rule: VertexRule, ends: Sequence[End], *, dimension: int | Fraction = 4, order: int = 0,
                  origin: str = "") -> GraphVertex:
        """The vertex of a :class:`~hepkit_renorm.rules.VertexRule` (its legs, value and outgoing momenta)."""
        return cls(tuple(rule.legs), rule.value, tuple(ends), frozenset(rule.outgoing), Fraction(dimension), order,
                   origin)

    @property
    def fields(self) -> tuple[str, ...]:
        return tuple(leg.field for leg in self.legs)


@dataclass(frozen=True)
class GraphPropagator:
    """``D_{XY}(k) = numerator / prod (k^2 - mass2)^power`` (module docstring); ``labels`` are the labels of the two
    ends in ``numerator`` (``X`` end first) and ``k = p(PROPAGATOR_MOMENTUM, mu)``."""

    fields: tuple[str, str]
    labels: tuple[tuple[Expression, ...], tuple[Expression, ...]]
    numerator: Expression
    denominators: tuple[tuple[Expression, int], ...]

    @classmethod
    def from_rule(cls, propagator: Propagator, substitutions: Mapping[Expression, Expression] | None = None,
                  ) -> GraphPropagator:
        """Split a :class:`~hepkit_renorm.rules.Propagator` into numerator and denominators.

        ``substitutions`` (symbol -> value) are applied first, e.g. ``{model.gauge_fixing.xi: 1}`` for the Feynman
        gauge, in which the ``k^mu k^nu`` terms cancel and a common ``k^2`` cancels between numerator and denominator.
        The denominator must be a polynomial in the full ``k^2`` (``spbar(0,0) + sphat(0,0)``); its roots are the
        masses.
        """
        value = propagator.value
        for sym, val in (substitutions or {}).items():
            value = value.replace(sym, val)
        numerator, denominators = split_denominators(value)
        return cls(propagator.fields, propagator.labels, numerator, denominators)


_K2 = S("hepkit_renorm_graphs::k2")


def split_denominators(value: Expression) -> tuple[Expression, tuple[tuple[Expression, int], ...]]:
    """``value = numerator / prod (k^2 - m^2)^n`` for an expression in the propagator momentum ``p(0, mu)``."""
    k = PROPAGATOR_MOMENTUM
    v = value.replace(tn.spbar(k, k), _K2 - tn.sphat(k, k)).replace(tn.sp(k, k), _K2).together()
    num, den_factors = Expression.num(1), []
    for f in (list(v) if v.get_type() == AtomType.Mul else [v]):
        if f.get_type() == AtomType.Pow:
            base, ex = list(f)
            if ex.get_type() == AtomType.Num and ex.to_canonical_string().startswith("-"):
                den_factors.append((base, -int(ex.to_canonical_string())))
                continue
        num = num * f
    const = Expression.num(1)
    dens: list[tuple[Expression, int]] = []
    pieces: list[tuple[Expression, int]] = []
    for base, power in den_factors:
        if tn.contains_head(base, tn.sphat) or tn.contains_head(base, tn.spbar) or tn.contains_head(base, tn.p):
            raise NotImplementedError(f"denominator {base} is not a function of the full k^2")
        fac = base.factor()
        for f in (list(fac) if fac.get_type() == AtomType.Mul else [fac]):
            if f.get_type() == AtomType.Pow and list(f)[1].get_type() == AtomType.Num:
                b, ex = list(f)
                pieces.append((b, power * int(ex.to_canonical_string())))
            else:
                pieces.append((f, power))
    for base, power in pieces:
        coeffs = {str(c[0]): c[1] for c in base.expand().coefficient_list(_K2)}
        if set(coeffs) == {"1"}:
            const = const * base ** power
            continue
        if set(coeffs) - {"1", str(_K2)} or str(_K2) not in coeffs:
            raise NotImplementedError(f"denominator factor {base} is not linear in k^2")
        lead = coeffs[str(_K2)]
        mass2 = (-coeffs.get("1", Expression.num(0)) / lead).expand()
        const = const * lead ** power
        dens.append((mass2, power))
    numerator = (num / const).replace(_K2, tn.sp(k, k)).expand()
    return numerator, tuple(dens)


@dataclass(frozen=True)
class GraphRecord:
    """One graph of ``Gamma_{legs}`` (module docstring)."""

    legs: tuple[Leg, ...]
    vertices: tuple[GraphVertex, ...]
    propagators: tuple[GraphPropagator, ...]
    symmetry_factor: Fraction = Fraction(1)
    name: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        object.__setattr__(self, "symmetry_factor", Fraction(self.symmetry_factor))
        seen: dict[End, tuple[int, int]] = {}
        for v, vert in enumerate(self.vertices):
            for j, end in enumerate(vert.ends):
                if end in seen:
                    raise ValueError(f"{end} is joined twice")
                if isinstance(end, ExternalEnd):
                    if not 0 <= end.leg < len(self.legs):
                        raise ValueError(f"external leg {end.leg} does not exist")
                    if self.legs[end.leg].field != vert.legs[j].field:
                        raise ValueError(f"vertex {v} leg {j}: field {vert.legs[j].field} joined to external "
                                         f"{self.legs[end.leg].field}")
                else:
                    if not 0 <= end.propagator < len(self.propagators) or end.side not in (0, 1):
                        raise ValueError(f"{end} does not exist")
                    want = self.propagators[end.propagator].fields[end.side]
                    if want != vert.legs[j].field:
                        raise ValueError(f"vertex {v} leg {j}: field {vert.legs[j].field} joined to a {want} end")
                seen[end] = (v, j)
        missing = [ExternalEnd(j) for j in range(len(self.legs)) if ExternalEnd(j) not in seen]
        missing += [PropagatorEnd(e, s) for e in range(len(self.propagators)) for s in (0, 1)
                    if PropagatorEnd(e, s) not in seen]
        if missing:
            raise ValueError(f"unjoined ends {missing}")
        object.__setattr__(self, "_attachment", seen)

    # ------------------------------------------------------------------ topology
    def attachment(self, end: End) -> tuple[int, int]:
        """``(vertex, leg position)`` joined to ``end``."""
        return self._attachment[end]  # type: ignore[attr-defined]

    @property
    def loops(self) -> int:
        """Number of loops ``E - V + 1`` (the graph must be connected)."""
        return len(self.propagators) - len(self.vertices) + 1

    def routing(self) -> tuple[dict[int, Fraction], ...]:
        """The momentum of every propagator (module docstring), as ``{momentum number: coefficient}``."""
        nv = len(self.vertices)
        adj: dict[int, list[tuple[int, int, int]]] = {v: [] for v in range(nv)}
        for e in range(len(self.propagators)):
            a, _ = self.attachment(PropagatorEnd(e, 0))
            b, _ = self.attachment(PropagatorEnd(e, 1))
            adj[a].append((e, b, 0))
            adj[b].append((e, a, 1))
        parent: dict[int, tuple[int, int] | None] = {0: None}
        order, stack = [], [0]
        while stack:
            v = stack.pop()
            order.append(v)
            for e, w, _ in adj[v]:
                if w not in parent:
                    parent[w] = (e, v)
                    stack.append(w)
        if len(parent) != nv:
            raise ValueError("the graph is not connected")
        tree = {pe[0] for pe in parent.values() if pe is not None}
        mom: dict[int, dict[int, Fraction]] = {}
        chords = [e for e in range(len(self.propagators)) if e not in tree]
        for l_, e in enumerate(chords):
            mom[e] = {LOOP_MOMENTUM + l_: Fraction(1)}
        # incoming momentum of a leg joined to propagator end (e, s): +k_e for s = 0, -k_e for s = 1
        for v in reversed(order):
            if parent[v] is None:
                continue
            e_up, _ = parent[v]
            total: dict[int, Fraction] = {}
            up_sign = 0
            for j, end in enumerate(self.vertices[v].ends):
                if isinstance(end, ExternalEnd):
                    num = self.legs[end.leg].momentum
                    total[num] = total.get(num, Fraction(0)) + 1
                elif end.propagator == e_up:
                    up_sign += 1 if end.side == 0 else -1
                else:
                    sgn = 1 if end.side == 0 else -1
                    for n, c in mom[end.propagator].items():
                        total[n] = total.get(n, Fraction(0)) + sgn * c
            if up_sign not in (1, -1):
                raise ValueError("a tadpole edge is a tree edge")
            mom[e_up] = {n: -c / up_sign for n, c in total.items() if c}
        return tuple(dict(sorted(mom[e].items())) for e in range(len(self.propagators)))

    def leg_momentum(self, vertex: int, position: int, routing: Sequence[Momentum] | None = None) -> dict[int, Fraction]:
        """The incoming momentum of a vertex leg."""
        end = self.vertices[vertex].ends[position]
        if isinstance(end, ExternalEnd):
            return {self.legs[end.leg].momentum: Fraction(1)}
        k = (routing or self.routing())[end.propagator]
        sgn = 1 if end.side == 0 else -1
        return {n: sgn * c for n, c in k.items()}

    # ------------------------------------------------------------------ derived numbers
    def grassmann_sign(self, table: FieldTable) -> int:
        """The Wick sign of the graph (module docstring)."""
        seq: list[tuple] = []
        for vert in self.vertices:
            for leg, end in zip(vert.legs, vert.ends):
                if table.odd(leg.field):
                    seq.append(("theta", end.leg) if isinstance(end, ExternalEnd) else ("end", end.propagator,
                                                                                         end.side))
        target: list[tuple] = []
        for e, prop in enumerate(self.propagators):
            odd = [table.odd(f) for f in prop.fields]
            if odd[0] != odd[1]:
                raise ValueError(f"propagator {e} joins an even and an odd field")
            if odd[0]:
                target += [("end", e, 0), ("end", e, 1)]
        target += [("theta", j) for j, leg in enumerate(self.legs) if table.odd(leg.field)]
        if sorted(seq) != sorted(target):
            raise ValueError("Grassmann-odd legs and odd ends do not match")
        pos = [seq.index(x) for x in target]
        inv = sum(1 for i in range(len(pos)) for j in range(i + 1, len(pos)) if pos[i] > pos[j])
        return -1 if inv % 2 else 1

    def superficial_degree(self, table: FieldTable) -> Fraction:
        """``omega(G)`` (module docstring)."""
        omega = Fraction(4 * self.loops)
        for vert in self.vertices:
            omega += vert.dimension - sum(table[f].dimension for f in vert.fields)
        for prop in self.propagators:
            omega += table[prop.fields[0]].dimension + table[prop.fields[1]].dimension - 4
        return omega

    def is_one_particle_irreducible(self) -> bool:
        """True if removing any single propagator leaves the graph connected."""
        for cut in range(len(self.propagators)):
            seen, stack = {0}, [0]
            while stack:
                v = stack.pop()
                for e in range(len(self.propagators)):
                    if e == cut:
                        continue
                    a, _ = self.attachment(PropagatorEnd(e, 0))
                    b, _ = self.attachment(PropagatorEnd(e, 1))
                    for x, y in ((a, b), (b, a)):
                        if x == v and y not in seen:
                            seen.add(y)
                            stack.append(y)
            if len(seen) != len(self.vertices):
                return False
        return True
