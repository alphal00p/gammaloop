"""Feynman graphs of a class of external legs, enumerated from the vertex rules of the order store.

The graphs of report 02_method (sec:m-algorithm steps 2 and 7, sec:m-beyond) are enumerated directly from the
package's rule sets: every vertex is a :class:`~hepkit_renorm.rules.VertexRule` of a store entry (``S_0``,
``Delta_1``, later the counterterms and breaking insertions of lower orders), every line a propagator of the tree
action. Nothing is typed in: the vertices, their Grassmann order, the propagators, the symmetry factors and the
Wick signs follow from the rules (:mod:`hepkit_renorm.graphs` holds the sign rule and the power counting). The
graphs are multiplet-level: a line or a vertex carries the full group and flavour structure of its fields, which
the loop back-ends contract as tensors (:mod:`hepkit_renorm.loops`).

Families of vertices (:class:`VertexFamily`, :func:`store_families`)
--------------------------------------------------------------------
A family is the rule set of one Laurent coefficient ``C_k`` of one store entry, with its loop order ``L_v`` (the
power of ``hbar`` it carries), its ``epsilon`` power ``k``, its mass dimension (``4 + g dim s`` for ghost number
``g``, :func:`~hepkit_renorm.classes.dimension_bound`) and whether it is a Slavnov-Taylor insertion (role
``breaking``: the graphs of ``Delta . Gamma`` are linear in the insertion). The tree entry gives the vertex classes
of :func:`~hepkit_renorm.rules.tree_rules` (at least three legs, or a BRST source); every other entry gives all
its classes, two-point ones included (counterterm and breaking insertions on lines).

Which graphs (:func:`enumerate_graphs`)
---------------------------------------
For external legs ``Phi_1 .. Phi_n`` (in the Grassmann order of the vertex function) and a total order ``l``:
connected graphs with ``L >= min_loops`` loops (or exactly ``loops``) whose vertices satisfy ``L + sum_v L_v = l``, with exactly
``insertions`` vertices of insertion families, one-particle irreducible unless ``one_particle_irreducible=False``.
Lines are propagators of the tree action, so BRST sources (no propagator) are external only; every external leg is
attached to a vertex leg of the same field (amputated graphs). At ``l = 1`` this gives the one-loop graphs of
``S_0`` (ghost number 0, ``insertions = 0``) and the graphs with exactly one ``Delta_1`` vertex (``insertions =
1``); the forest of their divergent subgraphs is empty.

Algorithm: half-edges are joined one at a time (first the external legs in order, then the lowest open vertex leg):
to an open leg of an existing vertex through a propagator, or to a new vertex of a family. Legs of one vertex with
the same field are interchangeable (the vertex functions are graded symmetric in them), so only the first open leg
of a field is tried. The search is finite: ``sum_v (d_v - 2) = 2 L - 2 + n`` (Euler) bounds the vertices of degree
``>= 3``; two-point and one-point vertices exist only in families of loop order ``>= 1`` or insertion families
(bounded by ``l`` and ``insertions``) or carry a source (bounded by the external sources). Each graph is generated
several times and kept once, by its canonical key.

Canonical key and symmetry factor
---------------------------------
A graph is determined by its vertices (family, rule) with the external legs they carry and the multiset of lines
``((vertex, field), (vertex, field))`` (the ``X`` end first; both ends sorted for a self-conjugate line). The
canonical key (:attr:`Topology.key`) is the smallest such description over the relabellings of equivalent vertices.
The symmetry factor is ``1/|Aut|``, ``Aut`` the group of permutations of vertices and of same-field legs of each
vertex that fix the external legs and map lines to lines (:func:`automorphisms`); with the fully symmetrised vertex
functions ``V`` this is the factor of Wick's theorem. An automorphism that permutes the Grassmann-odd legs by an
odd permutation would make the graph vanish; such graphs are reported (``vanishing``) and dropped.
"""

from __future__ import annotations

import itertools
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction

from symbolica import Expression

from . import rules as rl
from .classes import dimension_bound
from .fields import FieldTable
from .functional import LocalFunctional, standard_legs
from .graphs import ExternalEnd, GraphPropagator, GraphRecord, GraphVertex, PropagatorEnd
from .spec import Model
from .store import OrderStore, Role, StoreEntry

__all__ = [
    "Enumeration",
    "Topology",
    "VertexFamily",
    "automorphisms",
    "enumerate_graphs",
    "graph_propagators",
    "line_fields",
    "line_keys",
    "record",
    "store_families",
    "topology_key",
]

Slot = tuple[int, int]
"""A vertex leg ``(vertex, leg position)``."""


# ---------------------------------------------------------------------------------------------- families
@dataclass(frozen=True)
class VertexFamily:
    """The vertex rules of one Laurent coefficient of one store entry (module docstring).

    ``rules`` maps the leg class (canonical field order) to its rule; ``loop_order`` is ``L_v``, ``power`` the
    ``epsilon`` power of the coefficient, ``dimension`` the mass dimension of the functional, ``insertion`` marks a
    Slavnov-Taylor insertion (role ``breaking``).
    """

    name: str
    rules: Mapping[tuple[str, ...], rl.VertexRule]
    dimension: Fraction = Fraction(4)
    loop_order: int = 0
    power: int = 0
    insertion: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "dimension", Fraction(self.dimension))


def _ghost_number(table: FieldTable, functional: LocalFunctional) -> int:
    values = {sum(table[f].ghost_number for f in cls) for cls in functional.leg_classes()}
    if len(values) > 1:
        raise ValueError(f"functional of mixed ghost number {sorted(values)}")
    return values.pop() if values else 0


def _entry_rules(entry: StoreEntry, functional: LocalFunctional, model: Model) -> dict[tuple[str, ...], rl.VertexRule]:
    out = rl.tree_rules(functional, model)
    if entry.role != Role.TREE:
        simp = model.simplifier()
        for cls in functional.leg_classes():
            if cls not in out and not rl.is_vertex_class(model.fields, cls):
                rule = rl.vertex_rule(functional, cls, simp)
                if rule.value.expand() != Expression.num(0):
                    out[cls] = rule
    return out


def store_families(model: Model, store: OrderStore, order: int, bookkeeping: str = "A") -> list[VertexFamily]:
    """The vertex families available to the graphs of total order ``order``: every Laurent coefficient of every
    entry of loop order ``< order`` (module docstring) that belongs to both bookkeepings or to ``bookkeeping``
    (``"A"``: the insertions ``b_d deltaS``, ``b_d DeltaS``; ``"B"``: the remainder ``hatDelta``; report
    sec:m-beyond tab:m-bookkeeping; the two must never be inserted together). At ``order = 1`` these are ``S_0``
    and ``Delta_1``."""
    if bookkeeping not in ("A", "B"):
        raise ValueError("bookkeeping is 'A' or 'B'")
    out = []
    for entry in store.entries(max_loop_order=order - 1, bookkeeping=bookkeeping):
        for k in entry.powers:
            functional = entry.coefficient(k)
            g = _ghost_number(model.fields, functional)
            name = entry.name if entry.powers == (0,) else f"{entry.name}[eps^{k}]"
            out.append(VertexFamily(name, _entry_rules(entry, functional, model), dimension_bound(model, g),
                                    entry.loop_order, k, entry.role == Role.BREAKING))
    return out


def graph_propagators(model: Model, store: OrderStore, substitutions: Mapping[Expression, Expression] | None = None,
                      ) -> dict[tuple[str, str], GraphPropagator]:
    """The propagators of the tree action of the store, split for the back-ends; ``substitutions`` default to the
    Feynman gauge ``{xi: 1}``."""
    subs = {model.gauge_fixing.xi: Expression.num(1)} if substitutions is None else substitutions
    return {p.fields: GraphPropagator.from_rule(p, subs) for p in rl.propagators(store.tree_action(), model)}


# ---------------------------------------------------------------------------------------------- topologies
@dataclass(frozen=True)
class Topology:
    """A graph before evaluation: vertices ``(family index, leg class)``, ``externals[k]`` = the slot of external
    leg ``k``, ``lines`` = ``(X slot, Y slot, (X, Y))`` per propagator, and the derived numbers."""

    legs: tuple[str, ...]
    vertices: tuple[tuple[int, tuple[str, ...]], ...]
    externals: tuple[Slot, ...]
    lines: tuple[tuple[Slot, Slot, tuple[str, str]], ...]
    key: tuple
    symmetry_factor: Fraction
    loops: int
    order: int
    power: int

    def describe(self, families: Sequence[VertexFamily]) -> str:
        """A short readable form: ``[fam:fields @ext] ... | lines``."""
        vs = []
        for v, (f, cls) in enumerate(self.vertices):
            ext = [str(k) for k, s in enumerate(self.externals) if s[0] == v]
            vs.append(f"[{families[f].name}:{' '.join(cls)}{' @' + ','.join(ext) if ext else ''}]")
        ls = [f"{a[0]}-{b[0]}:{x}" + ("" if x == y else f"/{y}") for a, b, (x, y) in self.lines]
        return " ".join(vs) + " | " + " ".join(ls)


@dataclass
class Enumeration:
    """Result of :func:`enumerate_graphs`: the kept topologies, the vanishing ones (odd automorphism) and the
    families they refer to."""

    legs: tuple[str, ...]
    families: tuple[VertexFamily, ...]
    topologies: list[Topology] = field(default_factory=list)
    vanishing: list[Topology] = field(default_factory=list)

    def records(self, model: Model, propagators: Mapping[tuple[str, str], GraphPropagator]) -> list[GraphRecord]:
        """The graph records of the kept topologies (:func:`record`)."""
        return [record(t, self.families, model, propagators) for t in self.topologies]


def _partners(propagators: Iterable[tuple[str, str]]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for x, y in propagators:
        out.setdefault(x, set()).add(y)
        out.setdefault(y, set()).add(x)
    return out


def line_keys(propagators: Iterable[tuple[str, str]]) -> list[tuple[str, str]]:
    """The field pairs ``(X, Y)`` that name lines: one orientation per pair. A bosonic pair of different fields
    (``H``, ``H^dagger``) has both ``D_{XY}`` and ``D_{YX}`` (the same line read from either end); the first in the
    given order is kept, so that every line has one description and the canonical key is unique."""
    out: list[tuple[str, str]] = []
    for x, y in propagators:
        if (x, y) not in out and (y, x) not in out:
            out.append((x, y))
    return out


def line_fields(a: str, b: str, keys: Iterable[tuple[str, str]]) -> tuple[bool, tuple[str, str]]:
    """``(a is the X end, (X, Y))`` for a line between legs of fields ``a`` and ``b`` (``keys``: from
    :func:`line_keys`)."""
    if (a, b) in keys:
        return True, (a, b)
    if (b, a) in keys:
        return False, (b, a)
    raise KeyError(f"no propagator between {a} and {b}")


class _Search:
    def __init__(self, table: FieldTable, legs: Sequence[str], families: Sequence[VertexFamily],
                 propagators: Iterable[tuple[str, str]], order: int, insertions: int | None, min_loops: int,
                 loops: int | None, one_particle_irreducible: bool) -> None:
        self.table = table
        self.legs = tuple(legs)
        self.families = tuple(families)
        self.keys = set(line_keys(propagators))
        self.partners = _partners(self.keys)
        self.order, self.insertions = order, insertions
        self.min_loops = loops if loops is not None else min_loops
        self.exact_loops = loops
        self.opi = one_particle_irreducible
        self.n_sources = sum(1 for f in legs if table[f].is_source)
        by_field: dict[str, list[tuple[int, tuple[str, ...], int]]] = {}
        for fi, fam in enumerate(self.families):
            for cls in sorted(fam.rules, key=lambda c: [table.order(x) for x in c]):
                for j, f in enumerate(cls):
                    if j == cls.index(f):
                        by_field.setdefault(f, []).append((fi, cls, j))
        self.by_field = by_field
        # Euler: sum_v (d_v - 2) = 2 L - 2 + n with L <= order; one-point vertices lower
        # the sum by one each, and exist only in families of loop order >= 1 or in insertion families
        tadpoles = [fam for fam in self.families if any(len(c) == 1 for c in fam.rules)]
        slack = 0
        if any(fam.insertion for fam in tadpoles):
            slack += insertions if insertions is not None else order
        if any(not fam.insertion for fam in tadpoles):
            slack += order
        self.max_degree_sum = 2 * (loops if loops is not None else order) - 2 + len(self.legs) + slack
        self.degree_sum = self.orders = self.n_insertions = self.sources = 0
        self.vertices: list[tuple[int, tuple[str, ...]]] = []
        self.open: list[list[bool]] = []
        self.externals: list[Slot | None] = [None] * len(self.legs)
        self.lines: list[tuple[Slot, Slot]] = []
        self.found: dict[tuple, tuple] = {}

    # ------------------------------------------------------------------ bookkeeping
    def _budget_ok(self, fi: int, cls: tuple[str, ...]) -> bool:
        fam = self.families[fi]
        if self.orders + fam.loop_order + self.min_loops > self.order:
            return False
        if fam.insertion and self.insertions is not None and self.n_insertions + 1 > self.insertions:
            return False
        if self.degree_sum + len(cls) - 2 > self.max_degree_sum:
            return False
        return self.sources + self._n_sources(cls) <= self.n_sources

    def _n_sources(self, cls: tuple[str, ...]) -> int:
        return sum(1 for f in cls if self.table[f].is_source)

    def _add(self, fi: int, cls: tuple[str, ...]) -> int:
        fam = self.families[fi]
        self.vertices.append((fi, cls))
        self.open.append([True] * len(cls))
        self.degree_sum += len(cls) - 2
        self.orders += fam.loop_order
        self.n_insertions += fam.insertion
        self.sources += self._n_sources(cls)
        return len(self.vertices) - 1

    def _pop(self) -> None:
        fi, cls = self.vertices.pop()
        fam = self.families[fi]
        self.open.pop()
        self.degree_sum -= len(cls) - 2
        self.orders -= fam.loop_order
        self.n_insertions -= fam.insertion
        self.sources -= self._n_sources(cls)

    def _first_open(self, v: int, f: str, exclude: Slot | None = None) -> int | None:
        for j, (g, free) in enumerate(zip(self.vertices[v][1], self.open[v])):
            if free and g == f and (v, j) != exclude:
                return j
        return None

    def _dead(self, v: int, skip: int) -> bool:
        """True if a new vertex has a leg that can never be joined (no propagator, externals all placed)."""
        cls = self.vertices[v][1]
        return any(j != skip and f not in self.partners for j, f in enumerate(cls))

    # ------------------------------------------------------------------ search
    def run(self) -> None:
        self._step()

    def _step(self) -> None:
        if None in self.externals:
            self._place_external(self.externals.index(None))
            return
        slot = next(((v, j) for v, fr in enumerate(self.open) for j, free in enumerate(fr) if free), None)
        if slot is None:
            self._complete()
            return
        v, j = slot
        f = self.vertices[v][1][j]
        if f not in self.partners:
            return
        self.open[v][j] = False
        for g in sorted(self.partners[f], key=self.table.order):
            for w in range(len(self.vertices)):
                i = self._first_open(w, g, exclude=slot)
                if i is None:
                    continue
                self.open[w][i] = False
                self.lines.append(((v, j), (w, i)))
                self._step()
                self.lines.pop()
                self.open[w][i] = True
            for fi, cls, i in self.by_field.get(g, ()):
                if not self._budget_ok(fi, cls):
                    continue
                w = self._add(fi, cls)
                if not self._dead(w, i):
                    self.open[w][i] = False
                    self.lines.append(((v, j), (w, i)))
                    self._step()
                    self.lines.pop()
                self._pop()
        self.open[v][j] = True

    def _place_external(self, k: int) -> None:
        f = self.legs[k]
        for v in range(len(self.vertices)):
            j = self._first_open(v, f)
            if j is None:
                continue
            self.open[v][j] = False
            self.externals[k] = (v, j)
            self._step()
            self.externals[k] = None
            self.open[v][j] = True
        for fi, cls, j in self.by_field.get(f, ()):
            if not self._budget_ok(fi, cls):
                continue
            v = self._add(fi, cls)
            self.open[v][j] = False
            self.externals[k] = (v, j)
            self._step()
            self.externals[k] = None
            self._pop()

    def _complete(self) -> None:
        nv, ne = len(self.vertices), len(self.lines)
        loops = ne - nv + 1
        if loops < self.min_loops or (self.exact_loops is not None and loops != self.exact_loops):
            return
        orders = sum(self.families[f].loop_order for f, _ in self.vertices)
        if loops + orders != self.order:
            return
        n_ins = sum(1 for f, _ in self.vertices if self.families[f].insertion)
        if self.insertions is not None and n_ins != self.insertions:
            return
        ends = [(a[0], b[0]) for a, b in self.lines]
        if not _connected(nv, ends) or (self.opi and not all(_connected(nv, ends[:e] + ends[e + 1:])
                                                             for e in range(ne))):
            return
        lines = []
        for a, b in self.lines:
            fa, fb = self.vertices[a[0]][1][a[1]], self.vertices[b[0]][1][b[1]]
            a_is_x, xy = line_fields(fa, fb, self.keys)
            lines.append((a, b, xy) if a_is_x else (b, a, xy))
        key, order = _canonical(self.vertices, self.externals, lines)
        if key in self.found:
            return
        # relabel the vertices into the canonical order
        new = {old: n for n, old in enumerate(order)}
        verts = tuple(self.vertices[old] for old in order)
        exts = tuple((new[v], j) for v, j in self.externals)
        lns = tuple(sorted((((new[a[0]], a[1]), (new[b[0]], b[1]), xy) for a, b, xy in lines)))
        self.found[key] = (verts, exts, lns, loops, orders)


def _connected(nv: int, ends: Sequence[tuple[int, int]]) -> bool:
    seen, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for a, b in ends:
            for x, y in ((a, b), (b, a)):
                if x == v and y not in seen:
                    seen.add(y)
                    stack.append(y)
    return len(seen) == nv


def _vertex_labels(vertices, externals) -> list[tuple]:
    ext: dict[int, list[int]] = {}
    for k, (v, _) in enumerate(externals):
        ext.setdefault(v, []).append(k)
    return [(f, cls, tuple(sorted(ext.get(v, ())))) for v, (f, cls) in enumerate(vertices)]


def _line_key(lines, perm: Mapping[int, int], vertices) -> tuple:
    out = []
    for a, b, (x, y) in lines:
        ea = (perm[a[0]], vertices[a[0]][1][a[1]])
        eb = (perm[b[0]], vertices[b[0]][1][b[1]])
        out.append((ea, eb) if x != y else tuple(sorted((ea, eb))))
    return tuple(sorted(out))


def _relabellings(labels: Sequence[tuple]) -> tuple[list[int], Iterable[dict[int, int]]]:
    """Vertices sorted by label, and every bijection old -> new position that keeps the sorted labels."""
    base = sorted(range(len(labels)), key=lambda v: (repr(labels[v]), v))
    blocks: list[list[int]] = []
    for v in base:
        if blocks and labels[blocks[-1][0]] == labels[v]:
            blocks[-1].append(v)
        else:
            blocks.append([v])

    def gen():
        for choice in itertools.product(*(itertools.permutations(b) for b in blocks)):
            order = [v for block in choice for v in block]
            yield {old: n for n, old in enumerate(order)}
    return base, gen()


def _canonical(vertices, externals, lines) -> tuple[tuple, list[int]]:
    labels = _vertex_labels(vertices, externals)
    _, perms = _relabellings(labels)
    best, best_perm = None, None
    sorted_labels = tuple(repr(lab) for lab in sorted(labels, key=repr))
    for perm in perms:
        key = (sorted_labels, _line_key(lines, perm, vertices))
        if best is None or key < best:
            best, best_perm = key, perm
    order = sorted(best_perm, key=best_perm.get)
    return best, order


def topology_key(vertices: Sequence[tuple[int, tuple[str, ...]]], externals: Sequence[Slot],
                 lines: Sequence[tuple[Slot, Slot, tuple[str, str]]]) -> tuple:
    """The canonical key (module docstring) of a graph given by its vertices ``(family index, leg class)``, the
    slots of the external legs and its lines ``(X slot, Y slot, (X, Y))``; equal keys iff isomorphic graphs."""
    return _canonical(list(vertices), list(externals), list(lines))[0]


def automorphisms(top: Topology) -> list[tuple[dict[int, int], dict[Slot, Slot]]]:
    """The automorphisms of a topology as ``(vertex map, slot map)``: permutations of equivalent vertices and of
    same-field legs that fix the external legs and map the lines onto themselves (module docstring)."""
    labels = _vertex_labels(top.vertices, top.externals)
    _, perms = _relabellings(labels)
    target = _line_key(top.lines, {v: v for v in range(len(top.vertices))}, top.vertices)
    out = []
    for perm in perms:
        if _line_key(top.lines, perm, top.vertices) == target:
            out.extend((dict(perm), smap) for smap in _slot_maps(top, perm))
    return out


def _slot_maps(top: Topology, perm: Mapping[int, int]) -> list[dict[Slot, Slot]]:
    """The leg bijections over a vertex permutation that fix the external legs and preserve the lines."""
    line_set = Counter(_slot_pair(a, b, xy) for a, b, xy in top.lines)
    ext_slot = {s: k for k, s in enumerate(top.externals)}
    slots = [(v, j) for v in range(len(top.vertices)) for j in range(len(top.vertices[v][1]))]
    out: list[dict[Slot, Slot]] = []
    smap: dict[Slot, Slot] = {}
    used: set[Slot] = set()

    def extend(i: int) -> None:
        if i == len(slots):
            if Counter(_slot_pair(smap[a], smap[b], xy) for a, b, xy in top.lines) == line_set:
                out.append(dict(smap))
            return
        v, j = slots[i]
        w, f = perm[v], top.vertices[v][1][j]
        if (v, j) in ext_slot:
            cands = [top.externals[ext_slot[(v, j)]]]
        else:
            cands = [(w, k) for k, g in enumerate(top.vertices[w][1]) if g == f and (w, k) not in ext_slot]
        for c in cands:
            if c not in used:
                smap[(v, j)] = c
                used.add(c)
                extend(i + 1)
                used.discard(c)
                del smap[(v, j)]
    extend(0)
    return out


def _slot_pair(a: Slot, b: Slot, xy: tuple[str, str]) -> tuple:
    return (a, b) if xy[0] != xy[1] else tuple(sorted((a, b)))


def _odd_parity(top: Topology, table: FieldTable, smap: Mapping[Slot, Slot]) -> int:
    seq = [(v, j) for v in range(len(top.vertices)) for j, f in enumerate(top.vertices[v][1]) if table.odd(f)]
    pos = {s: n for n, s in enumerate(seq)}
    perm = [pos[smap[s]] for s in seq]
    inv = sum(1 for i in range(len(perm)) for k in range(i + 1, len(perm)) if perm[i] > perm[k])
    return -1 if inv % 2 else 1


def enumerate_graphs(model: Model, legs: Sequence[str], families: Sequence[VertexFamily],
                     propagators: Iterable[tuple[str, str]], *, order: int = 1, insertions: int | None = 0,
                     min_loops: int = 1, loops: int | None = None, one_particle_irreducible: bool = True,
                     ) -> Enumeration:
    """All graphs of the leg class ``legs`` (module docstring), deterministic order (by canonical key);
    ``loops`` restricts the loop number to one value (e.g. the one-loop graphs with counterterm insertions at
    ``order = 2``)."""
    table = model.fields
    search = _Search(table, legs, families, propagators, order, insertions, min_loops, loops,
                     one_particle_irreducible)
    search.run()
    out = Enumeration(tuple(legs), tuple(families))
    for key in sorted(search.found, key=repr):
        verts, exts, lns, loops, orders = search.found[key]
        power = sum(families[f].power for f, _ in verts)
        top = Topology(tuple(legs), verts, exts, lns, key, Fraction(1), loops, loops + orders, power)
        auts = automorphisms(top)
        top = Topology(top.legs, verts, exts, lns, key, Fraction(1, len(auts)), loops, loops + orders, power)
        if any(_odd_parity(top, table, smap) < 0 for _, smap in auts):
            out.vanishing.append(top)
        else:
            out.topologies.append(top)
    return out


def record(top: Topology, families: Sequence[VertexFamily], model: Model,
           propagators: Mapping[tuple[str, str], GraphPropagator], name: str = "") -> GraphRecord:
    """The :class:`~hepkit_renorm.graphs.GraphRecord` of a topology: vertices in topology order with their rules,
    propagators in line order, the symmetry factor; the Wick sign is derived by the record."""
    ends: dict[Slot, object] = {}
    for k, s in enumerate(top.externals):
        ends[s] = ExternalEnd(k)
    props = []
    for e, (x_slot, y_slot, xy) in enumerate(top.lines):
        ends[x_slot] = PropagatorEnd(e, 0)
        ends[y_slot] = PropagatorEnd(e, 1)
        props.append(propagators[xy])
    verts = []
    for v, (fi, cls) in enumerate(top.vertices):
        fam = families[fi]
        verts.append(GraphVertex.from_rule(fam.rules[cls], [ends[(v, j)] for j in range(len(cls))],
                                           dimension=fam.dimension, order=fam.loop_order, origin=fam.name))
    return GraphRecord(tuple(standard_legs(model.fields, top.legs)), tuple(verts), tuple(props),
                       symmetry_factor=top.symmetry_factor, name=name or top.describe(families))
