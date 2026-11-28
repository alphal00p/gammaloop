"""Graphs from the HepKIT generator on the written model, mapped back to the package's multiplet-level graphs.

This is the generator route of report 02_method step 2 (the second, independent enumeration next to
:func:`hepkit_renorm.diagrams.enumerate_graphs`). The written model (:mod:`hepkit_renorm.output`) encodes the
SU(2) and flavour slots in components, so a leg class of package fields becomes a set of *component tuples*
(:func:`component_tuples`: every choice of written particles for the legs). For each tuple HepKIT generates the
one-loop 1PI diagrams; every diagram is mapped back through the writer's :class:`~hepkit_renorm.output.ReverseMap`
(vertex -> store entry and package leg class, particle -> package field) to the canonical key of a multiplet-level
graph (:func:`~hepkit_renorm.diagrams.topology_key`). The union over the tuples is the set of multiplet-level
graphs that HepKIT sees.

Generation settings: today's workarounds (scout_generation notes, sec. 5.2, and the written model's coupling
orders):

* ``loops = L``, 1PI (default bridges), ``self_energy = tadpoles = zero_snails = None``;
* coupling orders: the breaking order (``DELTA``) equal to the number of insertions, the source order (``SRC``)
  equal to the number of external sources (every source vertex is linear in the sources, so no source can be
  internal);
* ``max_vertices = 2 L - 2 + n + (l - L)`` (Euler's formula for vertices of degree ``>= 3``; the two-point vertices
  of the tree model carry a source, which is external, and each counterterm vertex carries at least one unit of the
  order ``CT``), and a filter that rejects internal source lines;
* at order ``l > L`` the counterterm order ``CT = l - L``.

A multiplet-level graph is missing from HepKIT's set if no component assignment of its lines has non-zero
components at all its vertices (a group-theoretical zero such as ``Tr tau^a = 0``); such graphs must have vanishing
values. HepKIT's symmetry factors and fermion signs refer to the component encoding (commuting ghosts) and are not
used: the package derives them (:mod:`hepkit_renorm.diagrams`, :mod:`hepkit_renorm.graphs`). A written vertex is
mapped to the family of its store entry with ``eps`` power 0; an entry with several powers is written as one vertex
with a Laurent coupling, which a multiplet-level key cannot split (such entries raise ``KeyError``).
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction

from .diagrams import VertexFamily, line_fields, line_keys, topology_key
from .output import WrittenModel
from .spec import Model

__all__ = ["Generated", "component_tuples", "generated_keys"]


@dataclass(frozen=True)
class Generated:
    """One component diagram: the external particles, the particles of its internal lines and the symmetry factor
    ``1/|Aut|`` of the component graph (HepKIT's ``symmetry_factor`` is ``|Aut|``); it equals the multiplet-level
    factor when no line carries components."""

    particles: tuple[str, ...]
    lines: tuple[str, ...]
    symmetry_factor: Fraction


def component_tuples(written: WrittenModel, fields: Sequence[str]) -> list[tuple[str, ...]]:
    """Every choice of written particles (one per leg) for a class of package fields, in leg order; none if a
    field is not written (a field without vertices, such as the source of an abelian ghost)."""
    rmap = written.reverse_map
    options = []
    for f in fields:
        names = sorted(n for n, rec in rmap.particles.items() if rec.field == f)
        if not names:
            return []
        options.append(names)
    return list(itertools.product(*options))


def _process(written: WrittenModel, particles: Sequence[str]) -> tuple[list[str], list[str], dict[int, int]]:
    """Incoming and outgoing particle lists and the map HepKIT external index -> leg position. A leg of an
    antiparticle (negative PDG code) is the outgoing particle."""
    rmap = written.reverse_map
    inc, out = [], []
    for k, name in enumerate(particles):
        rec = rmap.particles[name]
        if rec.pdg_code < 0:
            out.append((k, rec.antiname))
        else:
            inc.append((k, name))
    order = inc + out
    return [n for _, n in inc], [n for _, n in out], {i: k for i, (k, _) in enumerate(order)}


def _diagram_key(diagram, written: WrittenModel, families: Sequence[VertexFamily], keys: set[tuple[str, str]],
                 index: Mapping[int, int]) -> tuple:
    rmap = written.reverse_map
    fam_of = {fam.name: i for i, fam in enumerate(families) if fam.power == 0}
    vid = {v.id: n for n, v in enumerate(diagram.vertices)}
    vertices, free = [], []
    for v in diagram.vertices:
        rec = rmap.vertex(v.interaction)
        if rec.entry not in fam_of:
            raise KeyError(f"vertex {v.interaction} of entry {rec.entry} is not in the families")
        vertices.append((fam_of[rec.entry], tuple(rec.fields)))
        free.append(list(range(len(rec.fields))))

    def take(v: int, field: str) -> tuple[int, int]:
        cls = vertices[v][1]
        for j in free[v]:
            if cls[j] == field:
                free[v].remove(j)
                return v, j
        raise ValueError(f"vertex {diagram.vertices[v].interaction} has no free {field} leg")

    def end_fields(e) -> tuple[str, str]:
        rec = rmap.particles[e.particle_name]
        return rmap.particles[rec.antiname].field, rec.field          # (field at the source, at the target)

    externals: dict[int, tuple[int, int]] = {}
    lines = []
    for e in diagram.edges:
        at_source, at_target = end_fields(e)
        if e.is_external:
            if e.source is None:
                externals[index[e.external_index]] = take(vid[e.target], at_target)
            else:
                externals[index[e.external_index]] = take(vid[e.source], at_source)
            continue
        a, b = take(vid[e.source], at_source), take(vid[e.target], at_target)
        a_is_x, xy = line_fields(at_source, at_target, keys)
        lines.append((a, b, xy) if a_is_x else (b, a, xy))
    return topology_key(vertices, [externals[k] for k in range(len(externals))], lines)


def generated_keys(written: WrittenModel, model: Model, fields: Sequence[str], families: Sequence[VertexFamily],
                   propagators: Iterable[tuple[str, str]], *, loops: int = 1, order: int | None = None,
                   insertions: int = 0, tuples: Sequence[Sequence[str]] | None = None) -> dict[tuple, list[Generated]]:
    """The canonical keys of the multiplet-level graphs generated by HepKIT for a class (module docstring), each
    with its component diagrams. ``order`` (default ``loops``) is the total order ``l``: the counterterm vertices
    carry ``l - loops`` units of the written loop order (``CT``). ``tuples`` defaults to :func:`component_tuples`."""
    hk_model = written.load()
    rmap = written.reverse_map
    orders = written.options.orders
    keys = set(line_keys(propagators))
    source_pdgs = {abs(rec.pdg_code) for rec in rmap.particles.values() if model.fields[rec.field].is_source}
    n_sources = sum(1 for f in fields if model.fields[f].is_source)

    def no_internal_sources(graph, _n) -> bool:
        nodes = graph.nodes()
        for u, v, _d, pdg in graph.edges():
            if abs(int(str(pdg))) in source_pdgs and int(str(nodes[u][1])) == 0 and int(str(nodes[v][1])) == 0:
                return False
        return True

    ct = (loops if order is None else order) - loops
    coupling_orders = {orders.breaking: insertions, orders.source: n_sources}
    if orders.loop in rmap.orders:
        coupling_orders[orders.loop] = ct
    options = {"loops": loops, "coupling_orders": coupling_orders,
               "max_vertices": 2 * loops - 2 + len(fields) + ct, "self_energy": None, "tadpoles": None,
               "zero_snails": None, "progress": None}
    if n_sources:
        options["filter"] = no_internal_sources
    out: dict[tuple, list[Generated]] = {}
    for particles in (tuples if tuples is not None else component_tuples(written, fields)):
        inc, outg, index = _process(written, particles)
        for d in hk_model.process(inc, outg).generate_diagrams(**options):
            lines = tuple(sorted(e.particle_name for e in d.internal_edges))
            out.setdefault(_diagram_key(d, written, families, keys, index), []).append(
                Generated(tuple(particles), lines, 1 / Fraction(str(d.symmetry_factor))))
    return out
