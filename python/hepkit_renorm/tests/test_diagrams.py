"""Graph enumeration from the order store (hepkit_renorm.diagrams): counts, symmetry factors, canonical keys, Wick
signs, families and the order-l interface.

References:

* ``tests/data/enumeration_counts.json``: labelled multiplet-level one-loop 1PI graph counts per assignment of the
  hand enumeration (``SCRATCHPAD/model/enumeration.json``, HepKIT on the hand-written model). The hand enumeration
  used component tuples up to relabelling of SU(2) adjoint values, which misses graphs that vanish by group theory
  in the chosen representative tuple (``Tr tau^a = 0``, ``f^{abc} {T^b, T^c} = 0``): ``EXTRA`` lists them; their
  values vanish (``test_class_results.py``), and HepKIT on all component tuples finds them
  (``test_hepkit_graphs.py``).
* Hand-derived Wick signs of the class results (``SCRATCHPAD/classes/<key>/final_result.json``, diagram tables): in
  the report's conventions (odd fields of a vertex in the order of its Lagrangian term, external legs in the order of
  the class header) every graph with an open odd line has the sign ``+1`` and every closed fermion or ghost loop
  ``-1``; with two external ghosts the graph with the first ghost on the ``Delta_1`` vertex has ``+1`` (g1_ccRF).

Enumerating every class of ghost number 0 and 1 takes about 15 s (module fixture, shared by the tests).
"""

from __future__ import annotations

import itertools
import json
from collections import Counter
from fractions import Fraction

import pytest
from conftest import DATA

from hepkit_renorm import LocalFunctional
from hepkit_renorm.action import gauge_kinetic
from hepkit_renorm.brst import BRST
from hepkit_renorm.classes import leg_classes
from hepkit_renorm.diagrams import (
    automorphisms,
    enumerate_graphs,
    graph_propagators,
    line_keys,
    record,
    store_families,
    topology_key,
)
from hepkit_renorm.functional import Density
from hepkit_renorm.store import OrderStore, Role, StoreEntry, tree_store

# graphs that the hand enumeration misses (group-theoretical zeros in its representative component tuples)
EXTRA = {("W", "W", "W"): 3, ("W", "W", "W", "B"): 9, ("Hdag", "H", "W"): 1}


@pytest.fixture(scope="module")
def setup(sm):
    brst = BRST(sm)
    store = tree_store(sm, brst)
    return brst, store, store_families(sm, store, 1), graph_propagators(sm, store)


@pytest.fixture(scope="module")
def graphs(sm, setup):
    """Enumeration of every derived class of ghost number 0 and 1 (canonical leg order)."""
    _, _, fams, props = setup
    out = {}
    for g in (0, 1):
        for cls in leg_classes(sm, g):
            out[cls.fields] = enumerate_graphs(sm, cls.fields, fams, props, insertions=g)
    return out


def _canonical(sm, fields):
    return tuple(sorted(fields, key=sm.fields.order))


def test_families(sm, setup):
    _, _, fams, props = setup
    assert [(f.name, f.dimension, f.loop_order, f.power, f.insertion) for f in fams] == [
        ("S0", 4, 0, 0, False), ("Delta1", 5, 0, 0, True)]
    s0 = fams[0].rules
    two_point = [c for c in s0 if len(c) == 2]
    assert two_point and all(any(sm.fields[f].is_source for f in c) for c in two_point)     # rho omega, ...
    assert ("G", "G") in props and ("Q", "Qbar") in props and ("gh_G", "agh_G") in props
    assert not any(sm.fields[f].is_source for pair in props for f in pair)                   # sources never propagate
    keys = line_keys(props)
    assert ("H", "Hdag") in keys and ("Hdag", "H") not in keys and ("Hdag", "H") in props      # one name per line


def test_counts_against_hand_enumeration(sm, setup, graphs):
    """V3/T3: every assignment of the eighteen classes has the hand enumeration's number of graphs, plus EXTRA."""
    _, _, fams, props = setup
    with open(DATA / "enumeration_counts.json") as fh:
        ref = json.load(fh)["classes"]
    total = 0
    for cls, rows in ref.items():
        for row in rows:
            fields = tuple(row["fields"])
            enum = graphs.get(_canonical(sm, fields))
            if enum is None:                        # the four hand-only classes with identically zero monomials
                enum = enumerate_graphs(sm, fields, fams, props, insertions=1)
            mine = len(enum.topologies)
            assert mine == row["n_multiplet_labelled"] + EXTRA.get(fields, 0), (cls, fields)
            assert not enum.vanishing
            total += mine
    assert total == 2975 + sum(EXTRA.values())


def test_symmetry_factors(sm, setup, graphs):
    _, _, fams, _ = setup
    gg = {t.describe(fams): t.symmetry_factor for t in graphs[("G", "G")].topologies}
    assert sorted(gg.values()) == [Fraction(1, 2)] * 2 + [Fraction(1)] * 4
    assert all(v == Fraction(1, 2) for k, v in gg.items() if "G G G" in k)       # gauge bubble, four-gluon seagull
    hh = Counter(t.symmetry_factor for t in graphs[("H", "Hdag")].topologies)
    assert hh == Counter({Fraction(1): 6, Fraction(1, 2): 2})                      # report g0_SS: c_W, c_B are 1/2
    ssss = Counter(t.symmetry_factor for t in graphs[("H", "H", "Hdag", "Hdag")].topologies)
    assert ssss[Fraction(1, 2)] > 0 and set(ssss) <= {Fraction(1), Fraction(1, 2)}
    for t in graphs[("G", "G", "G", "G")].topologies:                               # |Aut| = number of automorphisms
        assert len(automorphisms(t)) == 1 / t.symmetry_factor


def test_canonical_key_is_a_graph_invariant(sm, graphs):
    """The key does not depend on the order of the vertices or of same-field legs; different graphs differ."""
    for fields in [("G", "G", "G", "G"), ("gh_G", "Q", "Qbar"), ("W", "W", "H", "Hdag")]:
        tops = graphs[fields].topologies
        assert len({t.key for t in tops}) == len(tops)
        for t in tops[:6]:
            nv = len(t.vertices)
            for perm in itertools.islice(itertools.permutations(range(nv)), 0, 24):
                pos = {old: new for new, old in enumerate(perm)}
                verts = [t.vertices[old] for old in perm]
                slot = {(v, j): (pos[v], j) for v in range(nv) for j in range(len(t.vertices[v][1]))}
                # exchange the first two same-field legs of every vertex where possible
                for v in range(nv):
                    cls = t.vertices[v][1]
                    for f in set(cls):
                        js = [j for j, g in enumerate(cls) if g == f]
                        if len(js) >= 2:
                            slot[(v, js[0])], slot[(v, js[1])] = (pos[v], js[1]), (pos[v], js[0])
                exts = [slot[s] for s in t.externals]
                lines = [(slot[a], slot[b], xy) for a, b, xy in t.lines]
                assert topology_key(verts, exts, lines) == t.key


def _rank(table, f):
    """Position class of an odd field in the report's monomials: omegabar A omega, rho omega, sigma^dag omega H,
    omega Psibar Psi (Delta_1), omega Psibar R, Psibar ... Psi, H^dag omega sigma."""
    fi = table[f]
    if fi.kind == "antighost":
        return 0
    if fi.is_source:
        return 4 if table[fi.source_of].kind == "scalar_conj" else 0
    return {"ghost": 1, "fermion_bar": 2, "fermion": 3}[fi.kind]


def _report_order_sign(table, fields) -> int:
    odd = [j for j, f in enumerate(fields) if table.odd(f)]
    perm = sorted(odd, key=lambda j: (_rank(table, fields[j]), j))
    inv = sum(1 for i in range(len(perm)) for k in range(i + 1, len(perm)) if perm[i] > perm[k])
    return -1 if inv % 2 else 1


def test_wick_signs_reproduce_the_hand_signs(sm, setup, graphs):
    """The global Wick rule of graphs.py, converted to the report's leg and vertex orders, gives the hand signs:
    closed fermion and ghost loops -1, open lines +1, two external ghosts +1 with the first ghost on Delta_1."""
    _, _, fams, props = setup
    table = sm.fields
    checked = Counter()
    for fields, enum in graphs.items():
        g = sum(table[f].ghost_number for f in fields)
        for t in enum.topologies:
            sign = record(t, fams, sm, props).grassmann_sign(table) * _report_order_sign(table, fields)
            for _, cls in t.vertices:
                sign *= _report_order_sign(table, cls)
            closed = all(table.odd(xy[0]) for _, _, xy in t.lines)
            ghosts = sorted((j for j, f in enumerate(fields) if table[f].kind == "ghost"),
                            key=lambda j: (_rank(table, fields[j]), j))
            if closed:
                assert sign == -1, (fields, t.describe(fams))
                checked["closed"] += 1
            elif g == 1 and len(ghosts) == 2:
                ins = next(v for v, (fi, _) in enumerate(t.vertices) if fams[fi].insertion)
                at = [t.externals[j][0] == ins for j in ghosts]
                if any(at):
                    assert sign == (1 if at[0] else -1), (fields, t.describe(fams))
                    checked["two ghosts"] += 1
            else:
                assert sign == 1, (fields, t.describe(fams))
                checked["open"] += 1
    assert checked == Counter({"closed": 2031, "open": 788, "two ghosts": 34}), checked


def test_two_identical_ghosts_are_antisymmetric(sm, setup, graphs):
    """omega^a omega^b Rbar Psi: exchanging the two external ghosts maps graphs onto graphs with opposite sign."""
    _, _, fams, props = setup
    table = sm.fields
    fields = ("gh_G", "gh_G", "Q", "src_Q")
    tops = {t.key: t for t in graphs[fields].topologies}
    assert len(tops) == 4
    for t in tops.values():
        swapped = topology_key(t.vertices, [t.externals[1], t.externals[0]] + list(t.externals[2:]), t.lines)
        partner = tops[swapped]
        assert record(t, fams, sm, props).grassmann_sign(table) == -record(partner, fams, sm, props).grassmann_sign(table)


def test_order_two_with_counterterm_vertices(sm, setup):
    """The order-l interface: a one-loop counterterm entry (loop order 1) enters the graphs of order 2. At one loop
    with one counterterm vertex: the gluon self-energy graphs with the two-point counterterm on a gluon line, or a
    three- or four-point counterterm vertex replacing a tree vertex; at order 1 the entry is absent."""
    brst, _, _, props = setup
    kin = gauge_kinetic(sm)
    k_g = Density.from_monomials(sm.fields, [m for m in kin.monomials() if set(m.fields()) == {"G"}])
    store = tree_store(sm, brst)
    store.add(StoreEntry("DS1", 1, Role.CT_FINITE, {0: LocalFunctional.integrate(k_g)}))
    assert [f.name for f in store_families(sm, store, 1)] == ["S0", "Delta1"]
    fams = store_families(sm, store, 2)
    assert [(f.name, f.loop_order) for f in fams] == [("S0", 0), ("Delta1", 0), ("DS1", 1)]
    assert set(fams[2].rules) == {("G", "G"), ("G", "G", "G"), ("G", "G", "G", "G")}
    enum = enumerate_graphs(sm, ("G", "G"), fams, props, order=2, loops=1)
    by_ct = Counter(next(c for f, c in t.vertices if fams[f].name == "DS1") for t in enum.topologies)
    assert by_ct == Counter({("G", "G"): 2, ("G", "G", "G"): 2, ("G", "G", "G", "G"): 1})
    assert all(t.order == 2 and t.loops == 1 for t in enum.topologies)
    sym = sorted(t.symmetry_factor for t in enum.topologies)
    assert sym == [Fraction(1, 2)] * 4 + [Fraction(1)]
    two_loop = enumerate_graphs(sm, ("G", "G"), fams, props, order=2, loops=2)
    assert all(fams[f].loop_order == 0 for t in two_loop.topologies for f, _ in t.vertices)
    assert len(two_loop.topologies) == 61                     # = HepKIT (test_hepkit_graphs.py)


def test_store_round_trip_keeps_families(sm, setup):
    _, store, fams, _ = setup
    again = OrderStore.from_json(sm, store.to_json())
    assert [(f.name, sorted(f.rules)) for f in store_families(sm, again, 1)] == [(f.name, sorted(f.rules))
                                                                                  for f in fams]
