"""Graph records (hepkit_renorm.graphs): routing, Wick sign, power counting, propagator splitting, forests."""

from __future__ import annotations

from fractions import Fraction

import pytest
from symbolica import Expression

from hepkit_renorm import rules, standard_legs
from hepkit_renorm import tensors as tn
from hepkit_renorm.brst import BRST
from hepkit_renorm.graphs import (
    LOOP_MOMENTUM,
    ExternalEnd,
    GraphPropagator,
    GraphRecord,
    GraphVertex,
    PropagatorEnd,
    split_denominators,
)
from hepkit_renorm.loops import divergent_subgraphs

E, Pe = ExternalEnd, PropagatorEnd


@pytest.fixture(scope="module")
def tree(sm):
    action = BRST(sm).action
    props = {p.fields: p for p in rules.propagators(action, sm)}
    return rules.tree_rules(action, sm), props


def _prop(sm, props, *fields):
    return GraphPropagator.from_rule(props[fields], {sm.gauge_fixing.xi: Expression.num(1)})


def _fermion_loop(sm, tree):
    r, props = tree
    legs = tuple(standard_legs(sm.fields, ["G", "G"]))
    v = r[("G", "Q", "Qbar")]
    return GraphRecord(legs, (GraphVertex.from_rule(v, [E(0), Pe(0, 0), Pe(1, 1)]),
                              GraphVertex.from_rule(v, [E(1), Pe(1, 0), Pe(0, 1)])),
                       (_prop(sm, props, "Q", "Qbar"),) * 2)


def test_split_denominators(sm, tree):
    _, props = tree
    h = GraphPropagator.from_rule(props[("H", "Hdag")])
    assert h.denominators == ((sm.param("m2"), 1),)
    assert (h.numerator - tn.I * tn.delta(*h.labels[0], *h.labels[1])).expand() == Expression.num(0)
    g = _prop(sm, props, "G", "G")
    assert g.denominators == ((Expression.num(0), 1),)          # Feynman gauge: one k^2 cancels
    general = GraphPropagator.from_rule(props[("G", "G")])
    assert general.denominators == ((Expression.num(0), 2),)
    num, dens = split_denominators(tn.I / (tn.spbar(0, 0) + tn.sphat(0, 0) - sm.param("m2")) ** 2)
    assert dens == ((sm.param("m2"), 2),) and num == tn.I


def test_routing_and_numbers(sm, tree):
    g = _fermion_loop(sm, tree)
    routing = g.routing()
    assert g.loops == 1 and g.is_one_particle_irreducible()
    assert routing[1] == {LOOP_MOMENTUM: 1}                    # the chord carries the loop momentum
    assert routing[0] == {2: 1, LOOP_MOMENTUM: 1}             # conservation at the vertex of leg 2
    for v in range(2):                                        # sum of incoming momenta vanishes up to sum p = 0
        tot: dict[int, Fraction] = {}
        for j in range(3):
            for n, c in g.leg_momentum(v, j, routing).items():
                tot[n] = tot.get(n, 0) + c
        assert all(c == 0 for n, c in tot.items() if n == LOOP_MOMENTUM)
    assert g.grassmann_sign(sm.fields) == -1                  # closed fermion loop
    assert g.superficial_degree(sm.fields) == 2
    assert divergent_subgraphs(g, sm.fields) == []


def test_wick_sign_open_line_and_order(sm, tree):
    r, props = tree
    legs = tuple(standard_legs(sm.fields, ["Qbar", "Q"]))
    v = r[("G", "Q", "Qbar")]
    se = GraphRecord(legs, (GraphVertex.from_rule(v, [Pe(1, 0), E(1), Pe(0, 1)]),
                            GraphVertex.from_rule(v, [Pe(1, 1), Pe(0, 0), E(0)])),
                     (_prop(sm, props, "Q", "Qbar"), _prop(sm, props, "G", "G")))
    assert se.grassmann_sign(sm.fields) == 1
    swapped = GraphRecord(tuple(reversed(legs)), tuple(
        GraphVertex(vt.legs, vt.value, tuple(E(1 - e.leg) if isinstance(e, ExternalEnd) else e for e in vt.ends))
        for vt in se.vertices), se.propagators)
    assert swapped.grassmann_sign(sm.fields) == -1           # exchanging two odd legs flips the sign
    assert se.superficial_degree(sm.fields) == 1


def test_validation(sm, tree):
    r, props = tree
    legs = tuple(standard_legs(sm.fields, ["G", "G"]))
    v = r[("G", "Q", "Qbar")]
    with pytest.raises(ValueError, match="unjoined"):
        GraphRecord(legs, (GraphVertex.from_rule(v, [E(0), Pe(0, 0), Pe(0, 1)]),), (_prop(sm, props, "Q", "Qbar"),))
    with pytest.raises(ValueError, match="joined to a"):
        GraphRecord(legs, (GraphVertex.from_rule(v, [E(0), Pe(0, 1), Pe(1, 1)]),
                           GraphVertex.from_rule(v, [E(1), Pe(1, 0), Pe(0, 0)])), (_prop(sm, props, "Q", "Qbar"),) * 2)


def test_divergent_subgraphs_two_loops(sm, tree):
    """A two-loop graph (gluon exchanged inside the fermion loop of the gluon two-point function): the inner fermion
    self-energy is a divergent proper subgraph."""
    r, props = tree
    legs = tuple(standard_legs(sm.fields, ["G", "G"]))
    v = r[("G", "Q", "Qbar")]
    q = _prop(sm, props, "Q", "Qbar")
    # fermion loop v0 -> v2 -> v3 -> v1 -> v0 with a gluon exchanged between v2 and v3
    g = GraphRecord(legs, (GraphVertex.from_rule(v, [E(0), Pe(0, 0), Pe(1, 1)]),
                           GraphVertex.from_rule(v, [E(1), Pe(3, 0), Pe(0, 1)]),
                           GraphVertex.from_rule(v, [Pe(4, 0), Pe(1, 0), Pe(2, 1)]),
                           GraphVertex.from_rule(v, [Pe(4, 1), Pe(2, 0), Pe(3, 1)])),
                    (q, q, q, q, _prop(sm, props, "G", "G")))
    assert g.loops == 2
    subs = divergent_subgraphs(g, sm.fields)
    assert frozenset({2, 4}) in subs                         # the self-energy of the inner fermion line
    assert all(len(s) < 5 for s in subs)
