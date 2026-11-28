"""The generator route (hepkit_renorm.hepkit_graphs): HepKIT on the written tree model, mapped back to the package's
multiplet-level graphs, against the direct enumeration (hepkit_renorm.diagrams).

HepKIT sees a multiplet-level graph iff some component assignment of its lines has non-zero vertex components; the
graphs only the direct enumeration has are group-theoretical zeros (checked to vanish here) or scaleless tadpoles
of the one-point function ``omega^B`` (which HepKIT does not generate). Over all 253 classes of ghost number 0 and
1 the two enumerations differ only by these six graphs (run in SCRATCHPAD/automation/ms27_diagrams/x07_hkcross.py;
here a selection covering every sector, about 25 s).
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from hepkit_renorm import LocalFunctional
from hepkit_renorm.action import gauge_kinetic
from hepkit_renorm.brst import BRST
from hepkit_renorm.class_results import ClassComputer
from hepkit_renorm.diagrams import enumerate_graphs, graph_propagators, store_families
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.functional import Density
from hepkit_renorm.hepkit_graphs import component_tuples, generated_keys
from hepkit_renorm.output import WriterOptions, tree_model, write_model
from hepkit_renorm.store import Role, StoreEntry, tree_store

CLASSES = [  # (fields in the report's leg order, graphs only in the direct enumeration)
    (("G", "G"), 0), (("W", "W"), 0), (("Qbar", "Q"), 0), (("Hdag", "H"), 0), (("W", "W", "W"), 0),
    (("W", "H", "Hdag"), 1), (("Qbar", "D", "H"), 0), (("src_G", "gh_G"), 0), (("src_W", "gh_W", "W"), 0),
    (("gh_W", "Hdag", "src_Hdag"), 0), (("gh_G", "G"), 0), (("gh_G", "Qbar", "Q"), 0), (("gh_W", "Qbar", "Q"), 0),
    (("gh_B", "Qbar", "Q"), 0), (("gh_G", "gh_G", "src_Q", "Q"), 0), (("gh_W", "gh_B", "Qbar", "src_Qbar"), 0),
    (("gh_W", "Hdag", "H", "W"), 0), (("gh_B",), 5),
]


@pytest.fixture(scope="module")
def setup(sm):
    brst = BRST(sm)
    store = tree_store(sm, brst)
    return brst, store, store_families(sm, store, 1), graph_propagators(sm, store), tree_model(sm, brst=brst)


def test_component_tuples(sm, setup):
    written = setup[4]
    assert component_tuples(written, ("G", "G")) == [("G", "G")]
    assert len(component_tuples(written, ("W", "Q", "Qbar"))) == 12
    assert component_tuples(written, ("src_gh_B", "gh_B")) == []      # zeta_B has no vertex and is not written


@pytest.mark.parametrize("fields, extra", CLASSES, ids=[" ".join(c[0]) for c in CLASSES])
def test_generator_route_equals_direct_enumeration(sm, setup, fields, extra):
    _, _, fams, props, written = setup
    g = sum(sm.fields[f].ghost_number for f in fields)
    mine = {t.key: t for t in enumerate_graphs(sm, fields, fams, props, insertions=g).topologies}
    hk = generated_keys(written, sm, fields, fams, props, insertions=g)
    assert set(hk) <= set(mine)
    assert len(mine) - len(hk) == extra
    for key, diagrams in hk.items():            # component symmetry factors = multiplet ones when no line has components
        plain = [d for d in diagrams if all(not written.reverse_map.particles[p].components for p in d.lines)]
        assert all(d.symmetry_factor == mine[key].symmetry_factor for d in plain)


def test_graphs_missing_in_hepkit_vanish(sm, setup):
    """The group-theoretical zero of W H Hdag (W bubble on the W W H Hdag seagull: f^{abc} {T^b, T^c} = 0) has no
    pole: its Laurent data vanish on D = 6 samples."""
    _, store, fams, props, written = setup
    fields = ("W", "H", "Hdag")
    cc = ClassComputer(sm, store, families=fams, propagators=props)
    enum, recs = cc.records(fields)
    hk = generated_keys(written, sm, fields, fams, props)
    missing = [r for t, r in zip(enum.topologies, recs) if t.key not in hk]
    assert len(missing) == 1
    fp = Fingerprinter(sm, dimension=6, seed=3)
    smp = fp.samples(list(fields), 3)
    lau = cc.backend.evaluate(cc.backend.prepare(missing[0], sm), fp, smp)
    assert all(not np.asarray(c).any() for c in lau.coefficients.values())


def test_two_loop_sets_and_symmetry_factors(sm, setup):
    """The enumeration is not one-loop specific: the two-loop graphs of G G, Qbar Q and src_G gh_G equal HepKIT's,
    with the same symmetry factors where no line carries components."""
    _, store, _, props, written = setup
    fams = store_families(sm, store, 2)
    plain = 0
    for fields, count in ((("G", "G"), 61), (("Qbar", "Q"), 104), (("src_G", "gh_G"), 9)):
        mine = {t.key: t for t in enumerate_graphs(sm, fields, fams, props, order=2).topologies}
        hk = generated_keys(written, sm, fields, fams, props, loops=2)
        assert set(hk) == set(mine) and len(mine) == count
        for key, diagrams in hk.items():
            for d in diagrams:
                if all(not written.reverse_map.particles[p].components for p in d.lines):
                    assert d.symmetry_factor == mine[key].symmetry_factor
                    plain += 1
    assert plain >= 30
    gg = [t.symmetry_factor for t in enumerate_graphs(sm, ("G", "G"), fams, props, order=2).topologies]
    assert sorted(set(gg)) == [Fraction(1, 6), Fraction(1, 4), Fraction(1, 2), Fraction(1)]


def test_counterterm_insertions_against_hepkit(sm, setup):
    """Order 2, one loop, one counterterm vertex of loop order 1 (written with order CT = 1): the same graphs."""
    brst, _, _, props, _ = setup
    kin = gauge_kinetic(sm)
    k_g = Density.from_monomials(sm.fields, [m for m in kin.monomials() if set(m.fields()) == {"G"}])
    store = tree_store(sm, brst)
    store.add(StoreEntry("DS1", 1, Role.CT_FINITE, {0: LocalFunctional.integrate(k_g)}))
    written = write_model(sm, store, WriterOptions(substitutions={"xi": 1}))
    fams = store_families(sm, store, 2)
    for fields in (("G", "G"), ("G", "G", "G")):
        mine = {t.key for t in enumerate_graphs(sm, fields, fams, props, order=2, loops=1).topologies}
        hk = generated_keys(written, sm, fields, fams, props, loops=1, order=2)
        assert set(hk) == mine and mine
