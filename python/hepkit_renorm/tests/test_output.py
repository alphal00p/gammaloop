"""Gate V8 (first part): the model writer and the order store.

* The tree-level regularised model (``S_0`` with gauge fixing, ghosts and sources, and ``Delta_1``) is written from
  the specification, loads in the installed HepKIT, and every vertex rule and propagator, read back from its UFO
  strings by an independent converter in this file, equals the package's rule through the reverse map (random
  kinematics, parameters and polarisations; SU(2) components fixed), and equals the hand-written model
  ``tests/data/sm_unbroken_bmhv.json`` (origin ``SCRATCHPAD/model/build_model.py``) rule by rule, up to the
  Grassmann sign of the different leg orders.
* HepKIT generates the same one-loop diagram sets from the written model as from the hand-written one (vertices and
  internal lines mapped to package fields and components, symmetry factors), with the counts of
  ``SCRATCHPAD/model/ENUMERATION.md`` (per component tuple, from ``SCRATCHPAD/model/parts/<class>.json``).
* Counterterm entries (singular with Laurent coefficients, the evanescent tag and ``d``; finite barred ones) are
  written with their orders, roles and v2 metadata, and their values per power agree with the package.
"""

from __future__ import annotations

import itertools
import json
from collections import Counter

import numpy as np
import pytest
from component_values import ComponentEvaluator, Kinematics, model_parameters, unit_vector
from conftest import DATA
from symbolica import AtomType, Expression, S
from test_rules import hand_particle, written_order

from hepkit_renorm import tensors as tn
from hepkit_renorm.action import gauge_kinetic
from hepkit_renorm.brst import BRST
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.functional import Density, LocalFunctional, standard_legs
from hepkit_renorm.output import ReverseMap, WriterOptions, component_samples, tree_model, write_model
from hepkit_renorm.rules import propagators
from hepkit_renorm.store import EVANESCENT_TAG, OrderStore, Role, StoreEntry, tree_store

I = tn.I
EPS = S("test_output::eps")
HAND_PARAMS = {"gs": ("coupling", "G"), "gw": ("coupling", "W"), "gp": ("coupling", "B"), "yu": ("param", "y_u"),
               "yd": ("param", "y_d"), "ye": ("param", "y_e"), "lam": ("param", "lambda")}


# ---------------------------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def brst(sm):
    return BRST(sm)


@pytest.fixture(scope="module")
def written(sm, brst):
    return tree_model(sm, brst=brst)


@pytest.fixture(scope="module")
def hand():
    with open(DATA / "sm_unbroken_bmhv.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def ev(sm, field):
    return ComponentEvaluator(sm, field)


@pytest.fixture(scope="module")
def hk():
    from symbolica.community import hepkit

    return hepkit


# ---------------------------------------------------------------------------------------------- UFO -> package
def _short(e: Expression) -> str:
    return e.get_name().split("::")[-1]


def ufo_expression(text: str) -> Expression:
    s = text.replace("UFO::{}::", "").replace("1\U0001d456", "hkI").replace("\U0001d456", "hkI")
    return Expression.parse(s)


class UfoConverter:
    """Independent reader of the UFO strings: ``idx(1,k)`` / positive colour indices -> labels of leg ``k``,
    ``dummy(n)`` / negative colour indices -> fresh labels, ``bmhv::bar(x)`` / ``bmhv::hat(x)`` -> every occurrence
    of ``x`` projected (gammabar, gbar, p.gbar; gammahat, ghat, p.ghat), ``P(x)`` -> ``p(0, x)``."""

    def __init__(self, sm, legs, params, colour_factor="G"):
        self.sm, self.legs, self.params = sm, legs, params
        fac = sm.factors[colour_factor]
        self.spaces = {"lorentz": tn.LORENTZ, "spinor": tn.SPINOR, "adj": fac.adjoint_space(),
                       "fund": fac.rep_space("F")}
        self.colour_factor = colour_factor

    def leg_label(self, k: int, space: str) -> Expression:
        leg = self.legs[k - 1]
        for lab, s in zip(leg.labels, self.sm.fields[leg.field].slots):
            kind = s.space.kind
            if space in ("lorentz", "spinor") and kind == space:
                return lab
            if s.factor == self.colour_factor and space in (kind if kind != "adjoint" else "adj", "colour"):
                return lab
        raise KeyError((leg.field, space))

    def convert(self, text: str, colour: bool = False) -> Expression:
        """The package expression of a UFO string (markers apply to the product term they stand in)."""
        e = ufo_expression(text).expand()
        out = Expression.num(0)
        for term in (list(e) if e.get_type() == AtomType.Add else [e]):
            bar, hat = set(), set()
            kept = Expression.num(1)
            for fa in (list(term) if term.get_type() == AtomType.Mul else [term]):
                if fa.get_type() == AtomType.Fn and fa.get_name().startswith("bmhv::"):
                    (bar if _short(fa) == "bar" else hat).add(next(iter(fa)).to_canonical_string())
                else:
                    kept = kept * fa
            self.dummies: dict[str, Expression] = {}
            self.bar, self.hat, self.colour = bar, hat, colour
            out = out + self._conv(kept)
        return out.expand()

    def _lab(self, arg: Expression, space: str) -> Expression:
        key = arg.to_canonical_string()
        if arg.get_type() == AtomType.Fn and _short(arg) == "idx":
            return self.leg_label(int(list(arg)[1].to_canonical_string()), space)
        if arg.get_type() == AtomType.Num and int(key) > 0:
            return self.leg_label(int(key), space)
        if key not in self.dummies:
            self.dummies[key] = tn.fresh(self.spaces[space])
        return self.dummies[key]

    def _proj(self, arg: Expression):
        key = arg.to_canonical_string()
        if key in self.bar:
            return tn.gbar
        if key in self.hat:
            return tn.ghat
        return None

    def _conv(self, e: Expression) -> Expression:
        t = e.get_type()
        if t == AtomType.Add:
            return sum((self._conv(a) for a in e), Expression.num(0))
        if t == AtomType.Mul:
            out = Expression.num(1)
            for a in e:
                out = out * self._conv(a)
            return out
        if t == AtomType.Pow:
            b, x = list(e)
            if b.get_type() == AtomType.Fn and _short(b) == "P" and x == Expression.num(2):
                return tn.sp(0, 0)
            return self._conv(b) ** x
        if t == AtomType.Num:
            return e
        if t == AtomType.Var:
            name = _short(e)
            return I if name == "hkI" else self.params[name]
        name, args = _short(e), list(e)
        if self.colour:
            if name == "T":
                return tn.T(self._lab(args[0], "adj"), self._lab(args[2], "fund"), self._lab(args[1], "fund"))
            if name == "f":
                return tn.f(*(self._lab(a, "adj") for a in args))
            if name == "d":
                return tn.dsym(*(self._lab(a, "adj") for a in args))
            if name == "Identity":
                return tn.delta(*(self._lab(a, "colour") for a in args))
            raise ValueError(name)

        def lor(arg):
            lab = self._lab(arg, "lorentz")
            proj = self._proj(arg)
            if proj is None:
                return lab, Expression.num(1)
            new = tn.fresh(tn.LORENTZ)
            return new, proj(new, lab)

        if name == "Gamma":
            mu = self._lab(args[0], "lorentz")
            key = args[0].to_canonical_string()
            head = tn.gammabar if key in self.bar else tn.gammahat if key in self.hat else tn.gamma
            return head(mu, self._lab(args[1], "spinor"), self._lab(args[2], "spinor"))
        if name == "Metric":
            (a, pa), (b, pb) = lor(args[0]), lor(args[1])
            return tn.g(a, b) * pa * pb
        if name == "P":
            mu, pm = lor(args[0])
            k = 0 if len(args) == 1 else self.legs[int(args[1].to_canonical_string()) - 1].momentum
            return tn.p(k, mu) * pm
        if name == "PSlash":
            mu = tn.fresh(tn.LORENTZ)
            return tn.gamma(mu, self._lab(args[0], "spinor"), self._lab(args[1], "spinor")) * tn.p(0, mu)
        if name in ("ProjM", "ProjP", "Gamma5", "Identity"):
            head = {"ProjM": tn.PL, "ProjP": tn.PR, "Gamma5": tn.gamma5, "Identity": tn.delta}[name]
            return head(self._lab(args[0], "spinor"), self._lab(args[1], "spinor"))
        raise ValueError(f"unknown UFO head {name}")


def json_rule_expression(sm, rule: dict, couplings: dict, lorentz: dict, conv: UfoConverter) -> Expression:
    """sum_{c, l} coupling * colour * lorentz of a JSON vertex rule, in the package vocabulary."""
    total = Expression.num(0)
    for c, row in zip(rule["color_structures"], rule["couplings"]):
        colour = conv.convert(c, colour=True)
        for lname, cname in zip(rule["lorentz_structures"], row):
            if cname is None:
                continue
            coupling = conv.convert(couplings[cname]["expression"])
            total = total + colour * conv.convert(lorentz[lname]["structure"]) * coupling
    return total.expand()


def ufo_params(sm, rmap: ReverseMap) -> dict[str, Expression]:
    out = {name: rmap.symbol(sm, name) for name, (kind, _) in rmap.parameters.items()
           if kind in ("coupling", "parameter", "gauge_parameter")}
    out.update({"eps": EPS, "eta": EVANESCENT_TAG, "pi": Expression.PI, _short(Expression.PI): Expression.PI})
    return out


def polarisations(sm, ev, kin, legs, comps):
    """Random polarisations; the component slots of leg ``j`` on the basis vectors ``comps[j]``."""
    pol = {}
    for leg, cj in zip(legs, comps):
        fixed = dict(cj)
        for pos, (lab, _slot) in enumerate(zip(leg.labels, sm.fields[leg.field].slots)):
            n = ev.dim(lab)
            pol[lab.to_canonical_string()] = unit_vector(ev.F, n, fixed[pos]) if pos in fixed else kin.vector(n)
    return pol


def same_values(sm, ev, field, a, b, legs, comps, extra=(), seeds=(0, 1)) -> bool:
    names = model_parameters(sm) + list(extra)
    for seed in seeds:
        kin = Kinematics(field, len(legs), names, seed)
        pol = polarisations(sm, ev, kin, legs, comps)
        if ev.value(a, pol, kin.momenta, kin.params) != ev.value(b, pol, kin.momenta, kin.params):
            return False
    return True


def fraction_value(ev, expr: Expression, pol, kin) -> int:
    """Value of a rational function of momenta (numerator / denominator evaluated separately)."""
    e = expr.together()
    num, den = Expression.num(1), Expression.num(1)
    for fa in (list(e) if e.get_type() == AtomType.Mul else [e]):
        if fa.get_type() == AtomType.Pow and list(fa)[1].get_type() == AtomType.Num and \
                int(list(fa)[1].to_canonical_string()) < 0:
            den = den * next(iter(fa)) ** (-int(list(fa)[1].to_canonical_string()))
        else:
            num = num * fa
    F = ev.F
    return ev.value(num.expand(), pol, kin.momenta, kin.params) * F.inv(
        ev.value(den.expand(), pol, kin.momenta, kin.params)) % F.p


def by_name(items):
    return {x["name"]: x for x in items}


# ---------------------------------------------------------------------------------------------- store
def test_tree_store(sm, brst):
    store = tree_store(sm, brst)
    assert [(e.name, e.loop_order, e.role) for e in store] == [("S0", 0, Role.TREE), ("Delta1", 0, Role.BREAKING)]
    assert store.tree_action() == brst.action
    back = OrderStore.from_json(sm, store.to_json())
    assert [e.name for e in back] == ["S0", "Delta1"]
    assert back.entry("S0").coefficient(0) == brst.action and back.entry("Delta1").coefficient(0) == brst.breaking
    assert store.entries(roles=[Role.BREAKING])[0].name == "Delta1"


def test_store_validation(sm, brst):
    s0 = brst.action
    with pytest.raises(ValueError):
        StoreEntry("x", 1, Role.TREE, {0: s0})
    with pytest.raises(ValueError):
        StoreEntry("x", 0, Role.CT_SINGULAR, {-1: s0})
    with pytest.raises(ValueError):
        StoreEntry("x", 0, Role.TREE, {-1: s0})
    store = tree_store(sm, brst)
    with pytest.raises(ValueError):
        store.add(StoreEntry("S0", 0, Role.TREE, {0: s0}))


# ---------------------------------------------------------------------------------------------- tree model
def test_tree_model_content(sm, written, hand):
    d = written.data
    rmap = written.reverse_map
    assert len(d["vertex_rules"]) == len(hand["vertex_rules"]) == 174
    # the hand-written model also lists zeta_B (no vertex, f^{Bbc} = 0); fields without vertices are not written
    assert len(d["particles"]) == len(hand["particles"]) - 1
    assert {o["name"] for o in d["orders"]} == {"GG", "GW", "GB", "YUK", "POT", "DELTA", "SRC"}
    assert {p["particle"] for p in d["propagators"]} == {p["name"] for p in d["particles"] if p["propagating"]}
    assert all(p["ghost_number"] == 0 for p in d["particles"])
    assert all(not p["propagating"] for p in d["particles"] if p["field_kind"] == "brst_source")
    names = {p["name"] for p in d["particles"]}
    assert all(by_name(d["particles"])[p["antiname"]]["antiname"] == p["name"] for p in d["particles"])
    assert set(rmap.particles) == names
    # every vertex rule maps back to package legs, and the particle multisets are those of the hand model
    ours = Counter()
    for v in d["vertex_rules"]:
        rec = rmap.vertex(v["name"])
        assert list(rec.particles) == v["particles"]
        legs = []
        for p, f, comps in zip(v["particles"], rec.fields, rec.components):
            assert rmap.particle(p) == (f, dict(comps))
            legs.append((f, comps[0][1] if comps else None))
        ours[tuple(sorted(legs, key=repr))] += 1
        assert v["role"] in ("tree", "source", "breaking")
        assert (v["role"] == "breaking") == ("DELTA" in dict(rec.orders))
    theirs = Counter(tuple(sorted((hand_particle(p) for p in r["particles"]), key=repr)) for r in hand["vertex_rules"])
    assert ours == theirs
    assert sum(v["evanescent"] for v in d["vertex_rules"]) == 23  # the Delta_1 rules (eq:V-delta)
    assert ReverseMap.from_dict(d["hepkit_renorm"]["reverse_map"]).to_dict() == rmap.to_dict()


def test_tree_model_loads(written, hk):
    m = hk.Model.from_json(written.json_string())
    assert len(m.particles) == len(written.data["particles"])
    assert len(m.vertex_rules) == 174
    again = json.loads(m.to_json())
    assert {v["name"] for v in again["vertex_rules"]} == {v["name"] for v in written.data["vertex_rules"]}


def _compare_with_package(sm, written, ev, field, mutate=None):
    d = written.data
    couplings, lorentz = by_name(d["couplings"]), by_name(d["lorentz_structures"])
    params = ufo_params(sm, written.reverse_map)
    bad = []
    for v in d["vertex_rules"]:
        rec = written.reverse_map.vertex(v["name"])
        rule = written.package_rule(v["name"])
        conv = UfoConverter(sm, rule.legs, params)
        ours = json_rule_expression(sm, v, couplings, lorentz, conv)
        if mutate is not None:
            ours = mutate(v, ours)
        if not same_values(sm, ev, field, ours, (I * rule.value).expand(), rule.legs, rec.components):
            bad.append(v["name"])
    return bad


def test_tree_rules_equal_package_rules(sm, written, ev, field):
    """Every written vertex rule, read back from its strings, is i V of the package at its components."""
    assert _compare_with_package(sm, written, ev, field) == []


def test_tree_rules_negative_control(sm, written, ev, field):
    """The comparison sees a sign error in one rule and a transposed colour generator."""
    flipped = _compare_with_package(sm, written, ev, field,
                                    mutate=lambda v, e: -e if v["name"] == "S0_G_Q1_Q1bar" else e)
    assert flipped == ["S0_G_Q1_Q1bar"]
    a, b, c = S("ca_", "cb_", "cc_")

    def transpose(v, e):
        return e.replace(tn.T(a, b, c), tn.T(a, c, b))
    wrong = _compare_with_package(sm, written, ev, field, mutate=transpose)
    assert "S0_G_Q1_Q1bar" in wrong and all("Q" in n or "U" in n or "D" in n for n in wrong)


def test_tree_rules_equal_hand_written_rules(sm, written, hand, ev, field):
    """Each of the 174 hand-written rules equals the written rule with the same package legs: the hand value is
    i V in the report's written order of its legs, ours i V in the canonical leg order; they differ by the sign of
    the permutation of the odd legs."""
    d = written.data
    couplings, lorentz = by_name(d["couplings"]), by_name(d["lorentz_structures"])
    params = ufo_params(sm, written.reverse_map)
    hand_params = {k: (sm.coupling(r) if kind == "coupling" else sm.param(r)) for k, (kind, r) in HAND_PARAMS.items()}
    hand_coup = by_name(hand["couplings"])
    hand_lor = by_name(hand["lorentz_structures"])
    index = {}
    for v in d["vertex_rules"]:
        rec = written.reverse_map.vertex(v["name"])
        key = tuple(sorted(((f, c[0][1] if c else None) for f, c in zip(rec.fields, rec.components)), key=repr))
        index[key] = v
    bad = []
    for h in hand["vertex_rules"]:
        hlegs = [hand_particle(p) for p in h["particles"]]
        v = index[tuple(sorted(hlegs, key=repr))]
        rec = written.reverse_map.vertex(v["name"])
        legs = standard_legs(sm.fields, rec.fields)
        mine = [(f, c[0][1] if c else None) for f, c in zip(rec.fields, rec.components)]
        free = list(range(len(mine)))
        perm = []  # hand position -> our position
        for hl in hlegs:
            j = next(j for j in free if mine[j] == hl)
            free.remove(j)
            perm.append(j)
        hand_legs = [legs[j] for j in perm]
        order = [perm[i] for i in written_order(sm, [f for f, _ in hlegs])]
        odd = [j for j in order if sm.fields.odd(rec.fields[j])]
        sign = -1 if sum(1 for x, y in itertools.combinations(odd, 2) if x > y) % 2 else 1
        theirs = json_rule_expression(sm, h, hand_coup, hand_lor, UfoConverter(sm, hand_legs, hand_params))
        ours = json_rule_expression(sm, v, couplings, lorentz, UfoConverter(sm, legs, params))
        if not same_values(sm, ev, field, theirs, (sign * ours).expand(), legs, rec.components):
            bad.append(h["name"])
    assert bad == []


def test_propagators(sm, written, hand, brst, ev, field):
    """Written propagators equal the package's D (xi = 1) through the reverse map, and the hand-written ones."""
    d = written.data
    params = ufo_params(sm, written.reverse_map)
    package = {p.fields: p for p in propagators(brst.action, sm, sm.simplifier(expand_projectors=False))}
    hand_props = by_name(hand["propagators"])
    hand_names = {}
    for p in hand["particles"]:
        f, c = hand_particle(p["name"])
        hand_names[(f, c)] = p["name"]
    xi = sm.gauge_fixing.xi
    for p in d["propagators"]:
        rec = written.reverse_map.propagators[p["name"]]
        x, y = rec.fields
        prop = package[(x, y)]
        lx, ly = prop.labels
        ends = {1: ly, 2: lx} if not rec.reversed else {1: lx, 2: ly}
        legs = [type("L", (), {"field": y if not rec.reversed else x, "labels": ends[1]})(),
                type("L", (), {"field": x if not rec.reversed else y, "labels": ends[2]})()]
        conv = UfoConverter(sm, legs, params)
        mine = (conv.convert(p["numerator"]) / conv.convert(p["denominator"]))
        deltas = Expression.num(1)
        for la, lb, sl in zip(lx, ly, sm.fields[x].slots):
            if sl.space.kind not in ("lorentz", "spinor"):
                deltas = deltas * tn.delta(la, lb)
        a_, b_ = S("a_", "b_")
        target = prop.value.replace(xi, 1).replace(tn.sphat(a_, b_), tn.sp(a_, b_) - tn.spbar(a_, b_))
        target = target.replace(tn.gammahat(a_, b_, S("c_")), tn.gamma(a_, b_, S("c_")) - tn.gammabar(a_, b_, S("c_")))
        if rec.reversed:
            target = target.replace(tn.p(0, a_), -tn.p(0, a_))
        # the hand-written propagator of the same particle (M_H^2 = m2)
        field_, comps = written.reverse_map.particle(p["particle"])
        hname = hand_names[(field_, next(iter(comps.values())) if comps else None)]
        hp = hand_props[hname + "_prop"]
        hconv = UfoConverter(sm, legs, {"MH": sm.param("m2") ** (Expression.num(1) / 2)})
        theirs = hconv.convert(hp["numerator"]) / hconv.convert(hp["denominator"])
        names = model_parameters(sm)
        for seed in (0, 1):
            kin = Kinematics(field, 1, names, seed)
            pol = {lab.to_canonical_string(): kin.vector(ev.dim(lab)) for lab in list(lx) + list(ly)}
            ours = fraction_value(ev, mine * deltas, pol, kin)
            assert ours == fraction_value(ev, target, pol, kin), p["name"]
            assert fraction_value(ev, mine, pol, kin) == fraction_value(ev, theirs, pol, kin), p["name"]


# ---------------------------------------------------------------------------------------------- diagrams
# (hand incoming, hand outgoing, number of one-loop 1PI diagrams) for component tuples of SCRATCHPAD/model/parts
# (ENUMERATION.md: g0_VV, g0_VVV, g0_FF, g0_FFS, g0_ghost_sources, g1_cFF, g1_cFFV, g1_cSS_cSSV, g1_ccRF).
ENUMERATED = [
    (["G", "G"], [], 7), (["W1", "W1"], [], 13), (["W1", "W2"], [], 6), (["B", "B"], [], 11),
    (["W1", "W2", "W3"], [], 18), (["W1", "W1", "W2"], [], 0),
    (["Q1"], ["Q1"], 7), (["Q2"], ["Q1"], 0), (["U"], ["Q1", "H2"], 8), (["U"], ["Q1", "H1"], 0),
    (["rhoG", "ghG"], [], 1), (["ghG", "Q1"], ["RQ1"], 1), (["ghG", "Q2"], ["RQ1"], 0),
    (["ghG", "Q1"], ["Q1"], 9), (["ghG", "Q1", "G"], ["Q1"], 21),
    (["ghW1", "H2"], ["H1"], 3), (["ghW1", "H1"], ["H1"], 0),
    (["ghG", "ghG", "Q1"], ["RQ1"], 4),
]


def generate(model, inc, out, n_breaking, n_sources, source_pdgs):
    """One-loop 1PI diagrams with the options of SCRATCHPAD/model/enumerate_diagrams.py."""
    def no_internal_sources(graph, _n):
        nodes = graph.nodes()
        for u, v, _d, pdg in graph.edges():
            if abs(int(str(pdg))) in source_pdgs and int(str(nodes[u][1])) == 0 and int(str(nodes[v][1])) == 0:
                return False
        return True
    kw = {"loops": 1, "coupling_orders": {"DELTA": n_breaking, "SRC": n_sources}, "max_vertices": len(inc) + len(out),
              "self_energy": None, "tadpoles": None, "zero_snails": None, "progress": None}
    if n_sources:
        kw["filter"] = no_internal_sources
    return model.process(inc, out).generate_diagrams(**kw)


def test_diagram_sets_equal_hand_model(sm, written, hand, hk):
    """Same diagrams (package vertex legs, internal lines, symmetry factors) from both models, with the counts of
    the enumeration of the report."""
    ours_model, hand_model = written.load(), hk.Model.from_json(json.dumps(hand))
    rmap = written.reverse_map
    by_package = {(r.field, r.components[0][1] if r.components else None): name
                  for name, r in rmap.particles.items()}
    hand_rules = {r["name"]: r["particles"] for r in hand["vertex_rules"]}
    hand_src = {abs(p["pdg_code"]) for p in hand["particles"] if hand_particle(p["name"])[0].startswith("src_")}
    our_src = {abs(r.pdg_code) for r in rmap.particles.values() if sm.fields[r.field].is_source}

    def ours_of(name):
        f, c = rmap.particle(name)
        return f, next(iter(c.values())) if c else None

    def key(dg, legs_of, particle):
        vs = tuple(sorted(tuple(sorted((particle(p) for p in legs_of(v.interaction)), key=repr)) for v in dg.vertices))
        es = tuple(sorted((particle(e.particle_name.rstrip("~")) for e in dg.internal_edges), key=repr))
        return vs, es, str(dg.symmetry_factor)

    for inc, out, count in ENUMERATED:
        fields = [hand_particle(n) for n in inc + out]
        n_breaking = sum(sm.fields[f].ghost_number for f, _ in fields)  # one Delta_1 at ghost number one
        n_sources = sum(1 for f, _ in fields if sm.fields[f].is_source)
        mine = [by_package[hand_particle(n)] for n in inc], [by_package[hand_particle(n)] for n in out]
        theirs = generate(hand_model, inc, out, n_breaking, n_sources, hand_src)
        ours = generate(ours_model, *mine, n_breaking, n_sources, our_src)
        assert len(theirs) == len(ours) == count, (inc, out)
        assert Counter(key(dg, lambda n: hand_rules[n], hand_particle) for dg in theirs) == \
            Counter(key(dg, lambda n: rmap.vertex(n).particles, ours_of) for dg in ours), (inc, out)


# ---------------------------------------------------------------------------------------------- counterterms
@pytest.fixture(scope="module")
def ct_store(sm, brst):
    """S_0, Delta_1 and two synthetic one-loop counterterm entries built from the gluon kinetic term: a singular
    one ``(c d / eps + c' eta) K_G`` and a finite barred one ``c'' P[K_G]`` (stand-ins for deltaS^(1), DeltaS^(1))."""
    fields = sm.fields
    kin = gauge_kinetic(sm)
    k_g = Density.from_monomials(fields, [m for m in kin.monomials() if set(m.fields()) == {"G"}])
    gs = sm.coupling("G")
    c1 = gs**2 * Expression.num(5) / 3 / (16 * Expression.PI**2)
    c0 = gs**2 * EVANESCENT_TAG / 7
    singular = StoreEntry("dS1", 1, Role.CT_SINGULAR, {
        -1: LocalFunctional.integrate(k_g.scale(c1 * tn.d)), 0: LocalFunctional.integrate(k_g.scale(c0))})
    finite = StoreEntry("DS1", 1, Role.CT_FINITE, {
        0: LocalFunctional.integrate(k_g.physical(sm.simplifier()).scale(gs**2 / 11))})
    # vanishes at d = 6 only: the zero test must not drop it (it fingerprints in D = 4, 6 and 8)
    at_six = StoreEntry("dS1d", 1, Role.CT_SINGULAR, {0: LocalFunctional.integrate(k_g.scale(gs**2 * (tn.d - 6) / 13))})
    store = tree_store(sm, brst)
    for entry in (singular, finite, at_six):
        store.add(entry)
    return store


@pytest.fixture(scope="module")
def ct_written(sm, ct_store):
    return write_model(sm, ct_store, WriterOptions(substitutions={"xi": 1}))


def test_counterterm_encoding(sm, ct_written, ct_store):
    d = ct_written.data
    at_six = [v for v in d["vertex_rules"] if v["operator"] == "dS1d"]
    assert len(at_six) == 3
    lau = by_name(d["couplings"])[at_six[0]["couplings"][0][0]]["laurent"]["terms"]
    assert sorted(int(k) for k in lau) == [0, 1]  # (d - 6) = -2 - 2 eps
    rules = [v for v in d["vertex_rules"] if v["operator"] in ("dS1", "DS1")]
    assert {(v["operator"], tuple(v["particles"])) for v in rules} == {
        (op, ("G",) * n) for op in ("dS1", "DS1") for n in (2, 3, 4)}
    couplings = by_name(d["couplings"])
    for v in rules:
        rec = ct_written.reverse_map.vertex(v["name"])
        orders = dict(rec.orders)
        assert orders["CT"] == 1 and orders["GG"] == 2 + len(v["particles"]) - 2
        assert ("CTPOLE" in orders) == (v["operator"] == "dS1") and ("CTFIN" in orders) == (v["operator"] == "DS1")
        assert v["loop_order"] == 1 and v["role"] == ("ct_singular" if v["operator"] == "dS1" else "ct_finite")
        for row in v["couplings"]:
            for c in row:
                if c is None:
                    continue
                lau = couplings[c]["laurent"]
                assert lau["variable"] == "eps"
                powers = sorted(int(k) for k in lau["terms"])
                # d / eps = 4 / eps - 2: the singular entry has eps^-1 and eps^0 (with the tag), the finite eps^0
                assert powers == ([-1, 0] if v["operator"] == "dS1" else [0])
                if v["operator"] == "dS1":
                    assert "eta" in lau["terms"]["0"] and "pi" in lau["terms"]["-1"]
    finite_lor = {lz for v in rules if v["operator"] == "DS1" for lz in v["lorentz_structures"]}
    structures = by_name(d["lorentz_structures"])
    assert all("bmhv::bar(" in structures[lz]["structure"] for lz in finite_lor)
    assert all("Bar(" in structures[lz]["v2_overrides"]["structure"] for lz in finite_lor)
    ops = by_name(d["operators"])
    assert ops["dS1"]["laurent_powers"] == [-1, 0] and ops["dS1"]["role"] == "ct_singular"
    kinds = {o["name"]: o["kind"] for o in d["orders"]}
    assert kinds["CT"] == "loop_counterterm" and kinds["DELTA"] == "insertion" and kinds["GG"] == "coupling"
    # the store with the counterterms survives a JSON round trip
    back = OrderStore.from_json(sm, ct_store.to_json())
    assert back.entry("dS1").coefficient(-1) == ct_store.entry("dS1").coefficient(-1)
    # max_loop_order selects the entries of lower orders (step 1 of the order-l algorithm)
    tree_only = write_model(sm, ct_store, WriterOptions(substitutions={"xi": 1}, max_loop_order=0))
    assert len(tree_only.data["vertex_rules"]) == 174


def test_counterterm_values(sm, ct_written, ev, field):
    """Written counterterm couplings (sum_k eps^k c_k, d = 4 - 2 eps) equal the package's entries."""
    d = ct_written.data
    couplings, lorentz = by_name(d["couplings"]), by_name(d["lorentz_structures"])
    params = ufo_params(sm, ct_written.reverse_map)
    for v in d["vertex_rules"]:
        if v["operator"] not in ("dS1", "DS1", "dS1d"):
            continue
        rec = ct_written.reverse_map.vertex(v["name"])
        target = Expression.num(0)
        legs = None
        for k in rec.powers:
            rule = ct_written.package_rule(v["name"], k)
            legs = rule.legs
            target = target + (I * rule.value * EPS**k).replace(tn.d, 4 - 2 * EPS)
        ours = json_rule_expression(sm, v, couplings, lorentz, UfoConverter(sm, legs, params))
        extra = [EPS.get_name(), EVANESCENT_TAG.get_name(), Expression.PI.get_name()]
        assert same_values(sm, ev, field, ours, target.expand(), legs, rec.components, extra), v["name"]
        assert not same_values(sm, ev, field, ours, (2 * target).expand(), legs, rec.components, extra)


def test_counterterm_model_loads_and_generates(ct_written, hk):
    """Today's HepKIT loads the counterterm model and builds tree diagrams with the two-point counterterms (vertex
    bound 1: unbounded two-point chains crash the generator, scout bug B9); the one-loop gluon self-energy without
    counterterm insertions is unchanged."""
    m = hk.Model.from_json(ct_written.json_string())
    tree = m.process(["G"], ["G"]).generate_diagrams(loops=0, coupling_orders={"CT": 1, "DELTA": 0, "SRC": 0},
                                                       max_vertices=1, maximum_bridges=None, progress=None)
    assert sorted(dg.vertices[0].interaction for dg in tree) == ["DS1_G_G", "dS1_G_G", "dS1d_G_G"]
    loop = generate(m, ["G", "G"], [], 0, 0, set())
    assert len([dg for dg in loop if not any(v.interaction.startswith(("dS1", "DS1", "dS1d")) for v in dg.vertices)]) == 7


# ---------------------------------------------------------------------------------------------- v2 option
def test_v2_schema_option(sm, brst, written, hk):
    v2 = tree_model(sm, WriterOptions(schema="v2", substitutions={"xi": 1}), brst=brst)
    d = v2.data
    assert d["schema_version"] == 2 and "schema_version" not in written.data
    assert not any("v2_overrides" in x for x in d["particles"] + d["lorentz_structures"])
    ghosts = {p["name"]: p["ghost_number"] for p in d["particles"]}
    assert ghosts["gh_G"] == 1 and ghosts["gh_G~"] == -1 and ghosts["src_G"] == -1 and ghosts["src_gh_G"] == -2
    v1_ghosts = {p["name"]: p.get("v2_overrides", {}).get("ghost_number", 0) for p in written.data["particles"]}
    assert v1_ghosts == ghosts
    s1 = {x["name"]: x for x in written.data["lorentz_structures"]}
    for x in d["lorentz_structures"]:
        assert x["structure"] == s1[x["name"]].get("v2_overrides", {}).get("structure", s1[x["name"]]["structure"])
        assert "bmhv::" not in x["structure"]
    assert any("Hat(" in x["structure"] for x in d["lorentz_structures"])
    assert len(hk.Model.from_json(v2.json_string()).vertex_rules) == 174  # today's loader drops the v2 fields


def test_component_samples_and_package_rules(sm, written):
    """Every written vertex is non-zero at its components in D = 6 (the exact zero test of the writer), and a
    component the writer dropped (W1 W1 W1) vanishes."""
    fp = Fingerprinter(sm, dimension=6)
    for name, rec in written.reverse_map.vertices.items():
        rule = written.package_rule(name)
        assert np.any(fp.evaluate(rule.value, component_samples(fp, rec.fields, rec.components, 2)))
    rule = written.package_rule("S0_W1_W2_W3")
    w = [(1, 0)]
    assert not np.any(fp.evaluate(rule.value, component_samples(fp, rule.fields, [w, w, w], 2)))
