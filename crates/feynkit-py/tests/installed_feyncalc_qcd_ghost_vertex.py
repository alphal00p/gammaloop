"""Generated QCD ghost UV poles, including crossed antighost amplitudes.

All graph weights, momentum routing, UV expansion, color/tensor reduction and
IBP are shared APIs.  The only analytic integral input is I(1)|pole = M/eps;
physical loop integration restores i/(16*pi**2), and a4 = gs**2/(16*pi**2).
The SM model supplies the canonical UFO ghost rule; only the gluon gauge
parameter is specialized here.
References: FeynCalcExamples/QCD/OneLoop/{Renormalization,GhGl-Gh}.
"""

import json
from pathlib import Path

from symbolica import E, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
spec = json.loads(model.to_json())
for propagator in spec["propagators"]:
    if propagator["particle"] == "g":
        propagator["numerator"] = (
            "-1𝑖*(UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))"
            "-(1-ghost_uv::xi)*UFO::P(UFO::idx(1,1))*UFO::P(UFO::idx(1,2))"
            "/spenso::dot(UFO::P(spenso::mink(4)),UFO::P(spenso::mink(4))))"
        )
model = hep.Model.from_json(json.dumps(spec))
D, eps, M, mUV, x, integral, xi, s, dA, CA = S(
    "ghost_uv::D",
    "ghost_uv::eps",
    "ghost_uv::M",
    "ghost_uv::mUV",
    "ghost_uv::x",
    "ghost_uv::I",
    "ghost_uv::xi",
    "ghost_uv::s",
    "ghost_uv::dA",
    "ghost_uv::CA",
)
K, P, gs = S("gammalooprs::K", "gammalooprs::P", "UFO::G")
mink, coad, metric, hedge, f, cas = S(
    "spenso::mink",
    "spenso::coad",
    "spenso::g",
    "gammalooprs::hedge",
    "spenso::f",
    "spenso::cas",
)
idx, dim, a, b, c, arguments = S(
    "ghost_uv::idx_",
    "ghost_uv::dim_",
    "ghost_uv::a",
    "ghost_uv::b",
    "ghost_uv::c",
    "ghost_uv::arguments___",
)
den, ed, mom, ms, quad = S(
    "gammalooprs::denom",
    "ghost_uv::ed_",
    "ghost_uv::mom_",
    "ghost_uv::ms_",
    "ghost_uv::quad_",
)
zero, one = E("0"), E("1")
kinematics = hep.Kinematics(D, momenta=[K(0), P(0), P(1)]).with_scalar_product(
    P(0), P(0), s
)
vacuum = hep.Kinematics(D, momenta=[K(0)])
family = hep.IntegralFamily(
    [K(0)], [], [vacuum.scalar_product(K(0), K(0)) - M], kinematics=vacuum
)
reducer = hep.TensorReducer(D).with_integrated_vector(K(0, mink(D)))
color = f(coad(dA, a), coad(dA, b), coad(dA, c))
terms_by_diagram, trees, generated_diagrams = {}, {}, {}
for pdg in (9000005, -9000005):
    for kind, outgoing, loops, count in (
        ("tree", [21, pdg], 0, 1),
        ("vertex", [21, pdg], 1, 2),
        ("self", [pdg], 1, 1),
    ):
        diagrams = (
            model.process([pdg], outgoing, vertex_allow=["V_35", "V_36"])
            .generate_diagrams(
                loops=loops,
                max_vertices=len(outgoing) - 1 + 2 * loops,
                maximum_bridges=0,
                self_energy=None,
                tadpoles=None,
                zero_snails=None,
                numerator_grouping=None,
                progress=None,
            )
            .diagrams
        )
        assert len(diagrams) == count, (pdg, kind, len(diagrams))
        generated_diagrams[pdg, kind] = diagrams
        for diagram in diagrams:
            topology = tuple(sorted(vertex.interaction for vertex in diagram.vertices))
            weight = (
                diagram.overall_factor_expression(evaluate=True)
                * diagram.numerator_prefactor_expression()
            )
            # Compare amputated ghost kernels in canonical field order.
            # Preserve all internal-loop and graph symmetry factors.
            _ordering, _value = S(
                "feynkit_generator_factor::ExternalFermionOrderingSign",
                "ordering_value_",
            )
            _raw = diagram.overall_factor_expression()
            _external_ordering = (_raw / _raw.replace(_ordering(_value), one)).replace(
                _ordering(_value), _value
            )
            weight /= _external_ordering
            assert weight == one
            # Keep every external Lorentz/color slot open; no polarization or
            # contraction with external momentum can hide an unwanted tensor.
            numerator = model.expand_couplings(
                diagram.numerator_expression().to_expression()
            )
            for half in diagram.half_edges:
                if half.edge.data.is_external:
                    numerator = numerator.replace(
                        hedge(half.data, 1), [a, b, c][half.edge.data.external_index]
                    )
            numerator = numerator.replace(coad(8, idx), coad(dA, idx))
            if loops:
                numerator = diagram.uv_expansion(
                    mUV, numerator=numerator
                ).to_expression()
            numerator = diagram.momentum_basis().route_expression(numerator)
            numerator = (
                numerator.replace(mink(dim, idx), mink(D, idx))
                .replace(mink(dim), mink(D))
                .replace(mUV**2, M)
            )
            for match in list(numerator.match(den(ed, mom, ms, quad))):
                values = dict(match)
                formal = family.rewrite_numerator(values[quad], [x])
                assert formal == x
                numerator = numerator.replace(
                    den(values[ed], values[mom], values[ms], values[quad]), formal
                )
            tensor = (
                TensorExpression(numerator.expand())
                .simplify_color()
                .contract()
                .to_expression()
                .to_dots()
                .to_expression()
                .replace(cas(2, coad(dA)), CA)
            )
            scalar = family.rewrite_numerator(
                kinematics.apply(reducer.reduce(tensor)), [x]
            )
            scalar = (weight * scalar).together().expand()
            assert not scalar.matches(K(arguments))
            assert not scalar.matches(hedge(arguments))
            if not loops:
                expected_momentum = P(0, mink(D, b))
                if pdg < 0:
                    expected_momentum -= P(1, mink(D, b))
                assert (scalar + gs * expected_momentum * color).expand() == zero
                trees[pdg] = scalar
                continue
            terms = []
            for monomial, coefficient in scalar.coefficient_list(x):
                power = -int((monomial.derivative(x) * x / monomial).together())
                assert monomial == x**-power
                terms.append(([power], coefficient))
            key = (pdg, kind, topology)
            assert key not in terms_by_diagram
            terms_by_diagram[key] = terms

# Crossing exchanges ghost color slots and sends its incoming momentum to
# minus the outgoing antighost momentum; both sides retain the open gluon slot.
crossed_tree = (
    trees[9000005]
    .replace(P(0, mink(D, b)), -P(0, mink(D, b)) + P(1, mink(D, b)))
    .replace(color, f(coad(dA, c), coad(dA, b), coad(dA, a)))
)
assert (crossed_tree - trees[-9000005]).expand() == zero

targets = sorted({tuple(p) for terms in terms_by_diagram.values() for p, _ in terms})
solution = hep.IBPFamily(family, name="qcd_ghost_vertex").reduce_laporta(
    [list(target) for target in targets], max_depth=2
)
assert solution.residuals == [[1]]
vertex_poles = {9000005: zero, -9000005: zero}
self_poles = {}
for (pdg, kind, topology), terms in terms_by_diagram.items():
    reduced = sum(
        (
            coefficient * solution.reduce(powers, integral=integral)
            for powers, coefficient in terms
        ),
        zero,
    ).together()
    pole = (
        (reduced / gs**2)
        .replace(integral(1), M / eps)
        .replace(D, 4 - 2 * eps)
        .series(eps, 0, -1)
        .to_expression()
        .expand()
    )
    assert pole.derivative(M).expand() == zero
    assert pole.derivative(mUV).expand() == zero
    assert pole.coefficient(eps**-2) == zero
    assert not pole.matches(integral(arguments))
    if kind == "self":
        self_poles[pdg] = pole
        expected = CA * (xi - 3) * s * metric(coad(dA, a), coad(dA, b)) / (4 * eps)
        assert (pole - expected).expand() == zero
        assert pole.replace(s, 0) == zero  # No auxiliary ghost mass pole.
    else:
        assert topology in (("V_35", "V_35", "V_35"), ("V_35", "V_35", "V_36"))
        multiplicity = 3 if "V_36" in topology else 1
        expected_ratio = multiplicity * CA * xi / (8 * eps)
        tree = trees[pdg].replace(D, 4)
        physical_pole = (Symbol.I * pole).expand()
        assert (physical_pole - tree * expected_ratio).expand() == zero
        ratio = (physical_pole / tree).together()
        assert ratio == expected_ratio
        vertex_poles[pdg] += physical_pole
for pdg, pole in vertex_poles.items():
    assert (pole - trees[pdg].replace(D, 4) * CA * xi / (2 * eps)).expand() == zero

# The gluon field counterterm is an explicit reference input from the separate
# generated QCD renormalization calculation. The ghost kinetic pole is computed
# here; combining it with the vertex supplies an independent coupling check.
Nf = S("ghost_uv::Nf")
vertex_ratios = {
    pdg: (pole / trees[pdg].replace(D, 4)).together()
    for pdg, pole in vertex_poles.items()
}
ghost_ct = {
    pdg: -pole.coefficient(metric(coad(dA, a), coad(dA, b))).coefficient(s)
    for pdg, pole in self_poles.items()
}
gluon_ct = ((13 - 3 * xi) * CA - 4 * Nf) / (6 * eps)
coupling_ct = {
    pdg: (-vertex_ratios[pdg] - ghost_ct[pdg] - gluon_ct / 2).together()
    for pdg in trees
}
for counterterm in coupling_ct.values():
    assert (counterterm + (11 * CA - 2 * Nf) / (6 * eps)).together() == zero
    assert counterterm.derivative(xi).expand() == zero
print(
    "Generated ghost and antighost UV: full tensors, crossing, two IBP targets and coupling consistency passed",
    solution.stats,
)
