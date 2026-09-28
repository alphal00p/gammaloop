"""Generated phi3 renormalization and finite zero-momentum triangle.

Reference: FeynCalc Phi3/OneLoop/Renormalization, in four dimensions.
Native graph factors and propagators feed the shared tadpole IBP recurrence.
OneLOop evaluates its master; direct parameter integrals check finite terms.
"""

import json
import math

import numpy as np
from symbolica import E, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model.phi3()
particle = model.particle("phi")
mass, coupling = S("UFO::mass", "UFO::g")
d, M, p2, mu2, eps, k, coordinate, integral = S(
    "phi3::d",
    "phi3::M",
    "phi3::p2",
    "phi3::mu2",
    "phi3::eps",
    "phi3::k",
    "phi3::x",
    "phi3::I",
)
Q, K, P = S("gammalooprs::Q", "gammalooprs::K", "gammalooprs::P")
index, external = S("phi3::index_", "phi3::external_")
zero, one, pi = E("0"), E("1"), Symbol.PI
kinematics = hep.Kinematics(d, momenta=[k])
family = hep.IntegralFamily(
    [k], [], [kinematics.scalar_product(k, k) - M], kinematics=kinematics
)
ibp = hep.IBPFamily(family, name="phi3_tadpole")
recurrence = ibp.solve_parametric([True], max_depth=1)
assert recurrence.rules
options = {
    "maximum_bridges": 0,
    "self_energy": None,
    "tadpoles": None,
    "zero_snails": None,
    "numerator_grouping": None,
    "progress": None,
}
diagrams, inputs, coefficients = {}, {}, {}
for label, incoming, outgoing, loops in [
    ("self_energy", 1, 1, 1),
    ("vertex", 2, 1, 1),
    ("tree", 2, 1, 0),
]:
    generated = model.process(
        [particle] * incoming, [particle] * outgoing
    ).generate_diagrams(
        loops=loops, max_vertices=incoming + outgoing if loops else 1, **options
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    diagrams[label] = diagram
    coefficient = (
        model.expand_couplings(diagram.numerator_expression().to_expression())
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    coefficients[label] = coefficient
    if not loops:
        continue
    # Zero external momenta in the actual routed denominators. A nonzero mass
    # makes this Taylor coefficient infrared safe, including the finite triangle.
    generated_family = diagram.propagator_family(kinematics=hep.Kinematics(d))
    denominators = [
        den.replace(P(external, index), zero)
        .replace(K(0, index), k(index))
        .replace(mass**2, M)
        for den in generated_family.denominators
    ]
    powers = len(denominators)
    assert powers == incoming + outgoing
    assert all(
        family.rewrite_numerator(den, [coordinate]) == coordinate
        for den in denominators
    )
    inputs[label] = coefficient * integral(powers)
assert coefficients == {
    "self_energy": coupling**2 / 2,
    "vertex": coupling**3,
    "tree": -Symbol.I * coupling,
}

# Apply exactly one recurrence at each positive power, then back-substitute.
reductions = {1: integral(1)}
for power in (2, 3):
    terms = recurrence.reduce([power], integral=integral)
    reductions[power] = terms.replace(
        integral(power - 1), reductions[power - 1]
    ).together()
assert (reductions[2] / integral(1) - (d - 2) / (2 * M)).together() == zero
assert (reductions[3] / integral(1) - (d - 2) * (d - 4) / (8 * M**2)).together() == zero
laporta = ibp.reduce_laporta([[2], [3]], max_depth=2)
for power in (2, 3):
    assert (
        laporta.reduce([power], integral=integral) - reductions[power]
    ).together() == zero

master_reduction = oneloop.reduce(family, [1])
master = master_reduction.terms[0][1].to_expression(mu2)
master_pole = oneloop.get_expression(master, coefficient=-1)
assert master_pole == M
master_finite = oneloop.reduction_coefficients(master_reduction, mu2)[0]
finite_symbol = S("phi3::Afinite")
reduced, poles, finite = {}, {}, {}
for label, power in [("self_energy", 2), ("vertex", 3)]:
    reduced[label] = inputs[label].replace(integral(power), reductions[power])
    expanded = (
        reduced[label]
        .replace(d, 4 - 2 * eps)
        .replace(integral(1), M / eps + finite_symbol)
        .series(eps, 0, 0)
        .to_expression()
        .expand()
    )
    poles[label] = expanded.coefficient(eps**-1)
    finite[label] = (
        (expanded - poles[label] / eps).together().replace(finite_symbol, master_finite)
    )
assert poles == {"self_energy": coupling**2 / 2, "vertex": zero}
assert (finite["vertex"] + coupling**3 / (2 * M)).together() == zero

# Keep generic p^2 for the self-energy. Its actual graph routing gives B0,
# whose UV branch is independent of p^2; no field pole is hidden by p=0.
self_diagram = diagrams["self_energy"]
basis = self_diagram.momentum_basis()
routed_kinematics = hep.Kinematics(momenta=[K(0), P(0), P(1)]).with_scalar_product(
    P(0), P(0), p2
)
shifts = []
for edge in self_diagram.internal_edges:
    sign = basis.edge_signatures[edge.id].loops[0]
    assert abs(sign) == 1
    shift = (basis.route_expression(Q(edge.id)) / sign - K(0)).expand()
    shifts.append(shift)
invariant = routed_kinematics.scalar_product(
    shifts[0] - shifts[1], shifts[0] - shifts[1]
)
assert invariant == p2
bubble_kinematics = hep.Kinematics(d, momenta=[K(0), P(0)]).with_scalar_product(
    P(0), P(0), invariant
)
bubble_family = hep.IntegralFamily(
    [K(0)],
    [P(0)],
    [bubble_kinematics.scalar_product(q, q) - M for q in (K(0), K(0) - P(0))],
    kinematics=bubble_kinematics,
)
bubble_reduction = oneloop.reduce(bubble_family, [1, 1])
bubble_master = bubble_reduction.terms[0][1].to_expression(mu2)
bubble_pole = oneloop.select_branch(
    oneloop.get_expression(bubble_master, coefficient=-1),
    [Replacement(M, one), Replacement(p2, one)],
)
assert bubble_pole * coefficients["self_energy"] == poles["self_energy"]
assert poles["self_energy"].derivative(p2) == zero
bubble_finite = (
    coefficients["self_energy"]
    * oneloop.reduction_coefficients(bubble_reduction, mu2)[0]
)

# Expand the bare operators; the generated local diagrams determine the
# matching matrix. h counts loops, with 1/(16*pi^2) removed from each delta.
h, field, mass_ct, vertex = S("phi3::h", "phi3::field", "phi3::mass_ct", "phi3::vertex")
Zfield, Zmass, Zvertex = (1 + h * x for x in (field, mass_ct, vertex))
specification = json.loads(model.to_json())
specification["orders"].append({"name": "CT", "expansion_order": 1, "hierarchy": 1})
for label, valence, lorentz, factor in [
    ("kinetic", 2, "P(dummy(1),1)*P(dummy(1),1)", Symbol.I * (Zfield - 1)),
    ("mass", 2, "1", -Symbol.I * mass**2 * (Zfield * Zmass - 1)),
    ("cubic", 3, "1", -Symbol.I * coupling * (Zvertex * Zfield ** E("3/2") - 1)),
]:
    name = "CT_" + label
    coefficient = factor.series(h, 0, 1).to_expression().expand().coefficient(h)
    specification["lorentz_structures"].append(
        {"name": name, "spins": [1] * valence, "structure": lorentz}
    )
    specification["couplings"].append(
        {
            "name": name,
            "expression": repr(h * coefficient),
            "orders": [["SCALAR", int(valence == 3)], ["CT", 1]],
            "value": None,
        }
    )
    specification["vertex_rules"].append(
        {
            "name": name,
            "particles": [particle.name] * valence,
            "color_structures": ["1"],
            "lorentz_structures": [name],
            "couplings": [[name]],
        }
    )
ct_model = hep.Model.from_json(json.dumps(specification))
ct_diagrams, counterterms = {}, {}
self_kinematics = hep.Kinematics().with_scalar_product(P(0), P(0), p2)
for label, incoming, count in [("self_energy", 1, 2), ("vertex", 2, 1)]:
    generated = ct_model.process(
        [particle.name] * incoming, [particle.name]
    ).generate_diagrams(loops=0, max_vertices=1, coupling_orders={"CT": 1}, **options)
    ct_diagrams[label] = generated.diagrams
    assert len(generated.diagrams) == count
    amplitude = zero
    for diagram in generated.diagrams:
        numerator = ct_model.expand_couplings(
            diagram.numerator_expression().to_expression()
        )
        numerator = (
            TensorExpression(numerator)
            .contract()
            .to_expression()
            .to_dots()
            .to_expression()
        )
        numerator = diagram.momentum_basis().route_expression(numerator)
        amplitude += (
            self_kinematics.apply(numerator)
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
        )
    counterterms[label] = (
        (amplitude.replace(mass**2, M) / (Symbol.I * h)).together().expand()
    )
unknowns = [field, mass_ct, vertex]
rows = [
    counterterms["self_energy"].coefficient(p2),
    counterterms["self_energy"].replace(p2, zero),
    counterterms["vertex"],
]
ct_matrix = Matrix.from_linear(
    3, 3, [row.coefficient(x) for row in rows for x in unknowns]
)
solution = ct_matrix.solve(Matrix.vec([zero, -poles["self_energy"], -poles["vertex"]]))
residues = [solution[row, 0].to_expression() for row in range(3)]
assert residues == [zero, coupling**2 / (2 * M), zero]
log4pi, gamma_e = S("phi3::log4pi", "phi3::gamma_E")
constants = {}
for scheme, delta in [("MS", 1 / eps), ("MSbar", 1 / eps + log4pi - gamma_e)]:
    constants[scheme] = [1 + residue * delta / (16 * pi**2) for residue in residues]
    replacements = [
        Replacement(x, residue * delta)
        for x, residue in zip(unknowns, residues, strict=True)
    ]
    for label, pole in poles.items():
        remainder = (
            pole * (1 / eps + log4pi - gamma_e)
            + counterterms[label].replace_multiple(replacements)
        ).together()
        expected = zero if scheme == "MSbar" else pole * (log4pi - gamma_e)
        assert (remainder - expected).together() == zero

# Independent finite checks, including the zero-momentum bubble where the
# tadpole reduction's O(eps) coefficient changes the finite answer.
evaluator = oneloop.compile_native(
    [finite["self_energy"], finite["vertex"], bubble_finite], [M, coupling, mu2, p2]
)
nodes, weights = np.polynomial.legendre.leggauss(128)
nodes, weights = (nodes + 1) / 2, weights / 2
points = [
    (1.0, 0.5, 1.0, 0.0),
    (2.0, 1.0, 3.0, -4.0),
    (4.0, 2.0, 0.5, 3.0),
    (0.25, 0.1, 2.0, 0.5),
]
numeric_checks = []
for mv, gv, scale, pv in points:
    values = [
        complex(x)
        for x in evaluator.evaluate_complex(
            [[complex(x) for x in (mv, gv, scale, pv)]]
        )[0]
    ]
    reference = [
        -(gv**2) * math.log(mv / scale) / 2,
        -(gv**3) / (2 * mv),
        -(gv**2)
        * sum(
            float(w) * math.log((mv - pv * float(x) * (1 - float(x))) / scale)
            for x, w in zip(nodes, weights, strict=True)
        )
        / 2,
    ]
    assert all(
        abs(actual - expected) < 2e-11
        for actual, expected in zip(values, reference, strict=True)
    )
    numeric_checks.append((mv, gv, scale, pv, values, reference))
print(
    "Generated phi3 diagrams, native IBP, MS/MSbar constants and finite triangle passed",
    laporta.stats,
)
