"""Generated phi4 MS/MSbar counterterms and finite one-loop scattering.

References: FeynCalc Phi4/OneLoop/Renormalization and PhiPhi-PhiPhi.
Native graph routing determines all three bubble invariants. Independent
Feynman-parameter integrals check massive and massless finite coefficients.
"""

import cmath
import json
import math

import numpy as np
from symbolica import E, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model.phi4()
particle = model.particle("phi")
mass, coupling = S("UFO::mass", "UFO::lam")
d, M, s, t, u, mu2, p2, eps = S(
    "phi4_one::d",
    "phi4_one::M",
    "phi4_one::s",
    "phi4_one::t",
    "phi4_one::u",
    "phi4_one::mu2",
    "phi4_one::p2",
    "phi4_one::eps",
)
Q, K, P = S("gammalooprs::Q", "gammalooprs::K", "gammalooprs::P")
zero, one, pi = E("0"), E("1"), Symbol.PI
kinematics = hep.Kinematics.mandelstam([P(i) for i in range(4)], [M] * 4, [s, t, u])
routing_kinematics = hep.Kinematics(momenta=[K(0), *[P(i) for i in range(4)]])
options = {
    "maximum_bridges": 0,
    "self_energy": None,
    "tadpoles": None,
    "zero_snails": None,
    "numerator_grouping": None,
    "progress": None,
}
diagrams = {}
for label, legs, loops, vertices in [
    ("self_energy", 1, 1, 1),
    ("vertex", 2, 1, 2),
    ("tree", 2, 0, 1),
]:
    diagrams[label] = (
        model.process([particle] * legs, [particle] * legs)
        .generate_diagrams(loops=loops, max_vertices=vertices, **options)
        .diagrams
    )
assert [len(diagrams[k]) for k in ("self_energy", "vertex", "tree")] == [1, 3, 1]
tree = diagrams["tree"][0]
tree_amplitude = (
    model.expand_couplings(tree.numerator_expression().to_expression())
    * tree.overall_factor_expression(evaluate=True)
    * tree.numerator_prefactor_expression()
)
assert tree_amplitude == -Symbol.I * coupling

# Normalize each edge momentum to K(0)+shift. The graph's propagators, not
# labels assigned to its drawings, determine s, t and u.
channel_coefficients, channel_reductions = {}, {}
for diagram in diagrams["vertex"]:
    basis = diagram.momentum_basis()
    family = diagram.propagator_family(kinematics=kinematics)
    shifts = []
    for edge, denominator in zip(
        diagram.internal_edges, family.denominators, strict=True
    ):
        sign = basis.edge_signatures[edge.id].loops[0]
        assert abs(sign) == 1
        shift = (basis.route_expression(Q(edge.id)) / sign - K(0)).expand()
        reconstructed = (
            kinematics.apply(
                routing_kinematics.scalar_product(K(0) + shift, K(0) + shift)
            )
            - M
        )
        assert (denominator.replace(mass**2, M) - reconstructed).expand() == zero
        shifts.append(shift)
    invariant = kinematics.scalar_product(shifts[0] - shifts[1], shifts[0] - shifts[1])
    assert invariant in (s, t, u) and invariant not in channel_coefficients
    coefficient = (
        model.expand_couplings(diagram.numerator_expression().to_expression())
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    assert coefficient == coupling**2 / 2
    bubble_kinematics = hep.Kinematics(d, momenta=[K(0), P(0)]).with_scalar_product(
        P(0), P(0), invariant
    )
    bubble_family = hep.IntegralFamily(
        [K(0)],
        [P(0)],
        [bubble_kinematics.scalar_product(q, q) - M for q in (K(0), K(0) - P(0))],
        kinematics=bubble_kinematics,
    )
    reduction = oneloop.reduce(bubble_family, [1, 1])
    assert len(reduction.terms) == 1
    channel_coefficients[invariant] = coefficient
    channel_reductions[invariant] = reduction
assert set(channel_coefficients) == {s, t, u}

# OneLOop supplies Laurent coefficients in its common MSbar master convention.
# Extract the UV branch for nonzero mass; the finite calculation below also
# checks massless nonexceptional channels and the massive threshold.
loop_poles, finite_parts = {}, {}
vertex_finite, vertex_pole = zero, zero
for invariant, reduction in channel_reductions.items():
    coefficients = oneloop.reduction_coefficients(reduction, mu2)
    master = reduction.terms[0][1].to_expression(mu2)
    pole_expression = oneloop.get_expression(master, coefficient=-1)
    pole = oneloop.select_branch(
        pole_expression, [Replacement(invariant, one), Replacement(M, one)]
    )
    assert pole == one
    weight = channel_coefficients[invariant] / coupling**2
    vertex_pole += weight * pole
    vertex_finite += weight * coefficients[0]
self_diagram = diagrams["self_energy"][0]
self_numerator = (
    model.expand_couplings(self_diagram.numerator_expression().to_expression())
    * self_diagram.overall_factor_expression(evaluate=True)
    * self_diagram.numerator_prefactor_expression()
)
assert self_numerator == coupling / 2
self_family = self_diagram.propagator_family()
assert len(self_family.denominators) == 1
tadpole_kinematics = hep.Kinematics(d, momenta=[K(0)])
tadpole_family = hep.IntegralFamily(
    [K(0)],
    [],
    [tadpole_kinematics.scalar_product(K(0), K(0)) - M],
    kinematics=tadpole_kinematics,
)
self_reduction = oneloop.reduce(tadpole_family, [1])
self_master = self_reduction.terms[0][1].to_expression(mu2)
self_pole = (
    self_numerator / coupling * oneloop.get_expression(self_master, coefficient=-1)
)
self_finite = (
    self_numerator / coupling * oneloop.reduction_coefficients(self_reduction, mu2)[0]
)
assert (self_pole - M / 2).together() == zero
assert (vertex_pole - E("3/2")).together() == zero
loop_poles = {"self_energy": self_pole, "vertex": vertex_pole}
finite_parts = {"self_energy": self_finite, "vertex": vertex_finite}

# h is g/(16*pi^2). Counterterm rules are coefficients of the bare local
# operators, and actual generated numerators supply the matching matrix.
h, delta_field, delta_mass, delta_vertex = S(
    "phi4_one::h",
    "phi4_one::delta_field",
    "phi4_one::delta_mass",
    "phi4_one::delta_vertex",
)
Zfield, Zmass, Zvertex = (1 + h * x for x in (delta_field, delta_mass, delta_vertex))
specification = json.loads(model.to_json())
specification["orders"].append({"name": "CT", "expansion_order": 1, "hierarchy": 1})
for label, valence, lorentz, factor in [
    ("kinetic", 2, "P(dummy(1),1)*P(dummy(1),1)", Symbol.I * (Zfield - 1)),
    ("mass", 2, "1", -Symbol.I * mass**2 * (Zfield * Zmass - 1)),
    ("quartic", 4, "1", -Symbol.I * coupling * (Zvertex * Zfield**2 - 1)),
]:
    name = "CT_" + label
    specification["lorentz_structures"].append(
        {"name": name, "spins": [1] * valence, "structure": lorentz}
    )
    specification["couplings"].append(
        {
            "name": name,
            "expression": repr(factor.expand().coefficient(h) * h),
            "orders": [["SCALAR", int(valence == 4)], ["CT", 1]],
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
ct_diagrams, ct_amplitudes = {}, {}
self_kinematics = hep.Kinematics().with_scalar_product(P(0), P(0), p2)
for label, legs, count in [("self_energy", 1, 2), ("vertex", 2, 1)]:
    result = ct_model.process(
        [particle.name] * legs, [particle.name] * legs
    ).generate_diagrams(loops=0, max_vertices=1, coupling_orders={"CT": 1}, **options)
    assert len(result.diagrams) == count
    ct_diagrams[label] = result.diagrams
    amplitude = zero
    for diagram in result.diagrams:
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
    ct_amplitudes[label] = (
        (
            amplitude.replace(mass**2, M)
            / (Symbol.I * h * (coupling if legs == 2 else one))
        )
        .together()
        .expand()
    )
rows = [
    ct_amplitudes["self_energy"].coefficient(p2),
    (ct_amplitudes["self_energy"].replace(p2, zero) / M).together().expand(),
    ct_amplitudes["vertex"],
]
unknowns = [delta_field, delta_mass, delta_vertex]
ct_matrix = Matrix.from_linear(
    3, 3, [row.coefficient(x) for row in rows for x in unknowns]
)
rhs = Matrix.vec([zero, -loop_poles["self_energy"] / M, -loop_poles["vertex"]])
solved = ct_matrix.solve(rhs)
residues = [solved[i, 0].to_expression() for i in range(3)]
assert residues == [zero, E("1/2"), E("3/2")]
L, gamma_e = S("phi4_one::log4pi", "phi4_one::gamma_E")
renormalization_constants = {}
for scheme, subtraction in [("MS", 1 / eps), ("MSbar", 1 / eps + L - gamma_e)]:
    constants = [x * subtraction for x in residues]
    renormalization_constants[scheme] = constants
    replacement = [Replacement(x, y) for x, y in zip(unknowns, constants, strict=True)]
    for label in loop_poles:
        total = (
            loop_poles[label] * (1 / eps + L - gamma_e)
            + ct_amplitudes[label].replace_multiple(replacement)
        ).together()
        expected = zero if scheme == "MSbar" else loop_poles[label] * (L - gamma_e)
        assert (total - expected).together() == zero

# Independently check OneLOop's finite branch against the parameter integral
# -integral_0^1 log((M-q*x*(1-x)-i0)/mu2) dx. Above threshold integrate the two
# real logarithmic roots analytically, including the cut's imaginary length.
finite_evaluator = oneloop.compile_native([vertex_finite], [M, s, t, u, mu2])
points = [
    (0.0, 5.0, -1.0, 2.0),
    (0.0, 20.0, -6.0, 1.0),
    (1.0, 4.0, 0.0, 1.0),
    (1.0, 5.0, -0.25, 1.0),
    (1.0, 12.0, -3.0, 2.0),
    (4.0, 25.0, -4.0, 3.0),
]
nodes, weights = np.polynomial.legendre.leggauss(256)
numeric_checks = []
for mass_squared, sv, tv, scale in points:
    uv = 4 * mass_squared - sv - tv
    reference = 0j
    for qv in (sv, tv, uv):
        if mass_squared == 0:
            expected = 2 - cmath.log(complex(-qv / scale, -0.0))
        elif qv >= 4 * mass_squared:
            beta = math.sqrt(1 - 4 * mass_squared / qv)
            roots = [(1 - beta) / 2, (1 + beta) / 2]
            expected = (
                -math.log(qv / scale)
                - sum((1 - r) * math.log(1 - r) + r * math.log(r) - 1 for r in roots)
                + 1j * math.pi * beta
            )
        else:
            x = (nodes + 1) / 2
            expected = (
                -np.dot(weights, np.log((mass_squared - qv * x * (1 - x)) / scale)) / 2
            )
        actual = complex(oneloop.b0(qv, mass_squared, mass_squared, scale)[0])
        assert abs(actual - expected) < 2e-10, (qv, mass_squared, actual, expected)
        reference += expected / 2
    actual = complex(
        finite_evaluator.evaluate_complex(
            [[complex(v) for v in (mass_squared, sv, tv, uv, scale)]]
        )[0, 0]
    )
    assert abs(actual - reference) < 4e-10, (actual, reference)
    numeric_checks.append((mass_squared, sv, tv, uv, scale, actual, reference))

# In the physical massless region s>0,t<0,u<0, the Feynman prescription fixes
# the imaginary part. Differentiate before taking s->infinity at fixed t.
log, absolute = S("log", "abs")
physical_finite = zero
for invariant, sign in [(s, 1), (t, -1), (u, -1)]:
    massless_master = S("oneloopmaster::B0")(invariant, zero, zero, mu2)
    branch = oneloop.select_branch(
        oneloop.get_expression(massless_master, coefficient=0),
        [Replacement(invariant, E(str(sign))), Replacement(mu2, one)],
    )
    branch = branch.replace(
        absolute(-invariant / mu2),
        (sign * invariant / mu2).replace(u, -s - t).expand(),
    )
    physical_finite += channel_coefficients[invariant] / coupling**2 * branch
physical_finite = physical_finite.replace(u, -s - t)
expected_massless = (
    6 - log(s / mu2) - log(-t / mu2) - log(((s + t) / mu2).expand()) + Symbol.I * pi
) / 2
assert (physical_finite - expected_massless).expand() == zero
asymptotic_parameter = S("phi4_one::inverse_s")
logarithmic_coefficient = (
    (s * physical_finite.derivative(s))
    .replace(s, 1 / asymptotic_parameter)
    .series(asymptotic_parameter, 0, 0)
    .to_expression()
)
assert logarithmic_coefficient == -one
leading_amplitude = (
    Symbol.I * coupling**2 / (16 * pi**2) * logarithmic_coefficient * log(s / mu2)
)
assert (
    leading_amplitude + Symbol.I * coupling**2 * log(s / mu2) / (16 * pi**2)
).together() == zero
print(
    "Generated phi4 one-loop MS/MSbar constants, six finite scattering points and the high-energy logarithm passed"
)
