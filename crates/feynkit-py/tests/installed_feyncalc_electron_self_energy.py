"""Generated full massive off-shell electron self-energy, including finite terms.

Primary benchmark: https://arxiv.org/abs/hep-ph/0008171, Eqs. (2.19)--(2.21),
with its
xi_paper=1-xi and CF=1; its Appendix C.6 supplies the independent B0 formula.
Keep the native ordered external-state weight, then state the conversion to
an amputated Dirac kernel explicitly. No UV expansion or mass rearrangement.
"""

import json
from math import log, pi
from pathlib import Path

from symbolica import E, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
D, eps, xi, s = S("qed_full::D", "qed_full::eps", "qed_full::xi", "qed_full::s")
specification = json.loads(model.to_json())
for propagator in specification["propagators"]:
    if propagator["particle"] == "a":
        propagator["numerator"] = (
            "-1𝑖*(UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))"
            "-(1-qed_full::xi)*UFO::P(UFO::idx(1,1))*UFO::P(UFO::idx(1,2))"
            "/spenso::dot(UFO::P(spenso::mink(4)),UFO::P(spenso::mink(4))))"
        )
model = hep.Model.from_json(json.dumps(specification))
electron, photon = (model.particle_by_pdg(pdg) for pdg in (11, 22))
assert electron.mass_parameter == "Me"
vertices = [
    v
    for v in model.vertex_rules
    if sorted(v.particles) == sorted([electron.name, electron.antiname, photon.name])
]
assert len(vertices) == 1
result = model.generate_diagrams(
    [electron],
    [electron],
    loops=1,
    max_vertices=2,
    maximum_bridges=0,
    vertex_allow=vertices,
    self_energy=None,
    tadpoles=None,
    zero_snails=None,
    numerator_grouping=None,
    progress=None,
)
assert len(result.diagrams) == 1
diagram = result.diagrams[0]
K, P, mass, charge = S("gammalooprs::K", "gammalooprs::P", "UFO::Me", "UFO::ee")
mink, bis, gamma, metric = S(
    "spenso::mink", "spenso::bis", "spenso::gamma", "spenso::g"
)
index, wave, mu, args = S(
    "qed_full::index_", "qed_full::wave_", "qed_full::mu", "qed_full::args___"
)
ordering, value = S(
    "feynkit_generator_factor::ExternalFermionOrderingSign", "qed_full::value_"
)
zero, one = E("0"), E("1")
kinematics = hep.Kinematics(D, momenta=[K(0), P(0)]).with_scalar_product(P(0), P(0), s)
family = diagram.integral_family(kinematics=kinematics)
coordinates = S("qed_full::d0", "qed_full::d1")
q2 = kinematics.scalar_product(K(0), K(0))
pq = kinematics.scalar_product(P(0), K(0))
assert (family.denominators[0] - q2 + mass**2).expand() == zero
assert (family.denominators[1] - s - q2 + 2 * pq).expand() == zero
assert [e.particle_name for e in diagram.internal_edges] == [electron.name, photon.name]

# The UFO vertex is -i e gamma; electron and photon propagator numerators
# carry +i and -i. Their product gives the unweighted kernel -e^2 N.
# Both generated e->e and e->gamma e amplitudes have an additional external
# Wick-ordering factor -1. Establish that sign independently at tree level.
tree_result = model.generate_diagrams(
    [electron],
    [photon, electron],
    max_vertices=1,
    vertex_allow=vertices,
    numerator_grouping=None,
    progress=None,
)
assert len(tree_result.diagrams) == 1
tree_diagram = tree_result.diagrams[0]
tree_ports = {}
for edge in tree_diagram.external_edges:
    representation = mink if edge.external_index == 1 else bis
    tree_ports[edge.external_index] = dict(
        next(
            tree_diagram.projector_expression().match(
                wave(edge.id, representation(4, index)), max_level=0
            )
        )
    )[index]
tree_num = model.expand_couplings(
    tree_diagram.numerator_expression().to_expression()
).replace(mink(4, index), mink(D, index))
tree_probe = gamma(
    bis(4, tree_ports[0]), bis(4, tree_ports[2]), mink(D, tree_ports[1])
) / (4 * D)
tree_coupling = (
    TensorExpression((tree_num * tree_probe).expand())
    .simplify_gamma()
    .expand()
    .to_expression()
)
assert (tree_coupling + Symbol.I * charge).expand() == zero
raw_factor = diagram.overall_factor_expression()
external_ordering = (raw_factor / raw_factor.replace(ordering(value), one)).replace(
    ordering(value), value
)
assert external_ordering == -one
assert diagram.overall_factor_expression(evaluate=True) == external_ordering
assert tree_diagram.overall_factor_expression(evaluate=True) == external_ordering
assert (
    diagram.numerator_prefactor_expression()
    == tree_diagram.numerator_prefactor_expression()
    == one
)
assert (tree_coupling * external_ordering - Symbol.I * charge).expand() == zero

ports = {
    edge.external_index: dict(
        next(
            diagram.projector_expression().match(
                wave(edge.id, bis(4, index)), max_level=0
            )
        )
    )[index]
    for edge in diagram.external_edges
}
numerator = (
    model.expand_couplings(diagram.numerator_expression(in_lmb=True).to_expression())
    .replace(mink(4, index), mink(D, index))
    .replace(mink(4), mink(D))
)
probes = [
    gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
    * P(0, mink(D, mu))
    / (4 * s),
    metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass),
]
native_factor = (
    diagram.overall_factor_expression(evaluate=True)
    * diagram.numerator_prefactor_expression()
)
# Independent D-dimensional gamma contraction, k=p-q and q=K. The scalar
# coefficients multiply pslash and m; no external on-shell equation is used.
photon_square = s + q2 - 2 * pq
independent_numerators = [
    (3 - D - xi)
    + (D - 1 - xi) * (s - pq) / s
    - 2 * (1 - xi) * (s - pq) ** 2 / (s * photon_square),
    D - 1 + xi,
]
terms_by_structure, targets, traced_coefficients = [], set(), []
for projector, independent in zip(probes, independent_numerators, strict=True):
    trace = (
        TensorExpression((numerator * projector).expand())
        .simplify_gamma()
        .expand()
        .to_dots()
        .to_expression()
    )
    scalar = (kinematics.apply(trace) * native_factor / charge**2).together()
    assert (scalar - independent).together() == zero
    traced_coefficients.append(scalar)
    rational = family.rewrite_numerator(scalar, coordinates).together().expand()
    terms = []
    # The graph supplies one power of each denominator. Its longitudinal
    # numerator has another inverse photon square; retain that raised power.
    for monomial, coefficient in rational.coefficient_list(*coordinates):
        powers = tuple(
            1 - int((monomial.derivative(c) * c / monomial).together())
            for c in coordinates
        )
        assert (
            monomial
            - coordinates[0] ** (1 - powers[0]) * coordinates[1] ** (1 - powers[1])
        ).together() == zero
        assert not coefficient.matches(K(args))
        assert not coefficient.matches(P(args))
        assert all(coefficient.derivative(c) == zero for c in coordinates)
        targets.add(powers)
        terms.append((powers, coefficient))
    terms_by_structure.append(terms)
assert (1, 2) in targets
solution = hep.IBPFamily(family, name="qed_full").reduce_laporta(
    [list(p) for p in sorted(targets)], max_depth=2
)
assert {tuple(p) for p in solution.residuals} == {(1, 0), (1, 1)}
I, A, B = S("qed_full::I", "qed_full::A", "qed_full::B")
x0, x1 = S("qed_full::x0", "qed_full::x1")
U, F = family.symanzik([x0, x1])
assert (U - x0 - x1).expand() == zero
assert (F - mass**2 * x0 * U + s * x0 * x1).expand() == zero
# Hence I(1,0)=A0(m²), I(1,1)=B0(s; m²,0). Massless pinches are scaleless
# and are removed by the shared IBP solver, not assigned a numerical master.
master_rules = [Replacement(I(1, 0), A), Replacement(I(1, 1), B)]
native_coefficients = [
    sum(
        (
            coefficient * solution.reduce(list(powers), integral=I)
            for powers, coefficient in terms
        ),
        zero,
    )
    .replace_multiple(master_rules)
    .together()
    for terms in terms_by_structure
]
# Native ordered matrix element = i*a4*[V_native pslash + S_native m],
# a4=e²/(16π²). The amputated inverse-propagator insertion Γ₂=-iΣ is the
# ordered result divided by its explicitly established external Wick sign.
# Thus conventional Σ/a4 has the native coefficients, while Γ₂/(i*a4)
# has their opposites. No internal graph factor is changed or discarded.
amputated_coefficients = [
    (coefficient / external_ordering).together() for coefficient in native_coefficients
]
exact_reference = [xi * (D - 2) * ((s + mass**2) * B - A) / (2 * s), -(D - 1 + xi) * B]
for native, amputated, reference in zip(
    native_coefficients, amputated_coefficients, exact_reference, strict=True
):
    assert (amputated - reference).together() == zero
    assert (native + reference).together() == zero
    for master in (A, B):
        coefficient = amputated.expand().coefficient(master).together()
        assert coefficient.series(D, 4, 0).to_expression().replace(
            D, 4
        ) == coefficient.replace(D, 4)
        assert not coefficient.matches(I(args))
        assert not coefficient.matches(K(args))
        assert all(coefficient.derivative(c) == zero for c in coordinates)
assert amputated_coefficients[0].replace(xi, 0).together() == zero

Af, Bf, scale2, L = S("qed_full::Af", "qed_full::Bf", "qed_full::scale2", "qed_full::L")
laurent = [
    coefficient.replace(A, mass**2 / eps + Af)
    .replace(B, 1 / eps + Bf)
    .replace(D, 4 - 2 * eps)
    .series(eps, 0, 0)
    .to_expression()
    .expand()
    for coefficient in amputated_coefficients
]
expected_poles = [xi, -(3 + xi)]
expected_finite = [xi * ((s + mass**2) * Bf - Af) / s - xi, -(3 + xi) * Bf + 2]
finite = []
for expression, pole, reference in zip(
    laurent, expected_poles, expected_finite, strict=True
):
    assert expression.coefficient(eps**-2) == zero
    assert (expression.coefficient(eps**-1) - pole).together() == zero
    constant = dict(expression.coefficient_list(eps))[one]
    assert (constant - reference).together() == zero
    finite.append(constant)
# Exact s->0 limit: Bf=1-L+s/(2m²)+... . No evaluation of the projector's
# removable 1/s singularity at s=0. The on-shell value is finite but its
# derivative is infrared singular, so this is not an on-shell Z2 calculation.
zero_invariant = [
    expression.replace(Af, mass**2 * (1 - L))
    .replace(Bf, 1 - L + s / (2 * mass**2))
    .series(s, 0, 0)
    .to_expression()
    .expand()
    for expression in finite
]
assert (zero_invariant[0] - xi * (E("1/2") - L)).expand() == zero
assert (zero_invariant[1] + (3 + xi) * (1 - L) - 2).expand() == zero
on_shell = [
    expression.replace(s, mass**2)
    .replace(Af, mass**2 * (1 - L))
    .replace(Bf, 2 - L)
    .expand()
    for expression in finite
]
assert (sum(on_shell, zero) + 4 - 3 * L).expand() == zero
assert (sum(expected_poles, zero) + 3).expand() == zero
assert sum(on_shell, zero).derivative(xi) == zero

# Shared OneLOop callbacks. Preserve the D-dependent coefficients through the
# Laurent expansion: the finite rational terms are -xi and +2 respectively.
a_coefficients = oneloop.master_coefficients(S("oneloopmaster::A0")(mass**2, scale2))
b_coefficients = oneloop.master_coefficients(
    S("oneloopmaster::B0")(s, mass**2, 0, scale2)
)
finite_oneloop = [
    expression.replace(Af, a_coefficients[0]).replace(Bf, b_coefficients[0])
    for expression in finite
]
numeric_points = [
    (-10.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
    (-0.001, 1.0, 4.0),
    (0.001, 1.0, 1.0),
    (0.3, 1.0, 3.0),
    (0.99, 1.0, 1.0),
    (1.0, 1.0, 1.0),
    (1.01, 1.0, 1.0),
    (3.0, 1.0, 7.0),
    (100.0, 1.0, 1.0),
    (-3.0, 2.0, 7.0),
    (2.0, 2.0, 0.25),
    (4.0, 2.0, 3.0),
    (10.0, 2.0, 7.0),
    (1.0, 0.5, 2.0),
]
for invariant, mass_value, scale in numeric_points:
    ratio = invariant / mass_value**2
    mass_log = log(mass_value**2 / scale)
    log_cut = log(abs(1 - ratio)) - (1j * pi if ratio > 1 else 0j) if ratio != 1 else 0j
    b_reference = 2 - mass_log + (1 - ratio) / ratio * log_cut
    a_reference = mass_value**2 * (1 - mass_log)
    for gauge in (0.0, 1.0, 3.0, -0.5):
        point = {s: invariant, mass: mass_value, scale2: scale, xi: gauge}
        assert abs(complex(a_coefficients[1].evaluate(point)) - mass_value**2) < 1e-12
        assert abs(complex(b_coefficients[1].evaluate(point)) - 1) < 1e-12
        assert abs(complex(a_coefficients[2].evaluate(point))) < 1e-12
        assert abs(complex(b_coefficients[2].evaluate(point))) < 1e-12
        assert abs(complex(a_coefficients[0].evaluate(point)) - a_reference) < 3e-12
        assert abs(complex(b_coefficients[0].evaluate(point)) - b_reference) < 3e-12
        references = [
            gauge
            * ((invariant + mass_value**2) * b_reference - a_reference)
            / invariant
            - gauge,
            -(3 + gauge) * b_reference + 2,
        ]
        for expression, reference in zip(finite_oneloop, references, strict=True):
            evaluated = complex(expression.evaluate(point))
            # Separate A0/B0 evaluation amplifies scalar error by m²/|s|.
            tolerance = 3e-12 * max(
                1.0, abs(gauge) * (1 + mass_value**2 / abs(invariant)), 3 + abs(gauge)
            )
            assert abs(evaluated - reference) < tolerance, (
                point,
                evaluated,
                reference,
                tolerance,
            )
        if ratio > 1:
            assert complex(b_coefficients[0].evaluate(point)).imag > 0
# The scalar master at s=0 is regular; the exact projected expression uses its
# symbolic limit above. Check the separate master value without a 0/0 projector.
assert (
    abs(
        complex(b_coefficients[0].evaluate({s: 0.0, mass: 2.0, scale2: 7.0}))
        - (1 - log(4 / 7))
    )
    < 1e-12
)
print(
    "Electron self-energy: native phase, exact D, UV and finite terms, "
    "Landau and on-shell identities, 60 OneLOop points and s=0 limit passed"
)
