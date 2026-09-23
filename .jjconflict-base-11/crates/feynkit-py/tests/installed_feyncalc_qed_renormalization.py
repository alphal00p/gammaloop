"""Generated one-loop QED MS and MSbar counterterms in a symbolic covariant gauge.

Reference: https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Renormalization
The graph UV expansion, traces, tensor projection and native IBP determine the
bare poles. The tadpole pole and local counterterm operators are stated inputs.
Bare Z factors generate six local CT diagrams and their matching matrix through
the ordinary model/generator APIs; no subtraction forest is constructed.

The massless reference uses direct auxiliary-mass replacement (n=0 IRR),
validated independently against the shared full UV expansion and its explicit
mass-compensation terms. Its five constants reuse the generated CT matrix.
https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/RenormalizationMassless

MS/MSbar matching uses the same generated CT amplitudes and OneLOop measure.
https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Renormalization2
"""

import json
from pathlib import Path

from symbolica import E, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
D, epsilon, M, mUV, coordinate, integral, xi, s, Nf = S(
    "qed_ren::D",
    "qed_ren::eps",
    "qed_ren::M",
    "qed_ren::mUV",
    "qed_ren::d0",
    "qed_ren::I",
    "qed_ren::xi",
    "qed_ren::s",
    "qed_ren::Nf",
)
specification = json.loads(model.to_json())
for propagator in specification["propagators"]:
    if propagator["particle"] == "a":
        propagator["numerator"] = (
            "-1𝑖*(UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))-(1-qed_ren::xi)*UFO::P(UFO::idx(1,1))*UFO::P(UFO::idx(1,2))/spenso::dot(UFO::P(spenso::mink(4)),UFO::P(spenso::mink(4))))"
        )
model = hep.Model.from_json(json.dumps(specification))
electron, photon = (model.particle_by_pdg(pdg) for pdg in (11, 22))
vertices = [
    vertex
    for vertex in model.vertex_rules
    if sorted(vertex.particles)
    == sorted([electron.antiname, electron.name, photon.name])
]
assert len(vertices) == 1
K, P, mass, charge = S("gammalooprs::K", "gammalooprs::P", "UFO::Me", "UFO::ee")
mink, bis, gamma, metric = S(
    "spenso::mink", "spenso::bis", "spenso::gamma", "spenso::g"
)
index, dim, wave, mu, nu = S(
    "qed_ren::index_", "qed_ren::dim_", "qed_ren::wave_", "qed_ren::mu", "qed_ren::nu"
)
den, edge_, mom_, mass_, quad_ = S(
    "gammalooprs::denom",
    "qed_ren::edge_",
    "qed_ren::mom_",
    "qed_ren::mass_",
    "qed_ren::quad_",
)
ordering, value = S(
    "feynkit_generator_factor::ExternalFermionOrderingSign", "qed_ren::value_"
)
kinematics = hep.Kinematics(D, momenta=[K(0), P(0)]).with_scalar_product(P(0), P(0), s)
vacuum = hep.Kinematics(D, momenta=[K(0)])
family = hep.IntegralFamily(
    [K(0)], [], [vacuum.scalar_product(K(0), K(0)) - M], kinematics=vacuum
)
reducer = hep.TensorReducer(D).with_integrated_vector(K(0, mink(D)))
gmunu = metric(mink(D, mu), mink(D, nu))
ppmunu = P(0, mink(D, mu)) * P(0, mink(D, nu))
zero, one = E("0"), E("1")
parts, diagrams = {}, {}
for kind, incoming, outgoing, loops in [
    ("tree", [electron], [photon, electron], 0),
    ("electron", [electron], [electron], 1),
    ("photon", [photon], [photon], 1),
    ("vertex", [electron], [photon, electron], 1),
]:
    generated = hep.Generator(model).generate(
        hep.Process.amplitude(incoming, outgoing).with_loop_count(loops, loops),
        max_vertices=len(incoming) + len(outgoing) - 2 + 2 * loops,
        maximum_bridges=0,
        vertex_allow=vertices,
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    diagrams[kind] = diagram
    ports = {}
    for edge in diagram.external_edges:
        rep = (
            mink
            if kind == "photon"
            or (kind in ("tree", "vertex") and edge.external_index == 1)
            else bis
        )
        ports[edge.external_index] = dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, rep(4, index)), max_level=0
                )
            )
        )[index]
    numerator = model.expand_couplings(diagram.numerator_expression().to_expression())
    if loops:
        numerator = diagram.momentum_basis().route_expression(
            diagram.uv_expansion(mUV, numerator=numerator).to_expression()
        )
    numerator = (
        numerator.replace(mink(dim, index), mink(D, index))
        .replace(mink(dim), mink(D))
        .replace(mUV**2, M)
    )
    if kind == "photon":
        numerator = numerator.replace(mink(D, ports[0]), mink(D, mu)).replace(
            mink(D, ports[1]), mink(D, nu)
        )
    pattern = den(edge_, mom_, mass_, quad_)
    for match in list(numerator.match(pattern)):
        values = dict(match)
        formal = family.rewrite_numerator(values[quad_], [coordinate])
        assert formal == coordinate
        numerator = numerator.replace(
            den(values[edge_], values[mom_], values[mass_], values[quad_]), formal
        )
    factor = (
        diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    if kind == "electron":
        raw = diagram.overall_factor_expression()
        removed = (raw / raw.replace(ordering(value), one)).replace(
            ordering(value), value
        )
        assert removed == -one
        factor /= removed
        probes = [
            (
                "electron_p",
                gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
                * P(0, mink(D, mu))
                / (4 * s),
            ),
            ("electron_m", metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass)),
        ]
    elif kind in ("vertex", "tree"):
        probes = [
            (
                kind,
                gamma(bis(4, ports[0]), bis(4, ports[2]), mink(D, ports[1])) / (4 * D),
            )
        ]
    else:
        probes = [(kind, one)]
    for label, projector in probes:
        trace = (
            TensorExpression((numerator * projector).expand())
            .simplify_gamma()
            .expand()
            .to_dots()
            .to_expression()
        )
        scalar = family.rewrite_numerator(
            kinematics.apply(reducer.reduce(trace)), [coordinate]
        )
        scalar = (scalar * factor).together().expand()
        terms = []
        for monomial, coefficient in scalar.coefficient_list(coordinate):
            power = -int(
                (monomial.derivative(coordinate) * coordinate / monomial).together()
            )
            assert monomial == coordinate ** (-power)
            assert not coefficient.matches(K(S("qed_ren::args___")))
            terms.append(([power], coefficient))
        parts[label] = terms

tree = parts.pop("tree")[0][1]
assert tree == Symbol.I * charge
solution = hep.IBPFamily(family, name="qed_one_loop").reduce_laporta(
    sorted({tuple(p) for terms in parts.values() for p, c in terms}), max_depth=2
)
assert solution.residuals == [[1]]
assert abs(complex(oneloop.A0(1.0, 1.0)[1]) - 1) < 1e-12
reduced, uv_poles = {}, {}
for label, terms in parts.items():
    expression = sum(
        (c * solution.reduce(p, integral=integral) for p, c in terms), zero
    ).together()
    reduced[label] = expression
    normalized = (
        expression * (Symbol.I / tree if label == "vertex" else one) / charge**2
    )
    if label == "photon":
        normalized *= Nf
    pole = (
        normalized.replace(integral(1), M / epsilon)
        .replace(D, 4 - 2 * epsilon)
        .series(epsilon, 0, -1)
        .to_expression()
        .expand()
    )
    pole = pole.replace(mink(4, index), mink(D, index))
    uv_poles[label] = pole
    assert pole.derivative(M).expand() == zero
    assert pole.derivative(mass).expand() == zero
    assert pole.coefficient(epsilon**-2) == zero
# Compare the four UV structures with the published symbolic-gauge result.
expected = [
    xi / epsilon,
    -(xi + 3) / epsilon,
    -4 * Nf * (s * gmunu - ppmunu) / (3 * epsilon),
    xi / epsilon,
]
for label, reference in zip(
    ("electron_p", "electron_m", "photon", "vertex"), expected, strict=True
):
    assert (uv_poles[label] - reference).together() == zero

# Counterterm structures follow the local kinetic, mass and vertex operators.
# Solve for coefficients of Z=1+a4*deltaZ. Their coupling combinations come from
# expanding the bare factors, while actual generated diagrams supply the matrix.
# M=mUV². The reference auxiliary-mass operator is M*(ZAm²-1)*A²/2:
# its linear coefficient is 2*deltaZAm, unlike an additive mass-squared shift.
a4, delta_psi, delta_m, delta_A, delta_xi, delta_e, delta_Am = S(
    "qed_ren::a4",
    "qed_ren::deltaZpsi",
    "qed_ren::deltaZm",
    "qed_ren::deltaZA",
    "qed_ren::deltaZxi",
    "qed_ren::deltaZe",
    "qed_ren::deltaZAm",
)
unknowns = [delta_psi, delta_m, delta_A, delta_xi, delta_e, delta_Am]
Zpsi, Zm, ZA, Zxi, Ze, ZAm = (one + a4 * delta for delta in unknowns)
external_ordering = diagrams["tree"].overall_factor_expression(evaluate=True)
assert external_ordering == -one
assert (
    diagrams["electron"].overall_factor_expression(evaluate=True) == external_ordering
)
local_tree = tree / external_ordering
assert local_tree == -Symbol.I * charge
specification = json.loads(model.to_json())
specification["orders"].append({"name": "CT", "expansion_order": 1, "hierarchy": 1})
ffv = next(
    structure
    for structure in specification["lorentz_structures"]
    if structure["name"] == vertices[0].lorentz_structures[0]
)
for label, particles, spins, lorentz, coupling, qed_order in [
    (
        "ee_kinetic",
        [electron.antiname, electron.name],
        [2, 2],
        # Normalized JSON uses the spinor order of the imported SM FFV rule.
        # UFO momenta are incoming: leg 2 carries the electron's momentum.
        "Gamma(dummy(1),idx(1,1),idx(1,2))*P(dummy(1),2)",
        Symbol.I * (Zpsi - one),
        2,
    ),
    (
        "ee_mass",
        [electron.antiname, electron.name],
        [2, 2],
        "Identity(idx(1,1),idx(1,2))",
        -Symbol.I * mass * (Zpsi * Zm - one),
        2,
    ),
    (
        "aa_kinetic",
        [photon.name] * 2,
        [3, 3],
        (
            "Metric(idx(1,1),idx(1,2))*P(dummy(1),1)*P(dummy(1),1)"
            "-P(idx(1,1),1)*P(idx(1,2),1)"
        ),
        -Symbol.I * (ZA - one),
        2,
    ),
    (
        "aa_gauge",
        [photon.name] * 2,
        [3, 3],
        "P(idx(1,1),1)*P(idx(1,2),1)",
        -Symbol.I * (ZA / Zxi - one) / xi,
        2,
    ),
    (
        "aa_auxmass",
        [photon.name] * 2,
        [3, 3],
        "Metric(idx(1,1),idx(1,2))",
        Symbol.I * M * (ZAm**2 - one),
        2,
    ),
    (
        "eea",
        list(vertices[0].particles),
        list(ffv["spins"]),
        ffv["structure"],
        local_tree * (Zpsi * Ze * ZA.sqrt() - one),
        3,
    ),
]:
    specification["lorentz_structures"].append(
        {"name": "CT_L_" + label, "spins": spins, "structure": lorentz}
    )
    specification["couplings"].append(
        {
            "name": "CT_GC_" + label,
            "expression": repr(coupling.series(a4, 0, 1).to_expression()),
            "orders": [["QED", qed_order], ["CT", 1]],
            "value": None,
        }
    )
    specification["vertex_rules"].append(
        {
            "name": "CT_" + label,
            "particles": particles,
            "color_structures": ["1"],
            "lorentz_structures": ["CT_L_" + label],
            "couplings": [["CT_GC_" + label]],
        }
    )
ct_model = hep.Model.from_json(json.dumps(specification))
ct_electron, ct_photon = (ct_model.particle_by_pdg(pdg) for pdg in (11, 22))
ct_vertices = [
    vertex for vertex in ct_model.vertex_rules if vertex.name.startswith("CT_")
]
assert len(ct_vertices) == 6
ct_diagrams, ct_coefficients = {}, {}
for kind, incoming, outgoing, count, qed_order in [
    ("electron", [ct_electron], [ct_electron], 2, 2),
    ("photon", [ct_photon], [ct_photon], 3, 2),
    ("vertex", [ct_electron], [ct_photon, ct_electron], 1, 3),
]:
    # These are tree topologies at perturbative CT order one. Bound vertices
    # explicitly because arbitrarily many two-point insertions add no loops.
    options = {
        "loops": 0,
        "max_vertices": 1,
        "maximum_bridges": None,
        "vertex_allow": ct_vertices,
        "self_energy": None,
        "tadpoles": None,
        "zero_snails": None,
        "numerator_grouping": None,
        "progress": None,
    }
    generated = ct_model.generate_diagrams(
        incoming, outgoing, coupling_orders={"QED": qed_order, "CT": 1}, **options
    )
    assert len(generated.diagrams) == count
    assert not ct_model.generate_diagrams(
        incoming, outgoing, coupling_orders={"CT": 0}, **options
    ).diagrams
    ct_diagrams[kind] = generated.diagrams
    for diagram in generated.diagrams:
        assert len(diagram.internal_edges) == 0
        assert diagram.symmetry_factor == one
        assert diagram.numerator_prefactor_expression() == one
        factor = diagram.overall_factor_expression(evaluate=True)
        assert factor == (one if kind == "photon" else external_ordering)
        ports = {}
        for edge in diagram.external_edges:
            rep = (
                mink
                if kind == "photon" or (kind == "vertex" and edge.external_index == 1)
                else bis
            )
            ports[edge.external_index] = dict(
                next(
                    diagram.projector_expression().match(
                        wave(edge.id, rep(4, index)), max_level=0
                    )
                )
            )[index]
        numerator = ct_model.expand_couplings(
            diagram.numerator_expression(in_lmb=True).to_expression()
        ).replace(mink(4, index), mink(D, index))
        if kind == "electron":
            probes = [
                (
                    "electron_p",
                    gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
                    * P(0, mink(D, mu))
                    / (4 * s),
                ),
                ("electron_m", metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass)),
            ]
            # Convert only the named external-state convention, matching the
            # bare electron insertion above; retain the native generated factor.
            normalization = Symbol.I * a4 * external_ordering
        elif kind == "photon":
            numerator = numerator.replace(mink(D, ports[0]), mink(D, mu)).replace(
                mink(D, ports[1]), mink(D, nu)
            )
            probes = [(kind, one)]
            normalization = Symbol.I * a4
        else:
            probes = [
                (
                    kind,
                    gamma(bis(4, ports[0]), bis(4, ports[2]), mink(D, ports[1]))
                    / (4 * D),
                )
            ]
            normalization = a4 * tree
        for label, projector in probes:
            trace = (
                TensorExpression((numerator * projector).expand())
                .simplify_gamma()
                .expand()
                .to_dots()
                .to_expression()
            )
            coefficient = (
                (kinematics.apply(trace) * factor / normalization).together().expand()
            )
            ct_coefficients[label] = ct_coefficients.get(label, zero) + coefficient

ct_photon_g = ct_coefficients["photon"].coefficient(gmunu)
ct_photon_pp = ct_coefficients["photon"].coefficient(ppmunu)
assert (
    ct_coefficients["photon"] - ct_photon_g * gmunu - ct_photon_pp * ppmunu
).together() == zero
assert ct_photon_g.derivative(s).derivative(s) == zero
ct_rows = [
    ct_coefficients["electron_p"],
    ct_coefficients["electron_m"],
    ct_photon_g.coefficient(s),
    xi * ct_photon_pp,
    ct_coefficients["vertex"],
    ct_photon_g.replace(s, zero) / M,
]
# Unknown order: deltaZpsi, deltaZm, deltaZA, deltaZxi, deltaZe, deltaZAm.
# Extract every matrix entry from generated amplitudes; check linearity so no
# constant or nonlinear term can be silently dropped by differentiation.
ct_entries = [row.derivative(delta).together() for row in ct_rows for delta in unknowns]
assert all(
    entry.derivative(delta) == zero for entry in ct_entries for delta in unknowns
)
for position, row in enumerate(ct_rows):
    assert (
        row
        - sum(
            (ct_entries[6 * position + j] * delta for j, delta in enumerate(unknowns)),
            zero,
        )
    ).together() == zero
ct_matrix = Matrix.from_linear(6, 6, ct_entries)
photon_g = uv_poles["photon"].coefficient(gmunu)
photon_pp = uv_poles["photon"].coefficient(ppmunu)
rhs = Matrix.vec(
    [
        -uv_poles["electron_p"],
        -uv_poles["electron_m"],
        -photon_g.coefficient(s),
        -xi * photon_pp,
        -uv_poles["vertex"],
        -photon_g.replace(s, zero) / M,
    ]
)
deltas = ct_matrix.solve(rhs)
for row, reference in enumerate(
    [
        -xi / epsilon,
        -3 / epsilon,
        -4 * Nf / (3 * epsilon),
        -4 * Nf / (3 * epsilon),
        2 * Nf / (3 * epsilon),
        zero,
    ]
):
    assert (deltas[row, 0].to_expression() - reference).together() == zero
assert (deltas[0, 0].to_expression() + uv_poles["vertex"]).together() == zero
assert (
    deltas[4, 0].to_expression() + deltas[2, 0].to_expression() / 2
).together() == zero
residual = ct_matrix * deltas - rhs
assert all(residual[row, 0].to_expression().together() == zero for row in range(6))
print(
    "Generated QED one-loop: symbolic gauge, four IBP targets, six generated CT diagrams and Ward identities passed",
    solution.stats,
)

# The massless gallery uses a different, explicit IRR prescription: first set
# the physical electron mass to zero, then FCLoopAddAuxiliaryMass[..., -M, 0]
# directly replaces every massless denominator q² by q²-M. FourSeries then
# Taylor-expands external momenta at fixed M through the UV degree. The n=0
# choice omits the compensating mass terms of a UV-preserving rearrangement.
# Reference: https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/RenormalizationMassless
massless_specification = json.loads(model.to_json())
for particle in massless_specification["particles"]:
    if abs(particle["pdg_code"]) == 11:
        particle["mass"] = "ZERO"
massless_model = hep.Model.from_json(json.dumps(massless_specification))
massless_electron, massless_photon = (
    massless_model.particle_by_pdg(pdg) for pdg in (11, 22)
)
massless_vertices = [
    vertex
    for vertex in massless_model.vertex_rules
    if sorted(vertex.particles)
    == sorted(
        [massless_electron.antiname, massless_electron.name, massless_photon.name]
    )
]
assert len(massless_vertices) == 1
Q, dot, massless_scale, massless_a, massless_b, massless_args = S(
    "gammalooprs::Q",
    "spenso::dot",
    "qed_ren::massless_scale",
    "qed_ren::massless_a_",
    "qed_ren::massless_b_",
    "qed_ren::massless_args___",
)
massless_kinematics = hep.Kinematics(D, momenta=[K(0), P(0), P(1)]).with_scalar_product(
    P(0), P(0), s
)
massless_loop_square = vacuum.scalar_product(K(0), K(0))
massless_Qg, massless_Qpp = S("qed_ren::massless_Qg", "qed_ren::massless_Qpp")
massless_parts, massless_targets, massless_scalars, massless_diagrams = (
    {},
    set(),
    {},
    {},
)
massless_components, massless_component_scalars = {}, {}
massless_photon_coordinate = S("qed_ren::massless_photon_coordinate")
for kind, incoming, outgoing in [
    ("electron_p", [massless_electron], [massless_electron]),
    ("photon", [massless_photon], [massless_photon]),
    ("vertex", [massless_electron], [massless_photon, massless_electron]),
]:
    generated = massless_model.generate_diagrams(
        incoming,
        outgoing,
        loops=1,
        max_vertices=len(incoming) + len(outgoing),
        maximum_bridges=0,
        vertex_allow=massless_vertices,
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    massless_diagrams[kind] = diagram
    ports = {}
    for edge in diagram.external_edges:
        rep = (
            mink
            if kind == "photon" or (kind == "vertex" and edge.external_index == 1)
            else bis
        )
        ports[edge.external_index] = dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, rep(4, index)), max_level=0
                )
            )
        )[index]
    numerator = massless_model.expand_couplings(
        diagram.numerator_expression().to_expression()
    ).replace(mass, zero)
    # Preserve the primitive denominators used by the n=0 prescription.
    # Combining the Feynman term over a squared photon denominator first
    # would turn 1/q² into q²/(q²-M)² after massification and shift finite
    # terms. Only the longitudinal remainder has a genuine extra 1/q².
    photon_edges = [
        edge
        for edge in diagram.internal_edges
        if edge.particle_name == massless_photon.name
    ]
    assert len(photon_edges) <= 1
    components = [("fermion_loop", numerator, {})]
    if photon_edges:
        photon_edge = photon_edges[0]
        photon_square = dot(Q(photon_edge.id, mink(4)), Q(photon_edge.id, mink(4)))
        feynman_numerator = numerator.replace(xi, one)
        longitudinal_numerator = (
            ((numerator - feynman_numerator) * photon_square).together().expand()
        )
        assert longitudinal_numerator.replace(xi, one).expand() == zero
        assert (
            feynman_numerator + longitudinal_numerator / photon_square - numerator
        ).together() == zero
        components = [
            ("feynman", feynman_numerator, {}),
            ("longitudinal", longitudinal_numerator, {photon_edge.id: 2}),
        ]
        # At xi=1 the primitive numerator is unchanged and its photon has
        # the default single propagator; no q² has been multiplied into it.
        assert components[0][1] == numerator.replace(xi, one)
        assert components[0][2] == {}
    massless_components[kind] = components
    factor = (
        diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
    )
    if kind == "electron_p":
        assert factor == external_ordering
        factor /= external_ordering
        projector = (
            gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
            * P(0, mink(D, mu))
            / (4 * s)
        )
    elif kind == "vertex":
        projector = gamma(bis(4, ports[0]), bis(4, ports[2]), mink(D, ports[1])) / (
            4 * D
        )
    else:
        projector = one
    for component, numerator, powers in components:
        numerator = numerator.replace(mink(dim, index), mink(D, index)).replace(
            mink(dim), mink(D)
        )
        if kind == "photon":
            numerator = numerator.replace(mink(D, ports[0]), mink(D, mu)).replace(
                mink(D, ports[1]), mink(D, nu)
            )

        # Retain the shared full UV result as an independent contrast. Freeze its
        # tagged massive denominators before extracting the zeroth explicit mUV
        # coefficient; M in the family is held fixed, never replaced by zero.
        expanded = (
            diagram.momentum_basis()
            .route_expression(
                diagram.uv_expansion(
                    mUV, numerator=numerator, edge_powers=powers
                ).to_expression()
            )
            .replace(mink(dim, index), mink(D, index))
            .replace(mink(dim), mink(D))
        )
        pattern = den(edge_, mom_, mass_, quad_)
        for match in list(expanded.match(pattern)):
            values = dict(match)
            formal = family.rewrite_numerator(
                values[quad_].replace(mUV**2, M), [coordinate]
            )
            assert formal == coordinate
            expanded = expanded.replace(
                den(values[edge_], values[mom_], values[mass_], values[quad_]),
                coordinate,
            )
        assert not expanded.matches(pattern)
        variants = {
            "full_uv": expanded.replace(mUV**2, M),
            "selective_uv": expanded.replace(mUV, zero),
        }

        # Implement the primary prescription literally and independently: replace
        # massless graph denominators by q²-M, then Taylor-expand external momenta.
        routed = diagram.momentum_basis().route_expression(numerator)
        traced = (
            TensorExpression((routed * projector).expand())
            .simplify_gamma()
            .expand()
            .to_dots()
            .to_expression()
        )
        denominator = (
            diagram.denominator_expression(in_lmb=True, edge_powers=powers)
            .to_expression()
            .replace(mink(dim, index), mink(D, index))
            .replace(mink(dim), mink(D))
        )
        for match in list(denominator.match(pattern)):
            values = dict(match)
            assert values[mass_] == zero
            if photon_edges and values[edge_] == photon_edge.id:
                # Verify the actual graph denominator power independently of
                # the override map, before any auxiliary-mass replacement.
                photon_tag = den(
                    values[edge_], values[mom_], values[mass_], values[quad_]
                )
                checked_denominator = denominator.replace(
                    photon_tag, massless_photon_coordinate
                )
                actual_power = (
                    checked_denominator.derivative(massless_photon_coordinate)
                    * massless_photon_coordinate
                    / checked_denominator
                ).together()
                assert actual_power == (1 if component == "feynman" else 2)
            denominator = denominator.replace(
                den(values[edge_], values[mom_], values[mass_], values[quad_]),
                values[quad_] - M,
            )
        direct = massless_kinematics.apply(traced / denominator).replace(
            massless_loop_square, coordinate + M
        )
        direct = direct.replace_multiple(
            [
                Replacement(
                    dot(K(0, mink(D)), P(massless_a, mink(D))),
                    massless_scale * dot(K(0, mink(D)), P(massless_a, mink(D))),
                ),
                Replacement(
                    dot(P(massless_a, mink(D)), K(0, mink(D))),
                    massless_scale * dot(P(massless_a, mink(D)), K(0, mink(D))),
                ),
                Replacement(
                    dot(P(massless_a, mink(D)), P(massless_b, mink(D))),
                    massless_scale**2
                    * dot(P(massless_a, mink(D)), P(massless_b, mink(D))),
                ),
                Replacement(
                    P(massless_a, mink(D, index)),
                    massless_scale * P(massless_a, mink(D, index)),
                ),
                Replacement(s, massless_scale**2 * s),
            ]
        )
        # The electron's pslash/(4p²) projector lowers the external degree by one,
        # so projected order0 retains its degree1 self-energy. Photon needs order2;
        # the logarithmically divergent vertex needs only order0.
        variants["direct_irr"] = (
            direct.series(massless_scale, 0, 2 if kind == "photon" else 0)
            .to_expression()
            .replace(massless_scale, one)
        )
        for scheme, expression in variants.items():
            trace = (
                expression
                if scheme == "direct_irr"
                else TensorExpression((expression * projector).expand())
                .simplify_gamma()
                .expand()
                .to_dots()
                .to_expression()
            )
            scalar = family.rewrite_numerator(
                massless_kinematics.apply(reducer.reduce(trace)), [coordinate]
            )
            normalization = (
                factor
                / charge**2
                * (Symbol.I / tree if kind == "vertex" else one)
                * (Nf if kind == "photon" else one)
            )
            scalar = (normalization * scalar).together().expand()
            massless_component_scalars[scheme, kind, component] = scalar
            massless_scalars[scheme, kind] = (
                massless_scalars.get((scheme, kind), zero) + scalar
            ).expand()
        assert (
            massless_component_scalars["direct_irr", kind, component]
            - massless_component_scalars["selective_uv", kind, component]
        ).together() == zero
    for scheme in ("full_uv", "selective_uv", "direct_irr"):
        scalar = massless_scalars[scheme, kind]
        terms = []
        for monomial, coefficient in scalar.coefficient_list(coordinate):
            power = -int(
                (monomial.derivative(coordinate) * coordinate / monomial).together()
            )
            assert monomial == coordinate ** (-power)
            assert not coefficient.matches(K(massless_args))
            scalar_coefficient = coefficient.replace(gmunu, massless_Qg).replace(
                ppmunu, massless_Qpp
            )
            assert not scalar_coefficient.matches(P(massless_args))
            assert coefficient.derivative(coordinate) == zero
            massless_targets.add((power,))
            terms.append(([power], coefficient))
        massless_parts[scheme, kind] = terms
    # This is an exact integrand-level check after vacuum tensor reduction,
    # not merely agreement of the final UV pole with a supplied formula.
    assert (
        massless_scalars["direct_irr", kind] - massless_scalars["selective_uv", kind]
    ).together() == zero

massless_solution = hep.IBPFamily(family, name="qed_massless_irr").reduce_laporta(
    [list(target) for target in sorted(massless_targets)], max_depth=2
)
assert massless_solution.residuals == [[1]]
massless_integrated, massless_uv_poles = {}, {}
for label, terms in massless_parts.items():
    expression = sum(
        (
            coefficient * massless_solution.reduce(power, integral=integral)
            for power, coefficient in terms
        ),
        zero,
    ).together()
    massless_integrated[label] = expression
    massless_uv_poles[label] = (
        expression.replace(gmunu, massless_Qg)
        .replace(ppmunu, massless_Qpp)
        .replace(integral(1), M / epsilon)
        .replace(D, 4 - 2 * epsilon)
        .series(epsilon, 0, -1)
        .to_expression()
        .expand()
        .replace(massless_Qg, gmunu)
        .replace(massless_Qpp, ppmunu)
    )
    assert massless_uv_poles[label].coefficient(epsilon**-2) == zero
for scheme in ("direct_irr", "selective_uv", "full_uv"):
    assert (massless_uv_poles[scheme, "electron_p"] - xi / epsilon).together() == zero
    assert (massless_uv_poles[scheme, "vertex"] - xi / epsilon).together() == zero
    reference = (
        Nf
        * (
            -4 * (s * gmunu - ppmunu) / 3
            + (4 * M * gmunu if scheme != "full_uv" else zero)
        )
        / epsilon
    )
    assert (massless_uv_poles[scheme, "photon"] - reference).together() == zero
# Each full-UV propagator retains -M/(K²-M)² at this order. The sum of the
# two photon-bubble compensation terms reduces to -(D-2)²*Nf*A0(M)*g:
# its -4*M*Nf*g/epsilon pole cancels the direct IRR auxiliary photon mass.
assert (
    massless_integrated["full_uv", "photon"]
    - massless_integrated["direct_irr", "photon"]
    + (D - 2) ** 2 * Nf * integral(1) * gmunu
).together() == zero

# Reuse the already generated local CT amplitudes and matching matrix. There
# is no physical mass operator in the massless calculation: remove the
# electron_m row and deltaZm column, preserving all other generated entries.
massless_indices = [0, 2, 3, 4, 5]
massless_unknowns = [unknowns[position] for position in massless_indices]
massless_ct_matrix = Matrix.from_linear(
    5,
    5,
    [
        ct_matrix[row, column].to_expression()
        for row in massless_indices
        for column in massless_indices
    ],
)
assert massless_ct_matrix[4, 4].to_expression() == 2
assert all(
    ct_rows[position].derivative(delta_m) == zero for position in massless_indices
)
massless_photon_g = massless_uv_poles["direct_irr", "photon"].coefficient(gmunu)
massless_photon_pp = massless_uv_poles["direct_irr", "photon"].coefficient(ppmunu)
massless_rhs = Matrix.vec(
    [
        -massless_uv_poles["direct_irr", "electron_p"],
        -massless_photon_g.coefficient(s),
        -xi * massless_photon_pp,
        -massless_uv_poles["direct_irr", "vertex"],
        -massless_photon_g.replace(s, zero) / M,
    ]
)
massless_deltas = massless_ct_matrix.solve(massless_rhs)
for row, reference in enumerate(
    [
        -xi / epsilon,
        -4 * Nf / (3 * epsilon),
        -4 * Nf / (3 * epsilon),
        2 * Nf / (3 * epsilon),
        -2 * Nf / epsilon,
    ]
):
    assert (massless_deltas[row, 0].to_expression() - reference).together() == zero
assert (
    massless_deltas[0, 0].to_expression() + massless_uv_poles["direct_irr", "vertex"]
).together() == zero
assert (
    massless_deltas[3, 0].to_expression() + massless_deltas[1, 0].to_expression() / 2
).together() == zero
massless_residual = massless_ct_matrix * massless_deltas - massless_rhs
assert all(
    massless_residual[row, 0].to_expression().together() == zero for row in range(5)
)
# Cancel the complete generated tensor/Dirac structures, including the
# auxiliary-mass term; the matrix solution is not the only validation.
massless_ct_rules = [
    Replacement(delta, massless_deltas[position, 0].to_expression())
    for position, delta in enumerate(massless_unknowns)
]
for kind in ("electron_p", "photon", "vertex"):
    assert (
        massless_uv_poles["direct_irr", kind]
        + ct_coefficients[kind].replace_multiple(massless_ct_rules)
    ).together() == zero
print(
    "Massless QED: literal n=0 IRR, full-UV contrast, five generated CT coefficients and Ward identities passed",
    massless_solution.stats,
)

# Compare the reference's MS and MSbar conventions at the same scale. OneLOop
# divides by rGamma; (4*pi)^epsilon*rGamma=1+cDelta*epsilon+O(epsilon^2),
# with cDelta=log(4*pi)-EulerGamma. Keep that real constant symbolic through
# exact simplification. No finite part is inferred from a UV-expanded graph.
cDelta = S("qed_ren::cDelta")
assert (-2 / (D - 4)).replace(D, 4 - 2 * epsilon) == 1 / epsilon
counterterms_by_scheme = {
    "MS": [deltas[row, 0].to_expression().together() for row in range(6)],
    "MSbar": [
        (epsilon * deltas[row, 0].to_expression() * (1 / epsilon + cDelta)).together()
        for row in range(6)
    ],
}
massless_counterterms_by_scheme = {
    "MS": [massless_deltas[row, 0].to_expression().together() for row in range(5)],
    "MSbar": [
        (
            epsilon * massless_deltas[row, 0].to_expression() * (1 / epsilon + cDelta)
        ).together()
        for row in range(5)
    ],
}
for pole_set, delta_names, schemes in [
    (uv_poles, unknowns, counterterms_by_scheme),
    (
        {
            kind: massless_uv_poles["direct_irr", kind]
            for kind in ("electron_p", "photon", "vertex")
        },
        massless_unknowns,
        massless_counterterms_by_scheme,
    ),
]:
    for scheme, constants in schemes.items():
        for label, pole in pole_set.items():
            generated_ct = ct_coefficients[label]
            for delta, constant in zip(delta_names, constants, strict=True):
                generated_ct = generated_ct.replace(delta, constant)
            # These are the complete generated scalar/tensor CT coefficients,
            # retaining the external convention of the corresponding loop sector.
            converted_uv = pole * (1 + epsilon * cDelta if scheme == "MSbar" else one)
            assert (generated_ct + converted_uv).together() == zero

# Actual OneLOop values check the measure conversion independently of the
# counterterm solve. The remaining cDelta*P in MS is a finite scheme shift.
for master in [
    S("oneloopmaster::A0")(4, 7),
    S("oneloopmaster::B0")(-3, 4, 0, 7),
]:
    finite_master, pole_master, double_master = oneloop.get_expression(master)
    assert double_master == zero
    physical_master = (
        ((one + epsilon * cDelta) * (pole_master / epsilon + finite_master))
        .series(epsilon, 0, 0)
        .to_expression()
    )
    assert (
        physical_master - pole_master / epsilon - finite_master - cDelta * pole_master
    ).together().expand() == zero
    assert (
        physical_master - pole_master * (1 / epsilon + cDelta) - finite_master
    ).together().expand() == zero

print(
    "Generated QED CT amplitudes: MS/MSbar conversion and OneLOop measure checks passed"
)
