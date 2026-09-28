"""Generated scalar/pseudoscalar Yukawa renormalization in D=4-2eps.

References: FeynCalc YukawaS and YukawaPS / OneLoop / Renormalization.
The shared model, graph UV expansion, Spenso traces, tensor projection and
RustRed IBP calculate all two-, three- and four-point ultraviolet poles.
"""

import copy
import json

from symbolica import E, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import TensorExpression

# Reuse the built-in quartic scalar and Standard Model Dirac propagators.
# The only extra interaction is -g phi psibar psi, or -i g phi psibar gamma5 psi.
specification = json.loads(hep.Model.phi4().to_json())
standard = json.loads(hep.Model.standard_model().to_json())
specification["particles"] += [
    p for p in standard["particles"] if abs(p["pdg_code"]) == 11
]
fermion_names = {p["name"] for p in specification["particles"] if p["spin"] == 2}
specification["propagators"] += [
    p for p in standard["propagators"] if p["particle"] in fermion_names
]
specification["parameters"] += [p for p in standard["parameters"] if p["name"] == "Me"]
parameter = copy.deepcopy(specification["parameters"][2])
parameter.update(name="g", lhacode=[3])
specification["parameters"].append(parameter)
specification["orders"].append(
    {"name": "YUKAWA", "expansion_order": 99, "hierarchy": 1}
)

D, eps, M, mUV, coordinate, integral, s = S(
    "yukawa::D",
    "yukawa::eps",
    "yukawa::M",
    "yukawa::mUV",
    "yukawa::x",
    "yukawa::I",
    "yukawa::s",
)
K, P, mass, scalar_mass, g, lam = S(
    "gammalooprs::K",
    "gammalooprs::P",
    "UFO::Me",
    "UFO::mass",
    "UFO::g",
    "UFO::lam",
)
mink, bis, gamma, metric = S(
    "spenso::mink", "spenso::bis", "spenso::gamma", "spenso::g"
)
index, dim, wave, mu = S(
    "yukawa::index_", "yukawa::dim_", "yukawa::wave_", "yukawa::mu"
)
den, edge_, mom_, mass_, quad_ = S(
    "gammalooprs::denom",
    "yukawa::edge_",
    "yukawa::mom_",
    "yukawa::mass_",
    "yukawa::quad_",
)
zero, one = E("0"), E("1")
kinematics = hep.Kinematics(D, momenta=[K(0), P(0)]).with_scalar_product(P(0), P(0), s)
vacuum = hep.Kinematics(D, momenta=[K(0)])
family = hep.IntegralFamily(
    [K(0)], [], [vacuum.scalar_product(K(0), K(0)) - M], kinematics=vacuum
)
master_reduction = oneloop.reduce(family, [1])
master = master_reduction.terms[0][1].to_expression(one)
master_pole = oneloop.get_expression(master, coefficient=-1)
assert master_pole == M
reducer = hep.TensorReducer(D).with_integrated_vector(K(0, mink(D)))
options = {
    "maximum_bridges": 0,
    "self_energy": None,
    "tadpoles": None,
    "zero_snails": None,
    "numerator_grouping": None,
    "progress": None,
}
h, dpsi, dm, dphi, dM, dg, dlam = S(
    "yukawa::h",
    "yukawa::dpsi",
    "yukawa::dm",
    "yukawa::dphi",
    "yukawa::dM",
    "yukawa::dg",
    "yukawa::dlam",
)
unknowns = [dpsi, dm, dphi, dM, dg, dlam]
Zpsi, Zm, Zphi, ZM, Zg, Zlam = (1 + h * x for x in unknowns)
models, all_diagrams, all_poles, all_reductions = {}, {}, {}, {}
all_counterterms, all_matrices, all_constants = {}, {}, {}
for variant, lorentz, local_coupling in [
    ("Scalar", "Identity(idx(1,1),idx(1,2))", -Symbol.I * g),
    ("Pseudoscalar", "Gamma5(idx(1,1),idx(1,2))", g),
]:
    definition = copy.deepcopy(specification)
    definition["name"] = "yukawa_" + variant.lower()
    definition["lorentz_structures"].append(
        {"name": "YUKAWA", "spins": [2, 2, 1], "structure": lorentz}
    )
    definition["couplings"].append(
        {
            "name": "YUKAWA",
            "expression": repr(local_coupling),
            "orders": [["YUKAWA", 1]],
            "value": None,
        }
    )
    definition["vertex_rules"].append(
        {
            "name": "YUKAWA",
            "particles": ["e+", "e-", "phi"],
            "color_structures": ["1"],
            "lorentz_structures": ["YUKAWA"],
            "couplings": [["YUKAWA"]],
        }
    )
    model = hep.Model.from_json(json.dumps(definition))
    models[variant] = model
    fermion, scalar = model.particle("e-"), model.particle("phi")
    # Local operators generate the same six counterterms in both theories.
    # The scalar mass is squared; the fermion mass renormalizes linearly.
    ct_definition = json.loads(model.to_json())
    ct_definition["orders"].append({"name": "CT", "expansion_order": 1, "hierarchy": 1})
    for label, particles, structure, factor in [
        (
            "fermion_kinetic",
            ["e+", "e-"],
            "Gamma(dummy(1),idx(1,1),idx(1,2))*P(dummy(1),2)",
            Symbol.I * (Zpsi - 1),
        ),
        (
            "fermion_mass",
            ["e+", "e-"],
            "Identity(idx(1,1),idx(1,2))",
            -Symbol.I * mass * (Zpsi * Zm - 1),
        ),
        (
            "scalar_kinetic",
            ["phi"] * 2,
            "P(dummy(1),1)*P(dummy(1),1)",
            Symbol.I * (Zphi - 1),
        ),
        ("scalar_mass", ["phi"] * 2, "1", -Symbol.I * scalar_mass**2 * (Zphi * ZM - 1)),
        (
            "vertex",
            ["e+", "e-", "phi"],
            lorentz,
            local_coupling * (Zpsi * Zg * Zphi.sqrt() - 1),
        ),
        ("quartic", ["phi"] * 4, "1", -Symbol.I * lam * (Zlam * Zphi**2 - 1)),
    ]:
        name = "CT_" + label
        ct_definition["lorentz_structures"].append(
            {
                "name": name,
                "spins": [model.particle(p).spin for p in particles],
                "structure": structure,
            }
        )
        ct_definition["couplings"].append(
            {
                "name": name,
                "expression": repr(factor.series(h, 0, 1).to_expression()),
                "orders": [["CT", 1]],
                "value": None,
            }
        )
        ct_definition["vertex_rules"].append(
            {
                "name": name,
                "particles": particles,
                "color_structures": ["1"],
                "lorentz_structures": [name],
                "couplings": [[name]],
            }
        )
    ct_model = hep.Model.from_json(json.dumps(ct_definition))
    stage_parts, stage_diagrams = {}, {}
    for stage, stage_model, loops in [("bare", model, 1), ("ct", ct_model, 0)]:
        parts, diagrams = {}, {}
        for kind, incoming, outgoing, bare_count, ct_count in [
            ("fermion", [fermion], [fermion], 1, 2),
            ("scalar", [scalar], [scalar], 2, 2),
            ("vertex", [fermion], [scalar, fermion], 1, 1),
            ("quartic", [scalar] * 2, [scalar] * 2, 9, 1),
        ]:
            generated = stage_model.process(incoming, outgoing).generate_diagrams(
                loops=loops,
                max_vertices=len(incoming) + len(outgoing) if loops else 1,
                coupling_orders={"CT": 1} if not loops else None,
                **options,
            )
            diagrams[kind] = generated.diagrams
            assert len(generated.diagrams) == (bare_count if loops else ct_count), (
                variant,
                kind,
                len(generated.diagrams),
            )
            for diagram in generated.diagrams:
                ports = {}
                if kind in ("fermion", "vertex"):
                    for edge in diagram.external_edges:
                        if kind == "vertex" and edge.external_index == 1:
                            continue
                        ports[edge.external_index] = dict(
                            next(
                                diagram.projector_expression().match(
                                    wave(edge.id, bis(4, index)), max_level=0
                                )
                            )
                        )[index]
                numerator = stage_model.expand_couplings(
                    diagram.numerator_expression().to_expression()
                )
                if loops:
                    numerator = diagram.uv_expansion(
                        mUV, numerator=numerator
                    ).to_expression()
                numerator = diagram.momentum_basis().route_expression(numerator)
                numerator = (
                    numerator.replace(mink(dim, index), mink(D, index))
                    .replace(mink(dim), mink(D))
                    .replace(mUV**2, M)
                )
                for match in list(numerator.match(den(edge_, mom_, mass_, quad_))):
                    values = dict(match)
                    formal = family.rewrite_numerator(values[quad_], [coordinate])
                    assert formal == coordinate
                    numerator = numerator.replace(
                        den(values[edge_], values[mom_], values[mass_], values[quad_]),
                        formal,
                    )
                factor = (
                    diagram.overall_factor_expression(evaluate=True)
                    * diagram.numerator_prefactor_expression()
                )
                if not loops:
                    factor /= Symbol.I * h
                if kind == "fermion":
                    probes = [
                        (
                            "fermion_p",
                            gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
                            * P(0, mink(D, mu))
                            / (4 * s),
                        ),
                        (
                            "fermion_m",
                            metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass),
                        ),
                    ]
                elif kind == "vertex":
                    probe = (
                        metric(bis(4, ports[0]), bis(4, ports[2]))
                        if variant == "Scalar"
                        else TensorExpression.gamma5(4)(
                            ports[0], ports[2]
                        ).to_expression()
                    )
                    # Project both vertex phases onto the scalar convention;
                    # this removes the known i before the rational matrix solve.
                    vertex_phase = local_coupling / (-Symbol.I * g)
                    probes = [(kind, probe / (4 * vertex_phase))]
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
                    expression = family.rewrite_numerator(
                        kinematics.apply(reducer.reduce(trace)), [coordinate]
                    )
                    expression = (expression * factor).together().expand()
                    terms = parts.setdefault(label, [])
                    for monomial, coefficient in expression.coefficient_list(
                        coordinate
                    ):
                        power = -int(
                            (
                                monomial.derivative(coordinate) * coordinate / monomial
                            ).together()
                        )
                        assert monomial == coordinate ** (-power)
                        assert not coefficient.matches(K(S("yukawa::args___")))
                        terms.append(([power], coefficient))
        stage_parts[stage], stage_diagrams[stage] = parts, diagrams
    parts = stage_parts["bare"]
    all_diagrams[variant] = stage_diagrams
    ibp = hep.IBPFamily(family, name="yukawa_one_loop")
    solution = ibp.reduce_laporta(
        sorted({tuple(p) for terms in parts.values() for p, c in terms}), max_depth=2
    )
    assert solution.residuals == [[1]]
    reduced, poles = {}, {}
    for label, terms in parts.items():
        expression = sum(
            (c * solution.reduce(p, integral=integral) for p, c in terms), zero
        ).together()
        reduced[label] = expression
        poles[label] = (
            expression.replace(integral(1), master_pole / eps)
            .replace(D, 4 - 2 * eps)
            .series(eps, 0, -1)
            .to_expression()
            .expand()
        )
        assert poles[label].derivative(M).expand() == zero
    all_reductions[variant] = reduced
    all_poles[variant] = poles
    reference_poles = {
        "fermion_p": -(g**2) / (2 * eps),
        "fermion_m": (-1 if variant == "Scalar" else 1) * g**2 / eps,
        "scalar": (
            lam * scalar_mass**2 / 2
            + 2 * g**2 * s
            - (12 if variant == "Scalar" else 4) * g**2 * mass**2
        )
        / eps,
        "vertex": -(g**3) / eps,
        "quartic": (3 * lam**2 / 2 - 24 * g**4) / eps,
    }
    for label, reference in reference_poles.items():
        assert (poles[label] - reference).together() == zero
    # Counterterms and bare graphs retain the same native external-fermion
    # ordering. No independent phase adjustment is needed in the matching.
    counterterms = {}
    for label, terms in stage_parts["ct"].items():
        assert all(powers == [0] for powers, coefficient in terms)
        counterterms[label] = sum(
            (coefficient for powers, coefficient in terms), zero
        ).expand()
    ct_rows = [
        counterterms["fermion_p"],
        counterterms["fermion_m"],
        counterterms["scalar"].coefficient(s),
        counterterms["scalar"].replace(s, zero),
        counterterms["vertex"],
        counterterms["quartic"],
    ]
    loop_rows = [
        poles["fermion_p"],
        poles["fermion_m"],
        poles["scalar"].coefficient(s),
        poles["scalar"].replace(s, zero),
        poles["vertex"],
        poles["quartic"],
    ]
    entries = [row.coefficient(x) for row in ct_rows for x in unknowns]
    for row_index, row in enumerate(ct_rows):
        assert (
            row
            - sum(
                (entries[row_index * 6 + i] * x for i, x in enumerate(unknowns)), zero
            )
        ).expand() == zero
    matrix = Matrix.from_linear(6, 6, entries)
    solved = matrix.solve(Matrix.vec([-row for row in loop_rows]))
    constants = [solved[i, 0].to_expression().expand() for i in range(6)]
    expected = [
        -(g**2) / (2 * eps),
        (3 if variant == "Scalar" else -1) * g**2 / (2 * eps),
        -2 * g**2 / eps,
        (
            2 * g**2
            + lam / 2
            - (12 if variant == "Scalar" else 4) * g**2 * mass**2 / scalar_mass**2
        )
        / eps,
        5 * g**2 / (2 * eps),
        (4 * g**2 + 3 * lam / 2 - 24 * g**4 / lam) / eps,
    ]
    assert all(
        (actual - reference).together() == zero
        for actual, reference in zip(constants, expected, strict=True)
    ), (variant, constants)
    replacements = [
        Replacement(x, value) for x, value in zip(unknowns, constants, strict=True)
    ]
    for label, pole in poles.items():
        assert (
            pole + counterterms[label].replace_multiple(replacements)
        ).together() == zero
    all_counterterms[variant] = counterterms
    all_matrices[variant] = matrix
    all_constants[variant] = constants
    print(
        variant,
        "all six generated renormalization constants and pole cancellation passed",
        solution.stats,
        flush=True,
    )

# Check the finite scheme shift without supplying any renormalization constant.
# The additive quartic shift remains regular at lam=0, where Zlam alone does not.
log4pi, gamma_e = S("yukawa::log4pi", "yukawa::gamma_E")
scheme_constants, additive_shifts, beta_functions = {}, {}, {}
for variant, constants in all_constants.items():
    for scheme, delta in [("MS", 1 / eps), ("MSbar", 1 / eps + log4pi - gamma_e)]:
        shifts = [(constant * eps * delta).expand() for constant in constants]
        scheme_constants[variant, scheme] = [
            1 + shift / (16 * Symbol.PI**2) for shift in shifts
        ]
        replacements = [
            Replacement(x, value) for x, value in zip(unknowns, shifts, strict=True)
        ]
        for label, pole in all_poles[variant].items():
            remainder = (
                pole * (1 + eps * (log4pi - gamma_e))
                + all_counterterms[variant][label].replace_multiple(replacements)
            ).together()
            expected = zero if scheme == "MSbar" else pole * eps * (log4pi - gamma_e)
            assert (remainder - expected).together() == zero
    shifts = [(g * constants[4]).expand(), (lam * constants[5]).expand()]
    additive_shifts[variant] = shifts
    assert (shifts[1].replace(lam, zero) + 24 * g**4 / eps).together() == zero
    # Scale independence of mu^eps*g*Zg and mu^(2eps)*lam*Zlam.
    bare_couplings = [g + h * shifts[0], lam + h * shifts[1]]
    jacobian = Matrix.from_linear(
        2, 2, [value.derivative(x) for value in bare_couplings for x in (g, lam)]
    )
    flow = jacobian.solve(
        Matrix.vec([-eps * bare_couplings[0], -2 * eps * bare_couplings[1]])
    )
    beta = [
        flow[i, 0]
        .to_expression()
        .series(h, 0, 1)
        .to_expression()
        .series(eps, 0, 0)
        .to_expression()
        .expand()
        .coefficient(h)
        for i in range(2)
    ]
    assert (beta[0] - 5 * g**3).together() == zero
    assert (beta[1] - 3 * lam**2 - 8 * g**2 * lam + 48 * g**4).together() == zero
    beta_functions[variant] = [value / (16 * Symbol.PI**2) for value in beta]
assert beta_functions["Scalar"] == beta_functions["Pseudoscalar"]
print("Both schemes, the zero-quartic limit and derived beta functions passed")
