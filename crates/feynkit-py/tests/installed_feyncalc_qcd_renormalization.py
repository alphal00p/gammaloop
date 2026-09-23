"""Generated one-loop QCD UV poles and MS/MSbar counterterms in symbolic gauge.

Reference: https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/Renormalization
The shared graph expansion retains auxiliary mass corrections. Counterterm
operators and the analytic tadpole pole are explicit inputs. Bare-factor
expansion and the ordinary generator supply actual CT diagrams; their projected
numerators determine every entry of the matching matrix.
Reference MS/MSbar scope:
https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/Renormalization2
The separate massive and massless n=0 IRR calculations directly massify only
massless propagators and reproduce their distinct auxiliary gluon mass terms.
Massless reference:
https://feyncalc.github.io/FeynCalcExamples/QCD/OneLoop/RenormalizationMassless
"""

import json
from pathlib import Path

from symbolica import E, Matrix, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.hep import oneloop
from symbolica.community.spenso import (
    ColorCasimirSettings,
    Representation,
    TensorExpression,
)

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
D, epsilon, M, mUV, coordinate, integral, xi, s, Nf = S(
    "qcd_ren::D",
    "qcd_ren::eps",
    "qcd_ren::M",
    "qcd_ren::mUV",
    "qcd_ren::d0",
    "qcd_ren::I",
    "qcd_ren::xi",
    "qcd_ren::s",
    "qcd_ren::Nf",
)
Nc, dA, CF, CA = S("qcd_ren::Nc", "qcd_ren::dA", "qcd_ren::CF", "qcd_ren::CA")
specification = json.loads(model.to_json())
for propagator in specification["propagators"]:
    if propagator["particle"] == "g":
        propagator["numerator"] = (
            "-1𝑖*(UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))"
            "-(1-qcd_ren::xi)*UFO::P(UFO::idx(1,1))*UFO::P(UFO::idx(1,2))"
            "/spenso::dot(UFO::P(spenso::mink(4)),UFO::P(spenso::mink(4))))"
        )
model = hep.Model.from_json(json.dumps(specification))
K, P, mass, gs = S("gammalooprs::K", "gammalooprs::P", "UFO::MB", "UFO::G")
mink, bis, gamma, metric = S(
    "spenso::mink", "spenso::bis", "spenso::gamma", "spenso::g"
)
cof, coad, cas, trace_index = S(
    "spenso::cof", "spenso::coad", "spenso::cas", "spenso::idx"
)
index, dim, wave, mu, nu, left, right = S(
    "qcd_ren::index_",
    "qcd_ren::dim_",
    "qcd_ren::wave_",
    "qcd_ren::mu",
    "qcd_ren::nu",
    "qcd_ren::left_",
    "qcd_ren::right_",
)
den, edge_, mom_, mass_, quad_ = S(
    "gammalooprs::denom",
    "qcd_ren::edge_",
    "qcd_ren::mom_",
    "qcd_ren::mass_",
    "qcd_ren::quad_",
)
ordering, value, arguments = S(
    "feynkit_generator_factor::ExternalFermionOrderingSign",
    "qcd_ren::value_",
    "qcd_ren::arguments___",
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
quark, gluon, ghost = (model.particle_by_pdg(pdg) for pdg in (5, 21, 9000005))
interaction_rules = {}
for label, particles in [
    ("qqg", [quark.antiname, quark.name, gluon.name]),
    ("ccg", [ghost.antiname, ghost.name, gluon.name]),
    ("ggg", [gluon.name] * 3),
    ("gggg", [gluon.name] * 4),
]:
    interaction_rules[label] = [
        vertex
        for vertex in model.vertex_rules
        if sorted(vertex.particles) == sorted(particles)
    ]
    assert len(interaction_rules[label]) == 1
parts, diagrams, irr_inputs = {}, {}, {}
for kind, incoming, outgoing, loops, vertices, count in [
    ("tree", [quark], [gluon, quark], 0, interaction_rules["qqg"], 1),
    ("quark", [quark], [quark], 1, interaction_rules["qqg"], 1),
    ("ghost", [ghost], [ghost], 1, interaction_rules["ccg"], 1),
    ("gluon_loop", [gluon], [gluon], 1, interaction_rules["ggg"], 1),
    ("ghost_loop", [gluon], [gluon], 1, interaction_rules["ccg"], 1),
    ("quark_loop", [gluon], [gluon], 1, interaction_rules["qqg"], 1),
    ("tadpole", [gluon], [gluon], 1, interaction_rules["gggg"], 1),
    (
        "vertex",
        [quark],
        [gluon, quark],
        1,
        interaction_rules["qqg"] + interaction_rules["ggg"],
        2,
    ),
]:
    generated = model.generate_diagrams(
        incoming,
        outgoing,
        loops=loops,
        max_vertices=len(incoming) + len(outgoing) - 2 + 2 * loops,
        maximum_bridges=0,
        vertex_allow=vertices,
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(generated.diagrams) == count
    diagrams[kind] = generated.diagrams
    for number, diagram in enumerate(generated.diagrams):
        ports = {}
        if kind != "ghost":
            for edge in diagram.external_edges:
                is_gluon = incoming == [gluon] or (
                    kind in ("tree", "vertex") and edge.external_index == 1
                )
                rep = mink if is_gluon else bis
                ports[edge.external_index] = dict(
                    next(
                        diagram.projector_expression().match(
                            wave(edge.id, rep(4, index)), max_level=0
                        )
                    )
                )[index]
        numerator = model.expand_couplings(
            diagram.numerator_expression().to_expression()
        )
        if kind in ("tree", "vertex"):
            # Conjugate the tree's color tensor; the generated tree fixes the norm.
            color_projector = (
                TensorExpression.t(dA, Nc)(ports[1], ports[0], ports[2])
                .spenso_conjugate()
                .to_expression()
            )
            numerator = (
                numerator.replace(cof(3, index), cof(Nc, index)).replace(
                    coad(8, index), coad(dA, index)
                )
                * color_projector
            )
            settings = ColorCasimirSettings(rewrite_fundamental_dimension=False)
        else:
            particle = incoming[0]
            rep, numeric_dim, symbolic_dim = (
                (cof, 3, Nc) if kind == "quark" else (coad, 8, dA)
            )
            color_slots = [
                slot.dual().to_expression()
                for slot in TensorExpression(numerator).interface
                if slot.to_expression().matches(rep(numeric_dim, index))
            ]
            assert len(color_slots) == 2
            color_indices = dict(
                next(
                    metric(*color_slots).match(
                        particle.color_sum(left, right), max_level=0
                    )
                )
            )
            color_projector = particle.color_sum(
                color_indices[left], color_indices[right]
            )
            numerator = (
                (numerator * color_projector / symbolic_dim)
                .replace(cof(3, index), cof(Nc, index))
                .replace(coad(8, index), coad(dA, index))
            )
            settings = ColorCasimirSettings()
        numerator = (
            TensorExpression(numerator)
            .simplify_color()
            .to_color_casimir(
                fundamental=Representation.cof(Nc),
                adjoint=Representation.coad(dA),
                settings=settings,
            )
            .to_expression()
        )
        # CF and CA are display names for shared representation-aware invariants;
        # the quark-loop trace uses the conventional fundamental index TR=1/2.
        numerator = (
            numerator.replace(cas(2, cof(Nc)), CF)
            .replace(cas(2, coad(dA)), CA)
            .replace(trace_index(2, cof(Nc)), one / 2)
        )
        # Keep the unexpanded graph numerator for the independent direct IRR path.
        raw_numerator = numerator
        if loops:
            numerator = diagram.momentum_basis().route_expression(
                diagram.uv_expansion(mUV, numerator=numerator).to_expression()
            )
        numerator = (
            numerator.replace(mink(dim, index), mink(D, index))
            .replace(mink(dim), mink(D))
            .replace(mUV**2, M)
        )
        if incoming == [gluon]:
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
        if kind in ("quark", "ghost"):
            # Amputated Grassmann two-point kernels omit external ordering.
            # Closed-loop signs and all other graph factors remain included.
            raw = diagram.overall_factor_expression()
            removed = (raw / raw.replace(ordering(value), one)).replace(
                ordering(value), value
            )
            assert removed == -one
            factor /= removed
        if kind == "quark":
            probes = [
                (
                    "quark_p",
                    gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
                    * P(0, mink(D, mu))
                    / (4 * s),
                ),
                ("quark_m", metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass)),
            ]
        elif kind in ("tree", "vertex"):
            probes = [
                (
                    f"{kind}_{number}",
                    gamma(bis(4, ports[0]), bis(4, ports[2]), mink(D, ports[1]))
                    / (4 * D),
                )
            ]
        else:
            probes = [(kind, one)]
        if loops:
            irr_inputs[kind, number] = (
                diagram,
                raw_numerator,
                ports.copy(),
                probes,
                factor,
            )
        for label, projector in probes:
            traced = (
                TensorExpression((numerator * projector).expand())
                .simplify_gamma()
                .expand()
                .simplify_metrics()
                .to_dots()
                .to_expression()
            )
            scalar = family.rewrite_numerator(
                kinematics.apply(reducer.reduce(traced)), [coordinate]
            )
            scalar = (scalar * factor).together().expand()
            terms = []
            for monomial, coefficient in scalar.coefficient_list(coordinate):
                power = -int(
                    (monomial.derivative(coordinate) * coordinate / monomial).together()
                )
                assert monomial == coordinate ** (-power)
                assert not coefficient.matches(K(arguments))
                terms.append(([power], coefficient))
            parts[label] = terms

tree_terms = parts.pop("tree_0")
assert len(tree_terms) == 1 and tree_terms[0][0] == [0]
tree = tree_terms[0][1]
assert (tree + Symbol.I * gs * Nc * CF).together() == zero
solution = hep.IBPFamily(family, name="qcd_one_loop").reduce_laporta(
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
        expression * (Symbol.I / tree if label.startswith("vertex") else one) / gs**2
    )
    if label == "quark_loop":
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
uv_poles["gluon"] = sum(
    (
        uv_poles[label]
        for label in ("gluon_loop", "ghost_loop", "quark_loop", "tadpole")
    ),
    zero,
).expand()
uv_poles["vertex"] = (uv_poles["vertex_0"] + uv_poles["vertex_1"]).expand()
expected = [
    CF * xi / epsilon,
    -CF * (xi + 3) / epsilon,
    CA * (xi - 3) * s / (4 * epsilon),
    ((13 - 3 * xi) * CA - 4 * Nf) * (s * gmunu - ppmunu) / (6 * epsilon),
    (CF * xi + CA * (xi + 3) / 4) / epsilon,
]
for label, reference in zip(
    ("quark_p", "quark_m", "ghost", "gluon", "vertex"), expected, strict=True
):
    assert (uv_poles[label] - reference).together() == zero

# Local counterterm operators, with Z=1+a4*deltaZ, a4=gs²/(16*pi²).
# Unknowns: Zq, Zm, ZA, Zxi, Zc, Zg, ZAm, Zcm. The local operator basis
# is an explicit model input; all matching entries come from generated graphs.
qqg = interaction_rules["qqg"]
spec = json.loads(model.to_json())
vertex_definition = next(v for v in spec["vertex_rules"] if v["name"] == qqg[0].name)
# Auxiliary operators follow the reference model:
# M*(ZAm²-1)*A²/2 and M*(Zcm²-1)*cbar*c/2. Identical gluons supply
# the factor two in their local rule; the distinct ghost fields do not.
a4 = S("qcd_ct::a4")
unknowns = list(
    S(
        "qcd_ct::deltaZq",
        "qcd_ct::deltaZm",
        "qcd_ct::deltaZA",
        "qcd_ct::deltaZxi",
        "qcd_ct::deltaZc",
        "qcd_ct::deltaZg",
        "qcd_ct::deltaZAm",
        "qcd_ct::deltaZcm",
    )
)
Zq, Zm, ZA, Zxi, Zc, Zg, ZAm, Zcm = (one + a4 * delta for delta in unknowns)
spec["orders"].append({"name": "CT", "expansion_order": 1, "hierarchy": 1})
for label, particles, spins, color, lorentz, coupling in [
    (
        "qq_kinetic",
        [quark.antiname, quark.name],
        [2, 2],
        "Identity(1,2)",
        "Gamma(dummy(1),idx(1,1),idx(1,2))*P(dummy(1),2)",
        Symbol.I * (Zq - one),
    ),
    (
        "qq_mass",
        [quark.antiname, quark.name],
        [2, 2],
        "Identity(1,2)",
        "Identity(idx(1,1),idx(1,2))",
        -Symbol.I * mass * (Zq * Zm - one),
    ),
    (
        "gg_kinetic",
        [gluon.name] * 2,
        [3, 3],
        "Identity(1,2)",
        "Metric(idx(1,1),idx(1,2))*P(dummy(1),1)*P(dummy(1),1)-P(idx(1,1),1)*P(idx(1,2),1)",
        -Symbol.I * (ZA - one),
    ),
    (
        "gg_gauge",
        [gluon.name] * 2,
        [3, 3],
        "Identity(1,2)",
        "P(idx(1,1),1)*P(idx(1,2),1)",
        -Symbol.I * (ZA / Zxi - one) / xi,
    ),
    (
        "gg_auxmass",
        [gluon.name] * 2,
        [3, 3],
        "Identity(1,2)",
        "Metric(idx(1,1),idx(1,2))",
        Symbol.I * M * (ZAm**2 - one),
    ),
    (
        "ghost_kinetic",
        [ghost.antiname, ghost.name],
        [-1, -1],
        "Identity(1,2)",
        "P(dummy(1),2)*P(dummy(1),2)",
        Symbol.I * (Zc - one),
    ),
    (
        "ghost_auxmass",
        [ghost.antiname, ghost.name],
        [-1, -1],
        "Identity(1,2)",
        "1",
        Symbol.I * M * (Zcm**2 - one) / 2,
    ),
]:
    spec["lorentz_structures"].append(
        {"name": "CT_L_" + label, "spins": spins, "structure": lorentz}
    )
    spec["couplings"].append(
        {
            "name": "CT_GC_" + label,
            "expression": repr(coupling.series(a4, 0, 1).to_expression()),
            "orders": [["QCD", 2], ["CT", 1]],
            "value": None,
        }
    )
    spec["vertex_rules"].append(
        {
            "name": "CT_" + label,
            "particles": particles,
            "color_structures": [color],
            "lorentz_structures": ["CT_L_" + label],
            "couplings": [["CT_GC_" + label]],
        }
    )
# Copy every color/Lorentz slot of the actual model rule and dress its coupling.
ct_vertex = json.loads(json.dumps(vertex_definition))
ct_vertex["name"] = "CT_qqg"
for row, couplings in enumerate(ct_vertex["couplings"]):
    for col, coupling in enumerate(couplings):
        if coupling is None:
            continue
        dressed = model.expand_couplings(S("UFO::" + coupling)) * (
            Zq * Zg * ZA.sqrt() - one
        )
        name = f"CT_qqg_{row}_{col}"
        spec["couplings"].append(
            {
                "name": name,
                "expression": repr(dressed.series(a4, 0, 1).to_expression()),
                "orders": [["QCD", 3], ["CT", 1]],
                "value": None,
            }
        )
        couplings[col] = name
spec["vertex_rules"].append(ct_vertex)
ct_model = hep.Model.from_json(json.dumps(spec))
ct_quark, ct_gluon, ct_ghost = (
    ct_model.particle_by_pdg(pdg) for pdg in (5, 21, 9000005)
)
ct_vertices = [v for v in ct_model.vertex_rules if v.name.startswith("CT_")]
ct_diagrams, ct_coefficients = {}, {}
external_ordering = diagrams["tree"][0].overall_factor_expression(evaluate=True)
for kind, incoming, outgoing, count, qcd_order in [
    ("quark", [ct_quark], [ct_quark], 2, 2),
    ("gluon", [ct_gluon], [ct_gluon], 3, 2),
    ("ghost", [ct_ghost], [ct_ghost], 2, 2),
    ("vertex", [ct_quark], [ct_gluon, ct_quark], 1, 3),
]:
    # CT order is perturbative bookkeeping; two-point insertions have no loops.
    # Bound vertices explicitly to exclude arbitrarily long insertion chains.
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
        incoming, outgoing, coupling_orders={"QCD": qcd_order, "CT": 1}, **options
    )
    assert len(generated.diagrams) == count
    assert not ct_model.generate_diagrams(
        incoming, outgoing, coupling_orders={"CT": 0}, **options
    ).diagrams
    ct_diagrams[kind] = generated.diagrams
    for diagram in generated.diagrams:
        assert len(diagram.internal_edges) == 0
        assert diagram.symmetry_factor == one
        factor = (
            diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
        )
        assert factor == (
            external_ordering if kind in ("quark", "ghost", "vertex") else one
        )
        ports = {}
        if kind != "ghost":
            for edge in diagram.external_edges:
                rep = (
                    mink
                    if kind == "gluon"
                    or (kind == "vertex" and edge.external_index == 1)
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
        )
        if kind == "vertex":
            color_projector = (
                TensorExpression.t(dA, Nc)(ports[1], ports[0], ports[2])
                .spenso_conjugate()
                .to_expression()
            )
            numerator = (
                numerator.replace(cof(3, index), cof(Nc, index)).replace(
                    coad(8, index), coad(dA, index)
                )
                * color_projector
            )
            settings = ColorCasimirSettings(rewrite_fundamental_dimension=False)
        else:
            particle = incoming[0]
            rep, numeric_dim, symbolic_dim = (
                (cof, 3, Nc) if kind == "quark" else (coad, 8, dA)
            )
            color_slots = [
                slot.dual().to_expression()
                for slot in TensorExpression(numerator).interface
                if slot.to_expression().matches(rep(numeric_dim, index))
            ]
            assert len(color_slots) == 2, (kind, numerator)
            color_indices = dict(
                next(
                    metric(*color_slots).match(
                        particle.color_sum(left, right), max_level=0
                    )
                )
            )
            color_projector = particle.color_sum(
                color_indices[left], color_indices[right]
            )
            numerator = (
                (numerator * color_projector / symbolic_dim)
                .replace(cof(3, index), cof(Nc, index))
                .replace(coad(8, index), coad(dA, index))
            )
            settings = ColorCasimirSettings()
        numerator = (
            TensorExpression(numerator)
            .simplify_color()
            .to_color_casimir(
                fundamental=Representation.cof(Nc),
                adjoint=Representation.coad(dA),
                settings=settings,
            )
            .to_expression()
        )
        numerator = (
            numerator.replace(cas(2, cof(Nc)), CF)
            .replace(cas(2, coad(dA)), CA)
            .replace(trace_index(2, cof(Nc)), one / 2)
            .replace(mink(dim, index), mink(D, index))
            .replace(mink(dim), mink(D))
        )
        if kind == "quark":
            probes = [
                (
                    "quark_p",
                    gamma(bis(4, ports[0]), bis(4, ports[1]), mink(D, mu))
                    * P(0, mink(D, mu))
                    / (4 * s),
                ),
                ("quark_m", metric(bis(4, ports[0]), bis(4, ports[1])) / (4 * mass)),
            ]
            normalization = Symbol.I * a4 * external_ordering
        elif kind == "vertex":
            probes = [
                (
                    kind,
                    gamma(bis(4, ports[0]), bis(4, ports[2]), mink(D, ports[1]))
                    / (4 * D),
                )
            ]
            normalization = tree * a4
        else:
            if kind == "gluon":
                numerator = numerator.replace(mink(D, ports[0]), mink(D, mu)).replace(
                    mink(D, ports[1]), mink(D, nu)
                )
            probes = [(kind, one)]
            normalization = Symbol.I * a4
            if kind == "ghost":
                normalization *= external_ordering
        for label, projector in probes:
            trace = (
                TensorExpression((numerator * projector).expand())
                .simplify_gamma()
                .expand()
                .simplify_metrics()
                .to_dots()
                .to_expression()
            )
            coefficient = (
                (kinematics.apply(trace) * factor / normalization).together().expand()
            )
            ct_coefficients[label] = ct_coefficients.get(label, zero) + coefficient
ct_g = ct_coefficients["gluon"].coefficient(gmunu)
ct_pp = ct_coefficients["gluon"].coefficient(ppmunu)
ct_rows = [
    ct_coefficients["quark_p"],
    ct_coefficients["quark_m"],
    ct_g.coefficient(s),
    xi * ct_pp,
    ct_coefficients["ghost"].coefficient(s),
    ct_coefficients["vertex"],
    ct_g.replace(s, zero) / M,
    ct_coefficients["ghost"].replace(s, zero) / M,
]
entries = [r.derivative(z).together() for r in ct_rows for z in unknowns]
assert all(e.derivative(z) == zero for e in entries for z in unknowns)
for i, row in enumerate(ct_rows):
    assert (
        row - sum((entries[8 * i + j] * z for j, z in enumerate(unknowns)), zero)
    ).together() == zero
ct_matrix = Matrix.from_linear(8, 8, entries)
assert ct_matrix[6, 6].to_expression() == 2
assert ct_matrix[7, 7].to_expression() == 1
assert (ct_coefficients["gluon"] - ct_g * gmunu - ct_pp * ppmunu).together() == zero
assert ct_g.derivative(s).derivative(s) == zero
assert ct_coefficients["ghost"].derivative(s).derivative(s) == zero

gluon_g = uv_poles["gluon"].coefficient(gmunu)
gluon_pp = uv_poles["gluon"].coefficient(ppmunu)
rhs = Matrix.vec(
    [
        -uv_poles["quark_p"],
        -uv_poles["quark_m"],
        -gluon_g.coefficient(s),
        -xi * gluon_pp,
        -uv_poles["ghost"].coefficient(s),
        -uv_poles["vertex"],
        -gluon_g.replace(s, zero) / M,
        -uv_poles["ghost"].replace(s, zero) / M,
    ]
)
deltas = ct_matrix.solve(rhs)
references = [
    -CF * xi / epsilon,
    -3 * CF / epsilon,
    ((13 - 3 * xi) * CA - 4 * Nf) / (6 * epsilon),
    ((13 - 3 * xi) * CA - 4 * Nf) / (6 * epsilon),
    CA * (3 - xi) / (4 * epsilon),
    -(11 * CA - 2 * Nf) / (6 * epsilon),
    zero,
    zero,
]
for row, reference in enumerate(references):
    assert (deltas[row, 0].to_expression() - reference).together() == zero
residual = ct_matrix * deltas - rhs
assert all(residual[row, 0].to_expression().together() == zero for row in range(8))
assert deltas[1, 0].to_expression().derivative(xi).expand() == zero
assert deltas[5, 0].to_expression().derivative(xi).expand() == zero
print(
    "Generated QCD one-loop: symbolic gauge, four IBP targets, eight generated counterterms passed",
    solution.stats,
)

# The completed MS/MSbar reference only needs the local UV poles. Use the
# OneLOop measure conversion (4*pi)^eps*rGamma=1+cDelta*eps+O(eps²), where
# cDelta=log(4*pi)-EulerGamma, with D=4-2eps. Keep it symbolic in exact checks.
cDelta = S("qcd_ct::cDelta")
assert (-2 / (D - 4)).replace(D, 4 - 2 * epsilon) == 1 / epsilon
counterterms_by_scheme = {
    "MS": [deltas[row, 0].to_expression().together() for row in range(8)],
    "MSbar": [
        (epsilon * deltas[row, 0].to_expression() * (1 / epsilon + cDelta)).together()
        for row in range(8)
    ],
}
for scheme, constants in counterterms_by_scheme.items():
    for label, generated_ct in ct_coefficients.items():
        for delta, constant in zip(unknowns, constants, strict=True):
            generated_ct = generated_ct.replace(delta, constant)
        subtraction_part = uv_poles[label] * (
            one + epsilon * cDelta if scheme == "MSbar" else one
        )
        assert (generated_ct + subtraction_part).together() == zero
        # In one common conventional measure MS retains the finite cDelta*P.
        remainder = (
            uv_poles[label] * (one + epsilon * cDelta) + generated_ct
        ).together()
        expected_remainder = (
            epsilon * cDelta * uv_poles[label] if scheme == "MS" else zero
        )
        assert (remainder - expected_remainder).together() == zero
print("Generated QCD CT amplitudes: complete MS/MSbar tensor cancellation passed")


# Direct n=0 infrared rearrangement: massify only massless propagators, then
# Taylor-expand external momenta to the superficial divergence degree. This
# follows the original massive reference independently of graph.uv_expansion.
Q, dot, x, t, a, b = S(
    "gammalooprs::Q",
    "spenso::dot",
    "qcd_irr::x",
    "qcd_irr::t",
    "qcd_irr::a_",
    "qcd_irr::b_",
)
irr_kin = hep.Kinematics(D, momenta=[K(0), P(0), P(1)]).with_scalar_product(
    P(0), P(0), s
)
irr_results = {}
beyond_uv = S("qcd_irr::beyond_uv")
for scenario in ("massive", "massless"):
    massless = scenario == "massless"
    irr_scalars = {}
    # Match the reference's order of operations: set mq=0 before massification.
    # It is not the massless limit of an already rearranged massive integral.
    active_indices = [i for i in range(8) if not (massless and i == 1)]
    active_unknowns = [unknowns[i] for i in active_indices]
    for (kind, number), (diagram, raw, ports, probes, factor) in irr_inputs.items():
        if massless:
            raw = raw.replace(mass, zero)
        # Separate every primitive longitudinal 1/q² before changing propagators.
        # A common squared denominator would spuriously massify the Feynman term.
        tags = {
            edge.id: S(f"qcd_irr::gluon_{edge.id}")
            for edge in diagram.internal_edges
            if edge.particle_name == gluon.name
        }
        tagged = raw
        for edge_id, tag in tags.items():
            tagged = tagged.replace(dot(Q(edge_id, mink(4)), Q(edge_id, mink(4))), tag)
        components = (
            tagged.expand().coefficient_list(*tags.values())
            if tags
            else [(one, tagged)]
        )
        reconstructed = zero
        for monomial, numerator in components:
            powers = {}
            for edge_id, tag in tags.items():
                exponent = int((monomial.derivative(tag) * tag / monomial).together())
                assert exponent in (0, -1), (kind, exponent, monomial)
                powers[edge_id] = 1 - exponent
                monomial = monomial.replace(
                    tag, dot(Q(edge_id, mink(4)), Q(edge_id, mink(4)))
                )
            reconstructed += monomial * numerator
            if any(power == 2 for power in powers.values()):
                assert numerator.replace(xi, one).expand() == zero
            denominator = (
                diagram.denominator_expression(in_lmb=True, edge_powers=powers)
                .to_expression()
                .replace(mink(dim, index), mink(D, index))
                .replace(mink(dim), mink(D))
            )
            # Keep massive quark denominators intact; only massless ones get M.
            for match in list(denominator.match(den(edge_, mom_, mass_, quad_))):
                vals = dict(match)
                propagator_mass = (
                    vals[mass_].replace(mass, zero) if massless else vals[mass_]
                )
                quadratic = vals[quad_].replace(mass, zero) if massless else vals[quad_]
                denominator = denominator.replace(
                    den(vals[edge_], vals[mom_], vals[mass_], vals[quad_]),
                    quadratic - (M if propagator_mass == zero else zero),
                )
            numerator = (
                diagram.momentum_basis()
                .route_expression(numerator)
                .replace(mink(dim, index), mink(D, index))
                .replace(mink(dim), mink(D))
            )
            if kind in ("gluon_loop", "ghost_loop", "quark_loop", "tadpole"):
                numerator = numerator.replace(mink(D, ports[0]), mink(D, mu)).replace(
                    mink(D, ports[1]), mink(D, nu)
                )
            for label, projector in probes:
                if massless and label == "quark_m":
                    continue
                traced = (
                    TensorExpression((numerator * projector).expand())
                    .simplify_gamma()
                    .expand()
                    .simplify_metrics()
                    .to_dots()
                    .to_expression()
                )
                direct = irr_kin.apply(traced / denominator).replace(
                    vacuum.scalar_product(K(0), K(0)), x
                )
                direct = direct.replace_multiple(
                    [
                        Replacement(
                            dot(K(0, mink(D)), P(a, mink(D))),
                            t * dot(K(0, mink(D)), P(a, mink(D))),
                        ),
                        Replacement(
                            dot(P(a, mink(D)), K(0, mink(D))),
                            t * dot(P(a, mink(D)), K(0, mink(D))),
                        ),
                        Replacement(
                            dot(P(a, mink(D)), P(b, mink(D))),
                            t**2 * dot(P(a, mink(D)), P(b, mink(D))),
                        ),
                        Replacement(P(a, mink(D, index)), t * P(a, mink(D, index))),
                        Replacement(s, t**2 * s),
                    ]
                )
                # The pslash projector lowers the quark external degree by one.
                # Two-point gluon/ghost tensors need degree two; the vertex needs zero.
                degree = (
                    2
                    if kind
                    in ("gluon_loop", "ghost_loop", "quark_loop", "tadpole", "ghost")
                    else 0
                )
                taylor = (
                    direct.series(t, 0, degree + int(massless)).to_expression().expand()
                )
                direct = taylor.replace(t, one)
                if massless:
                    # The reference tags the first term beyond the divergence
                    # degree; prove that it contributes no UV pole after IBP.
                    direct += (beyond_uv - one) * taylor.coefficient(t ** (degree + 1))
                scalar = irr_kin.apply(reducer.reduce(direct)).replace(
                    vacuum.scalar_product(K(0), K(0)), x
                )
                scalar = (
                    (
                        scalar
                        * factor
                        / gs**2
                        * (Symbol.I / tree if kind == "vertex" else one)
                        * (Nf if kind == "quark_loop" else one)
                    )
                    .together()
                    .expand()
                )
                irr_scalars[label] = (irr_scalars.get(label, zero) + scalar).together()
        assert (reconstructed - raw).together() == zero
    irr_terms, irr_targets = {}, set()
    for label, scalar in irr_scalars.items():
        # Symbolica partial fractions separate the physical and auxiliary tadpoles.
        # Their principal parts must reconstruct the complete projected integrand.
        apart = scalar.apart(x)
        assert (apart - scalar).together() == zero
        reconstructed = zero
        terms = []
        for squared_mass in (M,) if massless else (M, mass**2):
            principal = (
                apart.replace(x, coordinate + squared_mass)
                .series(coordinate, 0, -1)
                .to_expression()
                .expand()
            )
            reconstructed += principal.replace(coordinate, x - squared_mass)
            for monomial, coefficient in principal.coefficient_list(coordinate):
                if coefficient == zero:
                    continue
                power = -int(
                    (monomial.derivative(coordinate) * coordinate / monomial).together()
                )
                assert monomial == coordinate ** (-power) and power > 0
                terms.append((squared_mass, power, coefficient))
                irr_targets.add((power,))
        assert (scalar - reconstructed).together() == zero, label
        irr_terms[label] = terms
    irr_solution = hep.IBPFamily(
        family, name="qcd_direct_irr_" + scenario
    ).reduce_laporta([list(p) for p in sorted(irr_targets)], max_depth=2)
    assert irr_solution.residuals == [[1]]
    irr_poles = {}
    for label, terms in irr_terms.items():
        expression = sum(
            (
                coefficient
                * irr_solution.reduce([power], integral=integral)
                .replace(M, squared_mass)
                .replace(integral(1), squared_mass / epsilon)
                for squared_mass, power, coefficient in terms
            ),
            zero,
        ).together()
        pole = (
            expression.replace(D, 4 - 2 * epsilon)
            .series(epsilon, 0, -1)
            .to_expression()
            .expand()
            .replace(mink(4, index), mink(D, index))
        )
        irr_poles[label] = pole.together().expand()
        assert pole.coefficient(epsilon**-2) == zero
        assert pole.derivative(mass).together() == zero
        assert pole.derivative(beyond_uv).together() == zero
    irr_poles["gluon"] = sum(
        (irr_poles[k] for k in ("gluon_loop", "ghost_loop", "quark_loop", "tadpole")),
        zero,
    ).expand()
    irr_poles["vertex"] = (irr_poles["vertex_0"] + irr_poles["vertex_1"]).expand()

    # Physical UV constants agree, but the auxiliary mass does depend on the
    # rearrangement prescription. The reference's quark loop retains its true mass.
    for label in (
        ["quark_p", "ghost", "vertex"]
        if massless
        else ["quark_p", "quark_m", "ghost", "vertex"]
    ):
        assert (irr_poles[label] - uv_poles[label]).together() == zero
    irr_auxiliary_pole = (
        (CA * (1 + 3 * xi) + (8 * Nf if massless else zero)) * M * gmunu / (4 * epsilon)
    )
    assert (
        irr_poles["gluon"] - uv_poles["gluon"] - irr_auxiliary_pole
    ).together() == zero
    if massless:
        assert (
            irr_poles["quark_loop"]
            - uv_poles["quark_loop"]
            - 2 * Nf * M * gmunu / epsilon
        ).together() == zero
    else:
        assert irr_poles["quark_loop"].derivative(M) == zero
    irr_gluon_g = irr_poles["gluon"].coefficient(gmunu)
    irr_gluon_pp = irr_poles["gluon"].coefficient(ppmunu)
    irr_rhs = Matrix.vec(
        [
            -irr_poles["quark_p"],
            -irr_poles.get("quark_m", zero),
            -irr_gluon_g.coefficient(s),
            -xi * irr_gluon_pp,
            -irr_poles["ghost"].coefficient(s),
            -irr_poles["vertex"],
            -irr_gluon_g.replace(s, zero) / M,
            -irr_poles["ghost"].replace(s, zero) / M,
        ]
    )
    active_matrix = Matrix.from_linear(
        len(active_indices),
        len(active_indices),
        [
            ct_matrix[i, j].to_expression()
            for i in active_indices
            for j in active_indices
        ],
    )
    active_rhs = Matrix.vec([irr_rhs[i, 0].to_expression() for i in active_indices])
    irr_deltas = active_matrix.solve(active_rhs)
    for row, original_row in enumerate(active_indices):
        reference = (
            -(CA * (1 + 3 * xi) + (8 * Nf if massless else zero)) / (8 * epsilon)
            if original_row == 6
            else deltas[original_row, 0].to_expression()
        )
        assert (irr_deltas[row, 0].to_expression() - reference).together() == zero
    irr_schemes = {}
    for scheme, measure in (("MS", one), ("MSbar", one + epsilon * cDelta)):
        constants = [
            (irr_deltas[row, 0].to_expression() * measure).together()
            for row in range(len(active_indices))
        ]
        irr_schemes[scheme] = constants
        irr_ct_rules = [
            Replacement(delta, constant)
            for delta, constant in zip(active_unknowns, constants, strict=True)
        ]
        for label, generated_ct in ct_coefficients.items():
            if massless and label == "quark_m":
                continue
            if massless:
                assert generated_ct.derivative(unknowns[1]) == zero
            actual_ct = generated_ct.replace_multiple(irr_ct_rules)
            assert (measure * irr_poles[label] + actual_ct).together() == zero
            common_measure_remainder = (
                (one + epsilon * cDelta) * irr_poles[label] + actual_ct
            ).together()
            expected_remainder = (
                epsilon * cDelta * irr_poles[label] if scheme == "MS" else zero
            )
            assert (common_measure_remainder - expected_remainder).together() == zero
    irr_results[scenario] = {
        "poles": irr_poles,
        "schemes": irr_schemes,
        "constants": irr_deltas,
        "indices": active_indices,
        "matrix": active_matrix,
        "solution": irr_solution,
        "residues": [
            (epsilon * irr_deltas[row, 0].to_expression()).together()
            for row in range(len(active_indices))
        ],
    }
    print(
        f"Direct {scenario} QCD IRR: {len(irr_targets)} IBP targets and {len(active_indices)} generated CT equations passed",
        irr_solution.stats,
    )

# Keep the massive result available for the notebook's original comparison.
irr_deltas = irr_results["massive"]["constants"]
irr_solution = irr_results["massive"]["solution"]
