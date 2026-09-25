"""Generated W and top widths with masses, complex CKM and color averaging.

References: FeynCalc EW/Tree/W-ElAnel, W-QiQjbar and Qt-QbW.
The top reference sums initial colors; the physical width averages them.
"""

import math

from symbolica import E, Replacement, S, Symbol
from symbolica.community import hep
from symbolica.community.spenso import CookSettings, TensorExpression

model = hep.Model.standard_model()
P = S("gammalooprs::P")
charge, sw, mw = S("UFO::ee", "UFO::sw", "UFO::MW")
me, mc, mb, mt = S("UFO::Me", "UFO::MC", "UFO::MB", "UFO::MT")
conjugate, wrapped = S("spenso::conj", "weak_decay::adjoint_index")
ports = S("weak_decay::parent", "weak_decay::first", "weak_decay::second")
wave, rep, index, left, right = S(
    "weak_decay::wave_",
    "weak_decay::rep_",
    "weak_decay::index_",
    "weak_decay::left_",
    "weak_decay::right_",
)
metric = S("spenso::g")
Vr, Vi, GF = S("weak_decay::Vr", "weak_decay::Vi", "weak_decay::GF")
H, x, y = (S("weak_decay::" + name, is_positive=True) for name in ("H", "x", "y"))
zero, one, pi = E("0"), E("1"), Symbol.PI
results, diagrams = {}, {}
for label, pdgs, masses, ckm_name in [
    ("W- leptons", (-24, 11, -12), (mw, me, zero), None),
    ("W+ leptons", (24, -11, 12), (mw, me, zero), None),
    ("W- quarks", (-24, -4, 5), (mw, mc, mb), "CKM2x3"),
    ("W+ quarks", (24, 4, -5), (mw, mc, mb), "CKM2x3"),
    ("W- light quarks", (-24, -2, 1), (mw, zero, zero), "CKM1x1"),
    ("W+ light quarks", (24, 2, -1), (mw, zero, zero), "CKM1x1"),
    ("Top", (6, 5, 24), (mt, mb, mw), "CKM3x3"),
    ("Antitop", (-6, -5, -24), (mt, mb, mw), "CKM3x3"),
]:
    generated = hep.Process(model, [pdgs[0]], list(pdgs[1:])).generate_diagrams(
        max_vertices=1, numerator_grouping=None, progress=None
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    diagrams[label] = diagram
    assert diagram.denominator_expression(dimension=4).to_expression() == one
    particles = [model.particle_by_pdg(pdg) for pdg in pdgs]
    numerator = model.expand_couplings(diagram.numerator_expression().to_expression())
    mixing = one
    if ckm_name is not None:
        # Specify a generic complex CKM element without imposing the embedded
        # model's default parameter values or Wolfenstein approximation.
        ckm = S("UFO::" + ckm_name)
        numerator = numerator.replace(
            S("UFO::complexconjugate")(ckm), Vr - Symbol.I * Vi
        ).replace(ckm, Vr + Symbol.I * Vi)
        mixing = Vr**2 + Vi**2
    for edge in diagram.external_edges:
        matches = list(
            diagram.projector_expression().match(
                wave(edge.id, rep(4, index)), max_level=0
            )
        )
        assert len(matches) == 1
        numerator = numerator.replace(
            dict(matches[0])[index], ports[edge.external_index]
        )
    operator = TensorExpression(
        (
            numerator
            * diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
        ).expand()
    )
    adjoint = operator.dirac_adjoint().expand().simplify_gamma0().to_expression()
    for real in (charge, sw, mw, me, mc, mb, mt, Vr, Vi):
        adjoint = adjoint.replace(conjugate(real), real)
    adjoint = TensorExpression(adjoint).wrap_indices(wrapped)

    is_top = abs(pdgs[0]) == 6
    vector_position = 2 if is_top else 0
    # Dirac adjunction exchanges the endpoints of the open fermion chain.
    if is_top:
        left_position, right_position = (1, 0) if pdgs[0] > 0 else (0, 1)
    else:
        left_position = next(i for i in (1, 2) if pdgs[i] > 0)
        right_position = 3 - left_position
    spin_projector = (
        particles[left_position].spin_sum(
            P(left_position),
            wrapped(ports[right_position]),
            ports[left_position],
            average=left_position == 0,
        )
        * particles[right_position].spin_sum(
            P(right_position),
            ports[right_position],
            wrapped(ports[left_position]),
            average=right_position == 0,
        )
        * particles[vector_position].spin_sum(
            P(vector_position),
            ports[vector_position],
            wrapped(ports[vector_position]),
            average=vector_position == 0,
        )
    )
    parent, first, second = masses
    kin = (
        hep.Kinematics()
        .with_scalar_product(P(0), P(0), parent**2)
        .with_scalar_product(P(1), P(1), first**2)
        .with_scalar_product(P(2), P(2), second**2)
        .with_scalar_product(P(1), P(2), (parent**2 - first**2 - second**2) / 2)
        .with_scalar_product(P(0), P(1), (parent**2 + first**2 - second**2) / 2)
        .with_scalar_product(P(0), P(2), (parent**2 - first**2 + second**2) / 2)
    )
    physical_first = x if first != zero else zero
    physical_second = y if second != zero else zero
    physical_masses = [Replacement(parent, H)]
    if first != zero:
        physical_masses.append(Replacement(first, x))
    if second != zero:
        physical_masses.append(Replacement(second, y))
    phase_kin = (
        hep.Kinematics()
        .with_scalar_product(P(0), P(0), H**2)
        .with_scalar_product(P(1), P(1), physical_first**2)
        .with_scalar_product(P(2), P(2), physical_second**2)
        .with_scalar_product(
            P(1), P(2), (H**2 - physical_first**2 - physical_second**2) / 2
        )
    )
    measure = phase_kin.two_body_phase_space(P(1), P(2))
    radicand = S("weak_decay::radicand_")
    for match in list(measure.match(radicand.sqrt())):
        value = dict(match)[radicand]
        measure = measure.replace(value.sqrt(), value.expand().sqrt())
    flux = phase_kin.flux(P(0))
    assert flux == 2 * H
    kallen = (
        H**4
        + physical_first**4
        + physical_second**4
        - 2 * H**2 * physical_first**2
        - 2 * H**2 * physical_second**2
        - 2 * physical_first**2 * physical_second**2
    ).expand()
    # The positive root applies on the open decay domain H > m1+m2.
    reference_measure = kallen.sqrt() / (32 * pi**2 * H**2)
    assert (measure - reference_measure).together() == zero

    for convention in ("Physical", "Gallery color sum"):
        color_projector = one
        for position, particle in enumerate(particles):
            if particle.color == 1:
                continue
            matches = []
            for slot in operator.structure.slots:
                original = slot.to_expression()
                if original.replace(ports[position], zero) == original:
                    continue
                closure = metric(
                    slot.dual().to_expression(),
                    original.replace(ports[position], wrapped(ports[position])),
                )
                match = next(
                    closure.match(particle.color_sum(left, right), max_level=0), None
                )
                if match is not None:
                    matches.append(dict(match))
            assert len(matches) == 1
            color_projector *= particle.color_sum(
                matches[0][left],
                matches[0][right],
                average=position == 0 and convention == "Physical",
            )
        contracted = (
            TensorExpression(
                operator.to_expression() * adjoint * spin_projector * color_projector,
                cook_indices=CookSettings.indices(),
            )
            .expand()
            .simplify_gamma()
            .expand()
            .simplify_gamma()
            .simplify_color()
            .simplify_metrics()
            .to_dots()
        )
        assert contracted.is_scalar
        squared = (
            kin.apply(contracted.to_expression())
            .replace(charge**2, 4 * E("2").sqrt() * GF * mw**2 * sw**2)
            .replace_multiple(physical_masses)
            .together()
        )
        color_factor = (
            (3 if convention == "Gallery color sum" else 1)
            if is_top
            else (3 if ckm_name is not None else 1)
        )
        if is_top:
            polynomial = (H**2 - x**2) ** 2 + y**2 * (H**2 + x**2) - 2 * y**4
            expected_squared = E("2").sqrt() * color_factor * GF * mixing * polynomial
        else:
            a, b = physical_first**2, physical_second**2
            polynomial = 2 * H**4 - H**2 * (a + b) - (a - b) ** 2
            expected_squared = (
                2 * E("2").sqrt() * color_factor * GF * mixing * polynomial / 3
            )
        assert (squared - expected_squared).together() == zero, (
            label,
            convention,
            squared,
        )
        width = (4 * pi * measure * squared / flux).together()
        expected_width = (
            color_factor
            * GF
            * mixing
            * kallen.sqrt()
            * polynomial
            * E("2").sqrt()
            / (pi * H**3 * (16 if is_top else 24))
        )
        assert (width - expected_width).together() == zero
        results[label, convention] = {
            "squared": squared,
            "width": width,
            "mixing": mixing,
            "first": physical_first,
            "second": physical_second,
        }
    ratio = 3 if is_top else 1
    assert (
        results[label, "Gallery color sum"]["width"]
        - ratio * results[label, "Physical"]["width"]
    ).together() == zero
    print(
        label, "massive squared amplitude and both color conventions passed", flush=True
    )

for negative, positive in [
    ("W- leptons", "W+ leptons"),
    ("W- quarks", "W+ quarks"),
    ("W- light quarks", "W+ light quarks"),
    ("Top", "Antitop"),
]:
    for convention in ("Physical", "Gallery color sum"):
        assert (
            results[negative, convention]["width"]
            - results[positive, convention]["width"]
        ).together() == zero

# Compare physical widths numerically with the dimensionless gallery formulas,
# including a non-real CKM element and points approaching the two-body threshold.
numeric_checks = []
for (label, convention), result in results.items():
    is_top = label in ("Top", "Antitop")
    ratios = (
        [(0.0, 0.46), (0.03, 0.46), (0.3, 0.69), (0.3, 0.7)]
        if is_top
        else [(0.0, 0.0), (0.01, 0.02), (0.2, 0.3), (0.45, 0.54)]
    )
    for r1, r2 in ratios:
        r1 = r1 if result["first"] != zero else 0.0
        r2 = r2 if result["second"] != zero else 0.0
        substitutions = [
            Replacement(H, E("100")),
            Replacement(x, E(str(100 * r1))),
            Replacement(y, E(str(100 * r2))),
            Replacement(GF, E("0.00001")),
            Replacement(Vr, E("0.8")),
            Replacement(Vi, E("0.3")),
        ]
        actual = complex(result["width"].replace_multiple(substitutions).evaluate({}))
        phase = math.sqrt(max(0.0, (1 - (r1 + r2) ** 2) * (1 - (r1 - r2) ** 2)))
        mixing = 1.0 if "leptons" in label else 0.8**2 + 0.3**2
        colors = (
            (3 if convention == "Gallery color sum" else 1)
            if is_top
            else (1 if "leptons" in label else 3)
        )
        shape = (
            (1 - r1**2) ** 2 + r2**2 * (1 + r1**2) - 2 * r2**4
            if is_top
            else 1 - (r1**2 + r2**2) / 2 - (r1**2 - r2**2) ** 2 / 2
        )
        expected = (
            colors
            * 1e-5
            * 100**3
            * mixing
            * phase
            * shape
            / ((8 if is_top else 6) * math.sqrt(2) * math.pi)
        )
        assert abs(actual - expected) < 2e-9, (
            label,
            convention,
            r1,
            r2,
            actual,
            expected,
        )
        numeric_checks.append((label, convention, r1, r2, actual, expected))
print("Charge conjugation, complex CKM and 64 finite-width comparisons passed")
