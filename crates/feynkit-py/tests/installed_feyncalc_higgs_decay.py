"""Generated massive, on-shell Higgs decays and their two-body widths.

References:
https://feyncalc.github.io/FeynCalcExamples/EW/Tree/H-FFbar
https://feyncalc.github.io/FeynCalcExamples/EW/Tree/H-WW
https://feyncalc.github.io/FeynCalcExamples/EW/Tree/H-ZZ

The vector channels describe a Higgs above the two-vector threshold, not the
off-shell WW* and ZZ* decays of the physical 125 GeV Higgs.
"""

from pathlib import Path

from symbolica import E, Expression, S
from symbolica.community import hep as fk
from symbolica.community.spenso import CookSettings, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
P = S("gammalooprs::P")
charge, sw, cw, mh, mw, mz, vev = S(
    "UFO::ee", "UFO::sw", "UFO::cw", "UFO::MH", "UFO::MW", "UFO::MZ", "UFO::vev"
)
conjugate, wrapped = S("spenso::conj", "higgs_decay::adjoint_index")
ports = S("higgs_decay::h", "higgs_decay::first", "higgs_decay::second")
wave, rep, index, left, right = S("wave_", "rep_", "index_", "left_", "right_")
metric = S("spenso::g")
pi = Expression.PI
alpha = S("higgs_decay::alpha")
physical_mass = S("higgs_decay::M", is_positive=True)
beta = S("higgs_decay::beta", is_positive=True)
mass_squared = physical_mass**2 * (1 - beta**2) / 4

# M>0 and 0<beta<1 select the physical, above-threshold two-body branch.
physical_kin = (
    fk.Kinematics()
    .with_scalar_product(P(0), P(0), physical_mass**2)
    .with_scalar_product(P(1), P(1), mass_squared)
    .with_scalar_product(P(2), P(2), mass_squared)
    .with_scalar_product(P(1), P(2), (physical_mass**2 - 2 * mass_squared) / 2)
)
measure = physical_kin.two_body_phase_space(P(1), P(2)).expand()
# Symbolica retains products under radicals; choose the positive root explicitly.
measure = measure.replace((physical_mass**4 * beta**2).sqrt(), physical_mass**2 * beta)
flux = physical_kin.flux(P(0))
assert (measure - beta / (32 * pi**2)).together() == E("0")
assert flux == 2 * physical_mass
assert model.particle_by_pdg(25).spin_sum(P(0), ports[0], wrapped(ports[0])) == E("1")

for pdgs, mass, yukawa, yukawa_mass in (
    ((11, -11), S("UFO::Me"), "ye", "yme"),
    ((4, -4), S("UFO::MC"), "yc", "ymc"),
    ((5, -5), S("UFO::MB"), "yb", "ymb"),
    ((-24, 24), mw, None, None),
    ((23, 23), mz, None, None),
):
    generated = fk.Process(model, [25], list(pdgs)).generate_diagrams(
        max_vertices=1, numerator_grouping=None, progress=None
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    assert diagram.denominator_expression(dimension=4).to_expression() == E("1")
    assert diagram.overall_factor_expression(evaluate=True) == E("1")
    numerator = model.expand_couplings(diagram.numerator_expression().to_expression())
    if yukawa is not None:
        numerator = numerator.replace(
            S(f"UFO::{yukawa}"), model.parameter(yukawa).expression
        )
        # UFO Yukawa masses and pole masses are distinct inputs. The gallery
        # comparison imposes their tree-level equality explicitly.
        numerator = numerator.replace(S(f"UFO::{yukawa_mass}"), mass)
    numerator = numerator.replace(vev, model.parameter("vev").expression)
    for edge in diagram.external_edges:
        if edge.external_index == 0:
            continue  # A scalar Higgs has no open spin or color index.
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
    for real in (charge, sw, cw, mh, mw, mz, mass):
        adjoint = adjoint.replace(conjugate(real), real)
    adjoint = TensorExpression(adjoint).wrap_indices(wrapped)

    color_projector = E("1")
    for position, particle_pdg in enumerate(pdgs, start=1):
        particle = model.particle_by_pdg(particle_pdg)
        if particle.color == 1:
            continue
        matches = []
        for slot in operator.structure.slots:
            original = slot.to_expression()
            if original.replace(ports[position], E("0")) == original:
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
        color_projector *= particle.color_sum(matches[0][left], matches[0][right])

    if yukawa is not None:
        # Dirac adjunction exchanges the open fermion-chain endpoints.
        spin_projector = model.particle_by_pdg(pdgs[0]).spin_sum(
            P(1), wrapped(ports[2]), ports[1]
        ) * model.particle_by_pdg(pdgs[1]).spin_sum(P(2), ports[2], wrapped(ports[1]))
    else:
        # Sum all three physical polarizations of each final vector, including
        # its longitudinal Proca term. There is no initial spin average.
        spin_projector = model.particle_by_pdg(pdgs[0]).spin_sum(
            P(1), ports[1], wrapped(ports[1])
        ) * model.particle_by_pdg(pdgs[1]).spin_sum(P(2), ports[2], wrapped(ports[2]))
    scalar = (
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
    assert scalar.is_scalar
    kin = (
        fk.Kinematics()
        .with_scalar_product(P(0), P(0), mh**2)
        .with_scalar_product(P(1), P(1), mass**2)
        .with_scalar_product(P(2), P(2), mass**2)
        .with_scalar_product(P(1), P(2), (mh**2 - 2 * mass**2) / 2)
    )
    squared = kin.apply(scalar.to_expression())
    if yukawa is not None:
        colors = abs(model.particle_by_pdg(pdgs[0]).color)
        expected = (
            colors * charge**2 * mass**2 * (mh**2 - 4 * mass**2) / (2 * mw**2 * sw**2)
        )
        assert squared.replace(mass, E("0")) == E("0")
    else:
        expected = (
            charge**2
            * (mh**4 - 4 * mh**2 * mass**2 + 12 * mass**4)
            / (4 * mw**2 * sw**2)
        )
        if pdgs == (23, 23):
            # This is the labeled amplitude: the identical-Z factor belongs to
            # phase space below, not to the amplitude or polarization sum.
            expected *= mw**4 / (mz**4 * cw**4)
    residual = (squared - expected).expand().replace(cw, (1 - sw**2).sqrt()).together()
    assert residual == E("0"), (pdgs, residual.format_plain())
    print(f"H -> {pdgs}: generated massive squared amplitude passed", flush=True)

    # Apply the explicit 1/2! only for the identical ZZ final state.
    symmetry = E("1/2") if pdgs == (23, 23) else E("1")
    width = symmetry * 4 * pi * measure * squared / flux
    width = width.replace(charge**2, 4 * pi * alpha)
    if yukawa is not None:
        expected_width = (
            colors * alpha * physical_mass * mass**2 * beta**3 / (8 * mw**2 * sw**2)
        )
    else:
        # Independent reference coefficients: ZZ has the identical-state 1/2.
        reference_denominator = 32 if pdgs == (23, 23) else 16
        expected_width = (
            alpha * physical_mass**3 * beta / (reference_denominator * mw**2 * sw**2)
        )
        expected_width *= (
            1
            - 4 * mass_squared / physical_mass**2
            + 12 * mass_squared**2 / physical_mass**4
        )
    if pdgs == (23, 23):
        # The gallery width uses the tree-level electroweak mass relation.
        width = width.replace(mw, cw * mz)
        expected_width = expected_width.replace(mw, cw * mz)
    width = width.replace(mh, physical_mass).replace(
        mass, physical_mass * (1 - beta**2).sqrt() / 2
    )
    expected_width = expected_width.replace(
        mass, physical_mass * (1 - beta**2).sqrt() / 2
    )
    residual = (
        (width - expected_width).expand().replace(cw, (1 - sw**2).sqrt()).together()
    )
    assert residual == E("0"), (pdgs, residual.format_plain())
    if yukawa is not None:
        assert width.replace(beta, E("1")).together() == E("0")
    print(
        f"H -> {pdgs}: gallery width and final-state normalization passed", flush=True
    )
