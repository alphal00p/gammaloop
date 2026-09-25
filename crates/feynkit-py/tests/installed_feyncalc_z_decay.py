"""Generated chiral Z decays, including physical polarization and decay widths.

Reference: https://feyncalc.github.io/FeynCalcExamples/EW/Tree/Z-FFbar
"""

from pathlib import Path

from symbolica import E, Expression, S
from symbolica.community import hep as fk
from symbolica.community.spenso import CookSettings, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
P = S("gammalooprs::P")
charge, sw, cw, mz = S("UFO::ee", "UFO::sw", "UFO::cw", "UFO::MZ")
conjugate, wrapped = S("spenso::conj", "z_decay::adjoint_index")
ports = S("z_decay::z", "z_decay::fermion", "z_decay::antifermion")
wave, rep, index, left, right = S("wave_", "rep_", "index_", "left_", "right_")
metric = S("spenso::g")
pi = Expression.PI

# Compare symbolic Dirac adjoints against the existing numerical Weyl matrices,
# independently multiplying gamma0 * matrix.conjugate().transpose() * gamma0.
matrices = (
    TensorExpression(E("spenso::gamma0(spenso::bis(4,i),spenso::bis(4,j))")),
    TensorExpression.gamma5(4)("i", "j"),
    TensorExpression.projm(4)("i", "j"),
    TensorExpression.projp(4)("i", "j"),
)
numeric_matrices = []
for matrix in matrices:
    network = matrix.to_network()
    network.execute()
    numeric_matrices.append([complex(value) for value in network.result_tensor()[:]])
gamma0 = numeric_matrices[0]
for matrix, values in zip(matrices, numeric_matrices, strict=True):
    adjoint = matrix.dirac_adjoint().simplify_gamma().undo_chain()
    network = adjoint.to_network()
    network.execute()
    actual = [complex(value) for value in network.result_tensor()[:]]
    expected_matrix = [
        sum(
            gamma0[4 * row + a] * values[4 * b + a].conjugate() * gamma0[4 * b + column]
            for a in range(4)
            for b in range(4)
        )
        for row in range(4)
        for column in range(4)
    ]
    assert actual == expected_matrix
print("Hermitian special matrices: all 64 Dirac-adjoint components passed", flush=True)

for pdg, mass, weak_isospin, electric_charge in (
    (12, E("0"), E("1/2"), E("0")),
    (11, S("UFO::Me"), E("-1/2"), E("-1")),
    (4, S("UFO::MC"), E("1/2"), E("2/3")),
    (5, S("UFO::MB"), E("-1/2"), E("-1/3")),
):
    particle = model.particle_by_pdg(pdg)
    assert particle.mass_expression == mass
    assert particle.weak_isospin == weak_isospin
    assert particle.charge == electric_charge
    generated = fk.Process(model, [23], [pdg, -pdg]).generate_diagrams(
        max_vertices=1, numerator_grouping=None, progress=None
    )
    assert len(generated.diagrams) == 1
    diagram = generated.diagrams[0]
    assert diagram.denominator_expression(dimension=4).to_expression() == E("1")
    numerator = model.expand_couplings(diagram.numerator_expression().to_expression())
    for edge in diagram.external_edges:
        matches = list(
            diagram.projector_expression().match(
                wave(edge.id, rep(4, index)), max_level=0
            )
        )
        assert len(matches) == 1
        # Both spin and color ports carry the external half-edge identifier.
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
    for real in (charge, sw, cw, mz, mass):
        adjoint = adjoint.replace(conjugate(real), real)
    adjoint = TensorExpression(adjoint).wrap_indices(wrapped)

    color_projector = E("1")
    for position, particle_pdg in ((1, pdg), (2, -pdg)):
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

    # The final fermions are summed, while the three initial Z polarizations
    # are averaged with the full Proca projector, including its longitudinal term.
    spin_projector = (
        model.particle_by_pdg(pdg).spin_sum(P(1), wrapped(ports[2]), ports[1])
        * model.particle_by_pdg(-pdg).spin_sum(P(2), ports[2], wrapped(ports[1]))
        * model.particle_by_pdg(23).spin_sum(
            P(0), ports[0], wrapped(ports[0]), average=True
        )
    )
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
        .with_scalar_product(P(0), P(0), mz**2)
        .with_scalar_product(P(1), P(1), mass**2)
        .with_scalar_product(P(2), P(2), mass**2)
        .with_scalar_product(P(1), P(2), (mz**2 - 2 * mass**2) / 2)
        .with_scalar_product(P(0), P(1), mz**2 / 2)
        .with_scalar_product(P(0), P(2), mz**2 / 2)
    )
    squared = kin.apply(scalar.to_expression())
    colors = abs(model.particle_by_pdg(pdg).color)
    cv, ca = weak_isospin - 2 * electric_charge * sw**2, weak_isospin
    expected = (
        colors
        * charge**2
        / (3 * sw**2 * cw**2)
        * (cv**2 * (mz**2 + 2 * mass**2) + ca**2 * (mz**2 - 4 * mass**2))
    )
    residual = (
        ((squared - expected) * sw**2 * cw**2 / charge**2)
        .expand()
        .replace(cw, (1 - sw**2).sqrt())
        .together()
    )
    assert residual == E("0"), (pdg, residual.format_plain())
    print(
        f"Z -> {pdg}, {-pdg}: generated massive chiral squared amplitude passed",
        flush=True,
    )

    # Express the physical branch with positive M and beta. The on-shell
    # relation m^2=M^2(1-beta^2)/4 describes 0<beta<=1 above threshold.
    physical_mass = S("z_decay::M", is_positive=True)
    beta = E("1") if mass == E("0") else S("z_decay::beta", is_positive=True)
    mass_squared = physical_mass**2 * (1 - beta**2) / 4
    physical_kin = (
        fk.Kinematics()
        .with_scalar_product(P(0), P(0), physical_mass**2)
        .with_scalar_product(P(1), P(1), mass_squared)
        .with_scalar_product(P(2), P(2), mass_squared)
        .with_scalar_product(P(1), P(2), (physical_mass**2 - 2 * mass_squared) / 2)
    )
    measure = physical_kin.two_body_phase_space(P(1), P(2)).expand()
    # Symbolica keeps products under radicals intact; select the positive
    # root explicitly using M>0 and beta>0.
    measure = measure.replace(
        (physical_mass**4 * beta**2).sqrt(), physical_mass**2 * beta
    )
    flux = physical_kin.flux(P(0))
    assert (measure - beta / (32 * pi**2)).together() == E("0")
    assert flux == 2 * physical_mass
    physical_squared = squared.replace(mz, physical_mass)
    if mass != E("0"):
        physical_squared = physical_squared.replace(
            mass, physical_mass * (1 - beta**2).sqrt() / 2
        )
    # f and fbar are distinct final particles, so no factorial is present.
    width = 4 * pi * measure * physical_squared / flux
    gf = S("z_decay::G_F")
    width = width.replace(
        charge**2, 4 * E("2").sqrt() * gf * physical_mass**2 * sw**2 * cw**2
    )
    expected_width = (
        colors
        * gf
        * physical_mass**3
        * beta
        * E("2").sqrt()
        / (12 * pi)
        * (
            cv**2 * (1 + 2 * mass_squared / physical_mass**2)
            + ca**2 * (1 - 4 * mass_squared / physical_mass**2)
        )
    )
    residual = (
        (width - expected_width).expand().replace(cw, (1 - sw**2).sqrt()).together()
    )
    assert residual == E("0"), (pdg, residual.format_plain())
    if mass != E("0"):
        massless_width = width.replace(beta, E("1"))
        expected_massless = (
            colors * gf * physical_mass**3 * (cv**2 + ca**2) * E("2").sqrt() / (12 * pi)
        )
        assert (massless_width - expected_massless).expand().replace(
            cw, (1 - sw**2).sqrt()
        ).together() == E("0")
    else:
        assert (
            width - gf * physical_mass**3 * E("2").sqrt() / (24 * pi)
        ).expand().replace(cw, (1 - sw**2).sqrt()).together() == E("0")

    if pdg == 11:
        # Isolate the vector and axial coefficients of the generated charged-lepton
        # current: at sw^2=1/4 its vector coupling vanishes, and the difference
        # from sw=0 removes the same axial contribution.
        normalized = (
            (squared * 3 * sw**2 * cw**2 / charge**2)
            .expand()
            .replace(cw, (1 - sw**2).sqrt())
            .expand()
        )
        axial = 4 * normalized.replace(sw, E("1/2"))
        vector = 4 * normalized.replace(sw, E("0")) - axial
        assert (axial - (mz**2 - 4 * mass**2)).together() == E("0")
        assert (vector - (mz**2 + 2 * mass**2)).together() == E("0")
        axial_threshold = (
            (beta * axial).replace(mz, physical_mass).replace(mass**2, mass_squared)
        )
        vector_threshold = (
            (beta * vector).replace(mz, physical_mass).replace(mass**2, mass_squared)
        )
        assert (axial_threshold - physical_mass**2 * beta**3).together() == E("0")
        assert (
            vector_threshold - physical_mass**2 * beta * (3 - beta**2) / 2
        ).together() == E("0")
    print(
        f"Z -> {pdg}, {-pdg}: gallery width and massless normalization passed",
        flush=True,
    )
