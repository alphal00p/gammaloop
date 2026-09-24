"""Generated massive quark annihilation into gluons with physical polarizations.

Reference: https://feyncalc.github.io/FeynCalcExamples/QCD/Tree/QQbar-GlGl
The ordinary amplitudes include all three diagrams and their interference.
"""

from pathlib import Path

from symbolica import E, Replacement, S
from symbolica.community import hep as fk
from symbolica.community.spenso import CookSettings, TensorExpression

model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
P = S("gammalooprs::P")
s, t, u, mass, gs = S("s", "t", "u", "UFO::MB", "UFO::G")
Nc, dA, cof, coad = S("spenso::Nc", "dA", "spenso::cof", "spenso::coad")
kinematics = fk.Kinematics.mandelstam(
    [P(0), P(1), P(2), P(3)], [mass**2, mass**2, E("0"), E("0")], [s, t, u]
)
generated = fk.Generator(model).generate(
    fk.Process.amplitude([5, -5], [21, 21]).with_loop_count(0, 0),
    max_vertices=2,
    maximum_bridges=None,
    vertex_allow=["V_76", "V_36"],
    numerator_grouping=None,
    progress=None,
)
assert len(generated.diagrams) == 3
ports = S("external_0", "external_1", "external_2", "external_3")
a, b, c, inverse, wave, rep, index = S(
    "a_", "b_", "c_", "inverse_", "wave_", "rep_", "index_"
)
conjugate, adjoint_index = S("spenso::conj", "adjoint_index")
amplitude = E("0")
denominators = []
for diagram in generated.diagrams:
    numerator = model.expand_couplings(
        diagram.numerator_expression(in_lmb=True).to_expression()
    )
    # Align both spin and color ports using each generated external wavefunction.
    # This preserves interference without depending on internal half-edge numbering.
    for edge in diagram.external_edges:
        match = dict(
            next(
                diagram.projector_expression().match(
                    wave(edge.id, rep(4, index)), max_level=0
                )
            )
        )
        numerator = numerator.replace(match[index], ports[edge.external_index])
    denominator = kinematics.apply(
        diagram.denominator_expression(dimension=4, in_lmb=True)
        .to_expression()
        .replace(S("gammalooprs::denom")(a, b, c, inverse), inverse)
    ).expand()
    denominators.append(denominator)
    amplitude += (
        numerator
        * diagram.overall_factor_expression(evaluate=True)
        * diagram.numerator_prefactor_expression()
        / denominator
    )
assert all(d in denominators for d in (s, t - mass**2, u - mass**2))
operator = TensorExpression(amplitude.expand())
assert len(operator.structure.slots) == 8
adjoint = operator.dirac_adjoint().expand().simplify_gamma0().to_expression()
adjoint = adjoint.replace(conjugate(P(a, b)), P(a, b))
for real in (mass, gs, s, t, u):
    adjoint = adjoint.replace(conjugate(real), real)
adjoint = TensorExpression(adjoint).wrap_indices(adjoint_index)

# Match the shared completeness tensor against the dual of each external slot.
# Matching selects color slots and determines the quark/antiquark orientation.
left, right, metric = S("left_", "right_", "spenso::g")
color_projector = E("1")
initial_colors = 1
for position, pdg in enumerate((5, -5, 21, 21)):
    particle = model.particle_by_pdg(pdg)
    if position < 2:
        initial_colors *= abs(particle.color)
    matches = []
    for slot in operator.structure.slots:
        original = slot.to_expression()
        if original.replace(ports[position], E("0")) == original:
            continue
        closure = metric(
            slot.dual().to_expression(),
            original.replace(ports[position], adjoint_index(ports[position])),
        )
        match = next(closure.match(particle.color_sum(left, right), max_level=0), None)
        if match is not None:
            matches.append(match)
    assert len(matches) == 1
    indices = dict(matches[0])
    color_projector *= particle.color_sum(
        indices[left], indices[right], average=position < 2
    )

# Spenso dimensions are integers or symbols. Keep dA symbolic until color
# contraction, then impose the SU(N) relation and convert scalar Casimirs.
generic = (
    (operator.to_expression() * adjoint * color_projector * initial_colors / Nc**2)
    .replace(cof(3, index), cof(Nc, index))
    .replace(coad(8, index), coad(dA, index))
)
colored = (
    TensorExpression(generic, cook_indices=CookSettings.indices())
    .simplify_color()
    .to_expression()
    .replace(dA, Nc**2 - 1)
)
colored = TensorExpression(colored).to_cof_dimension_invariants()

# Dirac adjunction exchanges the two endpoints of the open fermion chain.
# Reduce the fermion trace before introducing the physical gluon projectors.
spin_projector = model.particle_by_pdg(5).spin_sum(
    P(0), ports[0], adjoint_index(ports[1]), average=True
) * model.particle_by_pdg(-5).spin_sum(
    P(1), adjoint_index(ports[0]), ports[1], average=True
)
spin_summed = (
    TensorExpression(
        colored.to_expression() * spin_projector, cook_indices=CookSettings.indices()
    )
    .expand()
    .simplify_gamma()
    .expand()
    .simplify_gamma()
)
expected = (
    (Nc**2 - 1)
    * gs**4
    * (
        mass**4 * (3 * t**2 + 14 * t * u + 3 * u**2)
        - mass**2 * (t**3 + 7 * t**2 * u + 7 * t * u**2 + u**3)
        - 6 * mass**8
        + t * u * (t**2 + u**2)
    )
    * (
        -2 * Nc**2 * mass**2 * (t + u)
        + 2 * Nc**2 * mass**4
        + Nc**2 * (t**2 + u**2)
        - s**2
    )
    / (2 * Nc**3 * s**2 * (u - mass**2) ** 2 * (t - mass**2) ** 2)
)
results = []
for references in ((P(3), P(2)), (P(0), P(0))):
    polarizations = E("1")
    for position, reference in zip((2, 3), references):
        polarizations *= model.particle_by_pdg(21).spin_sum(
            P(position),
            ports[position],
            adjoint_index(ports[position]),
            reference=reference,
        )
    scalar = (
        TensorExpression(
            spin_summed.to_expression() * polarizations,
            cook_indices=CookSettings.indices(),
        )
        .expand()
        .simplify_metrics()
        .to_dots()
    )
    assert scalar.is_scalar
    squared = (
        kinematics.apply(scalar.to_expression())
        .replace(s, 2 * mass**2 - t - u)
        .together()
    )
    assert (squared - expected.replace(s, 2 * mass**2 - t - u)).together() == E("0")
    results.append(squared)
    print(
        f"Generated massive SU(N) quark annihilation into gluons: references={references} passed"
    )
assert results[0] == results[1]

massless = results[0].replace(mass, E("0")).replace(Nc, E("3"))
expected_massless = (
    E("32/27") * gs**4 * (t**2 + u**2) / (t * u)
    - E("8/3") * gs**4 * (t**2 + u**2) / s**2
)
assert (massless - expected_massless.replace(s, -t - u)).together() == E("0")
# Bose exchange acts on the labeled squared amplitude. A cross section integrated
# over both identical-gluon labels would additionally require a factor of 1/2!.
assert (
    results[0] - results[0].replace_multiple([Replacement(t, u), Replacement(u, t)])
).together() == E("0")
print(
    "Massive gauge-reference independence, Bose symmetry and massless SU(3) reference passed"
)
