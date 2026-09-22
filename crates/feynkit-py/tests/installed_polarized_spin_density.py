"""Massive physical-spin density matrices through shared public particle methods.

The spin vector is the boosted rest-frame spin direction for either charge.
Summing its two orientations must recover ordinary completeness; no extra
initial-state spin average belongs to a selected state.
"""

from pathlib import Path

from symbolica import E, S
from symbolica.community import hep as fk
from symbolica.community.spenso import GammaSimplifySettings, TensorExpression

# The fixture declares a nonzero tau mass; its muon mass parameter is zero.
model = fk.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
p, spin, mass = S("spin_density::p", "spin_density::s", "UFO::MTA")
i, j, k, slot = S(
    "spin_density::i", "spin_density::j", "spin_density::k", "spin_density::slot_"
)
bis = S("spenso::bis")
kinematics = (
    fk.Kinematics()
    .with_scalar_product(p, p, mass**2)
    .with_scalar_product(p, spin, E("0"))
    .with_scalar_product(spin, spin, E("-1"))
)

for pdg, sign, ket_name, bra_name in (
    (15, 1, "gammalooprs::u", "gammalooprs::ubar"),
    (-15, -1, "gammalooprs::v", "gammalooprs::vbar"),
):
    particle = model.particle_by_pdg(pdg)
    density = particle.spin_sum(p, i, j, spin_vector=spin)
    ordinary = particle.spin_sum(p, i, j)
    opposite = density.replace(spin(slot), -spin(slot))
    assert TensorExpression((density + opposite - ordinary).expand()).simplify_gamma(
        GammaSimplifySettings.canonical()
    ).expand().to_expression() == E("0")
    assert particle.spin_sum(p, i, j, spin_vector=None) == ordinary
    traced = (
        TensorExpression(density.replace(j, i))
        .simplify_gamma(GammaSimplifySettings.canonical())
        .to_expression()
    )
    assert (traced - sign * 2 * mass).expand() == E("0")

    # A fixed-spin outer product is rank one: rho^2 = (ubar u or vbar v) rho.
    # Both left and right Dirac equations also hold on the supplied mass shell.
    left_density = particle.spin_sum(p, i, k, spin_vector=spin)
    right_density = particle.spin_sum(p, k, j, spin_vector=spin)
    opposite_particle = model.particle_by_pdg(-pdg)
    left_dirac = opposite_particle.spin_sum(p, i, k)
    right_dirac = opposite_particle.spin_sum(p, k, j)
    for identity in (
        left_density * right_density - sign * 2 * mass * density,
        left_dirac * right_density,
        left_density * right_dirac,
    ):
        reduced = (
            TensorExpression(identity.expand())
            .simplify_gamma(GammaSimplifySettings.canonical())
            .expand()
            .to_dots()
            .to_expression()
        )
        assert kinematics.apply(reduced).expand() == E("0")

    ket, bra = S(ket_name, bra_name)
    pair = ket(7, bis(4, i)) * bra(7, bis(4, j))
    spectator = ket(8, bis(4, i)) * bra(8, bis(4, j))
    unpaired = ket(7, bis(4, i))
    # Test the substitution itself without multiplying separate terms that
    # reuse the same explicit spinor indices into a tensor network.
    output = particle.sum_spins(pair, p, edge=7, spin_vector=spin)
    assert TensorExpression((output - density).expand()).simplify_gamma(
        GammaSimplifySettings.canonical()
    ).expand().to_expression() == E("0")
    assert particle.sum_spins(output, p, edge=7, spin_vector=spin) == output
    assert particle.sum_spins(spectator, p, edge=7, spin_vector=spin) == spectator
    assert particle.sum_spins(unpaired, p, edge=7, spin_vector=spin) == unpaired
    other_ket = S("gammalooprs::v" if sign > 0 else "gammalooprs::u")
    wrong_species = other_ket(7, bis(4, i)) * bra(7, bis(4, j))
    assert (
        particle.sum_spins(wrong_species, p, edge=7, spin_vector=spin) == wrong_species
    )
    print(
        f"PDG {pdg}: completeness, trace, purity, Dirac equations and pair rules passed"
    )

# Both public paths validate the requested state even if there is no pair to
# replace. A symbolic vector name is required; its normalization is caller data.
for pdg, average, invalid_spin in (
    (15, True, spin),
    (-15, True, spin),
    (13, False, spin),
    (-13, False, spin),
    (12, False, spin),
    (-12, False, spin),
    (22, False, spin),
    (23, False, spin),
    (25, False, spin),
    (15, False, E("1")),
    (15, False, E("0")),
    (15, False, spin + p),
):
    particle = model.particle_by_pdg(pdg)
    for method, arguments, options in (
        (particle.spin_sum, (p, i, j), {}),
        (particle.sum_spins, (E("1"), p), {"edge": 7}),
    ):
        try:
            method(*arguments, **options, average=average, spin_vector=invalid_spin)
        except ValueError:
            pass
        else:
            raise AssertionError(
                f"Accepted invalid spin state: PDG {pdg}, average={average}, {invalid_spin}"
            )
print("Invalid polarized-state requests are rejected by both public methods")
