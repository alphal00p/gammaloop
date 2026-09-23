"""Full B -> eta_c topology input: algebraic partial fractions and completion.

Global topology minimization remains a separate check.
The input fixture records the pinned FeynCalc source and its digest.
"""

import json
from functools import reduce
from pathlib import Path

from symbolica import E, Replacement, S
from symbolica.community import hep

fixture = json.loads(
    (Path(__file__).parent / "fixtures/feyncalc_etac_topologies.json").read_text()
)
D, k1, k2, n, nb, gkin, meta, u0b = S(
    "etac::D",
    "etac::k1",
    "etac::k2",
    "etac::n",
    "etac::nb",
    "etac::gkin",
    "etac::meta",
    "etac::u0b",
    is_real=True,
)
kin = (
    hep.Kinematics(D, momenta=[k1, k2, n, nb])
    .with_scalar_product(n, n, E("0"))
    .with_scalar_product(nb, nb, E("0"))
    .with_scalar_product(n, nb, E("2"))
)
momenta = dict(zip(("k1", "k2", "n", "nb"), (k1, k2, n, nb)))
labels = S(*["etac::" + label for label in fixture["labels"]], is_real=True)
products = [
    kin.scalar_product(momenta[a], momenta[b]) for a, b in fixture["scalar_products"]
]
pool = []
for text in fixture["denominators"]:
    denominator = E(text)
    for label, product in zip(labels, products):
        denominator = denominator.replace(label, product)
    pool.append(denominator.expand())
assert len(pool) == len(set(pool)) == 89
assert len(fixture["integrals"]) == 251
assert len({tuple(sorted(row["propagators"])) for row in fixture["integrals"]}) == 248
# Two direct source checks fix the SFAD mass-term sign and the mixed linear form.
assert pool[0] == kin.scalar_product(k1, k1)
assert (
    pool[2]
    - kin.scalar_product(k1 + k2, k1 + k2)
    - 2 * gkin * meta * u0b * kin.scalar_product(k1 + k2, n)
    + meta * u0b * kin.scalar_product(k1 + k2, nb)
    + 2 * gkin * meta**2 * u0b**2
).expand() == E("0")

sectors = {}
input_rows = []
for number, row in enumerate(fixture["integrals"], 1):
    denominators = [pool[i] for i in row["propagators"]]
    family = hep.IntegralFamily([k1, k2], [n, nb], denominators, kinematics=kin)
    fractions = family.partial_fraction(row["powers"])
    original = reduce(
        lambda a, b: a * b,
        (d ** (-p) for d, p in zip(denominators, row["powers"])),
        E("1"),
    )
    reconstructed = E("0")
    for coefficient, powers in fractions:
        reconstructed += coefficient * reduce(
            lambda a, b: a * b,
            (d ** (-p) for d, p in zip(denominators, powers)),
            E("1"),
        )
        sector = family.sector(powers)
        assert sector.is_independent
        key = tuple(sorted(d.to_canonical_string() for d in sector.denominators))
        sectors.setdefault(key, sector)
    assert (original - reconstructed).together() == E("0"), number
    input_rows.append((family, fractions))
assert len(sectors) == 677

statistics = {
    "scaleless certificate": 0,
    "not detected": 0,
    "transverse certificate": 0,
}
completed_rows = []
for family in sectors.values():
    parameters = [S(f"etac::x{i}") for i in range(len(family.denominators))]
    U, F = family.symanzik(parameters)
    if U == E("0"):
        direction = family.scaleless_transverse_direction()
        assert direction is not None and any(w != E("0") for w in direction)
        assert all(w.is_real() is True for w in direction)
        # Independently shift each loop by w_i*r_perp. External products stay
        # fixed; z_i=r_perp.k_i and z2=r_perp^2 are independent formal symbols.
        loops = [k1, k2]
        z = S("etac_shift::z1", "etac_shift::z2")
        z_squared = S("etac_shift::squared")
        shifts = [
            Replacement(
                kin.scalar_product(ki, kj),
                kin.scalar_product(ki, kj)
                + direction[i] * z[j]
                + direction[j] * z[i]
                + direction[i] * direction[j] * z_squared,
            )
            for i, ki in enumerate(loops)
            for j, kj in enumerate(loops[i:], i)
        ]
        for denominator in family.denominators:
            assert (denominator.replace_multiple(shifts) - denominator).together() == E(
                "0"
            )
        status = "transverse certificate"
    else:
        weights = family.scaleless_scaling(parameters)
        status = "not detected" if weights is None else "scaleless certificate"
        if weights is not None:
            G = U + F
            assert (
                sum(
                    (w * x * G.derivative(x) for w, x in zip(weights, parameters)),
                    E("0"),
                )
                - G
            ).expand() == E("0")
    statistics[status] += 1
    if status == "not detected":
        # All surviving families must pass automatic self-mapping before any
        # global minimization can rely on the shift search. Fractional
        # light-cone offsets previously broke candidate reconstruction.
        mapping = family.find_mapping(family)
        assert mapping is not None
        for source, target in enumerate(mapping.denominator_map):
            assert (
                mapping.apply(family.denominators[source]) - family.denominators[target]
            ).together() == E("0")
    completed = family.complete(candidates=pool)
    assert completed.is_complete and completed.is_independent
    assert completed.denominators[: len(family.denominators)] == family.denominators
    assert all(d in pool for d in completed.denominators)
    assert completed.complete(candidates=pool).denominators == completed.denominators
    # Reconstruct every independent loop scalar product from the selected basis.
    names = [S(f"etac::d{i}") for i in range(len(completed.denominators))]
    rules = completed.scalar_product_rules(names)
    for scalar_product, expression in rules:
        for name, denominator in zip(names, completed.denominators):
            expression = expression.replace(name, denominator)
        assert (expression - scalar_product).together() == E("0")
    completed_rows.append((family, completed, status))
assert statistics == {
    "scaleless certificate": 112,
    "not detected": 434,
    "transverse certificate": 131,
}
print(
    "PASS: 251 exact partial fractions, 677 independent sectors, 112 parametric and 131 transverse certificates, 434 automatic self-mappings, and 677 preferred completions",
    flush=True,
)
print(
    "Pending: global topology minimization",
    flush=True,
)
