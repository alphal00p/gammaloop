"""Physical state counts, transverse projectors and sewing in arbitrary dimensions."""

from symbolica import E, S
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

model = hep.Model.standard_model()
p, n, i, j, D = S(
    "dim_spin::p", "dim_spin::n", "dim_spin::i", "dim_spin::j", "dim_spin::D"
)
mink, metric, ket, bra = S(
    "spenso::mink", "spenso::g", "gammalooprs::ϵ", "gammalooprs::ϵbar"
)
zero, one = E("0"), E("1")
for dimension in (E("4"), E("6"), D):
    for name, mass_squared, missing, reference in (
        ("g", zero, 2, n),
        ("Z", S("UFO::MZ") ** 2, 1, None),
    ):
        particle = model.particle(name)
        kin = hep.Kinematics(dimension).with_scalar_product(p, p, mass_squared)
        for average in (False, True):
            projector = particle.spin_sum(
                p, i, j, reference=reference, average=average, dimension=dimension
            )
            for vector in (p, n) if reference is not None else (p,):
                contraction = (
                    TensorExpression((projector * vector(mink(dimension, i))).expand())
                    .simplify_metrics()
                    .to_dots()
                )
                assert kin.apply(contraction.to_expression()).together() == zero
            trace = (
                TensorExpression(
                    (
                        projector * metric(mink(dimension, i), mink(dimension, j))
                    ).expand()
                )
                .simplify_metrics()
                .to_dots()
            )
            expected = -one if average else missing - dimension
            assert (kin.apply(trace.to_expression()) - expected).together() == zero
            pair = ket(7, mink(dimension, i)) * bra(7, mink(dimension, j))
            other = ket(8, mink(dimension, i))
            summed = particle.sum_spins(
                pair + other,
                p,
                edge=7,
                reference=reference,
                average=average,
                dimension=dimension,
            )
            assert (summed - projector - other).expand() == zero
            if dimension == 4:
                assert projector == particle.spin_sum(
                    p, i, j, reference=reference, average=average
                )
            # Fixed two-state averages are available independently of the Lorentz dimension.
            if name == "g" and average:
                fixed_two = (
                    particle.spin_sum(p, i, j, reference=n, dimension=dimension) / 2
                )
                assert (fixed_two - projector * (dimension - 2) / 2).together() == zero

for particle_name, invalid_dimension, options in (
    ("g", E("2"), {}),
    ("Z", E("1"), {}),
    ("e-", E("0"), {}),
    ("g", E("-1"), {}),
    ("g", E("3/2"), {}),
    ("g", D - 2, {}),
    ("ta-", D, {"spin_vector": n}),
):
    particle = model.particle(particle_name)
    for method, arguments, keywords in (
        (particle.spin_sum, (p, i, j), options),
        (particle.sum_spins, (one, p), {"edge": 7, **options}),
    ):
        try:
            method(*arguments, dimension=invalid_dimension, **keywords)
        except ValueError:
            pass
        else:
            raise AssertionError((particle_name, invalid_dimension, method))
print(
    "Dimensional spin sums: physical counts, transversality, sewing, defaults and invalid-state checks passed"
)
