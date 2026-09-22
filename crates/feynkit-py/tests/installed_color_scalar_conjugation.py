"""Scalar weights in compact color words reuse Symbolica conjugation."""

from symbolica import E, S
from symbolica.community.spenso import Representation, TensorExpression, _, chain, trace

z, w, form_factor = S("color_weight_z", "color_weight_w", "color_form_factor")
coefficients = (
    z**2,
    z**-1,
    (z + w) ** 2,
    z ** E("1/2"),
    z**w,
    form_factor(z),
    1 / form_factor(z),
    TensorExpression(form_factor(z)).spenso_conjugate().to_expression(),
)
for colors in (2, 3, 5):
    fundamental = Representation.cof(colors)
    generator = TensorExpression.t(colors**2 - 1, colors)
    for coefficient in coefficients:
        conjugate_coefficient = TensorExpression(coefficient).spenso_conjugate()
        word = chain(
            fundamental("i"),
            fundamental.dual()("j"),
            coefficient * generator("a", _, _),
            generator("b", _, _),
        )
        expected = chain(
            fundamental("j"),
            fundamental.dual()("i"),
            generator("b", _, _),
            conjugate_coefficient * generator("a", _, _),
        )
        conjugate = word.spenso_conjugate()
        assert conjugate.to_expression() == expected.to_expression()
        assert conjugate.spenso_conjugate().to_expression() == word.to_expression()
        assert word.dirac_adjoint().to_expression() == conjugate.to_expression()
        color_trace = trace(
            fundamental,
            coefficient * generator("a", _, _),
            generator("b", _, _),
            generator("c", _, _),
        )
        expected_trace = trace(
            fundamental,
            generator("c", _, _),
            generator("b", _, _),
            conjugate_coefficient * generator("a", _, _),
        )
        assert (
            color_trace.spenso_conjugate().to_expression()
            == expected_trace.to_expression()
        )
        assert (
            color_trace.spenso_conjugate().spenso_conjugate().to_expression()
            == color_trace.to_expression()
        )
    print(f"SU({colors}): weighted compact chain/trace conjugation and adjoint passed")

fundamental = Representation.cof(3)
generator = TensorExpression.t(8, 3)
weighted_word = chain(
    fundamental("i"),
    fundamental.dual()("j"),
    z**2 * generator("a", _, _),
    generator("b", _, _),
)
weighted_norm = TensorExpression(
    weighted_word.to_expression() * weighted_word.spenso_conjugate().to_expression()
)
# At z = 1 + 2i, |z^2|^2 = 25 and the unweighted SU(3) norm is 16/3.
network = (
    TensorExpression(weighted_norm.to_expression().replace(z, E("1+2𝑖")))
    .undo_chain()
    .to_network()
)
network.execute()
assert abs(complex(network.result_scalar()) - 400 / 3) < 1e-12
print("Independent weighted SU(3) matrix norm is 400/3")
