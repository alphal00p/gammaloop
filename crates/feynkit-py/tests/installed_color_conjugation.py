"""Compact color words use the same conjugation as explicit indexed networks."""

from symbolica import E
from symbolica.community.spenso import Representation, TensorExpression

for colors in (2, 3, 5):
    generator = TensorExpression.t(colors**2 - 1, colors)
    fundamental = Representation.cof(colors)
    explicit = generator("a", "i", "k") * generator("b", "k", "j")
    word = explicit.collect_chains(fundamental)
    conjugate = word.spenso_conjugate()
    expected = (generator("b", "j", "k") * generator("a", "k", "i")).collect_chains(
        fundamental
    )
    assert conjugate.to_expression() == expected.to_expression()
    assert conjugate.spenso_conjugate().to_expression() == word.to_expression()
    assert word.dirac_adjoint().to_expression() == conjugate.to_expression()
    assert conjugate.simplify_color().to_expression() == (
        explicit.spenso_conjugate().simplify_color().to_expression()
    )
    norm = (
        (word * conjugate)
        .simplify_color()
        .to_cof_dimension_invariants()
        .simplify_metrics()
    )
    assert norm.is_scalar
    assert norm.to_expression() == E(str((colors**2 - 1) ** 2)) / (4 * colors)

    loop = (
        generator("a", "i", "j") * generator("b", "j", "k") * generator("c", "k", "i")
    )
    trace = loop.collect_chains(fundamental)
    conjugate_trace = trace.spenso_conjugate()
    assert conjugate_trace.spenso_conjugate().to_expression() == trace.to_expression()
    assert conjugate_trace.simplify_color().to_expression() == (
        loop.spenso_conjugate().simplify_color().to_expression()
    )
    print(f"SU({colors}): compact/explicit color conjugation and norm passed")
