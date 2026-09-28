"""Compact color words use the same conjugation as explicit indexed networks."""

from symbolica import E, S
from symbolica.community.spenso import (
    AUTO,
    Representation,
    TensorExpression,
    chain,
)
from symbolica.community.spenso import (
    trace as tensor_trace,
)

for colors in (2, 3, 5):
    generator = TensorExpression.t(colors**2 - 1, colors)
    fundamental = Representation.cof(colors)
    explicit = generator("a", "i", "k") * generator("b", "k", "j")
    word = chain(
        fundamental("i"),
        fundamental.dual()("j"),
        generator("a", AUTO, AUTO),
        generator("b", AUTO, AUTO),
    )
    conjugate = word.spenso_conjugate()
    expected = chain(
        fundamental("j"),
        fundamental.dual()("i"),
        generator("b", AUTO, AUTO),
        generator("a", AUTO, AUTO),
    )
    assert conjugate.to_expression() == expected.to_expression()
    assert conjugate.spenso_conjugate().to_expression() == word.to_expression()
    assert word.dirac_adjoint().to_expression() == conjugate.to_expression()
    assert conjugate.simplify_color().to_expression().to_expression() == (
        explicit.spenso_conjugate().simplify_color().to_expression().to_expression()
    )
    norm = (
        (word * conjugate)
        .simplify_color()
        .to_expression()
        .to_cof_dimension_invariants()
        .contract()
        .to_expression()
    )
    assert norm.is_scalar
    assert norm.to_expression() == E(str((colors**2 - 1) ** 2)) / (4 * colors)

    loop = (
        generator("a", "i", "j") * generator("b", "j", "k") * generator("c", "k", "i")
    )
    trace = tensor_trace(
        fundamental, *(generator(a, AUTO, AUTO) for a in ("a", "b", "c"))
    )
    conjugate_trace = trace.spenso_conjugate()
    assert conjugate_trace.spenso_conjugate().to_expression() == trace.to_expression()
    assert conjugate_trace.simplify_color().to_expression().to_expression() == (
        loop.spenso_conjugate().simplify_color().to_expression().to_expression()
    )
    print(f"SU({colors}): compact/explicit color conjugation and norm passed")

    cyclic, sym, antisym, factors = S(
        "spenso::cyclic", "spenso::sym", "spenso::antisym", "factors__"
    )
    symmetric_trace = TensorExpression(
        trace.to_expression().replace(cyclic(factors), sym(factors))
    )
    antisymmetric_trace = TensorExpression(
        trace.to_expression().replace(cyclic(factors), antisym(factors))
    )
    # Cyclic invariance reduces the six permutations of three generators to
    # two orientations. Their half-sum is real and half-difference imaginary.
    assert (
        symmetric_trace.spenso_conjugate().to_expression()
        == symmetric_trace.to_expression()
    )
    assert (
        antisymmetric_trace.spenso_conjugate() + antisymmetric_trace
    ).simplify_color().to_expression().to_expression() == E("0")
    assert (
        symmetric_trace.simplify_color().to_expression().to_expression()
        == ((trace + conjugate_trace) / 2)
        .simplify_color()
        .to_expression()
        .to_expression()
    )
    antisymmetric = antisymmetric_trace.simplify_color().to_expression()
    expected_antisymmetric = (
        ((trace - conjugate_trace) / 2).simplify_color().to_expression()
    )
    # Factored half-differences need not have the same literal Atom form.
    assert (
        (antisymmetric - expected_antisymmetric)
        .simplify_color()
        .to_expression()
        .to_expression()
    ) == E("0")
    print(f"SU({colors}): symmetric and antisymmetric projector conjugation passed")
