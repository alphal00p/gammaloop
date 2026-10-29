"""Compact color words use the same conjugation as explicit indexed networks."""

from symbolica import E, S
from symbolica.community.tensor import (
    AUTO,
    FactorProjector,
    Representation,
    TensorExpression,
    chain,
)
from symbolica.community.tensor import (
    trace as tensor_trace,
)

for colors in (2, 3, 5):
    generator = TensorExpression.color_t(colors**2 - 1, colors)
    fundamental = Representation.cof(colors)
    explicit = generator("a", "i", "k") * generator("b", "k", "j")
    explicit_conjugate = explicit.dirac_adjoint()
    expected_explicit = generator("a", "k", "i") * generator("b", "j", "k")
    assert explicit.rank == explicit_conjugate.rank == 4
    assert explicit_conjugate.to_expression() == expected_explicit.to_expression()
    assert (
        explicit_conjugate.dirac_adjoint().to_expression() == explicit.to_expression()
    )
    # These k labels represent independent sums. Flattening them without scopes
    # must still reject four compatible occurrences rather than contract them.
    try:
        TensorExpression(explicit.to_expression() * explicit_conjugate.to_expression())
    except ValueError as error:
        assert "more than two compatible tensor ports" in str(error)
    else:
        raise AssertionError("unscoped independent dummy copies were accepted")
    ket = explicit.wrap_indices(S("color_norm_ket"), dummies_only=True)
    bra = explicit_conjugate.wrap_indices(S("color_norm_bra"), dummies_only=True)
    assert ket.structure.axes == explicit.structure.axes
    assert bra.structure.axes == explicit_conjugate.structure.axes
    explicit_norm = ket * bra
    color_settings = dict(
        gamma=False, color=True, color_substitute_cof_dimension_invariants=True
    )
    assert explicit_norm.rank == 0
    assert explicit_norm.simplify_algebra(**color_settings).contract(
        collect_chains=False, collect_traces=False
    ).to_expression() == E(str((colors**2 - 1) ** 2)) / (4 * colors)
    if colors == 3:
        # HEP's explicit SU(3) matrices are independent of the color identities.
        network = explicit_norm.to_network()
        network.execute()
        assert abs(complex(network.result_scalar()) - 16 / 3) < 1e-12
    word = chain(
        fundamental("i"),
        fundamental.dual()("j"),
        generator("a", AUTO, AUTO),
        generator("b", AUTO, AUTO),
    )
    conjugate = word.dirac_adjoint()
    expected = chain(
        fundamental("j"),
        fundamental.dual()("i"),
        generator("b", AUTO, AUTO),
        generator("a", AUTO, AUTO),
    )
    assert conjugate.to_expression() == expected.to_expression()
    assert conjugate.dirac_adjoint().to_expression() == word.to_expression()
    assert word.dirac_adjoint().to_expression() == conjugate.to_expression()
    assert conjugate.simplify_algebra(gamma=False, color=True).contract(
        collect_chains=False, collect_traces=False
    ).to_expression() == (
        explicit.dirac_adjoint()
        .simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
    )
    norm = (
        (word * conjugate)
        .simplify_algebra(
            gamma=False, color=True, color_substitute_cof_dimension_invariants=True
        )
        .contract(collect_chains=False, collect_traces=False)
        .contract(collect_chains=False, collect_traces=False)
    )
    assert norm.is_scalar
    assert norm.to_expression() == E(str((colors**2 - 1) ** 2)) / (4 * colors)

    loop = (
        generator("a", "i", "j") * generator("b", "j", "k") * generator("c", "k", "i")
    )
    trace = tensor_trace(
        fundamental, *(generator(a, AUTO, AUTO) for a in ("a", "b", "c"))
    )
    conjugate_trace = trace.dirac_adjoint()
    assert conjugate_trace.dirac_adjoint().to_expression() == trace.to_expression()
    assert conjugate_trace.simplify_algebra(gamma=False, color=True).contract(
        collect_chains=False, collect_traces=False
    ).to_expression() == (
        loop.dirac_adjoint()
        .simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
    )
    print(f"SU({colors}): compact/explicit color conjugation and norm passed")

    word_factors = [generator(a, AUTO, AUTO) for a in ("a", "b", "c")]
    symmetric_trace = tensor_trace(
        fundamental, FactorProjector.symmetric(*word_factors)
    )
    antisymmetric_trace = tensor_trace(
        fundamental, FactorProjector.antisymmetric(*word_factors)
    )
    # Cyclic invariance reduces the six permutations of three generators to
    # two orientations. Their half-sum is real and half-difference imaginary.
    assert (
        symmetric_trace.dirac_adjoint().to_expression()
        == symmetric_trace.to_expression()
    )
    assert (antisymmetric_trace.dirac_adjoint() + antisymmetric_trace).simplify_algebra(
        gamma=False, color=True
    ).contract(collect_chains=False, collect_traces=False).to_expression() == E("0")
    assert (
        symmetric_trace.simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
        == ((trace + conjugate_trace) / 2)
        .simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
    )
    antisymmetric = antisymmetric_trace.simplify_algebra(
        gamma=False, color=True
    ).contract(collect_chains=False, collect_traces=False)
    expected_antisymmetric = (
        ((trace - conjugate_trace) / 2)
        .simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
    )
    # Factored half-differences need not have the same literal Atom form.
    assert (
        (antisymmetric - expected_antisymmetric)
        .simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
    ) == E("0")
    print(f"SU({colors}): symmetric and antisymmetric projector conjugation passed")

# Scalar coefficients inside a compact chain conjugate before multiplication.
z = S("color_weight_z")
generator = TensorExpression.color_t(8, 3)
fundamental = Representation.cof(3)
weighted = chain(
    fundamental("i"),
    fundamental.dual()("j"),
    z**2 * generator("a", AUTO, AUTO),
    generator("b", AUTO, AUTO),
)
weighted_conjugate = weighted.dirac_adjoint()
assert weighted_conjugate.dirac_adjoint().to_expression() == weighted.to_expression()
weighted_norm = TensorExpression(
    (weighted * weighted_conjugate).to_expression().replace(z, E("1+2𝑖"))
)
weighted_network = weighted_norm.to_network()
weighted_network.execute()
assert abs(complex(weighted_network.result_scalar()) - 400 / 3) < 1e-12
print("SU(3): explicit scoped and weighted component norms passed")
