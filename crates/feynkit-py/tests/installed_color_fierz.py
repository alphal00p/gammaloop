"""FeynCalc trace/chain Fierz contractions through shared Idenso algebra.

Reference: https://feyncalc.github.io/FeynCalcBook/SUNSimplify.html
"""

from symbolica import E
from symbolica.community.tensor import AUTO, Representation, TensorExpression
from symbolica.community.tensor import (
    trace as tensor_trace,
)

for colors in (2, 3, 5):
    generator = TensorExpression.color_t(colors**2 - 1, colors)
    fundamental = Representation.cof(colors)
    adjoint = Representation.coad(colors**2 - 1)
    trace = tensor_trace(
        fundamental, *(generator(a, AUTO, AUTO) for a in ("b", "a", "c"))
    )
    mixed = TensorExpression(
        generator("a", "i", "j").to_expression() * trace.to_expression()
    )
    identity = TensorExpression.g(fundamental, fundamental.dual())("i", "j")
    expected = TensorExpression(
        (generator("c", "i", "k") * generator("b", "k", "j")).to_expression() / 2
        - identity.to_expression() * adjoint.g("b", "c").to_expression() / (4 * colors)
    )
    reduced = mixed.simplify_algebra(
        gamma=False, color=True, color_substitute_cof_dimension_invariants=True
    ).contract(collect_chains=False, collect_traces=False)
    assert reduced.to_expression().expand() == (
        expected.simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
        .expand()
    )
    assert (
        reduced.simplify_algebra(gamma=False, color=True)
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
        == reduced.to_expression()
    )
    disabled = dict(
        gamma=False, color=True, color_evaluate_traces=False, color_expand_fierz=False
    )
    preserved = mixed.simplify_algebra(**disabled).contract(
        collect_chains=False, collect_traces=False
    )
    assert preserved.to_expression() != reduced.to_expression()
    assert (
        preserved.simplify_algebra(
            gamma=False, color=True, color_substitute_cof_dimension_invariants=True
        )
        .contract(collect_chains=False, collect_traces=False)
        .to_expression()
        == reduced.to_expression()
    )
    # Close the fundamental endpoints against the opposite generator word.
    # This is the positive norm of the three-generator trace, summed over a,b,c.
    other_trace = trace.dirac_adjoint()
    closed = TensorExpression(trace.to_expression() * other_trace.to_expression())
    trace_norm = (
        closed.simplify_algebra(
            gamma=False, color=True, color_substitute_cof_dimension_invariants=True
        )
        .contract(collect_chains=False, collect_traces=False)
        .contract(collect_chains=False, collect_traces=False)
    )
    assert trace_norm.is_scalar
    assert trace_norm.to_expression() == E(str((colors**2 - 1) * (colors**2 - 2))) / (
        8 * colors
    )
    if colors == 3:
        # Check every free-index component of the mixed identity directly.
        # The dual fundamental metric must retain both tensor ports in the sum.
        difference = (mixed - expected).to_network()
        difference.execute()
        components = difference.result_tensor()[:]
        assert len(components) == colors**2 * (colors**2 - 1) ** 2
        assert all(abs(complex(value)) < 1e-12 for value in components)
        # Independently contract explicit SU(3) matrices, without applying
        # symbolic color simplification to the original trace product.
        network = closed.to_network()
        network.execute()
        assert abs(complex(network.result_scalar()) - 7 / 3) < 1e-12
    print(f"SU({colors}): mixed Fierz identity and trace norm passed")
