"""FeynCalc trace/chain Fierz contractions through shared Idenso algebra.

Reference: https://feyncalc.github.io/FeynCalcBook/SUNSimplify.html
"""

from symbolica import E
from symbolica.community.spenso import (
    ColorSimplifySettings,
    Representation,
    TensorExpression,
    TensorName,
)

for colors in (2, 3, 5):
    generator = TensorExpression.t(colors**2 - 1, colors)
    fundamental = Representation.cof(colors)
    adjoint = Representation.coad(colors**2 - 1)
    loop = (
        generator("b", "k", "l") * generator("a", "l", "m") * generator("c", "m", "k")
    )
    trace = loop.collect_chains(fundamental)
    mixed = TensorExpression(
        generator("a", "i", "j").to_expression() * trace.to_expression()
    )
    identity = TensorExpression(
        TensorName.g().to_expression()(
            fundamental("i").to_expression(), fundamental.dual()("j").to_expression()
        )
    )
    expected = TensorExpression(
        (generator("c", "i", "k") * generator("b", "k", "j")).to_expression() / 2
        - identity.to_expression() * adjoint.g("b", "c").to_expression() / (4 * colors)
    )
    reduced = mixed.simplify_color().to_cof_dimension_invariants()
    assert reduced.to_expression().expand() == (
        expected.simplify_color().to_expression().expand()
    )
    assert reduced.simplify_color().to_expression() == reduced.to_expression()
    disabled = ColorSimplifySettings(
        evaluate_traces=False, expand_cross_chain_fierz=False
    )
    preserved = mixed.simplify_color(disabled)
    assert preserved.to_expression() != reduced.to_expression()
    assert (
        preserved.simplify_color().to_cof_dimension_invariants().to_expression()
        == reduced.to_expression()
    )
    # Close the fundamental endpoints against the opposite generator word.
    # This is the positive norm of the three-generator trace, summed over a,b,c.
    other_trace = trace.spenso_conjugate()
    closed = TensorExpression(trace.to_expression() * other_trace.to_expression())
    trace_norm = (
        closed.simplify_color().to_cof_dimension_invariants().simplify_metrics()
    )
    assert trace_norm.is_scalar
    assert trace_norm.to_expression() == E(str((colors**2 - 1) * (colors**2 - 2))) / (
        8 * colors
    )
    if colors == 3:
        # Check every free-index component of the mixed identity directly.
        # The dual fundamental metric must retain both tensor ports in the sum.
        difference = (mixed - expected).undo_trace().undo_chain().to_network()
        difference.execute()
        components = difference.result_tensor()[:]
        assert len(components) == colors**2 * (colors**2 - 1) ** 2
        assert all(abs(complex(value)) < 1e-12 for value in components)
        # Independently contract explicit SU(3) matrices, without applying
        # symbolic color simplification to the original trace product.
        network = closed.undo_trace().undo_chain().to_network()
        network.execute()
        assert abs(complex(network.result_scalar()) - 7 / 3) < 1e-12
    print(f"SU({colors}): mixed Fierz identity and trace norm passed")
