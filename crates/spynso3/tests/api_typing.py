"""Static regressions: run `ty check` with the installed community Python environment."""

from collections.abc import Iterator
from typing import assert_type

from symbolica import Expression, S
from symbolica.community import spenso as sp


def check_types(
    expression: sp.TensorExpression,
    tensor: sp.Tensor,
    network: sp.TensorNetwork,
    scalar: Expression,
    representation: sp.Representation,
    library: sp.TensorLibrary,
    evaluator: sp.TensorEvaluator,
    compiled: sp.CompiledTensorEvaluator,
) -> None:
    settings = sp.DisplaySettings(tensor_view="matrix")
    assert_type(settings.tensor_view, str)
    assert_type(tensor.to_html(settings=settings), str)
    assert_type(representation("mu"), sp.Slot)
    assert_type(representation(1), sp.Slot)
    assert_type(representation(scalar), Expression | sp.Slot)
    assert_type(expression.interface, tuple[sp.Representation | sp.Slot, ...])
    assert_type(expression[0], list[int])
    assert_type(expression[[0, 1]], int)
    assert_type(expression[:], list[list[int]])
    assert_type(tensor[0], Expression | float | complex)
    assert_type(tensor[0, 1], Expression | float | complex)
    assert_type(tensor[:], list[Expression | float | complex])
    assert_type(iter(tensor), Iterator[Expression | float | complex])
    tensor[0, 1] = 2.0
    tensor[0] = 1j
    assert_type(expression + scalar, sp.TensorExpression)
    assert_type(expression * 2.0, sp.TensorExpression)
    assert_type(expression * 1j, sp.TensorExpression)
    assert_type(expression + tensor, sp.TensorNetwork)
    assert_type(expression * network, sp.TensorNetwork)
    assert_type(tensor * expression, sp.TensorNetwork)
    assert_type(network + 2.0, sp.TensorNetwork)
    assert_type(sp.dot(expression, expression), sp.TensorExpression)
    assert_type(sp.dot(tensor, expression), sp.TensorNetwork)
    assert_type(sp.dot(expression, network), sp.TensorNetwork)
    assert_type(expression.outer(expression), sp.TensorExpression)
    assert_type(expression.outer(tensor), sp.TensorNetwork)
    assert_type(expression.contract(expression, left=0, right=0), sp.TensorExpression)
    assert_type(expression.contract(network, left=0, right=0), sp.TensorNetwork)
    assert_type(
        expression.compose(expression, left=(0, 1), right=(0, 1)), sp.TensorExpression
    )
    assert_type(expression.compose(tensor, left=(0, 1), right=(0, 1)), sp.TensorNetwork)
    assert_type(expression.trace(), sp.TensorExpression)
    assert_type(tensor.trace(), sp.TensorNetwork)
    assert_type(sp.trace(representation, expression), sp.TensorExpression)
    assert_type(sp.trace(representation, tensor, expression), sp.TensorNetwork)
    assert_type(sp.trace(representation, expression, tensor), sp.TensorNetwork)
    assert_type(sp.trace(representation, tensor, 2.0), sp.TensorNetwork)
    # An arbitrary mixed variadic sequence conservatively retains both possibilities.
    assert_type(
        sp.trace(representation, expression, expression, tensor),
        sp.TensorExpression | sp.TensorNetwork,
    )
    start, end = representation("i"), representation("j")
    assert_type(sp.chain(start, end, expression), sp.TensorExpression)
    assert_type(sp.chain(start, end, tensor, expression), sp.TensorNetwork)
    assert_type(sp.chain(start, end, expression, tensor), sp.TensorNetwork)
    assert_type(expression.factor(), sp.TensorExpression)
    assert_type(expression.collect_num(), sp.TensorExpression)
    assert_type(expression.collect_factors(), sp.TensorExpression)
    assert_type(expression.simplify_gamma(), sp.TensorExpression)
    assert_type(expression.simplify_color(), sp.TensorExpression)
    assert_type(expression.cook_indices(sp.CookSettings.indices()), sp.TensorExpression)
    assert_type(library["A"], sp.TensorExpression)
    assert_type(tensor.structure(), sp.TensorExpression)
    assert_type(network.structure(), sp.TensorExpression)
    assert_type(evaluator.evaluate([[2.0]]), list[sp.Tensor])
    assert_type(evaluator.evaluate_complex([[2j]]), list[sp.Tensor])
    assert_type(compiled.evaluate_complex([[2j]]), list[sp.Tensor])
    policy: sp.ExecutionMode = sp.ExecutionMode.All
    assert_type(int(policy), int)
    settings = sp.GammaSimplifySettings(chain_ordering=sp.GammaChainOrdering.Canonical)
    assert_type(settings.chain_ordering, sp.GammaChainOrdering)
    name = sp.TensorName("typing::Jbar", print={"typst": "macron(J)"})
    assert_type(name(representation), sp.TensorExpression)
    function = sp.BroadcastFunction("typing::f")
    assert_type(function(expression), sp.TensorExpression)
    assert_type(function(tensor), sp.TensorNetwork)
    assert_type(function(S("typing::x")), Expression)
    callbacks = sp.TensorFunctionLibrary()
    callbacks.register(function, lambda value: value * value)
