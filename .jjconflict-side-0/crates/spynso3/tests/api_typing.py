"""Static regressions: run `ty check` with the installed community Python environment."""

from collections.abc import Iterator
from typing import Literal, assert_type

import linnet
from symbolica import Expression, S
from symbolica.community import tensor as sp


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
    assert_type(settings.tensor_view, Literal["interactive", "matrix"])
    assert_type(tensor.to_html(settings=settings), str)
    assert_type(representation("mu"), sp.Slot)
    assert_type(representation(1), sp.Slot)
    assert_type(representation(scalar), Expression | sp.Slot)
    assert_type(expression.structure.slots, tuple[sp.Representation | sp.Slot, ...])
    symbolic_group = sp.FactorProjector.symmetric(expression, expression)
    component_group = sp.FactorProjector.antisymmetric(expression, tensor)
    assert_type(symbolic_group, sp.FactorProjector[sp.TensorExpression])
    assert_type(component_group, sp.FactorProjector[sp.TensorNetwork])
    assert_type(
        sp.chain(representation("i"), representation("j"), symbolic_group),
        sp.TensorExpression,
    )
    assert_type(sp.trace(representation, component_group), sp.TensorNetwork)
    assert_type(
        sp.FactorProjector.cyclic(symbolic_group, component_group),
        sp.FactorProjector[sp.TensorNetwork],
    )
    assert_type(expression.expand_projectors(), sp.TensorExpression)
    assert_type(sp.TensorExpression.gamma0(4), sp.TensorExpression)
    assert_type(sp.TensorExpression.charge_conjugation(4), sp.TensorExpression)
    assert_type(sp.TensorExpression.levi_civita(representation), sp.TensorExpression)
    assert_type(sp.TensorPattern.chain(S("a_"), S("b_"), S("fs___")), sp.TensorPattern)
    assert_type(sp.Nc, Expression)
    assert_type(expression.structure, sp.TensorStructure)
    assert_type(tensor.structure, sp.TensorStructure)
    assert_type(network.structure, sp.TensorStructure)
    config = linnet.RenderConfig()
    assert_type(network.render(config=config), str)
    assert_type(network.to_linnest(config=config), str)
    assert_type(network.to_html(config=config), str)
    assert_type(expression.structure.name, sp.TensorName | None)
    assert_type(expression.structure.arguments, tuple[Expression, ...])
    assert_type(expression.structure.shape, tuple[int | Expression, ...])
    assert_type(representation.name, sp.RepresentationName)
    assert_type(representation.name.metric_sign(0), int)
    assert_type(representation("mu").index, Expression)
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
    assert_type(
        expression.contract_ports(expression, left=0, right=0), sp.TensorExpression
    )
    assert_type(expression.contract_ports(network, left=0, right=0), sp.TensorNetwork)
    assert_type(
        expression.compose(expression, left=(0, 1), right=(0, 1)), sp.TensorExpression
    )
    assert_type(expression.compose(tensor, left=(0, 1), right=(0, 1)), sp.TensorNetwork)
    assert_type(expression.contract(), sp.AliasedTensorExpression)
    assert_type(expression.contract().contraction_complete, bool)
    assert_type(expression.contract().expand(), sp.TensorExpression)
    assert_type(expression.trace(), sp.TensorExpression)
    assert_type(expression**2, sp.TensorExpression)
    assert_type(2**expression, sp.TensorExpression)
    assert_type(bool(expression), bool)
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
    assert_type(expression.to_expression().factor(), Expression)
    assert_type(expression.to_expression().collect_num(), Expression)
    assert_type(expression.to_expression().collect_factors(), Expression)
    assert_type(
        expression.reduce_algebra(
            sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
        ).contract(sp.ContractSettings(collect_chains=False, collect_traces=False)),
        sp.AliasedTensorExpression,
    )
    assert_type(
        expression.reduce_algebra(
            sp.AlgebraSettings(color=sp.ColorSimplifySettings())
        ).contract(sp.ContractSettings(collect_chains=False, collect_traces=False)),
        sp.AliasedTensorExpression,
    )
    assert_type(expression.replace(sp.TensorRule(expression, 0)), sp.TensorExpression)
    assert_type(expression.to_expression().derivative(scalar), Expression)
    assert_type(
        expression.reduce_algebra(sp.AlgebraSettings.hep()), sp.AliasedTensorExpression
    )
    assert_type(expression.reindex("mu", sp.AUTO), sp.TensorExpression)
    assert_type(expression.rename_indices({"mu": "nu"}), sp.TensorExpression)
    assert_type(expression.permute_axes([1, 0]), sp.TensorExpression)
    assert_type(tensor.permute_axes([1, 0]), sp.Tensor)
    assert_type(network.permute_axes([1, 0]), sp.TensorNetwork)
    assert_type(network.contract_ports(expression, left=0, right=0), sp.TensorNetwork)
    assert_type(tensor.rename_indices({representation("mu"): "nu"}), sp.Tensor)
    assert_type(expression("mu", intern="flattened"), sp.TensorExpression)
    intern_indices: Literal["indices"] = "indices"
    intern_flattened: Literal["flattened"] = "flattened"
    for intern in (None, intern_indices, intern_flattened):
        assert_type(sp.TensorExpression(expression, intern=intern), sp.TensorExpression)
        assert_type(expression("mu", intern=intern), sp.TensorExpression)
        assert_type(expression.index("mu", intern=intern), sp.TensorExpression)
        assert_type(expression.reindex("mu", intern=intern), sp.TensorExpression)
        assert_type(
            expression.rename_indices({"mu": "nu"}, intern=intern), sp.TensorExpression
        )
        assert_type(tensor("mu", intern=intern), sp.TensorNetwork)
        assert_type(network("mu", intern=intern), sp.TensorNetwork)
    assert_type(expression.to_tensor(library), sp.Tensor)
    assert_type(network.to_tensor(library), sp.Tensor)
    assert_type(network.step(library), sp.TensorNetwork)
    assert_type(network.status, sp.ExecutionStatus)
    assert_type(network.status.ready_operations, tuple[str, ...])
    assert_type(network.status.complete, bool)
    assert_type(tensor.dtype, type[float] | type[complex] | type[Expression])
    assert_type(tensor.storage, Literal["dense", "sparse"])
    assert_type(tensor.copy(), sp.Tensor)
    assert_type(tensor.to_sparse(), sp.Tensor)
    assert_type(tensor.to_dense(), sp.Tensor)
    assert_type(tensor.map_components(lambda value: value * 2), sp.Tensor)
    assert_type(expression.shape, tuple[int | Expression, ...])
    assert_type(network.shape, tuple[int | Expression, ...])
    assert_type(tensor.shape, tuple[int, ...])
    assert_type(
        sp.TensorExpression(expression.to_expression(), intern="flattened"),
        sp.TensorExpression,
    )
    assert_type(library["A"], sp.Tensor)
    assert_type(library[expression], sp.Tensor)
    assert_type(expression in library, bool)
    assert_type("A" in library, bool)
    assert_type(library.get(expression), sp.Tensor | None)
    assert_type(library.get("A", None), sp.Tensor | None)
    assert_type(library.get("A", tensor), sp.Tensor)
    fallback: list[int] = []
    assert_type(library.get("A", default=fallback), sp.Tensor | list[int])
    assert_type(library.keys(), list[sp.TensorExpression])
    assert_type(library.values(), list[sp.Tensor])
    assert_type(library.items(), list[tuple[sp.TensorExpression, sp.Tensor]])
    assert_type(iter(library), Iterator[sp.TensorExpression])
    assert_type(library.to_html(), str)
    assert_type(library.to_html(settings=settings), str)
    assert_type(tensor.expression(), sp.TensorExpression)
    assert_type(network.expression(), sp.TensorExpression)
    assert_type(evaluator.evaluate([[2.0]]), list[sp.Tensor])
    assert_type(evaluator.evaluate_complex([[2j]]), list[sp.Tensor])
    assert_type(compiled.evaluate_complex([[2j]]), list[sp.Tensor])
    assert_type(compiled.evaluate([[2.0]]), list[sp.Tensor])
    assert_type(evaluator.parameters, list[Expression])
    assert_type(compiled.parameters, list[Expression])
    assert_type(evaluator.output_shape, tuple[int, ...])
    assert_type(compiled.output_shape, tuple[int, ...])
    policy: sp.ExecutionMode = sp.ExecutionMode.All
    assert_type(int(policy), int)
    settings = sp.GammaSimplifySettings(
        chain_ordering=sp.GammaChainOrdering.Canonical, gamma0=True, conjugate=True
    )
    assert_type(settings.output, Literal["reduced", "chains"])
    assert_type(settings.gamma0, bool)
    assert_type(settings.conjugate, bool)
    assert_type(settings.chain_ordering, sp.GammaChainOrdering)
    name = sp.TensorName("typing::Jbar", print={"typst": "macron(J)"})
    assert_type(name(representation), sp.TensorExpression)
    function = sp.BroadcastFunction("typing::f")
    assert_type(function(expression), sp.TensorExpression)
    assert_type(function(tensor), sp.TensorNetwork)
    assert_type(function(S("typing::x")), Expression)
    callbacks = sp.TensorFunctionLibrary()
    callbacks.register(function, lambda value: value * value)
