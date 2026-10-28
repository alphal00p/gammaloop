"""Factorized contraction route for the existing ladder benchmark fixtures."""

from math import prod
from time import process_time_ns

from symbolica import E
from symbolica.community.tensor import TensorExpression


def prepare(case):
    # Import only for this route: the paired baseline predates TensorRule.
    from symbolica.community.tensor import TensorRule

    body = case.fixture["ladder_vertex_rule"] if hasattr(case, "fixture") else case.rule
    nodes = {next(iter(node)): node for node in case.source.to_expression()}
    patterns = (
        case.fixture["ladder_native_patterns"]
        if hasattr(case, "fixture")
        else case.patterns
    )
    return [
        (nodes[E(str(i))], TensorRule(pattern, body, rhs_cache_size=1000))
        for i, pattern in enumerate(patterns, 1)
    ]


def reduce(case, prepared, observer=None, *, order_policy="graph"):
    """Apply local rules, contract the factored graph, then explicitly expand."""
    start = process_time_ns() if observer is not None else 0
    factors = {}
    for i, (node, rule) in enumerate(prepared, 1):
        # The registered metric normalizer canonicalizes compact products while
        # the typed rule constructs its result; resolve only this local vertex.
        value = TensorExpression(node).replace(rule)
        factors[i] = value.to_expression()
    expression = prod(factors.values())
    actual = list(expression)
    if len(actual) != len(factors):
        raise ValueError("This ladder fixture no longer has one factor per vertex")
    order = (
        [actual.index(factors[i]) for i in case.order]
        if order_policy == "explicit"
        else None
    )
    value = TensorExpression(expression)
    replaced = process_time_ns() if observer is not None else 0
    result = (
        value.contract(collect_chains=False, collect_traces=False)
        if order is None
        else value.contract(order=order, collect_chains=False, collect_traces=False)
    )
    contracted = process_time_ns() if observer is not None else 0
    expanded = result.expand().to_expression()
    if hasattr(case, "replacements"):
        expanded = expanded.replace_multiple(case.replacements)
    materialized = process_time_ns() if observer is not None else 0
    if observer is not None:
        from tensor_benchmark_cases import completed

        observer.update(
            completion=completed(contract=result),
            stages=[
                {"phase": "rule_application_and_admission", "cpu_ns": replaced - start},
                {"phase": "contraction", "cpu_ns": contracted - replaced},
                {
                    "phase": "explicit_materialization",
                    "cpu_ns": materialized - contracted,
                },
            ],
            order_policy=order_policy,
            vertex_order=list(case.order) if order is not None else None,
            normalized_factor_order=order,
            factorized_bytes=result.get_byte_size(),
            final_terms=len(list(expanded.terms())),
        )
    return expanded
