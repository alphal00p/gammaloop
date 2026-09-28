"""Factorized contraction route for the existing ladder benchmark fixtures."""

from collections import Counter
from math import prod
from time import process_time_ns

from symbolica import E
from symbolica.community.spenso import TensorExpression


def alias_levels(result):
    """Count retained definitions/weight terms, outside every algebra clock."""
    depth = {}
    definitions = Counter()
    terms = Counter()
    for handle, body in result.aliases:
        expression = body.to_expression()
        pending = [expression]
        parents = []
        while pending:
            node = pending.pop()
            if node in depth:
                parents.append(depth[node])
            elif str(node.get_type()).rsplit(".", 1)[-1] in {"Add", "Mul", "Pow", "Fn"}:
                pending.extend(node)
        level = 1 + max(parents, default=0)
        depth[handle.to_expression()] = level
        definitions[level] += 1
        terms[level] += len(list(expression.terms()))
    return [
        {"depth": level, "definitions": count, "stored_weight_terms": terms[level]}
        for level, count in sorted(definitions.items())
    ]


def prepare(case):
    # Import only for this route: the paired baseline predates TensorRule.
    from symbolica.community.spenso import TensorRule

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


def reduce(case, prepared, observer=None):
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
    order = [actual.index(factors[i]) for i in case.order]
    value = TensorExpression(expression)
    replaced = process_time_ns() if observer is not None else 0
    result = value.contract(order=order)
    contracted = process_time_ns() if observer is not None else 0
    expanded = result.expand().to_expression()
    if hasattr(case, "replacements"):
        expanded = expanded.replace_multiple(case.replacements)
    materialized = process_time_ns() if observer is not None else 0
    if observer is not None:
        observer.update(
            stages=[
                {"phase": "rule_application_and_admission", "cpu_ns": replaced - start},
                {"phase": "contraction", "cpu_ns": contracted - replaced},
                {
                    "phase": "explicit_materialization",
                    "cpu_ns": materialized - contracted,
                },
            ],
            vertex_order=list(case.order),
            normalized_factor_order=order,
            aliased_bytes=result.get_byte_size(),
            definitions=len(result.aliases),
            alias_levels=alias_levels(result),
            final_terms=len(list(expanded.terms())),
        )
    return expanded
