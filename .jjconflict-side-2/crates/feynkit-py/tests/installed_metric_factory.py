"""Typed metrics support explicit dual representations and preserve index order."""

from symbolica import E, S
from symbolica.community.spenso import (
    PortPattern,
    Representation,
    TensorExpression,
    TensorName,
    TensorPattern,
)

for dimension in (2, 3, 5):
    fundamental = Representation.cof(dimension)
    dual = fundamental.dual()
    for left, right in ((fundamental, dual), (dual, fundamental)):
        metric = TensorExpression.g(left, right)
        indexed = metric("i", "j")
        expected = TensorName.g().to_expression()(
            left("i").to_expression(), right("j").to_expression()
        )
        assert indexed.to_expression() == expected
        assert metric("i", "i").is_scalar
        assert metric("i", "i").simplify_metrics().to_expression() == E(str(dimension))
        network = indexed.to_network()
        network.execute()
        values = network.result_tensor()[:]
        assert len(values) == dimension**2
        assert all(
            complex(value) == int(row == column)
            for (row, column), value in zip(
                (
                    (row, column)
                    for row in range(dimension)
                    for column in range(dimension)
                ),
                values,
                strict=True,
            )
        )
        i, j = S("metric_i_", "metric_j_")
        pattern = TensorPattern(
            TensorName.g(),
            ports=[PortPattern.exact(left, i), PortPattern.exact(right, j)],
        )
        assert indexed.to_expression().replace(pattern, E("17")) == E("17")

mink = Representation.mink(4)
assert (
    TensorExpression.g(mink)("mu", "nu").to_expression()
    == TensorExpression.g(mink, mink)("mu", "nu").to_expression()
)
assert TensorExpression.g(mink)("mu", "mu").simplify_metrics().to_expression() == E("4")
n, m = S("metric_N", "metric_M")
symbolic = Representation.cof(n)
assert (
    TensorExpression.g(symbolic, symbolic.dual())("i", "i")
    .simplify_metrics()
    .to_expression()
    == n
)
for left, right in (
    (Representation.cof(3), Representation.cof(4)),
    (Representation.cof(3), Representation.coad(3)),
    (symbolic, Representation.cof(m).dual()),
    (symbolic, Representation.cof(3).dual()),
):
    try:
        TensorExpression.g(left, right)
    except ValueError:
        pass
    else:
        raise AssertionError("mismatched metric representation spaces were accepted")
print(
    "Metric factory: dual port order, identity matrices, patterns and exact dimensions passed"
)
