"""Exact component checks for the physical gluon-ring ladder.

The component oracle uses the three-term cyclic vertex formula directly, not
GluonLadder.RULE_TERMS or its symbolic replacement rules. HEP-library contraction
then evaluates the original opaque-vertex graph, with exact rational components.
"""

from fractions import Fraction
from itertools import product

from fermion_ladder_validation import FermionRingComponents
from symbolica import E, Expression, Replacement
from symbolica.community.spenso import (
    Representation,
    Tensor,
    TensorExpression,
    TensorLibrary,
    TensorNetwork,
)


class GluonRingComponents:
    """Evaluate scalar results and original networks at exact four-vector inputs."""

    metric = (1, -1, -1, -1)
    # Share only the rational test inputs; no fermion/gamma algebra is used.
    assignments = FermionRingComponents.assignments

    def __init__(self, case):
        self.case = case
        assert case.loops == 4
        assert len(case.vertex_templates) == len(case.vertices) == 8

    def _four_dimensions(self, expression):
        if isinstance(self.case.dimension, int):
            return expression
        return expression.replace(self.case.dimension, E("4"))

    @classmethod
    def _vertex_component(cls, momenta, a, b, c):
        p, q, r = momenta
        return (
            (cls.metric[a] * (p[c] - q[c]) if a == b else 0)
            + (cls.metric[b] * (q[a] - r[a]) if b == c else 0)
            + (cls.metric[c] * (r[b] - p[b]) if c == a else 0)
        )

    def _library(self, independent):
        library = TensorLibrary.hep_lib_atom()
        metric = Tensor.sparse(TensorExpression.g(Representation.mink(4)), Expression)
        for axis, sign in enumerate(self.metric):
            metric[axis, axis] = E(str(sign))
        library.register(metric)
        for template, (routes, _ports) in zip(
            self.case.vertex_templates, self.case.vertices, strict=True
        ):
            momenta = [
                tuple(
                    sum(weight * independent[j][axis] for j, weight in enumerate(row))
                    for axis in range(4)
                )
                for row in routes
            ]
            assert all(
                sum(components) == 0 for components in zip(*momenta, strict=True)
            )
            vertex = Tensor.sparse(
                TensorExpression(self._four_dimensions(template)), Expression
            )
            for indices in product(range(4), repeat=3):
                component = self._vertex_component(momenta, *indices)
                if component:
                    vertex[indices] = E(str(component))
            library.register(vertex)
        return library

    def check(self, scalar_results, *, networks=None):
        """Check Idenso/FORM scalars against unchanged opaque vertex networks."""
        assert scalar_results
        if networks is None:
            networks = {
                "original compact network": self.case.source,
                "original explicit metric network": self.case.explicit,
            }
        rows = []
        for sample, values in enumerate(self.assignments, start=1):
            independent = [tuple(Fraction(x) for x in row) for row in values]
            substitutions = [
                Replacement(
                    symbol,
                    E(
                        str(
                            sum(
                                sign * x * y
                                for sign, x, y in zip(
                                    self.metric,
                                    independent[i],
                                    independent[j],
                                    strict=True,
                                )
                            )
                        )
                    ),
                )
                for (i, j), symbol in self.case.scalars.items()
            ]
            library = self._library(independent)
            observed = {}
            expected = None
            for route, expression in networks.items():
                if hasattr(expression, "to_expression"):
                    expression = expression.to_expression()
                network = TensorNetwork(
                    self._four_dimensions(expression), library=library
                )
                network.execute(library=library)
                value = network.result_scalar()
                assert value != E("0")
                if expected is None:
                    expected = value
                assert value == expected, (sample, route, value, expected)
                observed[route] = value.format_plain()
            assert expected is not None
            for route, expression in scalar_results.items():
                value = (
                    self._four_dimensions(expression)
                    .replace_multiple(substitutions)
                    .expand()
                )
                assert value == expected, (sample, route, value, expected)
                observed[route] = value.format_plain()
            rows.append(
                {
                    "sample": sample,
                    "independent_vectors": {
                        name: [str(x) for x in row]
                        for name, row in zip(
                            ("k1", "k2", "k3", "k4", "q"), independent, strict=True
                        )
                    },
                    "hep_vertex_oracle": expected.format_plain(),
                    "observed": observed,
                    "exact": True,
                }
            )
        return rows
