"""Exact four-dimensional checks for the routed fermion-ring ladder.

The matrix oracle reads gamma matrices from HEP-lib, but performs its contractions
with Python integers. It never calls simplify_gamma, trace4, or tracen. Finite
component assignments supplement (and do not replace) symbolic FORM equality.
"""

from fractions import Fraction
from itertools import combinations, product
from math import lcm

from symbolica import E, Expression, Replacement
from symbolica.community.spenso import (
    Representation,
    Tensor,
    TensorExpression,
    TensorLibrary,
    TensorNetwork,
)


class FermionRingComponents:
    """Compare scalar polynomials and original networks to a HEP-matrix oracle."""

    metric = (1, -1, -1, -1)
    # Rows are k1, k2, k3, k4, q. Three-loop checks omit k4. Components
    # are contravariant and span four dimensions, with no on-shell constraints.
    assignments = (
        ((2, 1, -1, 3), (1, -2, 3, 1), (3, 2, 1, -1), (2, -1, 2, 3), (1, 1, 2, 2)),
        (
            ("1/2", "-2/3", "3/5", "1/3"),
            ("-3/2", "1/5", "2/3", "4/5"),
            ("2/5", "3/2", "-1/3", "2/3"),
            ("-1/5", "4/3", "-3/2", "1/2"),
            ("4/3", "-1/2", "3/5", "-2/5"),
        ),
        (
            ("-2/3", "5/4", "1/2", "-3/2"),
            ("3/4", "-1/3", "-5/2", "2/3"),
            ("5/3", "1/4", "-2/3", "1/2"),
            ("3/2", "-5/4", "1/3", "-2/3"),
            ("-1/2", "2/3", "3/4", "5/4"),
        ),
    )

    def __init__(self, case):
        self.case = case
        if case.loops not in (3, 4):
            raise ValueError("Component assignments support three or four loops")
        assert len(case.routing) == len(case.momenta) == 2 * case.loops
        assert all(len(row) == case.loops + 1 for row in case.routing)
        components = TensorExpression.gamma(4).components(
            library=TensorLibrary.hep_lib_atom()
        )
        units = (
            (E("0"), (0, 0)),
            (E("1"), (1, 0)),
            (E("-1"), (-1, 0)),
            (E("1𝑖"), (0, 1)),
            (E("-1𝑖"), (0, -1)),
        )
        assert len(components) == 64
        exact = []
        for value in components:
            matches = [pair for atom, pair in units if value == atom]
            assert len(matches) == 1, ("non-Gaussian-unit HEP gamma component", value)
            exact.append(matches[0])
        # Public components use logical [spinor row, spinor column, Lorentz].
        self.gamma = [
            tuple(exact[16 * r + 4 * c + mu] for r in range(4) for c in range(4))
            for mu in range(4)
        ]
        identity = tuple((int(r == c), 0) for r in range(4) for c in range(4))
        for mu, nu in product(range(4), repeat=2):
            ab, ba = (
                self._multiply(self.gamma[mu], self.gamma[nu]),
                self._multiply(self.gamma[nu], self.gamma[mu]),
            )
            assert tuple((a[0] + b[0], a[1] + b[1]) for a, b in zip(ab, ba)) == tuple(
                (2 * self.metric[mu] * v[0] if mu == nu else 0, 0) for v in identity
            )
        assert sum(identity[5 * i][0] for i in range(4)) == 4

    @staticmethod
    def _multiply(left, right):
        result = []
        for r, c in product(range(4), repeat=2):
            real = imaginary = 0
            for k in range(4):
                a, b = left[4 * r + k]
                x, y = right[4 * k + c]
                real += a * x - b * y
                imaginary += a * y + b * x
            result.append((real, imaginary))
        return tuple(result)

    @staticmethod
    def _determinant(rows):
        rows = [list(row) for row in rows]
        result = Fraction(1)
        for c in range(4):
            pivot = next((r for r in range(c, 4) if rows[r][c]), None)
            if pivot is None:
                return Fraction(0)
            if pivot != c:
                rows[c], rows[pivot] = rows[pivot], rows[c]
                result = -result
            scale = rows[c][c]
            result *= scale
            rows[c] = [x / scale for x in rows[c]]
            for r in range(c + 1, 4):
                scale = rows[r][c]
                rows[r] = [a - scale * b for a, b in zip(rows[r], rows[c])]
        return result

    def routed(self, independent):
        return [
            tuple(
                sum(weight * independent[j][axis] for j, weight in enumerate(route))
                for axis in range(4)
            )
            for route in self.case.routing
        ]

    def matrix_value(self, independent):
        """Contract the vertex indices, keeping all arithmetic exactly integral."""
        denominator = lcm(*(x.denominator for row in independent for x in row))
        momenta = [
            tuple(int(x * denominator) for x in row) for row in self.routed(independent)
        ]
        slash = [
            tuple(
                (
                    sum(
                        self.metric[mu] * p[mu] * self.gamma[mu][entry][0]
                        for mu in range(4)
                    ),
                    sum(
                        self.metric[mu] * p[mu] * self.gamma[mu][entry][1]
                        for mu in range(4)
                    ),
                )
                for entry in range(16)
            )
            for p in momenta
        ]
        pair = [
            [self._multiply(self.gamma[mu], slash[i]) for mu in range(4)]
            for i in range(len(momenta))
        ]
        real = imaginary = 0
        for indices in product(range(4), repeat=self.case.loops):
            # The outgoing half reverses rung endpoints around the fermion ring.
            order = indices + indices[:1] + indices[:0:-1]
            word = pair[0][order[0]]
            for i, index in enumerate(order[1:], start=1):
                word = self._multiply(word, pair[i][index])
            sign = -1  # Closed fermion loop.
            for index in indices:
                sign *= self.metric[index]
            real += sign * sum(word[5 * i][0] for i in range(4))
            imaginary += sign * sum(word[5 * i][1] for i in range(4))
        assert imaginary == 0
        return Fraction(real, denominator ** len(momenta))

    def _four_dimensions(self, expression):
        if isinstance(self.case.dimension, int):
            return expression
        return expression.replace(self.case.dimension, E("4"))

    def network_value(self, expression, independent):
        """Evaluate an unsimplified gamma network with exact HEP-library inputs."""
        lorentz = Representation.mink(4)
        library = TensorLibrary.hep_lib_atom()
        metric = Tensor.sparse(TensorExpression.g(lorentz), Expression)
        for axis, sign in enumerate(self.metric):
            metric[axis, axis] = E(str(sign))
        library.register(metric)
        for momentum, values in zip(
            self.case.momenta, self.routed(independent), strict=True
        ):
            tensor = Tensor.sparse(
                TensorExpression(self._four_dimensions(momentum)), Expression
            )
            for axis, value in enumerate(values):
                tensor[axis] = E(str(value))
            library.register(tensor)
        if hasattr(expression, "to_expression"):
            expression = expression.to_expression()
        expression = self._four_dimensions(expression)
        network = TensorNetwork(expression, library=library)
        network.execute(library=library)
        return network.result_scalar()

    def check(self, scalar_results, *, networks=None):
        """Check named scalar results and original networks at three points.

        Results use the case's scalar-product basis, with D specialized to four.
        These finite checks supplement the separate exact symbolic FORM equality.
        """
        assert scalar_results
        if networks is None:
            networks = {
                "original compact network": self.case.source,
                "original explicit metric network": self.case.explicit,
            }
        labels = [f"k{i}" for i in range(1, self.case.loops + 1)] + ["q"]
        rows = []
        for sample, values in enumerate(self.assignments, start=1):
            independent = [
                tuple(Fraction(x) for x in row)
                for row in values[: self.case.loops] + values[-1:]
            ]
            minors = {
                indices: self._determinant([independent[i] for i in indices])
                for indices in combinations(range(len(independent)), 4)
            }
            spanning, determinant = next(
                (indices, value) for indices, value in minors.items() if value
            )
            expected_fraction = self.matrix_value(independent)
            assert (
                expected_fraction != 0
            )  # Ensures the overall closed-loop sign is tested.
            expected = E(str(expected_fraction))
            substitutions = []
            for (i, j), symbol in self.case.scalars.items():
                dot = sum(
                    sign * x * y
                    for sign, x, y in zip(self.metric, independent[i], independent[j])
                )
                substitutions.append(Replacement(symbol, E(str(dot))))
            observed = {}
            for route, expression in scalar_results.items():
                result = (
                    self._four_dimensions(expression)
                    .replace_multiple(substitutions)
                    .expand()
                )
                assert result == expected, (
                    sample,
                    route,
                    result.format_plain(),
                    str(expected_fraction),
                )
                observed[route] = result.format_plain()
            for route, expression in networks.items():
                result = self.network_value(expression, independent)
                assert result == expected, (
                    sample,
                    route,
                    result.format_plain(),
                    str(expected_fraction),
                )
                observed[route] = result.format_plain()
            rows.append(
                {
                    "sample": sample,
                    "independent_vectors": {
                        name: [str(x) for x in row]
                        for name, row in zip(labels, independent, strict=True)
                    },
                    "spanning_vectors": [labels[i] for i in spanning],
                    "spanning_determinant": str(determinant),
                    "matrix_oracle": str(expected_fraction),
                    "observed": observed,
                    "exact": True,
                }
            )
        return rows
