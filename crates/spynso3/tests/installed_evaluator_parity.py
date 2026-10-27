"""Expression/tensor evaluator parity and logical-layout regressions."""

import inspect
import unittest

import numpy as np
from symbolica import E, Evaluator, Expression, FunctionDefinition, S
from symbolica.community import tensor as sp


class EvaluatorParityTests(unittest.TestCase):
    def test_factory_and_evaluation_signatures_match(self):
        self.assertEqual(
            inspect.signature(sp.Tensor.evaluator),
            inspect.signature(Expression.evaluator),
        )
        for method in (
            "evaluate",
            "evaluate_complex",
            "jit_compile",
            "set_real_params",
            "compile",
        ):
            self.assertEqual(
                inspect.signature(getattr(sp.TensorEvaluator, method)),
                inspect.signature(getattr(Evaluator, method)),
            )

    def test_expression_algebra_signatures_match(self):
        for method in (
            "apart",
            "cancel",
            "collect",
            "collect_by_coefficient",
            "collect_factors",
            "collect_horner",
            "collect_num",
            "collect_symbol",
            "derivative",
            "expand",
            "expand_num",
            "factor",
            "map",
            "replace_multiple",
            "together",
        ):
            with self.subTest(method=method):
                self.assertEqual(
                    inspect.signature(getattr(sp.TensorExpression, method)),
                    inspect.signature(getattr(Expression, method)),
                )

    def test_functions_options_and_logical_layout(self):
        x, y, f = S("evaluator_parity::x", "evaluator_parity::y", "evaluator_parity::f")
        a = sp.TensorName("evaluator_parity::A")(
            7, sp.Representation.euc(2), sp.Representation.mink(3)
        )
        definition = FunctionDefinition(f, [y], y**2 + 1)
        tensor = sp.Tensor.dense(a, [f(x), x, 0, x**2, 2 * x, 3]).permute_axes([1, 0])
        for sparse in (False, True):
            source = tensor.to_sparse() if sparse else tensor
            for jit in (False, True):
                for direct in (False, True):
                    options = {
                        "functions": [definition],
                        "jit_compile": jit,
                        "direct_translation": direct,
                        "cpe_iterations": 1,
                        "n_cores": 1,
                        "jit_optimization_level": 1,
                        "max_horner_scheme_variables": 20,
                        "max_common_pair_cache_entries": 100,
                        "max_common_pair_distance": 10,
                    }
                    evaluator = source.evaluator([x], **options)
                    scalar = Expression.evaluator_multiple(source[:], [x], **options)
                    self.assertEqual(evaluator.output_shape, (3, 2))
                    self.assertIsInstance(evaluator.scalar_evaluator, Evaluator)
                    evaluator.set_real_params([0])
                    scalar.set_real_params([0])
                    points = np.array([[2.0], [-1.0]])
                    for method in ("evaluate", "evaluate_complex"):
                        expected = getattr(scalar, method)(points)
                        result = getattr(evaluator, method)(points)
                        for row, values in zip(result, expected, strict=True):
                            np.testing.assert_allclose(row[:], values)
                            self.assertEqual(row.structure, source.structure)
                            self.assertEqual(
                                row.expression().to_expression(),
                                source.expression().to_expression(),
                            )
                            self.assertEqual(
                                row.expression().arguments,
                                source.expression().arguments,
                            )
                    # Flat ArrayLike input has Symbolica's batch interpretation.
                    self.assertEqual(len(evaluator.evaluate([2.0, -1.0])), 2)
                    with self.assertRaises(ValueError):
                        evaluator.evaluate([[1.0, 2.0]])

    def test_zero_scalar_and_complex_constants(self):
        for expression, expected in ((E("0"), 0j), (E("3"), 3 + 0j), (E("1i"), 1j)):
            tensor = sp.TensorExpression(expression).to_tensor()
            for jit in (False, True):
                evaluator = tensor.evaluator([], jit_compile=jit)
                result = evaluator.evaluate_complex([[]])[0]
                self.assertEqual(result.shape, ())
                self.assertEqual(result[:], [expected])
                self.assertEqual(result.structure, tensor.structure)


if __name__ == "__main__":
    unittest.main()
