"""Typed factorized contraction against the installed Python extension."""

import unittest

from symbolica import E
from symbolica.community import spenso as sp


class FactorizedContractionTests(unittest.TestCase):
    def test_contraction_outputs_and_dag_evaluation(self):
        rep = sp.Representation.mink(4)
        p, q, r = [
            sp.TensorName.vector(f"api_contract::{name}") for name in ("p", "q", "r")
        ]
        g = sp.TensorName.g().to_expression()
        first = g(p(rep).to_expression(), r(rep).to_expression())
        second = g(q(rep).to_expression(), r(rep).to_expression())
        expected = first + second
        i = rep("i")
        source = sp.TensorExpression(
            (p(i).to_expression() + q(i).to_expression()) * r(i).to_expression()
        )
        uncertified = sp.AliasedTensorExpression.from_expression(source)
        self.assertFalse(uncertified.contraction_complete)
        with self.assertRaises(AttributeError):
            uncertified.contraction_complete = True
        for order in (None, [0, 1], [1, 0]):
            with self.subTest(order=order):
                result = source.contract(order=order)
                self.assertIsInstance(result, sp.AliasedTensorExpression)
                self.assertTrue(result.contraction_complete)
                self.assertGreater(len(result.aliases), 0)
                self.assertEqual(result.root.rank, 0)
                self.assertEqual(result.expand().to_expression(), expected)
                rerun = result.contract()
                self.assertTrue(rerun.contraction_complete)
                self.assertEqual(rerun.expand(), result.expand())
                self.assertEqual(rerun.root.structure, result.root.structure)
                self.assertEqual(
                    [
                        (handle.to_expression(), body.to_expression())
                        for handle, body in rerun.aliases
                    ],
                    [
                        (handle.to_expression(), body.to_expression())
                        for handle, body in result.aliases
                    ],
                )
                self.assertEqual(
                    result.to_expression().to_expression().expand(), expected
                )
                self.assertEqual(
                    source.contract(order=order).expand().to_expression(),
                    expected,
                )
                evaluator = result.evaluator([first, second], iterations=1, n_cores=1)
                self.assertEqual(evaluator.evaluate([[2.0, 3.0]])[0][0], 5.0)
        for order in ([0], [0, 0], [0, 2]):
            with self.assertRaisesRegex(
                ValueError, "every normalized top-level factor"
            ):
                source.contract(order=order)
        with self.assertRaisesRegex(TypeError, "output"):
            source.contract(output="invalid")
        with self.assertRaises(TypeError):
            source.contract(left=0, right=0)

    def test_contraction_preserves_unresolved_order_metadata_and_typed_zero(self):
        rep = sp.Representation.euc(2)
        other = sp.Representation.mink(4)
        opened = sp.TensorName("api_contract::T")(7, rep, other)
        permuted = opened.permute_axes([1, 0])
        for value in (opened, permuted, 0 * opened, 0 * permuted):
            with self.subTest(zero=value.to_expression() == E("0")):
                result = value.contract()
                self.assertEqual(result.root.structure, value.structure)
                self.assertEqual(result.to_expression().structure, value.structure)
                self.assertEqual(
                    result.to_expression().to_expression(), value.to_expression()
                )


if __name__ == "__main__":
    unittest.main()
