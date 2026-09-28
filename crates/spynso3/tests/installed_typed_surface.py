"""Standalone tensor values and explicit Symbolica boundaries, with the built extension."""

import copy
import operator
import unittest

from symbolica import E, Expression, S
from symbolica.community import spenso as sp


class TypedSurfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rep = sp.Representation.euc(3)
        cls.P = sp.TensorName("typed_surface::P")
        cls.Q = sp.TensorName("typed_surface::Q")
        cls.x, cls.y = S("typed_surface::x", "typed_surface::y")

    def test_explicit_symbolica_boundary_and_checked_reconstruction(self):
        source = self.P(self.x, self.rep("j"), self.rep("i")).permute_axes([1, 0])
        self.assertNotIsInstance(source, Expression)
        self.assertNotIn(Expression, sp.TensorExpression.__mro__)
        for name in (
            "replace_multiple",
            "map",
            "derivative",
            "factor",
            "together",
            "cancel",
            "apart",
        ):
            with self.subTest(name=name):
                self.assertFalse(hasattr(source, name))
        raw = source.to_expression().replace(self.x, self.y)
        changed = sp.TensorExpression(raw, structure=source.structure)
        self.assertEqual(changed.to_expression(), raw)
        self.assertEqual(changed.structure.slots, source.structure.slots)
        self.assertEqual(changed.structure.arguments, (self.y,))
        self.assertEqual(str(changed.structure.name), str(source.structure.name))
        zero = sp.TensorExpression(E("0"), structure=source.structure)
        self.assertEqual(zero.rank, source.rank)
        self.assertEqual(zero.structure.slots, source.structure.slots)
        with self.assertRaisesRegex(ValueError, "interface"):
            sp.TensorExpression(E("1"), structure=source.structure)

    def test_tensor_equality_hash_and_truth_preserve_logical_interface(self):
        source = self.P(self.rep("j"), self.rep("i"))
        same = copy.copy(source).with_name("typed_surface::PresentationOnly")
        self.assertEqual(source, same)
        self.assertEqual(hash(source), hash(same))
        self.assertEqual({source: "value"}[same], "value")
        self.assertNotEqual(source, source.permute_axes([1, 0]))
        self.assertNotEqual(source, source.to_expression())
        self.assertNotEqual(sp.TensorExpression(E("0")), 0 * source)
        self.assertTrue(source)
        self.assertTrue(sp.TensorExpression(self.x))
        self.assertFalse(0 * source)
        self.assertFalse(sp.TensorExpression(E("0")))
        self.assertEqual(sp.TensorExpression(source), source)

    def test_mixed_scalar_arithmetic_dispatches_in_both_orders(self):
        typed = sp.TensorExpression(self.x)
        for other in (2, self.y, sp.TensorExpression(self.y)):
            raw = (
                other.to_expression()
                if isinstance(other, sp.TensorExpression)
                else other
            )
            for operation in (
                operator.add,
                operator.sub,
                operator.mul,
                operator.truediv,
            ):
                for left, right, expected in (
                    (typed, other, operation(self.x, raw)),
                    (other, typed, operation(raw, self.x)),
                ):
                    with self.subTest(
                        operation=operation.__name__, reflected=left is other
                    ):
                        result = operation(left, right)
                        self.assertIsInstance(result, sp.TensorExpression)
                        self.assertEqual(result.to_expression(), expected)
        vector = self.P(self.rep("i"))
        self.assertEqual(self.x * vector, vector * self.x)
        quotient = vector.to_expression() / typed
        self.assertIsInstance(quotient, sp.TensorExpression)
        self.assertEqual(quotient, vector / self.x)
        self.assertEqual(quotient.structure.slots, vector.structure.slots)
        with self.assertRaisesRegex(ValueError, "denominator"):
            self.x / vector
        with self.assertRaises(ValueError):
            self.x + vector

    def test_scalar_tensor_denominators_keep_internal_indices_scoped(self):
        p = sp.TensorName.vector("typed_surface::denominator_p")
        q = sp.TensorName.vector("typed_surface::denominator_q")
        slot = self.rep("i")
        numerator = p(slot)
        denominator = numerator * q(slot)
        self.assertTrue(denominator.is_scalar)
        result = numerator / denominator
        self.assertEqual(result.structure.slots, numerator.structure.slots)
        self.assertEqual(
            (numerator.to_expression() / denominator).structure.slots,
            numerator.structure.slots,
        )
        self.assertEqual(
            (0 * numerator / denominator).structure.slots, numerator.structure.slots
        )
        self.assertFalse(0 * numerator / denominator)

    def test_selected_port_contraction_retains_declared_logical_order(self):
        left = self.P(self.rep("i"), self.rep("j")).permute_axes([1, 0])
        right = self.Q(self.rep("k"), self.rep("i")).permute_axes([1, 0])
        result = left.contract_ports(right, left=1, right=0)
        self.assertEqual(
            result.structure.slots, (left.structure.slots[0], right.structure.slots[1])
        )
        self.assertEqual(result.rank, 2)
        with self.assertRaises(TypeError):
            left.contract(right, left=1, right=0)

    def test_typed_powers_and_reflected_scalar_power(self):
        scalar = sp.TensorExpression(self.x + 1)
        for exponent in (0, 1, 2, -1, E("1/2"), self.y):
            result = scalar**exponent
            self.assertIsInstance(result, sp.TensorExpression)
            self.assertEqual(result.to_expression(), (self.x + 1) ** exponent)
        exponent = sp.TensorExpression(self.x)
        self.assertEqual((2**exponent).to_expression(), 2**self.x)
        vector = self.P(self.rep("i"))
        self.assertTrue((vector**2).is_scalar)
        self.assertFalse(0 * vector)
        self.assertTrue(((0 * vector) ** 2).is_scalar)
        for invalid in (3, -1, E("1/2"), self.x, vector):
            with self.subTest(exponent=invalid), self.assertRaises(ValueError):
                vector**invalid
        opened = self.P(self.rep)
        self.assertEqual((opened**2).rank, 0)
        self.assertEqual((opened**3).rank, 1)
        self.assertEqual(((0 * opened) ** 3).structure.slots, opened.structure.slots)

    def test_explicit_evaluator_uses_bare_and_aliased_symbolic_payloads(self):
        source = sp.TensorExpression((self.x + 1) ** 3 + self.y)
        aliased = sp.AliasedTensorExpression.from_expression(source)
        points = [[2.0, 5.0], [-1.0, 7.0]]
        for value in (source, aliased):
            evaluator = value.evaluator([self.x, self.y], iterations=1, n_cores=1)
            self.assertEqual(evaluator.evaluate(points).tolist(), [[32.0], [7.0]])
        self.assertEqual(source.to_expression(), (self.x + 1) ** 3 + self.y)

    def test_typed_rule_inputs_and_callback_outputs_use_existing_schedule(self):
        source = self.P(self.rep("i"))
        target = self.Q(self.rep("i"))
        rule = sp.TensorRule(source, target)
        self.assertEqual(source.replace(rule), target)
        self.assertEqual(source.replace([rule, sp.TensorRule(target, source)]), target)
        for cache_size, count in ((0, 4), (10, 2)):
            calls = []

            def rhs(matches, calls=calls):
                calls.append(matches)
                return target

            rule = sp.TensorRule(source, rhs, rhs_cache_size=cache_size)
            repeated = self.x * source + self.y * source
            for _ in range(2):
                self.assertEqual(
                    repeated.replace(rule), self.x * target + self.y * target
                )
            self.assertEqual(len(calls), count)
        zero = source.replace(sp.TensorRule(source, sp.TensorExpression(E("0"))))
        self.assertFalse(zero)
        self.assertEqual(zero.structure.slots, source.structure.slots)
        with self.assertRaisesRegex(ValueError, "interface"):
            source.replace(sp.TensorRule(source, sp.TensorExpression(E("1"))))


if __name__ == "__main__":
    unittest.main()
