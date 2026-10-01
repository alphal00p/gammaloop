"""Expression inheritance and tensor overrides, with the built extension."""

import copy
import operator
import pickle
import unittest

from symbolica import E, Expression, FunctionDefinition, S
from symbolica.community import tensor as sp


class TypedSurfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rep = sp.Representation.euc(3)
        cls.P = sp.TensorName("typed_surface::P")
        cls.Q = sp.TensorName("typed_surface::Q")
        cls.x, cls.y = S("typed_surface::x", "typed_surface::y")

    def test_expression_subclass_and_checked_reconstruction(self):
        source = self.P(self.x, self.rep("j"), self.rep("i")).permute_axes([1, 0])
        self.assertIsInstance(source, Expression)
        self.assertTrue(issubclass(sp.TensorExpression, Expression))
        self.assertIn(Expression, sp.TensorExpression.__mro__)
        raw = source.to_expression().replace(self.x, self.y)
        changed = sp.TensorExpression(raw, structure=source.structure)
        self.assertEqual(changed.to_expression(), raw)
        self.assertEqual(changed.structure.axes, source.structure.axes)
        self.assertEqual(changed.arguments, (self.y,))
        self.assertEqual(str(changed.name), str(source.name))
        zero = sp.TensorExpression(E("0"), structure=source.structure)
        self.assertEqual(zero.rank, source.rank)
        self.assertEqual(zero.structure.axes, source.structure.axes)
        with self.assertRaisesRegex(ValueError, "interface"):
            sp.TensorExpression(E("1"), structure=source.structure)

    def test_inherited_symbolica_api_uses_the_actual_expression(self):
        source = self.x * self.P(self.rep("i")) + self.y * self.Q(self.rep("i"))
        raw = source.to_expression()
        # These methods come from Expression itself, without forwarding wrappers.
        for name in ("coefficient", "terms", "match", "get_tags", "to_polynomial"):
            with self.subTest(method=name):
                self.assertIs(
                    getattr(sp.TensorExpression, name), getattr(Expression, name)
                )
        coefficient = source.coefficient(self.x)
        self.assertIs(type(coefficient), Expression)
        self.assertEqual(coefficient, self.P(self.rep("i")).to_expression())
        self.assertEqual(list(source.terms()), list(raw.terms()))
        pattern = S("typed_surface::term_")
        self.assertEqual(list(source.match(pattern)), list(raw.match(pattern)))
        self.assertEqual(
            self.P(self.rep("i")).get_tags(), raw.coefficient(self.x).get_tags()
        )
        function = S("typed_surface::f")
        self.assertEqual(function(source), function(raw))

    def test_base_payload_stays_synchronized_after_tensor_operations(self):
        source = self.P(self.rep("j"), self.rep("i")).permute_axes([1, 0])
        metric = sp.TensorExpression.g(self.rep)("i", "k")
        values = (
            copy.copy(source),
            source.with_name("typed_surface::Renamed"),
            source.permute_axes([1, 0]),
            source.rename_indices({"i": "k"}),
            self.P(self.rep)("i"),
            (self.x * source + self.y * source).collect_factors(),
            ((self.x**2 - 1) * source).factor(),
            (metric * source).contract(),
            (metric * source).simplify_algebra(),
            0 * source,
        )
        for value in values:
            with self.subTest(expression=value.to_expression()):
                raw = value.to_expression()
                # Bypass tensor equality to check the inherited native payload itself.
                self.assertTrue(Expression.__eq__(value, raw))
                self.assertEqual(
                    Expression.get_all_symbols(value), raw.get_all_symbols()
                )
                self.assertEqual(list(value.terms()), list(raw.terms()))
                self.assertEqual(
                    S("typed_surface::payload")(value), S("typed_surface::payload")(raw)
                )

    def test_tensor_equality_hash_and_truth_preserve_logical_interface(self):
        source = self.P(self.rep("j"), self.rep("i"))
        same = copy.copy(source).with_name("typed_surface::PresentationOnly")
        self.assertEqual(source, same)
        self.assertEqual(hash(source), hash(same))
        self.assertEqual({source: "value"}[same], "value")
        self.assertNotEqual(source, source.permute_axes([1, 0]))
        self.assertNotEqual(source, source.to_expression())
        self.assertNotEqual(source.to_expression(), source)
        self.assertNotEqual(sp.TensorExpression(E("0")), 0 * source)
        self.assertNotEqual(0 * source, E("0"))
        self.assertNotEqual(E("0"), 0 * source)
        self.assertTrue(source)
        self.assertTrue(sp.TensorExpression(self.x))
        self.assertFalse(0 * source)
        self.assertFalse(sp.TensorExpression(E("0")))
        self.assertEqual(sp.TensorExpression(source), source)

    def test_deepcopy_preserves_named_tensor_metadata_and_pickle_is_explicit(self):
        source = self.P(7, self.rep("j"), self.rep("i")).permute_axes([1, 0])
        for value in (source, 0 * source):
            value = value.with_name(sp.TensorName("typed_surface::CopyMetadata")(7))
            copied = copy.deepcopy(value)
            self.assertIsNot(copied, value)
            self.assertIsInstance(copied, sp.TensorExpression)
            self.assertEqual(copied, value)
            self.assertEqual(copied.structure, value.structure)
            self.assertEqual(copied.arguments, (E("7"),))
            self.assertEqual(copied.structure.axes, source.structure.axes)
            self.assertTrue(Expression.__eq__(copied, value.to_expression()))
            # Inherited Expression pickling must not silently drop tensor metadata.
            with self.assertRaises(TypeError):
                pickle.dumps(value)

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
        self.assertEqual(quotient.structure.axes, vector.structure.axes)
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
        self.assertEqual(result.structure.axes, numerator.structure.axes)
        self.assertEqual(
            (numerator.to_expression() / denominator).structure.axes,
            numerator.structure.axes,
        )
        self.assertEqual(
            (0 * numerator / denominator).structure.axes, numerator.structure.axes
        )
        self.assertFalse(0 * numerator / denominator)

    def test_selected_port_contraction_retains_declared_logical_order(self):
        left = self.P(self.rep("i"), self.rep("j")).permute_axes([1, 0])
        right = self.Q(self.rep("k"), self.rep("i")).permute_axes([1, 0])
        result = left.contract_ports(right, left=1, right=0)
        self.assertEqual(result.axes, (left.axes[0], right.axes[1]))
        self.assertEqual(result.rank, 2)
        with self.assertRaises(TypeError):
            left.contract(right, left=1, right=0)

    def test_network_binary_contraction_uses_contract_ports(self):
        left = self.P(self.rep("i"), self.rep("j")).to_network()
        right = self.Q(self.rep("k"), self.rep("i")).to_network()
        result = left.contract_ports(right, left=0, right=1)
        self.assertEqual(result.rank, 2)
        self.assertFalse(hasattr(left, "contract"))

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
        self.assertEqual(((0 * opened) ** 3).structure.axes, opened.structure.axes)

    def test_explicit_evaluator_uses_symbolic_payloads(self):
        source = sp.TensorExpression((self.x + 1) ** 3 + self.y)
        points = [[2.0, 5.0], [-1.0, 7.0]]
        for value in (source, source.simplify_algebra()):
            evaluator = value.evaluator([self.x, self.y], iterations=1, n_cores=1)
            self.assertEqual(evaluator.evaluate(points).tolist(), [[32.0], [7.0]])
        self.assertEqual(source.to_expression(), (self.x + 1) ** 3 + self.y)

    def test_inherited_evaluator_accepts_function_definitions_and_jit_settings(self):
        function = S("typed_surface::evaluate_f")
        source = sp.TensorExpression(function(self.x) + 3 * self.x)
        definition = FunctionDefinition(function, [self.y], self.y**2 + 1)
        evaluator = source.evaluator(
            [self.x], functions=[definition], jit_compile=False, n_cores=1
        )
        self.assertEqual(evaluator.evaluate([[2.0], [-1.0]]).tolist(), [[11.0], [-1.0]])

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
        self.assertEqual(zero.structure.axes, source.structure.axes)
        with self.assertRaisesRegex(ValueError, "interface"):
            source.replace(sp.TensorRule(source, sp.TensorExpression(E("1"))))


if __name__ == "__main__":
    unittest.main()
