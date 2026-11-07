"""Tensor subclasses remain scalar-checked at component storage boundaries."""

import unittest

from symbolica import E, Expression, S
from symbolica.community import tensor as sp


class TensorSubclassComponentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.descriptor = sp.TensorName("subclass_components::v")(
            sp.Representation.euc(2)
        )
        cls.x = S("subclass_components::x")

    def test_dense_rejects_rankful_values_including_numeric_typed_zero(self):
        for component in (self.descriptor, 0 * self.descriptor):
            with (
                self.subTest(component=component),
                self.assertRaisesRegex(TypeError, "components must be scalar"),
            ):
                sp.Tensor.dense(self.descriptor, [component, component])

    def test_dense_accepts_scalar_tensor_and_existing_raw_expression(self):
        tensor = sp.Tensor.dense(self.descriptor, [sp.TensorExpression(self.x), E("2")])
        self.assertEqual(tensor[:], [self.x, E("2")])

    def test_dense_and_sparse_assignment_preserve_component_scope(self):
        for tensor in (
            sp.Tensor.dense(self.descriptor, [self.x, E("0")]),
            sp.Tensor.sparse(self.descriptor, Expression),
        ):
            for component in (self.descriptor, 0 * self.descriptor):
                with (
                    self.subTest(storage=tensor.storage, component=component),
                    self.assertRaisesRegex(Exception, "components must be scalar"),
                ):
                    tensor[0] = component
            tensor[0] = sp.TensorExpression(self.x)
            tensor[1] = E("2")
            self.assertEqual(tensor[:], [self.x, E("2")])

    def test_map_rejects_rankful_callback_results_before_dtype_coercion(self):
        for tensor in (
            sp.Tensor.dense(self.descriptor, [self.x, E("0")]),
            sp.Tensor.sparse(self.descriptor, Expression),
        ):
            for dtype in (float, complex, Expression):
                for component in (self.descriptor, 0 * self.descriptor):
                    with (
                        self.subTest(dtype=dtype, component=component),
                        self.assertRaisesRegex(TypeError, "components must be scalar"),
                    ):
                        tensor.map_components(
                            lambda _, component=component: component, dtype=dtype
                        )
            result = tensor.map_components(
                lambda _: sp.TensorExpression(self.x), dtype=Expression
            )
            self.assertEqual(result[:], [self.x, self.x])


if __name__ == "__main__":
    unittest.main()
