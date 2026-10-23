"""Metadata boundaries and rich displays against the installed native extension."""

import re
import unittest

from symbolica import E, S
from symbolica.community.tensor import (
    Representation,
    RepresentationName,
    Slot,
    Tensor,
    TensorExpression,
    TensorLibrary,
    TensorName,
    TensorStructure,
)


class TensorMetadataTests(unittest.TestCase):
    def test_signature_is_independent_of_identity_and_layout(self):
        d, x = S("metadata_tests::D", "metadata_tests::x")
        lorentz, euclidean = Representation.mink(d), Representation.euc(3)
        name = TensorName("metadata_tests::A")
        expression = name(x, 7, lorentz("mu"), euclidean)
        structure = expression.structure
        self.assertIsInstance(structure, TensorStructure)
        self.assertEqual(structure.rank, 2)
        self.assertCountEqual(structure.shape, (d, 3))
        self.assertEqual(expression.arguments, (x, E("7")))
        self.assertEqual(expression.name.to_expression(), name.to_expression())
        self.assertCountEqual(structure.axes, (lorentz("mu"), euclidean))
        self.assertEqual(TensorStructure(structure.axes[::-1]), structure)
        self.assertEqual(expression.with_name("metadata_tests::B").structure, structure)
        self.assertFalse(hasattr(structure, "name"))
        self.assertFalse(hasattr(structure, "arguments"))
        self.assertFalse(hasattr(structure, "permute_axes"))
        with self.assertRaises(AttributeError):
            structure.rank = 9
        indexed = expression("nu")
        self.assertEqual(indexed.permute_axes([1, 0]).structure, indexed.structure)
        self.assertTrue(all(isinstance(axis, Slot) for axis in indexed.structure.axes))
        self.assertIn(euclidean, structure.axes)

    def test_homogeneous_accessors_preserve_dimensions_and_duality(self):
        d = S("metadata_accessors::D")
        spaces = [
            Representation.mink(d),
            Representation.cof(3).dual(),
            Representation.euc(2),
        ]
        opened = TensorStructure(spaces)
        self.assertCountEqual(opened.representations(), spaces)
        self.assertEqual(opened.representations(), list(opened.axes))
        self.assertEqual(TensorStructure(spaces[::-1]), opened)
        with self.assertRaisesRegex(ValueError, "axis 0 is unresolved"):
            opened.slots()
        slots = [space(index) for space, index in zip(spaces, (9, 3, 7))]
        explicit = TensorStructure(slots)
        self.assertCountEqual(explicit.slots(), slots)
        self.assertEqual(explicit.slots(), list(explicit.axes))
        self.assertEqual(TensorStructure(slots[::-1]), explicit)
        with self.assertRaisesRegex(ValueError, "axis 0 has an explicit index"):
            explicit.representations()
        e = Representation.euc(3)
        self.assertEqual(TensorStructure([e(9), e(2)]).slots(), [e(2), e(9)])

    def test_mixed_accessors_fail_without_materializing_or_discarding_axes(self):
        space = Representation.euc(3)
        axes = (space(7), space)
        mixed = TensorStructure(axes)
        for structure in (mixed, TensorStructure(axes[::-1])):
            with self.assertRaisesRegex(ValueError, "axis 1 is unresolved"):
                structure.slots()
            with self.assertRaisesRegex(ValueError, "axis 0 has an explicit index"):
                structure.representations()
            self.assertEqual(structure.axes, axes)
            self.assertEqual(structure, mixed)
        scalar = TensorStructure([])
        self.assertEqual(scalar.slots(), [])
        self.assertEqual(scalar.representations(), [])

    def test_leaf_sum_and_product_share_canonical_free_slots(self):
        e = Representation.euc(3)
        A, B = TensorName("metadata_order::A"), TensorName("metadata_order::B")
        left, transposed = A(e(2), e(9)), A(e(9), e(2))
        expressions = [
            left,
            transposed,
            left + transposed,
            transposed + left,
            left - transposed,
            left - left,
            B(7, e(9), e(2)),
            left.permute_axes([1, 0]),
        ]
        for value in expressions:
            self.assertEqual(value.structure, left.structure)
            self.assertEqual(value.structure.slots(), [e(2), e(9)])
            if value:
                self.assertEqual(
                    TensorExpression(value.to_expression()).structure, left.structure
                )
        self.assertNotEqual(left.to_expression(), transposed.to_expression())
        self.assertEqual((A(e(9), e(4)) * B(e(4), e(2))).structure, left.structure)
        self.assertEqual((B(e(4), e(2)) * A(e(9), e(4))).structure, left.structure)

    def test_unresolved_signature_retains_multiplicity_without_occurrence_order(self):
        e, m = Representation.euc(3), Representation.mink(4)
        A = TensorName("metadata_order::Open")
        first, other = A(m, e, m), A(e, m, m)
        self.assertEqual(first.structure, other.structure)
        self.assertEqual(first.axes, (m, e, m))
        self.assertEqual(other.axes, (e, m, m))
        self.assertCountEqual(first.structure.representations(), [e, m, m])
        # Interface observation must not change positional compatibility for operations.
        with self.assertRaisesRegex(ValueError, "interface"):
            first + other
        with self.assertRaisesRegex(ValueError, "interface"):
            A(e, m) + A(m, e)

    def test_supplied_signature_checks_axes_without_changing_argument_order(self):
        e, m = Representation.euc(3), Representation.mink(4)
        A = TensorName("metadata_order::Checked")
        original = A(m, e, m)
        signature = TensorStructure([e, m, m])
        checked = TensorExpression(original.to_expression(), structure=signature)
        self.assertEqual(checked.structure, signature)
        self.assertEqual(checked.axes, original.axes)
        self.assertEqual(checked.to_expression(), original.to_expression())
        zero = TensorExpression(0, structure=signature)
        self.assertEqual(zero.structure, signature)
        with self.assertRaisesRegex(ValueError, "interface"):
            TensorExpression(
                original.to_expression(), structure=TensorStructure([e, m])
            )

    def test_network_preserves_argument_permutations_under_a_canonical_signature(self):
        r = Representation.euc(2)
        A = TensorName("metadata_order::Matrix")(r, r)
        library = TensorLibrary()
        library.register(Tensor.dense(A, [1.0, 2.0, 3.0, 4.0]))
        i, j = r("i"), r("j")
        left, right = A(i, j), A(j, i)
        select_i = Tensor.dense(TensorName("metadata_order::select_i")(i), [1.0, 0.0])
        select_j = Tensor.dense(TensorName("metadata_order::select_j")(j), [0.0, 1.0])
        for value, expected in [
            (left, 2),
            (right, 3),
            (left + right, 5),
            (right + left, 5),
            (left - right, -1),
        ]:
            with self.subTest(expression=str(value)):
                self.assertEqual(value.structure, left.structure)
                network = value.to_network(library=library)
                self.assertEqual(network.structure, left.structure)
                network.execute(library=library)
                self.assertEqual(
                    network.result_tensor(library=library).structure, left.structure
                )
                selected = value.to_network(library=library) * select_i * select_j
                selected.execute(library=library)
                self.assertEqual(selected.result_tensor(library=library)[0], expected)

    def test_expression_descriptor_and_data_are_separate(self):
        m, e = Representation.mink(2), Representation.euc(3)
        expression = TensorName("metadata_tests::M")(m("mu"), e("i"))
        tensor = Tensor.dense(expression, [float(i) for i in range(6)])
        self.assertEqual(tensor.structure, expression.structure)
        self.assertEqual(
            tensor.expression().to_expression(), expression.to_expression()
        )
        tensor[0] = 17.0
        self.assertEqual(tensor.structure, expression.structure)
        self.assertEqual(
            tensor.expression().to_expression(), expression.to_expression()
        )
        self.assertEqual(tensor[0], 17.0)

    def test_network_execution_preserves_source_and_structure(self):
        expression = TensorExpression.dirac_gamma(4)("a", "b", "mu")
        library = TensorLibrary.hep_lib_atom()
        network = expression.to_network(library=library)
        source, metadata = network.expression(), network.structure
        network.execute(library=library)
        self.assertEqual(network.expression().to_expression(), source.to_expression())
        self.assertEqual(network.structure, metadata)
        self.assertEqual(
            network.result_tensor(library=library).structure.axes, metadata.axes
        )

    def test_scalar_unnamed_and_symbolic_dimension_displays(self):
        scalar = TensorExpression(0).structure
        self.assertEqual(scalar.shape, ())
        self.assertEqual(scalar.axes, ())
        self.assertIn("no external slots", scalar.to_html())
        rep = Representation.mink(S("metadata_tests::D"))
        self.assertIn("Dimension", rep.to_html())
        self.assertIn("(+, −, …)", rep.to_html())
        self.assertEqual(TensorStructure([rep]).shape, (S("metadata_tests::D"),))

    def test_representation_name_duality_and_metric(self):
        rep = Representation.mink(4)
        self.assertIsInstance(rep.name, RepresentationName)
        self.assertEqual(rep.name.duality, "self-dual")
        self.assertEqual([rep.name.metric_sign(i) for i in range(4)], [1, -1, -1, -1])
        self.assertEqual(rep.name.dual(), rep.name)
        self.assertIn("(+, −, −, −)", rep.to_html())
        self.assertNotIn("Inspect", rep.to_html())
        color = Representation.cof(3)
        self.assertEqual(color.name.duality, "base")
        self.assertEqual(color.dual().name.duality, "dual")
        self.assertEqual(color.name.dual(), color.dual().name)
        self.assertIn("dual pairing", color.to_html())
        self.assertIn("(+, +, +)", Representation.euc(3).to_html())

    def test_structure_display_uses_native_controls_and_distinct_groups(self):
        expression = TensorExpression.dirac_gamma(4)(1, 2, 1)
        original = expression.to_typst()
        first, second = expression.structure.to_html(), expression.structure.to_html()
        self.assertIn("TensorStructure", first)
        self.assertNotIn("TensorName", first)
        self.assertEqual(first.count('type="radio"'), 3)
        self.assertIn("Slot 0", first)
        self.assertNotIn("<script", first)
        group = r'name="(spenso-structure-\d+)"'
        self.assertNotEqual(re.search(group, first)[1], re.search(group, second)[1])
        self.assertEqual(expression.to_typst(), original)
        self.assertNotIn("data-spenso-metadata", expression.to_html())
        slot = Representation.mink(4)(1)
        self.assertIn("Index", slot.to_html())
        self.assertEqual(slot.index, 1)


if __name__ == "__main__":
    unittest.main()
