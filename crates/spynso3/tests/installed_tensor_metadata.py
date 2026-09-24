"""Metadata boundaries and rich displays against the installed native extension."""

import re
import unittest

from symbolica import E, S
from symbolica.community.spenso import (
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
    def test_mixed_ports_shape_arguments_and_snapshot(self):
        d, x = S("metadata_tests::D", "metadata_tests::x")
        lorentz, euclidean = Representation.mink(d), Representation.euc(3)
        name = TensorName("metadata_tests::A")
        expression = name(x, 7, lorentz("mu"), euclidean)
        structure = expression.structure
        self.assertIsInstance(structure, TensorStructure)
        self.assertEqual(structure.rank, 2)
        self.assertEqual(structure.shape, (d, 3))
        self.assertEqual(structure.arguments, (x, E("7")))
        self.assertEqual(structure.name.to_expression(), name.to_expression())
        self.assertIsInstance(structure.slots[0], Slot)
        self.assertIsInstance(structure.slots[1], Representation)
        self.assertEqual(structure.slots[0].index, S("spynso3::mu"))
        self.assertEqual(
            TensorStructure(structure.slots, name=name, arguments=[x, 7]), structure
        )
        with self.assertRaises(AttributeError):
            structure.rank = 9
        indexed = expression("nu")
        self.assertIsInstance(indexed.structure.slots[1], Slot)
        self.assertIsInstance(structure.slots[1], Representation)

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
        expression = TensorExpression.gamma(4)("a", "b", "mu")
        library = TensorLibrary.hep_lib_atom()
        network = expression.to_network(library=library)
        source, metadata = network.expression(), network.structure
        network.execute(library=library)
        self.assertEqual(network.expression().to_expression(), source.to_expression())
        self.assertEqual(network.structure, metadata)
        self.assertEqual(
            network.result_tensor(library=library).structure.slots, metadata.slots
        )

    def test_scalar_unnamed_and_symbolic_dimension_displays(self):
        scalar = TensorExpression(0).structure
        self.assertEqual(scalar.shape, ())
        self.assertEqual(scalar.slots, ())
        self.assertIsNone(scalar.name)
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
        expression = TensorExpression.gamma(4)(1, 2, 1)
        original = expression.to_typst()
        first, second = expression.structure.to_html(), expression.structure.to_html()
        self.assertIn("TensorStructure", first)
        self.assertIn("TensorName", first)
        self.assertEqual(first.count('type="radio"'), 4)
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
