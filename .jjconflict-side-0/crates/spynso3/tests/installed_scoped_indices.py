"""Scoped tensor copies preserve interfaces, contractions, and alphabet display."""

import unittest
import unicodedata
import xml.etree.ElementTree as ET

from symbolica import E, S
from symbolica.community.spenso import Representation, TensorExpression, TensorName


class ScopedIndicesTests(unittest.TestCase):
    def test_scopes_preserve_indices_and_contractions(self):
        rep = Representation.mink(4)
        p = TensorName("scope_test::p")(rep)
        q = TensorName("scope_test::q")(rep)
        index = E("gammalooprs::hedge(2,1)")
        bra, other = S("scope_test::bra", "scope_test::other")
        original = p(index)
        wrapped = original.wrap_indices(bra)
        self.assertIsInstance(wrapped, TensorExpression)
        self.assertEqual(wrapped.rank, 1)
        self.assertNotEqual(wrapped.structure.slots, original.structure.slots)
        self.assertEqual(
            TensorExpression(wrapped.to_expression()).structure.slots,
            wrapped.structure.slots,
        )
        self.assertEqual(wrapped.wrap_indices(bra), wrapped)
        self.assertEqual((original * wrapped).rank, 2)
        self.assertEqual((wrapped * q(index).wrap_indices(bra)).rank, 0)
        self.assertEqual((wrapped * q(index).wrap_indices(other)).rank, 2)
        self.assertEqual(
            (wrapped * q(index).wrap_indices(bra)).to_dots(),
            (p(index) * q(index)).to_dots(),
        )
        nested = wrapped.wrap_indices(other)
        self.assertEqual(
            TensorExpression(nested.to_expression()).structure.slots,
            nested.structure.slots,
        )
        zero = (0 * original).wrap_indices(bra)
        self.assertEqual(zero.rank, 1)
        self.assertEqual(zero.structure.slots, wrapped.structure.slots)

        for tensor in (original * wrapped, original * nested):
            html = tensor._repr_html_()
            assert isinstance(html, str)
            root = ET.fromstring(html[html.index("<math") : html.index("</math>") + 7])
            visible = unicodedata.normalize("NFKC", "".join(root.itertext()))
            self.assertIn("μ", visible)
            self.assertIn("′", visible)
            self.assertNotIn("hedge", visible)
            self.assertNotIn("index_scope", visible)
            self.assertNotIn("scope_test::bra", visible)


if __name__ == "__main__":
    unittest.main()
