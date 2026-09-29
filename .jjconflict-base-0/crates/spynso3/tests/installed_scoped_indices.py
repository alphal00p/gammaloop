"""Scoped tensor copies preserve interfaces, contractions, and alphabet display."""

import unicodedata
import unittest
import xml.etree.ElementTree as ET

from symbolica import E, S
from symbolica.community.spenso import Representation, TensorExpression, TensorName


class ScopedIndicesTests(unittest.TestCase):
    def test_dummy_only_scopes_keep_free_and_unresolved_ports(self):
        rep = Representation.euc(3)
        metric = TensorExpression.g(rep)
        p, q = TensorName.vector("scope_dummy::p"), TensorName.vector("scope_dummy::q")
        lhs, rhs = S("scope_dummy::lhs", "scope_dummy::rhs")
        left = metric("a", "b") * p(rep("a"))
        right = metric("a", "c") * q(rep("a"))
        scoped = left.wrap_indices(lhs, dummies_only=True)
        self.assertEqual(scoped.structure.slots, left.structure.slots)
        self.assertEqual(scoped.wrap_indices(lhs, dummies_only=True), scoped)
        result = (
            (scoped * right.wrap_indices(rhs, dummies_only=True))
            .contract(rank_one=False)
            .to_expression()
        )
        self.assertEqual(result, p(rep("b")) * q(rep("c")))
        for dummies_only in (False, True):
            self.assertEqual(
                metric.wrap_indices(lhs, dummies_only=dummies_only), metric
            )

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
