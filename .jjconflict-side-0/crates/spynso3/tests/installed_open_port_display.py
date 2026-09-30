"""Unresolved axes remain distinct from assigned and contracted index labels."""

import re
import unicodedata
import unittest
import xml.etree.ElementTree as ET

import typst
from symbolica import S
from symbolica.community.tensor import (
    AUTO,
    DisplaySettings,
    Representation,
    TensorExpression,
    TensorName,
    dot,
)


def mathml(document):
    return ET.fromstring(re.search(r"<math\b.*?</math>", document, re.DOTALL)[0])


def visible_nodes(node):
    if node.tag not in {"mphantom", "annotation", "annotation-xml"}:
        yield node
        for child in node:
            yield from visible_nodes(child)


def visible(node):
    text = "".join(child.text or "" for child in visible_nodes(node))
    return "".join(unicodedata.normalize("NFKC", text).split())


class OpenPortDisplayTests(unittest.TestCase):
    def assert_renderers_agree(self, expression, expected, settings=None):
        settings = settings or DisplaySettings()
        original, structure = expression.to_expression(), expression.structure
        source = expression.to_typst(settings=settings)
        native = typst.compile(f"$ {source} $".encode(), format="html").decode()
        rich = expression.to_html(settings=settings)
        self.assertEqual(visible(mathml(native)), expected)
        self.assertEqual(visible(mathml(rich)), expected)
        self.assertNotIn("open_index", source)
        self.assertNotIn("open_index", expression.to_latex(settings=settings))
        self.assertNotIn("open_index", expression.format_tensor())
        self.assertEqual(expression.to_expression(), original)
        self.assertEqual(expression.structure, structure)
        return list(visible_nodes(mathml(rich)))

    def test_factories_and_user_tensors_share_the_same_placeholder_convention(self):
        rep = Representation.euc(3)
        a = TensorName("open_display::A")(rep, rep)
        gamma = TensorExpression.dirac_gamma(4)
        for style in ("alphabet", "graph", "raw"):
            for expression, expected, axes in (
                (a, "A□0□1", ["0", "1"]),
                (a(AUTO, "j"), "A□j", ["0"]),
                (gamma, "γ□0□1□2", ["0", "1", "2"]),
                (gamma("a", "b", AUTO), "γab□", ["2"]),
            ):
                with self.subTest(style=style, expected=expected):
                    nodes = self.assert_renderers_agree(
                        expression, expected, DisplaySettings(index_style=style)
                    )
                    tips = [n for n in nodes if "data-spenso-open-axis" in n.attrib]
                    self.assertCountEqual(
                        [n.attrib["data-spenso-open-axis"] for n in tips], axes
                    )
                    self.assertTrue(all("spenso::" in n.attrib["title"] for n in tips))
        self.assertEqual(a(AUTO, "j")("i"), a("i", "j"))
        self.assertIn("<svg", gamma.to_svg())

    def test_axis_numbers_exclude_scalar_arguments_and_preserve_dual_rows(self):
        rep = Representation("open_display::V", 3, is_self_dual=False)
        a = TensorName("open_display::B")(S("x"), 7, rep.dual(), rep)
        nodes = self.assert_renderers_agree(a, "Bx7□0□1")
        tips = [n for n in nodes if "data-spenso-open-axis" in n.attrib]
        self.assertEqual(len(tips), 2)
        self.assertTrue(any("dual" in n.attrib["title"] for n in tips))
        self.assertEqual({n.attrib["data-spenso-open-axis"] for n in tips}, {"0", "1"})
        # Qualified display keeps the hollow square and adds representation data.
        settings = DisplaySettings(show_dimensions=True)
        for output in (
            a.to_html(settings=settings),
            typst.compile(
                f"$ {a.to_typst(settings=settings)} $".encode(), format="html"
            ).decode(),
        ):
            self.assertEqual(visible(mathml(output)).count("□"), 2)

    def test_compact_contractions_and_generated_alphabets_are_preserved(self):
        gamma = TensorExpression.dirac_gamma(4)
        p = TensorName.vector("open_display::p")(Representation.mink(4))
        self.assert_renderers_agree(p, "p□")
        self.assert_renderers_agree(dot(p, p), "p⋅p")
        self.assert_renderers_agree(gamma * gamma, "[γ□γ□]□0□1")
        indexed = (gamma * gamma).index("a", "c", "μ", "ν").undo_chain()
        self.assert_renderers_agree(indexed, "γabμγbcν")
        self.assertNotIn("□", indexed.format_tensor())
        self.assertEqual(indexed.rank, 4)


if __name__ == "__main__":
    unittest.main()
