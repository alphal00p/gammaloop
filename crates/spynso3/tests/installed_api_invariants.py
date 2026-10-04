"""Color invariant notation agrees across source and notebook renderers."""

import re
import unicodedata
import unittest
import xml.etree.ElementTree as ET

import typst
from symbolica import E, S
from symbolica.community.tensor import (
    DisplaySettings,
    Representation,
    TensorExpression,
    to_html,
    to_svg,
    to_typst,
)


def mathml(document):
    return ET.fromstring(re.search(r"<math\b.*?</math>", document, re.DOTALL).group())


def visible(node):
    if node.tag in {"mphantom", "annotation", "annotation-xml"}:
        return ""
    text = (node.text or "") + "".join(
        visible(child) + (child.tail or "") for child in node
    )
    return "".join(unicodedata.normalize("NFKC", text).split())


class InvariantDisplayTests(unittest.TestCase):
    def setUp(self):
        self.fundamental = Representation.cof(3)
        self.adjoint = Representation.coad(8)

    def assert_renderers_agree(self, expression, expected, settings):
        source = to_typst(expression, settings=settings)
        document = typst.compile(f"$ {source} $".encode(), format="html").decode()
        rich = to_html(expression, settings=settings)
        self.assertEqual(visible(mathml(document)), expected)
        self.assertEqual(visible(mathml(rich)), expected)
        tensor = TensorExpression(expression)
        self.assertEqual(tensor.to_typst(settings=settings), source)
        self.assertEqual(tensor.to_expression(), expression)
        return mathml(rich)

    def test_compact_and_explicit_invariants_agree_across_renderers(self):
        f, a = self.fundamental, self.adjoint
        for typed, name, compact, explicit in (
            (f.casimir(), "CF", "CF", "C2(F)"),
            (a.casimir(), "CA", "CA", "C2(A)"),
            (f.dynkin_index(), "TR", "TR", "I2(F)"),
        ):
            for style, expected in (("compact", compact), ("explicit", explicit)):
                settings = DisplaySettings(invariant_style=style)
                with self.subTest(name=name, style=style):
                    root = self.assert_renderers_agree(typed, expected, settings)
                    self.assertTrue(list(root.iter("msub")))
                    self.assertEqual(
                        TensorExpression(typed).to_latex(settings=settings),
                        "$$" + to_typst(typed, settings=settings) + "$$",
                    )
        self.assert_renderers_agree(S("spenso::Nc"), "Nc", DisplaySettings())
        self.assertEqual(DisplaySettings().invariant_style, "compact")
        with self.assertRaisesRegex(ValueError, "invariant_style"):
            DisplaySettings(invariant_style="unknown")

    def test_higher_degrees_and_representation_dimensions(self):
        f, a = self.fundamental, self.adjoint
        for expression, compact, dimensions in (
            (f.casimir(4), "C4(F)", "C4(F3)"),
            (a.dynkin_index(), "I2(A)", "I2(A8)"),
            (f.gram(4, a), "G4(F,A)", "G4(F3,A8)"),
            (f.casimir(), "CF", "C2(F3)"),
            (Representation.cof(5).casimir(), "CF", "C2(F5)"),
            (f.casimir(S("k") + 1), "C1+k(F)", "C1+k(F3)"),
        ):
            for show_dimensions, expected in ((False, compact), (True, dimensions)):
                with self.subTest(
                    expression=str(expression), dimensions=show_dimensions
                ):
                    self.assert_renderers_agree(
                        expression,
                        expected,
                        DisplaySettings(show_dimensions=show_dimensions),
                    )
        # A dual representation stays visibly distinct, including under scripts.
        root = self.assert_renderers_agree(
            f.dual().casimir(4), "C4(F̄)", DisplaySettings()
        )
        self.assertTrue(list(root.iter("mover")))

    def test_scalar_delimiters_are_independent_of_tensor_layout(self):
        expression = self.fundamental.gram(4, self.adjoint)
        for layout in ("ports", "schoonschip", "call"):
            settings = DisplaySettings(
                tensor_layout=layout, parentheses=False, commas=False
            )
            self.assertEqual(
                visible(mathml(to_html(expression, settings=settings))), "G4(F,A)"
            )
        for explicit in (False, True):
            settings = DisplaySettings(
                symbol_scripts=False,
                invariant_style="explicit" if explicit else "compact",
            )
            self.assert_renderers_agree(expression, "gram(4,F,A)", settings)
            self.assert_renderers_agree(
                self.fundamental.casimir(), "cas(2,F)" if explicit else "CF", settings
            )

    def test_tensor_spectators_powers_and_exact_identity_survive_display(self):
        coefficient = self.fundamental.casimir() ** 2
        gamma = TensorExpression.dirac_gamma(4)(1, 2, 1)
        tensor = coefficient * gamma
        before = tensor.to_expression()
        slots = tensor.structure.axes
        for style in ("compact", "explicit"):
            settings = DisplaySettings(invariant_style=style)
            rendered = tensor.to_html(settings=settings)
            self.assertIsNotNone(rendered)
            root = mathml(to_html(coefficient, settings=settings))
            self.assertTrue(list(root.iter("msup")) or list(root.iter("msubsup")))
            ET.fromstring(to_svg(coefficient, settings=settings))
            self.assertEqual(tensor.to_expression(), before)
            self.assertEqual(tensor.structure.axes, slots)

    def test_unrelated_functions_are_not_color_representations(self):
        expression = E("spenso::cas(2,invariant_display_tests::cof(3))")
        self.assert_renderers_agree(expression, "C2(cof(3))", DisplaySettings())
        expression = E(
            "invariant_display_tests::cas(2,invariant_display_tests::cof(3))"
        )
        self.assert_renderers_agree(expression, "cas(2,cof(3))", DisplaySettings())


if __name__ == "__main__":
    unittest.main()
