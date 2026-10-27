"""Bounded notebook previews, complete exports, and lazy backend caching."""

import re
import unittest
import xml.etree.ElementTree as ET

from symbolica import E, FormattedOutput, PrintMode
from symbolica.community import tensor as sp


class SymbolicDisplayTests(unittest.TestCase):
    def test_complete_arithmetic_uses_symbolica_layout(self):
        rep = sp.Representation.euc(3)
        cases = [
            (sp.TensorExpression(E(source)), expected)
            for source, expected in [
                ("a-b", ["a", "−", "b"]),
                ("x-3", ["−", "3", "+", "x"]),
                ("-x", ["−", "x"]),
                ("a-2*b", ["a", "−", "2", "b"]),
                ("x/y", ["x", "y"]),
            ]
        ]
        cases.append(
            (
                sp.TensorName("arithmetic_tests::A")(rep("i"))
                - sp.TensorName("arithmetic_tests::B")(rep("i")),
                ["A", "i", "−", "B", "i"],
            )
        )
        for value, expected in cases:
            before = value.to_expression()
            for render in (
                value.to_html,
                value._repr_html_,
                value.formatted()._repr_html_,
            ):
                html = render()
                root = ET.fromstring(
                    re.search(r"<math\b[^>]*>.*?</math>", html, re.DOTALL).group()
                )
                tokens = [
                    node.text
                    for node in root.iter()
                    if node.tag.rsplit("}", 1)[-1] in {"mi", "mn", "mo", "mtext"}
                    and node.text != "\u2062"
                ]
                self.assertEqual(tokens, expected, html)
                if expected == ["x", "y"]:
                    fraction = root.find(".//{*}mfrac")
                    self.assertIsNotNone(fraction)
                    self.assertEqual(
                        ["".join(child.itertext()) for child in fraction], ["x", "y"]
                    )
                self.assertNotIn("data-symbolica-marker", html)
                self.assertNotIn("Expression preview", html)
            self.assertEqual(value.to_expression(), before)

    def test_complete_export_keeps_a_declined_semantic_leaf(self):
        # The semantic fragment plus Symbolica's path wrapper exceeds 4 KiB,
        # although the source and its ordinary native call fit the viewer.
        name = "T" * 900
        labels = [f"i{index}" for index in range(100)]
        rep = sp.Representation.euc(3)
        value = sp.TensorName(f"declined_leaf_tests::{name}")(
            *(rep(label) for label in labels)
        )
        before = value.to_expression()
        for render in (value.to_html, value._repr_html_):
            html = render()
            root = ET.fromstring(
                re.search(r"<math\b[^>]*>.*?</math>", html, re.DOTALL).group()
            )
            self.assertIsNotNone(root.find(".//{*}mmultiscripts"))
            tokens = [
                node.text
                for node in root.iter()
                if node.tag.rsplit("}", 1)[-1] in {"mi", "mn", "mo"}
            ]
            self.assertEqual(tokens, [name, *labels])
            self.assertNotIn("spenso-complete-leaf", html)
            self.assertNotIn("data-symbolica-marker", html)
        self.assertEqual(value.to_expression(), before)

    def test_complete_arithmetic_fallback_retains_every_term(self):
        value = sp.TensorExpression(
            E("+".join(f"fallback_term_{i}" for i in range(600)))
        )
        html = value.to_html()
        root = ET.fromstring(
            re.search(r"<math\b[^>]*>.*?</math>", html, re.DOTALL).group()
        )
        tokens = [
            node.text for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "mi"
        ]
        self.assertCountEqual(tokens, [f"fallback_term_{i}" for i in range(600)])
        self.assertNotIn("data-symbolica-marker", html)
        self.assertNotIn("⋯", html)

    def assertBoundedNotebookPreview(self, value, *, full_text=False):
        html = value._repr_html_()
        self.assertIn("<math", html)
        self.assertIn("Expression preview", html)
        self.assertIn("data-symbolica-marker", html)
        self.assertLessEqual(len(html.encode()), 256 * 1024)

        output = value.formatted()
        self.assertEqual(output._repr_html_(), html)
        self.assertEqual(output._repr_html_(), html)
        if full_text:
            self.assertEqual(output._repr_latex_(), value.to_latex())
            self.assertEqual(value._repr_latex_(), value.to_latex())
            plain = value.format_tensor()
            self.assertEqual(output.format_plain(), plain)
        else:
            self.assertIsNone(output._repr_latex_())
            self.assertIsNone(output._repr_latex_())
            self.assertEqual(
                value._repr_latex_(),
                r"$$\text{Large tensor expression: use the paged viewer}$$",
            )
            plain = output.format_plain()
            self.assertEqual(
                plain,
                "Tensor expression (large; use the paged notebook viewer or explicit export)",
            )
        for rendered in (str(value), repr(value), str(output), repr(output)):
            self.assertEqual(rendered, plain)
        for rendered in (value, output):

            class Pretty:
                def __init__(self):
                    self.parts = []

                def text(self, text):
                    self.parts.append(text)

            pretty = Pretty()
            rendered._repr_pretty_(pretty, False)
            self.assertEqual(pretty.parts, [plain])
        return html

    def test_formatted_only_computes_each_requested_backend_once(self):
        modes = []

        def printer(expression, *, mode, **options):
            modes.append(mode)
            return {
                PrintMode.Symbolica: "Jbar",
                PrintMode.Latex: r"\bar{J}",
                PrintMode.Typst: "macron(J)",
                PrintMode.MathML: "<mover><mi>J</mi><mo>¯</mo></mover>",
            }.get(mode)

        value = sp.TensorName("full_display_tests::Current", print=printer)(
            sp.Representation.euc(2)
        )
        original = value.to_expression()
        modes.clear()
        formatted = value.formatted()
        self.assertIs(type(formatted), FormattedOutput)
        self.assertEqual(modes, [])

        plain = formatted.format_plain()
        self.assertIn("Jbar", plain)
        self.assertEqual(set(modes), {PrintMode.Symbolica})
        count = len(modes)
        self.assertEqual(str(formatted), plain)
        self.assertEqual(repr(formatted), plain)
        self.assertEqual(len(modes), count)

        modes.clear()
        latex = formatted._repr_latex_()
        self.assertIn(r"\bar{J}", latex)
        self.assertEqual(set(modes), {PrintMode.Latex})
        count = len(modes)
        self.assertEqual(formatted._repr_latex_(), latex)
        self.assertEqual(len(modes), count)

        modes.clear()
        html = formatted._repr_html_()
        self.assertIn("<mover><mi>J</mi><mo>¯</mo></mover>", html)
        self.assertEqual(set(modes), {PrintMode.MathML})
        count = len(modes)
        self.assertEqual(formatted._repr_html_(), html)
        self.assertEqual(len(modes), count)
        self.assertEqual(value.to_expression(), original)

    def test_free_formatted_owns_its_input_until_rendered(self):
        expression = E("(full_display_tests::x+1)^3")
        output = sp.formatted(expression)
        del expression
        self.assertIn("3", output._repr_latex_())
        self.assertIn("<msup", output._repr_html_())

    def test_declined_print_callback_caches_fallback_without_preventing_latex(self):
        for raises in (False, True):
            with self.subTest(raises=raises):
                modes = []

                def printer(expression, *, mode, modes=modes, raises=raises, **options):
                    modes.append(mode)
                    if mode == PrintMode.MathML:
                        if raises:
                            # Registered printers log Python errors and decline.
                            raise ValueError("intentional-print-decline")
                        return None
                    return r"\bar{J}" if mode == PrintMode.Latex else "J"

                name = "Raised" if raises else "Declined"
                value = sp.TensorName(f"full_display_tests::{name}", print=printer)(
                    7, sp.Representation.euc(2)("a")
                )
                original = value.to_expression()
                modes.clear()
                output = value.formatted()
                self.assertEqual(modes, [])
                html = output._repr_html_()
                self.assertIn("<math", html)
                self.assertIn(name, html)
                self.assertIn("<mn>7</mn>", html)
                self.assertEqual(set(modes), {PrintMode.MathML})
                count = len(modes)
                self.assertEqual(output._repr_html_(), html)
                self.assertEqual(len(modes), count)

                modes.clear()
                latex = output._repr_latex_()
                self.assertIn(r"\bar{J}", latex)
                self.assertEqual(set(modes), {PrintMode.Latex})
                count = len(modes)
                self.assertEqual(output._repr_latex_(), latex)
                self.assertEqual(len(modes), count)
                self.assertEqual(value.to_expression(), original)

    def test_export_limit_errors_propagate_and_retry_independently_of_latex(self):
        modes = []
        oversized = True

        def printer(expression, *, mode, **options):
            modes.append(mode)
            if mode == PrintMode.MathML:
                if oversized:
                    # Exceed the complete export's 16 MiB output budget.
                    return "<mtext>" + "x" * (16 * 1024 * 1024) + "</mtext>"
                return "<mover><mi>J</mi><mo>¯</mo></mover>"
            return r"\bar{J}" if mode == PrintMode.Latex else "J"

        value = sp.TensorName("full_display_tests::ExportLimit", print=printer)(
            sp.Representation.euc(2)
        )
        original = value.to_expression()
        modes.clear()
        output = value.formatted()
        self.assertEqual(modes, [])
        error = "complete MathML export failed: .*rendering limits"
        with self.assertRaisesRegex(ValueError, error):
            value.to_html()
        self.assertEqual(set(modes), {PrintMode.MathML})
        modes.clear()
        for _ in range(2):
            count = len(modes)
            with self.assertRaisesRegex(ValueError, error):
                output._repr_html_()
            self.assertEqual(set(modes), {PrintMode.MathML})
            self.assertGreater(len(modes), count)

        modes.clear()
        latex = output._repr_latex_()
        self.assertIn(r"\bar{J}", latex)
        self.assertEqual(set(modes), {PrintMode.Latex})
        count = len(modes)
        self.assertEqual(output._repr_latex_(), latex)
        self.assertEqual(len(modes), count)

        oversized = False
        modes.clear()
        html = output._repr_html_()
        self.assertIn("<mover><mi>J</mi><mo>¯</mo></mover>", html)
        self.assertEqual(set(modes), {PrintMode.MathML})
        count = len(modes)
        self.assertEqual(output._repr_html_(), html)
        self.assertEqual(output._repr_latex_(), latex)
        self.assertEqual(len(modes), count)
        self.assertEqual(value.to_expression(), original)

    def test_large_scalar_exports_print_every_term_and_notebooks_are_bounded(self):
        calls = []

        def printer(expression, *, mode, **options):
            calls.append(mode)

        head = sp.TensorName("full_display_tests::Term", print=printer)
        term_count = 448
        value = sp.TensorExpression(sum((head(i) for i in range(term_count)), start=0))
        original = value.to_expression()
        self.assertGreater(original.get_byte_size(), 4096)
        calls.clear()
        html = value.to_html()
        for i in range(term_count):
            self.assertIn(f"<mn>{i}</mn>", html)
        self.assertGreaterEqual(len(calls), term_count)
        self.assertEqual(set(calls), {PrintMode.MathML})
        calls.clear()
        latex = value.to_latex()
        self.assertIn(str(term_count - 1), latex)
        self.assertGreaterEqual(len(calls), term_count)
        self.assertEqual(set(calls), {PrintMode.Latex})
        output = sp.formatted(original)
        self.assertEqual(output._repr_html_(), html)
        self.assertEqual(output._repr_latex_(), latex)
        self.assertEqual(output.format_plain(), value.format_tensor())
        calls.clear()
        preview = self.assertBoundedNotebookPreview(value)
        self.assertGreater(len(calls), 0)
        self.assertLess(len(calls), term_count)
        self.assertEqual(set(calls), {PrintMode.MathML})
        self.assertLess(
            sum(f"<mn>{i}</mn>" in preview for i in range(term_count)), term_count
        )
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(value.rank, 0)

    def test_sum_header_boundary_keeps_both_complete_exported_terms(self):
        # Symbol ids and kernel encodings affect the exact byte boundary.
        # Vary an integer's encoded length one byte at a time so the sum header,
        # rather than either term, is what crosses the notebook byte limit.
        candidates = []
        found = False
        for argument_count in range(500, 800):
            arguments = ",".join(str(i) for i in range(argument_count))
            base = E(
                f"full_display_tests::f({arguments})+full_display_tests::g({arguments},0)"
            )
            if not 4064 <= base.get_byte_size() <= 4096:
                continue
            for padding in range(32):
                pad = 1 << (80 + 8 * padding)
                raw = E(
                    f"full_display_tests::f({arguments})+full_display_tests::g({arguments},0,{pad})"
                )
                sizes = (
                    raw.get_byte_size(),
                    sum(term.get_byte_size() for term in raw.terms()),
                )
                candidates.append(sizes)
                if sizes[0] > 4096 and sizes[1] <= 4096:
                    found = True
                    break
            if found:
                break
        self.assertTrue(
            found, f"no sum-header byte-boundary fixture found: {candidates}"
        )
        self.assertGreater(raw.get_byte_size(), 4096)
        self.assertLessEqual(sum(term.get_byte_size() for term in raw.terms()), 4096)
        value = sp.TensorExpression(raw)
        output = sp.formatted(raw)
        self.assertEqual(output.format_plain(), value.format_tensor())
        self.assertEqual(output._repr_latex_(), value.to_latex())
        html = value.to_html()
        self.assertEqual(output._repr_html_(), html)
        for i in range(argument_count):
            self.assertEqual(html.count(f"<mn>{i}</mn>"), 3 if i == 0 else 2)
        self.assertBoundedNotebookPreview(value, full_text=True)
        self.assertEqual(value.to_expression(), raw)

    def test_large_product_and_power_exports_render_every_factor(self):
        factor_count = 1024
        product = E("*".join(f"full_display_tests::x{i}" for i in range(factor_count)))
        self.assertGreater(product.get_byte_size(), 4096)
        for raw in (product, (product + 1) ** 7):
            with self.subTest(power=raw != product):
                value = sp.TensorExpression(raw)
                html = value.to_html()
                latex = value.to_latex()
                plain = value.format_tensor()
                output = sp.formatted(raw)
                self.assertEqual(html, output._repr_html_())
                self.assertEqual(plain, output.format_plain())
                self.assertEqual(latex, output._repr_latex_())
                for i in range(factor_count):
                    self.assertIn(f"x{i}", plain)
                    self.assertIn(f"x{i}", latex)
                    self.assertIn(f">x{i}<", html)
                self.assertBoundedNotebookPreview(value, full_text=True)
                if raw == product:
                    self.assertEqual(len(str(value)), 5033)
                self.assertEqual(value.to_expression(), raw)

    def test_large_display_keeps_logical_ports_metadata_and_status(self):
        m, e = sp.Representation.mink(4), sp.Representation.euc(3)
        tensor = sp.TensorName("full_display_tests::T")(m("a"), e("b"), m("c"))
        term_count = 2048
        coefficient = E(
            "+".join(f"full_display_tests::c{i}" for i in range(term_count))
        )
        value = (coefficient * tensor).with_name("full_display_tests::Named")
        original = value.to_expression()
        self.assertGreater(original.get_byte_size(), 4096)
        slots = [slot.to_expression() for slot in value.structure.axes]
        name, status = value.name.to_expression(), value.reduction_status
        html, latex, plain = value.to_html(), value.to_latex(), value.format_tensor()
        for i in range(term_count):
            self.assertIn(f">c{i}<", html)
            self.assertIn(f"c{i}", latex)
            self.assertIn(f"c{i}", plain)
        self.assertBoundedNotebookPreview(value)
        self.assertEqual(value.to_expression(), original)
        self.assertEqual([slot.to_expression() for slot in value.structure.axes], slots)
        self.assertEqual(value.name.to_expression(), name)
        self.assertEqual(value.reduction_status, status)
        self.assertEqual(value.rank, 3)

    def test_small_notation_and_typed_zero_keep_complete_display(self):
        r = sp.Representation.euc(3)
        a = sp.TensorName("full_display_tests::A")(r("i"), r("j"))
        zero = 0 * a
        for value in (a, zero):
            raw = value.to_expression()
            slots = [slot.to_expression() for slot in value.structure.axes]
            self.assertEqual(value._repr_html_(), value.to_html())
            self.assertEqual(value._repr_latex_(), value.to_latex())
            self.assertEqual(value.to_expression(), raw)
            self.assertEqual(
                [slot.to_expression() for slot in value.structure.axes], slots
            )

    def test_formatted_rejects_typst_notation_consistently(self):
        head = sp.TensorName("full_display_tests::FormattedNotation")(
            sp.Representation.euc(2)
        )
        tensor = sp.Tensor.dense(head, [1.0, 2.0])
        network = sp.TensorNetwork(tensor)
        for value in (head, tensor, network):
            with (
                self.subTest(kind=type(value).__name__),
                self.assertRaisesRegex(
                    ValueError, "does not accept Typst notation_source"
                ),
            ):
                value.formatted(notation_source="")

    def test_ready_formatted_output_remains_compatible(self):
        output = FormattedOutput("plain", "<em>rich</em>", "$$x$$")
        self.assertEqual(str(output), "plain")
        self.assertEqual(output._repr_html_(), "<em>rich</em>")
        self.assertEqual(output._repr_latex_(), "$$x$$")


if __name__ == "__main__":
    unittest.main()
