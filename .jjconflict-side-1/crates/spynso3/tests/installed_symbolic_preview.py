"""Complete notebook display and lazy backends against the installed extension."""

import unittest
from unittest.mock import patch

import typst
from symbolica import E, FormattedOutput, PrintMode
from symbolica.community import tensor as sp


class SymbolicDisplayTests(unittest.TestCase):
    def test_formatted_only_computes_each_requested_backend_once(self):
        modes = []

        def printer(expression, *, mode, **options):
            modes.append(mode)
            return {
                PrintMode.Symbolica: "Jbar",
                PrintMode.Latex: r"\bar{J}",
                PrintMode.Typst: "macron(J)",
            }.get(mode)

        value = sp.TensorName("full_display_tests::Current", print=printer)(
            sp.Representation.euc(2)
        )
        original = value.to_expression()
        modes.clear()
        with patch.object(typst, "compile", wraps=typst.compile) as compile_:
            formatted = value.formatted()
            self.assertIs(type(formatted), FormattedOutput)
            self.assertEqual(modes, [])
            compile_.assert_not_called()

            plain = formatted.format_plain()
            self.assertIn("Jbar", plain)
            self.assertEqual(set(modes), {PrintMode.Symbolica})
            count = len(modes)
            self.assertEqual(str(formatted), plain)
            self.assertEqual(repr(formatted), plain)
            self.assertEqual(len(modes), count)
            compile_.assert_not_called()

            modes.clear()
            latex = formatted._repr_latex_()
            self.assertIn(r"\bar{J}", latex)
            self.assertEqual(set(modes), {PrintMode.Latex})
            count = len(modes)
            self.assertEqual(formatted._repr_latex_(), latex)
            self.assertEqual(len(modes), count)
            compile_.assert_not_called()

            modes.clear()
            html = formatted._repr_html_()
            self.assertIn("<mover", html)
            self.assertEqual(set(modes), {PrintMode.Typst})
            count = len(modes)
            self.assertEqual(formatted._repr_html_(), html)
            self.assertEqual(len(modes), count)
            self.assertEqual(compile_.call_count, 1)
        self.assertEqual(value.to_expression(), original)

    def test_free_formatted_owns_its_input_until_rendered(self):
        expression = E("(full_display_tests::x+1)^3")
        with patch.object(typst, "compile", wraps=typst.compile) as compile_:
            output = sp.formatted(expression)
            del expression
            compile_.assert_not_called()
            self.assertIn("3", output._repr_latex_())
            compile_.assert_not_called()
            self.assertIn("<msup", output._repr_html_())
            self.assertEqual(compile_.call_count, 1)

    def test_missing_renderer_is_cached_without_preventing_latex(self):
        output = sp.TensorExpression(E("full_display_tests::x+1")).formatted()
        with patch.object(
            typst, "compile", side_effect=ImportError("no renderer")
        ) as compile_:
            self.assertIsNone(output._repr_html_())
            self.assertIsNone(output._repr_html_())
            self.assertEqual(compile_.call_count, 1)
        self.assertIn("x", output._repr_latex_())

    def test_large_scalar_display_prints_every_term(self):
        calls = []

        def printer(expression, *, mode, **options):
            calls.append(mode)

        head = sp.TensorName("full_display_tests::Term", print=printer)
        term_count = 448
        value = sp.TensorExpression(sum((head(i) for i in range(term_count)), start=0))
        original = value.to_expression()
        self.assertGreater(original.get_byte_size(), 4096)
        calls.clear()
        html = value._repr_html_()
        for i in range(term_count):
            self.assertIn(f"<mn>{i}</mn>", html)
        self.assertGreaterEqual(len(calls), term_count)
        calls.clear()
        latex = value._repr_latex_()
        self.assertIn(str(term_count - 1), latex)
        self.assertGreaterEqual(len(calls), term_count)
        self.assertEqual(latex, value.to_latex())
        output = value.formatted()
        self.assertEqual(output._repr_html_(), html)
        self.assertEqual(output._repr_latex_(), latex)
        self.assertEqual(output.format_plain(), value.format_tensor())
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(value.rank, 0)

    def test_sum_header_boundary_keeps_both_complete_terms(self):
        arguments = ",".join(str(i) for i in range(573))
        raw = E(
            f"full_display_tests::f({arguments})+full_display_tests::g({arguments},0)"
        )
        self.assertGreater(raw.get_byte_size(), 4096)
        self.assertLessEqual(sum(term.get_byte_size() for term in raw.terms()), 4096)
        value = sp.TensorExpression(raw)
        output = value.formatted()
        self.assertEqual(output.format_plain(), value.format_tensor())
        self.assertEqual(output._repr_latex_(), value.to_latex())
        self.assertEqual(output._repr_html_(), value.to_html())
        self.assertEqual(value.to_expression(), raw)

    def test_large_product_and_power_render_every_factor(self):
        class Pretty:
            def __init__(self):
                self.parts = []

            def text(self, text):
                self.parts.append(text)

        factor_count = 1024
        product = E("*".join(f"full_display_tests::x{i}" for i in range(factor_count)))
        self.assertGreater(product.get_byte_size(), 4096)
        for raw in (product, (product + 1) ** 7):
            with self.subTest(power=raw != product):
                value = sp.TensorExpression(raw)
                with patch.object(typst, "compile", wraps=typst.compile) as compile_:
                    html = value._repr_html_()
                    latex = value._repr_latex_()
                    plain = str(value.formatted())
                    free_html = sp.formatted(raw)._repr_html_()
                    pretty = Pretty()
                    value._repr_pretty_(pretty, False)
                    self.assertEqual(compile_.call_count, 2)
                self.assertEqual(html, free_html)
                self.assertEqual(plain, value.format_tensor())
                self.assertEqual(latex, value.to_latex())
                self.assertEqual("".join(pretty.parts), plain)
                for i in range(factor_count):
                    self.assertIn(f"x{i}", plain)
                    self.assertIn(f"x{i}", latex)
                    self.assertIn(f">x{i}<", html)
                self.assertEqual(value.to_expression(), raw)

    def test_large_display_keeps_logical_ports_metadata_and_status(self):
        m, e = sp.Representation.mink(4), sp.Representation.euc(3)
        tensor = sp.TensorName("full_display_tests::T")(m("a"), e("b"), m("c"))
        coefficient = E("+".join(f"full_display_tests::c{i}" for i in range(1024)))
        value = (coefficient * tensor).with_name("full_display_tests::Named")
        original = value.to_expression()
        self.assertGreater(original.get_byte_size(), 4096)
        slots = [slot.to_expression() for slot in value.structure.axes]
        name, status = value.name.to_expression(), value.reduction_status
        self.assertEqual(value._repr_html_(), value.to_html())
        self.assertEqual(value._repr_latex_(), value.to_latex())
        self.assertEqual(value.formatted().format_plain(), value.format_tensor())
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

    def test_ready_formatted_output_remains_compatible(self):
        output = FormattedOutput("plain", "<em>rich</em>", "$$x$$")
        self.assertEqual(str(output), "plain")
        self.assertEqual(output._repr_html_(), "<em>rich</em>")
        self.assertEqual(output._repr_latex_(), "$$x$$")


if __name__ == "__main__":
    unittest.main()
