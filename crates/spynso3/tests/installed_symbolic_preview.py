"""Bounded notebook display and lazy backends against the installed extension."""

import unittest
from unittest.mock import patch

import typst
from symbolica import E, FormattedOutput, PrintMode
from symbolica.community import tensor as sp


class SymbolicPreviewTests(unittest.TestCase):
    def test_formatted_only_computes_each_requested_backend_once(self):
        modes = []

        def printer(expression, *, mode, **options):
            modes.append(mode)
            return {
                PrintMode.Symbolica: "Jbar",
                PrintMode.Latex: r"\bar{J}",
                PrintMode.Typst: "macron(J)",
            }.get(mode)

        value = sp.TensorName("lazy_preview_tests::Current", print=printer)(
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
        expression = E("(lazy_preview_tests::x+1)^3")
        with patch.object(typst, "compile", wraps=typst.compile) as compile_:
            output = sp.formatted(expression)
            del expression
            compile_.assert_not_called()
            self.assertIn("3", output._repr_latex_())
            compile_.assert_not_called()
            self.assertIn("<msup", output._repr_html_())
            self.assertEqual(compile_.call_count, 1)

    def test_missing_renderer_is_cached_without_preventing_latex(self):
        output = sp.TensorExpression(E("lazy_preview_tests::x+1")).formatted()
        with patch.object(
            typst, "compile", side_effect=ImportError("no renderer")
        ) as compile_:
            self.assertIsNone(output._repr_html_())
            self.assertIsNone(output._repr_html_())
            self.assertEqual(compile_.call_count, 1)
        self.assertIn("x", output._repr_latex_())

    def test_large_scalar_preview_does_not_print_omitted_leaves(self):
        calls = []

        def printer(expression, *, mode, **options):
            calls.append(mode)

        head = sp.TensorName("lazy_preview_tests::Term", print=printer)
        value = sp.TensorExpression(sum((head(i) for i in range(1500)), start=0))
        original = value.to_expression()
        self.assertGreater(original.get_byte_size(), 4096)
        calls.clear()
        html = value._repr_html_()
        self.assertIn("<div data-spenso-preview>", html)
        self.assertLess(len(html), 30000)
        self.assertLess(len(calls), 50)
        calls.clear()
        latex = value._repr_latex_()
        self.assertLess(len(latex), 10000)
        self.assertLess(len(calls), 50)
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(value.rank, 0)

        # Explicit export must visit the complete expression, even after preview.
        calls.clear()
        full = value.to_latex()
        self.assertIn("1499", full)
        self.assertGreaterEqual(len(calls), 1500)
        self.assertEqual(value.to_expression(), original)

    def test_sum_header_boundary_still_omits_a_term_before_ellipsis(self):
        arguments = ",".join(str(i) for i in range(573))
        raw = E(
            f"lazy_preview_tests::f({arguments})+lazy_preview_tests::g({arguments},0)"
        )
        self.assertGreater(raw.get_byte_size(), 4096)
        self.assertLessEqual(sum(term.get_byte_size() for term in raw.terms()), 4096)
        value = sp.TensorExpression(raw)
        output = value.formatted()
        for preview in (
            output.format_plain(),
            output._repr_latex_(),
            output._repr_html_(),
        ):
            self.assertIn("1 of 2 terms shown", preview)
        self.assertEqual(value.to_expression(), raw)

    def test_large_single_product_and_power_do_not_bypass_preview(self):
        class Pretty:
            def __init__(self):
                self.parts = []

            def text(self, text):
                self.parts.append(text)

        product = E("*".join(f"lazy_preview_tests::x{i}" for i in range(1800)))
        for raw in (product, (product + 1) ** 7):
            with self.subTest(power=raw != product):
                value = sp.TensorExpression(raw)
                with patch.object(typst, "compile") as compile_:
                    html = value._repr_html_()
                    latex = value._repr_latex_()
                    plain = str(value.formatted())
                    free_html = sp.formatted(raw)._repr_html_()
                    pretty = Pretty()
                    value._repr_pretty_(pretty, False)
                    compile_.assert_not_called()
                self.assertIn("<div data-spenso-preview>", html)
                self.assertIn("<div data-spenso-preview>", free_html)
                self.assertLess(len(html), 10000)
                self.assertLess(len(latex), 1000)
                self.assertLess(len(plain), 1000)
                self.assertLess(len("".join(pretty.parts)), 1000)
                self.assertEqual(value.to_expression(), raw)

    def test_preview_keeps_logical_ports_metadata_and_status(self):
        m, e = sp.Representation.mink(4), sp.Representation.euc(3)
        tensor = sp.TensorName("lazy_preview_tests::T")(m("a"), e("b"), m("c"))
        coefficient = E("+".join(f"lazy_preview_tests::c{i}" for i in range(1500)))
        value = (coefficient * tensor).with_name("lazy_preview_tests::Named")
        original = value.to_expression()
        slots = [slot.to_expression() for slot in value.structure.slots]
        name, status = value.name.to_expression(), value.reduction_status
        with patch.object(typst, "compile") as compile_:
            self.assertIn("<div data-spenso-preview>", value._repr_html_())
            value._repr_latex_()
            value.formatted().format_plain()
            compile_.assert_not_called()
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(
            [slot.to_expression() for slot in value.structure.slots], slots
        )
        self.assertEqual(value.name.to_expression(), name)
        self.assertEqual(value.reduction_status, status)
        self.assertEqual(value.rank, 3)

    def test_small_notation_and_typed_zero_keep_complete_display(self):
        r = sp.Representation.euc(3)
        a = sp.TensorName("lazy_preview_tests::A")(r("i"), r("j"))
        zero = 0 * a
        for value in (a, zero):
            raw = value.to_expression()
            slots = [slot.to_expression() for slot in value.structure.slots]
            self.assertEqual(value._repr_html_(), value.to_html())
            self.assertEqual(value._repr_latex_(), value.to_latex())
            self.assertNotIn("<div data-spenso-preview>", value._repr_html_())
            self.assertEqual(value.to_expression(), raw)
            self.assertEqual(
                [slot.to_expression() for slot in value.structure.slots], slots
            )

    def test_ready_formatted_output_remains_compatible(self):
        output = FormattedOutput("plain", "<em>rich</em>", "$$x$$")
        self.assertEqual(str(output), "plain")
        self.assertEqual(output._repr_html_(), "<em>rich</em>")
        self.assertEqual(output._repr_latex_(), "$$x$$")


if __name__ == "__main__":
    unittest.main()
