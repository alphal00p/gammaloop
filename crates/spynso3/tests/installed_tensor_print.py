"""Run with the installed Spenso extension and its optional Typst renderer."""

import re
import unittest

from symbolica import PrintMode, S
from symbolica.community.spenso import (
    DisplaySettings,
    Representation,
    TensorExpression,
    TensorLibrary,
    TensorName,
    dot,
    format_tensor,
)


def mathml(expression):
    return "".join(re.findall(r"<math\b.*?</math>", expression.to_html(), re.DOTALL))


class TensorPrintTests(unittest.TestCase):
    def test_head_mapping_keeps_arguments_coordinates_and_styles(self):
        name = TensorName(
            "tensor_print_tests::Jbar",
            print={"typst": "macron(J)", "latex": r"\bar{J}"},
        )
        x = S("tensor_print_tests::x")
        tensor = name(x, 7, Representation.euc(2), Representation.euc(2))
        component = TensorExpression(tensor.components()[1])
        original = component.to_expression()
        self.assertEqual(component.to_typst(), "attach(macron(J)(x,7),t:0 comma 1)")
        self.assertEqual(
            component._repr_latex_(), r"$$\bar{J}\!\left(x,7\right)^{0,1}$$"
        )
        self.assertEqual(component.format_tensor(), "Jbar(x,7)^(0,1)")
        array = component.to_typst(settings=DisplaySettings(component_style="array"))
        self.assertEqual(array, "macron(J)(x,7) lr([0 comma 1])")
        self.assertIn("<mover", mathml(component))
        self.assertIn("<msup", mathml(component))
        self.assertIn("<svg", component.to_svg())
        self.assertEqual(component.to_expression(), original)
        # Repeated notebook evaluation must accept the same mapping in any order.
        repeated = TensorName(
            "tensor_print_tests::Jbar",
            print={"latex": r"\bar{J}", "typst": "macron(J)"},
        )
        self.assertEqual(repeated.to_expression(), name.to_expression())

    def test_current_and_concrete_components_keep_the_custom_name(self):
        spinor = Representation.bis(4)
        jbar = TensorName(
            "tensor_print_tests::CurrentJbar", print={"typst": "macron(J)"}
        )(spinor)
        j = TensorName("tensor_print_tests::J")(spinor)
        current = jbar("a") * TensorExpression.gamma(4)("a", "b", "mu") * j("b")
        original = current.to_expression()
        library = TensorLibrary.hep_lib_atom()
        network = current.to_network(library=library)
        network.execute(library=library)
        kernel = network.result_tensor(library=library)
        self.assertEqual(len(kernel), 4)
        self.assertIn("macron(J)", current.to_typst())
        self.assertIn("macron(J)", kernel.to_typst())
        self.assertIn("<mover", mathml(current))
        self.assertIn("<mover", mathml(kernel))
        self.assertIn("<svg", kernel.to_svg())
        self.assertEqual(current.to_expression(), original)

    def test_callable_owns_the_display_and_none_uses_tensor_notation(self):
        def printer(expression, *, mode, **options):
            if mode == PrintMode.Typst:
                return "cal(J)"
            return None

        tensor = TensorName("tensor_print_tests::Callback", print=printer)(
            Representation.euc(2)
        )
        component = TensorExpression(tensor.components()[0])
        self.assertEqual(component.to_typst(), "cal(J)")
        self.assertNotIn("<msup", mathml(component))
        self.assertIn("𝒥", mathml(component))
        self.assertEqual(component._repr_latex_(), "$$Callback^{0}$$")
        power = TensorExpression(component**2)
        self.assertIn("cal(J)", power.to_typst())
        self.assertIn("𝒥", mathml(power))
        self.assertIn("<mn>2</mn>", mathml(power))

        def defer(expression, **options):
            return None

        fallback = TensorName("tensor_print_tests::Fallback", print=defer)(
            Representation.euc(2)
        )
        default = TensorExpression(fallback.components()[0])
        ordinary = TensorName("tensor_print_defaults::Fallback")(Representation.euc(2))
        expected = TensorExpression(ordinary.components()[0])
        self.assertEqual(default.to_typst(), expected.to_typst())
        self.assertIn("<msup", mathml(default))

    def test_vector_names_and_plain_mapping(self):
        vector = TensorName.vector(
            "tensor_print_tests::Vector", print={"typst": "macron(p)", "plain": "pbar"}
        )(Representation.mink(4))
        component = vector.components()[0]
        self.assertEqual(format_tensor(component), "pbar^(0)")
        self.assertIn("macron(p)", dot(vector, vector).to_typst())
        self.assertIn("<mover", mathml(dot(vector, vector)))

    def test_invalid_print_specifications(self):
        with self.assertRaisesRegex(ValueError, "mapping keys"):
            TensorName("tensor_print_tests::BadKey", print={"svg": "J"})
        with self.assertRaises(TypeError):
            TensorName("tensor_print_tests::BadValue", print={"typst": 1})
        with self.assertRaisesRegex(TypeError, "print must be"):
            TensorName("tensor_print_tests::BadType", print="J")


if __name__ == "__main__":
    unittest.main()
