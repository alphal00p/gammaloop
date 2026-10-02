"""Run with the installed Spenso extension and its optional Typst renderer."""

import re
import tempfile
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import typst
from symbolica import E, PrintMode, S
from symbolica.community.tensor import (
    DisplaySettings,
    Representation,
    TensorExpression,
    TensorLibrary,
    TensorName,
    dot,
    format_tensor,
)


def mathml(expression):
    return "".join(
        re.findall(
            r"<math\b.*?</math>",
            expression.to_html(settings=DisplaySettings(tensor_view="matrix")),
            re.DOTALL,
        )
    )


class TensorPrintTests(unittest.TestCase):
    def test_grouped_factors_keep_their_mathml_parentheses(self):
        for source, expected_groups in (
            ("(1+1i)*x", 1),
            ("(-1+2i)*x", 1),
            ("(1/2+1i/3)*x", 1),
            ("x+(-1+2i)*y", 1),
            ("(1+1i)^x", 1),
            ("(2i)^x", 1),
            ("x-(a+b)", 1),
            ("-(-a+b)-(a-b+c)", 2),
            ("-(a+b)", 1),
            ("x-2*(a+b)", 1),
            ("x-a", 0),
        ):
            with self.subTest(source=source):
                expression = TensorExpression(E(source))
                original = expression.to_expression()
                root = ET.fromstring(mathml(expression))
                groups = [
                    node
                    for node in root.iter("mrow")
                    if len(node) >= 2
                    and node[0].tag == "mo"
                    and node[0].text == "("
                    and node[-1].tag == "mo"
                    and node[-1].text == ")"
                ]
                self.assertEqual(len(groups), expected_groups)
                self.assertEqual(expression.to_expression(), original)
                if "^" in source:
                    self.assertIs(root.find(".//msup")[0], groups[0])

    def test_pure_imaginary_factors_stay_compact(self):
        for source in ("1i*x", "-2i*x", "1i*x/2", "x-2i*y"):
            with self.subTest(source=source):
                root = ET.fromstring(mathml(TensorExpression(E(source))))
                operators = [node.text for node in root.iter("mo")]
                self.assertNotIn("(", operators)
                self.assertNotIn(")", operators)
                if source == "x-2i*y":
                    # Symbolica orders terms by symbol registration. Both exact
                    # orders retain the negative imaginary coefficient compactly.
                    visible = re.sub(
                        r"\s+",
                        "",
                        unicodedata.normalize("NFKC", "".join(root.itertext())),
                    )
                    self.assertIn(visible, ("x−2iy", "−2iy+x"))
                else:
                    self.assertNotIn("+", operators)

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
        current = jbar("a") * TensorExpression.dirac_gamma(4)("a", "b", "mu") * j("b")
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

    def test_portable_custom_cancellation_retains_its_html_stroke(self):
        assets = Path(__file__).resolve().parents[1] / "typst"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "render.typ").write_bytes((assets / "render.typ").read_bytes())
            # Exercise the portable custom-display AST independently of gamma's
            # specialized renderer, including Typst's lexical helper scope.
            (root / "notation.typ").write_text(
                (assets / "notation.typ").read_text()
                + '\n#let custom-cancel() = _display-node(("math", "cancel", (("symbol", "J"),)))\n'
            )
            source = root / "main.typ"
            source.write_text(
                '#import "notation.typ": custom-cancel\n$ #custom-cancel() $'
            )
            document = typst.compile(
                str(source), format="html", root=directory
            ).decode()
            svg = typst.compile(str(source), format="svg", root=directory).decode()
        root = ET.fromstring(re.search(r"<math\b.*?</math>", document, re.DOTALL)[0])
        strokes = [node for node in root.iter() if "data-spenso-cancel" in node.attrib]
        self.assertEqual(len(strokes), 1)
        self.assertEqual(strokes[0].attrib["data-spenso-cancel"], "updiagonal")
        self.assertEqual(
            unicodedata.normalize("NFKC", "".join(strokes[0].itertext())).strip(), "J"
        )
        self.assertIn("<svg", svg)
        self.assertNotIn("data-spenso-cancel", svg)

    def test_invalid_print_specifications(self):
        with self.assertRaisesRegex(ValueError, "mapping keys"):
            TensorName("tensor_print_tests::BadKey", print={"svg": "J"})
        with self.assertRaises(TypeError):
            TensorName("tensor_print_tests::BadValue", print={"typst": 1})
        with self.assertRaisesRegex(TypeError, "print must be"):
            TensorName("tensor_print_tests::BadType", print="J")


if __name__ == "__main__":
    unittest.main()
