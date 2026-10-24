"""Exercise direct MathML print mappings and optional Typst SVG exports."""

import importlib.util
import re
import subprocess
import sys
import tempfile
import textwrap
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

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
from symbolica.community.tensor_native import _compile_typst


def mathml(expression):
    sources = re.findall(
        r"<math\b.*?</math>",
        expression.to_html(settings=DisplaySettings(tensor_view="matrix")),
        re.DOTALL,
    )
    roots = [ET.fromstring(source) for source in sources]
    for root in roots:
        for node in root.iter():
            node.tag = node.tag.rsplit("}", 1)[-1]
    return "".join(ET.tostring(root, encoding="unicode") for root in roots)


class TensorPrintTests(unittest.TestCase):
    def test_portable_export_is_explicit_in_fresh_processes(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = str(Path(directory) / "plain.bin")
            producer = textwrap.dedent("""
                import sys
                from symbolica.community.tensor import TensorName, to_typst_project
                head = TensorName("portable_print::Head", print={"plain": "ProcessPortableLabel"}).to_expression()
                head.save(sys.argv[1])
                project = to_typst_project(head)
                assert b"spenso.notation" in project["tree.cbor"]
                assert b"ProcessPortableLabel" in project["tree.cbor"]
                assert b"ProcessPortableLabel" not in open(sys.argv[1], "rb").read()
            """)
            consumer = textwrap.dedent("""
                import sys
                from symbolica import Expression
                from symbolica.community.tensor import to_typst_project
                project = to_typst_project(Expression.load(sys.argv[1]))
                assert b"ProcessPortableLabel" not in project["tree.cbor"]
            """)
            subprocess.run([sys.executable, "-c", producer, filename], check=True)
            subprocess.run([sys.executable, "-c", consumer, filename], check=True)

    def test_print_mapping_validation_and_registration(self):
        mapping = {
            "plain": "J",
            "latex": r"\bar{J}",
            "typst": "overline(J)",
            "mathml": "<mover accent='true'><mi>J</mi><mo>¯</mo></mover>",
        }
        TensorName("validated_print::J", print=mapping)
        TensorName("validated_print::J", print=mapping)
        with self.assertRaisesRegex(ValueError, "different process-local notation"):
            TensorName("validated_print::J", print={"plain": "K"})
        for source, offending in [
            ("<img src=x onerror=alert(1)>", "MathML"),
            ("<img src='x' onerror='alert(1)'/>", "img"),
            ("<mi onclick='x'>J</mi>", "onclick"),
            ("<mi href='x'>J</mi>", "href"),
            ("<mi style='color:red'>J</mi>", "style"),
            ("<mi>&unknown;</mi>", "unknown"),
        ]:
            with (
                self.subTest(source=source),
                self.assertRaisesRegex(ValueError, offending),
            ):
                TensorName("validated_print::Rejected", print={"mathml": source})

    def test_saved_notation_is_process_local_and_old_tags_are_inert(self):
        with tempfile.TemporaryDirectory() as directory:
            filename = str(Path(directory) / "expression.bin")
            producer = textwrap.dedent("""
                import sys
                from symbolica import S
                from symbolica.community.tensor import TensorName
                head = TensorName("saved_print::SavedHead", print={"mathml": "<mi>LocalOnly</mi>"}).to_expression()
                legacy = S("saved_print::LegacyHead", tags=["spenso::tensor", "spenso::print-head:mathml:<img src=x onerror=boom>"])
                (head(1) + legacy(2)).save(sys.argv[1])
            """)
            consumer = textwrap.dedent("""
                import sys
                from symbolica import Expression
                from symbolica.community.tensor import TensorExpression
                html = TensorExpression(Expression.load(sys.argv[1])).to_html()
                assert "SavedHead" in html and "LegacyHead" in html, html
                assert "LocalOnly" not in html and "<img" not in html, html
            """)
            subprocess.run([sys.executable, "-c", producer, filename], check=True)
            subprocess.run([sys.executable, "-c", consumer, filename], check=True)

    def test_scalar_print_callbacks_do_not_require_tensor_tags(self):
        head = S(
            "tensor_print_tests::ScalarEnergy",
            is_scalar=True,
            print={"typst": "cal(E)", "mathml": "<mi>ℰ</mi>"},
        )
        expression = TensorExpression(head(27) ** 2)
        before = expression.to_expression()
        output = mathml(expression)
        self.assertIn("ℰ", output)
        self.assertNotIn("ScalarEnergy", output)
        self.assertIn("<mn>2</mn>", output)
        self.assertEqual(expression.to_expression(), before)

        deferred = S(
            "tensor_print_tests::DeferredScalar",
            print=lambda *args, **kwargs: None,
        )
        output = mathml(TensorExpression(deferred(27)))
        self.assertIn("DeferredScalar", output)
        self.assertIn("<mn>27</mn>", output)
        self.assertNotIn("<msub", output)

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
                    self.assertIn(groups[0], root.find(".//msup")[0].iter())

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
            print={
                "typst": "macron(J)",
                "latex": r"\bar{J}",
                "mathml": "<mover><mi>J</mi><mo>¯</mo></mover>",
            },
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
        self.assertEqual(component.to_expression(), original)
        # Repeated notebook evaluation must accept the same mapping in any order.
        repeated = TensorName(
            "tensor_print_tests::Jbar",
            print={
                "latex": r"\bar{J}",
                "typst": "macron(J)",
                "mathml": "<mover><mi>J</mi><mo>¯</mo></mover>",
            },
        )
        self.assertEqual(repeated.to_expression(), name.to_expression())

    def test_current_and_concrete_components_keep_the_custom_name(self):
        spinor = Representation.bis(4)
        jbar = TensorName(
            "tensor_print_tests::CurrentJbar",
            print={
                "typst": "macron(J)",
                "mathml": "<mover><mi>J</mi><mo>¯</mo></mover>",
            },
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
        self.assertEqual(current.to_expression(), original)

    def test_callable_owns_the_display_and_none_uses_tensor_notation(self):
        def printer(expression, *, mode, **options):
            if mode == PrintMode.Typst:
                return "cal(J)"
            if mode == PrintMode.MathML:
                return "<mi>𝒥</mi>"
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
            "tensor_print_tests::Vector",
            print={
                "typst": "macron(p)",
                "plain": "pbar",
                "mathml": "<mover><mi>p</mi><mo>¯</mo></mover>",
            },
        )(Representation.mink(4))
        component = vector.components()[0]
        self.assertEqual(format_tensor(component), "pbar^(0)")
        self.assertIn("macron(p)", dot(vector, vector).to_typst())
        self.assertIn("<mover", mathml(dot(vector, vector)))

    def test_custom_mathml_cancellation_retains_its_html_stroke(self):
        def printer(expression, *, mode, **options):
            if mode == PrintMode.MathML:
                return (
                    '<mrow data-spenso-cancel="updiagonal" '
                    'style="background-image:linear-gradient(to bottom right,transparent 49%,currentColor 50%,transparent 51%)">'
                    "<mi>J</mi></mrow>"
                )
            return None

        value = TensorName("tensor_print_tests::Cancelled", print=printer)(7)
        original = bytes(value.__getstate__())
        root = ET.fromstring(mathml(value))
        strokes = [node for node in root.iter() if "data-spenso-cancel" in node.attrib]
        self.assertEqual(len(strokes), 1)
        self.assertEqual(strokes[0].attrib["data-spenso-cancel"], "updiagonal")
        self.assertEqual(
            unicodedata.normalize("NFKC", "".join(strokes[0].itertext())).strip(), "J"
        )
        self.assertIn("linear-gradient", strokes[0].attrib["style"])
        self.assertIn("currentColor", strokes[0].attrib["style"])
        self.assertEqual(bytes(value.__getstate__()), original)

    @unittest.skipUnless(
        importlib.util.find_spec("typst") is not None,
        "optional real SVG export requires the typst package",
    )
    def test_portable_custom_cancellation_compiles_as_explicit_svg_source(self):
        assets = Path(__file__).resolve().parents[1] / "typst"
        # Exercise the portable custom-display AST independently of gamma's
        # specialized renderer, including Typst's lexical helper scope.
        notation = (
            (assets / "notation.typ").read_text()
            + '\n#let custom-cancel() = _display-node(("math", "cancel", (("symbol", "J"),)))\n'
        )
        source = '#import "notation.typ": custom-cancel\n$ #custom-cancel() $'
        svg = bytes(_compile_typst(source, "svg", notation)).decode()
        self.assertIn("<svg", svg)
        self.assertNotIn("data-spenso-cancel", svg)

    @unittest.skipUnless(
        importlib.util.find_spec("typst") is not None,
        "optional real SVG export requires the typst package",
    )
    def test_custom_names_compile_in_symbolic_and_concrete_svg_exports(self):
        spinor = Representation.bis(4)
        jbar = TensorName("tensor_print_svg::Jbar", print={"typst": "macron(J)"})(
            spinor
        )
        j = TensorName("tensor_print_svg::J")(spinor)
        current = jbar("a") * TensorExpression.dirac_gamma(4)("a", "b", "mu") * j("b")
        library = TensorLibrary.hep_lib_atom()
        kernel = current.to_tensor(library)
        component = TensorExpression(jbar.components()[0])
        for value in (component, current, kernel):
            with self.subTest(value=value):
                self.assertIn("macron(J)", value.to_typst())
                root = ET.fromstring(value.to_svg())
                self.assertEqual(root.tag.rsplit("}", 1)[-1], "svg")

    def test_invalid_print_specifications(self):
        with self.assertRaisesRegex(ValueError, "mapping keys"):
            TensorName("tensor_print_tests::BadKey", print={"svg": "J"})
        with self.assertRaises(TypeError):
            TensorName("tensor_print_tests::BadValue", print={"typst": 1})
        with self.assertRaisesRegex(TypeError, "print must be"):
            TensorName("tensor_print_tests::BadType", print="J")


if __name__ == "__main__":
    unittest.main()
