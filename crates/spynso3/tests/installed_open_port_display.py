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
    trace,
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
        for document in (native, rich):
            root = mathml(document)
            if isinstance(expected, tuple):
                # Indexed components commute; each factor must retain its ports.
                self.assertCountEqual([visible(factor) for factor in root], expected)
            else:
                self.assertEqual(visible(root), expected)
        self.assertNotIn("open_index", source)
        self.assertNotIn("open_index", expression.to_latex(settings=settings))
        self.assertNotIn("open_index", expression.format_tensor())
        self.assertEqual(expression.to_expression(), original)
        self.assertEqual(expression.structure, structure)
        return list(visible_nodes(mathml(rich)))

    def test_dot_products_keep_tight_contractions_and_separate_factors(self):
        space = Representation.mink(4)
        p, q, r, s = (
            TensorName.vector(f"dot_factor_spacing::{name}")(space)
            for name in ("p", "q", "r", "s")
        )
        single, product = p * q, (p * q) * (r * s)
        for value in (single, product):
            source = value.to_typst()
            native = typst.compile(f"$ {source} $".encode(), format="html").decode()
            for document in (native, value.to_html()):
                with self.subTest(product=value == product, native=document == native):
                    nodes = list(mathml(document))
                    dots = [i for i, node in enumerate(nodes) if node.text == "⋅"]
                    self.assertEqual(len(dots), 2 if value == product else 1)
                    for i in dots:
                        # Keep ordinary glyph spacing inside each contraction.
                        self.assertEqual(nodes[i].tag, "mi")
                        self.assertNotEqual(nodes[i - 1].tag, "mspace")
                        self.assertNotEqual(nodes[i + 1].tag, "mspace")
                    if value == product:
                        gaps = [
                            n for n in nodes[dots[0] + 1 : dots[1]] if n.tag == "mspace"
                        ]
                        self.assertEqual(len(gaps), 1)
                        self.assertEqual(gaps[0].attrib["width"], "0.12em")
                    else:
                        self.assertFalse(any(n.tag == "mspace" for n in nodes))

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
        scalar_product = dot(p, p)
        self.assert_renderers_agree(scalar_product, "p⋅p")
        self.assertNotIn("□", scalar_product.format_tensor())
        self.assertNotIn(r"\square", scalar_product.to_latex())
        self.assert_renderers_agree(gamma * gamma, "[γ□γ□]□0□1")
        chain = (gamma * gamma).index("a", "c", "μ", "ν")
        indexed = chain.undo_chain()
        self.assert_renderers_agree(indexed, ("γabμ", "γbcν"))
        self.assertEqual(indexed.contract().canonize(), chain.canonize())
        self.assertEqual(indexed.structure, chain.structure)
        self.assertNotIn("□", indexed.format_tensor())
        self.assertEqual(indexed.rank, 4)

    def test_slashes_keep_only_their_actual_spinor_ports(self):
        gamma = TensorExpression.dirac_gamma(4)
        p = TensorName.vector("open_slash_display::p")(Representation.mink(4))
        for left, right, expected, axes in (
            ("a", "b", "pab", []),
            (AUTO, AUTO, "p□0□1", ["0", "1"]),
            ("a", AUTO, "pa□", ["1"]),
        ):
            with self.subTest(left=left, right=right):
                slash = (gamma(left, right, "mu") * p("mu")).contract()
                self.assertEqual(slash.rank, 2)
                nodes = self.assert_renderers_agree(slash, expected)
                self.assertEqual(
                    [
                        node.attrib["data-spenso-open-axis"]
                        for node in nodes
                        if "data-spenso-open-axis" in node.attrib
                    ],
                    axes,
                )
                self.assertEqual(
                    sum("data-spenso-cancel" in node.attrib for node in nodes), 1
                )
                self.assertIn("cancel(p)", slash.to_typst())
                self.assertIn("slash(p)", slash.format_tensor())
                self.assertIn(r"\not{p}", slash.to_latex())
                self.assertEqual(slash.format_tensor().count("□"), len(axes))
                self.assertEqual(slash.to_latex().count(r"\square"), len(axes))
                self.assertEqual(
                    visible(
                        mathml(slash.to_html(settings=DisplaySettings.schoonschip()))
                    ),
                    expected,
                )

    def test_slash_chains_and_traces_retain_matrix_notation(self):
        gamma = TensorExpression.dirac_gamma(4)
        mink = Representation.mink(4)
        p, q = [
            TensorName.vector(f"open_slash_display::{name}")(mink)
            for name in ("p", "q")
        ]
        left = (gamma(AUTO, AUTO, "mu") * p("mu")).contract()
        right = (gamma(AUTO, AUTO, "nu") * q("nu")).contract()
        for expression, expected, rank in (
            (left * right, "[pq]□0□1", 2),
            ((left * right)("a", "b"), "[pq]ab", 2),
            (trace(Representation.bis(4), left, right), "Tr(pq)", 0),
        ):
            with self.subTest(expected=expected):
                nodes = self.assert_renderers_agree(expression, expected)
                self.assertEqual(expression.rank, rank)
                self.assertEqual(expression.to_typst().count("cancel("), 2)
                self.assertEqual(
                    sum("data-spenso-cancel" in node.attrib for node in nodes), 2
                )
                self.assertNotIn("⋅", visible(mathml(expression.to_html())))

    def test_contracted_momentum_difference_retains_each_slash_and_spinor_port(self):
        gamma = TensorExpression.dirac_gamma(4)
        mink = Representation.mink(4)
        p, q = [
            TensorName.vector(f"open_slash_display::{name}")(mink)
            for name in ("p", "q")
        ]
        difference = (gamma(AUTO, AUTO, "mu") * (p("mu") - q("mu"))).contract()
        nodes = self.assert_renderers_agree(difference, "p□0□1−q□0□1")
        self.assertEqual(difference.rank, 2)
        self.assertEqual(sum("data-spenso-cancel" in node.attrib for node in nodes), 2)
        self.assertEqual(difference.to_typst().count("cancel("), 2)
        self.assertEqual(difference.format_tensor().count("□"), 4)
        self.assertEqual(difference.to_latex().count(r"\square"), 4)

    def test_schoonschip_vectors_do_not_hide_uncontracted_tensor_ports(self):
        mink = Representation.mink(4)
        tensor = TensorName("open_slash_display::T")(mink, mink)
        p = TensorName.vector("open_slash_display::p")(mink)
        settings = DisplaySettings.schoonschip()
        for remaining, expected in (("nu", "Tpnu"), (AUTO, "Tp□")):
            with self.subTest(remaining=remaining):
                value = (tensor("mu", remaining) * p("mu")).contract()
                original, structure = value.to_expression(), value.structure
                self.assertEqual(
                    visible(mathml(value.to_html(settings=settings))), expected
                )
                self.assertEqual(value.rank, 1)
                self.assertEqual(value.format_tensor().count("□"), remaining is AUTO)
                self.assertEqual(value.to_latex().count(r"\square"), remaining is AUTO)
                self.assertEqual(value.to_expression(), original)
                self.assertEqual(value.structure, structure)

    def test_default_ports_preserve_existing_bra_and_ket_notation(self):
        mink = Representation.mink(4)
        tensor = TensorName("braket_display::T")(mink, mink)
        p = TensorName.vector("braket_display::p")(mink)
        rep = Representation("braket_display::V", 3, is_self_dual=False)
        bra = TensorName.vector("braket_display::a")(rep.dual())
        ket = TensorName.vector("braket_display::b")(rep)
        matrix = TensorName("braket_display::M")
        for remaining in ("ν", AUTO):
            index = "□" if remaining is AUTO else remaining
            cases = (
                ((tensor("mu", remaining) * p("mu")).contract(), "⟨p|T", index),
                # Exercise the existing compact representation directly so this
                # display test does not depend on contraction order or polarity.
                (
                    matrix(
                        bra.to_expression(),
                        ket.to_expression(),
                        mink if remaining is AUTO else mink(remaining),
                    ),
                    "⟨a|M",
                    index + "|b⟩",
                ),
            )
            for value, prefix, suffix in cases:
                with self.subTest(prefix=prefix, remaining=remaining):
                    original, structure = value.to_expression(), value.structure
                    source = value.to_typst()
                    native = typst.compile(
                        f"$ {source} $".encode(), format="html"
                    ).decode()
                    rich = value.to_html()
                    for document in (native, rich):
                        text = visible(mathml(document))
                        self.assertTrue(text.startswith(prefix), text)
                        self.assertTrue(text.endswith(suffix), text)
                        self.assertEqual(text.count("⟨"), 1)
                        self.assertEqual(text.count("⟩"), prefix == "⟨a|M")
                    open_ports = [
                        node
                        for node in visible_nodes(mathml(rich))
                        if "data-spenso-open-axis" in node.attrib
                    ]
                    self.assertEqual(len(open_ports), remaining is AUTO)
                    self.assertEqual(value.rank, 1)
                    self.assertEqual(value.to_expression(), original)
                    self.assertEqual(value.structure, structure)


if __name__ == "__main__":
    unittest.main()
