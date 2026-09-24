"""Exercise notebook payloads against the installed native Spenso extension."""

import json
import re
import unittest
from html.parser import HTMLParser

from symbolica import E, Expression, S
from symbolica.community.spenso import (
    DisplaySettings,
    Representation,
    Tensor,
    TensorExpression,
    TensorName,
)


class Frames(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.frames = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        if tag == "iframe":
            self.frames.append(dict(attrs))


def explorer(tensor, **kwargs):
    frames = Frames(tensor.to_html(**kwargs)).frames
    assert len(frames) == 1
    assert frames[0]["sandbox"] == "allow-scripts"
    document = frames[0]["srcdoc"]
    match = re.search(
        r'<script id="tensor-data" type="application/json">(.*?)</script>',
        document,
        re.DOTALL,
    )
    assert match is not None
    return json.loads(match[1]), document


class TensorExplorerTests(unittest.TestCase):
    def test_rank_three_payload_agrees_with_coordinate_indexing(self):
        # Mixed representations deliberately exercise canonical storage permutations.
        m, e = Representation.mink(2), Representation.euc(3)
        structure = TensorName("explorer_tests::T")(m(1), e(2), m(3))
        values = [(S("explorer_tests::x") + 1) ** i for i in range(12)]
        tensor = Tensor.dense(structure, values)
        payload, document = explorer(tensor)
        self.assertEqual(payload["shape"], [2, 3, 2])
        self.assertEqual(payload["axes"], ["mu", "j", "rho"])
        self.assertTrue(payload["complete"])
        self.assertEqual(len(payload["entries"]), 12)
        for coordinate, (bytes_, plain, html) in payload["entries"]:
            self.assertEqual(bytes_, tensor[coordinate].get_byte_size())
            self.assertTrue(plain)
            self.assertIn("<svg", html)
        self.assertIn("Matrix", document)
        self.assertIn("data-spenso-explorer", tensor._repr_html_())
        self.assertEqual(list(tensor), values)

    def test_numeric_size_is_independent_of_magnitude(self):
        structure = TensorName("explorer_tests::N")(Representation.euc(2))
        for values, expected in [([0.0, 1e300], 8), ([0j, 1e200 + 1e100j], 16)]:
            payload, _ = explorer(Tensor.dense(structure, values))
            self.assertEqual(
                [value[0] for _, value in payload["entries"]], [expected] * 2
            )

    def test_sparse_default_is_shared_and_display_does_not_densify(self):
        structure = TensorName("explorer_tests::S")(Representation.euc(1_000_000))
        tensor = Tensor.sparse(structure, Expression)
        tensor[999_999] = E("x+y")
        payload, document = explorer(tensor)
        self.assertTrue(payload["sparse"])
        self.assertTrue(payload["complete"])
        self.assertEqual(payload["stored"], 1)
        self.assertEqual(payload["entries"][0][0], [999_999])
        self.assertEqual(payload["default"][0], E("0").get_byte_size())
        self.assertLess(len(document), 100_000)
        again, _ = explorer(tensor)
        self.assertEqual(again["stored"], 1)

    def test_dense_and_sparse_previews_are_explicitly_incomplete(self):
        structure = TensorName("explorer_tests::Large")(Representation.euc(2000))
        dense = Tensor.dense(structure, [float(i) for i in range(2000)])
        payload, _ = explorer(dense)
        self.assertFalse(payload["complete"])
        self.assertEqual(len(payload["entries"]), 512)
        self.assertEqual(payload["entries"][0][0], [0])
        self.assertEqual(payload["entries"][-1][0], [1999])
        dense.to_sparse()
        payload, _ = explorer(dense)
        self.assertFalse(payload["complete"])
        self.assertLessEqual(len(payload["entries"]), 512)
        self.assertEqual(payload["stored"], 1999)

    def test_matrix_override_keeps_mathematical_output(self):
        structure = TensorName("explorer_tests::M")(
            Representation.euc(2), Representation.euc(2)
        )
        tensor = Tensor.dense(structure, [1.0, 2.0, 3.0, 4.0])
        before = tensor.to_typst()
        html = tensor.to_html(settings=DisplaySettings(tensor_view="matrix"))
        self.assertIn("<mtable", html)
        self.assertNotIn("<iframe", html)
        self.assertEqual(tensor.to_typst(), before)
        self.assertEqual(DisplaySettings().tensor_view, "interactive")
        with self.assertRaisesRegex(ValueError, "tensor_view"):
            DisplaySettings(tensor_view="invalid")

    def test_zero_rank_stays_scalar_and_zero_length_stays_empty(self):
        scalar = Tensor.dense(
            TensorExpression(E("1")).with_name("explorer_tests::Scalar"), [2.0]
        )
        self.assertNotIn("<iframe", scalar.to_html())
        empty = Tensor.dense(
            TensorName("explorer_tests::Empty")(Representation.euc(0)), []
        )
        payload, _ = explorer(empty)
        self.assertEqual(payload["shape"], [0])
        self.assertEqual(payload["entries"], [])
        self.assertTrue(payload["complete"])

    def test_component_text_cannot_close_the_json_script(self):
        name = TensorName(
            "explorer_tests::Escaping",
            print={"plain": "</script><script>bad()</script>"},
        )
        structure = name(Representation.euc(1))
        tensor = Tensor.dense(structure, structure.components())
        payload, document = explorer(tensor)
        self.assertIn("</script>", payload["entries"][0][1][1])
        self.assertNotIn("<script>bad()", document)


if __name__ == "__main__":
    unittest.main()
