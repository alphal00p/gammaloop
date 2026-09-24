"""Executable catalogue recipes and component snapshots in the installed host."""

# ruff: noqa: S102 -- execute the renderer's own examples to check they work.

import contextlib
import io
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
    TensorLibrary,
    TensorName,
)


class Catalogue(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.entries = []
        self.entry = None
        self.recipe = None
        self.code = None
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "section" and attrs.get("class") == "sl-panel":
            self.entry = {"recipes": {}, "frames": []}
            self.entries.append(self.entry)
        if self.entry is None:
            return
        if tag == "div" and "data-recipe" in attrs:
            self.recipe = attrs["data-recipe"]
        if tag == "pre":
            self.code = ""
        if tag == "iframe":
            self.entry["frames"].append(attrs)

    def handle_endtag(self, tag):
        if tag == "pre":
            self.entry["recipes"][self.recipe] = self.code
            self.code = None
        if tag == "section":
            self.entry = None

    def handle_data(self, data):
        if self.code is not None:
            self.code += data


class LibraryDisplayTests(unittest.TestCase):
    def test_name_only_lookup_keeps_printers_and_explicit_options_still_validate(self):
        original = TensorName(
            "library_display::reused",
            is_symmetric=True,
            print={"typst": "macron(R)"},
        )
        reused = TensorName("library_display::reused")
        self.assertEqual(reused.to_expression(), original.to_expression())
        self.assertEqual(
            reused.to_expression().get_attributes(),
            original.to_expression().get_attributes(),
        )
        with self.assertRaises(TypeError):
            TensorName("library_display::reused", is_symmetric=False)

    def test_empty_catalogue_explains_factories(self):
        library = TensorLibrary()
        html = library.to_html()
        self.assertIn("data-spenso-library", html)
        self.assertIn("No stored tensors", html)
        self.assertIn("9 dimension-dependent factories", html)
        self.assertEqual(Catalogue(html).entries, [])
        self.assertEqual(len(library), 0)

    def test_hep_recipes_retrieve_real_data_and_keep_symbolic_indices(self):
        library = TensorLibrary.hep_lib_atom()
        before = [key.structure for key in library]
        html = library._repr_html_()
        self.assertIsInstance(html, str)
        catalogue = Catalogue(html)
        self.assertEqual(len(catalogue.entries), len(library))
        self.assertNotIn('<details class="sl-python" open', html)
        for entry in catalogue.entries:
            namespace = {"library": library}
            exec(entry["recipes"]["0"], namespace)
            key, tensor = namespace["key"], namespace["tensor"]
            self.assertEqual(tensor[:], library[key][:])
            exec(entry["recipes"]["1"], namespace)
            expr = namespace["expr"]
            self.assertIsInstance(expr, TensorExpression)
            self.assertEqual(expr.rank, key.rank)
            self.assertNotIn("open_index", expr.to_typst())
            with contextlib.redirect_stdout(io.StringIO()) as output:
                exec(entry["recipes"]["2"], namespace)
            self.assertIn(expr.to_typst(), output.getvalue())
        self.assertEqual([key.structure for key in library], before)

    def test_registered_arguments_and_dual_representations_round_trip(self):
        rep = Representation("library_display::V", 2, is_self_dual=False)
        name = TensorName("library_display::A", print={"typst": "macron(A)"})
        x = S("library_display::x")
        library = TensorLibrary()
        signatures = [name(x + E("1/3"), n, rep.dual(), rep) for n in (7, 8)]
        for n, key in enumerate(signatures):
            library.register(Tensor.dense(key, [float(n), 2.0, 3.0, 4.0]))
        before = [library[key][:] for key in signatures]
        entries = Catalogue(library.to_html()).entries
        reconstructed = []
        for entry in entries:
            namespace = {"library": library}
            exec(entry["recipes"]["0"], namespace)
            reconstructed.append(namespace["key"].structure)
            self.assertEqual(namespace["tensor"][:], library[namespace["key"]][:])
            exec(entry["recipes"]["1"], namespace)
            self.assertEqual(namespace["expr"].rank, 2)
            self.assertIn("macron", namespace["expr"].to_typst())
        self.assertCountEqual(reconstructed, [key.structure for key in signatures])
        self.assertEqual([library[key][:] for key in signatures], before)

    def test_component_view_preserves_logical_order_and_sparse_storage(self):
        m, e = Representation.mink(2), Representation.euc(3)
        signature = TensorName("library_display::mixed")(m, e)
        tensor = Tensor.dense(signature, [float(i) for i in range(6)])
        library = TensorLibrary()
        library.register(tensor)
        entry = Catalogue(library.to_html()).entries[0]
        frame = entry["frames"][0]
        self.assertEqual(frame["sandbox"], "allow-scripts")
        payload = json.loads(
            re.search(
                r'<script id="tensor-data" type="application/json">(.*?)</script>',
                frame["srcdoc"],
                re.DOTALL,
            )[1]
        )
        self.assertEqual(payload["shape"], [2, 3])
        for indices, (_, plain, _) in payload["entries"]:
            self.assertEqual(float(plain), tensor[indices])
        sparse = Tensor.sparse(
            TensorName("library_display::sparse")(Representation.euc(1_000_000)),
            Expression,
        )
        sparse[999_999] = E("1/7")
        large = TensorLibrary()
        large.register(sparse)
        html = large.to_html()
        self.assertLess(len(html), 200_000)
        self.assertEqual(large[sparse.expression()][999_999], E("1/7"))

    def test_display_settings_and_independent_catalogues(self):
        library = TensorLibrary()
        key = TensorName("library_display::matrix")(
            Representation.euc(2), Representation.euc(2)
        )
        library.register(Tensor.dense(key, [1.0, 2.0, 3.0, 4.0]))
        first, second = library.to_html(), library.to_html()
        identity = r'<div data-spenso-library="(\d+)">'
        self.assertNotEqual(
            re.search(identity, first)[1], re.search(identity, second)[1]
        )
        static = library.to_html(settings=DisplaySettings(tensor_view="matrix"))
        self.assertNotIn("<iframe", static)
        self.assertIn("<mtable", static)

    def test_large_catalogue_is_bounded_without_hiding_mapping_entries(self):
        library = TensorLibrary()
        rep = Representation.euc(1)
        for n in range(40):
            key = TensorName("library_display::many")(n, rep)
            library.register(Tensor.dense(key, [float(n)]))
        html = library.to_html()
        self.assertLess(len(Catalogue(html).entries), 40)
        self.assertIn("of 40 tensors", html)
        self.assertEqual(len(library.keys()), 40)


if __name__ == "__main__":
    unittest.main()
