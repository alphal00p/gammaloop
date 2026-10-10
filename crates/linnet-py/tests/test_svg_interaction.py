"""Native SVG inspection preserves drawing and exports safe graph identities."""

import importlib.util
import json
import unittest
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

from linnet_under_test import graph as lp

SVG = "{http://www.w3.org/2000/svg}"
XLINK = "{http://www.w3.org/1999/xlink}"


class SvgInteractionTests(unittest.TestCase):
    def setUp(self):
        self.unsafe_node = 'left <script id="injected-node">alert(1)</script> & "x"'
        self.unsafe_edge = (
            'particle </script><script id="injected-edge">x</script> & "y"'
        )
        left = lp.Node(self.unsafe_node, label="left")
        right = lp.Node("right", label="right")
        self.graph = lp.Graph(
            nodes=[left, right, lp.Node("isolated")],
            edges=[
                lp.Edge(
                    source=left,
                    sink=right,
                    name="propagator",
                    label="p",
                    style={
                        "stroke": lp.Stroke(
                            paint=lp.Color("#356a9a"), thickness=lp.Length.pt(0.8)
                        )
                    },
                    extensions={"particle": self.unsafe_edge},
                ),
                lp.Edge(source=right, name="outgoing"),
            ],
            render_config=lp.RenderSettings(
                layouts=lp.LayoutSettings(seed=19, steps=20, label_steps=10)
            ),
        )

    def test_logical_ids_and_escaped_details_are_safe_and_keyboard_addressable(self):
        root = ET.fromstring(self.graph.to_svg(interactive=True))
        targets = root.findall(".//*[@data-linnet-kind]")
        identities = {
            (target.attrib["data-linnet-kind"], int(target.attrib["data-linnet-id"]))
            for target in targets
        }
        self.assertEqual(
            identities,
            {("node", 0), ("node", 1), ("node", 2), ("edge", 0), ("edge", 1)}
            | {("halfedge", half.index) for half in self.graph.half_edges()},
        )
        focusable = Counter(
            (target.attrib["data-linnet-kind"], int(target.attrib["data-linnet-id"]))
            for target in targets
            if target.attrib["tabindex"] == "0"
        )
        self.assertEqual(focusable, Counter(dict.fromkeys(identities, 1)))
        self.assertGreater(
            sum(target.attrib["data-linnet-kind"] == "edge" for target in targets),
            self.graph.n_edges,
        )
        for target in targets:
            self.assertEqual(target.attrib["role"], "button")
            self.assertEqual(target.attrib["aria-pressed"], "false")
            self.assertIsNone(target.get("href"))
            self.assertIsNone(target.get(XLINK + "href"))
            self.assertEqual(
                target.find(SVG + "title").text, target.attrib["aria-label"]
            )
            detail = json.loads(target.attrib["data-linnet-detail"])
            if (
                target.attrib["data-linnet-kind"] == "node"
                and target.attrib["data-linnet-id"] == "0"
            ):
                self.assertEqual(detail["name"], self.unsafe_node)
            if (
                target.attrib["data-linnet-kind"] == "edge"
                and target.attrib["data-linnet-id"] == "0"
            ):
                self.assertEqual(detail["particle"], self.unsafe_edge)
            if target.attrib["data-linnet-kind"] == "halfedge":
                half = self.graph.half_edge(int(target.attrib["data-linnet-id"]))
                self.assertEqual(detail["half-edge"], half.index)
                self.assertEqual(detail["node"], half.node.index)
                self.assertEqual(detail["edge"], half.edge.index)
                self.assertEqual(
                    detail["pair"], None if half.pair is None else half.pair.index
                )
        self.assertEqual(len(root.findall(".//" + SVG + "script")), 1)
        self.assertIsNone(root.find(".//*[@id='injected-node']"))
        self.assertIsNone(root.find(".//*[@id='injected-edge']"))

    def test_dangling_and_self_loop_half_edges_keep_native_identity(self):
        for shape in ("incoming", "outgoing", "self_loop"):
            with self.subTest(shape=shape):
                vertex = lp.Node("v")
                graph = lp.Graph(
                    nodes=[vertex],
                    edges=[
                        lp.Edge(
                            source=None if shape == "incoming" else vertex,
                            sink=None if shape == "outgoing" else vertex,
                            name=shape,
                        )
                    ],
                )
                root = ET.fromstring(graph.to_svg(interactive=True))
                targets = root.findall(".//*[@data-linnet-kind='halfedge']")
                self.assertEqual(
                    {int(target.get("data-linnet-id")) for target in targets},
                    {half.index for half in graph.half_edges()},
                )
                for target in targets:
                    half = graph.half_edge(int(target.get("data-linnet-id")))
                    detail = json.loads(target.get("data-linnet-detail"))
                    self.assertEqual(detail["edge"], half.edge.index)
                    self.assertEqual(detail["node"], half.node.index)
                    self.assertEqual(
                        detail["flow"],
                        "source" if half.flow == lp.Flow.Source else "sink",
                    )
                    if half.pair is None:
                        self.assertIsNone(detail["pair"])
                        self.assertEqual(
                            graph.subgraph(half_edges=[half.index]).half_edge_indices(),
                            graph.subgraph(edges=[half.edge.index]).half_edge_indices(),
                        )
                    else:
                        self.assertEqual(detail["pair"], half.pair.index)
                        self.assertNotEqual(detail["pair"], half.index)

    def test_native_drawing_contains_custom_paints_and_inspection_targets(self):
        for interactive in (False, True):
            root = ET.fromstring(self.graph.to_svg(interactive=interactive))
            self.assertTrue(
                any(
                    element.get("stroke", "").lower() == "#356a9a"
                    for element in root.iter()
                )
            )
            # Only interactive output carries inspection targets and their script.
            self.assertEqual(bool(root.findall(".//*[@data-linnet-kind]")), interactive)
            self.assertEqual(bool(root.findall(".//" + SVG + "a")), interactive)
            self.assertEqual(bool(root.findall(".//" + SVG + "script")), interactive)

    def assert_native_drawing_unchanged(self, drawing, interactive):
        # Also called by the Rust integration test at the embedded compiler boundary.
        drawing_root = ET.fromstring(drawing)
        interactive_root = ET.fromstring(interactive)
        self.assertFalse(drawing_root.findall(".//" + SVG + "a"))
        self.assertTrue(interactive_root.findall(".//*[@data-linnet-kind]"))

        # Compare the complete drawings after removing inspection-only hit
        # targets, their scripts/styles, and the groups Typst places them in.
        def remove_targets(parent):
            for child in list(parent):
                href = child.get("href", child.get(XLINK + "href", ""))
                if (
                    child.tag in {SVG + "script", SVG + "style"}
                    or "data-linnet-kind" in child.attrib
                    or (child.tag == SVG + "a" and href.startswith("#linnet-"))
                ):
                    parent.remove(child)
                    continue
                remove_targets(child)
                if child.tag == SVG + "g" and len(child) == 0:
                    parent.remove(child)

        for root in (drawing_root, interactive_root):
            root.attrib.pop("data-linnet-interactive", None)
            remove_targets(root)
        self.assertEqual(ET.tostring(drawing_root), ET.tostring(interactive_root))

    def test_static_and_interactive_output_share_the_drawing(self):
        typst_route = lp.RenderSettings(
            layouts=lp.LayoutSettings(impred_steps=1).then(impred_steps=1)
        )
        for config in (None, typst_route):
            with self.subTest(typst=config is not None):
                result = self.graph.render(config=config)
                self.assert_native_drawing_unchanged(
                    result.to_svg(), result.to_svg(interactive=True)
                )

    def test_svg_file_and_rich_representations_share_the_prepared_drawing(self):
        prepared = self.graph.render()
        expected = prepared.to_svg()
        with TemporaryDirectory() as directory:
            output = Path(directory) / "graph.svg"
            self.assertEqual(prepared.save(output), output)
            self.assertEqual(output.read_text(), expected)
        # HTML display is interactive; the SVG image representation is static.
        self.assertEqual(self.graph._repr_svg_(), expected)
        self.assertEqual(self.graph._repr_html_(), prepared.to_svg(interactive=True))
        self.assertIn('data-linnet-interactive="true"', self.graph._repr_html_())
        self.assertNotIn("data-linnet-interactive", self.graph._repr_svg_())
        selected = self.graph.subgraph(edges=["propagator"])
        self.assertIn('data-linnet-interactive="true"', selected._repr_html_())
        self.assertNotIn("data-linnet-interactive", selected._repr_svg_())

    @unittest.skipUnless(importlib.util.find_spec("marimo"), "Marimo is not installed")
    def test_marimo_prefers_scripted_html_and_wraps_it_for_execution(self):
        from marimo._output.formatting import try_format

        for value in (self.graph, self.graph.subgraph(edges=["propagator"])):
            with self.subTest(value=type(value).__name__):
                output = try_format(value)
                self.assertIsNone(output.exception)
                self.assertEqual(output.mimetype, "text/html")
                self.assertIn("<iframe", output.data)
                self.assertIn("data-linnet-interactive", output.data)


if __name__ == "__main__":
    unittest.main()
