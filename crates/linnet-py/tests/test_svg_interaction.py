"""Native SVG inspection preserves drawing and exports safe graph identities."""

import importlib.util
import json
import unittest
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import linnet as lp
import typst

SVG = "{http://www.w3.org/2000/svg}"
XLINK = "{http://www.w3.org/1999/xlink}"


class SvgInteractionTests(unittest.TestCase):
    def setUp(self):
        self.unsafe_node = 'left <script id="injected-node">alert(1)</script> & "x"'
        self.unsafe_edge = (
            'particle </script><script id="injected-edge">x</script> & "y"'
        )
        left = lp.node(self.unsafe_node, label="left")
        right = lp.node("right", label="right")
        self.graph = lp.build(
            left,
            right,
            lp.node("isolated"),
            lp.edge(
                lp.source(left),
                "propagator",
                lp.sink(right),
                label="p",
                style={
                    "stroke": lp.Stroke(
                        paint=lp.Color("#356a9a"), thickness=lp.Length.pt(0.8)
                    )
                },
                extensions={"particle": self.unsafe_edge},
            ),
            lp.edge(lp.source(right), "outgoing"),
            render_config=lp.RenderConfig(
                layouts=lp.LayoutOptions(seed=19, steps=20, label_steps=10)
            ),
        )

    def test_logical_ids_and_escaped_details_are_safe_and_keyboard_addressable(self):
        root = ET.fromstring(self.graph.to_svg())
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
                vertex = lp.node("v")
                endpoints = (
                    (lp.sink(vertex),)
                    if shape == "incoming"
                    else (lp.source(vertex),)
                    if shape == "outgoing"
                    else (lp.source(vertex), lp.sink(vertex))
                )
                graph = lp.build(
                    vertex,
                    lp.edge(endpoints[0], shape, *endpoints[1:]),
                )
                root = ET.fromstring(graph.to_svg())
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

    def test_native_drawing_geometry_and_paints_are_unchanged(self):
        compiled = []
        compile_native = typst.compile

        def capture_native(*args, **kwargs):
            result = compile_native(*args, **kwargs)
            compiled.append(result[0] if isinstance(result, list) else result)
            return result

        with patch.object(typst, "compile", side_effect=capture_native):
            decorated = self.graph.to_svg()
        native_root = ET.fromstring(compiled[0])
        decorated_root = ET.fromstring(decorated)
        self.assertTrue(
            any(
                element.get("stroke", "").lower() == "#356a9a"
                for element in native_root.iter()
            )
        )
        # Compare the complete native drawing after removing inspection-only
        # hit targets and scripts/styles from both documents.
        for root in (native_root, decorated_root):
            root.attrib.pop("data-linnet-interactive", None)
            for parent in list(root.iter()):
                for child in list(parent):
                    href = child.get("href", child.get(XLINK + "href", ""))
                    if (
                        child.tag in {SVG + "script", SVG + "style"}
                        or "data-linnet-kind" in child.attrib
                        or (child.tag == SVG + "a" and href.startswith("#linnet-"))
                    ):
                        parent.remove(child)
        self.assertEqual(ET.tostring(native_root), ET.tostring(decorated_root))

    def test_svg_file_and_rich_representations_share_the_prepared_drawing(self):
        prepared = self.graph.prepare_render()
        expected = prepared.to_svg()
        with TemporaryDirectory() as directory:
            output = Path(directory) / "graph.svg"
            self.assertEqual(prepared.render(output), output)
            self.assertEqual(output.read_text(), expected)
        self.assertEqual(self.graph._repr_html_(), self.graph._repr_svg_())
        self.assertIn('data-linnet-interactive="true"', self.graph._repr_html_())
        selected = self.graph.subgraph(edges=["propagator"])
        self.assertIn('data-linnet-interactive="true"', selected._repr_html_())

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
