"""Linnest renders executable networks without evaluating their source formula."""

import json
import re
import unittest
import xml.etree.ElementTree as ET

from symbolica.core import Expression

E = Expression.parse
from symbolica.community.tensor import (
    DiagramRender,
    LayoutSettings,
    RenderSettings,
    StrokeStyle,
    Representation,
    Tensor,
    TensorExpression,
    TensorName,
    TensorNetwork,
)


def graph_labels(network):
    labels = {}
    for node in ET.fromstring(network.render().to_svg()).iter():
        if node.attrib.get("data-linnet-kind") not in {"node", "edge"}:
            continue
        detail = json.loads(node.attrib["data-linnet-detail"])
        if "label-typst" in detail:
            labels[(node.attrib["data-linnet-kind"], node.attrib["data-linnet-id"])] = (
                detail["label-typst"]
            )
    return list(labels.values())


class NetworkDisplayTests(unittest.TestCase):
    def test_native_svg_keeps_the_expression_hierarchy(self):
        rep = Representation.euc(2)
        A = TensorName("network_layout_tests::A")(rep, rep, rep)
        p = TensorName.vector("network_layout_tests::p")(rep)
        q = TensorName.vector("network_layout_tests::q")(rep)
        network = (A("i", "j", "k") * p("i") * q("j") / (p * p)).to_network()
        before = network.to_dot()
        root = ET.fromstring(network.render().to_svg())
        centers, edges = {}, {}
        for element in root.iter():
            kind = element.attrib.get("data-linnet-kind")
            if kind not in {"node", "edge"}:
                continue
            detail = json.loads(element.attrib["data-linnet-detail"])
            if kind == "edge":
                edges[detail["edge"]] = detail
                continue
            transform = re.fullmatch(
                r"translate\(([-\d.]+) ([-\d.]+)\)",
                element.attrib.get("transform", ""),
            )
            rect = element.find("{http://www.w3.org/2000/svg}rect")
            if transform and rect is not None:
                centers[detail["node"]] = (
                    float(transform[1]) + float(rect.attrib["width"]) / 2,
                    float(transform[2]) + float(rect.attrib["height"]) / 2,
                )
        dependencies = [
            edge for edge in edges.values() if edge["title"] == "Expression dependency"
        ]
        self.assertTrue(dependencies)
        for edge in dependencies:
            # SVG y increases downward: every input belongs below its operator.
            self.assertGreater(centers[edge["source"]][1], centers[edge["sink"]][1])
        # A contraction must not turn siblings into successive dependency ranks.
        children = {}
        for edge in dependencies:
            children.setdefault(edge["sink"], []).append(edge["source"])
        for siblings in children.values():
            for sibling in siblings[1:]:
                self.assertAlmostEqual(
                    centers[sibling][1], centers[siblings[0]][1], delta=0.001
                )
        self.assertEqual(network.to_dot(), before)

    def test_execution_summary_tracks_remaining_work_without_executing(self):
        rep = Representation.euc(2)
        tensor = Tensor.dense(
            TensorName("network_status_tests::A")(rep, rep), [1.0, 2.0, 3.0, 4.0]
        )
        network = tensor("i", "j") * tensor("j", "k") * tensor("k", "l")

        def summary(value):
            before = value.to_dot()
            status = value.status
            html = value.to_html()
            self.assertEqual(value.to_dot(), before)
            text = re.search(r'class="spenso-network-status".*?</div>', html).group()
            self.assertIn(f'data-complete="{str(status.complete).lower()}"', text)
            for count, noun in (
                (status.nodes, "node"),
                (status.operations, "operation"),
                (status.contractions, "contraction"),
            ):
                self.assertIn(f"{count} {noun}" + ("" if count == 1 else "s"), text)
            self.assertIn("Graph reduced" if status.complete else "Pending", text)
            return status

        initial = summary(network)
        stepped = network.step()
        partial = summary(stepped)
        self.assertLess(partial.contractions, initial.contractions)
        self.assertFalse(partial.complete)
        stepped.execute()
        self.assertTrue(summary(stepped).complete)
        self.assertFalse(summary(network).complete)
        # A single stored tensor can still have a pending self-contraction.
        self.assertFalse(summary(tensor.reindex("i", "i")).complete)
        # A library leaf is structurally reduced without being materialized by display.
        self.assertTrue(summary(TensorExpression.dirac_gamma(4).to_network()).complete)

    def test_edge_slots_use_representation_alphabets(self):
        gamma = TensorExpression.dirac_gamma(4)(1, 2, 1).to_network()
        labels = graph_labels(gamma)
        self.assertCountEqual(labels, ["a", "b", "mu", "gamma"])
        ET.fromstring(gamma.render().to_svg())

        # Compound graph indices must be named together, not independently as mu.
        rep = Representation.mink(4)
        p = TensorName("network_slot_labels::p")(rep)
        q = TensorName("network_slot_labels::q")(rep)
        network = (
            p(E("gammalooprs::hedge(2,1)")) * q(E("gammalooprs::hedge(3,1)"))
        ).to_network()
        labels = graph_labels(network)
        self.assertIn("mu", labels)
        self.assertIn("nu", labels)
        ET.fromstring(network.render().to_svg())

    def test_registered_names_and_typed_inspection(self):
        spinor = Representation.bis(4)
        jbar = TensorName("network_details::Jbar", print={"typst": "macron(J)"})(spinor)
        gamma = TensorExpression.dirac_gamma(4)
        network = (jbar(1) * gamma(1, 2, 1)).to_network()
        self.assertIn("macron(J)", graph_labels(network))
        root = ET.fromstring(network.render().to_svg())
        details = [
            json.loads(node.attrib["data-linnet-detail"])
            for node in root.iter()
            if "data-linnet-detail" in node.attrib
        ]
        titles = {detail.get("title") for detail in details}
        self.assertTrue(
            {
                "Product",
                "Stored tensor",
                "Library tensor",
                "Tensor contraction",
                "Free tensor slot",
                "Expression dependency",
                "Network output",
            }
            <= titles
        )
        tensor = next(
            detail
            for detail in details
            if dict(detail.get("properties", [])).get("Tensor name")
            == "network_details::Jbar"
        )
        properties = dict(tensor["properties"])
        self.assertEqual(properties["Rank"], "1")
        self.assertEqual(properties["Components"], "Symbolic")
        self.assertIn("Rank: 1", tensor["summary"])
        contraction = next(
            detail for detail in details if detail.get("title") == "Tensor contraction"
        )
        self.assertEqual(dict(contraction["properties"])["Dimension"], "4")
        self.assertTrue(
            any(
                "Stored tensor" in (node.text or "")
                for node in root.iter("{http://www.w3.org/2000/svg}title")
            )
        )

    def test_brackets_protect_open_occurrences_and_disappear_after_indexing(self):
        rep = Representation.mink(4)
        p = TensorName("bracket_slots::p")(rep)
        q = TensorName("bracket_slots::q")(rep)
        dyad = p.outer(p)
        self.assertEqual(dyad.rank, 2)
        self.assertIn("bracket(", str(dyad.to_expression()))
        indexed = dyad("mu", "nu")
        self.assertEqual(indexed.rank, 2)
        self.assertEqual(
            indexed.to_expression(), p("mu").to_expression() * p("nu").to_expression()
        )
        for left, right in [(p, q), (q, p)]:
            product = left.outer(right)
            self.assertEqual(
                product("mu", "nu").to_expression(),
                left("mu").to_expression() * right("nu").to_expression(),
            )
        imported = TensorExpression(
            TensorNetwork.bracket()(p("mu").to_expression(), q("nu").to_expression())
        )
        self.assertEqual(
            imported.to_expression(), p("mu").to_expression() * q("nu").to_expression()
        )
        summed = (p("mu") + q("mu")).outer(p("nu"))
        self.assertEqual(
            summed.to_expression(),
            (p("mu").to_expression() + q("mu").to_expression())
            * p("nu").to_expression(),
        )

    def test_graph_tracks_execution_while_source_expression_is_preserved(self):
        rep = Representation.euc(2)
        tensor = Tensor.dense(
            TensorName("network_display_tests::A")(rep, rep), [1.0, 2.0, 3.0, 4.0]
        )
        network = tensor("i", "j") * tensor("j", "k")
        expression = network.expression().to_expression()
        before = network.to_dot()
        html = network.to_html()
        self.assertRegex(before, r'tree = "[0-9A-Za-z]+"')
        # Numeric-leading tree IDs must not become extra DOT nodes.
        node_ids = set(
            re.findall(r'data-linnet-kind="node" data-linnet-id="(\d+)"', html)
        )
        self.assertEqual(node_ids, {"0", "1", "2"})
        self.assertIn("TensorNetwork", html)
        self.assertIn('data-linnet-kind="node"', html)
        self.assertIn('data-linnet-kind="edge"', html)
        self.assertEqual(network.to_dot(), before)
        network.execute()
        self.assertNotEqual(network.to_dot(), before)
        self.assertIn("<svg", network._repr_html_())
        self.assertEqual(network.expression().to_expression(), expression)
        self.assertEqual(list(network.result_tensor()), [7.0, 10.0, 15.0, 22.0])

    def test_scalar_sum_and_library_nodes_render(self):
        rep = Representation.mink(4)
        p, q = (
            TensorName(f"network_display_tests::{name}")(rep) for name in ("p", "q")
        )
        for network in (
            TensorNetwork(E("2")),
            (p("mu") + q("mu")).to_network(),
            (p("mu") * q("mu")).to_network(),
            TensorExpression.dirac_gamma(4)("a", "b", "mu").to_network(),
        ):
            with self.subTest(network=repr(network)):
                svg = network.render().to_svg()
                self.assertIn("data-linnet-interactive", svg)
                self.assertIn("spenso-network-svg", svg)
                self.assertIn("light-dark", svg)
                # Every emitted default paint must participate in theme switching,
                # including operator backgrounds and Linnest's arrowheads.
                root = ET.fromstring(svg)
                styles = "".join(root.itertext())
                for node in root.iter():
                    for paint in ("fill", "stroke"):
                        color = node.get(paint, "")
                        if color.startswith("#"):
                            self.assertIn(
                                f'[{paint}="{color}"]{{{paint}:light-dark(', styles
                            )

    def test_configured_snapshot_has_rich_display_and_is_stable(self):
        network = TensorNetwork(E("2 + x"))
        config = RenderSettings(
            title="Configured graph",
            layout=LayoutSettings(layout_algo="dot"),
            edge_stroke=StrokeStyle(paint="#123456", thickness=2),
        )
        drawing = network.render(config=config)
        self.assertIsInstance(drawing, DiagramRender)
        svg, html = drawing.to_svg(), drawing.to_html()
        self.assertIn("#123456", svg)
        self.assertIn("spenso-network", html)
        self.assertEqual(drawing._repr_svg_(), svg)
        self.assertEqual(drawing._repr_html_(), html)
        self.assertEqual(drawing._mime_(), ("text/html", html))
        self.assertIn("#image(bytes(", drawing.to_linnest())
        network.execute()
        self.assertEqual(drawing.to_html(), html)
        self.assertEqual(drawing.to_svg(), svg)
        self.assertIn("layout=LayoutSettings", repr(config))
        with self.assertRaises(AttributeError):
            config.title = "changed"
        with self.assertRaises(ValueError):
            RenderSettings(node_radius=-1)
        with self.assertRaises(TypeError):
            network.render(config={"title": "obsolete"})

    def test_linnet_configuration_and_portable_source(self):
        network = TensorExpression.dirac_gamma(4)("a", "b", "mu").to_network()
        config = RenderSettings(title="Network preview")
        source = network.to_linnest(config=config)
        self.assertIn("#image(bytes(", source)
        self.assertNotIn("@preview", source)
        self.assertNotIn("digraph", source)
        self.assertNotIn("network-dot", source)
        # Portable source contains the complete SVG, including typeset glyphs.
        self.assertIn("<svg", source)
        self.assertIn("<svg", network.render(config=config).to_svg())
        self.assertIn("<math", network.expression().to_html())


if __name__ == "__main__":
    unittest.main()
