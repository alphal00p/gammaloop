"""Linnest renders executable networks without evaluating their source formula."""

import re
import unittest

import linnet
from symbolica import E
from symbolica.community.spenso import (
    Representation,
    Tensor,
    TensorExpression,
    TensorName,
    TensorNetwork,
)


class NetworkDisplayTests(unittest.TestCase):
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
            TensorExpression.gamma(4)("a", "b", "mu").to_network(),
        ):
            with self.subTest(network=repr(network)):
                svg = network.render()
                self.assertIn("data-linnet-interactive", svg)
                self.assertIn("spenso-network-svg", svg)
                self.assertIn("light-dark", svg)

    def test_linnet_configuration_and_portable_source(self):
        network = TensorExpression.gamma(4)("a", "b", "mu").to_network()
        config = linnet.RenderConfig(title="Network preview")
        source = network.to_linnest(config=config)
        self.assertIn("Network preview", source)
        self.assertIn("render-network", source)
        self.assertNotIn("digraph", source)
        self.assertNotIn("network-dot", source)
        # The entrypoint contains its data; Linnet supplies the shared assets.
        prepared = linnet.PreparedRender.from_sources({"main.typ": source.encode()})
        self.assertIn("<svg", prepared.to_svg())
        self.assertIn("<svg", network.render(config=config))
        self.assertIn("<math", network.expression().to_html())


if __name__ == "__main__":
    unittest.main()
