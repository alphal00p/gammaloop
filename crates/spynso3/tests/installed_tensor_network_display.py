"""Linnest renders executable networks without evaluating their source formula."""

import re
import unittest

import linnet
from symbolica import E
from symbolica.community.spenso import (
    Representation,
    Tensor,
    TensorExpression,
    TensorLibrary,
    TensorName,
    TensorNetwork,
)


class NetworkDisplayTests(unittest.TestCase):
    def test_spinor_current_is_nary_before_execution(self):
        spinor = Representation.bis(4)
        jbar = TensorName("nary_current::Jbar")(spinor)
        j = TensorName("nary_current::J")(spinor)
        gamma = TensorExpression.gamma(4)
        factors = [jbar(1), gamma(1, 2, 1), gamma(2, 3, 2), gamma(3, 4, 3), j(4)]
        current = factors[0]
        nested = current.to_expression()
        for factor in factors[1:]:
            current = current * factor
            nested = TensorNetwork.bracket()(nested, factor.to_expression())
        library = TensorLibrary.hep_lib_atom()
        network = current.to_network(library=library)
        self.assertEqual(network.to_dot().count('label = "∏"'), 1)
        self.assertEqual(
            list(current.to_expression()), [f.to_expression() for f in factors]
        )
        reference = TensorExpression(nested).to_network(library=library)
        network.execute(library=library)
        reference.execute(library=library)
        actual = list(network.result_tensor(library=library))
        expected = list(reference.result_tensor(library=library))
        self.assertEqual(len(actual), 64)
        self.assertTrue(
            all((a - b).expand() == E("0") for a, b in zip(actual, expected))
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
