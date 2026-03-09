import unittest

import linnet_py as lp


class TestLinnetPy(unittest.TestCase):
    def test_basic_graph_ops(self):
        dot = """
        digraph G {
          A;
          B;
          A -> B;
        }
        """
        g = lp.DotGraph.from_string(dot)

        print(g[lp.NodeIndex(0)].name)
        self.assertGreaterEqual(g.n_nodes(), 2)
        self.assertGreaterEqual(g.n_edges(), 1)
        self.assertGreaterEqual(g.n_hedges(), 2)

        sg = g.full_filter()
        edges = g.iter_edges_of(sg)
        nodes = g.iter_nodes_of(sg)

        self.assertGreaterEqual(len(edges), 1)
        self.assertGreaterEqual(len(nodes), 2)

        # Indexing: node -> vertex data
        node_idx, crown, node_data = nodes[0]
        fetched_node = g[node_idx]
        self.assertEqual(type(fetched_node), type(node_data))

        # Indexing: edge -> edge data
        edge_pair, edge_idx, edge_data = edges[0]
        fetched_edge = g[edge_idx]
        self.assertEqual(type(fetched_edge), type(edge_data.data))

        # Indexing: hedge -> hedge data (use first hedge from node neighbors)
        self.assertTrue(len(crown) > 0)
        fetched_hedge = g[crown[0]]
        self.assertTrue(hasattr(fetched_hedge, "statement"))

        comps = g.connected_components(sg)
        self.assertGreaterEqual(len(comps), 1)

    def test_builder(self):
        builder = lp.DotGraphBuilder()
        n0 = builder.add_node(lp.DotVertexData("A", None, None))
        n1 = builder.add_node(lp.DotVertexData("B", None, None))

        edge_data = lp.DotEdgeData({"label": "x"}, None, None)
        builder.add_edge(n0, n1, edge_data, lp.Orientation.default())

        g = builder.build()
        self.assertEqual(g.n_nodes(), 2)
        self.assertEqual(g.n_edges(), 1)

        sg = g.full_filter()
        edges = g.iter_edges_of(sg)
        self.assertEqual(len(edges), 1)


if __name__ == "__main__":
    unittest.main()
