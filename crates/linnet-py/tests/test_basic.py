import ast
import types
import unittest
from pathlib import Path

import linnet_py as lp


class TestLinnetPy(unittest.TestCase):
    def test_documented_surface_matches_runtime(self):
        stub_path = Path(__file__).resolve().parents[1] / "linnet_py.pyi"
        stub_text = stub_path.read_text(encoding="utf-8")
        docs_stub_path = (
            Path(__file__).resolve().parents[3]
            / "docs"
            / "api"
            / "python"
            / "linnet-py.pyi"
        )
        self.assertEqual(stub_text, docs_stub_path.read_text(encoding="utf-8"))
        stub = ast.parse(stub_text)
        export_assignments = [
            node
            for node in stub.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "__all__"
                for target in node.targets
            )
        ]
        self.assertEqual(
            len(export_assignments), 1, "stub must define one __all__ list"
        )
        documented = ast.literal_eval(export_assignments[0].value)
        self.assertIsInstance(documented, list)
        self.assertTrue(all(isinstance(name, str) for name in documented))

        self.assertCountEqual(lp.__all__, documented)
        runtime = {
            name
            for name in dir(lp)
            if not name.startswith("_")
            and not isinstance(getattr(lp, name), types.ModuleType)
        }
        self.assertCountEqual(runtime, documented)

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
