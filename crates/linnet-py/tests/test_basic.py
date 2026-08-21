import ast
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

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

    def test_from_file_accepts_pathlike(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "graph.dot"
            path.write_text("digraph G { A -> B; }", encoding="utf-8")

            graph = lp.DotGraph.from_file(path)

        self.assertEqual(graph.n_nodes(), 2)
        self.assertEqual(graph.n_edges(), 1)

    def test_typst_rendering(self):
        graph = lp.DotGraph.from_string(
            'digraph example { A -> B [label="p_0"]; B -> C [label="p_1"]; }'
        )

        with TemporaryDirectory() as directory:
            template = Path(directory) / "custom.typ"
            template.write_text(
                """
                #set page(width: auto, height: auto, margin: 2mm)
                #let source = read(sys.inputs.at("data-path"))
                #let marker = sys.inputs.at("marker")
                [#marker]
                #raw(source, lang: "dot")
                """,
                encoding="utf-8",
            )
            output = Path(directory) / "diagram.pdf"
            self.assertEqual(graph.render(output), output)
            self.assertTrue(output.read_bytes().startswith(b"%PDF-"))

            custom_output = Path(directory) / "custom.svg"
            graph.render(
                custom_output,
                template=template,
                inputs=types.MappingProxyType({"marker": "custom template"}),
            )
            self.assertIn("<svg", custom_output.read_text(encoding="utf-8"))

        svg = graph.to_svg()
        self.assertIn("<svg", svg)
        self.assertIn("</svg>", svg)
        self.assertIn("<svg", graph._repr_svg_())

    def test_indices_are_values(self):
        self.assertEqual(lp.Hedge(3), lp.Hedge(3))
        self.assertNotEqual(lp.Hedge(3), lp.Hedge(4))
        self.assertEqual(lp.NodeIndex(2), lp.NodeIndex(2))
        self.assertNotEqual(lp.NodeIndex(2), lp.NodeIndex(3))
        self.assertEqual(lp.EdgeIndex(1), lp.EdgeIndex(1))
        self.assertNotEqual(lp.EdgeIndex(1), lp.EdgeIndex(2))

        values = {
            lp.Hedge(3): "hedge",
            lp.NodeIndex(2): "node",
            lp.EdgeIndex(1): "edge",
        }
        self.assertEqual(values[lp.Hedge(3)], "hedge")
        self.assertEqual(values[lp.NodeIndex(2)], "node")
        self.assertEqual(values[lp.EdgeIndex(1)], "edge")

    def test_builder(self):
        builder = lp.DotGraphBuilder()
        n0 = builder.add_node("A", types.MappingProxyType({"shape": "circle"}), index=0)
        n1 = builder.add_node("B", {"color": "blue"})
        builder.add_edge(
            n0,
            n1,
            {"label": "x"},
            "reversed",
            local_statements=types.MappingProxyType({"weight": "2"}),
            edge_id=0,
            source_hedge=lp.DotHedgeData(port_label="out"),
        )
        builder.add_external_edge(
            n0,
            {"label": "incoming"},
            lp.Orientation.undirected(),
            lp.Flow.sink(),
            hedge=lp.DotHedgeData(statement="external"),
        )

        g = builder.build()
        self.assertEqual(g.n_nodes(), 2)
        self.assertEqual(g.n_edges(), 2)
        self.assertEqual(g[n0].name, "A")
        self.assertEqual(g[n0].index, lp.NodeIndex(0))
        self.assertEqual(g[n0].statements["shape"], "circle")
        self.assertEqual(g[n1].statements["color"], "blue")

        sg = g.full_filter()
        edges = g.iter_edges_of(sg)
        self.assertEqual(len(edges), 2)
        labels = {edge.data.statements["label"] for _, _, edge in edges}
        self.assertEqual(labels, {"x", "incoming"})

        paired = next(edge for edge in edges if edge[0].kind == "paired")
        external = next(edge for edge in edges if edge[0].kind == "unpaired")
        self.assertEqual(repr(paired[2].orientation), "Orientation.Reversed")
        self.assertEqual(paired[2].data.edge_id, lp.EdgeIndex(0))
        self.assertEqual(paired[2].data.local_statements["weight"], "2")
        self.assertEqual(repr(external[2].orientation), "Orientation.Undirected")
        self.assertEqual(repr(external[0].flow), "Flow.Sink")


if __name__ == "__main__":
    unittest.main()
