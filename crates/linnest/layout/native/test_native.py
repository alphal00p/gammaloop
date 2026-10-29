"""Check native decomposition identities and rotations independently of EC logic."""

import json
import os
import subprocess
import unittest
from pathlib import Path


class NativeDecompositionTests(unittest.TestCase):
    binary = Path(
        os.environ.get(
            "EC_SPQR_BINARY",
            Path(__file__).resolve().parents[4] / "target/ec-planarity-native/ec-spqr",
        )
    )

    @classmethod
    def request(cls, nodes, pairs, root_edge=None):
        request = {
            "nodes": nodes,
            "edges": [
                {"id": str(i), "source": a, "target": b}
                for i, (a, b) in enumerate(pairs)
            ],
        }
        if root_edge is not None:
            request["root_edge"] = root_edge
        result = subprocess.run(
            [str(cls.binary)],
            input=json.dumps(request) + "\n",
            text=True,
            capture_output=True,
            check=True,
            timeout=30,
        )
        return json.loads(result.stdout), request

    def assert_plane(self, nodes, edges, rotations):
        endpoints = {e["id"]: (e["source"], e["target"]) for e in edges}
        darts = {
            (node, edge) for node, rotation in rotations.items() for edge in rotation
        }
        self.assertEqual(len(darts), 2 * len(edges))
        self.assertEqual(set(rotations), set(nodes))
        remaining = darts.copy()
        faces = 0
        while remaining:
            start = next(iter(remaining))
            dart = start
            while True:
                self.assertIn(dart, remaining)
                remaining.remove(dart)
                node, edge = dart
                a, b = endpoints[edge]
                self.assertIn(node, (a, b))
                other = b if node == a else a
                rotation = rotations[other]
                dart = (other, rotation[(rotation.index(edge) + 1) % len(rotation)])
                if dart == start:
                    break
            faces += 1
        self.assertEqual(len(nodes) - len(edges) + faces, 2)

    def test_mixed_spr_decomposition_and_reconstruction(self):
        # K4 with one extra two-edge path across ab has S, P, and R skeletons.
        pairs = [
            ("a", "b"),
            ("a", "c"),
            ("a", "d"),
            ("b", "c"),
            ("b", "d"),
            ("c", "d"),
            ("a", "x"),
            ("x", "b"),
        ]
        output, request = self.request(["a", "b", "c", "d", "x"], pairs, "0")
        self.assertTrue(output["planar"])
        self.assertTrue(output["biconnected"])
        self.assertEqual({c["type"] for c in output["components"]}, {"S", "P", "R"})
        components = {c["id"]: c for c in output["components"]}
        edges = {e["id"]: e for c in components.values() for e in c["edges"]}
        self.assertEqual(
            sorted(
                e["real_edge"] for e in edges.values() if e["real_edge"] is not None
            ),
            sorted(e["id"] for e in request["edges"]),
        )
        for component in components.values():
            self.assert_plane(
                component["nodes"], component["edges"], component["rotation_cw"]
            )
            for edge in component["edges"]:
                if edge["twin"] is not None:
                    twin = edges[edge["twin"]["edge"]]
                    self.assertEqual(
                        twin["twin"], {"component": component["id"], "edge": edge["id"]}
                    )
                    self.assertEqual(
                        {edge["source"], edge["target"]},
                        {twin["source"], twin["target"]},
                    )

        def expand(component_id, parent_edge=None):
            component = components[component_id]
            rotations = {n: list(r) for n, r in component["rotation_cw"].items()}
            for edge in component["edges"]:
                twin = edge["twin"]
                if twin is None or edge["id"] == parent_edge:
                    continue
                child = expand(twin["component"], twin["edge"])
                for node, rotation in child.items():
                    if twin["edge"] in rotation:
                        pivot = rotation.index(twin["edge"])
                        replacement = rotation[pivot + 1 :] + rotation[:pivot]
                        pivot = rotations[node].index(edge["id"])
                        rotations[node][pivot : pivot + 1] = replacement
                    else:
                        self.assertNotIn(node, rotations)
                        rotations[node] = rotation
            return rotations

        rebuilt = {
            node: [edges[edge]["real_edge"] for edge in rotation]
            for node, rotation in expand(output["root"]).items()
        }
        for node, rotation in rebuilt.items():
            expected = output["rotation_cw"][node]
            pivot = rotation.index(expected[0])
            self.assertEqual(rotation[pivot:] + rotation[:pivot], expected)
        self.assert_plane(request["nodes"], request["edges"], rebuilt)

    def test_parallel_edges_and_root_identity(self):
        output, request = self.request(["u", "v"], [("u", "v")] * 4, "2")
        self.assertEqual(len(output["components"]), 1)
        component = output["components"][0]
        self.assertEqual(component["type"], "P")
        reference = next(
            e for e in component["edges"] if e["id"] == component["reference_edge"]
        )
        self.assertEqual(reference["real_edge"], "2")
        self.assert_plane(request["nodes"], request["edges"], output["rotation_cw"])

    def test_nonplanarity_is_explicit(self):
        output, _ = self.request(list("abcdef"), [(a, b) for a in "abc" for b in "def"])
        self.assertFalse(output["planar"])
        self.assertEqual(output["components"], [])

    def test_protocol_recovers_after_invalid_input(self):
        requests = [
            "{bad json",
            json.dumps({"nodes": ["a", "a"], "edges": []}),
            json.dumps(
                {
                    "nodes": ["a", "b", "c"],
                    "edges": [
                        {"id": "ab", "source": "a", "target": "b"},
                        {"id": "bc", "source": "b", "target": "c"},
                    ],
                }
            ),
        ]
        result = subprocess.run(
            [str(self.binary)],
            input="\n".join(requests) + "\n",
            text=True,
            capture_output=True,
            check=True,
            timeout=30,
        )
        responses = [json.loads(line) for line in result.stdout.splitlines()]
        self.assertEqual(len(responses), 3)
        self.assertIn("error", responses[0])
        self.assertIn("error", responses[1])
        self.assertTrue(responses[2]["planar"])
        self.assertFalse(responses[2]["biconnected"])


if __name__ == "__main__":
    unittest.main()
