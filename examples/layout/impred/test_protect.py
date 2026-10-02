"""Independent geometric stress tests of simultaneous protected movement."""

import json
import unittest
from pathlib import Path

import numpy as np
from init import initialize, validate_geometry
from solver import model_arrays, point_segments, protect, solve
from test_init import diagram


class ProtectionTests(unittest.TestCase):
    def test_existing_crossing_is_preserved_during_fixed_route_relaxation(self):
        state = {
            "positions": {
                "v:0": [-2, -1],
                "v:1": [2, 1],
                "v:2": [-2, 1],
                "v:3": [2, -1],
            },
            "node_ids": {str(i): f"v:{i}" for i in range(4)},
            "routes": {"0": ["v:0", "v:1"], "1": ["v:2", "v:3"]},
            "incoming_ids": [],
            "outgoing_ids": [],
            "report": {"planar": False},
        }
        result = solve(state, steps=100, flexible=False)
        self.assertEqual(len(result["report"]["geometry"]["crossings"]), 1)
        self.assertNotEqual(result["positions"], state["positions"])
        for snapshot in result["history"]:
            geometry = validate_geometry(snapshot["positions"], snapshot["routes"])
            self.assertEqual(len(geometry["crossings"]), 1)
            self.assertFalse(
                geometry["contacts"] or geometry["overlaps"] or geometry["degenerate"]
            )

    def test_arbitrary_proposals_preserve_all_old_separating_halfplanes(self):
        state = initialize(diagram(4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]))
        keys, _, p, _, segments, _, _ = model_arrays(state)
        random = np.random.default_rng(20260923)
        for _ in range(100):
            proposal = random.normal(size=p.shape) * random.uniform(0.01, 100)
            step, scales = protect(p, segments, proposal)
            _t, d, normal, valid = point_segments(p, segments)
            plane = np.einsum("nsi,ni->ns", normal, p) - d * 0.5
            moved_point = np.einsum("nsi,ni->ns", normal, p + step)
            self.assertTrue(np.all(moved_point[valid] >= plane[valid] - 1e-9))
            for endpoint in (0, 1):
                moved_endpoint = np.einsum(
                    "nsi,si->ns", normal, (p + step)[segments[:, endpoint]]
                )
                self.assertTrue(np.all(moved_endpoint[valid] <= plane[valid] + 1e-9))
            self.assertTrue(np.all(scales >= 0))
            self.assertTrue(np.all(scales <= 1))
            for fraction in (0.25, 0.5, 0.75, 1):
                positions = {
                    key: list(point) for key, point in zip(keys, p + fraction * step)
                }
                self.assertTrue(validate_geometry(positions, state["routes"])["valid"])

    def test_tiny_gap_does_not_tilt_halfplane_past_segment_endpoint(self):
        p = np.array(
            [[1.271456367, 3.781922547], [8.988221611, 5.545712136], [4.214, 0]]
        )
        vector = p[1] - p[0]
        perpendicular = np.array([-vector[1], vector[0]]) / np.linalg.norm(vector)
        p[2] = p[0] + 0.381 * vector + 1e-9 * perpendicular
        proposal = np.array([[0.1, 0.2], [-0.3, 0.1], [0.2, -0.7]])
        step, _ = protect(p, np.array([[0, 1]]), proposal)
        first, second = p[1] - p[0], p[2] - p[0]
        before = first[0] * second[1] - first[1] * second[0]
        after = p + step
        first, second = after[1] - after[0], after[2] - after[0]
        self.assertGreater((first[0] * second[1] - first[1] * second[0]) * before, 0)

    def test_adjacent_vertices_cannot_swap_or_collapse(self):
        p = np.array([[0.0, 0.0], [1.0, 0.0]])
        proposal = np.array([[10.0, 0.0], [-10.0, 0.0]])
        step, _ = protect(p, np.array([[0, 1]]), proposal)
        self.assertLess((p + step)[0, 0], (p + step)[1, 0])
        self.assertGreater(np.linalg.norm((p + step)[0] - (p + step)[1]), 1e-6)

    def test_b1_long_run_with_and_without_adaptive_bends(self):
        initial = initialize(
            json.loads(
                Path("/tmp/feynkit-dot-notebook/paper-b1/diagram.json").read_text()
            ),
            scale=2.4,
        )
        for flexible in (False, True):
            result = solve(initial, steps=1000, flexible=flexible)
            self.assertTrue(result["report"]["geometry"]["valid"])
            for epoch in result["history"]:
                self.assertTrue(
                    validate_geometry(epoch["positions"], epoch["routes"])["valid"]
                )

    def test_repeated_random_moves_keep_loops_parallel_edges_and_dangling_legs(self):
        initial = initialize(
            diagram(3, [(None, 0), (0, 1), (0, 1), (1, 1), (1, 2), (2, None)])
        )
        keys, _, p, _, segments, _, _ = model_arrays(initial)
        random = np.random.default_rng(1120)
        for _ in range(100):
            step, _ = protect(p, segments, random.normal(size=p.shape) * 0.15)
            for fraction in (0.5, 1):
                positions = {
                    key: list(point) for key, point in zip(keys, p + fraction * step)
                }
                self.assertTrue(
                    validate_geometry(positions, initial["routes"])["valid"]
                )
            p += step


if __name__ == "__main__":
    unittest.main(verbosity=2)
