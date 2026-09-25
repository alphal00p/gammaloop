"""Force regressions for flexible carriers and finite-segment applicability."""

import copy
import unittest
from itertools import pairwise

import numpy as np
from solver import (
    adapt,
    model_arrays,
    node_edge_forces,
    point_segments,
    segment_attraction_forces,
    solve,
)


class ForceTests(unittest.TestCase):
    def test_topology_scales_above_one_and_force_overflow(self):
        state = {
            "positions": {"v:0": [0.0, 0.0], "x:0": [-2.0, 0.0], "x:1": [2.0, 0.0]},
            "routes": {"0": ["v:0", "x:0"], "1": ["v:0", "x:1"]},
            "node_ids": {"0": "v:0"},
            "external_ids": {"0": "x:0", "1": "x:1"},
            "incoming_ids": ["x:0"],
            "outgoing_ids": ["x:1"],
            "report": {},
        }
        knobs = {
            "steps": 1,
            "target": 1.0,
            "flexible": False,
            "pull": 0.1,
            "repulsion": 0.0,
            "attraction": 0.0,
            "node_edge_strength": 0.0,
            "pull_scales": {"x:0": 0.5, "x:1": 2.0},
            "pull_balance": 2.0,
        }
        result = solve(state, **knobs)
        # With all other forces disabled, endpoint forces are -0.025 and +0.4;
        # their opposite reaction acts on the physical vertex.
        expected = {
            "v:0": [-0.375 * 0.045, 0],
            "x:0": [-2 - 0.025 * 0.045, 0],
            "x:1": [2 + 0.4 * 0.045, 0],
        }
        for key, position in expected.items():
            np.testing.assert_allclose(result["positions"][key], position, atol=1e-10)
        np.testing.assert_allclose(
            np.sum(list(result["positions"].values()), axis=0), [0, 0], atol=1e-10
        )
        for balance in (0.0, 1.0, 2.5):
            with self.subTest(balance=balance):
                solve(state, **{**knobs, "pull_balance": balance})
        with self.assertRaisesRegex(ValueError, "non-finite external force"):
            solve(state, **{**knobs, "pull_balance": 1e308})
        zero = solve(state, **{**knobs, "pull": 0.0, "pull_balance": 1e308})
        self.assertEqual(zero["positions"], state["positions"])
        for invalid in (0, -1, float("inf"), float("nan")):
            with (
                self.subTest(invalid=invalid),
                self.assertRaisesRegex(ValueError, "positive finite"),
            ):
                solve(state, **{**knobs, "pull_scales": {"x:0": 0.5, "x:1": invalid}})

    def test_internal_consecutive_and_external_connection_exemptions(self):
        state = {
            "positions": {
                "v:0": [0, 0],
                "b:0:a": [1, 1],
                "b:0:b": [2, 1],
                "b:0:c": [3, 1],
                "v:1": [4, 0],
                "v:2": [0, -2],
                "b:1:a": [5, 1],
                "b:1:b": [6, 1],
                "x:1": [7, 0],
            },
            "routes": {
                "0": ["v:0", "b:0:a", "b:0:b", "b:0:c", "v:1"],
                "1": ["v:1", "b:1:a", "b:1:b", "x:1"],
                "2": ["v:0", "v:2"],
            },
            "external_ids": {"1": "x:1"},
        }
        before = copy.deepcopy(state)
        _, ids, _, _, _, weights, _ = model_arrays(state, flexible_edges={"0", "1"})
        for route in (state["routes"]["0"], state["routes"]["1"]):
            for first, second in pairwise(route):
                self.assertEqual(weights[ids[first], ids[second]], 0)
                self.assertEqual(weights[ids[second], ids[first]], 0)
        self.assertAlmostEqual(weights[ids["b:0:a"], ids["b:0:c"]], (0.35 / 3) ** 2)
        self.assertAlmostEqual(weights[ids["b:0:a"], ids["b:1:a"]], 0.35 / 9)
        self.assertAlmostEqual(weights[ids["v:0"], ids["b:0:b"]], 0.35 / 3)
        external = [ids[key] for key in state["routes"]["1"]]
        np.testing.assert_array_equal(weights[np.ix_(external, external)], 0)
        self.assertEqual(weights[ids["v:0"], ids["v:1"]], 1)
        self.assertEqual(weights[ids["v:0"], ids["v:2"]], 1)
        np.testing.assert_array_equal(np.diag(weights), 0)
        self.assertEqual(state, before)

    def test_flexible_external_identity_persists_across_split_and_contraction(self):
        state = {
            "positions": {"v:0": [0, 0], "x:0": [4, 0], "v:probe": [0, 5]},
            "routes": {"0": ["v:0", "x:0"]},
            "external_ids": {"0": "x:0"},
            "flexible_edges": ["0"],
        }
        _, ids, _, _, _, weights, _ = model_arrays(state)
        self.assertEqual(weights[ids["v:0"], ids["x:0"]], 0)
        _, inserted, _ = adapt(state, 1, 0, max_points=1, remove=False)
        self.assertEqual(inserted, 1)
        _, ids, _, _, _, weights, _ = model_arrays(state)
        bend = state["routes"]["0"][1]
        self.assertEqual(weights[ids["v:0"], ids[bend]], 0)
        self.assertEqual(weights[ids[bend], ids["x:0"]], 0)
        self.assertEqual(weights[ids["v:0"], ids["x:0"]], 0)
        self.assertEqual(weights[ids[bend], ids["v:probe"]], 0.5)
        self.assertEqual(weights[ids["x:0"], ids["v:probe"]], 0.5)

        _, _, removed = adapt(state, 3, 0, max_points=0)
        self.assertEqual(removed, 1)
        _, ids, _, _, _, weights, _ = model_arrays(state)
        self.assertEqual(weights[ids["v:0"], ids["x:0"]], 0)
        self.assertEqual(weights[ids["x:0"], ids["v:probe"]], 1)

    def test_rigid_polyline_keeps_adjacent_repulsion(self):
        state = {
            "positions": {"v:0": [0, 0], "b:0:a": [1, 1], "v:1": [2, 0]},
            "routes": {"0": ["v:0", "b:0:a", "v:1"]},
        }
        _, ids, _, _, _, weights, _ = model_arrays(state)
        self.assertEqual(weights[ids["v:0"], ids["b:0:a"]], 0.35)
        self.assertEqual(weights[ids["v:1"], ids["b:0:a"]], 0.35)

    def test_internal_charge_budget_survives_actual_refinement(self):
        state = {
            "positions": {
                "v:0": [0, 0],
                "b:0:anchor": [4, 0],
                "v:1": [8, 0],
                "v:probe": [0, 8],
            },
            "routes": {"0": ["v:0", "b:0:anchor", "v:1"]},
            "anchor_ids": {"0": "b:0:anchor"},
            "flexible_edges": ["0"],
        }
        for phase in ("initial", "split", "contract"):
            with self.subTest(phase=phase):
                if phase == "split":
                    _, added, removed = adapt(state, 1, 0, max_points=3, remove=False)
                    self.assertEqual((added, removed), (2, 0))
                elif phase == "contract":
                    _, added, removed = adapt(state, 3, 0, max_points=0)
                    self.assertEqual((added, removed), (0, 2))
                _, ids, _, _, _, weights, _ = model_arrays(state)
                samples = [key for key in state["routes"]["0"] if key.startswith("b:")]
                observed = weights[[ids[key] for key in samples], ids["v:probe"]]
                np.testing.assert_allclose(observed, 0.35 / len(samples))
                self.assertAlmostEqual(float(observed.sum()), 0.35)
                self.assertIn("b:0:anchor", state["positions"])

    def test_loop_retains_nonadjacent_repulsion(self):
        route = ["v:0", "b:0:a", "b:0:anchor", "b:0:b", "v:0"]
        state = {
            "positions": dict(zip(route[:-1], ([0, 0], [1, 1], [2, 0], [1, -1]))),
            "routes": {"0": route},
            "anchor_ids": {"0": "b:0:anchor"},
            "flexible_edges": ["0"],
        }
        _, ids, _, _, _, weights, _ = model_arrays(state)
        for first, second in pairwise(route):
            self.assertEqual(weights[ids[first], ids[second]], 0)
        self.assertAlmostEqual(weights[ids["v:0"], ids["b:0:anchor"]], 0.35 / 3)
        self.assertAlmostEqual(weights[ids["b:0:a"], ids["b:0:b"]], (0.35 / 3) ** 2)

    def test_external_rigid_route_keeps_its_repulsion(self):
        state = {
            "positions": {"v:0": [0, 0], "b:0:a": [1, 1], "x:0": [2, 0]},
            "routes": {"0": ["v:0", "b:0:a", "x:0"]},
            "external_ids": {"0": "x:0"},
        }
        _, ids, _, _, _, weights, _ = model_arrays(state)
        self.assertEqual(weights[ids["v:0"], ids["x:0"]], 0.5)
        self.assertEqual(weights[ids["v:0"], ids["b:0:a"]], 0.5)
        self.assertEqual(weights[ids["b:0:a"], ids["x:0"]], 0.25)

    def test_subdivision_adds_serial_compliance(self):
        target, strength, total_length = 2.0, 1.7, 4.0
        for exponent in (1.0, 0.4):
            unsplit_tension = (
                strength * total_length * (total_length / target) ** exponent
            )
            for count in (1, 2, 3):
                for reverse in (False, True):
                    with self.subTest(exponent=exponent, count=count, reverse=reverse):
                        x = np.linspace(0, total_length, count + 1)
                        points = np.column_stack((x, np.zeros(count + 1)))
                        segments = np.array(list(pairwise(range(count + 1))))
                        if reverse:
                            segments = segments[:, ::-1]
                        knobs = {
                            "target": target,
                            "strength": strength,
                            "exponent": exponent,
                        }
                        force = segment_attraction_forces(points, segments, **knobs)
                        expected = np.zeros_like(points)
                        expected[0, 0] = unsplit_tension / count ** (exponent + 1)
                        expected[-1] = -expected[0]
                        np.testing.assert_allclose(force, expected, atol=1e-14)
                        # With the same external tension, n serial segments can
                        # support n times the total length of a single segment.
                        stretched = segment_attraction_forces(
                            points * count, segments, **knobs
                        )
                        expected[0, 0] = unsplit_tension
                        expected[-1] = -expected[0]
                        np.testing.assert_allclose(stretched, expected, atol=1e-14)

    def test_unequal_segments_move_knot_toward_equal_spacing(self):
        for exponent in (1.0, 0.4):
            for middle in (0.7, 3.1):
                with self.subTest(exponent=exponent, middle=middle):
                    points = np.array([[0.0, 0], [middle, 0], [4, 0]])
                    force = segment_attraction_forces(
                        points,
                        np.array([[0, 1], [1, 2]]),
                        target=2,
                        strength=1.7,
                        exponent=exponent,
                    )
                    self.assertGreater(force[1, 0] * (2 - middle), 0)
                    moved = middle + 0.01 * force[1, 0]
                    self.assertLess(abs(2 - moved), abs(2 - middle))
                    np.testing.assert_array_equal(force[:, 1], 0)
                    equal = segment_attraction_forces(
                        np.array([[0.0, 0], [2, 0], [4, 0]]),
                        np.array([[0, 1], [1, 2]]),
                        target=2,
                        strength=1.7,
                        exponent=exponent,
                    )
                    np.testing.assert_allclose(equal[1], [0, 0], atol=1e-14)

    def test_bent_route_uses_independent_segment_tensions(self):
        points = np.array([[0.0, 0], [2, 0], [2, 3]])
        for exponent in (1.0, 0.4):
            with self.subTest(exponent=exponent):
                force = segment_attraction_forces(
                    points,
                    np.array([[0, 1], [1, 2]]),
                    target=2,
                    strength=2.5,
                    exponent=exponent,
                )
                right = 2.5 * 3 * (3 / 2) ** exponent
                np.testing.assert_allclose(force, [[5, 0], [-5, right], [0, -right]])

    def test_attraction_is_negative_potential_gradient(self):
        points = np.array([[0.2, 0.5], [2.1, -0.4], [3.2, 1.7], [5.3, 0.8]])
        target, strength = 2.4, 1.7
        for exponent in (1.0, 0.4):
            for route in ([0, 1, 2, 3], [0, 1, 2, 0, 3]):
                with self.subTest(exponent=exponent, route=route):
                    segments = np.array(list(pairwise(route)))

                    def potential(p, exponent=exponent, route=route):
                        # Independent scalar oracle for a chain or a loop with
                        # a branch: every segment contributes its own energy.
                        return sum(
                            strength
                            * np.linalg.norm(p[b] - p[a]) ** (exponent + 2)
                            / ((exponent + 2) * target**exponent)
                            for a, b in pairwise(route)
                        )

                    force = segment_attraction_forces(
                        points,
                        segments,
                        target=target,
                        strength=strength,
                        exponent=exponent,
                    )
                    h = 1e-6
                    gradient = np.zeros_like(points)
                    for index in np.ndindex(points.shape):
                        offset = np.zeros_like(points)
                        offset[index] = h
                        gradient[index] = (
                            potential(points + offset) - potential(points - offset)
                        ) / (2 * h)
                    np.testing.assert_allclose(force, -gradient, atol=2e-8, rtol=2e-8)
                    np.testing.assert_allclose(force.sum(axis=0), [0, 0], atol=1e-14)
                    torque = points[:, 0] * force[:, 1] - points[:, 1] * force[:, 0]
                    self.assertAlmostEqual(float(torque.sum()), 0)

    def test_zero_length_segment_has_finite_zero_attraction(self):
        points = np.array([[1.0, 3], [1, 3], [4, 3]])
        for exponent in (1.0, 0.4):
            with self.subTest(exponent=exponent):
                knobs = {"target": 2, "strength": 2.5, "exponent": exponent}
                actual = segment_attraction_forces(
                    points,
                    np.array([[0, 1], [1, 2]]),
                    **knobs,
                )
                expected = segment_attraction_forces(
                    points, np.array([[1, 2]]), **knobs
                )
                np.testing.assert_array_equal(actual, expected)
                self.assertTrue(np.all(np.isfinite(actual)))

    def test_endpoint_projections_have_no_force_but_keep_guard_distances(self):
        segments = np.array([[0, 1]])
        for x in (-0.1, 0, 2, 2.1):
            with self.subTest(x=x):
                points = np.array([[0.0, 0], [2, 0], [x, 0.2]])
                data = point_segments(points, segments)
                t, distance, normal, valid = data
                self.assertIn(t[2, 0], (0, 1))
                self.assertTrue(valid[2, 0])
                self.assertGreater(distance[2, 0], 0)
                self.assertAlmostEqual(np.linalg.norm(normal[2, 0]), 1)
                forces = node_edge_forces(
                    segments, data, target=2, clearance=0.5, strength=3, exponent=2
                )
                np.testing.assert_array_equal(forces, np.zeros_like(points))
                # Force evaluation must not invalidate the finite-distance data
                # reused by crossing protection at exact/outside endpoints.
                self.assertTrue(valid[2, 0])
                self.assertGreater(distance[2, 0], 0)

    def test_interior_repulsion_has_balanced_endpoint_reactions(self):
        segments = np.array([[0, 1]])
        points = np.array([[0.0, 0], [2, 0], [0.5, 0.2]])
        data = point_segments(points, segments)
        forces = node_edge_forces(
            segments, data, target=2, clearance=0.5, strength=3, exponent=2
        )
        magnitude = 3 * 2 * (1 - 0.2) ** 2
        np.testing.assert_allclose(
            forces, [[0, -0.75 * magnitude], [0, -0.25 * magnitude], [0, magnitude]]
        )
        np.testing.assert_allclose(forces.sum(axis=0), [0, 0], atol=1e-14)
        rotation = np.array([[0.6, -0.8], [0.8, 0.6]])
        transformed = points @ rotation.T + [19, -37]
        actual = node_edge_forces(
            segments,
            point_segments(transformed, segments),
            target=2,
            clearance=0.5,
            strength=3,
            exponent=2,
        )
        np.testing.assert_allclose(actual, forces @ rotation.T, atol=1e-13)


class ParallelAttractionTests(unittest.TestCase):
    def test_physical_multiplicity_ignores_orientation_and_subdivision(self):
        routes = {
            "first": ["v:0", "b:first", "v:1"],
            "reverse": ["v:1", "b:reverse:0", "b:reverse:1", "v:0"],
            "bridge": ["v:1", "v:2"],
            "loop": ["v:2", "b:loop", "v:2"],
            "external": ["v:0", "b:external", "x:0"],
        }
        state = {
            "positions": {
                key: [float(i), 0.0]
                for i, key in enumerate(
                    dict.fromkeys(key for route in routes.values() for key in route)
                )
            },
            "routes": routes,
            "external_ids": {"external": "x:0"},
        }
        before = copy.deepcopy(state)
        baseline = model_arrays(state)
        balanced = model_arrays(state, parallel_attraction_balance=1.0)
        np.testing.assert_array_equal(baseline[-1], np.ones(10))
        np.testing.assert_array_equal(balanced[-1], [0.5] * 5 + [1.0] * 5)
        np.testing.assert_array_equal(balanced[5], baseline[5])
        self.assertEqual(state, before)
        state["routes"]["first"].insert(1, "b:extra")
        state["positions"]["b:extra"] = [0.3, 0.2]
        split = model_arrays(state, parallel_attraction_balance=1.0)
        np.testing.assert_array_equal(split[-1], [0.5] * 6 + [1.0] * 5)

    def test_full_balance_compensates_congruent_parallel_attractions(self):
        state = {
            "positions": {"v:0": [0.0, 0.0], "v:1": [2.0, 0.0]},
            "routes": {"a": ["v:0", "v:1"], "b": ["v:1", "v:0"]},
        }
        _, _, positions, _, segments, _, scales = model_arrays(
            state, parallel_attraction_balance=1.0
        )
        single = segment_attraction_forces(
            positions, segments[:1], target=1.3, strength=2.5, exponent=0.4
        )
        balanced = segment_attraction_forces(
            positions, segments, target=1.3, strength=2.5 * scales, exponent=0.4
        )
        original = segment_attraction_forces(
            positions, segments, target=1.3, strength=2.5, exponent=0.4
        )
        np.testing.assert_array_equal(balanced, single)
        np.testing.assert_array_equal(original, 2 * single)
        np.testing.assert_array_equal(balanced.sum(axis=0), [0, 0])

    def test_invalid_parallel_balance_is_rejected_before_solving(self):
        for balance in (-1, float("inf"), float("nan")):
            with (
                self.subTest(balance=balance),
                self.assertRaisesRegex(ValueError, "finite and nonnegative"),
            ):
                solve({}, parallel_attraction_balance=balance)


if __name__ == "__main__":
    unittest.main()
