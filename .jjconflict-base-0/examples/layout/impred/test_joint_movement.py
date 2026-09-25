"""Movement regressions for independent axes and raw-force group aggregation."""

import json
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from init import validate_geometry
from movement import MovementDOFs, _project_component, project_halfplanes
from solver import model_arrays, protect, separating_halfplanes


def records(p, x_groups=(), fixed=(), shifts=None):
    result = [
        {
            "kind": "edge",
            "id": index,
            "key": str(index),
            "initial_position": point.tolist(),
            "shift": None,
            "fixed_axes": [False, False],
            "constraints": {"x": "Free", "y": "Free"},
        }
        for index, point in enumerate(p)
    ]
    for group in x_groups:
        for index in group:
            result[index]["constraints"]["x"] = {"Grouped": [{"Edge": group[0]}, "Any"]}
    for index, axis in fixed:
        result[index]["fixed_axes"][axis] = True
        result[index]["constraints"]["xy"[axis]] = "Fixed"
    for index, shift in (shifts or {}).items():
        result[index]["shift"] = shift
    return result


class JointMovementTests(unittest.TestCase):
    def test_blocked_y_does_not_freeze_shared_x_or_other_free_y(self):
        p = np.array([[0.0, 0.0], [4.0, 0.0], [2.0, 1e-6], [2.0, 3.0]])
        metadata = records(
            p, x_groups=[(2, 3)], fixed=[(i, a) for i in (0, 1) for a in (0, 1)]
        )
        dofs = MovementDOFs(list(map(str, range(4))), metadata)
        proposal = np.array([[0.0, 0.0], [0.0, 0.0], [0.2, -0.2], [0.2, 0.2]])
        evidence = {}
        step, _ = protect(
            p, np.array([[0, 1]]), proposal, dofs=dofs, cap=0.2, certificate=evidence
        )
        self.assertGreater(step[2, 0], 0.1)
        self.assertEqual(step[2, 0], step[3, 0])
        self.assertGreater(step[3, 1], 0.05)
        self.assertGreaterEqual(step[2, 1], -1e-10)
        np.testing.assert_array_equal(step[:2], np.zeros((2, 2)))
        self.assertLess(evidence["separator_residual"], 1e-10)
        for fraction in (0.25, 0.5, 1):
            after = {
                str(i): point.tolist() for i, point in enumerate(p + fraction * step)
            }
            self.assertTrue(validate_geometry(after, {"edge": ["0", "1"]})["valid"])

    def test_group_averages_raw_force_before_trust_cap(self):
        p = np.array([[2.0, -2.0], [2.0, 2.0]])
        dofs = MovementDOFs(["0", "1"], records(p, x_groups=[(0, 1)]))
        proposal = np.array([[10.0, 1.0], [-1.0, 0.0]])
        step, _ = protect(p, np.empty((0, 2), dtype=int), proposal, dofs=dofs, cap=0.1)
        self.assertEqual(step[0, 0], step[1, 0])
        self.assertGreater(step[0, 0], 0.095)
        self.assertLessEqual(float(np.linalg.norm(step, axis=1).max()), 0.1 + 1e-10)
        # Averaging the old per-point saturated suggestions nearly cancels X.
        old = (
            proposal * np.minimum(1.0, 0.1 / np.linalg.norm(proposal, axis=1))[:, None]
        )
        self.assertLess(abs(float(old[:, 0].mean())), 0.001)

    def test_shifted_cycle_and_fixed_reference_elimination(self):
        p = np.array([[2.0, 0.0], [3.0, 2.0], [4.0, 4.0]])
        metadata = records(
            p,
            fixed=[(1, 1)],
            shifts={
                0: {"x": 1.0, "y": 0.0},
                1: {"x": 2.0, "y": 0.0},
                2: {"x": 3.0, "y": 0.0},
            },
        )
        for i in range(3):
            metadata[i]["constraints"]["x"] = {
                "Grouped": [{"Edge": (i + 1) % 3}, "PositiveOnly"]
            }
        dofs = MovementDOFs(["0", "1", "2"], metadata)
        step, _ = protect(
            p,
            np.empty((0, 2), dtype=int),
            np.array([[-100.0, 2.0], [-1.0, 4.0], [-1.0, -2.0]]),
            dofs=dofs,
            cap=2.0,
        )
        np.testing.assert_array_equal(step[:, 0], np.repeat(step[0, 0], 3))
        self.assertEqual(step[1, 1], 0.0)
        self.assertGreater(p[0, 0] + step[0, 0] - 1.0, 0.0)
        dofs.validate(p + step)
        # References are immutable even after a successful movement.
        invalid = p + step
        invalid[1, 1] += 0.1
        with self.assertRaisesRegex(ArithmeticError, "fixed reference"):
            dofs.validate(invalid)

    def test_group_reference_and_cycle_errors_are_explicit(self):
        p = np.array([[1.0, 0.0], [2.0, 1.0]])
        metadata = records(p)
        metadata[0]["constraints"]["x"] = {"Grouped": [{"Edge": 99}, "Any"]}
        with self.assertRaisesRegex(ValueError, "Unknown native group reference"):
            MovementDOFs(["0", "1"], metadata)
        metadata = records(p, x_groups=[(0, 1)])
        with self.assertRaisesRegex(ArithmeticError, "grouped coordinates"):
            MovementDOFs(["0", "1"], metadata)

    def test_frozen_m4_final_geometry_can_move_group_x_and_every_free_y(self):
        fixture = json.loads(
            (
                Path(__file__).parent / "fixtures-joint-movement/m4-blocked-final.json"
            ).read_text()
        )
        keys, ids, p, _, segments, _, _ = model_arrays(fixture)
        dofs = MovementDOFs(keys, fixture["constraints"])
        outgoing = [ids[key] for key in fixture["outgoing_ids"]]
        candidates = []
        right = np.zeros_like(p)
        right[outgoing, 0] = 0.001
        step, _ = protect(p, segments, right, dofs=dofs, cap=0.002)
        np.testing.assert_allclose(step, right, atol=1e-11)
        candidates.append(step)
        # Every free Y has a permitted direction. In particular x:5 must go
        # down next to edge 16, while x:17 can independently move either way.
        for key in fixture["outgoing_ids"]:
            alternatives = []
            for sign in (-1.0, 1.0):
                proposal = np.zeros_like(p)
                proposal[ids[key], 1] = sign * 0.001
                step, _ = protect(p, segments, proposal, dofs=dofs, cap=0.002)
                alternatives.append((abs(step[ids[key], 1]), step))
            retained, step = max(alternatives, key=lambda item: item[0])
            self.assertGreater(retained, 0.00099, key)
            candidates.append(step)
        for step in candidates:
            dofs.validate(p + step)
            for fraction in (0.25, 0.5, 1.0):
                positions = {
                    key: point.tolist() for key, point in zip(keys, p + fraction * step)
                }
                geometry = validate_geometry(positions, fixture["routes"])
                self.assertTrue(geometry["valid"], geometry)

    def test_frozen_m4_blocked_proposal_rotates_without_freezing_other_legs(self):
        fixture = json.loads(
            (
                Path(__file__).parent / "fixtures-joint-movement/m4-blocked-final.json"
            ).read_text()
        )
        keys, ids, p, _, segments, _, _ = model_arrays(fixture)
        dofs = MovementDOFs(keys, fixture["constraints"])
        proposal = np.array(
            [fixture["blocked_proposal"].get(key, [0.0, 0.0]) for key in keys]
        )
        points, normals, bounds = separating_halfplanes(p, segments)
        towards = np.einsum("ij,ij->i", normals, proposal[points])
        active = (points == ids["x:5"]) & (towards > 0)
        self.assertLess(float(np.min(bounds[active] / towards[active])), 1e-4)
        evidence = {}
        step, _ = protect(
            p, segments, proposal, dofs=dofs, cap=0.0048, certificate=evidence
        )
        self.assertGreater(step[ids["x:17"], 1], 0.003)
        self.assertLess(step[ids["x:5"], 1], -0.001)
        self.assertGreater(step[ids["x:5"], 0], proposal[ids["x:5"], 0] + 1e-6)
        self.assertEqual(step[ids["x:5"], 0], step[ids["x:17"], 0])
        self.assertLess(evidence["separator_residual"], 1e-11)
        after = {key: point.tolist() for key, point in zip(keys, p + step)}
        self.assertTrue(validate_geometry(after, fixture["routes"])["valid"])

    def test_frozen_thin_regions_have_primal_and_optimality_certificates(self):
        for path in sorted(
            (Path(__file__).parent / "fixtures-joint-movement").glob(
                "*-thin-region*.npz"
            )
        ):
            with self.subTest(path=path.name), np.load(path) as problem:
                dofs = MovementDOFs(range(len(problem["proposal"])))
                # The fixture freezes the already-eliminated native axis map.
                dofs.index = problem["indices"]
                dofs.weights = problem["weights"]
                dofs.count = len(dofs.weights)
                step, evidence = project_halfplanes(
                    problem["proposal"],
                    dofs,
                    problem["points"],
                    problem["normals"],
                    problem["bounds"],
                    scale=float(problem["scale"]),
                )
                self.assertTrue(np.isfinite(step).all())
                self.assertLess(
                    evidence["primal_residual"], 2e-9 * float(problem["scale"])
                )
                self.assertLess(evidence["stationarity_residual"], 2e-8)
                self.assertLess(evidence["complementarity_residual"], 2e-8)
                self.assertGreater(evidence["components"], 1)

    def test_grouped_m4_multiplier_amplification_gets_certified_refinement(self):
        from scipy import sparse

        path = (
            Path(__file__).parent / "fixtures-joint-movement/m4-grouped-component.npz"
        )
        with np.load(path) as problem:
            values, evidence = _project_component(
                sparse.csc_matrix(problem["matrix"]),
                problem["bounds"],
                problem["quadratic"],
                problem["desired"],
            )
        self.assertTrue(np.isfinite(values).all())
        self.assertIn("polished", evidence["status"])
        self.assertLess(evidence["primal_residual"], 2e-9)
        self.assertLess(evidence["stationarity_residual"], 2e-8)
        self.assertLess(evidence["complementarity_residual"], 2e-8)

    def test_second_grouped_m4_face_is_polished_without_weakening_certificates(self):
        from scipy import sparse

        path = (
            Path(__file__).parent
            / "fixtures-joint-movement/m4-grouped-second-component.npz"
        )
        with np.load(path) as problem:
            values, evidence = _project_component(
                sparse.csc_matrix(problem["matrix"]),
                problem["bounds"],
                problem["quadratic"],
                problem["desired"],
            )
        self.assertTrue(np.isfinite(values).all())
        self.assertIn("polished", evidence["status"])
        self.assertLess(evidence["primal_residual"], 2e-9)
        self.assertLess(evidence["stationarity_residual"], 2e-8)
        self.assertLess(evidence["complementarity_residual"], 2e-8)

    def test_solver_failure_is_not_substituted_with_zero_movement(self):
        p = np.array([[1.0, 0.0], [1.0, 1.0]])
        dofs = MovementDOFs(["0", "1"], records(p, x_groups=[(0, 1)]))
        with patch("movement.clarabel.DefaultSolver") as optimizer:
            optimizer.return_value.solve.return_value = SimpleNamespace(
                x=[float("nan")],
                z=[float("nan")],
                status="explicit test failure",
                iterations=1,
            )
            with self.assertRaisesRegex(ArithmeticError, "Uncertified joint movement"):
                protect(
                    p,
                    np.array([[0, 1]]),
                    np.array([[10.0, 0.0], [10.0, 0.0]]),
                    cap=0.5,
                    dofs=dofs,
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
