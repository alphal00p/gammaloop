"""Real-solver regression coverage for the interactive stage backend."""

import copy
import hashlib
import json
import math
import unittest
from concurrent.futures import CancelledError
from unittest.mock import patch

from playground_runner import DEFAULTS, OUT, PlaygroundRunner
from solver import solve


class PlaygroundRunnerTests(unittest.TestCase):
    def setUp(self):
        self.runner = PlaygroundRunner()

    def test_configuration_is_complete_and_returned_data_is_independent(self):
        config = self.runner.config()
        case_ids = [case["id"] for case in config["cases"]]
        required = {
            "paper-b1",
            "qed-compton",
            "qcd-radiation",
            "higgs-sunset",
            "electron-rainbow",
            "crossed-higgs-double-box",
            "multi-leg-six-higgs-theta",
            "multi-leg-eight-gluon-ladder",
            "multi-leg-six-gluon-incoming",
        }
        available = {
            case["id"]
            for manifest in OUT.glob("cases*.json")
            for case in json.loads(manifest.read_text())
            if (OUT / case["id"] / "constrained-initial.json").exists()
        }
        self.assertTrue(required <= set(case_ids))
        self.assertSetEqual(set(case_ids), available)
        self.assertEqual(len(case_ids), len(set(case_ids)))
        self.assertEqual(len(config["parameters"]), 13)
        self.assertEqual(config["defaults"], DEFAULTS)
        for case in config["cases"]:
            self.assertEqual(
                [stage["id"] for stage in case["stages"]], ["seed", "impred"]
            )
        self.assertTrue(
            all(
                stage["available"]
                for case in config["cases"]
                for stage in case["stages"]
            )
        )
        for spec in config["parameters"]:
            self.assertLessEqual(spec["min"], spec["default"])
            self.assertLessEqual(spec["default"], spec["max"])
        json.dumps(config, allow_nan=False)
        config["cases"][0]["stages"][0]["positions"].clear()
        config["defaults"].clear()
        self.assertTrue(self.runner.config()["defaults"])
        self.assertTrue(self.runner.config()["cases"][0]["stages"][0]["positions"])

    def test_validation_rejects_unknown_nonfinite_out_of_range_and_noninteger_values(
        self,
    ):
        for parameters in [
            None,
            [],
            {"unknown": 1},
            {"impred.pull": True},
            {"impred.pull": "1"},
            {"impred.pull": math.nan},
            {"impred.pull": math.inf},
            {"impred.pull": 10**400},
            {"impred.edge_clearance": 0},
            {"impred.steps": 100.5},
            {"impred.external_max_points": 1.5},
            {"impred.external_max_points": 4},
            {"impred.external_length_ratio": 1.5},
            {"impred.pull_balance": -0.1},
            {"impred.pull_attachment": -0.1},
            {"impred.pull_attachment": math.inf},
            {"impred.split_length_ratio": 0.1},
            {"impred.contract_chord_ratio": 0},
            {"impred.split_length_ratio": 1.5, "impred.contract_chord_ratio": 1.5},
            {"impred.split_length_ratio": 1.5, "impred.contract_chord_ratio": 2.0},
            {"impred.external_direction_strength": 96},
            {"springs.route-rest-scale": 1},
            {"anneal.cool": 0.85},
        ]:
            with (
                self.subTest(parameters=str(parameters)[:100]),
                self.assertRaises(ValueError),
            ):
                self.runner.validate("higgs-sunset", "impred", parameters)
        for diagram, stage in [
            ("../higgs-sunset", "impred"),
            ("higgs-sunset", "invalid"),
            ("higgs-sunset", "springs"),
            ("higgs-sunset", "anneal"),
        ]:
            with self.assertRaises(ValueError):
                self.runner.validate(diagram, stage, {})
        normalized = self.runner.validate(
            "higgs-sunset", "impred", {"impred.steps": 100.0}
        )
        self.assertIs(type(normalized["impred.steps"]), int)
        self.assertEqual(set(normalized), set(DEFAULTS))

    def test_default_live_stages_exactly_reproduce_saved_geometry_without_file_changes(
        self,
    ):
        folder = OUT / "higgs-sunset"
        names = [
            "constrained-initial.json",
            "constrained-initial.projected.graph.json",
            "constrained-impred.json",
            "constrained-native-force-only.json",
            "constrained-native-force-cleanup.json",
            "native-prepared.json",
            "native-prepared.bin",
        ]
        before = {
            name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
            for name in names
        }
        saved = next(
            case
            for case in self.runner.config()["cases"]
            if case["id"] == "higgs-sunset"
        )
        result = self.runner.run("higgs-sunset", "impred", {})
        for actual, expected in zip(result["case"]["stages"], saved["stages"]):
            self.assertTrue(actual["available"])
            for key in ("positions", "routes", "anchor_ids", "counts", "constraints"):
                self.assertEqual(actual[key], expected[key], (actual["id"], key))
        self.assertEqual(result["parameters"], DEFAULTS)
        self.assertEqual(
            before,
            {
                name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
                for name in names
            },
        )

    def test_force_knobs_change_actual_solver_geometry(self):
        params = {"impred.steps": 100}
        reference = self.runner.run("higgs-sunset", "impred", params)["case"]["stages"][
            1
        ]
        for key, value in [
            ("impred.repulsion", 0.4),
            ("impred.attraction", 2.0),
            ("impred.parallel_attraction_balance", 0.0),
        ]:
            changed = self.runner.run("higgs-sunset", "impred", {**params, key: value})[
                "case"
            ]["stages"][1]
            self.assertNotEqual(changed["positions"], reference["positions"])
            self.assertEqual(changed["routes"], reference["routes"])
            self.assertTrue(changed["constraints"]["valid"])
            self.assertEqual(changed["counts"], reference["counts"])

    def test_cache_reuses_matching_parameters_and_returns_independent_payloads(self):
        params = {"impred.steps": 100}
        with patch("playground_runner.solve", wraps=solve) as actual_solver:
            baseline = self.runner.run("qed-compton", "impred", params)
            changed = self.runner.run(
                "qed-compton", "impred", {**params, "impred.target": 3.0}
            )
            repeat = self.runner.run("qed-compton", "impred", params)
            self.assertEqual(actual_solver.call_count, 2)
        self.assertEqual(baseline["case"], repeat["case"])
        self.assertNotEqual(
            baseline["case"]["stages"][1]["positions"],
            changed["case"]["stages"][1]["positions"],
        )
        repeat["case"]["stages"][1]["positions"].clear()
        self.assertTrue(
            self.runner.run("qed-compton", "impred", params)["case"]["stages"][1][
                "positions"
            ]
        )

    def test_seed_selection_does_not_run_relaxation(self):
        with patch("playground_runner.solve", wraps=solve) as actual_solver:
            result = self.runner.run("qed-compton", "seed", {})
            actual_solver.assert_not_called()
        seed, impred = result["case"]["stages"]
        self.assertTrue(seed["available"])
        self.assertFalse(impred["available"])
        self.assertFalse(self.runner._cache)

    def test_all_available_cases_keep_constraints_bounded_routes_and_physical_ids(self):
        for case in self.runner.config()["cases"]:
            with self.subTest(case=case["id"]):
                result = self.runner.run(case["id"], "impred", {"impred.steps": 100})
                seed = result["case"]["stages"][0]
                for stage in result["case"]["stages"]:
                    self.assertTrue(stage["constraints"]["valid"])
                    self.assertLess(stage["constraints"]["max_residual"], 1e-9)
                    for edge, route in stage["routes"].items():
                        if edge in stage["external_ids"]:
                            self.assertLessEqual(len(route) - 2, 2)
                            self.assertEqual(route[0], seed["routes"][edge][0])
                            self.assertEqual(route[-1], seed["routes"][edge][-1])
                        else:
                            self.assertEqual(route, seed["routes"][edge])
                    self.assertEqual(stage["node_ids"], seed["node_ids"])
                    self.assertEqual(stage["edge_endpoints"], seed["edge_endpoints"])
                    self.assertTrue(
                        all(
                            math.isfinite(v)
                            for p in stage["positions"].values()
                            for v in p
                        )
                    )
                self.assertEqual(result["case"]["stages"][1]["counts"], seed["counts"])

    def test_balancing_and_external_route_controls_reach_the_solver(self):
        parameters = {
            "impred.steps": 100,
            "impred.pull": 3.6,
            "impred.pull_balance": 2.0,
            "impred.pull_attachment": 1.0,
            "impred.parallel_attraction_balance": 0.25,
            "impred.external_max_points": 1,
            "impred.split_length_ratio": 1.5,
            "impred.contract_chord_ratio": 0.75,
        }
        with patch("playground_runner.solve", wraps=solve) as actual_solver:
            result = self.runner.run(
                "multi-leg-eight-gluon-ladder", "impred", parameters
            )
        knobs = actual_solver.call_args.kwargs
        self.assertEqual(knobs["pull_balance"], 2.0)
        self.assertEqual(knobs["parallel_attraction_balance"], 0.25)
        self.assertNotIn("pull_attachment", knobs)
        self.assertNotIn("external_direction_strength", knobs)
        self.assertEqual(knobs["external_max_points"], 1)
        self.assertEqual(knobs["split_length_ratio"], 1.5)
        self.assertEqual(knobs["contract_chord_ratio"], 0.75)
        self.assertNotIn("external_length_ratio", knobs)
        self.assertTrue(any(scale < 1 for scale in knobs["pull_scales"].values()))
        stage = result["case"]["stages"][1]
        self.assertTrue(stage["constraints"]["valid"])
        self.assertTrue(
            all(len(stage["routes"][edge]) <= 3 for edge in stage["external_ids"])
        )

    def test_cancelled_work_is_not_published_and_solver_checks_during_execution(self):
        with self.assertRaises(CancelledError):
            self.runner.run("higgs-sunset", "impred", {}, cancelled=lambda: True)
        self.assertFalse(self.runner._cache)
        seed = json.loads((OUT / "higgs-sunset/constrained-initial.json").read_text())
        before = copy.deepcopy(seed)
        snapshots = []
        with self.assertRaises(CancelledError):
            solve(
                seed,
                steps=100,
                flexible=False,
                cancelled=lambda: bool(snapshots),
                progress=snapshots.append,
            )
        self.assertEqual(len(snapshots), 1)
        self.assertEqual(snapshots[0]["iteration"], 1)
        self.assertEqual(seed, before)

    def test_attachment_gain_uses_native_topology_and_preserves_source_graph(self):
        for diagram, distributed in [
            ("multi-leg-eight-gluon-ladder", True),
            ("xs-extra-higgs-scattering", False),
        ]:
            with self.subTest(diagram=diagram):
                graph = self.runner._sources[diagram]["graph"]
                scale_sets = []
                with patch("playground_runner.solve", wraps=solve) as actual_solver:
                    for gain in (1.0, 16.0):
                        self.runner.run(
                            diagram,
                            "impred",
                            {"impred.steps": 100, "impred.pull_attachment": gain},
                        )
                        knobs = actual_solver.call_args.kwargs
                        self.assertNotIn("pull_attachment", knobs)
                        scale_sets.append(knobs["pull_scales"])
                    self.assertEqual(actual_solver.call_count, 2)
                before, after = scale_sets
                self.assertEqual(before.keys(), after.keys())
                self.assertTrue(all(after[e] >= before[e] for e in before))
                self.assertEqual(before != after, distributed)
                self.assertEqual(self.runner._sources[diagram]["graph"], graph)


if __name__ == "__main__":
    unittest.main(verbosity=2)
