import copy
import json
import math
import os
import shlex
import signal
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_performance as producer


class ProducerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.repository = Path(__file__).resolve().parents[2]
        self.suite = producer.tomllib.loads(
            (self.repository / "benchmarks/performance/suite.toml").read_text()
        )
        self.fixture = producer.Fixture.from_suite(self.suite, "gl15", self.repository)
        self.value = {
            "process_id": 0,
            "integrand_name": "1L",
            "point": self.fixture.corpus["points"][0],
            "discrete_dim": [0, 0, 0],
            "target_seconds": self.fixture.config["bench_duration_seconds"],
            "warmup_samples": 10,
            "warmup_seconds": 0.01,
            "minimal_integrand": False,
            "momentum_space": False,
            "n_batches": 10,
            "total_samples": 50,
            "displayed_result": {"re": 1.0, "im": 0.5},
            "batches": [{"total": 0.001}] * 10,
            "summary": [{"category": "Total", "mean_seconds_per_sample": 0.001}],
        }
        self.generation_summary = {
            "peak_ram_bytes": 128 * 1024**2,
            "reports": [
                {
                    "graph_name": "GL15",
                    "stats": {
                        "evaluator_count": 4,
                        "total_time": {"secs": 1, "nanos": 0},
                        "evaluator_spenso_time": {"secs": 0, "nanos": 1},
                        "evaluator_symbolica_time": {"secs": 0, "nanos": 1},
                        "evaluator_compile_time": {"secs": 0, "nanos": 0},
                    },
                }
            ],
        }
        self.all_methods = producer.Fixture(
            self.fixture.card.replace(
                'orientations = "monte_carlo"', 'orientations = "summed"'
            ),
            self.fixture.corpus | {"discrete_dim": [0, 0]},
            self.fixture.config
            | {
                "evaluator_methods": list(producer.EVALUATOR_METHODS),
                "reference_method": "SingleParametric",
                "required_metrics": [
                    "generation_seconds",
                    "warm_runtime_seconds_per_sample",
                    "evaluator_build_seconds",
                    "generation_peak_ram_bytes",
                ],
            },
            self.repository,
        )
        self.construction_events = [
            {
                "stage": "evaluator_stack_new_start",
                "orientation_count": 6,
                "iterative_orientation_optimization": True,
                "summed": True,
                "summed_function_map": True,
            },
            {
                "stage": "evaluator_stack_new_iterative_done",
                "orientation_count": 6,
                "elapsed_ms": 2000,
            },
            {
                "stage": "evaluator_stack_new_summed_function_map_done",
                "orientation_count": 6,
                "elapsed_ms": 7000,
            },
            {
                "stage": "evaluator_stack_new_summed_done",
                "orientation_count": 6,
                "elapsed_ms": 5000,
            },
            {
                "stage": "evaluator_stack_new_single_parametric_done",
                "orientation_count": 6,
                "elapsed_ms": 3000,
            },
            {"stage": "evaluator_stack_new_done", "orientation_count": 6},
        ]

    def read_result(self, value):
        path = self.directory / "bench.json"
        path.write_text(json.dumps(value))
        return producer.benchmark_result(path, self.fixture, 0)

    def test_production_sample_timing_is_not_cli_or_batch_wall_time(self):
        self.assertEqual(self.read_result(self.value), (0.001, [1.0, 0.5], 50))

    def test_rejects_incomplete_or_changed_benchmark_work(self):
        changes = [
            ("process_id", 1),
            ("integrand_name", "other"),
            ("point", [0.1, 0.2, 0.3]),
            ("discrete_dim", [1, 0, 0]),
            ("target_seconds", 0.1),
            ("warmup_samples", 0),
            ("minimal_integrand", True),
            ("momentum_space", True),
            ("n_batches", 1),
            ("total_samples", 0),
            ("batches", []),
            ("summary", [{"category": "Total", "mean_seconds_per_sample": 0.01}]),
            ("displayed_result", {"re": 0.0, "im": 0.0}),
            ("displayed_result", {"re": math.inf, "im": 0.0}),
        ]
        for key, replacement in changes:
            with self.subTest(key=key, replacement=replacement):
                value = copy.deepcopy(self.value)
                value[key] = replacement
                with self.assertRaises(ValueError):
                    self.read_result(value)

    def test_expired_shared_budget_never_starts_another_cli(self):
        with patch.object(producer.subprocess, "Popen") as start:
            with self.assertRaises(ValueError):
                producer.run_cli(
                    ["unused"],
                    self.directory,
                    self.directory / "log",
                    time.monotonic() - 1,
                )
            start.assert_not_called()

    def test_cli_completion_keeps_output_and_reaps_the_process(self):
        start = producer.subprocess.Popen
        processes = []

        def record(*args, **kwargs):
            process = start(*args, **kwargs)
            processes.append(process)
            return process

        log = self.directory / "cli.log"
        with patch.object(producer.subprocess, "Popen", side_effect=record):
            producer.run_cli(
                [sys.executable, "-c", "print('finished')"],
                self.directory,
                log,
                time.monotonic() + 10,
            )
        self.assertEqual(log.read_text().strip(), "finished")
        self.assertEqual(processes[0].returncode, 0)
        with self.assertRaises(ProcessLookupError):
            os.kill(processes[0].pid, 0)

    def test_cli_nonzero_exit_retains_private_log_error(self):
        log = self.directory / "failed.log"
        with self.assertRaisesRegex(ValueError, "CLI exited with status 7"):
            producer.run_cli(
                [sys.executable, "-c", "raise SystemExit(7)"],
                self.directory,
                log,
                time.monotonic() + 10,
            )

    def test_cli_deadline_escalates_and_reaps_a_term_resistant_process(self):
        start = producer.subprocess.Popen
        processes = []

        def record(*args, **kwargs):
            process = start(*args, **kwargs)
            processes.append(process)
            return process

        log = self.directory / "timeout.log"
        command = (
            "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
            "print('ready', flush=True); time.sleep(60)"
        )
        with (
            patch.object(producer.subprocess, "Popen", side_effect=record),
            patch.object(producer.os, "killpg", wraps=os.killpg) as terminate,
            self.assertRaisesRegex(ValueError, "240-second budget"),
        ):
            producer.run_cli(
                [sys.executable, "-c", command],
                self.directory,
                log,
                time.monotonic() + 1,
            )
        self.assertEqual(log.read_text().strip(), "ready")
        self.assertEqual(processes[0].returncode, -signal.SIGKILL)
        self.assertEqual(
            [call.args[1] for call in terminate.call_args_list][:2],
            [signal.SIGTERM, signal.SIGKILL],
        )
        with self.assertRaises(ProcessLookupError):
            os.kill(processes[0].pid, 0)

    def test_native_stage_durations_allow_zero_and_reject_malformed_before_bench(self):
        for index, duration in enumerate(
            [
                {"secs": 0, "nanos": 0},
                {"secs": -1, "nanos": 0},
                {"secs": True, "nanos": 0},
                {"secs": 2**64, "nanos": 0},
                {"secs": 0, "nanos": 10**9},
                {"secs": 0, "nanos": 0.5},
            ]
        ):
            with self.subTest(duration=duration):
                directory = self.directory / str(index)

                def generate(*_, directory=directory, duration=duration, **_kwargs):
                    summary = directory / "state" / "generation_summary.json"
                    summary.parent.mkdir(exist_ok=True)
                    contents = copy.deepcopy(self.generation_summary)
                    contents["reports"][0]["stats"]["evaluator_compile_time"] = duration
                    producer.write_json(
                        summary,
                        contents,
                    )

                with (
                    patch.object(producer, "run_cli", side_effect=generate) as run,
                    patch.object(
                        producer,
                        "benchmark_result",
                        return_value=(0.001, [1.0, 0.5], 50),
                    ),
                ):
                    arguments = (
                        sys.executable,
                        self.repository,
                        directory,
                        self.fixture,
                        time.monotonic() + 10,
                    )
                    if index == 0:
                        _, _, metadata = producer.measure(*arguments)
                        self.assertEqual(
                            metadata["generation_evaluator_compile_seconds"], 0
                        )
                        self.assertEqual(run.call_count, 2)
                    else:
                        with self.assertRaises(ValueError):
                            producer.measure(*arguments)
                        self.assertEqual(run.call_count, 1)

    def test_fixed_corpus_cannot_silently_enable_point_cache(self):
        with self.assertRaises(ValueError):
            producer.Fixture(
                self.fixture.card.replace(
                    "enable_cache = false", "enable_cache = true"
                ),
                self.fixture.corpus,
                self.fixture.config,
                self.repository,
            )

    def test_one_fresh_generation_serves_all_methods_with_distinct_runtime_metrics(
        self,
    ):
        directory = self.directory / "all-methods"
        calls = []

        def run(command, repository, log, deadline, **kwargs):
            calls.append(command)
            if "--clean-state" in command:
                producer.write_json(
                    directory / "state" / "generation_summary.json",
                    self.generation_summary,
                )
                logfile = directory / "state" / "logs" / "gammalog-generation.jsonl"
                logfile.parent.mkdir()
                logfile.write_text(
                    "".join(
                        json.dumps(event) + "\n" for event in self.construction_events
                    )
                )
                self.assertIn("#summary", kwargs["logfile_filter"])
                return
            self.assertIn("--read-only-state", command)
            self.assertIn("-n", command)
            commands = producer.tomllib.loads(Path(command[-1]).read_text())["commands"]
            self.assertEqual(commands[-1], "quit --no-save-state")
            method = None
            seen = []
            for text in commands[:-1]:
                tokens = shlex.split(text)
                if tokens[0] == "set":
                    method = tokens[-1].split("=")[1]
                    seen.append(method)
                    continue
                self.assertEqual(tokens[0], "bench")
                index = sum(1 for path in directory.glob(f"bench-{method}-*.json"))
                value = copy.deepcopy(self.value)
                value["point"] = self.all_methods.corpus["points"][index]
                value["discrete_dim"] = self.all_methods.corpus["discrete_dim"]
                timing = (producer.EVALUATOR_METHODS.index(method) + 1) * 0.001
                value["batches"] = [{"total": timing}] * 10
                value["summary"] = [
                    {"category": "Total", "mean_seconds_per_sample": timing}
                ]
                producer.write_json(
                    Path(tokens[tokens.index("--json-output") + 1]), value
                )
            self.assertEqual(seen, list(producer.EVALUATOR_METHODS))

        with patch.object(producer, "run_cli", side_effect=run):
            metrics, outputs, metadata = producer.measure(
                sys.executable,
                self.repository,
                directory,
                self.all_methods,
                time.monotonic() + 10,
            )
        self.assertEqual(len(calls), 2)
        self.assertEqual(sum("--clean-state" in command for command in calls), 1)
        self.assertEqual(metrics["generation_peak_ram_bytes"], 128 * 1024**2)
        for index, method in enumerate(producer.EVALUATOR_METHODS):
            self.assertAlmostEqual(
                metrics[f"warm_runtime_seconds_per_sample.{method}"],
                (index + 1) * 0.001,
            )
            self.assertEqual(outputs[method], [1.0, 0.5] * 4)
            self.assertEqual(metadata["evaluator_build_counts"][method], 1)
        self.assertEqual(metrics["evaluator_build_seconds.Iterative"], 2.0)
        self.assertEqual(metrics["evaluator_build_seconds.SummedFunctionMap"], 7.0)

    def test_construction_timings_sum_complete_distinct_stacks(self):
        state = self.directory / "state"
        logs = state / "logs"
        logs.mkdir(parents=True)
        (logs / "gammalog-generation.jsonl").write_text(
            "".join(json.dumps(event) + "\n" for event in self.construction_events * 2)
        )
        metrics, metadata = self.all_methods.construction_metrics(
            state, time.monotonic() + 10
        )
        self.assertEqual(
            metrics,
            {
                "evaluator_build_seconds.Iterative": 4.0,
                "evaluator_build_seconds.SingleParametric": 6.0,
                "evaluator_build_seconds.Summed": 10.0,
                "evaluator_build_seconds.SummedFunctionMap": 14.0,
            },
        )
        self.assertEqual(metadata["construction_orientation_counts"], [6, 6])
        self.assertEqual(
            metadata["evaluator_build_counts"],
            dict.fromkeys(producer.EVALUATOR_METHODS, 2),
        )

    def test_missing_duplicate_or_incomplete_construction_is_rejected(self):
        changes = (
            lambda events: events.clear(),
            lambda events: events.pop(2),
            lambda events: events.insert(2, events[1]),
            lambda events: events.pop(),
            lambda events: events.pop(0),
            lambda events: events.insert(2, events[0]),
            lambda events: events[1].update(elapsed_ms=math.inf),
            lambda events: events[1].update(orientation_count=True),
            lambda events: events.insert(0, []),
        )
        for index, change in enumerate(changes):
            with self.subTest(change=index):
                events = copy.deepcopy(self.construction_events)
                change(events)
                state = self.directory / str(index)
                logs = state / "logs"
                logs.mkdir(parents=True)
                (logs / "gammalog-generation.jsonl").write_text(
                    "".join(json.dumps(event) + "\n" for event in events)
                )
                with self.assertRaises((TypeError, ValueError)):
                    self.all_methods.construction_metrics(state, time.monotonic() + 10)

    def test_invalid_peak_ram_stops_before_runtime(self):
        for index, peak in enumerate((0, -1, True, 0.5, 2**64)):
            with self.subTest(peak=peak):
                directory = self.directory / str(index)
                summary = self.generation_summary | {"peak_ram_bytes": peak}

                def generate(*_args, directory=directory, summary=summary, **_kwargs):
                    producer.write_json(
                        directory / "state" / "generation_summary.json", summary
                    )

                with (
                    patch.object(producer, "run_cli", side_effect=generate) as run,
                    self.assertRaisesRegex(ValueError, "peak RAM"),
                ):
                    producer.measure(
                        sys.executable,
                        self.repository,
                        directory,
                        self.fixture,
                        time.monotonic() + 10,
                    )
                self.assertEqual(run.call_count, 1)

    def test_confirmation_cannot_reset_the_shared_case_deadline(self):
        args = SimpleNamespace(
            baseline_commit="a" * 40,
            candidate_commit="b" * 40,
            baseline_binary=Path("before"),
            candidate_binary=Path("after"),
            input_root=self.repository,
            work_dir=self.directory / "private",
            output_dir=self.directory / "public",
        )
        now = [100.0]
        calls = []

        def measure(binary, repository, directory, fixture, deadline):
            calls.append((binary, deadline))
            if len(calls) == 10:
                now[0] = deadline + 0.01
            slower = binary == args.candidate_binary
            return (
                {
                    "generation_seconds": 13.0 if slower else 10.0,
                    "warm_runtime_seconds_per_sample.SingleParametric": 1e-5,
                },
                {"SingleParametric": [1.0, 0.5] * 4},
                {},
            )

        with (
            patch.object(producer.time, "monotonic", side_effect=lambda: now[0]),
            patch.object(producer, "measure", side_effect=measure),
        ):
            self.assertEqual(producer.run_case(args, self.suite, "gl15"), 2)
        self.assertEqual(len(calls), 10)
        self.assertEqual({deadline for _, deadline in calls}, {340.0})
        self.assertNotIn(
            "confirmation",
            json.loads((args.output_dir / "measurements.json").read_text()),
        )
        report = json.loads((args.output_dir / "report.json").read_text())
        self.assertEqual(report["invalid_input"], 1)

    def test_summed_scalar_corpus_has_nine_coordinates_without_discrete_selectors(self):
        scalar = producer.Fixture.from_suite(
            self.suite, "gl04-q1-squared-orientation-3d", self.repository
        )
        card = scalar.benchmark_card(self.directory)
        commands = producer.tomllib.loads(card.read_text())["commands"]
        for command in commands:
            tokens = shlex.split(command)
            if tokens[0] != "bench":
                continue
            self.assertNotIn("-d", tokens)
            self.assertEqual(tokens.index("--duration") - tokens.index("-x") - 1, 9)

    def test_graph_assets_are_part_of_the_shared_fixture_identity(self):
        corpus = self.fixture.corpus | {"inputs": ["graph.dot"]}
        asset = self.directory / "graph.dot"
        asset.write_text("first input")
        first = producer.Fixture(
            self.fixture.card, corpus, self.fixture.config, self.directory
        )
        asset.write_text("changed numerator")
        second = producer.Fixture(
            self.fixture.card, corpus, self.fixture.config, self.directory
        )
        self.assertNotEqual(
            first.hashes["fixture_sha256"], second.hashes["fixture_sha256"]
        )

    def test_scalar_corpus_rejects_malformed_dimensions_and_external_assets(self):
        for field, value in [
            ("continuous_dimension", 0),
            ("continuous_dimension", 4),
            ("continuous_dimension", 9),
            ("discrete_dim", [-1]),
            ("inputs", ["../graph.dot"]),
            ("inputs", ["/graph.dot"]),
        ]:
            with self.subTest(field=field, value=value):
                corpus = copy.deepcopy(self.fixture.corpus)
                corpus[field] = value
                with self.assertRaises(ValueError):
                    producer.Fixture(
                        self.fixture.card,
                        corpus,
                        self.fixture.config,
                        self.fixture.input_root,
                    )

    def test_registered_corpus_covers_three_nonzero_scalar_probes_in_all_routes(self):
        cases = producer.Fixture.cases(self.suite)
        expected = {"gl15"} | {
            f"{probe}-{route}"
            for probe in ("gl04-q1-squared", "gl16-q7-squared", "gl02-energy-quartic")
            for route in ("orientation-3d", "explicit-3d", "projected-4d")
        }
        self.assertEqual(
            cases.keys(),
            expected
            | {"gl16-evaluator-aliases", "gl1-evaluator-guards", "gl1-thermal"},
        )
        for name in expected:
            with self.subTest(case=name):
                fixture = producer.Fixture.from_suite(self.suite, name, self.repository)
                self.assertEqual(len(fixture.corpus["points"]), 4)
                if name != "gl15":
                    self.assertEqual(fixture.corpus["continuous_dimension"], 9)
                    self.assertEqual(fixture.corpus["discrete_dim"], [])
                    self.assertEqual(
                        fixture.config["correctness_tolerance"]["absolute"], 0
                    )

    def test_registered_evaluator_probes_build_all_methods_from_hashed_graphs(self):
        for name, dimension, absolute, graph in (
            ("gl16-evaluator-aliases", 9, 0.0, "scalar-3l-gl16-q7-squared.dot"),
            ("gl1-evaluator-guards", 6, 1e-14, "gluon-self-energy-gl1.dot"),
        ):
            with self.subTest(case=name):
                fixture = producer.Fixture.from_suite(self.suite, name, self.repository)
                self.assertEqual(len(fixture.corpus["points"]), 2)
                self.assertEqual(fixture.corpus["continuous_dimension"], dimension)
                self.assertEqual(fixture.corpus["discrete_dim"], [0, 0])
                self.assertEqual(
                    fixture.config["correctness_tolerance"]["absolute"], absolute
                )
                self.assertEqual(set(fixture.methods), set(producer.EVALUATOR_METHODS))
                self.assertEqual(fixture.config["reference_method"], "SingleParametric")
                self.assertEqual(
                    set(fixture.metrics),
                    {
                        "generation_seconds",
                        "generation_peak_ram_bytes",
                        *{
                            f"{metric}.{method}"
                            for metric in producer.MODE_METRICS
                            for method in producer.EVALUATOR_METHODS
                        },
                    },
                )
                card = producer.tomllib.loads(fixture.card)
                runtime = card["default_runtime_settings"]
                self.assertEqual(runtime["sampling"]["orientations"], "summed")
                self.assertEqual(
                    runtime["general"]["evaluator_method"], "SingleParametric"
                )
                evaluator = (
                    card["cli_settings"]["global"]
                    .get("generation", {})
                    .get("evaluator", {})
                    .copy()
                )
                # Cards can set generation options in TOML or via their ordered CLI commands.
                prefix = "global.generation.evaluator."
                for command in card["commands"]:
                    tokens = shlex.split(command)
                    if tokens[:3] == ["set", "global", "kv"]:
                        for setting in tokens[3:]:
                            key, _, value = setting.partition("=")
                            if key.startswith(prefix):
                                evaluator[key.removeprefix(prefix)] = (
                                    producer.tomllib.loads("value=" + value)["value"]
                                )
                for flag in (
                    "iterative_orientation_optimization",
                    "summed",
                    "summed_function_map",
                ):
                    self.assertIs(evaluator[flag], True)
                self.assertIs(evaluator["compile"], False)
                self.assertIs(evaluator["store_atom"], False)
                graph_input = "benchmarks/performance/" + graph
                self.assertIn(graph_input, fixture.corpus["inputs"])
                self.assertTrue(
                    any(
                        shlex.split(command)[:3] == ["import", "graphs", graph_input]
                        for command in card["commands"]
                    )
                )
                graph_path = self.repository / graph_input
                read_bytes = Path.read_bytes

                def changed_graph(path, graph_path=graph_path, read_bytes=read_bytes):
                    contents = read_bytes(path)
                    return contents + b"\n" if path == graph_path else contents

                with patch.object(Path, "read_bytes", changed_graph):
                    changed = producer.Fixture.from_suite(
                        self.suite, name, self.repository
                    )
                self.assertNotEqual(
                    fixture.hashes["fixture_sha256"], changed.hashes["fixture_sha256"]
                )

    def test_thermal_probe_keeps_nested_qcd_uv_forests(self):
        fixture = producer.Fixture.from_suite(
            self.suite, "gl1-thermal", self.repository
        )
        card = producer.tomllib.loads(fixture.card)
        generation = card["cli_settings"]["global"]["generation"]
        # Vacuum subtraction would drop the two-loop forest nesting both
        # one-loop counterterms of this QCD graph.
        self.assertEqual(
            generation["medium"],
            {"mode": "thermodynamic_equilibrium", "vacuum_subtraction": False},
        )
        self.assertIs(generation["uv"]["subtract_uv"], True)
        self.assertGreater(
            card["default_runtime_settings"]["general"]["inverse_temperature"], 0
        )
        self.assertIn(
            "benchmarks/performance/gluon-self-energy-gl1.dot", fixture.corpus["inputs"]
        )
        self.assertEqual(fixture.methods, ["SingleParametric", "SummedFunctionMap"])
        self.assertEqual(fixture.config["reference_method"], "SingleParametric")

    def test_suite_finishes_all_cases_and_propagates_any_failure(self):
        names = list(producer.Fixture.cases(self.suite))
        for failing_index in (0, len(names) // 2, len(names) - 1):
            for failure in (1, 2):
                with self.subTest(case=failing_index, failure=failure):
                    results = [0] * len(names)
                    results[failing_index] = failure
                    directory = self.directory / f"{failing_index}-{failure}"
                    arguments = [
                        "run_performance.py",
                        "--input-root",
                        str(self.repository),
                        "--baseline-binary",
                        sys.executable,
                        "--candidate-binary",
                        sys.executable,
                        "--baseline-commit",
                        "a" * 40,
                        "--candidate-commit",
                        "b" * 40,
                        "--work-dir",
                        str(directory / "private"),
                        "--output-dir",
                        str(directory / "public"),
                    ]
                    with (
                        patch.object(sys, "argv", arguments),
                        patch.object(producer, "run_case", side_effect=results) as run,
                    ):
                        self.assertEqual(producer.main(), failure)
                    self.assertEqual(
                        [call.args[2] for call in run.call_args_list], names
                    )
                    for call in run.call_args_list:
                        args, _, name = call.args
                        self.assertEqual(args.work_dir, directory / "private" / name)
                        self.assertEqual(args.output_dir, directory / "public" / name)


if __name__ == "__main__":
    unittest.main()
