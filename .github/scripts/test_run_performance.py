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
        ):
            with self.assertRaisesRegex(ValueError, "240-second budget"):
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

                def generate(*_, directory=directory, duration=duration):
                    summary = directory / "state" / "generation_summary.json"
                    summary.parent.mkdir(exist_ok=True)
                    producer.write_json(
                        summary,
                        {
                            "peak_ram_bytes": 1,
                            "reports": [
                                {
                                    "graph_name": "GL15",
                                    "stats": {
                                        "evaluator_count": 1,
                                        "total_time": {"secs": 1, "nanos": 0},
                                        "evaluator_spenso_time": {
                                            "secs": 0,
                                            "nanos": 1,
                                        },
                                        "evaluator_symbolica_time": {
                                            "secs": 0,
                                            "nanos": 1,
                                        },
                                        "evaluator_compile_time": duration,
                                    },
                                }
                            ],
                        },
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

    def test_summed_scalar_corpus_has_nine_coordinates_without_discrete_selectors(self):
        scalar = producer.Fixture.from_suite(
            self.suite, "gl04-q1-squared-orientation-3d", self.repository
        )
        card = scalar.benchmark_card(self.directory)
        commands = producer.tomllib.loads(card.read_text())["commands"]
        for command in commands[:-1]:
            tokens = shlex.split(command)
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
        self.assertEqual(cases.keys(), expected)
        for name in cases:
            with self.subTest(case=name):
                fixture = producer.Fixture.from_suite(self.suite, name, self.repository)
                self.assertEqual(len(fixture.corpus["points"]), 4)
                if name != "gl15":
                    self.assertEqual(fixture.corpus["continuous_dimension"], 9)
                    self.assertEqual(fixture.corpus["discrete_dim"], [])
                    self.assertEqual(
                        fixture.config["correctness_tolerance"]["absolute"], 0
                    )

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
