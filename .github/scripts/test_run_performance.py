import copy
import json
import math
import os
import shlex
import shutil
import signal
import sys
import tempfile
import time
import unittest
from argparse import Namespace
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
        self.hyperfine_read = patch(
            "gammabench.hyperfine.read_measurement",
            return_value={"process_wall_ns": 10_000_000_000, "returncode": 0},
        )
        self.hyperfine_read.start()
        self.addCleanup(self.hyperfine_read.stop)
        self.hyperfine_resolve = patch(
            "gammabench.hyperfine.resolve",
            return_value={
                "path": sys.executable,
                "sha256": producer.hashlib.sha256(
                    Path(sys.executable).read_bytes()
                ).hexdigest(),
                "version": "2.0.0",
            },
        )
        self.hyperfine_resolve.start()
        self.addCleanup(self.hyperfine_resolve.stop)
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
            self.assertRaisesRegex(ValueError, "declared budget"),
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
            if "--export-json" in command:
                command = shlex.split(command[-1])
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

    def test_prepared_case_hands_pins_and_original_budget_to_gammabench(self):
        args = SimpleNamespace(
            baseline_commit="a" * 40,
            candidate_commit="b" * 40,
            baseline_binary=Path(sys.executable),
            candidate_binary=Path(sys.executable),
            input_root=self.repository,
            work_dir=self.directory / "private",
            output_dir=self.directory / "public",
        )
        with patch("gammabench.gate.run", return_value=2) as gate:
            self.assertEqual(producer.run_case(args, self.suite, "gl15"), 2)
        prepared = json.loads(gate.call_args.args[0].read_text())
        policy = json.loads(gate.call_args.kwargs["policy_path"].read_text())
        self.assertEqual(policy["budget_seconds"], 240)
        self.assertEqual(policy["pairs"], 5)
        self.assertEqual(policy["thresholds"]["generation_seconds"], [0.20, 0.5])
        self.assertEqual(prepared["gate"]["baseline"]["commit"], "a" * 40)
        self.assertEqual(prepared["gate"]["candidate"]["commit"], "b" * 40)
        self.assertEqual(prepared["tasks"][0]["inputs"], prepared["tasks"][1]["inputs"])
        self.assertEqual(
            prepared["tasks"][0]["process_measurements"]["generation"]["metric"],
            "generation_seconds",
        )

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
            | {
                "gl16-evaluator-aliases",
                "gl1-evaluator-guards",
                "gl1-thermal",
                "paper-massless-soft-cff4d-arb",
                "gl04-q1-squared-ltd-vs-cff-4d",
            },
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

    @unittest.skipUnless(shutil.which("hyperfine"), "requires pinned Hyperfine")
    def test_all_thirteen_cases_run_through_real_hyperfine_and_seal(self):
        from gammabench.snapshot import load, verify

        self.hyperfine_read.stop()
        self.hyperfine_resolve.stop()
        binary = self.directory / "synthetic-cli"
        binary.write_text(
            "#!"
            + sys.executable
            + "\n"
            + r"""
import json, shlex, sys, tomllib
from pathlib import Path
args = sys.argv[1:]
fixture = json.loads((Path.cwd().parent / 'fixture.json').read_text())
state = Path(args[args.index('-s') + 1])
print('PRIVATE-SYNTHETIC-CLI-OUTPUT')
if '--clean-state' in args:
    state.mkdir()
    stats = {'evaluator_count': 1, **{name: {'secs': 0, 'nanos': 1} for name in
             ('total_time', 'evaluator_spenso_time', 'evaluator_symbolica_time', 'evaluator_compile_time')}}
    (state / 'generation_summary.json').write_text(json.dumps({'peak_ram_bytes': 128 * 1024**2,
        'reports': [{'graph_name': fixture['corpus']['graph_name'], 'stats': stats}]}))
    if 'evaluator_build_seconds' in fixture['config']['required_metrics']:
        logs = state / 'logs'
        logs.mkdir()
        events = [{'stage': 'evaluator_stack_new_start', 'orientation_count': 3,
                   'iterative_orientation_optimization': True, 'summed': True, 'summed_function_map': True}]
        events += [{'stage': 'evaluator_stack_new_' + method + '_done', 'orientation_count': 3, 'elapsed_ms': 1}
                   for method in ('iterative', 'single_parametric', 'summed', 'summed_function_map')]
        events += [{'stage': 'evaluator_stack_new_done', 'orientation_count': 3}]
        (logs / 'gammalog-test.jsonl').write_text(''.join(json.dumps(event) + '\n' for event in events))
else:
    card = tomllib.loads(Path(args[-1]).read_text())
    index = 0
    for command in card['commands']:
        tokens = shlex.split(command)
        if tokens[0] != 'bench': continue
        corpus, config = fixture['corpus'], fixture['config']
        value = {'process_id': corpus['process_id'], 'integrand_name': corpus['integrand'],
                 'point': corpus['points'][index % len(corpus['points'])], 'discrete_dim': corpus['discrete_dim'],
                 'target_seconds': config['bench_duration_seconds'], 'warmup_samples': 10, 'warmup_seconds': 0.001,
                 'minimal_integrand': False, 'momentum_space': False, 'n_batches': config['bench_batches'],
                 'total_samples': 100, 'displayed_result': {'re': 1.0, 'im': 0.5},
                 'batches': [{'total': 0.001}] * config['bench_batches'],
                 'summary': [{'category': 'Total', 'mean_seconds_per_sample': 0.001}]}
        Path(tokens[tokens.index('--json-output') + 1]).write_text(json.dumps(value))
        index += 1
"""
        )
        binary.chmod(0o755)
        arguments = [
            "run_performance.py",
            "--input-root",
            str(self.repository),
            "--baseline-binary",
            str(binary),
            "--candidate-binary",
            str(binary),
            "--baseline-commit",
            "a" * 40,
            "--candidate-commit",
            "b" * 40,
            "--work-dir",
            str(self.directory / "private"),
            "--output-dir",
            str(self.directory / "public"),
        ]
        main_cases = {
            name: value
            for name, value in producer.Fixture.cases(self.suite).items()
            if value[0]["family"] != "precise" and not value[0].get("comparison")
        }
        self.assertEqual(len(main_cases), 13)
        with (
            patch.object(sys, "argv", arguments),
            patch.object(producer.Fixture, "cases", return_value=main_cases),
        ):
            self.assertEqual(producer.main(), 0)
        for case in main_cases:
            snapshot = self.directory / "public" / case / "discovery"
            verify(snapshot)
            self.assertEqual(len(load(snapshot)["records"]), 10)
            self.assertEqual(
                len(list((snapshot / "records").glob("*-generation.json"))), 10
            )
            self.assertNotIn(
                "PRIVATE-SYNTHETIC-CLI-OUTPUT",
                "".join(p.read_text() for p in snapshot.rglob("*.json")),
            )

    def test_suite_finishes_all_cases_and_propagates_any_failure(self):
        names = [
            name
            for name, (workload, _) in producer.Fixture.cases(self.suite).items()
            if workload["family"] != "precise"
        ]
        self.assertEqual(len(names), 14)
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

    def test_ltd_cff_fixture_changes_only_representation_and_keeps_scalar_budget(self):
        name = "gl04-q1-squared-ltd-vs-cff-4d"
        fixtures = [
            producer.Fixture.from_suite(self.suite, name, self.repository, rep)
            for rep in ("cff", "ltd")
        ]
        cards = [producer.tomllib.loads(fixture.card) for fixture in fixtures]
        for fixture, card, rep in zip(fixtures, cards, ("cff", "ltd")):
            self.assertEqual(fixture.timeout, 240)
            self.assertEqual(
                fixture.config["correctness_tolerance"],
                {"absolute": 0.0, "relative": 1e-9},
            )
            self.assertEqual(len(fixture.corpus["points"]), 4)
            generation = card["cli_settings"]["global"]["generation"]
            self.assertEqual(generation.pop("three_dimensional_representations"), [rep])
            levels = card["default_runtime_settings"]["stability"]["levels"]
            self.assertEqual(len(levels), 1)
            self.assertEqual(levels[0]["precision"], "Double")
            self.assertEqual(levels[0].pop("three_dimensional_representation"), rep)
            self.assertFalse(
                card["default_runtime_settings"]["general"]["enable_cache"]
            )
            for setting in (
                "subtract_uv=true",
                "generate_integrated=true",
                "local_uv_cts_from_expanded_4d_integrands=true",
                "enable_thresholds=true",
                "summed_function_map=false",
            ):
                self.assertIn(setting, card["commands"][2])
        self.assertEqual(cards[0], cards[1])
        self.assertNotEqual(
            fixtures[0].hashes["fixture_sha256"], fixtures[1].hashes["fixture_sha256"]
        )
        for key in ("shared_fixture_sha256", "corpus_sha256", "config_sha256"):
            self.assertEqual(fixtures[0].hashes[key], fixtures[1].hashes[key])
        old = producer.Fixture.from_suite(
            self.suite, "gl04-q1-squared-projected-4d", self.repository
        )
        self.assertEqual(old.corpus, fixtures[0].corpus)
        self.assertIsNone(old.representation)
        self.assertNotIn("comparison", old.config)
        with self.assertRaises(ValueError):
            producer.Fixture.from_suite(self.suite, "gl15", self.repository, "ltd")

    def test_ltd_preparation_pins_one_candidate_and_keeps_the_scoped_policy(self):
        args = Namespace(
            baseline_commit="a" * 40,
            candidate_commit="b" * 40,
            baseline_binary=Path("never-run-old-main-for-ltd"),
            candidate_binary=Path(sys.executable),
            input_root=self.repository,
            work_dir=self.directory / "private",
            output_dir=self.directory / "public",
        )
        with patch("gammabench.gate.run", return_value=0) as run:
            self.assertEqual(
                producer.run_case(args, self.suite, "gl04-q1-squared-ltd-vs-cff-4d"), 0
            )
        definition = json.loads(run.call_args.args[0].read_text())
        policy = json.loads(run.call_args.kwargs["policy_path"].read_text())
        self.assertEqual(
            definition["tools"]["baseline"], definition["tools"]["candidate"]
        )
        self.assertEqual(
            definition["gate"]["baseline"], definition["gate"]["candidate"]
        )
        self.assertEqual(definition["gate"]["baseline"]["commit"], "b" * 40)
        self.assertEqual(policy["budget_seconds"], 240)
        self.assertEqual(policy["pairs"], 5)
        self.assertEqual(
            policy["thresholds"],
            {
                "generation_seconds": [0.5, 0.5],
                "warm_runtime_seconds_per_sample.SingleParametric": [0.5, 1e-6],
            },
        )
        self.assertEqual(
            policy["correctness_tolerance"], {"absolute": 0, "relative": 1e-9}
        )
        self.assertEqual(
            [c["index"] for c in policy["certified_zero"]["components"]], [1, 3, 5, 7]
        )
        for task in definition["tasks"]:
            self.assertEqual(
                task["argv"][task["argv"].index("--side") + 1], task["variant"]
            )
        fixture = json.loads(
            (run.call_args.args[0].parent / "fixture.json").read_text()
        )
        self.assertIn("$representation", fixture["card"])
        self.assertNotEqual(
            fixture["representation_hashes"]["cff"],
            fixture["representation_hashes"]["ltd"],
        )
        outputs = []
        for side, changed_inventory in [
            ("baseline", False),
            ("candidate", False),
            ("candidate", True),
        ]:
            work = self.directory / "capture" / "discovery" / str(len(outputs)) / side
            work.mkdir(parents=True)
            inventory = {"physical_cut": [1, 3] if changed_inventory else [1, 2]}
            metadata = {
                "representation_evidence": [
                    {"physical_inventory": inventory, "native_residue_count": 7}
                    for _ in fixture["corpus"]["points"]
                ]
            }
            arguments = [
                "adapter.py",
                "sample",
                "--fixture",
                str(run.call_args.args[0].parent / "fixture.json"),
                "--binary",
                sys.executable,
                "--work-dir",
                str(work),
                "--result",
                str(work / "result.json"),
                "--deadline",
                str(time.monotonic() + 240),
                "--side",
                side,
                "--commit",
                "b" * 40,
            ]
            with (
                patch.object(sys, "argv", arguments),
                patch.object(
                    producer,
                    "measure",
                    return_value=(
                        {
                            "generation_seconds": 1.0,
                            "warm_runtime_seconds_per_sample.SingleParametric": 0.001,
                        },
                        {"SingleParametric": [1.0, 1e-30] * 4},
                        metadata,
                    ),
                ),
            ):
                if changed_inventory:
                    with self.assertRaisesRegex(ValueError, "inventories differ"):
                        producer.sample_main()
                else:
                    self.assertEqual(producer.sample_main(), 0)
                    output = json.loads((work / "result.json").read_text())
                    self.assertEqual(
                        output["certificate"], policy["certified_zero"]["certificate"]
                    )
                    outputs.append(output)
        self.assertEqual(
            outputs[0]["metadata"]["physical_inventory_sha256_words"],
            outputs[1]["metadata"]["physical_inventory_sha256_words"],
        )

    def test_precise_preparation_stays_inside_one_deadline_and_is_reused_per_role(self):
        fixture = producer.Fixture.from_suite(
            self.suite, "paper-massless-soft-cff4d-arb", self.repository
        )
        prepared = self.directory / "suite"
        for relative in fixture.corpus["inputs"]:
            target = prepared / "inputs" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(self.repository / relative, target)
        fixture_path = prepared / "fixture.json"
        fixture_path.write_text(
            json.dumps(
                {
                    "card": fixture.card,
                    "corpus": fixture.corpus,
                    "config": fixture.config,
                }
            )
        )
        deadline = time.monotonic() + 7200
        paths = []

        def generate(binary, repository, directory, selected, actual_deadline):
            self.assertEqual(actual_deadline, deadline)
            directory.mkdir(parents=True)
            return (
                directory / "state",
                123.0,
                {"evaluator_count": 1, "generation_peak_ram_bytes": 1024},
            )

        def measure(binary, repository, directory, selected, actual_deadline, state):
            self.assertEqual(actual_deadline, deadline)
            paths.append(state)
            return (
                {
                    "evaluator_seconds_per_sample": 1.0,
                    "warm_runtime_seconds_per_sample.SingleParametric": 1.1,
                },
                {"SingleParametric": ["1e-411", "3e-158"]},
                {
                    "load_seconds": 0.1,
                    "model_warmup_seconds": 0.1,
                    "canonical_permutation": [0, 1, 2],
                    "point_bits": fixture.ray["canonical_point_bits"],
                },
            )

        results = []
        with (
            patch.object(producer, "generate", side_effect=generate) as generation,
            patch.object(producer, "measure", side_effect=measure),
        ):
            for repeat in range(2):
                for side, commit in [("baseline", "a" * 40), ("candidate", "b" * 40)]:
                    work = self.directory / "private" / "discovery" / str(repeat) / side
                    work.mkdir(parents=True)
                    output = work / "result.json"
                    arguments = [
                        "adapter.py",
                        "sample",
                        "--fixture",
                        str(fixture_path),
                        "--binary",
                        sys.executable,
                        "--work-dir",
                        str(work),
                        "--result",
                        str(output),
                        "--deadline",
                        str(deadline),
                        "--side",
                        side,
                        "--commit",
                        commit,
                        "--producer",
                        sys.executable,
                    ]
                    with patch.object(sys, "argv", arguments):
                        self.assertEqual(producer.sample_main(), 0)
                    results.append(json.loads(output.read_text()))
            self.assertEqual(generation.call_count, 2)
        self.assertEqual(paths[0], paths[2])
        self.assertEqual(paths[1], paths[3])
        self.assertNotEqual(paths[0][0], paths[1][0])
        self.assertEqual(
            [bool(r["metadata"]["preparation"]) for r in results],
            [True, True, False, False],
        )
        for result in results:
            self.assertEqual(
                result["outputs"]["SingleParametric"], ["1e-411", "3e-158"]
            )
            self.assertEqual(
                result["certificate"],
                fixture.qualification_policy()["certified_zero"]["certificate"],
            )

    def test_precise_suite_declares_slow_budget_and_native_clocks(self):
        from gammabench.gate import validate_policy

        args = Namespace(
            baseline_commit="a" * 40,
            candidate_commit="b" * 40,
            baseline_binary=Path(sys.executable),
            candidate_binary=Path(sys.executable),
            baseline_precise_producer=Path(sys.executable),
            candidate_precise_producer=Path(sys.executable),
            input_root=self.repository,
            work_dir=self.directory / "private",
            output_dir=self.directory / "public",
        )
        with patch("gammabench.gate.run", return_value=0) as run:
            self.assertEqual(
                producer.run_case(args, self.suite, "paper-massless-soft-cff4d-arb"), 0
            )
        definition = json.loads(run.call_args.args[0].read_text())
        policy = json.loads(run.call_args.kwargs["policy_path"].read_text())
        self.assertEqual(validate_policy(policy, definition), policy)
        self.assertEqual(
            (
                policy["budget_class"],
                policy["budget_seconds"],
                policy["correctness_format"],
            ),
            ("slow", 7200, "decimal"),
        )
        self.assertEqual(
            policy["correctness_tolerance"], {"absolute": 0, "relative": 1e-5}
        )
        self.assertNotIn("generation_seconds", policy["thresholds"])
        self.assertIs(policy["confirm_on_correctness_failure"], True)
        self.assertEqual(
            policy["thresholds"]["evaluator_seconds_per_sample"], [0.15, 1e-6]
        )
        for task in definition["tasks"]:
            self.assertNotIn("process_measurements", task)
            self.assertIn(task["variant"] + "-producer", task["tool_ids"])

    def test_ltd_timed_identity_requires_every_sample_and_native_threshold_inventory(
        self,
    ):
        self.fixture = producer.Fixture.from_suite(
            self.suite, "gl04-q1-squared-ltd-vs-cff-4d", self.repository, "ltd"
        )
        threshold = {
            "esurface_id": 12,
            "cut_boundary_edge_ids": [1, 2],
            "threshold_boundary_edge_ids": [4, 5, 6],
            "status": "can_become_pinched",
            "invariant_bound_is_applicable": True,
        }
        cut = {
            "cut_id": 0,
            "edge_ids": [1, 2],
            "raising_power": 1,
            "left_thresholds": [],
            "right_thresholds": [
                threshold,
                threshold | {"threshold_boundary_edge_ids": [4, 8]},
            ],
        }
        group = {
            "complete_residue_sum": True,
            "graphs": [{"name": "GL04", "native_residue_counts": [["ltd", 7]]}],
            "cuts": [
                cut,
                cut
                | {
                    "cut_id": 1,
                    "edge_ids": [1, 7, 8],
                    "right_thresholds": [
                        threshold
                        | {
                            "cut_boundary_edge_ids": [8],
                            "threshold_boundary_edge_ids": [5, 6],
                            "status": "potentially_existing",
                        }
                    ],
                },
                cut | {"cut_id": 2, "edge_ids": [3, 4], "right_thresholds": []},
            ],
            "loop_momentum_bases": [{"edge_ids": [1, 3, 5]}],
            "threshold_esurfaces": [{"edge_ids": [4, 5]}],
        }
        attempt = {
            "precision": "Double",
            "three_dimensional_representation": "ltd",
            "status": {"Stable": 2},
        }
        value = self.value | {
            "integrand_name": "numerator",
            "point": self.fixture.corpus["points"][0],
            "discrete_dim": [],
            "target_seconds": 2.0,
            "batches": [
                {
                    "total": 0.001,
                    "attempt_signatures": [{"samples": 5, "attempts": [attempt]}],
                    "displayed_result_bounds": [[1.0, 1.0], [0.5, 0.5]],
                    "invalid_result_samples": 0,
                }
                for _ in range(10)
            ],
            "integrand_info": {
                "process_id": 0,
                "integrand_name": "numerator",
                "kind": "cross_section",
                "generated_representations": ["ltd"],
                "generation_compilation": "Eager",
                "active_f64_backend": "eager",
                "graph_count": 1,
                "graph_group_count": 1,
                "graph_groups": [group],
            },
        }
        self.assertEqual(self.read_result(value), (0.001, [1.0, 0.5], 50))
        # Correct warmup and stable attempts must not hide a poisoned timed row.
        for bounds in (
            None,
            [],
            [[1.0]],
            [[1.0, 1.0], [0.5, 0.5, 0.5]],
            [[True, 1.0], [0.5, 0.5]],
            [[1.0, math.nan], [0.5, 0.5]],
            [[1.0, 1.0], [0.5, math.inf]],
            [[1.0, 1.0], [0.6, 0.5]],
            [[0.0, 1.0], [0.5, 0.5]],
            [[1.0, 1.01], [0.5, 0.5]],
            [[1.0, 1.0], [-0.5, 0.5]],
        ):
            changed = copy.deepcopy(value)
            changed["batches"][0]["displayed_result_bounds"] = bounds
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                self.read_result(changed)
        for invalid in (None, True, -1, 1, 0.0):
            changed = copy.deepcopy(value)
            changed["batches"][0]["invalid_result_samples"] = invalid
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                self.read_result(changed)
        for key in ("displayed_result_bounds", "invalid_result_samples"):
            changed = copy.deepcopy(value)
            del changed["batches"][0][key]
            with self.assertRaises(ValueError):
                self.read_result(changed)
        changed = copy.deepcopy(value)
        changed["batches"][0]["displayed_result_bounds"] = [
            [1.0 - 1e-10, 1.0 + 1e-10],
            [0.5, 0.5],
        ]
        self.assertEqual(self.read_result(changed), (0.001, [1.0, 0.5], 50))
        # A zero component has no allowance from the other component's scale.
        changed = copy.deepcopy(value)
        changed["displayed_result"]["im"] = 0.0
        for batch in changed["batches"]:
            batch["displayed_result_bounds"][1] = [0.0, 0.0]
        self.assertEqual(self.read_result(changed), (0.001, [1.0, 0.0], 50))
        changed["batches"][0]["displayed_result_bounds"][1][1] = 1e-300
        with self.assertRaises(ValueError):
            self.read_result(changed)
        for attempts in (
            [],
            [attempt | {"precision": "Quad"}],
            [attempt | {"three_dimensional_representation": "cff"}],
            [attempt, attempt],
            [attempt | {"status": {"Unstable": 2}}],
        ):
            changed = copy.deepcopy(value)
            changed["batches"][0]["attempt_signatures"][0]["attempts"] = attempts
            with self.subTest(attempts=attempts), self.assertRaises(ValueError):
                self.read_result(changed)
        for count in (0, True, 4):
            changed = copy.deepcopy(value)
            changed["batches"][0]["attempt_signatures"][0]["samples"] = count
            with self.assertRaises(ValueError):
                self.read_result(changed)
        changed = copy.deepcopy(value)
        changed["batches"][0]["attempt_signatures"][0]["samples"] = 4
        changed["batches"][1]["attempt_signatures"][0]["samples"] = 6
        with self.assertRaises(ValueError):
            self.read_result(changed)
        changed = copy.deepcopy(value)
        changed["total_samples"] = 53
        for batch in changed["batches"][:3]:
            batch["attempt_signatures"][0]["samples"] = 6
        self.assertEqual(self.read_result(changed), (0.001, [1.0, 0.5], 53))
        # The remainder belongs to the first batches, as in the native splitter.
        changed["batches"][0]["attempt_signatures"][0]["samples"] = 5
        changed["batches"][-1]["attempt_signatures"][0]["samples"] = 6
        with self.assertRaises(ValueError):
            self.read_result(changed)
        for key, wrong in (
            ("generated_representations", ["cff", "ltd"]),
            ("active_f64_backend", "symjit"),
        ):
            changed = copy.deepcopy(value)
            changed["integrand_info"][key] = wrong
            with self.assertRaises(ValueError):
                self.read_result(changed)
        for mutation in (
            "rows",
            "thresholds",
            "power",
            "cut",
            "survivor_threshold",
            "directive",
            "support",
        ):
            changed = copy.deepcopy(value)
            current = changed["integrand_info"]["graph_groups"][0]
            if mutation == "rows":
                current["graphs"][0]["native_residue_counts"] = [["ltd", 0]]
            elif mutation == "power":
                current["cuts"][0]["raising_power"] = 2
            elif mutation == "cut":
                current["cuts"][0]["edge_ids"] = [2, 7]
            elif mutation == "survivor_threshold":
                current["cuts"][2]["right_thresholds"] = [threshold]
            elif mutation == "directive":
                current["graphs"][0]["threshold_counterterm_directives"] = [
                    {"multiplier": "i"}
                ]
            elif mutation == "support":
                current["cuts"][0]["right_thresholds"][0][
                    "threshold_boundary_edge_ids"
                ] = [4, 5]
            else:
                for c in current["cuts"]:
                    c["right_thresholds"] = []
            with self.assertRaises(ValueError):
                self.read_result(changed)

    def test_gl04_zero_certificate_rejects_changed_source_and_configuration(self):
        name = "gl04-q1-squared-ltd-vs-cff-4d"
        fixture = producer.Fixture.from_suite(self.suite, name, self.repository)
        self.assertEqual(
            fixture.certified_imaginary_zero,
            producer.Fixture.GL04_IMAGINARY_ZERO,
        )
        other = copy.deepcopy(self.suite)
        other["representation_comparisons"]["other-ltd"] = {
            "case": "gl16-q7-squared-projected-4d",
            "comparison": "ltd_vs_cff",
        }
        strict = producer.Fixture.from_suite(other, "other-ltd", self.repository, "ltd")
        self.assertIsNone(strict.certified_imaginary_zero)

        for section, key, value in (
            ("scalar", "bench_duration_seconds", 3.0),
            ("scalar", "correctness_tolerance", {"absolute": 1e-8, "relative": 1e-9}),
        ):
            changed = copy.deepcopy(self.suite)
            changed[section][key] = value
            with self.assertRaises(ValueError):
                producer.Fixture.from_suite(changed, name, self.repository)
        changed = copy.deepcopy(self.suite)
        changed["scalar"]["corpus"]["points"][0][0] = 0.123
        with self.assertRaises(ValueError):
            producer.Fixture.from_suite(changed, name, self.repository)
        paths = ["benchmarks/performance/scalar.toml", *fixture.corpus["inputs"]]
        for relative in paths:
            target = self.directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((self.repository / relative).read_bytes())
        for relative in paths:
            target = self.directory / relative
            original = target.read_bytes()
            target.write_bytes(original + b"\n")
            with self.subTest(source=relative), self.assertRaises(ValueError):
                producer.Fixture.from_suite(self.suite, name, self.directory)
            target.write_bytes(original)

    def test_gl04_zero_certificate_checks_saved_phase_defaults_after_timing(self):
        fixture = producer.Fixture.from_suite(
            self.suite, "gl04-q1-squared-ltd-vs-cff-4d", self.repository, "ltd"
        )
        # Minimal emitted metadata: checks must bind native defaults as well as
        # the explicitly frozen input card, including restricted model widths.
        global_settings = """
[global.generation.threshold_subtraction]
enable_thresholds = true
[global.generation.uv]
softct = true
subtract_uv = true
generate_integrated = true
local_uv_cts_from_expanded_4d_integrands = true
[global.generation.uv.renormalization_prescription]
log_divergent = "MUV"
massive_power_divergent = "MUV"
massless_power_divergent = "MUV"
overrides = []
[global.generation.uv.vakint]
normalization = "MSbar"
additional_normalization = "1"
evaluation_methods = ["alphaloop", "matad", "fmft"]
"""
        runtime = """
[general]
m_uv = 20.0
mu_r = 3.0
renormalization_localization_scale = 1000.0
numerator_sampling_scale = 1.0
additional_param_values = []
disable_flux_factor = false
integral_unit = "auto"
[model]
[kinematics]
e_cm = 1.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[1.0, 0.0, 0.0, 0.0]]
helicities = [0]
[subtraction]
disable_threshold_subtraction = false
"""
        phase_values = {
            "ZERO": [0.0, 0.0],
            "lam": [1.0, 0.0],
            "mass_scalar_1": [0.1, 0.0],
            "width_scalar_1": [0.0, 0.0],
            "mass_scalar_2": [2.0, 0.0],
            "width_scalar_2": [0.0, 0.0],
        }
        model = {
            "parameters": [
                {"name": name, "value": value, "parameter_type": "real"}
                for name, value in phase_values.items()
            ],
            "couplings": [{"name": "SCALAR_COUPLING", "value": [0.0, -1.0]}],
        }
        files = {
            "global_settings.toml": global_settings,
            "default_runtime_settings.toml": runtime,
            "processes/cross_sections/scalar_gl04/numerator/integrand/settings.toml": "[model]"
            + runtime.split("[model]", 1)[1],
            "model.json": json.dumps(model),
            "model_parameters.json": json.dumps(
                {
                    name: phase_values[name]
                    for name in ("lam", "mass_scalar_1", "mass_scalar_2")
                }
            ),
            "generation_summary.json": json.dumps(
                {
                    "peak_ram_bytes": 1,
                    "reports": [
                        {
                            "graph_name": "GL04",
                            "stats": {
                                "evaluator_count": 1,
                                **{
                                    stage: {"secs": 0, "nanos": 1}
                                    for stage in (
                                        "total_time",
                                        "evaluator_spenso_time",
                                        "evaluator_symbolica_time",
                                        "evaluator_compile_time",
                                    )
                                },
                            },
                        }
                    ],
                }
            ),
        }
        mutations = [None]
        mutations.extend(
            (
                "global_settings.toml",
                global_settings,
                global_settings + "\n[global.generation.medium]\n" + medium,
            )
            for medium in [
                'mode = "thermodynamic_equilibrium"',
                'mode = "vacuum"\nvacuum_subtraction = true',
            ]
        )
        mutations.extend(
            ("global_settings.toml", old, new)
            for old, new in (
                ("softct = true", "softct = false"),
                ('"MSbar"', '"pySecDec"'),
                ('additional_normalization = "1"', 'additional_normalization = "i"'),
                ('log_divergent = "MUV"', 'log_divergent = "IR"'),
            )
        )
        mutations.extend(
            (path, old, new)
            for path in (
                "default_runtime_settings.toml",
                "processes/cross_sections/scalar_gl04/numerator/integrand/settings.toml",
            )
            for old, new in (
                ("e_cm = 1.0", "e_cm = 2.0"),
                ("[[1.0,", "[[2.0,"),
            )
        )
        mutations.extend(
            ("default_runtime_settings.toml", old, new)
            for old, new in (
                ('integral_unit = "auto"', 'integral_unit = "pb"'),
                ("m_uv = 20.0", "m_uv = -20.0"),
            )
        )
        path = "processes/cross_sections/scalar_gl04/numerator/integrand/settings.toml"
        mutations.append(
            (path, files[path], files[path] + "\n[general]\nm_uv = -20.0\n")
        )
        for name in ("lam", "mass_scalar_1", "width_scalar_1"):
            changed = copy.deepcopy(model)
            next(
                parameter
                for parameter in changed["parameters"]
                if parameter["name"] == name
            )["value"][1] = 0.1
            mutations.append(("model.json", files["model.json"], json.dumps(changed)))
        mutations.append(
            ("model_parameters.json", '"lam": [1.0, 0.0]', '"lam": [2.0, 0.0]')
        )
        external_width = json.loads(files["model_parameters.json"]) | {
            "width_scalar_1": [0.1, 0.0]
        }
        mutations.append(
            (
                "model_parameters.json",
                files["model_parameters.json"],
                json.dumps(external_width),
            )
        )
        for index, mutation in enumerate(mutations):
            directory = self.directory / str(index)
            contents = files.copy()
            if mutation:
                path, old, new = mutation
                self.assertIn(old, contents[path])
                contents[path] = contents[path].replace(old, new)

            def generate(*_, contents=contents, directory=directory, logfile_filter):
                self.assertEqual(logfile_filter, "off")
                for relative, text in contents.items():
                    path = directory / "state" / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(text)

            with (
                self.subTest(mutation=mutation),
                patch.object(producer, "run_cli", side_effect=generate),
            ):
                arguments = (
                    Path(sys.executable),
                    self.repository,
                    directory,
                    fixture,
                    time.monotonic() + 10,
                )
                if mutation:
                    with self.assertRaises(ValueError):
                        producer.generate(*arguments)
                else:
                    with patch.object(
                        producer.time, "monotonic", side_effect=[10.0, 12.0]
                    ):
                        _, seconds, _ = producer.generate(*arguments)
                        self.assertEqual(seconds, 10.0)

    def test_slow_family_keeps_fast_limits_and_original_soft_prescriptions(self):
        fixture = producer.Fixture.from_suite(
            self.suite, "paper-massless-soft-cff4d-arb", self.repository
        )
        self.assertEqual(fixture.timeout, 7200)
        self.assertEqual(
            fixture.config["required_metrics"],
            ["evaluator_seconds_per_sample", "warm_runtime_seconds_per_sample"],
        )
        original = producer.tomllib.loads(
            (
                self.repository
                / "tests/resources/run_cards/paper_figure_b1_double_triangle_soft_ir.toml"
            ).read_text()
        )
        card = producer.tomllib.loads(fixture.card)
        expected = original["cli_settings"]["global"]["generation"]["uv"] | {
            "local_uv_cts_from_expanded_4d_integrands": True,
            "softct": True,
            "subtract_uv": True,
            "final_integrand": "ThreeD",
        }
        self.assertEqual(card["cli_settings"]["global"]["generation"]["uv"], expected)
        self.assertEqual(
            card["cli_settings"]["global"]["generation"]["evaluator"],
            original["cli_settings"]["global"]["generation"]["evaluator"]
            | {
                "inline_numerator_functions": True,
                "summed": False,
                "summed_function_map": False,
                "compile": False,
                "store_atom": False,
            },
        )
        self.assertEqual(
            card["default_runtime_settings"]["kinematics"],
            original["default_runtime_settings"]["kinematics"],
        )
        self.assertEqual(
            card["default_runtime_settings"]["general"],
            original["default_runtime_settings"]["general"] | {"enable_cache": False},
        )
        self.assertFalse(
            card["cli_settings"]["global"]["generation"]["threshold_subtraction"][
                "enable_thresholds"
            ]
        )
        with self.assertRaises(ValueError):
            producer.Fixture(
                self.fixture.card,
                self.fixture.corpus,
                self.fixture.config | {"fixture_timeout_seconds": 7200},
                self.repository,
            )

    def test_precise_producer_contract_and_evaluator_timing(self):
        fixture = producer.Fixture.from_suite(
            self.suite, "paper-massless-soft-cff4d-arb", self.repository
        )
        evaluation = {
            "output": ["-5.33e-411", "-3.31e-158"],
            "evaluator_seconds": 38.0,
            "physical_seconds": 38.5,
            "attempts": [{"precision": "Arb", "status": {"Stable": 2}}],
        }
        value = {
            "commit": "a" * 40,
            "process_id": 0,
            "integrand": fixture.corpus["integrand"],
            "precision_bits": 1000,
            "certified_real_zero": producer.Fixture.PAPER_REAL_ZERO.copy(),
            "point_bits": [fixture.ray["canonical_point_bits"][i] for i in [1, 0, 2]],
            "canonical_point_bits": fixture.ray["canonical_point_bits"],
            "canonical_permutation": [1, 0, 2],
            "load_seconds": 2.0,
            "model_warmup_seconds": 0.01,
            "warmup_samples": 1,
            "measured_samples": 1,
            "evaluations": [
                evaluation | {"evaluator_seconds": 39.0, "physical_seconds": 40.0},
                evaluation,
            ],
        }
        path = self.directory / "precise.json"
        path.write_text(json.dumps(value))
        timing, output, _ = producer.precise_result(path, fixture, "a" * 40)
        self.assertEqual(timing, 38.0)
        self.assertEqual(output, evaluation["output"])
        with patch.object(
            producer,
            "run_cli",
            side_effect=lambda command, *_: Path(command[-1]).write_text(
                json.dumps(value)
            ),
        ):
            metrics, measured_output, _ = producer.measure(
                Path(sys.executable),
                self.repository,
                self.directory / "measurement",
                fixture,
                time.monotonic() + 10,
                (self.directory / "state", Path(sys.executable), "a" * 40),
            )
        self.assertEqual(
            metrics,
            {
                "evaluator_seconds_per_sample": 38.0,
                "warm_runtime_seconds_per_sample.SingleParametric": 38.5,
            },
        )
        self.assertEqual(measured_output, {"SingleParametric": evaluation["output"]})
        # A legitimate catalogue/runtime ordering change preserves physical input.
        permuted = value | {
            "canonical_permutation": [2, 0, 1],
            "point_bits": [fixture.ray["canonical_point_bits"][i] for i in [2, 0, 1]],
        }
        path.write_text(json.dumps(permuted))
        self.assertEqual(producer.precise_result(path, fixture, "a" * 40)[0], 38.0)
        changes = [
            ("precision_bits", 53),
            ("certified_real_zero", None),
            (
                "certified_real_zero",
                producer.Fixture.PAPER_REAL_ZERO | {"component": 1},
            ),
            ("commit", "b" * 40),
            ("point_bits", [[0, 0, 0]]),
            ("canonical_permutation", [0, 1, 2]),
            ("model_warmup_seconds", 0),
            ("warmup_samples", 0),
            ("measured_samples", 2),
            ("evaluations", [evaluation]),
        ]
        for key, wrong in changes:
            with self.subTest(key=key):
                path.write_text(json.dumps(value | {key: wrong}))
                with self.assertRaises(ValueError):
                    producer.precise_result(path, fixture, "a" * 40)
        for changed in (
            evaluation | {"output": [0.0, "-3.31e-158"]},
            evaluation | {"output": ["-5.33e-411", "0"]},
            evaluation | {"output": ["NaN", "-3.31e-158"]},
            evaluation | {"attempts": []},
            evaluation
            | {"attempts": [{"precision": "Double", "status": {"Stable": 2}}]},
            evaluation | {"attempts": [{"precision": "Arb", "status": {"Stable": 1}}]},
            evaluation
            | {
                "attempts": [
                    {
                        "precision": "Arb",
                        "status": {"Stable": 2},
                        "three_dimensional_representation": "ltd",
                    }
                ]
            },
        ):
            path.write_text(json.dumps(value | {"evaluations": [changed, changed]}))
            with self.assertRaises(ValueError):
                producer.precise_result(path, fixture, "a" * 40)

    def test_explicit_lane_selection_never_counts_the_slow_case_as_fast_coverage(self):
        cases = producer.Fixture.cases(self.suite)
        for lane, count in (("fast", 14), ("slow", 1), ("all", 15)):
            directory = self.directory / lane
            arguments = [
                "run_performance.py",
                "--lane",
                lane,
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
                patch.object(producer, "run_case", return_value=0) as run,
            ):
                self.assertEqual(producer.main(), 0)
            selected = [call.args[2] for call in run.call_args_list]
            self.assertEqual(len(selected), count)
            self.assertEqual(
                "paper-massless-soft-cff4d-arb" in selected, lane != "fast"
            )
            self.assertTrue(set(selected) <= cases.keys())


if __name__ == "__main__":
    unittest.main()
