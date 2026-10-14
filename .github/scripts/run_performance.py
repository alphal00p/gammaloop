#!/usr/bin/env python3
"""Measure paired fresh generation and warmed fixed-point CLI benchmarks."""

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import signal
import statistics
import subprocess
import sys
import time
from decimal import Decimal, Underflow, localcontext
import tomllib
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from string import Template
from typing import ClassVar

EVALUATOR_METHODS = ("Iterative", "SingleParametric", "Summed", "SummedFunctionMap")
MODE_METRICS = {
    "evaluator_build_seconds": (0.20, 0.1),
    "warm_runtime_seconds_per_sample": (0.15, 1e-6),
}
THRESHOLDS = {
    "generation_seconds": (0.20, 0.5),
    "generation_peak_ram_bytes": (0.20, 64 * 1024**2),
    "evaluator_seconds_per_sample": (0.15, 1e-6),
    **{
        f"{metric}.{method}": limits
        for metric, limits in MODE_METRICS.items()
        for method in EVALUATOR_METHODS
    },
}


def positive_number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a number")
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be positive and finite")
    return value


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def check_deadline(deadline):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise ValueError("The performance harness exceeded its declared budget")
    return remaining


class Fixture:
    @classmethod
    def paper_real_zero_bound(cls, certificate):
        if certificate != cls.PAPER_REAL_ZERO:
            raise ValueError("Unexpected paper source-zero certificate")
        with localcontext() as context:
            context.prec = 1024
            context.traps[Underflow] = True
            # Same dimensioned sqrt(epsilon) rule as the certified scalar tests.
            # VarFloat<1000>::epsilon() = 2^-999; this uncut density has E^-7.
            return (Decimal(2) ** -999).sqrt() * Decimal(300) ** -7

    def qualification_policy(self):
        if self.precise:
            certificate = self.PAPER_REAL_ZERO
            components = [
                {
                    "index": 0,
                    "absolute_bound": str(self.paper_real_zero_bound(certificate)),
                }
            ]
        elif self.certified_imaginary_zero is not None:
            certificate = self.GL04_IMAGINARY_ZERO
            components = [
                {"index": index, "absolute_bound": math.sqrt(2.0**-52)}
                for index in (1, 3, 5, 7)
            ]
        else:
            return {}
        # The framework owns numerical qualification. This adapter binds the
        # physical certificate to the frozen source, complete card and corpus.
        witness = {
            "id": certificate["id"],
            "sha256": hashlib.sha256(
                json.dumps(
                    {
                        "certificate": certificate,
                        "fixture": self.hashes.get(
                            "shared_fixture_sha256", self.hashes["fixture_sha256"]
                        ),
                        "corpus": self.hashes["corpus_sha256"],
                        "config": self.hashes["config_sha256"],
                    },
                    sort_keys=True,
                ).encode()
            ).hexdigest(),
        }
        return {
            "correctness_format": "decimal" if self.precise else "number",
            "certified_zero": {"certificate": witness, "components": components},
        }

    # This is the source predicate in precise_sample.rs, not a user-selected
    # absolute tolerance or a classification inferred from numerical output.
    PAPER_REAL_ZERO: ClassVar[dict[str, str | int]] = {
        "id": "paper-b1-local-4d-real-v1",
        "component": 0,
        "precision_bits": 1000,
        "energy_dimension": -7,
        "e_cm": 300,
    }
    # Fixed GL04 source: q1² cancels the unique simple massless edge 1.
    # Both threshold-bearing cuts contain that edge; the surviving cut [3,4]
    # has no thresholds. Its scalar phase is real, including MUV/MSbar vacuum
    # addbacks (i per integrated loop) and the remaining contours/LU factor.
    # Fixture checks the saved phase settings and actual cut catalogue as well.
    GL04_IMAGINARY_ZERO: ClassVar[dict[str, str | int]] = {
        "id": "gl04-q1-squared-imaginary-v1",
        "component": 1,
        "precision_bits": 53,
        "energy_dimension": 3,
        "e_cm": 1,
        "point_count": 4,
    }
    GL04_SOURCE_HASHES: ClassVar[dict[str, str]] = {
        "shared_fixture_sha256": "c8cf76ae153e8e1ae9a88e719ff59f6a6f41b64c6e42b43a94ca620ccf1991ac",
        "corpus_sha256": "91546c7de9b78fac35c1f29dd55d7e08033c13d8e4c40b8719ae2ae1cfbb497c",
        "config_sha256": "a1c9331b0dc3de9cdfd1e895d0c39c6bf5fe99981ec041cc8acdd9f11e65dbef",
    }
    GL04_REPRESENTATION_HASHES: ClassVar[dict[str, str]] = {
        "cff": "b21a3046667d70b9449f7119d7d654a2c8085c320a91bb654d932dbe34989673",
        "ltd": "543e27d418d9ccc03f369da9a847a1964e3c46b5eaf0d375df4624733d5a3aaa",
    }

    @staticmethod
    def cases(suite):
        cases = {
            f"{name}-{route}" if route else name: (workload, route)
            for name, workload in suite["workloads"].items()
            for route in (suite["routes"] if workload["family"] == "scalar" else [None])
        }
        for name, comparison in suite.get("representation_comparisons", {}).items():
            workload, route = cases[comparison["case"]]
            if (
                name in cases
                or set(comparison) != {"case", "comparison"}
                or comparison["comparison"] != "ltd_vs_cff"
                or workload["family"] != "scalar"
                or route != "projected-4d"
            ):
                raise ValueError("Invalid representation comparison")
            cases[name] = (workload | {"comparison": "ltd_vs_cff"}, route)
        if not cases or any(
            not re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in cases
        ):
            raise ValueError("Invalid performance case names")
        return cases

    @staticmethod
    def source(input_root, name):
        path = Path(name)
        source = (input_root / path).resolve()
        if (
            path.is_absolute()
            or ".." in path.parts
            or not source.is_relative_to(input_root)
        ):
            raise ValueError("Fixture inputs must stay in the shared repository")
        return source

    @classmethod
    def from_suite(cls, suite, name, input_root, representation=None):
        workload, route = cls.cases(suite)[name]
        family = suite[workload["family"]]
        config = suite["settings"] | {
            key: value for key, value in family.items() if key not in ("card", "corpus")
        }
        config["family"] = workload["family"]
        corpus = family["corpus"] | {
            "process": workload["process"],
            "graph_name": workload["graph_name"],
        }
        directory = Path("benchmarks/performance")
        card = cls.source(input_root, str(directory / family["card"])).read_text()
        if route:
            config |= {
                key: suite["routes"][route][key]
                for key in ("evaluator_methods", "reference_method")
            }
            graph = str(directory / workload["graph"])
            corpus["inputs"] = [graph, *corpus["inputs"]]
            values = suite["routes"][route] | {
                "process": corpus["process"],
                "graph": graph,
            }
            card = Template(card).substitute(
                {
                    key: str(value).lower() if type(value) is bool else value
                    for key, value in values.items()
                }
            )
        if workload.get("comparison") == "ltd_vs_cff":
            if representation not in (None, "cff", "ltd"):
                raise ValueError("Unknown representation")
            config["comparison"] = "ltd_vs_cff"
            # The shared card and every physical setting stay identical. These
            # two explicit declarations are the only representation differences.
            card += '\n[cli_settings.global.generation]\nthree_dimensional_representations = ["$representation"]\n'
            card += """
[default_runtime_settings.stability]
check_on_norm = true
escalate_if_exact_zero = false
min_abs_wgt_for_escalation = 0.0
loop_momenta_norm_escalation_factor = -1.0
rotation_axis = [{ type = "euler_angles", alpha = 0.1, beta = 0.2, gamma = 0.3 }]
levels = [{ precision = "Double", three_dimensional_representation = "$representation", required_precision_for_re = 1e-5, required_precision_for_im = 1e-5, ecm_relative_tolerance_for_re = 0.0, ecm_relative_tolerance_for_im = 0.0, escalate_for_large_weight_threshold = 0.9 }]
"""
            shared = cls(card, corpus, config, input_root)
            selected = representation or "cff"
            fixture = cls(
                Template(card).substitute(representation=selected),
                corpus,
                config,
                input_root,
            )
            fixture.representation = selected
            fixture.hashes["shared_fixture_sha256"] = shared.hashes["fixture_sha256"]
            if name == "gl04-q1-squared-ltd-vs-cff-4d":
                if (
                    any(
                        fixture.hashes[key] != expected
                        for key, expected in Fixture.GL04_SOURCE_HASHES.items()
                    )
                    or fixture.hashes["fixture_sha256"]
                    != Fixture.GL04_REPRESENTATION_HASHES[selected]
                ):
                    raise ValueError(
                        "GL04 source-zero certificate requires the frozen source/card/model/corpus"
                    )
                fixture.certified_imaginary_zero = Fixture.GL04_IMAGINARY_ZERO.copy()
            fixture.shared_card = card
            return fixture
        if representation is not None:
            raise ValueError("Representation overrides require a declared comparison")
        return cls(card, corpus, config, input_root)

    def __init__(self, card, corpus, config, input_root):
        self.card = card
        self.input_root = input_root.resolve()
        self.corpus = corpus
        self.config = config
        self.precise = config.get("family") == "precise"
        self.representation = None
        self.certified_imaginary_zero = None
        self.shared_card = None
        fixture_hash = hashlib.sha256(card.encode())
        inputs = self.corpus["inputs"]
        if not isinstance(inputs, list) or any(not isinstance(p, str) for p in inputs):
            raise ValueError("Fixture inputs must be relative paths")
        if len(set(inputs)) != len(inputs):
            raise ValueError("Fixture inputs must be distinct")
        for input_file in inputs:
            source = self.source(self.input_root, input_file)
            fixture_hash.update(b"\0" + input_file.encode() + b"\0")
            fixture_hash.update(source.read_bytes())
        self.hashes = {
            "fixture_sha256": fixture_hash.hexdigest(),
            "corpus_sha256": hashlib.sha256(
                json.dumps(corpus, sort_keys=True, allow_nan=False).encode()
            ).hexdigest(),
            "config_sha256": hashlib.sha256(
                json.dumps(config, sort_keys=True, allow_nan=False).encode()
            ).hexdigest(),
        }
        card = tomllib.loads(card)
        if card["default_runtime_settings"]["general"]["enable_cache"] is not False:
            raise ValueError("The production benchmark must disable point caching")
        if self.config["schema_version"] != 2 or self.config["pairs"] != 5:
            raise ValueError("The performance gate requires schema 2 and five pairs")
        self.methods = self.config["evaluator_methods"]
        if (
            not isinstance(self.methods, list)
            or not self.methods
            or any(not isinstance(method, str) for method in self.methods)
            or len(set(self.methods)) != len(self.methods)
            or not set(self.methods) <= set(EVALUATOR_METHODS)
            or self.config["reference_method"] not in self.methods
        ):
            raise ValueError("Invalid evaluator methods or reference")
        metrics = self.config["required_metrics"]
        if (
            not isinstance(metrics, list)
            or any(not isinstance(metric, str) for metric in metrics)
            or len(set(metrics)) != len(metrics)
            or not (
                {"evaluator_seconds_per_sample", "warm_runtime_seconds_per_sample"}
                if self.precise
                else {"generation_seconds", "warm_runtime_seconds_per_sample"}
            )
            <= set(metrics)
            or not set(metrics)
            <= {
                "generation_seconds",
                "generation_peak_ram_bytes",
                "evaluator_seconds_per_sample",
                *MODE_METRICS,
            }
        ):
            raise ValueError("Unexpected required performance metrics")
        self.metrics = [
            f"{metric}.{method}" if method else metric
            for metric in metrics
            for method in (self.methods if metric in MODE_METRICS else [None])
        ]
        if len(self.methods) > 1 and (
            card["default_runtime_settings"]["sampling"]["orientations"] != "summed"
        ):
            raise ValueError("Multiple evaluator methods require summed orientations")
        self.timeout = positive_number(
            self.config["fixture_timeout_seconds"], "Fixture timeout"
        )
        if self.timeout > (7200 if self.precise else 240):
            raise ValueError("The fixture timeout exceeds its family's budget")
        for label in ("process", "integrand", "graph_name"):
            if not re.fullmatch(r"[A-Za-z0-9_]+", self.corpus[label]):
                raise ValueError(f"Invalid corpus {label}")
        if type(self.corpus["process_id"]) is not int or self.corpus["process_id"] < 0:
            raise ValueError("Invalid process id")
        if self.precise:
            if self.methods != ["SingleParametric"] or metrics != [
                "evaluator_seconds_per_sample",
                "warm_runtime_seconds_per_sample",
            ]:
                raise ValueError(
                    "Precise measurements require the native SingleParametric metrics"
                )
            if self.config.get("correctness_format") != "decimal":
                raise ValueError(
                    "Precise measurements require decimal correctness values"
                )
            self.certificate = self.source(self.input_root, self.corpus["certificate"])
            if self.corpus["certificate"] not in inputs:
                raise ValueError("The physical certificate must be fingerprinted")
            self.ray = json.loads(self.certificate.read_text())
            if (
                self.ray["precision_bits"] != 1000
                or self.ray["process"] != self.corpus["process"]
            ):
                raise ValueError("Unexpected physical certificate")
            return
        positive_number(self.config["bench_duration_seconds"], "Bench duration")
        batches = self.config["bench_batches"]
        if isinstance(batches, bool) or not isinstance(batches, int) or batches < 2:
            raise ValueError("Bench requires at least two batches")
        dimension = self.corpus["continuous_dimension"]
        if type(dimension) is not int or dimension <= 0 or dimension % 3:
            raise ValueError("The corpus requires three coordinates per loop")
        dims = self.corpus["discrete_dim"]
        if not isinstance(dims, list) or any(type(v) is not int or v < 0 for v in dims):
            raise ValueError("Invalid discrete dimensions")
        points = self.corpus["points"]
        if len(points) < 2 or any(len(p) != dimension for p in points):
            raise ValueError(
                "The corpus requires multiple points of the declared dimension"
            )
        for point in points:
            if any(
                isinstance(v, bool)
                or not isinstance(v, (int, float))
                or not math.isfinite(v)
                or not 0 < v < 1
                for v in point
            ):
                raise ValueError("Invalid x-space point")
        if len({tuple(p) for p in points}) != len(points):
            raise ValueError("The fixed-point corpus must contain distinct points")

    def benchmark_card(self, directory):
        commands = []
        for method in self.methods:
            commands.append(
                f"set process -p {self.corpus['process']} -i {self.corpus['integrand']} "
                f"kv general.evaluator_method={method}"
            )
            for index, point in enumerate(self.corpus["points"]):
                point_args = " ".join(format(v, ".17g") for v in point)
                dims = self.corpus["discrete_dim"]
                selector = "-d " + " ".join(str(v) for v in dims) if dims else ""
                path = directory / f"bench-{method}-{index}.json"
                commands.append(
                    f"bench -p {self.corpus['process']} -i {self.corpus['integrand']} "
                    f"-x {point_args} {selector} "
                    f"--duration {self.config['bench_duration_seconds']} "
                    f"--n-batches {self.config['bench_batches']} "
                    f"--json-output {shlex.quote(str(path))}"
                )
        commands.append("quit --no-save-state")
        card = directory / "benchmark.toml"
        card.write_text("commands = " + json.dumps(commands, ensure_ascii=False) + "\n")
        return card

    def construction_metrics(self, state, deadline):
        stages = {
            "evaluator_stack_new_iterative_done": "Iterative",
            "evaluator_stack_new_single_parametric_done": "SingleParametric",
            "evaluator_stack_new_summed_done": "Summed",
            "evaluator_stack_new_summed_function_map_done": "SummedFunctionMap",
        }
        flags = {
            "iterative_orientation_optimization": "Iterative",
            "summed": "Summed",
            "summed_function_map": "SummedFunctionMap",
        }
        durations = dict.fromkeys(self.methods, 0.0)
        counts = dict.fromkeys(self.methods, 0)
        orientation_counts = []
        # Generation uses one worker. Validate each complete stack separately:
        # distinct cuts/dual groups legitimately repeat the same stage names.
        for path in sorted((state / "logs").glob("gammalog*.jsonl")):
            active = None
            with path.open() as stream:
                for line in stream:
                    check_deadline(deadline)
                    event = json.loads(line)
                    if not isinstance(event, dict):
                        raise TypeError("Invalid evaluator construction log record")
                    stage = event.get("stage")
                    if stage == "evaluator_stack_new_start":
                        if active is not None:
                            raise ValueError(
                                "Overlapping evaluator construction records"
                            )
                        orientations = event["orientation_count"]
                        if type(orientations) is not int or orientations <= 0:
                            raise ValueError("Invalid construction orientation count")
                        if any(type(event[flag]) is not bool for flag in flags):
                            raise ValueError("Invalid evaluator construction flags")
                        expected = {"SingleParametric"} | {
                            method for flag, method in flags.items() if event[flag]
                        }
                        active = (expected, set(), orientations)
                    elif stage in stages or stage == "evaluator_stack_new_done":
                        if active is None:
                            raise ValueError(
                                "Evaluator construction stage without a start"
                            )
                        if (
                            type(event["orientation_count"]) is not int
                            or event["orientation_count"] != active[2]
                        ):
                            raise ValueError(
                                "Construction orientation counts changed within a stack"
                            )
                        if stage == "evaluator_stack_new_done":
                            if active[0] != active[1]:
                                raise ValueError(
                                    "Incomplete evaluator construction records"
                                )
                            orientation_counts.append(active[2])
                            active = None
                            continue
                        method = stages[stage]
                        if method not in active[0] or method in active[1]:
                            raise ValueError(
                                "Unexpected or duplicate evaluator construction stage"
                            )
                        elapsed = (
                            positive_number(
                                event["elapsed_ms"], "Evaluator construction"
                            )
                            / 1000
                        )
                        active[1].add(method)
                        if method in durations:
                            durations[method] += elapsed
                            counts[method] += 1
            if active is not None:
                raise ValueError("Unfinished evaluator construction records")
        if any(count == 0 for count in counts.values()):
            raise ValueError("Missing requested evaluator construction timings")
        return (
            {
                f"evaluator_build_seconds.{method}": elapsed
                for method, elapsed in durations.items()
            },
            {
                "evaluator_build_counts": counts,
                "construction_orientation_counts": orientation_counts,
            },
        )


def run_cli(command, repository, log_path, deadline, *, logfile_filter="off"):
    check_deadline(deadline)
    environment = os.environ.copy()
    environment.pop("GL_ALL_LOG_FILTER", None)
    environment.update(
        {
            "RAYON_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "GL_LOGFILE_FILTER": logfile_filter,
        }
    )
    with log_path.open("w") as log:
        process = subprocess.Popen(
            command,
            cwd=repository,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=not os.environ.get("GAMMABENCH_PRIVATE_TRIAL"),
        )
        # POSIX wait(timeout=...) polls with sleeps of up to 50 ms. Wait for
        # the actual exit in a worker and enforce the shared deadline separately.
        with ThreadPoolExecutor(max_workers=1) as waiter:
            try:
                completion = waiter.submit(process.wait)
                status = completion.result(timeout=check_deadline(deadline))
            except TimeoutError:
                try:
                    os.kill(process.pid, signal.SIGTERM) if os.environ.get(
                        "GAMMABENCH_PRIVATE_TRIAL"
                    ) else os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    completion.result(timeout=5)
                except TimeoutError:
                    try:
                        os.kill(process.pid, signal.SIGKILL) if os.environ.get(
                            "GAMMABENCH_PRIVATE_TRIAL"
                        ) else os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    completion.result()
                raise ValueError("The performance harness exceeded its declared budget")
            finally:
                # A successful parent must not leave compiler or evaluator children.
                try:
                    os.kill(process.pid, signal.SIGKILL) if os.environ.get(
                        "GAMMABENCH_PRIVATE_TRIAL"
                    ) else os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
    if status:
        raise ValueError(f"CLI exited with status {status}; private log: {log_path}")


def benchmark_result(path, fixture, index, *, metadata=None):
    value = json.loads(path.read_text())
    if value["point"] != fixture.corpus["points"][index]:
        raise ValueError("Bench did not evaluate the requested corpus point")
    if value["discrete_dim"] != fixture.corpus["discrete_dim"]:
        raise ValueError("Bench changed the requested discrete dimensions")
    if (
        value["process_id"] != fixture.corpus["process_id"]
        or value["integrand_name"] != fixture.corpus["integrand"]
    ):
        raise ValueError("Bench changed the requested integrand")
    if value["target_seconds"] != fixture.config["bench_duration_seconds"]:
        raise ValueError("Bench changed the requested target duration")
    if value["warmup_samples"] != 10:
        raise ValueError("Bench did not complete the expected ten warmup samples")
    positive_number(value["warmup_seconds"], "Bench warmup duration")
    if value["minimal_integrand"] is not False or value["momentum_space"] is not False:
        raise ValueError("Bench must evaluate the production x-space path")
    if value["n_batches"] != fixture.config["bench_batches"]:
        raise ValueError("Bench did not complete all requested batches")
    if (
        type(value["total_samples"]) is not int
        or value["total_samples"] < value["n_batches"]
    ):
        raise ValueError("Bench did not complete a valid sample count")
    batches = value["batches"]
    if len(batches) != fixture.config["bench_batches"]:
        raise ValueError("Bench did not record all requested batches")
    batch_totals = [
        positive_number(batch["total"], "Bench batch Total") for batch in batches
    ]
    output = value["displayed_result"]
    numbers = [output["re"], output["im"]]
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in numbers):
        raise ValueError("Bench returned a nonfinite physical result")
    if math.hypot(*numbers) == 0:
        raise ValueError("The physical benchmark must not evaluate to zero")
    totals = [row for row in value["summary"] if row["category"] == "Total"]
    if len(totals) != 1:
        raise ValueError("Bench must provide exactly one Total summary")
    timing = positive_number(totals[0]["mean_seconds_per_sample"], "Bench Total")
    if not math.isclose(timing, statistics.mean(batch_totals), rel_tol=1e-12):
        raise ValueError("Bench Total disagrees with its per-sample batch timings")
    if fixture.representation is not None:
        representation = fixture.representation
        expected_attempt = [
            {
                "precision": "Double",
                "three_dimensional_representation": representation,
                "status": {"Stable": 2},
            }
        ]
        sample_counts = []
        for batch in batches:
            # Timed rows use the same fixed point and displayed normalization as
            # the warmup. Only this new representation gate enforces the bounds.
            bounds = batch.get("displayed_result_bounds")
            invalid = batch.get("invalid_result_samples")
            if (
                type(invalid) is not int
                or invalid != 0
                or not isinstance(bounds, list)
                or len(bounds) != 2
                or any(
                    not isinstance(component, list)
                    or len(component) != 2
                    or any(
                        type(endpoint) not in (int, float)
                        or not math.isfinite(endpoint)
                        for endpoint in component
                    )
                    or component[0] > component[1]
                    for component in bounds
                )
            ):
                raise ValueError(
                    "Timed displayed results require finite component bounds"
                )
            tolerance = fixture.config["correctness_tolerance"]
            if any(
                not math.isclose(
                    endpoint,
                    reference,
                    rel_tol=tolerance["relative"],
                    abs_tol=tolerance["absolute"],
                )
                for component, reference in zip(bounds, numbers)
                for endpoint in component
            ):
                raise ValueError("Timed displayed results disagree with the warmup")
            signatures = batch["attempt_signatures"]
            if not signatures or any(
                set(signature) != {"samples", "attempts"}
                or type(signature["samples"]) is not int
                or signature["samples"] <= 0
                or signature["attempts"] != expected_attempt
                for signature in signatures
            ):
                raise ValueError(
                    "Every timed sample must use one stable Double representation"
                )
            sample_counts.append(sum(signature["samples"] for signature in signatures))
        base, remainder = divmod(value["total_samples"], value["n_batches"])
        if sample_counts != [
            base + (index < remainder) for index in range(value["n_batches"])
        ]:
            raise ValueError(
                "Timed attempt identities do not match each native batch size"
            )
        info = value["integrand_info"]
        if (
            info["process_id"] != fixture.corpus["process_id"]
            or info["integrand_name"] != fixture.corpus["integrand"]
            or info["generated_representations"] != [representation]
            or info["kind"] != "cross_section"
            or info["generation_compilation"] != "Eager"
            or info["active_f64_backend"] != "eager"
            or info["graph_count"] != 1
            or info["graph_group_count"] != 1
            or len(info["graph_groups"]) != 1
        ):
            raise ValueError("Unexpected generated representation or evaluator backend")
        group = info["graph_groups"][0]
        if not group["complete_residue_sum"] or len(group["graphs"]) != 1:
            raise ValueError("The benchmark requires one complete residue sum")
        graph = group["graphs"][0]
        rows = graph["native_residue_counts"]
        if (
            graph["name"] != fixture.corpus["graph_name"]
            or len(rows) != 1
            or rows[0][0] != representation
            or type(rows[0][1]) is not int
            or rows[0][1] <= 0
        ):
            raise ValueError("Expected a nontrivial native residue-map catalogue")
        # Compare physical edge identities, not representation-local surface IDs
        # or source residue counts (which are intentionally different).
        cuts = []
        viable = 0
        for cut in group["cuts"]:
            physical_cut = {
                "edge_ids": sorted(cut["edge_ids"]),
                "raising_power": cut["raising_power"],
            }
            for side in ("left_thresholds", "right_thresholds"):
                thresholds = []
                for threshold in cut[side]:
                    thresholds.append(
                        {
                            "cut_boundary_edge_ids": sorted(
                                threshold["cut_boundary_edge_ids"]
                            ),
                            "threshold_boundary_edge_ids": sorted(
                                threshold["threshold_boundary_edge_ids"]
                            ),
                            "status": threshold["status"],
                            "invariant_bound_is_applicable": threshold[
                                "invariant_bound_is_applicable"
                            ],
                        }
                    )
                    viable += threshold["status"] in (
                        "can_become_pinched",
                        "potentially_existing",
                    )
                physical_cut[side] = sorted(
                    thresholds, key=lambda item: json.dumps(item, sort_keys=True)
                )
            cuts.append(physical_cut)
        if len(cuts) < 2 or not viable:
            raise ValueError(
                "The LTD/CFF corpus must retain cuts and viable thresholds"
            )
        physical_inventory = {
            "graph": graph["name"],
            "cuts": sorted(cuts, key=lambda item: json.dumps(item, sort_keys=True)),
            "loop_momentum_bases": sorted(
                sorted(basis["edge_ids"]) for basis in group["loop_momentum_bases"]
            ),
            "threshold_edges": sorted(
                sorted(threshold["edge_ids"])
                for threshold in group["threshold_esurfaces"]
            ),
            "threshold_directives": graph.get("threshold_counterterm_directives", []),
        }
        if fixture.certified_imaginary_zero is not None:
            # The numerator cancels every threshold-bearing cut. Never infer
            # this support certificate from small numerical components.
            expected = {
                (1, 2): (
                    (),
                    (
                        ((1, 2), (4, 5, 6), "can_become_pinched"),
                        ((1, 2), (4, 8), "can_become_pinched"),
                    ),
                ),
                (1, 7, 8): ((), (((8,), (5, 6), "potentially_existing"),)),
                (3, 4): ((), ()),
            }
            if (
                len(cuts) != 3
                or {tuple(cut["edge_ids"]) for cut in cuts} != expected.keys()
                or physical_inventory["threshold_directives"]
                or any(
                    cut["raising_power"] != 1
                    or tuple(
                        tuple(
                            sorted(
                                (
                                    tuple(threshold["cut_boundary_edge_ids"]),
                                    tuple(threshold["threshold_boundary_edge_ids"]),
                                    threshold["status"],
                                )
                                for threshold in cut[side]
                            )
                        )
                        for side in ("left_thresholds", "right_thresholds")
                    )
                    != expected[tuple(cut["edge_ids"])]
                    or any(
                        not threshold["invariant_bound_is_applicable"]
                        for side in ("left_thresholds", "right_thresholds")
                        for threshold in cut[side]
                    )
                    for cut in cuts
                )
            ):
                raise ValueError(
                    "GL04 source-zero certificate requires the original simple cuts and threshold support"
                )
        if metadata is not None:
            metadata.append(
                {
                    "representation": representation,
                    "timed_attempt_samples_per_batch": sample_counts,
                    "native_residue_count": rows[0][1],
                    "physical_inventory": physical_inventory,
                }
            )
    return timing, numbers, value["total_samples"]


def generate(binary, repository, directory, fixture, deadline):
    directory.mkdir(parents=True, exist_ok=False)
    state = directory / "state"
    generation_card = directory / "generation.toml"
    generation_card.write_text(fixture.card)
    from gammabench import hyperfine

    binding = hyperfine.resolve()
    generation_command = [
        str(binary),
        "--clean-state",
        "-s",
        str(state),
        str(generation_card),
    ]
    export = directory / "generation-hyperfine.json"
    run_cli(
        hyperfine.invocation(binding, generation_command, export),
        repository,
        directory / "generation.log",
        deadline,
        logfile_filter=(
            "[{#generation,#profile,#compile,#summary}]=debug"
            if "evaluator_build_seconds" in fixture.config["required_metrics"]
            else "off"
        ),
    )
    measured = hyperfine.read_measurement(export, generation_command)
    if measured["returncode"]:
        raise ValueError("Generation CLI failed")
    generation_seconds = measured["process_wall_ns"] / 1e9
    if fixture.certified_imaginary_zero is not None:
        # Check implicit native defaults after the generation timer. The frozen
        # card/model hashes alone cannot constrain a changed native default.
        generation = tomllib.loads((state / "global_settings.toml").read_text())[
            "global"
        ]["generation"]
        runtime = tomllib.loads((state / "default_runtime_settings.toml").read_text())
        overrides = tomllib.loads(
            (
                state
                / "processes/cross_sections/scalar_gl04/numerator/integrand/settings.toml"
            ).read_text()
        )
        general = runtime["general"] | overrides.get("general", {})
        uv = generation["uv"]
        model = json.loads((state / "model.json").read_text())
        parameters = {parameter["name"]: parameter for parameter in model["parameters"]}
        external = json.loads((state / "model_parameters.json").read_text())
        phase_values = {
            "ZERO": [0.0, 0.0],
            "lam": [1.0, 0.0],
            "mass_scalar_1": [0.1, 0.0],
            "width_scalar_1": [0.0, 0.0],
            "mass_scalar_2": [2.0, 0.0],
            "width_scalar_2": [0.0, 0.0],
        }
        if (
            any(
                parameters[name]["value"] != value
                or parameters[name]["parameter_type"] != "real"
                for name, value in phase_values.items()
            )
            or any(
                external.get(name, value) != value
                for name, value in phase_values.items()
            )
            or not {"lam", "mass_scalar_1", "mass_scalar_2"} <= external.keys()
            or [
                (coupling["name"], coupling["value"]) for coupling in model["couplings"]
            ]
            != [("SCALAR_COUPLING", [0.0, -1.0])]
            or any(
                uv[key] is not True
                for key in (
                    "softct",
                    "subtract_uv",
                    "generate_integrated",
                    "local_uv_cts_from_expanded_4d_integrands",
                )
            )
            or uv["renormalization_prescription"]
            != {
                "log_divergent": "MUV",
                "massive_power_divergent": "MUV",
                "massless_power_divergent": "MUV",
                "overrides": [],
            }
            or uv["vakint"]["normalization"] != "MSbar"
            or uv["vakint"]["additional_normalization"] != "1"
            or uv["vakint"]["evaluation_methods"] != ["alphaloop", "matad", "fmft"]
            or generation["threshold_subtraction"]["enable_thresholds"] is not True
            or any(
                general[key] != value
                for key, value in {
                    "m_uv": 20.0,
                    "mu_r": 3.0,
                    "renormalization_localization_scale": 1000.0,
                    "numerator_sampling_scale": 1.0,
                    "additional_param_values": [],
                    "disable_flux_factor": False,
                    "integral_unit": "auto",
                }.items()
            )
            or generation.get("medium", {}).get("mode", "vacuum") != "vacuum"
            or generation.get("medium", {}).get("vacuum_subtraction", False)
            or runtime["model"]
            or overrides.get("model")
            or any(
                settings.get("kinematics", {}).get("e_cm", 1.0) != 1.0
                for settings in (runtime, overrides)
            )
            or any(
                settings["kinematics"]["externals"]["type"] != "constant"
                or settings["kinematics"]["externals"]["data"]["momenta"]
                != [[1.0, 0.0, 0.0, 0.0]]
                or settings["kinematics"]["externals"]["data"]["helicities"] != [0]
                for settings in (runtime, overrides)
            )
            or runtime["subtraction"]["disable_threshold_subtraction"] is not False
            or overrides.get("subtraction", {}).get(
                "disable_threshold_subtraction", False
            )
            is not False
        ):
            raise ValueError(
                "GL04 source-zero certificate requires the original real model, MUV/MSbar and physical settings"
            )
    summaries = list(state.rglob("generation_summary.json"))
    if len(summaries) != 1:
        raise ValueError("Fresh generation must save one integrand summary")
    summary = json.loads(summaries[0].read_text())
    reports = summary["reports"]
    if len(reports) != 1 or reports[0]["graph_name"] != fixture.corpus["graph_name"]:
        raise ValueError("Fresh generation must build the requested graph")
    stats = reports[0]["stats"]
    count = stats["evaluator_count"]
    if type(count) is not int or count <= 0:
        raise ValueError("Fresh generation did not build any evaluators")
    peak_ram = summary["peak_ram_bytes"]
    if type(peak_ram) is not int or not 0 < peak_ram < 2**64:
        raise ValueError("Invalid whole-generation peak RAM")
    metrics = {
        "generation_seconds": generation_seconds,
        "generation_peak_ram_bytes": peak_ram,
    }
    construction = {}
    if "evaluator_build_seconds" in fixture.config["required_metrics"]:
        build_metrics, construction = fixture.construction_metrics(state, deadline)
        metrics.update(build_metrics)
    generation_timings = {}
    for stage in (
        "total_time",
        "evaluator_spenso_time",
        "evaluator_symbolica_time",
        "evaluator_compile_time",
    ):
        duration = stats[stage]
        seconds, nanos = duration["secs"], duration["nanos"]
        if (
            type(seconds) is not int
            or not 0 <= seconds < 2**64
            or type(nanos) is not int
            or not 0 <= nanos < 10**9
        ):
            raise ValueError(f"Invalid native generation duration: {stage}")
        generation_timings[f"generation_{stage.removesuffix('_time')}_seconds"] = (
            seconds + nanos / 10**9
        )
    return (
        state,
        generation_seconds,
        {
            "evaluator_count": count,
            "generation_peak_ram_bytes": peak_ram,
            **generation_timings,
            **construction,
            "metrics": metrics,
        },
    )


def precise_result(path, fixture, commit):
    from gammabench.gate_policy import PairedGate

    value = json.loads(path.read_text())
    expected = {
        "commit": commit,
        "process_id": fixture.corpus["process_id"],
        "integrand": fixture.corpus["integrand"],
        "precision_bits": 1000,
        "canonical_point_bits": fixture.ray["canonical_point_bits"],
        "warmup_samples": 1,
        "measured_samples": 1,
    }
    if any(value[key] != expected_value for key, expected_value in expected.items()):
        raise ValueError(
            "Precise producer changed revision, sample, routing or precision"
        )
    Fixture.paper_real_zero_bound(value["certified_real_zero"])
    permutation = value["canonical_permutation"]
    if (
        not isinstance(permutation, list)
        or any(type(i) is not int for i in permutation)
        or sorted(permutation) != [0, 1, 2]
    ):
        raise ValueError("Precise producer did not certify a loop-momentum permutation")
    if value["point_bits"] != [
        fixture.ray["canonical_point_bits"][i] for i in permutation
    ]:
        raise ValueError("Precise producer changed the mapped physical point")
    positive_number(value["load_seconds"], "Precise state load")
    positive_number(value["model_warmup_seconds"], "Resolved-model warmup")
    evaluations = value["evaluations"]
    if not isinstance(evaluations, list) or len(evaluations) != 2:
        raise ValueError(
            "Precise producer must record one warmup and one measured sample"
        )
    for evaluation in evaluations:
        evaluator = positive_number(
            evaluation["evaluator_seconds"], "Native evaluator time"
        )
        if (
            positive_number(evaluation["physical_seconds"], "Physical sample time")
            < evaluator
        ):
            raise ValueError("Evaluator time exceeds the complete physical call")
        attempts = evaluation["attempts"]
        if (
            len(attempts) != 1
            or attempts[0]["precision"] != "Arb"
            or attempts[0]["status"] != {"Stable": 2}
        ):
            raise ValueError("Precise producer changed attempts or stability probes")
        if attempts[0].get("three_dimensional_representation", "cff") != "cff":
            raise ValueError("Precise producer did not evaluate CFF")
        output = evaluation["output"]
        if not isinstance(output, list) or len(output) != 2:
            raise ValueError(
                "Precise output must retain both finite decimal components"
            )
        numbers = [
            PairedGate.number(component, decimal_output=True) for component in output
        ]
        if not numbers[1]:
            raise ValueError(
                "Certified paper sample lost its nonzero imaginary component"
            )
    if evaluations[0]["output"] != evaluations[1]["output"]:
        raise ValueError("Warmup and measured precise values disagree")
    return evaluations[1]["evaluator_seconds"], evaluations[1]["output"], value


def measure(binary, repository, directory, fixture, deadline, prepared=None):
    if fixture.precise:
        directory.mkdir(parents=True, exist_ok=False)
        state, producer, commit = prepared
        result = directory / "sample.json"
        run_cli(
            [str(producer), str(state), str(fixture.certificate), str(result)],
            repository,
            directory / "sample.log",
            deadline,
        )
        timing, output, metadata = precise_result(result, fixture, commit)
        return (
            {
                "evaluator_seconds_per_sample": timing,
                "warm_runtime_seconds_per_sample.SingleParametric": metadata[
                    "evaluations"
                ][1]["physical_seconds"],
            },
            {"SingleParametric": output},
            metadata,
        )
    state, _generation_seconds, generation_metadata = generate(
        binary, repository, directory, fixture, deadline
    )
    benchmark_card = fixture.benchmark_card(directory)
    run_cli(
        [str(binary), "--read-only-state", "-n", "-s", str(state), str(benchmark_card)],
        repository,
        directory / "benchmark.log",
        deadline,
    )
    metrics = generation_metadata.pop("metrics")
    outputs, sample_counts, per_point = {}, {}, {}
    representation_evidence = []
    for method in fixture.methods:
        per_point[method] = []
        timings, outputs[method], sample_counts[method] = [], [], []
        for index in range(len(fixture.corpus["points"])):
            check_deadline(deadline)
            timing, output, samples = benchmark_result(
                directory / f"bench-{method}-{index}.json",
                fixture,
                index,
                metadata=representation_evidence,
            )
            per_point[method].append(timing)
            timings.append(timing)
            outputs[method].extend(output)
            sample_counts[method].append(samples)
        metrics[f"warm_runtime_seconds_per_sample.{method}"] = statistics.mean(timings)
    check_deadline(deadline)
    return (
        {metric: metrics[metric] for metric in fixture.metrics},
        outputs,
        {
            "samples_per_point": sample_counts,
            "runtime_seconds_per_point": per_point,
            **(
                {"representation_evidence": representation_evidence}
                if fixture.representation
                else {}
            ),
            **generation_metadata,
        },
    )


def run_case(args, suite, name):
    """Prepare this workload; gammabench owns all scheduling and gate policy."""
    from gammabench import hyperfine
    from gammabench.gate import run
    from gammabench.snapshot import sha256

    fixture = Fixture.from_suite(suite, name, args.input_root)
    if fixture.representation:
        # The historical executable has no LTD. Both roles deliberately use the
        # same candidate binary; only their explicit representation differs.
        args = argparse.Namespace(
            **(
                vars(args)
                | {
                    "baseline_binary": args.candidate_binary,
                    "baseline_commit": args.candidate_commit,
                    "build": getattr(args, "build", {})
                    | {"baseline": getattr(args, "build", {}).get("candidate", {})},
                }
            )
        )
    prepared = args.work_dir.parent / (args.work_dir.name + "-suite")
    prepared.mkdir(parents=True, exist_ok=False)
    # Retain the complete public fixture inputs with their repository-relative
    # layout. Both binaries consume these exact bytes; generated states are private.
    (prepared / "inputs").mkdir()
    inputs = {}
    for index, relative in enumerate(fixture.corpus["inputs"]):
        target = prepared / "inputs" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(fixture.source(args.input_root, relative), target)
        inputs[f"input-{index}"] = str(target.relative_to(prepared))
    card = fixture.shared_card or fixture.card
    (prepared / "card.toml").write_text(card)
    representation_hashes = (
        {
            representation: Fixture.from_suite(
                suite, name, args.input_root, representation
            ).hashes["fixture_sha256"]
            for representation in ("cff", "ltd")
        }
        if fixture.representation
        else {}
    )
    write_json(
        prepared / "fixture.json",
        {
            "card": card,
            "corpus": fixture.corpus,
            "config": fixture.config,
            "representation_hashes": representation_hashes,
        },
    )
    shutil.copyfile(__file__, prepared / "adapter.py")
    inputs.update(card="card.toml", fixture="fixture.json")
    tools = {
        side: {
            "kind": "executable",
            "path": str(getattr(args, side + "_binary")),
            "sha256": sha256(getattr(args, side + "_binary")),
            "version": getattr(args, side + "_commit"),
        }
        for side in ("baseline", "candidate")
    }
    binding = hyperfine.resolve()
    tools["hyperfine"] = {
        "kind": "executable",
        "path": binding["path"],
        "sha256": binding["sha256"],
        "version": binding["version"],
    }
    tools["python"] = {
        "kind": "executable",
        "path": sys.executable,
        "sha256": sha256(Path(sys.executable)),
        "version": sys.version.split()[0],
    }
    if fixture.precise:
        for side in ("baseline", "candidate"):
            producer = getattr(args, side + "_precise_producer", None)
            if producer is None or not producer.is_file():
                raise ValueError(
                    "The slow precise case requires both source-bound producers"
                )
            tools[side + "-producer"] = {
                "kind": "executable",
                "path": str(producer),
                "sha256": sha256(producer),
                "version": getattr(args, side + "_commit"),
            }
    ci = getattr(args, "ci", {})
    contract = "fixed-point-integrand-v1"
    metrics = {
        metric: {
            "unit": "bytes"
            if metric.endswith("bytes")
            else "s/sample"
            if metric.startswith("warm_runtime")
            else "s",
            "direction": "lower",
            "result_field": "metrics." + metric,
            "scope": "Hyperfine 2.0.0 complete generation CLI process including import and state save"
            if metric == "generation_seconds"
            else "Native whole-generation peak RAM"
            if metric.endswith("bytes")
            else "Native evaluator construction stages summed over stacks"
            if metric.startswith("evaluator_build")
            else "Native 1000-bit evaluator time; numerical warmup excluded"
            if metric == "evaluator_seconds_per_sample"
            else "Complete native 1000-bit physical call; model/numerical warmup excluded"
            if fixture.precise
            else "Mean native warmed seconds per sample over the fixed corpus; ten warmup samples excluded",
        }
        for metric in fixture.metrics
    }
    definition = {
        "schema_version": 1,
        "runner": "paired-gate-v1",
        "id": "gammaloop-ci-" + name,
        "description": "GammaLoop CI " + name,
        "contracts": {
            contract: {
                "description": "Finite nonzero integrands agree within the declared tolerance across revisions and requested evaluators",
                "included": [
                    "Fresh generated state",
                    "Declared evaluator implementations",
                    "Fixed-point numerical validation",
                ],
                "excluded": ["CLI compilation", "Point caching", "Event processing"],
            }
        },
        "metadata": {
            "schema_version": 1,
            "values": {
                "implementation.tool": "gammaloop",
                "implementation.tool_label": "GammaLoop",
                "implementation.way": "native",
                "implementation.way_label": "Native",
            },
        },
        "metrics": metrics,
        "tools": tools,
        "gate_adapter": "adapter.py",
        "ci": ci,
        "build": getattr(args, "build", {}),
        "gate": {
            side: {
                "commit": getattr(args, side + "_commit"),
                "fixture_sha256": fixture.hashes.get(
                    "shared_fixture_sha256", fixture.hashes["fixture_sha256"]
                ),
                "corpus_sha256": fixture.hashes["corpus_sha256"],
                "config_sha256": fixture.hashes["config_sha256"],
            }
            for side in ("baseline", "candidate")
        },
        "tasks": [
            {
                "id": side,
                "case": name,
                "engine": "gammaloop",
                "variant": side,
                "contract": contract,
                "tool_ids": [side, "hyperfine", "python"],
                "inputs": inputs,
                "process_measurements": {
                    "generation": {
                        "path": "trial/generation-hyperfine.json",
                        "metric": "generation_seconds",
                        "argv": [
                            "{tool:" + side + "}",
                            "--clean-state",
                            "-s",
                            "{workdir}/trial/state",
                            "{workdir}/trial/generation.toml",
                        ],
                    }
                },
                "environment": {
                    "RAYON_NUM_THREADS": "1",
                    "OMP_NUM_THREADS": "1",
                    "GAMMABENCH_PRIVATE_TRIAL": "1",
                },
                "argv": [
                    "{tool:python}",
                    "{suite_dir}/adapter.py",
                    "sample",
                    "--fixture",
                    "{suite_dir}/fixture.json",
                    "--binary",
                    "{tool:" + side + "}",
                    "--work-dir",
                    "{workdir}",
                    "--result",
                    "{result}",
                    "--deadline",
                    "{deadline}",
                    "--side",
                    side,
                    "--commit",
                    getattr(args, side + "_commit"),
                ],
            }
            for side in ("baseline", "candidate")
        ],
    }
    if fixture.precise:
        for task in definition["tasks"]:
            side = task["variant"]
            task.pop("process_measurements")
            task["tool_ids"].append(side + "-producer")
            task["argv"].extend(["--producer", "{tool:" + side + "-producer}"])
    policy = {
        "schema_version": 1,
        "pairs": 5,
        "budget_seconds": fixture.timeout,
        "thresholds": {
            metric: [
                0.5 if fixture.representation else THRESHOLDS[metric][0],
                THRESHOLDS[metric][1],
            ]
            for metric in fixture.metrics
        },
        "output_groups": fixture.methods,
        "reference_group": fixture.config["reference_method"],
        "metric_families": list(MODE_METRICS),
        "correctness_tolerance": fixture.config["correctness_tolerance"],
        **({"budget_class": "slow"} if fixture.precise else {}),
        **fixture.qualification_policy(),
    }
    write_json(prepared / "suite.json", definition)
    write_json(prepared / "policy.json", policy)
    capture = (
        "gammaloop-"
        + args.candidate_commit
        + "-"
        + str(ci.get("run_id", "local"))
        + "-"
        + str(ci.get("run_attempt", 1))
        + "-"
        + name
    )
    return run(
        prepared / "suite.json",
        args.output_dir,
        policy_path=prepared / "policy.json",
        private_work_dir=args.work_dir,
        label=capture,
    )


def sample_main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--deadline", required=True, type=float)
    parser.add_argument("--side", choices=("baseline", "candidate"), required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--producer", type=Path)
    args = parser.parse_args(sys.argv[2:])
    value = json.loads(args.fixture.read_text())
    card = value["card"]
    representation = None
    if value["config"].get("comparison") == "ltd_vs_cff":
        representation = "cff" if args.side == "baseline" else "ltd"
        card = Template(card).substitute(representation=representation)
    fixture = Fixture(
        card, value["corpus"], value["config"], args.fixture.parent / "inputs"
    )
    if representation:
        shared = Fixture(
            value["card"], value["corpus"], value["config"], fixture.input_root
        )
        fixture.hashes["shared_fixture_sha256"] = shared.hashes["fixture_sha256"]
        if (
            fixture.hashes["fixture_sha256"]
            != value["representation_hashes"][representation]
            or fixture.hashes["fixture_sha256"]
            != Fixture.GL04_REPRESENTATION_HASHES[representation]
            or any(
                fixture.hashes[k] != expected
                for k, expected in Fixture.GL04_SOURCE_HASHES.items()
            )
        ):
            raise ValueError("Prepared LTD/CFF physical fixture identity changed")
        fixture.representation = representation
        fixture.certified_imaginary_zero = Fixture.GL04_IMAGINARY_ZERO.copy()
    # Gammabench creates private/phase/repeat/role. This cache belongs only to
    # the current fresh capture and shares each expensive precise preparation
    # across its five native samples and any independent confirmation.
    private = args.work_dir.parents[2]
    prepared = None
    preparation = {}
    if fixture.precise:
        if args.producer is None:
            raise ValueError("Missing source-bound precise producer")
        state_root = private / "prepared" / args.side
        receipt = state_root / "preparation.json"
        if not receipt.exists():
            state, seconds, details = generate(
                args.binary, fixture.input_root, state_root, fixture, args.deadline
            )
            preparation = {
                "generation_seconds": seconds,
                "generation_is_paired": 0,
                "evaluator_count": details["evaluator_count"],
                "generation_peak_ram_bytes": details["generation_peak_ram_bytes"],
            }
            write_json(
                receipt,
                {
                    "commit": args.commit,
                    "fixture_hashes": fixture.hashes,
                    "state": str(state),
                    "metadata": preparation,
                },
            )
        saved = json.loads(receipt.read_text())
        if saved["commit"] != args.commit or saved["fixture_hashes"] != fixture.hashes:
            raise ValueError("Precise state preparation identity changed")
        prepared = (Path(saved["state"]), args.producer, args.commit)
    metrics, outputs, metadata = measure(
        args.binary,
        fixture.input_root,
        args.work_dir / "trial",
        fixture,
        args.deadline,
        prepared,
    )
    if fixture.precise:
        # Raw native decimal output and source proof remain in trial/sample.json;
        # the public record retains exact decimal outputs and its certificate.
        metadata = {
            "preparation": preparation,
            "precision_bits": 1000,
            "load_seconds": metadata["load_seconds"],
            "model_warmup_seconds": metadata["model_warmup_seconds"],
            "canonical_permutation": metadata["canonical_permutation"],
            "point_bits": metadata["point_bits"],
        }
    if representation:
        evidence = metadata.pop("representation_evidence")
        if len(evidence) != len(fixture.corpus["points"]):
            raise ValueError("Missing point representation evidence")
        inventory = evidence[0]["physical_inventory"]
        if any(row["physical_inventory"] != inventory for row in evidence):
            raise ValueError(
                "Physical cut and threshold inventories differ across points"
            )
        saved_inventory = private / "physical-inventory.json"
        if saved_inventory.exists():
            if json.loads(saved_inventory.read_text()) != inventory:
                raise ValueError(
                    "LTD/CFF physical cut and threshold inventories differ"
                )
        else:
            write_json(saved_inventory, inventory)
        digest = hashlib.sha256(
            json.dumps(inventory, sort_keys=True).encode()
        ).hexdigest()
        metadata["physical_inventory_sha256_words"] = [
            int(digest[i : i + 8], 16) for i in range(0, 64, 8)
        ]
        metadata["native_residue_counts"] = [
            row["native_residue_count"] for row in evidence
        ]
        metadata["representation_is_ltd"] = int(representation == "ltd")
    result = {"metrics": metrics, "outputs": outputs, "metadata": metadata}
    qualification = fixture.qualification_policy()
    if qualification:
        result["certificate"] = qualification["certified_zero"]["certificate"]
    write_json(args.result, result)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for side in ("baseline", "candidate"):
        parser.add_argument(f"--{side}-binary", type=Path, required=True)
        parser.add_argument(f"--{side}-commit", required=True)
        parser.add_argument(f"--{side}-precise-producer", type=Path)
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Shared candidate repository supplying both revisions' inputs",
    )
    parser.add_argument(
        "--case", help="Run one manifest case, including the explicit slow precise gate"
    )
    parser.add_argument(
        "--lane",
        choices=("fast", "slow", "all"),
        default="fast",
        help="Default: thirteen revision comparisons and one LTD/CFF comparison; slow generation needs a >=30GB host",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="Private state/log directory, never uploaded",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Metadata and numeric report directory",
    )
    parser.add_argument("--ci-metadata", type=Path)
    parser.add_argument("--build-metadata", type=Path)
    args = parser.parse_args()
    args.ci = json.loads(args.ci_metadata.read_text()) if args.ci_metadata else {}
    args.build = (
        json.loads(args.build_metadata.read_text()) if args.build_metadata else {}
    )
    try:
        for side in ("baseline", "candidate"):
            if not re.fullmatch(r"[0-9a-f]{40}", getattr(args, side + "_commit")):
                raise ValueError("Commit identities must be complete hexadecimal SHA1s")
            binary = getattr(args, side + "_binary").resolve()
            if not binary.is_file() or not os.access(binary, os.X_OK):
                raise ValueError("Missing executable")
            setattr(args, side + "_binary", binary)
            producer = getattr(args, side + "_precise_producer")
            if producer is not None:
                producer = producer.resolve()
                if not producer.is_file() or not os.access(producer, os.X_OK):
                    raise ValueError("Missing precise performance producer")
                setattr(args, side + "_precise_producer", producer)
        args.input_root = args.input_root.resolve()
        if not args.input_root.is_dir():
            raise ValueError("Missing shared input repository")
        args.work_dir = args.work_dir.resolve()
        args.output_dir = args.output_dir.resolve()
        if args.work_dir.exists():
            raise ValueError(
                "Private work directory must be new; no generated state may be reused"
            )
        if args.work_dir.is_relative_to(
            args.output_dir
        ) or args.output_dir.is_relative_to(args.work_dir):
            raise ValueError(
                "Private state/logs must be separate from uploaded metadata"
            )
        if args.output_dir.exists():
            raise ValueError(
                "Uploaded metadata directory must be new; stale reports cannot be reused"
            )
        suite = tomllib.loads(
            (args.input_root / "benchmarks/performance/suite.toml").read_text()
        )
        cases = Fixture.cases(suite)
        names = (
            [args.case]
            if args.case
            else [
                name
                for name, (workload, _) in cases.items()
                if args.lane == "all"
                or (workload["family"] == "precise") == (args.lane == "slow")
            ]
        )
        if any(name not in cases for name in names):
            raise ValueError("Unknown performance case")
        args.work_dir.mkdir(parents=True, exist_ok=False)
        args.output_dir.mkdir(parents=True, exist_ok=False)
        status = 0
        captured_at = datetime.now(UTC).isoformat()
        for name in names:
            print(f"Measuring {name}", flush=True)
            case_args = argparse.Namespace(
                **(
                    vars(args)
                    | {
                        "work_dir": args.work_dir / name,
                        "output_dir": args.output_dir / name,
                    }
                )
            )
            try:
                case_status = run_case(case_args, suite, name)
            except (OSError, ValueError, KeyError, TypeError):
                case_status = 2
                case_args.output_dir.mkdir(parents=True, exist_ok=True)
                write_json(
                    case_args.output_dir / "failure.json",
                    {"schema_version": 1, "case": name, "status": 2},
                )
            status = max(status, case_status)
        write_json(
            args.output_dir / "run.json",
            {
                "schema_version": 1,
                "ci": args.ci,
                "cases": names,
                "status": status,
                "created_at": captured_at,
            },
        )
        return status
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Performance suite failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(sample_main() if sys.argv[1:2] == ["sample"] else main())
