"""Reproducible numerator benchmarks for fermionic and gluonic propagator ladders.

Run with the combined Symbolica/Spenso host, not the stock Symbolica wheel.
The notebook imports this driver so its timed route and CLI route stay identical.
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import math
import os
import platform
import re
import shutil
import subprocess
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from time import perf_counter_ns, process_time_ns

from symbolica import E, Replacement, S
from symbolica.community.spenso import (
    AUTO,
    GammaSimplifySettings,
    Representation,
    TensorExpression,
    TensorName,
    trace,
)


class FermionLadder:
    """Massless fermion ring with planar vector rungs and a metric projection.

    Both modes use tr(1)=4; only the Lorentz dimension changes. Momenta run
    clockwise along the fermion ring. The external momentum q enters vertex1
    and exits vertex loops+1. No denominators, couplings or color are included.
    """

    particle = "fermionic"
    dimension_symbol = S("fermion_ladder::D")

    def __init__(self, mode: str, loops: int = 3):
        if mode not in ("trace4", "tracen") or loops not in (3, 4):
            raise ValueError("Use trace4/tracen and three/four loops")
        self.mode, self.loops = mode, loops
        self.dimension = 4 if mode == "trace4" else self.dimension_symbol
        self.lorentz = Representation.mink(self.dimension)
        self.spin = Representation.bis(4)
        self.metric = TensorName.g().to_expression()
        self.momenta = [
            TensorName.vector(f"fermion_ladder::p{i}").to_expression()(
                self.lorentz.to_expression()
            )
            for i in range(1, 2 * loops + 1)
        ]
        upper = [tuple(int(i <= j) for i in range(loops)) + (0,) for j in range(loops)]
        self.routing = upper + [row[:-1] + (-1,) for row in reversed(upper)]
        gamma = TensorExpression.gamma(self.dimension)
        gamma_head = TensorName.gamma().to_expression()
        names = ["mu"] + [f"a{i}" for i in range(1, loops)]
        opposite = ["nu"] + [f"b{i}" for i in reversed(range(1, loops))]
        labels = {
            name: self.lorentz(S(f"fermion_ladder::{name}"))
            for name in names + opposite
        }

        def word(indices):
            factors = []
            for index, momentum in zip(indices, self.momenta, strict=True):
                factors.extend(
                    (
                        gamma(AUTO, AUTO, labels[index]),
                        TensorExpression(
                            gamma_head(
                                self.spin.to_expression(),
                                self.spin.to_expression(),
                                momentum,
                            )
                        ),
                    )
                )
            return trace(self.spin, *factors)

        # Sewing the metric pairs is input preparation, in both engines.
        self.source = -word(names + ["mu"] + list(reversed(names[1:])))
        projected = -word(names + opposite).to_expression()
        for left, right in [("mu", "nu")] + [
            (f"a{i}", f"b{i}") for i in range(1, loops)
        ]:
            projected *= self.metric(
                labels[left].to_expression(), labels[right].to_expression()
            )
        self.explicit = TensorExpression(projected)
        self.scalars = {
            (i, j): S(f"fermion_ladder::s{i + 1}{j + 1}")
            for i in range(loops + 1)
            for j in range(i, loops + 1)
        }
        # The three-loop 4D scalar is small enough to expand in one step. Larger
        # cases benefit from collecting each momentum substitution separately.
        self.settings = GammaSimplifySettings(
            expand_traces=mode == "tracen" or loops == 4
        )
        self.staged_routing = mode == "tracen" or loops == 4
        if self.staged_routing:
            self._prepare_polynomial_routing()
        else:
            self.replacements = []
            for i, left in enumerate(self.routing):
                for j in range(i, len(self.routing)):
                    right = self.routing[j]
                    routed = sum(
                        (
                            x * y * self.scalars[min(a, b), max(a, b)]
                            for a, x in enumerate(left)
                            for b, y in enumerate(right)
                            if x and y
                        ),
                        E("0"),
                    )
                    self.replacements.append(
                        Replacement(
                            self.metric(self.momenta[i], self.momenta[j]), routed
                        )
                    )
        self.form_order = "trace-first"
        self.strategy = {
            "expand_traces": self.settings.expand_traces,
            "scalar_expansion": (
                "grouped polynomial substitution, outer momentum pairs first"
                if self.staged_routing
                else "simultaneous dot substitution + expand(via_poly=True)"
            ),
            "form_order": self.form_order,
            "form_routing_collection": "after each momentum alias",
        }

    def _prepare_polynomial_routing(self):
        """Prepare a bilinear change of momentum basis, collected at each step."""
        count = self.loops + 1
        basis = [
            TensorName.vector(f"fermion_ladder::{name}").to_expression()(
                self.lorentz.to_expression()
            )
            for name in [f"k{i}" for i in range(1, count)] + ["q"]
        ]
        vectors = basis + self.momenta

        def dot(i, j):
            i, j = sorted((i, j))
            return (
                self.scalars[i, j] if j < count else self.metric(vectors[i], vectors[j])
            )

        self.polynomial_variables = (
            [self.dimension_symbol]
            + list(self.scalars.values())
            + [
                dot(i, j)
                for i in range(len(vectors))
                for j in range(max(i, count), len(vectors))
            ]
        )
        aliases = list(range(count, len(vectors)))
        order = [
            i
            for pair in zip(aliases[: self.loops], reversed(aliases[self.loops :]))
            for i in pair
        ]
        active = list(range(len(vectors)))
        self.routing_groups = []
        for i in order:
            row = self.routing[i - count]
            replacements = []
            for j in active:
                if i == j:
                    terms = (
                        a * b * dot(k, l)
                        for k, a in enumerate(row)
                        for l, b in enumerate(row)
                        if a and b
                    )
                else:
                    terms = (a * dot(k, j) for k, a in enumerate(row) if a)
                replacements.append(
                    sum(terms, E("0")).to_polynomial(vars=self.polynomial_variables)
                )
            self.routing_groups.append(([dot(i, j) for j in active], replacements))
            active.remove(i)

    def _substitute_momenta(self, polynomial):
        # All variables touching one momentum are substituted simultaneously.
        # Grouping their coefficients avoids rebuilding the same RHS for each
        # monomial. Collect each momentum step and emit an Atom only at the end.
        for variables, replacements in self.routing_groups:
            parts = []
            for powers, coefficient in polynomial.coefficient_list(variables):
                for exponent, replacement in zip(powers, replacements, strict=True):
                    if exponent:
                        coefficient *= replacement**exponent
                parts.append(coefficient)
            if not parts:
                return polynomial
            polynomial = sum(parts[1:], parts[0])
        return polynomial

    def reduce(self, source=None):
        """Trace, insert physical momentum routing, then fully collect the scalar."""
        if source is None:
            source = self.source
        traced = source.simplify_gamma(self.settings).to_expression()
        if self.staged_routing:
            polynomial = traced.to_polynomial(vars=self.polynomial_variables)
            return self._substitute_momenta(polynomial).to_expression()
        return traced.replace_multiple(self.replacements).expand(via_poly=True)

    def import_form(self, polynomial: str):
        positions = {f"k{i}": i for i in range(1, self.loops + 1)}
        positions["q"] = self.loops + 1

        def dot(match):
            i, j = sorted((positions[match[1]], positions[match[2]]))
            return f"fermion_ladder::s{i}{j}"

        text = re.sub(r"\b(k\d+|q)\.(k\d+|q)\b", dot, polynomial.strip().rstrip(";"))
        text = re.sub(r"\bD\b", "fermion_ladder::D", text)
        return E(text).expand()

    def phases(self, expected):
        """Separate diagnostic clock; never substituted for the complete timing."""
        start = process_time_ns()
        traced = self.source.simplify_gamma(self.settings).to_expression()
        after_trace = process_time_ns()
        if self.staged_routing:
            polynomial = traced.to_polynomial(vars=self.polynomial_variables)
            after_conversion = process_time_ns()
            polynomial = self._substitute_momenta(polynomial)
            after_routing = process_time_ns()
            result = polynomial.to_expression()
            end = process_time_ns()
            phases = {
                "polynomial_conversion_cpu_ns": after_conversion - after_trace,
                "polynomial_routing_cpu_ns": after_routing - after_conversion,
                "scalar_emission_cpu_ns": end - after_routing,
            }
        else:
            routed = traced.replace_multiple(self.replacements)
            after_routing = process_time_ns()
            result = routed.expand(via_poly=True)
            end = process_time_ns()
            phases = {
                "routing_cpu_ns": after_routing - after_trace,
                "expansion_cpu_ns": end - after_routing,
            }
        assert result == expected
        return {
            "trace_cpu_ns": after_trace - start,
            **phases,
            "expanded_trace_terms": len(list(traced.expand(via_poly=True).terms())),
            "final_terms": len(list(result.terms())),
        }

    def check_components(self, scalar_results, networks):
        from fermion_ladder_validation import FermionRingComponents

        return FermionRingComponents(self).check(scalar_results, networks=networks)


def _form_run(
    executable,
    script,
    directory,
    *,
    mode,
    loops,
    repeats=1,
    diagnostic=False,
    order=None,
):
    if order is None:
        order = "routing-first" if mode == "trace4" else "trace-first"
    output = Path(directory) / f"{mode}.txt"
    command = [
        executable,
        "-q",
        "-d",
        f"MODE={mode}",
        "-d",
        f"LOOPS={loops}",
        "-d",
        f"ORDER={order}",
        "-d",
        "COLLECT=1",
        "-d",
        "ROUTEGROUP=1",
        "-d",
        f"DIAGNOSTIC={int(diagnostic)}",
        "-d",
        f"REPEATS={repeats}",
        "-d",
        "BATCHES=1",
    ]
    if diagnostic:
        command += ["-d", f"OUTPUT={output}"]
    start = perf_counter_ns()
    run = subprocess.run(
        command + [str(script)],
        cwd=directory,
        text=True,
        capture_output=True,
        check=False,
        timeout=180,
    )
    wall = perf_counter_ns() - start
    if run.returncode:
        raise RuntimeError(f"FORM failed:\n{run.stdout}\n{run.stderr}")
    timers = re.findall(r"BATCH_CPU_MS=(\d+)", run.stdout)
    if not diagnostic and len(timers) != 1:
        raise RuntimeError(f"Expected one FORM body timer:\n{run.stdout}")
    return {
        "command": command + [str(script)],
        "source_sha256": hashlib.sha256(Path(script).read_bytes()).hexdigest(),
        "affinity": sorted(os.sched_getaffinity(0)),
        "term_counts": [int(value) for value in re.findall(r"TERMS=(\d+)", run.stdout)],
        "stdout": run.stdout,
        "stderr": run.stderr,
        "process_wall_ns": wall,
        "body_cpu_ns": int(timers[0]) * 1_000_000 if timers else None,
        "repeats": repeats,
        "polynomial": output.read_text() if diagnostic else None,
    }


def benchmark(
    form=None,
    *,
    loops=3,
    rounds=5,
    output=None,
    cpu=None,
    case_factory=FermionLadder,
    form_script=None,
):
    """Compare warm algebra bodies, excluding setup, checks and serialization."""

    def time_batch(case, repeats, expected):
        wall, cpu_start = perf_counter_ns(), process_time_ns()
        results = [case.reduce() for _ in range(repeats)]
        elapsed_cpu = process_time_ns() - cpu_start
        elapsed_wall = perf_counter_ns() - wall
        assert all(result == expected for result in results)
        return {
            "repeats": repeats,
            "body_cpu_ns": elapsed_cpu,
            "body_wall_ns": elapsed_wall,
        }

    if rounds < 3:
        raise ValueError("Use at least three measured rounds")
    if cpu is not None:
        os.sched_setaffinity(0, {cpu})
    form = form or os.environ.get("FORM_EXECUTABLE") or shutil.which("form")
    if not form:
        raise RuntimeError(
            "Set FORM_EXECUTABLE or pass --form to a native FORM executable"
        )
    form = str(Path(form).resolve())
    script = (
        Path(form_script)
        if form_script
        else Path(__file__).with_name("fermion_propagator_ladder.frm")
    )
    from symbolica import core

    report = {
        "measured_at_utc": datetime.now(UTC).isoformat(),
        "loops": loops,
        "particle": case_factory.particle,
        "internal_edges": 3 * loops - 1,
        "vertices": 2 * loops,
        "trace_unit": 4 if case_factory.particle == "fermionic" else None,
        "closed_fermion_loop_minus": case_factory.particle == "fermionic",
        "scope": "Numerator only; no denominators, couplings, color, integration or on-shell constraints",
        "timing": "Warm complete numerator reduction + momentum routing + expanded scalar collection; setup, checks, output and disposal excluded. Python call and result wrapping included. FORM internal CPU timer; process wall reported separately.",
        "rounds": rounds,
        "platform": platform.platform(),
        "cpu_affinity": sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else None,
        "core_sha256": hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),
        "driver_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "case_source_sha256": hashlib.sha256(
            Path(inspect.getfile(case_factory)).read_bytes()
        ).hexdigest(),
        "form_source_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "form_executable_sha256": hashlib.sha256(Path(form).read_bytes()).hexdigest(),
        "form_version": subprocess.run(
            [form, "-v"], text=True, capture_output=True, check=True
        ).stdout.strip(),
        "cases": [],
    }
    scalar_results, original_networks = {}, {}
    with tempfile.TemporaryDirectory(prefix="fermion-ladder-") as directory:
        for mode in ("trace4", "tracen"):
            setup_wall, setup_cpu = perf_counter_ns(), process_time_ns()
            case = case_factory(mode, loops)
            setup = {
                "wall_ns": perf_counter_ns() - setup_wall,
                "cpu_ns": process_time_ns() - setup_cpu,
            }
            start_wall, start_cpu = perf_counter_ns(), process_time_ns()
            expected = case.reduce()
            first = {
                "wall_ns": perf_counter_ns() - start_wall,
                "cpu_ns": process_time_ns() - start_cpu,
            }
            if case.explicit is not None:
                assert expected == case.reduce(case.explicit), (
                    "Explicit metric projection differs"
                )
                assert case.explicit.structure.rank == 0
            assert case.source.structure.rank == 0
            diagnostic = _form_run(
                form,
                script,
                directory,
                mode=mode,
                loops=loops,
                diagnostic=True,
                order=case.form_order,
            )
            reference = case.import_form(diagnostic["polynomial"])
            assert expected == reference, f"{mode}: exact FORM polynomial differs"
            diagnostic["polynomial_sha256"] = hashlib.sha256(
                diagnostic.pop("polynomial").encode()
            ).hexdigest()
            scalar_results[f"Idenso {mode}"] = expected
            scalar_results[f"FORM {mode}"] = reference
            original_networks[f"original compact {mode}"] = case.source
            if case.explicit is not None:
                original_networks[f"original explicit {mode}"] = case.explicit
            pilot_repeats = max(
                1, min(64, math.ceil(150_000_000 / max(first["cpu_ns"], 1)))
            )
            pilot = _form_run(
                form,
                script,
                directory,
                mode=mode,
                loops=loops,
                repeats=pilot_repeats,
                order=case.form_order,
            )
            form_repeats = min(
                8192,
                max(
                    1,
                    math.ceil(
                        200_000_000 / max(pilot["body_cpu_ns"], 1) * pilot_repeats
                    ),
                ),
            )
            idenso_repeats = max(
                1, min(512, math.ceil(150_000_000 / max(first["cpu_ns"], 1)))
            )
            samples = []
            for round_index in range(rounds):
                row = {"round": round_index + 1}
                order = (
                    ("idenso", "form") if round_index % 2 == 0 else ("form", "idenso")
                )
                for engine in order:
                    row[engine] = (
                        time_batch(case, idenso_repeats, expected)
                        if engine == "idenso"
                        else _form_run(
                            form,
                            script,
                            directory,
                            mode=mode,
                            loops=loops,
                            repeats=form_repeats,
                            order=case.form_order,
                        )
                    )
                if row["form"]["body_cpu_ns"] == 0:
                    raise RuntimeError("FORM batch is below the CPU timer resolution")
                samples.append(row)
            idenso_cpu = median(
                row["idenso"]["body_cpu_ns"] / idenso_repeats for row in samples
            )
            form_cpu = median(
                row["form"]["body_cpu_ns"] / form_repeats for row in samples
            )
            record = {
                "mode": mode,
                "dimension": "4D" if mode == "trace4" else "D",
                "routing": case.routing,
                "terms": len(list(expected.terms())),
                "strategy": case.strategy,
                "first_idenso": first,
                "input_and_rule_preparation": setup,
                "form_pilot": pilot,
                "samples": samples,
                "idenso_median_cpu_ms": idenso_cpu / 1e6,
                "idenso_median_wall_ms": median(
                    row["idenso"]["body_wall_ns"] / idenso_repeats for row in samples
                )
                / 1e6,
                "form_median_cpu_ms": form_cpu / 1e6,
                "cpu_ratio": idenso_cpu / form_cpu,
                "diagnostic": case.phases(expected),
                "exact_form_match": True,
                "explicit_projection_match": True
                if case.explicit is not None
                else None,
                "scalar_polynomial_sha256": hashlib.sha256(
                    expected.format_plain().encode()
                ).hexdigest(),
                "form_diagnostic": diagnostic,
            }
            report["cases"].append(record)
    assert (
        scalar_results["Idenso tracen"].replace(case.dimension_symbol, E("4")).expand()
        == scalar_results["Idenso trace4"]
    )
    report["D_to_four_exact"] = True
    report["component_checks"] = case.check_components(
        scalar_results, original_networks
    )
    for key, path in (
        ("driver_sha256", Path(__file__)),
        ("case_source_sha256", Path(inspect.getfile(case_factory))),
        ("form_source_sha256", script),
        ("core_sha256", Path(core.__file__)),
    ):
        assert hashlib.sha256(path.read_bytes()).hexdigest() == report[key], (
            f"Source changed during measurement: {path}"
        )
    report["within_factor_three"] = all(
        case["cpu_ratio"] <= 3 for case in report["cases"]
    )
    if output:
        Path(output).write_text(json.dumps(report, indent=2) + "\n")
    return report


def consolidation_child(args):
    """One interpreter/case/route. Setup and diagnostics never enter primary clocks."""
    from time import thread_time_ns

    from symbolica import core
    from tensor_benchmark_cases import HistoricalLadder, atom, describe, make_case
    from tensor_benchmark_diagnostics import (
        FactorProbe,
        keep_rest,
        partial_parse,
        validation_cost,
    )

    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    core_path = Path(core.__file__)
    before = hashlib.sha256(core_path.read_bytes()).hexdigest()
    if args.expected_core and before != args.expected_core:
        raise ValueError(f"Core identity changed: {before}")
    setup = process_time_ns()
    case = make_case(args.case)
    route = args.route
    if route == "scalar" and not isinstance(case, HistoricalLadder):
        raise ValueError(
            "The pure Symbolica recipe is defined only for the historical ladder"
        )
    if route == "factorized":
        from tensor_benchmark_factorized import prepare, reduce

        if not (args.case.startswith(("historical-", "gluon-"))):
            raise ValueError("The factorized route requires a gluonic ladder fixture")
        prepared = prepare(case)
        call = lambda: reduce(case, prepared)
    else:
        call = case.scalar if route == "scalar" else case.reduce
    setup = process_time_ns() - setup
    warmup = call()
    warmup_identity = hashlib.sha256(atom(warmup).format_plain().encode()).hexdigest()
    del warmup
    cpu_start, thread_start, wall_start = (
        process_time_ns(),
        thread_time_ns(),
        perf_counter_ns(),
    )
    outputs = [call() for _ in range(args.calls)]
    wall_ns = perf_counter_ns() - wall_start
    thread_ns, cpu_ns = thread_time_ns() - thread_start, process_time_ns() - cpu_start
    result = outputs[-1]
    assert (
        hashlib.sha256(atom(result).format_plain().encode()).hexdigest()
        == warmup_identity
    )
    assert all(value == result for value in outputs)
    # Conversion of scalar-reference notation is an oracle adapter, not algebra timing.
    result = case.convert_scalar(result) if route == "scalar" else atom(result)
    result_path = args.output.with_suffix(".expression")
    result_path.write_text(result.format_plain())
    record = {
        "input": describe(case.source),
        "case": args.case,
        "route": route,
        "calls": args.calls,
        "fixed_warmups": 1,
        "wall_ns": wall_ns,
        "process_ns": cpu_ns,
        "thread_ns": thread_ns,
        "setup_process_ns": setup,
        "core": str(core_path),
        "core_sha256": before,
        "affinity": sorted(os.sched_getaffinity(0)),
        "strategy": case.strategy,
        "result": describe(result),
        "expression": str(result_path),
        "scope": "Complete case.reduce / scalar recipe; Python dispatch and retained output list included; setup, checks, counts and disposal excluded",
        "resource": resource_snapshot(args.cpu),
    }
    if route == "factorized":
        phases = {}
        assert reduce(case, prepared, phases) == result
        record["phases"] = phases
        record["strategy"] = {
            "vertex_order": list(case.order),
            "rule_application": "reusable TensorRule per local vertex; factorized output",
            "contraction": "existing Rust component reducer with merged states",
            "materialization": "explicit polynomial forward pass over aliases",
        }
    else:
        record["phases"] = case.phases(result)
    if args.check:
        if hasattr(case, "explicit") and case.explicit is not None:
            assert result == case.reduce(case.explicit)
            record["explicit_projection_exact"] = True
        if isinstance(case, HistoricalLadder):
            assert describe(result)["terms"] == 9652
            typed = case.typed()
            assert atom(typed.schoonschip(case.settings)) == result
            record["scalar_rank_and_fixedpoint"] = typed.is_scalar
            assert record["scalar_rank_and_fixedpoint"]
        if args.form and hasattr(case, "import_form"):
            with tempfile.TemporaryDirectory(prefix="r3-form-check-") as directory:
                script = case_form_script(case, Path(directory))
                form = _form_run(
                    args.form,
                    script,
                    directory,
                    mode=case.mode,
                    loops=case.loops,
                    diagnostic=True,
                    order=case.form_order,
                )
                form_expression = case.import_form(form["polynomial"])
                literal = result == form_expression
                record["form_source"] = script.read_text()
                record["form_exact"] = literal
                if not literal and hasattr(case, "check_form_components"):
                    record["form_component_checks"] = case.check_form_components(
                        result, form_expression
                    )
                elif not literal:
                    args.output.write_text(
                        json.dumps({**record, "failed_form": form}, indent=2) + "\n"
                    )
                    raise AssertionError(
                        f"{args.case}: FORM identity differs; no component oracle is installed"
                    )
                form_path = args.output.with_suffix(".form.txt")
                form_path.write_text(form.pop("polynomial"))
                form["polynomial_path"] = str(form_path)
                form["polynomial_sha256"] = hashlib.sha256(
                    form_path.read_bytes()
                ).hexdigest()
                record["form_diagnostic"] = form
        if hasattr(case, "dimension_symbol") and case.mode == "tracen":
            other = make_case(args.case.replace("tracen", "trace4"))
            lowered = result.replace(case.dimension_symbol, E("4")).expand()
            strict_four = other.reduce()
            record["D_to_four_exact"] = lowered == strict_four
            if not record["D_to_four_exact"] and hasattr(case, "check_form_components"):
                record["D_to_four_component_checks"] = case.check_form_components(
                    lowered, strict_four
                )
            else:
                assert record["D_to_four_exact"]
        if hasattr(case, "check_components"):
            # Existing independent HEP evaluator, with the same exact input/output.
            record["components"] = case.check_components(
                {"result": result}, {"input": case.source}
            )
        if args.diagnostic:
            if args.diagnostic == "keep-rest":
                record["diagnostic"] = keep_rest()
                failures = []
                for row in record["diagnostic"]:
                    required = row.get("required_milestone", "M3")
                    milestone_required = args.gate_stage == "M3" or (
                        args.gate_stage in ("baseline", "M0") and required == "M0"
                    )
                    known_baseline = row.get("spinor_slots") == "AUTO" or row.get(
                        "case"
                    ) in ("vector_power", "foreign_color")
                    semantic_ok = row.get(
                        "exact_expanded", True
                    ) is not False and row.get("rerun", True)
                    shape_ok = row.get("preserves_outer_scalar_factor") is not False
                    row["gate"] = (
                        "pass"
                        if row.get("status") == "ok"
                        and semantic_ok
                        and (shape_ok or args.gate_stage != "M3")
                        else "known_baseline_limitation"
                        if args.gate_stage == "baseline" and known_baseline
                        else "pending_M3"
                        if not milestone_required and semantic_ok
                        else "fail"
                    )
                    if row["gate"] == "fail":
                        failures.append(row)
                record["diagnostic_acceptance"] = {
                    "stage": args.gate_stage,
                    "pass": not failures,
                    "failed_rows": len(failures),
                }
            elif not isinstance(case, HistoricalLadder):
                raise ValueError("This diagnostic requires a historical ladder case")
            elif args.diagnostic == "validation":
                record["diagnostic"] = validation_cost(case)
            elif args.diagnostic == "partial-parse":
                record["diagnostic"] = partial_parse(case, args.term_limit)
            else:
                value, diagnostic = FactorProbe(case).run(
                    "atom" if args.diagnostic == "factor-atom" else "polynomial"
                )
                assert value == result
                record["diagnostic"] = {**diagnostic, "exact_reference": True}
    assert hashlib.sha256(core_path.read_bytes()).hexdigest() == before
    args.output.write_text(json.dumps(record, indent=2) + "\n")
    if record.get("diagnostic_acceptance", {}).get("pass") is False:
        raise AssertionError(
            "Required keep-the-rest milestone gate failed; full outcomes retained"
        )
    return record


def resource_snapshot(cpu):
    import resource

    usage = resource.getrusage(resource.RUSAGE_SELF)
    data = {
        "time_ns": __import__("time").time_ns(),
        "load_average": os.getloadavg(),
        "max_rss": usage.ru_maxrss,
        "voluntary_switches": usage.ru_nvcsw,
        "involuntary_switches": usage.ru_nivcsw,
    }
    if cpu is not None:
        siblings = Path(
            f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"
        )
        data["smt_siblings"] = (
            siblings.read_text().strip() if siblings.exists() else None
        )
        sibling_ids = (data["smt_siblings"] or str(cpu)).split(",")
        names = {f"cpu{n}" for n in sibling_ids if n.isdigit()}
        data["cpu_stat"] = [
            line
            for line in Path("/proc/stat").read_text().splitlines()
            if line.split()[0] in names
        ]
    return data


def case_form_script(case, directory):
    if hasattr(case, "form_source"):
        path = directory / "case.frm"
        path.write_text(case.form_source())
        return path
    name = (
        "fermion_propagator_ladder.frm"
        if case.particle == "fermionic"
        else "gluon_propagator_ladder.frm"
    )
    return Path(__file__).with_name(name)


def consolidation_benchmark(args):
    """M0 and later milestones use this scheduler around the existing case bodies."""
    import sys

    from tensor_benchmark_cases import CASES, make_case

    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    if args.rounds < 3:
        raise ValueError("M0 requires at least three rounds")
    interpreters = args.interpreter or [f"baseline={sys.executable}"]
    variants = [value.split("=", 1) for value in interpreters]
    if len({label for label, _ in variants}) != len(variants):
        raise ValueError("Interpreter labels must be unique")
    ladder_routes = dict(value.split("=", 1) for value in (args.ladder_route or []))
    if set(ladder_routes) - {label for label, _ in variants} or any(
        route not in ("typed", "factorized") for route in ladder_routes.values()
    ):
        raise ValueError("Use --ladder-route LABEL=typed or LABEL=factorized")
    cases = args.cases.split(",") if args.cases else CASES
    if any(case not in CASES for case in cases):
        raise ValueError(f"Cases must be selected from {CASES}")
    directory = args.output.with_suffix("")
    directory.mkdir(parents=True, exist_ok=False)
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    sources = {
        str(path): sha(path)
        for path in [
            Path(__file__),
            *Path(__file__).parent.glob("tensor_benchmark_*.py"),
            *(
                Path(__file__).with_name(name)
                for name in (
                    "gluon_ladder.py",
                    "fermion_ladder_validation.py",
                    "gluon_ladder_validation.py",
                    "gamma_simplification.py",
                    "fermion_propagator_ladder.frm",
                    "gluon_propagator_ladder.frm",
                )
            ),
        ]
    }
    fixture = (
        Path(__file__).resolve().parents[2]
        / "crates/idenso/tests/fixtures/aa_aa_2l_gl16_integrated_uv_start_after_simplify_metrics.sym"
    )
    sources[str(fixture)] = sha(fixture)
    if args.form:
        sources[str(Path(args.form).resolve())] = sha(args.form)
    provenance = Path(__file__).with_name("tensor_benchmark_provenance.json")
    sources[str(provenance)] = sha(provenance)
    report = {
        "schema": 1,
        "provenance": json.loads(provenance.read_text()),
        "platform": platform.platform(),
        "milestone": args.milestone,
        "rounds": args.rounds,
        "created_utc": datetime.now(UTC).isoformat(),
        "sources": sources,
        "interpreters": dict(variants),
        "ladder_routes": ladder_routes,
        "cases": cases,
        "timing": "Fresh interpreter for each case/route/round; one fixed unmeasured warmup, checked and disposed before clocks; fixed measured calls; setup/checks/disposal excluded; no discarded measured samples or retries",
        "diagnostic_scope": "Separate check processes; phase clocks/counts never substituted for primary clocks",
        "scheduling": "Rotate interpreter order each round; raw scalar and FORM references separate from paired criterion",
        "host_scope": "Pinned among cooperating jobs, not an exclusive-host claim",
        "records": [],
        "checks": [],
        "form": [],
        "summary": {},
    }
    (directory / "protocol.json").write_text(
        json.dumps(
            {
                **report,
                "calls": args.calls,
                "fixed_warmups": 1,
                "expected_cores": dict(
                    value.split("=", 1) for value in (args.core or [])
                ),
                "form": args.form,
                "form_batch_target_cpu_ns": 200_000_000,
                "free4_calibration": "Start8192 copies, double until>=200ms, cap1048576; retain all trials then fix3-round copy count",
                "form_only": args.form_only,
                "cpu": args.cpu,
            },
            indent=2,
        )
        + "\n"
    )
    cores = dict(value.split("=", 1) for value in (args.core or []))
    report["expected_cores"] = dict(cores)

    def child(label, python, case, route, suffix, check=False, measurement=False):
        if route == "typed" and (case.startswith(("historical-", "gluon-"))):
            route = ladder_routes.get(label, "typed")
        output = directory / f"{label}-{case}-{route}-{suffix}.json"
        command = [
            python,
            str(Path(__file__).resolve()),
            "--suite",
            "consolidation",
            "--child",
            "--case",
            case,
            "--route",
            route,
            "--calls",
            str(args.calls),
            "--output",
            str(output),
        ]
        if args.cpu is not None:
            command += ["--cpu", str(args.cpu)]
        if label in cores:
            command += ["--expected-core", cores[label]]
        if check:
            command += [
                "--check",
                "--gate-stage",
                "baseline" if label == variants[0][0] else args.gate_stage,
            ]
            if args.form:
                command += ["--form", args.form]
            if args.diagnostic:
                command += [
                    "--diagnostic",
                    args.diagnostic,
                    "--term-limit",
                    str(args.term_limit),
                ]
        before = resource_snapshot(args.cpu)
        try:
            run = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as error:

            def decode(value):
                return (
                    value.decode(errors="replace")
                    if isinstance(value, bytes)
                    else (value or "")
                )

            run = subprocess.CompletedProcess(
                command,
                124,
                decode(error.stdout),
                decode(error.stderr) + "\nPreset timeout reached; no retry",
            )
        output.with_suffix(".stdout").write_text(run.stdout)
        output.with_suffix(".stderr").write_text(run.stderr)
        receipt = {
            "label": label,
            "case": case,
            "route": route,
            "command": command,
            "exit_code": run.returncode,
            "output": str(output),
            "before": before,
            "after": resource_snapshot(args.cpu),
        }
        report["checks" if check and not measurement else "records"].append(receipt)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        if run.returncode:
            raise RuntimeError(
                f"Child failed; retained output: {output.with_suffix('.stderr')}"
            )
        value = json.loads(output.read_text())
        receipt["record"] = value
        cores.setdefault(label, value["core_sha256"])
        assert cores[label] == value["core_sha256"]
        return value

    for case in cases:
        if not args.form_only:
            expected = None
            for label, python in variants:
                value = child(label, python, case, "typed", "check", check=True)
                expression = Path(value["expression"]).read_text()
                if expected is None:
                    expected = expression
                assert expression == expected, f"Interpreter output differs for {case}"
            if case.startswith("historical-"):
                value = child(*variants[0], case, "scalar", "check", check=True)
                assert Path(value["expression"]).read_text() == expected
            if args.check_only:
                continue
            for round_index in range(args.rounds):
                offset = round_index % len(variants)
                for label, python in variants[offset:] + variants[:offset]:
                    value = child(
                        label,
                        python,
                        case,
                        "typed",
                        f"round{round_index + 1}",
                        check=bool(args.diagnostic),
                        measurement=True,
                    )
                    assert Path(value["expression"]).read_text() == expected
                if case.startswith("historical-"):
                    value = child(
                        *variants[0], case, "scalar", f"round{round_index + 1}"
                    )
                    assert Path(value["expression"]).read_text() == expected
        if args.form and case != "production-aa-aa":
            reference = make_case(case)
            with tempfile.TemporaryDirectory(prefix="r3-form-") as temp:
                script = case_form_script(reference, Path(temp))
                diagnostic = None
                if args.form_only:
                    diagnostic = _form_run(
                        args.form,
                        script,
                        temp,
                        mode=reference.mode,
                        loops=reference.loops,
                        diagnostic=True,
                        order=reference.form_order,
                    )
                    polynomial = directory / f"{case}.form.txt"
                    polynomial.write_text(diagnostic.pop("polynomial"))
                    diagnostic.update(
                        polynomial_path=str(polynomial),
                        polynomial_sha256=sha(polynomial),
                    )
                pilot = _form_run(
                    args.form,
                    script,
                    temp,
                    mode=reference.mode,
                    loops=reference.loops,
                    repeats=1,
                    order=reference.form_order,
                )
                repeats = min(
                    8192,
                    max(
                        1, math.ceil(200_000_000 / max(pilot["body_cpu_ns"], 1_000_000))
                    ),
                )
                row = {
                    "case": case,
                    "diagnostic": diagnostic,
                    "pilot": pilot,
                    "calibration": [],
                    "samples": [],
                }
                report["form"].append(row)
                if case.startswith("free4-"):
                    repeats = 8192
                    while True:
                        calibration = _form_run(
                            args.form,
                            script,
                            temp,
                            mode=reference.mode,
                            loops=reference.loops,
                            repeats=repeats,
                            order=reference.form_order,
                        )
                        row["calibration"].append(calibration)
                        args.output.write_text(json.dumps(report, indent=2) + "\n")
                        if calibration["body_cpu_ns"] >= 200_000_000:
                            break
                        if repeats >= 1_048_576:
                            raise RuntimeError(
                                "FORM calibration exceeds the preset copy limit"
                            )
                        repeats *= 2
                row["fixed_repeats"] = repeats
                for _ in range(args.rounds):
                    sample = _form_run(
                        args.form,
                        script,
                        temp,
                        mode=reference.mode,
                        loops=reference.loops,
                        repeats=repeats,
                        order=reference.form_order,
                    )
                    row["samples"].append(sample)
                    if diagnostic is not None:
                        assert (
                            sample["term_counts"][-1] == diagnostic["term_counts"][-1]
                        )
                    args.output.write_text(json.dumps(report, indent=2) + "\n")
                    if sample["body_cpu_ns"] == 0:
                        raise RuntimeError(
                            "FORM body below timer resolution; sample retained, no retry"
                        )
        rows = [r for r in report["records"] if r["case"] == case]
        summary = {}
        for label, _ in [] if args.form_only else variants:
            own = [
                r["record"]
                for r in rows
                if r["label"] == label and r["route"] != "scalar"
            ]
            summary[label] = {
                clock: median(r[clock] / r["calls"] for r in own)
                for clock in ("wall_ns", "process_ns", "thread_ns")
            }
            if own and own[0]["route"] == "factorized":
                summary[label]["factorized"] = {
                    "phase_cpu_ns": {
                        phase: median(
                            next(
                                stage["cpu_ns"]
                                for stage in row["phases"]["stages"]
                                if stage["phase"] == phase
                            )
                            for row in own
                        )
                        for phase in (
                            "rule_application_and_admission",
                            "contraction",
                            "explicit_materialization",
                        )
                    },
                    "aliased_bytes": [row["phases"]["aliased_bytes"] for row in own],
                    "definitions": [row["phases"]["definitions"] for row in own],
                    "scope": "Separate diagnostic calls; primary end-to-end clocks above",
                }
            if label != variants[0][0]:
                base = [
                    r["record"]
                    for r in rows
                    if r["label"] == variants[0][0] and r["route"] != "scalar"
                ]
                summary[label]["paired_ratios"] = {
                    clock: [
                        c[clock] / b[clock] if b[clock] else None
                        for b, c in zip(base, own, strict=True)
                    ]
                    for clock in ("wall_ns", "process_ns", "thread_ns")
                }
                summary[label]["faster_pairs"] = {
                    clock: sum(value is not None and value < 1 for value in ratios)
                    for clock, ratios in summary[label]["paired_ratios"].items()
                }
                summary[label]["median_paired_ratio"] = {
                    clock: median(values)
                    if all(x is not None for x in values)
                    else None
                    for clock, values in summary[label]["paired_ratios"].items()
                }
        scalar = [r["record"] for r in rows if r["route"] == "scalar"]
        if scalar:
            summary["pure_symbolica_reference"] = {
                clock: median(r[clock] / r["calls"] for r in scalar)
                for clock in ("wall_ns", "process_ns", "thread_ns")
            }
        if report["form"] and report["form"][-1]["case"] == case:
            form_ns = median(
                r["body_cpu_ns"] / r["repeats"] for r in report["form"][-1]["samples"]
            )
            summary["FORM"] = {
                "body_cpu_ns": form_ns,
                "scope": "Independent batch expression; excludes startup/setup/checks/disposal",
            }
            for label, _ in [] if args.form_only else variants:
                summary[label]["FORM_cpu_ratio"] = (
                    summary[label]["process_ns"] / form_ns
                )
        if args.diagnostic:
            summary["diagnostics"] = {}
            for label, _ in variants:
                samples = [
                    r["record"]["diagnostic"]
                    for r in rows
                    if r["label"] == label and r["route"] == "typed"
                ]
                if samples and isinstance(samples[0], list):
                    medians = []
                    for observations in zip(*samples, strict=True):
                        row = {
                            k: v
                            for k, v in observations[0].items()
                            if k
                            in (
                                "stage",
                                "operation",
                                "case",
                                "route",
                                "terms",
                                "requested_terms",
                                "method",
                                "branches",
                                "spinor_slots",
                            )
                        }
                        for key in ("cpu_ns", "whole_over_individual"):
                            if all(key in o for o in observations):
                                row[key] = median(o[key] for o in observations)
                        row["statuses"] = [o.get("status", "ok") for o in observations]
                        medians.append(row)
                    summary["diagnostics"][label] = medians
                else:
                    summary["diagnostics"][label] = samples
        report["summary"][case] = summary
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    assert all(sha(path) == digest for path, digest in sources.items()), (
        "Source changed during cohort"
    )
    for receipt in report["checks"] + report["records"]:
        value = receipt["record"]
        assert sha(value["core"]) == value["core_sha256"], "Core changed during cohort"
    report["acceptance"] = {}
    if args.diagnostic == "partial-parse":
        for label, _ in variants[1:]:
            checks = [
                row
                for case in report["summary"].values()
                for row in case["diagnostics"][label]
                if row.get("operation") == "network_identical_termset"
            ]
            report["acceptance"][label] = {
                "gate": "R1 identical-termset whole/single median <= 2",
                "pass": bool(checks)
                and all(
                    row.get("whole_over_individual", float("inf")) <= 2
                    for row in checks
                ),
                "rows": checks,
            }
    if args.milestone == "M1" and "historical-early" in report["summary"]:
        for label, _ in variants:
            if ladder_routes.get(label) != "factorized":
                continue
            value = report["summary"]["historical-early"][label]["factorized"]
            clocks = value["phase_cpu_ns"]
            gates = {
                "contraction_under_50ms": clocks["contraction"] < 50_000_000,
                "materialization_under_150ms": clocks["explicit_materialization"]
                < 150_000_000,
                "aliases_under_100KB": max(value["aliased_bytes"]) < 100_000,
            }
            report["acceptance"][label] = {
                "gate": "M1 historical ladder; plan's reverse order = 5,4,6,3,7,2,8,1",
                "pass": all(gates.values()),
                "requirements": gates,
                "measurements": value,
            }
    report.update(core_identities=cores, source_drift=False, completed=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    if any(not gate["pass"] for gate in report["acceptance"].values()):
        raise AssertionError(
            "Required milestone acceptance failed; complete cohort retained"
        )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite", choices=("ladder", "consolidation"), default="ladder"
    )
    parser.add_argument(
        "--interpreter",
        action="append",
        help="LABEL=/path/to/release/python; repeat for paired milestones",
    )
    parser.add_argument(
        "--core",
        action="append",
        help="LABEL=expected core SHA256; binds release identity",
    )
    parser.add_argument("--milestone", default="M0")
    parser.add_argument(
        "--ladder-route",
        action="append",
        help="LABEL=typed or LABEL=factorized for historical/gluonic ladders only",
    )
    parser.add_argument("--cases", help="Comma-separated consolidation case names")
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument(
        "--form-only",
        action="store_true",
        help="Measure only the separate FORM reference; do not resample Python",
    )
    parser.add_argument("--gate-stage", choices=("baseline", "M0", "M3"), default="M0")
    parser.add_argument("--calls", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--diagnostic",
        choices=(
            "validation",
            "partial-parse",
            "keep-rest",
            "factor-atom",
            "factor-polynomial",
        ),
    )
    parser.add_argument("--term-limit", type=int, default=200)
    parser.add_argument("--child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--case", help=argparse.SUPPRESS)
    parser.add_argument(
        "--route",
        choices=("typed", "scalar", "factorized"),
        default="typed",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--expected-core", help=argparse.SUPPRESS)
    parser.add_argument("--check", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--form")
    parser.add_argument(
        "--particle", choices=("fermionic", "gluonic"), default="fermionic"
    )
    parser.add_argument("--loops", type=int, choices=(3, 4), default=3)
    parser.add_argument("--rounds", type=int)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.rounds is None:
        args.rounds = 3 if args.suite == "consolidation" else 5
    if args.suite == "consolidation":
        if args.form_only and (not args.form or args.check_only or args.diagnostic):
            parser.error("--form-only requires --form and excludes checks/diagnostics")
        if args.calls < 1 or args.term_limit < 1:
            parser.error("--calls and --term-limit must be positive")
        result = (
            consolidation_child(args) if args.child else consolidation_benchmark(args)
        )
        print(
            json.dumps(
                result.get("summary", {"case": args.case, "status": "ok"}), indent=2
            )
        )
        raise SystemExit(0)
    options = {}
    if args.particle == "gluonic":
        if args.loops != 4:
            parser.error("The gluonic comparison requires --loops 4")
        from gluon_ladder import GluonLadder

        options = {
            "case_factory": GluonLadder,
            "form_script": Path(__file__).with_name("gluon_propagator_ladder.frm"),
        }
    result = benchmark(
        args.form,
        loops=args.loops,
        rounds=args.rounds,
        output=args.output,
        cpu=args.cpu,
        **options,
    )
    for case in result["cases"]:
        print(
            f"{case['dimension']}: Idenso {case['idenso_median_cpu_ms']:.6f} ms, FORM {case['form_median_cpu_ms']:.6f} ms, {case['cpu_ratio']:.2f}x"
        )
