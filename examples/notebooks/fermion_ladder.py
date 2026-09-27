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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--form")
    parser.add_argument(
        "--particle", choices=("fermionic", "gluonic"), default="fermionic"
    )
    parser.add_argument("--loops", type=int, choices=(3, 4), default=3)
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
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
