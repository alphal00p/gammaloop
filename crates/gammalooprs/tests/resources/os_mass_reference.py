"""Numerically check one private, natively certified QCD mass-reference input.

The Rust fixture owns graph generation, tensor/color algebra, phases, source
weights and scaleless certificates. This consumer checks their finite scalar
serialization and evaluates the independent reference with a pinned backend.
"""

import argparse
import hashlib
import itertools
import json
import math
import multiprocessing
from pathlib import Path
import re
import subprocess
import sys
import time

import sympy as sp

EPS = sp.Symbol("eps")
D = sp.Symbol("d")
CHANNELS = {
    "D0": {"CF"},
    "D00": {"CFCA"},
    "D01": {"CFCA"},
    "D02": {"CFCA"},
    "D03": {"CFCA"},
    **{f"D0{i}": {"CFTF"} for i in range(4, 9)},
    "D09": {"CF2"},
    "D10": {"CFTF"},
    "D11": {"CF2", "CFCA"},
}
POLICY = {
    "generation_processes": 2,
    "form_threads": 2,
    "make_jobs": 2,
    "qmc_cputhreads": 2,
    "epsrel": 1e-7,
    "epsabs": 1e-99,
    "minn": 10000,
    "maxeval": 1000000,
    "seed": 1337,
    "comparison_bound": "1e-6*(1+abs(reference))",
}


def polynomial(text, symbols):
    """Accept only exact arithmetic in the named finite scalar coordinates."""
    assert isinstance(text, str) and re.fullmatch(r"[A-Za-z0-9_+*/^()\s-]+", text)
    assert set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text)) <= symbols.keys()
    value = sp.sympify(text.replace("^", "**"), locals=symbols)
    assert not value.atoms(sp.Float)
    return value


def prepare(source):
    loops = source["loops"]
    assert loops in (1, 2) and source["external_momenta"] == ["p"]
    momenta = source["loop_momenta"] + ["p"]
    variables = [entry["variable"] for entry in source["scalar_products"]]
    assert len(source["loop_momenta"]) == loops
    assert len(set(momenta + variables + ["d", "eps"])) == len(momenta + variables) + 2
    assert all(
        re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", name) for name in momenta + variables
    )
    symbols = {name: sp.Symbol(name) for name in momenta + variables + ["d", "eps"]}
    dots = {}
    for entry in source["scalar_products"]:
        pair = tuple(sorted((entry["left"], entry["right"])))
        assert all(name in momenta for name in pair) and pair not in dots
        dots[pair] = symbols[entry["variable"]]
    assert ("p", "p") not in dots
    dots[("p", "p")] = sp.Integer(1)
    assert set(dots) == {
        tuple(sorted(pair)) for pair in itertools.product(momenta, repeat=2)
    }
    propagators = []
    edges = set()
    for pole in source["propagators"]:
        assert pole["power"] == 1 and pole["mass_squared"] in ("0", "1")
        assert pole["edge"] not in edges
        edges.add(pole["edge"])
        momentum = polynomial(
            pole["momentum"], {name: symbols[name] for name in momenta}
        )
        linear = sp.Poly(momentum, *(symbols[name] for name in momenta))
        assert linear.total_degree() <= 1 and linear.coeff_monomial(1) == 0
        coefficients = {name: linear.coeff_monomial(symbols[name]) for name in momenta}
        assert all(coefficient.is_Rational for coefficient in coefficients.values())
        reconstructed = sum(
            coefficients[a] * coefficients[b] * dots[tuple(sorted((a, b)))]
            for a, b in itertools.product(momenta, repeat=2)
        ) - int(pole["mass_squared"])
        expected = polynomial(
            pole["polynomial"], {name: symbols[name] for name in variables}
        )
        assert sp.Poly(
            reconstructed - expected, *(symbols[name] for name in variables)
        ).is_zero
        propagators.append(
            f"({pole['momentum'].replace('^', '**')})**2-({pole['mass_squared']})"
        )
    rows = {}
    for row in source["color_rows"]:
        assert row["channel"] not in rows
        coefficient = polynomial(
            row["coefficient"], {name: symbols[name] for name in variables + ["d"]}
        )
        assert not sp.Poly(
            coefficient, D, *(symbols[name] for name in variables)
        ).is_zero
        rows[row["channel"]] = coefficient
    assert set(rows) == CHANNELS[source["name"]]
    return {
        "source": source,
        "symbols": symbols,
        "dots": dots,
        "rows": rows,
        "propagators": propagators,
    }


def backend_input(prepared, coefficient, powers=None):
    source = prepared["source"]
    variables = [entry["variable"] for entry in source["scalar_products"]]
    symbols = prepared["symbols"]
    coefficient = coefficient.subs(D, 4 - 2 * EPS)
    finite = sp.Poly(coefficient, *(symbols[name] for name in variables))
    pairs = {
        entry["variable"]: (entry["left"], entry["right"])
        for entry in source["scalar_products"]
    }
    terms, indices, rebuilt = [], [], sp.Integer(0)
    # Only the already contracted finite scalar coefficient is materialized.
    for exponents, weight in finite.terms():
        sp.Poly(weight, EPS)
        factors, monomial = [f"({sp.sstr(weight)})"], weight
        for name, exponent in zip(variables, exponents, strict=True):
            monomial *= symbols[name] ** exponent
            for _ in range(exponent):
                index = f"mu{len(indices)}"
                indices.append(index)
                left, right = pairs[name]
                factors.append(f"{left}({index})*{right}({index})")
        terms.append("*".join(factors))
        rebuilt += monomial
    assert sp.Poly(
        rebuilt - coefficient, *(symbols[name] for name in variables)
    ).is_zero
    return {
        "propagators": prepared["propagators"],
        "loop_momenta": source["loop_momenta"],
        "external_momenta": ["p"],
        "Lorentz_indices": indices,
        "numerator": "+".join(terms),
        "dimensionality": "4-2*eps",
        "powerlist": powers or [1] * len(prepared["propagators"]),
        "replacement_rules": [("p*p", "1")],
    }


def integrate(name, inputs, directory, order):
    import pySecDec as psd
    from pySecDec.integral_interface import IntegralLibrary, series_to_json

    assert psd.__version__ == "1.6.6"
    directory.mkdir()
    loops = len(inputs["loop_momenta"])
    (directory / "input.json").write_text(
        json.dumps(
            {"inputs": inputs, "policy": POLICY, "requested_order": order}, indent=2
        )
    )
    started = time.monotonic()
    # PySecDec writes relative paths during generation. The caller runs cases
    # serially and restores cwd after this concrete package is complete.
    import os

    original = Path.cwd()
    try:
        os.chdir(directory)
        integral = psd.loop_integral.LoopIntegralFromPropagators(**inputs)
        psd.loop_integral.loop_package(
            name="coefficient",
            loop_integral=integral,
            real_parameters=[],
            requested_orders=[order],
            contour_deformation=False,
            form_optimization_level=2,
            form_threads=2,
            processes=2,
            form_work_space="100M",
            decomposition_method="iterative",
            additional_prefactor=f"exp({loops}*EulerGamma*eps)",
        )
        generated = time.monotonic()
        subprocess.run(["make", "-j2"], cwd="coefficient", check=True)
        built = time.monotonic()
        library = directory / "coefficient/coefficient_pylink.so"
        evaluate = IntegralLibrary(str(library))
        evaluate.use_Qmc(
            transform="Korobov3",
            epsrel=1e-7,
            epsabs=1e-99,
            maxeval=1000000,
            minn=10000,
            cputhreads=2,
            seed=1337,
        )
        raw = evaluate()
        (directory / "raw-result.json").write_text(json.dumps(raw, indent=2))
        series = series_to_json(raw[2])
        assert series["regulators"] == ["eps"] and set(series["sums"]) == {"sum"}
        rows = {}
        for power, (value, error) in series["sums"]["sum"].items():
            assert len(power) == 1 and power[0] not in rows
            values = [value.real, value.imag, error.real, error.imag]
            assert all(math.isfinite(number) for number in values)
            rows[power[0]] = (value, complex(abs(error.real), abs(error.imag)))
        assert set(range(-loops, order + 1)) <= rows.keys()
        report = {
            "name": name,
            "generation_seconds": generated - started,
            "build_seconds": built - generated,
            "evaluation_seconds": time.monotonic() - built,
            "coefficients": {
                str(power): {
                    "value": [value.real, value.imag],
                    "error": [error.real, error.imag],
                }
                for power, (value, error) in sorted(rows.items())
            },
            "library_sha256": hashlib.sha256(library.read_bytes()).hexdigest(),
            "input_sha256": hashlib.sha256(
                (directory / "input.json").read_bytes()
            ).hexdigest(),
        }
        (directory / "summary.json").write_text(json.dumps(report, indent=2))
        return rows, report
    finally:
        os.chdir(original)


def compare(scope, actual, expected, checks):
    for power in sorted(set(actual) | set(expected)):
        value, error = actual.get(power, (0j, 0j))
        reference = complex(expected.get(power, 0))
        bound = 1e-6 * (1 + abs(reference))
        passed = (
            abs(value.real - reference.real) <= bound
            and abs(value.imag - reference.imag) <= bound
        )
        checks.append(
            {
                "scope": scope,
                "epsilon_power": power,
                "actual": [value.real, value.imag],
                "expected": [reference.real, reference.imag],
                "sum_reported_absolute_errors": [error.real, error.imag],
                "absolute_tolerance": bound,
                "passed": passed,
            }
        )


def combine(terms):
    result = {}
    for weight, rows in terms:
        for power, (value, error) in rows.items():
            old_value, old_error = result.get(power, (0j, 0j))
            result[power] = (
                old_value + weight * value,
                old_error + abs(weight) * error,
            )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source_path, output = args.input.resolve(), args.output.resolve()
    assert not output.exists(), "never overwrite an earlier reference run"
    data = json.loads(source_path.read_text())
    sources = data["sources"]
    assert len(sources) == 13 and {s["name"] for s in sources} == set(CHANNELS)
    assert all(s["loops"] == (1 if s["name"] == "D0" else 2) for s in sources)
    assert [s["name"] for s in sources if s["integral_zero_certified"]] == ["D03"]
    prepared = {source["name"]: prepare(source) for source in sources}
    light = [source for source in sources if source["light_flavor"] is not None]
    assert len(light) == 5 and {source["light_flavor"] for source in light} == set(
        range(1, 6)
    )
    assert {source["name"] for source in light} == {f"D0{i}" for i in range(4, 9)}
    for source in light[1:]:
        assert {
            k: v for k, v in source.items() if k not in ("name", "light_flavor")
        } == {k: v for k, v in light[0].items() if k not in ("name", "light_flavor")}
    crossed = prepared["D11"]
    assert sp.Poly(
        crossed["rows"]["CF2"] + 2 * crossed["rows"]["CFCA"],
        D,
        *(
            crossed["symbols"][entry["variable"]]
            for entry in crossed["source"]["scalar_products"]
        ),
    ).is_zero

    one = prepared["D0"]
    k = one["source"]["loop_momenta"][0]
    kk, pk = one["dots"][(k, k)], one["dots"][tuple(sorted((k, "p")))]
    poles = [
        polynomial(pole["polynomial"], one["symbols"])
        for pole in one["source"]["propagators"]
    ]
    assert (
        len(poles) == 2
        and sp.expand(poles[0] - (kk - 1)) == 0
        and sp.expand(poles[1] - (kk - 2 * pk + 1)) == 0
    )
    assert sp.expand(one["rows"]["CF"] - ((2 - D) * pk + D)) == 0
    # Independent fixed-external internal-M derivative, then M=m=1. The native
    # fixture certifies this identity before serializing the one-loop source.
    insertion = D * (kk - 1) + 2 * ((2 - D) * pk + D)
    cases = [
        ("native-one-loop", backend_input(one, one["rows"]["CF"]), 1),
        ("fixed-external-mass-insertion", backend_input(one, insertion, [2, 1]), 1),
    ]
    for name in ("D00", "D01", "D02", "D04", "D09", "D10", "D11"):
        channel = (
            "CF2"
            if name in ("D09", "D11")
            else "CFTF"
            if name in ("D04", "D10")
            else "CFCA"
        )
        cases.append(
            (name, backend_input(prepared[name], prepared[name]["rows"][channel]), 0)
        )
    assert len(cases) == 9
    f = {-1: sp.Integer(3), 0: sp.Integer(4), 1: 8 + sp.pi**2 / 4}
    g = {-1: sp.Integer(6), 0: sp.Integer(2), 1: 8 + sp.pi**2 / 2}
    bare = {
        "CF2": {
            -2: sp.Rational(27, 2),
            -1: sp.Rational(75, 4),
            0: sp.Rational(441, 8)
            + 29 * sp.pi**2 / 4
            + 12 * sp.zeta(3)
            - 8 * sp.pi**2 * sp.log(2),
        },
        "CFCA": {
            -2: sp.Rational(11, 2),
            -1: sp.Rational(91, 4),
            0: sp.Rational(605, 8)
            - 5 * sp.pi**2 / 12
            - 6 * sp.zeta(3)
            + 4 * sp.pi**2 * sp.log(2),
        },
        "CFTFnl": {-2: -2, -1: -7, 0: -sp.Rational(45, 2) - 5 * sp.pi**2 / 3},
        "CFTFnh": {-2: -2, -1: -7, 0: -sp.Rational(69, 2) + 7 * sp.pi**2 / 3},
    }
    zm = {
        "CF2": {
            -2: sp.Rational(9, 2),
            -1: sp.Rational(45, 4),
            0: sp.Rational(199, 8)
            - 17 * sp.pi**2 / 4
            - 12 * sp.zeta(3)
            + 8 * sp.pi**2 * sp.log(2),
        },
        "CFCA": {
            -2: sp.Rational(11, 2),
            -1: -sp.Rational(97, 12),
            0: -sp.Rational(1111, 24)
            + 4 * sp.pi**2 / 3
            + 6 * sp.zeta(3)
            - 4 * sp.pi**2 * sp.log(2),
        },
        "CFTFnl": {
            -2: -2,
            -1: sp.Rational(5, 3),
            0: sp.Rational(71, 6) + 4 * sp.pi**2 / 3,
        },
        "CFTFnh": {
            -2: -2,
            -1: sp.Rational(5, 3),
            0: sp.Rational(143, 6) - 8 * sp.pi**2 / 3,
        },
    }
    f_atom, g_atom = (
        sum(value * EPS**power for power, value in rows.items()) for rows in (f, g)
    )
    insertions = {
        "CF2": f_atom * g_atom,
        "CFCA": sp.Rational(11, 3) * f_atom / EPS,
        "CFTFnl": -sp.Rational(4, 3) * f_atom / EPS,
        "CFTFnh": -sp.Rational(4, 3) * f_atom / EPS,
    }
    for channel in bare:
        for power in (-2, -1, 0):
            assert (
                sp.simplify(
                    -bare[channel][power]
                    + sp.expand(insertions[channel]).coeff(EPS, power)
                    - zm[channel][power]
                )
                == 0
            )
    output.mkdir(parents=True)
    validation = {
        "input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_inventory": [source["name"] for source in sources],
        "native_sources": 13,
        "two_loop_sources": 12,
        "numeric_cases": [name for name, _, _ in cases],
        "scalar_and_denominator_serialization_exact": True,
        "exact_light_and_crossed_reuse": True,
        "retained_scaleless_zero": "D03",
        "fixed_external_anchor_checked": True,
        "reference_conversion_exact": True,
        "integration_performed": False,
        "policy": POLICY,
        "python": sys.version,
    }
    (output / "validation.json").write_text(json.dumps(validation, indent=2))
    if args.validate_only:
        return
    multiprocessing.set_start_method("fork")
    evaluations = {
        name: integrate(name, inputs, output / name, order)
        for name, inputs, order in cases
    }
    results = {name: result[0] for name, result in evaluations.items()}
    checks = []
    compare("one_loop", results["native-one-loop"], f, checks)
    compare(
        "fixed_external_mass_insertion",
        results["fixed-external-mass-insertion"],
        g,
        checks,
    )
    channels = {
        "CF2": combine([(1, results["D09"]), (1, results["D11"])]),
        "CFCA": combine(
            [
                (1, results["D00"]),
                (1, results["D01"]),
                (1, results["D02"]),
                (-0.5, results["D11"]),
            ]
        ),
        "CFTFnl": results["D04"],
        "CFTFnh": results["D10"],
    }
    reconstructed = {}
    for channel, rows in channels.items():
        compare("bare_" + channel, rows, bare[channel], checks)
        reconstructed[channel] = {
            power: (
                -value + complex(sp.expand(insertions[channel]).coeff(EPS, power)),
                error,
            )
            for power, (value, error) in rows.items()
        }
        compare("Zm_" + channel, reconstructed[channel], zm[channel], checks)
    weights = {
        "CF2": sp.Rational(16, 9),
        "CFCA": sp.Integer(4),
        "CFTFnl": sp.Rational(10, 3),
        "CFTFnh": sp.Rational(2, 3),
    }
    for scope, rows, reference in (
        ("SU3_bare", channels, bare),
        ("SU3_Zm", reconstructed, zm),
    ):
        expected = {
            power: sum(
                weights[channel] * reference[channel][power] for channel in weights
            )
            for power in (-2, -1, 0)
        }
        compare(
            scope,
            combine([(float(weights[channel]), rows[channel]) for channel in weights]),
            expected,
            checks,
        )
    summary = {
        **validation,
        "integration_performed": True,
        "one_loop_anchors_passed": all(
            row["passed"]
            for row in checks
            if row["scope"] in ("one_loop", "fixed_external_mass_insertion")
        ),
        "bare_four_color_reference_passed": all(
            row["passed"] for row in checks if row["scope"].startswith("bare_")
        ),
        "independent_insertions_reconstruct_Zm_reference": all(
            row["passed"] for row in checks if row["scope"].startswith("Zm_")
        ),
        "SU3_reference_passed": all(
            row["passed"] for row in checks if row["scope"].startswith("SU3_")
        ),
        "production_OS_forest_integrated_condition_evaluated": False,
        "checks": checks,
        "backend_cases": {name: result[1] for name, result in evaluations.items()},
        "insertion_source": "Known analytic one-loop f, g and beta0; not production counterterm integration",
        "uncertainty_policy": "Linear absolute error sums; no Gaussian confidence claim",
        "reference_assembly_note": "Zm reconstruction is the bare comparison with known insertions, not a second independent two-loop integration",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    assert all(row["passed"] for row in checks), (
        "independent mass-reference comparison failed"
    )


if __name__ == "__main__":
    main()
