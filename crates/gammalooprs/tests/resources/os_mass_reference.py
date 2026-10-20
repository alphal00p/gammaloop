"""Check certified QCD mass inputs and an optional separate D09 vector input.

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
import re
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp
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
# Retained run of this same pinned backend/policy after exact native D09 vector
# admission, independent of production counterterm construction. This measured
# finite datum is not an analytic constant or a four-dimensional IR proof.
D09_VECTOR_NUMERICAL_REFERENCE = {
    "finite_prepared": "-19.0743405536087186",
    "reported_absolute_error": "6.386262723980348e-7",
    "source_revision": "c9532785fbde412355d762bb81020624bb8799f3",
    "closed_receipt_sha256": "615bbd0efbd6e9fe0c9f9531ca885ac7ead2e857c57636cf05697909ecc99ecb",
    "normalization": "a^2*CF^2, amputated Sigma, p^2=m^2=mu_bar^2=1",
}
# Canonical JSON (sorted keys, separators=(',', ':')) of the qualified c953
# backend input dictionary, not its outer policy/requested-order envelope.
# The retained input file SHA is 1c3552b96fe27d6d9a696a2c99e80cbd43ce77750270b6de1b6fa51056d89440.
D09_VECTOR_BACKEND_INPUT_SHA256 = (
    "01ecb7d5c3f696f7c7d70344f846d7233a5c92075d5e5165e82b992a812d6b20"
)


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


def prepare_vector(payload, mass, gluon_mass_squared="0"):
    """Bind a separate A projection to the original D09 chart and full d trace."""
    assert payload["projection"] == "pslash/(4*Nc*p_squared)"
    assert payload["source_kernel_weight"] == "1"
    assert payload["external_ordering"] == "-1"
    assert payload["numerator_prefactor"] == "1"
    assert payload["propagator_prescription"] == "Feynman+i0"
    vector = prepare(payload["source"])
    source = vector["source"]
    assert source["name"] == "D09" and source["loops"] == 2
    assert source["light_flavor"] is None and not source["integral_zero_certified"]
    # Only the finite scalar numerator differs. Keep every original edge,
    # repeated line, signed momentum, mass, power and scalar-product binding.
    assert gluon_mass_squared in ("0", "1")
    if gluon_mass_squared == "0":
        assert {k: v for k, v in source.items() if k != "color_rows"} == {
            k: v for k, v in mass["source"].items() if k != "color_rows"
        }
    else:
        assert {
            k: v for k, v in source.items() if k not in ("color_rows", "propagators")
        } == {
            k: v
            for k, v in mass["source"].items()
            if k not in ("color_rows", "propagators")
        }
        for pole, original in zip(
            source["propagators"], mass["source"]["propagators"], strict=True
        ):
            expected = original.copy()
            if expected["edge"] in (3, 6):
                assert expected["mass_squared"] == "0"
                expected["mass_squared"] = "1"
            assert {k: v for k, v in pole.items() if k != "polynomial"} == {
                k: v for k, v in expected.items() if k != "polynomial"
            }
    l, k = source["loop_momenta"]
    p = vector["symbols"]["p"]
    loop_symbols = {name: vector["symbols"][name] for name in (l, k, "p")}
    l, k = (loop_symbols[name] for name in (l, k))
    for pole, (edge, momentum, mass_squared) in zip(
        source["propagators"],
        [
            (2, l, "1"),
            (3, k - l, gluon_mass_squared),
            (4, k, "1"),
            (5, k, "1"),
            (6, k - p, gluon_mass_squared),
        ],
        strict=True,
    ):
        assert pole["edge"] == edge and pole["mass_squared"] == mass_squared
        assert sp.Poly(
            polynomial(pole["momentum"], loop_symbols) - momentum, l, k, p
        ).is_zero
    dots = vector["dots"]
    kk = dots[(str(k), str(k))]
    kl = dots[tuple(sorted((str(k), str(l))))]
    pk = dots[tuple(sorted((str(k), "p")))]
    pl = dots[tuple(sorted((str(l), "p")))]
    # Tr(pslash N)/(4 Nc p^2), N=gamma_mu(k/+m)gamma_nu(l/+m)
    # gamma^nu(k/+m)gamma^mu. No external Dirac equation and no A+B inference.
    expected = (2 - D) ** 2 * (2 * kl * pk + (1 - kk) * pl) + 2 * D * (2 - D) * pk
    variables = [
        vector["symbols"][entry["variable"]] for entry in source["scalar_products"]
    ]
    assert sp.Poly(vector["rows"]["CF2"] - expected, D, *variables).is_zero
    return vector


def vector_protocol_controls(payload, mass, gluon_mass_squared="0"):
    """Reject concrete source/protocol errors before any backend package exists."""
    mutations = {
        "early_D4": lambda value: value["source"]["color_rows"][0].update(
            coefficient=re.sub(
                r"\bd\b", "4", value["source"]["color_rows"][0]["coefficient"]
            )
        ),
        "mass_projector": lambda value: value["source"].update(
            color_rows=mass["source"]["color_rows"]
        ),
        "duplicate_source_edge": lambda value: value["source"]["propagators"][1].update(
            edge=2
        ),
        "wrong_kernel_sign": lambda value: value.update(source_kernel_weight="-1"),
        "wrong_external_ordering": lambda value: value.update(external_ordering="1"),
        "wrong_prescription": lambda value: value.update(
            propagator_prescription="principal_value"
        ),
    }
    for mutate in mutations.values():
        changed = json.loads(json.dumps(payload))
        mutate(changed)
        try:
            prepare_vector(changed, mass, gluon_mass_squared)
        except AssertionError:
            continue
        raise AssertionError("malformed optional vector source was admitted")
    return list(mutations)


def regulated_vector_reference(
    vector,
    mass,
    source_path,
    output,
    validate_only,
    emitted_regulated_source=None,
    emitted_unregulated_source=None,
):
    """Independent dimensional reference; never the emitted CFF integral.

    A unit positive scalar mass changes only the original gluon poles. The
    unchanged native ordered numerator and source phase remain the anchor.
    """
    payload = json.loads(json.dumps(vector))
    symbols = prepare(payload["source"])["symbols"]
    for pole in payload["source"]["propagators"]:
        if pole["edge"] in (3, 6):
            assert pole["mass_squared"] == "0"
            pole["mass_squared"] = "1"
            pole["polynomial"] = sp.sstr(polynomial(pole["polynomial"], symbols) - 1)
    regulated = prepare_vector(payload, mass, "1")
    controls = vector_protocol_controls(payload, mass, "1")
    for edge in (3, 6):
        wrong = json.loads(json.dumps(payload))
        pole = next(p for p in wrong["source"]["propagators"] if p["edge"] == edge)
        pole["mass_squared"] = "0"
        pole["polynomial"] = sp.sstr(polynomial(pole["polynomial"], symbols) + 1)
        try:
            prepare_vector(wrong, mass, "1")
        except AssertionError:
            controls.append(f"missing_gluon_mass_edge_{edge}")
        else:
            raise AssertionError("a missing matched gluon regulator was admitted")

    # Delta=x^2+(1-x)*lambda^2 at p^2=m^2=1. H=exp(gamma*eps)*Gamma(eps).
    # J1=H*int (1-x)*Delta^-eps; J2=-eps*H*int x*(1-x)*Delta^(-1-eps).
    # The cubic Minkowski denominator supplies the minus sign in J2.
    # f_lambda=H*int [2+(2-2eps)*x]*Delta^-eps is the OS mass coefficient.
    # CT_child=(kslash-m)/eps-f_lambda*m, with fields in MS-bar.
    # P_lambda=(2-d)*(J1/eps-2*f_lambda*J2), d=4-2eps.
    # J1 through eps^1, f_lambda through eps^0 and J2 through eps^1 suffice;
    # retain f_lambda's eps^1 coefficient too as an enclosing-operation witness.
    with mp.workdps(80):
        delta = lambda x: x * x - x + 1
        integrands = {
            "L1": lambda x: (1 - x) * mp.log(delta(x)),
            "L2": lambda x: (1 - x) * mp.log(delta(x)) ** 2,
            "R0": lambda x: x * (1 - x) / delta(x),
            "R1": lambda x: x * (1 - x) * mp.log(delta(x)) / delta(x),
            "I0": lambda x: (2 + 2 * x) * mp.log(delta(x)),
            "f1_log": lambda x: (
                (1 + x) * mp.log(delta(x)) ** 2 + 2 * x * mp.log(delta(x))
            ),
        }
        moments = {name: mp.quad(fn, [0, 1]) for name, fn in integrands.items()}
        l1, l2, r0, r1, i0 = (moments[name] for name in ("L1", "L2", "R0", "R1", "I0"))
        assert abs(l1 - (-1 + mp.pi / (2 * mp.sqrt(3)))) < mp.mpf("1e-70")
        assert abs(r0 - (-1 + 2 * mp.pi / (3 * mp.sqrt(3)))) < mp.mpf("1e-70")
        assert abs(i0 - (-6 + mp.sqrt(3) * mp.pi)) < mp.mpf("1e-70")
        finite = -l2 - mp.zeta(2) / 2 + 12 * r1 + 4 * i0 * r0 - 2 * l1 + 16 * r0
        proper = {
            -2: mp.mpf(-1),
            -1: 11 - 7 * mp.pi / mp.sqrt(3),
            0: finite,
        }
        moment_report = {name: mp.nstr(value, 75) for name, value in moments.items()}
        moment_report["f1"] = mp.nstr(moments["f1_log"] + 3 * mp.zeta(2) / 2, 75)
        proper_report = {
            str(power): mp.nstr(value, 75) for power, value in proper.items()
        }
        epsilon_checks = []
        for epsilon in (mp.mpf("1e-12"), mp.mpf("2e-12")):
            h = mp.exp(mp.euler * epsilon) * mp.gamma(epsilon)
            j1 = h * mp.quad(
                lambda x, epsilon=epsilon: (1 - x) * delta(x) ** (-epsilon), [0, 1]
            )
            j2 = (
                -epsilon
                * h
                * mp.quad(
                    lambda x, epsilon=epsilon: x * (1 - x) * delta(x) ** (-1 - epsilon),
                    [0, 1],
                )
            )
            f = h * mp.quad(
                lambda x, epsilon=epsilon: (
                    (2 + (2 - 2 * epsilon) * x) * delta(x) ** (-epsilon)
                ),
                [0, 1],
            )
            direct = (-2 + 2 * epsilon) * (j1 / epsilon - 2 * f * j2)
            remainder = direct - sum(
                value * epsilon**power for power, value in proper.items()
            )
            assert abs(remainder) < 100 * epsilon
            epsilon_checks.append(
                {
                    "epsilon": mp.nstr(epsilon, 25),
                    "finite_subtracted_remainder": mp.nstr(remainder, 60),
                    "absolute_bound": mp.nstr(100 * epsilon, 25),
                }
            )
    inputs = backend_input(regulated, regulated["rows"]["CF2"])
    emitted_source_sha256 = None
    if emitted_regulated_source is not None:
        emitted = prepare_vector(emitted_regulated_source, mass, "1")
        assert backend_input(emitted, emitted["rows"]["CF2"]) == inputs
        emitted_source_sha256 = hashlib.sha256(
            json.dumps(
                emitted_regulated_source, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest()
    emitted_unregulated_sha256 = None
    emitted_unregulated_input_sha256 = None
    if emitted_unregulated_source is not None:
        original = prepare_vector(vector, mass, "0")
        emitted = prepare_vector(emitted_unregulated_source, mass, "0")
        emitted_inputs = backend_input(emitted, emitted["rows"]["CF2"])
        assert emitted_inputs == backend_input(original, original["rows"]["CF2"])
        emitted_unregulated_input_sha256 = hashlib.sha256(
            json.dumps(emitted_inputs, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        assert emitted_unregulated_input_sha256 == D09_VECTOR_BACKEND_INPUT_SHA256
        emitted_unregulated_sha256 = hashlib.sha256(
            json.dumps(
                emitted_unregulated_source, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest()
    output.mkdir(parents=True)
    validation = {
        "mode": "matched_positive_gluon_mass_D09_vector_only",
        "input_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "normalization": "a^2*CF^2, amputated Sigma, p^2=m^2=mu_bar^2=lambda^2=1",
        "regulated_source": payload,
        "backend_inputs": inputs,
        "emitted_regulated_source_bound": emitted_source_sha256 is not None,
        "emitted_regulated_source_sha256": emitted_source_sha256,
        "input_binding_status": (
            "BOUND_CURRENT_EMITTED_TYPED_SOURCE"
            if emitted_source_sha256 is not None
            else "REFERENCE_PENDING_CURRENT_EMITTED_INPUT_BINDING"
        ),
        "emitted_unregulated_source_bound": emitted_unregulated_sha256 is not None,
        "emitted_unregulated_source_sha256": emitted_unregulated_sha256,
        "emitted_unregulated_backend_input_sha256": emitted_unregulated_input_sha256,
        "unregulated_vector_reference_kind": "retained_pinned_backend_numerical_regression",
        "unregulated_vector_reference": D09_VECTOR_NUMERICAL_REFERENCE,
        "protocol_negative_controls": controls,
        "proper_moments": moment_report,
        "proper_coefficients": proper_report,
        "proper_direct_dimensional_parameter_checks": epsilon_checks,
        "proper_finite_formula": "-L2-zeta(2)/2+12*R1+4*I0*R0-2*L1+16*R0",
        "signed_parent_MS_poles": {"-2": "1/2", "-1": "-1/4", "0": "0"},
        "policy": POLICY,
        "integration_performed": False,
        "production_OS_forest_integrated_condition_evaluated": False,
    }
    (output / "validation.json").write_text(json.dumps(validation, indent=2))
    if validate_only:
        return
    multiprocessing.set_start_method("fork")
    rows, report = integrate("D09-vector-massive-gluon", inputs, output / "bare", 0)
    prepared = combine(
        [
            (1, rows),
            (1, {power: (complex(value), 0j) for power, value in proper.items()}),
        ]
    )
    complete = combine([(1, prepared), (1, {-2: (0.5 + 0j, 0j), -1: (-0.25 + 0j, 0j)})])
    checks = []
    compare(
        "D09_regulated_vector_signed_UV_pole_residual",
        {power: value for power, value in complete.items() if power < 0},
        {},
        checks,
    )
    summary = {
        **validation,
        "integration_performed": True,
        "backend_case": report,
        "checks": checks,
        "finite_vector_reference": {
            "value": [complete[0][0].real, complete[0][0].imag],
            "bare_reported_absolute_error": [complete[0][1].real, complete[0][1].imag],
        },
        "passed": all(row["passed"] for row in checks),
        "unregulated_limit_or_CFF_exchange_proved": False,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    assert summary["passed"], "matched-regulator reference UV poles do not cancel"


def vector_reference(rows, checks):
    # Independent proper insertion. Keep the mass and vector masters through
    # epsilon^1, and substitute d=4-2eps only after the full ordered trace.
    h = sp.exp(sp.EulerGamma * EPS) * sp.gamma(EPS)
    j1 = h / ((1 - 2 * EPS) * (2 - 2 * EPS))
    m2_j2 = h / (2 * (1 - 2 * EPS))
    f = h * (3 - 2 * EPS) / (1 - 2 * EPS)
    proper = {-2: sp.Integer(5), -1: sp.Integer(12), 0: 32 + 11 * sp.pi**2 / 12}
    exact = sp.series(
        (2 - (4 - 2 * EPS)) * (j1 / EPS - 2 * f * m2_j2), EPS, 0, 1
    ).removeO()
    assert (
        sp.simplify(exact - sum(value * EPS**power for power, value in proper.items()))
        == 0
    )
    prepared = combine(
        [
            (1, rows),
            (1, {power: (complex(value), 0j) for power, value in proper.items()}),
        ]
    )
    # Source-derived amputated convention: W=-(U-UO), i*Omega=+a^2*CF^2.
    # External ordering -1 stays outside. Signed-pole/source receipt:
    # 459298b98222213e8f615ab8212f0def89d43d9027e6885e090bfc6573da904e
    # Numeric source/phase receipt is retained above (615bbd0e...).
    uv = {-2: (0.5 + 0j, 0j), -1: (-0.25 + 0j, 0j)}
    residual = {
        power: value
        for power, value in combine([(1, prepared), (1, uv)]).items()
        if power < 0
    }
    # Each pole residual has reference zero: the unchanged bound is exactly
    # 1e-6 per real/imaginary component. Extra negative powers are checked too.
    compare("D09_vector_signed_UV_pole_residual", residual, {}, checks)
    compare(
        "D09_vector_finite_numerical_regression",
        {0: prepared[0]},
        {0: float(D09_VECTOR_NUMERICAL_REFERENCE["finite_prepared"])},
        checks,
    )
    finite = checks[-1]
    finite["reference_kind"] = "retained_pinned_backend_numerical_regression"
    finite["new_bare_reported_absolute_error"] = finite[
        "sum_reported_absolute_errors"
    ].copy()
    finite["reference_reported_absolute_error"] = [
        float(D09_VECTOR_NUMERICAL_REFERENCE["reported_absolute_error"]),
        0.0,
    ]
    finite["linear_sum_reported_absolute_errors"] = [
        actual + reference
        for actual, reference in zip(
            finite["new_bare_reported_absolute_error"],
            finite["reference_reported_absolute_error"],
            strict=True,
        )
    ]
    finite["uncertainty_policy"] = (
        "Linear absolute error sum; no Gaussian independence claim; comparison center bound unchanged"
    )
    return {
        "proper_exact": {str(power): sp.sstr(value) for power, value in proper.items()},
        "prepared": {
            str(power): {
                "value": [value.real, value.imag],
                "bare_error": [error.real, error.imag],
            }
            for power, (value, error) in sorted(prepared.items())
        },
        "finite_reference": D09_VECTOR_NUMERICAL_REFERENCE,
        "passed": all(
            row["passed"] for row in checks if row["scope"].startswith("D09_vector_")
        ),
        "ordinary_4D_IR_integrability_proved": False,
        "simultaneous_CFF_exchange_proved": False,
        "production_R2_term_added": False,
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
    parser.add_argument("--regulated-vector-only", action="store_true")
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
    vector = None
    vector_controls = []
    if "vector_source" in data:
        vector = prepare_vector(data["vector_source"], prepared["D09"])
        vector_controls = vector_protocol_controls(
            data["vector_source"], prepared["D09"]
        )
    light = [source for source in sources if source["light_flavor"] is not None]
    assert len(light) == 5 and {source["light_flavor"] for source in light} == set(
        range(1, 6)
    )
    assert {source["name"] for source in light} == {f"D0{i}" for i in range(4, 9)}
    for source in light[1:]:
        assert {
            k: v for k, v in source.items() if k not in ("name", "light_flavor")
        } == {k: v for k, v in light[0].items() if k not in ("name", "light_flavor")}
    if args.regulated_vector_only:
        assert vector is not None, "the native vector source is required"
        regulated_vector_reference(
            data["vector_source"],
            prepared["D09"],
            source_path,
            output,
            args.validate_only,
            data.get("emitted_regulated_source"),
            data.get("emitted_unregulated_source"),
        )
        return
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
        "optional_numeric_cases": ["D09-vector"] if vector is not None else [],
        "optional_vector_source_checked": vector is not None,
        "optional_vector_protocol_controls": vector_controls,
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
    vector_report = None
    vector_result = None
    if vector is not None:
        vector_rows, vector_report = integrate(
            "D09-vector",
            backend_input(vector, vector["rows"]["CF2"]),
            output / "D09-vector",
            0,
        )
        vector_result = vector_reference(vector_rows, checks)
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
        "optional_vector_backend_case": vector_report,
        "optional_vector_reference": vector_result,
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
