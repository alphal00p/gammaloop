"""Regenerate formal master-product coefficients from a saved q0 result.

The caller supplies its existing Symbolica API and canonical Vakint engine type.
No reduction, external-data loading, or master evaluation is performed. Masters
and their products remain formal; all Laurent remainders are coefficientwise.
"""

import hashlib
import json
import os
import re
import tempfile
import time
from fractions import Fraction
from pathlib import Path

from .artifacts import reconstruct_coefficient_table, require_raw_equal
from .series import compute_scalar_jet


class Q0Generator:
    """Produce exact and tracked 300-digit coefficients through epsilon**5.

    ``api`` is the already imported Symbolica module, and ``engine_type`` is the
    Vakint class belonging to that same native host. This object never imports a
    second Symbolica kernel. The optional independent oracle is a development
    check; mpmath is neither imported nor required during normal regeneration.
    """

    decimal_precision = 300
    guard_decimal_digits = 32
    cutoff = 5

    def __init__(self, api, engine_type):
        self.api = api
        self.engine_type = engine_type

    @staticmethod
    def _require(condition, message):
        if not condition:
            raise ValueError(message)

    @staticmethod
    def _sha(path):
        digest = hashlib.sha256()
        with Path(path).open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    @staticmethod
    def _write_json(path, value):
        with path.open("x", encoding="utf-8") as output:
            json.dump(value, output, indent=2, allow_nan=False)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())

    @staticmethod
    def _ratio_record(value):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }

    def _component_record(self, raw, precision, rounded, label):
        center = Fraction(*rounded.as_integer_ratio())
        radius = Fraction(*precision.absolute_error.as_integer_ratio()) + abs(
            center - Fraction(*raw.as_integer_ratio())
        )
        if precision.exact_zero:
            self._require(
                center == 0 and radius == 0,
                f"exact-zero claim changed during export: {label}",
            )
            achieved = None
        else:
            self._require(
                precision.significant_digits is not None
                and precision.significant_digits >= self.decimal_precision,
                f"tracked coefficient achieved fewer than 300 decimal digits: {label}",
            )
            self._require(
                abs(center) > radius
                and radius * 10**self.decimal_precision <= abs(center) - radius,
                f"final binary rounding does not establish 300 relative decimal digits: {label}",
            )
            achieved = self.decimal_precision
        return {
            "center_exact_ratio": self._ratio_record(center),
            "display_300_significant_digits": str(
                rounded.to_decimal(digits=self.decimal_precision)
            ),
            "display_is_not_additional_accuracy": True,
            "absolute_error_exact_ratio_including_export_rounding": self._ratio_record(
                radius
            ),
            "tracked_absolute_error": str(precision.absolute_error),
            "tracked_significant_digits": precision.significant_digits,
            "exported_reliable_decimal_digits_at_least": achieved,
            "exact_zero": precision.exact_zero,
            "error_quality": str(precision.quality),
            "float_binary_precision": rounded.precision,
        }

    def _integer_monomial_power(self, monomial, symbol):
        if monomial == self.api.E("1"):
            return 0
        if monomial == symbol:
            return 1
        self._require(
            monomial.get_type() == self.api.AtomType.Pow,
            "non-Laurent monomial in independent oracle",
        )
        base, exponent = list(monomial)
        self._require(base == symbol, "unexpected monomial base")
        text = exponent.to_canonical_string()
        self._require(
            re.fullmatch(r"-?[0-9]+", text) is not None,
            "noninteger monomial exponent",
        )
        return int(text)

    @staticmethod
    def _oracle_context(enabled):
        if not enabled:
            return None
        try:
            import mpmath
        except ImportError as error:
            raise ValueError(
                "independent validation was requested but optional mpmath is unavailable"
            ) from error
        context = mpmath.mp.clone()
        context.dps = 380
        return context, mpmath.__version__

    def _validate_oracle(self, oracle, verified_jets, rows, epsilon, gamma):
        if oracle is None:
            return {
                "status": "not_requested",
                "required_end_user_dependency": False,
            }
        context, version = oracle
        checked = 0
        for row in rows:
            key = row["master_product"]
            numbers = {value["epsilon_power"]: value for value in row["coefficients"]}
            for monomial, exact in verified_jets[key].coefficient_list(epsilon):
                power = self._integer_monomial_power(monomial, epsilon)
                expected = context.mpf(0)
                for monomial, rational in exact.coefficient_list(gamma):
                    exponent = self._integer_monomial_power(monomial, gamma)
                    rational = Fraction(rational.to_canonical_string())
                    expected += (
                        context.mpf(rational.numerator)
                        / rational.denominator
                        * context.euler**exponent
                    )
                number = numbers[power]
                for name, expected_component in (
                    ("real", expected),
                    ("imaginary", context.mpf(0)),
                ):
                    record = number[name]
                    center = record["center_exact_ratio"]
                    bound = record[
                        "absolute_error_exact_ratio_including_export_rounding"
                    ]
                    actual = context.mpf(center["numerator"]) / context.mpf(
                        center["denominator"]
                    )
                    radius = context.mpf(bound["numerator"]) / context.mpf(
                        bound["denominator"]
                    )
                    self._require(
                        abs(expected_component - actual) <= radius,
                        f"independent coefficient oracle outside enclosure: "
                        f"{key}, epsilon^{power}, {name}",
                    )
                checked += 1
        return {
            "status": "passed",
            "implementation": "mpmath",
            "version": version,
            "fixed_decimal_precision": 380,
            "real_coefficients_checked": checked,
            "imaginary_coefficients_checked": checked,
            "same_Vakint_formula_not_retried": True,
            "required_end_user_dependency": False,
        }

    def generate(self, input_path, output_path, validate_oracle=False):
        """Create a fresh artifact directory and return its path.

        ``input_path`` contains the saved retained expression. Each
        combined scalar coefficient is evaluated once at 300 requested digits
        with Vakint's fixed 32-digit guard allowance. Accuracy failures are
        errors, never retries. A failed run preserves its sibling staging
        directory for inspection and never publishes a partial final directory.
        """
        input_path, output_path = Path(input_path), Path(output_path)
        self._require(not output_path.exists(), "artifact output must be fresh")
        oracle = self._oracle_context(validate_oracle)
        started = time.monotonic()
        payload = input_path.read_bytes()
        text = payload.decode("utf-8")
        api, expression = self.api, self.api.E
        retained = expression(text)
        calls = re.findall(r"(vakint::\{\}::RustRedMaster\s*\([^()]*\))", text)
        self._require(
            calls and len(calls) == len(re.findall(r"\bRustRedMaster\s*\(", text)),
            "unrecognized retained master syntax or no formal masters",
        )
        masters_by_key = {}
        for call in calls:
            master = expression(call)
            masters_by_key[master.to_canonical_string()] = master
        masters = [masters_by_key[key] for key in sorted(masters_by_key)]
        epsilon, gamma = expression("vakint::ε"), expression("symbolica::γ")
        coupling = expression("UFO::GC_11")
        mass, scale = expression("vakint::mUV"), expression("vakint::mursq")
        normalization = expression("𝑒") ** (-5 * gamma * epsilon)
        unit = retained.replace(mass, 1).replace(scale, 1)
        scale_present = scale in retained.get_all_symbols(False)
        grouping = {}
        for monomial, coefficient in unit.coefficient_list(*masters, coupling):
            if coefficient == expression("0"):
                continue
            common = monomial.coefficient_list(coupling)
            self._require(
                len(common) == 1 and common[0][0] == coupling**10,
                "coupling is not homogeneous GC_11^10",
            )
            product = common[0][1]
            key = product.to_canonical_string()
            self._require(
                key not in grouping,
                "duplicate canonical product after coupling extraction",
            )
            grouping[key] = product, coefficient
        reconstructed = (
            sum(
                (product * coefficient for product, coefficient in grouping.values()),
                expression("0"),
            )
            * coupling**10
        )
        for _, coefficient in (unit - reconstructed).coefficient_list(
            *masters, coupling
        ):
            self._require(
                coefficient.expand() == expression("0"),
                "exact master-polynomial reconstruction failed",
            )
        allowed = {epsilon, gamma, expression("𝑒")}
        for key, (_, coefficient) in grouping.items():
            unknown = [
                symbol.to_canonical_string()
                for symbol in coefficient.get_all_symbols(False)
                if symbol not in allowed
            ]
            self._require(
                not unknown,
                f"unresolved coefficient parameters for {key}: {unknown}",
            )

        engine = self.engine_type(
            run_time_decimal_precision=self.decimal_precision,
            evaluation_order=[],
            epsilon_symbol=epsilon,
            mu_r_sq_symbol=scale,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(
            tempfile.mkdtemp(
                prefix=f".{output_path.name}.pending-", dir=output_path.parent
            )
        )
        rows, verified_jets = [], {}
        symbolic_truncated, pseudo = expression("0"), expression("0")
        # The same integer-only decimal-to-binary rule as checked Vakint export.
        bits = (
            self.decimal_precision * 3_321_928_094_887_362_348 + 10**18 - 1
        ) // 10**18 + 1
        for key in sorted(grouping):
            product, coefficient = grouping[key]
            exact_poly, exact_check = compute_scalar_jet(
                api, coefficient, epsilon, gamma, cutoff=self.cutoff
            )
            self._require(
                exact_check["status"] == "PASS"
                and exact_check["remainder_order_at_least"] >= self.cutoff + 1,
                "factorwise exact remainder proof failed",
            )
            verified_jets[key] = exact_poly
            numerical = engine.numerical_result_from_expression(exact_poly)
            self._require(
                numerical.requested_decimal_precision == self.decimal_precision,
                "tracked evaluation changed the requested precision",
            )
            numeric_poly = numerical.to_expression(epsilon)
            exported = []
            for power, item in sorted(numerical.coefficients.items()):
                value = item.numerical_value
                # Exact zero's storage precision is part of the native value.
                # Capping it to 998 bits would silently change raw Atom bytes.
                rounded = api.ComplexFloat(
                    value.real
                    if item.precision.real.exact_zero
                    else value.real.with_precision(precision=bits),
                    value.imag
                    if item.precision.imaginary.exact_zero
                    else value.imag.with_precision(precision=bits),
                )
                checked = item.to_expression()
                require_raw_equal(
                    api.Expression.num(rounded),
                    checked,
                    f"exported Float values and storage precisions: {key}, epsilon^{power}",
                )
                exported.append(
                    {
                        "epsilon_power": power,
                        "native_value_expression": checked.to_canonical_string(),
                        "real": self._component_record(
                            value.real,
                            item.precision.real,
                            rounded.real,
                            f"{key}, epsilon^{power}, real",
                        ),
                        "imaginary": self._component_record(
                            value.imag,
                            item.precision.imaginary,
                            rounded.imag,
                            f"{key}, epsilon^{power}, imaginary",
                        ),
                    }
                )
            warnings = [
                {
                    "epsilon_power": warning.epsilon_power,
                    "message": str(warning),
                    "master_product": key,
                }
                for warning in numerical.warnings
            ]
            rows.append(
                {
                    "master_product": key,
                    "normalization": normalization.to_canonical_string(),
                    "exact_laurent": exact_poly.to_canonical_string(),
                    "remainder": {
                        "coefficientwise": True,
                        "epsilon_power": self.cutoff + 1,
                    },
                    "exact_recurrence_check": exact_check,
                    "requested_decimal_precision": self.decimal_precision,
                    "fixed_guard_decimal_digits": self.guard_decimal_digits,
                    "automatic_retries": 0,
                    "coefficients": exported,
                    "warnings": warnings,
                }
            )
            symbolic_truncated += product * exact_poly
            pseudo += product * numeric_poly
        symbolic_truncated *= coupling**10
        pseudo *= coupling**10
        names = (
            ("q0-retained", retained),
            ("q0-coefficients-exact-through-epsilon5", symbolic_truncated),
            ("q0-coefficients-300digits-through-epsilon5", pseudo),
        )
        for name, value in names:
            native_path = stage / (name + ".atom")
            value.save(str(native_path))
            loaded = api.Expression.load(str(native_path))
            require_raw_equal(loaded, value, f"native Atom/state roundtrip: {name}")
            portable = stage / (name + ".expr")
            if name == "q0-retained":
                portable.write_bytes(payload)
            else:
                portable.write_text(value.to_canonical_string(), encoding="utf-8")
        table = {
            "schema": "vakint.q0-formal-master-product-coefficients.v1",
            "common_coupling_prefactor": (coupling**10).to_canonical_string(),
            "unit_mass": {
                "mUV": 1,
                "mu_r": 1,
                "mu_r_squared_symbol": scale.to_canonical_string(),
                "mu_r_symbol_present_in_saved_input": scale_present,
            },
            "dimension_substitution_already_present": True,
            "normalization_applied_again": False,
            "unknown_master_values_evaluated": False,
            "numerical_expr_format": (
                "human decimal display only; native Atom and exact-ratio JSON are authoritative"
            ),
            "cutoff": {
                "included_maximum_power": self.cutoff,
                "remainder_starts_at": self.cutoff + 1,
                "scope": (
                    "each formal master-PRODUCT scalar coefficient; "
                    "not the fully evaluated physical Laurent series"
                ),
            },
            "products": rows,
        }
        self._write_json(stage / "coefficient-table.json", table)
        require_raw_equal(
            reconstruct_coefficient_table(table, api=api),
            pseudo,
            "lossless JSON reconstruction of the whole numerical expression",
        )
        oracle_result = self._validate_oracle(
            oracle, verified_jets, rows, epsilon, gamma
        )
        source_paths = (
            Path(__file__),
            Path(__file__).with_name("series.py"),
            Path(__file__).with_name("artifacts.py"),
        )
        metadata = {
            "schema": "vakint.q0-coefficient-artifact-provenance.v1",
            "input": {
                "path": str(input_path.resolve()),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            },
            "regeneration": {
                "symbolica_api": getattr(api, "__name__", None),
                "symbolica_version": getattr(api, "__version__", None),
                "vakint_engine": (
                    f"{self.engine_type.__module__}.{self.engine_type.__qualname__}"
                ),
                "source_sha256": {path.name: self._sha(path) for path in source_paths},
            },
            "literal_master_call_count": len(calls),
            "master_product_count": len(rows),
            "individual_master_count": len(masters),
            "numerical_coefficient_count": sum(
                len(row["coefficients"]) for row in rows
            ),
            "exact_reconstruction": True,
            "exact_recurrence_every_product": True,
            "exact_algorithm": (
                "pole-stripped polynomial-unit finite jets; "
                "Fraction recurrence and independent native factorwise Series"
            ),
            "whole_scalar_series_or_rational_GCD_used": False,
            "lossless_JSON_native_export_equal": True,
            "human_display_only_artifacts": [
                "q0-coefficients-300digits-through-epsilon5.expr"
            ],
            "tracked_requested_decimal_digits": self.decimal_precision,
            "fixed_guard_decimal_digits": self.guard_decimal_digits,
            "automatic_retries": 0,
            "independent_oracle": oracle_result,
            "artifacts": {
                path.name: {"bytes": path.stat().st_size, "sha256": self._sha(path)}
                for path in stage.iterdir()
                if path.is_file()
            },
            "no_IBP_or_master_evaluation_or_cache_io": True,
            "elapsed_seconds": time.monotonic() - started,
        }
        self._require(
            all(
                value["bytes"] < 100 * 1024 * 1024
                for value in metadata["artifacts"].values()
            ),
            "artifact exceeds the repository 100 MiB file limit; no truncation performed",
        )
        self._write_json(stage / "provenance.json", metadata)
        self._require(
            not output_path.exists(), "final artifact directory appeared concurrently"
        )
        os.rename(stage, output_path)
        return output_path
