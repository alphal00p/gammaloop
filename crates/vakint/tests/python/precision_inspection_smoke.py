"""Exercise the installed, coherent Symbolica/Vakint extension without FORM."""

import copy
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    import vakint_hepkit_standalone as core
    from vakint_hepkit_standalone import (
        ComplexFloat,
        ConstantSymbol,
        E,
        FeynmanDiagram,
        Float,
        MasterSymbol,
        Vakint,
    )
else:
    import symbolica as core
    from symbolica import ComplexFloat, E, Float
    from symbolica.community.hepkit import FeynmanDiagram
    from symbolica.community.hepkit.vakint import ConstantSymbol, MasterSymbol, Vakint


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples"))
try:
    from q0_coefficients import CoefficientArtifactError, SavedQ0Result
    from q0_coefficients.artifacts import (
        reconstruct_coefficient_table,
        require_raw_equal,
    )
    from q0_coefficients.generate import Q0Generator
    from q0_coefficients.series import probe as q0_series_probe
finally:
    sys.path.pop(0)


class PrecisionInspectionSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = Vakint(evaluation_order=[], run_time_decimal_precision=32)

    def test_resource_settings_are_explicit_and_disabled_by_default(self):
        for seconds in (float("nan"), float("inf"), -1.0):
            with (
                self.subTest(seconds=seconds),
                self.assertRaisesRegex(ValueError, "finite and nonnegative"),
            ):
                Vakint(evaluation_order=[], max_elapsed_seconds=seconds)
        with self.assertRaisesRegex(ValueError, "elapsed time"):
            Vakint(evaluation_order=[], max_elapsed_seconds=0)
        with self.assertRaisesRegex(ValueError, "process RSS"):
            Vakint(evaluation_order=[], max_process_rss_bytes=1)
        ordinary = Vakint(
            evaluation_order=[],
            max_process_rss_bytes=None,
            max_elapsed_seconds=None,
        )
        record = ordinary.get_constant(ConstantSymbol("z3"), prec=32)
        self.assertGreaterEqual(record.precision.real.significant_digits, 32)

    def test_constant_record_preserves_native_precision_and_is_read_only(self):
        constant = self.engine.get_constant(ConstantSymbol("z3"), prec=32)
        self.assertEqual(constant.available_precision.kind, "ArbitraryPrecision")
        self.assertEqual(constant.requested_decimal_precision, 32)
        self.assertIsNone(constant.available_precision.max_decimal_precision)
        self.assertIsInstance(constant.numerical_value, Float)
        self.assertIsInstance(constant.precision.real.absolute_error, Float)
        self.assertGreaterEqual(constant.precision.real.significant_digits, 32)
        self.assertTrue(constant.precision.imaginary.exact_zero)
        self.assertTrue(constant.provenance.version)
        self.assertIsInstance(constant.definition, core.Expression)
        self.assertIn("Integral(", str(constant.definition))
        with self.assertRaises(AttributeError):
            constant.definition = constant.definition
        self.assertIn("Vakint constant", repr(constant))
        self.assertIn("Available precision", str(constant))
        with self.assertRaises(AttributeError):
            constant.numerical_value = E("0")
        with self.assertRaises(AttributeError):
            constant.precision.real.significant_digits = 1000
        with patch.dict(os.environ, {"NO_COLOR": "1"}):
            self.assertNotIn("\x1b[", repr(constant))

    def test_master_keeps_all_orders_and_native_graph(self):
        symbols = self.engine.master_symbols()
        self.assertTrue(symbols)
        master = self.engine.get_master(symbols[0], prec=24)
        self.assertEqual(master.requested_decimal_precision, 24)
        self.assertTrue(master.is_computed)
        self.assertIsInstance(master.graph, FeynmanDiagram)
        self.assertTrue(master.coefficients)
        self.assertEqual(list(master.coefficients), sorted(master.coefficients))
        for power, coefficient in master.coefficients.items():
            self.assertEqual(coefficient.epsilon_power, power)
            self.assertIsInstance(coefficient.numerical_value, Float)
            self.assertIn(
                coefficient.symbolic_expression.kind, ("Expression", "NumericalOnly")
            )
        self.assertIn("digraph", master.graph.to_dot())
        self.assertIn("Vakint master", repr(master))
        with self.assertRaises(AttributeError):
            master.coefficients = {}

    def test_notebook_display_uses_stored_records(self):
        constant = self.engine.get_constant(ConstantSymbol("z3"), prec=24)
        html = constant.to_html()
        self.assertEqual(html, constant._repr_html_())
        self.assertIn('class="vakint-card"', html)
        self.assertIn("Available precision", html)
        master = self.engine.get_master(self.engine.master_symbols()[0], prec=24)
        html = master.to_html()
        self.assertIn("<table>", html)
        self.assertIn("<svg", html)

    def test_unknown_and_invalid_identities_fail_explicitly(self):
        with self.assertRaisesRegex(ValueError, "unknown RustRed constant"):
            self.engine.get_constant(ConstantSymbol("missing_constant"))
        with self.assertRaisesRegex(ValueError, "unknown RustRed master"):
            self.engine.get_master(MasterSymbol("missing_master"))
        with self.assertRaisesRegex(ValueError, "ambiguous RustRed master"):
            self.engine.get_master(MasterSymbol("P1"))
        with self.assertRaises(ValueError):
            ConstantSymbol(E("2"))

    def test_uncertain_zero_remains_inspectable_but_cannot_become_expression(self):
        result, statistical_error = self.engine.numerical_evaluation(
            E("x-y", default_namespace="vakint"), {"x": 1.0, "y": 1.0}
        )
        self.assertIsNone(statistical_error)
        self.assertEqual(result.requested_decimal_precision, 32)
        self.assertIsInstance(result.coefficients[0].numerical_value, ComplexFloat)
        self.assertFalse(result.coefficients[0].precision.real.exact_zero)
        self.assertIsNone(result.coefficients[0].precision.real.significant_digits)
        self.assertTrue(result.warnings)
        with self.assertRaises(ValueError):
            result.to_expression()
        with self.assertRaises(ValueError):
            self.engine.numerical_result_to_expression(result)
        self.assertIn("Vakint numerical result", result.to_html())
        self.assertIn("Precision warning", result.to_html())

    def test_exact_cancellation_converts_normally(self):
        result, error = self.engine.numerical_evaluation(
            E("x-x", default_namespace="vakint"), {}
        )
        self.assertIsNone(error)
        self.assertEqual(result.to_expression(), E("0"))

    def test_inexact_expression_inputs_keep_independent_source_identities(self):
        for value in (
            E("0.123456789`20"),
            E("z3+0.123456789`20", default_namespace="vakint"),
        ):
            with self.subTest(value=str(value)):
                result, error = self.engine.numerical_evaluation(
                    E("a-b", default_namespace="vakint"), {"a": value, "b": value}
                )
                self.assertIsNone(error)
                self.assertFalse(result.coefficients[0].precision.real.exact_zero)
                self.assertTrue(result.warnings)
                with self.assertRaises(ValueError):
                    result.to_expression()
                exact, _ = self.engine.numerical_evaluation(
                    E("a-a", default_namespace="vakint"), {"a": value}
                )
                self.assertEqual(exact.to_expression(), E("0"))
        result, _ = self.engine.numerical_evaluation(
            E("dot(p(1),p(1))", default_namespace="vakint"),
            {},
            externals={1: (E("1`20"), E("1`20"), 0, 0)},
        )
        self.assertFalse(result.coefficients[0].precision.real.exact_zero)
        with self.assertRaises(ValueError):
            result.to_expression()

    def test_native_float_parameters_preserve_more_than_binary64(self):
        engine = Vakint(evaluation_order=[], run_time_decimal_precision=80)
        result, error = engine.numerical_evaluation(
            E("x+z", default_namespace="vakint"),
            {
                "x": Float("1.234567890123456789", decimal_digits=200),
                "z": ComplexFloat(
                    "2.345678901234567891", "3.456789012345678912", decimal_digits=200
                ),
            },
        )
        self.assertIsNone(error)
        self.assertGreaterEqual(
            result.coefficients[0].precision.real.significant_digits, 80
        )
        self.assertGreaterEqual(
            result.coefficients[0].precision.imaginary.significant_digits, 80
        )
        self.assertFalse(result.warnings)
        limited, _ = engine.numerical_evaluation(
            E("x", default_namespace="vakint"), {"x": 1.2345678901234567}
        )
        self.assertLessEqual(
            limited.coefficients[0].precision.real.significant_digits, 15
        )
        self.assertTrue(limited.warnings)

    def test_exact_rational_parameter_and_near_lightlike_external(self):
        result, error = self.engine.numerical_evaluation(
            E("muvsq^-2*(1+ε^-1)", default_namespace="vakint"),
            {"muvsq": E("3/2")},
        )
        self.assertIsNone(error)
        reference = self.engine.numerical_result_from_expression(
            E("4/9*(1+ε^-1)", default_namespace="vakint")
        )
        self.assertEqual(result.to_expression(), reference.to_expression())
        result, error = self.engine.numerical_evaluation(
            E("dot(p(1),p(1))", default_namespace="vakint"),
            {},
            externals={1: (10**80 + 1, 10**80, 0, 0)},
        )
        self.assertIsNone(error)
        reference = self.engine.numerical_result_from_expression(E(str(2 * 10**80 + 1)))
        self.assertEqual(result.to_expression(), reference.to_expression())
        self.assertGreaterEqual(
            result.coefficients[0].precision.real.significant_digits, 32
        )

    def test_complex_input_components_keep_independent_accuracy(self):
        engine = Vakint(evaluation_order=[], run_time_decimal_precision=80)
        value = ComplexFloat(
            Float("1.23456789", decimal_digits=20),
            Float("2.34567891", decimal_digits=200),
        )
        result, error = engine.numerical_evaluation(
            E("z", default_namespace="vakint"), {"z": value}
        )
        self.assertIsNone(error)
        self.assertLessEqual(
            result.coefficients[0].precision.real.significant_digits, 20
        )
        self.assertGreaterEqual(
            result.coefficients[0].precision.imaginary.significant_digits, 80
        )

    def test_native_external_components_keep_their_precision(self):
        engine = Vakint(evaluation_order=[], run_time_decimal_precision=80)
        result, error = engine.numerical_evaluation(
            E("dot(p(1),p(1))", default_namespace="vakint"),
            {},
            externals={
                1: (
                    Float("1.000000000000000000000000000001", decimal_digits=200),
                    Float("1", decimal_digits=200),
                    0,
                    0,
                )
            },
        )
        self.assertIsNone(error)
        self.assertGreaterEqual(
            result.coefficients[0].precision.real.significant_digits, 80
        )
        self.assertFalse(result.warnings)

    def test_numerical_input_validation_is_explicit(self):
        expression = E("x", default_namespace="vakint")
        with self.assertRaises(TypeError):
            self.engine.numerical_evaluation(expression, {"x": True})
        with self.assertRaisesRegex(ValueError, "duplicate canonical parameter"):
            self.engine.numerical_evaluation(expression, {"x": 1, "vakint::x": 2})
        with self.assertRaisesRegex(ValueError, "four components"):
            self.engine.numerical_evaluation(expression, {}, externals={1: (1, 0, 0)})

    def test_parameter_maps_cannot_override_registered_constants(self):
        expression = E("z3", default_namespace="vakint")
        for value in (E("2"), Float("2", decimal_digits=80)):
            with (
                self.subTest(value=str(value)),
                self.assertRaisesRegex(ValueError, "registered constant"),
            ):
                self.engine.numerical_evaluation(expression, {"z3": value})
        result, error = self.engine.numerical_evaluation(
            E("user_parameters::mass"), {"user_parameters::mass": E("3/2")}
        )
        self.assertIsNone(error)
        self.assertEqual(
            result.to_expression(),
            self.engine.numerical_result_from_expression(E("3/2")).to_expression(),
        )


class SavedQ0ArtifactSmoke(unittest.TestCase):
    """Inspect the accepted result offline; never regenerate the physical fixture."""

    @classmethod
    def setUpClass(cls):
        cls.saved = SavedQ0Result.load(api=core)

    def test_actual_bundle_is_lossless_without_regeneration(self):
        with patch.object(
            Q0Generator,
            "generate",
            side_effect=AssertionError("loading must not regenerate coefficients"),
        ):
            saved = SavedQ0Result.load(api=core)
            report = saved.verify()
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["numerical_evaluations"], 0)
        self.assertFalse(report["master_values_evaluated"])
        self.assertEqual(report["master_symbols"], 55)
        self.assertEqual(report["master_products"], 58)
        self.assertEqual(report["scalar_coefficients"], 449)
        self.assertEqual(report["epsilon_maximum"], 5)
        self.assertEqual(report["remainder_starts_at"], 6)
        self.assertEqual(report["component_storage_bits"], [998, 1104])
        self.assertEqual(report["minimum_tracked_nonzero_decimal_digits"], 330)
        self.assertEqual(report["minimum_exported_nonzero_decimal_digits"], 300)
        native = saved.expression("numerical", format="native")
        require_raw_equal(
            native,
            saved.expression("numerical", format="json"),
            "accepted numerical JSON/native precision parity",
        )
        require_raw_equal(
            native,
            reconstruct_coefficient_table(saved.coefficient_table, api=core),
            "public exact-ratio reconstruction",
        )
        self.assertEqual(
            hashlib.sha256(saved.display_text("retained").encode()).hexdigest(),
            "da909dc0b8f199a9ce46c27368bd2d1c88b277d72e765ab70f9a83680caa1212",
        )
        coupling = E("UFO::GC_11")
        self.assertEqual(
            [monomial for monomial, _ in native.coefficient_list(coupling)],
            [coupling**10],
        )

    def test_all_recorded_component_bounds_and_coefficientwise_orders(self):
        table = self.saved.coefficient_table
        self.assertEqual(table["cutoff"]["included_maximum_power"], 5)
        self.assertEqual(table["cutoff"]["remainder_starts_at"], 6)
        self.assertFalse(table["unknown_master_values_evaluated"])
        self.assertFalse(table["normalization_applied_again"])
        components = exact_zeros = 0
        for row in table["products"]:
            self.assertEqual(row["requested_decimal_precision"], 300)
            self.assertEqual(row["fixed_guard_decimal_digits"], 32)
            self.assertEqual(row["automatic_retries"], 0)
            self.assertEqual(
                row["remainder"], {"coefficientwise": True, "epsilon_power": 6}
            )
            self.assertEqual(row["exact_recurrence_check"]["status"], "PASS")
            self.assertGreaterEqual(
                row["exact_recurrence_check"]["remainder_order_at_least"], 6
            )
            for number in row["coefficients"]:
                self.assertLessEqual(number["epsilon_power"], 5)
                for part in ("real", "imaginary"):
                    component = number[part]
                    center = component["center_exact_ratio"]
                    center = Fraction(
                        int(center["numerator"]), int(center["denominator"])
                    )
                    error = component[
                        "absolute_error_exact_ratio_including_export_rounding"
                    ]
                    error = Fraction(int(error["numerator"]), int(error["denominator"]))
                    self.assertGreaterEqual(error, 0)
                    self.assertTrue(component["error_quality"])
                    components += 1
                    if component["exact_zero"]:
                        self.assertEqual(center, 0)
                        self.assertEqual(error, 0)
                        self.assertEqual(component["float_binary_precision"], 1104)
                        exact_zeros += 1
                    else:
                        self.assertEqual(component["float_binary_precision"], 998)
                        self.assertGreaterEqual(
                            component["tracked_significant_digits"], 330
                        )
                        self.assertEqual(
                            component["exported_reliable_decimal_digits_at_least"], 300
                        )
                        self.assertLessEqual(error * 10**300, abs(center) - error)
        self.assertEqual(components, 898)
        self.assertEqual(exact_zeros, 449)

    def test_metadata_properties_return_independent_copies(self):
        table = self.saved.coefficient_table
        provenance = self.saved.provenance
        before_table = copy.deepcopy(table)
        before_provenance = copy.deepcopy(provenance)
        table["products"][0]["coefficients"][0]["real"]["float_binary_precision"] = 2
        provenance["artifacts"].clear()
        self.assertEqual(self.saved.coefficient_table, before_table)
        self.assertEqual(self.saved.provenance, before_provenance)
        self.assertEqual(self.saved.summary()["scalar_coefficients"], 449)

    def test_formats_are_explicit_and_numeric_text_is_only_a_display(self):
        self.assertIn("RustRedMaster", self.saved.display_text())
        for variant, format, message in (
            ("unknown", "native", "unknown expression variant"),
            ("retained", "pickle", "unknown expression format"),
            ("retained", "json", "numerical scalar coefficients only"),
            ("numerical", "text", "human display"),
        ):
            with (
                self.subTest(variant=variant, format=format),
                self.assertRaisesRegex(CoefficientArtifactError, message),
            ):
                self.saved.expression(variant, format=format)

    def test_corrupted_native_payload_never_falls_back_to_other_formats(self):
        with tempfile.TemporaryDirectory(prefix="vakint-q0-corrupt-") as temporary:
            directory = Path(temporary) / "copy"
            shutil.copytree(self.saved.directory, directory)
            saved = SavedQ0Result.load(directory, api=core)
            path = directory / "q0-coefficients-300digits-through-epsilon5.atom"
            payload = bytearray(path.read_bytes())
            payload[-1] ^= 1
            path.write_bytes(payload)
            self.assertTrue(path.with_suffix(".expr").is_file())
            with self.assertRaisesRegex(CoefficientArtifactError, "checksum mismatch"):
                saved.expression("numerical", format="native")
            with (
                patch.object(
                    Q0Generator,
                    "generate",
                    side_effect=AssertionError(
                        "corruption must not trigger regeneration"
                    ),
                ),
                self.assertRaisesRegex(CoefficientArtifactError, "checksum mismatch"),
            ):
                SavedQ0Result.load(directory, api=core)

    def test_exact_zero_storage_precision_is_part_of_raw_atom_identity(self):
        real = Float.from_ratio(1, 8, precision=998)
        preserved = core.Expression.num(
            ComplexFloat(real, Float.from_ratio(0, 1, precision=1104))
        )
        shortened = core.Expression.num(
            ComplexFloat(real, Float.from_ratio(0, 1, precision=998))
        )
        self.assertEqual(preserved, shortened)  # Numeric equality is by exact value.
        with self.assertRaisesRegex(CoefficientArtifactError, "exact-zero storage"):
            require_raw_equal(preserved, shortened, "exact-zero storage changed")

    def test_factorwise_exact_controls_cover_poles_and_cancellation(self):
        report = q0_series_probe(core)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(len(report["cases"]), 7)
        self.assertEqual(report["numerical_evaluations"], 0)
        self.assertEqual(
            report["cases"]["sixth_order_pole"]["normalization_depth_required"], 11
        )
        cancelled = report["cases"]["cancelled_through_cutoff"]
        self.assertGreaterEqual(cancelled["valuation_lower_bound"], 6)
        self.assertFalse(cancelled["identically_zero_proven"])
        for case in report["cases"].values():
            self.assertGreaterEqual(case["remainder_order_at_least"], 6)

    def test_public_generator_handles_one_synthetic_master_product_without_ibps(self):
        first_product = self.saved.coefficient_table["products"][0]["master_product"]
        master_text = re.search(r"vakint::\{\}::RustRedMaster\([^()]*\)", first_product)
        self.assertIsNotNone(master_text)
        master = E(master_text.group())
        epsilon, gamma = E("vakint::ε"), E("symbolica::γ")
        mass, coupling = E("vakint::mUV"), E("UFO::GC_11")
        retained = (
            coupling**10
            * mass**2
            * (mass**2) ** (-5 * epsilon)
            * master**2
            * E("𝑒") ** (-5 * gamma * epsilon)
            * epsilon**-6
        )
        with tempfile.TemporaryDirectory(prefix="vakint-q0-synthetic-") as temporary:
            source = Path(temporary) / "synthetic-formal-input.expr"
            source.write_text(retained.to_canonical_string())
            output = Path(temporary) / "generated"
            result = Q0Generator(core, Vakint).generate(source, output)
            self.assertEqual(result, output)
            saved = SavedQ0Result.load(output, api=core)
            report = saved.verify()
            self.assertEqual(report["master_symbols"], 1)
            self.assertEqual(report["master_products"], 1)
            self.assertEqual(report["scalar_coefficients"], 12)
            self.assertEqual(report["epsilon_minimum"], -6)
            self.assertEqual(report["epsilon_maximum"], 5)
            self.assertEqual(report["remainder_starts_at"], 6)
            self.assertFalse(report["master_values_evaluated"])
            row = saved.coefficient_table["products"][0]
            require_raw_equal(
                E(row["master_product"]), master**2, "master multiplicity"
            )
            self.assertEqual(
                row["exact_recurrence_check"]["normalization_depth_required"], 11
            )
            with self.assertRaisesRegex(ValueError, "fresh"):
                Q0Generator(core, Vakint).generate(source, output)


@unittest.skipUnless(
    os.environ.get("VAKINT_TEST_EXTERNAL_DATA"),
    "set VAKINT_TEST_EXTERNAL_DATA to the certified four-loop checkout",
)
class LoadedMomentInspectionSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = Path(os.environ["VAKINT_TEST_EXTERNAL_DATA"])
        cls.engine = Vakint(evaluation_order=[], run_time_decimal_precision=16)
        cls.engine.load_external_data(cls.path, 4)
        manifest = json.loads((cls.path / "4L/manifest.json").read_text())
        cls.version = manifest["version"]
        metadata = next(row for row in manifest["masters"] if row["name"] == "P6")
        cls.master = cls.engine.get_master(MasterSymbol(metadata["master_id"]), prec=16)

    def test_discovery_and_all_four_finite_moment_records(self):
        symbols = set(self.engine.constant_symbols())
        for name in ("K5Moment(1)", "K5Moment(3)", "K5LogMoment(1)", "K5LogMoment(3)"):
            with self.subTest(name=name):
                symbol = ConstantSymbol(name)
                self.assertIn(symbol, symbols)
                record = self.engine.get_constant(symbol)
                self.assertEqual(record.requested_decimal_precision, 16)
                self.assertEqual(record.available_precision.kind, "FixedPrecision")
                self.assertEqual(
                    record.available_precision.max_decimal_precision, 20000
                )
                self.assertEqual(record.precision.real.quality.name, "Bounded")
                self.assertTrue(record.precision.imaginary.exact_zero)
                self.assertIsInstance(record.numerical_value, Float)
                self.assertEqual(record.provenance.version, self.version)

    def test_over_cap_lookup_preserves_source_limit_and_warning(self):
        record = self.engine.get_constant(ConstantSymbol("K5Moment(3)"), prec=20001)
        self.assertEqual(record.available_precision.max_decimal_precision, 20000)
        self.assertLessEqual(record.precision.real.significant_digits, 20000)
        self.assertTrue(
            any(
                warning.requested_decimal_precision == 20001
                and any(
                    source.available_decimal_precision == 20000
                    for source in warning.sources
                )
                for warning in record.warnings
            )
        )

    def test_master_formula_display_hides_only_the_snapshot_transport(self):
        coefficient = self.master.coefficients[0]
        symbolic = coefficient.symbolic_expression
        expression = symbolic.expression
        self.assertIn("RustRedConstant", str(expression))
        self.assertEqual(coefficient.available_precision.max_decimal_precision, 20000)
        self.assertIn("K5Moment", symbolic.to_html())
        self.assertNotIn("RustRedConstant", symbolic.to_html())
        self.assertNotIn("rustred_data_", symbolic.to_html())
        self.assertNotIn("RustRedConstant", repr(symbolic))
        self.assertEqual(expression, symbolic.expression)
        self.assertEqual(coefficient.provenance.version, self.version)
        self.assertIn("K5Moment", self.master.to_html())

    def test_retained_formula_requires_its_originating_snapshot(self):
        expression = self.master.coefficients[0].symbolic_expression.expression
        result, error = self.engine.numerical_evaluation(expression, {})
        self.assertIsNone(error)
        self.assertTrue(result.coefficients)
        other = Vakint(evaluation_order=[])
        with self.assertRaisesRegex(ValueError, "snapshot"):
            other.numerical_evaluation(expression, {})
        with self.assertRaisesRegex(ValueError, "registered constant"):
            self.engine.numerical_evaluation(
                E("K5Moment(1)", default_namespace="vakint"), {"K5Moment(1)": 1}
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
