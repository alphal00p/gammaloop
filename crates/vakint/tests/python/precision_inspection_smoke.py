"""Exercise the installed, coherent Symbolica/Vakint extension without FORM."""

import importlib.util
import json
import os
import unittest
from pathlib import Path
from unittest.mock import patch

from symbolica import ComplexFloat, E, Float
from symbolica.community.feynkit import FeynmanDiagram
from symbolica.community.hepkit.vakint import (
    ConstantSymbol,
    MasterSymbol,
    Vakint,
)


class PrecisionInspectionSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = Vakint(evaluation_order=[], run_time_decimal_precision=32)

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
        if importlib.util.find_spec("typst"):
            self.assertIn("<svg", html)
        else:
            self.assertIn("Graph preview unavailable", html)

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
