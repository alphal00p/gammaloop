"""Public API checks against one coherent Symbolica/Vakint extension.

Set VAKINT_TEST_EXTERNAL_DATA to a validated external-data checkout. These
tests never download or write into that checkout.
"""

import json
import os
import unittest
from pathlib import Path

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    from vakint_hepkit_standalone import E, Vakint, VakintEvaluationMethod
else:
    from symbolica import E
    from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

DATA = Path(os.environ["VAKINT_TEST_EXTERNAL_DATA"])
DATA_VERSION = json.loads((DATA / "manifest.json").read_text())["version"]


def engine(mode="symbolic"):
    return Vakint(
        evaluation_order=[VakintEvaluationMethod.new_rustred_method(master_mode=mode)],
        form_exe_path="/external-data-smoke-must-not-use-form",
        integral_normalization_factor="FMFTandMATAD",
        run_time_decimal_precision=80,
    )


class ExternalDataSmoke(unittest.TestCase):
    def test_local_checkout_load_is_cumulative_and_read_only(self):
        instance = engine()
        report = instance.load_official_external_data(2, local_artifacts_path=DATA)
        self.assertEqual(report["version"], DATA_VERSION)
        self.assertEqual([p["loop_count"] for p in report["loops"]], [1, 2])
        self.assertTrue(all(p["source"] == "local" for p in report["loops"]))
        self.assertTrue(all(p["cache_available"] for p in report["loops"]))
        self.assertFalse((DATA / ".vakint-external-data").exists())

    def test_prepared_cache_loading_can_be_disabled(self):
        instance, reference = engine(), engine()
        report = instance.load_official_external_data(
            2, local_artifacts_path=DATA, load_prepared_caches=False
        )
        default_report = reference.load_official_external_data(
            2, local_artifacts_path=DATA
        )
        self.assertFalse(report["load_prepared_caches"])
        self.assertTrue(default_report["load_prepared_caches"])
        self.assertEqual(report["version"], default_report["version"])
        self.assertEqual(report["path"], default_report["path"])
        for disabled, default in zip(
            report["loops"], default_report["loops"], strict=True
        ):
            self.assertTrue(disabled["cache_available"])
            self.assertTrue(default["cache_available"])
            self.assertFalse(disabled["cache_loaded"])
            self.assertTrue(default["cache_loaded"])
            self.assertEqual(
                {k: v for k, v in disabled.items() if k != "cache_loaded"},
                {k: v for k, v in default.items() if k != "cache_loaded"},
            )
            loop_count = disabled["loop_count"]
            self.assertFalse(instance.cache_statistics(loop_count)["attached"])
            self.assertTrue(reference.cache_statistics(loop_count)["attached"])
        integral = E("topo(I1L(1,1))", default_namespace="vakint")
        self.assertEqual(
            instance.evaluate_integral(integral), reference.evaluate_integral(integral)
        )

    def test_version_error_preserves_active_data(self):
        instance = engine()
        instance.load_official_external_data(2, local_artifacts_path=DATA)
        integral = E("topo(I1L(1,1))", default_namespace="vakint")
        before = instance.evaluate_integral(integral)
        with self.assertRaises(Exception) as error:
            instance.load_official_external_data(
                2, version="v999.0.0", local_artifacts_path=DATA
            )
        self.assertIn("version", str(error.exception))
        self.assertEqual(before, instance.evaluate_integral(integral))

    def test_symbolic_then_numeric_agrees_with_direct_numeric(self):
        symbolic, numerical = engine(), engine("numerical")
        for instance in (symbolic, numerical):
            instance.load_external_data(DATA, 2)
        integral = E("topo(I2L(1,1,1,1))", default_namespace="vakint")
        exact = symbolic.evaluate_integral(integral)
        direct = numerical.evaluate_integral(integral)
        left, error = symbolic.numerical_evaluation(exact, {"mursq": 1.0})
        self.assertIsNone(error)
        right, error = numerical.numerical_evaluation(direct, {"mursq": 1.0})
        self.assertIsNone(error)
        matches, detail = left.compare_to(right, relative_threshold=1e-60)
        self.assertTrue(matches, detail)

    def test_unknown_master_error_lists_identity(self):
        instance = engine()
        unknown = E("RustRedMaster(unknown_family,2,3)", default_namespace="vakint")
        with self.assertRaises(Exception) as error:
            instance.numerical_evaluation(unknown, {})
        self.assertIn("unknown_family", str(error.exception))
        self.assertIn("unresolved masters", str(error.exception))

    def test_unknown_master_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            VakintEvaluationMethod.new_rustred_method(master_mode="invalid")


if __name__ == "__main__":
    unittest.main(verbosity=2)
