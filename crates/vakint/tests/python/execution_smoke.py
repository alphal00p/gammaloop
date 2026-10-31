"""Worker cores and Ordered weight arithmetic of the Python engine.

Every arithmetic and every core count must return the same exact result as
the exact arithmetic on one core. The default arithmetic is
modular_polynomial.
"""

import os
import unittest

from symbolica import E
from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod


def engine(**overrides):
    options = {
        "evaluation_order": [
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        "form_exe_path": "/execution-python-test-must-not-use-form",
        "integral_normalization_factor": "FMFTandMATAD",
        "use_dot_product_notation": True,
    }
    options.update(overrides)
    return Vakint(**options)


def lower_loop_input():
    return E(
        "topo(I1L(muvsq,3))+dot(k(1),k(1))*topo(I2L(muvsq,2,1,1))"
        "+topo(I3L(muvsq,2,1,1,1,1,1))+dot(k(1),k(2))*topo(I3L(muvsq,1,1,1,2,1,1))",
        default_namespace="vakint",
    )


def available_cores():
    return len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 1


class ExecutionSettingsTests(unittest.TestCase):
    def test_defaults(self):
        obj = engine()
        self.assertGreaterEqual(obj.n_cores, 1)
        self.assertLessEqual(obj.n_cores, min(10, available_cores()))
        self.assertEqual(obj.ordered_weight_arithmetic, "modular_polynomial")

    def test_canonical_arithmetic_text(self):
        for text, canonical in (
            ("modular", "modular_sampled"),
            ("modular_sampled:threads=1", "modular_sampled:threads=1"),
            ("modular:batch=64", "modular_sampled:batch_points=64"),
            ("modular_polynomial:lanes=2", "modular_polynomial:lanes=2"),
        ):
            with self.subTest(text=text):
                obj = engine(n_cores=1, ordered_weight_arithmetic=text)
                self.assertEqual(obj.ordered_weight_arithmetic, canonical)
                self.assertEqual(obj.n_cores, 1)

    def test_invalid_values_raise_value_error(self):
        for options in (
            {"n_cores": 0},
            {"n_cores": available_cores() + 1},
            {"ordered_weight_arithmetic": "dfs"},
            {"ordered_weight_arithmetic": "exact:threads=2"},
            {"ordered_weight_arithmetic": "modular:threads=0"},
        ):
            with self.subTest(options=options), self.assertRaises(ValueError):
                engine(**options)


class ExecutionResultTests(unittest.TestCase):
    def test_every_arithmetic_and_core_count_agrees(self):
        expression = lower_loop_input()
        reference = engine(
            n_cores=1, ordered_weight_arithmetic="exact"
        ).evaluate_integral(expression)
        cores = sorted({1, min(4, available_cores())})
        for n_cores in cores:
            for arithmetic in (None, "exact", "modular_sampled", "modular_polynomial"):
                with self.subTest(n_cores=n_cores, arithmetic=arithmetic):
                    options = {"n_cores": n_cores}
                    if arithmetic is not None:
                        options["ordered_weight_arithmetic"] = arithmetic
                    result = engine(**options).evaluate_integral(expression)
                    self.assertEqual(result, reference)


if __name__ == "__main__":
    unittest.main()
