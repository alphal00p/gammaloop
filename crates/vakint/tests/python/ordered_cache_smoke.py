"""Ordered step cache through the Python module of one coherently built host.

Uses only the vendored lower-loop programs: no external data, no FORM.
"""

import tempfile
import unittest

from symbolica import E
from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

INTEGRALS = [
    "topo(I2L(muvsq,3,2,1))",
    "topo(I3L(muvsq,2,2,2,1,1,1))",
    "dot(k(1),k(2))*dot(k(2),k(3))*topo(I3L(muvsq,2,1,1,1,1,1))",
]


def engine(**overrides):
    options = {
        "evaluation_order": [
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        "form_exe_path": "/ordered-cache-python-test-must-not-use-form",
        "integral_normalization_factor": "MSbar",
        "use_dot_product_notation": True,
    }
    options.update(overrides)
    return Vakint(**options)


def integrals():
    return [E(text, default_namespace="vakint") for text in INTEGRALS]


def evaluate_all(vakint):
    return [str(vakint.evaluate(integral)) for integral in integrals()]


class OrderedCacheSmoke(unittest.TestCase):
    def test_saved_cache_replays_in_a_new_engine_with_identical_results(self):
        uncached = evaluate_all(engine(ordered_cache="off"))
        self.assertEqual(engine(ordered_cache="off").ordered_cache_statistics()["programs"], [])
        with tempfile.TemporaryDirectory() as directory:
            writer = engine(ordered_cache_directory=directory)
            self.assertEqual(evaluate_all(writer), uncached)
            cold = writer.ordered_cache_statistics()
            self.assertGreater(cold["inserts"], 0)
            self.assertEqual(cold["inserts"], cold["misses"])
            # Autosave after every evaluation leaves nothing to write.
            self.assertEqual(writer.save_ordered_cache()["records_written"], 0)
            self.assertEqual(
                sorted(program["program"] for program in cold["programs"]),
                ["2L-unit-mass-vacuum-k3", "3L-unit-mass-vacuum-k6"],
            )
            reader = engine(
                ordered_cache="read-only", ordered_cache_directory=directory
            )
            loaded = reader.load_ordered_cache()
            self.assertEqual(len(loaded["programs"]), 3)
            self.assertEqual(evaluate_all(reader), uncached)
            warm = reader.ordered_cache_statistics()
            self.assertEqual(warm["misses"], 0)
            self.assertEqual(warm["inserts"], 0)
            self.assertGreater(warm["hits"], 0)
            reader.clear_ordered_cache()
            self.assertEqual(reader.ordered_cache_statistics()["programs"], [])

    def test_warming_saves_a_batch_and_options_are_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            vakint = engine(
                ordered_cache_directory=directory,
                ordered_cache_autosave="never",
                ordered_cache_ram_bytes=1 << 30,
            )
            report = vakint.warm_ordered_cache(integrals())
            self.assertEqual(report["integrals"], 3)
            self.assertGreater(report["inserts"], 0)
            self.assertEqual(report["save"]["records_written"], report["inserts"])
            fresh = engine(ordered_cache="read-only", ordered_cache_directory=directory)
            evaluate_all(fresh)
            self.assertEqual(fresh.ordered_cache_statistics()["misses"], 0)
        with self.assertRaises(ValueError):
            engine(ordered_cache="sometimes")
        with self.assertRaises(ValueError):
            engine(ordered_cache_autosave="always")


if __name__ == "__main__":
    unittest.main()
