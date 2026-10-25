"""Run the same native-type assertions in standalone and Community processes.

Set VAKINT_TEST_HOST=standalone for the private ``vakint_hepkit_standalone``
all-in-one host. Otherwise use the installed, single-kernel Community wheel.
Neither mode imports a second native extension. Optional VAKINT_TEST_ARTIFACTS
selects a validated local external-data checkout for the 5L catalogue checks.
"""

import builtins
import hashlib
import json
import os
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    from vakint_hepkit_standalone import (
        ConstantSymbol,
        DiagramError,
        E,
        Expression,
        FeynmanDiagram,
        Float,
        IntegralFamily,
        MasterSymbol,
        TensorExpression,
        Vakint,
        VakintEvaluationMethod,
        integral_from_diagram,
    )
else:
    from symbolica import E, Expression, Float
    from symbolica.community.hepkit import DiagramError, FeynmanDiagram, IntegralFamily
    from symbolica.community.hepkit.vakint import (
        ConstantSymbol,
        MasterSymbol,
        Vakint,
        VakintEvaluationMethod,
        integral_from_diagram,
    )
    from symbolica.community.tensor import TensorExpression


class CommunityInteroperability(unittest.TestCase):
    def test_shared_native_records_and_ingress_without_module_lookup(self):
        if os.environ.get("VAKINT_TEST_HOST") == "standalone":
            self.assertFalse(
                any(
                    name == "symbolica" or name.startswith("symbolica.")
                    for name in sys.modules
                )
            )
        original_import = builtins.__import__

        def block_components(name, *args, **kwargs):
            if name == "symbolica" or name.startswith("symbolica."):
                raise AssertionError(f"native binding imported {name}")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", block_components):
            engine = Vakint(
                evaluation_order=[
                    VakintEvaluationMethod.new_rustred_method(master_mode="retain")
                ],
                form_exe_path="/native-interoperability-must-not-use-form",
            )
            constant = engine.get_constant(ConstantSymbol("z3"), prec=32)
            self.assertIsInstance(constant.numerical_value, Float)
            self.assertEqual(constant.available_precision.kind, "ArbitraryPrecision")
            self.assertIn("Vakint constant", constant.to_html())
            master = next(
                candidate
                for symbol in engine.master_symbols()
                if (candidate := engine.get_master(symbol, prec=24)).loop_count == 1
            )
            diagram = master.graph
            self.assertIsInstance(diagram, FeynmanDiagram)
            family = diagram.propagator_family()
            self.assertIsInstance(family, IntegralFamily)
            numerator = diagram.numerator_expression()
            self.assertIsInstance(numerator, TensorExpression)
            self.assertIsInstance(numerator.to_expression(), Expression)
            rules = diagram.loop_momentum_basis.momentum_replacements()
            self.assertEqual(
                numerator.to_expression().replace_multiple(rules),
                diagram.loop_momentum_basis.route_expression(numerator.to_expression()),
            )
            for edge in diagram.internal_edges:
                self.assertIsInstance(edge.propagator_power, int)
            integral = integral_from_diagram(diagram, family, E("1"))
            self.assertIsInstance(integral.to_expression(), Expression)
            self.assertIsInstance(
                engine.evaluate_integral(integral.to_expression()), Expression
            )
            self.assertIn("feynkit", diagram.to_html())
            self.assertIn("Vakint master", master.to_html())
            with self.assertRaises(ValueError):
                integral_from_diagram(diagram, family, E("1"), powers=[True])
            with self.assertRaisesRegex(DiagramError, "independent diagram"):
                integral_from_diagram(diagram.subgraph(), family, E("1"))
            with self.assertRaisesRegex(ValueError, "unknown RustRed constant"):
                engine.get_constant(ConstantSymbol("interop_unknown_constant"))
            with self.assertRaisesRegex(ValueError, "unknown RustRed master"):
                engine.get_master(MasterSymbol("interop_unknown_master"))

    @unittest.skipUnless(
        os.environ.get("VAKINT_TEST_ARTIFACTS"), "local 5L data not selected"
    )
    def test_uncomputed_five_loop_catalogue_keeps_shared_graph_types(self):
        artifacts = Path(os.environ["VAKINT_TEST_ARTIFACTS"])
        manifest = json.loads((artifacts / "5L/manifest.json").read_text())
        expected_cache = manifest["cache"]
        engine = Vakint(evaluation_order=[])
        report = engine.load_official_external_data(
            5, local_artifacts_path=os.environ["VAKINT_TEST_ARTIFACTS"]
        )
        self.assertEqual(
            [item["loop_count"] for item in report["loops"]], [1, 2, 3, 4, 5]
        )
        cache = engine.cache_statistics(5)
        self.assertEqual(cache["attached"], expected_cache["available"])
        self.assertFalse(cache["ram_is_process_rss"])
        self.assertIsNotNone(cache["saved_application"])
        if expected_cache["available"]:
            self.assertTrue(report["loops"][-1]["cache_loaded"])
            self.assertEqual(
                expected_cache["format"], "rustred.saved-application-cache-file.v1"
            )
            for relative in (expected_cache["path"], expected_cache["metadata"]):
                self.assertEqual(
                    hashlib.sha256(
                        (artifacts / "5L" / relative).read_bytes()
                    ).hexdigest(),
                    manifest["files"][relative]["sha256"],
                )
            metadata = json.loads(
                (artifacts / "5L" / expected_cache["metadata"]).read_text()
            )
            single_file = metadata["preparation_report"]["options"]["single_file"]
            self.assertEqual(single_file["binding"], expected_cache["binding"])
            self.assertEqual(
                cache["saved_application"]["binding"], expected_cache["binding"]
            )
            file = cache["saved_application"]["file"]
            self.assertEqual(file["records"], single_file["records"])
            for field in ("compressed_bytes", "decoded_bytes"):
                self.assertEqual(file[field], expected_cache["identity"][field])
        else:
            self.assertFalse(report["loops"][-1]["cache_loaded"])
            self.assertIsNone(cache["saved_application"]["file"])
        self.assertEqual(cache["saved_application"]["native"]["cached_integrals"], 0)
        binding = cache["saved_application"]["binding"]
        engine.detach_cache(5)
        detached = engine.cache_statistics(5)
        self.assertFalse(detached["attached"])
        self.assertEqual(detached["saved_application"]["binding"], binding)
        count = 0
        for symbol in engine.master_symbols():
            master = engine.get_master(symbol)
            if master.loop_count != 5:
                continue
            count += 1
            self.assertFalse(master.is_computed)
            self.assertEqual(master.coefficients, {})
            self.assertIsInstance(master.graph, FeynmanDiagram)
            self.assertIn("digraph", master.graph.to_dot())
        self.assertEqual(count, 640)


if __name__ == "__main__":
    unittest.main(verbosity=2)
