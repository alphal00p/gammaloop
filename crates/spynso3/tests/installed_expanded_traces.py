"""Trace-local expansion contracts against the installed community extension."""

import unittest

from symbolica import S
from symbolica.community import spenso as sp


class ExpandedTraceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dimension = S("expanded_trace_api::D")
        cls.spin = sp.Representation.bis(4)
        cls.scalar = sum(S("expanded_trace_api::x", "expanded_trace_api::y")) ** 8
        cls.expanded = sp.GammaSimplifySettings(expand_traces=True)

    def trace_word(self, dimension, labels):
        rep = sp.Representation.mink(dimension)
        gamma = sp.TensorExpression.gamma(dimension)
        return sp.trace(
            self.spin,
            *(gamma(sp.AUTO, sp.AUTO, rep(label)) for label in labels),
        )

    def test_settings_preserve_defaults_and_disabled_evaluation(self):
        self.assertFalse(sp.GammaSimplifySettings().expand_traces)
        self.assertFalse(sp.GammaSimplifySettings.canonical().expand_traces)
        self.assertTrue(self.expanded.expand_traces)
        self.assertIn("expand_traces=True", repr(self.expanded))
        trace = self.trace_word(self.dimension, [f"a{i}" for i in range(6)])
        self.assertEqual(
            trace.simplify_gamma(), trace.simplify_gamma(sp.GammaSimplifySettings())
        )
        self.assertEqual(
            trace.simplify_gamma(
                sp.GammaSimplifySettings(evaluate_traces=False, expand_traces=True)
            ),
            trace.simplify_gamma(sp.GammaSimplifySettings(evaluate_traces=False)),
        )

    def test_expansion_preserves_scalar_spectators_and_fixed_points(self):
        for dimension in [4, self.dimension]:
            with self.subTest(dimension=dimension):
                trace = self.trace_word(
                    dimension, ["mu", "a", "b", "mu", "c", "d", "e", "f"]
                )
                body = trace.simplify_gamma().to_expression().expand()
                source = self.scalar * trace
                output = source.simplify_gamma(self.expanded)
                self.assertEqual(output.to_expression(), self.scalar * body)
                self.assertEqual(output.simplify_gamma(self.expanded), output)
                self.assertEqual(
                    trace.simplify_gamma(self.expanded).to_expression(), body
                )

    def test_independent_traces_keep_independent_expansion_boundaries(self):
        left = self.trace_word(self.dimension, [f"left{i}" for i in range(6)])
        right = self.trace_word(self.dimension, [f"right{i}" for i in range(6)])
        expected = (
            self.scalar
            * left.simplify_gamma().to_expression().expand()
            * right.simplify_gamma().to_expression().expand()
        )
        source = self.scalar * left * right
        result = source.simplify_gamma(self.expanded)
        self.assertEqual(result.to_expression(), expected)
        self.assertEqual(result.simplify_gamma(self.expanded), result)

    def test_pipeline_forwards_trace_expansion_without_global_expansion(self):
        trace = self.trace_word(self.dimension, [f"pipeline{i}" for i in range(6)])
        source = self.scalar * trace
        pipeline = sp.SimplifySettings(metrics=False, gamma=self.expanded)
        self.assertFalse(pipeline.expand)
        self.assertTrue(pipeline.gamma.expand_traces)
        result = source.simplify(pipeline)
        self.assertEqual(result, source.simplify_gamma(self.expanded))
        self.assertEqual(result.simplify(pipeline), result)


if __name__ == "__main__":
    unittest.main()
