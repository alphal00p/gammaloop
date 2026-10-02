"""Explicit expansion after gamma identities in the installed extension."""

import unittest

from symbolica import S
from symbolica.community import tensor as sp


class ExpandedTraceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dimension = S("expanded_trace_api::D")
        cls.spin = sp.Representation.bis(4)
        cls.scalar = sum(S("expanded_trace_api::x", "expanded_trace_api::y")) ** 8

    def trace_word(self, dimension, labels):
        rep = sp.Representation.mink(dimension)
        gamma = sp.TensorExpression.dirac_gamma(dimension)
        return sp.trace(
            self.spin, *(gamma(sp.AUTO, sp.AUTO, rep(label)) for label in labels)
        )

    def test_settings_preserve_defaults_and_disabled_evaluation(self):
        defaults = {
            "gamma": True,
            "gamma_output": "reduced",
            "gamma_ordering": "repeated_pairs",
            "gamma_evaluate_traces": True,
            "gamma0": False,
            "gamma_conjugate": False,
            "gamma_expand_three_gamma_epsilon": False,
        }
        with self.assertRaises(TypeError):
            self.trace_word(4, ["a", "b"]).simplify_algebra(
                gamma=True, gamma_expand_traces=True
            )
        trace = self.trace_word(self.dimension, [f"a{i}" for i in range(6)])
        self.assertEqual(
            trace.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand(),
            trace.simplify_algebra(**defaults, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand(),
        )
        inert = trace.simplify_algebra(
            gamma=True, gamma_evaluate_traces=False, epsilon=True
        ).contract(collect_chains=False, collect_traces=False)
        self.assertEqual(inert, trace)
        self.assertEqual(inert.structure.axes, trace.structure.axes)

    def test_chain_output_is_distinct_from_disabled_trace_evaluation(self):
        chains = {"gamma": True, "gamma_output": "chains"}
        with self.assertRaises(ValueError):
            self.trace_word(4, ["a", "b"]).simplify_algebra(
                gamma=True, gamma_output="expanded"
            )
        gamma = sp.TensorExpression.dirac_gamma(4)
        source = gamma("chain_a", "chain_b", "chain_mu") * gamma(
            "chain_b", "chain_c", "chain_mu"
        )
        joined = source.simplify_algebra(**chains)
        reduced = source.simplify_algebra(
            gamma=True, gamma_evaluate_traces=False, epsilon=True
        ).contract(collect_chains=False, collect_traces=False)
        self.assertNotEqual(joined.to_expression(), reduced.to_expression())
        self.assertEqual(joined.structure.axes, source.structure.axes)
        self.assertEqual(
            joined.simplify_algebra(**chains).to_expression(), joined.to_expression()
        )
        self.assertEqual(
            joined.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand(),
            source.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand(),
        )
        for value in [source, joined]:
            self.assertFalse(hasattr(value, "collect_gamma_chains"))
            self.assertFalse(hasattr(value, "collect_chains"))

    def test_explicit_expansion_distributes_spectators_and_preserves_interface(self):
        for dimension in [4, self.dimension]:
            with self.subTest(dimension=dimension):
                trace = self.trace_word(
                    dimension, ["mu", "a", "b", "mu", "c", "d", "e", "f"]
                )
                body = (
                    trace.simplify_algebra(gamma=True, epsilon=True)
                    .contract(collect_chains=False, collect_traces=False)
                    .expand()
                    .to_expression()
                )
                source = self.scalar * trace
                result = source.simplify_algebra(gamma=True, epsilon=True).contract(
                    collect_chains=False, collect_traces=False
                )
                self.assertIsInstance(result, sp.TensorExpression)
                self.assertEqual(result.structure.axes, source.structure.axes)
                expanded = result.expand()
                self.assertEqual(
                    expanded.to_expression(), (self.scalar * body).expand()
                )
                self.assertEqual(
                    result.simplify_algebra(gamma=True, epsilon=True)
                    .contract(collect_chains=False, collect_traces=False)
                    .expand(),
                    expanded,
                )
                self.assertEqual(expanded.structure.axes, source.structure.axes)
                self.assertEqual(
                    result.to_expression().expand(),
                    expanded.to_expression(),
                )

    def test_independent_traces_reduce_without_expanding_spectators(self):
        left = self.trace_word(self.dimension, [f"left{i}" for i in range(6)])
        right = self.trace_word(self.dimension, [f"right{i}" for i in range(6)])
        expected = (
            self.scalar
            * left.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand()
            .to_expression()
            * right.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand()
            .to_expression()
        ).expand()
        source = self.scalar * left * right
        result = source.simplify_algebra(gamma=True, epsilon=True).contract(
            collect_chains=False, collect_traces=False
        )
        self.assertEqual(result.expand().to_expression(), expected)
        self.assertEqual(
            result.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand(),
            result.expand(),
        )
        self.assertEqual(result.structure.axes, source.structure.axes)

    def test_pipeline_preserves_factorization_and_requires_explicit_expansion(self):
        trace = self.trace_word(self.dimension, [f"pipeline{i}" for i in range(6)])
        source = self.scalar * trace
        pipeline = {"gamma": True}
        with self.assertRaises(TypeError):
            source.simplify_algebra(expand=True)
        capped = source.simplify_algebra(gamma=True, max_steps_per_domain=0)
        self.assertEqual(capped.reduction_status, sp.ReductionStatus.Capped)
        self.assertEqual(capped, source)
        result = source.simplify_algebra(**pipeline)
        self.assertIsInstance(result, sp.TensorExpression)
        self.assertEqual(
            result.expand(),
            source.simplify_algebra(gamma=True, epsilon=True)
            .contract(collect_chains=False, collect_traces=False)
            .expand(),
        )
        self.assertEqual(result.simplify_algebra(**pipeline).expand(), result.expand())
        zero = (0 * source).simplify_algebra(**pipeline)
        self.assertEqual(zero.expand().structure.axes, source.structure.axes)
        self.assertFalse(zero.expand())


if __name__ == "__main__":
    unittest.main()
