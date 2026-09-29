"""Explicit materialization after shared gamma identities in the installed extension."""

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
        gamma = sp.TensorExpression.gamma(dimension)
        return sp.trace(
            self.spin, *(gamma(sp.AUTO, sp.AUTO, rep(label)) for label in labels)
        )

    def test_settings_preserve_defaults_and_disabled_evaluation(self):
        defaults = sp.GammaSimplifySettings()
        self.assertEqual(defaults.output, "reduced")
        self.assertFalse(defaults.gamma0)
        self.assertFalse(defaults.conjugate)
        self.assertFalse(sp.GammaSimplifySettings.canonical().gamma0)
        selected = sp.GammaSimplifySettings(
            gamma0=True, conjugate=True, evaluate_traces=False
        )
        self.assertTrue(selected.gamma0)
        self.assertTrue(selected.conjugate)
        self.assertFalse(selected.evaluate_traces)
        self.assertIn("gamma0=True", repr(selected))
        with self.assertRaises(TypeError):
            sp.GammaSimplifySettings(expand_traces=True)
        trace = self.trace_word(self.dimension, [f"a{i}" for i in range(6)])
        self.assertEqual(
            trace.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand(),
            trace.reduce_algebra(sp.AlgebraSettings(gamma=defaults, epsilon=True))
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand(),
        )
        inert = trace.reduce_algebra(
            sp.AlgebraSettings(
                gamma=sp.GammaSimplifySettings(evaluate_traces=False), epsilon=True
            )
        ).contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
        self.assertEqual(inert.to_expression(), trace)
        self.assertEqual(inert.root.structure.slots, trace.structure.slots)

    def test_chain_output_is_distinct_from_disabled_trace_evaluation(self):
        chains = sp.GammaSimplifySettings(output="chains")
        self.assertEqual(chains.output, "chains")
        self.assertIn("output='chains'", repr(chains))
        with self.assertRaises(ValueError):
            sp.GammaSimplifySettings(output="expanded")
        gamma = sp.TensorExpression.gamma(4)
        source = gamma("chain_a", "chain_b", "chain_mu") * gamma(
            "chain_b", "chain_c", "chain_mu"
        )
        joined = source.reduce_algebra(sp.AlgebraSettings(gamma=chains))
        reduced = source.reduce_algebra(
            sp.AlgebraSettings(
                gamma=sp.GammaSimplifySettings(evaluate_traces=False), epsilon=True
            )
        ).contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
        self.assertNotEqual(joined.to_expression(), reduced.to_expression())
        self.assertEqual(joined.root.structure.slots, source.structure.slots)
        self.assertEqual(
            joined.reduce_algebra(sp.AlgebraSettings(gamma=chains)).to_expression(),
            joined.to_expression(),
        )
        self.assertEqual(
            joined.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand(),
            source.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
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
                    trace.reduce_algebra(
                        sp.AlgebraSettings(
                            gamma=sp.GammaSimplifySettings(), epsilon=True
                        )
                    )
                    .contract(
                        sp.ContractSettings(collect_chains=False, collect_traces=False)
                    )
                    .expand()
                    .to_expression()
                )
                source = self.scalar * trace
                result = source.reduce_algebra(
                    sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
                ).contract(
                    sp.ContractSettings(collect_chains=False, collect_traces=False)
                )
                self.assertIsInstance(result, sp.AliasedTensorExpression)
                self.assertEqual(result.root.structure.slots, source.structure.slots)
                expanded = result.expand()
                self.assertEqual(
                    expanded.to_expression(), (self.scalar * body).expand()
                )
                self.assertEqual(
                    result.reduce_algebra(
                        sp.AlgebraSettings(
                            gamma=sp.GammaSimplifySettings(), epsilon=True
                        )
                    )
                    .contract(
                        sp.ContractSettings(collect_chains=False, collect_traces=False)
                    )
                    .expand(),
                    expanded,
                )
                self.assertEqual(expanded.structure.slots, source.structure.slots)
                self.assertEqual(
                    result.to_expression().to_expression().expand(),
                    expanded.to_expression(),
                )

    def test_independent_traces_materialize_through_the_same_alias_owner(self):
        left = self.trace_word(self.dimension, [f"left{i}" for i in range(6)])
        right = self.trace_word(self.dimension, [f"right{i}" for i in range(6)])
        expected = (
            self.scalar
            * left.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand()
            .to_expression()
            * right.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand()
            .to_expression()
        ).expand()
        source = self.scalar * left * right
        result = source.reduce_algebra(
            sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
        ).contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
        self.assertEqual(result.expand().to_expression(), expected)
        self.assertEqual(
            result.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand(),
            result.expand(),
        )
        self.assertEqual(result.root.structure.slots, source.structure.slots)

    def test_pipeline_shares_alias_result_and_requires_explicit_expansion(self):
        trace = self.trace_word(self.dimension, [f"pipeline{i}" for i in range(6)])
        source = self.scalar * trace
        pipeline = sp.AlgebraSettings(gamma=sp.GammaSimplifySettings())
        self.assertTrue(pipeline.gamma.evaluate_traces)
        with self.assertRaises(TypeError):
            sp.AlgebraSettings(expand=True)
        capped = source.reduce_algebra(
            sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), max_passes=0)
        )
        self.assertEqual(capped.reduction_status, sp.ReductionStatus.Capped)
        self.assertEqual(capped.to_expression(), source)
        result = source.reduce_algebra(pipeline)
        self.assertIsInstance(result, sp.AliasedTensorExpression)
        self.assertEqual(
            result.expand(),
            source.reduce_algebra(
                sp.AlgebraSettings(gamma=sp.GammaSimplifySettings(), epsilon=True)
            )
            .contract(sp.ContractSettings(collect_chains=False, collect_traces=False))
            .expand(),
        )
        self.assertEqual(result.reduce_algebra(pipeline).expand(), result.expand())
        zero = (0 * source).reduce_algebra(pipeline)
        self.assertEqual(zero.expand().structure.slots, source.structure.slots)
        self.assertFalse(zero.expand())


if __name__ == "__main__":
    unittest.main()
