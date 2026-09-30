"""Typed factorized contraction against the installed Python extension."""

import unittest

from symbolica import E
from symbolica.community import tensor as sp


class FactorizedContractionTests(unittest.TestCase):
    def test_contraction_outputs_and_dag_evaluation(self):
        rep = sp.Representation.mink(4)
        p, q, r = [
            sp.TensorName.vector(f"api_contract::{name}") for name in ("p", "q", "r")
        ]
        g = sp.TensorName.g().to_expression()
        first = g(p(rep).to_expression(), r(rep).to_expression())
        second = g(q(rep).to_expression(), r(rep).to_expression())
        expected = first + second
        i = rep("i")
        source = sp.TensorExpression(
            (p(i).to_expression() + q(i).to_expression()) * r(i).to_expression()
        )
        uncertified = sp.AliasedTensorExpression.from_expression(source)
        self.assertFalse(uncertified.contraction_complete)
        with self.assertRaises(AttributeError):
            uncertified.contraction_complete = True
        for order in (None, [0, 1], [1, 0]):
            with self.subTest(order=order):
                result = source.contract(order=order)
                self.assertIsInstance(result, sp.AliasedTensorExpression)
                self.assertTrue(result.contraction_complete)
                self.assertGreater(len(result.aliases), 0)
                self.assertEqual(result.root.rank, 0)
                self.assertEqual(result.expand().to_expression(), expected)
                rerun = result.contract()
                self.assertTrue(rerun.contraction_complete)
                self.assertEqual(rerun.expand(), result.expand())
                self.assertEqual(rerun.root.structure, result.root.structure)
                self.assertEqual(
                    [
                        (handle.to_expression(), body.to_expression())
                        for handle, body in rerun.aliases
                    ],
                    [
                        (handle.to_expression(), body.to_expression())
                        for handle, body in result.aliases
                    ],
                )
                self.assertEqual(
                    result.to_expression().to_expression().expand(), expected
                )
                self.assertEqual(
                    source.contract(order=order).expand().to_expression(), expected
                )
                evaluator = result.evaluator([first, second], iterations=1, n_cores=1)
                self.assertEqual(evaluator.evaluate([[2.0, 3.0]])[0][0], 5.0)
        for order in ([0], [0, 0], [0, 2]):
            with self.assertRaisesRegex(
                ValueError, "every normalized top-level factor"
            ):
                source.contract(order=order)
        with self.assertRaisesRegex(TypeError, "output"):
            source.contract(output="invalid")
        with self.assertRaises(TypeError):
            source.contract(left=0, right=0)

    def test_independent_reduction_settings_and_notation(self):
        spin = sp.Representation.bis(4)
        lorentz = sp.Representation.mink(4)
        gamma = sp.TensorExpression.dirac_gamma(4)
        source = gamma(spin("i"), spin("j"), lorentz("mu")) * gamma(
            spin("j"), spin("i"), lorentz("nu")
        )
        untouched = source.contract(representations=[]).to_expression()
        self.assertEqual(untouched, source)
        collected = source.contract().to_expression()
        self.assertIn("trace", collected.to_expression().format_plain())
        unfolded = collected.undo_trace()
        self.assertNotIn("trace", unfolded.to_expression().format_plain())
        self.assertIn("chain", unfolded.to_expression().format_plain())
        explicit = unfolded.undo_chain()
        self.assertNotIn("chain", explicit.to_expression().format_plain())
        # Independent finite HEP components compare represented tensors, not printed dummy labels.
        self.assertEqual(explicit.components(), collected.components())
        self.assertEqual(source.reduce_algebra().to_expression(), source)
        reduced = source.reduce_algebra(gamma=True)
        self.assertEqual(reduced.reduction_status, sp.ReductionStatus.Complete)
        self.assertEqual(reduced.expand().components(), collected.components())
        for obsolete in (
            "simplify",
            "simplify_gamma",
            "simplify_color",
            "simplify_epsilon",
            "undo_all",
            "chainify",
        ):
            self.assertFalse(hasattr(source, obsolete))

    def test_keyword_configs_and_family_validation(self):
        import inspect

        tensor = sp.TensorExpression(3)
        for carrier in (tensor, sp.AliasedTensorExpression.from_expression(tensor)):
            for operation in (carrier.contract, carrier.reduce_algebra):
                signature = inspect.signature(operation)
                self.assertTrue(
                    all(
                        p.kind == inspect.Parameter.KEYWORD_ONLY
                        for p in signature.parameters.values()
                    )
                )
                with self.assertRaises(TypeError):
                    operation({})
                with self.assertRaises(TypeError):
                    operation(settings={})
                with self.assertRaises(TypeError):
                    operation(max_passes=0)
            for family, options in (
                (
                    "gamma",
                    dict(
                        gamma_output="reduced",
                        gamma_ordering="repeated_pairs",
                        gamma_evaluate_traces=False,
                        gamma0=False,
                        gamma_conjugate=False,
                        gamma_expand_three_gamma_epsilon=False,
                    ),
                ),
                (
                    "color",
                    dict(
                        color_evaluate_traces=False,
                        color_expand_fierz=False,
                        color_substitute_cof_dimension_invariants=False,
                    ),
                ),
            ):
                for option, value in options.items():
                    with self.subTest(option=option):
                        with self.assertRaisesRegex(
                            ValueError, f"{option} requires {family}=True"
                        ):
                            carrier.reduce_algebra(**{option: value})
                config = {family: True, **options}
                self.assertEqual(
                    carrier.reduce_algebra(**config).to_expression(), tensor
                )
            for option in ("gamma_ordering", "gamma_output"):
                with self.assertRaisesRegex(ValueError, option):
                    carrier.reduce_algebra(gamma=True, **{option: "invalid"})
            for representations in (
                None,
                [],
                [sp.Representation.mink(4)],
                [sp.Representation.mink(4).name],
            ):
                config = dict(representations=representations, collect_traces=False)
                self.assertEqual(carrier.contract(**config).to_expression(), tensor)
        for obsolete in (
            "ContractSettings",
            "AlgebraSettings",
            "GammaSimplifySettings",
            "ColorSimplifySettings",
            "GammaChainOrdering",
        ):
            self.assertFalse(hasattr(sp, obsolete), obsolete)

    def test_scalar_product_notation_is_distinct_from_index_contraction(self):
        rep = sp.Representation.mink(4)
        p = sp.TensorName.vector("notation_surface::p")
        q = sp.TensorName.vector("notation_surface::q")
        source = sp.TensorExpression(
            p(rep("mu")).to_expression() * q(rep("mu")).to_expression()
        )
        self.assertEqual(source.to_dots(), source)
        contracted = source.contract().to_expression()
        self.assertEqual(contracted.to_expression().get_name(), "spenso::g")
        self.assertEqual(contracted.to_dots().to_expression().get_name(), "spenso::dot")
        with self.assertRaisesRegex(ValueError, "final structural port"):
            sp.TensorExpression(p.to_expression()(q(rep).to_expression()))

    def test_unfolded_power_scopes_against_integer_component_oracle(self):
        rep = sp.Representation.euc(2)
        p = sp.TensorName.vector("unfolded_power_oracle::p")
        q = sp.TensorName.vector("unfolded_power_oracle::q")
        values_p, values_q = (2, 3), (5, 7)
        library = sp.TensorLibrary.hep_lib_atom()
        library.register(sp.Tensor.dense(p(rep), list(values_p)))
        library.register(sp.Tensor.dense(q(rep), list(values_q)))
        scalar = sp.dot(p(rep), q(rep))
        expected_dot = sum(a * b for a, b in zip(values_p, values_q, strict=True))
        for exponent in (2, 3):
            opened = ((scalar + 1) ** exponent).undo_dots()
            self.assertEqual(
                opened.components(library), [E(str((expected_dot + 1) ** exponent))]
            )
        first, second = scalar.undo_dots(), scalar.undo_dots()
        # Multiply raw atoms so tensor composition cannot repair a collision.
        product = sp.TensorExpression(first.to_expression() * second.to_expression())
        self.assertEqual(product.components(library), [E(str(expected_dot**2))])

    def test_contraction_preserves_unresolved_order_metadata_and_typed_zero(self):
        rep = sp.Representation.euc(2)
        other = sp.Representation.mink(4)
        opened = sp.TensorName("api_contract::T")(7, rep, other)
        permuted = opened.permute_axes([1, 0])
        for value in (opened, permuted, 0 * opened, 0 * permuted):
            with self.subTest(zero=value.to_expression() == E("0")):
                result = value.contract()
                self.assertEqual(result.root.structure, value.structure)
                self.assertEqual(result.to_expression().structure, value.structure)
                self.assertEqual(
                    result.to_expression().to_expression(), value.to_expression()
                )


if __name__ == "__main__":
    unittest.main()
