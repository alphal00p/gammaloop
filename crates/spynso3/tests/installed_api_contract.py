"""Typed factorized contraction against the installed Python extension."""

import unittest
from pathlib import Path

from symbolica import E, S
from symbolica.community import tensor as sp


class FactorizedContractionTests(unittest.TestCase):
    def test_contraction_outputs_and_evaluation(self):
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
        self.assertFalse(source.contraction_complete)
        with self.assertRaises(AttributeError):
            source.contraction_complete = True
        for order in (None, [0, 1], [1, 0]):
            with self.subTest(order=order):
                result = source.contract(order=order)
                self.assertIsInstance(result, sp.TensorExpression)
                self.assertTrue(result.contraction_complete)
                self.assertFalse(hasattr(result, "aliases"))
                self.assertEqual(result.rank, 0)
                self.assertEqual(result.expand().to_expression(), expected)
                rerun = result.contract()
                self.assertTrue(rerun.contraction_complete)
                self.assertEqual(rerun.expand(), result.expand())
                self.assertEqual(rerun.structure, result.structure)
                self.assertEqual(result.to_expression().expand(), expected)
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

    def test_vector_dot_preserves_closed_spinor_connections(self):
        from itertools import permutations

        lorentz = sp.Representation.mink(4)
        gamma = sp.TensorExpression.dirac_gamma(4)
        p = sp.TensorName.vector("closed_spinor_dot::p")(lorentz)
        q = sp.TensorName.vector("closed_spinor_dot::q")(lorentz)
        factors = [gamma("a", "b", "mu"), p("mu"), gamma("b", "a", "nu"), q("nu")]
        expected = sp.TensorExpression(
            E(
                "spenso::trace(spenso::bis(4),spenso::cyclic("
                "spenso::gamma(spenso::in,spenso::out,closed_spinor_dot::p(spenso::mink(4))),"
                "spenso::gamma(spenso::in,spenso::out,closed_spinor_dot::q(spenso::mink(4)))))"
            )
        )
        # Build the component oracle from explicitly indexed scalar atoms, so
        # intake cannot share the rank-one tensor multiplication that made the dot.
        indexed = factors[0].to_expression()
        for factor in factors[1:]:
            indexed = indexed * factor.to_expression()
        oracle = sp.TensorExpression(indexed).components()
        for order in permutations(range(4)):
            with self.subTest(order=order):
                source = factors[order[0]]
                for position in order[1:]:
                    source = source * factors[position]
                result = source.contract()
                self.assertEqual(result, expected)
                self.assertEqual(result.components(), oracle)
                self.assertEqual(result.reduction_status, sp.ReductionStatus.Complete)
                self.assertEqual(result.contract(), expected)
                self.assertEqual(
                    source.simplify_algebra(gamma=False, color=False), expected
                )
                capped = source.contract(max_steps_per_domain=0)
                self.assertEqual(capped, source)
                self.assertEqual(capped.contract(), expected)
                uncollected = source.contract(collect_traces=False)
                self.assertNotIn("trace", uncollected.to_expression().format_plain())
                self.assertEqual(uncollected.contract(), expected)
                self.assertEqual(source.contract(representations=[]), source)

    def test_independent_reduction_settings_and_notation(self):
        spin = sp.Representation.bis(4)
        lorentz = sp.Representation.mink(4)
        gamma = sp.TensorExpression.dirac_gamma(4)
        source = gamma(spin("i"), spin("j"), lorentz("mu")) * gamma(
            spin("j"), spin("i"), lorentz("nu")
        )
        untouched = source.contract(representations=[])
        self.assertEqual(untouched, source)
        collected = source.contract()
        self.assertIn("trace", collected.to_expression().format_plain())
        unfolded = collected.undo_trace()
        self.assertNotIn("trace", unfolded.to_expression().format_plain())
        self.assertIn("chain", unfolded.to_expression().format_plain())
        explicit = unfolded.undo_chain()
        self.assertNotIn("chain", explicit.to_expression().format_plain())
        # Independent finite HEP components compare represented tensors, not printed dummy labels.
        self.assertEqual(explicit.components(), collected.components())
        self.assertEqual(
            source.simplify_algebra(gamma=False, color=False, contract="none"), source
        )
        reduced = source.simplify_algebra(gamma=True)
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

    def test_gamma_and_color_identities_are_enabled_by_default(self):
        gamma = sp.TensorExpression.dirac_gamma(4)
        word = gamma("a", "b", "mu") * gamma("b", "a", "mu")
        reduced = word.simplify_algebra()
        self.assertEqual(reduced.to_expression(), E("16"))
        self.assertEqual(reduced.components(), word.components())
        self.assertEqual(reduced.reduction_status, sp.ReductionStatus.Complete)
        self.assertIn(
            "trace", word.simplify_algebra(gamma=False).to_expression().format_plain()
        )

        generator_trace = sp.TensorExpression.color_t(8, 3)("a", "i", "i")
        reduced = generator_trace.simplify_algebra()
        self.assertEqual(reduced.to_expression(), E("0"))
        self.assertEqual(reduced.structure, generator_trace.structure)
        self.assertEqual(reduced.components(), generator_trace.components())
        self.assertEqual(reduced.reduction_status, sp.ReductionStatus.Complete)
        self.assertNotEqual(
            generator_trace.simplify_algebra(color=False).to_expression(), E("0")
        )

    def test_fk0032_reduces_color_and_lorentz_components_without_expansion(self):
        from symbolica.community.hepkit import Model

        # Register the generated momentum and compound-index heads. The saved
        # original numerator avoids generating all 4970 four-loop diagrams.
        Model.qcd()
        fixture = (
            Path(__file__).resolve().parents[2]
            / "idenso/benches/data/fk0032-projected.expr"
        )
        original = E(fixture.read_text())
        na = S("fk0032::N_A")
        # Change the typed slots before admission and normalize the colour
        # projector by N_A. No SU(N) relation between N_A and C_A is assumed.
        generic = original.replace(
            E("spenso::coad(8,index_)"),
            E("spenso::coad(fk0032::N_A,index_)"),
        ) * (8 / na)
        lorentz = sp.Representation.mink(S("D"))
        k = sp.TensorName("gammalooprs::K")
        p = sp.TensorName("gammalooprs::P")
        k0, k1, k2, k3 = [k(index, lorentz) for index in range(4)]
        p0 = p(0, lorentz)

        def dot(left, right):
            return sp.dot(left, right).to_expression()

        # Contract the four independent Lorentz pairs directly. Colour is
        # independently fixed by exact Gell-Mann components and the adjoint
        # four-generator trace identity: N_A*C_A^4/24 - d_A^{abcd}d_A^{abcd}.
        factors = [
            dot(p0, p0) + dot(k0, p0) - dot(k1, p0) - dot(k0, k1),
            dot(k2, k2) - dot(k1, k2) - dot(k2, k3) + dot(k2, p0),
            dot(k0, k0) - dot(k0, k1) - dot(k0, k3) + dot(k0, p0),
            dot(k1, k3),
        ]
        kinematics = E("1")
        for factor in factors:
            kinematics *= factor
        color = E(
            "fk0032::N_A*spenso::cas(2,spenso::coad(fk0032::N_A))^4/24"
            "-spenso::gram(4,spenso::coad(fk0032::N_A),spenso::coad(fk0032::N_A))"
        )
        cases = (
            (na, generic, False),
            (E("8"), original, False),
            (E("8"), original, True),
        )
        for mode in ("fully", "dots"):
            for dimension, expression, explicit_su3 in cases:
                source = sp.TensorExpression(expression)
                options = {"contract": mode}
                if explicit_su3:
                    options["color_substitute_cof_dimension_invariants"] = True
                with self.subTest(
                    mode=mode, dimension=dimension, explicit_su3=explicit_su3
                ):
                    result = source.simplify_algebra(**options)
                    self.assertEqual(
                        result.reduction_status, sp.ReductionStatus.Complete
                    )
                    self.assertTrue(result.contraction_complete)
                    self.assertEqual(result.structure.slots(), source.structure.slots())
                    self.assertEqual(result.rank, 0)
                    self.assertEqual(result.simplify_algebra(**options), result)
                    # This only converts existing scalar-product notation.
                    raw = result.to_dots().to_expression()
                    expected = (
                        E("𝑖*UFO::G^8")
                        / dimension
                        * (-108 if explicit_su3 else color.replace(na, dimension))
                        * kinematics
                    )
                    self.assertEqual(
                        (raw - expected).collect_factors().factor(), E("0")
                    )
                    for factor in factors:
                        self.assertEqual(
                            raw.replace(factor, E("0")).replace(-factor, E("0")), E("0")
                        )
                    for head in ("spenso::f(", "spenso::chain(", "spenso::trace("):
                        self.assertNotIn(head, raw.format_plain())
                    if explicit_su3:
                        self.assertNotIn("spenso::gram(", raw.format_plain())
                    else:
                        self.assertIn("spenso::gram(", raw.format_plain())

    def test_algebra_contraction_modes(self):
        from symbolica import S

        color = sp.Representation.cof(S("algebra_modes::Nc"))
        lorentz = sp.Representation.mink(4)
        metric = sp.TensorExpression.g(color, color.dual())
        source = metric("i", "j") * metric("j", "i")
        nc = S("algebra_modes::Nc")
        full = source.simplify_algebra()
        self.assertEqual(full.to_expression(), nc)
        self.assertEqual(full.reduction_status, sp.ReductionStatus.Complete)
        self.assertEqual(full.simplify_algebra(), full)
        self.assertEqual(source.simplify_algebra(max_steps_per_domain=None), full)
        self.assertEqual(source.contract(max_steps_per_domain=None).to_expression(), nc)
        for mode, options in (
            ("none", {}),
            ("selected", {"representations": []}),
            ("selected", {"representations": [lorentz]}),
        ):
            self.assertEqual(
                source.simplify_algebra(
                    gamma=False, color=False, contract=mode, **options
                ),
                source,
            )
        selected = source.simplify_algebra(
            contract="selected", representations=[color.name]
        )
        self.assertEqual(selected.to_expression(), nc)
        capped = source.simplify_algebra(max_steps_per_domain=0)
        self.assertEqual(capped.reduction_status, sp.ReductionStatus.Capped)
        self.assertEqual(capped, source)
        self.assertEqual(capped.simplify_algebra().to_expression(), nc)
        for options in (
            {"contract": "invalid"},
            {"contract": "selected"},
            {"contract": "fully", "representations": []},
            {"contract": "minimal", "representations": []},
            {"contract": "none", "representations": [color]},
        ):
            with self.assertRaises(ValueError):
                source.simplify_algebra(**options)

        # Algebra prerequisites ignore the optional structural representation filter.
        gamma = sp.TensorExpression.dirac_gamma(4)
        spin = sp.Representation.bis(4)
        word = gamma("a", "b", "mu") * gamma("b", "c", "mu")
        expected = 4 * sp.TensorExpression.g(spin)("a", "c")
        for options in (
            {"contract": "none"},
            {"contract": "selected", "representations": []},
        ):
            result = word.simplify_algebra(gamma=True, **options)
            self.assertEqual(result.components(), expected.components())
        # The closed colour delta is independent of gamma work, but full cleanup
        # still produces Nc without enabling any generator or Fierz identities.
        mixed = source * gamma("a", "b", "mu") * gamma("b", "a", "mu")
        self.assertEqual(
            mixed.simplify_algebra(gamma=True, color=False).to_expression(), 16 * nc
        )
        unevaluated = (gamma("a", "b", "mu") * gamma("b", "a", "nu")).simplify_algebra(
            gamma=False, color=False
        )
        self.assertIn("trace", unevaluated.to_expression().format_plain())

        p = sp.TensorName.vector("algebra_modes::p")(lorentz)
        q = sp.TensorName.vector("algebra_modes::q")(lorentz)
        compact = (p("mu") * q("mu")).contract()
        dots = compact.simplify_algebra(contract="dots")
        self.assertEqual(dots, compact.to_dots())
        self.assertIn("dot", dots.to_expression().format_plain())
        self.assertEqual(
            dots.simplify_algebra(
                contract="dots", max_steps_per_domain=0
            ).reduction_status,
            sp.ReductionStatus.Complete,
        )
        pending = compact.simplify_algebra(contract="dots", max_steps_per_domain=0)
        if compact != dots:
            self.assertEqual(pending.reduction_status, sp.ReductionStatus.Capped)
            self.assertEqual(pending, compact)
        self.assertEqual(
            pending.simplify_algebra(contract="dots"),
            dots,
        )

    def test_minimal_contraction_keeps_sum_products_and_identity_prerequisites(self):
        lorentz = sp.Representation.mink(4)
        p, q = (
            sp.TensorName.vector(f"minimal_modes::{name}")(lorentz)
            for name in ("p", "q")
        )
        metric = sp.TensorExpression.g(lorentz)
        matrix = sp.TensorName("minimal_modes::A")(lorentz, lorentz)
        source = (metric("mu", "nu") + matrix("mu", "nu")) * (p("mu") + q("mu"))
        original = source.to_expression()
        for minimal in (
            source.contract(expand=False),
            source.simplify_algebra(gamma=False, color=False, contract="minimal"),
        ):
            self.assertEqual(minimal.to_expression(), original)
            self.assertEqual(minimal.reduction_status, sp.ReductionStatus.Complete)
            self.assertFalse(minimal.contraction_complete)
            self.assertEqual(minimal.contract(), source.contract())
        self.assertEqual(source.contract(), source.contract(expand=True))
        self.assertNotEqual(source.contract().to_expression(), original)
        self.assertEqual(source.to_expression(), original)
        monomial = metric("mu", "nu") * p("mu")
        self.assertEqual(monomial.contract(expand=False), p("nu"))

        gamma = sp.TensorExpression.dirac_gamma(4)
        word = gamma("a", "b", "rho") * gamma("b", "a", "rho")
        reduced = word.simplify_algebra(contract="minimal")
        self.assertEqual(reduced.to_expression(), E("16"))
        self.assertEqual(reduced.components(), word.components())
        self.assertEqual(reduced.reduction_status, sp.ReductionStatus.Complete)
        self.assertEqual(
            (word * source)
            .simplify_algebra(color=False, contract="minimal")
            .to_expression(),
            (16 * source).to_expression(),
        )

        fundamental = sp.Representation.cof(3)
        delta = sp.TensorExpression.g(fundamental, fundamental.dual())
        cycle = delta("i", "j") * delta("j", "i")
        dimension = cycle.simplify_algebra(gamma=False, color=False, contract="minimal")
        self.assertEqual(dimension.to_expression(), E("3"))
        self.assertEqual(dimension.reduction_status, sp.ReductionStatus.Complete)

    def test_full_contraction_preserves_disabled_color_identities(self):
        generator = sp.TensorExpression.color_t(8, 3)
        source = generator("a", "i", "j") * generator("a", "k", "l")
        expected = source.components()
        for options in (
            {"color": False},
            {"color": True, "color_expand_fierz": False},
            {
                "color": True,
                "color_expand_fierz": False,
                "contract": "selected",
                "representations": [sp.Representation.cof(3)],
            },
        ):
            result = source.simplify_algebra(**options)
            self.assertIn("t(", result.to_expression().format_plain())
            self.assertEqual(result.components(), expected)
            self.assertEqual(result.reduction_status, sp.ReductionStatus.Complete)
            rerun = result.simplify_algebra(**options)
            self.assertEqual(rerun.reduction_status, sp.ReductionStatus.Complete)
            self.assertEqual(rerun, result)

    def test_keyword_configs_and_family_validation(self):
        import inspect

        tensor = sp.TensorExpression(3)
        self.assertFalse(hasattr(tensor, "reduce_algebra"))
        defaults = inspect.signature(tensor.simplify_algebra).parameters
        self.assertIs(defaults["gamma"].default, True)
        self.assertIs(defaults["color"].default, True)
        self.assertIs(defaults["epsilon"].default, False)
        self.assertIs(
            inspect.signature(tensor.contract).parameters["expand"].default, True
        )
        for operation in (tensor.contract, tensor.simplify_algebra):
            signature = inspect.signature(operation)
            self.assertIsNone(signature.parameters["max_steps_per_domain"].default)
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
                {
                    "gamma_output": "reduced",
                    "gamma_ordering": "repeated_pairs",
                    "gamma_evaluate_traces": False,
                    "gamma0": False,
                    "gamma_conjugate": False,
                    "gamma_expand_three_gamma_epsilon": False,
                },
            ),
            (
                "color",
                {
                    "color_evaluate_traces": False,
                    "color_expand_fierz": False,
                    "color_substitute_cof_dimension_invariants": False,
                },
            ),
        ):
            for option, value in options.items():
                with (
                    self.subTest(option=option),
                    self.assertRaisesRegex(
                        ValueError, f"{option} requires {family}=True"
                    ),
                ):
                    tensor.simplify_algebra(**{family: False, option: value})
            config = {family: True, **options}
            self.assertEqual(tensor.simplify_algebra(**config), tensor)
        for option in ("gamma_ordering", "gamma_output"):
            with self.assertRaisesRegex(ValueError, option):
                tensor.simplify_algebra(gamma=True, **{option: "invalid"})
        for representations in (
            None,
            [],
            [sp.Representation.mink(4)],
            [sp.Representation.mink(4).name],
        ):
            config = {"representations": representations, "collect_traces": False}
            self.assertEqual(tensor.contract(**config), tensor)
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
        contracted = source.contract()
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
                self.assertEqual(result.structure, value.structure)
                self.assertEqual(result.to_expression(), value.to_expression())


if __name__ == "__main__":
    unittest.main()
