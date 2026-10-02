"""Tensor API contracts against the installed extension and native compiler."""

import copy
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from symbolica import E, Expression, Replacement, S, Symbol, T
from symbolica.community import tensor as sp


class TensorOperationsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rep = sp.Representation.euc(2)
        cls.other = sp.Representation.mink(3)
        cls.A = sp.TensorName("api_operations::A")
        cls.x, cls.y = S("api_operations::x", "api_operations::y")

    def matrix(self, symbolic=False):
        values = (
            [self.x, E("1"), E("0"), self.x + 1] if symbolic else [1.0, 2.0, 3.0, 4.0]
        )
        return sp.Tensor.dense(self.A(self.rep, self.rep), values)

    def test_expression_products_preserve_tensor_interfaces(self):
        vector = self.A(self.rep)("i")
        zero = sp.TensorExpression(0, structure=vector.structure)
        for tensor in (vector, zero):
            for result in (self.x * tensor, tensor * self.x, self.x.__mul__(tensor)):
                with self.subTest(tensor=tensor, result=result):
                    self.assertIsInstance(result, sp.TensorExpression)
                    self.assertEqual(result.axes, tensor.axes)
                    self.assertEqual(
                        result.to_expression(), self.x * tensor.to_expression()
                    )
        raw = vector.to_expression()
        self.assertEqual(raw.__mul__(vector), vector * vector)

    def test_tensor_product_supports_the_generic_extension_protocol(self):
        class Product:
            def __symbolica_rmul__(self, left: Expression) -> tuple[Expression, str]:
                return left, "product"

        tensor = self.A(self.rep)
        result = tensor * Product()
        self.assertIs(result[0], tensor)
        self.assertEqual(result[1], "product")

    def test_gamma_product_prefers_the_unresolved_spinor_partner(self):
        spinor, vector = sp.Representation.bis(4), sp.Representation.mink(4)
        Q = sp.TensorName.vector("preferred_ports::Q")(spinor)
        G = sp.TensorName.vector("preferred_ports::G")(vector)
        gamma = sp.TensorExpression.dirac_gamma(4)
        left = Symbol.I * gamma(1, sp.AUTO, 1) * Q(1) * G(1)
        result = left * gamma(2, sp.AUTO, 2)
        shared = "preferred_ports::spinor"
        expected = Symbol.I * gamma(1, shared, 1) * Q(1) * G(1) * gamma(2, shared, 2)
        self.assertEqual(result.axes, (spinor(2), vector(2)))
        library = sp.TensorLibrary.hep_lib_atom()
        library.register(sp.Tensor.dense(Q, [E(str(i)) for i in (1, 2, 3, 4)]))
        library.register(sp.Tensor.dense(G, [E(str(i)) for i in (5, 6, 7, 8)]))
        self.assertEqual(result.to_tensor(library)[:], expected.to_tensor(library)[:])

    def test_slot_constructor_reuses_supported_index_labels(self):
        source = self.A(self.rep)
        scope = S("api_operations::slot_scope")
        scoped = source("i").wrap_indices(scope).axes[0].index
        for rep, dimension, dual in (
            (self.rep, 2, False),
            (sp.Representation.cof(3), 3, True),
        ):
            representation = rep.dual() if dual else rep
            descriptor = self.A(representation)
            for label in (E("7"), E("-7"), E("gammalooprs::hedge(2,1)"), scoped):
                with self.subTest(label=label, dual=dual):
                    admitted = descriptor(label)
                    slot = sp.Slot(rep.name.name, dimension, label, dual=dual)
                    self.assertEqual(slot, admitted.axes[0])
                    self.assertEqual(self.A(slot), admitted)
                    self.assertEqual(slot.index, label)

        tensor = sp.Tensor.dense(source, [E("2"), E("3")])
        slot = sp.Slot(self.rep.name.name, 2, scoped)
        norm = (tensor(slot) * tensor(slot)).to_tensor()
        self.assertEqual(norm.scalar(), E("13"))

    def test_metric_helpers_reuse_supported_index_labels(self):
        source = self.A(self.rep)
        scoped = (
            source("i").wrap_indices(S("api_operations::metric_scope")).axes[0].index
        )
        reference = S("api_operations::metric_reference")
        for method in (self.rep.g, self.rep.flat, self.rep.id):
            original = method(reference, "j")
            for label in (E("7"), E("gammalooprs::hedge(2,1)"), scoped):
                with self.subTest(method=method.__name__, label=label):
                    expected = original.rename_indices({reference: label})
                    result = method(label, "j")
                    self.assertEqual(result, expected)
                    self.assertEqual(
                        result.to_expression(),
                        original.to_expression().replace(reference, label),
                    )
                    self.assertEqual(
                        [axis.index for axis in result.axes],
                        [label, self.rep("j").index],
                    )

    def test_numeric_index_admission_preserves_identity_and_contractions(self):
        source = self.A(self.rep)
        reference = source("reference")
        original = reference.axes[0].index
        maximum = 2 * sys.maxsize + 1
        labels = (
            7,
            2**31,
            -(2**31),
            2**32 + 7,
            sys.maxsize,
            sys.maxsize + 1,
            -sys.maxsize - 1,
            maximum,
            -maximum,
        )
        library = sp.TensorLibrary.hep_lib_atom()
        library.register(sp.Tensor.dense(source, [E("2"), E("3")]))
        for label in labels:
            if abs(label) > maximum:
                continue
            encoded = E(str(label))
            raw = reference.to_expression().replace(original, encoded)
            for value in (source(label), source(encoded), sp.TensorExpression(raw)):
                with self.subTest(label=label, expression=value.to_expression()):
                    self.assertEqual(value.axes[0].index, encoded)
                    self.assertEqual(value.to_expression(), raw)
                    self.assertEqual(sp.TensorExpression(raw).axes, value.axes)
                    self.assertEqual(self.rep(label).index, encoded)
                    product = value * source(7)
                    self.assertEqual(product.rank, 0 if label == 7 else 2)
                    components = product.to_tensor(library)
                    self.assertEqual(
                        components[:],
                        [E("13")] if label == 7 else [E("4"), E("6"), E("6"), E("9")],
                    )

    def test_numeric_index_admission_rejects_lossy_rationals_and_overflow(self):
        source = self.A(self.rep)
        reference = source("reference")
        original = reference.axes[0].index
        beyond_range = 2 * sys.maxsize + 2
        for label in (beyond_range, -beyond_range):
            with self.assertRaises(ValueError):
                source(label)
        for spelling in (
            "65537/2",
            "1/65536",
            "-1/2",
            "1+𝑖",
            str(beyond_range),
            str(-beyond_range),
        ):
            encoded = E(spelling)
            with self.subTest(label=spelling):
                with self.assertRaises(ValueError):
                    source(encoded)
                with self.assertRaises(ValueError):
                    sp.TensorExpression(
                        reference.to_expression().replace(original, encoded)
                    )

    def test_named_builtin_factories_preserve_heads_and_patterns(self):
        for api, dimensions, indices, head in (
            ("dirac_gamma", (4,), ("i", "j", "mu"), "spenso::gamma"),
            ("color_f", (8,), ("a", "b", "c"), "spenso::f"),
            ("color_t", (8, 3), ("a", "i", "j"), "spenso::t"),
        ):
            with self.subTest(api=api):
                name = getattr(sp.TensorName, api)()
                self.assertEqual(name.to_expression().get_name(), head)
                factory = getattr(sp.TensorExpression, api)
                value = factory(*dimensions)(*indices)
                self.assertEqual(value.rank, len(indices))
                self.assertEqual(value.to_expression().get_name(), head)
                wildcards = [S(f"builtin_api::{index}_") for index in indices]
                dim_patterns = [
                    S(f"builtin_api::d{k}_") for k in range(len(dimensions))
                ]
                pattern = getattr(sp.TensorPattern, api)(*dim_patterns, *wildcards)
                self.assertEqual(value.to_expression().replace(pattern, E("7")), E("7"))
                with self.assertRaisesRegex(TypeError, f"TensorExpression\\.{api}"):
                    name()
                old_name = head.removeprefix("spenso::")
                for owner in (sp.TensorName, sp.TensorExpression, sp.TensorPattern):
                    self.assertNotIn(old_name, owner.__dict__)

        structure_constant = sp.TensorExpression.color_f(8)
        self.assertEqual(
            structure_constant("b", "a", "c").to_expression(),
            -structure_constant("a", "b", "c").to_expression(),
        )
        generator = sp.TensorExpression.color_t(8, 3)("a", "i", "j")
        self.assertEqual(
            generator.axes[2].representation, sp.Representation.cof(3).dual()
        )

    def test_transforms_preserve_interface_and_refresh_identity(self):
        a = self.A(self.x, self.rep("i"), self.rep("j"))
        replaced = a.replace(self.x, self.y)
        self.assertIsInstance(replaced, sp.TensorExpression)
        self.assertEqual(replaced.arguments, (self.y,))
        self.assertEqual(replaced.structure.axes, a.structure.axes)
        self.assertEqual(
            a.map(T().replace(self.x, self.y)),
            replaced,
        )
        self.assertEqual(
            a.replace_multiple([Replacement(self.x, self.y)]),
            replaced,
        )
        # Formal derivatives of tensor-valued functions remain outside the
        # strict tensor grammar; scalar-prefactor derivatives preserve its ports.
        with self.assertRaisesRegex(ValueError, "formal tensor derivatives"):
            (a * self.x).derivative(self.x)
        weighted = a * self.y
        derivative = weighted.derivative(self.y)
        self.assertEqual(derivative.structure.axes, a.structure.axes)
        self.assertEqual(derivative.to_expression(), a.to_expression())
        zero = a.derivative(S("api_operations::absent"))
        self.assertEqual(zero.rank, 2)
        self.assertEqual(zero.to_expression(), E("0"))
        with self.assertRaisesRegex(ValueError, "interface"):
            sp.TensorExpression(E("1"), structure=a.structure)
        library = sp.TensorLibrary()
        library.register(
            sp.Tensor.dense(self.A(self.y, self.rep, self.rep), [5.0, 6.0, 7.0, 8.0])
        )
        self.assertEqual(replaced.to_tensor(library)[:], [5.0, 6.0, 7.0, 8.0])

    def test_symbolica_rewrites_preserve_ordered_ports_and_reject_interface_changes(
        self,
    ):
        source = self.A(self.x, self.rep("j"), self.rep("i")).permute_axes([1, 0])
        expected = self.A(self.y, self.rep("j"), self.rep("i")).permute_axes([1, 0])
        for changed, expected_result in (
            (source.replace(self.x, self.y), expected),
            (source.replace(self.x, lambda _: self.y), expected),
            (source.replace_multiple([Replacement(self.x, self.y)]), expected),
            (source.map(T().replace(self.x, self.y)), expected),
            ((self.y * source).derivative(self.y), source),
        ):
            self.assertIsInstance(changed, sp.TensorExpression)
            self.assertEqual(changed.structure.axes, source.structure.axes)
            self.assertTrue(Expression.__eq__(changed, changed.to_expression()))
            self.assertEqual(changed, expected_result)
            self.assertEqual(changed.arguments, expected_result.arguments)

        for replacement in (E("1"), self.A(self.rep("k"), self.rep("i"))):
            for operation in (
                lambda rhs: source.replace(source.to_expression(), rhs),
                lambda rhs: source.replace_multiple([Replacement(source, rhs)]),
                lambda rhs: source.map(T().replace(source, rhs)),
            ):
                with (
                    self.subTest(replacement=replacement),
                    self.assertRaisesRegex(ValueError, "interface"),
                ):
                    operation(replacement)
        # An atomic data key follows its function head, so a zero no longer has
        # that key. Explicit result metadata is retained across these rewrites.
        self.assertEqual(
            source.replace(source, 0).structure.axes, source.structure.axes
        )
        self.assertIsNone(source.replace(source, 0).name)
        source = source.with_name(sp.TensorName("api_operations::RewriteMetadata")(7))
        for zero in (
            source.replace(source.to_expression(), E("0")),
            source.replace_multiple([Replacement(source, E("0"))]),
            source.map(T().replace(source, E("0"))),
            source.derivative(self.y),
        ):
            self.assertIsInstance(zero, sp.TensorExpression)
            self.assertEqual(zero.to_expression(), E("0"))
            self.assertEqual(zero.structure, source.structure)

    def test_reusable_tensor_rules_share_checked_replacement(self):
        source = self.A(self.rep("i"))
        result = sp.TensorName("api_operations::RuleResult")(self.rep("i"))
        lhs, rhs = source.to_expression(), result.to_expression()
        rule = sp.TensorRule(lhs, rhs)
        for _ in range(2):
            changed = source.replace(rule)
            self.assertEqual(changed.to_expression(), rhs)
            self.assertEqual(changed.structure.axes, source.structure.axes)
        self.assertEqual(source.replace(sp.TensorRule(lhs, rhs)), result)
        with self.assertRaises(TypeError):
            source.replace(rule, rhs_cache_size=0)
        self.assertEqual(source.to_expression(), lhs)
        zero = source.replace(sp.TensorRule(lhs, E("0")))
        self.assertEqual(zero.to_expression(), E("0"))
        self.assertEqual(zero.structure.axes, source.structure.axes)

        repeated = self.x * source + self.y * source
        for cache_size, expected_calls in [(0, 4), (10, 2)]:
            calls = []

            def replace(_, calls=calls):
                calls.append(None)
                return rhs

            rule = sp.TensorRule(lhs, replace, rhs_cache_size=cache_size)
            self.assertEqual(calls, [])
            for _ in range(2):
                self.assertEqual(
                    repeated.replace(rule), self.x * result + self.y * result
                )
            self.assertEqual(len(calls), expected_calls)

    def test_nested_scalar_expression_uses_the_existing_symbolica_evaluator(self):
        value = sp.TensorExpression(((self.x + 1) ** 3 + 2) ** 2)
        expected = ((self.x + 1) ** 3 + 2) ** 2
        self.assertEqual(value.to_expression(), expected)
        self.assertEqual(value.expand().to_expression(), expected.expand())
        evaluator = value.evaluator([self.x], iterations=1, n_cores=1)
        self.assertEqual(evaluator.evaluate([[2.0]])[0][0], 841.0)

    def test_symbolica_algebra_returns_tensors_with_exact_requested_forms(self):
        x, y = self.x, self.y
        tensor = self.A(self.rep("j"), self.rep("i")).permute_axes([1, 0])
        a = tensor.to_expression()
        for method, arguments, source, expected in (
            ("expand_num", (), 3 * (x + y) * tensor, (3 * x + 3 * y) * a),
            ("collect_factors", (), x * tensor + y * tensor, (x + y) * a),
            ("factor", (), (x**2 - 1) * tensor, (x - 1) * (x + 1) * a),
            ("collect", (x,), x * tensor + y * x * tensor, x * (a + y * a)),
            ("collect_num", (), (2 * x + 4 * y) * tensor, 2 * (x + 2 * y) * a),
            (
                "collect_by_coefficient",
                (),
                2 * x * tensor + 2 * y * tensor,
                2 * (x * a + y * a),
            ),
            ("collect_symbol", (x,), x * tensor + y * x * tensor, x * (a + y * a)),
            ("collect_horner", ([x, y],), (x**2 + x * y) * tensor, x * (x + y) * a),
            ("together", (), tensor / x + tensor / y, (x * a + y * a) / (x * y)),
            ("cancel", (), (x**2 - 1) / (x - 1) * tensor, (x + 1) * a),
            ("apart", (x,), tensor / (x * (x + 1)), a / x - a / (x + 1)),
        ):
            with self.subTest(method=method):
                source = source.with_name("api_operations::AlgebraResult")
                before = source.to_expression()
                result = getattr(source, method)(*arguments)
                self.assertIsInstance(result, sp.TensorExpression)
                self.assertEqual(result.to_expression(), expected)
                self.assertEqual(result.structure, source.structure)
                self.assertEqual(result.structure.axes, tensor.structure.axes)
                self.assertEqual(source.to_expression(), before)

    def test_symbolica_algebra_preserves_unresolved_ports_metadata_and_typed_zero(self):
        tensor = self.A(7, self.rep("i"), self.other)
        self.assertIsInstance(tensor.axes[1], sp.Representation)
        for source in (tensor, 0 * tensor):
            source = source.with_name(
                sp.TensorName("api_operations::AlgebraMetadata")(7)
            )
            self.assertEqual(source.arguments, (E("7"),))
            for method, arguments in (
                ("factor", ()),
                ("expand_num", ()),
                ("collect", (self.x, self.y)),
                ("collect_num", ()),
                ("collect_factors", ()),
                ("collect_by_coefficient", ()),
                ("collect_symbol", (self.x,)),
                ("collect_horner", ()),
                ("together", ()),
                ("cancel", ()),
                ("apart", (self.x, self.y)),
            ):
                with self.subTest(method=method, zero=not bool(source)):
                    result = getattr(source, method)(*arguments)
                    self.assertIsInstance(result, sp.TensorExpression)
                    self.assertEqual(result.to_expression(), source.to_expression())
                    self.assertEqual(result.structure, source.structure)
                    self.assertEqual(result.rank, 2)

        opened = self.A(self.other, self.rep)
        self.assertEqual(opened.axes, (self.other, self.rep))
        library = sp.TensorLibrary()
        library.register(
            sp.Tensor.dense(
                self.A(self.other, self.rep), [E(str(i)) for i in range(1, 7)]
            )
        )
        for source, method, coefficient in (
            (self.x * opened + self.y * opened, "collect_factors", self.x + self.y),
            ((self.x**2 - 1) * opened, "factor", self.x**2 - 1),
        ):
            with self.subTest(method=method, changed=True):
                result = getattr(source, method)()
                self.assertEqual(result.structure.axes, opened.structure.axes)
                # Keep the declared component order as well as its shape.
                for actual, component in zip(
                    result.to_tensor(library)[:], range(1, 7), strict=True
                ):
                    self.assertEqual(
                        (actual - coefficient * component).expand(), E("0")
                    )

    def test_collection_callbacks_check_the_actual_replacement_interface(self):
        tensor = self.A(self.rep("j"), self.rep("i")).permute_axes([1, 0])
        source = (self.x * tensor).with_name("api_operations::CallbackResult")
        for method in ("collect", "collect_symbol"):
            with self.subTest(method=method):
                mapped = getattr(source, method)(
                    self.x, coeff_map=lambda value: 2 * value
                )
                self.assertEqual(mapped.to_expression(), 2 * source.to_expression())
                self.assertEqual(mapped.structure, source.structure)
                zero = getattr(source, method)(self.x, coeff_map=lambda value: E("0"))
                self.assertEqual(zero.to_expression(), E("0"))
                self.assertEqual(zero.structure, source.structure)
                for replacement in (
                    E("1"),
                    self.A(self.rep("k"), self.rep("i")).to_expression(),
                ):
                    with (
                        self.subTest(replacement=str(replacement)),
                        self.assertRaises(ValueError),
                    ):
                        getattr(source, method)(
                            self.x,
                            coeff_map=lambda value, replacement=replacement: (
                                replacement
                            ),
                        )
                with self.assertRaises(ValueError):
                    getattr(source, method)(
                        self.x,
                        key_map=lambda key: (
                            key * self.A(self.other("mu")).to_expression()
                        ),
                    )

        opened = self.A(self.other, self.rep)
        self.assertEqual(opened.axes, (self.other, self.rep))
        source = self.x * opened
        mapped = source.collect(self.x, coeff_map=lambda coefficient: 2 * coefficient)
        self.assertEqual(mapped.structure, source.structure)
        self.assertEqual(mapped.to_expression(), 2 * source.to_expression())

    def test_unchanged_algebra_retains_completed_reduction_status(self):
        source = sp.TensorExpression(7).simplify_algebra()
        self.assertEqual(source.reduction_status, sp.ReductionStatus.Complete)
        result = source.expand_num().collect_factors().factor().cancel()
        self.assertEqual(result, source)
        self.assertEqual(result.reduction_status, sp.ReductionStatus.Complete)

    def test_read_only_symbolica_inspection_preserves_tensor_value(self):
        tensor = self.A(self.x, self.rep("i"))
        value = (self.x + self.y) * tensor
        expression, structure = value.to_expression(), value.structure
        self.assertTrue(value.contains(self.x))
        self.assertTrue(self.x in value)
        self.assertFalse(S("api_operations::absent") in value)
        for method, arguments in (
            ("get_all_symbols", ()),
            ("get_all_symbols", (False,)),
            ("get_all_indeterminates", ()),
            ("get_all_indeterminates", (False,)),
            ("get_type", ()),
            ("get_byte_size", ()),
        ):
            with self.subTest(method=method, arguments=arguments):
                self.assertEqual(
                    getattr(value, method)(*arguments),
                    getattr(expression, method)(*arguments),
                )
        self.assertEqual(tensor.get_head(), self.A.to_expression())
        self.assertEqual(tensor.get_name(), "api_operations::A")
        self.assertTrue((0 * tensor).is_zero())
        self.assertTrue(sp.TensorExpression(1).is_one())
        self.assertFalse(value.is_expanded())
        self.assertTrue(value.expand().is_expanded(self.x))
        self.assertEqual(value.nterms(), 1)
        self.assertEqual(value.expand().nterms(), 2)
        self.assertEqual(value.to_expression(), expression)
        self.assertEqual(value.structure, structure)

    def test_scalar_tensor_denominator_accepts_an_open_expression(self):
        tensor = self.A(self.rep("j"), self.rep("i"))
        denominator = sp.TensorExpression(self.x)
        # A plain Expression on the left dispatches to TensorExpression.__rtruediv__.
        for numerator in (tensor.to_expression(), tensor, 0 * tensor):
            expected = sp.TensorExpression(numerator)
            quotient = denominator.__rtruediv__(numerator)
            self.assertIsInstance(quotient, sp.TensorExpression)
            self.assertEqual(quotient.structure.axes, expected.structure.axes)
            self.assertEqual(
                quotient.to_expression(), expected.to_expression() / self.x
            )
        quotient = tensor.to_expression() / denominator
        self.assertIsInstance(quotient, sp.TensorExpression)
        self.assertEqual(quotient, tensor / self.x)
        self.assertEqual((self.y / denominator).to_expression(), self.y / self.x)
        with self.assertRaisesRegex(ValueError, "denominator"):
            self.x / tensor

    def test_noop_transforms_preserve_complete_metadata_in_fresh_objects(self):
        atomic = self.A(7, self.rep("i"), self.other("mu"))
        opened = self.A(7, self.rep, self.other)
        composite = ((self.x + self.y) * atomic).with_name(
            "api_operations::NamedComposite"
        )
        zero = sp.TensorExpression(
            atomic.to_expression().derivative(S("api_operations::absent")),
            structure=atomic.structure,
        ).with_name("api_operations::NamedZero")
        permuted = atomic.permute_axes([1, 0])
        self.assertEqual(atomic.arguments, (E("7"),))
        self.assertEqual(zero.to_expression(), E("0"))
        self.assertEqual(zero.rank, 2)
        self.assertEqual(permuted.to_expression(), atomic.to_expression())
        self.assertEqual(permuted.structure, atomic.structure)

        for label, source in [
            ("atomic", atomic),
            ("open", opened),
            ("composite", composite),
            ("zero", zero),
            ("permuted", permuted),
        ]:
            expression, structure = source.to_expression(), source.structure
            for method in ["to_dots", "simplify_algebra"]:
                with self.subTest(case=label, method=method):
                    result = (
                        getattr(source, method)(gamma=True, epsilon=True)
                        if method == "simplify_algebra"
                        else source.to_dots()
                    )
                    self.assertIsNot(result, source)
                    self.assertIsInstance(result, sp.TensorExpression)
                    self.assertEqual(result.to_expression(), expression)
                    self.assertEqual(result.structure, structure)
                    rerun = (
                        getattr(result, method)(gamma=True, epsilon=True)
                        if method == "simplify_algebra"
                        else result.to_dots()
                    )
                    self.assertIsNot(rerun, result)
                    self.assertIsInstance(rerun, sp.TensorExpression)
                    self.assertEqual(rerun.to_expression(), expression)
                    self.assertEqual(rerun.structure, structure)
                    self.assertEqual(source.structure, structure)

    def test_simultaneous_renaming_and_explicit_reindex(self):
        a = self.A(self.rep("i"), self.rep("j"))
        swapped = a.rename_indices({"i": "j", "j": "i"})
        self.assertEqual(swapped.structure, a.structure)
        self.assertEqual(
            a.rename_indices({self.rep("i"): "k"}), a.reindex("k", sp.AUTO)
        )
        with self.assertRaises(KeyError):
            a.rename_indices({"missing": "k"})
        with self.assertRaises(ValueError):
            a.rename_indices({"i": "j"})
        self.assertEqual(a.reindex("k", "k").rank, 0)
        tensor = self.matrix().reindex("i", "j").result_tensor()
        renamed = tensor.rename_indices({"i": "j", "j": "i"})
        self.assertEqual(renamed[:], tensor[:])
        self.assertEqual(renamed.structure.axes, swapped.structure.axes)

    def test_axis_permutation_retains_data_and_library_identity(self):
        for reps, shape in [
            ((self.rep, self.rep), (2, 2)),
            ((self.rep, self.other), (2, 3)),
        ]:
            with self.subTest(shape=shape):
                expression = self.A(*reps)
                values = [float(i) for i in range(np.prod(shape))]
                tensor = sp.Tensor.dense(expression, values)
                expected = np.array(values).reshape(shape).T
                permuted = tensor.permute_axes([1, 0])
                self.assertEqual(permuted.shape, shape[::-1])
                np.testing.assert_array_equal(permuted.to_numpy(), expected)
                self.assertEqual(permuted.permute_axes([1, 0])[:], values)
                library = sp.TensorLibrary()
                library.register(tensor)
                expression_permuted = expression.permute_axes([1, 0])
                np.testing.assert_array_equal(
                    expression_permuted.to_tensor(library).to_numpy(), expected
                )
                network = expression.to_network(library=library).permute_axes([1, 0])
                np.testing.assert_array_equal(
                    network.to_tensor(library).to_numpy(), expected
                )
                indexed = expression("i", "j")
                self.assertEqual(
                    indexed.permute_axes([1, 0]).structure, indexed.structure
                )
                self.assertEqual(indexed.permute_axes([1, 0]).axes, indexed.axes[::-1])
                self.assertEqual(tensor[:], values)
        with self.assertRaises(ValueError):
            self.matrix().permute_axes([0, 0])

    def test_interning_is_shared_by_all_indexing_entry_points(self):
        payload = S("api_operations::edge")(7, 2)
        a = self.A(self.rep)
        tensor = sp.Tensor.dense(a, [1.0, 2.0])
        label = S("api_operations::raw_index")
        raw = a(label).to_expression().replace(label, payload)
        for mode in ("indices", "flattened"):
            expected = a(payload, intern=mode)
            self.assertEqual(sp.TensorExpression(raw, intern=mode), expected)
            for value in (a, tensor, a.to_network()):
                for indexed in (
                    value(payload, intern=mode),
                    value.index(payload, intern=mode),
                    value.reindex(payload, intern=mode),
                    value("i").rename_indices({"i": payload}, intern=mode),
                ):
                    expression = (
                        indexed
                        if isinstance(indexed, sp.TensorExpression)
                        else indexed.expression()
                    )
                    self.assertEqual(
                        expression.to_expression(), expected.to_expression()
                    )
                    self.assertEqual(expression.structure.axes, expected.structure.axes)
            self.assertEqual(sp.TensorExpression(expected, intern=mode), expected)
        self.assertNotEqual(
            a(payload, intern="indices"), a(payload, intern="flattened")
        )
        with self.assertRaises(ValueError):
            a(payload)

    def test_intern_defaults_validation_and_retired_cooking_surface(self):
        a = self.A(self.rep("i"))
        for source in (a, a.to_expression()):
            self.assertEqual(sp.TensorExpression(source), a)
            self.assertEqual(sp.TensorExpression(source, intern=None), a)
            self.assertEqual(
                sp.TensorExpression(source, structure=a.structure, intern=None), a
            )
        for invalid, error in (
            ("invalid", ValueError),
            (True, TypeError),
            (1, TypeError),
        ):
            with self.assertRaises(error):
                sp.TensorExpression(a, intern=invalid)
            with self.assertRaises(error):
                a.reindex("j", intern=invalid)
            # Reject invalid policies even when there are no indices to traverse.
            with self.assertRaises(error):
                sp.TensorExpression(1).index(intern=invalid)
            with self.assertRaises(error):
                a.rename_indices({}, intern=invalid)
        for name in ("CookSettings", "CookMode", "CookSourceFilter", "CookTagFilter"):
            self.assertFalse(hasattr(sp, name))
        self.assertFalse(
            any(name.startswith("cook") for name in dir(sp.TensorExpression))
        )
        with self.assertRaises(TypeError):
            sp.TensorExpression(a, cook_indices=None)

    def test_component_slices_and_independent_storage_conversions(self):
        tensor = self.matrix()
        self.assertEqual(tensor[::-1], [4.0, 3.0, 2.0, 1.0])
        self.assertEqual(tensor[-1], 4.0)
        self.assertEqual(tensor[:, -1], [2.0, 4.0])
        self.assertEqual(tensor[::-1, ::-1], [[4.0, 3.0], [2.0, 1.0]])
        self.assertEqual(tensor[0:0, :], [])
        with self.assertRaises(IndexError):
            _ = tensor[-5]
        with self.assertRaises(IndexError):
            _ = tensor[:, :, :]
        sparse = tensor.to_sparse()
        self.assertEqual(tensor.storage, "dense")
        self.assertEqual(sparse.storage, "sparse")
        sparse[-1, -1] = 9.0
        self.assertEqual(tensor[-1], 4.0)
        self.assertEqual(sparse[-1], 9.0)
        for independent in (tensor.copy(), copy.copy(tensor), sparse.to_dense()):
            independent[0] = 13.0
        self.assertEqual(tensor[0], 1.0)
        self.assertEqual(sparse[0], 1.0)

    def test_sparse_maps_materialize_changed_zero_before_contraction(self):
        tensor = sp.Tensor.sparse(self.A(self.rep), float)
        tensor[0] = 2.0
        scaled = tensor.map_components(lambda value: 3 * value)
        self.assertEqual(scaled.storage, "sparse")
        shifted = tensor.map_components(lambda value: value + 1)
        self.assertEqual(shifted.storage, "dense")
        self.assertEqual(shifted[:], [3.0, 1.0])
        self.assertEqual((shifted("i") * shifted("i")).to_tensor()[0], 10.0)
        symbolic = tensor.map_components(lambda value: self.x * value, dtype=Expression)
        self.assertIs(symbolic.dtype, Expression)
        self.assertEqual(symbolic[:], [2.0 * self.x, E("0")])
        complex_tensor = tensor.map_components(lambda value: value + 1j, dtype=complex)
        self.assertIs(complex_tensor.dtype, complex)
        self.assertEqual(
            complex_tensor.map_components(lambda value: value.conjugate())[:],
            [2 - 1j, -1j],
        )
        self.assertEqual(tensor[:], [2.0, 0.0])
        calls = []
        with self.assertRaises(TypeError):
            tensor.map_components(lambda value: calls.append(value), dtype=str)
        self.assertEqual(calls, [])

    def test_numpy_full_shape_order_dtype_and_copy(self):
        expression = self.A(self.rep, self.other)
        source = np.arange(6).reshape(3, 2).T
        tensor = sp.Tensor.from_numpy(expression, source)
        np.testing.assert_array_equal(tensor.to_numpy(), source)
        self.assertIs(tensor.dtype, float)
        source[0, 0] = 99
        self.assertEqual(tensor[0], 0.0)
        result = tensor.to_numpy()
        result[0, 0] = 77
        self.assertEqual(tensor[0], 0.0)
        complex_tensor = sp.Tensor.from_numpy(expression, source + 1j)
        self.assertIs(complex_tensor.dtype, complex)
        np.testing.assert_array_equal(complex_tensor.to_numpy(), source + 1j)
        with self.assertRaisesRegex(ValueError, "shape"):
            sp.Tensor.from_numpy(expression, np.zeros((3, 2)))
        with self.assertRaises(TypeError):
            sp.Tensor.from_numpy(expression, np.zeros((2, 3), dtype=object))
        with self.assertRaises(TypeError):
            self.matrix(symbolic=True).to_numpy()

    def test_network_progress_and_execution_copies(self):
        tensor = self.matrix()
        network = tensor("i", "j") * tensor("j", "k") * tensor("k", "l")
        status = network.status
        self.assertFalse(status.complete)
        self.assertEqual(status.contractions, 2)
        self.assertTrue(status.ready_operations)
        stepped = network.step()
        self.assertEqual(network.status.contractions, 2)
        self.assertLess(stepped.status.contractions, 2)
        np.testing.assert_array_equal(
            network.to_tensor().to_numpy(), np.linalg.matrix_power(tensor.to_numpy(), 3)
        )
        self.assertFalse(network.status.complete)
        stepped.execute()
        self.assertTrue(stepped.status.complete)
        self.assertEqual(stepped.result_tensor()[:], network.to_tensor()[:])
        self.assertFalse(status.complete)
        with self.assertRaises(AttributeError):
            status.complete = True
        trace = tensor.reindex("i", "i")
        self.assertFalse(trace.status.complete)
        self.assertTrue(trace.step().status.complete)
        self.assertEqual(trace.to_tensor()[0], 5.0)

    def test_rank_three_permutation_and_broadcast_execution(self):
        source = np.arange(12.0).reshape(2, 3, 2)
        tensor = sp.Tensor.from_numpy(self.A(self.rep, self.other, self.rep), source)
        permuted = tensor.permute_axes([2, 0, 1])
        np.testing.assert_array_equal(permuted.to_numpy(), source.transpose(2, 0, 1))
        square = sp.BroadcastFunction("api_operations::square")
        functions = sp.TensorFunctionLibrary()
        functions.register(square, lambda value: value * value)
        network = square(permuted)
        np.testing.assert_array_equal(
            network.to_tensor(function_library=functions).to_numpy(),
            source.transpose(2, 0, 1) ** 2,
        )
        self.assertFalse(network.status.complete)
        library = sp.TensorLibrary()
        library.register(tensor)
        expression = square(tensor.expression())
        np.testing.assert_array_equal(
            expression.to_tensor(library, function_library=functions).to_numpy(),
            source**2,
        )

    def test_renaming_rejects_dummy_capture_in_expressions_and_networks(self):
        tensor = self.matrix()
        expression = tensor.expression()("i", "j") * tensor.expression()("j", "k")
        network = tensor("i", "j") * tensor("j", "k")
        for value in (expression, network):
            with self.assertRaises(ValueError):
                value.rename_indices({"i": "j"})
            renamed = value.rename_indices({"i": "k", "k": "i"})
            self.assertEqual(renamed.structure, value.structure)

    def test_simplification_presets_and_expansion_control(self):
        a = self.A(self.rep("i"))
        expression = (self.x + 1) * (self.y + 1) * a
        self.assertEqual(expression.simplify_algebra(), expression)
        self.assertEqual(expression.simplify_algebra().expand(), expression.expand())
        metric = sp.TensorExpression.g(self.rep)("i", "j") * a
        self.assertEqual(
            metric.simplify_algebra(contract="none").expand(), metric.expand()
        )
        self.assertEqual(metric.simplify_algebra(), self.A(self.rep("j")))
        gamma = sp.TensorExpression.dirac_gamma(4)
        chain = gamma("a", "b", "mu") * gamma("b", "a", "nu")
        simplified = chain.simplify_algebra(gamma=True, color=True, epsilon=True)
        self.assertEqual(
            simplified.expand(),
            simplified.simplify_algebra(gamma=True, color=True, epsilon=True).expand(),
        )
        self.assertEqual(simplified.rank, 2)
        d = S("api_operations::D")
        symbolic_metric = sp.TensorExpression.g(sp.Representation.mink(d))("mu", "mu")
        self.assertEqual(symbolic_metric.simplify_algebra().expand().to_expression(), d)
        capped = chain.simplify_algebra(gamma=True, max_steps_per_domain=0)
        self.assertEqual(capped.reduction_status, sp.ReductionStatus.Capped)
        self.assertEqual(capped, chain)
        bounded = chain.simplify_algebra(gamma=True, max_steps_per_domain=1)
        self.assertIn(
            bounded.reduction_status,
            (sp.ReductionStatus.Complete, sp.ReductionStatus.Capped),
        )
        self.assertEqual(
            bounded.simplify_algebra(gamma=True, color=True, epsilon=True).expand(),
            simplified.expand(),
        )

    def test_generated_coefficient_collection_keeps_scalar_spectators(self):
        lorentz = sp.Representation.mink(4)
        gamma = sp.TensorExpression.dirac_gamma(4)
        p = sp.TensorName.vector("api_operations::coefficient_p")(lorentz)
        q = sp.TensorName.vector("api_operations::coefficient_q")(lorentz)
        spectator = (self.x + self.y) ** 12
        source = (
            spectator
            * (p("rho") + q("rho"))
            * gamma("a", "b", "rho")
            * gamma("b", "c", "mu")
            * (p("sigma") + q("sigma"))
            * gamma("c", "d", "sigma")
            * gamma("d", "a", "mu")
        ).with_name("api_operations::CollectedTrace")
        options = {"color": False, "contract": "dots"}
        collected = source.simplify_algebra(**options)
        nested = source.simplify_algebra(**options, collect_coefficients=False)
        self.assertEqual(
            collected, source.simplify_algebra(**options, collect_coefficients=True)
        )
        # Independent Clifford identity: Tr(r/ gamma_mu r/ gamma^mu) = -8 r^2.
        expected = -8 * (sp.dot(p, p) + 2 * sp.dot(p, q) + sp.dot(q, q))
        for result in (collected, nested):
            self.assertIsInstance(result, sp.TensorExpression)
            self.assertEqual(result.structure, source.structure)
            self.assertEqual(result.reduction_status, sp.ReductionStatus.Complete)
            self.assertTrue(bool(result.contains(spectator)))
            coefficient = result.to_expression().replace(spectator, E("1"))
            self.assertEqual(
                (coefficient - expected.to_expression()).collect_factors().factor(),
                E("0"),
            )
        self.assertEqual(collected.simplify_algebra(**options), collected)
        untouched = ((self.x + 1) * (self.y + 1) * self.A(self.rep("i"))).with_name(
            "api_operations::UnrelatedCoefficient"
        )
        for collect in (False, True):
            self.assertEqual(
                untouched.simplify_algebra(
                    gamma=False, color=False, collect_coefficients=collect
                ),
                untouched,
            )

    def test_interpreted_and_compiled_batch_contract(self):
        for values in (
            [self.x, self.y, self.x * self.y, E("0")],
            [self.x * E("1i"), self.y, E("0"), self.x],
        ):
            tensor = sp.Tensor.dense(self.A(self.rep, self.rep), values).to_sparse()
            evaluator = tensor.evaluator(
                {}, {}, [self.y, self.x], iterations=1, n_cores=1
            )
            self.assertEqual(evaluator.parameters, [self.y, self.x])
            self.assertEqual(evaluator.input_size, 2)
            self.assertEqual(evaluator.output_shape, (2, 2))
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                compiled = evaluator.compile(
                    "api_operations_eval",
                    str(root / "eval.cpp"),
                    str(root / "eval.so"),
                    inline_asm="none",
                    optimization_level=0,
                )
                self.assertEqual(compiled.parameters, evaluator.parameters)
                self.assertEqual(compiled.supports_real, evaluator.supports_real)
                for engine in (evaluator, compiled):
                    self.assertEqual(engine.evaluate_complex([]), [])
                    for batch in ([[1.0]], [[1.0, 2.0], [1.0, 2.0, 3.0]]):
                        with self.assertRaisesRegex(ValueError, "input row"):
                            engine.evaluate(batch)
                        with self.assertRaisesRegex(ValueError, "input row"):
                            engine.evaluate_complex(batch)
                    if engine.supports_real:
                        self.assertEqual(
                            engine.evaluate([[2.0, 3.0]])[0][:], [3.0, 2.0, 6.0, 0.0]
                        )
                    else:
                        with self.assertRaisesRegex(ValueError, "complex coefficients"):
                            engine.evaluate([[2.0, 3.0]])
                inputs = [[2 + 1j, 3 - 2j], [0j, 0j]]
                for a, b in zip(
                    evaluator.evaluate_complex(inputs),
                    compiled.evaluate_complex(inputs),
                    strict=True,
                ):
                    self.assertEqual(a[:], b[:])
                    self.assertEqual(a.structure, b.structure)
        constant = sp.Tensor.dense(self.A(self.rep), [1.0, 2.0]).evaluator(
            {}, {}, [], iterations=1, n_cores=1
        )
        self.assertEqual(constant.evaluate([[]])[0][:], [1.0, 2.0])


if __name__ == "__main__":
    unittest.main()
