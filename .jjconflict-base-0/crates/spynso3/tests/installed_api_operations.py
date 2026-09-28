"""Tensor API contracts against the installed extension and native compiler."""

import copy
import tempfile
import unittest
from pathlib import Path

import numpy as np
from symbolica import E, Expression, Replacement, S, T
from symbolica.community import spenso as sp


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

    def test_transforms_preserve_interface_and_refresh_identity(self):
        a = self.A(self.x, self.rep("i"), self.rep("j"))
        replaced = sp.TensorExpression(
            a.to_expression().replace(self.x, self.y), structure=a.structure
        )
        self.assertIsInstance(replaced, sp.TensorExpression)
        self.assertEqual(replaced.structure.arguments, (self.y,))
        self.assertEqual(replaced.structure.slots, a.structure.slots)
        self.assertEqual(
            sp.TensorExpression(
                a.to_expression().map(T().replace(self.x, self.y)),
                structure=a.structure,
            ),
            replaced,
        )
        self.assertEqual(
            sp.TensorExpression(
                a.to_expression().replace_multiple([Replacement(self.x, self.y)]),
                structure=a.structure,
            ),
            replaced,
        )
        # Formal derivatives of tensor-valued functions remain outside the
        # strict tensor grammar; scalar-prefactor derivatives preserve its ports.
        with self.assertRaisesRegex(ValueError, "not tagged as a tensor"):
            sp.TensorExpression(
                (a.to_expression() * self.x).derivative(self.x),
                structure=a.structure,
            )
        weighted = a * self.y
        derivative = sp.TensorExpression(
            weighted.to_expression().derivative(self.y),
            structure=weighted.structure,
        )
        self.assertEqual(derivative.structure.slots, a.structure.slots)
        self.assertEqual(derivative.to_expression(), a.to_expression())
        zero = sp.TensorExpression(
            a.to_expression().derivative(S("api_operations::absent")),
            structure=a.structure,
        )
        self.assertEqual(zero.rank, 2)
        self.assertEqual(zero.to_expression(), E("0"))
        with self.assertRaisesRegex(ValueError, "interface"):
            sp.TensorExpression(E("1"), structure=a.structure)
        library = sp.TensorLibrary()
        library.register(
            sp.Tensor.dense(self.A(self.y, self.rep, self.rep), [5.0, 6.0, 7.0, 8.0])
        )
        self.assertEqual(replaced.to_tensor(library)[:], [5.0, 6.0, 7.0, 8.0])

    def test_reusable_tensor_rules_share_checked_replacement(self):
        source = self.A(self.rep("i"))
        result = sp.TensorName("api_operations::RuleResult")(self.rep("i"))
        lhs, rhs = source.to_expression(), result.to_expression()
        rule = sp.TensorRule(lhs, rhs)
        for _ in range(2):
            changed = source.replace(rule)
            self.assertEqual(changed.to_expression(), rhs)
            self.assertEqual(changed.structure.slots, source.structure.slots)
        self.assertEqual(source.replace(sp.TensorRule(lhs, rhs)), result)
        with self.assertRaises(TypeError):
            source.replace(rule, rhs_cache_size=0)
        self.assertEqual(source.to_expression(), lhs)
        zero = source.replace(sp.TensorRule(lhs, E("0")))
        self.assertEqual(zero.to_expression(), E("0"))
        self.assertEqual(zero.structure.slots, source.structure.slots)

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

    def test_typed_aliases_retain_order_metadata_and_callback_errors(self):
        source = self.A(self.x, self.rep("j"), self.rep("i"))
        value = sp.AliasedTensorExpression.from_expression(source)
        self.assertEqual(value.root.structure.slots, source.structure.slots)
        self.assertEqual(value.to_expression(), source)
        handle, body = value.aliases[0]
        self.assertEqual(body.structure.slots, source.structure.slots)
        self.assertEqual(body.structure.arguments, (self.x,))
        self.assertEqual(handle.structure.slots, source.structure.slots)
        calls = []

        def identity(body):
            calls.append(body.structure.arguments)
            return body

        mapped = value.map_aliases(identity)
        self.assertEqual(calls, [(self.x,)])
        self.assertEqual(mapped.to_expression(), source)
        self.assertEqual(mapped.aliases[0][1].structure.arguments, (self.x,))
        zero = value.map_aliases(lambda body: 0 * body)
        self.assertEqual(zero.to_expression().to_expression(), E("0"))
        self.assertEqual(zero.to_expression().structure.slots, source.structure.slots)

        def failing(_):
            raise LookupError("alias callback sentinel")

        with self.assertRaisesRegex(LookupError, "alias callback sentinel"):
            value.map_aliases(failing)
        with self.assertRaisesRegex(ValueError, "interface"):
            value.map_aliases(lambda _: sp.TensorExpression(E("1")))

    def test_nested_scalar_aliases_use_the_existing_symbolica_evaluator(self):
        inner = sp.AliasedTensorExpression.from_expression(
            sp.TensorExpression((self.x + 1) ** 3)
        )
        outer = sp.TensorExpression(S("api_operations::outer_alias"))
        body = (inner.root + 2) ** 2
        value = sp.AliasedTensorExpression(outer, inner.aliases + [(outer, body)])
        expected = ((self.x + 1) ** 3 + 2) ** 2
        self.assertEqual(value.to_expression().to_expression(), expected)
        self.assertEqual(value.expand().to_expression(), expected.expand())
        evaluator = value.evaluator([self.x], iterations=1, n_cores=1)
        self.assertEqual(evaluator.evaluate([[2.0]])[0][0], 841.0)

    def test_scalar_tensor_denominator_accepts_an_open_expression(self):
        tensor = self.A(self.rep("j"), self.rep("i"))
        denominator = sp.TensorExpression(self.x)
        # A plain Expression on the left dispatches to TensorExpression.__rtruediv__.
        for numerator in (tensor.to_expression(), tensor, 0 * tensor):
            expected = sp.TensorExpression(numerator)
            quotient = denominator.__rtruediv__(numerator)
            self.assertIsInstance(quotient, sp.TensorExpression)
            self.assertEqual(quotient.structure.slots, expected.structure.slots)
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
        self.assertEqual(atomic.structure.arguments, (E("7"),))
        self.assertEqual(zero.to_expression(), E("0"))
        self.assertEqual(zero.rank, 2)
        self.assertEqual(permuted.to_expression(), atomic.to_expression())
        self.assertEqual(permuted.structure.slots, atomic.structure.slots[::-1])

        for label, source in [
            ("atomic", atomic),
            ("open", opened),
            ("composite", composite),
            ("zero", zero),
            ("permuted", permuted),
        ]:
            expression, structure = source.to_expression(), source.structure
            for method in ["to_dots", "simplify_gamma"]:
                with self.subTest(case=label, method=method):
                    result = getattr(source, method)()
                    self.assertIsNot(result, source)
                    typed = (
                        result.to_expression() if method == "simplify_gamma" else result
                    )
                    self.assertEqual(typed.to_expression(), expression)
                    self.assertEqual(typed.structure, structure)
                    rerun = getattr(result, method)()
                    self.assertIsNot(rerun, result)
                    typed_rerun = (
                        rerun.to_expression() if method == "simplify_gamma" else rerun
                    )
                    self.assertEqual(typed_rerun.to_expression(), expression)
                    self.assertEqual(typed_rerun.structure, structure)
                    self.assertEqual(source.structure, structure)

    def test_simultaneous_renaming_and_explicit_reindex(self):
        a = self.A(self.rep("i"), self.rep("j"))
        swapped = a.rename_indices({"i": "j", "j": "i"})
        self.assertEqual(swapped.structure.slots, a.structure.slots[::-1])
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
        self.assertEqual(renamed.structure.slots, swapped.structure.slots)

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
                self.assertEqual(permuted.structure.shape, shape[::-1])
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
                self.assertEqual(
                    expression.structure.permute_axes([1, 0]).shape, shape[::-1]
                )
                self.assertEqual(tensor[:], values)
        with self.assertRaises(ValueError):
            self.matrix().permute_axes([0, 0])

    def test_cooking_is_shared_by_all_indexing_entry_points(self):
        payload = S("api_operations::edge")(7, 2)
        settings = sp.CookSettings.reversible()
        a = self.A(self.rep)
        tensor = sp.Tensor.dense(a, [1.0, 2.0])
        expected = a(payload, cook_indices=settings).to_expression()
        for value in (a, tensor, a.to_network()):
            indexed = value(payload, cook_indices=settings)
            expression = (
                indexed
                if isinstance(indexed, sp.TensorExpression)
                else indexed.expression()
            )
            self.assertEqual(expression.to_expression(), expected)
            decoded_index = sp.TensorExpression(
                expression.structure.slots[0].index
            ).uncook(settings)
            self.assertEqual(decoded_index.to_expression(), payload)
        with self.assertRaises(ValueError):
            a(payload)
        with self.assertRaises(TypeError):
            a(payload, cook_indices=True)

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
            self.assertEqual(renamed.structure.slots, value.structure.slots[::-1])

    def test_simplification_presets_and_expansion_control(self):
        a = self.A(self.rep("i"))
        expression = (self.x + 1) * (self.y + 1) * a
        self.assertEqual(expression.simplify().to_expression(), expression)
        self.assertEqual(expression.simplify().expand(), expression.expand())
        metric = sp.TensorExpression.g(self.rep)("i", "j") * a
        self.assertEqual(metric.simplify().expand(), metric.contract().expand())
        gamma = sp.TensorExpression.gamma(4)
        chain = gamma("a", "b", "mu") * gamma("b", "a", "nu")
        simplified = chain.simplify(sp.SimplifySettings.hep())
        self.assertEqual(
            simplified.expand(), simplified.simplify(sp.SimplifySettings.hep()).expand()
        )
        self.assertEqual(simplified.root.rank, 2)
        d = S("api_operations::D")
        symbolic_metric = sp.TensorExpression.g(sp.Representation.mink(d))("mu", "mu")
        self.assertEqual(symbolic_metric.simplify().expand().to_expression(), d)
        with self.assertRaises(ValueError):
            sp.SimplifySettings(max_passes=0)
        with self.assertRaisesRegex(ValueError, "stabilize"):
            metric.simplify(sp.SimplifySettings(max_passes=1))

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
