"""Tensor API contracts against the installed extension and native compiler."""

import copy
import tempfile
import unittest
from pathlib import Path

import numpy as np
from symbolica import E, Expression, S, T
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
        replaced = a.replace(self.x, self.y)
        self.assertIsInstance(replaced, sp.TensorExpression)
        self.assertEqual(replaced.structure.arguments, (self.y,))
        self.assertEqual(replaced.structure.slots, a.structure.slots)
        self.assertEqual(a.map(T().replace(self.x, self.y)), replaced)
        self.assertEqual(a.replace_multiple([(self.x, self.y)]), replaced)
        weighted = a * self.x
        self.assertEqual(weighted.derivative(self.x).structure.slots, a.structure.slots)
        zero = a.derivative(S("api_operations::absent"))
        self.assertEqual(zero.rank, 2)
        self.assertEqual(zero.to_expression(), E("0"))
        with self.assertRaisesRegex(ValueError, "interface"):
            a.replace(a, 1)
        library = sp.TensorLibrary()
        library.register(sp.Tensor.dense(replaced, [5.0, 6.0, 7.0, 8.0]))
        self.assertEqual(replaced.to_tensor(library)[:], [5.0, 6.0, 7.0, 8.0])

    def test_noop_transforms_preserve_complete_metadata_in_fresh_objects(self):
        atomic = self.A(7, self.rep("i"), self.other("mu"))
        opened = self.A(7, self.rep, self.other)
        composite = ((self.x + self.y) * atomic).with_name(
            "api_operations::NamedComposite"
        )
        zero = atomic.derivative(S("api_operations::absent")).with_name(
            "api_operations::NamedZero"
        )
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
            for method in ["normalize_dots", "simplify_gamma"]:
                with self.subTest(case=label, method=method):
                    result = getattr(source, method)()
                    self.assertIsNot(result, source)
                    self.assertEqual(result.to_expression(), expression)
                    self.assertEqual(result.structure, structure)
                    rerun = getattr(result, method)()
                    self.assertIsNot(rerun, result)
                    self.assertEqual(rerun.to_expression(), expression)
                    self.assertEqual(rerun.structure, structure)
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
        expected = a(payload, cook_indices=settings).wrap_indices().to_expression()
        self.assertIn("edge", str(expected))
        for value in (a, tensor, a.to_network()):
            indexed = value(payload, cook_indices=settings)
            expression = (
                indexed
                if isinstance(indexed, sp.TensorExpression)
                else indexed.expression()
            )
            self.assertEqual(expression.wrap_indices().to_expression(), expected)
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
        self.assertEqual(symbolic[:], [2 * self.x, E("0")])
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
        self.assertEqual(expression.simplify(), expression)
        self.assertEqual(
            expression.simplify(sp.SimplifySettings(expand=True)), expression.expand()
        )
        metric = sp.TensorExpression.g(self.rep)("i", "j") * a
        self.assertEqual(metric.simplify(), metric.simplify_metrics())
        gamma = sp.TensorExpression.gamma(4)
        chain = gamma("a", "b", "mu") * gamma("b", "a", "nu")
        simplified = chain.simplify(sp.SimplifySettings.hep())
        self.assertEqual(simplified, simplified.simplify(sp.SimplifySettings.hep()))
        self.assertEqual(simplified.rank, 2)
        d = S("api_operations::D")
        symbolic_metric = sp.TensorExpression.g(sp.Representation.mink(d))("mu", "mu")
        self.assertEqual(symbolic_metric.simplify().to_expression(), d)
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
