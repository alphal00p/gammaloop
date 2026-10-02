"""Stored library data and signature lookup against the installed extension."""

import unittest

from symbolica import E, Expression, S
from symbolica.community.tensor import (
    Representation,
    Tensor,
    TensorExpression,
    TensorLibrary,
    TensorName,
)


class TensorLibraryTests(unittest.TestCase):
    def test_full_signatures_disambiguate_extra_arguments_and_namespaces(self):
        rep = Representation.euc(2)
        x = S("library_tests::x")
        name = TensorName("library_tests::A")
        signatures = [
            name(x, 7, rep),
            name(x, 8, rep),
            TensorName("library_other::A")(x, 7, rep),
        ]
        library = TensorLibrary()
        for number, signature in enumerate(signatures):
            library.register(Tensor.dense(signature, [float(number), 10.0]))
        for number, signature in enumerate(signatures):
            stored = library[signature]
            self.assertIsInstance(stored, Tensor)
            self.assertEqual(stored[:], [float(number), 10.0])
            self.assertEqual(stored.structure, signature.structure)
            self.assertEqual(
                stored.expression()("i").to_expression(),
                signature("i").to_expression(),
            )
        with self.assertRaisesRegex(KeyError, "ambiguous"):
            library[name]
        self.assertEqual(library["library_other::A"][0], 2.0)
        self.assertEqual(library[S("library_other::A")][0], 2.0)

    def test_mapping_discovery_and_snapshots(self):
        rep = Representation.euc(2)
        signature = TensorName("library_tests::snapshot")(rep)
        library = TensorLibrary()
        self.assertEqual(len(library), 0)
        self.assertEqual(library.keys(), [])
        self.assertEqual(library.values(), [])
        self.assertEqual(library.items(), [])
        library.register(Tensor.dense(signature, [1.0, 2.0]))
        keys = library.keys()
        iterator = iter(library)
        self.assertEqual(len(library), 1)
        self.assertEqual(keys[0].structure, signature.structure)
        key, stored = library.items()[0]
        self.assertEqual(library[key][:], stored[:])
        stored[0] = 9.0
        self.assertEqual(library[key][0], 1.0)
        value = library.values()[0]
        value[1] = 8.0
        self.assertEqual(library[key][1], 2.0)
        library.register(stored)
        self.assertEqual(library[key][0], 9.0)
        self.assertEqual(len(library), 1)
        other = TensorName("library_tests::another")(rep)
        library.register(Tensor.dense(other, [3.0, 4.0]))
        self.assertEqual(len(library), 2)
        self.assertEqual(len(keys), 1)
        self.assertEqual(len(list(iterator)), 1)
        self.assertEqual(
            [key.structure for key in library],
            [key.structure for key, _ in library.items()],
        )
        self.assertEqual(
            [value[:] for value in library.values()],
            [value[:] for _, value in library.items()],
        )

    def test_requested_axis_order_and_network_execution(self):
        m, e = Representation.mink(2), Representation.euc(3)
        name = TensorName("library_tests::mixed")
        library = TensorLibrary()
        library.register(Tensor.dense(name(m, e), [float(i) for i in range(6)]))
        self.assertEqual(library[name][:], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        transposed = library[name(e, m)]
        self.assertEqual(transposed.shape, (3, 2))
        self.assertEqual(transposed[:], [0.0, 3.0, 1.0, 4.0, 2.0, 5.0])
        network = transposed("i", "mu")
        network.execute()
        self.assertEqual(network.result_tensor()[:], transposed[:])

    def test_sparse_symbolic_complex_and_scalar_data(self):
        rep = Representation.euc(2)
        library = TensorLibrary()
        signature = TensorName("library_tests::sparse")(rep)
        sparse = Tensor.sparse(signature, Expression)
        x = S("library_tests::component")
        sparse[1] = x
        library.register(sparse)
        self.assertEqual(library[signature][:], [E("0"), x])
        complex_key = TensorName("library_tests::complex")(rep)
        library.register(Tensor.dense(complex_key, [1j, 2 + 3j]))
        self.assertEqual(library[complex_key][:], [1j, 2 + 3j])
        scalar = TensorName("library_tests::scalar")()
        library.register(Tensor.dense(scalar, [42.0]))
        self.assertEqual(library[scalar].structure.rank, 0)
        self.assertEqual(library[scalar][0], 42.0)

    def test_hep_entries_and_generic_factories(self):
        for library in [TensorLibrary.hep_lib(), TensorLibrary.hep_lib_atom()]:
            keys = library.keys()
            self.assertGreater(len(keys), 0)
            self.assertEqual(len(keys), len(library))
            for key, value in library.items():
                self.assertIsInstance(key, TensorExpression)
                self.assertEqual(value.structure, key.structure)
                # Verify that the enumerated signature can retrieve this value.
                self.assertEqual(library[key][:], value[:])  # noqa: PLR1733
            gamma = library[TensorExpression.dirac_gamma(4)]
            self.assertEqual(gamma.structure.shape, (4, 4, 4))
            self.assertEqual(gamma[:], library["spenso::gamma"][:])

        library = TensorLibrary()
        metric = library[TensorExpression.g(Representation.mink(2))]
        self.assertEqual(metric[:], [1.0, 0.0, 0.0, -1.0])
        self.assertEqual(len(library), 0)
        self.assertEqual(library.items(), [])
        with self.assertRaises(ValueError):
            library[TensorExpression.g(Representation.mink(S("library_tests::D")))]

    def test_explicit_symbolic_and_integer_component_data_remain_exact(self):
        rep = Representation.mink(4)
        momentum = TensorName.vector("library_tests::exact_momentum")(rep)
        values = [E("1/3"), E("2/7"), E("3/11"), E("5/13")]
        for tensor in (
            Tensor.dense(momentum, values),
            Tensor.sparse(momentum, Expression),
        ):
            for index, value in enumerate(values):
                tensor[index] = value
            self.assertIs(tensor.dtype, Expression)
            self.assertEqual(tensor[:], values)
            library = TensorLibrary.hep_lib_atom()
            library.register(tensor)
            self.assertEqual(library[momentum][:], values)
            square = momentum("mu") * momentum("mu")
            result = square.to_tensor(library)
            self.assertIs(result.dtype, Expression)
            self.assertEqual(result.scalar(), E("1/9-4/49-9/121-25/169"))

        large = 2**80 + 1
        vector = TensorName("library_tests::integer_vector")(Representation.euc(2))
        for values in ([large, -2], [E(str(large)), -2]):
            tensor = Tensor.dense(vector, values)
            self.assertIs(tensor.dtype, Expression)
            self.assertEqual(tensor[:], [E(str(large)), E("-2")])
            tensor[1] = 3
            self.assertEqual(tensor[1], E("3"))
            norm = (tensor("i") * tensor("i")).to_tensor()
            self.assertEqual(norm.scalar(), E(str(large**2 + 9)))
        sparse = Tensor.sparse(vector, Expression)
        sparse[0] = large
        self.assertEqual(sparse[:], [E(str(large)), E("0")])
        self.assertIs(Tensor.dense(vector, [1.0, 2.0]).dtype, float)
        self.assertIs(Tensor.dense(vector, [1j, 2j]).dtype, complex)

    def test_default_and_explicit_exact_hep_components_and_contractions(self):
        gamma = TensorExpression.dirac_gamma(4)
        generator = TensorExpression.color_t(8, 3)
        structure_constant = TensorExpression.color_f(8)
        metric = TensorExpression.g(Representation.mink(2))
        for library in (None, TensorLibrary.hep_lib_atom()):
            with self.subTest(library=library):
                components = [
                    descriptor.to_tensor(library)
                    for descriptor in (gamma, generator, structure_constant, metric)
                ]
                for tensor in components:
                    self.assertIs(tensor.dtype, Expression)
                    self.assertTrue(
                        all(isinstance(value, Expression) for value in tensor)
                    )
                gamma_data, generator_data, f_data, metric_data = components
                self.assertEqual(gamma_data[0, 2, 0], E("1"))
                self.assertEqual(gamma_data[0, 3, 2], E("-𝑖"))
                self.assertEqual(generator_data[0, 0, 1], E("1/2"))
                self.assertEqual(generator_data[7, 0, 0], E("sqrt(3)/6"))
                self.assertEqual(f_data[3, 4, 7], E("sqrt(3)/2"))
                self.assertEqual(metric_data[:], [E("1"), E("0"), E("0"), E("-1")])
                for source, expected in (
                    (
                        generator("a", "i", "j") * generator("a", "j", "i"),
                        E("4"),
                    ),
                    (structure_constant("a", "b", "c") ** 2, E("24")),
                    (gamma("a", "b", "mu") * gamma("b", "a", "mu"), E("16")),
                ):
                    result = source.to_tensor(library)
                    self.assertIs(result.dtype, Expression)
                    self.assertEqual(result.scalar(), expected)
                for value in ("0", "1"):
                    result = TensorExpression(E(value)).to_tensor(library)
                    self.assertIs(result.dtype, Expression)
                    self.assertEqual(result.scalar(), E(value))

        numerical = TensorLibrary.hep_lib()
        self.assertIs(generator.to_tensor(numerical).dtype, complex)
        self.assertIs(structure_constant.to_tensor(numerical).dtype, float)
        self.assertEqual(generator.to_tensor(numerical)[0, 0, 1], 0.5 + 0j)
        self.assertAlmostEqual(
            structure_constant.to_tensor(numerical)[3, 4, 7], 3**0.5 / 2
        )
        vector = TensorName("library_tests::numerical_vector")(Representation.euc(2))
        data = Tensor.dense(vector, [1.0, 2.0])
        numerical.register(data)
        self.assertIs(data.dtype, float)
        self.assertIs(numerical[vector].dtype, float)
        self.assertIs(vector("i").to_tensor(numerical).dtype, float)
        norm = (vector("i") * vector("i")).to_tensor(numerical)
        # Shared network scalars use Symbolica Atom; result wrapping preserves it.
        self.assertIs(norm.dtype, Expression)
        self.assertEqual(norm.scalar(), E("5.0"))
        self.assertEqual(
            norm.scalar().evaluate({}, decimal_digit_precision=100).real.precision, 53
        )

    def test_invalid_and_missing_keys(self):
        library = TensorLibrary()
        signature = TensorName("library_tests::missing")(Representation.euc(2))
        with self.assertRaises(KeyError):
            library[signature]
        with self.assertRaises(KeyError):
            library["library_tests::missing"]
        with self.assertRaises(ValueError):
            library[signature("i")]
        with self.assertRaises(TypeError):
            library[None]

    def test_membership_and_get_share_signature_resolution(self):
        rep = Representation.euc(2)
        name = TensorName("library_tests::mapping")
        signature = name(7, rep)
        library = TensorLibrary()
        library.register(Tensor.dense(signature, [1.0, 2.0]))
        for key in [
            signature,
            name(7, rep),
            library.keys()[0],
            name,
            S("library_tests::mapping"),
            "library_tests::mapping",
        ]:
            with self.subTest(key=str(key)):
                self.assertIn(key, library)
                stored = library.get(key)
                self.assertEqual(stored[:], [1.0, 2.0])
                stored[0] = 9.0
                self.assertEqual(library[key][0], 1.0)

        sentinel = object()
        for key in [
            name(8, rep),
            name(7, Representation.euc(3)),
            "library_tests::absent",
        ]:
            with self.subTest(missing=str(key)):
                self.assertNotIn(key, library)
                self.assertIsNone(library.get(key))
                self.assertIs(library.get(key, sentinel), sentinel)
                self.assertIs(library.get(key, default=sentinel), sentinel)
                with self.assertRaises(KeyError):
                    library[key]

        library.register(Tensor.dense(name(8, rep), [3.0, 4.0]))
        self.assertIn(signature, library)
        for key in [name, S("library_tests::mapping"), "library_tests::mapping"]:
            with self.subTest(ambiguous=str(key)):
                with self.assertRaisesRegex(KeyError, "ambiguous"):
                    self.assertIn(key, library)
                with self.assertRaisesRegex(KeyError, "ambiguous"):
                    library.get(key, sentinel)

    def test_factory_membership_and_get_preserve_errors_and_enumeration(self):
        library = TensorLibrary()
        metric = TensorExpression.g(Representation.mink(2))
        self.assertIn(metric, library)
        self.assertEqual(library.get(metric)[:], [1.0, 0.0, 0.0, -1.0])
        self.assertEqual(len(library), 0)
        self.assertEqual(library.keys(), [])
        # Name-only lookup cannot choose concrete dimensions for a factory.
        self.assertNotIn("spenso::g", library)
        self.assertIsNone(library.get("spenso::g"))
        for key, error in [
            (None, TypeError),
            (metric("mu", "nu"), ValueError),
            (
                TensorExpression.g(Representation.mink(S("library_tests::D"))),
                ValueError,
            ),
        ]:
            with self.subTest(invalid=str(key)):
                with self.assertRaises(error):
                    self.assertIn(key, library)
                with self.assertRaises(error):
                    library.get(key, object())


if __name__ == "__main__":
    unittest.main()
