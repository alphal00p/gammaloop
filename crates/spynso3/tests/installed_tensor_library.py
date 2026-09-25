"""Stored library data and signature lookup against the installed extension."""

import unittest

from symbolica import E, Expression, S
from symbolica.community.spenso import (
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
        self.assertEqual(transposed.structure.shape, (3, 2))
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
            gamma = library[TensorExpression.gamma(4)]
            self.assertEqual(gamma.structure.shape, (4, 4, 4))
            self.assertEqual(gamma[:], library["spenso::gamma"][:])

        library = TensorLibrary()
        metric = library[TensorExpression.g(Representation.mink(2))]
        self.assertEqual(metric[:], [1.0, 0.0, 0.0, -1.0])
        self.assertEqual(len(library), 0)
        self.assertEqual(library.items(), [])
        with self.assertRaises(ValueError):
            library[TensorExpression.g(Representation.mink(S("library_tests::D")))]

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
