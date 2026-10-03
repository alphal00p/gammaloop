"""Construction, rewrite patterns, and factor projectors in an installed host."""

import unittest

from symbolica import E, Expression, S
from symbolica.community.tensor import (
    AUTO,
    FactorProjector,
    Nc,
    PortPattern,
    Representation,
    Tensor,
    TensorExpression,
    TensorLibrary,
    TensorName,
    TensorNetwork,
    TensorPattern,
    chain,
    dot,
    trace,
)


class TensorVocabularyTests(unittest.TestCase):
    def test_named_factories_assign_distinct_open_ports(self):
        spin, mink = Representation.bis(4), Representation.mink(4)
        library = TensorLibrary.hep_lib_atom()
        for name, factory in (
            ("gamma0", TensorExpression.gamma0),
            ("charge_conjugation", TensorExpression.charge_conjugation),
        ):
            with self.subTest(name=name):
                tensor = factory(4)
                self.assertEqual(tensor.rank, 2)
                self.assertNotEqual(tensor.to_expression(), E("0"))
                indexed = tensor("a", "b")
                reference = library[f"spenso::{name}"].expression()("a", "b")
                self.assertEqual(indexed, reference)
                self.assertEqual(
                    indexed.to_tensor(library)[:], reference.to_tensor(library)[:]
                )

        epsilon = TensorExpression.levi_civita(mink)
        self.assertEqual(epsilon.rank, 4)
        indexed = epsilon("mu", "nu", "rho", "sigma")
        self.assertNotEqual(indexed.to_expression(), E("0"))
        self.assertEqual(
            epsilon("nu", "mu", "rho", "sigma").to_expression(),
            -indexed.to_expression(),
        )
        self.assertEqual(epsilon("mu", "mu", "rho", "sigma").to_expression(), E("0"))
        self.assertEqual(
            TensorExpression.levi_civita(
                Representation.euc(S("vocabulary::D")), rank=3
            ).rank,
            3,
        )
        # Repeated arguments are legitimately zero; factories express distinct ports.
        self.assertEqual(
            TensorName.charge_conjugation()(spin, spin).to_expression(), E("0")
        )
        self.assertEqual(
            TensorName.levi_civita()(mink, mink, mink, mink).to_expression(), E("0")
        )

    def test_scalar_constant_and_invariant_patterns(self):
        nc = Nc()
        self.assertIsInstance(nc, Expression)
        self.assertEqual(nc.get_name(), "spenso::Nc")
        self.assertEqual(nc.evaluate({}), 3)
        degree, rep = S("vocabulary::degree_", "vocabulary::rep_")
        for value, pattern in (
            (Representation.cof(3).casimir(), TensorPattern.casimir(degree, rep)),
            (
                Representation.coad(8).dynkin_index(),
                TensorPattern.dynkin_index(degree, rep),
            ),
        ):
            self.assertEqual(value.replace(pattern, 7), E("7"))

    def test_compact_patterns_and_reversed_contextual_ports(self):
        spin, mink = Representation.bis(4), Representation.mink(4)
        dimension, i, j, mu, fs, left, right = S(
            "vocabulary::dimension_",
            "vocabulary::i_",
            "vocabulary::j_",
            "vocabulary::mu_",
            "vocabulary::factors___",
            "vocabulary::left_",
            "vocabulary::right_",
        )
        gamma = TensorExpression.dirac_gamma(4)
        word = chain(
            spin("a"), spin("b"), gamma(AUTO, AUTO, "mu"), gamma(AUTO, AUTO, "nu")
        )
        closed = trace(spin, gamma(AUTO, AUTO, "mu"), gamma(AUTO, AUTO, "nu"))
        p = TensorName.vector("vocabulary::p")(mink)
        for value, pattern in (
            (dot(p, p), TensorPattern.dot(left, right)),
            (word, TensorPattern.chain(left, right, fs)),
            (
                closed,
                TensorPattern.trace(
                    PortPattern.exact(Representation.bis(dimension)), fs
                ),
            ),
            (
                TensorExpression.gamma0(4)("a", "b"),
                TensorPattern.gamma0(dimension, i, j),
            ),
            (
                TensorExpression.charge_conjugation(4)("a", "b"),
                TensorPattern.charge_conjugation(dimension, i, j),
            ),
            (
                TensorExpression.levi_civita(mink, rank=2)("a", "b"),
                TensorPattern.levi_civita(Representation.mink(dimension), i, j),
            ),
        ):
            with self.subTest(pattern=str(pattern)):
                self.assertEqual(value.to_expression().replace(pattern, 7), E("7"))
        incoming, outgoing = PortPattern.chain_in(), PortPattern.chain_out()
        forward = TensorPattern.dirac_gamma(dimension, incoming, outgoing, mu)
        reverse = TensorPattern.dirac_gamma(dimension, outgoing, incoming, mu)
        factor = gamma(AUTO, AUTO, "mu")
        reversed_factor = factor.permute_axes((1, 0, 2))
        reversed_word = chain(*reversed_factor.axes[:2], reversed_factor)
        self.assertEqual(
            reversed_word.to_expression().replace(
                TensorPattern.chain(left, right, reverse), 7
            ),
            E("7"),
        )
        self.assertNotEqual(
            reversed_word.to_expression().replace(
                TensorPattern.chain(left, right, forward), 7
            ),
            E("7"),
        )

    def test_symbolic_and_component_projectors_agree(self):
        rep = Representation.euc(2)
        a, b = (TensorName(f"vocabulary::{name}")(rep, rep) for name in ("A", "B"))
        data_a = Tensor.dense(a, [1.0, 2.0, 3.0, 4.0])
        data_b = Tensor.dense(b, [0.0, 1.0, -1.0, 0.0])
        library = TensorLibrary()
        library.register(data_a)
        library.register(data_b)
        for factory, expected in (
            (FactorProjector.symmetric, [0.5, 2.5, -2.5, 0.5]),
            (FactorProjector.antisymmetric, [-2.5, -1.5, -1.5, 2.5]),
            (FactorProjector.cyclic, [0.5, 2.5, -2.5, 0.5]),
        ):
            with self.subTest(factory=factory.__name__):
                group = factory(a, b)
                word = chain(rep("i"), rep("j"), group)
                self.assertIsInstance(word, TensorExpression)
                self.assertEqual(word.to_tensor(library)[:], expected)
                expanded = word.expand_projectors()
                self.assertEqual(expanded.structure.axes, word.structure.axes)
                self.assertEqual(expanded.to_tensor(library)[:], expected)
                self.assertEqual(expanded.expand_projectors(), expanded)
                network = chain(rep("i"), rep("j"), factory(data_a, data_b))
                self.assertIsInstance(network, TensorNetwork)
                network.execute()
                self.assertEqual(network.result_tensor()[:], expected)
        with self.assertRaises(ValueError):
            FactorProjector.symmetric()

    def test_projector_groups_remain_local_to_the_selected_factors(self):
        spin = Representation.bis(4)
        gamma = TensorExpression.dirac_gamma(4)
        factors = [gamma(AUTO, AUTO, name) for name in ("mu", "nu", "rho", "sigma")]
        word = chain(
            spin("a"),
            spin("b"),
            factors[0],
            FactorProjector.antisymmetric(*factors[1:3]),
            factors[3],
        )
        expected = (
            chain(spin("a"), spin("b"), *factors)
            - chain(
                spin("a"), spin("b"), factors[0], factors[2], factors[1], factors[3]
            )
        ) / 2
        self.assertEqual(
            (
                word.expand_projectors().to_expression() - expected.to_expression()
            ).expand(),
            E("0"),
        )
        projected_trace = trace(spin, FactorProjector.symmetric(*factors[:2]))
        fs, rep = S("vocabulary::fs___", "vocabulary::rep_")
        pattern = TensorPattern.trace(rep, TensorPattern.symmetric(fs))
        self.assertEqual(projected_trace.to_expression().replace(pattern, 7), E("7"))

    def test_projector_signs_and_scalar_coefficients_survive_composition(self):
        rep = Representation.euc(2)
        a, b = (
            TensorName(f"vocabulary::{name}")(rep, rep)
            for name in ("Asigned", "Bsigned")
        )
        data_a = Tensor.dense(a, [1.0, 2.0, 3.0, 4.0])
        data_b = Tensor.dense(b, [0.0, 1.0, -1.0, 0.0])
        library = TensorLibrary()
        library.register(data_a)
        library.register(data_b)
        for left, right in ((b, a), (-a, b), (data_b, data_a)):
            group = FactorProjector.antisymmetric(left, right)
            result = chain(rep("i"), rep("j"), group)
            if isinstance(result, TensorNetwork):
                result.execute()
                components = result.result_tensor()[:]
            else:
                components = result.to_tensor(library)[:]
            self.assertEqual(components, [2.5, 1.5, 1.5, -2.5])
        weighted = 2 * a.compose(b, left=(0, 1), right=(0, 1))
        self.assertEqual(trace(rep, weighted).to_tensor(library)[:], [2.0])
        for repeated in (a, data_a):
            result = chain(
                rep("i"), rep("j"), FactorProjector.antisymmetric(repeated, repeated)
            )
            if isinstance(result, TensorNetwork):
                result.execute()
                result = result.result_tensor()
            else:
                result = result.to_tensor(library)
            self.assertEqual(result.rank, 2)
            self.assertEqual(result[:], [0.0] * 4)

    def test_component_projectors_keep_spectator_axes_attached_to_factors(self):
        matrix_rep, spectator = Representation.euc(2), Representation.mink(2)
        a, b = (
            TensorName(f"vocabulary::{name}")(matrix_rep, matrix_rep, spectator)
            for name in ("A3", "B3")
        )
        values_a = [float(i + 1) for i in range(8)]
        values_b = [float(2 * i - 3) for i in range(8)]
        data_a, data_b = Tensor.dense(a, values_a), Tensor.dense(b, values_b)
        expected = [
            sum(
                (
                    values_a[(i * 2 + k) * 2 + mu] * values_b[(k * 2 + j) * 2 + nu]
                    - values_b[(i * 2 + k) * 2 + nu] * values_a[(k * 2 + j) * 2 + mu]
                )
                / 2
                for k in range(2)
            )
            for i in range(2)
            for j in range(2)
            for mu in range(2)
            for nu in range(2)
        ]
        projected = chain(
            matrix_rep("i"),
            matrix_rep("j"),
            FactorProjector.antisymmetric(data_a, data_b),
        )
        projected.execute()
        self.assertEqual(projected.result_tensor()[:], expected)

    def test_nested_component_groups_and_invalid_factors(self):
        rep = Representation.euc(2)
        a, b = (TensorName(f"vocabulary::{name}")(rep, rep) for name in ("An", "Bn"))
        data_a = Tensor.dense(a, [1.0, 2.0, 3.0, 4.0])
        data_b = Tensor.dense(b, [0.0, 1.0, -1.0, 0.0])
        library = TensorLibrary()
        library.register(data_a)
        library.register(data_b)
        symbolic = trace(
            rep, FactorProjector.symmetric(a, FactorProjector.antisymmetric(a, b))
        )
        concrete = trace(
            rep,
            FactorProjector.symmetric(
                data_a, FactorProjector.antisymmetric(data_a, data_b)
            ),
        )
        concrete.execute()
        self.assertEqual(concrete.result_tensor()[:], symbolic.to_tensor(library)[:])
        with self.assertRaisesRegex(ValueError, "individual matrix factors"):
            FactorProjector.symmetric(a.compose(b, left=(0, 1), right=(0, 1)), a)
        with self.assertRaises(ValueError):
            FactorProjector.symmetric(a, TensorName.vector("vocabulary::v")(rep))


if __name__ == "__main__":
    unittest.main()
