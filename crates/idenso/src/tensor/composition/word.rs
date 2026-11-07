use crate::tensor::SymbolicTensor;
use crate::tensor::composition;
use crate::tensor::composition::ExplicitIndexOccurrences;
use crate::tensor::inference::{InterfaceInference, TensorInferenceError, TensorInterface};
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};

impl SymbolicTensor<PartialStructure> {
    /// Build a variadic tensor word without discarding its operand interfaces.
    ///
    /// The head must carry [`SPENSO_TAG`]'s `tensor_word` tag. Its existing
    /// Symbolica attributes own normalization; for example `Flat`, `Linear`,
    /// and `Antisymmetric` describe an exterior word of formally odd operands.
    /// Field statistics and mixed-parity admission belong to the caller.
    /// Commuting coefficients must satisfy Symbolica's native `Scalar` property.
    /// Keep tensor-valued coefficients outside the word; they are not new odd
    /// operands and are not extracted by native multilinear normalization.
    ///
    /// Every unresolved operand occurrence receives independent explicit ports
    /// before normalization. Thus two copies of an unresolved vector do not
    /// become identical components and vanish. Written explicit indices retain
    /// their identities and contract through ordinary tensor interface rules.
    /// Native normalization may produce zero; the resulting tensor retains its
    /// planned interface. A detached zero Atom cannot recover that information.
    ///
    /// The empty word is the scalar identity. The constructor does not assign
    /// an adjoint convention: explicitly adjoint the operands and choose their
    /// order before constructing the result.
    /// Symbolic networks retain the word; component-shadow realization requires
    /// an explicit policy and is unsupported by the named component parser.
    ///
    /// # Errors
    /// Rejects invalid word heads, unpublished operands, unsupported operand
    /// grammar, excessive explicit index occurrences, and normalization that
    /// changes the planned tensor interface.
    ///
    /// # Example
    ///
    /// ```
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::inference::TensorInferenceError;
    /// use spenso::network::tags::SPENSO_TAG;
    /// use spenso::{mink, tensor_symbol};
    /// use symbolica::{function, symbol};
    /// let head = symbol!("word_example::exterior"; Flat, Linear, Antisymmetric;
    ///     tag = &SPENSO_TAG.tensor_word);
    /// let a = SymbolicTensor::infer(function!(tensor_symbol!("word_example::a"), mink!(4, 1)))?;
    /// let zero = SymbolicTensor::word(head, &[a.clone(), a])?;
    /// assert!(zero.is_zero());
    /// # Ok::<(), TensorInferenceError>(())
    /// ```
    pub fn word(head: Symbol, operands: &[Self]) -> Result<Self, TensorInferenceError> {
        if !head.has_tag(&SPENSO_TAG.tensor_word)
            || head.is_scalar()
            || head.has_tag(&SPENSO_TAG.tensor)
            || head.has_tag(&SPENSO_TAG.rank1)
            || head.has_tag(&SPENSO_TAG.representation)
            || head.has_tag(&SPENSO_TAG.broadcast)
        {
            return Err(TensorInferenceError::invalid(
                "a tensor word head must carry only the tensor-word structural role",
            ));
        }
        if operands.is_empty() {
            return Self::infer(Atom::one());
        }
        for operand in operands {
            operand.explicit_index_occurrences()?;
            if !Self::is_word_operand(head, operand.expression.as_view()) {
                return Err(TensorInferenceError::invalid(
                    "word operands must be individual generators or same-head words with scalar coefficients; keep commuting tensor coefficients outside the word",
                ));
            }
        }
        let dummies = Self::reserved_dummies(operands);
        let materialized = operands
            .iter()
            .map(|operand| {
                let replacements = operand
                    .structure
                    .open_positions()
                    .into_iter()
                    .map(|position| (position, dummies.fresh_index()))
                    .collect();
                operand.reindex_interface_ports(&replacements)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let occurrences =
            ExplicitIndexOccurrences::product(&materialized.iter().collect::<Vec<_>>())?;
        let interface = InterfaceInference::merge_explicit_interface_sequence(
            &materialized
                .iter()
                .map(|operand| operand.structure.clone())
                .collect::<Vec<_>>(),
        )?
        .into_algebra_result();
        let expression = FunctionBuilder::new(head)
            .add_args(
                materialized
                    .iter()
                    .map(|operand| operand.expression.as_view()),
            )
            .finish();
        Self::validate_interface(&expression, &interface)?;
        Self::from_validated_parts(expression, interface, Arc::new(occurrences)).finish_parts()
    }

    /// Check the operand grammar without rewriting or assigning field statistics.
    pub(in crate::tensor) fn is_word_operand(head: Symbol, value: AtomView<'_>) -> bool {
        if value.is_zero() {
            return true;
        }
        if value.is_scalar().is_true() {
            return false;
        }
        match value {
            AtomView::Var(_) => true,
            AtomView::Fun(function) => {
                function.get_symbol() == head
                    || (!composition::is_composite_head(function.get_symbol(), &SPENSO_TAG)
                        && !function.get_symbol().has_tag(&SPENSO_TAG.tensor_word))
            }
            AtomView::Add(sum) => sum.iter().all(|term| Self::is_word_operand(head, term)),
            AtomView::Mul(product) => {
                let mut operands = product
                    .iter()
                    .filter(|factor| !factor.is_scalar().is_true());
                matches!((operands.next(), operands.next()),
                    (Some(operand), None) if Self::is_word_operand(head, operand))
            }
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AdjointError;
    use crate::tensor::{SymbolicNetExt, SymbolicNetParse};
    use crate::test_support;
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::dimension::Dimension;
    use spenso::structure::partial::PartialIndex;
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use spenso::structure::slot::IsAbstractSlot;
    use std::collections::HashMap;
    use std::slice;
    use symbolica::{function, symbol};

    fn head() -> Symbol {
        test_support::test_initialize();
        symbol!("tensor_word_tests::exterior"; Flat, Linear, Antisymmetric;
            tag = &SPENSO_TAG.tensor_word)
    }

    fn component(name: &str, index: usize) -> SymbolicTensor<PartialStructure> {
        let port = Minkowski {}
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index));
        SymbolicTensor::infer(function!(SPENSO_TAG.tensor_symbol(name), port.to_atom())).unwrap()
    }

    fn unresolved(name: &str) -> SymbolicTensor<PartialStructure> {
        let rep = Minkowski {}.new_rep(4).cast::<LibraryRep>();
        SymbolicTensor::from_declared_parts(
            function!(SPENSO_TAG.tensor_symbol(name), rep.to_symbolic([])),
            PartialStructure::from_logical_slots([rep.slot(PartialIndex::open(0))]),
        )
        .unwrap()
    }

    #[test]
    fn tensor_word_reuses_native_exterior_normalization() {
        let head = head();
        let a = component("tensor_word_tests::a", 1);
        let b = component("tensor_word_tests::b", 2);
        let c = component("tensor_word_tests::c", 3);
        let ab = SymbolicTensor::word(head, &[a.clone(), b.clone()]).unwrap();
        let ba = SymbolicTensor::word(head, &[b.clone(), a.clone()]).unwrap();
        assert_eq!(ab.rank(), 2);
        assert_eq!(ab.expression, -ba.expression);
        assert!(
            SymbolicTensor::word(head, &[a.clone(), a.clone()])
                .unwrap()
                .is_zero()
        );
        let nested = SymbolicTensor::word(head, &[ab.clone(), c.clone()]).unwrap();
        let direct = SymbolicTensor::word(head, &[a, b, c.clone()]).unwrap();
        assert_eq!(nested.expression, direct.expression);
        assert_eq!(
            SymbolicTensor::word(head, &[c, ab]).unwrap().expression,
            direct.expression
        );
    }

    #[test]
    fn tensor_word_distinguishes_unresolved_occurrences_before_normalization() {
        let head = head();
        let a = unresolved("tensor_word_tests::unresolved");
        let word = SymbolicTensor::word(head, &[a.clone(), a]).unwrap();
        assert_eq!(word.rank(), 2);
        assert!(!word.is_zero());
        let ports = word.structure.logical_slots();
        assert_ne!(ports[0].aind, ports[1].aind);
        assert!(word.structure.open_positions().is_empty());
        assert_eq!(
            SymbolicTensor::infer(word.expression.clone())
                .unwrap()
                .rank(),
            2
        );
    }

    #[test]
    fn tensor_word_keeps_ranked_zero_compatible_with_a_nonzero_word() {
        let head = head();
        let a = component("tensor_word_tests::zero_a", 31);
        let b = component("tensor_word_tests::zero_b", 32);
        let zero_a = a.with_algebra_result(Atom::Zero).unwrap();
        let zero = SymbolicTensor::word(head, &[zero_a, b.clone()]).unwrap();
        let nonzero = SymbolicTensor::word(head, &[a, b]).unwrap();
        assert!(zero.is_zero());
        assert_eq!(zero.rank(), 2);
        assert!(zero.try_add(&nonzero).unwrap().same_value(&nonzero));
        assert!(nonzero.try_add(&zero).unwrap().same_value(&nonzero));
        assert_eq!(SymbolicTensor::infer(Atom::Zero).unwrap().rank(), 0);
    }

    #[test]
    fn tensor_word_preserves_native_real_and_complex_scalar_prefactors() {
        let head = head();
        let a = component("tensor_word_tests::weighted_a", 41);
        let b = component("tensor_word_tests::weighted_b", 42);
        let real = Atom::var(symbol!("tensor_word_tests::real"; Real, Scalar));
        let complex = Atom::var(symbol!("tensor_word_tests::complex"; Scalar)) + Atom::i();
        let base = SymbolicTensor::word(head, &[a.clone(), b.clone()]).unwrap();
        for coefficient in [real, complex, Atom::num((2, 3))] {
            let weighted = a
                .with_algebra_result(&coefficient * a.expression())
                .unwrap();
            let actual = SymbolicTensor::word(head, &[weighted, b.clone()]).unwrap();
            assert_eq!(actual.expression, coefficient * base.expression());
            assert_eq!(actual.rank(), 2);
        }
    }

    #[test]
    fn tensor_word_exposes_ports_to_reindexing_dimension_maps_and_literal_renaming() {
        let head = head();
        let a = component("tensor_word_tests::mapped_a", 51);
        let b = component("tensor_word_tests::mapped_b", 52);
        let word = SymbolicTensor::word(head, &[a, b]).unwrap();
        let renamed = word
            .reindex_interface_ports(&HashMap::from([
                (0, AbstractIndex::Normal(61)),
                (1, AbstractIndex::Normal(62)),
            ]))
            .unwrap();
        assert_eq!(renamed.rank(), 2);
        assert_ne!(word.expression, renamed.expression);
        let dimension = Atom::var(symbol!("tensor_word_tests::D"));
        let target_dimension = symbol!("tensor_word_tests::E");
        let mapped = renamed
            .with_lorentz_dimension(Dimension::try_from(dimension.as_view()).unwrap())
            .unwrap();
        assert_eq!(mapped.rank(), 2);
        assert!(mapped.expression.contains(&dimension));
        let literal = mapped
            .with_renamed_symbols(&HashMap::from([(
                dimension.get_symbol().unwrap(),
                target_dimension,
            )]))
            .unwrap();
        assert!(!literal.expression.contains(&dimension));
        assert!(literal.expression.contains_symbol(target_dimension));
        assert!(literal.expression.contains_symbol(head));
        let scoped = literal
            .wrap_indices(symbol!("tensor_word_tests::owned_scope"), false)
            .unwrap();
        assert_eq!(scoped.rank(), 2);
        assert!(
            scoped
                .structure
                .logical_slots()
                .iter()
                .all(|slot| matches!(slot.aind, PartialIndex::Explicit(AbstractIndex::Scoped(_))))
        );
        assert!(scoped.expression.contains_symbol(head));
    }

    #[test]
    fn tensor_word_rejects_unresolved_raw_ingress_and_implicit_adjoint() {
        let head = head();
        let a = unresolved("tensor_word_tests::raw_a");
        let b = unresolved("tensor_word_tests::raw_b");
        let raw = function!(head, a.expression(), b.expression());
        assert!(SymbolicTensor::infer(raw).is_err());
        let word = SymbolicTensor::word(head, &[a, b]).unwrap();
        assert!(matches!(
            word.dirac_adjoint(true),
            Err(AdjointError::ImplicitTensorWordAdjoint)
        ));
    }

    #[test]
    fn tensor_word_rejects_invalid_scope_and_preserves_zero_interface_mapping() {
        let head = head();
        let illegal = SymbolicTensor::new(
            function!(
                SPENSO_TAG.tensor_symbol("tensor_word_tests::illegal"),
                SPENSO_TAG.chain_in
            ),
            PartialStructure::from_logical_slots([]),
        );
        assert!(illegal.is_err());
        let mut invalidated = component("tensor_word_tests::invalidated", 70);
        invalidated.invalidate_proofs();
        assert!(SymbolicTensor::word(head, &[invalidated]).is_err());
        let a = component("tensor_word_tests::mapped_zero", 71);
        let zero = a.with_algebra_result(Atom::Zero).unwrap();
        let word = SymbolicTensor::word(head, &[zero]).unwrap();
        let dimension = Atom::var(symbol!("tensor_word_tests::zero_D"));
        let mapped = word
            .with_lorentz_dimension(Dimension::try_from(dimension.as_view()).unwrap())
            .unwrap();
        assert!(mapped.is_zero());
        assert_eq!(mapped.rank(), 1);
        assert_eq!(
            mapped.structure.logical_slots()[0].rep().dim.to_symbolic(),
            dimension
        );
    }

    #[test]
    fn tensor_word_stays_intact_when_contracting_a_commuting_tensor_coefficient() {
        let head = head();
        let a = component("tensor_word_tests::contract_a", 81);
        let b = component("tensor_word_tests::contract_b", 82);
        let word = SymbolicTensor::word(head, &[a, b]).unwrap();
        let metric = SymbolicTensor::infer(function!(
            ETS.metric,
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(81))
                .to_atom(),
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(82))
                .to_atom()
        ))
        .unwrap();
        let product = word.multiply(&metric).unwrap();
        assert!(product.is_scalar());
        let contracted = product.contract(Default::default()).unwrap();
        assert!(contracted.is_scalar());
        assert!(contracted.expression.contains_symbol(head));
        assert!(!contracted.expression.contains_symbol(ETS.metric));
    }

    #[test]
    fn tensor_word_survives_symbolic_network_execution() {
        let word = SymbolicTensor::word(
            head(),
            &[
                component("tensor_word_tests::network_a", 85),
                component("tensor_word_tests::network_b", 86),
            ],
        )
        .unwrap();
        let network = word
            .expression
            .parse_to_symbolic_net::<AbstractIndex>(&Default::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), 2);
        assert_eq!(network.simple_execute::<()>().unwrap(), word.expression);
    }

    #[test]
    fn tensor_word_extension_keeps_commuting_tensor_coefficients_outside() {
        let head = head();
        let a = component("tensor_word_tests::extended_a", 91);
        let b = component("tensor_word_tests::extended_b", 92);
        let c = component("tensor_word_tests::extended_c", 93);
        let coefficient = component("tensor_word_tests::coefficient", 91);
        let ab = SymbolicTensor::word(head, &[a.clone(), b.clone()]).unwrap();
        let weighted = coefficient.multiply(&ab).unwrap();
        assert!(SymbolicTensor::word(head, &[weighted, c.clone()]).is_err());
        let extended = SymbolicTensor::word(head, &[ab, c.clone()]).unwrap();
        let actual = coefficient.multiply(&extended).unwrap();
        let direct = SymbolicTensor::word(head, &[a, b, c]).unwrap();
        let expected = coefficient.multiply(&direct).unwrap();
        assert!(actual.same_value(&expected));
        assert_eq!(actual.rank(), 2);
        assert!(actual.expression.contains_symbol(head));
    }

    #[test]
    fn tensor_word_preserves_multilinearity_and_field_metadata() {
        let head = head();
        let a = component("tensor_word_tests::linear_a", 101);
        let b = component("tensor_word_tests::linear_b", 101);
        let jet = SymbolicTensor::infer(function!(
            SPENSO_TAG.tensor_symbol("tensor_word_tests::jet"),
            2,
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(102))
                .to_atom()
        ))
        .unwrap();
        let sum = a.try_add(&b).unwrap();
        let actual = SymbolicTensor::word(head, &[sum, jet.clone()]).unwrap();
        let expected = SymbolicTensor::word(head, &[a, jet.clone()])
            .unwrap()
            .try_add(&SymbolicTensor::word(head, &[b, jet.clone()]).unwrap())
            .unwrap();
        assert!(actual.same_value(&expected));
        assert!(actual.expression.contains(jet.expression()));
    }

    #[test]
    fn tensor_word_checks_native_callback_results_and_head_roles() {
        head();
        let a = component("tensor_word_tests::callback_a", 111);
        let erased = symbol!(
            "tensor_word_tests::erased",
            tag = &SPENSO_TAG.tensor_word,
            norm = |_, output| {
                **output = Atom::one();
            }
        );
        assert!(SymbolicTensor::word(erased, slice::from_ref(&a)).is_err());
        let annihilated = symbol!(
            "tensor_word_tests::annihilated",
            tag = &SPENSO_TAG.tensor_word,
            norm = |_, output| {
                **output = Atom::Zero;
            }
        );
        let zero = SymbolicTensor::word(annihilated, slice::from_ref(&a)).unwrap();
        assert!(zero.is_zero());
        assert_eq!(zero.rank(), 1);
        let conflicting = symbol!("tensor_word_tests::conflicting"; Scalar;
            tag = &SPENSO_TAG.tensor_word);
        assert!(SymbolicTensor::word(conflicting, &[a]).is_err());
        assert!(SymbolicTensor::infer(Atom::var(head())).is_err());
        assert!(
            SymbolicTensor::word(head(), &[])
                .unwrap()
                .expression
                .is_one()
        );
    }

    #[test]
    fn tensor_word_retains_hidden_index_occurrences_after_native_zero() {
        let head = head();
        let a = component("tensor_word_tests::history_a", 121);
        let zero = SymbolicTensor::word(head, &[a.clone(), a.clone()]).unwrap();
        assert!(zero.is_zero());
        assert!(zero.multiply(&a).is_err());
        assert!(SymbolicTensor::word(head, &[zero, a]).is_err());
    }
}
