use crate::AdjointError;
use crate::tensor::SymbolicTensor;
use crate::tensor::composition;
use crate::tensor::composition::{ExplicitIndexOccurrences, ProductPlan, TensorCompositionError};
use crate::tensor::inference::{InterfaceInference, TensorInferenceError, TensorInterface};
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::slot::{DualSlotTo, ParseableAind};
use std::mem;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::symbol;

#[cfg(test)]
mod adjoint_tests;

#[derive(Clone, Copy)]
pub(super) enum ProductAlgebra {
    Commuting,
    Word(Symbol),
}

impl ProductAlgebra {
    pub(super) fn expression(
        self,
        left: AtomView<'_>,
        right: AtomView<'_>,
    ) -> Result<Atom, TensorCompositionError> {
        match self {
            Self::Commuting => Ok(SymbolicTensor::bracket_product(left, right)),
            Self::Word(head) => SymbolicTensor::word_product_expression(head, left, right),
        }
    }

    pub(super) fn finish(
        self,
        expression: Atom,
        interface: PartialStructure,
        occurrences: Arc<ExplicitIndexOccurrences>,
        sources: &[&SymbolicTensor<PartialStructure>],
        operation: &str,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        let interface = match self {
            Self::Commuting => interface,
            Self::Word(_) => {
                let interface = composition::interface_in_atom_order(&expression, &interface)?;
                SymbolicTensor::validate_interface(&expression, &interface)?;
                interface
            }
        };
        SymbolicTensor::from_validated_parts(expression, interface, occurrences)
            .validate_rewrite(operation, sources)
    }
}

impl ProductPlan {
    /// Apply this positional product using a declared tensor-word algebra.
    ///
    /// Odd generators must already be marked by `SymbolicTensor::word`.
    /// Every summand may contain one same-head word, multiplied by any
    /// commuting, word-free tensor coefficient. The output combines the left
    /// word followed by the right word, leaving coefficients outside.
    /// Symbolica's existing `Flat`, `Linear`, and (when selected)
    /// `Antisymmetric` attributes own flattening, distribution, and signs.
    ///
    /// Alignment and contraction use this plan's original logical port
    /// positions. Surviving unresolved coefficient ports remain unresolved.
    /// The result's logical interface follows its normalized expression;
    /// collecting words can move explicit ports relative to coefficient ports.
    /// Existing explicit scopes are retained; independent bound indices must
    /// be scoped with the native `wrap_indices` API before composition.
    ///
    /// `Compose` is intentionally unsupported: use an explicit `Contract`
    /// plan for word endpoints and keep compact matrix chains word-free.
    /// This operation assigns no conjugation or adjoint convention.
    ///
    /// # Errors
    /// Rejects invalid heads, unpublished operands, unsupported word grammar,
    /// incompatible ports, and excessive index multiplicities.
    ///
    /// # Example
    /// ```
    /// use idenso::tensor::SymbolicTensor;
    /// use idenso::tensor::composition::{ProductPlan, TensorCompositionError};
    /// use spenso::network::tags::SPENSO_TAG;
    /// use spenso::{mink, tensor_symbol};
    /// use symbolica::{function, symbol};
    /// let head = symbol!("word_product_example::exterior"; Flat, Linear, Antisymmetric;
    ///     tag = &SPENSO_TAG.tensor_word);
    /// let a = SymbolicTensor::infer(function!(tensor_symbol!("word_product_example::a"), mink!(4, 1)))?;
    /// let b = SymbolicTensor::infer(function!(tensor_symbol!("word_product_example::b"), mink!(4, 2)))?;
    /// let a = SymbolicTensor::word(head, &[a])?;
    /// let b = SymbolicTensor::word(head, &[b])?;
    /// let product = ProductPlan::Outer.apply_word(head, &a, &b)?;
    /// assert!(product.same_value(&a.multiply_word(head, &b)?));
    /// assert!(a.pow_word(head, 2)?.is_zero());
    /// # Ok::<(), TensorCompositionError>(())
    /// ```
    pub fn apply_word(
        &self,
        head: Symbol,
        left: &SymbolicTensor<PartialStructure>,
        right: &SymbolicTensor<PartialStructure>,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorCompositionError> {
        SymbolicTensor::word(head, &[])?;
        if !head.is_flat() || !head.is_linear() {
            return Err(TensorCompositionError::InvalidResultInterface(
                "word multiplication requires a Flat, Linear tensor-word head".into(),
            ));
        }
        for operand in [left, right] {
            operand.explicit_index_occurrences()?;
            for term in operand.expression.expand_in(Atom::var(head)).terms() {
                SymbolicTensor::split_word_coefficient(
                    head,
                    term,
                    &mut Vec::new(),
                    &mut Vec::new(),
                    &mut false,
                )?;
            }
        }
        self.apply_algebra(left, right, ProductAlgebra::Word(head))
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Dirac adjoint in a tensor-word algebra with an explicit odd-generator map.
    ///
    /// Each term has the same grammar as [`Self::multiply_word`]: at most one
    /// selected word and a commuting, word-free tensor coefficient. The native
    /// Dirac adjoint conjugates that coefficient with fixed physical indices;
    /// this operation reverses the generator order and delegates each generator
    /// to `adjoint_generator`. Symbolica owns distribution and exterior signs.
    ///
    /// The callback supplies the field convention, for example a declared
    /// fermion/barred-fermion pairing. It must return an admissible, word-free
    /// operand with the source generator's dual representations, unchanged
    /// index labels, and unchanged logical axis order. Callback failures are
    /// propagated; no Majorana, ghost, or unknown-field policy is implicit.
    ///
    /// The source's logical interface is retained with dual representations.
    /// Unresolved ports are materialized before the operation, as with
    /// [`Self::dirac_adjoint`]. Ranked zero retains its interface and requires
    /// no generator callbacks. Ordinary `dirac_adjoint` continues to reject
    /// nonzero tensor words without an explicit policy.
    ///
    /// # Errors
    /// Propagates invalid word grammar, generator policy failures, incompatible
    /// dual ports, and invalid resulting tensor interfaces.
    ///
    /// # Example
    /// ```
    /// use idenso::AdjointError;
    /// use idenso::tensor::SymbolicTensor;
    /// use spenso::network::tags::SPENSO_TAG;
    /// use spenso::{mink, tensor_symbol};
    /// use symbolica::{function, symbol};
    /// let head = symbol!("word_adjoint_example::exterior"; Flat, Linear, Antisymmetric;
    ///     tag = &SPENSO_TAG.tensor_word);
    /// let a = SymbolicTensor::infer(function!(tensor_symbol!("word_adjoint_example::a"), mink!(4, 1)))?;
    /// let word = SymbolicTensor::word(head, &[a])?;
    /// // An explicit identity policy for a formally self-conjugate generator.
    /// let adjoint = word.dirac_adjoint_word(head, |generator| Ok(generator.clone()))?;
    /// assert!(adjoint.same_value(&word));
    /// # Ok::<(), AdjointError>(())
    /// ```
    pub fn dirac_adjoint_word(
        &self,
        head: Symbol,
        mut adjoint_generator: impl FnMut(&Self) -> Result<Self, AdjointError>,
    ) -> Result<Self, AdjointError> {
        let identity = Self::word(head, &[])?;
        // Reuse the product boundary's head and complete operand-grammar checks.
        ProductPlan::Scalar.apply_word(head, &identity, self)?;
        let source = self.materialized()?;
        let interface = PartialStructure::from_logical_slots(
            source
                .structure
                .logical_slots()
                .into_iter()
                .map(|slot| slot.dual()),
        );
        let mut occurrences = source
            .explicit_index_occurrences()?
            .mapped(|slot| Ok::<_, TensorInferenceError>(slot.dual()))?;
        if source.is_zero() {
            return Ok(Self::from_validated_parts(
                Atom::Zero,
                interface,
                Arc::new(occurrences),
            ));
        }
        let mut terms = Vec::new();
        for term in source.expression.expand_in(Atom::var(head)).terms() {
            let mut coefficient = Vec::new();
            let mut generators = Vec::new();
            Self::split_word_coefficient(
                head,
                term,
                &mut coefficient,
                &mut generators,
                &mut false,
            )?;
            let coefficient = coefficient
                .into_iter()
                .reduce(|left, right| Self::bracket_product(left.as_view(), right.as_view()))
                .unwrap_or_else(Atom::one);
            let coefficient = Self::infer(coefficient)?.dirac_adjoint(true)?;
            let generators = generators
                .into_iter()
                .rev()
                .map(|generator| {
                    let generator = Self::infer(generator)?;
                    let adjoint = adjoint_generator(&generator)?;
                    adjoint.explicit_index_occurrences()?;
                    let expected = generator
                        .structure
                        .logical_slots()
                        .into_iter()
                        .map(|slot| slot.dual())
                        .collect::<Vec<_>>();
                    if adjoint.structure.logical_slots() != expected {
                        return Err(TensorInferenceError::invalid(
                            "word-adjoint generator policy must preserve dual logical slots and index labels",
                        )
                        .into());
                    }
                    let mut word_free = true;
                    adjoint.expression.visitor(&mut |node| {
                        if matches!(node, AtomView::Fun(function)
                            if function.get_symbol().has_tag(&SPENSO_TAG.tensor_word))
                        {
                            word_free = false;
                        }
                        word_free
                    });
                    if !word_free {
                        return Err(TensorInferenceError::invalid(
                            "word-adjoint generator policy must return a word-free operand",
                        )
                        .into());
                    }
                    Ok(adjoint)
                })
                .collect::<Result<Vec<_>, AdjointError>>()?;
            let word = Self::word(head, &generators)?;
            let term = coefficient.multiply_word(head, &word)?;
            occurrences.branch_max(term.explicit_index_occurrences()?);
            terms.push(term.expression);
        }
        let value = Self::from_declared_parts(terms.into_iter().sum(), interface)?;
        occurrences.branch_max(value.explicit_index_occurrences()?);
        occurrences.validate()?;
        Ok(Self::from_validated_parts(
            value.expression,
            value.structure,
            Arc::new(occurrences),
        ))
    }

    /// Multiply marked tensor words using the ordinary automatic port plan.
    ///
    /// For a chosen endpoint contraction and an example, see [`ProductPlan::apply_word`].
    ///
    /// # Errors
    /// Returns planning or word-admission errors from [`ProductPlan::apply_word`].
    pub fn multiply_word(
        &self,
        head: Symbol,
        right: &Self,
    ) -> Result<Self, TensorCompositionError> {
        self.product_plan(right)?.apply_word(head, self, right)
    }

    /// A non-negative integer power in the selected word algebra.
    ///
    /// The native tensor-power owner first validates index multiplicities and
    /// the resulting interface. Multiplication then uses the same word route
    /// as an explicitly written product. Each repeated occurrence receives an
    /// independent bound-index scope while explicit free ports keep their
    /// identities; no separate nilpotency rule is used.
    /// See the example on [`ProductPlan::apply_word`].
    ///
    /// # Errors
    /// Rejects invalid word grammar or tensor powers with invalid index multiplicities.
    pub fn pow_word(&self, head: Symbol, exponent: u32) -> Result<Self, TensorCompositionError> {
        let scalar_exponent = Self::infer(Atom::num(exponent))?;
        // Each repeated bound contraction gets an independent scope below.
        // Validate only the public boundary with the native power owner; the
        // unsplit source history would falsely count its bound indices again.
        let boundary = Self::from_declared_parts(Atom::Zero, self.structure.clone())?;
        let expected = boundary.pow(&scalar_exponent)?;
        let identity = Self::word(head, &[])?;
        let initial = ProductPlan::Scalar.apply_word(head, &identity, self)?;
        let mut result = if exponent == 0 { identity } else { initial };
        let dummies = Self::reserved_dummies([self]);
        let mut used = self.expression.get_all_symbols(true);
        for index in self.explicit_indices()? {
            used.extend(index.to_atom().get_all_symbols(true));
        }
        for _ in 1..exponent {
            let scope = loop {
                let scope = dummies
                    .fresh_index()
                    .to_atom()
                    .get_symbol()
                    .expect("native dummy indices are symbols");
                if used.insert(scope) {
                    break scope;
                }
            };
            let copy = self.wrap_indices(scope, true)?;
            result = result.multiply_word(head, &copy)?;
        }
        let occurrences = Arc::clone(
            result
                .proofs
                .occurrences
                .get()
                .expect("published word history"),
        );
        ProductAlgebra::Word(head).finish(
            result.expression,
            expected.structure,
            occurrences,
            &[self],
            "pow_word",
        )
    }

    /// Extract the existing word and its commuting factors, without sorting
    /// generators or interpreting field statistics.
    fn split_word_coefficient(
        head: Symbol,
        term: AtomView<'_>,
        coefficient: &mut Vec<Atom>,
        arguments: &mut Vec<Atom>,
        found: &mut bool,
    ) -> Result<(), TensorCompositionError> {
        match term {
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    Self::split_word_coefficient(head, factor, coefficient, arguments, found)?;
                }
            }
            AtomView::Fun(function) if function.get_symbol() == SPENSO_TAG.bracket => {
                for factor in function.iter() {
                    Self::split_word_coefficient(head, factor, coefficient, arguments, found)?;
                }
            }
            AtomView::Fun(function) if function.get_symbol() == head => {
                if mem::replace(found, true) {
                    return Err(TensorCompositionError::InvalidResultInterface(
                        "each input term must contain at most one word; compose multiple words with multiply_word".into(),
                    ));
                }
                arguments.extend(function.iter().map(|argument| argument.to_owned()));
            }
            _ => {
                let mut hidden_word = false;
                term.visitor(&mut |node| {
                    if matches!(node, AtomView::Fun(function)
                        if function.get_symbol().has_tag(&SPENSO_TAG.tensor_word))
                    {
                        hidden_word = true;
                    }
                    !hidden_word
                });
                if hidden_word {
                    return Err(TensorCompositionError::InvalidResultInterface(
                        "word coefficients must be word-free; foreign heads, word powers, and opaque word-bearing wrappers are unsupported".into(),
                    ));
                }
                coefficient.push(term.to_owned());
            }
        }
        Ok(())
    }

    fn word_product_expression(
        head: Symbol,
        left: AtomView<'_>,
        right: AtomView<'_>,
    ) -> Result<Atom, TensorCompositionError> {
        // This temporary operation delegates additive distribution to
        // Symbolica. It never survives the typed construction boundary.
        let product_head = symbol!("idenso::private::word_product"; Linear);
        let product = FunctionBuilder::new(product_head)
            .add_arg(left.expand_in(Atom::var(head)))
            .add_arg(right.expand_in(Atom::var(head)))
            .finish();
        let mut error = None;
        let result = product.replace_map(|node, _, output| {
            let AtomView::Fun(function) = node else {
                return;
            };
            if function.get_symbol() != product_head || error.is_some() {
                return;
            }
            let mut coefficient = Vec::new();
            let mut arguments = Vec::new();
            for operand in function.iter() {
                if let Err(failure) = Self::split_word_coefficient(
                    head,
                    operand,
                    &mut coefficient,
                    &mut arguments,
                    &mut false,
                ) {
                    error = Some(failure);
                    return;
                }
            }
            let word = if arguments.is_empty() {
                Atom::one()
            } else {
                FunctionBuilder::new(head).add_args(arguments).finish()
            };
            let coefficient = coefficient
                .into_iter()
                .reduce(|left, right| Self::bracket_product(left.as_view(), right.as_view()))
                .unwrap_or_else(Atom::one);
            **output = Self::bracket_product(coefficient.as_view(), word.as_view());
        });
        match error {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }

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
    /// The nonzero result's logical ports follow the normalized word expression,
    /// including any operand permutation caused by its symmetry attributes.
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
        let interface = composition::interface_in_atom_order(&expression, &interface)?;
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
    use crate::dirac::AGS;
    use crate::representations::Bispinor;
    use crate::tensor::composition::{MatrixChannel, PortPair};
    use crate::tensor::{SymbolicNetExt, SymbolicNetParse};
    use crate::test_support;
    use crate::{AdjointError, bis};
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::dimension::Dimension;
    use spenso::structure::partial::PartialIndex;
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use spenso::structure::slot::{DummyAind, IsAbstractSlot, ParseableAind, SlotMatcher};
    use std::collections::{BTreeMap, HashMap};
    use std::slice;
    use symbolica::id::{Pattern, Replacement};
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

    fn canonical(tensor: &SymbolicTensor<PartialStructure>) -> Atom {
        let mut slots = SlotMatcher::default();
        let mut indices = BTreeMap::new();
        tensor.expression.visitor(&mut |part| {
            if let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(part) {
                indices.insert(slot.to_atom(), slot.rep());
                return false;
            }
            true
        });
        let canonical = tensor.expression.canonize_tensors(indices).unwrap();
        let replacements = canonical
            .dummy_indices
            .into_iter()
            .enumerate()
            .map(|(position, (source, representation))| {
                Replacement::new(
                    Pattern::Literal(source),
                    Pattern::Literal(
                        representation
                            .slot::<AbstractIndex, _>(AbstractIndex::Dummy(position))
                            .to_atom(),
                    ),
                )
            })
            .collect::<Vec<_>>();
        canonical.canonical_form.replace_multiple(&replacements)
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
    fn tensor_word_normalized_port_order_round_trips_independently_of_operand_order() {
        let head = head();
        let a = component("tensor_word_tests::order_a", 21);
        let b = component("tensor_word_tests::order_b", 22);
        for operands in [[a.clone(), b.clone()], [b, a]] {
            let word = SymbolicTensor::word(head, &operands).unwrap();
            let restored = SymbolicTensor::infer(word.expression.clone()).unwrap();
            assert!(word.same_value(&restored));
        }
        let open = unresolved("tensor_word_tests::order_unresolved");
        let word = SymbolicTensor::word(head, &[open.clone(), open]).unwrap();
        assert!(word.same_value(&SymbolicTensor::infer(word.expression.clone()).unwrap()));
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
    fn word_product_retains_cancelled_occurrences_through_native_plans() {
        let head = head();
        let a = component("word_product_history::a", 900);
        let zero = SymbolicTensor::word(head, &[a.clone(), a]).unwrap();
        assert!(zero.is_zero());
        let c = SymbolicTensor::word(head, &[component("word_product_history::c", 901)]).unwrap();
        let probe = component("word_product_history::probe", 900);
        for plan in [ProductPlan::Scalar, ProductPlan::Outer] {
            let result = plan.apply_word(head, &zero, &c).unwrap();
            assert!(matches!(
                result.multiply_word(head, &probe),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
            let contracted = ProductPlan::Contract(vec![PortPair { left: 0, right: 0 }])
                .apply_word(head, &result, &c)
                .unwrap();
            assert!(matches!(
                contracted.multiply_word(head, &probe),
                Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
            ));
        }
    }

    #[test]
    fn word_product_extends_through_a_contracted_tensor_coefficient() {
        let head = head();
        let a = component("word_product_tests::a", 201);
        let b = component("word_product_tests::b", 202);
        let c = component("word_product_tests::c", 203);
        let coefficient = component("word_product_tests::coefficient", 201);
        let ab = SymbolicTensor::word(head, &[a.clone(), b.clone()]).unwrap();
        let weighted = coefficient.multiply(&ab).unwrap();
        let c_word = SymbolicTensor::word(head, slice::from_ref(&c)).unwrap();
        let actual = weighted.multiply_word(head, &c_word).unwrap();
        let direct = SymbolicTensor::word(head, &[a, b, c]).unwrap();
        let expected = coefficient.multiply(&direct).unwrap();
        assert!(actual.same_value(&expected));
        assert_eq!(actual.rank(), 2);
    }

    #[test]
    fn word_product_merges_in_operand_order_and_preserves_exact_scalars() {
        let head = head();
        let a = component("word_product_tests::ordered_a", 211);
        let b = component("word_product_tests::ordered_b", 212);
        let a_word = SymbolicTensor::word(head, slice::from_ref(&a)).unwrap();
        let b_word = SymbolicTensor::word(head, slice::from_ref(&b)).unwrap();
        let forward = a_word.multiply_word(head, &b_word).unwrap();
        let reverse = b_word.multiply_word(head, &a_word).unwrap();
        assert_eq!(forward.expression, -&reverse.expression);
        assert_eq!(
            forward.expression,
            SymbolicTensor::word(head, &[a, b]).unwrap().expression
        );
        let scalar = SymbolicTensor::infer(Atom::num((2, 3)) + Atom::i()).unwrap();
        let weighted = scalar.multiply(&a_word).unwrap();
        assert_eq!(
            weighted.multiply_word(head, &b_word).unwrap().expression,
            scalar.expression * forward.expression
        );
    }

    #[test]
    fn word_product_uses_native_distribution_and_nilpotency() {
        let head = head();
        let a =
            SymbolicTensor::word(head, &[component("word_product_tests::linear_a", 221)]).unwrap();
        let b =
            SymbolicTensor::word(head, &[component("word_product_tests::linear_b", 221)]).unwrap();
        let c =
            SymbolicTensor::word(head, &[component("word_product_tests::linear_c", 222)]).unwrap();
        let sum = a.try_add(&b).unwrap();
        let coefficient = component("word_product_tests::linear_coefficient", 223);
        let weighted = coefficient.multiply(&sum).unwrap();
        let actual = weighted.multiply_word(head, &c).unwrap();
        let expected = coefficient
            .multiply(
                &a.multiply_word(head, &c)
                    .unwrap()
                    .try_add(&b.multiply_word(head, &c).unwrap())
                    .unwrap(),
            )
            .unwrap();
        let expected = SymbolicTensor::infer(expected.expression.expand()).unwrap();
        assert!(actual.same_value(&expected));
        assert!(a.pow_word(head, 2).unwrap().is_zero());
        let squared_sum = sum.pow_word(head, 2).unwrap();
        assert!(squared_sum.is_zero());
        assert!(squared_sum.is_scalar());
        assert!(a.pow_word(head, 0).unwrap().expression.is_one());
        assert!(a.pow_word(head, 3).is_err());
    }

    #[test]
    fn word_product_current_keeps_the_second_gamma_endpoint_open() {
        let head = head();
        let leaf = |name| {
            SymbolicTensor::infer(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_arg(Bispinor {}.new_rep(4).to_symbolic([]))
                    .finish(),
            )
            .unwrap()
        };
        let barred = SymbolicTensor::word(head, &[leaf("word_product_tests::barpsi")])
            .unwrap()
            .wrap_indices(symbol!("word_product_tests::bar_scope"), false)
            .unwrap();
        let psi = SymbolicTensor::word(head, &[leaf("word_product_tests::psi")])
            .unwrap()
            .wrap_indices(symbol!("word_product_tests::psi_scope"), false)
            .unwrap();
        let mu = Minkowski {}
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(231))
            .to_atom();
        let gamma = SymbolicTensor::infer(function!(AGS.gamma, bis!(4), bis!(4), &mu)).unwrap();
        let first = ProductPlan::Contract(vec![PortPair { left: 0, right: 0 }])
            .apply_word(head, &barred, &gamma)
            .unwrap();
        assert_eq!(first.rank(), 2);
        assert_eq!(first.structure.open_positions(), [0]);
        let current = ProductPlan::Contract(vec![PortPair { left: 0, right: 0 }])
            .apply_word(head, &first, &psi)
            .unwrap();
        assert_eq!(current.rank(), 1);
        let left = composition::port_atom(barred.structure.logical_slots()[0]);
        let right = composition::port_atom(psi.structure.logical_slots()[0]);
        let expected_gamma = SymbolicTensor::infer(function!(AGS.gamma, left, right, &mu)).unwrap();
        let expected_word = SymbolicTensor::word(head, &[barred, psi]).unwrap();
        let expected = expected_gamma.multiply(&expected_word).unwrap();
        assert!(current.same_value(&expected));
        assert!(
            SymbolicTensor::infer(current.expression.clone())
                .unwrap()
                .same_value(&current)
        );
    }

    #[test]
    fn word_product_keeps_open_even_coefficients_in_order() {
        let head = head();
        let word =
            SymbolicTensor::word(head, &[component("word_product_tests::open_word", 241)]).unwrap();
        let left = unresolved("word_product_tests::left_even");
        let right = unresolved("word_product_tests::right_even");
        let first = ProductPlan::Outer.apply_word(head, &left, &word).unwrap();
        let actual = ProductPlan::Outer.apply_word(head, &first, &right).unwrap();
        assert_eq!(actual.rank(), 3);
        assert_eq!(actual.structure.open_positions().len(), 2);
        assert!(
            SymbolicTensor::infer(actual.expression.clone())
                .unwrap()
                .same_value(&actual)
        );
        let positions = actual.structure.open_positions();
        let indexed = actual
            .reindex_interface_ports(&HashMap::from([
                (positions[0], AbstractIndex::Normal(242)),
                (positions[1], AbstractIndex::Normal(243)),
            ]))
            .unwrap();
        let expected = component("word_product_tests::left_even", 242)
            .multiply(&word)
            .unwrap()
            .multiply(&component("word_product_tests::right_even", 243))
            .unwrap();
        assert_eq!(indexed.expression, expected.expression);
    }

    #[test]
    fn word_product_rejects_lost_order_hidden_words_and_foreign_heads_even_with_zero() {
        let head = head();
        let a =
            SymbolicTensor::word(head, &[component("word_product_tests::invalid_a", 251)]).unwrap();
        let b =
            SymbolicTensor::word(head, &[component("word_product_tests::invalid_b", 252)]).unwrap();
        let lost_order = a.multiply(&b).unwrap();
        let zero = SymbolicTensor::infer(Atom::Zero).unwrap();
        assert!(lost_order.multiply_word(head, &zero).is_err());
        assert!(lost_order.pow_word(head, 0).is_err());
        let foreign_head = symbol!("word_product_tests::foreign"; Flat, Linear, Antisymmetric;
            tag = &SPENSO_TAG.tensor_word);
        let foreign = SymbolicTensor::word(
            foreign_head,
            &[component("word_product_tests::foreign_leaf", 253)],
        )
        .unwrap();
        assert!(a.multiply_word(head, &foreign).is_err());
        assert!(foreign.pow_word(head, 0).is_err());
        let closed = a
            .multiply(&component("word_product_tests::closing", 251))
            .unwrap();
        let hidden = SymbolicTensor::infer(function!(Symbol::EXP, closed.expression)).unwrap();
        assert!(hidden.multiply_word(head, &zero).is_err());
        assert!(hidden.pow_word(head, 0).is_err());
    }

    #[test]
    fn word_product_distributes_inhomogeneous_scalar_sums_without_private_heads() {
        let head = head();
        let a = component("word_product_tests::mixed_a", 254);
        let coefficient = component("word_product_tests::mixed_even", 254);
        let b = component("word_product_tests::mixed_b", 255);
        let a_word = SymbolicTensor::word(head, slice::from_ref(&a)).unwrap();
        let b_word = SymbolicTensor::word(head, slice::from_ref(&b)).unwrap();
        let scalar = SymbolicTensor::infer(Atom::num(2)).unwrap();
        let sum = coefficient
            .multiply(&a_word)
            .unwrap()
            .try_add(&scalar)
            .unwrap();
        let actual = sum.multiply_word(head, &b_word).unwrap();
        let expected = coefficient
            .multiply(&SymbolicTensor::word(head, &[a, b]).unwrap())
            .unwrap()
            .try_add(&scalar.multiply(&b_word).unwrap())
            .unwrap();
        assert!(actual.same_value(&expected));
        assert!(
            !actual
                .expression
                .contains_symbol(symbol!("idenso::private::word_product"; Linear))
        );
    }

    #[test]
    fn word_product_preserves_ranked_zero_and_owned_free_indices() {
        let head = head();
        let scope = symbol!("word_product_tests::free_scope");
        let a = SymbolicTensor::word(head, &[component("word_product_tests::zero_a", 261)])
            .unwrap()
            .wrap_indices(scope, false)
            .unwrap();
        let b = SymbolicTensor::word(head, &[component("word_product_tests::zero_b", 262)])
            .unwrap()
            .wrap_indices(scope, false)
            .unwrap();
        let zero = a.with_algebra_result(Atom::Zero).unwrap();
        let actual = zero.multiply_word(head, &b).unwrap();
        let nonzero = a.multiply_word(head, &b).unwrap();
        assert!(actual.is_zero());
        assert_eq!(actual.rank(), 2);
        assert!(actual.try_add(&nonzero).unwrap().same_value(&nonzero));
        assert_eq!(
            actual.structure.logical_slots(),
            nonzero.structure.logical_slots()
        );
    }

    #[test]
    fn word_product_honors_native_bound_index_scopes() {
        let head = head();
        let closed = SymbolicTensor::word(
            head,
            &[
                component("word_product_tests::closed_a", 271),
                component("word_product_tests::closed_b", 271),
            ],
        )
        .unwrap();
        assert!(closed.is_scalar());
        let left = closed
            .wrap_indices(symbol!("word_product_tests::left_bound"), true)
            .unwrap();
        let right = closed
            .wrap_indices(symbol!("word_product_tests::right_bound"), true)
            .unwrap();
        let result = left.multiply_word(head, &right).unwrap();
        assert!(result.is_scalar());
        assert!(!result.is_zero());
        assert!(
            SymbolicTensor::infer(result.expression.clone())
                .unwrap()
                .same_value(&result)
        );
        assert!(
            result
                .expression
                .contains_symbol(symbol!("word_product_tests::left_bound"))
        );
        assert!(
            result
                .expression
                .contains_symbol(symbol!("word_product_tests::right_bound"))
        );
    }

    #[test]
    fn word_power_reserves_scope_symbols_retained_only_in_zero_history() {
        let head = head();
        let a = component("word_power_hidden_scope::a", 920);
        let zero = SymbolicTensor::word(head, &[a.clone(), a]).unwrap();
        let AbstractIndex::Dummy(previous) = AbstractIndex::new_dummy() else {
            unreachable!()
        };
        let occupied_scope = AbstractIndex::Dummy(previous.checked_add(1).unwrap())
            .to_atom()
            .get_symbol()
            .unwrap();
        let zero = zero.wrap_indices(occupied_scope, true).unwrap();
        assert!(zero.is_zero());
        let squared = zero.pow_word(head, 2).unwrap();
        assert!(squared.is_zero());
        assert!(squared.is_scalar());
        assert_eq!(squared.explicit_indices().unwrap().count(), 2);
    }

    #[test]
    fn word_power_freshens_bound_indices_before_validating_free_port_multiplicities() {
        let head = head();
        let matrix = |name, indices: [usize; 2]| {
            SymbolicTensor::infer(
                FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                    .add_args(indices.map(|index| {
                        Minkowski {}
                            .new_rep(4)
                            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                            .to_atom()
                    }))
                    .finish(),
            )
            .unwrap()
        };
        let source = SymbolicTensor::word(
            head,
            &[
                matrix("word_power_mixed::a", [910, 911]),
                matrix("word_power_mixed::b", [910, 912]),
            ],
        )
        .unwrap();
        assert_eq!(source.rank(), 2);
        let copy = source
            .wrap_indices(symbol!("word_power_mixed::copy"), true)
            .unwrap();
        let explicit = source.multiply_word(head, &copy).unwrap();
        let squared = source.pow_word(head, 2).unwrap();
        assert!(!squared.is_zero());
        assert!(squared.is_scalar());
        assert_eq!(canonical(&squared), canonical(&explicit));
        assert!(source.pow_word(head, 3).is_err());
    }

    #[test]
    fn word_power_freshens_closed_contractions_before_exterior_normalization() {
        let head = head();
        let closed = SymbolicTensor::word(
            head,
            &[
                component("word_power_scopes::a", 281),
                component("word_power_scopes::b", 281),
            ],
        )
        .unwrap()
        .wrap_indices(
            AbstractIndex::Dummy(0).to_atom().get_symbol().unwrap(),
            true,
        )
        .unwrap();

        let left = closed
            .wrap_indices(symbol!("word_power_scopes::left"), true)
            .unwrap();
        let right = closed
            .wrap_indices(symbol!("word_power_scopes::right"), true)
            .unwrap();
        let squared = closed.pow_word(head, 2).unwrap();
        let explicit = left.multiply_word(head, &right).unwrap();
        assert!(!squared.is_zero());
        assert!(squared.is_scalar());
        assert_eq!(canonical(&squared), canonical(&explicit));
        let third = closed
            .wrap_indices(symbol!("word_power_scopes::third"), true)
            .unwrap();
        let cubed = closed.pow_word(head, 3).unwrap();
        assert!(!cubed.is_zero());
        assert_eq!(
            canonical(&cubed),
            canonical(&explicit.multiply_word(head, &third).unwrap())
        );
        assert!(closed.pow_word(head, 1).unwrap().same_value(&closed));
    }

    #[test]
    fn word_product_retains_even_compact_matrix_chains_as_coefficients() {
        let head = head();
        let matrix = SymbolicTensor::infer(function!(AGS.projp, bis!(4), bis!(4))).unwrap();
        let channel = MatrixChannel {
            input: 0,
            output: 1,
        };
        let coefficient = matrix.compose(&matrix, channel, channel).unwrap();
        let closed = |left, right, index| {
            SymbolicTensor::word(head, &[component(left, index), component(right, index)]).unwrap()
        };
        let ab = closed(
            "word_product_tests::chain_a",
            "word_product_tests::chain_b",
            281,
        );
        let cd = closed(
            "word_product_tests::chain_c",
            "word_product_tests::chain_d",
            282,
        );
        let first = coefficient.multiply_word(head, &ab).unwrap();
        let actual = first.multiply_word(head, &cd).unwrap();
        let expected = coefficient
            .multiply_word(head, &ab.multiply_word(head, &cd).unwrap())
            .unwrap();
        assert!(actual.same_value(&expected));
        assert_eq!(actual.structure.open_positions().len(), 2);
        assert!(actual.expression.contains_symbol(SPENSO_TAG.chain));
        assert!(
            !actual
                .expression
                .contains_symbol(symbol!("idenso::private::word_product"; Linear))
        );
        assert!(
            SymbolicTensor::infer(actual.expression.clone())
                .unwrap()
                .same_value(&actual)
        );
        assert!(
            ProductPlan::Compose(channel, channel)
                .apply_word(head, &matrix, &matrix)
                .is_err()
        );
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
