//! Collect exact terminal polynomials before materializing atoms.
//! Source scalar sums retain their boundaries in selected products and power bases.

use std::collections::HashMap;
use std::iter;
use std::mem;
#[cfg(test)]
use std::slice;
use std::sync::Arc;

use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::float::Complex;
use symbolica::domains::rational::{Q, Rational, RationalField};
#[cfg(test)]
use symbolica::parse;
use symbolica::poly::PolyVariable;
use symbolica::poly::polynomial::MultivariatePolynomial;
use symbolica_utils::PolynomialExpressionExt;

use crate::color::simplify::ColorAlgebraSimplifier;
#[cfg(feature = "reference-cases")]
use crate::reference_cases::timing;
#[cfg(feature = "reference-cases")]
use crate::reference_cases::timing::Phase;
use crate::tensor::inference::{InterfaceInference, TensorInferenceError};
use crate::tensor::simplification::observation::DomainObservations;
use crate::tensor::{SymbolicTensor, TensorProofs};

pub(crate) type CollectedPolynomial = MultivariatePolynomial<RationalField, u32>;

// Policy exception to treating dependencies' `#[doc(hidden)]` items as private:
// `PolyVariable::Temporary` is limited to this collector's opaque source sums,
// preserved color sums and imaginary unit. Public alternatives would register
// symbols (global symbol-table writes and permanently registered names) or repurpose
// `Function` or `Power` variables with a different meaning. These index-based variables
// all map back to their atoms before output is materialized; none escapes.
// An upstream public opaque-atom variable would remove this exception.
pub(crate) struct CollectedOutput<'a> {
    source: &'a DomainObservations,
    proofs: TensorProofs,
    atoms: [HashMap<Atom, CollectedPolynomial>; 2],
    opaque: Vec<Atom>,
    imaginary_unit: Option<PolyVariable>,
    preserve_sum: Option<fn(AtomView<'_>) -> bool>,
    scalar_sums: HashMap<Atom, CollectedPolynomial>,
    preserved_sums: HashMap<Atom, CollectedPolynomial>,
    variables: Arc<Vec<PolyVariable>>,
    tensors: HashMap<PolyVariable, Option<bool>>,
    refused: bool,
    emitted: HashMap<Atom, CollectedPolynomial>,
}

struct CollectedSector {
    selected: SymbolicTensor<PartialStructure>,
    polynomial: CollectedPolynomial,
}

impl<'a> CollectedOutput<'a> {
    pub(crate) fn new(source: &'a SymbolicTensor<PartialStructure>) -> Self {
        Self {
            source: source.reduction_observations(),
            proofs: source.proofs.clone(),
            atoms: Default::default(),
            opaque: Vec::new(),
            imaginary_unit: None,
            preserve_sum: None,
            scalar_sums: HashMap::new(),
            preserved_sums: HashMap::new(),
            variables: Arc::default(),
            tensors: HashMap::new(),
            refused: false,
            emitted: HashMap::new(),
        }
    }

    /// Keep domain sum factors as main variables, outside coefficient
    /// normalization. Domains use this to retain terminal decompositions whose
    /// product expansion would discard their fixed-point certificate.
    pub(crate) fn with_preserved_sums(mut self, preserve: fn(AtomView<'_>) -> bool) -> Self {
        self.preserve_sum = Some(preserve);
        self
    }

    /// Attach the boundary supplied by a leaf constructor before polynomial
    /// conversion can hide that leaf among generated coefficients.
    pub(crate) fn polynomial_with_interface(
        &mut self,
        value: AtomView<'_>,
        interface: &PartialStructure,
    ) -> CollectedPolynomial {
        let polynomial = self.polynomial(value, false);
        for variable in polynomial.variables().iter() {
            if !matches!(variable, PolyVariable::Temporary(_))
                && variable.to_atom().as_view() == value
            {
                self.tensors.insert(
                    variable.clone(),
                    Some(!interface.logical_slots().is_empty()),
                );
            }
        }
        polynomial
    }

    pub(crate) fn metric_polynomial(
        &mut self,
        metric: AtomView<'_>,
        endpoints: [AtomView<'_>; 2],
    ) -> CollectedPolynomial {
        if let Some(polynomial) = self.atoms[0].get(metric.get_data()) {
            return polynomial.clone();
        }
        let interface = self
            .proofs
            .metric_interface(endpoints)
            .expect("admitted trace endpoints establish their generated metric boundary");
        self.polynomial_with_interface(metric, &interface)
    }

    // Selected exterior coefficients can be collected across trace sectors.
    // Vector weights and power bases retain their original scalar-sum boundaries.
    pub(crate) fn polynomial(
        &mut self,
        value: AtomView<'_>,
        collect_sums: bool,
    ) -> CollectedPolynomial {
        self.polynomial_in(value, collect_sums, None)
    }

    fn polynomial_in(
        &mut self,
        value: AtomView<'_>,
        collect_sums: bool,
        row: Option<&DomainObservations>,
    ) -> CollectedPolynomial {
        if let Some(polynomial) = self.atoms[usize::from(collect_sums)].get(value.get_data()) {
            return polynomial.clone();
        }
        let mut polynomial = if !collect_sums
            && matches!(value, AtomView::Add(_))
            && self
                .source
                .region(value)
                .or_else(|| row.and_then(|observed| observed.region(value)))
                .is_some_and(|region| region.certifies_scalar_interface_syntax())
        {
            // The source or the materialized row observed this scalar sum
            // as a coefficient. Preserve its factorization.
            self.opaque_sum(value, false)
        } else {
            match value {
                AtomView::Add(sum) => sum
                    .iter()
                    .fold(CollectedPolynomial::new_zero(&Q), |sum, term| {
                        sum + self.polynomial_in(term, collect_sums, row)
                    }),
                AtomView::Mul(product) => {
                    product
                        .iter()
                        .fold(CollectedPolynomial::new_one(&Q), |product, factor| {
                            let factor = if matches!(factor, AtomView::Add(_))
                                && self.preserve_sum.is_some_and(|preserve| preserve(factor))
                            {
                                self.opaque_sum(factor, true)
                            } else {
                                self.polynomial_in(factor, collect_sums, row)
                            };
                            self.product(product, factor)
                        })
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    if let Ok(exponent) = i64::try_from(exponent)
                        && let Ok(exponent) = u32::try_from(exponent)
                    {
                        // Contraction can combine two source powers. Reuse the
                        // opaque base even when the resulting power is new.
                        let base = self.polynomial_in(base, false, row);
                        if exponent <= i32::MAX as u32
                            && (0..base.nvars()).all(|variable| {
                                base.degree(variable)
                                    .checked_mul(exponent)
                                    .is_some_and(|degree| degree <= i32::MAX as u32)
                            })
                        {
                            base.pow(exponent as usize)
                        } else {
                            self.refused = true;
                            CollectedPolynomial::new_zero(&Q)
                        }
                    } else {
                        CollectedPolynomial::new(
                            &Q,
                            Some(1),
                            Arc::new(vec![PolyVariable::Power(value.to_owned())]),
                        )
                        .monomial(Rational::one(), vec![1])
                    }
                }
                _ => match value.try_to_polynomial(&Q, None) {
                    Ok(polynomial) => polynomial,
                    Err(_) => {
                        if let Ok(number) = Complex::<Rational>::try_from(value) {
                            // Keep Q as the coefficient ring. The imaginary unit
                            // is an opaque numerical atom, never a registered symbol.
                            // Materialization normalizes its integer powers as well.
                            let variable = self.imaginary_unit.get_or_insert_with(|| {
                                let variable = PolyVariable::Temporary(self.opaque.len());
                                self.opaque.push(Atom::i());
                                variable
                            });
                            CollectedPolynomial::new_one(&Q).mul_coeff(number.re)
                                + CollectedPolynomial::new(
                                    &Q,
                                    Some(1),
                                    Arc::new(vec![variable.clone()]),
                                )
                                .monomial(number.im, vec![1])
                        } else {
                            // Rounded or unsupported coefficient domains retain
                            // established atom arithmetic, including its grouping.
                            self.refused = true;
                            CollectedPolynomial::new_zero(&Q)
                        }
                    }
                },
            }
        };
        if !self.refused {
            // Every temporary coefficient uses a prefix of the same ordering.
            // Unification can then extend exponent rows without re-sorting a
            // large polynomial every time another pairing introduces a variable.
            for variable in polynomial.variables().iter() {
                if !self.variables.contains(variable) {
                    Arc::make_mut(&mut self.variables).push(variable.clone());
                }
            }
            CollectedPolynomial::new(&Q, None, self.variables.clone())
                .unify_variables(&mut polynomial);
            self.atoms[usize::from(collect_sums)].insert(value.to_owned(), polynomial.clone());
        }
        polynomial
    }

    // An opaque sum has one variable in every row of this collection, so
    // independently completed rows can cancel before materializing the atom.
    fn opaque_sum(&mut self, value: AtomView<'_>, tensor: bool) -> CollectedPolynomial {
        let sums = if tensor {
            &mut self.preserved_sums
        } else {
            &mut self.scalar_sums
        };
        if let Some(polynomial) = sums.get(value.get_data()) {
            return polynomial.clone();
        }
        let variable = PolyVariable::Temporary(self.opaque.len());
        self.opaque.push(value.to_owned());
        if tensor {
            self.tensors.insert(variable.clone(), Some(true));
        }
        let polynomial = CollectedPolynomial::new(&Q, Some(1), Arc::new(vec![variable]))
            .monomial(Rational::one(), vec![1]);
        sums.insert(value.to_owned(), polynomial.clone());
        polynomial
    }

    pub(crate) fn product(
        &mut self,
        mut left: CollectedPolynomial,
        mut right: CollectedPolynomial,
    ) -> CollectedPolynomial {
        if self.refused {
            return CollectedPolynomial::new_zero(&Q);
        }
        // Reuse the owned inputs when aligning variables. Multiplication would
        // otherwise clone both just to unify their coefficient contexts.
        left.unify_variables(&mut right);
        // Symbolica stores u32 exponents but its arithmetic uses signed bounds.
        // Ordinary trace degrees are small; only unusually large input powers
        // need the more expensive per-variable degree check.
        if !left
            .exponents
            .iter()
            .chain(&right.exponents)
            .all(|&exponent| exponent <= i32::MAX as u32 / 2)
            && (0..left.nvars()).any(|variable| {
                u64::from(left.degree(variable)) + u64::from(right.degree(variable))
                    > i32::MAX as u64
            })
        {
            self.refused = true;
            return CollectedPolynomial::new_zero(&Q);
        }
        &left * &right
    }

    fn emit_coefficient(&self, polynomial: &CollectedPolynomial) -> Atom {
        let coefficient = polynomial
            .to_expression_with_variables(|variable| {
                Some(match variable {
                    PolyVariable::Temporary(index) => self.opaque[*index].clone(),
                    variable => variable.to_atom(),
                })
            })
            .expect("every temporary has its original scalar scope");
        // Distribute numerical coefficients using Symbolica's own normalizer.
        // This can cancel 2*(x+y)-2*x-2*y while leaving products and powers of
        // unrelated scalar sums factorized. Tensor factors are emitted outside.
        coefficient.expand_num()
    }

    /// A colour identity can emit both a source scalar sum and its individual
    /// terms. Test their exact polynomial relation without opening any other
    /// source factor or changing the grouping of a surviving result.
    fn cancels_with_source_sum(&self, value: &CollectedPolynomial) -> bool {
        if self.opaque.is_empty() || value.nterms() < 2 {
            return false;
        }
        let elementary =
            |atom| matches!(atom, AtomView::Num(_) | AtomView::Var(_) | AtomView::Fun(_));
        for (position, variable) in value.variables().iter().enumerate() {
            let PolyVariable::Temporary(index) = variable else {
                continue;
            };
            let AtomView::Add(sum) = self.opaque[*index].as_view() else {
                continue;
            };
            let degree = value.degree(position);
            if degree == 0
                || self.tensors.get(variable) == Some(&Some(true))
                || !sum.iter().all(|term| match term {
                    AtomView::Mul(product) => product.iter().all(elementary),
                    factor => elementary(factor),
                })
            {
                continue;
            }
            let Ok(expanded) = sum
                .as_view()
                .try_to_polynomial::<_, u32>(&Q, value.variables())
            else {
                continue;
            };
            // Only reconcile monomials in the existing basis. Bound the terms
            // per substitution by the size of the collected input and keep
            // exponent arithmetic within Symbolica's signed range.
            if expanded.variables() == value.variables()
                && expanded
                    .nterms()
                    .checked_pow(degree)
                    .is_some_and(|n| n <= value.nterms())
                && (0..value.nvars()).all(|i| {
                    u64::from(expanded.degree(i)) * u64::from(degree) + u64::from(value.degree(i))
                        <= i32::MAX as u64
                })
                && value.replace_with_poly(position, &expanded).is_zero()
            {
                return true;
            }
        }
        false
    }

    pub(crate) fn finish(&mut self, value: CollectedPolynomial) -> Option<Atom> {
        #[cfg(feature = "reference-cases")]
        let _phase = timing::scope(Phase::OutputNormalization);
        if self.refused {
            return None;
        }
        let mut tensors = Vec::new();
        for (position, variable) in value.variables().iter().enumerate() {
            if value.degree(position) == 0 {
                continue;
            }
            if matches!(variable, PolyVariable::Temporary(_)) {
                if self.tensors.get(variable) == Some(&Some(true)) {
                    tensors.push(position);
                }
                continue;
            }
            if (*self.tensors.entry(variable.clone()).or_insert_with(|| {
                let atom = variable.to_atom();
                self.proofs
                    .admitted_interface(atom.as_view())
                    .or_else(|| ColorAlgebraSimplifier::produced_interface(atom.as_view()))
                    .map(|interface| !interface.logical_slots().is_empty())
            }))? {
                tensors.push(position);
            }
        }
        let polynomial = value.to_polynomial_in(&tensors);
        // Both coefficients and main variables have custom atom boundaries.
        // Symbolica's coefficient mapper cannot emit temporary main variables.
        let variables = polynomial
            .variables()
            .iter()
            .map(|variable| match variable {
                PolyVariable::Temporary(index) => self.opaque[*index].clone(),
                variable => variable.to_atom(),
            })
            .collect::<Vec<_>>();
        let expression = Atom::add_many(polynomial.into_iter().map(|term| {
            Atom::mul_many(
                iter::once(self.emit_coefficient(term.coefficient)).chain(
                    variables
                        .iter()
                        .zip(term.exponents)
                        .filter(|(_, exponent)| **exponent != 0)
                        .map(|(variable, exponent)| variable.pow(i64::from(*exponent))),
                ),
            )
        }));
        // Numerical distribution can prove a cancellation against an opened
        // copy of a preserved sum. Accept that zero without changing any
        // surviving decomposition's arithmetic contexts.
        let expression = if (!self.preserved_sums.is_empty() && expression.expand_num().is_zero())
            || (!expression.is_zero() && self.cancels_with_source_sum(&value))
        {
            Atom::Zero
        } else {
            expression
        };
        let value = if expression.is_zero() {
            CollectedPolynomial::new_zero(&Q)
        } else {
            value
        };
        self.emitted.insert(expression.clone(), value);
        Some(expression)
    }

    pub(crate) fn begin(&mut self) {
        self.refused = false;
    }

    /// Complete a terminal atom result through the same arithmetic used by
    /// recurrence producers. Refusal leaves the caller's atom arithmetic intact.
    fn collect_atom(
        &mut self,
        expression: &Atom,
        row: Option<&DomainObservations>,
    ) -> Option<Atom> {
        self.begin();
        // A row can publish a newly scalar sum. Every enclosing power/product
        // can therefore have a different polynomial at the same expression bytes.
        let source_atoms = row.map(|_| mem::take(&mut self.atoms));
        let polynomial = self.polynomial_in(expression.as_view(), false, row);
        if let Some(cache) = source_atoms {
            self.atoms = cache;
        }
        self.finish(polynomial)
    }

    pub(crate) fn collect_tensor(
        &mut self,
        source: &SymbolicTensor<PartialStructure>,
    ) -> Option<Atom> {
        // This emission may have introduced leaves absent from the original
        // collector input. Use the actual producer's proof owner for the row.
        let previous = mem::replace(&mut self.proofs, source.proofs.clone());
        let observations = source.work.observations().map(AsRef::as_ref);
        let expression = self.collect_atom(&source.expression, observations);
        self.proofs = previous;
        expression
    }

    pub(crate) fn retain_product(&mut self, expression: &Atom) {
        if self.emitted.contains_key(expression) {
            return;
        }
        let polynomial = self.polynomial(expression.as_view(), false);
        if !self.refused {
            self.emitted.insert(expression.clone(), polynomial);
        }
    }

    /// Combine only certified generated sectors. The collector retains every
    /// unrelated spectator outside this frontier and does not open it here.
    pub(crate) fn collect(
        &mut self,
        states: &mut Vec<(
            SymbolicTensor<PartialStructure>,
            SymbolicTensor<PartialStructure>,
        )>,
    ) -> Result<(), TensorInferenceError> {
        #[cfg(feature = "reference-cases")]
        let _phase = timing::scope(Phase::CoefficientMerge);
        let mut generated: Vec<CollectedSector> = Vec::new();
        let mut retained = Vec::new();
        for (selected, coefficient) in mem::take(states) {
            let Some(polynomial) = self.emitted.get(&selected.expression) else {
                // Some terminal identities normalize to a metric before the
                // bound recurrence runs. Their selected coefficients still
                // receive the requested numerical-factor distribution.
                let expression = coefficient.expression.expand_num();
                let coefficient = if expression == coefficient.expression {
                    coefficient
                } else {
                    let observations = coefficient
                        .reduction_observations()
                        .updated(expression.as_view());
                    coefficient.with_identity_result(expression, Some(observations))?
                };
                retained.push((selected, coefficient));
                continue;
            };
            let facts = coefficient.reduction_observations();
            if facts.identity_candidates() != 0
                || !(facts.excludes_internal_connections()
                    || (coefficient.is_scalar() && facts.can_certify_scalar()))
                || !selected.reduction_observations().is_terminal_polynomial()
            {
                retained.push((selected, coefficient));
                continue;
            }
            let polynomial = polynomial.clone();
            let interface = InterfaceInference::merge_explicit_interface_sequence(&[
                selected.structure.clone(),
                coefficient.structure.clone(),
            ])?;
            if interface.logical_slots().len()
                != selected.structure.logical_slots().len()
                    + coefficient.structure.logical_slots().len()
            {
                retained.push((selected, coefficient));
                continue;
            }
            self.refused = false;
            let coefficient_polynomial = self.polynomial(coefficient.expression.as_view(), true);
            for variable in coefficient_polynomial.variables().iter() {
                if !matches!(variable, PolyVariable::Temporary(_))
                    && let Some(interface) = coefficient
                        .proofs
                        .admitted_interface(variable.to_atom().as_view())
                {
                    self.tensors.insert(
                        variable.clone(),
                        Some(!interface.logical_slots().is_empty()),
                    );
                }
            }
            if self.refused {
                retained.push((selected, coefficient));
                continue;
            }
            let polynomial = self.product(polynomial, coefficient_polynomial);
            if self.refused {
                retained.push((selected, coefficient));
                continue;
            }
            if let Some(position) = generated.iter().position(|previous| {
                InterfaceInference::additive_interfaces_match(
                    &previous.selected.structure,
                    &interface,
                )
            }) {
                let previous = &mut generated[position].polynomial;
                *previous = mem::replace(previous, CollectedPolynomial::new_zero(&Q)) + polynomial;
            } else {
                let selected = if coefficient.is_scalar() {
                    selected
                } else {
                    SymbolicTensor::collection_product(&[selected, coefficient])?
                };
                generated.push(CollectedSector {
                    selected,
                    polynomial,
                });
            }
        }
        self.refused = false;
        for CollectedSector {
            selected,
            polynomial,
        } in generated
        {
            let expression = self
                .finish(polynomial)
                .expect("certified polynomial variables retain their interfaces");
            let observations = DomainObservations::terminal_polynomial(&selected.structure, true);
            let selected = selected.with_identity_result(expression, Some(observations))?;
            let coefficient = SymbolicTensor::from_validated_parts(
                Atom::one(),
                PartialStructure::from_logical_slots([]),
                Arc::default(),
            );
            retained.push((selected, coefficient));
        }
        *states = retained;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generated_sum_factors_use_the_current_rows_observations() {
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let row = SymbolicTensor::infer(parse!("row_scope_weight * (row_scope_x + row_scope_y)"))
            .unwrap();
        let facts = row.reduction_observations().as_ref();
        for row_first in [true, false] {
            let mut output = CollectedOutput::new(&source);
            for use_row in [row_first, !row_first, row_first] {
                let result = output
                    .collect_atom(&row.expression, use_row.then_some(facts))
                    .unwrap();
                let expected = if use_row {
                    row.expression.clone()
                } else {
                    row.expression.expand()
                };
                assert_eq!(result, expected);
            }
            assert_eq!(output.opaque.len(), 1);
        }
    }

    #[test]
    fn exact_complex_products_collect_before_and_after_materialization() {
        let x = parse!("complex_collection_x");
        let y = parse!("complex_collection_y");
        for remainder in [0, 1] {
            let expression = (&x + Atom::i() * &y) * (&x - Atom::i() * &y) - x.pow(2) - y.pow(2)
                + Atom::num(remainder);
            let source = SymbolicTensor::infer(Atom::one()).unwrap();
            let mut output = CollectedOutput::new(&source);
            let result = output.collect_atom(&expression, None).unwrap();
            assert_eq!(result, Atom::num(remainder));
            assert_eq!(output.collect_atom(&result, None), Some(result));
        }
    }

    #[test]
    fn preserved_sum_factors_keep_scalar_weights_collectable() {
        let source = SymbolicTensor::infer(parse!(
            "(preserved_a+preserved_b)*(preserved_c+preserved_d)"
        ))
        .unwrap();
        let mut output = CollectedOutput::new(&source)
            .with_preserved_sums(|value| !matches!(value, AtomView::Num(_)));
        assert_eq!(
            output.collect_atom(&source.expression, None),
            Some(source.expression.clone())
        );
        let sum = parse!("preserved_a+preserved_b");
        let cancel = Atom::num(2) * &sum - parse!("2*preserved_a+2*preserved_b");
        assert_eq!(output.collect_atom(&cancel, None), Some(Atom::Zero));
    }

    #[test]
    fn unrepresentable_power_refuses_without_expanding_its_scalar_base() {
        let source = SymbolicTensor::infer(parse!("overflow_x + overflow_y")).unwrap();
        let mut output = CollectedOutput::new(&source);
        let expression = source.expression.pow(i64::from(i32::MAX) + 1);
        let polynomial = output.polynomial(expression.as_view(), false);
        assert!(output.refused);
        assert!(output.finish(polynomial).is_none());
        assert_eq!(output.opaque, slice::from_ref(&source.expression));
        assert!(output.emitted.is_empty());
    }

    #[test]
    fn product_checks_all_monomials_after_aligning_variables() {
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let mut output = CollectedOutput::new(&source);
        let left = CollectedPolynomial::from_coefficient_list(
            vec![Rational::one(), Rational::one()],
            vec![0, i32::MAX as u32, 1, 0],
            vec![PolyVariable::Temporary(0), PolyVariable::Temporary(1)].into(),
            &Q,
        );
        let right = CollectedPolynomial::from_coefficient_list(
            vec![Rational::one()],
            vec![1, 0],
            vec![PolyVariable::Temporary(1), PolyVariable::Temporary(0)].into(),
            &Q,
        );
        let polynomial = output.product(left, right);
        assert!(output.refused);
        assert!(output.finish(polynomial).is_none());
    }

    #[test]
    fn refused_finalizer_product_retains_the_exact_selected_state() {
        let source = SymbolicTensor::infer(parse!("retained_x + retained_y")).unwrap();
        let mut output = CollectedOutput::new(&source);
        let selected = SymbolicTensor::from_validated_parts(
            source.expression.pow(i64::from(i32::MAX)),
            PartialStructure::from_logical_slots([]),
            Arc::clone(source.proofs.occurrences.get().unwrap()),
        );
        selected
            .work
            .get_or_init()
            .observations
            .set(Arc::new(DomainObservations::terminal_polynomial(
                &selected.structure,
                false,
            )))
            .unwrap();
        output.retain_product(&selected.expression);
        assert!(!output.refused);
        let coefficient = SymbolicTensor::infer(source.expression.pow(2)).unwrap();
        let original = (selected, coefficient);
        let mut states = vec![original.clone()];
        output.collect(&mut states).unwrap();
        assert_eq!(states, [original]);
    }
}
