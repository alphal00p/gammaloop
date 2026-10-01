//! Project the open boundary of a loop-dependent tensor block. Its internal
//! contractions and vertex sums are coefficients, not projector monomials.

use super::external::ExternalProjector;
use super::*;
use idenso::shorthands::UndoShorthands;
use spenso::{
    network::parsing::ParseState,
    structure::{abstract_index::AbstractIndex, partial::PartialStructure, slot::ParseableAind},
};

impl TensorReducer {
    pub(super) fn has_tensor_sum(expression: AtomView<'_>) -> bool {
        let mut found = false;
        expression.visitor(&mut |value| {
            if matches!(value, AtomView::Add(_))
                && value.contains_symbol(Minkowski::selfless_symbol())
            {
                found = true;
            }
            !found
        });
        found
    }

    pub(super) fn reduce_blocks(
        &self,
        source: &SymbolicTensor<PartialStructure>,
        cooking: &CookSettings,
        external: Option<&ExternalProjector>,
    ) -> Result<TensorReduction, TensorReductionError> {
        // Refuse oversized selected powers before the collector instantiates
        // their independent dummy scopes. Scalar invariant powers stay opaque.
        let mut oversized_power = false;
        source.expression().visitor(&mut |value| {
            if let AtomView::Pow(power) = value
                && let AtomView::Add(sum) = power.get_base()
                && let Ok(exponent) = usize::try_from(power.get_exp())
            {
                oversized_power |= sum.iter().any(|term| {
                    self.parse_monomial(cooking.uncook(term).as_view())
                        .is_ok_and(|term| {
                            term.integrated.len().saturating_mul(exponent)
                                > OrthogonalWeingarten::MAX_RANK
                        })
                });
            }
            !oversized_power
        });
        if oversized_power {
            return Err(TensorReductionError::TensorAlgebra(
                "selected coefficient collection did not complete: tensor power exceeds the projector rank limit".into(),
            ));
        }
        // Only loop/spectator dots need open ports. Invariant dots, including
        // inverse Gram factors, must remain scalar coefficients.
        let mut admission_error = None;
        let indexed = source.expression().replace_map(|value, _, out| {
            let base = match value {
                AtomView::Pow(power) => power.get_base(),
                base => base,
            };
            if !matches!(base, AtomView::Fun(f) if f.get_symbol() == SPENSO_TAG.dot) {
                return;
            }
            match self.parse_monomial(cooking.uncook(value).as_view()) {
                Ok(term) if !term.integrated.is_empty() => {
                    match value.undo_dots::<AbstractIndex>() {
                        Ok(value) => **out = value,
                        Err(error) => {
                            admission_error =
                                Some(TensorReductionError::TensorAlgebra(error.to_string()))
                        }
                    }
                }
                Ok(_) => **out = value.to_owned(),
                Err(error) => admission_error = Some(error),
            }
        });
        if let Some(error) = admission_error.take() {
            return Err(error);
        }
        let indexed = SymbolicTensor::infer(indexed)
            .map_err(|e| TensorReductionError::TensorAlgebra(e.to_string()))?;
        let rows = indexed
            .coefficient_list_with(|value| {
                let value = cooking.uncook(value);
                match self.indexed_vector(value.as_view()) {
                    Ok(Some(vector)) => !self.is_integrated(&vector),
                    Err(error) => {
                        admission_error = Some(error);
                        false
                    }
                    Ok(None) => {
                        // Metrics rotate covariantly with the loop block. Other
                        // Lorentz tensors (including spin matrices) stay fixed.
                        matches!(value.as_view(), AtomView::Fun(f) if f.get_symbol() != ETS.metric && f.get_symbol() != SPENSO_TAG.dot)
                            && value.contains_symbol(Minkowski::selfless_symbol())
                    }
                }
            })
            .map_err(|e| TensorReductionError::TensorAlgebra(e.to_string()))?;
        if let Some(error) = admission_error {
            return Err(error);
        }
        let state = ParseState::<AbstractIndex>::default();
        state.reserve_indices(indexed.expression().as_view());
        let mut terms = Vec::new();
        let mut max_rank = 0;
        for (outside, inside) in rows {
            let boundary = PartialStructure::from_logical_slots(
                inside
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .filter(|slot| slot.rep().rep == LibraryRep::from(Minkowski {})),
            );
            let slots = boundary
                .slots()
                .map_err(|e| TensorReductionError::TensorAlgebra(e.to_string()))?;
            let rank = slots.len();
            max_rank = max_rank.max(rank);
            if rank > OrthogonalWeingarten::MAX_RANK {
                return Err(TensorReductionError::UnsupportedRank {
                    rank,
                    maximum: OrthogonalWeingarten::MAX_RANK,
                });
            }
            let slots = slots.iter().map(|slot| slot.to_atom()).collect::<Vec<_>>();
            let fresh = (0..rank)
                .map(|_| minkowski_slot(&self.dimension, &state.fresh_index().to_atom()))
                .map(|slot| cooking.cook(slot.as_view()))
                .collect::<Vec<_>>();
            // Inner contractions own distinct dummies even when retained inside
            // a factored sum; they must not reconnect to the outside projector.
            let renamed = inside.expression().replace_map(|value, _, out| {
                if let Some(position) = slots.iter().position(|slot| slot.as_view() == value) {
                    **out = fresh[position].clone();
                }
            });
            let contract = |expression: Atom| {
                SymbolicTensor::infer(expression)
                    .and_then(|value| {
                        value.contract(ContractSettings {
                            expand: false,
                            collect_chains: false,
                            collect_traces: false,
                            ..Default::default()
                        })
                    })
                    .and_then(|value| value.to_dots())
                    .map(|value| value.into_expression())
                    .map_err(|e| TensorReductionError::TensorAlgebra(e.to_string()))
            };
            let probe = SPENSO_TAG.rank_one_tensor_symbol("FeynKit::BlockProbe");
            let probes = (0..)
                .map(|i| {
                    FunctionBuilder::new(probe)
                        .add_arg(i)
                        .add_arg(
                            FunctionBuilder::new(Minkowski::selfless_symbol())
                                .add_arg(&self.dimension)
                                .finish(),
                        )
                        .finish()
                })
                .filter(|probe| !self.external_vectors.contains(probe))
                .take(rank)
                .collect::<Vec<_>>();
            let monomial = TensorMonomial {
                scalar: Atom::one(),
                integrated: probes.clone(),
                outside: slots
                    .iter()
                    .map(|slot| {
                        let slot = cooking.uncook(slot.as_view());
                        let AtomView::Fun(function) = slot.as_view() else {
                            unreachable!()
                        };
                        let index = function.iter().last().unwrap().to_owned();
                        Outside::Index { slot, index }
                    })
                    .collect(),
            };
            let projected = if let Some(external) = external {
                self.reduce_external_monomial(monomial, external)?
            } else if rank % 2 == 0 {
                self.reduce_monomial(monomial)?
            } else {
                Vec::new()
            };
            let fresh = fresh
                .iter()
                .map(|slot| cooking.uncook(slot.as_view()))
                .collect::<Vec<_>>();
            for term in projected {
                // Apply the existing projector to a general covariant block:
                // every probe occurrence denotes one of its open tensor ports.
                let projector = term.tensor.replace_map(|value, _, out| {
                    let AtomView::Fun(function) = value else {
                        return;
                    };
                    if function.get_symbol() != SPENSO_TAG.dot {
                        return;
                    }
                    let args = function.iter().collect::<Vec<_>>();
                    if args.len() != 2 {
                        return;
                    }
                    let positions = args
                        .iter()
                        .map(|arg| probes.iter().position(|probe| probe.as_view() == *arg))
                        .collect::<Vec<_>>();
                    **out = match positions.as_slice() {
                        [Some(left), Some(right)] => FunctionBuilder::new(ETS.metric)
                            .add_arg(&fresh[*left])
                            .add_arg(&fresh[*right])
                            .finish(),
                        [Some(position), None] | [None, Some(position)] => {
                            let other = if positions[0].is_none() {
                                args[0]
                            } else {
                                args[1]
                            };
                            let AtomView::Fun(slot) = fresh[*position].as_view() else {
                                unreachable!()
                            };
                            indexed_vector(
                                &other.to_owned(),
                                &slot.iter().last().unwrap().to_owned(),
                            )
                        }
                        _ => return,
                    };
                });
                let tensor =
                    contract(&renamed * outside.expression() * cooking.cook(projector.as_view()))?;
                if !tensor.is_zero() {
                    terms.push(TensorReductionTerm {
                        coefficient: term.coefficient,
                        tensor: cooking.uncook(tensor.as_view()),
                        integrated_orbit: None,
                        projector_orbit: None,
                    });
                }
                if terms.len() > self.output_term_limit {
                    return Err(TensorReductionError::OutputLimit {
                        terms: terms.len(),
                        limit: self.output_term_limit,
                    });
                }
            }
        }
        let fully_contracted = terms.is_empty()
            || source
                .structure()
                .logical_slots()
                .iter()
                .all(|slot| slot.rep().rep != LibraryRep::from(Minkowski {}));
        Ok(TensorReduction {
            terms,
            max_rank,
            fully_contracted,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{parse, symbol};

    #[test]
    fn loop_vector_sums_are_retained_under_a_rank_two_projector() {
        let k = spenso::vector_symbol!("block_test::k");
        let q = spenso::vector_symbol!("block_test::q");
        let p = spenso::vector_symbol!("block_test::p");
        let d = symbol!("block_test::D").to_atom();
        let mu = symbol!("block_test::mu").to_atom();
        let nu = symbol!("block_test::nu").to_atom();
        let vector = |head, index| {
            FunctionBuilder::new(head)
                .add_arg(minkowski_slot(&d, index))
                .finish()
        };
        let input = (vector(k, &mu) + vector(q, &mu))
            * (vector(k, &nu) + vector(q, &nu))
            * vector(p, &mu)
            * vector(p, &nu);
        let result = TensorReducer::new(d.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .with_output_term_limit(1)
            .reduce(input.as_view())
            .unwrap();
        assert_eq!(result.terms().len(), 1);
        assert!(result.is_fully_contracted());
        assert!(TensorReducer::has_tensor_sum(result.expression().as_view()));
        // Independently contract this small rank-two identity to scalar dots.
        let actual = SymbolicTensor::infer(result.expression())
            .unwrap()
            .contract(ContractSettings::default())
            .unwrap()
            .to_dots()
            .unwrap()
            .into_expression();
        let compact = |head| {
            FunctionBuilder::new(head)
                .add_arg(
                    FunctionBuilder::new(Minkowski::selfless_symbol())
                        .add_arg(&d)
                        .finish(),
                )
                .finish()
        };
        let [k, q, p] = [k, q, p].map(compact);
        let expected = (dot(&k, &k) + Atom::num(2) * dot(&k, &q) + dot(&q, &q)) * dot(&p, &p) / d;
        assert!(
            (&actual - &expected).expand().is_zero(),
            "actual: {actual}; expected: {expected}"
        );
    }

    #[test]
    fn mixed_loop_and_fixed_sums_keep_the_fixed_rank_two_term() {
        let k = spenso::vector_symbol!("block_mixed::k");
        let p = spenso::vector_symbol!("block_mixed::p");
        let input = parse!(
            "(block_mixed::k(spenso::mink(D,mu))+block_mixed::p(spenso::mink(D,mu)))
            *(block_mixed::k(spenso::mink(D,nu))+block_mixed::p(spenso::mink(D,nu)))"
        );
        let result = TensorReducer::new(parse!("D"))
            .with_integrated_head(k)
            .reduce(input.as_view())
            .unwrap();
        let expected = parse!(
            "spenso::dot(block_mixed::k(spenso::mink(D)),block_mixed::k(spenso::mink(D)))
            *spenso::g(spenso::mink(D,mu),spenso::mink(D,nu))/D
            +block_mixed::p(spenso::mink(D,mu))*block_mixed::p(spenso::mink(D,nu))"
        );
        assert!((result.expression() - expected).expand().is_zero());
        assert!(!result.is_fully_contracted());
        let _ = p;
    }
}
