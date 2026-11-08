//! Project the open boundary of a loop-dependent tensor block. Its internal
//! contractions and vertex sums are coefficients, not projector monomials.

use super::external::ExternalProjector;
use super::{
    Outside, TensorMonomial, TensorReducer, TensorReduction, TensorReductionError,
    TensorReductionTerm, indexed_vector,
};
use idenso::shorthands::UndoShorthands;
use idenso::tensor::composition;
use idenso::tensor::{ContractSettings, SymbolicTensor};
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::ParseState;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::representation::{BaseRepName, LibraryRep, Minkowski};
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::utils::resources;
use symbolica::utils::resources::ResourceBudget;

impl TensorReducer {
    pub(super) fn has_tensor_sum(expression: AtomView<'_>) -> Result<bool, TensorReductionError> {
        let mut found = false;
        let mut resource_error = None;
        expression.visitor(&mut |value| {
            if let Err(error) = resources::checkpoint() {
                resource_error = Some(error);
                return false;
            }
            if matches!(value, AtomView::Add(_))
                && value.contains_symbol(Minkowski::selfless_symbol())
            {
                found = true;
            }
            !found
        });
        if let Some(error) = resource_error {
            return Err(error.into());
        }
        Ok(found)
    }

    pub(super) fn reduce_blocks(
        &self,
        source: &SymbolicTensor<PartialStructure>,
        external: Option<&ExternalProjector>,
    ) -> Result<TensorReduction, TensorReductionError> {
        resources::checkpoint()?;
        // Only loop/spectator dots need open ports. Invariant dots, including
        // inverse Gram factors, must remain scalar coefficients.
        let mut admission_error = None;
        let indexed = source.expression().replace_map(|value, _, out| {
            if admission_error.is_some() {
                return;
            }
            if let Err(error) = resources::checkpoint() {
                admission_error = Some(TensorReductionError::Resource(error));
                return;
            }
            let base = match value {
                AtomView::Pow(power) => power.get_base(),
                base => base,
            };
            if !matches!(base, AtomView::Fun(f) if f.get_symbol() == SPENSO_TAG.dot) {
                return;
            }
            match self.parse_monomial(value) {
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
        let indexed = source.with_checked_expression(indexed)?;
        resources::checkpoint()?;
        let rows = indexed
            .coefficient_list_with(|value| {
                match self.indexed_vector(value) {
                    Ok(Some(vector)) => !self.is_integrated(&vector),
                    Err(error) => {
                        admission_error = Some(error);
                        false
                    }
                    Ok(None) => {
                        // Metrics rotate covariantly with the loop block. Other
                        // Lorentz tensors (including spin matrices) stay fixed.
                        matches!(value, AtomView::Fun(f) if f.get_symbol() != ETS.metric && f.get_symbol() != SPENSO_TAG.dot)
                            && value.contains_symbol(Minkowski::selfless_symbol())
                    }
                }
            })
            .map_err(|e| TensorReductionError::TensorAlgebra(e.to_string()))?;
        resources::checkpoint()?;
        if let Some(error) = admission_error {
            return Err(error);
        }
        let state = ParseState::<AbstractIndex>::default();
        state.reserve_indices(indexed.expression().as_view());
        let mut terms = Vec::new();
        let mut max_rank = 0;
        for (outside, inside) in rows {
            resources::checkpoint()?;
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
            let indices = (0..rank).map(|_| state.fresh_index()).collect::<Vec<_>>();
            let fresh = slots
                .iter()
                .zip(&indices)
                .map(|(slot, index)| slot.rep().slot::<AbstractIndex, _>(*index).to_atom())
                .collect::<Vec<_>>();
            let representations = slots
                .iter()
                .map(|slot| slot.rep().to_symbolic([]))
                .collect::<Vec<_>>();
            let slots = slots.iter().map(|slot| slot.to_atom()).collect::<Vec<_>>();
            // Inner contractions own distinct dummies even when retained inside
            // a factored sum; they must not reconnect to the outside projector.
            let replacements = inside
                .structure()
                .logical_slots()
                .into_iter()
                .enumerate()
                .filter_map(|(position, slot)| {
                    slots
                        .iter()
                        .position(|expected| expected == &composition::port_atom(slot))
                        .map(|axis| (position, indices[axis]))
                })
                .collect();
            let renamed = inside.reindex_interface_ports(&replacements)?;
            let probe = SPENSO_TAG.rank_one_tensor_symbol("FeynKit::BlockProbe");
            // Probe vectors use the same native representation as the block
            // ports; the raw scalar dimension belongs only to projector weights.
            let mut probe_index = 0_usize;
            let probes = representations
                .iter()
                .map(|representation| {
                    loop {
                        resources::checkpoint()?;
                        let candidate = FunctionBuilder::new(probe)
                            .add_arg(probe_index)
                            .add_arg(representation)
                            .finish();
                        probe_index = probe_index
                            .checked_add(1)
                            .ok_or(TensorReductionError::CounterOverflow)?;
                        if !self.external_vectors.contains(&candidate) {
                            break Ok::<_, TensorReductionError>(candidate);
                        }
                    }
                })
                .collect::<Result<Vec<_>, _>>()?;
            let monomial = TensorMonomial {
                scalar: Atom::one(),
                integrated: probes.clone(),
                outside: slots
                    .iter()
                    .map(|slot| {
                        let AtomView::Fun(function) = slot.as_view() else {
                            unreachable!()
                        };
                        let index = function.iter().last().unwrap().to_owned();
                        Outside::Index {
                            slot: slot.clone(),
                            index,
                        }
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
            for term in projected {
                resources::checkpoint()?;
                // Apply the existing projector to a general covariant block:
                // every probe occurrence denotes one of its open tensor ports.
                let projector = term.tensor.expression().replace_map(|value, _, out| {
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
                // The probe recipe becomes a new tensor projector with these
                // explicit boundary slots. Existing inside/outside bodies stay typed.
                let projector = SymbolicTensor::infer(projector)?;
                let tensor = renamed
                    .multiply(&outside)?
                    .multiply(&projector)?
                    .contract(ContractSettings {
                        expand: false,
                        collect_chains: false,
                        collect_traces: false,
                        ..Default::default()
                    })?
                    .to_dots()?;
                let tensor = source.with_checked_expression(tensor.expression().clone())?;
                if !tensor.is_zero() {
                    terms.push(TensorReductionTerm {
                        coefficient: term.coefficient,
                        tensor,
                        integrated_orbit: None,
                        projector_orbit: None,
                        budget: ResourceBudget::capture_current(),
                    });
                }
                if let Some(limit) = self.output_term_limit
                    && terms.len() > limit
                {
                    return Err(TensorReductionError::OutputLimit {
                        terms: terms.len(),
                        limit,
                    });
                }
            }
        }
        let fully_contracted = source
            .structure()
            .logical_slots()
            .iter()
            .all(|slot| slot.rep().rep != LibraryRep::from(Minkowski {}));
        Ok(TensorReduction {
            terms,
            zero: source.with_checked_expression(Atom::Zero)?,
            max_rank,
            fully_contracted,
            budget: ResourceBudget::capture_current(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reduction::{dot, minkowski_slot};
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
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert_eq!(result.terms().len(), 1);
        assert!(result.is_fully_contracted());
        assert!(
            TensorReducer::has_tensor_sum(result.expression().unwrap().into_expression().as_view())
                .unwrap()
        );
        // Independently contract this small rank-two identity to scalar dots.
        let actual = result
            .expression()
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
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let expected = parse!(
            "spenso::dot(block_mixed::k(spenso::mink(D)),block_mixed::k(spenso::mink(D)))
            *spenso::g(spenso::mink(D,mu),spenso::mink(D,nu))/D
            +block_mixed::p(spenso::mink(D,mu))*block_mixed::p(spenso::mink(D,nu))"
        );
        assert!(
            (result.expression().unwrap().into_expression() - expected)
                .expand()
                .is_zero()
        );
        assert!(!result.is_fully_contracted());
        let _ = p;
    }
}
