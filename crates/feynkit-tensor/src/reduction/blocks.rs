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
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::utils::resources;
use symbolica::utils::resources::ResourceBudget;

/// Row strategy of the block path.
///
/// The default opens every `dot(k, q)^n` with dummies keyed by `(k, q, copy)`,
/// so that equal tensor structures become equal collector rows.
/// `FEYNKIT_TENSORREDUCE_ROWS=legacy` restores the earlier path for
/// comparison: a fresh dummy per opened dot and one contraction per labelled
/// pair of pairings.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum RowMode {
    Legacy,
    Default,
}

impl RowMode {
    fn from_env() -> Self {
        static MODE: OnceLock<RowMode> = OnceLock::new();
        *MODE.get_or_init(
            || match std::env::var("FEYNKIT_TENSORREDUCE_ROWS").as_deref() {
                Ok("legacy") => RowMode::Legacy,
                _ => RowMode::Default,
            },
        )
    }
}

/// Whether a pooled dot label would occur more than twice in one product once
/// every nested sum is distributed.
///
/// Occurrences add up over the factors of a product (and inside function
/// arguments, conservatively), multiply under integer powers and take the
/// maximum over the summands of a sum.
fn reuses_pooled_label(
    expression: AtomView<'_>,
    pooled: &HashSet<Atom>,
) -> Result<bool, TensorReductionError> {
    type Counts = Vec<(usize, u64)>;
    fn add(into: &mut Counts, from: Counts, factor: u64) {
        for (label, count) in from {
            match into.iter_mut().find(|(seen, _)| *seen == label) {
                Some(entry) => entry.1 = entry.1.saturating_add(count.saturating_mul(factor)),
                None => into.push((label, count.saturating_mul(factor))),
            }
        }
    }
    fn visit<'a>(
        value: AtomView<'a>,
        labels: &HashMap<AtomView<'a>, usize>,
        clash: &mut bool,
    ) -> Result<Option<Counts>, TensorReductionError> {
        resources::checkpoint()?;
        let counts = match value {
            AtomView::Num(_) | AtomView::Var(_) => return Ok(None),
            AtomView::Fun(function) => {
                if function.get_symbol() == Minkowski::selfless_symbol()
                    && function.get_nargs() == 2
                {
                    let label = function.iter().last().unwrap();
                    return Ok(labels.get(&label).map(|id| vec![(*id, 1)]));
                }
                let mut counts = Counts::new();
                for argument in function.iter() {
                    if let Some(found) = visit(argument, labels, clash)? {
                        add(&mut counts, found, 1);
                    }
                }
                counts
            }
            AtomView::Mul(product) => {
                let mut counts = Counts::new();
                for factor in product.iter() {
                    if let Some(found) = visit(factor, labels, clash)? {
                        add(&mut counts, found, 1);
                    }
                }
                counts
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let Some(found) = visit(base, labels, clash)? else {
                    return Ok(None);
                };
                // A non-integer power of an indexed factor is not admitted;
                // count it as a clash so that the fresh-dummy path reports it.
                let factor = i64::try_from(exponent)
                    .ok()
                    .and_then(|exponent| u64::try_from(exponent).ok())
                    .unwrap_or(u64::MAX);
                let mut counts = Counts::new();
                add(&mut counts, found, factor);
                counts
            }
            AtomView::Add(sum) => {
                let mut counts = Counts::new();
                for term in sum.iter() {
                    for (label, count) in visit(term, labels, clash)?.unwrap_or_default() {
                        match counts.iter_mut().find(|(seen, _)| *seen == label) {
                            Some(entry) => entry.1 = entry.1.max(count),
                            None => counts.push((label, count)),
                        }
                    }
                }
                counts
            }
        };
        if counts.iter().any(|(_, count)| *count > 2) {
            *clash = true;
        }
        Ok((!counts.is_empty()).then_some(counts))
    }
    let labels = pooled
        .iter()
        .enumerate()
        .map(|(id, label)| (label.as_view(), id))
        .collect::<HashMap<_, _>>();
    let mut clash = false;
    visit(expression, &labels, &mut clash)?;
    Ok(clash)
}

/// Labelled block rows: the selected outside structure and its coefficient.
type BlockRows = Vec<(
    SymbolicTensor<PartialStructure>,
    SymbolicTensor<PartialStructure>,
)>;

impl TensorReducer {
    pub(super) fn row_mode(&self) -> RowMode {
        self.row_mode.unwrap_or_else(RowMode::from_env)
    }

    #[cfg(test)]
    pub(crate) fn with_row_mode(mut self, mode: RowMode) -> Self {
        self.row_mode = Some(mode);
        self
    }

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
        let mode = self.row_mode();
        let deterministic = match mode {
            RowMode::Legacy => None,
            RowMode::Default => self.open_block_dots(source, true)?,
        };
        let indexed = match deterministic {
            Some(indexed) => indexed,
            None => self
                .open_block_dots(source, false)?
                .expect("fresh dummies always admit the opened expression"),
        };
        let rows = self.block_rows(&indexed)?;
        self.project_block_rows(source, &indexed, rows, external)
    }

    /// Open every dot of an integrated vector with a spectator into explicit
    /// indexed vectors.
    ///
    /// With `deterministic`, `dot(k, q)^n` becomes
    /// `prod_j k(D(k,q,j)) q(D(k,q,j))` with one label per `(k, q, copy)` from
    /// a pool that avoids every label of `source`. The selected outside
    /// structures of the collector are then equal whenever their tensor
    /// structure is, so equal rows merge. The same label can occur twice in one
    /// product only if the same dot occurs at two nesting levels of it; such an
    /// expression fails the source check and this returns `None`, and the
    /// caller opens with fresh dummies instead. Without `deterministic` every
    /// occurrence gets a fresh dummy (`undo_dots`).
    pub(super) fn open_block_dots(
        &self,
        source: &SymbolicTensor<PartialStructure>,
        deterministic: bool,
    ) -> Result<Option<SymbolicTensor<PartialStructure>>, TensorReductionError> {
        resources::checkpoint()?;
        let labels = deterministic.then(|| {
            let state = ParseState::<AbstractIndex>::default();
            state.reserve_indices(source.expression().as_view());
            state
        });
        let mut pool: HashMap<(Atom, Atom, usize), Atom> = HashMap::new();
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
                    if let Some(labels) = &labels
                        && term.integrated.len() == term.outside.len()
                    {
                        let mut opened = term.scalar.clone();
                        let mut complete = true;
                        for (copy, (vector, outside)) in
                            term.integrated.iter().zip(&term.outside).enumerate()
                        {
                            let Outside::Vector(spectator) = outside else {
                                complete = false;
                                break;
                            };
                            let index = pool
                                .entry((vector.clone(), spectator.clone(), copy))
                                .or_insert_with(|| Atom::from(labels.fresh_index()));
                            opened = opened
                                * indexed_vector(vector, index)
                                * indexed_vector(spectator, index);
                        }
                        if complete {
                            **out = opened;
                            return;
                        }
                    }
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
        // Idenso scopes the dummies of a nested sum to its summands, so the
        // source check accepts `k(D) q(D) (y + x k(D) q(D))`. The collector
        // distributes that sum, after which `D` would occur four times.
        if !pool.is_empty() {
            let pooled = pool.into_values().collect::<HashSet<_>>();
            if reuses_pooled_label(indexed.as_view(), &pooled)? {
                return Ok(None);
            }
        }
        match source.with_checked_expression(indexed) {
            Ok(indexed) => Ok(Some(indexed)),
            Err(_) if deterministic => Ok(None),
            Err(error) => Err(error.into()),
        }
    }

    /// Split an opened block source into (selected outside, coefficient) rows.
    /// The collector merges rows whose selected structures are equal.
    pub(super) fn block_rows(
        &self,
        indexed: &SymbolicTensor<PartialStructure>,
    ) -> Result<BlockRows, TensorReductionError> {
        resources::checkpoint()?;
        let mut admission_error = None;
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
        Ok(rows)
    }

    fn project_block_rows(
        &self,
        source: &SymbolicTensor<PartialStructure>,
        indexed: &SymbolicTensor<PartialStructure>,
        rows: BlockRows,
        external: Option<&ExternalProjector>,
    ) -> Result<TensorReduction, TensorReductionError> {
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

    /// Fully contracted, expanded scalar form of a reduction.
    fn contracted(reduction: TensorReduction) -> Atom {
        reduction
            .into_expression()
            .unwrap()
            .contract(ContractSettings::default())
            .unwrap()
            .to_dots()
            .unwrap()
            .into_expression()
            .expand()
    }

    #[test]
    fn deterministic_dot_labels_merge_equal_outside_structures() {
        let _ = (
            spenso::vector_symbol!("block_rows::k"),
            spenso::vector_symbol!("block_rows::q"),
        );
        let input = parse!(
            "block_rows::x1*spenso::dot(block_rows::k(spenso::mink(D)),block_rows::q(spenso::mink(D)))*block_rows::k(spenso::mink(D,mu))
            +block_rows::x2*spenso::dot(block_rows::k(spenso::mink(D)),block_rows::q(spenso::mink(D)))*block_rows::k(spenso::mink(D,mu))
            +block_rows::x3*spenso::dot(block_rows::k(spenso::mink(D)),block_rows::q(spenso::mink(D)))*block_rows::k(spenso::mink(D,mu))"
        );
        let source = SymbolicTensor::infer(input).unwrap();
        let reducer =
            TensorReducer::new(parse!("D")).with_integrated_head(symbol!("block_rows::k"));
        let fresh = reducer.open_block_dots(&source, false).unwrap().unwrap();
        assert_eq!(reducer.block_rows(&fresh).unwrap().len(), 3);
        let keyed = reducer.open_block_dots(&source, true).unwrap().unwrap();
        let rows = reducer.block_rows(&keyed).unwrap();
        assert_eq!(rows.len(), 1, "one outside structure q(D(k,q,0))");

        let legacy = reducer
            .clone()
            .with_row_mode(RowMode::Legacy)
            .reduce_blocks(&source, None)
            .unwrap();
        let current = reducer
            .clone()
            .with_row_mode(RowMode::Default)
            .reduce_blocks(&source, None)
            .unwrap();
        assert_eq!(legacy.terms().len(), 3);
        assert_eq!(current.terms().len(), 1);
        let expected = parse!(
            "(block_rows::x1+block_rows::x2+block_rows::x3)
            *spenso::dot(block_rows::k(spenso::mink(D)),block_rows::k(spenso::mink(D)))
            *block_rows::q(spenso::mink(D,mu))/D"
        )
        .expand();
        assert_eq!(contracted(legacy), expected);
        assert_eq!(contracted(current), expected);
    }

    #[test]
    fn a_dot_at_two_nesting_levels_falls_back_to_fresh_dummies() {
        let _ = (
            spenso::vector_symbol!("block_clash::k"),
            spenso::vector_symbol!("block_clash::q"),
        );
        // One label for both dots would give k(D)^2 q(D)^2 after the collector
        // distributes the sum, i.e. dot(k,k) dot(q,q) q(mu) without the 1/D.
        let input = parse!(
            "spenso::dot(block_clash::k(spenso::mink(D)),block_clash::q(spenso::mink(D)))
            *(block_clash::y*block_clash::k(spenso::mink(D,mu))
              +block_clash::x*spenso::dot(block_clash::k(spenso::mink(D)),block_clash::q(spenso::mink(D)))*block_clash::q(spenso::mink(D,mu)))"
        );
        let source = SymbolicTensor::infer(input).unwrap();
        let reducer =
            TensorReducer::new(parse!("D")).with_integrated_head(symbol!("block_clash::k"));
        assert!(reducer.open_block_dots(&source, true).unwrap().is_none());
        let expected = parse!(
            "spenso::dot(block_clash::k(spenso::mink(D)),block_clash::k(spenso::mink(D)))/D
            *(block_clash::y*block_clash::q(spenso::mink(D,mu))
              +block_clash::x*spenso::dot(block_clash::q(spenso::mink(D)),block_clash::q(spenso::mink(D)))*block_clash::q(spenso::mink(D,mu)))"
        )
        .expand();
        for mode in [RowMode::Legacy, RowMode::Default] {
            let reduction = reducer
                .clone()
                .with_row_mode(mode)
                .reduce_blocks(&source, None)
                .unwrap();
            assert_eq!(contracted(reduction), expected, "{mode:?}");
        }
    }
}
