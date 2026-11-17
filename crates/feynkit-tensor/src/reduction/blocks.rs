//! Project the open boundary of a loop-dependent tensor block. Its internal
//! contractions and vertex sums are coefficients, not projector monomials.

use super::external::ExternalProjector;
use super::{
    Outside, TensorMonomial, TensorReducer, TensorReduction, TensorReductionError,
    TensorReductionTerm, coset_key_from_parts, indexed_vector, matching_coset_key,
    perfect_matching_count, perfect_matchings,
};
use crate::OrthogonalWeingarten;
use idenso::shorthands::UndoShorthands;
use idenso::tensor::composition;
use idenso::tensor::{ContractSettings, SymbolicTensor};
use spenso::network::library::symbolic::ETS;
use spenso::network::parsing::ParseState;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{
    PartialIndex, PartialSlot, PartialStructure, PartialStructureExt,
};
use spenso::structure::representation::{BaseRepName, LibraryRep, Minkowski};
use spenso::structure::slot::IsAbstractSlot;
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::utils::resources;
use symbolica::utils::resources::ResourceBudget;

/// Row strategy of the block path.
///
/// Both non-legacy modes open every `dot(k, q)^n` with dummies keyed by
/// `(k, q, copy)`, so that equal tensor structures become equal collector
/// rows, and project each row once per pairing sigma of its ports.
/// [`RowMode::Relabel`] (the default) builds the pieces of a row with
/// explicit ports by relabelling them, without tensor inference or
/// multiplication; other rows use [`RowMode::Grouped`]. The environment
/// variable `FEYNKIT_TENSORREDUCE_ROWS` selects `grouped` or `legacy` for
/// comparison. `legacy` is the earlier path: a fresh dummy per opened dot and
/// one contraction per labelled pair of pairings.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum RowMode {
    Legacy,
    Grouped,
    Relabel,
}

impl RowMode {
    fn from_env() -> Self {
        static MODE: OnceLock<RowMode> = OnceLock::new();
        *MODE.get_or_init(
            || match std::env::var("FEYNKIT_TENSORREDUCE_ROWS").as_deref() {
                Ok("legacy") => RowMode::Legacy,
                Ok("grouped") => RowMode::Grouped,
                _ => RowMode::Relabel,
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
            RowMode::Grouped | RowMode::Relabel => self.open_block_dots(source, true)?,
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
    /// structure is, so equal rows merge. If the same dot occurs at two nesting
    /// levels of one product its label would occur more than twice once the
    /// collector distributes the sum; this returns `None` then (and when the
    /// source check fails), and the caller opens with fresh dummies instead. Without `deterministic` every
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
        let mode = self.row_mode();
        // The Weingarten weights of one rank, for every pair of pairings.
        let mut tables: HashMap<usize, PairingTable> = HashMap::new();
        // Relabelling needs the source's explicit open labels; one shared pool
        // of pair labels serves every row, since rows are separate summands.
        let source_ports = explicit_ports(source.structure());
        let mut pair_labels: Vec<AbstractIndex> = Vec::new();
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
            if mode == RowMode::Relabel && external.is_none() {
                if rank % 2 == 1 {
                    continue;
                }
                if let Some(source_ports) = &source_ports
                    && let Some(plan) = RelabelPlan::new(source_ports, &outside, &inside)
                {
                    if !tables.contains_key(&rank) {
                        tables.insert(rank, self.pairing_table(rank)?);
                    }
                    while pair_labels.len() < rank / 2 {
                        pair_labels.push(state.fresh_index());
                    }
                    if self.project_relabelled_row(
                        source,
                        &outside,
                        &inside,
                        &plan,
                        &tables[&rank],
                        &pair_labels,
                        &mut terms,
                    )? {
                        continue;
                    }
                }
            }
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
            if mode == RowMode::Legacy || external.is_some() {
                self.project_row_pairs(
                    source,
                    &outside,
                    &renamed,
                    &slots,
                    &fresh,
                    &representations,
                    external,
                    &mut terms,
                )?;
                continue;
            }
            if rank % 2 == 1 {
                continue;
            }
            if !tables.contains_key(&rank) {
                tables.insert(rank, self.pairing_table(rank)?);
            }
            let table = &tables[&rank];
            // The projector term (sigma, tau) is W(sigma, tau) times metrics on
            // the block ports (sigma) and metrics on the outside slots (tau).
            // By bilinearity, sum_(sigma,tau) W (I.sigma)(O.tau) equals
            // sum_sigma (I.sigma)(O.sum_tau W tau): contract the large block
            // once per sigma instead of once per (sigma, tau).
            for (sigma_index, sigma) in table.pairings.iter().enumerate() {
                resources::checkpoint()?;
                let mut weights = Atom::Zero;
                for (tau_index, tau) in table.pairings.iter().enumerate() {
                    let weight = &table.weights[sigma_index][tau_index];
                    if !weight.is_zero() {
                        weights += weight * metric_pairing(&slots, tau);
                    }
                }
                if weights.is_zero() {
                    continue;
                }
                let block = renamed
                    .multiply(&SymbolicTensor::infer(metric_pairing(&fresh, sigma))?)?
                    .contract(block_contract_settings())?
                    .to_dots()?;
                let weighted = outside
                    .multiply(&SymbolicTensor::infer(weights)?)?
                    .contract(block_contract_settings())?;
                let tensor = block
                    .multiply(&weighted)?
                    .contract(block_contract_settings())?
                    .to_dots()?;
                let tensor = source.with_checked_expression(tensor.expression().clone())?;
                self.push_block_term(&mut terms, Atom::one(), tensor)?;
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

    fn push_block_term(
        &self,
        terms: &mut Vec<TensorReductionTerm>,
        coefficient: Atom,
        tensor: SymbolicTensor<PartialStructure>,
    ) -> Result<(), TensorReductionError> {
        if !tensor.is_zero() {
            terms.push(TensorReductionTerm {
                coefficient,
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
        Ok(())
    }

    /// Project one row by relabelling its ports; `false` if it cannot.
    ///
    /// The vacuum average of the row is
    /// `sum_(sigma,tau) W(sigma,tau) (I g_sigma)(O g_tau)`. Contracting a
    /// metric is a relabelling, so `B_sigma = I g_sigma` is `I` with both
    /// ports of every pair of sigma joined, and `O_tau = O g_tau` is `O` with,
    /// for every pair of tau, the leg on the outside renamed to the other port
    /// (two outside legs then contract; an outside leg paired with an open
    /// source slot takes over its label), or times the metric of two open
    /// source slots. Each piece has a known interface (`B_sigma` is a scalar,
    /// `O_tau` is `O` with these changes), so it is built without inferring a
    /// projector or multiplying tensors.
    ///
    /// If every summand of `I` holds each port on a direct vector factor,
    /// `B_sigma` is the sum of the other factors times
    /// `prod_(a,b in sigma) dot(v_a, v_b)`. Otherwise the ports are relabelled
    /// to one fresh label per pair and `B_sigma` is contracted. The row adds
    /// `B_sigma sum_tau W(sigma,tau) O_tau` for every sigma, each checked once
    /// against the source interface. Nothing is added when this returns
    /// `false` (a port spelling that does not occur in its piece).
    #[allow(clippy::too_many_arguments)]
    fn project_relabelled_row(
        &self,
        source: &SymbolicTensor<PartialStructure>,
        outside: &SymbolicTensor<PartialStructure>,
        inside: &SymbolicTensor<PartialStructure>,
        plan: &RelabelPlan,
        table: &PairingTable,
        pair_labels: &[AbstractIndex],
        terms: &mut Vec<TensorReductionTerm>,
    ) -> Result<bool, TensorReductionError> {
        let mut outsides = Vec::with_capacity(table.pairings.len());
        for tau in &table.pairings {
            resources::checkpoint()?;
            let Some(piece) = plan.outside(outside, tau) else {
                return Ok(false);
            };
            outsides.push(
                piece
                    .contract(block_contract_settings())?
                    .to_dots()?
                    .into_expression(),
            );
        }
        let monomials = self.port_monomials(plan, inside.expression().as_view())?;
        let mut row = Vec::with_capacity(table.pairings.len());
        for (sigma_index, sigma) in table.pairings.iter().enumerate() {
            resources::checkpoint()?;
            let mut weighted = Atom::Zero;
            for (tau_index, piece) in outsides.iter().enumerate() {
                let weight = &table.weights[sigma_index][tau_index];
                if !weight.is_zero() && !piece.is_zero() {
                    weighted += weight * piece;
                }
            }
            if weighted.is_zero() {
                continue;
            }
            let block = match &monomials {
                Some(monomials) => paired_monomials(monomials, sigma)?,
                None => {
                    let Some(block) = plan.inside(inside, sigma, pair_labels) else {
                        return Ok(false);
                    };
                    block
                        .contract(block_contract_settings())?
                        .to_dots()?
                        .into_expression()
                }
            };
            row.push(block * weighted);
        }
        for term in row {
            resources::checkpoint()?;
            let tensor = source.with_checked_expression(term)?;
            self.push_block_term(terms, Atom::one(), tensor)?;
        }
        Ok(true)
    }

    /// Split every summand of a row's inside into its port vectors (compact,
    /// by port) and the product of its other factors. `None` unless each port
    /// occurs exactly once per summand, as a direct vector factor, and no
    /// other factor holds a port (a metric, a sum or a function argument).
    fn port_monomials(
        &self,
        plan: &RelabelPlan,
        inside: AtomView<'_>,
    ) -> Result<Option<Vec<PortMonomial>>, TensorReductionError> {
        let summands = match inside {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            term => vec![term],
        };
        let holds_port = |factor: AtomView<'_>| {
            let mut found = false;
            factor.visitor(&mut |node| {
                if let AtomView::Fun(function) = node
                    && function.get_nargs() == 2
                    && function.get_symbol() == Minkowski::selfless_symbol()
                    && plan.port_atoms.iter().any(|port| port.as_view() == node)
                {
                    found = true;
                }
                !found
            });
            found
        };
        let mut monomials = Vec::with_capacity(summands.len());
        for term in summands {
            resources::checkpoint()?;
            let factors = match term {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                factor => vec![factor],
            };
            let mut vectors = vec![None; plan.port_atoms.len()];
            let mut rest = Vec::with_capacity(factors.len());
            for factor in factors {
                let port = match self.indexed_vector(factor) {
                    Ok(Some(vector)) => plan
                        .port_atoms
                        .iter()
                        .position(|port| *port == vector.slot)
                        .map(|port| (port, vector.compact)),
                    _ => None,
                };
                match port {
                    Some((port, compact)) if vectors[port].is_none() => {
                        vectors[port] = Some(compact)
                    }
                    Some(_) => return Ok(None),
                    None if holds_port(factor) => return Ok(None),
                    None => rest.push(factor),
                }
            }
            let Some(vectors) = vectors.into_iter().collect::<Option<Vec<_>>>() else {
                return Ok(None);
            };
            let rest = rest
                .into_iter()
                .fold(Atom::one(), |product, factor| product * factor.to_owned());
            monomials.push(PortMonomial { rest, vectors });
        }
        Ok(Some(monomials))
    }

    /// Weingarten weights of every (sigma, tau) pair of pairings of `rank`
    /// slots, with the reducer's pairing diagnostics.
    fn pairing_table(&self, rank: usize) -> Result<PairingTable, TensorReductionError> {
        resources::checkpoint()?;
        let pairing_count = perfect_matching_count(rank)?;
        if let Some(limit) = self.pairing_limit
            && pairing_count > limit
        {
            return Err(TensorReductionError::PairingLimit {
                rank,
                pairings: pairing_count,
                limit,
            });
        }
        if let Some(limit) = self.pairing_product_limit {
            let products = pairing_count
                .checked_mul(pairing_count)
                .ok_or(TensorReductionError::CounterOverflow)?;
            if products > limit {
                return Err(TensorReductionError::PairingProductLimit {
                    rank,
                    products,
                    limit,
                });
            }
        }
        let pairings = perfect_matchings(rank)?;
        let table = OrthogonalWeingarten::new(rank / 2, self.dimension.clone())?;
        let by_key = table
            .coefficients()
            .map(|(coset, coefficient)| {
                Ok((coset_key_from_parts(coset.parts())?, coefficient.clone()))
            })
            .collect::<Result<HashMap<_, _>, TensorReductionError>>()?;
        let weights = pairings
            .iter()
            .map(|sigma| {
                pairings
                    .iter()
                    .map(|tau| Ok(by_key[&matching_coset_key(sigma, tau)?].clone()))
                    .collect::<Result<Vec<_>, TensorReductionError>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(PairingTable { pairings, weights })
    }

    /// The previous row projection: one probe monomial reduced per row and
    /// one contraction of the whole row per labelled (sigma, tau) term. Used
    /// with an external projector and for `FEYNKIT_TENSORREDUCE_ROWS=legacy`.
    #[allow(clippy::too_many_arguments)]
    fn project_row_pairs(
        &self,
        source: &SymbolicTensor<PartialStructure>,
        outside: &SymbolicTensor<PartialStructure>,
        renamed: &SymbolicTensor<PartialStructure>,
        slots: &[Atom],
        fresh: &[Atom],
        representations: &[Atom],
        external: Option<&ExternalProjector>,
        terms: &mut Vec<TensorReductionTerm>,
    ) -> Result<(), TensorReductionError> {
        let rank = slots.len();
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
                        indexed_vector(&other.to_owned(), &slot.iter().last().unwrap().to_owned())
                    }
                    _ => return,
                };
            });
            // The probe recipe becomes a new tensor projector with these
            // explicit boundary slots. Existing inside/outside bodies stay typed.
            let projector = SymbolicTensor::infer(projector)?;
            let tensor = renamed
                .multiply(outside)?
                .multiply(&projector)?
                .contract(block_contract_settings())?
                .to_dots()?;
            let tensor = source.with_checked_expression(tensor.expression().clone())?;
            self.push_block_term(terms, term.coefficient, tensor)?;
        }
        Ok(())
    }
}

/// The Weingarten weights of one rank: `weights[sigma][tau]` for the
/// perfect matchings `pairings` (each as its partner array).
struct PairingTable {
    pairings: Vec<Vec<usize>>,
    weights: Vec<Vec<Atom>>,
}

/// The slot atoms of an interface whose slots are all explicitly labelled.
fn explicit_ports(structure: &PartialStructure) -> Option<HashSet<Atom>> {
    structure
        .logical_slots()
        .into_iter()
        .map(|slot| {
            matches!(slot.aind, PartialIndex::Explicit(_)).then(|| composition::port_atom(slot))
        })
        .collect()
}

/// The ports of a row whose projection can be built by relabelling.
///
/// It exists only if every slot of both row pieces is explicitly labelled,
/// every slot of the inside is a Minkowski port of one representation (so
/// that a projected inside is a scalar and its port labels are the only
/// labels to rename), and every port either connects to the outside or is an
/// open slot of the source, never both. Other rows take the grouped path.
pub(super) struct RelabelPlan {
    ports: Vec<PartialSlot>,
    port_atoms: Vec<Atom>,
    connects: Vec<bool>,
}

impl RelabelPlan {
    pub(super) fn new(
        source_ports: &HashSet<Atom>,
        outside: &SymbolicTensor<PartialStructure>,
        inside: &SymbolicTensor<PartialStructure>,
    ) -> Option<Self> {
        let ports = inside.structure().logical_slots();
        let outside_ports = explicit_ports(outside.structure())?;
        let representation = ports.first().map(|slot| slot.rep());
        if ports.iter().any(|slot| {
            !matches!(slot.aind, PartialIndex::Explicit(_))
                || slot.rep().rep != LibraryRep::from(Minkowski {})
                || Some(slot.rep()) != representation
        }) {
            return None;
        }
        let port_atoms = ports
            .iter()
            .map(|slot| composition::port_atom(*slot))
            .collect::<Vec<_>>();
        let connects = port_atoms
            .iter()
            .map(|port| {
                let connects = outside_ports.contains(port);
                (connects != source_ports.contains(port)).then_some(connects)
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Self {
            ports,
            port_atoms,
            connects,
        })
    }

    /// `B_sigma` before contraction: the inside with both ports of every pair
    /// of `sigma` renamed to one label of `pair_labels`. Its interface is
    /// empty. `None` if a port spelling does not occur in the inside.
    pub(super) fn inside(
        &self,
        inside: &SymbolicTensor<PartialStructure>,
        sigma: &[usize],
        pair_labels: &[AbstractIndex],
    ) -> Option<SymbolicTensor<PartialStructure>> {
        let mut renames = Vec::with_capacity(sigma.len());
        for (pair, left) in (0..sigma.len())
            .filter(|left| *left < sigma[*left])
            .enumerate()
        {
            let label = self.ports[left]
                .rep()
                .slot::<AbstractIndex, _>(pair_labels[pair])
                .to_atom();
            renames.push((self.port_atoms[left].clone(), label.clone()));
            renames.push((self.port_atoms[sigma[left]].clone(), label));
        }
        Some(SymbolicTensor::new(
            relabel(inside.expression().as_view(), &renames)?,
            PartialStructure::from_logical_slots([]),
        ))
    }

    /// `O_tau` before contraction: the outside times the metrics of `tau`, as
    /// a relabelling. Two outside legs contract (one is renamed to the other);
    /// an outside leg paired with an open source slot takes over that label;
    /// two open source slots get their metric. `None` if a port spelling does
    /// not occur in the outside.
    pub(super) fn outside(
        &self,
        outside: &SymbolicTensor<PartialStructure>,
        tau: &[usize],
    ) -> Option<SymbolicTensor<PartialStructure>> {
        let mut slots = outside.structure().logical_slots();
        let position = |slots: &[PartialSlot], port: &Atom| {
            slots
                .iter()
                .position(|slot| composition::port_atom(*slot) == *port)
        };
        let mut renames = Vec::new();
        let mut metrics = Atom::one();
        for left in (0..tau.len()).filter(|left| *left < tau[*left]) {
            let right = tau[left];
            match (self.connects[left], self.connects[right]) {
                (true, true) => {
                    renames.push((
                        self.port_atoms[right].clone(),
                        self.port_atoms[left].clone(),
                    ));
                    slots.remove(position(&slots, &self.port_atoms[left])?);
                    slots.remove(position(&slots, &self.port_atoms[right])?);
                }
                (true, false) | (false, true) => {
                    let (leg, open) = if self.connects[left] {
                        (left, right)
                    } else {
                        (right, left)
                    };
                    renames.push((self.port_atoms[leg].clone(), self.port_atoms[open].clone()));
                    let at = position(&slots, &self.port_atoms[leg])?;
                    slots[at] = self.ports[open];
                }
                (false, false) => {
                    metrics = metrics
                        * FunctionBuilder::new(ETS.metric)
                            .add_arg(&self.port_atoms[left])
                            .add_arg(&self.port_atoms[right])
                            .finish();
                    slots.push(self.ports[left]);
                    slots.push(self.ports[right]);
                }
            }
        }
        Some(SymbolicTensor::new(
            relabel(outside.expression().as_view(), &renames)? * metrics,
            PartialStructure::from_logical_slots(slots),
        ))
    }
}

/// One summand of a row's inside: its port vectors (compact, by port) and
/// the product of its other factors.
struct PortMonomial {
    rest: Atom,
    vectors: Vec<Atom>,
}

/// `B_sigma` of monomial-shaped insides:
/// `sum_i rest_i prod_(a < sigma[a]) dot(v_ia, v_i,sigma[a])`.
fn paired_monomials(
    monomials: &[PortMonomial],
    sigma: &[usize],
) -> Result<Atom, TensorReductionError> {
    let mut block = Atom::Zero;
    for monomial in monomials {
        resources::checkpoint()?;
        let mut term = monomial.rest.clone();
        for left in (0..sigma.len()).filter(|left| *left < sigma[*left]) {
            term = term
                * FunctionBuilder::new(SPENSO_TAG.dot)
                    .add_arg(&monomial.vectors[left])
                    .add_arg(&monomial.vectors[sigma[left]])
                    .finish();
        }
        block += term;
    }
    Ok(block)
}

/// Rename explicit Minkowski slot atoms (`from -> to`, all at once). `None`
/// if some `from` does not occur.
fn relabel(expression: AtomView<'_>, renames: &[(Atom, Atom)]) -> Option<Atom> {
    if renames.is_empty() {
        return Some(expression.to_owned());
    }
    let mut applied = vec![false; renames.len()];
    let relabelled = expression.replace_map(|value, _, out| {
        if let AtomView::Fun(function) = value
            && function.get_nargs() == 2
            && function.get_symbol() == Minkowski::selfless_symbol()
            && let Some(rename) = renames.iter().position(|(from, _)| from.as_view() == value)
        {
            applied[rename] = true;
            **out = renames[rename].1.clone();
        }
    });
    applied.iter().all(|applied| *applied).then_some(relabelled)
}

/// `prod_(a < pairing[a]) g(slots[a], slots[pairing[a]])`.
fn metric_pairing(slots: &[Atom], pairing: &[usize]) -> Atom {
    (0..pairing.len())
        .filter(|left| *left < pairing[*left])
        .fold(Atom::one(), |product, left| {
            product
                * FunctionBuilder::new(ETS.metric)
                    .add_arg(&slots[left])
                    .add_arg(&slots[pairing[left]])
                    .finish()
        })
}

fn block_contract_settings() -> ContractSettings<'static> {
    ContractSettings {
        expand: false,
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
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
        assert_eq!(legacy.terms().len(), 3);
        let expected = parse!(
            "(block_rows::x1+block_rows::x2+block_rows::x3)
            *spenso::dot(block_rows::k(spenso::mink(D)),block_rows::k(spenso::mink(D)))
            *block_rows::q(spenso::mink(D,mu))/D"
        )
        .expand();
        assert_eq!(contracted(legacy), expected);
        for mode in [RowMode::Grouped, RowMode::Relabel] {
            let current = reducer
                .clone()
                .with_row_mode(mode)
                .reduce_blocks(&source, None)
                .unwrap();
            assert_eq!(current.terms().len(), 1, "{mode:?}");
            assert_eq!(contracted(current), expected, "{mode:?}");
        }
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
        for mode in [RowMode::Legacy, RowMode::Grouped, RowMode::Relabel] {
            let reduction = reducer
                .clone()
                .with_row_mode(mode)
                .reduce_blocks(&source, None)
                .unwrap();
            assert_eq!(contracted(reduction), expected, "{mode:?}");
        }
    }

    #[test]
    fn per_sigma_grouping_equals_the_per_pair_projection() {
        let _ = (
            spenso::vector_symbol!("block_sigma::k"),
            spenso::vector_symbol!("block_sigma::q"),
            spenso::vector_symbol!("block_sigma::p"),
        );
        // Rank-four block: ports D(k,q) and rho connect to the outside
        // (q and p), mu and nu are open. The scalar sum keeps it a block.
        let input = parse!(
            "spenso::dot(block_sigma::k(spenso::mink(D)),block_sigma::q(spenso::mink(D)))
            *block_sigma::k(spenso::mink(D,mu))*block_sigma::k(spenso::mink(D,nu))
            *block_sigma::k(spenso::mink(D,rho))*block_sigma::p(spenso::mink(D,rho))
            *(block_sigma::x+block_sigma::y*spenso::dot(block_sigma::k(spenso::mink(D)),block_sigma::k(spenso::mink(D))))"
        );
        let source = SymbolicTensor::infer(input).unwrap();
        let reducer =
            TensorReducer::new(parse!("D")).with_integrated_head(symbol!("block_sigma::k"));
        let pairs = reducer
            .clone()
            .with_row_mode(RowMode::Legacy)
            .reduce_blocks(&source, None)
            .unwrap();
        assert_eq!(pairs.max_rank(), 4);
        assert_eq!(pairs.terms().len(), 9, "3 x 3 labelled (sigma, tau) terms");
        let pairs = contracted(pairs);
        assert!(!pairs.is_zero());
        for mode in [RowMode::Grouped, RowMode::Relabel] {
            let grouped = reducer
                .clone()
                .with_row_mode(mode)
                .reduce_blocks(&source, None)
                .unwrap();
            assert_eq!(grouped.terms().len(), 3, "one term per sigma ({mode:?})");
            assert_eq!(pairs, contracted(grouped), "{mode:?}");
        }
    }

    /// Small block inputs of inside rank 0, 2 or 4 whose ports mix open
    /// source slots, dummies to an outside vector and opened dots.
    fn random_block_inputs() -> Vec<(usize, Atom)> {
        let _ = (
            spenso::vector_symbol!("block_relabel::k"),
            spenso::vector_symbol!("block_relabel::q"),
            spenso::vector_symbol!("block_relabel::p"),
        );
        let vector = |head: &str, index: Option<&str>| {
            let slot = match index {
                Some(index) => format!("spenso::mink(D,block_relabel::{index})"),
                None => "spenso::mink(D)".to_owned(),
            };
            parse!(format!("block_relabel::{head}({slot})"))
        };
        let dot = |left: &str, right: &str| {
            FunctionBuilder::new(SPENSO_TAG.dot)
                .add_arg(vector(left, None))
                .add_arg(vector(right, None))
                .finish()
        };
        // A scalar sum keeps every input on the block path.
        let block = parse!(
            "block_relabel::x+block_relabel::y*spenso::dot(block_relabel::k(spenso::mink(D)),block_relabel::k(spenso::mink(D)))"
        );
        let mut seed = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = |bound: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % bound
        };
        let mut inputs = Vec::new();
        for case in 0..30 {
            let rank = [0, 2, 4][case % 3];
            let mut atom = block.clone();
            let mut placed = 0;
            let mut open = 0;
            while placed < rank {
                match next(7) {
                    // A metric between two open source slots (contraction path).
                    4 if rank - placed >= 2 => {
                        let (left, right) = (format!("e{}", open + 1), format!("e{}", open + 2));
                        open += 2;
                        atom = atom
                            * parse!(format!(
                                "spenso::g(spenso::mink(D,block_relabel::{left}),spenso::mink(D,block_relabel::{right}))"
                            ));
                        placed += 2;
                    }
                    // A sum carrying an open port (contraction path).
                    5 => {
                        open += 1;
                        let leg = vector("k", Some(&format!("e{open}")));
                        atom =
                            atom * (parse!("block_relabel::z") * leg.clone() + dot("k", "k") * leg);
                        placed += 1;
                    }
                    // An open source slot carried by the loop vector.
                    0 | 4 => {
                        open += 1;
                        atom = atom * vector("k", Some(&format!("e{open}")));
                        placed += 1;
                    }
                    // A dummy to an outside spectator.
                    1 => {
                        let index = format!("d{placed}");
                        let spectator = ["q", "p"][next(2) as usize];
                        atom = atom * vector("k", Some(&index)) * vector(spectator, Some(&index));
                        placed += 1;
                    }
                    // An opened dot, possibly squared (two copies).
                    _ => {
                        let spectator = ["q", "p"][next(2) as usize];
                        let power = if rank - placed >= 2 && next(2) == 0 {
                            2
                        } else {
                            1
                        };
                        atom = atom * dot("k", spectator).pow(Atom::num(power));
                        placed += power as usize;
                    }
                }
            }
            // Some outside slots are not ports at all.
            if next(2) == 0 {
                open += 1;
                atom = atom * vector("p", Some(&format!("e{open}")));
            }
            if rank == 0 {
                atom = atom * vector("q", Some("e9"));
            }
            inputs.push((rank, atom));
        }
        inputs
    }

    #[test]
    fn relabelled_pieces_have_the_inferred_interfaces() {
        let inputs = random_block_inputs();
        let reducer =
            TensorReducer::new(parse!("D")).with_integrated_head(symbol!("block_relabel::k"));
        let mut planned = [0_usize; 3];
        // Rows on the contraction path and on the monomial (dot) path.
        let mut paths = [0_usize; 2];
        for (rank, input) in inputs {
            let source = SymbolicTensor::infer(input.clone()).unwrap();
            let indexed = reducer.open_block_dots(&source, true).unwrap().unwrap();
            let source_ports = explicit_ports(source.structure()).unwrap();
            let table = reducer.pairing_table(rank).unwrap();
            let state = ParseState::<AbstractIndex>::default();
            state.reserve_indices(indexed.expression().as_view());
            let labels = (0..rank / 2)
                .map(|_| state.fresh_index())
                .collect::<Vec<_>>();
            for (outside, inside) in reducer.block_rows(&indexed).unwrap() {
                assert_eq!(inside.structure().logical_slots().len(), rank, "{input}");
                let plan = RelabelPlan::new(&source_ports, &outside, &inside)
                    .unwrap_or_else(|| panic!("explicit row of {input}"));
                planned[rank / 2] += 1;
                let monomial = reducer
                    .port_monomials(&plan, inside.expression().as_view())
                    .unwrap()
                    .is_some();
                paths[usize::from(monomial)] += 1;
                for pairing in &table.pairings {
                    for piece in [
                        plan.inside(&inside, pairing, &labels).unwrap(),
                        plan.outside(&outside, pairing).unwrap(),
                    ] {
                        // The same slots; the logical order of an all-explicit
                        // interface follows the spelling and is not compared.
                        let inferred = SymbolicTensor::infer(piece.expression().clone()).unwrap();
                        assert_eq!(
                            piece.structure().canonical(),
                            inferred.structure().canonical(),
                            "{} in {input}",
                            piece.expression()
                        );
                    }
                }
            }
            let legacy = contracted(
                reducer
                    .clone()
                    .with_row_mode(RowMode::Legacy)
                    .reduce_blocks(&source, None)
                    .unwrap(),
            );
            for mode in [RowMode::Grouped, RowMode::Relabel] {
                let reduction = reducer
                    .clone()
                    .with_row_mode(mode)
                    .reduce_blocks(&source, None)
                    .unwrap();
                for term in reduction.terms() {
                    assert_eq!(term.tensor().structure(), source.structure());
                }
                assert_eq!(contracted(reduction), legacy, "{mode:?} for {input}");
            }
        }
        assert!(planned.iter().all(|rows| *rows > 0), "{planned:?}");
        assert!(paths.iter().all(|rows| *rows > 0), "{paths:?}");
    }

    #[test]
    fn rows_with_unlabelled_ports_fall_back_to_the_grouped_path() {
        let _ = (
            spenso::vector_symbol!("block_open::k"),
            spenso::vector_symbol!("block_open::q"),
        );
        let input = parse!(
            "spenso::dot(block_open::k(spenso::mink(D)),block_open::q(spenso::mink(D)))
            *block_open::k(spenso::mink(D,mu))
            *(block_open::x+block_open::y*spenso::dot(block_open::k(spenso::mink(D)),block_open::k(spenso::mink(D))))"
        );
        let source = SymbolicTensor::infer(input).unwrap();
        let reducer =
            TensorReducer::new(parse!("D")).with_integrated_head(symbol!("block_open::k"));
        let indexed = reducer.open_block_dots(&source, true).unwrap().unwrap();
        let source_ports = explicit_ports(source.structure()).unwrap();
        let rows = reducer.block_rows(&indexed).unwrap();
        assert_eq!(rows.len(), 1);
        let (outside, inside) = &rows[0];
        assert!(RelabelPlan::new(&source_ports, outside, inside).is_some());
        // The same row with one unlabelled (open) port, on either piece.
        let unlabel = |tensor: &SymbolicTensor<PartialStructure>| {
            let mut slots = tensor.structure().logical_slots();
            slots[0].set_aind(PartialIndex::Open(spenso::structure::partial::OpenPortId(
                0,
            )));
            SymbolicTensor::new(
                tensor.expression().clone(),
                PartialStructure::from_logical_slots(slots),
            )
        };
        assert!(RelabelPlan::new(&source_ports, outside, &unlabel(inside)).is_none());
        assert!(RelabelPlan::new(&source_ports, &unlabel(outside), inside).is_none());
        // A source with an open slot disables relabelling for every row.
        assert!(explicit_ports(unlabel(&source).structure()).is_none());
        // A port that is neither an outside leg nor a source slot.
        assert!(RelabelPlan::new(&HashSet::new(), outside, inside).is_none());

        // An inside with a non-Minkowski slot (here a Euclidean index to an
        // outside tensor) keeps the row off the relabel path; the grouped
        // path gives the same projection as the legacy one.
        let (tensor, vector) = (
            spenso::tensor_symbol!("block_open::T"),
            spenso::tensor_symbol!("block_open::V"),
        );
        let input = source.expression().clone()
            * FunctionBuilder::new(tensor)
                .add_arg(parse!("spenso::mink(D,block_open::nu)"))
                .add_arg(parse!("spenso::euc(2,block_open::i)"))
                .finish()
            * FunctionBuilder::new(vector)
                .add_arg(parse!("spenso::euc(2,block_open::i)"))
                .finish();
        let source = SymbolicTensor::infer(input).unwrap();
        let source_ports = explicit_ports(source.structure()).unwrap();
        let indexed = reducer.open_block_dots(&source, true).unwrap().unwrap();
        let rows = reducer.block_rows(&indexed).unwrap();
        assert_eq!(rows.len(), 1);
        assert!(RelabelPlan::new(&source_ports, &rows[0].0, &rows[0].1).is_none());
        let legacy = reducer
            .clone()
            .with_row_mode(RowMode::Legacy)
            .reduce_blocks(&source, None)
            .unwrap();
        let relabel = reducer
            .with_row_mode(RowMode::Relabel)
            .reduce_blocks(&source, None)
            .unwrap();
        assert_eq!(relabel.terms().len(), 1);
        let legacy = contracted(legacy);
        assert!(!legacy.is_zero());
        assert_eq!(legacy, contracted(relabel));
    }
}
