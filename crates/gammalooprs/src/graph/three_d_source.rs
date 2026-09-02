use std::collections::{BTreeMap, BTreeSet};

use ahash::AHashSet;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, Flow, Hedge, HedgePair},
    subgraph::{Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike},
};
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    concrete_index::ExpandedIndex,
    representation::{LibraryRep, Minkowski, RepName},
};
use symbolica::domains::rational::Rational;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    function,
    graph::Graph as SymbolicaGraph,
    id::Replacement,
};
use three_dimensional_reps::{
    EnergyEdgeIndexMap, LinearSurface, MomentumSignature, ParsedGraph, ThreeDGraphSource,
    graph_io::{
        GraphIoError, ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge,
        ParsedGraphInternalEdge, initial_state_cut_external_alias,
    },
    utils::{rank_i64, solve_rational_system},
    validate_parsed_graph,
};

use crate::{
    cff::surface::{GammaLoopLinearEnergyExpr, LinearEnergyExpr},
    graph::{Graph, LMBext, LoopMomentumBasis},
    momentum::SignOrZero,
    numerator::energy_degree::{EnergyPowerAssignmentPlan, EquivalentEnergyCandidates},
    utils::{GS, W_, symbols::UvMomentumProvenanceRole},
    uv::uv_graph::UVE,
};

impl ThreeDGraphSource for Graph {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        GraphThreeDSource::new(self, &[])?.to_three_d_parsed_graph()
    }

    fn energy_edge_index_map(&self, parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        GraphThreeDSource::new(self, &[])
            .ok()?
            .energy_edge_index_map(parsed)
    }
}

pub(crate) struct GraphThreeDSource<'a> {
    graph: &'a Graph,
    contract_edges: AHashSet<EdgeIndex>,
    initial_state_cut_edges: AHashSet<EdgeIndex>,
    inner_loop_count: usize,
    outer_loop_edges: Vec<EdgeIndex>,
    edge_loop_coordinates: BTreeMap<EdgeIndex, Vec<Rational>>,
    parent_loop_coordinates: Vec<Vec<Rational>>,
    uv_edges: AHashSet<EdgeIndex>,
    // UV boundary incidence is provenance too. It is optional for the current
    // vacuum MUV expansion, while future non-vacuum schemes can retain source
    // crown hedges explicitly without inferring their attachment from momentum
    // balance.
    exact_uv_boundary_hedges: Vec<Hedge>,
    // Proper UV subgraphs use their own loop/external split. In particular, a
    // parent loop crossing the UV crown is an external source coordinate, not
    // an inactive loop to be discarded by the child residue map.
    exact_coordinate_lmb: Option<LoopMomentumBasis>,
    exact_sub_lmb_frame: Option<ExactUvSubLmbFrame>,
    exact_denominators: Option<&'a [FourDDenominator]>,
    exact_signatures: Vec<MomentumSignature>,
    // Exact incidence is inherited from the source graph. An unexpanded exact
    // source reuses the production parsed graph literally; rewritten terms
    // retain the deterministic exact topology built below.
    exact_parsed: Option<ParsedGraph>,
    // Kept in parsed-denominator order for diagnostics which inspect the
    // factor ordering. Semantic consumers use `exact_occurrences`, because a
    // source-preserving parsed graph can interleave cut carriers.
    exact_local_to_original_occurrence: Vec<usize>,
    // The single namespace boundary shared by exact numerator mapping, OSE
    // replacement, surface projection, and physical-owner projection.
    exact_occurrences: Vec<ExactParsedOccurrence>,
    exact_energy_edge_index_map: Option<EnergyEdgeIndexMap>,
}

/// Incidence frame used by an exact proper-UV source in sub-LMB coordinates.
/// A physical shell source retains the complete source-graph crown. A Taylor
/// vacuum keeps the same owner incidence and sub-LMB provenance after the UV
/// operator has set every external denominator shift to zero.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ExactUvSubLmbFrame {
    RetainedPhysicalCrown,
    TaylorVacuum,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ExactParsedOccurrence {
    local_edge_id: usize,
    energy_edge_id: usize,
    original_occurrence: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct FourDDenominator {
    pub(crate) source_edge: EdgeIndex,
    pub(crate) momentum: Atom,
    pub(crate) mass_squared: Atom,
    pub(crate) full_expr: Atom,
}

type ExactSourceEdgeCoordinates = (EdgeIndex, Vec<Rational>, Vec<(EdgeIndex, SignOrZero)>);

#[derive(Clone, Debug, PartialEq, Eq)]
struct ExactSourceOwnerOccurrence {
    energy_edge_id: usize,
    raw_momentum: Atom,
    raw_to_parsed_sign: i64,
}

/// The exact source coordinates needed to evaluate a parent numerator under
/// one source-local CFF energy map. The source graph is temporary; this owned
/// value deliberately survives it without introducing production orientation
/// IDs. Source coordinates follow the same carrier-occurrence maps as the
/// reference evaluator. Repeated denominators retain their immutable source
/// owner, and a factor-local energy-power assignment may choose only among the
/// serial occurrences reconstructed for that owner.
#[derive(Clone, Debug)]
pub(crate) struct ExactSourceEnergyMapper {
    inactive_loop_count: usize,
    parent_loop_coordinates: Vec<Vec<Rational>>,
    parent_loop_edges: Vec<EdgeIndex>,
    edge_coordinates: Vec<ExactSourceEdgeCoordinates>,
    source_edge_occurrences: BTreeMap<EdgeIndex, Vec<ExactSourceOwnerOccurrence>>,
    cut_alias_edges: AHashSet<EdgeIndex>,
    exact_ose_replacements: Vec<Replacement>,
}

impl ExactSourceEnergyMapper {
    pub(crate) fn exact_ose_replacements(&self) -> &[Replacement] {
        &self.exact_ose_replacements
    }

    /// Expose only the serial exact occurrences reconstructed from each source
    /// edge. This is the certification boundary for owner-local dispatch.
    pub(crate) fn equivalent_energy_candidates(
        &self,
        physical_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Result<EquivalentEnergyCandidates> {
        let physical_edges = physical_edges.into_iter().collect::<BTreeSet<_>>();
        Ok(EquivalentEnergyCandidates::try_from_source_occurrences(
            self.source_edge_occurrences
                .iter()
                .filter(|(edge, _)| physical_edges.contains(edge))
                .map(|(edge, occurrences)| {
                    (
                        *edge,
                        occurrences
                            .iter()
                            .map(|occurrence| occurrence.energy_edge_id)
                            .collect(),
                    )
                }),
        )?)
    }

    /// Replace only temporal components. Spatial components remain in the
    /// production graph's parent momentum basis, including abstract indices.
    /// Inactive source energies are placeholders while the affine maps are
    /// assembled and are set to zero only after the complete factorized atom
    /// has been rewritten.
    #[cfg(test)]
    pub(crate) fn map_numerator(
        &self,
        loop_energy_map: &[LinearEnergyExpr],
        edge_energy_map: &[LinearEnergyExpr],
        numerator: &Atom,
    ) -> Result<Atom> {
        let assignments = self
            .source_edge_occurrences
            .iter()
            .filter_map(|(edge, occurrences)| {
                occurrences
                    .first()
                    .map(|occurrence| (*edge, occurrence.energy_edge_id))
            })
            .collect();
        let mapped =
            self.map_numerator_factor(loop_energy_map, edge_energy_map, numerator, &assignments)?;
        Ok(self.set_inactive_loop_energies_to_zero(mapped))
    }

    /// Reconstruct the physical energy of each requested owner directly from
    /// the exact source's loop lift and affine external shift. This differs
    /// deliberately from numerator sampling: every serial occurrence of a
    /// generalized lower-sector contact may sample zero even though its
    /// contour still has a nonzero pole. Empty occurrence assignments prevent
    /// those contact samples from overriding the physical source coordinates;
    /// they also make this reconstruction invariant under `D(Q)` versus
    /// `D(-Q)` denominator routing.
    #[cfg(test)]
    pub(crate) fn map_physical_owner_loop_lift_energies(
        &self,
        loop_energy_map: &[LinearEnergyExpr],
        edge_energy_map: &[LinearEnergyExpr],
        physical_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Result<Vec<(EdgeIndex, Atom)>> {
        let assignments = BTreeMap::new();
        physical_edges
            .into_iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .filter(|edge| !self.cut_alias_edges.contains(edge))
            .map(|edge| {
                if !self
                    .edge_coordinates
                    .iter()
                    .any(|(source_edge, _, _)| *source_edge == edge)
                {
                    return Err(eyre::eyre!(
                        "physical source edge {} has no exact source coordinates",
                        usize::from(edge),
                    ));
                }
                Ok((
                    edge,
                    self.set_inactive_loop_energies_to_zero(self.map_numerator_factor(
                        loop_energy_map,
                        edge_energy_map,
                        &GS.emr_mom(edge, GS.cind(0)),
                        &assignments,
                    )?),
                ))
            })
            .collect()
    }

    pub(crate) fn map_planned_numerator(
        &self,
        loop_energy_map: &[LinearEnergyExpr],
        edge_energy_map: &[LinearEnergyExpr],
        plan: &EnergyPowerAssignmentPlan,
    ) -> Result<Atom> {
        let mapped = plan.map_factors(|factor, assignments| {
            self.map_numerator_factor(loop_energy_map, edge_energy_map, factor, assignments)
        })?;
        Ok(self.set_inactive_loop_energies_to_zero(mapped))
    }

    /// Map one factor using exactly the occurrence assignment which was used
    /// to bound that factor during generalized CFF generation.
    fn map_numerator_factor(
        &self,
        loop_energy_map: &[LinearEnergyExpr],
        edge_energy_map: &[LinearEnergyExpr],
        numerator: &Atom,
        assignments: &BTreeMap<EdgeIndex, usize>,
    ) -> Result<Atom> {
        let parent_loop_count = self.parent_loop_coordinates.len();
        let expected_active_loop_count = parent_loop_count - self.inactive_loop_count;
        if loop_energy_map.len() != expected_active_loop_count {
            return Err(eyre::eyre!(
                "exact source supplies {} active loop-energy maps, expected {expected_active_loop_count}",
                loop_energy_map.len(),
            ));
        }

        let loop_temporal = |loop_index: usize| {
            FunctionBuilder::new(GS.loop_mom)
                .add_arg(loop_index)
                .add_arg(GS.cind(0))
                .finish()
        };
        let inactive_energies = (0..self.inactive_loop_count)
            .map(|index| loop_temporal(parent_loop_count + index))
            .collect::<Vec<_>>();
        if let Some(edge) = self
            .parent_loop_edges
            .iter()
            .find(|edge| self.cut_alias_edges.contains(edge))
        {
            return Err(eyre::eyre!(
                "parent loop carrier {} is also an initial-state cut alias",
                usize::from(*edge),
            ));
        }
        // These are energies of the exact source coordinates themselves. In
        // particular, a carrier denominator may be shifted by an external
        // momentum, so its pole energy is not the corresponding loop energy.
        let active_source_energies = loop_energy_map.iter().map(|energy| {
            energy
                .to_atom_gs(&[])
                .replace_multiple(&self.exact_ose_replacements)
        });
        let source_energies = inactive_energies
            .iter()
            .cloned()
            .chain(active_source_energies)
            .collect::<Vec<_>>();
        let source_energy = |coordinates: &[Rational]| {
            coordinates.iter().zip(&source_energies).fold(
                Atom::Zero,
                |sum, (coefficient, energy)| {
                    if *coefficient == 0 {
                        sum
                    } else {
                        sum + Atom::num(coefficient.clone()) * energy
                    }
                },
            )
        };
        let mut selected_occurrences = BTreeMap::new();
        for (edge, occurrence) in assignments {
            // Cut momenta are external aliases, so occurrence-local pole
            // energies are neither required nor relevant for their numerator.
            if self.cut_alias_edges.contains(edge) {
                continue;
            }
            let source_occurrence = self
                .source_edge_occurrences
                .get(edge)
                .and_then(|occurrences| {
                    occurrences
                        .iter()
                        .find(|candidate| candidate.energy_edge_id == *occurrence)
                })
                .ok_or_else(|| {
                    eyre::eyre!(
                        "exact occurrence {occurrence} is not a certified source-owner candidate for physical edge {}",
                        usize::from(*edge),
                    )
                })?;
            let energy = edge_energy_map
                .get(*occurrence)
                .ok_or_else(|| {
                    eyre::eyre!(
                        "exact source occurrence {occurrence} for physical edge {} is outside its {} edge-energy maps",
                        usize::from(*edge),
                        edge_energy_map.len(),
                    )
                })?
                .to_atom_gs(&[])
                .replace_multiple(&self.exact_ose_replacements);
            // Each factor uses the same exact occurrence which owns its bound,
            // including a zero sample in lower-sector contact terms.
            selected_occurrences.insert(*edge, (source_occurrence, energy));
        }
        let external_shift = |terms: &[(EdgeIndex, SignOrZero)]| {
            terms.iter().fold(Atom::Zero, |sum, (edge, sign)| {
                let energy = crate::utils::external_energy_atom_from_index(*edge);
                match sign {
                    SignOrZero::Zero => sum,
                    SignOrZero::Plus => sum + energy,
                    SignOrZero::Minus => sum - energy,
                }
            })
        };

        let minkowski_symbol = LibraryRep::from(Minkowski {}).symbol();
        let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
        let mut tag_mapping_error = None;
        let mut has_local_provenance = false;
        let tagged_numerator = numerator.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
                return;
            }
            let Some((owner, role, hard)) = GS.uv_momentum_provenance_data(momentum.get(0))
            else {
                return;
            };
            // A factorized component maps only its own coordinate domain. Keep
            // remote tags intact so the later component can consume the same
            // unexpanded numerator. Membership in `edge_coordinates`, rather
            // than the selected denominator occurrences, is authoritative: a
            // pinched local owner still belongs to this component and must not
            // leak a spatial provenance tag to the evaluator.
            if !self
                .edge_coordinates
                .iter()
                .any(|(edge, _, _)| *edge == owner)
            {
                return;
            }
            has_local_provenance = true;
            let index = momentum.get(1);
            let concrete_temporal = usize::try_from(index).is_ok_and(|index| index == 0)
                || index.to_owned() == AIND_SYMBOLS.cind.call(Atom::Zero);
            let abstract_index = matches!(
                index,
                AtomView::Fun(function) if function.get_symbol() == minkowski_symbol
            );
            // A fixed factor remains on its canonical owner occurrence, which
            // keeps it inside the same generalized Laurent/contact functional
            // as its CFF rank bound. Only a denominator-derived factor can be
            // reassigned among derivative-created serial copies.
            if !concrete_temporal && !abstract_index {
                **output = GS.erase_uv_momentum_provenance(&Atom::from(momentum.to_owned()));
                return;
            }
            let Some((occurrence, energy)) = selected_occurrences.get(&owner) else {
                tag_mapping_error.get_or_insert_with(|| {
                    eyre::eyre!(
                        "tagged temporal momentum owned by edge {} has no planned exact occurrence",
                        usize::from(owner),
                    )
                });
                return;
            };
            let hard_to_raw_sign = if hard == occurrence.raw_momentum {
                1
            } else if hard == -occurrence.raw_momentum.clone() {
                -1
            } else {
                tag_mapping_error.get_or_insert_with(|| {
                    eyre::eyre!(
                        "tagged hard momentum `{}` for source edge {} is neither the exact denominator momentum `{}` nor its negative",
                        hard,
                        usize::from(owner),
                        occurrence.raw_momentum,
                    )
                });
                return;
            };
            // The oriented occurrence map already carries the incidence sign
            // of a fixed owner. A dispatched derivative-grown factor must also
            // cross from its selected parsed copy back through the raw
            // rewritten-denominator frame before reaching the stored hard
            // frame.
            let occurrence_to_hard_sign = hard_to_raw_sign
                * if role == UvMomentumProvenanceRole::DenominatorDerived {
                    occurrence.raw_to_parsed_sign
                } else {
                    1
                };
            let signed_energy = if occurrence_to_hard_sign == -1 {
                -energy.clone()
            } else {
                energy.clone()
            };
            if concrete_temporal {
                **output = signed_energy;
                return;
            }
            let erased = GS.erase_uv_momentum_provenance(&Atom::from(momentum.to_owned()));
            let spatial = erased.replace_map(|hard_view, _, hard_output| {
                let AtomView::Fun(hard_momentum) = hard_view else {
                    return;
                };
                if hard_momentum.get_symbol() == GS.emr_mom
                    && hard_momentum.get_nargs() == 2
                    && hard_momentum.get(1) == index
                {
                    **hard_output = FunctionBuilder::new(GS.emr_vec)
                        .add_arg(hard_momentum.get(0))
                        .add_arg(hard_momentum.get(1))
                        .finish();
                }
            });
            **output = spatial + signed_energy * GS.energy_delta(index.to_owned());
        });
        if let Some(error) = tag_mapping_error {
            return Err(error);
        }

        // An untagged numerator belongs to the unchanged physical graph and
        // must reproduce its ordinary edge-energy map exactly. Every tagged
        // factor has already consumed the exact occurrence selected by the
        // same immutable plan that supplied its CFF energy bound. Fixed
        // factors stay on the canonical occurrence; only denominator-derived
        // factors may be balanced over serial copies.
        let mut exact_edge_energies = BTreeMap::new();
        if !has_local_provenance {
            for (edge, (occurrence, energy)) in &selected_occurrences {
                let source_momentum = FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(*edge))
                    .finish();
                let source_to_raw_sign = if occurrence.raw_momentum == source_momentum {
                    Some(1)
                } else if occurrence.raw_momentum == -source_momentum {
                    Some(-1)
                } else {
                    None
                };
                if let Some(source_to_raw_sign) = source_to_raw_sign {
                    exact_edge_energies.insert(
                        *edge,
                        if source_to_raw_sign * occurrence.raw_to_parsed_sign == -1 {
                            -energy.clone()
                        } else {
                            energy.clone()
                        },
                    );
                }
            }
        }
        let mut replacements = Vec::new();
        for (edge, coordinates, external_terms) in &self.edge_coordinates {
            // Initial-state cut momenta are literal signed external aliases.
            // They remain aliases even if an exact denominator also happens
            // to use the same physical momentum literally.
            let energy = if self.cut_alias_edges.contains(edge) {
                source_energy(coordinates) + external_shift(external_terms)
            } else {
                exact_edge_energies
                    .get(edge)
                    .cloned()
                    .unwrap_or_else(|| source_energy(coordinates) + external_shift(external_terms))
            };
            replacements.push(Replacement::new(
                GS.emr_mom(*edge, AIND_SYMBOLS.cind.call(Atom::Zero))
                    .to_pattern(),
                energy.clone().to_pattern(),
            ));
            replacements.push(Replacement::new(
                GS.emr_mom(*edge, &mink_index).to_pattern(),
                (GS.emr_vec_index(*edge, &mink_index) + energy * GS.energy_delta(&mink_index))
                    .to_pattern(),
            ));
        }
        for (loop_index, (coordinates, loop_edge)) in self
            .parent_loop_coordinates
            .iter()
            .zip(&self.parent_loop_edges)
            .enumerate()
        {
            let energy = exact_edge_energies
                .get(loop_edge)
                .cloned()
                .unwrap_or_else(|| source_energy(coordinates));
            let loop_index_atom = Atom::num(loop_index as i64);
            replacements.push(Replacement::new(
                loop_temporal(loop_index).to_pattern(),
                energy.clone().to_pattern(),
            ));
            for spatial_index in 1..=3 {
                replacements.push(Replacement::new(
                    FunctionBuilder::new(GS.loop_mom)
                        .add_arg(loop_index_atom.clone())
                        .add_arg(AIND_SYMBOLS.cind.call(spatial_index))
                        .finish()
                        .to_pattern(),
                    GS.emr_mom(*loop_edge, AIND_SYMBOLS.cind.call(spatial_index))
                        .to_pattern(),
                ));
            }
            replacements.push(Replacement::new(
                FunctionBuilder::new(GS.loop_mom)
                    .add_arg(loop_index_atom)
                    .add_arg(mink_index.as_view())
                    .finish()
                    .to_pattern(),
                (GS.emr_vec_index(*loop_edge, &mink_index) + energy * GS.energy_delta(&mink_index))
                    .to_pattern(),
            ));
        }

        Ok(tagged_numerator.replace_multiple(&replacements))
    }

    fn set_inactive_loop_energies_to_zero(&self, atom: Atom) -> Atom {
        let parent_loop_count = self.parent_loop_coordinates.len();
        (0..self.inactive_loop_count).fold(atom, |atom, index| {
            atom.replace(
                FunctionBuilder::new(GS.loop_mom)
                    .add_arg(parent_loop_count + index)
                    .add_arg(GS.cind(0))
                    .finish(),
            )
            .with(Atom::Zero)
        })
    }
}

impl<'a> GraphThreeDSource<'a> {
    pub(crate) fn new(
        graph: &'a Graph,
        contract_edges: &[EdgeIndex],
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        let contract_edges = contract_edges.iter().copied().collect::<AHashSet<_>>();
        let initial_state_cut_edges = graph
            .iter_edges_of(&graph.initial_state_cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect::<AHashSet<_>>();
        let contracts_edge = |edge_id: EdgeIndex| {
            contract_edges.contains(&edge_id) && !initial_state_cut_edges.contains(&edge_id)
        };

        let mut contracted_filter: SuBitGraph = graph.empty_subgraph();
        for (pair, edge_id, edge_data) in graph.underlying.iter_edges() {
            if pair.is_paired() && !edge_data.data.is_dummy && contracts_edge(edge_id) {
                contracted_filter.add(pair);
            }
        }
        let contracted_subgraph =
            InternalSubGraph::cleaned_filter_optimist(contracted_filter, &graph.underlying);
        let inner_lmb = graph
            .try_compatible_sub_lmb(
                &contracted_subgraph,
                graph.dummy_stripped_external_flows_of(&contracted_subgraph),
                &graph.loop_momentum_basis,
            )
            .map_err(|error| {
                GraphIoError::Source(format!(
                    "contracted edges do not admit a parent-compatible inner loop basis: {error}"
                ))
            })?;

        let parent_loop_count = graph.loop_momentum_basis.loop_edges.len();
        let mut basis_rows = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .filter(|edge_id| inner_lmb.loop_edges.contains(edge_id))
            .map(|edge_id| {
                graph.loop_momentum_basis.edge_signatures[*edge_id]
                    .internal
                    .iter()
                    .map(|sign| sign_to_i32(*sign))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let inner_loop_count = basis_rows.len();
        if inner_loop_count != inner_lmb.loop_edges.len()
            || exact_integer_rank(&basis_rows) != inner_loop_count
        {
            return Err(GraphIoError::Source(format!(
                "parent-compatible inner loop basis has {} carriers but rank {} in the {}-dimensional production basis",
                inner_lmb.loop_edges.len(),
                exact_integer_rank(&basis_rows),
                parent_loop_count,
            )));
        }

        let mut outer_loop_edges = Vec::new();
        let mut basis_rank = inner_loop_count;
        for (pair, edge_id, edge_data) in graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_id, _)| *edge_id)
        {
            if !pair.is_paired()
                || edge_data.data.is_dummy
                || initial_state_cut_edges.contains(&edge_id)
                || contracts_edge(edge_id)
            {
                continue;
            }
            let row = graph.loop_momentum_basis.edge_signatures[edge_id]
                .internal
                .iter()
                .map(|sign| sign_to_i32(*sign))
                .collect::<Vec<_>>();
            let mut trial = basis_rows.clone();
            trial.push(row.clone());
            let trial_rank = exact_integer_rank(&trial);
            if trial_rank > basis_rank {
                basis_rows.push(row);
                outer_loop_edges.push(edge_id);
                basis_rank = trial_rank;
            }
            if basis_rank == parent_loop_count {
                break;
            }
        }
        if basis_rank != parent_loop_count {
            return Err(GraphIoError::Source(format!(
                "contracted source supplies loop-signature rank {basis_rank}, expected the production rank {parent_loop_count}"
            )));
        }

        let basis_matrix = (0..parent_loop_count)
            .map(|column| {
                basis_rows
                    .iter()
                    .map(|row| Rational::from(row[column]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let solve_coordinates = |row: Vec<SignOrZero>| {
            if row.is_empty() {
                return Ok(Vec::new());
            }
            let rhs = row
                .iter()
                .map(|sign| Rational::from(sign_to_i32(*sign)))
                .collect::<Vec<_>>();
            solve_rational_system(basis_matrix.clone(), rhs).ok_or_else(|| {
                GraphIoError::Source(
                    "failed to solve exact contracted-source loop coordinates".to_string(),
                )
            })
        };
        let edge_loop_coordinates = graph
            .loop_momentum_basis
            .edge_signatures
            .iter()
            .map(|(edge_id, signature)| {
                Ok((
                    edge_id,
                    solve_coordinates(signature.internal.iter().copied().collect())?,
                ))
            })
            .collect::<three_dimensional_reps::graph_io::Result<BTreeMap<_, _>>>()?;
        for (pair, edge_id, edge_data) in graph.underlying.iter_edges() {
            if pair.is_paired()
                && !edge_data.data.is_dummy
                && !contracts_edge(edge_id)
                && edge_loop_coordinates[&edge_id]
                    .iter()
                    .take(inner_loop_count)
                    .any(|coordinate| *coordinate != 0)
            {
                return Err(GraphIoError::Source(format!(
                    "surviving edge {} retains an inner-loop coordinate after contraction",
                    usize::from(edge_id)
                )));
            }
        }
        let parent_loop_coordinates = (0..parent_loop_count)
            .map(|loop_index| {
                let row = (0..parent_loop_count)
                    .map(|column| {
                        if column == loop_index {
                            SignOrZero::Plus
                        } else {
                            SignOrZero::Zero
                        }
                    })
                    .collect::<Vec<_>>();
                solve_coordinates(row)
            })
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;

        Ok(Self {
            graph,
            contract_edges,
            initial_state_cut_edges,
            inner_loop_count,
            outer_loop_edges,
            edge_loop_coordinates,
            parent_loop_coordinates,
            uv_edges: AHashSet::new(),
            exact_uv_boundary_hedges: Vec::new(),
            exact_coordinate_lmb: None,
            exact_sub_lmb_frame: None,
            exact_denominators: None,
            exact_signatures: Vec::new(),
            exact_parsed: None,
            exact_local_to_original_occurrence: Vec::new(),
            exact_occurrences: Vec::new(),
            exact_energy_edge_index_map: None,
        })
    }

    #[cfg(test)]
    pub(crate) fn from_exact_denominators(
        graph: &'a Graph,
        denominators: &'a [FourDDenominator],
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        Self::from_exact_denominators_in_uv_edges(graph, denominators, [])
    }

    #[cfg(test)]
    pub(crate) fn from_exact_denominators_in_uv_edges(
        graph: &'a Graph,
        denominators: &'a [FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        Self::from_exact_denominators_in_uv_edges_and_boundaries(graph, denominators, uv_edges, [])
    }

    /// Build an exact source from original graph incidence. `uv_boundary_hedges`
    /// are source-graph crown hedges retained by a non-vacuum UV prescription;
    /// their endpoint is never reconstructed from momentum signatures.
    pub(crate) fn from_exact_denominators_in_uv_edges_and_boundaries(
        graph: &'a Graph,
        denominators: &'a [FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        Self::from_exact_denominators_in_coordinates(
            graph,
            denominators,
            uv_edges,
            uv_boundary_hedges,
            None,
            None,
        )
    }

    /// Build a proper UV source in the loop/external coordinates of its
    /// subgraph LMB. Parent loops crossing the retained crown are external in
    /// this frame and therefore survive the child CFF and Taylor operation.
    pub(crate) fn from_exact_denominators_in_uv_sub_lmb(
        graph: &'a Graph,
        denominators: &'a [FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        sub_lmb: &LoopMomentumBasis,
        frame: ExactUvSubLmbFrame,
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        Self::from_exact_denominators_in_coordinates(
            graph,
            denominators,
            uv_edges,
            uv_boundary_hedges,
            Some(sub_lmb),
            Some(frame),
        )
    }

    fn from_exact_denominators_in_coordinates(
        graph: &'a Graph,
        denominators: &'a [FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        coordinate_lmb: Option<&LoopMomentumBasis>,
        sub_lmb_frame: Option<ExactUvSubLmbFrame>,
    ) -> three_dimensional_reps::graph_io::Result<Self> {
        if coordinate_lmb.is_some() != sub_lmb_frame.is_some() {
            return Err(GraphIoError::Source(
                "an exact UV sub-LMB source requires an explicit incidence frame".to_string(),
            ));
        }
        let uv_edges = uv_edges.into_iter().collect::<AHashSet<_>>();
        if sub_lmb_frame.is_some()
            && let Some(owner) = denominators
                .iter()
                .map(|denominator| denominator.source_edge)
                .find(|owner| !uv_edges.contains(owner))
        {
            return Err(GraphIoError::Source(format!(
                "exact UV sub-LMB denominator owner {} is outside the UV source edges",
                usize::from(owner),
            )));
        }
        let exact_uv_boundary_hedges = uv_boundary_hedges
            .into_iter()
            .sorted()
            .dedup()
            .collect::<Vec<_>>();
        let mut uv_subgraph: SuBitGraph = graph.empty_subgraph();
        for edge_id in &uv_edges {
            uv_subgraph.add(graph[edge_id].1);
        }
        let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
        for hedge in &exact_uv_boundary_hedges {
            if !crown.includes(hedge) {
                return Err(GraphIoError::Source(format!(
                    "exact UV boundary hedge {hedge} is not in the source UV subgraph crown",
                )));
            }
        }
        match sub_lmb_frame {
            Some(ExactUvSubLmbFrame::RetainedPhysicalCrown)
                if exact_uv_boundary_hedges
                    != crown.included_iter().sorted().collect::<Vec<_>>() =>
            {
                return Err(GraphIoError::Source(
                    "an exact physical UV sub-LMB source must retain its complete non-dummy crown"
                        .to_string(),
                ));
            }
            Some(ExactUvSubLmbFrame::TaylorVacuum) if !exact_uv_boundary_hedges.is_empty() => {
                return Err(GraphIoError::Source(
                    "an exact Taylor-vacuum UV sub-LMB source cannot retain crown hedges"
                        .to_string(),
                ));
            }
            _ => {}
        }
        let production_signatures = denominators
            .iter()
            .map(|denominator| {
                denominator.momentum_signature(graph, uv_edges.contains(&denominator.source_edge))
            })
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
        let coordinate_signatures = coordinate_lmb
            .map(|lmb| {
                denominators
                    .iter()
                    .map(|denominator| {
                        denominator.momentum_signature_in_lmb(
                            lmb,
                            uv_edges.contains(&denominator.source_edge),
                        )
                    })
                    .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()
            })
            .transpose()?;
        if matches!(sub_lmb_frame, Some(ExactUvSubLmbFrame::TaylorVacuum))
            && let Some(signatures) = &coordinate_signatures
        {
            let shifted = denominators
                .iter()
                .zip(signatures)
                .filter(|(_, signature)| {
                    signature
                        .external_signature
                        .iter()
                        .any(|coefficient| *coefficient != 0)
                })
                .map(|(denominator, signature)| {
                    (
                        denominator.source_edge,
                        denominator.momentum.clone(),
                        signature.external_signature.clone(),
                    )
                })
                .collect::<Vec<_>>();
            if !shifted.is_empty() {
                return Err(GraphIoError::Source(format!(
                    "an exact Taylor-vacuum UV sub-LMB source must have zero external denominator shifts; shifted (owner, momentum, signature): {shifted:?}"
                )));
            }
        }
        // Exact 4D sources must use the same coordinate frame as the ordinary
        // production CFF. The configured graph LMB is not necessarily the
        // canonical carrier basis selected by `GraphThreeDSource::new`.
        let mut production_source = Self::new(graph, &[])?;
        if coordinate_lmb.is_none() && uv_edges.is_empty() && exact_uv_boundary_hedges.is_empty() {
            let production_parsed = production_source.to_three_d_parsed_graph()?;
            let production_energy_map = production_source
                .energy_edge_index_map(&production_parsed)
                .expect("GammaLoop production sources provide an edge-index map");
            let original_occurrences = denominators
                .iter()
                .enumerate()
                .filter(|(_, denominator)| denominator.is_original_graph_denominator(graph))
                .map(|(original, denominator)| (denominator.source_edge, original))
                .collect::<BTreeMap<_, _>>();
            let retained_source_edges = production_energy_map
                .internal
                .values()
                .copied()
                .map(EdgeIndex)
                .filter(|edge| !production_source.initial_state_cut_edges.contains(edge))
                .collect::<BTreeSet<_>>();
            // A complete set of unchanged source denominators is the production
            // graph, not a reconstructed exact topology. Reuse that source
            // literally, including its unpaired-external-leg convention.
            // `vertex_external_balance_info` is not a zero oracle here because
            // production external coordinates may include a dependent momentum;
            // this identity path must not pass through exact-source validation.
            if original_occurrences.len() == denominators.len()
                && original_occurrences
                    .keys()
                    .copied()
                    .collect::<BTreeSet<_>>()
                    == retained_source_edges
            {
                let local_edges_by_source = production_energy_map
                    .internal
                    .iter()
                    .map(|(local, source)| (EdgeIndex(*source), *local))
                    .collect::<BTreeMap<_, _>>();
                let exact_signatures = denominators
                    .iter()
                    .map(|denominator| {
                        let signature = production_parsed.internal_edges
                            [local_edges_by_source[&denominator.source_edge]]
                            .signature
                            .clone();
                        let source_momentum = FunctionBuilder::new(GS.emr_mom)
                            .add_arg(usize::from(denominator.source_edge))
                            .finish();
                        if denominator.momentum == source_momentum {
                            signature
                        } else {
                            signature.negated()
                        }
                    })
                    .collect();
                let exact_occurrences = production_energy_map
                    .internal
                    .iter()
                    .filter_map(|(local_edge_id, energy_edge_id)| {
                        original_occurrences.get(&EdgeIndex(*energy_edge_id)).map(
                            |original_occurrence| ExactParsedOccurrence {
                                local_edge_id: *local_edge_id,
                                energy_edge_id: *energy_edge_id,
                                original_occurrence: *original_occurrence,
                            },
                        )
                    })
                    .collect::<Vec<_>>();
                production_source.exact_denominators = Some(denominators);
                production_source.exact_signatures = exact_signatures;
                production_source.exact_parsed = Some(production_parsed);
                production_source.exact_local_to_original_occurrence = exact_occurrences
                    .iter()
                    .map(|occurrence| occurrence.original_occurrence)
                    .collect();
                production_source.exact_occurrences = exact_occurrences;
                production_source.exact_energy_edge_index_map = Some(production_energy_map);
                return Ok(production_source);
            }
        }
        let parent_loop_count = graph.loop_momentum_basis.loop_edges.len();
        let integral_coordinate = |coordinate: &Rational| {
            (coordinate.denominator_ref().to_i64() == Some(1))
                .then(|| coordinate.numerator_ref().to_i64())
                .flatten()
                .and_then(|numerator| i32::try_from(numerator).ok())
        };
        let production_loop_rows = production_signatures
            .iter()
            .map(|signature| {
                (0..parent_loop_count)
                    .map(|production_column| {
                        let coordinate = signature
                            .loop_signature
                            .iter()
                            .zip(&production_source.parent_loop_coordinates)
                            .fold(Rational::from(0), |sum, (coefficient, parent_coordinates)| {
                                sum + Rational::from(*coefficient)
                                    * parent_coordinates[production_column].clone()
                            });
                        integral_coordinate(&coordinate).ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "exact 4D denominator has non-integral production CFF coordinate {coordinate}"
                            ))
                        })
                    })
                    .collect()
            })
            .collect::<three_dimensional_reps::graph_io::Result<Vec<Vec<_>>>>()?;
        let (
            inner_loop_count,
            outer_loop_edges,
            edge_loop_coordinates,
            parent_loop_coordinates,
            exact_signatures,
        ) = if let Some(sub_lmb) = coordinate_lmb {
            if sub_lmb
                .loop_edges
                .iter()
                .any(|edge| !uv_edges.contains(edge))
            {
                return Err(GraphIoError::Source(format!(
                    "exact UV sub-LMB loop carriers {:?} are not all in the UV source edges {:?}",
                    sub_lmb.loop_edges,
                    uv_edges.iter().copied().sorted().collect::<Vec<_>>(),
                )));
            }
            let exact_signatures = coordinate_signatures
                .expect("an explicit exact coordinate LMB has coordinate signatures");
            let source_loop_count = sub_lmb.loop_edges.len();
            let source_loop_rows = exact_signatures
                .iter()
                .map(|signature| signature.loop_signature.clone())
                .collect::<Vec<_>>();
            select_exact_loop_basis(&source_loop_rows, source_loop_count)?;
            let edge_loop_coordinates = sub_lmb
                .edge_signatures
                .iter()
                .map(|(edge, signature)| {
                    (
                        edge,
                        signature
                            .internal
                            .iter()
                            .map(|coefficient| Rational::from(sign_to_i32(*coefficient)))
                            .collect(),
                    )
                })
                .collect();
            let parent_loop_coordinates = (0..source_loop_count)
                .map(|row| {
                    (0..source_loop_count)
                        .map(|column| Rational::from(i32::from(row == column)))
                        .collect()
                })
                .collect();
            (
                0,
                sub_lmb.loop_edges.iter().copied().collect(),
                edge_loop_coordinates,
                parent_loop_coordinates,
                exact_signatures,
            )
        } else {
            let (outer_basis_rows, active_loop_columns) =
                select_exact_loop_basis(&production_loop_rows, parent_loop_count)?;
            let inactive_loop_columns = (0..parent_loop_count)
                .filter(|column| !active_loop_columns.contains(column))
                .collect::<Vec<_>>();
            let mut source_basis_rows = inactive_loop_columns
                .iter()
                .map(|inactive_column| {
                    (0..parent_loop_count)
                        .map(|column| i32::from(column == *inactive_column))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            source_basis_rows.extend(outer_basis_rows.iter().cloned());
            let basis_matrix = (0..parent_loop_count)
                .map(|column| {
                    source_basis_rows
                        .iter()
                        .map(|row| Rational::from(row[column]))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let solve_coordinates = |row: &[Rational]| {
                if row.is_empty() {
                    return Ok(Vec::new());
                }
                solve_rational_system(basis_matrix.clone(), row.to_vec()).ok_or_else(|| {
                    GraphIoError::Source(
                        "failed to solve canonical exact-source loop coordinates".to_string(),
                    )
                })
            };
            let edge_loop_coordinates = graph
                .loop_momentum_basis
                .edge_signatures
                .iter()
                .map(|(edge_id, _)| {
                    let production_coordinates = production_source
                        .edge_loop_coordinates
                        .get(&edge_id)
                        .expect("production source contains every graph edge");
                    Ok((edge_id, solve_coordinates(production_coordinates)?))
                })
                .collect::<three_dimensional_reps::graph_io::Result<BTreeMap<_, _>>>()?;
            let parent_loop_coordinates = production_source
                .parent_loop_coordinates
                .iter()
                .map(|coordinates| solve_coordinates(coordinates))
                .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
            let inner_loop_count = inactive_loop_columns.len();
            let exact_signatures = production_signatures
                .iter()
                .zip(&production_loop_rows)
                .map(|(signature, production_row)| {
                    let coordinates = solve_coordinates(
                        &production_row
                            .iter()
                            .map(|coefficient| Rational::from(*coefficient))
                            .collect::<Vec<_>>(),
                    )?;
                    if coordinates[..inner_loop_count]
                        .iter()
                        .any(|coordinate| *coordinate != 0)
                    {
                        return Err(GraphIoError::Source(
                            "exact 4D denominator retains an inactive-loop coordinate".to_string(),
                        ));
                    }
                    let loop_signature = coordinates[inner_loop_count..]
                        .iter()
                        .map(|coordinate| {
                            integral_coordinate(coordinate).ok_or_else(|| {
                                GraphIoError::Source(format!(
                                    "exact 4D denominator has non-integral source coordinate {coordinate}"
                                ))
                            })
                        })
                        .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
                    Ok(MomentumSignature {
                        loop_signature,
                        external_signature: signature.external_signature.clone(),
                    })
                })
                .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
            let outer_loop_edges = outer_basis_rows
                .iter()
                .map(|basis_row| {
                    production_loop_rows
                        .iter()
                        .zip(denominators)
                        .filter_map(|(row, denominator)| {
                            (row == basis_row
                                || row.iter().zip(basis_row).all(
                                    |(coefficient, basis_coefficient)| {
                                        coefficient.checked_neg() == Some(*basis_coefficient)
                                    },
                                ))
                            .then_some(denominator.source_edge)
                        })
                        .min()
                        .expect("every exact source basis row comes from a denominator")
                })
                .collect();
            (
                inner_loop_count,
                outer_loop_edges,
                edge_loop_coordinates,
                parent_loop_coordinates,
                exact_signatures,
            )
        };
        let occurrence_keys = denominators
            .iter()
            .zip(&exact_signatures)
            .map(|(denominator, signature)| {
                (
                    usize::from(denominator.source_edge),
                    uv_edges.contains(&denominator.source_edge),
                    signature.clone(),
                    denominator.mass_squared.to_canonical_string(),
                    denominator.momentum.to_canonical_string(),
                    denominator.full_expr.to_canonical_string(),
                )
            })
            .collect::<Vec<_>>();
        let mut exact_local_to_original_occurrence = (0..denominators.len()).collect::<Vec<_>>();
        exact_local_to_original_occurrence
            .sort_by(|left, right| occurrence_keys[*left].cmp(&occurrence_keys[*right]));
        let active_cograph_edges = denominators
            .iter()
            .filter(|denominator| !uv_edges.contains(&denominator.source_edge))
            .map(|denominator| denominator.source_edge)
            .collect::<AHashSet<_>>();
        // A proper UV child is uncut: parent Cutkosky carriers belong to the
        // factorized outer graph, even when the child LMB lists their physical
        // IDs among its external coordinates. A cut incident to the child is
        // already retained once, as part of the complete physical crown.
        let initial_state_cut_edges = if coordinate_lmb.is_none() {
            graph
                .iter_edges_of(&graph.initial_state_cut)
                .map(|(_, edge_id, _)| edge_id)
                .collect::<AHashSet<_>>()
        } else {
            AHashSet::new()
        };
        let contract_edges = graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge_id, edge)| {
                (pair.is_paired()
                    && !edge.data.is_dummy
                    && !active_cograph_edges.contains(&edge_id)
                    && !initial_state_cut_edges.contains(&edge_id))
                .then_some(edge_id)
            })
            .collect();

        let mut source = Self {
            graph,
            contract_edges,
            initial_state_cut_edges,
            inner_loop_count,
            outer_loop_edges,
            edge_loop_coordinates,
            parent_loop_coordinates,
            uv_edges,
            exact_uv_boundary_hedges,
            exact_coordinate_lmb: coordinate_lmb.cloned(),
            exact_sub_lmb_frame: sub_lmb_frame,
            exact_denominators: Some(denominators),
            exact_signatures,
            exact_parsed: None,
            exact_local_to_original_occurrence,
            exact_occurrences: Vec::new(),
            exact_energy_edge_index_map: None,
        };

        let (exact_parsed, mut exact_occurrences, mut exact_energy_edge_index_map) =
            source.build_exact_parsed_graph_with_occurrence_order()?;
        let exact_uv_domains = exact_occurrences
            .iter()
            .map(|occurrence| {
                source
                    .uv_edges
                    .contains(&denominators[occurrence.original_occurrence].source_edge)
            })
            .collect::<Vec<_>>();
        let (exact_parsed, canonical_to_current) =
            canonicalize_exact_parsed_graph(exact_parsed, &exact_occurrences, &exact_uv_domains)?;
        let mut current_to_canonical = vec![usize::MAX; canonical_to_current.len()];
        for (canonical, current) in canonical_to_current.into_iter().enumerate() {
            current_to_canonical[current] = canonical;
        }
        if current_to_canonical.contains(&usize::MAX) {
            return Err(GraphIoError::Source(
                "exact topology canonicalization did not return a complete edge permutation"
                    .to_string(),
            ));
        }

        let physical_edge_count = graph.underlying.n_edges();
        exact_occurrences.sort_by_key(|occurrence| current_to_canonical[occurrence.local_edge_id]);
        // A rewritten exact source has a canonical occurrence-local namespace;
        // physical owner IDs remain available through the provenance mapper.
        let mut next_synthetic_energy_edge = physical_edge_count;
        let mut canonical_energy_ids = BTreeMap::new();
        for occurrence in &mut exact_occurrences {
            occurrence.local_edge_id = current_to_canonical[occurrence.local_edge_id];
            canonical_energy_ids
                .entry(occurrence.energy_edge_id)
                .or_insert_with(|| {
                    let canonical = next_synthetic_energy_edge;
                    next_synthetic_energy_edge += 1;
                    canonical
                });
            occurrence.energy_edge_id = canonical_energy_ids[&occurrence.energy_edge_id];
        }
        exact_energy_edge_index_map.internal = exact_energy_edge_index_map
            .internal
            .into_iter()
            .map(|(current, energy)| {
                (
                    current_to_canonical[current],
                    canonical_energy_ids.get(&energy).copied().unwrap_or(energy),
                )
            })
            .collect();
        exact_energy_edge_index_map.orientation_edge_count = next_synthetic_energy_edge;
        source.exact_local_to_original_occurrence = exact_occurrences
            .iter()
            .map(|occurrence| occurrence.original_occurrence)
            .collect();
        source.exact_occurrences = exact_occurrences;
        source.exact_energy_edge_index_map = Some(exact_energy_edge_index_map);
        source.exact_parsed = Some(exact_parsed);
        Ok(source)
    }

    pub(crate) fn edge_loop_coordinates(&self, edge_id: EdgeIndex) -> Option<&[Rational]> {
        self.edge_loop_coordinates.get(&edge_id).map(Vec::as_slice)
    }

    fn coordinate_lmb(&self) -> &LoopMomentumBasis {
        self.exact_coordinate_lmb
            .as_ref()
            .unwrap_or(&self.graph.loop_momentum_basis)
    }

    /// Physical EMRs retained by the factorized outer graph of a proper UV
    /// child. They are external to its CFF numerator rank and stay T-inert.
    pub(crate) fn factorized_external_emr_edges(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.graph
            .underlying
            .iter_edges()
            .filter_map(move |(_, edge_id, edge)| {
                (self.exact_coordinate_lmb.is_some()
                    && !edge.data.is_dummy
                    && !self.uv_edges.contains(&edge_id))
                .then_some(edge_id)
            })
    }

    pub(crate) fn reconstructible_outer_loop_coordinates(
        &self,
        edge_id: EdgeIndex,
    ) -> Option<&[Rational]> {
        let coordinates = self.edge_loop_coordinates(edge_id)?;
        coordinates[..self.inner_loop_count]
            .iter()
            .all(|coordinate| *coordinate == 0)
            .then_some(&coordinates[self.inner_loop_count..])
    }

    fn outer_loop_signature(
        &self,
        edge_id: EdgeIndex,
    ) -> three_dimensional_reps::graph_io::Result<Vec<i32>> {
        self.edge_loop_coordinates(edge_id)
            .ok_or_else(|| {
                GraphIoError::Source(format!(
                    "missing contracted-source coordinates for edge {}",
                    usize::from(edge_id)
                ))
            })?
            .iter()
            .skip(self.inner_loop_count)
            .map(|coordinate| {
                let numerator = coordinate.numerator_ref().to_i64().ok_or_else(|| {
                    GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} is outside the i64 range"
                    ))
                })?;
                let denominator = coordinate.denominator_ref().to_i64().ok_or_else(|| {
                    GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} has an out-of-range denominator"
                    ))
                })?;
                if denominator != 1 {
                    return Err(GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} is not integral"
                    )));
                }
                i32::try_from(numerator).map_err(|_| {
                    GraphIoError::Source(format!(
                        "contracted-source coordinate {coordinate} is outside the i32 range"
                    ))
                })
            })
            .collect()
    }

    fn exact_source_routing_sign(
        &self,
        source_edge: EdgeIndex,
        uses_uv_loop_basis: bool,
        canonical_exact_signature: &MomentumSignature,
    ) -> three_dimensional_reps::graph_io::Result<i32> {
        let source_signature = self.outer_loop_signature(source_edge)?;
        if source_signature.len() != canonical_exact_signature.loop_signature.len() {
            return Err(GraphIoError::Source(format!(
                "source edge {} and its rewritten exact denominator use incompatible loop-coordinate dimensions",
                usize::from(source_edge),
            )));
        }
        let denominators = self
            .exact_denominators
            .expect("exact routing is requested only for exact denominators");
        let opposite_domain_rows = self
            .exact_local_to_original_occurrence
            .iter()
            .copied()
            .filter(|original| {
                self.uv_edges.contains(&denominators[*original].source_edge) != uses_uv_loop_basis
            })
            .map(|original| {
                self.exact_signatures[original]
                    .canonical_up_to_sign()
                    .0
                    .loop_signature
            })
            .collect::<Vec<_>>();
        let opposite_rank = exact_integer_rank(&opposite_domain_rows);
        let mut compatible_signs = Vec::new();
        // Source incidence remains authoritative. This rank test never infers
        // endpoints: it only determines whether the rewritten denominator is
        // routed with or against that incidence after quotienting out loop
        // directions belonging to the factorized opposite topology domain.
        for sign in [1_i32, -1] {
            let difference = canonical_exact_signature
                .loop_signature
                .iter()
                .zip(&source_signature)
                .map(|(exact, source)| {
                    exact
                        .checked_sub(sign.checked_mul(*source).ok_or_else(|| {
                            GraphIoError::Source(
                                "exact source routing coefficient overflow".to_string(),
                            )
                        })?)
                        .ok_or_else(|| {
                            GraphIoError::Source(
                                "exact source routing coefficient overflow".to_string(),
                            )
                        })
                })
                .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
            let mut rows_with_difference = opposite_domain_rows.clone();
            rows_with_difference.push(difference);
            if exact_integer_rank(&rows_with_difference) == opposite_rank {
                compatible_signs.push(sign);
            }
        }
        match compatible_signs.as_slice() {
            [sign] => Ok(*sign),
            _ => Err(GraphIoError::Source(format!(
                "source edge {} has no unique routing for rewritten exact signature {:?} modulo the factorized opposite-domain loop span",
                usize::from(source_edge),
                canonical_exact_signature.loop_signature,
            ))),
        }
    }

    fn contracts_edge(&self, edge_id: EdgeIndex) -> bool {
        self.contract_edges.contains(&edge_id) && !self.initial_state_cut_edges.contains(&edge_id)
    }

    pub(crate) fn physical_energy_edge_index_map(&self) -> Option<EnergyEdgeIndexMap> {
        let denominators = self.exact_denominators?;
        let parsed = self.exact_parsed_graph().ok()?;
        let mut internal = self
            .exact_occurrences
            .iter()
            .map(|occurrence| {
                (
                    occurrence.energy_edge_id,
                    usize::from(denominators[occurrence.original_occurrence].source_edge),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let energy_map = self.exact_energy_edge_index_map.as_ref()?;
        for (cut_edge, source_edge) in parsed
            .initial_state_cut_edges
            .iter()
            .zip(self.initial_state_cut_edges.iter().copied().sorted())
        {
            internal.insert(
                energy_map.internal[&cut_edge.edge_id],
                usize::from(source_edge),
            );
        }
        Some(EnergyEdgeIndexMap {
            internal,
            // Exact-source generation already remaps external coordinates to
            // physical graph edge IDs. This second map changes only the
            // occurrence-local internal-energy namespace.
            external: BTreeMap::new(),
            orientation_edge_count: self.graph.underlying.n_edges(),
        })
    }

    /// Project exact occurrence support onto every physical owner of the same
    /// rewritten denominator channel. Energy evaluation remains occurrence
    /// local, while Cutkosky support is neutral provenance: equal `D(Q)` and
    /// `D(-Q)` factors with equal mass and topology domain support every
    /// distinct source edge which instantiated that algebraic channel.
    pub(crate) fn physical_cut_support_edge_index_map(
        &self,
    ) -> Option<BTreeMap<usize, Vec<EdgeIndex>>> {
        let denominators = self.exact_denominators?;
        let channels = self
            .exact_occurrences
            .iter()
            .map(|occurrence| {
                let denominator = &denominators[occurrence.original_occurrence];
                let (signature, _) =
                    self.exact_signatures[occurrence.original_occurrence].canonical_up_to_sign();
                (
                    (
                        self.uv_edges.contains(&denominator.source_edge),
                        signature,
                        denominator.mass_squared.to_canonical_string(),
                    ),
                    denominator.source_edge,
                    occurrence.energy_edge_id,
                )
            })
            .collect::<Vec<_>>();
        let mut owners_by_channel = BTreeMap::<_, BTreeSet<EdgeIndex>>::new();
        for (channel, owner, _) in &channels {
            owners_by_channel
                .entry(channel.clone())
                .or_default()
                .insert(*owner);
        }

        let mut support = channels
            .iter()
            .map(|(channel, _, energy_edge_id)| {
                (
                    *energy_edge_id,
                    owners_by_channel[channel].iter().copied().collect(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        for (occurrence, owner) in self.physical_energy_edge_index_map()?.internal {
            support
                .entry(occurrence)
                .or_insert_with(|| vec![EdgeIndex(owner)]);
        }
        Some(support)
    }

    pub(crate) fn physical_linear_surface(&self, surface: &LinearSurface) -> Option<LinearSurface> {
        let denominators = self.exact_denominators?;
        let edge_map = self.physical_energy_edge_index_map()?;
        let all_physical = surface
            .expression
            .internal_terms
            .iter()
            .all(|(edge_id, _)| {
                let energy_edge_id = usize::from(*edge_id);
                if let Some(denominator) = self
                    .exact_occurrences
                    .iter()
                    .find(|occurrence| occurrence.energy_edge_id == energy_edge_id)
                    .and_then(|occurrence| denominators.get(occurrence.original_occurrence))
                {
                    !self.uv_edges.contains(&denominator.source_edge)
                        && denominator.is_original_graph_denominator(self.graph)
                } else {
                    edge_map.internal.contains_key(&energy_edge_id)
                }
            });
        all_physical.then(|| {
            let mut surface = surface.clone();
            surface.expression = surface
                .expression
                .remap_energy_edges(&edge_map.internal, &BTreeMap::new());
            surface
        })
    }

    fn energy_mapper(
        &self,
        inactive_loop_count: usize,
        exact_ose_replacements: Vec<Replacement>,
    ) -> ExactSourceEnergyMapper {
        let coordinate_lmb = self.coordinate_lmb();
        // A proper child maps only shell-owned EMRs. Crown and remote parent
        // EMRs stay as factorized, T-inert symbols until the outer graph maps
        // them; replacing an all-zero child row would incorrectly erase them.
        // A nested exact source retains the full enclosing scope solely so its
        // omitted prefix can be contracted. That prefix has zero rows in the
        // quotient LMB and is therefore outside this mapper's domain. A pinched
        // edge of the active quotient still has a nonzero row; when it has no
        // exact occurrence below, its ordinary quotient signature is retained.
        let edge_coordinates = self
            .graph
            .underlying
            .iter_edges()
            .filter(|(pair, edge_id, edge)| {
                pair.is_paired()
                    && !edge.data.is_dummy
                    && (self.exact_coordinate_lmb.is_none()
                        || (self.uv_edges.contains(edge_id)
                            && coordinate_lmb.edge_signatures[*edge_id]
                                .internal
                                .iter()
                                .any(|coefficient| *coefficient != SignOrZero::Zero)))
            })
            .map(|(_, edge_id, _)| {
                let signature = &coordinate_lmb.edge_signatures[edge_id];
                let taylor_signature = (self.exact_sub_lmb_frame
                    == Some(ExactUvSubLmbFrame::TaylorVacuum))
                .then(|| {
                    self.exact_denominators
                        .expect("a Taylor-vacuum source has exact denominators")
                        .iter()
                        .enumerate()
                        .find_map(|(occurrence, denominator)| {
                            (denominator.source_edge == edge_id)
                                .then(|| &self.exact_signatures[occurrence])
                        })
                })
                .flatten();
                let mut coordinates = taylor_signature.map_or_else(
                    || self.edge_loop_coordinates[&edge_id].clone(),
                    |signature| {
                        signature
                            .loop_signature
                            .iter()
                            .copied()
                            .map(Rational::from)
                            .collect()
                    },
                );
                let external_terms = if taylor_signature.is_some() {
                    Vec::new()
                } else if self.initial_state_cut_edges.contains(&edge_id) {
                    let alias_signature = MomentumSignature {
                        loop_signature: Vec::new(),
                        external_signature: (&signature.external)
                            .into_iter()
                            .map(sign_to_i32)
                            .collect(),
                    };
                    let (external_id, external_sign) = initial_state_cut_external_alias(
                        usize::from(edge_id),
                        &alias_signature,
                    )
                    .expect(
                        "exact source cut aliases were validated while building OSE replacements",
                    );
                    let external_edge = self
                        .coordinate_lmb()
                        .ext_edges
                        .iter()
                        .copied()
                        .nth(external_id)
                        .expect("validated cut alias refers to an existing external coordinate");
                    coordinates.fill(Rational::from(0));
                    vec![(
                        external_edge,
                        match external_sign {
                            -1 => SignOrZero::Minus,
                            1 => SignOrZero::Plus,
                            _ => unreachable!("initial-state cut aliases have unit sign"),
                        },
                    )]
                } else {
                    coordinate_lmb
                        .ext_edges
                        .iter()
                        .copied()
                        .zip(signature.external.iter().copied())
                        .collect()
                };
                (edge_id, coordinates, external_terms)
            })
            .collect::<Vec<_>>();
        let mut source_edge_occurrences =
            BTreeMap::<EdgeIndex, Vec<ExactSourceOwnerOccurrence>>::new();
        if let Some(denominators) = self.exact_denominators {
            let parsed = self
                .exact_parsed
                .as_ref()
                .expect("exact denominator source has a canonical parsed graph");
            for occurrence in &self.exact_occurrences {
                let local = occurrence.local_edge_id;
                let original = occurrence.original_occurrence;
                let denominator = &denominators[original];
                // Source incidence—not a literal momentum search—owns every
                // serial occurrence. Canonicalization may reverse the complete
                // routed denominator; retain that sign beside the raw momentum
                // so an odd tagged hard numerator can compose it explicitly.
                let raw_to_parsed_sign = if parsed.internal_edges[local].signature
                    == self.exact_signatures[original]
                {
                    1
                } else if parsed.internal_edges[local].signature
                    == self.exact_signatures[original].negated()
                {
                    -1
                } else {
                    unreachable!(
                        "canonical exact occurrence routing differs by more than an overall sign"
                    )
                };
                source_edge_occurrences
                    .entry(denominator.source_edge)
                    .or_default()
                    .push(ExactSourceOwnerOccurrence {
                        energy_edge_id: occurrence.energy_edge_id,
                        raw_momentum: denominator.momentum.clone(),
                        raw_to_parsed_sign,
                    });
            }
        }
        ExactSourceEnergyMapper {
            inactive_loop_count,
            parent_loop_coordinates: self.parent_loop_coordinates.clone(),
            parent_loop_edges: self.coordinate_lmb().loop_edges.iter().copied().collect(),
            edge_coordinates,
            source_edge_occurrences,
            cut_alias_edges: self.initial_state_cut_edges.clone(),
            exact_ose_replacements,
        }
    }

    pub(crate) fn exact_source_energy_mapper(&self) -> Option<ExactSourceEnergyMapper> {
        Some(self.energy_mapper(self.inner_loop_count, self.exact_ose_replacements()?))
    }

    pub(crate) fn exact_ose_replacements(&self) -> Option<Vec<Replacement>> {
        let denominators = self.exact_denominators?;
        let mut replacements = self
            .exact_occurrences
            .iter()
            .map(|occurrence| {
                let denominator = &denominators[occurrence.original_occurrence];
                // Keep unchanged physical propagators in GammaLoop's OSE
                // dialect so a later local-3D UV operation can still rescale
                // them. Rewritten exact denominators have no physical OSE
                // identity and therefore retain their literal square root.
                let on_shell_energy = if denominator.is_original_graph_denominator(self.graph) {
                    crate::utils::ose_atom_from_index(denominator.source_edge)
                } else {
                    denominator.on_shell_energy()
                };
                Replacement::new(
                    crate::utils::ose_atom_from_index(EdgeIndex(occurrence.energy_edge_id))
                        .to_pattern(),
                    on_shell_energy.to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let parsed = self.exact_parsed_graph().ok()?;
        let energy_map = self.exact_energy_edge_index_map.as_ref()?;
        replacements.extend(
            parsed
                .initial_state_cut_edges
                .iter()
                .zip(self.initial_state_cut_edges.iter().copied().sorted())
                .map(|(cut_edge, source_edge)| {
                    Replacement::new(
                        crate::utils::ose_atom_from_index(EdgeIndex(
                            energy_map.internal[&cut_edge.edge_id],
                        ))
                        .to_pattern(),
                        crate::utils::ose_atom_from_index(source_edge).to_pattern(),
                    )
                }),
        );
        Some(replacements)
    }

    pub(crate) fn exact_inverse_energy_product(&self) -> Option<Atom> {
        Some(
            self.exact_denominators?
                .iter()
                .map(|denominator| {
                    // Match the residue-map dialect above: every occurrence
                    // factor must enter a subsequent local-3D T operation in
                    // the same rescalable form as its pole energy.
                    let on_shell_energy = if denominator.is_original_graph_denominator(self.graph) {
                        crate::utils::ose_atom_from_index(denominator.source_edge)
                    } else {
                        denominator.on_shell_energy()
                    };
                    -Atom::num(2) * on_shell_energy
                })
                .reduce(|product, energy| product * energy)
                .map(|product| Atom::one() / product)
                .unwrap_or_else(Atom::one),
        )
    }

    pub(crate) fn active_loop_count(&self) -> usize {
        self.outer_loop_edges.len()
    }

    /// Rebuild the exact source with one representative for each immutable
    /// physical denominator owner. The original graph incidence is reused;
    /// no topology is inferred from momentum signatures.
    #[cfg(test)]
    pub(crate) fn physical_owner_skeleton(
        &self,
    ) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        let denominators = self.exact_denominators.ok_or_else(|| {
            GraphIoError::Source("a physical-owner skeleton requires exact denominators".into())
        })?;
        let mut representatives = BTreeMap::<EdgeIndex, FourDDenominator>::new();
        for denominator in denominators {
            match representatives.entry(denominator.source_edge) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(denominator.clone());
                }
                std::collections::btree_map::Entry::Occupied(mut entry)
                    if denominator.is_original_graph_denominator(self.graph)
                        && !entry.get().is_original_graph_denominator(self.graph) =>
                {
                    entry.insert(denominator.clone());
                }
                std::collections::btree_map::Entry::Occupied(_) => {}
            }
        }
        let representatives = representatives.into_values().collect::<Vec<_>>();
        physical_owner_skeleton_from_denominators(
            self.graph,
            &representatives,
            self.uv_edges.iter().copied(),
            self.exact_uv_boundary_hedges.iter().copied(),
            self.exact_coordinate_lmb.as_ref(),
            self.exact_sub_lmb_frame,
        )
    }

    pub(crate) fn contract_subgraph(&self) -> SuBitGraph {
        let mut subgraph: SuBitGraph = self.graph.empty_subgraph();
        for (pair, edge_id, _) in self.graph.underlying.iter_edges() {
            if pair.is_paired() && self.contracts_edge(edge_id) {
                subgraph.add(pair);
            }
        }
        subgraph
    }

    fn exact_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        self.exact_parsed.clone().ok_or_else(|| {
            GraphIoError::Source("exact parsed graph was not initialized".to_string())
        })
    }

    fn build_exact_parsed_graph_with_occurrence_order(
        &self,
    ) -> three_dimensional_reps::graph_io::Result<(
        ParsedGraph,
        Vec<ExactParsedOccurrence>,
        EnergyEdgeIndexMap,
    )> {
        let denominators = self
            .exact_denominators
            .expect("exact parsed source has exact denominators");
        let coordinate_lmb = self.coordinate_lmb();

        // Denominators absent from this additive 4D term are contracted only
        // in the source topology. The typed atom remains untouched and keeps
        // exact multiplicity and denominator expressions.
        let mut parent = (0..self.graph.n_nodes()).collect_vec();
        for (pair, edge_id, edge_data) in self.graph.underlying.iter_edges() {
            if edge_data.data.is_dummy || !self.contracts_edge(edge_id) {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                union_parent(
                    &mut parent,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink)),
                );
            }
        }

        let node_ids = self
            .graph
            .underlying
            .iter_nodes()
            .map(|(node_id, _, _)| node_id)
            .sorted()
            .collect::<Vec<_>>();
        let mut root_to_internal = BTreeMap::<usize, usize>::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            if !root_to_internal.contains_key(&root) {
                root_to_internal.insert(root, root_to_internal.len());
            }
        }
        let node_to_internal = node_ids
            .iter()
            .map(|node_id| {
                let root = find_parent(&mut parent, usize::from(*node_id));
                (*node_id, root_to_internal[&root])
            })
            .collect::<BTreeMap<_, _>>();
        // An exact term retains one energy ID per denominator occurrence, but
        // repeated occurrences with the same rewritten rational propagator
        // and physical incidence encode one powered line. Subdivide that
        // incidence before CFF recursion. UV and cograph occurrences use
        // separate node domains even if all other data agree.
        let mut occurrence_groups =
            BTreeMap::<EdgeIndex, (bool, MomentumSignature, String, i32, Vec<usize>)>::new();
        for (occurrence, original) in self
            .exact_local_to_original_occurrence
            .iter()
            .copied()
            .enumerate()
        {
            let denominator = &denominators[original];
            let source_edge = denominator.source_edge;
            let (_, pair) = self.graph[&source_edge];
            if !matches!(pair, HedgePair::Paired { .. }) {
                return Err(GraphIoError::Source(format!(
                    "4D denominator source edge {} is not a paired internal edge",
                    usize::from(source_edge),
                )));
            }
            let uses_uv_loop_basis = self.uv_edges.contains(&source_edge);
            let (signature, _) = self.exact_signatures[original].canonical_up_to_sign();
            if self
                .outer_loop_signature(source_edge)?
                .iter()
                .all(|coefficient| *coefficient == 0)
            {
                return Err(GraphIoError::Source(format!(
                    "loop-active exact denominator has loop-independent source edge {}",
                    usize::from(source_edge),
                )));
            }
            let source_routing_sign =
                self.exact_source_routing_sign(source_edge, uses_uv_loop_basis, &signature)?;
            let mass = denominator.mass_squared.to_canonical_string();
            match occurrence_groups.entry(source_edge) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert((
                        uses_uv_loop_basis,
                        signature,
                        mass,
                        source_routing_sign,
                        vec![occurrence],
                    ));
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    let (domain, group_signature, group_mass, _, members) = entry.get_mut();
                    if (*domain, &*group_signature, group_mass.as_str())
                        != (uses_uv_loop_basis, &signature, mass.as_str())
                    {
                        return Err(GraphIoError::Source(format!(
                            "source edge {} instantiates more than one rewritten rational denominator in one additive 4D term",
                            usize::from(source_edge),
                        )));
                    }
                    members.push(occurrence);
                }
            }
        }
        for (source_edge, (_, _, _, _, members)) in &mut occurrence_groups {
            members.sort_by_key(|occurrence| {
                let original = self.exact_local_to_original_occurrence[*occurrence];
                (
                    !denominators[original].is_original_graph_denominator(self.graph),
                    *occurrence,
                )
            });
            if members.is_empty() {
                return Err(GraphIoError::Source(format!(
                    "source edge {} has an empty exact denominator group",
                    usize::from(*source_edge),
                )));
            }
        }
        let active_uv_edges = occurrence_groups
            .iter()
            .filter(|(_, (uses_uv_loop_basis, _, _, _, _))| *uses_uv_loop_basis)
            .map(|(source_edge, _)| *source_edge)
            .collect::<BTreeSet<_>>();
        let mut uv_parent = (0..self.graph.n_nodes()).collect_vec();
        for (pair, edge_id, edge) in self.graph.underlying.iter_edges() {
            if edge.data.is_dummy
                || !self.uv_edges.contains(&edge_id)
                || active_uv_edges.contains(&edge_id)
                || self.initial_state_cut_edges.contains(&edge_id)
            {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                union_parent(
                    &mut uv_parent,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink)),
                );
            }
        }
        let mut uv_roots = BTreeSet::new();
        for source_edge in &active_uv_edges {
            let (_, HedgePair::Paired { source, sink }) = self.graph[source_edge] else {
                unreachable!("exact source edges were validated as paired")
            };
            for node in [self.graph.node_id(source), self.graph.node_id(sink)] {
                uv_roots.insert(find_parent(&mut uv_parent, usize::from(node)));
            }
        }
        let uv_root_to_internal = uv_roots
            .into_iter()
            .enumerate()
            .map(|(offset, root)| (root, root_to_internal.len() + offset))
            .collect::<BTreeMap<_, _>>();
        let uv_node_to_internal = node_ids
            .iter()
            .filter_map(|node| {
                let root = find_parent(&mut uv_parent, usize::from(*node));
                uv_root_to_internal
                    .get(&root)
                    .copied()
                    .map(|internal| (*node, internal))
            })
            .collect::<BTreeMap<_, _>>();
        let mut next_node = root_to_internal.len() + uv_root_to_internal.len();
        let mut power_chain_nodes = BTreeMap::new();
        let loop_names = self
            .outer_loop_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();
        let mut internal_edges = Vec::new();
        let mut exact_occurrences = Vec::with_capacity(denominators.len());
        let mut initial_state_cut_edges = Vec::new();
        let mut internal_energy_edges = BTreeMap::new();
        let mut next_synthetic_energy_edge = self.graph.underlying.n_edges();
        let mut group_id = 0;
        for (pair, edge_index, edge) in self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
        {
            if edge.data.is_dummy {
                continue;
            }
            let HedgePair::Paired { source, sink } = pair else {
                continue;
            };
            if self.initial_state_cut_edges.contains(&edge_index) {
                let signature = &coordinate_lmb.edge_signatures[edge_index];
                let local_edge_id = internal_edges.len();
                internal_edges.push(ParsedGraphInternalEdge {
                    edge_id: local_edge_id,
                    tail: node_to_internal[&self.graph.node_id(source)],
                    head: node_to_internal[&self.graph.node_id(sink)],
                    label: edge.data.name.value.clone(),
                    mass_key: Some(edge.data.particle.mass_atom().to_canonical_string()),
                    signature: MomentumSignature {
                        loop_signature: self.outer_loop_signature(edge_index)?,
                        external_signature: (&signature.external)
                            .into_iter()
                            .map(sign_to_i32)
                            .collect(),
                    },
                    had_pow: false,
                });
                let (external_id, external_sign) = initial_state_cut_external_alias(
                    usize::from(edge_index),
                    &internal_edges[local_edge_id].signature,
                )?;
                initial_state_cut_edges.push(ParsedGraphInitialStateCutEdge {
                    edge_id: local_edge_id,
                    external_id,
                    external_sign,
                });
                internal_energy_edges.insert(local_edge_id, usize::from(edge_index));
                continue;
            }

            let Some((uses_uv_loop_basis, canonical_signature, mass, source_routing_sign, members)) =
                occurrence_groups.remove(&edge_index)
            else {
                continue;
            };
            let node_map = if uses_uv_loop_basis {
                &uv_node_to_internal
            } else {
                &node_to_internal
            };
            let (incidences, auxiliary_nodes) = serial_power_chain_incidences(
                node_map[&self.graph.node_id(source)],
                node_map[&self.graph.node_id(sink)],
                1,
                members.len(),
                &mut next_node,
            );
            power_chain_nodes.extend(auxiliary_nodes.into_iter().enumerate().map(
                |(position, node)| {
                    (
                        format!("__gammaloop_exact_power_{group_id}_{position}"),
                        node,
                    )
                },
            ));
            let representative = &denominators[self.exact_local_to_original_occurrence[members[0]]];
            let signature = if source_routing_sign == 1 {
                canonical_signature
            } else {
                canonical_signature.negated()
            };
            // The rational denominator is even in its momentum.  Rewriting a
            // physical D(Q) as D(-Q') must therefore retain the source
            // particle's mass identity even when `Q'` is not the literal
            // momentum wrapper of this owner.  Momentum spelling still
            // controls occurrence provenance below; it cannot split repeated
            // denominators into different CFF mass classes.
            let mass_key = if representative.mass_squared == edge.data.mass_atom().pow(2) {
                edge.data.particle.mass_atom().to_canonical_string()
            } else {
                mass
            };
            for (position, (occurrence, (tail, head))) in
                members.into_iter().zip(incidences).enumerate()
            {
                let original_occurrence = self.exact_local_to_original_occurrence[occurrence];
                let denominator = &denominators[original_occurrence];
                let energy_edge_id =
                    if position == 0 && denominator.is_original_graph_denominator(self.graph) {
                        usize::from(edge_index)
                    } else {
                        let energy_edge_id = next_synthetic_energy_edge;
                        next_synthetic_energy_edge += 1;
                        energy_edge_id
                    };
                let local_edge_id = internal_edges.len();
                internal_edges.push(ParsedGraphInternalEdge {
                    edge_id: local_edge_id,
                    tail,
                    head,
                    label: if position == 0 {
                        edge.data.name.value.clone()
                    } else {
                        format!("{}__exact_power_{position}", edge.data.name.value)
                    },
                    mass_key: Some(mass_key.clone()),
                    signature: signature.clone(),
                    had_pow: false,
                });
                exact_occurrences.push(ExactParsedOccurrence {
                    local_edge_id,
                    energy_edge_id,
                    original_occurrence,
                });
                internal_energy_edges.insert(local_edge_id, energy_edge_id);
            }
            group_id += 1;
        }
        if !occurrence_groups.is_empty() {
            return Err(GraphIoError::Source(format!(
                "exact denominator groups retain unknown source edges {:?}",
                occurrence_groups
                    .keys()
                    .map(|edge| usize::from(*edge))
                    .collect::<Vec<_>>(),
            )));
        }

        let mut next_external_id = 10_000_000usize;
        let mut external_edges = self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_id, _)| *edge_id)
            .filter_map(|(pair, edge_id, edge)| match pair {
                HedgePair::Unpaired { hedge, flow } if !edge.data.is_dummy => {
                    let node = node_to_internal[&self.graph.node_id(hedge)];
                    let (source, destination) = match flow {
                        Flow::Source => (Some(node), None),
                        Flow::Sink => (None, Some(node)),
                    };
                    let incidence_sign = match flow {
                        Flow::Source => -1,
                        Flow::Sink => 1,
                    };
                    let external_edge = ParsedGraphExternalEdge {
                        edge_id: next_external_id,
                        source,
                        destination,
                        label: edge.data.name.value.clone(),
                        external_coefficients: (&coordinate_lmb.edge_signatures[edge_id].external)
                            .into_iter()
                            .map(|coefficient| {
                                let coefficient = sign_to_i32(coefficient);
                                if self.uv_edges.is_empty() {
                                    coefficient
                                } else {
                                    incidence_sign * coefficient
                                }
                            })
                            .collect(),
                    };
                    next_external_id += 1;
                    Some(external_edge)
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        for hedge in &self.exact_uv_boundary_hedges {
            let edge_id = self.graph.underlying[hedge];
            let node = self.graph.node_id(*hedge);
            let Some(&internal_node) = uv_node_to_internal.get(&node) else {
                return Err(GraphIoError::Source(format!(
                    "UV boundary hedge {hedge} of source edge {} is not incident to the retained UV source minor",
                    usize::from(edge_id),
                )));
            };
            let signature = &coordinate_lmb.edge_signatures[edge_id];
            if signature
                .internal
                .iter()
                .any(|coefficient| *coefficient != SignOrZero::Zero)
            {
                return Err(GraphIoError::Source(format!(
                    "UV boundary hedge {hedge} of source edge {} carries a loop coordinate in the exact source LMB",
                    usize::from(edge_id),
                )));
            }
            let (source, destination, incidence_sign) = match self.graph.underlying.flow(*hedge) {
                Flow::Source => (Some(internal_node), None, -1),
                Flow::Sink => (None, Some(internal_node), 1),
            };
            external_edges.push(ParsedGraphExternalEdge {
                edge_id: next_external_id,
                source,
                destination,
                label: self.graph.underlying[edge_id].name.value.clone(),
                external_coefficients: (&signature.external)
                    .into_iter()
                    .map(|coefficient| incidence_sign * sign_to_i32(coefficient))
                    .collect(),
            });
            next_external_id += 1;
        }
        // Contracted cograph components without any exact denominator are
        // multiplicative identities. Keeping their external-only vertices in
        // the causal graph would manufacture a zero E-surface. Retain every
        // disconnected component that contains an exact denominator (or cut
        // carrier), and only discard external incidence on zero-denominator
        // components.
        let active_nodes = internal_edges
            .iter()
            .flat_map(|edge| [edge.tail, edge.head])
            .collect::<BTreeSet<_>>();
        external_edges.retain(|edge| {
            edge.source
                .into_iter()
                .chain(edge.destination)
                .any(|node| active_nodes.contains(&node))
        });

        let parsed = ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges,
            loop_names,
            external_names: coordinate_lmb
                .ext_edges
                .iter()
                .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
                .collect(),
            node_name_to_internal: root_to_internal
                .into_iter()
                .map(|(root, node)| (format!("n{root}"), node))
                .chain(
                    uv_node_to_internal
                        .into_iter()
                        .map(|(node, internal)| (format!("uv{}", usize::from(node)), internal)),
                )
                .chain(power_chain_nodes)
                .collect(),
        };
        let validation = validate_parsed_graph(&parsed);
        if !validation.vertex_balance_violations.is_empty() {
            return Err(GraphIoError::Source(format!(
                "source-backed exact topology violates loop-momentum conservation at vertices {:?}",
                validation.vertex_balance_violations,
            )));
        }
        let initial_state_external_ids = parsed
            .initial_state_cut_edges
            .iter()
            .map(|edge| edge.external_id)
            .collect::<BTreeSet<_>>();
        // External coordinates can include a dependent momentum, so their
        // per-vertex residual is relative to the source graph rather than
        // absolutely zero. Contracted production nodes inherit the sum of their
        // source-node residuals. A proper UV source instead retains its complete
        // crown and restores original owner signatures on the same incidence;
        // serial dotted copies then cancel on their auxiliary nodes.
        let mut source_external_balance = BTreeMap::<usize, Vec<i32>>::new();
        if self.uv_edges.is_empty() {
            let production_parsed =
                GraphThreeDSource::new(self.graph, &[])?.to_three_d_parsed_graph()?;
            let production_external_balance =
                validate_parsed_graph(&production_parsed).vertex_external_balance_info;
            for node in &node_ids {
                let source_node = production_parsed
                    .node_name_to_internal
                    .get(&format!("n{}", usize::from(*node)))
                    .ok_or_else(|| {
                        GraphIoError::Source(format!(
                            "production source is missing original node {}",
                            usize::from(*node),
                        ))
                    })?;
                let Some(coefficients) = production_external_balance.get(source_node) else {
                    continue;
                };
                let contracted = source_external_balance
                    .entry(node_to_internal[node])
                    .or_insert_with(|| vec![0; parsed.external_names.len()]);
                for (contracted, source) in contracted.iter_mut().zip(coefficients) {
                    *contracted += source;
                }
            }
        } else if self.exact_coordinate_lmb.is_some()
            && self.exact_sub_lmb_frame == Some(ExactUvSubLmbFrame::RetainedPhysicalCrown)
        {
            let mut original_source = parsed.clone();
            for occurrence in &exact_occurrences {
                let source_edge = denominators[occurrence.original_occurrence].source_edge;
                original_source.internal_edges[occurrence.local_edge_id]
                    .signature
                    .external_signature = (&coordinate_lmb.edge_signatures[source_edge].external)
                    .into_iter()
                    .map(sign_to_i32)
                    .collect();
            }
            source_external_balance =
                validate_parsed_graph(&original_source).vertex_external_balance_info;
        }
        let unresolved_external_balance = active_nodes
            .into_iter()
            .filter_map(|node| {
                let actual = validation
                    .vertex_external_balance_info
                    .get(&node)
                    .cloned()
                    .unwrap_or_else(|| vec![0; parsed.external_names.len()]);
                let source = source_external_balance
                    .remove(&node)
                    .unwrap_or_else(|| vec![0; parsed.external_names.len()]);
                let mut coefficients = actual
                    .into_iter()
                    .zip(source)
                    .map(|(actual, source)| actual - source)
                    .collect::<Vec<_>>();
                for external_id in &initial_state_external_ids {
                    if let Some(coefficient) = coefficients.get_mut(*external_id) {
                        *coefficient = 0;
                    }
                }
                coefficients
                    .iter()
                    .any(|coefficient| *coefficient != 0)
                    .then_some((node, coefficients))
            })
            .collect::<BTreeMap<_, _>>();
        if !unresolved_external_balance.is_empty() {
            return Err(GraphIoError::Source(format!(
                "source-backed exact topology has external momentum imbalance {:?}; retain the corresponding source-graph UV boundary hedges explicitly",
                unresolved_external_balance,
            )));
        }
        let energy_edge_index_map = EnergyEdgeIndexMap {
            internal: internal_energy_edges,
            external: coordinate_lmb
                .ext_edges
                .iter()
                .enumerate()
                .map(|(external_id, edge_id)| (external_id, usize::from(*edge_id)))
                .collect(),
            orientation_edge_count: next_synthetic_energy_edge,
        };
        Ok((parsed, exact_occurrences, energy_edge_index_map))
    }
}

impl FourDDenominator {
    fn is_original_graph_denominator(&self, graph: &Graph) -> bool {
        let source_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(self.source_edge))
            .finish();
        (self.momentum == source_momentum || self.momentum == -source_momentum)
            && self.mass_squared == graph.underlying[self.source_edge].mass_atom().pow(2)
    }

    fn on_shell_energy(&self) -> Atom {
        let spatial_square = (1..=3).fold(Atom::Zero, |sum, spatial_index| {
            let component_index = Atom::from(ExpandedIndex::from_iter([spatial_index]));
            let component = self
                .momentum
                .replace(function!(GS.emr_mom, W_.i_))
                .with(function!(GS.emr_mom, W_.i_, component_index));
            sum + component.clone() * component
        });
        (&self.mass_squared + spatial_square).sqrt()
    }

    pub(crate) fn momentum_signature(
        &self,
        graph: &Graph,
        uses_uv_loop_basis: bool,
    ) -> three_dimensional_reps::graph_io::Result<MomentumSignature> {
        self.momentum_signature_in_lmb(&graph.loop_momentum_basis, uses_uv_loop_basis)
    }

    pub(crate) fn momentum_signature_in_lmb(
        &self,
        lmb: &LoopMomentumBasis,
        uses_uv_loop_basis: bool,
    ) -> three_dimensional_reps::graph_io::Result<MomentumSignature> {
        let mut loop_signature = vec![0; lmb.loop_edges.len()];
        let mut external_signature = vec![0; lmb.ext_edges.len()];
        // Provenance is ownership metadata, not part of the linear momentum.
        // Erase it only in this temporary topology-classification copy; the
        // stored denominator and factorized numerator remain fully tagged.
        let physical_momentum = GS.erase_uv_momentum_provenance(&self.momentum);
        accumulate_momentum_signature(
            lmb,
            physical_momentum.as_view(),
            1,
            uses_uv_loop_basis,
            &mut loop_signature,
            &mut external_signature,
        )
        .map_err(|error| {
            GraphIoError::Source(format!(
                "could not extract 4D denominator momentum signature for edge {} from `{}`: {error}",
                usize::from(self.source_edge),
                self.momentum
            ))
        })?;
        Ok(MomentumSignature {
            loop_signature,
            external_signature,
        })
    }

    #[cfg(test)]
    pub(crate) fn depends_on_loop(
        &self,
        graph: &Graph,
        uses_uv_loop_basis: bool,
    ) -> three_dimensional_reps::graph_io::Result<bool> {
        Ok(self
            .momentum_signature(graph, uses_uv_loop_basis)?
            .loop_signature
            .into_iter()
            .any(|coefficient| coefficient != 0))
    }
}

#[cfg(test)]
fn physical_owner_skeleton_from_denominators<'a>(
    graph: &'a Graph,
    denominators: &'a [FourDDenominator],
    uv_edges: impl IntoIterator<Item = EdgeIndex>,
    uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
    coordinate_lmb: Option<&LoopMomentumBasis>,
    sub_lmb_frame: Option<ExactUvSubLmbFrame>,
) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
    GraphThreeDSource::from_exact_denominators_in_coordinates(
        graph,
        denominators,
        uv_edges,
        uv_boundary_hedges,
        coordinate_lmb,
        sub_lmb_frame,
    )?
    .to_three_d_parsed_graph()
}

impl ThreeDGraphSource for GraphThreeDSource<'_> {
    fn to_three_d_parsed_graph(&self) -> three_dimensional_reps::graph_io::Result<ParsedGraph> {
        if self.exact_denominators.is_some() {
            return self.exact_parsed_graph();
        }
        // Keep contraction virtual here: linnet graph contraction mutates/deletes
        // edges, while the CFF expression must still map back to the original
        // EdgeIndex values used by GammaLoop's energy and surface caches.
        let mut parent = (0..self.graph.n_nodes()).collect_vec();
        for (pair, edge_id, edge_data) in self.graph.underlying.iter_edges() {
            if edge_data.data.is_dummy || !self.contracts_edge(edge_id) {
                continue;
            }
            if let HedgePair::Paired { source, sink } = pair {
                union_parent(
                    &mut parent,
                    usize::from(self.graph.node_id(source)),
                    usize::from(self.graph.node_id(sink)),
                );
            }
        }

        let node_ids = self
            .graph
            .underlying
            .iter_nodes()
            .map(|(node_id, _, _)| node_id)
            .sorted()
            .collect::<Vec<_>>();
        let mut root_to_internal = BTreeMap::<usize, usize>::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            if !root_to_internal.contains_key(&root) {
                root_to_internal.insert(root, root_to_internal.len());
            }
        }
        let mut node_to_internal = BTreeMap::new();
        for node_id in &node_ids {
            let root = find_parent(&mut parent, usize::from(*node_id));
            node_to_internal.insert(*node_id, root_to_internal[&root]);
        }

        let loop_names = self
            .outer_loop_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();
        let external_names = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .map(|edge_id| self.graph.underlying[*edge_id].name.value.clone())
            .collect::<Vec<_>>();

        let mut internal_edges = Vec::new();
        let mut external_edges = Vec::new();
        let mut initial_state_cut_edges = Vec::new();
        let mut next_external_id = 10_000_000usize;

        for (pair, edge_index, edge_data) in self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
        {
            if edge_data.data.is_dummy {
                continue;
            }
            let signature = &self.graph.loop_momentum_basis.edge_signatures[edge_index];
            let momentum_signature = MomentumSignature {
                loop_signature: self.outer_loop_signature(edge_index)?,
                external_signature: (&signature.external).into_iter().map(sign_to_i32).collect(),
            };
            let label = edge_data.data.name.value.clone();
            match pair {
                HedgePair::Paired { source, sink } => {
                    if self.contracts_edge(edge_index) {
                        continue;
                    }
                    let local_edge_id = internal_edges.len();
                    let tail = *node_to_internal
                        .get(&self.graph.node_id(source))
                        .ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "missing contracted source node mapping for edge {}",
                                usize::from(edge_index)
                            ))
                        })?;
                    let head =
                        *node_to_internal
                            .get(&self.graph.node_id(sink))
                            .ok_or_else(|| {
                                GraphIoError::Source(format!(
                                    "missing contracted sink node mapping for edge {}",
                                    usize::from(edge_index)
                                ))
                            })?;
                    internal_edges.push(ParsedGraphInternalEdge {
                        edge_id: local_edge_id,
                        tail,
                        head,
                        label,
                        mass_key: Some(edge_data.data.particle.mass_atom().to_canonical_string()),
                        signature: momentum_signature,
                        had_pow: false,
                    });
                    if self.initial_state_cut_edges.contains(&edge_index) {
                        let (external_id, external_sign) = initial_state_cut_external_alias(
                            usize::from(edge_index),
                            &internal_edges[local_edge_id].signature,
                        )?;
                        initial_state_cut_edges.push(ParsedGraphInitialStateCutEdge {
                            edge_id: local_edge_id,
                            external_id,
                            external_sign,
                        });
                    }
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let node = *node_to_internal
                        .get(&self.graph.node_id(hedge))
                        .ok_or_else(|| {
                            GraphIoError::Source(format!(
                                "missing contracted external node mapping for edge {}",
                                usize::from(edge_index)
                            ))
                        })?;
                    let (source, destination) = match flow {
                        Flow::Source => (Some(node), None),
                        Flow::Sink => (None, Some(node)),
                    };
                    external_edges.push(ParsedGraphExternalEdge {
                        edge_id: next_external_id,
                        source,
                        destination,
                        label,
                        external_coefficients: momentum_signature.external_signature,
                    });
                    next_external_id += 1;
                }
                HedgePair::Split { .. } => {
                    return Err(GraphIoError::Source(
                        "split edges are not supported when extracting GammaLoop Graph input"
                            .to_string(),
                    ));
                }
            }
        }

        Ok(ParsedGraph {
            internal_edges,
            external_edges,
            initial_state_cut_edges,
            loop_names,
            external_names,
            node_name_to_internal: root_to_internal
                .into_iter()
                .map(|(root, node)| (format!("n{root}"), node))
                .collect(),
        })
    }

    fn energy_edge_index_map(&self, _parsed: &ParsedGraph) -> Option<EnergyEdgeIndexMap> {
        if self.exact_denominators.is_some() {
            return self.exact_energy_edge_index_map.clone();
        }
        let internal = self
            .graph
            .underlying
            .iter_edges()
            .sorted_by_key(|(_, edge_index, _)| *edge_index)
            .filter_map(|(pair, edge_index, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !self.contracts_edge(edge_index))
                    .then_some(usize::from(edge_index))
            })
            .enumerate()
            .collect::<BTreeMap<_, _>>();

        let external = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .enumerate()
            .map(|(external_id, edge_id)| (external_id, usize::from(*edge_id)))
            .collect::<BTreeMap<_, _>>();

        Some(EnergyEdgeIndexMap {
            internal,
            external,
            orientation_edge_count: self.graph.underlying.n_edges(),
        })
    }
}

/// Canonically relabel an exact topology after its incidences have been
/// inherited from the source graph. This operation only deduplicates already
/// source-built topologies; momentum signatures never determine endpoints.
fn canonicalize_exact_parsed_graph(
    mut parsed: ParsedGraph,
    exact_occurrences: &[ExactParsedOccurrence],
    exact_uv_domains: &[bool],
) -> three_dimensional_reps::graph_io::Result<(ParsedGraph, Vec<usize>)> {
    let exact_count = exact_uv_domains.len();
    if exact_occurrences.len() != exact_count || parsed.internal_edges.len() < exact_count {
        return Err(GraphIoError::Source(format!(
            "exact topology has {} internal edges and {} occurrence records for {exact_count} denominator domains",
            parsed.internal_edges.len(),
            exact_occurrences.len(),
        )));
    }
    let mut exact_domain_by_current = BTreeMap::new();
    for (occurrence, domain) in exact_occurrences.iter().zip(exact_uv_domains) {
        if occurrence.local_edge_id >= parsed.internal_edges.len() {
            return Err(GraphIoError::Source(format!(
                "exact occurrence local edge {} is outside the {}-edge parsed topology",
                occurrence.local_edge_id,
                parsed.internal_edges.len(),
            )));
        }
        if exact_domain_by_current
            .insert(occurrence.local_edge_id, *domain)
            .is_some()
        {
            return Err(GraphIoError::Source(format!(
                "exact topology repeats local occurrence edge {}",
                occurrence.local_edge_id,
            )));
        }
    }
    let exact_current_edges = exact_domain_by_current.keys().copied().collect::<Vec<_>>();

    let source_node_names = parsed.node_name_to_internal.clone();
    let active_nodes = parsed
        .internal_edges
        .iter()
        .flat_map(|edge| [edge.tail, edge.head])
        .chain(
            parsed
                .external_edges
                .iter()
                .flat_map(|edge| edge.source.into_iter().chain(edge.destination)),
        )
        .collect::<BTreeSet<_>>();
    let mut source_node_to_graphica = BTreeMap::new();
    let mut canonicalizer = SymbolicaGraph::<(u8, Vec<i32>), _>::new();
    for node in active_nodes {
        source_node_to_graphica.insert(node, canonicalizer.add_node((0, Vec::new())));
    }

    let initial_cut_by_edge = parsed
        .initial_state_cut_edges
        .iter()
        .map(|edge| (edge.edge_id, (edge.external_id, edge.external_sign)))
        .collect::<BTreeMap<_, _>>();
    let mut canonical_routing_signs = Vec::with_capacity(parsed.internal_edges.len());
    for (position, edge) in parsed.internal_edges.iter().enumerate() {
        let domain = exact_domain_by_current
            .get(&position)
            .copied()
            .map(u8::from)
            .unwrap_or(2);
        // A propagator is even under Q -> -Q. Canonical node labels must
        // therefore use canonical Q together with the correspondingly routed
        // source incidence; external crown incidences remain genuinely
        // directed.
        // Fix Q/-Q in the generation-LMB frame: the first active coordinate
        // is positive. Unlike a purely lexicographic sign, this cannot flip
        // the common loop-coordinate orientation between different powers of
        // the same rational topology.
        let routing_sign = edge
            .signature
            .loop_signature
            .iter()
            .chain(&edge.signature.external_signature)
            .find(|coefficient| **coefficient != 0)
            .map(|coefficient| coefficient.signum())
            .unwrap_or(1);
        let canonical_signature = if routing_sign == 1 {
            edge.signature.clone()
        } else {
            edge.signature.negated()
        };
        canonical_routing_signs.push(routing_sign);
        let (tail, head) = if routing_sign == 1 {
            (edge.tail, edge.head)
        } else {
            (edge.head, edge.tail)
        };
        let canonical_cut = initial_cut_by_edge
            .get(&edge.edge_id)
            .map(|(external_id, external_sign)| (*external_id, routing_sign * external_sign));
        canonicalizer
            .add_edge(
                source_node_to_graphica[&tail],
                source_node_to_graphica[&head],
                true,
                (
                    domain,
                    canonical_signature,
                    edge.mass_key.clone(),
                    canonical_cut,
                ),
            )
            .map_err(|error| GraphIoError::Source(error.to_string()))?;
    }
    for edge in &parsed.external_edges {
        let external_node = canonicalizer.add_node((1, edge.external_coefficients.clone()));
        match (edge.source, edge.destination) {
            (Some(source), None) => {
                canonicalizer
                    .add_edge(
                        source_node_to_graphica[&source],
                        external_node,
                        true,
                        (
                            3,
                            MomentumSignature {
                                loop_signature: Vec::new(),
                                external_signature: edge.external_coefficients.clone(),
                            },
                            None,
                            None,
                        ),
                    )
                    .map_err(|error| GraphIoError::Source(error.to_string()))?;
            }
            (None, Some(destination)) => {
                canonicalizer
                    .add_edge(
                        external_node,
                        source_node_to_graphica[&destination],
                        true,
                        (
                            3,
                            MomentumSignature {
                                loop_signature: Vec::new(),
                                external_signature: edge.external_coefficients.clone(),
                            },
                            None,
                            None,
                        ),
                    )
                    .map_err(|error| GraphIoError::Source(error.to_string()))?;
            }
            _ => {
                return Err(GraphIoError::Source(format!(
                    "exact topology external edge {} must have exactly one internal endpoint",
                    edge.edge_id,
                )));
            }
        }
    }

    let canonical = canonicalizer.canonize();
    let mut source_nodes_by_canonical_position = source_node_to_graphica
        .iter()
        .map(|(source_node, graphica_node)| (*source_node, canonical.vertex_map[*graphica_node]))
        .collect::<Vec<_>>();
    source_nodes_by_canonical_position.sort_by_key(|(_, canonical_node)| *canonical_node);
    let source_to_canonical = source_nodes_by_canonical_position
        .into_iter()
        .enumerate()
        .map(|(canonical_node, (source_node, _))| (source_node, canonical_node))
        .collect::<BTreeMap<_, _>>();
    let mut reversed_internal_edges = BTreeSet::new();
    for (edge, routing_sign) in parsed
        .internal_edges
        .iter_mut()
        .zip(canonical_routing_signs)
    {
        if routing_sign == -1 {
            reversed_internal_edges.insert(edge.edge_id);
            std::mem::swap(&mut edge.tail, &mut edge.head);
            edge.signature = edge.signature.negated();
        }
        edge.tail = source_to_canonical[&edge.tail];
        edge.head = source_to_canonical[&edge.head];
    }
    // Graphica canonizes the source-built rational skeleton only after each
    // even propagator has been put in its canonical Q routing. Swapping both
    // incidence and momentum sign leaves the routed edge unchanged, keeps
    // every serial dotted chain balanced, and makes owner IDs irrelevant to
    // the algebraic topology.
    for cut_edge in &mut parsed.initial_state_cut_edges {
        if reversed_internal_edges.contains(&cut_edge.edge_id) {
            cut_edge.external_sign = cut_edge.external_sign.checked_neg().ok_or_else(|| {
                GraphIoError::Source(
                    "exact topology cut routing sign cannot be negated".to_string(),
                )
            })?;
        }
    }
    for edge in &mut parsed.external_edges {
        edge.source = edge.source.map(|node| source_to_canonical[&node]);
        edge.destination = edge.destination.map(|node| source_to_canonical[&node]);
    }
    parsed.external_edges.sort_by_key(|edge| {
        (
            edge.source,
            edge.destination,
            edge.external_coefficients.clone(),
            edge.edge_id,
        )
    });
    for (edge_id, edge) in parsed.external_edges.iter_mut().enumerate() {
        edge.edge_id = edge_id;
        edge.label = format!("__gammaloop_exact_external_{edge_id}");
    }

    let mut canonical_to_current = exact_current_edges;
    canonical_to_current.sort_by_key(|current| {
        let edge = &parsed.internal_edges[*current];
        (
            edge.tail,
            edge.head,
            exact_domain_by_current[current],
            edge.signature.clone(),
            edge.mass_key.clone(),
            edge.had_pow,
            *current,
        )
    });
    let cut_by_current = parsed
        .initial_state_cut_edges
        .iter()
        .map(|cut| (cut.edge_id, (cut.external_id, cut.external_sign)))
        .collect::<BTreeMap<_, _>>();
    let mut cut_current_edges = (0..parsed.internal_edges.len())
        .filter(|current| !exact_domain_by_current.contains_key(current))
        .collect::<Vec<_>>();
    cut_current_edges.sort_by_key(|current| {
        let edge = &parsed.internal_edges[*current];
        (
            edge.tail,
            edge.head,
            edge.signature.clone(),
            edge.mass_key.clone(),
            cut_by_current.get(current).copied(),
            *current,
        )
    });
    canonical_to_current.extend(cut_current_edges);
    let current_to_canonical = canonical_to_current
        .iter()
        .enumerate()
        .map(|(canonical, current)| (*current, canonical))
        .collect::<BTreeMap<_, _>>();
    let internal_edges = canonical_to_current
        .iter()
        .enumerate()
        .map(|(canonical, current)| {
            let mut edge = parsed.internal_edges[*current].clone();
            edge.edge_id = canonical;
            edge.label = if canonical < exact_count {
                format!("__gammaloop_exact_edge_{canonical}")
            } else {
                format!("__gammaloop_exact_cut_{}", canonical - exact_count)
            };
            edge
        })
        .collect::<Vec<_>>();
    for cut_edge in &mut parsed.initial_state_cut_edges {
        cut_edge.edge_id = current_to_canonical[&cut_edge.edge_id];
    }
    parsed.internal_edges = internal_edges;
    parsed.loop_names = (0..parsed.loop_names.len())
        .map(|index| format!("ell{index}"))
        .collect();
    parsed.node_name_to_internal = source_node_names
        .into_iter()
        .filter_map(|(name, source_node)| {
            source_to_canonical
                .get(&source_node)
                .copied()
                .map(|canonical_node| (name, canonical_node))
        })
        .collect();

    let validation = validate_parsed_graph(&parsed);
    if !validation.vertex_balance_violations.is_empty() {
        return Err(GraphIoError::Source(format!(
            "source-backed exact topology violates loop-momentum conservation at vertices {:?}",
            validation.vertex_balance_violations,
        )));
    }
    // The source builder has already checked external balance relative to the
    // original graph, including dependent production momenta and cut aliases.
    // Canonical relabeling preserves those incidences, so reinterpreting the
    // residual here as an absolute zero condition would reject valid sources.
    Ok((parsed, canonical_to_current))
}

fn find_parent(parent: &mut [usize], node: usize) -> usize {
    let parent_node = parent[node];
    if parent_node == node {
        node
    } else {
        let root = find_parent(parent, parent_node);
        parent[node] = root;
        root
    }
}

fn union_parent(parent: &mut [usize], left: usize, right: usize) {
    let left_root = find_parent(parent, left);
    let right_root = find_parent(parent, right);
    if left_root != right_root {
        parent[right_root] = left_root;
    }
}

fn serial_power_chain_incidences(
    tail: usize,
    head: usize,
    source_routing_sign: i32,
    power: usize,
    next_node: &mut usize,
) -> (Vec<(usize, usize)>, Vec<usize>) {
    if power == 0 {
        return (Vec::new(), Vec::new());
    }
    let auxiliary_nodes = (*next_node..*next_node + power.saturating_sub(1)).collect::<Vec<_>>();
    *next_node += auxiliary_nodes.len();
    let chain_nodes = std::iter::once(tail)
        .chain(auxiliary_nodes.iter().copied())
        .chain(std::iter::once(head))
        .collect::<Vec<_>>();
    let incidences = chain_nodes
        .windows(2)
        .map(|endpoints| {
            if source_routing_sign == 1 {
                (endpoints[0], endpoints[1])
            } else {
                (endpoints[1], endpoints[0])
            }
        })
        .collect();
    (incidences, auxiliary_nodes)
}

fn select_exact_loop_basis(
    rows: &[Vec<i32>],
    loop_count: usize,
) -> three_dimensional_reps::graph_io::Result<(Vec<Vec<i32>>, Vec<usize>)> {
    if rows.iter().any(|row| row.len() != loop_count) {
        return Err(GraphIoError::Source(
            "exact 4D denominator loop signatures have inconsistent dimensions".to_string(),
        ));
    }
    let active_rank = exact_integer_rank(rows);
    if active_rank == 0 {
        return Ok((Vec::new(), Vec::new()));
    }

    // Prefer the production carriers, but only accept denominator rows whose
    // active minor keeps all rewritten momentum signatures integral.
    let mut candidates = rows
        .iter()
        .filter(|row| row.iter().any(|coefficient| *coefficient != 0))
        .map(|row| {
            let mut canonical = row.clone();
            if canonical
                .iter()
                .find(|coefficient| **coefficient != 0)
                .is_some_and(|coefficient| *coefficient < 0)
            {
                for coefficient in &mut canonical {
                    *coefficient = coefficient.checked_neg().ok_or_else(|| {
                        GraphIoError::Source(
                            "exact 4D denominator loop coefficient cannot be canonicalized"
                                .to_string(),
                        )
                    })?;
                }
            }
            Ok(canonical)
        })
        .collect::<three_dimensional_reps::graph_io::Result<BTreeSet<_>>>()?
        .into_iter()
        .collect::<Vec<_>>();
    let unit_axis = |row: &[i32]| {
        row.iter().enumerate().find_map(|(axis, coefficient)| {
            (*coefficient == 1
                && row
                    .iter()
                    .enumerate()
                    .all(|(column, value)| column == axis || *value == 0))
            .then_some(axis)
        })
    };
    candidates.sort_by(|left, right| match (unit_axis(left), unit_axis(right)) {
        (Some(left), Some(right)) => left.cmp(&right),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => left.cmp(right),
    });

    for selected in candidates.iter().cloned().combinations(active_rank) {
        for columns in (0..loop_count).combinations(active_rank) {
            let matrix = selected
                .iter()
                .map(|row| {
                    columns
                        .iter()
                        .map(|column| Rational::from(row[*column]))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let inverse_is_integral = (0..active_rank).all(|rhs_index| {
                let rhs = (0..active_rank)
                    .map(|index| Rational::from(i32::from(index == rhs_index)))
                    .collect::<Vec<_>>();
                solve_rational_system(matrix.clone(), rhs).is_some_and(|solution| {
                    solution
                        .iter()
                        .all(|coordinate| coordinate.denominator_ref().to_i64() == Some(1))
                })
            });
            if inverse_is_integral {
                return Ok((selected, columns));
            }
        }
    }

    Err(GraphIoError::Source(format!(
        "exact 4D source rank {active_rank} has no unimodular denominator-carrier basis in the {loop_count}-dimensional production CFF coordinates"
    )))
}

fn accumulate_momentum_signature(
    lmb: &LoopMomentumBasis,
    view: AtomView<'_>,
    coefficient: i32,
    use_uv_loop_basis: bool,
    loop_signature: &mut [i32],
    external_signature: &mut [i32],
) -> std::result::Result<(), String> {
    match view {
        AtomView::Add(add) => {
            for term in add.iter() {
                accumulate_momentum_signature(
                    lmb,
                    term,
                    coefficient,
                    use_uv_loop_basis,
                    loop_signature,
                    external_signature,
                )?;
            }
        }
        AtomView::Mul(mul) => {
            let mut scalar = coefficient;
            let mut momentum_factor = None;
            for factor in mul.iter() {
                if let Ok(integer) = i64::try_from(factor) {
                    scalar = scalar
                        .checked_mul(i32::try_from(integer).map_err(|error| error.to_string())?)
                        .ok_or_else(|| {
                            "integer coefficient overflow in 4D denominator momentum".to_string()
                        })?;
                } else if momentum_factor.replace(factor).is_some() {
                    return Err(format!(
                        "expected a linear momentum expression, found product `{}`",
                        view.to_owned()
                    ));
                }
            }
            if let Some(momentum_factor) = momentum_factor {
                accumulate_momentum_signature(
                    lmb,
                    momentum_factor,
                    scalar,
                    use_uv_loop_basis,
                    loop_signature,
                    external_signature,
                )?;
            } else if scalar != 0 {
                return Err(format!(
                    "expected a momentum factor, found scalar `{}`",
                    view.to_owned()
                ));
            }
        }
        AtomView::Fun(function) => {
            if function.get_symbol() != GS.emr_mom || function.get_nargs() != 1 {
                return Err(format!(
                    "expected an EMR momentum, found `{}`",
                    view.to_owned()
                ));
            }
            let edge = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
                format!(
                    "EMR momentum has non-integer edge id `{}`",
                    function.get(0).to_owned()
                )
            })?);
            if use_uv_loop_basis
                && let Some(loop_index) = lmb
                    .loop_edges
                    .iter()
                    .position(|loop_edge| *loop_edge == edge)
            {
                loop_signature[loop_index] += coefficient;
                return Ok(());
            }
            let signature = lmb
                .edge_signatures
                .get(edge)
                .ok_or_else(|| format!("edge {edge} is not present in the loop-momentum basis"))?;
            for (loop_index, sign) in signature.internal.iter_enumerated() {
                loop_signature[usize::from(loop_index)] += coefficient * sign_to_i32(*sign);
            }
            for (external_index, sign) in signature.external.iter_enumerated() {
                external_signature[usize::from(external_index)] += coefficient * sign_to_i32(*sign);
            }
        }
        AtomView::Num(_) if view.to_owned().is_zero() => {}
        AtomView::Num(_) => {
            return Err(format!(
                "expected a momentum expression, found scalar `{}`",
                view.to_owned()
            ));
        }
        _ => {
            return Err(format!(
                "expected a linear momentum expression, found `{}`",
                view.to_owned()
            ));
        }
    }
    Ok(())
}

fn exact_integer_rank(rows: &[Vec<i32>]) -> usize {
    if rows.first().is_none_or(Vec::is_empty) {
        return 0;
    }
    rank_i64(
        &rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|coefficient| i64::from(*coefficient))
                    .collect()
            })
            .collect::<Vec<Vec<_>>>(),
    )
}

fn sign_to_i32(sign: SignOrZero) -> i32 {
    match sign {
        SignOrZero::Minus => -1,
        SignOrZero::Zero => 0,
        SignOrZero::Plus => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::expression::GammaLoopOrientationExpression,
        cff::surface::GammaLoopSurfaceCache,
        dot,
        graph::{feynman_graph::FeynmanGraph, parse::IntoGraph},
        initialisation::test_initialise,
        momentum::sample::LoopIndex,
        numerator::energy_degree::{EnergyPowerAnalyzer, EquivalentEnergyCandidates},
    };
    use color_eyre::Result;
    use linnet::half_edge::subgraph::SubSetOps;

    #[test]
    fn contracted_source_preserves_emr_ids_without_denominator_aliases() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            d -> e [id=3 lmb_id=1]
            e -> f [id=4]
            f -> d [id=5]
            c -> d [id=6]
        })?;
        let contracted_edges = [EdgeIndex::from(0), EdgeIndex::from(1), EdgeIndex::from(2)];
        let source = GraphThreeDSource::new(&graph, &contracted_edges)?;
        assert_eq!(source.inner_loop_count, 1);
        assert_eq!(source.outer_loop_edges.len(), 1);

        let inner_carrier = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .find(|edge_id| contracted_edges.contains(edge_id))
            .expect("contracted loop should retain a parent-basis carrier");
        assert!(
            source
                .reconstructible_outer_loop_coordinates(inner_carrier)
                .is_none()
        );
        let parsed = source.to_three_d_parsed_graph()?;
        let edge_map = source
            .energy_edge_index_map(&parsed)
            .expect("GammaLoop sources provide an edge-index map");
        assert!(
            !edge_map
                .internal
                .values()
                .any(|edge| *edge == usize::from(inner_carrier))
        );
        let outer_carrier = source.outer_loop_edges[0];
        assert!(
            edge_map
                .internal
                .values()
                .any(|edge| *edge == usize::from(outer_carrier))
        );

        let contracted_tree_edge = EdgeIndex::from(6);
        let tree_source = GraphThreeDSource::new(&graph, &[contracted_tree_edge])?;
        assert_eq!(tree_source.inner_loop_count, 0);
        assert_eq!(
            tree_source.reconstructible_outer_loop_coordinates(contracted_tree_edge),
            tree_source.edge_loop_coordinates(contracted_tree_edge),
        );
        let parsed = tree_source.to_three_d_parsed_graph()?;
        let edge_map = tree_source
            .energy_edge_index_map(&parsed)
            .expect("GammaLoop sources provide an edge-index map");
        assert!(
            !edge_map
                .internal
                .values()
                .any(|edge| *edge == usize::from(contracted_tree_edge)),
            "a shrunken edge must not be aliased to a surviving denominator"
        );
        Ok(())
    }

    #[test]
    fn exact_source_serializes_dotted_same_edge_occurrences() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            c -> d [id=3]
        })?;
        let edge = EdgeIndex::from(0);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .finish();
        let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
        let first = FourDDenominator {
            source_edge: edge,
            momentum: momentum.clone(),
            mass_squared: mass_squared.clone(),
            full_expr: Atom::var(symbolica::symbol!("three_d_source_test::first")),
        };
        let second = FourDDenominator {
            source_edge: edge,
            momentum,
            mass_squared,
            full_expr: Atom::var(symbolica::symbol!("three_d_source_test::second")),
        };
        let denominators = [first.clone(), first, second];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let parsed = source.to_three_d_parsed_graph()?;

        assert_eq!(parsed.internal_edges.len(), 3);
        assert_eq!(
            parsed
                .internal_edges
                .iter()
                .map(|edge| edge.edge_id)
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(
            parsed
                .node_name_to_internal
                .keys()
                .filter(|name| name.starts_with("__gammaloop_exact_power_"))
                .count(),
            2
        );
        let [first, second, third] = parsed.internal_edges.as_slice() else {
            panic!("the cubic exact propagator must retain three occurrences");
        };
        let mut in_degree = BTreeMap::<usize, usize>::new();
        let mut out_degree = BTreeMap::<usize, usize>::new();
        for edge in [first, second, third] {
            assert_ne!(edge.tail, edge.head);
            *out_degree.entry(edge.tail).or_default() += 1;
            *in_degree.entry(edge.head).or_default() += 1;
        }
        assert!(in_degree.values().all(|degree| *degree == 1));
        assert!(out_degree.values().all(|degree| *degree == 1));
        assert_eq!(
            parsed
                .internal_edges
                .iter()
                .flat_map(|edge| [edge.tail, edge.head])
                .collect::<BTreeSet<_>>()
                .len(),
            3,
            "a cubic self-loop power must be a three-edge cycle"
        );
        assert!(
            three_dimensional_reps::validate_parsed_graph(&parsed)
                .vertex_balance_violations
                .is_empty()
        );
        let offset = graph.underlying.n_edges();
        assert_eq!(
            source.energy_edge_index_map(&parsed).unwrap().internal,
            BTreeMap::from([(0, offset), (1, offset + 1), (2, offset + 2)])
        );
        assert_eq!(
            source.physical_energy_edge_index_map().unwrap().internal,
            BTreeMap::from([
                (offset, usize::from(edge)),
                (offset + 1, usize::from(edge)),
                (offset + 2, usize::from(edge)),
            ])
        );
        Ok(())
    }

    #[test]
    fn exact_source_owner_relabeling_preserves_residue_rank_and_cut_provenance() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_distinct_owners {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let distinct_owners = [EdgeIndex(0), EdgeIndex(1)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        // Relabel the second algebraically identical denominator with the first
        // physical owner. This changes provenance and the source-minor
        // realization, but not the rational function to be integrated.
        let mut relabeled_owner = distinct_owners.clone();
        relabeled_owner[1].source_edge = EdgeIndex(0);

        let options = graph.denominator_only_cff_3d_expression_options();
        let cutset = crate::graph::cuts::CutSet::empty(graph.n_hedges());
        let mut build = |denominators: &[FourDDenominator], numerator: &Atom| {
            let (parsed, physical_map) = {
                let source = GraphThreeDSource::from_exact_denominators(&graph, denominators)?;
                (
                    source.to_three_d_parsed_graph()?,
                    source
                        .physical_energy_edge_index_map()
                        .expect("exact source retains physical-energy provenance")
                        .internal,
                )
            };
            let (cff, _) =
                graph.cff_from_4d_denominators(denominators, &cutset, &options, numerator)?;
            let cut_support = cff
                .terms
                .values()
                .flat_map(|term| &term.orientations)
                .flat_map(|orientation| &orientation.orientation.variants)
                .flat_map(|variant| &variant.denominator_edges)
                .copied()
                .collect::<BTreeSet<_>>();
            let expression = cff
                .terms
                .values()
                .flat_map(|term| {
                    term.orientations.iter().map(|orientation| {
                        Ok(orientation.expression.clone()
                            * term.map_exact_source_numerator(&orientation.orientation)?)
                    })
                })
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .fold(Atom::Zero, |sum, orientation| sum + orientation)
                * Atom::num(cff.production_prefactor_factor());
            Ok::<_, color_eyre::Report>((parsed, physical_map, cut_support, expression))
        };

        let (distinct_parsed, distinct_map, distinct_support, distinct_expression) =
            build(&distinct_owners, &Atom::one())?;
        let (relabeled_parsed, relabeled_map, relabeled_support, relabeled_expression) =
            build(&relabeled_owner, &Atom::one())?;

        for parsed in [&distinct_parsed, &relabeled_parsed] {
            assert_eq!(parsed.internal_edges.len(), 2);
            assert!(three_dimensional_reps::validate_parsed_graph(parsed).ok);
            let nodes = parsed
                .internal_edges
                .iter()
                .flat_map(|edge| [edge.tail, edge.head])
                .collect::<BTreeSet<_>>();
            assert_eq!(
                parsed.internal_edges.len() + 1 - nodes.len(),
                1,
                "an algebraic squared propagator must retain one topological loop"
            );
            assert_eq!(parsed.loop_names.len(), 1);
        }

        // Owner aliases only rename equal on-shell energies. Once that neutral
        // alias is identified, both source realizations must give the same CFF
        // residue including the production prefactor convention.
        let canonicalize_owner_alias = |mut expression: Atom| {
            expression = expression
                .replace(GS.ose(EdgeIndex(1)))
                .with(GS.ose(EdgeIndex(0)));
            for component in 0..=3 {
                expression = expression
                    .replace(GS.emr_mom(EdgeIndex(1), GS.cind(component)))
                    .with(GS.emr_mom(EdgeIndex(0), GS.cind(component)));
            }
            expression
        };
        let residue_difference = canonicalize_owner_alias(distinct_expression)
            - canonicalize_owner_alias(relabeled_expression);
        for [qx, qy, qz, energy] in [
            [Atom::Zero, Atom::Zero, Atom::Zero, Atom::one()],
            [
                Atom::num(3) / Atom::num(4),
                Atom::Zero,
                Atom::Zero,
                Atom::num(5) / Atom::num(4),
            ],
            [
                Atom::num(4) / Atom::num(3),
                Atom::Zero,
                Atom::Zero,
                Atom::num(5) / Atom::num(3),
            ],
        ] {
            let fixed_difference = residue_difference
                .replace(GS.emr_mom(EdgeIndex(0), GS.cind(1)))
                .with(qx)
                .replace(GS.emr_mom(EdgeIndex(0), GS.cind(2)))
                .with(qy)
                .replace(GS.emr_mom(EdgeIndex(0), GS.cind(3)))
                .with(qz)
                .replace(GS.ose(EdgeIndex(0)))
                .with(energy)
                .together();
            assert_eq!(
                fixed_difference,
                Atom::Zero,
                "the exact CFF residue must not depend on which physical owner labels equal denominators"
            );
        }

        let distinct_numerator = GS.emr_mom(EdgeIndex(1), GS.cind(0)).pow(2);
        let relabeled_numerator = GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2);
        let (_, _, _, distinct_quadratic) = build(&distinct_owners, &distinct_numerator)?;
        let (_, _, _, relabeled_quadratic) = build(&relabeled_owner, &relabeled_numerator)?;
        let quadratic_difference = canonicalize_owner_alias(distinct_quadratic)
            - canonicalize_owner_alias(relabeled_quadratic);
        for [qx, qy, qz, energy] in [
            [Atom::Zero, Atom::Zero, Atom::Zero, Atom::one()],
            [
                Atom::num(3) / Atom::num(4),
                Atom::Zero,
                Atom::Zero,
                Atom::num(5) / Atom::num(4),
            ],
            [
                Atom::num(4) / Atom::num(3),
                Atom::Zero,
                Atom::Zero,
                Atom::num(5) / Atom::num(3),
            ],
        ] {
            let fixed_quadratic_difference = quadratic_difference
                .clone()
                .replace(GS.emr_mom(EdgeIndex(0), GS.cind(1)))
                .with(qx)
                .replace(GS.emr_mom(EdgeIndex(0), GS.cind(2)))
                .with(qy)
                .replace(GS.emr_mom(EdgeIndex(0), GS.cind(3)))
                .with(qz)
                .replace(GS.ose(EdgeIndex(0)))
                .with(energy)
                .together();
            assert!(
                fixed_quadratic_difference.is_zero(),
                "equal rational denominators must preserve an attached quadratic energy numerator across a physical-owner relabel: {fixed_quadratic_difference}; symbolic difference: {quadratic_difference}",
            );
        }

        assert_eq!(
            distinct_map.values().copied().collect::<BTreeSet<_>>(),
            BTreeSet::from([0, 1]),
            "distinct-owner physical-energy provenance must remain distinct"
        );
        assert_eq!(
            relabeled_map.values().copied().collect::<BTreeSet<_>>(),
            BTreeSet::from([0]),
            "owner relabeling must remain visible in physical-energy provenance"
        );
        assert_eq!(
            distinct_support,
            BTreeSet::from([EdgeIndex(0), EdgeIndex(1)]),
            "cut support must retain both distinct physical owners"
        );
        assert_eq!(
            relabeled_support,
            BTreeSet::from([EdgeIndex(0)]),
            "cut support must retain the relabeled physical owner without changing the residue"
        );
        Ok(())
    }

    #[test]
    fn exact_source_keeps_source_instantiated_domains_and_masses_separate() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_same_owner_domains {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=3]
            c -> a [id=4]
            d -> e [id=1 lmb_id=1]
            e -> f [id=5]
            f -> d [id=6]
            g -> h [id=2 lmb_id=2]
            h -> i [id=7]
            i -> g [id=8]
        })?;
        let cograph_edge = EdgeIndex(0);
        let other_cograph_edge = EdgeIndex(1);
        let uv_edge = EdgeIndex(2);
        let cograph_mass = graph.underlying[cograph_edge].particle.mass_atom().pow(2);
        let other_cograph_mass = Atom::num(2);
        let uv_mass = Atom::var(GS.m_uv_expansion).pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: cograph_edge,
                momentum: FunctionBuilder::new(GS.emr_mom).add_arg(0).finish(),
                mass_squared: cograph_mass.clone(),
                full_expr: Atom::one(),
            },
            FourDDenominator {
                source_edge: other_cograph_edge,
                momentum: FunctionBuilder::new(GS.emr_mom).add_arg(1).finish(),
                mass_squared: other_cograph_mass.clone(),
                full_expr: Atom::one(),
            },
            FourDDenominator {
                source_edge: uv_edge,
                momentum: FunctionBuilder::new(GS.emr_mom).add_arg(2).finish(),
                mass_squared: uv_mass.clone(),
                full_expr: Atom::one(),
            },
        ];
        let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &denominators,
            [uv_edge, EdgeIndex(7), EdgeIndex(8)],
        )?;
        let parsed = source.to_three_d_parsed_graph()?;
        let edge_for_mass = |mass: &Atom| {
            let mass_key = mass.to_canonical_string();
            parsed
                .internal_edges
                .iter()
                .find(|edge| edge.mass_key.as_deref() == Some(mass_key.as_str()))
                .expect("each exact mass retains one occurrence")
        };
        let cograph = edge_for_mass(&cograph_mass);
        let other_cograph = edge_for_mass(&other_cograph_mass);
        let uv = edge_for_mass(&uv_mass);

        assert_eq!(cograph.tail, cograph.head);
        assert_eq!(other_cograph.tail, other_cograph.head);
        assert_eq!(uv.tail, uv.head);
        assert_eq!(
            [cograph.tail, other_cograph.tail, uv.tail]
                .into_iter()
                .collect::<BTreeSet<_>>()
                .len(),
            3,
            "distinct source components and UV/cograph domains must remain separate after contraction"
        );
        assert!(
            parsed
                .node_name_to_internal
                .keys()
                .all(|name| !name.starts_with("__gammaloop_exact_power_")),
            "different masses and topology domains must not form a power chain"
        );
        let offset = graph.underlying.n_edges();
        let expected_physical_map = source
            .exact_local_to_original_occurrence
            .iter()
            .enumerate()
            .map(|(local, original)| {
                (
                    offset + local,
                    usize::from(denominators[*original].source_edge),
                )
            })
            .collect::<BTreeMap<_, _>>();
        assert_eq!(
            source.physical_energy_edge_index_map().unwrap().internal,
            expected_physical_map
        );
        Ok(())
    }

    #[test]
    fn exact_source_normalizes_opposite_spelling_inside_one_power_chain() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_opposite_same_owner {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let momentum = FunctionBuilder::new(GS.emr_mom).add_arg(0).finish();
        let repeated_mass = graph.underlying[EdgeIndex(0)].particle.mass_atom().pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: EdgeIndex(0),
                momentum: momentum.clone(),
                mass_squared: repeated_mass.clone(),
                full_expr: Atom::var(symbolica::symbol!("three_d_source_test::positive_routing")),
            },
            FourDDenominator {
                source_edge: EdgeIndex(0),
                momentum: -momentum.clone(),
                mass_squared: repeated_mass,
                full_expr: Atom::var(symbolica::symbol!("three_d_source_test::negative_routing")),
            },
            // Keep the physical endpoints distinct after exact-source
            // contraction and balance the effective opposite routing.
            FourDDenominator {
                source_edge: EdgeIndex(1),
                momentum,
                mass_squared: Atom::num(2),
                full_expr: Atom::one(),
            },
        ];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let parsed = source.to_three_d_parsed_graph()?;
        let owner_occurrences = source
            .exact_local_to_original_occurrence
            .iter()
            .enumerate()
            .filter_map(|(occurrence, original)| {
                (denominators[*original].source_edge == EdgeIndex(0)).then_some(occurrence)
            })
            .collect::<Vec<_>>();
        let [first_id, second_id] = owner_occurrences.as_slice() else {
            panic!("the repeated owner must retain two occurrence-local IDs");
        };
        let first = &parsed.internal_edges[*first_id];
        let second = &parsed.internal_edges[*second_id];
        assert_eq!(first.signature, second.signature);
        assert!(
            first.head == second.tail || second.head == first.tail,
            "Q and -Q copies of one even denominator must form one consistently routed power chain"
        );
        assert_eq!(
            parsed
                .node_name_to_internal
                .keys()
                .filter(|name| name.starts_with("__gammaloop_exact_power_"))
                .count(),
            1
        );
        assert!(
            three_dimensional_reps::validate_parsed_graph(&parsed)
                .vertex_balance_violations
                .is_empty(),
            "canonical denominator spelling with source-owned routing must preserve momentum balance"
        );
        let offset = graph.underlying.n_edges();
        assert_eq!(
            source.energy_edge_index_map(&parsed).unwrap().internal,
            BTreeMap::from([(0, offset), (1, offset + 1), (2, offset + 2),])
        );
        let expected_physical_map = source
            .exact_local_to_original_occurrence
            .iter()
            .enumerate()
            .map(|(local, original)| {
                (
                    offset + local,
                    usize::from(denominators[*original].source_edge),
                )
            })
            .collect::<BTreeMap<_, _>>();
        assert_eq!(
            source.physical_energy_edge_index_map().unwrap().internal,
            expected_physical_map
        );

        // D(Q) and D(-Q) are the same even denominator. Their canonical CFF
        // energy is therefore identical, while recovering the odd physical
        // numerator Q^0 must compose the literal and denominator->canonical
        // signs rather than inherit either sign separately.
        let mapper = source
            .exact_source_energy_mapper()
            .expect("exact source has an energy mapper");
        let mapped_energy = Atom::var(symbolica::symbol!(
            "three_d_source_test::canonical_odd_energy"
        ));
        let mut loop_energy = LinearEnergyExpr::zero();
        loop_energy.constant = mapped_energy.clone();
        let numerator = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let mut physical_energies = Vec::new();
        for original in [0, 1] {
            let local = source
                .exact_local_to_original_occurrence
                .iter()
                .position(|candidate| *candidate == original)
                .expect("both source-edge occurrences survive canonicalization");
            let occurrence = graph.underlying.n_edges() + local;
            let mut edge_energy_map =
                vec![LinearEnergyExpr::zero(); graph.underlying.n_edges() + denominators.len()];
            edge_energy_map[occurrence].constant = mapped_energy.clone();
            physical_energies.push(mapper.map_numerator_factor(
                std::slice::from_ref(&loop_energy),
                &edge_energy_map,
                &numerator,
                &BTreeMap::from([(EdgeIndex(0), occurrence)]),
            )?);
        }
        assert_eq!(
            physical_energies[0], physical_energies[1],
            "literal +Q and -Q denominator spellings must recover the same physical odd numerator energy",
        );
        assert!(
            physical_energies[0] == mapped_energy || physical_energies[0] == -mapped_energy.clone(),
            "canonical Q routing may choose either overall frame but must preserve the odd numerator magnitude",
        );
        Ok(())
    }

    #[test]
    fn exact_massless_surface_uses_normalized_edge_mass() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_massless_surface {
            edge [num=1 mass="UFO::ZERO"]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
        })?;
        let denominators = [EdgeIndex(0), EdgeIndex(1)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        assert_eq!(denominators[0].mass_squared, Atom::Zero);
        assert_ne!(
            denominators[0].mass_squared,
            graph.underlying[EdgeIndex(0)].particle.mass_atom().pow(2)
        );

        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let physical_edges = source
            .physical_energy_edge_index_map()
            .expect("exact source has a physical energy-edge projection");
        let occurrence = physical_edges
            .internal
            .iter()
            .find_map(|(occurrence, physical)| (*physical == 0).then_some(*occurrence))
            .expect("exact source contains the first physical edge");
        let surface = LinearSurface {
            kind: three_dimensional_reps::LinearSurfaceKind::Esurface,
            expression: three_dimensional_reps::LinearEnergyExpr::ose(EdgeIndex(occurrence), 1),
            origin: three_dimensional_reps::surface::SurfaceOrigin::Physical,
            numerator_only: false,
        };
        let physical = source
            .physical_linear_surface(&surface)
            .expect("massless original denominator remains a physical surface");

        assert_eq!(
            physical.expression,
            three_dimensional_reps::LinearEnergyExpr::ose(EdgeIndex(0), 1)
        );
        Ok(())
    }

    #[test]
    fn unexpanded_exact_source_is_the_production_source() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph exact_unexpanded_identity {
                num = 1
                edge [particle="scalar_1" num=1]
                node [num=1]

                ext0 [style=invis is_cut=0]
                v2 -> ext0 [id=0]
                ext0 -> v3
                v0 -> v1 [id=1 lmb_id=0]
                v0 -> v1 [id=2]
                v3 -> v0 [id=3 lmb_id=1]
                v1 -> v2 [id=4]
                v2 -> v3 [id=5]
            },
            "scalars"
        )?;
        let cut_edges = graph
            .get_edges_in_initial_state_cut()
            .into_iter()
            .collect::<BTreeSet<_>>();
        let denominators = graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, source_edge, edge)| {
                (pair.is_paired() && !edge.data.is_dummy && !cut_edges.contains(&source_edge)).then(
                    || FourDDenominator {
                        source_edge,
                        momentum: FunctionBuilder::new(GS.emr_mom)
                            .add_arg(usize::from(source_edge))
                            .finish(),
                        mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
                        full_expr: Atom::one(),
                    },
                )
            })
            .collect::<Vec<_>>();
        let production = GraphThreeDSource::new(&graph, &[])?;
        let exact = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let production_parsed = production.to_three_d_parsed_graph()?;
        let exact_parsed = exact.to_three_d_parsed_graph()?;

        assert_eq!(exact_parsed, production_parsed);
        assert_eq!(
            exact.energy_edge_index_map(&exact_parsed),
            production.energy_edge_index_map(&production_parsed),
        );
        let physical_edges: BTreeMap<_, _> = production
            .energy_edge_index_map(&production_parsed)
            .unwrap()
            .internal
            .into_values()
            .map(|edge_id| (edge_id, edge_id))
            .collect();
        assert_eq!(
            exact.physical_energy_edge_index_map().unwrap().internal,
            physical_edges,
        );
        let mapper = exact.exact_source_energy_mapper().unwrap();
        for occurrence in &exact.exact_occurrences {
            let owner = denominators[occurrence.original_occurrence].source_edge;
            assert_eq!(occurrence.energy_edge_id, usize::from(owner));
            assert_eq!(
                mapper.source_edge_occurrences[&owner],
                vec![ExactSourceOwnerOccurrence {
                    energy_edge_id: usize::from(owner),
                    raw_momentum: denominators[occurrence.original_occurrence]
                        .momentum
                        .clone(),
                    raw_to_parsed_sign: 1,
                }],
                "an empty-forest source must preserve numerator provenance and routing",
            );
        }
        Ok(())
    }

    #[test]
    fn unexpanded_massless_opposite_spelling_preserves_the_source_frame() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_massless_opposite_identity {
            edge [num=1 mass="UFO::ZERO"]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let production = GraphThreeDSource::new(&graph, &[])?;
        let production_parsed = production.to_three_d_parsed_graph()?;
        let mut denominators = (0..3)
            .map(|source_edge| {
                let source_edge = EdgeIndex(source_edge);
                FourDDenominator {
                    source_edge,
                    momentum: FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(source_edge))
                        .finish(),
                    mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
                    full_expr: Atom::one(),
                }
            })
            .collect::<Vec<_>>();
        denominators[1].momentum = -denominators[1].momentum.clone();
        denominators.reverse();
        let exact = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let exact_parsed = exact.to_three_d_parsed_graph()?;

        assert_eq!(exact_parsed, production_parsed);
        assert_eq!(
            exact.energy_edge_index_map(&exact_parsed),
            production.energy_edge_index_map(&production_parsed),
        );
        assert_eq!(
            exact
                .exact_source_energy_mapper()
                .unwrap()
                .source_edge_occurrences[&EdgeIndex(1)],
            vec![ExactSourceOwnerOccurrence {
                energy_edge_id: 1,
                raw_momentum: -FunctionBuilder::new(GS.emr_mom).add_arg(1).finish(),
                raw_to_parsed_sign: -1,
            }],
            "D(-Q) must compose its literal sign with the source routing sign without changing the odd physical numerator Q",
        );
        Ok(())
    }

    #[test]
    fn tagged_odd_hard_momentum_composes_the_d_neg_q_routing_sign() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_tagged_massless_opposite {
            edge [num=1 mass="UFO::ZERO"]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let mut denominators = (0..3)
            .map(|source_edge| {
                let source_edge = EdgeIndex(source_edge);
                FourDDenominator {
                    source_edge,
                    momentum: FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(source_edge))
                        .finish(),
                    mass_squared: Atom::Zero,
                    full_expr: Atom::one(),
                }
            })
            .collect::<Vec<_>>();
        let owner = EdgeIndex(1);
        denominators[usize::from(owner)].momentum =
            -denominators[usize::from(owner)].momentum.clone();
        let expected_energy = denominators[usize::from(owner)].on_shell_energy();
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let mapper = source.exact_source_energy_mapper().unwrap();
        let hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(owner))
            .finish();
        let tagged_hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_momentum_provenance_tag(
                Atom::num(usize::from(owner) as i64).as_view(),
                false,
                hard.as_view(),
            ))
            .add_arg(GS.cind(0))
            .finish();
        let edge_energy_map = (0..3)
            .map(|edge| LinearEnergyExpr::ose(EdgeIndex(edge), 1))
            .collect::<Vec<_>>();

        let mapped = mapper.map_numerator(
            &[LinearEnergyExpr::ose(owner, 1)],
            &edge_energy_map,
            &tagged_hard,
        )?;
        assert_eq!(
            mapped,
            crate::utils::ose_atom_from_index(owner),
            "rewriting D(Q) as D(-Q) must not reverse the fixed odd numerator Q",
        );
        assert_eq!(
            mapped
                .replace(crate::utils::ose_atom_from_index(owner))
                .with(expected_energy.clone()),
            expected_energy,
        );
        Ok(())
    }

    #[test]
    fn exact_full_source_preserves_production_carrier_orientations() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_production_basis {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
            c -> a [id=2 lmb_id=0]
            c -> d [id=3]
            d -> e [id=4]
            e -> c [id=5 lmb_id=1]
        })?;
        let production = GraphThreeDSource::new(&graph, &[])?;
        let configured = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        assert_ne!(production.outer_loop_edges, configured);

        let denominators = (0..6)
            .map(|edge| {
                let edge = EdgeIndex(edge);
                FourDDenominator {
                    source_edge: edge,
                    momentum: FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(edge))
                        .finish(),
                    mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                    full_expr: Atom::one(),
                }
            })
            .collect::<Vec<_>>();
        let exact = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        assert_eq!(exact.outer_loop_edges, production.outer_loop_edges);

        let production_parsed = production.to_three_d_parsed_graph()?;
        let exact_parsed = exact.to_three_d_parsed_graph()?;
        assert_eq!(exact_parsed.internal_edges.len(), denominators.len());
        for (exact_edge, original) in exact_parsed
            .internal_edges
            .iter()
            .zip(&exact.exact_local_to_original_occurrence)
        {
            let denominator = &denominators[*original];
            assert_eq!(
                exact_edge.signature,
                production_parsed.internal_edges[usize::from(denominator.source_edge)]
                    .signature
                    .clone(),
                "an unchanged exact denominator must preserve the production carrier orientation for odd numerators",
            );
        }

        let mut reversed_denominators = denominators.clone();
        reversed_denominators.reverse();
        let reversed = GraphThreeDSource::from_exact_denominators(&graph, &reversed_denominators)?;
        let reversed_parsed = reversed.to_three_d_parsed_graph()?;
        assert_eq!(reversed.outer_loop_edges, exact.outer_loop_edges);
        assert_eq!(reversed_parsed.loop_names, exact_parsed.loop_names);
        for (exact_edge, original) in reversed_parsed
            .internal_edges
            .iter()
            .zip(&reversed.exact_local_to_original_occurrence)
        {
            let denominator = &reversed_denominators[*original];
            assert_eq!(
                exact_edge.signature,
                production_parsed.internal_edges[usize::from(denominator.source_edge)]
                    .signature
                    .clone(),
                "denominator occurrence ordering must not change production carrier orientation",
            );
        }
        Ok(())
    }

    #[test]
    fn exact_shifted_factor_reversal_preserves_projected_affine_maps() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_shifted_reversal {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=2]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> outgoing [id=3]
        })?;
        let denominators = [EdgeIndex(0), EdgeIndex(1)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let options = graph.denominator_only_cff_3d_expression_options();
        let forward_source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let forward = graph
            .generate_3d_expression_for_4d_term(&forward_source, &options, &Atom::one(), None)?
            .0
            .expression
            .remap_energy_edge_indices(
                &forward_source
                    .physical_energy_edge_index_map()
                    .expect("exact source has a physical energy map"),
            );

        let reversed_denominators = [denominators[1].clone(), denominators[0].clone()];
        let reversed_source =
            GraphThreeDSource::from_exact_denominators(&graph, &reversed_denominators)?;
        let reversed = graph
            .generate_3d_expression_for_4d_term(&reversed_source, &options, &Atom::one(), None)?
            .0
            .expression
            .remap_energy_edge_indices(
                &reversed_source
                    .physical_energy_edge_index_map()
                    .expect("exact source has a physical energy map"),
            );

        assert_eq!(
            forward_source.to_three_d_parsed_graph()?,
            reversed_source.to_three_d_parsed_graph()?
        );
        assert_eq!(forward.orientations.len(), reversed.orientations.len());
        let mut unmatched = reversed.orientations.iter().collect::<Vec<_>>();
        for forward in &forward.orientations {
            let position = unmatched
                .iter()
                .position(|reversed| {
                    reversed.data.orientation == forward.data.orientation
                        && reversed.loop_energy_map == forward.loop_energy_map
                        && reversed.edge_energy_map == forward.edge_energy_map
                })
                .ok_or_else(|| {
                    eyre::eyre!(
                        "projected forward affine map {:?} is absent after factor reversal",
                        forward.edge_energy_map
                    )
                })?;
            unmatched.swap_remove(position);
        }
        assert!(unmatched.is_empty());
        Ok(())
    }

    #[test]
    fn unexpanded_external_source_reuses_production_frame_for_reversed_input() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_external_production_identity {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=2]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> outgoing [id=3]
        })?;
        let mut denominators = [EdgeIndex(1), EdgeIndex(0)].map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        denominators[0].momentum = -denominators[0].momentum.clone();

        let production_source = GraphThreeDSource::new(&graph, &[])?;
        let production_parsed = production_source.to_three_d_parsed_graph()?;
        let production_energy_map = production_source
            .energy_edge_index_map(&production_parsed)
            .expect("production source has a full energy map");
        let exact_source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let exact_parsed = exact_source.to_three_d_parsed_graph()?;

        assert_eq!(exact_parsed, production_parsed);
        assert_eq!(
            exact_source.energy_edge_index_map(&exact_parsed),
            Some(production_energy_map)
        );

        let factorized_numerator = (GS.emr_mom(EdgeIndex(0), GS.cind(0))
            + Atom::var(symbolica::symbol!("exact_external_identity::a")))
            * (GS.emr_mom(EdgeIndex(1), GS.cind(0))
                + Atom::var(symbolica::symbol!("exact_external_identity::b")));
        let (generated, mapper, plan, _) = graph.generate_3d_expression_for_4d_term(
            &exact_source,
            &graph.denominator_only_cff_3d_expression_options(),
            &factorized_numerator,
            None,
        )?;
        assert_eq!(
            mapper.source_edge_occurrences,
            BTreeMap::from([
                (
                    EdgeIndex(0),
                    vec![ExactSourceOwnerOccurrence {
                        energy_edge_id: 0,
                        raw_momentum: FunctionBuilder::new(GS.emr_mom).add_arg(0).finish(),
                        raw_to_parsed_sign: 1,
                    }],
                ),
                (
                    EdgeIndex(1),
                    vec![ExactSourceOwnerOccurrence {
                        energy_edge_id: 1,
                        raw_momentum: -FunctionBuilder::new(GS.emr_mom).add_arg(1).finish(),
                        raw_to_parsed_sign: -1,
                    }],
                ),
            ])
        );
        let orientation = generated
            .expression
            .orientations
            .first()
            .expect("production bubble has a causal orientation");
        let production_mapped = factorized_numerator
            .replace_multiple(crate::cff::expression::energy_map_replacements_gs(
                &orientation.edge_energy_map,
                &graph,
            ))
            .replace_multiple(
                exact_source
                    .exact_ose_replacements()
                    .expect("exact source has literal on-shell energies"),
            );
        let exact_mapped = mapper.map_planned_numerator(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            &plan,
        )?;
        assert_eq!(exact_mapped, production_mapped);
        assert!(matches!(exact_mapped.as_view(), AtomView::Mul(_)));
        Ok(())
    }

    #[test]
    fn exact_rank_deficient_source_keeps_the_complete_active_direction() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_rank_deficient {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
            c -> a [id=2 lmb_id=0]
            c -> d [id=3]
            d -> e [id=4]
            e -> c [id=5 lmb_id=1]
        })?;
        let parent_carriers = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let denominator = FourDDenominator {
            source_edge: EdgeIndex(0),
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(parent_carriers[0]))
                .finish()
                + FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(parent_carriers[1]))
                    .finish(),
            mass_squared: graph.underlying[EdgeIndex(0)].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };
        let raw_signature = denominator
            .momentum_signature(&graph, false)?
            .loop_signature;
        let denominators = [denominator];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        assert_eq!(source.inner_loop_count, 1);
        assert_eq!(source.active_loop_count(), 1);

        let source_coordinates = (0..2)
            .map(|source_column| {
                raw_signature.iter().enumerate().fold(
                    Rational::from(0),
                    |sum, (parent_column, coefficient)| {
                        sum + Rational::from(*coefficient)
                            * source.parent_loop_coordinates[parent_column][source_column].clone()
                    },
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(source_coordinates[0], Rational::from(0));
        assert_eq!(
            source_coordinates[1]
                .numerator_ref()
                .to_i64()
                .unwrap()
                .abs(),
            1
        );
        assert_eq!(source_coordinates[1].denominator_ref().to_i64(), Some(1));
        assert_eq!(
            source.exact_signatures[0].loop_signature,
            vec![source_coordinates[1].numerator_ref().to_i64().unwrap() as i32]
        );
        Ok(())
    }

    #[test]
    fn exact_rank_deficient_initial_cut_energy_uses_its_external_alias() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph exact_rank_deficient_initial_cut {
            edge [num=1 mass=1]
            node [num=1]

            c -> d [id=0 is_cut=0]
            a -> b [id=1]
            b -> c [id=2]
            c -> a [id=3 lmb_id=0]
            d -> e [id=4]
            e -> f [id=5]
            f -> d [id=6 lmb_id=1]
            f -> a [id=7]
        })?;
        let cut_edge = graph
            .get_edges_in_initial_state_cut()
            .into_iter()
            .exactly_one()
            .expect("test graph has one initial-state cut edge");
        let parent_carriers = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        assert_eq!(parent_carriers.len(), 2);

        // A cut carrier is a literal external alias even if a noncanonical
        // stored signature gives it a component outside this term's active span.
        graph.loop_momentum_basis.edge_signatures[cut_edge].internal[LoopIndex::from(1)] =
            SignOrZero::Plus;
        let denominator = FourDDenominator {
            source_edge: parent_carriers[0],
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(parent_carriers[0]))
                .finish(),
            mass_squared: graph.underlying[parent_carriers[0]]
                .particle
                .mass_atom()
                .pow(2),
            full_expr: Atom::one(),
        };
        let denominators = [denominator];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        assert_eq!(source.inner_loop_count, 1);
        assert_eq!(source.active_loop_count(), 1);

        let temporal_cut = GS.emr_mom(cut_edge, GS.cind(0));
        let signature = &graph.loop_momentum_basis.edge_signatures[cut_edge];
        let alias_signature = MomentumSignature {
            loop_signature: Vec::new(),
            external_signature: (&signature.external).into_iter().map(sign_to_i32).collect(),
        };
        let (external_id, external_sign) =
            initial_state_cut_external_alias(usize::from(cut_edge), &alias_signature)?;
        let external_edge = graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .copied()
            .nth(external_id)
            .expect("cut alias has a physical external carrier");
        let expected = match external_sign {
            -1 => -crate::utils::external_energy_atom_from_index(external_edge),
            1 => crate::utils::external_energy_atom_from_index(external_edge),
            _ => unreachable!("initial-state cut aliases have unit sign"),
        };
        let occurrence = graph.underlying.n_edges();
        let edge_energy_map = vec![LinearEnergyExpr::zero(); occurrence + 1];
        assert_eq!(
            source
                .exact_source_energy_mapper()
                .expect("exact source has an energy mapper")
                .map_numerator(&[LinearEnergyExpr::zero()], &edge_energy_map, &temporal_cut,)?,
            expected
        );
        Ok(())
    }

    #[test]
    fn exact_energy_mapper_preserves_factorized_affine_and_tensor_structure() -> Result<()> {
        test_initialise()?;
        let parent_edges = [EdgeIndex(3), EdgeIndex(5)];
        let mapped_edge = EdgeIndex(7);
        let exact_ose_edge = EdgeIndex(11);
        let external_edge = EdgeIndex(13);
        let exact_energy = Atom::var(symbolica::symbol!("three_d_source_test::exact_energy"));
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 1,
            parent_loop_coordinates: vec![
                vec![Rational::from(1), Rational::from(2)],
                vec![Rational::from(-1), Rational::from(3)],
            ],
            parent_loop_edges: parent_edges.to_vec(),
            edge_coordinates: vec![(
                mapped_edge,
                vec![Rational::from(0), Rational::from(-1)],
                vec![(external_edge, SignOrZero::Plus)],
            )],
            source_edge_occurrences: BTreeMap::new(),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: vec![Replacement::new(
                crate::utils::ose_atom_from_index(exact_ose_edge).to_pattern(),
                exact_energy.to_pattern(),
            )],
        };
        let active_map = LinearEnergyExpr {
            internal_terms: vec![(exact_ose_edge, Atom::num(-2))],
            external_terms: vec![(external_edge, Atom::num(3))],
            uniform_scale_coeff: Atom::num(-5),
            constant: Atom::num(7),
        };
        let active_energy = -Atom::num(2) * &exact_energy
            + Atom::num(3) * crate::utils::external_energy_atom_from_index(external_edge)
            - Atom::num(5) * Atom::var(GS.numerator_sampling_scale)
            + Atom::num(7);
        let loop_component = |loop_index: usize, index: &Atom| {
            FunctionBuilder::new(GS.loop_mom)
                .add_arg(loop_index)
                .add_arg(index)
                .finish()
        };

        // The inactive coordinates cancel before their placeholders are set to
        // zero, while every affine term in the active map remains explicit.
        let factor_a = Atom::var(symbolica::symbol!("three_d_source_test::factor_a"));
        let factor_b = Atom::var(symbolica::symbol!("three_d_source_test::factor_b"));
        let temporal_sum = loop_component(0, &GS.cind(0)) + loop_component(1, &GS.cind(0));
        let numerator =
            (temporal_sum + &factor_a) * (GS.emr_mom(mapped_edge, GS.cind(0)) + &factor_b);
        let mapped = mapper.map_numerator(std::slice::from_ref(&active_map), &[], &numerator)?;
        let expected = (Atom::num(5) * &active_energy + &factor_a)
            * (-&active_energy
                + crate::utils::external_energy_atom_from_index(external_edge)
                + &factor_b);
        assert_eq!(mapped, expected);
        let AtomView::Mul(product) = mapped.as_view() else {
            panic!("exact energy mapping should preserve the numerator product");
        };
        assert_eq!(
            product
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)))
                .count(),
            2
        );

        let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
        let abstract_loops = loop_component(0, &mink_index) + loop_component(1, &mink_index);
        assert_eq!(
            mapper.map_numerator(std::slice::from_ref(&active_map), &[], &abstract_loops)?,
            GS.emr_vec_index(parent_edges[0], &mink_index)
                + GS.emr_vec_index(parent_edges[1], &mink_index)
                + Atom::num(5) * &active_energy * GS.energy_delta(&mink_index)
        );
        assert_eq!(
            mapper.map_numerator(
                std::slice::from_ref(&active_map),
                &[],
                &GS.emr_mom(mapped_edge, &mink_index),
            )?,
            GS.emr_vec_index(mapped_edge, &mink_index)
                + (-&active_energy + crate::utils::external_energy_atom_from_index(external_edge))
                    * GS.energy_delta(&mink_index)
        );

        let concrete_spatial = loop_component(1, &GS.cind(2)) * GS.emr_mom(mapped_edge, GS.cind(3));
        assert_eq!(
            mapper.map_numerator(&[active_map], &[], &concrete_spatial)?,
            GS.emr_mom(parent_edges[1], GS.cind(2)) * GS.emr_mom(mapped_edge, GS.cind(3))
        );
        Ok(())
    }

    #[test]
    fn exact_energy_bounds_require_a_literal_occurrence() {
        let unique_edge = EdgeIndex(3);
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: Vec::new(),
            parent_loop_edges: Vec::new(),
            edge_coordinates: Vec::new(),
            source_edge_occurrences: BTreeMap::from([(
                unique_edge,
                vec![ExactSourceOwnerOccurrence {
                    energy_edge_id: 12,
                    raw_momentum: FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(unique_edge))
                        .finish(),
                    raw_to_parsed_sign: 1,
                }],
            )]),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: Vec::new(),
        };

        let candidates = mapper.equivalent_energy_candidates([unique_edge]).unwrap();
        let numerator = GS.emr_mom(unique_edge, GS.cind(0)).pow(2);
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([unique_edge])
            .plan_atom_assignment(&numerator, &candidates)
            .unwrap();
        assert_eq!(plan.energy_degree_bounds(), &[(12, 2)]);

        let missing_edge = EdgeIndex(5);
        let candidates = mapper.equivalent_energy_candidates([missing_edge]).unwrap();
        let missing = EnergyPowerAnalyzer::for_physical_emr_edges([missing_edge])
            .plan_atom_assignment(&GS.emr_mom(missing_edge, GS.cind(0)).pow(3), &candidates)
            .expect_err("an absent literal EMR occurrence must not use an affine fallback");
        assert!(
            missing
                .to_string()
                .contains("no certified equivalent exact-energy candidates")
        );
    }

    #[test]
    fn exact_energy_bounds_accept_equivalent_literal_occurrences() -> Result<()> {
        test_initialise()?;
        let repeated_edge = EdgeIndex(4);
        let occurrences = [EdgeIndex(13), EdgeIndex(14)];
        let energy = Atom::var(symbolica::symbol!(
            "three_d_source_test::equivalent_bound_energy"
        ));
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: Vec::new(),
            parent_loop_edges: Vec::new(),
            edge_coordinates: Vec::new(),
            source_edge_occurrences: BTreeMap::from([(
                repeated_edge,
                occurrences
                    .iter()
                    .copied()
                    .map(|occurrence| ExactSourceOwnerOccurrence {
                        energy_edge_id: usize::from(occurrence),
                        raw_momentum: FunctionBuilder::new(GS.emr_mom)
                            .add_arg(usize::from(repeated_edge))
                            .finish(),
                        raw_to_parsed_sign: 1,
                    })
                    .collect(),
            )]),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: occurrences
                .iter()
                .map(|occurrence| {
                    Replacement::new(
                        crate::utils::ose_atom_from_index(*occurrence).to_pattern(),
                        energy.to_pattern(),
                    )
                })
                .collect(),
        };

        let candidates = mapper.equivalent_energy_candidates([repeated_edge])?;
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([repeated_edge])
            .plan_atom_assignment(&GS.emr_mom(repeated_edge, GS.cind(0)).pow(4), &candidates)?;
        assert_eq!(
            plan.energy_degree_bounds(),
            &[(usize::from(occurrences[0]), 4)],
            "an original untagged factor remains fixed on the canonical source occurrence"
        );
        Ok(())
    }

    #[test]
    fn exact_energy_bounds_choose_the_lowest_certified_occurrence_deterministically() -> Result<()>
    {
        test_initialise()?;
        let repeated_edge = EdgeIndex(4);
        let occurrences = [EdgeIndex(13), EdgeIndex(14)];
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: Vec::new(),
            parent_loop_edges: Vec::new(),
            edge_coordinates: Vec::new(),
            source_edge_occurrences: BTreeMap::from([(
                repeated_edge,
                occurrences
                    .iter()
                    .rev()
                    .copied()
                    .map(|occurrence| ExactSourceOwnerOccurrence {
                        energy_edge_id: usize::from(occurrence),
                        raw_momentum: FunctionBuilder::new(GS.emr_mom)
                            .add_arg(usize::from(repeated_edge))
                            .finish(),
                        raw_to_parsed_sign: 1,
                    })
                    .collect(),
            )]),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: occurrences
                .iter()
                .map(|occurrence| {
                    Replacement::new(
                        crate::utils::ose_atom_from_index(*occurrence).to_pattern(),
                        Atom::one().to_pattern(),
                    )
                })
                .collect(),
        };

        let candidates = mapper.equivalent_energy_candidates([repeated_edge])?;
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([repeated_edge])
            .plan_atom_assignment(&GS.emr_mom(repeated_edge, GS.cind(0)).pow(4), &candidates)?;
        assert_eq!(
            plan.energy_degree_bounds(),
            &[(usize::from(occurrences[0]), 4)]
        );
        Ok(())
    }

    #[test]
    fn exact_energy_bounds_keep_original_on_canonical_occurrence_across_routing_signs() -> Result<()>
    {
        test_initialise()?;
        let repeated_edge = EdgeIndex(4);
        let occurrences = [EdgeIndex(13), EdgeIndex(14), EdgeIndex(15)];
        let energy = Atom::var(symbolica::symbol!(
            "three_d_source_test::signed_class_energy"
        ));

        for signs in [[1, -1, -1], [1, 1, -1]] {
            let mapper = ExactSourceEnergyMapper {
                inactive_loop_count: 0,
                parent_loop_coordinates: vec![vec![Rational::from(1)]],
                parent_loop_edges: vec![repeated_edge],
                edge_coordinates: vec![(repeated_edge, vec![Rational::from(1)], Vec::new())],
                source_edge_occurrences: BTreeMap::from([(
                    repeated_edge,
                    occurrences
                        .iter()
                        .zip(signs)
                        .map(|(occurrence, sign)| ExactSourceOwnerOccurrence {
                            energy_edge_id: usize::from(*occurrence),
                            raw_momentum: FunctionBuilder::new(GS.emr_mom)
                                .add_arg(usize::from(repeated_edge))
                                .finish(),
                            raw_to_parsed_sign: sign,
                        })
                        .collect(),
                )]),
                cut_alias_edges: AHashSet::new(),
                exact_ose_replacements: occurrences
                    .iter()
                    .map(|occurrence| {
                        Replacement::new(
                            crate::utils::ose_atom_from_index(*occurrence).to_pattern(),
                            energy.to_pattern(),
                        )
                    })
                    .collect(),
            };
            let candidates = mapper.equivalent_energy_candidates([repeated_edge])?;
            let numerator = GS.emr_mom(repeated_edge, GS.cind(0)).pow(2);
            let plan = EnergyPowerAnalyzer::for_physical_emr_edges([repeated_edge])
                .plan_atom_assignment(&numerator, &candidates)?;
            assert_eq!(
                plan.energy_degree_bounds(),
                &[(usize::from(occurrences[0]), 2)],
                "routing signs do not make an original factor eligible for cross-copy dispatch"
            );

            let mut edge_energy_map = vec![LinearEnergyExpr::zero(); 16];
            for (occurrence, sign) in occurrences.iter().zip(signs) {
                edge_energy_map[usize::from(*occurrence)] =
                    LinearEnergyExpr::ose(*occurrence, sign);
            }
            assert_eq!(
                mapper.map_planned_numerator(
                    &[LinearEnergyExpr::ose(occurrences[0], 1)],
                    &edge_energy_map,
                    &plan,
                )?,
                energy.pow(2)
            );
        }
        Ok(())
    }

    #[test]
    fn denominator_derived_dispatch_composes_each_occurrence_routing_sign() -> Result<()> {
        test_initialise()?;
        let owner = EdgeIndex(4);
        let occurrences = [EdgeIndex(13), EdgeIndex(14)];
        let energy = Atom::var(symbolica::symbol!(
            "three_d_source_test::derived_dispatch_energy"
        ));
        let hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(owner))
            .finish();
        let left_shift = Atom::var(symbolica::symbol!(
            "three_d_source_test::derived_dispatch_left"
        ));
        let right_shift = Atom::var(symbolica::symbol!(
            "three_d_source_test::derived_dispatch_right"
        ));
        let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
        let reference = GS.emr_mom(EdgeIndex(7), &mink_index);

        for hard_sign in [1, -1] {
            let immutable_hard = if hard_sign == 1 {
                hard.clone()
            } else {
                -hard.clone()
            };
            let tag = GS.uv_momentum_provenance_tag(
                Atom::num(usize::from(owner) as i64).as_view(),
                true,
                immutable_hard.as_view(),
            );
            let tagged_temporal = FunctionBuilder::new(GS.emr_mom)
                .add_arg(tag.as_view())
                .add_arg(GS.cind(0))
                .finish();
            let tagged_abstract = FunctionBuilder::new(GS.emr_mom)
                .add_arg(tag.as_view())
                .add_arg(mink_index.as_view())
                .finish();
            let temporal_numerator =
                (&tagged_temporal + &left_shift) * (&tagged_temporal + &right_shift);
            let abstract_numerator =
                function!(GS.dot, reference.as_view(), tagged_abstract.as_view());

            for (raw_signs, parsed_sign) in
                [([1, -1], 1), ([-1, 1], 1), ([1, -1], -1), ([-1, 1], -1)]
            {
                let mapper = ExactSourceEnergyMapper {
                    inactive_loop_count: 0,
                    parent_loop_coordinates: Vec::new(),
                    parent_loop_edges: Vec::new(),
                    edge_coordinates: vec![(owner, Vec::new(), Vec::new())],
                    source_edge_occurrences: BTreeMap::from([(
                        owner,
                        occurrences
                            .iter()
                            .zip(raw_signs)
                            .map(|(occurrence, raw_sign)| ExactSourceOwnerOccurrence {
                                energy_edge_id: usize::from(*occurrence),
                                raw_momentum: if raw_sign == 1 {
                                    hard.clone()
                                } else {
                                    -hard.clone()
                                },
                                // P = parsed_sign * H and R = raw_sign * H,
                                // hence P = (parsed_sign * raw_sign) * R.
                                raw_to_parsed_sign: parsed_sign * raw_sign,
                            })
                            .collect(),
                    )]),
                    cut_alias_edges: AHashSet::new(),
                    exact_ose_replacements: occurrences
                        .iter()
                        .map(|occurrence| {
                            Replacement::new(
                                crate::utils::ose_atom_from_index(*occurrence).to_pattern(),
                                energy.to_pattern(),
                            )
                        })
                        .collect(),
                };
                let candidates = mapper.equivalent_energy_candidates([owner])?;
                let mut edge_energy_map = vec![LinearEnergyExpr::zero(); 15];
                for occurrence in occurrences {
                    edge_energy_map[usize::from(occurrence)] =
                        LinearEnergyExpr::ose(occurrence, parsed_sign);
                }
                let expected_bounds = [
                    (usize::from(occurrences[0]), 1),
                    (usize::from(occurrences[1]), 1),
                ];
                let temporal_plan = EnergyPowerAnalyzer::for_physical_emr_edges([owner])
                    .plan_atom_assignment(&temporal_numerator, &candidates)?;
                assert_eq!(temporal_plan.energy_degree_bounds(), &expected_bounds);
                let signed_energy = Atom::num(hard_sign) * &energy;
                assert_eq!(
                    mapper.map_planned_numerator(&[], &edge_energy_map, &temporal_plan)?,
                    (&signed_energy + &left_shift) * (&signed_energy + &right_shift),
                    "each dispatched temporal factor must retain hard sign {hard_sign} for raw signs {raw_signs:?} and parsed sign {parsed_sign}",
                );

                let expected_vector = Atom::num(hard_sign) * GS.emr_vec_index(owner, &mink_index)
                    + Atom::num(hard_sign) * &energy * GS.energy_delta(&mink_index);
                for occurrence in occurrences {
                    let forced_candidates =
                        EquivalentEnergyCandidates::try_from_source_occurrences([(
                            owner,
                            vec![usize::from(occurrence)],
                        )])?;
                    let abstract_plan = EnergyPowerAnalyzer::for_physical_emr_edges([owner])
                        .plan_atom_assignment(&abstract_numerator, &forced_candidates)?;
                    assert_eq!(
                        abstract_plan.energy_degree_bounds(),
                        &[(usize::from(occurrence), 1)]
                    );
                    assert_eq!(
                        mapper.map_planned_numerator(&[], &edge_energy_map, &abstract_plan)?,
                        function!(GS.dot, reference.as_view(), expected_vector.as_view()),
                        "the abstract factor on occurrence {} must retain the same spatial and temporal hard sign",
                        usize::from(occurrence),
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn exact_energy_mapper_uses_source_loop_map_across_distinct_occurrence_pole_maps() -> Result<()>
    {
        test_initialise()?;
        let repeated_edge = EdgeIndex(4);
        let occurrences = [EdgeIndex(13), EdgeIndex(14)];
        let energy = Atom::var(symbolica::symbol!(
            "three_d_source_test::routing_reversed_energy"
        ));
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: vec![vec![Rational::from(1)]],
            parent_loop_edges: vec![repeated_edge],
            edge_coordinates: vec![(repeated_edge, vec![Rational::from(1)], Vec::new())],
            source_edge_occurrences: BTreeMap::from([(
                repeated_edge,
                occurrences
                    .iter()
                    .copied()
                    .map(|occurrence| ExactSourceOwnerOccurrence {
                        energy_edge_id: usize::from(occurrence),
                        raw_momentum: FunctionBuilder::new(GS.emr_mom)
                            .add_arg(usize::from(repeated_edge))
                            .finish(),
                        raw_to_parsed_sign: 1,
                    })
                    .collect(),
            )]),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: occurrences
                .iter()
                .map(|occurrence| {
                    Replacement::new(
                        crate::utils::ose_atom_from_index(*occurrence).to_pattern(),
                        energy.to_pattern(),
                    )
                })
                .collect(),
        };
        let mut edge_energy_map = vec![LinearEnergyExpr::zero(); 15];
        edge_energy_map[usize::from(occurrences[0])] = LinearEnergyExpr::ose(occurrences[0], 1);
        edge_energy_map[usize::from(occurrences[1])] = LinearEnergyExpr::ose(occurrences[1], -1);
        let loop_energy_map = [LinearEnergyExpr::ose(occurrences[0], 1)];
        let hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(repeated_edge))
            .finish();
        let fixed = FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_momentum_provenance_tag(
                Atom::num(usize::from(repeated_edge) as i64).as_view(),
                false,
                hard.as_view(),
            ))
            .add_arg(GS.cind(0))
            .finish();

        assert_eq!(
            mapper.map_numerator(&loop_energy_map, &edge_energy_map, &fixed)?,
            energy
        );

        edge_energy_map[usize::from(occurrences[1])] = LinearEnergyExpr::ose(occurrences[1], -2);
        assert_eq!(
            mapper.map_numerator(&loop_energy_map, &edge_energy_map, &fixed)?,
            energy
        );
        Ok(())
    }

    #[test]
    fn exact_energy_mapper_separates_fixed_owner_routing_from_derived_dispatch() -> Result<()> {
        test_initialise()?;
        let owners = [EdgeIndex(4), EdgeIndex(5)];
        let occurrences = [
            [EdgeIndex(13), EdgeIndex(15)],
            [EdgeIndex(14), EdgeIndex(16)],
        ];
        let source_energy = Atom::var(symbolica::symbol!(
            "three_d_source_test::fixed_source_energy"
        ));
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: vec![vec![Rational::from(1)]],
            parent_loop_edges: vec![owners[0]],
            edge_coordinates: owners
                .iter()
                .copied()
                .map(|owner| (owner, vec![Rational::from(1)], Vec::new()))
                .collect(),
            source_edge_occurrences: owners
                .iter()
                .copied()
                .zip(occurrences)
                .enumerate()
                .map(|(owner_position, (owner, occurrences))| {
                    (
                        owner,
                        occurrences
                            .into_iter()
                            .map(|occurrence| ExactSourceOwnerOccurrence {
                                energy_edge_id: usize::from(occurrence),
                                raw_momentum: FunctionBuilder::new(GS.emr_mom)
                                    .add_arg(usize::from(owner))
                                    .finish(),
                                // Mirror the GL24 T^0 triangle: the second
                                // even denominator is routed oppositely by
                                // canonicalization. Its oriented occurrence
                                // sample already carries that incidence sign
                                // for a fixed factor; only a dispatched
                                // derivative factor crosses this frame again.
                                raw_to_parsed_sign: if owner_position == 0 { 1 } else { -1 },
                            })
                            .collect(),
                    )
                })
                .collect(),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: Vec::new(),
        };
        let contacts = [
            Atom::var(symbolica::symbol!("three_d_source_test::left_contact")),
            Atom::var(symbolica::symbol!("three_d_source_test::right_contact")),
        ];
        let mut edge_energy_map = vec![LinearEnergyExpr::zero(); 17];
        for (owner_occurrences, contact) in occurrences.iter().zip(&contacts) {
            edge_energy_map[usize::from(owner_occurrences[0])].constant = contact.clone();
            edge_energy_map[usize::from(owner_occurrences[1])].constant = Atom::var(
                symbolica::symbol!("three_d_source_test::unselected_contact"),
            );
        }
        let mut source_loop_energy = LinearEnergyExpr::zero();
        source_loop_energy.constant = source_energy.clone();
        let loop_energy_map = [source_loop_energy];
        let tagged = |owner: EdgeIndex, role: UvMomentumProvenanceRole, index: &Atom| {
            let hard = FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(owner))
                .finish();
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(usize::from(owner) as i64).as_view(),
                    role,
                    hard.as_view(),
                ))
                .add_arg(index.as_view())
                .finish()
        };
        let candidates = mapper.equivalent_energy_candidates(owners)?;
        let fixed_temporal = tagged(
            owners[0],
            UvMomentumProvenanceRole::TaylorFixed,
            &GS.cind(0),
        ) * tagged(
            owners[1],
            UvMomentumProvenanceRole::TaylorFixed,
            &GS.cind(0),
        );
        let fixed_temporal_plan = EnergyPowerAnalyzer::for_physical_emr_edges(owners)
            .plan_atom_assignment(&fixed_temporal, &candidates)?;

        assert_eq!(
            fixed_temporal_plan.energy_degree_bounds(),
            &[
                (usize::from(occurrences[0][0]), 1),
                (usize::from(occurrences[1][0]), 1),
            ],
            "fixed factors must stay on the canonical occurrence of their own physical owner",
        );
        // Exercise every occurrence-orientation sign pair and every branch in
        // which either selected contact sample vanishes. A fixed numerator
        // remains in its canonical occurrence's Laurent functional without a
        // second parsed-incidence sign.
        for left_sign in [-1, 0, 1] {
            for right_sign in [-1, 0, 1] {
                edge_energy_map[usize::from(occurrences[0][0])].constant =
                    Atom::num(left_sign) * &contacts[0];
                edge_energy_map[usize::from(occurrences[1][0])].constant =
                    Atom::num(right_sign) * &contacts[1];
                let actual = mapper.map_planned_numerator(
                    &loop_energy_map,
                    &edge_energy_map,
                    &fixed_temporal_plan,
                )?;
                let expected = Atom::num(left_sign * right_sign) * &contacts[0] * &contacts[1];
                assert!(
                    (actual - expected).expand().is_zero(),
                    "fixed factors must consume owner-local samples without a second parsed-routing sign ({left_sign:+}, {right_sign:+})",
                );
            }
        }

        for (owner_occurrences, contact) in occurrences.iter().zip(&contacts) {
            edge_energy_map[usize::from(owner_occurrences[0])].constant = contact.clone();
        }

        let abstract_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
        let fixed_abstract = tagged(
            owners[0],
            UvMomentumProvenanceRole::TaylorFixed,
            &abstract_index,
        ) * tagged(
            owners[1],
            UvMomentumProvenanceRole::TaylorFixed,
            &abstract_index,
        );
        let fixed_abstract_plan = EnergyPowerAnalyzer::for_physical_emr_edges(owners)
            .plan_atom_assignment(&fixed_abstract, &candidates)?;
        let mapped_vectors = owners.map(|owner| GS.emr_vec_index(owner, &abstract_index));
        assert_eq!(
            mapper.map_planned_numerator(
                &loop_energy_map,
                &edge_energy_map,
                &fixed_abstract_plan,
            )?,
            (&mapped_vectors[0] + &contacts[0] * GS.energy_delta(&abstract_index))
                * (&mapped_vectors[1] + &contacts[1] * GS.energy_delta(&abstract_index)),
            "abstract fixed factors must keep their spatial part and consume the oriented owner-local sample once",
        );

        let derived_temporal = tagged(
            owners[0],
            UvMomentumProvenanceRole::DenominatorDerived,
            &GS.cind(0),
        ) * tagged(
            owners[1],
            UvMomentumProvenanceRole::DenominatorDerived,
            &GS.cind(0),
        );
        let derived_temporal_plan = EnergyPowerAnalyzer::for_physical_emr_edges(owners)
            .plan_atom_assignment(&derived_temporal, &candidates)?;
        assert_eq!(
            mapper.map_planned_numerator(
                &loop_energy_map,
                &edge_energy_map,
                &derived_temporal_plan,
            )?,
            -&contacts[0] * &contacts[1],
            "denominator-derived factors must retain occurrence-local contact samples and routing signs",
        );

        let derived_abstract = tagged(
            owners[0],
            UvMomentumProvenanceRole::DenominatorDerived,
            &abstract_index,
        ) * tagged(
            owners[1],
            UvMomentumProvenanceRole::DenominatorDerived,
            &abstract_index,
        );
        let derived_abstract_plan = EnergyPowerAnalyzer::for_physical_emr_edges(owners)
            .plan_atom_assignment(&derived_abstract, &candidates)?;
        assert_eq!(
            mapper.map_planned_numerator(
                &loop_energy_map,
                &edge_energy_map,
                &derived_abstract_plan,
            )?,
            (&mapped_vectors[0] + &contacts[0] * GS.energy_delta(&abstract_index))
                * (&mapped_vectors[1] - &contacts[1] * GS.energy_delta(&abstract_index)),
            "abstract denominator-derived factors must keep their spatial part and occurrence-local temporal samples",
        );
        Ok(())
    }

    #[test]
    fn exact_energy_mapper_uses_loop_maps_for_source_coordinates() -> Result<()> {
        test_initialise()?;
        let literal_edge = EdgeIndex(3);
        let reconstructed_edge = EdgeIndex(5);
        let carrier_occurrence = EdgeIndex(11);
        let literal_occurrence = EdgeIndex(12);
        let duplicate_occurrence = EdgeIndex(13);
        let duplicate_ose = EdgeIndex(14);
        let fallback = EdgeIndex(17);
        let external_edge = EdgeIndex(19);
        let literal_energy = Atom::var(symbolica::symbol!("three_d_source_test::literal_energy"));
        let mut mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: vec![vec![Rational::from(1)]],
            parent_loop_edges: vec![reconstructed_edge],
            edge_coordinates: vec![
                (literal_edge, vec![Rational::from(1)], Vec::new()),
                (reconstructed_edge, vec![Rational::from(1)], Vec::new()),
            ],
            source_edge_occurrences: BTreeMap::from([(
                literal_edge,
                vec![ExactSourceOwnerOccurrence {
                    energy_edge_id: usize::from(literal_occurrence),
                    raw_momentum: FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(literal_edge))
                        .finish(),
                    raw_to_parsed_sign: 1,
                }],
            )]),
            cut_alias_edges: AHashSet::new(),
            exact_ose_replacements: vec![
                Replacement::new(
                    crate::utils::ose_atom_from_index(literal_occurrence).to_pattern(),
                    literal_energy.to_pattern(),
                ),
                Replacement::new(
                    crate::utils::ose_atom_from_index(duplicate_occurrence).to_pattern(),
                    literal_energy.to_pattern(),
                ),
                Replacement::new(
                    crate::utils::ose_atom_from_index(duplicate_ose).to_pattern(),
                    literal_energy.to_pattern(),
                ),
            ],
        };
        let loop_energy_map = [LinearEnergyExpr::ose(fallback, 1)];
        let mut edge_energy_map =
            vec![LinearEnergyExpr::zero(); usize::from(duplicate_occurrence) + 1];
        edge_energy_map[usize::from(carrier_occurrence)] = LinearEnergyExpr {
            internal_terms: vec![(carrier_occurrence, Atom::num(2))],
            external_terms: vec![(external_edge, Atom::num(3))],
            uniform_scale_coeff: Atom::num(-5),
            constant: Atom::num(7),
        };
        edge_energy_map[usize::from(literal_occurrence)] = LinearEnergyExpr {
            internal_terms: vec![(literal_occurrence, Atom::num(3))],
            external_terms: vec![(external_edge, Atom::num(-2))],
            uniform_scale_coeff: Atom::num(4),
            constant: Atom::num(-1),
        };
        let loop_temporal = FunctionBuilder::new(GS.loop_mom)
            .add_arg(0)
            .add_arg(GS.cind(0))
            .finish();
        let factor_a = Atom::var(symbolica::symbol!("three_d_source_test::edge_factor"));
        let factor_b = Atom::var(symbolica::symbol!(
            "three_d_source_test::reconstructed_factor"
        ));
        let factor_c = Atom::var(symbolica::symbol!("three_d_source_test::loop_factor"));
        let numerator = (GS.emr_mom(literal_edge, GS.cind(0)) + &factor_a)
            * (GS.emr_mom(reconstructed_edge, GS.cind(0)) + &factor_b)
            * (loop_temporal + &factor_c);
        let mapped = mapper.map_numerator(&loop_energy_map, &edge_energy_map, &numerator)?;
        let source_sample = crate::utils::ose_atom_from_index(fallback);
        let literal_sample = Atom::num(3) * &literal_energy
            - Atom::num(2) * crate::utils::external_energy_atom_from_index(external_edge)
            + Atom::num(4) * Atom::var(GS.numerator_sampling_scale)
            - Atom::num(1);
        assert_eq!(
            mapped,
            (literal_sample.clone() + &factor_a)
                * (source_sample.clone() + &factor_b)
                * (source_sample.clone() + &factor_c)
        );

        edge_energy_map[usize::from(carrier_occurrence)] = LinearEnergyExpr::zero();
        assert_eq!(
            mapper.map_numerator(
                &loop_energy_map,
                &edge_energy_map,
                &GS.emr_mom(reconstructed_edge, GS.cind(0)).pow(2),
            )?,
            source_sample.clone().pow(2)
        );

        edge_energy_map[usize::from(carrier_occurrence)] =
            LinearEnergyExpr::ose(carrier_occurrence, 2);
        mapper
            .source_edge_occurrences
            .get_mut(&literal_edge)
            .unwrap()
            .push(ExactSourceOwnerOccurrence {
                energy_edge_id: usize::from(duplicate_occurrence),
                raw_momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(literal_edge))
                    .finish(),
                raw_to_parsed_sign: 1,
            });
        edge_energy_map[usize::from(duplicate_occurrence)] = LinearEnergyExpr {
            internal_terms: vec![(duplicate_ose, Atom::num(3))],
            external_terms: vec![(external_edge, Atom::num(-2))],
            uniform_scale_coeff: Atom::num(4),
            constant: Atom::num(-1),
        };
        assert_eq!(
            mapper.map_numerator(&loop_energy_map, &edge_energy_map, &numerator)?,
            (literal_sample.clone() + &factor_a)
                * (source_sample.clone() + &factor_b)
                * (source_sample + &factor_c)
        );

        edge_energy_map[usize::from(duplicate_occurrence)].internal_terms[0].1 = Atom::num(-3);
        let factorized_literal_numerator = (GS.emr_mom(literal_edge, GS.cind(0)) + &factor_a)
            * (GS.emr_mom(literal_edge, GS.cind(0)) + &factor_b);
        let candidates = mapper.equivalent_energy_candidates([literal_edge])?;
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([literal_edge])
            .plan_atom_assignment(&factorized_literal_numerator, &candidates)?;
        assert_eq!(
            plan.energy_degree_bounds(),
            &[(usize::from(literal_occurrence), 2)],
            "both original factors remain fixed on the canonical literal occurrence"
        );
        assert_eq!(
            mapper.map_planned_numerator(&loop_energy_map, &edge_energy_map, &plan)?,
            (literal_sample.clone() + factor_a) * (literal_sample + factor_b)
        );
        Ok(())
    }

    #[test]
    fn exact_energy_mapper_preserves_cut_external_aliases() -> Result<()> {
        test_initialise()?;
        let cut_edge = EdgeIndex(3);
        let loop_edge = EdgeIndex(5);
        let external_edge = EdgeIndex(7);
        let mapper = ExactSourceEnergyMapper {
            inactive_loop_count: 0,
            parent_loop_coordinates: vec![vec![Rational::from(1)]],
            parent_loop_edges: vec![loop_edge],
            edge_coordinates: vec![(
                cut_edge,
                vec![Rational::from(0)],
                vec![(external_edge, SignOrZero::Minus)],
            )],
            source_edge_occurrences: BTreeMap::from([(
                cut_edge,
                [1, 2]
                    .into_iter()
                    .map(|energy_edge_id| ExactSourceOwnerOccurrence {
                        energy_edge_id,
                        raw_momentum: FunctionBuilder::new(GS.emr_mom)
                            .add_arg(usize::from(cut_edge))
                            .finish(),
                        raw_to_parsed_sign: 1,
                    })
                    .collect(),
            )]),
            cut_alias_edges: AHashSet::from([cut_edge]),
            exact_ose_replacements: Vec::new(),
        };
        let edge_energy_map = [
            LinearEnergyExpr::zero(),
            LinearEnergyExpr::ose(EdgeIndex(1), 2),
            LinearEnergyExpr::ose(EdgeIndex(2), -3),
        ];

        assert_eq!(
            mapper.map_numerator(
                &[LinearEnergyExpr::zero()],
                &edge_energy_map,
                &GS.emr_mom(cut_edge, GS.cind(0)),
            )?,
            -crate::utils::external_energy_atom_from_index(external_edge)
        );
        Ok(())
    }

    #[test]
    fn exact_energy_mapper_uses_selected_fixed_owner_sample_with_uv_mass() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_literal_momentum {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let carrier = graph.loop_momentum_basis.loop_edges[LoopIndex::from(0)];
        let ownership_edge = EdgeIndex(1);
        assert_ne!(carrier, ownership_edge);
        let denominator = FourDDenominator {
            source_edge: ownership_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(carrier))
                .finish(),
            mass_squared: Atom::var(GS.m_uv_expansion).pow(2),
            full_expr: Atom::one(),
        };
        let expected_energy = denominator.on_shell_energy();
        let denominators = [denominator];
        let source = GraphThreeDSource::from_exact_denominators(&graph, &denominators)?;
        let mapper = source
            .exact_source_energy_mapper()
            .expect("exact UV source has an energy mapper");
        let occurrence = EdgeIndex(graph.underlying.n_edges());
        assert_eq!(
            mapper.source_edge_occurrences[&ownership_edge],
            vec![ExactSourceOwnerOccurrence {
                energy_edge_id: usize::from(occurrence),
                raw_momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(carrier))
                    .finish(),
                raw_to_parsed_sign: 1,
            }]
        );
        let mut edge_energy_map = vec![LinearEnergyExpr::zero(); usize::from(occurrence) + 1];
        edge_energy_map[usize::from(occurrence)] = LinearEnergyExpr::ose(occurrence, 2);

        let hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(carrier))
            .finish();
        let tagged_hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_momentum_provenance_tag(
                Atom::num(usize::from(ownership_edge) as i64).as_view(),
                false,
                hard.as_view(),
            ))
            .add_arg(GS.cind(0))
            .finish();
        assert_eq!(
            mapper.map_numerator(
                &[LinearEnergyExpr::ose(occurrence, 1)],
                &edge_energy_map,
                &tagged_hard.pow(2),
            )?,
            Atom::num(4) * expected_energy.pow(2),
            "the fixed factor must consume the selected occurrence sample 2*OSE before the UV-mass OSE replacement",
        );
        Ok(())
    }

    #[test]
    fn exact_uv_component_inherits_source_minor_and_rejects_wrong_provenance() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_source_minor {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
            c -> a [id=2 lmb_id=0]
            c -> d [id=3]
            d -> e [id=4]
            e -> c [id=5 lmb_id=1]
        })?;
        let uv_edges = (0..graph.underlying.n_edges())
            .map(EdgeIndex)
            .collect::<Vec<_>>();
        let uv_mass = Atom::var(GS.m_uv_expansion).pow(2);
        let denominators = uv_edges
            .iter()
            .copied()
            .map(|source_edge| FourDDenominator {
                source_edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(source_edge))
                    .finish(),
                mass_squared: uv_mass.clone(),
                full_expr: Atom::one(),
            })
            .collect::<Vec<_>>();
        let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &denominators,
            uv_edges.iter().copied(),
        )?;
        let parsed = source.to_three_d_parsed_graph()?;

        assert!(
            three_dimensional_reps::validate_parsed_graph(&parsed).ok,
            "the complete UV component must retain its source-graph momentum balance"
        );
        assert!(
            parsed
                .node_name_to_internal
                .keys()
                .all(|name| !name.starts_with("__gammaloop_exact_power_")),
            "unraised source edges must not acquire power-chain vertices"
        );
        for (occurrence, original) in source
            .exact_local_to_original_occurrence
            .iter()
            .copied()
            .enumerate()
        {
            let source_edge = denominators[original].source_edge;
            let (
                _,
                HedgePair::Paired {
                    source: tail,
                    sink: head,
                },
            ) = graph[&source_edge]
            else {
                unreachable!("the fixture contains only paired source edges")
            };
            let mut parsed_endpoints = [
                parsed.internal_edges[occurrence].tail,
                parsed.internal_edges[occurrence].head,
            ];
            parsed_endpoints.sort_unstable();
            let mut source_endpoints = [
                parsed.node_name_to_internal[&format!("uv{}", usize::from(graph.node_id(tail)))],
                parsed.node_name_to_internal[&format!("uv{}", usize::from(graph.node_id(head)))],
            ];
            source_endpoints.sort_unstable();
            assert_eq!(
                parsed_endpoints, source_endpoints,
                "each exact UV denominator must inherit its own source-edge incidence; canonical D(Q)=D(-Q) routing may reverse both momentum and direction"
            );
        }

        let omitted_edge = EdgeIndex(1);
        let contracted_denominators = denominators
            .iter()
            .filter(|denominator| denominator.source_edge != omitted_edge)
            .cloned()
            .collect::<Vec<_>>();
        let contracted = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &contracted_denominators,
            uv_edges.iter().copied(),
        )?
        .to_three_d_parsed_graph()?;
        let (
            _,
            HedgePair::Paired {
                source: tail,
                sink: head,
            },
        ) = graph[&omitted_edge]
        else {
            unreachable!("the omitted UV source edge is paired")
        };
        assert_eq!(
            contracted.node_name_to_internal[&format!("uv{}", usize::from(graph.node_id(tail)))],
            contracted.node_name_to_internal[&format!("uv{}", usize::from(graph.node_id(head)))],
            "an absent UV denominator is contracted in the known source minor"
        );
        assert!(three_dimensional_reps::validate_parsed_graph(&contracted).ok);

        let mut malformed = denominators.clone();
        malformed[0].momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let error = match GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph, &malformed, uv_edges,
        ) {
            Ok(_) => panic!("wrong provenance must not be repaired by changing source incidence"),
            Err(error) => error,
        };
        assert!(
            error
                .to_string()
                .contains("has no unique routing for rewritten exact signature"),
            "unexpected provenance-validation error: {error}"
        );
        Ok(())
    }

    #[test]
    fn exact_source_routes_a_hard_uv_row_modulo_the_soft_cograph_span() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_hard_soft_routing {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v5 [id=0]
            v0 -> v2 [id=1 lmb_id=0]
            v3 -> v0 [id=2]
            v0 -> v5 [id=3]
            v2 -> v1 [id=4]
            v1 -> v3 [id=5 lmb_id=1]
            v1 -> v4 [id=6]
            v2 -> v3 [id=7]
            v4 -> outgoing [id=8]
        })?;
        let cograph_edges = [EdgeIndex(1), EdgeIndex(2)];
        let uv_edges = [EdgeIndex(4), EdgeIndex(5), EdgeIndex(7)];
        let hard_carrier = EdgeIndex(5);
        let physical_denominators = cograph_edges.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let uv_mass = Atom::var(GS.m_uv_expansion).pow(2);
        let uv_denominators = uv_edges.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(hard_carrier))
                .finish(),
            mass_squared: uv_mass.clone(),
            full_expr: Atom::one(),
        });
        let denominators = physical_denominators
            .into_iter()
            .chain(uv_denominators)
            .collect::<Vec<_>>();
        let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &denominators,
            uv_edges,
        )?;

        let hard_soft_owner = uv_edges
            .into_iter()
            .find(|owner| {
                source.outer_loop_signature(*owner).is_ok_and(|row| {
                    row.iter().filter(|coefficient| **coefficient != 0).count() > 1
                })
            })
            .expect("one source UV edge carries both the hard and soft loop directions");
        let original = denominators
            .iter()
            .position(|denominator| denominator.source_edge == hard_soft_owner)
            .expect("the mixed source owner retains one rewritten denominator");
        let exact_signature = source.exact_signatures[original].canonical_up_to_sign().0;
        assert_eq!(
            exact_signature
                .loop_signature
                .iter()
                .filter(|coefficient| **coefficient != 0)
                .count(),
            1,
            "the UV Taylor rewrite retains only the hard loop direction",
        );
        let routing_sign =
            source.exact_source_routing_sign(hard_soft_owner, true, &exact_signature)?;
        let soft_rows = cograph_edges
            .into_iter()
            .map(|edge| source.outer_loop_signature(edge))
            .collect::<three_dimensional_reps::graph_io::Result<Vec<_>>>()?;
        let soft_rank = exact_integer_rank(&soft_rows);
        let source_row = source.outer_loop_signature(hard_soft_owner)?;
        let residual = |sign: i32| {
            exact_signature
                .loop_signature
                .iter()
                .zip(&source_row)
                .map(|(exact, original)| exact - sign * original)
                .collect::<Vec<_>>()
        };
        let mut compatible = soft_rows.clone();
        compatible.push(residual(routing_sign));
        assert_eq!(exact_integer_rank(&compatible), soft_rank);
        let mut incompatible = soft_rows;
        incompatible.push(residual(-routing_sign));
        assert!(exact_integer_rank(&incompatible) > soft_rank);
        assert!(
            three_dimensional_reps::validate_parsed_graph(&source.to_three_d_parsed_graph()?).ok,
            "routing in the quotient must preserve the source-minor endpoint balance",
        );
        Ok(())
    }

    #[test]
    fn exact_uv_triangle_reuses_only_three_and_four_edge_cff_topologies() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_triangle_cache {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v5 [id=0]
            v0 -> v2 [id=1 lmb_id=0]
            v3 -> v0 [id=2]
            v0 -> v5 [id=3]
            v2 -> v1 [id=4]
            v1 -> v3 [id=5 lmb_id=1]
            v1 -> v4 [id=6]
            v2 -> v3 [id=7]
            v4 -> outgoing [id=8]
        })?;
        let owners = [EdgeIndex(4), EdgeIndex(5), EdgeIndex(7)];
        let carrier = EdgeIndex(5);
        let momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(carrier))
            .finish();
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let base = owners
            .iter()
            .copied()
            .map(|source_edge| FourDDenominator {
                source_edge,
                momentum: if source_edge == EdgeIndex(7) {
                    -momentum.clone()
                } else {
                    momentum.clone()
                },
                mass_squared: mass_squared.clone(),
                full_expr: Atom::one(),
            })
            .collect::<Vec<_>>();
        let mut terms = vec![base.clone()];
        terms.extend((0..owners.len()).map(|dotted_owner| {
            let mut denominators = base.clone();
            denominators.push(base[dotted_owner].clone());
            denominators
        }));

        let options = graph.denominator_only_cff_3d_expression_options();
        let mut cache = crate::cff::generation::ExactCffGenerationCache::default();
        let mut edge_counts = Vec::new();
        for denominators in &terms {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                &graph,
                denominators,
                owners,
            )?;
            graph.register_3d_expression_for_4d_term(
                &source,
                &options,
                &Atom::one(),
                &mut cache,
            )?;
        }
        for denominators in terms {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                &graph,
                &denominators,
                owners,
            )?;
            edge_counts.push(source.to_three_d_parsed_graph()?.internal_edges.len());
            graph.generate_3d_expression_for_4d_term(
                &source,
                &options,
                &Atom::one(),
                Some(&mut cache),
            )?;
        }

        assert_eq!(edge_counts, vec![3, 4, 4, 4]);
        assert_eq!(
            cache.len(),
            2,
            "the numerator hit and all three denominator hits require only the three-edge and four-edge vacuum CFF topologies"
        );
        Ok(())
    }

    #[test]
    fn exact_uv_triangle_matches_direct_factorized_minkowski_cff_variants() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_triangle_source {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v1 [id=0]
            v0 -> v3 [id=1 lmb_id=0]
            v0 -> v5 [id=2]
            v1 -> v2 [id=3 lmb_id=1]
            v1 -> v3 [id=4]
            v2 -> v4 [id=5 lmb_id=2]
            v2 -> v4 [id=6]
            v3 -> v5 [id=7]
            v4 -> v5 [id=8]
            v0 -> outgoing [id=9]
        })?;
        let owners = [EdgeIndex(1), EdgeIndex(2), EdgeIndex(7)];
        let uv_filter = owners
            .into_iter()
            .map(|edge| graph.get_edge_subgraph(edge))
            .reduce(|left, right| left.union(&right))
            .expect("the diagnostic UV triangle has three edges");
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let sub_lmb = graph.try_compatible_sub_lmb(
            &uv_subgraph,
            graph.dummy_stripped_external_flows_of(&uv_subgraph),
            &graph.loop_momentum_basis,
        )?;
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let denominators = owners.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: sub_lmb.loop_atom::<Atom>(source_edge, GS.emr_mom, &[], true),
            mass_squared: mass_squared.clone(),
            full_expr: Atom::one(),
        });
        let exact_source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
            &graph,
            &denominators,
            owners,
            [],
            &sub_lmb,
            ExactUvSubLmbFrame::TaylorVacuum,
        )?;
        let mut direct_graph: Graph = dot!(digraph direct_physical_uv_triangle {
            edge [num=1 mass=1]
            node [num=1]

            v0 -> v3 [id=0 lmb_id=0]
            v0 -> v5 [id=1]
            v3 -> v5 [id=2]
        })?;
        let direct_owners = [EdgeIndex(0), EdgeIndex(1), EdgeIndex(2)];
        let direct_parsed =
            GraphThreeDSource::new(&direct_graph, &[])?.to_three_d_parsed_graph()?;
        let owner_skeleton = exact_source.physical_owner_skeleton()?;
        let exact_parsed = exact_source.to_three_d_parsed_graph()?;
        let incidence = |parsed: &ParsedGraph| {
            parsed
                .internal_edges
                .iter()
                .map(|edge| {
                    if edge.tail <= edge.head {
                        (edge.tail, edge.head)
                    } else {
                        (edge.head, edge.tail)
                    }
                })
                .sorted()
                .collect::<Vec<_>>()
        };
        assert_eq!(incidence(&direct_parsed), incidence(&owner_skeleton));
        assert_eq!(owner_skeleton, exact_parsed);
        assert_eq!(
            direct_parsed
                .internal_edges
                .iter()
                .map(|edge| edge.signature.loop_signature.clone())
                .collect::<Vec<_>>(),
            vec![vec![1], vec![-1], vec![1]],
            "the direct physical triangle carries q1=+Q, q2=-Q, q7=+Q",
        );

        let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
        let direct_numerator =
            GS.emr_mom(direct_owners[0], &mink_index) * GS.emr_mom(direct_owners[1], &mink_index);
        let tagged_momentum = |owner: EdgeIndex, hard: &Atom, index: &Atom| {
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(usize::from(owner) as i64).as_view(),
                    false,
                    hard.as_view(),
                ))
                .add_arg(index.as_view())
                .finish()
        };
        let exact_factors = [
            tagged_momentum(owners[0], &denominators[0].momentum, &mink_index),
            tagged_momentum(owners[1], &denominators[1].momentum, &mink_index),
        ];
        let exact_numerator = &exact_factors[0] * &exact_factors[1];
        let options = graph.denominator_only_cff_3d_expression_options();
        let direct_canonization =
            direct_graph.get_esurface_canonization(&direct_graph.loop_momentum_basis);
        let direct = direct_graph.generate_3d_expression_for_integrand(
            &[],
            &direct_canonization,
            &options,
            Some(&direct_numerator),
        )?;
        let (exact, mapper, plan, _) = graph.generate_3d_expression_for_4d_term(
            &exact_source,
            &options,
            &exact_numerator,
            None,
        )?;
        let owner_occurrences = owners.map(|owner| &mapper.source_edge_occurrences[&owner][0]);
        assert_eq!(
            owner_occurrences.map(|occurrence| occurrence.raw_momentum.clone()),
            std::array::from_fn(|index| denominators[index].momentum.clone()),
        );
        assert_eq!(
            owner_occurrences.map(|occurrence| occurrence.raw_to_parsed_sign),
            [1, -1, 1],
        );
        assert_eq!(direct.source_energy_degree_bounds, vec![(0, 1), (1, 1)]);
        assert_eq!(
            plan.energy_degree_bounds(),
            &[
                (owner_occurrences[0].energy_edge_id, 1),
                (owner_occurrences[1].energy_edge_id, 1),
            ],
        );

        let direct_energy_relabels = direct_owners
            .iter()
            .zip(&owner_occurrences)
            .map(|(direct_owner, occurrence)| {
                Replacement::new(
                    crate::utils::ose_atom_from_index(*direct_owner).to_pattern(),
                    crate::utils::ose_atom_from_index(EdgeIndex(occurrence.energy_edge_id))
                        .to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let direct_spatial_relabels = direct_owners[..2]
            .iter()
            .zip(&owners[..2])
            .map(|(direct_owner, physical_owner)| {
                Replacement::new(
                    GS.emr_vec_index(*direct_owner, &mink_index).to_pattern(),
                    GS.emr_vec_index(*physical_owner, &mink_index).to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let spatial_from = GS.emr_vec_index(owners[1], &mink_index);
        let spatial_to = -GS.emr_vec_index(owners[0], &mink_index);
        let direct_surface_replacements = direct.expression.surfaces.get_all_replacements_gs(&[]);
        let exact_surface_replacements = exact.expression.surfaces.get_all_replacements_gs(&[]);
        let mut direct_terms = direct
            .expression
            .orientations
            .iter()
            .map(|orientation| {
                let coefficient = orientation
                    .to_atom_gs()
                    .replace_multiple(&direct_surface_replacements)
                    .replace_multiple(&direct_energy_relabels)
                    .replace_multiple(mapper.exact_ose_replacements())
                    .together();
                let numerator = direct_numerator
                    .replace_multiple(orientation.energy_replacements_gs(&direct_graph))
                    .replace_multiple(&direct_spatial_relabels)
                    .replace_multiple(&direct_energy_relabels)
                    .replace_multiple(mapper.exact_ose_replacements())
                    .replace(spatial_from.clone())
                    .with(spatial_to.clone());
                (coefficient, numerator)
            })
            .collect::<Vec<_>>();
        let mut exact_terms = exact
            .expression
            .orientations
            .iter()
            .map(|orientation| {
                let coefficient = orientation
                    .to_atom_gs()
                    .replace_multiple(&exact_surface_replacements)
                    .replace_multiple(mapper.exact_ose_replacements())
                    .together();
                let numerator = mapper
                    .map_planned_numerator(
                        &orientation.loop_energy_map,
                        &orientation.edge_energy_map,
                        &plan,
                    )?
                    .replace(spatial_from.clone())
                    .with(spatial_to.clone());
                Ok((coefficient, numerator))
            })
            .collect::<Result<Vec<_>>>()?;
        assert_eq!(direct_terms.len(), exact_terms.len());
        for (direct_coefficient, direct_numerator) in direct_terms.drain(..) {
            let Some(position) = exact_terms.iter().position(|(coefficient, numerator)| {
                (&direct_coefficient - coefficient).together().is_zero()
                    && (&direct_numerator - numerator).expand().is_zero()
            }) else {
                panic!(
                    "direct triangle variant has no reconstructed match: coefficient={direct_coefficient}, numerator={direct_numerator}, remaining={exact_terms:?}",
                );
            };
            let _ = exact_terms.swap_remove(position);
        }
        assert!(exact_terms.is_empty());
        Ok(())
    }

    #[test]
    fn exact_sub_lmb_taylor_vacuum_keeps_dotted_bubble_owner_incidence() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_taylor_vacuum_bubble {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b -> c [id=3]
            c -> outgoing [id=4]
        })?;
        let uv_edges = [EdgeIndex(1), EdgeIndex(2)];
        let uv_filter = graph
            .get_edge_subgraph(uv_edges[0])
            .union(&graph.get_edge_subgraph(uv_edges[1]));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let crown = graph.dummy_stripped_external_flows_of(&uv_subgraph);
        let sub_lmb = graph.try_compatible_sub_lmb(
            &uv_subgraph,
            crown.clone(),
            &graph.loop_momentum_basis,
        )?;
        let physical_denominators = uv_edges.map(|source_edge| FourDDenominator {
            source_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(source_edge))
                .finish(),
            mass_squared: graph.underlying[source_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        });
        let physical_source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
            &graph,
            &physical_denominators,
            uv_edges,
            crown.included_iter(),
            &sub_lmb,
            ExactUvSubLmbFrame::RetainedPhysicalCrown,
        )?;
        let physical_parsed = physical_source.to_three_d_parsed_graph()?;
        assert!(!physical_parsed.external_edges.is_empty());
        assert!(physical_parsed.internal_edges.iter().any(|edge| {
            edge.signature
                .external_signature
                .iter()
                .any(|coefficient| *coefficient != 0)
        }));
        let carrier = FunctionBuilder::new(GS.emr_mom).add_arg(1).finish();
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: uv_edges[0],
                momentum: carrier.clone(),
                mass_squared: mass_squared.clone(),
                full_expr: Atom::one(),
            },
            FourDDenominator {
                source_edge: uv_edges[1],
                momentum: -carrier,
                mass_squared,
                full_expr: Atom::one(),
            },
        ];

        assert!(
            GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &denominators,
                uv_edges,
                [],
                &sub_lmb,
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
            )
            .is_err(),
            "the physical-source frame must continue to reject an omitted crown",
        );
        assert!(
            GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                &graph,
                &physical_denominators,
                uv_edges,
                [],
                &sub_lmb,
                ExactUvSubLmbFrame::TaylorVacuum,
            )
            .is_err(),
            "the Taylor-vacuum frame must reject a denominator that retains its physical crown shift",
        );
        let outside_owner = FourDDenominator {
            source_edge: EdgeIndex(3),
            momentum: FunctionBuilder::new(GS.emr_mom).add_arg(3).finish(),
            mass_squared: Atom::var(GS.m_uv_expansion).pow(2),
            full_expr: Atom::one(),
        };
        for (frame, boundaries) in [
            (
                ExactUvSubLmbFrame::RetainedPhysicalCrown,
                crown.included_iter().collect::<Vec<_>>(),
            ),
            (ExactUvSubLmbFrame::TaylorVacuum, Vec::new()),
        ] {
            assert!(
                GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                    &graph,
                    std::slice::from_ref(&outside_owner),
                    uv_edges,
                    boundaries,
                    &sub_lmb,
                    frame,
                )
                .is_err(),
                "both exact sub-LMB frames must reject owners outside the UV source edges",
            );
        }
        let source = GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
            &graph,
            &denominators,
            uv_edges,
            [],
            &sub_lmb,
            ExactUvSubLmbFrame::TaylorVacuum,
        )?;
        let parsed = source.to_three_d_parsed_graph()?;
        assert_eq!(parsed.internal_edges.len(), 2);
        assert!(parsed.external_edges.is_empty());
        assert!(parsed.initial_state_cut_edges.is_empty());
        assert!(parsed.internal_edges.iter().all(|edge| {
            edge.signature
                .external_signature
                .iter()
                .all(|coefficient| *coefficient == 0)
        }));
        assert_eq!(
            source
                .physical_energy_edge_index_map()
                .expect("the vacuum source retains physical owners")
                .internal
                .into_values()
                .collect::<BTreeSet<_>>(),
            uv_edges
                .into_iter()
                .map(usize::from)
                .collect::<BTreeSet<_>>(),
        );
        assert!(validate_parsed_graph(&parsed).ok);

        // The physical owner-2 momentum is q+p in this sub-LMB, whereas the
        // Taylor-vacuum denominator and its odd hard numerator retain only q.
        // The remote crown EMR belongs to the factorized outer graph and must
        // therefore pass through the child energy mapper literally.
        let hard = sub_lmb.loop_atom::<Atom>(uv_edges[1], GS.emr_mom, &[], true);
        let soft = sub_lmb.ext_atom(uv_edges[1], GS.emr_mom, &[W_.x___], true);
        assert!(!soft.is_zero());
        let tag = GS.uv_momentum_provenance_tag(
            Atom::num(usize::from(uv_edges[1]) as i64).as_view(),
            true,
            hard.as_view(),
        );
        let tagged_hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(tag)
            .add_arg(GS.cind(0))
            .finish();
        let crown_energy = GS.emr_mom(EdgeIndex(3), GS.cind(0));
        let analysis_numerator = &tagged_hard * &crown_energy;
        let options = graph.denominator_only_cff_3d_expression_options();
        let (generated, mapper, _, _) = graph.generate_3d_expression_for_4d_term(
            &source,
            &options,
            &analysis_numerator,
            None,
        )?;
        let orientation = generated
            .expression
            .orientations
            .first()
            .expect("the vacuum bubble has a CFF orientation");
        let left_occurrence = &mapper.source_edge_occurrences[&uv_edges[0]][0];
        let right_occurrence = &mapper.source_edge_occurrences[&uv_edges[1]][0];
        assert_eq!(left_occurrence.raw_momentum, denominators[0].momentum);
        assert_eq!(right_occurrence.raw_momentum, denominators[1].momentum);
        assert_eq!(
            [
                left_occurrence.raw_to_parsed_sign,
                right_occurrence.raw_to_parsed_sign,
            ],
            [1, -1],
            "canonicalizing the even D(-Q) denominator must retain the odd raw-to-parsed numerator sign",
        );
        let parsed_signature = |energy_edge_id| {
            let occurrence = source
                .exact_occurrences
                .iter()
                .find(|occurrence| occurrence.energy_edge_id == energy_edge_id)
                .expect("each owned energy edge has one exact occurrence");
            &parsed.internal_edges[occurrence.local_edge_id].signature
        };
        assert_eq!(
            parsed_signature(left_occurrence.energy_edge_id),
            parsed_signature(right_occurrence.energy_edge_id),
            "the source canonicalizer puts the two even bubble denominators in one routed energy frame",
        );
        let parsed_energy = |occurrence: &ExactSourceOwnerOccurrence| {
            orientation.edge_energy_map[occurrence.energy_edge_id]
                .to_atom_gs(&[])
                .replace_multiple(mapper.exact_ose_replacements())
        };
        assert_eq!(
            parsed_energy(right_occurrence),
            -parsed_energy(left_occurrence),
            "the canonically reversed incidence gives this representative opposite oriented parsed-edge samples before raw projection",
        );
        assert_ne!(
            parsed_energy(left_occurrence),
            parsed_energy(right_occurrence),
            "the first causal representative must retain its nonzero incidence-relative orientation signs",
        );
        let abstract_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::num(2)]);
        let tagged_fixed_hard = |owner: EdgeIndex, hard: &Atom, index: &Atom| {
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(usize::from(owner) as i64).as_view(),
                    false,
                    hard.as_view(),
                ))
                .add_arg(index.as_view())
                .finish()
        };
        let left_fixed = tagged_fixed_hard(uv_edges[0], &denominators[0].momentum, &abstract_index);
        let right_fixed =
            tagged_fixed_hard(uv_edges[1], &denominators[1].momentum, &abstract_index);
        let left_temporal = tagged_fixed_hard(uv_edges[0], &denominators[0].momentum, &GS.cind(0));
        let right_temporal = tagged_fixed_hard(uv_edges[1], &denominators[1].momentum, &GS.cind(0));
        let fixed_pair = &left_fixed * &right_fixed;
        let (pair_generated, pair_mapper, pair_plan, _) =
            graph.generate_3d_expression_for_4d_term(&source, &options, &fixed_pair, None)?;
        assert_eq!(
            pair_plan.energy_degree_bounds(),
            &[
                (left_occurrence.energy_edge_id, 1),
                (right_occurrence.energy_edge_id, 1),
            ],
            "the production analyzer must keep the two fixed owner factors on their own exact occurrences",
        );
        assert!(!pair_generated.expression.orientations.is_empty());
        let positive_right_energy =
            crate::utils::ose_atom_from_index(EdgeIndex(right_occurrence.energy_edge_id))
                .replace_multiple(pair_mapper.exact_ose_replacements());
        let sample_sign = |sample: &Atom| {
            [-1, 0, 1].into_iter().find(|sign| {
                (sample - Atom::num(*sign) * &positive_right_energy)
                    .expand()
                    .is_zero()
            })
        };
        let mut contact_samples = BTreeSet::new();
        for (orientation_id, pair_orientation) in
            pair_generated.expression.orientations.iter().enumerate()
        {
            let pair_mapped_left = pair_mapper.map_numerator(
                &pair_orientation.loop_energy_map,
                &pair_orientation.edge_energy_map,
                &left_fixed,
            )?;
            let pair_mapped_right = pair_mapper.map_numerator(
                &pair_orientation.loop_energy_map,
                &pair_orientation.edge_energy_map,
                &right_fixed,
            )?;
            let planned_pair = pair_mapper.map_planned_numerator(
                &pair_orientation.loop_energy_map,
                &pair_orientation.edge_energy_map,
                &pair_plan,
            )?;
            assert_eq!(
                planned_pair,
                &pair_mapped_left * &pair_mapped_right,
                "the production assignment plan must recombine both mapped factors unchanged in orientation {orientation_id}",
            );
            let mapped_left_temporal = pair_mapper.map_numerator(
                &pair_orientation.loop_energy_map,
                &pair_orientation.edge_energy_map,
                &left_temporal,
            )?;
            let mapped_right_temporal = pair_mapper.map_numerator(
                &pair_orientation.loop_energy_map,
                &pair_orientation.edge_energy_map,
                &right_temporal,
            )?;
            let is_contact = pair_orientation.variants.iter().any(|variant| {
                variant.origin.as_deref().is_some_and(|origin| {
                    origin.starts_with("bounded_degree_quadratic_recursive_contact:")
                })
            });
            if !is_contact {
                assert!(
                    !mapped_left_temporal.is_zero()
                        && !mapped_right_temporal.is_zero()
                        && (&mapped_left_temporal + &mapped_right_temporal)
                            .expand()
                            .is_zero(),
                    "nonzero fixed +Q and -Q owner samples must remain opposite in orientation {orientation_id}: left={mapped_left_temporal}, right={mapped_right_temporal}, directions={:?}, label={:?}, origins={:?}",
                    pair_orientation.data.orientation,
                    pair_orientation.data.label,
                    pair_orientation
                        .variants
                        .iter()
                        .filter_map(|variant| variant.origin.as_deref())
                        .collect::<Vec<_>>(),
                );
            } else {
                contact_samples.insert((
                    sample_sign(&mapped_left_temporal).ok_or_else(|| {
                        eyre::eyre!(
                            "left contact sample is not -E, 0, or +E in orientation {orientation_id}: {mapped_left_temporal}",
                        )
                    })?,
                    sample_sign(&mapped_right_temporal).ok_or_else(|| {
                        eyre::eyre!(
                            "right contact sample is not -E, 0, or +E in orientation {orientation_id}: {mapped_right_temporal}",
                        )
                    })?,
                ));
            }
        }
        assert_eq!(
            contact_samples,
            BTreeSet::from([(-1, 1), (0, 1), (1, 1)]),
            "the fixed-owner bubble must retain the deterministic e0=-,0,+ contact family while the untouched right owner stays at +E, including the (0,+E) lower-sector map",
        );

        // Canonicalizing the even D(-Q) denominator may change its oriented
        // parsed-edge contact sample. A fixed factor stays on its canonical
        // owner occurrence, whereas a denominator-derived factor may use a
        // serial copy and therefore composes that copy's routing sign. The
        // dedicated mapper-role regression exercises the latter dispatch.
        let mapped_hard = mapper.map_numerator(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            &tagged_hard,
        )?;
        let mapped_product = mapper.map_numerator(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            &analysis_numerator,
        )?;
        assert_eq!(mapped_product, mapped_hard * crown_energy);

        let remote_owner = EdgeIndex(3);
        let remote_hard = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(remote_owner))
            .finish();
        let remote_tag = GS.uv_momentum_provenance_tag(
            Atom::num(usize::from(remote_owner) as i64).as_view(),
            false,
            remote_hard.as_view(),
        );
        let mink_index = LibraryRep::from(Minkowski {}).to_symbolic([Atom::var(W_.a__)]);
        for index in [GS.cind(0), mink_index] {
            let remote_tagged = FunctionBuilder::new(GS.emr_mom)
                .add_arg(remote_tag.clone())
                .add_arg(index)
                .finish();
            assert_eq!(
                mapper.map_numerator(
                    &orientation.loop_energy_map,
                    &orientation.edge_energy_map,
                    &remote_tagged,
                )?,
                remote_tagged,
                "a remote temporal or abstract provenance tag must survive for its owning component"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_uv_source_retains_a_non_vacuum_two_point_shift() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_two_point {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b -> outgoing [id=3]
        })?;
        let loop_carrier = FunctionBuilder::new(GS.emr_mom).add_arg(1).finish();
        let mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let denominators = [
            FourDDenominator {
                source_edge: EdgeIndex(1),
                momentum: loop_carrier.clone(),
                mass_squared: mass_squared.clone(),
                full_expr: Atom::one(),
            },
            FourDDenominator {
                source_edge: EdgeIndex(2),
                momentum: FunctionBuilder::new(GS.emr_mom).add_arg(2).finish(),
                mass_squared,
                full_expr: Atom::one(),
            },
        ];
        let uv_boundary_hedges = [EdgeIndex(0), EdgeIndex(3)].map(|edge_id| {
            let (_, HedgePair::Unpaired { hedge, .. }) = graph[&edge_id] else {
                panic!("two-point boundary edge must be external")
            };
            hedge
        });
        let parsed = GraphThreeDSource::from_exact_denominators_in_uv_edges_and_boundaries(
            &graph,
            &denominators,
            [EdgeIndex(1), EdgeIndex(2)],
            uv_boundary_hedges,
        )?
        .to_three_d_parsed_graph()?;

        assert_eq!(parsed.internal_edges.len(), 2);
        assert_eq!(
            parsed
                .internal_edges
                .iter()
                .flat_map(|edge| [edge.tail, edge.head])
                .collect::<BTreeSet<_>>()
                .len(),
            2
        );
        assert_eq!(parsed.external_edges.len(), 2);
        assert!(parsed.external_edges.iter().all(|edge| {
            edge.external_coefficients
                .iter()
                .any(|coefficient| *coefficient != 0)
        }));
        let validation = three_dimensional_reps::validate_parsed_graph(&parsed);
        assert!(validation.ok);
        assert!(validation.vertex_external_balance_info.is_empty());
        Ok(())
    }

    #[test]
    fn exact_source_separates_loop_denominators_from_external_tree_factors() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph shifted_loop {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> a [id=4]
            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
            c -> d [id=3]
            d -> outgoing [id=5]
        })?;
        let loop_edge = graph.loop_momentum_basis.loop_edges[LoopIndex::from(0)];
        let tree_edge = EdgeIndex::from(3);
        let loop_denominator = FourDDenominator {
            source_edge: loop_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(loop_edge))
                .finish(),
            mass_squared: graph.underlying[loop_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };
        let tree_denominator = FourDDenominator {
            source_edge: tree_edge,
            momentum: FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(tree_edge))
                .finish(),
            mass_squared: graph.underlying[tree_edge].particle.mass_atom().pow(2),
            full_expr: Atom::one(),
        };

        assert!(loop_denominator.depends_on_loop(&graph, false)?);
        assert!(!tree_denominator.depends_on_loop(&graph, false)?);

        let external_edge = *graph
            .loop_momentum_basis
            .ext_edges
            .first()
            .expect("shifted loop graph has an external momentum carrier");
        let shifted = FourDDenominator {
            source_edge: loop_edge,
            momentum: &loop_denominator.momentum
                + FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(external_edge))
                    .finish(),
            mass_squared: Atom::var(GS.m_uv_expansion).pow(2),
            full_expr: Atom::one(),
        };
        let shifted_signature = shifted.momentum_signature(&graph, false)?;
        assert_eq!(shifted_signature.loop_signature, vec![1]);
        assert_eq!(
            shifted_signature.external_signature,
            (&graph.loop_momentum_basis.edge_signatures[external_edge].external)
                .into_iter()
                .map(sign_to_i32)
                .collect::<Vec<_>>()
        );
        Ok(())
    }

    #[test]
    fn pure_tree_exact_denominators_have_no_active_loop_source() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0]
            b -> c [id=1]
        })?;
        for edge in [EdgeIndex::from(0), EdgeIndex::from(1)] {
            let denominator = FourDDenominator {
                source_edge: edge,
                momentum: FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(edge))
                    .finish(),
                mass_squared: graph.underlying[edge].particle.mass_atom().pow(2),
                full_expr: Atom::one(),
            };
            assert!(!denominator.depends_on_loop(&graph, false)?);
        }
        Ok(())
    }
}
