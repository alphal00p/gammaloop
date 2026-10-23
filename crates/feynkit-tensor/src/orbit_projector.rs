//! Exact projector coefficients between repeated-vector contraction orbits.
//!
//! Let `G` be the Gram matrix on labeled perfect matchings,
//! `G[P,Q] = D^cycles(P,Q)`, and let `W = G^-1` be the orthogonal-Weingarten
//! matrix.  If `S_A` and `S_B` are the raw zero-or-one incidence matrices from
//! labeled matchings to internal and projector color orbits, respectively,
//! tensor reduction needs
//!
//! `C = S_A^T W S_B`.
//!
//! The full matching space is never constructed here.  Color symmetry implies
//! `G S_A = S_A H_A` for a small orbit Gram matrix `H_A`, and likewise for
//! `B`.  Consequently,
//!
//! `C = H_A^-T N = N H_B^-1`, where `N = S_A^T S_B`.
//!
//! We construct `N` by memoized matching of joint-color multiplicities and
//! construct `H` by memoized alternating cycles around a representative
//! matching.  The smaller of the two orbit systems is solved.  This is exactly
//! the same double-orbit sum of orthogonal-Weingarten weights, without the
//! factorial-squared labeled-pairing loop.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::mem;

use symbolica::atom::{Atom, AtomCore};
use symbolica::solve::SolveError;
use symbolica::symbol;
use symbolica::utils::resources::{self, ResourceLimitError};
use thiserror::Error;

/// One symmetry-reduced coefficient in a fully contracted tensor projector.
#[derive(Clone, Debug)]
pub(crate) struct OrbitProjectorTerm {
    internal_alpha: Box<[u16]>,
    projector_alpha: Box<[u16]>,
    internal_labeled_pairings: u128,
    projector_labeled_pairings: u128,
    coefficient: Atom,
}

impl OrbitProjectorTerm {
    /// Upper-triangular contraction exponents for integrated-vector classes.
    pub(crate) fn internal_alpha(&self) -> &[u16] {
        &self.internal_alpha
    }

    /// Upper-triangular contraction exponents for projector-vector classes.
    pub(crate) fn projector_alpha(&self) -> &[u16] {
        &self.projector_alpha
    }

    /// Number of labeled internal matchings in this orbit.
    pub(crate) fn internal_labeled_pairings(&self) -> u128 {
        self.internal_labeled_pairings
    }

    /// Number of labeled projector matchings in this orbit.
    pub(crate) fn projector_labeled_pairings(&self) -> u128 {
        self.projector_labeled_pairings
    }

    /// Exact sum of orthogonal-Weingarten weights between both orbits.
    pub(crate) fn coefficient(&self) -> &Atom {
        &self.coefficient
    }
}

/// Exact double-orbit projector built without labeled matching products.
#[derive(Clone, Debug)]
pub(crate) struct OrbitProjector {
    terms: Vec<OrbitProjectorTerm>,
}

impl OrbitProjector {
    /// Build coefficients for two colorings of the same tensor indices.
    ///
    /// Equal labels identify repeated vectors.  Labels need not be dense; they
    /// are canonicalized independently on both sides. An explicit `output_limit`
    /// bounds intermediate and final orbit maps, not labeled matchings.
    pub(crate) fn new(
        internal_colors: &[usize],
        projector_colors: &[usize],
        dimension: Atom,
        output_limit: Option<usize>,
    ) -> Result<Self, OrbitProjectorError> {
        resources::checkpoint()?;
        if internal_colors.len() != projector_colors.len() {
            return Err(OrbitProjectorError::ColorLengthMismatch {
                internal: internal_colors.len(),
                projector: projector_colors.len(),
            });
        }
        let rank = internal_colors.len();
        if !rank.is_multiple_of(2) {
            return Err(OrbitProjectorError::OddRank(rank));
        }
        let (internal_colors, internal_class_count) = dense_colors(internal_colors);
        let (projector_colors, projector_class_count) = dense_colors(projector_colors);
        resources::checkpoint()?;
        let intersections = MatchingIntersections::new(
            &internal_colors,
            internal_class_count,
            &projector_colors,
            projector_class_count,
            output_limit,
        )?;
        let internal_alphas = intersections
            .counts
            .keys()
            .map(|key| key.internal.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let projector_alphas = intersections
            .counts
            .keys()
            .map(|key| key.projector.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let term_bound = internal_alphas
            .len()
            .checked_mul(projector_alphas.len())
            .ok_or(OrbitProjectorError::CounterOverflow)?;
        if let Some(limit) = output_limit
            && term_bound > limit
        {
            return Err(OrbitProjectorError::OutputLimit {
                terms: term_bound,
                limit,
            });
        }
        let matrix_bytes = term_bound
            .checked_mul(mem::size_of::<u128>() + mem::size_of::<Atom>())
            .and_then(|bytes| u64::try_from(bytes).ok())
            .ok_or(OrbitProjectorError::CounterOverflow)?;
        resources::check_allocation(matrix_bytes)?;
        let internal_positions = orbit_positions(&internal_alphas);
        let projector_positions = orbit_positions(&projector_alphas);
        let mut incidence = vec![vec![0_u128; projector_alphas.len()]; internal_alphas.len()];
        for (key, count) in intersections.counts {
            resources::checkpoint()?;
            incidence[internal_positions[&key.internal]][projector_positions[&key.projector]] =
                count;
        }
        let internal_sizes = incidence
            .iter()
            .map(|row| {
                row.iter().try_fold(0_u128, |sum, count| {
                    resources::checkpoint()?;
                    sum.checked_add(*count)
                        .ok_or(OrbitProjectorError::CounterOverflow)
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let projector_sizes = (0..projector_alphas.len())
            .map(|column| {
                incidence.iter().try_fold(0_u128, |sum, row| {
                    resources::checkpoint()?;
                    sum.checked_add(row[column])
                        .ok_or(OrbitProjectorError::CounterOverflow)
                })
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut coefficients =
            vec![vec![Atom::Zero; projector_alphas.len()]; internal_alphas.len()];
        if internal_alphas.len() <= projector_alphas.len() {
            let gram = orbit_gram(
                &internal_alphas,
                internal_class_count,
                &dimension,
                output_limit,
            )?;
            for projector in 0..projector_alphas.len() {
                resources::checkpoint()?;
                let right_hand_side = incidence
                    .iter()
                    .map(|row| atom_from_count(row[projector]))
                    .collect::<Result<Vec<_>, _>>()?;
                let column = solve_transposed(&gram, &right_hand_side, rank)?;
                for (internal, coefficient) in column.into_iter().enumerate() {
                    resources::checkpoint()?;
                    coefficients[internal][projector] = coefficient.cancel();
                    resources::checkpoint()?;
                }
            }
        } else {
            let gram = orbit_gram(
                &projector_alphas,
                projector_class_count,
                &dimension,
                output_limit,
            )?;
            for internal in 0..internal_alphas.len() {
                resources::checkpoint()?;
                let right_hand_side = incidence[internal]
                    .iter()
                    .copied()
                    .map(atom_from_count)
                    .collect::<Result<Vec<_>, _>>()?;
                let row = solve_transposed(&gram, &right_hand_side, rank)?;
                for (projector, coefficient) in row.into_iter().enumerate() {
                    resources::checkpoint()?;
                    coefficients[internal][projector] = coefficient.cancel();
                    resources::checkpoint()?;
                }
            }
        }

        let mut terms = Vec::with_capacity(term_bound);
        for (internal, internal_alpha) in internal_alphas.iter().enumerate() {
            for (projector, projector_alpha) in projector_alphas.iter().enumerate() {
                resources::checkpoint()?;
                let coefficient = coefficients[internal][projector].clone();
                if coefficient.is_zero() {
                    continue;
                }
                terms.push(OrbitProjectorTerm {
                    internal_alpha: internal_alpha.clone().into_boxed_slice(),
                    projector_alpha: projector_alpha.clone().into_boxed_slice(),
                    internal_labeled_pairings: internal_sizes[internal],
                    projector_labeled_pairings: projector_sizes[projector],
                    coefficient,
                });
            }
        }

        resources::checkpoint()?;
        Ok(Self { terms })
    }

    pub(crate) fn terms(&self) -> &[OrbitProjectorTerm] {
        &self.terms
    }
}

/// Failures detected before a symmetry-reduced projector can be constructed.
#[derive(Debug, Error)]
pub(crate) enum OrbitProjectorError {
    #[error("{0}")]
    Resource(#[from] ResourceLimitError),
    #[error("internal colors have rank {internal}, but projector colors have rank {projector}")]
    ColorLengthMismatch { internal: usize, projector: usize },
    #[error("tensor rank {0} is odd; vacuum tensor reduction requires an even rank")]
    OddRank(usize),
    #[error("symmetry reduction produced {terms} orbit terms, above limit {limit}")]
    OutputLimit { terms: usize, limit: usize },
    #[error("orbit multiplicity or representation counter overflowed")]
    CounterOverflow,
    #[error("could not solve the rank-{rank} symmetry-reduced projector: {source}")]
    Solve {
        rank: usize,
        #[source]
        source: SolveError,
    },
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct OrbitPairKey {
    internal: Vec<u16>,
    projector: Vec<u16>,
}

struct MatchingIntersections {
    counts: BTreeMap<OrbitPairKey, u128>,
}

impl MatchingIntersections {
    fn new(
        internal: &[usize],
        internal_class_count: usize,
        projector: &[usize],
        projector_class_count: usize,
        limit: Option<usize>,
    ) -> Result<Self, OrbitProjectorError> {
        resources::checkpoint()?;
        let joint_len = internal_class_count
            .checked_mul(projector_class_count)
            .ok_or(OrbitProjectorError::CounterOverflow)?;
        resources::check_allocation(
            u64::try_from(joint_len).map_err(|_| OrbitProjectorError::CounterOverflow)?,
        )?;
        let mut joint_counts = vec![0_u8; joint_len];
        for (&internal, &projector) in internal.iter().zip(projector) {
            resources::checkpoint()?;
            let count = &mut joint_counts[internal * projector_class_count + projector];
            *count = count
                .checked_add(1)
                .ok_or(OrbitProjectorError::CounterOverflow)?;
        }
        let alpha_lengths = (
            triangle_len(internal_class_count)?,
            triangle_len(projector_class_count)?,
        );
        let mut cache = HashMap::new();
        let counts = matching_intersections(
            &joint_counts,
            internal_class_count,
            projector_class_count,
            alpha_lengths,
            limit,
            &mut cache,
        )?;
        Ok(Self { counts })
    }
}

fn matching_intersections(
    counts: &[u8],
    internal_class_count: usize,
    projector_class_count: usize,
    alpha_lengths: (usize, usize),
    limit: Option<usize>,
    cache: &mut HashMap<Vec<u8>, BTreeMap<OrbitPairKey, u128>>,
) -> Result<BTreeMap<OrbitPairKey, u128>, OrbitProjectorError> {
    resources::checkpoint()?;
    if let Some(cached) = cache.get(counts) {
        return Ok(cached.clone());
    }
    let Some(anchor) = counts.iter().position(|count| *count != 0) else {
        return Ok(BTreeMap::from([(
            OrbitPairKey {
                internal: vec![0; alpha_lengths.0],
                projector: vec![0; alpha_lengths.1],
            },
            1,
        )]));
    };

    let mut after_anchor = counts.to_vec();
    after_anchor[anchor] -= 1;
    let (anchor_internal, anchor_projector) = split_joint(anchor, projector_class_count);
    let mut result = BTreeMap::<OrbitPairKey, u128>::new();
    for partner in 0..after_anchor.len() {
        resources::checkpoint()?;
        let choices = after_anchor[partner];
        if choices == 0 {
            continue;
        }
        let mut remaining = after_anchor.clone();
        remaining[partner] -= 1;
        let (partner_internal, partner_projector) = split_joint(partner, projector_class_count);
        for (mut key, multiplicity) in matching_intersections(
            &remaining,
            internal_class_count,
            projector_class_count,
            alpha_lengths,
            limit,
            cache,
        )? {
            let internal = &mut key.internal
                [upper_triangle_index(anchor_internal, partner_internal, internal_class_count)];
            *internal = internal
                .checked_add(1)
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            let projector = &mut key.projector
                [upper_triangle_index(anchor_projector, partner_projector, projector_class_count)];
            *projector = projector
                .checked_add(1)
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            add_count(
                &mut result,
                key,
                multiplicity
                    .checked_mul(u128::from(choices))
                    .ok_or(OrbitProjectorError::CounterOverflow)?,
                limit,
            )?;
        }
    }
    cache.insert(counts.to_vec(), result.clone());
    Ok(result)
}

fn orbit_gram(
    alphas: &[Vec<u16>],
    class_count: usize,
    dimension: &Atom,
    limit: Option<usize>,
) -> Result<Vec<Vec<Atom>>, OrbitProjectorError> {
    resources::checkpoint()?;
    let matrix_bytes = alphas
        .len()
        .checked_mul(alphas.len())
        .and_then(|cells| cells.checked_mul(mem::size_of::<Atom>()))
        .and_then(|bytes| u64::try_from(bytes).ok())
        .ok_or(OrbitProjectorError::CounterOverflow)?;
    resources::check_allocation(matrix_bytes)?;
    let positions = orbit_positions(alphas);
    let mut gram = Vec::with_capacity(alphas.len());
    for representative in alphas {
        resources::checkpoint()?;
        let distribution = relative_matching_distribution(representative, class_count, limit)?;
        let mut row = vec![Atom::Zero; alphas.len()];
        for (key, multiplicity) in distribution {
            resources::checkpoint()?;
            let target = positions[&key.alpha];
            row[target] += atom_from_count(multiplicity)?
                * dimension.clone().pow(Atom::num(i64::from(key.cycles)));
            resources::checkpoint()?;
        }
        gram.push(row);
    }
    Ok(gram)
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct RelativeKey {
    alpha: Vec<u16>,
    cycles: u8,
}

fn relative_matching_distribution(
    representative: &[u16],
    class_count: usize,
    limit: Option<usize>,
) -> Result<BTreeMap<RelativeKey, u128>, OrbitProjectorError> {
    resources::checkpoint()?;
    let edge_types = upper_triangle_pairs(class_count)?;
    let counts = representative
        .iter()
        .copied()
        .map(|count| u8::try_from(count).map_err(|_| OrbitProjectorError::CounterOverflow))
        .collect::<Result<Vec<_>, _>>()?;
    let mut cache = HashMap::new();
    relative_components(
        &counts,
        &edge_types,
        triangle_len(class_count)?,
        limit,
        &mut cache,
    )
}

fn relative_components(
    counts: &[u8],
    edge_types: &[(usize, usize)],
    alpha_len: usize,
    limit: Option<usize>,
    cache: &mut HashMap<Vec<u8>, BTreeMap<RelativeKey, u128>>,
) -> Result<BTreeMap<RelativeKey, u128>, OrbitProjectorError> {
    resources::checkpoint()?;
    if let Some(cached) = cache.get(counts) {
        return Ok(cached.clone());
    }
    let Some(anchor_type) = counts.iter().position(|count| *count != 0) else {
        return Ok(BTreeMap::from([(
            RelativeKey {
                alpha: vec![0; alpha_len],
                cycles: 0,
            },
            1,
        )]));
    };

    let mut after_anchor = counts.to_vec();
    after_anchor[anchor_type] -= 1;
    let (anchor_entrance, anchor_exit) = edge_types[anchor_type];
    let mut component_options = BTreeMap::new();

    // A one-edge alternating cycle is the matching edge shared with the
    // representative.  Unlike longer cycles it has only one orientation.
    let mut shared_alpha = vec![0_u16; alpha_len];
    shared_alpha
        [upper_triangle_index(anchor_entrance, anchor_exit, edge_color_count(edge_types))] = 1;
    add_component_option(
        &mut component_options,
        after_anchor.clone(),
        shared_alpha,
        1,
        limit,
    )?;

    extend_component(
        &after_anchor,
        edge_types,
        anchor_entrance,
        anchor_exit,
        vec![0; alpha_len],
        1,
        limit,
        &mut component_options,
    )?;

    let mut result = BTreeMap::new();
    for (component, component_multiplicity) in component_options {
        resources::checkpoint()?;
        for (mut rest, rest_multiplicity) in
            relative_components(&component.remaining, edge_types, alpha_len, limit, cache)?
        {
            add_alpha(&mut rest.alpha, &component.alpha)?;
            rest.cycles = rest
                .cycles
                .checked_add(1)
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            let multiplicity = component_multiplicity
                .checked_mul(rest_multiplicity)
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            add_count(&mut result, rest, multiplicity, limit)?;
        }
    }
    cache.insert(counts.to_vec(), result.clone());
    Ok(result)
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct ComponentKey {
    remaining: Vec<u8>,
    alpha: Vec<u16>,
}

#[allow(clippy::too_many_arguments)]
fn extend_component(
    remaining: &[u8],
    edge_types: &[(usize, usize)],
    anchor_entrance: usize,
    last_exit: usize,
    alpha: Vec<u16>,
    path_multiplicity: u128,
    limit: Option<usize>,
    output: &mut BTreeMap<ComponentKey, u128>,
) -> Result<(), OrbitProjectorError> {
    resources::checkpoint()?;
    let class_count = edge_color_count(edge_types);
    for next_type in 0..remaining.len() {
        resources::checkpoint()?;
        let choices = remaining[next_type];
        if choices == 0 {
            continue;
        }
        let mut after_next = remaining.to_vec();
        after_next[next_type] -= 1;
        let (left, right) = edge_types[next_type];
        for (entrance, exit) in [(left, right), (right, left)] {
            resources::checkpoint()?;
            let multiplicity = path_multiplicity
                .checked_mul(u128::from(choices))
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            let mut extended_alpha = alpha.clone();
            let extended =
                &mut extended_alpha[upper_triangle_index(last_exit, entrance, class_count)];
            *extended = extended
                .checked_add(1)
                .ok_or(OrbitProjectorError::CounterOverflow)?;

            let mut closed_alpha = extended_alpha.clone();
            let closed =
                &mut closed_alpha[upper_triangle_index(exit, anchor_entrance, class_count)];
            *closed = closed
                .checked_add(1)
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            add_component_option(
                output,
                after_next.clone(),
                closed_alpha,
                multiplicity,
                limit,
            )?;
            extend_component(
                &after_next,
                edge_types,
                anchor_entrance,
                exit,
                extended_alpha,
                multiplicity,
                limit,
                output,
            )?;
        }
    }
    Ok(())
}

fn add_component_option(
    output: &mut BTreeMap<ComponentKey, u128>,
    remaining: Vec<u8>,
    alpha: Vec<u16>,
    multiplicity: u128,
    limit: Option<usize>,
) -> Result<(), OrbitProjectorError> {
    add_count(
        output,
        ComponentKey { remaining, alpha },
        multiplicity,
        limit,
    )
}

fn solve_transposed(
    gram: &[Vec<Atom>],
    right_hand_side: &[Atom],
    rank: usize,
) -> Result<Vec<Atom>, OrbitProjectorError> {
    resources::checkpoint()?;
    let variables = (0..gram.len())
        .map(|index| {
            resources::checkpoint()?;
            Ok(Atom::var(symbol!(&format!(
                "FeynKit::OrbitProjectorRank{rank}_{index}"
            ))))
        })
        .collect::<Result<Vec<_>, OrbitProjectorError>>()?;
    let equations = (0..gram.len())
        .map(|column| {
            resources::checkpoint()?;
            variables.iter().enumerate().try_fold(
                -right_hand_side[column].clone(),
                |equation, (row, variable)| {
                    resources::checkpoint()?;
                    Ok(equation + &gram[row][column] * variable)
                },
            )
        })
        .collect::<Result<Vec<_>, OrbitProjectorError>>()?;
    resources::checkpoint()?;
    let (matrix, rhs) = Atom::system_to_matrix::<u16, _, _>(&equations, &variables)
        .map_err(|source| OrbitProjectorError::Solve { rank, source })?;
    resources::checkpoint()?;
    let solution = matrix
        .solve(&rhs)
        .map_err(|error| OrbitProjectorError::Solve {
            rank,
            source: SolveError::Other(error.to_string()),
        })?;
    resources::checkpoint()?;
    solution
        .into_vec()
        .into_iter()
        .map(|value| {
            resources::checkpoint()?;
            let expression = value.to_expression();
            resources::checkpoint()?;
            Ok(expression)
        })
        .collect()
}

fn dense_colors(colors: &[usize]) -> (Vec<usize>, usize) {
    let labels = colors.iter().copied().collect::<BTreeSet<_>>();
    let positions = labels
        .into_iter()
        .enumerate()
        .map(|(dense, label)| (label, dense))
        .collect::<BTreeMap<_, _>>();
    (
        colors.iter().map(|color| positions[color]).collect(),
        positions.len(),
    )
}

fn orbit_positions(alphas: &[Vec<u16>]) -> BTreeMap<Vec<u16>, usize> {
    alphas
        .iter()
        .cloned()
        .enumerate()
        .map(|(position, alpha)| (alpha, position))
        .collect()
}

fn upper_triangle_pairs(size: usize) -> Result<Vec<(usize, usize)>, OrbitProjectorError> {
    let count = triangle_len(size)?;
    let bytes = count
        .checked_mul(mem::size_of::<(usize, usize)>())
        .and_then(|bytes| u64::try_from(bytes).ok())
        .ok_or(OrbitProjectorError::CounterOverflow)?;
    resources::check_allocation(bytes)?;
    let mut pairs = Vec::with_capacity(count);
    for left in 0..size {
        for right in left..size {
            resources::checkpoint()?;
            pairs.push((left, right));
        }
    }
    Ok(pairs)
}

fn edge_color_count(edge_types: &[(usize, usize)]) -> usize {
    edge_types.last().map_or(0, |(_, right)| right + 1)
}

fn triangle_len(size: usize) -> Result<usize, OrbitProjectorError> {
    size.checked_add(1)
        .and_then(|next| size.checked_mul(next))
        .map(|product| product / 2)
        .ok_or(OrbitProjectorError::CounterOverflow)
}

fn upper_triangle_index(left: usize, right: usize, size: usize) -> usize {
    let (left, right) = if left <= right {
        (left, right)
    } else {
        (right, left)
    };
    left * size - left * left.saturating_sub(1) / 2 + right - left
}

fn split_joint(joint: usize, projector_class_count: usize) -> (usize, usize) {
    (joint / projector_class_count, joint % projector_class_count)
}

fn add_alpha(target: &mut [u16], contribution: &[u16]) -> Result<(), OrbitProjectorError> {
    for (target, contribution) in target.iter_mut().zip(contribution) {
        resources::checkpoint()?;
        *target = target
            .checked_add(*contribution)
            .ok_or(OrbitProjectorError::CounterOverflow)?;
    }
    Ok(())
}

fn add_count<K: Ord>(
    output: &mut BTreeMap<K, u128>,
    key: K,
    multiplicity: u128,
    limit: Option<usize>,
) -> Result<(), OrbitProjectorError> {
    resources::checkpoint()?;
    let is_new = !output.contains_key(&key);
    if let Some(limit) = limit
        && is_new
        && output.len() >= limit
    {
        return Err(OrbitProjectorError::OutputLimit {
            terms: output.len().saturating_add(1),
            limit,
        });
    }
    let count = output.entry(key).or_default();
    *count = count
        .checked_add(multiplicity)
        .ok_or(OrbitProjectorError::CounterOverflow)?;
    Ok(())
}

fn atom_from_count(count: u128) -> Result<Atom, OrbitProjectorError> {
    resources::checkpoint()?;
    let coefficient = Atom::num(count);
    resources::checkpoint()?;
    Ok(coefficient)
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use crate::{CosetType, OrthogonalWeingarten};
    use symbolica::parse;

    #[test]
    fn joint_color_dp_counts_every_labeled_matching_once() {
        for rank in (0..=20).step_by(2) {
            let internal = (0..rank)
                .map(|position| usize::from(position >= 2))
                .collect::<Vec<_>>();
            let projector = (0..rank)
                .map(|position| usize::from(position >= rank / 2))
                .collect::<Vec<_>>();
            let (internal, internal_classes) = dense_colors(&internal);
            let (projector, projector_classes) = dense_colors(&projector);
            let intersections = MatchingIntersections::new(
                &internal,
                internal_classes,
                &projector,
                projector_classes,
                Some(10_000),
            )
            .unwrap();
            assert_eq!(
                intersections.counts.values().sum::<u128>(),
                perfect_matching_count(rank)
            );
        }
    }

    #[test]
    fn alternating_cycle_dp_counts_every_relative_matching_once() {
        for multiplicities in [vec![4], vec![2, 4], vec![2, 2, 4], vec![2, 18]] {
            let colors = multiplicities
                .iter()
                .enumerate()
                .flat_map(|(color, count)| std::iter::repeat_n(color, *count))
                .collect::<Vec<_>>();
            let intersections = MatchingIntersections::new(
                &colors,
                multiplicities.len(),
                &vec![0; colors.len()],
                1,
                Some(10_000),
            )
            .unwrap();
            let alphas = intersections
                .counts
                .keys()
                .map(|key| key.internal.clone())
                .collect::<BTreeSet<_>>();
            for alpha in alphas {
                let distribution =
                    relative_matching_distribution(&alpha, multiplicities.len(), Some(10_000))
                        .unwrap();
                assert_eq!(
                    distribution.values().sum::<u128>(),
                    perfect_matching_count(colors.len()),
                    "multiplicities={multiplicities:?}, representative={alpha:?}"
                );
            }
        }
    }

    #[test]
    fn orbit_gram_at_dimension_one_has_raw_orbit_sizes_in_every_row() {
        let colors = [0, 0, 1, 1, 1, 1];
        let intersections = MatchingIntersections::new(&colors, 2, &[0; 6], 1, Some(100)).unwrap();
        let alphas = intersections
            .counts
            .keys()
            .map(|key| key.internal.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let sizes = alphas
            .iter()
            .map(|alpha| {
                intersections
                    .counts
                    .iter()
                    .filter(|(key, _)| &key.internal == alpha)
                    .map(|(_, count)| count)
                    .sum::<u128>()
            })
            .collect::<Vec<_>>();
        let gram = orbit_gram(&alphas, 2, &Atom::one(), Some(100)).unwrap();
        for row in gram {
            assert_eq!(
                row,
                sizes
                    .iter()
                    .map(|size| Atom::num(*size))
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn dense_orbit_product_obeys_the_output_budget_before_solving() {
        let error =
            OrbitProjector::new(&[0, 0, 1, 1], &[0, 1, 0, 1], Atom::one(), Some(3)).unwrap_err();
        assert!(matches!(
            error,
            OrbitProjectorError::OutputLimit { terms: 4, limit: 3 }
        ));
    }

    #[test]
    fn orbit_representation_overflows_are_reported() {
        assert!(matches!(
            MatchingIntersections::new(&[0; 256], 1, &[0; 256], 1, None),
            Err(OrbitProjectorError::CounterOverflow)
        ));
        assert!(matches!(
            relative_matching_distribution(&[256], 1, None),
            Err(OrbitProjectorError::CounterOverflow)
        ));
        assert!(matches!(
            add_alpha(&mut [u16::MAX], &[1]),
            Err(OrbitProjectorError::CounterOverflow)
        ));
    }

    #[test]
    fn reduced_projector_matches_explicit_weingarten_sum_through_rank_eight() {
        let cases = [
            (vec![0, 0, 1, 1], vec![0, 1, 0, 1]),
            (vec![0, 0, 1, 1, 1, 1], vec![0, 0, 0, 0, 1, 1]),
            (vec![0, 0, 1, 1, 1, 1, 1, 1], vec![0, 0, 0, 0, 1, 1, 1, 1]),
            // This reversal exercises `N H_B^-1`: the projector side has
            // fewer contraction orbits than the internal side.
            (vec![0, 0, 0, 0, 1, 1, 1, 1], vec![0, 0, 1, 1, 1, 1, 1, 1]),
        ];
        for (internal, projector) in cases {
            compare_with_brute_force(&internal, &projector);
        }
    }

    #[test]
    #[ignore = "rank-ten explicit Weingarten cross-check exercises 893,025 matching pairs"]
    fn rank_ten_matches_explicit_weingarten_sum() {
        compare_with_brute_force(
            &[0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
            &[0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        );
    }

    #[test]
    fn rank_twenty_partial_symmetries_are_practical() {
        let started = Instant::now();
        let mut internal = vec![0; 2];
        internal.extend([1; 18]);
        let mut projector = vec![1; 18];
        projector.extend([0; 2]);
        let result = OrbitProjector::new(&internal, &projector, parse!("D"), Some(100)).unwrap();
        let elapsed = started.elapsed();

        assert_eq!(result.terms().len(), 4);
        assert_eq!(
            result
                .terms()
                .iter()
                .map(|term| term.internal_alpha())
                .collect::<BTreeSet<_>>()
                .len(),
            2
        );
        assert_eq!(
            result
                .terms()
                .iter()
                .map(|term| term.projector_alpha())
                .collect::<BTreeSet<_>>()
                .len(),
            2
        );
        assert_eq!(
            result
                .terms()
                .iter()
                .map(|term| term.internal_labeled_pairings())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([34_459_425, 620_269_650])
        );
        assert!(
            result
                .terms()
                .iter()
                .all(|term| term.coefficient().to_string().contains('D'))
        );
        eprintln!("rank-20 [2,18] x [2,18] orbit projector: {elapsed:?}");
    }

    fn compare_with_brute_force(internal: &[usize], projector: &[usize]) {
        let dimension = parse!("D");
        let reduced =
            OrbitProjector::new(internal, projector, dimension.clone(), Some(10_000)).unwrap();
        let reduced = reduced
            .terms()
            .iter()
            .map(|term| {
                (
                    (
                        term.internal_alpha().to_vec(),
                        term.projector_alpha().to_vec(),
                    ),
                    term.coefficient().clone(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let pairings = perfect_matchings(internal.len());
        let table = OrthogonalWeingarten::for_rank(internal.len(), dimension).unwrap();
        let mut brute = BTreeMap::<(Vec<u16>, Vec<u16>), Atom>::new();
        let mut pair_total = 0_u128;
        for left in &pairings {
            let internal_alpha = pairing_alpha(internal, left);
            for right in &pairings {
                let projector_alpha = pairing_alpha(projector, right);
                let coset = CosetType::between_matchings(left, right).unwrap();
                *brute
                    .entry((internal_alpha.clone(), projector_alpha))
                    .or_insert(Atom::Zero) += table.coefficient(&coset).unwrap();
                pair_total += 1;
            }
        }
        assert_eq!(
            pair_total,
            perfect_matching_count(internal.len()).pow(2),
            "the explicit double-orbit multiplicities must cover every matching pair"
        );
        let keys = reduced
            .keys()
            .chain(brute.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        for key in keys {
            let actual = reduced.get(&key).cloned().unwrap_or(Atom::Zero);
            let expected = brute.get(&key).cloned().unwrap_or(Atom::Zero);
            assert!(
                (actual.clone() - &expected).together().is_zero(),
                "rank {} orbit {key:?}: reduced={actual}, brute={expected}",
                internal.len()
            );
        }
    }

    fn pairing_alpha(colors: &[usize], pairing: &[usize]) -> Vec<u16> {
        let class_count = colors
            .iter()
            .copied()
            .max()
            .map_or(0, |maximum| maximum + 1);
        let mut alpha = vec![0; triangle_len(class_count).unwrap()];
        for left in 0..pairing.len() {
            let right = pairing[left];
            if left < right {
                alpha[upper_triangle_index(colors[left], colors[right], class_count)] += 1;
            }
        }
        alpha
    }

    fn perfect_matching_count(rank: usize) -> u128 {
        (1..rank).step_by(2).map(|factor| factor as u128).product()
    }

    fn perfect_matchings(rank: usize) -> Vec<Vec<usize>> {
        fn generate(mask: u32, partners: &mut [usize], output: &mut Vec<Vec<usize>>) {
            if mask == 0 {
                output.push(partners.to_vec());
                return;
            }
            let left = mask.trailing_zeros() as usize;
            let without_left = mask & !(1 << left);
            let mut candidates = without_left;
            while candidates != 0 {
                let right = candidates.trailing_zeros() as usize;
                partners[left] = right;
                partners[right] = left;
                generate(without_left & !(1 << right), partners, output);
                candidates &= !(1 << right);
            }
        }

        let mut output = Vec::with_capacity(perfect_matching_count(rank) as usize);
        generate((1_u32 << rank) - 1, &mut vec![0; rank], &mut output);
        output
    }
}
