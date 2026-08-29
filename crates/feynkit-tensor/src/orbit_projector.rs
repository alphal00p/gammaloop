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

use symbolica::atom::{Atom, AtomCore};
use symbolica::symbol;
use thiserror::Error;

use crate::OrthogonalWeingarten;

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
    /// are canonicalized independently on both sides.  `output_limit` bounds
    /// intermediate and final orbit maps, not labeled matchings.
    pub(crate) fn new(
        internal_colors: &[usize],
        projector_colors: &[usize],
        dimension: Atom,
        output_limit: usize,
    ) -> Result<Self, OrbitProjectorError> {
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
        if rank > OrthogonalWeingarten::MAX_RANK {
            return Err(OrbitProjectorError::UnsupportedRank {
                rank,
                maximum: OrthogonalWeingarten::MAX_RANK,
            });
        }

        let (internal_colors, internal_class_count) = dense_colors(internal_colors);
        let (projector_colors, projector_class_count) = dense_colors(projector_colors);
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
        if term_bound > output_limit {
            return Err(OrbitProjectorError::OutputLimit {
                terms: term_bound,
                limit: output_limit,
            });
        }
        let internal_positions = orbit_positions(&internal_alphas);
        let projector_positions = orbit_positions(&projector_alphas);
        let mut incidence = vec![vec![0_u128; projector_alphas.len()]; internal_alphas.len()];
        for (key, count) in intersections.counts {
            incidence[internal_positions[&key.internal]][projector_positions[&key.projector]] =
                count;
        }
        let internal_sizes = incidence
            .iter()
            .map(|row| row.iter().sum::<u128>())
            .collect::<Vec<_>>();
        let projector_sizes = (0..projector_alphas.len())
            .map(|column| incidence.iter().map(|row| row[column]).sum::<u128>())
            .collect::<Vec<_>>();

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
                let right_hand_side = incidence
                    .iter()
                    .map(|row| atom_from_count(row[projector]))
                    .collect::<Result<Vec<_>, _>>()?;
                let column = solve_transposed(&gram, &right_hand_side, rank)?;
                for (internal, coefficient) in column.into_iter().enumerate() {
                    coefficients[internal][projector] = coefficient.cancel();
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
                let right_hand_side = incidence[internal]
                    .iter()
                    .copied()
                    .map(atom_from_count)
                    .collect::<Result<Vec<_>, _>>()?;
                let row = solve_transposed(&gram, &right_hand_side, rank)?;
                for (projector, coefficient) in row.into_iter().enumerate() {
                    coefficients[internal][projector] = coefficient.cancel();
                }
            }
        }

        let mut terms = Vec::with_capacity(internal_alphas.len() * projector_alphas.len());
        for (internal, internal_alpha) in internal_alphas.iter().enumerate() {
            for (projector, projector_alpha) in projector_alphas.iter().enumerate() {
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

        Ok(Self { terms })
    }

    pub(crate) fn terms(&self) -> &[OrbitProjectorTerm] {
        &self.terms
    }
}

/// Failures detected before a symmetry-reduced projector can be constructed.
#[derive(Debug, Error)]
pub(crate) enum OrbitProjectorError {
    #[error("internal colors have rank {internal}, but projector colors have rank {projector}")]
    ColorLengthMismatch { internal: usize, projector: usize },
    #[error("tensor rank {0} is odd; vacuum tensor reduction requires an even rank")]
    OddRank(usize),
    #[error("tensor rank {rank} exceeds the supported maximum {maximum}")]
    UnsupportedRank { rank: usize, maximum: usize },
    #[error("symmetry reduction produced {terms} orbit terms, above limit {limit}")]
    OutputLimit { terms: usize, limit: usize },
    #[error("pairing multiplicity {0} does not fit a signed 64-bit coefficient")]
    MultiplicityOverflow(u128),
    #[error("pairing multiplicity overflowed the rank-20 exact counter")]
    CounterOverflow,
    #[error("could not solve the rank-{rank} symmetry-reduced projector: {source}")]
    Solve {
        rank: usize,
        #[source]
        source: symbolica::solve::SolveError,
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
        limit: usize,
    ) -> Result<Self, OrbitProjectorError> {
        let mut joint_counts = vec![0_u8; internal_class_count * projector_class_count];
        for (&internal, &projector) in internal.iter().zip(projector) {
            joint_counts[internal * projector_class_count + projector] += 1;
        }
        let alpha_lengths = (
            triangle_len(internal_class_count),
            triangle_len(projector_class_count),
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
    limit: usize,
    cache: &mut HashMap<Vec<u8>, BTreeMap<OrbitPairKey, u128>>,
) -> Result<BTreeMap<OrbitPairKey, u128>, OrbitProjectorError> {
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
            key.internal
                [upper_triangle_index(anchor_internal, partner_internal, internal_class_count)] +=
                1;
            key.projector[upper_triangle_index(
                anchor_projector,
                partner_projector,
                projector_class_count,
            )] += 1;
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
    limit: usize,
) -> Result<Vec<Vec<Atom>>, OrbitProjectorError> {
    let positions = orbit_positions(alphas);
    let mut gram = Vec::with_capacity(alphas.len());
    for representative in alphas {
        let distribution = relative_matching_distribution(representative, class_count, limit)?;
        let mut row = vec![Atom::Zero; alphas.len()];
        for (key, multiplicity) in distribution {
            let target = positions[&key.alpha];
            row[target] += atom_from_count(multiplicity)?
                * dimension.clone().pow(Atom::num(i64::from(key.cycles)));
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
    limit: usize,
) -> Result<BTreeMap<RelativeKey, u128>, OrbitProjectorError> {
    let edge_types = upper_triangle_pairs(class_count);
    let counts = representative
        .iter()
        .copied()
        .map(|count| u8::try_from(count).expect("rank twenty fits in u8"))
        .collect::<Vec<_>>();
    let mut cache = HashMap::new();
    relative_components(
        &counts,
        &edge_types,
        triangle_len(class_count),
        limit,
        &mut cache,
    )
}

fn relative_components(
    counts: &[u8],
    edge_types: &[(usize, usize)],
    alpha_len: usize,
    limit: usize,
    cache: &mut HashMap<Vec<u8>, BTreeMap<RelativeKey, u128>>,
) -> Result<BTreeMap<RelativeKey, u128>, OrbitProjectorError> {
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
        for (mut rest, rest_multiplicity) in
            relative_components(&component.remaining, edge_types, alpha_len, limit, cache)?
        {
            add_alpha(&mut rest.alpha, &component.alpha);
            rest.cycles += 1;
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
    limit: usize,
    output: &mut BTreeMap<ComponentKey, u128>,
) -> Result<(), OrbitProjectorError> {
    let class_count = edge_color_count(edge_types);
    for next_type in 0..remaining.len() {
        let choices = remaining[next_type];
        if choices == 0 {
            continue;
        }
        let mut after_next = remaining.to_vec();
        after_next[next_type] -= 1;
        let (left, right) = edge_types[next_type];
        for (entrance, exit) in [(left, right), (right, left)] {
            let multiplicity = path_multiplicity
                .checked_mul(u128::from(choices))
                .ok_or(OrbitProjectorError::CounterOverflow)?;
            let mut extended_alpha = alpha.clone();
            extended_alpha[upper_triangle_index(last_exit, entrance, class_count)] += 1;

            let mut closed_alpha = extended_alpha.clone();
            closed_alpha[upper_triangle_index(exit, anchor_entrance, class_count)] += 1;
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
    limit: usize,
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
    let variables = (0..gram.len())
        .map(|index| {
            Atom::var(symbol!(&format!(
                "FeynKit::OrbitProjectorRank{rank}_{index}"
            )))
        })
        .collect::<Vec<_>>();
    let equations = (0..gram.len())
        .map(|column| {
            variables.iter().enumerate().fold(
                -right_hand_side[column].clone(),
                |equation, (row, variable)| equation + &gram[row][column] * variable,
            )
        })
        .collect::<Vec<_>>();
    Atom::solve_linear_system::<u16, _, _>(&equations, &variables)
        .map_err(|source| OrbitProjectorError::Solve { rank, source })
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

fn upper_triangle_pairs(size: usize) -> Vec<(usize, usize)> {
    (0..size)
        .flat_map(|left| (left..size).map(move |right| (left, right)))
        .collect()
}

fn edge_color_count(edge_types: &[(usize, usize)]) -> usize {
    edge_types.last().map_or(0, |(_, right)| right + 1)
}

fn triangle_len(size: usize) -> usize {
    size * (size + 1) / 2
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

fn add_alpha(target: &mut [u16], contribution: &[u16]) {
    for (target, contribution) in target.iter_mut().zip(contribution) {
        *target += contribution;
    }
}

fn add_count<K: Ord>(
    output: &mut BTreeMap<K, u128>,
    key: K,
    multiplicity: u128,
    limit: usize,
) -> Result<(), OrbitProjectorError> {
    let is_new = !output.contains_key(&key);
    if is_new && output.len() == limit {
        return Err(OrbitProjectorError::OutputLimit {
            terms: output.len() + 1,
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
    let count =
        i64::try_from(count).map_err(|_| OrbitProjectorError::MultiplicityOverflow(count))?;
    Ok(Atom::num(count))
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
                10_000,
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
                10_000,
            )
            .unwrap();
            let alphas = intersections
                .counts
                .keys()
                .map(|key| key.internal.clone())
                .collect::<BTreeSet<_>>();
            for alpha in alphas {
                let distribution =
                    relative_matching_distribution(&alpha, multiplicities.len(), 10_000).unwrap();
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
        let intersections = MatchingIntersections::new(&colors, 2, &[0; 6], 1, 100).unwrap();
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
        let gram = orbit_gram(&alphas, 2, &Atom::one(), 100).unwrap();
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
        let error = OrbitProjector::new(&[0, 0, 1, 1], &[0, 1, 0, 1], Atom::one(), 3).unwrap_err();
        assert!(matches!(
            error,
            OrbitProjectorError::OutputLimit { terms: 4, limit: 3 }
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
        let result = OrbitProjector::new(&internal, &projector, parse!("D"), 100).unwrap();
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
        let reduced = OrbitProjector::new(internal, projector, dimension.clone(), 10_000).unwrap();
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
        let mut alpha = vec![0; triangle_len(class_count)];
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
