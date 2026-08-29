//! Orthogonal-Weingarten coefficients for Lorentz-invariant projectors.
//!
//! This implements the Collins--Matsumoto orthogonality recurrence.  At rank
//! `2k`, invariance under the hyperoctahedral group reduces the coefficient
//! system from `(2k - 1)!!` pairings to `p(k)` integer partitions.
//! This is the coefficient-generation counterpart of the orbit projectors in
//! Goode *et al.*, [arXiv:2408.05137](https://arxiv.org/abs/2408.05137); the
//! recurrence is Lemma 4.3 of Collins and Matsumoto,
//! [arXiv:1701.04493](https://arxiv.org/abs/1701.04493).

use std::{
    collections::BTreeMap,
    sync::{Mutex, OnceLock},
};

use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};
use thiserror::Error;

/// The coset type of two perfect matchings.
///
/// Superimposing two matchings gives disjoint alternating cycles.  Dividing
/// every cycle length by two produces an integer partition of the number of
/// pairs; that partition is the coset type.  Orthogonal-Weingarten
/// coefficients are constant on each such class.
///
/// # Examples
///
/// The crossed and uncrossed contractions of a rank-four vacuum tensor have
/// coset types `[2]` and `[1, 1]`, respectively.
///
/// ```
/// use feynkit_tensor::CosetType;
///
/// let crossed = CosetType::new([2])?;
/// assert_eq!(crossed.parts(), [2]);
/// # Ok::<(), feynkit_tensor::WeingartenError>(())
/// ```
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct CosetType {
    parts: Box<[usize]>,
}

impl CosetType {
    /// Construct a coset type from positive partition parts.
    ///
    /// The parts are canonicalized into non-increasing order.  The empty
    /// partition is the unique rank-zero coset type.  Coset types above the
    /// validated rank-20 target are rejected.
    ///
    /// # Examples
    ///
    /// ```
    /// use feynkit_tensor::CosetType;
    ///
    /// let orbit = CosetType::new([1, 3, 1])?;
    /// assert_eq!(orbit.parts(), [3, 1, 1]);
    /// # Ok::<(), feynkit_tensor::WeingartenError>(())
    /// ```
    pub fn new(parts: impl IntoIterator<Item = usize>) -> Result<Self, WeingartenError> {
        let mut parts = parts.into_iter().collect::<Vec<_>>();
        if parts.contains(&0) {
            return Err(WeingartenError::ZeroPartitionPart);
        }
        let pair_count = parts.iter().try_fold(0_usize, |total, part| {
            total
                .checked_add(*part)
                .ok_or(WeingartenError::RankOverflow)
        })?;
        if pair_count > OrthogonalWeingarten::MAX_RANK / 2 {
            return Err(WeingartenError::UnsupportedRank {
                rank: pair_count.saturating_mul(2),
                maximum: OrthogonalWeingarten::MAX_RANK,
            });
        }
        parts.sort_unstable_by(|left, right| right.cmp(left));
        Ok(Self {
            parts: parts.into_boxed_slice(),
        })
    }

    /// Infer the coset type of two matchings represented by partner arrays.
    ///
    /// A matching `[1, 0, 3, 2]`, for example, contains the pairs `(0, 1)`
    /// and `(2, 3)`.  Both arrays must be fixed-point-free involutions of the
    /// same even length.
    ///
    /// # Examples
    ///
    /// ```
    /// use feynkit_tensor::CosetType;
    ///
    /// let direct = [1, 0, 3, 2];
    /// let crossed = [2, 3, 0, 1];
    /// assert_eq!(CosetType::between_matchings(&direct, &crossed)?.parts(), [2]);
    /// # Ok::<(), feynkit_tensor::WeingartenError>(())
    /// ```
    pub fn between_matchings(left: &[usize], right: &[usize]) -> Result<Self, WeingartenError> {
        validate_matching(left)?;
        validate_matching(right)?;
        if left.len() != right.len() {
            return Err(WeingartenError::MatchingLengthMismatch {
                left: left.len(),
                right: right.len(),
            });
        }

        let mut seen = vec![false; left.len()];
        let mut parts = Vec::new();
        for start in 0..left.len() {
            if seen[start] {
                continue;
            }
            let mut stack = vec![start];
            let mut vertices = 0;
            while let Some(vertex) = stack.pop() {
                if seen[vertex] {
                    continue;
                }
                seen[vertex] = true;
                vertices += 1;
                stack.push(left[vertex]);
                stack.push(right[vertex]);
            }
            debug_assert_eq!(vertices % 2, 0);
            parts.push(vertices / 2);
        }
        Self::new(parts)
    }

    /// Integer-partition parts in canonical non-increasing order.
    pub fn parts(&self) -> &[usize] {
        &self.parts
    }

    /// Number of pairs, so the associated tensor rank is twice this value.
    pub fn pair_count(&self) -> usize {
        self.parts.iter().sum()
    }

    /// Number of matchings in this coset class relative to one fixed matching.
    /// This is the orbit-size formula
    /// `2^(k - length(lambda)) k! / z_lambda` of Eq. (2.39) in Goode *et al.*
    /// for a partition `lambda` of `k`.
    ///
    /// # Examples
    ///
    /// There are six rank-six pairings in the `[2, 1]` orbit.
    ///
    /// ```
    /// use feynkit_tensor::CosetType;
    ///
    /// assert_eq!(CosetType::new([2, 1])?.class_size(), 6);
    /// # Ok::<(), feynkit_tensor::WeingartenError>(())
    /// ```
    pub fn class_size(&self) -> u128 {
        let pair_count = self.pair_count();
        let mut numerator = 1_u128 << (pair_count - self.parts.len());
        numerator *= factorial(pair_count);

        let mut denominator = 1_u128;
        let mut offset = 0;
        while offset < self.parts.len() {
            let part = self.parts[offset];
            let multiplicity = self.parts[offset..]
                .iter()
                .take_while(|candidate| **candidate == part)
                .count();
            denominator *= (part as u128).pow(multiplicity as u32);
            denominator *= factorial(multiplicity);
            offset += multiplicity;
        }
        numerator / denominator
    }

    fn without_one(&self) -> Option<Self> {
        let one = self.parts.iter().rposition(|part| *part == 1)?;
        let mut parts = self.parts.to_vec();
        parts.remove(one);
        Some(Self {
            parts: parts.into_boxed_slice(),
        })
    }
}

/// Orthogonal-Weingarten coefficients at one even tensor rank.
///
/// Coefficients are exact Symbolica expressions in the caller-supplied
/// dimension.  Generation proceeds through systems of size `p(k)` at rank
/// `2k`; in particular, rank 20 has only 42 symmetry classes.
///
/// # Examples
///
/// Generate the Lorentz projector coefficients for a rank-four vacuum
/// numerator in dimensional regularization.
///
/// ```no_run
/// use feynkit_tensor::{CosetType, OrthogonalWeingarten};
/// use symbolica::parse;
///
/// let table = OrthogonalWeingarten::for_rank(4, parse!("4 - 2*eps"))?;
/// let crossed = table.coefficient(&CosetType::new([2])?)?;
/// println!("crossed projector coefficient = {crossed}");
/// # Ok::<(), feynkit_tensor::WeingartenError>(())
/// ```
#[derive(Clone, Debug)]
pub struct OrthogonalWeingarten {
    pair_count: usize,
    dimension: Atom,
    coefficients: BTreeMap<CosetType, Atom>,
}

impl OrthogonalWeingarten {
    /// Largest tensor rank supported by the native coefficient engine.
    pub const MAX_RANK: usize = 20;

    /// Generate all coefficients for `pair_count` metric pairs.
    ///
    /// `dimension` may itself be symbolic, such as `4 - 2*eps`.  A rank-zero
    /// table contains the empty coset type with coefficient one.  At most ten
    /// pairs (rank 20) are accepted. Keep mixed high-rank projectors symbolic:
    /// at fixed positive integer dimension `d < pair_count`, the labeled
    /// metric basis is redundant and this constructor returns
    /// [`WeingartenError::SingularDimension`]. The all-equal isotropic formula
    /// does not require this inverse.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feynkit_tensor::OrthogonalWeingarten;
    /// use symbolica::parse;
    ///
    /// let rank_six = OrthogonalWeingarten::new(3, parse!("D"))?;
    /// assert_eq!(rank_six.rank(), 6);
    /// # Ok::<(), feynkit_tensor::WeingartenError>(())
    /// ```
    pub fn new(pair_count: usize, dimension: Atom) -> Result<Self, WeingartenError> {
        if pair_count > Self::MAX_RANK / 2 {
            return Err(WeingartenError::UnsupportedRank {
                rank: pair_count.saturating_mul(2),
                maximum: Self::MAX_RANK,
            });
        }
        Self::validate_inverse_dimension(pair_count, &dimension)?;
        let canonical_dimension = canonical_dimension();
        let canonical_coefficients = canonical_coefficients(pair_count, &canonical_dimension)?;

        let coefficients = if dimension == canonical_dimension {
            canonical_coefficients
        } else {
            canonical_coefficients
                .into_iter()
                .map(|(coset_type, coefficient)| {
                    let coefficient = coefficient
                        .replace(canonical_dimension.to_pattern())
                        .with(dimension.clone());
                    (coset_type, coefficient)
                })
                .collect()
        };

        Ok(Self {
            pair_count,
            dimension,
            coefficients,
        })
    }

    /// Generate coefficients for an even tensor rank.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feynkit_tensor::OrthogonalWeingarten;
    /// use symbolica::parse;
    ///
    /// let projector = OrthogonalWeingarten::for_rank(12, parse!("4 - 2*eps"))?;
    /// assert_eq!(projector.coefficients().len(), 11);
    /// # Ok::<(), feynkit_tensor::WeingartenError>(())
    /// ```
    pub fn for_rank(rank: usize, dimension: Atom) -> Result<Self, WeingartenError> {
        if !rank.is_multiple_of(2) {
            return Err(WeingartenError::OddRank(rank));
        }
        Self::new(rank / 2, dimension)
    }

    /// Tensor rank covered by this table.
    pub fn rank(&self) -> usize {
        2 * self.pair_count
    }

    /// Number of metric pairs covered by this table.
    pub fn pair_count(&self) -> usize {
        self.pair_count
    }

    /// Dimension expression used in the generated coefficients.
    pub fn dimension(&self) -> &Atom {
        &self.dimension
    }

    /// Look up the coefficient of a coset type.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feynkit_tensor::{CosetType, OrthogonalWeingarten};
    /// use symbolica::parse;
    ///
    /// let projector = OrthogonalWeingarten::for_rank(6, parse!("D"))?;
    /// let coefficient = projector.coefficient(&CosetType::new([3])?)?;
    /// println!("fully crossed rank-six coefficient = {coefficient}");
    /// # Ok::<(), feynkit_tensor::WeingartenError>(())
    /// ```
    pub fn coefficient(&self, coset_type: &CosetType) -> Result<&Atom, WeingartenError> {
        if coset_type.pair_count() != self.pair_count {
            return Err(WeingartenError::CosetRankMismatch {
                table_rank: self.rank(),
                coset_rank: 2 * coset_type.pair_count(),
            });
        }
        Ok(&self.coefficients[coset_type])
    }

    /// Iterate over all coset types and exact coefficients.
    pub fn coefficients(&self) -> impl ExactSizeIterator<Item = (&CosetType, &Atom)> {
        self.coefficients.iter()
    }

    /// Coefficient of every metric pairing for a totally symmetric tensor.
    ///
    /// For an integrand containing one repeated loop momentum, this is the
    /// closed isotropic-tensor result
    /// `1 / (D (D + 2) ... (D + 2k - 2))`.  It avoids constructing the full
    /// Weingarten table and is therefore the preferred rank-20 fast path.
    ///
    /// # Examples
    ///
    /// A numerator with twenty powers of the same loop momentum uses this
    /// coefficient for every product of ten metrics.
    ///
    /// ```
    /// use feynkit_tensor::OrthogonalWeingarten;
    /// use symbolica::parse;
    ///
    /// let coefficient =
    ///     OrthogonalWeingarten::isotropic_pairing_coefficient(10, parse!("D"));
    /// assert!(coefficient.to_string().contains("D"));
    /// ```
    pub fn isotropic_pairing_coefficient(pair_count: usize, dimension: Atom) -> Atom {
        let denominator = (0..pair_count).fold(Atom::num(1), |denominator, offset| {
            denominator * (dimension.clone() + Atom::num(2 * offset as i64))
        });
        Atom::num(1) / denominator
    }

    pub(crate) fn validate_inverse_dimension(
        pair_count: usize,
        dimension: &Atom,
    ) -> Result<(), WeingartenError> {
        if pair_count == 0 {
            return Ok(());
        }
        if let Ok(fixed_dimension) = i64::try_from(dimension.as_view())
            && (fixed_dimension <= 0 || fixed_dimension < pair_count as i64)
        {
            return Err(WeingartenError::SingularDimension {
                rank: 2 * pair_count,
                dimension: dimension.clone(),
            });
        }
        Ok(())
    }
}

fn canonical_coefficients(
    pair_count: usize,
    dimension: &Atom,
) -> Result<BTreeMap<CosetType, Atom>, WeingartenError> {
    static TABLES: OnceLock<Mutex<Vec<BTreeMap<CosetType, Atom>>>> = OnceLock::new();
    let tables = TABLES.get_or_init(|| {
        Mutex::new(vec![BTreeMap::from([(
            CosetType {
                parts: Box::new([]),
            },
            Atom::num(1),
        )])])
    });
    let mut tables = tables
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    while tables.len() <= pair_count {
        let current_pairs = tables.len();
        let next = solve_level(current_pairs, dimension, tables.last().unwrap())?;
        tables.push(next);
    }
    Ok(tables[pair_count].clone())
}

/// Failures in coset-type validation or coefficient generation.
///
/// # Examples
///
/// An odd-rank vacuum numerator has no Lorentz-invariant metric projector.
///
/// ```
/// use feynkit_tensor::{OrthogonalWeingarten, WeingartenError};
/// use symbolica::parse;
///
/// let error = OrthogonalWeingarten::for_rank(3, parse!("D")).unwrap_err();
/// assert!(matches!(error, WeingartenError::OddRank(3)));
/// ```
#[derive(Debug, Error)]
pub enum WeingartenError {
    /// Tensor projectors only exist at even rank.
    #[error("tensor rank {0} is odd; vacuum tensor reduction requires an even rank")]
    OddRank(usize),
    /// The exact engine is bounded to the validated high-rank target.
    #[error("tensor rank {rank} exceeds the supported maximum {maximum}")]
    UnsupportedRank { rank: usize, maximum: usize },
    /// Summing user-provided partition parts overflowed the platform integer.
    #[error("coset-type tensor rank is too large to represent")]
    RankOverflow,
    /// Integer partitions cannot contain zero.
    #[error("coset-type partition parts must be positive")]
    ZeroPartitionPart,
    /// Matching arrays must describe the same vertex set.
    #[error("matching lengths differ: left has {left} vertices, right has {right}")]
    MatchingLengthMismatch { left: usize, right: usize },
    /// A partner array did not describe a perfect matching.
    #[error("invalid matching at vertex {vertex}: partner is {partner}")]
    InvalidMatching { vertex: usize, partner: usize },
    /// A coefficient was requested from a table at another rank.
    #[error("coset type has rank {coset_rank}, but this table has rank {table_rank}")]
    CosetRankMismatch {
        table_rank: usize,
        coset_rank: usize,
    },
    /// At low fixed integer dimension, metric pairings obey dimension-specific
    /// identities and the universal Gram matrix has no inverse.
    #[error(
        "rank-{rank} metric-pairing Gram matrix is singular in fixed dimension {dimension}; keep the dimension symbolic or use the all-equal isotropic projector"
    )]
    SingularDimension { rank: usize, dimension: Atom },
    /// Symbolica could not solve the reduced recurrence system.
    #[error("could not solve rank-{rank} orthogonal-Weingarten system: {source}")]
    Solve {
        rank: usize,
        #[source]
        source: symbolica::solve::SolveError,
    },
}

fn solve_level(
    pair_count: usize,
    dimension: &Atom,
    lower: &BTreeMap<CosetType, Atom>,
) -> Result<BTreeMap<CosetType, Atom>, WeingartenError> {
    let coset_types = partitions(pair_count);
    let variables = (0..coset_types.len())
        .map(|index| {
            Atom::var(symbol!(&format!(
                "FeynKit::OrthogonalWeingarten{pair_count}_{index}"
            )))
        })
        .collect::<Vec<_>>();
    let positions = coset_types
        .iter()
        .enumerate()
        .map(|(index, coset_type)| (coset_type.clone(), index))
        .collect::<BTreeMap<_, _>>();

    let equations = coset_types
        .iter()
        .enumerate()
        .map(|(row, coset_type)| {
            let matching = representative(coset_type);
            let mut equation = dimension * &variables[row];
            let pivot = 2 * pair_count - 2;
            for vertex in 0..pivot {
                let moved = relabel_transposition(&matching, vertex, pivot);
                let target = CosetType::between_matchings(&reference_matching(pair_count), &moved)
                    .expect("a relabelled perfect matching remains valid");
                equation += &variables[positions[&target]];
            }
            if let Some(lower_type) = coset_type.without_one() {
                equation -= &lower[&lower_type];
            }
            equation
        })
        .collect::<Vec<_>>();

    let solutions =
        Atom::solve_linear_system::<u16, _, _>(&equations, &variables).map_err(|source| {
            WeingartenError::Solve {
                rank: 2 * pair_count,
                source,
            }
        })?;

    Ok(coset_types
        .into_iter()
        .zip(solutions.into_iter().map(|solution| solution.cancel()))
        .collect())
}

fn canonical_dimension() -> Atom {
    Atom::var(symbol!("FeynKit::Dimension"))
}

fn partitions(total: usize) -> Vec<CosetType> {
    fn append(total: usize, maximum: usize, current: &mut Vec<usize>, output: &mut Vec<CosetType>) {
        if total == 0 {
            output.push(CosetType {
                parts: current.clone().into_boxed_slice(),
            });
            return;
        }
        for next in (1..=total.min(maximum)).rev() {
            current.push(next);
            append(total - next, next, current, output);
            current.pop();
        }
    }

    let mut output = Vec::new();
    append(total, total, &mut Vec::new(), &mut output);
    output
}

fn reference_matching(pair_count: usize) -> Vec<usize> {
    (0..2 * pair_count).map(|vertex| vertex ^ 1).collect()
}

fn representative(coset_type: &CosetType) -> Vec<usize> {
    let mut matching = vec![usize::MAX; 2 * coset_type.pair_count()];
    let mut start = 0;
    for &part in coset_type.parts() {
        for pair in 0..part {
            let right = start + 2 * pair + 1;
            let next_left = start + 2 * ((pair + 1) % part);
            matching[right] = next_left;
            matching[next_left] = right;
        }
        start += 2 * part;
    }
    matching
}

fn relabel_transposition(matching: &[usize], left: usize, right: usize) -> Vec<usize> {
    let transpose = |vertex| {
        if vertex == left {
            right
        } else if vertex == right {
            left
        } else {
            vertex
        }
    };
    let mut moved = vec![usize::MAX; matching.len()];
    for (vertex, partner) in matching.iter().copied().enumerate() {
        moved[transpose(vertex)] = transpose(partner);
    }
    moved
}

fn validate_matching(matching: &[usize]) -> Result<(), WeingartenError> {
    for (vertex, partner) in matching.iter().copied().enumerate() {
        if partner >= matching.len() || partner == vertex || matching[partner] != vertex {
            return Err(WeingartenError::InvalidMatching { vertex, partner });
        }
    }
    Ok(())
}

fn factorial(value: usize) -> u128 {
    (2..=value).map(|factor| factor as u128).product()
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use symbolica::{
        atom::{Atom, AtomCore},
        parse,
    };

    use super::{CosetType, OrthogonalWeingarten, partitions, reference_matching, representative};

    #[test]
    fn partition_system_sizes_reach_only_42_at_rank_twenty() {
        let counts = (0..=10)
            .map(|pair_count| partitions(pair_count).len())
            .collect::<Vec<_>>();
        assert_eq!(counts, [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42]);
        assert_eq!(
            partitions(10)
                .iter()
                .map(CosetType::class_size)
                .sum::<u128>(),
            654_729_075
        );
    }

    #[test]
    fn representatives_have_the_requested_coset_type() {
        for pair_count in 0..=10 {
            let reference = reference_matching(pair_count);
            for coset_type in partitions(pair_count) {
                assert_eq!(
                    CosetType::between_matchings(&reference, &representative(&coset_type)).unwrap(),
                    coset_type
                );
            }
        }
    }

    #[test]
    fn empty_coset_class_contains_the_empty_matching() {
        let empty = CosetType::new([]).unwrap();
        assert_eq!(empty.class_size(), 1);
        assert_eq!(CosetType::between_matchings(&[], &[]).unwrap(), empty);
    }

    #[test]
    fn malformed_odd_matching_is_rejected() {
        assert!(CosetType::between_matchings(&[1, 0, 2], &[1, 0, 2]).is_err());
    }

    #[test]
    fn ranks_above_twenty_are_rejected_before_combinatorics() {
        assert!(matches!(
            OrthogonalWeingarten::new(11, parse!("D")),
            Err(super::WeingartenError::UnsupportedRank {
                rank: 22,
                maximum: 20
            })
        ));
        assert!(matches!(
            CosetType::new([usize::MAX, 1]),
            Err(super::WeingartenError::RankOverflow)
        ));
        assert!(CosetType::new([11]).is_err());
    }

    #[test]
    fn rank_four_matches_the_closed_form() {
        let dimension = parse!("D");
        let table = OrthogonalWeingarten::new(2, dimension.clone()).unwrap();
        let denominator = dimension.clone()
            * (dimension.clone() - Atom::num(1))
            * (dimension.clone() + Atom::num(2));
        let diagonal = (dimension.clone() + Atom::num(1)) / denominator.clone();
        let crossed = -Atom::num(1) / denominator;

        assert!(
            (table.coefficient(&CosetType::new([1, 1]).unwrap()).unwrap() - diagonal)
                .together()
                .is_zero()
        );
        assert!(
            (table.coefficient(&CosetType::new([2]).unwrap()).unwrap() - crossed)
                .together()
                .is_zero()
        );
    }

    #[test]
    fn exceptional_fixed_dimensions_are_rejected_before_substitution() {
        assert!(matches!(
            OrthogonalWeingarten::for_rank(10, Atom::num(4)),
            Err(super::WeingartenError::SingularDimension {
                rank: 10,
                dimension
            }) if dimension == Atom::num(4)
        ));
        assert!(OrthogonalWeingarten::for_rank(8, Atom::num(4)).is_ok());

        // The totally symmetric projector has a closed isotropic form and
        // does not invert the redundant labeled-pairing basis.
        let isotropic = OrthogonalWeingarten::isotropic_pairing_coefficient(10, Atom::num(4));
        assert!(!isotropic.is_zero());
    }

    #[test]
    fn isotropic_fast_path_is_the_sum_over_coset_classes() {
        let dimension = parse!("D");
        for pair_count in 1..=5 {
            let table = OrthogonalWeingarten::new(pair_count, dimension.clone()).unwrap();
            let summed = table
                .coefficients()
                .fold(Atom::num(0), |sum, (kind, value)| {
                    sum + Atom::num(kind.class_size()) * value
                });
            let expected =
                OrthogonalWeingarten::isotropic_pairing_coefficient(pair_count, dimension.clone());
            assert!((summed - expected).together().is_zero());
        }
    }

    #[test]
    fn vakint_projector_tables_are_reproduced_through_rank_ten() {
        let epsilon = parse!("ep");
        let dimension = Atom::num(4) - Atom::num(2) * epsilon;
        for rank in [2, 4, 6, 8, 10] {
            let table = OrthogonalWeingarten::for_rank(rank, dimension.clone()).unwrap();
            let expected = vakint_coefficients(rank);
            let types = vakint_coset_order(rank / 2);
            assert_eq!(expected.len(), types.len());
            for (index, (expected, coset_type)) in expected.into_iter().zip(types).enumerate() {
                let actual = table.coefficient(&coset_type).unwrap();
                assert!(
                    (actual - &expected).together().is_zero(),
                    "VAKINT rank-{rank} coefficient {} ({:?}) differs: FeynKit={}, VAKINT={}",
                    index + 1,
                    coset_type.parts(),
                    actual,
                    expected
                );
            }
        }
    }

    #[test]
    fn rank_twenty_table_is_practical() {
        // Goode et al.'s `coefficients.prc` ancillary is expressed in their
        // transformed symmetric-`dd` orbit basis, rather than as raw
        // orthogonal-Weingarten coefficients.  The exact VAKINT projector
        // entries above are therefore the direct coefficient oracle here.
        let started = Instant::now();
        let table = OrthogonalWeingarten::new(10, parse!("D")).unwrap();
        let elapsed = started.elapsed();
        assert_eq!(table.coefficients().len(), 42);
        assert!(
            table
                .coefficients()
                .all(|(_, coefficient)| coefficient.cancel() == coefficient.clone())
        );
        let isotropic_sum = table
            .coefficients()
            .fold(Atom::num(0), |sum, (kind, coefficient)| {
                sum + Atom::num(kind.class_size()) * coefficient
            });
        let isotropic_closed = OrthogonalWeingarten::isotropic_pairing_coefficient(10, parse!("D"));
        assert!((isotropic_sum - isotropic_closed).together().is_zero());
        eprintln!("rank-20 orthogonal-Weingarten table: {elapsed:?}");
    }

    fn vakint_coset_order(pair_count: usize) -> Vec<CosetType> {
        let ascending_largest_part = match pair_count {
            1 => vec![vec![1]],
            2 => vec![vec![1, 1], vec![2]],
            3 => vec![vec![1, 1, 1], vec![2, 1], vec![3]],
            4 => vec![
                vec![1, 1, 1, 1],
                vec![2, 1, 1],
                vec![2, 2],
                vec![3, 1],
                vec![4],
            ],
            5 => vec![
                vec![1, 1, 1, 1, 1],
                vec![2, 1, 1, 1],
                vec![2, 2, 1],
                vec![3, 1, 1],
                vec![3, 2],
                vec![4, 1],
                vec![5],
            ],
            _ => unreachable!("VAKINT only preloads projector tables through rank ten"),
        };
        ascending_largest_part
            .into_iter()
            .map(|parts| CosetType::new(parts).unwrap())
            .collect()
    }

    fn vakint_coefficients(rank: usize) -> Vec<Atom> {
        // Exact MIT-licensed oracle transcribed from VAKINT's FORM tables
        // `tensorreduce.frm` and `pvtab10.h` at D = 4 - 2 ep.  Keeping the
        // compact coefficient data here lets this crate package independently
        // instead of reaching into a sibling workspace member at compile time.
        let expressions: &[&str] = match rank {
            2 => &["1/(4-2*ep)"],
            4 => &[
                "(5-2*ep)/(72-108*ep+52*ep^2-8*ep^3)",
                "-1/(72-108*ep+52*ep^2-8*ep^3)",
            ],
            6 => &[
                "(26-22*ep+4*ep^2)/(1152-3168*ep+3280*ep^2-1600*ep^3+368*ep^4-32*ep^5)",
                "-1/(192-464*ep+392*ep^2-136*ep^3+16*ep^4)",
                "2/(1152-3168*ep+3280*ep^2-1600*ep^3+368*ep^4-32*ep^5)",
            ],
            8 => &[
                "(287-278*ep+84*ep^2-8*ep^3)/(28800-125280*ep+199504*ep^2-159680*ep^3+71152*ep^4-17888*ep^5+2368*ep^6-128*ep^7)",
                "(-166+198*ep-72*ep^2+8*ep^3)/(57600-308160*ep+649568*ep^2-718368*ep^3+461664*ep^4-178080*ep^5+40512*ep^6-4992*ep^7+256*ep^8)",
                "(54-26*ep+4*ep^2)/(57600-308160*ep+649568*ep^2-718368*ep^3+461664*ep^4-178080*ep^5+40512*ep^6-4992*ep^7+256*ep^8)",
                "2/(1800-8280*ep+13864*ep^2-11016*ep^3+4432*ep^4-864*ep^5+64*ep^6)",
                "(-26+10*ep)/(57600-308160*ep+649568*ep^2-718368*ep^3+461664*ep^4-178080*ep^5+40512*ep^6-4992*ep^7+256*ep^8)",
            ],
            10 => &[
                "(1280-5628*ep+6324*ep^2-2728*ep^3+496*ep^4-32*ep^5)/(-1382400*ep+7626240*ep^2-16822272*ep^3+19839104*ep^4-13953408*ep^5+6120576*ep^6-1684608*ep^7+281856*ep^8-26112*ep^9+1024*ep^10)",
                "(-640+1412*ep-892*ep^2+208*ep^3-16*ep^4)/(-1382400*ep+7626240*ep^2-16822272*ep^3+19839104*ep^4-13953408*ep^5+6120576*ep^6-1684608*ep^7+281856*ep^8-26112*ep^9+1024*ep^10)",
                "(32-22*ep+4*ep^2)/(-138240*ep+734976*ep^2-1535232*ep^3+1676864*ep^4-1059968*ep^5+400064*ep^6-88448*ep^7+10496*ep^8-512*ep^9)",
                "(320-452*ep+160*ep^2-16*ep^3)/(-1382400*ep+7626240*ep^2-16822272*ep^3+19839104*ep^4-13953408*ep^5+6120576*ep^6-1684608*ep^7+281856*ep^8-26112*ep^9+1024*ep^10)",
                "(-160+60*ep-8*ep^2)/(-1382400*ep+7626240*ep^2-16822272*ep^3+19839104*ep^4-13953408*ep^5+6120576*ep^6-1684608*ep^7+281856*ep^8-26112*ep^9+1024*ep^10)",
                "(-16+10*ep)/(-138240*ep+734976*ep^2-1535232*ep^3+1676864*ep^4-1059968*ep^5+400064*ep^6-88448*ep^7+10496*ep^8-512*ep^9)",
                "(80-28*ep)/(-1382400*ep+7626240*ep^2-16822272*ep^3+19839104*ep^4-13953408*ep^5+6120576*ep^6-1684608*ep^7+281856*ep^8-26112*ep^9+1024*ep^10)",
            ],
            _ => unreachable!("VAKINT only preloads projector tables through rank ten"),
        };
        expressions
            .iter()
            .map(|expression| parse!(expression))
            .collect()
    }
}
