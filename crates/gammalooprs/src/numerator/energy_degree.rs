use std::{
    borrow::Borrow,
    collections::{BTreeMap, BTreeSet},
};

use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski},
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::domains::rational::Rational;
use thiserror::Error;
use three_dimensional_reps::utils::{rank_i64, solve_rational_system};

use crate::{
    graph::Graph,
    momentum::SignOrZero,
    utils::{
        GS,
        symbols::{UvDenominatorClassId, UvMomentumProvenanceRole},
    },
    uv::UltravioletGraph,
};

// Bound native CFF trials while sharing the same budget for hard and soft dispatch.
const ENERGY_ASSIGNMENT_PROPOSAL_BUDGET: usize = 3;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct EnergyPowerCapMap<K = EdgeIndex> {
    degrees: BTreeMap<K, usize>,
}

impl<K> Default for EnergyPowerCapMap<K> {
    fn default() -> Self {
        Self {
            degrees: BTreeMap::new(),
        }
    }
}

impl<K: Copy + Ord> EnergyPowerCapMap<K> {
    pub fn is_empty(&self) -> bool {
        self.degrees.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (K, usize)> + '_ {
        self.degrees.iter().map(|(edge, degree)| (*edge, *degree))
    }

    fn unit(edge: K) -> Self {
        let mut degrees = BTreeMap::new();
        degrees.insert(edge, 1);
        Self { degrees }
    }

    fn add_assign(&mut self, other: Self) -> Result<(), EnergyPowerAnalysisError> {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = slot
                .checked_add(degree)
                .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
        }
        Ok(())
    }

    pub(crate) fn max_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = (*slot).max(degree);
        }
    }

    fn scale(&mut self, exponent: usize) -> Result<(), EnergyPowerAnalysisError> {
        for degree in self.degrees.values_mut() {
            *degree = degree
                .checked_mul(exponent)
                .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
        }
        Ok(())
    }
}

impl EnergyPowerCapMap {
    pub fn into_generation_bounds(self) -> Vec<(usize, usize)> {
        self.degrees
            .into_iter()
            .map(|(edge, degree)| (usize::from(edge), degree))
            .collect()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum EnergyReference {
    Physical(EdgeIndex),
    UvClass(UvDenominatorClassId),
}

impl EnergyReference {
    fn missing_candidates(self, degree: usize) -> EnergyPowerAnalysisError {
        match self {
            Self::Physical(edge) => EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates {
                edge: edge.into(),
                degree,
            },
            Self::UvClass(class) => EnergyPowerAnalysisError::InvalidClassEnergyCandidates {
                class: class.0,
                detail: format!("no candidates for degree {degree}"),
            },
        }
    }
}

/// Exact-energy occurrences certified for one physical owner or UV class.
///
/// Untagged and fixed numerator factors stay on the explicitly designated base
/// occurrence. Denominator-derived factors are load-balanced only over the
/// derivative-created copies of that owner. Momentum-routing signs are restored
/// while mapping the numerator, and candidate sets belonging to distinct source
/// owners remain disjoint. Completed UV factors instead share their certified
/// denominator-class pool, without changing physical or soft-carrier rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum EnergyCandidateFamily {
    Unprovenanced(EdgeIndex),
    Fixed(EdgeIndex),
    DenominatorDerived(EdgeIndex),
    UvClass(UvDenominatorClassId),
}

impl EnergyCandidateFamily {
    fn reference(self) -> EnergyReference {
        match self {
            Self::Unprovenanced(edge) | Self::Fixed(edge) | Self::DenominatorDerived(edge) => {
                EnergyReference::Physical(edge)
            }
            Self::UvClass(class) => EnergyReference::UvClass(class),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EquivalentEnergyCandidates {
    by_family: BTreeMap<EnergyCandidateFamily, Vec<usize>>,
    // Exact complete wrappers rewritten by the source's certified affine lift.
    // Only these may jointly reference several fixed base occurrences.
    pub(crate) fixed_affine_blocks: BTreeSet<Atom>,
}

impl EquivalentEnergyCandidates {
    pub(crate) fn add_uv_classes(
        &mut self,
        groups: impl IntoIterator<Item = (UvDenominatorClassId, Vec<usize>)>,
    ) -> Result<(), EnergyPowerAnalysisError> {
        for (class, mut occurrences) in groups {
            occurrences.sort_unstable();
            if occurrences.is_empty()
                || occurrences.windows(2).any(|pair| pair[0] == pair[1])
                || self
                    .by_family
                    .contains_key(&EnergyCandidateFamily::UvClass(class))
                || self.by_family.iter().any(|(family, existing)| {
                    matches!(family, EnergyCandidateFamily::UvClass(_))
                        && existing
                            .iter()
                            .any(|occurrence| occurrences.binary_search(occurrence).is_ok())
                })
            {
                return Err(EnergyPowerAnalysisError::InvalidClassEnergyCandidates {
                    class: class.0,
                    detail: "empty, duplicate or overlapping class pool".into(),
                });
            }
            // Physical restrictions may refer to an occurrence in this class;
            // distinct certified pole classes themselves must remain disjoint.
            self.by_family
                .insert(EnergyCandidateFamily::UvClass(class), occurrences);
        }
        Ok(())
    }
    /// Build candidate sets from the serial denominator occurrences certified
    /// by source-edge topology reconstruction. Algebraic energy equality is
    /// not rediscovered here: members already share one immutable source owner.
    pub(crate) fn try_from_source_occurrences(
        groups: impl IntoIterator<Item = (EdgeIndex, Vec<usize>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        let mut by_family = BTreeMap::new();
        let mut candidate_owners = BTreeMap::<usize, EdgeIndex>::new();

        for (edge, mut candidates) in groups {
            if by_family.contains_key(&EnergyCandidateFamily::Fixed(edge)) {
                return Err(
                    EnergyPowerAnalysisError::DuplicatePhysicalEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            candidates.sort_unstable();
            if candidates.is_empty() {
                return Err(
                    EnergyPowerAnalysisError::EmptyEquivalentEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            for pair in candidates.windows(2) {
                if pair[0] == pair[1] {
                    return Err(
                        EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate {
                            edge: edge.into(),
                            candidate: pair[0],
                        },
                    );
                }
            }
            for candidate in &candidates {
                if let Some(first_edge) = candidate_owners.insert(*candidate, edge) {
                    return Err(
                        EnergyPowerAnalysisError::OverlappingEquivalentEnergyCandidates {
                            candidate: *candidate,
                            first_edge: first_edge.into(),
                            second_edge: edge.into(),
                        },
                    );
                }
            }

            // Serial copies can reduce the maximal denominator-derived degree.
            // The lowest canonical occurrence supplies the base when
            // the caller has not explicitly partitioned the serial copies.
            by_family.insert(
                EnergyCandidateFamily::Unprovenanced(edge),
                vec![candidates[0]],
            );
            by_family.insert(EnergyCandidateFamily::Fixed(edge), vec![candidates[0]]);
            by_family.insert(EnergyCandidateFamily::DenominatorDerived(edge), candidates);
        }

        Ok(Self {
            by_family,
            fixed_affine_blocks: BTreeSet::new(),
        })
    }

    /// Build production candidate sets from one retained source occurrence and
    /// its derivative-created serial copies. The base is explicit because
    /// exact-topology canonicalization may place it after a generated copy in
    /// the occurrence-local energy namespace.
    pub(crate) fn try_from_partitioned_source_occurrences(
        groups: impl IntoIterator<Item = (EdgeIndex, usize, Vec<usize>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        let mut groups = groups.into_iter().collect::<Vec<_>>();
        let mut candidates = Self::try_from_source_occurrences(groups.iter().map(
            |(edge, base, derivative_copies)| {
                (
                    *edge,
                    std::iter::once(*base)
                        .chain(derivative_copies.iter().copied())
                        .collect(),
                )
            },
        ))?;
        for (edge, base, derivative_copies) in &mut groups {
            derivative_copies.sort_unstable();
            candidates
                .by_family
                .insert(EnergyCandidateFamily::Unprovenanced(*edge), vec![*base]);
            candidates
                .by_family
                .insert(EnergyCandidateFamily::Fixed(*edge), vec![*base]);
            if derivative_copies.is_empty() {
                candidates
                    .by_family
                    .remove(&EnergyCandidateFamily::DenominatorDerived(*edge));
            } else {
                candidates.by_family.insert(
                    EnergyCandidateFamily::DenominatorDerived(*edge),
                    derivative_copies.clone(),
                );
            }
        }
        Ok(candidates)
    }

    fn get(&self, family: EnergyCandidateFamily) -> Option<&[usize]> {
        self.by_family.get(&family).map(Vec::as_slice)
    }
}

// Each factor retains its own signed physical/class-to-occurrence assignment.
pub(crate) type FactorEnergyAssignments = BTreeMap<EnergyReference, usize>;

/// One immutable assignment of factor-local energy dependencies to certified
/// exact occurrences. The same tree supplies generation bounds and numerator
/// substitutions; repeated assigned cycles remain compressed.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct EnergyPowerAssignmentPlan {
    expression: PlannedEnergyExpression,
    energy_degree_bounds: Vec<(usize, usize)>,
    placement_family: Option<(EnergyCandidateFamily, Vec<usize>)>,
}

impl EnergyPowerAssignmentPlan {
    pub(crate) fn visit_factors<E>(
        &self,
        mut visit: impl FnMut(&Atom, &FactorEnergyAssignments) -> Result<(), E>,
    ) -> Result<(), E> {
        self.expression.visit_factors(&mut visit)
    }

    pub(crate) fn certify_bounds(&self) -> Result<(), EnergyPowerAnalysisError> {
        let actual = self
            .expression
            .exact_degrees()?
            .into_iter()
            .collect::<Vec<_>>();
        if actual != self.energy_degree_bounds {
            return Err(EnergyPowerAnalysisError::InvalidAssignmentCertificate {
                detail: "declared capacities differ from the immutable assigned expression".into(),
            });
        }
        Ok(())
    }

    fn certify_roundtrip(&self, original: &Atom) -> Result<(), EnergyPowerAnalysisError> {
        let collapsed = self
            .map_factors(|factor, _| Ok::<_, std::convert::Infallible>(factor.clone()))
            .unwrap();
        if collapsed != *original
            && collapsed.expand_num().collect_factors() != original.expand_num().collect_factors()
        {
            return Err(EnergyPowerAnalysisError::InvalidAssignmentCertificate {
                detail: "factorized assignment does not collapse to its analyzed numerator".into(),
            });
        }
        self.certify_bounds()
    }

    pub(crate) fn accounted_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.expression.accounted_bytes()
            + self.energy_degree_bounds.capacity() * std::mem::size_of::<(usize, usize)>()
            + self
                .placement_family
                .as_ref()
                .map_or(0, |(_, occurrences)| {
                    occurrences.capacity() * std::mem::size_of::<usize>()
                })
    }

    pub(crate) fn energy_degree_bounds(&self) -> &[(usize, usize)] {
        &self.energy_degree_bounds
    }

    pub(crate) fn expression(&self) -> &PlannedEnergyExpression {
        &self.expression
    }

    pub(crate) fn map_factors<E>(
        &self,
        mut map: impl FnMut(&Atom, &FactorEnergyAssignments) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        self.expression
            .map(&mut |_, expression, assignments| map(expression, assignments))
    }

    pub(crate) fn prepare_factors<E>(
        &self,
        mut map: impl FnMut(&Atom, &FactorEnergyAssignments) -> Result<Atom, E>,
    ) -> Result<PlannedEnergyExpression, E> {
        self.expression.prepare_factors(&mut map)
    }

    /// The caller supplies native CFF scores before requesting the third plan.
    /// Equal ordered bounds need only one trial even if leaf assignments differ.
    pub(crate) fn placement_challenger(
        &self,
        seen_bounds: &[Vec<(usize, usize)>],
    ) -> Result<Option<Self>, EnergyPowerAnalysisError> {
        let Some((family, occurrences)) = &self.placement_family else {
            return Ok(None);
        };
        for reverse in [true, false] {
            let mut expression = self.expression.clone();
            expression.permute_family(*family, occurrences, reverse);
            let bounds = expression.exact_degrees()?.into_iter().collect::<Vec<_>>();
            if !seen_bounds.contains(&bounds) {
                return Ok(Some(Self {
                    expression,
                    energy_degree_bounds: bounds,
                    placement_family: self.placement_family.clone(),
                }));
            }
        }
        Ok(None)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) enum PlannedEnergyExpression {
    Factor {
        id: usize,
        expression: Atom,
        degrees: EnergyPowerCapMap<EnergyReference>,
        families: BTreeMap<EnergyReference, EnergyCandidateFamily>,
        assignments: FactorEnergyAssignments,
    },
    Add(Vec<Self>),
    Mul(Vec<Self>),
    MultilinearFunction {
        symbol: Symbol,
        arguments: Vec<Self>,
    },
    Repeat {
        base: Box<Self>,
        exponent: usize,
    },
}

#[derive(Clone, Default)]
struct EnergyAllocationState {
    offsets: BTreeMap<EnergyCandidateFamily, usize>,
    loads: EnergyPowerCapMap<usize>,
}

impl EnergyAllocationState {
    fn cycle_key(
        &self,
        candidates: &EquivalentEnergyCandidates,
        packed: Option<EnergyCandidateFamily>,
        dependencies: &BTreeMap<EnergyCandidateFamily, bool>,
    ) -> Vec<(EnergyCandidateFamily, usize, Vec<usize>)> {
        dependencies
            .iter()
            .map(|(family, greedy)| {
                let occurrences = &candidates.by_family[family];
                let offset = self.offsets.get(family).copied().unwrap_or(0);
                let offset = if packed == Some(*family) {
                    offset.min(occurrences.len())
                } else {
                    offset % occurrences.len()
                };
                let loads = occurrences
                    .iter()
                    .filter(|_| *greedy)
                    .map(|occurrence| self.loads.degrees.get(occurrence).copied().unwrap_or(0))
                    .collect::<Vec<_>>();
                let minimum = loads.iter().copied().min().unwrap_or(0);
                (
                    *family,
                    offset,
                    loads.into_iter().map(|load| load - minimum).collect(),
                )
            })
            .collect()
    }

    fn advance_cycles(
        &mut self,
        previous: &Self,
        cycles: usize,
    ) -> Result<(), EnergyPowerAnalysisError> {
        for (family, offset) in &mut self.offsets {
            let delta = *offset - previous.offsets.get(family).copied().unwrap_or(0);
            *offset = offset
                .checked_add(
                    delta
                        .checked_mul(cycles)
                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?,
                )
                .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
        }
        for (occurrence, load) in &mut self.loads.degrees {
            let delta = *load - previous.loads.degrees.get(occurrence).copied().unwrap_or(0);
            *load = load
                .checked_add(
                    delta
                        .checked_mul(cycles)
                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?,
                )
                .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
        }
        Ok(())
    }
}

impl PlannedEnergyExpression {
    fn visit_factors<E>(
        &self,
        visit: &mut impl FnMut(&Atom, &FactorEnergyAssignments) -> Result<(), E>,
    ) -> Result<(), E> {
        match self {
            Self::Factor {
                expression,
                assignments,
                ..
            } => visit(expression, assignments),
            Self::Add(children)
            | Self::Mul(children)
            | Self::MultilinearFunction {
                arguments: children,
                ..
            } => children
                .iter()
                .try_for_each(|child| child.visit_factors(visit)),
            Self::Repeat { base, .. } => base.visit_factors(visit),
        }
    }

    pub(crate) fn accounted_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + match self {
                Self::Factor {
                    expression,
                    degrees,
                    families,
                    assignments,
                    ..
                } => {
                    expression.as_view().get_byte_size()
                        + degrees.degrees.len()
                            * (std::mem::size_of::<(EnergyReference, usize)>() + 32)
                        + families.len()
                            * (std::mem::size_of::<(EnergyReference, EnergyCandidateFamily)>() + 32)
                        + assignments.len() * (std::mem::size_of::<(EnergyReference, usize)>() + 32)
                }
                Self::Add(children)
                | Self::Mul(children)
                | Self::MultilinearFunction {
                    arguments: children,
                    ..
                } => {
                    children.capacity() * std::mem::size_of::<Self>()
                        + children.iter().map(Self::accounted_bytes).sum::<usize>()
                }
                Self::Repeat { base, .. } => base.accounted_bytes(),
            }
    }

    /// Select independently at factor leaves, adding ranks in products and
    /// taking their componentwise maximum across sums. The rank Pareto frontier
    /// is exact; no numerator polynomial or Cartesian product of its monomials is
    /// ever constructed.
    fn soft_momentum_assignments(
        &self,
        alternatives: &BTreeMap<Atom, Vec<(Atom, EnergyPowerCapMap<EnergyReference>)>>,
    ) -> Result<
        Vec<(EnergyPowerCapMap<EnergyReference>, BTreeMap<usize, usize>)>,
        EnergyPowerAnalysisError,
    > {
        match self {
            Self::Factor {
                id,
                expression,
                degrees,
                ..
            } => Ok(alternatives
                .get(expression)
                .map(|choices| {
                    choices
                        .iter()
                        .enumerate()
                        .map(|(choice, (_, degrees))| {
                            (degrees.clone(), BTreeMap::from([(*id, choice)]))
                        })
                        .collect()
                })
                .unwrap_or_else(|| vec![(degrees.clone(), BTreeMap::new())])),
            Self::Add(children)
            | Self::Mul(children)
            | Self::MultilinearFunction {
                arguments: children,
                ..
            } => {
                let additive = matches!(self, Self::Add(_));
                children.iter().try_fold(
                    vec![(EnergyPowerCapMap::default(), BTreeMap::new())],
                    |states, child| {
                        let right = child.soft_momentum_assignments(alternatives)?;
                        let mut unique = BTreeMap::new();
                        for (left_degrees, left_assignments) in &states {
                            for (right_degrees, right_assignments) in &right {
                                let mut degrees = left_degrees.clone();
                                if additive {
                                    degrees.max_assign(right_degrees.clone());
                                } else {
                                    degrees.add_assign(right_degrees.clone())?;
                                }
                                unique.entry(degrees.degrees).or_insert_with(|| {
                                    let mut assignments = left_assignments.clone();
                                    assignments.extend(right_assignments);
                                    assignments
                                });
                            }
                        }
                        Ok(unique
                            .iter()
                            .filter(|(degrees, _)| {
                                !unique.keys().any(|other| {
                                    other != *degrees
                                        && other.iter().all(|(edge, degree)| {
                                            *degree <= degrees.get(edge).copied().unwrap_or(0)
                                        })
                                })
                            })
                            .map(|(degrees, assignments)| {
                                (
                                    EnergyPowerCapMap {
                                        degrees: degrees.clone(),
                                    },
                                    assignments.clone(),
                                )
                            })
                            .collect())
                    },
                )
            }
            Self::Repeat { base, exponent } => {
                let mut alternatives = base.soft_momentum_assignments(alternatives)?;
                for (degrees, _) in &mut alternatives {
                    degrees.scale(*exponent)?;
                }
                Ok(alternatives)
            }
        }
    }

    fn degrees(&self) -> Result<EnergyPowerCapMap<EnergyReference>, EnergyPowerAnalysisError> {
        match self {
            Self::Factor { degrees, .. } => Ok(degrees.clone()),
            Self::Add(terms) => {
                terms
                    .iter()
                    .try_fold(EnergyPowerCapMap::default(), |mut sum, term| {
                        sum.max_assign(term.degrees()?);
                        Ok(sum)
                    })
            }
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => factors
                .iter()
                .try_fold(EnergyPowerCapMap::default(), |mut product, factor| {
                    product.add_assign(factor.degrees()?)?;
                    Ok(product)
                }),
            Self::Repeat { base, exponent } => {
                let mut degrees = base.degrees()?;
                degrees.scale(*exponent)?;
                Ok(degrees)
            }
        }
    }

    fn allocate(
        &self,
        candidates: &EquivalentEnergyCandidates,
        state: &mut EnergyAllocationState,
        packed: Option<EnergyCandidateFamily>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        match self {
            Self::Factor {
                degrees, families, ..
            } => {
                let mut result = self.clone();
                let Self::Factor { assignments, .. } = &mut result else {
                    unreachable!()
                };
                for (reference, degree) in degrees.iter() {
                    let family = families[&reference];
                    let eligible = candidates
                        .get(family)
                        .ok_or_else(|| reference.missing_candidates(degree))?;
                    let offset = state.offsets.entry(family).or_default();
                    let candidate =
                        if degree == 1 && matches!(family, EnergyCandidateFamily::UvClass(_)) {
                            let position = if packed == Some(family) && *offset >= eligible.len() {
                                0
                            } else {
                                *offset % eligible.len()
                            };
                            eligible[position]
                        } else {
                            // Fixed families have one candidate. Indivisible blocks
                            // use least-loaded whole placement, never factorization.
                            *eligible
                                .iter()
                                .min_by_key(|candidate| {
                                    (
                                        state.loads.degrees.get(candidate).copied().unwrap_or(0),
                                        **candidate,
                                    )
                                })
                                .unwrap()
                        };
                    // Assignment insertion is production work: putting it
                    // inside debug_assert! would remove it from release builds.
                    assignments.insert(reference, candidate);
                    *offset = offset
                        .checked_add(degree)
                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
                    let load = state.loads.degrees.entry(candidate).or_default();
                    *load = load
                        .checked_add(degree)
                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
                }
                Ok(result)
            }
            Self::Add(terms) => {
                // Branches share their starting interval and combine maxima;
                // products concatenate intervals without expanding any sums.
                let initial = state.clone();
                let mut result = Vec::with_capacity(terms.len());
                for term in terms {
                    let mut branch = initial.clone();
                    result.push(term.allocate(candidates, &mut branch, packed)?);
                    state.loads.max_assign(branch.loads);
                    for (family, offset) in branch.offsets {
                        let current = state.offsets.entry(family).or_default();
                        *current = (*current).max(offset);
                    }
                }
                Ok(Self::Add(result))
            }
            Self::Mul(factors) => Ok(Self::Mul(
                factors
                    .iter()
                    .map(|factor| factor.allocate(candidates, state, packed))
                    .collect::<Result<_, _>>()?,
            )),
            Self::MultilinearFunction { symbol, arguments } => Ok(Self::MultilinearFunction {
                symbol: *symbol,
                arguments: arguments
                    .iter()
                    .map(|argument| argument.allocate(candidates, state, packed))
                    .collect::<Result<_, _>>()?,
            }),
            Self::Repeat { base, exponent } => {
                let mut dependencies = BTreeMap::new();
                base.allocation_dependencies(&mut dependencies);
                for family in dependencies.keys() {
                    if candidates.get(*family).is_none() {
                        return Err(family
                            .reference()
                            .missing_candidates(base.family_degrees()?.degrees[family]));
                    }
                }
                let mut contexts = BTreeMap::new();
                let mut factors = Vec::new();
                let mut repetition = 0;
                while repetition < *exponent {
                    let key = state.cycle_key(candidates, packed, &dependencies);
                    if let Some((previous_repetition, start, previous_state)) = contexts.get(&key) {
                        let period = repetition - previous_repetition;
                        let extra_cycles = (*exponent - repetition) / period;
                        if extra_cycles > 0 {
                            // Equal cyclic offsets and relative greedy loads
                            // determine the same next assignment. Whole cycles
                            // can therefore retain a factorized integer power.
                            let cycle = factors.split_off(*start);
                            factors.push(Self::Repeat {
                                base: Box::new(Self::Mul(cycle)),
                                exponent: extra_cycles + 1,
                            });
                            state.advance_cycles(previous_state, extra_cycles)?;
                            repetition += extra_cycles * period;
                            contexts.clear();
                            continue;
                        }
                    } else if contexts.len() < 64 {
                        contexts.insert(key, (repetition, factors.len(), state.clone()));
                    }
                    factors.push(base.allocate(candidates, state, packed)?);
                    repetition += 1;
                }
                Ok(Self::Mul(factors))
            }
        }
    }

    fn allocation_dependencies(&self, dependencies: &mut BTreeMap<EnergyCandidateFamily, bool>) {
        match self {
            Self::Factor {
                degrees, families, ..
            } => {
                for (reference, degree) in degrees.iter() {
                    let family = families[&reference];
                    *dependencies.entry(family).or_default() |=
                        degree != 1 || !matches!(family, EnergyCandidateFamily::UvClass(_));
                }
            }
            Self::Add(children)
            | Self::Mul(children)
            | Self::MultilinearFunction {
                arguments: children,
                ..
            } => {
                for child in children {
                    child.allocation_dependencies(dependencies);
                }
            }
            Self::Repeat { base, .. } => base.allocation_dependencies(dependencies),
        }
    }

    fn family_degrees(
        &self,
    ) -> Result<EnergyPowerCapMap<EnergyCandidateFamily>, EnergyPowerAnalysisError> {
        match self {
            Self::Factor {
                degrees, families, ..
            } => Ok(EnergyPowerCapMap {
                degrees: degrees
                    .iter()
                    .map(|(reference, degree)| (families[&reference], degree))
                    .collect(),
            }),
            Self::Add(children)
            | Self::Mul(children)
            | Self::MultilinearFunction {
                arguments: children,
                ..
            } => {
                let mut result = EnergyPowerCapMap::default();
                for child in children {
                    let degrees = child.family_degrees()?;
                    if matches!(self, Self::Add(_)) {
                        result.max_assign(degrees);
                    } else {
                        result.add_assign(degrees)?;
                    }
                }
                Ok(result)
            }
            Self::Repeat { base, exponent } => {
                let mut degrees = base.family_degrees()?;
                degrees.scale(*exponent)?;
                Ok(degrees)
            }
        }
    }

    fn permute_family(
        &mut self,
        family: EnergyCandidateFamily,
        occurrences: &[usize],
        reverse: bool,
    ) {
        match self {
            Self::Factor {
                families,
                assignments,
                ..
            } => {
                for (reference, assigned) in assignments {
                    if families[reference] == family {
                        let position = occurrences
                            .binary_search(assigned)
                            .expect("assignment belongs to its certified pool");
                        *assigned = occurrences[if reverse {
                            occurrences.len() - position - 1
                        } else {
                            (position + 1) % occurrences.len()
                        }];
                    }
                }
            }
            Self::Add(children)
            | Self::Mul(children)
            | Self::MultilinearFunction {
                arguments: children,
                ..
            } => {
                for child in children {
                    child.permute_family(family, occurrences, reverse);
                }
            }
            Self::Repeat { base, .. } => base.permute_family(family, occurrences, reverse),
        }
    }

    pub(crate) fn prepare_factors<E>(
        &self,
        map: &mut impl FnMut(&Atom, &FactorEnergyAssignments) -> Result<Atom, E>,
    ) -> Result<Self, E> {
        Ok(match self {
            Self::Factor {
                expression,
                assignments,
                ..
            } => {
                let mut result = self.clone();
                let Self::Factor {
                    expression: prepared,
                    ..
                } = &mut result
                else {
                    unreachable!()
                };
                *prepared = map(expression, assignments)?;
                result
            }
            Self::Add(children) => Self::Add(
                children
                    .iter()
                    .map(|child| child.prepare_factors(map))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Mul(children) => Self::Mul(
                children
                    .iter()
                    .map(|child| child.prepare_factors(map))
                    .collect::<Result<_, _>>()?,
            ),
            Self::MultilinearFunction { symbol, arguments } => Self::MultilinearFunction {
                symbol: *symbol,
                arguments: arguments
                    .iter()
                    .map(|child| child.prepare_factors(map))
                    .collect::<Result<_, _>>()?,
            },
            Self::Repeat { base, exponent } => Self::Repeat {
                base: Box::new(base.prepare_factors(map)?),
                exponent: *exponent,
            },
        })
    }

    fn exact_degrees(&self) -> Result<BTreeMap<usize, usize>, EnergyPowerAnalysisError> {
        match self {
            Self::Factor {
                id,
                expression,
                degrees,
                assignments,
                ..
            } => {
                let mut exact = BTreeMap::new();
                for (edge, degree) in degrees.iter() {
                    let candidate = assignments.get(&edge).copied().ok_or_else(|| {
                        EnergyPowerAnalysisError::MissingFactorEnergyAssignment {
                            factor_id: *id,
                            expression: expression.log_print(None),
                            edge: match edge {
                                EnergyReference::Physical(edge) => edge.into(),
                                EnergyReference::UvClass(class) => class.0,
                            },
                        }
                    })?;
                    let slot = exact.entry(candidate).or_insert(0usize);
                    *slot = slot
                        .checked_add(degree)
                        .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
                }
                Ok(exact)
            }
            Self::Add(terms) => {
                let mut envelope = BTreeMap::new();
                for term in terms {
                    for (candidate, degree) in term.exact_degrees()? {
                        let slot = envelope.entry(candidate).or_insert(0);
                        *slot = (*slot).max(degree);
                    }
                }
                Ok(envelope)
            }
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => {
                let mut product = BTreeMap::new();
                for factor in factors {
                    for (candidate, degree) in factor.exact_degrees()? {
                        let slot = product.entry(candidate).or_insert(0usize);
                        *slot = slot
                            .checked_add(degree)
                            .ok_or(EnergyPowerAnalysisError::EnergyDegreeOverflow)?;
                    }
                }
                Ok(product)
            }
            Self::Repeat { base, exponent } => {
                let mut degrees = EnergyPowerCapMap {
                    degrees: base.exact_degrees()?,
                };
                degrees.scale(*exponent)?;
                Ok(degrees.degrees)
            }
        }
    }

    pub(crate) fn map<E>(
        &self,
        map: &mut impl FnMut(usize, &Atom, &FactorEnergyAssignments) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        match self {
            Self::Factor {
                id,
                expression,
                assignments,
                ..
            } => map(*id, expression, assignments),
            Self::Add(terms) => terms
                .iter()
                .try_fold(Atom::Zero, |sum, term| Ok(sum + term.map(map)?)),
            Self::Mul(factors) => factors.iter().try_fold(Atom::one(), |product, factor| {
                Ok(product * factor.map(map)?)
            }),
            Self::MultilinearFunction { symbol, arguments } => {
                let mut builder = FunctionBuilder::new(*symbol);
                for argument in arguments {
                    let mapped = argument.map(map)?;
                    builder = builder.add_arg(mapped.as_view());
                }
                Ok(builder.finish())
            }
            Self::Repeat { base, exponent } => Ok(base.map(map)?.pow(Atom::num(*exponent as u64))),
        }
    }
}

#[derive(Debug, Error)]
pub enum EnergyPowerAnalysisError {
    #[error("invalid exact-energy assignment certificate: {detail}")]
    InvalidAssignmentCertificate { detail: String },
    #[error("invalid exact-energy candidates for UV class {class}: {detail}")]
    InvalidClassEnergyCandidates { class: usize, detail: String },
    #[error("numerator energy degree exceeds the supported integer range")]
    EnergyDegreeOverflow,
    #[error("physical EMR edge {edge} has more than one exact-energy candidate set")]
    DuplicatePhysicalEnergyCandidateSet { edge: usize },
    #[error("physical EMR edge {edge} has an empty exact-energy candidate set")]
    EmptyEquivalentEnergyCandidateSet { edge: usize },
    #[error(
        "physical EMR edge {edge} repeats exact-energy candidate {candidate} in one candidate set"
    )]
    DuplicateEquivalentEnergyCandidate { edge: usize, candidate: usize },
    #[error(
        "exact-energy candidate {candidate} is shared by physical EMR edges {first_edge} and {second_edge}; minimax assignment requires disjoint certified candidate sets"
    )]
    OverlappingEquivalentEnergyCandidates {
        candidate: usize,
        first_edge: usize,
        second_edge: usize,
    },
    #[error(
        "physical EMR energy-power degree {degree} for edge {edge} has no certified equivalent exact-energy candidates"
    )]
    MissingEquivalentEnergyCandidates { edge: usize, degree: usize },
    #[error(
        "factor {factor_id} (`{expression}`) depends on physical EMR edge {edge} but has no exact-energy assignment"
    )]
    MissingFactorEnergyAssignment {
        factor_id: usize,
        expression: String,
        edge: usize,
    },
    #[error(
        "positive denominator wrapper `{expression}` is owned by edge {owner} but its active energy depends on edge {dependent}"
    )]
    NonlocalPositiveDenominator {
        expression: String,
        owner: usize,
        dependent: usize,
    },
    #[error(
        "positive denominator wrapper `{expression}` mixes incompatible energy-provenance families for owner {owner}"
    )]
    MixedPositiveDenominatorProvenance { expression: String, owner: usize },
    #[error(
        "physical-EMR energy-power analysis cannot assign loop-basis expression `{expression}` to a source edge; normalize it to Q(edge, index) before CFF analysis"
    )]
    LoopMomentumInPhysicalEmrAnalysis { expression: String },
    #[error("loop momentum K({loop_index}, ...) has no carrier edge in this loop momentum basis")]
    MissingLoopCarrierEdge { loop_index: usize },
    #[error(
        "energy-power analysis expected an integer edge id in Q(edge, index), found `{argument}`"
    )]
    InvalidEmrEdgeArgument { argument: String },
    #[error("energy-power analysis found malformed UV momentum provenance `{argument}`")]
    InvalidUvMomentumProvenance { argument: String },
    #[error(
        "energy-power analysis expected an integer loop id in K(loop, index), found `{argument}`"
    )]
    InvalidLoopMomentumArgument { argument: String },
    #[error(
        "non-polynomial energy dependence encountered in numerator subexpression `{expression}` with exponent `{exponent}`"
    )]
    NonPolynomialEnergyPower {
        expression: String,
        exponent: String,
    },
    #[error(
        "energy dependence encountered in function without declared multilinearity: `{expression}`"
    )]
    OpaqueEnergyFunction { expression: String },
}

pub struct EnergyPowerAnalyzer {
    loop_edges: Option<Vec<EdgeIndex>>,
    internal_edges: Option<BTreeSet<EdgeIndex>>,
    active_uv_classes: Option<BTreeSet<UvDenominatorClassId>>,
    minkowski_symbol: Symbol,
}

#[derive(Clone, Copy)]
enum EnergyPowerAnalysisMode {
    StrictPolynomial,
    ConservativeUpperBound,
}

#[derive(Clone, Copy)]
enum LorentzIndexKind {
    Temporal,
    Spatial,
    Abstract,
}

impl EnergyPowerAnalyzer {
    pub fn new(loop_edges: impl IntoIterator<Item = EdgeIndex>) -> Self {
        Self {
            loop_edges: Some(loop_edges.into_iter().collect()),
            internal_edges: None,
            active_uv_classes: None,
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn with_internal_edges(
        loop_edges: impl IntoIterator<Item = EdgeIndex>,
        internal_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Self {
        Self {
            loop_edges: Some(loop_edges.into_iter().collect()),
            internal_edges: Some(internal_edges.into_iter().collect()),
            active_uv_classes: None,
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn for_physical_emr_edges(internal_edges: impl IntoIterator<Item = EdgeIndex>) -> Self {
        Self {
            loop_edges: None,
            internal_edges: Some(internal_edges.into_iter().collect()),
            active_uv_classes: None,
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn analyze_atom(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_physical_atom(expression, EnergyPowerAnalysisMode::StrictPolynomial)
    }

    pub fn analyze_atom_upper_bound(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_physical_atom(expression, EnergyPowerAnalysisMode::ConservativeUpperBound)
    }

    fn analyze_physical_atom(
        &self,
        expression: &Atom,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        Ok(EnergyPowerCapMap {
            degrees: self
                .analyze_view(expression.as_view(), mode)?
                .iter()
                .map(|(reference, degree)| match reference {
                    EnergyReference::Physical(edge) => Ok((edge, degree)),
                    EnergyReference::UvClass(class) => {
                        Err(EnergyPowerAnalysisError::InvalidClassEnergyCandidates {
                            class: class.0,
                            detail: "canonical class requires reference-aware analysis".into(),
                        })
                    }
                })
                .collect::<Result<_, _>>()?,
        })
    }

    pub(crate) fn analyze_reference_atom(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap<EnergyReference>, EnergyPowerAnalysisError> {
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::StrictPolynomial,
        )
    }

    #[cfg(test)]
    pub(crate) fn plan_atom_assignment(
        &self,
        expression: &Atom,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        Ok(self
            .plan_atom_assignment_proposals(expression, candidates)?
            .remove(0))
    }

    pub(crate) fn plan_atom_assignment_proposals(
        &self,
        expression: &Atom,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<Vec<EnergyPowerAssignmentPlan>, EnergyPowerAnalysisError> {
        let started = std::time::Instant::now();
        let expected_degrees = self.analyze_reference_atom(expression)?;
        let mut next_factor_id = 0;
        let planned = self.plan_view(
            expression.as_view(),
            &mut next_factor_id,
            true,
            &candidates.fixed_affine_blocks,
        )?;
        debug_assert_eq!(planned.degrees()?, expected_degrees);
        let family_degrees = planned.family_degrees()?;
        for (family, degree) in family_degrees.iter() {
            if candidates.get(family).is_none() {
                return Err(family.reference().missing_candidates(degree));
            }
        }
        let baseline = planned.allocate(candidates, &mut EnergyAllocationState::default(), None)?;
        let baseline_bounds = baseline.exact_degrees()?.into_iter().collect::<Vec<_>>();
        let placement_family = family_degrees
            .iter()
            .filter_map(|(family, degree)| {
                let occurrences = candidates.get(family)?;
                if occurrences.len() < 2 {
                    return None;
                }
                let maximum = occurrences
                    .iter()
                    .map(|occurrence| {
                        baseline_bounds
                            .iter()
                            .find(|(id, _)| id == occurrence)
                            .map_or(0, |(_, bound)| *bound)
                    })
                    .max()
                    .unwrap_or(0);
                (maximum > 1).then_some((
                    std::cmp::Reverse(maximum),
                    std::cmp::Reverse(degree.saturating_sub(occurrences.len())),
                    family,
                    occurrences.to_vec(),
                ))
            })
            .min()
            .map(|(_, _, family, occurrences)| (family, occurrences));
        // The old Cartesian/Pareto frontier is replaced by one baseline and a
        // single shape change. Native row counts select the incumbent before
        // its bounded placement reversal; rank never substitutes for counts.
        let mut plans = vec![EnergyPowerAssignmentPlan {
            expression: baseline,
            energy_degree_bounds: baseline_bounds,
            placement_family: placement_family.clone(),
        }];
        if let Some((family @ EnergyCandidateFamily::UvClass(_), _)) = placement_family {
            let packed = planned.allocate(
                candidates,
                &mut EnergyAllocationState::default(),
                Some(family),
            )?;
            let bounds = packed.exact_degrees()?.into_iter().collect::<Vec<_>>();
            if bounds != plans[0].energy_degree_bounds {
                plans.push(EnergyPowerAssignmentPlan {
                    expression: packed,
                    energy_degree_bounds: bounds,
                    placement_family: plans[0].placement_family.clone(),
                });
            }
        }
        let allocation_elapsed = started.elapsed();
        let certification_started = std::time::Instant::now();
        for plan in &plans {
            plan.certify_roundtrip(expression)?;
        }
        crate::debug_tags!(#generation, #uv, #local, #four_d, #cff, #profile;
            stage = "energy_assignment_preparation",
            allocation_ms = allocation_elapsed.as_secs_f64() * 1000.0,
            structural_certificate_ms = certification_started.elapsed().as_secs_f64() * 1000.0,
            proposals = plans.len(),
            "Prepared and certified factorized energy assignments"
        );
        Ok(plans)
    }

    fn plan_view(
        &self,
        expression: AtomView<'_>,
        next_factor_id: &mut usize,
        compress_powers: bool,
        fixed_affine_blocks: &BTreeSet<Atom>,
    ) -> Result<PlannedEnergyExpression, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Add(add) => Ok(PlannedEnergyExpression::Add(
                add.iter()
                    .map(|term| {
                        self.plan_view(term, next_factor_id, compress_powers, fixed_affine_blocks)
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            )),
            AtomView::Mul(mul) => Ok(PlannedEnergyExpression::Mul(
                mul.iter()
                    .map(|factor| {
                        self.plan_view(factor, next_factor_id, compress_powers, fixed_affine_blocks)
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            )),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let base_degrees =
                    self.analyze_view(base, EnergyPowerAnalysisMode::StrictPolynomial)?;
                if base_degrees.is_empty() {
                    return self.plan_factor(expression, next_factor_id);
                }
                if !self
                    .analyze_view(exponent, EnergyPowerAnalysisMode::StrictPolynomial)?
                    .is_empty()
                {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                }
                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                };
                if exponent < 0 {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_string(),
                    });
                }
                if exponent == 0 {
                    return self.plan_factor(expression, next_factor_id);
                }
                if compress_powers {
                    return Ok(PlannedEnergyExpression::Repeat {
                        base: Box::new(self.plan_view(
                            base,
                            next_factor_id,
                            true,
                            fixed_affine_blocks,
                        )?),
                        exponent: usize::try_from(exponent)
                            .map_err(|_| EnergyPowerAnalysisError::EnergyDegreeOverflow)?,
                    });
                }
                Ok(PlannedEnergyExpression::Mul(
                    (0..exponent)
                        .map(|_| {
                            self.plan_view(
                                base,
                                next_factor_id,
                                compress_powers,
                                fixed_affine_blocks,
                            )
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                ))
            }
            AtomView::Fun(function)
                if function.get_nargs() == 4 && function.get_symbol() == GS.den =>
            {
                let owner = self.emr_candidate_family(function.get(0))?.reference();
                let diagnostic_owner = match owner {
                    EnergyReference::Physical(edge) => edge.into(),
                    EnergyReference::UvClass(class) => class.0,
                };
                let expression = expression.to_owned();
                let degrees =
                    self.analyze_view(function.get(3), EnergyPowerAnalysisMode::StrictPolynomial)?;
                let mut value_families = BTreeSet::new();
                let mut family_error = None;
                let _ = function.get(3).to_owned().replace_map(|view, _, _| {
                    if family_error.is_some() {
                        return;
                    }
                    let AtomView::Fun(momentum) = view else {
                        return;
                    };
                    if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() != 2 {
                        return;
                    }
                    let family = match self.emr_candidate_family(momentum.get(0)) {
                        Ok(family) => family,
                        Err(error) => {
                            family_error = Some(error);
                            return;
                        }
                    };
                    if self.is_internal_reference(family.reference())
                        && !self
                            .component_degree(family.reference(), momentum.get(1))
                            .is_empty()
                    {
                        value_families.insert(family);
                    }
                });
                if let Some(error) = family_error {
                    return Err(error);
                }
                let fixed_affine = fixed_affine_blocks.contains(&expression);
                if (!fixed_affine && value_families.len() > 1)
                    || (fixed_affine
                        && value_families
                            .iter()
                            .any(|family| !matches!(family, EnergyCandidateFamily::Fixed(_))))
                {
                    return Err(
                        EnergyPowerAnalysisError::MixedPositiveDenominatorProvenance {
                            expression: expression.log_print(None),
                            owner: diagnostic_owner,
                        },
                    );
                }
                let mut families = BTreeMap::new();
                for (dependent, _) in degrees.iter() {
                    if !fixed_affine && dependent != owner {
                        return Err(EnergyPowerAnalysisError::NonlocalPositiveDenominator {
                            expression: expression.log_print(None),
                            owner: diagnostic_owner,
                            dependent: match dependent {
                                EnergyReference::Physical(edge) => edge.into(),
                                EnergyReference::UvClass(class) => class.0,
                            },
                        });
                    }
                    let family = value_families
                        .iter()
                        .copied()
                        .find(|family| family.reference() == dependent)
                        .ok_or_else(|| EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                            argument: expression.log_print(None),
                        })?;
                    families.insert(dependent, family);
                }
                let id = *next_factor_id;
                *next_factor_id += 1;
                // Conservatively keep a positive common-denominator
                // completion as one assignment unit. This gives its temporal
                // and spatial/mass pieces one occurrence frame without
                // expanding the wrapper. Separate wrappers may still be
                // load-balanced over serial copies of the same line, and the
                // mapper restores each selected routing sign. A source-certified
                // affine block stays equally opaque: its several fixed base
                // occurrences jointly supply its complete polynomial.
                Ok(PlannedEnergyExpression::Factor {
                    id,
                    expression,
                    degrees,
                    families,
                    assignments: BTreeMap::new(),
                })
            }
            AtomView::Fun(function)
                if function.get_symbol() == GS.dot
                    || (function.get_symbol().is_linear()
                        && function.get_symbol() != GS.emr_mom
                        && function.get_symbol() != GS.loop_mom
                        && function.get_symbol() != GS.emr_vec) =>
            {
                Ok(PlannedEnergyExpression::MultilinearFunction {
                    symbol: function.get_symbol(),
                    arguments: function
                        .iter()
                        .map(|argument| {
                            self.plan_view(
                                argument,
                                next_factor_id,
                                compress_powers,
                                fixed_affine_blocks,
                            )
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                })
            }
            _ => self.plan_factor(expression, next_factor_id),
        }
    }

    fn plan_factor(
        &self,
        expression: AtomView<'_>,
        next_factor_id: &mut usize,
    ) -> Result<PlannedEnergyExpression, EnergyPowerAnalysisError> {
        let id = *next_factor_id;
        *next_factor_id += 1;
        let degrees = self.analyze_view(expression, EnergyPowerAnalysisMode::StrictPolynomial)?;
        let mut families = BTreeMap::new();
        if !degrees.is_empty() {
            let AtomView::Fun(momentum) = expression else {
                unreachable!("factorized energy dependence must be one momentum component");
            };
            debug_assert_eq!(momentum.get_symbol(), GS.emr_mom);
            let family = self.emr_candidate_family(momentum.get(0))?;
            families.insert(family.reference(), family);
        }
        Ok(PlannedEnergyExpression::Factor {
            id,
            expression: expression.to_owned(),
            degrees,
            families,
            assignments: BTreeMap::new(),
        })
    }

    fn analyze_view(
        &self,
        expression: AtomView<'_>,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap<EnergyReference>, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Num(_) | AtomView::Var(_) => Ok(EnergyPowerCapMap::default()),
            AtomView::Add(add) => {
                let mut max_degree = EnergyPowerCapMap::default();
                for term in add.iter() {
                    max_degree.max_assign(self.analyze_view(term, mode)?);
                }
                Ok(max_degree)
            }
            AtomView::Mul(mul) => {
                let mut product_degree = EnergyPowerCapMap::default();
                for factor in mul.iter() {
                    product_degree.add_assign(self.analyze_view(factor, mode)?)?;
                }
                Ok(product_degree)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut base_degree = self.analyze_view(base, mode)?;
                if !self.analyze_view(exponent, mode)?.is_empty() {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                }
                if base_degree.is_empty() {
                    return Ok(base_degree);
                }

                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                };
                if exponent < 0 {
                    if matches!(mode, EnergyPowerAnalysisMode::ConservativeUpperBound) {
                        return Ok(EnergyPowerCapMap::default());
                    }
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_string(),
                    });
                }
                base_degree.scale(
                    usize::try_from(exponent)
                        .map_err(|_| EnergyPowerAnalysisError::EnergyDegreeOverflow)?,
                )?;
                Ok(base_degree)
            }
            AtomView::Fun(function) => {
                if function.get_symbol() == GS.emr_vec {
                    return Ok(EnergyPowerCapMap::default());
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.emr_mom {
                    let edge = self.emr_candidate_family(function.get(0))?.reference();
                    return Ok(if self.is_internal_reference(edge) {
                        self.component_degree(edge, function.get(1))
                    } else {
                        EnergyPowerCapMap::default()
                    });
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.loop_mom {
                    let Some(loop_edges) = &self.loop_edges else {
                        return Err(
                            EnergyPowerAnalysisError::LoopMomentumInPhysicalEmrAnalysis {
                                expression: Atom::from(function.to_owned()).log_print(None),
                            },
                        );
                    };
                    let loop_index = usize::try_from(function.get(0)).map_err(|_| {
                        EnergyPowerAnalysisError::InvalidLoopMomentumArgument {
                            argument: function.get(0).to_owned().log_print(None),
                        }
                    })?;
                    let edge = loop_edges
                        .get(loop_index)
                        .copied()
                        .ok_or(EnergyPowerAnalysisError::MissingLoopCarrierEdge { loop_index })?;
                    let edge = EnergyReference::Physical(edge);
                    return Ok(if self.is_internal_reference(edge) {
                        self.component_degree(edge, function.get(1))
                    } else {
                        EnergyPowerCapMap::default()
                    });
                }

                // `den(edge, momentum, mass, full_expr)` is a provenance
                // wrapper whose algebraic value is its final argument. A
                // positive wrapper can remain in a factorized numerator after
                // a local 4D Taylor expansion, so analyze that polynomial
                // without also counting its edge/momentum metadata.
                if function.get_nargs() == 4 && function.get_symbol() == GS.den {
                    return self.analyze_view(function.get(3), mode);
                }

                // Dot products are bilinear. Keep that composition explicit:
                // the Symbolica symbol is intentionally not globally declared
                // linear, while numerator energy degrees still add between its
                // two vector slots.
                if function.get_nargs() == 2 && function.get_symbol() == GS.dot {
                    let mut degree = self.analyze_view(function.get(0), mode)?;
                    degree.add_assign(self.analyze_view(function.get(1), mode)?)?;
                    return Ok(degree);
                }

                let argument_degrees = function
                    .iter()
                    .map(|argument| self.analyze_view(argument, mode))
                    .collect::<Result<Vec<_>, _>>()?;
                if argument_degrees.iter().any(|degree| !degree.is_empty())
                    && !function.get_symbol().is_linear()
                {
                    return Err(EnergyPowerAnalysisError::OpaqueEnergyFunction {
                        expression: Atom::from(function.to_owned()).log_print(None),
                    });
                }
                if function.get_symbol().is_linear() {
                    let mut degree = EnergyPowerCapMap::default();
                    for argument_degree in argument_degrees {
                        degree.add_assign(argument_degree)?;
                    }
                    Ok(degree)
                } else {
                    Ok(EnergyPowerCapMap::default())
                }
            }
        }
    }

    fn component_degree(
        &self,
        edge: EnergyReference,
        index: AtomView<'_>,
    ) -> EnergyPowerCapMap<EnergyReference> {
        match self.lorentz_index_kind(index) {
            LorentzIndexKind::Temporal | LorentzIndexKind::Abstract => {
                EnergyPowerCapMap::unit(edge)
            }
            LorentzIndexKind::Spatial => EnergyPowerCapMap::default(),
        }
    }

    fn emr_candidate_family(
        &self,
        argument: AtomView<'_>,
    ) -> Result<EnergyCandidateFamily, EnergyPowerAnalysisError> {
        if let Some(class) = GS.uv_class_data(argument) {
            return Ok(EnergyCandidateFamily::UvClass(class));
        }
        if let Ok(edge) = usize::try_from(argument) {
            return Ok(EnergyCandidateFamily::Unprovenanced(EdgeIndex(edge)));
        }

        if let Some((edge, role, _)) = GS.uv_momentum_provenance_data(argument) {
            return Ok(match role {
                UvMomentumProvenanceRole::DenominatorDerived => {
                    EnergyCandidateFamily::DenominatorDerived(edge)
                }
                UvMomentumProvenanceRole::TaylorFixed
                | UvMomentumProvenanceRole::PhysicalSourceFixed => {
                    EnergyCandidateFamily::Fixed(edge)
                }
                UvMomentumProvenanceRole::DenominatorDerivedSoft => {
                    EnergyCandidateFamily::Unprovenanced(edge)
                }
            });
        }
        match argument {
            AtomView::Fun(provenance) if provenance.get_symbol() == GS.uv_momentum_provenance => {
                Err(EnergyPowerAnalysisError::InvalidUvMomentumProvenance {
                    argument: argument.to_owned().log_print(None),
                })
            }
            _ => Err(EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                argument: argument.to_owned().log_print(None),
            }),
        }
    }

    fn lorentz_index_kind(&self, index: AtomView<'_>) -> LorentzIndexKind {
        if let Ok(index) = usize::try_from(index) {
            return if index == 0 {
                LorentzIndexKind::Temporal
            } else {
                LorentzIndexKind::Spatial
            };
        }

        let AtomView::Fun(function) = index else {
            return LorentzIndexKind::Spatial;
        };

        if function.get_symbol() == AIND_SYMBOLS.cind {
            return if function.get_nargs() == 1
                && usize::try_from(function.get(0)).is_ok_and(|index| index == 0)
            {
                LorentzIndexKind::Temporal
            } else {
                LorentzIndexKind::Spatial
            };
        }

        if function.get_symbol() == self.minkowski_symbol && matches!(function.get_nargs(), 1 | 2) {
            // `mink(4)` is the compact full-vector slot and `mink(4, n)` is an
            // abstract Lorentz-index label in the numerator algebra. Neither
            // selects a concrete spatial component; those use `cind(n)` and are
            // handled above. Both forms can contract onto the temporal
            // component and therefore contribute to the energy-power bound.
            return LorentzIndexKind::Abstract;
        }

        LorentzIndexKind::Spatial
    }

    fn is_internal_reference(&self, reference: EnergyReference) -> bool {
        match reference {
            EnergyReference::Physical(edge) => self
                .internal_edges
                .as_ref()
                .is_none_or(|edges| edges.contains(&edge)),
            EnergyReference::UvClass(class) => self
                .active_uv_classes
                .as_ref()
                .is_none_or(|classes| classes.contains(&class)),
        }
    }
}

impl Graph {
    #[cfg(test)]
    pub(crate) fn optimize_soft_momentum_routing(
        &self,
        numerator: &Atom,
        active_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> color_eyre::Result<Atom> {
        Ok(self
            .soft_momentum_routing_proposals(numerator, active_edges)?
            .remove(0))
    }

    /// Route only denominator-derived soft momenta once the untouched outer
    /// numerator is present. Original factors retain their literal owners.
    /// Every candidate is an exact off-shell linear identity in the production
    /// LMB, including the fixed external shift; on-shell energies play no role.
    pub(crate) fn soft_momentum_routing_proposals(
        &self,
        numerator: &Atom,
        active_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> color_eyre::Result<Vec<Atom>> {
        if !numerator.contains_symbol(GS.uv_momentum_provenance) {
            return Ok(vec![numerator.clone()]);
        }
        let mut soft_components = BTreeSet::new();
        let _ = numerator.replace_map(|view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() == 2
                && let Some((_, UvMomentumProvenanceRole::DenominatorDerivedSoft, _)) =
                    GS.uv_momentum_provenance_data(momentum.get(0))
            {
                soft_components.insert(view.to_owned());
                **output = view.to_owned();
            }
        });
        if soft_components.is_empty() {
            return Ok(vec![numerator.clone()]);
        }

        let active_edges = active_edges.into_iter().collect::<BTreeSet<_>>();
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(active_edges.iter().copied());
        let lmb = &self.loop_momentum_basis;
        let row = |edge: EdgeIndex| {
            lmb.edge_signatures[edge]
                .internal
                .iter()
                .map(|sign| match sign {
                    SignOrZero::Minus => -1,
                    SignOrZero::Zero => 0,
                    SignOrZero::Plus => 1,
                })
                .collect::<Vec<i64>>()
        };
        let rows = active_edges
            .iter()
            .map(|edge| (*edge, row(*edge)))
            .filter(|(_, row)| row.iter().any(|coefficient| *coefficient != 0))
            .collect::<Vec<_>>();
        let rank = rank_i64(&rows.iter().map(|(_, row)| row.clone()).collect::<Vec<_>>());
        let mut alternatives = BTreeMap::new();
        let mut planning_edges = active_edges.clone();
        for component in soft_components {
            let AtomView::Fun(momentum) = component.as_view() else {
                unreachable!()
            };
            let (owner, _, payload) = GS.uv_momentum_provenance_data(momentum.get(0)).unwrap();
            planning_edges.insert(owner);
            let original = GS.erase_uv_momentum_provenance(&component);
            // Soft tags name actual crown carriers. A nested Taylor operation
            // consumes its new hard part separately before reaching this stage.
            if payload
                != FunctionBuilder::new(GS.emr_mom)
                    .add_arg(usize::from(owner))
                    .finish()
            {
                return Err(eyre::eyre!(
                    "soft momentum provenance does not retain its literal crown carrier"
                ));
            }
            let target = row(owner);
            let indices = [momentum.get(1).to_owned()];
            let spatial = matches!(
                analyzer.lorentz_index_kind(momentum.get(1)),
                LorentzIndexKind::Spatial
            );
            let fixed_external = target.iter().all(|coefficient| *coefficient == 0);
            let mut candidates = BTreeSet::new();
            if spatial || fixed_external || active_edges.contains(&owner) {
                candidates.insert(original);
            }
            if spatial {
                // Concrete spatial components carry no energy rank and keep
                // their original routing.
            } else if fixed_external {
                candidates.insert(lmb.ext_atom(owner, GS.emr_mom, &indices, true));
            } else {
                // Every minimal support is linearly independent and extends to
                // a basis of the available row space. Enumerating those bases
                // therefore includes every nondominated energy support.
                for selected in rows.iter().combinations(rank) {
                    for columns in (0..target.len()).combinations(rank) {
                        let matrix = columns
                            .iter()
                            .map(|axis| {
                                selected
                                    .iter()
                                    .map(|(_, row)| Rational::from(row[*axis]))
                                    .collect()
                            })
                            .collect();
                        let rhs = columns
                            .iter()
                            .map(|axis| Rational::from(target[*axis]))
                            .collect();
                        let Some(coefficients) = solve_rational_system(matrix, rhs) else {
                            continue;
                        };
                        if !(0..target.len()).all(|axis| {
                            selected.iter().zip(&coefficients).fold(
                                Rational::from(0),
                                |sum, ((_, row), coefficient)| {
                                    sum + coefficient * &Rational::from(row[axis])
                                },
                            ) == target[axis]
                        }) {
                            continue;
                        }
                        let mut candidate = lmb.ext_atom(owner, GS.emr_mom, &indices, true);
                        for ((edge, _), coefficient) in selected.iter().zip(coefficients) {
                            if coefficient != 0 {
                                candidate += Atom::num(coefficient)
                                    * (GS.emr_mom(*edge, &indices[0])
                                        - lmb.ext_atom(*edge, GS.emr_mom, &indices, true));
                            }
                        }
                        candidates.insert(candidate);
                        // Different invertible minors of one basis give the
                        // same unique exact solution.
                        break;
                    }
                }
            }
            if candidates.is_empty() {
                return Err(eyre::eyre!(
                    "soft energy on crown edge {owner} has no exact routing through the active outer CFF edges {active_edges:?}"
                ));
            }
            alternatives.insert(
                component,
                candidates
                    .into_iter()
                    .map(|candidate| {
                        Ok((
                            candidate.clone(),
                            analyzer.analyze_reference_atom(&candidate)?,
                        ))
                    })
                    .collect::<Result<Vec<_>, EnergyPowerAnalysisError>>()?,
            );
        }

        // An unavailable crown carrier must still expose its power factors to
        // the planner; treating it as constant would silently lose the ranks
        // of the active edges onto which its energy has to be routed. Fixed
        // loads on these extra axes are identical in every assignment and do
        // not enter the final active-edge objective.
        let planned = EnergyPowerAnalyzer::for_physical_emr_edges(planning_edges).plan_view(
            numerator.as_view(),
            &mut 0,
            false,
            &BTreeSet::new(),
        )?;
        let mut proposals = planned.soft_momentum_assignments(&alternatives)?;
        proposals.sort_by_cached_key(|(degrees, assignments)| {
            let mut envelope = active_edges
                .iter()
                .map(|edge| {
                    degrees
                        .degrees
                        .get(&EnergyReference::Physical(*edge))
                        .copied()
                        .unwrap_or(0)
                })
                .collect::<Vec<_>>();
            envelope.sort_unstable_by(|left, right| right.cmp(left));
            (envelope, assignments.clone())
        });
        // The bounded search uses the certified basis/Pareto class only to
        // propose candidates. Actual native map count selects between them before
        // graph surface conversion; no global count optimum is claimed.
        proposals
            .into_iter()
            .unique_by(|(degrees, _)| {
                active_edges
                    .iter()
                    .map(|edge| {
                        degrees
                            .degrees
                            .get(&EnergyReference::Physical(*edge))
                            .copied()
                            .unwrap_or(0)
                    })
                    .collect::<Vec<_>>()
            })
            .take(ENERGY_ASSIGNMENT_PROPOSAL_BUDGET)
            .map(|(_, assignments)| {
                let optimized = planned.map(&mut |id, factor, _| {
                    Ok::<_, EnergyPowerAnalysisError>(assignments.get(&id).map_or_else(
                        || factor.clone(),
                        |choice| alternatives[factor][*choice].0.clone(),
                    ))
                })?;
                // Energy-independent powers/functions need no assignment traversal.
                // Erase their remaining soft metadata without changing their routing.
                Ok(optimized.replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() == 2
                        && let Some((_, UvMomentumProvenanceRole::DenominatorDerivedSoft, _)) =
                            GS.uv_momentum_provenance_data(momentum.get(0))
                    {
                        **output = GS.erase_uv_momentum_provenance(&view.to_owned());
                    }
                }))
            })
            .collect()
    }

    pub fn full_numerator_atom(&self) -> Atom {
        self.numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .expect("Graph numerator should be available")
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor
    }

    pub fn numerator_energy_power_caps(
        &self,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_numerator_energy_powers(EnergyPowerAnalysisMode::StrictPolynomial)
    }

    pub(crate) fn numerator_energy_power_upper_bounds(
        &self,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_numerator_energy_powers(EnergyPowerAnalysisMode::ConservativeUpperBound)
    }

    fn analyze_numerator_energy_powers(
        &self,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        // Keep edge-energy variables attached to their source denominator.
        // Rewriting them to a loop-momentum basis loses cancellation
        // information in non-coordinate UV directions, e.g. when q_i-q_j is
        // kept fixed while q_i and q_j are both large.
        let numerator = self.full_numerator_atom();
        self.analyze_numerator_atom_energy_powers(&numerator, mode)
    }

    fn analyze_numerator_atom_energy_powers(
        &self,
        numerator: &Atom,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        let internal_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy).then_some(edge)
            });

        EnergyPowerAnalyzer::for_physical_emr_edges(internal_edges)
            .analyze_physical_atom(numerator, mode)
    }

    pub fn automatic_numerator_energy_degree_bounds(
        &self,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        Ok(self
            .numerator_energy_power_upper_bounds()?
            .into_generation_bounds())
    }

    pub(crate) fn automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
        &self,
        numerators: impl IntoIterator<Item = impl Borrow<Atom>>,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        min_degree: usize,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(active_edges);
        // Branches are evaluated independently. Take the maximum of their
        // ranks without constructing a sum in which leading powers could cancel.
        let mut bounds = EnergyPowerCapMap::default();
        for numerator in numerators {
            bounds.max_assign(analyzer.analyze_atom(numerator.borrow())?);
        }
        Ok(bounds
            .iter()
            .filter_map(|(edge, degree)| (degree >= min_degree).then_some((edge.into(), degree)))
            .collect())
    }

    #[cfg(test)]
    pub(crate) fn plan_numerator_energy_assignment_in_atom_excluding(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        Ok(self
            .plan_numerator_energy_assignment_proposals_in_atom_excluding(
                numerator,
                excluded_edges,
                candidates,
            )?
            .remove(0))
    }

    pub(crate) fn plan_numerator_energy_assignment_proposals_in_atom_excluding(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<Vec<EnergyPowerAssignmentPlan>, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        let mut analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(active_edges);
        // Other components remain symbolic coefficients of this contour.
        analyzer.active_uv_classes = Some(
            candidates
                .by_family
                .keys()
                .filter_map(|family| match family {
                    EnergyCandidateFamily::UvClass(class) => Some(*class),
                    _ => None,
                })
                .collect(),
        );
        analyzer.plan_atom_assignment_proposals(numerator, candidates)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::{
        dot,
        graph::{Graph, LMBext, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::energy_degree::{
            EnergyCandidateFamily, EnergyPowerAnalysisError, EnergyPowerAnalyzer, EnergyReference,
            EquivalentEnergyCandidates,
        },
        utils::{
            GS, W_,
            symbols::{UvDenominatorClassId, UvMomentumProvenanceRole},
        },
    };
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::shadowing::symbolica_utils::ReplaceBuilderExt;
    use spenso::structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, Minkowski},
    };
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        function, symbol,
    };

    fn mink_index(label: &str) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::var(symbol!(label)).as_view())
            .finish()
    }

    fn mink_component(index: i64) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::num(index).as_view())
            .finish()
    }

    fn compact_minkowski_vector() -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .finish()
    }

    fn bounds(expression: Atom) -> Vec<(usize, usize)> {
        EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap()
            .into_generation_bounds()
    }

    fn bounds_with_internal_edges(expression: Atom) -> Vec<(usize, usize)> {
        EnergyPowerAnalyzer::with_internal_edges([EdgeIndex(7), EdgeIndex(9)], [EdgeIndex(3)])
            .analyze_atom(&expression)
            .unwrap()
            .into_generation_bounds()
    }

    #[test]
    fn emr_minkowski_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_index("mu"));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_numeric_minkowski_label_counts_as_abstract_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_component(1));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_numeric_zero_minkowski_label_counts_as_abstract_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_component(0));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn compact_emr_vector_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, compact_minkowski_vector());
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn dot_product_composes_compact_vector_energy_powers() {
        let q3 = function!(GS.emr_mom, 3, compact_minkowski_vector());
        let q7 = function!(GS.emr_mom, 7, compact_minkowski_vector());
        assert_eq!(bounds(function!(GS.dot, q3.clone(), q3)), vec![(3, 2)]);
        assert_eq!(bounds(function!(GS.dot, q7, q3)), vec![(3, 1), (7, 1)]);
    }

    #[test]
    fn positive_denominator_wrapper_analyzes_only_its_polynomial_value() {
        let q3 = function!(GS.emr_mom, 3, compact_minkowski_vector());
        let full_expr = function!(GS.dot, q3.clone(), q3.clone()) - Atom::one();
        let denominator = GS.den(3, q3, Atom::one(), full_expr);
        assert_eq!(bounds(denominator), vec![(3, 2)]);
    }

    #[test]
    fn emr_external_edge_minkowski_index_does_not_count_when_internal_edges_are_known() {
        let expression = function!(GS.emr_mom, 0, mink_index("mu"));
        assert!(bounds_with_internal_edges(expression).is_empty());
    }

    #[test]
    fn emr_temporal_concrete_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.call(0));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_spatial_concrete_index_does_not_count() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.call(2));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn emr_direct_temporal_component_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, 0);
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_direct_spatial_component_does_not_count() {
        let expression = function!(GS.emr_mom, 3, 1);
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn spatial_emr_vector_does_not_count() {
        let expression = function!(GS.emr_vec, 3, mink_index("mu"));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn lmb_energy_maps_to_carrier_edge() {
        let expression = function!(GS.loop_mom, 1, mink_index("mu"));
        assert_eq!(bounds(expression), vec![(9, 1)]);
    }

    #[test]
    fn physical_emr_analysis_rejects_loop_basis_energy() {
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(7)]);
        let loop_momentum = function!(GS.loop_mom, 0, mink_index("mu"));
        assert!(matches!(
            analyzer.analyze_atom(&loop_momentum),
            Err(EnergyPowerAnalysisError::LoopMomentumInPhysicalEmrAnalysis { .. })
        ));

        let emr_momentum = function!(GS.emr_mom, 7, mink_index("mu"));
        assert_eq!(
            analyzer
                .analyze_atom(&emr_momentum)
                .unwrap()
                .into_generation_bounds(),
            vec![(7, 1)]
        );
    }

    #[test]
    fn sums_take_edgewise_max_and_products_add() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let k1 = function!(GS.loop_mom, 1, mink_index("nu"));
        let expression = q3.clone() * k1.clone() + q3.pow(Atom::num(3));
        assert_eq!(bounds(expression), vec![(3, 3), (9, 1)]);
    }

    #[test]
    fn powers_of_affine_atoms_produce_edgewise_bounds() {
        let k0 = function!(GS.loop_mom, 0, mink_index("mu"));
        let k1 = function!(GS.loop_mom, 1, mink_index("nu"));
        let expression = (k0 + k1).pow(Atom::num(2));
        assert_eq!(bounds(expression), vec![(7, 2), (9, 2)]);
    }

    #[test]
    fn factorized_products_compose_energy_bounds_without_expansion() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression = (q3.clone() + Atom::num(1)) * (q3 + Atom::num(2));
        assert_eq!(bounds(expression), vec![(3, 2)]);
    }

    fn two_equivalent_candidates() -> EquivalentEnergyCandidates {
        EquivalentEnergyCandidates::try_from_source_occurrences([(EdgeIndex(3), vec![10, 11])])
            .unwrap()
    }

    fn uv_owned_q3(denominator_derived: bool, index: Atom) -> Atom {
        let hard_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_momentum_provenance_tag(
                Atom::num(3).as_view(),
                denominator_derived,
                hard_momentum.as_view(),
            ))
            .add_arg(index.as_view())
            .finish()
    }

    fn production_candidates_with_two_derived_copies() -> EquivalentEnergyCandidates {
        EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            EdgeIndex(3),
            10,
            vec![11, 12],
        )])
        .unwrap()
    }

    fn planned_candidates(expression: &Atom, candidates: &EquivalentEnergyCandidates) {
        let edge = EdgeIndex(3);
        let plans = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment_proposals(expression, candidates)
            .unwrap();
        assert!(!plans.is_empty());
        for plan in plans {
            let mapped = plan.map_factors(|factor, assignments| {
                let Some(candidate) = assignments.get(&EnergyReference::Physical(edge)) else {
                    return Ok::<_, ()>(factor.clone());
                };
                let _ = factor.replace_map(|view, _, _| {
                    let AtomView::Fun(momentum) = view else { return; };
                    if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 2 {
                        return;
                    }
                    let family = match GS.uv_momentum_provenance_data(momentum.get(0)) {
                        Some((owner, UvMomentumProvenanceRole::DenominatorDerived, _)) =>
                            EnergyCandidateFamily::DenominatorDerived(owner),
                        Some((owner, _, _)) => EnergyCandidateFamily::Fixed(owner),
                        None => EnergyCandidateFamily::Unprovenanced(edge),
                    };
                    assert!(candidates.get(family).unwrap().contains(candidate),
                        "every original or derivative-owned factor must stay in its certified occurrence domain");
                });
                Ok(GS.erase_uv_momentum_provenance(factor)
                    .replace(function!(GS.emr_mom, 3, W_.x___))
                    .with(function!(GS.emr_mom, *candidate, W_.x___)))
            }).unwrap();
            let bounds = plan
                .energy_degree_bounds()
                .iter()
                .copied()
                .collect::<BTreeMap<_, _>>();
            let actual = EnergyPowerAnalyzer::new([])
                .analyze_atom(&mapped)
                .unwrap()
                .into_generation_bounds();
            assert!(
                actual
                    .iter()
                    .all(|(edge, degree)| bounds.get(edge).is_some_and(|bound| bound >= degree)),
                "the returned capacity must certify the fully mapped numerator"
            );
            let mut reconstructed = mapped;
            for candidate in bounds.keys() {
                reconstructed = reconstructed
                    .replace(function!(GS.emr_mom, *candidate, W_.x___))
                    .with(function!(GS.emr_mom, 3, W_.x___));
            }
            assert!(
                (reconstructed - GS.erase_uv_momentum_provenance(expression))
                    .expand()
                    .is_zero(),
                "every proposed assignment must reconstruct the complete numerator value: {expression}"
            );
        }
    }

    #[test]
    fn assignment_plans_preserve_owner_domains_and_complete_factorized_values() {
        let fixed = uv_owned_q3(false, GS.cind(0));
        let derived = uv_owned_q3(true, GS.cind(0));
        let plain = function!(GS.emr_mom, 3, GS.cind(0));
        let vector = uv_owned_q3(true, compact_minkowski_vector());
        let a = Atom::var(symbol!("energy_assignment_test::a"));
        let b = Atom::var(symbol!("energy_assignment_test::b"));
        let linear = symbol!("energy_assignment_test::linear"; Linear);
        let denominator = function!(GS.den, 3, function!(GS.emr_mom, 3), 1, derived.pow(2));
        // The quadratic wrapper stays indivisible even below a fixed quartic.
        // The complete certificate permits any valid distribution of the
        // derivative-created factors, including additive capacity reuse.
        let expressions = [
            (&plain + Atom::one()) * (&plain + Atom::num(2)),
            function!(GS.dot, &vector, &vector),
            function!(linear, &derived) + function!(linear, &derived * &a),
            (&derived + &a).pow(2),
            denominator.pow(2) + derived.pow(2),
            fixed.pow(4) * &denominator * derived.pow(2),
            fixed.pow(4) * derived.pow(2),
            (&derived + Atom::one()) * (&derived + Atom::num(2)) * (&derived + Atom::num(3)),
            &fixed * &derived,
            &fixed * (&a * &derived + &b * &derived),
            &a * &fixed * &derived + &b * derived.pow(2),
            derived.pow(3) + derived.pow(2),
            derived.pow(2),
        ];
        for candidates in [
            two_equivalent_candidates(),
            production_candidates_with_two_derived_copies(),
            EquivalentEnergyCandidates::try_from_source_occurrences([(
                EdgeIndex(3),
                vec![10, 11, 12],
            )])
            .unwrap(),
        ] {
            for expression in &expressions {
                planned_candidates(expression, &candidates);
            }
        }
    }

    #[test]
    fn soft_taylor_routing_balances_the_outer_envelope_without_moving_originals()
    -> color_eyre::Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph soft_taylor_cograph {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2]
            b -> c [id=3 lmb_id=1]
            c -> d [id=4]
            d -> a [id=5]
            d -> outgoing [id=6]
        })?;
        let soft_component = |index: Atom| {
            let owner = Atom::num(4);
            let payload = function!(GS.emr_mom, 4);
            function!(
                GS.emr_mom,
                GS.uv_momentum_provenance_tag(
                    owner.as_view(),
                    UvMomentumProvenanceRole::DenominatorDerivedSoft,
                    payload.as_view(),
                ),
                index
            )
        };
        let soft = |index| soft_component(GS.cind(index));
        let fixed_quartic = GS.emr_mom(EdgeIndex(3), GS.cind(0)).pow(4);
        let fixed_linear = GS.emr_mom(EdgeIndex(4), GS.cind(0));
        let numerator = &fixed_quartic * &fixed_linear * soft(0).pow(3);
        let active = [EdgeIndex(3), EdgeIndex(4), EdgeIndex(5)];
        let optimized = graph.optimize_soft_momentum_routing(&numerator, active)?;
        assert_ne!(
            optimized.replace(fixed_quartic.clone()).with(Atom::one()),
            optimized,
            "the original quartic remains on its original edge",
        );
        let neutral = graph.normal_emr_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[W_.x___],
            |_| true,
        );
        assert!(
            (optimized.replace_multiple(&neutral).collect_factors()
                - GS.erase_uv_momentum_provenance(&numerator)
                    .replace_multiple(&neutral)
                    .collect_factors())
            .collect_factors()
            .is_zero(),
            "soft routing must be an exact off-shell identity"
        );
        assert!(!optimized.contains_symbol(GS.uv_momentum_provenance));
        assert_eq!(
            graph.optimize_soft_momentum_routing(&soft(1).pow(2), active)?,
            GS.emr_mom(EdgeIndex(4), GS.cind(1)).pow(2),
            "spatial powers lose metadata without being expanded or rerouted",
        );
        let surviving = [EdgeIndex(3), EdgeIndex(5)];
        let rerouted = graph.optimize_soft_momentum_routing(&soft(0).pow(2), surviving)?;
        assert!(
            !rerouted
                .replace(GS.emr_mom(EdgeIndex(4), GS.cind(0)))
                .matches(),
            "an unavailable crown carrier must be routed through active edges"
        );
        assert!(
            (rerouted.replace_multiple(&neutral).collect_factors()
                - GS.erase_uv_momentum_provenance(&soft(0).pow(2))
                    .replace_multiple(&neutral)
                    .collect_factors())
            .collect_factors()
            .is_zero()
        );
        assert!(
            graph
                .optimize_soft_momentum_routing(&soft(0), [EdgeIndex(1)])
                .unwrap_err()
                .to_string()
                .contains("no exact routing")
        );
        let vector_index = compact_minkowski_vector();
        let dot = function!(
            GS.dot,
            soft_component(vector_index.clone()),
            GS.emr_mom(EdgeIndex(0), vector_index)
        );
        // This diagnostic routing polynomial exercises abstract vector slots
        // and a factorized additive power without using numerator expansion.
        let factorized = &fixed_quartic * (dot + Atom::one()).pow(2);
        let routed = graph.optimize_soft_momentum_routing(&factorized, active)?;
        assert!(
            (routed.replace_multiple(&neutral).collect_factors()
                - GS.erase_uv_momentum_provenance(&factorized)
                    .replace_multiple(&neutral)
                    .collect_factors())
            .collect_factors()
            .is_zero(),
            "abstract-slot routing must obey the same exact linear identity"
        );
        // A dependent external carrier has two literal exact representations,
        // but both request the same all-zero active CFF capacity. A lone factor
        // bypasses composite-node collection; every retained proposal must still
        // preserve its complete external value.
        let external_owner = [EdgeIndex(0), EdgeIndex(6)]
            .into_iter()
            .find(|edge| {
                graph
                    .loop_momentum_basis
                    .ext_atom(*edge, GS.emr_mom, &[GS.cind(0)], true)
                    != GS.emr_mom(*edge, GS.cind(0))
            })
            .expect("one of the two external momenta is dependent");
        let owner = Atom::num(usize::from(external_owner));
        let payload = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(external_owner))
            .finish();
        let external_soft = function!(
            GS.emr_mom,
            GS.uv_momentum_provenance_tag(
                owner.as_view(),
                UvMomentumProvenanceRole::DenominatorDerivedSoft,
                payload.as_view(),
            ),
            GS.cind(0)
        );
        let literal = GS.erase_uv_momentum_provenance(&external_soft);
        let external =
            graph
                .loop_momentum_basis
                .ext_atom(external_owner, GS.emr_mom, &[GS.cind(0)], true);
        assert_ne!(
            literal, external,
            "the fixture must expose duplicate capacity through different atoms"
        );
        for source in [external_soft.clone(), (external_soft + Atom::one()).pow(8)] {
            let proposals = graph.soft_momentum_routing_proposals(&source, active)?;
            assert!(!proposals.is_empty());
            for proposal in proposals {
                assert!(
                    EnergyPowerAnalyzer::for_physical_emr_edges(active)
                        .analyze_atom(&proposal)?
                        .is_empty()
                );
                assert!(
                    (proposal.replace_multiple(&neutral).collect_factors()
                        - GS.erase_uv_momentum_provenance(&source)
                            .replace_multiple(&neutral)
                            .collect_factors())
                    .is_zero(),
                    "all zero-capacity routes must retain the complete factorized external coefficient"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn positive_denominator_completion_stays_in_one_occurrence_frame() {
        let q3 = uv_owned_q3(true, GS.cind(0));
        let spatial_mass = Atom::var(symbol!("energy_assignment_test::spatial_mass"));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg((q3.pow(2) - &spatial_mass).as_view())
            .finish();
        planned_candidates(&denominator, &two_equivalent_candidates());
    }

    #[test]
    fn positive_denominator_completion_obeys_embedded_emr_provenance() {
        let edge = EdgeIndex(3);
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = |energy: &Atom| {
            FunctionBuilder::new(GS.den)
                .add_arg(3)
                .add_arg(source_momentum.as_view())
                .add_arg(Atom::one().as_view())
                .add_arg(energy.clone().pow(2).as_view())
                .finish()
        };
        let assert_assignment = |energy: Atom, candidates: EquivalentEnergyCandidates| {
            planned_candidates(&denominator(&energy), &candidates);
        };

        let base_only = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            edge,
            10,
            Vec::new(),
        )])
        .unwrap();
        assert_assignment(function!(GS.emr_mom, 3, GS.cind(0)), base_only);

        assert_assignment(
            uv_owned_q3(false, GS.cind(0)),
            production_candidates_with_two_derived_copies(),
        );
        assert_assignment(
            uv_owned_q3(true, GS.cind(0)),
            production_candidates_with_two_derived_copies(),
        );
    }

    #[test]
    fn positive_denominator_completion_rejects_mixed_emr_provenance() {
        let edge = EdgeIndex(3);
        let fixed = uv_owned_q3(false, GS.cind(0));
        let derived = uv_owned_q3(true, GS.cind(0));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg((fixed * derived).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(
                &denominator,
                &production_candidates_with_two_derived_copies(),
            )
            .unwrap_err();

        assert!(matches!(
            error,
            EnergyPowerAnalysisError::MixedPositiveDenominatorProvenance { owner: 3, .. }
        ));
    }

    #[test]
    fn separate_positive_denominator_completions_may_use_separate_occurrences() {
        let q3 = uv_owned_q3(true, GS.cind(0));
        let source_momentum = FunctionBuilder::new(GS.emr_mom).add_arg(3).finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(source_momentum.as_view())
            .add_arg(Atom::one().as_view())
            .add_arg(q3.pow(2).as_view())
            .finish();
        let expression = denominator.pow(2);
        planned_candidates(&expression, &two_equivalent_candidates());
    }

    #[test]
    fn positive_denominator_completion_rejects_a_nonlocal_energy_owner() {
        let q4 = FunctionBuilder::new(GS.emr_mom)
            .add_arg(
                GS.uv_momentum_provenance_tag(
                    Atom::num(4).as_view(),
                    true,
                    FunctionBuilder::new(GS.emr_mom)
                        .add_arg(4)
                        .finish()
                        .as_view(),
                ),
            )
            .add_arg(GS.cind(0))
            .finish();
        let denominator = FunctionBuilder::new(GS.den)
            .add_arg(3)
            .add_arg(
                FunctionBuilder::new(GS.emr_mom)
                    .add_arg(3)
                    .finish()
                    .as_view(),
            )
            .add_arg(Atom::one().as_view())
            .add_arg(q4.pow(2).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(3), EdgeIndex(4)])
            .plan_atom_assignment(&denominator, &two_equivalent_candidates())
            .unwrap_err();

        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonlocalPositiveDenominator {
                owner: 3,
                dependent: 4,
                ..
            }
        ));
    }

    #[test]
    fn source_energy_candidates_fix_originals_and_offer_all_derived_occurrences() {
        let candidates = EquivalentEnergyCandidates::try_from_source_occurrences([
            (EdgeIndex(3), vec![12, 10, 11]),
            (EdgeIndex(7), vec![22, 21, 20]),
            (EdgeIndex(9), vec![33, 31, 32, 30]),
        ])
        .unwrap();

        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(3))),
            Some([10, 11, 12].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(7))),
            Some([20, 21, 22].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(9))),
            Some([30, 31, 32, 33].as_slice())
        );
        for (edge, canonical) in [(EdgeIndex(3), 10), (EdgeIndex(7), 20), (EdgeIndex(9), 30)] {
            assert_eq!(
                candidates.get(EnergyCandidateFamily::Unprovenanced(edge)),
                Some([canonical].as_slice())
            );
            assert_eq!(
                candidates.get(EnergyCandidateFamily::Fixed(edge)),
                Some([canonical].as_slice())
            );
        }
    }

    #[test]
    fn source_energy_candidates_keep_occurrences_disjoint_between_owners() {
        let overlapping = EquivalentEnergyCandidates::try_from_source_occurrences([
            (EdgeIndex(3), vec![10, 11, 12]),
            (EdgeIndex(7), vec![12]),
        ])
        .unwrap_err();
        assert!(matches!(
            overlapping,
            EnergyPowerAnalysisError::OverlappingEquivalentEnergyCandidates { .. }
        ));

        let duplicate = EquivalentEnergyCandidates::try_from_source_occurrences([(
            EdgeIndex(3),
            vec![10, 11, 11],
        )])
        .unwrap_err();
        assert!(matches!(
            duplicate,
            EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate { .. }
        ));
    }

    #[test]
    fn partitioned_source_candidates_do_not_infer_the_base_from_energy_id_order() {
        let candidates = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([
            (EdgeIndex(3), 12, vec![10, 11]),
            (EdgeIndex(7), 20, Vec::new()),
        ])
        .unwrap();

        assert_eq!(
            candidates.get(EnergyCandidateFamily::Unprovenanced(EdgeIndex(3))),
            Some([12].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::Fixed(EdgeIndex(3))),
            Some([12].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(3))),
            Some([10, 11].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::Fixed(EdgeIndex(7))),
            Some([20].as_slice())
        );
        assert_eq!(
            candidates.get(EnergyCandidateFamily::DenominatorDerived(EdgeIndex(7))),
            None,
            "an owner without a derivative-created serial copy must not accept a derivative-owned numerator factor"
        );
    }

    #[test]
    fn partitioned_source_candidates_reject_derived_factors_without_a_copy() {
        let edge = EdgeIndex(3);
        let candidates = EquivalentEnergyCandidates::try_from_partitioned_source_occurrences([(
            edge,
            12,
            Vec::new(),
        )])
        .unwrap();
        let derived = uv_owned_q3(true, mink_index("mu"));
        assert!(matches!(
            EnergyPowerAnalyzer::for_physical_emr_edges([edge])
                .plan_atom_assignment(&derived, &candidates),
            Err(EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates { .. })
        ));
    }

    #[test]
    fn assignment_certificate_rejects_understated_capacity() {
        let class = UvDenominatorClassId(0);
        let q = function!(GS.emr_mom, GS.uv_class_ref(class), GS.cind(0));
        let mut candidates = EquivalentEnergyCandidates::try_from_source_occurrences([]).unwrap();
        candidates.add_uv_classes([(class, vec![10, 11])]).unwrap();
        let mut plan = EnergyPowerAnalyzer::for_physical_emr_edges([])
            .plan_atom_assignment(&q.pow(3), &candidates)
            .unwrap();
        plan.certify_bounds().unwrap();
        plan.energy_degree_bounds[0].1 -= 1;
        assert!(matches!(
            plan.certify_bounds(),
            Err(EnergyPowerAnalysisError::InvalidAssignmentCertificate { .. })
        ));
    }

    #[test]
    fn canonical_class_cyclic_lifts_preserve_factorization_and_exact_bounds() {
        let class = UvDenominatorClassId(0);
        let reference = EnergyReference::UvClass(class);
        let q = function!(GS.emr_mom, GS.uv_class_ref(class), GS.cind(0));
        let expressions = [
            q.pow(5),
            (q.clone() + Atom::one()).pow(5),
            (q.pow(2) + Atom::one()) * (q.pow(3) + Atom::one()),
            q.pow(5) + q.pow(3),
            function!(GS.dot, &q, &q) * q.pow(3),
        ];
        let mut candidates = EquivalentEnergyCandidates::try_from_source_occurrences([]).unwrap();
        candidates
            .add_uv_classes([(class, (10..15).collect())])
            .unwrap();
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([]);
        for expression in expressions {
            let plans = analyzer
                .plan_atom_assignment_proposals(&expression, &candidates)
                .unwrap();
            assert_eq!(plans.len(), 1);
            let plan = &plans[0];
            assert_eq!(
                plan.energy_degree_bounds(),
                &[(10, 1), (11, 1), (12, 1), (13, 1), (14, 1)]
            );
            let mut prepared_leaves = 0;
            let prepared = plan
                .prepare_factors(|factor, assignments| {
                    prepared_leaves += 1;
                    Ok::<_, ()>(assignments.get(&reference).map_or_else(
                        || factor.clone(),
                        |occurrence| {
                            factor
                                .replace(q.clone())
                                .with(GS.emr_mom(EdgeIndex(*occurrence), GS.cind(0)))
                        },
                    ))
                })
                .unwrap();
            assert!(prepared_leaves > 0);
            let mapped = prepared
                .map(&mut |_, factor, _| Ok::<_, ()>(factor.clone()))
                .unwrap();
            let actual = EnergyPowerAnalyzer::new([])
                .analyze_atom(&mapped)
                .unwrap()
                .into_generation_bounds();
            assert_eq!(actual, plan.energy_degree_bounds());
            let collapsed = mapped.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() == 2
                    && usize::try_from(momentum.get(0)).is_ok_and(|id| (10..15).contains(&id))
                {
                    **output = q.clone();
                }
            });
            assert_eq!(collapsed.collect_factors(), expression.collect_factors());
        }
    }

    #[test]
    fn canonical_large_powers_keep_assignment_cycles_compressed() {
        let class = UvDenominatorClassId(2);
        let q = function!(GS.emr_mom, GS.uv_class_ref(class), GS.cind(0));
        let mut candidates = EquivalentEnergyCandidates::try_from_source_occurrences([]).unwrap();
        candidates
            .add_uv_classes([(class, (0..5).collect())])
            .unwrap();
        let plans = EnergyPowerAnalyzer::for_physical_emr_edges([])
            .plan_atom_assignment_proposals(&q.pow(500_000), &candidates)
            .unwrap();
        assert_eq!(
            plans[0].energy_degree_bounds(),
            &[
                (0, 100_000),
                (1, 100_000),
                (2, 100_000),
                (3, 100_000),
                (4, 100_000)
            ]
        );
        assert_eq!(
            plans[1].energy_degree_bounds(),
            &[(0, 499_996), (1, 1), (2, 1), (3, 1), (4, 1)]
        );
        for plan in plans {
            let mut leaves = 0;
            plan.prepare_factors(|factor, _| {
                leaves += 1;
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
            assert!(
                leaves <= 12,
                "cycles must not instantiate half a million factor assignments"
            );
        }
    }

    #[test]
    fn canonical_challengers_preserve_pools_and_opaque_blocks() {
        let class = UvDenominatorClassId(0);
        let q = function!(GS.emr_mom, GS.uv_class_ref(class), GS.cind(0));
        let mut candidates = EquivalentEnergyCandidates::try_from_source_occurrences([]).unwrap();
        candidates
            .add_uv_classes([(class, (10..15).collect())])
            .unwrap();
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([]);
        let plans = analyzer
            .plan_atom_assignment_proposals(&q.pow(7), &candidates)
            .unwrap();
        assert_eq!(plans.len(), 2);
        assert_eq!(
            plans[0].energy_degree_bounds(),
            &[(10, 2), (11, 2), (12, 1), (13, 1), (14, 1)]
        );
        assert_eq!(
            plans[1].energy_degree_bounds(),
            &[(10, 3), (11, 1), (12, 1), (13, 1), (14, 1)]
        );
        let seen = plans
            .iter()
            .map(|plan| plan.energy_degree_bounds().to_vec())
            .collect::<Vec<_>>();
        let third = plans[1].placement_challenger(&seen).unwrap().unwrap();
        assert_eq!(
            third.energy_degree_bounds(),
            &[(10, 1), (11, 1), (12, 1), (13, 1), (14, 3)]
        );

        let den = function!(
            GS.den,
            GS.uv_class_ref(class),
            &q,
            1,
            q.pow(2) - Atom::one()
        );
        let opaque = analyzer
            .plan_atom_assignment_proposals(&den.pow(2), &candidates)
            .unwrap();
        assert_eq!(opaque[0].energy_degree_bounds(), &[(10, 2), (11, 2)]);
        let _ = opaque[0]
            .map_factors(|factor, assignments| {
                assert_eq!(factor, &den);
                assert_eq!(assignments.len(), 1);
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
    }

    #[test]
    fn canonical_pool_validation_preserves_physical_restrictions_and_domains() {
        let class = UvDenominatorClassId(0);
        let other = UvDenominatorClassId(1);
        let mut candidates = production_candidates_with_two_derived_copies();
        candidates
            .add_uv_classes([(class, vec![10, 11, 12])])
            .unwrap();
        assert_eq!(
            candidates.get(EnergyCandidateFamily::Fixed(EdgeIndex(3))),
            Some([10].as_slice())
        );
        assert!(candidates.add_uv_classes([(other, vec![12, 13])]).is_err());
        let mut analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([]);
        analyzer.active_uv_classes = Some([class].into_iter().collect());
        let active = function!(GS.emr_mom, GS.uv_class_ref(class), GS.cind(0));
        let inactive = function!(GS.emr_mom, GS.uv_class_ref(other), GS.cind(0));
        let expression = active.pow(3) * inactive.pow(7);
        let plans = analyzer
            .plan_atom_assignment_proposals(&expression, &candidates)
            .unwrap();
        assert_eq!(
            plans[0].energy_degree_bounds(),
            &[(10, 1), (11, 1), (12, 1)]
        );
        let mapped = plans[0]
            .map_factors(|factor, _| Ok::<_, ()>(factor.clone()))
            .unwrap();
        assert_eq!(mapped, expression);
    }

    #[test]
    fn assignment_rejects_a_missing_certified_candidate_set() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let candidates =
            EquivalentEnergyCandidates::try_from_source_occurrences([(EdgeIndex(7), vec![12])])
                .unwrap();
        let error = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(3)])
            .plan_atom_assignment(&q3, &candidates)
            .unwrap_err();

        assert!(matches!(
            error,
            EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates { edge: 3, degree: 1 }
        ));
    }

    #[test]
    fn overflowing_polynomial_degrees_return_an_error() {
        let energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
        let large_degree = usize::MAX / 3 + 1;
        let exponent = Atom::num(large_degree as u64);
        let powered = energy.clone().pow(exponent.clone());
        let overflow_product = (1..=3).fold(Atom::one(), |product, shift| {
            product * (energy.clone() + Atom::num(shift)).pow(exponent.clone())
        });
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges([EdgeIndex(0)]);
        for expression in [(powered.clone() + Atom::one()).pow(3), overflow_product] {
            assert!(matches!(
                analyzer.analyze_atom(&expression),
                Err(EnergyPowerAnalysisError::EnergyDegreeOverflow)
            ));
            assert!(matches!(
                analyzer.analyze_atom_upper_bound(&expression),
                Err(EnergyPowerAnalysisError::EnergyDegreeOverflow)
            ));
        }
        assert_eq!(
            analyzer
                .analyze_atom(&(powered + Atom::one()).pow(2))
                .unwrap()
                .into_generation_bounds(),
            vec![(0, 2 * large_degree)]
        );
    }

    #[test]
    fn strict_analysis_rejects_negative_energy_powers() {
        let expression = function!(GS.emr_mom, 3, mink_index("mu")).pow(Atom::num(-1));
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonPolynomialEnergyPower { .. }
        ));
    }

    #[test]
    fn strict_analysis_ignores_energy_powers_outside_the_active_edge_set() {
        let analyzer =
            EnergyPowerAnalyzer::with_internal_edges([EdgeIndex(7), EdgeIndex(9)], [EdgeIndex(5)]);
        for expression in [
            function!(GS.emr_mom, 3, mink_index("mu")).pow(Atom::num(-1)),
            function!(GS.loop_mom, 0, mink_index("mu")).pow(Atom::num(-1)),
        ] {
            assert!(analyzer.analyze_atom(&expression).unwrap().is_empty());
        }
    }

    #[test]
    fn energy_dependent_exponents_are_rejected_even_for_constant_bases() {
        let exponent = function!(GS.emr_mom, 3, 0);
        let expression = Atom::num(2).pow(exponent);
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonPolynomialEnergyPower { .. }
        ));
    }

    #[test]
    fn strict_analysis_rejects_opaque_energy_functions() {
        let expression = FunctionBuilder::new(symbol!("opaque_energy_function"))
            .add_arg(function!(GS.emr_mom, 3, mink_index("mu")).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::OpaqueEnergyFunction { .. }
        ));
    }

    #[test]
    fn zero_arguments_do_not_hide_energy_in_opaque_functions() {
        let expression = FunctionBuilder::new(symbol!("opaque_energy_function_with_zero"))
            .add_arg(Atom::zero().as_view())
            .add_arg(function!(GS.emr_mom, 3, 0).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom_upper_bound(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::OpaqueEnergyFunction { .. }
        ));
    }

    #[test]
    fn declared_multilinear_functions_preserve_energy_bounds() {
        let expression = FunctionBuilder::new(symbol!("multilinear_energy_function"; Linear))
            .add_arg(function!(GS.emr_mom, 3, mink_index("mu")).as_view())
            .finish();
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn malformed_emr_edge_argument_is_reported() {
        let expression = FunctionBuilder::new(GS.emr_mom)
            .add_arg(Atom::var(symbol!("edge")).as_view())
            .add_arg(mink_index("mu").as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::InvalidEmrEdgeArgument { .. }
        ));
    }

    #[test]
    fn malformed_loop_argument_is_reported() {
        let expression = FunctionBuilder::new(GS.loop_mom)
            .add_arg(Atom::var(symbol!("loop_id")).as_view())
            .add_arg(mink_index("mu").as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::InvalidLoopMomentumArgument { .. }
        ));
    }
}
