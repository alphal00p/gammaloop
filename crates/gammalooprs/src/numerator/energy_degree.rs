use std::{
    collections::{BTreeMap, BTreeSet},
    ops::MulAssign,
};

use linnet::half_edge::involution::EdgeIndex;
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski},
};
use symbolica::atom::{Atom, AtomView, FunctionBuilder, Symbol};
use thiserror::Error;

use crate::{graph::Graph, utils::GS, uv::UltravioletGraph};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct EnergyPowerCapMap {
    degrees: BTreeMap<EdgeIndex, usize>,
}

impl EnergyPowerCapMap {
    pub fn is_empty(&self) -> bool {
        self.degrees.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (EdgeIndex, usize)> + '_ {
        self.degrees.iter().map(|(edge, degree)| (*edge, *degree))
    }

    pub fn into_generation_bounds(self) -> Vec<(usize, usize)> {
        self.degrees
            .into_iter()
            .map(|(edge, degree)| (usize::from(edge), degree))
            .collect()
    }

    fn unit(edge: EdgeIndex) -> Self {
        let mut degrees = BTreeMap::new();
        degrees.insert(edge, 1);
        Self { degrees }
    }

    fn add_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            *self.degrees.entry(edge).or_insert(0) += degree;
        }
    }

    fn max_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = (*slot).max(degree);
        }
    }

    fn scale(&mut self, exponent: usize) {
        for degree in self.degrees.values_mut() {
            degree.mul_assign(exponent);
        }
    }
}

/// Exact energy variables which have been certified to represent the same
/// algebraic energy for one physical EMR edge.
///
/// For each physical edge, candidates are restricted to the largest class
/// with one algebraic on-shell energy; the lowest occurrence fixes equal-size
/// ties. Momentum-routing signs are not part of that energy identity and are
/// restored while mapping the numerator. Candidate sets are also kept disjoint. This makes the minimax
/// assignment below separable by physical edge and avoids introducing a
/// general flow solver for independent load-balancing problems.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EquivalentEnergyCandidates {
    by_physical_edge: BTreeMap<EdgeIndex, Vec<usize>>,
}

impl EquivalentEnergyCandidates {
    fn try_new<K: Eq>(
        groups: impl IntoIterator<Item = (EdgeIndex, Vec<(usize, K)>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        let mut by_physical_edge = BTreeMap::new();
        let mut candidate_owners = BTreeMap::<usize, EdgeIndex>::new();

        for (edge, mut candidates) in groups {
            if by_physical_edge.contains_key(&edge) {
                return Err(
                    EnergyPowerAnalysisError::DuplicatePhysicalEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            candidates.sort_by_key(|(candidate, _)| *candidate);
            if candidates.is_empty() {
                return Err(
                    EnergyPowerAnalysisError::EmptyEquivalentEnergyCandidateSet {
                        edge: edge.into(),
                    },
                );
            }
            for pair in candidates.windows(2) {
                if pair[0].0 == pair[1].0 {
                    return Err(
                        EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate {
                            edge: edge.into(),
                            candidate: pair[0].0,
                        },
                    );
                }
            }
            for (candidate, _) in &candidates {
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

            let mut equivalent_classes = Vec::<Vec<(usize, K)>>::new();
            for candidate in candidates {
                if let Some(class) = equivalent_classes
                    .iter_mut()
                    .find(|class| class[0].1 == candidate.1)
                {
                    class.push(candidate);
                } else {
                    equivalent_classes.push(vec![candidate]);
                }
            }
            // More equivalent occurrences minimize the maximal assigned
            // degree. The lowest canonical occurrence fixes equal-size ties.
            equivalent_classes.sort_by(|left, right| {
                right
                    .len()
                    .cmp(&left.len())
                    .then_with(|| left[0].0.cmp(&right[0].0))
            });
            let candidate_ids = equivalent_classes
                .remove(0)
                .into_iter()
                .map(|(candidate, _)| candidate)
                .collect();
            by_physical_edge.insert(edge, candidate_ids);
        }

        Ok(Self { by_physical_edge })
    }

    /// Build candidate sets from the exact algebraic on-shell energies of
    /// denominator occurrences. This is the production certification boundary;
    /// arbitrary equality keys are accepted only by the focused unit tests.
    pub(crate) fn try_from_exact_energies(
        groups: impl IntoIterator<Item = (EdgeIndex, Vec<(usize, Atom)>)>,
    ) -> Result<Self, EnergyPowerAnalysisError> {
        Self::try_new(groups)
    }

    fn get(&self, edge: EdgeIndex) -> Option<&[usize]> {
        self.by_physical_edge.get(&edge).map(Vec::as_slice)
    }
}

/// One immutable assignment of factor-local physical energy dependencies to
/// certified equivalent exact energy variables.
///
/// The same plan owns both the bounds passed to CFF generation and the
/// factor-local assignments used when mapping the numerator. This prevents a
/// balanced bound from understating a numerator which is still evaluated
/// wholly through one canonical exact energy.
#[derive(Debug, Clone)]
pub(crate) struct EnergyPowerAssignmentPlan {
    expression: PlannedEnergyExpression,
    factor_assignments: BTreeMap<usize, BTreeMap<EdgeIndex, usize>>,
    energy_degree_bounds: Vec<(usize, usize)>,
}

impl EnergyPowerAssignmentPlan {
    pub(crate) fn energy_degree_bounds(&self) -> &[(usize, usize)] {
        &self.energy_degree_bounds
    }

    pub(crate) fn map_factors<E>(
        &self,
        mut map: impl FnMut(&Atom, &BTreeMap<EdgeIndex, usize>) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        self.expression.map(&self.factor_assignments, &mut map)
    }
}

#[derive(Debug, Clone)]
enum PlannedEnergyExpression {
    Factor {
        id: usize,
        expression: Atom,
        degrees: EnergyPowerCapMap,
    },
    Add(Vec<Self>),
    Mul(Vec<Self>),
    MultilinearFunction {
        symbol: Symbol,
        arguments: Vec<Self>,
    },
    Denominator {
        metadata: [Atom; 3],
        value: Box<Self>,
    },
}

impl PlannedEnergyExpression {
    fn degrees(&self) -> EnergyPowerCapMap {
        match self {
            Self::Factor { degrees, .. } => degrees.clone(),
            Self::Add(terms) => terms
                .iter()
                .fold(EnergyPowerCapMap::default(), |mut sum, term| {
                    sum.max_assign(term.degrees());
                    sum
                }),
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => factors
                .iter()
                .fold(EnergyPowerCapMap::default(), |mut product, factor| {
                    product.add_assign(factor.degrees());
                    product
                }),
            Self::Denominator { value, .. } => value.degrees(),
        }
    }

    fn assign_factor_slots(
        &self,
        offsets: &BTreeMap<EdgeIndex, usize>,
        slot_candidates: &BTreeMap<EdgeIndex, Vec<usize>>,
        factor_assignments: &mut BTreeMap<usize, BTreeMap<EdgeIndex, usize>>,
    ) -> Result<(), EnergyPowerAnalysisError> {
        match self {
            Self::Factor {
                id,
                expression,
                degrees,
            } => {
                for (edge, degree) in degrees.iter() {
                    if degree != 1 {
                        return Err(EnergyPowerAnalysisError::UnfactorizedEnergyPower {
                            expression: expression.log_print(None),
                            edge: edge.into(),
                            degree,
                        });
                    }
                    let slot = offsets.get(&edge).copied().unwrap_or(0);
                    factor_assignments
                        .entry(*id)
                        .or_default()
                        .insert(edge, slot_candidates[&edge][slot]);
                }
                Ok(())
            }
            Self::Add(terms) => terms.iter().try_for_each(|term| {
                term.assign_factor_slots(offsets, slot_candidates, factor_assignments)
            }),
            Self::Mul(factors)
            | Self::MultilinearFunction {
                arguments: factors, ..
            } => {
                let mut child_offsets = offsets.clone();
                for factor in factors {
                    factor.assign_factor_slots(
                        &child_offsets,
                        slot_candidates,
                        factor_assignments,
                    )?;
                    for (edge, degree) in factor.degrees().iter() {
                        *child_offsets.entry(edge).or_insert(0) += degree;
                    }
                }
                Ok(())
            }
            Self::Denominator { value, .. } => {
                value.assign_factor_slots(offsets, slot_candidates, factor_assignments)
            }
        }
    }

    fn map<E>(
        &self,
        assignments: &BTreeMap<usize, BTreeMap<EdgeIndex, usize>>,
        map: &mut impl FnMut(&Atom, &BTreeMap<EdgeIndex, usize>) -> Result<Atom, E>,
    ) -> Result<Atom, E> {
        match self {
            Self::Factor { id, expression, .. } => {
                let empty = BTreeMap::new();
                map(expression, assignments.get(id).unwrap_or(&empty))
            }
            Self::Add(terms) => {
                terms.iter().try_fold(
                    Atom::Zero,
                    |sum, term| Ok(sum + term.map(assignments, map)?),
                )
            }
            Self::Mul(factors) => factors.iter().try_fold(Atom::one(), |product, factor| {
                Ok(product * factor.map(assignments, map)?)
            }),
            Self::MultilinearFunction { symbol, arguments } => {
                let mut builder = FunctionBuilder::new(*symbol);
                for argument in arguments {
                    let mapped = argument.map(assignments, map)?;
                    builder = builder.add_arg(mapped.as_view());
                }
                Ok(builder.finish())
            }
            Self::Denominator { metadata, value } => {
                let mapped = value.map(assignments, map)?;
                Ok(FunctionBuilder::new(GS.den)
                    .add_arg(metadata[0].as_view())
                    .add_arg(metadata[1].as_view())
                    .add_arg(metadata[2].as_view())
                    .add_arg(mapped.as_view())
                    .finish())
            }
        }
    }
}

#[derive(Debug, Error)]
pub enum EnergyPowerAnalysisError {
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
        "energy-dependent numerator factor `{expression}` retains unsplit degree {degree} in physical EMR edge {edge}"
    )]
    UnfactorizedEnergyPower {
        expression: String,
        edge: usize,
        degree: usize,
    },
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
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn for_physical_emr_edges(internal_edges: impl IntoIterator<Item = EdgeIndex>) -> Self {
        Self {
            loop_edges: None,
            internal_edges: Some(internal_edges.into_iter().collect()),
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn analyze_atom(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::StrictPolynomial,
        )
    }

    pub fn analyze_atom_upper_bound(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::ConservativeUpperBound,
        )
    }

    pub(crate) fn plan_atom_assignment(
        &self,
        expression: &Atom,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        let expected_degrees = self.analyze_atom(expression)?;
        let mut next_factor_id = 0;
        let planned = self.plan_view(expression.as_view(), &mut next_factor_id)?;
        debug_assert_eq!(planned.degrees(), expected_degrees);

        let mut factor_assignments = BTreeMap::<usize, BTreeMap<EdgeIndex, usize>>::new();
        let mut exact_bounds = BTreeMap::<usize, usize>::new();
        let mut slot_candidates = BTreeMap::<EdgeIndex, Vec<usize>>::new();
        for (edge, degree) in expected_degrees.iter() {
            let exact_candidates = candidates.get(edge).ok_or(
                EnergyPowerAnalysisError::MissingEquivalentEnergyCandidates {
                    edge: edge.into(),
                    degree,
                },
            )?;
            // Disjoint equivalent-candidate sets reduce the exact minimax to
            // quotient/remainder balancing. Giving the remainder to canonical
            // candidate IDs first fixes the lexicographic tie deterministically.
            let quotient = degree / exact_candidates.len();
            let remainder = degree % exact_candidates.len();
            let mut slots = Vec::with_capacity(degree);
            for (candidate_index, &candidate) in exact_candidates.iter().enumerate() {
                let load = quotient + usize::from(candidate_index < remainder);
                if load != 0 {
                    let previous = exact_bounds.insert(candidate, load);
                    debug_assert!(
                        previous.is_none(),
                        "certified exact-energy candidate sets are disjoint"
                    );
                }
                slots.extend(std::iter::repeat_n(candidate, load));
            }
            slot_candidates.insert(edge, slots);
        }
        planned.assign_factor_slots(&BTreeMap::new(), &slot_candidates, &mut factor_assignments)?;

        Ok(EnergyPowerAssignmentPlan {
            expression: planned,
            factor_assignments,
            energy_degree_bounds: exact_bounds.into_iter().collect(),
        })
    }

    fn plan_view(
        &self,
        expression: AtomView<'_>,
        next_factor_id: &mut usize,
    ) -> Result<PlannedEnergyExpression, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Add(add) => Ok(PlannedEnergyExpression::Add(
                add.iter()
                    .map(|term| self.plan_view(term, next_factor_id))
                    .collect::<Result<Vec<_>, _>>()?,
            )),
            AtomView::Mul(mul) => Ok(PlannedEnergyExpression::Mul(
                mul.iter()
                    .map(|factor| self.plan_view(factor, next_factor_id))
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
                Ok(PlannedEnergyExpression::Mul(
                    (0..exponent)
                        .map(|_| self.plan_view(base, next_factor_id))
                        .collect::<Result<Vec<_>, _>>()?,
                ))
            }
            AtomView::Fun(function)
                if function.get_nargs() == 4 && function.get_symbol() == GS.den =>
            {
                Ok(PlannedEnergyExpression::Denominator {
                    metadata: [
                        function.get(0).to_owned(),
                        function.get(1).to_owned(),
                        function.get(2).to_owned(),
                    ],
                    value: Box::new(self.plan_view(function.get(3), next_factor_id)?),
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
                        .map(|argument| self.plan_view(argument, next_factor_id))
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
        Ok(PlannedEnergyExpression::Factor {
            id,
            expression: expression.to_owned(),
            degrees: self.analyze_view(expression, EnergyPowerAnalysisMode::StrictPolynomial)?,
        })
    }

    fn analyze_view(
        &self,
        expression: AtomView<'_>,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
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
                    product_degree.add_assign(self.analyze_view(factor, mode)?);
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
                base_degree.scale(exponent as usize);
                Ok(base_degree)
            }
            AtomView::Fun(function) => {
                if function.get_symbol() == GS.emr_vec {
                    return Ok(EnergyPowerCapMap::default());
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.emr_mom {
                    let edge = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
                        EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                            argument: function.get(0).to_owned().log_print(None),
                        }
                    })?);
                    return Ok(if self.is_internal_edge(edge) {
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
                    return Ok(if self.is_internal_edge(edge) {
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
                    degree.add_assign(self.analyze_view(function.get(1), mode)?);
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
                        degree.add_assign(argument_degree);
                    }
                    Ok(degree)
                } else {
                    Ok(EnergyPowerCapMap::default())
                }
            }
        }
    }

    fn component_degree(&self, edge: EdgeIndex, index: AtomView<'_>) -> EnergyPowerCapMap {
        match self.lorentz_index_kind(index) {
            LorentzIndexKind::Temporal | LorentzIndexKind::Abstract => {
                EnergyPowerCapMap::unit(edge)
            }
            LorentzIndexKind::Spatial => EnergyPowerCapMap::default(),
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

    fn is_internal_edge(&self, edge: EdgeIndex) -> bool {
        self.internal_edges
            .as_ref()
            .is_none_or(|internal_edges| internal_edges.contains(&edge))
    }
}

impl Graph {
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
            .analyze_view(numerator.as_view(), mode)
    }

    pub fn automatic_numerator_energy_degree_bounds(
        &self,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        Ok(self
            .numerator_energy_power_upper_bounds()?
            .into_generation_bounds())
    }

    pub(crate) fn automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
        &self,
        numerator: &Atom,
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
        Ok(EnergyPowerAnalyzer::for_physical_emr_edges(active_edges)
            .analyze_atom(numerator)?
            .iter()
            .filter_map(|(edge, degree)| (degree >= min_degree).then_some((edge.into(), degree)))
            .collect())
    }

    pub(crate) fn plan_numerator_energy_assignment_in_atom_excluding(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        candidates: &EquivalentEnergyCandidates,
    ) -> Result<EnergyPowerAssignmentPlan, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        EnergyPowerAnalyzer::for_physical_emr_edges(active_edges)
            .plan_atom_assignment(numerator, candidates)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        numerator::energy_degree::{
            EnergyPowerAnalysisError, EnergyPowerAnalyzer, EquivalentEnergyCandidates,
        },
        utils::GS,
    };
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, Minkowski},
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
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
        EquivalentEnergyCandidates::try_new([(
            EdgeIndex(3),
            vec![(10, "same exact energy"), (11, "same exact energy")],
        )])
        .unwrap()
    }

    fn planned_candidates(expression: &Atom) -> (Vec<(usize, usize)>, Vec<usize>, Atom) {
        let edge = EdgeIndex(3);
        let plan = EnergyPowerAnalyzer::for_physical_emr_edges([edge])
            .plan_atom_assignment(expression, &two_equivalent_candidates())
            .unwrap();
        let bounds = plan.energy_degree_bounds().to_vec();
        let mut assignments = Vec::new();
        let mapped = plan
            .map_factors(|factor, factor_assignments| {
                if let Some(candidate) = factor_assignments.get(&edge) {
                    assignments.push(*candidate);
                }
                Ok::<_, ()>(factor.clone())
            })
            .unwrap();
        (bounds, assignments, mapped)
    }

    #[test]
    fn minimax_assignment_balances_two_linear_factors() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression = (q3.clone() + Atom::num(1)) * (q3 + Atom::num(2));
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1), (11, 1)]);
        assert_eq!(assignments, vec![10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn minimax_assignment_balances_multilinear_function_slots() {
        let q3 = function!(GS.emr_mom, 3, compact_minkowski_vector());
        let expression = function!(GS.dot, q3.clone(), q3);
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1), (11, 1)]);
        assert_eq!(assignments, vec![10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn minimax_assignment_reuses_capacity_across_additive_branches() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let left = FunctionBuilder::new(symbol!("energy_assignment_test::left"; Linear))
            .add_arg(q3.as_view())
            .finish();
        let right = FunctionBuilder::new(symbol!("energy_assignment_test::right"; Linear))
            .add_arg(q3.as_view())
            .finish();
        let expression = left + right;
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1)]);
        assert_eq!(assignments, vec![10, 10]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn minimax_assignment_repeats_power_bases_without_expansion() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression =
            (q3 + Atom::var(symbol!("energy_assignment_test::shift"))).pow(Atom::num(2));
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 1), (11, 1)]);
        assert_eq!(assignments, vec![10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn minimax_assignment_uses_canonical_lexicographic_tie_break() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression =
            (q3.clone() + Atom::num(1)) * (q3.clone() + Atom::num(2)) * (q3 + Atom::num(3));
        let (bounds, assignments, mapped) = planned_candidates(&expression);

        assert_eq!(bounds, vec![(10, 2), (11, 1)]);
        assert_eq!(assignments, vec![10, 10, 11]);
        assert_eq!(mapped, expression);
    }

    #[test]
    fn exact_energy_candidates_select_the_largest_class_with_a_canonical_tie_break() {
        let candidates = EquivalentEnergyCandidates::try_new([
            (
                EdgeIndex(3),
                vec![(12, "negative"), (10, "positive"), (11, "negative")],
            ),
            (
                EdgeIndex(7),
                vec![(22, "negative"), (21, "positive"), (20, "positive")],
            ),
            (
                EdgeIndex(9),
                vec![
                    (33, "later"),
                    (31, "earlier"),
                    (32, "later"),
                    (30, "earlier"),
                ],
            ),
        ])
        .unwrap();

        assert_eq!(candidates.get(EdgeIndex(3)), Some([11, 12].as_slice()));
        assert_eq!(candidates.get(EdgeIndex(7)), Some([20, 21].as_slice()));
        assert_eq!(candidates.get(EdgeIndex(9)), Some([30, 31].as_slice()));
    }

    #[test]
    fn exact_energy_candidates_require_disjoint_occurrences_before_class_selection() {
        let overlapping = EquivalentEnergyCandidates::try_new([
            (
                EdgeIndex(3),
                vec![(10, "selected"), (11, "selected"), (12, "unselected")],
            ),
            (EdgeIndex(7), vec![(12, "other edge")]),
        ])
        .unwrap_err();
        assert!(matches!(
            overlapping,
            EnergyPowerAnalysisError::OverlappingEquivalentEnergyCandidates { .. }
        ));

        let duplicate = EquivalentEnergyCandidates::try_new([(
            EdgeIndex(3),
            vec![(10, "selected"), (11, "selected"), (11, "unselected")],
        )])
        .unwrap_err();
        assert!(matches!(
            duplicate,
            EnergyPowerAnalysisError::DuplicateEquivalentEnergyCandidate { .. }
        ));
    }

    #[test]
    fn assignment_rejects_a_missing_certified_candidate_set() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let candidates =
            EquivalentEnergyCandidates::try_new([(EdgeIndex(7), vec![(12, "unrelated energy")])])
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
