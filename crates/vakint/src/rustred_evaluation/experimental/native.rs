//! Thin steering of RustRed-owned search and per-integral application.

use std::collections::BTreeSet;
use std::sync::{Arc, Mutex};

use rustred::family::{IntegralFamily, IntegralKey};
use rustred::reduction::ReductionLimits;
use rustred::sector::{CoordinatePriority, CoordinatePriorityLimits, Mask, OrderingPolicy, zero};
use rustred::solver::{
    CandidateReducer, IntegralOrder, SectorConfig, SectorExecutor, SectorSolveOptions, SourceSystem,
};
use symbolica::atom::Atom;

use super::{CandidateReduction, CandidateScalarReduction};

/// Owns candidate rules, not a `ClosedArtifact` and not a family-closure proof.
#[derive(Debug)]
pub struct NativeCandidate<const N: usize> {
    reducer: Mutex<CandidateReducer<N>>,
    parent_momenta: Vec<Atom>,
    dimension: Atom,
    terminals: BTreeSet<IntegralKey>,
}

impl<const N: usize> NativeCandidate<N> {
    /// Solve an explicitly supplied family with the existing RustRed engine.
    ///
    /// The bounded experiment enumerates sectors only through arity 16.
    /// Parent momenta are the family's physical slots, not topology names.
    /// No algebra, source generation or rule application is implemented here.
    pub fn solve(
        family: &IntegralFamily,
        parent_momenta: Vec<Atom>,
        workers: usize,
        permutation: Option<[usize; N]>,
    ) -> Result<Self, String> {
        if N > 16
            || family.denominator_count() != N
            || parent_momenta.is_empty()
            || parent_momenta.len() > N
        {
            return Err("experimental candidate family/physical-slot arity is unsupported".into());
        }
        let root: [bool; N] = std::array::from_fn(|axis| axis < parent_momenta.len());
        let ordering = candidate_ordering(permutation)?;
        let analyzer = zero::Analyzer::try_unrestricted(family).map_err(|e| e.to_string())?;
        let mut certificates = Vec::new();
        let mut zero_sectors = Vec::new();
        let mut sectors = Vec::new();
        for bits in 0..(1usize << N) {
            let sector: [bool; N] = std::array::from_fn(|axis| bits & (1 << axis) != 0);
            let mask = Mask::try_new(sector).map_err(|e| e.to_string())?;
            match analyzer.analyze(&mask).map_err(|e| e.to_string())? {
                zero::Decision::ProvedZero(certificate) => {
                    zero_sectors.push(sector);
                    certificates.push(certificate);
                }
                zero::Decision::Inconclusive(_) => {
                    if sector
                        .iter()
                        .zip(root)
                        .all(|(&active, allowed)| !active || allowed)
                    {
                        sectors.push(sector);
                    }
                }
                zero::Decision::Excluded(_) => {
                    return Err("unrestricted sector was excluded".into());
                }
            }
        }
        sectors.sort_unstable();
        let config = SectorConfig {
            permutation,
            zero_sectors: Arc::from(zero_sectors),
            ..SectorConfig::default()
        };
        let source = SourceSystem::<N>::from_family(family).map_err(|e| e.to_string())?;
        let solutions = SectorExecutor::new(workers)
            .map_err(|e| e.to_string())?
            .map(
                &source,
                &sectors,
                &config,
                SectorSolveOptions::default(),
                |done| Ok::<_, String>((done.sector, done.solution)),
            )
            .map_err(|e| e.to_string())?;
        let reducer = CandidateReducer::try_new(
            family,
            root,
            ordering,
            solutions,
            certificates,
            ReductionLimits::default(),
        )
        .map_err(|e| e.to_string())?;
        let terminals = reducer.terminals().clone();
        Ok(Self {
            reducer: Mutex::new(reducer),
            parent_momenta,
            dimension: family.dimension().to_expression(),
            terminals,
        })
    }

    pub fn terminals(&self) -> &BTreeSet<IntegralKey> {
        &self.terminals
    }

    /// Drop pointwise memoized results without rerunning candidate search.
    /// This is useful for distinguishing a cold native application experiment
    /// from the earlier traversal used to discover its offline terminal set.
    pub fn clear_cache(&self) -> Result<(), String> {
        self.reducer
            .lock()
            .map_err(|_| "candidate reducer lock poisoned")?
            .clear_cache()
            .map_err(|error| error.to_string())
    }

    /// Cumulative formula applications; clearing the cache preserves counters.
    pub fn rule_applications(&self) -> Result<usize, String> {
        Ok(self
            .reducer
            .lock()
            .map_err(|_| "candidate reducer lock poisoned")?
            .statistics()
            .rule_applications())
    }
}

fn candidate_ordering<const N: usize>(
    permutation: Option<[usize; N]>,
) -> Result<OrderingPolicy, String> {
    let Some(slots) = permutation else {
        return Ok(OrderingPolicy::SpiredUncutV1);
    };
    // Solver permutations list slots in priority order; the durable ordering
    // descriptor instead stores the inverse rank for each slot.
    IntegralOrder::new([false; N], [false; N])
        .with_permutation(slots)
        .map_err(|error| error.to_string())?;
    let mut ranks = [0; N];
    for (rank, slot) in slots.into_iter().enumerate() {
        ranks[slot] = rank;
    }
    OrderingPolicy::try_spired_with_coordinate_priority(
        &CoordinatePriority::try_new(N, &ranks, CoordinatePriorityLimits::default())
            .map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_ordering_inverts_non_involutive_solver_permutation() {
        let expected = OrderingPolicy::try_spired_with_coordinate_priority(
            &CoordinatePriority::try_new(3, &[1, 2, 0], CoordinatePriorityLimits::default())
                .unwrap(),
        )
        .unwrap();
        assert_eq!(candidate_ordering(Some([2, 0, 1])).unwrap(), expected);
        assert_eq!(
            candidate_ordering::<3>(None).unwrap(),
            OrderingPolicy::SpiredUncutV1
        );
        assert!(candidate_ordering(Some([0, 0, 2])).is_err());
        assert!(candidate_ordering(Some([0, 1, 3])).is_err());
    }
}

impl<const N: usize> CandidateScalarReduction for NativeCandidate<N> {
    fn index_count(&self) -> usize {
        N
    }
    fn parent_momenta(&self) -> &[Atom] {
        &self.parent_momenta
    }
    fn dimension(&self) -> &Atom {
        &self.dimension
    }
    fn reduce_unit_mass(&self, target: &IntegralKey) -> Result<CandidateReduction, String> {
        let mut reducer = self
            .reducer
            .lock()
            .map_err(|_| "candidate reducer lock poisoned")?;
        let before = reducer.statistics().rule_applications();
        let decomposition = reducer
            .reduce_unit_mass(target)
            .map_err(|e| e.to_string())?;
        let applied_rules = reducer
            .statistics()
            .rule_applications()
            .checked_sub(before)
            .ok_or("candidate rule counter decreased")?;
        Ok(CandidateReduction {
            terms: decomposition
                .terms()
                .iter()
                .map(|(key, coefficient)| (key.clone(), coefficient.to_expression()))
                .collect(),
            applied_rules,
        })
    }
}
