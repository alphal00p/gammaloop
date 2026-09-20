//! Thin steering of RustRed-owned search and per-integral application.

use std::collections::BTreeSet;
use std::sync::{Arc, Mutex};

use rustred::family::{IntegralFamily, IntegralKey};
use rustred::persistence::BinaryIoLimits;
use rustred::reduction::ReductionLimits;
use rustred::reduction::terminal_normalization::{
    TerminalAliasPlan, TerminalNormalizationLimits, TerminalNormalizationPlan,
    VacuumParametricLimits,
};
use rustred::sector::{CoordinatePriority, CoordinatePriorityLimits, Mask, OrderingPolicy, zero};
use rustred::solver::{
    CandidateReducer, IntegralOrder, SectorConfig, SectorExecutor, SectorSolveOptions, SourceSystem,
};
use symbolica::atom::Atom;

use super::{CandidateReduction, CandidateScalarReduction};
use crate::symbols::S;
use crate::vakint_parse;

/// Owns candidate rules, not a `ClosedArtifact` and not a family-closure proof.
#[derive(Debug)]
pub struct NativeCandidate<const N: usize> {
    reducer: Mutex<CandidateReducer<N>>,
    family: Arc<IntegralFamily>,
    parent_momenta: Vec<Atom>,
    dimension: Atom,
    terminals: BTreeSet<IntegralKey>,
    output_terminals: BTreeSet<IntegralKey>,
}

impl<const N: usize> NativeCandidate<N> {
    /// Steer an existing RustRed candidate applier without running discovery.
    ///
    /// The caller may load it with RustRed's candidate-bundle loader. Its
    /// formulas remain explicitly uncertified; loading does not promote them
    /// to a closed artifact. Binding the supplied family here ensures scalar
    /// numerator lowering uses the same denominator coordinates as reduction.
    pub fn from_reducer(
        family: Arc<IntegralFamily>,
        parent_momenta: Vec<Atom>,
        reducer: CandidateReducer<N>,
    ) -> Result<Self, String> {
        if family.denominator_count() != N
            || parent_momenta.is_empty()
            || parent_momenta.len() > N
            || family.fingerprint() != reducer.family_fingerprint()
        {
            return Err("candidate reducer/family binding or physical-slot arity differs".into());
        }
        let terminals = reducer.terminals().clone();
        let output_terminals = reducer.canonical_terminals().clone();
        Ok(Self {
            reducer: Mutex::new(reducer),
            dimension: family.dimension().to_expression(),
            family,
            parent_momenta,
            terminals,
            output_terminals,
        })
    }

    /// Bind a fresh owner and prepare RustRed's exact vacuum terminal aliases.
    ///
    /// Raw declarations still bind the rule program. The plan is prepared
    /// once using native U-polynomial equality, then used by the core applier
    /// before memoization. A previously used reducer must have its point cache
    /// explicitly cleared by the caller;
    /// unlike `from_reducer`, this constructor changes output representatives.
    pub fn from_reducer_with_terminal_aliases(
        family: Arc<IntegralFamily>,
        parent_momenta: Vec<Atom>,
        reducer: CandidateReducer<N>,
    ) -> Result<Self, String> {
        let mut native = Self::from_reducer(family, parent_momenta, reducer)?;
        let reducer = native
            .reducer
            .get_mut()
            .map_err(|_| "candidate reducer lock poisoned")?;
        let aliases = TerminalAliasPlan::vacuum_parametric_equivalences(
            &native.family,
            &native.terminals,
            reducer.ordering(),
            VacuumParametricLimits::default(),
        )
        .map_err(|error| error.to_string())?;
        reducer
            .install_terminal_aliases(aliases)
            .map_err(|error| error.to_string())?;
        native.output_terminals = reducer.canonical_terminals().clone();
        Ok(native)
    }

    /// Load a generated finite output convention into a fresh native owner.
    ///
    /// RustRed reconstructs the exact terminal proof and compares the complete
    /// stored output once, then owns weighted application and memoization.
    /// The raw terminal set still binds the rule program. As with explicit
    /// aliases, the caller must clear an already populated point cache itself.
    /// Native bytes must come from the matching trusted RustRed/Symbolica stack;
    /// this changes neither candidate authority nor the plain constructor.
    pub fn from_reducer_with_terminal_normalization(
        family: Arc<IntegralFamily>,
        parent_momenta: Vec<Atom>,
        reducer: CandidateReducer<N>,
        bytes: &[u8],
    ) -> Result<Self, String> {
        let mut native = Self::from_reducer(family, parent_momenta, reducer)?;
        let reducer = native
            .reducer
            .get_mut()
            .map_err(|_| "candidate reducer lock poisoned")?;
        let plan = TerminalNormalizationPlan::decode_generated(
            bytes,
            &native.family,
            reducer.terminals(),
            reducer.ordering(),
            TerminalNormalizationLimits::default(),
            BinaryIoLimits::default(),
        )
        .map_err(|error| error.to_string())?;
        reducer
            .install_terminal_normalization(plan)
            .map_err(|error| error.to_string())?;
        native.output_terminals = reducer.canonical_terminals().clone();
        Ok(native)
    }

    /// Solve an explicitly supplied family with the existing RustRed engine.
    ///
    /// The bounded experiment enumerates sectors only through arity 16.
    /// Parent momenta are the family's physical slots, not topology names.
    /// No algebra, source generation or rule application is implemented here.
    pub fn solve(
        family: Arc<IntegralFamily>,
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
        let analyzer = zero::Analyzer::try_unrestricted(&family).map_err(|e| e.to_string())?;
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
        let source = SourceSystem::<N>::from_family(&family).map_err(|e| e.to_string())?;
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
            &family,
            root,
            ordering,
            solutions,
            certificates,
            ReductionLimits::default(),
        )
        .map_err(|e| e.to_string())?;
        Self::from_reducer(family, parent_momenta, reducer)
    }

    pub fn terminals(&self) -> &BTreeSet<IntegralKey> {
        &self.terminals
    }

    /// Exact possible outputs of the installed RustRed terminal convention.
    /// Offline values need only cover these keys; raw declarations remain
    /// available through `terminals()` and owned by the rule program.
    pub fn output_terminals(&self) -> &BTreeSet<IntegralKey> {
        &self.output_terminals
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
    fn candidate_binding_installs_aliases_without_changing_raw_catalog_keys() {
        use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
        use rustred::solver::{Integral, SectorSolution, SectorStats};

        crate::Vakint::initialize_vakint_symbols();
        // Only declare three equal dotted terminals. This fixture has no rules
        // and never runs IBP discovery; the core owns the equality proof.
        let family = Arc::new(
            Compiler::new(Limits::default())
                .unwrap()
                .compile_text(TextProject {
                    name: None,
                    parameters: None,
                    loop_momenta: vec!["k1".into(), "k2".into()],
                    external_momenta: Vec::new(),
                    dimension: "d".into(),
                    propagators: ["k1^2-1", "k2^2-1", "(k1+k2)^2-1"]
                        .into_iter()
                        .enumerate()
                        .map(|(slot, expression)| TextPropagator {
                            id: format!("D{}", slot + 1),
                            expression: expression.into(),
                            target_power: 1,
                            power_shift: None,
                        })
                        .collect(),
                    external_gram: Vec::new(),
                    numerator: None,
                })
                .unwrap()
                .into_lowered(LoweringLimits::default())
                .unwrap()
                .into_family(),
        );
        let solution = SectorSolution {
            rules: Vec::new(),
            finite_residuals: [[2, 1, 1], [1, 2, 1], [1, 1, 2]]
                .into_iter()
                .map(|powers| Integral::numeric(powers).unwrap())
                .collect(),
            stats: SectorStats::default(),
        };
        let mut reducer = CandidateReducer::try_new(
            &family,
            [true; 3],
            OrderingPolicy::SpiredUncutV1,
            [([true; 3], solution)],
            Vec::new(),
            ReductionLimits::default(),
        )
        .unwrap();
        let raw = reducer.terminals().clone();
        assert!(reducer.terminal_aliases().is_none());
        let momenta = ["k(1)", "k(2)", "k(1)+k(2)"]
            .into_iter()
            .map(|momentum| vakint_parse!(momentum).unwrap())
            .collect();
        let raw_key = raw.first().unwrap().clone();
        reducer.reduce_unit_mass(&raw_key).unwrap();
        // The existing constructor must still accept populated reducers and
        // preserve their raw output convention without implicit normalization.
        let raw_native = NativeCandidate::from_reducer(family.clone(), momenta, reducer).unwrap();
        assert!(
            raw_native
                .reducer
                .lock()
                .unwrap()
                .terminal_aliases()
                .is_none()
        );
        assert_eq!(
            raw_native.reduce_unit_mass(&raw_key).unwrap().terms,
            vec![(raw_key, Atom::num(1))]
        );
        raw_native.clear_cache().unwrap();
        let native = NativeCandidate::from_reducer_with_terminal_aliases(
            family.clone(),
            raw_native.parent_momenta,
            raw_native.reducer.into_inner().unwrap(),
        )
        .unwrap();
        assert_eq!(native.terminals(), &raw);
        assert_eq!(native.output_terminals().len(), 1);
        let (source, representative) = {
            let reducer = native.reducer.lock().unwrap();
            let plan = reducer.terminal_aliases().unwrap();
            assert_eq!(plan.raw_terminals(), &raw);
            assert_eq!(plan.statistics().verified_aliases, 2);
            assert_eq!(reducer.canonical_terminals().len(), 1);
            let (source, alias) = plan.aliases().first_key_value().unwrap();
            (source.clone(), alias.representative().clone())
        };
        let output = native.reduce_unit_mass(&source).unwrap();
        assert_eq!(output.terms, vec![(representative, Atom::num(1))]);
        assert_eq!(output.applied_rules, 0);
        assert_eq!(
            native.reduce_unit_mass(&source).unwrap().terms,
            output.terms
        );
        native.clear_cache().unwrap();
        assert_eq!(native.terminals(), &raw);
        assert_eq!(native.output_terminals().len(), 1);
        assert_eq!(
            native.reduce_unit_mass(&source).unwrap().terms,
            output.terms
        );
        let error = NativeCandidate::from_reducer_with_terminal_aliases(
            family,
            native.parent_momenta,
            native.reducer.into_inner().unwrap(),
        )
        .unwrap_err();
        assert!(error.contains("empty reduction cache"));
    }

    #[test]
    fn candidate_binding_loads_weighted_terminal_outputs_without_discovery() {
        use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
        use rustred::solver::{Integral, SectorSolution, SectorStats};
        use std::collections::BTreeMap;

        crate::Vakint::initialize_vakint_symbols();
        let family = Arc::new(
            Compiler::new(Limits::default())
                .unwrap()
                .compile_text(TextProject {
                    name: None,
                    parameters: None,
                    loop_momenta: vec!["k1".into(), "k2".into(), "k3".into()],
                    external_momenta: Vec::new(),
                    dimension: "d".into(),
                    propagators: [
                        "k1^2-1",
                        "k2^2-1",
                        "k3^2-1",
                        "(k1+k2+k3)^2-1",
                        "(k1+k2)^2-1",
                        "(k1+k3)^2-1",
                    ]
                    .into_iter()
                    .enumerate()
                    .map(|(slot, expression)| TextPropagator {
                        id: format!("D{}", slot + 1),
                        expression: expression.into(),
                        target_power: i64::from(slot < 4),
                        power_shift: None,
                    })
                    .collect(),
                    external_gram: Vec::new(),
                    numerator: None,
                })
                .unwrap()
                .into_lowered(LoweringLimits::default())
                .unwrap()
                .into_family(),
        );
        let scalar = [1, 1, 1, 1, 0, 0];
        let numerator = [1, 1, 1, 1, -1, 0];
        let mut declarations = vec![scalar, numerator];
        for slot in 0..4 {
            let mut pinch = scalar;
            pinch[slot] = 0;
            declarations.push(pinch);
        }
        let mut solutions = BTreeMap::new();
        for powers in declarations {
            let sector = powers.map(|power| power > 0);
            solutions
                .entry(sector)
                .or_insert_with(|| SectorSolution {
                    rules: Vec::new(),
                    finite_residuals: Vec::new(),
                    stats: SectorStats::default(),
                })
                .finite_residuals
                .push(Integral::numeric(powers).unwrap());
        }
        let mut reducer = CandidateReducer::try_new(
            &family,
            [true, true, true, true, false, false],
            OrderingPolicy::SpiredUncutV1,
            solutions,
            Vec::new(),
            ReductionLimits::default(),
        )
        .unwrap();
        let raw = reducer.terminals().clone();
        let plan = TerminalNormalizationPlan::vacuum_quadratic_numerators(
            &family,
            &raw,
            reducer.ordering(),
            TerminalNormalizationLimits::default(),
        )
        .unwrap();
        assert_eq!(plan.statistics().projected_numerators, 1);
        assert_eq!(plan.canonical_terminals().len(), 2);
        let bytes = plan.encode_native(BinaryIoLimits::default()).unwrap();
        let numerator = IntegralKey::try_new(numerator.map(i64::from).to_vec()).unwrap();
        let scalar = IntegralKey::try_new(scalar.map(i64::from).to_vec()).unwrap();
        let pinch = plan
            .canonical_terminals()
            .iter()
            .find(|key| **key != scalar)
            .unwrap()
            .clone();
        let momenta = ["k(1)", "k(2)", "k(3)", "k(1)+k(2)+k(3)"]
            .into_iter()
            .map(|value| vakint_parse!(value).unwrap())
            .collect();
        reducer.reduce_unit_mass(&numerator).unwrap();
        let raw_native = NativeCandidate::from_reducer(family.clone(), momenta, reducer).unwrap();
        assert_eq!(
            raw_native.reduce_unit_mass(&numerator).unwrap().terms,
            vec![(numerator.clone(), Atom::num(1))]
        );
        raw_native.clear_cache().unwrap();
        let native = NativeCandidate::from_reducer_with_terminal_normalization(
            family.clone(),
            raw_native.parent_momenta,
            raw_native.reducer.into_inner().unwrap(),
            &bytes,
        )
        .unwrap();
        assert_eq!(native.terminals(), &raw);
        let expected = BTreeMap::from([
            (scalar, vakint_parse!("1/3").unwrap()),
            (pinch, vakint_parse!("4/3").unwrap()),
        ]);
        assert_eq!(
            native.output_terminals(),
            &expected.keys().cloned().collect()
        );
        assert!(!native.output_terminals().contains(&numerator));
        // The four-line circuit symmetry gives D5 -> (D1+D2+D3+D4+1)/3
        // after integration. All four positive pinches have the same value.
        for _ in 0..2 {
            let result = native.reduce_unit_mass(&numerator).unwrap();
            assert_eq!(
                result.terms.into_iter().collect::<BTreeMap<_, _>>(),
                expected
            );
            assert_eq!(result.applied_rules, 0);
        }
        assert!(
            native
                .reducer
                .lock()
                .unwrap()
                .terminal_normalization()
                .is_some()
        );
        let error = NativeCandidate::from_reducer_with_terminal_normalization(
            family,
            native.parent_momenta,
            native.reducer.into_inner().unwrap(),
            &bytes,
        )
        .unwrap_err();
        assert!(error.contains("empty reduction cache"));
    }

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

    fn lower_scalar_numerator(
        &self,
        numerator: &Atom,
        base: &IntegralKey,
    ) -> Result<Vec<super::CandidateLoweredTerm>, String> {
        let loop_momenta = self
            .family
            .loop_momenta()
            .iter()
            .map(|label| {
                vakint_parse!(label.as_str())
                    .map_err(|error| format!("parse family loop momentum {label}: {error}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let service = rustred::scalar_numerator::FamilyScalarNumeratorService::try_new(
            &self.family,
            S.dot,
            loop_momenta,
            rustred::scalar_numerator::ScalarNumeratorLimits::default(),
        )
        .map_err(|error| error.to_string())?;
        let lowering = service
            .lower(numerator, base)
            .map_err(|error| error.to_string())?;
        Ok(lowering
            .terms()
            .iter()
            .map(|term| super::CandidateLoweredTerm {
                target: term.integral().clone(),
                coefficient: term.coefficient().to_expression(),
                scalar_spectator: term.scalar_spectator().clone(),
                common_mass_squared_power: term.common_mass_squared_power(),
            })
            .collect())
    }
}
