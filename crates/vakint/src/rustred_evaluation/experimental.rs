//! Opt-in, per-integral comparisons of unsealed RustRed candidate rules.
//!
//! This development seam is deliberately absent from `EvaluationMethod` and
//! the shipped artifact registry. Successful numerical comparisons certify
//! neither a whole family nor arbitrary indices. Rule selection, pointwise
//! guards, descent, memoization and unresolved-leaf errors belong to RustRed.
//! Vakint only consumes its existing match/routing witness and an explicitly
//! supplied, offline terminal-value catalog. There is no FORM fallback here.

pub mod catalog;
pub mod native;

use std::collections::BTreeMap;
use std::fmt::Debug;
use std::sync::Arc;

use rustred::family::IntegralKey;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use crate::fmft::FMFT;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    FMFTOptions, Vakint, VakintError, VakintExpression, VakintSettings, get_integer_from_atom,
    get_prop_with_id,
};

/// Exact output from RustRed's experimental per-target rule application.
#[derive(Debug)]
pub struct CandidateReduction {
    /// Unit-mass coefficients in the same Minkowski convention as the family.
    pub terms: Vec<(IntegralKey, Atom)>,
    /// New applications during this call. Zero can also mean a cache hit;
    /// terminal membership must be checked independently.
    pub applied_rules: usize,
}

/// Narrow steering boundary; implementations must delegate to RustRed.
///
/// Families use `q² - 1` denominators. Physical slots precede auxiliary slots,
/// which enter a numerator-free target with power zero. The parent descriptor
/// is input data, not a topology-name dispatch. A reducer must reject unknown
/// leaves rather than promoting them to new terminals.
pub trait CandidateScalarReduction: Debug + Send + Sync {
    fn index_count(&self) -> usize;
    fn parent_momenta(&self) -> &[Atom];
    fn dimension(&self) -> &Atom;
    fn reduce_unit_mass(&self, target: &IntegralKey) -> Result<CandidateReduction, String>;
    fn lower_scalar_numerator(
        &self,
        numerator: &Atom,
        base: &IntegralKey,
    ) -> Result<Vec<CandidateLoweredTerm>, String>;
}

/// One exact term emitted by RustRed's family-bound scalar-numerator service.
/// This adapter only transports the service result; it does not expand or
/// simplify scalar products itself.
#[derive(Clone, Debug)]
pub struct CandidateLoweredTerm {
    pub target: IntegralKey,
    pub coefficient: Atom,
    pub scalar_spectator: Atom,
    pub common_mass_squared_power: u32,
}

/// Observable finite-target coverage, separate from family closure.
#[derive(Debug)]
pub struct CandidateEvaluation {
    pub value: Atom,
    pub targets: usize,
    pub applied_rules: usize,
}

/// Existing matcher output bound to the supplied candidate parent descriptor.
#[derive(Debug)]
pub struct PreparedCandidateIntegral {
    pub target: IntegralKey,
    mass_squared: Atom,
    scalar_numerator: Atom,
    coefficient: Atom,
    common_mass_squared_power: u32,
}

/// Prepare actual matched target keys before an offline terminal catalog exists.
/// The native evaluator uses this same helper; preparation never guesses the
/// key or rematches against a separate graph/family implementation.
pub fn prepare_candidate_integrals(
    vakint: &Vakint,
    settings: &VakintSettings,
    input: AtomView,
    reducer: &dyn CandidateScalarReduction,
) -> Result<Vec<PreparedCandidateIntegral>, VakintError> {
    if reducer.parent_momenta().is_empty() || reducer.parent_momenta().len() > reducer.index_count()
    {
        return Err(candidate_error(
            "inconsistent candidate physical-slot arity",
        ));
    }
    let mut prepared = Vec::new();
    for mut term in VakintExpression::split_integrals(input)? {
        let mut matched = vakint
            .topologies
            .match_topologies_to_user_input(term.integral.as_view(), false)?
            .ok_or_else(|| candidate_error("input has no registered topology match"))?;
        matched.apply_replacement_rules()?;
        term.apply_numerator_replacement_rules(&matched, settings)?;
        let integral = matched.canonical_topology.get_integral();
        if integral.n_loops != 4 {
            return Err(candidate_error(
                "experimental FMFT finalization requires four loops",
            ));
        }
        integral.validate_parent_routing(reducer.parent_momenta())?;
        let canonical_loop_momenta = (1..=integral.n_loops)
            .map(|axis| function!(S.k, Atom::num(axis)))
            .collect::<Vec<_>>();
        let routed_numerator = super::numerator::route_to_parent(
            term.numerator.as_view(),
            &canonical_loop_momenta,
            integral
                .parent_routing
                .as_ref()
                .map(|(_, coordinates)| coordinates.as_ref()),
        );
        let family_numerator = (1..=integral.n_loops).fold(routed_numerator, |value, axis| {
            value
                .replace(function!(S.k, Atom::num(axis)).to_pattern())
                .with(
                    vk_parse!(format!("k{axis}"))
                        .expect("family loop label parses")
                        .to_pattern(),
                )
        });
        let expression = integral
            .canonical_expression
            .as_ref()
            .ok_or_else(|| candidate_error("matched canonical expression is absent"))?;
        let mut powers = vec![0; reducer.index_count()];
        let mut mass: Option<Atom> = None;
        for slot in 1..=integral.n_props {
            let Some(properties) = get_prop_with_id(expression.as_view(), slot) else {
                continue;
            };
            powers[slot - 1] = properties
                .get(&vk_symbol!("pow_"))
                .and_then(|power| get_integer_from_atom(power.as_view()))
                .ok_or_else(|| candidate_error("noninteger matched propagator power"))?;
            let current = properties
                .get(&vk_symbol!("mUVsq_"))
                .ok_or_else(|| candidate_error("matched propagator mass is absent"))?;
            if current.is_zero() || mass.as_ref().is_some_and(|mass| mass != current) {
                return Err(candidate_error(
                    "candidate input requires one nonzero common mass",
                ));
            }
            mass.get_or_insert_with(|| current.clone());
        }
        let base =
            IntegralKey::try_new(powers).map_err(|error| candidate_error(error.to_string()))?;
        let lowered = reducer
            .lower_scalar_numerator(&family_numerator, &base)
            .map_err(candidate_error)?;
        if lowered.is_empty() {
            return Err(candidate_error(
                "scalar-numerator lowering returned no terms",
            ));
        }
        let mass_squared = mass.ok_or_else(|| candidate_error("candidate input has no mass"))?;
        prepared.extend(lowered.into_iter().map(|term| PreparedCandidateIntegral {
            target: term.target,
            mass_squared: mass_squared.clone(),
            scalar_numerator: term.scalar_spectator,
            coefficient: term.coefficient,
            common_mass_squared_power: term.common_mass_squared_power,
        }));
    }
    Ok(prepared)
}

/// Development-only adapter with an explicit offline FMFT terminal catalog.
#[derive(Debug)]
pub struct ExperimentalRustRed {
    reducer: Arc<dyn CandidateScalarReduction>,
    /// Exact unit-mass, pre-normalization Minkowski FMFT PR expressions.
    terminals: BTreeMap<IntegralKey, Atom>,
    master_options: FMFTOptions,
}

impl ExperimentalRustRed {
    pub fn new(
        reducer: Arc<dyn CandidateScalarReduction>,
        terminals: BTreeMap<IntegralKey, Atom>,
    ) -> Result<Self, VakintError> {
        if reducer.parent_momenta().is_empty()
            || reducer.parent_momenta().len() > reducer.index_count()
            || terminals
                .keys()
                .any(|key| key.powers().len() != reducer.index_count())
        {
            return Err(candidate_error(
                "inconsistent candidate family/catalog arity",
            ));
        }
        for (key, value) in &terminals {
            let mut approximate = false;
            value.visitor(&mut |part| {
                if let AtomView::Num(number) = part {
                    approximate |= number.get_coeff_view().is_float();
                }
                !approximate
            });
            if approximate {
                return Err(candidate_error(format!(
                    "terminal {:?} contains an approximate coefficient; the offline PR catalog must be exact",
                    key.powers()
                )));
            }
        }
        Ok(Self {
            reducer,
            terminals,
            master_options: FMFTOptions::default(),
        })
    }

    /// Select the existing pure-Rust master finalization options.
    ///
    /// For example, retaining symbolic master coefficients permits inspection
    /// before approximate numerical substitution. This does not select a FORM
    /// scalar backend or change candidate reduction.
    pub fn with_master_options(mut self, options: FMFTOptions) -> Self {
        self.master_options = options;
        self
    }

    /// Consume the existing matcher once, then run a FORM-free scalar tail.
    ///
    /// Scalar-product lowering is delegated to RustRed's family-bound service;
    /// this adapter only applies the resulting integral-key decomposition.
    pub fn evaluate_integral(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<CandidateEvaluation, VakintError> {
        let terms = prepare_candidate_integrals(vakint, settings, input, self.reducer.as_ref())?;
        let mut value = Atom::Zero;
        let mut applied_rules = 0usize;
        for term in &terms {
            let reduction = self
                .reducer
                .reduce_unit_mass(&term.target)
                .map_err(candidate_error)?;
            applied_rules = applied_rules
                .checked_add(reduction.applied_rules)
                .ok_or_else(|| candidate_error("candidate application counter overflow"))?;
            let mut raw = Atom::Zero;
            for (master, coefficient) in reduction.terms {
                let terminal = self.terminals.get(&master).ok_or_else(|| {
                    candidate_error(format!(
                        "no supplied offline terminal value for {:?}",
                        master.powers()
                    ))
                })?;
                raw += coefficient * terminal;
            }
            // (sum(master)-sum(target)) + (2L-sum(master)) = 2L-sum(target).
            // The shared FMFT finalizer supplies the remaining (m²)^(-L ep).
            let exponent = term
                .target
                .powers()
                .iter()
                .try_fold(8i128, |sum, power| sum.checked_sub(i128::from(*power)))
                .and_then(|sum| i64::try_from(sum).ok())
                .ok_or_else(|| candidate_error("candidate mass exponent overflow"))?;
            let dimension = vk_parse!("4-2*ep").expect("fixed dimension parses");
            raw = raw
                .replace(self.reducer.dimension().to_pattern())
                .with(dimension.to_pattern())
                // Different declared residuals may project onto the same PR
                // master. Combine its exact rational coefficient before any
                // Laurent truncation, so intermediate poles cannot request
                // unknown master orders that cancel in the complete sum.
                .together();
            let exponent = exponent
                .checked_add(i64::from(term.common_mass_squared_power))
                .ok_or_else(|| candidate_error("lowered mass exponent overflow"))?;
            raw *= term.mass_squared.clone().pow(Atom::num(exponent))
                * &term.coefficient
                * &term.scalar_numerator;
            value += FMFT::with_settings(settings.clone()).finalize_native_reduced_masters(
                raw,
                &term.mass_squared,
                &self.master_options,
            )?;
        }
        Ok(CandidateEvaluation {
            value,
            targets: terms.len(),
            applied_rules,
        })
    }
}

fn candidate_error(detail: impl Into<String>) -> VakintError {
    VakintError::RustRedEvaluation(super::RustRedEvaluationError::Reduction {
        detail: format!("experimental candidate: {}", detail.into()),
    })
}
