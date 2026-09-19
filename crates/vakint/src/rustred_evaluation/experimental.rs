//! Opt-in, per-integral comparisons of unsealed RustRed candidate rules.
//!
//! This diagnostic seam accepts explicitly supplied candidate programs. The
//! shipped four-loop backend reuses its already-matched-input adapter.
//! Successful numerical comparisons certify neither a whole family nor
//! arbitrary indices. Rule selection, pointwise
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
    FMFTOptions, ReplacementRules, Vakint, VakintError, VakintExpression, VakintSettings,
    get_integer_from_atom, get_prop_with_id,
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

impl ReplacementRules {
    /// Ordered physical momenta from the existing four-loop match witness.
    ///
    /// This authenticates the retained defining parent and contracted routing
    /// before exposing a descriptor for finite candidate selection. It neither
    /// matches another topology nor selects a reducer from a topology name.
    pub fn candidate_parent_momenta(&mut self) -> Result<Vec<Atom>, VakintError> {
        self.apply_replacement_rules()?;
        let integral = self.canonical_topology.get_integral();
        if integral.n_loops != 4 {
            return Err(candidate_error("candidate parent requires four loops"));
        }
        let (parent, _) = integral
            .parent_routing
            .as_ref()
            .ok_or_else(|| candidate_error("matched family has no retained parent routing"))?;
        let momenta = (1..=integral.n_props)
            .map(|slot| {
                get_prop_with_id(parent.as_view(), slot)
                    .and_then(|properties| properties.get(&vk_symbol!("q_")).cloned())
                    .ok_or_else(|| {
                        candidate_error(format!("defining parent slot {slot} is absent"))
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        integral.validate_parent_routing(&momenta)?;
        Ok(momenta)
    }
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
        prepared.extend(prepare_matched_candidate_integral(
            term.numerator.as_view(),
            &matched,
            reducer,
        )?);
    }
    Ok(prepared)
}

/// Consume the match and simultaneous routing supplied by Vakint's dispatcher.
/// Unlike the diagnostic convenience entry point, this never rematches input.
fn prepare_matched_candidate_integral(
    numerator: AtomView,
    matched: &ReplacementRules,
    reducer: &dyn CandidateScalarReduction,
) -> Result<Vec<PreparedCandidateIntegral>, VakintError> {
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
        numerator,
        &canonical_loop_momenta,
        integral
            .parent_routing
            .as_ref()
            .map(|(_, coordinates)| coordinates.as_ref()),
    );
    // Parent routing expands vector sums at component level. Restore the
    // scalar-product representation before replacing vector labels, just
    // as the production materializer does. Otherwise k(i,mu) is not
    // matched by k(i), so loop components become opaque spectators in the
    // family-bound scalar lowerer.
    let scalar_numerator = Vakint::convert_to_dot_notation(routed_numerator.as_view());
    let family_numerator = (1..=integral.n_loops).fold(scalar_numerator, |value, axis| {
        value
            .replace(function!(S.k, Atom::num(axis)).to_pattern())
            .with(
                vk_parse!(format!("k{axis}"))
                    .expect("family loop label parses")
                    .to_pattern(),
            )
    });
    if family_numerator.contains_symbol(S.k) {
        return Err(candidate_error(
            "parent-routed numerator retains loop-vector components; tensor preprocessing must produce scalar products",
        ));
    }
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
    let base = IntegralKey::try_new(powers).map_err(|error| candidate_error(error.to_string()))?;
    let lowered = reducer
        .lower_scalar_numerator(&family_numerator, &base)
        .map_err(candidate_error)?;
    if lowered.is_empty() {
        return Err(candidate_error(
            "scalar-numerator lowering returned no terms",
        ));
    }
    let mass_squared = mass.ok_or_else(|| candidate_error("candidate input has no mass"))?;
    Ok(lowered
        .into_iter()
        .map(|term| PreparedCandidateIntegral {
            target: term.target,
            mass_squared: mass_squared.clone(),
            scalar_numerator: term.scalar_spectator,
            coefficient: term.coefficient,
            common_mass_squared_power: term.common_mass_squared_power,
        })
        .collect())
}

/// Unsealed-program adapter with an explicit offline FMFT terminal catalog.
///
/// Diagnostics supply their own program. The public four-loop backend uses
/// only the bundled, finite-tested programs selected from its structural
/// registry and passes the existing matcher witness directly.
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
        self.evaluate_prepared(settings, &terms, &self.master_options)
    }

    pub(super) fn evaluate_matched(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        matched: &ReplacementRules,
        substitute_masters: bool,
    ) -> Result<CandidateEvaluation, VakintError> {
        let terms = prepare_matched_candidate_integral(numerator, matched, self.reducer.as_ref())?;
        self.evaluate_prepared(
            settings,
            &terms,
            &FMFTOptions {
                expand_masters: substitute_masters,
                susbstitute_masters: substitute_masters,
            },
        )
    }

    fn evaluate_prepared(
        &self,
        settings: &VakintSettings,
        terms: &[PreparedCandidateIntegral],
        master_options: &FMFTOptions,
    ) -> Result<CandidateEvaluation, VakintError> {
        let mut contributions = BTreeMap::<Atom, Atom>::new();
        let mut applied_rules = 0usize;
        for term in terms {
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
            contributions
                .entry(term.mass_squared.clone())
                .and_modify(|sum| *sum += raw.clone())
                .or_insert(raw);
        }
        // Keep every exact contribution with the same common mass together.
        // Native FMFT finalization expands Laurent series and substitutes
        // finite master tables; doing that per target can expose an
        // unsupported pole before another target cancels it exactly.
        let value = finalize_candidate_contributions(settings, contributions, master_options)?;
        Ok(CandidateEvaluation {
            value,
            targets: terms.len(),
            applied_rules,
        })
    }
}

fn finalize_candidate_contributions(
    settings: &VakintSettings,
    contributions: BTreeMap<Atom, Atom>,
    options: &FMFTOptions,
) -> Result<Atom, VakintError> {
    let fmft = FMFT::with_settings(settings.clone());
    contributions
        .into_iter()
        .try_fold(Atom::Zero, |mut value, (mass, raw)| {
            value += fmft.finalize_native_reduced_masters(raw.together(), &mass, options)?;
            Ok(value)
        })
}

fn candidate_error(detail: impl Into<String>) -> VakintError {
    VakintError::RustRedEvaluation(super::RustRedEvaluationError::Reduction {
        detail: format!("experimental candidate: {}", detail.into()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_finalization_cancels_unknown_orders_before_expansion() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings {
            form_exe_path: "/definitely/not/a/form/executable".into(),
            number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: 25,
            integral_normalization_factor: crate::LoopNormalizationFactor::FMFTandMATAD,
            use_dot_product_notation: true,
            ..VakintSettings::default()
        };
        let mut contributions = BTreeMap::new();
        let mass = Atom::num(1);
        contributions.insert(mass.clone(), vk_parse!("PR9x").unwrap());
        contributions
            .entry(mass)
            .and_modify(|sum| *sum += vk_parse!("-PR9x").unwrap())
            .or_insert_with(|| vk_parse!("-PR9x").unwrap());

        let value =
            finalize_candidate_contributions(&settings, contributions, &FMFTOptions::default())
                .expect("exact cancellation should remove unsupported master orders");
        assert!(value.is_zero());
    }
}
