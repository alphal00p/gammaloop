//! Close a massless propagator integral with one massive propagator.
//!
//! Keep its exact homogeneity through the final radial integration: expanding
//! the p-integral first would turn the necessary epsilon-dependent power of
//! the hard momentum into logarithms and lose the simple Gamma integral.

use std::time::{Duration, Instant};

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::function;
use symbolica::transcendental::gamma;

use super::{ReducedPropagatorValue, external_head, forcer, laurent, scalar_numerator};
use crate::symbols::S;
use crate::utils::simplify_real;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{EvaluationMethod, RustRedMasterMode, Vakint, VakintError, VakintExpression};

/// The massive line used to close an already cut, scalar vacuum graph.
///
/// `momentum` is the massless graph's external `p(i)`. `evaluate` integrates
/// that momentum with an extra `(p² - mass_squared)^(-power)` and the same
/// per-loop measure as its p-integral. All other propagators must be massless.
/// The caller owns the graph cut and its momentum routing.
#[derive(Clone, Debug)]
pub struct MassiveClosure {
    pub momentum: Atom,
    pub mass_squared: Atom,
    pub power: i64,
}

impl MassiveClosure {
    fn validated_head(&self, settings: &crate::VakintSettings) -> Result<Atom, VakintError> {
        let (head, _) = external_head(&self.momentum)?;
        if self.mass_squared.contains_symbol(S.k)
            || self.mass_squared.contains_symbol(S.p)
            || self
                .mass_squared
                .contains_symbol(vk_symbol!(&settings.epsilon_symbol))
        {
            return Err(super::invalid(
                "the closing mass must be independent of momenta and epsilon",
            ));
        }
        if self.power <= 0 || self.mass_squared.is_zero() {
            return Err(super::invalid(
                "the closing propagator needs a positive power and a nonzero mass squared",
            ));
        }
        Ok(head)
    }

    /// Evaluate independent closures in bounded batches. Exact scalar
    /// reductions may share a FORCER invocation, but every closure keeps its
    /// own mass, power, external momentum and Laurent expansion. Spectator
    /// coefficients, including epsilon poles, remain attached to their
    /// original homogeneous channels until that expansion.
    pub fn evaluate_batch(
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        inputs: &[(Self, Atom)],
    ) -> Result<Vec<Atom>, VakintError> {
        let mut results = Vec::with_capacity(inputs.len());
        for chunk in inputs.chunks(32) {
            results
                .extend(settings.with_resource_budget(|| {
                    Self::evaluate_batch_inner(vakint, settings, chunk)
                })?);
        }
        Ok(results)
    }

    fn evaluate_batch_inner(
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        inputs: &[(Self, Atom)],
    ) -> Result<Vec<Atom>, VakintError> {
        let started = Instant::now();
        let mut preparation = Duration::ZERO;
        let mut finalization = Duration::ZERO;
        let mut combination = Duration::ZERO;
        let mut numerical = Duration::ZERO;
        let mut local_terms = 0_usize;
        let mut scaleless_terms = 0_usize;
        let mut source_terms = 0_usize;
        let mut cut_loops_max = 0_usize;
        let mode = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::Forcer(options) => Some(options.master_mode),
                _ => None,
            })
            .unwrap_or(RustRedMasterMode::Symbolic);
        let mut results = vec![Atom::Zero; inputs.len()];
        let mut heads = Vec::with_capacity(inputs.len());
        let mut requests = Vec::new();
        let mut owners = Vec::new();
        for (owner, (closure, input)) in inputs.iter().enumerate() {
            let preparation_started = Instant::now();
            let head = closure.validated_head(settings)?;
            heads.push(head.clone());
            if input.is_zero() {
                preparation += preparation_started.elapsed();
                continue;
            }
            settings.validate_execution()?;
            let mut topologies = Vec::new();
            input.visitor(&mut |part| {
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.topo
                {
                    topologies.push(part.to_owned());
                    return false;
                }
                true
            });
            topologies.sort();
            topologies.dedup();
            let coefficients = input.coefficient_list::<i16>(&topologies);
            source_terms += coefficients.len();
            preparation += preparation_started.elapsed();
            for (topology, numerator) in coefficients {
                let preparation_started = Instant::now();
                if topology == Atom::one() {
                    let values = [ReducedPropagatorValue {
                        scale_power: 0,
                        reduced: Vakint::convert_to_dot_notation(numerator.as_view()),
                    }];
                    preparation += preparation_started.elapsed();
                    let finalization_started = Instant::now();
                    let value = closure.finalize(&values, 0, &head, settings, mode)?;
                    finalization += finalization_started.elapsed();
                    let combination_started = Instant::now();
                    results[owner] += value;
                    combination += combination_started.elapsed();
                    local_terms += 1;
                    continue;
                }
                let mut expression = VakintExpression::try_from(topology * numerator)?;
                for term in &mut expression.0 {
                    if !term.normalize_propagator()? {
                        scaleless_terms += 1;
                        continue;
                    }
                    vakint.validate_settings(settings)?;
                    let external = term.propagator_external()?.ok_or_else(|| {
                        super::invalid("a massive closure requires a massless p-integral")
                    })?;
                    if external != head {
                        return Err(super::invalid(format!(
                            "the cut graph uses {external}, but the closing line uses {head}"
                        )));
                    }
                    term.numerator = scalar_numerator(term.numerator.as_view(), &head)?;
                    let matched = term.prepare_integral_specs(vakint, settings, None)?;
                    if !matches!(
                        vakint.evaluation_method(settings, &matched.canonical_topology),
                        Some(EvaluationMethod::Forcer(_))
                    ) {
                        return Err(super::invalid(
                            "massive closure of a p-integral requires the Forcer backend",
                        ));
                    }
                    let request = matched
                        .propagator_request(term.numerator.as_view())
                        .expect("a p-integral match has a propagator request");
                    if let Some(input) = forcer::forcer_input_with_epsilon(
                        request.family,
                        request.powers,
                        request.canonical_numerator,
                        &request.external,
                        vk_symbol!(&settings.epsilon_symbol),
                    )? && !input.parts.is_empty()
                    {
                        cut_loops_max = cut_loops_max.max(request.family.loops);
                        requests.push(input);
                        owners.push((owner, request.family.loops));
                    }
                }
                preparation += preparation_started.elapsed();
            }
        }
        let reduction_started = Instant::now();
        let values = forcer::reduce_inputs(&requests, settings)?;
        let reduction = reduction_started.elapsed();
        for ((owner, loops), values) in owners.iter().copied().zip(values) {
            let finalization_started = Instant::now();
            let value = inputs[owner]
                .0
                .finalize(&values, loops, &heads[owner], settings, mode)?;
            finalization += finalization_started.elapsed();
            let combination_started = Instant::now();
            results[owner] += value;
            combination += combination_started.elapsed();
        }
        for ((_, input), result) in inputs.iter().zip(&mut results) {
            if input.is_zero() {
                continue;
            }
            if mode == RustRedMasterMode::Numerical {
                let numerical_started = Instant::now();
                *result = vakint.resolve_rustred_numerical(std::mem::take(result), settings)?;
                numerical += numerical_started.elapsed();
            }
            let combination_started = Instant::now();
            *result = result.expand();
            combination += combination_started.elapsed();
        }
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            closure_requests = inputs.len(),
            source_terms = source_terms,
            local_terms = local_terms,
            scaleless_terms = scaleless_terms,
            propagator_terms = requests.len(),
            cut_loops_max = cut_loops_max,
            preparation_seconds = preparation.as_secs_f64(),
            forcer_reduction_seconds = reduction.as_secs_f64(),
            radial_master_laurent_seconds = finalization.as_secs_f64(),
            result_combination_seconds = combination.as_secs_f64(),
            numerical_resolution_seconds = numerical.as_secs_f64(),
            seconds = started.elapsed().as_secs_f64(),
            "Evaluated massive p-integral closure batch"
        );
        Ok(results)
    }

    /// Evaluate scalar p-integrals and loop-free powers of `dot(p,p)`, including
    /// sums of both. The configured Forcer backend reduces the p-integrals;
    /// loop-free terms need no external backend. Tensor dependence on the hard
    /// momentum is rejected rather than discarded.
    ///
    /// The Laurent depth counts the closing loop as well: a cut with `L` loops
    /// is expanded through `number_of_terms_in_epsilon_expansion - L - 2`.
    /// Forcer's master mode selects retained or expanded output; without a
    /// Forcer entry, loop-free terms use symbolic expansion.
    pub fn evaluate(
        &self,
        vakint: &Vakint,
        settings: &crate::VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            let started = Instant::now();
            let mut preparation = Duration::ZERO;
            let mut reduction = Duration::ZERO;
            let mut finalization = Duration::ZERO;
            let mut combination = Duration::ZERO;
            let mut numerical = Duration::ZERO;
            let mut local_terms = 0_usize;
            let mut scaleless_terms = 0_usize;
            let mut propagator_terms = 0_usize;
            let mut cut_loops_max = 0_usize;
            let preparation_started = Instant::now();
            let head = self.validated_head(settings)?;
            if input.is_zero() {
                return Ok(Atom::Zero);
            }
            settings.validate_execution()?;
            let mode = settings
                .evaluation_order
                .0
                .iter()
                .find_map(|method| match method {
                    EvaluationMethod::Forcer(options) => Some(options.master_mode),
                    _ => None,
                })
                .unwrap_or(RustRedMasterMode::Symbolic);

            // VakintExpression intentionally requires a topology in every term.
            // Separate the local remainder here so closing a tree or a local
            // counterterm also receives its final massive loop integration.
            let mut topologies = Vec::new();
            input.visitor(&mut |part| {
                if let AtomView::Fun(function) = part
                    && function.get_symbol() == S.topo
                {
                    topologies.push(part.to_owned());
                    return false;
                }
                true
            });
            topologies.sort();
            topologies.dedup();
            let coefficients = input.coefficient_list::<i16>(&topologies);
            let source_terms = coefficients.len();
            preparation += preparation_started.elapsed();
            let mut result = Atom::Zero;
            for (topology, numerator) in coefficients {
                if topology == Atom::one() {
                    let preparation_started = Instant::now();
                    let values = [ReducedPropagatorValue {
                        scale_power: 0,
                        reduced: Vakint::convert_to_dot_notation(numerator.as_view()),
                    }];
                    preparation += preparation_started.elapsed();
                    let finalization_started = Instant::now();
                    let integrated = self.finalize(&values, 0, &head, settings, mode)?;
                    finalization += finalization_started.elapsed();
                    let combination_started = Instant::now();
                    result += integrated;
                    combination += combination_started.elapsed();
                    local_terms += 1;
                    continue;
                }
                let preparation_started = Instant::now();
                let mut expression = VakintExpression::try_from(topology * numerator)?;
                preparation += preparation_started.elapsed();
                for term in &mut expression.0 {
                    let preparation_started = Instant::now();
                    if !term.normalize_propagator()? {
                        preparation += preparation_started.elapsed();
                        scaleless_terms += 1;
                        continue;
                    }
                    vakint.validate_settings(settings)?;
                    let external = term.propagator_external()?.ok_or_else(|| {
                        super::invalid("a massive closure requires a massless p-integral")
                    })?;
                    if external != head {
                        return Err(super::invalid(format!(
                            "the cut graph uses {external}, but the closing line uses {head}"
                        )));
                    }
                    term.numerator = scalar_numerator(term.numerator.as_view(), &head)?;
                    let matched = term.prepare_integral_specs(vakint, settings, None)?;
                    if !matches!(
                        vakint.evaluation_method(settings, &matched.canonical_topology),
                        Some(EvaluationMethod::Forcer(_))
                    ) {
                        return Err(super::invalid(
                            "massive closure of a p-integral requires the Forcer backend",
                        ));
                    }
                    let request = matched
                        .propagator_request(term.numerator.as_view())
                        .expect("a p-integral match has a propagator request");
                    preparation += preparation_started.elapsed();
                    cut_loops_max = cut_loops_max.max(request.family.loops);
                    let reduction_started = Instant::now();
                    let values = forcer::reduce_p(
                        request.family,
                        request.powers,
                        request.canonical_numerator,
                        &request.external,
                        settings,
                    )?;
                    reduction += reduction_started.elapsed();
                    let finalization_started = Instant::now();
                    let integrated =
                        self.finalize(&values, request.family.loops, &head, settings, mode)?;
                    finalization += finalization_started.elapsed();
                    let combination_started = Instant::now();
                    result += integrated;
                    combination += combination_started.elapsed();
                    propagator_terms += 1;
                }
            }
            if mode == RustRedMasterMode::Numerical {
                let numerical_started = Instant::now();
                result = vakint.resolve_rustred_numerical(result, settings)?;
                numerical += numerical_started.elapsed();
            }
            let combination_started = Instant::now();
            let result = result.expand();
            combination += combination_started.elapsed();
            tracing::info!(
                uv = true,
                rstar = true,
                timing = true,
                source_terms = source_terms,
                local_terms = local_terms,
                scaleless_terms = scaleless_terms,
                propagator_terms = propagator_terms,
                cut_loops_max = cut_loops_max,
                preparation_seconds = preparation.as_secs_f64(),
                forcer_reduction_seconds = reduction.as_secs_f64(),
                radial_master_laurent_seconds = finalization.as_secs_f64(),
                result_combination_seconds = combination.as_secs_f64(),
                numerical_resolution_seconds = numerical.as_secs_f64(),
                seconds = started.elapsed().as_secs_f64(),
                "Evaluated massive p-integral closure"
            );
            Ok(result)
        })
    }

    fn finalize(
        &self,
        values: &[ReducedPropagatorValue],
        cut_loops: usize,
        head: &Atom,
        settings: &crate::VakintSettings,
        mode: RustRedMasterMode,
    ) -> Result<Atom, VakintError> {
        let epsilon = vk_symbol!(&settings.epsilon_symbol);
        let square = function!(S.dot, head, head);
        let loops = cut_loops as i64 + 1;
        let target = settings.number_of_terms_in_epsilon_expansion - loops - 1;
        let mut result = Atom::Zero;
        for value in values {
            // The backend has already reduced a scalar numerator. Its PG(x,y)
            // factors are scalar master tags, not two-index vectors for the
            // generic tensor notation converter.
            let reduced = value.reduced.expand();
            for (monomial, coefficient) in
                reduced.coefficient_list::<i16>(std::slice::from_ref(&square))
            {
                let radial_power = if monomial == Atom::one() {
                    0
                } else if monomial == square {
                    1
                } else if let AtomView::Pow(power) = monomial.as_view()
                    && power.get_base_exp().0 == square.as_view()
                    && let Some(power) = crate::get_integer_from_atom(power.get_base_exp().1)
                {
                    power
                } else {
                    return Err(super::invalid(format!(
                        "the closing momentum needs integer powers of its square: {monomial}"
                    )));
                };
                let mut hard_dependence = coefficient.contains_symbol(S.k);
                coefficient.visitor(&mut |part| {
                    if let AtomView::Fun(function) = part
                        && function.get_symbol() == S.p
                        && let AtomView::Fun(external) = head.as_view()
                        && function.iter().next() == external.iter().next()
                    {
                        hard_dependence = true;
                    }
                    !hard_dependence
                });
                if hard_dependence {
                    return Err(VakintError::InvalidNumerator(format!(
                        "massive closure requires scalar dependence on the hard momentum: {coefficient}"
                    )));
                }
                super::compare::check_symbols(&coefficient)?;
                let scale_power = value.scale_power + radial_power;
                let sign = if (self.power + radial_power) % 2 == 0 {
                    1
                } else {
                    -1
                };
                let prefactor = self.radial_prefactor(settings, cut_loops, scale_power)?;
                let coefficient = Atom::num(sign) * coefficient;
                result += if mode == RustRedMasterMode::Retain {
                    simplify_real((prefactor * coefficient).as_view())
                } else {
                    laurent::expand(&prefactor, &coefficient, epsilon, target)?
                };
            }
        }
        if mode != RustRedMasterMode::Retain {
            crate::truncation_guard::ensure_known_orders(
                &result,
                epsilon,
                Some(target),
                "Forcer massive closure",
                crate::truncation_guard::plain_name,
            )?;
        }
        let mu = Atom::var(vk_symbol!(&settings.mu_r_sq_symbol));
        Ok(simplify_real(result.as_view())
            .replace(Atom::var(vk_symbol!("logMassiveClosure")).to_pattern())
            .with(function!(Symbol::LOG, &self.mass_squared).to_pattern())
            .replace(Atom::var(vk_symbol!("log_mu_sq")).to_pattern())
            .with(function!(Symbol::LOG, mu).to_pattern()))
    }

    fn radial_prefactor(
        &self,
        settings: &crate::VakintSettings,
        cut_loops: usize,
        scale_power: i64,
    ) -> Result<Atom, VakintError> {
        let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let loops = cut_loops as i64 + 1;
        let normalization = settings
            .get_staged_integral_normalization_factor_atom()?
            .replace(S.n_loops.to_pattern())
            .with(Atom::num(loops).to_pattern());
        let measure = vk_parse!(&format!(
            "(𝑖*𝜋^(2-({eps})))^{loops}\
             *(Gamma(1+{eps})*Gamma(1-{eps})^2/Gamma(2-2*{eps}))^{cut_loops}\
             *exp(-logMassiveClosure)^({eps}*{loops})",
            eps = settings.epsilon_symbol,
        ))
        .map_err(VakintError::SymbolicaError)?;
        // ∫ d^d q_E/π^(d/2) (q_E²)^(w-l ε)/(q_E²+M²)^a
        // = (M²)^(2+w-a-(l+1)ε) Γ(2+w-(l+1)ε)
        //   Γ(a-2-w+(l+1)ε)/(Γ(2-ε) Γ(a)).
        let gamma_ratio = gamma().call(Atom::num(2 + scale_power) - Atom::num(loops) * &epsilon)
            * gamma().call(Atom::num(self.power - 2 - scale_power) + Atom::num(loops) * &epsilon)
            / (gamma().call(Atom::num(2) - &epsilon) * gamma().call(Atom::num(self.power)));
        Ok(Vakint::normalize_gamma_functions(
            (normalization
                * measure
                * self
                    .mass_squared
                    .pow(Atom::num(2 + scale_power - self.power))
                * gamma_ratio)
                .as_view(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::LazyLock;

    use super::*;
    use crate::{EvaluationOrder, ForcerOptions, LoopNormalizationFactor, VakintSettings};

    static VAKINT: LazyLock<Vakint> = LazyLock::new(|| Vakint::new().unwrap());

    fn parse(text: &str) -> Atom {
        Vakint::initialize_vakint_symbols();
        vk_parse!(text).unwrap()
    }

    fn settings(terms: i64) -> VakintSettings {
        VakintSettings {
            evaluation_order: EvaluationOrder::empty(),
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            n_cores: 1,
            number_of_terms_in_epsilon_expansion: terms,
            use_dot_product_notation: true,
            ..Default::default()
        }
    }

    fn closure(power: i64) -> MassiveClosure {
        MassiveClosure {
            momentum: parse("p(1)"),
            mass_squared: Atom::one(),
            power,
        }
    }

    fn assert_equal(left: &Atom, right: &Atom) {
        assert!(
            (left - right).expand().is_zero(),
            "left: {left}\nright: {right}\ndifference: {}",
            (left - right).expand()
        );
    }

    #[test]
    fn massive_closure_preserves_tadpole_poles_finite_and_higher_orders() {
        let settings = settings(5);
        let value = closure(2)
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        let expected = parse("1/ε+𝜋^2*ε/12+𝜋^4*ε^3/160")
            - symbolica::transcendental::zeta().call(3) * parse("ε^2/3");
        // Both m² and μ² are one; no scale logarithms survive.
        let value = value
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::one().to_pattern());
        assert_equal(&value, &expected);
    }

    #[test]
    fn massive_closure_matches_the_independent_vacuum_backend() {
        let mut settings = settings(5);
        for power in [1, 2, 3, 4] {
            let value = closure(power)
                .evaluate(&VAKINT, &settings, Atom::one().as_view())
                .unwrap();
            settings.evaluation_order = EvaluationOrder(vec![EvaluationMethod::RustRed(
                crate::RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]);
            let reference = VAKINT
                .evaluate(&settings, parse(&format!("topo(I1L(1,{power}))")).as_view())
                .unwrap();
            settings.evaluation_order = EvaluationOrder::empty();
            assert_equal(&value, &reference);
        }
    }

    #[test]
    fn massive_closure_handles_mass_scaling_and_radial_numerators() {
        let settings = settings(4);
        let massive = MassiveClosure {
            mass_squared: parse("mass_squared"),
            ..closure(3)
        };
        // q²/(q²-m²)^3 = 1/(q²-m²)^2 + m²/(q²-m²)^3.
        let with_numerator = massive
            .evaluate(&VAKINT, &settings, parse("dot(p(1),p(1))").as_view())
            .unwrap();
        let lower_power = MassiveClosure {
            power: 2,
            ..massive.clone()
        }
        .evaluate(&VAKINT, &settings, Atom::one().as_view())
        .unwrap();
        let scalar = massive
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        assert_equal(
            &with_numerator,
            &(lower_power + &massive.mass_squared * scalar),
        );
    }

    #[test]
    fn massive_closure_keeps_the_massless_tree_denominator() {
        let settings = settings(3);
        let value = closure(1)
            .evaluate(&VAKINT, &settings, parse("1/dot(p(1),p(1))").as_view())
            .unwrap()
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::one().to_pattern());
        assert_equal(&value, &parse("1/ε+1+ε*(1+𝜋^2/12)"));
    }

    #[test]
    fn existing_tensor_prepass_projects_a_mixed_mass_vacuum_before_matching() {
        let settings = settings(3);
        let topology = parse("topo(prop(1,edge(1,2),k(1),1,1)*prop(2,edge(2,1),k(1),0,1))");
        let input = parse("dot(k(1),p(2))^2") * &topology;
        let projected = VAKINT.tensor_reduce(&settings, input.as_view()).unwrap();
        let expression = VakintExpression::try_from(projected).unwrap();
        assert_eq!(expression.0.len(), 1);
        assert_eq!(expression.0[0].integral, topology);
        assert!(
            (expression.0[0].numerator.clone() - parse("dot(k(1),k(1))*dot(p(2),p(2))/(4-2*ε)"))
                .together()
                .is_zero()
        );
    }

    #[test]
    fn massive_closure_rejects_open_hard_indices_and_unintegrated_loops() {
        for input in ["p(1,mu)", "dot(p(1),p(2))", "dot(k(1),k(1))"] {
            let error = closure(2)
                .evaluate(&VAKINT, &settings(3), parse(input).as_view())
                .unwrap_err();
            assert!(matches!(error, VakintError::InvalidNumerator(_)), "{error}");
        }
    }

    #[test]
    fn massive_closure_expands_the_exact_bubble_before_truncation() {
        let settings = settings(5);
        let value = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: parse("1/ε"),
                }],
                1,
                &parse("p(1)"),
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        // Direct successive one-loop integration, with both loop measures
        // e^(γε)/(iπ^(d/2)); Γ(2-2ε) cancels before the final expansion.
        let reference =
            parse("exp(2*EulerGamma*ε)*Gamma(1+ε)*Gamma(1-ε)^2*Gamma(2*ε)/(ε*Gamma(2-ε))");
        let reference =
            super::super::finalize::gamma_series(&reference, vk_symbol!("ε"), 2).unwrap();
        let value = value
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::one().to_pattern());
        assert_equal(&value, &reference);
    }

    #[test]
    fn massive_closure_retains_exact_gamma_dependence_when_requested() {
        let mut settings = settings(5);
        let expanded = closure(2)
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Retain,
            })]);
        let retained = closure(2)
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        assert!(retained.contains_symbol(gamma()));
        let restored = super::super::finalize::gamma_series(&retained, vk_symbol!("ε"), 3).unwrap();
        assert_equal(&restored, &expanded);
    }

    #[test]
    fn massive_closure_keeps_products_of_scalar_forcer_tags() {
        let settings = settings(5);
        let reduced = parse("PG(1,0)^2*PG(2,0)");
        let value = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced,
                }],
                2,
                &parse("p(1)"),
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        let gamma_value = super::super::finalize::pg_value(1, 0, vk_symbol!("ε")).pow(Atom::num(2))
            * super::super::finalize::pg_value(2, 0, vk_symbol!("ε"));
        let reference = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: gamma_value,
                }],
                2,
                &parse("p(1)"),
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        assert_equal(&value, &reference);
    }

    #[test]
    fn massive_closure_checks_master_depth_after_the_extra_loop_pole() {
        let error = closure(2)
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: parse("PMaster(no)"),
                }],
                3,
                &parse("p(1)"),
                &settings(13),
                RustRedMasterMode::Symbolic,
            )
            .unwrap_err();
        assert!(
            matches!(error, VakintError::MissingMasterOrders(_)),
            "{error}"
        );
    }

    #[test]
    fn closure_batch_matches_serial_with_masses_powers_poles_and_local_terms() {
        let mut settings = settings(5);
        if forcer::forcer_form_version(&settings.form_exe_path).is_err() {
            eprintln!("skipping live massive closure batch: FORM 5 is unavailable");
            return;
        }
        let directory = tempfile::tempdir().unwrap();
        settings.clean_tmp_dir = false;
        settings.temporary_directory = Some(directory.path().display().to_string());
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Symbolic,
            })]);
        let inputs = vec![
            (closure(2), parse("(1+ε)/ε*topo(P1L(p(1),1,1))+1/ε^2")),
            (
                MassiveClosure {
                    momentum: parse("p(3)"),
                    mass_squared: parse("mass_squared"),
                    power: 3,
                },
                parse("(2-ε)/ε^2*dot(k(1),p(3))^2*topo(P1L(p(3),2,1))"),
            ),
            (closure(1), parse("topo(P1L(p(1),1,0))")),
            (closure(4), parse("dot(p(1),p(1))/ε")),
            (closure(2), Atom::Zero),
        ];
        for master_mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            settings.evaluation_order =
                EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                    master_mode,
                })]);
            let actual = MassiveClosure::evaluate_batch(&VAKINT, &settings, &inputs).unwrap();
            for ((closure, input), actual) in inputs.iter().zip(actual) {
                let expected = closure
                    .evaluate(&VAKINT, &settings, input.as_view())
                    .unwrap();
                assert_equal(&actual, &expected);
            }
        }
    }

    #[test]
    fn massive_closure_runs_forcer_and_retains_local_terms() {
        let mut settings = settings(5);
        if forcer::forcer_form_version(&settings.form_exe_path).is_err() {
            eprintln!("skipping live massive closure: FORM 5 is unavailable");
            return;
        }
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
                master_mode: RustRedMasterMode::Symbolic,
            })]);
        let bubble = parse("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))");
        let massive = closure(2);
        let expected = massive
            .finalize(
                &[ReducedPropagatorValue {
                    scale_power: 0,
                    reduced: parse("1/ε"),
                }],
                1,
                &massive.momentum,
                &settings,
                RustRedMasterMode::Symbolic,
            )
            .unwrap();
        let local = massive
            .evaluate(&VAKINT, &settings, Atom::one().as_view())
            .unwrap();
        let value = massive
            .evaluate(&VAKINT, &settings, (&bubble + Atom::one()).as_view())
            .unwrap();
        assert_equal(&value, &(expected + local));
        let scaleless = parse("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,0))");
        assert!(
            massive
                .evaluate(&VAKINT, &settings, scaleless.as_view())
                .unwrap()
                .is_zero()
        );
    }
}
