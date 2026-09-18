//! Opt-in finite-target experiments, not four-loop production acceptance.
//!
//! The ignored offline oracle is deliberately separate from the native scalar
//! lane. Its input must be generated from the actual candidate terminal keys;
//! it neither decides which leaves are terminals nor supplies reduction rules.

#![cfg(feature = "experimental-rustred")]

#[path = "experimental_rustred_4l/input.rs"]
mod input;
mod test_utils;
#[path = "experimental_rustred_4l/timing.rs"]
mod timing;

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use rustred::family::IntegralKey;
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::rustred_evaluation::experimental::native::NativeCandidate;
use vakint::rustred_evaluation::experimental::{
    CandidateReduction, CandidateScalarReduction, ExperimentalRustRed, prepare_candidate_integrals,
};
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, Vakint, VakintSettings,
};
use vakint::{vakint_parse as vk_parse, vakint_symbol as vk_symbol};

/// The complete four-loop analytic/FMFT acceptance inventory currently used
/// by Vakint.  This is intentionally kept as an explicit list: a future
/// sealed four-loop RustRed artifact must exercise every one of these inputs,
/// not just the finite parent/dotted/pinch probes below.
const FOUR_LOOP_ANALYTIC_FMFT_CASES: [&str; 14] = [
    "test_integrate_4l_h",
    "test_integrate_4l_h_squared_mass",
    "test_integrate_4l_h_rank_4",
    "test_integrate_4l_h_rank_4_additional_symbols_numerator",
    "test_integrate_4l_PR9d_from_H",
    "test_integrate_4l_PR9d_from_X",
    "test_integrate_4l_PR9d_from_H_pinch",
    "test_integrate_4l_PR9d_from_FG",
    "test_integrate_4l_PR9d_from_FG_pinch",
    "test_integrate_4l_PR11d",
    "test_integrate_4l_clover",
    "test_integrate_4l_clover_with_non_unit_scales",
    "test_integrate_4l_dotted_clover",
    "test_integrate_4l_clover_with_numerator",
];

/// Inventory-only gate for the eventual complete four-loop lane.
///
/// This deliberately remains ignored until a sealed four-loop artifact and
/// terminal catalog are shipped.  It checks that the upstream FMFT cases are
/// still present and reports the precise production prerequisite instead of
/// pretending that finite candidate experiments constitute whole-family
/// acceptance.
#[test]
#[ignore = "requires a sealed four-loop RustRed artifact and terminal catalog"]
fn four_loop_rustred_fmft_inventory_requires_sealed_artifact() {
    let source = include_str!("integral_evaluation_analytic_tests.rs");
    for case in FOUR_LOOP_ANALYTIC_FMFT_CASES {
        assert!(
            source.contains(&format!("fn {case}()")),
            "upstream four-loop FMFT case disappeared: {case}"
        );
    }
    panic!(
        "four-loop RustRed lane is not enabled: {} FMFT cases are inventoried, \
         but no authenticated artifact/catalog is registered in Vakint",
        FOUR_LOOP_ANALYTIC_FMFT_CASES.len()
    );
}

/// A negative-boundary stub only; never a numerical reduction oracle.
#[derive(Debug)]
struct RejectingReducer {
    momenta: Vec<Atom>,
    dimension: Atom,
}

impl CandidateScalarReduction for RejectingReducer {
    fn index_count(&self) -> usize {
        10
    }
    fn parent_momenta(&self) -> &[Atom] {
        &self.momenta
    }
    fn dimension(&self) -> &Atom {
        &self.dimension
    }
    fn reduce_unit_mass(&self, _target: &IntegralKey) -> Result<CandidateReduction, String> {
        Err("sentinel: no candidate rules supplied".into())
    }
}

#[test]
fn experimental_catalog_rejects_wrong_arity_without_claiming_reduction() {
    Vakint::initialize_vakint_symbols();
    let reducer = Arc::new(RejectingReducer {
        momenta: vec![vk_parse!("k(1)").unwrap()],
        dimension: vk_parse!("candidate_dimension").unwrap(),
    });
    let terminals = BTreeMap::from([(IntegralKey::try_new([1; 9]).unwrap(), Atom::num(1))]);
    assert!(ExperimentalRustRed::new(reducer, terminals).is_err());
}

#[test]
fn experimental_catalog_requires_exact_coefficients_before_normalization() {
    Vakint::initialize_vakint_symbols();
    let reducer = Arc::new(RejectingReducer {
        momenta: vec![vk_parse!("k(1)").unwrap()],
        dimension: vk_parse!("candidate_dimension").unwrap(),
    });
    let key = IntegralKey::try_new([1; 10]).unwrap();
    let exact = BTreeMap::from([(key.clone(), vk_parse!("5/4*PR4").unwrap())]);
    assert!(ExperimentalRustRed::new(reducer.clone(), exact).is_ok());
    for expression in ["1.25`40*PR4", "f(1+1.25`40*x)"] {
        let approximate = BTreeMap::from([(key.clone(), vk_parse!(expression).unwrap())]);
        let error = ExperimentalRustRed::new(reducer.clone(), approximate).unwrap_err();
        assert!(error.to_string().contains("approximate coefficient"));
    }
}

/// Input lines: `ten comma-separated powers<TAB>full common-mass Vakint input`.
///
/// The caller constructs the second field with the retained parent routing and
/// auxiliary-denominator numerator for the first field, using symbolic
/// `muvsq` until after evaluation, when both mass scales are set to one.
/// This is an explicit
/// offline trust boundary: output records are oracle values, not certificates.
/// No guessed PR assignment, topology-name switch or missing-master promotion
/// occurs in the experimental native evaluator.
#[test]
#[ignore = "offline FMFT oracle: requires actual candidate terminal inputs and FORM"]
fn offline_fmft_candidate_terminal_catalog() {
    let source = std::env::var("VAKINT_4L_CANDIDATE_TERMINAL_INPUTS")
        .expect("supply candidate-derived terminal inputs, not guessed masters");
    let form = std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH")
        .expect("offline oracle requires an explicit FORM executable");
    let input = std::fs::read_to_string(&source).expect("read terminal input file");
    let settings = VakintSettings {
        form_exe_path: form,
        epsilon_symbol: "ep".into(),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {
            expand_masters: false,
            susbstitute_masters: false,
        })]),
        ..VakintSettings::default()
    };
    let vakint = Vakint::new().unwrap();
    let mut seen = BTreeSet::new();
    for (line_number, line) in input.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (key, expression) = line.split_once('\t').unwrap_or_else(|| {
            panic!(
                "line {} requires powers and input separated by TAB",
                line_number + 1
            )
        });
        let powers = key
            .split(',')
            .map(|power| power.trim().parse::<i64>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(
            powers.len(),
            10,
            "four-loop complete family requires ten powers"
        );
        assert!(seen.insert(powers), "duplicate terminal key");
        assert!(
            seen.len() <= 512,
            "offline terminal count exceeds explicit bound"
        );
        let expression = vk_parse!(expression).unwrap();
        let result = vakint
            .evaluate_integral(&settings, expression.as_view())
            .unwrap();
        let result = result
            .replace(vk_parse!("mursq").unwrap().to_pattern())
            .with(Atom::num(1).to_pattern())
            .replace(vk_parse!("muvsq").unwrap().to_pattern())
            .with(Atom::num(1).to_pattern())
            .together();
        assert!(
            !result.contains_symbol(vakint::symbols::S.topo),
            "unreduced integral"
        );
        assert!(
            !result.contains_symbol(vk_symbol!("muvsq")),
            "input must use unit mass"
        );
        println!("{key}\t{}", result.to_canonical_string());
    }
    assert!(!seen.is_empty(), "no candidate terminal inputs supplied");
}

/// A finite-target experiment with three deliberately separate phases:
/// RustRed candidate search/application, offline terminal evaluation, then a
/// frozen-catalog native scalar lane with an invalid FORM executable path.
/// It does not prove arbitrary-index or whole-family completeness.
#[test]
#[ignore = "experimental: generates actual four-loop candidates and uses an explicit offline FMFT oracle"]
fn candidate_parent_dotted_and_pinch_match_fmft() {
    let source = std::env::var("VAKINT_4L_CANDIDATE_PARENT_INPUT")
        .map(|path| std::fs::read_to_string(path).unwrap())
        .unwrap_or_else(|_| include_str!("inputs/experimental_four_loop_h.csv").into());
    run_candidate_finite_family("H", source);
}

/// Run the same finite-target candidate/FeynKit/FMFT comparison for every
/// supplied complete four-loop parent descriptor.  This remains an ignored,
/// offline experiment: each family gets its own RustRed candidate reducer and
/// terminal catalog, and no result is registered in the production artifact
/// registry.  Keeping the matrix here reuses exactly the same harness as the
/// historical H probe instead of silently testing a weaker direct evaluator.
#[test]
#[ignore = "experimental: generates four-loop candidates and uses an explicit offline FMFT oracle"]
fn candidate_all_four_loop_parents_match_fmft() {
    for (name, source) in [
        ("H", include_str!("inputs/experimental_four_loop_h.csv")),
        ("FG", include_str!("inputs/experimental_four_loop_fg.csv")),
        ("BMW", include_str!("inputs/experimental_four_loop_bmw.csv")),
        ("X", include_str!("inputs/experimental_four_loop_x.csv")),
    ] {
        run_candidate_finite_family(name, source.to_owned());
    }
}

fn run_candidate_finite_family(family_name: &'static str, source: String) {
    test_utils::run_multi_lane_acceptance(move || {
        Vakint::initialize_vakint_symbols();
        let form = std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH")
            .expect("explicit offline FMFT oracle executable required");
        let workers = std::env::var("VAKINT_4L_CANDIDATE_WORKERS")
            .map_or(1, |value| value.parse::<usize>().unwrap());
        assert!(
            (1..=6).contains(&workers),
            "experimental worker budget is 1..=6"
        );
        let dot_power = std::env::var("VAKINT_4L_CANDIDATE_DOT_POWER")
            .map_or(2, |value| value.parse::<i64>().expect("integer dot power"));
        assert!(dot_power >= 2, "dotted target power must be at least two");
        let parent = input::ParentInput::from_csv(&source);
        let started = std::time::Instant::now();
        let native = Arc::new(
            NativeCandidate::<10>::solve(
                &parent.family,
                parent.physical_momenta.clone(),
                workers,
                None,
            )
            .expect("RustRed candidate generation and pointwise bridge"),
        );
        println!(
            "candidate {family_name} search: {:?}; fixed residuals: {}",
            started.elapsed(),
            native.terminals().len()
        );
        let physical = parent.physical_momenta.len();
        let corner: Vec<i64> = (0..10).map(|axis| i64::from(axis < physical)).collect();
        let mut dotted = corner.clone();
        dotted[0] = dot_power;
        let mut pinch = corner.clone();
        pinch[2] = 0;
        let cases = [("parent", corner), ("dotted", dotted), ("pinch", pinch)];
        let vakint = Vakint::new().unwrap();
        let preparation_settings = VakintSettings {
            form_exe_path: "/candidate-preparation-must-not-use-form".into(),
            use_dot_product_notation: true,
            ..VakintSettings::default()
        };
        let mut reached = BTreeSet::new();
        let mut applications = 0usize;
        for (name, powers) in &cases {
            let input = parent.integral(powers);
            let canonical = vakint
                .to_canonical(&preparation_settings, input.as_view(), false)
                .unwrap();
            let scalar = vakint
                .tensor_reduce(&preparation_settings, canonical.as_view())
                .unwrap();
            let prepared = prepare_candidate_integrals(
                &vakint,
                &preparation_settings,
                scalar.as_view(),
                native.as_ref(),
            )
            .unwrap();
            for target in prepared {
                let result = native
                    .reduce_unit_mass(&target.target)
                    .unwrap_or_else(|error| panic!("candidate target {name} unresolved: {error}"));
                applications += result.applied_rules;
                println!(
                    "candidate {name} {:?}: {} new rule applications; {} reached terminals; declared-terminal target={}",
                    target.target.powers(),
                    result.applied_rules,
                    result.terms.len(),
                    native.terminals().contains(&target.target)
                );
                for (terminal, _) in result.terms {
                    assert!(
                        native.terminals().contains(&terminal),
                        "no undeclared master promotion"
                    );
                    reached.insert(terminal);
                }
            }
        }
        assert!(
            applications > 0,
            "terminal-only comparisons do not test IBP application"
        );
        assert!(
            reached.len() <= 512,
            "explicit offline terminal budget exceeded"
        );

        // The only FORM use in this experiment is independent oracle setup and
        // the legacy comparison lane. No oracle output becomes a candidate rule.
        let oracle_settings = VakintSettings {
            form_exe_path: form.clone(),
            epsilon_symbol: "ep".into(),
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions {
                expand_masters: false,
                susbstitute_masters: false,
            })]),
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            use_dot_product_notation: true,
            ..VakintSettings::default()
        };
        let mut catalog = BTreeMap::new();
        for terminal in reached {
            let oracle_input = parent.integral(terminal.powers());
            let value = vakint
                .evaluate(&oracle_settings, oracle_input.as_view())
                .unwrap_or_else(|error| {
                    panic!("offline FMFT terminal {:?}: {error}", terminal.powers())
                })
                .replace(vk_parse!("mursq").unwrap().to_pattern())
                .with(Atom::num(1).to_pattern())
                .replace(vk_parse!("muvsq").unwrap().to_pattern())
                .with(Atom::num(1).to_pattern())
                .together();
            println!(
                "offline terminal {:?}\t{}",
                terminal.powers(),
                value.to_canonical_string()
            );
            catalog.insert(terminal, value);
        }
        if std::env::var_os("VAKINT_ACCEPTANCE_EXACT_DIAGNOSTICS").is_some() {
            let options = FMFTOptions {
                expand_masters: true,
                susbstitute_masters: false,
            };
            let symbolic_evaluator = ExperimentalRustRed::new(native.clone(), catalog.clone())
                .unwrap()
                .with_master_options(options.clone());
            let mut symbolic_settings = oracle_settings.clone();
            symbolic_settings.number_of_terms_in_epsilon_expansion = 5;
            symbolic_settings.run_time_decimal_precision = 25;
            symbolic_settings.integral_normalization_factor = LoopNormalizationFactor::MSbar;
            symbolic_settings.evaluation_order =
                EvaluationOrder(vec![EvaluationMethod::FMFT(options)]);
            for (name, powers) in &cases {
                let input = parent.integral(powers);
                let oracle_symbolic = vakint
                    .evaluate(&symbolic_settings, input.as_view())
                    .unwrap();
                test_utils::log_laurent_diagnostics(
                    &format!("before numerical masters, FMFT {name}"),
                    oracle_symbolic.as_view(),
                    &symbolic_settings.epsilon_symbol,
                );
                let canonical = vakint
                    .to_canonical(&preparation_settings, input.as_view(), false)
                    .unwrap();
                let scalar = vakint
                    .tensor_reduce(&preparation_settings, canonical.as_view())
                    .unwrap();
                let mut native_settings = symbolic_settings.clone();
                native_settings.form_exe_path =
                    "/symbolic-candidate-diagnostics-must-not-use-form".into();
                let native_symbolic = symbolic_evaluator
                    .evaluate_integral(&vakint, &native_settings, scalar.as_view())
                    .unwrap();
                test_utils::log_laurent_diagnostics(
                    &format!("before numerical masters, native {name}"),
                    native_symbolic.value.as_view(),
                    &native_settings.epsilon_symbol,
                );
                let unit_mass = |value: &Atom| {
                    let mut approximate = false;
                    value.visitor(&mut |part| {
                        if let AtomView::Num(number) = part {
                            approximate |= number.get_coeff_view().is_float();
                        }
                        !approximate
                    });
                    assert!(!approximate, "exact diagnostic requires exact coefficients");
                    value
                        .replace(vk_parse!("mursq").unwrap().to_pattern())
                        .with(Atom::num(1).to_pattern())
                        .replace(vk_parse!("muvsq").unwrap().to_pattern())
                        .with(Atom::num(1).to_pattern())
                        .expand()
                        .together()
                };
                let oracle_exact = unit_mass(&oracle_symbolic);
                let native_exact = unit_mass(&native_symbolic.value);
                let epsilon = vk_parse!(native_settings.epsilon_symbol.as_str()).unwrap();
                for (lane, value) in [("FMFT", &oracle_exact), ("native", &native_exact)] {
                    let coefficients = value.coefficient_list::<i8>(&[epsilon.clone()]);
                    for exponent in -4..0 {
                        let monomial = epsilon.clone().pow(Atom::num(exponent));
                        let coefficient = coefficients
                            .iter()
                            .find(|(power, _)| *power == monomial)
                            .map_or(Atom::Zero, |(_, coefficient)| coefficient.expand());
                        eprintln!(
                            "EXACT_PRE_MASTER [{name}, {lane}] epsilon^{exponent}: zero={}; coefficient={}",
                            coefficient.is_zero(),
                            coefficient.to_canonical_string()
                        );
                    }
                }
                let difference = (oracle_exact - native_exact).expand().together();
                eprintln!(
                    "EXACT_PRE_MASTER [{name}] unit-mass difference zero={}; expression={}",
                    difference.is_zero(),
                    difference.to_canonical_string()
                );
            }
        }
        let evaluator = Arc::new(ExperimentalRustRed::new(native.clone(), catalog).unwrap());
        let lanes = [
            test_utils::EvaluationTestLane::oracle(
                EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
                test_utils::EvaluationTestInput::TensorReduced,
            ),
            test_utils::EvaluationTestLane::experimental_rustred_scalar(
                test_utils::EvaluationTestInput::TensorReduced,
                evaluator.clone(),
            ),
        ];
        let settings = VakintSettings {
            form_exe_path: form,
            run_time_decimal_precision: 25,
            number_of_terms_in_epsilon_expansion: 5,
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            ..VakintSettings::default()
        };
        let masses = vakint.params_from_f64(
            &settings,
            &std::collections::HashMap::from([("mursq".into(), 1.0), ("muvsq".into(), 1.0)]),
        );
        // Terminal discovery above necessarily populated the point cache.
        // Empty it so numerical parity also exercises actual recurrence
        // application under the harness's invalid scalar FORM path, not only
        // retrieval of values computed before the offline oracle ran.
        native.clear_cache().unwrap();
        let native_applications_before = native.rule_applications().unwrap();
        for (name, powers) in &cases {
            let integral = parent.integral(&powers);
            test_utils::compare_evaluations(
                settings.clone(),
                &lanes,
                test_utils::TensorPrepass::FeynKit,
                integral.as_view(),
                masses.clone(),
                std::collections::HashMap::default(),
                1e-20,
                10.0,
                true,
            );
            println!("finite-target numerical parity PASS: {family_name}/{name}; no family-closure claim");
        }
        let native_applications = native
            .rule_applications()
            .unwrap()
            .checked_sub(native_applications_before)
            .expect("candidate rule counter must not decrease");
        assert!(
            native_applications > 0,
            "invalid-FORM scalar comparisons must include cold candidate-rule application"
        );
        println!("invalid-FORM cold scalar lane: {native_applications} new rule applications");
        if let Ok(repeats) = std::env::var("VAKINT_4L_CANDIDATE_BENCH_REPEATS") {
            let inputs = cases
                .iter()
                .map(|(name, powers)| (*name, parent.integral(powers)))
                .collect::<Vec<_>>();
            timing::run(
                &vakint,
                &settings,
                native.as_ref(),
                evaluator.as_ref(),
                &inputs,
                repeats.parse().expect("integer benchmark repeat count"),
            );
        }
    });
}
