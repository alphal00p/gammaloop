//! Independent numerator-cancellation and graph-contraction comparisons.
//! These finite tests do not certify arbitrary-index closure.

use super::acceptance_inputs::ParentDescriptor;
use super::input::{NumeratorPropagatorPinch, ParentInput};
use super::test_utils::{self, EvaluationTestInput, EvaluationTestLane, TensorPrepass};
use super::{
    ExperimentalRustRed, NativeCandidate, OfflineTerminalCatalog, retained_parent_descriptor,
};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;
use symbolica::atom::{AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, Vakint, VakintSettings,
};

const PARENTS: [ParentDescriptor; 4] = [
    ParentDescriptor::H,
    ParentDescriptor::Fg,
    ParentDescriptor::Bmw,
    ParentDescriptor::X,
];

fn label(parent: ParentDescriptor) -> &'static str {
    match parent {
        ParentDescriptor::H => "H",
        ParentDescriptor::Fg => "FG",
        ParentDescriptor::Bmw => "BMW",
        ParentDescriptor::X => "X",
    }
}

fn cases(parent: ParentDescriptor) -> Vec<NumeratorPropagatorPinch> {
    let input = ParentInput::from_csv(parent.csv());
    // Each double contraction is a forest in all four input graphs. In
    // particular, do not cancel both FG edges 5 and 7: they carry its only
    // fourth-loop denominators and would give a vacuous scaleless zero.
    [
        ("D1", "P1", &[0][..]),
        ("D7", "P7", &[6][..]),
        ("D1D2", "P1P2", &[0, 1][..]),
        ("D5D6", "P5P6", &[4, 5][..]),
    ]
    .into_iter()
    .map(|(numerator, pinch, slots)| input.numerator_equals_propagators(slots, numerator, pinch))
    .collect()
}

#[test]
fn sixteen_inputs_are_expanded_and_have_registered_four_loop_witnesses() {
    Vakint::initialize_vakint_symbols();
    let mut count = 0;
    for parent in PARENTS {
        for pair in cases(parent) {
            // Both inputs go independently through the existing matcher.
            let first = retained_parent_descriptor(&pair.with_numerator);
            let second = retained_parent_descriptor(&pair.pinched);
            assert_eq!(first, parent);
            pair.with_numerator.visitor(&mut |node| {
                if let AtomView::Fun(fun) = node {
                    if fun.get_symbol() == vakint::symbols::S.dot {
                        for argument in fun.iter() {
                            assert!(matches!(argument, AtomView::Fun(f) if f.get_symbol() == vakint::symbols::S.k), "unexpanded scalar product: {node}");
                        }
                    }
                }
                true
            });
            println!(
                "{}/{} -> {} parent",
                label(parent),
                pair.pinch_name,
                label(second)
            );
            count += 1;
        }
    }
    assert_eq!(count, 16);
}

struct CandidateLane {
    native: Arc<NativeCandidate<10>>,
    evaluator: Arc<ExperimentalRustRed>,
}

impl CandidateLane {
    fn load(parent: ParentDescriptor, workers: usize) -> Self {
        let input = ParentInput::from_csv(parent.csv());
        let start = std::time::Instant::now();
        let native = Arc::new(
            NativeCandidate::solve(input.family.clone(), input.physical_momenta, workers, None)
                .unwrap(),
        );
        println!(
            "pinch candidate {}: {:?}; {} terminals",
            label(parent),
            start.elapsed(),
            native.terminals().len()
        );
        let directory = std::env::var_os("VAKINT_4L_CANDIDATE_CATALOG_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/inputs/experimental_four_loop_catalogs")
            });
        let path = directory.join(format!("{}.rrcat", label(parent).to_ascii_lowercase()));
        let catalog = OfflineTerminalCatalog::decode(
            &std::fs::read_to_string(&path).unwrap(),
            input.family.fingerprint(),
            10,
        )
        .unwrap();
        catalog.require_complete().unwrap();
        assert!(
            native
                .terminals()
                .iter()
                .all(|key| catalog.terms().contains_key(key)),
            "incomplete catalog {}",
            path.display()
        );
        let evaluator =
            Arc::new(ExperimentalRustRed::new(native.clone(), catalog.into_terms()).unwrap());
        Self { native, evaluator }
    }

    fn lane(&self) -> EvaluationTestLane {
        EvaluationTestLane::experimental_rustred_scalar(
            EvaluationTestInput::TensorReduced,
            self.evaluator.clone(),
        )
    }
}

#[test]
#[ignore = "release four-loop candidate generation and independent FMFT oracle"]
fn sixteen_numerator_propagator_pinches_match_fmft_and_rustred() {
    test_utils::run_multi_lane_acceptance(|| {
        Vakint::initialize_vakint_symbols();
        let form = std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH")
            .expect("explicit FMFT oracle executable required");
        let workers =
            std::env::var("VAKINT_4L_CANDIDATE_WORKERS").map_or(1, |v| v.parse::<usize>().unwrap());
        assert!((1..=6).contains(&workers));
        let filter = std::env::var("VAKINT_4L_CANDIDATE_FAMILY_FILTER").ok();
        assert!(
            filter
                .as_deref()
                .is_none_or(|f| PARENTS.iter().any(|p| label(*p) == f)),
            "unknown parent filter"
        );
        let mut candidates = BTreeMap::new();
        let mut passed = 0;
        for parent in PARENTS {
            if filter.as_deref().is_some_and(|f| label(parent) != f) {
                continue;
            }
            for pair in cases(parent) {
                let first = retained_parent_descriptor(&pair.with_numerator);
                let second = retained_parent_descriptor(&pair.pinched);
                for owner in [first, second] {
                    candidates
                        .entry(owner)
                        .or_insert_with(|| CandidateLane::load(owner, workers));
                }
                let first_native = &candidates[&first];
                let second_native = &candidates[&second];
                first_native.native.clear_cache().unwrap();
                second_native.native.clear_cache().unwrap();
                let before_first = first_native.native.rule_applications().unwrap();
                let before_second = second_native.native.rule_applications().unwrap();
                let oracle = EvaluationTestLane::oracle(
                    EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
                    EvaluationTestInput::TensorReduced,
                );
                let settings = VakintSettings {
                    form_exe_path: form.clone(),
                    run_time_decimal_precision: 32,
                    number_of_terms_in_epsilon_expansion: 5,
                    integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
                    ..VakintSettings::default()
                };
                // Non-unit scales make a dropped mass term/power observable.
                let masses = vakint::params_from_f64(
                    &[("muvsq".into(), 3.0), ("mursq".into(), 7.0)]
                        .into_iter()
                        .collect(),
                    32,
                );
                println!(
                    "pinch START {}/{} == {} (owners {}/{})",
                    label(parent),
                    pair.numerator_name,
                    pair.pinch_name,
                    label(first),
                    label(second)
                );
                test_utils::compare_input_pair(
                    settings,
                    &[
                        (oracle.clone(), oracle),
                        (first_native.lane(), second_native.lane()),
                    ],
                    TensorPrepass::FeynKit,
                    pair.with_numerator.as_view(),
                    pair.pinched.as_view(),
                    masses,
                    Default::default(),
                    1.0e-20,
                    0.0,
                );
                let mut applications =
                    first_native.native.rule_applications().unwrap() - before_first;
                if first != second {
                    applications +=
                        second_native.native.rule_applications().unwrap() - before_second;
                }
                println!(
                    "pinch PASS {}/{} == {}; {applications} applications",
                    label(parent),
                    pair.numerator_name,
                    pair.pinch_name
                );
                passed += 1;
            }
        }
        assert_eq!(passed, if filter.is_some() { 4 } else { 16 });
    });
}
