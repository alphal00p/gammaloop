//! Measure individual public stages on saved contraction-benchmark inputs.
//!
//! First run `metric_contraction_benchmark` to create input snapshots, then run
//! `cargo run -p idenso --profile dev-optim --example contraction_phase_benchmark
//! -- INPUT_DIRECTORY OUTPUT_DIRECTORY [samples=5] [target_batch_ms=8]`.
//!
//! Reads `.input` fixtures and `.gamma_input` original traces. Original and
//! already-simplified gamma expressions are reported separately: their execution
//! paths differ, and overlapping stage timings cannot be added to reconstruct a
//! full gamma call. Parsing, validation and output snapshots are outside timing;
//! timed calls include output destruction. Zero samples validates and snapshots
//! outputs without measuring. Symbols used by the source benchmark are registered
//! before parsing so tagged rank-one tensors retain their contraction semantics.
//! Three additional factored generic-tensor controls exercise a selected bispinor
//! representation and both orientations of an unrelated dual representation.
//!
//! Each method is called directly through the linked production API. Compare
//! output snapshots across builds; term counts alone do not establish equivalence.

use std::{hint::black_box, path::Path, time::Instant};

use idenso::{
    dirac::GammaSimplifier,
    epsilon::EpsilonSimplifier,
    representations::Bispinor,
    shorthands::{
        chain::Chain,
        schoonschip::{Schoonschip, SchoonschipSettings},
    },
};
use spenso::network::{parsing::AtomStructureExt, tags::SPENSO_TAG};
use symbolica::atom::{Atom, AtomCore, AtomView};

#[derive(Clone, Copy, Debug)]
enum Method {
    Clone,
    RepeatedIndices,
    NormalizeDots,
    MetricSchoonschip,
    Schoonschip,
    Epsilon,
    Chainify,
    CollectChains,
    CollectGammaChains,
    Gamma,
}
const METHODS: [Method; 10] = [
    Method::Clone,
    Method::RepeatedIndices,
    Method::NormalizeDots,
    Method::MetricSchoonschip,
    Method::Schoonschip,
    Method::Epsilon,
    Method::Chainify,
    Method::CollectChains,
    Method::CollectGammaChains,
    Method::Gamma,
];

enum Outcome {
    Atom(Atom),
    Boolean(bool),
}
impl Outcome {
    fn atom(self) -> Atom {
        match self {
            Self::Atom(atom) => atom,
            Self::Boolean(value) => Atom::num(i64::from(value)),
        }
    }
}

fn run(method: Method, view: AtomView<'_>) -> Outcome {
    use Method::*;
    if matches!(method, RepeatedIndices) {
        return Outcome::Boolean(view.has_repeated_explicit_indices());
    }
    Outcome::Atom(match method {
        Clone => view.to_owned(),
        NormalizeDots => view.normalize_dots(),
        MetricSchoonschip => view.schoonschip_with_settings(
            &SchoonschipSettings::default()
                .without_rank1_tensors()
                .with_chain_like_functions(),
        ),
        Schoonschip => view
            .schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions()),
        Epsilon => view.simplify_epsilon(),
        Chainify => view.chainify(Bispinor {}.into()),
        CollectChains => view.collect_chains(Bispinor {}.into()),
        CollectGammaChains => view.collect_gamma_chains(),
        Gamma => view.simplify_gamma(),
        RepeatedIndices => unreachable!(),
    })
}
fn batch(method: Method, input: &Atom, repetitions: usize) -> f64 {
    let start = Instant::now();
    for _ in 0..repetitions {
        drop(black_box(run(method, black_box(input.as_view()))));
    }
    start.elapsed().as_nanos() as f64 / repetitions as f64
}

fn main() {
    let arguments: Vec<_> = std::env::args().collect();
    let input_directory = Path::new(&arguments[1]);
    let output_directory = Path::new(&arguments[2]);
    let samples: usize = arguments.get(3).map_or(5, |s| s.parse().unwrap());
    let target_ms: f64 = arguments.get(4).map_or(8., |s| s.parse().unwrap());
    assert!(target_ms.is_finite() && target_ms > 0.);
    std::fs::create_dir_all(output_directory).unwrap();
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
    for name in ["cleanup_benchmark::V", "cleanup_benchmark::W"] {
        SPENSO_TAG.rank_one_tensor_symbol(name);
    }
    for i in 0..10 {
        SPENSO_TAG.rank_one_tensor_symbol(&format!("metric_benchmark::p{i}"));
    }
    let mut paths: Vec<_> = std::fs::read_dir(input_directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            matches!(
                path.extension().and_then(|s| s.to_str()),
                Some("input" | "gamma_input")
            )
        })
        .collect();
    paths.sort();
    let inputs = paths.into_iter().map(|path| {
        let case = path.file_stem().unwrap().to_str().unwrap().to_owned();
        let original_gamma = path.extension().unwrap() == "gamma_input";
        let has_gamma = input_directory.join(format!("{case}.gamma_input")).exists();
        let variant = if original_gamma {
            "gamma_source"
        } else if has_gamma {
            "gamma_terminal"
        } else {
            "fixture"
        };
        let source = std::fs::read_to_string(path).unwrap();
        (case, variant, has_gamma, source)
    });
    let controls = [
        (
            "generic_bispinor_factored",
            "((x+y)^8+generic((x+y)^8,bis(n^2-1,label(a)),metadata,bis(n^2-1,label(b))))^2",
        ),
        (
            "generic_other_representation_factored",
            "((x+y)^8+generic((x+y)^8,cof(n^2-1,label(a)),metadata,dind(cof(n^2-1,label(b)))))^2",
        ),
        (
            "generic_other_representation_reversed",
            "((x+y)^8+generic((x+y)^8,dind(cof(n^2-1,label(a))),metadata,cof(n^2-1,label(b))))^2",
        ),
    ]
    .into_iter()
    .map(|(name, source)| (name.to_owned(), "fixture", false, source.to_owned()));
    let mut cases = Vec::new();
    for (case, variant, has_gamma, source) in inputs.chain(controls) {
        let input = Atom::parse(
            &source,
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        let mut nodes = 0;
        input.as_view().visitor(&mut |_| {
            nodes += 1;
            true
        });
        std::fs::write(
            output_directory.join(format!("{case}.{variant}.input")),
            input.to_plain_string(),
        )
        .unwrap();
        for method in METHODS {
            if matches!(method, Method::Gamma) && !has_gamma {
                continue;
            }
            let output = run(method, input.as_view()).atom();
            if matches!(method, Method::Gamma) && variant == "gamma_terminal" {
                assert_eq!(output, input, "{case}: terminal gamma expression changed");
            }
            let changed = output != input;
            let output_terms = output.nterms();
            std::fs::write(
                output_directory.join(format!("{case}.{variant}.{method:?}.output")),
                output.to_plain_string(),
            )
            .unwrap();
            cases.push((
                case.clone(),
                variant,
                input.clone(),
                method,
                nodes,
                changed,
                output_terms,
            ));
        }
    }
    if samples == 0 {
        eprintln!("Validated {} case/method outputs", cases.len());
        return;
    }
    let repetitions: Vec<_> = cases
        .iter()
        .map(|(_, _, input, method, _, _, _)| {
            let estimate = batch(*method, input, 1).max(1.);
            ((target_ms * 1e6 / estimate).ceil() as usize).clamp(1, 100_000)
        })
        .collect();
    let mut timings = vec![Vec::new(); cases.len()];
    for sample in 0..samples {
        for offset in 0..cases.len() {
            let index = (sample + offset) % cases.len();
            let (_, _, input, method, _, _, _) = &cases[index];
            timings[index].push(batch(*method, input, repetitions[index]));
        }
    }
    for (index, (case, variant, input, method, nodes, changed, output_terms)) in
        cases.iter().enumerate()
    {
        let mut sorted = timings[index].clone();
        sorted.sort_by(f64::total_cmp);
        let median_ns = if samples.is_multiple_of(2) {
            (sorted[samples / 2 - 1] + sorted[samples / 2]) / 2.
        } else {
            sorted[samples / 2]
        };
        println!(
            "{{\"case\":{case:?},\"input_variant\":{variant:?},\"method\":\"{method:?}\",\"nodes\":{nodes},\"terms\":{},\"output_terms\":{output_terms},\"changed\":{changed},\"repetitions\":{},\"samples_ns\":{:?},\"median_ns\":{median_ns}}}",
            input.nterms(),
            repetitions[index],
            timings[index]
        );
    }
}
