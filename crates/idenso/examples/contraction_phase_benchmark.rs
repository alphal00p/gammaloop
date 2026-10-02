//! Measure public operations and isolated kernels on saved contraction inputs.
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
//! Typed gamma methods include admission in each call.
//! Public operations call the linked production API. `Chainify` isolates the
//! internal notation primitive through the `reference-cases` feature. Compare
//! output snapshots across builds; term counts alone do not establish equivalence.

use std::{hint::black_box, path::Path, time::Instant};

use idenso::tensor::{SymbolicTensor, contract::ContractSettings};
use spenso::network::{parsing::AtomStructureExt, tags::SPENSO_TAG};
use symbolica::atom::{Atom, AtomCore, AtomView};

#[derive(Clone, Copy, Debug)]
enum Method {
    Clone,
    RepeatedIndices,
    NotationToDots,
    MetricContract,
    Contract,
    Epsilon,
    Chainify,
    TypedGammaChains,
    TypedGamma,
}
const METHODS: [Method; 9] = [
    Method::Clone,
    Method::RepeatedIndices,
    Method::NotationToDots,
    Method::MetricContract,
    Method::Contract,
    Method::Epsilon,
    Method::Chainify,
    Method::TypedGammaChains,
    Method::TypedGamma,
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
        NotationToDots => SymbolicTensor::infer(view.to_owned())
            .unwrap()
            .to_dots()
            .unwrap()
            .into_expression(),
        MetricContract => SymbolicTensor::infer(view.to_owned())
            .unwrap()
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        Contract => SymbolicTensor::infer(view.to_owned())
            .unwrap()
            .contract(idenso::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        Epsilon => SymbolicTensor::infer(view.to_owned())
            .unwrap()
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                epsilon: true,
                ..Default::default()
            })
            .unwrap()
            .contract(ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        Chainify => idenso::reference_cases::chain_notation_phase(view),
        TypedGammaChains => idenso::tensor::SymbolicTensor::infer((view).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                gamma: Some(idenso::dirac::GammaSimplifySettings {
                    output: idenso::dirac::GammaOutput::Chains,
                    ..Default::default()
                }),
                epsilon: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        TypedGamma => idenso::tensor::SymbolicTensor::infer((view).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&idenso::tensor::AlgebraSettings {
                gamma: Some(idenso::dirac::GammaSimplifySettings::default()),
                epsilon: true,
                ..Default::default()
            })
            .unwrap()
            .contract(idenso::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
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
    if arguments
        .get(1)
        .is_some_and(|mode| mode == "--reduce-input")
    {
        reduce_input(&arguments[2..]);
        return;
    }
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
    SPENSO_TAG.tensor_symbol("spenso::generic");
    // Typed controls use atomic dimensions and one coherent open interface.
    // Their scalar spectator remains factored throughout every stage.
    let controls = [
        (
            "generic_bispinor_factored",
            "(x+y)^8*generic((x+y)^8,metadata,bis(Ns,a),bis(Ns,b))",
        ),
        (
            "generic_other_representation_factored",
            "(x+y)^8*generic((x+y)^8,metadata,cof(Nc,a),dind(cof(Nc,b)))",
        ),
        (
            "generic_other_representation_reversed",
            "(x+y)^8*generic((x+y)^8,metadata,dind(cof(Nc,a)),cof(Nc,b))",
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
            if matches!(method, Method::TypedGamma) && !has_gamma {
                continue;
            }
            let output = run(method, input.as_view()).atom();
            if matches!(method, Method::TypedGamma) && variant == "gamma_terminal" {
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
        #[cfg(feature = "reference-cases")]
        let (_, diagnostic) = idenso::reference_cases::timing::measure(
            &format!("phase/{case}/{variant}/{method:?}"),
            || drop(run(*method, input.as_view())),
        );
        #[cfg(not(feature = "reference-cases"))]
        let diagnostic = "null";
        println!(
            "{{\"case\":{case:?},\"input_variant\":{variant:?},\"method\":\"{method:?}\",\"nodes\":{nodes},\"terms\":{},\"output_terms\":{output_terms},\"changed\":{changed},\"repetitions\":{},\"samples_ns\":{:?},\"median_ns\":{median_ns},\"diagnostic\":{diagnostic}}}",
            input.nterms(),
            repetitions[index],
            timings[index]
        );
    }
}

/// Diagnose an exact production input without distributing its numerator.
/// `--reduce-input INPUT OUTPUT_DIRECTORY [gamma|algebra|algebra-dots] [collected|nested]`
/// records admission and each exact continuation separately. `algebra` matches
/// the public gamma/colour defaults with full contraction; `algebra-dots` also
/// forms canonical dots. `gamma` retains the trace-only diagnostic.
fn reduce_input(arguments: &[String]) {
    use idenso::tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus};

    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
    let output = Path::new(&arguments[1]);
    std::fs::create_dir_all(output).unwrap();
    let atom = if Path::new(&arguments[0])
        .extension()
        .is_some_and(|extension| extension == "raw")
    {
        Atom::import(&mut std::fs::File::open(&arguments[0]).unwrap(), None).unwrap()
    } else {
        for head in ["gammalooprs::K", "gammalooprs::Q", "gammalooprs::P"] {
            SPENSO_TAG.rank_one_tensor_symbol(head);
        }
        let _ = spenso::index_symbol!("gammalooprs::hedge");
        let _ = spenso::index_symbol!("gammalooprs::edge");
        let source = std::fs::read_to_string(&arguments[0]).unwrap();
        Atom::parse(
            &source,
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let (input, clocks) = idenso::reference_cases::timing::measure("admission", || {
        SymbolicTensor::infer(atom).unwrap()
    });
    std::fs::write(output.join("admission.json"), clocks.to_string()).unwrap();
    let mode = arguments.get(2).map(String::as_str);
    assert!(matches!(
        mode,
        None | Some("gamma" | "algebra" | "algebra-dots")
    ));
    let coefficients = arguments.get(3).map(String::as_str);
    assert!(matches!(coefficients, None | Some("collected" | "nested")));
    let settings = AlgebraSettings {
        gamma: mode.map(|_| Default::default()),
        color: matches!(mode, Some("algebra" | "algebra-dots")).then(Default::default),
        collect_coefficients: coefficients != Some("nested"),
        contract: if mode == Some("algebra") {
            AlgebraContraction::Fully
        } else {
            AlgebraContraction::Dots
        },
        ..Default::default()
    };
    let mut current = input;
    for step in 1..=8 {
        let start = Instant::now();
        let reduced = current.simplify_algebra(&settings).unwrap();
        let elapsed = start.elapsed().as_secs_f64();
        println!(
            "step={step} seconds={elapsed:.9} status={:?} atom_bytes={}",
            reduced.reduction_status(),
            reduced.expression().as_view().get_byte_size()
        );
        let (diagnostic, clocks) =
            idenso::reference_cases::timing::measure("production reduction", || {
                current.simplify_algebra(&settings).unwrap()
            });
        assert_eq!(diagnostic.expression(), reduced.expression());
        assert_eq!(diagnostic.reduction_status(), reduced.reduction_status());
        std::fs::write(output.join(format!("step-{step}.json")), clocks.to_string()).unwrap();
        std::fs::write(
            output.join(format!("step-{step}.expr")),
            reduced.expression().to_plain_string(),
        )
        .unwrap();
        if reduced.reduction_status() == ReductionStatus::Complete {
            let start = Instant::now();
            let rerun = reduced.simplify_algebra(&settings).unwrap();
            assert_eq!(rerun.reduction_status(), ReductionStatus::Complete);
            println!(
                "rerun_seconds={:.9} stable={}",
                start.elapsed().as_secs_f64(),
                rerun.expression() == reduced.expression()
            );
            break;
        }
        assert!(
            reduced.expression() != current.expression(),
            "unfinished result made no progress"
        );
        current = reduced;
    }
}
