//! Compare explicit-index discovery, syntactic inference, and partial parsing.
//!
//! Run on one CPU, for example:
//! `taskset -c 2 cargo run -p idenso --profile dev-optim --example
//! structure_matching_benchmark -- /tmp/structure-matching 7 20`.
//!
//! Arguments are an output directory, sample count, and target batch milliseconds.
//! JSONL on stdout includes every sample and the median. Exact inputs, inferred
//! structures (including canonicalization layouts), and every top-level term's
//! structure are saved separately for comparisons between builds. Input creation,
//! gamma simplification, validation, and snapshot formatting are outside timing.
//! Method order rotates between samples; result destruction is included.
//!
//! `infer_fast` uses only the first summand for structure discovery, but validates
//! chain/trace nesting throughout the expression.
//! `infer_all_terms` invokes Fast inference separately on every top-level term;
//! nested sums retain Fast inference's first-summand semantics. Partial parsing
//! preserves all summands and keeps shorthand opaque, without executing a network.

use std::{hint::black_box, path::Path, time::Instant};

use idenso::{
    dirac::GammaSimplifier, gamma, gamma5, representations::Bispinor, tensor::SymbolicNetParse,
};
use spenso::{
    network::parsing::{AtomStructureExt, ParseSettings, ShorthandParsing, StructureInferenceMode},
    structure::{
        OrderedStructure, TensorStructure,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomView},
    symbol,
};

struct Case {
    name: String,
    expression: Atom,
    repeated: bool,
}

impl Case {
    fn parse(name: &str, source: &str, repeated: bool) -> Self {
        Self {
            name: name.into(),
            expression: Atom::parse(
                source,
                "spenso",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap(),
            repeated,
        }
    }

    fn cases() -> Vec<Self> {
        let mut cases = vec![
            Self::parse("open_metric", "g(mink(4,a),mink(4,b))", false),
            Self::parse(
                "contracted_metrics",
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))",
                true,
            ),
            Self::parse(
                "dual_slots",
                "T(cof(3,a),dind(cof(3,b)))*U(cof(3,b),dind(cof(3,c)))",
                true,
            ),
            Self::parse("independent_sum", "T(mink(4,a))+U(mink(4,a))", false),
            Self::parse(
                "factored_open_sum",
                "(x*T(mink(4,a))+y*U(mink(4,a)))*V(mink(4,b))",
                false,
            ),
            Self::parse(
                "factored_contracted_sum",
                "(x*T(mink(4,a))+y*U(mink(4,a)))*V(mink(4,a),mink(4,b))",
                true,
            ),
            Self::parse("compact_inner_product", "g(P(mink(4)),Q(mink(4)))", false),
            Self::parse("compact_open_vector", "g(P(mink(4)),mink(4,a))", false),
            Self::parse("tensor_power", "T(mink(4,a),mink(4,b))^2", true),
        ];
        for heads in [64, 256, 1024] {
            let scalars = (0..heads)
                .map(|i| format!("f{i}(x{i},y{i})"))
                .collect::<Vec<_>>()
                .join("*");
            let name = if heads == 64 {
                "scalar_heavy".into()
            } else {
                format!("scalar_heavy_{heads}")
            };
            cases.push(Self::parse(
                &name,
                &format!("{scalars}*T(mink(4,a),mink(4,b))"),
                false,
            ));
        }
        let product = (0..64)
            .map(|i| format!("T{i}(mink(4,i{i}))"))
            .collect::<Vec<_>>()
            .join("*");
        cases.push(Self::parse("long_open_product", &product, false));
        cases.push(Self::parse(
            "long_contracted_product",
            &format!("{product}*U(mink(4,i63))"),
            true,
        ));
        let mink = Minkowski {}.new_rep(4);
        let spin = Bispinor {}.new_rep(4).to_symbolic([]);
        for length in [6, 8, 10, 12] {
            let mut factors = vec![gamma5!()];
            factors.extend((0..length).map(|i| {
                let index = mink.pattern(symbol!(&format!("structure_benchmark::mu{i}")));
                gamma!(index)
            }));
            cases.push(Self {
                name: format!("axial_{length}"),
                expression: trace!(&spin; factors).simplify_gamma(),
                repeated: false,
            });
        }
        cases
    }

    fn terms(&self) -> Vec<AtomView<'_>> {
        match self.expression.as_view() {
            AtomView::Add(sum) => sum.iter().collect(),
            expression => vec![expression],
        }
    }

    fn validate_and_snapshot(&self, directory: &Path, settings: &ParseSettings) {
        assert_eq!(
            self.expression.has_repeated_explicit_indices(),
            self.repeated,
            "{}",
            self.name
        );
        let infer = |expression: AtomView<'_>| {
            let structure = match expression
                .infer_structure::<OrderedStructure>(StructureInferenceMode::Fast)
            {
                Ok(structure) => structure,
                Err(error) => return format!("error: {error}\n"),
            };
            let slots = structure
                .canonical()
                .external_structure()
                .iter()
                .map(|slot| slot.to_atom().to_plain_string())
                .collect::<Vec<_>>();
            format!("slots: {slots:?}\nlayout: {:?}\n", structure.layout())
        };
        let mut snapshot = format!(
            "repeated: {}\nroot:\n{}",
            self.repeated,
            infer(self.expression.as_view())
        );
        for (number, term) in self.terms().into_iter().enumerate() {
            snapshot.push_str(&format!("term {number}:\n{}", infer(term)));
        }
        let network = self
            .expression
            .parse_to_symbolic_net::<AbstractIndex>(settings);
        snapshot.push_str(&format!(
            "partial_parse_error: {:?}\n",
            network.err().map(|error| error.to_string())
        ));
        std::fs::write(
            directory.join(format!("{}.input", self.name)),
            self.expression.to_plain_string(),
        )
        .unwrap();
        std::fs::write(directory.join(format!("{}.structure", self.name)), snapshot).unwrap();
    }
}

#[derive(Clone, Copy)]
enum Method {
    Repeated,
    InferFast,
    InferAllTerms,
    PartialParse,
}

impl Method {
    const ALL: [Self; 4] = [
        Self::Repeated,
        Self::InferFast,
        Self::InferAllTerms,
        Self::PartialParse,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::Repeated => "repeated_explicit_indices",
            Self::InferFast => "infer_fast",
            Self::InferAllTerms => "infer_all_terms",
            Self::PartialParse => "partial_parse_all_terms",
        }
    }

    fn run(self, case: &Case, settings: &ParseSettings) {
        let expression = black_box(case.expression.as_view());
        match self {
            Self::Repeated => {
                black_box(expression.has_repeated_explicit_indices());
            }
            Self::InferFast => {
                let _ = black_box(
                    expression.infer_structure::<OrderedStructure>(StructureInferenceMode::Fast),
                );
            }
            Self::InferAllTerms => match expression {
                AtomView::Add(sum) => {
                    for term in sum.iter() {
                        let _ = black_box(
                            term.infer_structure::<OrderedStructure>(StructureInferenceMode::Fast),
                        );
                    }
                }
                expression => {
                    let _ = black_box(
                        expression
                            .infer_structure::<OrderedStructure>(StructureInferenceMode::Fast),
                    );
                }
            },
            Self::PartialParse => {
                let _ = black_box(expression.parse_to_symbolic_net::<AbstractIndex>(settings));
            }
        }
    }

    fn batch(self, case: &Case, settings: &ParseSettings, repetitions: usize) -> f64 {
        let start = Instant::now();
        for _ in 0..repetitions {
            self.run(case, settings);
        }
        start.elapsed().as_secs_f64() * 1e9 / repetitions as f64
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let directory = Path::new(args.get(1).expect(
        "usage: structure_matching_benchmark OUTPUT_DIRECTORY [samples=7] [target_batch_ms=20]",
    ));
    let samples: usize = args.get(2).map_or(7, |s| s.parse().unwrap());
    let target_ms: f64 = args.get(3).map_or(20., |s| s.parse().unwrap());
    assert!(samples > 0 && target_ms.is_finite() && target_ms > 0.);
    std::fs::create_dir_all(directory).unwrap();
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
    let settings = ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..ParseSettings::default()
    };
    for case in Case::cases() {
        case.validate_and_snapshot(directory, &settings);
        let inference_ok = case
            .expression
            .infer_structure::<OrderedStructure>(StructureInferenceMode::Fast)
            .is_ok();
        let partial_parse_ok = case
            .expression
            .parse_to_symbolic_net::<AbstractIndex>(&settings)
            .is_ok();
        let repetitions = Method::ALL.map(|method| {
            method.run(&case, &settings);
            let ns = method.batch(&case, &settings, 3);
            (target_ms * 1e6 / ns).ceil().clamp(1., 1_000_000.) as usize
        });
        let mut timings: [Vec<f64>; 4] = std::array::from_fn(|_| Vec::new());
        for sample in 0..samples {
            for offset in 0..Method::ALL.len() {
                let index = (sample + offset) % Method::ALL.len();
                timings[index].push(Method::ALL[index].batch(&case, &settings, repetitions[index]));
            }
        }
        for (index, method) in Method::ALL.iter().enumerate() {
            let mut sorted = timings[index].clone();
            sorted.sort_by(f64::total_cmp);
            let median = if samples.is_multiple_of(2) {
                (sorted[samples / 2 - 1] + sorted[samples / 2]) / 2.
            } else {
                sorted[samples / 2]
            };
            println!(
                "{{\"case\":{:?},\"terms\":{},\"repeated\":{},\"inference_ok\":{inference_ok},\"partial_parse_ok\":{partial_parse_ok},\"method\":{:?},\"repetitions\":{},\"samples_ns\":{:?},\"median_ns\":{}}}",
                case.name,
                case.expression.nterms(),
                case.repeated,
                method.name(),
                repetitions[index],
                timings[index],
                median
            );
        }
    }
}
