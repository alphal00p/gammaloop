//! Shared metric/vector contraction and complete pipeline timings on several shapes.
//!
//! `cargo run -p idenso --profile dev-optim --example metric_contraction_benchmark
//! -- OUTPUT_DIRECTORY [samples=5] [target_batch_ms=15]`
//!
//! The isolated metric-only contractor is compiled from its production source so this example
//! can time the private operation without a new public API. The two public
//! Schoonschip methods and gamma simplification call the linked production library.
//! The repeated-index check is timed separately; it never guards the contractor.
//! Input creation, validation and output snapshots are outside timing. Timed calls
//! include output destruction. Settings with chain-like functions enabled also
//! cover contractions into ordered chains and cyclic/symmetric traces.
//!
//! JSONL records whether the operation changed its input. Exact input/output
//! snapshots permit before/after comparisons; different symbolic forms still need
//! an independent tensor-network equality check, not a term-count comparison.

use std::{hint::black_box, path::Path, time::Instant};

pub use idenso::W_;
#[cfg(test)]
use idenso::representations;
use idenso::{
    dirac::GammaSimplifier,
    epsilon::EpsilonSimplifier,
    gamma, gamma5,
    representations::Bispinor,
    shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
};
use spenso::{
    g,
    network::parsing::AtomStructureExt,
    network::tags::SPENSO_TAG,
    p, q,
    structure::representation::{Minkowski, RepName},
    trace,
};
use symbolica::{
    atom::{Atom, AtomView},
    function, symbol,
};

#[path = "../src/shorthands/schoonschip/slot_contraction.rs"]
mod slot_contraction;

struct Case {
    name: String,
    expression: Atom,
    gamma_input: Option<Atom>,
}

impl Case {
    fn parse(name: &str, source: &str) -> Self {
        Self {
            name: name.into(),
            expression: Atom::parse(
                source,
                "spenso",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap(),
            gamma_input: None,
        }
    }

    fn cases() -> Vec<Self> {
        // Register canonical rank-one heads before parsing the vector cases.
        SPENSO_TAG.rank_one_tensor_symbol("cleanup_benchmark::V");
        SPENSO_TAG.rank_one_tensor_symbol("cleanup_benchmark::W");
        let mut cases = vec![
            Self::parse("metric_hit", "g(mink(4,a),mink(4,b))*T(mink(4,b))"),
            Self::parse("metric_miss", "g(mink(4,a),mink(4,b))*T(mink(4,c))"),
            Self::parse("dimension_miss", "g(mink(4,a),mink(4,b))*T(mink(5,b))"),
            Self::parse("representation_miss", "g(mink(4,a),mink(4,b))*T(euc(4,b))"),
            Self::parse("dual_hit", "g(lor(4,a),dind(lor(4,b)))*T(lor(4,b))"),
            Self::parse(
                "dual_variance_miss",
                "g(lor(4,a),dind(lor(4,b)))*T(dind(lor(4,b)))",
            ),
            Self::parse(
                "two_metrics",
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*T(mink(4,c))",
            ),
            Self::parse(
                "closed_metric_loop",
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*g(mink(4,c),mink(4,a))",
            ),
            Self::parse(
                "epsilon_partner",
                "g(mink(4,a),mink(4,b))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
            ),
            Self::parse(
                "chain_partner",
                "g(mink(4,a),mink(4,b))*chain(bis(4,s),bis(4,t),T(mink(4,b),in,out),U(in,out))",
            ),
            Self::parse(
                "trace_partner",
                "g(mink(4,a),mink(4,b))*trace(bis(4),cyclic(T(mink(4,b),in,out),U(in,out)))",
            ),
            Self::parse(
                "sum_of_products",
                "x*g(mink(4,a),mink(4,b))*T(mink(4,b))+y*g(mink(4,a),mink(4,c))*U(mink(4,c))",
            ),
            Self::parse(
                "factored_sum_partner",
                "g(mink(4,a),mink(4,b))*(T(mink(4,b))+U(mink(4,b)))",
            ),
            Self::parse("compact_metric", "g(mink(4,a),P(mink(4)))*T(mink(4,a))"),
            Self::parse("no_metric", "T(mink(4,a))*U(mink(4,a))"),
            Self::parse(
                "vector_tensor",
                "cleanup_benchmark::V(mink(4,a))*T(mink(4,a))",
            ),
            Self::parse(
                "vector_miss",
                "cleanup_benchmark::V(mink(4,a))*T(mink(4,b))",
            ),
            Self::parse(
                "vector_dimension_miss",
                "cleanup_benchmark::V(mink(4,a))*T(mink(5,a))",
            ),
            Self::parse(
                "vector_dual",
                "cleanup_benchmark::V(lor(4,a))*T(dind(lor(4,a)))",
            ),
            Self::parse(
                "vector_variance_miss",
                "cleanup_benchmark::V(lor(4,a))*T(lor(4,a))",
            ),
            Self::parse(
                "vector_parameter",
                "cleanup_benchmark::V(a,mink(4,a))*T(mink(4,a))",
            ),
            Self::parse(
                "vector_dot",
                "cleanup_benchmark::V(mink(4,a))*cleanup_benchmark::W(mink(4,a))",
            ),
            Self::parse(
                "vector_chain",
                "cleanup_benchmark::V(mink(4,a))*chain(bis(4,s),bis(4,t),T(mink(4,a),in,out),U(in,out))",
            ),
            Self::parse(
                "vector_trace",
                "cleanup_benchmark::V(mink(4,a))*trace(bis(4),cyclic(T(mink(4,a),in,out),U(in,out)))",
            ),
            Self::parse(
                "vector_symmetric_trace",
                "cleanup_benchmark::V(mink(4,a))*trace(bis(4),sym(T(mink(4,a),in,out),U(in,out)))",
            ),
            Self::parse(
                "vector_epsilon",
                "cleanup_benchmark::V(mink(4,a))*epsilon(mink(4,a),mink(4,b),mink(4,c),mink(4,d))",
            ),
            Self::parse(
                "tagged_compact_metric",
                "g(mink(4,a),cleanup_benchmark::V(mink(4)))*T(mink(4,a))",
            ),
            Self::parse(
                "epsilon_square",
                "epsilon(mink(4,a),mink(4,b),mink(4,c),mink(4,d))^2",
            ),
            Self::parse(
                "epsilon_distinct_pair",
                "epsilon(mink(4,a),mink(4,b),mink(4,c),mink(4,d))*epsilon(mink(4,e),mink(4,f),mink(4,h),mink(4,j))",
            ),
            Self::parse(
                "epsilon_compact_spectator",
                "epsilon(mink(4,a),mink(4,b),mink(4,c),mink(4,d))*g(mink(4,e),cleanup_benchmark::V(mink(4)))",
            ),
        ];
        cases.push(Self::parse(
            "symmetric_trace_partner",
            "g(mink(4,a),mink(4,b))*trace(bis(4),sym(T(mink(4,b),in,out),U(in,out)))",
        ));
        let spectators = (0..48)
            .map(|i| format!("S{i}(x{i})"))
            .collect::<Vec<_>>()
            .join("*");
        cases.push(Self::parse(
            "many_spectators_hit",
            &format!("{spectators}*g(mink(4,a),mink(4,b))*T(mink(4,b))"),
        ));
        cases.push(Self::parse(
            "many_spectators_miss",
            &format!("{spectators}*g(mink(4,a),mink(4,b))*T(mink(4,c))"),
        ));
        let metrics = (0..8)
            .map(|i| format!("g(mink(4,i{i}),mink(4,i{}))", i + 1))
            .collect::<Vec<_>>()
            .join("*");
        cases.push(Self::parse(
            "eight_metrics",
            &format!("{metrics}*T(mink(4,i8))"),
        ));
        let mink = Minkowski {}.new_rep(4);
        let spin = Bispinor {}.new_rep(4).to_symbolic([]);
        for length in [6, 8, 10, 12] {
            let mut factors = vec![gamma5!()];
            factors.extend(
                (0..length)
                    .map(|i| gamma!(mink.pattern(symbol!(&format!("metric_benchmark::mu{i}"))))),
            );
            let input = trace!(&spin; factors);
            cases.push(Self {
                name: format!("axial_{length}"),
                expression: input.simplify_gamma(),
                gamma_input: Some(input),
            });
        }
        let indices: Vec<_> = (0..8)
            .map(|i| mink.pattern(symbol!(&format!("metric_benchmark::nu{i}"))))
            .collect();
        let free = trace!(&spin;indices.iter().map(|index|gamma!(index)));
        cases.push(Self {
            name: "free_trace_8".into(),
            expression: free.simplify_gamma(),
            gamma_input: Some(free),
        });
        let p = p!(mink.to_symbolic([]));
        let q = q!(mink.to_symbolic([]));
        let compact = trace!(&spin;(0..8).map(|i|gamma!(if i<2{&p}else{&q})));
        cases.push(Self {
            name: "compact_paired_8".into(),
            expression: compact.simplify_gamma(),
            gamma_input: Some(compact),
        });
        let a = mink.pattern(symbol!("metric_benchmark::a"));
        let b = mink.pattern(symbol!("metric_benchmark::b"));
        let c = mink.pattern(symbol!("metric_benchmark::c"));
        let d = mink.pattern(symbol!("metric_benchmark::d"));
        let momenta: Vec<_> = (0..10)
            .map(|i| {
                function!(
                    SPENSO_TAG.rank_one_tensor_symbol(&format!("metric_benchmark::p{i}")),
                    mink.to_symbolic([])
                )
            })
            .collect();
        let p = &momenta;
        let odd = [&a, &p[0], &a, &p[1], &p[2], &p[3], &p[4], &p[5]];
        let even = [
            &a, &p[0], &p[1], &p[2], &p[3], &a, &p[4], &p[5], &p[6], &p[7],
        ];
        let crossing = [
            &a, &p[0], &b, &p[1], &p[2], &a, &p[3], &p[4], &b, &p[5], &p[6], &p[7],
        ];
        for (name, factors, axial) in [
            ("order_sensitive_12", crossing.as_slice(), false),
            ("ordinary_odd_8", odd.as_slice(), false),
            ("ordinary_even4_10", even.as_slice(), false),
            ("axial_odd_8", odd.as_slice(), true),
            ("axial_even4_10", even.as_slice(), true),
            ("axial_crossing_12", crossing.as_slice(), true),
            (
                "ordinary_even2_8",
                [&a, &p[0], &p[1], &a, &p[2], &p[3], &p[4], &p[5]].as_slice(),
                false,
            ),
            (
                "ordinary_cyclic_12",
                [
                    &a, &p[0], &p[1], &p[2], &p[3], &p[4], &p[5], &p[6], &p[7], &p[8], &p[9], &a,
                ]
                .as_slice(),
                false,
            ),
        ] {
            let factors = axial
                .then(|| gamma5!())
                .into_iter()
                .chain(factors.iter().map(|p| gamma!(*p)));
            let input = trace!(&spin; factors);
            cases.push(Self {
                name: name.into(),
                expression: input.simplify_gamma(),
                gamma_input: Some(input),
            });
        }
        let external = g!(&a, &c)
            * g!(&b, &d)
            * trace!(&spin; [&a, &p[0], &b, &p[1], &p[2], &c, &p[3], &p[4], &d, &p[5], &p[6], &p[7]].map(|p| gamma!(p)));
        let free = trace!(&spin; (0..12).map(|i| gamma!(mink.pattern(symbol!(&format!("metric_benchmark::mu{i}"))))));
        for (name, input) in [("external_metrics_12", external), ("free_12", free)] {
            cases.push(Self {
                name: name.into(),
                expression: input.simplify_gamma(),
                gamma_input: Some(input),
            });
        }
        cases
    }
}

#[derive(Clone, Copy)]
enum Method {
    Repeated,
    MetricCore,
    GuardedMetricCore,
    MetricSettings,
    FullSchoonschip,
    EpsilonCleanup,
    FullGamma,
}
impl Method {
    const ALL: [Self; 7] = [
        Self::Repeated,
        Self::MetricCore,
        Self::GuardedMetricCore,
        Self::MetricSettings,
        Self::FullSchoonschip,
        Self::EpsilonCleanup,
        Self::FullGamma,
    ];
    fn name(self) -> &'static str {
        match self {
            Self::Repeated => "repeated_indices_only",
            Self::MetricCore => "metric_contractor_only",
            Self::GuardedMetricCore => "guarded_metric_contractor",
            Self::MetricSettings => "metric_schoonschip",
            Self::FullSchoonschip => "full_schoonschip",
            Self::EpsilonCleanup => "epsilon_cleanup",
            Self::FullGamma => "full_gamma",
        }
    }
    fn input(self, case: &Case) -> Option<AtomView<'_>> {
        match self {
            Self::FullGamma => case.gamma_input.as_ref().map(Atom::as_view),
            _ => Some(case.expression.as_view()),
        }
    }
    fn apply(self, expression: AtomView<'_>) -> Atom {
        match self {
            Self::Repeated => Atom::num(i64::from(expression.has_repeated_explicit_indices())),
            Self::MetricCore => slot_contraction::SlotContraction::run(expression, true, false),
            Self::GuardedMetricCore => {
                if expression.has_repeated_explicit_indices() {
                    slot_contraction::SlotContraction::run(expression, true, false)
                } else {
                    expression.to_owned()
                }
            }
            Self::MetricSettings => expression.schoonschip_with_settings(
                &SchoonschipSettings::default()
                    .without_rank1_tensors()
                    .with_chain_like_functions(),
            ),
            Self::FullSchoonschip => expression.schoonschip_with_settings(
                &SchoonschipSettings::default().with_chain_like_functions(),
            ),
            Self::EpsilonCleanup => expression.simplify_epsilon(),
            Self::FullGamma => expression.simplify_gamma(),
        }
    }
    fn run(self, expression: AtomView<'_>) {
        if matches!(self, Self::Repeated) {
            black_box(black_box(expression).has_repeated_explicit_indices());
        } else {
            let _ = black_box(self.apply(black_box(expression)));
        }
    }
    fn batch(self, expression: AtomView<'_>, repetitions: usize) -> f64 {
        let started = Instant::now();
        for _ in 0..repetitions {
            self.run(expression);
        }
        started.elapsed().as_secs_f64() * 1e9 / repetitions as f64
    }
}

fn main() {
    let arguments: Vec<_> = std::env::args().collect();
    let directory = Path::new(arguments.get(1).expect(
        "usage: metric_contraction_benchmark OUTPUT_DIRECTORY [samples=5] [target_batch_ms=15]",
    ));
    let samples: usize = arguments.get(2).map_or(5, |s| s.parse().unwrap());
    let target_ms: f64 = arguments.get(3).map_or(15., |s| s.parse().unwrap());
    assert!(samples > 0 && target_ms.is_finite() && target_ms > 0.);
    std::fs::create_dir_all(directory).unwrap();
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
    for case in Case::cases() {
        std::fs::write(
            directory.join(format!("{}.input", case.name)),
            case.expression.to_plain_string(),
        )
        .unwrap();
        if let Some(input) = &case.gamma_input {
            std::fs::write(
                directory.join(format!("{}.gamma_input", case.name)),
                input.to_plain_string(),
            )
            .unwrap();
        }
        let methods: Vec<_> = Method::ALL
            .into_iter()
            .filter(|m| m.input(&case).is_some())
            .collect();
        let mut changed = Vec::new();
        let mut repetitions = Vec::new();
        let mut output_terms = Vec::new();
        for &method in &methods {
            let input = method.input(&case).unwrap();
            let output = method.apply(input);
            changed.push(output.as_view() != input);
            output_terms.push(output.nterms());
            std::fs::write(
                directory.join(format!("{}.{}.output", case.name, method.name())),
                output.to_plain_string(),
            )
            .unwrap();
            if matches!(
                method,
                Method::MetricCore
                    | Method::GuardedMetricCore
                    | Method::MetricSettings
                    | Method::FullSchoonschip
                    | Method::EpsilonCleanup
            ) {
                assert_eq!(
                    method.apply(output.as_view()),
                    output,
                    "{} {} is not idempotent",
                    case.name,
                    method.name()
                );
            }
            let ns = method.batch(input, 1);
            repetitions.push((target_ms * 1e6 / ns).ceil().clamp(1., 500_000.) as usize);
        }
        let mut timings: Vec<Vec<f64>> = vec![Vec::new(); methods.len()];
        for sample in 0..samples {
            for offset in 0..methods.len() {
                let m = (sample + offset) % methods.len();
                timings[m].push(methods[m].batch(methods[m].input(&case).unwrap(), repetitions[m]));
            }
        }
        for (m, method) in methods.iter().enumerate() {
            let mut sorted = timings[m].clone();
            sorted.sort_by(f64::total_cmp);
            let median = if samples.is_multiple_of(2) {
                (sorted[samples / 2 - 1] + sorted[samples / 2]) / 2.
            } else {
                sorted[samples / 2]
            };
            println!(
                "{{\"case\":{:?},\"method\":{:?},\"terms\":{},\"output_terms\":{},\"changed\":{},\"repeated_indices\":{},\"repetitions\":{},\"samples_ns\":{:?},\"median_ns\":{median}}}",
                case.name,
                method.name(),
                method.input(&case).unwrap().nterms(),
                output_terms[m],
                changed[m],
                case.expression.has_repeated_explicit_indices(),
                repetitions[m],
                timings[m]
            );
        }
    }
}
