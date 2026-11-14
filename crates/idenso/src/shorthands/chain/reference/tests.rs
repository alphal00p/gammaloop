//! The pre-checked chainify must return exactly what the matcher-only
//! reference returns. Every chainify a test reaches is compared already;
//! these tests add a broad generated corpus, the stored workflow fixtures
//! and the simplification pipelines that call chainify.

use std::time::{Duration, Instant};

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    structure::representation::{LibraryRep, Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    symbol,
};

use super::{MatcherCounts, Mode, cross_checks, record, take_counts, with_mode};
use crate::{
    color::{CS, ColorSimplifySettings},
    cook::{CookMode, CookSettings},
    dirac::{AGS, GammaSimplifySettings},
    representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet},
    shorthands::chain::Chain,
    tensor::{AlgebraSettings, ContractSettings, SymbolicTensor},
    test_support::test_initialize,
};

/// Register every head with its attributes before an input names it.
fn initialize() {
    test_initialize();
    let _ = (CS.t, CS.f, CS.d, AGS.gamma, AGS.gamma0, ETS.metric, T.chain);
    for vector in ["p", "q", "Q", "K", "P"] {
        T.rank_one_tensor_symbol(&format!("prefilter::{vector}"));
    }
    let _ = (
        symbol!("prefilter::sym"; Symmetric),
        symbol!("prefilter::antisym"; Antisymmetric),
        symbol!("prefilter::cyclic"; Cyclesymmetric),
        symbol!("prefilter::den"; Scalar),
        symbol!("prefilter::scalar"; Scalar),
        // The stored workflows label ports by gammaloop's half-edge functions.
        symbol!("prefilter_gl::hedge", tags = [T.index.clone()]),
        symbol!("prefilter_gl::edge", tags = [T.index.clone()]),
    );
}

fn parse(source: &str) -> Atom {
    Atom::parse(
        source,
        "spenso",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap_or_else(|error| panic!("{error}: {source}"))
}

/// Momenta of the stored workflow inputs become rank-one vectors, as the
/// gammaloop momenta are in production. Their other symbols move to this
/// test's namespaces, so that no other test sees them declared differently.
fn parse_workflow(source: &str) -> Atom {
    parse(
        &source
            .replace("gammalooprs::Q(", "prefilter::Q(")
            .replace("gammalooprs::K(", "prefilter::K(")
            .replace("gammalooprs::P(", "prefilter::P(")
            .replace("gammalooprs::", "prefilter_gl::")
            .replace("subclass_preservation_benchmark::", "prefilter_gl::")
            .replace("UFO::", "prefilter_ufo::"),
    )
}

fn representations() -> Vec<LibraryRep> {
    vec![
        Bispinor {}.into(),
        Minkowski {}.into(),
        ColorFundamental {}.into(),
        LibraryRep::from(ColorFundamental {}).dual(),
        ColorAdjoint {}.into(),
        ColorSextet {}.into(),
        LibraryRep::from(ColorSextet {}).dual(),
    ]
}

/// Factors sharing index labels, so that products connect into partial,
/// open and closed chains of every representation.
const FACTORS: &[&str] = &[
    // Dirac matrices, slashed momenta, existing words and traces.
    "gamma(bis(4,i1),bis(4,i2),mink(4,mu1))",
    "gamma(bis(4,i2),bis(4,i3),mink(4,mu2))",
    "gamma(bis(4,i3),bis(4,i1),mink(4,mu1))",
    "gamma(bis(4,i2),bis(4,i1),prefilter::p(mink(4)))",
    "gamma(bis(4,i4),bis(4,i2),prefilter::Q(1,mink(4)))",
    "gamma(bis(D,j1),bis(D,j2),mink(D,nu1))",
    "gamma(bis(D,j2),bis(D,j1),prefilter::Q(2,mink(D)))",
    "gamma0(bis(4,i3),bis(4,i4))",
    "gamma5(bis(4,i3),bis(4,i4))",
    "projm(bis(4,i4),bis(4,i1))",
    "g(bis(4,i1),bis(4,i4))",
    "chain(bis(4,i1),bis(4,i2),gamma(in,out,mink(4,mu1)),gamma(in,out,mink(4,mu2)))",
    "chain(bis(4,i2),bis(4,i3),gamma(in,out,prefilter::p(mink(4))))",
    "chain(bis(4,i3),bis(4,i3),gamma(in,out,mink(4,mu3)))",
    "chain(bis(4,i4),bis(4,i1))",
    "trace(bis(4),gamma(in,out,mink(4,mu1)),gamma(in,out,mink(4,mu2)))",
    "prefilter::u(1,bis(4,i1))*prefilter::ubar(2,bis(4,i3))",
    // Lorentz metrics, vectors and inner products.
    "g(mink(4,mu1),mink(4,mu2))",
    "g(mink(D,nu1),mink(D,nu2))",
    "g(prefilter::p(mink(4)),prefilter::q(mink(4)))",
    "dot(prefilter::p(mink(4)),prefilter::q(mink(4)))",
    "dot(prefilter::Q(1,mink(D)),prefilter::Q(2,mink(D)))",
    "dot(prefilter::Q(1,mink(D)),prefilter::Q(1,mink(D)))",
    "prefilter::p(mink(4,mu1))",
    "prefilter::Q(2,mink(D,nu2))",
    "g(gamma(bis(4,i1),bis(4,i2),mink(4)),prefilter::p(mink(4)))",
    "dot(prefilter::q(mink(4)),gamma(bis(4,i2),bis(4,i3),mink(4)))",
    "dot(chain(bis(4,i3),bis(4,i4),gamma(in,out,mink(4))),prefilter::p(mink(4)))",
    "g(prefilter::p(mink(4)),prefilter::q(mink(4))*gamma(bis(4,i1),bis(4,i1),mink(4,mu3)))",
    "prefilter::den(dot(prefilter::Q(1,mink(D)),prefilter::Q(1,mink(D)))-prefilter::m^2)",
    "prefilter::den(dot(prefilter::Q(1,mink(D)),prefilter::Q(2,mink(D))))^-2",
    // Colour generators, structure constants and existing colour words.
    "t(coad(8,a1),cof(3,c1),dind(cof(3,c2)))",
    "t(coad(8,a2),cof(3,c2),dind(cof(3,c3)))",
    "t(coad(8,a1),cof(3,c3),dind(cof(3,c1)))",
    "t(coad(8,a1),cof(Nc,c1),dind(cof(Nc,c2)))",
    "f(coad(8,a1),coad(8,a2),coad(8,a3))",
    "f(coad(8,a3),coad(8,a2),coad(8,a4))",
    "d(coad(8,a1),coad(8,a2),coad(8,a3))",
    "g(coad(8,a3),coad(8,a4))",
    "g(cof(3,c1),dind(cof(3,c3)))",
    "chain(cof(3,c1),dind(cof(3,c2)),t(coad(8,a1),in,out))",
    "prefilter::sextet(cos(6,s1),dind(cos(6,s2)),cof(3,c1))",
    // Heads whose arguments the matcher permutes, and ordered heads. Both
    // implementations hang in the matcher on a cyclesymmetric head with an
    // endpoint pair, so those are checked without it below.
    "prefilter::sym(dind(cof(3,c2)),cof(3,c1))",
    "prefilter::sym(cof(3,c1),prefilter::x,dind(cof(3,c2)))",
    "prefilter::antisym(dind(cof(3,c2)),cof(3,c1),prefilter::x)",
    "prefilter::sym(bis(4,i1),mink(4,mu1),bis(4,i2))",
    "prefilter::sym(prefilter::p(mink(4)),prefilter::q(mink(4)))",
    "prefilter::plain(dind(cof(3,c2)),cof(3,c1))",
    "prefilter::plain(prefilter::x,cof(3,c1),prefilter::y,dind(cof(3,c2)),prefilter::z)",
    "prefilter::plain(cof(3,c1),dind(cof(4,c2)))",
    "prefilter::plain(cof(3,c1),cof(3,c2))",
    "prefilter::plain(in,bis(4,i1),bis(4,i2))",
    "prefilter::plain(bis(4,i1),out,bis(4,i2))",
    "prefilter::plain(prefilter::plain(cof(3,c1),dind(cof(3,c2))),bis(4,i3),bis(4,i4))",
    "prefilter::plain(cof(3),dind(cof(3)))",
    "prefilter::plain(bis(4),bis(4,i2))",
    "prefilter::plain(prefilter::p(mink(4)),prefilter::q(mink(4)))",
    "prefilter::plain(prefilter::p(mink(4)),mink(4,mu1))",
    "prefilter::scalar(gamma(bis(4,i1),bis(4,i2),mink(4,mu1)),bis(4,i3),bis(4,i4))",
    // Scalars.
    "prefilter::x+prefilter::y",
    "(prefilter::x+prefilter::y)^3",
    "prefilter::m^2",
];

/// Deterministic products, sums, powers and nested functions of the factors.
fn generated_corpus() -> Vec<Atom> {
    let factors = FACTORS
        .iter()
        .map(|source| parse(source))
        .collect::<Vec<_>>();
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = |bound: usize| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 33) as usize) % bound
    };
    let mut corpus = factors.clone();
    for _ in 0..600 {
        let length = 1 + next(4);
        let mut product = Atom::num(1);
        for _ in 0..length {
            product *= &factors[next(factors.len())];
        }
        let expression = match next(6) {
            0 => &product + &factors[next(factors.len())] * &factors[next(factors.len())],
            1 => product.pow(Atom::num(2)),
            2 => symbolica::function!(
                symbol!("prefilter::plain"),
                &product,
                &factors[next(factors.len())]
            ),
            3 => {
                symbolica::function!(symbol!("prefilter::scalar"), &product)
                    * &factors[next(factors.len())]
            }
            _ => product,
        };
        corpus.push(expression);
    }
    corpus
}

/// Inputs of the stored gammaloop workflows, split into their summands.
fn workflow_inputs() -> Vec<(&'static str, Atom)> {
    let sources = [
        (
            "aa_aa_2l_gl16",
            include_str!(
                "../../../../tests/fixtures/aa_aa_2l_gl16_integrated_uv_start_after_simplify_metrics.sym"
            ),
        ),
        (
            "aa_aa_2l_gl16_odd_trace",
            include_str!("../../../../tests/fixtures/aa_aa_2l_gl16_odd_trace_term.sym"),
        ),
        (
            "massive_fermion_ladder",
            include_str!("../../../../benches/data/massive-fermion-ladder.expr"),
        ),
        (
            "fk0032_projected",
            include_str!("../../../../benches/data/fk0032-projected.expr"),
        ),
    ];
    sources
        .into_iter()
        .map(|(name, source)| (name, parse_workflow(source)))
        .collect()
}

fn summands(expression: &Atom) -> Vec<Atom> {
    match expression.as_view() {
        AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect(),
        _ => vec![expression.clone()],
    }
}

fn assert_identical(expression: &Atom, representation: LibraryRep) {
    let chained = expression.chainify(representation);
    let expected = super::chainify(&expression.as_view(), representation);
    assert_eq!(
        chained.as_view().get_data(),
        expected.as_view().get_data(),
        "{representation:?}: {expression}"
    );
}

/// The pre-checks admit every binding the matcher can make: start before end
/// for ordered heads, either order for heads whose arguments it permutes,
/// whatever order normalization leaves. They are checked directly here for
/// cyclesymmetric heads too, where the matcher itself does not terminate.
#[test]
fn pre_checks_admit_every_binding_the_matcher_can_make() {
    use crate::shorthands::chain::{chainify_rule_may_match, has_chainify_endpoint_pair};
    use spenso::structure::slot::SlotMatcher;
    use symbolica::atom::representation::FunView;

    fn function_of(atom: &Atom) -> FunView<'_> {
        match atom.as_view() {
            AtomView::Fun(function) => function,
            AtomView::Mul(product) => product
                .iter()
                .find_map(|factor| match factor {
                    AtomView::Fun(function) => Some(function),
                    _ => None,
                })
                .unwrap(),
            other => panic!("no function in {other}"),
        }
    }
    initialize();
    let fundamental = LibraryRep::from(ColorFundamental {});
    let slots = &mut SlotMatcher::default();
    let mut pair =
        |source: &str| has_chainify_endpoint_pair(function_of(&parse(source)), fundamental, slots);
    for source in [
        "prefilter::plain(cof(3,c1),dind(cof(3,c2)))",
        "prefilter::plain(prefilter::x,cof(3,c1),prefilter::y,dind(cof(3,c2)))",
        "prefilter::plain(cof(3),dind(cof(3)))",
        "prefilter::plain(prefilter::K(1,cof(3)),dind(cof(3,c2)))",
    ] {
        assert!(pair(source), "{source}");
    }
    for source in [
        "prefilter::plain(dind(cof(3,c2)),cof(3,c1))",
        "prefilter::plain(cof(3,c1),dind(cof(4,c2)))",
        "prefilter::plain(cof(3,c1),cof(3,c2))",
        "prefilter::plain(cof(3,c1),prefilter::x)",
    ] {
        assert!(!pair(source), "{source}");
    }
    for head in ["sym", "antisym", "cyclic"] {
        for (first, second) in [
            ("cof(3,c1)", "dind(cof(3,c2))"),
            ("dind(cof(3,c2))", "cof(3,c1)"),
        ] {
            let source = format!("prefilter::{head}({first},prefilter::y,{second})");
            assert!(pair(&source), "{source}");
        }
    }
    // The rule's head filter.
    for source in [
        "g(cof(3,c1),dind(cof(3,c2)))",
        "chain(cof(3,c1),dind(cof(3,c2)))",
        "prefilter::scalar(cof(3,c1),dind(cof(3,c2)))",
    ] {
        assert!(
            !chainify_rule_may_match(parse(source).as_view(), fundamental, slots),
            "{source}"
        );
    }
    assert!(chainify_rule_may_match(
        parse("prefilter::plain(cof(3,c1),dind(cof(3,c2)))").as_view(),
        fundamental,
        slots
    ));
}

#[test]
fn pre_checked_chainify_matches_reference_on_generated_corpus() {
    initialize();
    assert_eq!(parse("dot(x,y)").as_fun_view().unwrap().get_symbol(), T.dot);
    assert_eq!(
        parse("g(x,y)").as_fun_view().unwrap().get_symbol(),
        ETS.metric
    );
    assert_eq!(parse("t(x,y)").as_fun_view().unwrap().get_symbol(), CS.t);
    let corpus = generated_corpus();
    let checks = cross_checks();
    take_counts();
    let mut changed = 0;
    for expression in &corpus {
        for representation in representations() {
            assert_identical(expression, representation);
            changed += usize::from(expression.chainify(representation) != *expression);
        }
    }
    let counts = take_counts();
    assert_eq!(
        cross_checks() - checks,
        2 * corpus.len() * representations().len()
    );
    // The corpus exercises both outcomes of both pre-checks and real rewrites.
    assert!(counts.rule.iter().all(|&count| count > 0), "{counts:?}");
    assert!(
        counts.selection.iter().all(|&count| count > 0),
        "{counts:?}"
    );
    assert!(changed > 100, "{changed}");
}

#[test]
fn pre_checked_chainify_matches_reference_in_chain_collection() {
    initialize();
    let corpus = generated_corpus();
    let checks = cross_checks();
    for expression in corpus.iter().step_by(3) {
        for representation in representations() {
            // Each collection chainifies its input; the hook compares it.
            for (traces, retain_singletons) in [(true, true), (true, false), (false, false)] {
                let _ = expression.collect_chains(representation, true, traces, retain_singletons);
            }
            let _ = expression.join_chains(representation);
        }
    }
    assert!(cross_checks() > checks);
}

#[test]
fn pre_checked_chainify_matches_reference_on_workflow_summands() {
    initialize();
    let checks = cross_checks();
    let mut summand_count = 0;
    for (_, expression) in workflow_inputs() {
        for summand in summands(&expression) {
            summand_count += 1;
            for representation in representations() {
                assert_identical(&summand, representation);
            }
        }
    }
    assert!(summand_count > 64, "{summand_count}");
    assert_eq!(
        cross_checks() - checks,
        summand_count * representations().len()
    );
}

/// R* simplification as `RStarIntegrand::simplify` runs it, followed by the
/// colour stage of `RStarIntegrand::finish_tensor`. Both collect matrix
/// notation through chainify.
fn simplify_like_rstar(expression: &Atom) -> Result<Atom, String> {
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let gamma = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::canonical()),
        epsilon: true,
        ..Default::default()
    };
    let error = |error: &dyn std::fmt::Display| error.to_string();
    let cooked = cooking
        .try_cook(expression.as_view())
        .map_err(|e| error(&e))?;
    let simplified = SymbolicTensor::infer_indexed(cooked)
        .map_err(|e| error(&e))?
        .simplify_algebra(&gamma)
        .map_err(|e| error(&e))?
        .to_dots()
        .map_err(|e| error(&e))?;
    let simplified = cooking.uncook(simplified.expression().as_view());
    let cooked = cooking
        .try_cook(simplified.as_view())
        .map_err(|e| error(&e))?;
    let collected = SymbolicTensor::infer_indexed(cooked)
        .map_err(|e| error(&e))?
        .contract(ContractSettings::default())
        .map_err(|e| error(&e))?
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            ..gamma
        })
        .map_err(|e| error(&e))?;
    Ok(cooking.uncook(collected.expression().as_view()))
}

fn simplification_inputs() -> Vec<(&'static str, Vec<Atom>)> {
    workflow_inputs()
        .into_iter()
        .map(|(name, expression)| {
            let terms = if name == "aa_aa_2l_gl16" {
                // A representative slice of the 64 closed nine-gamma traces.
                summands(&expression).into_iter().step_by(8).collect()
            } else {
                vec![expression]
            };
            (name, terms)
        })
        .collect()
}

#[test]
fn pre_checked_chainify_matches_reference_in_rstar_like_simplification() {
    initialize();
    let checks = cross_checks();
    for (name, terms) in simplification_inputs() {
        for term in terms {
            if let Err(error) = simplify_like_rstar(&term) {
                panic!("{name}: {error}");
            }
        }
    }
    // Not every reference case has explicitly indexed free ports.
    let simplified = crate::reference_cases::reference_cases()
        .iter()
        .filter(|case| simplify_like_rstar(&case.expression()).is_ok())
        .count();
    assert!(simplified > 0);
    assert!(cross_checks() - checks > 100, "{}", cross_checks() - checks);
}

fn time<R>(repetitions: usize, mut body: impl FnMut() -> R) -> Duration {
    let start = Instant::now();
    for _ in 0..repetitions {
        std::hint::black_box(body());
    }
    start.elapsed() / repetitions as u32
}

/// Times both implementations on the chainify inputs that R*-like
/// simplification of the stored workflows produces, then the whole
/// simplification with either implementation inside it.
/// `cargo test -p idenso --release chainify_pre_check_benchmark -- --ignored --nocapture`
#[test]
#[ignore = "timing benchmark"]
fn chainify_pre_check_benchmark() {
    initialize();
    let mut total = [Duration::ZERO; 4];
    for (name, terms) in simplification_inputs() {
        let (_, inputs) = with_mode(Mode::PreChecked, || {
            record(|| {
                for term in &terms {
                    let _ = simplify_like_rstar(term).unwrap();
                }
            })
        });
        let repetitions = 5;
        let reference = time(repetitions, || {
            for (input, representation) in &inputs {
                let _ = std::hint::black_box(super::chainify(&input.as_view(), *representation));
            }
        });
        take_counts();
        let pre_checked = with_mode(Mode::PreChecked, || {
            time(repetitions, || {
                for (input, representation) in &inputs {
                    let _ = std::hint::black_box(input.chainify(*representation));
                }
            })
        });
        let counts = take_counts();
        let per_run = |count: usize| count / repetitions;
        let MatcherCounts { rule, selection } = counts;
        let pipeline_reference = time(3, || {
            with_mode(Mode::Reference, || {
                terms.iter().map(simplify_like_rstar).count()
            })
        });
        let pipeline_pre_checked = time(3, || {
            with_mode(Mode::PreChecked, || {
                terms.iter().map(simplify_like_rstar).count()
            })
        });
        println!(
            "{name}: {} chainify inputs ({} bytes); chainify {reference:?} -> {pre_checked:?} \
             ({:.2}x); rule matcher skipped {} of {}, selection skipped {} of {}; \
             simplification {pipeline_reference:?} -> {pipeline_pre_checked:?} ({:.3}x)",
            inputs.len(),
            inputs
                .iter()
                .map(|(input, _)| input.as_view().get_data().len())
                .sum::<usize>(),
            reference.as_secs_f64() / pre_checked.as_secs_f64(),
            per_run(rule[0]),
            per_run(rule[0] + rule[1]),
            per_run(selection[0]),
            per_run(selection[0] + selection[1]),
            pipeline_reference.as_secs_f64() / pipeline_pre_checked.as_secs_f64(),
        );
        for (sum, value) in total.iter_mut().zip([
            reference,
            pre_checked,
            pipeline_reference,
            pipeline_pre_checked,
        ]) {
            *sum += value;
        }
    }
    println!(
        "total: chainify {:?} -> {:?} ({:.2}x); simplification {:?} -> {:?} ({:.3}x)",
        total[0],
        total[1],
        total[0].as_secs_f64() / total[1].as_secs_f64(),
        total[2],
        total[3],
        total[2].as_secs_f64() / total[3].as_secs_f64(),
    );
}
