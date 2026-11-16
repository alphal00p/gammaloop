use spenso::network::library::symbolic::ETS;
use spenso::network::tags::SPENSO_TAG as T;
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::SlotMatcher;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::parser::ParseSettings;
use symbolica::{function, symbol};

use super::has_endpoint_pair;
use crate::color::CS;
use crate::dirac::AGS;
use crate::representations::{Bispinor, ColorAdjoint, ColorFundamental, ColorSextet};
use crate::shorthands::chain::chainify_with;
use crate::test_support::test_initialize;

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
    Atom::parse(source, "spenso", ParseSettings::symbolica())
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

fn representations() -> [LibraryRep; 7] {
    [
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
            2 => function!(
                symbol!("prefilter::plain"),
                &product,
                &factors[next(factors.len())]
            ),
            3 => function!(symbol!("prefilter::scalar"), &product) * &factors[next(factors.len())],
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

fn run<const PRECHECK: bool>(expression: &Atom, rep: LibraryRep) -> (Atom, Vec<(Atom, Atom)>) {
    let mut updates = Vec::new();
    let result =
        chainify_with::<PRECHECK>(expression.as_view(), rep, None, None, &mut |old, new| {
            updates.push((old.to_owned(), new.to_owned()));
        });
    (result, updates)
}

/// Both the rewritten expression and the sequence of observations used to
/// publish materialized interfaces must agree with the matcher-only path.
#[test]
fn chainify_prechecks_preserve_rewrites_and_observations() {
    initialize();
    for input in generated_corpus().into_iter().chain(
        workflow_inputs()
            .into_iter()
            .flat_map(|(_, input)| summands(&input)),
    ) {
        for rep in representations() {
            assert_eq!(
                run::<true>(&input, rep),
                run::<false>(&input, rep),
                "{rep:?}: {input}"
            );
        }
    }
}

#[test]
fn endpoint_precheck_respects_dimension_and_head_symmetry() {
    initialize();
    let rep = LibraryRep::from(ColorFundamental {});
    let mut slots = SlotMatcher::default();
    for (head, expected) in [
        ("plain", false),
        ("sym", true),
        ("antisym", true),
        ("cyclic", true),
    ] {
        let expr = parse(&format!(
            "prefilter::{head}(dind(cof(3,k)),prefilter::x,cof(3,j))"
        ));
        // Symmetric normalization may reorder the arguments or extract a sign.
        let mut observed = None;
        expr.visitor(&mut |node| {
            if let AtomView::Fun(f) = node
                && f.get_symbol().get_name() == format!("prefilter::{head}")
            {
                observed = Some(has_endpoint_pair(f, rep, &mut slots));
                return false;
            }
            true
        });
        assert_eq!(observed, Some(expected), "{expr}");
    }
    let unequal = parse("prefilter::plain(cof(3,j),dind(cof(4,k)))");
    assert!(!has_endpoint_pair(
        unequal.as_fun_view().unwrap(),
        rep,
        &mut slots
    ));
}
