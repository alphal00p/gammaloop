use idenso::dirac::{AGS, GammaSimplifier, GammaSimplifySettings};
use spenso::network::tags::SPENSO_TAG;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, symbol,
};

fn parse(source: &str) -> Atom {
    Atom::parse(
        source,
        "spenso",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap()
}
fn initialize() {
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
}
fn trace(arguments: impl IntoIterator<Item = Atom>) -> Atom {
    spenso::trace!(parse("bis(4)"); arguments.into_iter().map(|a| {
        function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, a)
    }))
}
fn assert_expanded(input: &Atom, expected: &Atom) {
    let settings = GammaSimplifySettings::default().with_expanded_traces();
    let output = input.simplify_gamma_with(settings);
    assert_eq!(&output, expected);
    assert_eq!(output.simplify_gamma_with(settings), output);
    assert_eq!(
        input.simplify_gamma_with(settings.without_trace_evaluation()),
        input.simplify_gamma_with(GammaSimplifySettings::default().without_trace_evaluation()),
    );
}

#[test]
fn trace_expansion_preserves_spectators_and_independent_bodies() {
    initialize();
    let spectator = parse("(expanded_api::x+expanded_api::y)^8");
    for dimension in ["4", "D"] {
        let a = trace(
            ["mu", "a", "b", "mu", "c", "d", "e", "f"]
                .map(|i| parse(&format!("mink({dimension},expanded_api::{i})"))),
        );
        let b = trace(
            ["v0", "v1", "v2", "v3", "v4", "v5"]
                .map(|i| parse(&format!("mink({dimension},expanded_api::{i})"))),
        );
        let ea = a.simplify_gamma().expand();
        let eb = b.simplify_gamma().expand();
        assert_expanded(&a, &ea);
        assert_expanded(&(&spectator * &a), &(&spectator * &ea));
        assert_expanded(&(&spectator * &a * &b), &(&spectator * &ea * &eb));
    }
}

#[test]
fn four_dimensional_component_callbacks_are_expanded() {
    initialize();
    let vector = spenso::vector_symbol!(
        "expanded_api::sum_callback",
        norm = |value, out| {
            if matches!(value, AtomView::Fun(v) if matches!(v.iter().last(), Some(AtomView::Fun(slot)) if slot.get_nargs() == 2))
            {
                **out = parse("expanded_api::x+expanded_api::y");
            }
        }
    );
    let q = SPENSO_TAG.rank_one_tensor_symbol("expanded_api::q");
    let input = trace([
        parse("mink(4,a)"),
        function!(vector, parse("mink(4)")),
        parse("mink(4,b)"),
        function!(q, parse("mink(4)")),
    ]);
    assert_expanded(&input, &input.simplify_gamma().expand());
}

fn inner_trace() -> Atom {
    trace((0..6).map(|i| {
        let head = SPENSO_TAG.rank_one_tensor_symbol(&format!("expanded_api::inner_p{i}"));
        function!(head, parse("mink(D)"))
    }))
}
fn scalar(value: Atom) -> Atom {
    function!(symbol!("expanded_api::F"; Scalar), value)
}

#[test]
fn nested_traces_in_callbacks_and_metadata_keep_their_own_boundary() {
    initialize();
    let inner = inner_trace();
    let expanded = inner.simplify_gamma().expand();
    assert_ne!(expanded, inner.simplify_gamma());
    let vector = spenso::vector_symbol!(
        "expanded_api::trace_callback",
        norm = |value, out| {
            if matches!(value, AtomView::Fun(v) if matches!(v.iter().last(), Some(AtomView::Fun(slot)) if slot.get_nargs() == 2))
            {
                **out = scalar(inner_trace());
            }
        }
    );
    let q = SPENSO_TAG.rank_one_tensor_symbol("expanded_api::nested_q");
    let compact = function!(q, parse("mink(D)"));
    let callback = trace([
        parse("mink(D,a)"),
        function!(vector, parse("mink(D)")),
        parse("mink(D,b)"),
        compact.clone(),
    ]);
    let metadata = trace([
        function!(q, scalar(inner.clone()), parse("mink(D)")),
        compact.clone(),
    ]);
    let unit = spenso::trace!(function!(symbol!("spenso::bis"), scalar(inner.clone())); [
        function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, compact.clone()),
        function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, compact),
    ]);
    let opaque = |value: Atom| {
        spenso::trace!(parse("bis(4)"); [
            function!(symbol!("expanded_api::OpaqueMatrix"), SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out, scalar(value)),
        ])
    };
    let spectator = parse("(expanded_api::s+expanded_api::t)^8");
    for input in [callback, metadata, unit] {
        let expected = input
            .simplify_gamma()
            .replace(inner.to_pattern())
            .with(expanded.clone())
            .expand();
        assert_expanded(&input, &expected);
        assert_expanded(&(&spectator * input), &(&spectator * expected));
    }
    assert_expanded(&opaque(inner), &opaque(expanded));
}

#[test]
fn rewritten_and_collected_traces_expand_the_complete_body() {
    initialize();
    let spectator = parse("(expanded_api::s+expanded_api::t)^8");
    for heads in [
        vec![AGS.projp],
        vec![AGS.gamma5],
        vec![AGS.gamma5, AGS.gamma5],
    ] {
        let mut factors: Vec<_> = heads
            .into_iter()
            .map(|head| function!(head, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out))
            .collect();
        factors.extend((0..6).map(|i| {
            function!(
                AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                parse(&format!("mink(4,expanded_api::v{i})"))
            )
        }));
        let input = spenso::trace!(parse("bis(4)"); factors);
        let expected = input.simplify_gamma().expand();
        assert_expanded(&(&spectator * input), &(&spectator * expected));
    }
    let closed = Atom::mul_many((0..4).map(|i| {
        function!(
            AGS.gamma,
            parse(&format!("bis(4,expanded_api::s{i})")),
            parse(&format!("bis(4,expanded_api::s{})", (i + 1) % 4)),
            parse(&format!("mink(D,expanded_api::v{i})"))
        )
    }));
    let expected = closed.simplify_gamma().expand();
    assert_expanded(&(&spectator * closed), &(&spectator * expected));
}

#[test]
fn repeated_four_dimensional_indices_preserve_expanded_scalar_context() {
    initialize();
    let p_head = SPENSO_TAG.rank_one_tensor_symbol("expanded_four::p");
    let q_head = SPENSO_TAG.rank_one_tensor_symbol("expanded_four::q");
    let p = function!(p_head, parse("mink(4)"));
    let q = function!(q_head, parse("mink(4)"));
    let mu = parse("mink(4,expanded_four::mu)");
    let dot = spenso::g!(&p, &q);
    let spectator = parse("(expanded_four::s+expanded_four::t)^8");
    // gamma(mu) p/ gamma(mu) = -2 p/ in four dimensions.
    let odd = trace([mu.clone(), p.clone(), mu.clone(), q.clone()]);
    // gamma(mu) p/ q/ gamma(mu) = 4(p.q) times the identity.
    let even = trace([mu.clone(), p.clone(), q.clone(), mu, p, q]);
    for (input, expected) in [
        (odd, Atom::num(-8) * &dot),
        (even, Atom::num(16) * &dot * &dot),
    ] {
        assert_expanded(&input, &expected);
        let input = &spectator * input;
        let expected = &spectator * expected;
        assert_expanded(&input, &expected);
        assert!(matches!(expected.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == spectator.as_view())));
    }
}
