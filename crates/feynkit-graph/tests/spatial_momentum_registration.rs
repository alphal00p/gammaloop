use feynkit_graph::symbols;
use symbolica::atom::{Atom, AtomView};
use symbolica::parser::ParseSettings;

fn parse(expression: &str) -> Atom {
    Atom::parse(expression, "test", ParseSettings::default()).unwrap()
}

// A fresh binary ensures no earlier getter call registers Q3 for this test: the
// head must already carry its normalizer when a card or DOT string is parsed.
#[test]
fn standalone_parsing_normalizes_spatial_momentum_without_a_getter() {
    let time = parse("gammalooprs::Q3(0,spenso::cind(0))");
    assert!(time.is_zero(), "{time}");
    let spatial = parse("gammalooprs::Q3(0,spenso::cind(2))");
    assert_eq!(spatial, parse("gammalooprs::Q(0,spenso::cind(2))"));
    // A later getter call returns the parsed head instead of redefining it.
    let indexed = parse("gammalooprs::Q3(0,spenso::mink(4,7))");
    let AtomView::Fun(call) = indexed.as_view() else {
        panic!("an indexed Q3 stays a function: {indexed}");
    };
    assert_eq!(call.get_symbol(), symbols::spatial_momentum());
}
