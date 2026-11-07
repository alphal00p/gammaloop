#![cfg(feature = "symbolica")]

// Symbolica's symbol table is process global, and nextest runs every test in
// its own process: each test sees written dummy-shaped names before Spenso
// has defined the tagged generated symbol of the same spelling.

use spenso::network::parsing::{NetworkParse, ParseSettings, ParseState};
use spenso::structure::abstract_index::{AbstractIndex, DUMMY_INDEX_TAG};
use spenso::structure::representation::{Minkowski, RepName};
use spenso::structure::slot::{DummyAind, IsAbstractSlot, ParseableAind};
use spenso::{chain, chain_factor, mink, tensor, vector};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::try_parse;

fn tagged(atom: &Atom) -> bool {
    matches!(atom.as_view(), AtomView::Var(variable) if variable.get_symbol().has_tag(DUMMY_INDEX_TAG))
}

// A chain with a compact vector is lowered with the parser's first dummy.
fn lowered_chain() -> Atom {
    chain!(
        mink!(4, i),
        mink!(4, j),
        chain_factor!(factor, in, out, (vector!(compact_p, mink!(4))))
    )
}

#[test]
fn written_parser_dummy_name_does_not_block_lowering() {
    let written = try_parse!("spenso::d_1000000").unwrap();
    assert!(!tagged(&written));
    let explicit = tensor!(
        explicit,
        Minkowski {}.new_slot::<AbstractIndex, _, _>(
            4,
            AbstractIndex::try_from(written.as_view()).unwrap()
        )
    );
    let expression = lowered_chain() * &explicit;

    let network = expression
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    let dangling = network.graph.dangling_indices();
    assert!(
        dangling.iter().any(|slot| slot.aind.to_atom() == written),
        "the written index must remain external: {dangling:?}"
    );

    let generated = AbstractIndex::new_dummy_at(1_000_000).to_atom();
    assert!(tagged(&generated));
    assert_ne!(generated, written);
}

#[test]
fn indices_parsed_from_dummy_names_do_not_block_fresh_dummies() {
    let written = AbstractIndex::try_from("d_3").unwrap().to_atom();
    assert!(!tagged(&written));
    let state = ParseState::<AbstractIndex>::default();
    state.reserve_index(AbstractIndex::try_from(written.as_view()).unwrap());
    let fresh = (0..5)
        .map(|_| state.fresh_index().to_atom())
        .collect::<Vec<_>>();
    for (position, atom) in fresh.iter().enumerate() {
        assert!(tagged(atom), "{atom}");
        assert_ne!(atom, &written);
        assert!(!fresh[..position].contains(atom));
    }
    // A later lookup of the same counter is stable.
    assert_eq!(
        AbstractIndex::new_dummy_at(3).to_atom(),
        AbstractIndex::new_dummy_at(3).to_atom()
    );
}

#[test]
fn attributed_dummy_names_round_trip_through_printing_and_parsing() {
    // A canonical string carries the tag, so a fresh process re-creates the
    // generated symbol itself rather than a written one of the same spelling.
    let written = try_parse!("spenso::{spenso::dummy-index}::d_2000000").unwrap();
    assert!(tagged(&written));
    assert_eq!(AbstractIndex::new_dummy_at(2_000_000).to_atom(), written);

    let slot =
        Minkowski {}.new_slot::<AbstractIndex, _, _>(4, AbstractIndex::new_dummy_at(1_000_001));
    let original = tensor!(left, slot.to_atom()) * tensor!(right, slot.to_atom());
    let printed = original.to_canonical_string();
    assert!(printed.contains(DUMMY_INDEX_TAG), "{printed}");
    assert_eq!(try_parse!(&printed).unwrap(), original);
}
