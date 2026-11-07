#![cfg(feature = "symbolica")]

// The parser numbers its dummies from 1_000_000. This test needs the tagged
// `d_1000000` to be created first, by parsing its canonical string. Symbolica's
// symbol table is process global and `dummy_symbol_collisions` writes an
// untagged `d_1000000`, so this case lives in its own test binary: both
// `cargo test` and nextest run it in a fresh process.

use spenso::network::parsing::{NetworkParse, ParseSettings};
use spenso::structure::abstract_index::{AbstractIndex, DUMMY_INDEX_TAG};
use spenso::structure::representation::{Minkowski, RepName};
use spenso::structure::slot::{DummyAind, ParseableAind};
use spenso::{chain, chain_factor, mink, tensor, vector};
use symbolica::atom::{Atom, AtomView};
use symbolica::try_parse;

fn tagged(atom: &Atom) -> bool {
    matches!(atom.as_view(), AtomView::Var(variable) if variable.get_symbol().has_tag(DUMMY_INDEX_TAG))
}

#[test]
fn lowering_follows_a_parsed_tagged_first_parser_dummy() {
    let parsed = try_parse!("spenso::{spenso::dummy-index}::d_1000000").unwrap();
    assert!(tagged(&parsed));
    assert_eq!(AbstractIndex::new_dummy_at(1_000_000).to_atom(), parsed);
    // The parser's second dummy spelling is also written, without the tag.
    let written = try_parse!("spenso::d_1000001").unwrap();
    assert!(!tagged(&written));

    let explicit = tensor!(
        explicit,
        Minkowski {}.new_slot::<AbstractIndex, _, _>(
            4,
            AbstractIndex::try_from(written.as_view()).unwrap()
        )
    );
    let lowered = chain!(
        mink!(4, i),
        mink!(4, j),
        chain_factor!(first, in, out, (vector!(compact_p, mink!(4)))),
        chain_factor!(second, in, out, (vector!(compact_q, mink!(4))))
    ) * &explicit;
    let network = lowered
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    let dangling = network.graph.dangling_indices();
    assert!(
        dangling.iter().any(|slot| slot.aind.to_atom() == written),
        "the written index must remain external: {dangling:?}"
    );
    assert_eq!(AbstractIndex::new_dummy_at(1_000_000).to_atom(), parsed);
    let generated = AbstractIndex::new_dummy_at(1_000_001).to_atom();
    assert!(tagged(&generated));
    assert_ne!(generated, written);
}
