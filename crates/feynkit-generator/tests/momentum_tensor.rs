use symbolica::symbol;

#[test]
fn momentum_head_is_tagged_before_its_first_public_lookup() {
    // Link and initialize the generator crate without generating a diagram.
    let _run = feynkit_generator::Run::default();

    let momentum = symbol!("gammalooprs::Q");
    assert!(momentum.has_tag("spenso::tensor"));
    assert!(momentum.has_tag("spenso::rank1"));
}
