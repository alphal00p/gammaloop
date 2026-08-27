use spenso::network::tags::SPENSO_TAG;
use symbolica::symbol;

#[test]
fn momentum_head_is_tagged_before_its_first_public_lookup() {
    // Link and initialize the generator crate without generating a diagram.
    let _options = feynkit_generator::GenerationOptions::default();

    let momentum = symbol!("FeynKit::Momentum");
    assert!(momentum.has_tag(&SPENSO_TAG.tensor));
    assert!(momentum.has_tag(&SPENSO_TAG.rank1));
}
