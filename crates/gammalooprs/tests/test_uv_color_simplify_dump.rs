use std::{path::PathBuf, time::Instant};

use gammalooprs::initialisation::test_initialise;
use idenso::color::ColorSimplifier;
use symbolica::parse;

fn fixture_path(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/resources/uv_color_simplify")
        .join(file_name)
}

#[test]
#[ignore = "large UV scalar hotspot fixture; run explicitly when optimizing color simplification"]
fn loads_uv_scalar_profile_hotspot_dump_and_simplifies_color() {
    test_initialise().unwrap();

    let path = fixture_path("banana_topo29_dod4_term0.sym");
    let source = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));

    assert!(
        source.len() > 500_000,
        "expected the extracted dod4 S_GS⊛44 hotspot fixture, got only {} bytes from {}",
        source.len(),
        path.display()
    );

    let parse_started = Instant::now();
    let expr = parse!(&source);
    eprintln!(
        "parsed {} bytes from {} in {:?}",
        source.len(),
        path.display(),
        parse_started.elapsed()
    );

    let simplify_started = Instant::now();
    let simplified = expr.simplify_color();
    eprintln!(
        "color-simplified fixture in {:?}",
        simplify_started.elapsed()
    );

    assert!(
        !simplified.is_zero(),
        "color simplification unexpectedly reduced the UV hotspot fixture to zero"
    );
}
