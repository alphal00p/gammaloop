//! The frozen Forcer grid fixtures `tests/fixtures/forcer/p{1,2,3}l_grid.json`
//! and `p_anchors.json`. With `tests/fixtures/benchmark_p_integrals.json` they
//! total about 0.38 MB and these checks take about a second, so they all run in
//! the fast gate (`just test_vakint`).
//!
//! Each file must hold exactly the cases of its definition
//! (`vakint::p_integral::forcer_fixtures::{grid_cases, anchor_cases}`), and
//! every frozen output must parse, have a tag-basis form
//! (`vakint::p_integral::compare::tag_basis`) and finalize. With FORM 5,
//!
//! ```text
//! VAKINT_FORCER_FIXTURES=write VAKINT_FORCER_FIXTURES_JOBS=4 \
//!     cargo nextest run -p vakint --test p_integral_forcer_fixtures
//! ```
//!
//! regenerates the files (the anchor values are checked by the unit tests in
//! `src/p_integral/forcer/tests.rs`).

#[path = "support/p_fixtures.rs"]
mod p_fixtures;

use p_fixtures::*;
use vakint::p_integral::compare::tag_basis;
use vakint::p_integral::forcer_fixtures::ForcerFixture;
use vakint::{LoopNormalizationFactor, RustRedMasterMode, VakintSettings};

fn retain_settings() -> VakintSettings {
    VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    }
}

/// Every case parses, has a tag-basis form and finalizes exactly to a nonzero
/// value; scaleless cases have no output. A sector whose fixed rank-2 numerator
/// gives an exact 0 is excluded by construction (`fixtures::grid_cases`), so a
/// zero here is a change in Forcer's output, or its silent zero.
fn assert_complete(file: &str, fixture: &ForcerFixture) {
    let settings = retain_settings();
    let epsilon = vakint::vakint_symbol!(&settings.epsilon_symbol);
    for case in &fixture.cases {
        let values = frozen_values(case);
        assert_eq!(
            values.is_empty(),
            case.outputs.is_empty(),
            "{file}: {}",
            case.name
        );
        for (_, reduced) in &values {
            let basis = tag_basis(reduced, epsilon)
                .unwrap_or_else(|error| panic!("{file}: {}: {error}", case.name));
            assert_eq!(basis.is_empty(), reduced.is_zero(), "{file}: {}", case.name);
        }
        let value = finalized(case, &settings, RustRedMasterMode::Retain);
        assert_eq!(value.is_zero(), values.is_empty(), "{file}: {}", case.name);
    }
}

#[test]
fn the_anchor_fixture_is_complete() {
    let fixture = anchor_fixture();
    assert_complete("p_anchors.json", &fixture);
}

#[test]
fn the_one_loop_grid_fixture_is_complete() {
    assert_complete("p1l_grid.json", &grid_fixture(1));
}

#[test]
fn the_two_loop_grid_fixture_is_complete() {
    assert_complete("p2l_grid.json", &grid_fixture(2));
}

#[test]
fn the_three_loop_grid_fixture_is_complete() {
    assert_complete("p3l_grid.json", &grid_fixture(3));
}
