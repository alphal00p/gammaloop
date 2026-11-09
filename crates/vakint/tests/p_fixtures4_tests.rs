//! The frozen four-loop Forcer fixtures of the M3 plan (M3-V7), one parse test
//! per file: each file holds exactly the cases of its definition, carries its
//! provenance (Forcer commit, FORM version, run flags and seconds), and every
//! frozen output parses, has a tag-basis form and finalizes. With FORM 5,
//!
//! ```text
//! VAKINT_FORCER_FIXTURES=write VAKINT_FORCER_FIXTURES_WORKERS=4 \
//!     cargo nextest run -p vakint --test p_fixtures4_tests
//! ```
//!
//! regenerates the files (`tform -w4`, one FORM process at a time).

#[path = "support/p_fixtures.rs"]
mod p_fixtures;
#[path = "support/p_fixtures4.rs"]
mod p_fixtures4;

use p_fixtures::*;
use p_fixtures4::*;
use vakint::p_integral::compare::tag_basis;
use vakint::p_integral::forcer_fixtures::ForcerFixture;
use vakint::{LoopNormalizationFactor, RustRedMasterMode, VakintSettings};

/// The provenance block is complete and every case has a nonzero, parseable,
/// finalizable value (none of these cases is scaleless).
fn assert_frozen(file: &str, fixture: &ForcerFixture) {
    assert!(!fixture.form_version.is_empty(), "{file}: FORM version");
    let run = fixture.run.as_ref().expect("a run block");
    assert!(run.flags.contains("Forcer()"), "{file}: {run:?}");
    let settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let epsilon = vakint::vakint_symbol!(&settings.epsilon_symbol);
    for case in &fixture.cases {
        let values = frozen_values(case);
        assert!(!values.is_empty(), "{file}: {}", case.name);
        for (_, reduced) in &values {
            let basis = tag_basis(reduced, epsilon)
                .unwrap_or_else(|error| panic!("{file}: {}: {error}", case.name));
            assert!(!basis.is_empty(), "{file}: {}", case.name);
        }
        let value = finalized(case, &settings, RustRedMasterMode::Retain);
        assert!(!value.is_zero(), "{file}: {}", case.name);
    }
}

#[test]
fn the_four_loop_anchor_fixture_is_complete() {
    let fixture = anchor4_fixture();
    assert_frozen("p4l_anchors.json", &fixture);
    // A master anchor is the master integral itself: its Forcer value carries
    // its own master and no other.
    for case in &fixture.cases {
        if let Some(master) = &case.master {
            let output = case.outputs.concat();
            assert!(
                output.contains(&format!("Master({master})")),
                "{}",
                case.name
            );
            assert_eq!(output.matches("Master(").count(), 1, "{}", case.name);
        }
    }
    // la4_top is the battery's ladder, pinned to -70 zeta(7) by
    // forcer_four_loop_battery.
    let ladder = case(&fixture, "la4_top");
    let frozen = battery("p4l_battery.json");
    let pinned = case(&frozen, "la4_top");
    assert_eq!(
        (&ladder.family, &ladder.powers, &ladder.outputs),
        (&pinned.family, &pinned.powers, &pinned.outputs)
    );
}
