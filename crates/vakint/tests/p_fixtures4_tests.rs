//! The frozen four-loop Forcer fixtures of the M3 plan (M3-V7), one parse test
//! per file: each file holds exactly the cases of its definition, carries its
//! provenance (Forcer commit, FORM version, run flags and seconds), and every
//! frozen output parses, has a tag-basis form and finalizes. With FORM 5,
//!
//! ```text
//! VAKINT_FORCER_FIXTURES=write VAKINT_FORCER_FIXTURES_WORKERS=4 \
//!     cargo nextest run -p vakint --test p_fixtures4_tests -E 'test(<one test>)'
//! ```
//!
//! regenerates the file of that test (`tform -w4`, one FORM process at a
//! time) and its depths; `VAKINT_FORCER_FIXTURES=depths` recomputes only the
//! depths. Select one test: `p4l_anchors.json` is frozen for M3-T17. The
//! ignored `every_four_loop_grid_depth_is_the_finalizers` recomputes every
//! grid depth (`--run-ignored all`), and the Forcer programs of the grid and
//! the anchors are pinned by digest.

#[path = "support/p_fixtures.rs"]
mod p_fixtures;
#[path = "support/p_fixtures4.rs"]
mod p_fixtures4;

use std::collections::BTreeSet;

use p_fixtures::*;
use p_fixtures4::*;
use vakint::p_integral::compare::tag_basis;
use vakint::p_integral::forcer_fixtures::{
    ANCHOR4_PROVENANCE, CaseDepth, DEPTH_CAP_TERMS, Finalized, ForcerCase, ForcerFixture,
    GRID4_PROVENANCE, I7_TERMS, anchor_cases_4l, grid_cases_4l,
};
use vakint::{LoopNormalizationFactor, RustRedMasterMode, VakintSettings};

/// The four-loop benchmark tops whose Forcer value is an exact 0: the J11
/// top with `dot(k(1),p(1))*dot(k(2),p(1))`, undotted and with one dot (raw
/// output ` 0`). The grid keeps them as data.
const ZERO_GRID_CASES: [&str; 2] = ["bench_L4_d0_r2_P4L_J11", "bench_L4_d1_r2_P4L_J11"];

/// The provenance block is complete and every case has a parseable,
/// finalizable value, nonzero except for the cases `zeros` (none of these
/// cases is scaleless).
fn assert_frozen(file: &str, fixture: &ForcerFixture, zeros: &[&str]) {
    assert!(!fixture.form_version.is_empty(), "{file}: FORM version");
    let run = fixture.run.as_ref().expect("a run block");
    assert!(run.flags.contains("Forcer()"), "{file}: {run:?}");
    let settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let epsilon = vakint::vakint_symbol!(&settings.epsilon_symbol);
    for case in fixture.cases.iter().chain(&fixture.boundary) {
        let values = frozen_values(case);
        if zeros.contains(&case.name.as_str()) {
            assert_eq!(case.outputs, [" 0\n"], "{file}: {}", case.name);
            assert!(
                values.iter().all(|(_, v)| v.is_zero()),
                "{file}: {}",
                case.name
            );
            continue;
        }
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

/// The depth of a case is consistent: `ε^(terms-L-1)`, within the cap, a
/// limit exactly when below the cap, and a nonzero value.
fn assert_depth<'a>(file: &str, case: &'a ForcerCase) -> &'a CaseDepth {
    let depth = case
        .depth
        .as_ref()
        .unwrap_or_else(|| panic!("{file}: {}: no depth", case.name));
    let loops = case.loops().unwrap() as i64;
    assert_eq!(
        depth.order,
        depth.terms - loops - 1,
        "{file}: {}",
        case.name
    );
    assert!(
        (1..=DEPTH_CAP_TERMS).contains(&depth.terms),
        "{file}: {}",
        case.name
    );
    assert_eq!(
        depth.limit.is_empty(),
        depth.terms == DEPTH_CAP_TERMS,
        "{file}: {}",
        case.name
    );
    assert_nonzero_series(&depth.symbolic, depth.terms, &case.name);
    depth
}

/// `terms` parseable unit-scale Laurent coefficients, not all zero.
fn assert_nonzero_series(coefficients: &[String], terms: i64, name: &str) {
    assert_eq!(coefficients.len() as i64, terms, "{name}");
    let parsed = coefficients.iter().map(|c| parse(c)).collect::<Vec<_>>();
    assert!(parsed.iter().any(|c| !c.is_zero()), "{name}");
}

#[test]
fn the_four_loop_anchor_fixture_is_complete() {
    let fixture = anchor4_fixture();
    assert_frozen("p4l_anchors.json", &fixture, &[]);
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

#[test]
fn the_four_loop_grid_fixture_is_complete() {
    let file = "p4l_grid.json";
    let fixture = grid4_fixture();
    assert_frozen(file, &fixture, &ZERO_GRID_CASES);
    let count = |prefix: &str| {
        fixture
            .cases
            .iter()
            .filter(|case| case.name.starts_with(prefix))
            .count()
    };
    // 15 cells x 11 families; 2 sectors at rank 0 and 2 and 1 dotted sector
    // per family.
    assert_eq!(
        [count("bench_L4_"), count("sector_"), count("dotted_")],
        [165, 44, 11]
    );
    assert_eq!(fixture.cases.len(), 220);
    for case in &fixture.cases {
        if ZERO_GRID_CASES.contains(&case.name.as_str()) {
            // A zero finalizes at any depth: the cap, every coefficient 0.
            let depth = case.depth.as_ref().unwrap();
            assert_eq!(depth.terms, DEPTH_CAP_TERMS, "{}", case.name);
            assert!(depth.symbolic.iter().all(|c| c == "0"), "{}", case.name);
        } else {
            assert_depth(file, case);
        }
    }
    // The depths are the finalizer's (recomputed here for the dotted sectors;
    // forcer_four_loop_battery recomputes every battery depth).
    for case in fixture
        .cases
        .iter()
        .filter(|case| case.name.starts_with("dotted_"))
    {
        assert_eq!(&case.deepest().unwrap(), assert_depth(file, case));
    }
}

#[test]
fn the_four_loop_battery_fixture_is_complete() {
    let file = "p4l_battery.json";
    let fixture = battery4_fixture();
    assert_frozen(file, &fixture, &[]);
    // All eleven four-loop families, with ranks, dots and ISPs.
    let families = fixture
        .cases
        .iter()
        .map(|case| case.family.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(families.len(), 11, "{families:?}");
    assert!(families.iter().all(|family| family.starts_with("P4L_")));
    for kind in ["_top", "_isp", "_rank3", "_dot"] {
        let with = fixture
            .cases
            .iter()
            .filter(|case| case.name.ends_with(kind))
            .map(|case| case.family.as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(with, families, "{kind}");
    }
    // Every case has its depth and its 11-term value; the Oep cases need an
    // ALARM order of no1, no2 or lala at ε^6, or of haha, no6, cross or nono
    // at ε^7.
    let mut alarms = BTreeSet::new();
    for case in &fixture.cases {
        let depth = assert_depth(file, case);
        match case.eleven_terms.as_ref().expect("an 11-term value") {
            Finalized::Symbolic(value) => {
                assert!(depth.terms >= I7_TERMS, "{}", case.name);
                assert_nonzero_series(value, I7_TERMS, &case.name);
            }
            Finalized::MissingMasterOrders(missing) => {
                assert!(depth.terms < I7_TERMS, "{}", case.name);
                alarms.extend(
                    missing
                        .iter()
                        .map(|master| (master.master.clone(), master.first_missing_order)),
                );
            }
        }
    }
    let alarm = |master: &str, order: i64| alarms.contains(&(format!("PMaster({master})"), order));
    assert!(
        ["no1", "no2", "lala"].iter().all(|master| alarm(master, 6)),
        "{alarms:?}"
    );
    assert!(
        ["haha", "no6", "cross", "nono"]
            .iter()
            .all(|master| alarm(master, 7)),
        "{alarms:?}"
    );
    // M2's three-loop Oep boundary: the NO top finalizes through 11 terms
    // (ε^7) and no further (Master(no) is known through ε^7).
    let [boundary] = fixture.boundary.as_slice() else {
        panic!("one boundary case")
    };
    let depth = assert_depth(file, boundary);
    assert_eq!((depth.terms, depth.limit.len()), (I7_TERMS, 1));
    assert_eq!(
        boundary.eleven_terms,
        Some(Finalized::Symbolic(depth.symbolic.clone()))
    );
    // The whole-family blind spots, with the per-integral cases Forcer runs:
    // exactly the cases whose output carries the master (FORM wraps lines
    // anywhere, so the text is compared without whitespace).
    let nodes = fixture
        .forcer_blind
        .iter()
        .map(|blind| blind.node.as_str())
        .collect::<Vec<_>>();
    assert_eq!(nodes, ["T147", "T152", "T153"]);
    for blind in &fixture.forcer_blind {
        let master = format!("Master({})", blind.master);
        let carriers = fixture
            .cases
            .iter()
            .filter(|case| {
                let output = case.outputs.concat();
                let compact = output
                    .chars()
                    .filter(|c| !c.is_whitespace())
                    .collect::<String>();
                compact.contains(&master)
            })
            .map(|case| case.name.as_str())
            .collect::<Vec<_>>();
        assert!(!carriers.is_empty(), "{}", blind.node);
        assert_eq!(blind.per_integral, carriers, "{}", blind.node);
    }
}

/// Every stored grid depth is the finalizer's: the parse test recomputes the
/// dotted sectors only, this recomputes all 220 cases from the frozen outputs
/// (no FORM; about 90 s with `VAKINT_FORCER_FIXTURES_WORKERS=4`).
#[test]
#[ignore = "recomputes all 220 four-loop grid depths (minutes)"]
fn every_four_loop_grid_depth_is_the_finalizers() {
    let fixture = grid4_fixture();
    let mut recomputed = fixture.clone();
    recomputed
        .write_depths(false, workers())
        .unwrap_or_else(|error| panic!("p4l_grid.json: {error}"));
    for (stored, fresh) in fixture.cases.iter().zip(&recomputed.cases) {
        assert_eq!(stored.depth, fresh.depth, "{}", stored.name);
    }
}

/// `blake3` of the Forcer programs (name, program) of the grid and anchor
/// definitions, as `forcer_programs_are_byte_identical` keeps for the other
/// fixture files: the frozen outputs belong to these programs.
const PROGRAM_DIGESTS_4L: [(&str, &str); 2] = [
    (
        "p4l_grid.json",
        "176becab393e704d529b02346a8e9eb837209b170c087bef79f1d8fac6009ebe",
    ),
    (
        "p4l_anchors.json",
        "1ee821498d2db0349bd0e57eadfe3598177ca33a91a2d6affcd4a49735230b7f",
    ),
];

#[test]
fn the_four_loop_grid_and_anchor_programs_are_byte_identical() {
    let found = [
        ForcerFixture::unwritten(
            GRID4_PROVENANCE,
            grid_cases_4l(&benchmark_inputs_text()).unwrap(),
        ),
        ForcerFixture::unwritten(ANCHOR4_PROVENANCE, anchor_cases_4l()),
    ]
    .map(|fixture| fixture.program_digest().unwrap());
    assert_eq!(
        found,
        PROGRAM_DIGESTS_4L.map(|(_, digest)| digest.to_string()),
        "the Forcer programs of {:?} changed",
        PROGRAM_DIGESTS_4L.map(|(file, _)| file)
    );
}
