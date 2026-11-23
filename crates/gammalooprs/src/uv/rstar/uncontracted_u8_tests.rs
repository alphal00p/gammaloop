//! U8 (results/consolidation/U8-DESIGN.md of the beta project): the switch's
//! source and the default kinds, the drop / demotion of a DEFAULT switch by
//! an incompatible setting (a setting the user wrote wins over a default;
//! every written value keeps its refusals), and one identity for a default
//! run and the explicit list of the default kinds.
//!
//! The tests marked "own process" set process environment variables: run
//! them under cargo nextest (one process per test), as the existing
//! start-up refusal tests.

use eyre::{Result, ensure, eyre};

use super::{
    top_mass::IrQuotient,
    uncontracted::{
        self, DEFAULT_KINDS, Kind, Mode, Source, Switch, TEST_DEFAULT_VARIABLE, UNSET_IS_ON,
        unset_is_on,
        k5::{HoleCut, with_hole_cut},
        with_switch,
    },
    uncontracted_ledger::{self as ledger, Exceptions, Routes},
};

fn parse(value: Option<&str>) -> Switch {
    Switch::parse(value).unwrap_or_else(|error| panic!("{value:?}: {error}"))
}

fn init(switch: Switch, routes: Routes, pint6: bool) -> Result<ledger::Initialised, String> {
    ledger::initialise(switch, None, None, || routes, || pint6, None)
}

fn routes(uv_children_rstar: bool, ir_hard_child: bool) -> Routes {
    Routes {
        uv_children_rstar,
        ir_hard_child,
        ..Routes::PRODUCTION
    }
}

/// The default kinds are root, child and ir-hard (one list for the unset
/// switch and a bare `on` / `check`); surrogate must be listed (ruling
/// 18:28Z).
#[test]
fn default_kinds_are_root_child_ir_hard() {
    assert_eq!(
        DEFAULT_KINDS,
        [Kind::Root, Kind::Child, Kind::IrHard]
    );
    assert!(!Kind::Surrogate.in_default() && Kind::Surrogate.implemented());
    for kind in Kind::ALL {
        assert_eq!(kind.in_default(), DEFAULT_KINDS.contains(&kind), "{kind:?}");
    }
    assert_eq!(Switch::DEFAULT_ON.enabled_kinds(), DEFAULT_KINDS.to_vec());
    assert_eq!(parse(Some("on")).enabled_kinds(), DEFAULT_KINDS.to_vec());
    assert_eq!(parse(Some("check")).enabled_kinds(), DEFAULT_KINDS.to_vec());
    assert!(Switch::DEFAULT_ON.inert_kinds().is_empty());
}

/// Only a `:<kinds>` list makes kinds explicit (ruling D1); the unset
/// switch is the default source, a bare value the bare one.
#[test]
fn switch_source_explicit_bits() {
    let default = Switch::DEFAULT_ON;
    assert_eq!(default.mode, Mode::On);
    assert_eq!(default.source(), Source::Default);
    assert!(default.is_default() && default.explicit_kinds().is_empty());
    for bare in ["on", "1", "check", " on "] {
        let switch = parse(Some(bare));
        assert_eq!(switch.source(), Source::Bare, "{bare}");
        assert!(!switch.is_default(), "{bare}");
        assert!(switch.explicit_kinds().is_empty(), "{bare}");
    }
    let listed = parse(Some("on:root,child"));
    assert_eq!(listed.source(), Source::Listed);
    assert_eq!(listed.explicit_kinds(), vec![Kind::Root, Kind::Child]);
    assert!(listed.explicit(Kind::Root) && !listed.explicit(Kind::IrHard));
    assert_eq!(parse(Some("check:c")).explicit_kinds(), vec![Kind::IrHard]);
    // an inert listed kind is never binding (it is not enabled)
    assert_eq!(
        parse(Some("on:root,quotient")).explicit_kinds(),
        vec![Kind::Root]
    );
    let all = parse(Some("on:root,child,ir-hard"));
    assert_eq!(all.enabled_kinds(), default.enabled_kinds());
    assert_eq!(all.explicit_kinds(), DEFAULT_KINDS.to_vec());
    assert_ne!(all, default, "the source differs");
    // the fields appended to UNCONTRACTED_SWITCH
    assert_eq!(default.source_fields(), " source=default binding=none");
    assert_eq!(
        parse(Some("on")).source_fields(),
        " source=bare binding=none"
    );
    assert_eq!(listed.source_fields(), " source=listed binding=root,child");
}

/// The unset switch follows [`UNSET_IS_ON`]; an EMPTY value is off (ruling
/// D4), as are `0` and `off`, with or without blanks.
#[test]
fn unset_follows_the_default_and_empty_is_off() {
    // the test build keeps the unset switch off unless the own-process tests
    // ask for the production meaning (`unset_is_on`)
    assert_eq!(
        parse(None),
        if unset_is_on() {
            Switch::DEFAULT_ON
        } else {
            Switch::OFF
        }
    );
    assert!(!unset_is_on() || UNSET_IS_ON);
    for off in ["", " ", "0", "off", " off "] {
        assert_eq!(parse(Some(off)), Switch::OFF, "{off:?}");
    }
    // a malformed value is refused as before
    for bad in ["of", "on:", "default", "off:root"] {
        assert!(Switch::parse(Some(bad)).is_err(), "{bad}");
    }
}

/// Ruling D2: a default run and an explicit
/// `on:root,child,ir-hard` run have ONE stored identity (they
/// resume each other), and a bare `on` the same.
#[test]
fn default_identity_equals_the_explicit_default_list() {
    let identity = |switch: Switch| {
        with_hole_cut(HoleCut::Form, || {
            ledger::identity_of(switch, &Exceptions::default())
        })
    };
    let default = identity(Switch::DEFAULT_ON);
    assert_eq!(
        default.as_deref(),
        Some("on:root,child,ir-hard:strict")
    );
    assert_eq!(identity(parse(Some("on:root,child,ir-hard"))), default);
    assert_eq!(identity(parse(Some("on:a,c,b"))), default);
    // a list with surrogate is another identity (the hole-cut segment)
    assert_eq!(
        identity(parse(Some("on:root,child,ir-hard,surrogate"))).as_deref(),
        Some("on:root,child,ir-hard,surrogate:strict:hole_cut=form")
    );
    assert_eq!(identity(parse(Some("on"))), default);
    // a default with a dropped kind is the explicit list without it
    assert_eq!(
        identity(Switch::DEFAULT_ON.without(Kind::IrHard)),
        identity(parse(Some("on:root,child")))
    );
    assert_eq!(identity(Switch::OFF), None);
}

/// A default kind whose route a setting turns off is DROPPED with a notice
/// (`GAMMALOOP_IR_HARD=off`; `GAMMALOOP_UV_CHILDREN=massive` turns the
/// default off instead); a written value that enables the kind is refused
/// exactly as before.
#[test]
fn default_kind_route_off_is_dropped() -> Result<()> {
    let initialised = init(Switch::DEFAULT_ON, routes(true, false), false).map_err(|e| eyre!(e))?;
    ensure!(
        initialised.switch.enabled_kinds() == vec![Kind::Root, Kind::Child]
            && initialised.switch.is_default()
            && initialised.dropped == vec![(Kind::IrHard, "route_off")]
            && initialised.demoted.is_none(),
        "{initialised:?}"
    );
    // GAMMALOOP_UV_CHILDREN=massive (with or without the IR hard route) turns
    // a default OFF as a whole (beta's ruling 2026-11-23): see
    // `default_demoted_by_massive_children`
    for ir_hard_child in [true, false] {
        let initialised =
            init(Switch::DEFAULT_ON, routes(false, ir_hard_child), false).map_err(|e| eyre!(e))?;
        ensure!(
            initialised.switch == Switch::OFF
                && initialised.demoted.as_deref() == Some("uv_children:massive")
                && initialised.dropped.is_empty(),
            "{initialised:?}"
        );
    }
    ensure!(
        ledger::dropped_text(Kind::IrHard, "route_off").contains("GAMMALOOP_IR_HARD"),
        "notice"
    );
    // the production routes drop nothing
    let initialised = init(Switch::DEFAULT_ON, Routes::PRODUCTION, false).map_err(|e| eyre!(e))?;
    ensure!(
        initialised.switch == Switch::DEFAULT_ON && initialised.dropped.is_empty(),
        "{initialised:?}"
    );
    // written values keep their refusals (bare and listed)
    for value in ["on", "check", "on:ir-hard", "on:root,child,ir-hard,surrogate"] {
        let error = init(parse(Some(value)), routes(true, false), false).unwrap_err();
        ensure!(error.contains("GAMMALOOP_IR_HARD"), "{value}: {error}");
    }
    for value in ["on", "check:child"] {
        let error = init(parse(Some(value)), routes(false, true), false).unwrap_err();
        ensure!(
            error.contains("GAMMALOOP_UV_CHILDREN=massive"),
            "{value}: {error}"
        );
    }
    // a written list without the kind is not affected
    ensure!(
        init(parse(Some("on:root,child,surrogate")), routes(true, false), false)
            .is_ok_and(|initialised| initialised.dropped.is_empty()),
        "on:root,child,surrogate"
    );
    Ok(())
}

/// `GAMMALOOP_UV_CHILDREN=massive` turns a DEFAULT switch OFF as a whole
/// (beta's ruling 2026-11-23: a default must never fail a configuration that
/// succeeds with the switch off; the massive children's Δ make the other
/// kinds' records error), with the reason `uv_children:massive`; a WRITTEN
/// value keeps today's behaviour: one that enables child is refused at
/// start-up, one without child (e.g. `on:root`) runs.
#[test]
fn default_demoted_by_massive_children() -> Result<()> {
    let massive = routes(false, true);
    let initialised = init(Switch::DEFAULT_ON, massive.clone(), false).map_err(|e| eyre!(e))?;
    ensure!(
        initialised.switch == Switch::OFF
            && initialised.demoted.as_deref() == Some("uv_children:massive"),
        "{initialised:?}"
    );
    for value in ["on", "check", "on:child", "on:root,child,ir-hard"] {
        let error = init(parse(Some(value)), massive.clone(), false).unwrap_err();
        ensure!(
            error.contains("GAMMALOOP_UV_CHILDREN=massive"),
            "{value}: {error}"
        );
    }
    for value in ["on:root", "check:root,ir-hard"] {
        let initialised = init(parse(Some(value)), massive.clone(), false).map_err(|e| eyre!(e))?;
        ensure!(
            initialised.switch == parse(Some(value)) && initialised.demoted.is_none(),
            "{value}: {initialised:?}"
        );
    }
    Ok(())
}

/// `GAMMALOOP_IR_DIRAC_CLOSURE=off` drops ir-hard from a default switch; a
/// written value with ir-hard is refused.
#[test]
fn default_ir_dirac_off_drops_ir_hard() -> Result<()> {
    let off = Routes {
        ir_dirac: Ok(super::ir_dirac::IrDiracMode::Off),
        ..Routes::PRODUCTION
    };
    let initialised = init(Switch::DEFAULT_ON, off.clone(), false).map_err(|e| eyre!(e))?;
    ensure!(
        initialised.switch.enabled_kinds() == vec![Kind::Root, Kind::Child]
            && initialised.dropped == vec![(Kind::IrHard, "ir_dirac_off")],
        "{initialised:?}"
    );
    for value in ["on", "check:ir-hard"] {
        let error = init(parse(Some(value)), off.clone(), false).unwrap_err();
        ensure!(
            error.contains("GAMMALOOP_IR_DIRAC_CLOSURE=off"),
            "{value}: {error}"
        );
    }
    Ok(())
}

/// A quotient route other than `rstar` and `PINT6_CLOSURE_DUMP` would refuse
/// every kind: a default switch is turned OFF with a notice; written values
/// are refused as before.
#[test]
fn default_demoted_by_ir_quotient_and_pint6() -> Result<()> {
    for (route, name) in [
        (IrQuotient::Massive, "massive"),
        (IrQuotient::Off, "off"),
        (IrQuotient::Check, "check"),
        (IrQuotient::Auto, "auto"),
    ] {
        let quotient = Routes {
            ir_quotient: Ok(route),
            ..Routes::PRODUCTION
        };
        let initialised = init(Switch::DEFAULT_ON, quotient.clone(), false).map_err(|e| eyre!(e))?;
        ensure!(
            initialised.switch == Switch::OFF
                && initialised.demoted.as_deref() == Some(&format!("ir_quotient:{name}")[..]),
            "{name}: {initialised:?}"
        );
        for value in ["on", "on:root"] {
            let error = init(parse(Some(value)), quotient.clone(), false).unwrap_err();
            ensure!(
                error.contains(&format!("GAMMALOOP_IR_QUOTIENT={name} is refused")),
                "{value}: {error}"
            );
        }
    }
    // a malformed quotient value stays an error (it is refused by its own
    // reader too)
    let malformed = Routes {
        ir_quotient: Err("GAMMALOOP_IR_QUOTIENT must be one of ...".into()),
        ..Routes::PRODUCTION
    };
    ensure!(init(Switch::DEFAULT_ON, malformed, false).is_err(), "malformed");
    let initialised = init(Switch::DEFAULT_ON, Routes::PRODUCTION, true).map_err(|e| eyre!(e))?;
    ensure!(
        initialised.switch == Switch::OFF
            && initialised.demoted.as_deref() == Some("pint6_dump"),
        "{initialised:?}"
    );
    for value in ["on", "check:root"] {
        let error = init(parse(Some(value)), Routes::PRODUCTION, true).unwrap_err();
        ensure!(error.contains("PINT6_CLOSURE_DUMP"), "{value}: {error}");
    }
    // the other start-up refusals stay refusals for a default too: they are
    // the switch's own variables (STRICT=0, a bad exception list)
    ensure!(
        ledger::initialise(
            Switch::DEFAULT_ON,
            Some("0"),
            None,
            || Routes::PRODUCTION,
            || false,
            None
        )
        .is_err(),
        "STRICT=0"
    );
    ensure!(
        ledger::initialise(
            Switch::DEFAULT_ON,
            None,
            None,
            || Routes::PRODUCTION,
            || false,
            Some("/nonexistent")
        )
        .is_err(),
        "exceptions"
    );
    Ok(())
}

/// The per-diagram "nonzero counterterm with no record" rule keeps its
/// trigger (root, child and ir-hard enabled): the U7b list
/// `root,child,ir-hard` (the default kinds), lists with surrogate and the
/// bare `on` / `check` keep it as a failure; the DEFAULT switch has it as a
/// notice.
#[test]
fn the_diagram_rule_keeps_its_trigger() {
    let failures = |switch: Switch| {
        with_switch(switch, || {
            with_hole_cut(HoleCut::Form, || {
                let evidence = serde_json::json!({
                    "schema": ledger::EVIDENCE_SCHEMA,
                    "identity": ledger::current_identity(),
                    "kinds": {},
                });
                let mut tally = ledger::RunTally::new().expect("on");
                tally.add("g/000001", Some(&evidence), true);
                tally.failures()
            })
        })
    };
    let rule = |failures: &[String]| {
        failures
            .iter()
            .any(|failure| failure.contains("the counterterm is nonzero"))
    };
    for switch in [
        parse(Some("check:root,child,ir-hard")),
        parse(Some("on:root,child,ir-hard,surrogate")),
        parse(Some("on")),
        parse(Some("check")),
    ] {
        let failures = failures(switch);
        assert!(rule(&failures), "{switch:?}: {failures:?}");
    }
    // under the DEFAULT switch the rule is a notice (`the_diagram_rule_is_a_notice_under_the_default`)
    assert!(!rule(&failures(Switch::DEFAULT_ON)));
    // a partial list: the kind-level rule covers it
    let failures = failures(parse(Some("on:root,child")));
    assert!(!rule(&failures), "{failures:?}");
}

/// Beta's ruling of 2026-11-23 (U8 N): under a DEFAULT switch (the variable
/// unset) the run-level diagram rule is a NOTICE, never a failure: a
/// constructed diagram that meets the rule (nonzero counterterm, no record of
/// any kind taken, none all-declined or excepted) leaves the run's failures
/// empty and gives one notice naming the diagram, the rule and "default kind,
/// not binding"; the same evidence under a WRITTEN value (bare `on`, bare
/// `check`, a list of the same kinds) fails as before, byte for byte; a zero
/// counterterm meets no rule.
#[test]
fn the_diagram_rule_is_a_notice_under_the_default() {
    let tally_of = |switch: Switch, nonzero: bool| {
        with_switch(switch, || {
            with_hole_cut(HoleCut::Form, || {
                // every enabled kind expected its records, none was taken or
                // declined (the hooks missed them): no evidence, no excepted kind
                let mut kinds = serde_json::Map::new();
                for kind in switch.enabled_kinds() {
                    kinds.insert(
                        kind.as_str().to_owned(),
                        serde_json::json!({"expected": 1, "attempted": 0, "outcome": "ok"}),
                    );
                }
                let evidence = serde_json::json!({
                    "schema": ledger::EVIDENCE_SCHEMA,
                    "identity": ledger::current_identity(),
                    "kinds": kinds,
                });
                let mut tally = ledger::RunTally::new().expect("on");
                tally.add("g/000007", Some(&evidence), nonzero);
                (tally.failures(), tally.notices())
            })
        })
    };
    let rule = "g/000007: the counterterm is nonzero, and no record of any kind was";
    let (failures, notices) = tally_of(Switch::DEFAULT_ON, true);
    assert!(
        !failures.iter().any(|failure| failure.starts_with(rule)),
        "{failures:#?}"
    );
    assert_eq!(
        notices
            .iter()
            .filter(|notice| notice.starts_with(rule))
            .cloned()
            .collect::<Vec<_>>(),
        vec![
            "g/000007: the counterterm is nonzero, and no record of any kind was taken; default kind, not binding: the contracted value is used"
                .to_owned()
        ],
        "{notices:#?}"
    );
    for (value, verb) in [
        ("on", "taken"),
        ("check", "compared"),
        ("on:root,child,ir-hard", "taken"),
        ("check:root,child,ir-hard", "compared"),
    ] {
        let switch = parse(Some(value));
        assert_eq!(
            switch.enabled_kinds(),
            Switch::DEFAULT_ON.enabled_kinds(),
            "{value}"
        );
        let (failures, notices) = tally_of(switch, true);
        assert!(
            failures.contains(&format!(
                "g/000007: the counterterm is nonzero, and no record of any kind was {verb}"
            )),
            "{value}: {failures:#?}"
        );
        assert!(
            !notices.iter().any(|notice| notice.starts_with(rule)),
            "{value}: {notices:#?}"
        );
    }
    for switch in [Switch::DEFAULT_ON, parse(Some("on"))] {
        let (failures, notices) = tally_of(switch, false);
        assert!(
            !failures
                .iter()
                .chain(&notices)
                .any(|line| line.starts_with(rule)),
            "zero counterterm: {failures:#?} {notices:#?}"
        );
    }
}

/// Own process: the environment's switch with the variable unset and
/// `GAMMALOOP_IR_HARD=off` (the six-loop setting): with the default on,
/// root and child and no refusal; with it off, off.
#[test]
fn unset_switch_with_ir_hard_off_in_the_environment() -> Result<()> {
    ensure!(
        std::env::var_os("NEXTEST").is_some(),
        "this test sets process environment variables: run it under cargo nextest (one process per test)"
    );
    unsafe {
        std::env::remove_var(uncontracted::VARIABLE);
        std::env::set_var("GAMMALOOP_IR_HARD", "off");
        // the production meaning of an unset switch (`unset_is_on`)
        std::env::set_var(TEST_DEFAULT_VARIABLE, "1");
    }
    super::validate_uncontracted()?;
    let switch = uncontracted::switch();
    if UNSET_IS_ON {
        ensure!(
            switch.enabled_kinds() == vec![Kind::Root, Kind::Child]
                && switch.is_default(),
            "{switch:?}"
        );
        ensure!(
            super::uncontracted_source()?.as_deref() == Some("default"),
            "source"
        );
    } else {
        ensure!(switch == Switch::OFF, "{switch:?}");
    }
    Ok(())
}

/// Own process: an explicit list with `GAMMALOOP_IR_HARD=off` is still
/// refused at start-up.
#[test]
fn explicit_ir_hard_with_ir_hard_off_is_refused() -> Result<()> {
    ensure!(
        std::env::var_os("NEXTEST").is_some(),
        "this test sets process environment variables: run it under cargo nextest (one process per test)"
    );
    unsafe {
        std::env::set_var(uncontracted::VARIABLE, "on:ir-hard");
        std::env::set_var("GAMMALOOP_IR_HARD", "off");
    }
    let error = super::validate_uncontracted().expect_err("refused");
    ensure!(format!("{error}").contains("GAMMALOOP_IR_HARD"), "{error}");
    Ok(())
}

/// Own process: `uncontracted_demote_default` turns a default switch off
/// (the beta's `--grouping skeletons`) and leaves a written value alone.
#[test]
fn demote_default_turns_only_a_default_off() -> Result<()> {
    ensure!(
        std::env::var_os("NEXTEST").is_some(),
        "this test sets process environment variables: run it under cargo nextest (one process per test)"
    );
    unsafe {
        std::env::remove_var(uncontracted::VARIABLE);
        std::env::set_var(TEST_DEFAULT_VARIABLE, "1");
    }
    let demoted = super::uncontracted_demote_default("grouping_skeletons")?;
    ensure!(demoted == UNSET_IS_ON, "unset: {demoted}");
    super::validate_uncontracted()?;
    ensure!(uncontracted::switch() == Switch::OFF, "off after demotion");
    ensure!(super::uncontracted_identity()?.is_none(), "identity");
    Ok(())
}

/// Own process: a written value is never demoted.
#[test]
fn demote_default_leaves_a_written_value() -> Result<()> {
    ensure!(
        std::env::var_os("NEXTEST").is_some(),
        "this test sets process environment variables: run it under cargo nextest (one process per test)"
    );
    unsafe { std::env::set_var(uncontracted::VARIABLE, "on:root") };
    ensure!(
        !super::uncontracted_demote_default("grouping_skeletons")?,
        "a written value was demoted"
    );
    super::validate_uncontracted()?;
    ensure!(
        uncontracted::switch().enabled_kinds() == vec![Kind::Root],
        "on:root"
    );
    unsafe { std::env::set_var(uncontracted::VARIABLE, "bogus") };
    ensure!(
        super::uncontracted_demote_default("grouping_skeletons").is_err(),
        "a malformed value is the start-up error"
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// G-4: the run-level rule binds only the kinds named in a list (ruling D1)
// ---------------------------------------------------------------------------

/// A run of one diagram under `switch` that takes (on) / compares (check)
/// one record of each kind in `kinds`: its tally's failures, notices and
/// default kinds without a record.
fn tally_of(switch: Switch, kinds: &[Kind]) -> Result<(Vec<String>, Vec<String>, Vec<&'static str>)> {
    with_switch(switch, || {
        with_hole_cut(HoleCut::Form, || {
            let mut tally = ledger::RunTally::new().ok_or_else(|| eyre!("no tally"))?;
            let opened = ledger::DiagramLedger::open("u8-g4", "000001")
                .ok_or_else(|| eyre!("no ledger"))?;
            for kind in kinds {
                ledger::expect(*kind);
                ledger::record_value(*kind, switch.mode);
            }
            let evidence = opened.close(true)?;
            tally.add("u8-g4/000001", evidence.as_ref(), true);
            Ok((tally.failures(), tally.notices(), tally.default_no_record()))
        })
    })
}

const U7B: [Kind; 3] = [Kind::Root, Kind::Child, Kind::IrHard];

/// A default (or bare) kind with no record over the run is a notice, not a
/// failure.
#[test]
fn run_tally_default_kind_without_record_is_a_notice() -> Result<()> {
    for switch in [Switch::DEFAULT_ON, parse(Some("on")), parse(Some("check"))] {
        let (failures, notices, no_record) = tally_of(switch, &[Kind::Root, Kind::Child])?;
        ensure!(failures.is_empty(), "{switch:?}: {failures:#?}");
        ensure!(no_record == vec!["ir-hard"], "{switch:?}: {no_record:?}");
        let verb = if switch.mode == Mode::Check {
            "compared"
        } else {
            "taken"
        };
        ensure!(
            notices
                == vec![format!(
                    "kind ir-hard: no record {verb} over the run (0 expected, 0 declined); default kind, not binding"
                )],
            "{switch:?}: {notices:#?}"
        );
    }
    // the measured L1 case: root records only
    let (failures, _, no_record) = tally_of(Switch::DEFAULT_ON, &[Kind::Root])?;
    ensure!(
        failures.is_empty() && no_record == vec!["child", "ir-hard"],
        "{failures:#?} {no_record:?}"
    );
    Ok(())
}

/// A listed kind with no record fails the run with the existing message,
/// byte for byte, and gives no notice.
#[test]
fn run_tally_listed_kind_without_record_fails() -> Result<()> {
    let (failures, notices, no_record) =
        tally_of(parse(Some("on:root,child,ir-hard,surrogate")), &U7B)?;
    ensure!(
        failures
            == vec![
                "kind surrogate: no record taken over the run (0 expected, 0 declined); a run-level exception (*<TAB>surrogate<TAB>none_compared<TAB>reason in GAMMALOOP_UNCONTRACTED_EXCEPTIONS) is required"
                    .to_owned()
            ],
        "{failures:#?}"
    );
    ensure!(notices.is_empty() && no_record.is_empty(), "{notices:#?}");
    let (failures, ..) = tally_of(parse(Some("check:root,child,ir-hard,surrogate")), &U7B)?;
    ensure!(
        failures.len() == 1 && failures[0].starts_with("kind surrogate: no record compared"),
        "{failures:#?}"
    );
    Ok(())
}

/// A list binds its own kinds only: `on:root,surrogate` with root records
/// only fails on surrogate and gives no notice.
#[test]
fn run_tally_mixed() -> Result<()> {
    let (failures, notices, no_record) = tally_of(parse(Some("on:root,surrogate")), &[Kind::Root])?;
    ensure!(
        failures.len() == 1 && failures[0].starts_with("kind surrogate: no record taken"),
        "{failures:#?}"
    );
    ensure!(notices.is_empty() && no_record.is_empty(), "{notices:#?}");
    // every kind with a record: nothing either way
    let (failures, notices, _) = tally_of(Switch::DEFAULT_ON, &DEFAULT_KINDS)?;
    ensure!(failures.is_empty() && notices.is_empty(), "{failures:#?} {notices:#?}");
    Ok(())
}

/// A run-level exception covers a default kind too: no notice.
#[test]
fn run_tally_exception_covers_a_default_kind() -> Result<()> {
    let exceptions = ledger::parse_exceptions("*\tir-hard\tnone_compared\tno IR correction in this set\n")
        .map_err(|error| eyre!(error))?;
    let (failures, notices, no_record) = ledger::with_exceptions(exceptions.clone(), || {
        tally_of(Switch::DEFAULT_ON, &[Kind::Root, Kind::Child])
    })?;
    ensure!(
        failures.is_empty() && notices.is_empty() && no_record.is_empty(),
        "{failures:#?} {notices:#?}"
    );
    let (failures, ..) = ledger::with_exceptions(exceptions, || {
        tally_of(parse(Some("on:root,child,ir-hard")), &[Kind::Root, Kind::Child])
    })?;
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

// ---------------------------------------------------------------------------
// G-6: the all-declined coverage rules bind only the kinds named in a list
// (ruling 18:11Z, in the sense of D1)
// ---------------------------------------------------------------------------

/// One diagram under `switch` whose one root record declined by design
/// (`code`, frame K1): the ledger's result (Err: the diagram failed) and its
/// lines.
fn one_declined_root(switch: Switch, code: &str) -> (Result<Option<serde_json::Value>>, Vec<String>) {
    with_switch(switch, || {
        with_hole_cut(HoleCut::Form, || {
            let _ = ledger::take_lines();
            let opened = ledger::DiagramLedger::open("u8-g6", "000001").expect("a ledger");
            ledger::record_shaped(Kind::Root, "K1", Some(code), switch.mode);
            let closed = opened.close(true);
            (closed, ledger::take_lines())
        })
    })
}

/// A default / bare kind whose every record is a designed decline (the
/// one-record uvind_head case of helper VH 2's L1, on a root record): the
/// diagram passes, outcome `declined`, the notice; an explicit list fails it
/// `uncovered` as before; an exception naming the code makes it `excepted`.
#[test]
fn all_declined_binds_listed_kinds_only() -> Result<()> {
    for switch in [Switch::DEFAULT_ON, parse(Some("on")), parse(Some("check"))] {
        let (closed, lines) = one_declined_root(switch, "uvind_head");
        let evidence = closed.map_err(|error| eyre!("{switch:?}: {error}"))?.ok_or_else(|| eyre!("no evidence"))?;
        ensure!(
            evidence["outcome"] == "ok" && evidence["kinds"]["root"]["outcome"] == "declined",
            "{switch:?}: {evidence}"
        );
        ensure!(
            lines.iter().any(|line| line.starts_with("UNCONTRACTED_NOTICE ")
                && line.contains("kind root: all 1 records declined by design (uvind_head:1); contracted value used")),
            "{switch:?}: {lines:#?}"
        );
        // the run tally: covered (no diagram-rule failure)
        let failures = with_switch(switch, || {
            with_hole_cut(HoleCut::Form, || {
                let mut tally = ledger::RunTally::new().expect("on");
                tally.add("u8-g6/000001", Some(&evidence), true);
                tally.failures()
            })
        });
        ensure!(
            !failures.iter().any(|failure| failure.contains("the counterterm is nonzero")),
            "{switch:?}: {failures:#?}"
        );
    }
    for value in ["on:root,child,ir-hard", "check:root"] {
        let (closed, _) = one_declined_root(parse(Some(value)), "uvind_head");
        let error = closed.expect_err(value).to_string();
        ensure!(
            error.contains("outcome=uncovered") && error.contains("all 1 records declined"),
            "{value}: {error}"
        );
    }
    let exceptions = ledger::parse_exceptions("u8-g6\troot\tuvind_head\ta constructed uvind record\n")
        .map_err(|error| eyre!(error))?;
    let (closed, _) = ledger::with_exceptions(exceptions, || {
        one_declined_root(parse(Some("on:root,child,ir-hard")), "uvind_head")
    });
    let evidence = closed?.ok_or_else(|| eyre!("no evidence"))?;
    ensure!(evidence["kinds"]["root"]["outcome"] == "excepted", "{evidence}");
    Ok(())
}

/// The setting-wide class rule: a default kind's records of one class all
/// declined with a setting's code are a notice, not a failure (diagram and
/// run); a listed kind fails as before.
#[test]
fn setting_wide_classes_bind_listed_kinds_only() -> Result<()> {
    let code = uncontracted::SETTING_CODES[0];
    let (closed, lines) = one_declined_root(Switch::DEFAULT_ON, code);
    let evidence = closed?.ok_or_else(|| eyre!("no evidence"))?;
    ensure!(evidence["outcome"] == "ok", "{evidence}");
    ensure!(lines.iter().any(|line| line.starts_with("UNCONTRACTED_NOTICE ")), "{lines:#?}");
    let (failures, notices) = with_switch(Switch::DEFAULT_ON, || {
        with_hole_cut(HoleCut::Form, || {
            let mut tally = ledger::RunTally::new().expect("on");
            tally.add("u8-g6/000001", Some(&evidence), true);
            (tally.failures(), tally.notices())
        })
    });
    ensure!(
        !failures.iter().any(|failure| failure.contains("-shaped records over the run declined")),
        "{failures:#?}"
    );
    ensure!(
        notices.iter().any(|notice| notice.contains("K1-shaped records over the run declined with a setting's code")),
        "{notices:#?}"
    );
    let (closed, _) = one_declined_root(parse(Some("on:root")), code);
    ensure!(closed.is_err(), "a listed kind passed");
    Ok(())
}
