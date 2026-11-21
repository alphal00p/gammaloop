//! A value WITHOUT a new-feature marker (kind closure, insertions, a v3
//! direct record: `FactorizedNumerator::has_new_features`) is validated and
//! evaluated by exactly the base's rules. These tests pin the base's values
//! and refusals where the series once differed, next to the differential
//! harness (direct_differential_tests.rs):
//! - the v1 bubble with 0, 6 and 7 external spectators: the base's integrand
//!   and value, as printed by the base build;
//! - a v1 text keeps the base reader: `insertions n` other than 0 is
//!   `Reserved`; a marked value is written as v3, an unmarked one as v1;
//! - a declared section count never allocates: huge counts of the v3
//!   sections are refused as text errors (the process does not abort);
//! - a port that occurs only in the projector is refused for marked values
//!   too.
use std::path::PathBuf;

use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;
use vakint::p_integral::factorized::{
    DirectRecord, FACTORIZED_SCHEMA_V3, Insertion, evaluate_factorized, vakint_input,
};
use vakint::p_integral::{FACTORIZED_SCHEMA, FactorizedError, FactorizedNumerator};
use vakint::refpath;

fn fixture(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn rank(n: usize) -> String {
    fixture(&format!(
        "tests/fixtures/direct_differential/rank{n}-v1.txt"
    ))
}

fn direct_v3() -> String {
    fixture("src/p_integral/factorized/fixtures/wh-k5-direct.record.txt")
}

fn closure_v3() -> String {
    fixture("src/p_integral/factorized/fixtures/wh-k1-root.record.txt")
}

fn marked(value: &FactorizedNumerator) -> FactorizedNumerator {
    FactorizedNumerator {
        record: Some(DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
        ..value.clone()
    }
}

/// The base build's output for the reviewer's three records (integrand of
/// `vakint_input`, value of `evaluate_factorized`, settings
/// `refpath::vakint_settings(1)`, prefactor `prefactor.total()`), printed by
/// the base dea282615.
const BASE: &[(usize, &str, &str)] = &[
    (
        0,
        "topo(prop(2,edge(0,1),k(1),0,1)*prop(3,edge(0,1),-p(1)-k(1),0,1))",
        "1/ε",
    ),
    (
        6,
        "p(1,hedge(0,1))*p(1,hedge(0,2))*p(1,hedge(0,3))*p(1,hedge(0,4))*p(1,hedge(0,5))*p(1,hedge(0,6))*topo(prop(2,edge(0,1),k(1),0,1)*prop(3,edge(0,1),-p(1)-k(1),0,1))",
        "p(1,hedge(0,1))*p(1,hedge(0,2))*p(1,hedge(0,3))*p(1,hedge(0,4))*p(1,hedge(0,5))*p(1,hedge(0,6))/ε",
    ),
    (
        7,
        "p(1,hedge(0,1))*p(1,hedge(0,2))*p(1,hedge(0,3))*p(1,hedge(0,4))*p(1,hedge(0,5))*p(1,hedge(0,6))*p(1,hedge(0,7))*topo(prop(2,edge(0,1),k(1),0,1)*prop(3,edge(0,1),-p(1)-k(1),0,1))",
        "p(1,hedge(0,1))*p(1,hedge(0,2))*p(1,hedge(0,3))*p(1,hedge(0,4))*p(1,hedge(0,5))*p(1,hedge(0,6))*p(1,hedge(0,7))/ε",
    ),
];

#[test]
fn the_spectator_bubbles_keep_the_base_values() {
    Vakint::initialize_vakint_symbols();
    let vakint = Vakint::new().unwrap();
    let settings = refpath::vakint_settings(1);
    for (n, integrand, result) in BASE {
        let value = FactorizedNumerator::from_text(&rank(*n)).unwrap();
        assert!(
            !value.has_new_features(),
            "rank {n}: a v1 record is unmarked"
        );
        assert_eq!(value.ports.len(), *n);
        let prefactor = value.prefactor.total();
        let (input, _) = vakint_input(&settings, &value, &prefactor)
            .unwrap_or_else(|e| panic!("rank {n}: the base accepted it: {e}"));
        assert_eq!(input.to_string(), *integrand, "rank {n}: integrand");
        let got = evaluate_factorized(&vakint, &settings, &value, &prefactor)
            .unwrap_or_else(|e| panic!("rank {n}: the base evaluated it: {e}"));
        assert_eq!(got.value.to_string(), *result, "rank {n}: value");
    }
}

#[test]
fn a_v1_text_keeps_the_base_reader() {
    let text = rank(0);
    assert!(FactorizedNumerator::parse_text(&text).is_ok());
    for count in ["1", "x", "18446744073709551615", "100000000000"] {
        let edited = text.replace("insertions\t0", &format!("insertions\t{count}"));
        assert_ne!(edited, text);
        match FactorizedNumerator::parse_text(&edited) {
            Err(FactorizedError::Reserved { feature }) => {
                assert!(feature.starts_with("insertion vertices"), "{feature}")
            }
            other => panic!("v1 `insertions {count}`: {other:?}"),
        }
    }
    // the same count in a v3 record is a section (here missing)
    let v3 = marked(&FactorizedNumerator::from_text(&text).unwrap())
        .to_text()
        .unwrap();
    assert!(v3.starts_with(&format!("schema\t{FACTORIZED_SCHEMA_V3}\n")));
    assert!(FactorizedNumerator::parse_text(&v3).is_ok());
    let v3 = v3.replace("insertions\t0", "insertions\t1");
    assert!(matches!(
        FactorizedNumerator::parse_text(&v3),
        Err(FactorizedError::TextSectionCount { .. })
    ));
}

#[test]
fn a_marked_value_is_written_as_v3_and_an_unmarked_one_as_v1() {
    let schema = |text: &str| text.lines().next().unwrap().to_string();
    let v1 = FactorizedNumerator::from_text(&rank(6)).unwrap();
    let written = v1.to_text().unwrap();
    assert_eq!(schema(&written), format!("schema\t{FACTORIZED_SCHEMA}"));
    assert_eq!(FactorizedNumerator::from_text(&written).unwrap(), v1);
    // marked by a record
    let written = marked(&v1).to_text().unwrap();
    assert_eq!(schema(&written), format!("schema\t{FACTORIZED_SCHEMA_V3}"));
    // marked by an insertion only (no record): v3 as well, never v1
    let mut with_insertion = FactorizedNumerator::parse_text(&rank(0)).unwrap();
    with_insertion.insertions.push(Insertion {
        vertex: 1,
        hedges: vec![0, 3],
        cycles: (0, 0),
        coefficient: Atom::num(1),
        factor: Atom::num(1),
    });
    assert!(with_insertion.has_new_features());
    let written = with_insertion.to_text().unwrap();
    assert_eq!(schema(&written), format!("schema\t{FACTORIZED_SCHEMA_V3}"));
    assert!(written.contains("\ninsertions\t1\n"), "{written}");
    let back = FactorizedNumerator::parse_text(&written).unwrap();
    assert_eq!(back.insertions, with_insertion.insertions);
    // the real v3 records
    for text in [direct_v3(), closure_v3()] {
        let value = FactorizedNumerator::from_text(&text).unwrap();
        assert!(value.has_new_features());
        assert_eq!(
            schema(&value.to_text().unwrap()),
            format!("schema\t{FACTORIZED_SCHEMA_V3}")
        );
    }
}

/// `text` with the line starting `key\t` replaced by `key\t<count>` (or a
/// `key\t<count>` line inserted after the line starting `after`, and the
/// `lines` count raised by one: the refusal must come from the section).
fn with_count(text: &str, key: &str, after: Option<&str>, count: &str) -> String {
    let mut out = Vec::new();
    let mut done = false;
    for line in text.lines() {
        match after {
            None if line.starts_with(&format!("{key}\t")) && !done => {
                out.push(format!("{key}\t{count}"));
                done = true;
            }
            Some(a) if line.starts_with(a) && !done => {
                out.push(line.to_string());
                out.push(format!("{key}\t{count}"));
                done = true;
            }
            // one line more: the record's line count follows
            Some(_) if line.starts_with("lines\t") => {
                let n: usize = line["lines\t".len()..].parse().unwrap();
                out.push(format!("lines\t{}", n + 1));
            }
            _ => out.push(line.to_string()),
        }
    }
    assert!(done, "{key}");
    out.join("\n") + "\n"
}

/// A refusal by the section reader itself (not the record's line count).
fn is_text_error(error: &FactorizedError) -> bool {
    matches!(
        error,
        FactorizedError::TextSyntax { .. } | FactorizedError::TextSectionCount { .. }
    )
}

/// A declared count of a new section is never an allocation size (A-11):
/// each huge count is a text refusal; the control with the true count reads.
#[test]
fn huge_counts_of_the_new_sections_are_refused_not_allocated() {
    let direct = direct_v3();
    let closure = closure_v3();
    assert!(FactorizedNumerator::parse_text(&direct).is_ok());
    assert!(FactorizedNumerator::parse_text(&closure).is_ok());
    // (record, key, inserted after, true count)
    let sections: [(&str, &str, &str, Option<&str>, &str); 6] = [
        ("direct", &direct, "insertions", None, "1"),
        ("direct", &direct, "open", None, "2"),
        ("direct", &direct, "weights", Some("original_cycles\t"), "0"),
        ("closure", &closure, "softs", None, "0"),
        (
            "closure",
            &closure,
            "weights",
            Some("original_cycles\t"),
            "0",
        ),
        ("closure", &closure, "open", None, "0"),
    ];
    for (name, text, key, after, truth) in sections {
        let control = with_count(text, key, after, truth);
        if after.is_none() {
            assert_eq!(control, *text, "{name} {key}: the control is the record");
        }
        assert!(
            FactorizedNumerator::parse_text(&control).is_ok(),
            "{name} {key}: the true count reads"
        );
        for count in ["18446744073709551615", "100000000000"] {
            let huge = with_count(text, key, after, count);
            match FactorizedNumerator::parse_text(&huge) {
                Err(error) if is_text_error(&error) => {}
                other => panic!("{name} `{key} {count}`: {other:?}"),
            }
        }
    }
    // kind direct refuses soft legs before reading them
    for count in ["18446744073709551615", "100000000000"] {
        let huge = with_count(&direct, "softs", None, count);
        match FactorizedNumerator::parse_text(&huge) {
            Err(FactorizedError::TextSyntax { detail, .. }) => {
                assert!(detail.starts_with("kind direct with soft legs"), "{detail}")
            }
            other => panic!("direct `softs {count}`: {other:?}"),
        }
    }
}

/// A port that occurs only in the projector, never in a factor (A-1): the
/// base refuses it as `Ports`; a marked value too. Without the absent port
/// (the defect removed) both validate.
#[test]
fn a_port_only_in_the_projector_is_refused() {
    let text = rank(1);
    let value = FactorizedNumerator::parse_text(&text).unwrap();
    let port = value.ports[0].clone();
    let absent = symbolica::try_parse!(
        "spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(99,0))"
    )
    .unwrap();
    let q = |slot: &Atom| {
        symbolica::try_parse!(&format!(
            "gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}}::Q(0,{})",
            slot.to_canonical_string()
        ))
        .unwrap()
    };
    let mut closed = value.clone();
    closed.projector = Some(vakint::p_integral::factorized::Projector {
        tensor: q(&port),
        p2_power: -1,
    });
    closed.validate().unwrap();
    marked(&closed).validate().unwrap();
    let mut defect = closed.clone();
    defect.ports.push(absent.clone());
    defect.ports.sort_by_key(|p| p.to_canonical_string());
    defect.projector.as_mut().unwrap().tensor = q(&port) * q(&absent);
    for (what, v) in [("unmarked", defect.clone()), ("marked", marked(&defect))] {
        match v.validate() {
            Err(FactorizedError::Ports { detail }) => {
                assert!(detail.starts_with("port "), "{what}: {detail}");
                assert!(
                    detail.contains("does not occur in the factors"),
                    "{what}: {detail}"
                );
            }
            other => panic!("{what}: {other:?}"),
        }
    }
}

/// PRE-EXISTING (the base aborted here too): a declared count of the base's
/// own sections (loops, externals, edges, vertices, ports) is no allocation
/// size either; each huge count of a v1 text is a section refusal.
#[test]
fn huge_counts_of_the_base_sections_are_refused_not_allocated() {
    let text = rank(1);
    assert!(FactorizedNumerator::parse_text(&text).is_ok());
    let edits: [(&str, &str, Box<dyn Fn(&str) -> String>); 5] = [
        (
            "loops",
            "momenta\tloops=1\texternals=1",
            Box::new(|n| format!("momenta\tloops={n}\texternals=1")),
        ),
        (
            "externals",
            "momenta\tloops=1\texternals=1",
            Box::new(|n| format!("momenta\tloops=1\texternals={n}")),
        ),
        ("edges", "edges\t4", Box::new(|n| format!("edges\t{n}"))),
        (
            "vertices",
            "vertices\t2",
            Box::new(|n| format!("vertices\t{n}")),
        ),
        ("ports", "ports\t1", Box::new(|n| format!("ports\t{n}"))),
    ];
    for (name, line, edit) in &edits {
        assert_eq!(text.matches(&format!("\n{line}\n")).count(), 1, "{name}");
        for count in ["18446744073709551615", "100000000000"] {
            let huge = text.replace(&format!("\n{line}\n"), &format!("\n{}\n", edit(count)));
            match FactorizedNumerator::parse_text(&huge) {
                Err(error) if is_text_error(&error) => {}
                other => panic!("v1 {name} {count}: {other:?}"),
            }
        }
    }
}
