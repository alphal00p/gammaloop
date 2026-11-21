//! Every designed refusal of the ports entry, of the port rules of a marked
//! value and of the v3 direct-record text: each test builds a MINIMAL input
//! (a record that is accepted, with exactly one defect added), asserts the
//! refusal's VARIANT and DETAIL PREFIX, and asserts that the same input
//! without the defect is ACCEPTED (so a deleted rule makes the test fail).
//! Inputs: the vendored v1 bubble with `n` external spectators
//! (tests/fixtures/direct_differential/rank<n>-v1.txt, marked with an empty
//! v3 direct record), the beta's v3 kind-direct records k5-direct and K6
//! quark.
use std::path::PathBuf;

use symbolica::atom::Atom;
use vakint::Vakint;
use vakint::p_integral::factorized::{DirectRecord, evaluate_factorized_ports};
use vakint::p_integral::{FactorizedError, FactorizedNumerator};
use vakint::refpath;

const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const BIS: &str = "spenso::{spenso::representation,spenso::self_dual}::bis";
const DIM: &str = "gammalooprs::{}::dim";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const Q: &str = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
const G: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
const GAMMA: &str = "spenso::{linear,spenso::tensor}::gamma";

fn read(name: &str) -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn k5() -> String {
    read("src/p_integral/factorized/fixtures/wh-k5-direct.record.txt")
}

fn engine() -> Vakint {
    Vakint::initialize_vakint_symbols();
    Vakint::new().unwrap()
}

fn slot(h: usize, k: usize) -> String {
    format!("{MINK}({DIM},{HEDGE}({h},{k}))")
}

fn bis(h: usize, k: usize) -> String {
    format!("{BIS}(4,{HEDGE}({h},{k}))")
}

/// `text` with the unique occurrence of `from` replaced by `to`.
fn edit(text: &str, from: &str, to: &str) -> String {
    assert_eq!(text.matches(from).count(), 1, "not unique: {from:?}");
    text.replace(from, to)
}

/// `text` with its `lines` count raised by `n` (for added lines).
fn more_lines(text: &str, n: usize) -> String {
    text.lines()
        .map(|l| match l.strip_prefix("lines\t") {
            Some(c) => format!("lines\t{}", c.parse::<usize>().unwrap() + n),
            None => l.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

/// The one-loop bubble with external leg 0 and vertex 1 carrying `factor`,
/// the given ports and projector, MARKED (an empty v3 direct record).
fn bubble(factor: &str, ports: &[String], projector: Option<&str>) -> FactorizedNumerator {
    let text = read("tests/fixtures/direct_differential/rank0-v1.txt");
    let mut lines = Vec::new();
    for line in text.lines() {
        if line == "vertex\t1\trule=V_BQQ\t1" {
            lines.push(format!("vertex\t1\trule=V_BQQ\t{factor}"));
        } else if line == "ports\t0" {
            lines.push(format!("ports\t{}", ports.len()));
            lines.extend(ports.iter().map(|p| format!("port\t{p}")));
        } else if line == "projector\tnone" {
            lines.push(match projector {
                None => line.to_string(),
                Some(t) => format!("projector\tp2_power=0\t{t}"),
            });
        } else if line.starts_with("lines\t") {
            lines.push(format!("lines\t{}", 28 + ports.len()));
        } else {
            lines.push(line.to_string());
        }
    }
    let text = lines.join("\n") + "\n";
    let value = FactorizedNumerator::parse_text(&text).unwrap_or_else(|e| panic!("{e}\n{text}"));
    FactorizedNumerator {
        record: Some(DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
        ..value
    }
}

fn ports_entry(vakint: &Vakint, value: &FactorizedNumerator) -> Result<(), FactorizedError> {
    let settings = refpath::vakint_settings(value.loops() as i64 + 1);
    evaluate_factorized_ports(vakint, &settings, value, &Atom::num(1)).map(|_| ())
}

/// `Ports { detail }` with `detail` starting with `prefix`.
fn assert_ports(result: Result<(), FactorizedError>, prefix: &str, contains: &str) {
    match result {
        Err(FactorizedError::Ports { detail }) => {
            assert!(
                detail.starts_with(prefix) && detail.contains(contains),
                "{detail}"
            )
        }
        other => panic!("expected Ports({prefix}..), got {other:?}"),
    }
}

/// `TextSyntax { detail }` with `detail` starting with `prefix`.
fn assert_text(result: Result<FactorizedNumerator, FactorizedError>, prefix: &str) {
    match result {
        Err(FactorizedError::TextSyntax { detail, .. }) => {
            assert!(detail.starts_with(prefix), "{detail}")
        }
        other => panic!("expected TextSyntax({prefix}..), got {other:?}"),
    }
}

// ----------------------------------------------------------------------------
// The ports entry (evaluate.rs check_open_ports)

/// An OPEN spinor port (an open spinor chain: two bispinor slots no
/// projector closes) through the ports entry. Twin: the same bubble with a
/// Lorentz port only.
#[test]
fn refusal_ports_entry_open_spinor_port() {
    let vakint = engine();
    let p = slot(0, 1);
    let good = bubble(&format!("{Q}(0,{p})"), std::slice::from_ref(&p), None);
    ports_entry(&vakint, &good).unwrap();
    let mut ports = vec![p.clone(), bis(99, 0), bis(99, 1)];
    ports.sort();
    let bad = bubble(
        &format!(
            "{Q}(0,{p})*{GAMMA}({},{},{})*{Q}(0,{})",
            bis(99, 0),
            bis(99, 1),
            slot(0, 7),
            slot(0, 7)
        ),
        &ports,
        None,
    );
    bad.validate().unwrap();
    assert_ports(
        ports_entry(&vakint, &bad),
        "open port ",
        "is not a Lorentz slot",
    );
}

/// The CLOSED bispinor pair (the K6 quark record: a compact-slash gamma
/// projector closes the two bispinor ports) is accepted; without its
/// projector the fermion line through the two amputated legs is open.
#[test]
fn refusal_k6_quark_without_its_projector_is_an_open_line() {
    let vakint = engine();
    let text = read("tests/fixtures/wh_ports/k6-quark.input.txt");
    let good = FactorizedNumerator::from_text(&text).unwrap();
    assert!(good.projector.is_some());
    let settings = refpath::vakint_settings(good.loops() as i64 + 1);
    evaluate_factorized_ports(&vakint, &settings, &good, &Atom::num(1)).unwrap();
    let mut bad = good.clone();
    bad.projector = None;
    match bad.validate() {
        Err(FactorizedError::OpenLine { .. }) => {}
        other => panic!("expected OpenLine, got {other:?}"),
    }
}

/// An open Lorentz slot with an INTEGER index (a projector slot that is not
/// a port). Twin: the same open slot with a named index.
#[test]
fn refusal_ports_entry_integer_index() {
    let vakint = engine();
    let p = slot(0, 1);
    let factor = format!("{Q}(0,{p})");
    let named = format!("{MINK}({DIM},ports_refusal::{{}}::nu)");
    let good = bubble(
        &factor,
        std::slice::from_ref(&p),
        Some(&format!("{G}({p},{named})")),
    );
    ports_entry(&vakint, &good).unwrap();
    let bad = bubble(
        &factor,
        std::slice::from_ref(&p),
        Some(&format!("{G}({p},{MINK}({DIM},7))")),
    );
    bad.validate().unwrap();
    assert_ports(
        ports_entry(&vakint, &bad),
        "open port ",
        "has an integer index",
    );
}

// ----------------------------------------------------------------------------
// Port rules of a marked value (factorized.rs check_contraction)

/// A port that is not a slot.
#[test]
fn refusal_port_not_a_slot() {
    let vakint = engine();
    let p = slot(0, 1);
    let factor = format!("{Q}(0,{p})");
    let good = bubble(&factor, std::slice::from_ref(&p), None);
    ports_entry(&vakint, &good).unwrap();
    let bad = bubble(
        &factor,
        &[p.clone(), "ports_refusal::{}::x".to_string()],
        None,
    );
    assert_ports(
        ports_entry(&vakint, &bad),
        "ports_refusal::{}::x is not a Lorentz or spinor slot",
        "",
    );
}

/// A port listed twice.
#[test]
fn refusal_port_listed_twice() {
    let vakint = engine();
    let p = slot(0, 1);
    let factor = format!("{Q}(0,{p})");
    let good = bubble(&factor, std::slice::from_ref(&p), None);
    ports_entry(&vakint, &good).unwrap();
    let bad = bubble(&factor, &[p.clone(), p.clone()], None);
    assert_ports(ports_entry(&vakint, &bad), "a port is listed twice", "");
}

/// A free index of the projector that occurs in the factors but is not a
/// port. Twin: the same index listed as a port (the projector closes both).
#[test]
fn refusal_projector_index_in_the_factors_not_a_port() {
    let vakint = engine();
    let (a, b) = (slot(0, 1), slot(0, 2));
    let factor = format!("{Q}(0,{a})*{Q}(0,{b})");
    let projector = format!("{G}({a},{b})");
    let mut both = vec![a.clone(), b.clone()];
    both.sort();
    let good = bubble(&factor, &both, Some(&projector));
    ports_entry(&vakint, &good).unwrap();
    let bad = bubble(&factor, std::slice::from_ref(&a), Some(&projector));
    assert_ports(
        ports_entry(&vakint, &bad),
        "the projector's free index ",
        "occurs in the factors but is not a port",
    );
}

/// A port that occurs only in the projector, never in a factor (A-1).
/// Twin: the projector without that port.
#[test]
fn refusal_port_only_in_the_projector() {
    let vakint = engine();
    let (a, absent) = (slot(0, 1), slot(99, 0));
    let factor = format!("{Q}(0,{a})");
    let good = bubble(
        &factor,
        std::slice::from_ref(&a),
        Some(&format!("{Q}(0,{a})")),
    );
    ports_entry(&vakint, &good).unwrap();
    let mut ports = vec![a.clone(), absent.clone()];
    ports.sort();
    let bad = bubble(
        &factor,
        &ports,
        Some(&format!("{Q}(0,{a})*{Q}(0,{absent})")),
    );
    assert_ports(
        ports_entry(&vakint, &bad),
        "port ",
        "does not occur in the factors",
    );
}

// ----------------------------------------------------------------------------
// The v3 direct record (text.rs read_direct_section, the open section;
// factorized.rs validate_contract with original_cycles)

/// The given open-slot list differs from the slots the record leaves open
/// (one slot exchanged; same line count).
#[test]
fn refusal_v3_direct_wrong_open_slot() {
    let text = k5();
    FactorizedNumerator::from_text(&text).unwrap();
    let bad = edit(
        &text,
        &format!("open_slot\t{}\n", slot(1, 1)),
        &format!("open_slot\t{}\n", slot(1, 9)),
    );
    assert_text(FactorizedNumerator::parse_text(&bad), "open slots ");
}

/// The `open` count differs from its lines (one fewer declared).
#[test]
fn refusal_v3_direct_open_count() {
    let text = k5();
    let bad = edit(&text, "open\t2\n", "open\t1\n");
    assert_text(
        FactorizedNumerator::parse_text(&bad),
        "expected `prefactor`, found `open_slot`",
    );
}

/// The closure section of a direct record: `closure`, `taylor` must be
/// `none`, `softs` 0.
#[test]
fn refusal_v3_direct_closure_section() {
    let text = k5();
    FactorizedNumerator::from_text(&text).unwrap();
    let bad = edit(&text, "closure\tnone\n", "closure\tedge=2\tpower=1\tR\n");
    assert_text(
        FactorizedNumerator::parse_text(&bad),
        "kind direct with `closure` not none",
    );
    let bad = edit(&text, "taylor\tnone\n", "taylor\tdegree=0\n");
    assert_text(
        FactorizedNumerator::parse_text(&bad),
        "kind direct with `taylor` not none",
    );
    let bad = edit(&text, "softs\t0\n", "softs\t1\n");
    assert_text(
        FactorizedNumerator::parse_text(&bad),
        "kind direct with soft legs",
    );
}

/// A weight's `applied_by` other than vakint / beta. Twin: the same weight
/// applied by the beta.
#[test]
fn refusal_v3_direct_weight_applied_by() {
    let text = k5();
    let with = |by: &str| {
        more_lines(
            &edit(
                &text,
                "original_cycles\tfermion=0\tghost=1\n",
                &format!(
                    "original_cycles\tfermion=0\tghost=1\nweights\t1\nweight\tloop_sign\tapplied_by={by}\t-1\n"
                ),
            ),
            2,
        )
    };
    FactorizedNumerator::from_text(&with("beta")).unwrap();
    assert_text(
        FactorizedNumerator::parse_text(&with("nobody")),
        "applied_by=nobody",
    );
}

/// `depth` must be `none` or a number (one field).
#[test]
fn refusal_v3_direct_depth() {
    let text = k5();
    FactorizedNumerator::from_text(&text).unwrap();
    let bad = edit(&text, "depth\t2\n", "depth\t2\t3\n");
    assert_text(FactorizedNumerator::parse_text(&bad), "depth<TAB><n>");
    let bad = edit(&text, "depth\t2\n", "depth\ttwo\n");
    assert_text(FactorizedNumerator::parse_text(&bad), "bad depth \"two\"");
}

/// With `original_cycles` the weights are the caller's: loop_sign 1,
/// nf_power 0, prefactor_overall 1 (each broken alone).
#[test]
fn refusal_v3_direct_weights_with_original_cycles() {
    let text = k5();
    FactorizedNumerator::from_text(&text).unwrap();
    for (from, to) in [
        (
            "prefactor\tloop_sign=1\tnf_power=0\n",
            "prefactor\tloop_sign=-1\tnf_power=0\n",
        ),
        (
            "prefactor\tloop_sign=1\tnf_power=0\n",
            "prefactor\tloop_sign=1\tnf_power=1\n",
        ),
        ("prefactor_overall\t1\n", "prefactor_overall\t2\n"),
    ] {
        match FactorizedNumerator::from_text(&edit(&text, from, to)) {
            Err(FactorizedError::NotScalar { location, detail }) => {
                assert_eq!(location, "prefactor");
                assert!(
                    detail.starts_with("with original_cycles the weights are the caller's"),
                    "{detail}"
                );
            }
            other => panic!("{to}: {other:?}"),
        }
    }
}

/// The record and its insertions carry more cycles than `original_cycles`
/// (the k5-direct record carries one ghost cycle).
#[test]
fn refusal_v3_direct_cycles_above_original() {
    let text = k5();
    FactorizedNumerator::from_text(&text).unwrap();
    let bad = edit(
        &text,
        "original_cycles\tfermion=0\tghost=1\n",
        "original_cycles\tfermion=0\tghost=0\n",
    );
    match FactorizedNumerator::from_text(&bad) {
        Err(FactorizedError::Graph { detail }) => {
            assert!(
                detail.starts_with("the record and its insertions carry "),
                "{detail}"
            )
        }
        other => panic!("{other:?}"),
    }
}

/// Insertions on a direct record: the insertion vertex's own factor must be
/// 1; an `ins:` vertex carries exactly two boundary half-edges.
#[test]
fn refusal_v3_direct_insertion_rules() {
    let text = k5();
    let good = FactorizedNumerator::from_text(&text).unwrap();
    assert_eq!(good.insertions.len(), 1);
    let bad = edit(
        &text,
        "vertex\t2\trule=ins:scalar_insertion_0\t1\n",
        "vertex\t2\trule=ins:scalar_insertion_0\t2\n",
    );
    match FactorizedNumerator::from_text(&bad) {
        Err(FactorizedError::Insertion { vertex: 2, detail }) => {
            assert!(
                detail.starts_with("the vertex factor must be 1"),
                "{detail}"
            )
        }
        other => panic!("{other:?}"),
    }
    let bad = edit(
        &text,
        "insertion\t2\thedges=5,6\t",
        "insertion\t2\thedges=5\t",
    );
    match FactorizedNumerator::from_text(&bad) {
        Err(FactorizedError::Insertion { vertex: 2, detail }) => {
            assert!(
                detail.starts_with("a surrogate coefficient vertex (rule "),
                "{detail}"
            )
        }
        other => panic!("{other:?}"),
    }
}
