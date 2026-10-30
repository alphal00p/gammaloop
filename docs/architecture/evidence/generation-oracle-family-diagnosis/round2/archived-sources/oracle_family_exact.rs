//! Exact factor-preserving comparison for bounded oracle-family expression pairs.

use std::env;
use std::fs;

use serde::Serialize;
use serde_json::Value;
use symbolica::atom::{Atom, AtomCore};
use symbolica::parser::ParseSettings;

#[derive(Serialize)]
struct ResultRow {
    parsed_equal: bool,
    canonical_equal: bool,
    sign_rejected: bool,
    coefficient_rejected: bool,
    power_rejected: bool,
}

fn main() {
    let arguments = env::args().collect::<Vec<_>>();
    assert_eq!(arguments.len(), 3, "pairs.json results.json");
    let pairs: Vec<Value> = serde_json::from_slice(&fs::read(&arguments[1]).unwrap()).unwrap();
    let mut results = Vec::new();
    for pair in pairs {
        let historical = pair
            .get("historical")
            .and_then(Value::as_str)
            .expect("historical expression");
        let current = pair
            .get("current")
            .and_then(Value::as_str)
            .expect("current expression");
        let historical =
            Atom::parse(historical, "oracle_family", ParseSettings::default()).unwrap();
        let current = Atom::parse(current, "oracle_family", ParseSettings::default()).unwrap();
        let row = ResultRow {
            parsed_equal: historical == current,
            canonical_equal: historical.to_canonical_string() == current.to_canonical_string(),
            sign_rejected: -&historical != current,
            coefficient_rejected: &historical * Atom::num(2) != current,
            power_rejected: historical.pow(2) != current,
        };
        assert!(row.parsed_equal && row.canonical_equal);
        assert!(row.sign_rejected && row.coefficient_rejected && row.power_rejected);
        results.push(row);
    }
    fs::write(&arguments[2], serde_json::to_vec(&results).unwrap()).unwrap();
    println!("{} exact factorized expression pairs", results.len());
}
