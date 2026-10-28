use symbolica::parser::ParseSettings;
use serde_json::Value;
use std::env;
use std::fs;
use symbolica::atom::{Atom, AtomCore};
fn main() {
    let paths = env::args().collect::<Vec<_>>();
    let pairs: Vec<Value> = serde_json::from_slice(&fs::read(&paths[1]).unwrap()).unwrap();
    let mut results = Vec::new();
    for pair in pairs {
        let before = Atom::parse(pair["baseline"].as_str().unwrap(), "frozen_local", ParseSettings::default()).unwrap();
        let after = Atom::parse(pair["current"].as_str().unwrap(), "frozen_local", ParseSettings::default()).unwrap();
        let parsed_equal = before == after;
        let factored_equal = before.collect_factors() == after.collect_factors();
        assert!(!before.is_zero() && !after.is_zero());
        assert_ne!(before, -&after, "sign corruption must be rejected");
        assert_ne!(&before * Atom::num(2), after, "coefficient corruption must be rejected");
        results.push((parsed_equal, factored_equal));
    }
    fs::write(&paths[2], serde_json::to_vec(&results).unwrap()).unwrap();
    println!("{} exact local pairs: parsed {}, factored {}", results.len(), results.iter().filter(|r| r.0).count(), results.iter().filter(|r| r.1).count());
}
