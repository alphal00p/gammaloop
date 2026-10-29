//! Diagnostic exact tensor admission for the frozen grouping-frame captures.

use std::env;
use std::fs;
use std::path::Path;

use feynkit_graph::NumeratorInterface;
use feynkit_model::Model;
use serde_json::Value;
use symbolica::atom::Atom;
use symbolica::parser::ParseSettings;

fn diagram(path: &Path, name: &str) -> Value {
    let data: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    data["diagrams"]
        .as_array()
        .unwrap()
        .iter()
        .find(|diagram| diagram["name"] == name)
        .unwrap()
        .clone()
}

#[test]
#[ignore = "diagnostic frozen producer comparison"]
fn d65_prefactor_constructor_admission() {
    const OLD_Q: &str = "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q";
    const CURRENT_Q: &str =
        "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q";
    const OLD_VERTEX: &str = "gammalooprs::{spenso::index,spenso::index-label:v}::vertex";
    const CURRENT_VERTEX: &str = "gammalooprs::{spenso::index}::vertex";
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/frame-certificate");
    let old = diagram(&root.join("historical-frozen/comparison-qed.json"), "D65");
    let current = diagram(&root.join("current-frozen/comparison-qed.json"), "D65");
    let old = old["payload"]["numerator_prefactor"].as_str().unwrap();
    let current = current["payload"]["numerator_prefactor"].as_str().unwrap();
    // Bind only the presentation attributes removed from the same graph constructors.
    // The parser rejects duplicate symbol names with different tag sets in one process.
    assert_eq!(old.matches(OLD_Q).count(), 4);
    let old = old
        .replace(OLD_Q, CURRENT_Q)
        .replace(OLD_VERTEX, CURRENT_VERTEX);
    assert_eq!(old, current);
    let old = Atom::parse(&old, "frame certificate old", ParseSettings::default()).unwrap();
    let current = Atom::parse(
        current,
        "frame certificate current",
        ParseSettings::default(),
    )
    .unwrap();
    let old = NumeratorInterface::admit(old, None).unwrap();
    let current = NumeratorInterface::admit(current, None).unwrap();
    assert!(old.same_value(&current));
}

#[test]
#[ignore = "diagnostic effective-model fingerprint verifier"]
fn model_fingerprint_from_file() {
    let path = env::var("FEYNKIT_FRAME_MODEL_INPUT").unwrap();
    let definition = fs::read_to_string(path).unwrap();
    let model = Model::from_json(&definition).unwrap();
    println!("FRAME_MODEL_FINGERPRINT={}", model.fingerprint());
}
