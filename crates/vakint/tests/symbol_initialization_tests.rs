//! Public initialization must not parse a numeric literal as a Symbol.

use spenso::network::tags::SPENSO_TAG;
use symbolica::atom::{Atom, AtomCore};
use symbolica::domains::float::{Complex, Float};
use vakint::{Vakint, VakintSettings, symbols::S, vakint_parse};

#[test]
fn exported_macros_initialize_momentum_tags_in_a_fresh_process() {
    const MODE: &str = "VAKINT_COLD_MACRO_TEST";
    if let Ok(mode) = std::env::var(MODE) {
        // Each branch runs before any Vakint instance or explicit symbol
        // initialization, regardless of other tests' ordering or concurrency.
        let head = match mode.as_str() {
            "parse_default" => {
                vakint_parse!("k(1,101)*p(2,101)").unwrap();
                S.k
            }
            "parse_explicit_namespace" => {
                vakint_parse!("k(1,101)*p(2,101)", "vakint").unwrap();
                S.p
            }
            "symbol_unqualified" => vakint::vakint_symbol!("k"),
            "symbol_qualified" => vakint::vakint_symbol!("vakint::p"),
            _ => panic!("unknown cold initialization case"),
        };
        assert!(head.has_tag(&SPENSO_TAG.tensor));
        assert!(head.has_tag(&SPENSO_TAG.rank1));
        assert_eq!(S.mom, S.k);
        Vakint::initialize_vakint_symbols();
        Vakint::initialize_vakint_symbols();
        Vakint::new().unwrap();
        return;
    }
    for mode in [
        "parse_default",
        "parse_explicit_namespace",
        "symbol_unqualified",
        "symbol_qualified",
    ] {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "exported_macros_initialize_momentum_tags_in_a_fresh_process",
                "--nocapture",
            ])
            .env(MODE, mode)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "cold initialization case {mode} failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }
}

#[test]
fn public_initialization_preserves_namespaced_imaginary_placeholder() {
    Vakint::initialize_vakint_symbols();
    Vakint::initialize_vakint_symbols();
    assert_eq!(S.cmplx_i.get_namespace(), "vakint");
    let placeholder = Atom::var(S.cmplx_i);
    assert_eq!(vakint_parse!("vakint::𝑖").unwrap(), placeholder);

    let numeric_i = vakint_parse!("𝑖").unwrap();
    assert_ne!(numeric_i, placeholder);
    assert_eq!(numeric_i.pow(Atom::num(2)), Atom::num(-1));

    // The historical FORM-facing symbolic placeholder still becomes the
    // numerical imaginary unit only at the existing evaluator boundary.
    let settings = VakintSettings::default();
    let (_, complex) =
        Vakint::get_constants_map(&settings, &Default::default(), &Default::default(), None)
            .unwrap();
    let precision = settings.get_binary_precision();
    assert_eq!(
        complex.get(&placeholder),
        Some(&Complex::new(
            Float::with_val(precision, 0),
            Float::with_val(precision, 1),
        )),
    );
}

#[test]
fn momentum_heads_retain_tensor_tags_and_form_identities() {
    Vakint::initialize_vakint_symbols();
    for (head, name) in [(S.k, "k"), (S.p, "p")] {
        assert_eq!(head.get_namespace(), "vakint");
        assert_eq!(head.get_stripped_name(), name);
        assert!(head.has_tag(&SPENSO_TAG.tensor));
        assert!(head.has_tag(&SPENSO_TAG.rank1));
        assert!(!S.should_symbol_be_escaped_in_form(&head));
    }
    assert_eq!(S.mom, S.k);
    let indexed = vakint_parse!("k(1,101)*p(1,101)").unwrap();
    let reparsed = vakint_parse!(&indexed.to_canonical_string()).unwrap();
    assert_eq!(indexed, reparsed);
    assert!(reparsed.get_all_symbols(true).contains(&S.k));
    assert!(reparsed.get_all_symbols(true).contains(&S.p));
}
