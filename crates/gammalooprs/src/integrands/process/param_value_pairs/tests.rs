use std::env;
use std::fs::{self, File};
use std::io::Cursor;
use std::path::Path;
use std::process::{self, Command};
use std::sync::Arc;

use bincode::config;
use symbolica::coefficient::Coefficient;
use symbolica::domains::finite_field::{FiniteFieldCore, Zp64};
use symbolica::domains::float::{Complex as SymbolicaComplex, Float as SymbolicaFloat};
use symbolica::state::State;
use symbolica::symbol;

use super::*;
use crate::GammaLoopContextContainer;
use crate::utils::load_generic_model;

impl ParamValuePairs {
    fn coefficient_fixture(domain: &str, noise: bool) -> Atom {
        if domain.starts_with("field") {
            let field = Zp64::new(if noise { 23 } else { 7 });
            let element = field.to_element(if domain == "field-zero" { 0 } else { 3 });
            let mut atom = Atom::Zero;
            atom.to_num(Coefficient::from_finite_field(field, element));
            return atom;
        }
        let namespaces = if noise {
            [
                "parameter_domain_noise_left",
                "parameter_domain_noise_right",
            ]
        } else {
            ["parameter_domain_left", "parameter_domain_right"]
        };
        let left = symbol!(&format!("{}::x", namespaces[0]));
        let right = symbol!(&format!("{}::x", namespaces[1]));
        let expression = Atom::parse(
            format!(
                "({}::x+2*{}::x+7)/({}::x*{}::x+3)",
                namespaces[0], namespaces[1], namespaces[0], namespaces[1]
            ),
            GAMMALOOP_NAMESPACE,
            ParseSettings::default(),
        )
        .unwrap();
        let atom = expression.set_coefficient_ring(Arc::new(vec![left.into(), right.into()]));
        if domain == "polynomial-zero" {
            let AtomView::Num(number) = atom.as_view() else {
                panic!("Expected polynomial coefficient")
            };
            let CoefficientView::RationalPolynomial(polynomial) = number.get_coeff_view() else {
                panic!("Expected polynomial domain")
            };
            let mut polynomial = polynomial.deserialize();
            polynomial.numerator =
                MultivariatePolynomial::new(&Z, None, Arc::clone(polynomial.numerator.variables()));
            let mut zero = Atom::Zero;
            zero.to_num(Coefficient::RationalPolynomial(polynomial));
            zero
        } else {
            atom
        }
    }

    fn assert_coefficient_fresh_process_roundtrip(domain: &str, test: &str) {
        const MODE: &str = "GAMMALOOP_PARAMETER_DOMAIN_MODE";
        const DIRECTORY: &str = "GAMMALOOP_PARAMETER_DOMAIN_DIRECTORY";
        if let Ok(mode) = env::var(MODE) {
            let directory = env::var(DIRECTORY).unwrap();
            let directory = Path::new(&directory);
            let model = load_generic_model("scalars");
            if mode == "read" {
                // Register unrelated fields, variable lists and symbols first.
                // Native indices must change without changing binding identity.
                let _ = Self::coefficient_fixture("field", true);
                let _ = Self::coefficient_fixture("polynomial", true);
                let _ = symbol!("parameter_domain_right::x", "parameter_domain_left::x");
                let _ = symbol!("parameter_domain::b", "parameter_domain::a");
                let state_map = State::import(
                    &mut Cursor::new(fs::read(directory.join("state")).unwrap()),
                    None,
                )
                .unwrap();
                let encoded = fs::read(directory.join("bindings")).unwrap();
                let (restored, consumed): (Self, _) = bincode::decode_from_slice_with_context(
                    &encoded,
                    config::standard(),
                    GammaLoopContextContainer {
                        state_map: &state_map,
                        model: &model,
                    },
                )
                .unwrap();
                assert_eq!(consumed, encoded.len());
                assert_eq!(restored.value_range, 11..13);
                let expected = Self::coefficient_fixture(domain, false);
                assert_eq!(restored.params[0], expected);
                assert_eq!(
                    restored.params[1],
                    symbol!("parameter_domain::binding").call_args([&expected])
                );
                match (domain, restored.params[0].as_view()) {
                    (domain, AtomView::Num(number)) if domain.starts_with("field") => {
                        assert!(matches!(
                            number.get_coeff_view(),
                            CoefficientView::FiniteField(..)
                        ));
                        assert_eq!(
                            expected.to_canonical_string(),
                            if domain == "field-zero" {
                                "[0%7]"
                            } else {
                                "[3%7]"
                            }
                        );
                    }
                    (domain, AtomView::Num(number)) if domain.starts_with("polynomial") => {
                        assert!(matches!(
                            number.get_coeff_view(),
                            CoefficientView::RationalPolynomial(..)
                        ));
                    }
                    _ => panic!("The native coefficient domain must survive"),
                }
                if domain == "field-antisymmetric" {
                    let (a, b) = symbol!("parameter_domain::a", "parameter_domain::b");
                    let anti = symbol!("parameter_domain::anti"; Antisymmetric);
                    let field = Zp64::new(7);
                    let element = field.to_element(4);
                    let expected = Atom::num(Coefficient::from_finite_field(field, element))
                        * anti.call_args([Atom::var(b), Atom::var(a)]);
                    assert_eq!(
                        restored.params[2], expected,
                        "Antisymmetric swap must negate within GF(7)"
                    );
                }
                assert_ne!(
                    fs::read(directory.join("native")).unwrap(),
                    restored.params[0].as_view().get_data()
                );
                assert_eq!(
                    bincode::encode_to_vec(&restored, config::standard()).unwrap(),
                    encoded,
                    "Remapped domains must preserve the exact resume fingerprint"
                );
                return;
            }
            let _ = symbol!("parameter_domain::a", "parameter_domain::b");
            let coefficient = Self::coefficient_fixture(domain, false);
            let mut pairs = Self {
                value_range: 11..13,
                params: vec![
                    coefficient.clone(),
                    symbol!("parameter_domain::binding").call_args([&coefficient]),
                ],
            };
            if domain == "field-antisymmetric" {
                let (a, b) = symbol!("parameter_domain::a", "parameter_domain::b");
                let anti = symbol!("parameter_domain::anti"; Antisymmetric);
                pairs
                    .params
                    .push(&coefficient * anti.call_args([Atom::var(a), Atom::var(b)]));
            }
            let encoded = bincode::encode_to_vec(&pairs, config::standard()).unwrap();
            fs::write(directory.join("bindings"), encoded).unwrap();
            fs::write(directory.join("native"), coefficient.as_view().get_data()).unwrap();
            State::export(&mut File::create(directory.join("state")).unwrap()).unwrap();
            return;
        }
        let directory =
            env::temp_dir().join(format!("gammaloop-parameter-{domain}-{}", process::id()));
        fs::create_dir(&directory).unwrap();
        for mode in ["write", "read"] {
            let output = Command::new(env::current_exe().unwrap())
                .args(["--exact", test, "--nocapture"])
                .env(MODE, mode)
                .env(DIRECTORY, &directory)
                .output()
                .unwrap();
            if !output.status.success() {
                fs::remove_dir_all(&directory).unwrap();
                panic!(
                    "{domain} {mode} child failed:\n{}\n{}",
                    String::from_utf8_lossy(&output.stdout),
                    String::from_utf8_lossy(&output.stderr)
                );
            }
        }
        fs::remove_dir_all(directory).unwrap();
    }
}

#[test]
fn parameter_pairs_field_zero_roundtrip_in_fresh_process() {
    ParamValuePairs::assert_coefficient_fresh_process_roundtrip(
        "field-zero",
        "integrands::process::param_value_pairs::tests::parameter_pairs_field_zero_roundtrip_in_fresh_process",
    );
}

#[test]
fn parameter_pairs_polynomial_zero_roundtrip_in_fresh_process() {
    ParamValuePairs::assert_coefficient_fresh_process_roundtrip(
        "polynomial-zero",
        "integrands::process::param_value_pairs::tests::parameter_pairs_polynomial_zero_roundtrip_in_fresh_process",
    );
}

#[test]
fn parameter_pairs_field_antisymmetric_roundtrip_in_fresh_process() {
    ParamValuePairs::assert_coefficient_fresh_process_roundtrip(
        "field-antisymmetric",
        "integrands::process::param_value_pairs::tests::parameter_pairs_field_antisymmetric_roundtrip_in_fresh_process",
    );
}

#[test]
fn parameter_pairs_finite_field_roundtrip_in_fresh_process() {
    ParamValuePairs::assert_coefficient_fresh_process_roundtrip(
        "field",
        "integrands::process::param_value_pairs::tests::parameter_pairs_finite_field_roundtrip_in_fresh_process",
    );
}

#[test]
fn parameter_pairs_rational_polynomial_roundtrip_in_fresh_process() {
    ParamValuePairs::assert_coefficient_fresh_process_roundtrip(
        "polynomial",
        "integrands::process::param_value_pairs::tests::parameter_pairs_rational_polynomial_roundtrip_in_fresh_process",
    );
}
#[test]
fn parameter_pairs_preserve_finite_field_signs_and_domains() {
    let coefficient = |prime, value| {
        let field = Zp64::new(prime);
        let element = field.to_element(value);
        Atom::num(Coefficient::from_finite_field(field, element))
    };
    let values = [coefficient(7, 3), coefficient(7, 4), coefficient(11, 3)];
    assert_roundtrip(&values.iter().collect());
    let encodings: Vec<_> = values
        .iter()
        .map(|value| {
            bincode::encode_to_vec(ParamValuePairs::from_iter([value]), config::standard()).unwrap()
        })
        .collect();
    assert_ne!(encodings[0], encodings[1]);
    assert_ne!(encodings[0], encodings[2]);
    assert_ne!(encodings[1], encodings[2]);
}

#[test]
fn parameter_pairs_preserve_polynomial_signs() {
    let polynomial = ParamValuePairs::coefficient_fixture("polynomial", false);
    assert_roundtrip(&[polynomial.clone(), -&polynomial].into_iter().collect());
}

#[test]
fn parameter_pairs_roundtrip_preserves_namespace_order_and_content() {
    let (left, right) = symbol!("pair_namespace_left::p", "pair_namespace_right::p");
    let mut pairs = ParamValuePairs {
        value_range: 4..6,
        params: vec![
            left.call_args([Atom::num(0)]),
            right.call_args([Atom::num(0)]),
        ],
    };
    let encoded = bincode::encode_to_vec(&pairs, config::standard()).unwrap();
    let mut state = Vec::new();
    State::export(&mut state).unwrap();
    let state_map = State::import(&mut Cursor::new(state), None).unwrap();
    let model = load_generic_model("scalars");
    let (restored, consumed): (ParamValuePairs, _) = bincode::decode_from_slice_with_context(
        &encoded,
        config::standard(),
        GammaLoopContextContainer {
            state_map: &state_map,
            model: &model,
        },
    )
    .unwrap();
    assert_eq!(consumed, encoded.len());
    assert_eq!(restored.value_range, pairs.value_range);
    assert_eq!(restored.params, pairs.params);
    assert_eq!(
        bincode::encode_to_vec(&restored, config::standard()).unwrap(),
        encoded
    );

    pairs.params.reverse();
    assert_ne!(
        bincode::encode_to_vec(&pairs, config::standard()).unwrap(),
        encoded
    );
    pairs.params.reverse();
    pairs.value_range = 5..7;
    assert_ne!(
        bincode::encode_to_vec(&pairs, config::standard()).unwrap(),
        encoded
    );
    pairs.value_range = 4..6;
    pairs.params[0] = right.call_args([Atom::num(0)]);
    assert_ne!(
        bincode::encode_to_vec(&pairs, config::standard()).unwrap(),
        encoded
    );
    pairs.params[0] = left.call_args([Atom::num(1)]);
    assert_ne!(
        bincode::encode_to_vec(&pairs, config::standard()).unwrap(),
        encoded
    );
}

fn assert_roundtrip(pairs: &ParamValuePairs) -> Vec<u8> {
    let encoded = bincode::encode_to_vec(pairs, config::standard()).unwrap();
    let mut state = Vec::new();
    State::export(&mut state).unwrap();
    let state_map = State::import(&mut Cursor::new(state), None).unwrap();
    let model = load_generic_model("scalars");
    let (restored, consumed): (ParamValuePairs, _) = bincode::decode_from_slice_with_context(
        &encoded,
        config::standard(),
        GammaLoopContextContainer {
            state_map: &state_map,
            model: &model,
        },
    )
    .unwrap();
    assert_eq!(consumed, encoded.len());
    assert_eq!(restored.value_range, pairs.value_range);
    for (before, after) in pairs.params.iter().zip(&restored.params) {
        // Equality alone does not establish component precision. Native
        // bytes include both float precisions, including imaginary zero.
        assert_eq!(
            before.as_view().get_data(),
            after.as_view().get_data(),
            "{before} != {after}"
        );
    }
    assert_eq!(
        bincode::encode_to_vec(&restored, config::standard()).unwrap(),
        encoded
    );
    encoded
}

#[test]
fn parameter_pairs_preserve_precision_and_registration_independent_symmetry() {
    const CHILD: &str = "GAMMALOOP_PARAMETER_PAIR_REGISTRATION_ORDER";
    const TEST: &str = "integrands::process::param_value_pairs::tests::parameter_pairs_preserve_precision_and_registration_independent_symmetry";
    const MARKER: &str = "PORTABLE_PARAMETER_BYTES=";
    if let Ok(order) = env::var(CHILD) {
        let names = if order == "forward" {
            ["a", "b", "c"]
        } else {
            ["c", "b", "a"]
        };
        for name in names {
            let _ = symbol!(&format!("pair_portability::{name}"));
        }
        let (a, b) = symbol!("pair_portability::a", "pair_portability::b");
        let symmetric = symbol!("pair_portability::symmetric"; Symmetric);
        let antisymmetric = symbol!("pair_portability::antisymmetric"; Antisymmetric);
        let linear = symbol!("pair_portability::linear"; Linear);
        let ordered = symbol!("pair_portability::ordered");
        let a = Atom::var(a);
        let b = Atom::var(b);
        let anti = antisymmetric.call_args([&a, &b]);
        assert_eq!(anti, -antisymmetric.call_args([&b, &a]));
        let exact = Atom::parse(
            "123456789012345678901234567890123456789/99887766554433221100998877665544332211",
            GAMMALOOP_NAMESPACE,
            ParseSettings::default(),
        )
        .unwrap();
        let real = SymbolicaFloat::with_val(173, 1) / SymbolicaFloat::with_val(173, 3);
        let float = Atom::num(SymbolicaComplex::new(
            real,
            SymbolicaFloat::with_val(211, 0),
        ));
        let complex = Atom::num(SymbolicaComplex::new(
            SymbolicaFloat::with_val(157, 2) / SymbolicaFloat::with_val(157, 7),
            SymbolicaFloat::with_val(199, 5) / SymbolicaFloat::with_val(199, 11),
        ));
        // The preferred placeholder namespace is deliberately occupied by
        // an original symbol. It must survive without becoming a leaf.
        let collision = Atom::var(symbol!("gammaloop_parameter_float_0::original"));
        let pairs: ParamValuePairs = [
            exact,
            float.clone(),
            -&float,
            complex.clone(),
            -&complex,
            symmetric.call_args([&a, &b]),
            anti.clone(),
            antisymmetric.call_args([&anti, &a]),
            anti.clone().pow(2),
            linear.call_args([&anti]),
            ordered.call_args([&a, &b]),
            ordered.call_args([&b, &a]),
            symmetric.call_args([&float, &complex]),
            antisymmetric.call_args([&float, &complex]),
            ordered.call_args([&collision, &float]),
        ]
        .into_iter()
        .collect();
        let bytes = assert_roundtrip(&pairs);
        let mut hex = String::new();
        for byte in bytes {
            write!(hex, "{byte:02x}").unwrap();
        }
        println!("{MARKER}{hex}");

        // An unrelated ambient definition cannot silently lend semantic
        // attributes to an internal leaf name. Fail explicitly rather
        // than choosing an ambient-dependent alternative wire name.
        let new_float = Atom::num(SymbolicaComplex::new(
            SymbolicaFloat::with_val(149, 9),
            SymbolicaFloat::with_val(223, 0),
        ));
        let negative = -&new_float;
        let representative = if new_float.as_view().get_data() <= negative.as_view().get_data() {
            &new_float
        } else {
            &negative
        };
        let mut conflicting_name = "gammaloop_parameter_float_0::n".to_owned();
        for byte in representative.as_view().get_data() {
            write!(conflicting_name, "{byte:02x}").unwrap();
        }
        let _ = symbol!(&conflicting_name; Symmetric);
        assert!(PortableParameter::from_atom(&new_float).is_err());
        return;
    }
    let mut outputs = Vec::new();
    for order in ["forward", "reverse"] {
        let output = Command::new(env::current_exe().unwrap())
            .args(["--exact", TEST, "--nocapture"])
            .env(CHILD, order)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "child {order} failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout).unwrap();
        outputs.push(
            stdout
                .lines()
                .find_map(|line| line.split_once(MARKER).map(|(_, bytes)| bytes.to_owned()))
                .expect("child must execute the roundtrip test"),
        );
    }
    assert_eq!(
        outputs[0], outputs[1],
        "registration order must not change parameter identity"
    );
}
