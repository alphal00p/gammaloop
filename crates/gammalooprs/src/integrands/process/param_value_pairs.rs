//! Ordered parameter bindings with portable symbolic identity in full-state archives.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;
use std::num::ParseIntError;
use std::ops::Range;
use std::sync::Arc;

use bincode::config;
use bincode::de::Decoder;
use bincode::enc::Encoder;
use bincode::error::{DecodeError, EncodeError};
use bincode::{Decode, Encode};
use symbolica::atom::SymbolBuilder;
use symbolica::coefficient::{Coefficient, CoefficientView};
use symbolica::domains::finite_field::{FiniteFieldCore, Zp64};
use symbolica::domains::integer::{Integer, IntegerRing, Z};
use symbolica::domains::rational_polynomial::RationalPolynomial;
use symbolica::parser::ParseSettings;
use symbolica::poly::PolyVariable;
use symbolica::poly::polynomial::MultivariatePolynomial;
use symbolica::prelude::{Atom, AtomCore, AtomOrView, AtomView, Symbol};
use symbolica::wrap_symbol;

use crate::{GAMMALOOP_NAMESPACE, GammaLoopContext};

/// Parameters in binding-slot order and the corresponding range of numeric values.
#[derive(Clone, Default)]
pub struct ParamValuePairs {
    pub value_range: Range<usize>,
    pub params: Vec<Atom>,
}

/// Canonical symbolic text with lossless, domain-aware coefficient leaves.
///
/// Symbolica's canonical printer handles full namespaces, symbol attributes and
/// antisymmetric permutation signs. Float printing loses precision, and printed
/// finite-field and rational-polynomial coefficients cannot be parsed as atoms.
/// Leaf identities therefore preserve exact values and domains without using
/// process-local field, symbol or variable-list indices.
#[derive(bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
struct PortableParameter {
    expression: String,
    coefficients: Vec<(String, PortableCoefficient)>,
}

#[derive(bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
enum PortableCoefficient {
    Native(Atom),
    FiniteField {
        value: u64,
        prime: u64,
    },
    Polynomial {
        variables: Vec<PortableVariable>,
        numerator: Vec<(String, Vec<u16>)>,
        denominator: Vec<(String, Vec<u16>)>,
    },
}

#[derive(bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = GammaLoopContext)]
enum PortableVariable {
    Symbol(String),
    Function(String, Box<PortableParameter>),
    Power(Box<PortableParameter>),
    // The pinned kernel's private, index-only variables contain no symbol IDs.
    // Preserve them through its public codec without accessing hidden variants.
    Internal(PolyVariable),
}

impl PortableVariable {
    fn from_variable(variable: &PolyVariable) -> Result<Self, EncodeError> {
        Ok(match variable {
            PolyVariable::Symbol(symbol) => Self::Symbol(Atom::var(*symbol).to_canonical_string()),
            PolyVariable::Function(symbol, expression) => Self::Function(
                Atom::var(*symbol).to_canonical_string(),
                Box::new(PortableParameter::from_atom(expression)?),
            ),
            PolyVariable::Power(expression) => {
                Self::Power(Box::new(PortableParameter::from_atom(expression)?))
            }
            variable => Self::Internal(variable.clone()),
        })
    }

    fn symbol(text: &str) -> Result<Symbol, DecodeError> {
        match PortableParameter::parse(text)?.as_view() {
            AtomView::Var(variable) => Ok(variable.get_symbol()),
            _ => Err(DecodeError::Other(
                "Polynomial variable symbol is not a symbol",
            )),
        }
    }

    fn into_variable(self) -> Result<PolyVariable, DecodeError> {
        Ok(match self {
            Self::Symbol(symbol) => PolyVariable::Symbol(Self::symbol(&symbol)?),
            Self::Function(symbol, expression) => {
                PolyVariable::Function(Self::symbol(&symbol)?, expression.into_atom()?)
            }
            Self::Power(expression) => PolyVariable::Power(expression.into_atom()?),
            Self::Internal(variable) => match variable {
                PolyVariable::Symbol(_) | PolyVariable::Function(..) | PolyVariable::Power(_) => {
                    return Err(DecodeError::Other(
                        "Internal polynomial variable contains symbolic identity",
                    ));
                }
                variable => variable,
            },
        })
    }
}

impl PortableCoefficient {
    fn terms(polynomial: &MultivariatePolynomial<IntegerRing, u16>) -> Vec<(String, Vec<u16>)> {
        (0..polynomial.nterms())
            .map(|index| {
                (
                    polynomial
                        .coefficient_back(polynomial.nterms() - index - 1)
                        .to_string(),
                    polynomial.exponents(index).to_vec(),
                )
            })
            .collect()
    }

    fn from_atom(atom: Atom) -> Result<Self, EncodeError> {
        let AtomView::Num(number) = atom.as_view() else {
            return Err(EncodeError::Other("Parameter coefficient is not numeric"));
        };
        match number.get_coeff_view() {
            CoefficientView::Float(..)
            | CoefficientView::Indeterminate
            | CoefficientView::Infinity(_) => Ok(Self::Native(atom)),
            CoefficientView::FiniteField(..) => {
                // The public printer exposes the field modulus; the native index
                // and its lookup are private to Symbolica. Do not parse this as an atom.
                let text = atom.to_canonical_string();
                let (value, prime) = text
                    .strip_prefix('[')
                    .and_then(|text| text.strip_suffix(']'))
                    .and_then(|text| text.split_once('%'))
                    .ok_or(EncodeError::Other(
                        "Invalid finite-field coefficient spelling",
                    ))?;
                Ok(Self::FiniteField {
                    value: value.parse().map_err(|error: ParseIntError| {
                        EncodeError::OtherString(error.to_string())
                    })?,
                    prime: prime.parse().map_err(|error: ParseIntError| {
                        EncodeError::OtherString(error.to_string())
                    })?,
                })
            }
            CoefficientView::RationalPolynomial(polynomial) => {
                let polynomial = polynomial.deserialize();
                Ok(Self::Polynomial {
                    variables: polynomial
                        .numerator
                        .variables()
                        .iter()
                        .map(PortableVariable::from_variable)
                        .collect::<Result<_, _>>()?,
                    numerator: Self::terms(&polynomial.numerator),
                    denominator: Self::terms(&polynomial.denominator),
                })
            }
            _ => Err(EncodeError::Other(
                "Exact rational coefficients use canonical expression text",
            )),
        }
    }

    fn identity(&self) -> Result<Cow<'_, [u8]>, EncodeError> {
        match self {
            Self::Native(atom) => Ok(Cow::Borrowed(atom.as_view().get_data())),
            _ => bincode::encode_to_vec(self, config::standard()).map(Cow::Owned),
        }
    }

    fn polynomial(
        terms: Vec<(String, Vec<u16>)>,
        variables: Arc<Vec<PolyVariable>>,
    ) -> Result<MultivariatePolynomial<IntegerRing, u16>, DecodeError> {
        let mut polynomial = MultivariatePolynomial::new(&Z, Some(terms.len()), variables);
        for (coefficient, exponents) in terms {
            if exponents.len() != polynomial.nvars() {
                return Err(DecodeError::Other(
                    "Polynomial coefficient has an invalid exponent shape",
                ));
            }
            polynomial.append_monomial(
                coefficient
                    .parse::<Integer>()
                    .map_err(|error| DecodeError::OtherString(error.to_string()))?,
                &exponents,
            );
        }
        Ok(polynomial)
    }

    fn into_atom(self) -> Result<Atom, DecodeError> {
        match self {
            Self::Native(atom) => {
                if !matches!(atom.as_view(), AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Float(..) | CoefficientView::Indeterminate | CoefficientView::Infinity(_)))
                {
                    return Err(DecodeError::Other(
                        "Native parameter leaf contains state-dependent identity",
                    ));
                }
                Ok(atom)
            }
            Self::FiniteField { value, prime } => {
                if prime < 3 || prime.is_multiple_of(2) || value >= prime {
                    return Err(DecodeError::Other(
                        "Invalid finite-field parameter value or modulus",
                    ));
                }
                let field = Zp64::new(prime);
                let element = field.to_element(value);
                let mut atom = Atom::Zero;
                atom.to_num(Coefficient::from_finite_field(field, element));
                Ok(atom)
            }
            Self::Polynomial {
                variables,
                numerator,
                denominator,
            } => {
                let variables = Arc::new(
                    variables
                        .into_iter()
                        .map(PortableVariable::into_variable)
                        .collect::<Result<_, _>>()?,
                );
                let numerator = Self::polynomial(numerator, Arc::clone(&variables))?;
                let denominator = Self::polynomial(denominator, variables)?;
                if denominator.is_zero() {
                    return Err(DecodeError::Other(
                        "Polynomial parameter denominator is zero",
                    ));
                }
                let mut atom = Atom::Zero;
                atom.to_num(Coefficient::RationalPolynomial(RationalPolynomial {
                    numerator,
                    denominator,
                }));
                Ok(atom)
            }
        }
    }
}

impl PortableParameter {
    fn from_atom(parameter: &Atom) -> Result<Self, EncodeError> {
        let mut namespaces = BTreeSet::new();
        parameter.visitor(&mut |part| {
            if let Some(symbol) = part.get_symbol() {
                namespaces.insert(symbol.get_namespace());
            }
            true
        });
        // Arbitrary user parameters may themselves use our preferred namespace.
        // Choose the first absent namespace rather than reserving user symbols.
        let namespace = (0..)
            .map(|index| format!("gammaloop_parameter_float_{index}"))
            .find(|candidate| !namespaces.contains(candidate.as_str()))
            .expect("Finitely many source namespaces leave an available placeholder namespace");
        let mut coefficients = BTreeMap::new();
        let mut symbol_error = None;
        let skeleton = parameter.replace_map(|part, _, out| {
            let AtomView::Num(number) = part else {
                return;
            };
            if matches!(
                number.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            ) {
                return;
            }
            let positive = part.to_owned();
            // Canonicalize the sign orbit so normalizing an antisymmetric
            // function can cancel a minus sign without changing leaf identity.
            let representation = (|| {
                // Negate in the coefficient domain: rational -1 times a
                // finite-field atom is intentionally indeterminate in Symbolica.
                let mut negative = Atom::Zero;
                negative.to_num(-number.get_coeff_view().to_owned());
                let positive = PortableCoefficient::from_atom(positive)?;
                let negative = PortableCoefficient::from_atom(negative)?;
                let negate = positive.identity()?.as_ref() > negative.identity()?.as_ref();
                Ok::<_, EncodeError>(if !negate {
                    (positive, false)
                } else {
                    (negative, true)
                })
            })();
            let (representative, negate) = match representation {
                Ok(value) => value,
                Err(error) => {
                    symbol_error = Some(error);
                    return;
                }
            };
            let kind = if matches!(representative, PortableCoefficient::Native(_)) {
                "n"
            } else {
                "c"
            };
            let mut name = format!("{namespace}::{kind}");
            let identity = match representative.identity() {
                Ok(identity) => identity,
                Err(error) => {
                    symbol_error = Some(error);
                    return;
                }
            };
            for byte in identity.iter() {
                write!(name, "{byte:02x}").expect("Formatting into a String cannot fail");
            }
            drop(identity);
            let placeholder = match SymbolBuilder::new(wrap_symbol!(&name))
                .with_attributes(&[])
                .build()
            {
                Ok(symbol) => Atom::var(symbol),
                Err(error) => {
                    symbol_error = Some(EncodeError::OtherString(error.to_string()));
                    return;
                }
            };
            coefficients.insert(name, representative);
            **out = if negate { -placeholder } else { placeholder };
        });
        if let Some(error) = symbol_error {
            return Err(error);
        }
        Ok(Self {
            expression: skeleton.to_canonical_string(),
            coefficients: coefficients.into_iter().collect(),
        })
    }

    fn parse(expression: &str) -> Result<Atom, DecodeError> {
        Atom::parse(expression, GAMMALOOP_NAMESPACE, ParseSettings::default())
            .map_err(DecodeError::OtherString)
    }

    fn into_atom(self) -> Result<Atom, DecodeError> {
        let expression = Self::parse(&self.expression)?;
        let coefficients: BTreeMap<_, _> = self
            .coefficients
            .into_iter()
            .map(|(name, coefficient)| Ok((name, coefficient.into_atom()?)))
            .collect::<Result<_, DecodeError>>()?;
        Ok(Self::restore_coefficients(&expression, &coefficients))
    }

    fn restore_coefficients(expression: &Atom, coefficients: &BTreeMap<String, Atom>) -> Atom {
        expression.replace_map(|part, _, out| {
            // Parsing can permute antisymmetric arguments and introduce a
            // rational minus sign. Absorb that sign in the field before normalizing
            // the product; mixing rational -1 with the field loses its value.
            if let AtomView::Mul(product) = part {
                let field = product.iter().find_map(|factor| {
                    let AtomView::Var(variable) = factor else {
                        return None;
                    };
                    let value = coefficients.get(variable.get_symbol().get_name())?;
                    let AtomView::Num(number) = value.as_view() else {
                        return None;
                    };
                    matches!(number.get_coeff_view(), CoefficientView::FiniteField(..))
                        .then_some((factor, number.get_coeff_view().to_owned()))
                });
                if let Some((placeholder, coefficient)) = field
                    && product
                        .iter()
                        .any(|factor| factor == Atom::num(-1).as_view())
                {
                    let mut restored = Atom::Zero;
                    restored.to_num(-coefficient);
                    for factor in product.iter() {
                        if factor != placeholder && factor != Atom::num(-1).as_view() {
                            restored *=
                                Self::restore_coefficients(&factor.to_owned(), coefficients);
                        }
                    }
                    **out = restored;
                    return;
                }
            }
            if let AtomView::Var(variable) = part
                && let Some(value) = coefficients.get(variable.get_symbol().get_name())
            {
                **out = value.clone();
            }
        })
    }
}

impl Encode for ParamValuePairs {
    fn encode<E: Encoder>(&self, encoder: &mut E) -> Result<(), EncodeError> {
        self.value_range.encode(encoder)?;
        self.params
            .iter()
            .map(PortableParameter::from_atom)
            .collect::<Result<Vec<_>, _>>()?
            .encode(encoder)
    }
}

impl<C: GammaLoopContext> Decode<C> for ParamValuePairs {
    fn decode<D: Decoder<Context = C>>(decoder: &mut D) -> Result<Self, DecodeError> {
        let value_range = Decode::decode(decoder)?;
        let params = Vec::<PortableParameter>::decode(decoder)?
            .into_iter()
            .map(PortableParameter::into_atom)
            .collect::<Result<_, _>>()?;
        Ok(Self {
            value_range,
            params,
        })
    }
}
impl ParamValuePairs {
    pub fn validate(&self) {
        assert_eq!(
            self.value_range.len(),
            self.params.len(),
            "Number of values and parameters must match"
        );
    }

    pub fn default_from_symbol(param: Symbol) -> Self {
        Self {
            value_range: 0..1,
            params: vec![Atom::var(param)],
        }
    }
}

impl<'a, A: Into<AtomOrView<'a>>> FromIterator<A> for ParamValuePairs {
    fn from_iter<T: IntoIterator<Item = A>>(iter: T) -> Self {
        let params: Vec<Atom> = iter.into_iter().map(|a| a.into().into_owned()).collect();
        let len = params.len();
        Self {
            value_range: 0..len,
            params,
        }
    }
}

#[cfg(test)]
mod tests {
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
                let CoefficientView::RationalPolynomial(polynomial) = number.get_coeff_view()
                else {
                    panic!("Expected polynomial domain")
                };
                let mut polynomial = polynomial.deserialize();
                polynomial.numerator = MultivariatePolynomial::new(
                    &Z,
                    None,
                    Arc::clone(polynomial.numerator.variables()),
                );
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
                bincode::encode_to_vec(ParamValuePairs::from_iter([value]), config::standard())
                    .unwrap()
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
            let representative = if new_float.as_view().get_data() <= negative.as_view().get_data()
            {
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
}
