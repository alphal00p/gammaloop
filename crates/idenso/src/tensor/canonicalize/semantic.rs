use spenso::structure::{
    representation::{LibraryRep, LibrarySlot, Representation},
    slot::{AbsInd, IsAbstractSlot, ParseableAind},
};
use symbolica::{
    atom::{Atom, AtomView, Symbol, UserData, UserDataKey},
    coefficient::CoefficientView,
    domains::integer::IntegerRing,
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
};

/// Process-independent identity for a Symbolica symbol used in graph colors.
///
/// Symbolica's numeric symbol IDs depend on interning order. Canonical graph
/// colors instead use the qualified name and the semantic flags that can alter
/// normalization or tensor interpretation. Symbol user data is retained as a
/// complete recursive payload rather than collapsed to process-local identity.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct SemanticSymbolKey {
    name: String,
    wildcard_level: u8,
    attributes: [bool; 8],
    tags: Vec<String>,
    payload: SemanticUserDataKey,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SemanticUserDataMapKey {
    Integer(i64),
    String(String),
    Atom(Box<SemanticAtomKey>),
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SemanticUserDataKey {
    None,
    Atom(Box<SemanticAtomKey>),
    Integer(i64),
    String(String),
    List(Vec<Self>),
    Map(Vec<(SemanticUserDataMapKey, Self)>),
    Serialized(Vec<u8>),
    Unknown(String),
}

impl SemanticUserDataMapKey {
    fn new(key: &UserDataKey) -> Self {
        match key {
            UserDataKey::Integer(value) => Self::Integer(*value),
            UserDataKey::String(value) => Self::String(value.clone()),
            UserDataKey::Atom(value) => Self::Atom(Box::new(SemanticAtomKey::new(value.as_view()))),
        }
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        match self {
            Self::Integer(value) => {
                hasher.update(&[0]);
                hasher.update(&value.to_le_bytes());
            }
            Self::String(value) => {
                hasher.update(&[1]);
                update_digest_bytes(hasher, value.as_bytes());
            }
            Self::Atom(value) => {
                hasher.update(&[2]);
                value.update_digest(hasher);
            }
        }
    }
}

impl SemanticUserDataKey {
    fn new(data: &UserData) -> Self {
        match data {
            UserData::None => Self::None,
            UserData::Atom(value) => Self::Atom(Box::new(SemanticAtomKey::new(value.as_view()))),
            UserData::Integer(value) => Self::Integer(*value),
            UserData::String(value) => Self::String(value.clone()),
            UserData::List(values) => Self::List(values.iter().map(Self::new).collect()),
            UserData::Map(values) => {
                let mut values = values
                    .iter()
                    .map(|(key, value)| (SemanticUserDataMapKey::new(key), Self::new(value)))
                    .collect::<Vec<_>>();
                values.sort();
                Self::Map(values)
            }
            UserData::Serialized(value) => Self::Serialized(value.clone()),
            value => Self::Unknown(format!("{value:?}")),
        }
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        match self {
            Self::None => {
                hasher.update(&[0]);
            }
            Self::Atom(value) => {
                hasher.update(&[1]);
                value.update_digest(hasher);
            }
            Self::Integer(value) => {
                hasher.update(&[2]);
                hasher.update(&value.to_le_bytes());
            }
            Self::String(value) => {
                hasher.update(&[3]);
                update_digest_bytes(hasher, value.as_bytes());
            }
            Self::List(values) => {
                hasher.update(&[4]);
                update_digest_len(hasher, values.len());
                for value in values {
                    value.update_digest(hasher);
                }
            }
            Self::Map(values) => {
                hasher.update(&[5]);
                update_digest_len(hasher, values.len());
                for (key, value) in values {
                    key.update_digest(hasher);
                    value.update_digest(hasher);
                }
            }
            Self::Serialized(value) => {
                hasher.update(&[6]);
                update_digest_bytes(hasher, value);
            }
            Self::Unknown(value) => {
                hasher.update(&[u8::MAX]);
                update_digest_bytes(hasher, value.as_bytes());
            }
        }
    }
}

impl From<Symbol> for SemanticSymbolKey {
    fn from(symbol: Symbol) -> Self {
        let mut tags = symbol.get_tags().to_vec();
        tags.sort();
        Self {
            name: symbol.get_name().to_owned(),
            wildcard_level: symbol.get_wildcard_level(),
            attributes: [
                symbol.is_symmetric(),
                symbol.is_antisymmetric(),
                symbol.is_cyclesymmetric(),
                symbol.is_linear(),
                symbol.is_scalar(),
                symbol.is_real(),
                symbol.is_integer(),
                symbol.is_positive(),
            ],
            tags,
            payload: SemanticUserDataKey::new(symbol.get_data()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SemanticPolyVariableKey {
    Symbol(SemanticSymbolKey),
    Function(SemanticSymbolKey, Box<SemanticAtomKey>),
    Power(Box<SemanticAtomKey>),
    Temporary(usize),
}

impl SemanticPolyVariableKey {
    fn new(variable: &PolyVariable) -> Self {
        match variable {
            PolyVariable::Symbol(symbol) => Self::Symbol((*symbol).into()),
            PolyVariable::Function(symbol, atom) => Self::Function(
                (*symbol).into(),
                Box::new(SemanticAtomKey::new(atom.as_view())),
            ),
            PolyVariable::Power(atom) => {
                Self::Power(Box::new(SemanticAtomKey::new(atom.as_view())))
            }
            PolyVariable::Temporary(temporary) => Self::Temporary(*temporary),
        }
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        match self {
            Self::Symbol(symbol) => {
                hasher.update(&[0]);
                symbol.update_digest(hasher);
            }
            Self::Function(symbol, atom) => {
                hasher.update(&[1]);
                symbol.update_digest(hasher);
                atom.update_digest(hasher);
            }
            Self::Power(atom) => {
                hasher.update(&[2]);
                atom.update_digest(hasher);
            }
            Self::Temporary(temporary) => {
                hasher.update(&[3]);
                update_digest_len(hasher, *temporary);
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct SemanticPolynomialTermKey {
    powers: Vec<(SemanticPolyVariableKey, u16)>,
    coefficient: String,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct SemanticPolynomialKey(Vec<SemanticPolynomialTermKey>);

impl SemanticPolynomialKey {
    fn new(polynomial: &MultivariatePolynomial<IntegerRing, u16>) -> Self {
        let variables = polynomial
            .variables
            .iter()
            .map(SemanticPolyVariableKey::new)
            .collect::<Vec<_>>();
        let width = variables.len();
        let mut terms = polynomial
            .coefficients
            .iter()
            .enumerate()
            .map(|(term, coefficient)| {
                let powers = variables
                    .iter()
                    .cloned()
                    .zip(&polynomial.exponents[term * width..(term + 1) * width])
                    .filter_map(|(variable, exponent)| {
                        (*exponent != 0).then_some((variable, *exponent))
                    })
                    .collect::<Vec<_>>();
                let mut powers = powers;
                powers.sort();
                SemanticPolynomialTermKey {
                    powers,
                    coefficient: coefficient.to_string(),
                }
            })
            .collect::<Vec<_>>();
        terms.sort();
        Self(terms)
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        update_digest_len(hasher, self.0.len());
        for term in &self.0 {
            update_digest_len(hasher, term.powers.len());
            for (variable, exponent) in &term.powers {
                variable.update_digest(hasher);
                hasher.update(&exponent.to_le_bytes());
            }
            update_digest_bytes(hasher, term.coefficient.as_bytes());
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct SemanticRationalPolynomialKey {
    numerator: SemanticPolynomialKey,
    denominator: SemanticPolynomialKey,
}

impl SemanticRationalPolynomialKey {
    fn new(
        polynomial: symbolica::domains::rational_polynomial::RationalPolynomial<IntegerRing>,
    ) -> Self {
        Self {
            numerator: SemanticPolynomialKey::new(&polynomial.numerator),
            denominator: SemanticPolynomialKey::new(&polynomial.denominator),
        }
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        self.numerator.update_digest(hasher);
        self.denominator.update_digest(hasher);
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum SemanticNumberKey {
    Indeterminate,
    Infinity(Option<(String, String)>),
    Rational { real: String, imaginary: String },
    Float { real: Vec<u8>, imaginary: Vec<u8> },
    FiniteField(String),
    RationalPolynomial(SemanticRationalPolynomialKey),
}

impl SemanticNumberKey {
    fn new(number: CoefficientView<'_>) -> Self {
        match number {
            CoefficientView::Indeterminate => Self::Indeterminate,
            CoefficientView::Infinity(direction) => Self::Infinity(
                direction
                    .map(|(real, imag)| (real.to_rat().to_string(), imag.to_rat().to_string())),
            ),
            CoefficientView::Natural(real, real_den, imaginary, imaginary_den) => Self::Rational {
                real: symbolica::domains::rational::Rational::from_int_unchecked(real, real_den)
                    .to_string(),
                imaginary: symbolica::domains::rational::Rational::from_int_unchecked(
                    imaginary,
                    imaginary_den,
                )
                .to_string(),
            },
            CoefficientView::Large(real, imaginary) => Self::Rational {
                real: real.to_rat().to_string(),
                imaginary: imaginary.to_rat().to_string(),
            },
            CoefficientView::Float(real, imaginary) => Self::Float {
                real: real.to_float().serialize(),
                imaginary: imaginary.to_float().serialize(),
            },
            CoefficientView::FiniteField(_, _) => Self::FiniteField(number.to_owned().to_string()),
            CoefficientView::RationalPolynomial(polynomial) => Self::RationalPolynomial(
                SemanticRationalPolynomialKey::new(polynomial.deserialize()),
            ),
        }
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        match self {
            Self::Indeterminate => {
                hasher.update(&[0]);
            }
            Self::Infinity(direction) => {
                hasher.update(&[1]);
                if let Some((real, imaginary)) = direction {
                    hasher.update(&[1]);
                    update_digest_bytes(hasher, real.as_bytes());
                    update_digest_bytes(hasher, imaginary.as_bytes());
                } else {
                    hasher.update(&[0]);
                }
            }
            Self::Rational { real, imaginary } => {
                hasher.update(&[2]);
                update_digest_bytes(hasher, real.as_bytes());
                update_digest_bytes(hasher, imaginary.as_bytes());
            }
            Self::Float { real, imaginary } => {
                hasher.update(&[3]);
                update_digest_bytes(hasher, real);
                update_digest_bytes(hasher, imaginary);
            }
            Self::FiniteField(value) => {
                hasher.update(&[4]);
                update_digest_bytes(hasher, value.as_bytes());
            }
            Self::RationalPolynomial(value) => {
                hasher.update(&[5]);
                value.update_digest(hasher);
            }
        }
    }
}

/// Stable recursive Atom identity shared by graph colors and proof keys.
///
/// Sum and Product children are sorted by this key instead of Symbolica's
/// interner-sensitive internal term order. Function arguments remain ordered.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum SemanticAtomKey {
    Number(SemanticNumberKey),
    Variable(SemanticSymbolKey),
    Function {
        head: SemanticSymbolKey,
        arguments: Vec<Self>,
    },
    Power {
        base: Box<Self>,
        exponent: Box<Self>,
    },
    Product(Vec<Self>),
    Sum(Vec<Self>),
}

impl SemanticAtomKey {
    pub(super) fn new(value: AtomView<'_>) -> Self {
        match value {
            AtomView::Num(number) => Self::Number(SemanticNumberKey::new(number.get_coeff_view())),
            AtomView::Var(variable) => Self::Variable(variable.get_symbol().into()),
            AtomView::Fun(function) => Self::Function {
                head: function.get_symbol().into(),
                arguments: function.iter().map(Self::new).collect(),
            },
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::Power {
                    base: Box::new(Self::new(base)),
                    exponent: Box::new(Self::new(exponent)),
                }
            }
            AtomView::Mul(product) => {
                let mut factors = product.iter().map(Self::new).collect::<Vec<_>>();
                factors.sort();
                Self::Product(factors)
            }
            AtomView::Add(sum) => {
                let mut terms = sum.iter().map(Self::new).collect::<Vec<_>>();
                terms.sort();
                Self::Sum(terms)
            }
        }
    }

    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        match self {
            Self::Number(number) => {
                hasher.update(&[0]);
                number.update_digest(hasher);
            }
            Self::Variable(symbol) => {
                hasher.update(&[1]);
                symbol.update_digest(hasher);
            }
            Self::Function { head, arguments } => {
                hasher.update(&[2]);
                head.update_digest(hasher);
                update_digest_len(hasher, arguments.len());
                for argument in arguments {
                    argument.update_digest(hasher);
                }
            }
            Self::Power { base, exponent } => {
                hasher.update(&[3]);
                base.update_digest(hasher);
                exponent.update_digest(hasher);
            }
            Self::Product(factors) => {
                hasher.update(&[4]);
                update_digest_len(hasher, factors.len());
                for factor in factors {
                    factor.update_digest(hasher);
                }
            }
            Self::Sum(terms) => {
                hasher.update(&[5]);
                update_digest_len(hasher, terms.len());
                for term in terms {
                    term.update_digest(hasher);
                }
            }
        }
    }
}

impl SemanticSymbolKey {
    fn update_digest(&self, hasher: &mut blake3::Hasher) {
        update_digest_bytes(hasher, self.name.as_bytes());
        hasher.update(&[self.wildcard_level]);
        for attribute in self.attributes {
            hasher.update(&[u8::from(attribute)]);
        }
        update_digest_len(hasher, self.tags.len());
        for tag in &self.tags {
            update_digest_bytes(hasher, tag.as_bytes());
        }
        self.payload.update_digest(hasher);
    }
}

fn update_digest_len(hasher: &mut blake3::Hasher, len: usize) {
    hasher.update(&(len as u128).to_le_bytes());
}

fn update_digest_bytes(hasher: &mut blake3::Hasher, bytes: &[u8]) {
    update_digest_len(hasher, bytes.len());
    hasher.update(bytes);
}

/// Stable digest of the complete semantic Atom key.
///
/// This digest is suitable for generated symbol names, but equality and graph
/// colors continue to retain and compare the complete recursive key.
pub(crate) fn semantic_atom_digest(value: AtomView<'_>) -> blake3::Hash {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"idenso::semantic-atom-key::v1");
    SemanticAtomKey::new(value).update_digest(&mut hasher);
    hasher.finalize()
}

pub(super) fn representation_key(representation: Representation<LibraryRep>) -> SemanticAtomKey {
    SemanticAtomKey::new(representation.to_symbolic([Atom::Zero]).as_view())
}

pub(super) fn slot_key<Aind: AbsInd + ParseableAind>(slot: LibrarySlot<Aind>) -> SemanticAtomKey
where
    LibrarySlot<Aind>: IsAbstractSlot<Aind = Aind>,
{
    SemanticAtomKey::new(slot.to_atom().as_view())
}

#[cfg(test)]
mod tests {
    use std::hash::{Hash, Hasher};

    use symbolica::{
        atom::{Atom, AtomCore, NamespacedSymbol, SymbolBuilder, UserData},
        coefficient::Coefficient,
        domains::float::Float,
        domains::integer::Z,
        parse,
    };

    use super::SemanticAtomKey;

    #[test]
    fn commutative_children_use_the_shared_stable_order() {
        let left = parse!("f(z)+x*y+y*x", default_namespace = "semantic_key");
        let right = parse!("y*x+f(z)+x*y", default_namespace = "semantic_key");

        assert_eq!(
            SemanticAtomKey::new(left.as_view()),
            SemanticAtomKey::new(right.as_view())
        );
    }

    #[test]
    fn ordered_function_arguments_remain_distinct() {
        let left = parse!("f(x,y)", default_namespace = "semantic_key");
        let right = parse!("f(y,x)", default_namespace = "semantic_key");

        assert_ne!(
            SemanticAtomKey::new(left.as_view()),
            SemanticAtomKey::new(right.as_view())
        );
    }

    #[test]
    fn qualified_names_are_semantic_key_data() {
        let first = parse!("f(x)", default_namespace = "semantic_key_first");
        let second = parse!("f(x)", default_namespace = "semantic_key_second");

        assert_ne!(
            SemanticAtomKey::new(first.as_view()),
            SemanticAtomKey::new(second.as_view())
        );
    }

    #[test]
    fn multiprecision_float_identity_does_not_use_lossy_display() {
        let first = Atom::num(Float::with_val(8, 1.0));
        let second = Atom::num(Float::with_val(8, 1.0078125));

        assert_ne!(first, second);
        assert_eq!(first.to_string(), second.to_string());
        assert_ne!(
            SemanticAtomKey::new(first.as_view()),
            SemanticAtomKey::new(second.as_view())
        );
        assert_ne!(
            super::semantic_atom_digest(first.as_view()),
            super::semantic_atom_digest(second.as_view())
        );
    }

    #[test]
    fn rational_polynomial_identity_retains_qualified_variable_names() {
        let first_polynomial = parse!("a::x+1", default_namespace = "semantic_rat_poly")
            .to_rational_polynomial::<_, _, u16>(&Z, &Z, None);
        let second_polynomial = parse!("b::x+1", default_namespace = "semantic_rat_poly")
            .to_rational_polynomial::<_, _, u16>(&Z, &Z, None);
        let first = Atom::num(Coefficient::RationalPolynomial(first_polynomial));
        let second = Atom::num(Coefficient::RationalPolynomial(second_polynomial));

        assert_ne!(first, second);
        assert_eq!(first.to_string(), second.to_string());
        assert_ne!(
            SemanticAtomKey::new(first.as_view()),
            SemanticAtomKey::new(second.as_view())
        );
        assert_ne!(
            super::semantic_atom_digest(first.as_view()),
            super::semantic_atom_digest(second.as_view())
        );
    }

    #[test]
    fn hash_collision_never_substitutes_for_full_key_equality() {
        #[derive(Default)]
        struct CollidingHasher;

        impl Hasher for CollidingHasher {
            fn finish(&self) -> u64 {
                0
            }

            fn write(&mut self, _bytes: &[u8]) {}
        }

        let left = SemanticAtomKey::new(
            parse!("f(x)", default_namespace = "semantic_key_collision").as_view(),
        );
        let right = SemanticAtomKey::new(
            parse!("g(x)", default_namespace = "semantic_key_collision").as_view(),
        );
        let hash = |key: &SemanticAtomKey| {
            let mut hasher = CollidingHasher;
            key.hash(&mut hasher);
            hasher.finish()
        };

        assert_eq!(hash(&left), hash(&right));
        assert_ne!(left, right);
    }

    #[test]
    fn symbol_keys_retain_the_complete_cooked_atom_payload() {
        let payload = parse!(
            "outer(ns::inner(x+y),z^2)",
            default_namespace = "semantic_key_payload"
        );
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(
            "semantic_key_payload::cooked_payload",
        ))
        .with_user_data(UserData::Atom(payload.clone()))
        .build()
        .unwrap();
        let SemanticAtomKey::Variable(key) = SemanticAtomKey::new(Atom::var(symbol).as_view())
        else {
            panic!("expected a variable semantic key");
        };

        assert_eq!(
            key.payload,
            super::SemanticUserDataKey::Atom(Box::new(SemanticAtomKey::new(payload.as_view())))
        );
    }

    #[test]
    fn symbol_keys_retain_recursive_user_data_with_stable_map_order() {
        let payload = parse!("f(x+y)", default_namespace = "semantic_user_data");
        let first = UserData::Map(
            [
                (
                    symbolica::atom::UserDataKey::String("nested".to_owned()),
                    UserData::List(vec![UserData::Integer(7), UserData::Atom(payload.clone())]),
                ),
                (
                    symbolica::atom::UserDataKey::Integer(3),
                    UserData::Serialized(vec![0, 1, 255]),
                ),
            ]
            .into_iter()
            .collect(),
        );
        let second = UserData::Map(
            [
                (
                    symbolica::atom::UserDataKey::Integer(3),
                    UserData::Serialized(vec![0, 1, 255]),
                ),
                (
                    symbolica::atom::UserDataKey::String("nested".to_owned()),
                    UserData::List(vec![UserData::Integer(7), UserData::Atom(payload)]),
                ),
            ]
            .into_iter()
            .collect(),
        );

        assert_eq!(
            super::SemanticUserDataKey::new(&first),
            super::SemanticUserDataKey::new(&second)
        );
    }

    #[test]
    fn semantic_digest_uses_commutative_key_order() {
        let left = parse!("f(x+y,x*y)", default_namespace = "semantic_digest");
        let right = parse!("f(y+x,y*x)", default_namespace = "semantic_digest");

        assert_eq!(
            super::semantic_atom_digest(left.as_view()),
            super::semantic_atom_digest(right.as_view())
        );
    }
}
