//! Exact quantum numbers at the JSON boundary. Fractional values must be strings;
//! accepting a rounded float such as 0.6666666666666666 would change the charge.

use serde::{Deserialize, Deserializer, Serialize, Serializer, de};
use symbolica::domains::{integer::Integer, rational::Rational};

#[derive(Deserialize)]
#[serde(untagged)]
enum Input {
    Integer(i64),
    Unsigned(u64),
    Fraction(String),
    Float(f64),
}

impl Input {
    fn rational<E: de::Error>(self) -> Result<Rational, E> {
        match self {
            Self::Integer(n) => Ok(n.into()),
            Self::Unsigned(n) => Ok(n.into()),
            Self::Fraction(s) => {
                let (n, d) = s.split_once('/').unwrap_or((&s, "1"));
                let n: Integer = n.trim().parse().map_err(E::custom)?;
                let d: Integer = d.trim().parse().map_err(E::custom)?;
                if d.is_zero() {
                    return Err(E::custom("quantum number denominator cannot be zero"));
                }
                Ok(Rational::from((n, d)))
            }
            // UFO JSON also uses floating-point literals for integral charges.
            Self::Float(n) if n.fract() == 0.0 && n.abs() <= 9_007_199_254_740_991.0 => {
                Rational::try_from(n).map_err(E::custom)
            }
            Self::Float(_) => Err(E::custom(
                "fractional or large floating-point quantum numbers require an exact string, e.g. \"2/3\"",
            )),
        }
    }
}

pub(super) fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Rational, D::Error> {
    Input::deserialize(d)?.rational()
}

pub(super) fn serialize<S: Serializer>(value: &Rational, s: S) -> Result<S::Ok, S::Error> {
    s.collect_str(value)
}

pub(super) mod optional {
    use super::*;

    pub(crate) fn deserialize<'de, D: Deserializer<'de>>(
        d: D,
    ) -> Result<Option<Rational>, D::Error> {
        Option::<Input>::deserialize(d)?
            .map(Input::rational)
            .transpose()
    }

    pub(crate) fn serialize<S: Serializer>(
        value: &Option<Rational>,
        s: S,
    ) -> Result<S::Ok, S::Error> {
        value.as_ref().map(ToString::to_string).serialize(s)
    }
}
