use std::{fmt::Display, str::FromStr};

use linnet::num_traits::{Sign, SignOrZero};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A fixed external-state helicity in GammaLoop's integer normalization.
///
/// Fermion helicities use `-1` and `+1` rather than half-integers. Summation
/// policies belong to event or runtime configuration and are deliberately not
/// represented by this fixed-helicity primitive.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Helicity(SignOrZero);

/// An integer or textual value is not a fixed helicity.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum HelicityError {
    #[error("invalid helicity {value}; expected -1, 0, or 1")]
    InvalidInteger { value: i8 },
    #[error("invalid helicity '{value}'; expected minus, zero, or plus")]
    InvalidText { value: String },
}

impl Helicity {
    pub const MINUS: Self = Self(SignOrZero::Minus);
    pub const ZERO: Self = Self(SignOrZero::Zero);
    pub const PLUS: Self = Self(SignOrZero::Plus);

    /// Construct a fixed helicity from a validated sign or zero.
    pub const fn new(value: SignOrZero) -> Self {
        Self(value)
    }

    /// Return the underlying sign representation.
    pub const fn sign(self) -> SignOrZero {
        self.0
    }

    /// Return the integer convention used by serialized kinematics.
    pub const fn integer(self) -> i8 {
        match self.0 {
            SignOrZero::Minus => -1,
            SignOrZero::Zero => 0,
            SignOrZero::Plus => 1,
        }
    }

    pub const fn is_minus(self) -> bool {
        matches!(self.0, SignOrZero::Minus)
    }

    pub const fn is_zero(self) -> bool {
        matches!(self.0, SignOrZero::Zero)
    }

    pub const fn is_plus(self) -> bool {
        matches!(self.0, SignOrZero::Plus)
    }
}

impl TryFrom<i8> for Helicity {
    type Error = HelicityError;

    fn try_from(value: i8) -> Result<Self, Self::Error> {
        SignOrZero::try_from(value)
            .map(Self)
            .map_err(|_| HelicityError::InvalidInteger { value })
    }
}

impl From<SignOrZero> for Helicity {
    fn from(value: SignOrZero) -> Self {
        Self(value)
    }
}

impl From<Sign> for Helicity {
    fn from(value: Sign) -> Self {
        Self(value.into())
    }
}

impl From<Helicity> for SignOrZero {
    fn from(value: Helicity) -> Self {
        value.0
    }
}

impl Display for Helicity {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(match self.0 {
            SignOrZero::Minus => "-",
            SignOrZero::Zero => "0",
            SignOrZero::Plus => "+",
        })
    }
}

impl FromStr for Helicity {
    type Err = HelicityError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "-1" | "-" | "minus" => Ok(Self::MINUS),
            "0" | "zero" => Ok(Self::ZERO),
            "1" | "+1" | "+" | "plus" => Ok(Self::PLUS),
            _ => Err(HelicityError::InvalidText {
                value: value.to_owned(),
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validates_integer_helicities() {
        assert_eq!(Helicity::try_from(-1).unwrap(), Helicity::MINUS);
        assert_eq!(Helicity::try_from(0).unwrap(), Helicity::ZERO);
        assert_eq!(Helicity::try_from(1).unwrap(), Helicity::PLUS);
        assert_eq!(
            Helicity::try_from(2).unwrap_err(),
            HelicityError::InvalidInteger { value: 2 }
        );
    }

    #[test]
    fn serde_uses_the_integer_convention() {
        assert_eq!(serde_json::to_string(&Helicity::MINUS).unwrap(), "-1");
        assert_eq!(
            serde_json::from_str::<Helicity>("0").unwrap(),
            Helicity::ZERO
        );
        assert!(serde_json::from_str::<Helicity>("2").is_err());
    }

    #[test]
    fn parses_and_displays_named_helicities() {
        assert_eq!("minus".parse(), Ok(Helicity::MINUS));
        assert_eq!("+1".parse(), Ok(Helicity::PLUS));
        assert_eq!(Helicity::ZERO.to_string(), "0");
    }
}
