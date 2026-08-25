use std::{
    fmt::Display,
    ops::{Add, Neg},
};

use linnet::num_traits::SignOrZero;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A linear momentum combination whose coefficients are restricted to `-1`, `0`, and `1`.
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Signature {
    coefficients: Vec<SignOrZero>,
}

/// Invalid signature data or application.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum SignatureError {
    /// Integer coefficients must be signs or zero.
    #[error("invalid signature coefficient {value}; expected -1, 0, or 1")]
    InvalidCoefficient {
        /// Invalid integer value.
        value: i8,
    },
    /// A signature must have one coefficient per basis vector.
    #[error("signature expects {expected} basis vectors, received {actual}")]
    BasisLength {
        /// Number of signature coefficients.
        expected: usize,
        /// Number of supplied basis vectors.
        actual: usize,
    },
    /// Adding two sign-only signatures would produce a coefficient outside the domain.
    #[error("signature coefficient at index {index} would be outside -1, 0, or 1")]
    CoefficientOverflow {
        /// Coefficient position.
        index: usize,
    },
}

impl Signature {
    /// Construct a signature from validated sign coefficients.
    pub fn new(coefficients: impl IntoIterator<Item = SignOrZero>) -> Self {
        Self {
            coefficients: coefficients.into_iter().collect(),
        }
    }

    /// Validate integer coefficients and construct a signature.
    pub fn try_from_integers(
        coefficients: impl IntoIterator<Item = i8>,
    ) -> Result<Self, SignatureError> {
        coefficients
            .into_iter()
            .map(|value| {
                SignOrZero::try_from(value)
                    .map_err(|_| SignatureError::InvalidCoefficient { value })
            })
            .collect()
    }

    /// Number of coefficients.
    pub fn len(&self) -> usize {
        self.coefficients.len()
    }

    /// Whether there are no coefficients.
    pub fn is_empty(&self) -> bool {
        self.coefficients.is_empty()
    }

    /// Borrow a coefficient by position.
    pub fn get(&self, index: usize) -> Option<SignOrZero> {
        self.coefficients.get(index).copied()
    }

    /// Iterate over the coefficients.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = SignOrZero> + '_ {
        self.coefficients.iter().copied()
    }

    /// Export coefficients in the format used by graph and loop-momentum tools.
    pub fn integer_coefficients(&self) -> Vec<isize> {
        self.iter().map(sign_value).collect()
    }

    /// Return an overall-sign canonical form whose first nonzero coefficient is positive.
    pub fn canonicalized(&self) -> Self {
        match self.iter().find(|coefficient| coefficient.is_sign()) {
            Some(SignOrZero::Minus) => Self::new(self.iter().map(Neg::neg)),
            _ => self.clone(),
        }
    }

    /// Add two sign-only signatures, rejecting coefficients with magnitude greater than one.
    pub fn checked_add(&self, rhs: &Self) -> Result<Self, SignatureError> {
        if self.len() != rhs.len() {
            return Err(SignatureError::BasisLength {
                expected: self.len(),
                actual: rhs.len(),
            });
        }

        self.iter()
            .zip(rhs.iter())
            .enumerate()
            .map(|(index, (left, right))| match (left, right) {
                (SignOrZero::Zero, sign) | (sign, SignOrZero::Zero) => Ok(sign),
                (SignOrZero::Plus, SignOrZero::Minus) | (SignOrZero::Minus, SignOrZero::Plus) => {
                    Ok(SignOrZero::Zero)
                }
                _ => Err(SignatureError::CoefficientOverflow { index }),
            })
            .collect()
    }

    /// Apply the signature to a basis.
    ///
    /// An all-zero signature returns `None`, avoiding a separate zero-construction
    /// trait for arbitrary basis values.
    pub fn apply<B>(&self, basis: &[B]) -> Result<Option<B>, SignatureError>
    where
        B: Clone + Neg<Output = B> + Add<Output = B>,
    {
        if self.len() != basis.len() {
            return Err(SignatureError::BasisLength {
                expected: self.len(),
                actual: basis.len(),
            });
        }

        Ok(self
            .iter()
            .zip(basis.iter())
            .filter_map(|(coefficient, value)| match coefficient {
                SignOrZero::Zero => None,
                SignOrZero::Plus => Some(value.clone()),
                SignOrZero::Minus => Some(-value.clone()),
            })
            .reduce(Add::add))
    }
}

impl FromIterator<SignOrZero> for Signature {
    fn from_iter<I: IntoIterator<Item = SignOrZero>>(iter: I) -> Self {
        Self::new(iter)
    }
}

impl IntoIterator for Signature {
    type Item = SignOrZero;
    type IntoIter = std::vec::IntoIter<SignOrZero>;

    fn into_iter(self) -> Self::IntoIter {
        self.coefficients.into_iter()
    }
}

impl<'a> IntoIterator for &'a Signature {
    type Item = SignOrZero;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, SignOrZero>>;

    fn into_iter(self) -> Self::IntoIter {
        self.coefficients.iter().copied()
    }
}

impl Display for Signature {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for coefficient in &self.coefficients {
            write!(formatter, "{coefficient}")?;
        }
        Ok(())
    }
}

/// Coefficients of one momentum in loop and external momentum bases.
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct MomentumSignature {
    /// Loop-momentum coefficients.
    pub loops: Signature,
    /// External-momentum coefficients.
    pub external: Signature,
}

impl MomentumSignature {
    /// Construct a combined momentum signature.
    pub const fn new(loops: Signature, external: Signature) -> Self {
        Self { loops, external }
    }

    /// Apply both parts to their corresponding bases.
    pub fn apply<B>(
        &self,
        loop_momenta: &[B],
        external_momenta: &[B],
    ) -> Result<Option<B>, SignatureError>
    where
        B: Clone + Neg<Output = B> + Add<Output = B>,
    {
        let loop_part = self.loops.apply(loop_momenta)?;
        let external_part = self.external.apply(external_momenta)?;
        Ok(match (loop_part, external_part) {
            (Some(loop_part), Some(external_part)) => Some(loop_part + external_part),
            (Some(loop_part), None) => Some(loop_part),
            (None, Some(external_part)) => Some(external_part),
            (None, None) => None,
        })
    }

    /// Export `(loop, external)` integer coefficients.
    pub fn integer_coefficients(&self) -> (Vec<isize>, Vec<isize>) {
        (
            self.loops.integer_coefficients(),
            self.external.integer_coefficients(),
        )
    }

    /// Return a canonical form up to one overall sign across both bases.
    pub fn canonicalized(&self) -> Self {
        let first = self
            .loops
            .iter()
            .chain(self.external.iter())
            .find(|coefficient| coefficient.is_sign());
        if first == Some(SignOrZero::Minus) {
            Self {
                loops: Signature::new(self.loops.iter().map(Neg::neg)),
                external: Signature::new(self.external.iter().map(Neg::neg)),
            }
        } else {
            self.clone()
        }
    }

    /// Format the linear combination using `k_i` and `p_i` basis labels.
    pub fn format_momentum(&self) -> String {
        let terms = self
            .loops
            .iter()
            .enumerate()
            .map(|(index, sign)| (sign, format!("k_{index}")))
            .chain(
                self.external
                    .iter()
                    .enumerate()
                    .map(|(index, sign)| (sign, format!("p_{index}"))),
            );

        let mut result = String::new();
        for (sign, label) in terms.filter(|(sign, _)| sign.is_sign()) {
            if result.is_empty() {
                if sign == SignOrZero::Minus {
                    result.push('-');
                }
            } else if sign == SignOrZero::Plus {
                result.push('+');
            } else {
                result.push('-');
            }
            result.push_str(&label);
        }
        if result.is_empty() {
            result.push('0');
        }
        result
    }
}

fn sign_value(sign: SignOrZero) -> isize {
    match sign {
        SignOrZero::Zero => 0,
        SignOrZero::Plus => 1,
        SignOrZero::Minus => -1,
    }
}

#[cfg(test)]
mod tests {
    use crate::FourMomentum;

    use super::*;

    #[test]
    fn applies_to_scalar_and_momentum_bases() {
        let signature = Signature::new([
            SignOrZero::Plus,
            SignOrZero::Minus,
            SignOrZero::Zero,
            SignOrZero::Plus,
        ]);
        assert_eq!(signature.apply(&[1, 2, 3, 4]).unwrap(), Some(3));

        let basis = [
            FourMomentum::from_args(1, 1, 0, 0),
            FourMomentum::from_args(1, 0, 1, 0),
            FourMomentum::from_args(1, 0, 0, 1),
            FourMomentum::from_args(1, 1, 1, 1),
        ];
        assert_eq!(
            signature.apply(&basis).unwrap(),
            Some(FourMomentum::from_args(1, 2, 0, 1))
        );
    }

    #[test]
    fn all_zero_signature_has_no_value() {
        let signature = Signature::new([SignOrZero::Zero; 3]);
        assert_eq!(signature.apply(&[1, 2, 3]).unwrap(), None);
    }

    #[test]
    fn checked_add_rejects_non_sign_coefficients() {
        let plus = Signature::new([SignOrZero::Plus]);
        assert_eq!(
            plus.checked_add(&plus),
            Err(SignatureError::CoefficientOverflow { index: 0 })
        );
    }

    #[test]
    fn combined_canonicalization_uses_one_overall_sign() {
        let signature = MomentumSignature::new(
            Signature::new([SignOrZero::Zero, SignOrZero::Minus]),
            Signature::new([SignOrZero::Plus]),
        );
        let canonical = signature.canonicalized();
        assert_eq!(canonical.format_momentum(), "k_1-p_0");
    }
}
