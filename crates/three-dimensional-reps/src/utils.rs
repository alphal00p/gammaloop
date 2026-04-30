//! Local helpers built from Symbolica and linnet primitives.
//!
//! Keep this module small and focused. Helpers that turn out to be generally
//! useful should be reported in the top-level implementation log and moved
//! upstream into Symbolica or linnet later.

use std::{
    fmt,
    ops::{Deref, DerefMut},
};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::involution::EdgeIndex;
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::{Q, RationalField},
    parse,
    tensors::matrix::Matrix,
};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Encode, Decode)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct StringSerializedAtom(pub Atom);

impl Deref for StringSerializedAtom {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for StringSerializedAtom {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for StringSerializedAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Serialize for StringSerializedAtom {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.0.to_canonical_string().serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for StringSerializedAtom {
    fn deserialize<D>(deserializer: D) -> Result<StringSerializedAtom, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(StringSerializedAtom(parse!(String::deserialize(
            deserializer
        )?)))
    }
}

pub(crate) mod serde_atom {
    use serde::{Deserialize, Serialize};
    use symbolica::{
        atom::{Atom, AtomCore},
        parse,
    };

    pub(crate) fn serialize<S>(atom: &Atom, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        atom.to_canonical_string().serialize(serializer)
    }

    pub(crate) fn deserialize<'de, D>(deserializer: D) -> Result<Atom, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(parse!(String::deserialize(deserializer)?))
    }
}

pub(crate) mod serde_atom_terms {
    use super::EdgeIndex;
    use serde::{Deserialize, Serialize};
    use symbolica::{
        atom::{Atom, AtomCore},
        parse,
    };

    pub(crate) fn serialize<S>(
        terms: &[(EdgeIndex, Atom)],
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        terms
            .iter()
            .map(|(edge, atom)| (*edge, atom.to_canonical_string()))
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub(crate) fn deserialize<'de, D>(deserializer: D) -> Result<Vec<(EdgeIndex, Atom)>, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(Vec::<(EdgeIndex, String)>::deserialize(deserializer)?
            .into_iter()
            .map(|(edge, atom)| (edge, parse!(atom)))
            .collect())
    }
}

pub(crate) use symbolica::domains::rational::Rational;

pub(crate) trait RationalExt {
    fn pow_usize(&self, exponent: usize) -> Rational;
    fn to_i64_pair(&self) -> Option<(i64, i64)>;
    fn signum_i32(&self) -> i32;
}

impl RationalExt for Rational {
    fn pow_usize(&self, exponent: usize) -> Rational {
        self.pow(u64::try_from(exponent).expect("rational exponent exceeds u64"))
    }

    fn to_i64_pair(&self) -> Option<(i64, i64)> {
        Some((
            self.numerator_ref().to_i64()?,
            self.denominator_ref().to_i64()?,
        ))
    }

    fn signum_i32(&self) -> i32 {
        if self.is_zero() {
            0
        } else if self.is_negative() {
            -1
        } else {
            1
        }
    }
}

pub(crate) fn rational_from_usize(value: usize) -> Rational {
    Rational::from(value)
}

pub(crate) fn rational_pow_i64(base: i64, exponent: usize) -> Rational {
    Rational::from(base).pow_usize(exponent)
}

pub(crate) fn rank_i64(rows: &[Vec<i64>]) -> usize {
    rational_matrix_from_i64(rows).map_or(0, |matrix| matrix.rank())
}

pub(crate) fn rank_rational(rows: &[Vec<Rational>]) -> usize {
    rational_matrix_from_rows(rows).map_or(0, |matrix| matrix.rank())
}

pub(crate) fn solve_rational_system(
    matrix: Vec<Vec<Rational>>,
    rhs: Vec<Rational>,
) -> Option<Vec<Rational>> {
    let n = matrix.len();
    if rhs.len() != n || matrix.iter().any(|row| row.len() != n) {
        return None;
    }

    let lhs = rational_matrix_from_rows(&matrix)?;
    let rhs =
        Matrix::from_nested_vec(rhs.into_iter().map(|value| vec![value]).collect(), Q).ok()?;
    let solution = lhs.solve(&rhs).ok()?;
    Some(
        solution
            .row_iter()
            .map(|row| {
                row.first()
                    .expect("matrix solve output is a vector")
                    .clone()
            })
            .collect(),
    )
}

pub(crate) fn determinant_i32(matrix: &[Vec<i32>]) -> Rational {
    if matrix.is_empty() {
        return Rational::one();
    }
    if matrix.iter().any(|row| row.len() != matrix.len()) {
        return Rational::zero();
    }
    rational_matrix_from_i32(matrix)
        .and_then(|matrix| matrix.det().ok())
        .unwrap_or_else(Rational::zero)
}

pub(crate) fn determinant_i32_is_nonzero(matrix: &[Vec<i32>]) -> bool {
    !determinant_i32(matrix).is_zero()
}

pub(crate) fn determinant_i32_signum(matrix: &[Vec<i32>]) -> i32 {
    determinant_i32(matrix).signum_i32()
}

pub(crate) fn factorial(value: usize) -> Rational {
    (1..=value).fold(Rational::one(), |acc, item| acc * rational_from_usize(item))
}

pub(crate) fn binomial(n: usize, k: usize) -> Rational {
    if k > n {
        return Rational::zero();
    }
    let k = k.min(n - k);
    (0..k).fold(Rational::one(), |acc, item| {
        acc * rational_from_usize(n - item) / rational_from_usize(item + 1)
    })
}

pub(crate) fn rising(value: usize, order: usize) -> Rational {
    (0..order).fold(Rational::one(), |acc, offset| {
        acc * rational_from_usize(value + offset)
    })
}

pub(crate) fn multi_factorial(values: &[usize]) -> Rational {
    values
        .iter()
        .fold(Rational::one(), |acc, value| acc * factorial(*value))
}

pub(crate) fn multiindices_leq(bounds: &[usize]) -> Vec<Vec<usize>> {
    if bounds.is_empty() {
        return vec![Vec::new()];
    }
    let mut out = vec![Vec::new()];
    for bound in bounds {
        let mut next = Vec::new();
        for prefix in &out {
            for value in 0..=*bound {
                let mut item = prefix.clone();
                item.push(value);
                next.push(item);
            }
        }
        out = next;
    }
    out
}

fn rational_matrix_from_i64(rows: &[Vec<i64>]) -> Option<Matrix<RationalField>> {
    rational_matrix_from_rows(
        &rows
            .iter()
            .map(|row| row.iter().map(|value| Rational::from(*value)).collect())
            .collect::<Vec<Vec<_>>>(),
    )
}

fn rational_matrix_from_i32(rows: &[Vec<i32>]) -> Option<Matrix<RationalField>> {
    rational_matrix_from_rows(
        &rows
            .iter()
            .map(|row| row.iter().map(|value| Rational::from(*value)).collect())
            .collect::<Vec<Vec<_>>>(),
    )
}

fn rational_matrix_from_rows(rows: &[Vec<Rational>]) -> Option<Matrix<RationalField>> {
    if rows.is_empty() {
        return None;
    }
    Matrix::from_nested_vec(rows.to_vec(), Q).ok()
}
