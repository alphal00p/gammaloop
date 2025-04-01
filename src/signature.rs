use crate::momentum::{FourMomentum, SignOrZero, ThreeMomentum};
use crate::momentum_sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta,
};
use crate::utils::{FloatLike, Length, F};
use bincode::{BorrowDecode, Decode, Encode};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use spenso::contraction::RefZero;
use std::fmt::Display;
use std::ops::{AddAssign, Index, IndexMut, Neg, SubAssign};
use typed_index_collections::TiVec;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, PartialOrd, Eq, Ord, Hash)]
pub struct SignatureLike<T: From<usize>>(TiVec<T, SignOrZero>);
pub type LoopSignature = SignatureLike<LoopIndex>;
pub type ExternalSignature = SignatureLike<ExternalIndex>;

// manual implementations because TiVec is not Encode/Devode
impl<T: Encode + From<usize>> Encode for SignatureLike<T> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        self.0.raw.encode(encoder)
    }
}

impl<'de, Context, T: Decode<Context> + From<usize>> BorrowDecode<'de, Context>
    for SignatureLike<T>
{
    fn borrow_decode<D: bincode::de::Decoder>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(SignatureLike(Vec::decode(decoder)?.into()))
    }
}

impl<Context, T: Decode<Context> + From<usize>> Decode<Context> for SignatureLike<T> {
    fn decode<D: bincode::de::Decoder>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(SignatureLike(Vec::decode(decoder)?.into()))
    }
}

#[derive(
    Debug, Clone, Serialize, Deserialize, PartialEq, PartialOrd, Eq, Ord, Hash, Encode, Decode,
)]
pub struct LoopExtSignature {
    pub internal: LoopSignature,
    pub external: ExternalSignature,
}

impl From<(Vec<isize>, Vec<isize>)> for LoopExtSignature {
    fn from(value: (Vec<isize>, Vec<isize>)) -> Self {
        Self {
            internal: LoopSignature::from_iter(value.0),
            external: ExternalSignature::from_iter(value.1),
        }
    }
}

impl<T> Index<T> for SignatureLike<T>
where
    usize: From<T>,
    T: From<usize>,
{
    type Output = SignOrZero;

    fn index(&self, index: T) -> &Self::Output {
        &self.0[index]
    }
}

impl<T> IndexMut<T> for SignatureLike<T>
where
    usize: From<T>,
    T: From<usize>,
{
    fn index_mut(&mut self, index: T) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T> Display for SignatureLike<T>
where
    T: From<usize>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for sign in &self.0 {
            write!(f, "{}", sign)?;
        }
        Ok(())
    }
}

impl<T> FromIterator<SignOrZero> for SignatureLike<T>
where
    T: From<usize>,
{
    fn from_iter<I: IntoIterator<Item = SignOrZero>>(iter: I) -> Self {
        SignatureLike(iter.into_iter().collect())
    }
}

impl<T> FromIterator<i8> for SignatureLike<T>
where
    T: From<usize>,
{
    fn from_iter<I: IntoIterator<Item = i8>>(iter: I) -> Self {
        SignatureLike(
            iter.into_iter()
                .map(|x| match x {
                    0 => SignOrZero::Zero,
                    1 => SignOrZero::Plus,
                    -1 => SignOrZero::Minus,
                    _ => panic!("Invalid value for Signature"),
                })
                .collect(),
        )
    }
}

impl<T> FromIterator<isize> for SignatureLike<T>
where
    T: From<usize>,
{
    fn from_iter<I: IntoIterator<Item = isize>>(iter: I) -> Self {
        SignatureLike(
            iter.into_iter()
                .map(|x| match x {
                    0 => SignOrZero::Zero,
                    1 => SignOrZero::Plus,
                    -1 => SignOrZero::Minus,
                    _ => panic!("Invalid value for Signature"),
                })
                .collect(),
        )
    }
}

impl<T> From<Vec<i8>> for SignatureLike<T>
where
    T: From<usize>,
{
    fn from(value: Vec<i8>) -> Self {
        SignatureLike::from_iter(value)
    }
}

impl<T> IntoIterator for SignatureLike<T>
where
    T: From<usize>,
{
    type Item = SignOrZero;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a SignatureLike<T>
where
    T: From<usize>,
{
    type Item = SignOrZero;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, Self::Item>>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().copied()
    }
}

impl<T> SignatureLike<T>
where
    T: From<usize> + Copy,
    usize: From<T>,
{
    pub fn iter_enumerated(&self) -> impl Iterator<Item = (T, &SignOrZero)> {
        self.0.iter_enumerated()
    }
    pub fn validate_basis<B>(&self, basis: &[B]) -> bool {
        self.len() == basis.len()
    }

    pub fn sum(&mut self, other: &Self) {
        for (i, sign) in other.iter_enumerated() {
            match (self[i], sign) {
                (SignOrZero::Zero, SignOrZero::Zero) => self.0[i] = SignOrZero::Zero,
                (SignOrZero::Zero, SignOrZero::Plus) => self.0[i] = SignOrZero::Plus,
                (SignOrZero::Zero, SignOrZero::Minus) => self.0[i] = SignOrZero::Minus,
                (SignOrZero::Plus, SignOrZero::Zero) => self.0[i] = SignOrZero::Plus,
                (SignOrZero::Plus, SignOrZero::Plus) => panic!("cannot add two positive signs"),
                (SignOrZero::Plus, SignOrZero::Minus) => self.0[i] = SignOrZero::Zero,
                (SignOrZero::Minus, SignOrZero::Zero) => self.0[i] = SignOrZero::Minus,
                (SignOrZero::Minus, SignOrZero::Plus) => self.0[i] = SignOrZero::Zero,
                (SignOrZero::Minus, SignOrZero::Minus) => panic!("cannot add two negative signs"),
            }
        }
    }

    pub fn panic_validate_basis<B>(&self, basis: &[B]) {
        if !self.validate_basis(basis) {
            panic!(
                "Invalid basis for Signature, expected length {}, got length {}",
                self.len(),
                basis.len()
            );
        }
    }

    pub fn to_momtrop_format(&self) -> Vec<isize> {
        self.0
            .iter()
            .map(|x| match x {
                SignOrZero::Zero => 0,
                SignOrZero::Plus => 1,
                SignOrZero::Minus => -1,
            })
            .collect()
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> std::slice::Iter<SignOrZero> {
        self.0.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn apply<B>(&self, basis: &[B]) -> B
    where
        B: RefZero + Clone + Neg<Output = B> + AddAssign<B>,
    {
        // self.panic_validate_basis(basis);
        let mut result = basis[0].ref_zero();
        for (&sign, t) in self.0.iter().zip(basis.iter().cloned()) {
            result += sign * t;
        }
        result
    }

    pub fn apply_typed<O, I, V>(&self, basis: &V) -> O
    where
        V: Index<I, Output = O>,
        O: RefZero + Neg<Output = O> + AddAssign<O> + Clone,
        I: From<usize>,
        usize: From<I>,
    {
        let mut result = basis[I::from(0)].ref_zero();
        for (&sign, i) in self.0.iter().zip(0..) {
            result += sign * basis[I::from(i)].clone();
        }

        result
    }

    pub fn apply_iter<I, O>(&self, basis: I) -> Option<O>
    where
        I: IntoIterator,
        I::Item: RefZero<O>,
        O: Clone + SubAssign<I::Item> + AddAssign<I::Item>,
    {
        let mut basis_iter = basis.into_iter();
        let mut signature_iter = self.into_iter();

        while let (Some(sign), Some(item)) = (signature_iter.next(), basis_iter.next()) {
            if sign.is_sign() {
                // Initialize the result based on the first non-zero sign
                let mut result = item.ref_zero();
                match sign {
                    SignOrZero::Zero => {
                        panic!("unreachable");
                        // return None;
                    }
                    SignOrZero::Plus => {
                        result += item;
                    }
                    SignOrZero::Minus => {
                        result -= item;
                    }
                }

                // Continue processing the rest of the iterator
                while let (Some(sign), Some(item)) = (signature_iter.next(), basis_iter.next()) {
                    match sign {
                        SignOrZero::Zero => {}
                        SignOrZero::Plus => {
                            result += item;
                        }
                        SignOrZero::Minus => {
                            result -= item;
                        }
                    }
                }

                return Some(result);
            }
        }

        // Return None if no non-zero sign was found
        None
    }

    pub fn label_with(&self, label: &str) -> String {
        let mut result = String::new();
        let mut first = true;
        for (i, sign) in self.0.iter().enumerate() {
            if !first {
                result.push_str(&sign.to_string());
            } else {
                first = false;
            }
            if sign.is_sign() {
                result.push_str(&format!("{}_{}", label, i));
            }
        }
        result
    }
}

#[test]
fn test_signature() {
    use crate::momentum::FourMomentum;
    let sig = LoopSignature::from_iter(vec![SignOrZero::Plus, SignOrZero::Minus]);
    let basis: Vec<i32> = vec![1, 2];
    assert_eq!(sig.apply(&basis), 1 - 2);
    assert_eq!(sig.apply_iter(basis.iter()), Some(-1));

    let basis: [FourMomentum<i32>; 4] = [
        FourMomentum::from_args(1, 1, 0, 0),
        FourMomentum::from_args(1, 0, 1, 0),
        FourMomentum::from_args(1, 0, 0, 1),
        FourMomentum::from_args(1, 1, 1, 1),
    ];

    let sig = LoopSignature::from_iter(vec![
        SignOrZero::Plus,
        SignOrZero::Minus,
        SignOrZero::Zero,
        SignOrZero::Plus,
    ]);

    assert_eq!(sig.apply(&basis), FourMomentum::from_args(1, 2, 0, 1));
    let sig = ExternalSignature::from_iter(vec![
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Zero,
    ]);
    assert_eq!(sig.apply_iter(basis.iter()), None);
    let sig = ExternalSignature::from_iter(vec![
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Zero,
        SignOrZero::Minus,
    ]);
    assert_eq!(sig.apply_iter(basis.iter()), Some(-basis[3]));
}

impl LoopExtSignature {
    pub fn compute_momentum_untyped<'a, 'b: 'a, T>(
        &self,
        loop_moms: &'a [T],
        external_moms: &'b [T],
    ) -> T
    where
        T: RefZero + Clone + Neg<Output = T> + AddAssign<T>,
    {
        if loop_moms.is_empty() {
            return self.external.apply(external_moms);
        }
        if external_moms.is_empty() {
            return self.internal.apply(loop_moms);
        }
        let mut res = self.internal.apply(loop_moms);
        res += self.external.apply(external_moms);
        res
    }

    pub fn compute_momentum<L, E, M>(&self, loop_momenta: &L, external_momenta: &E) -> M
    where
        M: RefZero + Clone + Neg<Output = M> + AddAssign<M>,
        L: Index<LoopIndex, Output = M> + Length,
        E: Index<ExternalIndex, Output = M> + Length,
    {
        if loop_momenta.is_empty() {
            return self.external.apply_typed(external_momenta);
        }

        if external_momenta.is_empty() {
            return self.internal.apply_typed(loop_momenta);
        }

        let mut res = self.internal.apply_typed(loop_momenta);
        res += self.external.apply_typed(external_momenta);
        res
    }

    pub fn to_momtrop_format(&self) -> (Vec<isize>, Vec<isize>) {
        (
            self.internal.to_momtrop_format(),
            self.external.to_momtrop_format(),
        )
    }

    /// Usefull for debugging
    pub fn format_momentum(&self) -> String {
        let mut res = String::new();
        let mut first = true;

        for (i, sign) in (&self.internal).into_iter().enumerate() {
            if !first {
                res.push_str(&sign.to_string());
            } else {
                first = false;
            }
            if sign.is_sign() {
                res.push_str(&format!("k_{}", i));
            }
        }

        for (i, sign) in (&self.external).into_iter().enumerate() {
            if !first {
                res.push_str(&sign.to_string());
            } else {
                first = false;
            }
            if sign.is_sign() {
                res.push_str(&format!("l_{}", i));
            }
        }
        res
    }

    #[allow(unused)]
    pub fn compute_four_momentum_from_three<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) -> FourMomentum<F<T>> {
        let loop_moms = loop_moms
            .iter()
            .map(|m| m.clone().into_on_shell_four_momentum(None))
            .collect::<TiVec<LoopIndex, _>>();

        self.compute_momentum(&loop_moms, external_moms)
    }

    pub fn compute_three_momentum_from_four<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) -> ThreeMomentum<F<T>> {
        let external_moms: ExternalThreeMomenta<F<T>> =
            external_moms.iter().map(|m| m.spatial.clone()).collect();
        self.compute_momentum(loop_moms, &external_moms)
    }
}
