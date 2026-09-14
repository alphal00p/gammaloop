#![allow(dead_code)]

use crate::momentum::sample::{
    ExternalFourMomenta, ExternalIndex, ExternalThreeMomenta, LoopIndex, LoopMomenta,
};
use crate::momentum::{FourMomentum, SignOrZero, ThreeMomentum};
use crate::utils::{F, FloatLike, Length};
use bincode::{BorrowDecode, Decode, Encode};
use serde::{Deserialize, Serialize};
use spenso::algebra::algebraic_traits::RefZero;
use std::fmt::Display;
use std::ops::{Add, AddAssign, Index, IndexMut, Neg, SubAssign};
use symbolica::atom::{Atom, AtomOrView, FunctionBuilder, Symbol};
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

impl LoopExtSignature {
    pub(crate) fn swap_loops(&mut self, i: LoopIndex, j: LoopIndex) {
        // println!("i{i},j{j}");÷
        if !self.internal.is_empty() {
            self.internal.0.swap(i, j);
        }
    }

    pub(crate) fn put_loop_to_ext(&mut self, l: LoopIndex) {
        let a = self.internal.0.remove(l);
        self.external.0.push(a);
    }
    pub(crate) fn swap_external(&mut self, i: ExternalIndex, j: ExternalIndex) {
        self.external.0.swap(i, j);
    }

    pub(crate) fn loop_atom<'a, I>(
        &self,
        mom_symbol: Symbol,
        additional_args: &'a [I],
        id_map: impl Fn(LoopIndex) -> Atom,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.internal.atom(mom_symbol, additional_args, id_map)
    }

    pub(crate) fn ext_atom<'a, I>(
        &self,
        mom_symbol: Symbol,
        additional_args: &'a [I],
        id_map: impl Fn(ExternalIndex) -> Atom,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        self.external.atom(mom_symbol, additional_args, id_map)
    }

    pub(crate) fn equality_up_to_sign(&self, other: &Self) -> bool {
        self == other
            || (self.internal.len() == other.internal.len()
                && self.external.len() == other.external.len()
                && self
                    .internal
                    .iter()
                    .zip(other.internal.iter())
                    .all(|(a, b)| *a == -*b)
                && self
                    .external
                    .iter()
                    .zip(other.external.iter())
                    .all(|(a, b)| *a == -*b))
    }

    /// Sampling may identify energies in an exact represented spatial frame.
    /// This does not change structural routing equality used by raised edges.
    pub(crate) fn spatial_equality_up_to_sign<T: FloatLike>(
        &self,
        other: &Self,
        externals: &ExternalThreeMomenta<F<T>>,
    ) -> eyre::Result<bool> {
        use crate::integrands::process::sampling_maps::SamplingEvaluationError;
        use rug::{Float, float::Round};
        const PRECISION: u32 = 2048;
        if self.internal.len() != other.internal.len()
            || self.external.len() != other.external.len()
            || self.external.len() != externals.len()
        {
            return Err(eyre::eyre!(
                "spatial route comparison requires equal routing dimensions and the complete external frame"
            ));
        }
        if externals
            .iter()
            .flat_map(|p| [&p.px, &p.py, &p.pz])
            .any(|value| !value.0.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "spatial route comparison",
                detail: "nonfinite original external component".into(),
            }
            .into());
        }
        if self.equality_up_to_sign(other) {
            return Ok(true);
        }
        let mut uncertain = false;
        for sign in [1, -1] {
            if !self
                .internal
                .iter()
                .zip(other.internal.iter())
                .all(|(a, b)| *a == if sign == 1 { *b } else { -*b })
            {
                continue;
            }
            // Differences can be +/-2, which SignatureLike cannot represent.
            let coefficients = self
                .external
                .to_momtrop_format()
                .into_iter()
                .zip(other.external.to_momtrop_format())
                .map(|(a, b)| a - sign * b)
                .collect::<Vec<_>>();
            let mut sums: [[Float; 2]; 3] = std::array::from_fn(|_| {
                [Float::with_val(PRECISION, 0), Float::with_val(PRECISION, 0)]
            });
            for (&coefficient, momentum) in coefficients.iter().zip(externals) {
                if coefficient == 0 {
                    continue;
                }
                for (sum, component) in
                    sums.iter_mut()
                        .zip([&momentum.px, &momentum.py, &momentum.pz])
                {
                    let (mut lo, mut hi) = component.0.mpfr_enclosure(PRECISION);
                    if coefficient < 0 {
                        std::mem::swap(&mut lo, &mut hi);
                    }
                    let lo = Float::with_val_round(PRECISION, &lo * coefficient, Round::Down).0;
                    let hi = Float::with_val_round(PRECISION, &hi * coefficient, Round::Up).0;
                    sum[0] = Float::with_val_round(PRECISION, &sum[0] + &lo, Round::Down).0;
                    sum[1] = Float::with_val_round(PRECISION, &sum[1] + &hi, Round::Up).0;
                }
            }
            if sums.iter().flatten().any(|value| !value.is_finite()) {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "spatial route comparison",
                    detail: "nonfinite directed external sum".into(),
                }
                .into());
            }
            if sums.iter().any(|[lo, hi]| lo > &0 || hi < &0) {
                continue;
            }
            if sums.iter().all(|[lo, hi]| lo == &0 && hi == &0) {
                return Ok(true);
            }
            uncertain = true;
        }
        if uncertain {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: format!("spatial routing equality is unresolved at {PRECISION}-bit directed precision: {self:?} versus {other:?}"),
            }
            .into());
        }
        Ok(false)
    }
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

impl<T: From<usize>> Default for SignatureLike<T> {
    fn default() -> Self {
        SignatureLike(TiVec::new())
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

impl<T> AddAssign<()> for SignatureLike<T>
where
    T: From<usize> + Copy,
    usize: From<T>,
{
    fn add_assign(&mut self, _rhs: ()) {
        self.0.push(SignOrZero::Plus);
    }
}

impl<T> SignatureLike<T>
where
    T: From<usize> + Copy,
    usize: From<T>,
{
    pub(crate) fn atom<'a, I>(
        &self,
        mom_symbol: Symbol,
        additional_args: &'a [I],
        id_map: impl Fn(T) -> Atom,
    ) -> Atom
    where
        &'a I: Into<AtomOrView<'a>>,
    {
        let mut rep = Atom::Zero;
        for (l, s) in self.iter_enumerated() {
            let mom = FunctionBuilder::new(mom_symbol)
                .add_arg(id_map(l))
                .add_args(additional_args)
                .finish();
            // println!("mom: {mom} {s}");

            rep += *s * mom;
        }

        // println!("rep{rep}");

        rep
    }

    pub(crate) fn iter_enumerated(&self) -> impl Iterator<Item = (T, &SignOrZero)> {
        self.0.iter_enumerated()
    }
    pub(crate) fn validate_basis<B>(&self, basis: &[B]) -> bool {
        self.len() == basis.len()
    }

    pub(crate) fn sum(&mut self, other: &Self) {
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

    pub(crate) fn panic_validate_basis<B>(&self, basis: &[B]) {
        if !self.validate_basis(basis) {
            panic!(
                "Invalid basis for Signature, expected length {}, got length {}",
                self.len(),
                basis.len()
            );
        }
    }

    pub(crate) fn to_momtrop_format(&self) -> Vec<isize> {
        self.0
            .iter()
            .map(|x| match x {
                SignOrZero::Zero => 0,
                SignOrZero::Plus => 1,
                SignOrZero::Minus => -1,
            })
            .collect()
    }
    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&'_ self) -> std::slice::Iter<'_, SignOrZero> {
        self.0.iter()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub(crate) fn apply<B>(&self, basis: &[B]) -> B
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

    pub(crate) fn try_apply<B>(&self, basis: &[B]) -> Option<B>
    where
        B: Clone + Neg<Output = B> + Add<B, Output = B>,
    {
        self.0
            .iter()
            .zip(basis.iter().cloned())
            .filter_map(|(sign, t)| match sign {
                SignOrZero::Zero => None,
                SignOrZero::Plus => Some(t),
                SignOrZero::Minus => Some(-t),
            })
            .reduce(|sum, t| sum + t)
    }

    pub(crate) fn apply_typed<O, I, V>(&self, basis: &V) -> O
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

    pub(crate) fn apply_iter<I, O>(&self, basis: I) -> Option<O>
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

    pub(crate) fn label_with(&self, label: &str) -> String {
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

    /// Canonization function to compare two signatures up to an overall sign,
    /// If the first nonzero entry is positive, it will return itself,
    /// otherwise it will return the negative of itself.
    pub(crate) fn first_abs(&self) -> Self {
        let sign = self.iter().find(|x| x.is_sign());

        if let Some(sign) = sign {
            if sign.is_positive() {
                self.clone()
            } else {
                self.iter().map(|x| -*x).collect()
            }
        } else {
            self.clone()
        }
    }

    pub(crate) fn pop(&mut self) -> Option<SignOrZero> {
        self.0.pop()
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
    pub(crate) fn compute_momentum_untyped<'a, 'b: 'a, T>(
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

    pub(crate) fn try_compute_momentum<'a, 'b: 'a, T>(
        &self,
        loop_moms: &'a [T],
        external_moms: &'b [T],
    ) -> Option<T>
    where
        T: Clone + Neg<Output = T> + Add<T, Output = T>,
    {
        let loop_part = self.internal.try_apply(loop_moms);
        let external_part = self.external.try_apply(external_moms);

        match (loop_part, external_part) {
            (Some(l), Some(e)) => Some(l + e),
            (Some(l), None) => Some(l),
            (None, Some(e)) => Some(e),
            (None, None) => None,
        }
    }

    pub(crate) fn compute_momentum<L, E, M>(&self, loop_momenta: &L, external_momenta: &E) -> M
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

    pub(crate) fn to_momtrop_format(&self) -> (Vec<isize>, Vec<isize>) {
        (
            self.internal.to_momtrop_format(),
            self.external.to_momtrop_format(),
        )
    }

    /// Usefull for debugging
    pub(crate) fn format_momentum(&self) -> String {
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
    pub(crate) fn compute_four_momentum_from_three<T: FloatLike>(
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

    pub(crate) fn compute_three_momentum_from_four<T: FloatLike>(
        &self,
        loop_moms: &LoopMomenta<F<T>>,
        external_moms: &ExternalFourMomenta<F<T>>,
    ) -> ThreeMomentum<F<T>> {
        let external_moms: ExternalThreeMomenta<F<T>> =
            external_moms.iter().map(|m| m.spatial.clone()).collect();
        self.compute_momentum(loop_moms, &external_moms)
    }
}

#[cfg(test)]
mod tests {
    use super::LoopExtSignature;

    #[test]
    fn affine_signature_equality_uses_one_global_sign() {
        let a = LoopExtSignature::from((vec![1], vec![1]));
        let independent_reversal = LoopExtSignature::from((vec![-1], vec![1]));
        assert!(!a.equality_up_to_sign(&independent_reversal));
        assert!(a.equality_up_to_sign(&LoopExtSignature::from((vec![-1], vec![-1]))));
        for (internal, external) in [
            (vec![0, 1], vec![0, 0]),
            (vec![0, 0], vec![-1, 1]),
            (vec![0, 0], vec![0, 0]),
            (vec![], vec![]),
        ] {
            let reversed = LoopExtSignature::from((
                internal.iter().map(|value| -*value).collect(),
                external.iter().map(|value| -*value).collect(),
            ));
            assert!(LoopExtSignature::from((internal, external)).equality_up_to_sign(&reversed));
        }
        assert!(
            !LoopExtSignature::from((vec![1], vec![]))
                .equality_up_to_sign(&LoopExtSignature::from((vec![], vec![1])))
        );
        assert!(!a.equality_up_to_sign(&LoopExtSignature::from((vec![1, 0], vec![1]))));
    }

    #[test]
    fn spatial_routing_equality_certifies_native_external_frames() {
        use crate::{
            integrands::process::sampling_maps::SamplingEvaluationError,
            momentum::{ThreeMomentum, sample::ExternalThreeMomenta},
            utils::{ArbPrec, F, FloatLike, QuadFloat},
        };
        fn check<T: FloatLike>() {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let left = LoopExtSignature::from((vec![1], vec![0, 0]));
            let right = LoopExtSignature::from((vec![1], vec![1, 1]));
            let beam = ThreeMomentum::new(one.from_i64(2), zero.clone(), one.from_i64(3));
            let cm = ExternalThreeMomenta::from_iter([beam.clone(), -beam]);
            assert!(!left.equality_up_to_sign(&right));
            assert!(left.spatial_equality_up_to_sign(&right, &cm).unwrap());
            let reversed = LoopExtSignature::from((vec![-1], vec![-1, -1]));
            assert!(left.spatial_equality_up_to_sign(&reversed, &cm).unwrap());

            // A common boost and even a tiny nonzero represented displacement
            // break the fixed-frame identity; neither has a norm tolerance.
            for displacement in [one.clone(), one.epsilon().square()] {
                let mut shifted = cm.clone();
                for momentum in &mut shifted {
                    // Preserve the tiny original input on a zero-valued axis.
                    momentum.py += &displacement;
                }
                assert!(!left.spatial_equality_up_to_sign(&right, &shifted).unwrap());
            }
            let a = LoopExtSignature::from((vec![1], vec![1, 0]));
            let wrong_sign = LoopExtSignature::from((vec![-1], vec![1, 0]));
            assert!(!a.spatial_equality_up_to_sign(&wrong_sign, &cm).unwrap());
            // With matching internal signs, the external difference is +2.
            let doubled = LoopExtSignature::from((vec![1], vec![-1, 0]));
            assert!(!a.spatial_equality_up_to_sign(&doubled, &cm).unwrap());

            let big = &one / one.epsilon().square();
            let medium = big.sqrt();
            let values = [big.clone(), medium.clone(), one.clone(), -medium, -big];
            assert_eq!(
                values.iter().fold(zero.clone(), |sum, value| sum + value),
                zero
            );
            let cancelled = ExternalThreeMomenta::from_iter(
                values.map(|x| ThreeMomentum::new(x, zero.clone(), zero.clone())),
            );
            let empty = LoopExtSignature::from((vec![1], vec![0; 5]));
            let sum = LoopExtSignature::from((vec![1], vec![1; 5]));
            assert!(!empty.spatial_equality_up_to_sign(&sum, &cancelled).unwrap());

            // These original binary64 values span more than the fixed 2048-bit
            // certificate budget, also when represented exactly in Quad/Arb.
            let huge = F(T::from_f64_exact_binary(2.0_f64.powi(1023)));
            let tiny = F(T::from_f64_exact_binary(f64::from_bits(1)));
            let unresolved = ExternalThreeMomenta::from_iter(
                [huge.clone(), tiny, -huge]
                    .map(|x| ThreeMomentum::new(x, zero.clone(), zero.clone())),
            );
            let empty = LoopExtSignature::from((vec![1], vec![0; 3]));
            let sum = LoopExtSignature::from((vec![1], vec![1; 3]));
            let error = empty
                .spatial_equality_up_to_sign(&sum, &unresolved)
                .unwrap_err();
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));

            let mut nonfinite = cm.clone();
            nonfinite.iter_mut().next().unwrap().px = F(T::from_f64_exact_binary(f64::INFINITY));
            // Input validation precedes even the structural fast path.
            let error = left
                .spatial_equality_up_to_sign(&left, &nonfinite)
                .unwrap_err();
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::Unrepresentable { .. })
            ));
            assert!(
                left.spatial_equality_up_to_sign(&right, &ExternalThreeMomenta::<F<T>>::default())
                    .is_err()
            );
        }
        check::<f64>();
        check::<QuadFloat>();
        check::<ArbPrec>();
    }

    #[test]
    fn raised_edge_comparison_preserves_complete_affine_routing() {
        use crate::{
            dot,
            graph::{Graph, parse::from_dot::IntoGraph},
            initialisation::test_initialise,
        };
        use linnet::half_edge::involution::EdgeIndex;
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph raised_affine_signature {
            ext [style=invis]
            edge [num=1 mass=0]
            node [num=1]
            ext -> a [id=0]
            a -> b [id=1]
            b -> c [id=2]
            c -> a [id=3]
            ext -> c [id=4]
        })
        .unwrap();
        let groups = graph.get_raised_edge_groups();
        assert!(groups.iter().any(|group| group.len() >= 2));
        for group in groups {
            for pair in group.windows(2) {
                assert!(graph.loop_momentum_basis.edges_are_raised(pair[0], pair[1]));
            }
        }
        // Deliberately supplied affine rows isolate the caller contract; they
        // are not a new physical routing assignment for this graph.
        let mut basis = graph.loop_momentum_basis.clone();
        basis.edge_signatures[EdgeIndex(1)] = LoopExtSignature::from((vec![1], vec![1, 0]));
        basis.edge_signatures[EdgeIndex(2)] = LoopExtSignature::from((vec![-1], vec![1, 0]));
        assert!(!basis.edges_are_raised(EdgeIndex(1), EdgeIndex(2)));
        basis.edge_signatures[EdgeIndex(2)] = LoopExtSignature::from((vec![-1], vec![-1, 0]));
        assert!(basis.edges_are_raised(EdgeIndex(1), EdgeIndex(2)));
    }
}
