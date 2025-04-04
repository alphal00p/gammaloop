use crate::cff::expression::CFFFloat;
use crate::momentum::{FourMomentum, Signature, ThreeMomentum};
use crate::numerator::NumeratorEvaluateFloat;
use crate::SamplingSettings;
use crate::{ParameterizationMapping, ParameterizationMode, Settings, MAX_LOOP};
use bincode::{Decode, Encode};
use colored::Colorize;

use itertools::{izip, Itertools};
use rand::Rng;
use ref_ops::{RefAdd, RefDiv, RefMul, RefNeg, RefRem, RefSub};
use rug::float::{Constant, ParseFloatError};
use rug::ops::{CompleteRound, Pow};
use rug::Float;
use serde::{Deserialize, Deserializer, Serialize};
use spenso::complex::SymbolicaComplex;
use spenso::{
    complex::{Complex, R},
    contraction::{RefOne, RefZero},
    upgrading_arithmetic::TrySmallestUpgrade,
};
use symbolica::atom::Symbol;
use symbolica::domains::float::{
    ConstructibleFloat, NumericalFloatLike, RealNumberLike, SingleFloat,
};
use symbolica::domains::integer::Integer;
use symbolica::evaluate::CompiledEvaluatorFloat;
use symbolica::symb;

use statrs::function::gamma::{gamma, gamma_lr, gamma_ur};
use std::cmp::{Ord, Ordering};
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, Sub, SubAssign};
use std::str::FromStr;
use std::sync::LazyLock;
use std::time::Duration;
use symbolica::domains::float::Real;
use symbolica::domains::rational::Rational;
// use symbolica::domains::Field;
use symbolica::numerical_integration::Sample;
use typed_index_collections::TiSlice;

#[allow(unused_imports)]
use log::{debug, info};
use symbolica::atom::Atom;
use symbolica::printer::{AtomPrinter, PrintOptions};

use git_version::git_version;
pub const GIT_VERSION: &str = git_version!(fallback = "unavailable");
pub const VERSION: &str = "0.0.1";

#[allow(unused)]
const MAX_DIMENSION: usize = MAX_LOOP * 3;

pub const PINCH_TEST_THRESHOLD: f64 = 1e-10;

pub const LEFT: usize = 0;
pub const RIGHT: usize = 1;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum Side {
    LEFT = 0,
    RIGHT = 1,
}

impl From<Side> for usize {
    fn from(val: Side) -> Self {
        val as usize
    }
}

pub mod sorted_vectorize {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use std::iter::FromIterator;

    pub fn serialize<'a, T, K, V, S>(target: T, ser: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
        T: IntoIterator<Item = (&'a K, &'a V)>,
        K: Serialize + PartialOrd + 'a,
        V: Serialize + PartialOrd + 'a,
    {
        let mut container: Vec<_> = target.into_iter().collect::<Vec<_>>();
        container.sort_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap());
        serde::Serialize::serialize(&container, ser)
    }

    pub fn deserialize<'de, T, K, V, D>(des: D) -> Result<T, D::Error>
    where
        D: Deserializer<'de>,
        T: FromIterator<(K, V)>,
        K: Deserialize<'de>,
        V: Deserialize<'de>,
    {
        let container: Vec<_> = serde::Deserialize::deserialize(des)?;
        Ok(T::from_iter(container))
    }
}

pub trait FloatConvertFrom<U> {
    fn convert_from(x: &U) -> Self;
}

//     fn convert_from(x: &f128::f128) -> f64 {
//         (*x).to_f64().unwrap()
//     }
// }

// impl FloatConvertFrom<f128::f128> for f128::f128 {
//     fn convert_from(x: &f128::f128) -> f128::f128 {
//         *x
//     }
// }

// impl FloatConvertFrom<f64> for f128::f128 {
//     fn convert_from(x: &f64) -> f128::f128 {
//         f128::f128::from_f64(*x).unwrap()
//     }
// }

#[derive(Debug, Clone, PartialEq, PartialOrd, Encode, Decode)]
pub struct VarFloat<const N: u32> {
    #[bincode(with_serde)]
    float: rug::Float,
}

impl<const N: u32> Serialize for VarFloat<N> {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let string = self.float.to_string();
        string.serialize(serializer)
    }
}

impl<'de, const N: u32> Deserialize<'de> for VarFloat<N> {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let string: String = serde::Deserialize::deserialize(deserializer)?;
        let val: Self = string
            .parse()
            .unwrap_or_else(|_| panic!("failed to parse arb prec from string: {}", string));

        Ok(val)
    }
}

impl<const N: u32> FromStr for VarFloat<N> {
    type Err = ParseFloatError;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        let float = rug::Float::parse(s)?;
        Ok(Self {
            float: rug::Float::with_val(N, float),
        })
    }
}

impl<const N: u32> Rem<&VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn rem(self, rhs: &VarFloat<N>) -> Self::Output {
        (&self.float % &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> From<Float> for VarFloat<N> {
    fn from(x: Float) -> Self {
        VarFloat {
            float: rug::Float::with_val(N, x),
        }
    }
}

impl<const N: u32> From<&Rational> for VarFloat<N> {
    fn from(x: &Rational) -> Self {
        let n = x.numerator();

        let n = match n {
            Integer::Double(f) => Float::with_val(N, f),
            Integer::Large(f) => Float::with_val(N, f),
            Integer::Natural(f) => Float::with_val(N, f),
        };

        let d = x.denominator();

        let d = match d {
            Integer::Double(f) => Float::with_val(N, f),
            Integer::Large(f) => Float::with_val(N, f),
            Integer::Natural(f) => Float::with_val(N, f),
        };

        let r = n / d;

        VarFloat {
            float: rug::Float::with_val(N, r),
        }
    }
}

impl<const N: u32> std::ops::Mul for VarFloat<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        (self.float * rhs.float).into()
    }
}

impl<const N: u32> std::ops::Mul<&VarFloat<N>> for VarFloat<N> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        (self.float * &rhs.float).into()
    }
}

impl<const N: u32> std::ops::Mul<VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn mul(self, rhs: VarFloat<N>) -> Self::Output {
        (&self.float * &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::Mul<&VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn mul(self, rhs: &VarFloat<N>) -> Self::Output {
        (&self.float * &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::MulAssign for VarFloat<N> {
    fn mul_assign(&mut self, rhs: Self) {
        self.float *= rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::MulAssign<&VarFloat<N>> for VarFloat<N> {
    fn mul_assign(&mut self, rhs: &Self) {
        self.float *= &rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::Add for VarFloat<N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        (self.float + rhs.float).into()
    }
}

impl<const N: u32> std::ops::Add<&VarFloat<N>> for VarFloat<N> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        (self.float + &rhs.float).into()
    }
}

impl<const N: u32> std::ops::Add<VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn add(self, rhs: VarFloat<N>) -> Self::Output {
        (&self.float + &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::Add<&VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn add(self, rhs: &VarFloat<N>) -> Self::Output {
        (&self.float + &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::AddAssign for VarFloat<N> {
    fn add_assign(&mut self, rhs: Self) {
        self.float += rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::AddAssign<&VarFloat<N>> for VarFloat<N> {
    fn add_assign(&mut self, rhs: &Self) {
        self.float += &rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::Sub for VarFloat<N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        (self.float - rhs.float).into()
    }
}

impl<const N: u32> std::ops::Sub<&VarFloat<N>> for VarFloat<N> {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        (self.float - &rhs.float).into()
    }
}

impl<const N: u32> std::ops::Sub<VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn sub(self, rhs: VarFloat<N>) -> Self::Output {
        (&self.float - &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::Sub<&VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn sub(self, rhs: &VarFloat<N>) -> Self::Output {
        (&self.float - &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::SubAssign for VarFloat<N> {
    fn sub_assign(&mut self, rhs: Self) {
        self.float -= rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::SubAssign<&VarFloat<N>> for VarFloat<N> {
    fn sub_assign(&mut self, rhs: &Self) {
        self.float -= &rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> Div for VarFloat<N> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        (self.float / rhs.float).into()
    }
}

impl<const N: u32> Div<&VarFloat<N>> for VarFloat<N> {
    type Output = Self;

    fn div(self, rhs: &Self) -> Self::Output {
        (self.float / &rhs.float).into()
    }
}

impl<const N: u32> Div<VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn div(self, rhs: VarFloat<N>) -> Self::Output {
        (&self.float / &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> Div<&VarFloat<N>> for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn div(self, rhs: &VarFloat<N>) -> Self::Output {
        (&self.float / &rhs.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::ops::DivAssign for VarFloat<N> {
    fn div_assign(&mut self, rhs: Self) {
        self.float /= rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::DivAssign<&VarFloat<N>> for VarFloat<N> {
    fn div_assign(&mut self, rhs: &Self) {
        self.float /= &rhs.float;
        self.float.set_prec(N);
    }
}

impl<const N: u32> std::ops::Neg for VarFloat<N> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        (-self.float).into()
    }
}

impl<const N: u32> std::ops::Neg for &VarFloat<N> {
    type Output = VarFloat<N>;

    fn neg(self) -> Self::Output {
        (-&self.float).complete(N as i64).into()
    }
}

impl<const N: u32> std::fmt::Display for VarFloat<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.float)
    }
}

impl<const N: u32> std::fmt::LowerExp for VarFloat<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:e}", self.float)
    }
}

impl<const N: u32> RefZero for VarFloat<N> {
    fn ref_zero(&self) -> Self {
        Self::new_zero()
    }
}

impl<const N: u32> R for VarFloat<N> {}

impl<const N: u32> RefOne for VarFloat<N> {
    fn ref_one(&self) -> Self {
        self.one()
    }
}

impl<const N: u32> NumericalFloatLike for VarFloat<N> {
    fn mul_add(&self, a: &Self, b: &Self) -> Self {
        (&self.float * &a.float + &b.float)
            .complete(N as i64)
            .into()
    }

    // fn norm(&self) -> Self {
    //     self.float.clone().abs().into()
    // }

    fn from_i64(&self, a: i64) -> Self {
        VarFloat {
            float: Float::with_val(N, a),
        }
    }

    fn from_usize(&self, a: usize) -> Self {
        VarFloat {
            float: Float::with_val(N, a),
        }
    }

    fn get_precision(&self) -> u32 {
        N
    }

    fn zero(&self) -> Self {
        Self::new_zero()
    }

    fn one(&self) -> Self {
        self.from_i64(1)
    }

    fn new_zero() -> Self {
        VarFloat {
            float: Float::new(N),
        }
    }

    fn inv(&self) -> Self {
        self.float.clone().recip().into()
    }

    fn pow(&self, e: u64) -> Self {
        rug::ops::Pow::pow(&self.float, e).complete(N as i64).into()
    }

    fn sample_unit<R: Rng + ?Sized>(&self, rng: &mut R) -> Self {
        let f: f64 = rng.gen();
        Float::with_val(N, f).into()
    }

    fn neg(&self) -> Self {
        (-self.float.clone()).into()
    }

    fn get_epsilon(&self) -> f64 {
        2.0f64.powi(-(N as i32))
    }

    fn fixed_precision(&self) -> bool {
        true
    }
}

impl<const N: u32> SingleFloat for VarFloat<N> {
    fn is_finite(&self) -> bool {
        self.float.is_finite()
    }

    fn is_one(&self) -> bool {
        self.float == 1.
    }

    fn is_zero(&self) -> bool {
        self.float == 0.
    }

    fn from_rational(&self, rat: &Rational) -> Self {
        rat.into()
    }
}
impl<const N: u32> RealNumberLike for VarFloat<N> {
    fn to_f64(&self) -> f64 {
        self.float.to_f64()
    }

    fn round_to_nearest_integer(&self) -> Integer {
        self.float.clone().round().to_integer().unwrap().into()
    }

    fn to_usize_clamped(&self) -> usize {
        self.float
            .to_integer()
            .unwrap()
            .to_usize()
            .unwrap_or(usize::MAX)
    }
}

impl<const N: u32> Real for VarFloat<N> {
    fn i(&self) -> Option<Self> {
        None
    }
    #[inline(always)]
    fn pi(&self) -> Self {
        Float::with_val(N, rug::float::Constant::Pi).into()
    }

    #[inline(always)]
    fn e(&self) -> Self {
        self.one().exp()
    }

    #[inline(always)]
    fn euler(&self) -> Self {
        Float::with_val(N, rug::float::Constant::Euler).into()
    }

    #[inline(always)]
    fn phi(&self) -> Self {
        (self.one() + self.from_i64(5).sqrt()) / Self::from_f64(2.)
    }
    fn atan2(&self, x: &Self) -> Self {
        self.float.clone().atan2(&x.float).into()
    }

    fn powf(&self, e: &Self) -> Self {
        self.float.clone().pow(e.float.clone()).into()
    }

    fn log(&self) -> Self {
        self.float.ln_ref().complete(N as i64).into()
    }
    fn norm(&self) -> Self {
        self.float.clone().abs().into()
    }

    delegate! {
        #[into]
        to self.float.clone(){
            fn sqrt(&self) -> Self;
            fn exp(&self) -> Self;
            fn sin(&self) -> Self;
            fn cos(&self) -> Self;
            fn tan(&self) -> Self;
            fn asin(&self) -> Self;
            fn acos(&self) -> Self;
            fn sinh(&self) -> Self;
            fn cosh(&self) -> Self;
            fn tanh(&self) -> Self;
            fn asinh(&self) -> Self;
            fn acosh(&self) -> Self;
            fn atanh(&self) -> Self;
        }
    }
}

impl FloatLike for VarFloat<113> {
    fn E(&self) -> Self {
        Self::E()
    }

    fn PIHALF(&self) -> Self {
        Self::PIHALF()
    }

    fn SQRT_2(&self) -> Self {
        Self::from_f64(2.0).sqrt()
    }

    fn SQRT_2_HALF(&self) -> Self {
        Self::from_f64(2.0).sqrt() / Self::from_f64(2.0)
    }

    fn rem_euclid(&self, rhs: &Self) -> Self {
        let r = self.ref_rem(rhs);
        if r < r.zero() {
            r + rhs
        } else {
            r
        }
    }

    fn FRAC_1_PI(&self) -> Self {
        Self::FRAC_1_PI()
    }

    fn PI(&self) -> Self {
        Self::PI()
    }

    fn TAU(&self) -> Self {
        Self::TAU()
    }

    fn from_f64(x: f64) -> Self {
        VarFloat::from_f64(x)
    }

    fn into_f64(&self) -> f64 {
        self.to_f64()
    }

    fn is_nan(&self) -> bool {
        self.float.is_nan()
    }

    fn is_infinite(&self) -> bool {
        self.float.is_infinite()
    }

    fn floor(&self) -> Self {
        self.float.clone().floor().into()
    }
}

impl<const N: u32> VarFloat<N> {
    fn one() -> Self {
        VarFloat {
            float: rug::Float::with_val(N, 1.0),
        }
    }
    #[allow(non_snake_case)]
    fn E() -> Self {
        Self::one().exp()
    }

    #[allow(non_snake_case)]
    fn PIHALF() -> Self {
        Self::PI() / Self::from_f64(2.0)
    }

    #[allow(non_snake_case)]
    fn PI() -> Self {
        VarFloat {
            float: rug::Float::with_val(N, Constant::Pi),
        }
    }

    #[allow(non_snake_case)]
    fn TAU() -> Self {
        let mut tau = Self::PI() + Self::PI();
        tau.float.set_prec(N);
        tau
    }

    #[allow(non_snake_case)]
    fn FRAC_1_PI() -> Self {
        Self::PI().inv()
    }

    pub fn from_f64(x: f64) -> Self {
        VarFloat {
            float: rug::Float::with_val(N, x),
        }
    }

    pub fn to_f64(&self) -> f64 {
        self.float.to_f64()
    }
}

impl<const N: u32> Default for VarFloat<N> {
    fn default() -> Self {
        VarFloat {
            float: rug::Float::with_val(N, 0.0),
        }
    }
}

impl PrecisionUpgradable for f128 {
    type Higher = f128;
    type Lower = f64;

    fn higher(&self) -> Self::Higher {
        self.clone()
    }

    fn lower(&self) -> Self::Lower {
        self.to_f64()
    }
}

impl<T: Real + PrecisionUpgradable, H: Real, L: Real> PrecisionUpgradable for Complex<T>
where
    T: PrecisionUpgradable<Higher = H, Lower = L>,
{
    type Higher = Complex<H>;
    type Lower = Complex<L>;

    fn higher(&self) -> Self::Higher {
        Complex::new(self.re.higher(), self.im.higher())
    }

    fn lower(&self) -> Self::Lower {
        Complex::new(self.re.lower(), self.im.lower())
    }
}

// #[allow(non_camel_case_types)]
// pub type f256 = VarFloat<243>;

pub trait PrecisionUpgradable {
    type Higher;
    type Lower;

    fn higher(&self) -> Self::Higher;
    fn lower(&self) -> Self::Lower;
}

pub trait FloatLike:
    Real
    +R
    + PartialOrd
    + RealNumberLike
    + for<'a> RefAdd<&'a Self, Output = Self>
    // + for<'a> RefMutAdd<&'a Self, Output = Self>
    + RefAdd<Self, Output = Self>
    // + RefMutAdd<Self, Output = Self>
    + for<'a> RefMul<&'a Self, Output = Self>
    // + for<'a> RefMutMul<&'a Self, Output = Self>
    + RefMul<Self, Output = Self>
    // + RefMutMul<Self, Output = Self>
    + for<'a> RefSub<&'a Self, Output = Self>
    // + for<'a> RefMutSub<&'a Self, Output = Self>
    + RefSub<Self, Output = Self>
    // + RefMutSub<Self, Output = Self>
    + for<'a> RefDiv<&'a Self, Output = Self>
    // + for<'a> RefMutDiv<&'a Self, Output = Self> f64 doesn't have RefMutDiv
    + RefDiv<Self, Output = Self>
    + for<'a> RefRem<&'a Self, Output = Self>
    // + RefMutDiv<Self, Output = Self>
    + RefNeg<Output = Self>
    + RefZero
    + RefOne
    // + RefMutNeg<Output = Self> f64 doesn't have RefMutNeg
    + PrecisionUpgradable
    + Serialize
    + Display
    + CFFFloat<Self>
    + NumeratorEvaluateFloat
    {

    #[allow(non_snake_case)]
    fn PI(&self) -> Self;
    #[allow(non_snake_case)]
    fn E(&self) -> Self;
    #[allow(non_snake_case)]
    fn TAU(&self) -> Self;
    #[allow(non_snake_case)]
    fn SQRT_2(&self) -> Self;
    #[allow(non_snake_case)]
    fn SQRT_2_HALF(&self) -> Self;
    #[allow(non_snake_case)]
    fn PIHALF(&self) -> Self;
    #[allow(non_snake_case)]
    fn FRAC_1_PI(&self) -> Self;

    fn from_f64(x: f64) -> Self;

    #[allow(clippy::wrong_self_convention)]
    fn into_f64(&self) -> f64; // for inverse gamma in tropical sampling

    fn is_nan(&self) -> bool;

    fn is_infinite(&self) -> bool;

    fn floor(&self) -> Self;

    fn square(&self) -> Self {
        self.pow(2)
    }

    fn powi(&self, n: i32) -> Self {
        let absn = n.unsigned_abs() as u64;
        if n.is_negative() {
            self.pow(absn).inv()
        } else {
            self.pow(absn)
        }
    }

    fn epsilon(&self) -> Self {
        Self::from_f64(f64::EPSILON)
    }

    fn less_than_epsilon(&self) -> bool {
        self < &self.epsilon()
    }

    fn positive(&self) -> bool {
        self > &self.zero()
    }

    fn max_value(&self) -> Self {
        Self::from_f64(f64::MAX)
    }

    fn min_value(&self) -> Self {
        Self::from_f64(f64::MIN)
    }

    fn ln(&self) -> Self {
        panic!("ln not implemented for {:?}", self);
        // self.log() //FIXME
    }

    fn rem_euclid(&self, rhs: &Self) -> Self;
}

#[derive(
    Debug, Clone, PartialEq, PartialOrd, Copy, Default, Serialize, Deserialize, Encode, Decode, Hash,
)]
pub struct F<T: FloatLike>(pub T);

impl<T: FloatLike> R for F<T> {}

impl<T: FloatLike> Rem<&F<T>> for &F<T> {
    type Output = F<T>;

    fn rem(self, rhs: &F<T>) -> Self::Output {
        F(self.0.ref_rem(&rhs.0))
    }
}

impl<T: FloatLike> RefZero<F<T>> for &F<T> {
    fn ref_zero(&self) -> F<T> {
        F(self.0.ref_zero())
    }
}

impl<T: FloatLike> RealNumberLike for F<T> {
    delegate! {
        to self.0{
            fn to_usize_clamped(&self)->usize;
            fn to_f64(&self)->f64;
            fn round_to_nearest_integer(&self)->Integer;
        }
    }
}

impl<T: FloatLike> SingleFloat for F<T> {
    delegate! {
        to self.0{
            fn is_zero(&self)->bool;
            fn is_one(&self)->bool;
            fn is_finite(&self)->bool;
        }
    }
    fn from_rational(&self, rat: &Rational) -> Self {
        F(self.0.from_rational(rat))
    }
}

impl<T: FloatLike> PrecisionUpgradable for F<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    type Higher = F<T::Higher>;
    type Lower = F<T::Lower>;

    fn higher(&self) -> Self::Higher {
        F(self.0.higher())
    }

    fn lower(&self) -> Self::Lower {
        F(self.0.lower())
    }
}

impl<T: FloatLike> RefZero for F<T> {
    fn ref_zero(&self) -> Self {
        F(self.0.zero())
    }
}

impl<T: FloatLike> RefOne for F<T> {
    fn ref_one(&self) -> Self {
        F(self.0.one())
    }
}

impl<T: FloatLike> TrySmallestUpgrade<F<T>> for F<T> {
    type LCM = F<T>;
    fn try_upgrade(&self) -> Option<std::borrow::Cow<Self::LCM>> {
        Some(std::borrow::Cow::Borrowed(self))
    }
}

impl<T: FloatLike> TrySmallestUpgrade<F<T>> for Complex<F<T>> {
    type LCM = Complex<F<T>>;
    fn try_upgrade(&self) -> Option<std::borrow::Cow<Self::LCM>> {
        Some(std::borrow::Cow::Borrowed(self))
    }
}

// impl<T:FloatLike> TrySmallestUpgrade<Complex<F<T>>> for F<T> {
//     type LCM = Complex<F<T>>;
//     fn try_upgrade(&self) -> Option<std::borrow::Cow<Self::LCM>> {
//         Some(std::borrow::Cow::Borrowed(self))
//     }
// }

impl<'a, T: FloatLike> From<&'a Rational> for F<T> {
    fn from(x: &'a Rational) -> Self {
        F(T::from_f64(x.to_f64()))
    }
}

impl<T: FloatLike> From<T> for F<T> {
    fn from(x: T) -> Self {
        F(x)
    }
}

impl<T: FloatLike> std::fmt::Display for F<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl<T: FloatLike> std::fmt::LowerExp for F<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:e}", self.0)
    }
}

impl<T: FloatLike> NumericalFloatLike for F<T> {
    fn mul_add(&self, a: &Self, b: &Self) -> Self {
        F(self.0.mul_add(&a.0, &b.0))
    }
    fn new_zero() -> Self {
        F(T::new_zero())
    }
    fn sample_unit<R: Rng + ?Sized>(&self, rng: &mut R) -> Self {
        F(self.0.sample_unit(rng))
    }
    fn neg(&self) -> Self {
        F(self.0.ref_neg())
    }

    delegate! {
        #[into]
        to self.0{
            fn zero(&self) -> Self;
            fn one(&self) -> Self;
            // fn norm(&self) -> Self;
            fn from_usize(&self, x: usize) -> Self;
            fn from_i64(&self, x: i64) -> Self;
            fn pow(&self, n: u64) -> Self;
            fn inv(&self) -> Self;
            fn get_precision(&self) -> u32;
            fn get_epsilon(&self) -> f64;
            fn fixed_precision(&self) -> bool;
        }
    }
}

impl<T: FloatLike + ConstructibleFloat> ConstructibleFloat for F<T> {
    fn new_from_i64(a: i64) -> Self {
        F(T::new_from_i64(a))
    }

    fn new_from_usize(a: usize) -> Self {
        F(T::new_from_usize(a))
    }

    fn new_one() -> Self {
        F(T::new_one())
    }

    fn new_sample_unit<R: Rng + ?Sized>(rng: &mut R) -> Self {
        F(T::new_sample_unit(rng))
    }
}

impl<T: FloatLike> Real for F<T> {
    fn atan2(&self, x: &Self) -> Self {
        F(self.0.atan2(&x.0))
    }

    fn i(&self) -> Option<Self> {
        None
    }

    fn powf(&self, e: &Self) -> Self {
        F(self.0.powf(&e.0))
    }

    fn norm(&self) -> Self {
        F(self.0.norm())
    }

    delegate! {
        #[into]
        to self.0{
            fn e(&self)->Self;
            fn phi(&self)->Self;
            fn euler(&self)->Self;
            fn pi(&self)->Self;
            fn sqrt(&self) -> Self;
            fn log(&self) -> Self;
            fn exp(&self) -> Self;
            fn sin(&self) -> Self;
            fn cos(&self) -> Self;
            fn tan(&self) -> Self;
            fn asin(&self) -> Self;
            fn acos(&self) -> Self;
            fn sinh(&self) -> Self;
            fn cosh(&self) -> Self;
            fn tanh(&self) -> Self;
            fn asinh(&self) -> Self;
            fn acosh(&self) -> Self;
            fn atanh(&self) -> Self;

        }
    }
}

use delegate::delegate;

impl<T: FloatLike> F<T> {
    pub fn max(self, other: F<T>) -> F<T> {
        if self < other {
            other
        } else {
            self
        }
    }

    pub fn negate(&mut self) {
        self.0 = -self.0.clone();
    }

    pub fn from_ff64(x: F<f64>) -> Self {
        F(T::from_f64(x.0))
    }

    pub fn higher(&self) -> F<T::Higher>
    where
        T::Higher: FloatLike,
    {
        F(self.0.higher())
    }

    pub fn lower(&self) -> F<T::Lower>
    where
        T::Lower: FloatLike,
    {
        F(self.0.lower())
    }

    pub fn from_f64(x: f64) -> Self {
        F(T::from_f64(x))
    }

    pub fn into_ff64(&self) -> F<f64> {
        F(self.0.into_f64())
    }

    pub fn abs(&self) -> Self {
        F(self.0.norm())
    }

    pub fn i(&self) -> Complex<Self> {
        Complex::new(self.zero(), self.one())
    }

    pub fn log10(&self) -> Self {
        self.ln()
    }

    pub fn complex_sqrt(&self) -> Complex<Self> {
        if self.positive() {
            Complex::new(self.sqrt(), self.zero())
        } else {
            Complex::new(self.zero(), (-self).sqrt())
        }
    }

    pub fn rem_euclid(&self, rhs: &Self) -> Self {
        F(self.0.rem_euclid(&rhs.0))
    }

    delegate! {
        #[into]
        to self.0 {
            #[allow(non_snake_case)]
            pub fn PI(&self) -> Self;
            #[allow(non_snake_case)]
            pub fn E(&self) -> Self;
            #[allow(non_snake_case)]
            pub fn TAU(&self) -> Self;
            #[allow(non_snake_case)]
            pub fn PIHALF(&self) -> Self;
            #[allow(non_snake_case)]
            pub fn SQRT_2(&self) -> Self;
            #[allow(non_snake_case)]
            pub fn SQRT_2_HALF(&self) -> Self;
            #[allow(non_snake_case)]
            pub fn FRAC_1_PI(&self) -> Self;
            pub fn into_f64(&self) -> f64;
            pub fn square(&self) -> Self;
            pub fn powi(&self, n: i32) -> Self;
            pub fn epsilon(&self) -> Self;
            pub fn less_than_epsilon(&self) -> bool;
            pub fn positive(&self) -> bool;
            pub fn max_value(&self) -> Self;
            pub fn min_value(&self) -> Self;
            pub fn ln(&self) -> Self;
            pub fn is_nan(&self) -> bool;
            pub fn is_infinite(&self) -> bool;
            pub fn floor(&self) -> Self;
        }
    }
}

impl CompiledEvaluatorFloat for F<f64> {
    fn evaluate(
        eval: &mut symbolica::evaluate::CompiledEvaluator,
        args: &[Self],
        out: &mut [Self],
    ) {
        // cast to f64
        let args_f64: Vec<f64> = args.iter().map(|x| x.0).collect_vec();
        let mut out_f64 = out.iter().map(|x| x.0).collect_vec();

        eval.evaluate_double(&args_f64, &mut out_f64);

        // write the result to out
        out.iter_mut()
            .zip(out_f64)
            .for_each(|(out_ff64, out_f64)| *out_ff64 = F(out_f64));
    }
}

impl<T: FloatLike> Add<F<T>> for F<T> {
    type Output = F<T>;
    fn add(self, rhs: F<T>) -> Self::Output {
        F(self.0 + rhs.0)
    }
}

impl<T: FloatLike> Add<&F<T>> for F<T> {
    type Output = F<T>;
    fn add(self, rhs: &F<T>) -> Self::Output {
        F(self.0 + &rhs.0)
    }
}

impl<T: FloatLike> Add<&F<T>> for &F<T> {
    type Output = F<T>;
    fn add(self, rhs: &F<T>) -> Self::Output {
        F(self.0.ref_add(&rhs.0))
    }
}

impl<T: FloatLike> Add<F<T>> for &F<T> {
    type Output = F<T>;
    fn add(self, rhs: F<T>) -> Self::Output {
        F(self.0.ref_add(rhs.0))
    }
}

impl<T: FloatLike> AddAssign<&F<T>> for F<T> {
    fn add_assign(&mut self, rhs: &F<T>) {
        self.0 += &rhs.0;
    }
}

impl<T: FloatLike> AddAssign<F<T>> for F<T> {
    fn add_assign(&mut self, rhs: F<T>) {
        self.0 += rhs.0;
    }
}

impl<T: FloatLike> Sub<F<T>> for F<T> {
    type Output = F<T>;
    fn sub(self, rhs: F<T>) -> Self::Output {
        F(self.0 - rhs.0)
    }
}

impl<T: FloatLike> Sub<&F<T>> for F<T> {
    type Output = F<T>;
    fn sub(self, rhs: &F<T>) -> Self::Output {
        F(self.0 - &rhs.0)
    }
}

impl<T: FloatLike> Sub<&F<T>> for &F<T> {
    type Output = F<T>;
    fn sub(self, rhs: &F<T>) -> Self::Output {
        F(self.0.ref_sub(&rhs.0))
    }
}

impl<T: FloatLike> Sub<F<T>> for &F<T> {
    type Output = F<T>;
    fn sub(self, rhs: F<T>) -> Self::Output {
        F(self.0.ref_sub(rhs.0))
    }
}

impl<T: FloatLike> SubAssign<&F<T>> for F<T> {
    fn sub_assign(&mut self, rhs: &F<T>) {
        self.0 -= &rhs.0;
    }
}

impl<T: FloatLike> SubAssign<F<T>> for F<T> {
    fn sub_assign(&mut self, rhs: F<T>) {
        self.0 -= rhs.0;
    }
}

impl<T: FloatLike> Mul<F<T>> for F<T> {
    type Output = F<T>;
    fn mul(self, rhs: F<T>) -> Self::Output {
        F(self.0 * rhs.0)
    }
}

impl<T: FloatLike> Mul<&F<T>> for F<T> {
    type Output = F<T>;
    fn mul(self, rhs: &F<T>) -> Self::Output {
        F(self.0 * &rhs.0)
    }
}

impl<T: FloatLike> Mul<&F<T>> for &F<T> {
    type Output = F<T>;
    fn mul(self, rhs: &F<T>) -> Self::Output {
        F(self.0.ref_mul(&rhs.0))
    }
}

impl<T: FloatLike> Mul<F<T>> for &F<T> {
    type Output = F<T>;
    fn mul(self, rhs: F<T>) -> Self::Output {
        F(self.0.ref_mul(rhs.0))
    }
}

impl<T: FloatLike> MulAssign<&F<T>> for F<T> {
    fn mul_assign(&mut self, rhs: &F<T>) {
        self.0 *= &rhs.0;
    }
}

impl<T: FloatLike> MulAssign<F<T>> for F<T> {
    fn mul_assign(&mut self, rhs: F<T>) {
        self.0 *= rhs.0;
    }
}

impl<T: FloatLike> Div<F<T>> for F<T> {
    type Output = F<T>;
    fn div(self, rhs: F<T>) -> Self::Output {
        F(self.0 / rhs.0)
    }
}

impl<T: FloatLike> Div<&F<T>> for F<T> {
    type Output = F<T>;
    fn div(self, rhs: &F<T>) -> Self::Output {
        F(self.0 / &rhs.0)
    }
}

impl<T: FloatLike> Div<&F<T>> for &F<T> {
    type Output = F<T>;
    fn div(self, rhs: &F<T>) -> Self::Output {
        F(self.0.ref_div(&rhs.0))
    }
}

impl<T: FloatLike> Div<F<T>> for &F<T> {
    type Output = F<T>;
    fn div(self, rhs: F<T>) -> Self::Output {
        F(self.0.ref_div(rhs.0))
    }
}

impl<T: FloatLike> DivAssign<&F<T>> for F<T> {
    fn div_assign(&mut self, rhs: &F<T>) {
        self.0 /= &rhs.0;
    }
}

impl<T: FloatLike> DivAssign<F<T>> for F<T> {
    fn div_assign(&mut self, rhs: F<T>) {
        self.0 /= rhs.0;
    }
}

impl<T: FloatLike> Neg for F<T> {
    type Output = F<T>;
    fn neg(self) -> Self::Output {
        F(-self.0)
    }
}

impl<T: FloatLike> Neg for &F<T> {
    type Output = F<T>;
    fn neg(self) -> Self::Output {
        F(self.0.ref_neg())
    }
}

pub trait RefDefault {
    fn default(&self) -> Self;
}

impl<T: FloatLike> RefDefault for T {
    fn default(&self) -> Self {
        self.zero()
    }
}

impl<T: FloatLike> RefDefault for F<T> {
    fn default(&self) -> Self {
        F(self.0.default())
    }
}
impl PrecisionUpgradable for f64 {
    type Higher = f128;
    type Lower = f64;

    fn higher(&self) -> Self::Higher {
        f128::from_f64(*self)
    }

    fn lower(&self) -> Self::Lower {
        *self
    }
}

impl FloatLike for f64 {
    fn PI(&self) -> Self {
        std::f64::consts::PI
    }

    fn SQRT_2(&self) -> Self {
        std::f64::consts::SQRT_2
    }

    fn SQRT_2_HALF(&self) -> Self {
        std::f64::consts::SQRT_2 / 2.0
    }

    fn PIHALF(&self) -> Self {
        std::f64::consts::PI / 2.0
    }

    fn E(&self) -> Self {
        std::f64::consts::E
    }

    fn TAU(&self) -> Self {
        std::f64::consts::TAU
    }

    fn FRAC_1_PI(&self) -> Self {
        std::f64::consts::FRAC_1_PI
    }

    fn from_f64(x: f64) -> Self {
        x
    }

    fn into_f64(&self) -> f64 {
        *self
    }

    fn is_nan(&self) -> bool {
        f64::is_nan(*self)
    }

    fn is_infinite(&self) -> bool {
        f64::is_infinite(*self)
    }

    fn floor(&self) -> Self {
        f64::floor(*self)
    }

    fn rem_euclid(&self, rhs: &Self) -> Self {
        f64::rem_euclid(*self, *rhs)
    }
}
impl From<F<f64>> for f64 {
    fn from(value: F<f64>) -> Self {
        value.0
    }
}

impl From<F<f64>> for Rational {
    fn from(value: F<f64>) -> Self {
        value.0.into()
    }
}

#[allow(non_camel_case_types)]
pub type f128 = VarFloat<113>;

/// An iterator which iterates two other iterators simultaneously
#[derive(Clone, Debug)]
#[must_use = "iterator adaptors are lazy and do nothing unless consumed"]
pub struct ZipEq<I, J> {
    a: I,
    b: J,
}

/// An iterator which iterates two other iterators simultaneously and checks
/// if the sizes are equal in debug mode.
#[allow(unused)]
pub fn zip_eq<I, J>(i: I, j: J) -> ZipEq<I::IntoIter, J::IntoIter>
where
    I: IntoIterator,
    J: IntoIterator,
{
    ZipEq {
        a: i.into_iter(),
        b: j.into_iter(),
    }
}

impl<I, J> Iterator for ZipEq<I, J>
where
    I: Iterator,
    J: Iterator,
{
    type Item = (I::Item, J::Item);

    fn next(&mut self) -> Option<Self::Item> {
        match (self.a.next(), self.b.next()) {
            (None, None) => None,
            (Some(a), Some(b)) => Some((a, b)),
            (None, Some(_)) => {
                #[cfg(debug_assertions)]
                panic!("Unequal length of iterators; first iterator finished first");
                #[cfg(not(debug_assertions))]
                None
            }
            (Some(_), None) => {
                #[cfg(debug_assertions)]
                panic!("Unequal length of iterators; second iterator finished first");
                #[cfg(not(debug_assertions))]
                None
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let sa = self.a.size_hint();
        let sb = self.b.size_hint();
        (sa.0.min(sb.0), sa.1.zip(sb.1).map(|(ua, ub)| ua.min(ub)))
    }
}

impl<I, J> ExactSizeIterator for ZipEq<I, J>
where
    I: ExactSizeIterator,
    J: ExactSizeIterator,
{
}

pub fn parse_python_expression(expression: &str) -> Atom {
    let processed_string = String::from(expression)
        .replace("**", "^")
        .replace("cmath.sqrt", "sqrt")
        .replace("cmath.pi", "pi")
        .replace("math.sqrt", "sqrt")
        .replace("math.pi", "pi");
    Atom::parse(processed_string.as_str())
        .map_err(|e| {
            format!(
                "Failed to parse expression : '{}'\nError: {}",
                processed_string, e
            )
        })
        .unwrap()
}

pub fn to_str_expression(expression: &Atom) -> String {
    format!(
        "{}",
        AtomPrinter::new_with_options(
            expression.as_view(),
            PrintOptions {
                pretty_matrix: true,
                terms_on_new_line: false,
                color_top_level_sum: false,
                color_builtin_symbols: false,
                print_finite_field: false,
                explicit_rational_polynomial: false,
                symmetric_representation_for_finite_field: false,
                number_thousands_separator: None,
                multiplication_operator: '*',
                square_brackets_for_function: false,
                num_exp_as_superscript: false,
                latex: false,
                double_star_for_exponentiation: false,
                precision: None,
            },
        )
    )
}

/// Format a mean ± sdev as mean(sdev) with the correct number of digits.
/// Based on the Python package gvar.
pub fn format_uncertainty(mean: F<f64>, sdev: F<f64>) -> String {
    let mean = mean.0;
    let sdev = sdev.0;

    fn ndec(x: f64, offset: usize) -> i32 {
        let mut ans = (offset as f64 - x.log10()) as i32;
        if ans > 0 && x * 10.0.powi(ans) >= [0.5, 9.5, 99.5][offset] {
            ans -= 1;
        }
        if ans < 0 {
            0
        } else {
            ans
        }
    }
    let v = mean;
    let dv = sdev.abs();

    // special cases
    if v.is_nan() || dv.is_nan() {
        format!("{:e} ± {:e}", v, dv)
    } else if dv.is_infinite() {
        format!("{:e} ± inf", v)
    } else if v.is_zero() && !(1e-4..1e5).contains(&dv) {
        if dv.is_zero() {
            "0(0)".to_owned()
        } else {
            let e = format!("{:.1e}", dv);
            let mut ans = e.split('e');
            let e1 = ans.next().unwrap();
            let e2 = ans.next().unwrap();
            "0.0(".to_owned() + e1 + ")e" + e2
        }
    } else if v.is_zero() {
        if dv >= 9.95 {
            format!("0({:.0})", dv)
        } else if dv >= 0.995 {
            format!("0.0({:.1})", dv)
        } else {
            let ndecimal = ndec(dv, 2);
            format!(
                "{:.*}({:.0})",
                ndecimal as usize,
                v,
                dv * (10.).powi(ndecimal)
            )
        }
    } else if dv.is_zero() {
        let e = format!("{:e}", v);
        let mut ans = e.split('e');
        let e1 = ans.next().unwrap();
        let e2 = ans.next().unwrap();
        if e2 != "0" {
            e1.to_owned() + "(0)e" + e2
        } else {
            e1.to_owned() + "(0)"
        }
    } else if dv > 1e4 * v.abs() {
        format!("{:.1e} ± {:.2e}", v, dv)
    } else if v.abs() >= 1e6 || v.abs() < 1e-5 {
        // exponential notation for large |self.mean|
        let exponent = v.abs().log10().floor();
        let fac = (10.0).powf(&exponent);
        let mantissa = format_uncertainty(F(v / fac), F(dv / fac));
        let e = format!("{:.0e}", fac);
        let mut ee = e.split('e');
        mantissa + "e" + ee.nth(1).unwrap()
    }
    // normal cases
    else if dv >= 9.95 {
        if v.abs() >= 9.5 {
            format!("{:.0}({:.0})", v, dv)
        } else {
            let ndecimal = ndec(v.abs(), 1);
            format!("{:.*}({:.*})", ndecimal as usize, v, ndecimal as usize, dv)
        }
    } else if dv >= 0.995 {
        if v.abs() >= 0.95 {
            format!("{:.1}({:.1})", v, dv)
        } else {
            let ndecimal = ndec(v.abs(), 1);
            format!("{:.*}({:.*})", ndecimal as usize, v, ndecimal as usize, dv)
        }
    } else {
        let ndecimal = ndec(v.abs(), 1).max(ndec(dv, 2));
        format!(
            "{:.*}({:.0})",
            ndecimal as usize,
            v,
            dv * (10.).powi(ndecimal)
        )
    }
}

/// Compare two slices, selecting on length first
#[allow(unused)]
pub fn compare_slice<T: Ord>(slice1: &[T], slice2: &[T]) -> Ordering {
    match slice1.len().cmp(&slice2.len()) {
        Ordering::Equal => (),
        non_eq => return non_eq,
    }

    let l = slice1.len();
    // Slice to the loop iteration range to enable bound check
    // elimination in the compiler
    let lhs = &slice1[..l];
    let rhs = &slice2[..l];

    for i in 0..l {
        match lhs[i].cmp(&rhs[i]) {
            Ordering::Equal => (),
            non_eq => return non_eq,
        }
    }

    Ordering::Equal
}

pub trait Signum {
    fn multiply_sign(&self, sign: i8) -> Self;
}

// impl Signum for f128::f128 {
//     #[inline]
//     fn multiply_sign(&self, sign: i8) -> f128::f128 {
//         match sign {
//             1 => *self,
//             0 => f128::f128::zero(),
//             -1 => self.neg(),
//             _ => unreachable!("Sign should be -1,0,1"),
//         }
//     }
// }

impl Signum for f64 {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> f64 {
        match sign {
            1 => *self,
            0 => self.zero(),
            -1 => -self,
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl Signum for f32 {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> Self {
        match sign {
            1 => *self,
            0 => 0.0,
            -1 => self.neg(),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: FloatLike> Signum for Complex<F<T>> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> Complex<F<T>> {
        match sign {
            1 => self.clone(),
            0 => self.ref_zero(),
            -1 => -self.clone(),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: FloatLike> Signum for FourMomentum<F<T>> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> FourMomentum<F<T>> {
        match sign {
            1 => self.clone(),
            0 => self.ref_zero(),
            -1 => -self,
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

#[allow(unused)]
#[inline]
/// Invert with better precision
pub fn finv<T: FloatLike>(c: Complex<F<T>>) -> Complex<F<T>> {
    let norm = c.norm_squared();
    c.conj() / norm
}

#[allow(unused)]
#[inline]
pub fn powi<T: FloatLike>(c: Complex<F<T>>, n: i32) -> Complex<F<T>> {
    if n.is_negative() {
        let u = -n as u64;
        finv(c.pow(u))
    } else {
        let u = n as u64;
        c.pow(u)
    }
}

#[allow(unused)]
pub fn evaluate_signature<T>(signature: &[i8], momenta: &[FourMomentum<F<T>>]) -> FourMomentum<F<T>>
where
    T: FloatLike,
{
    let mut momentum = momenta[0].zero();
    for (&sign, mom) in zip_eq(signature, momenta) {
        match sign {
            0 => {}
            1 => momentum += mom,
            -1 => momentum -= mom,
            _ => {
                #[cfg(debug_assertions)]
                panic!("Sign should be -1,0,1")
            }
        }
    }

    momentum
}

#[allow(unused)]
#[inline]
pub fn pinch_dampening_function<T: FloatLike>(
    dampening_arg: F<T>,
    delta_t: F<T>,
    powers: (u64, u64),
    multiplier: f64,
) -> F<T> {
    // Make sure the function is even in t-tstar
    assert!(powers.1 % 2 == 0);
    let a = dampening_arg.pow(powers.0);
    &a / (&a + F::<T>::from_f64(multiplier) * delta_t.pow(powers.1))
}

pub fn h<T: FloatLike>(
    t: F<T>,
    tstar: Option<F<T>>,
    sigma: Option<F<T>>,
    h_function_settings: &crate::HFunctionSettings,
) -> F<T> {
    let sqrt_pi = t.PI().sqrt();
    let sig = if let Some(s) = sigma {
        s
    } else {
        F::<T>::from_f64(h_function_settings.sigma)
    };
    let power = h_function_settings.power;
    match h_function_settings.function {
        crate::HFunction::Exponential => {
            (-(t.square()) / (sig.square())).exp() * F::<T>::from_f64(2_f64) / (sqrt_pi * &sig)
        }
        crate::HFunction::PolyExponential => {
            // Result of \int_0^{\infty} dt (t/sigma)^{-p} exp(2-t^2/sigma^2-sigma^2/t^2)
            let normalisation = match power {
                None | Some(0) => sqrt_pi * &sig / F::<T>::from_f64(2_f64),
                Some(1) => F::<T>::from_f64(0.841_568_215_070_771_4) * &sig,
                Some(3) => F::<T>::from_f64(1.033_476_847_068_688_6) * &sig,
                Some(4) => F::<T>::from_f64(1.329_340_388_179_137) * &sig,
                Some(6) => F::<T>::from_f64(2.880_237_507_721_463_7) * &sig,
                Some(7) => F::<T>::from_f64(4.783_566_971_347_609) * &sig,
                Some(9) => F::<T>::from_f64(16.225_745_976_182_285) * &sig,
                Some(10) => F::<T>::from_f64(32.735_007_058_911_25) * &sig,
                Some(12) => F::<T>::from_f64(155.837_465_922_583_42) * &sig,
                Some(13) => F::<T>::from_f64(364.658_500_356_566_04) * &sig,
                Some(15) => F::<T>::from_f64(2_257.637_553_015_473) * &sig,
                Some(16) => F::<T>::from_f64(5_939.804_418_537_864) * &sig,
                _ => panic!(
                    "Value {} of power in poly exponential h function not supported",
                    power.unwrap()
                ),
            };
            let prefactor = match power {
                None | Some(0) => normalisation.inv(),
                Some(p) => (&t / &sig).powi(-(p as i32)) / normalisation,
            };
            prefactor
                * (F::<T>::from_f64(2_f64)
                    - (t.square()) / (sig.square())
                    - (sig.square()) / (t.square()))
                .exp()
        }
        crate::HFunction::PolyLeftRightExponential => {
            // Result of \int_0^{\infty} dt (t/sigma)^{-p} exp( -((t^2/sigma^2 +1)/ (t/sigma) -2) )
            let normalisation = match power {
                None | Some(0) => F::<T>::from_f64(2.066_953_694_137_377) * &sig,
                Some(1) => F::<T>::from_f64(1.683_136_430_141_542_8) * &sig,
                Some(3) => F::<T>::from_f64(3.750_090_124_278_92) * &sig,
                Some(4) => F::<T>::from_f64(9.567_133_942_695_218) * &sig,
                Some(6) => F::<T>::from_f64(139.373_101_752_153_5) * &sig,
                Some(7) => F::<T>::from_f64(729.317_000_713_132_1) * &sig,
                Some(9) => F::<T>::from_f64(32_336.242_742_929_753) * &sig,
                Some(10) => F::<T>::from_f64(263_205.217_049_469) * &sig,
                Some(12) => F::<T>::from_f64(2.427_503_717_893_097_5e7) * &sig,
                Some(13) => F::<T>::from_f64(2.694_265_921_644_289e8) * &sig,
                Some(15) => F::<T>::from_f64(9.040_742_057_760_125e12) * &sig,
                Some(16) => F::<T>::from_f64(1.452_517_480_246_491_3e14) * &sig,
                _ => panic!(
                    "Value {} of power in poly exponential h function not supported",
                    power.unwrap()
                ),
            };

            // println!("normalisation: {}", normalisation);
            // println!("t: {}", t);
            // println!("sig: {}", sig);
            // println!("power: {:?}", power);

            let prefactor = match power {
                None | Some(0) => normalisation.inv(),
                Some(p) => (&t / &sig).powi(-(p as i32)) / normalisation,
            };

            // println!("prefactor: {}", prefactor);
            prefactor
                * (F::<T>::from_f64(2_f64) - ((t.square()) / (sig.square()) + t.one()) / (t / sig))
                    .exp()
        }
        crate::HFunction::ExponentialCT => {
            let delta_t_sq = (tstar.clone().unwrap() - &t).square();
            let tstar_sq = tstar.unwrap().square();
            // info!("dampener: {}", dampener);
            // info!("delta_t_sq: {}", delta_t_sq);
            // info!("tstar_sq: {}", tstar_sq);
            // info!(
            //     "Exp arg: {}",
            //     -sig.inv() * (delta_t_sq / tstar_sq + sig * sig * (dampener * dampener))
            // );
            // info!(
            //     "result: {}",
            //     (-sig.inv() * (delta_t_sq / tstar_sq + sig * sig * (dampener * dampener))).exp()
            // );
            if h_function_settings.enabled_dampening {
                let dampener = delta_t_sq.clone() / (delta_t_sq.clone() - &tstar_sq);
                (-sig.inv() * (delta_t_sq.clone() / tstar_sq + sig.square() * (dampener.square())))
                    .exp()
            } else {
                (-sig.inv() * (delta_t_sq / tstar_sq)).exp()
            }
        }
    }
}

/// Calculate the determinant of any complex-valued input matrix using LU-decomposition.
/// Original C-code by W. Gong and D.E. Soper.
#[allow(unused)]
pub fn determinant<T: FloatLike>(bb: &[Complex<F<T>>], dimension: usize) -> Complex<F<T>> {
    let one = bb[0].re.one();
    let zero = one.zero();
    // Define matrix related variables.
    let mut determinant = Complex::new(one.clone(), zero.clone());
    let mut indx = [0; MAX_DIMENSION];
    let mut d = 1; // initialize parity parameter

    // Inintialize the matrix to be decomposed with the transferred matrix b.
    let mut aa = bb.to_vec();

    // Define parameters used in decomposition.
    let mut imax = 0;
    let mut flag = 1;
    let mut dumc;
    let mut sum;

    let mut aamax;
    let mut dumr;
    let mut vv = vec![zero.clone(); MAX_DIMENSION];

    // Get the implicit scaling information.
    for i in 0..dimension {
        aamax = zero.clone();
        for j in 0..dimension {
            let r = aa[i * dimension + j].norm_squared();
            if r > aamax {
                aamax = r;
            }
        }
        // Set a flag to check if the determinant is zero.
        if aamax.is_zero() {
            flag = 0;
        }
        // Save the scaling.
        vv[i] = aamax.inv();
    }
    if flag == 1 {
        for j in 0..dimension {
            for i in 0..j {
                sum = aa[i * dimension + j].clone();
                for k in 0..i {
                    sum -= &aa[i * dimension + k] * &aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
            }
            //Initialize for the search for largest pivot element.
            aamax = zero.clone();
            for i in j..dimension {
                sum = aa[i * dimension + j].clone();
                for k in 0..j {
                    sum -= &aa[i * dimension + k] * &aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum.clone();
                // Figure of merit for the pivot.
                dumr = &vv[i] * sum.norm_squared();
                // Is it better than the best so far?
                if dumr >= aamax {
                    imax = i;
                    aamax = dumr;
                }
            }
            // See if we need to interchange rows.
            if j != imax {
                for k in 0..dimension {
                    dumc = aa[imax * dimension + k].clone();
                    aa[imax * dimension + k] = aa[j * dimension + k].clone();
                    aa[j * dimension + k] = dumc.clone();
                }
                // Change the parity of d.
                d = -d;
                // Interchange the scale factor.
                vv[imax] = vv[j].clone();
            }
            indx[j] = imax;
            if j + 1 != dimension {
                dumc = aa[j * dimension + j].inv();
                for i in j + 1..dimension {
                    aa[i * dimension + j] *= &dumc;
                }
            }
        }
    }
    // Calculate the determinant using the decomposed matrix.
    if flag == 0 {
        determinant = Complex::new(zero.clone(), zero.clone());
    } else {
        // Multiply the diagonal elements.
        for diagonal in 0..dimension {
            determinant *= &aa[diagonal * dimension + diagonal];
        }
        determinant *= &zero.from_i64(d);
    }
    determinant
}

#[allow(unused)]
pub fn next_combination_with_replacement(state: &mut [usize], max_entry: usize) -> bool {
    for i in (0..state.len()).rev() {
        if state[i] < max_entry {
            state[i] += 1;
            for j in i + 1..state.len() {
                state[j] = state[i]
            }
            return true;
        }
    }
    false
}

pub fn compute_loop_part<T: FloatLike>(
    loop_signature: &Signature,
    loop_moms: &[ThreeMomentum<F<T>>],
) -> ThreeMomentum<F<T>> {
    loop_signature.apply(loop_moms)
}

pub fn compute_shift_part<T: FloatLike>(
    external_signature: &Signature,
    external_moms: &[FourMomentum<F<T>>],
) -> FourMomentum<F<T>> {
    external_signature.apply(external_moms)
}

pub fn compute_t_part_of_shift_part<T: FloatLike>(
    external_signature: &Signature,
    external_moms: &[FourMomentum<F<T>>],
) -> F<T> {
    // external_signature.panic_validate_basis(external_moms);
    external_signature
        .apply_iter(external_moms.iter().map(|m| m.temporal.value.clone()))
        .unwrap_or(external_moms[0].temporal.value.zero())
}

// Bilinear form for E-surface defined as sqrt[(k+p1)^2+m1sq] + sqrt[(k+p2)^2+m2sq] + e_shift
// The Bilinear system then reads 4 k.a.k + 4 k.n + C = 0
#[allow(unused, clippy::type_complexity)]
pub fn one_loop_e_surface_bilinear_form<T: FloatLike>(
    p1: &[F<T>; 3],
    p2: &[F<T>; 3],
    m1_sq: F<T>,
    m2_sq: F<T>,
    e_shift: F<T>,
) -> ([[F<T>; 3]; 3], [F<T>; 3], F<T>) {
    let zero = e_shift.zero();
    let two = zero.from_i64(2);
    let e_shift_sq = e_shift.square();
    let p1_sq = p1[0].square() + p1[1].square() + p1[2].square();
    let p2_sq = p2[0].square() + p2[1].square() + p2[2].square();

    let zeros = [zero.clone(), zero.clone(), zero.clone()];
    let mut a = [zeros.clone(), zeros.clone(), zeros.clone()];
    a[0][0] = (&p1[0] - &p2[0] - &e_shift) * (&p2[0] - &p1[0] - &e_shift);
    a[0][1] = (&p1[0] - &p2[0]) * (&p2[1] - &p1[1]);
    a[1][0] = a[0][1].clone();
    a[0][2] = (&p1[0] - &p2[0]) * (&p2[2] - &p1[2]);
    a[2][0] = a[0][2].clone();
    a[1][1] = (&p1[1] - &p2[1] - &e_shift) * (&p2[1] - &p1[1] - &e_shift);
    a[1][2] = (&p1[1] - &p2[1]) * (&p2[2] - &p1[2]);
    a[2][1] = a[1][2].clone();
    a[2][2] = (&p1[2] - &p2[2] - &e_shift) * (&p2[2] - &p1[2] - &e_shift);

    let mut b = zeros.clone();
    b[0] =
        (&p2[0] - &p1[0]) * (&m1_sq - &m2_sq + &p1_sq - &p2_sq) + &e_shift_sq * (&p1[0] + &p2[0]);
    b[1] =
        (&p2[1] - &p1[1]) * (&m1_sq - &m2_sq + &p1_sq - &p2_sq) + &e_shift_sq * (&p1[1] + &p2[1]);
    b[2] =
        (&p2[2] - &p1[2]) * (&m1_sq - &m2_sq + &p1_sq - &p2_sq) + &e_shift_sq * (&p1[2] + &p2[2]);

    let c = -&e_shift_sq * &e_shift_sq + &two * &e_shift_sq * (&m1_sq + &m2_sq + &p1_sq + &p2_sq)
        - (&m1_sq - &m2_sq + &p1_sq - &p2_sq) * (&m1_sq - &m2_sq + &p1_sq - &p2_sq);

    (a, b, c)
}

#[allow(unused)]
pub fn one_loop_e_surface_exists<T: FloatLike>(
    p1: &[F<T>; 3],
    p2: &[F<T>; 3],
    m1_sq: F<T>,
    m2_sq: F<T>,
    e_shift: F<T>,
) -> (bool, bool) {
    let p_norm_sq = (&p1[0] - &p2[0]) * (&p1[0] - &p2[0])
        + (&p1[1] - &p2[1]) * (&p1[1] - &p2[1])
        + (&p1[2] - &p2[2]) * (&p1[2] - &p2[2]);

    // /!\ In alphaLoop this should be done without numerical check but purely symbolically, or at least
    // mul_unit must make sure there is no weird transition from non-pinched to pinched for the 2->2 massless E-surface sandwich,
    // i.e. such cases should be *both* existing and pinched!
    if e_shift > F::<T>::from_f64(PINCH_TEST_THRESHOLD) {
        return (false, false);
    }
    let test = (e_shift.square() - &p_norm_sq)
        - (m1_sq.sqrt() + m2_sq.sqrt()) * (m1_sq.sqrt() + m2_sq.sqrt());
    if test.norm() < F::<T>::from_f64(PINCH_TEST_THRESHOLD) {
        (false, true)
    } else if test < e_shift.zero() {
        (false, false)
    } else {
        (true, false)
    }
}

use color_eyre::Result;
use eyre::eyre;
#[allow(unused)]
use std::fmt::LowerExp;

pub trait ApproxEq<U: LowerExp, T: LowerExp>: LowerExp {
    fn approx_eq(&self, other: &U, tolerance: &T) -> bool;

    fn approx_eq_slice(lhs: &[Self], rhs: &[U], tolerance: &T) -> bool
    where
        Self: Sized,
    {
        lhs.iter()
            .zip_eq(rhs)
            .all(|(l, r)| l.approx_eq(r, tolerance))
    }

    fn approx_eq_iterator<'a, I, J>(lhs: I, rhs: J, tolerance: &'a T) -> bool
    where
        Self: Sized + 'a,
        U: 'a,
        I: IntoIterator<Item = &'a Self>,
        J: IntoIterator<Item = &'a U>,
    {
        lhs.into_iter()
            .zip_eq(rhs)
            .all(|(l, r)| l.approx_eq(r, tolerance))
    }

    fn assert_approx_eq(&self, other: &U, tolerance: &T) {
        assert!(
            self.approx_eq(other, tolerance),
            "assert_approx_eq failed: \n{:+e} != \n{:+e} with tolerance {:+e}",
            self,
            other,
            tolerance
        )
    }
    fn approx_eq_res(&self, other: &U, tolerance: &T) -> Result<()> {
        if self.approx_eq(other, tolerance) {
            Ok(())
        } else {
            Err(eyre!(
                "assert_approx_eq failed: \n{:+e} != \n{:+e} with tolerance {:+e}",
                self,
                other,
                tolerance
            ))
        }
    }
}

// pub trait ApproxEqable: Real+PartialOrd+ for<'a> RefSub<&'a Self,Output = Self>+IsZero{}

// impl<T: Real+PartialOrd+ for<'a> RefSub<&'a T,Output = T>+IsZero> ApproxEqable for T{}

impl<T: FloatLike> ApproxEq<F<T>, F<T>> for F<T> {
    fn approx_eq(&self, other: &F<T>, tolerance: &F<T>) -> bool {
        if other.is_zero() {
            self.norm() < tolerance.clone()
        } else {
            ((self.ref_sub(other)) / other).norm() < tolerance.clone()
        }
    }
}

impl<T: FloatLike> ApproxEq<Complex<F<T>>, F<T>> for Complex<F<T>> {
    fn approx_eq(&self, other: &Complex<F<T>>, tolerance: &F<T>) -> bool {
        if !self.norm().re.approx_eq(&other.norm().re, tolerance) {
            return false;
        } else if self.norm().is_zero() || other.norm().is_zero() {
            return true;
        }
        let two_pi = self.re.PI() + self.re.PI();
        let arg_self = self.arg().rem_euclid(&two_pi);
        let arg_other = other.arg().rem_euclid(&two_pi);
        if !arg_self.approx_eq(&arg_other, tolerance) {
            return false;
        }
        true
    }
    fn approx_eq_res(&self, other: &Complex<F<T>>, tolerance: &F<T>) -> Result<()> {
        if !self.norm().re.approx_eq(&other.norm().re, tolerance) {
            return Err(eyre!(
                "Norms are not approximately equal: \n{:+e} != \n{:+e} with tolerance {:+e}",
                &self.norm().re,
                other.norm().re,
                tolerance
            ));
        } else if self.norm().is_zero() || other.norm().is_zero() {
            return Ok(());
        }

        let two_pi = self.re.PI() + self.re.PI();
        let arg_self = self.arg().rem_euclid(&two_pi);
        let arg_other = other.arg().rem_euclid(&two_pi);
        // let arg_diff = (&self.arg() - &other.arg()).rem_euclid(&two_pi);
        // let arg_zero = self.re.zero();
        if !arg_self.approx_eq(&arg_other, tolerance) {
            return  Err(eyre!(
                "Phases are not approximately equal: \n{:+e} - \n{:+e}= \n{:+e}!=0 with tolerance {:+e}",
                arg_self, arg_other,&arg_self-&arg_other, tolerance
            ));
        }
        Ok(())
    }
}

impl<T: FloatLike> ApproxEq<F<T>, F<T>> for Complex<F<T>> {
    fn approx_eq(&self, other: &F<T>, tolerance: &F<T>) -> bool {
        self.re.approx_eq(other, tolerance) && self.im.approx_eq(tolerance, tolerance)
    }
    fn approx_eq_res(&self, other: &F<T>, tolerance: &F<T>) -> Result<()> {
        if self.im.approx_eq(tolerance, tolerance) {
            return Err(eyre!(
                "Non-zero imaginary part: \n{:+e} with tolerance {:+e}",
                &self.im,
                tolerance
            ));
        }
        if !self.re.approx_eq(other, tolerance) {
            return Err(eyre!(
                "Real parts are not approximately equal: \n{:+e} != \n{:+e} with tolerance {:+e}",
                &self.re,
                other,
                tolerance
            ));
        }
        Ok(())
    }
}

impl<T: FloatLike> ApproxEq<Complex<F<T>>, F<T>> for F<T> {
    fn approx_eq(&self, other: &Complex<F<T>>, tolerance: &F<T>) -> bool {
        other.re.approx_eq(self, tolerance) && other.im.approx_eq(tolerance, tolerance)
    }

    fn approx_eq_res(&self, other: &Complex<F<T>>, tolerance: &F<T>) -> Result<()> {
        if other.im.approx_eq(tolerance, tolerance) {
            return Err(eyre!(
                "Non-zero imaginary part: \n{:+e} with tolerance {:+e}",
                &other.im,
                tolerance
            ));
        }
        if !other.re.approx_eq(self, tolerance) {
            return Err(eyre!(
                "Real parts are not approximately equal: \n{:+e} != \n{:+e} with tolerance {:+e}",
                &other.re,
                self,
                tolerance
            ));
        }
        Ok(())
    }
}

#[allow(unused)]
pub fn one_loop_eval_e_surf<T: FloatLike>(
    k: &[F<T>; 3],
    p1: &[F<T>; 3],
    p2: &[F<T>; 3],
    m1_sq: F<T>,
    m2_sq: F<T>,
    e_shift: F<T>,
) -> F<T> {
    ((&k[0] + &p1[0]) * (&k[0] + &p1[0])
        + (&k[1] + &p1[1]) * (&k[1] + &p1[1])
        + (&k[2] + &p1[2]) * (&k[2] + &p1[2])
        + m1_sq)
        .sqrt()
        + ((&k[0] + &p2[0]) * (&k[0] + &p2[0])
            + (&k[1] + &p2[1]) * (&k[1] + &p2[1])
            + (&k[2] + &p2[2]) * (&k[2] + &p2[2])
            + m2_sq)
            .sqrt()
        + e_shift
}

#[allow(unused)]
pub fn one_loop_eval_e_surf_k_derivative<T: FloatLike>(
    k: &[F<T>; 3],
    p1: &[F<T>; 3],
    p2: &[F<T>; 3],
    m1_sq: F<T>,
    m2_sq: F<T>,
) -> [F<T>; 3] {
    let e1 = ((&k[0] + &p1[0]) * (&k[0] + &p1[0])
        + (&k[1] + &p1[1]) * (&k[1] + &p1[1])
        + (&k[2] + &p1[2]) * (&k[2] + &p1[2])
        + m1_sq)
        .sqrt();
    let e2 = ((&k[0] + &p2[0]) * (&k[0] + &p2[0])
        + (&k[1] + &p2[1]) * (&k[1] + &p2[1])
        + (&k[2] + &p2[2]) * (&k[2] + &p2[2])
        + m2_sq)
        .sqrt();
    [
        (&k[0] + &p1[0]) / &e1 + (&k[0] + &p2[0]) / &e2,
        (&k[1] + &p1[1]) / &e1 + (&k[1] + &p2[1]) / &e2,
        (&k[2] + &p1[2]) / &e1 + (&k[2] + &p2[2]) / &e2,
    ]
}

#[allow(unused)]
pub fn one_loop_get_e_surf_t_scaling<T: FloatLike>(
    k: &[F<T>; 3],
    p1: &[F<T>; 3],
    p2: &[F<T>; 3],
    m1_sq: F<T>,
    m2_sq: F<T>,
    e_shift: F<T>,
) -> [F<T>; 2] {
    let zero = e_shift.zero();
    let one = zero.one();
    let (a, b, c_coef) = one_loop_e_surface_bilinear_form(p1, p2, m1_sq, m2_sq, e_shift);
    let mut a_coef = zero.clone();
    for i in 0..=2 {
        for j in 0..=2 {
            a_coef += &k[i] * &a[i][j] * &k[j];
        }
    }
    a_coef *= zero.from_i64(4);
    let mut b_coef = zero.clone();
    for i in 0..=2 {
        b_coef += &k[i] * &b[i];
    }
    b_coef *= zero.from_i64(4);
    let discr = b_coef.square() - zero.from_i64(4) * &a_coef * &c_coef;
    if discr < zero {
        [zero.clone(), zero.clone()]
    } else {
        [
            (-&b_coef + discr.sqrt()) / (zero.from_i64(2) * &a_coef),
            (-&b_coef - discr.sqrt()) / (zero.from_i64(2) * &a_coef),
        ]
    }
}

pub fn box_muller<T: FloatLike>(x1: F<T>, x2: F<T>) -> (F<T>, F<T>) {
    let r = (-x1.from_i64(2) * x1.log()).sqrt();
    let theta = r.from_i64(2) * r.PI() * x2;
    (r.clone() * theta.cos(), r * theta.sin())
}

pub fn compute_surface_and_volume<T: FloatLike>(n_dim: usize, radius: F<T>) -> (F<T>, F<T>) {
    let mut surface = radius.from_i64(2);
    let one = radius.one();
    let mut volume = one.clone();
    for i in 1..n_dim + 1 {
        (surface, volume) = (
            one.from_i64(2) * one.PI() * volume,
            surface / one.from_i64(i as i64),
        );
    }
    (
        surface * radius.pow(n_dim as u64),
        volume * radius.pow(n_dim as u64),
    )
}

pub fn get_n_dim_for_n_loop_momenta(
    settings: &Settings,
    n_loop_momenta: usize,
    force_radius: bool,
    n_edges: Option<usize>, // for tropical parameterization, we need to know the number of edges
) -> usize {
    if matches!(
        settings.sampling,
        SamplingSettings::DiscreteGraphs(crate::DiscreteGraphSamplingSettings::TropicalSampling(_))
    ) {
        let tropical_part = 2 * n_edges.expect("No tropical subgraph table generated, please run without tropical sampling or regenerate with tables") - 1;
        let d_l = 3 * n_loop_momenta;
        return if d_l % 2 == 1 {
            tropical_part + d_l + 1
        } else {
            tropical_part + d_l
        };
    }
    match settings.parameterization.mode {
        ParameterizationMode::HyperSphericalFlat => {
            // Because we use Box-Muller, we need to have an even number of angular dimensions
            let mut n_dim = 3 * n_loop_momenta;
            if n_dim % 2 == 1 {
                n_dim += 1;
            }
            // Then if the radius is not forced, then we need to add mul_unit more dimension
            if !force_radius {
                n_dim += 1;
            }
            n_dim
        }
        ParameterizationMode::HyperSpherical
        | ParameterizationMode::Cartesian
        | ParameterizationMode::Spherical => {
            if force_radius {
                3 * n_loop_momenta - 1
            } else {
                3 * n_loop_momenta
            }
        }
    }
}

pub fn global_parameterize<T: FloatLike>(
    x: &[F<T>],
    e_cm_squared: F<T>,
    settings: &Settings,
    force_radius: bool,
) -> (Vec<[F<T>; 3]>, F<T>) {
    let zero = e_cm_squared.zero();
    let one = zero.one();
    match settings.parameterization.mode {
        ParameterizationMode::HyperSpherical | ParameterizationMode::HyperSphericalFlat => {
            let e_cm =
                e_cm_squared.sqrt() * F::<T>::from_f64(settings.parameterization.shifts[0].0);
            let mut jac = one.clone();
            // rescale the input to the desired range
            let mut x_r = Vec::with_capacity(x.len());
            if !force_radius {
                x_r.push(x[0].clone());
            } else {
                let lo = F::<T>::from_f64(settings.parameterization.input_rescaling[0][0].0);
                let hi = F::<T>::from_f64(settings.parameterization.input_rescaling[0][0].1);
                x_r.push(&lo + &x[0] * (&hi - &lo));
                jac *= &hi - &lo;
            }
            let lo = F::<T>::from_f64(settings.parameterization.input_rescaling[0][1].0);
            let hi = F::<T>::from_f64(settings.parameterization.input_rescaling[0][1].1);
            x_r.push(&lo + &x[1] * (&hi - &lo));
            jac *= &hi - &lo;
            for xi in &x[2..] {
                let lo = F::<T>::from_f64(settings.parameterization.input_rescaling[0][2].0);
                let hi = F::<T>::from_f64(settings.parameterization.input_rescaling[0][2].1);
                x_r.push(&lo + xi * (&hi - &lo));
                jac *= &hi - &lo;
            }

            let radius: F<T> = if force_radius {
                x[0].clone()
            } else {
                match settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        // r = e_cm * ln(1 + b*x/(1-x))
                        let b = F::<T>::from_f64(settings.parameterization.b);
                        let radius = &e_cm * (&one + &b * &x_r[0] / (&one - &x_r[0])).log();
                        jac *= &e_cm * &b / (&one - &x_r[0]) / (&one + &x_r[0] * (&b - &one));
                        radius
                    }
                    ParameterizationMapping::Linear => {
                        // r = e_cm * b * x/(1-x)
                        let b = F::<T>::from_f64(settings.parameterization.b);
                        let radius = &e_cm * &b * &x_r[0] / (&one - &x_r[0]);
                        jac *= (&e_cm * &b + &radius).powi(2) / &e_cm / &b;
                        radius
                    }
                }
            };
            match settings.parameterization.mode {
                ParameterizationMode::HyperSpherical => {
                    let phi = zero.from_i64(2) * zero.PI() * &x_r[1];
                    jac *= zero.from_i64(2) * zero.PI();

                    let mut cos_thetas = Vec::with_capacity(x.len() - 2);
                    let mut sin_thetas = Vec::with_capacity(x.len() - 2);

                    for (i, xi) in x_r[2..].iter().enumerate() {
                        let cos_theta = -&one + zero.from_i64(2) * xi;
                        jac *= zero.from_i64(2);
                        let sin_theta = (&one - cos_theta.square()).sqrt();
                        if i > 0 {
                            jac *= sin_theta.pow(i as u64);
                        }
                        cos_thetas.push(cos_theta);
                        sin_thetas.push(sin_theta);
                    }

                    let mut concatenated_vecs = Vec::with_capacity(x.len() / 3);
                    let mut base = radius.clone();
                    for (cos_theta, sin_theta) in cos_thetas.iter().zip(sin_thetas.iter()) {
                        concatenated_vecs.push(&base * cos_theta);
                        base *= sin_theta;
                    }
                    concatenated_vecs.push(&base * phi.cos());
                    concatenated_vecs.push(&base * phi.sin());

                    jac *= radius.pow((x.len() - 1) as u64); // hyperspherical coords

                    (
                        concatenated_vecs
                            .chunks(3)
                            .map(|v| [v[0].clone(), v[1].clone(), v[2].clone()])
                            .collect(),
                        jac,
                    )
                }
                ParameterizationMode::HyperSphericalFlat => {
                    // As we will use Box Muller we expect an even number of random variables
                    assert!(x_r[1..].len() % 2 == 0);
                    let mut normal_distributed_xs = vec![];
                    for x_pair in x_r[1..].chunks(2) {
                        let (z1, z2) = box_muller(x_pair[0].clone(), x_pair[1].clone());
                        normal_distributed_xs.push(z1);
                        normal_distributed_xs.push(z2);
                    }
                    // ignore the last variable generated if we had to pad to get the 3*n_loop_momenta
                    if normal_distributed_xs.len() % 3 != 0 {
                        normal_distributed_xs.pop();
                    }
                    let curr_norm = normal_distributed_xs[..]
                        .iter()
                        .map(|x| x.square())
                        .reduce(|acc, e| acc + &e)
                        .unwrap_or(zero.clone())
                        .sqrt();
                    let surface =
                        compute_surface_and_volume(normal_distributed_xs.len() - 1, radius.clone())
                            .0;
                    jac *= surface;
                    let rescaling_factor = &radius / &curr_norm;
                    (
                        normal_distributed_xs
                            .chunks(3)
                            .map(|v| {
                                [
                                    &v[0] * &rescaling_factor,
                                    &v[1] * &rescaling_factor,
                                    &v[2] * &rescaling_factor,
                                ]
                            })
                            .collect(),
                        jac,
                    )
                }
                _ => unreachable!(),
            }
        }
        ParameterizationMode::Cartesian | ParameterizationMode::Spherical => {
            if force_radius {
                panic!("Cannot force radius for non-hyperspherical parameterization.");
            }
            let mut jac = one.clone();
            let mut vecs = Vec::with_capacity(x.len() / 3);
            for (i, xi) in x.chunks(3).enumerate() {
                let (vec_i, jac_i) = parameterize3d(xi, e_cm_squared.clone(), i, settings);
                vecs.push(vec_i);
                jac *= jac_i;
            }
            (vecs, jac)
        }
    }
}

#[allow(unused)]
pub fn global_inv_parameterize<T: FloatLike>(
    moms: &[ThreeMomentum<F<T>>],
    e_cm_squared: F<T>,
    settings: &Settings,
    force_radius: bool,
) -> (Vec<F<T>>, F<T>) {
    let one = e_cm_squared.one();
    let zero = one.zero();
    if matches!(
        settings.sampling,
        SamplingSettings::DiscreteGraphs(crate::DiscreteGraphSamplingSettings::TropicalSampling(_))
    ) {
        panic!("Trying to inverse parameterize a tropical parametrization.")
    }
    match settings.parameterization.mode {
        ParameterizationMode::HyperSpherical => {
            let e_cm =
                e_cm_squared.sqrt() * F::<T>::from_f64(settings.parameterization.shifts[0].0);
            let mut inv_jac = one.clone();
            let mut xs = Vec::with_capacity(moms.len() * 3);

            let cartesian_xs = moms
                .iter()
                .flat_map(|lv| lv.clone().into_iter())
                .collect::<Vec<F<T>>>();

            let mut k_r_sq = cartesian_xs
                .iter()
                .map(|xi| xi.square())
                .reduce(|acc, e| acc + &e)
                .unwrap_or(zero.clone());
            // cover the degenerate case
            if k_r_sq.is_zero() {
                return (vec![zero.clone(); cartesian_xs.len()], zero);
            }
            let k_r = k_r_sq.sqrt();
            if force_radius {
                xs.push(k_r.clone());
            } else {
                match settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        let b = F::<T>::from_f64(settings.parameterization.b);
                        let x1 = &one - &b / (-&one + &b + (&k_r / &e_cm).exp());
                        inv_jac /= e_cm * &b / (&one - &x1) / (&one + &x1 * (&b - &one));
                        xs.push(x1);
                    }
                    ParameterizationMapping::Linear => {
                        let b = F::<T>::from_f64(settings.parameterization.b);
                        inv_jac /= (&e_cm * &b + &k_r).powi(2) / &e_cm / &b;
                        xs.push(&k_r / (&e_cm * &b + &k_r));
                    }
                }
            };

            let y = cartesian_xs[cartesian_xs.len() - 2].clone();
            let x = cartesian_xs[cartesian_xs.len() - 1].clone();
            let xphi = if x < zero {
                &one + F::<T>::from_f64(0.5) * zero.FRAC_1_PI() * x.atan2(&y)
            } else {
                F::<T>::from_f64(0.5) * zero.FRAC_1_PI() * x.atan2(&y)
            };
            xs.push(xphi);
            inv_jac /= F::<T>::from_f64(2.) * zero.PI();

            for (i, x) in cartesian_xs[..cartesian_xs.len() - 2].iter().enumerate() {
                xs.push(F::<T>::from_f64(0.5) * (&one + x / k_r_sq.sqrt()));
                inv_jac /= F::<T>::from_f64(2.);
                if i > 0 {
                    inv_jac /= (&one - (x * x / &k_r_sq)).sqrt().powi(i as i32);
                }
                k_r_sq -= x * x;
            }

            inv_jac /= k_r.powi((cartesian_xs.len() - 1) as i32);

            let lo = F::<T>::from_f64(settings.parameterization.input_rescaling[0][0].0);
            let hi = F::<T>::from_f64(settings.parameterization.input_rescaling[0][0].1);
            xs[0] = (&xs[0] - &lo) / (&hi - &lo);
            inv_jac /= (&hi - &lo);

            let lo = F::<T>::from_f64(settings.parameterization.input_rescaling[0][1].0);
            let hi = F::<T>::from_f64(settings.parameterization.input_rescaling[0][1].1);
            xs[1] = (&xs[1] - &lo) / (&hi - &lo);
            inv_jac /= &hi - &lo;

            let lo = F::<T>::from_f64(settings.parameterization.input_rescaling[0][2].0);
            let hi = F::<T>::from_f64(settings.parameterization.input_rescaling[0][2].1);
            for x in &mut xs[2..] {
                *x -= &lo / &hi - &lo;
                inv_jac /= (&hi - &lo);
            }
            (xs, inv_jac)
        }
        ParameterizationMode::HyperSphericalFlat => {
            panic!("Inverse of flat hyperspherical sampling is not available since it is not bijective.");
        }
        ParameterizationMode::Cartesian | ParameterizationMode::Spherical => {
            if force_radius {
                panic!("Cannot force radius for non-hyperspherical parameterization.");
            }
            let mut inv_jac = one;
            let mut xs = Vec::with_capacity(moms.len() * 3);
            for (i, mom) in moms.iter().enumerate() {
                let (xs_i, inv_jac_i) = inv_parametrize3d(mom, e_cm_squared.clone(), i, settings);
                xs.extend(xs_i);
                inv_jac *= inv_jac_i;
            }
            (xs, inv_jac)
        }
    }
}

/// Map a vector in the unit hypercube to the infinite hypercube.
/// Also compute the Jacobian.
pub fn parameterize3d<T: FloatLike>(
    x: &[F<T>],
    e_cm_squared: F<T>,
    loop_index: usize,
    settings: &Settings,
) -> ([F<T>; 3], F<T>) {
    let zero = e_cm_squared.zero();
    let one = zero.one();
    let e_cm =
        e_cm_squared.sqrt() * F::<T>::from_f64(settings.parameterization.shifts[loop_index].0);
    let mut l_space = [zero.clone(), zero.clone(), zero.clone()];
    let mut jac = one.clone();

    // rescale the input to the desired range
    let mut x_r = [zero.clone(), zero.clone(), zero.clone()];
    for (xd, xi, &(lo, hi)) in izip!(
        &mut x_r,
        x,
        &settings.parameterization.input_rescaling[loop_index]
    ) {
        let lo = F::<T>::from_f64(lo);
        let hi = F::<T>::from_f64(hi);
        *xd = &lo + xi * (&hi - &lo);
        jac *= &hi - &lo;
    }

    match settings.parameterization.mode {
        ParameterizationMode::Cartesian => match settings.parameterization.mapping {
            ParameterizationMapping::Log => {
                for i in 0..3 {
                    let x = x_r[i].clone();
                    l_space[i] = &e_cm * (&x / (&one - &x)).log();
                    jac *= &e_cm / (&x - &x * &x);
                }
            }
            ParameterizationMapping::Linear => {
                for i in 0..3 {
                    let x = x_r[i].clone();
                    l_space[i] = &e_cm * (&one / (&one - &x) - &one / &x);
                    jac *= &e_cm * (&one / (&x * &x) + &one / ((&one - &x) * (&one - &x)));
                }
            }
        },
        ParameterizationMode::Spherical => {
            let radius = match settings.parameterization.mapping {
                ParameterizationMapping::Log => {
                    // r = &e_cm * ln(1 + b*&x/(1-&x))
                    let x = x_r[0].clone();
                    let b = F::<T>::from_f64(settings.parameterization.b);
                    let radius = &e_cm * (&one + &b * &x / (&one - &x)).log();
                    jac *= &e_cm * &b / (&one - &x) / (&one + &x * (&b - &one));

                    radius
                }
                ParameterizationMapping::Linear => {
                    // r = &e_cm * b * x/(1-x)
                    let b = F::<T>::from_f64(settings.parameterization.b);
                    let radius = &e_cm * &b * &x_r[0] / (&one - &x_r[0]);
                    jac *= (&e_cm * &b + &radius).powi(2) / &e_cm / &b;
                    radius
                }
            };
            let phi = F::<T>::from_f64(2.) * zero.PI() * &x_r[1];
            jac *= F::<T>::from_f64(2.) * zero.PI();

            let cos_theta = -&one + F::<T>::from_f64(2.) * &x_r[2]; // out of range
            jac *= F::<T>::from_f64(2.);
            let sin_theta = (&one - cos_theta.square()).sqrt();

            l_space[0] = &radius * &sin_theta * phi.cos();
            l_space[1] = &radius * &sin_theta * phi.sin();
            l_space[2] = &radius * &cos_theta;

            jac *= radius.square(); // spherical coord
        }
        _ => {
            panic!(
                "Inappropriate parameterization mapping specified for parameterize: {:?}.",
                settings.parameterization.mode.clone()
            );
        }
    }

    // add a shift such that k=l is harder to be picked up by integrators such as cuhre
    l_space[0] += &e_cm * F::<T>::from_f64(settings.parameterization.shifts[loop_index].1);
    l_space[1] += &e_cm * F::<T>::from_f64(settings.parameterization.shifts[loop_index].2);
    l_space[2] += &e_cm * F::<T>::from_f64(settings.parameterization.shifts[loop_index].3);

    (l_space, jac)
}

pub fn inv_parametrize3d<T: FloatLike>(
    mom: &ThreeMomentum<F<T>>,
    e_cm_squared: F<T>,
    loop_index: usize,
    settings: &Settings,
) -> ([F<T>; 3], F<T>) {
    let one = e_cm_squared.one();
    let zero = one.zero();
    if settings.parameterization.mode != ParameterizationMode::Spherical {
        panic!("Inverse mapping is only implemented for spherical coordinates");
    }

    let mut jac = one.clone();
    let e_cm =
        e_cm_squared.sqrt() * F::<T>::from_f64(settings.parameterization.shifts[loop_index].0);

    let x = &mom.px - &e_cm * F::<T>::from_f64(settings.parameterization.shifts[loop_index].1);
    let y = &mom.py - &e_cm * F::<T>::from_f64(settings.parameterization.shifts[loop_index].2);
    let z = &mom.pz - &e_cm * F::<T>::from_f64(settings.parameterization.shifts[loop_index].3);

    let k_r_sq = x.square() + y.square() + z.square();
    let k_r = k_r_sq.sqrt();

    let x2 = if y < zero {
        &one + F::<T>::from_f64(0.5) * zero.FRAC_1_PI() * y.atan2(&x)
    } else {
        F::<T>::from_f64(0.5) * zero.FRAC_1_PI() * y.atan2(&x)
    };

    // cover the degenerate case
    if k_r_sq.is_zero() {
        return ([zero.clone(), x2.clone(), zero.clone()], zero.clone());
    }

    let x1 = match settings.parameterization.mapping {
        ParameterizationMapping::Log => {
            let b = F::<T>::from_f64(settings.parameterization.b);
            let x1 = &one - &b / (-&one + &b + (&k_r / &e_cm).exp());
            jac /= &e_cm * &b / (&one - &x1) / (&one + &x1 * (&b - &one));
            x1
        }
        ParameterizationMapping::Linear => {
            let b = F::<T>::from_f64(settings.parameterization.b);
            jac /= (&e_cm * &b + &k_r).powi(2) / &e_cm / &b;
            &k_r / (&e_cm * &b + &k_r)
        }
    };

    let x3 = F::<T>::from_f64(0.5) * (&one + &z / &k_r);

    jac /= F::<T>::from_f64(2.) * zero.PI();
    jac /= F::<T>::from_f64(2.);
    jac /= k_r.square();

    let mut x = [x1, x2, x3];
    for (xi, &(lo, hi)) in x
        .iter_mut()
        .zip_eq(&settings.parameterization.input_rescaling[loop_index])
    {
        *xi = (xi.clone() - F::<T>::from_f64(lo)) / F::<T>::from_f64(hi - lo);
        jac /= F::<T>::from_f64(hi - lo);
    }

    (x, jac)
}

pub const MINUTE: usize = 60;
pub const HOUR: usize = 3_600;
pub const DAY: usize = 86_400;
pub const WEEK: usize = 604_800;
pub fn format_wdhms(seconds: usize) -> String {
    let mut compound_duration = vec![];
    if seconds == 0 {
        compound_duration.push("0s".to_string());
        return compound_duration.join(" ");
    }

    let mut sec = seconds % WEEK;
    // weeks
    let ws = seconds / WEEK;
    if ws != 0 {
        compound_duration.push(format!("{ws}w"));
    }

    // days
    let ds = sec / DAY;
    sec %= DAY;
    if ds != 0 {
        compound_duration.push(format!("{ds}d"));
    }

    // hours
    let hs = sec / HOUR;
    sec %= HOUR;
    if hs != 0 {
        compound_duration.push(format!("{hs}h"));
    }

    // minutes
    let ms = sec / MINUTE;
    sec %= MINUTE;
    if ms != 0 {
        compound_duration.push(format!("{ms}m"));
    }

    // seconds
    if sec != 0 {
        compound_duration.push(format!("{sec}s"));
    }

    compound_duration.join(" ")
}

pub fn format_wdhms_from_duration(duration: Duration) -> String {
    format_wdhms(duration.as_secs() as usize)
}

#[allow(unused)]
pub fn inverse_gamma_lr(a: f64, p: f64, n_iter: usize) -> f64 {
    // this algorithm is taken from https://dl.acm.org/doi/pdf/10.1145/22721.23109

    // get an estimate for x0 to start newton iterations.
    let q = 1.0 - p;

    if (1.0 - 1.0e-8..=1.0 + 1.0e-8).contains(&a) {
        return -q.ln();
    }

    let gamma_a = gamma(a);
    let b = q * gamma_a;
    let c = 0.577_215_664_901_532_9;

    let mut x0 = 0.5;
    if a < 1.0 {
        if b > 0.6 || (b >= 0.45 && a >= 0.3) {
            let u = if b * q > 10e-8 {
                (p * gamma(a + 1.0)).powf(a.recip())
            } else {
                (-q / a - c).exp()
            };
            x0 = u / (1.0 - u / (a + 1.0));
        } else if a < 0.3 && (0.35..=0.6).contains(&b) {
            let t = (-c - b).exp();
            let u = t * t.exp();
            x0 = t * u.exp();
        } else if (0.15..=0.35).contains(&b) || ((0.15..0.45).contains(&b) && a >= 0.3) {
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            x0 = y - (1.0 - a) * y.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln();
        } else if 0.01 < b && b < 0.15 {
            let y = -b.ln();
            let u = y - (1.0 - a) * y.ln();
            x0 = y
                - (1.0 - a) * u.ln()
                - ((u * u + 2.0 * (3.0 - a) * u + (2.0 - a) * (3.0 - a))
                    / (u * u + (5.0 - a) * u + 2.0))
                    .ln();
        } else if b <= 0.01 {
            let y = -b.ln();
            let c1 = (a - 1.0) * y.ln();
            let c2 = (a - 1.0) * (1.0 + c1);
            let c3 = (a - 1.0) * (-0.5 * c1 * c1 + (a - 2.0) * c1 + (3.0 * a - 5.0) * 0.5);
            let c4 = (a - 1.0)
                * (1.0 / 3.0 * c1 * c1 * c1 - (3.0 * a - 5.0) * 0.5 * c1 * c1
                    + (a * a - 6.0 * a + 7.0) * c1
                    + (11.0 * a * a - 46.0 * a + 47.0) / 6.0);
            let c5 = (a - 1.0)
                * (-0.25 * c1 * c1 * c1 * c1
                    + (11.0 * a - 7.0) / 6.0 * c1 * c1 * c1
                    + (-3.0 * a * a - 13.0) * c1 * c1
                    + (2.0 * a * a * a - 25.0 * a * a + 72.0 * a - 61.0) * 0.5 * c1
                    + (25.0 * a * a * a - 195.0 * a * a + 477.0 * a - 379.0) / 12.0);
            x0 = y + c1 + c2 / (y) + c3 / (y * y) + c4 / (y * y * y) + c5 / (y * y * y * y);

            if b <= 1.0e-28 {
                return x0;
            }
        }
    } else {
        let pref;
        let tau;
        if p < 0.5 {
            pref = -1.0;
            tau = p;
        } else {
            pref = 1.0;
            tau = q;
        }
        let t = (-2.0 * tau.ln()).sqrt();

        let a_0 = 3.31125922108741;
        let a_1 = 11.6616720288968;
        let a_2 = 4.28342155967104;
        let a_3 = 0.213623493715853;

        let b_1 = 6.61053765625462;
        let b_2 = 6.40691597760039;
        let b_3 = 1.27364489782223;
        let b_4 = 3.611_708_101_884_203e-2;

        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;

        let numerator = a_0 + a_1 * t + a_2 * t2 + a_3 * t3;
        let denominator = 1.0 + b_1 * t + b_2 * t2 + b_3 * t3 + b_4 * t4;

        let s = pref * (t - numerator / denominator);
        let s2 = s * s;
        let s3 = s * s2;
        let s4 = s * s3;
        let s5 = s * s4;

        let a_sqrt = a.sqrt();

        let w = a + s * a_sqrt + (s2 - 1.0) / 3.0 + (s3 - 7.0 * s) / (36.0 * a_sqrt)
            - (3.0 * s4 + 7.0 * s2 - 16.0) / (810.0 * a)
            + (9.0 * s5 + 256.0 * s3 - 433.0 * s) / (38880.0 * a * a_sqrt);

        if a >= 500.0 && (1.0 - w / a).abs() < 1.0e-6 {
            return w;
        } else if p > 0.5 {
            if w < 3.0 * a {
                x0 = w;
            } else {
                let d = 2f64.max(a * (a - 1.0));
                if b > 10f64.powf(-d) {
                    let u = -b.ln() + (a - 1.0) * w.ln() - (1.0 + (1.0 - a) / (1.0 + w)).ln();
                    x0 = -b.ln() + (a - 1.0) * u.ln() - (1.0 + (1.0 - a) / (1.0 + u)).ln();
                } else {
                    let y = -b.ln();
                    let c1 = (a - 1.0) * y.ln();
                    let c2 = (a - 1.0) * (1.0 + c1);
                    let c3 = (a - 1.0) * (-0.5 * c1 * c1 + (a - 2.0) * c1 + (3.0 * a - 5.0) * 0.5);
                    let c4 = (a - 1.0)
                        * (1.0 / 3.0 * c1 * c1 * c1 - (3.0 * a - 5.0) * 0.5 * c1 * c1
                            + (a * a - 6.0 * a + 7.0) * c1
                            + (11.0 * a * a - 46.0 * a + 47.0) / 6.0);
                    let c5 = (a - 1.0)
                        * (-0.25 * c1 * c1 * c1 * c1
                            + (11.0 * a - 7.0) / 6.0 * c1 * c1 * c1
                            + (-3.0 * a * a - 13.0) * c1 * c1
                            + (2.0 * a * a * a - 25.0 * a * a + 72.0 * a - 61.0) * 0.5 * c1
                            + (25.0 * a * a * a - 195.0 * a * a + 477.0 * a - 379.0) / 12.0);
                    x0 = y + c1 + c2 / (y) + c3 / (y * y) + c4 / (y * y * y) + c5 / (y * y * y * y);
                }
            }
        } else {
            // this part is heavily simplified from the paper, if any issues occur this estimate
            // will need more refinement.
            let v = (p * gamma(a + 1.0)).ln();
            x0 = ((v + w) / a).exp();
        }
    }

    // start iteration
    let mut x_n = x0;
    for _ in 0..n_iter {
        let r = x_n.powf(a - 1.0) * (-x_n).exp() / gamma_a;
        if x_n <= 0. {
            x_n = 1.0e-16;
        }
        let t_n = if p <= 0.5 {
            (gamma_lr(a, x_n) - p) / r
        } else {
            -(gamma_ur(a, x_n) - q) / r
        };
        let w_n = (a - 1.0 - x_n) / 2.0;

        let h_n = if t_n.abs() <= 0.1 && (w_n * t_n).abs() <= 0.1 {
            t_n + w_n * t_n * t_n
        } else {
            t_n
        };

        x_n -= h_n;
    }

    x_n
}

#[allow(unused)]
pub fn inv_3x3_sig_matrix(mat: [[isize; 3]; 3]) -> [[isize; 3]; 3] {
    let denom = -mat[0][2] * mat[1][1] * mat[2][0]
        + mat[0][1] * mat[1][2] * mat[2][0]
        + mat[0][2] * mat[1][0] * mat[2][1]
        - mat[0][0] * mat[1][2] * mat[2][1]
        - mat[0][1] * mat[1][0] * mat[2][2]
        + mat[0][0] * mat[1][1] * mat[2][2];
    if denom != 1 && denom != -1 {
        panic!("Non invertible signature matrix.");
    }
    let mut inv_mat = [[0; 3]; 3];
    inv_mat[0][0] = (-mat[1][2] * mat[2][1] + mat[1][1] * mat[2][2]) * denom;
    inv_mat[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * denom;
    inv_mat[0][2] = (-mat[0][2] * mat[1][1] + mat[0][1] * mat[1][2]) * denom;
    inv_mat[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * denom;
    inv_mat[1][1] = (-mat[0][2] * mat[2][0] + mat[0][0] * mat[2][2]) * denom;
    inv_mat[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) * denom;
    inv_mat[2][0] = (-mat[1][1] * mat[2][0] + mat[1][0] * mat[2][1]) * denom;
    inv_mat[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) * denom;
    inv_mat[2][2] = (-mat[0][1] * mat[1][0] + mat[0][0] * mat[1][1]) * denom;

    inv_mat
}

pub fn print_banner() {
    info!(
        "\n{}{}\n",
        r"                                        _
                                       | |
   __ _  __ _ _ __ ___  _ __ ___   __ _| |     ___   ___  _ __
  / _` |/ _` | '_ ` _ \| '_ ` _ \ / _` | |    / _ \ / _ \| '_ \
 | (_| | (_| | | | | | | | | | | | (_| | |___| (_) | (_) | |_) |
  \__, |\__,_|_| |_| |_|_| |_| |_|\__,_|______\___/ \___/| .__/
   __/ |                                                 | |
"
        .to_string()
        .bold()
        .blue(),
        format!(
            r#"  |___/    {}                    |_|    "#,
            format!("{:-26}", GIT_VERSION).green(),
        )
        .bold()
        .blue(),
    );
}

#[allow(unused)]
pub fn format_for_compare_digits(x: F<f64>, y: F<f64>) -> (String, String) {
    let mut string_x = format!("{:.16e}", x);
    let mut string_y = format!("{:.16e}", y);

    #[allow(clippy::comparison_chain)]
    if string_x.len() > string_y.len() {
        for _ in 0..(string_x.len() - string_y.len()) {
            string_y.push(' ');
        }
    } else if string_y.len() > string_x.len() {
        for _ in 0..(string_y.len() - string_x.len()) {
            string_x.push(' ');
        }
    }

    let string_vec = string_x
        .chars()
        .zip(string_y.chars())
        .map(|(char_x, char_y)| {
            if char_x == char_y {
                (char_x.to_string().green(), char_y.to_string().green())
            } else {
                (char_x.to_string().red(), char_y.to_string().red())
            }
        })
        .collect_vec();

    let string_x = string_vec.iter().map(|(x, _)| x).join("");
    let string_y = string_vec.iter().map(|(_, y)| y).join("");

    (string_x, string_y)
}

#[allow(unused)]
pub fn format_evaluation_time(time: Duration) -> String {
    let time_secs = time.as_secs_f64();
    if time_secs < 1e-6 {
        format!("{} ns", time.as_nanos())
    } else if time_secs < 1e-3 {
        format!("{:.2} µs", (time.as_nanos() as f64) / 1000.)
    } else if time_secs < 1.0 {
        format!("{:.2} ms", (time.as_micros() as f64) / 1000.)
    } else {
        format!("{:.2} s", (time.as_millis() as f64) / 1000.)
    }
}

pub fn format_evaluation_time_from_f64(time: f64) -> String {
    format_evaluation_time(Duration::from_secs_f64(time))
}

pub fn format_sample(sample: &Sample<F<f64>>) -> String {
    match sample {
        Sample::Continuous(_, xs) => {
            let xs_point = xs.iter().map(|x| format!("{:.16}", x)).join(", ");
            format!("xs: [{}]", xs_point)
        }
        Sample::Discrete(_, graph_index, Some(nested_sample)) => match nested_sample.as_ref() {
            Sample::Continuous(_, xs) => {
                let xs_point = xs.iter().map(|x| format!("{:.16}", x)).join(", ");
                format!("graph: {}, xs: [{}]", graph_index, xs_point)
            }
            Sample::Discrete(_, channel_index, Some(nested_cont_sample)) => {
                match nested_cont_sample.as_ref() {
                    Sample::Continuous(_, xs) => {
                        let xs_point = xs.iter().map(|x| format!("{:.16}", x)).join(", ");
                        format!(
                            "graph: {}, channel: {}, xs: [{}]",
                            graph_index, channel_index, xs_point
                        )
                    }
                    _ => String::from("N/A"),
                }
            }
            _ => String::from("N/A"),
        },
        _ => String::from("N/A"),
    }
}

pub fn view_list_diff_typed<K, T: PartialEq + std::fmt::Debug>(
    vec1: &TiSlice<K, T>,
    vec2: &TiSlice<K, T>,
) -> String {
    let mut result = String::new();

    result.push_str("elements of vec1 that are not in vec2:\n");

    vec1.iter()
        .filter(|vec1_element| !vec2.contains(vec1_element))
        .for_each(|vec1_element| result.push_str(&format!("{:#?}\n", vec1_element)));

    result.push_str("elements of vec2 that are not in vec1:\n");

    vec2.iter()
        .filter(|vec2_element| !vec1.contains(vec2_element))
        .for_each(|vec2_element| result.push_str(&format!("{:#?}", vec2_element)));

    result
}

pub fn into_complex_ff64<T: FloatLike>(c: &Complex<F<T>>) -> Complex<F<f64>> {
    Complex::new(c.re.into_ff64(), c.im.into_ff64())
}

#[test]
fn complex_compare() {
    let ltd = Complex::new(0.11773583919739394, -0.22157450463964778).map(F);

    let cff = Complex::new(0.11773583589023458, -0.22157450446824836).map(F);

    ltd.approx_eq_res(&cff, &F(0.00000001)).unwrap();
}

pub struct GammaloopSymbols {
    pub ubar: Symbol,
    pub vbar: Symbol,
    pub v: Symbol,
    pub u: Symbol,
    pub color_wrap: Symbol,
    pub epsilon: Symbol,
    pub epsilonbar: Symbol,
    pub x_: Symbol,
    pub y_: Symbol,
    pub z_: Symbol,
    pub a_: Symbol,
    pub b_: Symbol,
    pub c_: Symbol,
    pub d_: Symbol,
    pub e_: Symbol,
    pub f_: Symbol,
    pub g_: Symbol,
    pub h_: Symbol,

    pub x__: Symbol,
    pub y__: Symbol,
    pub z__: Symbol,
    pub a__: Symbol,
    pub b__: Symbol,
    pub c__: Symbol,
    pub d__: Symbol,
    pub e__: Symbol,
    pub f__: Symbol,
    pub g__: Symbol,
    pub h__: Symbol,

    pub x___: Symbol,
    pub y___: Symbol,
    pub z___: Symbol,
    pub a___: Symbol,
    pub b___: Symbol,
    pub c___: Symbol,
    pub d___: Symbol,
    pub e___: Symbol,
    pub f___: Symbol,
    pub g___: Symbol,
    pub h___: Symbol,
    pub dim: Symbol,
    pub coeff: Symbol,
}

pub static GS: LazyLock<GammaloopSymbols> = LazyLock::new(|| GammaloopSymbols {
    ubar: symb!("ubar"),
    vbar: symb!("vbar"),
    dim: symb!("dim"),
    v: symb!("v"),
    u: symb!("u"),
    epsilon: symb!("ϵ"),
    color_wrap: symb!("color"),
    epsilonbar: symb!("ϵbar"),
    coeff: symb!("coef"),
    x_: symb!("x_"),
    y_: symb!("y_"),
    z_: symb!("z_"),
    a_: symb!("a_"),
    b_: symb!("b_"),
    c_: symb!("c_"),
    d_: symb!("d_"),
    e_: symb!("e_"),
    f_: symb!("f_"),
    g_: symb!("g_"),
    h_: symb!("h_"),
    x__: symb!("x__"),
    y__: symb!("y__"),
    z__: symb!("z__"),
    a__: symb!("a__"),
    b__: symb!("b__"),
    c__: symb!("c__"),
    d__: symb!("d__"),
    e__: symb!("e__"),
    f__: symb!("f__"),
    g__: symb!("g__"),
    h__: symb!("h__"),
    x___: symb!("x___"),
    y___: symb!("y___"),
    z___: symb!("z___"),
    a___: symb!("a___"),
    b___: symb!("b___"),
    c___: symb!("c___"),
    d___: symb!("d___"),
    e___: symb!("e___"),
    f___: symb!("f___"),
    g___: symb!("g___"),
    h___: symb!("h___"),
});

/// Checks if two lists are permutations of eachother, and establish a map between indices
pub fn is_permutation<T: PartialEq>(left: &[T], right: &[T]) -> Option<PermutationMap> {
    if left.len() != right.len() {
        return None;
    }

    let mut left_to_right = Vec::with_capacity(left.len());
    for elem_in_left in left.iter() {
        let option_position = right
            .iter()
            .enumerate()
            .position(|(right_index, elem_in_right)| {
                elem_in_right == elem_in_left && !left_to_right.contains(&right_index)
            });

        if let Some(position) = option_position {
            left_to_right.push(position);
        } else {
            return None;
        }
    }

    let mut right_to_left = vec![0; left.len()];
    for (index, left_to_right_elem) in left_to_right.iter().enumerate() {
        right_to_left[*left_to_right_elem] = index
    }

    Some(PermutationMap {
        left_to_right,
        right_to_left,
    })
}

#[derive(Clone, Debug)]
pub struct PermutationMap {
    left_to_right: Vec<usize>,
    right_to_left: Vec<usize>,
}

impl PermutationMap {
    pub fn left_to_right(&self, left_index: usize) -> usize {
        self.left_to_right[left_index]
    }

    pub fn right_to_left(&self, right_index: usize) -> usize {
        self.right_to_left[right_index]
    }
}

#[test]
fn test_is_permutation() {
    let a = ["a", "b"];
    let b = ["a", "c"];

    assert!(is_permutation(&a, &b).is_none());

    let a = ["a", "b", "b", "c", "d"];
    let b = ["d", "b", "a", "c", "b"];

    let permutation_map = is_permutation(&a, &b).unwrap();

    for ind in 0..5 {
        assert_eq!(a[ind], b[permutation_map.left_to_right[ind]]);
        assert_eq!(b[ind], a[permutation_map.right_to_left[ind]]);
    }
}

impl<T: FloatLike> momtrop::float::MomTropFloat for F<T> {
    #[inline]
    fn PI(&self) -> Self {
        self.PI()
    }

    #[inline]
    fn abs(&self) -> Self {
        self.abs()
    }

    #[inline]
    fn cos(&self) -> Self {
        <F<T> as Real>::cos(self)
    }

    #[inline]
    fn exp(&self) -> Self {
        <F<T> as Real>::exp(self)
    }

    #[inline]
    fn one(&self) -> Self {
        <F<T> as NumericalFloatLike>::one(self)
    }

    #[inline]
    fn from_f64(&self, value: f64) -> Self {
        F::from_f64(value)
    }

    #[inline]
    fn from_isize(&self, value: isize) -> Self {
        Self(self.0.from_i64(value as i64))
    }

    #[inline]
    fn inv(&self) -> Self {
        <F<T> as NumericalFloatLike>::inv(self)
    }

    #[inline]
    fn ln(&self) -> Self {
        <F<T> as Real>::log(self)
    }

    #[inline]
    fn powf(&self, power: &Self) -> Self {
        <F<T> as Real>::powf(self, power)
    }

    #[inline]
    fn sin(&self) -> Self {
        <F<T> as Real>::sin(self)
    }

    #[inline]
    fn sqrt(&self) -> Self {
        <F<T> as Real>::sqrt(self)
    }

    #[inline]
    fn zero(&self) -> Self {
        <F<T> as NumericalFloatLike>::zero(self)
    }

    #[inline]
    fn to_f64(&self) -> f64 {
        self.into_f64()
    }
}
