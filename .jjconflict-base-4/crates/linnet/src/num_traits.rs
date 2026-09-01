//! # Generic Numerical Traits
//!
//! This module provides generic traits for numerical types, specifically for
//! obtaining "zero" and "one" values in a way that can be generic over
//! both owned and borrowed contexts, especially for types that might not
//! directly implement `num_traits::Zero` or `num_traits::One` in all desired forms,
//! or when specific reference-based behavior is needed.
//!
//! ## Traits:
//!
//! - **`RefZero<T = Self>`**:
//!   Provides a method `ref_zero(&self) -> T`. This is designed to return a "zero"
//!   value of type `T`. Implementations are provided for standard numeric primitives
//!   (e.g., `f32`, `i64`) and conditionally for types from the `symbolica` crate
//!   when the `symbolica` feature is enabled. The `Borrow<T>` supertrait allows
//!   flexibility in how the zero value is obtained or represented.
//!
//! - **`RefOne`**:
//!   Provides a method `ref_one(&self) -> Self`. This is designed to return a "one"
//!   value of the implementor's type. Similar to `RefZero`, implementations are
//!   provided for standard numeric primitives and conditionally for `symbolica` types.
//!
//! These traits are primarily used internally within the `linnet` library to abstract
//! over the creation of zero and one values for generic calculations, potentially
//! in contexts involving external libraries like `symbolica`.

use duplicate::duplicate;
use std::{
    borrow::Borrow,
    fmt::Display,
    ops::{Mul, Neg},
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SignError {
    #[error("Invalid value for Sign")]
    InvalidValue,
    #[error("Zero is not a valid value for Sign")]
    ZeroValue,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, PartialOrd, Ord, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(serde_repr::Serialize_repr, serde_repr::Deserialize_repr)
)]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[repr(i8)]
pub enum SignOrZero {
    Zero = 0,
    Plus = 1,
    Minus = -1,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(serde_repr::Serialize_repr, serde_repr::Deserialize_repr)
)]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
#[repr(i8)]
pub enum Sign {
    Positive = 1,
    Negative = -1,
}

impl From<Sign> for SignOrZero {
    fn from(sign: Sign) -> Self {
        match sign {
            Sign::Positive => Self::Plus,
            Sign::Negative => Self::Minus,
        }
    }
}

pub trait Pow<E> {
    fn pow(&self, exponent: E) -> Self;
}

macro_rules! impl_pow_for_sign {
    ($type:ty, $exponent:ty, $positive:expr) => {
        impl Pow<$exponent> for $type {
            fn pow(&self, exponent: $exponent) -> Self {
                match exponent {
                    0 => $positive,
                    exponent => {
                        if exponent % 2 == 0 {
                            $positive
                        } else {
                            *self
                        }
                    }
                }
            }
        }
    };
}

impl_pow_for_sign!(Sign, i32, Sign::Positive);
impl_pow_for_sign!(Sign, u32, Sign::Positive);
impl_pow_for_sign!(Sign, i64, Sign::Positive);
impl_pow_for_sign!(Sign, u64, Sign::Positive);
impl_pow_for_sign!(Sign, isize, Sign::Positive);
impl_pow_for_sign!(Sign, usize, Sign::Positive);
impl_pow_for_sign!(Sign, i128, Sign::Positive);
impl_pow_for_sign!(Sign, u128, Sign::Positive);

impl_pow_for_sign!(SignOrZero, i32, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, u32, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, i64, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, u64, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, isize, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, usize, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, i128, SignOrZero::Plus);
impl_pow_for_sign!(SignOrZero, u128, SignOrZero::Plus);

impl TryFrom<SignOrZero> for Sign {
    type Error = SignError;

    fn try_from(value: SignOrZero) -> Result<Self, Self::Error> {
        match value {
            SignOrZero::Zero => Err(SignError::ZeroValue),
            SignOrZero::Plus => Ok(Sign::Positive),
            SignOrZero::Minus => Ok(Sign::Negative),
        }
    }
}

impl TryFrom<i8> for SignOrZero {
    type Error = SignError;

    fn try_from(value: i8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(SignOrZero::Zero),
            1 => Ok(SignOrZero::Plus),
            -1 => Ok(SignOrZero::Minus),
            _ => Err(SignError::InvalidValue),
        }
    }
}

impl Display for SignOrZero {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SignOrZero::Zero => write!(f, "."),
            SignOrZero::Plus => write!(f, "+"),
            SignOrZero::Minus => write!(f, "-"),
        }
    }
}

impl SignOrZero {
    pub fn is_zero(&self) -> bool {
        matches!(self, SignOrZero::Zero)
    }

    pub fn is_sign(&self) -> bool {
        matches!(self, SignOrZero::Plus | SignOrZero::Minus)
    }

    pub fn is_positive(&self) -> bool {
        matches!(self, SignOrZero::Plus)
    }

    pub fn is_negative(&self) -> bool {
        matches!(self, SignOrZero::Minus)
    }
}

pub trait RefZero<T = Self>: Borrow<T> {
    fn ref_zero(&self) -> T;
}

pub trait RefOne {
    fn ref_one(&self) -> Self;
}

#[cfg(feature = "symbolica")]
impl RefZero for symbolica::atom::Atom {
    fn ref_zero(&self) -> Self {
        symbolica::atom::Atom::num(0)
    }
}

#[cfg(feature = "symbolica")]
impl RefOne for symbolica::atom::Atom {
    fn ref_one(&self) -> Self {
        symbolica::atom::Atom::num(1)
    }
}

#[cfg(feature = "symbolica")]
impl RefZero for symbolica::domains::integer::Integer {
    fn ref_zero(&self) -> Self {
        symbolica::domains::integer::Integer::zero()
    }
}

#[cfg(feature = "symbolica")]
impl RefZero for symbolica::domains::rational::Rational {
    fn ref_zero(&self) -> Self {
        symbolica::domains::rational::Rational::zero()
    }
}

#[cfg(feature = "symbolica")]
impl<T: RefOne + symbolica::domains::float::Real + RefZero> RefOne
    for symbolica::domains::float::Complex<T>
{
    fn ref_one(&self) -> Self {
        Self::new(self.re.ref_one(), self.im.ref_zero())
    }
}

#[cfg(feature = "symbolica")]
impl<T: symbolica::domains::float::Real + RefZero> RefZero
    for symbolica::domains::float::Complex<T>
{
    fn ref_zero(&self) -> Self {
        Self::new(self.re.ref_zero(), self.im.ref_zero())
    }
}

impl RefOne for Sign {
    fn ref_one(&self) -> Self {
        Sign::Positive
    }
}

impl RefZero for SignOrZero {
    fn ref_zero(&self) -> Self {
        SignOrZero::Zero
    }
}

impl RefOne for SignOrZero {
    fn ref_one(&self) -> Self {
        SignOrZero::Plus
    }
}

impl<T: Neg<Output = T>> Mul<T> for Sign {
    type Output = T;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Sign::Positive => rhs,
            Sign::Negative => -rhs,
        }
    }
}

impl Neg for Sign {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Sign::Positive => Sign::Negative,
            Sign::Negative => Sign::Positive,
        }
    }
}

impl<T: Neg<Output = T> + RefZero> Mul<T> for SignOrZero {
    type Output = T;

    fn mul(self, rhs: T) -> Self::Output {
        match self {
            SignOrZero::Plus => rhs,
            SignOrZero::Minus => -rhs,
            SignOrZero::Zero => rhs.ref_zero(),
        }
    }
}

impl Neg for SignOrZero {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            SignOrZero::Plus => SignOrZero::Minus,
            SignOrZero::Minus => SignOrZero::Plus,
            SignOrZero::Zero => SignOrZero::Zero,
        }
    }
}

// impl<T: num::Zero> RefZero for T { future impls Grrr
//     fn zero(&self) -> Self {
//         T::zero()
//     }
// }

duplicate! {
    [types zero one;
        [f32] [0.0] [1.];
        [f64] [0.0] [1.];
        [i8] [0] [1];
        [i16] [0] [1] ;
        [i32] [0] [1];
        [i64] [0] [1];
        [i128] [0] [1];
        [u8] [0] [1];
        [u16] [0] [1];
        [u32] [0] [1];
        [u64] [0] [1];
        [u128] [0] [1];
        ]

    impl RefZero for types{
        fn ref_zero(&self)-> Self{
            zero
        }
    }

    impl RefZero<types> for &types{
        fn ref_zero(&self)-> types{
            zero
        }
    }

    impl RefOne for types{
        fn ref_one(&self)-> Self{
            one
        }
    }
}
