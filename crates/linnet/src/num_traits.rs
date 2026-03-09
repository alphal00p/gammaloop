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
use std::borrow::Borrow;
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
