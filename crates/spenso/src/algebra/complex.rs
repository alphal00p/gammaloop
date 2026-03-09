use std::{
    fmt::{Debug, Display, LowerExp},
    ops::{Add, Neg, Sub},
};

use crate::network::Ref;

use duplicate::duplicate;
use enum_try_as_inner::EnumTryAsInner;
use num::{Float, One, Zero};
#[cfg(feature = "python")]
use pyo3::types::{PyComplex, PyComplexMethods};
use ref_ops::{RefAdd, RefDiv, RefMul, RefSub};
use serde::{Deserialize, Serialize};
#[cfg(feature = "shadowing")]
use symbolica::domains::{
    float::{Complex as SymComplex, Constructible, FloatLike, Real},
    rational::Rational,
};

use crate::algebra::algebraic_traits::{RefOne, RefZero};

pub trait R {}
duplicate! {
    [t;
    [f32];
    [f64];
    [i8];
    [i16];
    [i32];
    [i64];
    [i128];
    [u8];
    [u16];
    [u32];
    [u64];
    [u128];
    ]
    impl R for t {}
}
#[derive(
    Copy,
    Clone,
    PartialEq,
    Serialize,
    Deserialize,
    Hash,
    Eq,
    PartialOrd,
    Ord,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]

pub struct Complex<T> {
    pub re: T,
    pub im: T,
}

#[cfg(feature = "python_stubgen")]
impl pyo3_stub_gen::PyStubType for Complex<f64> {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo::builtin("complex")
    }
}

#[cfg(feature = "python")]
impl<'a> pyo3::FromPyObject<'_, 'a> for Complex<f64> {
    type Error = pyo3::PyErr;

    fn extract(ob: pyo3::Borrowed<'_, 'a, pyo3::PyAny>) -> Result<Self, Self::Error> {
        if let Ok(a) = ob.cast::<PyComplex>() {
            Ok(Complex::new(a.real(), a.imag()))
        } else if let Ok(a) = ob.extract::<f64>() {
            Ok(Complex::new(a, 0.))
        } else {
            Err(pyo3::exceptions::PyValueError::new_err(
                "Not a valid complex number",
            ))
        }
    }
}

impl<T> Complex<T> {
    pub fn as_ref(&self) -> Complex<&T> {
        Complex {
            re: &self.re,
            im: &self.im,
        }
    }
}

impl<T: Clone> From<&Complex<T>> for Complex<T> {
    fn from(value: &Complex<T>) -> Self {
        value.clone()
    }
}

impl<T> Ref for Complex<T> {
    type Ref<'a>
        = &'a Complex<T>
    where
        Self: 'a;
    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

pub mod add;
pub mod add_assign;
pub mod div;
pub mod div_assign;
pub mod mul;
pub mod mul_assign;
pub mod neg;
pub mod sub;
pub mod sub_assign;

#[cfg(feature = "shadowing")]
pub mod symbolica_traits;

#[cfg(feature = "shadowing")]
impl RefZero for Rational {
    fn ref_zero(&self) -> Self {
        self.zero()
    }
}

#[cfg(feature = "shadowing")]
impl From<f64> for Complex<Rational> {
    fn from(re: f64) -> Self {
        Complex {
            re: Rational::try_from(re).unwrap(),
            im: Rational::zero(),
        }
    }
}

#[cfg(feature = "shadowing")]
impl From<Complex<f64>> for Complex<Rational> {
    fn from(value: Complex<f64>) -> Self {
        Complex {
            re: Rational::try_from(value.re).unwrap(),
            im: Rational::try_from(value.im).unwrap(),
        }
    }
}

impl<T: RefZero> RefZero for Complex<T> {
    fn ref_zero(&self) -> Self {
        Complex::new(self.re.ref_zero(), self.im.ref_zero())
    }
}

impl<T: RefZero + RefOne> RefOne for Complex<T> {
    fn ref_one(&self) -> Self {
        Complex::new(self.re.ref_one(), self.im.ref_zero())
    }
}

impl<T> From<T> for Complex<T>
where
    T: RefZero,
{
    fn from(re: T) -> Self {
        Complex {
            im: re.ref_zero(),
            re,
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Real> From<Complex<T>> for SymComplex<T> {
    fn from(complex: Complex<T>) -> Self {
        SymComplex::new(complex.re, complex.im)
    }
}

#[cfg(feature = "shadowing")]
impl<T: Real> From<SymComplex<T>> for Complex<T> {
    fn from(complex: SymComplex<T>) -> Self {
        Complex::new(complex.re, complex.im)
    }
}

impl<T: Default> Default for Complex<T> {
    fn default() -> Self {
        Complex {
            re: T::default(),
            im: T::default(),
        }
    }
}

pub trait SymbolicaComplex {
    type R;
    fn arg(&self) -> Self::R;
}

pub trait NumTraitComplex {
    type R;
    fn arg(&self) -> Self::R;
}

#[cfg(feature = "shadowing")]
impl<T: Real> SymbolicaComplex for Complex<T> {
    type R = T;
    fn arg(&self) -> T {
        self.im.atan2(&self.re)
    }
}

impl<T: Float> NumTraitComplex for Complex<T> {
    type R = T;
    fn arg(&self) -> T {
        self.im.atan2(self.re)
    }
}

impl<T: Zero> Zero for Complex<T> {
    fn zero() -> Self {
        Complex {
            re: T::zero(),
            im: T::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }

    fn set_zero(&mut self) {
        self.re.set_zero();
        self.im.set_zero();
    }
}

impl<T: Zero + One + PartialEq + Sub<Output = T> + Clone> One for Complex<T> {
    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        self.re.is_one() && self.im.is_zero()
    }

    fn one() -> Self {
        Complex {
            re: T::one(),
            im: T::zero(),
        }
    }

    fn set_one(&mut self) {
        self.re.set_one();
        self.im.set_zero();
    }
}

pub trait FloatDerived<T: Float> {
    fn norm(&self) -> T;

    fn to_polar_coordinates(self) -> (T, T);
    fn from_polar_coordinates(r: T, phi: T) -> Self;
}

impl<T: Float + for<'a> RefMul<&'a T, Output = T> + Add<T, Output = T>> FloatDerived<T>
    for Complex<T>
{
    fn norm(&self) -> T {
        self.norm_squared().sqrt()
    }

    #[inline]
    fn to_polar_coordinates(self) -> (T, T)
    where
        T: num::Float,
    {
        (self.norm_squared().sqrt(), self.arg())
    }

    #[inline]
    fn from_polar_coordinates(r: T, phi: T) -> Complex<T> {
        Complex::new(r * phi.cos(), r * phi.sin())
    }
}

impl<T> Complex<T> {
    pub fn map_ref<U>(&self, f: impl Fn(&T) -> U) -> Complex<U> {
        Complex {
            re: f(&self.re),
            im: f(&self.im),
        }
    }

    pub fn map<U>(self, f: impl Fn(T) -> U) -> Complex<U> {
        Complex {
            re: f(self.re),
            im: f(self.im),
        }
    }

    pub fn map_mut(&mut self, mut f: impl FnMut(&mut T)) {
        f(&mut self.re);
        f(&mut self.im);
    }

    pub fn map_ref_mut<U>(&mut self, mut f: impl FnMut(&mut T) -> U) -> Complex<U> {
        Complex {
            re: f(&mut self.re),
            im: f(&mut self.im),
        }
    }
    #[inline]
    pub fn new(re: T, im: T) -> Complex<T> {
        Complex { re, im }
    }

    pub fn new_re(re: T) -> Complex<T>
    where
        T: RefZero,
    {
        Complex {
            im: re.ref_zero(),
            re,
        }
    }

    pub fn new_im(im: T) -> Complex<T>
    where
        T: RefZero,
    {
        Complex {
            re: im.ref_zero(),
            im,
        }
    }

    #[cfg(feature = "shadowing")]
    pub fn new_zero() -> Self
    where
        T: Constructible,
    {
        Complex {
            re: T::new_zero(),
            im: T::new_zero(),
        }
    }

    #[cfg(feature = "shadowing")]
    #[inline]
    pub fn new_i() -> Self
    where
        T: Constructible,
    {
        Complex {
            re: T::new_zero(),
            im: T::new_one(),
        }
    }

    pub fn ref_i(&self) -> Complex<T>
    where
        T: RefOne + RefZero,
    {
        Complex {
            re: self.re.ref_zero(),
            im: self.im.ref_one(),
        }
    }

    pub fn pow<'a>(&'a self, e: u64) -> Self
    where
        T: RefOne
            + RefZero
            + for<'c> RefMul<&'c T, Output = T>
            + for<'c> RefAdd<&'c T, Output = T>
            + for<'c> RefSub<&'c T, Output = T>,
    {
        // TODO: use binary exponentiation
        let mut r = self.ref_one();
        for _ in 0..e {
            r *= self;
        }
        r
    }

    pub fn conj(&self) -> Complex<T>
    where
        T: Clone + Neg<Output = T>,
    {
        Complex::new(self.re.clone(), -self.im.clone())
    }

    #[inline]
    pub fn i() -> Complex<T>
    where
        T: num::Zero + num::One,
    {
        Complex {
            re: T::zero(),
            im: T::one(),
        }
    }

    #[inline]
    pub fn norm_squared(&self) -> T
    where
        T: for<'a> RefMul<&'a T, Output = T> + Add<T, Output = T>,
    {
        (self.re.ref_mul(&self.re)) + (self.im.ref_mul(&self.im))
    }

    pub fn inv(&self) -> Self
    where
        T: for<'a> RefMul<&'a T, Output = T>
            + Add<T, Output = T>
            + for<'a> RefDiv<&'a T, Output = T>
            + Neg<Output = T>,
    {
        let n = self.norm_squared();
        Complex::new(self.re.ref_div(&n), -self.im.ref_div(&n))
    }
}

impl<T: Display> std::fmt::Display for Complex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({}+{}i)", self.re, self.im))
    }
}

impl<T: Debug> std::fmt::Debug for Complex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({:?}+{:?}i)", self.re, self.im))
    }
}

impl<T: LowerExp> std::fmt::LowerExp for Complex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("({:e}+{:e}i)", self.re, self.im))
    }
}

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum RealOrComplexRef<'a, T> {
    Real(&'a T),
    Complex(&'a Complex<T>),
}

#[derive(Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum RealOrComplexMut<'a, T> {
    Real(&'a mut T),
    Complex(&'a mut Complex<T>),
}

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum RealOrComplex<T> {
    Real(T),
    Complex(Complex<T>),
}

impl<T: RefOne> RefOne for RealOrComplex<T> {
    fn ref_one(&self) -> Self {
        match self {
            RealOrComplex::Real(r) => RealOrComplex::Real(r.ref_one()),
            RealOrComplex::Complex(c) => RealOrComplex::Real(c.re.ref_one()),
        }
    }
}

impl<T: RefZero> From<RealOrComplex<T>> for Complex<T> {
    fn from(value: RealOrComplex<T>) -> Self {
        match value {
            RealOrComplex::Real(r) => {
                let z = r.ref_zero();
                Complex::new(r, z)
            }
            RealOrComplex::Complex(c) => c,
        }
    }
}

impl<T> Ref for RealOrComplex<T> {
    type Ref<'a>
        = RealOrComplexRef<'a, T>
    where
        Self: 'a;
    fn refer(&self) -> Self::Ref<'_> {
        match self {
            RealOrComplex::Real(r) => RealOrComplexRef::Real(r),
            RealOrComplex::Complex(c) => RealOrComplexRef::Complex(c),
        }
    }
}

impl<T: Default> RealOrComplex<T> {
    pub fn to_complex(self) -> Complex<T> {
        match self {
            RealOrComplex::Real(r) => Complex::new(r, T::default()),
            RealOrComplex::Complex(c) => c,
        }
    }
}

impl<T: RefZero> RealOrComplex<T> {
    pub fn zero(&self) -> Self {
        match self {
            RealOrComplex::Real(r) => RealOrComplex::Real(r.ref_zero()),
            RealOrComplex::Complex(c) => RealOrComplex::Complex(c.ref_zero()),
        }
    }

    pub fn to_complex_mut(&mut self) {
        if self.is_real() {
            let old = std::mem::replace(self, self.zero());

            if let RealOrComplex::Real(re) = old {
                *self = RealOrComplex::Complex(Complex::new_re(re));
            }
        }
    }
}

impl<T: Display> Display for RealOrComplex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RealOrComplex::Complex(c) => c.fmt(f),
            RealOrComplex::Real(r) => r.fmt(f),
        }
    }
}
