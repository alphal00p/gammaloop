use duplicate::duplicate;
use num::{traits::WrappingAdd, Complex};
use std::ops::{Add, Mul};

use arbitrary_int::u2;

// use crate::tensor::{UpgradingAdd, UpgradingMul};

use super::SmallestUpgrade;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fouroot {
    r: u2,
}

impl Fouroot {
    pub fn one() -> Self {
        Fouroot { r: u2::new(0) }
    }
    pub fn i() -> Self {
        Fouroot { r: u2::new(1) }
    }
    pub fn n_one() -> Self {
        Fouroot { r: u2::new(2) }
    }
    pub fn n_i() -> Self {
        Fouroot { r: u2::new(3) }
    }
}

duplicate! {
    [
  U T;
  [ Fouroot ]    [Complex<f64>];
  [ &Fouroot]    [Complex<f64>];
  [ Fouroot ]    [Complex<i64>];
  [ &Fouroot]    [Complex<i64>];
  [ Fouroot ]    [Complex<f32>];
  [ &Fouroot]    [Complex<f32>];
  [ Fouroot ]    [Complex<i32>];
  [ &Fouroot]    [Complex<i32>];
]
impl SmallestUpgrade<U> for T {
    type LCM = T;
    fn upgrade(self) -> Self::LCM {
        self
    }
}
impl<'a> SmallestUpgrade<U> for &'a T {
    type LCM = T;
    fn upgrade(self) -> Self::LCM {
        *self
    }
}
}

impl<T> SmallestUpgrade<Complex<T>> for Fouroot
where
    T: From<i8>,
{
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        let one = u2::new(0);
        let i = u2::new(1);
        let minus_one = u2::new(2);
        let minus_i = u2::new(3);

        match self.r {
            x if x == one => Complex::new(T::from(1), T::from(0)),
            x if x == i => Complex::new(T::from(0), T::from(1)),
            x if x == minus_one => Complex::new(T::from(-1), T::from(0)),
            x if x == minus_i => Complex::new(T::from(0), T::from(-1)),
            _ => unreachable!(),
        }
    }
}

impl<T> SmallestUpgrade<&Complex<T>> for Fouroot
where
    T: From<i8>,
{
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<T>>::upgrade(self)
    }
}

duplicate! {
    [
  F U;
  [ Fouroot ]    [Complex<T>];
  [ Fouroot]    [&Complex<T>];
]

impl<'a, T> SmallestUpgrade<U> for &'a Fouroot
where
    T: From<i8>,
{
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<T>>::upgrade(*self)
    }
}
}

// This is why we don't just use Into. And why we need the SmallestUpgrade trait.

impl SmallestUpgrade<f64> for Fouroot {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<f64>>::upgrade(self)
    }
}

impl SmallestUpgrade<Fouroot> for f64 {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        self.into()
    }
}

impl<'a> SmallestUpgrade<Fouroot> for &'a Fouroot {
    type LCM = Fouroot;
    fn upgrade(self) -> Self::LCM {
        *self
    }
}

impl SmallestUpgrade<&Fouroot> for Fouroot {
    type LCM = Fouroot;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

// Arithmetic:

impl Add<Fouroot> for Fouroot {
    type Output = Complex<i8>;
    fn add(self, rhs: Fouroot) -> Self::Output {
        SmallestUpgrade::<Complex<i8>>::upgrade(self) + SmallestUpgrade::<Complex<i8>>::upgrade(rhs)
    }
}

impl<T, Out, M> Add<T> for Fouroot
where
    T: SmallestUpgrade<Fouroot, LCM = M>,
    Fouroot: SmallestUpgrade<T, LCM = M>,

    M: Add<M, Output = Out>,
{
    type Output = Out;
    fn add(self, rhs: T) -> Self::Output {
        self.upgrade() + rhs.upgrade()
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Mul<Fouroot> for Fouroot {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Fouroot {
            r: self.r.wrapping_add(&rhs.r),
        }
    }
}

impl<T, Out, M> Mul<T> for Fouroot
where
    T: SmallestUpgrade<Fouroot, LCM = M>,
    Fouroot: SmallestUpgrade<T, LCM = M>,

    M: Mul<M, Output = Out>,
{
    type Output = Out;
    fn mul(self, rhs: T) -> Self::Output {
        self.upgrade() * rhs.upgrade()
    }
}

duplicate! {
    [
        F ;
        [Complex<f64>];
        [Complex<i32>];
    ]

    impl<Out> Add<Fouroot> for F
    where
        Fouroot: SmallestUpgrade<F>,
        F: SmallestUpgrade<Fouroot>,

        <Fouroot as SmallestUpgrade<F>>::LCM: Add<<Fouroot as SmallestUpgrade<F>>::LCM, Output = Out>,
    {
        type Output = Out;
        fn add(self, rhs: Fouroot) -> Self::Output {
            <F as SmallestUpgrade<Fouroot>>::upgrade(self)
                + <Fouroot as SmallestUpgrade<F>>::upgrade(rhs)
        }
    }

    impl<Out> Add<&Fouroot> for F
    where
        Fouroot: SmallestUpgrade<F>,
        F: SmallestUpgrade<Fouroot>,

        <Fouroot as SmallestUpgrade<F>>::LCM: Add<<Fouroot as SmallestUpgrade<F>>::LCM, Output = Out>,
    {
        type Output = Out;
        fn add(self, rhs: &Fouroot) -> Self::Output {
            <F as SmallestUpgrade<Fouroot>>::upgrade(self)
                + <Fouroot as SmallestUpgrade<F>>::upgrade(*rhs)
        }
    }

    impl<Out> Mul<Fouroot> for F
    where
        Fouroot: SmallestUpgrade<F>,
        F: SmallestUpgrade<Fouroot>,

        <Fouroot as SmallestUpgrade<F>>::LCM: Mul<<Fouroot as SmallestUpgrade<F>>::LCM, Output = Out>,
    {
        type Output = Out;
        fn mul(self, rhs: Fouroot) -> Self::Output {
            <F as SmallestUpgrade<Fouroot>>::upgrade(self)
                * <Fouroot as SmallestUpgrade<F>>::upgrade(rhs)
        }
    }

    impl<Out> Mul<&Fouroot> for F
    where
        Fouroot: SmallestUpgrade<F>,
        F: SmallestUpgrade<Fouroot>,

        <Fouroot as SmallestUpgrade<F>>::LCM: Mul<<Fouroot as SmallestUpgrade<F>>::LCM, Output = Out>,
    {
        type Output = Out;
        fn mul(self, rhs: &Fouroot) -> Self::Output {
            <F as SmallestUpgrade<Fouroot>>::upgrade(self)
                * <Fouroot as SmallestUpgrade<F>>::upgrade(*rhs)
        }
    }
}
#[test]
fn test_upgrading_mul() {
    let a = Complex::new(1, 0);
    let b = Fouroot::i();
    let c = a + b;
    assert_eq!(c, Complex::new(1, 1));

    let af = Complex::new(1.0, 0.0);
    let bf = Fouroot::i();

    let _cf = af + bf;

    let _aa = bf + 1.;
    // assert_eq!(cf, Complex::new(0.0, 1.0));

    let i = Fouroot::i();

    let f = i + i;

    let _ff = i + af;
    assert_eq!(f, Complex::new(0, 2));

    let f = Fouroot::i() * Fouroot::i() * Fouroot::i() * Fouroot::i();

    assert_eq!(f, Fouroot::one());
}
