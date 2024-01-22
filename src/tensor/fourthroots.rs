use num::{traits::WrappingAdd, Complex};
use std::{
    ops::{Add, Mul},
    process::Output,
};

use arbitrary_int::u2;

use crate::tensor::{UpgradingAdd, UpgradingMul};

use super::SmallestUpgrade;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Fouroot {
    r: u2,
}

impl Fouroot {
    pub fn one() -> Self {
        Fouroot { r: u2::new(0) }
    }
    pub fn i() -> Self {
        Fouroot { r: u2::new(1) }
    }
    pub fn m_one() -> Self {
        Fouroot { r: u2::new(2) }
    }
    pub fn m_i() -> Self {
        Fouroot { r: u2::new(3) }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl Mul for Fouroot {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Fouroot {
            r: self.r.wrapping_add(&rhs.r),
        }
    }
}

impl<'a> SmallestUpgrade<Fouroot> for &'a Fouroot {
    type LCM = &'a Fouroot;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl SmallestUpgrade<&Fouroot> for Fouroot {
    type LCM = Fouroot;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl<T> SmallestUpgrade<Fouroot> for Complex<T> {
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl<T> SmallestUpgrade<&Fouroot> for Complex<T> {
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl<'a, T> SmallestUpgrade<&Fouroot> for &'a Complex<T> {
    type LCM = &'a Complex<T>;
    fn upgrade(self) -> Self::LCM {
        self
    }
}

impl<'a, T> SmallestUpgrade<Fouroot> for &'a Complex<T> {
    type LCM = &'a Complex<T>;
    fn upgrade(self) -> Self::LCM {
        self
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

impl<'a, T> SmallestUpgrade<&Complex<T>> for &'a Fouroot
where
    T: From<i8>,
{
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<T>>::upgrade(*self)
    }
}

impl<'a, T> SmallestUpgrade<Complex<T>> for &'a Fouroot
where
    T: From<i8>,
{
    type LCM = Complex<T>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<T>>::upgrade(*self)
    }
}

impl SmallestUpgrade<f64> for Fouroot {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<f64>>::upgrade(self)
    }
}

impl SmallestUpgrade<Fouroot> for f64 {
    type LCM = Complex<f64>;
    fn upgrade(self) -> Self::LCM {
        SmallestUpgrade::<Complex<f64>>::upgrade(self)
    }
}

impl Add<Fouroot> for Fouroot {
    type Output = Complex<i8>;
    fn add(self, rhs: Fouroot) -> Self::Output {
        SmallestUpgrade::<Complex<i8>>::upgrade(self) + SmallestUpgrade::<Complex<i8>>::upgrade(rhs)
    }
}

impl Add<&Fouroot> for Fouroot {
    type Output = Complex<i8>;
    fn add(self, rhs: &Fouroot) -> Self::Output {
        SmallestUpgrade::<Complex<i8>>::upgrade(self)
            + SmallestUpgrade::<Complex<i8>>::upgrade(*rhs)
    }
}

impl Add<&Fouroot> for &Fouroot {
    type Output = Complex<i8>;
    fn add(self, rhs: &Fouroot) -> Self::Output {
        SmallestUpgrade::<Complex<i8>>::upgrade(*self)
            + SmallestUpgrade::<Complex<i8>>::upgrade(*rhs)
    }
}

impl Add<Fouroot> for &Fouroot {
    type Output = Complex<i8>;
    fn add(self, rhs: Fouroot) -> Self::Output {
        SmallestUpgrade::<Complex<i8>>::upgrade(*self)
            + SmallestUpgrade::<Complex<i8>>::upgrade(rhs)
    }
}

#[test]
fn test_upgrading_mul() {
    let a = Complex::new(1, 0);
    let b = Fouroot::i();
    let c = a.up_mul(b);
    assert_eq!(c, Complex::new(0, 1));

    let af = Complex::new(1.0, 0.0);
    let bf = Fouroot::i();

    let cf = af.up_mul(&bf);

    let aa = bf.up_add(1.);
    assert_eq!(cf, Complex::new(0.0, 1.0));

    let i = Fouroot::i();

    let f = i.up_add(&i);

    let ff = i.up_add(&af);
    assert_eq!(f, Complex::new(0, 2));

    let f = Fouroot::i() * Fouroot::i() * Fouroot::i() * Fouroot::i();

    assert_eq!(f, Fouroot::one());
}
