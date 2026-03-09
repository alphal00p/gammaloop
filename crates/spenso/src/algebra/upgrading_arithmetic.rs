use std::borrow::Cow;

use duplicate::duplicate;
use ref_ops::{RefAdd, RefMul, RefSub};
#[cfg(feature = "shadowing")]
use std::ops::{Add, Div, Mul, Sub};

#[cfg(feature = "shadowing")]
use crate::shadowing::symbolica_utils::SerializableAtom;

#[cfg(feature = "shadowing")]
use symbolica::{
    atom::Atom,
    domains::{
        float::{Complex as SymbolicaComplex, Real},
        integer::Integer,
        rational::Rational,
    },
};

use crate::{
    algebra::algebraic_traits::RefZero,
    algebra::complex::{Complex, R},
};
// #[derive(Copy, Clone, PartialEq)]
// pub struct Complex<T: Scalar> {
//     re: T,
//     im: T,
// }

// impl<T: Scalar> Complex<T> {
//     #[inline]
//     pub fn new(re: T, im: T) -> Complex<T> {
//         Complex { re, im }
//     }

//     #[inline]
//     pub fn i() -> Complex<T>
//     where
//         T: num::Zero + num::One,
//     {
//         Complex {
//             re: T::zero(),
//             im: T::one(),
//         }
//     }

//     #[inline]
//     pub fn norm_squared(&self) -> T
//     where
//         for<'a, 'b> &'a T: Mul<&'b T, Output = T>,
//         T: Add<T, Output = T>,
// #[derive(Copy, Clone, PartialEq)]
// pub struct Complex<T: Scalar> {
//     re: T,
//     im: T,
// }

// impl<T: Scalar> Complex<T> {
//     #[inline]
//     pub fn new(re: T, im: T) -> Complex<T> {
//         Complex { re, im }
//     }

//     #[inline]
//     pub fn i() -> Complex<T>
//     where
//         T: num::Zero + num::One,
//     {
//         Complex {
//             re: T::zero(),
//             im: T::one(),
//         }
//     }

//     #[inline]
//     pub fn norm_squared(&self) -> T
//     where
//         for<'a, 'b> &'a T: Mul<&'b T, Output = T>,
//         T: Add<T, Output = T>,
//     {
//         &self.re * &self.re + &self.im * &self.im
//     }

//     #[inline]
//     pub fn arg(&self) -> T
//     where
//         T: num::Float,
//     {
//         self.im.atan2(self.re)
//     }

//     #[inline]
//     pub fn to_polar_coordinates(self) -> (T, T)
//     where
//         T: num::Float,
//         for<'a, 'b> &'a T: Mul<&'b T, Output = T>,
//         T: Add<T, Output = T>,
//     {
//         (self.norm_squared().sqrt(), self.arg())
//     }

//     #[inline]
//     pub fn from_polar_coordinates(r: T, phi: T) -> Complex<T>
//     where
//         T: num::Float,
//     {
//         Complex::new(r * phi.cos(), r * phi.sin())
//     }
// }

// impl<'a, 'b, T: Scalar> Add<&'a Complex<T>> for &'b Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn add(self, rhs: &'a Complex<T>) -> Self::Output {
//         Complex::new(&self.re + &rhs.re, &self.im + &rhs.im)
//     }
// }

// impl<'a, T: Scalar> Add<&'a Complex<T>> for Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn add(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self + rhs
//     }
// }

// impl<'a, T: Scalar> Add<Complex<T>> for &'a Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn add(self, rhs: Complex<T>) -> Self::Output {
//         self + &rhs
//     }
// }

// impl<T: Scalar> Add for Complex<T> {
//     type Output = Self;

//     #[inline]
//     fn add(self, rhs: Self) -> Self::Output {
//         &self + &rhs
//     }
// }

// impl<T: Scalar> AddAssign for Complex<T>
// where
//     for<'a> T: AddAssign<&'a T>,
// {
//     #[inline]
//     fn add_assign(&mut self, rhs: Self) {
//         self.add_assign(&rhs)
//     }
// }

// impl<T: Scalar> AddAssign<&Complex<T>> for Complex<T>
// where
//     for<'a> T: AddAssign<&'a T>,
// {
//     #[inline]
//     fn add_assign(&mut self, rhs: &Self) {
//         self.re += &rhs.re;
//         self.im += &rhs.im;
//     }
// }

// impl<'a, 'b, T: Scalar> Sub<&'a Complex<T>> for &'b Complex<T>
// where
//     for<'c, 'd> &'c T: Sub<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn sub(self, rhs: &'a Complex<T>) -> Self::Output {
//         Complex::new(&self.re - &rhs.re, &self.im - &rhs.im)
//     }
// }

// impl<'a, T: Scalar> Sub<&'a Complex<T>> for Complex<T>
// where
//     for<'c, 'd> &'c T: Sub<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn sub(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self - rhs
//     }
// }

// impl<'a, T: Scalar> Sub<Complex<T>> for &'a Complex<T>
// where
//     for<'c, 'd> &'c T: Sub<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn sub(self, rhs: Complex<T>) -> Self::Output {
//         self - &rhs
//     }
// }

// impl<T: Scalar> Sub for Complex<T>
// where
//     for<'a, 'b> &'a T: Sub<&'b T, Output = T>,
// {
//     type Output = Self;

//     #[inline]
//     fn sub(self, rhs: Self) -> Self::Output {
//         &self - &rhs
//     }
// }

// impl<T: Scalar> SubAssign for Complex<T>
// where
//     for<'a> T: SubAssign<&'a T>,
// {
//     #[inline]
//     fn sub_assign(&mut self, rhs: Self) {
//         self.sub_assign(&rhs)
//     }
// }

// impl<T: Scalar> SubAssign<&Complex<T>> for Complex<T>
// where
//     for<'a> T: SubAssign<&'a T>,
// {
//     #[inline]
//     fn sub_assign(&mut self, rhs: &Self) {
//         self.re -= &rhs.re;
//         self.im -= &rhs.im;
//     }
// }

// impl<'a, 'b, T: Scalar> Mul<&'a Complex<T>> for &'b Complex<T>
// where
//     for<'c, 'd> &'c T: Mul<&'d T, Output = T>,
//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn mul(self, rhs: &'a Complex<T>) -> Self::Output {
//         Complex::new(
//             &self.re * &rhs.re - &self.im * &rhs.im,
//             &self.re * &rhs.im + &self.im * &rhs.re,
//         )
//     }
// }

// impl<'a, T: Scalar> Mul<&'a Complex<T>> for Complex<T>
// where
//     for<'c, 'd> &'c T: Mul<&'d T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn mul(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self * rhs
//     }
// }

// impl<'a, T: Scalar> Mul<Complex<T>> for &'a Complex<T>
// where
//     for<'c, 'd> &'c T: Mul<&'d T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn mul(self, rhs: Complex<T>) -> Self::Output {
//         self * &rhs
//     }
// }

// impl<T: Scalar> Mul for Complex<T>
// where
//     for<'a, 'b> &'a T: Mul<&'b T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Self;

//     #[inline]
//     fn mul(self, rhs: Self) -> Self::Output {
//         &self * &rhs
//     }
// }

// impl<T: Scalar> MulAssign for Complex<T>
// where
//     for<'a> T: MulAssign<&'a T>,
// {
//     #[inline]
//     fn mul_assign(&mut self, rhs: Self) {
//         self.mul_assign(rhs)
//     }
// }

// impl<T: Scalar> MulAssign<&Complex<T>> for Complex<T>
// where
//     for<'a, 'b> &'a T: Mul<&'b T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     #[inline]
//     fn mul_assign(&mut self, rhs: &Self) {
//         let res = Mul::mul(&*self, rhs);
//         *self = res;
//     }
// }

// impl<'a, 'b, T: Scalar> Div<&'a Complex<T>> for &'b Complex<T>
// where
//     for<'c, 'd> &'c T: Div<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn div(self, rhs: &'a Complex<T>) -> Self::Output {
//         let n = rhs.norm_squared();
//         let re = self.re * rhs.re + self.im * rhs.im;
//         let im = self.im * rhs.re - self.re * rhs.im;
//         Complex::new(re / n, im / n)
//     }
// }

// impl<'a, T: Scalar> Div<&'a Complex<T>> for Complex<T>
// where
//     for<'c, 'd> &'c T: Div<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn div(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self / rhs
//     }
// }

// impl<'a, T: Scalar> Div<Complex<T>> for &'a Complex<T>
// where
//     for<'c, 'd> &'c T: Div<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn div(self, rhs: Complex<T>) -> Self::Output {
//         self / &rhs
//     }
// }

// impl<T: Scalar> Div for Complex<T>
// where
//     for<'a, 'b> &'a T: Div<&'b T, Output = T>,
// {
//     type Output = Self;

//     #[inline]
//     fn div(self, rhs: Self) -> Self::Output {
//         &self / &rhs
//     }
// }

// impl<T: Scalar> DivAssign for Complex<T>
// where
//     for<'a> T: DivAssign<&'a T>,
// {
//     #[inline]
//     fn div_assign(&mut self, rhs: Self) {
//         self.div_assign(&rhs)
//     }
// }

// impl<T: Scalar> DivAssign<&Complex<T>> for Complex<T>
// where
//     for<'a> T: DivAssign<&'a T>,
// {
//     #[inline]
//     fn div_assign(&mut self, rhs: &Self) {
//         self.re /= &rhs.re;
//         self.im /= &rhs.im;
//     }
// }

// impl<T: Scalar> Neg for Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn neg(self) -> Complex<T> {
//         Complex::new(-self.re, -self.im)
//     }
// }

// impl<T: Scalar> std::fmt::Display for Complex<T> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         f.write_fmt(format_args!("({}+{}i)", self.re, self.im))
//     }
// }

// impl<T: Scalar> std::fmt::Debug for Complex<T> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         f.write_fmt(format_args!("({:?}+{:?}i)", self.re, self.im))
//     }
// }
//         Complex::new(r * phi.cos(), r * phi.sin())
//     }
// }

// impl<'a, 'b, T: Scalar> Add<&'a Complex<T>> for &'b Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn add(self, rhs: &'a Complex<T>) -> Self::Output {
//         Complex::new(&self.re + &rhs.re, &self.im + &rhs.im)
//     }
// }

// impl<'a, T: Scalar> Add<&'a Complex<T>> for Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn add(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self + rhs
//     }
// }

// impl<'a, T: Scalar> Add<Complex<T>> for &'a Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn add(self, rhs: Complex<T>) -> Self::Output {
//         self + &rhs
//     }
// }

// impl<T: Scalar> Add for Complex<T> {
//     type Output = Self;

//     #[inline]
//     fn add(self, rhs: Self) -> Self::Output {
//         &self + &rhs
//     }
// }

// impl<T: Scalar> AddAssign for Complex<T>
// where
//     for<'a> T: AddAssign<&'a T>,
// {
//     #[inline]
//     fn add_assign(&mut self, rhs: Self) {
//         self.add_assign(&rhs)
//     }
// }

// impl<T: Scalar> AddAssign<&Complex<T>> for Complex<T>
// where
//     for<'a> T: AddAssign<&'a T>,
// {
//     #[inline]
//     fn add_assign(&mut self, rhs: &Self) {
//         self.re += &rhs.re;
//         self.im += &rhs.im;
//     }
// }

// impl<'a, 'b, T: Scalar> Sub<&'a Complex<T>> for &'b Complex<T>
// where
//     for<'c, 'd> &'c T: Sub<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn sub(self, rhs: &'a Complex<T>) -> Self::Output {
//         Complex::new(&self.re - &rhs.re, &self.im - &rhs.im)
//     }
// }

// impl<'a, T: Scalar> Sub<&'a Complex<T>> for Complex<T>
// where
//     for<'c, 'd> &'c T: Sub<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn sub(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self - rhs
//     }
// }

// impl<'a, T: Scalar> Sub<Complex<T>> for &'a Complex<T>
// where
//     for<'c, 'd> &'c T: Sub<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn sub(self, rhs: Complex<T>) -> Self::Output {
//         self - &rhs
//     }
// }

// impl<T: Scalar> Sub for Complex<T>
// where
//     for<'a, 'b> &'a T: Sub<&'b T, Output = T>,
// {
//     type Output = Self;

//     #[inline]
//     fn sub(self, rhs: Self) -> Self::Output {
//         &self - &rhs
//     }
// }

// impl<T: Scalar> SubAssign for Complex<T>
// where
//     for<'a> T: SubAssign<&'a T>,
// {
//     #[inline]
//     fn sub_assign(&mut self, rhs: Self) {
//         self.sub_assign(&rhs)
//     }
// }

// impl<T: Scalar> SubAssign<&Complex<T>> for Complex<T>
// where
//     for<'a> T: SubAssign<&'a T>,
// {
//     #[inline]
//     fn sub_assign(&mut self, rhs: &Self) {
//         self.re -= &rhs.re;
//         self.im -= &rhs.im;
//     }
// }

// impl<'a, 'b, T: Scalar> Mul<&'a Complex<T>> for &'b Complex<T>
// where
//     for<'c, 'd> &'c T: Mul<&'d T, Output = T>,
//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn mul(self, rhs: &'a Complex<T>) -> Self::Output {
//         Complex::new(
//             &self.re * &rhs.re - &self.im * &rhs.im,
//             &self.re * &rhs.im + &self.im * &rhs.re,
//         )
//     }
// }

// impl<'a, T: Scalar> Mul<&'a Complex<T>> for Complex<T>
// where
//     for<'c, 'd> &'c T: Mul<&'d T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn mul(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self * rhs
//     }
// }

// impl<'a, T: Scalar> Mul<Complex<T>> for &'a Complex<T>
// where
//     for<'c, 'd> &'c T: Mul<&'d T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn mul(self, rhs: Complex<T>) -> Self::Output {
//         self * &rhs
//     }
// }

// impl<T: Scalar> Mul for Complex<T>
// where
//     for<'a, 'b> &'a T: Mul<&'b T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     type Output = Self;

//     #[inline]
//     fn mul(self, rhs: Self) -> Self::Output {
//         &self * &rhs
//     }
// }

// impl<T: Scalar> MulAssign for Complex<T>
// where
//     for<'a> T: MulAssign<&'a T>,
// {
//     #[inline]
//     fn mul_assign(&mut self, rhs: Self) {
//         self.mul_assign(rhs)
//     }
// }

// impl<T: Scalar> MulAssign<&Complex<T>> for Complex<T>
// where
//     for<'a, 'b> &'a T: Mul<&'b T, Output = T>,

//     T: Sub<T, Output = T> + Add<T, Output = T>,
// {
//     #[inline]
//     fn mul_assign(&mut self, rhs: &Self) {
//         let res = Mul::mul(&*self, rhs);
//         *self = res;
//     }
// }

// impl<'a, 'b, T: Scalar> Div<&'a Complex<T>> for &'b Complex<T>
// where
//     for<'c, 'd> &'c T: Div<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn div(self, rhs: &'a Complex<T>) -> Self::Output {
//         let n = rhs.norm_squared();
//         let re = self.re * rhs.re + self.im * rhs.im;
//         let im = self.im * rhs.re - self.re * rhs.im;
//         Complex::new(re / n, im / n)
//     }
// }

// impl<'a, T: Scalar> Div<&'a Complex<T>> for Complex<T>
// where
//     for<'c, 'd> &'c T: Div<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn div(self, rhs: &'a Complex<T>) -> Self::Output {
//         &self / rhs
//     }
// }

// impl<'a, T: Scalar> Div<Complex<T>> for &'a Complex<T>
// where
//     for<'c, 'd> &'c T: Div<&'d T, Output = T>,
// {
//     type Output = Complex<T>;

//     #[inline]
//     fn div(self, rhs: Complex<T>) -> Self::Output {
//         self / &rhs
//     }
// }

// impl<T: Scalar> Div for Complex<T>
// where
//     for<'a, 'b> &'a T: Div<&'b T, Output = T>,
// {
//     type Output = Self;

//     #[inline]
//     fn div(self, rhs: Self) -> Self::Output {
//         &self / &rhs
//     }
// }

// impl<T: Scalar> DivAssign for Complex<T>
// where
//     for<'a> T: DivAssign<&'a T>,
// {
//     #[inline]
//     fn div_assign(&mut self, rhs: Self) {
//         self.div_assign(&rhs)
//     }
// }

// impl<T: Scalar> DivAssign<&Complex<T>> for Complex<T>
// where
//     for<'a> T: DivAssign<&'a T>,
// {
//     #[inline]
//     fn div_assign(&mut self, rhs: &Self) {
//         self.re /= &rhs.re;
//         self.im /= &rhs.im;
//     }
// }

// impl<T: Scalar> Neg for Complex<T> {
//     type Output = Complex<T>;

//     #[inline]
//     fn neg(self) -> Complex<T> {
//         Complex::new(-self.re, -self.im)
//     }
// }

// impl<T: Scalar> std::fmt::Display for Complex<T> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         f.write_fmt(format_args!("({}+{}i)", self.re, self.im))
//     }
// }

// impl<T: Scalar> std::fmt::Debug for Complex<T> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         f.write_fmt(format_args!("({:?}+{:?}i)", self.re, self.im))
//     }
// }

pub trait SmallestUpgrade<T> {
    type LCM;
    type Order;
    fn upgrade(self) -> Self::LCM;
}

pub struct LessThan;
pub struct GreaterThan;
pub struct Equal;

pub trait TryFromUpgrade<T> {
    fn try_from_upgrade(value: &T) -> Option<Self>
    where
        Self: Sized;
}

pub trait TrySmallestUpgrade<T> {
    type LCM: Clone;
    // type Order;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>;
}

// pub trait SmallerThan<T> {}
// pub trait LargerThan<T> {}

// pub trait Incomparable<T> {}

// impl<T, U> LargerThan<U> for T where U: SmallerThan<T> {}

// impl<T, U> TrySmallestUpgrade<T> for U
// where
//     T: TrySmallestUpgrade<U, LCM = U, Order = LessThan>,
//     U: Clone,
// {
//     type LCM = U;
//

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>> {
//         Some(Cow::Borrowed(self))
//     }
// }

// impl<U> TrySmallestUpgrade<U> for U {
//     type LCM = U;
//

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>> {
//         Some(Cow::Borrowed(self))
//     }
// }

#[cfg(feature = "shadowing")]
impl<O: Real, U: Real, T: Real> TrySmallestUpgrade<SymbolicaComplex<T>> for SymbolicaComplex<U>
where
    U: TrySmallestUpgrade<T, LCM = O>,
    T: TrySmallestUpgrade<U, LCM = O>,
{
    type LCM = SymbolicaComplex<O>;
    // type Order = LessThan;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        let re = self.re.try_upgrade()?;
        let im = self.im.try_upgrade()?;
        Some(Cow::Owned(SymbolicaComplex::new(
            re.into_owned(),
            im.into_owned(),
        )))
    }
}

impl<O: R + Clone, U: R, T: R> TrySmallestUpgrade<Complex<T>> for Complex<U>
where
    U: TrySmallestUpgrade<T, LCM = O>,
    T: TrySmallestUpgrade<U, LCM = O>,
{
    type LCM = Complex<O>;
    // type Order = LessThan;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        let re = self.re.try_upgrade()?;
        let im = self.im.try_upgrade()?;
        Some(Cow::Owned(Complex::new(re.into_owned(), im.into_owned())))
    }
}

impl<T, U> TryFromUpgrade<T> for U
where
    T: TrySmallestUpgrade<U, LCM = U>,
    U: Clone,
{
    fn try_from_upgrade(value: &T) -> Option<Self> {
        let cow = value.try_upgrade()?;
        Some(cow.into_owned())
    }
}

pub trait TryIntoUpgrade<T> {
    fn try_into_upgrade(&self) -> Option<T>;
}

impl<T, U> TryIntoUpgrade<U> for T
where
    U: TryFromUpgrade<T>,
{
    fn try_into_upgrade(&self) -> Option<U> {
        U::try_from_upgrade(self)
    }
}

// duplicate! {
//     [ num;
//     [f64] ;
//     [i32] ;]

// impl TrySmallestUpgrade<num> for Complex<num> {
//     type LCM = Complex<num>;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Borrowed(self))
//     }
// }

// impl TrySmallestUpgrade<Complex<num>> for num {
//     type LCM = Complex<num>;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Owned(Complex::from(*self)))
//     }
// }
// }

// impl<T, U> TrySmallestUpgrade<U> for T
// where
//     T: Borrow<T>,
//     U: Borrow<T>,
// {
//     type LCM = T;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//     where
//         Self::LCM: Clone,
//     {
//         Some(Cow::Borrowed(self))
//     }
// } can't do this because of future impls GRR.

// duplicate! {
//     [smaller larger;
//     [f64] [Atom];
//     [Atom] [Atom];
//     [Complex<f64>] [Atom];
//     [i32] [Atom];]
// impl TrySmallestUpgrade<smaller> for larger {
//     type LCM = larger;
//
//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Borrowed(self))
//     }
// }

// impl<'a> TrySmallestUpgrade<&'a smaller> for larger {
//     type LCM = larger;
//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Borrowed(self))
//     }
// }

// impl<'a,'b> TrySmallestUpgrade<&'a smaller> for &'b larger {
//     type LCM = larger;
//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Borrowed(*self))
//     }
// }

// impl<'b> TrySmallestUpgrade<smaller> for &'b larger {
//     type LCM = larger;
//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
//         where
//             Self::LCM: Clone {
//         Some(Cow::Borrowed(*self))
//     }
// }

// }

duplicate! {
    [equal;
    [i8];
    [i16];
    [i32];
    [f64];
    // [Complex<f64>];
    ]

impl TrySmallestUpgrade<equal> for equal {
    type LCM = equal;

    fn try_upgrade(&self) -> Option<Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}

}
#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Atom> for Atom {
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        Some(Cow::Borrowed(self))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Integer> for Integer {
    type LCM = Integer;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        Some(Cow::Borrowed(self))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Complex<Integer>> for Complex<Integer> {
    type LCM = Complex<Integer>;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        Some(Cow::Borrowed(self))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Rational> for Rational {
    type LCM = Rational;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        Some(Cow::Borrowed(self))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Complex<Rational>> for Complex<Rational> {
    type LCM = Complex<Rational>;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        Some(Cow::Borrowed(self))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<SerializableAtom> for SerializableAtom {
    type LCM = SerializableAtom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>>
    where
        Self::LCM: Clone,
    {
        Some(Cow::Borrowed(self))
    }
}

duplicate! {
    [smaller larger;
    [f32] [f64];
    [i32] [f64];]

impl TrySmallestUpgrade<larger> for smaller {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Owned(larger::from(*self)))
    }
}

impl TrySmallestUpgrade<smaller> for larger {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}

}

// impl<T> TrySmallestUpgrade<Complex<T>> for T
// where
//     T: RefZero + Clone,
// {
//     type LCM = Complex<T>;

//     fn try_upgrade(&self) -> Option<Cow<Self::LCM>> {
//         let new = Complex::new(self.clone(), self.ref_zero());
//         Some(Cow::Owned(new))
//     }
// }

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Atom> for f64 {
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>> {
        let symrat = Atom::num(*self);

        Some(Cow::Owned(symrat))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<SerializableAtom> for f64 {
    type LCM = SerializableAtom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>> {
        let symrat = Atom::num(*self);

        Some(Cow::Owned(symrat.into()))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<Atom> for i32 {
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>> {
        let symnum = Atom::num(*self);

        Some(Cow::Owned(symnum))
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<SerializableAtom> for i32 {
    type LCM = SerializableAtom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>> {
        let symnum = Atom::num(*self).into();

        Some(Cow::Owned(symnum))
    }
}

#[cfg(feature = "shadowing")]
impl Add for SerializableAtom {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        SerializableAtom::from(self.0 + rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Add<&SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn add(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 + &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Add<&SerializableAtom> for SerializableAtom {
    type Output = SerializableAtom;
    fn add(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 + &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Add<SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn add(self, rhs: SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 + &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Sub for SerializableAtom {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        SerializableAtom::from(self.0 - rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Sub<&SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn sub(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 - &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Sub<&SerializableAtom> for SerializableAtom {
    type Output = SerializableAtom;
    fn sub(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 - &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Sub<SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn sub(self, rhs: SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 - &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Mul for SerializableAtom {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        SerializableAtom::from(self.0 * rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Mul<&SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn mul(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 * &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Mul<&SerializableAtom> for SerializableAtom {
    type Output = SerializableAtom;
    fn mul(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 * &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Mul<SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn mul(self, rhs: SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 * &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Div for SerializableAtom {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        SerializableAtom::from(self.0 / rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Div<&SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn div(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 / &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Div<&SerializableAtom> for SerializableAtom {
    type Output = SerializableAtom;
    fn div(self, rhs: &SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 / &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl Div<SerializableAtom> for &SerializableAtom {
    type Output = SerializableAtom;
    fn div(self, rhs: SerializableAtom) -> Self::Output {
        SerializableAtom::from(&self.0 / &rhs.0)
    }
}

#[cfg(feature = "shadowing")]
impl TrySmallestUpgrade<SerializableAtom> for Complex<f64> {
    type LCM = SerializableAtom;

    fn try_upgrade(&self) -> Option<Cow<'_, Self::LCM>> {
        let real: Cow<'_, SerializableAtom> =
            <f64 as TrySmallestUpgrade<SerializableAtom>>::try_upgrade(&self.re)?;
        let imag: Cow<'_, SerializableAtom> =
            <f64 as TrySmallestUpgrade<SerializableAtom>>::try_upgrade(&self.im)?;

        let i = SerializableAtom::from(Atom::i());
        let symrat = (i * imag.as_ref()) + real.as_ref();

        Some(Cow::Owned(symrat))
    }
}

#[cfg(feature = "shadowing")]
duplicate! {
    [smaller larger;
[f64] [Atom];
[f64] [SerializableAtom];
[i32] [Atom];
[i32] [SerializableAtom];
// [Complex<f64>] [Atom];
[Complex<f64>] [SerializableAtom];
]

impl TrySmallestUpgrade<smaller> for larger {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}


}

duplicate! {
    [smaller larger;

[f64] [Complex<f64>];

]

impl TrySmallestUpgrade<smaller> for larger {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
        Some(Cow::Borrowed(self))
    }
}

impl TrySmallestUpgrade<larger> for smaller {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
               let z = self.ref_zero();
                Some(Cow::Owned(Complex::new(*self, z)))
    }
}


}

pub trait FallibleMul<T = Self> {
    type Output;
    fn mul_fallible(&self, rhs: &T) -> Option<Self::Output>;
}

impl<T, U> FallibleMul<T> for U
where
    U: TrySmallestUpgrade<T>,
    T: TrySmallestUpgrade<U, LCM = <U as TrySmallestUpgrade<T>>::LCM>,
    U::LCM: Clone,
    U::LCM: for<'b> RefMul<&'b U::LCM, Output = U::LCM>,
{
    type Output = U::LCM;

    fn mul_fallible(&self, rhs: &T) -> Option<Self::Output> {
        let lhs = self.try_upgrade()?;
        let rhs = rhs.try_upgrade()?;
        Some(lhs.as_ref().ref_mul(rhs.as_ref()))
    }
}

pub trait FallibleAdd<T> {
    type Output;
    fn add_fallible(&self, rhs: &T) -> Option<Self::Output>;
}

impl<T, U> FallibleAdd<T> for U
where
    U: TrySmallestUpgrade<T>,
    T: TrySmallestUpgrade<U, LCM = U::LCM>,
    U::LCM: for<'b> RefAdd<&'b U::LCM, Output = U::LCM>,
{
    type Output = U::LCM;

    fn add_fallible(&self, rhs: &T) -> Option<Self::Output> {
        let lhs = self.try_upgrade()?;
        let rhs = rhs.try_upgrade()?;
        Some(lhs.as_ref().ref_add(rhs.as_ref()))
    }
}

pub trait FallibleSub<T> {
    type Output;
    fn sub_fallible(&self, rhs: &T) -> Option<Self::Output>;
}

impl<T, U> FallibleSub<T> for U
where
    U: TrySmallestUpgrade<T>,
    T: TrySmallestUpgrade<U, LCM = <U as TrySmallestUpgrade<T>>::LCM>,
    U::LCM: Clone,
    U::LCM: for<'b> RefSub<&'b U::LCM, Output = U::LCM>,
{
    type Output = U::LCM;

    fn sub_fallible(&self, rhs: &T) -> Option<Self::Output> {
        let lhs = self.try_upgrade()?;
        let rhs = rhs.try_upgrade()?;
        Some(lhs.as_ref().ref_sub(rhs.as_ref()))
    }
}

pub trait FallibleAddAssign<T = Self> {
    fn add_assign_fallible(&mut self, rhs: &T);
}

impl<T, U> FallibleAddAssign<T> for U
where
    U: TrySmallestUpgrade<T, LCM = U>,
    T: TrySmallestUpgrade<U, LCM = U>,
    U::LCM: Clone,
    U::LCM: for<'b> RefAdd<&'b U::LCM, Output = U::LCM>,
{
    fn add_assign_fallible(&mut self, rhs: &T) {
        let lhs = self.try_upgrade().unwrap();
        let rhs = rhs.try_upgrade().unwrap();
        let out = lhs.as_ref().ref_add(rhs.as_ref());
        *self = out;
    }
}

pub trait FallibleSubAssign<T> {
    fn sub_assign_fallible(&mut self, rhs: &T);
}

impl<T, U> FallibleSubAssign<T> for U
where
    U: TrySmallestUpgrade<T, LCM = U>,
    T: TrySmallestUpgrade<U, LCM = U>,
    U::LCM: Clone,
    U::LCM: for<'b> RefSub<&'b U::LCM, Output = U::LCM>,
{
    fn sub_assign_fallible(&mut self, rhs: &T) {
        let lhs = self.try_upgrade().unwrap();
        let rhs = rhs.try_upgrade().unwrap();
        let out = lhs.as_ref().ref_sub(rhs.as_ref());
        *self = out;
    }
}

#[cfg(test)]
mod test {
    #[cfg(feature = "shadowing")]
    use ahash::{HashMap, HashMapExt};

    #[cfg(feature = "shadowing")]
    use symbolica::{atom::Atom, parse};

    #[cfg(feature = "shadowing")]
    use crate::algebra::complex::Complex;

    use super::*;

    #[test]
    fn i32_arithmetic() {
        let a: i32 = 4;
        let b: i32 = 4;
        let mut c = a.mul_fallible(&b).unwrap();
        c.add_assign_fallible(&a);
        c.sub_assign_fallible(&b);
        let d = b.sub_fallible(&a);
        let e = a.add_fallible(&b);
        assert_eq!(c, 16);
        assert_eq!(d, Some(0));
        assert_eq!(e, Some(8));
    }

    #[test]
    #[cfg(feature = "shadowing")]
    fn test_fallible_mul() {
        use symbolica::atom::AtomCore;

        let a: i32 = 4;
        let b: f64 = 4.;
        let mut c: f64 = a.mul_fallible(&b).unwrap();
        c.add_assign_fallible(&a);
        let d: Option<f64> = b.mul_fallible(&a);
        let a: &i32 = &a;
        let e: Option<f64> = a.mul_fallible(&b);
        assert_eq!(c, 20.);
        assert_eq!(d, Some(16.));
        assert_eq!(e, Some(16.));

        let a = parse!("a(2)");
        let b = &parse!("b(1)");

        let mut f = a.mul_fallible(&4.).unwrap();
        f.add_assign_fallible(b);

        let i = Atom::i();

        f.add_assign_fallible(&i);

        let function_map = HashMap::new();

        let mut const_map = HashMap::new();
        const_map.insert(i.as_view(), Complex::<f64>::new(0., 1.).into());

        const_map.insert(a.as_view(), Complex::<f64>::new(3., 1.).into());

        const_map.insert(b.as_view(), Complex::<f64>::new(3., 1.).into());

        let ev: symbolica::domains::float::Complex<f64> =
            f.evaluate(|r| r.into(), &const_map, &function_map).unwrap();

        println!("{}", ev);
        // print!("{}", f.unwrap());

        let g = Complex::new(0.1, 3.);

        let mut h = a.sub_fallible(&g).unwrap();

        h.add_assign_fallible(&a);
        let _f = a.mul_fallible(&a);

        // Atom::default();

        println!("{}", h);
    }
}

// impl<T, U> SmallestUpgrade<U> for T
// where
//     U: From<T>,
// {
//     type LCM = U;
//     fn upgrade(self) -> Self::LCM {
//         U::from(self)
//     }
// }

// impl<T, U> SmallestUpgrade<Up<T>> for Up<U>
// where
//     T: From<U>,
// {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         T::from(self.up)
//     }
// } We can't do this because of possible future impls, means that any specialization is forbidden

// impl<T> SmallestUpgrade<T> for T {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// } We don't want this, so that we can specialize binary operations

// impl<T, U> SmallestUpgrade<T> for U
// where
//     T: SmallestUpgrade<U, LCM = U>,
// {
//     type LCM = U;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// } This should work but doesn't

// impl<T, U> SmallestUpgrade<Up<T>> for Up<U>
// where
//     T: From<U>,
// {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         T::from(self.up)
//     }
// } We can't do this because of possible future impls, means that any specialization is forbidden

// impl<T> SmallestUpgrade<T> for T {
//     type LCM = T;
//     fn upgrade(self) -> Self::LCM {
//         self
//     }
// } We don't want this, so that we can specialize binary operations
