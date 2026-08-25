use std::ops::{Add, Div, Sub};

use ref_ops::{RefDiv, RefMul};

use crate::algebra::{
    algebraic_traits::RefZero,
    complex::{RealOrComplex, RealOrComplexRef},
};

use super::Complex;

impl<'a, T> Div<&'a Complex<T>> for &Complex<T>
where
    T: Clone
        + Sub<T, Output = T>
        + Div<T, Output = T>
        + for<'r> RefMul<&'r T, Output = T>
        + Add<T, Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: &'a Complex<T>) -> Self::Output {
        let n = rhs.norm_squared();
        let re = self.re.ref_mul(&rhs.re) + self.im.ref_mul(&rhs.im);
        let im = self.im.ref_mul(&rhs.re) - self.re.ref_mul(&rhs.im);
        Complex::new(re / n.clone(), im / n)
    }
}

impl<'a, T> Div<&'a T> for &Complex<T>
where
    T: for<'c> RefDiv<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: &'a T) -> Self::Output {
        Complex::new(self.re.ref_div(rhs), self.im.ref_div(rhs))
    }
}

impl<'a, T> Div<&'a Complex<T>> for Complex<T>
where
    T: Clone
        + Sub<T, Output = T>
        + Div<T, Output = T>
        + for<'b> RefMul<&'b T, Output = T>
        + Add<T, Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: &'a Complex<T>) -> Self::Output {
        &self / rhs
    }
}

impl<'a, T> Div<&'a T> for Complex<T>
where
    T: for<'c> Div<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: &'a T) -> Self::Output {
        Complex::new(self.re.div(rhs), self.im.div(rhs))
    }
}

impl<T> Div<Complex<T>> for &Complex<T>
where
    T: Clone
        + Sub<T, Output = T>
        + Div<T, Output = T>
        + for<'b> RefMul<&'b T, Output = T>
        + Add<T, Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: Complex<T>) -> Self::Output {
        self / &rhs
    }
}

impl<T> Div<T> for &Complex<T>
where
    T: RefDiv<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Complex::new(self.re.ref_div(rhs.clone()), self.im.ref_div(rhs))
    }
}

impl<T> Div for Complex<T>
where
    T: Clone
        + Sub<T, Output = T>
        + Div<T, Output = T>
        + for<'a> RefMul<&'a T, Output = T>
        + Add<T, Output = T>,
{
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        &self / &rhs
    }
}

impl<T> Div<T> for Complex<T>
where
    T: Div<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        Complex::new(self.re.div(rhs.clone()), self.im.div(rhs))
    }
}

impl<T, U, Out> Div<RealOrComplex<T>> for RealOrComplex<U>
where
    U: Div<T, Output = Out> + RefZero,
    Complex<U>: Div<Complex<T>, Output = Complex<Out>> + Div<T, Output = Complex<Out>>,
{
    type Output = RealOrComplex<Out>;
    fn div(self, rhs: RealOrComplex<T>) -> RealOrComplex<Out> {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplex::Complex(b)) => RealOrComplex::Complex(a / b),
            (RealOrComplex::Real(a), RealOrComplex::Real(b)) => RealOrComplex::Real(a / b),
            (RealOrComplex::Complex(a), RealOrComplex::Real(b)) => RealOrComplex::Complex(a / b),
            (mut a, b) => {
                a.to_complex_mut();
                a / b
            }
        }
    }
}

impl<T, U, Out> Div<&RealOrComplex<T>> for RealOrComplex<U>
where
    U: for<'a> Div<&'a T, Output = Out> + RefZero,
    Complex<U>: for<'a> Div<&'a Complex<T>, Output = Complex<Out>>
        + for<'a> Div<&'a T, Output = Complex<Out>>,
{
    type Output = RealOrComplex<Out>;
    fn div(self, rhs: &RealOrComplex<T>) -> RealOrComplex<Out> {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplex::Complex(b)) => RealOrComplex::Complex(a / b),
            (RealOrComplex::Real(a), RealOrComplex::Real(b)) => RealOrComplex::Real(a / b),
            (RealOrComplex::Complex(a), RealOrComplex::Real(b)) => RealOrComplex::Complex(a / b),
            (mut a, b) => {
                a.to_complex_mut();
                a / b
            }
        }
    }
}

impl<T, U, Out> Div<RealOrComplexRef<'_, T>> for RealOrComplex<U>
where
    U: for<'a> Div<&'a T, Output = Out> + RefZero,
    Complex<U>: for<'a> Div<&'a Complex<T>, Output = Complex<Out>>
        + for<'a> Div<&'a T, Output = Complex<Out>>,
{
    type Output = RealOrComplex<Out>;
    fn div(self, rhs: RealOrComplexRef<'_, T>) -> RealOrComplex<Out> {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplexRef::Complex(b)) => {
                RealOrComplex::Complex(a / b)
            }
            (RealOrComplex::Real(a), RealOrComplexRef::Real(b)) => RealOrComplex::Real(a / b),
            (RealOrComplex::Complex(a), RealOrComplexRef::Real(b)) => RealOrComplex::Complex(a / b),
            (mut a, b) => {
                a.to_complex_mut();
                a / b
            }
        }
    }
}

impl<T, U, Out> Div<&RealOrComplex<T>> for Complex<U>
where
    Complex<U>: for<'a> Div<&'a Complex<T>, Output = Complex<Out>>
        + for<'a> Div<&'a T, Output = Complex<Out>>,
{
    type Output = Complex<Out>;
    fn div(self, rhs: &RealOrComplex<T>) -> Complex<Out> {
        match rhs {
            RealOrComplex::Complex(b) => self / b,
            RealOrComplex::Real(a) => self / a,
        }
    }
}
impl<T, U, Out> Div<RealOrComplexRef<'_, T>> for Complex<U>
where
    Complex<U>: for<'a> Div<&'a Complex<T>, Output = Complex<Out>>
        + for<'a> Div<&'a T, Output = Complex<Out>>,
{
    type Output = Complex<Out>;
    fn div(self, rhs: RealOrComplexRef<'_, T>) -> Complex<Out> {
        match rhs {
            RealOrComplexRef::Complex(b) => self / b,
            RealOrComplexRef::Real(a) => self / a,
        }
    }
}

#[cfg(test)]
mod test {
    // use crate::algebra::complex::Complex;

    // fn div() {
    //     let a = Complex::new(1.0, 2.0);
    //     a / 2.0;
    // }
}
