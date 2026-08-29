use std::ops::{Add, Mul, MulAssign};

use ref_ops::{RefAdd, RefMul, RefSub};

use crate::algebra::algebraic_traits::RefZero;

use super::{Complex, RealOrComplex, RealOrComplexRef};

impl<T> MulAssign for Complex<T>
where
    T: Mul<T, Output = T> + Add<T, Output = T> + std::ops::Sub<T, Output = T> + Clone,
{
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

impl<T> MulAssign<T> for Complex<T>
where
    T: MulAssign<T> + Clone,
{
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        self.re *= rhs.clone();
        self.im *= rhs;
    }
}

impl<T> MulAssign<&Complex<T>> for Complex<T>
where
    T: for<'c> RefMul<&'c T, Output = T>
        + for<'c> RefAdd<&'c T, Output = T>
        + for<'c> RefSub<&'c T, Output = T>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: &Self) {
        *self = &*self * rhs;
    }
}

impl<T> MulAssign<&T> for Complex<T>
where
    T: for<'c> RefMul<&'c T, Output = T>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: &T) {
        self.re = self.re.ref_mul(rhs);
        self.im = self.im.ref_mul(rhs);
    }
}

impl<T, U> MulAssign<&RealOrComplex<T>> for RealOrComplex<U>
where
    U: for<'a> MulAssign<&'a T> + RefZero,
    Complex<U>: for<'a> MulAssign<&'a Complex<T>> + for<'a> MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: &RealOrComplex<T>) {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplex::Complex(b)) => {
                *a *= b;
            }
            (RealOrComplex::Real(a), RealOrComplex::Real(b)) => {
                *a *= b;
            }
            (RealOrComplex::Complex(a), RealOrComplex::Real(b)) => {
                *a *= b;
            }
            (a, b) => {
                a.to_complex_mut();
                *a *= b;
            }
        }
    }
}

impl<T, U> MulAssign<RealOrComplexRef<'_, T>> for RealOrComplex<U>
where
    U: for<'a> MulAssign<&'a T> + RefZero,
    Complex<U>: for<'a> MulAssign<&'a Complex<T>> + for<'a> MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: RealOrComplexRef<'_, T>) {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplexRef::Complex(b)) => {
                *a *= b;
            }
            (RealOrComplex::Real(a), RealOrComplexRef::Real(b)) => {
                *a *= b;
            }
            (RealOrComplex::Complex(a), RealOrComplexRef::Real(b)) => {
                *a *= b;
            }
            (a, b) => {
                a.to_complex_mut();
                *a *= b;
            }
        }
    }
}

impl<T, U> MulAssign<RealOrComplexRef<'_, T>> for Complex<U>
where
    U: for<'a> MulAssign<&'a T> + RefZero,
    Complex<U>: for<'a> MulAssign<&'a Complex<T>> + for<'a> MulAssign<&'a T>,
{
    fn mul_assign(&mut self, rhs: RealOrComplexRef<'_, T>) {
        match rhs {
            RealOrComplexRef::Complex(b) => {
                *self *= b;
            }
            RealOrComplexRef::Real(b) => {
                self.re *= b;
            }
        }
    }
}
