use std::ops::AddAssign;

use ref_ops::RefAdd;

use crate::algebra::algebraic_traits::RefZero;

use super::{Complex, RealOrComplex, RealOrComplexRef};

impl<T> AddAssign for Complex<T>
where
    for<'a> T: AddAssign<&'a T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign(&rhs)
    }
}

impl<T> AddAssign<T> for Complex<T>
where
    T: AddAssign<T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        self.re += rhs;
    }
}

impl<T> AddAssign<&Complex<T>> for Complex<T>
where
    for<'a> T: AddAssign<&'a T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: &Self) {
        self.re += &rhs.re;
        self.im += &rhs.im;
    }
}

impl<T> AddAssign<&T> for Complex<T>
where
    T: for<'a> RefAdd<&'a T, Output = T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: &T) {
        self.re = self.re.ref_add(rhs);
    }
}

impl<T, U> AddAssign<&RealOrComplex<T>> for RealOrComplex<U>
where
    U: for<'a> AddAssign<&'a T> + RefZero,
    Complex<U>: for<'a> AddAssign<&'a Complex<T>> + for<'a> AddAssign<&'a T>,
{
    fn add_assign(&mut self, rhs: &RealOrComplex<T>) {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplex::Complex(b)) => {
                *a += b;
            }
            (RealOrComplex::Real(a), RealOrComplex::Real(b)) => {
                *a += b;
            }
            (RealOrComplex::Complex(a), RealOrComplex::Real(b)) => {
                *a += b;
            }
            (a, b) => {
                a.to_complex_mut();
                *a += b;
            }
        }
    }
}

impl<T, U> AddAssign<RealOrComplexRef<'_, T>> for RealOrComplex<U>
where
    U: for<'a> AddAssign<&'a T> + RefZero,
    Complex<U>: for<'a> AddAssign<&'a Complex<T>> + for<'a> AddAssign<&'a T>,
{
    fn add_assign(&mut self, rhs: RealOrComplexRef<'_, T>) {
        match (self, rhs) {
            (RealOrComplex::Complex(a), RealOrComplexRef::Complex(b)) => {
                *a += b;
            }
            (RealOrComplex::Real(a), RealOrComplexRef::Real(b)) => {
                *a += b;
            }
            (RealOrComplex::Complex(a), RealOrComplexRef::Real(b)) => {
                *a += b;
            }
            (a, b) => {
                a.to_complex_mut();
                *a += b;
            }
        }
    }
}

impl<T, U> AddAssign<RealOrComplexRef<'_, T>> for Complex<U>
where
    U: for<'a> AddAssign<&'a T> + RefZero,
    Complex<U>: for<'a> AddAssign<&'a Complex<T>> + for<'a> AddAssign<&'a T>,
{
    fn add_assign(&mut self, rhs: RealOrComplexRef<'_, T>) {
        match rhs {
            RealOrComplexRef::Complex(b) => {
                *self += b;
            }
            RealOrComplexRef::Real(b) => {
                self.re += b;
            }
        }
    }
}
