use std::ops::Neg;

use super::{Complex, RealOrComplex, RealOrComplexRef};

impl<T> Neg for Complex<T>
where
    T: Neg<Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn neg(self) -> Complex<T> {
        Complex::new(-self.re, -self.im)
    }
}

impl<'a, T> Neg for &'a Complex<T>
where
    &'a T: Neg<Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn neg(self) -> Complex<T> {
        Complex::new(-&self.re, -&self.im)
    }
}

impl<T: Neg<Output = T>> Neg for RealOrComplex<T> {
    type Output = RealOrComplex<T>;

    fn neg(self) -> Self::Output {
        match self {
            RealOrComplex::Real(r) => RealOrComplex::Real(-r),
            RealOrComplex::Complex(c) => RealOrComplex::Complex(-c),
        }
    }
}

impl<'a, T> Neg for RealOrComplexRef<'a, T>
where
    &'a T: Neg<Output = T>,
{
    type Output = RealOrComplex<T>;

    fn neg(self) -> Self::Output {
        match self {
            RealOrComplexRef::Real(r) => RealOrComplex::Real(-r),
            RealOrComplexRef::Complex(c) => RealOrComplex::Complex(-c),
        }
    }
}
