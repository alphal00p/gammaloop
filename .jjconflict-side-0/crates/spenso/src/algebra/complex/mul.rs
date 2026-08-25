use std::ops::{Add, Mul, Sub};

use ref_ops::{RefAdd, RefMul, RefSub};

use super::Complex;

impl<'a, T> Mul<&'a Complex<T>> for &Complex<T>
where
    T: for<'c> RefMul<&'c T, Output = T>
        + for<'c> RefAdd<&'c T, Output = T>
        + for<'c> RefSub<&'c T, Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: &'a Complex<T>) -> Self::Output {
        Complex::new(
            self.re.ref_mul(&rhs.re).ref_sub(&self.im.ref_mul(&rhs.im)),
            self.re.ref_mul(&rhs.im).ref_add(&self.im.ref_mul(&rhs.re)),
        )
    }
}

impl<'a, T> Mul<&'a T> for &Complex<T>
where
    T: for<'c> RefMul<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: &'a T) -> Self::Output {
        Complex::new(self.re.ref_mul(rhs), self.im.ref_mul(rhs))
    }
}

impl<'a, T> Mul<&'a Complex<T>> for Complex<T>
where
    for<'b> T: Mul<&'b T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: &'a Complex<T>) -> Self::Output {
        Complex::new(
            self.re.clone() * &rhs.re - self.im.clone() * &rhs.im,
            self.re.clone() * &rhs.im + self.im.clone() * &rhs.re,
        )
    }
}

impl<'a, T> Mul<&'a T> for Complex<T>
where
    T: for<'c> RefMul<&'c T, Output = T>,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: &'a T) -> Self::Output {
        Complex::new(self.re.ref_mul(rhs), self.im.ref_mul(rhs))
    }
}

impl<T> Mul<Complex<T>> for &Complex<T>
where
    T: Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: Complex<T>) -> Self::Output {
        self.clone() * rhs
    }
}

impl<T> Mul<T> for &Complex<T>
where
    T: RefMul<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Complex::new(self.re.ref_mul(rhs.clone()), self.im.ref_mul(rhs))
    }
}

impl<T> Mul for Complex<T>
where
    T: Mul<T, Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Clone,
{
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Complex::new(
            self.re.clone() * rhs.re.clone() - self.im.clone() * rhs.im.clone(),
            self.re.clone() * rhs.im.clone() + self.im.clone() * rhs.re.clone(),
        )
    }
}

impl<T> Mul<T> for Complex<T>
where
    T: Mul<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Complex::new(self.re * rhs.clone(), self.im * rhs)
    }
}
