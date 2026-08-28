use std::ops::Sub;

use ref_ops::RefSub;

use super::Complex;

impl<'a, T> Sub<&'a Complex<T>> for &Complex<T>
where
    T: Sub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: &'a Complex<T>) -> Self::Output {
        self.clone() - rhs.clone()
    }
}

impl<'a, T> Sub<&'a T> for &Complex<T>
where
    T: for<'c> RefSub<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: &'a T) -> Self::Output {
        Complex {
            re: self.re.ref_sub(rhs),
            im: self.im.clone(),
        }
    }
}

impl<'a, T> Sub<&'a Complex<T>> for Complex<T>
where
    T: Sub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: &'a Complex<T>) -> Self::Output {
        self - rhs.clone()
    }
}

impl<'a, T> Sub<&'a T> for Complex<T>
where
    T: for<'c> RefSub<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: &'a T) -> Self::Output {
        Complex {
            re: self.re.ref_sub(rhs),
            im: self.im,
        }
    }
}

impl<T> Sub<Complex<T>> for &Complex<T>
where
    T: Sub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: Complex<T>) -> Self::Output {
        self.clone() - rhs
    }
}

impl<T> Sub<T> for &Complex<T>
where
    T: RefSub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        Complex {
            re: self.re.ref_sub(rhs),
            im: self.im.clone(),
        }
    }
}

impl<T> Sub for Complex<T>
where
    T: Sub<T, Output = T>,
{
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Complex::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl<T> Sub<T> for Complex<T>
where
    T: Sub<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        Complex {
            re: self.re - rhs,
            im: self.im,
        }
    }
}
