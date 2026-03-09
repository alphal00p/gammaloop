use std::ops::Add;

use ref_ops::RefAdd;

use super::Complex;

impl<'a, T> Add<&'a Complex<T>> for &Complex<T>
where
    T: Add<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: &'a Complex<T>) -> Self::Output {
        self.clone() + rhs.clone()
    }
}

impl<'a, T> Add<&'a T> for &Complex<T>
where
    T: for<'c> RefAdd<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: &'a T) -> Self::Output {
        Complex {
            re: self.re.ref_add(rhs),
            im: self.im.clone(),
        }
    }
}

impl<'a, T> Add<&'a Complex<T>> for Complex<T>
where
    T: Add<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: &'a Complex<T>) -> Self::Output {
        self + rhs.clone()
    }
}

impl<'a, T> Add<&'a T> for Complex<T>
where
    T: for<'c> RefAdd<&'c T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: &'a T) -> Self::Output {
        Complex {
            re: self.re.ref_add(rhs),
            im: self.im.clone(),
        }
    }
}

impl<T> Add<Complex<T>> for &Complex<T>
where
    T: Add<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: Complex<T>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<T> Add<T> for &Complex<T>
where
    T: RefAdd<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        Complex {
            re: self.re.ref_add(rhs),
            im: self.im.clone(),
        }
    }
}

impl<T> Add for Complex<T>
where
    T: Add<T, Output = T>,
{
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Complex::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl<T> Add<T> for Complex<T>
where
    T: Add<T, Output = T> + Clone,
{
    type Output = Complex<T>;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        Complex {
            re: self.re.add(rhs),
            im: self.im,
        }
    }
}
