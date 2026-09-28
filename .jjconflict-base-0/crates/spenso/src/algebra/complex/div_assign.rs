use std::ops::{Add, Div, DivAssign, Sub};

use ref_ops::{RefDiv, RefMul};

use super::Complex;

impl<T> DivAssign for Complex<T>
where
    T: Clone
        + Sub<T, Output = T>
        + Div<T, Output = T>
        + for<'a> RefMul<&'a T, Output = T>
        + Add<T, Output = T>,
{
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        self.div_assign(&rhs)
    }
}

impl<T> DivAssign<T> for Complex<T>
where
    T: DivAssign<T> + Clone,
{
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        self.re /= rhs.clone();
        self.im /= rhs;
    }
}

impl<T> DivAssign<&Complex<T>> for Complex<T>
where
    T: Clone
        + Sub<T, Output = T>
        + Div<T, Output = T>
        + for<'a> RefMul<&'a T, Output = T>
        + Add<T, Output = T>,
{
    #[inline]
    fn div_assign(&mut self, rhs: &Self) {
        *self = &*self / rhs;
    }
}

impl<T> DivAssign<&T> for Complex<T>
where
    T: for<'a> RefDiv<&'a T, Output = T>,
{
    #[inline]
    fn div_assign(&mut self, rhs: &T) {
        self.re = self.re.ref_div(rhs);
        self.im = self.im.ref_div(rhs);
    }
}
