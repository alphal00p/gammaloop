use std::ops::SubAssign;

use ref_ops::RefSub;

use super::Complex;

impl<T> SubAssign for Complex<T>
where
    for<'a> T: SubAssign<&'a T>,
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.sub_assign(&rhs)
    }
}

impl<T> SubAssign<T> for Complex<T>
where
    T: SubAssign<T>,
{
    #[inline]
    fn sub_assign(&mut self, rhs: T) {
        self.re -= rhs;
    }
}

impl<T> SubAssign<&Complex<T>> for Complex<T>
where
    for<'a> T: SubAssign<&'a T>,
{
    #[inline]
    fn sub_assign(&mut self, rhs: &Self) {
        self.re -= &rhs.re;
        self.im -= &rhs.im;
    }
}

impl<T> SubAssign<&T> for Complex<T>
where
    T: for<'a> RefSub<&'a T, Output = T>,
{
    #[inline]
    fn sub_assign(&mut self, rhs: &T) {
        self.re = self.re.ref_sub(rhs);
    }
}
