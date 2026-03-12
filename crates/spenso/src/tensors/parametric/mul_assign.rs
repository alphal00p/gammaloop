use std::ops::MulAssign;

use symbolica::atom::{Atom, AtomView};

use super::{to_param::ToAtom, ConcreteOrParam, ConcreteOrParamRef};

// Scalars

impl<T, U> MulAssign<ConcreteOrParam<U>> for ConcreteOrParam<T>
where
    Atom: MulAssign<U> + MulAssign<Atom>,
    T: ToAtom + MulAssign<U>,
{
    fn mul_assign(&mut self, rhs: ConcreteOrParam<U>) {
        match (self, rhs) {
            (ConcreteOrParam::Concrete(c1), ConcreteOrParam::Concrete(c2)) => *c1 *= c2,
            (ConcreteOrParam::Param(a), ConcreteOrParam::Concrete(c)) => *a *= c,
            (ConcreteOrParam::Param(a1), ConcreteOrParam::Param(a2)) => *a1 *= a2,
            (c, a) => {
                c.to_param();
                *c *= a
            }
        }
    }
}

impl<T, U> MulAssign<&ConcreteOrParam<U>> for ConcreteOrParam<T>
where
    Atom: for<'a> MulAssign<&'a U> + for<'a> MulAssign<&'a Atom>,
    T: ToAtom + for<'a> MulAssign<&'a U>,
{
    fn mul_assign(&mut self, rhs: &ConcreteOrParam<U>) {
        match (self, rhs) {
            (ConcreteOrParam::Concrete(c1), ConcreteOrParam::Concrete(c2)) => *c1 *= c2,
            (ConcreteOrParam::Param(a), ConcreteOrParam::Concrete(c)) => *a *= c,
            (ConcreteOrParam::Param(a1), ConcreteOrParam::Param(a2)) => *a1 *= a2,
            (c, a) => {
                c.to_param();
                *c *= a
            }
        }
    }
}

impl<T, U> MulAssign<ConcreteOrParamRef<'_, U>> for ConcreteOrParam<T>
where
    Atom: MulAssign<U> + for<'a> MulAssign<AtomView<'a>>,
    T: ToAtom + MulAssign<U>,
{
    fn mul_assign(&mut self, rhs: ConcreteOrParamRef<'_, U>) {
        match (self, rhs) {
            (ConcreteOrParam::Concrete(c1), ConcreteOrParamRef::Concrete(c2)) => *c1 *= c2,
            (ConcreteOrParam::Param(a), ConcreteOrParamRef::Concrete(c)) => *a *= c,
            (ConcreteOrParam::Param(a1), ConcreteOrParamRef::Param(a2)) => *a1 *= a2,
            (c, a) => {
                c.to_param();
                *c *= a
            }
        }
    }
}

// Atom

impl<T> MulAssign<ConcreteOrParamRef<'_, T>> for Atom
where
    Atom: MulAssign<T> + for<'a> MulAssign<AtomView<'a>>,
{
    fn mul_assign(&mut self, rhs: ConcreteOrParamRef<'_, T>) {
        match rhs {
            ConcreteOrParamRef::Concrete(c) => *self *= c,
            ConcreteOrParamRef::Param(a) => *self *= a,
        }
    }
}

// impl<T> MulAssign<RealOrComplexRef<'_, T>> for Atom
// where
//     Atom: for<'a> MulAssign<&'a T> + for<'a> MulAssign<&'a Complex<T>>,
// {
//     fn mul_assign(&mut self, rhs: RealOrComplexRef<'_, T>) {
//         match rhs {
//             RealOrComplexRef::Real(c) => *self *= c,
//             RealOrComplexRef::Complex(a) => *self *= a,
//         }
//     }
// }

// impl<T> MulAssign<&RealOrComplex<T>> for Atom
// where
//     Atom: for<'a> MulAssign<&'a T> + for<'a> MulAssign<&'a Complex<T>>,
// {
//     fn mul_assign(&mut self, rhs: &RealOrComplex<T>) {
//         match rhs {
//             RealOrComplex::Real(c) => *self *= c,
//             RealOrComplex::Complex(a) => *self *= a,
//         }
//     }
// }

// impl<T> MulAssign<RealOrComplex<T>> for Atom
// where
//     Atom: MulAssign<T> + MulAssign<Complex<T>>,
// {
//     fn mul_assign(&mut self, rhs: RealOrComplex<T>) {
//         match rhs {
//             RealOrComplex::Real(c) => *self *= c,
//             RealOrComplex::Complex(a) => *self *= a,
//         }
//     }
// }

// impl<T: Into<Coefficient>> MulAssign<Complex<T>> for Atom {
//     fn mul_assign(&mut self, rhs: Complex<T>) {
//         *self = &*self * (Atom::i() * rhs.im + rhs.re);
//     }
// }

// impl<T> MulAssign<&Complex<T>> for Atom
// where
//     for<'a> &'a T: Into<Coefficient>,
// {
//     fn mul_assign(&mut self, rhs: &Complex<T>) {
//         *self = &*self + &rhs.re + Atom::i() + &rhs.im;
//     }
// }

impl<T> MulAssign<ConcreteOrParam<T>> for Atom
where
    Atom: MulAssign<T> + MulAssign<Atom>,
{
    fn mul_assign(&mut self, rhs: ConcreteOrParam<T>) {
        match rhs {
            ConcreteOrParam::Concrete(c) => *self *= c,
            ConcreteOrParam::Param(a) => *self *= a,
        }
    }
}
