use std::ops::Div;

use symbolica::atom::{Atom, AtomView};

use crate::tensors::parametric::{to_param::ToAtom, ConcreteOrParam, ConcreteOrParamRef};

impl<T, U, Out> Div<ConcreteOrParamRef<'_, U>> for ConcreteOrParam<T>
where
    Atom: Div<U, Output = Atom> + for<'a> Div<AtomView<'a>, Output = Atom>,
    T: ToAtom + Div<U, Output = Out>,
{
    type Output = ConcreteOrParam<Out>;
    fn div(self, rhs: ConcreteOrParamRef<'_, U>) -> ConcreteOrParam<Out> {
        match (self, rhs) {
            (ConcreteOrParam::Concrete(c1), ConcreteOrParamRef::Concrete(c2)) => {
                ConcreteOrParam::Concrete(c1 / c2)
            }
            (ConcreteOrParam::Param(a), ConcreteOrParamRef::Concrete(c)) => {
                ConcreteOrParam::Param(a / c)
            }
            (ConcreteOrParam::Param(a1), ConcreteOrParamRef::Param(a2)) => {
                ConcreteOrParam::Param(a1 / a2)
            }
            (mut c, a) => {
                c.to_param();
                c / a
            }
        }
    }
}

// impl<T: Into<Coefficient>> Div<Complex<T>> for Atom {
//     type Output = Atom;
//     fn div(self, rhs: Complex<T>) -> Atom {
//         self / (Atom::num(rhs.re) + Atom::i() * rhs.im)
//     }
// }

// impl<T> Div<&Complex<T>> for Atom
// where
//     for<'a> &'a T: Into<Coefficient>,
// {
//     type Output = Atom;
//     fn div(self, rhs: &Complex<T>) -> Atom {
//         self / (Atom::num(&rhs.re) + Atom::i() * &rhs.im)
//     }
// }

impl<T, U, Out> Div<ConcreteOrParam<T>> for ConcreteOrParam<U>
where
    Atom: Div<T, Output = Atom> + Div<Atom, Output = Atom>,
    U: Div<T, Output = Out> + ToAtom,
{
    type Output = ConcreteOrParam<Out>;
    fn div(self, rhs: ConcreteOrParam<T>) -> Self::Output {
        match (self, rhs) {
            (ConcreteOrParam::Concrete(c1), ConcreteOrParam::Concrete(c2)) => {
                ConcreteOrParam::Concrete(c1 / c2)
            }
            (ConcreteOrParam::Param(a), ConcreteOrParam::Concrete(c)) => {
                ConcreteOrParam::Param(a / c)
            }
            (ConcreteOrParam::Param(a1), ConcreteOrParam::Param(a2)) => {
                ConcreteOrParam::Param(a1 / a2)
            }
            (mut c, a) => {
                c.to_param();
                c / a
            }
        }
    }
}
