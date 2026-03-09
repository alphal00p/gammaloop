use crate::{
    algebra::algebraic_traits::RefZero,
    algebra::complex::{Complex, RealOrComplex},
    algebra::upgrading_arithmetic::{FallibleMul, TrySmallestUpgrade},
    algebra::ScalarMul,
    shadowing::symbolica_utils::SerializableAtom,
    structure::TensorStructure,
    tensors::complex::RealOrComplexTensor,
    tensors::data::DataTensor,
};
use duplicate::duplicate;
use symbolica::{atom::Atom, domains::rational::Rational};

use super::{ConcreteOrParam, MixedTensor, ParamOrComposite, ParamOrConcrete, ParamTensor};

impl<S: TensorStructure + Clone> ScalarMul<SerializableAtom> for ParamTensor<S> {
    type Output = ParamTensor<S>;
    fn scalar_mul(&self, rhs: &SerializableAtom) -> Option<Self::Output> {
        Some(ParamTensor {
            tensor: self.tensor.scalar_mul(&rhs.0)?,
            param_type: ParamOrComposite::Composite,
        })
    }
}

impl<S: TensorStructure + Clone, T> ScalarMul<T> for ParamTensor<S>
where
    DataTensor<Atom, S>: ScalarMul<T, Output = DataTensor<Atom, S>>,
{
    type Output = ParamTensor<S>;
    fn scalar_mul(&self, rhs: &T) -> Option<Self::Output> {
        Some(ParamTensor::composite(self.tensor.scalar_mul(rhs)?))
    }
}

impl<T> TrySmallestUpgrade<Complex<T>> for Atom
where
    Atom: TrySmallestUpgrade<T, LCM = Atom>,
{
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_, Self::LCM>> {
        Some(std::borrow::Cow::Borrowed(self))
    }
}

impl<T> TrySmallestUpgrade<Atom> for Complex<T>
where
    T: TrySmallestUpgrade<Atom, LCM = Atom>,
{
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_, Self::LCM>> {
        Some(std::borrow::Cow::Owned(
            self.re.try_upgrade()?.as_ref() + self.im.try_upgrade()?.as_ref() * Atom::i(),
        ))
    }
}

impl<T> TrySmallestUpgrade<RealOrComplex<T>> for Atom
where
    Atom: TrySmallestUpgrade<T, LCM = Atom>,
{
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_, Self::LCM>> {
        Some(std::borrow::Cow::Borrowed(self))
    }
}

impl<T> TrySmallestUpgrade<Atom> for RealOrComplex<T>
where
    T: TrySmallestUpgrade<Atom, LCM = Atom>,
{
    type LCM = Atom;

    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_, Self::LCM>> {
        match self {
            RealOrComplex::Real(r) => r.try_upgrade(),
            RealOrComplex::Complex(c) => Some(std::borrow::Cow::Owned(
                c.re.try_upgrade()?.as_ref() + c.im.try_upgrade()?.as_ref() * Atom::i(),
            )),
        }
    }
}

impl<T> TrySmallestUpgrade<ConcreteOrParam<T>> for Atom
where
    Atom: TrySmallestUpgrade<T>,
{
    type LCM = <Atom as TrySmallestUpgrade<T>>::LCM;

    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_, Self::LCM>> {
        <Atom as TrySmallestUpgrade<T>>::try_upgrade(self)
    }
}

impl<T> TrySmallestUpgrade<Atom> for ConcreteOrParam<T>
where
    T: TrySmallestUpgrade<Atom>,
    Atom: TrySmallestUpgrade<T, LCM = T::LCM>,
{
    type LCM = <T as TrySmallestUpgrade<Atom>>::LCM;

    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_, Self::LCM>> {
        match self {
            ConcreteOrParam::Param(a) => a.try_upgrade(),
            ConcreteOrParam::Concrete(c) => c.try_upgrade(),
        }
    }
}

duplicate! {
    [smaller larger;

[Rational][Complex<Rational>]

]

impl TrySmallestUpgrade<smaller> for larger {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
        Some(std::borrow::Cow::Borrowed(self))
    }
}

impl TrySmallestUpgrade<larger> for smaller {
    type LCM = larger;



    fn try_upgrade(&self) -> Option<std::borrow::Cow<'_,Self::LCM>>
        where
            Self::LCM: Clone {
               let z = self.ref_zero();
                Some(std::borrow::Cow::Owned(Complex::new(self.clone(), z)))
    }
}


}

// impl<S: TensorStructure + Clone, T> ScalarMul<RealOrComplex<T>> for ParamTensor<S>
// where
//     ParamTensor<S>:
//         ScalarMul<T, Output = ParamTensor<S>> + ScalarMul<Complex<T>, Output = ParamTensor<S>>,
// {
//     type Output = ParamTensor<S>;
//     fn scalar_mul(&self, rhs: &RealOrComplex<T>) -> Option<Self::Output> {
//         Some(match rhs {
//             RealOrComplex::Complex(c) => self.scalar_mul(c)?,
//             RealOrComplex::Real(c) => self.scalar_mul(c)?,
//         })
//     }
// }

// impl<S: TensorStructure + Clone, T> ScalarMul<ConcreteOrParam<T>> for ParamTensor<S>
// where
//     ParamTensor<S>: ScalarMul<T, Output = ParamTensor<S>>,
// {
//     type Output = ParamTensor<S>;
//     fn scalar_mul(&self, rhs: &ConcreteOrParam<T>) -> Option<Self::Output> {
//         Some(match rhs {
//             ConcreteOrParam::Param(a) => ParamTensor::composite(self.tensor.scalar_mul(a)?),
//             ConcreteOrParam::Concrete(c) => self.scalar_mul(c)?,
//         })
//     }
// }

impl<T, I> ScalarMul<Atom> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
    T: Clone + FallibleMul<Atom, Output = Atom>,
    Complex<T>: FallibleMul<Atom, Output = Atom>,
{
    type Output = MixedTensor<T, I>;
    fn scalar_mul(&self, rhs: &Atom) -> Option<Self::Output> {
        match self {
            MixedTensor::Param(a) => Some(MixedTensor::Param(a.scalar_mul(rhs)?)),

            MixedTensor::Concrete(RealOrComplexTensor::Real(a)) => Some(MixedTensor::Param(
                ParamTensor::composite(a.scalar_mul(rhs)?),
            )),
            MixedTensor::Concrete(RealOrComplexTensor::Complex(a)) => Some(MixedTensor::Param(
                ParamTensor::composite(a.scalar_mul(rhs)?),
            )),
        }
    }
}

impl<U, S> ScalarMul<Atom> for RealOrComplexTensor<U, S>
where
    DataTensor<U, S>: ScalarMul<Atom, Output = DataTensor<Atom, S>>,

    DataTensor<Complex<U>, S>: ScalarMul<Atom, Output = DataTensor<Atom, S>>,
    S: TensorStructure + Clone,
{
    type Output = ParamTensor<S>;

    fn scalar_mul(&self, rhs: &Atom) -> Option<Self::Output> {
        Some(match self {
            RealOrComplexTensor::Real(rt) => ParamTensor::composite(rt.scalar_mul(rhs)?),
            RealOrComplexTensor::Complex(rt) => ParamTensor::composite(rt.scalar_mul(rhs)?),
        })
    }
}

impl<T, U, S, Out> ScalarMul<ConcreteOrParam<T>> for ParamOrConcrete<U, S>
where
    ParamTensor<S>:
        ScalarMul<T, Output = ParamTensor<S>> + ScalarMul<Atom, Output = ParamTensor<S>>,
    U: ScalarMul<T, Output = Out> + ScalarMul<Atom, Output = ParamTensor<S>>,
    S: TensorStructure + Clone,
{
    type Output = ParamOrConcrete<Out, S>;

    fn scalar_mul(&self, rhs: &ConcreteOrParam<T>) -> Option<Self::Output> {
        Some(match (self, rhs) {
            (ParamOrConcrete::Param(rt), ConcreteOrParam::Param(rs)) => {
                ParamOrConcrete::Param(rt.scalar_mul(rs)?)
            }
            (ParamOrConcrete::Concrete(rt), ConcreteOrParam::Param(rs)) => {
                ParamOrConcrete::Param(rt.scalar_mul(rs)?)
            }
            (ParamOrConcrete::Param(rt), ConcreteOrParam::Concrete(rs)) => {
                ParamOrConcrete::Param(rt.scalar_mul(rs)?)
            }
            (ParamOrConcrete::Concrete(rt), ConcreteOrParam::Concrete(rs)) => {
                ParamOrConcrete::Concrete(rt.scalar_mul(rs)?)
            }
        })
    }
}

impl<T, I> ScalarMul<SerializableAtom> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
    T: Clone + FallibleMul<Atom, Output = Atom>,
    Complex<T>: FallibleMul<Atom, Output = Atom>,
{
    type Output = MixedTensor<T, I>;
    fn scalar_mul(&self, rhs: &SerializableAtom) -> Option<Self::Output> {
        self.scalar_mul(&rhs.0)
    }
}

impl<I, T> ScalarMul<T> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
    T: FallibleMul<T, Output = T> + Clone,
    Atom: FallibleMul<T, Output = Atom>,
    Complex<T>: FallibleMul<T, Output = Complex<T>>,
{
    type Output = MixedTensor<T, I>;
    fn scalar_mul(&self, rhs: &T) -> Option<Self::Output> {
        match self {
            MixedTensor::Param(a) => Some(MixedTensor::Param(ParamTensor::composite(
                a.tensor.scalar_mul(rhs)?,
            ))),
            MixedTensor::Concrete(RealOrComplexTensor::Real(a)) => Some(MixedTensor::Concrete(
                RealOrComplexTensor::Real(a.scalar_mul(rhs)?),
            )),
            MixedTensor::Concrete(RealOrComplexTensor::Complex(a)) => Some(MixedTensor::Concrete(
                RealOrComplexTensor::Complex(a.scalar_mul(rhs)?),
            )),
        }
    }
}

impl<I, T> ScalarMul<Complex<T>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
    T: FallibleMul<Complex<T>, Output = Complex<T>> + Clone,
    Atom: FallibleMul<Complex<T>, Output = Atom>,
    Complex<T>: FallibleMul<Complex<T>, Output = Complex<T>>,
{
    type Output = MixedTensor<T, I>;
    fn scalar_mul(&self, rhs: &Complex<T>) -> Option<Self::Output> {
        match self {
            MixedTensor::Param(a) => Some(MixedTensor::Param(ParamTensor::composite(
                a.tensor.scalar_mul(rhs)?,
            ))),
            MixedTensor::Concrete(RealOrComplexTensor::Real(a)) => Some(MixedTensor::Concrete(
                RealOrComplexTensor::Complex(a.scalar_mul(rhs)?),
            )),
            MixedTensor::Concrete(RealOrComplexTensor::Complex(a)) => Some(MixedTensor::Concrete(
                RealOrComplexTensor::Complex(a.scalar_mul(rhs)?),
            )),
        }
    }
}
