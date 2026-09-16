use std::ops::Neg;

use crate::{
    algebra::algebraic_traits::IsZero,
    algebra::complex::Complex,
    algebra::upgrading_arithmetic::{FallibleSub, TrySmallestUpgrade},
    iterators::IteratableTensor,
    structure::{HasStructure, TensorStructure},
    tensors::complex::RealOrComplexTensor,
    tensors::data::{DataTensor, DenseTensor, GetTensorData, SparseTensor},
};

impl<T, U, I, Out> FallibleSub<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: FallibleSub<T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(&self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(
            self.structure().same_external(rhs.structure()),
            "{} != {}",
            self.structure().string_rep(),
            rhs.structure().string_rep()
        );
        let self_to_rhs = rhs.structure().find_permutation(self.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter_expanded()
            .map(|(indices, u)| {
                let permuted_indices = indices.apply_permutation(&self_to_rhs);
                let t = rhs.get_ref(&permuted_indices).unwrap();
                u.sub_fallible(t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<T, U, I, Out> FallibleSub<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: FallibleSub<T, Output = Out>,
    I: TensorStructure + Clone,
    T: TrySmallestUpgrade<U, LCM = Out>,
    Out: Neg<Output = Out> + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(&self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(
            self.structure().same_external(rhs.structure()),
            "dense sparse: {} != {}",
            self.structure().string_rep(),
            rhs.structure().string_rep()
        );
        let rhs_to_self = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .iter_expanded()
            .map(|(indices, t)| {
                let permuted_indices = indices.apply_permutation(&rhs_to_self);
                let u = self.get_ref(&permuted_indices);
                if let Ok(u) = u {
                    u.sub_fallible(t)
                } else {
                    Some(t.try_upgrade().unwrap().into_owned().neg())
                }
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<T, U, I, Out> FallibleSub<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: FallibleSub<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    I: TensorStructure + Clone,
    T: Clone,
    Out: Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(&self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(
            self.structure().same_external(rhs.structure()),
            "sparse dense:{} != {}",
            self.structure().string_rep(),
            rhs.structure().string_rep()
        );
        let self_to_rhs = rhs.structure().find_permutation(self.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter_expanded()
            .map(|(indices, u)| {
                let permuted_indices = indices.apply_permutation(&self_to_rhs);
                let t = rhs.get_ref(&permuted_indices);
                if let Ok(t) = t {
                    u.sub_fallible(t)
                } else {
                    Some(u.try_upgrade().unwrap().into_owned())
                }
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<T, U, I, Out> FallibleSub<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: FallibleSub<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    I: TensorStructure + Clone,
    T: Clone + TrySmallestUpgrade<U, LCM = Out>,
    Out: IsZero + Clone + Neg<Output = Out>,
{
    type Output = SparseTensor<Out, I>;
    fn sub_fallible(&self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(
            self.structure().same_external(rhs.structure()),
            "sparse sparse:{} != {}",
            self.structure().string_rep(),
            rhs.structure().string_rep()
        );
        let structure = self.structure().clone();
        let mut data = SparseTensor::empty(structure, self.zero.sub_fallible(&rhs.zero)?);
        let self_to_rhs = rhs.structure().find_permutation(self.structure()).unwrap();

        for (indices, u) in self.iter_expanded() {
            let permuted_indices = indices.apply_permutation(&self_to_rhs);
            let t = rhs.get_ref(&permuted_indices);
            if let Ok(t) = t {
                data.smart_set(&indices, u.sub_fallible(t)?).unwrap();
            } else {
                data.smart_set(&indices, u.try_upgrade().unwrap().into_owned())
                    .unwrap();
            }
        }

        for (i, t) in rhs.iter_expanded() {
            let permuted_indices = i.apply_inverse_permutation(&self_to_rhs);

            if self.get_ref(&permuted_indices).is_err() {
                data.smart_set(&i, t.try_upgrade().unwrap().into_owned().neg())
                    .unwrap();
            }
        }

        Some(data)
    }
}

impl<T, U, Out, I> FallibleSub<DataTensor<T, I>> for DataTensor<U, I>
where
    U: FallibleSub<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    T: Clone + TrySmallestUpgrade<U, LCM = Out>,
    Out: IsZero + Clone + Neg<Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DataTensor<Out, I>;
    fn sub_fallible(&self, rhs: &DataTensor<T, I>) -> Option<Self::Output> {
        // println!("sub_fallible");
        match (self, rhs) {
            (DataTensor::Dense(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.sub_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Sparse(a.sub_fallible(b)?))
            }
            (DataTensor::Dense(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Dense(a.sub_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.sub_fallible(b)?))
            }
        }
    }
}

impl<T: Clone, S: TensorStructure + Clone> FallibleSub<RealOrComplexTensor<T, S>>
    for RealOrComplexTensor<T, S>
where
    T: FallibleSub<T, Output = T>
        + TrySmallestUpgrade<T, LCM = T>
        + IsZero
        + Neg<Output = T>
        + Clone
        + FallibleSub<Complex<T>, Output = Complex<T>>
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>,
    Complex<T>: FallibleSub<T, Output = Complex<T>>
        + FallibleSub<Complex<T>, Output = Complex<T>>
        + IsZero
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>
        + TrySmallestUpgrade<T, LCM = Complex<T>>
        + Neg<Output = Complex<T>>
        + Clone,
{
    type Output = RealOrComplexTensor<T, S>;

    fn sub_fallible(&self, rhs: &RealOrComplexTensor<T, S>) -> Option<Self::Output> {
        match (self, rhs) {
            (RealOrComplexTensor::Real(s), RealOrComplexTensor::Real(r)) => {
                Some(RealOrComplexTensor::Real(s.sub_fallible(r)?))
            }
            (RealOrComplexTensor::Real(s), RealOrComplexTensor::Complex(r)) => {
                Some(RealOrComplexTensor::Complex(s.sub_fallible(r)?))
            }
            (RealOrComplexTensor::Complex(s), RealOrComplexTensor::Real(r)) => {
                Some(RealOrComplexTensor::Complex(s.sub_fallible(r)?))
            }
            (RealOrComplexTensor::Complex(s), RealOrComplexTensor::Complex(r)) => {
                Some(RealOrComplexTensor::Complex(s.sub_fallible(r)?))
            }
        }
    }
}
#[cfg(feature = "shadowing")]
use crate::tensors::parametric::{MixedTensor, ParamTensor};
#[cfg(feature = "shadowing")]
use symbolica::atom::Atom;
#[cfg(feature = "shadowing")]
impl<S: TensorStructure + Clone> FallibleSub<ParamTensor<S>> for ParamTensor<S> {
    type Output = ParamTensor<S>;
    fn sub_fallible(&self, rhs: &ParamTensor<S>) -> Option<Self::Output> {
        let t = self.tensor.sub_fallible(&rhs.tensor)?;
        Some(ParamTensor::composite(t))
    }
}

#[cfg(feature = "shadowing")]
impl<I, T> FallibleSub<MixedTensor<T, I>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
    T: FallibleSub<T, Output = T>
        + FallibleSub<Atom, Output = Atom>
        + TrySmallestUpgrade<T, LCM = T>
        + IsZero
        + Clone
        + Neg<Output = T>
        + FallibleSub<Complex<T>, Output = Complex<T>>
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>
        + TrySmallestUpgrade<Atom, LCM = Atom>,
    Complex<T>: FallibleSub<T, Output = Complex<T>>
        + FallibleSub<Complex<T>, Output = Complex<T>>
        + FallibleSub<Atom, Output = Atom>
        + IsZero
        + Neg<Output = Complex<T>>
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>
        + TrySmallestUpgrade<T, LCM = Complex<T>>
        + TrySmallestUpgrade<Atom, LCM = Atom>
        + Clone,
    Atom: TrySmallestUpgrade<T, LCM = Atom>
        + TrySmallestUpgrade<Complex<T>, LCM = Atom>
        + FallibleSub<T, Output = Atom>
        + FallibleSub<Complex<T>, Output = Atom>,
{
    type Output = MixedTensor<T, I>;
    fn sub_fallible(&self, rhs: &MixedTensor<T, I>) -> Option<Self::Output> {
        match (self, rhs) {
            (MixedTensor::Param(a), MixedTensor::Param(b)) => {
                Some(MixedTensor::Param(a.sub_fallible(b)?))
            }
            (MixedTensor::Param(s), MixedTensor::Concrete(o)) => match o {
                RealOrComplexTensor::Real(o) => Some(MixedTensor::Param(ParamTensor::composite(
                    s.tensor.sub_fallible(o)?,
                ))),
                RealOrComplexTensor::Complex(o) => Some(MixedTensor::Param(
                    ParamTensor::composite(s.tensor.sub_fallible(o)?),
                )),
            },
            (MixedTensor::Concrete(s), MixedTensor::Param(o)) => match s {
                RealOrComplexTensor::Real(s) => Some(MixedTensor::Param(ParamTensor::composite(
                    s.sub_fallible(&o.tensor)?,
                ))),
                RealOrComplexTensor::Complex(s) => Some(MixedTensor::Param(
                    ParamTensor::composite(s.sub_fallible(&o.tensor)?),
                )),
            },
            (MixedTensor::Concrete(s), MixedTensor::Concrete(o)) => {
                Some(MixedTensor::Concrete(s.sub_fallible(o)?))
            }
        }
    }
}
