use crate::{
    algebra::algebraic_traits::IsZero,
    algebra::complex::Complex,
    algebra::upgrading_arithmetic::{FallibleAdd, TrySmallestUpgrade},
    iterators::IteratableTensor,
    structure::{HasStructure, TensorStructure},
    tensors::complex::RealOrComplexTensor,
    tensors::data::{DataTensor, DenseTensor, GetTensorData, SparseTensor},
};

impl<T, U, I, Out> FallibleAdd<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: FallibleAdd<T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(&self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(
            self.structure().same_external(rhs.structure()),
            "dense dense: {}!={}",
            self.structure().string_rep(),
            rhs.structure().string_rep()
        );
        // Makes rhs into self ,when applied.
        let self_to_rhs = rhs.structure().find_permutation(self.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter_expanded()
            .map(|(indices, u)| {
                let permuted_indices = indices.apply_permutation(&self_to_rhs);
                let t = rhs.get_ref(&permuted_indices).unwrap();
                u.add_fallible(t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<T, U, I, Out> FallibleAdd<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: FallibleAdd<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    I: TensorStructure + Clone,
    T: Clone,
    Out: Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(&self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let self_to_rhs = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter_expanded()
            .map(|(indices, u)| {
                let permuted_indices = indices.apply_permutation(&self_to_rhs);
                let t = rhs.get_ref(&permuted_indices);
                if let Ok(t) = t {
                    u.add_fallible(t)
                } else {
                    Some(u.try_upgrade().unwrap().into_owned())
                }
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<T, U, I, Out> FallibleAdd<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: FallibleAdd<T, Output = Out>,
    I: TensorStructure + Clone,
    T: TrySmallestUpgrade<U, LCM = Out>,

    Out: Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(&self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let rhs_to_self = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .iter_expanded()
            .map(|(indices, t)| {
                let permuted_indices = indices.apply_permutation(&rhs_to_self);
                let u = self.get_ref(&permuted_indices);
                if let Ok(u) = u {
                    u.add_fallible(t)
                } else {
                    Some(t.try_upgrade().unwrap().into_owned())
                }
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<T, U, I, Out> FallibleAdd<SparseTensor<T, I>> for SparseTensor<U, I>
where
    I: TensorStructure + Clone,
    U: FallibleAdd<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    T: Clone + TrySmallestUpgrade<U, LCM = Out>,
    Out: IsZero + Clone,
{
    type Output = SparseTensor<Out, I>;
    fn add_fallible(&self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));

        let rhs_to_self = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();
        let mut data = SparseTensor::empty(structure, self.zero.add_fallible(&rhs.zero)?);
        for (indices, u) in self.iter_expanded() {
            let permuted_indices = indices.apply_inverse_permutation(&rhs_to_self);
            let t = rhs.get_ref(&permuted_indices);
            if let Ok(t) = t {
                data.smart_set(&indices, u.add_fallible(t)?).unwrap();
            } else {
                data.smart_set(&indices, u.try_upgrade().unwrap().into_owned())
                    .unwrap();
            }
            // println!("{:?}", t);
            // data.smart_set(&indices, u.add_fallible(&t)?).unwrap();
        }

        for (i, t) in rhs.iter_expanded() {
            let permuted_indices = i.apply_permutation(&rhs_to_self);
            if self.get_ref(&permuted_indices).is_err() {
                data.smart_set(&permuted_indices, t.try_upgrade().unwrap().into_owned())
                    .unwrap();
            }
        }

        Some(data)
    }
}

impl<T, U, Out, I> FallibleAdd<DataTensor<T, I>> for DataTensor<U, I>
where
    U: FallibleAdd<T, Output = Out> + TrySmallestUpgrade<T, LCM = Out>,
    T: Clone + TrySmallestUpgrade<U, LCM = Out>,
    Out: IsZero + Clone,
    I: TensorStructure + Clone,
{
    type Output = DataTensor<Out, I>;
    fn add_fallible(&self, rhs: &DataTensor<T, I>) -> Option<Self::Output> {
        match (self, rhs) {
            (DataTensor::Dense(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.add_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Sparse(a.add_fallible(b)?))
            }
            (DataTensor::Dense(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Dense(a.add_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.add_fallible(b)?))
            }
        }
    }
}

impl<T: Clone, S: TensorStructure + Clone> FallibleAdd<RealOrComplexTensor<T, S>>
    for RealOrComplexTensor<T, S>
where
    T: FallibleAdd<T, Output = T>
        + TrySmallestUpgrade<T, LCM = T>
        + IsZero
        + Clone
        + FallibleAdd<Complex<T>, Output = Complex<T>>
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>,
    Complex<T>: FallibleAdd<T, Output = Complex<T>>
        + FallibleAdd<Complex<T>, Output = Complex<T>>
        + IsZero
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>
        + TrySmallestUpgrade<T, LCM = Complex<T>>
        + Clone,
{
    type Output = RealOrComplexTensor<T, S>;

    fn add_fallible(&self, rhs: &RealOrComplexTensor<T, S>) -> Option<Self::Output> {
        match (self, rhs) {
            (RealOrComplexTensor::Real(s), RealOrComplexTensor::Real(r)) => {
                Some(RealOrComplexTensor::Real(s.add_fallible(r)?))
            }
            (RealOrComplexTensor::Real(s), RealOrComplexTensor::Complex(r)) => {
                Some(RealOrComplexTensor::Complex(s.add_fallible(r)?))
            }
            (RealOrComplexTensor::Complex(s), RealOrComplexTensor::Real(r)) => {
                Some(RealOrComplexTensor::Complex(s.add_fallible(r)?))
            }
            (RealOrComplexTensor::Complex(s), RealOrComplexTensor::Complex(r)) => {
                Some(RealOrComplexTensor::Complex(s.add_fallible(r)?))
            }
        }
    }
}

#[cfg(feature = "shadowing")]
use crate::tensors::parametric::{MixedTensor, ParamTensor};
#[cfg(feature = "shadowing")]
use symbolica::atom::Atom;

#[cfg(feature = "shadowing")]
impl<S: TensorStructure + Clone> FallibleAdd<ParamTensor<S>> for ParamTensor<S> {
    type Output = ParamTensor<S>;
    fn add_fallible(&self, rhs: &ParamTensor<S>) -> Option<Self::Output> {
        let t = self.tensor.add_fallible(&rhs.tensor)?;
        Some(ParamTensor::composite(t))
    }
}

#[cfg(feature = "shadowing")]
impl<I, T> FallibleAdd<MixedTensor<T, I>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
    T: FallibleAdd<T, Output = T>
        + FallibleAdd<Atom, Output = Atom>
        + TrySmallestUpgrade<T, LCM = T>
        + IsZero
        + Clone
        + FallibleAdd<Complex<T>, Output = Complex<T>>
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>
        + TrySmallestUpgrade<Atom, LCM = Atom>,
    Complex<T>: FallibleAdd<T, Output = Complex<T>>
        + FallibleAdd<Complex<T>, Output = Complex<T>>
        + FallibleAdd<Atom, Output = Atom>
        + IsZero
        + TrySmallestUpgrade<Complex<T>, LCM = Complex<T>>
        + TrySmallestUpgrade<T, LCM = Complex<T>>
        + TrySmallestUpgrade<Atom, LCM = Atom>
        + Clone,
    Atom: TrySmallestUpgrade<T, LCM = Atom>
        + TrySmallestUpgrade<Complex<T>, LCM = Atom>
        + FallibleAdd<T, Output = Atom>
        + FallibleAdd<Complex<T>, Output = Atom>,
{
    type Output = MixedTensor<T, I>;
    fn add_fallible(&self, rhs: &MixedTensor<T, I>) -> Option<Self::Output> {
        match (self, rhs) {
            (MixedTensor::Param(a), MixedTensor::Param(b)) => {
                Some(MixedTensor::Param(a.add_fallible(b)?))
            }
            (MixedTensor::Param(s), MixedTensor::Concrete(o)) => match o {
                RealOrComplexTensor::Real(o) => Some(MixedTensor::Param(ParamTensor::composite(
                    s.tensor.add_fallible(o)?,
                ))),
                RealOrComplexTensor::Complex(o) => Some(MixedTensor::Param(
                    ParamTensor::composite(s.tensor.add_fallible(o)?),
                )),
            },
            (MixedTensor::Concrete(s), MixedTensor::Param(o)) => match s {
                RealOrComplexTensor::Real(s) => Some(MixedTensor::Param(ParamTensor::composite(
                    s.add_fallible(&o.tensor)?,
                ))),
                RealOrComplexTensor::Complex(s) => Some(MixedTensor::Param(
                    ParamTensor::composite(s.add_fallible(&o.tensor)?),
                )),
            },
            (MixedTensor::Concrete(s), MixedTensor::Concrete(o)) => {
                Some(MixedTensor::Concrete(s.add_fallible(o)?))
            }
        }
    }
}
