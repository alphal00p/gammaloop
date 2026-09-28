use symbolica::atom::Atom;

use crate::{
    algebra::complex::{Complex, RealOrComplex},
    structure::{HasStructure, ScalarStructure, TensorStructure},
    tensors::complex::RealOrComplexTensor,
    tensors::data::{DataTensor, SparseTensor, StorageTensor},
};

use super::{ParamOrComposite, ParamOrConcrete, ParamTensor};

pub trait ToParam: HasStructure {
    fn to_param(self) -> ParamTensor<Self::Structure>;
}

pub trait ToAtom {
    fn to_atom(self) -> Atom;
}

impl ToAtom for f64 {
    fn to_atom(self) -> Atom {
        Atom::num(self)
    }
}

impl ToAtom for Complex<f64> {
    fn to_atom(self) -> Atom {
        Atom::num(self.re) + Atom::num(self.im) * Atom::i()
    }
}

impl<T> ToAtom for RealOrComplex<T>
where
    T: ToAtom,
{
    fn to_atom(self) -> Atom {
        match self {
            RealOrComplex::Real(r) => r.to_atom(),
            RealOrComplex::Complex(c) => c.re.to_atom() + c.im.to_atom() * Atom::i(),
        }
    }
}

impl<I: TensorStructure + Clone, D: ToAtom> ToParam for DataTensor<D, I> {
    fn to_param(self) -> ParamTensor<Self::Structure> {
        ParamTensor {
            tensor: self.map_data(ToAtom::to_atom),
            param_type: ParamOrComposite::Composite,
        }
    }
}

impl<I: TensorStructure + Clone, D: ToAtom + Clone> ToParam for RealOrComplexTensor<D, I> {
    fn to_param(self) -> ParamTensor<Self::Structure> {
        ParamTensor {
            tensor: match self {
                RealOrComplexTensor::Real(r) => r.map_data(ToAtom::to_atom),
                RealOrComplexTensor::Complex(c) => {
                    c.map_data(|a| a.re.to_atom() + a.im.to_atom() * Atom::i())
                }
            },
            param_type: ParamOrComposite::Composite,
        }
    }
}

impl<I: TensorStructure> ToParam for ParamTensor<I> {
    fn to_param(self) -> ParamTensor<Self::Structure> {
        self
    }
}

impl<T: ToParam + HasStructure<Structure = I>, I: ScalarStructure + TensorStructure>
    ParamOrConcrete<T, I>
{
    pub fn to_param(&mut self) {
        if matches!(self, ParamOrConcrete::Concrete(_)) {
            let old = std::mem::replace(
                self,
                ParamOrConcrete::Param(ParamTensor::param(DataTensor::Sparse(
                    SparseTensor::empty(I::scalar_structure(), Atom::Zero),
                ))),
            );

            if let ParamOrConcrete::Concrete(r) = old {
                *self = ParamOrConcrete::Param(r.to_param());
            }
        }
    }
}
