use std::ops::AddAssign;

use symbolica::atom::{Atom, AtomView};

use crate::{
    algebra::complex::Complex,
    structure::{HasStructure, ScalarStructure, TensorStructure},
    tensors::complex::RealOrComplexTensor,
    tensors::data::DataTensor,
};

use super::{
    to_param::{ToAtom, ToParam},
    ConcreteOrParam, ConcreteOrParamRef, ParamOrComposite, ParamOrConcrete, ParamTensor,
};

// Scalars:

impl<T, U> AddAssign<ConcreteOrParamRef<'_, U>> for ConcreteOrParam<T>
where
    Atom: AddAssign<U> + for<'a> AddAssign<AtomView<'a>>,
    T: ToAtom + AddAssign<U>,
{
    fn add_assign(&mut self, rhs: ConcreteOrParamRef<'_, U>) {
        match (self, rhs) {
            (ConcreteOrParam::Concrete(c1), ConcreteOrParamRef::Concrete(c2)) => *c1 += c2,
            (ConcreteOrParam::Param(a), ConcreteOrParamRef::Concrete(c)) => *a += c,
            (ConcreteOrParam::Param(a1), ConcreteOrParamRef::Param(a2)) => *a1 += a2,
            (c, a) => {
                c.to_param();
                *c += a
            }
        }
    }
}

// impl<T: Into<Coefficient>> AddAssign<Complex<T>> for Atom {
//     fn add_assign(&mut self, rhs: Complex<T>) {
//         *self = &*self + rhs.re + Atom::i() * rhs.im;
//     }
// }

// impl<T> AddAssign<&Complex<T>> for Atom
// where
//     for<'a> &'a T: Into<Coefficient>,
// {
//     fn add_assign(&mut self, rhs: &Complex<T>) {
//         *self = &*self + &rhs.re + Atom::i() * &rhs.im;
//     }
// }

impl<T> AddAssign<ConcreteOrParam<T>> for Atom
where
    Atom: AddAssign<T> + AddAssign<Atom>,
{
    fn add_assign(&mut self, rhs: ConcreteOrParam<T>) {
        match rhs {
            ConcreteOrParam::Concrete(c) => *self += c,
            ConcreteOrParam::Param(a) => *self += a,
        }
    }
}

impl<T> AddAssign<ConcreteOrParamRef<'_, T>> for Atom
where
    Atom: AddAssign<T> + for<'a> AddAssign<AtomView<'a>>,
{
    fn add_assign(&mut self, rhs: ConcreteOrParamRef<'_, T>) {
        match rhs {
            ConcreteOrParamRef::Concrete(c) => *self += c,
            ConcreteOrParamRef::Param(a) => *self += a,
        }
    }
}

// impl<T> AddA/

// Tensors

impl<I, D> AddAssign<&RealOrComplexTensor<D, I>> for ParamTensor<I>
where
    DataTensor<Atom, I>:
        for<'a> AddAssign<&'a DataTensor<D, I>> + for<'a> AddAssign<&'a DataTensor<Complex<D>, I>>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &RealOrComplexTensor<D, I>) {
        match rhs {
            RealOrComplexTensor::Real(r) => self.tensor += r,
            RealOrComplexTensor::Complex(c) => self.tensor += c,
        }
        self.param_type = ParamOrComposite::Composite;
    }
}

impl<I, D> AddAssign<RealOrComplexTensor<D, I>> for ParamTensor<I>
where
    DataTensor<Atom, I>: AddAssign<DataTensor<D, I>> + AddAssign<DataTensor<Complex<D>, I>>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: RealOrComplexTensor<D, I>) {
        match rhs {
            RealOrComplexTensor::Real(r) => self.tensor += r,
            RealOrComplexTensor::Complex(c) => self.tensor += c,
        }
        self.param_type = ParamOrComposite::Composite;
    }
}

impl<I, D> AddAssign<DataTensor<D, I>> for ParamTensor<I>
where
    DataTensor<Atom, I>: AddAssign<DataTensor<D, I>>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: DataTensor<D, I>) {
        self.tensor += rhs;
        self.param_type = ParamOrComposite::Composite;
    }
}

impl<T, U, I> AddAssign<&ParamOrConcrete<T, I>> for ParamOrConcrete<U, I>
where
    U: for<'a> AddAssign<&'a T> + ToParam + HasStructure<Structure = I>,
    ParamTensor<I>: for<'a> AddAssign<&'a T> + for<'a> AddAssign<&'a ParamTensor<I>>,
    I: TensorStructure + Clone + ScalarStructure,
{
    fn add_assign(&mut self, rhs: &ParamOrConcrete<T, I>) {
        match (self, rhs) {
            (ParamOrConcrete::Param(p), ParamOrConcrete::Param(q)) => {
                *p += q;
            }
            (ParamOrConcrete::Param(p), ParamOrConcrete::Concrete(q)) => {
                *p += q;
            }
            (ParamOrConcrete::Concrete(p), ParamOrConcrete::Concrete(q)) => {
                *p += q;
            }
            (p, q) => {
                p.to_param();
                *p += q;
            }
        }
    }
}

impl<T, U, I> AddAssign<ParamOrConcrete<T, I>> for ParamOrConcrete<U, I>
where
    U: AddAssign<T> + ToParam + HasStructure<Structure = I>,
    ParamTensor<I>: AddAssign<T> + AddAssign<ParamTensor<I>>,
    I: TensorStructure + Clone + ScalarStructure,
{
    fn add_assign(&mut self, rhs: ParamOrConcrete<T, I>) {
        match (self, rhs) {
            (ParamOrConcrete::Param(p), ParamOrConcrete::Param(q)) => {
                *p += q;
            }
            (ParamOrConcrete::Param(p), ParamOrConcrete::Concrete(q)) => {
                *p += q;
            }
            (ParamOrConcrete::Concrete(p), ParamOrConcrete::Concrete(q)) => {
                *p += q;
            }
            (p, q) => {
                p.to_param();
                *p += q;
            }
        }
    }
}

impl<I> AddAssign<&ParamTensor<I>> for ParamTensor<I>
where
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &ParamTensor<I>) {
        self.tensor += &rhs.tensor;

        self.param_type = ParamOrComposite::Composite;
    }
}

impl<I> AddAssign<ParamTensor<I>> for ParamTensor<I>
where
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: ParamTensor<I>) {
        self.tensor += rhs.tensor;
        self.param_type = ParamOrComposite::Composite;
    }
}

impl<I, D> AddAssign<&DataTensor<D, I>> for ParamTensor<I>
where
    DataTensor<Atom, I>: for<'a> AddAssign<&'a DataTensor<D, I>>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &DataTensor<D, I>) {
        self.tensor += rhs;
        self.param_type = ParamOrComposite::Composite;
    }
}
