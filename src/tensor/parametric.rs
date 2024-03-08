use ahash::{AHashMap, HashMap};
use enum_try_as_inner::EnumTryAsInner;

use symbolica::{
    domains::float::Complex,
    evaluate::EvaluationFn,
    representations::{Atom, AtomView, Symbol},
};

use super::{
    Contract, DataIterator, DataTensor, DenseTensor, HasName, HistoryStructure, Slot, SparseTensor,
    StructureContract, TensorStructure, TracksCount, VecStructure,
};

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum MixedTensor<T: TensorStructure = VecStructure> {
    Float(DataTensor<f64, T>),
    Complex(DataTensor<Complex<f64>, T>),
    Symbolic(DataTensor<Atom, T>),
}

impl<'a, I: TensorStructure + Clone + 'a> MixedTensor<I> {
    pub fn evaluate_float<'b>(&mut self, const_map: &'b HashMap<AtomView<'a>, f64>)
    where
        'b: 'a,
    {
        let content = match self {
            MixedTensor::Symbolic(x) => Some(x),
            _ => None,
        };

        if let Some(x) = content {
            *self = MixedTensor::Float(x.evaluate(const_map));
        }
    }

    pub fn evaluate_complex<'b>(&mut self, const_map: &'b HashMap<AtomView<'a>, Complex<f64>>)
    where
        'b: 'a,
    {
        let content = match self {
            MixedTensor::Symbolic(x) => Some(x),
            _ => None,
        };

        if let Some(x) = content {
            *self = MixedTensor::Complex(x.evaluate(const_map));
        }
    }
}

impl<I> DataTensor<Atom, I>
where
    I: Clone + TensorStructure,
{
    pub fn evaluate<'a, 'b, T>(&self, const_map: &'b HashMap<AtomView<'a>, T>) -> DataTensor<T, I>
    where
        T: symbolica::domains::float::Real
            + for<'c> std::convert::From<&'c symbolica::domains::rational::Rational>,
        'a: 'b,
    {
        match self {
            DataTensor::Dense(x) => DataTensor::Dense(x.evaluate(const_map)),
            DataTensor::Sparse(x) => DataTensor::Sparse(x.evaluate(const_map)),
        }
    }
}

impl<I> SparseTensor<Atom, I>
where
    I: Clone,
{
    pub fn evaluate<'a, T>(&self, const_map: &HashMap<AtomView<'a>, T>) -> SparseTensor<T, I>
    where
        T: symbolica::domains::float::Real
            + for<'d> std::convert::From<&'d symbolica::domains::rational::Rational>,
    {
        let fn_map: HashMap<_, EvaluationFn<_>> = HashMap::default();
        let mut cache = HashMap::default();
        let structure = self.structure.clone();
        let data = self
            .elements
            .iter()
            .map(|(idx, x)| {
                (
                    *idx,
                    x.as_view().evaluate::<T>(const_map, &fn_map, &mut cache),
                )
            })
            .collect::<AHashMap<_, _>>();

        SparseTensor {
            elements: data,
            structure,
        }
    }
}

impl<I> DenseTensor<Atom, I>
where
    I: Clone,
{
    pub fn evaluate<'a, T>(&'a self, const_map: &HashMap<AtomView<'a>, T>) -> DenseTensor<T, I>
    where
        T: symbolica::domains::float::Real
            + for<'b> std::convert::From<&'b symbolica::domains::rational::Rational>,
    {
        let fn_map: HashMap<_, EvaluationFn<_>> = HashMap::default();
        let mut cache = HashMap::default();
        let structure = self.structure.clone();
        let data = self
            .data
            .iter()
            .map(|x| x.as_view().evaluate::<T>(const_map, &fn_map, &mut cache))
            .collect::<Vec<_>>();

        DenseTensor { data, structure }
    }

    pub fn append_const_map<'a, 'b, T>(
        &'a self,
        data: &DenseTensor<T, I>,
        const_map: &mut HashMap<AtomView<'b>, T>,
    ) where
        I: TensorStructure,
        T: Copy,
        'a: 'b,
    {
        for ((i, a), (j, v)) in self.flat_iter().zip(data.flat_iter()) {
            assert_eq!(i, j);
            const_map.insert(a.as_view(), *v);
        }
    }
}

impl<T> TensorStructure for MixedTensor<T>
where
    T: TensorStructure,
{
    type Structure = T;

    fn structure(&self) -> &Self::Structure {
        match self {
            MixedTensor::Float(t) => t.structure(),
            MixedTensor::Complex(t) => t.structure(),
            MixedTensor::Symbolic(t) => t.structure(),
        }
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        match self {
            MixedTensor::Float(t) => t.mut_structure(),
            MixedTensor::Complex(t) => t.mut_structure(),
            MixedTensor::Symbolic(t) => t.mut_structure(),
        }
    }
    fn external_structure(&self) -> &[Slot] {
        match self {
            MixedTensor::Float(t) => t.external_structure(),
            MixedTensor::Complex(t) => t.external_structure(),
            MixedTensor::Symbolic(t) => t.external_structure(),
        }
    }
}

impl<T> HasName for MixedTensor<T>
where
    T: HasName + TensorStructure,
{
    type Name = T::Name;

    fn name(&self) -> Option<std::borrow::Cow<'_, <T as HasName>::Name>> {
        match self {
            MixedTensor::Float(t) => t.name(),
            MixedTensor::Complex(t) => t.name(),
            MixedTensor::Symbolic(t) => t.name(),
        }
    }

    fn set_name(&mut self, name: &Self::Name) {
        match self {
            MixedTensor::Float(t) => t.set_name(name),
            MixedTensor::Complex(t) => t.set_name(name),
            MixedTensor::Symbolic(t) => t.set_name(name),
        }
    }
}

impl<T> TracksCount for MixedTensor<T>
where
    T: TracksCount + TensorStructure,
{
    fn contractions_num(&self) -> usize {
        match self {
            MixedTensor::Float(t) => t.contractions_num(),
            MixedTensor::Complex(t) => t.contractions_num(),
            MixedTensor::Symbolic(t) => t.contractions_num(),
        }
    }
}

pub type MixedTensors = MixedTensor<HistoryStructure<Symbol>>;

impl<I> From<DenseTensor<f64, I>> for MixedTensor<I>
where
    I: TensorStructure,
{
    fn from(other: DenseTensor<f64, I>) -> Self {
        MixedTensor::<I>::Float(DataTensor::Dense(other))
    }
}

impl<I> From<SparseTensor<f64, I>> for MixedTensor<I>
where
    I: TensorStructure,
{
    fn from(other: SparseTensor<f64, I>) -> Self {
        MixedTensor::<I>::Float(DataTensor::Sparse(other))
    }
}

impl<I> From<DenseTensor<Complex<f64>, I>> for MixedTensor<I>
where
    I: TensorStructure,
{
    fn from(other: DenseTensor<Complex<f64>, I>) -> Self {
        MixedTensor::<I>::Complex(DataTensor::Dense(other))
    }
}

impl<I> From<SparseTensor<Complex<f64>, I>> for MixedTensor<I>
where
    I: TensorStructure,
{
    fn from(other: SparseTensor<Complex<f64>, I>) -> Self {
        MixedTensor::<I>::Complex(DataTensor::Sparse(other))
    }
}

impl<I> From<DenseTensor<Atom, I>> for MixedTensor<I>
where
    I: TensorStructure,
{
    fn from(other: DenseTensor<Atom, I>) -> Self {
        MixedTensor::<I>::Symbolic(DataTensor::Dense(other))
    }
}

impl<I> From<SparseTensor<Atom, I>> for MixedTensor<I>
where
    I: TensorStructure,
{
    fn from(other: SparseTensor<Atom, I>) -> Self {
        MixedTensor::<I>::Symbolic(DataTensor::Sparse(other))
    }
}

impl<I> Contract<MixedTensor<I>> for MixedTensor<I>
where
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = MixedTensor<I>;
    fn contract(&self, other: &MixedTensor<I>) -> Option<Self::LCM> {
        match (self, other) {
            (MixedTensor::<I>::Float(s), MixedTensor::<I>::Float(o)) => {
                Some(MixedTensor::<I>::Float(s.contract(o)?))
            }
            (MixedTensor::<I>::Float(s), MixedTensor::<I>::Complex(o)) => {
                Some(MixedTensor::<I>::Complex(s.contract(o)?))
            }
            (MixedTensor::<I>::Float(s), MixedTensor::<I>::Symbolic(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract(o)?))
            }
            (MixedTensor::<I>::Complex(s), MixedTensor::<I>::Float(o)) => {
                Some(MixedTensor::<I>::Complex(s.contract(o)?))
            }
            (MixedTensor::<I>::Complex(s), MixedTensor::<I>::Complex(o)) => {
                Some(MixedTensor::<I>::Complex(s.contract(o)?))
            }
            (MixedTensor::<I>::Complex(s), MixedTensor::<I>::Symbolic(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract(o)?))
            }
            (MixedTensor::<I>::Symbolic(s), MixedTensor::<I>::Float(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract(o)?))
            }
            (MixedTensor::<I>::Symbolic(s), MixedTensor::<I>::Complex(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract(o)?))
            }
            (MixedTensor::<I>::Symbolic(s), MixedTensor::<I>::Symbolic(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract(o)?))
            }
        }
    }
}
