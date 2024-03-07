use ahash::{AHashMap, HashMap};
use enum_try_as_inner::EnumTryAsInner;

use num::Complex;
use smartstring::alias::String;
use symbolica::{
    domains::{
        float::NumericalFloatLike,
        rational::{Rational, RationalField},
    },
    evaluate::EvaluationFn,
    poly::{evaluate::InstructionEvaluator, polynomial::MultivariatePolynomial, Variable},
    representations::{AsAtomView, Atom, AtomView, FunctionBuilder, Symbol},
    state::{State, Workspace},
};

use super::{
    Contract, DataIterator, DataTensor, DenseTensor, HasName, HistoryStructure, SetTensorData,
    Slot, SparseTensor, StructureContract, SymbolicAdd, SymbolicAddAssign, SymbolicInto,
    SymbolicMul, SymbolicNeg, SymbolicStructureContract, SymbolicSub, SymbolicSubAssign,
    SymbolicZero, TensorStructure, TracksCount,
};

impl<'a, T, I> SparseTensor<T, I>
where
    for<'d> &'d T: SymbolicInto,
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> SparseTensor<Atom, I> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            let _ = result.set(&index, value.into_sym().unwrap());
        }
        result
    }
}

impl<'a, I> DenseTensor<Atom, I>
where
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    pub fn symbolic_zeros(structure: I) -> DenseTensor<Atom, I> {
        let result_data = vec![0; structure.size()];

        DenseTensor {
            data: result_data.iter().map(|&x| Atom::new_num(x)).collect(),
            structure,
        }
    }

    pub fn symbolic_labels(
        label: &str,
        structure: I,
        _ws: &'a Workspace,
        state: &'a mut State,
    ) -> DenseTensor<Atom, I> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let indices_str = index
                .into_iter()
                .map(|index| index.to_string().into())
                .collect::<Vec<String>>()
                .join("_");

            let value = Atom::parse(&format!("{}_{}", label, indices_str), state).unwrap();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn numbered_labeled(
        number: usize,
        label: Symbol,
        structure: I,
        _ws: &'a Workspace,
        _state: &'a State,
    ) -> DenseTensor<Atom, I> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let mut value_builder = FunctionBuilder::new(label);
            value_builder = value_builder.add_arg(Atom::new_num(number as i64).as_atom_view());

            for i in index {
                value_builder = value_builder.add_arg(Atom::new_num(i as i64).as_atom_view());
            }
            // Atom::parse(&format!("{}_{}_{}", label, indices_si)).unwrap();

            let value = value_builder.finish();

            data.push(value);
        }
        DenseTensor { data, structure }
    }
}

impl<'a, T, I> DenseTensor<T, I>
where
    for<'d> &'d T: SymbolicInto,
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> DenseTensor<Atom, I> {
        let mut result = DenseTensor::symbolic_zeros(self.structure.clone());
        for (index, value) in self.iter() {
            let _ = result.set(&index, value.into_sym().unwrap());
        }
        result
    }
}

pub trait FromStucture: Sized {
    fn from_structure(
        structure: HistoryStructure<String>,
        state: &mut State,
        ws: &Workspace,
    ) -> Option<Self>;
}

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum MixedTensor<T: TensorStructure> {
    Float(DataTensor<f64, T>),
    Complex(DataTensor<Complex<f64>, T>),
    Symbolic(DataTensor<Atom, T>),
}

impl<'a, I: TensorStructure + Clone + 'a> MixedTensor<I> {
    pub fn evaluate<'b>(&mut self, const_map: &'b HashMap<AtomView<'a>, f64>)
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

    pub fn append_const_map<'a, T>(
        &'a self,
        data: &DenseTensor<T, I>,
        const_map: &mut HashMap<AtomView<'a>, T>,
    ) where
        I: TensorStructure,
        T: Copy,
    {
        for ((i, a), (j, v)) in self.flat_iter().zip(data.flat_iter()) {
            assert_eq!(i, j);
            const_map.insert(a.as_view(), *v);
        }
    }

    pub fn to_evaluator<'a, N>(
        &'a self,
        _var_map: &mut HashMap<AtomView<'a>, Variable>,
        _state: &State,
    ) -> DenseTensor<InstructionEvaluator<N>, I>
    where
        N: NumericalFloatLike + for<'b> std::convert::From<&'b Rational>,
    {
        let structure = self.structure.clone();
        let data = self
            .data
            .iter()
            .map(|x| {
                let poly: MultivariatePolynomial<_, u8> = x
                    .as_view()
                    .to_polynomial_with_conversion(&RationalField::new());

                // x.as_view().evaluate(var_map, function_map, cache);

                // println!("{}", poly.printer(state));

                let (h, _ops, _) = poly.optimize_horner_scheme(4000);
                let mut i = h.to_instr(20);

                i.fuse_operations();

                for _ in 0..100_000 {
                    if !i.common_pair_elimination() {
                        break;
                    }
                    i.fuse_operations();
                }

                i.to_output(poly.var_map.as_ref().unwrap().to_vec(), true)
                    .convert::<N>()
                    .evaluator()
            })
            .collect::<Vec<_>>();

        DenseTensor { data, structure }
    }
}

impl<I, N> DenseTensor<InstructionEvaluator<N>, I>
where
    I: Clone,
    N: NumericalFloatLike + for<'a> std::convert::From<&'a Rational> + Copy,
{
    pub fn evaluate(&self, param: &[N]) -> DenseTensor<N, I> {
        let structure = self.structure.clone();
        let data = self
            .data
            .iter()
            .enumerate()
            .map(|(i, x)| x.clone().evaluate(&[param[i]])[0])
            .collect::<Vec<_>>();
        DenseTensor { data, structure }
    }
}

#[test]

fn test_evaluator() {
    let state = State::get_global_state().write().unwrap();
    let _ws = Workspace::new();
    let structure = crate::tensor::NamedStructure::from_integers(
        &[
            (crate::tensor::AbstractIndex(4), crate::tensor::Dimension(5)),
            (crate::tensor::AbstractIndex(5), crate::tensor::Dimension(4)),
        ],
        "r",
    );
    let p = crate::tensor::Shadowable::shadow(structure).unwrap();

    let mut var_map = AHashMap::new();
    let a: DenseTensor<InstructionEvaluator<f64>, crate::tensor::NamedStructure> =
        p.to_evaluator(&mut var_map, &state);

    // for (k, v) in var_map {
    //     println!("{} {:?}", k.printer(&state), v);
    // }

    let b = a.evaluate(&[
        1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0,
        4.0, 5.0,
    ]);

    println!("{:?}", b);
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
