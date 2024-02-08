use crate::tensor::IntoId;
use ahash::AHashMap;
use enum_dispatch::enum_dispatch;
use enum_try_as_inner::EnumTryAsInner;
use indexmap::IndexMap;
use num::Complex;
use smartstring::alias::String;
use symbolica::{
    representations::{AsAtomView, Atom, FunctionBuilder, Identifier},
    state::{State, Workspace},
};

use rustc_hash::FxHasher;
type FxIndexMap<K, V> = IndexMap<K, V, FxHasher>;
use nohash_hasher::BuildNoHashHasher;
type NHIndexMap<K, V> = IndexMap<K, V, BuildNoHashHasher<K>>;
use super::{
    ConcreteIndex, DenseTensor, Expr, HasTensorStructure, NumTensor, SparseTensor, SymbolicAdd,
    SymbolicAddAssign, SymbolicInto, SymbolicMul, SymbolicNeg, SymbolicSub, SymbolicSubAssign,
    SymbolicZero, TensorSkeleton,
};
use intmap::IntMap;

pub trait SymbolicInternalContract {
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self;
}

impl<T, I> SymbolicInternalContract for DenseTensor<T, I>
where
    T: for<'a> SymbolicAddAssign<&'a T>
        + for<'b> SymbolicSubAssign<&'b T>
        + SymbolicNeg
        + Clone
        + std::fmt::Debug,
    I: Clone,
{
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self {
        let mut result: DenseTensor<T, I> = self.clone();
        for trace in self.traces() {
            let mut new_structure = self.structure().clone();
            new_structure.trace(trace[0], trace[1]);

            let mut new_result: DenseTensor<T, I> =
                DenseTensor::from_data_coerced(&self.data, new_structure).unwrap();
            for (idx, t) in result.iter_symbolic_trace(trace, state, ws) {
                new_result.set(&idx, t);
            }
            result = new_result;
        }
        result
    }
}

impl<T, I> SymbolicInternalContract for SparseTensor<T, I>
where
    T: for<'a> SymbolicAddAssign<&'a T>
        + for<'b> SymbolicSubAssign<&'b T>
        + SymbolicNeg
        + Clone
        + std::fmt::Debug,
    I: Clone,
{
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self {
        let trace = self.traces()[0];

        // println!("trace {:?}", trace);
        let mut new_structure = self.structure().clone();
        new_structure.trace(trace[0], trace[1]);

        let mut new_result = SparseTensor::empty(new_structure);
        for (idx, t) in self.iter_symbolic_trace(trace, state, ws)
        // .filter(|(_, t)| *t != T::default())
        {
            new_result.set(&idx, t).unwrap();
        }

        if new_result.traces().is_empty() {
            new_result
        } else {
            new_result.internal_contract_sym(ws, state)
        }
    }
}

pub trait SymbolicContract<T> {
    type LCM;
    fn contract_sym(&self, other: &T, state: &State, ws: &Workspace) -> Option<Self::LCM>;
}

impl<T, U, I, Out> SymbolicContract<DenseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: SymbolicAdd<&'b T, Output = Out>
        + SymbolicMul<&'b T, Output = Out>
        + SymbolicSub<&'b T, Output = Out>,
    for<'a, 'b> &'a Out: SymbolicAdd<&'b Out, Output = Out> + SymbolicSub<&'b Out, Output = Out>,
    Out: SymbolicZero
        + Clone
        + std::fmt::Debug
        + for<'a> SymbolicAddAssign<&'a Out>
        + SymbolicNeg
        + for<'b> SymbolicSubAssign<&'b Out>,
    I: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((_, i, j)) = self.structure().match_index(other.structure()) {
            let dimension_of_contraction = self.shape()[i];
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let metric = self.structure()[i].representation.negative();
            let mut result_index = 0;
            for fiber_a in self.iter_fibers(i) {
                for fiber_b in other.iter_fibers(j) {
                    for i in 0..dimension_of_contraction {
                        if metric[i] {
                            result_data[result_index] = result_data[result_index].sub_sym(
                                &fiber_a[i].mul_sym(fiber_b[i], ws, state)?,
                                ws,
                                state,
                            )?;
                        } else {
                            result_data[result_index] = result_data[result_index].add_sym(
                                &fiber_a[i].mul_sym(fiber_b[i], ws, state)?,
                                ws,
                                state,
                            )?;
                        }
                    }
                    result_index += 1;
                }
            }

            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract_sym(ws, state));
            }
        }
        None
    }
}

impl<T, U, I, Out> SymbolicContract<DenseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: SymbolicAdd<&'b T, Output = Out>
        + SymbolicMul<&'b T, Output = Out>
        + SymbolicSub<&'b T, Output = Out>,
    for<'a, 'b> &'a Out: SymbolicAdd<&'b Out, Output = Out> + SymbolicSub<&'b Out, Output = Out>,
    Out: SymbolicZero
        + Clone
        + std::fmt::Debug
        + for<'a> SymbolicAddAssign<&'a Out>
        + SymbolicNeg
        + for<'b> SymbolicSubAssign<&'b Out>,
    I: Clone,
    U: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((_, i, j)) = self.structure().match_index(other.structure()) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let metric = self.get_ith_metric(i).unwrap();
            let mut result_index = 0;
            for (skipped, nonzeros, fiber_a) in self.iter_fibers(i) {
                result_index += skipped;
                for fiber_b in other.iter_fibers(j) {
                    for (i, k) in nonzeros.iter().enumerate() {
                        if metric[*k] {
                            result_data[result_index] = result_data[result_index].sub_sym(
                                &fiber_a[i].mul_sym(fiber_b[*k], ws, state)?,
                                ws,
                                state,
                            )?;
                        } else {
                            result_data[result_index] = result_data[result_index].add_sym(
                                &fiber_a[i].mul_sym(fiber_b[*k], ws, state)?,
                                ws,
                                state,
                            )?;
                        }
                    }
                    result_index += 1;
                }
            }

            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract_sym(ws, state));
            }
        }
        None
    }
}

#[allow(clippy::if_same_then_else)]
impl<T, U, Out, I> SymbolicContract<SparseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: SymbolicAdd<&'b T, Output = Out>
        + SymbolicMul<&'b T, Output = Out>
        + SymbolicSub<&'b T, Output = Out>,
    for<'a, 'b> &'a Out: SymbolicAdd<&'b Out, Output = Out> + SymbolicSub<&'b Out, Output = Out>,
    Out: SymbolicZero
        + Clone
        + std::fmt::Debug
        + for<'a> SymbolicAddAssign<&'a Out>
        + SymbolicNeg
        + for<'b> SymbolicSubAssign<&'b Out>,
    I: Clone,
    T: Clone,
    U: Clone,
{
    type LCM = SparseTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &SparseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((_, i, j)) = self.structure().match_index(other.structure()) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = AHashMap::default();
            let metric = self.get_ith_metric(i).unwrap();
            let one = 1;
            let stride_other = *final_structure
                .strides()
                .get(self.structure().order() - 2)
                .unwrap_or(&one);

            let mut result_index = 0;
            for (skipped_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                result_index += skipped_a;
                for (skipped_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {
                    result_index += skipped_b * stride_other;
                    let mut value = Out::zero(state, ws);
                    let mut nonzero = false;
                    for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                        nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)) // Only store the positions
                    }) {
                        nonzero = true;
                        if metric[x] {
                            value = value.sub_sym(
                                &fiber_a[i].mul_sym(fiber_b[j], ws, state)?,
                                ws,
                                state,
                            )?;
                        } else {
                            value = value.add_sym(
                                &fiber_a[i].mul_sym(fiber_b[j], ws, state)?,
                                ws,
                                state,
                            )?;
                        }
                    }

                    if nonzero
                    // && value.as_atom_view() != zero.as_atom_view() TODO impl generic
                    {
                        result_data.insert(result_index, value);
                    }
                    result_index += 1;
                }
            }

            let result = SparseTensor {
                elements: result_data,
                structure: final_structure,
            };

            if result.traces().is_empty() {
                return Some(result);
            } else {
                return Some(result.internal_contract_sym(ws, state)); //.internal_contract() need to impl for symbolic
            }
        }
        None
    }
}

impl<T, U, I, Out> SymbolicContract<SparseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a T: SymbolicAdd<&'b U, Output = Out>
        + SymbolicMul<&'b U, Output = Out>
        + SymbolicSub<&'b U, Output = Out>,
    for<'a, 'b> &'a Out: SymbolicAdd<&'b Out, Output = Out> + SymbolicSub<&'b Out, Output = Out>,
    Out: SymbolicZero
        + Clone
        + SymbolicAdd<Out, Output = Out>
        + SymbolicSub<Out, Output = Out>
        + std::fmt::Debug
        + for<'a> SymbolicAddAssign<&'a Out>
        + SymbolicNeg
        + for<'b> SymbolicSubAssign<&'b Out>,
    I: Clone,
    T: Clone,
{
    type LCM = DenseTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &SparseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        other.contract_sym(self, state, ws)
    }
}

impl<I> SparseTensor<Atom, I>
where
    I: Clone,
{
    pub fn builder<'a>(&self, state: &'a State, ws: &'a Workspace) -> SparseTensor<Expr<'a>, I> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(&index, value.builder(state, ws)).unwrap();
        }
        result
    }
}

impl<I> DenseTensor<Atom, I>
where
    I: Clone,
{
    pub fn builder<'a>(&self, state: &'a State, ws: &'a Workspace) -> DenseTensor<Expr<'a>, I> {
        let mut result = DenseTensor::neutral_builder(self.structure.clone(), ws, state);
        for (index, value) in self.iter() {
            result.set(&index, value.builder(state, ws));
        }
        result
    }
}

// impl<T: ConvertableToSymbolic> SparseTensor<T> {
//     pub fn to_symbolic<'a>(&self, ws: &'a Workspace, state: &'a State) -> SparseTensor<Atom> {
//         let mut result = SparseTensor::empty(self.structure.clone());
//         for (index, value) in self.iter() {
//             result
//                 .set(index, value.to_symbolic(ws, state).unwrap())
//                 .unwrap();
//         }
//         result
//     }

//     pub fn to_symbolic_builder<'a>(
//         &self,
//         ws: &'a Workspace,
//         state: &'a State,
//     ) -> SparseTensor<Expr<'a>> {
//         let mut result = SparseTensor::empty(self.structure.clone());
//         for (index, value) in self.iter() {
//             result
//                 .set(
//                     index,
//                     value.to_symbolic(ws, state).unwrap().builder(state, ws),
//                 )
//                 .unwrap();
//         }
//         result
//     }
// }

impl<'a, T, I> SparseTensor<T, I>
where
    for<'d> &'d T: SymbolicInto,
    I: Clone,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> SparseTensor<Atom, I> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            let _ = result.set(&index, value.into_sym(ws, state).unwrap());
        }
        result
    }

    pub fn to_symbolic_builder<'b>(
        &'b self,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> SparseTensor<Expr<'a>, I> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            let _ = result.set(
                &index,
                value.into_sym(ws, state).unwrap().builder(state, ws),
            );
        }
        result
    }
}

impl<'a, I> DenseTensor<Expr<'a>, I> {
    pub fn symbolic_zeros(structure: TensorSkeleton<I>) -> DenseTensor<Atom, I> {
        let result_data = vec![0; structure.size()];

        DenseTensor {
            data: result_data.iter().map(|&x| Atom::new_num(x)).collect(),
            structure,
        }
    }

    pub fn neutral_builder(
        structure: TensorSkeleton<I>,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>, I> {
        let zero = Atom::new_num(0).builder(state, ws);
        let result_data = vec![zero; structure.size()];
        DenseTensor {
            data: result_data,
            structure,
        }
    }
    pub fn symbolic_labels(
        label: &str,
        structure: TensorSkeleton<I>,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> DenseTensor<Atom, I> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let indices_str = index
                .into_iter()
                .map(|index| index.to_string().into())
                .collect::<Vec<String>>()
                .join("_");

            let value = Atom::parse(&format!("{}_{}", label, indices_str), state, ws).unwrap();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn numbered_labeled(
        number: usize,
        label: Identifier,
        structure: TensorSkeleton<I>,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Atom, I> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let mut value_builder = FunctionBuilder::new(label, state, ws);
            value_builder = value_builder.add_arg(Atom::new_num(number as i64).as_atom_view());

            for i in index {
                value_builder = value_builder.add_arg(Atom::new_num(i as i64).as_atom_view());
            }
            // Atom::parse(&format!("{}_{}_{}", label, indices_str, i), state, ws).unwrap();

            let value = value_builder.finish().into_atom();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn numbered_labeled_builder(
        number: usize,
        label: Identifier,
        structure: TensorSkeleton<I>,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>, I> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let mut value_builder = FunctionBuilder::new(label, state, ws);
            value_builder = value_builder.add_arg(Atom::new_num(number as i64).as_atom_view());

            for i in index {
                value_builder = value_builder.add_arg(Atom::new_num(i as i64).as_atom_view());
            }
            // Atom::parse(&format!("{}_{}_{}", label, indices_str, i), state, ws).unwrap();

            let value = value_builder.finish();

            data.push(value);
        }
        DenseTensor { data, structure }
    }

    pub fn finish(self) -> DenseTensor<Atom, I> {
        DenseTensor {
            data: self.data.into_iter().map(|x| x.into_atom()).collect(),
            structure: self.structure,
        }
    }
}

impl<'a, T, I> DenseTensor<T, I>
where
    for<'d> &'d T: SymbolicInto,
    I: Clone,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> DenseTensor<Atom, I> {
        let mut result = DenseTensor::symbolic_zeros(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(&index, value.into_sym(ws, state).unwrap());
        }
        result
    }

    pub fn to_symbolic_builder<'b>(
        &'b self,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> DenseTensor<Expr<'a>, I> {
        let mut result = DenseTensor::neutral_builder(self.structure.clone(), ws, state);
        for (index, value) in self.iter() {
            result.set(
                &index,
                value.into_sym(ws, state).unwrap().builder(state, ws),
            );
        }
        result
    }
}

impl<T, U, I, Out> SymbolicContract<NumTensor<T, I>> for NumTensor<U, I>
where
    for<'a, 'b> &'a U: SymbolicAdd<&'b T, Output = Out>
        + SymbolicMul<&'b T, Output = Out>
        + SymbolicSub<&'b T, Output = Out>,
    for<'a, 'b> &'a T: SymbolicAdd<&'b U, Output = Out>
        + SymbolicMul<&'b U, Output = Out>
        + SymbolicSub<&'b U, Output = Out>,
    for<'a, 'b> &'a Out: SymbolicAdd<&'b Out, Output = Out> + SymbolicSub<&'b Out, Output = Out>,
    Out: SymbolicZero
        + Clone
        + std::fmt::Debug
        + for<'a> SymbolicAddAssign<&'a Out>
        + SymbolicNeg
        + for<'b> SymbolicSubAssign<&'b Out>,
    I: Clone,
    T: Clone,
    U: Clone,
{
    type LCM = NumTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &NumTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        match (self, other) {
            (NumTensor::Dense(s), NumTensor::Dense(o)) => {
                Some(NumTensor::Dense(s.contract_sym(o, state, ws)?))
            }
            (NumTensor::Dense(s), NumTensor::Sparse(o)) => {
                Some(NumTensor::Dense(o.contract_sym(s, state, ws)?))
            }
            (NumTensor::Sparse(s), NumTensor::Dense(o)) => {
                Some(NumTensor::Dense(s.contract_sym(o, state, ws)?))
            }
            (NumTensor::Sparse(s), NumTensor::Sparse(o)) => {
                Some(NumTensor::Sparse(s.contract_sym(o, state, ws)?))
            }
        }
    }
}

pub trait FromStucture: Sized {
    fn from_structure(
        structure: TensorSkeleton<String>,
        state: &mut State,
        ws: &Workspace,
    ) -> Option<Self>;
}

#[derive(Debug, Clone, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum MixedTensor<T> {
    Float(NumTensor<f64, T>),
    Complex(NumTensor<Complex<f64>, T>),
    Symbolic(NumTensor<Atom, T>),
}

impl<I> HasTensorStructure for MixedTensor<I>
where
    I: Clone,
{
    type Name = I;
    fn structure(&self) -> &TensorSkeleton<Self::Name> {
        match self {
            MixedTensor::Float(t) => t.structure(),
            MixedTensor::Complex(t) => t.structure(),
            MixedTensor::Symbolic(t) => t.structure(),
        }
    }

    fn mut_structure(&mut self) -> &mut TensorSkeleton<Self::Name> {
        match self {
            MixedTensor::Float(t) => t.mut_structure(),
            MixedTensor::Complex(t) => t.mut_structure(),
            MixedTensor::Symbolic(t) => t.mut_structure(),
        }
    }
}

pub type MixedTensors = MixedTensor<Identifier>;

impl<I> From<DenseTensor<f64, I>> for MixedTensor<I> {
    fn from(other: DenseTensor<f64, I>) -> Self {
        MixedTensor::<I>::Float(NumTensor::Dense(other))
    }
}

impl<I> From<SparseTensor<f64, I>> for MixedTensor<I> {
    fn from(other: SparseTensor<f64, I>) -> Self {
        MixedTensor::<I>::Float(NumTensor::Sparse(other))
    }
}

impl<I> From<DenseTensor<Complex<f64>, I>> for MixedTensor<I> {
    fn from(other: DenseTensor<Complex<f64>, I>) -> Self {
        MixedTensor::<I>::Complex(NumTensor::Dense(other))
    }
}

impl<I> From<SparseTensor<Complex<f64>, I>> for MixedTensor<I> {
    fn from(other: SparseTensor<Complex<f64>, I>) -> Self {
        MixedTensor::<I>::Complex(NumTensor::Sparse(other))
    }
}

impl<I> From<DenseTensor<Atom, I>> for MixedTensor<I> {
    fn from(other: DenseTensor<Atom, I>) -> Self {
        MixedTensor::<I>::Symbolic(NumTensor::Dense(other))
    }
}

impl<I> From<SparseTensor<Atom, I>> for MixedTensor<I> {
    fn from(other: SparseTensor<Atom, I>) -> Self {
        MixedTensor::<I>::Symbolic(NumTensor::Sparse(other))
    }
}

impl<I> SymbolicContract<MixedTensor<I>> for MixedTensor<I>
where
    I: Clone,
{
    type LCM = MixedTensor<I>;
    fn contract_sym(
        &self,
        other: &MixedTensor<I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        match (self, other) {
            (MixedTensor::<I>::Float(s), MixedTensor::<I>::Float(o)) => {
                Some(MixedTensor::<I>::Float(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Float(s), MixedTensor::<I>::Complex(o)) => {
                Some(MixedTensor::<I>::Complex(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Float(s), MixedTensor::<I>::Symbolic(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Complex(s), MixedTensor::<I>::Float(o)) => {
                Some(MixedTensor::<I>::Complex(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Complex(s), MixedTensor::<I>::Complex(o)) => {
                Some(MixedTensor::<I>::Complex(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Complex(s), MixedTensor::<I>::Symbolic(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Symbolic(s), MixedTensor::<I>::Float(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Symbolic(s), MixedTensor::<I>::Complex(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensor::<I>::Symbolic(s), MixedTensor::<I>::Symbolic(o)) => {
                Some(MixedTensor::<I>::Symbolic(s.contract_sym(o, state, ws)?))
            }
        }
    }
}
