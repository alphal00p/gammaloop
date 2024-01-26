use ahash::AHashMap;
use enum_dispatch::enum_dispatch;
use hyperdual::Num;
use num::Complex;
use rustc_hash::FxHashMap;
use std::{collections::HashMap, process::Output};
use symbolica::{
    representations::{
        default::Linear, number::Number, AsAtomView, Atom, FunctionBuilder, Identifier,
    },
    state::{ResettableBuffer, State, Workspace},
};

use super::{
    ConcreteIndex, Contract, DenseTensor, Expr, HasTensorStructure, NumTensor, SmallestUpgrade,
    SparseTensor, SymbolicAdd, SymbolicAddAssign, SymbolicInto, SymbolicMul, SymbolicNeg,
    SymbolicSub, SymbolicSubAssign, SymbolicZero, TensorSkeleton, VecSlotExtension,
};

pub trait SymbolicInternalContract {
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self;
}

impl<T> SymbolicInternalContract for DenseTensor<T>
where
    T: for<'a> SymbolicAddAssign<&'a T>
        + for<'b> SymbolicSubAssign<&'b T>
        + SymbolicNeg
        + Clone
        + std::fmt::Debug,
{
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self {
        let mut result = self.clone();
        for trace in self.traces() {
            let mut new_structure = self.structure().clone();
            new_structure.trace(trace[0], trace[1]);

            let mut new_result = DenseTensor::from_data_coerced(&self.data, new_structure).unwrap();
            for (idx, t) in result.iter_symbolic_trace(trace, state, ws) {
                new_result.set(&idx, t);
            }
            result = new_result;
        }
        result
    }
}

impl<T> SymbolicInternalContract for SparseTensor<T>
where
    T: for<'a> SymbolicAddAssign<&'a T>
        + for<'b> SymbolicSubAssign<&'b T>
        + SymbolicNeg
        + Clone
        + std::fmt::Debug,
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

impl<T, U, Out> SymbolicContract<DenseTensor<T>> for DenseTensor<U>
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
{
    type LCM = DenseTensor<Out>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            let dimension_of_contraction = self.shape()[i];
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let metric = self.structure()[i].representation.negative();
            for (index_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();
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

impl<T, U, Out> SymbolicContract<DenseTensor<T>> for SparseTensor<U>
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
{
    type LCM = DenseTensor<Out>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let metric = self.get_ith_metric(i).unwrap();
            for (index_a, nonzeros, fiber_a) in self.iter_fibers(i) {
                for (index_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = final_structure
                        .flat_index(
                            &index_a[..i]
                                .iter()
                                .chain(&index_a[i + 1..])
                                .chain(&index_b[..j])
                                .chain(&index_b[j + 1..])
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .unwrap();
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
impl<T, U, Out> SymbolicContract<SparseTensor<T>> for SparseTensor<U>
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
{
    type LCM = SparseTensor<Out>;
    fn contract_sym(
        &self,
        other: &SparseTensor<T>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = AHashMap::new();
            let metric = self.get_ith_metric(i).unwrap();
            for (index_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = index_a[..i]
                        .iter()
                        .chain(&index_a[i + 1..])
                        .chain(&index_b[..j])
                        .chain(&index_b[j + 1..])
                        .cloned()
                        .collect::<Vec<_>>();
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

impl<T, U, Out> SymbolicContract<SparseTensor<T>> for DenseTensor<U>
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
{
    type LCM = DenseTensor<Out>;
    fn contract_sym(
        &self,
        other: &SparseTensor<T>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        other.contract_sym(self, state, ws)
    }
}

impl SparseTensor<Atom> {
    pub fn builder<'a>(&self, state: &'a State, ws: &'a Workspace) -> SparseTensor<Expr<'a>> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(index, value.builder(state, ws)).unwrap();
        }
        result
    }
}

impl DenseTensor<Atom> {
    pub fn builder<'a>(&self, state: &'a State, ws: &'a Workspace) -> DenseTensor<Expr<'a>> {
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

impl<'a, T> SparseTensor<T>
where
    for<'d> &'d T: SymbolicInto,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> SparseTensor<Atom> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(&index, value.into_sym(ws, state).unwrap());
        }
        result
    }

    pub fn to_symbolic_builder<'b>(
        &'b self,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> SparseTensor<Expr<'a>> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result.set(
                &index,
                value.into_sym(ws, state).unwrap().builder(state, ws),
            );
        }
        result
    }
}

impl<'a> DenseTensor<Expr<'a>> {
    pub fn symbolic_zeros(structure: TensorSkeleton) -> DenseTensor<Atom> {
        let result_data = vec![0; structure.size()];

        DenseTensor {
            data: result_data.iter().map(|&x| Atom::new_num(x)).collect(),
            structure,
        }
    }

    pub fn neutral_builder(
        structure: TensorSkeleton,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>> {
        let zero = Atom::new_num(0).builder(state, ws);
        let result_data = vec![zero; structure.size()];
        DenseTensor {
            data: result_data,
            structure,
        }
    }
    pub fn symbolic_labels(
        label: &str,
        structure: TensorSkeleton,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> DenseTensor<Atom> {
        let mut data = vec![];
        for index in structure.index_iter() {
            let indices_str = index
                .into_iter()
                .map(|index| index.to_string())
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
        structure: TensorSkeleton,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Atom> {
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
        structure: TensorSkeleton,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>> {
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

    pub fn finish(self) -> DenseTensor<Atom> {
        DenseTensor {
            data: self.data.into_iter().map(|x| x.into_atom()).collect(),
            structure: self.structure,
        }
    }
}

impl<'a, T> DenseTensor<T>
where
    for<'d> &'d T: SymbolicInto,
{
    pub fn to_symbolic<'c: 'a, 'b>(
        &'b self,
        ws: &'c Workspace,
        state: &'c mut State,
    ) -> DenseTensor<Atom> {
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
    ) -> DenseTensor<Expr<'a>> {
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

impl<T, U, Out> SymbolicContract<NumTensor<T>> for NumTensor<U>
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
{
    type LCM = NumTensor<Out>;
    fn contract_sym(
        &self,
        other: &NumTensor<T>,
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

#[enum_dispatch(HasTensorStructure)]
#[derive(Debug, Clone)]
pub enum MixedTensors {
    Float(NumTensor<f64>),
    Complex(NumTensor<Complex<f64>>),
    Symbolic(NumTensor<Atom>),
}

impl From<DenseTensor<f64>> for MixedTensors {
    fn from(other: DenseTensor<f64>) -> Self {
        MixedTensors::Float(NumTensor::Dense(other))
    }
}

impl From<SparseTensor<f64>> for MixedTensors {
    fn from(other: SparseTensor<f64>) -> Self {
        MixedTensors::Float(NumTensor::Sparse(other))
    }
}

impl From<DenseTensor<Complex<f64>>> for MixedTensors {
    fn from(other: DenseTensor<Complex<f64>>) -> Self {
        MixedTensors::Complex(NumTensor::Dense(other))
    }
}

impl From<SparseTensor<Complex<f64>>> for MixedTensors {
    fn from(other: SparseTensor<Complex<f64>>) -> Self {
        MixedTensors::Complex(NumTensor::Sparse(other))
    }
}

impl From<DenseTensor<Atom>> for MixedTensors {
    fn from(other: DenseTensor<Atom>) -> Self {
        MixedTensors::Symbolic(NumTensor::Dense(other))
    }
}

impl From<SparseTensor<Atom>> for MixedTensors {
    fn from(other: SparseTensor<Atom>) -> Self {
        MixedTensors::Symbolic(NumTensor::Sparse(other))
    }
}

impl SymbolicContract<MixedTensors> for MixedTensors {
    type LCM = MixedTensors;
    fn contract_sym(
        &self,
        other: &MixedTensors,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        match (self, other) {
            (MixedTensors::Float(s), MixedTensors::Float(o)) => {
                Some(MixedTensors::Float(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Float(s), MixedTensors::Complex(o)) => {
                Some(MixedTensors::Complex(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Float(s), MixedTensors::Symbolic(o)) => {
                Some(MixedTensors::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Complex(s), MixedTensors::Float(o)) => {
                Some(MixedTensors::Complex(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Complex(s), MixedTensors::Complex(o)) => {
                Some(MixedTensors::Complex(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Complex(s), MixedTensors::Symbolic(o)) => {
                Some(MixedTensors::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Symbolic(s), MixedTensors::Float(o)) => {
                Some(MixedTensors::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Symbolic(s), MixedTensors::Complex(o)) => {
                Some(MixedTensors::Symbolic(s.contract_sym(o, state, ws)?))
            }
            (MixedTensors::Symbolic(s), MixedTensors::Symbolic(o)) => {
                Some(MixedTensors::Symbolic(s.contract_sym(o, state, ws)?))
            }
        }
    }
}
