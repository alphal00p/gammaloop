use std::sync::Arc;

use ahash::{AHashMap, HashMap, HashMapExt};
use arbitrary_int::Number;
use enum_try_as_inner::EnumTryAsInner;

use num::Complex;
use smartstring::alias::String;
use symbolica::{
    domains::{
        float::NumericalFloatLike,
        rational::{Rational, RationalField},
    },
    poly::{evaluate::InstructionEvaluator, polynomial::MultivariatePolynomial, Variable},
    representations::{AsAtomView, Atom, AtomView, FunctionBuilder, Identifier},
    state::{State, Workspace},
};

use crate::tensor::{NamedStructure, Shadowable};

use super::{
    DataTensor, DenseTensor, HasName, HistoryStructure, SetTensorData, Slot, SparseTensor,
    SymbolicAdd, SymbolicAddAssign, SymbolicInto, SymbolicMul, SymbolicNeg,
    SymbolicStructureContract, SymbolicSub, SymbolicSubAssign, SymbolicZero, TensorStructure,
    TracksCount,
};

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
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self {
        let mut result: DenseTensor<T, I> = self.clone();
        for trace in self.traces() {
            let mut new_structure = self.structure.clone();
            new_structure.trace_sym(trace[0], trace[1], state, ws);

            let mut new_result: DenseTensor<T, I> =
                DenseTensor::from_data_coerced(&self.data, new_structure).unwrap();
            for (idx, t) in result.iter_symbolic_trace(trace, state, ws) {
                let _ = new_result.set(&idx, t);
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
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    fn internal_contract_sym(&self, ws: &Workspace, state: &State) -> Self {
        let trace = self.traces()[0];

        // println!("trace {:?}", trace);
        let mut new_structure = self.structure.clone();
        new_structure.trace_sym(trace[0], trace[1], state, ws);

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
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    type LCM = DenseTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &DenseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure =
                        self.structure
                            .merge_at_sym(&other.structure, (i, j), state, ws);
                    let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);
                    let dimension = self_iter.fiber_dimension;
                    let metric = self.external_structure()[i].representation.negative();

                    for fiber_a in self_iter.by_ref() {
                        for fiber_b in other_iter.by_ref() {
                            for i in 0..dimension {
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
                        other_iter.reset();
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    return Some(result);
                }
                // println!("SparseTensor DenseTensor");
                let (permutation, self_matches, other_matches) =
                    self.structure().match_indices(other.structure()).unwrap();

                let mut final_structure = self.structure.clone();
                final_structure.merge_sym(&other.structure, state, ws);

                // Initialize result tensor with default values
                let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
                let mut result_index = 0;
                let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);

                let mut other_iter = other.iter_multi_fibers(&other_matches);
                while let Some(fiber_a) = selfiter.next() {
                    for fiber_b in other_iter.by_ref() {
                        for k in 0..fiber_a.len() {
                            if fiber_a[k].1 {
                                result_data[result_index] = result_data[result_index].sub_sym(
                                    &fiber_a[k].0.mul_sym(fiber_b[selfiter.map[k]], ws, state)?,
                                    ws,
                                    state,
                                )?;
                            } else {
                                result_data[result_index] = result_data[result_index].add_sym(
                                    &fiber_a[k].0.mul_sym(fiber_b[selfiter.map[k]], ws, state)?,
                                    ws,
                                    state,
                                )?;
                            }
                        }
                        result_index += 1;
                    }
                    other_iter.reset();
                }
                let result: DenseTensor<Out, I> = DenseTensor {
                    data: result_data,
                    structure: final_structure,
                };

                return Some(result);
            }
            return other.contract_sym(self, state, ws);
        }
        let mut final_structure = self.structure.clone();
        final_structure.merge_sym(&other.structure, state, ws);

        let mut result_data = vec![Out::zero(state, ws); final_structure.size()];

        let stride = other.size();
        for (i, u) in self.iter_flat() {
            for (j, t) in other.iter_flat() {
                result_data[i * stride + j] = u.mul_sym(t, ws, state).unwrap();
            }
        }
        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };
        Some(result)
    }
}

impl<T, U, I, Out> SymbolicContract<DenseTensor<T, I>> for SparseTensor<U, I>
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
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    type LCM = DenseTensor<Out, I>;
    #[allow(clippy::too_many_lines)]
    fn contract_sym(
        &self,
        other: &DenseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure =
                        self.structure
                            .merge_at_sym(&other.structure, (i, j), state, ws);
                    let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
                    let metric = self.get_ith_metric(i).unwrap();
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    for (skipped, nonzeros, fiber_a) in self_iter.by_ref() {
                        result_index += skipped;
                        for fiber_b in other_iter.by_ref() {
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
                        other_iter.reset();
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                } else {
                    // println!("SparseTensor DenseTensor");
                    let (permutation, self_matches, other_matches) =
                        self.structure().match_indices(other.structure()).unwrap();

                    let mut final_structure = self.structure.clone();
                    final_structure.merge_sym(&other.structure, state, ws);
                    let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
                    let mut result_index = 0;

                    let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);
                    let mut other_iter = other.iter_multi_fibers(&other_matches);

                    while let Some((skipped, nonzeros, fiber_a)) = selfiter.next() {
                        result_index += skipped;
                        for fiber_b in other_iter.by_ref() {
                            for (i, k) in nonzeros.iter().enumerate() {
                                if fiber_a[i].1 {
                                    result_data[result_index] = result_data[result_index].sub_sym(
                                        &fiber_a[i].0.mul_sym(
                                            fiber_b[selfiter.map[*k]],
                                            ws,
                                            state,
                                        )?,
                                        ws,
                                        state,
                                    )?;
                                } else {
                                    result_data[result_index] = result_data[result_index].add_sym(
                                        &fiber_a[i].0.mul_sym(
                                            fiber_b[selfiter.map[*k]],
                                            ws,
                                            state,
                                        )?,
                                        ws,
                                        state,
                                    )?;
                                }
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }
                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                }
            } else {
                other.contract_sym(self, state, ws)
            }
        } else {
            let mut final_structure = self.structure.clone();
            final_structure.merge_sym(&other.structure, state, ws);

            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let stride = other.size();

            for (i, u) in self.iter_flat() {
                for (j, t) in other.iter_flat() {
                    result_data[i * stride + j] = u.mul_sym(t, ws, state).unwrap();
                }
            }
            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };
            Some(result)
        }
    }
}

#[allow(clippy::if_same_then_else)]
impl<T, U, Out, I> SymbolicContract<SparseTensor<T, I>> for SparseTensor<U, I>
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
    I: TensorStructure + Clone + SymbolicStructureContract,
    T: Clone,
    U: Clone,
{
    type LCM = SparseTensor<Out, I>;
    #[allow(clippy::too_many_lines)]
    fn contract_sym(
        &self,
        other: &SparseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure =
                        self.structure
                            .merge_at_sym(&other.structure, (i, j), state, ws);
                    let mut result_data = AHashMap::default();
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.external_structure()[i].representation.negative();

                    let stride = *final_structure
                        .strides()
                        .get(i.saturating_sub(1))
                        .unwrap_or(&1);

                    for (skipped_a, nonzeros_a, fiber_a) in self_iter.by_ref() {
                        result_index += skipped_a;
                        for (skipped_b, nonzeros_b, fiber_b) in other_iter.by_ref() {
                            result_index += skipped_b * stride;
                            let mut value = Out::zero(state, ws);
                            let mut nonzero = false;

                            for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                                nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x))
                            }) {
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
                                nonzero = true;
                            }

                            if nonzero {
                                // && value.as_atom_view() != zero.as_atom_view() TODO impl generic
                                result_data.insert(result_index, value);
                            }
                            result_index += 1;
                        }
                        other_iter.reset();
                    }

                    let result = SparseTensor {
                        elements: result_data,
                        structure: final_structure,
                    };

                    return Some(result);
                }
                let (permutation, self_matches, other_matches) =
                    self.structure().match_indices(other.structure()).unwrap();

                let mut final_structure = self.structure.clone();
                final_structure.merge_sym(&other.structure, state, ws);
                let mut result_data = AHashMap::default();
                let one = if let Some(o) = final_structure.strides().first() {
                    *o
                } else {
                    1
                };
                let stride = *final_structure
                    .strides()
                    .get(i.checked_sub(1)?)
                    .unwrap_or(&one);

                let mut result_index = 0;

                let mut self_iter = self.iter_multi_fibers_metric(&self_matches, permutation);
                let mut other_iter = other.iter_multi_fibers(&other_matches);
                while let Some((skipped_a, nonzeros_a, fiber_a)) = self_iter.next() {
                    result_index += skipped_a;
                    for (skipped_b, nonzeros_b, fiber_b) in other_iter.by_ref() {
                        result_index += skipped_b * stride;
                        let mut value = Out::zero(state, ws);
                        let mut nonzero = false;

                        for (i, j) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                            nonzeros_b
                                .binary_search(&self_iter.map[x])
                                .ok()
                                .map(|j| (i, j))
                        }) {
                            if fiber_a[i].1 {
                                value = value.sub_sym(
                                    &fiber_a[i].0.mul_sym(fiber_b[j], ws, state)?,
                                    ws,
                                    state,
                                )?;
                            } else {
                                value = value.add_sym(
                                    &fiber_a[i].0.mul_sym(fiber_b[j], ws, state)?,
                                    ws,
                                    state,
                                )?;

                                nonzero = true;
                            }
                        }

                        if nonzero {
                            result_data.insert(result_index, value);
                        }
                        result_index += 1;
                    }
                    other_iter.reset();
                }

                let result = SparseTensor {
                    elements: result_data,
                    structure: final_structure,
                };

                return Some(result);
            }
            return other.contract_sym(self, state, ws);
        }
        let mut final_structure = self.structure.clone();
        final_structure.merge_sym(&other.structure, state, ws);

        let mut result_data = AHashMap::default();
        let stride = other.size();
        for (i, u) in self.iter_flat() {
            for (j, t) in other.iter_flat() {
                result_data.insert(i * stride + j, u.mul_sym(t, ws, state).unwrap());
            }
        }
        let result = SparseTensor {
            elements: result_data,
            structure: final_structure,
        };
        Some(result)
    }
}

impl<T, U, I, Out> SymbolicContract<SparseTensor<T, I>> for DenseTensor<U, I>
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
    I: TensorStructure + Clone + SymbolicStructureContract,
{
    type LCM = DenseTensor<Out, I>;
    #[allow(clippy::too_many_lines)]
    fn contract_sym(
        &self,
        other: &SparseTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    let final_structure =
                        self.structure
                            .merge_at_sym(&other.structure, (i, j), state, ws);
                    let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
                    let mut result_index = 0;

                    let mut self_iter = self.iter_fiber(i);
                    let mut other_iter = other.iter_fiber(j);

                    let metric = self.external_structure()[i].representation.negative();

                    let stride = *final_structure
                        .strides()
                        .get(i.saturating_sub(1))
                        .unwrap_or(&1);

                    for fiber_a in self_iter.by_ref() {
                        for (skipped, nonzeros, fiber_b) in other_iter.by_ref() {
                            result_index += skipped * stride;
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
                        other_iter.reset();
                    }

                    let result = DenseTensor {
                        data: result_data,
                        structure: final_structure,
                    };

                    Some(result)
                } else {
                    let (permutation, self_matches, other_matches) =
                        self.structure().match_indices(other.structure()).unwrap();

                    let mut final_structure = self.structure.clone();
                    final_structure.merge_sym(&other.structure, state, ws);

                    let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
                    let mut result_index = 0;

                    let one = if let Some(o) = final_structure.strides().first() {
                        *o
                    } else {
                        1
                    };
                    let stride = *final_structure
                        .strides()
                        .get(i.saturating_sub(1))
                        .unwrap_or(&one);

                    let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);

                    while let Some(fiber_a) = selfiter.next() {
                        for (skipped, nonzeros, fiber_b) in other.iter_multi_fibers(&other_matches)
                        {
                            result_index += skipped * stride;
                            for (i, k) in nonzeros.iter().enumerate() {
                                if fiber_a[i].1 {
                                    result_data[result_index] = result_data[result_index].sub_sym(
                                        &fiber_a[i].0.mul_sym(
                                            fiber_b[selfiter.map[*k]],
                                            ws,
                                            state,
                                        )?,
                                        ws,
                                        state,
                                    )?;
                                } else {
                                    result_data[result_index] = result_data[result_index].add_sym(
                                        &fiber_a[i].0.mul_sym(
                                            fiber_b[selfiter.map[*k]],
                                            ws,
                                            state,
                                        )?,
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

                    Some(result)
                }
            } else {
                other.contract_sym(self, state, ws)
            }
        } else {
            let mut final_structure = self.structure.clone();
            final_structure.merge_sym(&other.structure, state, ws);

            let mut result_data = vec![Out::zero(state, ws); final_structure.size()];
            let stride = other.size();

            for (i, u) in self.iter_flat() {
                for (j, t) in other.iter_flat() {
                    result_data[i * stride + j] = u.mul_sym(t, ws, state).unwrap();
                }
            }
            let result = DenseTensor {
                data: result_data,
                structure: final_structure,
            };
            Some(result)
        }
    }
}

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
            let _ = result.set(&index, value.into_sym(ws, state).unwrap());
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
        structure: I,
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

            let value = Atom::new_from_view(&value_builder.finish().as_atom_view());

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
            let _ = result.set(&index, value.into_sym(ws, state).unwrap());
        }
        result
    }
}

impl<T, U, I, Out> SymbolicContract<DataTensor<T, I>> for DataTensor<U, I>
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
    I: TensorStructure + Clone + SymbolicStructureContract,
    T: Clone,
    U: Clone,
{
    type LCM = DataTensor<Out, I>;
    fn contract_sym(
        &self,
        other: &DataTensor<T, I>,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        match (self, other) {
            (DataTensor::Dense(s), DataTensor::Dense(o)) => {
                Some(DataTensor::Dense(s.contract_sym(o, state, ws)?))
            }
            (DataTensor::Dense(s), DataTensor::Sparse(o)) => {
                Some(DataTensor::Dense(o.contract_sym(s, state, ws)?))
            }
            (DataTensor::Sparse(s), DataTensor::Dense(o)) => {
                Some(DataTensor::Dense(s.contract_sym(o, state, ws)?))
            }
            (DataTensor::Sparse(s), DataTensor::Sparse(o)) => {
                Some(DataTensor::Sparse(s.contract_sym(o, state, ws)?))
            }
        }
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

impl<I> DenseTensor<Atom, I>
where
    I: Clone,
{
    pub fn to_evaluator<'a, N>(
        &'a self,
        var_map: &mut HashMap<AtomView<'a>, Variable>,
        state: &State,
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
                    .to_polynomial_with_map(&RationalField::new(), var_map);

                println!("{}", poly.printer(state));

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
    let mut state = State::new();
    let ws = Workspace::new();
    let structure = NamedStructure::from_integers(&[(4, 5), (5, 4)], "r");
    let p = structure.shadow(&mut state, &ws).unwrap();

    let mut var_map = HashMap::new();
    let a: DenseTensor<InstructionEvaluator<f64>, NamedStructure> =
        p.to_evaluator(&mut var_map, &state);

    for (k, v) in var_map {
        println!("{} {:?}", k.printer(&state), v);
    }

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

    fn name(&self) -> Option<&Self::Name> {
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

pub type MixedTensors = MixedTensor<HistoryStructure<Identifier>>;

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

impl<I> SymbolicContract<MixedTensor<I>> for MixedTensor<I>
where
    I: TensorStructure + Clone + SymbolicStructureContract,
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
