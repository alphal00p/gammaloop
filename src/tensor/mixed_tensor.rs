use std::collections::BTreeMap;
use symbolica::{
    representations::{number::Number, AsAtomView, Atom},
    state::{State, Workspace},
};

use super::{
    ContractableWithDense, ContractableWithSparse, DenseTensor, Expr, HasTensorStructure,
    SparseTensor, TensorStructure, VecSlotExtension,
};

impl<'a> ContractableWithDense<Expr<'a>> for SparseTensor<Expr<'a>> {
    fn contract_with_dense(&self, other: &DenseTensor<Expr<'a>>) -> Option<DenseTensor<Expr<'a>>> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let state = self.iter().next().unwrap().1.state;
            let workspace = self.iter().next().unwrap().1.workspace;
            let zero = workspace.new_num(0);

            let neutral_summand = zero.builder(state, workspace);
            let mut result_data = vec![neutral_summand; final_structure.size()];
            let metric = self.structure()[i].representation.negative();

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
                            result_data[result_index] =
                                -(fiber_a[i].clone() * fiber_b[*k]) + &result_data[result_index];
                        } else {
                            result_data[result_index] =
                                (fiber_a[i].clone() * fiber_b[*k]) + &result_data[result_index];
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
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<'a> ContractableWithDense<Expr<'a>> for DenseTensor<Expr<'a>> {
    fn contract_with_dense(&self, other: &DenseTensor<Expr<'a>>) -> Option<DenseTensor<Expr<'a>>> {
        if let Some((i, j)) = self.match_index(other) {
            let dimension_of_contraction = self.shape()[i];
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let state = self.iter().next().unwrap().1.state;
            let workspace = self.iter().next().unwrap().1.workspace;
            let zero = workspace.new_num(0);
            let neutral_summand = zero.builder(state, workspace);
            let mut result_data = vec![neutral_summand; final_structure.size()];
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
                            result_data[result_index] -= &(fiber_a[i].clone() * fiber_b[i]);
                        } else {
                            result_data[result_index] += &(fiber_a[i].clone() * fiber_b[i]);
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
                return Some(result.internal_contract());
            }
        }
        None
    }
}

impl<'a> ContractableWithSparse<Expr<'a>> for DenseTensor<Expr<'a>> {
    fn contract_with_sparse(
        &self,
        other: &SparseTensor<Expr<'a>>,
    ) -> Option<DenseTensor<Expr<'a>>> {
        other.contract_with_dense(self)
    }
}

impl<'a> ContractableWithSparse<Expr<'a>> for SparseTensor<Expr<'a>> {
    fn contract_with_sparse(
        &self,
        other: &SparseTensor<Expr<'a>>,
    ) -> Option<SparseTensor<Expr<'a>>> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let state = self.iter().next().unwrap().1.state;
            let workspace = self.iter().next().unwrap().1.workspace;
            let zero = workspace.new_num(0);
            let neutral_summand = zero.builder(state, workspace);
            let mut result_data = BTreeMap::new();
            let metric = self.structure()[i].representation.negative();

            for (index_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {
                    let result_index = index_a[..i]
                        .iter()
                        .chain(&index_a[i + 1..])
                        .chain(&index_b[..j])
                        .chain(&index_b[j + 1..])
                        .cloned()
                        .collect::<Vec<_>>();
                    let mut value = neutral_summand.clone();
                    let mut nonzero = false;
                    for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                        nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)) // Only store the positions
                    }) {
                        // Adjust indices for fetching from the other tensor
                        if metric[x] {
                            value -= &(fiber_a[i].clone() * fiber_b[j]);
                        } else {
                            value += &(fiber_a[i].clone() * fiber_b[j]);
                        }

                        nonzero = true;
                    }

                    if nonzero && value.as_atom_view() != zero.as_atom_view() {
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
                return Some(result.internal_contract());
            }
        }
        None
    }
}

pub trait ConvertableToSymbolic {
    fn to_symbolic<'a>(&self, ws: &'a Workspace, state: &'a State) -> Option<Expr<'a>>;
}

impl ConvertableToSymbolic for f64 {
    fn to_symbolic<'a>(&self, ws: &'a Workspace, state: &'a State) -> Option<Expr<'a>> {
        let rugrat = rug::Rational::from_f64(*self)?;
        let natrat = symbolica::rings::rational::Rational::from_large(rugrat);
        let symrat = ws.new_num(Number::from(natrat));

        return Some(symrat.builder(state, ws));
    }
}

impl ConvertableToSymbolic for Complex<f64> {
    fn to_symbolic<'a>(&self, ws: &'a Workspace, state: &'a State) -> Option<Expr<'a>> {
        let real = self.re.to_symbolic(ws, state)?;
        let imag = self.im.to_symbolic(ws, state)?;
        let zero = ws.new_num(0);
        let one = ws.new_num(1);
        let i = ws.new_num(Number::from(
            symbolica::rings::rational::Rational::from_large(rug::Rational::from_f64(1.0).unwrap()),
        ));
        let symrat = zero + &real + &(i * &imag);

        return Some(symrat.builder(state, ws));
    }
}

impl<T: ConvertableToSymbolic> SparseTensor<T> {
    pub fn to_symbolic<'a>(&self, ws: &'a Workspace, state: &'a State) -> SparseTensor<Expr<'a>> {
        let mut result = SparseTensor::empty(self.structure.clone());
        for (index, value) in self.iter() {
            result
                .set(index, value.to_symbolic(ws, state).unwrap())
                .unwrap();
        }
        result
    }
}

impl<'a> DenseTensor<Expr<'a>> {
    pub fn symbolic_zeros(
        structure: TensorStructure,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'a>> {
        let zero = ws.new_num(0);
        let neutral_summand = zero.builder(state, ws);
        let result_data = vec![neutral_summand; structure.size()];
        DenseTensor {
            data: result_data,
            structure,
        }
    }

    pub fn symbolic_labels(
        label: &str,
        structure: &TensorStructure,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> Vec<Atom> {
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
        data
    }
}

impl<T: ConvertableToSymbolic> DenseTensor<T> {
    pub fn to_symbolic<'a: 'b, 'b>(
        &'b self,
        ws: &'a Workspace,
        state: &'a State,
    ) -> DenseTensor<Expr<'b>> {
        let mut result = DenseTensor::symbolic_zeros(self.structure.clone(), ws, state);
        for (index, value) in self.iter() {
            result.set(&index, value.to_symbolic(ws, state).unwrap());
        }
        result
    }
}
