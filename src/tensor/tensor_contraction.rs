use super::*;

impl<T> DenseTensor<T>
where
    T: Default
        + Clone
        + std::ops::AddAssign
        + std::ops::SubAssign
        + std::ops::Mul<Output = T>
        + Copy
        + std::cmp::PartialEq
        + std::fmt::Debug,
{
    pub fn internal_contract(&self) -> Self {
        let mut result = self.clone();
        for trace in self.traces() {
            let new_structure = self
                .structure()
                .clone()
                .into_iter()
                .enumerate()
                .filter(|&(i, _)| !trace.contains(&i))
                .map(|(_, x)| x)
                .collect();

            let mut new_result = DenseTensor::default(new_structure);
            for (idx, t) in result.iter_trace(trace) {
                new_result.set(&idx, t);
            }
            result = new_result;
        }
        result
    }

    pub fn contract_with_dense(&self, other: &Self) -> Option<Self> {
        if let Some((i, j)) = self.match_index(other) {
            // println!("{},{}", i, j);
            let self_shape = self.shape();

            let dimension_of_contraction = self_shape[i];
            let metric = self.structure()[i].representation.negative();

            let final_structure = self.structure().merge_at(other.structure(), (i, j));

            // Initialize result tensor with default values
            let mut result_data = vec![T::default(); final_structure.size()];

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

                    for k in 0..dimension_of_contraction {
                        // Adjust indices for fetching from the other tensor
                        if metric[k] {
                            result_data[result_index] -= *fiber_a[k] * *fiber_b[k];
                        } else {
                            result_data[result_index] += *fiber_a[k] * *fiber_b[k];
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

impl<T> SparseTensor<T>
where
    T: Default
        + Clone
        + std::ops::AddAssign
        + std::ops::SubAssign
        + std::ops::Mul<Output = T>
        + Copy
        + std::cmp::PartialEq
        + std::fmt::Debug,
{
    pub fn internal_contract(&self) -> Self {
        let mut result = self.clone();
        for trace in self.traces() {
            // println!("trace {:?}", trace);
            let new_structure = self
                .structure()
                .clone()
                .into_iter()
                .enumerate()
                .filter(|&(i, _)| !trace.contains(&i))
                .map(|(_, x)| x)
                .collect();

            let mut new_result = SparseTensor::empty(new_structure);
            for (idx, t) in result.iter_trace(trace).filter(|(_, t)| *t != T::default()) {
                new_result.set(&idx, t).unwrap();
            }
            result = new_result;
        }
        result
    }
}

pub trait ContractableWithSparse<T> {
    fn contract_with_sparse(&self, other: &SparseTensor<T>) -> Option<Self>
    where
        Self: std::marker::Sized;
}

pub trait ContractableWithDense<T> {
    fn contract_with_dense(&self, other: &DenseTensor<T>) -> Option<DenseTensor<T>>;
}

macro_rules! contract_with_dense_impl {
    ($t:ty,$u:ty) => {
        impl ContractableWithDense<$t> for SparseTensor<$u> {
            fn contract_with_dense(&self, other: &DenseTensor<$t>) -> Option<DenseTensor<$t>> {
                if let Some((i, j)) = self.match_index(other) {
                    let final_structure = self.structure().merge_at(other.structure(), (i, j));
                    let mut result_data = vec![<$t>::default(); final_structure.size()];

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
                                // Adjust indices for fetching from the other tensor
                                if metric[*k] {
                                    result_data[result_index] -= *fiber_a[i] * *fiber_b[*k];
                                } else {
                                    result_data[result_index] += *fiber_a[i] * *fiber_b[*k];
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
        impl ContractableWithDense<$t> for DenseTensor<$u> {
            fn contract_with_dense(&self, other: &DenseTensor<$t>) -> Option<DenseTensor<$t>> {
                if let Some((i, j)) = self.match_index(other) {
                    // println!("{},{}", i, j);
                    let self_shape = self.shape();

                    let dimension_of_contraction = self_shape[i];
                    let metric = self.structure()[i].representation.negative();

                    let final_structure = self.structure().merge_at(other.structure(), (i, j));

                    // Initialize result tensor with default values
                    let mut result_data = vec![<$t>::default(); final_structure.size()];

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

                            for k in 0..dimension_of_contraction {
                                // Adjust indices for fetching from the other tensor
                                if metric[k] {
                                    result_data[result_index] -= *fiber_a[k] * *fiber_b[k];
                                } else {
                                    result_data[result_index] += *fiber_a[k] * *fiber_b[k];
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
    };
}

contract_with_dense_impl!(f64, f64);
contract_with_dense_impl!(num::Complex<f64>, num::Complex<f64>);

contract_with_dense_impl!(num::Complex<f64>, f64);

macro_rules! contract_with_sparse_impl {
    ($t:ty,$u:ty) => {
        impl ContractableWithSparse<$u> for SparseTensor<$t> {
            fn contract_with_sparse(&self, other: &SparseTensor<$u>) -> Option<Self> {
                if let Some((i, j)) = self.match_index(other) {
                    let final_structure = self.structure().merge_at(other.structure(), (i, j));
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

                            let mut value = <$t>::default();
                            let mut nonzero = false;
                            for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                                nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)) // Only store the positions
                            }) {
                                // Adjust indices for fetching from the other tensor
                                if metric[x] {
                                    value -= *fiber_a[i] * *fiber_b[j];
                                } else {
                                    value += *fiber_a[i] * *fiber_b[j];
                                }

                                nonzero = true;
                            }

                            if nonzero && value != <$t>::default() {
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
        impl ContractableWithSparse<$u> for DenseTensor<$t> {
            fn contract_with_sparse(&self, other: &SparseTensor<$u>) -> Option<Self> {
                other.contract_with_dense(self)
            }
        }
    };
}

contract_with_sparse_impl!(f64, f64);
contract_with_sparse_impl!(num::Complex<f64>, num::Complex<f64>);

contract_with_sparse_impl!(num::Complex<f64>, f64);
