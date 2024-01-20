use num::Complex;

use super::*;
use std::ops::Neg;
impl<T> DenseTensor<T>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + Neg<Output = T>
        + Clone
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

            let mut new_result = DenseTensor::from_data_coerced(&self.data, new_structure).unwrap();
            for (idx, t) in result.iter_trace(trace) {
                new_result.set(&idx, t);
            }
            result = new_result;
        }
        result
    }
}

impl<T> SparseTensor<T>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + std::ops::Neg<Output = T>
        + Clone,
{
    pub fn internal_contract(&self) -> Self {
        let trace = self.traces()[0];

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
        for (idx, t) in self.iter_trace(trace) {
            //filter(|(_, t)| *t != T::default())
            new_result.set(&idx, t).unwrap();
        }

        if new_result.traces().is_empty() {
            new_result
        } else {
            new_result.internal_contract()
        }
    }
}

// pub trait NumTensor {}

// impl<T> NumTensor for DenseTensor<T> {}
// impl<T> NumTensor for SparseTensor<T> {}
// pub trait Densest<A: NumTensor, B: NumTensor> {}
// pub trait Contractable<T> {
//     fn contract<C: NumTensor + Densest<SparseTensor<T>, Self>>(
//         &self,
//         other: &SparseTensor<T>,
//     ) -> Option<C>;
// }

// pub fn mul<'a, 'b, T, U>(right: &'a T, left: &'b U) -> <&'b U as SmallestUpgrade<&'b T>>::LCM
// where
//     &'b U: SmallestUpgrade<&'b T>,
//     &'a T: SmallestUpgrade<&'b U, LCM = <&'b U as SmallestUpgrade<&'b T>>::LCM>,
//     <&'b U as SmallestUpgrade<&'b T>>::LCM: std::ops::Mul<
//         <&'b U as SmallestUpgrade<&'b T>>::LCM,
//         Output = <&'b U as SmallestUpgrade<&'b T>>::LCM,
//     >,
// {
//     left.upgrade() * right.upgrade()
// }

pub fn mul<T, U>(right: T, left: U) -> U::LCM
where
    U: SmallestUpgrade<T>,
    T: SmallestUpgrade<U>,
    U::LCM: std::ops::Mul<T::LCM, Output = U::LCM>,
{
    left.upgrade() * right.upgrade()
}
pub trait Contract<T> {
    type LCM;
    fn contract(&self, other: &T) -> Option<Self::LCM>;
}

impl<T, U> Contract<DenseTensor<T>> for DenseTensor<U>
where
    T: SmallestUpgrade<U> + Copy,
    U: SmallestUpgrade<T, LCM = T::LCM> + Copy,
    T::LCM: std::ops::AddAssign<T::LCM>
        + std::ops::SubAssign<T::LCM>
        + for<'a> std::ops::AddAssign<&'a T::LCM>
        + for<'b> std::ops::SubAssign<&'b T::LCM>
        + std::fmt::Debug
        + Neg<Output = T::LCM>
        + Default
        + Clone
        + std::ops::Mul<T::LCM, Output = T::LCM>,
{
    type LCM = DenseTensor<T::LCM>;
    fn contract(&self, other: &DenseTensor<T>) -> Option<Self::LCM> {
        if let Some((i, j)) = self.match_index(other) {
            // println!("{},{}", i, j);
            let self_shape = self.shape();

            let dimension_of_contraction = self_shape[i];
            let metric = self.structure()[i].representation.negative();

            let final_structure = self.structure().merge_at(other.structure(), (i, j));

            // Initialize result tensor with default values
            let mut result_data = vec![T::LCM::default(); final_structure.size()];

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
                            result_data[result_index] -= mul(*fiber_a[k], *fiber_b[k]);
                        } else {
                            result_data[result_index] += mul(*fiber_a[k], *fiber_b[k]);
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

impl<T, U> Contract<DenseTensor<T>> for SparseTensor<U>
where
    T: SmallestUpgrade<U> + Copy,
    U: SmallestUpgrade<T, LCM = T::LCM> + Copy,
    T::LCM: std::ops::AddAssign<T::LCM>
        + std::ops::SubAssign<T::LCM>
        + for<'a> std::ops::AddAssign<&'a T::LCM>
        + for<'b> std::ops::SubAssign<&'b T::LCM>
        + std::fmt::Debug
        + Neg<Output = T::LCM>
        + Default
        + Clone
        + std::ops::Mul<T::LCM, Output = T::LCM>,
{
    type LCM = DenseTensor<T::LCM>;
    fn contract(&self, other: &DenseTensor<T>) -> Option<Self::LCM> {
        // println!("Contracting SparseTensor DenseTensor");
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![<T::LCM>::default(); final_structure.size()];

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
                            result_data[result_index] -= mul(*fiber_a[i], *fiber_b[*k]);
                        } else {
                            result_data[result_index] += mul(*fiber_a[i], *fiber_b[*k]);
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

impl<T, U> Contract<SparseTensor<T>> for SparseTensor<U>
where
    T: SmallestUpgrade<U> + Copy,
    U: SmallestUpgrade<T, LCM = T::LCM> + Copy,
    T::LCM: std::ops::AddAssign<T::LCM>
        + std::ops::SubAssign<T::LCM>
        + for<'a> std::ops::AddAssign<&'a T::LCM>
        + for<'b> std::ops::SubAssign<&'b T::LCM>
        + std::fmt::Debug
        + Neg<Output = T::LCM>
        + std::cmp::PartialOrd
        + Default
        + Clone
        + std::ops::Mul<T::LCM, Output = T::LCM>,
{
    type LCM = SparseTensor<T::LCM>;
    fn contract(&self, other: &SparseTensor<T>) -> Option<Self::LCM> {
        // println!("Contracting SparseTensor SparseTensor");
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

                    let mut value = <T::LCM>::default();
                    let mut nonzero = false;
                    for (i, j, x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                        nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)) // Only store the positions
                    }) {
                        // Adjust indices for fetching from the other tensor
                        if metric[x] {
                            value -= mul(*fiber_a[i], *fiber_b[j]);
                        } else {
                            value += mul(*fiber_a[i], *fiber_b[j]);
                        }

                        nonzero = true;
                    }

                    if nonzero && value != <T::LCM>::default() {
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

impl<T, U> Contract<SparseTensor<T>> for DenseTensor<U>
where
    T: SmallestUpgrade<U> + Copy,
    U: SmallestUpgrade<T, LCM = T::LCM> + Copy,
    T::LCM: std::ops::AddAssign<T::LCM>
        + std::ops::SubAssign<T::LCM>
        + for<'a> std::ops::AddAssign<&'a T::LCM>
        + for<'b> std::ops::SubAssign<&'b T::LCM>
        + std::fmt::Debug
        + Neg<Output = T::LCM>
        + Default
        + Clone
        + std::ops::Mul<T::LCM, Output = T::LCM>,
{
    type LCM = DenseTensor<T::LCM>;
    fn contract(&self, other: &SparseTensor<T>) -> Option<Self::LCM> {
        other.contract(self)
    }
}

enum UpgradableType {
    F64(f64),
    Complex(Complex<f64>),
}

enum ContractableTensor {
    Dense(DenseTensor<UpgradableType>),
    Sparse(SparseTensor<UpgradableType>),
}

struct TensorNetwork {
    tensors: Vec<ContractableTensor>,
}
