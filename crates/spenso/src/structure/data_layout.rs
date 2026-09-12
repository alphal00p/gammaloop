use thiserror::Error;

use super::{CanonicalLayout, Canonicalized, TensorStructure};

/// Errors raised while translating between logical row-major data and tensor storage.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum TensorDataLayoutError {
    #[error("tensor data requires concrete representation dimensions")]
    NonConcreteDimension,
    #[error("tensor dimensions overflow addressable storage")]
    Overflow,
    #[error("expected {expected} tensor indices, got {actual}")]
    RankMismatch { expected: usize, actual: usize },
    #[error("index {index} out of bounds for axis {axis} of size {dimension}")]
    IndexOutOfBounds {
        axis: usize,
        index: usize,
        dimension: usize,
    },
    #[error("flat index {index} out of bounds for tensor storage of size {size}")]
    FlatIndexOutOfBounds { index: usize, size: usize },
    #[error("tensor data has {actual} elements, expected {expected}")]
    DataLength { expected: usize, actual: usize },
}

/// Bidirectional row-major coordinate mapping between logical and canonical storage axes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TensorDataLayout {
    logical_shape: Vec<usize>,
    logical_strides: Vec<usize>,
    storage_axes: Vec<usize>,
    storage_shape: Vec<usize>,
    storage_strides: Vec<usize>,
    size: usize,
}

impl TensorDataLayout {
    pub fn from_canonicalized<S: TensorStructure>(
        canonicalized: &Canonicalized<S>,
    ) -> Result<Self, TensorDataLayoutError> {
        let storage_shape = canonicalized
            .canonical()
            .external_dims_iter()
            .map(|dimension| {
                usize::try_from(dimension).map_err(|_| TensorDataLayoutError::NonConcreteDimension)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_storage_shape(&storage_shape, canonicalized.layout())
    }

    pub fn from_storage_shape(
        storage_shape: &[usize],
        layout: &CanonicalLayout,
    ) -> Result<Self, TensorDataLayoutError> {
        if storage_shape.len() != layout.order() {
            return Err(TensorDataLayoutError::RankMismatch {
                expected: layout.order(),
                actual: storage_shape.len(),
            });
        }

        let logical_shape = layout.canonical_to_logical(storage_shape);
        let logical_axes = (0..logical_shape.len()).collect::<Vec<_>>();
        let storage_axes = layout.logical_to_canonical(&logical_axes);
        let logical_strides = Self::strides(&logical_shape)?;
        let storage_strides = Self::strides(storage_shape)?;
        let size = logical_shape.iter().try_fold(1usize, |size, &dimension| {
            size.checked_mul(dimension)
                .ok_or(TensorDataLayoutError::Overflow)
        })?;

        Ok(Self {
            logical_shape,
            logical_strides,
            storage_axes,
            storage_shape: storage_shape.to_vec(),
            storage_strides,
            size,
        })
    }

    fn strides(shape: &[usize]) -> Result<Vec<usize>, TensorDataLayoutError> {
        let mut strides = vec![1usize; shape.len()];
        for axis in (0..shape.len().saturating_sub(1)).rev() {
            strides[axis] = strides[axis + 1]
                .checked_mul(shape[axis + 1])
                .ok_or(TensorDataLayoutError::Overflow)?;
        }
        Ok(strides)
    }

    pub fn logical_shape(&self) -> &[usize] {
        &self.logical_shape
    }

    pub fn storage_shape(&self) -> &[usize] {
        &self.storage_shape
    }

    pub fn size(&self) -> usize {
        self.size
    }

    fn expanded_from_flat(
        &self,
        flat: usize,
        strides: &[usize],
    ) -> Result<Vec<usize>, TensorDataLayoutError> {
        if flat >= self.size {
            return Err(TensorDataLayoutError::FlatIndexOutOfBounds {
                index: flat,
                size: self.size,
            });
        }
        let mut remainder = flat;
        Ok(strides
            .iter()
            .map(|&stride| {
                let index = remainder / stride;
                remainder %= stride;
                index
            })
            .collect())
    }

    fn validate_logical_expanded(&self, logical: &[usize]) -> Result<(), TensorDataLayoutError> {
        if logical.len() != self.logical_shape.len() {
            return Err(TensorDataLayoutError::RankMismatch {
                expected: self.logical_shape.len(),
                actual: logical.len(),
            });
        }
        for (axis, (&index, &dimension)) in logical.iter().zip(&self.logical_shape).enumerate() {
            if index >= dimension {
                return Err(TensorDataLayoutError::IndexOutOfBounds {
                    axis,
                    index,
                    dimension,
                });
            }
        }
        Ok(())
    }

    fn validate_storage_expanded(&self, storage: &[usize]) -> Result<(), TensorDataLayoutError> {
        if storage.len() != self.storage_shape.len() {
            return Err(TensorDataLayoutError::RankMismatch {
                expected: self.storage_shape.len(),
                actual: storage.len(),
            });
        }
        for (axis, (&index, &dimension)) in storage.iter().zip(&self.storage_shape).enumerate() {
            if index >= dimension {
                return Err(TensorDataLayoutError::IndexOutOfBounds {
                    axis,
                    index,
                    dimension,
                });
            }
        }
        Ok(())
    }

    pub fn logical_expanded_to_storage_expanded(
        &self,
        logical: &[usize],
    ) -> Result<Vec<usize>, TensorDataLayoutError> {
        self.validate_logical_expanded(logical)?;
        Ok(self
            .storage_axes
            .iter()
            .map(|&logical_axis| logical[logical_axis])
            .collect())
    }

    pub fn storage_expanded_to_logical_expanded(
        &self,
        storage: &[usize],
    ) -> Result<Vec<usize>, TensorDataLayoutError> {
        self.validate_storage_expanded(storage)?;
        let mut logical = vec![0; storage.len()];
        for (storage_axis, &logical_axis) in self.storage_axes.iter().enumerate() {
            logical[logical_axis] = storage[storage_axis];
        }
        Ok(logical)
    }

    pub fn logical_expanded_to_storage_flat(
        &self,
        logical: &[usize],
    ) -> Result<usize, TensorDataLayoutError> {
        Ok(self
            .logical_expanded_to_storage_expanded(logical)?
            .iter()
            .zip(&self.storage_strides)
            .map(|(&index, &stride)| index * stride)
            .sum())
    }

    pub fn logical_flat_to_storage_flat(
        &self,
        logical_flat: usize,
    ) -> Result<usize, TensorDataLayoutError> {
        let logical = self.expanded_from_flat(logical_flat, &self.logical_strides)?;
        self.logical_expanded_to_storage_flat(&logical)
    }

    pub fn storage_flat_to_logical_expanded(
        &self,
        storage_flat: usize,
    ) -> Result<Vec<usize>, TensorDataLayoutError> {
        let storage = self.expanded_from_flat(storage_flat, &self.storage_strides)?;
        self.storage_expanded_to_logical_expanded(&storage)
    }

    pub fn storage_flat_to_logical_flat(
        &self,
        storage_flat: usize,
    ) -> Result<usize, TensorDataLayoutError> {
        let logical = self.storage_flat_to_logical_expanded(storage_flat)?;
        Ok(logical
            .iter()
            .zip(&self.logical_strides)
            .map(|(&index, &stride)| index * stride)
            .sum())
    }

    pub fn reorder_to_storage<T>(&self, logical: Vec<T>) -> Result<Vec<T>, TensorDataLayoutError> {
        if logical.len() != self.size {
            return Err(TensorDataLayoutError::DataLength {
                expected: self.size,
                actual: logical.len(),
            });
        }
        let mut storage = std::iter::repeat_with(|| None)
            .take(self.size)
            .collect::<Vec<_>>();
        for (logical_flat, value) in logical.into_iter().enumerate() {
            storage[self.logical_flat_to_storage_flat(logical_flat)?] = Some(value);
        }
        Ok(storage
            .into_iter()
            .map(|value| value.expect("axis permutations are bijective"))
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::TensorDataLayout;
    use crate::structure::{
        Canonicalized, IndexLess, TensorStructure,
        abstract_index::AbstractIndex,
        concrete_index::FlatIndex,
        dimension::Dimension,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        representation::{ExtendibleReps, LibraryRep, RepName},
    };
    use crate::tensors::data::{DenseTensor, SparseTensor};

    #[test]
    fn mixed_representation_data_stays_in_logical_row_major_order() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let interface = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::Explicit(AbstractIndex::Normal(0))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(1))),
        ]);
        let layout = TensorDataLayout::from_canonicalized(&interface).unwrap();

        assert_eq!(layout.logical_shape(), [2, 3]);
        assert_eq!(layout.storage_axes, vec![1, 0]);
        assert_eq!(
            layout.reorder_to_storage(vec![0, 1, 2, 3, 4, 5]).unwrap(),
            vec![0, 3, 1, 4, 2, 5]
        );
        assert_eq!(
            (0..6)
                .map(|index| layout.logical_flat_to_storage_flat(index).unwrap())
                .collect::<Vec<_>>(),
            vec![0, 2, 4, 1, 3, 5]
        );
    }

    #[test]
    fn explicit_index_sorting_is_part_of_the_layout_mapping() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(8))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(3))),
        ]);
        let layout = TensorDataLayout::from_canonicalized(&interface).unwrap();

        assert_eq!(layout.storage_axes, vec![1, 0]);
        assert_eq!(
            layout.reorder_to_storage(vec![0, 1, 2, 3]).unwrap(),
            vec![0, 2, 1, 3]
        );
    }

    #[test]
    fn mixed_open_and_explicit_ports_keep_logical_data_order() {
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(2));
        let interface = PartialStructure::from_logical_slots([
            euc.slot(PartialIndex::open(0)),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(3))),
        ]);
        let layout = TensorDataLayout::from_canonicalized(&interface).unwrap();

        assert_eq!(layout.storage_axes, vec![1, 0]);
        assert_eq!(
            layout.reorder_to_storage(vec![0, 1, 2, 3]).unwrap(),
            vec![0, 2, 1, 3]
        );
    }

    #[test]
    fn rank_zero_layout_contains_one_scalar_element() {
        let interface = PartialStructure::from_logical_slots([]);
        let layout = TensorDataLayout::from_canonicalized(&interface).unwrap();

        assert!(layout.logical_shape().is_empty());
        assert!(layout.storage_axes.is_empty());
        assert_eq!(layout.size(), 1);
        assert_eq!(layout.logical_flat_to_storage_flat(0).unwrap(), 0);
        assert_eq!(layout.reorder_to_storage(vec![7]).unwrap(), vec![7]);
    }

    #[test]
    fn both_canonicalization_stages_and_inverse_mapping_are_consistent() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let interface = PartialStructure::from_logical_slots([
            mink.slot(PartialIndex::Explicit(AbstractIndex::Normal(0))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(9))),
            euc.slot(PartialIndex::Explicit(AbstractIndex::Normal(3))),
        ]);
        let layout = TensorDataLayout::from_canonicalized(&interface).unwrap();
        let logical = (0..18).collect::<Vec<_>>();
        let expected_storage = vec![0, 9, 3, 12, 6, 15, 1, 10, 4, 13, 7, 16, 2, 11, 5, 14, 8, 17];

        assert_eq!(layout.logical_shape(), [2, 3, 3]);
        assert_eq!(layout.storage_axes, vec![2, 1, 0]);
        assert_eq!(
            layout.reorder_to_storage(logical.clone()).unwrap(),
            expected_storage
        );
        for logical_flat in 0..layout.size() {
            let storage_flat = layout.logical_flat_to_storage_flat(logical_flat).unwrap();
            assert_eq!(
                layout.storage_flat_to_logical_flat(storage_flat).unwrap(),
                logical_flat
            );
        }
    }

    #[test]
    fn representation_input_order_and_late_index_order_are_separate_stages() {
        let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(2));
        let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
        let input = vec![mink, euc, euc]
            .into_iter()
            .collect::<Canonicalized<IndexLess<LibraryRep, AbstractIndex>>>();
        let layout = TensorDataLayout::from_canonicalized(&input).unwrap();
        let representation_storage = layout.reorder_to_storage((0..18).collect()).unwrap();
        assert_eq!(
            representation_storage,
            vec![0, 9, 1, 10, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17]
        );

        let storage_indices = input.layout().logical_to_canonical(&[
            AbstractIndex::Normal(50),
            AbstractIndex::Normal(9),
            AbstractIndex::Normal(3),
        ]);
        let indexless = input.into_canonical();
        let dense = DenseTensor::from_storage_data(representation_storage, indexless.clone())
            .unwrap()
            .reindex_storage(&storage_indices)
            .unwrap()
            .apply();
        assert_eq!(
            dense.data,
            vec![0, 9, 3, 12, 6, 15, 1, 10, 4, 13, 7, 16, 2, 11, 5, 14, 8, 17]
        );

        let sparse = SparseTensor {
            elements: HashMap::from([(FlatIndex::from(5), 7), (FlatIndex::from(14), 11)]),
            zero: 0,
            structure: indexless,
        }
        .reindex_storage(&storage_indices)
        .unwrap()
        .apply();
        assert_eq!(
            sparse.elements,
            HashMap::from([(FlatIndex::from(13), 7), (FlatIndex::from(10), 11)])
        );
        assert!(!sparse.elements.contains_key(&FlatIndex::from(5)));
        assert!(!sparse.elements.contains_key(&FlatIndex::from(14)));
    }
}
