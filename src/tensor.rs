use std::{
    borrow::Cow,
    collections::{ BTreeMap, HashMap},
    iter::FromIterator,
    ops::Bound::Included,
    usize,
};


type AbstractIndex = usize;
type Dimension = usize;
type ConcreteIndex = usize;
type Position = usize;

#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash)]
pub enum Signature {
    Euclidean(Dimension),
    Lorentz(Dimension),
}

impl Signature {
    pub fn negative(&self) -> Vec<bool> {
        match self {
            Signature::Euclidean(value) => vec![false; *value],
            Signature::Lorentz(value) => std::iter::once(false)
                .chain(std::iter::repeat(true).take(*value - 1))
                .collect::<Vec<_>>(),
        }
    }
}

impl From<Dimension> for Signature {
    fn from(value: Dimension) -> Self {
        Signature::Euclidean(value)
    }
}

impl<'a> std::iter::FromIterator<&'a Signature> for Vec<Dimension> {
    fn from_iter<T: IntoIterator<Item = &'a Signature>>(iter: T) -> Self {
        iter.into_iter()
            .map(|&rep| -> Dimension { (&rep).into() })
            .collect()
    }
}

impl From<&Signature> for Dimension {
    fn from(rep: &Signature) -> Self {
        match rep {
            Signature::Euclidean(value) => *value,
            Signature::Lorentz(value) => *value,
        }
    }
}

impl From<Signature> for Dimension {
    fn from(rep: Signature) -> Self {
        match rep {
            Signature::Euclidean(value) => value,
            Signature::Lorentz(value) => value,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone, Hash)]
pub struct Slot {
    index: AbstractIndex,
    signature: Signature,
}

impl From<(AbstractIndex, Signature)> for Slot {
    fn from(value: (AbstractIndex, Signature)) -> Self {
        Slot {
            index: value.0,
            signature: value.1,
        }
    }
}

pub type TensorStructure = Vec<Slot>;

pub trait VecSlotExt {
    fn from_idxsing(indices: &[AbstractIndex], dims: &[Signature]) -> Self;
    fn from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self;
    fn match_index(&self, other: &Self) -> Option<(Position, Position)>;
    fn traces(&self) -> Vec<Vec<usize>>;
    fn merge_at(&self, other: &Self, positions: (Position, Position)) -> Self;
    fn shape(&self) -> Vec<Dimension>;
    fn order(&self) -> usize;
    fn strides_row_major(&self) -> Vec<usize>;
    fn strides_column_major(&self) -> Vec<usize>;
    fn strides(&self) -> Vec<usize>;
    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String>;
    fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String>;
    fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String>;
    fn size(&self) -> usize;
}

impl VecSlotExt for TensorStructure {
    fn from_idxsing(indices: &[AbstractIndex], signatures: &[Signature]) -> Self {
        indices
            .iter()
            .zip(signatures.iter())
            .map(|(&index, &dim)| Slot::from((index, dim)))
            .collect()
    }

    fn from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self {
        indices
            .iter()
            .zip(dims.iter())
            .map(|(&index, &dim)| Slot::from((index, Signature::Euclidean(dim))))
            .collect()
    }

    fn match_index(&self, other: &Self) -> Option<(Position, Position)> {
        for (i, slot_a) in self.iter().enumerate().rev() {
            for (j, slot_b) in other.iter().enumerate() {
                if slot_a == slot_b {
                    return Some((i, j));
                }
            }
        }
        None
    }

    fn traces(&self) -> Vec<Vec<usize>> {
        let mut positions = HashMap::new();

        // Track the positions of each element
        for (index, &value) in self.iter().enumerate() {
            positions.entry(value).or_insert_with(Vec::new).push(index);
        }

        // Collect only the positions of repeated elements
        positions
            .into_iter()
            .filter_map(|(_, indices)| {
                if indices.len() > 1 {
                    Some(indices)
                } else {
                    None
                }
            })
            .collect()
    }

    fn merge_at(&self, other: &Self, positions: (Position, Position)) -> Self {
        let mut slots_b = other.clone();
        let mut slots_a = self.clone();

        slots_a.remove(positions.0);
        slots_b.remove(positions.1);

        slots_a.append(&mut slots_b);
        slots_a
    }

    fn shape(&self) -> Vec<Dimension> {
        self.iter().map(|slot| &slot.signature).collect()
    }

    fn order(&self) -> usize {
        //total valence (or misnamed : rank)
        self.len()
    }

    fn strides_column_major(&self) -> Vec<usize> {
        let mut strides: Vec<usize> = vec![1; self.order()];

        if self.order() == 0 {
            return strides;
        }

        for i in 0..self.order() - 1 {
            strides[i + 1] = strides[i] * usize::from(self[i].signature);
        }

        strides
    }

    fn strides_row_major(&self) -> Vec<usize> {
        let mut strides = vec![1; self.order()];
        if self.order() == 0 {
            return strides;
        }

        for i in (0..self.order() - 1).rev() {
            strides[i] = strides[i + 1] * usize::from(self[i + 1].signature);
        }

        strides
    }

    fn strides(&self) -> Vec<usize> {
        self.strides_row_major()
    }

    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String> {
        if indices.len() != self.order() {
            return Err("Mismatched order".to_string());
        }

        for (i, &dim_len) in self.iter().map(|slot| &slot.signature).enumerate() {
            if indices[i] >= usize::from(dim_len) {
                return Err(format!(
                    "Index {} out of bounds for dimension {} of size {}",
                    indices[i],
                    i,
                    usize::from(dim_len)
                ));
            }
        }
        Ok(())
    }

    fn flat_index(&self, indices: &[usize]) -> Result<usize, String> {
        let strides = self.strides();
        self.verify_indices(indices)?;

        let mut idx = 0;
        for (i, &index) in indices.iter().enumerate() {
            idx += index * strides[i];
        }
        Ok(idx)
    }

    fn expanded_index(&self, flat_index: usize) -> Result<Vec<usize>, String> {
        let mut indices = vec![];
        let mut index = flat_index;
        for &stride in self.strides().iter() {
            indices.push(index / stride);
            index %= stride;
        }
        if index == 0 {
            Ok(indices)
        } else {
            Err(format!("Index {} out of bounds", flat_index))
        }
    }

    fn size(&self) -> usize {
        self.shape().iter().product()
    }
}

pub trait HasTensorStructure {
    fn structure(&self) -> &Vec<Slot>;
    // inline
    fn order(&self) -> usize {
        self.structure().order()
    }

    fn shape(&self) -> Vec<usize> {
        self.structure().shape()
    }

    fn size(&self) -> usize {
        self.structure().size()
    }

    fn verify_indices(&self, indices: &[usize]) -> Result<(), String> {
        self.structure().verify_indices(indices)
    }

    fn strides(&self) -> Vec<usize> {
        self.structure().strides()
    }

    fn flat_index(&self, indices: &[usize]) -> Result<usize, String> {
        self.structure().flat_index(indices)
    }

    fn expanded_index(&self, flat_index: usize) -> Result<Vec<usize>, String> {
        self.structure().expanded_index(flat_index)
    }

    fn match_index(&self, other: &dyn HasTensorStructure) -> Option<(usize, usize)> {
        self.structure().match_index(other.structure())
    }

    fn traces(&self) -> Vec<Vec<usize>> {
        self.structure().traces()
    }
}

#[derive(Debug, Clone)]
pub struct SparseTensor<T> {
    elements: BTreeMap<Vec<usize>, T>,
    structure: Vec<Slot>,
}

impl<T> HasTensorStructure for SparseTensor<T> {
    fn structure(&self) -> &Vec<Slot> {
        &self.structure
    }
}

impl<T: PartialEq + Default + Clone> SparseTensor<T> {
    pub fn empty(structure: TensorStructure) -> Self {
        SparseTensor {
            elements: BTreeMap::new(),
            structure,
        }
    }

    pub fn empty_from_integers(indices: &[AbstractIndex], dims: &[Dimension]) -> Self {
        let structure = TensorStructure::from_integers(indices, dims);
        SparseTensor {
            elements: BTreeMap::new(),
            structure,
        }
    }

    pub fn from_data(
        data: &[(Vec<ConcreteIndex>, T)],
        indices: &[AbstractIndex],
    ) -> Result<Self, String> {
        let mut dimensions = vec![0; indices.len()];
        for (index, _) in data {
            if index.len() != indices.len() {
                return Err("Mismatched order".to_string());
            }
            for (i, &idx) in index.iter().enumerate() {
                if idx >= dimensions[i] {
                    dimensions[i] = idx + 1;
                }
            }
        }
        Ok(SparseTensor {
            elements: BTreeMap::from_iter(data.iter().cloned()),
            structure: TensorStructure::from_integers(indices, &dimensions),
        })
    }

    pub fn set(&mut self, indices: &[usize], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        self.elements.insert(indices.to_vec(), value);
        Ok(())
    }

    pub fn get(&self, indices: &[usize]) -> Result<Cow<T>, String> {
        self.verify_indices(indices)?;
        // if the index is in the bTree return the value, else return default, lazily allocating the default
        Ok(match self.elements.get(indices) {
            Some(value) => Cow::Borrowed(value),
            None => Cow::Owned(T::default()),
        })
    }
}

pub struct SparseTensorTraceIterator<'a, T> {
    tensor: &'a SparseTensor<T>,
    trace_indices: Vec<Position>,
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T> SparseTensorTraceIterator<'a, T> {
    fn new(tensor: &'a SparseTensor<T>, trace_indices: Vec<Position>) -> Self {
        assert!(trace_indices.len() >= 2, "Invalid trace indices");
        //trace positions must point to the same dimension
        assert!(
            trace_indices
                .iter()
                .map(|&pos| tensor.structure()[pos].signature)
                .collect::<Vec<_>>()
                .iter()
                .all(|&sig| sig == tensor.structure()[trace_indices[0]].signature),
            "Trace indices must point to the same dimension"
        );
        SparseTensorTraceIterator {
            tensor,
            trace_indices,
            current_indices: vec![0; tensor.order()],
            done: false,
        }
    }

    fn increment_indices(&mut self) -> bool {
        for (i, index) in self
            .current_indices
            .iter_mut()
            .enumerate()
            .rev()
            .filter(|(pos, _)| !self.trace_indices.contains(pos))
        {
            *index += 1;
            // If the index goes beyond the shape boundary, wrap around to 0
            if index >= &mut self.tensor.shape()[i] {
                *index = 0;
                continue; // carry over to the next dimension
            }
            return true; // We've successfully found the next combination
        }
        false // No more combinations left
    }
}

impl<'a, T> Iterator for SparseTensorTraceIterator<'a, T>
where
    T: Clone + Copy + Default + std::ops::AddAssign + std::ops::SubAssign + std::cmp::PartialEq,
{
    type Item = (Vec<ConcreteIndex>, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let trace_dimension = self.tensor.structure()[self.trace_indices[0]].signature;
        let trace_sign = trace_dimension.negative();
        let mut trace = T::default();

        for i in 0..trace_dimension.into() {
            let mut indices = self.current_indices.clone();
            for &pos in self.trace_indices.iter() {
                indices[pos] = i;
            }
            if let Some(value) = self.tensor.elements.get(&indices) {
                if trace_sign[i] {
                    trace -= *value;
                } else {
                    trace += *value;
                }
            }
        }

        //make a vector withouth the trace indices
        let trace_indices: Vec<usize> = self
            .current_indices
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(i, _)| !self.trace_indices.contains(&i))
            .map(|(_, x)| x)
            .collect();

        self.done = !self.increment_indices();

        Some((trace_indices, trace))
    }
}

pub struct SparseTensorFiberIterator<'a, T> {
    tensor: &'a SparseTensor<T>,
    fiber_index: Position,
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T> SparseTensorFiberIterator<'a, T> {
    fn new(tensor: &'a SparseTensor<T>, fiber_index: usize) -> Self {
        assert!(fiber_index < tensor.order(), "Invalid fiber index");

        SparseTensorFiberIterator {
            tensor,
            fiber_index,
            current_indices: vec![0; tensor.order()],
            done: false,
        }
    }

    fn increment_indices(&mut self) -> bool {
        for (i, index) in self
            .current_indices
            .iter_mut()
            .enumerate()
            .rev()
            .filter(|(pos, _)| *pos != self.fiber_index)
        {
            *index += 1;
            // If the index goes beyond the shape boundary, wrap around to 0
            if index >= &mut self.tensor.shape()[i] {
                *index = 0;
                continue; // carry over to the next dimension
            }
            return true; // We've successfully foun+155 m / -164 md the next combination
        }
        false // No more combinations left
    }
}

impl<'a, T> Iterator for SparseTensorFiberIterator<'a, T>
where
    T: Clone + Default,
{
    type Item = (Vec<ConcreteIndex>, Vec<Position>, Vec<&'a T>);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let mut lower_bound = self.current_indices.clone();
        let mut upper_bound = self.current_indices.clone();

        // Set the range for the varying dimension to cover all indices
        lower_bound[self.fiber_index] = 0;
        upper_bound[self.fiber_index] = self.tensor.shape()[self.fiber_index];

        let range = self
            .tensor
            .elements
            .range((Included(lower_bound), Included(upper_bound)));

        let mut nonzeros = Vec::new();

        let mut values: Vec<&'a T> = Vec::new();

        for (indices, value) in range.filter(|(key, _)| {
            // Ensure that the difference with start_key is at the same position
            for (i, index) in key.iter().enumerate() {
                if self.fiber_index != i && index != &self.current_indices[i] {
                    return false;
                }
            }
            true
        }) {
            nonzeros.push(indices[self.fiber_index]);
            values.push(value);
        }

        // The upper bound of the range (exclusive)

        // Prepare a vector to hold the combined values
        let fiber_indices = self.current_indices.clone();

        self.done = !self.increment_indices();

        // Check if there are any elements in the range
        if !values.is_empty() {
            Some((fiber_indices, nonzeros, values))
        } else {
            self.next()
        }
    }
}

impl<T: Clone + Default> SparseTensor<T> {
    pub fn iter_fibers(&self, fiber_index: usize) -> SparseTensorFiberIterator<T> {
        SparseTensorFiberIterator::new(self, fiber_index)
    }

    pub fn iter_trace(&self, trace_indices: Vec<Position>) -> SparseTensorTraceIterator<T> {
        SparseTensorTraceIterator::new(self, trace_indices)
    }
}

#[derive(Debug, Clone)]
pub struct DenseTensor<T> {
    data: Vec<T>,
    structure: TensorStructure,
}

impl<T> HasTensorStructure for DenseTensor<T> {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}

impl<T: Default + Clone> DenseTensor<T> {
    pub fn default(structure: TensorStructure) -> Self {
        DenseTensor {
            data: vec![T::default(); structure.size()],
            structure,
        }
    }

    pub fn default_from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self {
        let structure = TensorStructure::from_integers(indices, dims);
        DenseTensor {
            data: vec![T::default(); structure.size()],
            structure,
        }
    }

    pub fn from_data(data: &[T], structure: TensorStructure) -> Result<Self, String> {
        if data.len() != structure.size() {
            return Err("Data length does not match shape".to_string());
        }
        Ok(DenseTensor {
            data: data.to_vec(),
            structure,
        })
    }

    pub fn set(&mut self, indices: &[usize], value: T) {
        let idx = self.flat_index(&indices);
        if let Ok(i) = idx {
            self.data[i] = value;
        }
    }

    pub fn get_linear(&self, index: usize) -> Option<&T> {
        self.data.get(index)
    }

    pub fn get(&self, indices: &[usize]) -> Option<&T> {
        if let Ok(idx) = self.flat_index(&indices) {
            Some(&self.data[idx])
        } else {
            None
        }
    }
}

pub struct TensorIterator<'a, T> {
    tensor: &'a DenseTensor<T>,
    current_flat_index: usize,
}

impl<'a, T> TensorIterator<'a, T> {
    fn new(tensor: &'a DenseTensor<T>) -> Self {
        TensorIterator {
            tensor,
            current_flat_index: 0,
        }
    }
}

impl<'a, T: Clone> Iterator for TensorIterator<'a, T> {
    type Item = (Vec<usize>, T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.tensor.expanded_index(self.current_flat_index) {
            let value = self.tensor.data[self.current_flat_index].clone();

            self.current_flat_index += 1;

            Some((indices, value))
        } else {
            None
        }
    }
}

pub struct DenseTensorTraceIterator<'a, T> {
    tensor: &'a DenseTensor<T>,
    trace_indices: Vec<Position>,
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T> DenseTensorTraceIterator<'a, T> {
    fn new(tensor: &'a DenseTensor<T>, trace_indices: Vec<Position>) -> Self {
        assert!(trace_indices.len() >= 2, "Invalid trace indices");
        //trace positions must point to the same dimension
        assert!(
            trace_indices
                .iter()
                .map(|&pos| tensor.structure()[pos].signature)
                .collect::<Vec<_>>()
                .iter()
                .all(|&sig| sig == tensor.structure()[trace_indices[0]].signature),
            "Trace indices must point to the same dimension"
        );
        DenseTensorTraceIterator {
            tensor,
            trace_indices,
            current_indices: vec![0; tensor.order()],
            done: false,
        }
    }

    fn increment_indices(&mut self) -> bool {
        for (i, index) in self
            .current_indices
            .iter_mut()
            .enumerate()
            .rev()
            .filter(|(pos, _)| !self.trace_indices.contains(pos))
        {
            *index += 1;
            // If the index goes beyond the shape boundary, wrap around to 0
            if index >= &mut self.tensor.shape()[i] {
                *index = 0;
                continue; // carry over to the next dimension
            }
            return true; // We've successfully found the next combination
        }
        false // No more combinations left
    }
}

impl<'a, T> Iterator for DenseTensorTraceIterator<'a, T>
where
    T: Clone + Copy + Default + std::ops::AddAssign + std::ops::SubAssign + std::cmp::PartialEq,
{
    type Item = (Vec<ConcreteIndex>, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let trace_dimension = self.tensor.structure()[self.trace_indices[0]].signature;
        let trace_sign = trace_dimension.negative();
        let mut trace = T::default();

        for i in 0..trace_dimension.into() {
            let mut indices = self.current_indices.clone();
            for &pos in self.trace_indices.iter() {
                indices[pos] = i;
            }
            if trace_sign[i] {
                trace -= *self.tensor.get(&indices).unwrap();
            } else {
                trace += *self.tensor.get(&indices).unwrap();
            }
        }

        //make a vector withouth the trace indices
        let trace_indices: Vec<usize> = self
            .current_indices
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(i, _)| !self.trace_indices.contains(&i))
            .map(|(_, x)| x)
            .collect();

        self.done = !self.increment_indices();

        Some((trace_indices, trace))
    }
}
pub struct DenseTensorFiberIterator<'a, T> {
    tensor: &'a DenseTensor<T>,
    strides: Vec<usize>,
    fixedindex: usize,
    linear_start: usize,
    current_fiber: usize,
    total_fibers: usize,
}

impl<'a, T> DenseTensorFiberIterator<'a, T> {
    fn new(tensor: &'a DenseTensor<T>, fixedindex: usize) -> Self {
        assert!(fixedindex < tensor.order(), "Invalid fixedindex");

        let fiber_length = tensor.shape()[fixedindex];
        let total_fibers = tensor.size() / fiber_length;
        let strides = tensor.strides();

        DenseTensorFiberIterator {
            tensor,
            strides,
            fixedindex,
            linear_start: 0,
            current_fiber: 0,
            total_fibers,
        }
    }

    fn update_linear_start(&mut self) {
        let mut expanded_index = self.tensor.expanded_index(self.linear_start).unwrap();

        for (i, index) in expanded_index
            .iter_mut()
            .enumerate()
            .rev()
            .filter(|(pos, _)| *pos != self.fixedindex)
        {
            *index += 1;
            // If the index goes beyond the shape boundary, wrap around to 0
            if index >= &mut self.tensor.shape()[i] {
                *index = 0;
                continue; // carry over to the next dimension
            }
            break; // We've successfully foun+155 m / -164 md the next combination
        }
        self.linear_start = self.tensor.flat_index(&expanded_index).unwrap();
    }
}

impl<'a, T> Iterator for DenseTensorFiberIterator<'a, T>
where
    T: Clone + Default,
{
    type Item = (Vec<usize>, Vec<&'a T>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_fiber >= self.total_fibers {
            return None;
        }

        // Determine start index for the current fiber

        let mut fiberdata = Vec::with_capacity(self.tensor.shape()[self.fixedindex]);

        let fiberindices = self.tensor.expanded_index(self.linear_start).unwrap();

        for i in 0..self.tensor.shape()[self.fixedindex] {
            let linear_index = self.linear_start + i * self.strides[self.fixedindex];
            fiberdata.push(self.tensor.get_linear(linear_index).unwrap());
        }

        self.update_linear_start();
        self.current_fiber += 1;
        // Determine end index for the current fiber

        Some((fiberindices, fiberdata))
    }
}

impl<T: Clone> DenseTensor<T> {
    // ... [Other methods] ...

    pub fn iter(&self) -> TensorIterator<T> {
        TensorIterator::new(self)
    }

    pub fn iter_fibers(&self, fixedindex: usize) -> DenseTensorFiberIterator<T> {
        DenseTensorFiberIterator::new(self, fixedindex)
    }

    pub fn iter_trace(&self, trace_indices: Vec<Position>) -> DenseTensorTraceIterator<T> {
        DenseTensorTraceIterator::new(self, trace_indices)
    }
}

pub enum NumTensor<T> {
    Dense(DenseTensor<T>),
    Sparse(SparseTensor<T>),
}

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
    pub fn contract_with_sparse(&self, other: &SparseTensor<T>) -> Option<Self> {
        other.contract_with_dense(self)
    }

    pub fn contract_with_dense(&self, other: &Self) -> Option<Self> {
        if let Some((i, j)) = self.match_index(other) {
            // println!("{},{}", i, j);
            let self_shape = self.shape();

            let dimension_of_contraction = self_shape[i];
            let metric = self.structure()[i].signature.negative();

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
    pub fn contract_with_dense(&self, other: &DenseTensor<T>) -> Option<DenseTensor<T>> {
        if let Some((i, j)) = self.match_index(other) {
            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = vec![T::default(); final_structure.size()];

            let metric = self.structure()[i].signature.negative();

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

    pub fn contract_with_sparse(&self, other: &Self) -> Option<Self> {
        if let Some((i, j)) = self.match_index(other) {

            let final_structure = self.structure().merge_at(other.structure(), (i, j));
            let mut result_data = BTreeMap::new();


            let metric = self.structure()[i].signature.negative();

            for (index_a, nonzeros_a, fiber_a) in self.iter_fibers(i) {
                for (index_b, nonzeros_b, fiber_b) in other.iter_fibers(j) {


                    let result_index = index_a[..i]
                        .iter()
                        .chain(&index_a[i + 1..])
                        .chain(&index_b[..j])
                        .chain(&index_b[j + 1..])
                        .cloned()
                        .collect::<Vec<_>>();

                    let mut value = T::default();
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

                    if nonzero {

                        if value != T::default() {
                            result_data.insert(result_index, value);
                        }
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


#[allow(dead_code)]
struct SymbolicTensor {
    structure: TensorStructure,
    expression: String,
}

impl HasTensorStructure for SymbolicTensor {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}

#[allow(dead_code)]
enum Tensor<T> {
    Num(NumTensor<T>),
    Symbolic(SymbolicTensor),
}

#[cfg(test)]
mod tests;
