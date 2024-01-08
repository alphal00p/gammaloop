use std::ops::{AddAssign, Neg, SubAssign};

use super::*;

pub struct TensorStructureIndexIterator<'a> {
    structure: &'a TensorStructure,
    current_flat_index: usize,
}

impl<'a> Iterator for TensorStructureIndexIterator<'a> {
    type Item = Vec<ConcreteIndex>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.structure.expanded_index(self.current_flat_index) {
            self.current_flat_index += 1;

            Some(indices)
        } else {
            None
        }
    }
}

impl<'a> TensorStructureIndexIterator<'a> {
    pub fn new(structure: &'a TensorStructure) -> Self {
        TensorStructureIndexIterator {
            structure,
            current_flat_index: 0,
        }
    }
}

pub struct SparseTensorIterator<'a, T> {
    iter: std::collections::btree_map::Iter<'a, Vec<usize>, T>,
}

impl<'a, T> SparseTensorIterator<'a, T> {
    fn new(tensor: &'a SparseTensor<T>) -> Self {
        SparseTensorIterator {
            iter: tensor.elements.iter(),
        }
    }
}

impl<'a, T> Iterator for SparseTensorIterator<'a, T> {
    type Item = (&'a Vec<usize>, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl<'a, T> IntoIterator for &'a SparseTensor<T> {
    type Item = (&'a Vec<usize>, &'a T);
    type IntoIter = SparseTensorIterator<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        SparseTensorIterator::new(self)
    }
}

pub struct SparseTensorTraceIterator<'a, T> {
    tensor: &'a SparseTensor<T>,
    trace_indices: [Position; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T> SparseTensorTraceIterator<'a, T> {
    fn new(tensor: &'a SparseTensor<T>, trace_indices: [Position; 2]) -> Self {
        //trace positions must point to the same dimension
        assert!(
            trace_indices
                .iter()
                .map(|&pos| tensor.structure()[pos].representation)
                .collect::<Vec<_>>()
                .iter()
                .all(|&sig| sig == tensor.structure()[trace_indices[0]].representation),
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
        // Filter out the trace indices
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
    T: for<'c> std::ops::AddAssign<&'c T>
        + for<'b> std::ops::SubAssign<&'b T>
        + std::ops::Neg<Output = T>
        + Clone,
{
    type Item = (Vec<ConcreteIndex>, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let trace_dimension = self.tensor.structure()[self.trace_indices[0]].representation;
        let trace_sign = trace_dimension.negative();
        let mut iter = trace_sign.iter().enumerate();
        let mut indices = self.current_indices.clone();
        let (i, mut sign) = iter.next().unwrap(); //First element (to eliminate the need for default)

        indices[self.trace_indices[0]] = i;
        indices[self.trace_indices[1]] = i;

        // Data might not exist at that concrete index position, we advance it till it does, and if not we skip

        while self.tensor.is_empty_at(&indices) {
            let Some((i, signint)) = iter.next() else {
                self.done = !self.increment_indices();
                return self.next(); // skip
            };
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            sign = signint;
        }

        let value = (*self.tensor.elements.get(&indices).unwrap()).clone(); //Should now be safe to unwrap
        let mut trace = if *sign { value.neg() } else { value };

        for (i, sign) in iter {
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            if let Some(value) = self.tensor.elements.get(&indices) {
                if *sign {
                    trace -= value;
                } else {
                    trace += value;
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

impl<'a, T> Iterator for SparseTensorFiberIterator<'a, T> {
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

impl<T> SparseTensor<T> {
    pub fn iter_fibers(&self, fiber_index: usize) -> SparseTensorFiberIterator<T> {
        SparseTensorFiberIterator::new(self, fiber_index)
    }

    pub fn iter_trace(&self, trace_indices: [Position; 2]) -> SparseTensorTraceIterator<T> {
        SparseTensorTraceIterator::new(self, trace_indices)
    }
}
impl<T> SparseTensor<T> {
    pub fn iter(&self) -> SparseTensorIterator<T> {
        SparseTensorIterator::new(self)
    }
}

pub struct DenseTensorIterator<'a, T> {
    tensor: &'a DenseTensor<T>,
    current_flat_index: usize,
}

impl<'a, T> DenseTensorIterator<'a, T> {
    fn new(tensor: &'a DenseTensor<T>) -> Self {
        DenseTensorIterator {
            tensor,
            current_flat_index: 0,
        }
    }
}

impl<'a, T> Iterator for DenseTensorIterator<'a, T> {
    type Item = (Vec<usize>, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.tensor.expanded_index(self.current_flat_index) {
            let value = self.tensor.get_linear(self.current_flat_index).unwrap();

            self.current_flat_index += 1;

            Some((indices, value))
        } else {
            None
        }
    }
}

impl<'a, T> IntoIterator for &'a DenseTensor<T> {
    type Item = (Vec<usize>, &'a T);
    type IntoIter = DenseTensorIterator<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        DenseTensorIterator::new(self)
    }
}

impl<T> IntoIterator for DenseTensor<T> {
    type Item = (Vec<usize>, T);
    type IntoIter = DenseTensorIntoIterator<T>;

    fn into_iter(self) -> Self::IntoIter {
        DenseTensorIntoIterator::new(self)
    }
}

pub struct DenseTensorIntoIterator<T> {
    tensor: DenseTensor<T>,
    current_flat_index: usize,
}

impl<T> DenseTensorIntoIterator<T> {
    fn new(tensor: DenseTensor<T>) -> Self {
        DenseTensorIntoIterator {
            tensor,
            current_flat_index: 0,
        }
    }
}

impl<T> Iterator for DenseTensorIntoIterator<T> {
    type Item = (Vec<usize>, T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.tensor.expanded_index(self.current_flat_index) {
            let indices = indices.clone();
            let value = self.tensor.data.remove(self.current_flat_index);

            self.current_flat_index += 1;

            Some((indices, value))
        } else {
            None
        }
    }
}

pub struct DenseTensorTraceIterator<'a, T> {
    tensor: &'a DenseTensor<T>,
    trace_indices: [Position; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T> DenseTensorTraceIterator<'a, T> {
    fn new(tensor: &'a DenseTensor<T>, trace_indices: [Position; 2]) -> Self {
        assert!(trace_indices.len() >= 2, "Invalid trace indices");
        //trace positions must point to the same dimension
        assert!(
            trace_indices
                .iter()
                .map(|&pos| tensor.structure()[pos].representation)
                .collect::<Vec<_>>()
                .iter()
                .all(|&sig| sig == tensor.structure()[trace_indices[0]].representation),
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
    T: for<'c> AddAssign<&'c T>
        + for<'b> SubAssign<&'b T>
        + Neg<Output = T>
        + Clone
        + std::fmt::Debug,
{
    type Item = (Vec<ConcreteIndex>, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let trace_dimension = self.tensor.structure()[self.trace_indices[0]].representation;
        let trace_sign = trace_dimension.negative();

        let mut iter = trace_sign.iter().enumerate();
        let mut indices = self.current_indices.clone();
        let (_, sign) = iter.next().unwrap(); //First sign

        for &pos in self.trace_indices.iter() {
            indices[pos] = 0;
        }

        let value = self.tensor.get(&indices).unwrap().clone();

        let mut trace = if *sign { value.neg() } else { value };

        for (i, sign) in iter {
            for &pos in self.trace_indices.iter() {
                indices[pos] = i;
            }

            if let Some(value) = self.tensor.get(&indices) {
                if *sign {
                    trace -= value;
                } else {
                    trace += value;
                }
            }
        }

        //make a vector without the trace indices
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
        let mut expanded_index = self
            .tensor
            .expanded_index(self.linear_start)
            .unwrap()
            .clone();

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

impl<'a, T> Iterator for DenseTensorFiberIterator<'a, T> {
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

impl<T> DenseTensor<T> {
    // ... [Other methods] ...

    pub fn iter(&self) -> DenseTensorIterator<T> {
        DenseTensorIterator::new(self)
    }

    pub fn iter_fibers(&self, fixedindex: usize) -> DenseTensorFiberIterator<T> {
        DenseTensorFiberIterator::new(self, fixedindex)
    }

    pub fn iter_trace(&self, trace_indices: [Position; 2]) -> DenseTensorTraceIterator<T> {
        DenseTensorTraceIterator::new(self, trace_indices)
    }
}
