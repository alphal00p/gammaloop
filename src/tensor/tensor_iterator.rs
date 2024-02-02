use std::ops::{AddAssign, Neg, SubAssign};

use serde_yaml::mapping::Iter;
use symbolica::state::{State, Workspace};

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

pub struct SparseTensorIterator<'a, T, N> {
    iter: intmap::iter::Iter<'a, u64, T>,
    structure: &'a TensorSkeleton<N>,
}

impl<'a, T, N> SparseTensorIterator<'a, T, N> {
    fn new(tensor: &'a SparseTensor<T, N>) -> Self {
        SparseTensorIterator {
            iter: tensor.elements.iter(),
            structure: &tensor.structure,
        }
    }
}

impl<'a, T, N> Iterator for SparseTensorIterator<'a, T, N> {
    type Item = (Vec<ConcreteIndex>, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((k, v)) = self.iter.next() {
            let indices = self.structure.expanded_index(*k as usize).unwrap();
            Some((indices, v))
        } else {
            None
        }
    }
}

// impl<'a, T, I> IntoIterator for &'a SparseTensor<T, I> {
//     type Item = (&'a Vec<ConcreteIndex>, &'a T);
//     type IntoIter = SparseTensorIterator<'a, T>;

//     fn into_iter(self) -> Self::IntoIter {
//         SparseTensorIterator::new(self)
//     }
// }

pub struct SparseTensorTraceIterator<'a, T, I> {
    tensor: &'a SparseTensor<T, I>,
    trace_indices: [usize; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T, I> SparseTensorTraceIterator<'a, T, I> {
    fn new(tensor: &'a SparseTensor<T, I>, trace_indices: [usize; 2]) -> Self {
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

impl<'a, T, I> Iterator for SparseTensorTraceIterator<'a, T, I>
where
    T: for<'c> std::ops::AddAssign<&'c T>
        + for<'b> std::ops::SubAssign<&'b T>
        + std::ops::Neg<Output = T>
        + Clone,
    I: Clone,
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

        // Data might not exist at that concrete index usize, we advance it till it does, and if not we skip

        while self.tensor.is_empty_at(&indices) {
            let Some((i, signint)) = iter.next() else {
                self.done = !self.increment_indices();
                return self.next(); // skip
            };
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            sign = signint;
        }

        let value = (*self.tensor.get(&indices).unwrap()).clone(); //Should now be safe to unwrap
        let mut trace = if *sign { value.neg() } else { value };

        for (i, sign) in iter {
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            if let Ok(value) = self.tensor.get(&indices) {
                if *sign {
                    trace -= value;
                } else {
                    trace += value;
                }
            }
        }

        //make a vector withouth the trace indices
        let trace_indices: Vec<ConcreteIndex> = self
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

pub struct SparseTensorSymbolicTraceIterator<'a, 'b, T, I> {
    tensor: &'a SparseTensor<T, I>,
    trace_indices: [usize; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
    state: &'b State,
    ws: &'b Workspace,
}

impl<'a, 'b, T, I> SparseTensorSymbolicTraceIterator<'a, 'b, T, I> {
    fn new(
        tensor: &'a SparseTensor<T, I>,
        trace_indices: [usize; 2],
        state: &'b State,
        ws: &'b Workspace,
    ) -> Self {
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
        SparseTensorSymbolicTraceIterator {
            tensor,
            trace_indices,
            current_indices: vec![0; tensor.order()],
            done: false,
            state,
            ws,
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

impl<'a, 'b, T, I> Iterator for SparseTensorSymbolicTraceIterator<'a, 'b, T, I>
where
    T: for<'c> SymbolicAddAssign<&'c T> + for<'d> SymbolicSubAssign<&'d T> + SymbolicNeg + Clone,
    I: Clone,
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

        // Data might not exist at that concrete index usize, we advance it till it does, and if not we skip

        while self.tensor.is_empty_at(&indices) {
            let Some((i, signint)) = iter.next() else {
                self.done = !self.increment_indices();
                return self.next(); // skip
            };
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            sign = signint;
        }

        let value = (*self.tensor.get(&indices).unwrap()).clone(); //Should now be safe to unwrap
        let mut trace = if *sign {
            value.neg_sym(self.ws, self.state)
        } else {
            value
        };

        for (i, sign) in iter {
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            if let Ok(value) = self.tensor.get(&indices) {
                if *sign {
                    trace.sub_assign_sym(value, self.ws, self.state);
                } else {
                    trace.add_assign_sym(value, self.ws, self.state);
                }
            }
        }

        //make a vector withouth the trace indices
        let trace_indices: Vec<ConcreteIndex> = self
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

pub struct SparseTensorFiberIterator<'a, T, I> {
    tensor: &'a SparseTensor<T, I>,
    fiber_index: usize,
    fiber_stride: usize,
    fiber_dimension: usize,
    left_stride: Option<usize>,
    right_stride: Option<usize>,
    skipped: usize,
    max: usize,
}

impl<'a, T, I> SparseTensorFiberIterator<'a, T, I> {
    fn new(tensor: &'a SparseTensor<T, I>, fiber_position: usize) -> Self {
        assert!(fiber_position < tensor.order(), "Invalid fiber index");

        let fiber_stride = tensor.strides()[fiber_position];
        let left_stride = tensor.strides().get(fiber_position - 1).map(|x| *x);
        let right_stride = tensor.strides().get(fiber_position + 1).map(|x| *x);
        let max = tensor.size() - fiber_stride * (tensor.shape()[fiber_position] - 1);
        let fiber_dimension = tensor.shape()[fiber_position];

        SparseTensorFiberIterator {
            tensor,
            fiber_index: 0,
            fiber_stride,
            fiber_dimension,
            left_stride,
            right_stride,
            skipped: 0,
            max,
        }
    }

    fn update_linear_start(&mut self) {
        match (self.left_stride, self.right_stride) {
            (Some(l), Some(r)) => {
                self.fiber_index += 1;
                if self.fiber_index % self.fiber_stride == 0 {
                    self.fiber_index += l - self.fiber_stride;
                }
            }
            (Some(l), None) => {
                self.fiber_index += l;
            }
            (None, Some(r)) => {
                self.fiber_index += 1;
            }
            (None, None) => {
                self.fiber_index += 1;
            }
        }
    }

    // fn increment_indices(&mut self) -> bool {
    //     for (i, index) in self
    //         .current_indices
    //         .iter_mut()
    //         .enumerate()
    //         .rev()
    //         .filter(|(pos, _)| *pos != self.fiber_index)
    //     {
    //         *index += 1;
    //         // If the index goes beyond the shape boundary, wrap around to 0
    //         if index >= &mut self.tensor.shape()[i] {
    //             *index = 0;
    //             continue; // carry over to the next dimension
    //         }
    //         return true; // We've successfully found the next combination
    //     }
    //     false // No more combinations left
    // }
}

impl<'a, T, I> Iterator for SparseTensorFiberIterator<'a, T, I>
where
    T: Clone,
    I: Clone,
{
    type Item = (usize, Vec<usize>, Vec<&'a T>);
    fn next(&mut self) -> Option<Self::Item> {
        if self.fiber_index >= self.max {
            return None;
        }

        // let range = self
        //     .tensor
        //     .elements
        //     .range((Included(lower_bound), Included(upper_bound)));

        let mut nonzeros = Vec::new();

        let mut values: Vec<&'a T> = Vec::new();

        // for (indices, value) in range.filter(|(key, _)| {
        //     // Ensure that the difference with start_key is at the same usize
        //     for (i, index) in key.iter().enumerate() {
        //         if self.fiber_index != i && index != &self.current_indices[i] {
        //             return false;
        //         }
        //     }
        //     true
        // }) {
        //     nonzeros.push(indices[self.fiber_index]);
        //     values.push(value);
        // }

        for i in 0..self.fiber_dimension {
            if let Some(v) = self
                .tensor
                .elements
                .get((self.fiber_index + i * self.fiber_stride) as u64)
            {
                nonzeros.push(i);
                values.push(v);
            }
        }
        self.update_linear_start();

        // Check if there are any elements in the range
        if !values.is_empty() {
            let skipped = self.skipped;
            self.skipped = 0;
            Some((skipped, nonzeros, values))
        } else {
            self.skipped += 1;
            self.next()
        }
    }
}

impl<T, I> SparseTensor<T, I> {
    pub fn iter_fibers(&self, fiber_index: usize) -> SparseTensorFiberIterator<T, I> {
        SparseTensorFiberIterator::new(self, fiber_index)
    }

    pub fn iter_trace(&self, trace_indices: [usize; 2]) -> SparseTensorTraceIterator<T, I> {
        SparseTensorTraceIterator::new(self, trace_indices)
    }

    pub fn iter(&self) -> SparseTensorIterator<T, I> {
        SparseTensorIterator::new(self)
    }

    pub fn iter_symbolic_trace<'a, 'b>(
        &'a self,
        trace_indices: [usize; 2],
        state: &'b State,
        ws: &'b Workspace,
    ) -> SparseTensorSymbolicTraceIterator<'a, 'b, T, I> {
        SparseTensorSymbolicTraceIterator::new(self, trace_indices, state, ws)
    }
}

pub struct DenseTensorIterator<'a, T, I> {
    tensor: &'a DenseTensor<T, I>,
    current_flat_index: usize,
}

impl<'a, T, I> DenseTensorIterator<'a, T, I> {
    fn new(tensor: &'a DenseTensor<T, I>) -> Self {
        DenseTensorIterator {
            tensor,
            current_flat_index: 0,
        }
    }
}

impl<'a, T, I> Iterator for DenseTensorIterator<'a, T, I> {
    type Item = (Vec<ConcreteIndex>, &'a T);

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

impl<'a, T, I> IntoIterator for &'a DenseTensor<T, I> {
    type Item = (Vec<ConcreteIndex>, &'a T);
    type IntoIter = DenseTensorIterator<'a, T, I>;

    fn into_iter(self) -> Self::IntoIter {
        DenseTensorIterator::new(self)
    }
}

impl<T, I> IntoIterator for DenseTensor<T, I> {
    type Item = (Vec<ConcreteIndex>, T);
    type IntoIter = DenseTensorIntoIterator<T, I>;

    fn into_iter(self) -> Self::IntoIter {
        DenseTensorIntoIterator::new(self)
    }
}

pub struct DenseTensorIntoIterator<T, I> {
    tensor: DenseTensor<T, I>,
    current_flat_index: usize,
}

impl<T, I> DenseTensorIntoIterator<T, I> {
    fn new(tensor: DenseTensor<T, I>) -> Self {
        DenseTensorIntoIterator {
            tensor,
            current_flat_index: 0,
        }
    }
}

impl<T, I> Iterator for DenseTensorIntoIterator<T, I> {
    type Item = (Vec<ConcreteIndex>, T);

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

pub struct DenseTensorTraceIterator<'a, T, I> {
    tensor: &'a DenseTensor<T, I>,
    trace_indices: [usize; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T, I> DenseTensorTraceIterator<'a, T, I> {
    fn new(tensor: &'a DenseTensor<T, I>, trace_indices: [usize; 2]) -> Self {
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

impl<'a, T, I> Iterator for DenseTensorTraceIterator<'a, T, I>
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
        let trace_indices: Vec<ConcreteIndex> = self
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

pub struct DenseTensorSymbolicTraceIterator<'a, 'b, T, I> {
    tensor: &'a DenseTensor<T, I>,
    trace_indices: [usize; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
    state: &'b State,
    ws: &'b Workspace,
}

impl<'a, 'b, T, I> DenseTensorSymbolicTraceIterator<'a, 'b, T, I> {
    fn new(
        tensor: &'a DenseTensor<T, I>,
        trace_indices: [usize; 2],
        state: &'b State,
        ws: &'b Workspace,
    ) -> Self {
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
        DenseTensorSymbolicTraceIterator {
            tensor,
            trace_indices,
            current_indices: vec![0; tensor.order()],
            done: false,
            state,
            ws,
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

impl<'a, 'b, T, I> Iterator for DenseTensorSymbolicTraceIterator<'a, 'b, T, I>
where
    T: Clone
        + for<'c> SymbolicAddAssign<&'c T>
        + for<'d> SymbolicSubAssign<&'d T>
        + SymbolicNeg
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

        let mut trace = if *sign {
            value.neg_sym(self.ws, self.state)
        } else {
            value
        };

        for (i, sign) in iter {
            for &pos in self.trace_indices.iter() {
                indices[pos] = i;
            }

            if let Some(value) = self.tensor.get(&indices) {
                if *sign {
                    trace.sub_assign_sym(value, self.ws, self.state);
                } else {
                    trace.add_assign_sym(value, self.ws, self.state);
                }
            }
        }

        //make a vector without the trace indices
        let trace_indices: Vec<ConcreteIndex> = self
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

pub struct DenseTensorFiberIterator<'a, T, I> {
    tensor: &'a DenseTensor<T, I>,
    fiber_stride: usize,
    left_stride: Option<usize>,
    right_stride: Option<usize>,
    fiber_indices: Vec<usize>,
    fixedindex: usize,
    linear_start: usize,
    current_fiber: usize,
    total_fibers: usize,
    max: usize,
    done: bool,
}

impl<'a, T, I> DenseTensorFiberIterator<'a, T, I> {
    fn new(tensor: &'a DenseTensor<T, I>, fixedindex: usize) -> Self {
        assert!(fixedindex < tensor.order(), "Invalid fixedindex");

        let fiber_length = tensor.shape()[fixedindex];
        let total_fibers = tensor.size() / fiber_length;
        let fiber_stride = tensor.strides()[fixedindex];
        let left_stride = tensor.strides().get(fixedindex - 1).map(|x| *x);
        let right_stride = tensor.strides().get(fixedindex + 1).map(|x| *x);
        let max = tensor.size() - fiber_stride * (fiber_length - 1);

        DenseTensorFiberIterator {
            tensor,
            fiber_stride,
            left_stride,
            right_stride,
            fiber_indices: vec![0; tensor.order()],
            fixedindex,
            linear_start: 0,
            current_fiber: 0,
            total_fibers,
            max,
            done: 0 == max,
        }
    }

    fn update_linear_start(&mut self) {
        match (self.left_stride, self.right_stride) {
            (Some(l), Some(r)) => {
                self.linear_start += 1;
                if self.linear_start % self.fiber_stride == 0 {
                    self.linear_start += l - self.fiber_stride;
                }
            }
            (Some(l), None) => {
                self.linear_start += l;
            }
            (None, Some(r)) => {
                self.linear_start += 1;
            }
            (None, None) => {
                self.linear_start += 1;
            }
        }
    }
}

impl<'a, T, I> Iterator for DenseTensorFiberIterator<'a, T, I> {
    type Item = Vec<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Determine start index for the current fiber

        let mut fiberdata = Vec::with_capacity(self.tensor.shape()[self.fixedindex]);

        for i in 0..self.tensor.shape()[self.fixedindex] {
            let linear_index = self.linear_start + i * self.fiber_stride;
            fiberdata.push(self.tensor.get_linear(linear_index).unwrap());
        }

        self.update_linear_start();
        self.done = self.linear_start >= self.max;

        // Determine end index for the current fiber

        Some(fiberdata)
    }
}

impl<T, I> DenseTensor<T, I> {
    // ... [Other methods] ...

    pub fn iter(&self) -> DenseTensorIterator<T, I> {
        DenseTensorIterator::new(self)
    }

    pub fn iter_fibers(&self, fixedindex: usize) -> DenseTensorFiberIterator<T, I> {
        DenseTensorFiberIterator::new(self, fixedindex)
    }

    pub fn iter_trace(&self, trace_indices: [usize; 2]) -> DenseTensorTraceIterator<T, I> {
        DenseTensorTraceIterator::new(self, trace_indices)
    }

    pub fn iter_symbolic_trace<'a, 'b>(
        &'a self,
        trace_indices: [usize; 2],
        state: &'b State,
        ws: &'b Workspace,
    ) -> DenseTensorSymbolicTraceIterator<'a, 'b, T, I> {
        DenseTensorSymbolicTraceIterator::new(self, trace_indices, state, ws)
    }
}
