use std::ops::{AddAssign, Neg, SubAssign};

use super::*;
use ahash::AHashSet;

use permutation::Permutation;

use symbolica::state::{State, Workspace};
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

pub struct TensorSkeletonFiberIterator {
    pub varying_fiber_index: usize,
    pub increment: usize,
    pub stride: Option<usize>,
    pub shift: Option<usize>,
    pub max: usize,
}

impl TensorSkeletonFiberIterator {
    pub fn new<N>(skeleton: &TensorSkeleton<N>, fiber_position: usize) -> Self {
        assert!(fiber_position < skeleton.order(), "Invalid fiber index");

        let strides = skeleton.strides();
        let fiber_stride = strides[fiber_position];
        let dim = skeleton.shape()[fiber_position];

        let max = skeleton.size() - fiber_stride * (dim - 1) - 1;

        let mut stride = None;
        let mut shift = None;
        let mut increment = 1;

        if fiber_position == skeleton.order() - 1 {
            increment = *strides.get(skeleton.order() - 2).unwrap_or(&1);
        } else if fiber_position != 0 {
            shift = Some(strides[fiber_position - 1]);
            stride = Some(strides[fiber_position]);
        }

        TensorSkeletonFiberIterator {
            varying_fiber_index: 0,
            increment,
            stride,
            shift,
            max,
        }
    }

    pub fn reset(&mut self) {
        self.varying_fiber_index = 0;
    }
}

impl Iterator for TensorSkeletonFiberIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.varying_fiber_index > self.max {
            return None;
        }
        let ret = self.varying_fiber_index;

        self.varying_fiber_index += self.increment;

        if let Some(s) = self.stride {
            if self.varying_fiber_index % s == 0 {
                self.varying_fiber_index += self.shift.unwrap() - s;
            }
        }

        Some(ret)
    }
}
#[derive(Debug)]
pub struct TensorSkeletonMultiFiberIterator {
    pub varying_fiber_index: usize,
    pub increment: usize,
    pub fixed_strides: Vec<usize>,
    pub shifts: Vec<usize>,
    pub max: usize,
}

#[derive(Debug)]
pub struct TensorSkeletonMultiFiberMetricIterator {
    pub iterator: TensorSkeletonMultiFiberIterator,
    pos: usize,
    pub reps: Vec<Representation>,
    pub indices: Vec<usize>,
    pub has_neg: AHashSet<usize>,
    map: Vec<usize>,
    permutation: Permutation,
}

impl TensorSkeletonMultiFiberIterator {
    pub fn new<N>(
        skeleton: &TensorSkeleton<N>,
        fiber_positions: &[bool],
    ) -> TensorSkeletonMultiFiberIterator {
        let strides = skeleton.strides();
        let dims = skeleton.shape();
        let order = skeleton.order();
        let mut max = 0;

        // max -= 1;

        let mut increment = 1;

        let mut fixed_strides = vec![];
        let mut shifts = vec![];

        let mut before = true;
        let mut has_seen_stride = false;
        let mut first = true;

        for pos in (0..order).rev() {
            let is_fixed = fiber_positions[pos];

            if is_fixed && !before && !first {
                has_seen_stride = true;
                fixed_strides.push(strides[pos]);
            }
            if !is_fixed && before && has_seen_stride {
                shifts.push(strides[pos]);
            }

            if !is_fixed {
                max += (dims[pos] - 1) * strides[pos];
                if first {
                    increment = strides[pos];
                    first = false;
                }
            }

            before = is_fixed;
        }

        if fixed_strides.len() > shifts.len() {
            fixed_strides.pop();
        }

        TensorSkeletonMultiFiberIterator {
            varying_fiber_index: 0,
            increment,
            fixed_strides,
            shifts,
            max,
        }
    }

    pub fn new_conjugate<N>(
        skeleton: &TensorSkeleton<N>,
        fiber_positions: &[bool],
    ) -> (Self, Self) {
        let strides = skeleton.strides();
        let dims = skeleton.shape();
        let order = skeleton.order();
        let mut max = 0;

        // max -= 1;

        let mut increment = 1;

        // println!("{:?}{}", fiber_positions, increment);
        let mut fixed_strides = vec![];
        let mut fixed_strides_conj = vec![];
        let mut shifts = vec![];
        let mut shifts_conj = vec![];

        let mut before = true;
        let mut has_seen_stride = false;
        let mut has_seen_stride_conj = false;
        let mut first = true;
        let mut first_conj = true;
        let mut increment_conj = 1;

        let mut max_conj = 0;

        for pos in (0..order).rev() {
            let is_fixed = fiber_positions[pos];

            if is_fixed && !before {
                if !first {
                    has_seen_stride = true;
                    fixed_strides.push(strides[pos]);
                }

                if has_seen_stride_conj {
                    shifts_conj.push(strides[pos]);
                }
            }
            if !is_fixed && before {
                if has_seen_stride {
                    shifts.push(strides[pos]);
                }

                if !first_conj {
                    fixed_strides_conj.push(strides[pos]);
                    has_seen_stride_conj = true;
                }
            }

            if is_fixed {
                max_conj += (dims[pos] - 1) * strides[pos];
                if first_conj {
                    increment_conj = strides[pos];
                    first_conj = false;
                }
            } else {
                max += (dims[pos] - 1) * strides[pos];
                if first {
                    increment = strides[pos];
                    first = false;
                }
            }

            before = is_fixed;
        }

        if fixed_strides.len() > shifts.len() {
            fixed_strides.pop();
        }

        if fixed_strides_conj.len() > shifts_conj.len() {
            fixed_strides_conj.pop();
        }

        (
            TensorSkeletonMultiFiberIterator {
                varying_fiber_index: 0,
                increment: increment_conj,
                fixed_strides: fixed_strides_conj,
                shifts: shifts_conj,
                max: max_conj,
            },
            TensorSkeletonMultiFiberIterator {
                varying_fiber_index: 0,
                increment,
                fixed_strides,
                shifts,
                max,
            },
        )
    }

    pub fn reset(&mut self) {
        self.varying_fiber_index = 0;
    }
}

impl Iterator for TensorSkeletonMultiFiberIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.varying_fiber_index > self.max {
            return None;
        }
        let ret = self.varying_fiber_index;

        self.varying_fiber_index += self.increment;

        for (i, s) in self.fixed_strides.iter().enumerate() {
            if self.varying_fiber_index % s == 0 {
                self.varying_fiber_index += self.shifts[i] - s;
            } else {
                break;
            }
        }

        Some(ret)
    }
}

impl TensorSkeletonMultiFiberMetricIterator {
    pub fn new<N>(
        skeleton: &TensorSkeleton<N>,
        fiber_positions: &[bool],
        permutation: Permutation,
    ) -> TensorSkeletonMultiFiberMetricIterator {
        // for f in fiber_positions {
        //     filter[*f] = true;
        // }
        let iterator = TensorSkeletonMultiFiberIterator::new(skeleton, fiber_positions);

        let mut f = fiber_positions.iter();
        let mut reps = skeleton.reps();
        reps.retain(|_| !*f.next().unwrap());
        let capacity = reps.iter().map(usize::from).product();
        let indices = vec![0; reps.len()];
        // println!("Reps : {:?}", reps);

        TensorSkeletonMultiFiberMetricIterator {
            iterator,
            reps,
            indices,
            has_neg: AHashSet::with_capacity(capacity),
            pos: 0,
            permutation,
            map: Vec::new(),
        }
    }

    pub fn new_conjugates<N>(
        skeleton: &TensorSkeleton<N>,
        fiber_positions: &[bool],
        permutation: Permutation,
    ) -> (Self, TensorSkeletonMultiFiberIterator) {
        let iters = TensorSkeletonMultiFiberIterator::new_conjugate(skeleton, fiber_positions);
        let mut f = fiber_positions.iter();
        let mut reps = skeleton.reps();
        reps.retain(|_| *f.next().unwrap());
        let capacity = reps.iter().map(usize::from).product();
        let indices = vec![0; reps.len()];
        // println!("Reps : {:?}", reps);
        (
            TensorSkeletonMultiFiberMetricIterator {
                iterator: iters.0,
                reps,
                indices,
                has_neg: AHashSet::with_capacity(capacity),
                pos: 0,
                permutation,
                map: Vec::new(),
            },
            iters.1,
        )
    }

    fn increment_indices(&mut self) {
        let mut carry = true;
        for (i, r) in self.indices.iter_mut().rev().zip(self.reps.iter().rev()) {
            if carry {
                *i += 1;
                carry = *i == usize::from(r);
                *i %= usize::from(r);
            }
        }
    }

    fn has_neg(&mut self, i: usize) -> bool {
        if self.has_neg.get(&i).is_some() {
            return true;
        }
        let mut neg = false;
        for (i, r) in self
            .indices
            .iter()
            .zip(self.reps.iter())
            .filter(|(_, r)| !matches!(r, Representation::Euclidean(_)))
        {
            neg ^= r.is_neg(*i);
        }

        // let neg = self
        //     .reps
        //     .last()
        //     .unwrap()
        //     .is_neg(*self.indices.last().unwrap());
        if neg {
            self.has_neg.insert(i);
        }
        neg
    }

    fn generate_map(&mut self) {
        // println!("{:?}{:?}", self.map.len(), self.pos);
        if self.map.len() <= self.pos {
            let mut stride = 1;
            // println!("{:?}", self.indices);
            // println!("{:?}", self.reps);
            // println!("{:?}", self.permutation);
            let mut ind = 0;
            for (i, d) in self
                .permutation
                .apply_slice(self.indices.as_slice())
                .iter()
                .zip(self.permutation.apply_slice(self.reps.as_slice()).iter())
                .rev()
            {
                ind += stride * i;
                stride *= usize::from(d);
            }
            self.map.push(ind);
            // println!("{:?}", self.map)
        }
    }

    pub fn reset(&mut self) -> &[usize] {
        self.iterator.reset();
        self.indices = vec![0; self.reps.len()];
        self.pos = 0;
        // println!("reset!");
        &self.map
    }
}

impl Iterator for TensorSkeletonMultiFiberMetricIterator {
    type Item = (usize, usize, bool);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.iterator.next() {
            let pos = self.pos;
            self.generate_map();
            self.pos += 1;
            let neg = self.has_neg(i);
            self.increment_indices();
            Some((pos, i, neg))
        } else {
            None
        }
    }
}

#[test]
fn construct() {
    let a = TensorSkeleton::from_idxsing(
        &[
            (1, Representation::Euclidean(2)), //0
            (3, Representation::Euclidean(2)), //1     inc
                                               // (33, Representation::Lorentz(4)),  //2
                                               // (23, Representation::Lorentz(3)),  //3
                                               // (22, Representation::Lorentz(3)),  //4     inc
                                               // (35, Representation::Lorentz(3)),  //5
                                               // (42, Representation::Lorentz(3)),  //6     inc
        ],
        "t",
    );

    let fixed = [false, true];
    let free: Vec<bool> = fixed.iter().map(|x| !x).collect();
    let fi = 0;

    let fiber = TensorSkeletonFiberIterator::new(&a, fi);

    let fiber_iter = TensorSkeletonMultiFiberIterator::new(&a, &free);
    let mut free_iter = TensorSkeletonMultiFiberIterator::new(&a, &fixed);

    println!("{:?}", free_iter);
    println!("{:?}", fiber_iter);

    println!("single fiber");
    for f in fiber {
        // println!("{:?}", a.expanded_index(f));
        for i in 0..2 {
            println!("{:?}", a.expanded_index(f + i * a.strides()[fi]));
        }
    }

    println!("multi fiber");
    for f in fiber_iter {
        for i in free_iter.by_ref() {
            println!("{:?}", a.expanded_index(f + i));
        }
        free_iter.reset();
    }

    let iters = TensorSkeletonMultiFiberIterator::new_conjugate(&a, &fixed);
    println!("{:?}", iters.0);
    println!("{:?}", iters.1);
    let t = TensorSkeletonMultiFiberIterator::new(&a, &fixed);
    // let tt = TensorSkeletonMultiFiberMetricIterator::new(&a, &fixed);
    println!("{:?}", a.strides());
    // println!("{:?}", a.strides_column_major());
    println!("{:?}", t.fixed_strides);
    println!("{:?}", t.shifts);
    println!("{:?}", a.expanded_index(t.max));
    println!("{:?}", t.increment);
    // println!("{:?}", a.order());

    println!("{:?}", a.expanded_index(t.max));

    let _reps = a.reps();

    // for f in tt {
    //     println!("{:?},{:?}", a.expanded_index(f.1), f.2);
    // }
}

pub struct TensorFiberIterator<'a, T, N>
where
    T: HasTensorStructure<Name = N>,
{
    tensor: &'a T,
    fiber_iter: TensorSkeletonFiberIterator,
    skipped: usize,
    pub fiber_dimension: usize,
    increment: usize,
}

impl<'a, T, N> TensorFiberIterator<'a, T, N>
where
    T: HasTensorStructure<Name = N>,
{
    pub fn new(tensor: &'a T, fiber_position: usize) -> Self {
        let fiber_iter = TensorSkeletonFiberIterator::new(tensor.structure(), fiber_position);
        let increment = tensor.strides()[fiber_position];

        TensorFiberIterator {
            tensor,
            fiber_iter,
            skipped: 0,
            fiber_dimension: tensor.shape()[fiber_position],
            increment,
        }
    }

    pub fn reset(&mut self) {
        self.fiber_iter.reset();
        self.skipped = 0;
    }
}

impl<'a, T, N> Iterator for TensorFiberIterator<'a, SparseTensor<T, N>, N>
where
    N: Clone,
{
    type Item = (usize, Vec<ConcreteIndex>, Vec<&'a T>);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.fiber_iter.next() {
            let mut out = Vec::new();
            let mut nonzeros = Vec::new();
            for i in 0..self.fiber_dimension {
                if let Some(v) = self.tensor.elements.get(&(s + i * self.increment)) {
                    nonzeros.push(i);
                    out.push(v);
                }
            }
            if !out.is_empty() {
                let skipped = self.skipped;
                self.skipped = 0;
                Some((skipped, nonzeros, out))
            } else {
                self.skipped += 1;
                self.next()
            }
        } else {
            None
        }
    }
}

impl<'a, T, N> Iterator for TensorFiberIterator<'a, DenseTensor<T, N>, N>
where
    N: Clone,
{
    type Item = Vec<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(s) = self.fiber_iter.next() {
            let mut out = Vec::with_capacity(self.fiber_dimension);
            for i in 0..self.fiber_dimension {
                if let Some(v) = self.tensor.get_linear(s + i * self.increment) {
                    out.push(v);
                }
            }

            Some(out)
        } else {
            None
        }
    }
}

pub struct TensorMultiFiberMetricIterator<'a, T, N>
where
    T: HasTensorStructure<Name = N>,
{
    tensor: &'a T,
    fiber_iter: TensorSkeletonMultiFiberMetricIterator,
    free_iter: TensorSkeletonMultiFiberIterator,
    skipped: usize,
    pub map: Vec<usize>,
    capacity: usize,
}

impl<'a, T, N> TensorMultiFiberMetricIterator<'a, T, N>
where
    T: HasTensorStructure<Name = N>,
{
    pub fn new(tensor: &'a T, fiber_positions: &[bool], permutation: Permutation) -> Self {
        let iters = TensorSkeletonMultiFiberMetricIterator::new_conjugates(
            tensor.structure(),
            fiber_positions,
            permutation,
        );

        // let free: Vec<bool> = fiber_positions.iter().map(|x| !x).collect();#

        // let free: Vec<bool> = fiber_positions.iter().map(|x| !x).collect();

        // let fiber_iter = TensorSkeletonMultiFiberMetricIterator::new(tensor.structure(), &free);
        // let free_iter = TensorSkeletonMultiFiberIterator::new(tensor.structure(), fiber_positions);

        let mut f = fiber_positions.iter();
        let mut reps = tensor.shape();
        reps.retain(|_| !*f.next().unwrap());
        let capacity = reps.iter().product();
        TensorMultiFiberMetricIterator {
            tensor,
            map: vec![],
            fiber_iter: iters.0,
            free_iter: iters.1,
            skipped: 0,
            capacity,
        }
    }

    pub fn reset(&mut self) {
        self.fiber_iter.reset();
        self.free_iter.reset();
        self.skipped = 0;
    }
}

impl<'a, T, N> Iterator for TensorMultiFiberMetricIterator<'a, SparseTensor<T, N>, N>
where
    N: Clone,
{
    type Item = (usize, Vec<ConcreteIndex>, Vec<(&'a T, bool)>);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.free_iter.next() {
            let mut out = Vec::new();
            let mut nonzeros = Vec::new();
            for (pos, ind, met) in self.fiber_iter.by_ref() {
                if let Some(v) = self.tensor.get_linear(ind + i) {
                    out.push((v, met));
                    nonzeros.push(pos);
                    // println!("hi");
                }
            }
            if self.map.is_empty() {
                self.map = self.fiber_iter.reset().to_owned();
            } else {
                self.fiber_iter.reset();
            }
            if !out.is_empty() {
                let skipped = self.skipped;
                self.skipped = 0;
                Some((skipped, nonzeros, out))
            } else {
                self.skipped += 1;
                self.next()
            }
        } else {
            None
        }
    }
}

impl<'a, T, N> Iterator for TensorMultiFiberMetricIterator<'a, DenseTensor<T, N>, N>
where
    N: Clone,
{
    type Item = Vec<(&'a T, bool)>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.free_iter.next() {
            let mut out = Vec::with_capacity(self.capacity);
            for (_, ind, met) in self.fiber_iter.by_ref() {
                if let Some(v) = self.tensor.get_linear(ind + i) {
                    out.push((v, met));
                }
            }
            if self.map.is_empty() {
                self.map = self.fiber_iter.reset().to_owned();
            } else {
                self.fiber_iter.reset();
            }

            self.fiber_iter.reset();
            if out.is_empty() {
                return self.next();
            }
            Some(out)
        } else {
            None
        }
    }
}

pub struct TensorMultiFiberIterator<'a, T, N>
where
    T: HasTensorStructure<Name = N>,
{
    tensor: &'a T,
    fiber_iter: TensorSkeletonMultiFiberIterator,
    free_iter: TensorSkeletonMultiFiberIterator,
    skipped: usize,
}

impl<'a, T, N> TensorMultiFiberIterator<'a, T, N>
where
    T: HasTensorStructure<Name = N>,
{
    pub fn new(tensor: &'a T, fiber_positions: &[bool]) -> Self {
        let iters =
            TensorSkeletonMultiFiberIterator::new_conjugate(tensor.structure(), fiber_positions);
        TensorMultiFiberIterator {
            tensor,
            fiber_iter: iters.0,
            free_iter: iters.1,
            skipped: 0,
        }
    }

    pub fn reset(&mut self) {
        self.fiber_iter.reset();
        self.free_iter.reset();
        self.skipped = 0;
    }
}

impl<'a, T, N> Iterator for TensorMultiFiberIterator<'a, SparseTensor<T, N>, N>
where
    N: Clone,
{
    type Item = (usize, Vec<ConcreteIndex>, Vec<&'a T>);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.free_iter.next() {
            let mut out = Vec::new();
            let mut nonzeros = Vec::new();
            for (pos, ind) in self.fiber_iter.by_ref().enumerate() {
                if let Some(v) = self.tensor.get_linear(ind + i) {
                    out.push(v);
                    nonzeros.push(pos);
                }
            }
            self.fiber_iter.reset();
            if !out.is_empty() {
                let skipped = self.skipped;
                self.skipped = 0;
                Some((skipped, nonzeros, out))
            } else {
                self.skipped += 1;
                self.next()
            }
        } else {
            None
        }
    }
}

impl<'a, T, N> Iterator for TensorMultiFiberIterator<'a, DenseTensor<T, N>, N>
where
    N: Clone,
{
    type Item = Vec<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.free_iter.next() {
            let mut out = Vec::new();
            for ind in self.fiber_iter.by_ref() {
                if let Some(v) = self.tensor.get_linear(ind + i) {
                    out.push(v);
                }
            }
            self.fiber_iter.reset();
            Some(out)
        } else {
            None
        }
    }
}

pub struct SparseTensorIterator<'a, T, N> {
    iter: std::collections::hash_map::Iter<'a, usize, T>,
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
            let indices = self.structure.expanded_index(*k).unwrap();
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
        let left_stride = tensor.strides().get(fiber_position - 1).copied();
        let right_stride = tensor.strides().get(fiber_position + 1).copied();
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
            (Some(l), Some(_r)) => {
                self.fiber_index += 1;
                if self.fiber_index % self.fiber_stride == 0 {
                    self.fiber_index += l - self.fiber_stride;
                }
            }
            (Some(l), None) => {
                self.fiber_index += l;
            }
            (None, Some(_r)) => {
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
                .get(&(self.fiber_index + i * self.fiber_stride))
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

    pub fn iter_fiber(&self, fiber_index: usize) -> TensorFiberIterator<Self, I> {
        TensorFiberIterator::new(self, fiber_index)
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

    pub fn iter_multi_fibers_metric(
        &self,
        fiber_positions: &[bool],
        permutation: Permutation,
    ) -> TensorMultiFiberMetricIterator<Self, I> {
        TensorMultiFiberMetricIterator::new(self, fiber_positions, permutation)
    }

    pub fn iter_multi_fibers(&self, fiber_positions: &[bool]) -> TensorMultiFiberIterator<Self, I> {
        TensorMultiFiberIterator::new(self, fiber_positions)
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

impl<T, I> DenseTensor<T, I> {
    // ... [Other methods] ...

    pub fn iter(&self) -> DenseTensorIterator<T, I> {
        DenseTensorIterator::new(self)
    }

    pub fn iter_fiber(&self, fixedindex: usize) -> TensorFiberIterator<Self, I> {
        TensorFiberIterator::new(self, fixedindex)
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

    pub fn iter_multi_fibers_metric(
        &self,
        fiber_positions: &[bool],
        permutation: Permutation,
    ) -> TensorMultiFiberMetricIterator<Self, I> {
        TensorMultiFiberMetricIterator::new(self, fiber_positions, permutation)
    }

    pub fn iter_multi_fibers(&self, fiber_positions: &[bool]) -> TensorMultiFiberIterator<Self, I> {
        TensorMultiFiberIterator::new(self, fiber_positions)
    }
}
