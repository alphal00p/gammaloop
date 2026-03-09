//! # Permutations
//!
//! This module provides a `Permutation` struct and associated functionalities
//! for representing and working with permutations of a sequence of integers `0..n`.
//!
//! ## Key Features:
//!
//! - **Representation**: A `Permutation` is stored by its direct mapping (`map[i]` is
//!   the image of `i`) and its inverse mapping.
//! - **Construction**:
//!   - Identity permutation: `Permutation::id(n)`.
//!   - From a mapping vector: `Permutation::from_map(vec![...])`.
//!   - From disjoint cycles: `Permutation::from_disjoint_cycles(&[vec![...]])`.
//!   - From potentially overlapping cycles with specified composition order:
//!     `Permutation::from_cycles_ordered(&[vec![...]], order)`.
//! - **Basic Operations**:
//!   - Inverse: `p.inverse()`.
//!   - Composition: `p1.compose(&p2)` (applies `p2` then `p1`).
//!   - Apply to slices: `p.apply_slice(data)` (returns a new `Vec`),
//!     `p.apply_slice_in_place(data_mut)`.
//!   - Power: `p.pow(k)`.
//!   - Sign: `p.sign()` (+1 for even, -1 for odd).
//!   - Check for identity: `p.is_identity()`.
//! - **Cycle Utilities**:
//!   - Find cycle decomposition: `p.find_cycles()`.
//!   - Convert to transpositions: `p.transpositions()`.
//! - **Iterators**:
//!   - `p.iter_slice(data)`: Iterates over `data` in the order defined by `p.map`.
//!   - `p.iter_slice_inv(data)`: Iterates over `data` in the order defined by `p.inv`.
//!   - Mutable versions: `p.iter_slice_mut(data_mut)`, `p.iter_slice_inv_mut(data_mut)`.
//! - **Ranking and Unranking**:
//!   - Myrvold & Ruskey's ranking/unranking algorithms (`myrvold_ruskey_rank1`, `myrvold_ruskey_unrank1`, etc.).
//! - **Sorting**:
//!   - `Permutation::sort(slice)`: Returns the permutation that sorts the slice.
//! - **Graph Integration**:
//!   - The `HedgeGraphExt` and `PermutationExt` traits suggest integration with
//!     graph structures, likely for tasks like finding graph automorphisms or
//!     canonical labellings by permuting vertices/hedges.
//!
//! The module is designed to be comprehensive for permutation-related tasks,
//! including those relevant to combinatorial algorithms and graph theory.

use std::{
    collections::BTreeMap, fmt, iter::FusedIterator, marker::PhantomData, ops::Index, ptr::NonNull,
};

use crate::half_edge::{
    involution::Hedge,
    subgraph::{BaseSubgraph, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
    HedgeGraph, NodeIndex,
};
use ahash::AHashSet;
use thiserror::Error;

use crate::half_edge::involution::Flow;

/// A permutation of `0..n`, with the ability to apply itself (or its inverse) to slices.
///
/// # Examples
///
/// ```
/// use linnet::permutation::Permutation;
///
/// // Create a permutation that maps 0->2, 1->0, 2->1, 3->3
/// let p = Permutation::from_map(vec![2, 0, 1, 3]);
///
/// // Apply the permutation to a slice
/// let data = vec![10, 20, 30, 40];
/// let permuted = p.apply_slice(&data);
/// assert_eq!(permuted, vec![20, 30, 10, 40]);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "bincode", derive(bincode::Encode, bincode::Decode))]
pub struct Permutation {
    map: Vec<usize>,
    inv: Vec<usize>,
}

/// Implement ordering comparisons for permutations based on their `map` field.
impl PartialOrd for Permutation {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.map.partial_cmp(&other.map)
    }
}

impl Permutation {
    // --------------------------------------------------------------------------------------------
    // Basic Constructors and Accessors
    // --------------------------------------------------------------------------------------------

    /// Creates the identity permutation of length `n`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::id(4);
    /// assert_eq!(p.apply_slice(&[10,20,30,40]), vec![10,20,30,40]);
    /// ```
    pub fn id(n: usize) -> Self {
        Permutation {
            map: (0..n).collect(),
            inv: (0..n).collect(),
        }
    }

    /// Creates a permutation from a mapping vector.
    /// The `map` vector states where index `i` is sent: `map[i]` is the image of `i`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// assert_eq!(p.apply_slice(&[10,20,30]), vec![20,30,10]);
    /// ```
    pub fn from_map(map: Vec<usize>) -> Self {
        let mut inv = vec![0; map.len()];
        for (i, &j) in map.iter().enumerate() {
            inv[j] = i;
        }
        Permutation { map, inv }
    }

    /// Creates a permutation from a inverse mapping vector.
    /// The `inv` vector states that index `i` is actually inv[i]
    pub fn from_inv(inv: Vec<usize>) -> Self {
        let mut map = vec![0; inv.len()];
        for (i, &j) in inv.iter().enumerate() {
            map[j] = i;
        }
        Permutation { map, inv }
    }

    /// Creates a permutation from (from,to) pairs, minimizing swaps among all valid completions.
    pub fn from_mappings<I>(mappings: I, size: usize) -> Result<Self, String>
    where
        I: IntoIterator<Item = (usize, usize)>,
    {
        if size == 0 {
            return Ok(Self::id(0));
        }

        // Identity to start from
        let mut map: Vec<usize> = (0..size).collect();

        // forced_sources[s] = true iff s has an explicitly forced mapping
        let mut forced_sources = vec![false; size];

        // owner_of_target[t] = Some(s) if *any* assignment currently owns t (forced or identity we keep)
        let mut owner_of_target: Vec<Option<usize>> = vec![None; size];

        // forced_owner_of_target[t] = Some(s) iff a *forced* mapping s -> t exists
        let mut forced_owner_of_target: Vec<Option<usize>> = vec![None; size];

        // 1) Apply & validate forced mappings in O(1) each
        for (from, to) in mappings {
            if from >= size || to >= size {
                return Err(format!(
                    "Index out of bounds: mapping ({from}, {to}) for size {size}",
                ));
            }
            if let Some(prev_from) = forced_owner_of_target[to] {
                if prev_from != from {
                    return Err(format!(
                        "Target {to} is mapped from multiple sources ({prev_from} and {from})",
                    ));
                }
            }
            if forced_sources[from] && map[from] != to {
                return Err(format!(
                    "Source {} has conflicting mappings ({} vs {})",
                    from, map[from], to
                ));
            }

            map[from] = to;
            forced_sources[from] = true;
            forced_owner_of_target[to] = Some(from);
            owner_of_target[to] = Some(from);
        }

        // 2) Keep identities whenever possible; collect displaced (unforced whose identity is taken)
        let mut displaced_sources: Vec<usize> = Vec::new();
        for s in 0..size {
            if forced_sources[s] {
                continue; // fixed by input
            }
            if owner_of_target[s].is_none() {
                // identity slot s is free -> keep s->s
                map[s] = s;
                owner_of_target[s] = Some(s);
            } else {
                // identity slot s is already taken by a forced mapping -> s is displaced
                displaced_sources.push(s);
            }
        }

        // 3) For each displaced source s, assign it to the unique free target reachable
        //    by walking through forced owners: t0 = s; while forced_owner_of_target[t] = Some(u) { t = u; }
        //    Use a small memo to avoid re-walking shared tails.
        let mut sink_cache: Vec<Option<usize>> = vec![None; size];

        for s in displaced_sources {
            let mut path: Vec<usize> = Vec::new();
            let mut t = s;

            // Walk until we hit a cached sink or a free target
            loop {
                if let Some(sink) = sink_cache[t] {
                    // propagate cached sink to the path and assign
                    for x in path.drain(..) {
                        sink_cache[x] = Some(sink);
                    }
                    map[s] = sink;
                    owner_of_target[sink] = Some(s);
                    break;
                }
                path.push(t);

                match forced_owner_of_target[t] {
                    Some(owner) => {
                        // advance along the forced-owner chain
                        t = owner;
                        continue;
                    }
                    None => {
                        // t is a free target (no forced owner) -> this is the sink
                        for x in path.drain(..) {
                            sink_cache[x] = Some(t);
                        }
                        map[s] = t;
                        owner_of_target[t] = Some(s);
                        break;
                    }
                }
            }
        }

        // 4) Final uniqueness check
        let mut seen = vec![false; size];
        for &to in &map {
            if to >= size || std::mem::replace(&mut seen[to], true) {
                return Err("Failed to create valid permutation".to_string());
            }
        }

        Ok(Self::from_map(map))
    }

    /// Creates a permutation from a list of sources mapping to the prefix targets 0..k-1.
    ///
    /// The j-th element of `sources` is forced to map to target `j`.
    /// Remaining indices are chosen to (1) keep identities when possible and
    /// (2) minimize the number of swaps (maximize the number of cycles).
    ///
    /// # Arguments
    /// * `sources` - iterator of distinct indices in `0..size`, length = k
    /// * `size`    - total permutation size (n)
    ///
    /// # Returns
    /// A swaps-optimal completion satisfying `sources[j] -> j`.
    ///
    /// # Example
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_prefix_sources([5, 0, 3], 6).unwrap();
    /// // forced: 5->0, 0->1, 3->2  (targets 0,1,2)
    /// assert_eq!(p[5], 0);
    /// assert_eq!(p[0], 1);
    /// assert_eq!(p[3], 2);
    /// // others chosen swaps-optimally
    /// ```
    pub fn from_prefix_sources<I>(sources: I, size: usize) -> Result<Self, String>
    where
        I: IntoIterator<Item = usize>,
    {
        if size == 0 {
            return Ok(Self::id(0));
        }

        // Collect sources and validate in O(n+k)
        let src_vec: Vec<usize> = sources.into_iter().collect();
        let k = src_vec.len();
        if k > size {
            return Err(format!("Too many sources: k={k} > size={size}"));
        }

        let mut seen_source = vec![false; size];
        for &s in &src_vec {
            if s >= size {
                return Err(format!("Source {s} out of bounds for size {size}"));
            }
            if std::mem::replace(&mut seen_source[s], true) {
                return Err(format!("Duplicate source index {s}"));
            }
        }

        // Identity base
        let mut map: Vec<usize> = (0..size).collect();

        // Bookkeeping
        let mut forced_sources = vec![false; size]; // s is in src_vec
        let mut owner_of_target: Vec<Option<usize>> = vec![None; size]; // who currently owns target t

        // forced_owner_of_target[t] = Some(s) only for t in 0..k with s = src_vec[t]
        // (This defines the fixed "owner graph" we walk to find sinks.)
        let mut forced_owner_of_target: Vec<Option<usize>> = vec![None; size];

        // 1) Apply forced prefix mappings: src_vec[j] -> j
        for (t, &s) in (0..k).zip(src_vec.iter()) {
            // note: t is guaranteed < size since k <= size
            if let Some(prev) = owner_of_target[t] {
                // Should never happen here because each prefix target appears exactly once
                if prev != s {
                    return Err(format!("Target {t} already owned by {prev}"));
                }
            }
            map[s] = t;
            forced_sources[s] = true;
            owner_of_target[t] = Some(s);
            forced_owner_of_target[t] = Some(s);
        }

        // 2) Keep identities whenever possible; collect displaced
        let mut displaced: Vec<usize> = Vec::new();
        for i in 0..size {
            if forced_sources[i] {
                continue;
            }
            if owner_of_target[i].is_none() {
                // identity target i is free -> keep i->i
                map[i] = i;
                owner_of_target[i] = Some(i);
            } else {
                // identity target taken by some forced mapping -> i is displaced
                displaced.push(i);
            }
        }

        // 3) Route each displaced source to its unique sink via forced-owner chains.
        //    Cache sinks so shared tails are walked once.
        let mut sink_cache: Vec<Option<usize>> = vec![None; size];

        for s in displaced {
            let mut path: Vec<usize> = Vec::new();
            let mut t = s;

            loop {
                if let Some(cached) = sink_cache[t] {
                    // Assign cached sink; mark ownership
                    for x in path.drain(..) {
                        sink_cache[x] = Some(cached);
                    }
                    // Safety: cached is free by construction of earlier assignments
                    map[s] = cached;
                    owner_of_target[cached] = Some(s);
                    break;
                }

                path.push(t);

                match forced_owner_of_target[t] {
                    Some(owner) => {
                        // Follow the forced owner chain
                        t = owner;
                    }
                    None => {
                        // We reached a target not in the forced prefix (or outside 0..k)
                        // This is the unique sink for the whole component.
                        // It is guaranteed free at this point because:
                        // - If t == i for some identity i kept earlier, that i cannot lie on the
                        //   forced-owner chain (no forced owner on t), hence we would not be
                        //   displacing s into an already-owned t in this component.
                        for x in path.drain(..) {
                            sink_cache[x] = Some(t);
                        }
                        map[s] = t;
                        owner_of_target[t] = Some(s);
                        break;
                    }
                }
            }
        }

        // 4) Verify permutation validity
        let mut seen = vec![false; size];
        for &to in &map {
            if to >= size || std::mem::replace(&mut seen[to], true) {
                return Err("Construction failed: duplicate or OOB target".into());
            }
        }

        Ok(Self::from_map(map))
    }

    /// Returns the internal mapping as a slice.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// assert_eq!(p.map(), &[2, 0, 1]);
    /// ```
    // -- ADDED
    pub fn map(&self) -> &[usize] {
        &self.map
    }

    /// Returns the inverse mapping as a slice.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// assert_eq!(p.inv(), &[1, 2, 0]);
    /// ```
    // -- ADDED
    pub fn inv(&self) -> &[usize] {
        &self.inv
    }

    // --------------------------------------------------------------------------------------------
    // Basic Operations
    // --------------------------------------------------------------------------------------------

    /// Returns the inverse of the permutation.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let inv = p.inverse();
    /// assert_eq!(inv.apply_slice(&[10,20,30]), vec![30, 10, 20]);
    /// ```
    pub fn inverse(&self) -> Self {
        Permutation {
            map: self.inv.clone(),
            inv: self.map.clone(),
        }
    }

    /// Applies `self` to a slice, returning a new `Vec<T>` in permuted order.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let data = vec![10, 20, 30];
    /// assert_eq!(p.apply_slice(&data), vec![20, 30, 10]);
    /// ```
    pub fn apply_slice<T: Clone, S>(&self, slice: S) -> Vec<T>
    where
        S: AsRef<[T]>,
    {
        let s = slice.as_ref();
        self.inv.iter().map(|&idx| s[idx].clone()).collect()
    }

    /// Applies the inverse of `self` to a slice, returning a new `Vec<T>` in permuted order.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let data = vec![10, 20, 30];
    /// assert_eq!(p.apply_slice_inv(&data), vec![30, 10, 20]);
    /// ```
    pub fn apply_slice_inv<T: Clone, S>(&self, slice: S) -> Vec<T>
    where
        S: AsRef<[T]>,
    {
        let s = slice.as_ref();
        self.map.iter().map(|&idx| s[idx].clone()).collect()
    }

    /// Applies `self` in-place to the provided slice by using transpositions
    /// derived from the cycle decomposition.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let mut data = vec![10, 20, 30];
    /// p.apply_slice_in_place(&mut data);
    /// assert_eq!(data, vec![20, 30, 10]);
    /// ```
    pub fn apply_slice_in_place<T: Clone, S>(&self, slice: &mut S)
    where
        S: AsMut<[T]>,
    {
        let transpositions = self.transpositions();
        for (i, j) in transpositions.iter().rev() {
            slice.as_mut().swap(*i, *j);
        }
    }

    /// Applies the inverse of `self` in-place to the provided slice by using transpositions
    /// derived from the cycle decomposition.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let mut data = vec![10, 20, 30];
    /// p.apply_slice_in_place_inv(&mut data);
    /// assert_eq!(data, vec![30, 10, 20]);
    /// ```
    pub fn apply_slice_in_place_inv<T, S>(&self, slice: &mut S)
    where
        S: AsMut<[T]>,
    {
        let transpositions = self.transpositions();
        for (i, j) in transpositions {
            slice.as_mut().swap(i, j);
        }
    }

    /// Composes `self` with another permutation `other`, returning a new permutation:
    /// `(self ◦ other)(i) = self.map[other.map[i]]`.
    pub fn compose(&self, other: &Self) -> Self {
        let map = other.map.iter().map(|&i| self.map[i]).collect();
        Self::from_map(map)
    }

    // --------------------------------------------------------------------------------------------
    // Iterating
    // --------------------------------------------------------------------------------------------

    /// Returns an iterator that yields references to the elements of the slice
    /// in the order specified by the permutation.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let data = vec![10, 20, 30];
    /// let mut iter = p.iter_slice(&data);
    /// assert_eq!(iter.next(), Some(&30));
    /// assert_eq!(iter.next(), Some(&10));
    /// assert_eq!(iter.next(), Some(&20));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter_slice<'a, T>(&'a self, slice: &'a [T]) -> PermutationMapIter<'a, T> {
        PermutationMapIter::new(slice, &self.map)
    }

    /// Returns an iterator that yields references to the elements of the slice
    /// in the order specified by the inverse of the permutation.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]); // map: 0->2, 1->0, 2->1
    ///                                             // inv: 0->1, 1->2, 2->0
    /// let data = vec![10, 20, 30]; // p.apply_slice_inv(&data) would be [20, 30, 10]
    /// let mut iter = p.iter_slice_inv(&data);
    /// assert_eq!(iter.next(), Some(&20)); // data[inv[0]] = data[1]
    /// assert_eq!(iter.next(), Some(&30)); // data[inv[1]] = data[2]
    /// assert_eq!(iter.next(), Some(&10)); // data[inv[2]] = data[0]
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter_slice_inv<'a, T>(&'a self, slice: &'a [T]) -> PermutationMapIter<'a, T> {
        PermutationMapIter::new(slice, &self.inv)
    }

    /// Returns a mutable iterator that yields mutable references to the elements of the slice
    /// in the order specified by the permutation's `map`.
    ///
    /// The iterator's length is `self.len()`.
    ///
    /// # Panics
    ///
    /// Panics if `slice.len()` is less than `self.len()`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let mut data = vec![10, 20, 30];
    /// for (i, val_ref) in p.iter_slice_mut(&mut data).enumerate() {
    ///     *val_ref += i * 100; // 0->30, 1->10, 2->20
    /// }
    /// // data was permuted to [10,20,30] then modified in permuted order
    /// // Original indices: data[map[0]=2] (30) gets +0 -> 30
    /// //                   data[map[1]=0] (10) gets +100 -> 110
    /// //                   data[map[2]=1] (20) gets +200 -> 220
    /// assert_eq!(data, vec![110, 220, 30]);
    /// ```
    pub fn iter_slice_mut<'a, T>(&'a self, slice: &'a mut [T]) -> PermutationMapIterMut<'a, T> {
        PermutationMapIterMut::new(slice, &self.map)
    }

    /// Returns a mutable iterator that yields mutable references to the elements of the slice
    /// in the order specified by the permutation's `inv` (inverse) map.
    ///
    /// The iterator's length is `self.len()`.
    ///
    /// # Panics
    ///
    /// Panics if `slice.len()` is less than `self.len()`.
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]); // inv: [1, 2, 0]
    /// let mut data = vec![10, 20, 30];
    /// for (i, val_ref) in p.iter_slice_inv_mut(&mut data).enumerate() {
    ///     *val_ref += i * 100; // 0->data[1]=20, 1->data[2]=30, 2->data[0]=10
    /// }
    /// // Original indices: data[inv[0]=1] (20) gets +0 -> 20
    /// //                   data[inv[1]=2] (30) gets +100 -> 130
    /// //                   data[inv[2]=0] (10) gets +200 -> 210
    /// assert_eq!(data, vec![210, 20, 130]);
    /// ```
    pub fn iter_slice_inv_mut<'a, T>(&'a self, slice: &'a mut [T]) -> PermutationMapIterMut<'a, T> {
        PermutationMapIterMut::new(slice, &self.inv)
    }

    // --------------------------------------------------------------------------------------------
    // Sorting Utilities
    // --------------------------------------------------------------------------------------------

    /// Given a slice of items that implement `Ord`, returns the permutation that sorts them
    /// in ascending order.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let data = vec![30, 10, 20, 40];
    /// let perm = Permutation::sort(&data);
    /// // perm.map should be [1, 2, 0, 3]
    /// assert_eq!(perm.apply_slice(&data), vec![10, 20, 30, 40]);
    /// ```
    pub fn sort<T, S>(slice: S) -> Permutation
    where
        T: Ord,
        S: AsRef<[T]>,
    {
        let s = slice.as_ref();
        let mut permutation: Vec<usize> = (0..s.len()).collect();
        permutation.sort_by_key(|&i| &s[i]);
        Self::from_inv(permutation)
    }

    pub fn sort_by_key<T, S, F, O>(slice: S, mut key: F) -> Permutation
    where
        O: Ord,
        S: AsRef<[T]>,
        F: FnMut(&T) -> O,
    {
        let s = slice.as_ref();
        let mut permutation: Vec<usize> = (0..s.len()).collect();
        permutation.sort_by_key(|&i| key(&s[i]));
        Self::from_inv(permutation)
    }

    // --------------------------------------------------------------------------------------------
    // Cycles and Transpositions
    // --------------------------------------------------------------------------------------------

    /// Returns the cycle decomposition of `self` as a `Vec` of cycles,
    /// each cycle represented as a `Vec<usize>`.
    /// Each cycle lists the indices of a single cycle, e.g. `[0, 2, 1]` means `0->2, 2->1, 1->0`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1, 3]);
    /// let cycles = p.find_cycles();
    /// // cycles might be [[0, 2, 1], [3]]
    /// assert_eq!(cycles.len(), 2);
    /// ```
    pub fn find_cycles(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.map.len()];
        let mut cycles = Vec::new();
        for i in 0..self.map.len() {
            if visited[i] {
                continue;
            }
            let mut cycle = Vec::new();
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                cycle.push(j);
                j = self.map[j];
            }
            if !cycle.is_empty() {
                cycles.push(cycle);
            }
        }
        cycles
    }

    /// Converts a single cycle to a list of transpositions that produce that cycle.
    /// This is a helper method and typically not used directly.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let cycle = vec![0, 2, 1];
    /// let transpositions = Permutation::cycle_to_transpositions(&cycle);
    /// // cycle 0->2,2->1,1->0 can be built from swaps (0,1) and (0,2)
    /// assert_eq!(transpositions, vec![(0, 1), (0, 2)]);
    /// ```
    pub fn cycle_to_transpositions(cycle: &[usize]) -> Vec<(usize, usize)> {
        let mut transpositions = Vec::new();
        for i in (1..cycle.len()).rev() {
            transpositions.push((cycle[0], cycle[i]));
        }
        transpositions
    }

    /// Returns the list of transpositions for `self`, by decomposing it into cycles
    /// and then converting each cycle to transpositions.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 0, 1]);
    /// let transpositions = p.transpositions();
    /// assert_eq!(transpositions, vec![(0, 1), (0, 2)]);
    /// ```
    pub fn transpositions(&self) -> Vec<(usize, usize)> {
        let cycles = self.find_cycles();
        let mut transpositions = Vec::new();
        for cycle in cycles {
            transpositions.extend(Self::cycle_to_transpositions(&cycle));
        }
        transpositions
    }

    // --------------------------------------------------------------------------------------------
    // Myrvold & Ruskey Ranking/Unranking
    // --------------------------------------------------------------------------------------------

    /// Computes the rank of the permutation in the Myrvold & Ruskey "Rank1" ordering.
    ///
    /// This is a recursive implementation. For permutations of size `n`, it removes
    /// the position of `n-1` from the permutation, multiplies the result by `n`,
    /// and adds the index of `n-1`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 1, 3, 0]);
    /// assert_eq!(p.myrvold_ruskey_rank1(), 12);
    /// ```
    pub fn myrvold_ruskey_rank1(mut self) -> usize {
        let n = self.map.len();
        if self.map.len() == 1 {
            return 0;
        }

        let s = self.map[n - 1];
        self.map.swap_remove(self.inv[n - 1]);
        self.inv.swap_remove(s);

        s + n * self.myrvold_ruskey_rank1()
    }

    /// Unranks a permutation of size `n` from its Myrvold & Ruskey "Rank1" index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::myrvold_ruskey_unrank1(4, 12);
    /// assert_eq!(p.map(), &[2, 1, 3, 0]);
    /// ```
    pub fn myrvold_ruskey_unrank1(n: usize, mut rank: usize) -> Self {
        let mut p = (0..n).collect::<Vec<_>>();
        for i in (1..=n).rev() {
            let j = rank % i;
            rank /= i;
            p.swap(i - 1, j);
        }
        Permutation::from_map(p)
    }

    fn factorial(n: usize) -> usize {
        (1..=n).product()
    }

    /// Computes the rank of the permutation in the Myrvold & Ruskey "Rank2" ordering.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![2, 1, 3, 0]);
    /// // Suppose it has rank = R. We can test we get p back by unranking R.
    /// let rank = p.clone().myrvold_ruskey_rank2();
    /// let q = Permutation::myrvold_ruskey_unrank2(4, rank);
    /// assert_eq!(q, p);
    /// ```
    pub fn myrvold_ruskey_rank2(mut self) -> usize {
        let n = self.map.len();
        if n == 1 {
            return 0;
        }
        let s = self.map[n - 1];
        self.map.swap_remove(self.inv[n - 1]);
        self.inv.swap_remove(s);
        s * Self::factorial(n - 1) + self.myrvold_ruskey_rank2()
    }

    /// Unranks a permutation of size `n` from its Myrvold & Ruskey "Rank2" index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::myrvold_ruskey_unrank2(4, 1);
    /// assert_eq!(p.map(), &[2, 1, 3, 0]);
    /// ```
    pub fn myrvold_ruskey_unrank2(n: usize, mut rank: usize) -> Self {
        let mut p = (0..n).collect::<Vec<_>>();
        for i in (1..=n).rev() {
            let s = rank / (Self::factorial(i - 1));
            p.swap(i - 1, s);
            rank %= Self::factorial(i - 1);
        }
        Permutation::from_map(p)
    }

    // --------------------------------------------------------------------------------------------
    // Additional Suggested Methods
    // --------------------------------------------------------------------------------------------

    /// Checks if this permutation is the identity permutation (i.e., does nothing).
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::id(4);
    /// assert!(p.is_identity());
    ///
    /// let q = Permutation::from_map(vec![1,0,2,3]);
    /// assert!(!q.is_identity());
    /// ```
    // -- ADDED
    pub fn is_identity(&self) -> bool {
        self.map.iter().enumerate().all(|(i, &m)| i == m)
    }

    /// Returns the sign (+1 or -1) of the permutation,
    /// indicating whether it is an even (+1) or odd (-1) permutation.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![1,0,3,2]);
    /// assert_eq!(p.sign(), 1); // even
    ///
    /// let q = Permutation::from_map(vec![2,1,0]);
    /// assert_eq!(q.sign(), -1); // odd
    /// ```
    // -- ADDED
    pub fn sign(&self) -> i8 {
        // Count inversions or use transpositions
        let mut sign = 1i8;
        for cycle in self.find_cycles() {
            // Each cycle of length k contributes (k-1) to the total parity
            let k = cycle.len();
            if k > 1 && (k - 1) % 2 == 1 {
                sign = -sign;
            }
        }
        sign
    }

    /// Computes the k-th power of the permutation (composition with itself k times).
    /// For k = 0, it returns the identity of the same size.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let p = Permutation::from_map(vec![1, 2, 0]);
    /// // p^2 maps 0->p(1)=2, 1->p(2)=0, 2->p(0)=1 => [2,0,1]
    /// let p2 = p.pow(2);
    /// assert_eq!(p2.map(), &[2, 0, 1]);
    /// ```
    // -- ADDED
    pub fn pow(&self, k: usize) -> Self {
        if k == 0 {
            return Permutation::id(self.map.len());
        }
        let mut result = Permutation::id(self.map.len());
        let mut base = self.clone();
        let mut exp = k;

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.compose(&base);
            }
            base = base.compose(&base);
            exp /= 2;
        }
        result
    }
}

/// Specifies the direction for reading cycle compositions
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum CycleOrder {
    /// Apply rightmost cycles first: (a b)(c d) applies (c d) then (a b)
    LastFirst,
    /// Apply leftmost cycles first: (a b)(c d) applies (a b) then (c d)
    FirstFirst,
}

impl fmt::Display for Permutation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // First show cycle notation
        let cycles = self.find_cycles();
        let mut first = true;
        for cycle in cycles {
            if cycle.len() > 1 {
                // Only show non-trivial cycles
                if !first {
                    write!(f, " ")?;
                }
                write!(f, "(")?;
                for (i, &x) in cycle.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{x}")?;
                }
                write!(f, ")")?;
                first = false;
            }
        }
        if first {
            // If no cycles were printed (identity permutation)
            write!(f, "()")?;
        }

        // Then show one-line notation
        write!(f, " [")?;
        for (i, &x) in self.map.iter().enumerate() {
            if i > 0 {
                write!(f, " ")?;
            }
            write!(f, "{x}")?;
        }
        write!(f, "]")
    }
}

impl Permutation {
    /// Creates a permutation from a set of disjoint cycles.
    /// Each cycle should be a Vec<usize> representing indices in the cycle.
    /// The cycles must be disjoint (no shared elements).
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let cycles = vec![vec![0, 1, 2], vec![3, 4]];
    /// let p = Permutation::from_disjoint_cycles(&cycles).unwrap();
    /// assert_eq!(p.map(), &[1, 2, 0, 4, 3]);
    ///
    /// // Error if cycles are not disjoint
    /// let invalid = vec![vec![0, 1], vec![1, 2]];
    /// assert!(Permutation::from_disjoint_cycles(&invalid).is_err());
    /// ```
    pub fn from_disjoint_cycles(cycles: &[Vec<usize>]) -> Result<Self, String> {
        // Find the size of the permutation (maximum index + 1)
        let n = cycles
            .iter()
            .flat_map(|cycle| cycle.iter())
            .max()
            .map(|&max| max + 1)
            .unwrap_or(0);

        // Check for duplicates across cycles
        let mut seen = vec![false; n];
        for cycle in cycles {
            for &idx in cycle {
                if seen[idx] {
                    return Err("Cycles are not disjoint".to_string());
                }
                seen[idx] = true;
            }
        }

        // Create the permutation map
        let mut map = (0..n).collect::<Vec<_>>();
        for cycle in cycles {
            if cycle.len() <= 1 {
                continue;
            }

            // Map each element to the next element in the cycle
            for i in 0..cycle.len() {
                let from = cycle[i];
                let to = cycle[(i + 1) % cycle.len()];
                map[from] = to;
            }
        }

        Ok(Permutation::from_map(map))
    }
    /// Creates a permutation from any set of cycles with specified reading order.
    pub fn from_cycles_ordered(cycles: &[Vec<usize>], order: CycleOrder) -> Self {
        // Find the size of the permutation
        let n = cycles
            .iter()
            .flat_map(|cycle| cycle.iter())
            .max()
            .map(|&max| max + 1)
            .unwrap_or(0);

        // Start with identity permutation
        let mut result = Permutation::id(n);

        // Choose iteration order based on CycleOrder
        let cycle_iter: Box<dyn Iterator<Item = &Vec<usize>>> = match order {
            CycleOrder::LastFirst => Box::new(cycles.iter().rev()),
            CycleOrder::FirstFirst => Box::new(cycles.iter()),
        };

        // Apply cycles in specified order
        for cycle in cycle_iter {
            if cycle.len() <= 1 {
                continue;
            }

            // Create a single cycle permutation
            let mut cycle_map = (0..n).collect::<Vec<_>>();
            for i in 0..cycle.len() {
                let from = cycle[i];
                let to = cycle[(i + 1) % cycle.len()];
                cycle_map[from] = to;
            }
            let cycle_perm = Permutation::from_map(cycle_map);
            result = cycle_perm.compose(&result);
        }

        result
    }

    /// Creates a permutation from cycles using right-to-left reading order (default).
    /// Equivalent to `from_cycles_ordered(cycles, CycleOrder::RightToLeft)`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use linnet::permutation::Permutation;
    /// let cycles = vec![vec![0, 1, 2], vec![1, 2]];
    /// let p = Permutation::from_cycles(&cycles);
    /// assert_eq!(p.map(), &[1, 0, 2]);
    /// ```
    pub fn from_cycles(cycles: &[Vec<usize>]) -> Self {
        Self::from_cycles_ordered(cycles, CycleOrder::LastFirst)
    }

    pub fn length(&self) -> usize {
        self.map.len()
    }

    pub fn generate_all(generators: &[Permutation]) -> Result<Vec<Permutation>, PermutationError> {
        let size = if let Some(generator) = generators.first() {
            generator.length()
        } else {
            return Err(PermutationError::EmptyGenerators);
        };
        let mut all = AHashSet::new();
        let mut stack = vec![Permutation::id(size)];
        for g in generators {
            if g.length() != size {
                return Err(PermutationError::InvalidGeneratorLength);
            }
        }

        while let Some(current) = stack.pop() {
            for g in generators {
                let next = current.compose(g);
                if all.insert(next.clone()) {
                    stack.push(next);
                }
            }
        }

        Ok(all.drain().collect())
    }
}

#[derive(Error, Debug)]
pub enum PermutationError {
    #[error("Invalid generator length")]
    InvalidGeneratorLength,

    #[error("Empty generators")]
    EmptyGenerators,
}

pub trait HedgeGraphExt {
    fn hedges_between(&self, a: NodeIndex, b: NodeIndex) -> SuBitGraph;

    fn permute_subgraph<S: SubSetLike>(&self, subgraph: &S, hedge_perm: &Permutation)
        -> SuBitGraph;

    fn orientation_ord(&self, hedge: Hedge) -> u8;
}

pub trait PermutationExt<Orderer: Ord = ()> {
    // fn from_disjoint_cycles(cycles: &[Vec<usize>]) -> Result<Self, String>;
    //
    type Output;
    type Edges;

    fn permute_vertices(
        &self,
        perm: &Permutation,
        ord: &impl Fn(&Self::Edges) -> Orderer,
    ) -> Vec<Self::Output>;

    fn sort_by(&self, hedge: Hedge, ord: &impl Fn(&Self::Edges) -> Orderer) -> impl Ord;

    fn sort_by_perm(
        &self,
        hedge: Hedge,
        perm: &Permutation,
        ord: &impl Fn(&Self::Edges) -> Orderer,
    ) -> impl Ord;
}

impl<E, V, H> HedgeGraphExt for HedgeGraph<E, V, H> {
    fn hedges_between(&self, a: NodeIndex, b: NodeIndex) -> SuBitGraph {
        let a_ext = InternalSubGraph::cleaned_filter_optimist(
            SuBitGraph::from_hedge_iter(self.iter_crown(a), self.n_hedges()),
            self,
        )
        .filter;
        let b_ext = InternalSubGraph::cleaned_filter_optimist(
            SuBitGraph::from_hedge_iter(self.iter_crown(b), self.n_hedges()),
            self,
        )
        .filter;
        a_ext.intersection(&b_ext)
    }

    fn permute_subgraph<S: SubSetLike>(
        &self,
        subgraph: &S,
        hedge_perm: &Permutation,
    ) -> SuBitGraph {
        let mut permuted_subgraph: SuBitGraph = self.empty_subgraph();

        for h in subgraph.included_iter() {
            permuted_subgraph.add(Hedge(hedge_perm[h.0]));
        }
        permuted_subgraph
    }

    fn orientation_ord(&self, hedge: Hedge) -> u8 {
        match self.superficial_hedge_orientation(hedge) {
            Some(Flow::Sink) => 1,
            Some(Flow::Source) => 2,
            None => 0,
        }
    }
}

impl<E, V, H, O: Ord> PermutationExt<O> for HedgeGraph<E, V, H> {
    type Output = Permutation;
    type Edges = E;

    fn permute_vertices(
        &self,
        perm: &Permutation,
        ord: &impl Fn(&Self::Edges) -> O,
    ) -> Vec<Self::Output> {
        // fn generator_pair(shift: usize, n: usize, total_size: usize) -> [Permutation; 2] {
        //     let mut map = (0..total_size).collect::<Vec<_>>();
        //     for i in 0..n {
        //         map[i + shift] = (i + 1) % n;
        //     }
        //     [
        //         Permutation::from_map(map),
        //         Permutation::from_cycles(&[vec![shift, shift + n - 1], vec![total_size]]),
        //     ]
        // }
        //
        fn generator_pair(n: usize) -> [Permutation; 2] {
            let mut map = Vec::new();
            for i in 0..n {
                map.push((i + 1) % n)
            }
            [
                Permutation::from_map(map),
                Permutation::from_cycles(&[vec![0, n - 1]]),
            ]
        }
        let transpositions = perm.map();
        let mut perms = Vec::new();

        let n = self.n_hedges();
        let mut map = (0..n).collect::<Vec<_>>();

        for (from, &to) in transpositions.iter().enumerate() {
            let from_hairs = self.iter_crown(NodeIndex(from));
            let mut from_hedges = from_hairs.collect::<Vec<_>>();
            let to_hairs = self.iter_crown(NodeIndex(to));
            let mut to_hedges = to_hairs.collect::<Vec<_>>();

            let mut to = BTreeMap::new();

            to_hedges.iter().for_each(|&hedge| {
                to.entry(self.sort_by_perm(hedge, perm, ord))
                    .or_insert_with(Vec::new)
                    .push(hedge);
            });
            to_hedges = to
                .values()
                .flat_map(|values| {
                    if values.len() > 1
                        && self.node_id(self.inv(values[0])).0 <= self.node_id(values[0]).0
                    //insure no double counting of inter-edge permutation
                    {
                        let gen_pair = generator_pair(values.len());

                        let cycle = gen_pair[0].apply_slice(values);
                        let trans = gen_pair[1].apply_slice(values);

                        // println!("{values:?}{cycle:?}{trans:?}");

                        let mut cycle_map = (0..n).collect::<Vec<_>>();
                        let mut trans_map = (0..n).collect::<Vec<_>>();
                        for i in 0..values.len() {
                            cycle_map[cycle[i].0] = values[i].0;
                            trans_map[trans[i].0] = values[i].0;
                        }

                        perms.push(Permutation::from_map(cycle_map));
                        perms.push(Permutation::from_map(trans_map));
                    }
                    values
                })
                .cloned()
                .collect();

            // let mut to = BTreeMap::new();

            // to_hedges.iter().for_each(|&hedge| {
            //     to.entry(self.sort_by(hedge, ord))
            //         .or_insert_with(Vec::new)
            //         .push(hedge);
            // });

            from_hedges.sort_by(|a, b| self.sort_by(*a, ord).cmp(&self.sort_by(*b, ord)));

            // to_hedges.sort_by(|a, b| {
            //     self.sort_by_perm(*a, perm, ord)
            //         .cmp(&self.sort_by_perm(*b, perm, ord))
            // });

            // println!("from:{from_hedges:?}");

            for (from_hedge, to_hedge) in from_hedges.iter().zip(to_hedges.iter()) {
                // println!("{from_hedge}->{to_hedge}");
                map[from_hedge.0] = to_hedge.0;
            }
        }
        let mut maps = vec![];
        let map = Permutation::from_map(map);
        // println!("Map:{map}");
        match Permutation::generate_all(&perms) {
            Ok(a) => {
                for p in a {
                    // println!("Permutation:{p}");
                    // println!("Composition:{}", map.compose(&p));
                    maps.push(map.compose(&p));
                }
            }
            Err(PermutationError::EmptyGenerators) => maps.push(map),
            Err(e) => panic!("Error generating permutations: {e}"),
        }
        maps
    }

    fn sort_by(&self, hedge: Hedge, ord: &impl Fn(&Self::Edges) -> O) -> impl Ord {
        (
            self.involved_node_id(hedge)
                .unwrap_or(self.node_id(hedge))
                .0,
            self.orientation_ord(hedge),
            ord(self.get_edge_data(hedge)),
        )
    }

    fn sort_by_perm(
        &self,
        hedge: Hedge,
        perm: &Permutation,
        ord: &impl Fn(&Self::Edges) -> O,
    ) -> impl Ord {
        (
            perm.inv()[self
                .involved_node_id(hedge)
                .unwrap_or(self.node_id(hedge))
                .0],
            self.orientation_ord(hedge),
            ord(self.get_edge_data(hedge)),
        )
    }
}

#[derive(Debug)]
pub struct PermutationMapIterMut<'a, T: 'a> {
    ptr: NonNull<T>,
    map_indices: &'a [usize],
    current_map_idx: usize,
    len: usize,
    _marker: PhantomData<&'a mut T>,
}

impl<'a, T: 'a> PermutationMapIterMut<'a, T> {
    fn new(slice: &'a mut [T], map_indices: &'a [usize]) -> Self {
        let slice_len = slice.len();
        let perm_len = map_indices.len();
        assert!(
            slice_len >= perm_len,
            "Slice length ({slice_len}) must be at least the permutation length ({perm_len}) for mutable iteration via map.",
        );
        // This assertion also implicitly ensures that all indices in map_indices are valid
        // if the permutation itself is valid (i.e., all its map indices are < perm_len).
        PermutationMapIterMut {
            ptr: NonNull::new(slice.as_mut_ptr()).expect("Slice pointer cannot be null"),
            map_indices,
            current_map_idx: 0,
            len: perm_len,
            _marker: PhantomData,
        }
    }
}

impl<'a, T> Iterator for PermutationMapIterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_map_idx >= self.map_indices.len() {
            None
        } else {
            let perm_idx = self.map_indices[self.current_map_idx];
            self.current_map_idx += 1;
            self.len -= 1;
            // SAFETY:
            // - `ptr` is valid and non-null for the lifetime 'a.
            // - `perm_idx` is an index from the permutation's map.
            // - `PermutationMapIterMut::new` asserts that `slice.len()` is sufficient,
            //   and that all `map_indices` are valid for a slice of length `map_indices.len()`.
            // - Each call to `next` uses a unique index from `map_indices`,
            //   yielding a mutable reference to a distinct element.
            unsafe { Some(&mut *self.ptr.as_ptr().add(perm_idx)) }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}

impl<T> ExactSizeIterator for PermutationMapIterMut<'_, T> {
    fn len(&self) -> usize {
        self.len
    }
}

impl<T> FusedIterator for PermutationMapIterMut<'_, T> {}

/// An iterator over a slice, ordered by a permutation's `map`.
#[derive(Debug, Clone)] // Added Clone
pub struct PermutationMapIter<'a, T: 'a> {
    slice: &'a [T],
    map_indices: &'a [usize],
    current_map_idx: usize,
}

impl<'a, T: 'a> PermutationMapIter<'a, T> {
    fn new(slice: &'a [T], map_indices: &'a [usize]) -> Self {
        // Assertion to ensure slice is long enough for all indices in map_indices.
        // This relies on the Permutation's map_indices being valid for its own length.
        if !map_indices.is_empty() {
            assert!(
                map_indices.iter().all(|&idx| idx < slice.len()),
                "Slice length ({}) is too short for all indices in permutation map. Max index in map: {:?}",
                slice.len(),
                map_indices.iter().max()
            );
        }

        PermutationMapIter {
            slice,
            map_indices,
            current_map_idx: 0,
        }
    }
}

impl<'a, T> Iterator for PermutationMapIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_map_idx >= self.map_indices.len() {
            None
        } else {
            let perm_idx = self.map_indices[self.current_map_idx];
            // println!("{perm_idx}{}", self.current_map_idx);
            self.current_map_idx += 1;
            // Safety: new() asserts that all map_indices are valid for the slice.
            Some(&self.slice[perm_idx])
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.map_indices.len() - self.current_map_idx;
        (len, Some(len))
    }
}

impl<T> ExactSizeIterator for PermutationMapIter<'_, T> {
    fn len(&self) -> usize {
        self.map_indices.len() - self.current_map_idx
    }
}

impl<T> FusedIterator for PermutationMapIter<'_, T> {}

impl Index<usize> for Permutation {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.map()[index]
    }
}

#[cfg(test)]
mod tests {

    use crate::half_edge::{builder::HedgeGraphBuilder, involution::Flow};

    use super::*;

    #[test]
    fn test_from_disjoint_cycles() {
        // Test disjoint cycles
        let cycles = vec![vec![0, 3, 2], vec![1, 4]];
        let p = Permutation::from_disjoint_cycles(&cycles).unwrap();
        assert_eq!(p.map(), &[3, 4, 0, 2, 1]);

        // Test single cycle
        let cycles = vec![vec![0, 1, 2]];
        let p = Permutation::from_disjoint_cycles(&cycles).unwrap();
        assert_eq!(p.map(), &[1, 2, 0]);

        // Test non-disjoint cycles
        let cycles = vec![vec![0, 1], vec![1, 2]];
        assert!(Permutation::from_disjoint_cycles(&cycles).is_err());

        // Test empty cycles
        let cycles: Vec<Vec<usize>> = vec![];
        let p = Permutation::from_disjoint_cycles(&cycles).unwrap();
        assert!(p.map().is_empty());

        // Test single element cycles
        let cycles = vec![vec![0]];
        let p = Permutation::from_disjoint_cycles(&cycles).unwrap();
        assert_eq!(p.map(), &[0]);
    }

    #[test]
    fn test_from_cycles_ordered() {
        let cycles = vec![vec![0, 1, 2], vec![1, 2]];

        // Right to left: (0 1 2)(1 2)
        // First (1 2): [0 2 1]
        // Then (0 1 2): [1 0 2]
        let p = Permutation::from_cycles_ordered(&cycles, CycleOrder::LastFirst);
        assert_eq!(p.map(), &[1, 0, 2]);

        // Left to right: (1 2)(0 1 2)
        // First (0 1 2): [1 2 0]
        // Then (1 2): [0 2 1]
        let p = Permutation::from_cycles_ordered(&cycles, CycleOrder::FirstFirst);
        assert_eq!(p.map(), &[2, 1, 0]);
    }

    #[test]
    fn test_cycle_composition_both_orders() {
        // Test (1 2)(0 1)
        let cycles = vec![vec![1, 2], vec![0, 1]];

        // Right to left: apply (0 1) then (1 2)
        let p = Permutation::from_cycles_ordered(&cycles, CycleOrder::LastFirst);
        assert_eq!(p.map(), &[2, 0, 1]);

        // Left to right: apply (1 2) then (0 1)
        let p = Permutation::from_cycles_ordered(&cycles, CycleOrder::FirstFirst);
        assert_eq!(p.map(), &[1, 2, 0]);
    }

    #[test]
    fn test_default_order() {
        let cycles = vec![vec![0, 1, 2], vec![1, 2]];
        let p1 = Permutation::from_cycles(&cycles);
        let p2 = Permutation::from_cycles_ordered(&cycles, CycleOrder::LastFirst);
        assert_eq!(p1, p2);
    }

    #[test]
    fn test_disjoint_cycles_order_invariant() {
        let cycles = vec![vec![0, 1], vec![2, 3]];
        let p1 = Permutation::from_cycles_ordered(&cycles, CycleOrder::LastFirst);
        let p2 = Permutation::from_cycles_ordered(&cycles, CycleOrder::FirstFirst);
        // Order shouldn't matter for disjoint cycles
        assert_eq!(p1, p2);
    }

    #[test]
    fn test_from_cycles() {
        // Test single cycle
        let cycles = vec![vec![0, 1, 2]];
        let p = Permutation::from_cycles(&cycles);
        assert_eq!(p.map(), &[1, 2, 0]);

        // Test overlapping cycles
        let cycles = vec![vec![0, 1, 2], vec![1, 2]];
        let p = Permutation::from_cycles(&cycles);
        assert_eq!(p.map(), &[1, 0, 2]);

        // Test multiple disjoint cycles
        let cycles = vec![vec![0, 1], vec![2, 3]];
        let p = Permutation::from_cycles(&cycles);
        assert_eq!(p.map(), &[1, 0, 3, 2]);

        // Test composition order
        let cycles = vec![vec![0, 1], vec![1, 2]]; // (0 1)(1 2) = (0 1 2)
        let p = Permutation::from_cycles(&cycles);
        assert_eq!(p.map(), &[1, 2, 0]);
    }

    #[test]
    fn test_cycle_composition_properties() {
        // Test that (0 1 2)(1 2) = (0 1)()
        let p = Permutation::from_cycles(&[vec![0, 1, 2], vec![1, 2]]);
        let q = Permutation::from_cycles(&[vec![1, 0], vec![2]]);
        assert_eq!(p, q);

        // Test that (0 1)(1 2) = (0 1 2)
        let p = Permutation::from_cycles(&[vec![0, 1], vec![1, 2]]);
        let q = Permutation::from_cycles(&[vec![1, 2, 0]]);
        assert_eq!(p, q);
    }

    #[test]
    fn test_myrvold_ruskey_rank1() {
        let p = Permutation::from_map(vec![2, 1, 3, 0]);
        assert_eq!(p.myrvold_ruskey_rank1(), 12);
        for i in 0..=23 {
            assert_eq!(
                i,
                Permutation::myrvold_ruskey_rank1(Permutation::myrvold_ruskey_unrank1(4, i))
            );
        }
    }

    #[test]
    fn test_myrvold_ruskey_rank2() {
        let p = Permutation::myrvold_ruskey_unrank2(4, 1);
        assert_eq!(p.map, vec![2, 1, 3, 0]);
        for i in 0..=23 {
            assert_eq!(
                i,
                Permutation::myrvold_ruskey_rank2(Permutation::myrvold_ruskey_unrank2(4, i))
            );
        }
    }

    #[test]
    fn test_apply_slice() {
        let map = vec![2, 1, 3, 0];
        let pp = Permutation::sort(&map);
        let p = Permutation::from_map(map);
        assert_eq!(p, pp);
        let data = vec![10, 20, 30, 40];
        let permuted = p.apply_slice(&data);
        assert_eq!(permuted, vec![40, 20, 10, 30]);
    }

    #[test]
    fn test_iter_slice_concrete_type_properties() {
        let p = Permutation::from_map(vec![1, 0, 2]);
        let data = vec![10, 20, 30];
        let mut iter = p.iter_slice(&data);

        assert_eq!(iter.len(), 3);
        assert_eq!(iter.size_hint(), (3, Some(3)));
        assert_eq!(iter.next(), Some(&20));
        assert_eq!(iter.len(), 2);
        assert_eq!(iter.size_hint(), (2, Some(2)));

        let iter_clone = iter.clone();
        assert_eq!(iter.collect::<Vec<_>>(), vec![&10, &30]); // Consumes iter
        assert_eq!(iter_clone.collect::<Vec<_>>(), vec![&10, &30]); // iter_clone is independent
    }

    #[test]
    fn test_iter_slice_inv_concrete_type_properties() {
        let p = Permutation::from_map(vec![1, 0, 2]); // inv: [1, 0, 2]
        let data = vec![10, 20, 30];
        let mut iter = p.iter_slice_inv(&data);

        assert_eq!(iter.len(), 3);
        assert_eq!(iter.size_hint(), (3, Some(3)));
        assert_eq!(iter.next(), Some(&20));
        assert_eq!(iter.len(), 2);
        assert_eq!(iter.size_hint(), (2, Some(2)));

        let iter_clone = iter.clone();
        // original iter: after next(), it points to data[p.inv[1]=0] (10), then data[p.inv[2]=2] (30)
        assert_eq!(iter.collect::<Vec<_>>(), vec![&10, &30]);
        assert_eq!(iter_clone.collect::<Vec<_>>(), vec![&10, &30]);
    }

    #[test]
    #[should_panic(expected = "Slice length (2) is too short for all indices in permutation map.")]
    fn test_iter_slice_panic_too_short() {
        let p = Permutation::from_map(vec![0, 1, 2]); // Needs slice of length at least 3
        let data = vec![0, 0]; // Too short
        let _iter = p.iter_slice(&data);
    }

    #[test]
    #[should_panic(expected = "Slice length (3) is too short for all indices in permutation map.")]
    fn test_iter_slice_panic_map_idx_out_of_bounds() {
        // This case should ideally be caught by Permutation::from_map validation
        // but if map somehow contains an index larger than its own length, this test covers iter_slice.
        // However, current Permutation::from_map generates inv based on map.len(),
        // so map itself won't have indices > map.len().
        // This test ensures that if map_indices somehow had an index pointing beyond slice.len(), it would panic.
        let p = Permutation {
            map: vec![0, 5, 1],
            inv: vec![0, 2, 1],
        }; // map[1] = 5
        let data = vec![10, 20, 30]; // slice.len() = 3, but map[1] is 5.
        let _iter = p.iter_slice(&data);
    }

    #[test]
    fn test_apply_slice_inv() {
        let p = Permutation::from_map(vec![2, 1, 3, 0]);
        let data = vec![10, 20, 30, 40];
        let permuted = p.apply_slice_inv(&data);
        assert_eq!(permuted, vec![30, 20, 40, 10]);
    }

    #[test]
    fn test_find_cycles() {
        let p = Permutation::from_map(vec![2, 0, 1, 3]);
        let cycles = p.find_cycles();
        assert_eq!(cycles, vec![vec![0, 2, 1], vec![3]]);
    }

    #[test]
    fn test_cycle_to_transpositions() {
        let cycle = vec![0, 2, 1];
        let transpositions = Permutation::cycle_to_transpositions(&cycle);
        assert_eq!(transpositions, vec![(0, 1), (0, 2)]);
    }

    #[test]
    fn test_transpositions() {
        let p = Permutation::sort([2, 1, 0, 3]);

        let transpositions = p.transpositions();
        println!("{:?}", transpositions.len());

        let p = Permutation::from_map(vec![2, 0, 1, 3]);
        let transpositions = p.transpositions();
        assert_eq!(transpositions, vec![(0, 1), (0, 2)]);
    }

    #[test]
    fn test_apply_slice_in_place() {
        let p = Permutation::from_map(vec![2, 0, 1, 3]);
        let mut data = vec![10, 20, 30, 40];

        let a = p.apply_slice(&data);
        p.apply_slice_in_place(&mut data);
        assert_eq!(data, vec![20, 30, 10, 40]);
        assert_eq!(a, data);
    }

    #[test]
    fn test_apply_slice_in_place_inv() {
        let p = Permutation::from_map(vec![2, 0, 1, 3]);
        let mut data = vec![10, 20, 30, 40];
        let a = p.apply_slice_inv(&data);
        p.apply_slice_in_place_inv(&mut data);
        assert_eq!(data, vec![30, 10, 20, 40]);
        assert_eq!(data, a);
    }

    #[test]
    fn test_sort() {
        let data = vec![30, 10, 20, 40];
        let perm = Permutation::sort(&data);
        let sorted_data = perm.apply_slice(&data);
        assert_eq!(sorted_data, vec![10, 20, 30, 40]);
    }

    #[test]
    fn test_sort_inverse() {
        let data = vec![30, 10, 20, 40];
        let perm = Permutation::sort(&data);
        let sorted_data = perm.apply_slice(&data);
        assert_eq!(sorted_data, vec![10, 20, 30, 40]);

        let inv_perm = perm.inverse();
        let original_data = inv_perm.apply_slice(&sorted_data);
        assert_eq!(original_data, data);
    }

    #[test]
    fn test_iter_slice() {
        let p = Permutation::from_map(vec![2, 0, 1, 3]);
        let data = vec!["a", "b", "c", "d"];
        let mut iter = p.iter_slice(&data);
        assert_eq!(iter.next(), Some(&"c")); // map[0] = 2 -> data[2]
        assert_eq!(iter.next(), Some(&"a")); // map[1] = 0 -> data[0]
        assert_eq!(iter.next(), Some(&"b")); // map[2] = 1 -> data[1]
        assert_eq!(iter.next(), Some(&"d")); // map[3] = 3 -> data[3]
        assert_eq!(iter.next(), None);

        let p_id = Permutation::id(3);
        let data_id = vec![100, 200, 300];
        let collected: Vec<&i32> = p_id.iter_slice(&data_id).collect();
        assert_eq!(collected, vec![&100, &200, &300]);
    }

    #[test]
    fn test_iter_slice_inv() {
        let p = Permutation::from_map(vec![2, 0, 1, 3]); // inv: [1, 2, 0, 3]
        let data = vec!["a", "b", "c", "d"];
        let mut iter = p.iter_slice_inv(&data);
        assert_eq!(iter.next(), Some(&"b")); // inv[0] = 1 -> data[1]
        assert_eq!(iter.next(), Some(&"c")); // inv[1] = 2 -> data[2]
        assert_eq!(iter.next(), Some(&"a")); // inv[2] = 0 -> data[0]
        assert_eq!(iter.next(), Some(&"d")); // inv[3] = 3 -> data[3]
        assert_eq!(iter.next(), None);

        let p_id = Permutation::id(3);
        let data_id = vec![100, 200, 300];
        let collected: Vec<&i32> = p_id.iter_slice_inv(&data_id).collect();
        assert_eq!(collected, vec![&100, &200, &300]);
    }

    #[test]
    fn test_is_identity() {
        let p = Permutation::id(5);
        assert!(p.is_identity());

        let q = Permutation::from_map(vec![1, 0, 2]);
        assert!(!q.is_identity());
    }

    #[test]
    fn test_sign() {
        // Even permutation
        let p = Permutation::from_map(vec![1, 0, 3, 2]);
        assert_eq!(p.sign(), 1);

        // Odd permutation
        let q = Permutation::from_map(vec![2, 1, 0]);
        assert_eq!(q.sign(), -1);
    }

    #[test]
    fn test_pow() {
        let p = Permutation::from_map(vec![1, 2, 0]);
        // p^1 = p
        assert_eq!(p.pow(1), p);

        // p^2 = 0->2,1->0,2->1 => [2,0,1]
        let p2 = p.pow(2);
        assert_eq!(p2.map(), &[2, 0, 1]);

        // p^3 = identity
        let p3 = p.pow(3);
        assert_eq!(p3, Permutation::id(3));
    }

    #[test]
    fn test_compose() {
        // p1: 0→1, 1→2, 2→0 (a cycle of length 3)
        let p1 = Permutation::from_map(vec![1, 2, 0]);

        // p2: 0→2, 1→0, 2→1 (the inverse of p1)
        let p2 = Permutation::from_map(vec![2, 0, 1]);

        // By definition in your code, compose(self, other) = x ↦ other(self(x)).
        //
        // So p2 ∘ p1 means apply p1 first, then p2. Since p2 is the inverse of p1,
        // their composition should be the identity permutation.
        let c1 = p1.compose(&p2);
        assert_eq!(c1, Permutation::id(3), "Expected p2 ∘ p1 = identity");

        // Likewise, p1 ∘ p2 = identity
        let c2 = p2.compose(&p1);
        assert_eq!(c2, Permutation::id(3), "Expected p1 ∘ p2 = identity");

        {
            let p1 = Permutation::from_map(vec![1, 0, 2]); // (0 1)
            let p2 = Permutation::from_map(vec![0, 2, 1]); // (1 2)

            let c1 = p1.compose(&p2); // (0 1)(1 2) = (0 1 2)
            let c2 = p2.compose(&p1); // (1 2)(0 1) = (0 2 1)

            assert_ne!(
                c1, c2,
                "Different order of composition should give different results"
            );
            assert_eq!(c1.map(), &[1, 2, 0]);
            assert_eq!(c2.map(), &[2, 0, 1]);
        }

        // Test associativity: (a ∘ b) ∘ c = a ∘ (b ∘ c)
        {
            let p1 = Permutation::from_map(vec![1, 2, 0]);
            let p2 = Permutation::from_map(vec![2, 0, 1]);
            let p3 = Permutation::from_map(vec![0, 2, 1]);

            let c1 = p1.compose(&p2).compose(&p3);
            let c2 = p1.compose(&p2.compose(&p3));

            assert_eq!(c1, c2, "Composition should be associative");
        }

        // Test composition of disjoint cycles
        {
            let p1 = Permutation::from_map(vec![1, 0, 2, 3]); // (0 1)
            let p2 = Permutation::from_map(vec![0, 1, 3, 2]); // (2 3)

            let c = p1.compose(&p2);
            assert_eq!(
                c.map(),
                &[1, 0, 3, 2],
                "Disjoint cycles should compose independently"
            );
        }

        // Test composition order with multiple elements
        {
            let p1 = Permutation::from_map(vec![1, 2, 3, 0]); // (0 1 2 3)
            let p2 = Permutation::from_map(vec![3, 2, 1, 0]); // (0 3)(1 2)

            let c = p1.compose(&p2);
            assert_eq!(
                c.map(),
                &[0, 3, 2, 1],
                "Complex composition should work correctly"
            );
        }
    }

    #[test]
    fn test_permute_graph() {
        let mut triangle = HedgeGraphBuilder::new();

        let a = triangle.add_node(());
        let b = triangle.add_node(());
        let c = triangle.add_node(());

        triangle.add_edge(a, b, (), false);
        triangle.add_edge(b, c, (), false);
        triangle.add_edge(c, a, (), false);

        let graph: HedgeGraph<(), (), ()> = triangle.build();
        let perm = Permutation::from_cycles(&[vec![0, 2], vec![1]]); //permutes a and c
        let perm2 = Permutation::from_cycles(&[vec![0, 1, 2]]); //permutes a and c

        let a_b_edge = graph.hedges_between(a, b);
        let b_c_edge = graph.hedges_between(b, c);
        let c_a_edge = graph.hedges_between(c, a);

        let hedge_perm = graph.permute_vertices(&perm, &|_| ());
        let hedge_perm2 = graph.permute_vertices(&perm2, &|_| ());
        let permuted_b_c_edge = graph.permute_subgraph(&b_c_edge, &hedge_perm[0]);
        let permuted_b_c_edge2 = graph.permute_subgraph(&b_c_edge, &hedge_perm2[0]);
        let permuted_a_b_edge = graph.permute_subgraph(&a_b_edge, &hedge_perm[0]);
        let permuted_a_b_edge2 = graph.permute_subgraph(&a_b_edge, &hedge_perm2[0]);
        let permuted_c_a_edge = graph.permute_subgraph(&c_a_edge, &hedge_perm[0]);
        let permuted_c_a_edge2 = graph.permute_subgraph(&c_a_edge, &hedge_perm2[0]);

        assert_eq!(hedge_perm.len(), 1);

        println!(
            "//a-c\n{}\n//a-b\n{}\n//b-c\n{}",
            graph.dot(&c_a_edge),
            graph.dot(&a_b_edge),
            graph.dot(&b_c_edge)
        );
        println!(
            "//permuted a-b\n{}\n//permuted b-c\n{}\n//permuted c-a\n{}",
            graph.dot(&permuted_a_b_edge),
            graph.dot(&permuted_b_c_edge),
            graph.dot(&permuted_c_a_edge)
        );
        println!(
            "//permuted a-b\n{}\n//permuted b-c\n{}\n//permuted c-a\n{}",
            graph.dot(&permuted_a_b_edge2),
            graph.dot(&permuted_b_c_edge2),
            graph.dot(&permuted_c_a_edge2)
        );

        assert_eq!(a_b_edge, permuted_b_c_edge);
        assert_eq!(b_c_edge, permuted_a_b_edge);
        assert_eq!(c_a_edge, permuted_c_a_edge);

        // println!(
        //     "{}\n{}",
        //     graph.dot(&cleaned_subgraph),
        //     graph.dot(&permuted_subgraph)
        // )
        // let permuted_subgraph = hedgeperm.apply_slice(cleaned_subgraph.filter);
    }

    #[test]
    fn test_permute_graph_double_edges() {
        let mut triangle = HedgeGraphBuilder::new();

        let a = triangle.add_node(());
        let b = triangle.add_node(());
        let c = triangle.add_node(());

        triangle.add_edge(a, b, (), false);

        triangle.add_external_edge(a, (), false, Flow::Sink);
        triangle.add_edge(b, c, (), false);
        triangle.add_external_edge(c, (), false, Flow::Source);
        triangle.add_edge(c, a, (), false);
        triangle.add_edge(c, a, (), false);
        let graph: HedgeGraph<(), (), ()> = triangle.build();
        let perm = Permutation::from_cycles(&[vec![0, 2], vec![1]]); //permutes a and c

        let a_b_edge = graph.hedges_between(a, b);
        let b_c_edge = graph.hedges_between(b, c);
        let c_a_edge = graph.hedges_between(c, a);

        let hedge_perm = graph.permute_vertices(&perm, &|_| ());
        let permuted_b_c_edge = graph.permute_subgraph(&b_c_edge, &hedge_perm[0]);
        let permuted_a_b_edge = graph.permute_subgraph(&a_b_edge, &hedge_perm[0]);
        let permuted_c_a_edge = graph.permute_subgraph(&c_a_edge, &hedge_perm[0]);

        assert_eq!(hedge_perm.len(), 2);
        assert_eq!(a_b_edge, permuted_b_c_edge);
        assert_eq!(b_c_edge, permuted_a_b_edge);
        assert_eq!(c_a_edge, permuted_c_a_edge);

        println!(
            "//a-c\n{}\n//a-b\n{}\n//b-c\n{}",
            graph.dot(&c_a_edge),
            graph.dot(&a_b_edge),
            graph.dot(&b_c_edge)
        );
        println!(
            "//permuted a-b\n{}\n//permuted b-c\n{}\n//permuted c-a\n{}",
            graph.dot(&permuted_a_b_edge),
            graph.dot(&permuted_b_c_edge),
            graph.dot(&permuted_c_a_edge)
        );
        // let permuted_subgraph = hedgeperm.apply_slice(a.filter);
    }

    #[test]
    fn test_permute_double_triangle() {
        let mut doubletriangle = HedgeGraphBuilder::new();

        let a = doubletriangle.add_node(());
        let b = doubletriangle.add_node(());
        let c = doubletriangle.add_node(());
        let d = doubletriangle.add_node(());

        doubletriangle.add_edge(a, b, (), false);

        doubletriangle.add_edge(b, c, (), false);
        doubletriangle.add_edge(c, d, (), false);
        doubletriangle.add_edge(d, a, (), false);
        doubletriangle.add_edge(c, a, (), false);
        let graph: HedgeGraph<(), (), ()> = doubletriangle.build();
        let perm = Permutation::from_cycles(&[vec![0, 2], vec![3]]); //permutes a and c

        let a_b_edge = graph.hedges_between(a, b);
        let b_c_edge = graph.hedges_between(b, c);
        let c_a_edge = graph.hedges_between(c, a);

        let hedge_perm = graph.permute_vertices(&perm, &|_| ());
        let permuted_b_c_edge = graph.permute_subgraph(&b_c_edge, &hedge_perm[0]);
        let permuted_a_b_edge = graph.permute_subgraph(&a_b_edge, &hedge_perm[0]);
        let permuted_c_a_edge = graph.permute_subgraph(&c_a_edge, &hedge_perm[0]);

        assert_eq!(hedge_perm.len(), 1);
        assert_eq!(a_b_edge, permuted_b_c_edge);
        assert_eq!(b_c_edge, permuted_a_b_edge);
        assert_eq!(c_a_edge, permuted_c_a_edge);

        println!(
            "//a-c\n{}\n//a-b\n{}\n//b-c\n{}",
            graph.dot(&c_a_edge),
            graph.dot(&a_b_edge),
            graph.dot(&b_c_edge)
        );
        println!(
            "//permuted a-b\n{}\n//permuted b-c\n{}\n//permuted c-a\n{}",
            graph.dot(&permuted_a_b_edge),
            graph.dot(&permuted_b_c_edge),
            graph.dot(&permuted_c_a_edge)
        );
        // let permuted_subgraph = hedgeperm.apply_slice(a.filter);
    }

    #[test]
    fn cyclic_doubled_triangle() {
        let mut doubledtriangle = HedgeGraphBuilder::new();

        let a = doubledtriangle.add_node(());
        let b = doubledtriangle.add_node(());
        let c = doubledtriangle.add_node(());

        doubledtriangle.add_edge(a, b, (), false);
        doubledtriangle.add_edge(a, b, (), false);
        doubledtriangle.add_edge(b, c, (), false);
        doubledtriangle.add_edge(b, c, (), false);
        doubledtriangle.add_edge(c, a, (), false);
        doubledtriangle.add_edge(c, a, (), false);
        let graph: HedgeGraph<(), (), ()> = doubledtriangle.build();
        let perm = Permutation::from_cycles(&[vec![0, 1, 2]]); //single cycle

        let a_b_edge = graph.hedges_between(a, b);
        let b_c_edge = graph.hedges_between(b, c);
        let c_a_edge = graph.hedges_between(c, a);

        let hedge_perm = graph.permute_vertices(&perm, &|_| ());
        let permuted_b_c_edge = graph.permute_subgraph(&b_c_edge, &hedge_perm[0]);
        let permuted_a_b_edge = graph.permute_subgraph(&a_b_edge, &hedge_perm[0]);
        let permuted_c_a_edge = graph.permute_subgraph(&c_a_edge, &hedge_perm[0]);

        assert_eq!(hedge_perm.len(), 8); //2^3

        println!(
            "//a-c\n{}\n//a-b\n{}\n//b-c\n{}",
            graph.dot(&c_a_edge),
            graph.dot(&a_b_edge),
            graph.dot(&b_c_edge)
        );
        println!(
            "//permuted a-b\n{}\n//permuted b-c\n{}\n//permuted c-a\n{}",
            graph.dot(&permuted_a_b_edge),
            graph.dot(&permuted_b_c_edge),
            graph.dot(&permuted_c_a_edge)
        );
        assert_eq!(a_b_edge, permuted_c_a_edge);
        assert_eq!(c_a_edge, permuted_b_c_edge);
        assert_eq!(b_c_edge, permuted_a_b_edge);
        // let permuted_subgraph = hedgeperm.apply_slice(a.filter);
    }

    #[test]
    fn cycle_permutation() {
        let mut cycle = HedgeGraphBuilder::new();

        let a = cycle.add_node(());
        let b = cycle.add_node(());

        cycle.add_edge(a, b, (), false);
        cycle.add_edge(b, a, (), false);

        let graph: HedgeGraph<(), (), ()> = cycle.build();
        let perm = Permutation::from_cycles(&[vec![0, 1]]); //permutes a and b

        let h = Hedge(0);
        let mut h_sub: SuBitGraph = graph.empty_subgraph();
        h_sub.add(h);

        let hedge_perm = graph.permute_vertices(&perm, &|_| ());
        let permuted_h = graph.permute_subgraph(&h_sub, &hedge_perm[0]);
        println!("n perms{}", hedge_perm.len());

        println!("//origninal:\n{}", graph.dot(&h_sub));
        println!("//permuted:\n{}", graph.dot(&permuted_h));
        // let permuted_subgraph = hedgeperm.apply_slice(a.filter);
    }

    #[test]
    fn triple_loop_two_nodes() {
        let mut cycle = HedgeGraphBuilder::new();

        let a = cycle.add_node(());
        let b = cycle.add_node(());

        cycle.add_edge(a, b, (), false);
        cycle.add_edge(b, a, (), false);
        cycle.add_edge(b, a, (), false);

        let graph: HedgeGraph<(), (), ()> = cycle.build();
        let perm = Permutation::from_cycles(&[vec![0, 1]]); //permutes a and b

        let h = Hedge(0);
        let mut h_sub: SuBitGraph = graph.empty_subgraph();
        h_sub.add(h);

        let hedge_perm = graph.permute_vertices(&perm, &|_| ());
        let permuted_h = graph.permute_subgraph(&h_sub, &hedge_perm[0]);
        println!("n perms{}", hedge_perm.len());

        println!("//origninal:\n{}", graph.dot(&h_sub));
        println!("//permuted:\n{}", graph.dot(&permuted_h));
        // let permuted_subgraph = hedgeperm.apply_slice(a.filter);
    }

    #[test]
    fn generate_all() {
        fn generator_pair(n: usize) -> [Permutation; 2] {
            let mut map = Vec::new();
            for i in 0..n {
                map.push((i + 1) % n)
            }
            [
                Permutation::from_map(map),
                Permutation::from_cycles(&[vec![0, n - 1]]),
            ]
        }

        let all_2 = Permutation::generate_all(&generator_pair(2)).unwrap();

        let all_3 = Permutation::generate_all(&generator_pair(3)).unwrap();
        let all_4 = Permutation::generate_all(&generator_pair(4)).unwrap();
        assert_eq!(all_2.len(), 2);
        assert_eq!(all_3.len(), 6);
        assert_eq!(all_4.len(), 24);
    }

    #[test]
    fn test_apply_in_place_same_as_ref() {
        // Create a few different permutations for testing
        let permutations = [
            Permutation::from_map(vec![2, 0, 1, 3]),
            Permutation::from_map(vec![3, 2, 1, 0]),
            Permutation::from_map(vec![1, 0, 3, 2]),
        ];

        for perm in &permutations {
            // Test data
            let data = vec![10, 20, 30, 40];

            // Apply using by-reference method
            let result_by_ref = perm.apply_slice(&data);

            // Apply using in-place method
            let mut data_copy = data.clone();
            perm.apply_slice_in_place(&mut data_copy);

            // Results should be the same
            assert_eq!(result_by_ref, data_copy);

            // Also test inverse operations
            let result_inv_by_ref = perm.apply_slice_inv(&data);

            let mut data_copy = data.clone();
            perm.apply_slice_in_place_inv(&mut data_copy);

            // Inverse results should be the same
            assert_eq!(result_inv_by_ref, data_copy);
        }
    }

    #[test]
    fn test_apply_in_place_inverse_relation() {
        // Create a permutation
        let perm = Permutation::from_map(vec![2, 0, 3, 1]);

        // Apply the permutation and then its inverse should give back the original data
        let data = vec![10, 20, 30, 40];

        // Using by-reference methods
        let permuted = perm.apply_slice(&data);
        let unpermuted = perm.apply_slice_inv(&permuted);
        assert_eq!(data, unpermuted);

        // Using in-place methods
        let mut data_copy = data.clone();
        perm.apply_slice_in_place(&mut data_copy);
        perm.apply_slice_in_place_inv(&mut data_copy);
        assert_eq!(data, data_copy);
    }

    #[test]
    fn test_validate_with_sort() {
        // Create data that's out of order
        let data = vec![30, 10, 40, 20];

        // Get sorting permutation
        let perm = Permutation::sort(&data);

        // Sort using by-reference method
        let sorted_by_ref = perm.apply_slice(&data);
        assert_eq!(sorted_by_ref, vec![10, 20, 30, 40]);

        // Sort using in-place method
        let mut data_copy = data.clone();
        perm.apply_slice_in_place(&mut data_copy);
        assert_eq!(data_copy, vec![10, 20, 30, 40]);

        // Verify that unsort works too
        let mut sorted = vec![10, 20, 30, 40];
        perm.apply_slice_in_place_inv(&mut sorted);
        assert_eq!(sorted, data);

        let unsorted_by_ref = perm.apply_slice_inv([10, 20, 30, 40]);
        assert_eq!(unsorted_by_ref, data);
    }

    #[test]
    fn test_apply_equivalence_random_data() {
        // Test with different sized permutations
        for size in 1..10 {
            // Create a random permutation
            let map: Vec<usize> = (0..size).collect();
            let mut perm_map = map.clone();

            // Simple shuffle algorithm for testing
            for i in 0..size {
                let j = (i * 7 + 3) % size;
                perm_map.swap(i, j);
            }

            let perm = Permutation::from_map(perm_map);

            // Create test data
            let test_data: Vec<i32> = (0..size as i32).collect();

            // Apply using by-reference methods
            let result_by_ref = perm.apply_slice(&test_data);

            // Apply using in-place methods
            let mut data_copy = test_data.clone();
            perm.apply_slice_in_place(&mut data_copy);

            // Results should match
            assert_eq!(result_by_ref, data_copy);

            // Test inverse operations
            let result_inv_by_ref = perm.apply_slice_inv(&test_data);

            let mut data_copy = test_data.clone();
            perm.apply_slice_in_place_inv(&mut data_copy);

            assert_eq!(result_inv_by_ref, data_copy);
        }
    }

    #[test]
    fn test_from_mappings() {
        // Test the example from the docstring
        let mappings = [(5, 0), (0, 1)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 6).unwrap();
        assert_eq!(p.map(), &[1, 5, 2, 3, 4, 0]);
        assert_eq!(p[5], 0);
        assert_eq!(p[0], 1);

        // Test simple swap
        let mappings = [(0, 1), (1, 0)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 3).unwrap();
        assert_eq!(p.map(), &[1, 0, 2]);

        // Test single mapping
        let mappings = [(0, 2)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 3).unwrap();
        // Should have 0->2, and need to place 1 and 2 somewhere
        // 2 is used by 0, so 2 can't stay at position 2
        // Available targets are 0 and 1
        // Position 1 is not forced, so it can map to available target
        // Position 2 is not forced, so it can map to available target
        assert_eq!(p[0], 2);
        // The exact arrangement of 1 and 2 may vary, but should be valid
        let mut targets: Vec<_> = p.map().to_vec();
        targets.sort();
        assert_eq!(targets, vec![0, 1, 2]);

        // Test cycle
        let mappings = [(0, 1), (1, 2), (2, 0)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 3).unwrap();
        assert_eq!(p.map(), &[1, 2, 0]);

        // Test conflicting mappings - should fail
        let mappings = [(0, 1), (1, 1)];
        let result = Permutation::from_mappings(mappings.iter().copied(), 2);
        assert!(result.is_err());

        // Test out of bounds - should fail
        let mappings = [(0, 5)];
        let result = Permutation::from_mappings(mappings.iter().copied(), 3);
        assert!(result.is_err());

        // Test empty mappings - should be identity
        let mappings: Vec<(usize, usize)> = vec![];
        let p = Permutation::from_mappings(mappings.iter().copied(), 4).unwrap();
        assert_eq!(p, Permutation::id(4));

        // Test size 0
        let mappings: Vec<(usize, usize)> = vec![];
        let _p = Permutation::from_mappings(mappings.iter().copied(), 0).unwrap();
        // assert_eq!(p.map(), &[]);
    }

    #[test]
    fn test_from_mappings_complex() {
        // Test more complex displacement
        let mappings = [(0, 3), (1, 0)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 4).unwrap();
        assert_eq!(p[0], 3);
        assert_eq!(p[1], 0);
        // Positions 2 and 3 need to be mapped somewhere
        // Target 3 is taken by position 0
        // Target 0 is taken by position 1
        // Available targets are 1 and 2
        // So positions 2 and 3 should map to 1 and 2 in some order
        let mut targets: Vec<_> = p.map().to_vec();
        targets.sort();
        assert_eq!(targets, vec![0, 1, 2, 3]);

        // Test that the permutation is valid
        assert!(p.map().iter().all(|&target| target < 4));
        let mut seen = [false; 4];
        for &target in p.map() {
            assert!(!seen[target]);
            seen[target] = true;
        }
    }

    #[test]
    fn test_from_mappings_user_example() {
        // Test the specific example from the user's request:
        // Going from invalid map [1,1,2,3,4,0] to valid [1,5,2,3,4,0]
        // This means we want: 5->0 and 0->1
        let mappings = [(5, 0), (0, 1)];
        let p = Permutation::from_mappings(mappings.iter().copied(), 6).unwrap();

        // Should result in [1, 5, 2, 3, 4, 0]
        assert_eq!(p.map(), &[1, 5, 2, 3, 4, 0]);

        // Verify the specific mappings
        assert_eq!(p[5], 0); // 5 -> 0 ✓
        assert_eq!(p[0], 1); // 0 -> 1 ✓

        // Verify other elements stay in place where possible
        assert_eq!(p[2], 2); // 2 -> 2 (identity preserved)
        assert_eq!(p[3], 3); // 3 -> 3 (identity preserved)
        assert_eq!(p[4], 4); // 4 -> 4 (identity preserved)

        // Element 1 had to be displaced because its identity target (1)
        // was taken by element 0, so it maps to the available target 5
        assert_eq!(p[1], 5); // 1 -> 5 (displaced to available target)
    }
}

#[cfg(test)]
mod swaps_optimality_tests {
    use super::*;

    /// Return the minimum number of swaps to realize permutation p from identity:
    /// n - number_of_cycles.
    fn min_swaps_of_perm(p: &Permutation) -> usize {
        let n = p.length();
        let cycles = p.find_cycles();
        let num_cycles = cycles.len();
        n - num_cycles
    }

    /// Check whether permutation p satisfies a set of forced mappings.
    fn satisfies_forced(p: &Permutation, forced: &[(usize, usize)]) -> bool {
        forced.iter().all(|&(f, t)| p[f] == t)
    }

    /// Simple next_permutation over a slice of usize in-place; returns false when wrapped.
    fn next_permutation(a: &mut [usize]) -> bool {
        // Find longest non-increasing suffix
        if a.len() < 2 {
            return false;
        }
        let mut i = a.len() - 2;
        while a[i] >= a[i + 1] {
            if i == 0 {
                return false;
            }
            i -= 1;
        }
        // Find rightmost successor to pivot
        let mut j = a.len() - 1;
        while a[j] <= a[i] {
            j -= 1;
        }
        a.swap(i, j);
        // Reverse suffix
        a[i + 1..].reverse();
        true
    }

    /// Exhaustively compute the minimum possible swap-count over all permutations
    /// of size `n` that satisfy `forced`. (Brute force = (n!) scan; only use for small n.)
    fn brute_min_swaps_given_forced(
        n: usize,
        forced: &[(usize, usize)],
    ) -> Option<(usize, Vec<usize>)> {
        if n == 0 {
            return Some((0, vec![]));
        }
        let mut map: Vec<usize> = (0..n).collect(); // identity
        let mut best: Option<(usize, Vec<usize>)> = None;

        loop {
            // candidate permutation
            let cand = Permutation::from_map(map.clone());
            if satisfies_forced(&cand, forced) {
                let cost = min_swaps_of_perm(&cand);
                match &best {
                    None => best = Some((cost, map.clone())),
                    Some((best_cost, _)) if cost < *best_cost => best = Some((cost, map.clone())),
                    _ => {}
                }
            }
            if !next_permutation(&mut map) {
                break;
            }
        }
        best
    }

    #[test]
    fn from_mappings_is_swaps_optimal_on_counterexample() {
        // Classic counterexample where naive zip can make a single 4-cycle (3 swaps)
        // but the optimal completion yields two 2-cycles (2 swaps).
        // Forced: 0->1, 2->3 on n=4.
        let n = 4;
        let forced = vec![(0, 1), (2, 3)];

        let p = Permutation::from_mappings(forced.iter().copied(), n)
            .expect("from_mappings should succeed");

        // Our constructor must achieve the brute-force minimal swap count.
        let (min_cost, _witness) =
            brute_min_swaps_given_forced(n, &forced).expect("At least one completion must exist");
        let cost = min_swaps_of_perm(&p);
        assert_eq!(
            cost, min_cost,
            "from_mappings is not swaps-optimal on the counterexample"
        );

        // Strengthen: for this particular instance, the true optimum is 2.
        assert_eq!(
            min_cost, 2,
            "Unexpected brute-force minimum for the counterexample"
        );
    }

    #[test]
    fn from_mappings_handles_long_forced_chains_optimally() {
        // Forced chain that displaces a block:
        // 0->1, 1->2, 3->4 with n=6.
        // The optimal wiring should maximize the number of cycles among the remaining
        // (displaced) sources.
        let n = 6;
        let forced = vec![(0, 1), (1, 2), (3, 4)];
        let p = Permutation::from_mappings(forced.iter().copied(), n).unwrap();

        let (min_cost, _witness) = brute_min_swaps_given_forced(n, &forced).unwrap();
        let cost = min_swaps_of_perm(&p);
        assert_eq!(
            cost,
            min_cost,
            "Non-optimal swaps for long forced chain: got {cost}, want {min_cost}. p={:?}",
            p.map()
        );

        // Sanity: still respects all forced pairs
        for (f, t) in forced {
            assert_eq!(p[f], t);
        }
    }

    #[test]
    fn from_mappings_preserves_zero_swaps_when_possible() {
        // No forced pairs: must return identity (0 swaps).
        let n = 7;
        let forced: Vec<(usize, usize)> = vec![];
        let p = Permutation::from_mappings(forced.iter().copied(), n).unwrap();
        assert!(p.is_identity());
        assert_eq!(min_swaps_of_perm(&p), 0);

        // Forced identity only: still 0 swaps.
        let forced = [(0, 0), (3, 3), (6, 6)];
        let p = Permutation::from_mappings(forced.iter().copied(), n).unwrap();
        assert_eq!(min_swaps_of_perm(&p), 0);
    }

    #[test]
    fn from_mappings_extends_forced_2cycle_optimally() {
        // Forced 2-cycle plus one extra forced edge that displaces someone else.
        // n=5, forced: (0,1), (1,0) and (3,2).
        // Optimal completion pairs (2,3) into a 2-cycle and leaves 4 fixed, for min swaps.
        let n = 5;
        let forced = vec![(0, 1), (1, 0), (3, 2)];

        let p = Permutation::from_mappings(forced.iter().copied(), n).unwrap();
        let (min_cost, _witness) = brute_min_swaps_given_forced(n, &forced).unwrap();
        let cost = min_swaps_of_perm(&p);
        assert_eq!(
            cost, min_cost,
            "Should close remaining component(s) into max cycles"
        );
        for (f, t) in forced {
            assert_eq!(p[f], t);
        }
    }

    #[test]
    fn from_mappings_user_case_is_swaps_optimal() {
        // The user’s earlier scenario: n=6, forced (5->0),(0->1).
        let n = 6;
        let forced = vec![(5, 0), (0, 1)];

        let p = Permutation::from_mappings(forced.iter().copied(), n).unwrap();
        assert_eq!(p[5], 0);
        assert_eq!(p[0], 1);

        let (min_cost, _witness) = brute_min_swaps_given_forced(n, &forced).unwrap();
        let cost = min_swaps_of_perm(&p);
        assert_eq!(cost, min_cost, "User case should be swaps-optimal too");
    }

    #[test]
    fn from_mappings_exhaustive_small_n() {
        // Fully adversarial: for n=5 we scan a *suite* of forced sets and
        // verify the constructor attains the brute-force minimum swaps each time.
        // (We don't iterate over *all* subsets to keep test time reasonable.)
        let n = 5;

        let forced_sets: Vec<Vec<(usize, usize)>> = vec![
            vec![],                       // empty
            vec![(0, 1)],                 // single edge
            vec![(0, 1), (2, 3)],         // two disjoint forced edges
            vec![(0, 1), (1, 2)],         // small forced chain
            vec![(0, 1), (1, 0)],         // forced 2-cycle
            vec![(0, 2), (2, 4), (4, 0)], // forced 3-cycle
            vec![(0, 2), (2, 4), (1, 3)], // 3-chain + extra edge
            vec![(4, 0)],                 // target at opposite end
        ];

        for forced in forced_sets {
            let p = Permutation::from_mappings(forced.iter().copied(), n)
                .expect("from_mappings should succeed on these cases");

            // Must satisfy constraints
            assert!(
                satisfies_forced(&p, &forced),
                "Does not satisfy forced pairs: {forced:?}",
            );

            // Must achieve global min swaps among all completions
            let (min_cost, _witness) = brute_min_swaps_given_forced(n, &forced)
                .expect("No valid completion in brute force (unexpected)");
            let cost = min_swaps_of_perm(&p);
            assert_eq!(
                cost,
                min_cost,
                "Not swaps-optimal for forced={forced:?}. Got {cost}, want {min_cost}. p={:?}",
                p.map()
            );
        }
    }

    #[test]
    fn from_mappings_rejects_conflicts() {
        // Conflicting forced targets for same to-value
        let n = 4;
        let forced = [(0, 1), (2, 1)];
        assert!(Permutation::from_mappings(forced.iter().copied(), n).is_err());

        // Out of bounds
        let forced = [(0, 5)];
        assert!(Permutation::from_mappings(forced.iter().copied(), n).is_err());
    }
}
