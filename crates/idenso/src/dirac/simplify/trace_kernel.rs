//! Dimension-generic traces reuse factored pairing polynomials within one call.
//! Short four-dimensional traces use more compact integer recipes generated
//! from the same three-gamma identity as open chains. No user atoms or dummy
//! indices escape into the global recipe cache.

use std::{
    collections::{BTreeMap, HashMap},
    sync::{LazyLock, OnceLock},
};

use itertools::Itertools;
use spenso::g;
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{THREE_GAMMA_METRIC_TERMS, is_four_dimension, is_minkowski_slot};
use crate::epsilon::epsilon4;

mod factored;
mod sparse;

#[derive(Clone, Copy)]
pub(super) enum TraceOutput {
    Expanded,
    Factored,
}

pub(super) struct TraceKernel<const N: usize> {
    // For axial traces, the first four entries belong to epsilon; the rest
    // are sorted metric pairs. Ordinary traces contain metric pairs only.
    terms: Vec<([u8; N], i32)>,
    factored: OnceLock<FactoredKernel>,
}

struct FactoredKernel {
    polynomial: factored::FactoredTrace,
    epsilons: Vec<[u8; 4]>,
}

#[derive(Clone, Default)]
struct Monomial {
    metrics: Vec<[u8; 2]>,
    epsilons: Vec<[u8; 4]>,
}

impl<const N: usize> TraceKernel<N> {
    fn generate(axial: bool) -> Self {
        let mut terms = BTreeMap::new();
        Self::reduce_triples(0, 1, 1, Monomial::default(), axial, &mut terms);
        Self {
            terms: terms.into_iter().filter(|(_, c)| *c != 0).collect(),
            factored: OnceLock::new(),
        }
    }

    fn reduce_triples(
        current: u8,
        position: u8,
        sign: i32,
        monomial: Monomial,
        axial: bool,
        terms: &mut BTreeMap<[u8; N], i32>,
    ) {
        if usize::from(position) == N - 1 {
            if monomial.epsilons.len() % 2 == usize::from(axial) {
                let mut monomial = monomial;
                monomial.metrics.push([current, position]);
                // An axial odd prefix is gamma(current) gamma5. Moving its
                // gamma5 past the final gamma adds the remaining minus sign.
                monomial.collect::<N>(if axial { -sign } else { sign }, terms);
            }
            return;
        }
        let triple = [current, position, position + 1];
        for (coefficient, [a, b], remaining) in THREE_GAMMA_METRIC_TERMS {
            let mut next = monomial.clone();
            next.metrics.push([triple[a], triple[b]]);
            Self::reduce_triples(
                triple[remaining],
                position + 2,
                sign * coefficient,
                next,
                axial,
                terms,
            );
        }
        let dummy = N as u8 + position / 2;
        let mut next = monomial;
        next.epsilons.push([triple[0], triple[1], triple[2], dummy]);
        Self::reduce_triples(dummy, position + 2, -sign, next, axial, terms);
    }

    fn factored(&self, axial: bool) -> &FactoredKernel {
        self.factored.get_or_init(|| {
            let mut epsilons = Vec::new();
            let terms = self
                .terms
                .iter()
                .map(|(recipe, coefficient)| {
                    let mut factors = Vec::with_capacity(N / 2);
                    if axial {
                        let key: [u8; 4] = recipe[..4].try_into().unwrap();
                        // Recipes are sorted, so each epsilon's terms are a
                        // contiguous group. Cache only its integer positions.
                        if epsilons.last() != Some(&key) {
                            epsilons.push(key);
                        }
                        factors.push(N * N + epsilons.len() - 1);
                    }
                    factors.extend(
                        recipe[if axial { 4 } else { 0 }..]
                            .as_chunks::<2>()
                            .0
                            .iter()
                            .map(|pair| usize::from(pair[0]) * N + usize::from(pair[1])),
                    );
                    (factors, *coefficient)
                })
                .collect();
            FactoredKernel {
                polynomial: factored::FactoredTrace::new(terms, N * N + epsilons.len()),
                epsilons,
            }
        })
    }

    fn evaluate(&self, indices: &[AtomView<'_>], axial: bool, output: TraceOutput) -> Atom {
        debug_assert_eq!(indices.len(), N);
        // Each scalar product is built once and shared by all template terms.
        let mut metrics: Vec<_> = (0..N)
            .flat_map(|a| (0..N).map(move |b| g!(indices[a], indices[b])))
            .collect();
        if matches!(output, TraceOutput::Factored) && N >= 10 {
            let factored = self.factored(axial);
            metrics.extend(factored.epsilons.iter().map(|&[a, b, c, d]| {
                epsilon4(
                    indices[a as usize],
                    indices[b as usize],
                    indices[c as usize],
                    indices[d as usize],
                )
            }));
            return factored.polynomial.evaluate(
                factored.polynomial.root,
                &metrics,
                Atom::num(4).as_view(),
            );
        }
        let max_coefficient = self.terms.iter().map(|(_, c)| c.abs()).max().unwrap_or(0);
        let coefficients: Vec<_> = (-max_coefficient..=max_coefficient)
            .map(|coefficient| Atom::num(4 * i64::from(coefficient)))
            .collect();
        let mut last_epsilon_key = [u8::MAX; 4];
        let mut last_epsilon = Atom::Zero;
        Atom::add_many(self.terms.iter().map(|(indices_recipe, coefficient)| {
            let epsilon = if axial {
                let key: [u8; 4] = indices_recipe[..4].try_into().unwrap();
                if key != last_epsilon_key {
                    last_epsilon_key = key;
                    let [a, b, c, d] = key;
                    last_epsilon = epsilon4(
                        indices[a as usize],
                        indices[b as usize],
                        indices[c as usize],
                        indices[d as usize],
                    );
                }
                Some(last_epsilon.as_view())
            } else {
                None
            };
            let pairs = &indices_recipe[if axial { 4 } else { 0 }..];
            Atom::mul_many(
                std::iter::once(coefficients[(*coefficient + max_coefficient) as usize].as_view())
                    .chain(epsilon)
                    .chain(
                        pairs
                            .as_chunks::<2>()
                            .0
                            .iter()
                            .map(|p| metrics[usize::from(p[0]) * N + usize::from(p[1])].as_view()),
                    ),
            )
        }))
    }
}

/// Build the ordinary pairing polynomial without intermediate trace nodes.
/// Repeated subwords share their result for this evaluation only.
struct PairingTrace<'a, Output: TraceAlgebra = FactoredOutput> {
    width: usize,
    arguments: Vec<AtomView<'a>>,
    metrics: Vec<Option<Atom>>,
    trace_unit: AtomView<'a>,
    subwords: HashMap<Vec<usize>, Output::Value>,
    output: Output,
    canonical_arguments: bool,
    summed_indices: Vec<usize>,
}

/// Scalar output operations shared by both representations. The word evaluator
/// below owns contraction order and identities, independently of Atom assembly.
trait TraceAlgebra {
    type Value: Clone;
    fn unit(&mut self, trace_unit: AtomView<'_>) -> Self::Value;
    fn scale(&mut self, coefficient: i32, value: Self::Value) -> Self::Value;
    fn metric_product(
        &mut self,
        metric: AtomView<'_>,
        coefficient: i32,
        value: Self::Value,
    ) -> Self::Value;
    fn dimension_product(
        &mut self,
        dimension: AtomView<'_>,
        constant: i32,
        linear: i32,
        value: Self::Value,
    ) -> Self::Value;
    fn add(&mut self, left: Self::Value, right: Self::Value) -> Self::Value;
    fn subtract(&mut self, left: Self::Value, right: Self::Value) -> Self::Value;
    fn sum(&mut self, values: Vec<Self::Value>) -> Self::Value;
    fn finish(self, value: Self::Value) -> Atom;
}

struct FactoredOutput;

impl TraceAlgebra for FactoredOutput {
    type Value = Atom;
    fn unit(&mut self, trace_unit: AtomView<'_>) -> Atom {
        trace_unit.to_owned()
    }
    fn scale(&mut self, coefficient: i32, value: Atom) -> Atom {
        match coefficient {
            1 => value,
            -1 => -value,
            _ => Atom::num(coefficient) * value,
        }
    }
    fn metric_product(&mut self, metric: AtomView<'_>, coefficient: i32, value: Atom) -> Atom {
        if coefficient == 1 {
            metric * value.as_view()
        } else {
            Atom::num(coefficient) * metric * value
        }
    }
    fn dimension_product(
        &mut self,
        dimension: AtomView<'_>,
        constant: i32,
        linear: i32,
        value: Atom,
    ) -> Atom {
        // Preserve the current factored grouping, including rounded domains.
        let coefficient = match (constant, linear) {
            (0, 1) => dimension.to_owned(),
            (_, 1) => dimension - Atom::num(-constant).as_view(),
            (_, -1) => Atom::num(constant) - dimension,
            _ => unreachable!(),
        };
        coefficient * value
    }
    fn add(&mut self, left: Atom, right: Atom) -> Atom {
        left + right
    }
    fn subtract(&mut self, left: Atom, right: Atom) -> Atom {
        left - right
    }
    fn sum(&mut self, values: Vec<Atom>) -> Atom {
        Atom::add_many(values)
    }
    fn finish(self, value: Atom) -> Atom {
        value
    }
}

impl<'a, Output: TraceAlgebra> PairingTrace<'a, Output> {
    fn metric(&mut self, left: usize, right: usize) -> AtomView<'_> {
        // Summed endpoints often disappear before any metric involving them is
        // needed. In particular, do not construct indexed vector components for
        // a dummy that a word identity is about to eliminate.
        let (left, right) = (left.min(right), left.max(right));
        self.metrics[left * self.width + right]
            .get_or_insert_with(|| g!(self.arguments[left], self.arguments[right]))
            .as_view()
    }

    fn supports_sparse_output(&mut self, word: &[usize]) -> bool {
        // Extend only contracted words: a free 16-factor trace already has over
        // two million pairings. Sparse emission also bounds its row buffers.
        if !self.canonical_arguments
            || (word.len() > 14 && (word.len() > 16 || self.summed_indices.is_empty()))
            || i64::try_from(self.trace_unit).is_err()
        {
            return false;
        }
        // Mixed free slots and compact vectors can normalize an indexed vector
        // through arbitrary callbacks. The factored trace-local fallback owns
        // that domain; homogeneous pairs retain normalized metric functions.
        let mut free_explicit = false;
        let mut compact = false;
        for (index, &argument) in self.arguments.iter().enumerate() {
            if self.summed_indices.contains(&index) {
                continue;
            }
            if is_minkowski_slot(argument) {
                free_explicit = true;
            } else {
                compact = true;
            }
        }
        if free_explicit && compact {
            return false;
        }
        // Symbolic dimensions become polynomial variables; literal four folds
        // into the integer coefficient. Other dimensions retain Atom arithmetic.
        self.summed_indices.first().copied().is_none_or(|index| {
            let dimension = self.metric(index, index);
            matches!(dimension, AtomView::Var(_)) || is_four_dimension(dimension)
        })
    }

    fn metric_product(
        &mut self,
        left: usize,
        right: usize,
        coefficient: i32,
        value: Output::Value,
    ) -> Output::Value {
        let (left, right) = (left.min(right), left.max(right));
        let metric = self.metrics[left * self.width + right]
            .get_or_insert_with(|| g!(self.arguments[left], self.arguments[right]));
        self.output
            .metric_product(metric.as_view(), coefficient, value)
    }

    fn with_output<Other: TraceAlgebra>(self, output: Other) -> PairingTrace<'a, Other> {
        debug_assert!(self.subwords.is_empty());
        PairingTrace {
            width: self.width,
            arguments: self.arguments,
            metrics: self.metrics,
            trace_unit: self.trace_unit,
            subwords: HashMap::new(),
            output,
            canonical_arguments: self.canonical_arguments,
            summed_indices: self.summed_indices,
        }
    }

    fn evaluate(&mut self, remaining: &[usize]) -> Output::Value {
        let Some(&first) = remaining.first() else {
            return self.output.unit(self.trace_unit);
        };
        if let Some(result) = self.subwords.get(remaining) {
            return result.clone();
        }
        let result = self
            .contract_word(remaining)
            .or_else(|| self.reduce_scalar_word(remaining))
            .unwrap_or_else(|| self.pair(remaining, first));
        // A long trace must not retain an unbounded number of large subwords.
        // These limits cover all reusable subwords through length fourteen.
        if remaining.len() <= 14 && self.subwords.len() < 1024 {
            self.subwords.insert(remaining.to_vec(), result.clone());
        }
        result
    }

    /// Eliminate a summed explicit pair before evaluating the free pairing
    /// polynomial. Reordering the interior keeps every other summed index in
    /// its word, so each branch can contract independently without an external
    /// tensor metric or a new symbolic trace node.
    fn contract_word(&mut self, word: &[usize]) -> Option<Output::Value> {
        if self.summed_indices.is_empty() {
            return None;
        }
        let mut best = None;
        for (left, &index) in word.iter().enumerate() {
            if !self.summed_indices.contains(&index) {
                continue;
            }
            let Some(right) = word[left + 1..].iter().position(|&i| i == index) else {
                continue;
            };
            let right = left + 1 + right;
            let four_dimensional = i64::try_from(self.metric(index, index)) == Ok(4);
            for (start, distance) in [(left, right - left), (right, word.len() - right + left)] {
                let priority = if distance == 1 {
                    0
                } else if distance == 2 || (four_dimensional && distance.is_multiple_of(2)) {
                    1
                } else {
                    2
                };
                let score = (priority, distance);
                if best.is_none_or(|(_, _, previous)| score < previous) {
                    best = Some((start, distance, score));
                }
            }
        }
        let (start, distance, _) = best?;
        let mut rotated = word.to_vec();
        rotated.rotate_left(start);
        let interior = &rotated[1..distance];
        let outside = &rotated[distance + 1..];
        let dimension = self.metric(rotated[0], rotated[0]).to_owned();
        let position = self
            .summed_indices
            .iter()
            .position(|&i| i == rotated[0])
            .unwrap();
        let contracted = self.summed_indices.swap_remove(position);
        // Every child omits this pair. Suspend its classification while solving
        // those children so free subwords skip the repeated-pair scan entirely.
        // Memo keys remain valid: the removed ID is absent from every child.
        let result = match interior.len() {
            0 => {
                let trace = self.evaluate(outside);
                Some(
                    self.output
                        .dimension_product(dimension.as_view(), 0, 1, trace),
                )
            }
            1 => {
                let reduced: Vec<_> = interior.iter().chain(outside).copied().collect();
                let trace = self.evaluate(&reduced);
                Some(
                    self.output
                        .dimension_product(dimension.as_view(), 2, -1, trace),
                )
            }
            n if n % 2 == 1 && i64::try_from(dimension.as_view()) == Ok(4) => {
                let reduced: Vec<_> = interior.iter().rev().chain(outside).copied().collect();
                let trace = self.evaluate(&reduced);
                Some(self.output.scale(-2, trace))
            }
            2 if interior.iter().all(|i| !self.summed_indices.contains(i)) => {
                // gamma(mu) a/ b/ gamma(mu) = 4(a.b) + (D-4)a/b/.
                // A metric with another summed index must stay in the word;
                // the general permutation recurrence below owns that case.
                let outside_trace = self.evaluate(outside);
                let first = self.metric_product(interior[0], interior[1], 4, outside_trace);
                let coefficient = dimension.as_view() - Atom::num(4).as_view();
                if coefficient.is_zero() {
                    Some(first)
                } else {
                    let reduced: Vec<_> = interior.iter().chain(outside).copied().collect();
                    let trace = self.evaluate(&reduced);
                    let second = self
                        .output
                        .dimension_product(dimension.as_view(), -4, 1, trace);
                    Some(self.output.add(first, second))
                }
            }
            n => {
                // Anticommute the left endpoint through the interior, then
                // substitute its explicit index in each anticommutator term:
                // (-1)^n(D-2) A + 2 sum_i (-1)^i A_without_i gamma(i).
                // The final i=n-1 term is included in the D-2 coefficient.
                let unchanged: Vec<_> = interior.iter().chain(outside).copied().collect();
                let coefficient = if n % 2 == 0 {
                    dimension.as_view() - Atom::num(2).as_view()
                } else {
                    Atom::num(2) - dimension.as_view()
                };
                let mut terms = Vec::with_capacity(n);
                if !coefficient.is_zero() {
                    let trace = self.evaluate(&unchanged);
                    let (constant, linear) = if n % 2 == 0 { (-2, 1) } else { (2, -1) };
                    terms.push(self.output.dimension_product(
                        dimension.as_view(),
                        constant,
                        linear,
                        trace,
                    ));
                }
                for moved in 0..n - 1 {
                    let reduced: Vec<_> = interior[..moved]
                        .iter()
                        .chain(&interior[moved + 1..])
                        .chain(std::iter::once(&interior[moved]))
                        .chain(outside)
                        .copied()
                        .collect();
                    let trace = self.evaluate(&reduced);
                    terms.push(
                        self.output
                            .scale(if moved % 2 == 0 { 2 } else { -2 }, trace),
                    );
                }
                Some(self.output.sum(terms))
            }
        };
        self.summed_indices.push(contracted);
        result
    }

    /// Compact slashes square to scalar norms. Keeping these identities inside
    /// the word recurrence avoids constructing and revisiting shorter traces.
    fn reduce_scalar_word(&mut self, word: &[usize]) -> Option<Output::Value> {
        if !self.canonical_arguments {
            return None;
        }
        for start in 0..word.len() {
            if word[start] == word[(start + 1) % word.len()] {
                let rest: Vec<_> = word
                    .iter()
                    .cycle()
                    .skip(start + 2)
                    .take(word.len() - 2)
                    .copied()
                    .collect();
                let trace = self.evaluate(&rest);
                return Some(self.metric_product(word[start], word[start], 1, trace));
            }
        }
        if word.len() >= 4 {
            for start in 0..word.len() {
                let p = word[start];
                let q = word[(start + 1) % word.len()];
                if p != word[(start + 2) % word.len()] {
                    continue;
                }
                // p/ q/ p/ = 2 (p.q) p/ - p^2 q/. All summed explicit
                // pairs were eliminated first: p is compact and a free
                // explicit q can occur only in this coefficient.
                let mut rest: Vec<_> = std::iter::once(p)
                    .chain(
                        word.iter()
                            .cycle()
                            .skip(start + 3)
                            .take(word.len() - 3)
                            .copied(),
                    )
                    .collect();
                let first = self.evaluate(&rest);
                rest[0] = q;
                let second = self.evaluate(&rest);
                let first = self.metric_product(p, q, 2, first);
                let second = self.metric_product(p, p, 1, second);
                return Some(self.output.subtract(first, second));
            }
        }
        None
    }

    fn pair(&mut self, remaining: &[usize], first: usize) -> Output::Value {
        let mut terms = Vec::with_capacity(remaining.len() - 1);
        for partner in 1..remaining.len() {
            let rest: Vec<_> = remaining[1..partner]
                .iter()
                .chain(&remaining[partner + 1..])
                .copied()
                .collect();
            let subword = self.evaluate(&rest);
            let product = self.metric_product(first, remaining[partner], 1, subword);
            terms.push(
                self.output
                    .scale(if partner % 2 == 1 { 1 } else { -1 }, product),
            );
        }
        self.output.sum(terms)
    }
}

/// `canonical_arguments` admits only canonical compact vectors and explicit
/// slots with compatible dimensions and at most two occurrences. It enables
/// word contractions; unrestricted arguments retain the plain pairing formula.
pub(super) fn evaluate_generic(
    indices: &[AtomView<'_>],
    trace_unit: AtomView<'_>,
    canonical_arguments: bool,
    output: TraceOutput,
) -> Atom {
    if indices.len() % 2 == 1 {
        return Atom::Zero;
    }
    // Ordered argument identities, including their dimension and metadata,
    // determine a subtrace. Equal arguments can share states across positions.
    let mut arguments = Vec::new();
    let word: Vec<_> = indices
        .iter()
        .map(|index| {
            arguments
                .iter()
                .position(|known| known == index)
                .unwrap_or_else(|| {
                    arguments.push(*index);
                    arguments.len() - 1
                })
        })
        .collect();
    let width = arguments.len();
    let summed_indices = if canonical_arguments {
        arguments
            .iter()
            .enumerate()
            .filter_map(|(index, &argument)| {
                (is_minkowski_slot(argument) && word.iter().filter(|&&i| i == index).count() == 2)
                    .then_some(index)
            })
            .collect()
    } else {
        Vec::new()
    };
    let mut trace = PairingTrace {
        width,
        arguments,
        metrics: vec![None; width * width],
        trace_unit,
        subwords: HashMap::new(),
        output: FactoredOutput,
        canonical_arguments,
        summed_indices,
    };
    if matches!(output, TraceOutput::Expanded) && trace.supports_sparse_output(&word) {
        let mut trace = trace.with_output(sparse::SparseOutput::default());
        let result = trace.evaluate(&word);
        trace.output.finish(result)
    } else {
        let result = trace.evaluate(&word);
        match output {
            TraceOutput::Factored => result,
            TraceOutput::Expanded => result.expand(),
        }
    }
}

impl Monomial {
    fn collect<const N: usize>(mut self, mut sign: i32, terms: &mut BTreeMap<[u8; N], i32>) {
        // Eliminate temporary metric-linked indices before expanding epsilon
        // pairs, so two epsilons sharing a dummy use a 3x3 determinant.
        while let Some(position) = self
            .metrics
            .iter()
            .position(|p| p.iter().any(|&i| usize::from(i) >= N))
        {
            let [mut a, mut b] = self.metrics.swap_remove(position);
            if a == b {
                sign *= 4;
                continue;
            }
            if usize::from(a) < N {
                std::mem::swap(&mut a, &mut b);
            }
            for index in self
                .metrics
                .iter_mut()
                .flatten()
                .chain(self.epsilons.iter_mut().flatten())
            {
                if *index == a {
                    *index = b;
                }
            }
        }
        if self
            .epsilons
            .iter()
            .any(|e| (0..4).any(|i| e[i + 1..].contains(&e[i])))
        {
            return;
        }
        if self.epsilons.len() > 1 {
            let shared = self.epsilons.iter().enumerate().find_map(|(i, left)| {
                self.epsilons
                    .iter()
                    .enumerate()
                    .skip(i + 1)
                    .find_map(|(j, right)| {
                        left.iter().enumerate().find_map(|(a, index)| {
                            right.iter().position(|x| x == index).map(|b| (i, j, a, b))
                        })
                    })
            });
            let (i, j) = shared.map_or((0, 1), |(i, j, _, _)| (i, j));
            let right = self.epsilons.remove(j);
            let left = self.epsilons.remove(i);
            let (left, right) = if let Some((_, _, a, b)) = shared {
                if (a + b) % 2 == 1 {
                    sign = -sign;
                }
                (
                    left.into_iter()
                        .enumerate()
                        .filter_map(|(i, x)| (i != a).then_some(x))
                        .collect::<Vec<_>>(),
                    right
                        .into_iter()
                        .enumerate()
                        .filter_map(|(i, x)| (i != b).then_some(x))
                        .collect::<Vec<_>>(),
                )
            } else {
                (left.to_vec(), right.to_vec())
            };
            for permutation in (0..left.len()).permutations(left.len()) {
                let inversions = permutation
                    .iter()
                    .enumerate()
                    .map(|(i, a)| permutation[i + 1..].iter().filter(|b| a > *b).count())
                    .sum::<usize>();
                let mut next = self.clone();
                next.metrics.extend(
                    left.iter()
                        .enumerate()
                        .map(|(i, &a)| [a, right[permutation[i]]]),
                );
                next.collect::<N>(if inversions % 2 == 0 { sign } else { -sign }, terms);
            }
            return;
        }
        let mut key = [0; N];
        let offset = if let Some(epsilon) = self.epsilons.first_mut() {
            for i in 0..4 {
                for j in i + 1..4 {
                    if epsilon[i] > epsilon[j] {
                        sign = -sign;
                    }
                }
            }
            epsilon.sort_unstable();
            key[..4].copy_from_slice(epsilon);
            4
        } else {
            0
        };
        for pair in &mut self.metrics {
            pair.sort_unstable();
        }
        self.metrics.sort_unstable();
        for (destination, index) in key[offset..].iter_mut().zip(self.metrics.iter().flatten()) {
            debug_assert!(usize::from(*index) < N);
            *destination = *index;
        }
        *terms.entry(key).or_default() += sign;
    }
}

// Separate statics per arity and parity: Rust statics inside a generic function
// would otherwise be shared across monomorphizations. Odd traces short-circuit
// before constructing a kernel. Tables are initialized only when needed.
macro_rules! short_trace_dispatch {
    ($($length:literal),* $(,)?) => {
        pub(super) fn evaluate(indices: &[AtomView<'_>], axial: bool, output: TraceOutput) -> Option<Atom> {
            match indices.len() {
                1 | 3 | 5 | 7 | 9 | 11 | 13 => Some(Atom::Zero),
                $($length => {
                    static ORDINARY: LazyLock<TraceKernel<$length>> = LazyLock::new(|| TraceKernel::generate(false));
                    static AXIAL: LazyLock<TraceKernel<$length>> = LazyLock::new(|| TraceKernel::generate(true));
                    Some(if axial { AXIAL.evaluate(indices, true, output) } else { ORDINARY.evaluate(indices, false, output) })
                },)*
                _ => None,
            }
        }

    };
}

short_trace_dispatch!(2, 4, 6, 8, 10, 12, 14);

#[cfg(test)]
pub(super) mod tests {
    use super::*;

    // Independent Clifford-algebra oracle in a Euclidean orthonormal basis.
    // Store all basis blades and multiply by each vector, without trace
    // recurrences or epsilon contraction identities from the implementation.
    pub(in crate::dirac::simplify) fn clifford_trace<const D: usize>(
        vectors: &[[i64; D]],
        axial: bool,
    ) -> i64 {
        assert!(!axial || D == 4);
        let mut blades = vec![0; 1 << D];
        blades[0] = 1;
        for vector in vectors {
            let mut next = vec![0; 1 << D];
            for (mask, coefficient) in blades.into_iter().enumerate() {
                for (axis, component) in vector.iter().enumerate() {
                    let sign = if (mask >> (axis + 1)).count_ones() % 2 == 0 {
                        1
                    } else {
                        -1
                    };
                    next[mask ^ (1 << axis)] += sign * coefficient * component;
                }
            }
            blades = next;
        }
        4 * blades[if axial { 15 } else { 0 }]
    }

    fn determinant(rows: [[i64; 4]; 4]) -> i64 {
        (0..4)
            .permutations(4)
            .map(|p| {
                let inversions = (0..4)
                    .map(|i| (i + 1..4).filter(|&j| p[i] > p[j]).count())
                    .sum::<usize>();
                let sign = if inversions % 2 == 0 { 1 } else { -1 };
                sign * (0..4).map(|i| rows[i][p[i]]).product::<i64>()
            })
            .sum()
    }

    fn check<const N: usize>(ordinary_terms: usize, axial_terms: usize) {
        for (axial, expected_terms) in [(false, ordinary_terms), (true, axial_terms)] {
            let kernel = TraceKernel::<N>::generate(axial);
            assert_eq!(kernel.terms.len(), expected_terms);
            // Terminal free-trace evaluation relies on each original index
            // appearing once, across the epsilon and all metric pairs.
            for (recipe, _) in &kernel.terms {
                let mut indices = *recipe;
                indices.sort_unstable();
                assert_eq!(indices, std::array::from_fn(|index| index as u8));
            }
            if N >= 10 {
                let factored = kernel.factored(axial);
                let mut expected = BTreeMap::<Vec<usize>, i32>::new();
                for (recipe, coefficient) in &kernel.terms {
                    let mut factors: Vec<_> = recipe[if axial { 4 } else { 0 }..]
                        .as_chunks::<2>()
                        .0
                        .iter()
                        .map(|pair| usize::from(pair[0]) * N + usize::from(pair[1]))
                        .collect();
                    if axial {
                        let epsilon: [u8; 4] = recipe[..4].try_into().unwrap();
                        factors.push(
                            N * N
                                + factored
                                    .epsilons
                                    .iter()
                                    .position(|key| *key == epsilon)
                                    .unwrap(),
                        );
                    }
                    factors.sort_unstable();
                    *expected.entry(factors).or_default() += coefficient;
                }
                assert_eq!(factored.polynomial.coefficient_map(), expected);
            }
            let mut seed = 17u64;
            for sample in 0..8 {
                let vectors: [[i64; 4]; N] = std::array::from_fn(|_| {
                    std::array::from_fn(|_| {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((seed >> 32) % 5) as i64 - 2
                    })
                });
                let result: i64 = kernel
                    .terms
                    .iter()
                    .map(|(recipe, coefficient)| {
                        let (epsilon, pairs) = if axial {
                            (
                                determinant(std::array::from_fn(|i| {
                                    vectors[usize::from(recipe[i])]
                                })),
                                &recipe[4..],
                            )
                        } else {
                            (1, &recipe[..])
                        };
                        4 * i64::from(*coefficient)
                            * epsilon
                            * pairs
                                .as_chunks::<2>()
                                .0
                                .iter()
                                .map(|pair| {
                                    vectors[usize::from(pair[0])]
                                        .iter()
                                        .zip(vectors[usize::from(pair[1])])
                                        .map(|(a, b)| a * b)
                                        .sum::<i64>()
                                })
                                .product::<i64>()
                    })
                    .sum();
                assert_eq!(
                    result,
                    clifford_trace(&vectors, axial),
                    "length {N}, axial {axial}, sample {sample}"
                );
            }
        }
    }

    #[test]
    fn short_kernels_match_independent_clifford_products() {
        check::<2>(1, 0);
        check::<4>(3, 1);
        check::<6>(15, 6);
        check::<8>(105, 33);
        check::<10>(693, 180);
        check::<12>(4383, 1029);
        check::<14>(26931, 6042);
    }

    fn check_generic<const N: usize>() {
        // Six independent axes distinguish the generic formula from a 4D
        // reduction. Keep the formal unit trace at four in both evaluations.
        let mut seed = 23u64;
        for sample in 0..8 {
            let vectors: [[i64; 6]; N] = std::array::from_fn(|_| {
                std::array::from_fn(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    ((seed >> 32) % 5) as i64 - 2
                })
            });
            let unit = Atom::num(4);
            let mut trace = PairingTrace {
                width: N,
                arguments: Vec::new(),
                metrics: (0..N)
                    .flat_map(|a| {
                        (0..N).map(move |b| {
                            Atom::num(
                                vectors[a]
                                    .iter()
                                    .zip(vectors[b])
                                    .map(|(a, b)| a * b)
                                    .sum::<i64>(),
                            )
                        })
                    })
                    .map(Some)
                    .collect(),
                trace_unit: unit.as_view(),
                subwords: HashMap::new(),
                output: FactoredOutput,
                canonical_arguments: false,
                summed_indices: Vec::new(),
            };
            let result = trace.evaluate(&(0..N).collect::<Vec<_>>());
            assert_eq!(
                result,
                Atom::num(clifford_trace(&vectors, false)),
                "generic length {N}, sample {sample}"
            );
            assert!(trace.subwords.len() <= 1024);
        }
    }

    #[test]
    fn generic_kernels_match_independent_six_dimensional_clifford_products() {
        check_generic::<2>();
        check_generic::<4>();
        check_generic::<6>();
        check_generic::<8>();
        check_generic::<10>();
        check_generic::<12>();
        check_generic::<14>();
        check_generic::<16>();
    }

    #[test]
    fn compact_word_reductions_match_six_dimensional_clifford_products() {
        let words: Vec<Vec<usize>> = vec![
            (0..14).map(|i| usize::from(i >= 2)).collect(),
            (0..14).map(|i| i % 2).collect(),
            vec![0, 1, 2, 3, 4, 5, 0, 1],
            vec![0, 1, 2, 3, 4, 5, 1, 0],
            vec![0, 1, 2, 3, 4, 5, 0, 5],
        ];
        let mut seed = 31u64;
        for sample in 0..8 {
            let vectors: [[i64; 6]; 6] = std::array::from_fn(|_| {
                std::array::from_fn(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    ((seed >> 32) % 5) as i64 - 2
                })
            });
            let unit = Atom::num(4);
            let metrics: Vec<_> = (0..6)
                .flat_map(|a| {
                    (0..6).map(move |b| {
                        Atom::num(
                            vectors[a]
                                .iter()
                                .zip(vectors[b])
                                .map(|(a, b)| a * b)
                                .sum::<i64>(),
                        )
                    })
                })
                .collect();
            for word in &words {
                let mut trace = PairingTrace {
                    width: 6,
                    arguments: Vec::new(),
                    metrics: metrics.iter().cloned().map(Some).collect(),
                    trace_unit: unit.as_view(),
                    subwords: HashMap::new(),
                    output: FactoredOutput,
                    canonical_arguments: true,
                    summed_indices: Vec::new(),
                };
                let factors: Vec<_> = word.iter().map(|&index| vectors[index]).collect();
                assert_eq!(
                    trace.evaluate(word),
                    Atom::num(clifford_trace(&factors, false)),
                    "word {word:?}, sample {sample}"
                );
            }
        }
    }

    fn check_summed_word<const D: usize>() {
        use super::super::DiracSimplifier;
        use crate::gamma;
        use spenso::{
            network::{library::symbolic::ETS, tags::SPENSO_TAG},
            structure::representation::{Minkowski, RepName},
            trace,
        };
        use symbolica::{
            atom::{AtomCore, FunctionBuilder},
            symbol,
        };

        let reps = crate::test_support::test_initialize();
        let mink = Minkowski {}.new_rep(D);
        let spin = reps.bis4.to_symbolic([]);
        let arguments: Vec<_> = (0..12)
            .map(|i| {
                FunctionBuilder::new(
                    SPENSO_TAG.rank_one_tensor_symbol(&format!("idenso::summed_trace::p{i}")),
                )
                .add_arg(mink.to_symbolic([]))
                .finish()
            })
            .chain([
                mink.pattern(symbol!("summed_trace_a")),
                mink.pattern(symbol!("summed_trace_b")),
            ])
            .collect();
        let mut words = vec![
            vec![12, 0, 1, 2, 3, 12],
            vec![12, 0, 12, 1, 2, 3],
            vec![12, 0, 13, 13, 12, 1, 2, 3],
            vec![12, 0, 13, 1, 13, 2, 12, 3],
            vec![12, 0, 13, 1, 12, 2, 13, 3],
            vec![12, 0, 13, 1, 2, 12, 3, 4, 13, 5, 6, 7],
            vec![12, 0, 1, 0, 1, 12, 0, 1],
        ];
        for interior in 2..=5 {
            words.push(
                std::iter::once(12)
                    .chain(0..interior)
                    .chain(std::iter::once(12))
                    .chain(interior..2 * interior)
                    .collect(),
            );
        }
        let mut seed = 43u64;
        for sample in 0..4 {
            let vectors: [[i64; D]; 12] = std::array::from_fn(|_| {
                std::array::from_fn(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    ((seed >> 32) % 5) as i64 - 2
                })
            });
            for word in &words {
                let input = trace!(&spin; word.iter().map(|&i| gamma!(&arguments[i])));
                let reduced = DiracSimplifier::new(&crate::dirac::GammaSimplifySettings::default())
                    .evaluate_terminal_trace::<false>(input.as_view())
                    .expect("scalar word contractions must finish within the terminal evaluator");
                reduced.visitor(&mut |node| {
                    assert!(
                        node != arguments[12].as_view() && node != arguments[13].as_view(),
                        "summed index remains in word {word:?}"
                    );
                    true
                });
                let result = reduced.replace_map(|atom, _, out| {
                    if let AtomView::Fun(metric) = atom
                        && metric.get_symbol() == ETS.metric
                    {
                        let args: Vec<_> = metric
                            .iter()
                            .map(|arg| {
                                arguments[..12]
                                    .iter()
                                    .position(|vector| vector.as_view() == arg)
                                    .expect("every explicit summed index must have been eliminated")
                            })
                            .collect();
                        **out = Atom::num(
                            vectors[args[0]]
                                .iter()
                                .zip(vectors[args[1]])
                                .map(|(a, b)| a * b)
                                .sum::<i64>(),
                        );
                    }
                });
                let mut expected = 0;
                for a in 0..D {
                    for b in 0..if word.contains(&13) { D } else { 1 } {
                        let factors: Vec<_> = word
                            .iter()
                            .map(|&i| match i {
                                12 => std::array::from_fn(|axis| i64::from(axis == a)),
                                13 => std::array::from_fn(|axis| i64::from(axis == b)),
                                _ => vectors[i],
                            })
                            .collect();
                        expected += clifford_trace(&factors, false);
                    }
                }
                assert_eq!(
                    result,
                    Atom::num(expected),
                    "D={D}, sample {sample}, word {word:?}"
                );
            }
        }
    }

    #[test]
    fn summed_trace_words_match_independent_clifford_products() {
        check_summed_word::<4>();
        check_summed_word::<6>();
    }

    #[test]
    fn long_axial_fallback_matches_a_clifford_component() {
        use super::super::{DiracFactor, DiracSimplifier};
        use crate::{
            dirac::GammaSimplifier, epsilon::EpsilonSimplifier, gamma, gamma5,
            shorthands::schoonschip::Schoonschip,
        };
        use spenso::{network::library::symbolic::ETS, p, q, vector};
        use symbolica::{atom::AtomCore, symbol};

        let reps = crate::test_support::test_initialize();
        let spin = reps.bis4.to_symbolic([]);
        let indices: Vec<_> = (0..16)
            .map(|i| reps.mink4.pattern(symbol!(&format!("long_axial_mu{i}"))))
            .collect();
        let factors: Vec<_> = [gamma5!()]
            .into_iter()
            .chain(indices.iter().map(|mu| gamma!(mu)))
            .collect();
        let parsed: Vec<_> = factors
            .iter()
            .map(|factor| DiracFactor::parse(factor.as_view()))
            .collect();
        // Force one long-trace step before projection so the fallback is
        // exercised. Projecting first would only test repeated-pair reductions.
        let mut reduced =
            DiracSimplifier::simplify_single_gamma5_trace(spin.as_view(), &parsed, 0).unwrap();
        let basis = [
            p!(reps.mink4.to_symbolic([])),
            q!(reps.mink4.to_symbolic([])),
            vector!(long_axial_r, reps.mink4.to_symbolic([])),
            vector!(long_axial_s, reps.mink4.to_symbolic([])),
        ];
        let axes = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 0, 1, 1];
        for (index, axis) in indices.iter().zip(axes) {
            reduced = reduced
                .replace(index.to_pattern())
                .with(basis[axis].to_pattern());
        }
        let result = reduced
            .simplify_gamma()
            .expand()
            .simplify_epsilon()
            .schoonschip()
            .replace_map(|atom, _, out| {
                if let AtomView::Fun(f) = atom
                    && (f.get_symbol() == ETS.metric
                        || f.get_symbol() == *crate::epsilon::EPSILON_SYMBOL)
                {
                    let axes: Vec<_> = f
                        .iter()
                        .map(|arg| {
                            basis
                                .iter()
                                .position(|vector| vector.as_view() == arg)
                                .expect("all dummy indices must contract")
                        })
                        .collect();
                    **out = if axes.len() == 2 {
                        Atom::num(i64::from(axes[0] == axes[1]))
                    } else {
                        Atom::num(determinant(std::array::from_fn(|i| {
                            std::array::from_fn(|j| i64::from(axes[i] == j))
                        })))
                    };
                }
            });
        let vectors: Vec<[i64; 4]> = axes
            .iter()
            .map(|&axis| std::array::from_fn(|i| i64::from(i == axis)))
            .collect();
        assert_eq!(result, Atom::num(clifford_trace(&vectors, true)));
        assert_eq!(result, Atom::num(4));
    }
}
