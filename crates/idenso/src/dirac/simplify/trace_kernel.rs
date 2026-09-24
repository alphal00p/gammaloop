//! Dimension-generic traces reuse factored pairing polynomials within one call.
//! Short four-dimensional traces use more compact integer recipes generated
//! from the same three-gamma identity as open chains. No user atoms or dummy
//! indices escape into the global recipe cache.

use std::{collections::{BTreeMap, HashMap}, sync::{LazyLock, OnceLock}};

use itertools::Itertools;
use spenso::g;
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::THREE_GAMMA_METRIC_TERMS;
use crate::epsilon::epsilon4;

mod factored;

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
                epsilon4(indices[a as usize], indices[b as usize], indices[c as usize], indices[d as usize])
            }));
            return factored.polynomial.evaluate(&metrics, Atom::num(4).as_view());
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
                    last_epsilon = epsilon4(indices[a as usize], indices[b as usize], indices[c as usize], indices[d as usize]);
                }
                Some(last_epsilon.as_view())
            } else {
                None
            };
            let pairs = &indices_recipe[if axial { 4 } else { 0 }..];
            Atom::mul_many(std::iter::once(coefficients[(*coefficient + max_coefficient) as usize].as_view())
                .chain(epsilon)
                .chain(
                pairs
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|p| metrics[usize::from(p[0]) * N + usize::from(p[1])].as_view()),
            ))
        }))
    }
}

/// Build the ordinary pairing polynomial without intermediate trace nodes.
/// Repeated subwords share their result for this evaluation only.
struct PairingTrace<'a> {
    width: usize,
    metrics: Vec<Atom>,
    trace_unit: AtomView<'a>,
    subwords: HashMap<Vec<usize>, Atom>,
}

impl PairingTrace<'_> {
    fn evaluate(&mut self, remaining: &[usize]) -> Atom {
        let Some(&first) = remaining.first() else {
            return self.trace_unit.to_owned();
        };
        if let Some(result) = self.subwords.get(remaining) {
            return result.clone();
        }
        let mut terms = Vec::with_capacity(remaining.len() - 1);
        for partner in 1..remaining.len() {
            let rest: Vec<_> = remaining[1..partner]
                .iter()
                .chain(&remaining[partner + 1..])
                .copied()
                .collect();
            let subword = self.evaluate(&rest);
            let product = self.metrics[first * self.width + remaining[partner]].as_view()
                * subword.as_view();
            terms.push(if partner % 2 == 1 { product } else { -product });
        }
        let result = Atom::add_many(terms);
        // A long trace must not retain an unbounded number of large subwords.
        // These limits cover all reusable subwords through length fourteen.
        if remaining.len() <= 14 && self.subwords.len() < 1024 {
            self.subwords.insert(remaining.to_vec(), result.clone());
        }
        result
    }
}

pub(super) fn evaluate_generic(indices: &[AtomView<'_>], trace_unit: AtomView<'_>) -> Atom {
    if indices.len() % 2 == 1 {
        return Atom::Zero;
    }
    let width = indices.len();
    let metrics = (0..width)
        .flat_map(|a| (0..width).map(move |b| g!(indices[a], indices[b])))
        .collect();
    PairingTrace {
        width,
        metrics,
        trace_unit,
        subwords: HashMap::new(),
    }
    .evaluate(&(0..width).collect::<Vec<_>>())
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
mod tests {
    use super::*;

    // Independent Clifford-algebra oracle in a Euclidean orthonormal basis.
    // Store all basis blades and multiply by each vector, without trace
    // recurrences or epsilon contraction identities from the implementation.
    fn clifford_trace<const D: usize>(vectors: &[[i64; D]], axial: bool) -> i64 {
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
                metrics: (0..N)
                    .flat_map(|a| {
                        (0..N).map(move |b| {
                            Atom::num(vectors[a].iter().zip(vectors[b]).map(|(a, b)| a * b).sum::<i64>())
                        })
                    })
                    .collect(),
                trace_unit: unit.as_view(),
                subwords: HashMap::new(),
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
