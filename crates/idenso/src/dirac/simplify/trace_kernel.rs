//! Short trace recipes: signed pairings in arbitrary dimension, and the more
//! compact four-dimensional reduction from the same three-gamma identity as
//! open chains. Only integer index recipes are cached; no user atoms or dummy
//! indices escape into the global cache.

use std::{collections::BTreeMap, sync::LazyLock};

use itertools::Itertools;
use spenso::g;
use symbolica::atom::{Atom, AtomView};

use super::THREE_GAMMA_METRIC_TERMS;
use crate::epsilon::epsilon4;

pub(super) struct TraceKernel<const N: usize> {
    // For axial traces, the first four entries belong to epsilon; the rest
    // are sorted metric pairs. Ordinary traces contain metric pairs only.
    terms: Vec<([u8; N], i32)>,
}

#[derive(Clone, Default)]
struct Monomial {
    metrics: Vec<[u8; 2]>,
    epsilons: Vec<[u8; 4]>,
}

impl<const N: usize> TraceKernel<N> {
    fn generate_pairings() -> Self {
        let mut terms = Vec::with_capacity((1..N).step_by(2).product());
        Self::pairings((1 << N) - 1, 0, &mut [0; N], 1, &mut terms);
        Self { terms }
    }

    fn pairings(
        remaining: u16,
        position: usize,
        recipe: &mut [u8; N],
        mut sign: i32,
        terms: &mut Vec<([u8; N], i32)>,
    ) {
        if remaining == 0 {
            terms.push((*recipe, sign));
            return;
        }
        let first = remaining.trailing_zeros();
        let remaining = remaining & !(1 << first);
        recipe[position] = first as u8;
        let mut partners = remaining;
        while partners != 0 {
            let partner = partners.trailing_zeros();
            recipe[position + 1] = partner as u8;
            Self::pairings(
                remaining & !(1 << partner),
                position + 2,
                recipe,
                sign,
                terms,
            );
            partners &= partners - 1;
            sign = -sign;
        }
    }

    fn generate(axial: bool) -> Self {
        let mut terms = BTreeMap::new();
        Self::reduce_triples(0, 1, 1, Monomial::default(), axial, &mut terms);
        Self {
            terms: terms.into_iter().filter(|(_, c)| *c != 0).collect(),
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

    fn evaluate(&self, indices: &[AtomView<'_>], axial: bool, trace_unit: AtomView<'_>) -> Atom {
        debug_assert_eq!(indices.len(), N);
        // Each scalar product is built once and shared by all template terms.
        let metrics: Vec<_> = (0..N)
            .flat_map(|a| (0..N).map(move |b| g!(indices[a], indices[b])))
            .collect();
        Atom::add_many(self.terms.iter().map(|(indices_recipe, coefficient)| {
            let mut factors = Vec::with_capacity(N / 2 + 1);
            factors.push(Atom::num(i64::from(*coefficient)) * trace_unit);
            let pairs = if axial {
                let [a, b, c, d]: [u8; 4] = indices_recipe[..4].try_into().unwrap();
                factors.push(epsilon4(
                    indices[a as usize],
                    indices[b as usize],
                    indices[c as usize],
                    indices[d as usize],
                ));
                &indices_recipe[4..]
            } else {
                &indices_recipe[..]
            };
            factors.extend(
                pairs
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|p| metrics[usize::from(p[0]) * N + usize::from(p[1])].clone()),
            );
            Atom::mul_many(factors)
        }))
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
        pub(super) fn evaluate(indices: &[AtomView<'_>], axial: bool) -> Option<Atom> {
            let trace_unit = Atom::num(4);
            match indices.len() {
                1 | 3 | 5 | 7 | 9 | 11 | 13 => Some(Atom::Zero),
                $($length => {
                    static ORDINARY: LazyLock<TraceKernel<$length>> = LazyLock::new(|| TraceKernel::generate(false));
                    static AXIAL: LazyLock<TraceKernel<$length>> = LazyLock::new(|| TraceKernel::generate(true));
                    Some(if axial { AXIAL.evaluate(indices, true, trace_unit.as_view()) } else { ORDINARY.evaluate(indices, false, trace_unit.as_view()) })
                },)*
                _ => None,
            }
        }

        pub(super) fn evaluate_generic(indices: &[AtomView<'_>], trace_unit: AtomView<'_>) -> Option<Atom> {
            match indices.len() {
                1 | 3 | 5 | 7 | 9 | 11 | 13 => Some(Atom::Zero),
                $($length => {
                    static GENERIC: LazyLock<TraceKernel<$length>> = LazyLock::new(TraceKernel::generate_pairings);
                    Some(GENERIC.evaluate(indices, false, trace_unit))
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

    fn check_generic<const N: usize>(expected_terms: usize) {
        let kernel = TraceKernel::<N>::generate_pairings();
        assert_eq!(kernel.terms.len(), expected_terms);
        for (recipe, coefficient) in &kernel.terms {
            assert!(coefficient.abs() == 1);
            let mut indices = *recipe;
            indices.sort_unstable();
            assert_eq!(indices, std::array::from_fn(|index| index as u8));
        }
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
            let result: i64 = kernel
                .terms
                .iter()
                .map(|(recipe, coefficient)| {
                    4 * i64::from(*coefficient)
                        * recipe
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
                clifford_trace(&vectors, false),
                "generic length {N}, sample {sample}"
            );
        }
    }

    #[test]
    fn generic_kernels_match_independent_six_dimensional_clifford_products() {
        check_generic::<2>(1);
        check_generic::<4>(3);
        check_generic::<6>(15);
        check_generic::<8>(105);
        check_generic::<10>(945);
        check_generic::<12>(10395);
        check_generic::<14>(135135);
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
