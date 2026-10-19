use spenso::{
    network::{
        library::symbolic::ETS,
        parsing::{AtomStructureExt, StrictTensorFilter},
    },
    structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        representation::{LibraryRep, RepName, Representation},
        slot::Slot,
    },
    tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
};

#[cfg(test)]
use crate::shorthands::schoonschip::{SchoonschipSettings, SchoonschipWithSettings};

spenso::symbolica_init_lazy_static! {
    /// Symbolica-level Levi-Civita symbol. The `Antisymmetric` attribute lets
    /// Symbolica canonicalize argument order and annihilate repeated arguments.
    ///
    /// In four Minkowski dimensions its normalization is fixed by the identities
    /// that produce it: the axial trace `Tr[γ5 γ^μ γ^ν γ^ρ γ^σ] = 4 ε^{μνρσ}`,
    /// the three-gamma expansion with `-ε^{μνρσ} γ_σ γ5`, and the pair rule
    /// `ε(a…) ε(b…) = det g(a_i, b_j)`. With the metric `(+,−,−,−)` and
    /// `γ5 = i γ^0 γ^1 γ^2 γ^3` (`γ5² = 1`, `diag(−1,−1,1,1)` in the Weyl basis),
    /// they require `ε^{0123} = −i`: the tensor is `−i` times the real
    /// Levi-Civita symbol with `ε^{0123} = +1`. [`contract_epsilon4`] evaluates
    /// it on explicit components.
    pub static EPSILON_SYMBOL, EPSILON_SYMBOL_INNER: Symbol =
        || tensor_symbol!("spenso::epsilon"; Antisymmetric);
}

/// Builds an antisymmetric epsilon tensor with arbitrary rank.
pub fn epsilon<'a, F>(factors: impl IntoIterator<Item = F>) -> Atom
where
    F: Into<AtomOrView<'a>>,
{
    factors
        .into_iter()
        .fold(FunctionBuilder::new(*EPSILON_SYMBOL), |builder, factor| {
            builder.add_arg(factor)
        })
        .finish()
}

/// Builds an antisymmetric epsilon tensor with arbitrary rank.
///
/// Arguments are converted through `spenso::shadowing::IntoAtom`, so tests
/// and rewrite code can pass slots, atoms, atom views, symbols, and integers
/// without spelling out the conversion.
#[macro_export]
macro_rules! epsilon {
    (; $factors:expr $(,)?) => {
        $crate::epsilon::epsilon(($factors).into_iter())
    };
    ($($arg:expr),* $(,)?) => {{
        let builder = symbolica::atom::FunctionBuilder::new(*$crate::epsilon::EPSILON_SYMBOL);
        $(
            let builder = builder.add_arg(spenso::shadowing::IntoAtom::into_atom($arg));
        )*
        builder.finish()
    }};
}

/// Builds a four-dimensional Levi-Civita tensor.
pub fn epsilon4<'a, 'b, 'c, 'd>(
    mu: impl Into<AtomOrView<'a>>,
    nu: impl Into<AtomOrView<'b>>,
    rho: impl Into<AtomOrView<'c>>,
    sigma: impl Into<AtomOrView<'d>>,
) -> Atom {
    FunctionBuilder::new(*EPSILON_SYMBOL)
        .add_arg(mu)
        .add_arg(nu)
        .add_arg(rho)
        .add_arg(sigma)
        .finish()
}

/// Contracts the four-dimensional Minkowski epsilon with four covectors given
/// by components, `ε^{ν1ν2ν3ν4} r1[ν1] r2[ν2] r3[ν3] r4[ν4] = ε^{0123} det(r)`
/// with `ε^{0123} = −i` (see [`EPSILON_SYMBOL`]).
///
/// A row holds the covariant components `a_ν = g_νμ a^μ` of a contracted
/// vector, or a unit row `δ_ν^μ` that keeps the contravariant index `μ`. For
/// four vectors with contravariant components, `ε(a,b,c,d) = i det[a;b;c;d]`.
pub fn contract_epsilon4(rows: &[[Atom; 4]; 4]) -> Atom {
    -Atom::i() * EpsilonSimplifierPass::determinant(4, |i, j| rows[i][j].clone())
}

/// Simplifies epsilon-metric contractions and epsilon-pair contractions.
#[cfg(test)]
pub(crate) trait EpsilonSimplifier {
    fn simplify_epsilon(&self) -> Atom;
}

#[cfg(test)]
impl EpsilonSimplifier for Atom {
    fn simplify_epsilon(&self) -> Atom {
        self.as_view().simplify_epsilon()
    }
}

#[cfg(test)]
impl EpsilonSimplifier for AtomView<'_> {
    fn simplify_epsilon(&self) -> Atom {
        EpsilonSimplifierPass::run(*self)
    }
}

pub(crate) struct EpsilonSimplifierPass;

impl EpsilonSimplifierPass {
    /// Apply epsilon identities once. The shared tensor scheduler contracts
    /// emitted metrics at the replacement boundary.
    pub(crate) fn step(expr: AtomView<'_>) -> Atom {
        #[cfg(feature = "reference-cases")]
        let _phase = crate::reference_cases::timing::scope(
            crate::reference_cases::timing::Phase::EpsilonKernel,
        );
        let _ = *EPSILON_SYMBOL;
        expr.replace_map(|term, _context, out| {
            let rewritten = match term {
                AtomView::Pow(_) => Self::simplify_power(term),
                AtomView::Mul(_) => Self::simplify_pair_product(term),
                _ => None,
            };
            if let Some(rewritten) = rewritten {
                **out = rewritten;
            }
        })
    }

    /// Runs the epsilon identities to a fixed point, then lets metric
    /// simplification consume any determinants that collapsed to traces.
    #[cfg(test)]
    fn run(expr: AtomView<'_>) -> Atom {
        // Force Symbolica to register epsilon as antisymmetric before any
        // builders below create new epsilon nodes.
        let _ = *EPSILON_SYMBOL;
        if !Self::contains_epsilon(expr) {
            return expr.to_owned();
        }

        let settings = SchoonschipSettings::default();
        let reducer = SchoonschipWithSettings {
            settings: &settings,
        };
        let mut current = reducer.run(expr);

        loop {
            let next = Self::step(current.as_view());
            if next == current {
                return next;
            }

            // Epsilon identities introduce determinants whose metrics may now
            // contract. Unchanged visits need no further tensor normalization.
            let next = reducer.run(next.as_view());
            if next == current {
                return next;
            }

            current = next;
        }
    }

    #[cfg(test)]
    fn contains_epsilon(expr: AtomView<'_>) -> bool {
        match expr {
            AtomView::Fun(f) => {
                f.get_symbol() == *EPSILON_SYMBOL || f.iter().any(Self::contains_epsilon)
            }
            AtomView::Add(add) => add.iter().any(Self::contains_epsilon),
            AtomView::Mul(mul) => mul.iter().any(Self::contains_epsilon),
            AtomView::Pow(pow) => pow.iter().any(Self::contains_epsilon),
            AtomView::Num(_) | AtomView::Var(_) => false,
        }
    }

    /// Applies the repeated-epsilon identity to powers:
    /// `epsilon(a_1,...,a_n)^k -> det(g(a_i,a_j)) epsilon(a_1,...,a_n)^(k-2)`.
    ///
    /// This handles the common canonical form where Symbolica has already
    /// merged an epsilon product into a positive integer power.
    fn simplify_power(expr: AtomView<'_>) -> Option<Atom> {
        let AtomView::Pow(pow) = expr else {
            return None;
        };
        let (base, exponent) = pow.get_base_exp();
        let exponent = Self::positive_integer(exponent)?;
        if exponent < 2 {
            return None;
        }
        let (epsilon_args, representation) = Self::epsilon_args(base)?;
        if !representation.rep.is_self_dual() {
            return None;
        }

        Some(
            Self::metric_determinant(&epsilon_args, &epsilon_args)
                * Self::power(base, exponent - 2),
        )
    }

    /// Applies the epsilon-pair contraction:
    /// `epsilon(a_1,...,a_n) epsilon(b_1,...,b_n) -> det(g(a_i,b_j))`.
    ///
    /// Extra multiplicative factors are preserved around the determinant.
    fn simplify_pair_product(expr: AtomView<'_>) -> Option<Atom> {
        let AtomView::Mul(product) = expr else {
            return None;
        };
        let mut epsilons = product
            .iter()
            .enumerate()
            .filter(|(_, factor)| {
                matches!(factor, AtomView::Fun(fun) if fun.get_symbol() == *EPSILON_SYMBOL)
            })
            .peekable();

        while let Some((left_index, left)) = epsilons.next() {
            // A lone epsilon needs no representation parsing or owned arguments.
            epsilons.peek()?;
            let Some((left_args, left_rep)) = Self::epsilon_args(left) else {
                continue;
            };

            for (right_index, right) in epsilons.clone() {
                let Some((right_args, right_rep)) = Self::epsilon_args(right) else {
                    continue;
                };
                // A determinant pairs a representation with its dual. In
                // particular, two SU(3) triplet epsilons cannot contract.
                if left_args.len() != right_args.len() || left_rep.dual() != right_rep {
                    continue;
                }

                let determinant = Self::metric_determinant(&left_args, &right_args);
                let spectators =
                    Atom::mul_many(product.iter().enumerate().filter_map(|(index, factor)| {
                        (index != left_index && index != right_index).then_some(factor)
                    }));
                return Some(spectators * determinant);
            }
        }

        None
    }

    fn positive_integer(expr: AtomView<'_>) -> Option<i64> {
        let AtomView::Num(number) = expr else {
            return None;
        };
        let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
            return None;
        };

        (value > 0).then_some(value)
    }

    fn power(base: AtomView<'_>, exponent: i64) -> Atom {
        match exponent {
            0 => Atom::num(1),
            1 => base.to_owned(),
            _ => base.to_owned().pow(Atom::num(exponent)),
        }
    }

    /// Emits the determinant `det(g(left_i,right_j))`.
    fn metric_determinant(left: &[Atom], right: &[Atom]) -> Atom {
        Self::determinant(right.len(), |i, j| ETS.metric_literal(&left[i], &right[j]))
    }

    /// Sums `sign(π) Π_i entry(i, π(i))` over all permutations `π`.
    fn determinant(size: usize, entry: impl Fn(usize, usize) -> Atom) -> Atom {
        let mut permutation = (0..size).collect::<Vec<_>>();
        let mut sum = Atom::Zero;
        Self::determinant_terms(&entry, 0, &mut permutation, 1, &mut sum);
        sum
    }

    fn determinant_terms(
        entry: &impl Fn(usize, usize) -> Atom,
        position: usize,
        permutation: &mut [usize],
        sign: i64,
        sum: &mut Atom,
    ) {
        if position == permutation.len() {
            let term = (0..permutation.len()).fold(Atom::num(sign), |product, i| {
                product * entry(i, permutation[i])
            });
            *sum += term;
            return;
        }

        for i in position..permutation.len() {
            permutation.swap(position, i);
            let next_sign = if i == position { sign } else { -sign };
            Self::determinant_terms(entry, position + 1, permutation, next_sign, sum);
            permutation.swap(position, i);
        }
    }

    fn epsilon_args(epsilon: AtomView<'_>) -> Option<(Vec<Atom>, Representation<LibraryRep>)> {
        let AtomView::Fun(f) = epsilon else {
            return None;
        };
        if f.get_symbol() != *EPSILON_SYMBOL {
            return None;
        }

        let mut representation = None;
        let mut args = Vec::with_capacity(f.get_nargs());
        for arg in f.iter() {
            let rep = Representation::<LibraryRep>::try_from(arg)
                .ok()
                .or_else(|| {
                    // Compact vectors can have ordinary Symbolica heads and
                    // scalar labels, but expose exactly one bare representation.
                    let AtomView::Fun(vector) = arg else {
                        return None;
                    };
                    if vector.get_symbol().is_scalar() {
                        return None;
                    }
                    let mut reps = vector.iter().filter_map(|argument| {
                        Representation::<LibraryRep>::try_from(argument)
                            .ok()
                            .map(|rep| (argument, rep))
                    });
                    let (rep_arg, rep) = reps.next()?;
                    if reps.next().is_some()
                        || Slot::<LibraryRep, AbstractIndex>::try_from(rep_arg).is_ok()
                        || vector.iter().any(|argument| {
                            argument != rep_arg
                                && (argument.is_tensorial(StrictTensorFilter::Tagged)
                                    || argument.is_tensorial(StrictTensorFilter::ContainsReps))
                        })
                    {
                        return None;
                    }
                    Some(rep)
                })?;
            if representation.is_some_and(|previous| previous != rep) {
                return None;
            }
            representation = Some(rep);
            args.push(arg.to_owned());
        }
        let representation = representation?;
        if !representation.rep.is_self_dual()
            && representation.dim != Dimension::Concrete(args.len())
        {
            return None;
        }
        Some((args, representation))
    }
}

#[cfg(test)]
mod test {
    use insta::assert_snapshot;
    use spenso::{
        g, mink,
        network::{library::symbolic::ETS, parsing::AtomStructureExt, tags::SPENSO_TAG},
        p, q,
        structure::partial::PartialStructure,
        trace,
    };
    use symbolica::{
        atom::{Atom, AtomCore, AtomView, FunctionBuilder},
        parse,
    };
    use symbolica_utils::AtomPrintExt;

    use crate::{
        epsilon as eps, gamma, gamma5,
        tensor::{AlgebraSettings, SymbolicTensor},
        test_support::test_initialize,
    };

    use super::{EPSILON_SYMBOL, EpsilonSimplifier, contract_epsilon4, epsilon4};

    #[test]
    fn color_epsilon_dual_pair_contracts_to_positive_norm() {
        test_initialize();
        let epsilon = parse!(
            "epsilon(cof(3,i),cof(3,j),cof(3,k))",
            default_namespace = "spenso"
        );
        for (barred, expected) in [
            ("epsilon(dind(cof(3,i)),dind(cof(3,j)),dind(cof(3,k)))", "6"),
            (
                "epsilon(dind(cof(3,i)),dind(cof(3,j)),dind(cof(3,l)))",
                "2*g(cof(3,k),dind(cof(3,l)))",
            ),
            (
                "epsilon(dind(cof(3,i)),dind(cof(3,l)),dind(cof(3,m)))",
                "g(cof(3,j),dind(cof(3,l)))*g(cof(3,k),dind(cof(3,m)))-g(cof(3,j),dind(cof(3,m)))*g(cof(3,k),dind(cof(3,l)))",
            ),
        ] {
            let product = &epsilon * parse!(barred, default_namespace = "spenso");
            assert_eq!(
                product.simplify_epsilon().expand(),
                parse!(expected, default_namespace = "spenso").expand()
            );
        }
    }

    #[test]
    fn epsilon_pairs_require_dual_homogeneous_spaces() {
        test_initialize();
        for source in [
            "epsilon(cof(3,i),cof(3,j),cof(3,k))^2",
            "epsilon(dind(cof(3,i)),dind(cof(3,j)),dind(cof(3,k)))^2",
            "epsilon(cof(3,i),cof(3,j),cof(3,k))*epsilon(cof(3,l),cof(3,m),cof(3,n))",
            "epsilon(cof(3,i),cof(3,j),cof(3,k))*epsilon(mink(3,l),mink(3,m),mink(3,n))",
            "epsilon(mink(3,i),mink(3,j))*epsilon(mink(4,l),mink(4,m))",
            "epsilon(mink(D,i),mink(D,j))*epsilon(mink(E,l),mink(E,m))",
            "epsilon(cof(3,i),cof(3,j),cof(3,k))*epsilon(dind(cof(4,l)),dind(cof(4,m)),dind(cof(4,n)))",
            "epsilon(cof(4,i),cof(4,j),cof(4,k))*epsilon(dind(cof(4,l)),dind(cof(4,m)),dind(cof(4,n)))",
            "epsilon(mink(4,i),cof(3,j))*epsilon(mink(4,l),dind(cof(3,m)))",
            "epsilon(mink(4,i),A(mink(4,j),mink(4,k)))*epsilon(mink(4,l),B(mink(4,m),mink(4,n)))",
            "epsilon(mink(4,i),A(mink(4),mink(4)))*epsilon(mink(4,l),B(mink(4),mink(4)))",
        ] {
            let expression = parse!(source, default_namespace = "spenso");
            assert_eq!(expression.simplify_epsilon(), expression, "{source}");
        }
    }

    #[test]
    fn epsilon_pair_with_contracted_vectors_retains_lorentz_identity() {
        test_initialize();
        let p = p!(mink!(4));
        let q = q!(mink!(4));
        let expression = eps!(mink!(4, i), &p) * eps!(mink!(4, i), &q);
        let expected = 3 * g!(p, q);
        assert_eq!(expression.simplify_epsilon(), expected);
        for (source, expected) in [
            (
                "epsilon(mink(4,i),plain_p(mink(4)))*epsilon(mink(4,i),plain_q(mink(4)))",
                "3*g(plain_p(mink(4)),plain_q(mink(4)))",
            ),
            (
                "epsilon(mink(4,i),P(label,mink(4)))*epsilon(mink(4,i),Q(other_label,mink(4)))",
                "3*g(P(label,mink(4)),Q(other_label,mink(4)))",
            ),
        ] {
            let expression = parse!(source, default_namespace = "spenso");
            assert_eq!(
                expression.simplify_epsilon(),
                parse!(expected, default_namespace = "spenso")
            );
        }
    }

    #[test]
    fn contracted_epsilon4_reproduces_the_axial_trace() {
        // With unit vectors e_μ, e̸_0 = γ^0 and e̸_k = -γ^k, so
        // γ5 = i γ^0 γ^1 γ^2 γ^3 = -i e̸_0 e̸_1 e̸_2 e̸_3. Both traces are
        // simplified by Idenso and evaluated on contravariant components.
        let r = test_initialize();
        let spin = r.bis4.to_symbolic([]);
        let components: [[i64; 4]; 8] = [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [3, 1, -2, 5],
            [2, -1, 4, 1],
            [7, 3, 1, -2],
            [1, 5, -3, 2],
        ];
        let vectors: Vec<_> = (0..components.len())
            .map(|i| {
                let name = format!("idenso::epsilon_components::v{i}");
                FunctionBuilder::new(SPENSO_TAG.rank_one_tensor_symbol(&name))
                    .add_arg(r.mink4.to_symbolic([]))
                    .finish()
            })
            .collect();
        let lowered = |argument: AtomView| -> [Atom; 4] {
            let vector = vectors
                .iter()
                .position(|v| v.as_view() == argument)
                .unwrap();
            std::array::from_fn(|mu| {
                let sign = if mu == 0 { 1 } else { -1 };
                Atom::num(sign * components[vector][mu])
            })
        };
        let evaluate = |factors: Vec<Atom>| {
            SymbolicTensor::<PartialStructure>::infer(trace!(&spin; factors))
                .unwrap()
                .simplify_algebra(&AlgebraSettings::hep())
                .unwrap()
                .into_expression()
                .replace_map(|part, _, out| {
                    let AtomView::Fun(f) = part else {
                        return;
                    };
                    let arguments: Vec<_> = f.iter().collect();
                    if f.get_symbol() == *EPSILON_SYMBOL {
                        **out = contract_epsilon4(&std::array::from_fn(|i| lowered(arguments[i])));
                    } else if f.get_symbol() == ETS.metric {
                        let vector = vectors.iter().position(|v| v.as_view() == arguments[1]);
                        **out = lowered(arguments[0])
                            .into_iter()
                            .zip(components[vector.unwrap()])
                            .fold(Atom::Zero, |sum, (a, b)| sum + a * Atom::num(b));
                    }
                })
                .expand()
        };
        let slashed = |range: std::ops::Range<usize>| range.map(|i| gamma!(&vectors[i]));
        let axial = evaluate(std::iter::once(gamma5!()).chain(slashed(4..8)).collect());
        let explicit = evaluate(slashed(0..8).collect());
        assert!(!axial.is_zero());
        assert_eq!(axial, -Atom::i() * explicit);
        // ε(a,b,c,d) = i det[a;b;c;d] for contravariant components.
        assert_eq!(
            contract_epsilon4(&std::array::from_fn(|i| lowered(vectors[i].as_view()))),
            Atom::i()
        );
    }

    #[test]
    fn epsilon_symbol_is_antisymmetric() {
        let _ = test_initialize();
        let mu = mink!(4, mu);
        let nu = mink!(4, nu);
        let rho = mink!(4, rho);
        let sigma = mink!(4, sigma);

        assert_snapshot!(epsilon4(&nu, &mu, &rho, &sigma).to_bare_ordered_string(), @"-1*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))");
        assert!(epsilon4(&mu, &mu, &rho, &sigma).is_zero());
    }

    #[test]
    fn epsilon_macro_builds_arbitrary_rank() {
        let _ = test_initialize();
        let expr = eps!(mink!(4, a), mink!(4, b), mink!(4, c));

        assert_snapshot!(
            expr.to_bare_ordered_string(),
            @"epsilon(mink(4,a),mink(4,b),mink(4,c))"
        );
    }

    #[test]
    fn epsilon_metric_contraction_replaces_one_slot() {
        let _ = test_initialize();
        let mu = mink!(4, mu);
        let nu = mink!(4, nu);
        let rho = mink!(4, rho);
        let sigma = mink!(4, sigma);
        let p = p!(mink!(4));
        let expr = g!(&mu, &p) * epsilon4(&mu, &nu, &rho, &sigma);

        assert_snapshot!(expr.simplify_epsilon().to_bare_ordered_string(), @"-1*epsilon(mink(4,nu),mink(4,rho),mink(4,sigma),p(mink(4)))");
    }

    #[test]
    fn epsilon_pair_expands_to_metric_determinant() {
        let _ = test_initialize();
        let a = mink!(4, a);
        let b = mink!(4, b);
        let c = mink!(4, c);
        let d = mink!(4, d);
        let expr = epsilon4(&a, &b, &c, &d) * epsilon4(&a, &b, &c, &d);

        assert_snapshot!(&expr.simplify_epsilon().to_bare_ordered_string(), @"24");
    }

    #[test]
    fn epsilon_pair_with_two_indices_is_metric_determinant() {
        let a = mink!(4, a);
        let b = mink!(4, b);
        let c = mink!(4, c);
        let d = mink!(4, d);
        let left = eps!(a, b);
        let right = eps!(c, d);

        assert_snapshot!((left * right).simplify_epsilon().to_bare_ordered_string(), @"-1*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,c))+g(mink(4,a),mink(4,c))*g(mink(4,b),mink(4,d))");
    }

    #[test]
    fn epsilon_powers_reach_a_fixed_point() {
        test_initialize();
        let epsilon = eps!(mink!(4, a), mink!(4, b), mink!(4, c), mink!(4, d));
        for exponent in [2, 3, 4, 5] {
            let expression = epsilon.clone().pow(Atom::num(exponent));
            let expected = Atom::num(24_i64.pow((exponent / 2) as u32))
                * epsilon.clone().pow(Atom::num(exponent % 2));
            let result = expression.simplify_epsilon();
            assert_eq!(result, expected);
            assert_eq!(result.simplify_epsilon(), result);
        }
        for exponent in ["-1", "1/2", "n"] {
            let expression = epsilon.clone().pow(parse!(exponent));
            assert_eq!(expression.simplify_epsilon(), expression);
        }
    }

    #[test]
    fn epsilon_cleanup_preserves_lone_epsilon_and_compact_spectator_normalization() {
        test_initialize();
        let tags = &spenso::network::tags::SPENSO_TAG;
        tags.rank_one_tensor_symbol("spenso::epsilon_cleanup_p");
        tags.rank_one_tensor_symbol("spenso::epsilon_cleanup_q");
        for (source, expected) in [
            (
                "g(mink(4,a),mink(4,x))*g(mink(4,x),mink(4,y))*epsilon(mink(4,a),mink(4,b),mink(4,c),mink(4,d))",
                "epsilon(mink(4,y),mink(4,b),mink(4,c),mink(4,d))",
            ),
            (
                "g(mink(4,a),epsilon_cleanup_p(mink(4)))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
                "epsilon_cleanup_p(mink(4,a))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
            ),
            (
                "epsilon_cleanup_p(epsilon_cleanup_q(mink(4)))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
                "g(epsilon_cleanup_p(mink(4)),epsilon_cleanup_q(mink(4)))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
            ),
            (
                "g(mink(4,a),mink(4,x))*epsilon(mink(4,a),mink(4,b))*epsilon(mink(4,x),mink(4,c))",
                "3*g(mink(4,b),mink(4,c))",
            ),
        ] {
            let input = parse!(source, default_namespace = "spenso");
            let expected = parse!(expected, default_namespace = "spenso");
            let result = input.simplify_epsilon();
            assert_eq!(result, expected, "{source}");
            assert_eq!(result.simplify_epsilon(), result, "{source}");
        }
    }

    #[test]
    fn distinct_epsilon_pairs_simplify_without_expanding_surrounding_factors() {
        test_initialize();
        let pair = parse!(
            "epsilon(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d))",
            default_namespace = "spenso"
        );
        assert!(!pair.has_repeated_explicit_indices());
        // Preserve the determinant's own extracted antisymmetry sign rather
        // than distributing it through a surrounding factorized expression.
        let determinant = pair.simplify_epsilon();
        assert_ne!(determinant, pair);
        let spectator = parse!(
            "epsilon(cof(3,i),cof(3,j),cof(3,k))*(x+y)^5",
            default_namespace = "spenso"
        );
        let input = &spectator * &pair;
        assert_eq!(input.simplify_epsilon(), &spectator * &determinant);

        let scalar = parse!("z", default_namespace = "spenso");
        let input = (scalar.clone() + pair).pow(Atom::num(2));
        let expected = (scalar + determinant).pow(Atom::num(2));
        assert_eq!(input.simplify_epsilon(), expected);
        assert_eq!(expected.simplify_epsilon(), expected);
    }
}
