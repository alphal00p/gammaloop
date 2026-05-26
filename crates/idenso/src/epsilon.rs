use std::sync::LazyLock;

use spenso::g;
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
};

use crate::shorthands::schoonschip::Schoonschip;

/// Symbolica-level Levi-Civita symbol. The `Antisymmetric` attribute lets
/// Symbolica canonicalize argument order and annihilate repeated arguments.
pub static EPSILON_SYMBOL: LazyLock<Symbol> =
    LazyLock::new(|| spenso::tensor_symbol!("spenso::epsilon"; Antisymmetric));

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

/// Simplifies epsilon-metric contractions and epsilon-pair contractions.
pub trait EpsilonSimplifier {
    fn simplify_epsilon(&self) -> Atom;
}

impl EpsilonSimplifier for Atom {
    fn simplify_epsilon(&self) -> Atom {
        self.as_view().simplify_epsilon()
    }
}

impl EpsilonSimplifier for AtomView<'_> {
    fn simplify_epsilon(&self) -> Atom {
        EpsilonSimplifierPass::run(*self)
    }
}

struct EpsilonSimplifierPass;

impl EpsilonSimplifierPass {
    /// Runs the epsilon identities to a fixed point, then lets metric
    /// simplification consume any determinants that collapsed to traces.
    fn run(expr: AtomView<'_>) -> Atom {
        // Force Symbolica to register epsilon as antisymmetric before any
        // builders below create new epsilon nodes.
        let _ = *EPSILON_SYMBOL;
        if !Self::contains_epsilon(expr) {
            return expr.to_owned();
        }

        let mut current = expr.schoonschip();

        loop {
            let next = current.replace_map(|term, _context, out| {
                if let Some(rewritten) = Self::simplify_power(term.schoonschip().as_view())
                    .or_else(|| Self::simplify_pair_product(term))
                {
                    **out = rewritten;
                }
            });
            // .collect_metrics()
            // .schoonschip();

            if next == current {
                return next;
            }

            current = next;
        }
    }

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
        let epsilon_args = Self::epsilon_args(base)?;
        let exponent = Self::positive_integer(exponent)?;
        if exponent < 2 {
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
        let factors = Self::multiplicative_factors(expr);
        if factors.len() < 2 {
            return None;
        }

        for (left_index, left) in factors.iter().enumerate() {
            let Some(left_args) = Self::epsilon_args(*left) else {
                continue;
            };

            for (right_index, right) in factors.iter().enumerate().skip(left_index + 1) {
                let Some(right_args) = Self::epsilon_args(*right) else {
                    continue;
                };
                if left_args.len() != right_args.len() || left_args.is_empty() {
                    continue;
                }

                let determinant = Self::metric_determinant(&left_args, &right_args);
                let mut excluded = vec![false; factors.len()];
                excluded[left_index] = true;
                excluded[right_index] = true;
                return Some(Self::product_excluding(&factors, &excluded) * determinant);
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

    /// Emits the determinant `det(g(left_i,right_j))` by summing over all
    /// signed permutations of the right indices.
    fn metric_determinant(left: &[Atom], right: &[Atom]) -> Atom {
        let mut permutation = (0..right.len()).collect::<Vec<_>>();
        let mut sum = Atom::Zero;
        Self::determinant_terms(left, right, 0, &mut permutation, 1, &mut sum);
        sum
    }

    fn determinant_terms(
        left: &[Atom],
        right: &[Atom],
        position: usize,
        permutation: &mut [usize],
        sign: i64,
        sum: &mut Atom,
    ) {
        if position == permutation.len() {
            let term = left
                .iter()
                .enumerate()
                .fold(Atom::num(sign), |product, (i, arg)| {
                    product * g!(arg.clone(), right[permutation[i]].clone())
                });
            *sum += term;
            return;
        }

        for i in position..permutation.len() {
            permutation.swap(position, i);
            let next_sign = if i == position { sign } else { -sign };
            Self::determinant_terms(left, right, position + 1, permutation, next_sign, sum);
            permutation.swap(position, i);
        }
    }

    fn epsilon_args(epsilon: AtomView<'_>) -> Option<Vec<Atom>> {
        let AtomView::Fun(f) = epsilon else {
            return None;
        };
        if f.get_symbol() != *EPSILON_SYMBOL {
            return None;
        }

        Some(f.iter().map(|arg| arg.to_owned()).collect())
    }

    fn multiplicative_factors<'a>(expr: AtomView<'a>) -> Vec<AtomView<'a>> {
        match expr {
            AtomView::Mul(mul) => mul.iter().collect(),
            _ => vec![expr],
        }
    }

    fn product_excluding(factors: &[AtomView<'_>], excluded: &[bool]) -> Atom {
        factors
            .iter()
            .enumerate()
            .filter(|(index, _)| !excluded[*index])
            .fold(Atom::num(1), |product, (_, factor)| {
                product * factor.to_owned()
            })
    }
}

#[cfg(test)]
mod test {
    use insta::assert_snapshot;
    use spenso::{g, mink, p, shadowing::symbolica_utils::AtomCoreExt};

    use crate::{epsilon as eps, test_support::test_initialize};

    use super::{EpsilonSimplifier, epsilon4};

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

        assert_snapshot!(expr.simplify_epsilon().to_bare_ordered_string(), @"(g(mink(4,a),mink(4,b)))^2*(g(mink(4,c),mink(4,d)))^2+(g(mink(4,a),mink(4,b)))^2*-1*g(mink(4,c),mink(4,c))*g(mink(4,d),mink(4,d))+(g(mink(4,a),mink(4,c)))^2*(g(mink(4,b),mink(4,d)))^2+(g(mink(4,a),mink(4,c)))^2*-1*g(mink(4,b),mink(4,b))*g(mink(4,d),mink(4,d))+(g(mink(4,a),mink(4,d)))^2*(g(mink(4,b),mink(4,c)))^2+(g(mink(4,a),mink(4,d)))^2*-1*g(mink(4,b),mink(4,b))*g(mink(4,c),mink(4,c))+(g(mink(4,b),mink(4,c)))^2*-1*g(mink(4,a),mink(4,a))*g(mink(4,d),mink(4,d))+(g(mink(4,b),mink(4,d)))^2*-1*g(mink(4,a),mink(4,a))*g(mink(4,c),mink(4,c))+(g(mink(4,c),mink(4,d)))^2*-1*g(mink(4,a),mink(4,a))*g(mink(4,b),mink(4,b))+-2*g(mink(4,a),mink(4,b))*g(mink(4,a),mink(4,c))*g(mink(4,b),mink(4,d))*g(mink(4,c),mink(4,d))+-2*g(mink(4,a),mink(4,b))*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,c))*g(mink(4,c),mink(4,d))+-2*g(mink(4,a),mink(4,c))*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,c))*g(mink(4,b),mink(4,d))+2*g(mink(4,a),mink(4,a))*g(mink(4,b),mink(4,c))*g(mink(4,b),mink(4,d))*g(mink(4,c),mink(4,d))+2*g(mink(4,a),mink(4,b))*g(mink(4,a),mink(4,c))*g(mink(4,b),mink(4,c))*g(mink(4,d),mink(4,d))+2*g(mink(4,a),mink(4,b))*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,d))*g(mink(4,c),mink(4,c))+2*g(mink(4,a),mink(4,c))*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,b))*g(mink(4,c),mink(4,d))+g(mink(4,a),mink(4,a))*g(mink(4,b),mink(4,b))*g(mink(4,c),mink(4,c))*g(mink(4,d),mink(4,d))");
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
}
