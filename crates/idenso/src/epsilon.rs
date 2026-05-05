use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    structure::abstract_index::AIND_SYMBOLS,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::CoefficientView,
    function, symbol,
};

use crate::metric::MetricSimplifier;

/// Symbolica-level Levi-Civita symbol. The `Antisymmetric` attribute lets
/// Symbolica canonicalize argument order and annihilate repeated arguments.
pub static EPSILON_SYMBOL: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("spenso::epsilon"; Antisymmetric));

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

/// Builds a four-dimensional Levi-Civita tensor.
pub fn epsilon4(mu: &Atom, nu: &Atom, rho: &Atom, sigma: &Atom) -> Atom {
    epsilon([mu, nu, rho, sigma])
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
        simplify_epsilon_impl(*self)
    }
}

fn simplify_epsilon_impl(expr: AtomView<'_>) -> Atom {
    // Force Symbolica to register epsilon as antisymmetric before any builders
    // below create new epsilon nodes.
    let _ = *EPSILON_SYMBOL;
    let mut current = expr.to_owned();

    loop {
        let next = current
            .replace_map(|term, _context, out| {
                if let Some(rewritten) = simplify_epsilon_metric_product(term)
                    .or_else(|| simplify_epsilon_power(term))
                    .or_else(|| simplify_epsilon_pair_product(term))
                {
                    **out = rewritten;
                }
            })
            .expand()
            .simplify_metrics();

        if next == current {
            return next;
        }

        current = next;
    }
}

fn simplify_epsilon_metric_product(expr: AtomView<'_>) -> Option<Atom> {
    let factors = multiplicative_factors(expr);
    if factors.len() < 2 {
        return None;
    }

    for (metric_index, metric) in factors.iter().enumerate() {
        let Some((left, right)) = metric_args(metric.as_view()) else {
            continue;
        };

        for (old, new) in [(&left, &right), (&right, &left)] {
            if !is_slot(old) {
                continue;
            }

            for (epsilon_index, epsilon) in factors.iter().enumerate() {
                if epsilon_index == metric_index {
                    continue;
                }

                let Some(rewritten_epsilon) =
                    replace_first_epsilon_arg(epsilon.as_view(), old, new)
                else {
                    continue;
                };

                return Some(product_with_replaced_factor(
                    &factors,
                    metric_index,
                    epsilon_index,
                    rewritten_epsilon,
                ));
            }
        }
    }

    None
}

fn simplify_epsilon_power(expr: AtomView<'_>) -> Option<Atom> {
    let AtomView::Pow(pow) = expr else {
        return None;
    };
    let (base, exponent) = pow.get_base_exp();
    let epsilon_args = epsilon_args(base)?;
    let exponent = positive_integer(exponent)?;
    if exponent < 2 {
        return None;
    }

    // Antisymmetric products often normalize as epsilon(...)^2. Split off one
    // epsilon pair and leave any remaining power untouched.
    Some(metric_determinant(&epsilon_args, &epsilon_args) * power(base, exponent - 2))
}

fn simplify_epsilon_pair_product(expr: AtomView<'_>) -> Option<Atom> {
    let factors = multiplicative_factors(expr);
    if factors.len() < 2 {
        return None;
    }

    for (left_index, left) in factors.iter().enumerate() {
        let Some(left_args) = epsilon_args(left.as_view()) else {
            continue;
        };

        for (right_index, right) in factors.iter().enumerate().skip(left_index + 1) {
            let Some(right_args) = epsilon_args(right.as_view()) else {
                continue;
            };
            if left_args.len() != right_args.len() || left_args.is_empty() {
                continue;
            }

            let determinant = metric_determinant(&left_args, &right_args);
            let mut excluded = vec![false; factors.len()];
            excluded[left_index] = true;
            excluded[right_index] = true;
            return Some(product_excluding(&factors, &excluded) * determinant);
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

fn metric_determinant(left: &[Atom], right: &[Atom]) -> Atom {
    let mut permutation = (0..right.len()).collect::<Vec<_>>();
    let mut sum = Atom::Zero;
    determinant_terms(left, right, 0, &mut permutation, 1, &mut sum);
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
                product * function!(ETS.metric, arg.clone(), right[permutation[i]].clone())
            });
        *sum += term;
        return;
    }

    for i in position..permutation.len() {
        permutation.swap(position, i);
        let next_sign = if i == position { sign } else { -sign };
        determinant_terms(left, right, position + 1, permutation, next_sign, sum);
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

fn replace_first_epsilon_arg(epsilon: AtomView<'_>, old: &Atom, new: &Atom) -> Option<Atom> {
    let AtomView::Fun(f) = epsilon else {
        return None;
    };
    if f.get_symbol() != *EPSILON_SYMBOL {
        return None;
    }

    let mut replaced = false;
    let mut builder = FunctionBuilder::new(*EPSILON_SYMBOL);
    for arg in f.iter() {
        if !replaced && arg == old.as_view() {
            builder = builder.add_arg(new.clone());
            replaced = true;
        } else {
            builder = builder.add_arg(arg);
        }
    }

    replaced.then(|| builder.finish())
}

fn metric_args(metric: AtomView<'_>) -> Option<(Atom, Atom)> {
    let AtomView::Fun(f) = metric else {
        return None;
    };
    if f.get_symbol() != ETS.metric || f.get_nargs() != 2 {
        return None;
    }

    let args = f.iter().map(|arg| arg.to_owned()).collect::<Vec<_>>();
    Some((args[0].clone(), args[1].clone()))
}

fn is_slot(atom: &Atom) -> bool {
    is_slot_view(atom.as_view())
}

fn is_slot_view(atom: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = atom else {
        return false;
    };

    if f.get_symbol().has_tag(&T.representation) {
        return f.get_nargs() == 2;
    }

    f.get_symbol() == AIND_SYMBOLS.dind
        && f.get_nargs() == 1
        && f.iter().next().is_some_and(is_slot_view)
}

fn multiplicative_factors(expr: AtomView<'_>) -> Vec<Atom> {
    match expr {
        AtomView::Mul(mul) => mul.iter().map(|factor| factor.to_owned()).collect(),
        _ => vec![expr.to_owned()],
    }
}

fn product_with_replaced_factor(
    factors: &[Atom],
    removed_index: usize,
    replaced_index: usize,
    replacement: Atom,
) -> Atom {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != removed_index)
        .fold(Atom::num(1), |product, (index, factor)| {
            if index == replaced_index {
                product * &replacement
            } else {
                product * factor
            }
        })
}

fn product_excluding(factors: &[Atom], excluded: &[bool]) -> Atom {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| !excluded[*index])
        .fold(Atom::num(1), |product, (_, factor)| product * factor)
}

#[cfg(test)]
fn mink4(index: impl Into<Atom>) -> Atom {
    use spenso::structure::representation::{Minkowski, RepName};

    Minkowski {}.to_symbolic([Atom::num(4), index.into()])
}

#[cfg(test)]
mod test {
    use insta::assert_snapshot;
    use spenso::{
        network::library::symbolic::ETS,
        shadowing::symbolica_utils::AtomCoreExt,
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::{atom::Atom, atom::FunctionBuilder, function, symbol};

    use super::{EpsilonSimplifier, epsilon4, mink4};

    #[test]
    fn epsilon_symbol_is_antisymmetric() {
        let mu = mink4(symbol!("mu"));
        let nu = mink4(symbol!("nu"));
        let rho = mink4(symbol!("rho"));
        let sigma = mink4(symbol!("sigma"));

        assert_snapshot!(epsilon4(&nu, &mu, &rho, &sigma).to_bare_ordered_string(), @"-1*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))");
        assert!(epsilon4(&mu, &mu, &rho, &sigma).is_zero());
    }

    #[test]
    fn epsilon_metric_contraction_replaces_one_slot() {
        let mu = mink4(symbol!("mu"));
        let nu = mink4(symbol!("nu"));
        let rho = mink4(symbol!("rho"));
        let sigma = mink4(symbol!("sigma"));
        let p = function!(
            symbol!("spenso::p"),
            Minkowski {}.to_symbolic([Atom::num(4)])
        );
        let expr = function!(ETS.metric, &mu, &p) * epsilon4(&mu, &nu, &rho, &sigma);

        assert_snapshot!(expr.simplify_epsilon().to_bare_ordered_string(), @"-1*epsilon(mink(4,nu),mink(4,rho),mink(4,sigma),p(mink(4)))");
    }

    #[test]
    fn epsilon_pair_expands_to_metric_determinant() {
        let a = mink4(symbol!("a"));
        let b = mink4(symbol!("b"));
        let c = mink4(symbol!("c"));
        let d = mink4(symbol!("d"));
        let expr = epsilon4(&a, &b, &c, &d) * epsilon4(&a, &b, &c, &d);

        assert_snapshot!(expr.simplify_epsilon().to_bare_ordered_string(), @"24");
    }

    #[test]
    fn epsilon_pair_with_two_indices_is_metric_determinant() {
        let a = mink4(symbol!("a"));
        let b = mink4(symbol!("b"));
        let c = mink4(symbol!("c"));
        let d = mink4(symbol!("d"));
        let left = FunctionBuilder::new(*super::EPSILON_SYMBOL)
            .add_arg(a)
            .add_arg(b)
            .finish();
        let right = FunctionBuilder::new(*super::EPSILON_SYMBOL)
            .add_arg(c)
            .add_arg(d)
            .finish();

        assert_snapshot!((left * right).simplify_epsilon().to_bare_ordered_string(), @"-1*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,c))+g(mink(4,a),mink(4,c))*g(mink(4,b),mink(4,d))");
    }
}
