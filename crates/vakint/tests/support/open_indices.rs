//! Shared helpers of the D6 open-index tests (helper VH, lane
//! d6-open-indices): close the open indices of a result with a closing tensor
//! (metric cycles give `d = 4 - 2ε`, vector chains `dot`s), read Laurent
//! coefficients, and compare order by order.
#![allow(dead_code)]

use std::collections::BTreeMap;

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use vakint::vakint_symbol;

/// The open index carried by a vector or metric argument: a raw index, or
/// the index of a Spenso Minkowski slot `mink(D, index)`.
pub fn index_of(argument: AtomView<'_>) -> Atom {
    if let AtomView::Fun(slot) = argument
        && slot.get_symbol().get_stripped_name() == "mink"
        && slot.get_nargs() == 2
    {
        return slot.iter().nth(1).unwrap().to_owned();
    }
    argument.to_owned()
}

/// Contract every repeated index of `expression` (a sum of products of
/// scalars, metrics `g(a, b)` and vectors `p(i, a)`), in `d = 4 - 2ε`: a metric
/// cycle gives `d`, a chain between two vectors their `dot`. Every index must
/// occur exactly twice in a product.
pub fn close_indices(expression: &Atom, epsilon: Symbol) -> Atom {
    let dimension = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let expanded = expression.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
        _ => vec![expanded.clone()],
    };
    let mut total = Atom::Zero;
    for term in terms {
        let mut factors = Vec::new();
        let push = |factor: AtomView<'_>, factors: &mut Vec<Atom>| {
            if let AtomView::Pow(power) = factor {
                let (base, exponent) = power.get_base_exp();
                if let Ok(n) = i64::try_from(exponent)
                    && n > 0
                    && is_indexed(base)
                {
                    for _ in 0..n {
                        factors.push(base.to_owned());
                    }
                    return;
                }
            }
            factors.push(factor.to_owned());
        };
        match term.as_view() {
            AtomView::Mul(product) => product.iter().for_each(|f| push(f, &mut factors)),
            other => push(other, &mut factors),
        }
        let mut scalar = Atom::num(1);
        // index -> endpoints: (factor position, which end)
        let mut metrics = Vec::new();
        let mut vectors = Vec::new();
        for factor in factors {
            match factor.as_view() {
                AtomView::Fun(function)
                    if function.get_symbol().get_stripped_name() == "g"
                        && function.get_nargs() == 2 =>
                {
                    let mut arguments = function.iter();
                    metrics.push((
                        index_of(arguments.next().unwrap()),
                        index_of(arguments.next().unwrap()),
                    ));
                }
                AtomView::Fun(function)
                    if function.get_symbol().get_stripped_name() == "p"
                        && function.get_nargs() == 2 =>
                {
                    let mut arguments = function.iter();
                    let id = arguments.next().unwrap().to_owned();
                    vectors.push((id, index_of(arguments.next().unwrap())));
                }
                _ => scalar = scalar * factor,
            }
        }
        // Walk chains from every vector, then the remaining metric cycles.
        let mut used = vec![false; metrics.len()];
        let mut vector_used = vec![false; vectors.len()];
        for start in 0..vectors.len() {
            if vector_used[start] {
                continue;
            }
            vector_used[start] = true;
            let mut index = vectors[start].1.clone();
            loop {
                if let Some(other) =
                    (0..vectors.len()).find(|&v| !vector_used[v] && vectors[v].1 == index)
                {
                    vector_used[other] = true;
                    let head = |id: &Atom| {
                        FunctionBuilder::new(vakint_symbol!("p"))
                            .add_arg(id)
                            .finish()
                    };
                    scalar = scalar
                        * FunctionBuilder::new(vakint_symbol!("dot"))
                            .add_arg(head(&vectors[start].0))
                            .add_arg(head(&vectors[other].0))
                            .finish();
                    break;
                }
                let Some(next) = (0..metrics.len())
                    .find(|&m| !used[m] && (metrics[m].0 == index || metrics[m].1 == index))
                else {
                    panic!("index {index} occurs once in {term}");
                };
                used[next] = true;
                index = if metrics[next].0 == index {
                    metrics[next].1.clone()
                } else {
                    metrics[next].0.clone()
                };
            }
        }
        for start in 0..metrics.len() {
            if used[start] {
                continue;
            }
            used[start] = true;
            let first = metrics[start].0.clone();
            let mut index = metrics[start].1.clone();
            while index != first {
                let next = (0..metrics.len())
                    .find(|&m| !used[m] && (metrics[m].0 == index || metrics[m].1 == index))
                    .unwrap_or_else(|| panic!("open metric chain in {term}"));
                used[next] = true;
                index = if metrics[next].0 == index {
                    metrics[next].1.clone()
                } else {
                    metrics[next].0.clone()
                };
            }
            scalar = scalar * &dimension;
        }
        total = total + scalar;
    }
    total
}

pub fn is_indexed(view: AtomView<'_>) -> bool {
    matches!(view, AtomView::Fun(function)
        if (function.get_symbol().get_stripped_name() == "g"
            || function.get_symbol().get_stripped_name() == "p")
            && function.get_nargs() == 2)
}

/// The coefficients of the Laurent series `value` in `ε`, by power.
pub fn epsilon_coefficients(value: &Atom, epsilon: Symbol) -> BTreeMap<i64, Atom> {
    let expanded = value.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
        view if view.is_zero() => vec![],
        _ => vec![expanded.clone()],
    };
    let mut coefficients = BTreeMap::<i64, Atom>::new();
    for term in terms {
        let factors = match term.as_view() {
            AtomView::Mul(product) => product.iter().map(|f| f.to_owned()).collect::<Vec<_>>(),
            _ => vec![term.clone()],
        };
        let mut order = 0;
        let mut rest = Atom::num(1);
        for factor in factors {
            match factor.as_view() {
                AtomView::Var(variable) if variable.get_symbol() == epsilon => order += 1,
                AtomView::Pow(power)
                    if matches!(power.get_base_exp().0, AtomView::Var(v) if v.get_symbol() == epsilon) =>
                {
                    order += i64::try_from(power.get_base_exp().1).unwrap();
                }
                _ => rest = rest * factor,
            }
        }
        let entry = coefficients.entry(order).or_insert(Atom::Zero);
        *entry = &*entry + rest;
    }
    coefficients
        .into_iter()
        .map(|(order, value)| (order, value.expand()))
        .filter(|(_, value)| !value.is_zero())
        .collect()
}

/// `open` contracted with `closing` equals `scalar` through the last order
/// of `scalar` (or, for an exact zero, through the last complete order of the
/// contraction).
pub fn assert_closed_case(open: &Atom, closing: &Atom, scalar: &Atom, epsilon: Symbol, what: &str) {
    let closed = close_indices(&(open * closing), epsilon);
    let last = last_order(scalar, &closed, epsilon);
    assert_series_through(&closed, scalar, last, what, epsilon);
}

/// The last order to compare: that of `reference`, or, when it is an exact
/// zero, the last order `closed` determines exactly (all but its last).
pub fn last_order(reference: &Atom, closed: &Atom, epsilon: Symbol) -> i64 {
    match epsilon_coefficients(reference, epsilon).keys().next_back() {
        Some(last) => *last,
        None => epsilon_coefficients(closed, epsilon)
            .keys()
            .next_back()
            .map_or(0, |last| last - 1),
    }
}

/// `left` and `right` agree order by order in `ε` through `last`.
pub fn assert_series_through(left: &Atom, right: &Atom, last: i64, what: &str, epsilon: Symbol) {
    let left = epsilon_coefficients(left, epsilon);
    let right = epsilon_coefficients(right, epsilon);
    let first = left.keys().chain(right.keys()).copied().min().unwrap_or(0);
    for order in first..=last {
        let l = left.get(&order).cloned().unwrap_or(Atom::Zero);
        let r = right.get(&order).cloned().unwrap_or(Atom::Zero);
        let difference = (&l - &r).expand();
        assert!(
            difference.is_zero(),
            "{what}: order ε^{order}: closed open-index result\n  {l}\nreference (scalar path)\n  {r}\ndifference {difference}"
        );
    }
}
