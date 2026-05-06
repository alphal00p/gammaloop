use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    structure::abstract_index::AIND_SYMBOLS,
};

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};

use super::utils::multiplicative_factors;

pub(super) fn simplify_chain_like_metric_products(expr: AtomView<'_>) -> Atom {
    let mut current = expr.to_owned();

    loop {
        let next = current.replace_map(|term, _context, out| {
            if let Some(rewritten) = simplify_chain_like_metric_product(term) {
                **out = rewritten;
            }
        });

        if next == current {
            return next;
        }

        current = next;
    }
}

fn simplify_chain_like_metric_product(expr: AtomView<'_>) -> Option<Atom> {
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

            for (target_index, target) in factors.iter().enumerate() {
                if target_index == metric_index {
                    continue;
                }

                let Some(rewritten_target) =
                    replace_first_in_chain_like(target.as_view(), old, new)
                else {
                    continue;
                };

                return Some(product_with_replaced_factor(
                    &factors,
                    metric_index,
                    target_index,
                    rewritten_target,
                ));
            }
        }
    }

    None
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

pub(super) fn is_slot(atom: &Atom) -> bool {
    is_slot_view(atom.as_view())
}

fn is_slot_view(atom: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = atom else {
        return false;
    };

    if f.get_symbol().has_tag(&T.representation) {
        return f.get_nargs() == 2;
    }

    // Downstairs indices wrap the actual slot as `dind(slot)`.
    f.get_symbol() == AIND_SYMBOLS.dind
        && f.get_nargs() == 1
        && f.iter().next().is_some_and(is_slot_view)
}

fn replace_first_in_chain_like(chain_like: AtomView<'_>, old: &Atom, new: &Atom) -> Option<Atom> {
    let AtomView::Fun(f) = chain_like else {
        return None;
    };

    if f.get_symbol() != T.chain && f.get_symbol() != T.trace {
        return None;
    }

    let mut replaced = false;
    let mut builder = FunctionBuilder::new(f.get_symbol());
    for arg in f.iter() {
        if !replaced && let Some(new_arg) = replace_first_exact_atom(arg, old, new) {
            builder = builder.add_arg(new_arg);
            replaced = true;
        } else {
            builder = builder.add_arg(arg);
        }
    }

    replaced.then(|| builder.finish())
}

fn replace_first_exact_atom(expr: AtomView<'_>, old: &Atom, new: &Atom) -> Option<Atom> {
    if expr == old.as_view() {
        return Some(new.clone());
    }

    let AtomView::Fun(f) = expr else {
        return None;
    };

    let mut replaced = false;
    let mut builder = FunctionBuilder::new(f.get_symbol());
    for arg in f.iter() {
        if !replaced && let Some(new_arg) = replace_first_exact_atom(arg, old, new) {
            builder = builder.add_arg(new_arg);
            replaced = true;
        } else {
            builder = builder.add_arg(arg);
        }
    }

    replaced.then(|| builder.finish())
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
