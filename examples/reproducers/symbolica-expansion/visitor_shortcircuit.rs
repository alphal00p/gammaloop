//! Show the difference between subtree pruning and terminating an Atom search.

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};

fn is_target(value: AtomView<'_>, target: Symbol) -> bool {
    matches!(value, AtomView::Var(variable) if variable.get_symbol() == target)
}

fn is_opaque(value: AtomView<'_>, opaque: Symbol) -> bool {
    matches!(value, AtomView::Fun(function) if function.get_symbol() == opaque)
}

fn visitor_search(value: AtomView<'_>, target: Symbol, opaque: Symbol) -> (bool, usize) {
    let mut found = false;
    let mut visits = 0;
    value.visitor(&mut |node| {
        visits += 1;
        if found {
            return false;
        }
        found = is_target(node, target);
        !found && !is_opaque(node, opaque)
    });
    (found, visits)
}

fn shortcircuit_search(
    value: AtomView<'_>,
    target: Symbol,
    opaque: Symbol,
    visits: &mut usize,
) -> bool {
    *visits += 1;
    if is_target(value, target) {
        return true;
    }
    if is_opaque(value, opaque) {
        return false;
    }
    match value {
        AtomView::Fun(function) => function
            .iter()
            .any(|child| shortcircuit_search(child, target, opaque, visits)),
        AtomView::Mul(product) => product
            .iter()
            .any(|child| shortcircuit_search(child, target, opaque, visits)),
        AtomView::Add(sum) => sum
            .iter()
            .any(|child| shortcircuit_search(child, target, opaque, visits)),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            shortcircuit_search(base, target, opaque, visits)
                || shortcircuit_search(exponent, target, opaque, visits)
        }
        AtomView::Num(_) | AtomView::Var(_) => false,
    }
}

fn main() {
    let (scope, leaf, opaque, target, x) = symbolica::symbol!(
        "visitor_probe::scope",
        "visitor_probe::leaf",
        "visitor_probe::opaque",
        "visitor_probe::target",
        "visitor_probe::x"
    );
    let target_atom = Atom::var(target);
    let leaf_atom = FunctionBuilder::new(leaf).add_arg(Atom::var(x)).finish();
    let hidden = FunctionBuilder::new(opaque).add_arg(&target_atom).finish();

    for width in [8, 256, 8192] {
        for position in ["first", "last", "absent", "opaque"] {
            let mut args = vec![leaf_atom.as_view(); width];
            match position {
                "first" => args.insert(0, target_atom.as_view()),
                "last" => args.push(target_atom.as_view()),
                "opaque" => args.insert(0, hidden.as_view()),
                _ => {}
            }
            let source = FunctionBuilder::new(scope).add_args(args).finish();
            let (found, visitor_visits) = visitor_search(source.as_view(), target, opaque);
            let mut short_visits = 0;
            let short_found =
                shortcircuit_search(source.as_view(), target, opaque, &mut short_visits);
            assert_eq!(found, matches!(position, "first" | "last"));
            assert_eq!(short_found, found);
            if position == "first" {
                assert_eq!(visitor_visits, width + 2);
                assert_eq!(short_visits, 2);
            } else {
                assert_eq!(visitor_visits, short_visits);
            }
            println!(
                "{{\"width\":{width},\"position\":\"{position}\",\"found\":{found},\"visitor_visits\":{visitor_visits},\"shortcircuit_visits\":{short_visits}}}"
            );
        }
    }
}
