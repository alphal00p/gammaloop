//! Linearisation of momentum sums after a substitution
//! (`GAMMALOOP_LINEARIZE_MOMENTA`, on when unset: production since r5g,
//! results/numerator-routing; `0`, `off`, `false` or empty: inert).
//!
//! Routing substitutes every edge momentum `Q(e, ...)` by its basis sum. A
//! Feynman-rule combination `Q_a - Q_b` then becomes a nested sum such as
//! `Q_5-Q_8+Q_11-(-Q_9-Q_11)` that nothing collapses. Here every sum that
//! directly contains a numeric multiple of another sum is flattened with
//! `expand_num` (numbers distributed, never products of sums), so the example
//! becomes `Q_5+2Q_11-Q_8+Q_9`. Sums without such a child, e.g. the factors of
//! `c*(Q_1+Q_2)*(Q_3+Q_4)`, are left as they are: nothing is expanded or
//! rewritten that is already linear. The value is unchanged.

use std::sync::LazyLock;

use symbolica::atom::{Atom, AtomCore, AtomView};

/// Whether `GAMMALOOP_LINEARIZE_MOMENTA` is on (unset, `1`, `on`, `true`).
pub(crate) fn enabled() -> bool {
    static ON: LazyLock<bool> =
        LazyLock::new(
            || match std::env::var("GAMMALOOP_LINEARIZE_MOMENTA").as_deref() {
                Err(std::env::VarError::NotPresent) | Ok("1" | "on" | "true") => true,
                Err(std::env::VarError::NotUnicode(_)) | Ok("" | "0" | "off" | "false") => false,
                Ok(other) => {
                    panic!("GAMMALOOP_LINEARIZE_MOMENTA={other:?}: expected 0/off or 1/on")
                }
            },
        );
    *ON
}

/// `c * (sum)` with a numeric `c` (the shape a substituted sum takes inside
/// another sum).
fn is_scaled_sum(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Mul(product) => {
            product.get_nargs() == 2
                && product.iter().any(|f| matches!(f, AtomView::Num(_)))
                && product.iter().any(|f| matches!(f, AtomView::Add(_)))
        }
        _ => false,
    }
}

/// A sum with a nested (scaled) sum among its terms.
fn is_nested_sum(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Add(sum) => sum
            .iter()
            .any(|t| is_scaled_sum(t) || matches!(t, AtomView::Add(_))),
        _ => false,
    }
}

/// Whether `view` contains a nested sum anywhere (also inside function
/// arguments).
pub(crate) fn has_nested_sum(view: AtomView<'_>) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        if found {
            return false;
        }
        if is_nested_sum(part) {
            found = true;
            return false;
        }
        true
    });
    found
}

/// Flatten the nested sums of `view`; `None` when there is none.
pub(crate) fn linearize_view(view: AtomView<'_>) -> Option<Atom> {
    if !has_nested_sum(view) {
        return None;
    }
    Some(view.replace_map(|term, _, out| {
        if is_nested_sum(term) {
            **out = term.expand_num();
        }
    }))
}

/// [`linearize_view`] when [`enabled`], by value.
pub(crate) fn linearize(atom: Atom) -> Atom {
    if !enabled() {
        return atom;
    }
    linearize_view(atom.as_view()).unwrap_or(atom)
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse;

    #[test]
    fn flattens_nested_sums_only() {
        let nested = parse!("q5-q8+q11-(-q9-q11)");
        assert!(has_nested_sum(nested.as_view()));
        let flat = linearize_view(nested.as_view()).unwrap();
        assert_eq!(flat, parse!("q5+2*q11-q8+q9"));
        // products of sums stay products
        let product = parse!("3*(a+b)*(c+d)*f(x)");
        assert!(linearize_view(product.as_view()).is_none());
        // a nested sum inside a product of sums: only that sum is flattened
        let mixed = parse!("(a-(b+c))*(c+d)");
        let out = linearize_view(mixed.as_view()).unwrap();
        assert_eq!(out, parse!("(a-b-c)*(c+d)"));
        // inside a function argument
        let inside = parse!("f(a-2*(b-c),e)");
        assert_eq!(
            linearize_view(inside.as_view()).unwrap(),
            parse!("f(a-2*b+2*c,e)")
        );
    }
}
