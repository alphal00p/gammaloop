use symbolica::atom::{Atom, AtomCore, AtomView};

pub(super) const TRACE_SCHOONSCHIP: bool = false;
pub(super) fn is_sum(expr: &Atom) -> bool {
    matches!(expr.as_view(), AtomView::Add(_))
}

pub(super) fn expression_size(expr: &Atom) -> (usize, usize) {
    (expr.as_view().get_byte_size(), expr.nterms())
}

pub(super) fn trace_sum_contractions() -> bool {
    std::env::var_os("IDENSO_TRACE_SUM_CONTRACTIONS").is_some()
}

pub(super) fn disable_direct_sum_contractions() -> bool {
    std::env::var_os("IDENSO_DISABLE_DIRECT_SUM_CONTRACTIONS").is_some()
}

pub(super) fn trace_contraction_ordering() -> bool {
    std::env::var_os("IDENSO_TRACE_CONTRACTION_ORDERING").is_some()
}

pub(super) fn trace_direct_sum_terms() -> bool {
    std::env::var_os("IDENSO_TRACE_DIRECT_SUM_TERMS").is_some()
}

pub(super) fn trace_direct_sum_term_expressions() -> bool {
    std::env::var_os("IDENSO_TRACE_DIRECT_SUM_TERM_EXPRESSIONS").is_some()
}

pub(super) fn trace_schoonschip_patterns() -> bool {
    std::env::var_os("IDENSO_TRACE_SCHOONSCHIP_PATTERNS").is_some()
}

pub(super) fn trace_schoonschip_pattern_misses() -> bool {
    std::env::var_os("IDENSO_TRACE_SCHOONSCHIP_PATTERN_MISSES").is_some()
}
fn distribute_expanded_left_sum(expanded_left: &Atom, right: &Atom) -> Atom {
    expanded_left.terms().fold(Atom::Zero, |sum, term| {
        let term = term.to_owned();
        sum + &term * right
    })
}

fn distribute_expanded_right_sum(left: &Atom, expanded_right: &Atom) -> Atom {
    expanded_right.terms().fold(Atom::Zero, |sum, term| {
        let term = term.to_owned();
        sum + left * &term
    })
}

pub(super) fn distribute_smallest_expanded_sum_side(left: &Atom, right: &Atom) -> Atom {
    match (is_sum(left), is_sum(right)) {
        (true, true) => {
            if expression_size(left) <= expression_size(right) {
                distribute_expanded_left_sum(&left.expand(), right)
            } else {
                distribute_expanded_right_sum(left, &right.expand())
            }
        }
        (true, false) => distribute_expanded_left_sum(&left.expand(), right),
        (false, true) => distribute_expanded_right_sum(left, &right.expand()),
        (false, false) => left * right,
    }
}

pub(super) fn multiplicative_factors(expr: AtomView<'_>) -> Vec<Atom> {
    match expr {
        AtomView::Mul(mul) => mul.iter().map(|factor| factor.to_owned()).collect(),
        _ => vec![expr.to_owned()],
    }
}

pub(super) fn product_excluding(factors: &[Atom], excluded: &[bool]) -> Atom {
    factors
        .iter()
        .enumerate()
        .filter(|(index, _)| !excluded[*index])
        .fold(Atom::num(1), |product, (_, factor)| product * factor)
}
