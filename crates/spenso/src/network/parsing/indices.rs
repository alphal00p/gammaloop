//! Syntactic index occurrences, before structure inference removes contractions.

use crate::structure::slot::{SlotMatch, SlotMatcher};
use symbolica::atom::AtomView;

/// Store common symbol IDs separately from arbitrary borrowed index payloads.
struct Indices<'a> {
    symbols: Vec<u32>,
    others: Vec<AtomView<'a>>,
}

#[derive(Default)]
pub(super) struct RepeatedIndices {
    slots: SlotMatcher,
}

impl RepeatedIndices {
    pub(super) fn contains(expr: AtomView<'_>) -> bool {
        let mut scanner = Self::default();
        let mut seen = Indices {
            symbols: Vec::with_capacity(16),
            others: Vec::new(),
        };
        // Top-level alternatives do not need their possible-index sets merged.
        match expr {
            AtomView::Add(sum) => sum.iter().any(|term| {
                seen.symbols.clear();
                seen.others.clear();
                scanner.visit(term, &mut seen)
            }),
            _ => scanner.visit(expr, &mut seen),
        }
    }

    fn visit<'a>(&mut self, expr: AtomView<'a>, seen: &mut Indices<'a>) -> bool {
        match expr {
            AtomView::Num(_) | AtomView::Var(_) => false,
            AtomView::Fun(fun) => {
                match self.slots.classify(expr) {
                    SlotMatch::Explicit(slot) => match slot.index() {
                        AtomView::Var(var) => {
                            let index = var.get_symbol_id();
                            if seen.symbols.contains(&index) {
                                return true;
                            }
                            seen.symbols.push(index);
                            false
                        }
                        index => {
                            if seen.others.contains(&index) {
                                return true;
                            }
                            seen.others.push(index);
                            false
                        }
                    },
                    // Dimensions and index payloads are opaque, including malformed slots.
                    SlotMatch::Opaque => false,
                    SlotMatch::Other => fun.iter().any(|arg| self.visit(arg, seen)),
                }
            }
            AtomView::Mul(product) => product.iter().any(|factor| self.visit(factor, seen)),
            AtomView::Add(sum) => {
                let outer_symbols = seen.symbols.len();
                let outer_others = seen.others.len();
                let mut possible_symbols = Vec::new();
                let mut possible_others = Vec::new();
                for term in sum.iter() {
                    if self.visit(term, seen) {
                        return true;
                    }
                    for &index in &seen.symbols[outer_symbols..] {
                        if !possible_symbols.contains(&index) {
                            possible_symbols.push(index);
                        }
                    }
                    for &index in &seen.others[outer_others..] {
                        if !possible_others.contains(&index) {
                            possible_others.push(index);
                        }
                    }
                    seen.symbols.truncate(outer_symbols);
                    seen.others.truncate(outer_others);
                }
                // Later product factors can meet an index from any summand.
                seen.symbols.extend(possible_symbols);
                seen.others.extend(possible_others);
                false
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let before = seen.symbols.len() + seen.others.len();
                if self.visit(base, seen) {
                    return true;
                }
                if seen.symbols.len() + seen.others.len() > before
                    && match i64::try_from(exponent) {
                        Ok(n) => n.unsigned_abs() > 1,
                        Err(_) => matches!(exponent, AtomView::Num(number)
                            if number.get_coeff_view().is_integer()),
                    }
                {
                    return true;
                }
                self.visit(exponent, seen)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::network::parsing::AtomStructureExt;
    use symbolica::{atom::Atom, parser::ParseSettings};
    fn parse(source: &str) -> Atom {
        Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap()
    }
    #[test]
    fn repeated_explicit_indices_follow_product_and_sum_scopes() {
        crate::structure::representation::initialize();
        let cases = [
            ("g(mink(4,a),mink(4,b))", false),
            ("x*g(mink(4,a),mink(4,b))", false),
            ("g(mink(4,a),mink(4,b))*g(mink(4,c),mink(4,d))", false),
            ("g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c))", false),
            (
                "g(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d),mink(4,e),mink(4,f))",
                false,
            ),
            // Metric traces normalize during construction, before the scan.
            ("g(mink(4,a),mink(4,a))", false),
            ("Unknown(mink(4,a),mink(4,a))", true),
            ("g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))", true),
            (
                "g(mink(4,a),mink(4,b))*epsilon(mink(4,b),mink(4,c),mink(4,d),mink(4,e))",
                true,
            ),
            ("g(mink(4,a),mink(4,b))^2", true),
            ("f(g(mink(4,a),mink(4,a)))", false),
            ("f(Unknown(mink(4,a),mink(4,a)))", true),
            (
                "g(mink(4,a),mink(4,b))*(g(mink(4,b),mink(4,c))+g(mink(4,b),mink(4,d)))",
                true,
            ),
            ("g(mink(4,a),P(mink(4)))", false),
            ("g(P(mink(4)),Q(mink(4)))", false),
            (
                "epsilon(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d))",
                false,
            ),
            ("g(mink(4,a),mink(4,b))*g(mink(5,a),mink(5,b))", true),
            ("Unknown(mink(4,a))*Another(mink(4,b))", false),
            ("Unknown(mink(4,a))*Another(mink(4,a))", true),
            ("f(Unknown(mink(4,a)))", false),
            ("Unknown(mink(4,a))*(A(mink(4,a))+B(mink(4,b)))", true),
            ("Unknown(mink(4,c))*(A(mink(4,a))+B(mink(4,b)))", false),
            (
                "(A(mink(4,a))+B(mink(4,b)))*(C(mink(4,a))+D(mink(4,c)))",
                true,
            ),
            (
                "(A(mink(4,a))+B(mink(4,b)))*(C(mink(4,c))+D(mink(4,d)))",
                false,
            ),
            ("f(A(mink(4,a))+B(mink(4,a)))", false),
            ("T(mink(4,a))+U(mink(4,b),mink(4,b))", true),
            (
                "T(mink(4,a))*(U(mink(4,b))+V(mink(4,c)))*W(mink(4,c))",
                true,
            ),
            (
                "T(mink(4,a))*(U(mink(4,b))+V(mink(4,c)))*W(mink(4,d))",
                false,
            ),
            ("f(mink(4,1))*g(mink(4,1),mink(4,2))", true),
            ("f(1)*h(1)*P(mink(4,a))", false),
            ("T(mink(4,a))^(-2)", true),
            ("T(mink(4,a))^n", false),
            ("(A(mink(4,a))+B(mink(4,b)))^2", true),
            ("T(mink(4,f(a)))*U(mink(4,f(a)))", true),
            ("T(mink(4,f(a)))*U(mink(4,a))", false),
            ("T(mink(f(mink(4,a))))*U(mink(4,a))", false),
            ("T(mink(4,f(mink(4,a))))*U(mink(4,a))", false),
            ("T(mink(4,a,b))*U(mink(4,a))", false),
            ("uind(mink(4,a))*U(mink(4,a))", true),
            ("dind(mink(4,a))*U(mink(4,a))", true),
            ("dind(mink(4,a),extra)*U(mink(4,a))", false),
            ("dind(dind(mink(4,a),extra))*U(mink(4,a))", false),
            ("T(mink(4,a))^18446744073709551616", true),
            ("T(mink(4,a))^(-18446744073709551616)", true),
        ];
        for (source, expected) in cases {
            let expr = parse(source);
            assert_eq!(expr.has_repeated_explicit_indices(), expected, "{source}");
        }
    }
}
