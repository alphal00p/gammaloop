use spenso::{
    network::library::symbolic::ETS,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::{Pattern, Replacement},
    symbol,
};

use crate::{rep_symbols::RS, representations::Bispinor};

pub trait SelectiveExpand {
    fn expand_in_patterns(&self, pats: &[Pattern]) -> Vec<(Atom, Atom)>;

    fn expand_metrics(&self) -> Vec<(Atom, Atom)> {
        let metric_pat = function!(ETS.metric, RS.a__).to_pattern();
        let id_pat = function!(ETS.metric, RS.a__).to_pattern();

        self.expand_in_patterns(&[metric_pat, id_pat])
    }

    fn expand_bis(&self) -> Vec<(Atom, Atom)> {
        let bis = Bispinor {};
        let bis_pat = function!(RS.f_, RS.a___, bis.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[bis_pat])
    }

    fn expand_mink(&self) -> Vec<(Atom, Atom)> {
        let mink = Minkowski {};
        let mink_pat = function!(RS.f_, RS.a___, mink.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[mink_pat])
    }

    fn expand_mink_bis(&self) -> Vec<(Atom, Atom)> {
        let mink = Minkowski {};
        let mink_pat = function!(RS.f_, RS.a___, mink.to_symbolic([RS.b__]), RS.c___).to_pattern();

        let bis = Bispinor {};
        let bis_pat = function!(RS.f_, RS.a___, bis.to_symbolic([RS.b__]), RS.c___).to_pattern();

        self.expand_in_patterns(&[mink_pat, bis_pat])
    }
}

impl SelectiveExpand for Atom {
    fn expand_in_patterns(&self, pats: &[Pattern]) -> Vec<(Atom, Atom)> {
        self.as_view().expand_in_patterns(pats)
    }
}

impl SelectiveExpand for AtomView<'_> {
    fn expand_in_patterns(&self, pats: &[Pattern]) -> Vec<(Atom, Atom)> {
        let mut coefs = Vec::new();

        // A (x+y)(z*B+x*C) => A(x*x*C+y*x*C+y*z*B+y*z*B)
        for p in pats {
            for m in self.pattern_match(p, None, None) {
                let matched = p.replace_wildcards(&m).unwrap();
                let coef = match matched.as_view() {
                    AtomView::Mul(mul) if mul.has_coefficient() => {
                        let mut coef = Atom::num(1);
                        for factor in mul.iter() {
                            if !matches!(factor, AtomView::Num(_)) {
                                coef *= factor.to_owned();
                            }
                        }
                        coef
                    }
                    _ => matched,
                };
                if !coefs.iter().any(|known| known == &coef) {
                    coefs.push(coef);
                }
            }
        }

        // Symbolica's coefficient collector requires polynomial variables.
        // Tensor factors can contain polynomial dimensions, so collect through placeholders.
        let replacement_pairs = coefs
            .iter()
            .enumerate()
            .map(|(i, coef)| {
                (
                    coef,
                    Atom::var(symbol!(format!("idenso::selective_expand_{i}"))),
                )
            })
            .collect::<Vec<_>>();

        let to_placeholders = replacement_pairs
            .iter()
            .map(|(coef, placeholder)| Replacement::new(coef.to_pattern(), placeholder.clone()))
            .collect::<Vec<_>>();
        let from_placeholders = replacement_pairs
            .iter()
            .map(|(coef, placeholder)| Replacement::new(placeholder.to_pattern(), (*coef).clone()))
            .collect::<Vec<_>>();
        let placeholders = replacement_pairs
            .iter()
            .map(|(_, placeholder)| placeholder.clone())
            .collect::<Vec<_>>();

        self.replace_multiple(&to_placeholders)
            .coefficient_list::<i8>(&placeholders)
            .into_iter()
            .map(|(key, coefficient)| (key.replace_multiple(&from_placeholders), coefficient))
            .collect()
    }
}
