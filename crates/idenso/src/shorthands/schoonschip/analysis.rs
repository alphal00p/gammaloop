//! Share the syntactic checks that decide whether normalization has any work.

use spenso::{
    network::{parsing::AtomStructureExt, tags::SPENSO_TAG},
    structure::slot::SlotMatch,
};
use symbolica::atom::{AtomView, Symbol};

/// Call-local observations of one expression, invalidated by any rewrite.
/// A compound opaque payload leaves head observations incomplete; consumers
/// must then defer to the individual simplifiers.
pub(crate) struct SimplificationCandidates<const N: usize> {
    pub(crate) repeated_indices: bool,
    pub(crate) brackets: bool,
    pub(crate) dots: bool,
    pub(crate) symbols: [bool; N],
    pub(crate) complete: bool,
}

impl<const N: usize> SimplificationCandidates<N> {
    pub(crate) fn scan(expression: AtomView<'_>, symbols: [Symbol; N]) -> Self {
        let symbols = symbols.map(|symbol| symbol.get_id());
        let bracket = SPENSO_TAG.bracket.get_id();
        let rank_one = &SPENSO_TAG.rank1;
        let mut rank_one_heads = [None; 16];
        let mut candidates = Self {
            repeated_indices: false,
            brackets: false,
            dots: false,
            symbols: [false; N],
            complete: true,
        };
        candidates.repeated_indices =
            expression.has_repeated_explicit_indices_with_observer(|node, slot| match node {
                AtomView::Fun(function) => {
                    let id = function.get_symbol_id();
                    candidates.observe_symbol(id, bracket, &symbols);
                    let entry = &mut rank_one_heads[(id.wrapping_mul(0x9e37_79b9) >> 28) as usize];
                    let tagged = match *entry {
                        Some((cached, tagged)) if cached == id => tagged,
                        _ => {
                            let tagged = function.get_symbol().has_tag(rank_one);
                            *entry = Some((id, tagged));
                            tagged
                        }
                    };
                    candidates.dots |= tagged;
                    match slot {
                        SlotMatch::Explicit(slot) => {
                            // A variance wrapper can hide the representation
                            // head from the observer's outer function node.
                            if slot.representation().wrapper().is_some() {
                                candidates.observe_symbol(
                                    slot.representation().head().get_id(),
                                    bracket,
                                    &symbols,
                                );
                            }
                            if !candidates.observe_leaf(slot.dimension(), bracket, &symbols)
                                || !candidates.observe_leaf(slot.index(), bracket, &symbols)
                            {
                                candidates.allow_all();
                            }
                        }
                        SlotMatch::Opaque => {
                            if function.iter().any(|argument| {
                                !candidates.observe_leaf(argument, bracket, &symbols)
                            }) {
                                candidates.allow_all();
                            }
                        }
                        SlotMatch::Other => {}
                    }
                }
                AtomView::Var(variable) => {
                    candidates.observe_symbol(variable.get_symbol_id(), bracket, &symbols);
                }
                AtomView::Pow(_) => candidates.dots = true,
                _ => {}
            });
        candidates
    }

    pub(crate) fn normalized(&self) -> bool {
        !self.repeated_indices && !self.brackets && !self.dots
    }

    fn observe_symbol(&mut self, id: u32, bracket: u32, symbols: &[u32; N]) {
        if id == bracket {
            self.brackets = true;
        }
        if symbols.contains(&id) {
            for (seen, &symbol) in self.symbols.iter_mut().zip(symbols) {
                *seen |= id == symbol;
            }
        }
    }

    fn observe_leaf(&mut self, node: AtomView<'_>, bracket: u32, symbols: &[u32; N]) -> bool {
        match node {
            AtomView::Var(variable) => {
                self.observe_symbol(variable.get_symbol_id(), bracket, symbols);
                true
            }
            AtomView::Num(_) => true,
            _ => false,
        }
    }

    fn allow_all(&mut self) {
        self.brackets = true;
        self.dots = true;
        self.complete = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shorthands::{bracket::BracketNormalizer, schoonschip::Schoonschip};
    use symbolica::{
        atom::{Atom, AtomCore},
        parser::ParseSettings,
    };

    #[test]
    fn completed_observations_preserve_heads_and_normalization_boundaries() {
        crate::test_support::test_initialize();
        let _ = spenso::p!(spenso::mink!(4));
        let heads = [
            SPENSO_TAG.chain,
            SPENSO_TAG.trace,
            SPENSO_TAG.bracket,
            *crate::epsilon::EPSILON_SYMBOL,
        ];
        for source in [
            "g(mink(D,a),mink(D,b))",
            "g(mink(D,a),mink(D,b))*(g(mink(D,c),mink(D,d))+g(mink(D,c),mink(D,e)))",
            "T(dind(lor(4,a)))",
            "T(mink(spenso::chain,a))",
            "T(mink(4,spenso::trace))",
            "T(mink(4,a,spenso::bracket))",
            "T(dind(mink(4)))",
            "T(mink(spenso::trace(4),a))",
            "T(mink(4,f(spenso::epsilon)))",
            "g(p(mink(4)),p(mink(4)))",
            "p(mink(4,a))^2",
            "(x+y)^6*g(mink(D,a),mink(D,b))",
            "bracket(p(mink(4,a))^2)",
            "g(mink(4,a),mink(4,b))*T(mink(4,b))*later(spenso::trace)",
            "scope(mink(4,a),mink(4,a),bracket(p(mink(4,b))^2))",
            "scope(mink(4,a),mink(4,a),p(q(mink(4))))",
            "scope(mink(4,a),mink(4,a),T(mink(4,f(spenso::epsilon))))",
        ] {
            let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
            let candidates = SimplificationCandidates::scan(expression.as_view(), heads);
            assert_eq!(
                candidates.repeated_indices,
                expression.has_repeated_explicit_indices(),
                "{source}"
            );
            if candidates.complete {
                assert_eq!(
                    candidates.symbols,
                    heads.map(|head| expression.contains_symbol(head)),
                    "{source}"
                );
            }
            if !candidates.dots {
                assert_eq!(expression.normalize_dots(), expression, "{source}");
            }
            if !candidates.brackets {
                assert_eq!(
                    BracketNormalizer::normalize(expression.as_view()),
                    expression,
                    "{source}"
                );
            }
        }
        // A wrapped representation is still visible to the old full head scan.
        use spenso::structure::representation::LibraryRep;
        let representation: LibraryRep = crate::representations::Bispinor {}.into();
        let head = representation.symbol();
        let expression =
            Atom::parse("T(dind(bis(4,a)))", "spenso", ParseSettings::symbolica()).unwrap();
        let candidates = SimplificationCandidates::scan(expression.as_view(), [head]);
        assert!(!candidates.complete || candidates.symbols == [expression.contains_symbol(head)]);
    }
}
