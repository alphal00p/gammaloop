use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
};

use crate::{
    utils::{GS, W_, symbolica_ext::CallSymbol},
    uv::UVgenerationSettings,
};

#[derive(Clone, Copy, Debug)]
pub(crate) enum UvOperation {
    Approx,
    Integrate,
    Series,
    Truncate,
}

impl UvOperation {
    fn symbol(self) -> Symbol {
        match self {
            Self::Approx => GS.uv_approx,
            Self::Integrate => GS.uv_integrate,
            Self::Series => GS.uv_series,
            Self::Truncate => GS.uv_truncate,
        }
    }

    fn apply(self, history: Atom) -> Atom {
        self.symbol().f(history)
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct UvMarker {
    enabled: bool,
    keep: bool,
}

impl UvMarker {
    pub(crate) fn new(settings: &UVgenerationSettings) -> Self {
        Self {
            enabled: settings.add_marker,
            keep: settings.keep_marker,
        }
    }

    fn label(subgraph: &SuBitGraph) -> Atom {
        if subgraph.is_empty() {
            Atom::Zero
        } else {
            Atom::var(subgraph.symbol())
        }
    }

    pub(crate) fn subgraph(current: &SuBitGraph, given: &SuBitGraph) -> Atom {
        function!(GS.uv_subgraph, Self::label(current), Self::label(given))
    }

    fn normalize(atom: &Atom) -> Atom {
        atom.replace(function!(GS.ct_marker, W_.a_).pow(Atom::var(W_.b_)))
            .when(W_.b_.filter(|exponent| exponent.is_integer()))
            .with(function!(
                GS.ct_marker,
                Atom::var(W_.a_).pow(Atom::var(W_.b_))
            ))
            .replace(function!(GS.ct_marker, W_.a_) * function!(GS.ct_marker, W_.b_))
            .repeat()
            .with(function!(GS.ct_marker, W_.a_ * W_.b_))
    }

    pub(crate) fn apply(
        self,
        operation: UvOperation,
        current: &SuBitGraph,
        given: &SuBitGraph,
        atom: &Atom,
    ) -> Atom {
        if !self.enabled || atom.is_zero() {
            return atom.clone();
        }

        let shadow = Self::subgraph(current, given);
        atom.map_terms_single_core(|term| {
            let normalized = Self::normalize(&term.to_owned());
            // Approximation enters a subgraph; later operations wrap that same history.
            let history = if matches!(operation, UvOperation::Approx) {
                &shadow * W_.a_
            } else {
                Atom::var(W_.a_)
            };
            let rewritten = normalized
                .replace(function!(GS.ct_marker, W_.a_))
                .with(function!(GS.ct_marker, operation.apply(history)));
            if rewritten == normalized {
                normalized * function!(GS.ct_marker, operation.apply(shadow.clone()))
            } else {
                rewritten
            }
        })
    }

    pub(crate) fn prefix(self, current: &SuBitGraph, given: &SuBitGraph, atom: &Atom) -> Atom {
        if !self.enabled || atom.is_zero() {
            return atom.clone();
        }

        let shadow = Self::subgraph(current, given);
        atom.map_terms_single_core(|term| {
            let normalized = Self::normalize(&term.to_owned());
            let rewritten = normalized
                .replace(function!(GS.ct_marker, W_.a_))
                .with(function!(GS.ct_marker, &shadow * W_.a_));
            if rewritten == normalized {
                normalized * function!(GS.ct_marker, shadow.clone())
            } else {
                rewritten
            }
        })
    }

    pub(crate) fn finish(self, atom: &Atom) -> Atom {
        if self.enabled && !self.keep {
            atom.replace(function!(GS.ct_marker, W_.a_))
                .with(Atom::one())
        } else {
            atom.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
    use symbolica::{atom::AtomCore, function, symbol};

    use super::*;

    fn marker(keep: bool) -> UvMarker {
        UvMarker::new(&UVgenerationSettings {
            add_marker: true,
            keep_marker: keep,
            ..Default::default()
        })
    }

    #[test]
    fn operations_nest_marker_history() {
        let graph = SuBitGraph::empty(1);
        let marker = marker(true);
        let mut actual = Atom::num(2);
        for operation in [
            UvOperation::Approx,
            UvOperation::Integrate,
            UvOperation::Series,
            UvOperation::Truncate,
        ] {
            actual = marker.apply(operation, &graph, &graph, &actual);
        }

        let subgraph = UvMarker::subgraph(&graph, &graph);
        let history = UvOperation::Truncate.apply(
            UvOperation::Series
                .apply(UvOperation::Integrate.apply(UvOperation::Approx.apply(subgraph))),
        );
        assert_eq!(actual, Atom::num(2) * function!(GS.ct_marker, history));
    }

    #[test]
    fn marker_products_merge_and_strip() {
        let a = Atom::var(symbol!("a"));
        let b = Atom::var(symbol!("b"));
        let product = function!(GS.ct_marker, &a) * function!(GS.ct_marker, &b);

        assert_eq!(
            UvMarker::normalize(&product),
            function!(GS.ct_marker, &a * b)
        );
        assert_eq!(
            UvMarker::normalize(&(function!(GS.ct_marker, &a) * function!(GS.ct_marker, &a))),
            function!(GS.ct_marker, a.pow(2))
        );
        assert_eq!(marker(false).finish(&product), Atom::one());
    }

    #[test]
    fn every_term_receives_operation_and_subgraph_history() {
        let graph = SuBitGraph::empty(1);
        let marker = marker(true);
        let a = Atom::var(symbol!("a"));
        let mixed = Atom::num(2) + Atom::num(3) * function!(GS.ct_marker, a.clone());
        let subgraph = UvMarker::subgraph(&graph, &graph);

        assert_eq!(
            marker.apply(UvOperation::Approx, &graph, &graph, &mixed),
            Atom::num(2) * function!(GS.ct_marker, UvOperation::Approx.apply(subgraph.clone()))
                + Atom::num(3) * function!(GS.ct_marker, UvOperation::Approx.apply(&subgraph * &a))
        );
        assert_eq!(
            marker.prefix(&graph, &graph, &mixed),
            Atom::num(2) * function!(GS.ct_marker, subgraph.clone())
                + Atom::num(3) * function!(GS.ct_marker, subgraph * a)
        );
    }

    #[test]
    fn uv_heads_have_spenso_shorthands() {
        for (symbol, name) in [
            (GS.uv_subgraph, "gammalooprs::uv::subgraph"),
            (GS.uv_approx, "gammalooprs::uv::Approx"),
            (GS.uv_integrate, "gammalooprs::uv::Integrate"),
            (GS.uv_series, "gammalooprs::uv::Series"),
            (GS.uv_truncate, "gammalooprs::uv::Truncate"),
        ] {
            assert_eq!(symbol.get_name(), name);
        }

        let current = Atom::var(symbol!("S_47"));
        let given = Atom::var(symbol!("S_3s"));
        let subgraph = function!(GS.uv_subgraph, current, given);
        let operations = function!(
            GS.uv_truncate,
            function!(
                GS.uv_series,
                function!(GS.uv_integrate, function!(GS.uv_approx, Atom::one()))
            )
        );
        let settings = SpensoPrintSettings::compact();
        let canonical = (&subgraph * &operations).to_canonical_string();
        let marked = function!(GS.ct_marker, operations.clone());

        for name in ["subgraph", "Approx", "Integrate", "Series", "Truncate"] {
            assert!(canonical.contains(name));
        }
        assert!(!canonical.contains('⊛'));
        assert!(!canonical.contains("K["));
        assert!(marked.to_canonical_string().contains("CT"));
        assert_eq!(
            subgraph.printer(settings.nice_symbolica()).to_string(),
            "S_47⊛3s"
        );
        assert_eq!(
            operations.printer(settings.nice_symbolica()).to_string(),
            "Tr(Σ(⟨K[1]⟩))"
        );
        assert_eq!(
            marked.printer(settings.nice_symbolica()).to_string(),
            nu_ansi_term::Color::Magenta
                .paint("Tr(Σ(⟨K[1]⟩))")
                .to_string()
        );
        assert_eq!(
            marked.printer(settings.typst_symbolica()).to_string(),
            "Tr(Σ(⟨K[1]⟩))"
        );
    }
}
