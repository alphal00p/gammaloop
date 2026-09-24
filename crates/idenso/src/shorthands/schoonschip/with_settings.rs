use spenso::network::parsing::AtomStructureExt;
use symbolica::atom::{Atom, AtomView};

use crate::shorthands::bracket::BracketNormalizer;

use super::{api::Schoonschip, settings::SchoonschipSettings, slot_contraction::SlotContraction};

pub(crate) struct SchoonschipWithSettings<'a> {
    pub(crate) settings: &'a SchoonschipSettings,
}

impl SchoonschipWithSettings<'_> {
    pub(crate) fn run(&self, view: AtomView<'_>) -> Atom {
        let mut current = view.to_owned();
        loop {
            let next = self.apply_once(current.as_view());
            if next == current {
                return next;
            }
            current = next;
        }
    }

    fn apply_once(&self, view: AtomView<'_>) -> Atom {
        let normalized = BracketNormalizer::normalize(view).normalize_dots();
        // Metric and vector contractions both consume two explicit occurrences.
        // Dot normalization stays outside this guard: compact vector rewrites
        // can require no repeated explicit index at all.
        if !normalized.has_repeated_explicit_indices() {
            return normalized;
        }
        let contracted = SlotContraction::run(
            normalized.as_view(),
            self.settings.simplify_chain_like_functions,
            self.settings.schoonschip_rank1_tensors,
        );
        if contracted == normalized {
            return contracted;
        }
        // Only a contraction can introduce new dot or bracket simplifications
        // here. Initial normalization changes are still covered by run's fixed
        // point, including powers that turn a bracket's payload into a scalar.
        BracketNormalizer::normalize(contracted.normalize_dots().as_view())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parser::ParseSettings;

    #[test]
    fn unchanged_contractions_preserve_normalization_fixed_points() {
        crate::test_support::test_initialize();
        let _ = spenso::p!(spenso::mink!(4));
        let _ = spenso::q!(spenso::mink!(4));
        let sources = [
            "g(mink(4,a),mink(4,b))*epsilon(mink(4,c),mink(4,d),mink(4,e),mink(4,f))",
            "bracket(p(mink(4,a))^2)",
            "bracket(g(mink(4,a),mink(4,b))^2)",
            "bracket(p(mink(4,a)),p(mink(4,a)))",
            "bracket(bracket(g(mink(4,a),mink(4,b))^2),p(mink(4,c)))",
            "bracket(p(mink(4,a))^2+q(mink(4,a))^2)^2",
            "g(mink(4,a),mink(4,b))*unknown(mink(4,b))",
            "p(mink(4,a))*q(mink(4,a))",
            "unknown(mink(4,a))*another(mink(4,a))",
            "bracket(g(mink(4,a),mink(4,b)),unknown(mink(4,b)))",
        ];
        for (mode, settings) in [
            SchoonschipSettings::default(),
            SchoonschipSettings::default().with_chain_like_functions(),
            SchoonschipSettings::default().without_rank1_tensors(),
        ]
        .into_iter()
        .enumerate()
        {
            for source in sources {
                let expression = Atom::parse(source, "spenso", ParseSettings::symbolica()).unwrap();
                let mut expected = expression.clone();
                loop {
                    // Preserve the previous eager cleanup schedule as an oracle
                    // for scalar powers and unresolved indexed bracket scopes.
                    let normalized =
                        BracketNormalizer::normalize(expected.as_view()).normalize_dots();
                    let contracted = if normalized.has_repeated_explicit_indices() {
                        SlotContraction::run(
                            normalized.as_view(),
                            settings.simplify_chain_like_functions,
                            settings.schoonschip_rank1_tensors,
                        )
                    } else {
                        normalized
                    };
                    let next = BracketNormalizer::normalize(contracted.normalize_dots().as_view());
                    if next == expected {
                        break;
                    }
                    expected = next;
                }
                let actual = expression.schoonschip_with_settings(&settings);
                assert_eq!(actual, expected, "{source}, mode {mode}");
                assert_eq!(
                    actual.schoonschip_with_settings(&settings),
                    actual,
                    "{source}"
                );
            }
        }
    }
}
