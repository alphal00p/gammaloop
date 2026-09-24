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
            let normalized = BracketNormalizer::normalize(current.as_view()).normalize_dots();
            // Dot normalization stays outside this guard: compact vector
            // rewrites can require no repeated explicit index at all.
            if normalized.has_repeated_explicit_indices() {
                let contracted = SlotContraction::run(
                    normalized.as_view(),
                    self.settings.simplify_chain_like_functions,
                    self.settings.schoonschip_rank1_tensors,
                );
                if contracted != normalized {
                    let next = BracketNormalizer::normalize(contracted.normalize_dots().as_view());
                    // The contractor already reaches its fixed point. If
                    // cleanup changes nothing, no new contraction or
                    // normalization can be exposed by another full pass.
                    if next == contracted || next == current {
                        return next;
                    }
                    current = next;
                    continue;
                }
            }
            // Initial normalization can expose another bracket identity even
            // without a contraction, such as bracket(p(mu)^2) becoming scalar.
            if normalized == current {
                return normalized;
            }
            current = normalized;
        }
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
            "(x+y)^6*g(mink(4,a),mink(4,b))*(unknown(mink(4,b))+another(mink(4,b)))",
            "(g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c)))*epsilon(mink(4,a),mink(4,d),mink(4,e),mink(4,f))",
            "p(mink(4,a))*(p(mink(4,a))+q(mink(4,a)))",
            "bracket(g(mink(4,a),mink(4,b))*p(mink(4,a))*q(mink(4,b)))",
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
