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
        let simplified = if normalized.has_repeated_explicit_indices() {
            SlotContraction::run(
                normalized.as_view(),
                self.settings.simplify_chain_like_functions,
                self.settings.schoonschip_rank1_tensors,
            )
        } else {
            normalized
        };
        BracketNormalizer::normalize(simplified.normalize_dots().as_view())
    }
}
