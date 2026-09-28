use symbolica::atom::{Atom, AtomView};

use super::settings::SchoonschipSettings;
use crate::shorthands::schoonschip::{DotNormalizer, with_settings::SchoonschipWithSettings};

pub(crate) trait Schoonschip {
    fn schoonschip(&self) -> Atom;

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom;

    fn normalize_dots(&self) -> Atom;

    fn to_dots(&self) -> Atom;
}

impl Schoonschip for Atom {
    fn schoonschip(&self) -> Atom {
        self.as_view().schoonschip()
    }

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom {
        self.as_view().schoonschip_with_settings(settings)
    }

    fn normalize_dots(&self) -> Atom {
        self.as_view().normalize_dots()
    }

    fn to_dots(&self) -> Atom {
        self.as_view().to_dots()
    }
}

impl Schoonschip for AtomView<'_> {
    fn normalize_dots(&self) -> Atom {
        super::normalize_dots::DotNormalizer::run(*self)
    }

    fn schoonschip(&self) -> Atom {
        self.schoonschip_with_settings(&SchoonschipSettings::default())
    }

    fn schoonschip_with_settings(&self, settings: &SchoonschipSettings) -> Atom {
        SchoonschipWithSettings { settings }.run(*self, &mut Vec::new())
    }

    fn to_dots(&self) -> Atom {
        // Resolve explicit vector pairs before tensor collection can absorb a
        // registered momentum into an ordinary, untagged vector head.
        let explicit = crate::shorthands::metric::to_dots_impl(*self);
        let simplified = explicit
            .schoonschip_with_settings(&SchoonschipSettings::default().with_rank1_tensors());
        let simplified = DotNormalizer::metric_shorthand_to_dot(simplified.as_view());
        // Explicit representation slots also identify vectors whose heads
        // were created as ordinary Symbolica symbols without rank-one tags.
        crate::shorthands::metric::to_dots_impl(simplified.as_view())
    }
}
