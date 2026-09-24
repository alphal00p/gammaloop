mod analysis;
mod api;
mod contraction;
mod normalize_dots;
mod settings;
mod slot_contraction;
mod utils;
mod with_settings;

#[cfg(test)]
mod test;

pub(crate) use analysis::SimplificationCandidates;
pub use api::Schoonschip;
pub use contraction::Schoonschipify;
pub(crate) use normalize_dots::DotNormalizer;
pub use settings::{SchoonschipContractionOrder, SchoonschipSettings, SchoonschipTraversal};
