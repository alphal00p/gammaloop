mod api;
mod contraction;
mod normalize_dots;
mod settings;
mod utils;
mod with_settings;

#[cfg(test)]
mod test;

pub use api::Schoonschip;
pub use contraction::Schoonschipify;
pub(crate) use normalize_dots::DotNormalizer;
pub use settings::{SchoonschipContractionOrder, SchoonschipSettings, SchoonschipTraversal};
