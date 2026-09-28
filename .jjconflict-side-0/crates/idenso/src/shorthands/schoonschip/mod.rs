mod analysis;
#[cfg(test)]
mod api;
mod normalize_dots;
mod settings;
mod slot_contraction;
mod with_settings;

#[cfg(test)]
mod test;

pub(crate) use analysis::SimplificationCandidates;
// The legacy raw pipeline exists solely as an independent unit-test oracle.
#[cfg(test)]
pub(crate) use api::Schoonschip;
pub(crate) use normalize_dots::DotNormalizer;
pub(crate) use settings::SchoonschipSettings;
pub(crate) use slot_contraction::{ContractionStatus, FactorizedContraction, SlotContraction};
pub(crate) use with_settings::SchoonschipWithSettings;
