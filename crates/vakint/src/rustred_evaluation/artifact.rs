use std::sync::LazyLock;

use rustred::foundry::artifact::ClosedArtifact;

use super::RustRedEvaluationError;

const K1_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k1.rr");
const K3_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k3.rr");

// The embedded files and the pinned RustRed revision are one atomic build-time
// contract. Decode only the current RustRed schema: obsolete artifact schemas
// intentionally have no Vakint migration, compatibility reader, or fallback.
static K1_ARTIFACT: LazyLock<Result<ClosedArtifact, String>> =
    LazyLock::new(|| ClosedArtifact::decode_durable(K1_BYTES).map_err(|error| error.to_string()));
static K3_ARTIFACT: LazyLock<Result<ClosedArtifact, String>> =
    LazyLock::new(|| ClosedArtifact::decode_durable(K3_BYTES).map_err(|error| error.to_string()));

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ArtifactFamily {
    UnitMassVacuumK1,
    UnitMassVacuumK3,
}

impl ArtifactFamily {
    pub(super) const fn name(self) -> &'static str {
        match self {
            Self::UnitMassVacuumK1 => "unit-mass-vacuum-k1",
            Self::UnitMassVacuumK3 => "unit-mass-vacuum-k3",
        }
    }

    pub(super) fn artifact(self) -> Result<&'static ClosedArtifact, RustRedEvaluationError> {
        let loaded = match self {
            Self::UnitMassVacuumK1 => &*K1_ARTIFACT,
            Self::UnitMassVacuumK3 => &*K3_ARTIFACT,
        };
        loaded
            .as_ref()
            .map_err(|detail| RustRedEvaluationError::ArtifactLoad {
                family: self.name(),
                detail: detail.clone(),
            })
    }

    pub(super) fn validate_root_powers(self, powers: &[i64]) -> Result<(), RustRedEvaluationError> {
        let artifact = self.artifact()?;
        if powers.len() != artifact.arity() {
            return Err(RustRedEvaluationError::InvalidMatchedFamily {
                detail: format!(
                    "matched family has {} powers but artifact {} requires {}",
                    powers.len(),
                    self.name(),
                    artifact.arity()
                ),
            });
        }
        for (position, (&value, bounds)) in powers
            .iter()
            .zip(artifact.supported_root_power_bounds())
            .enumerate()
        {
            if !bounds.contains(value) {
                return Err(RustRedEvaluationError::InvalidMatchedFamily {
                    detail: format!(
                        "power {value} at position {position} is outside artifact {} root domain [{}, {}]",
                        self.name(),
                        bounds.lower(),
                        bounds.upper()
                    ),
                });
            }
        }
        Ok(())
    }
}
