use std::sync::LazyLock;

use rustred::foundry::artifact::{ArtifactPersistenceError, ClosedArtifact};

use super::RustRedEvaluationError;

const K1_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k1.rr");
const K3_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k3.rr");

fn decode_current_artifact(bytes: &[u8]) -> Result<ClosedArtifact, ArtifactPersistenceError> {
    ClosedArtifact::decode_durable(bytes)
}

// The embedded files and the pinned RustRed revision are one atomic build-time
// contract. Decode only the current RustRed schema: obsolete artifact schemas
// intentionally have no Vakint migration, compatibility reader, or fallback.
static K1_ARTIFACT: LazyLock<Result<ClosedArtifact, String>> =
    LazyLock::new(|| decode_current_artifact(K1_BYTES).map_err(|error| error.to_string()));
static K3_ARTIFACT: LazyLock<Result<ClosedArtifact, String>> =
    LazyLock::new(|| decode_current_artifact(K3_BYTES).map_err(|error| error.to_string()));

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

#[cfg(test)]
mod tests {
    use rustred::foundry::artifact::ArtifactPersistenceError;

    use super::{K1_BYTES, decode_current_artifact};

    #[test]
    fn embedded_artifact_contract_rejects_obsolete_schemas() {
        let schema_offset = b"RRIBP\0\r\n".len();
        for schema in [1_u32, 2] {
            let mut obsolete = K1_BYTES.to_vec();
            obsolete[schema_offset..schema_offset + std::mem::size_of::<u32>()]
                .copy_from_slice(&schema.to_le_bytes());
            assert!(matches!(
                decode_current_artifact(&obsolete),
                Err(ArtifactPersistenceError::UnsupportedSchema { actual })
                    if actual == schema
            ));
        }
    }
}
