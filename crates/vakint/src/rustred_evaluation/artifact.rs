use std::sync::LazyLock;

use rustred::foundry::artifact::{ArtifactPersistenceError, ArtifactSchemaVersion, ClosedArtifact};

use super::RustRedEvaluationError;
use super::terminal::{TerminalCatalog, TerminalManifest, TerminalSource};

const K1_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k1.rr");
const K3_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k3.rr");

const K1_ALGORITHM_ID: &str = "rustred.generated.one-loop-unit-mass-tadpole.v1";
const K1_FAMILY_FINGERPRINT: &str = "rustred-integral-family-v2;N37:rustred-one-loop-unit-mass-tadpole-v1;L1;1:k;E0;P1;1:d;Q1,0,1;M;RY1,1;I+1;X1;Y1,1;I+1;X0;D1;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;G0;U1;RY1,0;Y1,1;I+1;X0;";
const K3_ALGORITHM_ID: &str = "rustred.generated.two-loop-unit-mass-sunset.v1";
const K3_FAMILY_FINGERPRINT: &str = "rustred-integral-family-v2;N36:rustred-two-loop-unit-mass-sunset-v1;L2;2:k1;2:k2;E0;P1;1:d;Q2,0,3;M;RY1,1;I+1;X1;Y1,1;I+1;X0;D3;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I+2;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;G0;U3;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;";

fn decode_current_artifact(bytes: &[u8]) -> Result<ClosedArtifact, ArtifactPersistenceError> {
    ClosedArtifact::decode_durable(bytes)
}

// The embedded files and the pinned RustRed revision are one atomic build-time
// contract. Decode only the current RustRed schema: obsolete artifact schemas
// intentionally have no Vakint migration, compatibility reader, or fallback.
struct FamilyAssets {
    artifact: ClosedArtifact,
    terminals: TerminalCatalog,
}

fn load_assets(
    bytes: &[u8],
    terminal_manifest: &TerminalManifest<'_>,
) -> Result<FamilyAssets, String> {
    let artifact = decode_current_artifact(bytes).map_err(|error| error.to_string())?;
    let terminals = TerminalCatalog::compile(&artifact, terminal_manifest)?;
    Ok(FamilyAssets {
        artifact,
        terminals,
    })
}

static K1_ASSETS: LazyLock<Result<FamilyAssets, String>> = LazyLock::new(|| {
    let sources = [TerminalSource::exact_matad_basis(
        &[1],
        "-Gam(1,1)/(ep*(ep-1))",
    )];
    load_assets(
        K1_BYTES,
        &TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            K1_ALGORITHM_ID,
            K1_FAMILY_FINGERPRINT,
            &sources,
        ),
    )
});
static K3_ASSETS: LazyLock<Result<FamilyAssets, String>> = LazyLock::new(|| {
    let sources = [
        TerminalSource::exact_matad_basis(&[0, 1, 1], "(Gam(1,1)/(ep*(ep-1)))^2"),
        TerminalSource::exact_matad_basis(&[1, 1, 1], "-miT111"),
    ];
    load_assets(
        K3_BYTES,
        &TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            K3_ALGORITHM_ID,
            K3_FAMILY_FINGERPRINT,
            &sources,
        ),
    )
});

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
            Self::UnitMassVacuumK1 => &*K1_ASSETS,
            Self::UnitMassVacuumK3 => &*K3_ASSETS,
        };
        loaded
            .as_ref()
            .map(|assets| &assets.artifact)
            .map_err(|detail| RustRedEvaluationError::ArtifactLoad {
                family: self.name(),
                detail: detail.clone(),
            })
    }

    pub(super) fn terminals(self) -> Result<&'static TerminalCatalog, RustRedEvaluationError> {
        let loaded = match self {
            Self::UnitMassVacuumK1 => &*K1_ASSETS,
            Self::UnitMassVacuumK3 => &*K3_ASSETS,
        };
        loaded
            .as_ref()
            .map(|assets| &assets.terminals)
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
