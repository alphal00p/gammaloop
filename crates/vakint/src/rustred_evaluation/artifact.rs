use std::sync::LazyLock;

use rustred::foundry::artifact::{ArtifactPersistenceError, ArtifactSchemaVersion, ClosedArtifact};

use super::terminal::{TerminalCatalog, TerminalManifest, TerminalSource};
use super::RustRedEvaluationError;

const K1_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k1.rr");
const K3_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k3.rr");
const K6_BYTES: &[u8] = include_bytes!("../../data/rustred/unit_mass_vacuum_k6.rr");

const K1_ALGORITHM_ID: &str = "rustred.generated.one-loop-unit-mass-tadpole.v1";
const K1_FAMILY_FINGERPRINT: &str = "rustred-integral-family-v2;N37:rustred-one-loop-unit-mass-tadpole-v1;L1;1:k;E0;P1;1:d;Q1,0,1;M;RY1,1;I+1;X1;Y1,1;I+1;X0;D1;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;G0;U1;RY1,0;Y1,1;I+1;X0;";
const K3_ALGORITHM_ID: &str = "rustred.generated.two-loop-unit-mass-sunset.v1";
const K3_FAMILY_FINGERPRINT: &str = "rustred-integral-family-v2;N36:rustred-two-loop-unit-mass-sunset-v1;L2;2:k1;2:k2;E0;P1;1:d;Q2,0,3;M;RY1,1;I+1;X1;Y1,1;I+1;X0;D3;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I+2;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;G0;U3;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;";

const K6_ALGORITHM_ID: &str = "rustred.source-port-original-domain.v1";
const K6_FAMILY_FINGERPRINT: &str = "rustred-integral-family-v2;N41:rustred-three-loop-unit-mass-vacuum-k6-v1;L3;2:k1;2:k2;2:k3;E0;P1;1:d;Q3,0,6;M;RY1,1;I+1;X1;Y1,1;I+1;X0;D6;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I-2;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I-2;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I-1;X0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;RY1,1;I-2;X0;Y1,1;I+1;X0;RY1,1;I+1;X0;Y1,1;I+1;X0;G0;U6;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;RY1,0;Y1,1;I+1;X0;";

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
            ArtifactSchemaVersion::V5,
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
            ArtifactSchemaVersion::V5,
            K3_ALGORITHM_ID,
            K3_FAMILY_FINGERPRINT,
            &sources,
        ),
    )
});

static K6_ASSETS: LazyLock<Result<FamilyAssets, String>> = LazyLock::new(|| {
    load_assets(
        K6_BYTES,
        &TerminalManifest::new(
            ArtifactSchemaVersion::V5,
            K6_ALGORITHM_ID,
            K6_FAMILY_FINGERPRINT,
            &super::terminal::k6::SOURCES,
        ),
    )
});

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ArtifactFamily {
    UnitMassVacuumK1,
    UnitMassVacuumK3,
    UnitMassVacuumK6,
}

/// Build-time registry view used by the matcher and future artifact loaders.
///
/// The registry deliberately contains only sealed assets.  A four-loop family
/// must be added here together with its durable bytes, authenticated family
/// fingerprint, and terminal manifest; there is no placeholder artifact and
/// no topology-name dispatch.  Keeping this lookup separate from matching
/// makes that addition an atomic registry change rather than a change to the
/// scalar reducer.
pub(super) const fn shipped_family_for_loop_count(loop_count: usize) -> Option<ArtifactFamily> {
    match loop_count {
        1 => Some(ArtifactFamily::UnitMassVacuumK1),
        2 => Some(ArtifactFamily::UnitMassVacuumK3),
        3 => Some(ArtifactFamily::UnitMassVacuumK6),
        _ => None,
    }
}

impl ArtifactFamily {
    pub(super) const fn name(self) -> &'static str {
        match self {
            Self::UnitMassVacuumK1 => "unit-mass-vacuum-k1",
            Self::UnitMassVacuumK3 => "unit-mass-vacuum-k3",
            Self::UnitMassVacuumK6 => "unit-mass-vacuum-k6",
        }
    }

    pub(super) fn artifact(self) -> Result<&'static ClosedArtifact, RustRedEvaluationError> {
        let loaded = match self {
            Self::UnitMassVacuumK1 => &*K1_ASSETS,
            Self::UnitMassVacuumK3 => &*K3_ASSETS,
            Self::UnitMassVacuumK6 => &*K6_ASSETS,
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
            Self::UnitMassVacuumK6 => &*K6_ASSETS,
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
    use rustred::foundry::artifact::{
        derive_one_loop_unit_mass_tadpole, derive_two_loop_unit_mass_sunset,
        ArtifactPersistenceError,
    };

    use super::{
        decode_current_artifact, shipped_family_for_loop_count, ArtifactFamily, K1_BYTES, K3_BYTES,
    };

    #[test]
    fn registry_exposes_only_sealed_assets_until_four_loop_artifact_exists() {
        assert_eq!(
            shipped_family_for_loop_count(1),
            Some(ArtifactFamily::UnitMassVacuumK1)
        );
        assert_eq!(
            shipped_family_for_loop_count(2),
            Some(ArtifactFamily::UnitMassVacuumK3)
        );
        assert_eq!(
            shipped_family_for_loop_count(3),
            Some(ArtifactFamily::UnitMassVacuumK6)
        );
        assert_eq!(shipped_family_for_loop_count(4), None);
    }

    #[test]
    fn embedded_artifacts_match_current_rustred_generators() {
        let k1 = derive_one_loop_unit_mass_tadpole()
            .expect("derive the registered K=1 artifact")
            .encode_durable()
            .expect("encode the registered K=1 artifact");
        let k3 = derive_two_loop_unit_mass_sunset()
            .expect("derive the registered K=3 artifact")
            .encode_durable()
            .expect("encode the registered K=3 artifact");

        assert_eq!(k1, K1_BYTES);
        assert_eq!(k3, K3_BYTES);
    }

    #[test]
    fn embedded_artifact_contract_rejects_obsolete_schemas() {
        let schema_offset = b"RRIBP\0\r\n".len();
        for schema in [1_u32, 2, 3, 4] {
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
