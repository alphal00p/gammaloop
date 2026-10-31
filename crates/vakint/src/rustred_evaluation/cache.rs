//! Persistent reduction-cache boundary.
//!
//! Vakint reduces only with RustRed's ordered scheduler. The depth-first
//! completed-result caches (prepared bundles, running caches and their
//! transport) were removed together with depth-first reduction: a completed
//! terminal image replaces a frontier key and would skip every grouping and
//! cancellation below it. A persistent Ordered cache stores weight-independent
//! one-step records instead and attaches at this boundary, bound to the
//! application identity below.
use std::{fs, path::Path};

use sha2::{Digest, Sha256};

use super::{RustRedEvaluationError, data::RustRedData};
use crate::VakintError;

pub(super) fn operation(error: rustred_app::AppError) -> RustRedEvaluationError {
    RustRedEvaluationError::CacheOperation { source: error }
}

fn failure(error: impl std::fmt::Display) -> VakintError {
    operation(rustred_app::AppError::new(
        rustred_app::AppErrorKind::Input,
        error.to_string(),
    ))
    .into()
}

/// An application identity deliberately excludes cache inventory, release names
/// and numerical master coefficients: none changes a raw Route/Apply record.
pub(super) fn saved_application_digest(
    data: &RustRedData,
    path: &Path,
) -> Result<[u8; 32], VakintError> {
    let manifest = if path.is_dir() {
        path.join("application.json")
    } else {
        path.to_owned()
    };
    let application = fs::read(manifest).map_err(failure)?;
    let mut digest = Sha256::new();
    digest.update(b"vakint-saved-application-cache-v1\0");
    digest.update((application.len() as u64).to_le_bytes());
    digest.update(application);
    if let Some(pack) = data.loops.get(&5) {
        fn field(digest: &mut Sha256, bytes: &[u8]) {
            digest.update((bytes.len() as u64).to_le_bytes());
            digest.update(bytes);
        }
        field(&mut digest, pack.manifest.dimension.as_bytes());
        field(&mut digest, pack.manifest.normalization.as_bytes());
        let mut families = pack
            .manifest
            .families
            .iter()
            .map(|family| {
                family["fingerprint"]
                    .as_str()
                    .ok_or_else(|| failure("family fingerprint is absent"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        families.sort_unstable();
        digest.update((families.len() as u64).to_le_bytes());
        for family in families {
            field(&mut digest, family.as_bytes());
        }
        let mut masters = pack
            .manifest
            .masters
            .iter()
            .map(|master| {
                let family = master["family_sha256"]
                    .as_str()
                    .ok_or_else(|| failure("master family identity is absent"))?;
                let powers = master["powers"]
                    .as_array()
                    .ok_or_else(|| failure("master powers are absent"))?
                    .iter()
                    .map(|power| {
                        power
                            .as_i64()
                            .ok_or_else(|| failure("invalid master power"))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok((family, powers))
            })
            .collect::<Result<Vec<_>, VakintError>>()?;
        masters.sort_unstable();
        digest.update((masters.len() as u64).to_le_bytes());
        for (family, powers) in masters {
            field(&mut digest, family.as_bytes());
            digest.update((powers.len() as u64).to_le_bytes());
            for power in powers {
                digest.update(power.to_le_bytes());
            }
        }
    }
    Ok(digest.finalize().into())
}
