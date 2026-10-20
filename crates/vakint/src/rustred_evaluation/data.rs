//! Validated immutable data snapshots shared by one Vakint engine.

use std::{
    collections::BTreeMap,
    fs,
    path::{Component, Path, PathBuf},
    sync::{Arc, LazyLock},
};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use super::{
    RustRedEvaluationError,
    artifact::{ArtifactFamily, ArtifactSnapshot},
};
use crate::{Vakint, VakintError};

mod vendored {
    include!(concat!(env!("OUT_DIR"), "/vakint_vendored_data.rs"));
}

pub(super) const RELEASE_SCHEMA: &str = "vakint.external-data.v1";
const LOOP_SCHEMA: &str = "vakint.external-data.loop.v1";

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(super) struct ReleaseManifest {
    pub schema: String,
    pub version: String,
    pub loop_packs: Vec<LoopPack>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(super) struct LoopPack {
    pub loop_count: usize,
    pub manifest: String,
    pub manifest_sha256: String,
    pub archive: String,
    pub archive_sha256: Option<String>,
    pub archive_bytes: Option<u64>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(super) struct FileDigest {
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Clone, Debug, Deserialize)]
pub(super) struct LoopManifest {
    pub schema: String,
    pub version: String,
    pub loop_count: usize,
    pub normalization: String,
    pub dimension: String,
    pub families: Vec<serde_json::Value>,
    pub masters: Vec<serde_json::Value>,
    pub cache: CacheManifest,
    pub coverage: String,
    pub files: BTreeMap<String, FileDigest>,
}

#[derive(Clone, Debug, Deserialize)]
pub(super) struct CacheManifest {
    pub available: bool,
    pub path: Option<String>,
}

#[derive(Clone, Debug)]
enum Source {
    Vendored,
    Directory(PathBuf),
    Override {
        original: Arc<Source>,
        path: String,
        bytes: Arc<[u8]>,
    },
}

#[derive(Clone, Debug)]
pub(super) struct LoopData {
    pub manifest: LoopManifest,
    source: Source,
}

pub(super) fn failure(detail: impl std::fmt::Display) -> RustRedEvaluationError {
    RustRedEvaluationError::ExternalData {
        detail: detail.to_string(),
    }
}

pub(super) fn relative_path(path: &str) -> Result<&Path, RustRedEvaluationError> {
    let value = Path::new(path);
    if path.is_empty()
        || path.contains('\\')
        || !value
            .components()
            .all(|p| matches!(p, Component::Normal(_)))
    {
        return Err(failure(format!("invalid package-relative path {path:?}")));
    }
    Ok(value)
}

pub(super) fn checked_bytes(bytes: &[u8], expected: &str) -> Result<(), RustRedEvaluationError> {
    let actual = format!("{:x}", Sha256::digest(bytes));
    if actual != expected {
        return Err(failure(format!(
            "SHA256 mismatch: expected {expected}, found {actual}"
        )));
    }
    Ok(())
}

pub(super) fn validate_loop_directory(
    root: &Path,
    pack: &LoopPack,
    version: &str,
) -> Result<(), RustRedEvaluationError> {
    LoopData::load(Source::Directory(root.to_owned()), pack, version).map(|_| ())
}

impl Source {
    fn read(&self, path: &str) -> Result<Vec<u8>, RustRedEvaluationError> {
        let relative = relative_path(path)?;
        match self {
            Self::Override {
                original,
                path: overridden,
                bytes,
            } => {
                if path == overridden {
                    Ok(bytes.to_vec())
                } else {
                    original.read(path)
                }
            }
            Self::Vendored => vendored::FILES
                .iter()
                .find(|(name, _)| *name == path)
                .map(|(_, bytes)| bytes.to_vec())
                .ok_or_else(|| failure(format!("missing vendored data file {path}"))),
            Self::Directory(root) => {
                // Symlinks would allow a developer checkout to escape the same
                // file boundary that the release extractor enforces.
                let mut current = root.clone();
                for part in relative.components() {
                    current.push(part.as_os_str());
                    if fs::symlink_metadata(&current)
                        .map_err(failure)?
                        .file_type()
                        .is_symlink()
                    {
                        return Err(failure(format!(
                            "symlink in package path {}",
                            current.display()
                        )));
                    }
                }
                fs::read(current).map_err(failure)
            }
        }
    }
}

impl LoopData {
    fn load(
        source: Source,
        pack: &LoopPack,
        version: &str,
    ) -> Result<Self, RustRedEvaluationError> {
        let bytes = source.read(&pack.manifest)?;
        checked_bytes(&bytes, &pack.manifest_sha256)?;
        let manifest: LoopManifest = serde_json::from_slice(&bytes).map_err(failure)?;
        if manifest.schema != LOOP_SCHEMA
            || manifest.loop_count != pack.loop_count
            || manifest.version != version
        {
            return Err(failure(format!(
                "inconsistent {}L package manifest",
                pack.loop_count
            )));
        }
        if manifest.normalization != "FMFTandMATAD" || manifest.dimension != "4-2*epsilon" {
            return Err(failure(
                "unsupported master-data normalization or dimension",
            ));
        }
        let data = Self { manifest, source };
        for path in data.manifest.files.keys() {
            data.read(path)?;
        }
        let coverage: serde_json::Value =
            serde_json::from_slice(&data.read(&data.manifest.coverage)?).map_err(failure)?;
        if coverage["schema"].as_str() != Some("vakint.external-data.coverage.v1")
            || coverage["loop_count"].as_u64() != Some(pack.loop_count as u64)
        {
            return Err(failure("incompatible campaign coverage manifest"));
        }
        Ok(data)
    }

    pub(super) fn read(&self, path: &str) -> Result<Vec<u8>, RustRedEvaluationError> {
        relative_path(path)?;
        let digest = self.manifest.files.get(path).ok_or_else(|| {
            failure(format!(
                "unlisted {}L data file {path}",
                self.manifest.loop_count
            ))
        })?;
        let bytes = self
            .source
            .read(&format!("{}L/{path}", self.manifest.loop_count))?;
        if bytes.len() as u64 != digest.bytes {
            return Err(failure(format!(
                "incorrect size of {}L/{path}",
                self.manifest.loop_count
            )));
        }
        checked_bytes(&bytes, &digest.sha256)?;
        Ok(bytes)
    }

    pub(super) fn path(&self, relative: &str) -> Result<Option<PathBuf>, RustRedEvaluationError> {
        let relative = relative_path(relative)?;
        Ok(match &self.source {
            Source::Vendored => None,
            Source::Directory(root) => Some(
                root.join(format!("{}L", self.manifest.loop_count))
                    .join(relative),
            ),
            Source::Override { .. } => None,
        })
    }
}

#[derive(Clone)]
pub(crate) struct RustRedData {
    pub(super) version: String,
    pub(super) loops: BTreeMap<usize, Arc<LoopData>>,
    lower: BTreeMap<usize, ArtifactSnapshot>,
    pub(super) four_loop: Option<Arc<super::four_loop::FourLoopPrograms>>,
    pub(super) masters: Arc<super::masters::MasterCatalog>,
}

impl std::fmt::Debug for RustRedData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RustRedData")
            .field("version", &self.version)
            .field("loops", &self.loops.keys())
            .finish_non_exhaustive()
    }
}

impl RustRedData {
    fn load(
        source: Source,
        manifest: ReleaseManifest,
        max_loop_count: usize,
    ) -> Result<Arc<Self>, RustRedEvaluationError> {
        if manifest.schema != RELEASE_SCHEMA || !(1..=5).contains(&max_loop_count) {
            return Err(failure(
                "unsupported external-data schema or loop count (expected 1 through 5)",
            ));
        }
        let mut loops = BTreeMap::new();
        for count in 1..=max_loop_count {
            let mut matches = manifest
                .loop_packs
                .iter()
                .filter(|pack| pack.loop_count == count);
            let pack = matches.next().ok_or_else(|| {
                failure(format!("release {} has no {count}L pack", manifest.version))
            })?;
            if matches.next().is_some() || pack.manifest != format!("{count}L/manifest.json") {
                return Err(failure(format!("ambiguous {count}L pack manifest")));
            }
            loops.insert(
                count,
                Arc::new(LoopData::load(source.clone(), pack, &manifest.version)?),
            );
        }
        // Vendored lower-loop defaults remain available after a 1L/2L load.
        if max_loop_count < 3 {
            for (&count, data) in &Self::vendored()?.loops {
                loops.entry(count).or_insert_with(|| Arc::clone(data));
            }
        }
        let mut lower = BTreeMap::new();
        for (&count, data) in loops.range(..=3) {
            let family = super::artifact::shipped_family_for_loop_count(count).unwrap();
            lower.insert(
                count,
                family.load_snapshot(&data.read("ibps/rules.rrbin")?)?,
            );
        }
        let four_loop = loops
            .get(&4)
            .map(|data| super::four_loop::FourLoopPrograms::load(data).map(Arc::new))
            .transpose()?;
        let masters = Arc::new(super::masters::MasterCatalog::load(
            &loops,
            &manifest.version,
        )?);
        Ok(Arc::new(Self {
            version: manifest.version,
            loops,
            lower,
            four_loop,
            masters,
        }))
    }

    pub(super) fn from_directory(
        path: &Path,
        version: &str,
        max_loop_count: usize,
    ) -> Result<Arc<Self>, RustRedEvaluationError> {
        let source = Source::Directory(path.canonicalize().map_err(failure)?);
        let manifest: ReleaseManifest =
            serde_json::from_slice(&source.read("manifest.json")?).map_err(failure)?;
        if version != "latest" && manifest.version != version {
            return Err(failure(format!(
                "requested version {version}, but local data declares {}",
                manifest.version
            )));
        }
        Self::load(source, manifest, max_loop_count)
    }

    pub(super) fn vendored() -> Result<Arc<Self>, RustRedEvaluationError> {
        static DATA: LazyLock<Result<Arc<RustRedData>, String>> = LazyLock::new(|| {
            let source = Source::Vendored;
            let manifest =
                serde_json::from_slice(&source.read("manifest.json").map_err(|e| e.to_string())?)
                    .map_err(|e| e.to_string())?;
            RustRedData::load(source, manifest, 3).map_err(|e| e.to_string())
        });
        DATA.as_ref().cloned().map_err(failure)
    }

    pub(super) fn artifact(
        &self,
        family: ArtifactFamily,
    ) -> Result<ArtifactSnapshot, RustRedEvaluationError> {
        self.lower
            .get(&family.loop_count())
            .cloned()
            .ok_or_else(|| failure("missing lower-loop IBP program"))
    }
}

impl Vakint {
    pub(crate) fn resolve_rustred_numerical(
        &self,
        expression: symbolica::atom::Atom,
        settings: &crate::VakintSettings,
    ) -> Result<symbolica::atom::Atom, VakintError> {
        self.rustred_data()?
            .masters
            .resolve_numerical(expression, settings)
    }

    pub(crate) fn has_rustred_loop(&self, count: usize) -> bool {
        self.rustred_data
            .as_ref()
            .is_some_and(|data| data.loops.contains_key(&count))
    }

    pub(crate) fn rustred_data(&self) -> Result<Arc<RustRedData>, RustRedEvaluationError> {
        self.rustred_data
            .clone()
            .map(Ok)
            .unwrap_or_else(RustRedData::vendored)
    }

    /// Load all loop packs up to the requested maximum from an existing checkout.
    /// Validation and cache attachment finish before this engine is replaced.
    pub fn load_external_data(
        &mut self,
        path: impl AsRef<Path>,
        max_loop_count: usize,
    ) -> Result<(), VakintError> {
        let data = RustRedData::from_directory(path.as_ref(), "latest", max_loop_count)?;
        self.activate_external_data(data)
    }

    /// Load native IBPs from a local pack or a sealed 1L–3L rules file. A bare
    /// closed program retains the active master catalog and drops that loop's
    /// prepared cache, whose identity belonged to the preceding program.
    pub fn load_rustred_artifacts(
        &mut self,
        loop_count: usize,
        path: impl AsRef<Path>,
    ) -> Result<(), VakintError> {
        let path = path.as_ref();
        for root in path.ancestors().take(4) {
            if root.join("manifest.json").is_file() {
                let bytes = fs::read(root.join("manifest.json")).map_err(failure)?;
                if serde_json::from_slice::<ReleaseManifest>(&bytes).is_ok() {
                    return self.load_external_data(root, loop_count);
                }
            }
        }
        let family = super::artifact::shipped_family_for_loop_count(loop_count)
            .ok_or_else(|| failure("4L/5L IBPs require a local external-data pack; the explicit saved-application API also accepts a bare 5L application"))?;
        let path = if path.is_dir() {
            path.join("rules.rrbin")
        } else {
            path.to_owned()
        };
        let bytes = fs::read(&path).map_err(failure)?;
        let snapshot = family.load_snapshot(&bytes)?;
        let mut data = (*self.rustred_data()?).clone();
        let mut pack = (*data.loops[&loop_count]).clone();
        let digest = format!("{:x}", Sha256::digest(&bytes));
        pack.manifest.files.insert(
            "ibps/rules.rrbin".into(),
            FileDigest {
                sha256: digest.clone(),
                bytes: bytes.len() as u64,
            },
        );
        pack.manifest.cache = CacheManifest {
            available: false,
            path: None,
        };
        let original = match pack.source {
            Source::Override { original, .. } => original,
            source => Arc::new(source),
        };
        pack.source = Source::Override {
            original,
            path: format!("{loop_count}L/ibps/rules.rrbin"),
            bytes: bytes.into(),
        };
        data.loops.insert(loop_count, Arc::new(pack));
        data.lower.insert(loop_count, snapshot);
        data.version = format!("local-{}", &digest[..16]);
        self.activate_external_data(Arc::new(data))
    }

    pub(super) fn activate_external_data(
        &mut self,
        data: Arc<RustRedData>,
    ) -> Result<(), VakintError> {
        let mut staged = self.clone();
        staged.rustred_data = Some(Arc::clone(&data));
        staged.rustred_lower_caches = Arc::default();
        staged.rustred_application = Arc::default();
        for (&count, pack) in &data.loops {
            if pack.manifest.cache.available {
                let relative = pack
                    .manifest
                    .cache
                    .path
                    .as_deref()
                    .ok_or_else(|| failure("available cache has no path"))?;
                let path = pack
                    .path(relative)?
                    .ok_or_else(|| failure("vendored data must not contain caches"))?;
                staged.load_cache(count, path, &Default::default())?;
            }
        }
        if data.loops.contains_key(&5) {
            staged.saved_application(&crate::VakintSettings::default())?;
        }
        *self = staged;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RustRedMasterMode, VakintSettings, vakint_parse};

    #[test]
    #[ignore = "requires the authenticated v0.1.0 5L data pack; no IBP application is loaded"]
    fn authentic_five_loop_declared_terminals_remain_uncomputed() {
        Vakint::initialize_vakint_symbols();
        let source = Source::Directory(
            super::super::cache_tests::external_path()
                .canonicalize()
                .unwrap(),
        );
        let release: ReleaseManifest =
            serde_json::from_slice(&source.read("manifest.json").unwrap()).unwrap();
        assert_eq!(release.schema, RELEASE_SCHEMA);
        assert_eq!(release.version, "v0.1.0");
        let pack = release
            .loop_packs
            .iter()
            .find(|pack| pack.loop_count == 5)
            .unwrap();
        assert_eq!(
            pack.manifest_sha256,
            "44d047c2bd280e24ecda4730096bdc0fa449afd0f6fcba90f1002d886ed06f51"
        );
        // Authenticate the real 5L pack and its payloads without constructing
        // an IBP application, attaching caches, or loading lower-loop packs.
        let data = Arc::new(LoopData::load(source, pack, &release.version).unwrap());
        assert_eq!(data.manifest.masters.len(), 905);
        for terminal in &data.manifest.masters {
            assert_eq!(terminal["kind"], "uncomputed-terminal");
            assert_eq!(terminal["status"], "not_computed");
            assert!(terminal["coefficients"].as_array().unwrap().is_empty());
            assert!(
                terminal["available_epsilon_powers"]
                    .as_array()
                    .unwrap()
                    .is_empty()
            );
        }
        let terminal = &data.manifest.masters[0];
        let identity = terminal["master_id"].as_str().unwrap();
        assert_eq!(
            identity,
            "L5_44b030daa664fca0c71f4bb930624b098348f6184bfe46ca3ad83b86ef71df3b_m1_m1_m1_p0_p1_p1_p0_p0_p1_p0_p0_p1_p0_p1_p1"
        );
        let powers: Vec<i64> = serde_json::from_value(terminal["powers"].clone()).unwrap();
        assert_eq!(powers, [-1, -1, -1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1]);
        let family = data
            .manifest
            .families
            .iter()
            .find(|family| family["sha256"] == terminal["family_sha256"])
            .unwrap();
        let graph = data
            .path(terminal["parent_graph"].as_str().unwrap())
            .unwrap()
            .unwrap();
        assert!(graph.is_file());
        assert_eq!(graph.extension().unwrap(), "dot");
        let loops = BTreeMap::from([(5, Arc::clone(&data))]);
        let catalog = super::super::masters::MasterCatalog::load(&loops, &release.version).unwrap();
        let master = catalog.symbol(family["fingerprint"].as_str().unwrap(), &powers);
        let settings = VakintSettings {
            epsilon_symbol: "ep".into(),
            number_of_terms_in_epsilon_expansion: 6,
            run_time_decimal_precision: 32,
            ..Default::default()
        };
        let retained = catalog
            .finalize(master.clone(), &settings, 5, RustRedMasterMode::Retain)
            .unwrap();
        assert_eq!(retained, master);
        let direct_error = catalog
            .finalize(master.clone(), &settings, 5, RustRedMasterMode::Numerical)
            .unwrap_err()
            .to_string();
        let retained_error = catalog
            .resolve_numerical(retained, &settings)
            .unwrap_err()
            .to_string();
        for error in [&direct_error, &retained_error] {
            assert!(error.contains(identity), "{error}");
            assert!(
                error.contains("declared terminal has no computed coefficients"),
                "{error}"
            );
            assert!(error.contains(&release.version), "{error}");
            assert!(error.contains(graph.to_str().unwrap()), "{error}");
            assert!(!error.contains("undeclared residual master"), "{error}");
        }
        // The cancellation is meromorphic and requires exact combination of
        // rational epsilon coefficients before inspecting the unavailable data.
        let prefactor = vakint_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
        let cancelling = prefactor * master;
        assert!(!cancelling.is_zero());
        assert!(
            catalog
                .finalize(
                    cancelling.clone(),
                    &settings,
                    5,
                    RustRedMasterMode::Numerical
                )
                .unwrap()
                .is_zero()
        );
        assert!(
            catalog
                .resolve_numerical(cancelling, &settings)
                .unwrap()
                .is_zero()
        );
        eprintln!(
            "Authenticated {} declared/uncomputed terminals in {}; selected {identity}, powers {powers:?}, graph {}; Retain and both exact-cancellation routes passed. Direct error: {direct_error}; retained-catalog error: {retained_error}",
            data.manifest.masters.len(),
            release.version,
            graph.display(),
        );
    }
}
