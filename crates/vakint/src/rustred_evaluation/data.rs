//! Validated immutable data snapshots shared by one Vakint engine.

use std::{
    collections::BTreeMap,
    fs,
    path::{Component, Path, PathBuf},
    sync::Arc,
};

use once_cell::sync::OnceCell;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use symbolica::utils::resources;

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
pub(crate) struct FileDigest {
    pub sha256: String,
    pub bytes: u64,
}

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct FileReference {
    pub path: String,
    pub sha256: String,
}

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct LoopManifest {
    pub schema: String,
    pub version: String,
    pub loop_count: usize,
    pub normalization: String,
    pub dimension: String,
    pub families: Vec<serde_json::Value>,
    pub masters: Vec<serde_json::Value>,
    #[serde(default)]
    pub constant_bindings: Option<FileReference>,
    pub cache: CacheManifest,
    pub coverage: String,
    pub files: BTreeMap<String, FileDigest>,
}

#[derive(Clone, Debug, Deserialize)]
pub(crate) struct CacheManifest {
    pub available: bool,
    pub path: Option<String>,
    #[serde(default)]
    pub format: Option<String>,
    #[serde(default)]
    pub identity: Option<rustred_app::SavedApplicationCacheFileIdentity>,
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
pub(crate) struct LoopData {
    pub manifest: LoopManifest,
    source: Source,
    manifest_sha256: String,
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
        resources::checkpoint()?;
        let relative = relative_path(path)?;
        match self {
            Self::Override {
                original,
                path: overridden,
                bytes,
            } => {
                if path == overridden {
                    resources::check_allocation(bytes.len().try_into().map_err(failure)?)?;
                    Ok(bytes.to_vec())
                } else {
                    original.read(path)
                }
            }
            Self::Vendored => {
                let (_, bytes) = vendored::FILES
                    .iter()
                    .find(|(name, _)| *name == path)
                    .ok_or_else(|| failure(format!("missing vendored data file {path}")))?;
                resources::check_allocation(bytes.len().try_into().map_err(failure)?)?;
                Ok(bytes.to_vec())
            }
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
                resources::check_allocation(fs::metadata(&current).map_err(failure)?.len())?;
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
        let data = Self {
            manifest,
            source,
            manifest_sha256: pack.manifest_sha256.clone(),
        };
        for path in data.manifest.files.keys() {
            resources::checkpoint()?;
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

    pub(crate) fn read(&self, path: &str) -> Result<Vec<u8>, RustRedEvaluationError> {
        relative_path(path)?;
        let digest = self.manifest.files.get(path).ok_or_else(|| {
            failure(format!(
                "unlisted {}L data file {path}",
                self.manifest.loop_count
            ))
        })?;
        resources::check_allocation(digest.bytes)?;
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

    pub(crate) fn path(&self, relative: &str) -> Result<Option<PathBuf>, RustRedEvaluationError> {
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
    pub(crate) masters: Arc<super::masters::MasterCatalog>,
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
    /// Recognize a checkout or an installed pack before considering a bare IBP path.
    pub(super) fn release_root(path: &Path) -> Result<Option<PathBuf>, RustRedEvaluationError> {
        for root in path.ancestors().take(4) {
            let manifest = root.join("manifest.json");
            if manifest.is_file() {
                let bytes = fs::read(manifest).map_err(failure)?;
                if serde_json::from_slice::<ReleaseManifest>(&bytes).is_ok() {
                    return Ok(Some(root.to_owned()));
                }
            }
        }
        Ok(None)
    }

    fn load_loops(
        source: &Source,
        manifest: &ReleaseManifest,
        max_loop_count: usize,
    ) -> Result<BTreeMap<usize, Arc<LoopData>>, RustRedEvaluationError> {
        resources::checkpoint()?;
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
        Ok(loops)
    }

    fn load(
        source: Source,
        manifest: ReleaseManifest,
        max_loop_count: usize,
    ) -> Result<Arc<Self>, RustRedEvaluationError> {
        let loops = Self::load_loops(&source, &manifest, max_loop_count)?;
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
        let data = Arc::new(Self {
            version: manifest.version,
            loops,
            lower,
            four_loop,
            masters,
        });
        resources::checkpoint()?;
        Ok(data)
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

    /// Recognize only the already authenticated pack's exact IBP path.
    pub(super) fn directory_root_for_ibps(
        &self,
        loop_count: usize,
        path: &Path,
    ) -> Option<PathBuf> {
        let pack = self.loops.get(&loop_count)?;
        let Source::Directory(root) = &pack.source else {
            return None;
        };
        let expected = root.join(format!("{loop_count}L/ibps"));
        if path == expected {
            return Some(root.clone());
        }
        // Resolve parent aliases without following a replaced terminal symlink.
        // The ordinary pack validator must reject that symlink, not classify
        // its target as an unrelated bare application.
        let parent = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty());
        let candidate = parent
            .unwrap_or_else(|| Path::new("."))
            .canonicalize()
            .ok()?
            .join(path.file_name()?);
        (candidate == expected).then(|| root.clone())
    }

    /// Authenticate a mutable checkout before preserving this snapshot's native owners.
    pub(super) fn matches_directory(
        &self,
        path: &Path,
        max_loop_count: usize,
    ) -> Result<bool, RustRedEvaluationError> {
        let root = path.canonicalize().map_err(failure)?;
        if self.loops.len() != max_loop_count
            || self.loops.values().any(
                |pack| !matches!(&pack.source, Source::Directory(directory) if directory == &root),
            )
        {
            return Ok(false);
        }
        let source = Source::Directory(root);
        let manifest: ReleaseManifest =
            serde_json::from_slice(&source.read("manifest.json")?).map_err(failure)?;
        // This is the same complete schema/checksum/coverage validation used by
        // ordinary loading. It deliberately precedes an identity reuse decision.
        let loops = Self::load_loops(&source, &manifest, max_loop_count)?;
        Ok(self.version == manifest.version
            && loops.iter().all(|(count, pack)| {
                self.loops
                    .get(count)
                    .is_some_and(|current| current.manifest_sha256 == pack.manifest_sha256)
            }))
    }

    pub(super) fn vendored() -> Result<Arc<Self>, RustRedEvaluationError> {
        static DATA: OnceCell<Arc<RustRedData>> = OnceCell::new();
        Self::initialize_vendored(&DATA)
    }

    fn initialize_vendored(
        cell: &OnceCell<Arc<Self>>,
    ) -> Result<Arc<Self>, RustRedEvaluationError> {
        // A constrained first caller must not poison offline construction for
        // every later engine. Only a fully validated snapshot is shared.
        resources::checkpoint()?;
        let data = cell.get_or_try_init(|| {
            let source = Source::Vendored;
            let manifest =
                serde_json::from_slice(&source.read("manifest.json")?).map_err(failure)?;
            Self::load(source, manifest, 3)
        })?;
        resources::checkpoint()?;
        Ok(Arc::clone(data))
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
        self.require_bare_master_coefficients(&expression)?;
        let data = self.rustred_data()?;
        crate::precision::materialize(
            expression.as_view(),
            settings,
            &std::collections::HashMap::new(),
            Some(&data.masters),
        )
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
        self.activate_external_data(data, true)
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
        if let Some(root) = RustRedData::release_root(path)? {
            return self.load_external_data(root, loop_count);
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
            format: None,
            identity: None,
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
        self.activate_external_data(Arc::new(data), true)
    }

    pub(super) fn activate_external_data(
        &mut self,
        data: Arc<RustRedData>,
        load_prepared_caches: bool,
    ) -> Result<(), VakintError> {
        let mut staged = self.clone();
        staged.rustred_data = Some(Arc::clone(&data));
        staged.rustred_lower_caches = Arc::default();
        staged.rustred_application = Arc::default();
        for (&count, pack) in &data.loops {
            if load_prepared_caches && pack.manifest.cache.available {
                let relative = pack
                    .manifest
                    .cache
                    .path
                    .as_deref()
                    .ok_or_else(|| failure("available cache has no path"))?;
                let path = pack
                    .path(relative)?
                    .ok_or_else(|| failure("vendored data must not contain caches"))?;
                let mut options = super::cache::RustRedCacheLoadOptions::default();
                if count == 5 {
                    if pack.manifest.cache.format.as_deref()
                        != Some("rustred.saved-application-cache-file.v1")
                    {
                        return Err(failure("unsupported five-loop prepared-cache format").into());
                    }
                    options.file_identity =
                        Some(pack.manifest.cache.identity.ok_or_else(|| {
                            failure("five-loop prepared cache has no authenticated file identity")
                        })?);
                }
                staged.load_cache(count, path, &options)?;
            }
        }
        if data.loops.contains_key(&5) {
            staged.saved_application(&crate::VakintSettings::default())?;
        }
        resources::checkpoint()?;
        *self = staged;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{RustRedMasterMode, VakintSettings, vakint_parse};
    use std::collections::BTreeSet;

    fn vendored_directory() -> tempfile::TempDir {
        let directory = tempfile::tempdir().unwrap();
        for (relative, bytes) in vendored::FILES {
            let path = directory.path().join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, bytes).unwrap();
        }
        directory
    }

    #[test]
    fn known_directory_ibps_path_still_requires_valid_root_manifest() {
        let directory = vendored_directory();
        let data = RustRedData::from_directory(directory.path(), "latest", 3).unwrap();
        let ibps = directory.path().join("3L/ibps").canonicalize().unwrap();
        let root = data.directory_root_for_ibps(3, &ibps).unwrap();
        #[cfg(unix)]
        let root_manifest = fs::read(root.join("manifest.json")).unwrap();
        assert_eq!(root, directory.path().canonicalize().unwrap());
        assert!(
            data.directory_root_for_ibps(3, &directory.path().join("1L/masters"))
                .is_none()
        );
        assert!(
            RustRedData::vendored()
                .unwrap()
                .directory_root_for_ibps(3, &ibps)
                .is_none()
        );
        fs::remove_file(root.join("manifest.json")).unwrap();
        assert_eq!(data.directory_root_for_ibps(3, &ibps), Some(root.clone()));
        assert!(data.matches_directory(&root, 3).is_err());
        #[cfg(unix)]
        {
            fs::write(root.join("manifest.json"), &root_manifest).unwrap();
            let outside = tempfile::tempdir().unwrap();
            fs::rename(&ibps, directory.path().join("3L/original-ibps")).unwrap();
            std::os::unix::fs::symlink(outside.path(), &ibps).unwrap();
            assert_eq!(
                RustRedData::release_root(&ibps).unwrap(),
                Some(root.clone())
            );
            assert_eq!(data.directory_root_for_ibps(3, &ibps), Some(root.clone()));
            let parent_alias = directory.path().join("three-loop-alias");
            std::os::unix::fs::symlink(ibps.parent().unwrap(), &parent_alias).unwrap();
            let relative_terminal = parent_alias.join("ibps");
            assert_eq!(
                data.directory_root_for_ibps(3, &relative_terminal),
                Some(root.clone())
            );
            assert!(
                data.directory_root_for_ibps(3, &ibps.canonicalize().unwrap())
                    .is_none()
            );
            assert!(
                data.matches_directory(&root, 3)
                    .unwrap_err()
                    .to_string()
                    .contains("symlink")
            );
            fs::remove_file(root.join("manifest.json")).unwrap();
            assert!(RustRedData::release_root(&ibps).unwrap().is_none());
            assert_eq!(data.directory_root_for_ibps(3, &ibps), Some(root.clone()));
            assert!(data.matches_directory(&root, 3).is_err());
        }
        fs::write(root.join("manifest.json"), b"{broken-json").unwrap();
        assert_eq!(data.directory_root_for_ibps(3, &ibps), Some(root.clone()));
        assert!(data.matches_directory(&root, 3).is_err());
    }

    #[test]
    fn directory_reuse_checks_content_identity_and_rejects_corruption() {
        let directory = vendored_directory();
        let data = RustRedData::from_directory(directory.path(), "latest", 3).unwrap();
        let mut engine = Vakint::new().unwrap();
        engine.rustred_data = Some(Arc::clone(&data));
        let caches = Arc::clone(&engine.rustred_lower_caches);
        let application = Arc::clone(&engine.rustred_application);
        assert!(data.matches_directory(directory.path(), 3).unwrap());
        assert!(
            !RustRedData::vendored()
                .unwrap()
                .matches_directory(directory.path(), 3)
                .unwrap()
        );
        let different_directory = tempfile::tempdir().unwrap();
        assert!(
            !data
                .matches_directory(different_directory.path(), 3)
                .unwrap()
        );

        let coverage_path = directory
            .path()
            .join("1L")
            .join(&data.loops[&1].manifest.coverage);
        let coverage = fs::read(&coverage_path).unwrap();
        fs::write(&coverage_path, b"corrupt payload").unwrap();
        assert!(data.matches_directory(directory.path(), 3).is_err());
        assert!(engine.load_external_data(directory.path(), 3).is_err());
        assert!(Arc::ptr_eq(&data, &engine.rustred_data().unwrap()));
        assert!(Arc::ptr_eq(&caches, &engine.rustred_lower_caches));
        assert!(Arc::ptr_eq(&application, &engine.rustred_application));
        fs::write(&coverage_path, &coverage).unwrap();
        assert!(data.matches_directory(directory.path(), 3).unwrap());
        #[cfg(unix)]
        {
            let outside = tempfile::tempdir().unwrap();
            let outside_path = outside.path().join("coverage.json");
            fs::write(&outside_path, &coverage).unwrap();
            fs::remove_file(&coverage_path).unwrap();
            std::os::unix::fs::symlink(&outside_path, &coverage_path).unwrap();
            assert!(data.matches_directory(directory.path(), 3).is_err());
            fs::remove_file(&coverage_path).unwrap();
            fs::write(&coverage_path, &coverage).unwrap();
        }

        let manifest_path = directory.path().join("1L/manifest.json");
        let mut manifest = fs::read(&manifest_path).unwrap();
        manifest.push(b'\n');
        fs::write(&manifest_path, &manifest).unwrap();
        let release_path = directory.path().join("manifest.json");
        let mut release: ReleaseManifest =
            serde_json::from_slice(&fs::read(&release_path).unwrap()).unwrap();
        release
            .loop_packs
            .iter_mut()
            .find(|pack| pack.loop_count == 1)
            .unwrap()
            .manifest_sha256 = format!("{:x}", Sha256::digest(&manifest));
        fs::write(&release_path, serde_json::to_vec(&release).unwrap()).unwrap();
        // A valid new raw manifest is a different snapshot even when its parsed
        // mathematical payload happens to be unchanged.
        assert!(!data.matches_directory(directory.path(), 3).unwrap());

        let mut overridden = (*data).clone();
        let mut pack = (*overridden.loops[&1]).clone();
        pack.source = Source::Override {
            original: Arc::new(pack.source),
            path: format!("1L/{}", pack.manifest.coverage),
            bytes: coverage.into(),
        };
        overridden.loops.insert(1, Arc::new(pack));
        assert!(!overridden.matches_directory(directory.path(), 3).unwrap());
    }

    #[test]
    fn failed_vendored_initialization_remains_retryable_and_warm_hits_check_resources() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        use resources::{ResourceBudget, ResourceLimitError, ResourceLimits};

        let cell = OnceCell::new();
        let calls = AtomicUsize::new(0);
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1),
                ..Default::default()
            },
            move || {
                Ok(if calls.fetch_add(1, Ordering::Relaxed) == 0 {
                    0
                } else {
                    2
                })
            },
        );
        assert!(matches!(
            budget.run(|| RustRedData::initialize_vendored(&cell)),
            Err(RustRedEvaluationError::ResourceLimit(
                ResourceLimitError::Memory { .. }
            ))
        ));
        assert!(cell.get().is_none());

        let first = RustRedData::initialize_vendored(&cell).unwrap();
        let second = RustRedData::initialize_vendored(&cell).unwrap();
        assert!(Arc::ptr_eq(&first, &second));
        assert_eq!(first.loops.keys().copied().collect::<Vec<_>>(), [1, 2, 3]);

        let expired = ResourceBudget::new(ResourceLimits {
            max_elapsed: Some(std::time::Duration::ZERO),
            ..Default::default()
        });
        assert!(matches!(
            expired.run(|| RustRedData::initialize_vendored(&cell)),
            Err(RustRedEvaluationError::ResourceLimit(
                ResourceLimitError::Elapsed { .. }
            ))
        ));
        assert!(Arc::ptr_eq(cell.get().unwrap(), &first));
    }

    #[test]
    fn finite_constant_payloads_are_validated_and_participate_in_snapshot_identity() {
        use symbolica::atom::AtomCore;
        let original = RustRedData::vendored().unwrap();
        let valid = serde_json::json!({
            "symbol": "vakint::K5Moment(1)", "value": "1.25",
            "significant_digits": 30, "absolute_error": "1e-31", "error_quality": "bounded",
            "definition": "certified test moment", "provenance": {"certificate": "test-proof"}
        });
        let load = |rows: Vec<serde_json::Value>, tamper: bool| {
            let bytes = serde_json::to_vec(&rows).unwrap();
            let digest = format!("{:x}", Sha256::digest(&bytes));
            let mut pack = (*original.loops[&1]).clone();
            pack.manifest.constant_bindings = Some(FileReference {
                path: "constants.json".into(),
                sha256: if tamper {
                    "0".repeat(64)
                } else {
                    digest.clone()
                },
            });
            pack.manifest.files.insert(
                "constants.json".into(),
                FileDigest {
                    sha256: digest,
                    bytes: bytes.len() as u64,
                },
            );
            pack.source = Source::Override {
                original: Arc::new(pack.source),
                path: "1L/constants.json".into(),
                bytes: bytes.into(),
            };
            super::super::masters::MasterCatalog::load(
                &BTreeMap::from([(1, Arc::new(pack))]),
                "test-constants",
            )
        };
        let first = load(vec![valid.clone()], false).unwrap();
        let symbol = vakint_parse!("K5Moment(1)").unwrap();
        assert_eq!(first.constants[&symbol].significant_digits, 30);
        let mut changed = valid.clone();
        changed["value"] = "1.5".into();
        let second = load(vec![changed], false).unwrap();
        assert_ne!(
            first.coefficient_snapshot_tag(),
            second.coefficient_snapshot_tag()
        );
        let retained = first.prepare_evaluation(symbol).unwrap().expression;
        assert!(retained.contains_symbol(crate::vakint_symbol!("RustRedConstant")));
        assert!(
            second
                .prepare_evaluation(retained)
                .unwrap_err()
                .to_string()
                .contains("another data snapshot")
        );
        assert!(load(vec![valid.clone()], true).is_err());
        assert!(load(vec![valid.clone(), valid.clone()], false).is_err());
        for (field, value) in [
            ("significant_digits", serde_json::json!(0)),
            ("absolute_error", serde_json::json!("-1")),
            ("value", serde_json::json!("NaN")),
            ("definition", serde_json::json!("")),
            ("provenance", serde_json::json!({})),
            ("symbol", serde_json::json!("z3")),
            ("symbol", serde_json::json!("PR4fin")),
            ("symbol", serde_json::json!("ep")),
            ("symbol", serde_json::json!("p(1)")),
            ("symbol", serde_json::json!("NumericalSource()")),
            ("symbol", serde_json::json!("symbolica::log(2)")),
        ] {
            let mut invalid = valid.clone();
            invalid[field] = value;
            assert!(load(vec![invalid], false).is_err(), "{field}");
        }
    }

    #[test]
    fn vendored_master_inspection_returns_native_graph_and_all_orders() {
        let engine = Vakint::new().unwrap();
        let data = engine.rustred_data().unwrap();
        let metadata = &data.loops[&1].manifest.masters[0];
        let symbol = metadata["master_id"].as_str().unwrap().parse().unwrap();
        let master = engine.get_master(&symbol, 20).unwrap();
        assert_eq!(master.loop_count(), 1);
        assert_eq!(
            master.coefficients().keys().copied().collect::<Vec<_>>(),
            vec![-1, 0, 1, 2, 3, 4]
        );
        assert_eq!(master.requested_decimal_precision(), 20);
        master.graph().validate().unwrap();
        let dot = master.graph().to_dot().unwrap();
        assert!(dot.contains("lmb_id=0"));
        let restored =
            feynkit_graph::FeynmanDiagram::from_dot(master.graph().model_arc(), &dot).unwrap();
        assert_eq!(restored.to_dot().unwrap(), dot);
        assert!(master.coefficients().values().all(|coefficient| matches!(
            coefficient.symbolic_expression(),
            crate::inspection::MasterSymbolicExpression::Expression(_)
        )));
    }

    #[test]
    fn external_inspection_preserves_dots_and_uncomputed_numerators() {
        Vakint::initialize_vakint_symbols();
        let source = Source::Directory(super::super::cache_tests::external_path_through(5));
        let release: ReleaseManifest =
            serde_json::from_slice(&source.read("manifest.json").unwrap()).unwrap();
        // Graph lookup does not need to construct a reduction application or
        // attach its caches. Authenticate the same pack boundary directly.
        for (loop_count, name) in [(4, "P9"), (5, "P1")] {
            let pack = release
                .loop_packs
                .iter()
                .find(|pack| pack.loop_count == loop_count)
                .unwrap();
            let data = Arc::new(LoopData::load(source.clone(), pack, &release.version).unwrap());
            let metadata = data
                .manifest
                .masters
                .iter()
                .find(|master| master["name"] == name)
                .unwrap();
            let symbol = metadata["master_id"].as_str().unwrap().parse().unwrap();
            let catalog = Arc::new(
                super::super::masters::MasterCatalog::load(
                    &BTreeMap::from([(loop_count, Arc::clone(&data))]),
                    &release.version,
                )
                .unwrap(),
            );
            let master = catalog.get_master(symbol, 16).unwrap();
            master.graph().validate().unwrap();
            let dot = master.graph().to_dot().unwrap();
            let restored =
                feynkit_graph::FeynmanDiagram::from_dot(master.graph().model_arc(), &dot).unwrap();
            assert_eq!(restored.to_dot().unwrap(), dot);
            if loop_count == 4 {
                assert!(
                    master
                        .graph()
                        .edges()
                        .any(|(_, _, edge)| edge.propagator_power == 2)
                );
                assert!(master.is_computed());
            } else {
                assert!(!master.graph().numerator().expression().is_one());
                assert!(!master.is_computed());
                assert!(master.coefficients().is_empty());
                #[cfg(feature = "feynkit-ingress")]
                {
                    use symbolica::atom::AtomCore;

                    let dimension = symbolica::atom::Atom::var(feynkit_graph::symbols::dimension());
                    let kinematics =
                        feynkit_kinematics::Kinematics::in_dimension(&dimension).unwrap();
                    let family = master.graph().integral_family(&kinematics, &[]).unwrap();
                    let integral = crate::VakintExpression::from_diagram(
                        master.graph(),
                        &family,
                        master.graph().numerator().expression(),
                        &crate::diagram_integral::DiagramIntegralOptions::default(),
                    )
                    .unwrap();
                    let integral = symbolica::atom::Atom::from(integral);
                    assert!(!integral.is_zero());
                    assert!(integral.get_all_symbols(true).iter().all(|symbol| {
                        !symbol.get_name().starts_with("spenso::")
                            && symbol != &feynkit_graph::symbols::loop_momentum()
                            && symbol != &feynkit_graph::symbols::external_momentum()
                    }));
                }
            }
        }
    }

    #[test]
    #[ignore = "requires the authenticated 640-master 5L pack; no IBP application is loaded"]
    fn authentic_five_loop_declared_terminals_remain_uncomputed() {
        Vakint::initialize_vakint_symbols();
        let source = Source::Directory(
            super::super::cache_tests::external_path_through(5)
                .canonicalize()
                .unwrap(),
        );
        let release: ReleaseManifest =
            serde_json::from_slice(&source.read("manifest.json").unwrap()).unwrap();
        assert_eq!(release.schema, RELEASE_SCHEMA);
        let pack = release
            .loop_packs
            .iter()
            .find(|pack| pack.loop_count == 5)
            .unwrap();
        // Authenticate the real 5L pack and its payloads without constructing
        // an IBP application, attaching caches, or loading lower-loop packs.
        let data = Arc::new(LoopData::load(source, pack, &release.version).unwrap());
        assert_eq!(data.manifest.masters.len(), 640);
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
            "L5_44b030daa664fca0c71f4bb930624b098348f6184bfe46ca3ad83b86ef71df3b_m4_p1_p1_p0_p0_p1_p1_p1_m4_p1_p1_p0_p1_m1_p0"
        );
        let powers: Vec<i64> = serde_json::from_value(terminal["powers"].clone()).unwrap();
        assert_eq!(powers, [-4, 1, 1, 0, 0, 1, 1, 1, -4, 1, 1, 0, 1, -1, 0]);
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
        let symbolic_error = catalog
            .finalize(master.clone(), &settings, 5, RustRedMasterMode::Symbolic)
            .unwrap_err()
            .to_string();
        let retained_error = catalog
            .resolve_numerical(retained, &settings)
            .unwrap_err()
            .to_string();
        for error in [&direct_error, &symbolic_error, &retained_error] {
            assert!(error.contains(identity), "{error}");
            assert!(
                error.contains("declared terminal has no computed coefficients"),
                "{error}"
            );
            assert!(error.contains(&release.version), "{error}");
            assert!(error.contains(graph.to_str().unwrap()), "{error}");
            assert!(!error.contains("undeclared residual master"), "{error}");
        }
        // The common-family graph is a powered P640 representative, not a
        // certified parent for an arbitrary residual in that family. Reuse the
        // authenticated metadata's maximal-parent references without inferring
        // another routing from positive slots alone.
        let parent_references = data
            .manifest
            .masters
            .iter()
            .map(|master| {
                let relative = master["parent_graph"].as_str().unwrap();
                data.read(relative).unwrap();
                data.path(relative).unwrap().unwrap()
            })
            .collect::<BTreeSet<_>>();
        let expected_parents = [30527, 30699, 31740, 32745]
            .map(|sector| {
                data.path(&format!("graphs/parent-{sector}.dot"))
                    .unwrap()
                    .unwrap()
            })
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(parent_references, expected_parents);
        let family_graph = data
            .path(family["graph"].as_str().unwrap())
            .unwrap()
            .unwrap();
        assert_eq!(
            fs::read(&family_graph).unwrap(),
            data.read("masters/P640/graph.dot").unwrap()
        );
        let mut dotted_powers = powers.clone();
        dotted_powers[1] += 1;
        let unknown = catalog.symbol(family["fingerprint"].as_str().unwrap(), &dotted_powers);
        assert!(!catalog.masters.contains_key(&unknown));
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            let error = catalog
                .finalize(unknown.clone(), &settings, 5, mode)
                .unwrap_err()
                .to_string();
            assert!(error.contains("undeclared residual master"), "{error}");
            assert!(error.contains("no unique parent inferred"), "{error}");
            for parent in &parent_references {
                assert_eq!(
                    error.matches(parent.to_str().unwrap()).count(),
                    1,
                    "{error}"
                );
            }
            let positions = parent_references
                .iter()
                .map(|parent| error.find(parent.to_str().unwrap()).unwrap())
                .collect::<Vec<_>>();
            assert!(
                positions.windows(2).all(|pair| pair[0] < pair[1]),
                "{error}"
            );
            assert!(!error.contains(family_graph.to_str().unwrap()), "{error}");
            let cancel_unknown =
                vakint_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap() * &unknown;
            assert!(
                catalog
                    .finalize(cancel_unknown, &settings, 5, mode)
                    .unwrap()
                    .is_zero()
            );
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
