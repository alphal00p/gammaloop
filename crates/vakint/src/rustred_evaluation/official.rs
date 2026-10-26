//! Explicit, checksummed installation of official per-loop release assets.

use std::{
    collections::BTreeSet,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    time::Duration,
};

use flate2::read::MultiGzDecoder;
use serde::Serialize;
use sha2::{Digest, Sha256};

use super::{
    RustRedEvaluationError,
    data::{self, LoopPack, ReleaseManifest, RustRedData},
};
use crate::{Vakint, VakintError};

const RELEASES: &str = "https://github.com/alphal00p/vakint-external-data/releases/download";
const MAX_MANIFEST_BYTES: u64 = 16 * 1024 * 1024;

#[derive(Clone, Debug)]
pub struct OfficialExternalDataOptions {
    pub version: String,
    /// None resolves the current working directory when the loader is called.
    pub local_artifacts_path: Option<PathBuf>,
    pub refresh: bool,
    /// Attach advertised native prepared caches. Disabling this still validates
    /// every package checksum, but never decodes or imports its cache backing.
    pub load_prepared_caches: bool,
    /// Timeout for each HTTP request. Native authentication, copying and reuse
    /// have no default compute deadline; caller resource budgets still apply.
    pub timeout: Duration,
    pub max_archive_bytes: u64,
    pub max_extracted_bytes: u64,
    pub max_entries: usize,
}

impl Default for OfficialExternalDataOptions {
    fn default() -> Self {
        Self {
            version: "latest".into(),
            local_artifacts_path: None,
            refresh: false,
            load_prepared_caches: true,
            timeout: Duration::from_secs(3600),
            max_archive_bytes: 2 * 1024 * 1024 * 1024,
            max_extracted_bytes: 8 * 1024 * 1024 * 1024,
            max_entries: 250_000,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ExternalDataLoopReport {
    pub loop_count: usize,
    pub path: PathBuf,
    pub cache_available: bool,
    /// Whether this activation imported the available prepared cache.
    pub cache_loaded: bool,
    pub cache_path: Option<PathBuf>,
    pub cache_format: Option<String>,
    /// Cache transport is separate from installing the ordinary loop pack.
    pub cache_source: Option<String>,
    /// "local", "reused", or "downloaded".
    pub source: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct ExternalDataLoadReport {
    pub version: String,
    pub path: PathBuf,
    pub loops: Vec<ExternalDataLoopReport>,
    pub load_prepared_caches: bool,
}

fn version_component(version: &str) -> Result<(), RustRedEvaluationError> {
    if version.is_empty()
        || version == "."
        || version == ".."
        || !version
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"._-".contains(&c))
    {
        return Err(data::failure(
            "version must be a release tag, not a filesystem path",
        ));
    }
    Ok(())
}

/// Distinguish a checkout from a directory that merely contains another
/// application's manifest. Recognizable incomplete checkouts must never turn
/// into download destinations when their release manifest is lost or damaged.
fn local_dataset(base: &Path) -> Result<bool, RustRedEvaluationError> {
    let manifest_path = base.join("manifest.json");
    let has_manifest = manifest_path.try_exists().map_err(data::failure)?;
    if has_manifest {
        let bytes = read_bounded(
            File::open(&manifest_path).map_err(data::failure)?,
            MAX_MANIFEST_BYTES,
        )?;
        let manifest: serde_json::Value = serde_json::from_slice(&bytes).map_err(|error| {
            data::failure(format!("invalid {}: {error}", manifest_path.display()))
        })?;
        let external_schema = manifest["schema"].as_str().is_some_and(|schema| {
            schema == "vakint.external-data" || schema.starts_with("vakint.external-data.")
        });
        if external_schema || manifest.get("loop_packs").is_some() {
            return Ok(true);
        }
    }
    for count in 1..=5 {
        if base
            .join(format!("{count}L/manifest.json"))
            .try_exists()
            .map_err(data::failure)?
        {
            return Err(data::failure(format!(
                "incomplete external-data checkout at {}: {count}L/manifest.json exists but {}",
                base.display(),
                if has_manifest {
                    "the root manifest does not declare external data"
                } else {
                    "the root manifest.json is missing"
                }
            )));
        }
    }
    Ok(false)
}

impl Vakint {
    /// Install and activate all official packs from 1L through `max_loop_count`.
    /// Existing repository roots are always read locally, including when refresh
    /// is requested. Ordinary engine construction never calls this method.
    pub fn load_official_external_data(
        &mut self,
        max_loop_count: usize,
        options: &OfficialExternalDataOptions,
    ) -> Result<ExternalDataLoadReport, VakintError> {
        self.load_official_external_data_from(max_loop_count, options, RELEASES)
    }

    fn load_official_external_data_from(
        &mut self,
        max_loop_count: usize,
        options: &OfficialExternalDataOptions,
        releases: &str,
    ) -> Result<ExternalDataLoadReport, VakintError> {
        version_component(&options.version)?;
        if !(1..=5).contains(&max_loop_count)
            || options.timeout.is_zero()
            || options.max_archive_bytes == 0
            || options.max_extracted_bytes == 0
            || options.max_entries == 0
        {
            return Err(data::failure(
                "loop count must be 1 through 5 and transfer limits must be positive",
            )
            .into());
        }
        let base = options
            .local_artifacts_path
            .clone()
            .map(Ok)
            .unwrap_or_else(std::env::current_dir)
            .map_err(data::failure)?;
        if local_dataset(&base)? {
            let loaded = RustRedData::from_directory(&base, &options.version, max_loop_count)?;
            let report = report(
                &loaded,
                base,
                max_loop_count,
                &BTreeSet::new(),
                &BTreeSet::new(),
                true,
                options.load_prepared_caches,
            )?;
            self.activate_external_data(loaded, options.load_prepared_caches)?;
            return Ok(report);
        }
        let cache = base.join(".vakint-external-data");
        fs::create_dir_all(&cache).map_err(data::failure)?;
        let lock = OpenOptions::new()
            .create(true)
            .truncate(false)
            .read(true)
            .write(true)
            .open(cache.join("install.lock"))
            .map_err(data::failure)?;
        lock.lock().map_err(data::failure)?;
        let latest = cache.join("latest");
        let cached_version = if options.version == "latest" && !options.refresh && latest.exists() {
            Some(
                fs::read_to_string(&latest)
                    .map_err(data::failure)?
                    .trim()
                    .to_owned(),
            )
        } else if options.version != "latest" {
            Some(options.version.clone())
        } else {
            None
        };
        if let Some(version) = &cached_version {
            version_component(version)?;
        }
        let local_manifest = cached_version
            .as_ref()
            .map(|v| cache.join(v).join("manifest.json"));
        let tls = ureq::native_tls::TlsConnector::new().map_err(data::failure)?;
        let agent = ureq::AgentBuilder::new()
            .tls_connector(std::sync::Arc::new(tls))
            .timeout(options.timeout)
            .redirects(8)
            .build();
        let bytes = if let Some(path) = local_manifest.filter(|p| p.is_file()) {
            read_bounded(File::open(path).map_err(data::failure)?, MAX_MANIFEST_BYTES)?
        } else {
            let tag = if options.version == "latest" {
                cached_version.as_deref().unwrap_or("data-latest")
            } else {
                &options.version
            };
            download_manifest(&agent, &format!("{releases}/{tag}/manifest.json"))?
        };
        let manifest: ReleaseManifest = serde_json::from_slice(&bytes).map_err(data::failure)?;
        validate_manifest(&manifest, max_loop_count)?;
        if options.version != "latest" && manifest.version != options.version {
            return Err(data::failure(format!(
                "release manifest declares {}, requested {}",
                manifest.version, options.version
            ))
            .into());
        }
        if cached_version
            .as_ref()
            .is_some_and(|v| v != &manifest.version)
        {
            return Err(data::failure(
                "installed manifest does not match its immutable release directory",
            )
            .into());
        }
        let destination = cache.join(&manifest.version);
        fs::create_dir_all(&destination).map_err(data::failure)?;
        let installed_manifest = destination.join("manifest.json");
        if installed_manifest.exists() {
            if fs::read(&installed_manifest).map_err(data::failure)? != bytes {
                return Err(data::failure(
                    "immutable release manifest changed; use a new release version",
                )
                .into());
            }
        } else {
            atomic_write(&installed_manifest, &bytes)?;
        }
        let mut downloaded = BTreeSet::new();
        for count in 1..=max_loop_count {
            let pack = manifest
                .loop_packs
                .iter()
                .find(|p| p.loop_count == count)
                .unwrap();
            if destination.join(format!("{count}L")).exists() {
                // A corrupt installation is an error, never an invitation to
                // overwrite user data or silently skip checksum validation.
                super::data::validate_loop_directory(&destination, pack, &manifest.version)?;
                continue;
            }
            install_pack(
                &agent,
                releases,
                &manifest.version,
                pack,
                &destination,
                &cache,
                options,
            )?;
            downloaded.insert(count);
        }
        let loaded = RustRedData::from_directory(&destination, &manifest.version, max_loop_count)?;
        let mut staged = self.clone();
        // Authenticate any selected NativePack against this new registry, never
        // the engine's preceding release. No backing is activated yet.
        staged.activate_external_data(loaded.clone(), false)?;
        let mut downloaded_caches = BTreeSet::new();
        if options.load_prepared_caches {
            for (&count, snapshot) in &loaded.runtime_cache_snapshots {
                snapshot.descriptor.validate_program(&staged)?;
                let path = loaded.runtime_cache_path(count)?;
                if path.try_exists().map_err(data::failure)? {
                    continue;
                }
                let parent = path.parent().unwrap();
                match fs::symlink_metadata(parent) {
                    Ok(metadata) if !metadata.is_dir() || metadata.file_type().is_symlink() => {
                        return Err(
                            data::failure("invalid official cache installation directory").into(),
                        );
                    }
                    Ok(_) => {}
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                        fs::create_dir(parent).map_err(data::failure)?;
                    }
                    Err(error) => return Err(data::failure(error).into()),
                }
                let transfer = snapshot
                    .descriptor
                    .official_options(
                        snapshot.reference.sha256()?,
                        snapshot.reference.manifest_bytes,
                        options.timeout,
                    )
                    .map_err(super::cache::operation)?;
                let basename = Path::new(&snapshot.reference.manifest)
                    .file_name()
                    .unwrap()
                    .to_str()
                    .unwrap();
                staged.download_cache(
                    count,
                    &format!("{releases}/{}/{basename}", manifest.version),
                    &path,
                    &transfer,
                )?;
                downloaded_caches.insert(count);
            }
            staged.attach_external_caches()?;
        }
        let report = report(
            &loaded,
            destination,
            max_loop_count,
            &downloaded,
            &downloaded_caches,
            false,
            options.load_prepared_caches,
        )?;
        if options.version == "latest" {
            atomic_write(&latest, manifest.version.as_bytes())?;
        }
        *self = staged;
        Ok(report)
    }
}

fn report(
    data: &RustRedData,
    path: PathBuf,
    maximum: usize,
    downloaded: &BTreeSet<usize>,
    downloaded_caches: &BTreeSet<usize>,
    local: bool,
    load_prepared_caches: bool,
) -> Result<ExternalDataLoadReport, RustRedEvaluationError> {
    Ok(ExternalDataLoadReport {
        version: data.version.clone(),
        load_prepared_caches,
        loops: (1..=maximum)
            .map(|n| {
                let pack = &data.loops[&n];
                let selected = data.runtime_cache_snapshots.contains_key(&n);
                let available = selected || pack.manifest.cache.available;
                let cache_path = if selected {
                    Some(data.runtime_cache_path(n)?)
                } else {
                    pack.manifest
                        .cache
                        .path
                        .as_deref()
                        .map(|p| pack.path(p))
                        .transpose()?
                        .flatten()
                };
                Ok(ExternalDataLoopReport {
                    loop_count: n,
                    path: path.join(format!("{n}L")),
                    cache_available: available,
                    cache_loaded: load_prepared_caches && available,
                    cache_path,
                    cache_format: if selected {
                        Some(super::cache::transfer::NativePackManifest::FORMAT.into())
                    } else {
                        pack.manifest.cache.format.clone()
                    },
                    cache_source: (load_prepared_caches && available).then(|| {
                        if local {
                            "local"
                        } else if downloaded_caches.contains(&n)
                            || (!selected && downloaded.contains(&n))
                        {
                            "downloaded"
                        } else {
                            "reused"
                        }
                        .into()
                    }),
                    source: if local {
                        "local"
                    } else if downloaded.contains(&n) {
                        "downloaded"
                    } else {
                        "reused"
                    }
                    .into(),
                })
            })
            .collect::<Result<Vec<_>, RustRedEvaluationError>>()?,
        path,
    })
}

fn validate_manifest(
    manifest: &ReleaseManifest,
    maximum: usize,
) -> Result<(), RustRedEvaluationError> {
    version_component(&manifest.version)?;
    if manifest.schema != data::RELEASE_SCHEMA
        || manifest.version == "latest"
        || manifest.version == "data-latest"
    {
        return Err(data::failure(
            "release manifest requires a supported schema and an immutable version",
        ));
    }
    for count in 1..=maximum {
        let packs = manifest
            .loop_packs
            .iter()
            .filter(|p| p.loop_count == count)
            .collect::<Vec<_>>();
        if packs.len() != 1 || packs[0].manifest != format!("{count}L/manifest.json") {
            return Err(data::failure(format!(
                "release has no unique {count}L pack"
            )));
        }
        let pack = packs[0];
        if data::relative_path(&pack.archive)?.components().count() != 1
            || pack.archive_sha256.as_ref().is_none_or(|h| h.len() != 64)
            || pack.archive_bytes.is_none()
        {
            return Err(data::failure(format!(
                "release {count}L archive lacks identity/size"
            )));
        }
    }
    Ok(())
}

fn download_manifest(agent: &ureq::Agent, url: &str) -> Result<Vec<u8>, RustRedEvaluationError> {
    let response = agent
        .get(url)
        .set("Accept-Encoding", "identity")
        .call()
        .map_err(data::failure)?;
    read_bounded(response.into_reader(), MAX_MANIFEST_BYTES)
}

fn read_bounded(reader: impl Read, limit: u64) -> Result<Vec<u8>, RustRedEvaluationError> {
    let mut bytes = Vec::new();
    reader
        .take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(data::failure)?;
    if bytes.len() as u64 > limit {
        return Err(data::failure("external-data manifest exceeds size limit"));
    }
    Ok(bytes)
}

fn atomic_write(path: &Path, bytes: &[u8]) -> Result<(), RustRedEvaluationError> {
    let mut file =
        tempfile::NamedTempFile::new_in(path.parent().unwrap()).map_err(data::failure)?;
    file.write_all(bytes).map_err(data::failure)?;
    file.as_file().sync_all().map_err(data::failure)?;
    file.persist(path).map_err(data::failure)?;
    Ok(())
}

fn install_pack(
    agent: &ureq::Agent,
    releases: &str,
    version: &str,
    pack: &LoopPack,
    destination: &Path,
    cache: &Path,
    options: &OfficialExternalDataOptions,
) -> Result<(), RustRedEvaluationError> {
    let expected_bytes = pack
        .archive_bytes
        .ok_or_else(|| data::failure("archive size is absent"))?;
    if expected_bytes > options.max_archive_bytes {
        return Err(data::failure(
            "release archive exceeds configured size limit",
        ));
    }
    let scratch = tempfile::Builder::new()
        .prefix(".install-")
        .tempdir_in(cache)
        .map_err(data::failure)?;
    let archive_path = scratch.path().join("pack.tar.gz");
    let mut archive = File::create(&archive_path).map_err(data::failure)?;
    let mut reader = agent
        .get(&format!("{releases}/{version}/{}", pack.archive))
        .set("Accept-Encoding", "identity")
        .call()
        .map_err(data::failure)?
        .into_reader();
    let mut digest = Sha256::new();
    let mut total = 0u64;
    let mut buffer = [0u8; 64 * 1024];
    loop {
        let n = reader.read(&mut buffer).map_err(data::failure)?;
        if n == 0 {
            break;
        }
        total = total
            .checked_add(n as u64)
            .ok_or_else(|| data::failure("archive length overflow"))?;
        if total > expected_bytes || total > options.max_archive_bytes {
            return Err(data::failure("download exceeds declared archive size"));
        }
        digest.update(&buffer[..n]);
        archive.write_all(&buffer[..n]).map_err(data::failure)?;
    }
    if total != expected_bytes || Some(format!("{:x}", digest.finalize())) != pack.archive_sha256 {
        return Err(data::failure("release archive size or SHA256 mismatch"));
    }
    drop(archive);
    let extracted = scratch.path().join("extracted");
    fs::create_dir(&extracted).map_err(data::failure)?;
    extract_pack(
        File::open(archive_path).map_err(data::failure)?,
        &extracted,
        pack.loop_count,
        options,
    )?;
    super::data::validate_loop_directory(&extracted, pack, version)?;
    fs::rename(
        extracted.join(format!("{}L", pack.loop_count)),
        destination.join(format!("{}L", pack.loop_count)),
    )
    .map_err(data::failure)?;
    Ok(())
}

fn extract_pack(
    reader: impl Read,
    destination: &Path,
    count: usize,
    options: &OfficialExternalDataOptions,
) -> Result<(), RustRedEvaluationError> {
    // Bound the complete inflated stream, including GNU/PAX extension records
    // consumed internally by tar before an ordinary entry is yielded.
    let inflated_limit = options
        .max_extracted_bytes
        .checked_add(1)
        .ok_or_else(|| data::failure("extracted size limit overflow"))?;
    let mut inflated = MultiGzDecoder::new(reader).take(inflated_limit);
    let mut archive = tar::Archive::new(&mut inflated);
    let mut seen = BTreeSet::new();
    let mut total = 0u64;
    for (index, entry) in archive.entries().map_err(data::failure)?.enumerate() {
        if index >= options.max_entries {
            return Err(data::failure("archive entry limit exceeded"));
        }
        let mut entry = entry.map_err(data::failure)?;
        let path = entry.path().map_err(data::failure)?.into_owned();
        let name = path
            .to_str()
            .ok_or_else(|| data::failure("non-UTF8 archive path"))?;
        data::relative_path(name)?;
        if path.components().next().unwrap().as_os_str()
            != std::ffi::OsStr::new(&format!("{count}L"))
            || !seen.insert(path.clone())
        {
            return Err(data::failure(
                "duplicate archive entry or wrong loop directory",
            ));
        }
        let kind = entry.header().entry_type();
        if !kind.is_file() && !kind.is_dir() {
            return Err(data::failure(
                "release archives may contain only regular files and directories",
            ));
        }
        total = total
            .checked_add(entry.size())
            .ok_or_else(|| data::failure("extracted size overflow"))?;
        if total > options.max_extracted_bytes {
            return Err(data::failure("extracted size limit exceeded"));
        }
        if !entry.unpack_in(destination).map_err(data::failure)? {
            return Err(data::failure("archive path escapes extraction directory"));
        }
    }
    if inflated.limit() == 0 {
        return Err(data::failure("inflated archive size limit exceeded"));
    }
    Ok(())
}

#[cfg(test)]
#[path = "official_tests.rs"]
mod tests;
