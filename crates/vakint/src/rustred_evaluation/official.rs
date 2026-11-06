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
/// Bounds of the `data-latest` check that an existing `latest` installation makes on
/// every load (one attempt, no retry, a small manifest behind a redirect): connecting,
/// and the whole request. The caller's `timeout` and `read_idle_timeout` lower both
/// (the read-idle bound applies to the whole check).
const LATEST_CHECK_CONNECT: Duration = Duration::from_secs(5);
const LATEST_CHECK_TOTAL: Duration = Duration::from_secs(20);

#[derive(Clone, Debug)]
pub struct OfficialExternalDataOptions {
    /// A release tag, or "latest" (the default): the release the `data-latest` alias
    /// names. An existing "latest" installation checks the alias on each load (one
    /// bounded attempt, no retry) and moves to the release it names; when the alias
    /// cannot be read, the cached release is used offline.
    pub version: String,
    /// None resolves the current working directory when the loader is called.
    pub local_artifacts_path: Option<PathBuf>,
    /// With "latest": resolve `data-latest` with the ordinary transfer options and
    /// report a failure instead of falling back to the cached release.
    pub refresh: bool,
    /// Optional timeout for each HTTP request; `None` (the default) sets none, so a
    /// large download on a slow link is never cut. Native authentication, copying
    /// and reuse have no compute deadline; caller resource budgets still apply.
    pub timeout: Option<Duration>,
    /// Read-idle detector, a progress check and not a transfer deadline: connecting, or a read
    /// of the response, that receives no byte for this long fails with a transport error, and
    /// the transfer is retried. `None` disables it. Default 300 s.
    pub read_idle_timeout: Option<Duration>,
    /// Retries of a transfer after a transport failure (refused or reset connection, read idle
    /// beyond `read_idle_timeout`, body ended before its declared size). A retry counter, not a
    /// limit on the data: a size or SHA256 mismatch and an HTTP error status are never retried.
    /// Default 5.
    pub transfer_retries: usize,
    /// Optional caps on the downloaded archive, the extracted bytes and the archive
    /// entries; `None` (the default) sets none. Every download is already bounded
    /// by the size and SHA256 that its release manifest declares, and extraction
    /// starts only after that check. `Some(0)` is rejected.
    pub max_archive_bytes: Option<u64>,
    pub max_extracted_bytes: Option<u64>,
    pub max_entries: Option<usize>,
    /// Workers verifying the five-loop routes while the application loads,
    /// usually the [`crate::VakintSettings::n_cores`] of later evaluations;
    /// one verifies on the calling thread. Defaults to [`crate::default_n_cores`].
    pub n_cores: usize,
}

impl Default for OfficialExternalDataOptions {
    fn default() -> Self {
        Self {
            version: "latest".into(),
            local_artifacts_path: None,
            refresh: false,
            timeout: None,
            read_idle_timeout: Some(Duration::from_secs(300)),
            transfer_retries: 5,
            max_archive_bytes: None,
            max_extracted_bytes: None,
            max_entries: None,
            n_cores: crate::default_n_cores(),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ExternalDataLoopReport {
    pub loop_count: usize,
    pub path: PathBuf,
    /// "local", "reused", or "downloaded".
    pub source: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct ExternalDataLoadReport {
    pub version: String,
    pub path: PathBuf,
    pub loops: Vec<ExternalDataLoopReport>,
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
            || options.timeout.is_some_and(|timeout| timeout.is_zero())
            || options.read_idle_timeout.is_some_and(|idle| idle.is_zero())
            || options.max_archive_bytes == Some(0)
            || options.max_extracted_bytes == Some(0)
            || options.max_entries == Some(0)
        {
            return Err(data::failure(
                "loop count must be 1 through 5 and transfer limits must be positive",
            )
            .into());
        }
        crate::execution::validate_n_cores(options.n_cores)?;
        let base = options
            .local_artifacts_path
            .clone()
            .map(Ok)
            .unwrap_or_else(std::env::current_dir)
            .map_err(data::failure)?;
        if local_dataset(&base)? {
            let loaded = RustRedData::from_directory(&base, &options.version, max_loop_count)?;
            let report = report(&loaded, base, max_loop_count, &BTreeSet::new(), true);
            self.activate_external_data(loaded, options.n_cores)?;
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
        let tls =
            std::sync::Arc::new(ureq::native_tls::TlsConnector::new().map_err(data::failure)?);
        let mut builder = ureq::AgentBuilder::new()
            .tls_connector(tls.clone())
            .redirects(8);
        if let Some(timeout) = options.timeout {
            builder = builder.timeout(timeout);
        }
        if let Some(idle) = options.read_idle_timeout {
            builder = builder.timeout_connect(idle).timeout_read(idle);
        }
        let agent = builder.build();
        // An existing `latest` installation follows the `data-latest` alias: when it
        // names another valid release, that release is installed and activated as with
        // `refresh` (a failure is reported and leaves the engine and the marker as they
        // were). When the alias cannot be read in time or names no valid release, the
        // cached release is used without network, as before.
        let followed = match &cached_version {
            Some(cached) if options.version == "latest" && !options.refresh => {
                moved_data_latest(&tls, releases, options, max_loop_count, cached)
            }
            _ => None,
        };
        let cached_version = cached_version.filter(|_| followed.is_none());
        let local_manifest = cached_version
            .as_ref()
            .map(|v| cache.join(v).join("manifest.json"));
        let bytes = if let Some(bytes) = followed {
            bytes
        } else if let Some(path) = local_manifest.filter(|p| p.is_file()) {
            read_bounded(File::open(path).map_err(data::failure)?, MAX_MANIFEST_BYTES)?
        } else {
            let tag = if options.version == "latest" {
                cached_version.as_deref().unwrap_or("data-latest")
            } else {
                &options.version
            };
            download_manifest(
                &agent,
                &format!("{releases}/{tag}/manifest.json"),
                options.transfer_retries,
            )?
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
        staged.activate_external_data(loaded.clone(), options.n_cores)?;
        let report = report(&loaded, destination, max_loop_count, &downloaded, false);
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
    local: bool,
) -> ExternalDataLoadReport {
    ExternalDataLoadReport {
        version: data.version.clone(),
        loops: (1..=maximum)
            .map(|n| ExternalDataLoopReport {
                loop_count: n,
                path: path.join(format!("{n}L")),
                source: if local {
                    "local"
                } else if downloaded.contains(&n) {
                    "downloaded"
                } else {
                    "reused"
                }
                .into(),
            })
            .collect(),
        path,
    }
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

/// Outcome of one transfer attempt: a transport failure is retried, anything else is final.
enum Transfer {
    Transport(RustRedEvaluationError),
    Final(RustRedEvaluationError),
}

fn call_error(error: ureq::Error) -> Transfer {
    match error {
        ureq::Error::Status(..) => Transfer::Final(data::failure(error)),
        ureq::Error::Transport(_) => Transfer::Transport(data::failure(error)),
    }
}

#[cfg(not(test))]
const RETRY_BACKOFF_MS: u64 = 1000;
#[cfg(test)]
const RETRY_BACKOFF_MS: u64 = 1;

/// Run `attempt` again after a transport failure, at most `retries` more times, waiting
/// 1, 2, 4, ... s (at most 60 s) between attempts.
fn retry_transfer<T>(
    retries: usize,
    mut attempt: impl FnMut() -> Result<T, Transfer>,
) -> Result<T, RustRedEvaluationError> {
    let mut done = 0usize;
    loop {
        match attempt() {
            Ok(value) => return Ok(value),
            Err(Transfer::Final(error)) => return Err(error),
            Err(Transfer::Transport(error)) if done >= retries => return Err(error),
            Err(Transfer::Transport(_)) => {
                done += 1;
                std::thread::sleep(Duration::from_millis(
                    (RETRY_BACKOFF_MS << (done - 1).min(6)).min(60_000),
                ));
            }
        }
    }
}

fn download_manifest(
    agent: &ureq::Agent,
    url: &str,
    retries: usize,
) -> Result<Vec<u8>, RustRedEvaluationError> {
    retry_transfer(retries, || {
        let response = agent
            .get(url)
            .set("Accept-Encoding", "identity")
            .call()
            .map_err(call_error)?;
        let mut bytes = Vec::new();
        response
            .into_reader()
            .take(MAX_MANIFEST_BYTES + 1)
            .read_to_end(&mut bytes)
            .map_err(|error| Transfer::Transport(data::failure(error)))?;
        if bytes.len() as u64 > MAX_MANIFEST_BYTES {
            return Err(Transfer::Final(data::failure(
                "external-data manifest exceeds size limit",
            )));
        }
        Ok(bytes)
    })
}

/// The `data-latest` manifest if it names a valid release other than `cached`, from one
/// bounded attempt without retry. Any failure (transport, HTTP status, size, an invalid
/// manifest or one without the requested loop packs) gives `None`: the cached release
/// stays in use.
fn moved_data_latest(
    tls: &std::sync::Arc<ureq::native_tls::TlsConnector>,
    releases: &str,
    options: &OfficialExternalDataOptions,
    max_loop_count: usize,
    cached: &str,
) -> Option<Vec<u8>> {
    // ureq's whole-request timeout takes precedence over a read timeout, so the
    // caller's read-idle bound lowers the whole-check bound instead.
    let shortest = |limit: Duration| {
        [options.timeout, options.read_idle_timeout]
            .into_iter()
            .flatten()
            .fold(limit, Duration::min)
    };
    let agent = ureq::AgentBuilder::new()
        .tls_connector(tls.clone())
        .redirects(8)
        .timeout_connect(shortest(LATEST_CHECK_CONNECT))
        .timeout(shortest(LATEST_CHECK_TOTAL))
        .build();
    let bytes =
        download_manifest(&agent, &format!("{releases}/data-latest/manifest.json"), 0).ok()?;
    let manifest: ReleaseManifest = serde_json::from_slice(&bytes).ok()?;
    validate_manifest(&manifest, max_loop_count).ok()?;
    (manifest.version != cached).then_some(bytes)
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
    if options
        .max_archive_bytes
        .is_some_and(|cap| expected_bytes > cap)
    {
        return Err(data::failure(
            "release archive exceeds configured size limit",
        ));
    }
    let scratch = tempfile::Builder::new()
        .prefix(".install-")
        .tempdir_in(cache)
        .map_err(data::failure)?;
    let archive_path = scratch.path().join("pack.tar.gz");
    // One attempt downloads the whole archive again; a transport failure (including a body
    // that ends before its declared size) is retried, a size overrun or SHA256 mismatch is final.
    retry_transfer(options.transfer_retries, || {
        let mut archive =
            File::create(&archive_path).map_err(|error| Transfer::Final(data::failure(error)))?;
        let mut reader = agent
            .get(&format!("{releases}/{version}/{}", pack.archive))
            .set("Accept-Encoding", "identity")
            .call()
            .map_err(call_error)?
            .into_reader();
        let mut digest = Sha256::new();
        let mut total = 0u64;
        let mut buffer = [0u8; 64 * 1024];
        loop {
            let n = reader
                .read(&mut buffer)
                .map_err(|error| Transfer::Transport(data::failure(error)))?;
            if n == 0 {
                break;
            }
            total = total
                .checked_add(n as u64)
                .ok_or_else(|| Transfer::Final(data::failure("archive length overflow")))?;
            if total > expected_bytes || options.max_archive_bytes.is_some_and(|cap| total > cap) {
                return Err(Transfer::Final(data::failure(
                    "download exceeds declared archive size",
                )));
            }
            digest.update(&buffer[..n]);
            archive
                .write_all(&buffer[..n])
                .map_err(|error| Transfer::Final(data::failure(error)))?;
        }
        if total < expected_bytes {
            return Err(Transfer::Transport(data::failure(
                "release archive ended before its declared size",
            )));
        }
        if Some(format!("{:x}", digest.finalize())) != pack.archive_sha256 {
            return Err(Transfer::Final(data::failure(
                "release archive size or SHA256 mismatch",
            )));
        }
        Ok(())
    })?;
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
    let inflated_limit = match options.max_extracted_bytes {
        Some(cap) => cap
            .checked_add(1)
            .ok_or_else(|| data::failure("extracted size limit overflow"))?,
        None => u64::MAX,
    };
    let mut inflated = MultiGzDecoder::new(reader).take(inflated_limit);
    let mut archive = tar::Archive::new(&mut inflated);
    let mut seen = BTreeSet::new();
    let mut total = 0u64;
    for (index, entry) in archive.entries().map_err(data::failure)?.enumerate() {
        if options.max_entries.is_some_and(|cap| index >= cap) {
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
        if options.max_extracted_bytes.is_some_and(|cap| total > cap) {
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
