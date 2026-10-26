//! Explicit transport of trusted, checksummed native cache snapshots.
//! RustRed validates program bindings and streams every referenced chunk.
use super::{operation, program};
use crate::{Vakint, VakintError};
use flate2::read::MultiGzDecoder;
use rustred_app::{
    AppError, AppErrorKind, SavedApplicationCacheBundle, SavedApplicationCacheBundleWriter,
    SavedApplicationCacheCopyStatistics, SavedApplicationCachePackLimits,
};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeSet,
    fs::{self, File, OpenOptions},
    io::{Read, Seek, SeekFrom, Write},
    path::{Component, Path, PathBuf},
    sync::Arc,
    time::{Duration, Instant},
};

#[path = "cache_transfer_multipart.rs"]
mod multipart;
pub(in crate::rustred_evaluation) use multipart::Manifest as NativePackManifest;

const BUFFER_BYTES: usize = 64 * 1024;
// Conventional tar writers pad the final record to twenty 512-byte blocks.
const TAR_PADDING_BYTES: u64 = 10 * 1024;

/// Explicit logical-byte bounds on transport/staging, separate from reduction
/// RAM and filesystem allocated blocks or metadata.
/// The checksum must come from a trusted source independently of the input.
#[derive(Clone, Debug)]
pub struct RustRedCacheTransferOptions {
    multipart_manifest: bool,
    /// Explicit transfers retain their complete-operation deadline. Official
    /// loading keeps its existing per-HTTP timeout and ambient resource budget.
    complete_operation_timeout: bool,
    pub expected_sha256: [u8; 32],
    /// Enable entry-count diagnostics for ordinary single-archive transfers.
    /// Multipart manifests always enforce their aggregate entry allowance.
    pub diagnostics_enabled: bool,
    /// Maximum authenticated multipart manifest bytes. Raise explicitly for
    /// large snapshots; ambient memory limits still govern allocation.
    pub max_manifest_bytes: u64,
    pub max_archive_bytes: u64,
    pub max_extracted_bytes: u64,
    /// Counts all tar entries, including directories; also bounds path metadata.
    pub max_entries: usize,
    /// Archive + extracted files + verified native snapshot, concurrently.
    pub max_temporary_bytes: u64,
    /// HTTP timeout; multipart mode also enforces one complete operation deadline.
    pub timeout: Duration,
    pub storage: SavedApplicationCachePackLimits,
}

impl RustRedCacheTransferOptions {
    pub fn new(expected_sha256: [u8; 32]) -> Self {
        Self {
            multipart_manifest: false,
            complete_operation_timeout: true,
            expected_sha256,
            diagnostics_enabled: false,
            max_manifest_bytes: 16 * 1024 * 1024,
            max_archive_bytes: u64::MAX,
            max_extracted_bytes: u64::MAX,
            max_entries: 100_000,
            max_temporary_bytes: u64::MAX,
            timeout: Duration::from_secs(3600),
            storage: SavedApplicationCachePackLimits::default(),
        }
    }
    pub fn with_diagnostics(mut self) -> Self {
        self.diagnostics_enabled = true;
        self
    }

    /// The required checksum authenticates a multipart transport manifest.
    /// Its parts contain unchanged native files, not a different cache codec.
    pub fn multipart_manifest(expected_sha256: [u8; 32]) -> Self {
        Self {
            multipart_manifest: true,
            ..Self::new(expected_sha256)
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RustRedCacheTransferReport {
    pub output: PathBuf,
    /// Archive bytes, or manifest plus all archive parts in multipart mode.
    pub archive_bytes: u64,
    pub extracted_bytes: u64,
    pub entries: usize,
    /// Authenticated input checksum: archive or multipart manifest, respectively.
    pub sha256: [u8; 32],
    pub copied: SavedApplicationCacheCopyStatistics,
}

impl Vakint {
    fn native_cache_bindings(
        &self,
        loop_count: usize,
    ) -> Result<Vec<rustred::persistence::application_cache::ApplicationCacheBinding>, VakintError>
    {
        if loop_count == 5 {
            let application = self.five_loop_cache_application()?;
            let prepared = application
                .prepared
                .lock()
                .map_err(|_| operation(input_error("saved application lock poisoned")))?;
            Ok(vec![prepared.cache_binding(application.cache_digest)])
        } else if loop_count == 4 {
            Ok(super::super::four_loop::cache::programs(self)?
                .iter()
                .map(|program| program.binding().clone())
                .collect())
        } else {
            Ok(vec![program(self, loop_count)?.1.binding().clone()])
        }
    }

    /// Import a gzip tar, or a multipart manifest selected by the options.
    /// Local multipart archives are regular sibling files named by the manifest.
    /// Verify the required input checksum and native files before publishing an
    /// independent snapshot at a NEW destination. This does not attach a cache;
    /// preparation checkpoints are not copied. Supported on Linux/macOS.
    pub fn import_cache(
        &self,
        loop_count: usize,
        archive: impl AsRef<Path>,
        output: impl AsRef<Path>,
        options: &RustRedCacheTransferOptions,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        let transfer =
            Transfer::new(self, loop_count, output.as_ref(), options).map_err(operation)?;
        if options.multipart_manifest {
            return transfer.import_multipart(archive.as_ref());
        }
        let input = File::open(archive)
            .map_err(input_error)
            .map_err(operation)?;
        transfer.receive(input)
    }

    /// Download only when explicitly called. HTTP(S) transport uses native TLS
    /// verification and a required independently trusted input checksum.
    /// Ordinary archive mode refuses redirects. Explicit multipart mode requires
    /// HTTPS and permits at most eight redirects, including release asset URLs.
    /// Identity encoding is requested; checksums cover the resulting input bytes.
    /// Loading the published cache is a separate call.
    pub fn download_cache(
        &self,
        loop_count: usize,
        url: &str,
        output: impl AsRef<Path>,
        options: &RustRedCacheTransferOptions,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        let transfer =
            Transfer::new(self, loop_count, output.as_ref(), options).map_err(operation)?;
        if options.multipart_manifest {
            return transfer.download_multipart(url);
        }
        if !url.starts_with("https://") && !url.starts_with("http://") {
            return Err(operation(input_error("cache URL must use HTTP or HTTPS")).into());
        }
        let tls = ureq::native_tls::TlsConnector::new()
            .map_err(input_error)
            .map_err(operation)?;
        let agent = ureq::AgentBuilder::new()
            .tls_connector(Arc::new(tls))
            .redirects(0)
            .timeout(options.timeout)
            .build();
        let response = agent
            .get(url)
            .set("Accept-Encoding", "identity")
            .call()
            .map_err(input_error)
            .map_err(operation)?;
        if response.status() != 200
            || response
                .header("Content-Encoding")
                .is_some_and(|value| value != "identity")
        {
            return Err(operation(input_error("expected unencoded HTTP 200 archive bytes")).into());
        }
        if let Some(length) = response.header("Content-Length") {
            let length = length
                .parse::<u64>()
                .map_err(input_error)
                .map_err(operation)?;
            if length > options.max_archive_bytes.min(options.max_temporary_bytes) {
                return Err(operation(limit_error("archive exceeds transfer allowance")).into());
            }
        }
        transfer.receive(response.into_reader())
    }
}

struct Transfer<'a> {
    engine: &'a Vakint,
    loop_count: usize,
    output: &'a Path,
    parent: &'a Path,
    options: &'a RustRedCacheTransferOptions,
    started: Instant,
}

struct TransferCensus {
    archive_bytes: u64,
    held_archive_bytes: u64,
    extracted_bytes: u64,
    entries: usize,
    sha256: [u8; 32],
    expected_native_inventory: Option<(usize, u64)>,
}

impl<'a> Transfer<'a> {
    fn new(
        engine: &'a Vakint,
        loop_count: usize,
        output: &'a Path,
        options: &'a RustRedCacheTransferOptions,
    ) -> Result<Self, AppError> {
        checkpoint()?;
        if !(1..=5).contains(&loop_count) {
            return Err(input_error(
                "cache transfer supports one through five loops",
            ));
        }
        if !cfg!(any(target_os = "linux", target_os = "macos")) {
            return Err(input_error(
                "atomic cache publication requires Linux or macOS",
            ));
        }
        if output.file_name().is_none()
            || options.timeout.is_zero()
            || (options.multipart_manifest && options.max_manifest_bytes == 0)
            || options.max_archive_bytes == 0
            || options.max_extracted_bytes == 0
            || ((options.diagnostics_enabled || options.multipart_manifest)
                && options.max_entries == 0)
            || options.max_temporary_bytes == 0
        {
            return Err(input_error(
                "a named destination and positive transfer limits are required",
            ));
        }
        match fs::symlink_metadata(output) {
            Ok(_) => return Err(input_error("cache transfer destination already exists")),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(input_error(error)),
        }
        let parent = output
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        Ok(Self {
            engine,
            loop_count,
            output,
            parent,
            options,
            started: Instant::now(),
        })
    }

    fn tick(&self) -> Result<(), AppError> {
        checkpoint()?;
        if self.options.multipart_manifest
            && self.options.complete_operation_timeout
            && self.started.elapsed() >= self.options.timeout
        {
            return Err(limit_error("complete multipart transfer deadline exceeded"));
        }
        Ok(())
    }

    fn bindings(
        &self,
    ) -> Result<Vec<rustred::persistence::application_cache::ApplicationCacheBinding>, VakintError>
    {
        self.engine.native_cache_bindings(self.loop_count)
    }

    fn receive(&self, mut input: impl Read) -> Result<RustRedCacheTransferReport, VakintError> {
        checkpoint().map_err(operation)?;
        let scratch = tempfile::Builder::new()
            .prefix(".vakint-cache-transfer-")
            .tempdir_in(self.parent)
            .map_err(input_error)
            .map_err(operation)?;
        let mut archive = tempfile::tempfile_in(scratch.path())
            .map_err(input_error)
            .map_err(operation)?;
        let (archive_bytes, sha256) = self.spool(&mut input, &mut archive).map_err(operation)?;
        if sha256 != self.options.expected_sha256 {
            return Err(operation(input_error("cache archive SHA-256 mismatch")).into());
        }
        archive
            .seek(SeekFrom::Start(0))
            .map_err(input_error)
            .map_err(operation)?;
        let extracted = scratch.path().join("extracted");
        fs::create_dir(&extracted)
            .map_err(input_error)
            .map_err(operation)?;
        let (extracted_bytes, entries) = self
            .extract(archive, &extracted, archive_bytes)
            .map_err(operation)?;
        self.finish(
            &extracted,
            scratch.path(),
            TransferCensus {
                archive_bytes,
                held_archive_bytes: archive_bytes,
                extracted_bytes,
                entries,
                sha256,
                expected_native_inventory: None,
            },
        )
    }

    fn finish(
        &self,
        extracted: &Path,
        scratch: &Path,
        census: TransferCensus,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        self.tick().map_err(operation)?;
        let bindings = self.bindings()?;
        let bundle = SavedApplicationCacheBundle::open_generated(
            &extracted.join("cache"),
            &bindings,
            self.options.storage,
        )
        .map_err(operation)?;
        let bases = bindings
            .iter()
            .map(|binding| bundle.open_reader(binding))
            .collect::<Result<Vec<_>, _>>()
            .map_err(operation)?;
        let mut limits = self.options.storage;
        limits.max_disk_bytes = limits.max_disk_bytes.min(
            self.options
                .max_temporary_bytes
                .checked_sub(census.held_archive_bytes)
                .and_then(|n| n.checked_sub(census.extracted_bytes))
                .ok_or_else(|| operation(limit_error("native copy staging allowance exhausted")))?,
        );
        let verified = scratch.join("verified");
        let copy = || {
            SavedApplicationCacheBundleWriter::copy_bases_to_new(
                &verified,
                &bindings,
                &bases.iter().collect::<Vec<_>>(),
                limits,
            )
        };
        let (writer, copied) =
            if self.options.multipart_manifest && self.options.complete_operation_timeout {
                self.tick().map_err(operation)?;
                let remaining = self.options.timeout.saturating_sub(self.started.elapsed());
                symbolica::utils::resources::ResourceBudget::new(
                    symbolica::utils::resources::ResourceLimits {
                        max_process_rss_bytes: None,
                        max_elapsed: Some(remaining),
                    },
                )
                .run(copy)
            } else {
                copy()
            }
            .map_err(operation)?;
        drop(writer);
        if census
            .expected_native_inventory
            .is_some_and(|(chunks, bytes)| {
                copied.copied_chunks != chunks || copied.compressed_bytes != bytes
            })
        {
            return Err(operation(input_error(
                "native referenced chunks differ from transport inventory",
            ))
            .into());
        }
        self.tick().map_err(operation)?;
        let report = RustRedCacheTransferReport {
            output: self.output.to_path_buf(),
            archive_bytes: census.archive_bytes,
            extracted_bytes: census.extracted_bytes,
            entries: census.entries,
            sha256: census.sha256,
            copied,
        };
        self.publish(&verified).map_err(operation)?;
        Ok(report)
    }

    fn spool(&self, input: &mut impl Read, output: &mut File) -> Result<(u64, [u8; 32]), AppError> {
        let allowance = self
            .options
            .max_archive_bytes
            .min(self.options.max_temporary_bytes);
        let mut bytes = 0u64;
        let mut hash = Sha256::new();
        let mut buffer = [0u8; BUFFER_BYTES];
        loop {
            checkpoint()?;
            let read = input.read(&mut buffer).map_err(input_error)?;
            if read == 0 {
                break;
            }
            bytes = bytes
                .checked_add(read as u64)
                .ok_or_else(|| limit_error("archive byte count overflow"))?;
            if bytes > allowance {
                return Err(limit_error("archive exceeds transfer allowance"));
            }
            output.write_all(&buffer[..read]).map_err(input_error)?;
            hash.update(&buffer[..read]);
        }
        output.flush().map_err(input_error)?;
        Ok((bytes, hash.finalize().into()))
    }

    fn extract(
        &self,
        archive: File,
        output: &Path,
        archive_bytes: u64,
    ) -> Result<(u64, usize), AppError> {
        let mut state = multipart::Extraction::default();
        self.extract_into(archive, output, archive_bytes, &mut state, None)?;
        Ok((state.bytes, state.entries))
    }

    fn publish(&self, verified: &Path) -> Result<(), AppError> {
        checkpoint()?;
        #[cfg(any(target_os = "linux", target_os = "macos"))]
        rustix::fs::renameat_with(
            rustix::fs::CWD,
            verified,
            rustix::fs::CWD,
            self.output,
            rustix::fs::RenameFlags::NOREPLACE,
        )
        .map_err(input_error)?;
        #[cfg(not(any(target_os = "linux", target_os = "macos")))]
        return Err(input_error(
            "atomic cache publication requires Linux or macOS",
        ));
        // A sync failure after rename is reported without deleting the published
        // snapshot. Such an error must not be interpreted as durable success.
        File::open(self.parent)
            .and_then(|parent| parent.sync_all())
            .map_err(input_error)
    }
}

fn input_error(error: impl std::fmt::Display) -> AppError {
    AppError::new(AppErrorKind::Input, format!("cache transfer: {error}"))
}
fn limit_error(message: &str) -> AppError {
    AppError::new(AppErrorKind::Limit, format!("cache transfer: {message}"))
}

#[cfg(test)]
#[path = "cache_transfer_tests.rs"]
pub(in crate::rustred_evaluation) mod tests;

fn checkpoint() -> Result<(), AppError> {
    symbolica::utils::resources::checkpoint().map_err(AppError::resource)
}
