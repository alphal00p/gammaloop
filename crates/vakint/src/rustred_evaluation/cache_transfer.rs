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
    time::Duration,
};

const BUFFER_BYTES: usize = 64 * 1024;
// Conventional tar writers pad the final record to twenty 512-byte blocks.
const TAR_PADDING_BYTES: u64 = 10 * 1024;

/// Explicit logical-byte bounds on transport/staging, separate from reduction
/// RAM and filesystem allocated blocks or metadata.
/// The checksum must come from a trusted source independently of the archive.
#[derive(Clone, Debug)]
pub struct RustRedCacheTransferOptions {
    pub expected_sha256: [u8; 32],
    pub max_archive_bytes: u64,
    pub max_extracted_bytes: u64,
    /// Counts all tar entries, including directories; also bounds path metadata.
    pub max_entries: usize,
    /// Archive + extracted files + verified native snapshot, concurrently.
    pub max_temporary_bytes: u64,
    pub timeout: Duration,
    pub storage: SavedApplicationCachePackLimits,
}

impl RustRedCacheTransferOptions {
    pub fn new(expected_sha256: [u8; 32]) -> Self {
        Self {
            expected_sha256,
            max_archive_bytes: 16_000_000_000,
            max_extracted_bytes: 16_000_000_000,
            max_entries: 100_000,
            max_temporary_bytes: 48_000_000_000,
            timeout: Duration::from_secs(3600),
            storage: SavedApplicationCachePackLimits::default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RustRedCacheTransferReport {
    pub output: PathBuf,
    pub archive_bytes: u64,
    pub extracted_bytes: u64,
    pub entries: usize,
    pub sha256: [u8; 32],
    pub copied: SavedApplicationCacheCopyStatistics,
}

impl Vakint {
    /// Import a gzip-compressed tar with regular files/directories below
    /// `cache/`. Verify its required checksum before extracting, then publish
    /// an independent native snapshot at a NEW destination. No cache is loaded
    /// and preparation checkpoints are not copied. Supported on Linux/macOS.
    pub fn import_cache(
        &self,
        loop_count: usize,
        archive: impl AsRef<Path>,
        output: impl AsRef<Path>,
        options: &RustRedCacheTransferOptions,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        let transfer = Transfer::new(loop_count, output.as_ref(), options).map_err(operation)?;
        let input = File::open(archive)
            .map_err(input_error)
            .map_err(operation)?;
        transfer.receive(input)
    }

    /// Download only when explicitly called. HTTP(S) transport uses native TLS
    /// verification and a required independently trusted archive checksum.
    /// Redirects are disabled; supply the final URL of the `.tar.gz` bytes.
    /// Identity encoding is requested. The checksum and archive-byte bound
    /// apply after any HTTP content decoding performed by the transport.
    /// Loading the published cache is a separate call.
    pub fn download_cache(
        &self,
        loop_count: usize,
        url: &str,
        output: impl AsRef<Path>,
        options: &RustRedCacheTransferOptions,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        let transfer = Transfer::new(loop_count, output.as_ref(), options).map_err(operation)?;
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
    loop_count: usize,
    output: &'a Path,
    parent: &'a Path,
    options: &'a RustRedCacheTransferOptions,
}

impl<'a> Transfer<'a> {
    fn new(
        loop_count: usize,
        output: &'a Path,
        options: &'a RustRedCacheTransferOptions,
    ) -> Result<Self, AppError> {
        if !(1..=4).contains(&loop_count) {
            return Err(input_error(
                "cache transfer supports one through four loops",
            ));
        }
        if !cfg!(any(target_os = "linux", target_os = "macos")) {
            return Err(input_error(
                "atomic cache publication requires Linux or macOS",
            ));
        }
        if output.file_name().is_none()
            || options.timeout.is_zero()
            || options.max_archive_bytes == 0
            || options.max_extracted_bytes == 0
            || options.max_entries == 0
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
            loop_count,
            output,
            parent,
            options,
        })
    }

    fn receive(&self, mut input: impl Read) -> Result<RustRedCacheTransferReport, VakintError> {
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
        let bindings = if self.loop_count == 4 {
            super::super::four_loop::cache::programs()?
                .iter()
                .map(|program| program.binding().clone())
                .collect::<Vec<_>>()
        } else {
            vec![program(self.loop_count)?.1.binding().clone()]
        };
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
        limits.max_disk_bytes = limits
            .max_disk_bytes
            .min(self.options.max_temporary_bytes - archive_bytes - extracted_bytes);
        let verified = scratch.path().join("verified");
        let (writer, copied) = SavedApplicationCacheBundleWriter::copy_bases_to_new(
            &verified,
            &bindings,
            &bases.iter().collect::<Vec<_>>(),
            limits,
        )
        .map_err(operation)?;
        drop(writer);
        self.publish(&verified).map_err(operation)?;
        Ok(RustRedCacheTransferReport {
            output: self.output.to_path_buf(),
            archive_bytes,
            extracted_bytes,
            entries,
            sha256,
            copied,
        })
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
        let allowance = self
            .options
            .max_extracted_bytes
            .min(self.options.max_temporary_bytes - archive_bytes);
        let tar_allowance = (self.options.max_entries as u64)
            .checked_mul(1024)
            .and_then(|overhead| allowance.checked_add(overhead))
            .and_then(|value| value.checked_add(TAR_PADDING_BYTES))
            .ok_or_else(|| limit_error("tar allowance overflow"))?;
        let mut decoded = MultiGzDecoder::new(archive).take(
            tar_allowance
                .checked_add(1)
                .ok_or_else(|| limit_error("tar allowance overflow"))?,
        );
        let mut tar = tar::Archive::new(&mut decoded);
        let mut seen = BTreeSet::new();
        let mut bytes = 0u64;
        for entry in tar.entries().map_err(input_error)?.raw(true) {
            let mut entry = entry.map_err(input_error)?;
            if seen.len() >= self.options.max_entries {
                return Err(limit_error("archive entry count exceeds allowance"));
            }
            let kind = entry.header().entry_type();
            if !kind.is_file() && !kind.is_dir() {
                return Err(input_error(
                    "cache archive permits only regular files and directories",
                ));
            }
            let path = entry.path().map_err(input_error)?.into_owned();
            let mut components = path.components();
            if components.next() != Some(Component::Normal("cache".as_ref()))
                || components.any(|component| !matches!(component, Component::Normal(_)))
                || !seen.insert(path.clone())
            {
                return Err(input_error("invalid or duplicate cache archive path"));
            }
            let size = entry.size();
            if kind.is_dir() && size != 0 {
                return Err(input_error("archive directory contains data"));
            }
            bytes = bytes
                .checked_add(size)
                .ok_or_else(|| limit_error("extracted byte count overflow"))?;
            if bytes > allowance {
                return Err(limit_error("extracted cache exceeds transfer allowance"));
            }
            let destination = output.join(path);
            if kind.is_dir() {
                fs::create_dir_all(&destination).map_err(input_error)?;
            } else {
                fs::create_dir_all(
                    destination
                        .parent()
                        .ok_or_else(|| input_error("invalid file path"))?,
                )
                .map_err(input_error)?;
                let mut file = OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(destination)
                    .map_err(input_error)?;
                let copied = std::io::copy(&mut entry, &mut file).map_err(input_error)?;
                if copied != size {
                    return Err(input_error("truncated cache archive file"));
                }
            }
        }
        drop(tar);
        // Consume the gzip trailer to verify its CRC and reject data hidden
        // beyond the tar terminator. Padding is bounded just like tar headers.
        let mut buffer = [0u8; BUFFER_BYTES];
        loop {
            let read = decoded.read(&mut buffer).map_err(input_error)?;
            if read == 0 {
                break;
            }
            if buffer[..read].iter().any(|&byte| byte != 0) {
                return Err(input_error("data follows cache tar terminator"));
            }
        }
        if decoded.limit() == 0 {
            return Err(limit_error("decoded tar exceeds allowance"));
        }
        Ok((bytes, seen.len()))
    }

    fn publish(&self, verified: &Path) -> Result<(), AppError> {
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
mod tests;
