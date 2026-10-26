//! Multipart is transport metadata around unchanged native bundle files.
use super::*;
use rustred::persistence::application_cache::ApplicationCacheBinding;
use serde::{
    Deserialize, Deserializer,
    de::{Error, MapAccess, Visitor},
};
use std::collections::BTreeMap;

const SCHEMA: &str = "vakint.nativepack-transport.v1";
const MAX_PART_BYTES: u64 = 1_800_000_000;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Binding {
    artifact_digest: [u8; 32],
    family_fingerprint: String,
    physical_arity: usize,
    layer: String,
}

impl Binding {
    fn native(&self) -> ApplicationCacheBinding {
        ApplicationCacheBinding {
            artifact_digest: self.artifact_digest,
            family_fingerprint: self.family_fingerprint.clone(),
            physical_arity: self.physical_arity,
            layer: self.layer.clone(),
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct FileIdentity {
    bytes: u64,
    sha256: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Part {
    kind: String,
    archive: String,
    archive_bytes: u64,
    archive_sha256: String,
    members: Vec<String>,
    #[serde(default)]
    url: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(in crate::rustred_evaluation) struct Manifest {
    schema: String,
    loop_count: usize,
    binding: Binding,
    generation: u64,
    stored_record_references: u64,
    #[serde(deserialize_with = "unique_files")]
    files: BTreeMap<String, FileIdentity>,
    parts: Vec<Part>,
}

fn unique_files<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<BTreeMap<String, FileIdentity>, D::Error> {
    struct Files;
    impl<'de> Visitor<'de> for Files {
        type Value = BTreeMap<String, FileIdentity>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter.write_str("a native file inventory with unique paths")
        }

        fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
            let mut files = BTreeMap::new();
            while let Some(path) = map.next_key::<String>()? {
                match files.entry(path) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        entry.insert(map.next_value()?);
                    }
                    std::collections::btree_map::Entry::Occupied(_) => {
                        return Err(A::Error::custom("duplicate multipart file inventory path"));
                    }
                }
            }
            Ok(files)
        }
    }
    deserializer.deserialize_map(Files)
}

fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// Return true for native chunks, false for native indexes.
fn native_path(value: &str) -> Result<bool, AppError> {
    let pieces = value.split('/').collect::<Vec<_>>();
    match pieces.as_slice() {
        ["cache", "bundle.rrindex"] => Ok(false),
        ["cache", "bindings", binding, "cache.rrindex"] if hex(binding) => Ok(false),
        ["cache", "bindings", binding, "chunks", chunk]
            if hex(binding) && chunk.strip_suffix(".rrchunk.gz").is_some_and(hex) =>
        {
            Ok(true)
        }
        _ => Err(input_error(
            "multipart inventory contains a non-native file path",
        )),
    }
}

impl Manifest {
    pub(in crate::rustred_evaluation) const FORMAT: &str = SCHEMA;

    pub(in crate::rustred_evaluation) fn parse(
        bytes: &[u8],
        loop_count: usize,
    ) -> Result<Self, AppError> {
        let manifest: Self = serde_json::from_slice(bytes).map_err(input_error)?;
        manifest.validate(loop_count)?;
        Ok(manifest)
    }

    pub(in crate::rustred_evaluation) fn validate_program(
        &self,
        engine: &Vakint,
    ) -> Result<(), VakintError> {
        let bindings = engine.native_cache_bindings(self.loop_count)?;
        if bindings.len() != 1 {
            return Err(operation(input_error(
                "multipart descriptors require one native binding; this loop pack needs multiple bindings, so use its ordinary prepared-cache bundle",
            )).into());
        }
        if bindings != [self.binding.native()] {
            return Err(operation(input_error(
                "multipart namespace differs from the authenticated program",
            ))
            .into());
        }
        Ok(())
    }

    /// Official aggregate allowances follow the authenticated inventory rather
    /// than the unrelated size of one ordinary loop archive.
    pub(in crate::rustred_evaluation) fn official_options(
        &self,
        expected_sha256: [u8; 32],
        manifest_bytes: u64,
        timeout: Duration,
    ) -> Result<RustRedCacheTransferOptions, AppError> {
        let overflow = || limit_error("official cache inventory size overflow");
        let extracted = self.files.values().try_fold(0u64, |sum, file| {
            sum.checked_add(file.bytes).ok_or_else(overflow)
        })?;
        let archives = self.parts.iter().try_fold(manifest_bytes, |sum, part| {
            sum.checked_add(part.archive_bytes).ok_or_else(overflow)
        })?;
        // The verified native copy coexists with extracted files and one held
        // archive. This deliberately conservative grant includes every archive.
        let temporary = extracted
            .checked_mul(2)
            .and_then(|n| n.checked_add(archives))
            .ok_or_else(overflow)?;
        let mut options = RustRedCacheTransferOptions::multipart_manifest(expected_sha256);
        options.complete_operation_timeout = false;
        options.max_manifest_bytes = manifest_bytes;
        options.max_archive_bytes = archives;
        options.max_extracted_bytes = extracted;
        options.max_entries = self.files.len();
        options.max_temporary_bytes = temporary;
        options.timeout = timeout;
        options.storage.max_disk_bytes = extracted;
        Ok(options)
    }

    /// Reuse authenticates every installed native file without downloading,
    /// decompressing, importing coefficient state or changing the snapshot.
    pub(in crate::rustred_evaluation) fn validate_directory(
        &self,
        directory: &Path,
    ) -> Result<(), AppError> {
        let parent = directory
            .parent()
            .ok_or_else(|| input_error("official cache snapshot needs a parent directory"))?;
        let metadata = fs::symlink_metadata(parent).map_err(input_error)?;
        if !metadata.is_dir() || metadata.file_type().is_symlink() {
            return Err(input_error(
                "invalid official cache snapshot parent directory",
            ));
        }
        let mut buffer = [0u8; BUFFER_BYTES];
        for (relative, identity) in &self.files {
            checkpoint()?;
            let mut path = directory.to_owned();
            if !fs::symlink_metadata(&path).map_err(input_error)?.is_dir() {
                return Err(input_error(
                    "official cache snapshot is not a regular directory",
                ));
            }
            for component in Path::new(relative.strip_prefix("cache/").unwrap()).components() {
                path.push(component);
                if fs::symlink_metadata(&path)
                    .map_err(input_error)?
                    .file_type()
                    .is_symlink()
                {
                    return Err(input_error("symlink in official cache snapshot"));
                }
            }
            let metadata = fs::symlink_metadata(&path).map_err(input_error)?;
            if !metadata.is_file() || metadata.len() != identity.bytes {
                return Err(input_error("installed cache file length/type differs"));
            }
            let mut input = File::open(path).map_err(input_error)?;
            let mut hash = Sha256::new();
            let mut count = 0u64;
            loop {
                checkpoint()?;
                let read = input.read(&mut buffer).map_err(input_error)?;
                if read == 0 {
                    break;
                }
                count = count
                    .checked_add(read as u64)
                    .ok_or_else(|| limit_error("cache length overflow"))?;
                if count > identity.bytes {
                    return Err(input_error(
                        "installed cache grew beyond its authenticated length",
                    ));
                }
                hash.update(&buffer[..read]);
            }
            if count != identity.bytes || format!("{:x}", hash.finalize()) != identity.sha256 {
                return Err(input_error("installed cache file SHA256/length differs"));
            }
        }
        let binding = self.binding.native();
        let bundle = SavedApplicationCacheBundle::open_generated(
            directory,
            std::slice::from_ref(&binding),
            Default::default(),
        )?;
        let reader = bundle.open_reader(&binding)?;
        if reader.generation() != self.generation
            || reader.statistics().records as u64 != self.stored_record_references
        {
            return Err(input_error(
                "installed native generation/reference census differs from descriptor",
            ));
        }
        Ok(())
    }

    fn validate(&self, loop_count: usize) -> Result<(), AppError> {
        if self.schema != SCHEMA
            || self.loop_count != loop_count
            || self.files.is_empty()
            || !self.files.contains_key("cache/bundle.rrindex")
            || self.generation == 0
        {
            return Err(input_error(
                "incompatible multipart native snapshot manifest",
            ));
        }
        let mut members = BTreeSet::new();
        let mut archives = BTreeSet::new();
        let mut index_parts = 0;
        let mut namespaces = BTreeSet::new();
        for (path, identity) in &self.files {
            native_path(path)?;
            if let Some(namespace) = path.strip_prefix("cache/bindings/") {
                namespaces.insert(namespace.split('/').next().unwrap());
            }
            if !hex(&identity.sha256) {
                return Err(input_error("invalid multipart file checksum"));
            }
        }
        if namespaces.len() != 1
            || !self.files.contains_key(&format!(
                "cache/bindings/{}/cache.rrindex",
                namespaces
                    .iter()
                    .next()
                    .ok_or_else(|| input_error("missing native namespace"))?
            ))
        {
            return Err(input_error(
                "multipart descriptor requires one complete native namespace",
            ));
        }
        for part in &self.parts {
            let index = match part.kind.as_str() {
                "indexes" => {
                    index_parts += 1;
                    true
                }
                "chunks" => false,
                _ => return Err(input_error("unknown multipart container kind")),
            };
            if part.archive.is_empty()
                || part.archive.contains(['/', '\\'])
                || !part.archive.ends_with(".tar.gz")
                || !archives.insert(&part.archive)
                || part.archive_bytes == 0
                || part.archive_bytes > MAX_PART_BYTES
                || !hex(&part.archive_sha256)
                || part.members.is_empty()
            {
                return Err(input_error("invalid multipart archive identity"));
            }
            for path in &part.members {
                if !self.files.contains_key(path)
                    || !members.insert(path)
                    || native_path(path)? == index
                {
                    return Err(input_error(
                        "multipart members overlap or differ from their kind/inventory",
                    ));
                }
            }
        }
        if index_parts != 1 || members.len() != self.files.len() {
            return Err(input_error(
                "multipart requires one complete index container and an exact file partition",
            ));
        }
        Ok(())
    }
}

#[derive(Default)]
pub(super) struct Extraction {
    seen: BTreeSet<PathBuf>,
    pub(super) bytes: u64,
    pub(super) entries: usize,
}

impl Transfer<'_> {
    pub(super) fn admit_manifest_growth(
        &self,
        length: usize,
        capacity: usize,
        additional: usize,
    ) -> Result<(), AppError> {
        self.tick()?;
        let next = length
            .checked_add(additional)
            .ok_or_else(|| limit_error("multipart manifest length overflow"))?;
        let count = u64::try_from(next).map_err(input_error)?;
        if count > self.options.max_manifest_bytes
            || count > self.options.max_archive_bytes
            || count > self.options.max_temporary_bytes
        {
            return Err(limit_error("multipart manifest exceeds allowance"));
        }
        if next > capacity {
            symbolica::utils::resources::check_allocation(
                u64::try_from(next - capacity).map_err(input_error)?,
            )
            .map_err(AppError::resource)?;
        }
        Ok(())
    }

    pub(super) fn import_multipart(
        &self,
        path: &Path,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        self.multipart_budget(|| self.import_multipart_inner(path))
    }

    fn multipart_budget(
        &self,
        operation: impl FnOnce() -> Result<RustRedCacheTransferReport, VakintError>,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        self.tick().map_err(super::operation)?;
        if !self.options.complete_operation_timeout {
            return operation();
        }
        symbolica::utils::resources::ResourceBudget::new(
            symbolica::utils::resources::ResourceLimits {
                max_process_rss_bytes: None,
                max_elapsed: Some(self.options.timeout.saturating_sub(self.started.elapsed())),
            },
        )
        .run(operation)
    }

    fn import_multipart_inner(
        &self,
        path: &Path,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        if !fs::symlink_metadata(path)
            .map_err(input_error)
            .map_err(operation)?
            .is_file()
        {
            return Err(operation(input_error(
                "local multipart descriptor must be a regular file",
            ))
            .into());
        }
        let input = File::open(path).map_err(input_error).map_err(operation)?;
        self.receive_multipart(input, |part| {
            let path = path.parent().unwrap_or(Path::new(".")).join(&part.archive);
            if !fs::symlink_metadata(&path).map_err(input_error)?.is_file() {
                return Err(input_error(
                    "local multipart archive must be a regular file",
                ));
            }
            Ok(Box::new(File::open(path).map_err(input_error)?) as Box<dyn Read>)
        })
    }

    pub(super) fn download_multipart(
        &self,
        url: &str,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        self.multipart_budget(|| self.download_multipart_inner(url))
    }

    fn download_multipart_inner(
        &self,
        url: &str,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        let response = self.multipart_response(url).map_err(operation)?;
        self.receive_multipart(response.into_reader(), |part| {
            let target = if let Some(url) = &part.url {
                url.clone()
            } else {
                ureq::get(url)
                    .request_url()
                    .map_err(input_error)?
                    .as_url()
                    .join(&part.archive)
                    .map_err(input_error)?
                    .to_string()
            };
            Ok(self.multipart_response(&target)?.into_reader())
        })
    }

    fn multipart_response(&self, url: &str) -> Result<ureq::Response, AppError> {
        self.tick()?;
        let parsed = ureq::get(url).request_url().map_err(input_error)?;
        let loopback_test = cfg!(test)
            && parsed.scheme() == "http"
            && matches!(parsed.host(), "127.0.0.1" | "localhost" | "[::1]");
        if parsed.scheme() != "https" && !loopback_test
            || !parsed.as_url().username().is_empty()
            || parsed.as_url().password().is_some()
        {
            return Err(input_error(
                "multipart downloads require HTTPS without URL credentials",
            ));
        }
        let tls = ureq::native_tls::TlsConnector::new().map_err(input_error)?;
        let response = ureq::AgentBuilder::new()
            .tls_connector(Arc::new(tls))
            .https_only(!loopback_test)
            .redirects(8)
            .timeout(if self.options.complete_operation_timeout {
                self.options.timeout.saturating_sub(self.started.elapsed())
            } else {
                self.options.timeout
            })
            .build()
            .get(url)
            .set("Accept-Encoding", "identity")
            .call()
            .map_err(input_error)?;
        if response.status() != 200
            || response
                .header("Content-Encoding")
                .is_some_and(|v| v != "identity")
        {
            return Err(input_error("expected unencoded HTTP 200 multipart bytes"));
        }
        self.tick()?;
        Ok(response)
    }

    fn receive_multipart(
        &self,
        mut input: impl Read,
        mut part_reader: impl FnMut(&Part) -> Result<Box<dyn Read>, AppError>,
    ) -> Result<RustRedCacheTransferReport, VakintError> {
        self.tick().map_err(operation)?;
        let mut bytes = Vec::new();
        let mut block = [0u8; BUFFER_BYTES];
        loop {
            self.tick().map_err(operation)?;
            let read = input
                .read(&mut block)
                .map_err(input_error)
                .map_err(operation)?;
            if read == 0 {
                break;
            }
            self.admit_manifest_growth(bytes.len(), bytes.capacity(), read)
                .map_err(operation)?;
            bytes
                .try_reserve_exact(read)
                .map_err(input_error)
                .map_err(operation)?;
            bytes.extend_from_slice(&block[..read]);
        }
        let sha256: [u8; 32] = Sha256::digest(&bytes).into();
        if sha256 != self.options.expected_sha256 {
            return Err(operation(input_error("multipart manifest SHA-256 mismatch")).into());
        }
        let manifest = Manifest::parse(&bytes, self.loop_count).map_err(operation)?;
        // Bind before fetching any large part, under the SAME operation budget.
        // Saved-program authentication is never replaced by descriptor metadata.
        manifest.validate_program(self.engine)?;
        let bindings = [manifest.binding.native()];
        let mut archive_bytes = bytes.len() as u64;
        drop(bytes);
        let mut declared_archives = archive_bytes;
        let mut declared_files = 0u64;
        for file in manifest.files.values() {
            declared_files = declared_files
                .checked_add(file.bytes)
                .ok_or_else(|| operation(limit_error("multipart file census overflow")))?;
        }
        for part in &manifest.parts {
            declared_archives = declared_archives
                .checked_add(part.archive_bytes)
                .ok_or_else(|| operation(limit_error("multipart archive census overflow")))?;
        }
        if declared_archives > self.options.max_archive_bytes
            || declared_files > self.options.max_extracted_bytes
            || manifest.files.len() > self.options.max_entries
        {
            return Err(operation(limit_error("aggregate multipart allowance exceeded")).into());
        }
        let scratch = tempfile::Builder::new()
            .prefix(".vakint-cache-transfer-")
            .tempdir_in(self.parent)
            .map_err(input_error)
            .map_err(operation)?;
        let extracted = scratch.path().join("extracted");
        fs::create_dir(&extracted)
            .map_err(input_error)
            .map_err(operation)?;
        let mut state = Extraction::default();
        for part in &manifest.parts {
            self.tick().map_err(operation)?;
            let mut reader = part_reader(part).map_err(operation)?;
            let mut archive = tempfile::tempfile_in(scratch.path())
                .map_err(input_error)
                .map_err(operation)?;
            let mut hash = Sha256::new();
            let mut part_bytes = 0u64;
            loop {
                self.tick().map_err(operation)?;
                let read = reader
                    .read(&mut block)
                    .map_err(input_error)
                    .map_err(operation)?;
                if read == 0 {
                    break;
                }
                part_bytes = part_bytes
                    .checked_add(read as u64)
                    .ok_or_else(|| operation(limit_error("multipart part length overflow")))?;
                if part_bytes > part.archive_bytes
                    || state
                        .bytes
                        .checked_add(part_bytes)
                        .is_none_or(|n| n > self.options.max_temporary_bytes)
                {
                    return Err(operation(limit_error(
                        "multipart part exceeds declared/staging allowance",
                    ))
                    .into());
                }
                archive
                    .write_all(&block[..read])
                    .map_err(input_error)
                    .map_err(operation)?;
                hash.update(&block[..read]);
            }
            if part_bytes != part.archive_bytes
                || format!("{:x}", hash.finalize()) != part.archive_sha256
            {
                return Err(
                    operation(input_error("multipart archive length or SHA-256 mismatch")).into(),
                );
            }
            archive_bytes = archive_bytes
                .checked_add(part_bytes)
                .ok_or_else(|| operation(limit_error("multipart aggregate archive overflow")))?;
            archive
                .seek(SeekFrom::Start(0))
                .map_err(input_error)
                .map_err(operation)?;
            self.extract_into(
                archive,
                &extracted,
                part_bytes,
                &mut state,
                Some((&manifest, part)),
            )
            .map_err(operation)?;
        }
        let bundle = SavedApplicationCacheBundle::open_generated(
            &extracted.join("cache"),
            &bindings,
            self.options.storage,
        )
        .map_err(operation)?;
        let pack = bundle.open_reader(&bindings[0]).map_err(operation)?;
        if pack.generation() != manifest.generation
            || pack.statistics().records as u64 != manifest.stored_record_references
        {
            return Err(operation(input_error(
                "multipart native generation/reference census differs",
            ))
            .into());
        }
        let chunks = manifest
            .files
            .iter()
            .filter(|(path, _)| native_path(path).unwrap_or(false));
        let (count, compressed) = chunks.fold((0usize, 0u64), |(count, bytes), (_, file)| {
            (count + 1, bytes + file.bytes)
        });
        self.finish(
            &extracted,
            scratch.path(),
            TransferCensus {
                archive_bytes,
                held_archive_bytes: 0,
                extracted_bytes: state.bytes,
                entries: state.entries,
                sha256,
                expected_native_inventory: Some((count, compressed)),
            },
        )
    }

    pub(super) fn extract_into(
        &self,
        archive: File,
        output: &Path,
        archive_bytes: u64,
        state: &mut Extraction,
        expected: Option<(&Manifest, &Part)>,
    ) -> Result<(), AppError> {
        let allowance = self.options.max_extracted_bytes.min(
            self.options
                .max_temporary_bytes
                .checked_sub(archive_bytes)
                .ok_or_else(|| limit_error("archive exceeds staging allowance"))?,
        );
        let tar_allowance = if let Some((manifest, part)) = expected {
            part.members
                .iter()
                .try_fold(0u64, |n, path| n.checked_add(manifest.files[path].bytes))
                .and_then(|n| n.checked_add((part.members.len() as u64).checked_mul(1024)?))
                .and_then(|n| n.checked_add(TAR_PADDING_BYTES))
                .ok_or_else(|| limit_error("multipart tar allowance overflow"))?
        } else if self.options.diagnostics_enabled && allowance != u64::MAX {
            (self.options.max_entries as u64)
                .checked_mul(1024)
                .and_then(|n| n.checked_add(allowance))
                .and_then(|n| n.checked_add(TAR_PADDING_BYTES))
                .ok_or_else(|| limit_error("tar allowance overflow"))?
        } else {
            u64::MAX
        };
        let mut decoded = MultiGzDecoder::new(archive).take(tar_allowance.saturating_add(1));
        let mut local_members = BTreeSet::new();
        {
            let mut tar = tar::Archive::new(&mut decoded);
            for entry in tar.entries().map_err(input_error)?.raw(true) {
                let mut entry = entry.map_err(input_error)?;
                self.tick()?;
                if (self.options.diagnostics_enabled || expected.is_some())
                    && state.entries >= self.options.max_entries
                {
                    return Err(limit_error(
                        "aggregate archive entry count exceeds allowance",
                    ));
                }
                let kind = entry.header().entry_type();
                if !kind.is_file() && (!kind.is_dir() || expected.is_some()) {
                    return Err(input_error(
                        "cache archive permits only regular files and directories",
                    ));
                }
                let path = entry.path().map_err(input_error)?.into_owned();
                let mut components = path.components();
                if components.next() != Some(Component::Normal("cache".as_ref()))
                    || components.any(|component| !matches!(component, Component::Normal(_)))
                    || path.to_string_lossy().contains('\\')
                    || !state.seen.insert(path.clone())
                {
                    return Err(input_error("invalid or duplicate cache archive path"));
                }
                let identity = if let Some((manifest, part)) = expected {
                    let name = path
                        .to_str()
                        .ok_or_else(|| input_error("non-UTF8 multipart path"))?;
                    if !part.members.iter().any(|member| member == name)
                        || !local_members.insert(name.to_owned())
                    {
                        return Err(input_error(
                            "archive member not assigned to this multipart container",
                        ));
                    }
                    Some(&manifest.files[name])
                } else {
                    None
                };
                let size = entry.size();
                if kind.is_dir() && size != 0 || identity.is_some_and(|file| file.bytes != size) {
                    return Err(input_error("archive member size differs from inventory"));
                }
                state.bytes = state
                    .bytes
                    .checked_add(size)
                    .ok_or_else(|| limit_error("extracted byte overflow"))?;
                if state.bytes > allowance {
                    return Err(limit_error("extracted cache exceeds aggregate allowance"));
                }
                state.entries += 1;
                let destination = output.join(&path);
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
                    let mut hash = Sha256::new();
                    let mut copied = 0u64;
                    let mut buffer = [0u8; BUFFER_BYTES];
                    loop {
                        self.tick()?;
                        let read = entry.read(&mut buffer).map_err(input_error)?;
                        if read == 0 {
                            break;
                        }
                        file.write_all(&buffer[..read]).map_err(input_error)?;
                        if identity.is_some() {
                            hash.update(&buffer[..read]);
                        }
                        copied = copied
                            .checked_add(read as u64)
                            .ok_or_else(|| limit_error("member length overflow"))?;
                    }
                    if copied != size
                        || identity
                            .is_some_and(|file| file.sha256 != format!("{:x}", hash.finalize()))
                    {
                        return Err(input_error(
                            "cache archive member length or SHA-256 mismatch",
                        ));
                    }
                }
            }
        }
        let mut buffer = [0u8; BUFFER_BYTES];
        loop {
            self.tick()?;
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
        if let Some((_, part)) = expected
            && local_members.len() != part.members.len()
        {
            return Err(input_error("multipart container omits declared members"));
        }
        Ok(())
    }
}
