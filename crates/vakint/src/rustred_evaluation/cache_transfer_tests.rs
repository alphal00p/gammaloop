use super::*;
use crate::{RustRedCachePreparationOptions, RustRedCacheRoot, RustRedCacheWorkload};
use flate2::{Compression, write::GzEncoder};
use rustred::family::IntegralKey;
use rustred_app::CachePreparationStatus;
use serde_json::{Value, json};
use std::{collections::BTreeMap, io::Cursor, sync::atomic::AtomicBool, time::Instant};

pub(in crate::rustred_evaluation) struct Fixture {
    scratch: tempfile::TempDir,
    pub(in crate::rustred_evaluation) engine: Vakint,
    pub(in crate::rustred_evaluation) source: PathBuf,
    archive: PathBuf,
    options: RustRedCacheTransferOptions,
}

impl Fixture {
    fn new() -> Self {
        let scratch = tempfile::tempdir().unwrap();
        let engine = Vakint::new().unwrap();
        let source = scratch.path().join("source");
        let options = RustRedCachePreparationOptions {
            workload: RustRedCacheWorkload::Explicit(vec![RustRedCacheRoot {
                namespace: 0,
                integral: IntegralKey::try_new([4]).unwrap(),
            }]),
            ..Default::default()
        };
        let report = engine
            .prepare_cache(1, &source, &options, &AtomicBool::new(false), |_| Ok(()))
            .unwrap();
        assert_eq!(report.status, CachePreparationStatus::Complete);
        assert_eq!(report.completed_roots, 1);
        assert_eq!(report.newly_stored_roots, 1);
        assert_eq!(report.newly_stored_intermediates, 2);
        let archive = scratch.path().join("cache.tar.gz");
        let mut fixture = Self {
            scratch,
            engine,
            source,
            archive,
            options: RustRedCacheTransferOptions::new([0; 32]),
        };
        fixture.repack();
        fixture
    }

    fn repack(&mut self) {
        let encoder = GzEncoder::new(File::create(&self.archive).unwrap(), Compression::fast());
        let mut tar = tar::Builder::new(encoder);
        Self::append(&mut tar, &self.source, Path::new("cache"));
        tar.into_inner().unwrap().finish().unwrap();
        self.options.expected_sha256 = Sha256::digest(fs::read(&self.archive).unwrap()).into();
    }

    fn append<W: Write>(tar: &mut tar::Builder<W>, source: &Path, path: &Path) {
        let metadata = fs::metadata(source).unwrap();
        let mut header = tar::Header::new_ustar();
        header.set_mode(0o600);
        if metadata.is_dir() {
            header.set_entry_type(tar::EntryType::Directory);
            header.set_size(0);
            tar.append_data(&mut header, path, std::io::empty())
                .unwrap();
            let mut children = fs::read_dir(source)
                .unwrap()
                .map(|entry| entry.unwrap())
                .collect::<Vec<_>>();
            children.sort_by_key(|entry| entry.file_name());
            for child in children {
                Self::append(tar, &child.path(), &path.join(child.file_name()));
            }
        } else {
            header.set_entry_type(tar::EntryType::Regular);
            header.set_size(metadata.len());
            tar.append_data(&mut header, path, File::open(source).unwrap())
                .unwrap();
        }
    }

    fn output(&self) -> PathBuf {
        self.scratch.path().join("imported")
    }

    fn assert_clean_failure(&self, options: &RustRedCacheTransferOptions) {
        assert!(
            self.engine
                .import_cache(1, &self.archive, self.output(), options)
                .is_err()
        );
        assert!(!self.output().exists());
        self.assert_no_staging();
    }

    fn assert_no_staging(&self) {
        assert!(fs::read_dir(self.scratch.path()).unwrap().all(|entry| {
            !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with(".vakint-cache-transfer-")
        }));
    }

    pub(in crate::rustred_evaluation) fn assert_imported_records(&self, output: &Path) {
        let (_, program) = program(&self.engine, 1).unwrap();
        let bindings = [program.binding().clone()];
        let imported =
            SavedApplicationCacheBundle::open_generated(output, &bindings, Default::default())
                .unwrap();
        let original = SavedApplicationCacheBundle::open_generated(
            &self.source,
            &bindings,
            Default::default(),
        )
        .unwrap();
        let mut imported = imported.open_reader(program.binding()).unwrap();
        let mut original = original.open_reader(program.binding()).unwrap();
        for power in [2, 3, 4] {
            let key = rustred::persistence::application_cache::ApplicationCacheKey::guarded_root(
                &IntegralKey::try_new([power]).unwrap(),
            );
            let expected = original.lookup(&key, program.context()).unwrap();
            assert!(expected.is_some(), "prepared root/descendant is absent");
            assert_eq!(imported.lookup(&key, program.context()).unwrap(), expected);
        }
    }
}

#[test]
fn import_verifies_native_snapshot_and_does_not_attach_or_copy_preparation_checkpoint() {
    let fixture = Fixture::new();
    let report = fixture
        .engine
        .import_cache(1, &fixture.archive, fixture.output(), &fixture.options)
        .unwrap();
    assert_eq!(report.sha256, fixture.options.expected_sha256);
    assert_eq!(
        report.archive_bytes,
        fs::metadata(&fixture.archive).unwrap().len()
    );
    // I4 preparation retains its completed I2/I3 descendants as well; import
    // must preserve those shortcuts, not merely the explicitly selected root.
    assert_eq!(report.copied.preserved_record_references, 3);
    assert_eq!(report.copied.transfer_buffer_bytes, BUFFER_BYTES);
    assert!(!fixture.output().join("preparation.rrcheckpoint").exists());
    assert!(!fixture.engine.cache_statistics(1).unwrap().attached);
    fixture
        .engine
        .load_cache(1, fixture.output(), &Default::default())
        .unwrap();
    assert!(fixture.engine.cache_statistics(1).unwrap().attached);
    fixture.assert_imported_records(&fixture.output());
    fixture.assert_no_staging();
}

#[test]
fn checksum_binding_and_native_chunk_corruption_refuse_without_publication() {
    let mut fixture = Fixture::new();
    let mut wrong_hash = fixture.options.clone();
    wrong_hash.expected_sha256[0] ^= 1;
    fixture.assert_clean_failure(&wrong_hash);
    assert!(
        fixture
            .engine
            .import_cache(2, &fixture.archive, fixture.output(), &fixture.options)
            .is_err()
    );
    assert!(!fixture.output().exists());
    let bindings = fs::read_dir(fixture.source.join("bindings"))
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let chunk = fs::read_dir(bindings.join("chunks"))
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let mut data = fs::read(&chunk).unwrap();
    let last = data.len() - 1;
    data[last] ^= 1;
    fs::write(chunk, data).unwrap();
    fixture.repack();
    fixture.assert_clean_failure(&fixture.options);
}

#[test]
fn transport_extraction_entry_and_combined_staging_limits_are_independent() {
    let fixture = Fixture::new();
    for kind in 0..4 {
        let mut options = fixture.options.clone();
        match kind {
            0 => options.max_archive_bytes = fs::metadata(&fixture.archive).unwrap().len() - 1,
            1 => options.max_extracted_bytes = 0,
            2 => {
                options.diagnostics_enabled = true;
                options.max_entries = 0;
            }
            _ => options.max_temporary_bytes = fs::metadata(&fixture.archive).unwrap().len(),
        }
        fixture.assert_clean_failure(&options);
    }
    // Permit archive and extraction exactly, leaving no grant for the verified
    // native copy. Its own preflight must refuse before publication.
    let report = fixture
        .engine
        .import_cache(1, &fixture.archive, fixture.output(), &fixture.options)
        .unwrap();
    let output = fixture.scratch.path().join("no-copy-budget");
    let mut options = fixture.options.clone();
    options.max_temporary_bytes = report.archive_bytes + report.extracted_bytes;
    assert!(
        fixture
            .engine
            .import_cache(1, &fixture.archive, &output, &options)
            .is_err()
    );
    assert!(!output.exists());
    fixture.assert_no_staging();
}

#[test]
fn entry_count_template_is_inactive_until_diagnostics_are_enabled() {
    let fixture = Fixture::new();
    let mut options = fixture.options.clone();
    assert!(!options.diagnostics_enabled);
    options.max_entries = 0;
    let report = fixture
        .engine
        .import_cache(1, &fixture.archive, fixture.output(), &options)
        .unwrap();
    assert!(report.entries > 0);
    fixture.assert_no_staging();
}

#[test]
fn exact_entry_and_payload_limits_allow_standard_ustar_record_padding() {
    let fixture = Fixture::new();
    let mut tar = tar::Builder::new(Vec::new());
    Fixture::append(&mut tar, &fixture.source, Path::new("cache"));
    let mut bytes = tar.into_inner().unwrap();
    bytes.resize(
        bytes.len().div_ceil(TAR_PADDING_BYTES as usize) * TAR_PADDING_BYTES as usize,
        0,
    );
    let mut tar = tar::Archive::new(Cursor::new(&bytes));
    let sizes = tar
        .entries()
        .unwrap()
        .map(|entry| entry.unwrap().size())
        .collect::<Vec<_>>();
    let mut gzip = GzEncoder::new(Vec::new(), Compression::fast());
    gzip.write_all(&bytes).unwrap();
    let bytes = gzip.finish().unwrap();
    let mut options = fixture.options.clone();
    options.expected_sha256 = Sha256::digest(&bytes).into();
    options.diagnostics_enabled = true;
    options.max_entries = sizes.len();
    options.max_extracted_bytes = sizes.iter().sum();
    let output = fixture.output();
    let report = Transfer::new(&fixture.engine, 1, &output, &options)
        .unwrap()
        .receive(Cursor::new(bytes))
        .unwrap();
    assert_eq!(report.entries, options.max_entries);
    assert_eq!(report.extracted_bytes, options.max_extracted_bytes);
    fixture.assert_no_staging();
}

#[test]
fn existing_destination_and_racing_publication_are_never_overwritten() {
    let fixture = Fixture::new();
    fs::create_dir(fixture.output()).unwrap();
    let sentinel = fixture.output().join("keep");
    fs::write(&sentinel, b"user data").unwrap();
    assert!(
        fixture
            .engine
            .import_cache(1, &fixture.archive, fixture.output(), &fixture.options)
            .is_err()
    );
    assert_eq!(fs::read(&sentinel).unwrap(), b"user data");
    // Exercise the final atomic boundary separately: a destination can appear
    // after initial admission while the download and verification are running.
    let race_output = fixture.scratch.path().join("race");
    let transfer = Transfer::new(&fixture.engine, 1, &race_output, &fixture.options).unwrap();
    let staged = fixture.scratch.path().join("staged");
    fs::create_dir(&staged).unwrap();
    fs::create_dir(&race_output).unwrap();
    assert!(transfer.publish(&staged).is_err());
    assert!(staged.is_dir());
    fixture.assert_no_staging();
}

fn hostile_archive(entries: &[(&str, tar::EntryType)], trailing: usize) -> Vec<u8> {
    let mut tar = tar::Builder::new(Vec::new());
    for &(path, kind) in entries {
        let mut header = tar::Header::new_ustar();
        header.set_entry_type(kind);
        header.set_size(0);
        header.set_mode(0o600);
        let raw = header.as_mut_bytes();
        raw[..100].fill(0);
        raw[..path.len()].copy_from_slice(path.as_bytes());
        header.set_cksum();
        tar.append(&header, std::io::empty()).unwrap();
    }
    let mut tar = tar.into_inner().unwrap();
    tar.resize(tar.len() + trailing, 0);
    let mut gzip = GzEncoder::new(Vec::new(), Compression::fast());
    gzip.write_all(&tar).unwrap();
    gzip.finish().unwrap()
}

#[test]
fn traversal_links_sparse_extensions_duplicates_and_padding_bombs_are_rejected() {
    let fixture = Fixture::new();
    for entries in [
        vec![("cache/../escape", tar::EntryType::Regular)],
        vec![("/cache/escape", tar::EntryType::Regular)],
        vec![("wrong/file", tar::EntryType::Regular)],
        vec![("cache/link", tar::EntryType::Symlink)],
        vec![("cache/link", tar::EntryType::Link)],
        vec![("cache/sparse", tar::EntryType::GNUSparse)],
        vec![("cache/pax", tar::EntryType::XHeader)],
        vec![
            ("cache/file", tar::EntryType::Regular),
            ("cache/file", tar::EntryType::Regular),
        ],
    ] {
        let bytes = hostile_archive(&entries, 0);
        let mut options = fixture.options.clone();
        options.expected_sha256 = Sha256::digest(&bytes).into();
        let output = fixture.output();
        assert!(
            Transfer::new(&fixture.engine, 1, &output, &options)
                .unwrap()
                .receive(Cursor::new(bytes))
                .is_err()
        );
        assert!(!output.exists());
        fixture.assert_no_staging();
    }
    let bytes = hostile_archive(&[], 100_000);
    let mut options = fixture.options.clone();
    options.expected_sha256 = Sha256::digest(&bytes).into();
    options.diagnostics_enabled = true;
    options.max_entries = 1;
    options.max_extracted_bytes = 1;
    let output = fixture.output();
    assert!(
        Transfer::new(&fixture.engine, 1, &output, &options)
            .unwrap()
            .receive(Cursor::new(bytes))
            .is_err()
    );
    fixture.assert_no_staging();
}

#[test]
fn truncated_gzip_with_matching_transport_checksum_is_rejected() {
    let fixture = Fixture::new();
    let mut bytes = fs::read(&fixture.archive).unwrap();
    bytes.truncate(bytes.len() - 5);
    let mut options = fixture.options.clone();
    options.expected_sha256 = Sha256::digest(&bytes).into();
    let output = fixture.output();
    assert!(
        Transfer::new(&fixture.engine, 1, &output, &options)
            .unwrap()
            .receive(Cursor::new(bytes))
            .is_err()
    );
    fixture.assert_no_staging();
}

#[test]
fn explicit_loopback_http_download_matches_local_import_and_refuses_large_response() {
    let fixture = Fixture::new();
    for too_large in [false, true] {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let address = listener.local_addr().unwrap();
        let bytes = fs::read(&fixture.archive).unwrap();
        let server = std::thread::spawn(move || {
            let deadline = Instant::now() + Duration::from_secs(10);
            let (mut stream, _) = loop {
                match listener.accept() {
                    Ok(connection) => break connection,
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        assert!(
                            Instant::now() < deadline,
                            "loopback download did not connect within ten seconds"
                        );
                        std::thread::sleep(Duration::from_millis(10));
                    }
                    Err(error) => panic!("loopback server accept failed: {error}"),
                }
            };
            stream
                .set_read_timeout(Some(Duration::from_secs(10)))
                .unwrap();
            let mut request = [0u8; 4096];
            let mut received = 0;
            while !request[..received]
                .windows(4)
                .any(|part| part == b"\r\n\r\n")
            {
                assert!(
                    received < request.len(),
                    "loopback request headers too large"
                );
                let count = stream.read(&mut request[received..]).unwrap();
                assert_ne!(count, 0, "loopback request ended before its headers");
                received += count;
            }
            let _ = write!(
                stream,
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                bytes.len()
            );
            let _ = stream.write_all(&bytes);
        });
        let mut options = fixture.options.clone();
        options.timeout = Duration::from_secs(10);
        if too_large {
            options.max_archive_bytes = 1;
        }
        let output =
            fixture
                .scratch
                .path()
                .join(if too_large { "too-large" } else { "downloaded" });
        let result = fixture.engine.download_cache(
            1,
            &format!("http://{address}/cache.tar.gz"),
            &output,
            &options,
        );
        server.join().unwrap();
        if too_large {
            assert!(result.is_err());
            assert!(!output.exists());
        } else {
            assert_eq!(result.unwrap().copied.preserved_record_references, 3);
        }
        fixture.assert_no_staging();
    }
}

pub(in crate::rustred_evaluation) struct MultipartFixture {
    pub(in crate::rustred_evaluation) fixture: Fixture,
    pub(in crate::rustred_evaluation) descriptor: PathBuf,
    pub(in crate::rustred_evaluation) manifest: Value,
}

impl MultipartFixture {
    pub(in crate::rustred_evaluation) fn new() -> Self {
        let fixture = Fixture::new();
        let mut files = BTreeMap::new();
        Self::inventory(&fixture.source, Path::new("cache"), &mut files);
        assert!(files.len() >= 3);
        let (_, program) = program(&fixture.engine, 1).unwrap();
        let binding = program.binding();
        let bundle = SavedApplicationCacheBundle::open_generated(
            &fixture.source,
            std::slice::from_ref(binding),
            Default::default(),
        )
        .unwrap();
        let reader = bundle.open_reader(binding).unwrap();
        let indexes = files
            .keys()
            .filter(|path| path.ends_with(".rrindex"))
            .cloned()
            .collect::<Vec<_>>();
        let chunks = files
            .keys()
            .filter(|path| path.ends_with(".rrchunk.gz"))
            .cloned()
            .collect::<Vec<_>>();
        assert_eq!(indexes.len(), 2);
        assert!(!chunks.is_empty());
        let manifest = json!({
            "schema":"vakint.nativepack-transport.v1", "loop_count":1,
            "binding":{"artifact_digest":binding.artifact_digest,
                "family_fingerprint":binding.family_fingerprint,
                "physical_arity":binding.physical_arity,"layer":binding.layer},
            "generation":reader.generation(),
            "stored_record_references":reader.statistics().records,
            "files":files,
            "parts":[
                {"kind":"indexes","archive":"indexes.tar.gz","members":indexes},
                {"kind":"chunks","archive":"chunks.tar.gz","members":chunks}
            ]
        });
        let descriptor = fixture.scratch.path().join("transport.json");
        let mut multipart = Self {
            fixture,
            descriptor,
            manifest,
        };
        for part in 0..2 {
            multipart.replace_part(part, multipart.archive_part(part));
        }
        multipart.write_manifest();
        multipart
    }

    fn inventory(source: &Path, relative: &Path, files: &mut BTreeMap<String, Value>) {
        if source.is_dir() {
            for entry in fs::read_dir(source).unwrap() {
                let entry = entry.unwrap();
                Self::inventory(&entry.path(), &relative.join(entry.file_name()), files);
            }
        } else {
            let name = source.file_name().unwrap().to_string_lossy();
            if name.ends_with(".rrindex") || name.ends_with(".rrchunk.gz") {
                let bytes = fs::read(source).unwrap();
                files.insert(
                    relative.to_str().unwrap().to_owned(),
                    Self::identity(&bytes),
                );
            }
        }
    }

    fn identity(bytes: &[u8]) -> Value {
        json!({"bytes":bytes.len(),"sha256":format!("{:x}",Sha256::digest(bytes))})
    }

    pub(in crate::rustred_evaluation) fn part_path(&self, part: usize) -> PathBuf {
        self.fixture
            .scratch
            .path()
            .join(self.manifest["parts"][part]["archive"].as_str().unwrap())
    }

    fn archive_part(&self, part: usize) -> Vec<u8> {
        let mut archive = tar::Builder::new(Vec::new());
        for member in self.manifest["parts"][part]["members"].as_array().unwrap() {
            let member = member.as_str().unwrap();
            let source = self
                .fixture
                .source
                .join(member.strip_prefix("cache/").unwrap());
            let mut header = tar::Header::new_ustar();
            header.set_entry_type(tar::EntryType::Regular);
            header.set_mode(0o600);
            header.set_size(fs::metadata(&source).unwrap().len());
            archive
                .append_data(&mut header, member, File::open(source).unwrap())
                .unwrap();
        }
        Self::gzip(&archive.into_inner().unwrap())
    }

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut gzip = GzEncoder::new(Vec::new(), Compression::fast());
        gzip.write_all(bytes).unwrap();
        gzip.finish().unwrap()
    }

    fn replace_part(&mut self, part: usize, bytes: Vec<u8>) {
        fs::write(self.part_path(part), &bytes).unwrap();
        self.manifest["parts"][part]["archive_bytes"] = bytes.len().into();
        self.manifest["parts"][part]["archive_sha256"] =
            format!("{:x}", Sha256::digest(&bytes)).into();
    }

    fn write_manifest(&self) {
        fs::write(
            &self.descriptor,
            serde_json::to_vec(&self.manifest).unwrap(),
        )
        .unwrap();
    }

    fn options(&self) -> RustRedCacheTransferOptions {
        self.write_manifest();
        RustRedCacheTransferOptions::multipart_manifest(
            Sha256::digest(fs::read(&self.descriptor).unwrap()).into(),
        )
    }

    fn assert_failure(&self, options: &RustRedCacheTransferOptions) {
        assert!(
            self.fixture
                .engine
                .import_cache(1, &self.descriptor, self.fixture.output(), options)
                .is_err()
        );
        assert!(!self.fixture.output().exists());
        assert!(!self.fixture.engine.cache_statistics(1).unwrap().attached);
        self.fixture.assert_no_staging();
    }
}

#[test]
fn multipart_manifest_above_default_size_imports_with_explicit_larger_grant() {
    const MIB: u64 = 1024 * 1024;
    let multipart = MultipartFixture::new();
    let descriptor_bytes = 17 * MIB;
    let mut remaining = descriptor_bytes - fs::metadata(&multipart.descriptor).unwrap().len();
    let mut block = [b' '; BUFFER_BYTES];
    {
        let mut descriptor = OpenOptions::new()
            .append(true)
            .open(&multipart.descriptor)
            .unwrap();
        // JSON permits trailing whitespace. Stream padding and hashing so the
        // fixture does not retain another manifest-sized allocation.
        while remaining > 0 {
            let count = remaining.min(block.len() as u64) as usize;
            descriptor.write_all(&block[..count]).unwrap();
            remaining -= count as u64;
        }
    }
    assert_eq!(
        fs::metadata(&multipart.descriptor).unwrap().len(),
        descriptor_bytes
    );
    let mut descriptor = File::open(&multipart.descriptor).unwrap();
    let mut hash = Sha256::new();
    loop {
        let count = descriptor.read(&mut block).unwrap();
        if count == 0 {
            break;
        }
        hash.update(&block[..count]);
    }
    let mut options = RustRedCacheTransferOptions::multipart_manifest(hash.finalize().into());
    assert_eq!(options.max_manifest_bytes, 16 * MIB);
    multipart.assert_failure(&options);
    options.max_manifest_bytes = 32 * MIB;
    let output = multipart.fixture.output();
    let report = multipart
        .fixture
        .engine
        .import_cache(1, &multipart.descriptor, &output, &options)
        .unwrap();
    assert_eq!(report.sha256, options.expected_sha256);
    assert_eq!(report.copied.preserved_record_references, 3);
    assert!(
        !multipart
            .fixture
            .engine
            .cache_statistics(1)
            .unwrap()
            .attached
    );
    multipart.fixture.assert_imported_records(&output);
    multipart.fixture.assert_no_staging();
}

#[test]
fn multipart_manifest_growth_checks_ambient_allocation_headroom() {
    use symbolica::utils::resources::{ResourceBudget, ResourceLimitError, ResourceLimits};

    let scratch = tempfile::tempdir().unwrap();
    let engine = Vakint::new().unwrap();
    let output = scratch.path().join("imported");
    let options = RustRedCacheTransferOptions::multipart_manifest([0; 32]);
    let transfer = Transfer::new(&engine, 1, &output, &options).unwrap();
    let budget = ResourceBudget::with_memory_monitor(
        ResourceLimits {
            max_process_rss_bytes: Some(64),
            max_elapsed: None,
        },
        || Ok(64),
    );
    budget
        .run(|| transfer.admit_manifest_growth(4, 8, 4))
        .unwrap();
    let error = budget
        .run(|| transfer.admit_manifest_growth(4, 8, 5))
        .unwrap_err();
    assert_eq!(error.kind(), AppErrorKind::Limit);
    assert_eq!(
        error.resource_limit_error(),
        Some(&ResourceLimitError::Memory {
            rss_bytes: 64,
            additional_bytes: 1,
            limit_bytes: 64,
        })
    );
    ResourceBudget::with_memory_monitor(
        ResourceLimits {
            max_process_rss_bytes: Some(65),
            max_elapsed: None,
        },
        || Ok(64),
    )
    .run(|| transfer.admit_manifest_growth(4, 8, 5))
    .unwrap();
}

#[test]
fn multipart_manifest_exact_byte_grant_succeeds_and_one_byte_less_refuses() {
    let multipart = MultipartFixture::new();
    let mut options = multipart.options();
    let descriptor_bytes = fs::metadata(&multipart.descriptor).unwrap().len();
    options.max_manifest_bytes = descriptor_bytes - 1;
    multipart.assert_failure(&options);
    options.max_manifest_bytes = descriptor_bytes;
    let report = multipart
        .fixture
        .engine
        .import_cache(
            1,
            &multipart.descriptor,
            multipart.fixture.output(),
            &options,
        )
        .unwrap();
    assert_eq!(report.copied.preserved_record_references, 3);
    multipart.fixture.assert_no_staging();
}

#[test]
fn multipart_manifest_raw_duplicate_inventory_key_refuses_before_opening_parts() {
    let multipart = MultipartFixture::new();
    let mut descriptor = serde_json::to_string(&multipart.manifest).unwrap();
    let marker = "\"files\":{";
    let position = descriptor.find(marker).unwrap() + marker.len();
    let identity =
        serde_json::to_string(&multipart.manifest["files"]["cache/bundle.rrindex"]).unwrap();
    // Keep both identical keys in the raw JSON; Value would normalize them.
    descriptor.insert_str(position, &format!("\"cache/bundle.rrindex\":{identity},"));
    fs::write(&multipart.descriptor, descriptor.as_bytes()).unwrap();
    let options = RustRedCacheTransferOptions::multipart_manifest(
        Sha256::digest(descriptor.as_bytes()).into(),
    );
    for part in 0..multipart.manifest["parts"].as_array().unwrap().len() {
        fs::remove_file(multipart.part_path(part)).unwrap();
    }
    let error = multipart
        .fixture
        .engine
        .import_cache(
            1,
            &multipart.descriptor,
            multipart.fixture.output(),
            &options,
        )
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("duplicate multipart file inventory path"),
        "unexpected failure before the absent parts: {error}"
    );
    assert!(!multipart.fixture.output().exists());
    assert!(
        !multipart
            .fixture
            .engine
            .cache_statistics(1)
            .unwrap()
            .attached
    );
    multipart.fixture.assert_no_staging();
}

#[test]
fn multipart_offline_import_preserves_native_records_without_attaching() {
    let mut multipart = MultipartFixture::new();
    // Local descriptor import must use adjacent verified parts, without fetching
    // the optional remote locations supplied for explicit download callers.
    for part in multipart.manifest["parts"].as_array_mut().unwrap() {
        part["url"] = "http://127.0.0.1:9/not-used.tar.gz".into();
    }
    let options = multipart.options();
    let report = multipart
        .fixture
        .engine
        .import_cache(
            1,
            &multipart.descriptor,
            multipart.fixture.output(),
            &options,
        )
        .unwrap();
    let archive_bytes = fs::metadata(&multipart.descriptor).unwrap().len()
        + multipart.manifest["parts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|part| part["archive_bytes"].as_u64().unwrap())
            .sum::<u64>();
    let files = multipart.manifest["files"].as_object().unwrap();
    assert_eq!(report.archive_bytes, archive_bytes);
    assert_eq!(report.entries, files.len());
    assert_eq!(
        report.extracted_bytes,
        files
            .values()
            .map(|file| file["bytes"].as_u64().unwrap())
            .sum::<u64>()
    );
    assert_eq!(report.copied.preserved_record_references, 3);
    assert_eq!(report.sha256, options.expected_sha256);
    assert!(
        !multipart
            .fixture
            .engine
            .cache_statistics(1)
            .unwrap()
            .attached
    );
    assert!(
        !multipart
            .fixture
            .output()
            .join("preparation.rrcheckpoint")
            .exists()
    );
    multipart
        .fixture
        .assert_imported_records(&multipart.fixture.output());
    multipart.fixture.assert_no_staging();
}

#[test]
fn multipart_descriptor_part_and_native_integrity_are_independent() {
    let mut multipart = MultipartFixture::new();
    let mut options = multipart.options();
    options.expected_sha256[0] ^= 1;
    multipart.assert_failure(&options);

    let original = multipart.manifest.clone();
    for field in ["archive_sha256", "archive_bytes"] {
        multipart.manifest = original.clone();
        multipart.manifest["parts"][1][field] = if field == "archive_sha256" {
            "00".repeat(32).into()
        } else {
            (multipart.manifest["parts"][1][field].as_u64().unwrap() + 1).into()
        };
        multipart.assert_failure(&multipart.options());
    }
    multipart.manifest = original;
    let member = multipart.manifest["parts"][1]["members"][0]
        .as_str()
        .unwrap()
        .to_owned();
    let source = multipart
        .fixture
        .source
        .join(member.strip_prefix("cache/").unwrap());
    let mut bytes = fs::read(&source).unwrap();
    let last = bytes.len() - 1;
    bytes[last] ^= 1;
    fs::write(source, &bytes).unwrap();
    multipart.manifest["files"][&member] = MultipartFixture::identity(&bytes);
    multipart.replace_part(1, multipart.archive_part(1));
    // All transport checksums now agree with the changed file. Native's
    // committed chunk descriptor must still reject its changed payload.
    multipart.assert_failure(&multipart.options());
}

#[test]
fn multipart_binding_census_and_partition_contract_are_enforced() {
    let mut multipart = MultipartFixture::new();
    let original = multipart.manifest.clone();
    for case in 0..10 {
        multipart.manifest = original.clone();
        match case {
            0 => {
                multipart.manifest["binding"]["artifact_digest"][0] =
                    (original["binding"]["artifact_digest"][0].as_u64().unwrap() ^ 1).into()
            }
            1 => {
                multipart.manifest["generation"] =
                    (original["generation"].as_u64().unwrap() + 1).into()
            }
            2 => multipart.manifest["stored_record_references"] = 4.into(),
            3 => multipart.manifest["parts"][0]["kind"] = "chunks".into(),
            4 => multipart.manifest["parts"][1]["kind"] = "indexes".into(),
            5 => {
                let duplicate = multipart.manifest["parts"][0]["members"][0].clone();
                multipart.manifest["parts"][1]["members"]
                    .as_array_mut()
                    .unwrap()
                    .push(duplicate);
            }
            6 => {
                multipart.manifest["parts"][1]["members"]
                    .as_array_mut()
                    .unwrap()
                    .pop();
            }
            7 => {
                let member = multipart.manifest["parts"][1]["members"][0]
                    .as_str()
                    .unwrap()
                    .to_owned();
                let file = &mut multipart.manifest["files"][&member];
                file["bytes"] = (file["bytes"].as_u64().unwrap() + 1).into();
            }
            8 => multipart.manifest["parts"][1]["archive"] = "../outside.tar.gz".into(),
            _ => {
                let member = multipart.manifest["parts"][1]["members"][0]
                    .as_str()
                    .unwrap()
                    .to_owned();
                multipart.manifest["files"][&member]["sha256"] = "00".repeat(32).into();
            }
        }
        multipart.assert_failure(&multipart.options());
    }
}

#[test]
fn multipart_unreferenced_native_indexes_are_not_silently_discarded() {
    let mut multipart = MultipartFixture::new();
    let index = multipart.manifest["parts"][0]["members"]
        .as_array()
        .unwrap()
        .iter()
        .map(|member| member.as_str().unwrap())
        .find(|member| member.contains("/bindings/"))
        .unwrap()
        .to_owned();
    let namespace = index.split('/').nth(2).unwrap();
    let mut other_namespace = namespace.as_bytes().to_vec();
    other_namespace[0] = if other_namespace[0] == b'0' {
        b'1'
    } else {
        b'0'
    };
    let extra = format!(
        "cache/bindings/{}/cache.rrindex",
        String::from_utf8(other_namespace).unwrap()
    );
    let source = multipart
        .fixture
        .source
        .join(index.strip_prefix("cache/").unwrap());
    let destination = multipart
        .fixture
        .source
        .join(extra.strip_prefix("cache/").unwrap());
    fs::create_dir_all(destination.parent().unwrap()).unwrap();
    fs::write(destination, fs::read(source).unwrap()).unwrap();
    multipart.manifest["files"][&extra] = multipart.manifest["files"][&index].clone();
    multipart.manifest["parts"][0]["members"]
        .as_array_mut()
        .unwrap()
        .push(extra.into());
    multipart.replace_part(0, multipart.archive_part(0));
    multipart.assert_failure(&multipart.options());
}

#[test]
fn multipart_tar_paths_trailing_data_and_crc_fail_after_transport_authentication() {
    let mut multipart = MultipartFixture::new();
    let original = fs::read(multipart.part_path(1)).unwrap();
    for case in 0..7 {
        let bytes = match case {
            0 => hostile_archive(&[("cache/../escape", tar::EntryType::Regular)], 0),
            1 => hostile_archive(&[("cache/link", tar::EntryType::Symlink)], 0),
            2 => {
                let mut decoded = Vec::new();
                MultiGzDecoder::new(Cursor::new(&original))
                    .read_to_end(&mut decoded)
                    .unwrap();
                decoded.push(1);
                MultipartFixture::gzip(&decoded)
            }
            3 => {
                let mut bytes = original.clone();
                let crc = bytes.len() - 8;
                bytes[crc] ^= 1;
                bytes
            }
            4 => {
                let mut bytes = original.clone();
                bytes.extend(hostile_archive(
                    &[("cache/hidden", tar::EntryType::Regular)],
                    0,
                ));
                bytes
            }
            5 => {
                // The descriptor's partition stays disjoint, but the later
                // chunk archive actually tries to repeat an earlier index.
                let duplicate = multipart.manifest["parts"][0]["members"][0].clone();
                multipart.manifest["parts"][1]["members"]
                    .as_array_mut()
                    .unwrap()
                    .push(duplicate);
                let bytes = multipart.archive_part(1);
                multipart.manifest["parts"][1]["members"]
                    .as_array_mut()
                    .unwrap()
                    .pop();
                bytes
            }
            _ => original[..original.len() - 5].to_vec(),
        };
        multipart.replace_part(1, bytes);
        multipart.assert_failure(&multipart.options());
        assert!(!multipart.fixture.scratch.path().join("escape").exists());
    }
}

#[test]
fn multipart_limits_cover_all_parts_with_entry_limits_active_by_default() {
    let multipart = MultipartFixture::new();
    let options = multipart.options();
    assert!(!options.diagnostics_enabled);
    let report = multipart
        .fixture
        .engine
        .import_cache(
            1,
            &multipart.descriptor,
            multipart.fixture.scratch.path().join("baseline"),
            &options,
        )
        .unwrap();
    for case in 0..4 {
        let mut limited = options.clone();
        match case {
            0 => limited.max_archive_bytes = report.archive_bytes - 1,
            1 => limited.max_extracted_bytes = report.extracted_bytes - 1,
            2 => limited.max_entries = report.entries - 1,
            _ => limited.max_temporary_bytes = report.extracted_bytes,
        }
        multipart.assert_failure(&limited);
    }
    let mut exact = options;
    exact.max_archive_bytes = report.archive_bytes;
    exact.max_extracted_bytes = report.extracted_bytes;
    exact.max_entries = report.entries;
    let output = multipart.fixture.output();
    let exact_report = multipart
        .fixture
        .engine
        .import_cache(1, &multipart.descriptor, &output, &exact)
        .unwrap();
    assert_eq!(exact_report.entries, report.entries);
    assert_eq!(exact_report.extracted_bytes, report.extracted_bytes);
    multipart.fixture.assert_no_staging();
}

#[test]
fn multipart_expired_whole_transfer_refuses_before_publication() {
    let multipart = MultipartFixture::new();
    let mut options = multipart.options();
    options.timeout = Duration::from_secs(1);
    let output = multipart.fixture.output();
    let mut transfer = Transfer::new(&multipart.fixture.engine, 1, &output, &options).unwrap();
    transfer.started = Instant::now() - Duration::from_secs(2);
    assert!(transfer.import_multipart(&multipart.descriptor).is_err());
    assert!(!output.exists());
    multipart.fixture.assert_no_staging();
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
#[test]
fn multipart_descriptor_refuses_links_and_fifos_before_opening() {
    let multipart = MultipartFixture::new();
    let options = multipart.options();
    let link = multipart
        .fixture
        .scratch
        .path()
        .join("descriptor-link.json");
    std::os::unix::fs::symlink(&multipart.descriptor, &link).unwrap();
    assert!(
        multipart
            .fixture
            .engine
            .import_cache(1, &link, multipart.fixture.output(), &options)
            .is_err()
    );

    let fifo = multipart.fixture.scratch.path().join("descriptor.fifo");
    rustix::fs::mkfifoat(
        rustix::fs::CWD,
        &fifo,
        rustix::fs::Mode::RUSR | rustix::fs::Mode::WUSR,
    )
    .unwrap();
    let stop = Arc::new(AtomicBool::new(false));
    let writer_stop = Arc::clone(&stop);
    let writer_path = fifo.clone();
    let bytes = fs::read(&multipart.descriptor).unwrap();
    // If a regression opens the FIFO, supply a valid descriptor so the test
    // fails instead of leaving a blocked reader until the test-run deadline.
    let writer = std::thread::spawn(move || {
        let deadline = Instant::now() + Duration::from_secs(2);
        while !writer_stop.load(std::sync::atomic::Ordering::Relaxed) {
            if let Ok(fd) = rustix::fs::open(
                &writer_path,
                rustix::fs::OFlags::WRONLY | rustix::fs::OFlags::NONBLOCK,
                rustix::fs::Mode::empty(),
            ) {
                let _ = File::from(fd).write_all(&bytes);
                return true;
            }
            assert!(Instant::now() < deadline, "FIFO admission did not return");
            std::thread::sleep(Duration::from_millis(2));
        }
        false
    });
    let result =
        multipart
            .fixture
            .engine
            .import_cache(1, &fifo, multipart.fixture.output(), &options);
    stop.store(true, std::sync::atomic::Ordering::Relaxed);
    assert!(!writer.join().unwrap(), "descriptor FIFO was opened");
    assert!(result.is_err());
    assert!(!multipart.fixture.output().exists());
    multipart.fixture.assert_no_staging();
}

#[test]
fn multipart_http_download_verifies_parts_and_bounds_redirects() {
    let multipart = MultipartFixture::new();
    let mut options = multipart.options();
    options.timeout = Duration::from_secs(10);
    for case in 0..3 {
        let redirect = case == 1;
        let redirect_loop = case == 2;
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        listener.set_nonblocking(true).unwrap();
        let address = listener.local_addr().unwrap();
        let mut responses = Vec::new();
        if redirect {
            responses.push((
                "/transport.json".to_owned(),
                b"HTTP/1.1 302 Found\r\nLocation: /accepted.json\r\nContent-Length: 0\r\nConnection: close\r\n\r\n".to_vec(),
            ));
        }
        if redirect_loop {
            for _ in 0..10 {
                responses.push((
                    "/transport.json".to_owned(),
                    b"HTTP/1.1 302 Found\r\nLocation: /transport.json\r\nContent-Length: 0\r\nConnection: close\r\n\r\n".to_vec(),
                ));
            }
        }
        let descriptor_path = if redirect {
            "/accepted.json"
        } else {
            "/transport.json"
        };
        let descriptor = fs::read(&multipart.descriptor).unwrap();
        let parts = multipart.manifest["parts"].as_array().unwrap();
        for (path, body) in std::iter::once((descriptor_path.to_owned(), descriptor)).chain(
            parts.iter().enumerate().map(|(index, part)| {
                (
                    format!("/{}", part["archive"].as_str().unwrap()),
                    fs::read(multipart.part_path(index)).unwrap(),
                )
            }),
        ) {
            let mut response = format!(
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len(),
            )
            .into_bytes();
            response.extend(body);
            if !redirect_loop {
                responses.push((path, response));
            }
        }
        let stop = Arc::new(AtomicBool::new(false));
        let server_stop = Arc::clone(&stop);
        let server = std::thread::spawn(move || {
            let deadline = Instant::now() + Duration::from_secs(10);
            let mut requests = 0;
            for (path, response) in responses {
                let mut stream = loop {
                    if server_stop.load(std::sync::atomic::Ordering::Relaxed) {
                        return requests;
                    }
                    match listener.accept() {
                        Ok((stream, _)) => break stream,
                        Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                            assert!(
                                Instant::now() < deadline,
                                "multipart HTTP request timed out"
                            );
                            std::thread::sleep(Duration::from_millis(2));
                        }
                        Err(error) => panic!("multipart HTTP accept failed: {error}"),
                    }
                };
                stream
                    .set_read_timeout(Some(Duration::from_secs(10)))
                    .unwrap();
                let mut request = Vec::new();
                let mut block = [0; 1024];
                while !request.windows(4).any(|part| part == b"\r\n\r\n") {
                    let read = stream.read(&mut block).unwrap();
                    assert_ne!(read, 0);
                    request.extend_from_slice(&block[..read]);
                    assert!(request.len() <= 4096);
                }
                assert!(String::from_utf8_lossy(&request).starts_with(&format!("GET {path} ")));
                requests += 1;
                let _ = stream.write_all(&response);
            }
            requests
        });
        let output = multipart.fixture.scratch.path().join(match case {
            0 => "download",
            1 => "redirect",
            _ => "redirect-loop",
        });
        let result = multipart.fixture.engine.download_cache(
            1,
            &format!("http://{address}/transport.json"),
            &output,
            &options,
        );
        stop.store(true, std::sync::atomic::Ordering::Relaxed);
        let requests = server.join().unwrap();
        if redirect_loop {
            assert!(result.is_err());
            assert!(
                (2..=9).contains(&requests),
                "redirect limit did not bound requests: {requests}"
            );
            assert!(!output.exists());
        } else {
            assert_eq!(result.unwrap().copied.preserved_record_references, 3);
            assert_eq!(requests, 1 + parts.len() + usize::from(redirect));
            multipart.fixture.assert_imported_records(&output);
        }
        assert!(
            !multipart
                .fixture
                .engine
                .cache_statistics(1)
                .unwrap()
                .attached
        );
        multipart.fixture.assert_no_staging();
    }
}

#[test]
fn official_multipart_uses_http_timeout_without_a_default_compute_deadline() {
    let fixture = MultipartFixture::new();
    let manifest = NativePackManifest::parse(&fs::read(&fixture.descriptor).unwrap(), 1).unwrap();
    let options = manifest
        .official_options(
            [0; 32],
            fs::metadata(&fixture.descriptor).unwrap().len(),
            Duration::from_secs(1),
        )
        .unwrap();
    assert!(!options.complete_operation_timeout);
    let output = fixture.fixture.output();
    let mut transfer = Transfer::new(&fixture.fixture.engine, 1, &output, &options).unwrap();
    transfer.started = Instant::now().checked_sub(Duration::from_secs(2)).unwrap();
    transfer.tick().unwrap();
    let explicit = RustRedCacheTransferOptions::multipart_manifest([0; 32]);
    assert!(explicit.complete_operation_timeout);
    let output = fixture.fixture.output();
    let mut explicit_transfer =
        Transfer::new(&fixture.fixture.engine, 1, &output, &explicit).unwrap();
    explicit_transfer.started = Instant::now()
        .checked_sub(Duration::from_secs(3601))
        .unwrap();
    assert!(explicit_transfer.tick().is_err());
    let error = symbolica::utils::resources::ResourceBudget::new(
        symbolica::utils::resources::ResourceLimits {
            max_process_rss_bytes: None,
            max_elapsed: Some(Duration::ZERO),
        },
    )
    .run(|| transfer.tick())
    .unwrap_err();
    assert_eq!(error.kind(), AppErrorKind::Limit);
}
