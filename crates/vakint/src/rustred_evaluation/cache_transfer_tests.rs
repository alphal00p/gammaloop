use super::*;
use crate::{RustRedCachePreparationOptions, RustRedCacheRoot, RustRedCacheWorkload};
use flate2::{Compression, write::GzEncoder};
use rustred::family::IntegralKey;
use rustred_app::CachePreparationStatus;
use std::{io::Cursor, sync::atomic::AtomicBool, time::Instant};

struct Fixture {
    scratch: tempfile::TempDir,
    engine: Vakint,
    source: PathBuf,
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
    let (_, program) = program(&fixture.engine, 1).unwrap();
    let bundle = SavedApplicationCacheBundle::open_generated(
        &fixture.output(),
        &[program.binding().clone()],
        Default::default(),
    )
    .unwrap();
    let source = SavedApplicationCacheBundle::open_generated(
        &fixture.source,
        &[program.binding().clone()],
        Default::default(),
    )
    .unwrap();
    let mut imported = bundle.open_reader(program.binding()).unwrap();
    let mut original = source.open_reader(program.binding()).unwrap();
    for power in [2, 3, 4] {
        let key = rustred::persistence::application_cache::ApplicationCacheKey::guarded_root(
            &IntegralKey::try_new([power]).unwrap(),
        );
        let expected = original.lookup(&key, program.context()).unwrap();
        assert!(expected.is_some(), "prepared root/descendant is absent");
        assert_eq!(imported.lookup(&key, program.context()).unwrap(), expected);
    }
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
