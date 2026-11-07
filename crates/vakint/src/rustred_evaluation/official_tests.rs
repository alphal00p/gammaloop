use super::*;
use std::{
    collections::BTreeMap,
    net::TcpListener,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread,
};

fn copy_tree(source: &Path, destination: &Path) {
    fs::create_dir_all(destination).unwrap();
    for entry in fs::read_dir(source).unwrap() {
        let entry = entry.unwrap();
        let target = destination.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

fn fixture() -> tempfile::TempDir {
    let directory = tempfile::tempdir().unwrap();
    copy_tree(
        &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data"),
        directory.path(),
    );
    directory
}

fn fixture_manifest(root: &Path) -> ReleaseManifest {
    serde_json::from_slice(&fs::read(root.join("manifest.json")).unwrap()).unwrap()
}

#[test]
fn incomplete_checkouts_fail_without_network_or_local_mutation() {
    let directory = fixture();
    let root_manifest = directory.path().join("manifest.json");
    let loop_manifest = directory.path().join("1L/manifest.json");
    let original_loop_manifest = fs::read(&loop_manifest).unwrap();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    listener.set_nonblocking(true).unwrap();
    let releases = format!("http://{}", listener.local_addr().unwrap());
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(directory.path().to_owned()),
        timeout: Some(Duration::from_millis(100)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    for contents in [
        None,
        Some("{broken-json"),
        Some(r#"{"schema":"vakint.external-data.v99","version":"v0.1.0","loop_packs":[]}"#),
        Some(r#"{"loop_packs":[]}"#),
        Some(r#"{"schema":"unrelated.project.v1","name":"other application"}"#),
    ] {
        if let Some(contents) = contents {
            fs::write(&root_manifest, contents).unwrap();
        } else {
            fs::remove_file(&root_manifest).unwrap();
        }
        let error = engine
            .load_official_external_data_from(1, &options, &releases)
            .unwrap_err()
            .to_string();
        if contents.is_none() {
            assert!(error.contains("root manifest.json is missing"), "{error}");
        }
        assert!(
            matches!(listener.accept(), Err(error) if error.kind() == std::io::ErrorKind::WouldBlock)
        );
        assert!(!directory.path().join(".vakint-external-data").exists());
        assert_eq!(fs::read(&loop_manifest).unwrap(), original_loop_manifest);
        assert_eq!(fs::read_to_string(&root_manifest).ok().as_deref(), contents);
        assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    }
}

#[test]
fn checkout_is_cumulative_offline_and_version_checked_transactionally() {
    let directory = fixture();
    let mut engine = Vakint::new().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(directory.path().to_owned()),
        ..Default::default()
    };
    let report = engine
        .load_official_external_data_from(2, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(
        report
            .loops
            .iter()
            .map(|p| p.loop_count)
            .collect::<Vec<_>>(),
        [1, 2]
    );
    assert!(report.loops.iter().all(|p| p.source == "local"));
    assert!(!directory.path().join(".vakint-external-data").exists());
    let vendored = RustRedData::vendored().unwrap();
    let two_loop = engine.rustred_data().unwrap();
    assert_eq!(
        two_loop.loops.keys().copied().collect::<Vec<_>>(),
        [1, 2, 3]
    );
    assert!(Arc::ptr_eq(&two_loop.loops[&3], &vendored.loops[&3]));
    engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap();
    let one_loop = engine.rustred_data().unwrap();
    assert_eq!(
        one_loop.loops.keys().copied().collect::<Vec<_>>(),
        [1, 2, 3]
    );
    assert!(!Arc::ptr_eq(&one_loop.loops[&1], &vendored.loops[&1]));
    for count in [2, 3] {
        assert!(Arc::ptr_eq(
            &one_loop.loops[&count],
            &vendored.loops[&count]
        ));
    }
    let active = engine.rustred_data().unwrap();
    assert!(
        engine
            .load_official_external_data(
                2,
                &OfficialExternalDataOptions {
                    version: "v999.0.0".into(),
                    ..options.clone()
                }
            )
            .is_err()
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    fs::write(directory.path().join("1L/masters/P1/symbolic.json"), b"[]").unwrap();
    assert!(engine.load_official_external_data(2, &options).is_err());
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
}

fn release_assets(root: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut manifest = fixture_manifest(root);
    let mut assets = BTreeMap::new();
    for pack in &mut manifest.loop_packs {
        let gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        let mut tar = tar::Builder::new(gz);
        let label = format!("{}L", pack.loop_count);
        tar.append_dir_all(&label, root.join(&label)).unwrap();
        let bytes = tar.into_inner().unwrap().finish().unwrap();
        pack.archive_sha256 = Some(format!("{:x}", Sha256::digest(&bytes)));
        pack.archive_bytes = Some(bytes.len() as u64);
        assets.insert(format!("/{}/{}", manifest.version, pack.archive), bytes);
    }
    assets.insert(
        format!("/{}/manifest.json", manifest.version),
        serde_json::to_vec(&manifest).unwrap(),
    );
    assets
}

struct Server {
    url: String,
    stop: Arc<AtomicBool>,
    thread: Option<thread::JoinHandle<()>>,
}
impl Server {
    fn start(assets: BTreeMap<String, Vec<u8>>) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let url = format!("http://{}", listener.local_addr().unwrap());
        listener.set_nonblocking(true).unwrap();
        let stop = Arc::new(AtomicBool::new(false));
        let stopped = Arc::clone(&stop);
        let latest_manifest = assets
            .keys()
            .rfind(|path| path.ends_with("/manifest.json"))
            .unwrap()
            .clone();
        let thread = thread::spawn(move || {
            while !stopped.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((mut stream, _)) => {
                        stream
                            .set_read_timeout(Some(Duration::from_secs(2)))
                            .unwrap();
                        let mut request = [0u8; 4096];
                        let n = stream.read(&mut request).unwrap();
                        let text = String::from_utf8_lossy(&request[..n]);
                        let path = text.split_whitespace().nth(1).unwrap();
                        if path == "/data-latest/manifest.json" {
                            write!(stream, "HTTP/1.1 302 Found\r\nLocation: {latest_manifest}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n").unwrap();
                        } else if let Some(bytes) = assets.get(path) {
                            write!(stream, "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n", bytes.len()).unwrap();
                            stream.write_all(bytes).unwrap();
                        } else {
                            stream.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n").unwrap();
                        }
                    }
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        thread::sleep(Duration::from_millis(5))
                    }
                    Err(error) => panic!("fixture server: {error}"),
                }
            }
        });
        Self {
            url,
            stop,
            thread: Some(thread),
        }
    }
}
impl Drop for Server {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        self.thread.take().unwrap().join().unwrap();
    }
}

#[test]
fn unrelated_manifest_and_empty_loop_directory_allow_managed_installation() {
    let source = fixture();
    let version = fixture_manifest(source.path()).version;
    let server = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let unrelated =
        br#"{"schema":"unrelated.project.v1","version":"1","name":"other application"}"#;
    fs::write(destination.path().join("manifest.json"), unrelated).unwrap();
    fs::create_dir(destination.path().join("1L")).unwrap();
    let mut engine = Vakint::new().unwrap();
    let report = engine
        .load_official_external_data_from(
            1,
            &OfficialExternalDataOptions {
                local_artifacts_path: Some(destination.path().to_owned()),
                timeout: Some(Duration::from_secs(3)),
                ..Default::default()
            },
            &server.url,
        )
        .unwrap();
    assert_eq!(report.loops[0].source, "downloaded");
    assert_eq!(report.version, version);
    assert_eq!(
        report.path,
        destination
            .path()
            .join(".vakint-external-data")
            .join(&version)
    );
    assert_eq!(
        fs::read(destination.path().join("manifest.json")).unwrap(),
        unrelated
    );
    assert!(
        fs::read_dir(destination.path().join("1L"))
            .unwrap()
            .next()
            .is_none()
    );
}

#[test]
fn official_download_follows_redirect_then_reuses_offline_without_extraction() {
    let source = fixture();
    let server = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Some(Duration::from_secs(3)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let first = engine
        .load_official_external_data_from(2, &options, &server.url)
        .unwrap();
    assert!(first.loops.iter().all(|p| p.source == "downloaded"));
    let program = first.path.join("1L/ibps/rules.rrbin");
    let modified = fs::metadata(&program).unwrap().modified().unwrap();
    drop(server);
    let second = engine
        .load_official_external_data_from(2, &options, "http://127.0.0.1:1")
        .unwrap();
    assert!(second.loops.iter().all(|p| p.source == "reused"));
    assert_eq!(
        fs::metadata(&program).unwrap().modified().unwrap(),
        modified
    );
    let active = engine.rustred_data().unwrap();
    assert!(
        engine
            .load_official_external_data_from(
                2,
                &OfficialExternalDataOptions {
                    refresh: true,
                    ..options
                },
                "http://127.0.0.1:1"
            )
            .is_err()
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
}

#[test]
fn fixture_release_bypasses_stale_latest_and_reuses_immutable_data_offline() {
    let source = fixture();
    let required = super::super::test_data::external_data_fixture::required_version();
    assert_eq!(fixture_manifest(source.path()).version, required);
    let server = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let cache = destination.path().join(".vakint-external-data");
    fs::create_dir_all(&cache).unwrap();
    let stale = b"v0.1.0";
    assert_ne!(required.as_bytes(), stale);
    fs::write(cache.join("latest"), stale).unwrap();
    let options = OfficialExternalDataOptions {
        version: required.clone(),
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Some(Duration::from_secs(3)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let first = engine
        .load_official_external_data_from(1, &options, &server.url)
        .unwrap();
    assert_eq!(first.version, required);
    assert_eq!(first.path, cache.join(&required));
    assert_eq!(first.loops[0].source, "downloaded");
    assert_eq!(fs::read(cache.join("latest")).unwrap(), stale);
    let program = first.path.join("1L/ibps/rules.rrbin");
    let modified = fs::metadata(&program).unwrap().modified().unwrap();
    drop(server);

    let reused = engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(reused.path, first.path);
    assert_eq!(reused.version, required);
    assert_eq!(reused.loops[0].source, "reused");
    assert_eq!(
        fs::metadata(&program).unwrap().modified().unwrap(),
        modified
    );
    assert_eq!(fs::read(cache.join("latest")).unwrap(), stale);
}

#[test]
fn extraction_refuses_links_wrong_loop_and_resource_overflow() {
    for (path, kind, max_bytes) in [
        ("1L/link", tar::EntryType::Symlink, 100),
        ("2L/value", tar::EntryType::Regular, 100),
        ("1L/value", tar::EntryType::Regular, 2),
    ] {
        let mut header = tar::Header::new_gnu();
        header.set_path(path).unwrap();
        header.set_mode(0o644);
        header.set_entry_type(kind);
        let contents: &[u8] = if kind.is_symlink() {
            header.set_link_name("/outside").unwrap();
            b""
        } else {
            b"123"
        };
        header.set_size(contents.len() as u64);
        header.set_cksum();
        let gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        let mut archive = tar::Builder::new(gz);
        archive.append(&header, contents).unwrap();
        let bytes = archive.into_inner().unwrap().finish().unwrap();
        let destination = tempfile::tempdir().unwrap();
        assert!(
            extract_pack(
                bytes.as_slice(),
                destination.path(),
                1,
                &OfficialExternalDataOptions {
                    max_extracted_bytes: Some(max_bytes),
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn versions_and_package_paths_cannot_escape_the_installation() {
    for version in ["", "..", "../v1", "/tmp/v1", "v1/x", "v1?tag=x"] {
        assert!(version_component(version).is_err());
    }
    for path in ["../a", "/tmp/a", "a/../../b", "a\\b"] {
        assert!(data::relative_path(path).is_err());
    }
}

#[test]
fn interrupted_pack_installation_is_transactional_and_resumes_completed_packs() {
    let source = fixture();
    let version = fixture_manifest(source.path()).version;
    let complete_assets = release_assets(source.path());
    let mut interrupted_assets = complete_assets.clone();
    let second_pack = interrupted_assets
        .iter_mut()
        .find(|(name, _)| name.ends_with("-2L.tar.gz"))
        .unwrap()
        .1;
    second_pack.truncate(second_pack.len() / 2);
    let interrupted = Server::start(interrupted_assets);
    let destination = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Some(Duration::from_secs(3)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    assert!(
        engine
            .load_official_external_data_from(2, &options, &interrupted.url)
            .is_err()
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    let installations = destination.path().join(".vakint-external-data");
    let version_root = installations.join(&version);
    assert!(version_root.join("1L").is_dir());
    assert!(!version_root.join("2L").exists());
    assert!(!installations.join("latest").exists());
    assert!(fs::read_dir(&installations).unwrap().all(|entry| {
        !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".install-")
    }));
    drop(interrupted);
    let complete = Server::start(complete_assets);
    let report = engine
        .load_official_external_data_from(2, &options, &complete.url)
        .unwrap();
    assert_eq!(report.version, version);
    assert_eq!(report.loops[0].source, "reused");
    assert_eq!(report.loops[1].source, "downloaded");
}

#[test]
fn refresh_switches_all_requested_loops_together_and_preserves_the_old_installation() {
    let source = fixture();
    let original_version = fixture_manifest(source.path()).version;
    let refreshed_version = format!("{original_version}-refresh");
    let first = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Some(Duration::from_secs(3)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let old = engine
        .load_official_external_data_from(2, &options, &first.url)
        .unwrap();
    assert_eq!(old.version, original_version);
    retag_release(source.path(), &refreshed_version);
    let newer = Server::start(release_assets(source.path()));
    // Offline the existing installation keeps its cached release (with the release
    // server reachable it would follow data-latest: see the test below).
    let unchanged = engine
        .load_official_external_data_from(2, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(unchanged.version, original_version);
    assert_eq!(unchanged.path, old.path);
    let refreshed = engine
        .load_official_external_data_from(
            3,
            &OfficialExternalDataOptions {
                refresh: true,
                ..options.clone()
            },
            &newer.url,
        )
        .unwrap();
    assert_eq!(refreshed.version, refreshed_version);
    assert_eq!(refreshed.loops.len(), 3);
    assert!(
        refreshed
            .loops
            .iter()
            .all(|pack| pack.source == "downloaded")
    );
    assert!(
        engine
            .rustred_data()
            .unwrap()
            .loops
            .values()
            .all(|pack| pack.manifest.version == refreshed_version)
    );
    assert!(old.path.join("1L/ibps/rules.rrbin").is_file());
    let pinned = engine
        .load_official_external_data_from(
            2,
            &OfficialExternalDataOptions {
                version: original_version.clone(),
                ..options
            },
            "http://127.0.0.1:1",
        )
        .unwrap();
    assert_eq!(pinned.version, original_version);
    assert_eq!(pinned.path, old.path);
    assert!(pinned.loops.iter().all(|pack| pack.source == "reused"));
}

/// Re-declare the release in `root` as `version` (its root and loop-pack manifests and
/// the archive names), as the next published release would be.
fn retag_release(root: &Path, version: &str) {
    let mut manifest = fixture_manifest(root);
    manifest.version = version.to_owned();
    for pack in &mut manifest.loop_packs {
        let path = root.join(&pack.manifest);
        let mut payload: serde_json::Value =
            serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        payload["version"] = serde_json::json!(manifest.version);
        let bytes = serde_json::to_vec(&payload).unwrap();
        pack.manifest_sha256 = format!("{:x}", Sha256::digest(&bytes));
        pack.archive = format!(
            "vakint-external-data-{}-{}L.tar.gz",
            manifest.version, pack.loop_count
        );
        fs::write(path, bytes).unwrap();
    }
    fs::write(
        root.join("manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
}

#[test]
fn cached_latest_follows_a_moved_data_latest_and_falls_back_offline() {
    let source = fixture();
    let original_version = fixture_manifest(source.path()).version;
    let moved_version = format!("{original_version}-next");
    let first = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let installations = destination.path().join(".vakint-external-data");
    let marker = installations.join("latest");
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Some(Duration::from_secs(3)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let old = engine
        .load_official_external_data_from(2, &options, &first.url)
        .unwrap();
    assert_eq!(old.version, original_version);
    assert_eq!(fs::read_to_string(&marker).unwrap(), original_version);
    // data-latest still names the cached release: reused, nothing downloaded
    let same = engine
        .load_official_external_data_from(2, &options, &first.url)
        .unwrap();
    assert_eq!(same.path, old.path);
    assert!(same.loops.iter().all(|pack| pack.source == "reused"));
    drop(first);

    // A new release is published and data-latest moves to it: the next ordinary
    // load (no refresh) installs and activates it and moves the marker.
    retag_release(source.path(), &moved_version);
    let moved = Server::start(release_assets(source.path()));
    let followed = engine
        .load_official_external_data_from(2, &options, &moved.url)
        .unwrap();
    assert_eq!(followed.version, moved_version);
    assert_eq!(followed.path, installations.join(&moved_version));
    assert!(
        followed
            .loops
            .iter()
            .all(|pack| pack.source == "downloaded")
    );
    assert_eq!(fs::read_to_string(&marker).unwrap(), moved_version);
    let active = engine.rustred_data().unwrap();
    assert!((1..=2).all(|count| active.loops[&count].manifest.version == moved_version));
    assert!(old.path.join("1L/ibps/rules.rrbin").is_file());
    drop(moved);

    // A data-latest manifest that is no valid release is not followed.
    let broken = Server::start(BTreeMap::from([(
        "/broken/manifest.json".to_owned(),
        b"{}".to_vec(),
    )]));
    let kept = engine
        .load_official_external_data_from(2, &options, &broken.url)
        .unwrap();
    assert_eq!(kept.version, moved_version);
    assert!(kept.loops.iter().all(|pack| pack.source == "reused"));
    drop(broken);
    // Offline: the cached release, without network.
    let offline = engine
        .load_official_external_data_from(2, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(offline.version, moved_version);
    assert_eq!(offline.path, followed.path);
    assert!(offline.loops.iter().all(|pack| pack.source == "reused"));
    // An explicit version stays pinned and leaves the marker alone.
    let pinned = engine
        .load_official_external_data_from(
            2,
            &OfficialExternalDataOptions {
                version: original_version.clone(),
                ..options
            },
            "http://127.0.0.1:1",
        )
        .unwrap();
    assert_eq!(pinned.version, original_version);
    assert_eq!(pinned.path, old.path);
    assert_eq!(fs::read_to_string(&marker).unwrap(), moved_version);
}

#[test]
fn data_latest_check_is_one_bounded_attempt() {
    let source = fixture();
    let version = fixture_manifest(source.path()).version;
    let server = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let mut options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Some(Duration::from_secs(3)),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    engine
        .load_official_external_data_from(1, &options, &server.url)
        .unwrap();
    drop(server);
    // A release server that accepts connections and never answers.
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    listener.set_nonblocking(true).unwrap();
    let stop = Arc::new(AtomicBool::new(false));
    let accepted = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let handle = {
        let (stop, accepted) = (Arc::clone(&stop), Arc::clone(&accepted));
        thread::spawn(move || {
            let mut open = Vec::new();
            while !stop.load(Ordering::Relaxed) {
                match listener.accept() {
                    Ok((stream, _)) => {
                        accepted.fetch_add(1, Ordering::Relaxed);
                        open.push(stream);
                    }
                    Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                        thread::sleep(Duration::from_millis(5))
                    }
                    Err(error) => panic!("silent fixture server: {error}"),
                }
            }
        })
    };
    // The read-idle bound limits the whole check; the retry counter does not apply.
    options.read_idle_timeout = Some(Duration::from_millis(300));
    options.transfer_retries = 5;
    let started = std::time::Instant::now();
    let report = engine
        .load_official_external_data_from(1, &options, &url)
        .unwrap();
    let elapsed = started.elapsed();
    thread::sleep(Duration::from_millis(50));
    stop.store(true, Ordering::Relaxed);
    handle.join().unwrap();
    assert_eq!(report.version, version);
    assert_eq!(report.loops[0].source, "reused");
    assert_eq!(accepted.load(Ordering::Relaxed), 1, "one attempt, no retry");
    assert!(elapsed < Duration::from_secs(3), "check took {elapsed:?}");
}

#[test]
fn extraction_bounds_extension_metadata_before_tar_yields_file_entries() {
    let gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    let mut archive = tar::Builder::new(gz);
    let comment = vec![b'x'; 16 * 1024];
    archive
        .append_pax_extensions([("comment", comment.as_slice())])
        .unwrap();
    let mut header = tar::Header::new_gnu();
    header.set_path("1L/value").unwrap();
    header.set_size(1);
    header.set_mode(0o644);
    header.set_cksum();
    archive.append(&header, &b"1"[..]).unwrap();
    let bytes = archive.into_inner().unwrap().finish().unwrap();
    let destination = tempfile::tempdir().unwrap();
    assert!(
        extract_pack(
            bytes.as_slice(),
            destination.path(),
            1,
            &OfficialExternalDataOptions {
                max_extracted_bytes: Some(8 * 1024),
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(!destination.path().join("1L/value").exists());
}

// The bytes are package-authenticated but deliberately not a native cache.
// Vakint never reads prepared depth-first caches, so loading must not decode,
// import or even checksum them.
fn advertise_undecodable_cache(root: &Path, loops: usize) -> PathBuf {
    let path = root.join(format!("{loops}L/cache/cache.json"));
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    let bytes = b"authenticated bytes, not a native cache";
    fs::write(&path, bytes).unwrap();
    let manifest_path = root.join(format!("{loops}L/manifest.json"));
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&manifest_path).unwrap()).unwrap();
    manifest["cache"] = if loops == 5 {
        serde_json::json!({"available":true,"path":"cache/cache.json",
            "format":"rustred.saved-application-cache-file.v1",
            "identity":{"compressed_bytes":bytes.len(),"compressed_blake3":vec![0u8;32],
                "decoded_bytes":0,"decoded_blake3":vec![0u8;32]}})
    } else {
        serde_json::json!({"available":true,"path":"cache"})
    };
    manifest["files"]["cache/cache.json"] = serde_json::json!({
        "bytes":bytes.len(),"sha256":format!("{:x}",Sha256::digest(bytes))});
    let bytes = serde_json::to_vec(&manifest).unwrap();
    fs::write(&manifest_path, &bytes).unwrap();
    let release_path = root.join("manifest.json");
    let mut release: ReleaseManifest =
        serde_json::from_slice(&fs::read(&release_path).unwrap()).unwrap();
    release
        .loop_packs
        .iter_mut()
        .find(|pack| pack.loop_count == loops)
        .unwrap()
        .manifest_sha256 = format!("{:x}", Sha256::digest(bytes));
    fs::write(release_path, serde_json::to_vec(&release).unwrap()).unwrap();
    path
}

#[test]
fn advertised_prepared_caches_are_neither_decoded_nor_validated_at_load() {
    let directory = fixture();
    let payload = advertise_undecodable_cache(directory.path(), 1);
    let mut engine = Vakint::new().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(directory.path().to_owned()),
        ..Default::default()
    };
    let report = engine.load_official_external_data(1, &options).unwrap();
    assert_eq!(report.loops[0].source, "local");
    assert!(engine.rustred_application.lock().unwrap().is_none());
    // Tampered cache bytes are irrelevant to an ordered engine: loading skips
    // them, while an explicit read still authenticates every listed file.
    fs::write(&payload, b"tampered cache bytes").unwrap();
    engine.load_official_external_data(1, &options).unwrap();
    let error = engine.rustred_data().unwrap().loops[&1]
        .read("cache/cache.json")
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("incorrect size") || error.contains("SHA256 mismatch"),
        "{error}"
    );
    // Every other package file is still authenticated at load.
    fs::write(directory.path().join("1L/masters/P1/symbolic.json"), b"[]").unwrap();
    let active = engine.rustred_data().unwrap();
    assert!(engine.load_official_external_data(1, &options).is_err());
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
}

#[test]
#[ignore = "requires the authenticated 5L application pack; no IBP reduction"]
fn advertised_prepared_caches_are_ignored_through_five_loops() {
    let source = super::super::test_data::external_path_through(5);
    let directory = tempfile::tempdir().unwrap();
    fs::copy(
        source.join("manifest.json"),
        directory.path().join("manifest.json"),
    )
    .unwrap();
    for loops in 1..=5 {
        copy_tree(
            &source.join(format!("{loops}L")),
            &directory.path().join(format!("{loops}L")),
        );
    }
    // Both ends of the supported loop range advertise undecodable, correctly
    // checksummed cache payloads. Scientific programs/catalogues are unchanged.
    advertise_undecodable_cache(directory.path(), 1);
    advertise_undecodable_cache(directory.path(), 5);
    let mut engine = Vakint::new().unwrap();
    let report = engine
        .load_official_external_data(
            5,
            &OfficialExternalDataOptions {
                local_artifacts_path: Some(directory.path().to_owned()),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(report.loops.len(), 5);
    assert!(engine.loaded_application().is_ok());
}

#[test]
fn official_transfer_options_have_no_default_limits() {
    // No arbitrary limits by default: no request timeout, no archive, extraction or
    // entry cap (were 3600 s, 2 GiB, 8 GiB and 250000 entries).
    let options = OfficialExternalDataOptions::default();
    assert_eq!(options.timeout, None);
    // progress-based detector and retry counter instead of a whole-transfer timeout
    assert_eq!(options.read_idle_timeout, Some(Duration::from_secs(300)));
    assert_eq!(options.transfer_retries, 5);
    assert_eq!(options.max_archive_bytes, None);
    assert_eq!(options.max_extracted_bytes, None);
    assert_eq!(options.max_entries, None);
}

/// Serves `assets`, one thread per connection; the first request for each archive sends its
/// headers and half of the body, then stalls longer than the client's read-idle detector.
fn start_stalling_server(
    assets: BTreeMap<String, Vec<u8>>,
) -> (String, Arc<AtomicBool>, thread::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let url = format!("http://{}", listener.local_addr().unwrap());
    listener.set_nonblocking(true).unwrap();
    let stop = Arc::new(AtomicBool::new(false));
    let stopped = Arc::clone(&stop);
    // like Server::start: the default version "latest" asks for /data-latest/manifest.json, answered by a redirect
    let latest_manifest = assets
        .keys()
        .rfind(|path| path.ends_with("/manifest.json"))
        .unwrap()
        .clone();
    let assets = Arc::new(assets);
    let stalled = Arc::new(std::sync::Mutex::new(std::collections::BTreeSet::new()));
    let handle = thread::spawn(move || {
        while !stopped.load(Ordering::Relaxed) {
            match listener.accept() {
                Ok((mut stream, _)) => {
                    let (assets, stalled) = (Arc::clone(&assets), Arc::clone(&stalled));
                    let latest_manifest = latest_manifest.clone();
                    thread::spawn(move || {
                        stream.set_nonblocking(false).unwrap();
                        stream
                            .set_read_timeout(Some(Duration::from_secs(2)))
                            .unwrap();
                        let mut request = [0u8; 4096];
                        let Ok(n) = stream.read(&mut request) else {
                            return;
                        };
                        let text = String::from_utf8_lossy(&request[..n]);
                        let Some(path) = text.split_whitespace().nth(1).map(str::to_owned) else {
                            return;
                        };
                        if path == "/data-latest/manifest.json" {
                            let _ = write!(
                                stream,
                                "HTTP/1.1 302 Found\r\nLocation: {latest_manifest}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
                            );
                            return;
                        }
                        match assets.get(&path) {
                            Some(bytes)
                                if path.ends_with(".tar.gz")
                                    && stalled.lock().unwrap().insert(path.clone()) =>
                            {
                                let _ = write!(
                                    stream,
                                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                    bytes.len()
                                );
                                let _ = stream.write_all(&bytes[..bytes.len() / 2]);
                                let _ = stream.flush();
                                thread::sleep(Duration::from_millis(800));
                            }
                            Some(bytes) => {
                                let _ = write!(
                                    stream,
                                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                                    bytes.len()
                                );
                                let _ = stream.write_all(bytes);
                            }
                            None => {
                                let _ = stream.write_all(b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n");
                            }
                        }
                    });
                }
                Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                    thread::sleep(Duration::from_millis(5))
                }
                Err(error) => panic!("stalling fixture server: {error}"),
            }
        }
    });
    (url, stop, handle)
}

#[test]
fn stalled_transfer_is_detected_by_read_idle_and_retried() {
    // NL2: no whole-transfer deadline; a body that stops arriving trips the read-idle detector
    // and the transfer is retried (the second request of each archive is served in full).
    let source = fixture();
    let (url, stop, handle) = start_stalling_server(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        read_idle_timeout: Some(Duration::from_millis(200)),
        transfer_retries: 2,
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let report = engine
        .load_official_external_data_from(2, &options, &url)
        .unwrap();
    assert!(report.loops.iter().all(|p| p.source == "downloaded"));
    // without retries the same stall is a transport error
    let again = tempfile::tempdir().unwrap();
    let (url2, stop2, handle2) = start_stalling_server(release_assets(source.path()));
    let no_retry = OfficialExternalDataOptions {
        local_artifacts_path: Some(again.path().to_owned()),
        read_idle_timeout: Some(Duration::from_millis(200)),
        transfer_retries: 0,
        ..Default::default()
    };
    assert!(
        Vakint::new()
            .unwrap()
            .load_official_external_data_from(2, &no_retry, &url2)
            .is_err()
    );
    for (s, h) in [(stop, handle), (stop2, handle2)] {
        s.store(true, Ordering::Relaxed);
        h.join().unwrap();
    }
}
