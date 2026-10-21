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
        timeout: Duration::from_millis(100),
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
    assert!(
        report
            .loops
            .iter()
            .all(|p| p.source == "local" && !p.cache_available)
    );
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

#[test]
fn authenticated_pack_with_wrong_cache_binding_preserves_active_data_and_cache() {
    let directory = fixture();
    let prepared = tempfile::tempdir().unwrap();
    let cache = prepared.path().join("cache");
    let mut engine = Vakint::new().unwrap();
    engine
        .prepare_cache(
            1,
            &cache,
            &crate::RustRedCachePreparationOptions {
                workload: crate::RustRedCacheWorkload::Explicit(vec![crate::RustRedCacheRoot {
                    namespace: 0,
                    integral: rustred::family::IntegralKey::try_new([2]).unwrap(),
                }]),
                ..Default::default()
            },
            &AtomicBool::new(false),
            |_| Ok(()),
        )
        .unwrap();
    engine.load_cache(1, &cache, &Default::default()).unwrap();
    let active_data = engine.rustred_data().unwrap();
    let active_caches = Arc::clone(&engine.rustred_lower_caches);
    let statistics = engine.cache_statistics(1).unwrap();
    assert!(statistics.attached);

    // This is a valid native 1L cache with valid package checksums. Only the
    // later attachment to the declared 2L program must reject its identity.
    let loop_root = directory.path().join("2L");
    copy_tree(&cache, &loop_root.join("cache"));
    let loop_manifest_path = loop_root.join("manifest.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&fs::read(&loop_manifest_path).unwrap()).unwrap();
    manifest["cache"] = serde_json::json!({"available": true, "path": "cache"});
    let mut pending = vec![loop_root.join("cache")];
    while let Some(path) = pending.pop() {
        for entry in fs::read_dir(path).unwrap() {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_dir() {
                pending.push(entry.path());
            } else {
                let bytes = fs::read(entry.path()).unwrap();
                let relative = entry.path().strip_prefix(&loop_root).unwrap().to_owned();
                manifest["files"][relative.to_str().unwrap()] = serde_json::json!({
                    "sha256": format!("{:x}", Sha256::digest(&bytes)), "bytes": bytes.len()
                });
            }
        }
    }
    let bytes = serde_json::to_vec(&manifest).unwrap();
    fs::write(&loop_manifest_path, &bytes).unwrap();
    let release_path = directory.path().join("manifest.json");
    let mut release: ReleaseManifest =
        serde_json::from_slice(&fs::read(&release_path).unwrap()).unwrap();
    let pack = release
        .loop_packs
        .iter_mut()
        .find(|p| p.loop_count == 2)
        .unwrap();
    pack.manifest_sha256 = format!("{:x}", Sha256::digest(&bytes));
    data::validate_loop_directory(directory.path(), pack, &release.version).unwrap();
    fs::write(release_path, serde_json::to_vec(&release).unwrap()).unwrap();
    let error = engine
        .load_official_external_data_from(
            2,
            &OfficialExternalDataOptions {
                local_artifacts_path: Some(directory.path().to_owned()),
                ..Default::default()
            },
            "http://127.0.0.1:1",
        )
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("cache bundle program bindings differ"),
        "{error}"
    );
    assert!(Arc::ptr_eq(&active_data, &engine.rustred_data().unwrap()));
    assert!(Arc::ptr_eq(&active_caches, &engine.rustred_lower_caches));
    assert_eq!(statistics, engine.cache_statistics(1).unwrap());
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
                timeout: Duration::from_secs(3),
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
        timeout: Duration::from_secs(3),
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
    let required = super::super::cache_tests::external_data_fixture::required_version();
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
        timeout: Duration::from_secs(3),
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
                    max_extracted_bytes: max_bytes,
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
        timeout: Duration::from_secs(3),
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
    let mut manifest = fixture_manifest(source.path());
    let original_version = manifest.version.clone();
    let refreshed_version = format!("{original_version}-refresh");
    let first = Server::start(release_assets(source.path()));
    let destination = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(destination.path().to_owned()),
        timeout: Duration::from_secs(3),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let old = engine
        .load_official_external_data_from(2, &options, &first.url)
        .unwrap();
    assert_eq!(old.version, original_version);
    manifest.version = refreshed_version.clone();
    for pack in &mut manifest.loop_packs {
        let path = source.path().join(&pack.manifest);
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
        source.path().join("manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    let newer = Server::start(release_assets(source.path()));
    let unchanged = engine
        .load_official_external_data_from(2, &options, &newer.url)
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
                max_extracted_bytes: 8 * 1024,
                ..Default::default()
            }
        )
        .is_err()
    );
    assert!(!destination.path().join("1L/value").exists());
}
