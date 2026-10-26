//! Tiny real native records exercise the same official path used by 5L.
use super::*;
use crate::rustred_evaluation::cache::transfer::tests::MultipartFixture;
use serde_json::{Value, json};

struct NativeRelease {
    directory: tempfile::TempDir,
    native: MultipartFixture,
}

impl NativeRelease {
    fn new() -> Self {
        let release = Self {
            directory: fixture(),
            native: MultipartFixture::new(),
        };
        release.write_descriptor(&release.native.manifest);
        release
    }

    fn write_descriptor(&self, descriptor: &Value) {
        let root = self.directory.path();
        let relative = "cache/transport.json";
        fs::create_dir_all(root.join("1L/cache")).unwrap();
        let bytes = serde_json::to_vec(descriptor).unwrap();
        fs::write(root.join("1L").join(relative), &bytes).unwrap();
        let identity =
            json!({"sha256": format!("{:x}", Sha256::digest(&bytes)), "bytes": bytes.len()});
        let loop_path = root.join("1L/manifest.json");
        let mut pack: Value = serde_json::from_slice(&fs::read(&loop_path).unwrap()).unwrap();
        pack["files"][relative] = identity.clone();
        let loop_bytes = serde_json::to_vec(&pack).unwrap();
        fs::write(loop_path, &loop_bytes).unwrap();
        let mut manifest = fixture_manifest(root);
        manifest
            .loop_packs
            .iter_mut()
            .find(|p| p.loop_count == 1)
            .unwrap()
            .manifest_sha256 = format!("{:x}", Sha256::digest(&loop_bytes));
        manifest.runtime_cache_snapshots = vec![data::RuntimeCacheReference {
            loop_count: 1,
            manifest: "1L/cache/transport.json".into(),
            manifest_sha256: identity["sha256"].as_str().unwrap().into(),
            manifest_bytes: bytes.len() as u64,
        }];
        fs::write(
            root.join("manifest.json"),
            serde_json::to_vec(&manifest).unwrap(),
        )
        .unwrap();
    }

    fn assets(&self) -> BTreeMap<String, Vec<u8>> {
        let root = self.directory.path();
        let mut assets = release_assets(root);
        let version = fixture_manifest(root).version;
        assets.insert(
            format!("/{version}/transport.json"),
            fs::read(root.join("1L/cache/transport.json")).unwrap(),
        );
        for (n, part) in self.native.manifest["parts"]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            assets.insert(
                format!("/{version}/{}", part["archive"].as_str().unwrap()),
                fs::read(self.native.part_path(n)).unwrap(),
            );
        }
        assets
    }

    fn native_path(&self, root: &Path) -> PathBuf {
        root.join("cache-snapshots").join(
            &fixture_manifest(self.directory.path()).runtime_cache_snapshots[0].manifest_sha256,
        )
    }
}

#[test]
fn official_nativepack_installs_then_reuses_without_network_or_static_cache() {
    let release = NativeRelease::new();
    let server = Server::start(release_assets(release.directory.path()));
    let installation = tempfile::tempdir().unwrap();
    let mut engine = Vakint::new().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(installation.path().into()),
        ..Default::default()
    };
    let disabled = OfficialExternalDataOptions {
        load_prepared_caches: false,
        ..options.clone()
    };
    let cold = engine
        .load_official_external_data_from(1, &disabled, &server.url)
        .unwrap();
    assert!(cold.loops[0].cache_available);
    assert!(!cold.loops[0].cache_loaded);
    assert!(!engine.cache_statistics(1).unwrap().attached);
    assert!(!cold.loops[0].cache_path.as_ref().unwrap().exists());
    // Parts are deliberately unavailable during the disabled load above: a
    // separate server now exposes them for the first real cache installation.
    drop(server);
    let server = Server::start(release.assets());
    let loaded = engine
        .load_official_external_data_from(1, &options, &server.url)
        .unwrap();
    let item = &loaded.loops[0];
    assert_eq!(item.source, "reused");
    assert_eq!(item.cache_source.as_deref(), Some("downloaded"));
    assert_eq!(
        item.cache_format.as_deref(),
        Some("vakint.nativepack-transport.v1")
    );
    assert!(item.cache_loaded && engine.cache_statistics(1).unwrap().attached);
    release
        .native
        .fixture
        .assert_imported_records(item.cache_path.as_ref().unwrap());
    let expected = release.native_path(&loaded.path);
    assert_eq!(item.cache_path.as_ref(), Some(&expected));
    drop(server);
    let reused = engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(reused.loops[0].cache_source.as_deref(), Some("reused"));
    assert_eq!(reused.loops[0].cache_path, item.cache_path);
    let no_cache = engine
        .load_official_external_data_from(1, &disabled, "http://127.0.0.1:1")
        .unwrap();
    assert!(!no_cache.loops[0].cache_loaded);
    assert!(!engine.cache_statistics(1).unwrap().attached);
    assert!(expected.exists());
}

#[test]
fn official_nativepack_checkout_never_downloads_or_mutates_missing_snapshot() {
    let release = NativeRelease::new();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(release.directory.path().into()),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    let error = engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("selected official native cache is missing"),
        "{error}"
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    assert!(!release.directory.path().join("cache-snapshots").exists());
    assert!(
        !release
            .directory
            .path()
            .join(".vakint-external-data")
            .exists()
    );
    let target = release.native_path(release.directory.path());
    copy_tree(&release.native.fixture.source, &target);
    let loaded = engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(loaded.loops[0].cache_source.as_deref(), Some("local"));
    assert!(loaded.loops[0].cache_loaded);
    let active = engine.rustred_data().unwrap();
    let caches = Arc::clone(&engine.rustred_lower_caches);
    let index = target.join("bundle.rrindex");
    let mut bytes = fs::read(&index).unwrap();
    bytes[0] ^= 1;
    fs::write(index, bytes).unwrap();
    assert!(
        engine
            .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
            .is_err()
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    assert!(Arc::ptr_eq(&caches, &engine.rustred_lower_caches));
}

#[test]
fn official_nativepack_failed_download_preserves_engine_and_latest_marker() {
    let release = NativeRelease::new();
    let mut assets = release.assets();
    let part = assets
        .keys()
        .find(|p| p.ends_with("/chunks.tar.gz"))
        .unwrap()
        .clone();
    assets.get_mut(&part).unwrap().pop();
    let server = Server::start(assets);
    let installation = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(installation.path().into()),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    let cache = Arc::clone(&engine.rustred_lower_caches);
    let error = engine
        .load_official_external_data_from(1, &options, &server.url)
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("archive length or SHA-256 mismatch"),
        "{error}"
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    assert!(Arc::ptr_eq(&cache, &engine.rustred_lower_caches));
    assert!(
        !installation
            .path()
            .join(".vakint-external-data/latest")
            .exists()
    );
    let version = fixture_manifest(release.directory.path()).version;
    let root = installation
        .path()
        .join(".vakint-external-data")
        .join(version);
    assert!(!release.native_path(&root).exists());
    assert!(
        fs::read_dir(root.join("cache-snapshots"))
            .unwrap()
            .next()
            .is_none()
    );
}

#[test]
fn official_nativepack_binding_and_reference_failures_are_transactional() {
    let release = NativeRelease::new();
    let mut descriptor = release.native.manifest.clone();
    descriptor["binding"]["family_fingerprint"] = "wrong-program".into();
    release.write_descriptor(&descriptor);
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    let disabled = OfficialExternalDataOptions {
        local_artifacts_path: Some(release.directory.path().into()),
        load_prepared_caches: false,
        ..Default::default()
    };
    let error = engine
        .load_official_external_data_from(1, &disabled, "http://127.0.0.1:1")
        .unwrap_err()
        .to_string();
    assert!(error.contains("authenticated program"), "{error}");
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    release.write_descriptor(&release.native.manifest);
    let root = release.directory.path().join("manifest.json");
    let original = fixture_manifest(release.directory.path());
    let mut duplicate = original.clone();
    duplicate
        .runtime_cache_snapshots
        .push(duplicate.runtime_cache_snapshots[0].clone());
    fs::write(&root, serde_json::to_vec(&duplicate).unwrap()).unwrap();
    let error = engine
        .load_official_external_data_from(1, &disabled, "http://127.0.0.1:1")
        .unwrap_err()
        .to_string();
    assert!(error.contains("multiple selected"), "{error}");
    let mut mismatch = original;
    mismatch.runtime_cache_snapshots[0].manifest_bytes += 1;
    fs::write(root, serde_json::to_vec(&mismatch).unwrap()).unwrap();
    let error = engine
        .load_official_external_data_from(1, &disabled, "http://127.0.0.1:1")
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("differs from loop file inventory"),
        "{error}"
    );
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
}

#[test]
fn official_nativepack_inventory_grants_scale_past_old_aggregate_limits() {
    let release = NativeRelease::new();
    let mut value = release.native.manifest.clone();
    let files = value["files"].as_object_mut().unwrap();
    for file in files.values_mut() {
        file["bytes"] = (40u64 * 1024 * 1024 * 1024).into();
    }
    let parsed = super::super::super::cache::transfer::NativePackManifest::parse(
        &serde_json::to_vec(&value).unwrap(),
        1,
    )
    .unwrap();
    let grants = parsed
        .official_options([1; 32], 2048, Duration::from_secs(60))
        .unwrap();
    assert!(grants.max_extracted_bytes > 32 * 1024 * 1024 * 1024);
    assert!(grants.max_temporary_bytes > grants.max_extracted_bytes * 2);
    assert_eq!(grants.max_manifest_bytes, 2048);
    for file in value["files"].as_object_mut().unwrap().values_mut() {
        file["bytes"] = u64::MAX.into();
    }
    let parsed = super::super::super::cache::transfer::NativePackManifest::parse(
        &serde_json::to_vec(&value).unwrap(),
        1,
    )
    .unwrap();
    assert!(
        parsed
            .official_options([1; 32], 2048, Duration::from_secs(60))
            .is_err()
    );
}

#[test]
fn official_nativepack_reuse_checks_native_generation_and_reference_census() {
    let release = NativeRelease::new();
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(release.directory.path().into()),
        ..Default::default()
    };
    for field in ["generation", "stored_record_references"] {
        let mut descriptor = release.native.manifest.clone();
        descriptor[field] = (descriptor[field].as_u64().unwrap() + 1).into();
        release.write_descriptor(&descriptor);
        let target = release.native_path(release.directory.path());
        copy_tree(&release.native.fixture.source, &target);
        let error = engine
            .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("generation/reference census differs"),
            "{field}: {error}"
        );
        assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
        assert!(
            !release
                .directory
                .path()
                .join(".vakint-external-data")
                .exists()
        );
    }
}

#[test]
fn official_nativepack_reuse_rejects_a_symlinked_parent_without_http() {
    let release = NativeRelease::new();
    let outside = tempfile::tempdir().unwrap();
    let target = release.native_path(release.directory.path());
    let outside_target = outside.path().join(target.file_name().unwrap());
    copy_tree(&release.native.fixture.source, &outside_target);
    std::os::unix::fs::symlink(outside.path(), target.parent().unwrap()).unwrap();
    let mut engine = Vakint::new().unwrap();
    let active = engine.rustred_data().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(release.directory.path().into()),
        ..Default::default()
    };
    let error = engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap_err()
        .to_string();
    assert!(error.contains("snapshot parent directory"), "{error}");
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    assert!(target.parent().unwrap().is_symlink());
    assert!(
        !release
            .directory
            .path()
            .join(".vakint-external-data")
            .exists()
    );
}

#[test]
#[ignore = "requires an authenticated 640-master 5L TEST copy with the preserved 1089-record selected NativePack; no IBP reduction"]
fn official_five_loop_selected_running_snapshot_attaches_with_lower_loop_caches() {
    let root = PathBuf::from(
        std::env::var("VAKINT_TEST_OFFICIAL_NATIVEPACK_DATA")
            .expect("provide the prepared isolated 5L selected-snapshot test copy"),
    );
    let mut engine = Vakint::new().unwrap();
    let report = engine
        .load_official_external_data(
            5,
            &OfficialExternalDataOptions {
                local_artifacts_path: Some(root),
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(report.loops.len(), 5);
    for item in &report.loops {
        assert!(item.cache_available && item.cache_loaded);
        assert_eq!(item.cache_source.as_deref(), Some("local"));
        assert!(engine.cache_statistics(item.loop_count).unwrap().attached);
    }
    assert_eq!(
        report.loops[4].cache_format.as_deref(),
        Some("vakint.nativepack-transport.v1")
    );
    let backing = engine
        .cache_statistics(5)
        .unwrap()
        .saved_application
        .unwrap();
    assert!(backing.file.is_none());
    assert_eq!(backing.pack.unwrap().records, 1089);
    let data = engine.rustred_data().unwrap();
    assert_eq!(data.loops[&5].manifest.masters.len(), 640);
    let snapshot = &data.runtime_cache_snapshots[&5];
    let descriptor: Value =
        serde_json::from_slice(&fs::read(report.path.join(&snapshot.reference.manifest)).unwrap())
            .unwrap();
    assert_eq!(descriptor["generation"], 19);
    assert_eq!(descriptor["stored_record_references"], 1089);
    assert_eq!(
        report.loops[4].cache_path.as_ref(),
        Some(&data.runtime_cache_path(5).unwrap())
    );
    println!(
        "{}",
        json!({"status":"PASS_CURRENT_FIVE_LOOP_SELECTED_NATIVEPACK_ATTACHMENT", "catalogue_entries":640, "stored_record_references":1089, "generation":19, "lower_loop_caches_attached":true, "ibp_reduction_performed":false})
    );
}

#[test]
fn official_nativepack_selection_is_part_of_registry_directory_identity() {
    let release = NativeRelease::new();
    let root = release.directory.path();
    let data = RustRedData::from_directory(root, "latest", 3).unwrap();
    assert!(data.matches_directory(root, 3).unwrap());
    let mut manifest = fixture_manifest(root);
    manifest.runtime_cache_snapshots[0].manifest_sha256 = "0".repeat(64);
    fs::write(
        root.join("manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    assert!(data.matches_directory(root, 3).is_err());
    manifest.runtime_cache_snapshots.clear();
    fs::write(
        root.join("manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    assert!(!data.matches_directory(root, 3).unwrap());
}

#[test]
fn official_nativepack_failed_refresh_preserves_old_cache_and_offline_alias() {
    let release = NativeRelease::new();
    let server = Server::start(release.assets());
    let installation = tempfile::tempdir().unwrap();
    let options = OfficialExternalDataOptions {
        local_artifacts_path: Some(installation.path().into()),
        ..Default::default()
    };
    let mut engine = Vakint::new().unwrap();
    let old = engine
        .load_official_external_data_from(1, &options, &server.url)
        .unwrap();
    drop(server);
    let active = engine.rustred_data().unwrap();
    let caches = Arc::clone(&engine.rustred_lower_caches);
    let latest_path = installation.path().join(".vakint-external-data/latest");
    let marker = fs::read(&latest_path).unwrap();
    let mut manifest = fixture_manifest(release.directory.path());
    manifest.version.push_str("-new-selected-cache");
    for pack in &mut manifest.loop_packs {
        let path = release.directory.path().join(&pack.manifest);
        let mut payload: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        payload["version"] = manifest.version.clone().into();
        let bytes = serde_json::to_vec(&payload).unwrap();
        fs::write(path, &bytes).unwrap();
        pack.manifest_sha256 = format!("{:x}", Sha256::digest(bytes));
        pack.archive = format!(
            "vakint-external-data-{}-{}L.tar.gz",
            manifest.version, pack.loop_count
        );
    }
    fs::write(
        release.directory.path().join("manifest.json"),
        serde_json::to_vec(&manifest).unwrap(),
    )
    .unwrap();
    let mut assets = release.assets();
    let part = assets
        .keys()
        .find(|p| p.ends_with("/chunks.tar.gz"))
        .unwrap()
        .clone();
    assets.get_mut(&part).unwrap().pop();
    let server = Server::start(assets);
    let refresh = OfficialExternalDataOptions {
        refresh: true,
        ..options.clone()
    };
    assert!(
        engine
            .load_official_external_data_from(1, &refresh, &server.url)
            .is_err()
    );
    drop(server);
    assert_eq!(fs::read(latest_path).unwrap(), marker);
    assert!(Arc::ptr_eq(&active, &engine.rustred_data().unwrap()));
    assert!(Arc::ptr_eq(&caches, &engine.rustred_lower_caches));
    let reused = engine
        .load_official_external_data_from(1, &options, "http://127.0.0.1:1")
        .unwrap();
    assert_eq!(reused.version, old.version);
    assert_eq!(reused.loops[0].cache_path, old.loops[0].cache_path);
    assert_eq!(reused.loops[0].cache_source.as_deref(), Some("reused"));
}
