//! Shared external-data fixtures for unit tests that need the 4L/5L packs.
use crate::Vakint;
use std::path::PathBuf;

#[cfg(all(unix, not(target_arch = "wasm32")))]
#[path = "../../tests/support/external_data.rs"]
pub(super) mod external_data_fixture;

/// An explicit checkout overrides the official release used by these tests.
/// Otherwise install/reuse the vendored manifest's compatible immutable release
/// in the caller's current working directory, independently of cached `latest`.
/// Ordinary engine construction remains offline and vendors only 1L–3L.
pub(super) fn external_path() -> PathBuf {
    external_path_through(4)
}

pub(super) fn external_path_through(max_loop_count: usize) -> PathBuf {
    if let Some(path) = std::env::var_os("VAKINT_TEST_EXTERNAL_DATA") {
        return PathBuf::from(path);
    }
    #[cfg(all(unix, not(target_arch = "wasm32")))]
    {
        static PATHS: std::sync::Mutex<std::collections::BTreeMap<usize, PathBuf>> =
            std::sync::Mutex::new(std::collections::BTreeMap::new());
        let mut paths = PATHS.lock().unwrap();
        if let Some((_, path)) = paths.range(max_loop_count..).next() {
            return path.clone();
        }
        let mut engine = Vakint::new().expect("initialize the official-data test loader");
        let options = super::OfficialExternalDataOptions {
            version: external_data_fixture::required_version(),
            ..Default::default()
        };
        let report = engine
            .load_official_external_data(max_loop_count, &options)
            .expect("install or reuse official external data for native asset tests");
        eprintln!(
            "Vakint native test data: version={} path={}",
            report.version,
            report.path.display()
        );
        paths.insert(max_loop_count, report.path.clone());
        report.path
    }
    #[cfg(not(all(unix, not(target_arch = "wasm32"))))]
    {
        let _ = max_loop_count;
        panic!(
            "automatic official-data fixtures require a native Unix target; set VAKINT_TEST_EXTERNAL_DATA to use a local checkout"
        )
    }
}

/// Share immutable program ownership (validated 1L–4L data) between tests while
/// giving each test a fresh engine.
pub(super) fn external_engine() -> Vakint {
    static DATA: std::sync::OnceLock<std::sync::Arc<super::data::RustRedData>> =
        std::sync::OnceLock::new();
    let data = DATA.get_or_init(|| {
        super::data::RustRedData::from_directory(&external_path(), "latest", 4)
            .expect("external test data must validate and load")
    });
    let mut engine = Vakint::new().unwrap();
    engine.rustred_data = Some(std::sync::Arc::clone(data));
    engine
}
