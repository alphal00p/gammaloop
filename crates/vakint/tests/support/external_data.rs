//! The vendored manifest also pins the compatible immutable test release.

pub(crate) fn required_version() -> String {
    #[derive(serde::Deserialize)]
    struct ManifestVersion {
        version: String,
    }
    serde_json::from_slice::<ManifestVersion>(include_bytes!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/manifest.json"
    )))
    .expect("vendored data must declare the compatible test release")
    .version
}
