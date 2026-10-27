//! Shared immutable file or Native-index identity for cache collectors.
use eyre::{Result, ensure};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::Read,
    path::Path,
};

#[derive(Debug, PartialEq, Eq)]
pub(super) enum CacheIdentity {
    File(String),
    Native(Value),
}

impl CacheIdentity {
    pub(super) fn read(path: &Path) -> Result<Self> {
        // Preserve the established single-file behavior, including file symlinks.
        if fs::metadata(path)?.is_file() {
            return Ok(Self::File(digest(path)?));
        }
        ensure!(
            fs::symlink_metadata(path)?.is_dir(),
            "cache must be a regular file or native directory"
        );
        let bundle = path.join("bundle.rrindex");
        let pack = path.join("cache.rrindex");
        // Dangling index symlinks are present invalid inputs, not a fallback.
        let bundle = bundle.is_symlink() || bundle.try_exists()?;
        let pack = pack.is_symlink() || pack.try_exists()?;
        ensure!(
            bundle != pack,
            "native cache needs exactly one bundle or pack index"
        );
        let mut indexes = BTreeMap::new();
        if bundle {
            Self::index(path, "bundle.rrindex", &mut indexes)?;
            let bindings = path.join("bindings");
            ensure!(
                fs::symlink_metadata(&bindings)?.is_dir(),
                "native bindings must be a regular directory"
            );
            for entry in fs::read_dir(bindings)? {
                let entry = entry?;
                let name = entry.file_name();
                let name = name
                    .to_str()
                    .ok_or_else(|| eyre::eyre!("invalid native namespace"))?;
                ensure!(
                    name.len() == 64
                        && name
                            .bytes()
                            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                        && entry.file_type()?.is_dir(),
                    "native namespace must be a lowercase hash directory"
                );
                Self::index(
                    path,
                    &format!("bindings/{name}/cache.rrindex"),
                    &mut indexes,
                )?;
            }
            ensure!(indexes.len() > 1, "native bundle has no child indexes");
        } else {
            Self::index(path, "cache.rrindex", &mut indexes)?;
        }
        // Keep the aggregate stable even if a host enables serde_json/preserve_order.
        let mut document = BTreeMap::from([
            ("schema", json!("vakint.native-cache-index-identity.v1")),
            ("layout", json!(if bundle { "bundle" } else { "pack" })),
            ("indexes", json!(indexes)),
        ]);
        let sha256 = format!("{:x}", Sha256::digest(serde_json::to_vec(&document)?));
        document.insert("sha256", sha256.into());
        Ok(Self::Native(json!(document)))
    }

    fn index(root: &Path, relative: &str, indexes: &mut BTreeMap<String, Value>) -> Result<()> {
        let path = root.join(relative);
        let metadata = fs::symlink_metadata(&path)?;
        ensure!(
            metadata.is_file(),
            "native index must be a regular file: {relative}"
        );
        indexes.insert(
            relative.into(),
            json!({"bytes":metadata.len(),"sha256":digest(&path)?}),
        );
        Ok(())
    }
}

pub(super) fn digest(path: &Path) -> Result<String> {
    let mut input = File::open(path)?;
    let mut digest = Sha256::new();
    let mut block = [0; 65536];
    loop {
        let count = input.read(&mut block)?;
        if count == 0 {
            break;
        }
        digest.update(&block[..count]);
    }
    Ok(format!("{:x}", digest.finalize()))
}
