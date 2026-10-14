//! External assets belong to the host; consumers retain the expected hashes.
//! The same interface reads files in native/Pyodide filesystems or accepts bytes
//! supplied by an async host. Missing or damaged assets never enter the cache.

use std::{
    collections::BTreeMap,
    fs,
    io::{self, Read},
    path::{Path, PathBuf},
    sync::{Arc, Mutex, OnceLock},
};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

#[derive(Default)]
pub struct Data {
    directories: Vec<PathBuf>,
    bytes: BTreeMap<(&'static str, &'static str), Arc<[u8]>>,
}

impl Data {
    fn shared() -> &'static Mutex<Self> {
        static DATA: OnceLock<Mutex<Data>> = OnceLock::new();
        DATA.get_or_init(|| {
            Mutex::new(Self {
                directories: std::env::var_os("SYMBOLICA_DATA_DIRECTORY")
                    .map(|paths| std::env::split_paths(&paths).collect())
                    .unwrap_or_default(),
                ..Self::default()
            })
        })
    }

    /// Select an extracted data bundle without reading its payloads.
    pub fn set_directory(directory: impl Into<PathBuf>) {
        Self::set_directories(vec![directory.into()]);
    }

    /// Search separately installed resource bundles without copying their files.
    pub fn set_directories(directories: Vec<PathBuf>) {
        Self::shared().lock().unwrap().directories = directories;
    }
}

/// A consumer's immutable expectation, generated from its own source assets.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct Asset {
    pub name: &'static str,
    pub sha256: &'static str,
    pub size: u64,
    pub description: &'static str,
}

impl Asset {
    /// Check a dataset before entering an infallible lazy mathematical cache.
    pub fn prepare(assets: &[Self], prefix: &str) -> io::Result<()> {
        for asset in assets.iter().filter(|asset| asset.name.starts_with(prefix)) {
            asset.bytes()?;
        }
        Ok(())
    }

    fn validate(&self, bytes: &[u8]) -> io::Result<()> {
        if bytes.len() as u64 != self.size || format!("{:x}", Sha256::digest(bytes)) != self.sha256
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "data asset {} differs from the compiled expectation",
                    self.name
                ),
            ));
        }
        Ok(())
    }

    /// Verify host-supplied bytes before any mathematical/native decoding.
    pub fn install(&self, bytes: Vec<u8>) -> io::Result<()> {
        self.validate(&bytes)?;
        Data::shared()
            .lock()
            .unwrap()
            .bytes
            .insert((self.name, self.sha256), bytes.into());
        Ok(())
    }

    pub fn is_loaded(&self) -> bool {
        Data::shared()
            .lock()
            .unwrap()
            .bytes
            .contains_key(&(self.name, self.sha256))
    }

    /// Missing data can be supplied later and the same request retried.
    pub fn bytes(&self) -> io::Result<Arc<[u8]>> {
        let mut data = Data::shared().lock().unwrap();
        if let Some(bytes) = data.bytes.get(&(self.name, self.sha256)) {
            return Ok(bytes.clone());
        }
        let file = data.directories.iter().find_map(|root| {
            match fs::File::open(root.join(self.name)) {
                Err(error) if error.kind() == io::ErrorKind::NotFound => None,
                result => Some(result),
            }
        }).unwrap_or_else(|| Err(io::Error::new(io::ErrorKind::NotFound, format!(
            "data asset {} is unavailable; install its resource wheel or configure a data directory", self.name
        ))))?;
        let mut bytes = Vec::new();
        file.take(self.size + 1).read_to_end(&mut bytes)?;
        self.validate(&bytes)?;
        let bytes: Arc<[u8]> = bytes.into();
        data.bytes.insert((self.name, self.sha256), bytes.clone());
        Ok(bytes)
    }

    pub fn text(&self) -> io::Result<Arc<str>> {
        let bytes = self.bytes()?;
        let text = std::str::from_utf8(&bytes)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        Ok(Arc::from(text))
    }

    /// Build-script entry point: emit expectations, never payloads, into Rust.
    /// Dependency-owned sources can provide their bytes through `external`.
    pub fn generate(
        root: &Path,
        output: &Path,
        external: &BTreeMap<String, &[u8]>,
    ) -> io::Result<()> {
        #[derive(Deserialize)]
        struct Definition {
            constant: String,
            name: String,
            path: String,
            description: String,
            package: Option<String>,
        }
        #[derive(Deserialize)]
        struct Manifest {
            schema: u32,
            assets: Vec<Definition>,
        }
        let manifest_path = root.join("data.json");
        println!("cargo:rerun-if-changed={}", manifest_path.display());
        let manifest: Manifest = serde_json::from_slice(&fs::read(manifest_path)?)?;
        if manifest.schema != 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unsupported asset manifest schema",
            ));
        }
        let mut generated = String::new();
        let mut constants = Vec::new();
        for entry in manifest.assets {
            let bytes = if entry.package.is_some() {
                external
                    .get(&entry.name)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::NotFound,
                            format!("missing build asset {}", entry.name),
                        )
                    })?
                    .to_vec()
            } else {
                let path = root.join(&entry.path);
                println!("cargo:rerun-if-changed={}", path.display());
                fs::read(path)?
            };
            let hash = format!("{:x}", Sha256::digest(&bytes));
            generated.push_str(&format!(
                "pub const {}: symbolica_data::Asset = symbolica_data::Asset {{ name: {:?}, sha256: {:?}, size: {}, description: {:?} }};\n",
                entry.constant, entry.name, hash, bytes.len(), entry.description
            ));
            constants.push(entry.constant);
        }
        generated.push_str(&format!(
            "pub const ALL: &[symbolica_data::Asset] = &[{}];\n",
            constants.join(",")
        ));
        fs::write(output, generated)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failed_load_and_corrupt_install_do_not_prevent_a_verified_retry() {
        let hash = format!("{:x}", Sha256::digest(b"verified payload"));
        let asset = Asset {
            name: "test/retry",
            sha256: Box::leak(hash.into_boxed_str()),
            size: 16,
            description: "retry fixture",
        };
        assert!(asset.bytes().is_err());
        assert!(asset.install(b"damaged payload!".to_vec()).is_err());
        assert!(!asset.is_loaded());
        asset.install(b"verified payload".to_vec()).unwrap();
        assert_eq!(&*asset.bytes().unwrap(), b"verified payload");
    }
}
