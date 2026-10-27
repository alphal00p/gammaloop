//! Explicit single-file admission for the measurement collectors.
use std::{
    ffi::{OsStr, OsString},
    io,
    path::Path,
};

use serde_json::{Value, json};
use vakint::{RustRedCacheFileOptions, RustRedCacheLoadOptions};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct CacheFileArgs {
    options: Option<RustRedCacheFileOptions>,
}

impl CacheFileArgs {
    pub(super) fn recognizes(flag: &str) -> bool {
        matches!(
            flag,
            "--cache-max-compressed-bytes"
                | "--cache-max-decoded-bytes"
                | "--cache-max-index-bytes"
                | "--cache-decoded-directory"
        )
    }

    fn byte_grant(flag: &str, value: &OsStr) -> io::Result<usize> {
        value
            .to_str()
            .and_then(|value| value.parse().ok())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("{flag} requires an unsigned byte count fitting this platform"),
                )
            })
    }

    pub(super) fn set(&mut self, flag: &str, value: OsString) -> io::Result<()> {
        match flag {
            "--cache-max-compressed-bytes" => {
                self.options
                    .get_or_insert_default()
                    .limits
                    .max_compressed_bytes = Self::byte_grant(flag, &value)?;
            }
            "--cache-max-decoded-bytes" => {
                let options = self.options.get_or_insert_default();
                options.limits = options
                    .limits
                    .with_decoded_bytes(Self::byte_grant(flag, &value)?)
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
            }
            "--cache-max-index-bytes" => {
                self.options
                    .get_or_insert_default()
                    .limits
                    .chunk
                    .max_index_owned_bytes = Self::byte_grant(flag, &value)?;
            }
            "--cache-decoded-directory" => {
                self.options.get_or_insert_default().decoded_directory = Some(value.into());
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown cache file option {flag}"),
                ));
            }
        }
        Ok(())
    }

    pub(super) fn normalize(&mut self, cache: Option<&Path>) -> io::Result<()> {
        let Some(options) = self.options.as_mut() else {
            return Ok(());
        };
        if cache.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cache file grants and --cache-decoded-directory require --cache",
            ));
        }
        if let Some(directory) = &options.decoded_directory {
            let resolved = directory.canonicalize()?;
            if !resolved.is_dir() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--cache-decoded-directory must name an existing directory",
                ));
            }
            options.decoded_directory = Some(resolved);
        }
        // Paths must be representable in the durable measurement identity.
        self.document()?;
        Ok(())
    }

    pub(super) fn load_options(&self) -> RustRedCacheLoadOptions {
        RustRedCacheLoadOptions {
            file_options: self.options.clone(),
            ..Default::default()
        }
    }

    pub(super) fn document(&self) -> io::Result<Value> {
        let options = self.options.clone().unwrap_or_default();
        let directory = serde_json::to_value(&options.decoded_directory)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        Ok(json!({
            "explicit": self.options.is_some(),
            "max_compressed_bytes": options.limits.max_compressed_bytes,
            "max_decoded_bytes": options.limits.max_decoded_bytes(),
            "max_index_bytes": options.limits.chunk.max_index_owned_bytes,
            "decoded_directory": directory,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{collections::BTreeMap, fs};

    #[test]
    fn omitted_grants_keep_public_compact_defaults_and_explicit_grants_need_cache() {
        let mut omitted = CacheFileArgs::default();
        omitted.normalize(None).unwrap();
        assert!(omitted.load_options().file_options.is_none());
        let defaults = RustRedCacheFileOptions::default();
        let document = omitted.document().unwrap();
        assert_eq!(document["explicit"], false);
        assert_eq!(
            document["max_compressed_bytes"],
            defaults.limits.max_compressed_bytes
        );
        assert_eq!(
            document["max_decoded_bytes"],
            defaults.limits.max_decoded_bytes()
        );
        for (flag, value) in [
            ("--cache-max-compressed-bytes", "0"),
            ("--cache-max-decoded-bytes", "1048576"),
            ("--cache-max-index-bytes", "0"),
            ("--cache-decoded-directory", "unopened-directory"),
        ] {
            let mut args = CacheFileArgs::default();
            args.set(flag, value.into()).unwrap();
            assert!(
                args.normalize(None)
                    .unwrap_err()
                    .to_string()
                    .contains("require --cache")
            );
        }
        assert!(!CacheFileArgs::recognizes("--cache-max-records"));
    }

    #[test]
    fn invalid_bytes_and_native_decoded_bounds_error_without_widening_zero_grants() {
        for value in ["-1", "1.5", "", "not-bytes", &format!("{}0", usize::MAX)] {
            assert!(
                CacheFileArgs::default()
                    .set("--cache-max-index-bytes", value.into())
                    .is_err()
            );
        }
        assert!(
            CacheFileArgs::default()
                .set("--cache-max-decoded-bytes", "0".into())
                .is_err()
        );
        let mut args = CacheFileArgs::default();
        args.set("--cache-max-compressed-bytes", "0".into())
            .unwrap();
        args.set("--cache-max-index-bytes", "0".into()).unwrap();
        args.normalize(Some(Path::new("cache.gz"))).unwrap();
        let options = args.load_options().file_options.unwrap();
        assert_eq!(options.limits.max_compressed_bytes, 0);
        assert_eq!(options.limits.chunk.max_index_owned_bytes, 0);
        assert!(!options.limits.chunk.diagnostics_enabled);
    }

    #[test]
    fn decoded_directory_is_absolute_and_requires_an_existing_directory() -> io::Result<()> {
        let scratch = tempfile::tempdir()?;
        fs::create_dir(scratch.path().join("decoded"))?;
        let mut args = CacheFileArgs::default();
        args.set(
            "--cache-decoded-directory",
            scratch.path().join("decoded/../decoded").into_os_string(),
        )?;
        args.normalize(Some(Path::new("cache.gz")))?;
        assert_eq!(
            args.load_options().file_options.unwrap().decoded_directory,
            Some(scratch.path().join("decoded").canonicalize()?)
        );
        fs::write(scratch.path().join("regular-file"), [])?;
        for invalid in [
            scratch.path().join("missing"),
            scratch.path().join("regular-file"),
        ] {
            let mut args = CacheFileArgs::default();
            args.set("--cache-decoded-directory", invalid.into_os_string())?;
            assert!(args.normalize(Some(Path::new("cache.gz"))).is_err());
        }
        Ok(())
    }

    #[test]
    fn real_native_file_checks_each_byte_grant_and_shares_one_lazy_import()
    -> Result<(), Box<dyn std::error::Error>> {
        use rustred::{
            algebra::CoefficientContext,
            family::IntegralKey,
            persistence::application_cache::{
                ApplicationCacheBinding, ApplicationCacheKey, ApplicationCacheRecord,
            },
        };
        use rustred_app::{AppErrorKind, SavedApplicationCacheFile};

        let scratch = tempfile::tempdir()?;
        let decoded_directory = scratch.path().join("decoded");
        fs::create_dir(&decoded_directory)?;
        let context = CoefficientContext::try_new(["collector_cache_file_control_d"])?;
        let d = context.parameter("collector_cache_file_control_d").unwrap();
        let guard = context.try_sub(&d, &context.integer(4), Default::default())?;
        let record = ApplicationCacheRecord {
            key: ApplicationCacheKey::route(&IntegralKey::try_new([2])?),
            terms: BTreeMap::from([(
                IntegralKey::try_new([1])?,
                context.try_div(&context.one(), &guard, Default::default())?,
            )]),
            nonzero_conditions: vec![guard],
        };
        let binding = ApplicationCacheBinding {
            artifact_digest: [77; 32],
            family_fingerprint: "collector-file-options-control".into(),
            physical_arity: 1,
            layer: "routed-before-terminal-substitution-v1".into(),
        };
        let (encoded, identity) = SavedApplicationCacheFile::encode_generated(
            &binding,
            std::slice::from_ref(&record),
            &context,
            Default::default(),
        )?;
        let cache = scratch.path().join("small.rrcache.gz");
        fs::write(&cache, encoded)?;
        for (flag, value) in [
            (
                "--cache-max-compressed-bytes",
                identity.compressed_bytes - 1,
            ),
            ("--cache-max-decoded-bytes", identity.decoded_bytes - 1),
            ("--cache-max-index-bytes", 0),
        ] {
            let mut args = CacheFileArgs::default();
            args.set(flag, value.to_string().into())?;
            args.set(
                "--cache-decoded-directory",
                decoded_directory.clone().into_os_string(),
            )?;
            args.normalize(Some(&cache))?;
            let error = SavedApplicationCacheFile::open_with_options(
                &cache,
                &binding,
                Some(&identity),
                &args.load_options().file_options.unwrap(),
            )
            .err()
            .unwrap();
            assert_eq!(error.kind(), AppErrorKind::Limit, "{flag}: {error}");
        }
        let mut args = CacheFileArgs::default();
        args.set(
            "--cache-decoded-directory",
            decoded_directory.into_os_string(),
        )?;
        args.normalize(Some(&cache))?;
        let mut reader = SavedApplicationCacheFile::open_with_options(
            &cache,
            &binding,
            Some(&identity),
            &args.load_options().file_options.unwrap(),
        )?;
        assert_eq!(reader.statistics()?.state_imports, 0);
        assert_eq!(reader.lookup(&record.key, &context)?, Some(record.clone()));
        let mut fork = reader.fork_reader();
        assert_eq!(fork.lookup(&record.key, &context)?, Some(record));
        assert_eq!(reader.statistics()?.state_imports, 1);
        assert_eq!(fork.statistics()?.state_imports, 1);
        Ok(())
    }
}
