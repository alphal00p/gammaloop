use std::fs;
use std::path::PathBuf;

use tempfile::TempDir;

pub(crate) struct TestDirectory(pub(crate) PathBuf, Option<TempDir>);

impl TestDirectory {
    pub(crate) fn new() -> Self {
        let directory = tempfile::tempdir().unwrap();
        Self(fs::canonicalize(directory.path()).unwrap(), Some(directory))
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        self.1.take().unwrap().close().unwrap();
    }
}
