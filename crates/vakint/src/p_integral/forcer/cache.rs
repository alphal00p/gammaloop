//! Exact raw FORM outputs, before coefficients, epsilon symbols or master
//! expansions are attached. Embedded Forcer sources are fixed for this
//! process; executable identity and the complete program distinguish runs.

use std::{
    collections::HashMap,
    ffi::OsString,
    path::{Path, PathBuf},
    process::Command,
    sync::{Arc, Mutex},
    time::SystemTime,
};

use once_cell::sync::OnceCell;

use crate::VakintError;

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(super) struct ForcerRunKey {
    executable: PathBuf,
    modified: SystemTime,
    executable_bytes: u64,
    #[cfg(unix)]
    file_identity: (u64, u64, i64, i64),
    arguments: Vec<OsString>,
    program: String,
    locals: Vec<String>,
}

impl ForcerRunKey {
    pub(super) fn new(command: &Command, program: &str, locals: &[&str]) -> Option<Self> {
        let executable = Path::new(command.get_program());
        let executable = if executable.components().count() == 1 {
            std::env::split_paths(&std::env::var_os("PATH")?)
                .map(|directory| directory.join(executable))
                .find(|path| {
                    std::fs::metadata(path).is_ok_and(|metadata| {
                        #[cfg(unix)]
                        {
                            use std::os::unix::fs::PermissionsExt;
                            metadata.is_file() && metadata.permissions().mode() & 0o111 != 0
                        }
                        #[cfg(not(unix))]
                        metadata.is_file()
                    })
                })?
        } else {
            executable.to_owned()
        };
        let executable = executable.canonicalize().ok()?;
        let metadata = std::fs::metadata(&executable).ok()?;
        Some(Self {
            executable,
            modified: metadata.modified().ok()?,
            executable_bytes: metadata.len(),
            #[cfg(unix)]
            file_identity: {
                use std::os::unix::fs::MetadataExt;
                (
                    metadata.dev(),
                    metadata.ino(),
                    metadata.ctime(),
                    metadata.ctime_nsec(),
                )
            },
            arguments: command.get_args().map(OsString::from).collect(),
            program: program.to_owned(),
            locals: locals.iter().map(|local| (*local).to_owned()).collect(),
        })
    }

    fn bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.executable.as_os_str().len()
            + self
                .arguments
                .iter()
                .map(|argument| argument.len())
                .sum::<usize>()
            + self.program.len()
            + self.locals.iter().map(String::len).sum::<usize>()
    }
}

type Output = Arc<Vec<String>>;

struct Entry {
    output: Arc<OnceCell<Output>>,
    bytes: usize,
    used: u64,
}

#[derive(Default)]
struct State {
    entries: HashMap<ForcerRunKey, Entry>,
    bytes: usize,
    clock: u64,
    hits: u64,
    misses: u64,
    evictions: u64,
}

impl State {
    fn remove(&mut self, key: &ForcerRunKey) {
        if let Some(entry) = self.entries.remove(key) {
            self.bytes -= entry.bytes;
            self.evictions += 1;
        }
    }

    fn make_room(&mut self, bytes: usize, max_bytes: usize, max_entries: usize) -> bool {
        while self.bytes.saturating_add(bytes) > max_bytes || self.entries.len() >= max_entries {
            // Pending runs keep their cells so concurrent callers share them.
            let oldest = self
                .entries
                .iter()
                .filter(|(_, entry)| entry.output.get().is_some())
                .min_by_key(|(_, entry)| entry.used)
                .map(|(key, _)| key.clone());
            let Some(oldest) = oldest else { return false };
            self.remove(&oldest);
        }
        true
    }
}

pub(super) struct ForcerRunCache {
    state: Mutex<State>,
    max_bytes: usize,
    max_entries: usize,
}

impl Default for ForcerRunCache {
    fn default() -> Self {
        Self::new(64 << 20, 4096)
    }
}

impl ForcerRunCache {
    fn new(max_bytes: usize, max_entries: usize) -> Self {
        Self {
            state: Mutex::new(State::default()),
            max_bytes,
            max_entries,
        }
    }

    pub(super) fn run(
        &self,
        key: Option<ForcerRunKey>,
        compute: impl FnOnce() -> Result<Vec<String>, VakintError>,
    ) -> Result<Vec<String>, VakintError> {
        let Some(key) = key else { return compute() };
        let cell = {
            let mut state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            state.clock += 1;
            let used = state.clock;
            if let Some(entry) = state.entries.get_mut(&key) {
                entry.used = used;
                let cell = entry.output.clone();
                state.hits += 1;
                cell
            } else {
                state.misses += 1;
                let bytes = key.bytes() + std::mem::size_of::<Entry>();
                if !state.make_room(bytes, self.max_bytes, self.max_entries) {
                    drop(state);
                    return compute();
                }
                let cell = Arc::new(OnceCell::new());
                state.entries.insert(
                    key.clone(),
                    Entry {
                        output: cell.clone(),
                        bytes,
                        used,
                    },
                );
                state.bytes += bytes;
                cell
            }
        };
        // OnceCell shares successful concurrent requests. Its fallible
        // initializer leaves errors uncached; the registry lock is released.
        let result = cell.get_or_try_init(|| {
            let output = Arc::new(compute()?);
            let bytes = output
                .iter()
                .map(|value| value.capacity() + std::mem::size_of::<String>())
                .sum();
            let mut state = self
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            if state
                .entries
                .get(&key)
                .is_some_and(|entry| Arc::ptr_eq(&entry.output, &cell))
            {
                if state.make_room(bytes, self.max_bytes, usize::MAX) {
                    state.entries.get_mut(&key).unwrap().bytes += bytes;
                    state.bytes += bytes;
                } else {
                    state.remove(&key);
                }
            }
            Ok::<_, VakintError>(output)
        });
        let mut state = self
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if result.is_err()
            && state
                .entries
                .get(&key)
                .is_some_and(|entry| Arc::ptr_eq(&entry.output, &cell))
        {
            state.remove(&key);
        }
        log::debug!(target: "vakint::p_integral::forcer_cache",
            "raw Forcer cache hits={} misses={} evictions={} entries={} resident_bytes={}",
            state.hits, state.misses, state.evictions, state.entries.len(), state.bytes);
        drop(state);
        result.map(|output| output.as_ref().clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{
        Barrier,
        atomic::{AtomicUsize, Ordering},
    };

    fn key(program: &str) -> ForcerRunKey {
        ForcerRunKey::new(
            &Command::new(std::env::current_exe().unwrap()),
            program,
            &["F1"],
        )
        .unwrap()
    }

    #[test]
    fn raw_cache_shares_concurrent_identical_runs() {
        let cache = ForcerRunCache::default();
        let calls = AtomicUsize::new(0);
        let barrier = Barrier::new(8);
        std::thread::scope(|scope| {
            for _ in 0..8 {
                scope.spawn(|| {
                    barrier.wait();
                    let output = cache
                        .run(Some(key("same")), || {
                            calls.fetch_add(1, Ordering::SeqCst);
                            Ok(vec!["rat(1,ep)".into()])
                        })
                        .unwrap();
                    assert_eq!(output, ["rat(1,ep)"]);
                });
            }
        });
        assert_eq!(calls.load(Ordering::SeqCst), 1);
        let state = cache.state.lock().unwrap();
        assert_eq!((state.hits, state.misses), (7, 1));
    }

    #[test]
    fn raw_cache_keeps_different_runs_parallel_and_retries_errors() {
        let cache = ForcerRunCache::default();
        let barrier = Barrier::new(2);
        std::thread::scope(|scope| {
            for program in ["left", "right"] {
                let cache = &cache;
                let barrier = &barrier;
                scope.spawn(move || {
                    cache
                        .run(Some(key(program)), || {
                            barrier.wait();
                            Ok(vec![program.into()])
                        })
                        .unwrap();
                });
            }
        });
        assert!(
            cache
                .run(Some(key("error")), || Err(
                    VakintError::InvalidIntegralFormat("failed".into())
                ))
                .is_err()
        );
        assert_eq!(
            cache
                .run(Some(key("error")), || Ok(vec!["success".into()]))
                .unwrap(),
            ["success"]
        );
    }

    #[test]
    fn raw_cache_bounds_retained_outputs_and_distinguishes_run_inputs() {
        let cache = ForcerRunCache::new(4096, 2);
        for program in ["first", "second", "third", "fourth"] {
            cache
                .run(Some(key(program)), || Ok(vec!["x".repeat(512)]))
                .unwrap();
            let state = cache.state.lock().unwrap();
            assert!(state.bytes <= 4096);
            assert!(state.entries.len() <= 2);
        }
        cache
            .run(Some(key("oversized")), || Ok(vec!["x".repeat(8192)]))
            .unwrap();
        assert!(
            !cache
                .state
                .lock()
                .unwrap()
                .entries
                .contains_key(&key("oversized"))
        );
        let mut other = key("same");
        other.arguments.push("-w2".into());
        assert_ne!(key("same"), other);
        other = key("same");
        other.locals.push("F2".into());
        assert_ne!(key("same"), other);
        assert_ne!(key("same"), key("different"));
    }

    #[test]
    fn raw_cache_key_tracks_executable_replacement() {
        let directory = tempfile::tempdir().unwrap();
        let executable = directory.path().join("form");
        std::fs::write(&executable, "first executable").unwrap();
        let command = Command::new(&executable);
        let first = ForcerRunKey::new(&command, "program", &["F1"]).unwrap();
        std::fs::write(&executable, "replacement executable").unwrap();
        assert_ne!(
            first,
            ForcerRunKey::new(&command, "program", &["F1"]).unwrap()
        );
    }
}
