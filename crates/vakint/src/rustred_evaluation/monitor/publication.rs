//! The existing atomic snapshot sink, serviced by one bounded latest-value writer.
//! Blocking file/stream I/O never owns the terminal or the native observer lock.

use std::{
    collections::HashMap,
    io,
    path::{Path, PathBuf},
    sync::{
        Arc, Condvar, Mutex, OnceLock,
        atomic::{AtomicBool, AtomicU64, Ordering},
    },
    thread,
    time::{Duration, Instant},
};

use super::{LatestJsonFile, RustRedMonitorOptions, SnapshotWriter};

#[derive(Default)]
struct Destination {
    generation: AtomicU64,
    commit: Mutex<()>,
}

pub(super) struct FileLease {
    destination: Arc<Destination>,
    generation: u64,
}

impl FileLease {
    pub(super) fn next_owner() -> u64 {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        NEXT.fetch_add(1, Ordering::Relaxed) + 1
    }

    pub(super) fn with_owner(path: &Path, owner: u64) -> Self {
        static DESTINATIONS: OnceLock<Mutex<HashMap<PathBuf, Arc<Destination>>>> = OnceLock::new();
        // Resolve an existing parent on the I/O owner, never on the UI loop.
        let parent = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let key = match (parent.canonicalize(), path.file_name()) {
            (Ok(parent), Some(name)) => parent.join(name),
            _ if path.is_absolute() => path.to_owned(),
            _ => std::env::current_dir().unwrap_or_default().join(path),
        };
        let owners = DESTINATIONS.get_or_init(Default::default);
        let mut owners = owners
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        // Keep a tiny generation tombstone per requested destination even when
        // its current writer has finished. An older first registration may
        // still be delayed in I/O; weak ownership would forget the high-water
        // token and wrongly let that writer become current again. These are
        // destination identities only, never snapshots or native cache state.
        let destination = Arc::clone(owners.entry(key).or_default());
        // Registration may happen much later on a blocked I/O thread. Keep
        // request creation order, never late writer arrival order.
        destination.generation.fetch_max(owner, Ordering::AcqRel);
        let generation = owner;
        Self {
            destination,
            generation,
        }
    }
    pub(super) fn commit(&self) -> std::sync::MutexGuard<'_, ()> {
        self.destination
            .commit
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }
    pub(super) fn is_current(&self) -> bool {
        self.destination.generation.load(Ordering::Acquire) == self.generation
    }
}

#[derive(Default)]
struct Queue {
    pending: Option<serde_json::Value>,
    closing: bool,
    done: bool,
    error: Option<String>,
}

#[derive(Default)]
struct Shared {
    queue: Mutex<Queue>,
    ready: Condvar,
    active: AtomicBool,
}

/// One pending document, one current writer, no unbounded frame backlog.
pub(super) struct SnapshotPublisher {
    shared: Arc<Shared>,
}

impl SnapshotPublisher {
    pub(super) fn new(
        options: &RustRedMonitorOptions,
        stderr_json: bool,
        owner: u64,
    ) -> io::Result<Self> {
        let options = options.clone();
        let mut file = None;
        let mut stream = SnapshotWriter::new(io::stderr());
        Self::spawn(move |document, active| {
            let file = file.get_or_insert_with(|| LatestJsonFile::from_owner(&options, owner));
            file.write_if(document, active)?;
            if stderr_json && active.load(Ordering::Acquire) {
                stream.write(document)?;
            }
            Ok(())
        })
    }

    pub(super) fn spawn(
        mut write: impl FnMut(&serde_json::Value, &AtomicBool) -> io::Result<()> + Send + 'static,
    ) -> io::Result<Self> {
        let shared = Arc::new(Shared::default());
        shared.active.store(true, Ordering::Release);
        let worker = Arc::clone(&shared);
        thread::Builder::new()
            .name("vakint-monitor-publisher".into())
            .spawn(move || {
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    loop {
                        let document = {
                            let mut queue = worker
                                .queue
                                .lock()
                                .unwrap_or_else(std::sync::PoisonError::into_inner);
                            while queue.pending.is_none() && !queue.closing {
                                queue = worker
                                    .ready
                                    .wait(queue)
                                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                            }
                            match queue.pending.take() {
                                Some(document) => document,
                                None => return Ok(()),
                            }
                        };
                        write(&document, &worker.active)?;
                    }
                }))
                .unwrap_or_else(|_| Err(io::Error::other("monitor publisher panicked")));
                let mut queue = worker
                    .queue
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                queue.error = result.err().map(|error| error.to_string());
                queue.done = true;
                worker.ready.notify_all();
            })?;
        Ok(Self { shared })
    }

    pub(super) fn publish(&self, document: serde_json::Value) {
        let mut queue = self
            .shared
            .queue
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if !queue.closing && !queue.done {
            queue.pending = Some(document);
            self.shared.ready.notify_one();
        }
    }

    pub(super) fn error(&self) -> Option<String> {
        self.shared
            .queue
            .try_lock()
            .ok()
            .and_then(|queue| queue.error.clone())
    }

    fn close(&self, allowance: Duration) -> bool {
        let deadline = Instant::now() + allowance;
        let mut queue = self
            .shared
            .queue
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        queue.closing = true;
        self.shared.ready.notify_one();
        while !queue.done {
            let remaining = deadline.saturating_duration_since(Instant::now());
            if remaining.is_zero() {
                self.shared.active.store(false, Ordering::Release);
                queue.pending = None;
                return false;
            }
            let (next, _) = self
                .shared
                .ready
                .wait_timeout(queue, remaining)
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            queue = next;
        }
        true
    }

    pub(super) fn finish(&self) {
        // A diagnostic sink receives a short drain allowance, never an
        // unbounded join. Its delayed atomic commit is revoked on timeout.
        if !self.close(Duration::from_millis(100)) {
            log::warn!(
                "RustRed final monitor publication did not drain; pending telemetry cancelled, scientific result unchanged"
            );
        }
        if let Some(error) = self.error() {
            log::warn!("RustRed monitor publication failed: {error}");
        }
    }
}

impl Drop for SnapshotPublisher {
    fn drop(&mut self) {
        let _ = self.close(Duration::ZERO);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;

    #[test]
    fn stalled_writer_is_revoked_before_atomic_commit_and_new_owner_survives() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let options = RustRedMonitorOptions {
            progress_snapshot: Some(path.clone()),
            ..Default::default()
        };
        let alias = directory.path().join("alias");
        #[cfg(unix)]
        std::os::unix::fs::symlink(directory.path(), &alias).unwrap();
        #[cfg(unix)]
        let old_options = RustRedMonitorOptions {
            progress_snapshot: Some(alias.join("latest.json")),
            ..options.clone()
        };
        #[cfg(not(unix))]
        let old_options = options.clone();
        let old = LatestJsonFile::from_options(&old_options);
        let (entered, ready) = mpsc::channel();
        let (release, resume) = mpsc::channel();
        let publisher = SnapshotPublisher::spawn(move |document, active| {
            old.write_with_check(
                |file| {
                    serde_json::to_writer(file, document)?;
                    entered.send(()).unwrap();
                    resume.recv().unwrap();
                    Ok(())
                },
                active,
            )
        })
        .unwrap();
        publisher.publish(serde_json::json!({"old": 1}));
        ready.recv_timeout(Duration::from_secs(2)).unwrap();
        let started = Instant::now();
        assert!(!publisher.close(Duration::from_millis(100)));
        assert!(started.elapsed() < Duration::from_secs(1));
        let current = LatestJsonFile::from_options(&options);
        current.write(&serde_json::json!({"new": 2})).unwrap();
        release.send(()).unwrap();
        assert!(publisher.close(Duration::from_secs(2)));
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&std::fs::read(path).unwrap()).unwrap(),
            serde_json::json!({"new": 2})
        );
    }

    #[cfg(unix)]
    #[test]
    fn active_superseded_alias_writer_cannot_overwrite_a_newer_owner() {
        let directory = tempfile::tempdir().unwrap();
        let alias = directory.path().join("alias");
        std::os::unix::fs::symlink(directory.path(), &alias).unwrap();
        let options = RustRedMonitorOptions {
            progress_snapshot: Some(directory.path().join("latest.json")),
            ..Default::default()
        };
        let old = LatestJsonFile::from_options(&RustRedMonitorOptions {
            progress_snapshot: Some(alias.join("latest.json")),
            ..options.clone()
        });
        let (entered, ready) = mpsc::channel();
        let (release, resume) = mpsc::channel();
        let writer = thread::spawn(move || {
            old.write_with_check(
                |file| {
                    use std::io::Write;
                    file.write_all(b"old")?;
                    entered.send(()).unwrap();
                    resume.recv().unwrap();
                    Ok(())
                },
                &AtomicBool::new(true),
            )
        });
        ready.recv_timeout(Duration::from_secs(2)).unwrap();
        LatestJsonFile::from_options(&options)
            .write(&serde_json::json!({"current": true}))
            .unwrap();
        release.send(()).unwrap();
        writer.join().unwrap().unwrap();
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(
                &std::fs::read(options.progress_snapshot.unwrap()).unwrap()
            )
            .unwrap(),
            serde_json::json!({"current": true})
        );
    }

    #[test]
    fn late_first_registration_cannot_supersede_a_newer_request() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("latest.json");
        let options = RustRedMonitorOptions {
            progress_snapshot: Some(path.clone()),
            ..Default::default()
        };
        let old_owner = FileLease::next_owner();
        let new_owner = FileLease::next_owner();
        LatestJsonFile::from_owner(&options, new_owner)
            .write(&serde_json::json!({"new": true}))
            .unwrap();
        LatestJsonFile::from_owner(&options, old_owner)
            .write(&serde_json::json!({"old": true}))
            .unwrap();
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(&std::fs::read(path).unwrap()).unwrap(),
            serde_json::json!({"new": true})
        );
    }

    #[test]
    fn slow_sink_coalesces_one_pending_document_and_retains_final_error() {
        let (entered, ready) = mpsc::channel();
        let (release, resume) = mpsc::channel();
        let values = Arc::new(Mutex::new(Vec::new()));
        let written = Arc::clone(&values);
        let publisher = SnapshotPublisher::spawn(move |document, _active| {
            if document == &serde_json::json!(0) {
                entered.send(()).unwrap();
                resume.recv().unwrap();
            }
            written.lock().unwrap().push(document.clone());
            if document == &serde_json::json!(99) {
                return Err(io::Error::other("synthetic sink failure"));
            }
            Ok(())
        })
        .unwrap();
        publisher.publish(serde_json::json!(0));
        ready.recv_timeout(Duration::from_secs(2)).unwrap();
        for value in 1..100 {
            publisher.publish(serde_json::json!(value));
        }
        release.send(()).unwrap();
        assert!(publisher.close(Duration::from_secs(2)));
        assert_eq!(
            *values.lock().unwrap(),
            vec![serde_json::json!(0), serde_json::json!(99)]
        );
        assert_eq!(publisher.error().as_deref(), Some("synthetic sink failure"));
    }
}
