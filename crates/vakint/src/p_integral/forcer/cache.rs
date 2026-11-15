//! Exact raw FORM outputs, before coefficients, epsilon symbols or master
//! expansions are attached. Embedded Forcer sources are fixed for this
//! process; executable identity and the complete program distinguish runs.

use std::{
    collections::HashMap,
    ffi::OsString,
    path::{Path, PathBuf},
    process::Command,
    sync::{Arc, Condvar, Mutex},
    time::SystemTime,
};

use once_cell::sync::OnceCell;

use crate::{VakintError, stage_timing::StageSpan};

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

#[derive(Default)]
struct Cell {
    output: OnceCell<Output>,
    running: Mutex<bool>,
    ready: Condvar,
}

impl Cell {
    fn claim(self: &Arc<Self>) -> Option<Reservation> {
        let mut running = self
            .running
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if *running || self.output.get().is_some() {
            return None;
        }
        *running = true;
        Some(Reservation(self.clone()))
    }

    fn wait(&self) {
        let mut running = self
            .running
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        while *running && self.output.get().is_none() {
            // Includes the condvar's mutex reacquisition, but no ready hits or
            // work performed by this caller to produce its own requests.
            let _stage = StageSpan::new("p_integral.forcer.cache.pending_wait");
            running = self
                .ready
                .wait(running)
                .unwrap_or_else(std::sync::PoisonError::into_inner);
        }
    }
}

// Failed or panicking owners release their reservations, allowing a waiting
// caller to retry. A successful value is installed before waking its peers.
struct Reservation(Arc<Cell>);

impl Drop for Reservation {
    fn drop(&mut self) {
        *self
            .0
            .running
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = false;
        self.0.ready.notify_all();
    }
}

struct Entry {
    output: Arc<Cell>,
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
                .filter(|(_, entry)| entry.output.output.get().is_some())
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

/// Own only unfinished requests. On failure or unwind, completed independent
/// results remain reusable and peers can retry the released reservations.
struct PendingBatch<'a> {
    cache: &'a ForcerRunCache,
    keys: &'a [Option<ForcerRunKey>],
    reservations: Vec<Option<Reservation>>,
    invalid_publication: bool,
}

impl PendingBatch<'_> {
    fn publish(&mut self, index: usize, output: Vec<String>) -> Result<(), VakintError> {
        let _stage = StageSpan::new("p_integral.forcer.cache.publish");
        let Some(Some(reservation)) = self.reservations.get(index) else {
            self.invalid_publication = true;
            return Err(VakintError::InvalidGenericExpression(format!(
                "FORCER batch published unowned or already published request {index}"
            )));
        };
        let output = Arc::new(output);
        let bytes = output
            .iter()
            .map(|value| value.capacity() + std::mem::size_of::<String>())
            .sum();
        if let Some(key) = &self.keys[index] {
            let mut state = self
                .cache
                .state
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            if state
                .entries
                .get(key)
                .is_some_and(|entry| Arc::ptr_eq(&entry.output, &reservation.0))
            {
                if state.make_room(bytes, self.cache.max_bytes, usize::MAX) {
                    state.entries.get_mut(key).unwrap().bytes += bytes;
                    state.bytes += bytes;
                } else {
                    state.remove(key);
                }
            }
        }
        // Publish even if retention failed: callers already holding this cell
        // must receive its value, independently of later map replacements.
        reservation
            .0
            .output
            .set(output)
            .expect("exclusive FORCER reservation");
        drop(self.reservations[index].take());
        Ok(())
    }
}

impl Drop for PendingBatch<'_> {
    fn drop(&mut self) {
        if self.reservations.iter().all(Option::is_none) {
            return;
        }
        let mut state = self
            .cache
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        for (key, reservation) in self.keys.iter().zip(&self.reservations) {
            if let (Some(key), Some(reservation)) = (key, reservation)
                && state
                    .entries
                    .get(key)
                    .is_some_and(|entry| Arc::ptr_eq(&entry.output, &reservation.0))
            {
                state.remove(key);
            }
        }
        // The state lock is released before the fields drop and wake peers.
    }
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
        let mut compute = Some(compute);
        let mut values = self.run_batch(&[key], |_, publish| {
            publish(0, compute.take().expect("one reservation per request")()?)
        })?;
        Ok(values.remove(0))
    }

    /// Reserve each logical program independently, execute this caller's
    /// misses together, publishing each completed group immediately. Finish
    /// owned work before waiting for peers: reversed overlapping batches must
    /// never wait for each other while holding unpublished owned requests.
    pub(super) fn run_batch(
        &self,
        keys: &[Option<ForcerRunKey>],
        mut compute: impl FnMut(
            &[usize],
            &mut dyn FnMut(usize, Vec<String>) -> Result<(), VakintError>,
        ) -> Result<(), VakintError>,
    ) -> Result<Vec<Vec<String>>, VakintError> {
        let lookup_stage = StageSpan::new("p_integral.forcer.cache.map_lookup");
        let cells = keys
            .iter()
            .map(|key| {
                let Some(key) = key else {
                    return Arc::new(Cell::default());
                };
                let mut state = self
                    .state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                state.clock += 1;
                let used = state.clock;
                if let Some(entry) = state.entries.get_mut(key) {
                    entry.used = used;
                    let cell = entry.output.clone();
                    state.hits += 1;
                    cell
                } else {
                    state.misses += 1;
                    let bytes = key.bytes() + std::mem::size_of::<Entry>();
                    if !state.make_room(bytes, self.max_bytes, self.max_entries) {
                        return Arc::new(Cell::default());
                    }
                    let cell = Arc::new(Cell::default());
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
            })
            .collect::<Vec<_>>();
        drop(lookup_stage);
        while cells.iter().any(|cell| cell.output.get().is_none()) {
            let mut pending = PendingBatch {
                cache: self,
                keys,
                reservations: cells.iter().map(Cell::claim).collect(),
                invalid_publication: false,
            };
            let indices = pending
                .reservations
                .iter()
                .enumerate()
                .filter_map(|(index, reservation)| reservation.as_ref().map(|_| index))
                .collect::<Vec<_>>();
            if !indices.is_empty() {
                compute(&indices, &mut |index, output| {
                    pending.publish(index, output)
                })?;
                if pending.invalid_publication || pending.reservations.iter().any(Option::is_some) {
                    return Err(VakintError::InvalidGenericExpression(
                        "FORCER batch must publish each owned request exactly once".into(),
                    ));
                }
            }
            drop(pending);
            if let Some(cell) = cells.iter().find(|cell| cell.output.get().is_none()) {
                cell.wait();
            }
        }
        if log::log_enabled!(target: "vakint::p_integral::forcer_cache", log::Level::Debug) {
            let (hits, misses, evictions, entries, bytes) = {
                let state = self
                    .state
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                (
                    state.hits,
                    state.misses,
                    state.evictions,
                    state.entries.len(),
                    state.bytes,
                )
            };
            log::debug!(target: "vakint::p_integral::forcer_cache",
                "raw Forcer cache hits={hits} misses={misses} evictions={evictions} entries={entries} resident_bytes={bytes}");
        }
        let _stage = StageSpan::new("p_integral.forcer.cache.output_clone");
        Ok(cells
            .iter()
            .map(|cell| cell.output.get().unwrap().as_ref().clone())
            .collect())
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

    #[test]
    fn batches_preserve_request_hits_across_order_and_mixed_misses() {
        let cache = ForcerRunCache::default();
        cache.run(Some(key("a")), || Ok(vec!["a".into()])).unwrap();
        let programs = ["b", "a", "c", "b"];
        let keys = programs.map(|program| Some(key(program)));
        let values = cache
            .run_batch(&keys, |missing, publish| {
                assert_eq!(missing, [0, 2]);
                for &index in missing.iter().rev() {
                    publish(index, vec![programs[index].into()])?;
                }
                Ok(())
            })
            .unwrap();
        assert_eq!(values, programs.map(|program| vec![program.to_string()]));
        let values = cache
            .run_batch(&[Some(key("c")), Some(key("b")), Some(key("a"))], |_, _| {
                panic!("reordering a batch must preserve every request cache hit")
            })
            .unwrap();
        assert_eq!(values, [vec!["c"], vec!["b"], vec!["a"]]);
        let state = cache.state.lock().unwrap();
        assert_eq!((state.hits, state.misses), (5, 3));
    }

    #[test]
    fn batch_publishes_owned_misses_before_waiting_for_peer() {
        let cache = ForcerRunCache::default();
        let (peer_started_tx, peer_started_rx) = std::sync::mpsc::channel();
        let (release_peer_tx, release_peer_rx) = std::sync::mpsc::channel();
        let (own_started_tx, own_started_rx) = std::sync::mpsc::channel();
        std::thread::scope(|scope| {
            let cache = &cache;
            scope.spawn(move || {
                cache
                    .run(Some(key("peer")), || {
                        peer_started_tx.send(()).unwrap();
                        release_peer_rx.recv().unwrap();
                        Ok(vec!["peer".into()])
                    })
                    .unwrap()
            });
            peer_started_rx.recv().unwrap();
            let worker = scope.spawn(move || {
                cache
                    .run_batch(
                        &[Some(key("peer")), Some(key("own"))],
                        |missing, publish| {
                            assert_eq!(missing, [1]);
                            publish(1, vec!["own".into()])?;
                            own_started_tx.send(()).unwrap();
                            Ok(())
                        },
                    )
                    .unwrap()
            });
            // The peer remains pending: batching must still run and publish
            // this caller's miss before waiting for it.
            let published = own_started_rx.recv_timeout(std::time::Duration::from_secs(5));
            let output = published
                .is_ok()
                .then(|| cache.run(Some(key("own")), || panic!("owned miss was not shared")));
            release_peer_tx.send(()).unwrap();
            published.unwrap();
            assert_eq!(output.unwrap().unwrap(), ["own"]);
            assert_eq!(worker.join().unwrap(), [vec!["peer"], vec!["own"]]);
        });
    }

    #[test]
    fn failed_batch_releases_all_request_reservations() {
        let cache = ForcerRunCache::default();
        let keys = [Some(key("left")), Some(key("right"))];
        assert!(
            cache
                .run_batch(&keys, |_, _| Err(VakintError::InvalidIntegralFormat(
                    "failed batch".into()
                )))
                .is_err()
        );
        let values = cache
            .run_batch(&keys, |missing, publish| {
                assert_eq!(missing, [0, 1]);
                publish(0, vec!["left".into()])?;
                publish(1, vec!["right".into()])
            })
            .unwrap();
        assert_eq!(values, [vec!["left"], vec!["right"]]);
    }

    #[test]
    fn batch_publishes_early_group_while_later_group_is_blocked() {
        let cache = ForcerRunCache::default();
        let (published_tx, published_rx) = std::sync::mpsc::channel();
        let (release_tx, release_rx) = std::sync::mpsc::channel();
        let (peer_done_tx, peer_done_rx) = std::sync::mpsc::channel();
        std::thread::scope(|scope| {
            let cache = &cache;
            let owner = scope.spawn(move || {
                cache.run_batch(&[Some(key("early")), Some(key("late"))], |_, publish| {
                    publish(0, vec!["early".into()])?;
                    published_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                    publish(1, vec!["late".into()])
                })
            });
            published_rx.recv().unwrap();
            let peer = scope.spawn(move || {
                let output = cache.run(Some(key("early")), || panic!("early result recomputed"));
                peer_done_tx.send(output).unwrap();
            });
            let early = peer_done_rx.recv_timeout(std::time::Duration::from_secs(5));
            // Always unblock the owner before asserting, including on timeout.
            release_tx.send(()).unwrap();
            assert_eq!(
                owner.join().unwrap().unwrap(),
                [vec!["early"], vec!["late"]]
            );
            peer.join().unwrap();
            assert_eq!(early.unwrap().unwrap(), ["early"]);
        });
    }

    #[test]
    fn partial_error_and_panic_keep_successes_and_release_unpublished_requests() {
        for panics in [false, true] {
            let cache = ForcerRunCache::default();
            let keys = [Some(key("complete")), Some(key("unfinished"))];
            let mut unfinished = None;
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                cache.run_batch(&keys, |_, publish| {
                    unfinished = Some(
                        cache.state.lock().unwrap().entries[keys[1].as_ref().unwrap()]
                            .output
                            .clone(),
                    );
                    publish(0, vec!["complete".into()])?;
                    if panics {
                        panic!("later FORM group failed");
                    }
                    Err(VakintError::InvalidIntegralFormat(
                        "later group failed".into(),
                    ))
                })
            }));
            if panics {
                assert!(result.is_err());
            } else {
                assert!(result.unwrap().is_err());
            }
            let state = cache.state.lock().unwrap();
            assert!(state.entries.contains_key(keys[0].as_ref().unwrap()));
            assert!(!state.entries.contains_key(keys[1].as_ref().unwrap()));
            drop(state);
            // A peer retaining the removed cell can claim it after either
            // normal error cleanup or unwinding; it never remains running.
            assert!(unfinished.unwrap().claim().is_some());
            let values = cache
                .run_batch(&keys, |missing, publish| {
                    assert_eq!(missing, [1]);
                    publish(1, vec!["retried".into()])
                })
                .unwrap();
            assert_eq!(values, [vec!["complete"], vec!["retried"]]);
        }
    }

    #[test]
    fn batch_rejects_missing_duplicate_and_unowned_publications() {
        for invalid in ["missing", "duplicate", "unowned", "out of range"] {
            let cache = ForcerRunCache::default();
            cache
                .run(Some(key("ready")), || Ok(vec!["ready".into()]))
                .unwrap();
            let keys = [Some(key("owned")), Some(key("ready")), Some(key("pending"))];
            let result = cache.run_batch(&keys, |_, publish| {
                if invalid == "missing" {
                    return Ok(());
                }
                publish(0, vec!["owned".into()])?;
                let index = match invalid {
                    "duplicate" => 0,
                    "unowned" => 1,
                    _ => 3,
                };
                // A callback cannot hide an invalid emission by ignoring it.
                assert!(publish(index, vec!["invalid".into()]).is_err());
                publish(2, vec!["pending".into()])?;
                Ok(())
            });
            assert!(result.is_err(), "{invalid}");
            cache
                .run(Some(key("pending")), || Ok(vec!["pending".into()]))
                .unwrap();
        }
    }

    #[test]
    fn mixed_cached_and_uncached_requests_preserve_empty_outputs_and_order() {
        let cache = ForcerRunCache::default();
        cache
            .run(Some(key("ready")), || Ok(vec!["ready".into()]))
            .unwrap();
        let keys = [
            None,
            Some(key("ready")),
            Some(key("empty")),
            None,
            Some(key("empty")),
        ];
        let outputs = cache
            .run_batch(&keys, |missing, publish| {
                assert_eq!(missing, [0, 2, 3]);
                publish(3, vec!["second uncached".into()])?;
                publish(2, Vec::new())?;
                publish(0, vec!["first uncached".into()])
            })
            .unwrap();
        assert_eq!(
            outputs,
            [
                vec!["first uncached"],
                vec!["ready"],
                vec![],
                vec!["second uncached"],
                vec![]
            ]
        );
    }

    #[test]
    fn oversized_publication_reaches_callers_holding_the_evicted_cell() {
        let cache = ForcerRunCache::new(4096, 2);
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (release_tx, release_rx) = std::sync::mpsc::channel();
        let output = "x".repeat(8192);
        std::thread::scope(|scope| {
            let cache = &cache;
            let output = &output;
            let owner = scope.spawn(move || {
                cache.run(Some(key("large")), || {
                    started_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                    Ok(vec![output.clone()])
                })
            });
            started_rx.recv().unwrap();
            let cell = cache.state.lock().unwrap().entries[&key("large")]
                .output
                .clone();
            release_tx.send(()).unwrap();
            owner.join().unwrap().unwrap();
            cell.wait();
            assert_eq!(
                cell.output.get().unwrap().as_ref(),
                std::slice::from_ref(output)
            );
            let state = cache.state.lock().unwrap();
            assert!(!state.entries.contains_key(&key("large")));
            assert!(state.bytes <= 4096);
        });
    }
}
