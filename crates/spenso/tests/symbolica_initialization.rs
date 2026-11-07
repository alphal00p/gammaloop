#![cfg(feature = "symbolica")]

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, RecvTimeoutError, SyncSender};
use std::sync::{Mutex, OnceLock};
use std::thread;
use std::time::Duration;

use spenso::network::library::symbolic::ETS;
use symbolica::initialize;
use symbolica::state::State;

struct InitializerGate {
    entered: SyncSender<()>,
    release: Mutex<Receiver<()>>,
}

static GATE: OnceLock<InitializerGate> = OnceLock::new();
static CALLBACK_RETURNED: AtomicBool = AtomicBool::new(false);

// A downstream initializer is allowed to access already warmed bundles without
// Spenso's local marker. Such a reentrant access must not certify completion.
initialize!(
    || {
        let _ = ETS.metric;
        let gate = GATE
            .get()
            .expect("test installs the initializer gate first");
        gate.entered.send(()).unwrap();
        gate.release
            .lock()
            .unwrap()
            .recv_timeout(Duration::from_secs(10))
            .expect("test releases the initializer");
        CALLBACK_RETURNED.store(true, Ordering::Release);
    },
    "spenso",
);

// This binary intentionally has one test: Symbolica initialization is process
// global and cannot be reset between lifecycle scenarios.
#[test]
fn public_bundle_access_waits_for_all_initializers() {
    let (entered_tx, entered_rx) = mpsc::sync_channel(0);
    let (release_tx, release_rx) = mpsc::sync_channel(0);
    assert!(
        GATE.set(InitializerGate {
            entered: entered_tx,
            release: Mutex::new(release_rx),
        })
        .is_ok()
    );

    let initializer = thread::spawn(|| {
        let _ = State::is_builtin("__spenso_lifecycle_test__");
    });
    entered_rx
        .recv_timeout(Duration::from_secs(10))
        .expect("downstream callback has read the warmed bundle");

    let (started_tx, started_rx) = mpsc::sync_channel(0);
    let (finished_tx, finished_rx) = mpsc::sync_channel(0);
    let reader = thread::spawn(move || {
        started_tx.send(()).unwrap();
        let _ = ETS.metric;
        let completed = CALLBACK_RETURNED.load(Ordering::Acquire);
        finished_tx.send(completed).unwrap();
    });
    started_rx.recv().unwrap();
    let premature = finished_rx.recv_timeout(Duration::from_millis(200));
    // Always release and join the initializer before reporting the failure.
    release_tx.send(()).unwrap();
    initializer.join().unwrap();
    assert!(matches!(premature, Err(RecvTimeoutError::Timeout)));
    assert!(
        finished_rx
            .recv_timeout(Duration::from_secs(10))
            .expect("reader resumes after global initialization completes")
    );
    reader.join().unwrap();
}
