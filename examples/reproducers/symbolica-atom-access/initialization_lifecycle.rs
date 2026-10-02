//! Run each mode in a fresh process against matching patched Symbolica/Spenso.
//! No global-state reset is used: it cannot reset Symbolica's initializer Once.

use std::{
    panic::{AssertUnwindSafe, catch_unwind},
    sync::{
        LazyLock, Mutex, OnceLock,
        atomic::{AtomicBool, AtomicUsize, Ordering},
        mpsc::{Receiver, RecvTimeoutError, SyncSender, sync_channel},
    },
    thread,
    time::Duration,
};

use spenso::{
    network::library::symbolic::ETS,
    symbolica_init::{SymbolicaInitLazy, in_symbolica_initializer},
};
use symbolica::{initialize, state::State};

#[derive(Clone, Copy, Debug)]
enum Mode {
    BundleFirst,
    StateFirst,
    Concurrent,
    InitializerPanic,
    BundlePanic,
    MarkerUnwind,
}

struct Gate {
    entered: SyncSender<()>,
    release: Mutex<Receiver<()>>,
}

static MODE: OnceLock<Mode> = OnceLock::new();
static GATE: OnceLock<Gate> = OnceLock::new();
static CALLBACK_COUNT: AtomicUsize = AtomicUsize::new(0);
static CALLBACK_RETURNED: AtomicBool = AtomicBool::new(false);
static BUNDLE_ATTEMPTS: AtomicUsize = AtomicUsize::new(0);

static PANICKING_INNER: LazyLock<()> = LazyLock::new(|| {
    BUNDLE_ATTEMPTS.fetch_add(1, Ordering::Relaxed);
    panic!("intentional independent bundle initializer panic");
});
static PANICKING: SymbolicaInitLazy<()> = SymbolicaInitLazy::new(&PANICKING_INNER);

// Deliberately do not enter Spenso's initializer marker. This is a legal
// downstream Symbolica callback, and the dependency has already warmed ETS.
initialize!(
    || {
        CALLBACK_COUNT.fetch_add(1, Ordering::Relaxed);
        assert!(!State::is_initialized());
        let _ = ETS.metric;
        assert!(!State::is_initialized());

        match MODE.get().expect("mode set before any Symbolica access") {
            Mode::Concurrent => {
                let gate = GATE.get().unwrap();
                gate.entered.send(()).unwrap();
                gate.release
                    .lock()
                    .unwrap()
                    .recv_timeout(Duration::from_secs(10))
                    .expect("main thread releases initializer");
            }
            Mode::InitializerPanic => panic!("intentional registered initializer panic"),
            _ => {}
        }
        CALLBACK_RETURNED.store(true, Ordering::Release);
    },
    "spenso",
);

fn start_through_state() {
    let _ = State::is_builtin("__spenso_completion_lifecycle_probe__");
}

fn concurrent_access() {
    let (entered_tx, entered_rx) = sync_channel(0);
    let (release_tx, release_rx) = sync_channel(0);
    assert!(
        GATE.set(Gate {
            entered: entered_tx,
            release: Mutex::new(release_rx),
        })
        .is_ok()
    );

    let initializer = thread::spawn(start_through_state);
    entered_rx
        .recv_timeout(Duration::from_secs(10))
        .expect("unmarked initializer entered after reading ETS");
    assert!(!State::is_initialized());

    let (started_tx, started_rx) = sync_channel(0);
    let (finished_tx, finished_rx) = sync_channel(0);
    let reader = thread::spawn(move || {
        started_tx.send(()).unwrap();
        let _ = ETS.metric;
        assert!(CALLBACK_RETURNED.load(Ordering::Acquire));
        assert!(State::is_initialized());
        finished_tx.send(()).unwrap();
    });
    started_rx.recv().unwrap();
    // This negative timeout is a concurrency regression check, not a timing
    // benchmark. The initializer cannot return before main releases it.
    assert!(matches!(
        finished_rx.recv_timeout(Duration::from_millis(200)),
        Err(RecvTimeoutError::Timeout)
    ));

    release_tx.send(()).unwrap();
    initializer.join().unwrap();
    finished_rx
        .recv_timeout(Duration::from_secs(10))
        .expect("reader resumes after every initializer completes");
    reader.join().unwrap();
}

fn main() {
    let mode = match std::env::args().nth(1).as_deref() {
        Some("bundle-first") => Mode::BundleFirst,
        Some("state-first") => Mode::StateFirst,
        Some("concurrent") => Mode::Concurrent,
        Some("initializer-panic") => Mode::InitializerPanic,
        Some("bundle-panic") => Mode::BundlePanic,
        Some("marker-unwind") => Mode::MarkerUnwind,
        _ => panic!(
            "mode: bundle-first | state-first | concurrent | initializer-panic | \
             bundle-panic | marker-unwind"
        ),
    };
    MODE.set(mode).unwrap();
    assert!(!State::is_initialized());

    match mode {
        Mode::BundleFirst => {
            let _ = ETS.metric;
        }
        Mode::StateFirst => start_through_state(),
        Mode::Concurrent => concurrent_access(),
        Mode::InitializerPanic => {
            assert!(catch_unwind(start_through_state).is_err());
            assert!(!State::is_initialized());
            assert!(
                catch_unwind(|| {
                    let _ = ETS.metric;
                })
                .is_err()
            );
            assert!(!State::is_initialized());
            assert!(!CALLBACK_RETURNED.load(Ordering::Acquire));
            assert_eq!(CALLBACK_COUNT.load(Ordering::Relaxed), 1);
            println!("PASS {mode:?}");
            return;
        }
        Mode::BundlePanic => {
            start_through_state();
            for _ in 0..2 {
                assert!(
                    catch_unwind(AssertUnwindSafe(|| {
                        let _ = &*PANICKING;
                    }))
                    .is_err()
                );
                assert!(State::is_initialized());
            }
            assert_eq!(BUNDLE_ATTEMPTS.load(Ordering::Relaxed), 1);
        }
        Mode::MarkerUnwind => {
            assert!(
                catch_unwind(|| {
                    in_symbolica_initializer(|| {
                        in_symbolica_initializer(|| panic!("intentional marker unwind"));
                    });
                })
                .is_err()
            );
            assert!(!State::is_initialized());
            let _ = ETS.metric;
        }
    }

    assert!(State::is_initialized());
    assert!(CALLBACK_RETURNED.load(Ordering::Acquire));
    assert_eq!(CALLBACK_COUNT.load(Ordering::Relaxed), 1);
    // Repeated public access remains valid after completion.
    for _ in 0..4 {
        let _ = ETS.metric;
    }
    println!("PASS {mode:?}");
}
