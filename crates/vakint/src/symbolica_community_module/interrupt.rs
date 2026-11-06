//! One interruptible runner for the potentially long native calls.
//!
//! The native work runs on a worker thread with the GIL released. The calling
//! thread checks Python signals every 25 ms (Ctrl-C, marimo's Interrupt,
//! `monitor="tui"`'s re-sent SIGINT) and delivers progress. An interruption or
//! a failing progress callback cancels the request ([`crate::Cancellation`]):
//! in [`Scope::Everywhere`] every Symbolica resource checkpoint of the worker
//! (and of the workers it dispatches) then fails; in [`Scope::Coarse`] (data
//! loads, whose dense Symbolica checkpoints made a per-checkpoint budget cost
//! 28% of a four-loop load) Vakint's own checkpoints do (per file, program and
//! pack). The worker is joined, and the ORIGINAL exception is raised.
//!
//! The cancellation scope costs a little at every checkpoint, so it is
//! installed only when something can cancel the call: a caller on Python's
//! main thread (the only thread Python delivers signals to) or a progress
//! callback. Other callers (a marimo app-mode kernel thread, a user thread)
//! run without it, on the calling thread when no progress is consumed. The
//! engine's `interrupt` setting ("fine", "coarse" or "off") selects the scope.
//! The request threads are kept for later calls (see `pool`). A coarse call
//! without a progress consumer (a data load in a script) stays on the calling
//! thread like the unlicensed path below: a load allocates gigabytes, which is
//! about 10% slower on another thread's allocator arena.
//!
//! Restricted (unlicensed) Symbolica allows one Symbolica thread per user, and
//! the calling thread already is it. Then the work stays on the calling thread
//! (still without the GIL); a ticker marks a SIGINT check as due every 25 ms,
//! and the next checkpoint consumes a pending SIGINT and cancels. The signal is
//! re-delivered after the work returned, so Python's own handler raises
//! (KeyboardInterrupt, or marimo's interruption). No Python code runs inside a
//! checkpoint; live progress needs a license (the final snapshot is delivered).

use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use std::time::Duration;

use pyo3::{
    PyErr, PyResult, Python,
    exceptions::{PyException, PyRuntimeError, PyValueError},
    types::PyAnyMethods,
};
use symbolica::license::LicenseManager;
use symbolica::utils::resources::ResourceBudget;

use super::{progress::Progress, vakint_to_python_error};
use crate::{Cancellation, ResourceLimits, VakintError};

/// Period of the calling thread's signal checks.
const SIGNAL_POLL: Duration = Duration::from_millis(25);

/// Smallest stack of the request thread. The work used to run on the calling
/// thread (at least 8 MiB on Python's main thread); Symbolica's normalization
/// recurses deeply on nested input. Only reserved, not committed.
const MIN_REQUEST_STACK_BYTES: usize = 64 << 20;
/// Request thread stack when the caller's stack limit is unlimited.
const UNLIMITED_REQUEST_STACK_BYTES: usize = 1 << 30;

/// Where an interruption stops the native work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Scope {
    /// At every resource checkpoint (reductions, evaluations).
    Everywhere,
    /// At Vakint's own checkpoints only (data loads: per file, program and
    /// pack; a long native load step in between finishes first).
    Coarse,
}

/// The engine's `interrupt` setting.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum InterruptMode {
    /// At every resource checkpoint (data loads: at Vakint's checkpoints).
    #[default]
    Fine,
    /// At Vakint's own checkpoints only (per data file, program step and
    /// expression term), at no cost to the native checkpoints.
    Coarse,
    /// Never during a call: a Ctrl-C raises once the call returned.
    Off,
}

impl InterruptMode {
    pub(super) fn parse(value: &str) -> PyResult<Self> {
        match value {
            "fine" => Ok(Self::Fine),
            "coarse" => Ok(Self::Coarse),
            "off" => Ok(Self::Off),
            _ => Err(PyValueError::new_err(format!(
                "interrupt must be 'fine', 'coarse' or 'off', not {value:?}"
            ))),
        }
    }

    pub(super) fn as_str(self) -> &'static str {
        match self {
            Self::Fine => "fine",
            Self::Coarse => "coarse",
            Self::Off => "off",
        }
    }
}

/// Per-engine settings of the runner.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Options {
    /// The engine's time and memory limits.
    pub limits: ResourceLimits,
    pub mode: InterruptMode,
}

/// Run `operation` under the engine's resource limits, interruptibly.
///
/// The operation must not touch Python objects: it runs while the calling
/// thread does not hold the GIL. Its error is converted on the calling thread.
pub(super) fn run<T: Send>(
    py: Python<'_>,
    options: Options,
    progress: Progress,
    operation: impl FnOnce() -> Result<T, VakintError> + Send,
) -> PyResult<T> {
    run_in(py, Scope::Everywhere, options, progress, operation)
}

/// [`run`] for data loads: interruptible at Vakint's coarse checkpoints.
pub(super) fn run_load<T: Send>(
    py: Python<'_>,
    options: Options,
    progress: Progress,
    operation: impl FnOnce() -> Result<T, VakintError> + Send,
) -> PyResult<T> {
    run_in(py, Scope::Coarse, options, progress, operation)
}

fn enter<R>(cancellation: &Cancellation, scope: Option<Scope>, operation: impl FnOnce() -> R) -> R {
    match scope {
        Some(Scope::Everywhere) => cancellation.run(operation),
        Some(Scope::Coarse) => cancellation.run_coarse(operation),
        None => operation(),
    }
}

/// Whether Python delivers signals to this thread (its main thread).
fn on_signal_thread(py: Python<'_>) -> bool {
    let main = || -> PyResult<bool> {
        let threading = py.import("threading")?;
        let main = threading.call_method0("main_thread")?.getattr("ident")?;
        let current = threading.call_method0("get_ident")?;
        main.eq(current)
    };
    main().unwrap_or(false)
}

/// Whether `error` is an ordinary exception (not KeyboardInterrupt,
/// SystemExit or another BaseException that must propagate).
fn is_ordinary(py: Python<'_>, error: &PyErr) -> bool {
    error.is_instance_of::<PyException>(py)
}

fn run_in<T: Send>(
    py: Python<'_>,
    scope: Scope,
    options: Options,
    mut progress: Progress,
    operation: impl FnOnce() -> Result<T, VakintError> + Send,
) -> PyResult<T> {
    // A Ctrl-C pending before the call raises without starting any work.
    if let Err(error) = py.check_signals() {
        let _ = progress.finish(py, Some(&error), true);
        return Err(error);
    }
    // Only a signal (main thread) or a raising callback can cancel a call;
    // without either the cancellation scope would cost without use.
    let signals = options.mode != InterruptMode::Off && on_signal_thread(py);
    let callback = options.mode != InterruptMode::Off && progress.is_callback();
    let scope = (signals || callback).then_some(match options.mode {
        InterruptMode::Coarse => Scope::Coarse,
        _ => scope,
    });
    let limits = options.limits;
    let request = || {
        // The engine's budget starts with the request, inside the
        // cancellation scope; its errors are unchanged.
        ResourceBudget::new(limits).run(|| {
            crate::execution::checkpoint()?;
            let result = operation();
            let checked = crate::execution::checkpoint();
            result.and_then(|value| checked.map(|()| value).map_err(VakintError::from))
        })
    };
    let (outcome, interruption) = if scope.is_none() && !progress.is_active() {
        // Nothing can interrupt the call and nothing watches it: run it here,
        // without the GIL, as an uninterruptible call always did.
        (Ok(py.detach(request)), None)
    } else if scope == Some(Scope::Coarse) && !progress.is_active() {
        // Coarse checks only and nothing to deliver (a data load in a
        // script): stay on this thread, whose allocator a load of gigabytes
        // runs about 10% faster on; its own checkpoints consume a SIGINT.
        let (outcome, interrupted) = py.detach(|| on_caller(scope, request));
        let interruption = interrupted.then(|| redeliver_sigint(py));
        (Ok(outcome), interruption)
    } else {
        let cancellation = Cancellation::default();
        let watched = py.detach(|| watch(&cancellation, &mut progress, signals, scope, request));
        match watched {
            Watched::Done(joined, interruption) => (joined, interruption),
            Watched::Restricted(request) => {
                // Unlicensed: the work stays on this (Symbolica's) thread.
                let (outcome, interrupted) =
                    py.detach(|| on_caller(scope.filter(|_| signals), request));
                let interruption = interrupted.then(|| redeliver_sigint(py));
                (Ok(outcome), interruption)
            }
            Watched::SpawnFailed(error) => {
                let error = PyRuntimeError::new_err(format!(
                    "cannot start the Vakint request thread: {error}"
                ));
                let _ = progress.finish(py, Some(&error), false);
                return Err(error);
            }
        }
    };
    let outcome = match outcome {
        Ok(result) => result,
        Err(panic) => {
            // Close the indicator; the panic itself propagates unchanged.
            let error = PyRuntimeError::new_err("the Vakint request panicked");
            let _ = progress.finish(py, Some(&error), false);
            std::panic::resume_unwind(panic)
        }
    };
    // An interruption wins over the worker's result (usually the error its
    // cancelled checkpoint returned).
    if let Some(error) = interruption {
        let _ = progress.finish(py, Some(&error), true);
        return Err(error);
    }
    if let Err(error) = py.check_signals() {
        let _ = progress.finish(py, Some(&error), true);
        return Err(error);
    }
    // MissingMasterOrdersError attributes need the GIL: convert here.
    match outcome.map_err(vakint_to_python_error) {
        Ok(value) => {
            // The result is computed: a failing final snapshot (callback or
            // indicator) is reported, not allowed to discard it. A Ctrl-C
            // arriving there still propagates, as in plain Python.
            if let Err(error) = progress.finish(py, None, false) {
                if !is_ordinary(py, &error) {
                    return Err(error);
                }
                error.write_unraisable(py, None);
            }
            Ok(value)
        }
        Err(error) => {
            let _ = progress.finish(py, Some(&error), false);
            Err(error)
        }
    }
}

type Joined<T> = std::thread::Result<Result<T, VakintError>>;

/// Outcome of [`watch`].
enum Watched<T, R> {
    /// The worker ran the request; the interruption that cancelled it, if any.
    Done(Joined<T>, Option<PyErr>),
    /// The worker may not use Symbolica (no license on it): the request is
    /// handed back untouched, to run on the calling thread.
    Restricted(R),
    SpawnFailed(std::io::Error),
}

/// Whether threads other than the caller may use Symbolica. Asked on the
/// request thread itself, which holds no library unlock of the caller: only a
/// process-wide license counts. A license, once seen, stays.
fn symbolica_threads_allowed() -> bool {
    static LICENSED: AtomicBool = AtomicBool::new(false);
    if LICENSED.load(Ordering::Acquire) {
        return true;
    }
    let licensed = LicenseManager::is_licensed();
    if licensed {
        LICENSED.store(true, Ordering::Release);
    }
    licensed
}

/// Stack of the request thread: at least the caller's stack limit (the work
/// used to run on the calling thread), at least [`MIN_REQUEST_STACK_BYTES`]
/// and at least `RUST_MIN_STACK`.
fn request_stack_bytes() -> usize {
    #[cfg(all(unix, not(target_arch = "wasm32")))]
    let limit = match rustix::process::getrlimit(rustix::process::Resource::Stack).current {
        Some(bytes) => usize::try_from(bytes).unwrap_or(UNLIMITED_REQUEST_STACK_BYTES),
        None => UNLIMITED_REQUEST_STACK_BYTES,
    };
    #[cfg(not(all(unix, not(target_arch = "wasm32"))))]
    let limit = 0;
    let configured = std::env::var("RUST_MIN_STACK")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    limit
        .clamp(MIN_REQUEST_STACK_BYTES, UNLIMITED_REQUEST_STACK_BYTES)
        .max(configured)
}

/// Run `request` on the calling thread (which must not hold the GIL). With a
/// scope, cancel at the first checkpoint after a SIGINT arrived. Returns the
/// result and whether a SIGINT cancelled it (the SIGINT is then consumed).
fn on_caller<T: Send>(
    scope: Option<Scope>,
    request: impl FnOnce() -> Result<T, VakintError> + Send,
) -> (Result<T, VakintError>, bool) {
    use std::sync::mpsc::{RecvTimeoutError, sync_channel};
    if scope.is_none() {
        return (request(), false);
    }
    let due = Arc::new(AtomicBool::new(false));
    let interrupted = Arc::new(AtomicBool::new(false));
    let caller = std::thread::current().id();
    let poll = {
        let (due, interrupted) = (Arc::clone(&due), Arc::clone(&interrupted));
        move || {
            if !due.load(Ordering::Relaxed) || std::thread::current().id() != caller {
                return false;
            }
            due.store(false, Ordering::Relaxed);
            // Consumes a pending SIGINT without running Python code here
            // (main thread only, else 0); re-delivered by the caller.
            let pending = Python::attach(|_| unsafe { pyo3::ffi::PyOS_InterruptOccurred() } != 0);
            if pending {
                interrupted.store(true, Ordering::Relaxed);
            }
            pending
        }
    };
    let cancellation = Cancellation::with_poll(poll);
    let outcome = std::thread::scope(|threads| {
        let (stop, stopped) = sync_channel::<()>(0);
        // The ticker never uses Symbolica, so it needs no thread permit.
        let ticker = Arc::clone(&due);
        threads.spawn(move || {
            while let Err(RecvTimeoutError::Timeout) = stopped.recv_timeout(SIGNAL_POLL) {
                ticker.store(true, Ordering::Relaxed);
            }
        });
        let outcome = enter(&cancellation, scope, request);
        drop(stop);
        outcome
    });
    (outcome, interrupted.load(Ordering::Relaxed))
}

/// Deliver a SIGINT that a checkpoint consumed to Python's own handler, at a
/// safe point on the calling thread, and return the exception it raised.
fn redeliver_sigint(py: Python<'_>) -> PyErr {
    // Async-signal-safe; the handler runs at the next signal check.
    unsafe { pyo3::ffi::PyErr_SetInterrupt() };
    py.check_signals().err().unwrap_or_else(|| {
        // The installed handler did not raise: the request is cancelled
        // nevertheless, so report it as an interruption.
        pyo3::exceptions::PyKeyboardInterrupt::new_err("Vakint request interrupted")
    })
}

/// Run `request` on a pooled request thread inside `cancellation`'s `scope`;
/// poll signals (when `signals`) and progress until it returns.
#[cfg(not(target_arch = "wasm32"))]
fn watch<T: Send, R: FnOnce() -> Result<T, VakintError> + Send>(
    cancellation: &Cancellation,
    progress: &mut Progress,
    signals: bool,
    scope: Option<Scope>,
    request: R,
) -> Watched<T, R> {
    let mut interruption = None;
    let ran = pool::run(
        move || {
            if !symbolica_threads_allowed() {
                return Err(request);
            }
            Ok(enter(cancellation, scope, request))
        },
        || {
            if interruption.is_some() {
                return;
            }
            let polled = Python::attach(|py| {
                if signals {
                    py.check_signals()?;
                }
                progress.deliver(py)
            });
            if let Err(error) = polled {
                interruption = Some(error);
                if scope.is_some() {
                    cancellation.cancel();
                    // Show that the stop is under way (a native step
                    // without checkpoints finishes first).
                    Python::attach(|py| progress.stopping(py));
                }
            }
        },
    );
    match ran {
        Ok(Ok(Ok(joined))) => Watched::Done(Ok(joined), interruption),
        Ok(Ok(Err(request))) => Watched::Restricted(request),
        Ok(Err(panic)) => Watched::Done(Err(panic), interruption),
        Err(error) => Watched::SpawnFailed(error),
    }
}

/// Persistent request threads. A thread spawned for each call starts with a
/// cold allocator arena and empty thread-local Symbolica workspaces, which
/// measurably slowed every call (about 3% of a four-loop evaluation and 8%
/// of its data load); an idle request thread is kept for the next call
/// instead (one per concurrent caller, at most [`MAX_IDLE`]).
#[cfg(not(target_arch = "wasm32"))]
mod pool {
    use std::panic::{AssertUnwindSafe, catch_unwind};
    use std::sync::mpsc::{self, RecvTimeoutError};
    use std::sync::{Mutex, PoisonError};

    use super::{SIGNAL_POLL, request_stack_bytes};

    /// Idle request threads kept for later calls.
    const MAX_IDLE: usize = 4;

    type Job = Box<dyn FnOnce() + Send + 'static>;

    /// An idle request thread, owned by the process that spawned it (a
    /// forked child has no such thread).
    struct Worker {
        pid: u32,
        jobs: mpsc::Sender<Job>,
    }

    static IDLE: Mutex<Vec<Worker>> = Mutex::new(Vec::new());

    fn spawn() -> std::io::Result<Worker> {
        let (jobs, received) = mpsc::channel::<Job>();
        std::thread::Builder::new()
            .name("vakint-request".into())
            .stack_size(request_stack_bytes())
            // A job never unwinds (it catches its work's panic); the thread
            // ends when its worker is dropped.
            .spawn(move || received.into_iter().for_each(|job| job()))?;
        Ok(Worker {
            pid: std::process::id(),
            jobs,
        })
    }

    fn take() -> std::io::Result<Worker> {
        let pid = std::process::id();
        let idle = {
            let mut idle = IDLE.lock().unwrap_or_else(PoisonError::into_inner);
            idle.retain(|worker| worker.pid == pid);
            idle.pop()
        };
        idle.map_or_else(spawn, Ok)
    }

    fn give_back(worker: Worker) {
        let mut idle = IDLE.lock().unwrap_or_else(PoisonError::into_inner);
        if idle.len() < MAX_IDLE {
            idle.push(worker);
        }
    }

    /// The work of one call and what it must release, in this order, when it
    /// ran or was dropped unrun: the work (and everything it borrows), then
    /// the completion signal.
    struct Packet<F> {
        work: Option<F>,
        done: Option<mpsc::Sender<()>>,
    }

    impl<F> Drop for Packet<F> {
        fn drop(&mut self) {
            drop(self.work.take());
            drop(self.done.take());
        }
    }

    /// Waits, also while unwinding, until the job's packet is gone.
    struct Completion(mpsc::Receiver<()>);

    impl Completion {
        /// Poll `tick` every [`SIGNAL_POLL`] until the job is gone.
        fn wait(&self, mut tick: impl FnMut()) {
            while let Err(RecvTimeoutError::Timeout) = self.0.recv_timeout(SIGNAL_POLL) {
                tick();
            }
        }
    }

    impl Drop for Completion {
        fn drop(&mut self) {
            // Never sent: returns when the packet's sender is dropped.
            let _ = self.0.recv();
        }
    }

    /// Run `work` on a request thread, calling `tick` on this thread every
    /// [`SIGNAL_POLL`] until it returns: [`std::thread::scope`] with one
    /// pooled thread. Its panic is returned, not resumed.
    pub(super) fn run<O: Send>(
        work: impl FnOnce() -> O + Send,
        tick: impl FnMut(),
    ) -> std::io::Result<std::thread::Result<O>> {
        let slot: Mutex<Option<std::thread::Result<O>>> = Mutex::new(None);
        let (done, finished) = mpsc::channel::<()>();
        let completion = Completion(finished);
        let job: Box<dyn FnOnce() + Send + '_> = {
            let slot = &slot;
            let mut packet = Packet {
                work: Some(work),
                done: Some(done),
            };
            Box::new(move || {
                if let Some(work) = packet.work.take() {
                    let outcome = catch_unwind(AssertUnwindSafe(work));
                    *slot.lock().unwrap_or_else(PoisonError::into_inner) = Some(outcome);
                }
                drop(packet);
            })
        };
        // SAFETY: the job borrows `slot` and whatever `work` borrows, all of
        // which outlive this call: `completion` blocks this frame, on return
        // and on unwind, until the packet (and with it every borrow) is
        // dropped, whether the job ran or was dropped unrun by a thread that
        // exited. Only the type's lifetime changes.
        let job: Job = unsafe { std::mem::transmute::<Box<dyn FnOnce() + Send + '_>, Job>(job) };
        let mut worker = take()?;
        if let Err(mpsc::SendError(job)) = worker.jobs.send(job) {
            // That thread is gone: hand the job to a fresh one.
            worker = match spawn() {
                Ok(fresh) => fresh,
                Err(error) => {
                    drop(job);
                    drop(completion);
                    return Err(error);
                }
            };
            if let Err(mpsc::SendError(job)) = worker.jobs.send(job) {
                drop(job);
            }
        }
        completion.wait(tick);
        drop(completion);
        let outcome = slot.into_inner().unwrap_or_else(PoisonError::into_inner);
        match outcome {
            Some(Ok(value)) => {
                give_back(worker);
                Ok(Ok(value))
            }
            // A thread whose work panicked is not reused.
            Some(Err(panic)) => Ok(Err(panic)),
            None => Err(std::io::Error::other("the Vakint request thread exited")),
        }
    }
}

/// Browser kernels have no threads: run inline and poll signals only before
/// and after the request.
#[cfg(target_arch = "wasm32")]
fn watch<T: Send, R: FnOnce() -> Result<T, VakintError> + Send>(
    cancellation: &Cancellation,
    _progress: &mut Progress,
    _signals: bool,
    scope: Option<Scope>,
    request: R,
) -> Watched<T, R> {
    Watched::Done(
        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            enter(cancellation, scope, request)
        })),
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::{PyAnyMethods, PyDict, PyDictMethods};
    use std::ffi::CString;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::time::Instant;

    /// Threads of this process.
    fn thread_count() -> usize {
        std::fs::read_dir("/proc/self/task").map_or(0, |tasks| tasks.count())
    }

    /// Checkpoint until the request is cancelled (or a minute passed).
    fn spin_until_cancelled(entered: &AtomicBool) -> Result<(), VakintError> {
        entered.store(true, Ordering::Release);
        let started = Instant::now();
        while started.elapsed() < Duration::from_secs(60) {
            // Vakint's checkpoint: observed by both scopes.
            crate::execution::checkpoint()?;
            std::thread::sleep(Duration::from_micros(200));
        }
        Ok(())
    }

    /// Wait until the threads of this process are back to `baseline`
    /// (an exited thread can stay listed for a moment after its join).
    fn assert_threads_return_to(baseline: usize) {
        let started = Instant::now();
        while thread_count() > baseline && started.elapsed() < Duration::from_secs(2) {
            std::thread::sleep(Duration::from_millis(10));
        }
        assert!(
            thread_count() <= baseline,
            "{} > {baseline}",
            thread_count()
        );
    }

    /// Python handles signals on the thread that initialized it, and libtest
    /// runs several tests per process. A signal test therefore runs alone:
    /// under nextest (one process per test) directly, otherwise in a child
    /// process that runs only this test.
    fn on_python_main_thread(test: &str, body: impl FnOnce(Python<'_>)) {
        let alone = std::env::var_os("VAKINT_INTERRUPT_TEST_CHILD").is_some()
            || std::env::var("NEXTEST_EXECUTION_MODE").is_ok_and(|mode| mode == "process-per-test");
        if !alone {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args([test, "--exact", "--test-threads=1", "--nocapture"])
                .env("VAKINT_INTERRUPT_TEST_CHILD", "1")
                .status()
                .unwrap();
            assert!(status.success(), "{test} failed in its child process");
            return;
        }
        Python::initialize();
        Python::attach(|py| {
            // Python's own handler (the embedded interpreter installs none);
            // this also fails unless this thread initialized Python.
            py.run(
                c"import signal\nsignal.signal(signal.SIGINT, signal.default_int_handler)",
                None,
                None,
            )
            .expect("signal tests run on the thread that initialized Python");
            body(py)
        });
    }

    #[test]
    fn interrupt_reraises_original_exception_and_joins_worker() {
        on_python_main_thread(
            "symbolica_community_module::interrupt::tests::interrupt_reraises_original_exception_and_joins_worker",
            |py| {
                // The first call leaves an idle request thread for later calls.
                run(py, Options::default(), Progress::none(), || Ok(())).unwrap();
                let baseline = thread_count();
                let entered = Arc::new(AtomicBool::new(false));
                let signalled = Arc::new(std::sync::OnceLock::new());
                let trigger = {
                    let (entered, signalled) = (Arc::clone(&entered), Arc::clone(&signalled));
                    std::thread::spawn(move || {
                        while !entered.load(Ordering::Acquire) {
                            std::thread::sleep(Duration::from_millis(1));
                        }
                        std::thread::sleep(Duration::from_millis(100));
                        signalled.set(Instant::now()).unwrap();
                        // Async-signal-safe; needs no thread state.
                        unsafe { pyo3::ffi::PyErr_SetInterrupt() };
                    })
                };
                let result = run(py, Options::default(), Progress::none(), || {
                    spin_until_cancelled(&entered)
                });
                let returned = Instant::now();
                trigger.join().unwrap();
                let error = result.unwrap_err();
                assert!(
                    error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py),
                    "{error}"
                );
                let latency = returned.duration_since(*signalled.get().unwrap());
                assert!(latency < Duration::from_millis(300), "{latency:?}");
                // The worker (and the trigger) are joined.
                assert_threads_return_to(baseline);
                // Loads: the coarse scope stops at Vakint's checkpoints only.
                let entered = AtomicBool::new(false);
                let started = Instant::now();
                let interrupter = std::thread::spawn(|| {
                    std::thread::sleep(Duration::from_millis(150));
                    unsafe { pyo3::ffi::PyErr_SetInterrupt() };
                });
                let error = run_load(py, Options::default(), Progress::none(), || {
                    // Symbolica's own checkpoints are not observed here.
                    for _ in 0..1000 {
                        symbolica::utils::resources::checkpoint()?;
                    }
                    spin_until_cancelled(&entered)
                })
                .unwrap_err();
                interrupter.join().unwrap();
                assert!(error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py));
                assert!(started.elapsed() < Duration::from_millis(600));
                // The engine budget's own errors are unchanged without a signal.
                let error = run(
                    py,
                    Options {
                        limits: ResourceLimits {
                            max_elapsed: Some(Duration::ZERO),
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                    Progress::none(),
                    || Ok(()),
                )
                .unwrap_err();
                assert!(error.to_string().contains("time budget"), "{error}");
            },
        );
    }

    #[test]
    fn restricted_symbolica_path_consumes_sigint_at_a_checkpoint_and_redelivers_it() {
        on_python_main_thread(
            "symbolica_community_module::interrupt::tests::restricted_symbolica_path_consumes_sigint_at_a_checkpoint_and_redelivers_it",
            |py| {
                let entered = Arc::new(AtomicBool::new(false));
                let signalled = Arc::new(std::sync::OnceLock::new());
                let trigger = {
                    let (entered, signalled) = (Arc::clone(&entered), Arc::clone(&signalled));
                    std::thread::spawn(move || {
                        while !entered.load(Ordering::Acquire) {
                            std::thread::sleep(Duration::from_millis(1));
                        }
                        std::thread::sleep(Duration::from_millis(100));
                        signalled.set(Instant::now()).unwrap();
                        unsafe { pyo3::ffi::PyErr_SetInterrupt() };
                    })
                };
                let caller = std::thread::current().id();
                let (outcome, interrupted) = py.detach(|| {
                    on_caller(Some(Scope::Everywhere), || {
                        // The work stays on the calling thread.
                        assert_eq!(std::thread::current().id(), caller);
                        spin_until_cancelled(&entered)
                    })
                });
                let returned = Instant::now();
                trigger.join().unwrap();
                assert!(interrupted && outcome.is_err());
                let latency = returned.duration_since(*signalled.get().unwrap());
                assert!(latency < Duration::from_millis(300), "{latency:?}");
                // Python's own handler raises at the re-delivery, not before.
                let error = redeliver_sigint(py);
                assert!(error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py));
                py.check_signals().unwrap();
                // Without a signal the work runs to its end.
                let (outcome, interrupted) = py.detach(|| {
                    on_caller(Some(Scope::Coarse), || {
                        std::thread::sleep(Duration::from_millis(80));
                        crate::execution::checkpoint()?;
                        Ok(5)
                    })
                });
                assert!(!interrupted);
                assert_eq!(outcome.unwrap(), 5);
            },
        );
    }

    #[test]
    fn pending_interrupt_before_start_raises_without_running() {
        on_python_main_thread(
            "symbolica_community_module::interrupt::tests::pending_interrupt_before_start_raises_without_running",
            |py| {
                unsafe { pyo3::ffi::PyErr_SetInterrupt() };
                let called = AtomicBool::new(false);
                let error = run(py, Options::default(), Progress::none(), || {
                    called.store(true, Ordering::Relaxed);
                    Ok(())
                })
                .unwrap_err();
                assert!(error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py));
                assert!(!called.load(Ordering::Relaxed));
                // Nothing pending afterwards: the next call runs.
                assert_eq!(
                    run(py, Options::default(), Progress::none(), || Ok(7)).unwrap(),
                    7
                );
            },
        );
    }

    #[test]
    fn progress_callback_error_cancels_and_keeps_identity() {
        Python::initialize();
        Python::attach(|py| {
            let locals = PyDict::new(py);
            py.run(
                c"
class Stop(Exception):
    pass
calls = []
original = Stop('callback failed')
def callback(snapshot):
    calls.append(snapshot)
    if len(calls) == 2:
        raise original
",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();
            let callback = locals.get_item("callback").unwrap().unwrap().unbind();
            let progress = Progress::resolve(
                py,
                Some(&callback),
                "external_data_loading",
                "loading test data",
                Duration::ZERO,
            )
            .unwrap();
            let entered = AtomicBool::new(false);
            let started = Instant::now();
            let error = run(py, Options::default(), progress, || {
                spin_until_cancelled(&entered)
            })
            .unwrap_err();
            assert!(started.elapsed() < Duration::from_secs(5));
            let original = locals.get_item("original").unwrap().unwrap();
            assert!(error.value(py).is(&original));
            // Never called again after it raised.
            locals.set_item("error", error.value(py)).unwrap();
            py.run(
                c"
assert len(calls) == 2, calls
assert all(c.status == 'running' and c.stage == 'data_load' for c in calls), calls
assert calls[0].activity == 'loading test data'
assert calls[1].elapsed_seconds >= calls[0].elapsed_seconds + 0.2
",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();

            // A successful call ends with one final snapshot.
            let finals = PyDict::new(py);
            py.run(c"seen = []", Some(&finals), Some(&finals)).unwrap();
            let append = py
                .eval(c"seen.append", Some(&finals), Some(&finals))
                .unwrap()
                .unbind();
            let progress =
                Progress::resolve(py, Some(&append), "loading", "quick", Duration::ZERO).unwrap();
            assert_eq!(run(py, Options::default(), progress, || Ok(3)).unwrap(), 3);
            py.run(
                c"assert [s.status for s in seen][-1:] == ['succeeded'], seen",
                Some(&finals),
                Some(&finals),
            )
            .unwrap();
        });
    }

    #[test]
    fn auto_progress_uses_marimo_status_only_in_notebooks() {
        on_python_main_thread(
            "symbolica_community_module::interrupt::tests::auto_progress_uses_marimo_status_only_in_notebooks",
            |py| {
                let locals = PyDict::new(py);
                let code = |text: &str| CString::new(text).unwrap();
                let exec = |text: &str| py.run(&code(text), Some(&locals), Some(&locals));
                exec(
                    r#"
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock
context = MagicMock()
indicator = context.__enter__.return_value
marimo = SimpleNamespace(
    running_in_notebook=MagicMock(return_value=True),
    status=SimpleNamespace(
        spinner=MagicMock(return_value=context),
        progress_bar=MagicMock(return_value=context)),
    output=SimpleNamespace(_output=SimpleNamespace(flush=MagicMock())))
previous = sys.modules.get("marimo")
sys.modules["marimo"] = marimo
"#,
                )
                .unwrap();
                let result = (|| -> PyResult<()> {
                    let auto = pyo3::types::PyString::new(py, "auto").into_any().unbind();
                    let resolve = |stage| {
                        Progress::resolve(py, Some(&auto), stage, "testing", Duration::ZERO)
                    };
                    // A fast call never shows an indicator.
                    run(py, Options::default(), resolve("loading")?, || Ok(()))?;
                    exec("assert not marimo.status.spinner.called")?;
                    // A slow call shows a spinner after the delay, then exits it.
                    run(py, Options::default(), resolve("loading")?, || {
                        std::thread::sleep(Duration::from_millis(900));
                        Ok(())
                    })?;
                    exec(
                        r#"
assert marimo.status.spinner.call_count == 1, marimo.status.spinner.call_args_list
assert marimo.status.spinner.call_args.args[0] == "Vakint · testing"
assert context.__exit__.call_count == 1
context.__exit__.assert_called_with(None, None, None)
context.reset_mock(); marimo.status.spinner.reset_mock()
"#,
                    )?;
                    // A known Ordered phase total uses a bar.
                    let mut progress = resolve("")?;
                    let sink = progress.native_sink().unwrap();
                    run(py, Options::default(), progress, || {
                        let mut model = serde_json::json!({
                            "status": "Running", "elapsed_seconds": 1.0,
                            "stage": {"kind": "native", "id": "native", "activity": "ordered reduction"},
                            "ordered": {"phase": {"index": 2, "total": 5, "kind": "apply"}},
                        });
                        std::thread::sleep(Duration::from_millis(600));
                        sink.store(model.clone());
                        std::thread::sleep(Duration::from_millis(300));
                        model["ordered"]["phase"]["index"] = 3.into();
                        sink.store(model);
                        std::thread::sleep(Duration::from_millis(300));
                        Ok(())
                    })?;
                    exec(
                        r#"
assert marimo.status.progress_bar.call_count == 1, marimo.status.progress_bar.call_args_list
assert marimo.status.progress_bar.call_args.kwargs["total"] == 5
increments = [c.args[0] for c in indicator.update.call_args_list if c.args and isinstance(c.args[0], int)]
assert increments == [2, 1], indicator.update.call_args_list
context.reset_mock(); marimo.status.spinner.reset_mock(); marimo.status.progress_bar.reset_mock()
"#,
                    )?;
                    // An interrupt finishes the indicator with the original error.
                    let entered = AtomicBool::new(false);
                    let interrupter = std::thread::spawn(|| {
                        std::thread::sleep(Duration::from_millis(800));
                        unsafe { pyo3::ffi::PyErr_SetInterrupt() };
                    });
                    let error = run(py, Options::default(), resolve("loading")?, || {
                        spin_until_cancelled(&entered)
                    })
                    .unwrap_err();
                    interrupter.join().unwrap();
                    assert!(error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py));
                    locals.set_item("error", error.value(py))?;
                    exec(
                        r#"
assert context.__exit__.call_count == 1
assert context.__exit__.call_args.args[:2] == (KeyboardInterrupt, error)
# The indicator said the stop was under way before it closed.
assert any(c.args[:1] == ("Vakint · stopping",) for c in indicator.update.call_args_list), indicator.update.call_args_list
context.reset_mock(); marimo.status.spinner.reset_mock()
"#,
                    )?;
                    // Outside a notebook "auto" is no consumer at all.
                    exec("marimo.running_in_notebook.return_value = False")?;
                    assert!(!resolve("loading")?.is_active());
                    Ok(())
                })();
                exec(
                    r#"
if previous is None:
    del sys.modules["marimo"]
else:
    sys.modules["marimo"] = previous
"#,
                )
                .unwrap();
                result.unwrap();
            },
        );
    }

    #[test]
    fn runner_keeps_python_responsive_and_converts_errors_on_the_caller() {
        Python::initialize();
        Python::attach(|py| {
            // Another Python thread runs while the request holds no GIL.
            let locals = PyDict::new(py);
            py.run(
                c"
import threading, time
ticks = []
stop = threading.Event()
def tick():
    while not stop.is_set():
        ticks.append(time.monotonic())
        time.sleep(0.01)
thread = threading.Thread(target=tick)
thread.start()
",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();
            let counted = Arc::new(AtomicUsize::new(0));
            run(py, Options::default(), Progress::none(), || {
                std::thread::sleep(Duration::from_millis(500));
                counted.fetch_add(1, Ordering::Relaxed);
                Ok(())
            })
            .unwrap();
            py.run(
                c"stop.set(); thread.join()\nassert len(ticks) >= 20, len(ticks)",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();
            assert_eq!(counted.load(Ordering::Relaxed), 1);
            let error = run::<()>(py, Options::default(), Progress::none(), || {
                Err(VakintError::EvaluationError("native failure".into()))
            })
            .unwrap_err();
            assert!(error.is_instance_of::<pyo3::exceptions::PyValueError>(py));
            assert!(error.to_string().contains("Vakint error | "), "{error}");
        });
    }

    /// Recursion that needs about `depth` x 64 KiB of stack.
    #[inline(never)]
    fn deep(depth: usize) -> usize {
        let frame = std::hint::black_box([depth as u8; 64 << 10]);
        if depth == 0 {
            frame[0] as usize
        } else {
            deep(depth - 1) + std::hint::black_box(frame[depth % frame.len()]) as usize
        }
    }

    #[test]
    fn request_thread_has_at_least_the_callers_stack() {
        // The caller polls through Python while the work runs.
        Python::initialize();
        assert!(request_stack_bytes() >= MIN_REQUEST_STACK_BYTES);
        // About 20 MiB of stack: ten times Rust's default for a new thread.
        let cancellation = Cancellation::default();
        let mut progress = Progress::none();
        let watched = watch(
            &cancellation,
            &mut progress,
            false,
            Some(Scope::Everywhere),
            || Ok(deep(320)),
        );
        match watched {
            Watched::Done(Ok(Ok(value)), None) => assert!(value > 0),
            Watched::Restricted(_) => {} // unlicensed: the caller's stack is used
            _ => panic!("deep recursion failed on the request thread"),
        }
    }

    #[test]
    fn interrupt_modes_parse_and_print() {
        Python::initialize();
        for mode in [
            InterruptMode::Fine,
            InterruptMode::Coarse,
            InterruptMode::Off,
        ] {
            assert_eq!(InterruptMode::parse(mode.as_str()).unwrap(), mode);
        }
        assert_eq!(InterruptMode::default(), InterruptMode::Fine);
        Python::attach(|py| {
            let error = InterruptMode::parse("sometimes").unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
        });
    }

    #[test]
    fn uncancellable_calls_run_on_the_caller_without_a_scope() {
        Python::initialize();
        // A thread Python delivers no signals to, as a marimo app kernel.
        std::thread::spawn(|| {
            Python::attach(|py| {
                let caller = std::thread::current().id();
                let ran_here = run(py, Options::default(), Progress::none(), || {
                    // No cancellation scope: Vakint's checkpoints see nothing.
                    assert!(Cancellation::current().is_none());
                    Ok(std::thread::current().id() == caller)
                })
                .unwrap();
                assert!(ran_here);
                // A callback can cancel: then the call runs on the request
                // thread inside a scope.
                let callback = py
                    .eval(c"lambda snapshot: None", None, None)
                    .unwrap()
                    .unbind();
                let progress =
                    Progress::resolve(py, Some(&callback), "loading", "test", Duration::ZERO)
                        .unwrap();
                let scoped = run(py, Options::default(), progress, || {
                    Ok(Cancellation::current().is_some())
                })
                .unwrap();
                assert!(scoped || !LicenseManager::is_licensed());
                // "off" never installs a scope, even with a callback.
                let progress =
                    Progress::resolve(py, Some(&callback), "loading", "test", Duration::ZERO)
                        .unwrap();
                let options = Options {
                    mode: InterruptMode::Off,
                    ..Default::default()
                };
                let scoped = run(py, options, progress, || {
                    Ok(Cancellation::current().is_some())
                })
                .unwrap();
                assert!(!scoped);
            });
        })
        .join()
        .unwrap();
    }

    #[test]
    fn a_failing_final_snapshot_keeps_the_result() {
        Python::initialize();
        Python::attach(|py| {
            let locals = PyDict::new(py);
            py.run(
                c"
import sys
reported = []
previous_hook = sys.unraisablehook
sys.unraisablehook = reported.append
def callback(snapshot):
    if snapshot.status == 'succeeded':
        raise RuntimeError('display failed')
",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();
            let callback = locals.get_item("callback").unwrap().unwrap().unbind();
            let progress =
                Progress::resolve(py, Some(&callback), "loading", "final", Duration::ZERO).unwrap();
            let value = run(py, Options::default(), progress, || Ok(11));
            py.run(
                c"sys.unraisablehook = previous_hook",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();
            assert_eq!(value.unwrap(), 11);
            py.run(
                c"assert len(reported) == 1 and 'display failed' in str(reported[0].exc_value), reported",
                Some(&locals),
                Some(&locals),
            )
            .unwrap();
        });
    }

    #[test]
    fn request_threads_are_reused_and_a_panicking_one_is_retired() {
        // Alone in its process: other tests would share the pool.
        on_python_main_thread(
            "symbolica_community_module::interrupt::tests::request_threads_are_reused_and_a_panicking_one_is_retired",
            |_| pool_reuse_checks(),
        );
    }

    fn pool_reuse_checks() {
        let names = || -> Vec<String> {
            std::fs::read_dir("/proc/self/task")
                .map(|tasks| {
                    tasks
                        .filter_map(|task| {
                            std::fs::read_to_string(task.ok()?.path().join("comm")).ok()
                        })
                        .map(|name| name.trim().to_owned())
                        .collect()
                })
                .unwrap_or_default()
        };
        let request_threads = || {
            names()
                .iter()
                .filter(|name| *name == "vakint-request")
                .count()
        };
        let first = pool::run(|| std::thread::current().id(), || {})
            .unwrap()
            .unwrap();
        let second = pool::run(|| std::thread::current().id(), || {})
            .unwrap()
            .unwrap();
        assert_eq!(first, second, "an idle request thread is reused");
        assert_ne!(first, std::thread::current().id());
        let idle = request_threads();
        // Borrowed data is safe to use: the call waits for the work.
        let mut counted = 0;
        pool::run(|| counted += 1, || {}).unwrap().unwrap();
        assert_eq!(counted, 1);
        // A panic is returned, and its thread is not reused.
        let panicked: std::thread::Result<()> =
            pool::run(|| panic!("request panicked"), || {}).unwrap();
        assert!(panicked.is_err());
        let started = Instant::now();
        while request_threads() >= idle && started.elapsed() < Duration::from_secs(2) {
            std::thread::sleep(Duration::from_millis(10));
        }
        assert!(request_threads() < idle, "{} >= {idle}", request_threads());
        let third = pool::run(|| std::thread::current().id(), || {})
            .unwrap()
            .unwrap();
        assert_ne!(third, first);
        // Ticks run on this thread while the work runs.
        let mut ticks = 0;
        pool::run(
            || std::thread::sleep(Duration::from_millis(120)),
            || ticks += 1,
        )
        .unwrap()
        .unwrap();
        assert!(ticks >= 2, "{ticks}");
    }
}
