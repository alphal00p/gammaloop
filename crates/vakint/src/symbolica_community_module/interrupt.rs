//! One interruptible runner for the potentially long native calls.
//!
//! The native work runs on one scoped worker thread with the GIL released.
//! The calling thread checks Python signals every 25 ms (Ctrl-C, marimo's
//! Interrupt, `monitor="tui"`'s re-sent SIGINT) and delivers progress. An
//! interruption or a failing progress callback cancels the request
//! ([`crate::Cancellation`]): in [`Scope::Everywhere`] every Symbolica
//! resource checkpoint of the worker (and of the workers it dispatches) then
//! fails; in [`Scope::Coarse`] (data loads, whose dense Symbolica checkpoints
//! made a per-checkpoint budget cost 28% of a four-loop load) Vakint's own
//! checkpoints do (per file, program and pack). The worker is joined, and the
//! ORIGINAL exception is raised. Python handles signals on its main thread only; other callers
//! (a marimo app-mode kernel thread, a user thread) cannot be interrupted but
//! keep the interpreter responsive while the call runs.
//!
//! Restricted (unlicensed) Symbolica allows one Symbolica thread per user, and
//! the calling thread already is it. Then the work stays on the calling thread
//! (still without the GIL); a ticker marks a SIGINT check as due every 25 ms,
//! and the next checkpoint consumes a pending SIGINT and cancels. The signal is
//! re-delivered after the work returned, so Python's own handler raises
//! (KeyboardInterrupt, or marimo's interruption). No Python code runs inside a
//! checkpoint; live progress needs a license (the final snapshot is delivered).

use std::sync::{
    Arc, OnceLock,
    atomic::{AtomicBool, Ordering},
};
use std::time::Duration;

use pyo3::{PyErr, PyResult, Python, exceptions::PyRuntimeError};
use symbolica::license::LicenseManager;
use symbolica::utils::resources::ResourceBudget;

use super::{progress::Progress, vakint_to_python_error};
use crate::{Cancellation, ResourceLimits, VakintError};

/// Period of the calling thread's signal checks.
const SIGNAL_POLL: Duration = Duration::from_millis(25);

/// Where an interruption stops the native work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Scope {
    /// At every resource checkpoint (reductions, evaluations).
    Everywhere,
    /// At Vakint's own checkpoints only (data loads: per file, program and
    /// pack; a long native load step in between finishes first).
    Coarse,
}

/// Run `operation` under the engine's resource `limits`, interruptibly.
///
/// The operation must not touch Python objects: it runs while the calling
/// thread does not hold the GIL. Its error is converted on the calling thread.
pub(super) fn run<T: Send>(
    py: Python<'_>,
    limits: ResourceLimits,
    progress: Progress,
    operation: impl FnOnce() -> Result<T, VakintError> + Send,
) -> PyResult<T> {
    run_in(py, Scope::Everywhere, limits, progress, operation)
}

/// [`run`] for data loads: interruptible at Vakint's coarse checkpoints.
pub(super) fn run_load<T: Send>(
    py: Python<'_>,
    limits: ResourceLimits,
    progress: Progress,
    operation: impl FnOnce() -> Result<T, VakintError> + Send,
) -> PyResult<T> {
    run_in(py, Scope::Coarse, limits, progress, operation)
}

fn enter<R>(cancellation: &Cancellation, scope: Scope, operation: impl FnOnce() -> R) -> R {
    match scope {
        Scope::Everywhere => cancellation.run(operation),
        Scope::Coarse => cancellation.run_coarse(operation),
    }
}

fn run_in<T: Send>(
    py: Python<'_>,
    scope: Scope,
    limits: ResourceLimits,
    mut progress: Progress,
    operation: impl FnOnce() -> Result<T, VakintError> + Send,
) -> PyResult<T> {
    // A Ctrl-C pending before the call raises without starting any work.
    if let Err(error) = py.check_signals() {
        let _ = progress.finish(py, Some(&error), true);
        return Err(error);
    }
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
    let (outcome, interruption) = if symbolica_threads_allowed() {
        let cancellation = Cancellation::default();
        py.detach(|| {
            watch(&cancellation, &mut progress, || {
                enter(&cancellation, scope, request)
            })
        })
    } else {
        let (outcome, interrupted) = py.detach(|| on_caller(scope, request));
        let interruption = interrupted.then(|| redeliver_sigint(py));
        (Ok(Ok(outcome)), interruption)
    };
    let outcome = match outcome {
        Ok(Ok(result)) => result,
        Ok(Err(panic)) => std::panic::resume_unwind(panic),
        Err(error) => {
            return Err(PyRuntimeError::new_err(format!(
                "cannot start the Vakint request thread: {error}"
            )));
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
            progress.finish(py, None, false)?;
            Ok(value)
        }
        Err(error) => {
            let _ = progress.finish(py, Some(&error), false);
            Err(error)
        }
    }
}

type Joined<T> = std::thread::Result<Result<T, VakintError>>;

/// Whether threads other than the caller may use Symbolica (a license).
/// Licensing cannot change once Symbolica is in use, so it is read once.
fn symbolica_threads_allowed() -> bool {
    static ALLOWED: OnceLock<bool> = OnceLock::new();
    *ALLOWED.get_or_init(LicenseManager::is_licensed)
}

/// Run `request` on the calling thread (which must not hold the GIL),
/// cancelling at the first checkpoint after a SIGINT arrived. Returns the
/// result and whether a SIGINT cancelled it (the SIGINT is then consumed).
fn on_caller<T: Send>(
    scope: Scope,
    request: impl FnOnce() -> Result<T, VakintError> + Send,
) -> (Result<T, VakintError>, bool) {
    use std::sync::mpsc::{RecvTimeoutError, sync_channel};
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

/// Run `request` on a worker; poll signals and progress until it returns.
#[cfg(not(target_arch = "wasm32"))]
fn watch<T: Send>(
    cancellation: &Cancellation,
    progress: &mut Progress,
    request: impl FnOnce() -> Result<T, VakintError> + Send,
) -> (std::io::Result<Joined<T>>, Option<PyErr>) {
    use std::sync::mpsc::{RecvTimeoutError, sync_channel};
    let mut interruption = None;
    let joined = std::thread::scope(|scope| {
        // The worker drops `done` when it returns or unwinds.
        let (done, finished) = sync_channel::<()>(0);
        let worker = std::thread::Builder::new()
            .name("vakint-request".into())
            .spawn_scoped(scope, move || {
                let _done = done;
                request()
            })?;
        // Never sent: the worker returning (or unwinding) disconnects.
        while let Err(RecvTimeoutError::Timeout) = finished.recv_timeout(SIGNAL_POLL) {
            if interruption.is_none() {
                let polled = Python::attach(|py| {
                    py.check_signals()?;
                    progress.deliver(py)
                });
                if let Err(error) = polled {
                    interruption = Some(error);
                    cancellation.cancel();
                }
            }
        }
        Ok(worker.join())
    });
    (joined, interruption)
}

/// Browser kernels have no threads: run inline and poll signals only before
/// and after the request.
#[cfg(target_arch = "wasm32")]
fn watch<T: Send>(
    _cancellation: &Cancellation,
    _progress: &mut Progress,
    request: impl FnOnce() -> Result<T, VakintError> + Send,
) -> (std::io::Result<Joined<T>>, Option<PyErr>) {
    (Ok(Ok(request())), None)
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
                let result = run(py, ResourceLimits::default(), Progress::none(), || {
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
                let error = run_load(py, ResourceLimits::default(), Progress::none(), || {
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
                    ResourceLimits {
                        max_elapsed: Some(Duration::ZERO),
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
                    on_caller(Scope::Everywhere, || {
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
                    on_caller(Scope::Coarse, || {
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
                let error = run(py, ResourceLimits::default(), Progress::none(), || {
                    called.store(true, Ordering::Relaxed);
                    Ok(())
                })
                .unwrap_err();
                assert!(error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py));
                assert!(!called.load(Ordering::Relaxed));
                // Nothing pending afterwards: the next call runs.
                assert_eq!(
                    run(py, ResourceLimits::default(), Progress::none(), || Ok(7)).unwrap(),
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
            let error = run(py, ResourceLimits::default(), progress, || {
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
            assert_eq!(
                run(py, ResourceLimits::default(), progress, || Ok(3)).unwrap(),
                3
            );
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
                    run(
                        py,
                        ResourceLimits::default(),
                        resolve("loading")?,
                        || Ok(()),
                    )?;
                    exec("assert not marimo.status.spinner.called")?;
                    // A slow call shows a spinner after the delay, then exits it.
                    run(py, ResourceLimits::default(), resolve("loading")?, || {
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
                    run(py, ResourceLimits::default(), progress, || {
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
                    let error = run(py, ResourceLimits::default(), resolve("loading")?, || {
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
            run(py, ResourceLimits::default(), Progress::none(), || {
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
            let error = run::<()>(py, ResourceLimits::default(), Progress::none(), || {
                Err(VakintError::EvaluationError("native failure".into()))
            })
            .unwrap_err();
            assert!(error.is_instance_of::<pyo3::exceptions::PyValueError>(py));
            assert!(error.to_string().contains("Vakint error | "), "{error}");
        });
    }
}
