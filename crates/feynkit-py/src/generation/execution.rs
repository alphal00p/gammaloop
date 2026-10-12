use super::{PyGenerationProgress, PyParticle, PyRequest, PyRun, error};
use crate::progress::MarimoProgress;
use feynkit_generator::{
    ExternalLabel, GenerationControl, GenerationError, GenerationProgress, Request, Run, Stage,
};
use feynkit_model::{Model, ParticleId};
#[cfg(not(target_arch = "wasm32"))]
use linnet::half_edge::generate::batch::Hook;
use linnet::half_edge::generate::batch::{Reply, Request as BatchRequest};
use linnet_py::{CancelCheck, GenerationCallbacks};
use pyo3::{exceptions::PyTypeError, prelude::*, types::PyTuple};
use std::{
    collections::VecDeque,
    ops::ControlFlow,
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
    },
    time::{Duration, Instant},
};

/// Run generation with callbacks and signals on the calling Python thread.
pub(super) fn execute<T: Send>(
    py: Python<'_>,
    model: Arc<Model>,
    hooks: Option<&PyRequest>,
    controls: Option<&PyRun>,
    work: impl FnOnce(Option<Request>, Run) -> Result<T, GenerationError> + Send,
) -> PyResult<T> {
    let mut request = hooks.map(|v| v.inner.clone());
    #[cfg(not(target_arch = "wasm32"))]
    let thread_safe = controls.is_some_and(|v| v.thread_safe);
    #[cfg(not(target_arch = "wasm32"))]
    let worker_hooks = hooks.map(|v| v.worker_hooks).unwrap_or_default();
    let mut run = controls.map_or_else(Run::default, |v| v.inner.clone());
    if hooks.is_some_and(|v| v.check_heredity) {
        run = run.check_heredity(true);
    }
    let progress = controls.and_then(|v| v.progress.as_ref().map(|p| p.clone_ref(py)));
    let prune_skeletons = hooks.and_then(|v| v.prune_skeletons.as_ref().map(|p| p.clone_ref(py)));
    let prune_dressings = hooks.and_then(|v| v.prune_dressings.as_ref().map(|p| p.clone_ref(py)));
    let select = hooks.and_then(|v| v.select.as_ref().map(|p| p.clone_ref(py)));
    py.check_signals()?;
    let automatic = progress.as_ref().is_some_and(|value| {
        value
            .bind(py)
            .extract::<String>()
            .is_ok_and(|s| s == "auto")
    });
    let progress = if automatic { None } else { progress };
    if progress
        .as_ref()
        .is_some_and(|callback| !callback.bind(py).is_callable())
    {
        return Err(PyTypeError::new_err(
            "progress must be 'auto', callable, or None",
        ));
    }
    let callbacks = GenerationCallbacks::new(py, prune_skeletons, prune_dressings, select)?;
    let mask = callbacks.mask();
    let has_hooks = mask.skeletons || mask.dressings || mask.select;
    #[cfg(not(target_arch = "wasm32"))]
    let on_worker = move |hook| {
        thread_safe
            || match hook {
                Hook::Skeleton => worker_hooks.skeletons,
                Hook::Dressing => worker_hooks.dressings,
                Hook::Select => worker_hooks.select,
            }
    };
    #[cfg(not(target_arch = "wasm32"))]
    let has_worker_hooks = (mask.skeletons && on_worker(Hook::Skeleton))
        || (mask.dressings && on_worker(Hook::Dressing))
        || (mask.select && on_worker(Hook::Select));
    let cancelled = Arc::new(AtomicBool::new(false));
    let check = cancelled.clone();
    run = run.cancellation_check(move || check.load(Ordering::Relaxed));
    let marimo_progress = if automatic {
        MarimoProgress::detect(py, "Generating diagrams", "Starting generation")?
            .map(|display| Arc::new(Mutex::new(display)))
    } else {
        None
    };
    let display = marimo_progress.clone();
    let started = Instant::now();
    let interruption = Arc::new(Mutex::new(None));
    let snapshots = Arc::new(Mutex::new(VecDeque::<GenerationProgress>::new()));
    let delivery = Mutex::new((Instant::now(), None));
    let pending = snapshots.clone();
    let has_progress = progress.is_some() || display.is_some();
    let progress_token = run.cancellation().clone();
    let deliver_progress = Arc::new(move |force: bool| -> PyResult<()> {
        let mut pending = pending.lock().unwrap();
        let mut delivery = delivery.lock().unwrap();
        let changed = pending.front().is_some_and(|p| Some(p.stage) != delivery.1);
        let finished = pending
            .back()
            .is_some_and(|p| p.total.as_ref() == Some(&p.completed));
        if !force
            && !changed
            && !finished
            && pending.len() < 2
            && delivery.0.elapsed() < Duration::from_millis(150)
        {
            return Ok(());
        }
        let updates = std::mem::take(&mut *pending);
        if let Some(last) = updates.back() {
            *delivery = (Instant::now(), Some(last.stage));
        }
        drop(pending);
        drop(delivery);
        if has_progress {
            Python::attach(|py| -> PyResult<()> {
                for inner in updates {
                    if let Some(display) = &display {
                        let title = match inner.stage {
                            Stage::Skeletons => "Enumerating skeletons",
                            Stage::Dressing => "Assigning interactions",
                            Stage::Prepare => "Preparing physical diagrams",
                            Stage::Numerators => "Constructing numerators",
                            Stage::Selection => "Selecting diagrams",
                            Stage::ZeroColor => "Filtering zero-color diagrams",
                            Stage::GroupingPreparation => "Preparing numerators for grouping",
                            Stage::GroupingSamples => "Sampling numerators for grouping",
                            Stage::GroupingComparison => "Comparing numerators",
                            Stage::Grouping => "Finalizing diagram groups",
                            Stage::Done if progress_token.is_cancelled() => {
                                "Diagram generation cancelled"
                            }
                            Stage::Done => "Diagram generation complete",
                        };
                        let counts = match (inner.stage, &inner.total) {
                            (Stage::Done, _) => format!("{} diagrams retained", inner.completed),
                            (_, Some(total)) => format!("{} / {total} processed", inner.completed),
                            (_, None) => format!("{} processed", inner.completed),
                        };
                        let counts = if inner.zero_numerator_count.bits() > 0 {
                            format!(
                                "{counts} · {} zero graphs filtered",
                                inner.zero_numerator_count
                            )
                        } else {
                            counts
                        };
                        let subtitle =
                            format!("{counts} · {:.1}s elapsed", started.elapsed().as_secs_f64());
                        display.lock().unwrap().update(
                            py,
                            title,
                            &subtitle,
                            (&inner.completed).try_into().unwrap_or(usize::MAX),
                            inner
                                .total
                                .as_ref()
                                .map(|n| n.try_into().unwrap_or(usize::MAX)),
                        )?;
                    }
                    if let Some(callback) = &progress {
                        let reply = callback.call1(py, (PyGenerationProgress { inner },))?;
                        if !reply.is_none(py)
                            && !reply.extract::<bool>(py).map_err(|_| {
                                PyTypeError::new_err("progress must return None or bool")
                            })?
                        {
                            progress_token.cancel();
                            break;
                        }
                    }
                }
                Ok(())
            })?;
        }
        Ok(())
    });
    #[cfg(not(target_arch = "wasm32"))]
    let (progress_notifications, progress_receiver) =
        std::sync::mpsc::sync_channel::<std::sync::mpsc::SyncSender<GenerationControl>>(1);
    #[cfg(not(target_arch = "wasm32"))]
    let progress_cancellation = run.cancellation().clone();
    if has_progress {
        #[cfg(target_arch = "wasm32")]
        let (deliver, error) = (deliver_progress.clone(), interruption.clone());
        run = run.progress(move |snapshot| {
            // These owner boundaries must reach Python before native completion is
            // finalized or construction/selection starts. Other updates coalesce.
            let synchronize = snapshot.stage == Stage::Done
                || (matches!(snapshot.stage, Stage::Numerators | Stage::Selection)
                    && snapshot.completed.bits() == 0);
            let mut pending = snapshots.lock().unwrap();
            if pending
                .back()
                .is_some_and(|last| last.stage == snapshot.stage)
            {
                *pending.back_mut().unwrap() = snapshot;
            } else {
                pending.push_back(snapshot);
            }
            drop(pending);
            #[cfg(not(target_arch = "wasm32"))]
            if synchronize {
                let (reply, response) = std::sync::mpsc::sync_channel(1);
                if progress_notifications.send(reply).is_err() {
                    return GenerationControl::Cancel;
                }
                return response.recv().unwrap_or(GenerationControl::Cancel);
            }
            #[cfg(target_arch = "wasm32")]
            if error.lock().unwrap().is_none() {
                let result = deliver(synchronize);
                if let Err(exception) = result {
                    *error.lock().unwrap() = Some(exception);
                    return GenerationControl::Cancel;
                }
            }
            GenerationControl::Continue
        });
    }
    // The key table belongs to one compiled plan. Processes with alternative external
    // states may use several plans, so refresh it only when their sorted keys differ.
    let keys = Mutex::new(
        None::<(
            Vec<ParticleId>,
            Vec<ExternalLabel>,
            Py<PyTuple>,
            Py<PyTuple>,
        )>,
    );
    let callback_model = Arc::clone(&model);
    let token = run.cancellation().clone();
    let callback_flag = Arc::clone(&cancelled);
    let callback_cancelled: CancelCheck =
        Arc::new(move || token.is_cancelled() || callback_flag.load(Ordering::Relaxed));
    let apply_batch =
        move |requests: Vec<BatchRequest<ParticleId, ExternalLabel>>| -> PyResult<Option<Vec<Reply>>> {
            Python::attach(|py| {
                let Some(first) = requests.first() else {
                    return Ok(Some(Vec::new()));
                };
                let mut keys = keys.lock().unwrap();
                if keys.as_ref().is_none_or(|(fields, labels, ..)| {
                    fields != first.fields() || labels != first.labels()
                }) {
                    let fields = first
                        .fields()
                        .iter()
                        .map(|&id| PyParticle::new(id, Arc::clone(&callback_model)));
                    let labels = first.labels().iter().map(|label| match *label {
                        ExternalLabel::Exact { state, index } => (state.as_str(), Some(index)),
                        ExternalLabel::WithinState(state) => (state.as_str(), None),
                    });
                    *keys = Some((
                        first.fields().to_vec(),
                        first.labels().to_vec(),
                        PyTuple::new(py, fields)?.unbind(),
                        PyTuple::new(py, labels)?.unbind(),
                    ));
                }
                let (_, _, fields, labels) = keys.as_ref().unwrap();
                callbacks.dispatch(
                    py,
                    requests.into_iter().map(BatchRequest::into_ranked).collect(),
                    fields.bind(py),
                    labels.bind(py),
                    &callback_cancelled,
                )
            })
        };
    #[cfg(not(target_arch = "wasm32"))]
    let result = {
        let (requests, receiver) = std::sync::mpsc::sync_channel(1);
        if has_hooks {
            request = request.map(|request| {
                request.deferred_hooks(mask, move |batch| {
                    let (reply, response) = std::sync::mpsc::sync_channel(1);
                    if requests.send((batch, reply)).is_err() {
                        return ControlFlow::Break(());
                    }
                    response.recv().unwrap_or(ControlFlow::Break(()))
                })
            });
        }
        let interruption = interruption.clone();
        let deliver_progress = deliver_progress.clone();
        py.detach(move || {
            std::thread::scope(|scope| {
                // At most one callback worker belongs to this execution. The native
                // Request retains only channel handles; no Python references escape in a set.
                type Work = (
                    Vec<BatchRequest<ParticleId, ExternalLabel>>,
                    std::sync::mpsc::SyncSender<PyResult<Option<Vec<Reply>>>>,
                );
                let (callback_requests, callback_receiver) =
                    std::sync::mpsc::sync_channel::<Work>(1);
                let callback_worker = has_worker_hooks.then(|| {
                    let apply_batch = &apply_batch;
                    scope.spawn(move || {
                        while let Ok((batch, reply)) = callback_receiver.recv() {
                            let _ = reply.send(apply_batch(batch));
                        }
                    })
                });
                let worker = scope.spawn(|| work(request, run));
                // Native workers finish their wave before the generation coordinator
                // sends a batch here. Only that coordinator awaits a combined reply.
                let service_caller = || {
                    let progress_reply = progress_receiver.try_recv().ok();
                    if interruption.lock().unwrap().is_none() {
                        let result = Python::attach(|py| {
                            let result = py
                                .check_signals()
                                .and_then(|()| deliver_progress(progress_reply.is_some()));
                            // Publish cancellation before releasing the GIL to the predicate
                            // worker, so it cannot start another callback after a signal.
                            if result.is_err() {
                                cancelled.store(true, Ordering::Relaxed);
                            }
                            result
                        });
                        if let Err(error) = result {
                            *interruption.lock().unwrap() = Some(error);
                        }
                    }
                    if let Some(reply) = progress_reply {
                        let control = if interruption.lock().unwrap().is_some()
                            || progress_cancellation.is_cancelled()
                        {
                            GenerationControl::Cancel
                        } else {
                            GenerationControl::Continue
                        };
                        let _ = reply.send(control);
                    }
                };
                while !worker.is_finished() {
                    service_caller();
                    if let Ok((batch, reply)) = receiver.recv_timeout(Duration::from_millis(25)) {
                        let response = if interruption.lock().unwrap().is_some() {
                            ControlFlow::Break(())
                        } else {
                            let applied = (|| -> PyResult<Option<Vec<Reply>>> {
                                // Preserve native order while keeping each unmarked hook on
                                // the caller. Both routes use the same scoped batch dispatcher.
                                let mut replies = Vec::with_capacity(batch.len());
                                let mut pending = batch.into_iter().peekable();
                                while let Some(first) = pending.next() {
                                    let use_worker = on_worker(first.hook());
                                    let mut batch = vec![first];
                                    while let Some(next) = pending
                                        .next_if(|next| on_worker(next.hook()) == use_worker)
                                    {
                                        batch.push(next);
                                    }
                                    let reply = if use_worker {
                                        let (done, response) = std::sync::mpsc::sync_channel(1);
                                        if callback_requests.send((batch, done)).is_err() {
                                            Err(pyo3::exceptions::PyRuntimeError::new_err(
                                                "callback worker stopped",
                                            ))
                                        } else {
                                            loop {
                                                match response.recv_timeout(Duration::from_millis(25)) {
                                                    Ok(result) => break result,
                                                    Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                                                        service_caller()
                                                    }
                                                    Err(
                                                        std::sync::mpsc::RecvTimeoutError::Disconnected,
                                                    ) => {
                                                        break Err(
                                                            pyo3::exceptions::PyRuntimeError::new_err(
                                                                "callback worker stopped",
                                                            ),
                                                        );
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        apply_batch(batch)
                                    }?;
                                    let Some(mut reply) = reply else {
                                        return Ok(None);
                                    };
                                    replies.append(&mut reply);
                                }
                                Ok(Some(replies))
                            })();
                            match applied {
                                Ok(Some(replies)) => ControlFlow::Continue(replies),
                                Ok(None) => ControlFlow::Break(()),
                                Err(error) => {
                                    interruption.lock().unwrap().get_or_insert(error);
                                    cancelled.store(true, Ordering::Relaxed);
                                    ControlFlow::Break(())
                                }
                            }
                        };
                        let _ = reply.send(response);
                    }
                }
                drop(callback_requests);
                if let Some(callback_worker) = callback_worker {
                    callback_worker
                        .join()
                        .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
                }
                // Native heredity diagnostics panic; keep that failure inside this binding.
                worker.join()
            })
        })
    };
    #[cfg(target_arch = "wasm32")]
    let result = {
        // Browser kernels service the same joined batch directly on their calling thread.
        let signal_error = interruption.clone();
        let signal_cancelled = Arc::clone(&cancelled);
        run = run.cancellation_check(move || {
            if signal_error.lock().unwrap().is_none() {
                if let Err(error) = Python::attach(|py| py.check_signals()) {
                    *signal_error.lock().unwrap() = Some(error);
                    signal_cancelled.store(true, Ordering::Relaxed);
                }
            }
            signal_cancelled.load(Ordering::Relaxed)
        });
        if has_hooks {
            let error = interruption.clone();
            request = request.map(|request| {
                request.deferred_hooks(mask, move |batch| match apply_batch(batch) {
                    Ok(Some(replies)) => ControlFlow::Continue(replies),
                    Ok(None) => ControlFlow::Break(()),
                    Err(exception) => {
                        *error.lock().unwrap() = Some(exception);
                        cancelled.store(true, Ordering::Relaxed);
                        ControlFlow::Break(())
                    }
                })
            });
        }
        py.detach(|| std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(request, run))))
    };
    let result = (|| {
        if let Some(error) = interruption.lock().unwrap().take() {
            return Err(error);
        }
        deliver_progress(true)?;
        py.check_signals()?;
        result
            .map_err(|panic| {
                let message = panic
                    .downcast_ref::<String>()
                    .map(String::as_str)
                    .or_else(|| panic.downcast_ref::<&str>().copied())
                    .unwrap_or("generation panicked");
                error::GenerationError::new_err(message.to_owned())
            })?
            .map_err(error::generation)
    })();
    if let Some(display) = marimo_progress {
        display
            .lock()
            .unwrap()
            .finish(py, result.as_ref().err(), "Diagram generation stopped")?;
    }
    result
}
