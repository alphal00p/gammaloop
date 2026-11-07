//! Live progress of long calls: `progress="auto"`, a callback, or `None`.
//!
//! The native monitor thread stores its latest model in a
//! [`RustRedProgressSink`]; the calling Python thread reads it between signal
//! checks and delivers [`VakintProgress`] snapshots to a marimo status
//! indicator or to a callback. Calls without a native monitor (data loads)
//! report their stage and elapsed time only.

use std::time::{Duration, Instant};

use feynkit_py::MarimoProgress;
use pyo3::{
    IntoPyObjectExt,
    exceptions::PyTypeError,
    prelude::*,
    types::{PyDict, PyList, PyString},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;

use crate::RustRedProgressSink;
use crate::rustred_evaluation::monitor::ProgressMetrics;

/// The `progress` keyword of the long calls: "auto" (the default), a
/// callable or None.
pub struct ProgressArg(Option<Py<PyAny>>);

impl ProgressArg {
    /// The default, "auto".
    pub(super) fn auto() -> Self {
        Self(Some(Python::attach(|py| {
            PyString::new(py, "auto").into_any().unbind()
        })))
    }

    pub(super) fn get(&self) -> Option<&Py<PyAny>> {
        self.0.as_ref()
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for ProgressArg {
    type Error = PyErr;

    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        Ok(Self((!value.is_none()).then(|| value.to_owned().unbind())))
    }
}

/// The stub generator renders the default through Python's repr.
impl<'py> IntoPyObject<'py> for ProgressArg {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = std::convert::Infallible;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        Ok(self.0.unwrap_or_else(|| py.None()).into_bound(py))
    }
}

#[cfg(feature = "python_stubgen")]
impl pyo3_stub_gen::PyStubType for ProgressArg {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo {
            name:
                "typing.Literal['auto'] | collections.abc.Callable[[VakintProgress], None] | None"
                    .into(),
            import: ["typing".into(), "collections.abc".into()]
                .into_iter()
                .collect(),
        }
    }
}

/// A marimo indicator appears only for calls that take longer than this.
const AUTO_DISPLAY_DELAY: Duration = Duration::from_millis(500);
/// Shortest period between two deliveries of the same stage.
const MIN_DELIVERY_PERIOD: Duration = Duration::from_millis(250);

/// One immutable progress snapshot of a running (or finished) Vakint call.
///
/// Passed to a ``progress`` callback of ``Vakint.evaluate`` and the other
/// long calls, on the calling thread, at most once per
/// ``progress_interval_seconds`` (at least 0.25 s apart), on every stage
/// change and once at the end with the final ``status``. Snapshots stay valid
/// after the call. Counts are absent (None) when the current stage has no
/// such measurement; ``phase`` is a position in the Ordered schedule, not a
/// fraction of the remaining time.
///
/// Examples
/// --------
/// >>> from symbolica import E
/// >>> from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod
/// >>> engine = Vakint(evaluation_order=[VakintEvaluationMethod.new_rustred_method(master_mode="retain")])
/// >>> snapshots = []
/// >>> integral = E("topo(I1L(1,2))", default_namespace="vakint")
/// >>> _ = engine.evaluate(integral, progress=snapshots.append)
/// >>> snapshots[-1].status
/// 'succeeded'
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VakintProgress",
    module = "symbolica.community.hepkit.vakint",
    frozen
)]
pub struct VakintProgress {
    metrics: ProgressMetrics,
    model: serde_json::Value,
}

impl super::ModuleInit for VakintProgress {}

impl VakintProgress {
    fn new(metrics: ProgressMetrics, model: serde_json::Value) -> Self {
        Self { metrics, model }
    }
}

fn to_python(py: Python<'_>, value: &serde_json::Value) -> PyResult<Py<PyAny>> {
    match value {
        serde_json::Value::Null => Ok(py.None()),
        serde_json::Value::Bool(value) => value.into_py_any(py),
        serde_json::Value::Number(number) => match (number.as_u64(), number.as_i64()) {
            (Some(value), _) => value.into_py_any(py),
            (None, Some(value)) => value.into_py_any(py),
            (None, None) => number.as_f64().unwrap_or(f64::NAN).into_py_any(py),
        },
        serde_json::Value::String(text) => text.into_py_any(py),
        serde_json::Value::Array(items) => {
            let list = PyList::empty(py);
            for item in items {
                list.append(to_python(py, item)?)?;
            }
            Ok(list.into_any().unbind())
        }
        serde_json::Value::Object(fields) => {
            let dict = PyDict::new(py);
            for (key, item) in fields {
                dict.set_item(key, to_python(py, item)?)?;
            }
            Ok(dict.into_any().unbind())
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl VakintProgress {
    /// "running", "succeeded", "failed", "partial" or "interrupted".
    #[getter]
    fn status(&self) -> &str {
        &self.metrics.status
    }

    /// Seconds since the call started.
    #[getter]
    fn elapsed_seconds(&self) -> f64 {
        self.metrics.elapsed_seconds
    }

    /// Stage kind: "setup", "data_load", "ingress", "lowering", "native",
    /// "post_processing", "output" or "other".
    #[getter]
    fn stage(&self) -> &str {
        &self.metrics.stage
    }

    /// What the call is doing now, in words (native activity or caller stage).
    #[getter]
    fn activity(&self) -> &str {
        &self.metrics.activity
    }

    /// Ordered batches of the evaluation: (completed, total, label).
    #[getter]
    fn batch(&self) -> Option<(usize, usize, String)> {
        self.metrics.batch.clone()
    }

    /// Ordered schedule position: (running phase index, phase total).
    #[getter]
    fn phase(&self) -> Option<(usize, usize)> {
        self.metrics.phase
    }

    /// Kind of the running Ordered phase: "route", "apply", "reduce" or "done".
    #[getter]
    fn phase_kind(&self) -> Option<&str> {
        self.metrics.phase_kind.as_deref()
    }

    /// Integral keys processed by the native reduction request.
    #[getter]
    fn processed_keys(&self) -> Option<u64> {
        self.metrics.processed_keys
    }

    /// Input keys of the native reduction request.
    #[getter]
    fn input_keys(&self) -> Option<u64> {
        self.metrics.input_keys
    }

    /// Recent (about 30 s) or whole-request key rate.
    #[getter]
    fn keys_per_second(&self) -> Option<f64> {
        self.metrics.keys_per_second
    }

    /// Keys on the Ordered frontier now.
    #[getter]
    fn frontier_keys(&self) -> Option<u64> {
        self.metrics.frontier_keys
    }

    /// Largest Ordered frontier of the request.
    #[getter]
    fn peak_frontier_keys(&self) -> Option<u64> {
        self.metrics.peak_frontier_keys
    }

    /// Ordered step-cache hits of the request.
    #[getter]
    fn step_cache_hits(&self) -> Option<u64> {
        self.metrics.step_cache_hits
    }

    /// Ordered step-cache misses (computed steps) of the request.
    #[getter]
    fn step_cache_misses(&self) -> Option<u64> {
        self.metrics.step_cache_misses
    }

    /// Resident memory of this process.
    #[getter]
    fn rss_bytes(&self) -> Option<u64> {
        self.metrics.rss_bytes
    }

    /// Peak resident memory of this process.
    #[getter]
    fn peak_rss_bytes(&self) -> Option<u64> {
        self.metrics.peak_rss_bytes
    }

    /// One line with the available metrics, e.g.
    /// "Apply 38/83 · 6312 keys/s · frontier 12.3k (peak 174.9k) · ...".
    #[getter]
    fn summary(&self) -> String {
        self.metrics.summary()
    }

    /// The full monitor model (schema ``vakint.monitor.v2``, as in the
    /// ``model`` object of a ``progress_snapshot`` file). Calls without a
    /// native monitor give the status, elapsed time and stage only. The final
    /// snapshot's ``"status"`` is the call's outcome: "Succeeded", "Failed"
    /// or "Interrupted" (the last one is not a ``vakint.monitor.v2`` status:
    /// the native monitor never sees the interruption).
    #[gen_stub(override_return_type(
        type_repr = "dict[str, typing.Any]",
        imports = ("typing")
    ))]
    fn to_dict(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        to_python(py, &self.model)
    }

    fn __repr__(&self) -> String {
        self.summary()
    }
}

enum Consumer {
    None,
    /// A marimo status indicator, created once the call is slow enough.
    Marimo(Option<MarimoProgress>),
    Callback(Py<PyAny>),
}

/// The progress consumer of one call and its delivery state.
pub(super) struct Progress {
    consumer: Consumer,
    sink: Option<RustRedProgressSink>,
    /// Stage of calls without a native monitor.
    stage_id: &'static str,
    activity: String,
    period: Duration,
    started: Instant,
    seen: u64,
    latest: Option<serde_json::Value>,
    /// A native model arrived that was not delivered yet.
    unshown: bool,
    delivered: Option<(Instant, String)>,
    /// The callback raised: never call it again.
    callback_failed: bool,
}

/// Whether an already imported marimo runs a notebook (never imports it).
fn marimo_notebook(py: Python<'_>) -> PyResult<bool> {
    let modules = py.import("sys")?.getattr("modules")?;
    let Some(marimo) = modules
        .cast::<PyDict>()?
        .get_item("marimo")?
        .filter(|module| !module.is_none())
    else {
        return Ok(false);
    };
    // An unexpected marimo (no such function) is not a notebook.
    Ok(marimo
        .call_method0("running_in_notebook")
        .and_then(|running| running.is_truthy())
        .unwrap_or(false))
}

impl Progress {
    /// No consumer (internal calls and the keyword's `None`).
    pub(super) fn none() -> Self {
        Self {
            consumer: Consumer::None,
            sink: None,
            stage_id: "",
            activity: String::new(),
            period: MIN_DELIVERY_PERIOD,
            started: Instant::now(),
            seen: 0,
            latest: None,
            unshown: false,
            delivered: None,
            callback_failed: false,
        }
    }

    /// Resolve the `progress` keyword: "auto" (a marimo indicator inside a
    /// notebook, else nothing), a callable or None. `stage_id` and
    /// `activity` describe calls without a native monitor.
    pub(super) fn resolve(
        py: Python<'_>,
        value: Option<&Py<PyAny>>,
        stage_id: &'static str,
        activity: impl Into<String>,
        interval: Duration,
    ) -> PyResult<Self> {
        let consumer = match value.map(|value| value.bind(py)) {
            None => Consumer::None,
            Some(value) if value.is_callable() => Consumer::Callback(value.clone().unbind()),
            Some(value) if value.extract::<String>().is_ok_and(|value| value == "auto") => {
                if marimo_notebook(py)? {
                    Consumer::Marimo(None)
                } else {
                    Consumer::None
                }
            }
            Some(_) => {
                return Err(PyTypeError::new_err(
                    "progress must be 'auto', a callable or None",
                ));
            }
        };
        Ok(Self {
            consumer,
            stage_id,
            activity: activity.into(),
            period: interval.max(MIN_DELIVERY_PERIOD),
            ..Self::none()
        })
    }

    /// Whether anything consumes progress (the native monitor is needed).
    pub(super) fn is_active(&self) -> bool {
        !matches!(self.consumer, Consumer::None)
    }

    /// Whether a callback consumes progress (its exception cancels the call).
    pub(super) fn is_callback(&self) -> bool {
        matches!(self.consumer, Consumer::Callback(_))
    }

    /// The call is being stopped: say so on the marimo indicator (a native
    /// step without checkpoints may take seconds to finish). Callbacks get
    /// the final "interrupted" snapshot only.
    pub(super) fn stopping(&mut self, py: Python<'_>) {
        let (completed, total) = self
            .latest
            .as_ref()
            .and_then(ProgressMetrics::from_model)
            .and_then(|metrics| metrics.phase)
            .map_or((0, None), |(index, total)| (index.min(total), Some(total)));
        if let Consumer::Marimo(Some(display)) = &mut self.consumer {
            let _ = display.update(
                py,
                "Vakint · stopping",
                "Interrupted: stopping at the next checkpoint (a native step without \
                 checkpoints finishes first)",
                completed,
                total,
            );
        }
    }

    /// The sink a monitored call publishes into, when progress is consumed.
    pub(super) fn native_sink(&mut self) -> Option<RustRedProgressSink> {
        if !self.is_active() {
            return None;
        }
        Some(self.sink.get_or_insert_with(Default::default).clone())
    }

    /// The current metrics and model; `None` until a native model arrived.
    /// Calls without a native monitor always have a new (caller) snapshot.
    fn current(&mut self) -> Option<(ProgressMetrics, serde_json::Value)> {
        match &self.sink {
            Some(sink) => {
                if let Some((sequence, model)) = sink.latest_since(self.seen) {
                    self.seen = sequence;
                    self.latest = Some(model);
                    self.unshown = true;
                }
                let model = self.latest.clone()?;
                let metrics = ProgressMetrics::from_model(&model)?;
                Some((metrics, model))
            }
            None => {
                self.unshown = true;
                let metrics = ProgressMetrics::caller(
                    self.stage_id,
                    &self.activity,
                    self.started.elapsed().as_secs_f64(),
                );
                let model = metrics.caller_model();
                Some((metrics, model))
            }
        }
    }

    fn show(
        &mut self,
        py: Python<'_>,
        metrics: ProgressMetrics,
        model: serde_json::Value,
    ) -> PyResult<()> {
        match &mut self.consumer {
            Consumer::None => Ok(()),
            Consumer::Callback(callback) => {
                if self.callback_failed {
                    return Ok(());
                }
                let snapshot = Py::new(py, VakintProgress::new(metrics, model))?;
                callback
                    .call1(py, (snapshot,))
                    .map(|_| ())
                    .inspect_err(|_| {
                        self.callback_failed = true;
                    })
            }
            Consumer::Marimo(display) => {
                let title = format!("Vakint · {}", metrics.activity);
                let subtitle = metrics.summary();
                if display.is_none() {
                    *display = MarimoProgress::detect(py, &title, &subtitle)?;
                }
                let Some(display) = display else {
                    return Ok(());
                };
                let (completed, total) = match metrics.phase {
                    Some((index, total)) => (index.min(total), Some(total)),
                    None => (0, None),
                };
                display.update(py, &title, &subtitle, completed, total)
            }
        }
    }

    /// Deliver the latest snapshot if it is new and due.
    pub(super) fn deliver(&mut self, py: Python<'_>) -> PyResult<()> {
        if !self.is_active() {
            return Ok(());
        }
        if matches!(self.consumer, Consumer::Marimo(None))
            && self.started.elapsed() < AUTO_DISPLAY_DELAY
        {
            return Ok(());
        }
        let Some((metrics, model)) = self.current() else {
            return Ok(());
        };
        let key = format!(
            "{}/{}/{}",
            metrics.stage,
            metrics.stage_id,
            metrics.phase_kind.as_deref().unwrap_or_default()
        );
        let due = match &self.delivered {
            None => true,
            Some((at, last)) => *last != key || at.elapsed() >= self.period,
        };
        if !self.unshown || !due {
            return Ok(());
        }
        self.unshown = false;
        self.delivered = Some((Instant::now(), key));
        let shown = self.show(py, metrics, model);
        match shown {
            // A failing marimo indicator is a display problem, not a reason to
            // stop the call: report it and stop showing (a KeyboardInterrupt
            // raised while its Python code ran still propagates).
            Err(error)
                if matches!(self.consumer, Consumer::Marimo(_))
                    && error.is_instance_of::<pyo3::exceptions::PyException>(py) =>
            {
                error.write_unraisable(py, None);
                if let Consumer::Marimo(Some(display)) = &mut self.consumer {
                    let _ = display.finish(py, None, "Vakint");
                }
                self.consumer = Consumer::None;
                Ok(())
            }
            shown => shown,
        }
    }

    /// Final snapshot and display cleanup. `error` is the call's error (the
    /// original exception is never replaced); `interrupted` marks a signal or
    /// callback interruption.
    pub(super) fn finish(
        &mut self,
        py: Python<'_>,
        error: Option<&PyErr>,
        interrupted: bool,
    ) -> PyResult<()> {
        let status = match (error, interrupted) {
            (None, _) => "succeeded",
            (Some(_), true) => "interrupted",
            (Some(_), false) => "failed",
        };
        let delivered = if self.is_active() {
            // The monitor thread has finished: its final model is stored.
            self.seen = 0;
            let (mut metrics, mut model) = match self.current() {
                Some((metrics, model)) => (metrics, model),
                None => {
                    let metrics = ProgressMetrics::caller(
                        self.stage_id,
                        &self.activity,
                        self.started.elapsed().as_secs_f64(),
                    );
                    let model = metrics.caller_model();
                    (metrics, model)
                }
            };
            metrics.status = status.into();
            if let Some(object) = model.as_object_mut() {
                let mut label = status.to_owned();
                label[..1].make_ascii_uppercase();
                object.insert("status".into(), label.into());
            }
            if matches!(self.consumer, Consumer::Callback(_)) {
                self.show(py, metrics, model)
            } else {
                Ok(())
            }
        } else {
            Ok(())
        };
        let cleanup = match &mut self.consumer {
            Consumer::Marimo(Some(display)) => display.finish(py, error, "Vakint stopped"),
            _ => Ok(()),
        };
        self.consumer = Consumer::None;
        if error.is_some() {
            // Cleanup never replaces the original error.
            return Ok(());
        }
        delivered.and(cleanup)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keyword_resolution_rejects_unknown_values_and_auto_needs_a_notebook() {
        Python::initialize();
        Python::attach(|py| {
            let auto = PyString::new(py, "auto").into_any().unbind();
            // Rust test binaries never import marimo: "auto" is exactly None.
            let progress = Progress::resolve(py, Some(&auto), "", "", Duration::ZERO).unwrap();
            assert!(!progress.is_active());
            let none = Progress::resolve(py, None, "", "", Duration::ZERO).unwrap();
            assert!(!none.is_active());
            let bad = PyString::new(py, "tui").into_any().unbind();
            let error = Progress::resolve(py, Some(&bad), "", "", Duration::ZERO)
                .err()
                .unwrap();
            assert!(error.is_instance_of::<PyTypeError>(py));
            let number = 3_i64.into_pyobject(py).unwrap().into_any().unbind();
            assert!(Progress::resolve(py, Some(&number), "", "", Duration::ZERO).is_err());
            let callback = py
                .eval(c"lambda snapshot: None", None, None)
                .unwrap()
                .unbind();
            let mut progress =
                Progress::resolve(py, Some(&callback), "", "", Duration::from_secs(2)).unwrap();
            assert!(progress.is_active());
            assert_eq!(progress.period, Duration::from_secs(2));
            assert!(progress.native_sink().is_some());
        });
    }

    #[test]
    fn snapshots_convert_the_model_to_python_values() {
        Python::initialize();
        Python::attach(|py| {
            let metrics = ProgressMetrics::caller("external_data_loading", "loading", 1.5);
            let model = serde_json::json!({
                "status": "Running", "elapsed_seconds": 1.5, "stages": [1, -2, 0.5, null, true],
                "stage": {"kind": "data_load"},
            });
            let snapshot = Py::new(py, VakintProgress::new(metrics, model)).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("p", snapshot).unwrap();
            py.run(
                cr#"
d = p.to_dict()
assert d["stages"] == [1, -2, 0.5, None, True], d
assert d["stage"]["kind"] == "data_load"
assert p.status == "running" and p.stage == "data_load" and p.activity == "loading"
assert p.phase is None and p.frontier_keys is None and p.batch is None
assert repr(p) == p.summary == "loading · 1.5s", repr(p)
"#,
                None,
                Some(&locals),
            )
            .unwrap();
        });
    }
}
