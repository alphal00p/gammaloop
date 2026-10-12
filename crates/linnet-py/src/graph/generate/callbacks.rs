//! Scoped, batched Python dispatch. Every predicate receives one typed view, validates
//! its bool immediately, and expires that view in a finally block.

use std::sync::{atomic::AtomicBool, Arc};

use linnet::half_edge::generate::batch::{Hook, HookMask, Reply, Request};
use pyo3::{
    class::gc::{PyTraverseError, PyVisit},
    exceptions::PyRuntimeError,
    prelude::*,
    types::PyTuple,
};

use pyo3::sync::PyOnceLock;

static DISPATCH: PyOnceLock<Py<PyAny>> = PyOnceLock::new();

/// Native cancellation controls, checked between individual predicates and after a batch.
#[doc(hidden)]
pub type CancelCheck = Arc<dyn Fn() -> bool + Send + Sync>;

use super::{
    cleared,
    views::{PyDiagramView, PyDressingView, PySkeletonView},
    Objects,
};

/// Exactly one Python owner traverses the actual callable references. Plans and executions own
/// ordinary Py references to this holder; no native Arc hides a Python reference from GC.
#[pyclass]
pub(super) struct Hooks {
    skeletons: Option<Py<PyAny>>,
    dressings: Option<Py<PyAny>>,
    select: Option<Py<PyAny>>,
}
impl Hooks {
    pub(super) fn new(
        py: Python<'_>,
        skeletons: Option<Py<PyAny>>,
        dressings: Option<Py<PyAny>>,
        select: Option<Py<PyAny>>,
    ) -> PyResult<Option<Py<Self>>> {
        for callback in [&skeletons, &dressings, &select].into_iter().flatten() {
            if !callback.bind(py).is_callable() {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "generation hooks must be callable",
                ));
            }
        }
        if skeletons.is_none() && dressings.is_none() && select.is_none() {
            return Ok(None);
        }
        Py::new(
            py,
            Self {
                skeletons,
                dressings,
                select,
            },
        )
        .map(Some)
    }
    pub(super) fn mask(&self) -> HookMask {
        HookMask {
            skeletons: self.skeletons.is_some(),
            dressings: self.dressings.is_some(),
            select: self.select.is_some(),
        }
    }
    fn callback(&self, py: Python<'_>, hook: Hook) -> PyResult<Py<PyAny>> {
        match hook {
            Hook::Skeleton => &self.skeletons,
            Hook::Dressing => &self.dressings,
            Hook::Select => &self.select,
        }
        .as_ref()
        .map(|c| c.clone_ref(py))
        .ok_or_else(cleared)
    }
}
#[pymethods]
impl Hooks {
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.skeletons)?;
        visit.call(&self.dressings)?;
        visit.call(&self.select)
    }
    fn __clear__(&mut self) {
        self.skeletons = None;
        self.dressings = None;
        self.select = None;
    }
}

enum View {
    Skeleton(Py<PySkeletonView>),
    Dressing(Py<PyDressingView>),
    Diagram(Py<PyDiagramView>),
}
impl View {
    fn activate(&self, py: Python<'_>) -> PyResult<()> {
        match self {
            Self::Skeleton(v) => v.borrow(py).activate(),
            Self::Dressing(v) => v.borrow(py).activate(),
            Self::Diagram(v) => v.borrow(py).activate(),
        }
    }
    fn object(&self, py: Python<'_>) -> Py<PyAny> {
        match self {
            Self::Skeleton(v) => v.clone_ref(py).into_any(),
            Self::Dressing(v) => v.clone_ref(py).into_any(),
            Self::Diagram(v) => v.clone_ref(py).into_any(),
        }
    }
    fn expire(&self, py: Python<'_>) {
        match self {
            Self::Skeleton(v) => v.borrow_mut(py).expire(),
            Self::Dressing(v) => v.borrow_mut(py).expire(),
            Self::Diagram(v) => v.borrow_mut(py).expire(),
        }
    }
    fn traverse(&self, visit: &PyVisit<'_>) -> Result<(), PyTraverseError> {
        match self {
            Self::Skeleton(v) => visit.call(v),
            Self::Dressing(v) => visit.call(v),
            Self::Diagram(v) => visit.call(v),
        }
    }
}

// Exception tracebacks can retain this owned request after its callback worker exits.
// Like the typed views, its native snapshots and Python handles are safe to release
// on the thread that later collects that traceback. PyO3 enforces its normal
// thread-safety bounds.
#[pyclass]
struct Pending {
    request: Option<Arc<Request<u32, u32>>>,
    view: Option<View>,
    callback: Option<Py<PyAny>>,
    reply: Option<Reply>,
    cancelled: Option<CancelCheck>,
}
impl Pending {
    fn new(
        py: Python<'_>,
        request: Request<u32, u32>,
        objects: &Objects,
        hooks: &Hooks,
        cancelled: CancelCheck,
    ) -> PyResult<Self> {
        let callback = hooks.callback(py, request.hook())?;
        let request = Arc::new(request);
        let valid = Arc::new(AtomicBool::new(false));
        let view = match request.hook() {
            Hook::Skeleton => View::Skeleton(Py::new(
                py,
                PySkeletonView::from_request(Arc::clone(&request), valid, objects.clone_ref(py)),
            )?),
            Hook::Dressing => View::Dressing(Py::new(
                py,
                PyDressingView::from_request(Arc::clone(&request), valid, objects.clone_ref(py)),
            )?),
            Hook::Select => View::Diagram(Py::new(
                py,
                PyDiagramView::from_request(Arc::clone(&request), valid, objects.clone_ref(py)),
            )?),
        };
        Ok(Self {
            request: Some(request),
            view: Some(view),
            callback: Some(callback),
            reply: None,
            cancelled: Some(cancelled),
        })
    }
}
#[pymethods]
impl Pending {
    fn activate(&self, py: Python<'_>) -> PyResult<()> {
        self.view.as_ref().ok_or_else(cleared)?.activate(py)
    }
    #[getter]
    fn cancelled(&self) -> PyResult<bool> {
        self.cancelled
            .as_ref()
            .map(|check| check())
            .ok_or_else(cleared)
    }
    #[getter]
    fn callback(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.callback
            .as_ref()
            .map(|c| c.clone_ref(py))
            .ok_or_else(cleared)
    }
    #[getter]
    fn view(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.view.as_ref().map(|v| v.object(py)).ok_or_else(cleared)
    }
    fn resolve(&mut self, result: &Bound<'_, PyAny>) -> PyResult<()> {
        if self.reply.is_some() {
            return Err(PyRuntimeError::new_err(
                "generation verdict was already supplied",
            ));
        }
        // Preserve PyO3's bool extraction (including supported numpy bools), not Python truthiness.
        let keep = result.extract::<bool>()?;
        self.reply = Some(self.request.as_ref().ok_or_else(cleared)?.reply(keep));
        Ok(())
    }
    fn expire(&mut self, py: Python<'_>) {
        if let Some(view) = self.view.take() {
            view.expire(py);
        }
        self.request = None;
        self.callback = None;
        self.cancelled = None;
    }
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.callback)?;
        if let Some(view) = &self.view {
            view.traverse(&visit)?;
        }
        Ok(())
    }
    fn __clear__(&mut self) {
        // GC only clears unreachable Pending objects; ordinary and exceptional dispatch use
        // expire() while the GIL is held so retained views lose their request snapshots too.
        self.request = None;
        self.view = None;
        self.callback = None;
        self.reply = None;
        self.cancelled = None;
    }
}

pub(super) fn dispatch(
    py: Python<'_>,
    requests: Vec<Request<u32, u32>>,
    objects: &Objects,
    hooks: &Py<Hooks>,
    cancelled: &CancelCheck,
) -> PyResult<Option<Vec<Reply>>> {
    let hooks = hooks.borrow(py);
    let pending = requests
        .into_iter()
        .map(|request| {
            Py::new(
                py,
                Pending::new(py, request, objects, &hooks, Arc::clone(cancelled))?,
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    let batch = PyTuple::new(py, pending.iter().map(|p| p.clone_ref(py)))?;
    // Python invokes predicates directly; native result extraction occurs before the next.
    let dispatcher = DISPATCH.get_or_try_init(py, || {
        let globals = pyo3::types::PyDict::new(py);
        py.run(c"def dispatch(batch):\n    for request in batch:\n        if request.cancelled:\n            break\n        try:\n            request.activate()\n            request.resolve(request.callback(request.view))\n        finally:\n            request.expire()\n", Some(&globals), None)?;
        Ok::<_,PyErr>(globals.get_item("dispatch")?.expect("defined dispatcher").unbind())
    })?;
    let result = dispatcher.call1(py, (batch,));
    // Also invalidate requests not visited after a failing predicate; no further callbacks run.
    for request in &pending {
        request.borrow_mut(py).expire(py);
    }
    result?;
    if cancelled() {
        return Ok(None);
    }
    pending
        .into_iter()
        .map(|request| {
            request
                .borrow_mut(py)
                .reply
                .take()
                .ok_or_else(|| PyRuntimeError::new_err("missing generation verdict"))
        })
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

/// Scoped transport shared with eager physics generation. Keep this Rust owner local
/// to the invocation; its ordinary Py reference owns the GC-visible callable holder.
#[doc(hidden)]
pub struct GenerationCallbacks {
    hooks: Py<Hooks>,
    mask: HookMask,
    calling_thread: std::thread::ThreadId,
}

impl GenerationCallbacks {
    pub fn new(
        py: Python<'_>,
        skeletons: Option<Py<PyAny>>,
        dressings: Option<Py<PyAny>>,
        select: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        let hooks = match Hooks::new(py, skeletons, dressings, select)? {
            Some(hooks) => hooks,
            None => Py::new(
                py,
                Hooks {
                    skeletons: None,
                    dressings: None,
                    select: None,
                },
            )?,
        };
        let mask = hooks.borrow(py).mask();
        Ok(Self {
            hooks,
            mask,
            calling_thread: std::thread::current().id(),
        })
    }

    pub fn mask(&self) -> HookMask {
        self.mask
    }

    pub fn dispatch(
        &self,
        py: Python<'_>,
        requests: Vec<Request<u32, u32>>,
        fields: &Bound<'_, PyTuple>,
        labels: &Bound<'_, PyTuple>,
        cancelled: &CancelCheck,
    ) -> PyResult<Option<Vec<Reply>>> {
        let objects = Objects {
            fields: fields.clone().unbind(),
            labels: labels.clone().unbind(),
            graph_thread: Some(self.calling_thread),
        };
        dispatch(py, requests, &objects, &self.hooks, cancelled)
    }
}
