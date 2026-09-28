//! Interoperate with the installed Linnet extension, not a second set of PyO3 classes.

use std::{collections::BTreeSet, sync::Mutex};

use linnet::half_edge::{
    involution::{Flow, Hedge, HedgePair, Orientation},
    subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
};
use pyo3::{
    PyTraverseError, PyVisit,
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
    types::{PyDict, PyTuple},
};

use crate::graph::PyFeynmanDiagram;

#[derive(Default)]
pub(crate) struct LinnetCache(Mutex<Option<Py<LinnetCacheHolder>>>);

// Each physics wrapper owns one counted Python reference to this holder. The
// holder owns the graph reference once, so shared refreshes remain visible to
// every view without reporting one graph reference repeatedly to Python's GC.
#[pyclass(frozen, name = "_LinnetCache", module = "symbolica.community.feynkit")]
#[derive(Default)]
struct LinnetCacheHolder(Mutex<Option<LinnetExport>>);

struct LinnetExport {
    graph: Py<PyAny>,
    revision: u64,
    /// Indexed by the exported half-edge ID; entries are native half-edge IDs.
    native_hedges: Vec<usize>,
    exported_hedges: Vec<usize>,
}

impl Clone for LinnetCache {
    fn clone(&self) -> Self {
        Python::attach(|py| {
            Self(Mutex::new(Some(
                self.holder(py).expect("initialize shared Linnet cache"),
            )))
        })
    }
}

impl LinnetCache {
    fn holder(&self, py: Python<'_>) -> PyResult<Py<LinnetCacheHolder>> {
        let cached = self
            .0
            .lock()
            .map_err(|_| PyRuntimeError::new_err("Linnet cache poisoned"))?
            .as_ref()
            .map(|holder| holder.clone_ref(py));
        if let Some(holder) = cached {
            return Ok(holder);
        }
        // Allocation may trigger GC traversal of this wrapper; keep it outside
        // the mutex, including on the first clone before any graph is exported.
        let candidate = Py::new(py, LinnetCacheHolder::default())?;
        let mut cache = self
            .0
            .lock()
            .map_err(|_| PyRuntimeError::new_err("Linnet cache poisoned"))?;
        let holder = cache
            .get_or_insert_with(|| candidate.clone_ref(py))
            .clone_ref(py);
        drop(cache);
        Ok(holder)
    }

    pub(crate) fn traverse(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Ok(cache) = self.0.lock()
            && let Some(holder) = &*cache
        {
            visit.call(holder)?;
        }
        Ok(())
    }

    pub(crate) fn clear(&self) {
        let previous = self.0.lock().ok().and_then(|mut cache| cache.take());
        drop(previous);
    }

    pub(crate) fn graph(&self, py: Python<'_>, diagram: &PyFeynmanDiagram) -> PyResult<Py<PyAny>> {
        self.holder(py)?.borrow(py).graph(py, diagram)
    }

    pub(crate) fn selection(
        &self,
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        selection: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<SuBitGraph> {
        self.holder(py)?
            .borrow(py)
            .selection(py, diagram, selection)
    }

    pub(crate) fn export_selection(
        &self,
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        selection: &SuBitGraph,
        isolated: &BTreeSet<usize>,
    ) -> PyResult<Py<PyAny>> {
        self.holder(py)?
            .borrow(py)
            .export_selection(py, diagram, selection, isolated)
    }
}

#[pymethods]
impl LinnetCacheHolder {
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        if let Ok(cache) = self.0.lock()
            && let Some(export) = &*cache
        {
            visit.call(&export.graph)?;
        }
        Ok(())
    }

    fn __clear__(&self) {
        let previous = self.0.lock().ok().and_then(|mut cache| cache.take());
        drop(previous);
    }
}

impl LinnetCacheHolder {
    fn graph(&self, py: Python<'_>, diagram: &PyFeynmanDiagram) -> PyResult<Py<PyAny>> {
        // Never hold this lock across Python calls: GC can traverse the diagram
        // while a Python factory allocates or a payload destructor runs.
        let cached = {
            let cache = self
                .0
                .lock()
                .map_err(|_| PyRuntimeError::new_err("Linnet cache poisoned"))?;
            cache
                .as_ref()
                .map(|export| (export.graph.clone_ref(py), export.revision))
        };
        if let Some((graph, revision)) = cached {
            let current: u64 = graph
                .bind(py)
                .call_method0("full_subgraph")?
                .getattr("revision")?
                .extract()?;
            if current == revision {
                return Ok(graph);
            }
        }
        let export = LinnetExport::build(py, diagram)?;
        let graph = export.graph.clone_ref(py);
        let previous = self
            .0
            .lock()
            .map_err(|_| PyRuntimeError::new_err("Linnet cache poisoned"))?
            .replace(export);
        drop(previous);
        Ok(graph)
    }

    fn selection(
        &self,
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        selection: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<SuBitGraph> {
        let Some(selection) = selection else {
            return Ok(diagram.inner.underlying().full_filter());
        };
        let module = py.import("linnet")?;
        if !selection.is_instance(&module.getattr("Subgraph")?)? {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "subgraph must be a linnet.Subgraph from diagram.to_linnet()",
            ));
        }
        let graph = self.graph(py, diagram)?;
        if !selection.getattr("graph")?.is(graph.bind(py)) {
            return Err(PyValueError::new_err(
                "subgraph belongs to a different diagram or topology revision",
            ));
        }
        let revision: u64 = selection.getattr("revision")?.extract()?;
        let hedges: Vec<usize> = selection.call_method0("half_edge_indices")?.extract()?;
        let cache = self
            .0
            .lock()
            .map_err(|_| PyRuntimeError::new_err("Linnet cache poisoned"))?;
        let export = cache.as_ref().expect("export initialized above");
        if revision != export.revision {
            return Err(PyValueError::new_err(
                "subgraph belongs to a different topology revision",
            ));
        }
        let mut native = SuBitGraph::empty(diagram.inner.underlying().n_hedges());
        for hedge in hedges {
            native.add(Hedge(export.native_hedges[hedge]));
        }
        Ok(native)
    }

    fn export_selection(
        &self,
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        selection: &SuBitGraph,
        isolated: &BTreeSet<usize>,
    ) -> PyResult<Py<PyAny>> {
        let graph = self.graph(py, diagram)?;
        let hedges: Vec<_> = {
            let cache = self
                .0
                .lock()
                .map_err(|_| PyRuntimeError::new_err("Linnet cache poisoned"))?;
            let export = cache.as_ref().expect("export initialized above");
            selection
                .included_iter()
                .map(|hedge| export.exported_hedges[hedge.0])
                .collect()
        };
        let kwargs = PyDict::new(py);
        kwargs.set_item("half_edges", hedges)?;
        kwargs.set_item("nodes", isolated.iter().copied().collect::<Vec<_>>())?;
        Ok(graph
            .bind(py)
            .call_method("subgraph", (), Some(&kwargs))?
            .unbind())
    }
}

impl LinnetExport {
    fn build(py: Python<'_>, diagram: &PyFeynmanDiagram) -> PyResult<Self> {
        let module = py.import("linnet")?;
        // A selection retains its parent's topology and stable element IDs.
        // Export the complete owner graph even when the caller is a view.
        let diagram = diagram.whole();
        let native = diagram.inner.underlying();
        let mut items = Vec::new();
        for vertex in diagram.vertices() {
            let kwargs = PyDict::new(py);
            kwargs.set_item("data", Py::new(py, vertex)?)?;
            items.push(module.getattr("node")?.call((), Some(&kwargs))?.unbind());
        }
        let edges = diagram.edges();
        for (pair, edge, data) in native.iter_edges() {
            let endpoint = |hedge: Hedge, flow: Flow| -> PyResult<Py<PyAny>> {
                let kwargs = PyDict::new(py);
                kwargs.set_item("data", hedge.0)?;
                Ok(module
                    .getattr(match flow {
                        Flow::Source => "source",
                        Flow::Sink => "sink",
                    })?
                    .call((native.node_id(hedge).0,), Some(&kwargs))?
                    .unbind())
            };
            let (first, second) = match pair {
                HedgePair::Paired { source, sink } => (
                    endpoint(source, Flow::Source)?,
                    Some(endpoint(sink, Flow::Sink)?),
                ),
                HedgePair::Unpaired { hedge, flow } => (endpoint(hedge, flow)?, None),
                HedgePair::Split { .. } => unreachable!("whole diagrams have no split edges"),
            };
            let kwargs = PyDict::new(py);
            kwargs.set_item("data", Py::new(py, edges[edge.0].clone())?)?;
            kwargs.set_item(
                "orientation",
                module
                    .getattr("Orientation")?
                    .getattr(match data.orientation {
                        Orientation::Default => "Default",
                        Orientation::Reversed => "Reversed",
                        Orientation::Undirected => "Undirected",
                    })?,
            )?;
            items.push(
                module
                    .getattr("edge")?
                    .call((first, format!("e{}", edge.0), second), Some(&kwargs))?
                    .unbind(),
            );
        }
        let kwargs = PyDict::new(py);
        kwargs.set_item("name", diagram.inner.name())?;
        let graph = module
            .getattr("build")?
            .call(PyTuple::new(py, items)?, Some(&kwargs))?;
        let mut native_hedges = vec![0; native.n_hedges()];
        let mut exported_hedges = vec![0; native.n_hedges()];
        for hedge in graph.call_method0("half_edges")?.try_iter()? {
            let hedge = hedge?;
            let exported: usize = hedge.getattr("index")?.extract()?;
            let native: usize = hedge.getattr("data")?.extract()?;
            native_hedges[exported] = native;
            exported_hedges[native] = exported;
        }
        let revision = graph
            .call_method0("full_subgraph")?
            .getattr("revision")?
            .extract()?;
        Ok(Self {
            graph: graph.unbind(),
            revision,
            native_hedges,
            exported_hedges,
        })
    }
}
