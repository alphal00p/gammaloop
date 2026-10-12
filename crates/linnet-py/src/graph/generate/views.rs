//! Typed read-only views. Borrowed callbacks expire before native search resumes; owned
//! skeleton views retain compact data and the Python key tuples, without a plan or execution.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use linnet::half_edge::generate::{self, VertexSet};
use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyIndexError, PyReferenceError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

use super::{Objects, PyEdgeSignature, PyGraph};

/// How far a typed generation view has been built.
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyclass_enum
)]
#[pyclass(from_py_object, frozen, eq, eq_int, hash, name = "GenerationStage")]
#[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
#[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum PyGenerationStage {
    /// An untyped topology with vertices, edges, attached legs and any remaining open stubs.
    Topology,
    /// A complete topology with fields assigned to some of its edges.
    Fields,
    /// Every edge has a field and every vertex has a rule.
    Complete,
}

impl From<generate::Stage> for PyGenerationStage {
    fn from(stage: generate::Stage) -> Self {
        match stage {
            generate::Stage::Topology => Self::Topology,
            generate::Stage::Fields => Self::Fields,
            generate::Stage::Complete => Self::Complete,
        }
    }
}

fn expired() -> PyErr {
    PyReferenceError::new_err("a callback generation view is only valid during its callback")
}

// Native typed views deliberately expose the same read-only methods. Generate their Python
// counterparts together so validation, expiry and payload conversion cannot diverge by stage.
macro_rules! typed_view {
    (
        $(#[$doc:meta])* $python:ident, $name:literal, $native:ident, $source:ident, $request:ident
        $(, $owned:ident($payload:ty) => $accessor:ident)* $(,)?
    ) => {
        enum $source {
            Request(Arc<generate::batch::Request<u32, u32>>),
            $($owned($payload),)*
        }

        $(#[$doc])*
        #[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
        #[pyclass(name = $name)]
        #[cfg_attr(feature = "community", pyo3(module = "symbolica.community.graph"))]
        #[cfg_attr(not(feature = "community"), pyo3(module = "linnet_py"))]
        pub(super) struct $python {
            source: Option<$source>,
            valid: Arc<AtomicBool>,
            objects: Option<Objects>,
        }

        impl $python {
            /// Own a compact callback request; the dispatcher expires its shared validity flag.
            pub(super) fn from_request(
                request: Arc<generate::batch::Request<u32, u32>>,
                valid: Arc<AtomicBool>,
                objects: Objects,
            ) -> Self {
                Self {
                    source: Some($source::Request(request)),
                    valid,
                    objects: Some(objects),
                }
            }

            /// Release the compact request and its Python keys as the callback ends.
            pub(super) fn expire(&mut self) {
                self.valid.store(false, Ordering::Release);
                self.source = None;
                self.objects = None;
            }

            pub(super) fn activate(&self) -> PyResult<()> {
                self.source.as_ref().ok_or_else(expired)?;
                self.valid.store(true, Ordering::Release);
                Ok(())
            }

            fn objects(&self) -> PyResult<&Objects> {
                self.objects.as_ref().ok_or_else(expired)
            }

            fn with<R>(
                &self,
                f: impl for<'v> FnOnce(&'v generate::$native<'v, u32, u32>) -> PyResult<R>,
            ) -> PyResult<R> {
                if !self.valid.load(Ordering::Acquire) {
                    return Err(expired());
                }
                match self.source.as_ref().ok_or_else(expired)? {
                    $source::Request(request) => {
                        let view = request.$request().ok_or_else(expired)?;
                        f(&view)
                    }
                    $($source::$owned(owner) => f(&owner.$accessor()),)*
                }
            }

            fn check_vertex(view: &generate::$native<'_, u32, u32>, vertex: usize) -> PyResult<()> {
                if vertex < view.n_vertices() {
                    Ok(())
                } else {
                    Err(PyIndexError::new_err(format!(
                        "vertex {vertex} is out of range for {} vertices",
                        view.n_vertices(),
                    )))
                }
            }

            fn at<R>(
                &self,
                vertex: usize,
                f: impl for<'v> FnOnce(&'v generate::$native<'v, u32, u32>, usize) -> PyResult<R>,
            ) -> PyResult<R> {
                self.with(|view| {
                    Self::check_vertex(view, vertex)?;
                    f(view, vertex)
                })
            }
        }

        #[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
        #[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
        #[pymethods]
        impl $python {
            #[getter]
            fn stage(&self) -> PyResult<PyGenerationStage> {
                self.with(|view| Ok(view.stage().into()))
            }

            /// The loop order being generated.
            #[getter]
            fn loop_order(&self) -> PyResult<usize> {
                self.with(|view| Ok(view.loop_order()))
            }

            #[getter]
            fn n_vertices(&self) -> PyResult<usize> {
                self.with(|view| Ok(view.n_vertices()))
            }

            /// Independent cycles of the edges placed so far.
            #[getter]
            fn loops(&self) -> PyResult<usize> {
                self.with(|view| Ok(view.loops()))
            }

            #[getter]
            fn unplaced_legs(&self) -> PyResult<usize> {
                self.with(|view| Ok(view.unplaced_legs()))
            }

            #[getter]
            fn is_topology_complete(&self) -> PyResult<bool> {
                self.with(|view| Ok(view.is_topology_complete()))
            }

            #[getter]
            fn is_complete(&self) -> PyResult<bool> {
                self.with(|view| Ok(view.is_complete()))
            }

            fn degree(&self, vertex: usize) -> PyResult<usize> {
                self.at(vertex, |view, v| Ok(view.degree(v)))
            }

            fn open_stubs(&self, vertex: usize) -> PyResult<usize> {
                self.at(vertex, |view, v| Ok(view.open_stubs(v)))
            }

            fn self_loops(&self, vertex: usize) -> PyResult<usize> {
                self.at(vertex, |view, v| Ok(view.self_loops(v)))
            }

            /// Number of edges between distinct vertices; use self_loops for equal ends.
            fn multiplicity(&self, u: usize, v: usize) -> PyResult<usize> {
                self.with(|view| {
                    Self::check_vertex(view, u)?;
                    Self::check_vertex(view, v)?;
                    if u == v {
                        return Err(PyValueError::new_err(
                            "multiplicity needs distinct vertices; use self_loops",
                        ));
                    }
                    Ok(view.multiplicity(u, v))
                })
            }

            /// Attached external labels and their field/direction at this vertex.
            #[gen_stub(override_return_type(type_repr = "list[tuple[_Key, EdgeSignature]]"))]
            fn legs(&self, py: Python<'_>, vertex: usize) -> PyResult<Vec<(Py<PyAny>, PyEdgeSignature)>> {
                self.at(vertex, |view, vertex| {
                    let objects = self.objects()?;
                    view.legs(vertex)
                        .map(|(&label, stub)| Ok((objects.label(py, label)?, objects.signature(py, stub)?)))
                        .collect()
                })
            }

            /// Internal edges `(u, v, stub)`; stub is None before its field is assigned.
            fn edges(&self, py: Python<'_>) -> PyResult<Vec<(usize, usize, Option<PyEdgeSignature>)>> {
                self.with(|view| {
                    let objects = self.objects()?;
                    view.edges()
                        .map(|edge| {
                            let [u, v] = edge.vertices;
                            Ok((u, v, edge.stub.map(|stub| objects.signature(py, stub)).transpose()?))
                        })
                        .collect()
                })
            }

            /// Exact matching rules once all vertex stubs have fields; otherwise None.
            fn rules(&self, vertex: usize) -> PyResult<Option<Vec<u32>>> {
                self.at(vertex, |view, v| Ok(view.rules(v).map(<[u32]>::to_vec)))
            }

            /// Compatible rules in declaration order, including rules for partial assignments.
            fn domains(&self, vertex: usize) -> PyResult<Vec<u32>> {
                self.at(vertex, |view, v| Ok(view.domains(v)))
            }

            /// Inclusive named weight bounds once the skeleton is complete and feasible.
            #[gen_stub(override_return_type(type_repr = "dict[str, tuple[int, int]] | None"))]
            fn order_bounds(&self, py: Python<'_>) -> PyResult<Option<Py<PyDict>>> {
                self.with(|view| {
                    let Some(bounds) = view.order_bounds() else { return Ok(None) };
                    let result = PyDict::new(py);
                    for bound in bounds {
                        result.set_item(bound.name, (bound.range[0], bound.range[1]))?;
                    }
                    Ok(Some(result.unbind()))
                })
            }

            /// Stable untyped topology identity; None while edges or external labels are missing.
            #[gen_stub(override_return_type(type_repr = "tuple[typing.Any, ...] | None", imports=("typing")))]
            fn skeleton_id(&self, py: Python<'_>) -> PyResult<Option<Py<PyTuple>>> {
                self.with(|view| {
                    view.skeleton_id()
                        .map(|form| self.objects()?.skeleton_form(py, form))
                        .transpose()
                })
            }

            /// Bridges with a final side, as `(outside, inside, side_vertices)`.
            fn bridges(&self) -> PyResult<Vec<(usize, usize, Vec<usize>)>> {
                self.with(|view| Ok(view.bridges().into_iter().map(|bridge| {
                    let [outside, inside] = bridge.ends;
                    (outside, inside, bridge.side.iter().collect())
                }).collect()))
            }

            /// Independent cycles on a nonempty connected set of vertices.
            fn loops_in(&self, side: Vec<usize>) -> PyResult<usize> {
                self.with(|view| {
                    for &vertex in &side {
                        Self::check_vertex(view, vertex)?;
                    }
                    let side = side.into_iter().collect::<VertexSet>();
                    if !view.is_connected(side) {
                        return Err(PyValueError::new_err(
                            "loops_in needs a non-empty set of vertices connected by placed edges",
                        ));
                    }
                    Ok(view.loops_in(side))
                })
            }

            /// Materialize the current state, omitting open slots and using None for unknown keys.
            fn graph(&self, py: Python<'_>) -> PyResult<PyGraph> {
                self.with(|view| self.objects()?.graph(py, &view.graph(), |&field| field))
            }

            /// Materialize a fully assigned graph; None before the Complete stage.
            fn complete_graph(&self, py: Python<'_>) -> PyResult<Option<PyGraph>> {
                self.with(|view| {
                    view.complete_graph()
                        .map(|graph| self.objects()?.graph(py, graph, |&field| Some(field)))
                        .transpose()
                })
            }

            #[gen_stub(skip)]
            fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
                if let Some(objects) = &self.objects {
                    objects.traverse(&visit)?;
                }
                Ok(())
            }

            #[gen_stub(skip)]
            fn __clear__(&mut self) {
                self.expire();
            }
        }
    };
}

typed_view!(
    /// An untyped topology, borrowed during a callback or retained as compact owned data.
    PySkeletonView, "SkeletonView", SkeletonView, SkeletonSource, skeleton,
    Snapshot(generate::SkeletonSnapshot<u32, u32>) => view,
    Diagram(Arc<generate::Diagram<u32, u32>>) => skeleton,
);
typed_view!(
    /// A scoped complete topology with fields and rules being assigned.
    PyDressingView, "DressingView", DressingView, DressingSource, dressing,
);
typed_view!(
    /// A scoped fully assigned diagram supplied to final selection.
    PyDiagramView, "DiagramView", DiagramView, DiagramSource, diagram,
);

impl PySkeletonView {
    pub(super) fn from_snapshot(
        snapshot: generate::SkeletonSnapshot<u32, u32>,
        objects: Objects,
    ) -> Self {
        Self {
            source: Some(SkeletonSource::Snapshot(snapshot)),
            valid: Arc::new(AtomicBool::new(true)),
            objects: Some(objects),
        }
    }

    pub(super) fn from_diagram(
        diagram: Arc<generate::Diagram<u32, u32>>,
        objects: Objects,
    ) -> Self {
        Self {
            source: Some(SkeletonSource::Diagram(diagram)),
            valid: Arc::new(AtomicBool::new(true)),
            objects: Some(objects),
        }
    }
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyGenerationStage>()?;
    module.add_class::<PySkeletonView>()?;
    module.add_class::<PyDressingView>()?;
    module.add_class::<PyDiagramView>()?;
    Ok(())
}
