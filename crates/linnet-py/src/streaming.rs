//! Incremental native force layout for notebook previews.

use std::collections::BTreeMap;

use linnest::ForceLayoutStream;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// A position snapshot from one continuous force-layout run.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(frozen, name = "LayoutFrame", get_all)]
pub(crate) struct PyLayoutFrame {
    nodes: Vec<(f64, f64)>,
    edges: Vec<(f64, f64)>,
    iteration: usize,
    done: bool,
    max_movement: f64,
}

/// A synchronous iterator that advances the native force solver in batches.
///
/// Topology is fixed for the lifetime of the stream. Frames contain only
/// coordinates and progress, so notebook viewers can retain their SVG elements.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(unsendable, name = "LayoutStream")]
pub(crate) struct PyLayoutStream {
    stream: ForceLayoutStream,
    every: usize,
    started: bool,
    finished: bool,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLayoutStream {
    /// Start a force-layout preview from a single DOT graph.
    ///
    /// This preview uses native graph geometry, without Typst label measurement
    /// or final edge-label relaxation. A new stream restarts from the same seed.
    #[staticmethod]
    #[pyo3(signature = (dot, *, every=4, steps=200, epochs=8, seed=1, step=0.02, cool=0.85, spring_strength=1.0, repulsion=1.5, length_scale=1.0, depth_scale=1.0, flattening_end=0.5, delta=0.1, early_tolerance=0.000001))]
    #[allow(clippy::too_many_arguments)]
    fn from_dot(
        dot: &str,
        every: usize,
        steps: usize,
        epochs: usize,
        seed: u64,
        step: f64,
        cool: f64,
        spring_strength: f64,
        repulsion: f64,
        length_scale: f64,
        depth_scale: f64,
        flattening_end: f64,
        delta: f64,
        early_tolerance: f64,
    ) -> PyResult<Self> {
        if every == 0 {
            return Err(PyValueError::new_err("every must be greater than zero"));
        }
        let mut settings = BTreeMap::new();
        for (key, value) in [
            ("step", step),
            ("cool", cool),
            ("k-spring", spring_strength),
            ("beta", repulsion),
            ("length-scale", length_scale),
            ("depth-scale", depth_scale),
            ("flattening-end", flattening_end),
            ("delta", delta),
            ("early-tol", early_tolerance),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(PyValueError::new_err(format!(
                    "{key} must be a non-negative finite number"
                )));
            }
            settings.insert(key, value.to_string());
        }
        if !(0.0..=1.0).contains(&cool) || !(0.0..=1.0).contains(&flattening_end) {
            return Err(PyValueError::new_err(
                "cool and flattening_end must be between zero and one",
            ));
        }
        if length_scale == 0.0 {
            return Err(PyValueError::new_err("length_scale must be positive"));
        }
        settings.insert("steps", steps.to_string());
        settings.insert("epochs", epochs.to_string());
        settings.insert("seed", seed.to_string());
        settings.insert("layout-algo", "force".to_owned());
        let mut bytes = Vec::new();
        ciborium::ser::into_writer(&settings, &mut bytes)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let stream = ForceLayoutStream::from_dot(dot, &bytes).map_err(PyValueError::new_err)?;
        Ok(Self {
            stream,
            every,
            started: false,
            finished: false,
        })
    }

    /// Node names in the same stable index order as every frame's coordinates.
    #[getter]
    fn node_names(&self) -> Vec<String> {
        self.stream.node_names()
    }

    /// Source and sink node indices; a missing endpoint denotes a dangling edge.
    #[getter]
    fn endpoints(&self) -> Vec<(Option<usize>, Option<usize>)> {
        self.stream.endpoints()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    #[gen_stub(override_return_type(type_repr = "LayoutFrame"))]
    fn __next__(&mut self) -> Option<PyLayoutFrame> {
        if self.finished {
            return None;
        }
        let count = if self.started { self.every } else { 0 };
        self.started = true;
        let frame = self.stream.step(count);
        self.finished = frame.done;
        Some(PyLayoutFrame {
            nodes: frame
                .nodes
                .into_iter()
                .map(|point| (point.x, point.y))
                .collect(),
            edges: frame
                .edges
                .into_iter()
                .map(|point| (point.x, point.y))
                .collect(),
            iteration: frame.iteration,
            done: frame.done,
            max_movement: frame.max_movement,
        })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyLayoutFrame>()?;
    module.add_class::<PyLayoutStream>()?;
    Ok(())
}
