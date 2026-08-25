use feynkit_cff::{
    CffGenerator, CffGraph, CffOptions, CffReport, CffResult, EdgeOrientation,
    OrientationExpression, Surface, SurfaceCache, SurfaceId,
};
use pyo3::{prelude::*, types::PyModule};
use symbolica::{api::python::PythonExpression, atom::Atom, symbol};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{error, graph::PyFeynmanDiagram};

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffSurface",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffSurface {
    id: SurfaceId,
    surface: Option<Surface>,
}

impl PyCffSurface {
    fn new(id: SurfaceId, cache: &SurfaceCache) -> Self {
        Self {
            id,
            surface: cache.get(id),
        }
    }

    fn symbol_name_for(id: SurfaceId) -> Option<String> {
        match id {
            SurfaceId::Energy(id) => Some(format!("feynkit::E{}", id.index())),
            SurfaceId::H(id) => Some(format!("feynkit::H{}", id.index())),
            SurfaceId::Unit | SurfaceId::Infinite => None,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCffSurface {
    #[getter]
    fn kind(&self) -> &'static str {
        match self.id {
            SurfaceId::Energy(_) => "energy",
            SurfaceId::H(_) => "h",
            SurfaceId::Unit => "unit",
            SurfaceId::Infinite => "infinite",
        }
    }

    #[getter]
    fn index(&self) -> Option<usize> {
        match self.id {
            SurfaceId::Energy(id) => Some(id.index()),
            SurfaceId::H(id) => Some(id.index()),
            SurfaceId::Unit | SurfaceId::Infinite => None,
        }
    }

    #[getter]
    fn symbol_name(&self) -> Option<String> {
        Self::symbol_name_for(self.id)
    }

    #[getter]
    fn positive_energies(&self) -> Vec<usize> {
        match &self.surface {
            Some(Surface::Energy(surface)) => {
                surface.energies.iter().map(|edge| edge.index()).collect()
            }
            Some(Surface::H(surface)) => surface
                .positive_energies
                .iter()
                .map(|edge| edge.index())
                .collect(),
            Some(Surface::Unit | Surface::Infinite) | None => Vec::new(),
        }
    }

    #[getter]
    fn negative_energies(&self) -> Vec<usize> {
        match &self.surface {
            Some(Surface::H(surface)) => surface
                .negative_energies
                .iter()
                .map(|edge| edge.index())
                .collect(),
            Some(Surface::Energy(_) | Surface::Unit | Surface::Infinite) | None => Vec::new(),
        }
    }

    #[getter]
    fn external_shift(&self) -> Vec<(usize, i64)> {
        match &self.surface {
            Some(Surface::Energy(surface)) => surface.external_shift.iter(),
            Some(Surface::H(surface)) => surface.external_shift.iter(),
            Some(Surface::Unit | Surface::Infinite) | None => return Vec::new(),
        }
        .map(|(edge, coefficient)| (edge.index(), coefficient))
        .collect()
    }

    #[getter]
    fn vertices(&self) -> Vec<usize> {
        match &self.surface {
            Some(Surface::Energy(surface)) => surface.vertices.iter(),
            Some(Surface::H(surface)) => surface.vertices.iter(),
            Some(Surface::Unit | Surface::Infinite) | None => return Vec::new(),
        }
        .map(|vertex| vertex.index())
        .collect()
    }

    fn __repr__(&self) -> String {
        self.symbol_name().unwrap_or_else(|| self.kind().to_owned())
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffOrientation",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffOrientation {
    inner: OrientationExpression,
    surfaces: SurfaceCache,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCffOrientation {
    #[getter]
    fn id(&self) -> usize {
        self.inner.id.index()
    }

    #[getter]
    fn edge_orientations(&self) -> Vec<(usize, &'static str)> {
        self.inner
            .orientation
            .iter()
            .map(|(edge, orientation)| {
                let orientation = match orientation {
                    EdgeOrientation::Default => "default",
                    EdgeOrientation::Reversed => "reversed",
                    EdgeOrientation::Undirected => "undirected",
                };
                (edge.index(), orientation)
            })
            .collect()
    }

    fn denominator_products(&self) -> Vec<Vec<PyCffSurface>> {
        self.inner
            .denominator_products()
            .into_iter()
            .map(|term| {
                term.into_iter()
                    .map(|surface| PyCffSurface::new(surface, &self.surfaces))
                    .collect()
            })
            .collect()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffReport",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffReport {
    inner: CffReport,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCffReport {
    #[getter]
    fn candidate_orientations(&self) -> usize {
        self.inner.candidate_orientations
    }
    #[getter]
    fn acyclic_orientations(&self) -> usize {
        self.inner.acyclic_orientations
    }
    #[getter]
    fn unfolded_terms(&self) -> usize {
        self.inner.unfolded_terms
    }
    #[getter]
    fn interned_surfaces(&self) -> usize {
        self.inner.interned_surfaces
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffResult",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffResult {
    inner: CffResult,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCffResult {
    #[getter]
    fn report(&self) -> PyCffReport {
        PyCffReport {
            inner: self.inner.report,
        }
    }

    #[getter]
    fn orientations(&self) -> Vec<PyCffOrientation> {
        self.inner
            .expression
            .orientations()
            .iter()
            .cloned()
            .map(|inner| PyCffOrientation {
                inner,
                surfaces: self.inner.surfaces.clone(),
            })
            .collect()
    }

    #[getter]
    fn surfaces(&self) -> Vec<PyCffSurface> {
        let energy = (0..self.inner.surfaces.energy_surfaces().len()).map(|index| {
            PyCffSurface::new(
                SurfaceId::Energy(feynkit_cff::EnergySurfaceId(index)),
                &self.inner.surfaces,
            )
        });
        let h = (0..self.inner.surfaces.h_surfaces().len()).map(|index| {
            PyCffSurface::new(
                SurfaceId::H(feynkit_cff::HSurfaceId(index)),
                &self.inner.surfaces,
            )
        });
        energy.chain(h).collect()
    }

    /// Convert the typed denominator forest into a native Symbolica expression.
    ///
    /// Surface variables are named `feynkit::E<N>` and `feynkit::H<N>`; their
    /// definitions remain available through `surfaces`.
    fn to_expression(&self) -> PythonExpression {
        let mut expression = Atom::Zero;
        for orientation in self.inner.expression.orientations() {
            for product in orientation.denominator_products() {
                let mut term = Atom::num(1);
                for surface in product {
                    match surface {
                        SurfaceId::Unit => {}
                        SurfaceId::Infinite => {
                            term = Atom::Zero;
                            break;
                        }
                        id => {
                            let Some(name) = PyCffSurface::symbol_name_for(id) else {
                                continue;
                            };
                            term /= Atom::var(symbol!(&name));
                        }
                    }
                }
                expression += term;
            }
        }
        PythonExpression { expr: expression }
    }

    fn __len__(&self) -> usize {
        self.inner.expression.unfolded_term_count()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffGenerator",
    module = "symbolica.community.feynkit",
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffGenerator {
    inner: CffGenerator,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCffGenerator {
    #[new]
    #[pyo3(signature = (max_orientations=None))]
    fn new(max_orientations: Option<usize>) -> Self {
        let options = max_orientations.map_or_else(CffOptions::default, |maximum| {
            CffOptions::default().with_max_orientations(maximum)
        });
        Self {
            inner: CffGenerator::new(options),
        }
    }

    fn fix_orientation(&mut self, edge: usize, reversed: bool) {
        let orientation = if reversed {
            EdgeOrientation::Reversed
        } else {
            EdgeOrientation::Default
        };
        let options = self
            .inner
            .options()
            .clone()
            .with_fixed_orientation(feynkit_cff::EdgeId::new(edge), orientation);
        self.inner = CffGenerator::new(options);
    }

    fn contract_edge(&mut self, edge: usize) {
        let options = self
            .inner
            .options()
            .clone()
            .with_contracted_edge(feynkit_cff::EdgeId::new(edge));
        self.inner = CffGenerator::new(options);
    }

    fn mark_initial_state_edge(&mut self, edge: usize) {
        let options = self
            .inner
            .options()
            .clone()
            .with_initial_state_edge(feynkit_cff::EdgeId::new(edge));
        self.inner = CffGenerator::new(options);
    }

    fn generate(&self, py: Python<'_>, diagram: &PyFeynmanDiagram) -> PyResult<PyCffResult> {
        let generator = self.inner.clone();
        let diagram = diagram.inner.clone();
        py.detach(move || {
            let graph = CffGraph::try_from(&diagram)?;
            generator.generate(&graph)
        })
        .map(|inner| PyCffResult { inner })
        .map_err(error::cff)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCffSurface>()?;
    module.add_class::<PyCffOrientation>()?;
    module.add_class::<PyCffReport>()?;
    module.add_class::<PyCffResult>()?;
    module.add_class::<PyCffGenerator>()?;
    Ok(())
}
