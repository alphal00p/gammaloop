use std::collections::BTreeMap;

use feynkit_cff::{
    CffGenerator, CffOptions, CffReport, CffResult, EdgeOrientation, FeynmanDiagramCffExt,
    OrientationExpression, Surface, SurfaceCache, SurfaceId,
};
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};
use symbolica::{api::python::PythonExpression, atom::Atom, symbol};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{error, graph::PyFeynmanDiagram};

/// A denominator surface appearing in a Cross-Free Family representation.
///
/// Surfaces identify the combinations of on-shell energies that can occur in
/// loop-energy denominators.  They are obtained from a ``CffResult`` rather
/// than constructed directly.
///
/// Examples
/// --------
/// >>> surface = next(iter(result.surfaces))
/// >>> print(surface, surface.positive_energies, surface.external_shift)
///
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
    /// Return the surface category: energy, h, unit, or infinite.
    #[getter]
    fn kind(&self) -> &'static str {
        match self.id {
            SurfaceId::Energy(_) => "energy",
            SurfaceId::H(_) => "h",
            SurfaceId::Unit => "unit",
            SurfaceId::Infinite => "infinite",
        }
    }

    /// Return the index of an energy or H surface.
    ///
    /// Raises :class:`CffError` for the special unit or infinite sentinels.
    #[getter]
    fn index(&self) -> PyResult<usize> {
        match self.id {
            SurfaceId::Energy(id) => Ok(id.index()),
            SurfaceId::H(id) => Ok(id.index()),
            SurfaceId::Unit | SurfaceId::Infinite => Err(error::CffError::new_err(format!(
                "the {} CFF sentinel has no surface index",
                self.kind()
            ))),
        }
    }

    /// Return the Symbolica variable name assigned to this surface.
    ///
    /// Raises :class:`CffError` for the special unit or infinite sentinels,
    /// which are not denominator variables.
    #[getter]
    fn symbol_name(&self) -> PyResult<String> {
        Self::symbol_name_for(self.id).ok_or_else(|| {
            error::CffError::new_err(format!(
                "the {} CFF sentinel has no Symbolica variable",
                self.kind()
            ))
        })
    }

    /// Return edge IDs whose on-shell energies enter with positive sign.
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

    /// Return edge IDs whose on-shell energies enter with negative sign.
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

    /// Return external edge IDs and their integer shift coefficients.
    #[getter]
    fn external_shift(&self) -> Vec<(usize, i64)> {
        match &self.surface {
            Some(Surface::Energy(surface)) => surface.external_shift.iter(),
            Some(Surface::H(surface)) => surface.external_shift.iter(),
            Some(Surface::Unit | Surface::Infinite) | None => return Vec::new(),
        }
        .map(|(edge, coefficient)| (edge.index(), *coefficient))
        .collect()
    }

    /// Return the dense CFF vertex IDs enclosed by this surface.
    ///
    /// CFF generation remaps the diagram's internal vertices to a compact
    /// zero-based range. These IDs therefore index the contracted CFF graph;
    /// they are not the original ``diagram.vertices`` identifiers.
    #[getter]
    fn vertices(&self) -> Vec<usize> {
        match &self.surface {
            Some(Surface::Energy(surface)) => surface.vertex_set.iter(),
            Some(Surface::H(surface)) => surface.vertex_set.iter(),
            Some(Surface::Unit | Surface::Infinite) | None => return Vec::new(),
        }
        .map(|vertex| vertex.index())
        .collect()
    }

    /// Return the surface's symbolic name, or its category for special surfaces.
    ///
    /// Examples
    /// --------
    /// >>> print(surface)  # Symbolica denominator variable, for example feynkit::E0
    ///
    fn __repr__(&self) -> String {
        Self::symbol_name_for(self.id).unwrap_or_else(|| self.kind().to_owned())
    }
}

/// One acyclic energy-flow orientation contributing to a CFF expression.
///
/// An orientation records the direction chosen for every diagram edge and the
/// products of denominator surfaces generated by that choice.
///
/// Examples
/// --------
/// >>> orientation = next(iter(result.orientations))
/// >>> for product in orientation.denominator_products():
/// ...     print([surface.symbol_name for surface in product])
///
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
    /// Return this orientation's stable index within the CFF result.
    #[getter]
    fn id(&self) -> usize {
        self.inner.id.index()
    }

    /// Return each edge ID and its selected orientation.
    #[getter]
    fn edge_orientations(&self) -> Vec<(usize, &'static str)> {
        self.inner
            .data
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

    /// Expand this orientation into products of denominator surfaces.
    ///
    /// Examples
    /// --------
    /// >>> products = orientation.denominator_products()
    /// >>> denominator_variables = [
    /// ...     [surface.symbol_name for surface in product] for product in products
    /// ... ]
    ///
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

/// Diagnostic counts collected while constructing a Cross-Free Family.
///
/// Use this report to compare the number of candidate energy flows with the
/// acyclic flows and denominator terms that survive.
///
/// Examples
/// --------
/// >>> report = result.report
/// >>> report.candidate_orientations >= report.acyclic_orientations
/// True
///
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
    /// Return the number of candidate edge orientations considered.
    #[getter]
    fn candidate_orientations(&self) -> usize {
        self.inner.candidate_orientations
    }
    /// Return the number of candidate orientations that are acyclic.
    #[getter]
    fn acyclic_orientations(&self) -> usize {
        self.inner.acyclic_orientations
    }
    /// Return the total number of unfolded denominator terms.
    #[getter]
    fn unfolded_terms(&self) -> usize {
        self.inner.unfolded_terms
    }
    /// Return the number of unique denominator surfaces in the result.
    #[getter]
    fn interned_surfaces(&self) -> usize {
        self.inner.interned_surfaces
    }

    /// Return a concise constructor-style CFF generation summary.
    ///
    /// Examples
    /// --------
    /// >>> print(result.report)
    ///
    fn __repr__(&self) -> String {
        format!(
            "CffReport(candidate_orientations={}, acyclic_orientations={}, unfolded_terms={}, interned_surfaces={})",
            self.inner.candidate_orientations,
            self.inner.acyclic_orientations,
            self.inner.unfolded_terms,
            self.inner.interned_surfaces,
        )
    }

    /// Render CFF generation statistics as a compact HTML table.
    ///
    /// Examples
    /// --------
    /// Leave ``result.report`` as the final expression in a notebook cell.
    ///
    fn _repr_html_(&self) -> String {
        format!(
            "<div class=\"feynkit-cff-report\" style=\"display:inline-block;max-width:100%;overflow-x:auto\">\
             <strong>CFF report</strong>\
             <table style=\"border-collapse:collapse;margin-top:.25rem\"><tbody>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">candidate orientations</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">acyclic orientations</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">unfolded terms</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             <tr><th style=\"padding:.2rem .65rem;text-align:left\">interned surfaces</th><td style=\"padding:.2rem .65rem;text-align:right\">{}</td></tr>\
             </tbody></table></div>",
            self.inner.candidate_orientations,
            self.inner.acyclic_orientations,
            self.inner.unfolded_terms,
            self.inner.interned_surfaces,
        )
    }

    /// Write a concise CFF summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "...".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// The Cross-Free Family representation of a Feynman diagram.
///
/// A result bundles the energy-flow orientations, their denominator surfaces,
/// generation statistics, and conversion to a native Symbolica expression.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> result = diagram.build_cff()
/// >>> expression = result.to_expression()
///
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
    /// Return generation statistics for this CFF result.
    #[getter]
    fn report(&self) -> PyCffReport {
        PyCffReport {
            inner: self.inner.report,
        }
    }

    /// Return the acyclic energy-flow orientations in this result.
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

    /// Return all unique energy and H surfaces in this result.
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
    ///
    /// Examples
    /// --------
    /// >>> expression = result.to_expression()
    /// >>> expression  # native Symbolica display in a notebook
    ///
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

    /// Return the number of unfolded denominator terms.
    ///
    /// Examples
    /// --------
    /// >>> denominator_term_count = len(result)
    ///
    fn __len__(&self) -> usize {
        self.inner.expression.unfolded_term_count()
    }

    /// Return a concise summary of the CFF expression and its surfaces.
    ///
    /// Examples
    /// --------
    /// >>> print(result)
    ///
    fn __repr__(&self) -> String {
        format!(
            "CffResult(orientations={}, terms={}, surfaces={})",
            self.inner.expression.orientations().len(),
            self.inner.expression.unfolded_term_count(),
            self.inner.surfaces.energy_surfaces().len() + self.inner.surfaces.h_surfaces().len(),
        )
    }

    /// Render the CFF report and its native Symbolica expression as HTML.
    ///
    /// The expression fragment comes from ``Expression._repr_html_`` so its
    /// Symbolica formatting is preserved in notebook output.
    ///
    /// Examples
    /// --------
    /// Leave ``result`` as the final expression in a notebook cell.
    ///
    fn _repr_html_(&self) -> PyResult<String> {
        let report = PyCffReport {
            inner: self.inner.report,
        }
        ._repr_html_();
        let expression = self.to_expression()._repr_html_()?;
        Ok(format!(
            "<section class=\"feynkit-cff-result\" style=\"max-width:100%\">\
             <h3 style=\"margin:.25rem 0\">Cross-free family</h3>{report}\
             <div style=\"margin-top:.65rem\"><strong>Expression</strong>\
             <div style=\"margin-top:.25rem;max-width:100%;overflow-x:auto\">{expression}</div>\
             </div></section>"
        ))
    }

    /// Write a summary with Symbolica's native expression formatting.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        if cycle {
            pretty.call_method1("text", ("...",))?;
            return Ok(());
        }

        pretty.call_method1(
            "text",
            (format!(
                "CffResult(orientations={}, terms={}, surfaces={}, expression=",
                self.inner.expression.orientations().len(),
                self.inner.expression.unfolded_term_count(),
                self.inner.surfaces.energy_surfaces().len()
                    + self.inner.surfaces.h_surfaces().len(),
            ),),
        )?;
        self.to_expression()._repr_pretty_(pretty, false)?;
        pretty.call_method1("text", (")",))?;
        Ok(())
    }
}

/// Generate Cross-Free Family energy denominators for Feynman diagrams.
///
/// The generator can constrain edge orientations or contractions before
/// enumerating the acyclic energy flows of a loop diagram.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> cff = fk.CffGenerator(max_orientations=10_000)
/// >>> result = cff.generate(diagram)
///
/// Parameters
/// ----------
/// max_orientations : int, optional
///     Maximum number of candidate edge orientations to inspect.
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
    /// Create a configurable Cross-Free Family generator.
    ///
    /// Examples
    /// --------
    /// >>> generator = CffGenerator(max_orientations=1000)
    ///
    /// Parameters
    /// ----------
    /// max_orientations : int, optional
    ///     Maximum number of candidate orientations to inspect.
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

    /// Fix an edge to its stored or reversed direction.
    ///
    /// Examples
    /// --------
    /// >>> generator.fix_orientation(0, reversed=True)
    ///
    /// Parameters
    /// ----------
    /// edge : int
    ///     Diagram edge ID.
    /// reversed : bool
    ///     Select the reversed direction when true.
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

    /// Contract an edge before generating denominator surfaces.
    ///
    /// Examples
    /// --------
    /// >>> generator.contract_edge(2)
    ///
    /// Parameters
    /// ----------
    /// edge : int
    ///     Diagram edge ID to contract.
    fn contract_edge(&mut self, edge: usize) {
        let options = self
            .inner
            .options()
            .clone()
            .with_contracted_edge(feynkit_cff::EdgeId::new(edge));
        self.inner = CffGenerator::new(options);
    }

    /// Mark an edge as belonging to the initial state.
    ///
    /// Examples
    /// --------
    /// >>> generator.mark_initial_state_edge(0)
    ///
    /// Parameters
    /// ----------
    /// edge : int
    ///     Diagram edge ID to classify as initial state.
    fn mark_initial_state_edge(&mut self, edge: usize) {
        let options = self
            .inner
            .options()
            .clone()
            .with_initial_state_edge(feynkit_cff::EdgeId::new(edge));
        self.inner = CffGenerator::new(options);
    }

    /// Generate a Cross-Free Family representation for a diagram.
    ///
    /// Examples
    /// --------
    /// >>> result = generator.generate(diagram)
    /// >>> result.to_expression()
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     Diagram whose energy-flow orientations are enumerated.
    fn generate(&self, py: Python<'_>, diagram: &PyFeynmanDiagram) -> PyResult<PyCffResult> {
        let options = self.inner.options().clone();
        let diagram = diagram.inner.clone();
        py.detach(move || diagram.build_cff(options))
            .map(|inner| PyCffResult { inner })
            .map_err(error::cff)
    }
}

pub(crate) fn build_cff_for_diagram(
    py: Python<'_>,
    diagram: &PyFeynmanDiagram,
    max_orientations: Option<usize>,
    fixed_orientations: Option<BTreeMap<usize, bool>>,
    contracted_edges: Option<Vec<usize>>,
    initial_state_edges: Option<Vec<usize>>,
) -> PyResult<PyCffResult> {
    let mut options = max_orientations.map_or_else(CffOptions::default, |maximum| {
        CffOptions::default().with_max_orientations(maximum)
    });
    for (edge, reversed) in fixed_orientations.unwrap_or_default() {
        options = options.with_fixed_orientation(
            feynkit_cff::EdgeId::new(edge),
            if reversed {
                EdgeOrientation::Reversed
            } else {
                EdgeOrientation::Default
            },
        );
    }
    for edge in contracted_edges.unwrap_or_default() {
        options = options.with_contracted_edge(feynkit_cff::EdgeId::new(edge));
    }
    for edge in initial_state_edges.unwrap_or_default() {
        options = options.with_initial_state_edge(feynkit_cff::EdgeId::new(edge));
    }

    let diagram = diagram.inner.clone();
    py.detach(move || diagram.build_cff(options))
        .map(|inner| PyCffResult { inner })
        .map_err(error::cff)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCffSurface>()?;
    module.add_class::<PyCffOrientation>()?;
    module.add_class::<PyCffReport>()?;
    module.add_class::<PyCffResult>()?;
    module.add_class::<PyCffGenerator>()?;
    Ok(())
}
