mod display;
mod kernel;
use kernel::{CffKernel, CffOrientationData, CffTerm};

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{Arc, OnceLock},
};

use feynkit_cff::{
    CffExpression, CffOptions, CffReport, CffResult, CutPropagator, EdgeOrientation,
    FeynmanDiagramCffExt, SurfaceId, SurfacePole,
};
use linnet::half_edge::{involution::HedgePair, subgraph::SuBitGraph};
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};
use spynso3::expression::TensorExpression;
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::{Atom, AtomCore, AtomView, Symbol},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{
    energy::{PyEnergySurface, PyOnShellEnergy, PySurfaceFactor},
    error,
    graph::PyFeynmanDiagram,
};

/// One acyclic energy-flow orientation in a CFF representation.
///
/// ``families`` contains its cross-free denominator products.
/// ``to_expression()`` sums their contributions, including the common
/// energy prefactors and any supplied numerator. Edge directions are relative to the
/// originating diagram's stored orientation.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> result = diagram.integrate_energy(method="cff")
/// >>> orientation = result.orientations[0]
/// >>> expression = orientation.to_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffOrientation",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffOrientation {
    inner: CffOrientationData,
    parent: PyCffRepresentation,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCffOrientation {
    /// Return this orientation's stable index within the CFF result.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffOrientation`` class example:
    ///
    /// >>> orientation_ids = [item.id for item in result.orientations]
    #[getter]
    fn id(&self) -> usize {
        self.inner.id
    }

    /// Edge signs relative to the diagram's stored arrows.
    ///
    /// +1 agrees with the stored direction, -1 reverses it, and None denotes
    /// an undirected edge. These are orientation signs, not surface coefficients.
    ///
    /// Examples
    /// --------
    /// >>> signs = orientation.edge_signs
    /// >>> assert all(sign in (-1, 1, None) for sign in signs.values())
    #[getter]
    fn edge_signs(&self) -> BTreeMap<usize, Option<i8>> {
        self.inner
            .directions
            .iter()
            .map(|(edge, direction)| {
                (
                    *edge,
                    match direction {
                        EdgeOrientation::Default => Some(1),
                        EdgeOrientation::Reversed => Some(-1),
                        EdgeOrientation::Undirected => None,
                    },
                )
            })
            .collect()
    }

    /// Summarize this orientation without displaying its parent representation.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(orientation)
    fn __repr__(&self) -> String {
        format!(
            "CffOrientation(id={}, families={})",
            self.id(),
            self.families().len()
        )
    }

    /// Return each edge ID and its selected orientation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffOrientation`` class example:
    ///
    /// >>> directions = dict(orientation.edge_orientations)
    #[getter]
    fn edge_orientations(&self) -> Vec<(usize, &'static str)> {
        self.inner
            .directions
            .iter()
            .map(|(edge, orientation)| {
                let orientation = match orientation {
                    EdgeOrientation::Default => "default",
                    EdgeOrientation::Reversed => "reversed",
                    EdgeOrientation::Undirected => "undirected",
                };
                (*edge, orientation)
            })
            .collect()
    }

    /// Cross-free families contributing to this orientation.
    ///
    /// Examples
    /// --------
    /// >>> families = orientation.families
    #[getter]
    fn families(&self) -> Vec<PyCrossFreeFamily> {
        self.inner
            .terms
            .iter()
            .cloned()
            .enumerate()
            .map(|(id, term)| PyCrossFreeFamily {
                parent: self.parent.clone(),
                orientation: self.id(),
                id,
                term,
            })
            .collect()
    }
    /// Sum this orientation's families, including the common energy prefactor.
    ///
    /// Examples
    /// --------
    /// Using the setup in the class example:
    ///
    /// >>> expression = orientation.to_expression()
    ///
    /// Parameters
    /// ----------
    /// expand_surfaces : bool, optional
    ///     Replace stored surface placeholders by affine energy definitions.
    ///     Defaults to True. On-shell energies remain symbolic in either case.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     This orientation's contribution to the CFF energy integral.
    #[pyo3(signature=(*, expand_surfaces=true))]
    fn to_expression(
        &self,
        py: Python<'_>,
        expand_surfaces: bool,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(
            py,
            self.parent
                .expression(self.inner.atom(), expand_surfaces)
                .expr,
            None,
        )
    }
    /// Inspect this orientation's family sum and fixed energy-flow graph.
    ///
    /// Examples
    /// --------
    /// >>> html = orientation._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.parent
            .explorer_html(py, display::CffScope::Orientation(self.id()))
    }
    /// Display this orientation in Marimo.
    ///
    /// Examples
    /// --------
    /// >>> view = orientation._display_()
    fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py
            .import("marimo")?
            .call_method1("iframe", (self._repr_html_(py)?,))?
            .unbind())
    }
}

/// One cross-free family of denominator surfaces in a CFF orientation.
///
/// ``factors`` records signed surfaces and their multiplicities. Its expression
/// is the product of ``coefficient``, ``numerator`` and the inverse of every
/// entry in ``factors``. Summing these expressions reproduces the parent
/// orientation. Generalized numerator sampling can produce several weighted
/// contributions with the same geometric family but different ``energy_map``s.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> result = diagram.integrate_energy(method="cff")
/// >>> family = result.orientations[0].families[0]
/// >>> expression = family.to_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CrossFreeFamily",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCrossFreeFamily {
    parent: PyCffRepresentation,
    orientation: usize,
    id: usize,
    term: CffTerm,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCrossFreeFamily {
    /// Summarize this family's orientation, index and denominator factors.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(family)
    fn __repr__(&self) -> String {
        format!(
            "CrossFreeFamily(orientation={}, id={}, factors={})",
            self.orientation,
            self.id,
            self.factors().len()
        )
    }
    /// Family index within its orientation.
    ///
    /// Examples
    /// --------
    /// >>> index = family.id
    #[getter]
    fn id(&self) -> usize {
        self.id
    }
    /// Signed denominator factors, grouping repeated surfaces into powers.
    ///
    /// Examples
    /// --------
    /// >>> denominators = [f.to_expression() for f in family.factors]
    #[getter]
    fn factors(&self) -> Vec<PySurfaceFactor> {
        let mut powers = BTreeMap::new();
        for id in &self.term.path {
            *powers.entry(*id).or_insert(0) += 1;
        }
        powers
            .into_iter()
            .map(|(id, power)| PySurfaceFactor {
                surface: self.parent.inner.surfaces[id].clone(),
                sign: 1,
                power,
            })
            .collect()
    }
    /// Signed coefficient and on-shell factors multiplying this family's numerator.
    ///
    /// Examples
    /// --------
    /// >>> coefficient = family.coefficient
    #[getter]
    fn coefficient(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(
            py,
            &self.parent.normalization * self.term.prefactor(),
            None,
        )
    }
    /// Evaluated numerator, including any generated numerator surface factors.
    ///
    /// Examples
    /// --------
    /// >>> numerator = family.numerator
    #[getter]
    fn numerator(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        let value = self
            .term
            .numerator_atom(&self.parent.inner.surfaces)
            .replace_multiple(
                self.parent
                    .inner
                    .surfaces
                    .iter()
                    .map(|s| symbolica::id::Replacement::new(s.atom(false), s.atom(true))),
            );
        TensorExpression::from_atom_interface(py, value, None)
    }
    /// Edge-energy substitutions used to evaluate this contribution's numerator.
    ///
    /// Keys are physical edge IDs; values replace Q(edge, cind(0)). Empty for
    /// topology-only scalar families. These sampling maps need not be on-shell cuts.
    ///
    /// Examples
    /// --------
    /// >>> substitutions = family.energy_map
    #[getter]
    fn energy_map(&self) -> BTreeMap<usize, PythonExpression> {
        self.term
            .energy_map
            .iter()
            .map(|(e, a)| (*e, PythonExpression { expr: a.clone() }))
            .collect()
    }
    /// Return this family's contribution, including the common energy prefactor.
    ///
    /// Examples
    /// --------
    /// Using the setup in the class example:
    ///
    /// >>> expression = family.to_expression()
    ///
    /// Parameters
    /// ----------
    /// expand_surfaces : bool, optional
    ///     Replace stored surface placeholders by affine energy definitions.
    ///     Defaults to True. On-shell energies remain symbolic in either case.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     The complete family contribution, including its evaluated numerator and energy factors.
    #[pyo3(signature=(*, expand_surfaces=true))]
    fn to_expression(
        &self,
        py: Python<'_>,
        expand_surfaces: bool,
    ) -> PyResult<Py<TensorExpression>> {
        let denominator = self.term.atom(&self.parent.inner.surfaces);
        TensorExpression::from_atom_interface(
            py,
            self.parent.expression(denominator, expand_surfaces).expr,
            None,
        )
    }
    /// Inspect only this family's contribution and surface circlings.
    ///
    /// Examples
    /// --------
    /// >>> html = family._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.parent.explorer_html(
            py,
            display::CffScope::Family {
                orientation: self.orientation,
                family: self.id,
            },
        )
    }
    /// Display this family in Marimo.
    ///
    /// Examples
    /// --------
    /// >>> view = family._display_()
    fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py
            .import("marimo")?
            .call_method1("iframe", (self._repr_html_(py)?,))?
            .unbind())
    }
}

/// Generation counts for a CFF representation.
///
/// Compare candidate energy flows with the acyclic flows and denominator
/// products retained in the result.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> result = diagram.integrate_energy(method="cff")
/// >>> report = result.report
/// >>> assert report.candidate_orientations >= report.acyclic_orientations
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffReport",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffReport {
    inner: CffReport,
    generalized: bool,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCffReport {
    /// Candidate search count, or None for generalized numerator generation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> assert report.candidate_orientations >= report.acyclic_orientations
    #[getter]
    fn candidate_orientations(&self) -> Option<usize> {
        (!self.generalized).then_some(self.inner.candidate_orientations)
    }
    /// Acyclic search count, or None for generalized numerator generation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> assert report.acyclic_orientations == len(result.orientations)
    #[getter]
    fn acyclic_orientations(&self) -> Option<usize> {
        (!self.generalized).then_some(self.inner.acyclic_orientations)
    }
    /// Return the total number of unfolded denominator terms.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> term_count = report.unfolded_terms
    #[getter]
    fn unfolded_terms(&self) -> usize {
        self.inner.unfolded_terms
    }
    /// Return the number of unique denominator surfaces in the result.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> assert report.interned_surfaces == len(result.surfaces)
    #[getter]
    fn interned_surfaces(&self) -> usize {
        self.inner.interned_surfaces
    }

    /// Return a concise constructor-style CFF generation summary.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> print(result.report)
    fn __repr__(&self) -> String {
        format!(
            "CffReport(candidate_orientations={}, acyclic_orientations={}, unfolded_terms={}, interned_surfaces={})",
            self.candidate_orientations()
                .map_or_else(|| "None".into(), |n| n.to_string()),
            self.acyclic_orientations()
                .map_or_else(|| "None".into(), |n| n.to_string()),
            self.inner.unfolded_terms,
            self.inner.interned_surfaces,
        )
    }

    /// Render CFF generation statistics as a compact HTML table.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(report)
    fn _repr_html_(&self) -> String {
        crate::display::record_html(
            "CffReport",
            "CFF report",
            &[
                (
                    "candidate orientations",
                    self.candidate_orientations()
                        .map_or_else(|| "—".into(), |n| n.to_string()),
                ),
                (
                    "acyclic orientations",
                    self.acyclic_orientations()
                        .map_or_else(|| "—".into(), |n| n.to_string()),
                ),
                ("unfolded terms", self.inner.unfolded_terms.to_string()),
                (
                    "interned surfaces",
                    self.inner.interned_surfaces.to_string(),
                ),
            ],
        )
    }

    /// Write a concise CFF summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffReport`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(report)
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

/// A loop-energy integral in the cross-free family (CFF) representation.
///
/// Create with ``diagram.integrate_energy(method="cff")``. The result groups
/// acyclic energy-flow orientations into cross-free families of denominators.
/// With :math:`E_e=\sqrt{\boldsymbol{q}_e^2+m_e^2}`, its scalar expression is
///
/// $$
/// I_{\mathrm{CFF}}=\left(\prod_e\frac{-1}{2E_e}\right)
/// \sum_o\sum_{F\in\mathcal{F}_o}\prod_{s\in F}\frac{1}{s}.
/// $$
///
/// Repeated surfaces appear with their multiplicity. Family and orientation
/// expressions include their energy prefactors, so they sum to the result.
/// The measure is :math:`d\ell^0/(2\pi i)` per loop, with contours closed below,
/// as in ``LtdRepresentation``. The scalar formula above applies without a
/// numerator. Supply ``numerator=...`` to ``integrate_energy`` to generate
/// bounded polynomial energy numerators, including quadratic and higher powers.
/// Each generalized family retains its own coefficient, numerator evaluation
/// map and denominator powers. Couplings, graph weights and the spatial measure
/// remain separate.
///
/// Both representations use ``gammalooprs::OSE(e)`` for :math:`E_e` and
/// ``gammalooprs::Q(e, spenso::cind(0))`` for external energies; ``e`` is a
/// physical diagram edge ID. ``on_shell_energies`` supplies the routed square
/// roots. Impose external momentum conservation when comparing CFF and LTD.
///
/// Displaying the object opens the orientation, family and surface explorer in
/// Marimo or an HTML notebook. ``len(result)`` counts denominator products.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> result = diagram.integrate_energy(method="cff")
/// >>> expression = result.to_expression()
/// >>> assert result.report.acyclic_orientations == len(result.orientations)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffRepresentation",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffRepresentation {
    inner: Arc<CffKernel>,
    numerator: Atom,
    energy_degree_bounds: BTreeMap<usize, usize>,
    generalized: bool,
    source: PyFeynmanDiagram,
    energy_edges: Vec<usize>,
    owner: Arc<()>,
    normalization: Atom,
    diagram: Arc<feynkit_graph::FeynmanDiagram>,
    drawing: Arc<OnceLock<String>>,
    pole_order: Option<usize>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCffRepresentation {
    /// Numerator supplied to energy integration; one for scalar generation.
    ///
    /// Examples
    /// --------
    /// >>> numerator = result.numerator
    #[getter]
    fn numerator(&self) -> PythonExpression {
        PythonExpression {
            expr: self.numerator.clone(),
        }
    }

    /// Inferred polynomial energy-degree bounds keyed by physical internal edge ID.
    ///
    /// Omitted edges have degree zero. Bounds are computed before applying routing.
    ///
    /// Examples
    /// --------
    /// >>> bounds = result.energy_degree_bounds
    #[getter]
    fn energy_degree_bounds(&self) -> BTreeMap<usize, usize> {
        self.energy_degree_bounds.clone()
    }

    /// Return generation statistics for this CFF result.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> assert result.report.acyclic_orientations == len(result.orientations)
    #[getter]
    fn report(&self) -> PyCffReport {
        PyCffReport {
            inner: self.inner.report,
            generalized: self.generalized,
        }
    }

    /// Return the acyclic energy-flow orientations in this result.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> products = [item.to_expression() for item in result.orientations]
    #[getter]
    fn orientations(&self) -> Vec<PyCffOrientation> {
        self.inner
            .orientations
            .iter()
            .cloned()
            .map(|inner| PyCffOrientation {
                inner,
                parent: self.clone(),
            })
            .collect()
    }

    /// Return all unique energy and H surfaces in this result.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> energies = [surface.to_expression() for surface in result.surfaces]
    #[getter]
    fn surfaces(&self) -> Vec<PyEnergySurface> {
        self.inner.surfaces.clone()
    }

    /// Return the scalar loop-energy integral as a Symbolica expression.
    ///
    /// Includes contour signs and on-shell energy factors for
    /// :math:`d\ell^0/(2\pi i)` per loop, with contours closed below. A supplied
    /// numerator is included; couplings, graph weights and the spatial measure are separate.
    /// ``expand_surfaces`` expands affine denominators while leaving
    /// ``gammalooprs::OSE(e)`` symbolic. Use ``on_shell_energies`` for their routed
    /// definitions or custom evaluation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the class example:
    ///
    /// >>> expression = result.to_expression()
    /// >>> compact = result.to_expression(expand_surfaces=False)
    ///
    /// Parameters
    /// ----------
    /// expand_surfaces : bool, optional
    ///     Replace surface placeholders by their affine energy definitions.
    ///     Defaults to True. False keeps placeholders local to this representation.
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Rank-zero tensor summing all orientation contributions and prefactors.
    #[pyo3(signature = (*, expand_surfaces=true))]
    fn to_expression(
        &self,
        py: Python<'_>,
        expand_surfaces: bool,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(
            py,
            self.expression(self.inner.to_atom(), expand_surfaces).expr,
            None,
        )
    }
    /// Originating diagram or selected subgraph, retaining its routing and owner.
    ///
    /// Examples
    /// --------
    /// >>> routing = result.diagram.loop_momentum_basis
    #[getter]
    fn diagram(&self) -> PyFeynmanDiagram {
        self.source.clone()
    }
    /// Shared symbolic energies and routed square roots, keyed by physical edge ID.
    ///
    /// Examples
    /// --------
    /// >>> definitions = [(e.symbol, e.to_expression()) for e in result.on_shell_energies.values()]
    #[getter]
    fn on_shell_energies(&self) -> PyResult<BTreeMap<usize, PyOnShellEnergy>> {
        PyOnShellEnergy::for_diagram(&self.diagram, self.energy_edges.iter().copied())
    }

    /// Group equivalent energy surfaces after identifying raised propagator edges.
    /// ``edge_representatives`` maps repeated edges to their canonical edge.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> groups = result.raised_surface_groups({3: 2})
    /// >>> [group.max_order for group in groups]
    ///
    /// Parameters
    /// ----------
    /// edge_representatives : dict[int, int], optional
    ///     Repeated propagator edge IDs mapped to their canonical representative.
    #[pyo3(signature = (edge_representatives=None))]
    fn raised_surface_groups(
        &self,
        edge_representatives: Option<BTreeMap<usize, usize>>,
    ) -> PyResult<Vec<PyCffSurfaceGroup>> {
        let representatives = edge_representatives.unwrap_or_default();
        let replacements: Vec<_> = representatives
            .iter()
            .map(|(edge, representative)| {
                symbolica::id::Replacement::new(
                    feynkit_cff::symbols::on_shell_atom(linnet::half_edge::involution::EdgeIndex(
                        *edge,
                    )),
                    feynkit_cff::symbols::on_shell_atom(linnet::half_edge::involution::EdgeIndex(
                        *representative,
                    )),
                )
            })
            .collect();
        let mut groups = BTreeMap::<Atom, Vec<usize>>::new();
        for (id, surface) in self
            .inner
            .surfaces
            .iter()
            .enumerate()
            .filter(|(_, s)| s.kind() == "E")
        {
            groups
                .entry(surface.atom(true).replace_multiple(&replacements))
                .or_default()
                .push(id);
        }
        Ok(groups
            .into_values()
            .filter_map(|ids| {
                let max_order = self
                    .inner
                    .orientations
                    .iter()
                    .flat_map(|o| &o.terms)
                    .map(|term| {
                        term.path
                            .iter()
                            .filter(|id| ids.contains(id))
                            .count()
                            .saturating_sub(
                                term.numerator_factors
                                    .iter()
                                    .filter(|id| ids.contains(id))
                                    .count(),
                            )
                    })
                    .max()
                    .unwrap_or(0);
                (max_order > 0).then(|| PyCffSurfaceGroup {
                    surfaces: ids
                        .iter()
                        .map(|id| self.inner.surfaces[*id].clone())
                        .collect(),
                    ids,
                    max_order,
                    owner: Arc::clone(&self.owner),
                })
            })
            .collect())
    }

    /// Return coefficients of each inverse surface power, indexed from order one.
    /// These are pole coefficients, before analytic residue derivatives.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> coefficients = result.pole_coefficients(result.raised_surface_groups()[0])
    /// >>> [coefficient.to_expression() for coefficient in coefficients]
    ///
    /// Parameters
    /// ----------
    /// group : CffSurfaceGroup
    ///     A raised-surface group belonging to this result.
    fn pole_coefficients(&self, group: &PyCffSurfaceGroup) -> PyResult<Vec<Self>> {
        self.validate_group(group)?;
        Ok((1..=group.max_order)
            .map(|order| {
                let mut result = self.clone();
                let inner = Arc::make_mut(&mut result.inner);
                for orientation in &mut inner.orientations {
                    orientation.terms.retain_mut(|term| {
                        let count = term.path.iter().filter(|id| group.ids.contains(id)).count();
                        let numerator_count = term
                            .numerator_factors
                            .iter()
                            .filter(|id| group.ids.contains(id))
                            .count();
                        if count != order + numerator_count {
                            return false;
                        }
                        term.path.retain(|id| !group.ids.contains(id));
                        term.numerator_factors.retain(|id| !group.ids.contains(id));
                        true
                    });
                    orientation.expression = orientation
                        .terms
                        .iter()
                        .fold(Atom::Zero, |sum, term| sum + term.atom(&inner.surfaces));
                }
                result.pole_order = Some(order);
                result
            })
            .collect())
    }

    /// Evaluate all pole-order contributions to a residue in an explicit variable.
    /// ``surface`` must be the group's energy surface expressed in that variable;
    /// ``coefficient`` is the complete remaining coefficient, including any factors
    /// whose derivatives must act. The supplied root is assumed to be a simple zero.
    ///
    /// Examples
    /// --------
    /// Using the setup in ``CffRepresentation``, illustrate a simple pole locally
    /// parameterized by ``surface=t`` with constant remaining coefficient:
    ///
    /// >>> from symbolica import E, S
    /// >>> t = S("t")
    /// >>> group = result.raised_surface_groups()[0]
    /// >>> residue = result.residue(group, variable=t, root=E("0"), surface=t, coefficient=E("1"))
    ///
    /// Parameters
    /// ----------
    /// group : CffSurfaceGroup
    ///     A raised-surface group belonging to this result.
    /// variable : Expression
    ///     Independent integration variable.
    /// root : Expression
    ///     Simple zero of the surface, independent of variable.
    /// surface : Expression
    ///     Energy surface expressed in the integration variable.
    /// coefficient : Expression
    ///     Complete remaining coefficient to differentiate.
    /// replacements : list[tuple[Expression, Expression]], optional
    ///     Route all energy dependence to the integration variable before differentiating.
    #[pyo3(signature = (group, *, variable, root, surface, coefficient, replacements=None))]
    #[allow(clippy::too_many_arguments)] // Keep Python's residue coordinates explicit keyword arguments.
    fn residue(
        &self,
        py: Python<'_>,
        group: &PyCffSurfaceGroup,
        variable: ConvertibleToExpression,
        root: ConvertibleToExpression,
        surface: ConvertibleToExpression,
        coefficient: ConvertibleToExpression,
        replacements: Option<Vec<(ConvertibleToExpression, ConvertibleToExpression)>>,
    ) -> PyResult<Py<TensorExpression>> {
        let variable = expression_variable(variable)?;
        let root = root.to_expression().expr;
        let surface = surface.to_expression().expr;
        let coefficient = coefficient.to_expression().expr;
        let replacements = replacements
            .unwrap_or_default()
            .into_iter()
            .map(|(from, to)| {
                symbolica::id::Replacement::new(from.to_expression().expr, to.to_expression().expr)
            })
            .collect::<Vec<_>>();
        let mut result = Atom::Zero;
        for (index, term) in self.pole_coefficients(group)?.iter().enumerate() {
            let remaining = (coefficient.clone()
                * term.expression(term.inner.to_atom(), true).expr)
                .replace_multiple(&replacements);
            result += SurfacePole {
                surface: surface.clone(),
                order: index + 1,
            }
            .residue(&remaining, variable, &root)
            .map_err(error::cff)?;
        }
        TensorExpression::from_atom_interface(py, result, None)
    }

    /// Return the number of unfolded denominator terms.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> denominator_term_count = len(result)
    fn __len__(&self) -> usize {
        self.inner.term_count()
    }

    /// Return a concise summary of the CFF expression and its surfaces.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> print(result)
    fn __repr__(&self) -> String {
        let summary = format!(
            "CffRepresentation(orientations={}, terms={}, surfaces={})",
            self.inner.orientations.len(),
            self.inner.term_count(),
            self.inner.surfaces.len(),
        );
        match self.pole_order {
            Some(order) => format!("{summary} [pole coefficient, order {order}]"),
            None => summary,
        }
    }

    /// Explore orientations, factored families and surface regions on the native graph.
    ///
    /// Arrowhead clicks select another retained orientation; they never change the result.
    /// Shift-click surface factors to inspect several circlings together. The displayed
    /// expression includes any supplied numerator and excludes the spatial measure.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(result)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.explorer_html(py, display::CffScope::Representation)
    }

    /// Return Marimo's interactive presentation in a script-enabled iframe.
    ///
    /// Marimo uses this hook automatically when displaying the result.
    /// Other notebook frontends use ``_repr_html_`` without requiring Marimo.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> presentation = result._display_()
    fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let html = self.explorer_html(py, display::CffScope::Representation)?;
        Ok(py
            .import("marimo")?
            .call_method1("iframe", (html,))?
            .unbind())
    }

    /// Write a summary with Symbolica's native expression formatting.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffRepresentation`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(result)
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
                "CffRepresentation(orientations={}, terms={}, surfaces={}, expression=",
                self.inner.orientations.len(),
                self.inner.term_count(),
                self.inner.surfaces.len(),
            ),),
        )?;
        self.expression(self.inner.to_atom(), false)
            ._repr_pretty_(pretty, false)?;
        pretty.call_method1("text", (")",))?;
        Ok(())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCrossFreeFamily>()?;
    module.add_class::<PyCffOrientation>()?;
    module.add_class::<PyCffReport>()?;
    module.add_class::<PyCffRepresentation>()?;
    module.add_class::<PyCffSurfaceGroup>()?;
    module.add_class::<PyCutPropagator>()?;
    Ok(())
}

impl PyCffRepresentation {
    fn expression(&self, denominator: Atom, expand_surfaces: bool) -> PythonExpression {
        let mut expr = denominator * &self.normalization;
        if expand_surfaces {
            expr = expr.replace_multiple(
                self.inner
                    .surfaces
                    .iter()
                    .map(|s| symbolica::id::Replacement::new(s.atom(false), s.atom(true))),
            );
        }
        PythonExpression { expr }
    }
    pub(crate) fn from_diagram(
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        max_orientations: Option<usize>,
        fixed_orientations: Option<BTreeMap<usize, bool>>,
        contracted_edges: Option<Vec<usize>>,
        initial_state_edges: Option<Vec<usize>>,
        numerator: Option<ConvertibleToExpression>,
    ) -> PyResult<Self> {
        if let Some(numerator) = numerator {
            diagram.require_complete()?;
            if max_orientations.is_some()
                || fixed_orientations.is_some()
                || contracted_edges.is_some()
                || initial_state_edges.is_some()
            {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "numerator-aware CFF currently requires an unconstrained complete diagram",
                ));
            }
            let source = diagram.clone();
            let diagram = Arc::clone(&diagram.inner);
            let numerator = numerator.to_expression().expr;
            return py.detach(move || {
                let (inner, normalization, energy_degree_bounds) =
                    CffKernel::from_numerator(&diagram, &numerator)?;
                let energy_edges = diagram
                    .edges()
                    .filter(|(_, _, e)| e.external.is_none())
                    .map(|(id, _, _)| id.0)
                    .collect();
                Ok(Self {
                    inner: Arc::new(inner),
                    source,
                    energy_edges,
                    owner: Arc::new(()),
                    normalization,
                    diagram,
                    drawing: Arc::new(OnceLock::new()),
                    pole_order: None,
                    numerator,
                    energy_degree_bounds,
                    generalized: true,
                })
            });
        }
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

        let selection = diagram.selection();
        Self::build(py, diagram, options, selection)
    }

    fn validate_group(&self, group: &PyCffSurfaceGroup) -> PyResult<()> {
        if !Arc::ptr_eq(&self.owner, &group.owner) {
            return Err(error::CffError::new_err(
                "surface group belongs to a different CFF result",
            ));
        }
        Ok(())
    }

    fn build(
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        options: CffOptions,
        selection: SuBitGraph,
    ) -> PyResult<Self> {
        let source = diagram.clone();
        let diagram = diagram.inner.clone();
        py.detach(move || {
            let graph = diagram.underlying();
            let initial_edges = graph
                .iter_edges()
                .filter_map(|(pair, edge, data)| {
                    (data.data.external.is_some() && matches!(pair, HedgePair::Paired { .. }))
                        .then_some(edge)
                })
                .collect::<BTreeSet<_>>();
            let mut energies = Vec::new();
            let mut energy_edges = Vec::new();
            for (pair, edge, data) in graph.iter_edges_of(&selection) {
                if let HedgePair::Paired { .. } = pair {
                    if initial_edges.contains(&edge) || data.data.is_dummy {
                        continue;
                    }
                    if !options.contracted_edges().contains(&edge) {
                        energy_edges.push(edge.0);
                        energies.push(feynkit_cff::symbols::on_shell_atom(edge));
                    }
                }
            }
            // CFF denominators use the source-frame sign per propagator.
            // dq0/(2*pi*i) contours need no additional phase or spatial measure.
            let normalization = CffExpression::inverse_energy_product(energies);
            let inner = diagram
                .build_cff_subgraph(&selection, options)
                .map_err(error::cff)?;
            Ok(Self {
                inner: Arc::new(CffKernel::from_topology(inner)),
                numerator: Atom::num(1),
                energy_degree_bounds: BTreeMap::new(),
                generalized: false,
                source,
                energy_edges,
                owner: Arc::new(()),
                normalization,
                diagram,
                drawing: Arc::new(OnceLock::new()),
                pole_order: None,
            })
        })
    }
}

/// Equivalent energy surfaces and their maximum simultaneous pole order.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> model = hep.Model.phi4()
/// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
/// >>> result = process.generate_diagrams(loops=1)
/// >>> diagram = result.diagrams[0]
/// >>> result = diagram.integrate_energy(method="cff")
/// >>> group = result.raised_surface_groups()[0]
/// >>> coefficients = result.pole_coefficients(group)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffSurfaceGroup",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffSurfaceGroup {
    ids: Vec<usize>,
    surfaces: Vec<PyEnergySurface>,
    max_order: usize,
    owner: Arc<()>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCffSurfaceGroup {
    /// Summarize the equivalent surfaces and maximum pole order.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(group)
    fn __repr__(&self) -> String {
        format!(
            "CffSurfaceGroup(surfaces={}, max_order={})",
            self.ids.len(),
            self.max_order()
        )
    }

    /// Render equivalent surface definitions and their maximum pole order.
    ///
    /// Examples
    /// --------
    /// >>> html = group._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        let definitions = self
            .surfaces()
            .iter()
            .map(|s| s.definition_html(py))
            .collect::<PyResult<Vec<_>>>()?
            .join("<br>");
        Ok(crate::display::record_html(
            "CffSurfaceGroup",
            "Equivalent surfaces",
            &[
                ("Definitions", definitions),
                ("Maximum pole order", self.max_order().to_string()),
            ],
        ))
    }
    /// Highest inverse surface power occurring on one CFF branch.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffSurfaceGroup`` class example:
    ///
    /// >>> highest_pole = group.max_order
    /// >>> assert highest_pole >= 1
    #[getter]
    fn max_order(&self) -> usize {
        self.max_order
    }
    /// Canonical energy surfaces identified by this group.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CffSurfaceGroup`` class example:
    ///
    /// >>> equivalent_surfaces = group.surfaces
    #[getter]
    fn surfaces(&self) -> Vec<PyEnergySurface> {
        self.surfaces.clone()
    }
}

fn expression_variable(value: ConvertibleToExpression) -> PyResult<Symbol> {
    match value.to_expression().expr.as_view() {
        AtomView::Var(variable) => Ok(variable.get_symbol()),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "integration variable must be a Symbolica symbol",
        )),
    }
}

/// An oriented generalized cut distribution for a possibly raised propagator.
/// ``orientation`` selects q0=+E or q0=-E; ``prescription`` is the sign of i0.
/// The default normalization is -2 pi i. The positive ``on_shell_energy`` includes
/// the mass. Raised cuts act on the entire coefficient, including the uncut factor.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hepkit as hep
/// >>> q0, energy = S("q0", "energy")
/// >>> cut = hep.CutPropagator(q0, energy, power=2)
/// >>> residue = cut.apply(q0**2, q0)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CutPropagator",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCutPropagator {
    pub(crate) inner: CutPropagator,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyCutPropagator {
    /// Construct a reciprocal-propagator cut with explicit orientation and i0 sign.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CutPropagator`` class example:
    ///
    /// >>> cut = hep.CutPropagator(q0, energy, power=3, orientation=-1)
    /// >>> residue = cut.apply(q0**3, q0)
    ///
    /// Parameters
    /// ----------
    /// energy : Expression
    ///     Time component of the propagator momentum.
    /// on_shell_energy : Expression
    ///     Positive square root of spatial momentum squared plus mass squared.
    /// power : int
    ///     Positive propagator power.
    /// orientation : int
    ///     Select the positive (+1) or negative (-1) energy root.
    /// prescription : int
    ///     Sign of the imaginary infinitesimal, +1 or -1.
    /// normalization : Expression, optional
    ///     Replacement prefactor, default -2 pi i.
    #[new]
    #[pyo3(signature = (energy, on_shell_energy, *, power=1, orientation=1, prescription=1, normalization=None))]
    fn new(
        energy: ConvertibleToExpression,
        on_shell_energy: ConvertibleToExpression,
        power: usize,
        orientation: i8,
        prescription: i8,
        normalization: Option<ConvertibleToExpression>,
    ) -> PyResult<Self> {
        let inner = CutPropagator {
            energy: energy.to_expression().expr,
            on_shell_energy: on_shell_energy.to_expression().expr,
            power,
            orientation,
            prescription,
            normalization: normalization.map_or_else(
                || -2 * Atom::var(Symbol::PI) * Atom::i(),
                |value| value.to_expression().expr,
            ),
        };
        inner.validate().map_err(error::cff)?;
        Ok(Self { inner })
    }

    /// Emit a covariant CutDelta or energy-space generalized Delta expression.
    /// Delta(n,x) acts as f^(n-1)(0)/(n-1)!; it is not the usual delta derivative.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CutPropagator`` class example:
    ///
    /// >>> cut.to_expression(covariant=False)
    ///
    /// Parameters
    /// ----------
    /// covariant : bool
    ///     Emit a mass-shell CutDelta; otherwise show the energy-space Delta.
    #[pyo3(signature = (*, covariant=true))]
    fn to_expression(&self, covariant: bool) -> PyResult<PythonExpression> {
        Ok(PythonExpression {
            expr: if covariant {
                self.inner.to_atom()
            } else {
                self.inner.to_energy_atom()
            }
            .map_err(error::cff)?,
        })
    }

    /// Apply the complete derivative action and then set the energy on shell.
    /// Route the integrand to this independent energy variable before applying.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``CutPropagator`` class example:
    ///
    /// >>> cut.apply(q0**2, q0)
    ///
    /// Parameters
    /// ----------
    /// coefficient : Expression
    ///     All remaining factors of the integrand.
    /// variable : Expression
    ///     Independent energy symbol, equal to this cut's energy.
    fn apply(
        &self,
        coefficient: ConvertibleToExpression,
        variable: ConvertibleToExpression,
    ) -> PyResult<PythonExpression> {
        let variable = expression_variable(variable)?;
        Ok(PythonExpression {
            expr: self
                .inner
                .apply(&coefficient.to_expression().expr, variable)
                .map_err(error::cff)?,
        })
    }
}
