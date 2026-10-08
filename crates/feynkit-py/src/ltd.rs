mod display;

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    sync::{Arc, OnceLock},
};

use feynkit_cff::generalized::{
    Generate3DExpressionOptions, GeneratedThreeDExpression, ParsedGraph, RepresentationMode,
    generate_3d_expression, graph_io::EnergyEdgeIndexMap,
};
use feynkit_graph::FeynmanDiagram;
use pyo3::prelude::*;
use spynso3::expression::TensorExpression;
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, AtomCore},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{
    energy::{PyEnergySurface, PyOnShellEnergy, PySurfaceFactor},
    error,
    graph::PyFeynmanDiagram,
};
use feynkit_cff::generalized::{LinearEnergyExpr, OrientationID};
use linnet::half_edge::involution::{EdgeIndex, Orientation};

/// A scalar loop-energy integral in the loop-tree duality (LTD) representation.
///
/// Create with ``diagram.integrate_energy(method="ltd")``. Each residue puts
/// edges complementary to a spanning tree on shell at
/// :math:`q_e^0=\sigma_e E_e`, where :math:`\sigma_e\in\{-1,+1\}` and
/// :math:`E_e=\sqrt{\boldsymbol{q}_e^2+m_e^2}`. The ordered loop-momentum
/// basis and successive contours determine the retained pole assignments and
/// coefficients. A spanning tree alone does not determine the signs.
///
/// ``to_expression()`` sums the residues, including signed coefficients and
/// on-shell factors. The measure is :math:`d\ell^0/(2\pi i)` per loop, with
/// contours closed below, as in ``CffRepresentation``. Numerators, couplings,
/// graph weights and the spatial integration measure are excluded. Repeated
/// propagators are handled by higher-order residues.
///
/// Both representations use ``gammalooprs::OSE(e)`` for :math:`E_e` and
/// ``gammalooprs::Q(e, spenso::cind(0))`` for external energies; ``e`` is a
/// physical diagram edge ID. ``on_shell_energies`` supplies the routed square
/// roots. Impose external momentum conservation when comparing CFF and LTD.
///
/// Displaying the object opens the residue and surface-pair explorer in Marimo
/// or an HTML notebook. ``len(ltd)`` counts retained signed residues.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> ltd = diagram.integrate_energy(method="ltd")
/// >>> expression = ltd.to_expression()
/// >>> assert ltd.report.residues == len(ltd.residues)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LtdRepresentation",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLtdRepresentation {
    inner: Arc<GeneratedThreeDExpression>,
    parsed: Arc<ParsedGraph>,
    indices: Arc<EnergyEdgeIndexMap>,
    diagram: Arc<FeynmanDiagram>,
    source: PyFeynmanDiagram,
    drawing: Arc<OnceLock<String>>,
}

impl PyLtdRepresentation {
    fn physical_energy(&self, expression: &LinearEnergyExpr) -> LinearEnergyExpr {
        expression
            .clone()
            .remap_energy_edges(&self.indices.internal, &self.indices.external)
    }
    fn residue_expression(&self, id: usize, expand_surfaces: bool) -> PythonExpression {
        let mut expr = self.inner.expression.orientations[OrientationID(id)].to_atom();
        if expand_surfaces {
            expr = self
                .inner
                .expression
                .surfaces
                .substitute_energies(&expr, &[]);
        }
        let replacements = self
            .indices
            .internal
            .iter()
            .map(|(local, physical)| {
                symbolica::id::Replacement::new(
                    feynkit_cff::symbols::on_shell_atom(EdgeIndex(*local)),
                    feynkit_cff::symbols::on_shell_atom(EdgeIndex(*physical)),
                )
            })
            .chain(self.indices.external.iter().map(|(local, physical)| {
                symbolica::id::Replacement::new(
                    feynkit_cff::symbols::external_energy_atom(EdgeIndex(*local)),
                    feynkit_cff::symbols::external_energy_atom(EdgeIndex(*physical)),
                )
            }))
            .collect::<Vec<_>>();
        PythonExpression {
            expr: expr.replace_multiple(replacements),
        }
    }
    pub(crate) fn from_diagram(py: Python<'_>, diagram: &PyFeynmanDiagram) -> PyResult<Self> {
        diagram.require_complete()?;
        let source = diagram.clone();
        let diagram = Arc::clone(&diagram.inner);
        py.detach(move || {
            let crate::energy::integration::EnergyGraph { parsed, indices } =
                crate::energy::integration::EnergyGraph::new(&diagram)?;
            let inner = generate_3d_expression(
                &parsed,
                &Generate3DExpressionOptions {
                    representation: RepresentationMode::Ltd,
                    energy_degree_bounds: Some(Vec::new()),
                    ..Default::default()
                },
            )
            .map_err(|e| error::FeynkitError::new_err(e.to_string()))?;
            Ok(Self {
                inner: Arc::new(inner),
                parsed: Arc::new(parsed),
                indices: Arc::new(indices),
                diagram,
                source,
                drawing: Arc::new(OnceLock::new()),
            })
        })
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyLtdRepresentation {
    /// Signed residues with tree, pole-sign and surface-pair inspection.
    ///
    /// Examples
    /// --------
    /// >>> residue = ltd.residues[0]
    #[getter]
    fn residues(&self) -> Vec<PyLtdResidue> {
        (0..self.inner.expression.orientations.len())
            .map(|id| PyLtdResidue {
                parent: self.clone(),
                id,
            })
            .collect()
    }
    /// All stored affine surfaces, with physical diagram edge IDs.
    ///
    /// Examples
    /// --------
    /// >>> definitions = [s.to_expression() for s in ltd.surfaces]
    #[getter]
    fn surfaces(&self) -> Vec<PyEnergySurface> {
        self.inner
            .expression
            .surfaces
            .linear_surface_cache
            .iter_enumerated()
            .map(|(id, s)| {
                PyEnergySurface::from_linear(id, s.kind, &self.physical_energy(&s.expression))
                    .with_provenance(s.origin, s.numerator_only)
            })
            .collect()
    }
    /// Originating diagram, retaining the ordered momentum routing and owner.
    ///
    /// Examples
    /// --------
    /// >>> routing = ltd.diagram.loop_momentum_basis
    #[getter]
    fn diagram(&self) -> PyFeynmanDiagram {
        self.source.clone()
    }
    /// Shared symbolic energies and routed square roots, keyed by physical edge ID.
    ///
    /// Examples
    /// --------
    /// >>> definitions = [(e.symbol, e.to_expression()) for e in ltd.on_shell_energies.values()]
    #[getter]
    fn on_shell_energies(&self) -> PyResult<BTreeMap<usize, PyOnShellEnergy>> {
        PyOnShellEnergy::for_diagram(&self.diagram, self.indices.internal.values().copied())
    }
    /// Counts of retained residues, distinct trees, unfolded terms and surfaces.
    ///
    /// Examples
    /// --------
    /// >>> assert ltd.report.residues == len(ltd.residues)
    #[getter]
    fn report(&self) -> PyLtdReport {
        PyLtdReport {
            residues: self.inner.expression.orientations.len(),
            trees: self
                .residues()
                .iter()
                .map(|r| r.cut_edges())
                .collect::<BTreeSet<_>>()
                .len(),
            unfolded_terms: self.inner.expression.num_unfolded_terms(),
            interned_surfaces: self.inner.expression.surfaces.linear_surface_cache.len(),
        }
    }

    /// Return the scalar loop-energy integral as a Symbolica expression.
    ///
    /// Includes contour signs and on-shell energy factors for
    /// :math:`d\ell^0/(2\pi i)` per loop, with contours closed below. Numerators,
    /// couplings, graph weights and the spatial measure are excluded.
    /// ``expand_surfaces`` expands affine denominators while leaving
    /// ``gammalooprs::OSE(e)`` symbolic. Use ``on_shell_energies`` for their routed
    /// definitions or custom evaluation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the class example:
    ///
    /// >>> expression = ltd.to_expression()
    /// >>> compact = ltd.to_expression(expand_surfaces=False)
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
    ///     Rank-zero tensor summing all residue contributions and prefactors.
    #[pyo3(signature = (*, expand_surfaces=true))]
    fn to_expression(
        &self,
        py: Python<'_>,
        expand_surfaces: bool,
    ) -> PyResult<Py<TensorExpression>> {
        let expression = self.residues().iter().fold(Atom::Zero, |sum, r| {
            sum + self.residue_expression(r.id, expand_surfaces).expr
        });
        TensorExpression::from_atom_interface(py, expression, None)
    }

    /// Return the number of signed residues.
    ///
    /// Examples
    /// --------
    /// Using the setup in ``LtdRepresentation``:
    ///
    /// >>> residue_count = len(ltd)
    fn __len__(&self) -> usize {
        self.inner.expression.orientations.len()
    }

    /// Summarize the native LTD result.
    ///
    /// Examples
    /// --------
    /// Using the setup in ``LtdRepresentation``:
    ///
    /// >>> print(ltd)
    fn __repr__(&self) -> String {
        format!(
            "LtdRepresentation(residues={}, terms={}, surfaces={})",
            self.__len__(),
            self.inner.expression.num_unfolded_terms(),
            self.inner.expression.surfaces.linear_surface_cache.len()
        )
    }

    /// Explore native residues, surface pairs, and region-relative pole signs.
    ///
    /// Click cuts to navigate residues, tree edges to switch their pair, and
    /// hover cuts to trace the cycle fixing their signs in the ordered routing.
    ///
    /// Examples
    /// --------
    /// Using the setup in ``LtdRepresentation``:
    ///
    /// >>> from IPython.display import display
    /// >>> display(ltd)
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.explorer_html(py, None)
    }

    /// Embed the explorer in Marimo's script-enabled notebook frame.
    ///
    /// Examples
    /// --------
    /// Using the setup in ``LtdRepresentation``:
    ///
    /// >>> presentation = ltd._display_()
    fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py
            .import("marimo")?
            .call_method1("iframe", (self.explorer_html(py, None)?,))?
            .unbind())
    }
}

/// One retained signed residue in an LTD representation.
///
/// The cut edges obey :math:`q_e^0=\sigma_e E_e`, with signs relative to the
/// stored momentum arrows. The diagram's ordered routing and contours determine
/// these assignments; several residues may share a spanning tree. Higher-order
/// poles can contribute several additive terms to one residue.
///
/// ``energy_map`` gives every internal energy at the residue. Substituting these
/// energies into an uncut propagator gives the ``surface_pairs`` factors
/// :math:`q_e^0-E_e` and :math:`q_e^0+E_e`.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> ltd = diagram.integrate_energy(method="ltd")
/// >>> residue = ltd.residues[0]
/// >>> signs = residue.pole_signs
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LtdResidue",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLtdResidue {
    parent: PyLtdRepresentation,
    id: usize,
}
impl PyLtdResidue {
    fn local_cuts(&self) -> BTreeSet<usize> {
        self.parent.inner.expression.orientations[OrientationID(self.id)]
            .variants
            .iter()
            .flat_map(|v| v.half_edges.iter().map(|e| e.0))
            .collect()
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyLtdResidue {
    /// Summarize this residue's cut edges and pole signs.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(residue)
    fn __repr__(&self) -> String {
        format!(
            "LtdResidue(id={}, pole_signs={:?})",
            self.id,
            self.pole_signs()
        )
    }
    /// Stable residue index within this representation.
    ///
    /// Examples
    /// --------
    /// >>> index = residue.id
    #[getter]
    fn id(&self) -> usize {
        self.id
    }
    /// Physical IDs of the on-shell edges complementary to the spanning tree.
    ///
    /// Examples
    /// --------
    /// >>> cuts = residue.cut_edges
    #[getter]
    fn cut_edges(&self) -> Vec<usize> {
        self.local_cuts()
            .iter()
            .map(|e| self.parent.indices.internal[e])
            .collect()
    }
    /// Physical IDs of uncut internal edges in the spanning tree.
    ///
    /// Examples
    /// --------
    /// >>> tree = residue.tree_edges
    #[getter]
    fn tree_edges(&self) -> Vec<usize> {
        let cuts = self.local_cuts();
        self.parent
            .indices
            .internal
            .iter()
            .filter_map(|(e, p)| (!cuts.contains(e)).then_some(*p))
            .collect()
    }
    /// Cut-edge pole signs relative to the diagram's stored momentum arrows.
    ///
    /// Each entry fixes :math:`q_e^0=\sigma_e E_e`, with
    /// :math:`\sigma_e\in\{-1,+1\}`. The ordered routing and successive energy
    /// contours determine these signs for the entire residue. Selecting another
    /// surface factor does not change them.
    ///
    /// Examples
    /// --------
    /// >>> assignments = residue.pole_signs
    ///
    /// Returns
    /// -------
    /// dict[int, int]
    ///     Physical cut-edge IDs mapped to their pole signs.
    #[getter]
    fn pole_signs(&self) -> BTreeMap<usize, i8> {
        let r = &self.parent.inner.expression.orientations[OrientationID(self.id)];
        self.local_cuts()
            .iter()
            .map(|e| {
                (
                    self.parent.indices.internal[e],
                    match r.data.orientation[EdgeIndex(*e)] {
                        Orientation::Default => 1,
                        Orientation::Reversed => -1,
                        Orientation::Undirected => unreachable!("retained cut has a pole sign"),
                    },
                )
            })
            .collect()
    }
    /// Internal energies :math:`q_e^0` evaluated at this residue.
    ///
    /// Examples
    /// --------
    /// >>> energies = residue.energy_map
    ///
    /// Returns
    /// -------
    /// dict[int, Expression]
    ///     Physical internal-edge IDs mapped to affine combinations of ``OSE`` and
    ///     external ``Q`` symbols. Cut entries equal :math:`\sigma_e E_e`.
    #[getter]
    fn energy_map(&self) -> BTreeMap<usize, PythonExpression> {
        self.parent.inner.expression.orientations[OrientationID(self.id)]
            .edge_energy_map
            .iter()
            .enumerate()
            .map(|(e, energy)| {
                (
                    self.parent.indices.internal[&e],
                    PythonExpression {
                        expr: self.parent.physical_energy(energy).to_atom(&[]),
                    },
                )
            })
            .collect()
    }
    /// The two factors of each surviving uncut propagator.
    ///
    /// Each pair preserves :math:`q_e^0-E_e` and :math:`q_e^0+E_e` after evaluating
    /// :math:`q_e^0` at the residue. Each factor's sign relates it to its stored
    /// surface; it does not change the cut-edge pole assignment.
    ///
    /// Examples
    /// --------
    /// >>> pairs = residue.surface_pairs
    ///
    /// Returns
    /// -------
    /// dict[int, SurfacePair]
    ///     Physical tree-edge IDs mapped to their ordered minus/plus factors.
    #[getter]
    fn surface_pairs(&self) -> BTreeMap<usize, PySurfacePair> {
        let cache = &self.parent.inner.expression.surfaces.linear_surface_cache;
        let lookup: HashMap<_, _> = cache
            .iter_enumerated()
            .map(|(id, s)| (s.expression.clone(), id))
            .collect();
        let cuts = self.local_cuts();
        self.parent.inner.expression.orientations[OrientationID(self.id)]
            .edge_energy_map
            .iter()
            .enumerate()
            .filter(|(e, _)| !cuts.contains(e))
            .filter_map(|(edge, energy)| {
                let factor = |sign| {
                    let raw =
                        (energy.clone() + LinearEnergyExpr::ose(EdgeIndex(edge), sign)).canonical();
                    let (id, sign) = lookup
                        .get(&raw)
                        .map(|id| (*id, 1))
                        .or_else(|| lookup.get(&(-raw)).map(|id| (*id, -1)))?;
                    let stored = &cache[id];
                    Some(PySurfaceFactor {
                        surface: PyEnergySurface::from_linear(
                            id,
                            stored.kind,
                            &self.parent.physical_energy(&stored.expression),
                        ),
                        sign,
                        power: 1,
                    })
                };
                Some((
                    self.parent.indices.internal[&edge],
                    PySurfacePair {
                        edge_id: self.parent.indices.internal[&edge],
                        minus: factor(-1)?,
                        plus: factor(1)?,
                    },
                ))
            })
            .collect()
    }
    /// Return this residue's contribution, including its signed coefficient.
    ///
    /// Examples
    /// --------
    /// Using the setup in the class example:
    ///
    /// >>> expression = residue.to_expression()
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
    ///     The scalar residue, including multiplicities and on-shell energy factors.
    #[pyo3(signature=(*, expand_surfaces=true))]
    fn to_expression(
        &self,
        py: Python<'_>,
        expand_surfaces: bool,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(
            py,
            self.parent
                .residue_expression(self.id, expand_surfaces)
                .expr,
            None,
        )
    }
    /// Inspect this residue's expression, fixed cuts and surface pairs.
    ///
    /// Examples
    /// --------
    /// >>> html = residue._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.parent.explorer_html(py, Some(self.id))
    }
    /// Display this residue in Marimo.
    ///
    /// Examples
    /// --------
    /// >>> view = residue._display_()
    fn _display_(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py
            .import("marimo")?
            .call_method1("iframe", (self._repr_html_(py)?,))?
            .unbind())
    }
}

/// The two linear factors of an uncut propagator at an LTD residue.
///
/// ``minus`` represents :math:`q_e^0-E_e` and ``plus`` represents
/// :math:`q_e^0+E_e`, with :math:`q_e^0` taken from ``residue.energy_map``.
/// Their product is :math:`(q_e^0)^2-E_e^2`. Each ``SurfaceFactor`` retains
/// its sign relative to the stored affine surface. Switching factors changes
/// only the coefficient of :math:`E_e`, not the cut-edge pole assignments.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> ltd = diagram.integrate_energy(method="ltd")
/// >>> residue = ltd.residues[0]
/// >>> pair = next(iter(residue.surface_pairs.values()))
/// >>> denominator = pair.minus.to_expression() * pair.plus.to_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SurfacePair",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PySurfacePair {
    edge_id: usize,
    minus: PySurfaceFactor,
    plus: PySurfaceFactor,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySurfacePair {
    /// Identify both denominator factors for this tree edge.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(pair)
    fn __repr__(&self) -> String {
        format!(
            "SurfacePair(edge={}, minus={}, plus={})",
            self.edge_id,
            self.minus.atom(false),
            self.plus.atom(false)
        )
    }

    /// Render the minus/plus factors evaluated at the owning residue.
    ///
    /// Examples
    /// --------
    /// >>> html = pair._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(crate::display::record_html(
            "SurfacePair",
            &format!("Tree edge {}", self.edge_id),
            &[
                ("q⁰ − E", self.minus.definition_html(py)?),
                ("q⁰ + E", self.plus.definition_html(py)?),
            ],
        ))
    }
    /// The :math:`q_e^0-E_e` denominator factor evaluated at this residue.
    ///
    /// Examples
    /// --------
    /// >>> factor = pair.minus
    ///
    /// Returns
    /// -------
    /// SurfaceFactor
    ///     The stored affine surface together with its occurrence sign.
    #[getter]
    fn minus(&self) -> PySurfaceFactor {
        self.minus.clone()
    }
    /// The :math:`q_e^0+E_e` denominator factor evaluated at this residue.
    ///
    /// Examples
    /// --------
    /// >>> factor = pair.plus
    ///
    /// Returns
    /// -------
    /// SurfaceFactor
    ///     The stored affine surface together with its occurrence sign.
    #[getter]
    fn plus(&self) -> PySurfaceFactor {
        self.plus.clone()
    }
}

/// Generation counts for an LTD representation.
///
/// Distinguish retained pole assignments from spanning trees: several residues
/// may share a tree, and a higher-order residue may contain several terms.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> ltd = diagram.integrate_energy(method="ltd")
/// >>> report = ltd.report
/// >>> assert report.residues == len(ltd)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LtdReport",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLtdReport {
    residues: usize,
    trees: usize,
    unfolded_terms: usize,
    interned_surfaces: usize,
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyLtdReport {
    /// Summarize the retained residues, spanning trees, terms and surfaces.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(report)
    fn __repr__(&self) -> String {
        format!(
            "LtdReport(residues={}, trees={}, unfolded_terms={}, interned_surfaces={})",
            self.residues, self.trees, self.unfolded_terms, self.interned_surfaces
        )
    }

    /// Render LTD generation statistics as a compact table.
    ///
    /// Examples
    /// --------
    /// >>> html = report._repr_html_()
    fn _repr_html_(&self) -> String {
        crate::display::record_html(
            "LtdReport",
            "LTD report",
            &[
                ("residues", self.residues.to_string()),
                ("spanning trees", self.trees.to_string()),
                ("unfolded terms", self.unfolded_terms.to_string()),
                ("interned surfaces", self.interned_surfaces.to_string()),
            ],
        )
    }
    /// Number of retained signed pole assignments.
    ///
    /// Examples
    /// --------
    /// >>> count = report.residues
    #[getter]
    fn residues(&self) -> usize {
        self.residues
    }
    /// Number of distinct spanning trees represented by the residues.
    ///
    /// Examples
    /// --------
    /// >>> count = report.trees
    #[getter]
    fn trees(&self) -> usize {
        self.trees
    }
    /// Number of additive denominator products, including higher-order residues.
    ///
    /// Examples
    /// --------
    /// >>> count = report.unfolded_terms
    #[getter]
    fn unfolded_terms(&self) -> usize {
        self.unfolded_terms
    }
    /// Number of distinct stored affine energy surfaces.
    ///
    /// Examples
    /// --------
    /// >>> count = report.interned_surfaces
    #[getter]
    fn interned_surfaces(&self) -> usize {
        self.interned_surfaces
    }
}
