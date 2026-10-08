//! Shared symbolic energy definitions and denominator inspection for CFF and LTD.
pub(crate) mod integration;

use std::collections::BTreeMap;

use feynkit_cff::generalized::{LinearEnergyExpr, LinearSurfaceID, LinearSurfaceKind};
use feynkit_cff::{Surface, SurfaceCache, SurfaceId};
use feynkit_graph::FeynmanDiagram;
use linnet::half_edge::involution::EdgeIndex;
use pyo3::prelude::*;
#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use spynso3::expression::TensorExpression;
use symbolica::{
    api::python::PythonExpression,
    atom::{Atom, AtomCore},
    symbol,
};

use crate::{
    display::{expression_html, record_html},
    error,
};

/// A symbolic on-shell energy shared by CFF and LTD.
///
/// ``symbol`` is ``gammalooprs::OSE(e)``, with physical diagram edge ID ``e``.
/// ``to_expression()`` returns the positive-root definition
///
/// $$
/// E_e=\sqrt{\boldsymbol{q}_e^2+m_e^2}.
/// $$
///
/// The spatial momentum is routed through the diagram's loop-momentum basis.
/// Its ``K(i, cind(a))`` and ``P(i, cind(a))`` coordinates use basis positions
/// ``i`` and Cartesian components ``a`` equal to 1, 2 or 3. These positions are
/// not physical edge IDs. Surface expansion leaves ``OSE`` symbolic; substitute
/// this definition for spatial evaluation, or supply your own energy evaluator.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> rep = diagram.integrate_energy(method="cff")
/// >>> energy = next(iter(rep.on_shell_energies.values()))
/// >>> replacement = (energy.symbol, energy.to_expression())
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "OnShellEnergy",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyOnShellEnergy {
    edge_id: usize,
    definition: Atom,
}
impl PyOnShellEnergy {
    pub(crate) fn for_diagram(
        diagram: &FeynmanDiagram,
        edges: impl IntoIterator<Item = usize>,
    ) -> PyResult<BTreeMap<usize, Self>> {
        let basis = diagram.loop_momentum_basis();
        edges
            .into_iter()
            .map(|edge_id| {
                let edge = diagram
                    .edges()
                    .find(|(id, _, _)| id.0 == edge_id)
                    .ok_or_else(|| error::DiagramError::new_err("on-shell energy edge is absent"))?
                    .2;
                let mass = diagram
                    .model()
                    .particle_by_id(edge.particle)
                    .map_err(error::model)?
                    .symbolic_mass(diagram.model());
                let squared = (1..=3).fold(mass.pow(2), |sum, axis| {
                    let component = feynkit_graph::symbols::momentum()
                        .call((edge_id, symbol!("spenso::cind").call(axis)));
                    sum + basis.route_expression(&component).pow(2)
                });
                Ok((
                    edge_id,
                    Self {
                        edge_id,
                        definition: squared.sqrt(),
                    },
                ))
            })
            .collect()
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyOnShellEnergy {
    /// Physical diagram edge owning this energy.
    ///
    /// Examples
    /// --------
    /// >>> edge_id = energy.edge_id
    #[getter]
    fn edge_id(&self) -> usize {
        self.edge_id
    }
    /// Canonical ``gammalooprs::OSE(edge_id)``, shared between CFF and LTD.
    ///
    /// Examples
    /// --------
    /// >>> symbol = energy.symbol
    #[getter]
    pub(crate) fn symbol(&self) -> PythonExpression {
        PythonExpression {
            expr: feynkit_cff::symbols::on_shell_atom(EdgeIndex(self.edge_id)),
        }
    }
    /// Return the routed definition :math:`\sqrt{\boldsymbol{q}_e^2+m_e^2}`.
    ///
    /// Examples
    /// --------
    /// >>> definition = energy.to_expression()
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Positive square root in the diagram's spatial K/P basis coordinates.
    fn to_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.definition.clone(), None)
    }
    /// Display the symbolic energy and its definition.
    ///
    /// Examples
    /// --------
    /// >>> print(energy)
    fn __repr__(&self) -> String {
        format!("{} = {}", self.symbol().expr, self.definition)
    }

    /// Render the energy symbol with its expandable routed square root.
    ///
    /// Examples
    /// --------
    /// >>> html = energy._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(record_html(
            "OnShellEnergy",
            &format!("Edge {}", self.edge_id),
            &[
                ("Energy", expression_html(py, self.symbol())?),
                (
                    "Definition",
                    format!(
                        "<details><summary>Routed square root</summary>{}</details>",
                        expression_html(
                            py,
                            PythonExpression {
                                expr: self.definition.clone()
                            }
                        )?
                    ),
                ),
            ],
        ))
    }
}

/// A stored affine energy denominator in a CFF or LTD representation.
///
/// Its definition has the form
///
/// $$
/// s=c_0+\sum_e c_e E_e+\sum_a d_a Q_a^0.
/// $$
///
/// ``energy_coefficients`` stores :math:`c_e`, ``external_shift`` stores
/// :math:`d_a`, and ``constant`` stores :math:`c_0`. Both coefficient mappings
/// use physical diagram edge IDs. Nonzero energy coefficients of one sign
/// define an E-surface; mixed signs define an H-surface.
///
/// ``symbol`` and ``index`` are local to the representation and surface category.
/// The ``OSE`` symbols inside the definition are shared between representations.
/// A ``SurfaceFactor`` records an occurrence's sign and multiplicity without
/// changing this stored definition.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> rep = diagram.integrate_energy(method="cff")
/// >>> surface = rep.surfaces[0]
/// >>> expression = surface.to_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EnergySurface",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEnergySurface {
    pub(crate) cff_id: Option<SurfaceId>,
    index: usize,
    kind: &'static str,
    symbol: Atom,
    energy_coefficients: BTreeMap<usize, Atom>,
    external_shift: BTreeMap<usize, Atom>,
    constant: Atom,
    vertices: Vec<usize>,
    origin: &'static str,
    numerator_only: bool,
}
impl PyEnergySurface {
    pub(crate) fn definition_html(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "{} = {}",
            expression_html(py, self.symbol())?,
            expression_html(
                py,
                PythonExpression {
                    expr: self.atom(true)
                }
            )?
        ))
    }
    pub(crate) fn from_cff(id: SurfaceId, cache: &SurfaceCache) -> Self {
        let surface = cache.get(id).expect("surface belongs to the CFF arena");
        let (index, kind, symbol) = match id {
            SurfaceId::Energy(id) => (id.index(), "E", Atom::from(id)),
            SurfaceId::H(id) => (id.index(), "H", Atom::from(id)),
            _ => unreachable!("unit and infinite are not energy surfaces"),
        };
        let mut energy_coefficients = BTreeMap::new();
        let (external_shift, vertices) = match surface {
            Surface::Energy(s) => {
                for e in s.energies {
                    *energy_coefficients.entry(e.index()).or_insert(Atom::Zero) += Atom::num(1);
                }
                (s.external_shift, s.vertex_set)
            }
            Surface::H(s) => {
                for e in s.positive_energies {
                    *energy_coefficients.entry(e.index()).or_insert(Atom::Zero) += Atom::num(1);
                }
                for e in s.negative_energies {
                    *energy_coefficients.entry(e.index()).or_insert(Atom::Zero) -= Atom::num(1);
                }
                (s.external_shift, s.vertex_set)
            }
            _ => unreachable!("unit and infinite are not energy surfaces"),
        };
        Self {
            origin: "physical",
            numerator_only: false,
            cff_id: Some(id),
            index,
            kind,
            symbol,
            energy_coefficients,
            external_shift: external_shift
                .iter()
                .map(|(e, c)| (e.index(), Atom::num(*c)))
                .collect(),
            constant: Atom::Zero,
            vertices: vertices.iter().map(|v| v.index()).collect(),
        }
    }
    pub(crate) fn from_linear(
        id: LinearSurfaceID,
        kind: LinearSurfaceKind,
        expression: &LinearEnergyExpr,
    ) -> Self {
        Self {
            origin: "physical",
            numerator_only: false,
            cff_id: None,
            index: id.0,
            kind: if kind == LinearSurfaceKind::Hsurface {
                "H"
            } else {
                "E"
            },
            symbol: Atom::from(id),
            energy_coefficients: expression
                .internal_terms
                .iter()
                .map(|(e, c)| (e.0, Atom::num(c.clone())))
                .collect(),
            external_shift: expression
                .external_terms
                .iter()
                .map(|(e, c)| (e.0, Atom::num(c.clone())))
                .collect(),
            constant: Atom::num(expression.constant.clone()),
            vertices: Vec::new(),
        }
    }
    pub(crate) fn with_provenance(
        mut self,
        origin: feynkit_cff::generalized::surface::SurfaceOrigin,
        numerator_only: bool,
    ) -> Self {
        self.origin = match origin {
            feynkit_cff::generalized::surface::SurfaceOrigin::Physical => "physical",
            feynkit_cff::generalized::surface::SurfaceOrigin::Helper => "helper",
        };
        self.numerator_only = numerator_only;
        self
    }
    pub(crate) fn with_vertices(mut self, vertices: Vec<usize>) -> Self {
        self.vertices = vertices;
        self
    }
    pub(crate) fn atom(&self, expand: bool) -> Atom {
        if !expand {
            return self.symbol.clone();
        }
        self.energy_coefficients
            .iter()
            .fold(self.constant.clone(), |sum, (e, c)| {
                sum + c * feynkit_cff::symbols::on_shell_atom(EdgeIndex(*e))
            })
            + self.external_shift.iter().fold(Atom::Zero, |sum, (e, c)| {
                sum + c * feynkit_cff::symbols::external_energy_atom(EdgeIndex(*e))
            })
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PyEnergySurface {
    /// Whether this is a physical surface or an algebraic helper ("physical" or "helper").
    ///
    /// Helper factors in generalized CFF need not correspond to a graph region.
    ///
    /// Examples
    /// --------
    /// >>> provenance = surface.origin
    #[getter]
    pub(crate) fn origin(&self) -> &'static str {
        self.origin
    }

    /// Whether this affine expression is used only as a numerator factor.
    ///
    /// Examples
    /// --------
    /// >>> numerator_only = surface.numerator_only
    #[getter]
    pub(crate) fn numerator_only(&self) -> bool {
        self.numerator_only
    }

    /// Index within its surface category and representation.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.index
    #[getter]
    pub(crate) fn index(&self) -> usize {
        self.index
    }
    /// ``"E"`` for consistent nonzero energy signs, ``"H"`` for mixed signs.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.kind
    #[getter]
    pub(crate) fn kind(&self) -> &'static str {
        self.kind
    }
    /// Result-local Symbolica placeholder for this stored surface.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.symbol
    #[getter]
    pub(crate) fn symbol(&self) -> PythonExpression {
        PythonExpression {
            expr: self.symbol.clone(),
        }
    }
    /// Exact signed OSE coefficients, keyed by physical edge ID.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.energy_coefficients
    #[getter]
    pub(crate) fn energy_coefficients(&self) -> BTreeMap<usize, PythonExpression> {
        self.energy_coefficients
            .iter()
            .map(|(e, c)| (*e, PythonExpression { expr: c.clone() }))
            .collect()
    }
    /// Exact coefficients of external ``gammalooprs::Q(edge, spenso::cind(0))``.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.external_shift
    #[getter]
    pub(crate) fn external_shift(&self) -> BTreeMap<usize, PythonExpression> {
        self.external_shift
            .iter()
            .map(|(e, c)| (*e, PythonExpression { expr: c.clone() }))
            .collect()
    }
    /// Energy-independent additive term.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.constant
    #[getter]
    pub(crate) fn constant(&self) -> PythonExpression {
        PythonExpression {
            expr: self.constant.clone(),
        }
    }
    /// CFF region vertices; LTD regions belong to individual surface occurrences.
    ///
    /// Examples
    /// --------
    /// >>> value = surface.vertices
    #[getter]
    pub(crate) fn vertices(&self) -> Vec<usize> {
        self.vertices.clone()
    }
    /// Return the affine definition :math:`c_0+\sum_e c_e E_e+\sum_a d_a Q_a^0`.
    ///
    /// Examples
    /// --------
    /// >>> expression = surface.to_expression()
    ///
    /// Returns
    /// -------
    /// TensorExpression
    ///     Signed energy combination in physical ``OSE`` and external ``Q`` symbols.
    fn to_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.atom(true), None)
    }
    /// Show the stored surface and its affine definition.
    ///
    /// Examples
    /// --------
    /// >>> print(surface)
    fn __repr__(&self) -> String {
        format!("{} = {}", self.symbol, self.atom(true))
    }

    /// Render the stored surface's mathematical definition and E/H category.
    ///
    /// Examples
    /// --------
    /// >>> html = surface._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(record_html(
            "EnergySurface",
            &format!("{}-surface {}", self.kind, self.index),
            &[("Definition", self.definition_html(py)?)],
        ))
    }
}

/// A signed denominator factor referring to a stored energy surface.
///
/// For stored surface :math:`s`, occurrence sign :math:`\epsilon\in\{-1,+1\}`
/// and positive multiplicity :math:`p`, ``to_expression()`` returns
/// :math:`(\epsilon s)^p`. This is the denominator, before inversion.
/// ``sign`` relates the occurrence to the stored surface; it is not an LTD
/// cut-edge pole sign.
///
/// Examples
/// --------
/// >>> from symbolica.community import hepkit as hep
/// >>> process = hep.Model.phi3().process(["phi"], ["phi", "phi"])
/// >>> diagrams = process.generate_diagrams(
/// ...     loops=1, max_vertices=3, maximum_bridges=0, progress=None
/// ... )
/// >>> diagram = diagrams.diagrams[0]
/// >>> rep = diagram.integrate_energy(method="cff")
/// >>> factor = rep.orientations[0].families[0].factors[0]
/// >>> denominator = factor.to_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "SurfaceFactor",
    module = "symbolica.community.hepkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PySurfaceFactor {
    pub(crate) surface: PyEnergySurface,
    pub(crate) sign: i8,
    pub(crate) power: usize,
}
impl PySurfaceFactor {
    pub(crate) fn atom(&self, expand_surfaces: bool) -> Atom {
        (Atom::num(self.sign as i64) * self.surface.atom(expand_surfaces)).pow(self.power as i64)
    }

    pub(crate) fn definition_html(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "{} = {}",
            expression_html(
                py,
                PythonExpression {
                    expr: self.atom(false)
                }
            )?,
            expression_html(
                py,
                PythonExpression {
                    expr: self.atom(true)
                }
            )?
        ))
    }
}
#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[spenso_macros::track_usage(crate::record_usage)]
#[pymethods]
impl PySurfaceFactor {
    /// Show the signed denominator and its multiplicity.
    ///
    /// Examples
    /// --------
    /// >>> text = repr(factor)
    fn __repr__(&self) -> String {
        format!(
            "SurfaceFactor(surface={}, sign={}, power={})",
            self.surface.symbol, self.sign, self.power
        )
    }

    /// Render this denominator before inversion, preserving its sign and power.
    ///
    /// Examples
    /// --------
    /// >>> html = factor._repr_html_()
    fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        Ok(record_html(
            "SurfaceFactor",
            "Denominator",
            &[("Factor", self.definition_html(py)?)],
        ))
    }
    /// Stored surface referenced by this factor.
    ///
    /// Examples
    /// --------
    /// >>> value = factor.surface
    #[getter]
    fn surface(&self) -> PyEnergySurface {
        self.surface.clone()
    }
    /// Overall sign relative to the stored surface, +1 or -1.
    ///
    /// Examples
    /// --------
    /// >>> value = factor.sign
    #[getter]
    fn sign(&self) -> i8 {
        self.sign
    }
    /// Positive multiplicity in this denominator product.
    ///
    /// Examples
    /// --------
    /// >>> value = factor.power
    #[getter]
    fn power(&self) -> usize {
        self.power
    }
    /// Return the signed denominator factor, including its power.
    ///
    /// Examples
    /// --------
    /// Using the setup in the class example:
    ///
    /// >>> expression = factor.to_expression()
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
    ///     The denominator :math:`(\epsilon s)^p`, before inversion.
    #[pyo3(signature=(*, expand_surfaces=true))]
    pub(crate) fn to_expression(
        &self,
        py: Python<'_>,
        expand_surfaces: bool,
    ) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.atom(expand_surfaces), None)
    }
}
