use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use feynkit_cff::{
    CffExpression, CffGenerator, CffOptions, CffReport, CffResult, CutPropagator, EdgeOrientation,
    FeynmanDiagramCffExt, OrientationExpression, RaisedEnergySurfaceGroup, Surface, SurfaceCache,
    SurfaceId, SurfacePole,
};
use linnet::half_edge::{
    involution::HedgePair,
    subgraph::{ModifySubSet, SuBitGraph},
};
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::{Atom, AtomCore, AtomView, Symbol},
};

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
    owner: Arc<()>,
}

impl PyCffSurface {
    fn new(id: SurfaceId, cache: &SurfaceCache, owner: &Arc<()>) -> Self {
        Self {
            id,
            surface: cache.get(id),
            owner: Arc::clone(owner),
        }
    }

    fn symbol_name_for(id: SurfaceId) -> Option<String> {
        match id {
            SurfaceId::Energy(id) => Some(Atom::from(id).to_canonical_string()),
            SurfaceId::H(id) => Some(Atom::from(id).to_canonical_string()),
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

    /// Return the original diagram vertex IDs enclosed by this surface.
    ///
    /// Contracted CFF vertices retain the identities of all interaction
    /// vertices they contain, including for selected subgraphs.
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
    /// >>> print(surface)  # Symbolica denominator variable, for example feynkit_cff::η(0)
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
    owner: Arc<()>,
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
                    .map(|surface| PyCffSurface::new(surface, &self.surfaces, &self.owner))
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
    owner: Arc<()>,
    normalization: Atom,
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
                owner: Arc::clone(&self.owner),
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
                &self.owner,
            )
        });
        let h = (0..self.inner.surfaces.h_surfaces().len()).map(|index| {
            PyCffSurface::new(
                SurfaceId::H(feynkit_cff::HSurfaceId(index)),
                &self.inner.surfaces,
                &self.owner,
            )
        });
        energy.chain(h).collect()
    }

    /// Convert to the canonical eta/H denominator expression.
    ///
    /// ``expand_surfaces`` substitutes on-shell/external energies. ``normalized``
    /// additionally includes the -1/(2 E) factors and GammaLoop loop measure;
    /// it implies ``expand_surfaces``. Numerators and global weights stay separate.
    ///
    /// Examples
    /// --------
    /// >>> result.to_expression(normalized=True)
    ///
    /// Parameters
    /// ----------
    /// expand_surfaces : bool
    ///     Substitute canonical on-shell and external energies.
    /// normalized : bool
    ///     Include the energy products and spatial loop measure.
    #[pyo3(signature = (*, expand_surfaces=false, normalized=false))]
    fn to_expression(&self, expand_surfaces: bool, normalized: bool) -> PythonExpression {
        let mut expr = self.inner.expression.to_atom();
        if expand_surfaces || normalized {
            expr = expr.replace_multiple(self.inner.surfaces.replacements_with(
                feynkit_cff::symbols::on_shell_atom,
                feynkit_cff::symbols::external_energy_atom,
            ));
        }
        if normalized {
            expr *= &self.normalization;
        }
        PythonExpression { expr }
    }

    /// Expand one surface belonging to this result into canonical energy symbols.
    ///
    /// Examples
    /// --------
    /// >>> result.surface_expression(result.surfaces[0])
    ///
    /// Parameters
    /// ----------
    /// surface : CffSurface
    ///     A surface obtained from this result.
    fn surface_expression(&self, surface: &PyCffSurface) -> PyResult<PythonExpression> {
        if !Arc::ptr_eq(&self.owner, &surface.owner) {
            return Err(error::CffError::new_err(
                "surface belongs to a different CFF result",
            ));
        }
        let surface = surface
            .surface
            .as_ref()
            .ok_or_else(|| error::CffError::new_err("surface is absent from the result arena"))?;
        Ok(PythonExpression {
            expr: surface.to_atom_with(
                feynkit_cff::symbols::on_shell_atom,
                feynkit_cff::symbols::external_energy_atom,
            ),
        })
    }

    /// Group equivalent energy surfaces after identifying raised propagator edges.
    /// ``edge_representatives`` maps repeated edges to their canonical edge.
    ///
    /// Examples
    /// --------
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
        let groups = self
            .inner
            .expression
            .raised_energy_surfaces(&self.inner.surfaces, |edge| {
                feynkit_cff::EdgeId::new(
                    representatives
                        .get(&edge.index())
                        .copied()
                        .unwrap_or(edge.index()),
                )
            })
            .map_err(error::cff)?;
        Ok(groups
            .groups
            .into_iter()
            .map(|inner| PyCffSurfaceGroup {
                inner,
                surfaces: self.inner.surfaces.clone(),
                owner: Arc::clone(&self.owner),
            })
            .collect())
    }

    /// Return coefficients of each inverse surface power, indexed from order one.
    /// These are pole coefficients, before analytic residue derivatives.
    ///
    /// Examples
    /// --------
    /// >>> coefficients = result.pole_coefficients(result.raised_surface_groups()[0])
    /// >>> [coefficient.to_expression() for coefficient in coefficients]
    ///
    /// Parameters
    /// ----------
    /// group : CffSurfaceGroup
    ///     A raised-surface group belonging to this result.
    fn pole_coefficients(&self, group: &PyCffSurfaceGroup) -> PyResult<Vec<Self>> {
        self.validate_group(group)?;
        Ok(self
            .inner
            .select_energy_surface_residue(&group.inner)
            .into_iter()
            .map(|inner| Self {
                inner,
                owner: Arc::clone(&self.owner),
                normalization: self.normalization.clone(),
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
    /// >>> result.residue(group, variable=t, root=t_star, surface=eta, coefficient=numerator)
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
    /// normalized : bool
    ///     Include the generated CFF normalization in the coefficient.
    /// replacements : list[tuple[Expression, Expression]], optional
    ///     Route all energy dependence to the integration variable before differentiating.
    #[pyo3(signature = (group, *, variable, root, surface, coefficient, normalized=false, replacements=None))]
    #[allow(clippy::too_many_arguments)] // Keep Python's residue coordinates explicit keyword arguments.
    fn residue(
        &self,
        group: &PyCffSurfaceGroup,
        variable: ConvertibleToExpression,
        root: ConvertibleToExpression,
        surface: ConvertibleToExpression,
        coefficient: ConvertibleToExpression,
        normalized: bool,
        replacements: Option<Vec<(ConvertibleToExpression, ConvertibleToExpression)>>,
    ) -> PyResult<PythonExpression> {
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
            let remaining = (coefficient.clone() * term.to_expression(true, normalized).expr)
                .replace_multiple(&replacements);
            result += SurfacePole {
                surface: surface.clone(),
                order: index + 1,
            }
            .residue(&remaining, variable, &root)
            .map_err(error::cff)?;
        }
        Ok(PythonExpression { expr: result })
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
        let expression = self.to_expression(false, false)._repr_html_()?;
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
        self.to_expression(false, false)
            ._repr_pretty_(pretty, false)?;
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
    /// subgraph : linnet.Subgraph, optional
    ///     Graph-bound selection from diagram.to_linnet().
    #[pyo3(signature = (diagram, *, subgraph=None))]
    fn generate(
        &self,
        py: Python<'_>,
        diagram: &PyFeynmanDiagram,
        #[gen_stub(override_type(type_repr="linnet.Subgraph | None", imports=("linnet")))]
        subgraph: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyCffResult> {
        let selection = diagram.selection(py, subgraph)?;
        PyCffResult::build(py, diagram, self.inner.options().clone(), selection)
    }
}

pub(crate) fn build_cff_for_diagram(
    py: Python<'_>,
    diagram: &PyFeynmanDiagram,
    max_orientations: Option<usize>,
    fixed_orientations: Option<BTreeMap<usize, bool>>,
    contracted_edges: Option<Vec<usize>>,
    initial_state_edges: Option<Vec<usize>>,
    subgraph: Option<SuBitGraph>,
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

    let selection = subgraph.unwrap_or_else(|| diagram.inner.underlying().full_filter());
    PyCffResult::build(py, diagram, options, selection)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCffSurface>()?;
    module.add_class::<PyCffOrientation>()?;
    module.add_class::<PyCffReport>()?;
    module.add_class::<PyCffResult>()?;
    module.add_class::<PyCffGenerator>()?;
    module.add_class::<PyCffSurfaceGroup>()?;
    module.add_class::<PyCutPropagator>()?;
    Ok(())
}

impl PyCffResult {
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
            let mut selected_internal: SuBitGraph = graph.empty_subgraph();
            let mut contracted: SuBitGraph = graph.empty_subgraph();
            let mut energies = Vec::new();
            for (pair, edge, data) in graph.iter_edges_of(&selection) {
                if let HedgePair::Paired { source, sink } = pair {
                    if initial_edges.contains(&edge) || data.data.is_dummy {
                        continue;
                    }
                    selected_internal.add(source);
                    selected_internal.add(sink);
                    if options.contracted_edges().contains(&edge) {
                        contracted.add(source);
                        contracted.add(sink);
                    } else {
                        energies.push(feynkit_cff::symbols::on_shell_atom(edge));
                    }
                }
            }
            let loops = graph
                .cyclotomatic_number(&selected_internal)
                .saturating_sub(graph.cyclotomatic_number(&contracted));
            let normalization = CffExpression::measure_normalization(loops)
                * CffExpression::inverse_energy_product(energies);
            let inner = diagram
                .build_cff_subgraph(&selection, options)
                .map_err(error::cff)?;
            Ok(Self {
                inner,
                owner: Arc::new(()),
                normalization,
            })
        })
    }
}

/// Equivalent energy surfaces and their maximum simultaneous pole order.
///
/// Examples
/// --------
/// >>> group = result.raised_surface_groups()[0]
/// >>> print(group.max_order, group.surfaces)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CffSurfaceGroup",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCffSurfaceGroup {
    inner: RaisedEnergySurfaceGroup,
    surfaces: SurfaceCache,
    owner: Arc<()>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCffSurfaceGroup {
    #[getter]
    /// Highest inverse surface power occurring on one CFF branch.
    fn max_order(&self) -> usize {
        self.inner.max_occurrence
    }
    #[getter]
    /// Canonical energy surfaces identified by this group.
    fn surfaces(&self) -> Vec<PyCffSurface> {
        self.inner
            .surface_ids
            .iter()
            .map(|id| PyCffSurface::new(SurfaceId::Energy(*id), &self.surfaces, &self.owner))
            .collect()
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
/// >>> cut = CutPropagator(q0, E, power=2)
/// >>> cut.apply(q0**2, q0)
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CutPropagator",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCutPropagator {
    pub(crate) inner: CutPropagator,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCutPropagator {
    /// Construct a reciprocal-propagator cut with explicit orientation and i0 sign.
    ///
    /// Examples
    /// --------
    /// >>> cut = CutPropagator(q0, E, power=3, orientation=-1)
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
