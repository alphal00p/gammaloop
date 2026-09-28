use feynkit_graph::{IntegralFamily, IntegralMapping, PropagatorMapping};
use pyo3::prelude::*;
use symbolica::api::python::PythonExpression;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{error, graph::PyFeynmanDiagram, kinematics::PyKinematics};

/// An ordered set of inverse propagators sharing loop momenta and external kinematics.
///
/// Pass denominators such as ``k.k - m2``, not their reciprocals. An integral
/// with powers ``[a1, a2, ...]`` has integrand ``1/(D1**a1 * D2**a2 * ...)``:
/// zero powers omit a denominator and negative powers put it in the numerator.
/// All power vectors and denominator labels follow the original input order.
///
/// Declare independent external momenta and their scalar products through
/// ``Kinematics``. ``is_independent`` checks the denominator basis for linear
/// dependence; ``is_complete`` checks whether it spans every loop scalar product.
/// For L loops and E independent external momenta this space has
/// L*(L+1)/2 + L*E coordinates. ``complete()`` returns a new family with missing
/// coordinates appended; give these auxiliary slots nonpositive integral powers.
/// Use ``partial_fraction()`` before completion when denominators are dependent.
///
/// Use this class to rewrite scalar numerators, compare momentum routings, test
/// scalelessness, and construct Symanzik polynomials. Reduction is provided by
/// ``hep.IBPFamily`` or, for one loop, ``hep.oneloop.reduce``. No integration
/// prescription or loop-measure normalization is inferred by this container.
///
/// Examples
/// --------
/// Build a massless bubble with external virtuality ``p.p = s``. The method
/// examples below reuse this family and its symbols.
///
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> D, k, p, s = S("D", "k", "p", "s")
/// >>> d1, d2, x1, x2 = S("d1", "d2", "x1", "x2")
/// >>> kin = hep.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)
/// >>> denominators = [kin.scalar_product(k, k), kin.scalar_product(k-p, k-p)]
/// >>> family = hep.IntegralFamily([k], [p], denominators, kinematics=kin)
/// >>> assert family.rank == 2 and family.is_complete and family.is_independent
/// >>> rewritten = family.rewrite_numerator(kin.scalar_product(k, p), [d1, d2])
/// >>> assert (rewritten - (d1 - d2 + s)/2).expand() == E("0")
/// >>> U, F = family.symanzik([x1, x2])
/// >>> assert U == x1 + x2
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "IntegralFamily",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyIntegralFamily {
    pub(crate) inner: IntegralFamily,
}

impl PyIntegralFamily {
    /// Borrow the shared native family when calling a reduction backend.
    pub fn as_family(&self) -> &IntegralFamily {
        &self.inner
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyIntegralFamily {
    /// Build a complete family from the diagram's routed internal propagators.
    ///
    /// Graph propagators retain their edge order and model masses. Preferred
    /// independent dot products are appended in order when they increase the
    /// rank; automatic scalar products fill any remaining directions. Auxiliary
    /// denominators have nonpositive powers when representing the original
    /// integral. Tree diagrams and dependent graph propagators raise DiagramError.
    /// Use ``diagram.propagator_family()`` to extract dependent propagators for
    /// partial fractioning before completing the resulting families.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hep
    /// >>> model = hep.Model.phi4()
    /// >>> process = model.process(["phi", "phi"], ["phi", "phi"])
    /// >>> result = process.generate_diagrams(loops=1)
    /// >>> diagram = result.diagrams[0]
    /// >>> family = hep.IntegralFamily.from_diagram(diagram)
    /// >>> assert family.is_complete
    ///
    /// Parameters
    /// ----------
    /// diagram : FeynmanDiagram
    ///     A complete diagram with its chosen loop-momentum routing.
    /// independent_dot_products : list[Expression] or None, optional
    ///     Preferred auxiliary inverse propagators in the routed momenta.
    ///     None selects a suitable basis automatically.
    /// kinematics : Kinematics or None, optional
    ///     External assumptions and dimension; defaults to the diagram's
    ///     symbolic dimension with no on-shell assumptions.
    #[staticmethod]
    #[pyo3(signature = (diagram, independent_dot_products=None, *, kinematics=None))]
    pub(crate) fn from_diagram(
        diagram: &PyFeynmanDiagram,
        independent_dot_products: Option<Vec<PythonExpression>>,
        kinematics: Option<&PyKinematics>,
    ) -> PyResult<Self> {
        diagram.require_complete()?;
        let default = feynkit_kinematics::Kinematics::in_dimension(&symbolica::atom::Atom::var(
            feynkit_graph::symbols::dimension(),
        ))
        .expect("the shared dimension is a Symbolica symbol");
        let kinematics = kinematics.map(|kin| &kin.inner).unwrap_or(&default);
        let candidates = independent_dot_products
            .unwrap_or_default()
            .into_iter()
            .map(|product| product.expr)
            .collect::<Vec<_>>();
        Ok(Self {
            inner: IntegralFamily::from_diagram(&diagram.inner, kinematics, &candidates)
                .map_err(error::diagram)?,
        })
    }

    /// Compute the independent loop scalar products and denominator rank.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S, E
    /// >>> from symbolica.community import hep
    /// >>> D, k, p, s = S("D", "k", "p", "s")
    /// >>> d1, d2, x1, x2 = S("d1", "d2", "x1", "x2")
    /// >>> kin = hep.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)
    /// >>> denominators = [kin.scalar_product(k, k), kin.scalar_product(k-p, k-p)]
    /// >>> family = hep.IntegralFamily([k], [p], denominators, kinematics=kin)
    /// >>> assert family.rank == 2
    ///
    /// Parameters
    /// ----------
    /// loop_momenta : list[Expression]
    ///     Distinct unindexed integrated momentum names.
    /// external_momenta : list[Expression]
    ///     Independent unindexed external momentum names.
    /// denominators : list[Expression]
    ///     Ordered inverse propagators, not their reciprocals.
    /// kinematics : Kinematics | None
    ///     External assumptions and dimension; defaults to unconstrained 4D.
    #[new]
    #[pyo3(signature = (loop_momenta, external_momenta, denominators, *, kinematics=None))]
    fn new(
        loop_momenta: Vec<PythonExpression>,
        external_momenta: Vec<PythonExpression>,
        denominators: Vec<PythonExpression>,
        kinematics: Option<&PyKinematics>,
    ) -> PyResult<Self> {
        let default = feynkit_kinematics::Kinematics::new();
        let kinematics = kinematics.map(|kin| &kin.inner).unwrap_or(&default);
        Ok(Self {
            inner: IntegralFamily::new(
                loop_momenta.into_iter().map(|p| p.expr).collect(),
                external_momenta.into_iter().map(|p| p.expr).collect(),
                denominators.into_iter().map(|d| d.expr).collect(),
                kinematics,
            )
            .map_err(error::integral_family)?,
        })
    }

    /// Integrated momentum names in their original order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.loop_momenta == [k]
    #[getter]
    fn loop_momenta(&self) -> Vec<PythonExpression> {
        self.inner
            .loop_momenta()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Scoped kinematics, including the family momentum declarations.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.kinematics.scalar_product(p, p) == s
    #[getter]
    fn kinematics(&self) -> PyKinematics {
        PyKinematics {
            inner: self.inner.kinematics().clone(),
        }
    }

    /// Independent external momentum names in their original order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.external_momenta == [p]
    #[getter]
    fn external_momenta(&self) -> Vec<PythonExpression> {
        self.inner
            .external_momenta()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Loop-loop and loop-external products spanning the numerator space.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert len(family.scalar_products) == 2
    #[getter]
    fn scalar_products(&self) -> Vec<PythonExpression> {
        self.inner
            .scalar_products()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Ordered inverse propagators, including any auxiliary completion terms.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.denominators == denominators
    /// >>> weighted = list(zip(family.denominators, [1, 2]))
    #[getter]
    fn denominators(&self) -> Vec<PythonExpression> {
        self.inner
            .denominators()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Number of independent affine forms in the loop scalar products.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.rank == 2
    #[getter]
    fn rank(&self) -> usize {
        self.inner.rank()
    }

    /// Whether the inverse propagators span every loop scalar product.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.is_complete
    /// >>> assert family.rank == len(family.scalar_products)
    #[getter]
    fn is_complete(&self) -> bool {
        self.inner.is_complete()
    }

    /// Whether no denominator can be eliminated by an affine relation.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> assert family.is_independent
    /// >>> assert family.rank == len(family.denominators)
    #[getter]
    fn is_independent(&self) -> bool {
        self.inner.is_independent()
    }

    /// Return a compact summary of the family and its scalar-product rank.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> print(family)
    fn __repr__(&self) -> String {
        format!(
            "IntegralFamily(loops={}, external_momenta={}, denominators={}, rank={}/{}, complete={}, independent={})",
            self.inner.loop_momenta().len(),
            self.inner.external_momenta().len(),
            self.inner.denominators().len(),
            self.inner.rank(),
            self.inner.scalar_products().len(),
            self.inner.is_complete(),
            self.inner.is_independent(),
        )
    }

    /// Display ordered inverse propagators and family metadata in a notebook.
    ///
    /// Expressions use Symbolica's native HTML printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> from IPython.display import display
    /// >>> display(family)
    fn _repr_html_(&self) -> PyResult<String> {
        let dimension = self.inner.kinematics().dimension().to_symbolic();
        let mut metadata = String::new();
        for (label, expressions) in [
            ("Loop momenta", self.inner.loop_momenta()),
            ("External momenta", self.inner.external_momenta()),
            ("Dimension", std::slice::from_ref(&dimension)),
        ] {
            let values = expressions
                .iter()
                .map(|expr| PythonExpression::from(expr.clone())._repr_html_())
                .collect::<PyResult<Vec<_>>>()?
                .join(", ");
            metadata.push_str(&format!(
                "<div class=\"feynkit-family-metadata\"><strong>{label}:</strong> {}</div>",
                if values.is_empty() {
                    "<em>none</em>"
                } else {
                    &values
                },
            ));
        }
        let mut rows = String::new();
        for (index, denominator) in self.inner.denominators().iter().enumerate() {
            let expression = PythonExpression::from(denominator.clone())._repr_html_()?;
            rows.push_str(&format!(
                "<tr><th scope=\"row\" style=\"padding:.3rem .65rem;text-align:right\">\
                 D<sub>{}</sub></th><td style=\"padding:.3rem .65rem;text-align:left\">\
                 {expression}</td></tr>",
                index + 1,
            ));
        }
        if rows.is_empty() {
            rows.push_str("<tr><td colspan=\"2\"><em>No inverse propagators</em></td></tr>");
        }
        Ok(format!(
            "<section class=\"feynkit-integral-family\" style=\"max-width:100%;overflow-x:auto\">\
             <style>.feynkit-integral-family > .feynkit-family-metadata > div {{ display:inline; }}</style>\
             <strong>Integral family</strong>{metadata}\
             <div>Rank: {} / {} &middot; {} &middot; {}</div>\
             <table style=\"border-collapse:collapse;margin-top:.4rem\">\
             <caption style=\"text-align:left\">Ordered inverse propagators ({})</caption>\
             <tbody>{rows}</tbody></table></section>",
            self.inner.rank(),
            self.inner.scalar_products().len(),
            if self.inner.is_complete() {
                "complete"
            } else {
                "incomplete"
            },
            if self.inner.is_independent() {
                "independent"
            } else {
                "dependent"
            },
            self.inner.denominators().len(),
        ))
    }

    /// Write the summary and ordered denominators using Symbolica's text printer.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> from IPython.lib.pretty import pretty
    /// >>> text = pretty(family)
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     IPython's pretty printer, providing a ``text`` method.
    /// cycle : bool
    ///     Whether this family is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        let mut text = if cycle {
            "IntegralFamily(...)".to_owned()
        } else {
            self.__repr__()
        };
        if !cycle {
            for (index, denominator) in self.inner.denominators().iter().enumerate() {
                let expression = PythonExpression::from(denominator.clone()).__str__()?;
                text.push_str(&format!("\n  D{} = {expression}", index + 1));
            }
        }
        pretty.call_method1("text", (text,))?;
        Ok(())
    }

    /// Complete a basis, preferring supplied inverse propagators.
    ///
    /// Original propagators retain their positions. Dependent families must
    /// first be partial-fractioned. Added propagators carry nonpositive powers
    /// when used to represent numerator factors in an IBP integral list.
    /// Candidates are tried in order after applying this family's kinematics.
    /// Redundant candidates are skipped; bare scalar products fill any missing
    /// directions. All candidates must be affine in the loop scalar products.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> incomplete = hep.IntegralFamily([k], [p], denominators[:1], kinematics=kin)
    /// >>> completed = incomplete.complete()
    /// >>> assert not incomplete.is_complete and completed.is_complete
    /// >>> assert completed.denominators[0] == denominators[0]
    /// >>> powers = [1, 0]  # the appended slot is absent from the original integral
    ///
    /// Parameters
    /// ----------
    /// candidates : list[Expression] | None
    ///     Preferred auxiliary inverse propagators in this family's momentum
    ///     coordinates. Defaults to using only bare scalar products.
    #[pyo3(signature = (*, candidates=None))]
    fn complete(&self, candidates: Option<Vec<PythonExpression>>) -> PyResult<Self> {
        let candidates = candidates
            .unwrap_or_default()
            .into_iter()
            .map(|x| x.expr)
            .collect::<Vec<_>>();
        Ok(Self {
            inner: self
                .inner
                .complete(&candidates)
                .map_err(error::integral_family)?,
        })
    }

    /// Partial-fraction dependent propagators while preserving family order.
    ///
    /// Both affine mass shifts and homogeneous dependencies are supported.
    /// The identity is algebraic: no momentum shifts or scaleless-term removal
    /// are applied. Impose exceptional kinematics before building the family;
    /// generic external invariants may appear in the returned coefficients.
    ///
    /// Examples
    /// --------
    /// Using the symbols and kinematics in ``IntegralFamily``:
    ///
    /// >>> m2 = S("m2")
    /// >>> dependent = hep.IntegralFamily([k], [], [kin.scalar_product(k, k),
    /// ...     kin.scalar_product(k, k) - m2], kinematics=kin)
    /// >>> terms = dependent.partial_fraction([1, 1])
    /// >>> assert all(sum(power > 0 for power in powers) == 1 for coefficient, powers in terms)
    ///
    /// Parameters
    /// ----------
    /// powers : list[int]
    ///     Signed powers in family order; negative powers are numerator factors.
    /// max_states : int
    ///     Maximum number of intermediate exponent vectors; defaults to 100000.
    ///
    /// Returns
    /// -------
    /// list[tuple[Expression, list[int]]]
    ///     Coefficients and powers whose positive-power denominators are independent.
    #[pyo3(signature = (powers, *, max_states=100_000))]
    fn partial_fraction(
        &self,
        powers: Vec<i32>,
        max_states: usize,
    ) -> PyResult<Vec<(PythonExpression, Vec<i32>)>> {
        Ok(self
            .inner
            .partial_fraction(&powers, max_states)
            .map_err(error::integral_family)?
            .into_iter()
            .map(|(coefficient, powers)| (coefficient.into(), powers))
            .collect())
    }

    /// Verify explicit source-loop images in another family's coordinates.
    ///
    /// Real affine images must have loop determinant +1 or -1 and map every
    /// source propagator to a distinct equal target propagator. External names,
    /// dimensions and scalar-product assumptions must agree. Additional target
    /// propagators are allowed for subtopology embeddings.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> mapping = family.mapping_to(family, [k])
    /// >>> assert mapping.map_powers([1, 2]) == [1, 2]
    ///
    /// Parameters
    /// ----------
    /// target : IntegralFamily
    ///     Family in whose coordinates the images are expressed.
    /// loop_images : list[Expression]
    ///     One image per source loop momentum, in source order.
    ///
    /// Returns
    /// -------
    /// IntegralMapping | None
    ///     Verified mapping, or None when the images fail the equivalence checks.
    fn mapping_to(
        &self,
        target: &PyIntegralFamily,
        loop_images: Vec<PythonExpression>,
    ) -> PyResult<Option<PyIntegralMapping>> {
        let images = loop_images.into_iter().map(|p| p.expr).collect::<Vec<_>>();
        self.inner
            .mapping_to(&target.inner, &images)
            .map(|mapping| mapping.map(|inner| PyIntegralMapping { inner }))
            .map_err(error::integral_family)
    }

    /// Find a verified affine loop-momentum shift into the target family.
    ///
    /// Search derives candidates from independent quadratic propagators, then
    /// checks every propagator, including eikonal and auxiliary entries. It
    /// includes loop mixtures, reversals and external shifts. It does not test
    /// parametric identities that have no affine loop-momentum map.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> target = hep.IntegralFamily([k], [p], list(reversed(denominators)), kinematics=kin)
    /// >>> mapping = family.find_mapping(target)
    /// >>> assert mapping is not None
    /// >>> transformed = mapping.apply(kin.scalar_product(k, p))
    ///
    /// Parameters
    /// ----------
    /// target : IntegralFamily
    ///     Target family with the same external kinematics and loop count.
    /// max_candidates : int
    ///     Candidate budget; exhaustion raises IntegralFamilyError.
    ///
    /// Returns
    /// -------
    /// IntegralMapping | None
    ///     Verified mapping, or None if no supported candidate matches.
    #[pyo3(signature = (target, *, max_candidates=100_000))]
    fn find_mapping(
        &self,
        target: &PyIntegralFamily,
        max_candidates: usize,
    ) -> PyResult<Option<PyIntegralMapping>> {
        self.inner
            .find_mapping(&target.inner, max_candidates)
            .map(|mapping| mapping.map(|inner| PyIntegralMapping { inner }))
            .map_err(error::integral_family)
    }

    /// Group families and return verified maps to their representatives.
    ///
    /// Returns one ``(target_index, mapping)`` per input family in input order.
    /// Indices refer to the original list. Representatives map to themselves;
    /// every other mapping goes directly to a retained representative.
    /// The first compatible representative wins. Symbolica canonizes each
    /// Symanzik pair once, then native affine-shift search verifies all merges.
    /// A parametric equivalence without a verified loop map remains separate.
    ///
    /// Families need compatible external kinematics and nonsingular quadratic
    /// forms. Unsupported shift searches and exhausted budgets raise errors.
    /// Different propagator counts stay separate. Select sectors and remove
    /// certified scaleless integrals before grouping; complete bases afterwards.
    /// This does not exchange external momenta or infer integration prescriptions.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> families = [family, family]
    /// >>> mappings = hep.IntegralFamily.find_mappings(families)
    /// >>> assert [target for target, mapping in mappings] == [0, 0]
    ///
    /// Parameters
    /// ----------
    /// families : list[IntegralFamily]
    ///     Ordered families; earlier entries are preferred as representatives.
    /// max_candidates : int
    ///     Affine-shift candidate budget per pair; exhaustion raises an error.
    ///
    /// Returns
    /// -------
    /// list[tuple[int, IntegralMapping]]
    ///     Target index and verified map for each input, including representatives.
    #[staticmethod]
    #[pyo3(signature = (families, *, max_candidates=100_000))]
    fn find_mappings(
        families: Vec<PyIntegralFamily>,
        max_candidates: usize,
    ) -> PyResult<Vec<(usize, PyIntegralMapping)>> {
        let families = families.into_iter().map(|f| f.inner).collect::<Vec<_>>();
        IntegralFamily::find_mappings(&families, max_candidates)
            .map(|mappings| {
                mappings
                    .into_iter()
                    .map(|(target, inner)| (target, PyIntegralMapping { inner }))
                    .collect()
            })
            .map_err(error::integral_family)
    }

    /// Select the positive-power propagators of an integral's sector.
    ///
    /// Zero and negative powers are omitted. Loop variables and external
    /// kinematics are retained. This selects sector support; it does not remove
    /// numerator factors algebraically from the original integrand.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> sector = family.sector([1, 0])
    /// >>> assert sector.denominators == denominators[:1]
    ///
    /// Parameters
    /// ----------
    /// powers : list[int]
    ///     One signed propagator power per family denominator.
    fn sector(&self, powers: Vec<i32>) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.sector(&powers).map_err(error::integral_family)?,
        })
    }

    /// Find an unconstrained transverse loop direction proving scalelessness.
    ///
    /// Returns a nonzero real list ``w`` in loop order such that each inverse
    /// propagator is invariant under ``k_i -> k_i + w_i*r_perp`` for any
    /// vector orthogonal to the external span. The corresponding unrestricted
    /// transverse integral vanishes in dimensional regularization, including
    /// polynomial numerators. Apply ``sector(powers)`` first to exclude
    /// numerator-only entries.
    ///
    /// Requires a nonsingular external Gram matrix. Symbolic dimension is
    /// interpreted generically; a concrete dimension must exceed the external
    /// basis size. None means no certificate was found. A vanishing Symanzik U
    /// alone is not sufficient, and complex loop directions are not accepted.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> empty_sector = family.sector([0, 0])
    /// >>> direction = empty_sector.scaleless_transverse_direction()
    /// >>> assert direction is not None  # no denominators constrain the loop momentum
    ///
    /// Returns
    /// -------
    /// list[Expression] | None
    ///     Verified real loop direction, or no certificate.
    fn scaleless_transverse_direction(&self) -> PyResult<Option<Vec<PythonExpression>>> {
        self.inner
            .scaleless_transverse_direction()
            .map(|direction| direction.map(|v| v.into_iter().map(Into::into).collect()))
            .map_err(error::integral_family)
    }

    /// Find parameter weights proving a sector scaleless in dimensional regularization.
    ///
    /// For G=U+F, the returned weights satisfy sum(w_i*x_i*dG/dx_i)=G.
    /// None means this criterion did not detect scalelessness, not that the
    /// integral is nonzero. Every denominator is treated as present; use
    /// ``sector(powers)`` first to select positive-power entries. Singular
    /// quadratic loop forms raise IntegralFamilyError: their algebraic U/F
    /// polynomials do not establish this parametric scaling certificate.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> tadpole = family.sector([1, 0])
    /// >>> weights = tadpole.scaleless_scaling([x1])
    /// >>> assert weights is not None  # massless tadpole vanishes in dimensional regularization
    ///
    /// Parameters
    /// ----------
    /// parameters : list[Expression]
    ///     One distinct new symbol or labeled call per sector denominator.
    ///
    /// Returns
    /// -------
    /// list[Expression] | None
    ///     Exact scaling weights in parameter order, or no certificate.
    fn scaleless_scaling(
        &self,
        parameters: Vec<PythonExpression>,
    ) -> PyResult<Option<Vec<PythonExpression>>> {
        let parameters = parameters.into_iter().map(|x| x.expr).collect::<Vec<_>>();
        self.inner
            .scaleless_scaling(&parameters)
            .map(|weights| weights.map(|weights| weights.into_iter().map(Into::into).collect()))
            .map_err(error::integral_family)
    }

    /// Find a parameter permutation identifying both Symanzik polynomials.
    ///
    /// Symbolica canonizes polynomial incidence graphs, preserving coefficients,
    /// powers and the distinction between U and F. This can identify families
    /// without a loop shift at fixed external momenta. The result supplies no
    /// momentum or tensor-numerator substitution and does not check contours or
    /// propagator prescriptions. Both families must have equal denominator
    /// counts and the same external kinematics. Singular quadratic forms are
    /// rejected because their U/F polynomials can discard physical parameters.
    ///
    /// Examples
    /// --------
    /// Using the bubble setup in ``IntegralFamily``:
    ///
    /// >>> target = hep.IntegralFamily([k], [p], list(reversed(denominators)), kinematics=kin)
    /// >>> mapping = family.parametric_mapping(target, [x1, x2])
    /// >>> assert mapping is not None
    /// >>> mapped_powers = mapping.map_powers([1, 2])
    ///
    /// Parameters
    /// ----------
    /// target : IntegralFamily
    ///     Family whose U and F polynomials are compared.
    /// parameters : list[Expression]
    ///     One distinct new symbol or labeled call per source denominator.
    ///
    /// Returns
    /// -------
    /// PropagatorMapping | None
    ///     Formal parameter permutation, or None if the polynomials differ.
    fn parametric_mapping(
        &self,
        target: &PyIntegralFamily,
        parameters: Vec<PythonExpression>,
    ) -> PyResult<Option<PyPropagatorMapping>> {
        let parameters = parameters.into_iter().map(|x| x.expr).collect::<Vec<_>>();
        self.inner
            .parametric_mapping(&target.inner, &parameters)
            .map(|mapping| mapping.map(|inner| PyPropagatorMapping { inner }))
            .map_err(error::integral_family)
    }

    /// Compute the Symanzik polynomials U and F in propagator order.
    ///
    /// Uses Minkowski inverse propagators: k^2-m^2 gives U=x and F=m^2*x^2.
    /// For a weighted denominator k.M.k + 2 k.Q + J, this returns
    /// U=det(M) and F=Q.adj(M).Q-U*J, using Symbolica determinants and cofactors.
    /// Singular quadratic forms are accepted as algebraic polynomial data;
    /// they do not establish a Gaussian integration formula or scalelessness.
    /// This prepares polynomials without performing integration.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> U, F = family.symanzik([x1, x2])
    /// >>> assert U == x1 + x2  # two standard one-loop propagators
    ///
    /// Parameters
    /// ----------
    /// parameters : list[Expression]
    ///     One distinct new symbol or labeled call per inverse propagator.
    ///
    /// Returns
    /// -------
    /// tuple[Expression, Expression]
    ///     First and second Symanzik polynomials, respectively.
    fn symanzik(
        &self,
        parameters: Vec<PythonExpression>,
    ) -> PyResult<(PythonExpression, PythonExpression)> {
        let parameters = parameters.into_iter().map(|x| x.expr).collect::<Vec<_>>();
        let (u, f) = self
            .inner
            .symanzik(&parameters)
            .map_err(error::integral_family)?;
        Ok((PythonExpression { expr: u }, PythonExpression { expr: f }))
    }

    /// Solve loop scalar products in terms of inverse-propagator labels.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> rules = family.scalar_product_rules([d1, d2])
    /// >>> assert len(rules) == 2
    ///
    /// Parameters
    /// ----------
    /// labels : list[Expression]
    ///     One distinct symbol or labeled call per denominator, in family order.
    ///
    /// Returns
    /// -------
    /// list[tuple[Expression, Expression]]
    ///     Simultaneous replacement pairs for an independent, complete family.
    fn scalar_product_rules(
        &self,
        labels: Vec<PythonExpression>,
    ) -> PyResult<Vec<(PythonExpression, PythonExpression)>> {
        let labels = labels
            .into_iter()
            .map(|label| label.expr)
            .collect::<Vec<_>>();
        Ok(self
            .inner
            .scalar_product_rules(&labels)
            .map_err(error::integral_family)?
            .into_iter()
            .map(|(product, value)| (product.into(), value.into()))
            .collect())
    }

    /// Rewrite a numerator in inverse-propagator variables and expand it.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralFamily`` class example:
    ///
    /// >>> numerator = kin.scalar_product(k, p)**2
    /// >>> reduced = family.rewrite_numerator(numerator, [d1, d2])
    ///
    /// Parameters
    /// ----------
    /// numerator : Expression
    ///     Scalar numerator after tensor reduction and momentum routing.
    /// labels : list[Expression]
    ///     One distinct symbol or labeled call per denominator, in family order.
    fn rewrite_numerator(
        &self,
        numerator: &PythonExpression,
        labels: Vec<PythonExpression>,
    ) -> PyResult<PythonExpression> {
        let labels = labels
            .into_iter()
            .map(|label| label.expr)
            .collect::<Vec<_>>();
        self.inner
            .rewrite_numerator(&numerator.expr, &labels)
            .map(Into::into)
            .map_err(error::integral_family)
    }
}

/// A verified real loop-momentum shift and propagator embedding.
///
/// Obtain a mapping from ``IntegralFamily.find_mapping`` or ``mapping_to``.
/// The loop determinant has absolute value one, so the loop integration measure
/// is unchanged. Numerator substitutions use the same scalar-product rules
/// that verified the propagator identities.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> D, k, p, s = S("D", "k", "p", "s")
/// >>> d1, d2, x1, x2 = S("d1", "d2", "x1", "x2")
/// >>> kin = hep.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)
/// >>> denominators = [kin.scalar_product(k, k), kin.scalar_product(k-p, k-p)]
/// >>> family = hep.IntegralFamily([k], [p], denominators, kinematics=kin)
/// >>> m2 = S("m2")
/// >>> denominators = [denominators[0] - m2, denominators[1]]
/// >>> source = hep.IntegralFamily([k], [p], denominators, kinematics=kin)
/// >>> target = hep.IntegralFamily([k], [p], list(reversed(denominators)), kinematics=kin)
/// >>> mapping = source.find_mapping(target)
/// >>> assert mapping is not None
/// >>> target_powers = mapping.map_powers([1, 2])
/// >>> assert target_powers == [2, 1]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "IntegralMapping",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyIntegralMapping {
    inner: IntegralMapping,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyIntegralMapping {
    /// Source loop momentum names paired with their target-coordinate images.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralMapping`` class example:
    ///
    /// >>> loop_substitutions = mapping.momentum_rules
    #[getter]
    fn momentum_rules(&self) -> Vec<(PythonExpression, PythonExpression)> {
        self.inner
            .momentum_rules()
            .iter()
            .cloned()
            .map(|(p, q)| (p.into(), q.into()))
            .collect()
    }

    /// Zero-based target denominator index for each source denominator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralMapping`` class example:
    ///
    /// >>> mapped_powers = mapping.map_powers([1, 2])
    /// >>> slot_map = mapping.denominator_map
    #[getter]
    fn denominator_map(&self) -> Vec<usize> {
        self.inner.denominator_map().to_vec()
    }

    /// Reorder signed powers and insert zeros for unused target propagators.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``IntegralMapping`` class example:
    ///
    /// >>> powers = mapping.map_powers([1, 2])
    ///
    /// Parameters
    /// ----------
    /// powers : list[int]
    ///     One signed power per source denominator, in source order.
    fn map_powers(&self, powers: Vec<i32>) -> PyResult<Vec<i32>> {
        self.inner
            .map_powers(&powers)
            .map_err(error::integral_family)
    }

    /// Apply the verified scalar-product substitutions simultaneously.
    ///
    /// Examples
    /// --------
    /// Using the setup in ``IntegralMapping``:
    ///
    /// >>> transformed = mapping.apply(kin.scalar_product(k, p))
    ///
    /// Parameters
    /// ----------
    /// expression : Expression
    ///     Scalar expression in compact Spenso dot notation; contract tensors first.
    fn apply(&self, expression: &PythonExpression) -> PythonExpression {
        self.inner.apply(&expression.expr).into()
    }
}

/// A formal propagator permutation determined from Symanzik polynomials.
///
/// Unlike a momentum mapping, this supplies no tensor-numerator substitution
/// or check of integration prescriptions. Signed propagator powers can be
/// reordered; scalar numerators must first be expressed in family coordinates.
///
/// Examples
/// --------
/// >>> from symbolica import S, E
/// >>> from symbolica.community import hep
/// >>> D, k, p, s = S("D", "k", "p", "s")
/// >>> d1, d2, x1, x2 = S("d1", "d2", "x1", "x2")
/// >>> kin = hep.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)
/// >>> denominators = [kin.scalar_product(k, k), kin.scalar_product(k-p, k-p)]
/// >>> family = hep.IntegralFamily([k], [p], denominators, kinematics=kin)
/// >>> m2 = S("m2")
/// >>> denominators = [denominators[0] - m2, denominators[1]]
/// >>> source = hep.IntegralFamily([k], [p], denominators, kinematics=kin)
/// >>> target = hep.IntegralFamily([k], [p], list(reversed(denominators)), kinematics=kin)
/// >>> mapping = source.parametric_mapping(target, [x1, x2])
/// >>> assert mapping is not None
/// >>> target_powers = mapping.map_powers([1, 2])
/// >>> assert target_powers == [2, 1]
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "PropagatorMapping",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyPropagatorMapping {
    inner: PropagatorMapping,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyPropagatorMapping {
    /// Zero-based target denominator index for each source denominator.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``PropagatorMapping`` class example:
    ///
    /// >>> slot_map = mapping.denominator_map
    /// >>> mapped_powers = mapping.map_powers([1, 2])
    #[getter]
    fn denominator_map(&self) -> Vec<usize> {
        self.inner.denominator_map().to_vec()
    }

    /// Reorder signed propagator powers into target-family order.
    ///
    /// Examples
    /// --------
    /// Using the setup in the ``PropagatorMapping`` class example:
    ///
    /// >>> reordered = mapping.map_powers([1, -2])
    ///
    /// Parameters
    /// ----------
    /// powers : list[int]
    ///     One signed power per source denominator.
    fn map_powers(&self, powers: Vec<i32>) -> PyResult<Vec<i32>> {
        self.inner
            .map_powers(&powers)
            .map_err(error::integral_family)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyIntegralFamily>()?;
    module.add_class::<PyIntegralMapping>()?;
    module.add_class::<PyPropagatorMapping>()
}
