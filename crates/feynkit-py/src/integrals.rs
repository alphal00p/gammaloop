use feynkit_graph::{IntegralFamily, IntegralMapping, PropagatorMapping};
use pyo3::prelude::*;
use symbolica::api::python::PythonExpression;

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{error, kinematics::PyKinematics};

/// An ordered family of affine inverse propagators in loop scalar products.
///
/// Construct denominators with ``Kinematics.scalar_product``. External momenta
/// must form an independent basis. Rank, completion, and numerator rewriting
/// use Symbolica's exact linear algebra. Verified momentum shifts are available
/// through ``mapping_to`` and ``find_mapping``. Parametric scaling certificates
/// detect scaleless sectors in dimensional regularization. Integration
/// prescriptions and IBP reduction are separate operations.
///
/// Examples
/// --------
/// >>> kin = fk.Kinematics(D, momenta=[k, p]).with_scalar_product(p, p, s)
/// >>> denominators = [kin.scalar_product(k, k), kin.scalar_product(k-p, k-p)]
/// >>> family = fk.IntegralFamily([k], [p], denominators, kinematics=kin)
/// >>> assert family.is_complete and family.is_independent
/// >>> reduced = family.rewrite_numerator(kin.scalar_product(k, p)**2, [d1, d2])
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyIntegralFamily {
    /// Compute the independent loop scalar products and denominator rank.
    ///
    /// Examples
    /// --------
    /// >>> family = fk.IntegralFamily([k], [p], denominators, kinematics=kin)
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
    #[getter]
    fn kinematics(&self) -> PyKinematics {
        PyKinematics {
            inner: self.inner.kinematics().clone(),
        }
    }

    /// Independent external momentum names in their original order.
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
    #[getter]
    fn rank(&self) -> usize {
        self.inner.rank()
    }

    /// Whether the inverse propagators span every loop scalar product.
    #[getter]
    fn is_complete(&self) -> bool {
        self.inner.is_complete()
    }

    /// Whether no denominator can be eliminated by an affine relation.
    #[getter]
    fn is_independent(&self) -> bool {
        self.inner.is_independent()
    }

    /// Return a compact summary of the family and its scalar-product rank.
    ///
    /// Examples
    /// --------
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
    /// Leave ``family`` as the final expression in a notebook cell.
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
    /// IPython uses this representation when rich HTML output is unavailable.
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
    /// >>> completed = family.complete()
    /// >>> assert completed.is_complete
    /// >>> completed = family.complete(candidates=other_family.denominators)
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
    /// >>> terms = family.partial_fraction([1, 1])
    /// >>> for coefficient, powers in terms:
    /// ...     print(coefficient, powers)
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
    /// >>> mapping = source.mapping_to(target, [l - p])
    /// >>> assert mapping is not None
    /// >>> target_powers = mapping.map_powers([1, 2])
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
    /// >>> mapping = source.find_mapping(target)
    /// >>> if mapping is not None:
    /// ...     transformed = mapping.apply(scalar_numerator)
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

    /// Select the positive-power propagators of an integral's sector.
    ///
    /// Zero and negative powers are omitted. Loop variables and external
    /// kinematics are retained. This selects sector support; it does not remove
    /// numerator factors algebraically from the original integrand.
    ///
    /// Examples
    /// --------
    /// >>> sector = family.sector([1, 2, -1])
    /// >>> assert len(sector.denominators) == 2
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
    /// >>> sector = family.sector([1, 1, 0])
    /// >>> weights = sector.scaleless_scaling([x1, x2])
    /// >>> if weights is not None:
    /// ...     print("Scaleless in dimensional regularization", weights)
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
    /// >>> mapping = source.parametric_mapping(target, [x1, x2])
    /// >>> if mapping is not None:
    /// ...     target_powers = mapping.map_powers([1, 2])
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
/// >>> mapping = source.find_mapping(target)
/// >>> assert mapping is not None
/// >>> print(mapping.momentum_rules, mapping.denominator_map)
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
    #[getter]
    fn denominator_map(&self) -> Vec<usize> {
        self.inner.denominator_map().to_vec()
    }

    /// Reorder signed powers and insert zeros for unused target propagators.
    ///
    /// Examples
    /// --------
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
    /// >>> transformed_numerator = mapping.apply(numerator)
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
/// >>> mapping = source.parametric_mapping(target, [x1, x2])
/// >>> if mapping is not None:
/// ...     print(mapping.denominator_map, mapping.map_powers([1, 2]))
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
    #[getter]
    fn denominator_map(&self) -> Vec<usize> {
        self.inner.denominator_map().to_vec()
    }

    /// Reorder signed propagator powers into target-family order.
    ///
    /// Examples
    /// --------
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
