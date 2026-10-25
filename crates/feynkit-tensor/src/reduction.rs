//! Spenso-facing covariant tensor reduction of loop numerators.
//!
//! The reducer consumes admitted symbolic tensors whose declared heads use
//! Spenso slots, for example
//! `K(1, spenso::mink(D, mu))`.  Integrated vectors are selected by their head
//! or by their compact indexed-free form `K(1, spenso::mink(D))`.  Results use
//! the same ecosystem notation: scalar contractions are emitted as
//! `spenso::dot(...)` and uncontracted Lorentz pairs as `spenso::g(...)`.
//! Explicit Spenso metric chains are first contracted with Idenso.  In
//! particular, users may pass the natural
//! `g(mu,nu) K(mu) K(nu)` representation without a separate normalization
//! step; genuinely free indices are retained in the output metric slots.
//!
//! Repeated vectors are represented by contraction orbits.  For vector
//! multiplicities `N_i`, an orbit stores the upper-triangular exponents
//! `alpha_ij` in
//! `prod_(i <= j) dot(p_i, p_j)^alpha_ij`.  This is the symmetry reduction used
//! by modern high-rank tensor projectors: a rank-20 product of one loop vector
//! projected on one external vector has one orbit instead of `19!!` labeled
//! pairings.  The internal/projector orbit construction follows Sections
//! 2.6--2.7 of Goode *et al.*,
//! [arXiv:2408.05137](https://arxiv.org/abs/2408.05137).
//!
//! Partially symmetric, fully contracted projectors also avoid the Cartesian
//! product of labeled pairings.  If `S_A` and `S_B` are raw matching-to-orbit
//! incidence matrices, `G` is the metric-pairing Gram matrix, and `W=G^-1`,
//! the desired coefficients are `C=S_A^T W S_B`.  On an orbit space,
//! `G S_A=S_A H_A`, so exactly
//! `C=H_A^-T N=N H_B^-1`, with `N=S_A^T S_B`.  FeynKit builds `N` by a
//! joint-color dynamic program, builds `H` from alternating cycles, and solves
//! whichever side has fewer contraction orbits.  Only residual free-index
//! output uses the explicitly budgeted labeled-pairing path.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::mem;

use feynkit_graph::FeynmanDiagram;
use idenso::tensor::{ContractSettings, SymbolicTensor};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use spenso::{
    shadowing::TensorCollectFilter,
    structure::{
        partial::{PartialStructure, PartialStructureExt},
        representation::{BaseRepName, LibraryRep, Minkowski},
        slot::IsAbstractSlot,
    },
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::utils::resources;
use symbolica::utils::resources::{ResourceBudget, ResourceLimitError, ResourceLimits};
use thiserror::Error;

mod blocks;
mod external;

use crate::{
    OrthogonalWeingarten, WeingartenError,
    orbit_projector::{OrbitProjector, OrbitProjectorError},
};

/// One invariant contraction class of repeated vectors.
///
/// If `vectors = [k1, k2]`, then `alpha = [a11, a12, a22]` represents
/// `dot(k1,k1)^a11 dot(k1,k2)^a12 dot(k2,k2)^a22`.  The degree constraints are
/// `2 a_ii + sum_(j != i) a_ij = N_i`, where `N_i` is the corresponding entry
/// of [`Self::multiplicities`].
///
/// # Examples
///
/// A rank-six numerator containing two powers of one loop momentum and four
/// of another has only two contraction orbits, rather than 15 labeled
/// pairings.
///
/// ```no_run
/// use feynkit_tensor::TensorReducer;
/// use symbolica::{parse, symbol};
///
/// let _ = (spenso::vector_symbol!("k"), spenso::vector_symbol!("q"), spenso::vector_symbol!("p"));
/// let numerator = parse!(
///     "k(spenso::mink(D,m1))*k(spenso::mink(D,m2))
///      *q(spenso::mink(D,m3))*q(spenso::mink(D,m4))
///      *q(spenso::mink(D,m5))*q(spenso::mink(D,m6))
///      *p(spenso::mink(D,m1))*p(spenso::mink(D,m2))
///      *p(spenso::mink(D,m3))*p(spenso::mink(D,m4))
///      *p(spenso::mink(D,m5))*p(spenso::mink(D,m6))"
/// );
/// let reduced = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .with_integrated_head(symbol!("q"))
///     .reduce(&idenso::tensor::SymbolicTensor::infer(numerator)?)?;
/// let orbit = reduced.terms()[0].integrated_orbit().unwrap();
/// println!("alpha = {:?}", orbit.alpha());
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct ContractionOrbit {
    vectors: Box<[Atom]>,
    multiplicities: Box<[usize]>,
    alpha: Box<[u16]>,
    labeled_pairings: u128,
}

impl ContractionOrbit {
    /// Distinct compact vectors in canonical order.
    pub fn vectors(&self) -> &[Atom] {
        &self.vectors
    }

    /// Number of occurrences of every distinct vector.
    pub fn multiplicities(&self) -> &[usize] {
        &self.multiplicities
    }

    /// Upper-triangular `alpha_ij` entries in row-major order.
    pub fn alpha(&self) -> &[u16] {
        &self.alpha
    }

    /// Number of labeled perfect matchings represented by this orbit.
    pub fn labeled_pairings(&self) -> u128 {
        self.labeled_pairings
    }

    /// Convert the orbit to a product of standard Spenso dot products.
    ///
    /// # Errors
    /// Returns a resource error if a configured process-memory or elapsed-time
    /// allowance is exceeded during materialization.
    pub fn invariant(&self) -> Result<Atom, TensorReductionError> {
        resources::checkpoint()?;
        let mut invariant = Atom::one();
        let mut alpha = 0;
        for left in 0..self.vectors.len() {
            for right in left..self.vectors.len() {
                resources::checkpoint()?;
                let exponent = self.alpha[alpha];
                alpha += 1;
                if exponent == 0 {
                    continue;
                }
                let dot = dot(&self.vectors[left], &self.vectors[right]);
                invariant *= if exponent == 1 {
                    dot
                } else {
                    dot.pow(Atom::num(i64::from(exponent)))
                };
            }
        }
        resources::checkpoint()?;
        Ok(invariant)
    }
}

/// One compact contribution to a tensor reduction.
///
/// A term keeps the exact projector coefficient separate from its scalar
/// invariant so callers can inspect symmetry orbits before materializing the
/// product with [`Self::expression`].
///
/// # Examples
///
/// ```no_run
/// use feynkit_tensor::TensorReducer;
/// use symbolica::{parse, symbol};
///
/// let _ = (spenso::vector_symbol!("k"), spenso::vector_symbol!("q"), spenso::vector_symbol!("p"));
/// let numerator = parse!(
///     "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))
///      *p(spenso::mink(D,mu))*p(spenso::mink(D,nu))"
/// );
/// let reduction = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .reduce(&idenso::tensor::SymbolicTensor::infer(numerator)?)?;
/// let scalar_graph_numerator = reduction.terms()[0].expression()?;
/// println!("{scalar_graph_numerator}");
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug)]
pub struct TensorReductionTerm {
    coefficient: Atom,
    tensor: SymbolicTensor<PartialStructure>,
    integrated_orbit: Option<ContractionOrbit>,
    projector_orbit: Option<ContractionOrbit>,
    budget: ResourceBudget,
}

impl TensorReductionTerm {
    /// Exact rational coefficient, including orbit multiplicities.
    pub fn coefficient(&self) -> &Atom {
        &self.coefficient
    }

    /// Scalar invariant or residual metric tensor multiplying the coefficient.
    pub fn tensor(&self) -> &SymbolicTensor<PartialStructure> {
        &self.tensor
    }

    /// Symmetry orbit of contractions among integrated vectors, when retained.
    pub fn integrated_orbit(&self) -> Option<&ContractionOrbit> {
        self.integrated_orbit.as_ref()
    }

    /// Symmetry orbit of the external projector vectors, when fully contracted.
    pub fn projector_orbit(&self) -> Option<&ContractionOrbit> {
        self.projector_orbit.as_ref()
    }

    /// Apply the scalar weight without discarding the retained tensor interface.
    pub fn expression(&self) -> Result<SymbolicTensor<PartialStructure>, TensorReductionError> {
        self.budget.run(|| {
            resources::checkpoint()?;
            let expression = self
                .tensor
                .with_checked_expression(&self.coefficient * self.tensor.expression())?;
            resources::checkpoint()?;
            Ok(expression)
        })
    }
}

/// Symmetry-aware result of reducing a tensor expression.
///
/// # Examples
///
/// Materialize compact symmetry orbits as a standard Spenso/Symbolica scalar
/// expression suitable for the next stage of a QFT calculation.
///
/// ```no_run
/// use feynkit_tensor::TensorReducer;
/// use symbolica::{parse, symbol};
///
/// let _ = (spenso::vector_symbol!("k"), spenso::vector_symbol!("q"), spenso::vector_symbol!("p"));
/// let numerator = parse!(
///     "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))
///      *p(spenso::mink(D,mu))*p(spenso::mink(D,nu))"
/// );
/// let scalar = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .reduce(&idenso::tensor::SymbolicTensor::infer(numerator)?)?
///     .into_expression()?;
/// println!("{scalar}");
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug)]
pub struct TensorReduction {
    terms: Vec<TensorReductionTerm>,
    zero: SymbolicTensor<PartialStructure>,
    max_rank: usize,
    fully_contracted: bool,
    budget: ResourceBudget,
}

impl TensorReduction {
    /// Compact reduction terms before their final sum is materialized.
    pub fn terms(&self) -> &[TensorReductionTerm] {
        &self.terms
    }

    /// Largest residual projector rank after direct Einstein contractions.
    pub fn max_rank(&self) -> usize {
        self.max_rank
    }

    /// Whether the result contains no dangling Minkowski slots.
    pub fn is_fully_contracted(&self) -> bool {
        self.fully_contracted
    }

    /// Materialize all terms as one standard Spenso/Symbolica expression.
    ///
    /// # Examples
    ///
    /// The returned atom includes every exact projector coefficient and is
    /// ready to use as a scalar-graph numerator.
    ///
    /// ```no_run
    /// use feynkit_tensor::TensorReduction;
    /// use symbolica::atom::Atom;
    ///
    /// fn scalar_numerator(reduction: &TensorReduction) -> Result<Atom, feynkit_tensor::TensorReductionError> {
    ///     Ok(reduction.expression()?.into_expression())
    /// }
    /// ```
    pub fn expression(&self) -> Result<SymbolicTensor<PartialStructure>, TensorReductionError> {
        self.budget.run(|| {
            resources::checkpoint()?;
            let expression = self.terms.iter().try_fold(self.zero.clone(), |sum, term| {
                resources::checkpoint()?;
                Ok::<_, TensorReductionError>(sum.try_add(&term.expression()?)?)
            })?;
            resources::checkpoint()?;
            Ok(expression)
        })
    }

    /// Consume the compact result and materialize its expression.
    pub fn into_expression(self) -> Result<SymbolicTensor<PartialStructure>, TensorReductionError> {
        self.budget.run(|| {
            resources::checkpoint()?;
            let expression = self.terms.into_iter().try_fold(self.zero, |sum, term| {
                resources::checkpoint()?;
                Ok::<_, TensorReductionError>(sum.try_add(&term.expression()?)?)
            })?;
            resources::checkpoint()?;
            Ok(expression)
        })
    }
}

/// A tensor reducer configured by dimension and integrated-vector selectors.
///
/// The primary high-rank path is a fully contracted expression.  For example,
/// selecting `K` as integrated in
/// `K(mu) K(nu) p(mu) p(nu)` directly produces a scalar polynomial in
/// `dot(K,K)` and `dot(p,p)`.  Internal and projector repetitions are both used
/// to avoid generating equivalent labeled pairings.
///
/// # Examples
///
/// Reduce a one-loop rank-two vacuum numerator in dimensional regularization.
/// The metric is normalized automatically; no preliminary Idenso call is
/// required.
///
/// ```no_run
/// use feynkit_tensor::TensorReducer;
/// use symbolica::{parse, symbol};
///
/// let _ = (spenso::vector_symbol!("k"), spenso::vector_symbol!("q"), spenso::vector_symbol!("p"));
/// let numerator = parse!(
///     "spenso::g(spenso::mink(D,mu),spenso::mink(D,nu))
///      *k(spenso::mink(D,mu))*k(spenso::mink(D,nu))"
/// );
/// let reduced = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .reduce(&idenso::tensor::SymbolicTensor::infer(numerator)?)?;
/// println!("{}", reduced.expression()?);
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug)]
pub struct TensorReducer {
    dimension: Atom,
    integrated_heads: BTreeSet<String>,
    integrated_vectors: BTreeSet<Atom>,
    external_vectors: BTreeSet<Atom>,
    pairing_limit: Option<u128>,
    pairing_product_limit: Option<u128>,
    output_term_limit: Option<usize>,
    resource_limits: ResourceLimits,
}

impl TensorReducer {
    /// Construct a reducer in `dimension` with no integrated vectors selected.
    ///
    /// `dimension` must exactly match the source dimension of every input
    /// `spenso::mink(dimension,index)` slot, including reversibly encoded
    /// compound dimensions. Thus generated four-dimensional
    /// FeynKit rules use `4`, while a symbolic R-operation numerator whose
    /// slots carry `D` uses `D`. Mixed high-rank projectors should retain a
    /// symbolic dimension until after reduction: a fixed positive integer
    /// dimension below half the tensor rank makes the universal metric basis
    /// singular. The all-equal isotropic fast path remains available there.
    ///
    /// Count diagnostics and process-memory/elapsed-time limits are disabled
    /// by default. Configure resources with [`Self::with_resources`]. Matching
    /// counters and exponent representations remain checked for overflow.
    ///
    /// Add selectors with [`Self::with_integrated_head`] or
    /// [`Self::with_integrated_vector`].
    pub fn new(dimension: Atom) -> Self {
        Self {
            dimension,
            integrated_heads: BTreeSet::new(),
            integrated_vectors: BTreeSet::new(),
            external_vectors: BTreeSet::new(),
            pairing_limit: None,
            pairing_product_limit: None,
            output_term_limit: None,
            resource_limits: ResourceLimits::default(),
        }
    }

    /// Construct a reducer for every canonical diagram `gammalooprs::Q` vector.
    ///
    /// This constructor stores the qualified head name without registering a
    /// bare Symbolica symbol, so FeynKit's tensor tags remain initialization
    /// order independent.
    pub fn feynkit(dimension: Atom) -> Self {
        Self::new(dimension).with_integrated_head(feynkit_graph::momentum_symbol())
    }

    /// Select every vector with this Symbolica function head.
    pub fn with_integrated_head(mut self, head: Symbol) -> Self {
        self.integrated_heads.insert(head.get_name().to_owned());
        self
    }

    /// Select every vector with a qualified Symbolica head name.
    pub fn with_integrated_head_name(mut self, head: impl Into<String>) -> Self {
        self.integrated_heads.insert(head.into());
        self
    }

    /// Admit one exact compact source vector, such as `K(1,spenso::mink(D))`.
    /// Shared source admission normalizes encoded dimensions and index labels
    /// before the selector is compared with retained native tensor expressions.
    pub fn with_integrated_vector(mut self, vector: Atom) -> Result<Self, TensorReductionError> {
        let vector = feynkit_graph::NumeratorInterface::admit(vector, None)?;
        self.integrated_vectors.insert(vector.into_expression());
        Ok(self)
    }

    /// Add an independent external momentum appearing in the denominators.
    ///
    /// Pass compact Spenso vectors such as `p(spenso::mink(D))`. The reducer
    /// resolves loop vectors along this basis and averages only their
    /// transverse components, reusing the vacuum projector in `D - n`
    /// dimensions. The external Gram matrix must be invertible. Include an
    /// auxiliary vector for a degenerate basis, such as one lightlike momentum.
    /// The resulting scalar products are left for integral-family reduction.
    /// External vectors take precedence over integrated-head selectors.
    pub fn with_external_vector(mut self, vector: Atom) -> Result<Self, TensorReductionError> {
        let vector = feynkit_graph::NumeratorInterface::admit(vector, None)?;
        self.external_vectors.insert(vector.into_expression());
        Ok(self)
    }

    /// Opt into a labeled-pairing diagnostic for the generic/free-index path.
    ///
    /// Disabled by default. Symmetric contraction-orbit fast paths do not
    /// consume this diagnostic allowance.
    pub fn with_pairing_limit(mut self, limit: u128) -> Self {
        self.pairing_limit = Some(limit);
        self
    }

    /// Opt into a matching-product diagnostic for residual-free-index output.
    /// Disabled by default.
    pub fn with_pairing_product_limit(mut self, limit: u128) -> Self {
        self.pairing_product_limit = Some(limit);
        self
    }

    /// Opt into a diagnostic bound on materialized invariant terms.
    /// Disabled by default.
    pub fn with_output_term_limit(mut self, limit: usize) -> Self {
        self.output_term_limit = Some(limit);
        self
    }

    /// Configure shared process-memory and elapsed-time limits for a request.
    /// Nested resource scopes retain any stricter limits set by their caller.
    ///
    /// ```
    /// use feynkit_tensor::TensorReducer;
    /// use symbolica::parse;
    /// use symbolica::utils::resources::ResourceLimits;
    ///
    /// let reducer = TensorReducer::new(parse!("D"))
    ///     .with_resources(ResourceLimits::default());
    /// assert_eq!(reducer.resource_limits(), ResourceLimits::default());
    /// ```
    pub fn with_resources(mut self, limits: ResourceLimits) -> Self {
        self.resource_limits = limits;
        self
    }

    /// Return the process-memory and elapsed-time configuration.
    pub fn resource_limits(&self) -> ResourceLimits {
        self.resource_limits
    }

    /// Dimension used by Spenso slots and projector coefficients.
    pub fn dimension(&self) -> &Atom {
        &self.dimension
    }

    /// Admit a newly constructed numerator at an external source boundary.
    ///
    /// The expression must already use native Spenso representation syntax,
    /// with compound dimensions and indices encoded by their source adapter.
    /// This method does not cook source notation or reduce a tensor. Native
    /// admission remains authoritative; rejected inputs receive the reducer's
    /// established index and power diagnostics when its parser can refine them.
    pub fn admit(
        &self,
        expression: Atom,
    ) -> Result<SymbolicTensor<PartialStructure>, TensorReductionError> {
        ResourceBudget::new(self.resource_limits).run(|| {
            resources::checkpoint()?;
            let source = SymbolicTensor::infer(expression.clone())
                .map_err(|error| self.input_error(expression.as_view(), error))?;
            resources::checkpoint()?;
            Ok(source)
        })
    }

    fn input_error(
        &self,
        expression: AtomView<'_>,
        error: idenso::tensor::inference::TensorInferenceError,
    ) -> TensorReductionError {
        self.parse_monomial(expression)
            .err()
            .unwrap_or_else(|| TensorReductionError::TensorAlgebra(error.to_string()))
    }

    /// Reduce a tensor expression to compact invariant terms.
    ///
    /// Explicit Spenso metric contractions are normalized first with Idenso.
    /// This contracts chains such as `g(mu,nu) K(mu) K(nu)` while preserving
    /// any residual free index in its full `spenso::mink(D,index)` slot.  Each
    /// selected Lorentz sector is then treated independently; its scalar and
    /// foreign-representation coefficient remains factored. Compact Spenso dots between an
    /// integrated vector and a spectator are projected directly, including
    /// nonnegative integer powers, without allocating dummy indices. Dots
    /// between integrated vectors or with declared external basis vectors stay
    /// scalar invariants. Negative or noninteger powers of spectator contractions
    /// are not polynomial tensor numerators and are rejected.
    /// Odd-rank vacuum tensors vanish;
    /// already-contracted integrated pairs become Spenso dots before the
    /// projector is built. A single index occurrence inside an otherwise
    /// opaque tensor is retained as an external projector leg; expressions
    /// with more than two occurrences of one Minkowski index are rejected as
    /// ambiguous rather than silently dropping a contraction.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feynkit_tensor::TensorReducer;
    /// use symbolica::{parse, symbol};
    ///
    /// let _ = (spenso::vector_symbol!("k"), spenso::vector_symbol!("q"), spenso::vector_symbol!("p"));
    /// let numerator = parse!(
    ///     "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))
    ///      *p(spenso::mink(D,mu))*p(spenso::mink(D,nu))"
    /// );
    /// let reduction = TensorReducer::new(parse!("D"))
    ///     .with_integrated_head(symbol!("k"))
    ///     .reduce(&idenso::tensor::SymbolicTensor::infer(numerator)?)?;
    /// assert!(reduction.is_fully_contracted());
    /// # Ok::<(), feynkit_tensor::TensorReductionError>(())
    /// ```
    pub fn reduce(
        &self,
        source: &SymbolicTensor<PartialStructure>,
    ) -> Result<TensorReduction, TensorReductionError> {
        ResourceBudget::new(self.resource_limits).run(|| {
            let reduction = self.reduce_scoped(source)?;
            resources::checkpoint()?;
            Ok(reduction)
        })
    }

    fn reduce_scoped(
        &self,
        source: &SymbolicTensor<PartialStructure>,
    ) -> Result<TensorReduction, TensorReductionError> {
        resources::checkpoint()?;
        if self.integrated_heads.is_empty() && self.integrated_vectors.is_empty() {
            return Err(TensorReductionError::NoIntegratedVectorsSelected);
        }
        let external = self.external_projector()?;
        // The shared collector distributes only selected Lorentz structures.
        // Scalar and foreign-representation coefficients retain their original
        // factorization, and a capped selection cannot publish partial rows.
        let normalized = source
            .contract(ContractSettings {
                rank_one: false,
                collect_chains: false,
                collect_traces: false,
                expand: false,
                ..Default::default()
            })
            // Contract validates through the shared cached owner; diagnostic
            // refinement never accepts a source that native validation rejects.
            .map_err(|error| self.input_error(source.expression().as_view(), error))?;
        resources::checkpoint()?;
        let normalized = normalized
            .to_dots()
            .map_err(|error| TensorReductionError::TensorAlgebra(error.to_string()))?;
        resources::checkpoint()?;
        if Self::has_tensor_sum(normalized.expression().as_view())? {
            return self.reduce_blocks(&normalized, external.as_ref());
        }
        let rows = normalized
            .coefficient_list(TensorCollectFilter::Reps([LibraryRep::from(Minkowski {})]))
            .map_err(|error| TensorReductionError::TensorAlgebra(error.to_string()))?;
        resources::checkpoint()?;
        if let Some(limit) = self.output_term_limit
            && rows.len() > limit
        {
            return Err(TensorReductionError::OutputLimit {
                terms: rows.len(),
                limit,
            });
        }

        let mut terms = Vec::new();
        let mut max_rank = 0;
        for (selected, coefficient) in rows {
            resources::checkpoint()?;
            let monomial = self.parse_monomial(selected.expression().as_view())?;
            max_rank = max_rank.max(monomial.integrated.len());
            let projected = if let Some(external) = &external {
                self.reduce_external_monomial(monomial, external)?
            } else if monomial.integrated.len() % 2 == 0 {
                self.reduce_monomial(monomial)?
            } else {
                Vec::new()
            };
            for mut term in projected {
                resources::checkpoint()?;
                let tensor = coefficient.multiply(&term.tensor)?;
                term.tensor = source.with_checked_expression(tensor.expression().clone())?;
                terms.push(term);
            }
            if let Some(limit) = self.output_term_limit
                && terms.len() > limit
            {
                return Err(TensorReductionError::OutputLimit {
                    terms: terms.len(),
                    limit,
                });
            }
        }
        let zero = source.with_checked_expression(Atom::Zero)?;
        let fully_contracted = zero
            .structure()
            .logical_slots()
            .iter()
            .all(|slot| slot.rep().rep != LibraryRep::from(Minkowski {}));
        resources::checkpoint()?;
        Ok(TensorReduction {
            terms,
            zero,
            max_rank,
            fully_contracted,
            budget: ResourceBudget::capture_current(),
        })
    }

    fn parse_monomial(&self, term: AtomView<'_>) -> Result<TensorMonomial, TensorReductionError> {
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            factor => vec![factor],
        };
        let mut scalar_factors = Vec::new();
        let mut indexed_vectors = Vec::new();
        let mut opaque_factors = Vec::new();
        let mut integrated = Vec::new();
        let mut outside = Vec::new();

        for factor in factors {
            resources::checkpoint()?;
            let (base, exponent) = match factor {
                AtomView::Pow(power) => (power.get_base(), power.get_exp().to_owned()),
                base => (base, Atom::one()),
            };
            if let AtomView::Fun(function) = base
                && function.get_symbol() == SPENSO_TAG.dot
                && function.get_nargs() == 2
            {
                // Reuse the indexed-vector parser solely to validate compact
                // slots and classify the two endpoints. No dummy index enters
                // the numerator: the contraction is already explicit in dot.
                let vectors = function
                    .iter()
                    .map(|argument| {
                        let compact = argument.to_owned();
                        let indexed = indexed_vector(&compact, &Atom::Zero);
                        self.indexed_vector(indexed.as_view())
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if let [Some(left), Some(right)] = vectors.as_slice() {
                    let pair = match (self.is_integrated(left), self.is_integrated(right)) {
                        (true, false) => Some((left, right)),
                        (false, true) => Some((right, left)),
                        _ => None,
                    };
                    if let Some((internal, external)) = pair
                        && !self.external_vectors.contains(&external.compact)
                    {
                        let count = i64::try_from(exponent.as_view())
                            .ok()
                            .and_then(|count| usize::try_from(count).ok())
                            .ok_or(TensorReductionError::InvalidVectorPower(exponent))?;
                        let bytes = count
                            .checked_mul(mem::size_of::<Atom>() + mem::size_of::<Outside>())
                            .and_then(|bytes| u64::try_from(bytes).ok())
                            .ok_or(TensorReductionError::CounterOverflow)?;
                        resources::check_allocation(bytes)?;
                        integrated.try_reserve(count).map_err(|error| {
                            TensorReductionError::TensorAlgebra(error.to_string())
                        })?;
                        outside.try_reserve(count).map_err(|error| {
                            TensorReductionError::TensorAlgebra(error.to_string())
                        })?;
                        integrated.extend(std::iter::repeat_n(internal.compact.clone(), count));
                        outside.extend(std::iter::repeat_n(
                            Outside::Vector(external.compact.clone()),
                            count,
                        ));
                        continue;
                    }
                }
            }
            if let AtomView::Pow(power) = factor
                && let Some(vector) = self.indexed_vector(power.get_base())?
            {
                let exponent = i64::try_from(power.get_exp()).map_err(|_| {
                    TensorReductionError::InvalidVectorPower(power.get_exp().to_owned())
                })?;
                let exponent = usize::try_from(exponent).map_err(|_| {
                    TensorReductionError::InvalidVectorPower(power.get_exp().to_owned())
                })?;
                let integrated = self.is_integrated(&vector);
                if exponent > 2 {
                    return Err(TensorReductionError::AmbiguousIndex {
                        index: vector.index,
                        integrated: if integrated { exponent } else { 0 },
                        outside: if integrated { 0 } else { exponent },
                    });
                }
                indexed_vectors.extend(std::iter::repeat_n(vector, exponent));
                continue;
            }
            if let Some(vector) = self.indexed_vector(factor)? {
                indexed_vectors.push(vector);
            } else {
                let factor = factor.to_owned();
                opaque_factors.push(factor.clone());
                scalar_factors.push(factor);
            }
        }

        let mut by_index: BTreeMap<Atom, IndexOccupants> = BTreeMap::new();
        for vector in indexed_vectors {
            resources::checkpoint()?;
            let integrated = self.is_integrated(&vector);
            let occupants =
                by_index
                    .entry(vector.index.clone())
                    .or_insert_with(|| IndexOccupants {
                        slot: vector.slot.clone(),
                        integrated: Vec::new(),
                        outside: Vec::new(),
                    });
            if integrated {
                occupants.integrated.push(vector.compact);
            } else {
                occupants.outside.push(vector.compact);
            }
        }

        for (index, occupants) in by_index {
            resources::checkpoint()?;
            let opaque_occurrences =
                opaque_factors
                    .iter()
                    .try_fold(0_usize, |occurrences, factor| {
                        Ok::<usize, TensorReductionError>(occurrences.saturating_add(
                            count_minkowski_index(factor.as_view(), &self.dimension, &index, 1)?,
                        ))
                    })?;
            let total_occurrences = occupants
                .integrated
                .len()
                .saturating_add(occupants.outside.len())
                .saturating_add(opaque_occurrences);
            if opaque_occurrences > 0 && total_occurrences > 2 {
                return Err(TensorReductionError::AmbiguousMinkowskiIndex {
                    index,
                    occurrences: total_occurrences,
                });
            }
            match (
                occupants.integrated.as_slice(),
                occupants.outside.as_slice(),
            ) {
                ([], [vector]) => scalar_factors.push(indexed_vector(vector, &index)),
                ([], [left, right]) => scalar_factors.push(dot(left, right)),
                ([left, right], []) => scalar_factors.push(dot(left, right)),
                ([vector], []) => {
                    integrated.push(vector.clone());
                    outside.push(Outside::Index {
                        slot: occupants.slot.clone(),
                        index,
                    });
                }
                ([internal], [external]) => {
                    integrated.push(internal.clone());
                    outside.push(Outside::Vector(external.clone()));
                }
                (internal, external) => {
                    return Err(TensorReductionError::AmbiguousIndex {
                        index,
                        integrated: internal.len(),
                        outside: external.len(),
                    });
                }
            }
        }

        Ok(TensorMonomial {
            scalar: scalar_factors
                .into_iter()
                .fold(Atom::one(), |product, factor| product * factor),
            integrated,
            outside,
        })
    }

    fn is_integrated(&self, vector: &IndexedVector) -> bool {
        !self.external_vectors.contains(&vector.compact)
            && (self.integrated_heads.contains(vector.head.get_name())
                || self.integrated_vectors.contains(&vector.compact))
    }

    fn indexed_vector(
        &self,
        factor: AtomView<'_>,
    ) -> Result<Option<IndexedVector>, TensorReductionError> {
        let AtomView::Fun(function) = factor else {
            return Ok(None);
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let Some(AtomView::Fun(representation)) = arguments.last() else {
            return Ok(None);
        };
        if representation.get_symbol().get_namespace() != "spenso"
            || representation.get_symbol().get_stripped_name() != "mink"
            || representation.get_nargs() != 2
        {
            return Ok(None);
        }
        if arguments[..arguments.len() - 1]
            .iter()
            .any(|argument| contains_representation(argument))
        {
            return Ok(None);
        }
        let representation_arguments = representation.iter().collect::<Vec<_>>();
        let dimension = representation_arguments[0]
            .get_symbol()
            .and_then(idenso::CookSettings::dimension_payload)
            .unwrap_or_else(|| representation_arguments[0].to_owned());
        if dimension != self.dimension {
            return Err(TensorReductionError::DimensionMismatch {
                expected: self.dimension.clone(),
                found: dimension,
            });
        }

        let compact_representation = FunctionBuilder::new(representation.get_symbol())
            .add_arg(representation_arguments[0])
            .finish();
        let compact = arguments[..arguments.len() - 1]
            .iter()
            .fold(
                FunctionBuilder::new(function.get_symbol()),
                |builder, argument| builder.add_arg(*argument),
            )
            .add_arg(compact_representation)
            .finish();
        Ok(Some(IndexedVector {
            head: function.get_symbol(),
            compact,
            index: representation_arguments[1].to_owned(),
            slot: representation.as_view().to_owned(),
        }))
    }

    fn reduce_monomial(
        &self,
        monomial: TensorMonomial,
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        resources::checkpoint()?;
        let rank = monomial.integrated.len();
        if rank == 0 {
            return Ok(vec![TensorReductionTerm {
                coefficient: Atom::one(),
                tensor: SymbolicTensor::infer(monomial.scalar)?,
                integrated_orbit: None,
                projector_orbit: None,
                budget: ResourceBudget::capture_current(),
            }]);
        }
        let pair_count = rank / 2;
        let all_integrated_equal = monomial
            .integrated
            .windows(2)
            .all(|pair| pair[0] == pair[1]);
        let projector_vectors = monomial
            .outside
            .iter()
            .map(|outside| match outside {
                Outside::Vector(vector) => Some(vector.clone()),
                Outside::Index { .. } => None,
            })
            .collect::<Option<Vec<_>>>();

        if all_integrated_equal {
            return self.reduce_isotropic_internal(monomial, projector_vectors);
        }
        if let Some(projector_vectors) = &projector_vectors
            && projector_vectors.windows(2).all(|pair| pair[0] == pair[1])
        {
            return self.reduce_isotropic_projector(monomial, projector_vectors);
        }
        if let Some(projector_vectors) = projector_vectors {
            return self.reduce_orbit_projector(monomial, projector_vectors);
        }

        // Residual free indices are not repeated-vector color classes, so the
        // orbit projector cannot compact them.  Keep this explicitly bounded
        // labeled path only for tensor-valued output.
        if let Some(limit) = self.pairing_limit {
            let pairing_count = perfect_matching_count(rank)?;
            if pairing_count > limit {
                return Err(TensorReductionError::PairingLimit {
                    rank,
                    pairings: pairing_count,
                    limit,
                });
            }
        }
        if let Some(limit) = self.pairing_product_limit {
            let pairing_count = perfect_matching_count(rank)?;
            let pairing_products = pairing_count
                .checked_mul(pairing_count)
                .ok_or(TensorReductionError::CounterOverflow)?;
            if pairing_products > limit {
                return Err(TensorReductionError::PairingProductLimit {
                    rank,
                    products: pairing_products,
                    limit,
                });
            }
        }

        let pairings = perfect_matchings(rank)?;
        let internal_classes = classified_pairings(&monomial.integrated, &pairings)?;
        let outside_classes = classified_outside(&monomial.outside, &pairings)?;
        let output_classes = internal_classes
            .orbits
            .len()
            .checked_mul(outside_classes.atoms.len())
            .ok_or(TensorReductionError::CounterOverflow)?;
        if let Some(limit) = self.output_term_limit
            && output_classes > limit
        {
            return Err(TensorReductionError::OutputLimit {
                terms: output_classes,
                limit,
            });
        }

        let table = OrthogonalWeingarten::new(pair_count, self.dimension.clone())?;
        resources::checkpoint()?;
        let mut coefficient_by_key = BTreeMap::new();
        let mut coset_index = BTreeMap::new();
        for (index, (coset, coefficient)) in table.coefficients().enumerate() {
            coset_index.insert(coset_key_from_parts(coset.parts())?, index);
            coefficient_by_key.insert(index, coefficient.clone());
        }
        let coefficient_count = coefficient_by_key.len();
        let count_len = output_classes
            .checked_mul(coefficient_count)
            .ok_or(TensorReductionError::CounterOverflow)?;
        let bytes = count_len
            .checked_mul(mem::size_of::<u128>())
            .and_then(|bytes| u64::try_from(bytes).ok())
            .ok_or(TensorReductionError::CounterOverflow)?;
        resources::check_allocation(bytes)?;
        let mut counts = vec![0_u128; count_len];
        for (left_index, left) in pairings.iter().enumerate() {
            resources::checkpoint()?;
            let internal = internal_classes.class_of_pairing[left_index];
            for (right_index, right) in pairings.iter().enumerate() {
                resources::checkpoint()?;
                let outside = outside_classes.class_of_pairing[right_index];
                let key = matching_coset_key(left, right)?;
                let coefficient = coset_index[&key];
                let offset = (internal * outside_classes.atoms.len() + outside) * coefficient_count
                    + coefficient;
                counts[offset] = counts[offset]
                    .checked_add(1)
                    .ok_or(TensorReductionError::CounterOverflow)?;
            }
        }

        let mut terms = Vec::new();
        for (integrated_index, integrated_orbit) in internal_classes.orbits.into_iter().enumerate()
        {
            resources::checkpoint()?;
            for outside_index in 0..outside_classes.atoms.len() {
                let mut coefficient = Atom::Zero;
                for coefficient_index in 0..coefficient_count {
                    let offset = (integrated_index * outside_classes.atoms.len() + outside_index)
                        * coefficient_count
                        + coefficient_index;
                    let count = counts[offset];
                    if count != 0 {
                        coefficient += atom_from_count(count)?
                            * coefficient_by_key[&coefficient_index].clone();
                    }
                }
                if coefficient.is_zero() {
                    continue;
                }
                terms.push(TensorReductionTerm {
                    coefficient,
                    tensor: SymbolicTensor::infer(
                        monomial.scalar.clone()
                            * integrated_orbit.invariant()?
                            * outside_classes.atoms[outside_index].clone(),
                    )?,
                    integrated_orbit: Some(integrated_orbit.clone()),
                    projector_orbit: outside_classes.orbits[outside_index].clone(),
                    budget: ResourceBudget::capture_current(),
                });
            }
        }
        Ok(terms)
    }

    fn reduce_orbit_projector(
        &self,
        monomial: TensorMonomial,
        projector_vectors: Vec<Atom>,
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        OrthogonalWeingarten::validate_inverse_dimension(
            monomial.integrated.len() / 2,
            &self.dimension,
        )?;
        let (integrated_vectors, integrated_multiplicities, integrated_colors) =
            classified_colors(&monomial.integrated);
        let (projector_vectors, projector_multiplicities, projector_colors) =
            classified_colors(&projector_vectors);
        let projector = OrbitProjector::new(
            &integrated_colors,
            &projector_colors,
            self.dimension.clone(),
            self.output_term_limit,
        )?;

        projector
            .terms()
            .iter()
            .map(|term| {
                let integrated_orbit = ContractionOrbit {
                    vectors: integrated_vectors.clone().into_boxed_slice(),
                    multiplicities: integrated_multiplicities.clone().into_boxed_slice(),
                    alpha: term.internal_alpha().to_vec().into_boxed_slice(),
                    labeled_pairings: term.internal_labeled_pairings(),
                };
                let projector_orbit = ContractionOrbit {
                    vectors: projector_vectors.clone().into_boxed_slice(),
                    multiplicities: projector_multiplicities.clone().into_boxed_slice(),
                    alpha: term.projector_alpha().to_vec().into_boxed_slice(),
                    labeled_pairings: term.projector_labeled_pairings(),
                };
                Ok(TensorReductionTerm {
                    coefficient: term.coefficient().clone(),
                    tensor: SymbolicTensor::infer(
                        monomial.scalar.clone()
                            * integrated_orbit.invariant()?
                            * projector_orbit.invariant()?,
                    )?,
                    integrated_orbit: Some(integrated_orbit),
                    projector_orbit: Some(projector_orbit),
                    budget: ResourceBudget::capture_current(),
                })
            })
            .collect()
    }

    fn reduce_isotropic_internal(
        &self,
        monomial: TensorMonomial,
        projector_vectors: Option<Vec<Atom>>,
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        let pair_count = monomial.integrated.len() / 2;
        let integrated_orbit = repeated_orbit(monomial.integrated[0].clone(), pair_count)?;
        let coefficient = OrthogonalWeingarten::isotropic_pairing_coefficient(
            pair_count,
            self.dimension.clone(),
        )?;
        let internal = integrated_orbit.invariant()?;

        if let Some(projector_vectors) = projector_vectors {
            let orbits = contraction_orbits(&projector_vectors, self.output_term_limit)?;
            return orbits
                .into_iter()
                .map(|projector_orbit| {
                    Ok(TensorReductionTerm {
                        coefficient: coefficient.clone()
                            * atom_from_count(projector_orbit.labeled_pairings)?,
                        tensor: SymbolicTensor::infer(
                            monomial.scalar.clone()
                                * internal.clone()
                                * projector_orbit.invariant()?,
                        )?,
                        integrated_orbit: Some(integrated_orbit.clone()),
                        projector_orbit: Some(projector_orbit),
                        budget: ResourceBudget::capture_current(),
                    })
                })
                .collect();
        }

        let outside = outside_pairing_atoms(
            &monomial.outside,
            self.pairing_limit,
            self.output_term_limit,
        )?;
        outside
            .into_iter()
            .map(|(tensor, multiplicity)| {
                Ok(TensorReductionTerm {
                    coefficient: coefficient.clone() * atom_from_count(multiplicity)?,
                    tensor: SymbolicTensor::infer(
                        monomial.scalar.clone() * internal.clone() * tensor,
                    )?,
                    integrated_orbit: Some(integrated_orbit.clone()),
                    projector_orbit: None,
                    budget: ResourceBudget::capture_current(),
                })
            })
            .collect()
    }

    fn reduce_isotropic_projector(
        &self,
        monomial: TensorMonomial,
        projector_vectors: &[Atom],
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        let pair_count = monomial.integrated.len() / 2;
        let projector_orbit = repeated_orbit(projector_vectors[0].clone(), pair_count)?;
        let coefficient = OrthogonalWeingarten::isotropic_pairing_coefficient(
            pair_count,
            self.dimension.clone(),
        )?;
        let projector = projector_orbit.invariant()?;
        contraction_orbits(&monomial.integrated, self.output_term_limit)?
            .into_iter()
            .map(|integrated_orbit| {
                Ok(TensorReductionTerm {
                    coefficient: coefficient.clone()
                        * atom_from_count(integrated_orbit.labeled_pairings)?,
                    tensor: SymbolicTensor::infer(
                        monomial.scalar.clone() * integrated_orbit.invariant()? * projector.clone(),
                    )?,
                    integrated_orbit: Some(integrated_orbit),
                    projector_orbit: Some(projector_orbit.clone()),
                    budget: ResourceBudget::capture_current(),
                })
            })
            .collect()
    }
}

/// Tensor-reduction methods added to FeynKit diagrams without wrapping them.
///
/// # Examples
///
/// Reduce the numerator already stored on a generated diagram while keeping
/// FeynKit's native graph type.
///
/// ```no_run
/// use feynkit_graph::FeynmanDiagram;
/// use feynkit_tensor::{FeynmanDiagramTensorExt, TensorReducer, TensorReductionError};
/// use symbolica::{atom::Atom, parse};
///
/// fn scalar_numerator(diagram: &FeynmanDiagram) -> Result<Atom, TensorReductionError> {
///     diagram
///         .reduce_tensor_numerator(&TensorReducer::feynkit(parse!("4")))
///         .and_then(|reduction| reduction.into_expression())
///         .map(idenso::tensor::SymbolicTensor::into_expression)
/// }
///
/// fn scalar_graphs(
///     diagram: &FeynmanDiagram,
/// ) -> Result<Vec<FeynmanDiagram>, TensorReductionError> {
///     diagram.reduce_tensor_graphs(&TensorReducer::feynkit(parse!("4")))
/// }
/// ```
pub trait FeynmanDiagramTensorExt {
    /// Reduce the numerator and projector, selecting exactly the graph's internal edge momenta.
    /// Four-dimensional Lorentz slots are promoted to `dimension` before reduction.
    /// External momenta remain projector vectors; the scalar numerator prefactor is separate.
    /// As with `TensorReducer`, reduction requires at least one integrated vector.
    fn tensor_reduce(&self, dimension: Atom) -> Result<TensorReduction, TensorReductionError>;

    /// Reduce a selected region using an explicitly supplied projector and its
    /// paired internal momenta; boundary and external-state momenta stay fixed.
    /// A supplied expression replaces the selected numerator and projector,
    /// allowing reduction after tensor contractions or a UV expansion.
    fn tensor_reduce_of(
        &self,
        subgraph: &SuBitGraph,
        dimension: Atom,
        projector: &SymbolicTensor<PartialStructure>,
        expression: Option<&SymbolicTensor<PartialStructure>>,
    ) -> Result<TensorReduction, TensorReductionError>;

    /// Reduce the finalized diagram numerator and external-state projector.
    ///
    /// The request-wide numerator prefactor is intentionally kept separate
    /// from the returned tensor reduction.
    fn reduce_tensor_numerator(
        &self,
        reducer: &TensorReducer,
    ) -> Result<TensorReduction, TensorReductionError>;

    /// Split a fully contracted tensor numerator into scalar contributions.
    ///
    /// Every nonzero compact reduction term becomes one cloned diagram whose
    /// numerator contains both the exact tensor-projector coefficient and its
    /// scalar invariant. The external-state projector is consumed and reset to
    /// one, while the scalar numerator prefactor is retained. Topology IDs are
    /// deliberately preserved. Derived names are deterministic and distinct:
    /// `name.tensor[0]`, `name.tensor[1]`, and so on. Interaction and particle
    /// assignments remain available as provenance; the original numerator
    /// fragments are replaced by the scalar contribution. A numerator with
    /// residual free indices is rejected; use
    /// [`Self::reduce_tensor_numerator`] to retain its metric-tensor output.
    fn reduce_tensor_graphs(
        &self,
        reducer: &TensorReducer,
    ) -> Result<Vec<FeynmanDiagram>, TensorReductionError>;
}

impl FeynmanDiagramTensorExt for FeynmanDiagram {
    fn tensor_reduce(&self, dimension: Atom) -> Result<TensorReduction, TensorReductionError> {
        self.tensor_reduce_of(
            &self.underlying().full_filter(),
            dimension,
            self.projector(),
            None,
        )
    }

    fn tensor_reduce_of(
        &self,
        subgraph: &SuBitGraph,
        dimension: Atom,
        projector: &SymbolicTensor<PartialStructure>,
        expression: Option<&SymbolicTensor<PartialStructure>>,
    ) -> Result<TensorReduction, TensorReductionError> {
        let representation = FunctionBuilder::new(Minkowski::selfless_symbol())
            .add_arg(&dimension)
            .finish();
        let mut reducer = TensorReducer::new(dimension.clone());
        for (pair, edge, data) in self.underlying().iter_edges_of(subgraph) {
            if pair.is_paired() && data.data.external.is_none() && !data.data.is_dummy {
                let vector = FunctionBuilder::new(feynkit_graph::momentum_symbol())
                    .add_arg(edge.0)
                    .add_arg(&representation)
                    .finish();
                reducer = reducer.with_integrated_vector(vector)?;
            }
        }
        let projected = if let Some(expression) = expression {
            expression.clone()
        } else {
            let numerator = if *subgraph == self.underlying().full_filter() {
                self.numerator().clone()
            } else {
                self.numerator_of(subgraph, &SuBitGraph::empty(self.underlying().n_hedges()))?
            };
            numerator.multiply(projector)?
        };
        let projected = projected.map_representation_dimensions(|representation, current| {
            if representation == LibraryRep::from(Minkowski {}) && current == Atom::num(4).as_view()
            {
                dimension.clone()
            } else {
                current.to_owned()
            }
        })?;
        reducer.reduce(&projected)
    }

    fn reduce_tensor_numerator(
        &self,
        reducer: &TensorReducer,
    ) -> Result<TensorReduction, TensorReductionError> {
        ResourceBudget::new(reducer.resource_limits()).run(|| {
            resources::checkpoint()?;
            let projected_numerator = self.numerator().multiply(self.projector())?;
            resources::checkpoint()?;
            reducer.reduce(&projected_numerator)
        })
    }

    fn reduce_tensor_graphs(
        &self,
        reducer: &TensorReducer,
    ) -> Result<Vec<FeynmanDiagram>, TensorReductionError> {
        ResourceBudget::new(reducer.resource_limits()).run(|| {
            resources::checkpoint()?;
            let reduction = self.reduce_tensor_numerator(reducer)?;
            if !reduction.is_fully_contracted() {
                return Err(TensorReductionError::ResidualFreeIndices);
            }
            let graphs = reduction
                .terms
                .into_iter()
                .map(|term| term.expression())
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .filter(|numerator| !numerator.is_zero())
                .enumerate()
                .map(|(index, numerator)| {
                    resources::checkpoint()?;
                    let diagram = self
                        .clone()
                        .with_name(format!("{}.tensor[{index}]", self.name()))
                        .with_projector(SymbolicTensor::infer(Atom::one())?)
                        .with_numerator(numerator)
                        .map_err(TensorReductionError::Diagram)?;
                    resources::checkpoint()?;
                    Ok(diagram)
                })
                .collect::<Result<Vec<_>, TensorReductionError>>()?;
            resources::checkpoint()?;
            Ok(graphs)
        })
    }
}

/// Failures while parsing or materializing a tensor reduction.
#[derive(Debug, Error)]
pub enum TensorReductionError {
    /// A configured shared process-memory or elapsed-time allowance was exceeded.
    #[error("{0}")]
    Resource(#[from] ResourceLimitError),
    /// External basis entries must be compact Lorentz vectors.
    #[error("invalid compact external vector: {0}")]
    InvalidExternalVector(Atom),
    /// A nondegenerate external span is required for covariant projection.
    #[error(
        "cannot invert the external Gram matrix; choose an independent nondegenerate basis: {0}"
    )]
    ExternalGram(String),
    /// A reduction requires at least one loop/integrated-momentum selector.
    #[error("no integrated vectors selected; add a head or exact-vector selector")]
    NoIntegratedVectorsSelected,
    /// Scalar graph splitting cannot retain residual metric indices.
    #[error(
        "tensor reduction has residual free indices; use reduce_tensor_numerator for tensor-valued output"
    )]
    ResidualFreeIndices,
    /// A selected vector used a Lorentz dimension different from the reducer.
    #[error("Minkowski slot has dimension {found}, expected {expected}")]
    DimensionMismatch { expected: Atom, found: Atom },
    /// Tensor factors carrying Minkowski indices must have a non-negative
    /// integer exponent.
    #[error("tensor factor has unsupported exponent {0}")]
    InvalidVectorPower(Atom),
    /// More than one Einstein contraction was attached to an index.
    #[error("index {index} has {integrated} integrated and {outside} outside vector occurrences")]
    AmbiguousIndex {
        index: Atom,
        integrated: usize,
        outside: usize,
    },
    /// Generic/free-index expansion exceeded its labeled-pairing budget.
    #[error("rank-{rank} reduction needs {pairings} pairings, above limit {limit}")]
    PairingLimit {
        rank: usize,
        pairings: u128,
        limit: u128,
    },
    /// The generic projector would require too many relative matching pairs.
    #[error("rank-{rank} reduction needs {products} pairing products, above limit {limit}")]
    PairingProductLimit {
        rank: usize,
        products: u128,
        limit: u128,
    },
    /// Symmetry reduction still produced too many distinct invariant terms.
    #[error("tensor reduction produces at least {terms} terms, above limit {limit}")]
    OutputLimit { terms: usize, limit: usize },
    /// Checked tensor admission, contraction, or complete coefficient collection failed.
    #[error("tensor reduction algebra failed: {0}")]
    TensorAlgebra(String),
    /// A multiplicity did not fit an explicitly requested coefficient boundary.
    #[error("pairing multiplicity {0} does not fit the requested coefficient boundary")]
    MultiplicityOverflow(u128),
    /// A native matching or exponent representation cannot hold this rank.
    #[error("tensor rank {rank} exceeds matching representation capacity {maximum}")]
    UnsupportedRank { rank: usize, maximum: usize },
    /// The orthogonal-Weingarten coefficient engine failed.
    #[error(transparent)]
    Weingarten(WeingartenError),
    /// A scalar contribution could not be represented by a valid native graph.
    #[error(transparent)]
    Diagram(#[from] feynkit_graph::DiagramError),
    /// More than one Einstein contraction reused the same Minkowski index.
    #[error("Minkowski index {index} occurs {occurrences} times in one tensor term")]
    AmbiguousMinkowskiIndex { index: Atom, occurrences: usize },
    /// A matching multiplicity, allocation size, or exponent counter overflowed.
    #[error("tensor matching arithmetic exceeds its native representation")]
    CounterOverflow,
    /// Internal and projector color lists must describe the same tensor rank.
    #[error("internal colors have rank {internal}, but projector colors have rank {projector}")]
    ProjectorRankMismatch { internal: usize, projector: usize },
}

impl From<WeingartenError> for TensorReductionError {
    fn from(error: WeingartenError) -> Self {
        match error {
            WeingartenError::Resource(error) => Self::Resource(error),
            error => Self::Weingarten(error),
        }
    }
}

impl From<idenso::tensor::inference::TensorInferenceError> for TensorReductionError {
    fn from(error: idenso::tensor::inference::TensorInferenceError) -> Self {
        Self::TensorAlgebra(error.to_string())
    }
}

impl From<idenso::tensor::composition::TensorCompositionError> for TensorReductionError {
    fn from(error: idenso::tensor::composition::TensorCompositionError) -> Self {
        Self::TensorAlgebra(error.to_string())
    }
}

impl From<OrbitProjectorError> for TensorReductionError {
    fn from(error: OrbitProjectorError) -> Self {
        match error {
            OrbitProjectorError::Resource(error) => Self::Resource(error),
            OrbitProjectorError::ColorLengthMismatch {
                internal,
                projector,
            } => Self::ProjectorRankMismatch {
                internal,
                projector,
            },
            OrbitProjectorError::OddRank(rank) => Self::Weingarten(WeingartenError::OddRank(rank)),
            OrbitProjectorError::OutputLimit { terms, limit } => Self::OutputLimit { terms, limit },
            OrbitProjectorError::CounterOverflow => Self::CounterOverflow,
            OrbitProjectorError::Solve { rank, source } => {
                Self::Weingarten(WeingartenError::Solve { rank, source })
            }
        }
    }
}

#[derive(Clone, Debug)]
struct IndexedVector {
    head: Symbol,
    compact: Atom,
    index: Atom,
    slot: Atom,
}

struct IndexOccupants {
    slot: Atom,
    integrated: Vec<Atom>,
    outside: Vec<Atom>,
}

struct TensorMonomial {
    scalar: Atom,
    integrated: Vec<Atom>,
    outside: Vec<Outside>,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum Outside {
    Vector(Atom),
    Index { slot: Atom, index: Atom },
}

struct ClassifiedPairings {
    orbits: Vec<ContractionOrbit>,
    class_of_pairing: Vec<usize>,
}

struct ClassifiedOutside {
    atoms: Vec<Atom>,
    orbits: Vec<Option<ContractionOrbit>>,
    class_of_pairing: Vec<usize>,
}

fn contains_representation(argument: &AtomView<'_>) -> bool {
    matches!(argument, AtomView::Fun(function) if function.get_symbol().has_tag(&SPENSO_TAG.representation))
}

fn count_minkowski_index(
    expression: AtomView<'_>,
    dimension: &Atom,
    index: &Atom,
    multiplicity: usize,
) -> Result<usize, TensorReductionError> {
    resources::checkpoint()?;
    if let AtomView::Fun(function) = expression
        && function.get_symbol().get_namespace() == "spenso"
        && function.get_symbol().get_stripped_name() == "mink"
    {
        let arguments = function.iter().collect::<Vec<_>>();
        if let [slot_dimension, slot_index] = arguments.as_slice()
            && *slot_index == index.as_view()
        {
            let source_dimension = slot_dimension
                .get_symbol()
                .and_then(idenso::CookSettings::dimension_payload)
                .unwrap_or_else(|| slot_dimension.to_owned());
            if source_dimension == *dimension {
                return Ok(multiplicity);
            }
        }
    }
    if let AtomView::Pow(power) = expression {
        let base_occurrences =
            count_minkowski_index(power.get_base(), dimension, index, multiplicity)?;
        if base_occurrences == 0 {
            return Ok(0);
        }
        let exponent = i64::try_from(power.get_exp())
            .map_err(|_| TensorReductionError::InvalidVectorPower(power.get_exp().to_owned()))?;
        let exponent = usize::try_from(exponent)
            .map_err(|_| TensorReductionError::InvalidVectorPower(power.get_exp().to_owned()))?;
        return Ok(base_occurrences.saturating_mul(exponent));
    }
    expression
        .children()
        .try_fold(0_usize, |occurrences, child| {
            Ok(occurrences.saturating_add(count_minkowski_index(
                child,
                dimension,
                index,
                multiplicity,
            )?))
        })
}

#[cfg(test)]
fn contains_metric_pair(
    expression: AtomView<'_>,
    dimension: &Atom,
    left: &Atom,
    right: &Atom,
) -> bool {
    if let AtomView::Fun(function) = expression
        && function.get_symbol() == ETS.metric
    {
        let arguments = function.iter().collect::<Vec<_>>();
        let left_slot = minkowski_slot(dimension, left);
        let right_slot = minkowski_slot(dimension, right);
        if arguments.as_slice() == [left_slot.as_view(), right_slot.as_view()]
            || arguments.as_slice() == [right_slot.as_view(), left_slot.as_view()]
        {
            return true;
        }
    }
    expression
        .children()
        .any(|child| contains_metric_pair(child, dimension, left, right))
}

fn indexed_vector(compact: &Atom, index: &Atom) -> Atom {
    let AtomView::Fun(function) = compact.as_view() else {
        return compact.clone();
    };
    let arguments = function.iter().collect::<Vec<_>>();
    let Some(AtomView::Fun(representation)) = arguments.last() else {
        return compact.clone();
    };
    let mut indexed_representation = FunctionBuilder::new(representation.get_symbol());
    for argument in representation.iter() {
        indexed_representation = indexed_representation.add_arg(argument);
    }
    indexed_representation = indexed_representation.add_arg(index);
    arguments[..arguments.len() - 1]
        .iter()
        .fold(
            FunctionBuilder::new(function.get_symbol()),
            |builder, argument| builder.add_arg(*argument),
        )
        .add_arg(indexed_representation.finish())
        .finish()
}

fn dot(left: &Atom, right: &Atom) -> Atom {
    FunctionBuilder::new(SPENSO_TAG.dot)
        .add_arg(left)
        .add_arg(right)
        .finish()
}

fn outside_pair(left: &Outside, right: &Outside) -> Atom {
    match (left, right) {
        (Outside::Vector(left), Outside::Vector(right)) => dot(left, right),
        (Outside::Index { slot: left, .. }, Outside::Index { slot: right, .. }) => {
            FunctionBuilder::new(ETS.metric)
                .add_arg(left)
                .add_arg(right)
                .finish()
        }
        (Outside::Vector(vector), Outside::Index { index, .. })
        | (Outside::Index { index, .. }, Outside::Vector(vector)) => indexed_vector(vector, index),
    }
}

#[cfg(test)]
fn minkowski_slot(dimension: &Atom, index: &Atom) -> Atom {
    FunctionBuilder::new(Minkowski::selfless_symbol())
        .add_arg(dimension)
        .add_arg(index)
        .finish()
}

fn repeated_orbit(
    vector: Atom,
    pair_count: usize,
) -> Result<ContractionOrbit, TensorReductionError> {
    let rank = pair_count
        .checked_mul(2)
        .ok_or(TensorReductionError::CounterOverflow)?;
    let exponent = u16::try_from(pair_count).map_err(|_| TensorReductionError::CounterOverflow)?;
    Ok(ContractionOrbit {
        vectors: vec![vector].into_boxed_slice(),
        multiplicities: vec![rank].into_boxed_slice(),
        alpha: vec![exponent].into_boxed_slice(),
        labeled_pairings: perfect_matching_count(rank)?,
    })
}

fn classified_colors(colors: &[Atom]) -> (Vec<Atom>, Vec<usize>, Vec<usize>) {
    let mut counts = BTreeMap::<Atom, usize>::new();
    for color in colors {
        *counts.entry(color.clone()).or_default() += 1;
    }
    let vectors = counts.keys().cloned().collect::<Vec<_>>();
    let multiplicities = counts.values().copied().collect::<Vec<_>>();
    let positions = vectors
        .iter()
        .enumerate()
        .map(|(index, vector)| (vector, index))
        .collect::<BTreeMap<_, _>>();
    let classes = colors.iter().map(|color| positions[color]).collect();
    (vectors, multiplicities, classes)
}

fn contraction_orbits(
    colors: &[Atom],
    limit: Option<usize>,
) -> Result<Vec<ContractionOrbit>, TensorReductionError> {
    let (vectors, multiplicities, _) = classified_colors(colors);
    let mut cache = HashMap::new();
    let multiplicities_u16 = multiplicities
        .iter()
        .map(|count| u16::try_from(*count).map_err(|_| TensorReductionError::CounterOverflow))
        .collect::<Result<Vec<_>, _>>()?;
    let alpha = orbit_alpha(&multiplicities_u16, &mut cache, limit)?;
    Ok(alpha
        .into_iter()
        .map(|(alpha, labeled_pairings)| ContractionOrbit {
            vectors: vectors.clone().into_boxed_slice(),
            multiplicities: multiplicities.clone().into_boxed_slice(),
            alpha: alpha.into_boxed_slice(),
            labeled_pairings,
        })
        .collect())
}

fn orbit_alpha(
    counts: &[u16],
    cache: &mut HashMap<Vec<u16>, BTreeMap<Vec<u16>, u128>>,
    limit: Option<usize>,
) -> Result<BTreeMap<Vec<u16>, u128>, TensorReductionError> {
    resources::checkpoint()?;
    if let Some(cached) = cache.get(counts) {
        return Ok(cached.clone());
    }
    let alpha_len = counts
        .len()
        .checked_add(1)
        .and_then(|next| counts.len().checked_mul(next))
        .map(|size| size / 2)
        .ok_or(TensorReductionError::CounterOverflow)?;
    let Some(left) = counts.iter().position(|count| *count != 0) else {
        return Ok(BTreeMap::from([(vec![0; alpha_len], 1)]));
    };

    let mut after_left = counts.to_vec();
    after_left[left] -= 1;
    let mut result = BTreeMap::<Vec<u16>, u128>::new();
    for right in left..counts.len() {
        resources::checkpoint()?;
        let choices = after_left[right];
        if choices == 0 {
            continue;
        }
        let mut remaining = after_left.clone();
        remaining[right] -= 1;
        for (mut alpha, multiplicity) in orbit_alpha(&remaining, cache, limit)? {
            let exponent = &mut alpha[upper_triangle_index(left, right, counts.len())];
            *exponent = exponent
                .checked_add(1)
                .ok_or(TensorReductionError::CounterOverflow)?;
            let addition = multiplicity
                .checked_mul(u128::from(choices))
                .ok_or(TensorReductionError::CounterOverflow)?;
            let count = result.entry(alpha).or_default();
            *count = count
                .checked_add(addition)
                .ok_or(TensorReductionError::CounterOverflow)?;
            if let Some(limit) = limit
                && result.len() > limit
            {
                return Err(TensorReductionError::OutputLimit {
                    terms: result.len(),
                    limit,
                });
            }
        }
    }
    cache.insert(counts.to_vec(), result.clone());
    Ok(result)
}

fn upper_triangle_index(left: usize, right: usize, size: usize) -> usize {
    left * size - left * left.saturating_sub(1) / 2 + right - left
}

fn perfect_matching_count(rank: usize) -> Result<u128, TensorReductionError> {
    (1..rank).step_by(2).try_fold(1_u128, |count, factor| {
        resources::checkpoint()?;
        count
            .checked_mul(factor as u128)
            .ok_or(TensorReductionError::CounterOverflow)
    })
}

fn perfect_matchings(rank: usize) -> Result<Vec<Vec<usize>>, TensorReductionError> {
    fn generate(
        mask: u32,
        partners: &mut [usize],
        output: &mut Vec<Vec<usize>>,
    ) -> Result<(), TensorReductionError> {
        resources::checkpoint()?;
        if mask == 0 {
            output.push(partners.to_vec());
            return Ok(());
        }
        let left = mask.trailing_zeros() as usize;
        let without_left = mask & !(1 << left);
        let mut candidates = without_left;
        while candidates != 0 {
            let right = candidates.trailing_zeros() as usize;
            partners[left] = right;
            partners[right] = left;
            generate(without_left & !(1 << right), partners, output)?;
            candidates &= !(1 << right);
        }
        Ok(())
    }

    // The labeled enumerator uses one bit per vertex. Orbit/isotropic paths
    // do not use this representation and are not subject to this boundary.
    if rank > u32::BITS as usize {
        return Err(TensorReductionError::UnsupportedRank {
            rank,
            maximum: u32::BITS as usize,
        });
    }
    let count = usize::try_from(perfect_matching_count(rank)?)
        .map_err(|_| TensorReductionError::CounterOverflow)?;
    let bytes_per_matching = rank
        .checked_mul(mem::size_of::<usize>())
        .and_then(|bytes| bytes.checked_add(mem::size_of::<Vec<usize>>()))
        .ok_or(TensorReductionError::CounterOverflow)?;
    let bytes = count
        .checked_mul(bytes_per_matching)
        .and_then(|bytes| u64::try_from(bytes).ok())
        .ok_or(TensorReductionError::CounterOverflow)?;
    resources::check_allocation(bytes)?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(count)
        .map_err(|error| TensorReductionError::TensorAlgebra(error.to_string()))?;
    let mask = if rank == u32::BITS as usize {
        u32::MAX
    } else {
        (1_u32 << rank) - 1
    };
    generate(mask, &mut vec![0; rank], &mut output)?;
    Ok(output)
}

fn classified_pairings(
    colors: &[Atom],
    pairings: &[Vec<usize>],
) -> Result<ClassifiedPairings, TensorReductionError> {
    let mut counts = BTreeMap::<Atom, usize>::new();
    for color in colors {
        *counts.entry(color.clone()).or_default() += 1;
    }
    let vectors = counts.keys().cloned().collect::<Vec<_>>();
    let multiplicities = counts.values().copied().collect::<Vec<_>>();
    let vector_index = vectors
        .iter()
        .enumerate()
        .map(|(index, vector)| (vector.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let mut class_by_alpha = BTreeMap::<Vec<u16>, usize>::new();
    let mut class_of_pairing = Vec::with_capacity(pairings.len());
    let mut class_counts = Vec::<u128>::new();
    for pairing in pairings {
        let alpha = pairing_alpha(colors, &vector_index, pairing);
        let next = class_by_alpha.len();
        let class = *class_by_alpha.entry(alpha).or_insert(next);
        if class == class_counts.len() {
            class_counts.push(0);
        }
        class_counts[class] += 1;
        class_of_pairing.push(class);
    }
    let mut orbits = class_by_alpha
        .into_iter()
        .map(|(alpha, class)| {
            (
                class,
                ContractionOrbit {
                    vectors: vectors.clone().into_boxed_slice(),
                    multiplicities: multiplicities.clone().into_boxed_slice(),
                    alpha: alpha.into_boxed_slice(),
                    labeled_pairings: class_counts[class],
                },
            )
        })
        .collect::<Vec<_>>();
    orbits.sort_unstable_by_key(|(class, _)| *class);
    Ok(ClassifiedPairings {
        orbits: orbits.into_iter().map(|(_, orbit)| orbit).collect(),
        class_of_pairing,
    })
}

fn pairing_alpha(
    colors: &[Atom],
    color_index: &BTreeMap<Atom, usize>,
    pairing: &[usize],
) -> Vec<u16> {
    let mut alpha = vec![0; color_index.len() * (color_index.len() + 1) / 2];
    for left in 0..pairing.len() {
        let right = pairing[left];
        if left >= right {
            continue;
        }
        let left = color_index[&colors[left]];
        let right = color_index[&colors[right]];
        let (left, right) = if left <= right {
            (left, right)
        } else {
            (right, left)
        };
        alpha[upper_triangle_index(left, right, color_index.len())] += 1;
    }
    alpha
}

fn classified_outside(
    outside: &[Outside],
    pairings: &[Vec<usize>],
) -> Result<ClassifiedOutside, TensorReductionError> {
    if let Some(vectors) = outside
        .iter()
        .map(|outside| match outside {
            Outside::Vector(vector) => Some(vector.clone()),
            Outside::Index { .. } => None,
        })
        .collect::<Option<Vec<_>>>()
    {
        let classified = classified_pairings(&vectors, pairings)?;
        let atoms = classified
            .orbits
            .iter()
            .map(ContractionOrbit::invariant)
            .collect::<Result<Vec<_>, _>>()?;
        return Ok(ClassifiedOutside {
            atoms,
            orbits: classified.orbits.into_iter().map(Some).collect(),
            class_of_pairing: classified.class_of_pairing,
        });
    }

    let mut class_by_atom = BTreeMap::<Atom, usize>::new();
    let mut class_of_pairing = Vec::with_capacity(pairings.len());
    for pairing in pairings {
        let atom = outside_pairing_atom(outside, pairing);
        let next = class_by_atom.len();
        class_of_pairing.push(*class_by_atom.entry(atom).or_insert(next));
    }
    let mut atoms = vec![Atom::Zero; class_by_atom.len()];
    for (atom, class) in class_by_atom {
        atoms[class] = atom;
    }
    Ok(ClassifiedOutside {
        orbits: vec![None; atoms.len()],
        atoms,
        class_of_pairing,
    })
}

fn outside_pairing_atoms(
    outside: &[Outside],
    pairing_limit: Option<u128>,
    output_limit: Option<usize>,
) -> Result<BTreeMap<Atom, u128>, TensorReductionError> {
    let pairings = perfect_matching_count(outside.len())?;
    if let Some(limit) = pairing_limit
        && pairings > limit
    {
        return Err(TensorReductionError::PairingLimit {
            rank: outside.len(),
            pairings,
            limit,
        });
    }
    let mut atoms = BTreeMap::<Atom, u128>::new();
    for pairing in perfect_matchings(outside.len())? {
        resources::checkpoint()?;
        let count = atoms
            .entry(outside_pairing_atom(outside, &pairing))
            .or_default();
        *count = count
            .checked_add(1)
            .ok_or(TensorReductionError::CounterOverflow)?;
        if let Some(limit) = output_limit
            && atoms.len() > limit
        {
            return Err(TensorReductionError::OutputLimit {
                terms: atoms.len(),
                limit,
            });
        }
    }
    Ok(atoms)
}

fn outside_pairing_atom(outside: &[Outside], pairing: &[usize]) -> Atom {
    (0..pairing.len())
        .filter(|left| *left < pairing[*left])
        .fold(Atom::one(), |product, left| {
            product * outside_pair(&outside[left], &outside[pairing[left]])
        })
}

fn coset_key_from_parts(parts: &[usize]) -> Result<u64, TensorReductionError> {
    let initial = u64::try_from(parts.len()).map_err(|_| TensorReductionError::CounterOverflow)?;
    parts.iter().try_fold(initial, |key, part| {
        let part = u64::try_from(*part).map_err(|_| TensorReductionError::CounterOverflow)?;
        if part >= 32 {
            return Err(TensorReductionError::CounterOverflow);
        }
        key.checked_mul(32)
            .and_then(|key| key.checked_add(part))
            .ok_or(TensorReductionError::CounterOverflow)
    })
}

fn matching_coset_key(left: &[usize], right: &[usize]) -> Result<u64, TensorReductionError> {
    resources::checkpoint()?;
    let mut seen = 0_u32;
    if left.len() > u32::BITS as usize {
        return Err(TensorReductionError::UnsupportedRank {
            rank: left.len(),
            maximum: u32::BITS as usize,
        });
    }
    let all = if left.len() == u32::BITS as usize {
        u32::MAX
    } else {
        (1_u32 << left.len()) - 1
    };
    let mut part_multiplicities = [0_u8; u32::BITS as usize / 2 + 1];
    let mut part_count = 0_u64;
    while seen != all {
        let start = (!seen & all).trailing_zeros() as usize;
        let mut current = start;
        let mut size = 0;
        loop {
            seen |= 1 << current;
            let partner = left[current];
            seen |= 1 << partner;
            current = right[partner];
            size += 1;
            if current == start {
                break;
            }
        }
        part_multiplicities[size] = part_multiplicities[size]
            .checked_add(1)
            .ok_or(TensorReductionError::CounterOverflow)?;
        part_count += 1;
    }
    let mut key = part_count;
    for size in (1..part_multiplicities.len()).rev() {
        for _ in 0..part_multiplicities[size] {
            key = key
                .checked_mul(32)
                .and_then(|key| key.checked_add(size as u64))
                .ok_or(TensorReductionError::CounterOverflow)?;
        }
    }
    Ok(key)
}

fn atom_from_count(count: u128) -> Result<Atom, TensorReductionError> {
    resources::checkpoint()?;
    let count = Atom::num(count);
    resources::checkpoint()?;
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CosetType;
    use feynkit_graph::{DiagramEdge, DiagramVertex, NumeratorInterface};
    use feynkit_model::Model;
    use idenso::shorthands::UndoShorthands;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::{broadcast_symbol, tensor_symbol, vector_symbol};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::time::Duration;
    use symbolica::symbol;

    fn vectors() -> (Symbol, Symbol, Symbol, Symbol) {
        (
            vector_symbol!("feynkit_tensor_test::k"),
            vector_symbol!("feynkit_tensor_test::q"),
            vector_symbol!("feynkit_tensor_test::p"),
            vector_symbol!("feynkit_tensor_test::r"),
        )
    }

    fn compact(head: Symbol, label: i64, dimension: &Atom) -> Atom {
        FunctionBuilder::new(head)
            .add_arg(label)
            .add_arg(
                FunctionBuilder::new(symbol!("spenso::mink"))
                    .add_arg(dimension)
                    .finish(),
            )
            .finish()
    }

    fn indexed(head: Symbol, label: i64, dimension: &Atom, index: &Atom) -> Atom {
        FunctionBuilder::new(head)
            .add_arg(label)
            .add_arg(
                FunctionBuilder::new(symbol!("spenso::mink"))
                    .add_arg(dimension)
                    .add_arg(index)
                    .finish(),
            )
            .finish()
    }

    fn metric(dimension: &Atom, left: &Atom, right: &Atom) -> Atom {
        FunctionBuilder::new(ETS.metric)
            .add_arg(minkowski_slot(dimension, left))
            .add_arg(minkowski_slot(dimension, right))
            .finish()
    }

    #[test]
    fn count_diagnostics_are_disabled_by_default_and_explicit_when_requested() {
        let (k, q, _, _) = vectors();
        let dimension = Atom::var(symbol!("tensor_diagnostic_dimension"));
        let mu = Atom::var(symbol!("tensor_diagnostic_mu"));
        let nu = Atom::var(symbol!("tensor_diagnostic_nu"));
        let source =
            SymbolicTensor::infer(indexed(k, 1, &dimension, &mu) * indexed(q, 1, &dimension, &nu))
                .unwrap();
        let reducer = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .with_integrated_head(q);
        assert_eq!(reducer.pairing_limit, None);
        assert_eq!(reducer.pairing_product_limit, None);
        assert_eq!(reducer.output_term_limit, None);
        assert_eq!(reducer.resource_limits(), ResourceLimits::default());
        assert_eq!(reducer.reduce(&source).unwrap().max_rank(), 2);
        assert!(matches!(
            reducer.clone().with_pairing_limit(0).reduce(&source),
            Err(TensorReductionError::PairingLimit { limit: 0, .. })
        ));
        assert!(matches!(
            reducer
                .clone()
                .with_pairing_product_limit(0)
                .reduce(&source),
            Err(TensorReductionError::PairingProductLimit { limit: 0, .. })
        ));
        assert!(matches!(
            reducer.with_output_term_limit(0).reduce(&source),
            Err(TensorReductionError::OutputLimit { limit: 0, .. })
        ));
    }

    #[test]
    fn matching_and_multiplicity_representation_overflows_are_errors() {
        assert!(matches!(
            perfect_matchings(34),
            Err(TensorReductionError::UnsupportedRank {
                rank: 34,
                maximum: 32
            })
        ));
        assert!(perfect_matching_count(56).is_ok());
        assert!(matches!(
            perfect_matching_count(58),
            Err(TensorReductionError::CounterOverflow)
        ));
        assert!(matches!(
            coset_key_from_parts(&[1; 13]),
            Err(TensorReductionError::CounterOverflow)
        ));
    }

    #[test]
    fn disabled_diagnostics_do_not_compute_pairing_product_counters() {
        let (k, q, _, _) = vectors();
        let dimension = Atom::var(symbol!("tensor_free_rank_dimension"));
        let integrated = (0..36)
            .map(|index| compact(if index % 2 == 0 { k } else { q }, 1, &dimension))
            .collect();
        let outside = (0..36)
            .map(|index| {
                let index = Atom::var(symbol!(&format!("tensor_free_rank_index_{index}")));
                Outside::Index {
                    slot: minkowski_slot(&dimension, &index),
                    index,
                }
            })
            .collect();
        let result = TensorReducer::new(dimension).reduce_monomial(TensorMonomial {
            scalar: Atom::one(),
            integrated,
            outside,
        });
        // 35!! squared overflows u128, but this path does not request a product
        // diagnostic. Its real limitation is the labeled iterator's u32 mask.
        assert!(matches!(
            result,
            Err(TensorReductionError::UnsupportedRank {
                rank: 36,
                maximum: 32
            })
        ));
    }

    #[test]
    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    fn elapsed_resource_failure_precedes_tensor_work() {
        let (k, _, _, _) = vectors();
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let reducer = TensorReducer::new(Atom::num(4))
            .with_integrated_head(k)
            .with_resources(ResourceLimits {
                max_elapsed: Some(Duration::ZERO),
                ..ResourceLimits::default()
            });
        assert!(matches!(
            reducer.reduce(&source),
            Err(TensorReductionError::Resource(
                ResourceLimitError::Elapsed { .. }
            ))
        ));
    }

    #[test]
    fn outer_memory_budget_is_preserved_until_result_materialization() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::num(4);
        let source = SymbolicTensor::infer(
            dot(&compact(k, 1, &dimension), &compact(p, 1, &dimension)).pow(2),
        )
        .unwrap();
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);
        let measured = Arc::new(AtomicUsize::new(1));
        let monitor = Arc::clone(&measured);
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(4096),
                ..ResourceLimits::default()
            },
            move || {
                u64::try_from(monitor.load(Ordering::Relaxed)).map_err(|error| error.to_string())
            },
        );
        let result = budget.run(|| reducer.reduce(&source)).unwrap();
        measured.store(4097, Ordering::Relaxed);
        assert!(matches!(
            result.terms()[0].expression(),
            Err(TensorReductionError::Resource(
                ResourceLimitError::Memory { .. }
            ))
        ));
        assert!(matches!(
            result.into_expression(),
            Err(TensorReductionError::Resource(
                ResourceLimitError::Memory { .. }
            ))
        ));
        assert!(matches!(
            budget.run(|| reducer.reduce(&source)),
            Err(TensorReductionError::Resource(
                ResourceLimitError::Memory { .. }
            ))
        ));
    }

    #[test]
    fn public_reduce_checks_after_block_work_and_preserves_domain_errors() {
        let (k, q, p, _) = vectors();
        let dimension = symbol!("tensor_late_block::D").to_atom();
        let mu = symbol!("tensor_late_block::mu").to_atom();
        let nu = symbol!("tensor_late_block::nu").to_atom();
        let source = SymbolicTensor::infer(
            (indexed(k, 1, &dimension, &mu) + indexed(q, 1, &dimension, &mu))
                * (indexed(k, 1, &dimension, &nu) + indexed(q, 1, &dimension, &nu))
                * indexed(p, 1, &dimension, &mu)
                * indexed(p, 1, &dimension, &nu),
        )
        .unwrap();
        assert!(TensorReducer::has_tensor_sum(source.expression().as_view()).unwrap());
        let reducer = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q);
        // Warm the coefficient cache before calibrating the same steady-state work.
        reducer.reduce_scoped(&source).unwrap();
        let work_checks = Arc::new(AtomicUsize::new(0));
        let monitor = Arc::clone(&work_checks);
        let measure = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1_000_000),
                ..ResourceLimits::default()
            },
            move || {
                monitor.fetch_add(1, Ordering::Relaxed);
                Ok(0)
            },
        );
        // Measure the core work separately so refusal occurs only after its
        // final contraction, at the public successful-result boundary.
        measure.run(|| reducer.reduce_scoped(&source)).unwrap();
        let checks = work_checks.load(Ordering::Relaxed);
        let calls = AtomicUsize::new(0);
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1_000_000),
                ..ResourceLimits::default()
            },
            move || {
                Ok(if calls.fetch_add(1, Ordering::Relaxed) >= checks {
                    2_000_000
                } else {
                    0
                })
            },
        );
        assert!(matches!(
            budget.run(|| reducer.reduce(&source)),
            Err(TensorReductionError::Resource(
                ResourceLimitError::Memory { .. }
            ))
        ));
        assert!(reducer.reduce(&source).unwrap().is_fully_contracted());

        let calls = Arc::new(AtomicUsize::new(0));
        let monitor = Arc::clone(&calls);
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1_000_000),
                ..ResourceLimits::default()
            },
            move || {
                Ok(if monitor.fetch_add(1, Ordering::Relaxed) == 0 {
                    0
                } else {
                    2_000_000
                })
            },
        );
        assert!(matches!(
            budget.run(|| TensorReducer::new(dimension).reduce(&source)),
            Err(TensorReductionError::NoIntegratedVectorsSelected)
        ));
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn missing_requested_memory_monitor_is_a_typed_failure() {
        let (k, _, _, _) = vectors();
        let source = SymbolicTensor::infer(Atom::one()).unwrap();
        let reducer = TensorReducer::new(Atom::num(4)).with_integrated_head(k);
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(4096),
                ..ResourceLimits::default()
            },
            || Err("test monitor unavailable".into()),
        );
        assert!(matches!(
            budget.run(|| reducer.reduce(&source)),
            Err(TensorReductionError::Resource(
                ResourceLimitError::MonitorUnavailable { .. }
            ))
        ));
    }

    #[test]
    fn source_admission_preserves_domain_errors_without_accepting_invalid_tensors() {
        let (k, q, p, _) = vectors();
        let dimension = Atom::num(4);
        let mu = Atom::var(symbol!("feynkit_tensor_test::admission_mu"));
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        let vector = indexed(k, 1, &dimension, &mu);
        let source = reducer.admit(vector.clone()).unwrap();
        assert_eq!(source.rank(), 1);
        assert_eq!(source.expression(), &vector);
        let odd = reducer.reduce(&source).unwrap().into_expression().unwrap();
        assert!(odd.is_zero());
        assert_eq!(odd.structure(), source.structure());

        let malformed = vector * indexed(p, 1, &dimension, &mu) * indexed(q, 1, &dimension, &mu);
        assert!(SymbolicTensor::infer(malformed.clone()).is_err());
        assert!(matches!(
            reducer.admit(malformed),
            Err(TensorReductionError::AmbiguousIndex {
                index,
                integrated: 1,
                outside: 2,
            }) if index == mu
        ));
        let invalid_rank_one = FunctionBuilder::new(
            SPENSO_TAG.rank_one_tensor_symbol("feynkit_tensor_test::invalid_admission_vector"),
        )
        .add_arg(1)
        .finish();
        assert!(
            matches!(
                reducer.admit(invalid_rank_one),
                Err(TensorReductionError::TensorAlgebra(message))
                    if message.contains("requires exactly one final structural port")
            ),
            "a successful monomial parse must never override failed native admission"
        );
    }

    #[test]
    fn encoded_dimensions_keep_native_slots_and_scalar_projector_coefficients() {
        let (k, q, _, _) = vectors();
        let epsilon = Atom::var(symbol!("feynkit_tensor_test::encoded_epsilon"));
        let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let mu = Atom::num(0);
        let nu = Atom::num(1);
        let cooking = idenso::CookSettings::indices()
            .with_mode(idenso::CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let expected = cooking
            .try_cook(
                (metric(&dimension, &mu, &nu)
                    * dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
                    / &dimension)
                    .as_view(),
            )
            .unwrap();
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        for second in [
            indexed(k, 1, &dimension, &nu),
            indexed(k, 1, &dimension, &nu) + indexed(q, 1, &dimension, &nu),
        ] {
            let source = reducer
                .admit(
                    cooking
                        .try_cook((indexed(k, 1, &dimension, &mu) * second).as_view())
                        .unwrap(),
                )
                .unwrap();
            let result = reducer.reduce(&source).unwrap().into_expression().unwrap();
            assert_eq!(
                result.structure().logical_slots(),
                source.structure().logical_slots()
            );
            assert_eq!(result.rank(), 2);
            assert!((result.expression() - &expected).together().is_zero());
            let zero = source.with_checked_expression(Atom::Zero).unwrap();
            let reduced_zero = reducer.reduce(&zero).unwrap().into_expression().unwrap();
            assert!(reduced_zero.is_zero());
            assert_eq!(reduced_zero.structure(), source.structure());
            assert!(matches!(
                TensorReducer::new(Atom::num(6) - Atom::num(2) * &epsilon)
                    .with_integrated_head(k)
                    .reduce(&source),
                Err(TensorReductionError::DimensionMismatch { .. })
            ));
        }
        let opaque = spenso::tensor_symbol!("feynkit_tensor_test::encoded_opaque");
        let opaque = FunctionBuilder::new(opaque)
            .add_arg(minkowski_slot(&dimension, &mu))
            .finish();
        let opaque = cooking.try_cook(opaque.as_view()).unwrap();
        assert_eq!(
            count_minkowski_index(opaque.as_view(), &dimension, &mu, 1).unwrap(),
            1
        );
    }

    #[test]
    fn reduction_keeps_declared_spectator_axes_and_zero_interfaces() {
        use idenso::representations::Bispinor;
        use spenso::structure::representation::RepName;
        let (k, _, _, _) = vectors();
        let dimension = Atom::num(4);
        let mu = Atom::var(symbol!("retained_reduction_mu"));
        let nu = Atom::var(symbol!("retained_reduction_nu"));
        let left = Atom::var(symbol!("retained_reduction_left"));
        let right = Atom::var(symbol!("retained_reduction_right"));
        let spectator = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("retained_reduction_spin"))
            .add_arg(Bispinor {}.new_rep(4).to_symbolic([left]))
            .add_arg(Bispinor {}.new_rep(4).to_symbolic([right]))
            .finish();
        let source = SymbolicTensor::infer(
            &spectator * indexed(k, 1, &dimension, &mu) * indexed(k, 1, &dimension, &nu),
        )
        .unwrap();
        let mut slots = source.structure().logical_slots();
        slots.reverse();
        let source = SymbolicTensor::from_declared_parts(
            source.expression().clone(),
            PartialStructure::from_logical_slots(slots),
        )
        .unwrap();
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        let reduced = reducer.reduce(&source).unwrap();
        let result = reduced.expression().unwrap();
        assert_eq!(result.structure(), source.structure());
        assert_eq!(result.rank(), 4);
        assert!(result.expression().contains(&spectator));
        for term in reduced.terms() {
            assert_eq!(term.tensor().structure(), source.structure());
        }
        let vector = compact(k, 1, &dimension);
        let expected =
            &spectator * dot(&vector, &vector) * metric(&dimension, &mu, &nu) / &dimension;
        assert!((result.expression() - expected).together().is_zero());

        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let result = reducer.reduce(&zero).unwrap().into_expression().unwrap();
        assert!(result.is_zero());
        assert_eq!(result.structure(), source.structure());
        let odd = SymbolicTensor::infer(&spectator * indexed(k, 1, &dimension, &mu)).unwrap();
        let reduced = reducer.reduce(&odd).unwrap();
        assert!(!reduced.is_fully_contracted());
        let result = reduced.into_expression().unwrap();
        assert!(result.is_zero());
        assert_eq!(result.structure(), odd.structure());
    }

    #[test]
    fn closed_single_vector_moments_above_rank_32_are_exact() {
        let (k, q, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_isotropic_extended"));
        let kc = compact(k, 1, &dimension);
        let qc = compact(q, 1, &dimension);
        let pc = compact(p, 1, &dimension);
        let reducer = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q);
        for rank in [32, 34, 36] {
            let input = dot(&kc, &pc).pow(rank);
            let actual = reducer
                .reduce(&SymbolicTensor::infer(input).unwrap())
                .unwrap();
            let multiplicity = (1..rank)
                .step_by(2)
                .fold(Atom::one(), |value, factor| value * Atom::num(factor));
            let expected = multiplicity
                * dot(&kc, &kc).pow(rank / 2)
                * dot(&pc, &pc).pow(rank / 2)
                * OrthogonalWeingarten::isotropic_pairing_coefficient(
                    rank as usize / 2,
                    dimension.clone(),
                )
                .unwrap();
            assert!(actual.is_fully_contracted());
            assert!(
                (actual.into_expression().unwrap().into_expression() - expected)
                    .together()
                    .is_zero(),
                "closed isotropic rank {rank} differs from the exact moment"
            );
        }
        // The same existing isotropic-projector path also handles mixed
        // integrated vectors above the former generic rank-20 admission cap.
        let mixed = dot(&kc, &pc).pow(2) * dot(&qc, &pc).pow(20);
        let mixed = reducer
            .reduce(&SymbolicTensor::infer(mixed).unwrap())
            .unwrap();
        assert!(mixed.is_fully_contracted());
        assert_eq!(mixed.max_rank(), 22);
        assert_eq!(mixed.terms().len(), 2);
        assert_eq!(
            mixed
                .terms()
                .iter()
                .map(|term| { term.integrated_orbit().unwrap().labeled_pairings() })
                .sum::<u128>(),
            perfect_matching_count(22).unwrap()
        );
    }

    #[test]
    fn weighted_compact_dots_use_the_isotropic_second_moment() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("weighted_dot_D"));
        let k_vector = compact(k, 1, &dimension);
        let p_vector = compact(p, 1, &dimension);
        let inner = dot(&k_vector, &p_vector);
        let coefficient = (Atom::var(symbol!("weighted_dot_coefficient")) + Atom::one()).pow(20);
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        // Rotational invariance gives <k_mu k_nu> = k^2 g_mu_nu / D.
        let expected =
            &coefficient * dot(&k_vector, &k_vector) * dot(&p_vector, &p_vector) / &dimension;
        let mu = Atom::var(symbol!("weighted_dot_mu"));
        let components = indexed(k, 1, &dimension, &mu) * indexed(p, 1, &dimension, &mu);
        for scalar in [inner, components] {
            let input = dot(&(&coefficient * scalar * &k_vector), &p_vector);
            let result = reducer
                .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
                .unwrap()
                .into_expression()
                .unwrap()
                .into_expression();
            assert_eq!(result, expected);
        }
    }

    fn diagram_with_numerator(name: &str, numerator: Atom) -> FeynmanDiagram {
        diagram_with_payload(name, numerator, Atom::one(), Atom::one())
    }

    fn diagram_with_payload(
        name: &str,
        numerator: Atom,
        numerator_prefactor: Atom,
        projector: Atom,
    ) -> FeynmanDiagram {
        let model = Model::from_json(include_str!(
            "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        ))
        .unwrap();
        let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
        let particle = model.particle_id("scalar_0").unwrap();
        let mut left = DiagramVertex::interaction("left", rule);
        let numerator = NumeratorInterface::admit(numerator, None).unwrap();
        left.numerator = numerator.clone();
        let mut builder = FeynmanDiagram::builder(model, name)
            .numerator(numerator)
            .numerator_prefactor(NumeratorInterface::admit(numerator_prefactor, None).unwrap())
            .projector(NumeratorInterface::admit(projector, None).unwrap());
        let left = builder.add_vertex(left);
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        for _ in 0..3 {
            builder
                .add_edge(left, right, DiagramEdge::new(particle, false))
                .unwrap();
        }
        let diagram = builder.build().unwrap();
        diagram.validate().unwrap();
        diagram
    }

    fn brute_projected_expression(
        integrated: &[Atom],
        projector: &[Atom],
        dimension: &Atom,
    ) -> Atom {
        let pairings = perfect_matchings(integrated.len()).unwrap();
        let table = OrthogonalWeingarten::for_rank(integrated.len(), dimension.clone()).unwrap();
        let mut result = Atom::Zero;
        for left in &pairings {
            let internal = (0..left.len())
                .filter(|index| *index < left[*index])
                .fold(Atom::one(), |product, index| {
                    product * dot(&integrated[index], &integrated[left[index]])
                });
            for right in &pairings {
                let outside = (0..right.len())
                    .filter(|index| *index < right[*index])
                    .fold(Atom::one(), |product, index| {
                        product * dot(&projector[index], &projector[right[index]])
                    });
                let coset = CosetType::between_matchings(left, right).unwrap();
                result += table.coefficient(&coset).unwrap() * &internal * outside;
            }
        }
        result
    }

    #[test]
    fn reduction_requires_an_integrated_vector_selector() {
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_no_selector"));
        let error = TensorReducer::new(dimension)
            .reduce(&SymbolicTensor::infer((Atom::one().as_view()).to_owned()).unwrap())
            .unwrap_err();
        assert!(matches!(
            error,
            TensorReductionError::NoIntegratedVectorsSelected
        ));
    }

    #[test]
    fn allocation_free_coset_keys_match_the_public_definition() {
        let pairings = perfect_matchings(6).unwrap();
        for left in &pairings {
            for right in &pairings {
                let coset = CosetType::between_matchings(left, right).unwrap();
                assert_eq!(
                    matching_coset_key(left, right).unwrap(),
                    coset_key_from_parts(coset.parts()).unwrap()
                );
            }
        }
    }

    #[test]
    fn exact_compound_dimension_selectors_preserve_unselected_vectors_and_open_ports() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::num(4)
            - Atom::num(2) * Atom::var(symbol!("feynkit_tensor_test::selector_epsilon"));
        let [mu, nu] = [Atom::num(7), Atom::num(11)];
        let k0 = compact(k, 0, &dimension);
        let k1 = compact(k, 1, &dimension);
        let unselected = indexed(k, 1, &dimension, &mu) * indexed(k, 1, &dimension, &nu);
        let source = NumeratorInterface::admit(
            indexed(k, 0, &dimension, &mu) * indexed(k, 0, &dimension, &nu) + &unselected,
            None,
        )
        .unwrap();
        let source = SymbolicTensor::from_declared_parts(
            source.expression().clone(),
            PartialStructure::from_logical_slots(
                source.structure().logical_slots().into_iter().rev(),
            ),
        )
        .unwrap();
        let reducer = TensorReducer::new(dimension.clone())
            .with_integrated_vector(k0.clone())
            .unwrap();
        let result = reducer.reduce(&source).unwrap().into_expression().unwrap();
        let expected = NumeratorInterface::admit(
            dot(&k0, &k0) * metric(&dimension, &mu, &nu) / &dimension + &unselected,
            None,
        )
        .unwrap();
        assert_eq!(
            result.structure().logical_slots(),
            source.structure().logical_slots()
        );
        assert!(
            (result.expression() - expected.expression())
                .together()
                .is_zero()
        );
        let head = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&source)
            .unwrap()
            .into_expression()
            .unwrap();
        let expected_head = NumeratorInterface::admit(
            (dot(&k0, &k0) + dot(&k1, &k1)) * metric(&dimension, &mu, &nu) / &dimension,
            None,
        )
        .unwrap();
        assert!(
            (head.expression() - expected_head.expression())
                .together()
                .is_zero()
        );
        assert_ne!(
            result.expression(),
            head.expression(),
            "the exact selector must leave k(1) outside"
        );
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let zero_result = reducer.reduce(&zero).unwrap().into_expression().unwrap();
        assert!(zero_result.is_zero());
        assert_eq!(zero_result.structure(), source.structure());

        let vector = NumeratorInterface::admit(indexed(k, 0, &dimension, &mu), None).unwrap();
        let external = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_external_vector(k1.clone())
            .unwrap();
        let actual = external.reduce(&vector).unwrap().into_expression().unwrap();
        let expected = NumeratorInterface::admit(
            dot(&k0, &k1) / dot(&k1, &k1) * indexed(k, 1, &dimension, &mu),
            None,
        )
        .unwrap();
        assert_eq!(actual.structure(), vector.structure());
        assert!(
            (actual.expression() - expected.expression())
                .together()
                .is_zero()
        );
    }

    #[test]
    fn graph_compound_dimension_preserves_selected_edge_and_declared_ports() {
        use linnet::half_edge::{involution::EdgeIndex, subgraph::ModifySubSet};
        let momentum = feynkit_graph::momentum_symbol();
        let [mu, nu] = [Atom::num(7), Atom::num(11)];
        let four = Atom::num(4);
        let numerator = indexed(momentum, 0, &four, &mu) * indexed(momentum, 0, &four, &nu)
            + indexed(momentum, 1, &four, &mu) * indexed(momentum, 1, &four, &nu);
        let diagram = diagram_with_numerator("compound-selected-edge", numerator);
        let source = SymbolicTensor::from_declared_parts(
            diagram.numerator().expression().clone(),
            PartialStructure::from_logical_slots(
                diagram
                    .numerator()
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .rev(),
            ),
        )
        .unwrap();
        let mut selected: SuBitGraph = diagram.underlying().empty_subgraph();
        selected.add(diagram.underlying()[&EdgeIndex(0)].1);
        let dimension = Atom::num(4)
            - Atom::num(2) * Atom::var(symbol!("feynkit_tensor_test::graph_selector_epsilon"));
        let result = diagram
            .tensor_reduce_of(
                &selected,
                dimension.clone(),
                diagram.projector(),
                Some(&source),
            )
            .unwrap()
            .into_expression()
            .unwrap();
        let k = compact(momentum, 0, &dimension);
        let expected = NumeratorInterface::admit(
            dot(&k, &k) * metric(&dimension, &mu, &nu) / &dimension
                + indexed(momentum, 1, &dimension, &mu) * indexed(momentum, 1, &dimension, &nu),
            None,
        )
        .unwrap();
        let expected_ports = expected
            .structure()
            .logical_slots()
            .into_iter()
            .rev()
            .collect::<Vec<_>>();
        assert_eq!(result.structure().logical_slots(), expected_ports);
        assert!(
            (result.expression() - expected.expression())
                .together()
                .is_zero()
        );
        let zero = source.with_checked_expression(Atom::Zero).unwrap();
        let zero = diagram
            .tensor_reduce_of(&selected, dimension, diagram.projector(), Some(&zero))
            .unwrap()
            .into_expression()
            .unwrap();
        assert!(zero.is_zero());
        assert_eq!(zero.structure().logical_slots(), expected_ports);
        assert_eq!(source.rank(), 2);
        assert!(
            source
                .structure()
                .logical_slots()
                .iter()
                .all(|slot| slot.dim().to_symbolic() == four)
        );
    }

    #[test]
    fn graph_tensor_reduce_selects_internal_edges_and_promotes_lorentz_slots() {
        let model = Model::from_json(include_str!(
            "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        ))
        .unwrap();
        let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
        let particle = model.particle_id("scalar_0").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "internal-selection");
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let mut incoming = DiagramEdge::new(particle, false);
        incoming.external = Some(feynkit_graph::ExternalLeg {
            name: "in".into(),
            index: 0,
            state: feynkit_graph::ExternalState::Incoming,
            connection: 0,
        });
        let internal = builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        let outside = builder.add_edge(None, left, incoming).unwrap();
        let mut outgoing = DiagramEdge::new(particle, false);
        outgoing.external = Some(feynkit_graph::ExternalLeg {
            name: "out".into(),
            index: 1,
            state: feynkit_graph::ExternalState::Outgoing,
            connection: 1,
        });
        builder.add_edge(right, None, outgoing).unwrap();
        let momentum = feynkit_graph::momentum_symbol();
        let d = Atom::var(symbol!("graph_reduce_D"));
        let mu = Atom::var(symbol!("graph_reduce_mu"));
        let nu = Atom::var(symbol!("graph_reduce_nu"));
        let four = Atom::num(4);
        let numerator = indexed(momentum, internal.0 as i64, &four, &mu)
            * indexed(momentum, internal.0 as i64, &four, &nu);
        let projector = indexed(momentum, outside.0 as i64, &four, &mu)
            * indexed(momentum, outside.0 as i64, &four, &nu);
        let diagram = builder
            .numerator(NumeratorInterface::admit(numerator, None).unwrap())
            .projector(NumeratorInterface::admit(projector, None).unwrap())
            .numerator_prefactor(SymbolicTensor::infer(Atom::num(7)).unwrap())
            .build()
            .unwrap();
        let result = diagram
            .tensor_reduce(d.clone())
            .unwrap()
            .into_expression()
            .unwrap()
            .into_expression();
        let k = compact(momentum, internal.0 as i64, &d);
        let p = compact(momentum, outside.0 as i64, &d);
        let expected = dot(&k, &k) * dot(&p, &p) / &d;
        assert_eq!(result, expected);
        assert_eq!(diagram.numerator_prefactor().expression(), &Atom::num(7));

        use linnet::half_edge::{
            involution::{EdgeIndex, HedgePair},
            subgraph::ModifySubSet,
        };
        let pair = diagram.underlying()[&EdgeIndex(internal.0)].1;
        let mut selected: SuBitGraph = diagram.underlying().empty_subgraph();
        selected.add(pair);
        let partial = diagram
            .tensor_reduce_of(
                &selected,
                d.clone(),
                diagram.projector(),
                Some(&diagram.numerator().multiply(diagram.projector()).unwrap()),
            )
            .unwrap()
            .into_expression()
            .unwrap()
            .into_expression();
        assert_eq!(partial, expected);
        let product = diagram.numerator().multiply(diagram.projector()).unwrap();
        let expanded = product
            .with_checked_expression(Atom::num(3) + Atom::num(2) * product.expression())
            .unwrap();
        let reduced_expansion = diagram
            .tensor_reduce_of(&selected, d.clone(), diagram.projector(), Some(&expanded))
            .unwrap()
            .into_expression()
            .unwrap()
            .into_expression();
        assert_eq!(reduced_expansion, Atom::num(3) + Atom::num(2) * &expected);
        let HedgePair::Paired { source, .. } = pair else {
            panic!("internal propagator")
        };
        let mut boundary: SuBitGraph = diagram.underlying().empty_subgraph();
        boundary.add(source);
        assert!(matches!(
            diagram.tensor_reduce_of(&boundary, d, diagram.projector(), None),
            Err(TensorReductionError::NoIntegratedVectorsSelected)
        ));
    }

    #[test]
    fn fully_contracted_reduction_splits_into_named_scalar_graphs() {
        let (k, q, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_graph_split"));
        let mut numerator = Atom::one();
        for position in 0..4 {
            let index = Atom::var(symbol!(format!(
                "feynkit_tensor_test::graph_split_mu{position}"
            )));
            let integrated = if position < 2 { k } else { q };
            numerator *= indexed(integrated, 1, &dimension, &index);
            numerator *= indexed(p, 1, &dimension, &index);
        }
        let diagram = diagram_with_numerator("vacuum", numerator);
        let reducer = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .with_integrated_head(q);
        let reduction = diagram.reduce_tensor_numerator(&reducer).unwrap();
        let expected_numerators = reduction
            .terms()
            .iter()
            .map(TensorReductionTerm::expression)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let graphs = diagram.reduce_tensor_graphs(&reducer).unwrap();

        assert_eq!(graphs.len(), 2);
        assert_eq!(
            graphs
                .iter()
                .map(|graph| graph.numerator().clone())
                .collect::<Vec<_>>(),
            expected_numerators
        );
        for (index, graph) in graphs.iter().enumerate() {
            assert_eq!(graph.id(), diagram.id());
            assert_eq!(graph.name(), format!("vacuum.tensor[{index}]"));
            graph.validate().unwrap();
            let json = graph.to_json().unwrap();
            let dot = graph.to_dot().unwrap();
            FeynmanDiagram::from_json(graph.model_arc(), &json).unwrap();
            FeynmanDiagram::from_dot(graph.model_arc(), &dot).unwrap();
        }
    }

    #[test]
    fn diagram_reduction_consumes_projector_and_retains_scalar_prefactor() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_projected_graph"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::projected_graph_mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::projected_graph_nu"));
        let prefactor = Atom::var(symbol!("feynkit_tensor_test::projected_graph_prefactor"));
        let numerator = indexed(k, 1, &dimension, &mu) * indexed(k, 1, &dimension, &nu);
        let projector = indexed(p, 1, &dimension, &mu) * indexed(p, 1, &dimension, &nu);
        let diagram = diagram_with_payload(
            "projected-vacuum",
            numerator,
            prefactor.clone(),
            projector.clone(),
        );
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(p, 1, &dimension))
            / dimension;

        let reduction = diagram.reduce_tensor_numerator(&reducer).unwrap();
        assert_eq!(reduction.expression().unwrap().into_expression(), expected);
        let graphs = diagram.reduce_tensor_graphs(&reducer).unwrap();
        assert_eq!(graphs.len(), 1);
        let graph = &graphs[0];
        assert_eq!(graph.projector().expression(), &Atom::one());
        assert_eq!(graph.numerator_prefactor().expression(), &prefactor);
        assert_eq!(graph.numerator().expression(), &expected);
        assert_eq!(diagram.projector().expression(), &projector);
        graph.validate().unwrap();
        FeynmanDiagram::from_json(graph.model_arc(), &graph.to_json().unwrap()).unwrap();
        FeynmanDiagram::from_dot(graph.model_arc(), &graph.to_dot().unwrap()).unwrap();
    }

    #[test]
    fn scalar_graph_split_rejects_residual_free_indices() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_graph_free"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::graph_free_mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::graph_free_nu"));
        let numerator = indexed(k, 1, &dimension, &mu) * indexed(k, 1, &dimension, &nu);
        let diagram = diagram_with_numerator("tensor-vacuum", numerator);
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);

        let tensor_result = diagram.reduce_tensor_numerator(&reducer).unwrap();
        assert!(!tensor_result.is_fully_contracted());
        assert!(matches!(
            diagram.reduce_tensor_graphs(&reducer),
            Err(TensorReductionError::ResidualFreeIndices)
        ));
    }

    #[test]
    fn rank_two_contracted_projection_uses_spenso_dots() {
        let (k, _, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::nu"));
        let input = indexed(k, 1, &dimension, &mu)
            * indexed(k, 2, &dimension, &nu)
            * indexed(p, 1, &dimension, &mu)
            * indexed(r, 1, &dimension, &nu);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 2, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(r, 1, &dimension))
            / dimension;
        assert_eq!(result.expression().unwrap().into_expression(), expected);
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn squared_dot_expansion_preserves_rank_two_angular_average() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_squared_dot"));
        let loop_vector = compact(k, 1, &dimension);
        let external_vector = compact(p, 1, &dimension);
        let input = dot(&loop_vector, &external_vector)
            .pow(Atom::num(2))
            .undo_dots::<AbstractIndex>()
            .unwrap();
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        // Independent contractions in (k.p)^2 give <k_mu k_nu> = k^2 g_mu_nu / D.
        let expected =
            dot(&loop_vector, &loop_vector) * dot(&external_vector, &external_vector) / dimension;
        assert_eq!(result.expression().unwrap().into_expression(), expected);
    }

    #[test]
    fn compact_spectator_dots_match_explicit_projection_through_rank_eight() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_compact"));
        let kc = compact(k, 1, &dimension);
        let pc = compact(p, 1, &dimension);
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        let mut even_expected = Atom::one();
        for rank in 1..=8 {
            let input = dot(&kc, &pc).pow(Atom::num(rank));
            let compact_result = reducer
                .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
                .unwrap();
            let explicit = input.undo_dots::<AbstractIndex>().unwrap();
            let explicit_result = reducer
                .reduce(&SymbolicTensor::infer((explicit.as_view()).to_owned()).unwrap())
                .unwrap();
            assert_eq!(
                compact_result.expression().unwrap().into_expression(),
                explicit_result.expression().unwrap().into_expression()
            );
            let expected = if rank % 2 == 0 {
                even_expected *= Atom::num(rank - 1) * dot(&kc, &kc) * dot(&pc, &pc)
                    / (&dimension + Atom::num(rank - 2));
                even_expected.clone()
            } else {
                Atom::Zero
            };
            assert!(
                (compact_result.expression().unwrap().into_expression() - expected)
                    .together()
                    .is_zero()
            );
            assert_eq!(compact_result.max_rank(), rank as usize);
        }
    }

    #[test]
    fn compact_dots_mix_with_free_indices_and_preserve_scalar_weights() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_compact_free"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::mu_compact_free"));
        let kc = compact(k, 1, &dimension);
        let pc = compact(p, 1, &dimension);
        let reducer = TensorReducer::new(dimension.clone()).with_integrated_head(k);
        let weight = (Atom::one() + dot(&kc, &kc)).pow(Atom::num(-2));
        let input = &weight * dot(&kc, &pc) * indexed(k, 1, &dimension, &mu);
        let expected = weight * dot(&kc, &kc) * indexed(p, 1, &dimension, &mu) / dimension;
        let result = reducer
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert!(
            (result.expression().unwrap().into_expression() - expected)
                .together()
                .is_zero()
        );
        assert!(!result.is_fully_contracted());
    }

    #[test]
    fn compact_dot_external_projection_is_idempotent() {
        let (k, q, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_compact_external"));
        let kc = compact(k, 1, &dimension);
        let qc = compact(q, 1, &dimension);
        let pc = compact(p, 1, &dimension);
        let rc = compact(r, 1, &dimension);
        let reducer = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .with_external_vector(pc.clone())
            .unwrap();
        let input = dot(&kc, &rc) * dot(&qc, &rc);
        let parallel = dot(&kc, &pc) * dot(&qc, &pc) * dot(&rc, &pc).pow(Atom::num(2))
            / dot(&pc, &pc).pow(Atom::num(2));
        let transverse = (dot(&kc, &qc) - dot(&kc, &pc) * dot(&qc, &pc) / dot(&pc, &pc))
            * (dot(&rc, &rc) - dot(&rc, &pc).pow(Atom::num(2)) / dot(&pc, &pc))
            / (dimension - Atom::one());
        let result = reducer
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap()
            .expression()
            .unwrap()
            .into_expression();
        assert!((&result - parallel - transverse).together().is_zero());
        let repeated = reducer
            .reduce(&SymbolicTensor::infer((result.as_view()).to_owned()).unwrap())
            .unwrap()
            .expression()
            .unwrap()
            .into_expression();
        assert!((result - repeated).together().is_zero());
        // Rational longitudinal weights are invariants of the transverse average.
        let invariant = dot(&kc, &pc).pow(Atom::num(-1)) * dot(&kc, &qc);
        assert_eq!(
            reducer
                .reduce(&SymbolicTensor::infer((invariant.as_view()).to_owned()).unwrap())
                .unwrap()
                .expression()
                .unwrap()
                .into_expression(),
            invariant
        );
    }

    #[test]
    fn compact_spectator_dots_enforce_power_and_dimension_validity() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_compact_limits"));
        let kc = compact(k, 1, &dimension);
        let pc = compact(p, 1, &dimension);
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);
        for exponent in [Atom::num(-1), Atom::num(1) / Atom::num(2)] {
            let input = dot(&kc, &pc).pow(exponent);
            assert!(matches!(
                reducer.reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap()),
                Err(TensorReductionError::InvalidVectorPower(_))
            ));
        }
        for rank in [21, 59] {
            let input = dot(&kc, &pc).pow(rank);
            assert!(
                reducer
                    .reduce(&SymbolicTensor::infer(input).unwrap())
                    .unwrap()
                    .into_expression()
                    .unwrap()
                    .is_zero()
            );
        }
        let wrong_dimension = dot(&kc, &compact(p, 1, &Atom::num(3)));
        assert!(matches!(
            reducer.parse_monomial(wrong_dimension.as_view()),
            Err(TensorReductionError::DimensionMismatch { .. })
        ));
    }

    #[test]
    fn rank_six_two_vector_symmetry_has_two_internal_orbits() {
        let (k, q, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D6"));
        let indices = (0..6)
            .map(|index| Atom::var(symbol!(format!("feynkit_tensor_test::mu{index}"))))
            .collect::<Vec<_>>();
        let mut input = Atom::one();
        for (position, index) in indices.iter().enumerate() {
            let head = if position < 2 { k } else { q };
            input *= indexed(head, 1, &dimension, index);
            input *= indexed(p, 1, &dimension, index);
        }
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert_eq!(result.terms().len(), 2);
        let multiplicities = result
            .terms()
            .iter()
            .map(|term| term.integrated_orbit().unwrap().labeled_pairings())
            .collect::<BTreeSet<_>>();
        assert_eq!(multiplicities, BTreeSet::from([3, 12]));
    }

    #[test]
    fn rank_four_mixed_projection_matches_the_inverse_gram_matrix() {
        let (k, q, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D4"));
        let indices = (0..4)
            .map(|index| Atom::var(symbol!(format!("feynkit_tensor_test::nu{index}"))))
            .collect::<Vec<_>>();
        let internal = [k, k, q, q];
        let external = [p, p, r, r];
        let mut input = Atom::one();
        for position in 0..4 {
            input *= indexed(internal[position], 1, &dimension, &indices[position]);
            input *= indexed(external[position], 1, &dimension, &indices[position]);
        }
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();

        let k = compact(k, 1, &dimension);
        let q = compact(q, 1, &dimension);
        let p = compact(p, 1, &dimension);
        let r = compact(r, 1, &dimension);
        let denominator = dimension.clone()
            * (dimension.clone() - Atom::one())
            * (dimension.clone() + Atom::num(2));
        let diagonal = (dimension.clone() + Atom::one()) / denominator.clone();
        let crossed = -Atom::one() / denominator;
        let aa_bb = dot(&k, &k) * dot(&q, &q);
        let ab_ab = dot(&k, &q).pow(Atom::num(2));
        let pp_rr = dot(&p, &p) * dot(&r, &r);
        let pr_pr = dot(&p, &r).pow(Atom::num(2));
        let expected = diagonal.clone() * &aa_bb * &pp_rr
            + Atom::num(2) * crossed.clone() * &aa_bb * &pr_pr
            + Atom::num(2) * crossed.clone() * &ab_ab * &pp_rr
            + Atom::num(2) * (diagonal + crossed) * ab_ab * pr_pr;
        let difference = (result.expression().unwrap().into_expression() - &expected).together();
        assert!(
            difference.is_zero(),
            "result: {}\nexpected: {expected}\ndifference: {difference}",
            result.expression().unwrap().into_expression()
        );
    }

    #[test]
    fn rank_twenty_repeated_projection_stays_one_term() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D20"));
        let mut input = Atom::one();
        for position in 0..20 {
            let index = Atom::var(symbol!(format!("feynkit_tensor_test::rho{position}")));
            input *= indexed(k, 1, &dimension, &index);
            input *= indexed(p, 1, &dimension, &index);
        }
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert_eq!(result.max_rank(), 20);
        assert_eq!(result.terms().len(), 1);
        assert_eq!(
            result.terms()[0]
                .projector_orbit()
                .unwrap()
                .labeled_pairings(),
            654_729_075
        );
        let k = compact(k, 1, &dimension);
        let p = compact(p, 1, &dimension);
        let denominator = (0..10).fold(Atom::one(), |product, offset| {
            product * (dimension.clone() + Atom::num(2 * offset))
        });
        let expected = Atom::num(654_729_075_i64)
            * dot(&k, &k).pow(Atom::num(10))
            * dot(&p, &p).pow(Atom::num(10))
            / denominator;
        assert!(
            (result.expression().unwrap().into_expression() - expected)
                .together()
                .is_zero()
        );
    }

    #[test]
    fn rank_twenty_partial_internal_and_projector_symmetries_are_practical() {
        let started = std::time::Instant::now();
        let (k, q, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D20_partial"));
        let mut input = Atom::one();
        for position in 0..20 {
            let index = Atom::var(symbol!(format!(
                "feynkit_tensor_test::partial_rho{position}"
            )));
            let integrated = if position < 2 { k } else { q };
            let projector = if position < 18 { p } else { r };
            input *= indexed(integrated, 1, &dimension, &index);
            input *= indexed(projector, 1, &dimension, &index);
        }
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let elapsed = started.elapsed();

        // Rank 20 has 654,729,075 labeled pairings. Reaching this compact
        // four-term result proves that both partial symmetries used the orbit
        // projector without relying on a default work-count cutoff.
        assert_eq!(result.max_rank(), 20);
        assert_eq!(result.terms().len(), 4);
        assert!(result.is_fully_contracted());
        let internal_sizes = result
            .terms()
            .iter()
            .map(|term| term.integrated_orbit().unwrap().labeled_pairings())
            .collect::<BTreeSet<_>>();
        let projector_sizes = result
            .terms()
            .iter()
            .map(|term| term.projector_orbit().unwrap().labeled_pairings())
            .collect::<BTreeSet<_>>();
        assert_eq!(internal_sizes, BTreeSet::from([34_459_425, 620_269_650]));
        assert_eq!(projector_sizes, internal_sizes);

        // Setting every scalar invariant to one sums every entry of
        // S_A^T W S_B.  Every W row has the closed isotropic row sum.
        let coefficient_sum = result
            .terms()
            .iter()
            .fold(Atom::Zero, |sum, term| sum + term.coefficient().clone());
        let expected_sum = Atom::num(654_729_075_i64)
            * OrthogonalWeingarten::isotropic_pairing_coefficient(10, dimension).unwrap();
        assert!((coefficient_sum - expected_sum).together().is_zero());
        eprintln!("public rank-20 [2,18] x [2,18] reduction: {elapsed:?}");
    }

    #[test]
    fn isotropic_fast_paths_match_labeled_projectors_at_low_rank() {
        let (k, q, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_fast_oracle"));
        let s = vector_symbol!("feynkit_tensor_test::s_fast_oracle");
        let indices = (0..6)
            .map(|position| {
                Atom::var(symbol!(format!(
                    "feynkit_tensor_test::fast_oracle_mu{position}"
                )))
            })
            .collect::<Vec<_>>();

        let internal_compact = vec![compact(k, 1, &dimension); 6];
        let projector_compact = [p, p, r, r, s, s].map(|head| compact(head, 1, &dimension));
        let mut isotropic_internal_input = Atom::one();
        for position in 0..6 {
            isotropic_internal_input *= indexed(k, 1, &dimension, &indices[position]);
            isotropic_internal_input *= indexed(
                [p, p, r, r, s, s][position],
                1,
                &dimension,
                &indices[position],
            );
        }
        let isotropic_internal = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(
                &SymbolicTensor::infer((isotropic_internal_input.as_view()).to_owned()).unwrap(),
            )
            .unwrap()
            .expression()
            .unwrap()
            .into_expression();
        let brute_internal =
            brute_projected_expression(&internal_compact, &projector_compact, &dimension);
        assert!((isotropic_internal - brute_internal).together().is_zero());

        let internal_heads = [k, k, q, q, q, q];
        let integrated_compact = internal_heads.map(|head| compact(head, 1, &dimension));
        let projector_compact = vec![compact(p, 1, &dimension); 6];
        let mut isotropic_projector_input = Atom::one();
        for position in 0..6 {
            isotropic_projector_input *=
                indexed(internal_heads[position], 1, &dimension, &indices[position]);
            isotropic_projector_input *= indexed(p, 1, &dimension, &indices[position]);
        }
        let isotropic_projector = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .reduce(
                &SymbolicTensor::infer((isotropic_projector_input.as_view()).to_owned()).unwrap(),
            )
            .unwrap()
            .expression()
            .unwrap()
            .into_expression();
        let brute_projector =
            brute_projected_expression(&integrated_compact, &projector_compact, &dimension);
        assert!((isotropic_projector - brute_projector).together().is_zero());
    }

    #[test]
    fn rank_twelve_multiloop_projection_uses_orbits() {
        let (k, q, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D12"));
        let third = vector_symbol!("feynkit_tensor_test::l");
        let mut input = Atom::one();
        for position in 0..12 {
            let index = Atom::var(symbol!(format!("feynkit_tensor_test::tau{position}")));
            let internal = match position / 4 {
                0 => k,
                1 => q,
                _ => third,
            };
            input *= indexed(internal, 1, &dimension, &index);
            input *= indexed(p, 1, &dimension, &index);
        }
        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .with_integrated_head(q)
            .with_integrated_head(third)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert!(result.terms().len() < 100);
        assert_eq!(
            result
                .terms()
                .iter()
                .map(|term| term.integrated_orbit().unwrap().labeled_pairings())
                .sum::<u128>(),
            10_395
        );
    }

    #[test]
    fn rank_twelve_partial_symmetries_use_orbit_projector() {
        let (k, q, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D12_generic"));
        let mut input = Atom::one();
        for position in 0..12 {
            let index = Atom::var(symbol!(format!(
                "feynkit_tensor_test::generic_mu{position}"
            )));
            let internal = if position < 6 { k } else { q };
            let external = if position % 2 == 0 { p } else { r };
            input *= indexed(internal, 1, &dimension, &index);
            input *= indexed(external, 1, &dimension, &index);
        }
        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .with_integrated_head(q)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let materialized = result.expression().unwrap().into_expression();
        eprintln!(
            "generic rank-12 output: {} compact terms, {} canonical characters",
            result.terms().len(),
            materialized.to_canonical_string().len()
        );
        assert_eq!(result.max_rank(), 12);
        assert_eq!(result.terms().len(), 16);
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn already_contracted_integrated_vectors_become_a_dot() {
        let (k, q, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::Dc"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::muc"));
        let input = indexed(k, 1, &dimension, &mu) * indexed(q, 2, &dimension, &mu);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .with_integrated_head(q)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert_eq!(
            result.expression().unwrap().into_expression(),
            dot(&compact(k, 1, &dimension), &compact(q, 2, &dimension))
        );
    }

    #[test]
    fn external_only_contractions_use_spenso_dots() {
        let (k, _, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_external_dot"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::external_dot_mu"));
        let input = indexed(p, 1, &dimension, &mu) * indexed(r, 2, &dimension, &mu);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();

        assert_eq!(
            result.expression().unwrap().into_expression(),
            dot(&compact(p, 1, &dimension), &compact(r, 2, &dimension))
        );
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn external_only_free_indices_prevent_scalar_graph_splitting() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_external_free"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::external_free_mu"));
        let input = indexed(p, 1, &dimension, &mu);
        let diagram = diagram_with_numerator("external-free", input.clone());
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);

        let result = reducer
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert!(!result.is_fully_contracted());
        assert!(matches!(
            diagram.reduce_tensor_graphs(&reducer),
            Err(TensorReductionError::ResidualFreeIndices)
        ));
    }

    #[test]
    fn external_only_triple_index_is_ambiguous() {
        let (k, _, p, r) = vectors();
        let third = vector_symbol!("feynkit_tensor_test::s_external_triple");
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_external_triple"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::external_triple_mu"));
        let input = indexed(p, 1, &dimension, &mu)
            * indexed(r, 1, &dimension, &mu)
            * indexed(third, 1, &dimension, &mu);
        assert!(SymbolicTensor::infer(input.clone()).is_err());
        assert!(matches!(
            TensorReducer::new(dimension)
                .with_integrated_head(k)
                .parse_monomial(input.as_view())
                .map(|_| ()),
            Err(TensorReductionError::AmbiguousIndex {
                integrated: 0,
                outside: 3,
                ..
            })
        ));
    }

    #[test]
    fn opaque_free_minkowski_slots_prevent_scalar_graph_splitting() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_opaque_free"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::opaque_free_mu"));
        let opaque = FunctionBuilder::new(spenso::broadcast_symbol!("feynkit_tensor_test::opaque"))
            .add_arg(indexed(p, 1, &dimension, &mu))
            .finish();
        let diagram = diagram_with_numerator("opaque-free", opaque.clone());
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);

        let result = reducer
            .reduce(&SymbolicTensor::infer((opaque.as_view()).to_owned()).unwrap())
            .unwrap();
        assert!(!result.is_fully_contracted());
        assert!(matches!(
            diagram.reduce_tensor_graphs(&reducer),
            Err(TensorReductionError::ResidualFreeIndices)
        ));
    }

    #[test]
    fn ambiguous_opaque_consumers_cannot_share_a_recognized_vector_index() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_opaque_consumer"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::opaque_consumer_mu"));
        let opaque = FunctionBuilder::new(spenso::broadcast_symbol!(
            "feynkit_tensor_test::opaque_consumer"
        ))
        .add_arg(indexed(p, 1, &dimension, &mu))
        .finish();
        let input = indexed(k, 1, &dimension, &mu) * indexed(k, 2, &dimension, &mu) * opaque;
        assert!(SymbolicTensor::infer(input.clone()).is_err());
        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .parse_monomial(input.as_view())
            .map(|_| ());
        assert!(
            matches!(
                result,
                Err(TensorReductionError::AmbiguousMinkowskiIndex { occurrences: 3, .. })
            ),
            "unexpected ambiguous opaque-consumer result: {result:?}"
        );
    }

    #[test]
    fn opaque_external_tensor_legs_are_preserved_by_the_projector() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_opaque_projector"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::opaque_projector_mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::opaque_projector_nu"));
        let tensor = spenso::tensor_symbol!("feynkit_tensor_test::opaque_projector_tensor");
        let tensor_mu = FunctionBuilder::new(tensor)
            .add_arg(minkowski_slot(&dimension, &mu))
            .finish();
        let loop_squared = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension));
        let input = indexed(k, 1, &dimension, &mu)
            * indexed(k, 1, &dimension, &nu)
            * &tensor_mu
            * indexed(p, 1, &dimension, &nu);
        let expected = loop_squared
            * dot(
                &FunctionBuilder::new(tensor)
                    .add_arg(
                        FunctionBuilder::new(symbol!("spenso::mink"))
                            .add_arg(&dimension)
                            .finish(),
                    )
                    .finish(),
                &compact(p, 1, &dimension),
            )
            / &dimension;

        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let difference = (result.expression().unwrap().into_expression() - expected).together();
        assert!(
            difference.is_zero(),
            "opaque external tensor projector differs by {difference}"
        );
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn repeated_opaque_tensor_slot_is_reported_as_ambiguous() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_opaque_power"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::opaque_power_mu"));
        let spin = Atom::var(symbol!("feynkit_tensor_test::opaque_power_spin"));
        let tensor = FunctionBuilder::new(spenso::tensor_symbol!(
            "feynkit_tensor_test::opaque_power_tensor"
        ))
        .add_arg(minkowski_slot(&dimension, &spin))
        .add_arg(minkowski_slot(&dimension, &mu))
        .finish();
        let input = indexed(k, 1, &dimension, &mu) * tensor.pow(Atom::num(2));

        assert!(SymbolicTensor::infer(input.clone()).is_err());
        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .parse_monomial(input.as_view())
            .map(|_| ());
        assert!(
            matches!(
                result,
                Err(TensorReductionError::AmbiguousMinkowskiIndex { occurrences: 3, .. })
            ),
            "unexpected repeated opaque-slot result: {result:?}"
        );
    }

    #[test]
    fn symbolic_opaque_tensor_powers_are_rejected() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_opaque_symbolic_power"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::opaque_symbolic_power_mu"));
        let spin = Atom::var(symbol!("feynkit_tensor_test::opaque_symbolic_power_spin"));
        let exponent = Atom::var(symbol!("feynkit_tensor_test::opaque_symbolic_power_n"));
        let tensor = FunctionBuilder::new(spenso::tensor_symbol!(
            "feynkit_tensor_test::opaque_symbolic_power_tensor"
        ))
        .add_arg(minkowski_slot(&dimension, &spin))
        .add_arg(minkowski_slot(&dimension, &mu))
        .finish();
        let input = indexed(k, 1, &dimension, &mu) * tensor.pow(exponent);

        assert!(SymbolicTensor::infer(input.clone()).is_err());
        assert!(matches!(
            TensorReducer::new(dimension)
                .with_integrated_head(k)
                .parse_monomial(input.as_view())
                .map(|_| ()),
            Err(TensorReductionError::InvalidVectorPower(_))
        ));
    }

    #[test]
    fn parser_rejects_large_vector_powers_before_materialization() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_large_power"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::large_power_mu"));
        let input = indexed(k, 1, &dimension, &mu).pow(Atom::num(1_000_000));
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);
        let result = reducer.parse_monomial(input.as_view());
        assert!(matches!(
            &result,
            Err(TensorReductionError::AmbiguousIndex {
                integrated: 1_000_000,
                outside: 0,
                ..
            })
        ));
    }

    #[test]
    fn scalar_spectators_do_not_consume_the_selected_term_budget() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_expansion_budget"));
        let variables = (0..6)
            .map(|index| {
                Atom::var(symbol!(&format!(
                    "feynkit_tensor_test::expansion_budget_{index}"
                )))
            })
            .collect::<Vec<_>>();
        let input = (variables[0].clone() + &variables[1])
            * (variables[2].clone() + &variables[3])
            * (variables[4].clone() + &variables[5]);
        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .with_output_term_limit(7)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert_eq!(result.terms().len(), 1);
        assert_eq!(result.expression().unwrap().into_expression(), input);
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn reduction_keeps_factored_scalar_and_foreign_tensor_spectators() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_factored_coefficient"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::coefficient_mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::coefficient_nu"));
        let x = Atom::var(symbol!("feynkit_tensor_test::coefficient_x"));
        let y = Atom::var(symbol!("feynkit_tensor_test::coefficient_y"));
        let foreign = spenso::tensor!(
            spenso::tensor_symbol!("feynkit_tensor_test::coefficient_foreign"),
            spenso::euc!(3, 93201)
        );
        let spectator = (&x + &y).pow(30) * foreign;
        let input = &spectator
            * indexed(k, 1, &dimension, &mu)
            * indexed(k, 1, &dimension, &nu)
            * indexed(p, 1, &dimension, &mu)
            * indexed(p, 1, &dimension, &nu);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let expected = spectator
            * dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(p, 1, &dimension))
            / dimension;
        assert_eq!(result.expression().unwrap().into_expression(), expected);
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn incomplete_selected_tensor_collection_is_refused() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_refused_collection"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::refused_mu"));
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&calls);
        let hook = tensor_symbol!(
            "feynkit_tensor_test::refused_collection_hook",
            norm = move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        );
        let value = FunctionBuilder::new(hook)
            .add_arg(minkowski_slot(&dimension, &mu))
            .finish();
        let input = FunctionBuilder::new(broadcast_symbol!(
            "feynkit_tensor_test::refused_collection_broadcast"
        ))
        .add_arg(value)
        .finish();
        let source = SymbolicTensor::infer(input).unwrap();
        let before = source.expression().as_view();
        calls.store(0, Ordering::Relaxed);
        // Rebuilding this selected scope would replay a user normalizer. The
        // existing collector must refuse it without publishing partial rows.
        let collection_error = source
            .coefficient_list(TensorCollectFilter::Reps([LibraryRep::from(Minkowski {})]))
            .unwrap_err()
            .to_string();
        assert!(collection_error.contains("did not complete"));
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        let error = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .reduce(&source)
            .unwrap_err();
        assert!(
            matches!(error, TensorReductionError::TensorAlgebra(ref message)
            if message == &collection_error),
            "{error}"
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
        assert_eq!(source.expression().as_view(), before);
    }

    #[test]
    fn finite_selected_tensor_polynomial_completes_without_diagnostic_limits() {
        let (k, q, p, r) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_completed_collection"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::completed_mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::completed_nu"));
        let left = indexed(k, 1, &dimension, &mu) * indexed(p, 1, &dimension, &mu);
        let right = indexed(q, 1, &dimension, &nu) * indexed(r, 1, &dimension, &nu);
        let source = SymbolicTensor::infer((left + right).pow(2)).unwrap();
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&source)
            .unwrap();
        // The odd k term vanishes; <(k.p)^2> = k^2 p^2 / D.
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(p, 1, &dimension))
            / &dimension
            + dot(&compact(q, 1, &dimension), &compact(r, 1, &dimension)).pow(2);
        assert!(
            (result.expression().unwrap().into_expression() - expected)
                .together()
                .is_zero()
        );
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn explicit_metric_contractions_are_normalized_automatically() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::Dg"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::mug"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::nug"));
        let input = metric(&dimension, &mu, &nu)
            * indexed(k, 1, &dimension, &mu)
            * indexed(k, 1, &dimension, &nu);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        assert_eq!(
            result.expression().unwrap().into_expression(),
            dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
        );
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn products_of_generated_vertex_sums_are_sewn_by_the_shared_collector() {
        let (k, _, p, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::D_vertex_sums"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::vertex_sums_mu"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::vertex_sums_nu"));
        let a = Atom::var(symbol!("feynkit_tensor_test::vertex_sums_a"));
        let b = Atom::var(symbol!("feynkit_tensor_test::vertex_sums_b"));
        let c = Atom::var(symbol!("feynkit_tensor_test::vertex_sums_c"));
        let d = Atom::var(symbol!("feynkit_tensor_test::vertex_sums_d"));
        let left = metric(&dimension, &mu, &a) * indexed(k, 1, &dimension, &a)
            + metric(&dimension, &mu, &b) * indexed(k, 1, &dimension, &b);
        let right = metric(&dimension, &nu, &c) * indexed(k, 1, &dimension, &c)
            + metric(&dimension, &nu, &d) * indexed(k, 1, &dimension, &d);
        let input = left * right * indexed(p, 1, &dimension, &mu) * indexed(p, 1, &dimension, &nu);

        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap()
            .expression()
            .unwrap()
            .into_expression();
        let expected = Atom::num(4)
            * dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(p, 1, &dimension))
            / dimension;
        assert_eq!(result, expected);
    }

    #[test]
    fn metric_normalization_preserves_residual_free_indices() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::Dgr"));
        let contracted = Atom::var(symbol!("feynkit_tensor_test::mugr"));
        let left_free = Atom::var(symbol!("feynkit_tensor_test::rhogr"));
        let right_free = Atom::var(symbol!("feynkit_tensor_test::nugr"));
        let input = metric(&dimension, &contracted, &left_free)
            * indexed(k, 1, &dimension, &contracted)
            * indexed(k, 1, &dimension, &right_free);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * metric(&dimension, &left_free, &right_free)
            / dimension.clone();
        let expression = result.expression().unwrap().into_expression();
        assert_eq!(expression, expected);
        assert!(
            count_minkowski_index(expression.as_view(), &dimension, &left_free, 1).unwrap() > 0
        );
        assert!(
            count_minkowski_index(expression.as_view(), &dimension, &right_free, 1).unwrap() > 0
        );
        assert!(!result.is_fully_contracted());
    }

    #[test]
    fn free_indices_emit_dimensioned_spenso_metrics() {
        let (k, _, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::Df"));
        let mu = Atom::var(symbol!("feynkit_tensor_test::muf"));
        let nu = Atom::var(symbol!("feynkit_tensor_test::nuf"));
        let input = indexed(k, 1, &dimension, &mu) * indexed(k, 1, &dimension, &nu);
        let result = TensorReducer::new(dimension.clone())
            .with_integrated_head(k)
            .reduce(&SymbolicTensor::infer((input.as_view()).to_owned()).unwrap())
            .unwrap();
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * FunctionBuilder::new(ETS.metric)
                .add_arg(minkowski_slot(&dimension, &mu))
                .add_arg(minkowski_slot(&dimension, &nu))
                .finish()
            / dimension.clone();
        let expression = result.expression().unwrap().into_expression();
        assert_eq!(expression, expected);
        assert!(contains_metric_pair(
            expression.as_view(),
            &dimension,
            &mu,
            &nu
        ));
        assert!(!result.is_fully_contracted());
    }

    #[test]
    fn contraction_orbit_multiplicities_sum_to_all_pairings() {
        let (k, q, _, _) = vectors();
        let dimension = Atom::var(symbol!("feynkit_tensor_test::Do"));
        let colors = vec![
            compact(k, 1, &dimension),
            compact(k, 1, &dimension),
            compact(q, 1, &dimension),
            compact(q, 1, &dimension),
            compact(q, 1, &dimension),
            compact(q, 1, &dimension),
        ];
        let orbits = contraction_orbits(&colors, Some(100)).unwrap();
        assert_eq!(orbits.len(), 2);
        assert_eq!(
            orbits
                .iter()
                .map(ContractionOrbit::labeled_pairings)
                .sum::<u128>(),
            15
        );
    }
}
