//! Spenso-facing tensor reduction of vacuum numerators.
//!
//! The reducer consumes ordinary Symbolica atoms whose rank-one tensors use
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

use feynkit_graph::FeynmanDiagram;
use idenso::shorthands::metric::MetricSimplifier;
use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use spenso::structure::representation::{BaseRepName, Minkowski};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use thiserror::Error;

use crate::{
    OrthogonalWeingarten, WeingartenError,
    orbit_projector::{OrbitProjector, OrbitProjectorError},
};

const DEFAULT_PAIRING_LIMIT: u128 = 20_000;
const DEFAULT_PAIRING_PRODUCT_LIMIT: u128 = 120_000_000;
const DEFAULT_OUTPUT_TERM_LIMIT: usize = 100_000;

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
///     .reduce(numerator.as_view())?;
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
    pub fn invariant(&self) -> Atom {
        let mut invariant = Atom::one();
        let mut alpha = 0;
        for left in 0..self.vectors.len() {
            for right in left..self.vectors.len() {
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
        invariant
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
/// let numerator = parse!(
///     "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))
///      *p(spenso::mink(D,mu))*p(spenso::mink(D,nu))"
/// );
/// let reduction = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .reduce(numerator.as_view())?;
/// let scalar_graph_numerator = reduction.terms()[0].expression();
/// println!("{scalar_graph_numerator}");
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug)]
pub struct TensorReductionTerm {
    coefficient: Atom,
    tensor: Atom,
    integrated_orbit: Option<ContractionOrbit>,
    projector_orbit: Option<ContractionOrbit>,
}

impl TensorReductionTerm {
    /// Exact rational coefficient, including orbit multiplicities.
    pub fn coefficient(&self) -> &Atom {
        &self.coefficient
    }

    /// Scalar invariant or residual metric tensor multiplying the coefficient.
    pub fn tensor(&self) -> &Atom {
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

    /// Materialize this contribution as an ordinary Symbolica expression.
    pub fn expression(&self) -> Atom {
        &self.coefficient * &self.tensor
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
/// let numerator = parse!(
///     "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))
///      *p(spenso::mink(D,mu))*p(spenso::mink(D,nu))"
/// );
/// let scalar = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .reduce(numerator.as_view())?
///     .into_expression();
/// println!("{scalar}");
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug)]
pub struct TensorReduction {
    terms: Vec<TensorReductionTerm>,
    max_rank: usize,
    fully_contracted: bool,
}

impl TensorReduction {
    /// Compact reduction terms before their final sum is expanded.
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
    /// fn scalar_numerator(reduction: &TensorReduction) -> Atom {
    ///     reduction.expression()
    /// }
    /// ```
    pub fn expression(&self) -> Atom {
        self.terms
            .iter()
            .fold(Atom::Zero, |sum, term| sum + term.expression())
    }

    /// Consume the compact result and materialize its expression.
    pub fn into_expression(self) -> Atom {
        self.terms
            .into_iter()
            .fold(Atom::Zero, |sum, term| sum + term.coefficient * term.tensor)
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
/// let numerator = parse!(
///     "spenso::g(spenso::mink(D,mu),spenso::mink(D,nu))
///      *k(spenso::mink(D,mu))*k(spenso::mink(D,nu))"
/// );
/// let reduced = TensorReducer::new(parse!("D"))
///     .with_integrated_head(symbol!("k"))
///     .reduce(numerator.as_view())?;
/// println!("{}", reduced.expression());
/// # Ok::<(), feynkit_tensor::TensorReductionError>(())
/// ```
#[derive(Clone, Debug)]
pub struct TensorReducer {
    dimension: Atom,
    integrated_heads: BTreeSet<String>,
    integrated_vectors: BTreeSet<Atom>,
    pairing_limit: u128,
    pairing_product_limit: u128,
    output_term_limit: usize,
}

impl TensorReducer {
    /// Construct a reducer in `dimension` with no integrated vectors selected.
    ///
    /// `dimension` must exactly match the first argument of every input
    /// `spenso::mink(dimension,index)` slot.  Thus generated four-dimensional
    /// FeynKit rules use `4`, while a symbolic R-operation numerator whose
    /// slots carry `D` uses `D`. Mixed high-rank projectors should retain a
    /// symbolic dimension until after reduction: a fixed positive integer
    /// dimension below half the tensor rank makes the universal metric basis
    /// singular. The all-equal isotropic fast path remains available there.
    ///
    /// Add selectors with [`Self::with_integrated_head`] or
    /// [`Self::with_integrated_vector`].
    pub fn new(dimension: Atom) -> Self {
        Self {
            dimension,
            integrated_heads: BTreeSet::new(),
            integrated_vectors: BTreeSet::new(),
            pairing_limit: DEFAULT_PAIRING_LIMIT,
            pairing_product_limit: DEFAULT_PAIRING_PRODUCT_LIMIT,
            output_term_limit: DEFAULT_OUTPUT_TERM_LIMIT,
        }
    }

    /// Construct a reducer for every `FeynKit::Momentum` vector.
    ///
    /// This constructor stores the qualified head name without registering a
    /// bare Symbolica symbol, so FeynKit's tensor tags remain initialization
    /// order independent.
    pub fn feynkit(dimension: Atom) -> Self {
        Self::new(dimension).with_integrated_head_name("FeynKit::Momentum")
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

    /// Select one exact compact vector, such as `K(1,spenso::mink(D))`.
    pub fn with_integrated_vector(mut self, vector: Atom) -> Self {
        self.integrated_vectors.insert(vector);
        self
    }

    /// Bound labeled pairings generated by the generic/free-index path.
    ///
    /// Symmetric contraction-orbit fast paths do not consume this budget.
    pub fn with_pairing_limit(mut self, limit: u128) -> Self {
        self.pairing_limit = limit;
        self
    }

    /// Bound the labeled matching product in the residual-free-index path.
    pub fn with_pairing_product_limit(mut self, limit: u128) -> Self {
        self.pairing_product_limit = limit;
        self
    }

    /// Bound the number of distinct materialized invariant terms.
    pub fn with_output_term_limit(mut self, limit: usize) -> Self {
        self.output_term_limit = limit;
        self
    }

    /// Dimension used by Spenso slots and projector coefficients.
    pub fn dimension(&self) -> &Atom {
        &self.dimension
    }

    /// Distribute additive factors into independently reducible summands.
    ///
    /// The same output-term limit configured on this reducer is enforced
    /// while distributing, so callers can attach term-local tensor metadata
    /// without first performing an unbounded symbolic expansion.
    ///
    /// # Examples
    ///
    /// ```
    /// use feynkit_tensor::TensorReducer;
    /// use symbolica::parse;
    ///
    /// let reducer = TensorReducer::new(parse!("D"));
    /// let terms = reducer.distribute_summands(parse!("(a+b)*(c+d)").as_view())?;
    /// assert_eq!(terms.len(), 4);
    /// # Ok::<(), feynkit_tensor::TensorReductionError>(())
    /// ```
    pub fn distribute_summands(
        &self,
        expression: AtomView<'_>,
    ) -> Result<Vec<Atom>, TensorReductionError> {
        bounded_expand_summands(expression, self.output_term_limit)
    }

    /// Reduce a tensor expression to compact invariant terms.
    ///
    /// Explicit Spenso metric contractions are normalized first with Idenso.
    /// This contracts chains such as `g(mu,nu) K(mu) K(nu)` while preserving
    /// any residual free index in its full `spenso::mink(D,index)` slot.  Each
    /// summand is then treated independently.  Odd-rank vacuum tensors vanish;
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
    /// let numerator = parse!(
    ///     "k(spenso::mink(D,mu))*k(spenso::mink(D,nu))
    ///      *p(spenso::mink(D,mu))*p(spenso::mink(D,nu))"
    /// );
    /// let reduction = TensorReducer::new(parse!("D"))
    ///     .with_integrated_head(symbol!("k"))
    ///     .reduce(numerator.as_view())?;
    /// assert!(reduction.is_fully_contracted());
    /// # Ok::<(), feynkit_tensor::TensorReductionError>(())
    /// ```
    pub fn reduce(
        &self,
        expression: AtomView<'_>,
    ) -> Result<TensorReduction, TensorReductionError> {
        if self.integrated_heads.is_empty() && self.integrated_vectors.is_empty() {
            return Err(TensorReductionError::NoIntegratedVectorsSelected);
        }
        // Generated Feynman rules commonly contain products of vertex sums.
        // Distribute those sums before metric simplification so Idenso sees
        // each complete tensor monomial and can sew all propagator indices.
        // The local distributor enforces the configured term budget while it
        // works, rather than allocating an unbounded global `expand()` first.
        let distributed = bounded_expand_summands(expression, self.output_term_limit)?;
        let mut summands = Vec::with_capacity(distributed.len());
        for summand in distributed {
            let normalized = summand.simplify_metrics().metric_shorthand_to_dot();
            let normalized = bounded_expand_summands(normalized.as_view(), self.output_term_limit)?;
            let term_count = summands.len().checked_add(normalized.len()).ok_or(
                TensorReductionError::ExpansionLimit {
                    terms: usize::MAX,
                    limit: self.output_term_limit,
                },
            )?;
            if term_count > self.output_term_limit {
                return Err(TensorReductionError::ExpansionLimit {
                    terms: term_count,
                    limit: self.output_term_limit,
                });
            }
            summands.extend(normalized);
        }

        let mut terms = Vec::new();
        let mut max_rank = 0;
        for summand in summands {
            let monomial = self.parse_monomial(summand.as_view())?;
            max_rank = max_rank.max(monomial.integrated.len());
            if monomial.integrated.len() > OrthogonalWeingarten::MAX_RANK {
                return Err(TensorReductionError::UnsupportedRank {
                    rank: monomial.integrated.len(),
                    maximum: OrthogonalWeingarten::MAX_RANK,
                });
            }
            if monomial.integrated.len() % 2 == 1 {
                continue;
            }
            terms.extend(self.reduce_monomial(monomial)?);
            if terms.len() > self.output_term_limit {
                return Err(TensorReductionError::OutputLimit {
                    terms: terms.len(),
                    limit: self.output_term_limit,
                });
            }
        }
        let fully_contracted = terms.iter().try_fold(true, |fully_contracted, term| {
            Ok::<bool, TensorReductionError>(
                fully_contracted && !has_dangling_minkowski_indices(&term.tensor, &self.dimension)?,
            )
        })?;

        Ok(TensorReduction {
            terms,
            max_rank,
            fully_contracted,
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
        let mut integrated_rank = 0_usize;

        for factor in factors {
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
                if exponent > OrthogonalWeingarten::MAX_RANK {
                    return Err(TensorReductionError::UnsupportedRank {
                        rank: exponent,
                        maximum: OrthogonalWeingarten::MAX_RANK,
                    });
                }
                if integrated {
                    integrated_rank = integrated_rank.checked_add(exponent).ok_or(
                        TensorReductionError::UnsupportedRank {
                            rank: usize::MAX,
                            maximum: OrthogonalWeingarten::MAX_RANK,
                        },
                    )?;
                    if integrated_rank > OrthogonalWeingarten::MAX_RANK {
                        return Err(TensorReductionError::UnsupportedRank {
                            rank: integrated_rank,
                            maximum: OrthogonalWeingarten::MAX_RANK,
                        });
                    }
                }
                indexed_vectors.extend(std::iter::repeat_n(vector, exponent));
                continue;
            }
            if let Some(vector) = self.indexed_vector(factor)? {
                if self.is_integrated(&vector) {
                    integrated_rank += 1;
                    if integrated_rank > OrthogonalWeingarten::MAX_RANK {
                        return Err(TensorReductionError::UnsupportedRank {
                            rank: integrated_rank,
                            maximum: OrthogonalWeingarten::MAX_RANK,
                        });
                    }
                }
                indexed_vectors.push(vector);
            } else {
                let factor = factor.to_owned();
                opaque_factors.push(factor.clone());
                scalar_factors.push(factor);
            }
        }

        let mut by_index: BTreeMap<Atom, IndexOccupants> = BTreeMap::new();
        for vector in indexed_vectors {
            let integrated = self.is_integrated(&vector);
            let occupants = by_index.entry(vector.index.clone()).or_default();
            if integrated {
                occupants.integrated.push(vector.compact);
            } else {
                occupants.outside.push(vector.compact);
            }
        }

        let mut integrated = Vec::new();
        let mut outside = Vec::new();
        for (index, occupants) in by_index {
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
                        slot: minkowski_slot(&self.dimension, &index),
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
        self.integrated_heads.contains(vector.head.get_name())
            || self.integrated_vectors.contains(&vector.compact)
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
        if representation_arguments[0] != self.dimension.as_view() {
            return Err(TensorReductionError::DimensionMismatch {
                expected: self.dimension.clone(),
                found: representation_arguments[0].to_owned(),
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
        }))
    }

    fn reduce_monomial(
        &self,
        monomial: TensorMonomial,
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        let rank = monomial.integrated.len();
        if rank == 0 {
            return Ok(vec![TensorReductionTerm {
                coefficient: Atom::one(),
                tensor: monomial.scalar,
                integrated_orbit: None,
                projector_orbit: None,
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
        let pairing_count = perfect_matching_count(rank);
        if pairing_count > self.pairing_limit {
            return Err(TensorReductionError::PairingLimit {
                rank,
                pairings: pairing_count,
                limit: self.pairing_limit,
            });
        }
        let pairing_products = pairing_count.saturating_mul(pairing_count);
        if pairing_products > self.pairing_product_limit {
            return Err(TensorReductionError::PairingProductLimit {
                rank,
                products: pairing_products,
                limit: self.pairing_product_limit,
            });
        }

        let pairings = perfect_matchings(rank);
        let internal_classes = classified_pairings(&monomial.integrated, &pairings)?;
        let outside_classes = classified_outside(&monomial.outside, &pairings)?;
        let output_classes = internal_classes.orbits.len() * outside_classes.atoms.len();
        if output_classes > self.output_term_limit {
            return Err(TensorReductionError::OutputLimit {
                terms: output_classes,
                limit: self.output_term_limit,
            });
        }

        let table = OrthogonalWeingarten::new(pair_count, self.dimension.clone())?;
        let mut coefficient_by_key = BTreeMap::new();
        let mut coset_index = BTreeMap::new();
        for (index, (coset, coefficient)) in table.coefficients().enumerate() {
            coset_index.insert(coset_key_from_parts(coset.parts()), index);
            coefficient_by_key.insert(index, coefficient.clone());
        }
        let coefficient_count = coefficient_by_key.len();
        let mut counts =
            vec![
                0_u128;
                internal_classes.orbits.len() * outside_classes.atoms.len() * coefficient_count
            ];
        for (left_index, left) in pairings.iter().enumerate() {
            let internal = internal_classes.class_of_pairing[left_index];
            for (right_index, right) in pairings.iter().enumerate() {
                let outside = outside_classes.class_of_pairing[right_index];
                let key = matching_coset_key(left, right);
                let coefficient = coset_index[&key];
                let offset = (internal * outside_classes.atoms.len() + outside) * coefficient_count
                    + coefficient;
                counts[offset] += 1;
            }
        }

        let mut terms = Vec::new();
        for (integrated_index, integrated_orbit) in internal_classes.orbits.into_iter().enumerate()
        {
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
                    tensor: monomial.scalar.clone()
                        * integrated_orbit.invariant()
                        * outside_classes.atoms[outside_index].clone(),
                    integrated_orbit: Some(integrated_orbit.clone()),
                    projector_orbit: outside_classes.orbits[outside_index].clone(),
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

        Ok(projector
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
                TensorReductionTerm {
                    coefficient: term.coefficient().clone(),
                    tensor: monomial.scalar.clone()
                        * integrated_orbit.invariant()
                        * projector_orbit.invariant(),
                    integrated_orbit: Some(integrated_orbit),
                    projector_orbit: Some(projector_orbit),
                }
            })
            .collect())
    }

    fn reduce_isotropic_internal(
        &self,
        monomial: TensorMonomial,
        projector_vectors: Option<Vec<Atom>>,
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        let pair_count = monomial.integrated.len() / 2;
        let coefficient =
            OrthogonalWeingarten::isotropic_pairing_coefficient(pair_count, self.dimension.clone());
        let integrated_orbit = repeated_orbit(monomial.integrated[0].clone(), pair_count)?;
        let internal = integrated_orbit.invariant();

        if let Some(projector_vectors) = projector_vectors {
            let orbits = contraction_orbits(&projector_vectors, self.output_term_limit)?;
            return orbits
                .into_iter()
                .map(|projector_orbit| {
                    Ok(TensorReductionTerm {
                        coefficient: coefficient.clone()
                            * atom_from_count(projector_orbit.labeled_pairings)?,
                        tensor: monomial.scalar.clone()
                            * internal.clone()
                            * projector_orbit.invariant(),
                        integrated_orbit: Some(integrated_orbit.clone()),
                        projector_orbit: Some(projector_orbit),
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
                    tensor: monomial.scalar.clone() * internal.clone() * tensor,
                    integrated_orbit: Some(integrated_orbit.clone()),
                    projector_orbit: None,
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
        let coefficient =
            OrthogonalWeingarten::isotropic_pairing_coefficient(pair_count, self.dimension.clone());
        let projector_orbit = repeated_orbit(projector_vectors[0].clone(), pair_count)?;
        let projector = projector_orbit.invariant();
        contraction_orbits(&monomial.integrated, self.output_term_limit)?
            .into_iter()
            .map(|integrated_orbit| {
                Ok(TensorReductionTerm {
                    coefficient: coefficient.clone()
                        * atom_from_count(integrated_orbit.labeled_pairings)?,
                    tensor: monomial.scalar.clone()
                        * integrated_orbit.invariant()
                        * projector.clone(),
                    integrated_orbit: Some(integrated_orbit),
                    projector_orbit: Some(projector_orbit.clone()),
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
///         .map(|reduction| reduction.into_expression())
/// }
///
/// fn scalar_graphs(
///     diagram: &FeynmanDiagram,
/// ) -> Result<Vec<FeynmanDiagram>, TensorReductionError> {
///     diagram.reduce_tensor_graphs(&TensorReducer::feynkit(parse!("4")))
/// }
/// ```
pub trait FeynmanDiagramTensorExt {
    /// Reduce the finalized diagram numerator and external-state projector.
    ///
    /// The request-wide numerator prefactor is scalar by construction and is
    /// intentionally kept separate from the returned tensor reduction.
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
    fn reduce_tensor_numerator(
        &self,
        reducer: &TensorReducer,
    ) -> Result<TensorReduction, TensorReductionError> {
        let projected_numerator = self.numerator() * self.projector();
        reducer.reduce(projected_numerator.as_view())
    }

    fn reduce_tensor_graphs(
        &self,
        reducer: &TensorReducer,
    ) -> Result<Vec<FeynmanDiagram>, TensorReductionError> {
        let reduction = self.reduce_tensor_numerator(reducer)?;
        if !reduction.is_fully_contracted() {
            return Err(TensorReductionError::ResidualFreeIndices);
        }
        reduction
            .terms
            .into_iter()
            .map(|term| term.expression())
            .filter(|numerator| !numerator.is_zero())
            .enumerate()
            .map(|(index, numerator)| {
                self.clone()
                    .with_name(format!("{}.tensor[{index}]", self.name()))
                    .with_projector(Atom::one())
                    .with_numerator(numerator)
                    .map_err(TensorReductionError::Diagram)
            })
            .collect()
    }
}

/// Failures while parsing or materializing a tensor reduction.
#[derive(Debug, Error)]
pub enum TensorReductionError {
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
    /// Distributing factored tensor sums exceeded the materialization budget.
    #[error("tensor numerator expansion produces at least {terms} terms, above limit {limit}")]
    ExpansionLimit { terms: usize, limit: usize },
    /// A combinatorial multiplicity did not fit Symbolica's integer boundary.
    #[error("pairing multiplicity {0} does not fit a signed 64-bit coefficient")]
    MultiplicityOverflow(u128),
    /// The native matching representation is intentionally bounded to the
    /// high-rank target validated by this API.
    #[error("tensor rank {rank} exceeds the supported maximum {maximum}")]
    UnsupportedRank { rank: usize, maximum: usize },
    /// The orthogonal-Weingarten coefficient engine failed.
    #[error(transparent)]
    Weingarten(#[from] WeingartenError),
    /// A scalar contribution could not be represented by a valid native graph.
    #[error(transparent)]
    Diagram(#[from] feynkit_graph::DiagramError),
    /// More than one Einstein contraction reused the same Minkowski index.
    #[error("Minkowski index {index} occurs {occurrences} times in one tensor term")]
    AmbiguousMinkowskiIndex { index: Atom, occurrences: usize },
    /// A symmetry-reduced multiplicity exceeded the exact rank-20 counter.
    #[error("pairing multiplicity overflowed the rank-20 exact counter")]
    CounterOverflow,
    /// Internal and projector color lists must describe the same tensor rank.
    #[error("internal colors have rank {internal}, but projector colors have rank {projector}")]
    ProjectorRankMismatch { internal: usize, projector: usize },
}

impl From<OrbitProjectorError> for TensorReductionError {
    fn from(error: OrbitProjectorError) -> Self {
        match error {
            OrbitProjectorError::ColorLengthMismatch {
                internal,
                projector,
            } => Self::ProjectorRankMismatch {
                internal,
                projector,
            },
            OrbitProjectorError::OddRank(rank) => Self::Weingarten(WeingartenError::OddRank(rank)),
            OrbitProjectorError::UnsupportedRank { rank, maximum } => {
                Self::UnsupportedRank { rank, maximum }
            }
            OrbitProjectorError::OutputLimit { terms, limit } => Self::OutputLimit { terms, limit },
            OrbitProjectorError::MultiplicityOverflow(count) => Self::MultiplicityOverflow(count),
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
}

#[derive(Default)]
struct IndexOccupants {
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

fn bounded_expand_summands(
    expression: AtomView<'_>,
    limit: usize,
) -> Result<Vec<Atom>, TensorReductionError> {
    if limit == 0 {
        return Err(TensorReductionError::ExpansionLimit { terms: 1, limit });
    }
    match expression {
        AtomView::Add(sum) => {
            let mut terms = Vec::new();
            for summand in sum.iter() {
                let expanded = bounded_expand_summands(summand, limit)?;
                let term_count = terms.len().checked_add(expanded.len()).ok_or(
                    TensorReductionError::ExpansionLimit {
                        terms: usize::MAX,
                        limit,
                    },
                )?;
                if term_count > limit {
                    return Err(TensorReductionError::ExpansionLimit {
                        terms: term_count,
                        limit,
                    });
                }
                terms.extend(expanded);
            }
            Ok(terms)
        }
        AtomView::Mul(product) => {
            let mut terms = vec![Atom::one()];
            for factor in product.iter() {
                let expanded = bounded_expand_summands(factor, limit)?;
                terms = bounded_term_product(&terms, &expanded, limit)?;
            }
            Ok(terms)
        }
        AtomView::Pow(power) => {
            let Ok(exponent) = i64::try_from(power.get_exp()) else {
                return Ok(vec![expression.to_owned()]);
            };
            let Ok(mut exponent) = usize::try_from(exponent) else {
                return Ok(vec![expression.to_owned()]);
            };
            let mut factor = bounded_expand_summands(power.get_base(), limit)?;
            if exponent == 0 {
                return Ok(vec![Atom::one()]);
            }
            if factor.len() == 1 {
                return Ok(vec![
                    factor
                        .pop()
                        .expect("one expanded base term")
                        .pow(power.get_exp().to_owned()),
                ]);
            }

            let mut terms = vec![Atom::one()];
            while exponent != 0 {
                if exponent & 1 == 1 {
                    terms = bounded_term_product(&terms, &factor, limit)?;
                }
                exponent >>= 1;
                if exponent != 0 {
                    factor = bounded_term_product(&factor, &factor, limit)?;
                }
            }
            Ok(terms)
        }
        _ => Ok(vec![expression.to_owned()]),
    }
}

fn bounded_term_product(
    left: &[Atom],
    right: &[Atom],
    limit: usize,
) -> Result<Vec<Atom>, TensorReductionError> {
    let term_count =
        left.len()
            .checked_mul(right.len())
            .ok_or(TensorReductionError::ExpansionLimit {
                terms: usize::MAX,
                limit,
            })?;
    if term_count > limit {
        return Err(TensorReductionError::ExpansionLimit {
            terms: term_count,
            limit,
        });
    }
    let mut terms = Vec::with_capacity(term_count);
    for left in left {
        for right in right {
            terms.push(left * right);
        }
    }
    Ok(terms)
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
    if let AtomView::Fun(function) = expression
        && function.get_symbol().get_namespace() == "spenso"
        && function.get_symbol().get_stripped_name() == "mink"
    {
        let arguments = function.iter().collect::<Vec<_>>();
        if arguments.as_slice() == [dimension.as_view(), index.as_view()] {
            return Ok(multiplicity);
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

fn has_dangling_minkowski_indices(
    expression: &Atom,
    dimension: &Atom,
) -> Result<bool, TensorReductionError> {
    fn collect(
        expression: AtomView<'_>,
        dimension: &Atom,
        multiplicity: usize,
        counts: &mut BTreeMap<Atom, usize>,
    ) -> Result<(), TensorReductionError> {
        if let AtomView::Fun(function) = expression
            && function.get_symbol().get_namespace() == "spenso"
            && function.get_symbol().get_stripped_name() == "mink"
            && function.get_nargs() == 2
        {
            let arguments = function.iter().collect::<Vec<_>>();
            if arguments[0] != dimension.as_view() {
                return Err(TensorReductionError::DimensionMismatch {
                    expected: dimension.clone(),
                    found: arguments[0].to_owned(),
                });
            }
            let count = counts.entry(arguments[1].to_owned()).or_default();
            *count = count.saturating_add(multiplicity);
            return Ok(());
        }
        if let AtomView::Pow(power) = expression
            && let Ok(exponent) = i64::try_from(power.get_exp())
            && let Ok(exponent) = usize::try_from(exponent)
        {
            return collect(
                power.get_base(),
                dimension,
                multiplicity.saturating_mul(exponent),
                counts,
            );
        }
        for child in expression.children() {
            collect(child, dimension, multiplicity, counts)?;
        }
        Ok(())
    }

    let mut counts = BTreeMap::new();
    collect(expression.as_view(), dimension, 1, &mut counts)?;
    for (index, occurrences) in &counts {
        if *occurrences > 2 {
            return Err(TensorReductionError::AmbiguousMinkowskiIndex {
                index: index.clone(),
                occurrences: *occurrences,
            });
        }
    }
    Ok(counts.values().any(|occurrences| *occurrences == 1))
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
    Ok(ContractionOrbit {
        vectors: vec![vector].into_boxed_slice(),
        multiplicities: vec![2 * pair_count].into_boxed_slice(),
        alpha: vec![
            u16::try_from(pair_count).map_err(|_| TensorReductionError::OutputLimit {
                terms: pair_count,
                limit: usize::from(u16::MAX),
            })?,
        ]
        .into_boxed_slice(),
        labeled_pairings: perfect_matching_count(2 * pair_count),
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
    limit: usize,
) -> Result<Vec<ContractionOrbit>, TensorReductionError> {
    let (vectors, multiplicities, _) = classified_colors(colors);
    let mut cache = HashMap::new();
    let multiplicities_u16 = multiplicities
        .iter()
        .map(|count| {
            u16::try_from(*count).map_err(|_| TensorReductionError::UnsupportedRank {
                rank: colors.len(),
                maximum: usize::from(u16::MAX),
            })
        })
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
    limit: usize,
) -> Result<BTreeMap<Vec<u16>, u128>, TensorReductionError> {
    if let Some(cached) = cache.get(counts) {
        return Ok(cached.clone());
    }
    let alpha_len = counts.len() * (counts.len() + 1) / 2;
    let Some(left) = counts.iter().position(|count| *count != 0) else {
        return Ok(BTreeMap::from([(vec![0; alpha_len], 1)]));
    };

    let mut after_left = counts.to_vec();
    after_left[left] -= 1;
    let mut result = BTreeMap::<Vec<u16>, u128>::new();
    for right in left..counts.len() {
        let choices = after_left[right];
        if choices == 0 {
            continue;
        }
        let mut remaining = after_left.clone();
        remaining[right] -= 1;
        for (mut alpha, multiplicity) in orbit_alpha(&remaining, cache, limit)? {
            alpha[upper_triangle_index(left, right, counts.len())] += 1;
            *result.entry(alpha).or_default() += multiplicity * u128::from(choices);
            if result.len() > limit {
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

fn perfect_matching_count(rank: usize) -> u128 {
    (1..rank).step_by(2).map(|factor| factor as u128).product()
}

fn perfect_matchings(rank: usize) -> Vec<Vec<usize>> {
    fn generate(mask: u32, partners: &mut [usize], output: &mut Vec<Vec<usize>>) {
        if mask == 0 {
            output.push(partners.to_vec());
            return;
        }
        let left = mask.trailing_zeros() as usize;
        let without_left = mask & !(1 << left);
        let mut candidates = without_left;
        while candidates != 0 {
            let right = candidates.trailing_zeros() as usize;
            partners[left] = right;
            partners[right] = left;
            generate(without_left & !(1 << right), partners, output);
            candidates &= !(1 << right);
        }
    }

    let mut output = Vec::with_capacity(perfect_matching_count(rank) as usize);
    generate((1_u32 << rank) - 1, &mut vec![0; rank], &mut output);
    output
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
            .collect();
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
    pairing_limit: u128,
    output_limit: usize,
) -> Result<BTreeMap<Atom, u128>, TensorReductionError> {
    let pairings = perfect_matching_count(outside.len());
    if pairings > pairing_limit {
        return Err(TensorReductionError::PairingLimit {
            rank: outside.len(),
            pairings,
            limit: pairing_limit,
        });
    }
    let mut atoms = BTreeMap::<Atom, u128>::new();
    for pairing in perfect_matchings(outside.len()) {
        *atoms
            .entry(outside_pairing_atom(outside, &pairing))
            .or_default() += 1;
        if atoms.len() > output_limit {
            return Err(TensorReductionError::OutputLimit {
                terms: atoms.len(),
                limit: output_limit,
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

fn coset_key_from_parts(parts: &[usize]) -> u64 {
    parts
        .iter()
        .fold(parts.len() as u64, |key, part| (key << 5) | *part as u64)
}

fn matching_coset_key(left: &[usize], right: &[usize]) -> u64 {
    let mut seen = 0_u32;
    let all = (1_u32 << left.len()) - 1;
    let mut part_multiplicities = [0_u8; OrthogonalWeingarten::MAX_RANK / 2 + 1];
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
        part_multiplicities[size] += 1;
        part_count += 1;
    }
    let mut key = part_count;
    for size in (1..part_multiplicities.len()).rev() {
        for _ in 0..part_multiplicities[size] {
            key = (key << 5) | size as u64;
        }
    }
    key
}

fn atom_from_count(count: u128) -> Result<Atom, TensorReductionError> {
    let count =
        i64::try_from(count).map_err(|_| TensorReductionError::MultiplicityOverflow(count))?;
    Ok(Atom::num(count))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CosetType;
    use feynkit_graph::{DiagramEdge, DiagramVertex};
    use feynkit_model::Model;
    use spenso::vector_symbol;
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
        left.numerator = numerator.clone();
        let mut builder = FeynmanDiagram::builder(model, name)
            .numerator(numerator)
            .numerator_prefactor(numerator_prefactor)
            .projector(projector);
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
        let pairings = perfect_matchings(integrated.len());
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
            .reduce(Atom::one().as_view())
            .unwrap_err();
        assert!(matches!(
            error,
            TensorReductionError::NoIntegratedVectorsSelected
        ));
    }

    #[test]
    fn allocation_free_coset_keys_match_the_public_definition() {
        let pairings = perfect_matchings(6);
        for left in &pairings {
            for right in &pairings {
                let coset = CosetType::between_matchings(left, right).unwrap();
                assert_eq!(
                    matching_coset_key(left, right),
                    coset_key_from_parts(coset.parts())
                );
            }
        }
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
            .collect::<Vec<_>>();
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
        assert!((reduction.expression() - &expected).together().is_zero());
        let graphs = diagram.reduce_tensor_graphs(&reducer).unwrap();
        assert_eq!(graphs.len(), 1);
        let graph = &graphs[0];
        assert_eq!(graph.projector(), &Atom::one());
        assert_eq!(graph.numerator_prefactor(), &prefactor);
        assert!((graph.numerator() - expected).together().is_zero());
        assert_eq!(diagram.projector(), &projector);
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
            .reduce(input.as_view())
            .unwrap();
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 2, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(r, 1, &dimension))
            / dimension;
        assert_eq!(result.expression().expand(), expected.expand());
        assert!(result.is_fully_contracted());
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
            .reduce(input.as_view())
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
            .reduce(input.as_view())
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
        let difference = (result.expression() - &expected).together();
        assert!(
            difference.is_zero(),
            "result: {}\nexpected: {expected}\ndifference: {difference}",
            result.expression()
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
            .reduce(input.as_view())
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
        assert!((result.expression() - expected).together().is_zero());
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
            .reduce(input.as_view())
            .unwrap();
        let elapsed = started.elapsed();

        // The default labeled budget is only 20,000, whereas rank 20 has
        // 654,729,075 pairings.  Reaching this compact four-term result proves
        // that both partial symmetries used the orbit projector.
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
            * OrthogonalWeingarten::isotropic_pairing_coefficient(10, dimension);
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
            .reduce(isotropic_internal_input.as_view())
            .unwrap()
            .expression();
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
            .reduce(isotropic_projector_input.as_view())
            .unwrap()
            .expression();
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
            .reduce(input.as_view())
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
            .reduce(input.as_view())
            .unwrap();
        let materialized = result.expression();
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
            .reduce(input.as_view())
            .unwrap();
        assert_eq!(
            result.expression(),
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
            .reduce(input.as_view())
            .unwrap();

        assert_eq!(
            result.expression(),
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

        let result = reducer.reduce(input.as_view()).unwrap();
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
        assert!(matches!(
            TensorReducer::new(dimension)
                .with_integrated_head(k)
                .reduce(input.as_view()),
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
        let opaque = FunctionBuilder::new(symbol!("feynkit_tensor_test::opaque"))
            .add_arg(indexed(p, 1, &dimension, &mu))
            .finish();
        let diagram = diagram_with_numerator("opaque-free", opaque.clone());
        let reducer = TensorReducer::new(dimension).with_integrated_head(k);

        let result = reducer.reduce(opaque.as_view()).unwrap();
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
        let opaque = FunctionBuilder::new(symbol!("feynkit_tensor_test::opaque_consumer"))
            .add_arg(indexed(p, 1, &dimension, &mu))
            .finish();
        let input = indexed(k, 1, &dimension, &mu) * indexed(k, 2, &dimension, &mu) * opaque;
        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .reduce(input.as_view());
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
        let tensor = symbol!("feynkit_tensor_test::opaque_projector_tensor");
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
            .reduce(input.as_view())
            .unwrap();
        let difference = (result.expression() - expected).together();
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
        let tensor = FunctionBuilder::new(symbol!("feynkit_tensor_test::opaque_power_tensor"))
            .add_arg(minkowski_slot(&dimension, &spin))
            .add_arg(minkowski_slot(&dimension, &mu))
            .finish();
        let input = indexed(k, 1, &dimension, &mu) * tensor.pow(Atom::num(2));

        let result = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .reduce(input.as_view());
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
        let tensor =
            FunctionBuilder::new(symbol!("feynkit_tensor_test::opaque_symbolic_power_tensor"))
                .add_arg(minkowski_slot(&dimension, &spin))
                .add_arg(minkowski_slot(&dimension, &mu))
                .finish();
        let input = indexed(k, 1, &dimension, &mu) * tensor.pow(exponent);

        assert!(matches!(
            TensorReducer::new(dimension)
                .with_integrated_head(k)
                .reduce(input.as_view()),
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
            Err(TensorReductionError::UnsupportedRank {
                rank: 1_000_000,
                maximum: 20
            })
        ));
    }

    #[test]
    fn factored_sum_expansion_obeys_the_output_budget_while_distributing() {
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
        let error = TensorReducer::new(dimension)
            .with_integrated_head(k)
            .with_output_term_limit(7)
            .reduce(input.as_view())
            .unwrap_err();
        assert!(matches!(
            error,
            TensorReductionError::ExpansionLimit { terms: 8, limit: 7 }
        ));
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
            .reduce(input.as_view())
            .unwrap();
        assert_eq!(
            result.expression(),
            dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
        );
        assert!(result.is_fully_contracted());
    }

    #[test]
    fn products_of_generated_vertex_sums_are_expanded_before_metric_sewing() {
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
            .reduce(input.as_view())
            .unwrap()
            .expression();
        let expected = Atom::num(4)
            * dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * dot(&compact(p, 1, &dimension), &compact(p, 1, &dimension))
            / dimension;
        assert!((result - expected).together().is_zero());
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
            .reduce(input.as_view())
            .unwrap();
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * metric(&dimension, &left_free, &right_free)
            / dimension.clone();
        let expression = result.expression().expand();
        assert_eq!(expression, expected.expand());
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
            .reduce(input.as_view())
            .unwrap();
        let expected = dot(&compact(k, 1, &dimension), &compact(k, 1, &dimension))
            * FunctionBuilder::new(ETS.metric)
                .add_arg(minkowski_slot(&dimension, &mu))
                .add_arg(minkowski_slot(&dimension, &nu))
                .finish()
            / dimension.clone();
        let expression = result.expression();
        assert_eq!(expression.expand(), expected.expand());
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
        let orbits = contraction_orbits(&colors, 100).unwrap();
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
