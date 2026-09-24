//! Symbolica-validated descriptions of graph-aware sampling maps.
//!
//! The graph-independent map language and the reusable ordinary-map numerical
//! kernel live here.  Graph-aware surface and cut maps are compiled by the
//! process layer. Parsing is structural: a complete Symbolica expression is
//! checked recursively, so a valid nested call cannot hide an invalid root or
//! sibling.

use std::sync::{Arc, Mutex};

use super::sampling_context::{PreparedSurfaceStatus, SamplingMapContext};
use super::sampling_evaluator::{SamplingDualValue, SamplingExpressionEvaluator};
use crate::momentum::ThreeMomentum;
use crate::momentum::sample::LoopMomenta;
use crate::settings::runtime::{
    ParameterizationMapping, ParameterizationMode, ParameterizationSettings, SamplingRadialProfile,
};
use crate::utils::newton_solver::safeguarded_newton_iteration_and_derivative;
use crate::utils::{F, FloatLike, global_inv_parameterize, global_parameterize};
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomView},
    symbol, try_parse,
};

/// The geometric role of a resolved sampling map.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum SamplingMapDefinition {
    /// A loop-momentum basis channel, identified by its ordered edge list.
    Lmb(Vec<usize>),
    /// An energy surface in the supplied edge subspace.
    Surface(Vec<usize>),
    /// A physical Cutkosky cut.
    Cut(Vec<usize>),
    /// A soft channel around one routed edge.
    Soft(usize),
    /// A qualified collinear channel between two routed edges.
    Collinear(usize, usize),
    /// Coordinates left for a complement map.
    Complement(Vec<usize>),
    /// Independent maps on a direct-sum block.
    Product(Vec<Self>),
    /// Joint coordinates for compatible constraints.
    Intersect(Vec<Self>),
    /// Ordered conditional maps.  Later maps may depend on earlier data.
    Then(Vec<Self>),
    /// Prepare physical cut phase-space data before mapping the remaining maps.
    PhaseSpace(Box<Self>),
    /// Resolve a target on the left side of a prepared cut.
    Left(Box<Self>),
    /// Resolve a target on the right side of a prepared cut.
    Right(Box<Self>),
    /// Resolve a target in the kinematics of an explicit physical host cut.
    AtCut { cut: Vec<usize>, map: Box<Self> },
    /// Assign coordinates to a child; its LMB descriptor consumes no dimensions.
    Block { edges: Vec<usize>, map: Box<Self> },
}

impl SamplingMapDefinition {
    /// Parse and structurally validate one complete Symbolica map expression.
    pub fn parse(source: &str) -> Result<Self> {
        let atom = try_parse!(source.trim())
            .map_err(|error| eyre!("failed to parse sampling map `{source}`: {error}"))?;
        Self::from_atom(atom.as_view())
    }

    /// Build a map from an already parsed Symbolica atom.
    pub fn from_atom(atom: AtomView<'_>) -> Result<Self> {
        let AtomView::Fun(function) = atom else {
            return Err(eyre!(
                "sampling map must be a constructor call, got `{atom}`"
            ));
        };
        let symbol = function.get_symbol();
        let name = symbol.get_name();
        let arguments = function.iter().collect::<Vec<_>>();
        let short_name = name.rsplit("::").next().unwrap_or(name);
        // `try_parse!` resolves unqualified constructors in the current
        // crate namespace (`gammalooprs::surface`), while canonical atoms
        // emitted by `to_atom` use the dedicated sampling-map namespace.
        // Accept both forms but reject unrelated user namespaces.
        if name != short_name
            && name != format!("gammalooprs::{short_name}")
            && name != format!("gammalooprs::sampling_map::{short_name}")
        {
            return Err(eyre!(
                "sampling-map constructor `{name}` is outside the reserved sampling-map namespace"
            ));
        }

        match short_name {
            "lmb" => Self::ordered_edges(arguments, "lmb").map(Self::Lmb),
            "surface" => Self::edges(arguments, "surface").map(Self::Surface),
            "cut" => Self::edges(arguments, "cut").map(Self::Cut),
            "complement" => Self::edges(arguments, "complement").map(Self::Complement),
            "soft" => Self::one_edge(arguments, "soft").map(Self::Soft),
            "collinear" => {
                Self::two_edges(arguments, "collinear").map(|(a, b)| Self::Collinear(a, b))
            }
            "product" => Self::maps(arguments, "product").map(Self::Product),
            "intersect" => Self::maps(arguments, "intersect").map(Self::Intersect),
            "then" => Self::maps(arguments, "then").map(Self::Then),
            "phase_space" => {
                Self::one_map(arguments, "phase_space").map(|map| Self::PhaseSpace(Box::new(map)))
            }
            "left" => Self::one_map(arguments, "left").map(|map| Self::Left(Box::new(map))),
            "right" => Self::one_map(arguments, "right").map(|map| Self::Right(Box::new(map))),
            "at_cut" | "block" => {
                let [descriptor, map] = arguments.as_slice() else {
                    return Err(eyre!("{short_name} expects exactly a descriptor and a map"));
                };
                let descriptor = Self::from_atom(*descriptor)?;
                let map = Box::new(Self::from_atom(*map)?);
                match (short_name, descriptor) {
                    ("at_cut", Self::Cut(cut)) => Ok(Self::AtCut { cut, map }),
                    ("block", Self::Lmb(edges)) => Ok(Self::Block { edges, map }),
                    _ => Err(eyre!(
                        "{short_name} requires {} as its first argument",
                        if short_name == "at_cut" {
                            "cut(...)"
                        } else {
                            "lmb(...)"
                        }
                    )),
                }
            }
            _ => Err(eyre!(
                "unknown sampling-map constructor `{name}`; expected lmb, surface, cut, soft, collinear, complement, product, intersect, then, phase_space, left, right, at_cut or block"
            )),
        }
    }

    /// Return a canonical Symbolica expression for this definition.
    pub fn to_atom(&self) -> Atom {
        match self {
            Self::Lmb(edges) => call(
                "lmb",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Surface(edges) => call(
                "surface",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Cut(edges) => call(
                "cut",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Soft(edge) => call("soft", [Atom::num(*edge as i64)]),
            Self::Collinear(a, b) => {
                call("collinear", [Atom::num(*a as i64), Atom::num(*b as i64)])
            }
            Self::Complement(edges) => call(
                "complement",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Product(maps) => call("product", maps.iter().map(Self::to_atom)),
            Self::Intersect(maps) => call("intersect", maps.iter().map(Self::to_atom)),
            Self::Then(maps) => call("then", maps.iter().map(Self::to_atom)),
            Self::PhaseSpace(map) => call("phase_space", [map.to_atom()]),
            Self::Left(map) => call("left", [map.to_atom()]),
            Self::Right(map) => call("right", [map.to_atom()]),
            Self::AtCut { cut, map } => {
                call("at_cut", [Self::Cut(cut.clone()).to_atom(), map.to_atom()])
            }
            Self::Block { edges, map } => {
                call("block", [Self::Lmb(edges.clone()).to_atom(), map.to_atom()])
            }
        }
    }

    /// Physical equations are distinct from the active coordinate block. A
    /// supported intersection keeps its two energy sets in normal-coordinate
    /// order; their union would lose the equations and their shared energy.
    pub fn energy_edge_sets(&self) -> Vec<&[usize]> {
        match self {
            Self::Surface(edges) | Self::Cut(edges) => vec![edges],
            Self::Intersect(maps)
                if maps.len() == 2
                    && maps[0] != maps[1]
                    && maps.iter().all(|map| matches!(map, Self::Surface(_))) =>
            {
                maps.iter().flat_map(Self::energy_edge_sets).collect()
            }
            Self::PhaseSpace(map)
            | Self::Left(map)
            | Self::Right(map)
            | Self::AtCut { map, .. }
            | Self::Block { map, .. } => map.energy_edge_sets(),
            _ => Vec::new(),
        }
    }

    pub fn host_cut(&self) -> Option<&[usize]> {
        match self {
            Self::AtCut { cut, .. } => Some(cut),
            Self::PhaseSpace(map) => match map.as_ref() {
                Self::Cut(edges) => Some(edges),
                _ => None,
            },
            Self::Left(map) | Self::Right(map) | Self::Block { map, .. } => map.host_cut(),
            _ => None,
        }
    }

    pub fn cut_side(&self) -> Option<super::SamplingCutSide> {
        match self {
            Self::Left(_) => Some(super::SamplingCutSide::Left),
            Self::Right(_) => Some(super::SamplingCutSide::Right),
            Self::AtCut { map, .. } | Self::Block { map, .. } => map.cut_side(),
            _ => None,
        }
    }

    pub fn is_phase_space(&self) -> bool {
        match self {
            Self::PhaseSpace(map) => matches!(map.as_ref(), Self::Cut(_)),
            Self::AtCut { map, .. } | Self::Block { map, .. } => map.is_phase_space(),
            _ => false,
        }
    }

    fn edges(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Vec<usize>> {
        if arguments.is_empty() {
            return Err(eyre!("{constructor} expects at least one edge"));
        }
        let mut edges = arguments
            .into_iter()
            .map(|argument| parse_edge(argument, constructor))
            .collect::<Result<Vec<_>>>()?;
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "{constructor} contains a duplicate edge in {edges:?}"
            ));
        }
        Ok(edges)
    }

    fn ordered_edges(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Vec<usize>> {
        if arguments.is_empty() {
            return Err(eyre!("{constructor} expects at least one edge"));
        }
        let edges = arguments
            .into_iter()
            .map(|argument| parse_edge(argument, constructor))
            .collect::<Result<Vec<_>>>()?;
        let mut sorted_edges = edges.clone();
        sorted_edges.sort_unstable();
        if sorted_edges.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "{constructor} contains a duplicate edge in {edges:?}"
            ));
        }
        Ok(edges)
    }

    fn one_edge(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<usize> {
        match arguments.as_slice() {
            [edge] => parse_edge(*edge, constructor),
            _ => Err(eyre!(
                "{constructor} expects exactly one edge, got {} arguments",
                arguments.len()
            )),
        }
    }

    fn two_edges(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<(usize, usize)> {
        match arguments.as_slice() {
            [a, b] => {
                let a = parse_edge(*a, constructor)?;
                let b = parse_edge(*b, constructor)?;
                if a == b {
                    return Err(eyre!("{constructor} requires two distinct edges, got {a}"));
                }
                Ok((a, b))
            }
            _ => Err(eyre!(
                "{constructor} expects exactly two edges, got {} arguments",
                arguments.len()
            )),
        }
    }

    fn maps(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Vec<Self>> {
        if arguments.is_empty() {
            return Err(eyre!("{constructor} expects at least one map"));
        }
        arguments.into_iter().map(Self::from_atom).collect()
    }

    fn one_map(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Self> {
        match arguments.as_slice() {
            [map] => Self::from_atom(*map),
            _ => Err(eyre!(
                "{constructor} expects exactly one map, got {} arguments",
                arguments.len()
            )),
        }
    }
}

/// How a numerical kernel establishes the domain of its push-forward density.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingSupport {
    /// The map covers the full raw domain, including absent-surface fallback data.
    Full,
    /// The map covers a certified subset and needs a full-support sibling.
    Restricted,
}

/// Contract for the Jacobian supplied by a map kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingJacobian {
    /// Exact determinant of the selected forward map.
    ExactForward,
    /// Exact determinant obtained by implicit differentiation of a solved root.
    ExactImplicit,
    /// A positive partition score only; never the selected channel's integration Jacobian.
    ProxyOnly,
}

/// Static domain, prerequisite, and Jacobian contract of a compiled map kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SamplingMapContract {
    pub support: SamplingSupport,
    /// An earlier conditional map must prepare data before this map is valid.
    /// Supplying those prerequisites does not enlarge a restricted domain.
    pub requires_context: bool,
    /// The original-draw phase must freeze this map's proposal decisions.
    pub requires_proposal_policy: bool,
    pub jacobian: SamplingJacobian,
}

/// A numerical kernel for the ordinary global loop-momentum maps.
///
/// This is deliberately a small, graph-independent layer.  LMB selection and
/// reinterpretation stay in the process sampler; this kernel owns only the
/// push-forward map on the selected loop coordinates and its exact determinant.
/// Keeping the forward and inverse calls together also gives acceptance tests a
/// residual without having to duplicate the global parameterisation logic.
#[derive(Clone, Debug)]
pub struct SamplingMapKernel {
    definition: SamplingMapDefinition,
    settings: ParameterizationSettings,
    e_cm: f64,
    dimensions: usize,
}

/// Values returned by one forward or inverse map evaluation.
#[derive(Clone, Debug)]
pub struct SamplingMapPoint<T: FloatLike> {
    pub coordinates: Vec<F<T>>,
    pub loop_momenta: LoopMomenta<F<T>>,
    /// The positive exact forward determinant `|d p / d x|`.
    pub jacobian: F<T>,
    /// The corresponding exact inverse determinant.
    pub inverse_jacobian: F<T>,
    /// Infinity norm of the independent round-trip residual.
    pub residual: F<T>,
}

/// Numerical sampling failures which may be retried from the original source
/// at higher precision. Structural metadata and capability errors deliberately
/// keep their existing error types and must not enter the rescue stack.
#[derive(Clone, Debug)]
pub enum SamplingEvaluationError {
    Unrepresentable {
        operation: &'static str,
        detail: String,
    },
    UncertifiedRoot {
        detail: String,
    },
    UncertifiedOverlap {
        detail: String,
    },
    UncertainGeometry {
        detail: String,
    },
}

impl std::fmt::Display for SamplingEvaluationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unrepresentable { operation, detail } => write!(
                formatter,
                "sampling {operation} is not representable at the current precision: {detail}"
            ),
            Self::UncertifiedRoot { detail } => {
                write!(
                    formatter,
                    "sampling root is not certified at the current precision: {detail}"
                )
            }
            Self::UncertifiedOverlap { detail } => {
                write!(formatter, "threshold overlap is not certified: {detail}")
            }
            Self::UncertainGeometry { detail } => {
                write!(
                    formatter,
                    "sampling geometry is uncertain at the current precision: {detail}"
                )
            }
        }
    }
}

impl std::error::Error for SamplingEvaluationError {}

/// Result of evaluating one graph-independent sampling component.
///
/// The point is represented as a flat vector in the master raw coordinate
/// frame.  A composition therefore never has to reinterpret a child's local
/// coordinates.  `diagnostics` carries branch/support information to callers;
/// it is deliberately data rather than logging so acceptance tests can audit
/// every branch taken by a map.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingMapEvaluation<T: FloatLike = f64> {
    pub coordinates: Vec<T>,
    pub point: Vec<T>,
    pub jacobian: T,
    pub inverse_jacobian: T,
    pub residual: T,
    pub support: SamplingSupport,
    pub diagnostics: Vec<String>,
}

impl<T: FloatLike> SamplingMapEvaluation<T> {
    /// Preserve the composed Jacobian-pair checks when only its inverse density
    /// is requested. Reciprocal overflow is still a numerical failure.
    pub(crate) fn validate_inverse_density(density: T) -> Result<T> {
        let inverse = F(density);
        let forward = inverse.inv();
        if [&inverse, &forward]
            .iter()
            .any(|value| !value.0.is_finite() || **value <= value.zero())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "composed inverse density",
                detail: "nonfinite or nonpositive Jacobian pair".into(),
            }
            .into());
        }
        Ok(inverse.0)
    }

    /// Certify the complete numerical result after composition or affine
    /// routing, where individually finite child determinants may overflow.
    pub(crate) fn validate(self, operation: &'static str) -> Result<Self> {
        if self
            .coordinates
            .iter()
            .chain(&self.point)
            .any(|value| !value.is_finite())
            || !self.residual.is_finite()
            || self.residual < self.residual.zero()
            || [&self.jacobian, &self.inverse_jacobian]
                .iter()
                .any(|value| !value.is_finite() || **value <= value.zero())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation,
                detail: format!(
                    "non-finite coordinates/residual or nonpositive Jacobian pair ({}, {})",
                    self.jacobian, self.inverse_jacobian
                ),
            }
            .into());
        }
        Ok(self)
    }
}

/// A graph-independent map component which can participate in a product or
/// ordered conditional composition. The borrowed context carries previously
/// produced raw coordinates for `then` maps (empty for independent products)
/// together with the original-draw proposal access and stable compiled path.
pub trait SamplingMapComponent<T: FloatLike = f64>:
    std::fmt::Debug + Send + Sync + dyn_clone::DynClone
{
    fn dimensions(&self) -> usize;
    fn output_dimensions(&self) -> usize;
    fn contract(&self) -> SamplingMapContract;
    fn name(&self) -> &'static str;
    fn forward(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>>;
    /// Return None only for certified exclusion from the map's support. An
    /// uncertain boundary, non-finite density, or failed solve remains an error.
    /// Regular branches must have a unique inverse under the map's coordinate
    /// convention; overlapping inverse branches require complete density
    /// accounting before they can enter this component interface.
    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>>;

    /// Exact supplied-point density and support, without requiring a diagnostic
    /// forward reconstruction. Components with cheap inverses reuse them; costly
    /// radial inverses override this while retaining the same root and density.
    fn inverse_density(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<T>> {
        self.inverse(point, context)
            .map(|evaluation| evaluation.map(|evaluation| evaluation.inverse_jacobian))
    }
}

dyn_clone::clone_trait_object!(<T> SamplingMapComponent<T> where T: FloatLike);

/// An exact affine change of coordinates between two equal-dimensional frames.
///
/// The matrix is stored row-major and the map is `point = matrix * coordinates
/// + translation`.  This is the graph-independent numerical owner for routing
///   a selected LMB into a parent frame; graph code supplies the matrix and
///   translation after resolving the relevant edge signatures and external data.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingMapAffine<T: FloatLike = f64> {
    matrix: Vec<T>,
    inverse_matrix: Vec<T>,
    translation: Vec<T>,
    dimension: usize,
    determinant: T,
}

impl<T: FloatLike> SamplingMapAffine<T> {
    /// Construct an affine map from a square matrix and translation vector.
    ///
    /// Partial pivoting is used while constructing the inverse.  Singular and
    /// numerically rank-deficient matrices are rejected before a map can enter
    /// a sampling composition, since they have no valid push-forward density.
    pub fn new(matrix: Vec<Vec<T>>, translation: Vec<T>) -> Result<Self> {
        let dimension = matrix.len();
        if dimension == 0 {
            return Err(eyre!("affine sampling map requires a non-empty matrix"));
        }
        if translation.len() != dimension {
            return Err(eyre!(
                "affine sampling map translation has dimension {}, expected {}",
                translation.len(),
                dimension
            ));
        }
        if translation.iter().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "affine sampling map translation must contain only finite values"
            ));
        }
        if matrix.iter().any(|row| row.len() != dimension) {
            return Err(eyre!(
                "affine sampling map matrix must be square ({} rows of length {})",
                dimension,
                dimension
            ));
        }
        if matrix.iter().flatten().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "affine sampling map matrix must contain only finite values"
            ));
        }

        let zero = F(matrix[0][0].zero());
        let one = zero.one();
        let mut augmented = vec![vec![zero.clone(); 2 * dimension]; dimension];
        let mut scale = zero.clone();
        for (row_index, row) in matrix.iter().enumerate() {
            for (column_index, value) in row.iter().cloned().enumerate() {
                augmented[row_index][column_index] = F(value.clone());
                scale = scale.max(F(value).abs());
            }
            augmented[row_index][dimension + row_index] = one.clone();
        }
        if scale == zero {
            return Err(eyre!("affine sampling map matrix is singular"));
        }
        let pivot_tolerance = &scale * scale.epsilon() * scale.from_usize(dimension * 32);
        let mut determinant = one.clone();
        let mut row_sign = one.clone();
        for pivot_column in 0..dimension {
            let (pivot_row, pivot_abs) = (pivot_column..dimension)
                .map(|row| (row, augmented[row][pivot_column].abs()))
                .max_by(|(_, left), (_, right)| {
                    left.partial_cmp(right).unwrap_or(std::cmp::Ordering::Equal)
                })
                .expect("non-empty pivot range");
            if pivot_abs <= pivot_tolerance || !pivot_abs.0.is_finite() {
                return Err(eyre!(
                    "affine sampling map matrix is singular or rank-deficient at pivot {} (|pivot| = {}, tolerance = {})",
                    pivot_column,
                    pivot_abs,
                    pivot_tolerance
                ));
            }
            if pivot_row != pivot_column {
                augmented.swap(pivot_row, pivot_column);
                row_sign = -row_sign;
            }
            let pivot = augmented[pivot_column][pivot_column].clone();
            determinant *= &pivot;
            let inverse_pivot = &one / pivot;
            for value in &mut augmented[pivot_column] {
                *value *= &inverse_pivot;
            }
            // Keep a private copy so that the pivot row can be read while
            // each other row is updated through the mutable matrix.
            let pivot_values = augmented[pivot_column].clone();
            for (row, augmented_row) in augmented.iter_mut().enumerate().take(dimension) {
                if row == pivot_column {
                    continue;
                }
                let factor = augmented_row[pivot_column].clone();
                if factor == zero {
                    continue;
                }
                for (column, pivot_value) in pivot_values.iter().enumerate().take(2 * dimension) {
                    let term = &factor * pivot_value;
                    augmented_row[column] -= term;
                }
            }
        }
        determinant *= row_sign;
        if !determinant.0.is_finite() || determinant == zero {
            return Err(eyre!(
                "affine sampling map determinant is not finite and non-zero: {determinant}"
            ));
        }
        let inverse_matrix = augmented
            .into_iter()
            .flat_map(|row| row.into_iter().skip(dimension).map(|value| value.0))
            .collect::<Vec<_>>();
        let determinant = determinant.abs();
        if !determinant.0.is_finite() || determinant <= zero {
            return Err(eyre!(
                "affine sampling map Jacobian is not finite and positive: {determinant}"
            ));
        }
        Ok(Self {
            matrix: matrix.into_iter().flatten().collect(),
            inverse_matrix,
            translation,
            dimension,
            determinant: determinant.0,
        })
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn determinant(&self) -> T {
        self.determinant.clone()
    }

    pub fn matrix(&self) -> &[T] {
        &self.matrix
    }

    pub fn inverse_matrix(&self) -> &[T] {
        &self.inverse_matrix
    }

    pub fn translation(&self) -> &[T] {
        &self.translation
    }

    /// Invert the complete affine frame with its exact determinant. The shared
    /// component adapter adds support metadata without another numerical path.
    pub fn inverse(&self, point: &[T], _context: &[T]) -> Result<SamplingMapEvaluation<T>> {
        self.validate_dimension(point, "point")?;
        let translated = point
            .iter()
            .zip(&self.translation)
            .map(|(value, shift)| (F(value.clone()) - F(shift.clone())).0)
            .collect::<Vec<_>>();
        let coordinates = self.apply(
            &self.inverse_matrix,
            &translated,
            &vec![self.determinant.zero(); self.dimension],
        );
        let recovered = self.apply(&self.matrix, &coordinates, &self.translation);
        SamplingMapEvaluation {
            coordinates,
            point: point.to_vec(),
            jacobian: self.determinant.clone(),
            inverse_jacobian: (F(self.determinant.one()) / F(self.determinant.clone())).0,
            residual: max_coordinate_residual_scalar(&recovered, point),
            support: SamplingSupport::Full,
            diagnostics: vec![format!("affine dimension={}", self.dimension)],
        }
        .validate("affine map")
    }

    fn apply(&self, matrix: &[T], input: &[T], translation: &[T]) -> Vec<T> {
        (0..self.dimension)
            .map(|row| {
                matrix[row * self.dimension..(row + 1) * self.dimension]
                    .iter()
                    .zip(input)
                    .fold(F(translation[row].clone()), |sum, (coefficient, value)| {
                        sum + F(coefficient.clone()) * F(value.clone())
                    })
                    .0
            })
            .collect()
    }

    fn validate_dimension(&self, values: &[T], role: &str) -> Result<()> {
        if values.len() != self.dimension {
            return Err(eyre!(
                "affine sampling map {role} has dimension {}, expected {}",
                values.len(),
                self.dimension
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "affine sampling map {role} must contain only finite values"
            ));
        }
        Ok(())
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for SamplingMapAffine<T> {
    fn dimensions(&self) -> usize {
        self.dimension
    }

    fn output_dimensions(&self) -> usize {
        self.dimension
    }

    fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: SamplingSupport::Full,
            requires_context: false,
            requires_proposal_policy: false,
            jacobian: SamplingJacobian::ExactForward,
        }
    }

    fn name(&self) -> &'static str {
        "affine"
    }

    fn forward(
        &self,
        coordinates: &[T],
        _context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        self.validate_dimension(coordinates, "input")?;
        let point = self.apply(&self.matrix, coordinates, &self.translation);
        let translated = point
            .iter()
            .zip(&self.translation)
            .map(|(value, shift)| (F(value.clone()) - F(shift.clone())).0)
            .collect::<Vec<_>>();
        let recovered = self.apply(
            &self.inverse_matrix,
            &translated,
            &vec![self.determinant.zero(); self.dimension],
        );
        SamplingMapEvaluation {
            coordinates: coordinates.to_vec(),
            point,
            jacobian: self.determinant.clone(),
            inverse_jacobian: (F(self.determinant.one()) / F(self.determinant.clone())).0,
            residual: max_coordinate_residual_scalar(coordinates, &recovered),
            support: SamplingSupport::Full,
            diagnostics: vec![format!("affine dimension={}", self.dimension)],
        }
        .validate("affine map")
    }

    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.inverse(point, context.previous).map(Some)
    }
}

/// Report from integrating a normalized Gaussian through one sampling map.
///
/// This is intentionally independent of a process integrand. It exercises the
/// map's complete unit-cube-to-raw push-forward and its exact Jacobian, making
/// it suitable for acceptance tests of ordinary, surface and composed maps.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingMapAcceptanceReport {
    pub sample_count: usize,
    pub finite_sample_count: usize,
    pub normalization: f64,
    pub normalization_stderr: f64,
    pub jacobian_min: f64,
    pub jacobian_max: f64,
}

impl SamplingMapAcceptanceReport {
    pub fn normalized_gaussian(
        map: &dyn SamplingMapComponent,
        sample_count: usize,
        width: f64,
        center: &[f64],
    ) -> Result<Self> {
        if sample_count == 0 {
            return Err(eyre!(
                "sampling acceptance harness needs at least one sample"
            ));
        }
        if !width.is_finite() || width <= 0.0 {
            return Err(eyre!(
                "sampling acceptance Gaussian width must be positive and finite"
            ));
        }
        if center.len() != map.output_dimensions() {
            return Err(eyre!(
                "sampling acceptance Gaussian has dimension {}, expected {}",
                center.len(),
                map.output_dimensions()
            ));
        }
        if center.iter().any(|component| !component.is_finite()) {
            return Err(eyre!("sampling acceptance Gaussian centre must be finite"));
        }
        let dimension = map.output_dimensions() as f64;
        let gaussian_normalization =
            (2.0 * std::f64::consts::PI * width * width).powf(-0.5 * dimension);
        let mut report = Self {
            sample_count,
            finite_sample_count: 0,
            normalization: 0.0,
            normalization_stderr: 0.0,
            jacobian_min: f64::INFINITY,
            jacobian_max: f64::NEG_INFINITY,
        };
        let mut square_sum = 0.0;
        for sample in 1..=sample_count {
            let coordinates = (0..map.dimensions())
                .map(|dimension| halton(sample, prime(dimension)))
                .collect::<Vec<_>>();
            let evaluation = map.forward(&coordinates, &mut SamplingMapContext::detached(&[]))?;
            let value = evaluation
                .point
                .iter()
                .zip(center)
                .map(|(point, centre)| (point - centre).powi(2))
                .sum::<f64>();
            let weight =
                gaussian_normalization * (-0.5 * value / width.powi(2)).exp() * evaluation.jacobian;
            report.jacobian_min = report.jacobian_min.min(evaluation.jacobian);
            report.jacobian_max = report.jacobian_max.max(evaluation.jacobian);
            if !weight.is_finite() {
                continue;
            }
            report.finite_sample_count += 1;
            report.normalization += weight;
            square_sum += weight * weight;
        }
        if report.finite_sample_count > 0 {
            // The cube average is over every requested sample.  Dividing by
            // only the finite subset would silently renormalize a map with
            // failed branches and could make an incomplete channel appear
            // normalized.  Non-finite samples therefore contribute zero to
            // this diagnostic while `finite_sample_count` exposes the loss.
            let count = report.sample_count as f64;
            report.normalization /= count;
            report.normalization_stderr =
                ((square_sum / count - report.normalization.powi(2)).max(0.0) / count).sqrt();
        }
        Ok(report)
    }
}

fn prime(index: usize) -> u64 {
    const PRIMES: [u64; 16] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53];
    PRIMES.get(index).copied().unwrap_or(59 + 2 * index as u64)
}

fn halton(mut index: usize, base: u64) -> f64 {
    let mut fraction = 1.0;
    let mut value = 0.0;
    while index > 0 {
        fraction /= base as f64;
        value += fraction * (index as u64 % base) as f64;
        index /= base as usize;
    }
    value
}

/// A block-triangular composition of graph-independent map components.
///
/// Products split both input and output blocks and multiply the exact child
/// determinants.  Ordered compositions pass all previous output blocks as
/// context to each child; their derivative is block triangular, so the exact
/// determinant is still the product of the diagonal child determinants.
/// Cloning recurses into each child's warmed evaluator buffers; immutable
/// geometry callbacks retain their existing shared ownership.
#[derive(Clone)]
pub struct SamplingMapComposition<T: FloatLike = f64> {
    kind: SamplingCompositionKind,
    children: Vec<Box<dyn SamplingMapComponent<T>>>,
    dimensions: usize,
    output_dimensions: usize,
}

/// A product or ordered map whose child output blocks are embedded into an
/// explicitly supplied master-frame permutation.  This is the exact generic
/// construction used for a surface subspace together with its complement:
/// each child remains in its native coordinates, while this wrapper only
/// permutes complete three-momentum blocks.  The permutation has unit
/// absolute determinant, so the product Jacobian is unchanged. A conditional
/// frame transform first prepares the actual preceding outputs, then converts
/// the active physical output back to raw coordinates with its exact determinant.
pub type SamplingMapContextTransform<T> = Arc<
    dyn Fn(&mut SamplingMapContext<'_, T>) -> Result<(Vec<T>, SamplingMapAffine<T>)>
        + Send
        + Sync
        + 'static,
>;

#[derive(Clone)]
pub struct SamplingMapEmbedding<T: FloatLike = f64> {
    map: Box<SamplingMapComposition<T>>,
    /// Product-output index -> master-frame output index.
    output_indices: Vec<usize>,
    context_transform: Option<SamplingMapContextTransform<T>>,
}

impl<T: FloatLike> std::fmt::Debug for SamplingMapEmbedding<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("SamplingMapEmbedding")
            .field("map", &self.map)
            .field("output_indices", &self.output_indices)
            .field("conditional_frame", &self.context_transform.is_some())
            .finish()
    }
}

impl<T: FloatLike> SamplingMapEmbedding<T> {
    /// Build an embedded map from either a direct product or an ordered
    /// conditional composition.  The latter preserves the block context
    /// passed by `then(...)` while this wrapper only permutes its complete
    /// output into the master frame.
    pub fn from_composition(
        map: SamplingMapComposition<T>,
        output_indices: Vec<usize>,
    ) -> Result<Self> {
        if output_indices.len() != map.output_dimensions {
            return Err(eyre!(
                "sampling-map embedding supplies {} output indices, expected {}",
                output_indices.len(),
                map.output_dimensions
            ));
        }
        let mut seen = vec![false; map.output_dimensions];
        for index in &output_indices {
            let Some(slot) = seen.get_mut(*index) else {
                return Err(eyre!(
                    "sampling-map embedding output index {index} is outside master frame 0..{}",
                    map.output_dimensions
                ));
            };
            if *slot {
                return Err(eyre!(
                    "sampling-map embedding repeats master-frame output index {index}"
                ));
            }
            *slot = true;
        }
        Ok(Self {
            map: Box::new(map),
            output_indices,
            context_transform: None,
        })
    }

    /// Prepare only declared prior coordinates. The host must certify that
    /// omitted later coordinates cannot change this cut or target equation.
    pub fn with_context_transform(mut self, transform: SamplingMapContextTransform<T>) -> Self {
        self.context_transform = Some(transform);
        self
    }

    fn transformed_context(
        &self,
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<(Vec<T>, SamplingMapAffine<T>)>> {
        self.context_transform
            .as_ref()
            .map(|transform| {
                if context.previous.iter().any(|value| !value.is_finite()) {
                    return Err(eyre!(
                        "conditional frame needs finite raw prerequisite coordinates"
                    ));
                }
                let (physical, frame) = transform(context)?;
                if physical.iter().any(|value| !value.is_finite()) {
                    return Err(SamplingEvaluationError::Unrepresentable {
                        operation: "conditional physical prerequisites",
                        detail: "native frame preparation produced nonfinite coordinates"
                            .to_owned(),
                    }
                    .into());
                }
                if frame.dimension() != self.output_dimensions() {
                    return Err(eyre!(
                        "conditional frame must supply an active affine map of dimension {}",
                        self.output_dimensions()
                    ));
                }
                Ok((physical, frame))
            })
            .transpose()
    }

    /// Build an embedded direct product.  `output_indices` must be a complete
    /// permutation of `0..sum(child.output_dimensions())`; duplicate or
    /// missing frame coordinates are rejected before any numerical use.
    pub fn product(
        children: Vec<Box<dyn SamplingMapComponent<T>>>,
        output_indices: Vec<usize>,
    ) -> Result<Self> {
        let map = SamplingMapComposition::product(children)?;
        Self::from_composition(map, output_indices)
    }

    pub fn output_indices(&self) -> &[usize] {
        &self.output_indices
    }

    fn embed_point(&self, product_point: &[T]) -> Vec<T> {
        let mut point = vec![product_point[0].zero(); product_point.len()];
        for (product_index, master_index) in self.output_indices.iter().enumerate() {
            point[*master_index] = product_point[product_index].clone();
        }
        point
    }

    fn unembed_point(&self, master_point: &[T]) -> Vec<T> {
        self.output_indices
            .iter()
            .map(|master_index| master_point[*master_index].clone())
            .collect()
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for SamplingMapEmbedding<T> {
    fn dimensions(&self) -> usize {
        self.map.dimensions()
    }

    fn output_dimensions(&self) -> usize {
        self.map.output_dimensions()
    }

    fn contract(&self) -> SamplingMapContract {
        let mut contract = self.map.contract();
        contract.requires_context |= self.context_transform.is_some();
        contract
    }

    fn name(&self) -> &'static str {
        "embedded_product"
    }

    fn forward(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        let prepared = self.transformed_context(context)?;
        let previous = prepared
            .as_ref()
            .map_or(context.previous, |(physical, _)| physical);
        let mut evaluation = self
            .map
            .forward(coordinates, &mut context.reborrow(previous, Some(0)))?;
        if let Some((_, frame)) = prepared {
            let transformed =
                frame.forward(&evaluation.point, &mut context.reborrow(&[], Some(1)))?;
            evaluation.point = transformed.point;
            evaluation.jacobian = (F(evaluation.jacobian) * F(transformed.jacobian)).0;
            evaluation.inverse_jacobian =
                (F(evaluation.inverse_jacobian) * F(transformed.inverse_jacobian)).0;
            evaluation.residual = F(evaluation.residual).max(F(transformed.residual)).0;
            evaluation.diagnostics.extend(transformed.diagnostics);
        }
        evaluation.point = self.embed_point(&evaluation.point);
        evaluation.support = self.contract().support;
        evaluation
            .diagnostics
            .push("unit-determinant master-frame embedding".to_owned());
        evaluation.validate("embedded conditional map")
    }

    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        if point.len() != self.output_dimensions() {
            return Err(eyre!(
                "sampling-map embedding inverse received output dimension {}, expected {}",
                point.len(),
                self.output_dimensions()
            ));
        }
        let raw = self.unembed_point(point);
        let prepared = self.transformed_context(context)?;
        let transformed = prepared
            .as_ref()
            .map(|(_, frame)| frame.inverse(&raw, &[]))
            .transpose()?;
        let previous = prepared
            .as_ref()
            .map_or(context.previous, |(physical, _)| physical);
        let Some(mut evaluation) = self.map.inverse(
            transformed
                .as_ref()
                .map_or(raw.as_slice(), |mapped| &mapped.coordinates),
            &mut context.reborrow(previous, Some(0)),
        )?
        else {
            return Ok(None);
        };
        if let Some(transformed) = transformed {
            evaluation.jacobian = (F(evaluation.jacobian) * F(transformed.jacobian)).0;
            evaluation.inverse_jacobian =
                (F(evaluation.inverse_jacobian) * F(transformed.inverse_jacobian)).0;
            evaluation.residual = F(evaluation.residual).max(F(transformed.residual)).0;
            evaluation.diagnostics.extend(transformed.diagnostics);
        }
        // The inverse's point is the user-supplied master-frame point.
        evaluation.point = point.to_vec();
        evaluation.support = self.contract().support;
        evaluation
            .diagnostics
            .push("unit-determinant master-frame embedding".to_owned());
        evaluation
            .validate("embedded conditional inverse")
            .map(Some)
    }

    fn inverse_density(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<T>> {
        if point.len() != self.output_dimensions() {
            return Err(eyre!(
                "sampling-map embedding inverse received output dimension {}, expected {}",
                point.len(),
                self.output_dimensions()
            ));
        }
        let raw = self.unembed_point(point);
        let prepared = self.transformed_context(context)?;
        let transformed = prepared
            .as_ref()
            .map(|(_, frame)| frame.inverse(&raw, &[]))
            .transpose()?;
        let previous = prepared
            .as_ref()
            .map_or(context.previous, |(physical, _)| physical);
        let Some(density) = self.map.inverse_density(
            transformed
                .as_ref()
                .map_or(raw.as_slice(), |mapped| &mapped.coordinates),
            &mut context.reborrow(previous, Some(0)),
        )?
        else {
            return Ok(None);
        };
        let density = if let Some(frame) = transformed {
            (F(density) * F(frame.inverse_jacobian)).0
        } else {
            density
        };
        SamplingMapEvaluation::validate_inverse_density(density).map(Some)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SamplingCompositionKind {
    Product,
    Then,
}

impl<T: FloatLike> std::fmt::Debug for SamplingMapComposition<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("SamplingMapComposition")
            .field("kind", &self.kind)
            .field("children", &self.children.len())
            .field("dimensions", &self.dimensions)
            .field("output_dimensions", &self.output_dimensions)
            .finish()
    }
}

impl<T: FloatLike> SamplingMapComposition<T> {
    pub fn product(children: Vec<Box<dyn SamplingMapComponent<T>>>) -> Result<Self> {
        Self::new(SamplingCompositionKind::Product, children)
    }

    pub fn then(children: Vec<Box<dyn SamplingMapComponent<T>>>) -> Result<Self> {
        Self::new(SamplingCompositionKind::Then, children)
    }

    fn new(
        kind: SamplingCompositionKind,
        children: Vec<Box<dyn SamplingMapComponent<T>>>,
    ) -> Result<Self> {
        if children.is_empty() {
            return Err(eyre!(
                "sampling-map composition requires at least one child"
            ));
        }
        for (index, child) in children.iter().enumerate() {
            if child.dimensions() == 0 || child.output_dimensions() == 0 {
                return Err(eyre!(
                    "sampling-map composition child {index} ({}) has zero dimensions (input {}, output {})",
                    child.name(),
                    child.dimensions(),
                    child.output_dimensions()
                ));
            }
            if child.dimensions() != child.output_dimensions() {
                return Err(eyre!(
                    "sampling-map composition child {index} ({}) is not square (input {}, output {}); exact block Jacobians require bijective blocks",
                    child.name(),
                    child.dimensions(),
                    child.output_dimensions()
                ));
            }
        }
        let dimensions = children.iter().map(|child| child.dimensions()).sum();
        let output_dimensions = children.iter().map(|child| child.output_dimensions()).sum();
        Ok(Self {
            kind,
            children,
            dimensions,
            output_dimensions,
        })
    }

    pub fn is_product(&self) -> bool {
        self.kind == SamplingCompositionKind::Product
    }

    pub fn is_then(&self) -> bool {
        self.kind == SamplingCompositionKind::Then
    }

    pub fn children(&self) -> &[Box<dyn SamplingMapComponent<T>>] {
        &self.children
    }

    fn validate_input(&self, length: usize, kind: &str) -> Result<()> {
        if length != self.dimensions {
            return Err(eyre!(
                "sampling-map {} received dimension {}, expected {}",
                kind,
                length,
                self.dimensions
            ));
        }
        Ok(())
    }

    fn validate_output(&self, length: usize) -> Result<()> {
        if length != self.output_dimensions {
            return Err(eyre!(
                "sampling-map inverse received output dimension {}, expected {}",
                length,
                self.output_dimensions
            ));
        }
        Ok(())
    }

    fn evaluate_forward(
        &self,
        coordinates: &[T],
        initial_context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        self.validate_input(coordinates.len(), "composition")?;
        let mut offset = 0;
        let mut context = if self.is_then() {
            initial_context.previous.to_vec()
        } else {
            Vec::new()
        };
        let mut evaluations = Vec::with_capacity(self.children.len());
        for (index, child) in self.children.iter().enumerate() {
            let end = offset + child.dimensions();
            let child_coordinates = &coordinates[offset..end];
            let evaluation = child.forward(
                child_coordinates,
                &mut initial_context
                    .reborrow(if self.is_then() { &context } else { &[] }, Some(index)),
            )?;
            if self.is_then() {
                context.extend_from_slice(&evaluation.point);
            }
            evaluations.push(evaluation);
            offset = end;
        }
        let mut evaluation = combine_evaluations(evaluations)?;
        evaluation.support = self.contract().support;
        Ok(evaluation)
    }

    fn evaluate_inverse(
        &self,
        point: &[T],
        initial_context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.validate_output(point.len())?;
        let mut offset = 0;
        let mut context = if self.is_then() {
            initial_context.previous.to_vec()
        } else {
            Vec::new()
        };
        let mut evaluations = Vec::with_capacity(self.children.len());
        // In a triangular map, prior output blocks are known before later
        // blocks are inverted, so inverse traversal remains in declaration
        // order even though the determinant is block triangular.
        for (index, child) in self.children.iter().enumerate() {
            let end = offset + child.output_dimensions();
            let child_point = &point[offset..end];
            let Some(evaluation) = child.inverse(
                child_point,
                &mut initial_context
                    .reborrow(if self.is_then() { &context } else { &[] }, Some(index)),
            )?
            else {
                return Ok(None);
            };
            if self.is_then() {
                context.extend_from_slice(child_point);
            }
            evaluations.push(evaluation);
            offset = end;
        }
        let mut evaluation = combine_evaluations(evaluations)?;
        evaluation.support = self.contract().support;
        Ok(Some(evaluation))
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for SamplingMapComposition<T> {
    fn dimensions(&self) -> usize {
        self.dimensions
    }

    fn output_dimensions(&self) -> usize {
        self.output_dimensions
    }

    fn contract(&self) -> SamplingMapContract {
        let mut contract = self.children.iter().map(|child| child.contract()).fold(
            SamplingMapContract {
                support: SamplingSupport::Full,
                requires_context: false,
                requires_proposal_policy: false,
                jacobian: SamplingJacobian::ExactForward,
            },
            combine_contracts,
        );
        if self.is_then() {
            // Earlier output blocks supply later children's declared context.
            // A conditional first block still needs external preparation;
            // satisfying either dependency never expands restricted coverage.
            contract.requires_context = self.children[0].contract().requires_context;
        }
        contract
    }

    fn name(&self) -> &'static str {
        match self.kind {
            SamplingCompositionKind::Product => "product",
            SamplingCompositionKind::Then => "then",
        }
    }

    fn forward(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        self.evaluate_forward(coordinates, context)
    }

    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        self.evaluate_inverse(point, context)
    }

    fn inverse_density(
        &self,
        point: &[T],
        initial_context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<T>> {
        self.validate_output(point.len())?;
        let mut offset = 0;
        let mut context = if self.is_then() {
            initial_context.previous.to_vec()
        } else {
            Vec::new()
        };
        let mut density = F(point[0].one());
        for (index, child) in self.children.iter().enumerate() {
            let end = offset + child.output_dimensions();
            let child_point = &point[offset..end];
            let Some(child_density) = child.inverse_density(
                child_point,
                &mut initial_context
                    .reborrow(if self.is_then() { &context } else { &[] }, Some(index)),
            )?
            else {
                return Ok(None);
            };
            if self.is_then() {
                context.extend_from_slice(child_point);
            }
            density *= F(child_density);
            offset = end;
        }
        SamplingMapEvaluation::validate_inverse_density(density.0).map(Some)
    }
}

pub(crate) fn combine_contracts(
    left: SamplingMapContract,
    right: SamplingMapContract,
) -> SamplingMapContract {
    let support = match (left.support, right.support) {
        (SamplingSupport::Restricted, _) | (_, SamplingSupport::Restricted) => {
            SamplingSupport::Restricted
        }
        _ => SamplingSupport::Full,
    };
    let jacobian = match (left.jacobian, right.jacobian) {
        (SamplingJacobian::ProxyOnly, _) | (_, SamplingJacobian::ProxyOnly) => {
            SamplingJacobian::ProxyOnly
        }
        (SamplingJacobian::ExactImplicit, _) | (_, SamplingJacobian::ExactImplicit) => {
            SamplingJacobian::ExactImplicit
        }
        _ => SamplingJacobian::ExactForward,
    };
    SamplingMapContract {
        support,
        requires_context: left.requires_context || right.requires_context,
        requires_proposal_policy: left.requires_proposal_policy || right.requires_proposal_policy,
        jacobian,
    }
}

fn combine_evaluations<T: FloatLike>(
    evaluations: Vec<SamplingMapEvaluation<T>>,
) -> Result<SamplingMapEvaluation<T>> {
    let mut coordinates = Vec::new();
    let mut point = Vec::new();
    let one = F(evaluations[0].jacobian.one());
    let mut jacobian = one.clone();
    let mut inverse_jacobian = one.clone();
    let mut residual = one.zero();
    let mut support = SamplingSupport::Full;
    let mut diagnostics = Vec::new();
    for evaluation in evaluations {
        coordinates.extend(evaluation.coordinates);
        point.extend(evaluation.point);
        jacobian *= F(evaluation.jacobian);
        inverse_jacobian *= F(evaluation.inverse_jacobian);
        residual = residual.max(F(evaluation.residual));
        support = match (support, evaluation.support) {
            (SamplingSupport::Restricted, _) | (_, SamplingSupport::Restricted) => {
                SamplingSupport::Restricted
            }
            _ => SamplingSupport::Full,
        };
        diagnostics.extend(evaluation.diagnostics);
    }
    SamplingMapEvaluation {
        coordinates,
        point,
        jacobian: jacobian.0,
        inverse_jacobian: inverse_jacobian.0,
        residual: residual.0,
        support,
        diagnostics,
    }
    .validate("map composition")
}

/// A graph-independent spherical radial proxy around an energy-surface centre.
///
/// The first unit-cube coordinate is mapped to a non-negative radius and the
/// remaining coordinates are uniform hyperspherical angles.  When a regular
/// surface radius is supplied, the radial interval is split at that radius;
/// this keeps the threshold at a deterministic coordinate while retaining
/// full support for all of `R^dimension`.  An absent surface (`None`) uses the
/// same compactification with the threshold radius set to zero.  The map is
/// deliberately independent of graph topology: process code supplies the
/// centre and the classified radius after its kinematics preparation.
///
/// This map has an exact determinant for the spherical proposal it defines,
/// including its rootless fallback. A physical multi-loop E-surface is
/// generally direction-dependent, so this proposal does not by itself certify
/// focusing on that physical threshold. The implicit radial-root kernel
/// supplies the actual energy equation's directional radius when needed;
/// physical focusing and exact proposal-density accounting are distinct.
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceRadialMap<T: FloatLike = f64> {
    dimension: usize,
    center: Vec<T>,
    threshold_radius: Option<T>,
    beta: f64,
    power: f64,
}

/// A full-support radial chart whose target surface can vary with direction.
///
/// The evaluator returns the residual and its derivative with respect to the
/// radial coordinate for a fixed unit direction.  The chart solves one
/// monotone positive root per direction, then uses the same signed power
/// profile as [`SurfaceRadialMap`] on both sides of that root.  The radial
/// derivative is therefore the exact implicit contribution to the determinant;
/// angular derivatives of the direction-dependent root are tangential and do
/// not change the wedge-product determinant.
///
/// A direction with no regular positive root uses the compactified absent-fibre
/// branch.  This preserves full raw support while exposing the branch in the
/// diagnostics returned by the map.  The caller must provide an evaluator with
/// a single increasing regular root on every direction where it expects exact
/// surface targeting.
///
/// Both the root equation and the chart centre may depend on a preceding
/// ordered-map context. In that mode the map is conditional by itself, while
/// a `then` composition whose earlier blocks are full-support can expose a
/// full-support chart: context derivatives are off-diagonal and do not alter
/// the product of active-block determinants.
pub type ImplicitSurfaceRadialEvaluator<T = f64> =
    Arc<dyn Fn(&[T], T) -> Result<(T, T)> + Send + Sync + 'static>;

/// Context-aware variant used by conditional surface charts. The context is
/// the output of an earlier `then(...)`/complement map and may contain the
/// sampled spectator loop momenta needed to solve a partial-subspace surface.
pub type ImplicitSurfaceRadialContextEvaluator<T = f64> =
    Arc<dyn Fn(&[T], T, &[T], &[T]) -> Result<(T, T)> + Send + Sync + 'static>;

/// Complement-only centre and existence classification of a directional chart.
///
/// The context is the output of an earlier ordered map (typically a sampled
/// complement block). Keeping this preparation separate from the scalar
/// surface evaluator makes the dependency explicit: derivatives with respect
/// to the active radial coordinates remain block diagonal, while centre and
/// root changes with the complement occupy only the off-diagonal Jacobian
/// block. The vector has the active dimension and is passed to each radial
/// equation evaluation, so an optimizer is never invoked inside the root solve.
pub type ImplicitSurfaceContextPreparer<T = f64> =
    Arc<dyn Fn(&[T]) -> Result<(Vec<T>, PreparedSurfaceStatus<T>)> + Send + Sync + 'static>;

/// Native fit and worker-local buffers for the optional LU-scale proposal.
/// The physical surface and the scalar inversion remain owned by the same
/// implicit radial map; this record carries no second channel enumeration.
struct LuHProfile<T: FloatLike> {
    settings: SamplingRadialProfile,
    program: Mutex<SamplingExpressionEvaluator>,
    log_scale: F<T>,
    shape: F<T>,
    max_occurrence: usize,
}

impl<T: FloatLike> Clone for LuHProfile<T> {
    fn clone(&self) -> Self {
        Self {
            settings: self.settings.clone(),
            program: Mutex::new(
                self.program
                    .lock()
                    .expect("LU profile evaluator lock poisoned")
                    .clone(),
            ),
            log_scale: self.log_scale.clone(),
            shape: self.shape.clone(),
            max_occurrence: self.max_occurrence,
        }
    }
}

#[derive(Clone)]
pub struct ImplicitSurfaceRadialMap<T: FloatLike = f64> {
    dimension: usize,
    center: Vec<T>,
    beta: f64,
    power: f64,
    evaluator: ImplicitSurfaceRadialEvaluator<T>,
    context_evaluator: Option<ImplicitSurfaceRadialContextEvaluator<T>>,
    context_preparer: Option<ImplicitSurfaceContextPreparer<T>>,
    prepared_status: Option<PreparedSurfaceStatus<T>>,
    root_tolerance: f64,
    lu_h_profile: Option<LuHProfile<T>>,
}

impl<T: FloatLike> std::fmt::Debug for ImplicitSurfaceRadialMap<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ImplicitSurfaceRadialMap")
            .field("dimension", &self.dimension)
            .field("center", &self.center)
            .field("beta", &self.beta)
            .field("power", &self.power)
            .field("context_dependent", &self.context_evaluator.is_some())
            .field("center_context_dependent", &self.context_preparer.is_some())
            .field("root_tolerance", &self.root_tolerance)
            .finish_non_exhaustive()
    }
}

impl<T: FloatLike> ImplicitSurfaceRadialMap<T> {
    pub fn new(
        dimension: usize,
        center: Vec<T>,
        beta: f64,
        power: f64,
        evaluator: ImplicitSurfaceRadialEvaluator<T>,
    ) -> Result<Self> {
        if dimension < 2 {
            return Err(eyre!(
                "implicit surface radial map requires dimension at least two"
            ));
        }
        if center.len() != dimension {
            return Err(eyre!(
                "implicit surface radial-map centre has dimension {}, expected {dimension}",
                center.len()
            ));
        }
        if center.iter().any(|value| !value.is_finite()) {
            return Err(eyre!("implicit surface radial-map centre must be finite"));
        }
        if !beta.is_finite() || beta <= 0.0 {
            return Err(eyre!(
                "implicit surface radial-map compactification scale must be positive and finite"
            ));
        }
        if !power.is_finite() || power <= 0.0 {
            return Err(eyre!(
                "implicit surface radial-map radial power must be positive and finite"
            ));
        }
        Ok(Self {
            dimension,
            center,
            beta,
            power,
            evaluator,
            context_evaluator: None,
            context_preparer: None,
            prepared_status: None,
            root_tolerance: 1.0e-11,
            lu_h_profile: None,
        })
    }

    /// Bind one named proposal after cloning the shared physical geometry.
    /// Only an origin-centred full LU dilation has t*=R/r; shifted or conditional
    /// charts must not silently use this auxiliary-scale construction.
    pub(crate) fn with_lu_h_profile(
        mut self,
        mut program: SamplingExpressionEvaluator,
        settings: &SamplingRadialProfile,
        max_occurrence: usize,
    ) -> Result<Self> {
        if self.center.iter().any(|value| *value != value.zero())
            || self.context_evaluator.is_some()
            || self.context_preparer.is_some()
            || max_occurrence == 0
        {
            return Err(eyre!(
                "LU-h sampling requires an origin-centred full-frame physical cut and a positive residue order"
            ));
        }
        if program.parameter_count() != 5
            || program.output_count() != 5
            || program.derivative_parameters() != [0]
        {
            return Err(eyre!(
                "invalid compiled LU-h profile program: expected five inputs/outputs and active derivative column [0], got {} inputs, {} outputs and {:?}",
                program.parameter_count(),
                program.output_count(),
                program.derivative_parameters(),
            ));
        }
        let zero = self.center[0].zero();
        let fit = program.evaluate(&vec![zero; 5])?;
        let log_scale = fit[0].re.clone();
        let shape = fit[1].re.clone();
        if !log_scale.0.is_finite() || !shape.0.is_finite() || shape <= shape.zero() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "native LU h profile fit",
                detail: "finite log-scale and strictly positive shape cannot be represented"
                    .to_owned(),
            }
            .into());
        }
        self.lu_h_profile = Some(LuHProfile {
            settings: settings.clone(),
            program: Mutex::new(program),
            log_scale,
            shape,
            max_occurrence,
        });
        Ok(self)
    }

    pub fn lu_h_max_occurrence(&self) -> Option<usize> {
        self.lu_h_profile
            .as_ref()
            .map(|profile| profile.max_occurrence)
    }

    pub fn with_root_tolerance(mut self, root_tolerance: f64) -> Result<Self> {
        if !root_tolerance.is_finite() || root_tolerance <= 0.0 {
            return Err(eyre!(
                "implicit surface radial-map root tolerance must be positive and finite"
            ));
        }
        self.root_tolerance = root_tolerance;
        Ok(self)
    }

    /// Attach a conditional evaluator while retaining the same exact radial
    /// determinant. Such a map must be composed after the map producing its
    /// context and therefore advertises conditional support.
    pub fn with_context_evaluator(
        mut self,
        evaluator: ImplicitSurfaceRadialContextEvaluator<T>,
    ) -> Self {
        self.context_evaluator = Some(evaluator);
        self
    }

    /// Attach a context-dependent centre to the chart.  The map remains
    /// exact in its active coordinates because the centre depends only on
    /// already sampled context, and therefore contributes no extra diagonal
    /// determinant factor in an ordered composition.
    pub fn with_context_preparer(mut self, evaluator: ImplicitSurfaceContextPreparer<T>) -> Self {
        self.context_preparer = Some(evaluator);
        self
    }

    /// Freeze a prepared complement at binding time. This is useful for a
    /// full-active graph surface: center finding runs once during warmup,
    /// while the same scalar radial engine retains its explicit classification.
    pub(crate) fn freeze_context(mut self, context: &[T]) -> Result<Self> {
        let (center, status) = self.prepare_context(context)?;
        if let Some(evaluator) = self.context_evaluator.take() {
            let center = center.clone();
            let context = context.to_vec();
            self.evaluator =
                Arc::new(move |direction, radius| evaluator(direction, radius, &center, &context));
        }
        self.center = center;
        self.prepared_status = status;
        self.context_preparer = None;
        Ok(self)
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn center(&self) -> &[T] {
        &self.center
    }

    pub fn beta(&self) -> f64 {
        self.beta
    }

    pub fn power(&self) -> f64 {
        self.power
    }

    pub(crate) fn prepare_context(
        &self,
        context: &[T],
    ) -> Result<(Vec<T>, Option<PreparedSurfaceStatus<T>>)> {
        if context.iter().any(|value| !value.is_finite()) {
            return Err(eyre!("implicit surface context must contain finite values"));
        }
        let (center, status) = if let Some(evaluator) = &self.context_preparer {
            let (center, status) = evaluator(context)?;
            let status = match status {
                PreparedSurfaceStatus::Existing { normalized_margin } => {
                    PreparedSurfaceStatus::existing(normalized_margin)?
                }
                PreparedSurfaceStatus::Pinched { normalized_margin } => {
                    PreparedSurfaceStatus::pinched(normalized_margin)?
                }
                PreparedSurfaceStatus::Absent { reason } => PreparedSurfaceStatus::absent(reason)?,
            };
            (center, Some(status))
        } else {
            (self.center.clone(), self.prepared_status.clone())
        };
        if center.len() != self.dimension {
            return Err(eyre!(
                "implicit surface radial-map centre has dimension {}, expected {}",
                center.len(),
                self.dimension
            ));
        }
        if center.iter().any(|value| !value.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "implicit surface context centre",
                detail: "derived centre contains nonfinite values".to_owned(),
            }
            .into());
        }
        Ok((center, status))
    }

    pub fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: SamplingSupport::Full,
            requires_context: self.context_evaluator.is_some() || self.context_preparer.is_some(),
            requires_proposal_policy: false,
            jacobian: SamplingJacobian::ExactImplicit,
        }
    }

    pub fn forward(&self, coordinates: &[T]) -> Result<SamplingMapEvaluation<T>> {
        self.forward_with_context(coordinates, &[])
    }

    pub fn forward_with_context(
        &self,
        coordinates: &[T],
        context: &[T],
    ) -> Result<SamplingMapEvaluation<T>> {
        self.validate_coordinates(coordinates)?;
        let (direction, angular_jacobian) = direction_from_coordinates(coordinates)?;
        let (center, status) = self.prepare_context(context)?;
        let root = self.root_for_direction(&direction, &center, context, status.as_ref())?;
        let (radius, radial_jacobian) =
            self.radius_from_coordinate(coordinates[0].clone(), root.clone())?;
        let radius = F(radius);
        if self.lu_h_profile.is_none()
            && self.power != 1.0
            && root
                .as_ref()
                .is_some_and(|(threshold, _)| &radius.0 == threshold)
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "implicit radial map",
                detail: "radius rounds to the singular threshold seam".to_owned(),
            }
            .into());
        }
        let point = center
            .iter()
            .zip(&direction)
            .map(|(center, direction)| (F(center.clone()) + &radius * F(direction.clone())).0)
            .collect::<Vec<_>>();
        let jacobian =
            F(radial_jacobian) * F(angular_jacobian) * radius.powi(self.dimension as i32 - 1);
        let inverse_jacobian = jacobian.clone().inv();
        if !radius.0.is_finite()
            || radius <= radius.zero()
            || point.iter().any(|component| !component.is_finite())
            || [&jacobian, &inverse_jacobian]
                .iter()
                .any(|value| !value.0.is_finite() || **value <= value.zero())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "implicit radial map",
                detail: format!("point or positive Jacobian pair ({jacobian}, {inverse_jacobian})"),
            }
            .into());
        }
        Ok(SamplingMapEvaluation {
            coordinates: coordinates.to_vec(),
            point,
            jacobian: jacobian.0,
            inverse_jacobian: inverse_jacobian.0,
            residual: radius.zero().0,
            support: self.contract().support,
            diagnostics: {
                let mut diagnostics = vec![if self
                    .lu_h_profile
                    .as_ref()
                    .is_some_and(|profile| profile.settings.broad_fraction == 1.0)
                {
                    "lu_h:broad_only_root_not_required".to_owned()
                } else if root.is_some() {
                    "implicit_surface:regular_root".to_owned()
                } else if matches!(status, Some(PreparedSurfaceStatus::Pinched { .. })) {
                    "implicit_surface:pinched_fallback".to_owned()
                } else {
                    "implicit_surface:absent_fallback".to_owned()
                }];
                if let Some(profile) = &self.lu_h_profile {
                    diagnostics.push(format!(
                        "lu_h:log_scale={},shape={},broad_fraction={},max_occurrence={}",
                        profile.log_scale,
                        profile.shape,
                        profile.settings.broad_fraction,
                        profile.max_occurrence
                    ));
                }
                diagnostics
            },
        })
    }

    pub fn inverse(&self, point: &[T]) -> Result<SamplingMapEvaluation<T>> {
        self.inverse_with_context(point, &[])
    }

    pub fn inverse_with_context(
        &self,
        point: &[T],
        context: &[T],
    ) -> Result<SamplingMapEvaluation<T>> {
        self.inverse_evaluation(point, context, true)
    }

    fn inverse_evaluation(
        &self,
        point: &[T],
        context: &[T],
        reconstruct: bool,
    ) -> Result<SamplingMapEvaluation<T>> {
        if point.len() != self.dimension {
            return Err(eyre!(
                "implicit surface radial-map inverse received dimension {}, expected {}",
                point.len(),
                self.dimension
            ));
        }
        if point.iter().any(|component| !component.is_finite()) {
            return Err(eyre!(
                "implicit surface radial-map inverse requires finite momenta"
            ));
        }
        let (center, status) = self.prepare_context(context)?;
        let displacement = point
            .iter()
            .zip(&center)
            .map(|(point, center)| F(point.clone()) - F(center.clone()))
            .collect::<Vec<_>>();
        let radius = displacement
            .iter()
            .map(F::square)
            .fold(displacement[0].zero(), |sum, value| sum + value)
            .sqrt();
        if !radius.0.is_finite() || radius <= radius.zero() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "implicit radial inverse",
                detail: "radius is non-finite or at the chart centre".to_owned(),
            }
            .into());
        }
        let direction = displacement
            .iter()
            .map(|value| (value / &radius).0)
            .collect::<Vec<_>>();
        let root = self.root_for_direction(&direction, &center, context, status.as_ref())?;
        let (coordinate, radial_jacobian) =
            self.coordinate_from_radius(radius.0.clone(), root.clone())?;
        let mut coordinates = coordinates_from_direction(&direction)?;
        coordinates[0] = coordinate;
        self.validate_coordinates(&coordinates).map_err(|error| {
            SamplingEvaluationError::Unrepresentable {
                operation: "implicit radial inverse coordinates",
                detail: error.to_string(),
            }
        })?;
        // The angular factor and supplied-radius derivative define the density.
        // Foreign partition queries need neither a second inverse-CDF solve
        // nor its Cartesian roundtrip, which remain available on full inverses.
        let (recovered_direction, angular) = direction_from_coordinates(&coordinates)?;
        let residual = if reconstruct {
            // Reuse the prepared fiber and root for the roundtrip diagnostic;
            // re-entering forward would repeat its complement-only optimization.
            let (recovered_radius, _) =
                self.radius_from_coordinate(coordinates[0].clone(), root)?;
            let reconstructed = center
                .iter()
                .zip(&recovered_direction)
                .map(|(center, direction)| {
                    (F(center.clone()) + F(recovered_radius.clone()) * F(direction.clone())).0
                })
                .collect::<Vec<_>>();
            max_coordinate_residual_scalar(point, &reconstructed)
        } else {
            // Private placeholder, discarded by the scalar density query; full
            // inverses always calculate the actual reconstruction residual.
            radius.zero().0
        };
        if !residual.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "implicit surface radial-map inverse residual is not representable at the current precision".to_owned() }.into());
        }
        // Evaluate the density at the supplied momentum. Reprojecting its
        // recovered cube coordinate can move it on a sharply varying map.
        let jacobian = F(radial_jacobian) * F(angular) * radius.powi(self.dimension as i32 - 1);
        let inverse_jacobian = jacobian.inv();
        if [&jacobian, &inverse_jacobian]
            .iter()
            .any(|value| !value.0.is_finite() || **value <= value.zero())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "implicit radial inverse density",
                detail: "positive Jacobian pair is not representable at the supplied momentum"
                    .to_owned(),
            }
            .into());
        }
        Ok(SamplingMapEvaluation {
            coordinates,
            point: point.to_vec(),
            jacobian: jacobian.0,
            inverse_jacobian: inverse_jacobian.0,
            residual,
            support: self.contract().support,
            diagnostics: vec!["implicit_surface:inverse_supplied_radius".to_owned()],
        })
    }

    fn validate_coordinates(&self, coordinates: &[T]) -> Result<()> {
        if coordinates.len() != self.dimension {
            return Err(eyre!(
                "implicit surface radial-map forward received {}, expected {} coordinates",
                coordinates.len(),
                self.dimension
            ));
        }
        if coordinates.iter().any(|coordinate| {
            !coordinate.is_finite()
                || *coordinate <= coordinate.zero()
                || *coordinate >= coordinate.one()
        }) {
            return Err(eyre!(
                "implicit surface radial-map coordinates must be finite and strictly inside the unit cube"
            ));
        }
        Ok(())
    }

    fn evaluate(&self, direction: &[T], radius: T, center: &[T], context: &[T]) -> Result<(T, T)> {
        if let Some(evaluator) = &self.context_evaluator {
            evaluator(direction, radius, center, context)
        } else {
            (self.evaluator)(direction, radius)
        }
    }

    fn root_for_direction(
        &self,
        direction: &[T],
        center: &[T],
        context: &[T],
        status: Option<&PreparedSurfaceStatus<T>>,
    ) -> Result<Option<(T, T)>> {
        if self
            .lu_h_profile
            .as_ref()
            .is_some_and(|profile| profile.settings.broad_fraction == 1.0)
        {
            // The normalized raw broad law has no R dependence at all.
            return Ok(None);
        }
        if status.is_some_and(PreparedSurfaceStatus::uses_full_support_fallback) {
            return Ok(None);
        }
        let zero = F(direction[0].zero());
        let (origin_value, origin_derivative) =
            self.evaluate(direction, zero.0.clone(), center, context)?;
        let origin_value = F(origin_value);
        if !origin_value.0.is_finite() || !origin_derivative.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "implicit root origin",
                detail: "evaluator returned non-finite origin data".to_owned(),
            }
            .into());
        }
        let scale = origin_value.abs().max(zero.one());
        let requested_tolerance = F::<T>::from_f64(self.root_tolerance);
        let native_tolerance = zero.epsilon() * zero.from_i64(64);
        let native_relative_tolerance = if requested_tolerance < native_tolerance {
            requested_tolerance.clone()
        } else {
            native_tolerance
        };
        let residual_tolerance = &native_relative_tolerance * &scale;
        if origin_value.abs() <= residual_tolerance {
            return Err(SamplingEvaluationError::UncertainGeometry {
                detail: format!("center sign is not certified at the current root tolerance: value={origin_value}, tolerance={residual_tolerance}"),
            }.into());
        }
        if origin_value > zero {
            if status.is_some_and(PreparedSurfaceStatus::is_existing) {
                return Err(SamplingEvaluationError::UncertainGeometry {
                    detail: "prepared existing fiber center is not strictly interior".to_owned(),
                }
                .into());
            }
            return Ok(None);
        }
        // An LU-h proposal uses R(direction) only as a positive focusing scale:
        // any deterministic positive approximation defines the same normalized
        // radial CDF construction. Its angular derivatives do not enter the
        // spherical determinant. Honor the declared root accuracy for that law,
        // identically in forward and inverse; physical host roots and singular
        // threshold-shell roots retain their native-precision requirements.
        // Keep origin classification above unchanged, and retain a strict
        // interior bracket even for a cut very close to production threshold.
        let relative_tolerance = if self.lu_h_profile.is_some() {
            let interior_tolerance = origin_value.abs() / (&scale * zero.from_i64(64));
            if requested_tolerance < interior_tolerance {
                requested_tolerance
            } else {
                interior_tolerance
            }
        } else {
            native_relative_tolerance
        };
        // Reuse the native bracket-preserving solver. Callback failures keep
        // their original structural/numerical classification across its scalar
        // interface; no failed equation evaluation becomes a rootless branch.
        let callback_error = std::cell::RefCell::new(None);
        let result = safeguarded_newton_iteration_and_derivative(
            &zero,
            &F::<T>::from_f64(self.beta.max(1.0)),
            |radius| match self.evaluate(direction, radius.0.clone(), center, context) {
                Ok((value, derivative)) => (F(value), F(derivative)),
                Err(error) => {
                    *callback_error.borrow_mut() = Some(error);
                    (F::<T>::from_f64(f64::NAN), F::<T>::from_f64(f64::NAN))
                }
            },
            &(&relative_tolerance / zero.epsilon()),
            2048,
            96,
            &scale,
        );
        if let Some(error) = callback_error.into_inner() {
            return Err(error);
        }
        let result = result.map_err(|error| SamplingEvaluationError::UncertifiedRoot {
            detail: format!("implicit surface radial root failed residual certification: {error}"),
        })?;
        Ok(Some((result.solution.0, result.derivative_at_solution.0)))
    }

    fn lu_h_values(&self, y: &F<T>, log_root: &F<T>) -> Result<Vec<SamplingDualValue<T>>> {
        let profile = self.lu_h_profile.as_ref().expect("bound LU profile");
        let log_beta = F::<T>::from_f64(self.beta).ln();
        let one = y.one();
        let focus_sign = if y <= &profile.log_scale {
            one.clone()
        } else {
            -&one
        };
        let broad_sign = if y <= &(log_root - &log_beta) {
            one.clone()
        } else {
            -one
        };
        profile
            .program
            .lock()
            .map_err(|_| eyre!("LU profile evaluator lock poisoned"))?
            .evaluate_with_derivatives(&[
                y.0.clone(),
                log_root.0.clone(),
                log_beta.0,
                focus_sign.0,
                broad_sign.0,
            ])
    }

    fn radius_from_coordinate(&self, coordinate: T, root: Option<(T, T)>) -> Result<(T, T)> {
        let u = F(coordinate);
        if let Some(profile) = &self.lu_h_profile {
            let one = u.one();
            let beta = F::<T>::from_f64(self.beta);
            let Some((root, _)) = root else {
                // Keep the CDF convention u=G(t) even for the root-independent
                // power-one fallback: r=beta*(1-u)/u.
                return Ok(((&beta * (&one - &u) / &u).0, (beta / u.square()).0));
            };
            let log_root = F(root).ln();
            let log_odds = u.ln() - (&one - &u).ln();
            let focused = &profile.log_scale + &log_odds / &profile.shape;
            let broad = &log_root - beta.ln() + log_odds;
            let y = if profile.settings.broad_fraction == 0.0 || focused == broad {
                focused
            } else {
                let lower = if focused < broad {
                    focused.clone()
                } else {
                    broad.clone()
                };
                let upper = focused.max(broad);
                let upper_tail = u > &one / u.from_i64(2);
                let tail = if upper_tail { &one - &u } else { u.clone() };
                let equation = |y: &F<T>| -> Result<(F<T>, F<T>)> {
                    let values = self.lu_h_values(y, &log_root)?;
                    let value = &values[if upper_tail { 3 } else { 2 }];
                    Ok(if upper_tail {
                        (
                            &one - &value.value.re / &tail,
                            -&value.derivatives[0].re / &tail,
                        )
                    } else {
                        (
                            &value.value.re / &tail - &one,
                            &value.derivatives[0].re / &tail,
                        )
                    })
                };
                let tolerance = u.from_i64(64);
                let maximum_residual = u.epsilon() * &tolerance;
                let (lower_value, _) = equation(&lower)?;
                let (upper_value, _) = equation(&upper)?;
                // Coincident component quantiles and an endpoint root are
                // legitimate; the shared solver otherwise needs strict signs.
                if lower_value.abs() <= maximum_residual {
                    lower
                } else if upper_value.abs() <= maximum_residual {
                    upper
                } else {
                    let callback_error = std::cell::RefCell::new(None);
                    let result = safeguarded_newton_iteration_and_derivative(
                        &lower,
                        &upper,
                        |y| match equation(y) {
                            Ok(values) => values,
                            Err(error) => {
                                *callback_error.borrow_mut() = Some(error);
                                (F::<T>::from_f64(f64::NAN), F::<T>::from_f64(f64::NAN))
                            }
                        },
                        &tolerance,
                        2048,
                        0,
                        &one,
                    );
                    if let Some(error) = callback_error.into_inner() {
                        return Err(error);
                    }
                    result
                        .map_err(|error| SamplingEvaluationError::UncertifiedRoot {
                            detail: format!("LU-h mixture inverse CDF failed: {error}"),
                        })?
                        .solution
                }
            };
            let values = self.lu_h_values(&y, &log_root)?;
            let radius = F(values[4].value.re.0.exp());
            let derivative = &values[2].derivatives[0].re;
            if !derivative.0.is_finite() || derivative <= &u.zero() {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "LU-h inverse CDF",
                    detail: "positive CDF derivative is not representable".to_owned(),
                }
                .into());
            }
            // dr/du = r*|d(log r)/dy|/(dG/dy), obtained from the same eager
            // dual program as G. No derivative of Newton iterations is used.
            let jacobian = &radius * values[4].derivatives[0].re.abs() / derivative;
            return Ok((radius.0, jacobian.0));
        }
        let threshold = root.map(|(radius, _)| radius).unwrap_or_else(|| u.0.zero());
        let (radius, jacobian) = SurfaceRadialMap::<T>::radius_from_coordinate(
            &u,
            F(threshold),
            F::from_f64(self.beta),
            F::from_f64(self.power),
        );
        Ok((radius.0, jacobian.0))
    }

    fn coordinate_from_radius(&self, radius: T, root: Option<(T, T)>) -> Result<(T, T)> {
        let radius = F(radius);
        let beta = F::<T>::from_f64(self.beta);
        if self.lu_h_profile.is_some() {
            let Some((root, _)) = root else {
                return Ok((
                    (&beta / (&beta + &radius)).0,
                    ((&beta + &radius).square() / beta).0,
                ));
            };
            let log_root = F(root).ln();
            let y = &log_root - radius.ln();
            let values = self.lu_h_values(&y, &log_root)?;
            let derivative = &values[2].derivatives[0].re;
            if !derivative.0.is_finite() || derivative <= &radius.zero() {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "LU-h density",
                    detail: "positive CDF derivative is not representable".to_owned(),
                }
                .into());
            }
            return Ok((values[2].value.re.0.clone(), (radius / derivative).0));
        }
        let threshold = F(root
            .map(|(radius, _)| radius)
            .unwrap_or_else(|| radius.0.zero()));
        let power = F::<T>::from_f64(self.power);
        let (coordinate, jacobian) =
            SurfaceRadialMap::<T>::coordinate_from_radius(&radius, threshold, beta, power);
        Ok((coordinate.0, jacobian.0))
    }
}

fn direction_from_coordinates<T: FloatLike>(coordinates: &[T]) -> Result<(Vec<T>, T)> {
    if coordinates.len() < 2 {
        return Err(eyre!(
            "spherical direction requires at least two coordinates"
        ));
    }
    let dimension = coordinates.len();
    let zero = F(coordinates[0].zero());
    let one = zero.one();
    let two = zero.from_i64(2);
    let mut direction = vec![zero.clone(); dimension];
    let mut angular_jacobian = zero.TAU();
    let mut base = one.clone();
    for (i, coordinate) in coordinates[2..].iter().enumerate() {
        let cos_theta = -&one + &two * F(coordinate.clone());
        if cos_theta <= -&one || cos_theta >= one {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "spherical direction",
                detail: "polar coordinate rounds to a singular spherical chart boundary".to_owned(),
            }
            .into());
        }
        let sin_theta = (&one - cos_theta.square()).sqrt();
        angular_jacobian *= &two;
        let angular_power = dimension - 3 - i;
        if angular_power > 0 {
            angular_jacobian *= sin_theta.powi(angular_power as i32);
        }
        direction[i] = &base * cos_theta;
        base *= sin_theta;
    }
    let phi = zero.TAU() * F(coordinates[1].clone());
    direction[dimension - 2] = &base * F(phi.0.cos());
    direction[dimension - 1] = &base * F(phi.0.sin());
    Ok((
        direction.into_iter().map(|value| value.0).collect(),
        angular_jacobian.0,
    ))
}

fn coordinates_from_direction<T: FloatLike>(direction: &[T]) -> Result<Vec<T>> {
    if direction.len() < 2 {
        return Err(eyre!(
            "spherical direction requires at least two components"
        ));
    }
    let dimension = direction.len();
    let zero = F(direction[0].zero());
    let one = zero.one();
    let two = zero.from_i64(2);
    let mut coordinates = vec![zero.clone(); dimension];
    let mut base = one.clone();
    for i in 0..dimension - 2 {
        let cos_theta = F(direction[i].clone()) / &base;
        if base <= zero || !cos_theta.0.is_finite() || cos_theta <= -&one || cos_theta >= one {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "spherical inverse",
                detail: "direction lies on a singular spherical chart boundary".to_owned(),
            }
            .into());
        }
        coordinates[2 + i] = (&one + &cos_theta) / &two;
        base *= (&one - cos_theta.square()).sqrt();
    }
    let mut phi = F(direction[dimension - 1].atan2(&direction[dimension - 2]));
    if phi < zero {
        phi += zero.TAU();
    }
    coordinates[1] = phi / zero.TAU();
    if !coordinates[1].0.is_finite() || coordinates[1] <= zero || coordinates[1] >= one {
        return Err(SamplingEvaluationError::Unrepresentable {
            operation: "spherical inverse",
            detail: "direction lies on the spherical azimuth seam".to_owned(),
        }
        .into());
    }
    Ok(coordinates.into_iter().map(|value| value.0).collect())
}

/// Values returned by [`SurfaceRadialMap::forward`] and `inverse`.
#[derive(Clone, Debug)]
pub struct SurfaceRadialPoint<T: FloatLike> {
    pub coordinates: Vec<F<T>>,
    pub point: Vec<F<T>>,
    pub radius: F<T>,
    /// Exact positive determinant of the forward map in the unit cube frame.
    pub jacobian: F<T>,
    /// Its reciprocal (where the determinant is non-zero).
    pub inverse_jacobian: F<T>,
    /// Infinity norm of the independent coordinate/point round trip residual.
    pub residual: F<T>,
}

impl<T: FloatLike> SurfaceRadialMap<T> {
    /// Construct a radial map in `dimension >= 2` dimensions.
    ///
    /// `threshold_radius = Some(r*)` places the surface at `r*`; `None` is the
    /// full-support absent-fibre fallback. `beta` is a positive compactification
    /// scale and `power` controls concentration on both sides of the threshold.
    /// For power greater than one the density diverges as
    /// `|r - threshold_radius|^(1 / power - 1)` at a nonzero threshold.
    /// Ordinary soft channels retain separate coverage of the origin.
    pub fn new(
        dimension: usize,
        center: Vec<T>,
        threshold_radius: Option<T>,
        beta: f64,
        power: f64,
    ) -> Result<Self> {
        if dimension < 2 {
            return Err(eyre!("surface radial map requires dimension at least two"));
        }
        if center.len() != dimension {
            return Err(eyre!(
                "surface radial-map centre has dimension {}, expected {dimension}",
                center.len()
            ));
        }
        if center.iter().any(|component| !component.is_finite()) {
            return Err(eyre!("surface radial-map centre must be finite"));
        }
        if threshold_radius
            .as_ref()
            .is_some_and(|radius| !radius.is_finite() || radius < &radius.zero())
        {
            return Err(eyre!(
                "surface radial-map threshold radius must be finite and non-negative"
            ));
        }
        if !beta.is_finite() || beta <= 0.0 {
            return Err(eyre!(
                "surface radial-map compactification scale beta must be positive and finite"
            ));
        }
        if !power.is_finite() || power <= 0.0 {
            return Err(eyre!(
                "surface radial-map radial power must be positive and finite"
            ));
        }
        Ok(Self {
            dimension,
            center,
            threshold_radius,
            beta,
            power,
        })
    }

    /// Construct the absent-fibre full-support fallback directly.
    pub fn absent(dimension: usize, center: Vec<T>, beta: f64, power: f64) -> Result<Self> {
        Self::new(dimension, center, None, beta, power)
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn center(&self) -> &[T] {
        &self.center
    }

    pub fn threshold_radius(&self) -> Option<T> {
        self.threshold_radius.clone()
    }

    pub fn beta(&self) -> f64 {
        self.beta
    }

    pub fn power(&self) -> f64 {
        self.power
    }

    pub fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: SamplingSupport::Full,
            requires_context: false,
            requires_proposal_policy: false,
            jacobian: SamplingJacobian::ExactForward,
        }
    }

    /// Map a unit-cube point to a point around the surface centre.
    pub fn forward(&self, coordinates: &[F<T>]) -> Result<SurfaceRadialPoint<T>> {
        self.validate_coordinates(coordinates)?;
        let (threshold, beta, power) = self.radial_parameters();
        let (radius, radial_jacobian) =
            Self::radius_from_coordinate(&coordinates[0], threshold.clone(), beta, power);
        if self.power != 1.0 && threshold > threshold.zero() && radius == threshold {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "surface radial-map radius rounds to the singular threshold seam at the current precision".to_owned() }.into());
        }
        let one = radius.one();
        let two = radius.from_i64(2);
        let phi = radius.TAU() * &coordinates[1];
        let mut direction = vec![radius.zero(); self.dimension];
        let mut angular_jacobian = radius.TAU();
        let mut base = radius.clone();
        for (i, coordinate) in coordinates[2..].iter().enumerate() {
            let cos_theta = -&one + &two * coordinate;
            if cos_theta <= -&one || cos_theta >= one {
                return Err(SamplingEvaluationError::Unrepresentable {
                    operation: "surface radial direction",
                    detail: "polar coordinate rounds to a singular spherical chart boundary"
                        .to_owned(),
                }
                .into());
            }
            let sin_theta = (&one - cos_theta.square()).sqrt();
            angular_jacobian *= &two;
            // Coordinates are uniform in cos(theta), so this exponent is one
            // lower than for the theta parameterisation.
            let angular_power = self.dimension - 3 - i;
            if angular_power > 0 {
                angular_jacobian *= sin_theta.powi(angular_power as i32);
            }
            direction[i] = &base * &cos_theta;
            base *= sin_theta;
        }
        direction[self.dimension - 2] = &base * F(phi.0.cos());
        direction[self.dimension - 1] = &base * F(phi.0.sin());

        let mut point = Vec::with_capacity(self.dimension);
        for (component, centre) in direction.iter().zip(&self.center) {
            point.push(component + F(centre.clone()));
        }
        let jacobian =
            radial_jacobian * angular_jacobian * radius.powi((self.dimension - 1) as i32);
        let inverse_jacobian = jacobian.clone().inv();
        if !radius.0.is_finite()
            || radius <= radius.zero()
            || point.iter().any(|component| !component.0.is_finite())
            || [&jacobian, &inverse_jacobian]
                .iter()
                .any(|value| !value.0.is_finite() || **value <= value.zero())
        {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "surface radial-map point or positive Jacobian pair is not representable at the current precision".to_owned() }.into());
        }
        Ok(SurfaceRadialPoint {
            coordinates: coordinates.to_vec(),
            point,
            radius,
            jacobian,
            inverse_jacobian,
            residual: F::from_f64(0.0),
        })
    }

    /// Invert a point in the map's full-support domain.
    pub fn inverse(&self, point: &[F<T>]) -> Result<SurfaceRadialPoint<T>> {
        if point.len() != self.dimension {
            return Err(eyre!(
                "surface radial-map inverse received dimension {}, expected {}",
                point.len(),
                self.dimension
            ));
        }
        if point.iter().any(|component| !component.0.is_finite()) {
            return Err(eyre!("surface radial-map inverse requires finite momenta"));
        }
        let mut displacement = Vec::with_capacity(self.dimension);
        for (component, centre) in point.iter().zip(&self.center) {
            displacement.push(component - F(centre.clone()));
        }
        let mut radius_squared = displacement[0].square();
        for component in &displacement[1..] {
            radius_squared += component.square();
        }
        let radius = radius_squared.sqrt();
        let zero = radius.zero();
        if displacement.iter().all(|component| component == &zero) {
            return Err(eyre!(
                "surface radial-map inverse is undefined at its centre"
            ));
        }
        if !radius.0.is_finite() || radius <= zero {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "surface radial inverse radius",
                detail: "finite nonzero displacement has an unrepresentable positive norm"
                    .to_owned(),
            }
            .into());
        }
        let (threshold, beta, power) = self.radial_parameters();
        let (radial_coordinate, radial_jacobian) =
            Self::coordinate_from_radius(&radius, threshold, beta, power);
        let one = radius.one();
        let two = radius.from_i64(2);
        let mut coordinates = vec![zero.clone(); self.dimension];
        let mut base = radius.clone();
        let mut angular_jacobian = radius.TAU();
        for i in 0..self.dimension - 2 {
            let cos_theta = &displacement[i] / &base;
            coordinates[2 + i] = (&one + &cos_theta) / &two;
            let sine = (&one - cos_theta.square()).sqrt();
            angular_jacobian *= &two;
            let angular_power = self.dimension - 3 - i;
            if angular_power > 0 {
                angular_jacobian *= sine.powi(angular_power as i32);
            }
            base *= sine;
        }
        let mut phi = F(displacement[self.dimension - 1]
            .0
            .atan2(&displacement[self.dimension - 2].0));
        if phi < zero {
            phi += radius.TAU();
        }
        coordinates[0] = radial_coordinate;
        coordinates[1] = phi / radius.TAU();
        self.validate_coordinates(&coordinates).map_err(|error| {
            SamplingEvaluationError::Unrepresentable {
                operation: "surface radial inverse coordinates",
                detail: error.to_string(),
            }
        })?;
        let mapped = self.forward(&coordinates)?;
        let residual = max_coordinate_residual(point, &mapped.point);
        if !residual.0.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "surface radial-map inverse residual is not representable at the current precision".to_owned() }.into());
        }
        // Use the actual supplied radius and angles for its density; forward
        // reconstruction above is only an independent round-trip diagnostic.
        let jacobian = radial_jacobian * angular_jacobian * radius.powi(self.dimension as i32 - 1);
        let inverse_jacobian = jacobian.inv();
        if [&jacobian, &inverse_jacobian]
            .iter()
            .any(|value| !value.0.is_finite() || **value <= value.zero())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "surface radial inverse density",
                detail: "positive Jacobian pair is not representable at the supplied momentum"
                    .to_owned(),
            }
            .into());
        }
        Ok(SurfaceRadialPoint {
            coordinates,
            point: point.to_vec(),
            radius,
            jacobian,
            inverse_jacobian,
            residual,
        })
    }

    fn validate_coordinates(&self, coordinates: &[F<T>]) -> Result<()> {
        if coordinates.len() != self.dimension {
            return Err(eyre!(
                "surface radial-map forward received dimension {}, expected {}",
                coordinates.len(),
                self.dimension
            ));
        }
        let zero = coordinates[0].zero();
        let one = zero.one();
        if coordinates
            .iter()
            .any(|coordinate| coordinate.is_nan() || coordinate <= &zero || coordinate >= &one)
        {
            return Err(eyre!(
                "surface radial-map coordinates must be finite and strictly inside the unit cube"
            ));
        }
        Ok(())
    }

    fn radial_parameters(&self) -> (F<T>, F<T>, F<T>) {
        (
            F(self
                .threshold_radius
                .clone()
                .unwrap_or_else(|| self.center[0].zero())),
            F::from_f64(self.beta),
            F::from_f64(self.power),
        )
    }

    fn radius_from_coordinate(
        coordinate: &F<T>,
        threshold: F<T>,
        beta: F<T>,
        power: F<T>,
    ) -> (F<T>, F<T>) {
        let one = coordinate.one();
        let split = &threshold / (&threshold + &beta);
        if split > coordinate.zero() && coordinate < &split {
            let fraction = coordinate / &split;
            let radius = &threshold * Self::power_complement(&fraction, &power);
            let jacobian = &threshold * &power * (&one - &fraction).powf(&(&power - &one)) / &split;
            (radius, jacobian)
        } else {
            let odds = (coordinate - &split) / (&one - coordinate);
            let radius = &threshold + &beta * odds.powf(&power);
            // d[(u-s)/(1-u)]/du = (1-s)/(1-u)^2; the branch
            // interval shrinks when the threshold occupies a finite radius.
            let jacobian = &beta * &power * odds.powf(&(&power - &one)) * (&one - &split)
                / (&one - coordinate).square();
            (radius, jacobian)
        }
    }

    /// Inverse coordinate and positive dr/du evaluated at the supplied radius.
    /// The density depends on its signed distance, not the rounded recovered u.
    fn coordinate_from_radius(
        radius: &F<T>,
        threshold: F<T>,
        beta: F<T>,
        power: F<T>,
    ) -> (F<T>, F<T>) {
        let one = radius.one();
        let split = &threshold / (&threshold + &beta);
        let alpha = &one - power.inv();
        if split > radius.zero() && radius < &threshold {
            let coordinate =
                &split * Self::power_complement(&(radius / &threshold), &(&one / &power));
            let jacobian =
                &power * (&threshold + &beta) * ((&threshold - radius) / &threshold).powf(&alpha);
            (coordinate, jacobian)
        } else {
            let distance = (radius - &threshold) / &beta;
            let z = distance.powf(&(&one / &power));
            let coordinate = (&z + &split) / (&one + &z);
            let jacobian =
                &power * (&threshold + &beta) * distance.powf(&alpha) * (&one + &z).square();
            (coordinate, jacobian)
        }
    }

    /// Evaluate `1 - (1 - fraction)^power` without cancellation at the origin.
    /// The existing native hyperbolic operations preserve small arguments;
    /// direct subtraction is safe once its result is bounded away from zero.
    /// These native operations are not yet registered as eager primitives.
    fn power_complement(fraction: &F<T>, power: &F<T>) -> F<T> {
        let one = fraction.one();
        if power == &one {
            return fraction.clone();
        }
        let two = fraction.from_i64(2);
        let half_epsilon = fraction.epsilon() / &two;
        let logarithm = if fraction < &half_epsilon {
            // The first-order limit is already accurate to native precision
            // and avoids underflow from halving the smallest subnormals.
            -fraction
        } else if fraction < &(&one / &two) {
            let argument = fraction / (&two - fraction);
            -&two * F(argument.0.atanh())
        } else {
            (&one - fraction).ln()
        };
        let exponent = power * logarithm;
        if -&exponent < half_epsilon {
            -exponent
        } else if exponent > -&one {
            let half = &exponent / &two;
            -&two * F(half.0.exp()) * F(half.0.sinh())
        } else {
            &one - F(exponent.0.exp())
        }
    }
}

impl SamplingMapKernel {
    /// Construct a kernel for one ordinary global map.
    ///
    /// Composite surface maps are intentionally rejected here: they need a
    /// prepared graph context and are implemented by the process-level map
    /// compiler.  This kernel covers the reusable LMB shell/cartesian base.
    pub fn new(
        definition: SamplingMapDefinition,
        settings: ParameterizationSettings,
        e_cm: f64,
        n_loop_momenta: usize,
    ) -> Result<Self> {
        if !e_cm.is_finite() || e_cm <= 0.0 {
            return Err(eyre!(
                "sampling-map centre-of-mass energy must be positive and finite"
            ));
        }
        if n_loop_momenta == 0 {
            return Err(eyre!(
                "sampling-map kernel requires at least one loop momentum"
            ));
        }
        if !settings.b.is_finite() || settings.b <= 0.0 {
            return Err(eyre!("sampling-map scale b must be positive and finite"));
        }
        if !settings.power.is_finite() || settings.power <= 0.0 {
            return Err(eyre!(
                "sampling-map radial power must be positive and finite"
            ));
        }
        if matches!(
            &settings.mode,
            ParameterizationMode::HyperSphericalFlat
                | ParameterizationMode::SphericalProductCommonRadial
                | ParameterizationMode::MomentumSpace
        ) {
            return Err(eyre!(
                "sampling-map kernel requires a bijective ordinary map; {:?} is not invertible",
                settings.mode
            ));
        }
        if matches!(&settings.mode, ParameterizationMode::Cartesian)
            && matches!(&settings.mapping, ParameterizationMapping::Power)
        {
            return Err(eyre!(
                "power radial mapping is not defined for cartesian sampling"
            ));
        }
        if !matches!(
            &definition,
            SamplingMapDefinition::Lmb(_) | SamplingMapDefinition::Complement(_)
        ) {
            return Err(eyre!(
                "ordinary sampling-map kernel accepts only lmb(...) or complement(...); prepared composite maps require a graph context"
            ));
        }
        if matches!(&definition, SamplingMapDefinition::Lmb(edges) | SamplingMapDefinition::Complement(edges) if edges.is_empty())
        {
            return Err(eyre!(
                "ordinary sampling-map kernel requires a non-empty edge definition"
            ));
        }
        if let SamplingMapDefinition::Lmb(edges) = &definition
            && edges.len() != n_loop_momenta
        {
            return Err(eyre!(
                "lmb definition has {} edges, but the kernel has {} loop coordinates",
                edges.len(),
                n_loop_momenta
            ));
        }
        Ok(Self {
            definition,
            settings,
            e_cm,
            dimensions: 3 * n_loop_momenta,
        })
    }

    pub fn definition(&self) -> &SamplingMapDefinition {
        &self.definition
    }

    pub fn settings(&self) -> &ParameterizationSettings {
        &self.settings
    }

    pub fn dimensions(&self) -> usize {
        self.dimensions
    }

    pub fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: SamplingSupport::Full,
            requires_context: false,
            requires_proposal_policy: false,
            jacobian: SamplingJacobian::ExactForward,
        }
    }

    pub fn forward<T: FloatLike>(&self, coordinates: &[F<T>]) -> Result<SamplingMapPoint<T>> {
        self.validate_coordinates(coordinates)?;
        let e_cm = F::<T>::from_f64(self.e_cm);
        let (raw, jacobian) = global_parameterize(coordinates, e_cm.clone(), &self.settings);
        if !jacobian.0.is_finite()
            || jacobian <= jacobian.zero()
            || raw
                .iter()
                .flatten()
                .any(|component| !component.0.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "sampling-map forward point or positive Jacobian is not representable at the current precision".to_owned() }.into());
        }
        let loop_momenta = LoopMomenta(
            raw.into_iter()
                .map(|p| ThreeMomentum::new(p[0].clone(), p[1].clone(), p[2].clone()))
                .collect(),
        );
        let (inverse_coordinates, inverse_jacobian) =
            global_inv_parameterize(&loop_momenta.0, e_cm, &self.settings);
        self.validate_coordinates(&inverse_coordinates)
            .map_err(|error| SamplingEvaluationError::Unrepresentable {
                operation: "ordinary forward round-trip inverse",
                detail: error.to_string(),
            })?;
        if !inverse_jacobian.0.is_finite() || inverse_jacobian <= inverse_jacobian.zero() {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "sampling-map inverse Jacobian must be finite and positive at the current precision".to_owned() }.into());
        }
        let residual = max_coordinate_residual(coordinates, &inverse_coordinates);
        Ok(SamplingMapPoint {
            coordinates: coordinates.to_vec(),
            loop_momenta,
            jacobian,
            inverse_jacobian,
            residual,
        })
    }

    pub fn inverse<T: FloatLike>(
        &self,
        loop_momenta: &LoopMomenta<F<T>>,
    ) -> Result<SamplingMapPoint<T>> {
        if loop_momenta.0.len() * 3 != self.dimensions {
            return Err(eyre!(
                "sampling-map inverse received {} loop momenta, expected {}",
                loop_momenta.0.len(),
                self.dimensions / 3
            ));
        }
        if loop_momenta
            .0
            .iter()
            .flat_map(|momentum| [&momentum.px, &momentum.py, &momentum.pz])
            .any(|component| !component.0.is_finite())
        {
            return Err(eyre!("sampling-map inverse requires finite momenta"));
        }
        let e_cm = F::<T>::from_f64(self.e_cm);
        let (coordinates, inverse_jacobian) =
            global_inv_parameterize(&loop_momenta.0, e_cm.clone(), &self.settings);
        self.validate_coordinates(&coordinates).map_err(|error| {
            SamplingEvaluationError::Unrepresentable {
                operation: "ordinary inverse",
                detail: error.to_string(),
            }
        })?;
        if !inverse_jacobian.0.is_finite() || inverse_jacobian <= inverse_jacobian.zero() {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "sampling-map inverse Jacobian must be finite and positive at the current precision".to_owned() }.into());
        }
        let (mapped, jacobian) = global_parameterize(&coordinates, e_cm, &self.settings);
        let residual = max_momentum_residual(&loop_momenta.0, &mapped);
        if !jacobian.0.is_finite()
            || jacobian <= jacobian.zero()
            || mapped
                .iter()
                .flatten()
                .any(|component| !component.0.is_finite())
            || !residual.0.is_finite()
        {
            return Err(SamplingEvaluationError::Unrepresentable { operation: "radial map", detail: "sampling-map inverse point, positive Jacobian or residual is not representable at the current precision".to_owned() }.into());
        }
        Ok(SamplingMapPoint {
            coordinates,
            loop_momenta: loop_momenta.clone(),
            jacobian,
            inverse_jacobian,
            residual,
        })
    }

    fn validate_coordinates<T: FloatLike>(&self, coordinates: &[F<T>]) -> Result<()> {
        if coordinates.len() != self.dimensions {
            return Err(eyre!(
                "sampling-map received {} coordinates, expected {}",
                coordinates.len(),
                self.dimensions
            ));
        }
        let zero = coordinates[0].zero();
        let one = zero.one();
        if coordinates
            .iter()
            .any(|coordinate| coordinate.is_nan() || coordinate <= &zero || coordinate >= &one)
        {
            return Err(eyre!(
                "sampling-map coordinates must be finite and strictly inside the unit cube at the current precision"
            ));
        }
        Ok(())
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for SamplingMapKernel {
    fn dimensions(&self) -> usize {
        self.dimensions()
    }

    fn output_dimensions(&self) -> usize {
        self.dimensions()
    }

    fn contract(&self) -> SamplingMapContract {
        self.contract()
    }

    fn name(&self) -> &'static str {
        "lmb"
    }

    fn forward(
        &self,
        coordinates: &[T],
        _context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        let coordinates = coordinates.iter().cloned().map(F).collect::<Vec<_>>();
        let evaluation = SamplingMapKernel::forward(self, &coordinates)?;
        let point = evaluation
            .loop_momenta
            .0
            .iter()
            .flat_map(|momentum| {
                [
                    momentum.px.0.clone(),
                    momentum.py.0.clone(),
                    momentum.pz.0.clone(),
                ]
            })
            .collect();
        Ok(SamplingMapEvaluation {
            coordinates: evaluation
                .coordinates
                .into_iter()
                .map(|value| value.0)
                .collect(),
            point,
            jacobian: evaluation.jacobian.0,
            inverse_jacobian: evaluation.inverse_jacobian.0,
            residual: evaluation.residual.0,
            support: self.contract().support,
            diagnostics: Vec::new(),
        })
    }

    fn inverse(
        &self,
        point: &[T],
        _context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        if point.len() != self.dimensions() {
            return Err(eyre!(
                "lmb sampling-map inverse received output dimension {}, expected {}",
                point.len(),
                self.dimensions()
            ));
        }
        let loop_momenta = LoopMomenta(
            point
                .as_chunks::<3>()
                .0
                .iter()
                .map(|components| {
                    ThreeMomentum::new(
                        F(components[0].clone()),
                        F(components[1].clone()),
                        F(components[2].clone()),
                    )
                })
                .collect(),
        );
        let evaluation = SamplingMapKernel::inverse(self, &loop_momenta)?;
        let mapped_point = evaluation
            .loop_momenta
            .0
            .iter()
            .flat_map(|momentum| {
                [
                    momentum.px.0.clone(),
                    momentum.py.0.clone(),
                    momentum.pz.0.clone(),
                ]
            })
            .collect();
        Ok(Some(SamplingMapEvaluation {
            coordinates: evaluation
                .coordinates
                .into_iter()
                .map(|value| value.0)
                .collect(),
            point: mapped_point,
            jacobian: evaluation.jacobian.0,
            inverse_jacobian: evaluation.inverse_jacobian.0,
            residual: evaluation.residual.0,
            support: self.contract().support,
            diagnostics: Vec::new(),
        }))
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for SurfaceRadialMap<T> {
    fn dimensions(&self) -> usize {
        self.dimension()
    }

    fn output_dimensions(&self) -> usize {
        self.dimension()
    }

    fn contract(&self) -> SamplingMapContract {
        self.contract()
    }

    fn name(&self) -> &'static str {
        "surface"
    }

    fn forward(
        &self,
        coordinates: &[T],
        _context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        let coordinates = coordinates.iter().cloned().map(F).collect::<Vec<_>>();
        let evaluation = SurfaceRadialMap::forward(self, &coordinates)?;
        Ok(SamplingMapEvaluation {
            coordinates: evaluation
                .coordinates
                .into_iter()
                .map(|value| value.0)
                .collect(),
            point: evaluation.point.into_iter().map(|value| value.0).collect(),
            jacobian: evaluation.jacobian.0,
            inverse_jacobian: evaluation.inverse_jacobian.0,
            residual: evaluation.residual.0,
            support: self.contract().support,
            diagnostics: Vec::new(),
        })
    }

    fn inverse(
        &self,
        point: &[T],
        _context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        let point_f = point.iter().cloned().map(F).collect::<Vec<_>>();
        let evaluation = SurfaceRadialMap::inverse(self, &point_f)?;
        Ok(Some(SamplingMapEvaluation {
            coordinates: evaluation
                .coordinates
                .into_iter()
                .map(|value| value.0)
                .collect(),
            point: evaluation.point.into_iter().map(|value| value.0).collect(),
            jacobian: evaluation.jacobian.0,
            inverse_jacobian: evaluation.inverse_jacobian.0,
            residual: evaluation.residual.0,
            support: self.contract().support,
            diagnostics: Vec::new(),
        }))
    }
}

impl<T: FloatLike> SamplingMapComponent<T> for ImplicitSurfaceRadialMap<T> {
    fn dimensions(&self) -> usize {
        self.dimension
    }

    fn output_dimensions(&self) -> usize {
        self.dimension
    }

    fn contract(&self) -> SamplingMapContract {
        self.contract()
    }

    fn name(&self) -> &'static str {
        "implicit_surface"
    }

    fn forward(
        &self,
        coordinates: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<SamplingMapEvaluation<T>> {
        ImplicitSurfaceRadialMap::forward_with_context(self, coordinates, context.previous)
    }

    fn inverse(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<SamplingMapEvaluation<T>>> {
        ImplicitSurfaceRadialMap::inverse_with_context(self, point, context.previous).map(Some)
    }

    fn inverse_density(
        &self,
        point: &[T],
        context: &mut SamplingMapContext<'_, T>,
    ) -> Result<Option<T>> {
        self.inverse_evaluation(point, context.previous, false)
            .map(|evaluation| Some(evaluation.inverse_jacobian))
    }
}

fn max_coordinate_residual<T: FloatLike>(a: &[F<T>], b: &[F<T>]) -> F<T> {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .max_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or_else(|| F::<T>::from_f64(0.0))
}

fn max_coordinate_residual_scalar<T: FloatLike>(a: &[T], b: &[T]) -> T {
    a.iter()
        .zip(b)
        .map(|(x, y)| (F(x.clone()) - F(y.clone())).abs())
        .fold(F(a[0].zero()), F::max)
        .0
}

fn max_momentum_residual<T: FloatLike>(a: &[ThreeMomentum<F<T>>], b: &[[F<T>; 3]]) -> F<T> {
    a.iter()
        .zip(b)
        .flat_map(|(x, y)| {
            [
                (x.px.clone() - &y[0]).abs(),
                (x.py.clone() - &y[1]).abs(),
                (x.pz.clone() - &y[2]).abs(),
            ]
        })
        .max_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or_else(|| F::<T>::from_f64(0.0))
}

fn parse_edge(atom: AtomView<'_>, constructor: &str) -> Result<usize> {
    let value = i64::try_from(atom)
        .map_err(|_| eyre!("{constructor} expects non-negative integer edge IDs, got `{atom}`"))?;
    usize::try_from(value)
        .map_err(|_| eyre!("{constructor} expects non-negative integer edge IDs, got {value}"))
}

fn call<I>(name: &str, arguments: I) -> Atom
where
    I: IntoIterator<Item = Atom>,
{
    symbol!(&format!("gammalooprs::sampling_map::{name}")).call_args(arguments)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lu_h_proposal_radius_stops_at_declared_accuracy_without_changing_density() {
        std::thread::Builder::new()
            .name("lu-h-proposal-root-test".to_owned())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                use crate::settings::runtime::{HFunctionSettings, SamplingRadialProfile};
                use std::sync::atomic::{AtomicUsize, Ordering};

                crate::initialisation::test_initialise().unwrap();
                let profile = SamplingRadialProfile {
                    broad_fraction: 0.0,
                    ..Default::default()
                };
                let program = SamplingExpressionEvaluator::new_lu_h_profile(
                    &profile,
                    &HFunctionSettings::default(),
                )
                .unwrap();
                fn check<T: FloatLike>(
                    program: &SamplingExpressionEvaluator,
                    profile: &SamplingRadialProfile,
                ) {
                    let one = F::<T>::default().one();
                    let calls = Arc::new(AtomicUsize::new(0));
                    let count = Arc::clone(&calls);
                    let shell = ImplicitSurfaceRadialMap::new(
                        3,
                        vec![one.zero().0; 3],
                        1.0,
                        2.0,
                        Arc::new(move |_, radius: T| {
                            count.fetch_add(1, Ordering::Relaxed);
                            let radius = F(radius);
                            Ok((
                                (radius.square() - radius.from_i64(2)).0,
                                (&radius * radius.from_i64(2)).0,
                            ))
                        }),
                    )
                    .unwrap()
                    .with_root_tolerance(1e-9)
                    .unwrap();
                    let proposal = shell
                        .clone()
                        .with_lu_h_profile(program.clone(), profile, 1)
                        .unwrap();
                    let direction = [one.0.clone(), one.zero().0, one.zero().0];
                    let (root, _) = proposal
                        .root_for_direction(&direction, proposal.center(), &[], None)
                        .unwrap()
                        .unwrap();
                    let proposal_calls = calls.swap(0, Ordering::Relaxed);
                    let residual = (F(root).square() - one.from_i64(2)).abs();
                    // This witness deliberately uses an inexact focusing radius;
                    // the actual map inverse must still satisfy the strict law.
                    assert!(residual < F::<T>::from_f64(2e-9));
                    assert!(residual > one.epsilon() * one.from_i64(1024));
                    let (strict_root, _) = shell
                        .root_for_direction(&direction, shell.center(), &[], None)
                        .unwrap()
                        .unwrap();
                    assert!(calls.load(Ordering::Relaxed) > proposal_calls);
                    assert!(
                        (F(strict_root).square() - one.from_i64(2)).abs()
                            <= one.epsilon() * one.from_i64(128)
                    );
                    for radial in [0.01, 0.3, 0.9] {
                        let coordinates =
                            [radial, 0.31, 0.64].map(|value| F::<T>::from_f64(value).0);
                        let forward = proposal.forward(&coordinates).unwrap();
                        let inverse = proposal.inverse(&forward.point).unwrap();
                        assert!(
                            (F(forward.jacobian) * F(inverse.inverse_jacobian) - &one).abs()
                                < F::<T>::from_f64(1e-13)
                        );
                    }
                }
                check::<crate::utils::SamplingFloat>(&program, &profile);
                check::<crate::utils::ArbPrec>(&program, &profile);
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn lu_h_radial_profile_has_native_exact_density_and_inverse() {
        // Eager/dual compilation needs a larger development-build stack, as
        // in the generated graph fixtures; this leaves all numeric checks intact.
        std::thread::Builder::new()
            .name("lu-h-profile-test".to_owned())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                crate::initialisation::test_initialise().unwrap();

                use crate::settings::runtime::{HFunctionSettings, SamplingRadialProfile};
                use crate::utils::{ArbPrec, QuadFloat};

                fn check<T: FloatLike>(
                    program: &SamplingExpressionEvaluator,
                    profile: &SamplingRadialProfile,
                ) {
                    let zero = F::<T>::from_f64(0.0);
                    let map = ImplicitSurfaceRadialMap::new(
                        3,
                        vec![zero.0.clone(); 3],
                        1.0,
                        2.0,
                        Arc::new(|_, radius: T| {
                            let r = F(radius);
                            Ok(((&r - r.from_i64(2)).0, r.one().0))
                        }),
                    )
                    .unwrap()
                    .with_lu_h_profile(program.clone(), profile, 2)
                    .unwrap();
                    assert_eq!(map.lu_h_max_occurrence(), Some(2));
                    for radial in [1e-8, 0.07, 0.5, 0.91, 1.0 - 1e-8] {
                        let coordinates =
                            [radial, 0.31, 0.64].map(|value| F::<T>::from_f64(value).0);
                        let forward = map.forward(&coordinates).unwrap();
                        let inverse = map.inverse(&forward.point).unwrap();
                        assert_eq!(
                            map.inverse_density(
                                &forward.point,
                                &mut SamplingMapContext::detached(&[]),
                            )
                            .unwrap(),
                            Some(inverse.inverse_jacobian.clone()),
                        );
                        let tolerance = F::<T>::from_f64(1e-7);
                        assert!(
                            (F(forward.jacobian.clone()) * F(inverse.inverse_jacobian)
                                - zero.one())
                            .abs()
                                < tolerance
                        );
                        assert!(
                            (F(inverse.coordinates[0].clone()) - F(coordinates[0].clone())).abs()
                                < F::<T>::from_f64(1e-12)
                        );
                        let radius = forward
                            .point
                            .iter()
                            .map(|value| F(value.clone()).square())
                            .fold(zero.clone(), |sum, term| sum + term)
                            .sqrt();
                        let t = radius.from_i64(2) / &radius;
                        let fit = map.lu_h_profile.as_ref().unwrap();
                        let z = (t.ln() - &fit.log_scale).0.exp();
                        let z = F(z).powf(&fit.shape);
                        let focused = &fit.shape * &z / (&t * (zero.one() + &z).square());
                        let broad = zero.from_i64(2) / (zero.from_i64(2) + &t).square();
                        let epsilon = F::<T>::from_f64(profile.broad_fraction);
                        let density = (zero.one() - &epsilon) * focused + epsilon * broad;
                        let expected =
                            zero.from_i64(2) * zero.TAU() * radius.powi(3) / (&t * density);
                        assert!((F(forward.jacobian) / expected - zero.one()).abs() < tolerance);
                    }
                }
                for broad_fraction in [0.0, 0.02, 1.0] {
                    let profile = SamplingRadialProfile {
                        broad_fraction,
                        ..Default::default()
                    };
                    let program = SamplingExpressionEvaluator::new_lu_h_profile(
                        &profile,
                        &HFunctionSettings::default(),
                    )
                    .unwrap();
                    check::<f64>(&program, &profile);
                    check::<QuadFloat>(&program, &profile);
                    check::<ArbPrec>(&program, &profile);
                }
                // The two component quantiles coincide everywhere. A strict-sign
                // scalar bracket is unnecessary and would reject this exact solution.
                let profile = SamplingRadialProfile {
                    scale: Some(2.0),
                    shape: Some(1.0),
                    ..Default::default()
                };
                let program = SamplingExpressionEvaluator::new_lu_h_profile(
                    &profile,
                    &HFunctionSettings::default(),
                )
                .unwrap();
                check::<f64>(&program, &profile);
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn lu_h_radial_profile_rejects_incompatible_derivative_columns() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let parameters = ["y", "r", "b", "f", "s"]
            .map(|name| Atom::var(symbol!(&format!("lu_column_contract::{name}"))));
        for columns in [&[][..], &[1], &[1, 0, 2, 3, 4]] {
            // The fit itself is finite with positive shape. Reject missing,
            // reordered or extra derivatives before derivative[0] can mean
            // another input or carry unused runtime work.
            let program = SamplingExpressionEvaluator::new(
                [
                    "0",
                    "1",
                    "lu_column_contract::y",
                    "1-lu_column_contract::y",
                    "-lu_column_contract::y",
                ]
                .map(|expression| try_parse!(expression).unwrap()),
                parameters.clone(),
                columns,
            )?;
            let error = ImplicitSurfaceRadialMap::new(
                3,
                vec![0.0; 3],
                1.0,
                1.0,
                Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
            )?
            .with_lu_h_profile(program, &SamplingRadialProfile::default(), 1)
            .unwrap_err();
            assert!(error.to_string().contains("active derivative column"));
        }
        Ok(())
    }

    #[test]
    fn lu_h_radial_profile_uses_root_independent_broad_and_absent_laws() {
        // Eager/dual compilation needs a larger development-build stack, as
        // in the generated graph fixtures; this leaves all numeric checks intact.
        std::thread::Builder::new()
            .name("lu-h-profile-test".to_owned())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                crate::initialisation::test_initialise().unwrap();

                use crate::settings::runtime::{HFunctionSettings, SamplingRadialProfile};
                for broad_only in [false, true] {
                    let profile = SamplingRadialProfile {
                        broad_fraction: if broad_only { 1.0 } else { 0.02 },
                        ..Default::default()
                    };
                    let program = SamplingExpressionEvaluator::new_lu_h_profile(
                        &profile,
                        &HFunctionSettings::default(),
                    )
                    .unwrap();
                    let map = ImplicitSurfaceRadialMap::new(
                        3,
                        vec![0.0; 3],
                        3.0,
                        2.5,
                        Arc::new(move |_, _| {
                            assert!(
                                !broad_only,
                                "broad-only proposal must never evaluate its cut root"
                            );
                            Ok((1.0, 0.0))
                        }),
                    )
                    .unwrap()
                    .with_lu_h_profile(program, &profile, 1)
                    .unwrap();
                    for u in [0.1, 0.5, 0.9] {
                        let forward = map.forward(&[u, 0.31, 0.64]).unwrap();
                        let radius = forward
                            .point
                            .iter()
                            .map(|value| value * value)
                            .sum::<f64>()
                            .sqrt();
                        assert!((radius / (3.0 * (1.0 - u) / u) - 1.0).abs() < 1e-13);
                        let expected_j =
                            4.0 * std::f64::consts::PI * radius * radius * 3.0 / (u * u);
                        assert!((forward.jacobian / expected_j - 1.0).abs() < 1e-13);
                        assert!(
                            (map.inverse(&forward.point).unwrap().coordinates[0] - u).abs() < 1e-14
                        );
                    }
                }
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn lu_h_radial_profile_jacobian_matches_full_finite_difference() {
        // Eager/dual compilation needs a larger development-build stack, as
        // in the generated graph fixtures; this leaves all numeric checks intact.
        std::thread::Builder::new()
            .name("lu-h-profile-test".to_owned())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                crate::initialisation::test_initialise().unwrap();

                use crate::settings::runtime::{HFunctionSettings, SamplingRadialProfile};
                let profile = SamplingRadialProfile::default();
                let program = SamplingExpressionEvaluator::new_lu_h_profile(
                    &profile,
                    &HFunctionSettings::default(),
                )
                .unwrap();
                let map = ImplicitSurfaceRadialMap::new(
                    3,
                    vec![0.0; 3],
                    1.7,
                    2.0,
                    Arc::new(|direction, radius| Ok((radius - (2.0 + 0.4 * direction[0]), 1.0))),
                )
                .unwrap()
                .with_lu_h_profile(program, &profile, 1)
                .unwrap();
                for u in [0.15, 0.5, 0.87] {
                    let x = [u, 0.31, 0.64];
                    let step = 1e-5;
                    let mut jac = [[0.0; 3]; 3];
                    for col in 0..3 {
                        let mut plus = x;
                        plus[col] += step;
                        let mut minus = x;
                        minus[col] -= step;
                        let a = map.forward(&plus).unwrap();
                        let b = map.forward(&minus).unwrap();
                        for (row, entries) in jac.iter_mut().enumerate() {
                            entries[col] = (a.point[row] - b.point[row]) / (2.0 * step);
                        }
                    }
                    let determinant = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
                        - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
                        + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
                    assert!(
                        (map.forward(&x).unwrap().jacobian / determinant.abs() - 1.0).abs() < 2e-7
                    );
                }
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn radial_profiles_include_the_outer_branch_interval_in_the_jacobian() {
        for threshold in [None, Some(2.0)] {
            for power in [1.0, 1.5, 2.0, 3.0] {
                let explicit = |coordinate| {
                    SurfaceRadialMap::radius_from_coordinate(
                        &F(coordinate),
                        F(threshold.unwrap_or(0.0)),
                        F(2.0),
                        F(power),
                    )
                };
                let implicit = ImplicitSurfaceRadialMap::new(
                    3,
                    vec![0.0; 3],
                    2.0,
                    power,
                    Arc::new(|_: &[f64], radius: f64| Ok((radius - 2.0, 1.0))),
                )
                .unwrap();
                for coordinate in [0.1, 0.8] {
                    let step = 1.0e-6;
                    let (radius, jacobian) = explicit(coordinate);
                    let plus = explicit(coordinate + step).0;
                    let minus = explicit(coordinate - step).0;
                    let finite_difference = (plus.0 - minus.0) / (2.0 * step);
                    assert!((jacobian.0 / finite_difference - 1.0).abs() < 1.0e-8);
                    let (implicit_radius, implicit_jacobian) = implicit
                        .radius_from_coordinate(coordinate, threshold.map(|radius| (radius, 1.0)))
                        .unwrap();
                    assert!((implicit_radius - radius.0).abs() < 1.0e-12);
                    assert!((implicit_jacobian / finite_difference - 1.0).abs() < 1.0e-8);

                    let map =
                        SurfaceRadialMap::new(6, vec![0.1; 6], threshold, 2.0, power).unwrap();
                    assert_eq!(map.contract().jacobian, SamplingJacobian::ExactForward);
                    let coordinates = [coordinate, 0.27, 0.61, 0.39, 0.72, 0.58].map(F);
                    let mapped = map.forward(&coordinates).unwrap();
                    let mut matrix = vec![vec![0.0; 6]; 6];
                    for axis in 0..6 {
                        let mut plus = coordinates;
                        let mut minus = coordinates;
                        plus[axis].0 += step;
                        minus[axis].0 -= step;
                        let plus = map.forward(&plus).unwrap();
                        let minus = map.forward(&minus).unwrap();
                        for (component, row) in matrix.iter_mut().enumerate() {
                            row[axis] =
                                (plus.point[component].0 - minus.point[component].0) / (2.0 * step);
                        }
                    }
                    let determinant = SamplingMapAffine::new(matrix, vec![0.0; 6])
                        .unwrap()
                        .determinant();
                    assert!((determinant / mapped.jacobian.0 - 1.0).abs() < 1.0e-4);
                    let inverse = map.inverse(&mapped.point).unwrap();
                    assert!(
                        (mapped.jacobian * inverse.inverse_jacobian - F(1.0)).abs() < F(1.0e-10)
                    );
                }
            }
        }
    }

    #[test]
    fn surface_radial_profiles_focus_both_signed_normal_limits() {
        for power in [1.5_f64, 2.0, 3.0] {
            let alpha = 1.0 - 1.0 / power;
            let implicit = ImplicitSurfaceRadialMap::new(
                3,
                vec![0.0; 3],
                2.0,
                power,
                Arc::new(|_, radius| Ok((radius - 3.0, 1.0))),
            )
            .unwrap();
            for sign in [-1.0, 1.0] {
                let (scale, interval) = if sign < 0.0 {
                    (3.0_f64, 0.6)
                } else {
                    (2.0_f64, 0.4)
                };
                let limiting_weight = power * scale.powf(1.0 / power) / interval;
                for distance in [1.0e-3, 3.0e-4, 1.0e-4] {
                    let coordinate = 0.6 + sign * distance;
                    let (radius, jacobian) = SurfaceRadialMap::radius_from_coordinate(
                        &F(coordinate),
                        F(3.0),
                        F(2.0),
                        F(power),
                    );
                    let delta = radius.0 - 3.0;
                    assert!(delta * sign > 0.0);
                    // A |delta|^-alpha normal singularity has a bounded
                    // radial weight on BOTH sides when alpha = 1 - 1/p.
                    let normal_weight = jacobian.0 / delta.abs().powf(alpha);
                    assert!(
                        (normal_weight / limiting_weight - 1.0).abs() < 0.02,
                        "p={power}, sign={sign}, u={coordinate}, normal weight={normal_weight}, limit={limiting_weight}",
                    );
                    let mapped = implicit.forward(&[coordinate, 0.27, 0.61]).unwrap();
                    let radial_jacobian =
                        mapped.jacobian / (4.0 * std::f64::consts::PI * radius.0.powi(2));
                    assert!((radial_jacobian / jacobian.0 - 1.0).abs() < 1.0e-12);
                }
            }
        }
    }

    #[test]
    fn surface_radial_profiles_match_eager_branch_derivatives() {
        use crate::integrands::process::sampling_evaluator::SamplingExpressionEvaluator;

        crate::initialisation::test_initialise().unwrap();
        let coordinate = try_parse!("radial_profile_test::u").unwrap();
        let power = try_parse!("radial_profile_test::p").unwrap();
        // Eager algebra is an independent branch oracle at ordinary points.
        // The production owner uses stable native primitives near endpoints.
        for (expression, coordinates) in [
            (
                "3*(1-(1-radial_profile_test::u/(3/5))^radial_profile_test::p)",
                [0.13, 0.43],
            ),
            (
                "3+2*((radial_profile_test::u-3/5)/(1-radial_profile_test::u))^radial_profile_test::p",
                [0.71, 0.89],
            ),
        ] {
            let mut evaluator = SamplingExpressionEvaluator::new(
                [try_parse!(expression).unwrap()],
                [coordinate.clone(), power.clone()],
                &[0, 1],
            )
            .unwrap();
            for power in [1.0, 1.5, 2.0, 3.0] {
                for coordinate in coordinates {
                    let (radius, jacobian) = SurfaceRadialMap::radius_from_coordinate(
                        &F(coordinate),
                        F(3.0),
                        F(2.0),
                        F(power),
                    );
                    let result = evaluator
                        .evaluate_with_derivatives(&[coordinate, power])
                        .unwrap();
                    assert!((result[0].value.re.0 / radius.0 - 1.0).abs() < 1.0e-12);
                    assert!((result[0].derivatives[0].re.0 / jacobian.0 - 1.0).abs() < 1.0e-12);
                }
            }
        }
    }

    #[test]
    fn native_composition_rejects_combined_overflow_and_recovers() {
        fn evaluate<T: FloatLike>() -> Result<SamplingMapEvaluation<T>> {
            let one = F::<T>::default().one();
            let scale = one.from_i64(10).powi(160);
            let make = || SamplingMapAffine::new(vec![vec![scale.0.clone()]], vec![one.zero().0]);
            let first = make()?;
            let second = make()?;
            let coordinate = (&one / one.from_i64(4)).0;
            assert!(
                first
                    .forward(
                        std::slice::from_ref(&coordinate),
                        &mut SamplingMapContext::detached(&[])
                    )?
                    .jacobian
                    .is_finite()
            );
            assert!(
                second
                    .forward(
                        std::slice::from_ref(&coordinate),
                        &mut SamplingMapContext::detached(&[])
                    )?
                    .jacobian
                    .is_finite()
            );
            SamplingMapComposition::product(vec![Box::new(first), Box::new(second)])?.forward(
                &[coordinate.clone(), coordinate],
                &mut SamplingMapContext::detached(&[]),
            )
        }
        let error = evaluate::<f64>().unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());
        // QuadFloat extends the mantissa with two f64 components; its exponent
        // range remains binary64, so this range failure requires Arb precision.
        let quad_error = evaluate::<crate::utils::QuadFloat>().unwrap_err();
        assert!(
            quad_error
                .downcast_ref::<SamplingEvaluationError>()
                .is_some()
        );
        let arb = evaluate::<crate::utils::ArbPrec>().unwrap();
        assert!(!F(arb.jacobian.clone()).is_infinite());
        let one = F(arb.jacobian.clone()).one();
        assert!(
            (F(arb.jacobian) * F(arb.inverse_jacobian) - &one).abs()
                < one.epsilon() * one.from_i64(64)
        );
    }

    #[test]
    fn inverse_density_rejects_composed_jacobian_range_failures() {
        for scale in [1.0e-160, 1.0e160] {
            let child = SamplingMapAffine::new(vec![vec![scale]], vec![0.0]).unwrap();
            let map =
                SamplingMapComposition::product(vec![Box::new(child.clone()), Box::new(child)])
                    .unwrap();
            let point = [scale / 4.0; 2];
            let inverse = map
                .inverse(&point, &mut SamplingMapContext::detached(&[]))
                .unwrap_err();
            let density = map
                .inverse_density(&point, &mut SamplingMapContext::detached(&[]))
                .unwrap_err();
            assert!(inverse.downcast_ref::<SamplingEvaluationError>().is_some());
            assert!(density.downcast_ref::<SamplingEvaluationError>().is_some());
            for nonfinite in [f64::NAN, f64::INFINITY] {
                let point = [nonfinite, scale / 4.0];
                let inverse = map
                    .inverse(&point, &mut SamplingMapContext::detached(&[]))
                    .unwrap_err();
                let density = map
                    .inverse_density(&point, &mut SamplingMapContext::detached(&[]))
                    .unwrap_err();
                assert_eq!(density.to_string(), inverse.to_string());
            }
        }
    }

    #[test]
    fn surface_radial_profiles_preserve_tiny_native_coordinates() {
        fn check<T: FloatLike>(decimal_exponent: i32) {
            let one = F::<T>::default().one();
            let tiny = &one / one.from_usize(10).powi(decimal_exponent);
            let tolerance = one.epsilon() * one.from_usize(1024);
            for power in [1.0, 2.0, 3.0] {
                let map = SurfaceRadialMap::new(
                    3,
                    vec![one.zero().0; 3],
                    Some(one.from_usize(3).0),
                    2.0,
                    power,
                )
                .unwrap();
                let coordinates = [
                    tiny.clone(),
                    &one / one.from_usize(3),
                    &one / one.from_usize(2),
                ];
                let mapped = map.forward(&coordinates).unwrap();
                let expected = &tiny * one.from_usize(5 * power as usize);
                assert!((&mapped.radius / expected - &one).abs() < tolerance);
                let inverse = map.inverse(&mapped.point).unwrap();
                assert!((&inverse.coordinates[0] / &tiny - &one).abs() < tolerance);
                assert!((&mapped.jacobian * &inverse.inverse_jacobian - &one).abs() < tolerance);
            }
            assert_eq!(SurfaceRadialMap::power_complement(&tiny, &one), tiny);
        }
        check::<f64>(100);
        check::<crate::utils::QuadFloat>(50);
        check::<crate::utils::ArbPrec>(400);
        let subnormal = F(f64::from_bits(1));
        assert_eq!(
            SurfaceRadialMap::power_complement(&subnormal, &F(2.0)),
            F(f64::from_bits(2))
        );
        let implicit = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            2.0,
            2.0,
            Arc::new(|_, radius| Ok((radius - 3.0, 1.0))),
        )
        .unwrap();
        let coordinates = [1.0e-100, 0.27, 0.61];
        let mapped = implicit.forward(&coordinates).unwrap();
        let inverse = implicit.inverse(&mapped.point).unwrap();
        assert!((inverse.coordinates[0] / coordinates[0] - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn surface_radial_profiles_reject_rounded_threshold_seams() {
        let explicit = SurfaceRadialMap::new(3, vec![0.0; 3], Some(3.0), 2.0, 2.0).unwrap();
        let implicit = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            2.0,
            2.0,
            Arc::new(|_, radius| Ok((radius - 3.0, 1.0))),
        )
        .unwrap();
        for coordinate in [0.6 - 1.0e-10, 0.6, 0.6 + 1.0e-10] {
            let error = explicit
                .forward(&[F(coordinate), F(0.27), F(0.61)])
                .unwrap_err();
            assert!(error.to_string().contains("threshold seam"), "{error}");
            let error = implicit.forward(&[coordinate, 0.27, 0.61]).unwrap_err();
            assert!(error.to_string().contains("threshold seam"), "{error}");
        }
        // Power one has no singular threshold seam and remains invertible.
        let regular = SurfaceRadialMap::new(3, vec![0.0; 3], Some(3.0), 2.0, 1.0).unwrap();
        let mapped = regular.forward(&[F(0.6), F(0.27), F(0.61)]).unwrap();
        assert!(regular.inverse(&mapped.point).is_ok());
    }

    #[test]
    fn surface_radial_profiles_preserve_gaussian_normalization_and_second_moment() {
        let count = 4096;
        let gaussian_normalization = (2.0 * std::f64::consts::PI).powf(-1.5);
        for threshold in [None, Some(2.0)] {
            for power in [1.0, 1.5, 2.0, 3.0] {
                let map = SurfaceRadialMap::new(3, vec![0.0; 3], threshold, 2.0, power).unwrap();
                let mut normalization = 0.0;
                let mut second_moment = 0.0;
                for sample in 0..count {
                    // In 3D the uniform-cos(theta) angular determinant is
                    // exactly 4pi. A radial Gaussian's angular integral is
                    // therefore exact at any nonsingular angular point.
                    let coordinate = (sample as f64 + 0.5) / count as f64;
                    let mapped = map.forward(&[F(coordinate), F(0.27), F(0.61)]).unwrap();
                    let radius_squared = mapped
                        .point
                        .iter()
                        .map(|component| component.0.powi(2))
                        .sum::<f64>();
                    let weight =
                        gaussian_normalization * (-0.5 * radius_squared).exp() * mapped.jacobian.0;
                    normalization += weight / count as f64;
                    second_moment += weight * radius_squared / count as f64;
                }
                assert!(
                    (normalization - 1.0).abs() < 1.0e-5,
                    "R={threshold:?}, p={power}, norm={normalization}"
                );
                assert!(
                    (second_moment - 3.0).abs() < 5.0e-5,
                    "R={threshold:?}, p={power}, moment={second_moment}"
                );
            }
        }
    }

    #[test]
    fn affine_map_round_trips_with_translation_and_exact_determinant() {
        let map = SamplingMapAffine::new(vec![vec![2.0, 1.0], vec![1.0, 3.0]], vec![0.5, -1.0])
            .expect("invertible affine map");
        assert_eq!(map.dimension(), 2);
        assert!((map.determinant() - 5.0).abs() < 1.0e-14);
        assert_eq!(map.contract().support, SamplingSupport::Full);
        assert_eq!(map.contract().jacobian, SamplingJacobian::ExactForward);

        let forward = map
            .forward(&[0.2, -0.4], &mut SamplingMapContext::detached(&[]))
            .unwrap();
        assert!(forward.residual < 1.0e-14, "{}", forward.residual);
        assert!((forward.jacobian - 5.0).abs() < 1.0e-14);
        assert!((forward.inverse_jacobian - 0.2).abs() < 1.0e-14);

        let inverse = map.inverse(&forward.point, &[]).unwrap();
        assert!(inverse.residual < 1.0e-14, "{}", inverse.residual);
        for (actual, expected) in inverse.coordinates.iter().zip([0.2, -0.4]) {
            assert!((actual - expected).abs() < 1.0e-14);
        }
    }

    #[test]
    fn affine_map_rejects_malformed_and_singular_matrices() {
        assert!(SamplingMapAffine::<f64>::new(vec![], vec![]).is_err());
        assert!(SamplingMapAffine::new(vec![vec![1.0, 0.0]], vec![0.0]).is_err());
        assert!(SamplingMapAffine::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]], vec![0.0]).is_err());
        assert!(
            SamplingMapAffine::new(vec![vec![1.0, 2.0], vec![2.0, 4.0]], vec![0.0, 0.0],).is_err()
        );
    }

    #[test]
    fn conditional_affine_embedding_uses_native_prerequisites_and_supplied_inverse_point() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        fn check<T: FloatLike>() -> Result<()> {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let physical = ImplicitSurfaceRadialMap::new(
                3,
                vec![zero.0.clone(); 3],
                2.0,
                1.0,
                Arc::new(|_, r: T| {
                    let r = F(r);
                    Ok(((&r - r.from_i64(2)).0, r.one().0))
                }),
            )?
            .with_context_preparer(Arc::new(|prior: &[T]| {
                if prior.len() != 3 {
                    return Err(eyre!("test fiber requires three fixed coordinates"));
                }
                Ok((prior.to_vec(), PreparedSurfaceStatus::existing(None)?))
            }));
            let calls = Arc::new(AtomicUsize::new(0));
            let transform: SamplingMapContextTransform<T> = Arc::new({
                let calls = calls.clone();
                move |context: &mut SamplingMapContext<'_, T>| {
                    let prior = context.previous;
                    calls.fetch_add(1, Ordering::Relaxed);
                    if prior.len() != 3 {
                        return Err(eyre!("missing prior block"));
                    }
                    let x = F(prior[0].clone());
                    let tau = x.one() + x.square();
                    let inverse = tau.one() / &tau;
                    let matrix = (0..3)
                        .map(|i| {
                            (0..3)
                                .map(|j| {
                                    if i == j {
                                        inverse.0.clone()
                                    } else {
                                        inverse.zero().0
                                    }
                                })
                                .collect()
                        })
                        .collect();
                    let shift = (&inverse.one() - &inverse) / inverse.from_i64(3);
                    Ok((
                        prior.iter().map(|x| (F(x.clone()) * &tau).0).collect(),
                        SamplingMapAffine::new(matrix, vec![shift.0; 3])?,
                    ))
                }
            });
            let active = SamplingMapEmbedding::from_composition(
                SamplingMapComposition::then(vec![Box::new(physical.clone())])?,
                (0..3).collect(),
            )?
            .with_context_transform(transform.clone());
            let prior = SurfaceRadialMap::new(3, vec![zero.0.clone(); 3], None, 2.0, 1.0)?;
            let full = SamplingMapEmbedding::from_composition(
                SamplingMapComposition::then(vec![Box::new(prior), Box::new(active.clone())])?,
                (0..6).collect(),
            )?;
            let cube = [17, 31, 63, 29, 43, 71].map(|n| (one.from_i64(n) / one.from_i64(100)).0);
            let mapped = full.forward(&cube, &mut SamplingMapContext::detached(&[]))?;
            assert_eq!(calls.load(Ordering::Relaxed), 1);
            let inverse = full
                .inverse(&mapped.point, &mut SamplingMapContext::detached(&[]))?
                .expect("full-support inverse");
            assert_eq!(calls.load(Ordering::Relaxed), 2);
            let tolerance = one.epsilon().sqrt() * one.from_i64(64);
            assert!(
                (F(mapped.jacobian.clone()) * F(inverse.inverse_jacobian) - &one).abs() < tolerance
            );
            for (actual, expected) in inverse.coordinates.iter().zip(&cube) {
                assert!((F(actual.clone()) - F(expected.clone())).abs() < tolerance);
            }
            // A foreign supplied point is pulled back using its own prior
            // block, not the selected map's scale or a reconstructed point.
            let foreign_prior = mapped.point[..3]
                .iter()
                .map(|x| (F(x.clone()) / one.from_i64(2)).0)
                .collect::<Vec<_>>();
            let foreign = active
                .inverse(
                    &mapped.point[3..],
                    &mut SamplingMapContext::detached(&foreign_prior),
                )?
                .expect("full-support inverse");
            let (fixed, frame) = transform(&mut SamplingMapContext::detached(&foreign_prior))?;
            let pullback = frame.inverse(&mapped.point[3..], &[])?;
            let expected = physical.inverse_with_context(&pullback.coordinates, &fixed)?;
            assert!(
                (F(foreign.inverse_jacobian)
                    / (F(expected.inverse_jacobian) * F(pullback.inverse_jacobian))
                    - &one)
                    .abs()
                    < tolerance
            );
            let step = &one / one.from_i64(1000000);
            let mut matrix = vec![vec![0.0; 6]; 6];
            for column in 0..6 {
                let mut plus = cube.clone();
                let mut minus = cube.clone();
                plus[column] = (F(plus[column].clone()) + &step).0;
                minus[column] = (F(minus[column].clone()) - &step).0;
                let plus = full.forward(&plus, &mut SamplingMapContext::detached(&[]))?;
                let minus = full.forward(&minus, &mut SamplingMapContext::detached(&[]))?;
                for (row, values) in matrix.iter_mut().enumerate() {
                    values[column] = ((F(plus.point[row].clone()) - F(minus.point[row].clone()))
                        / (&step * step.from_i64(2)))
                    .into_f64();
                }
            }
            assert!(
                (SamplingMapAffine::new(matrix, vec![0.0; 6])?.determinant()
                    / F(mapped.jacobian).into_f64()
                    - 1.0)
                    .abs()
                    < 2.0e-6
            );
            Ok(())
        }
        check::<f64>().unwrap();
        check::<crate::utils::QuadFloat>().unwrap();
        check::<crate::utils::ArbPrec>().unwrap();
    }

    #[test]
    fn conditional_affine_embedding_distinguishes_numeric_and_structural_context_failure() {
        let make = || {
            SamplingMapEmbedding::product(
                vec![Box::new(
                    SurfaceRadialMap::new(3, vec![0.0; 3], None, 2.0, 1.0).unwrap(),
                )],
                (0..3).collect(),
            )
            .unwrap()
        };
        let numerical = make().with_context_transform(Arc::new(|_| {
            Ok((
                vec![f64::INFINITY],
                SamplingMapAffine::new(
                    vec![
                        vec![1.0, 0.0, 0.0],
                        vec![0.0, 1.0, 0.0],
                        vec![0.0, 0.0, 1.0],
                    ],
                    vec![0.0; 3],
                )?,
            ))
        }));
        let error = numerical
            .forward(&[0.2, 0.3, 0.6], &mut SamplingMapContext::detached(&[1.0]))
            .unwrap_err();
        assert!(matches!(
            error.downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::Unrepresentable { .. })
        ));
        let structural = make().with_context_transform(Arc::new(|_| {
            Ok((
                vec![1.0],
                SamplingMapAffine::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]], vec![0.0; 2])?,
            ))
        }));
        let error = structural
            .forward(&[0.2, 0.3, 0.6], &mut SamplingMapContext::detached(&[1.0]))
            .unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_none());
        assert!(error.to_string().contains("dimension 3"));
        for (center, numeric) in [(vec![f64::INFINITY; 3], true), (vec![0.0; 2], false)] {
            let map = ImplicitSurfaceRadialMap::new(
                3,
                vec![0.0; 3],
                2.0,
                1.0,
                Arc::new(|_, r| Ok((r - 2.0, 1.0))),
            )
            .unwrap()
            .with_context_preparer(Arc::new(move |_| {
                Ok((center.clone(), PreparedSurfaceStatus::existing(None)?))
            }));
            let error = map.forward(&[0.2, 0.3, 0.6]).unwrap_err();
            assert_eq!(
                error.downcast_ref::<SamplingEvaluationError>().is_some(),
                numeric
            );
        }
    }

    #[test]
    fn parses_host_and_block_descriptors_with_symbolica() {
        let map =
            SamplingMapDefinition::parse("block(lmb(7,3),at_cut(cut(9,2),left(surface(6,4))))")
                .unwrap();
        assert_eq!(
            SamplingMapDefinition::from_atom(map.to_atom().as_view()).unwrap(),
            map
        );
        assert_eq!(map.host_cut(), Some([2, 9].as_slice()));
        assert_eq!(map.energy_edge_sets(), vec![[4, 6].as_slice()]);
        let joint = SamplingMapDefinition::parse(
            "block(lmb(7),at_cut(cut(9,2),intersect(surface(6,4),surface(6,3))))",
        )
        .unwrap();
        assert_eq!(
            SamplingMapDefinition::from_atom(joint.to_atom().as_view()).unwrap(),
            joint
        );
        assert_eq!(joint.host_cut(), Some([2, 9].as_slice()));
        assert_eq!(
            joint.energy_edge_sets(),
            vec![[4, 6].as_slice(), [3, 6].as_slice()]
        );
        for invalid in [
            "block(surface(1),surface(2))",
            "at_cut(lmb(1),surface(2))",
            "block(lmb(1,1),surface(2))",
            "at_cut(cut(1))",
        ] {
            assert!(SamplingMapDefinition::parse(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn parses_nested_definition_and_round_trips() {
        let source =
            "then(phase_space(cut(10, 2)), product(left(surface(12, 4)), right(complement(7, 8))))";
        let parsed = SamplingMapDefinition::parse(source).expect("valid map");
        assert_eq!(
            parsed,
            SamplingMapDefinition::Then(vec![
                SamplingMapDefinition::PhaseSpace(Box::new(SamplingMapDefinition::Cut(vec![
                    2, 10
                ]))),
                SamplingMapDefinition::Product(vec![
                    SamplingMapDefinition::Left(Box::new(SamplingMapDefinition::Surface(vec![
                        4, 12
                    ]))),
                    SamplingMapDefinition::Right(Box::new(SamplingMapDefinition::Complement(
                        vec![7, 8]
                    ))),
                ]),
            ])
        );
        let reparsed =
            SamplingMapDefinition::from_atom(parsed.to_atom().as_view()).expect("round trip");
        assert_eq!(parsed, reparsed);
    }

    #[test]
    fn rejects_unknown_root_and_invalid_nested_map() {
        for source in [
            "foo(surface(1))",
            "product(surface(1), foo(2))",
            "surface(1) + soft(2)",
            "user_defined::surface(1)",
        ] {
            assert!(
                SamplingMapDefinition::parse(source).is_err(),
                "accepted {source}"
            );
        }
    }

    #[test]
    fn rejects_bad_arities_and_duplicate_edges() {
        for source in [
            "surface()",
            "soft(1, 2)",
            "collinear(3, 3)",
            "then()",
            "phase_space(cut(1), cut(2))",
            "surface(2, 2)",
            "surface(-1)",
        ] {
            assert!(
                SamplingMapDefinition::parse(source).is_err(),
                "accepted {source}"
            );
        }
        assert!(SamplingMapDefinition::parse("lmb(1, 3, 1)").is_err());
    }

    #[test]
    fn preserves_parent_lmb_order_but_canonicalizes_edge_sets() {
        assert_eq!(
            SamplingMapDefinition::parse("lmb(4, 2, 7)").expect("valid lmb"),
            SamplingMapDefinition::Lmb(vec![4, 2, 7])
        );
        assert_eq!(
            SamplingMapDefinition::parse("surface(7, 2, 4)").expect("valid surface"),
            SamplingMapDefinition::Surface(vec![2, 4, 7])
        );
    }

    #[test]
    fn ordinary_kernel_round_trips_spherical_coordinates() {
        let settings = ParameterizationSettings {
            mode: ParameterizationMode::Spherical,
            ..Default::default()
        };
        let kernel =
            SamplingMapKernel::new(SamplingMapDefinition::Lmb(vec![4, 2]), settings, 173.0, 2)
                .unwrap();
        let point = kernel
            .forward(&[F(0.17), F(0.23), F(0.31), F(0.43), F(0.59), F(0.71)])
            .unwrap();
        assert!(point.jacobian.0.is_finite() && point.jacobian > F(0.0));
        assert!(point.inverse_jacobian.0.is_finite() && point.inverse_jacobian > F(0.0));
        assert!((point.jacobian * point.inverse_jacobian - F(1.0)).abs() < F(1.0e-12));
        assert!(point.residual < F(1.0e-12));
        let inverse = kernel.inverse(&point.loop_momenta).unwrap();
        assert!(inverse.residual < F(1.0e-12));
        assert_eq!(inverse.coordinates.len(), kernel.dimensions());
    }

    #[test]
    fn ordinary_kernel_supports_cartesian_shell_and_rejects_non_bijective_map() {
        let coordinates = [
            0.05, 0.2, 0.5, 0.4, 0.7, 0.9, 0.13, 0.87, 0.6, 0.31, 0.69, 0.95,
        ]
        .map(F);
        for mapping in [
            ParameterizationMapping::Log,
            ParameterizationMapping::Linear,
        ] {
            let kernel = SamplingMapKernel::new(
                SamplingMapDefinition::Complement(vec![0, 1, 2, 3]),
                ParameterizationSettings {
                    mode: ParameterizationMode::Cartesian,
                    mapping,
                    ..Default::default()
                },
                42.2,
                4,
            )
            .unwrap();
            let point = kernel.forward(&coordinates).unwrap();
            assert!(point.jacobian.0.is_finite() && point.jacobian > F(0.0));
            assert!(point.residual < F(1.0e-12));
            let inverse = kernel.inverse(&point.loop_momenta).unwrap();
            for (expected, actual) in coordinates.iter().zip(&inverse.coordinates) {
                assert!((expected - actual).abs() < F(1.0e-12));
            }
            assert!((point.jacobian * inverse.inverse_jacobian - F(1.0)).abs() < F(1.0e-12));

            // Independent derivatives of the forward map certify the inverse density.
            let step = 1.0e-6;
            let mut numerical_jacobian = 1.0;
            for dimension in 0..coordinates.len() {
                let mut lower = coordinates;
                let mut upper = coordinates;
                lower[dimension].0 -= step;
                upper[dimension].0 += step;
                let lower = kernel.forward(&lower).unwrap().loop_momenta;
                let upper = kernel.forward(&upper).unwrap().loop_momenta;
                let component = dimension % 3;
                let lower = &lower.0[dimension / 3];
                let upper = &upper.0[dimension / 3];
                let lower = [&lower.px, &lower.py, &lower.pz][component].0;
                let upper = [&upper.px, &upper.py, &upper.pz][component].0;
                numerical_jacobian *= (upper - lower) / (2.0 * step);
            }
            assert!((numerical_jacobian * inverse.inverse_jacobian.0 - 1.0).abs() < 1.0e-8);
        }
        assert!(
            SamplingMapKernel::new(
                SamplingMapDefinition::Lmb(vec![0]),
                ParameterizationSettings {
                    mode: ParameterizationMode::Cartesian,
                    mapping: ParameterizationMapping::Power,
                    ..Default::default()
                },
                42.2,
                1,
            )
            .is_err()
        );
        assert!(
            SamplingMapKernel::new(
                SamplingMapDefinition::Lmb(vec![0]),
                ParameterizationSettings {
                    mode: ParameterizationMode::HyperSphericalFlat,
                    ..Default::default()
                },
                42.2,
                1,
            )
            .is_err()
        );
    }

    #[test]
    fn ordinary_kernel_reports_dimension_and_support_errors() {
        let kernel = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![0]),
            ParameterizationSettings::default(),
            42.2,
            1,
        )
        .unwrap();
        assert!(kernel.forward(&[F(0.2), F(0.4)]).is_err());
        assert!(kernel.forward(&[F(0.0), F(0.4), F(0.7)]).is_err());
        assert!(
            SamplingMapKernel::new(
                SamplingMapDefinition::Surface(vec![0]),
                ParameterizationSettings::default(),
                42.2,
                1,
            )
            .is_err()
        );
    }

    #[test]
    fn sampling_maps_reject_unrepresentable_tails_without_clipping() {
        let kernel = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![0]),
            ParameterizationSettings {
                mode: ParameterizationMode::Cartesian,
                mapping: ParameterizationMapping::Log,
                ..Default::default()
            },
            42.2,
            1,
        )
        .unwrap();
        assert!(
            kernel
                .forward(&[F(f64::from_bits(1)), F(0.4), F(0.7)])
                .is_err()
        );
        assert!(
            kernel
                .inverse(&LoopMomenta::from_iter([ThreeMomentum::new(
                    F(4220.0),
                    F(10.0),
                    F(20.0),
                )]))
                .is_err()
        );
        let surface = SurfaceRadialMap::new(3, vec![0.0; 3], Some(1.0), 1.0, 1.0).unwrap();
        assert!(
            surface
                .forward(&[F(f64::from_bits(1)), F(0.3), F(0.4)])
                .is_err()
        );
        let implicit = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            1.0e-104,
            1.0,
            Arc::new(|_, radius| Ok((radius + 1.0, 1.0))),
        )
        .unwrap();
        assert!(implicit.forward(&[0.5, 0.3, 0.4]).is_err());
    }

    #[test]
    fn implicit_surface_inverse_preserves_resolvable_tail_and_angular_coordinates() {
        let map = ImplicitSurfaceRadialMap::new(
            2,
            vec![0.0; 2],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius + 1.0, 1.0))),
        )
        .unwrap();
        for coordinates in [[1.0e-16, 0.23], [0.3, 1.0e-16]] {
            let forward = map.forward(&coordinates).unwrap();
            let inverse = map.inverse(&forward.point).unwrap();
            for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
                assert!((actual / expected - 1.0).abs() < 1.0e-12);
            }
            assert!((forward.jacobian * inverse.inverse_jacobian - 1.0).abs() < 1.0e-12);
        }
        // The azimuth seam is a measure-zero chart boundary, not a nearby
        // interior point onto which the inverse may silently move its input.
        assert!(map.inverse(&[1.0, 0.0]).is_err());
        let polar_map = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius + 1.0, 1.0))),
        )
        .unwrap();
        assert!(polar_map.forward(&[0.3, 0.4, 1.0e-20]).is_err());
    }

    #[test]
    fn inverse_rounded_endpoints_are_retryable_but_original_endpoints_are_not() {
        let kernel = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![0]),
            ParameterizationSettings {
                mode: ParameterizationMode::Spherical,
                ..Default::default()
            },
            1.0,
            1,
        )
        .unwrap();
        let surface = SurfaceRadialMap::absent(3, vec![0.0; 3], 1.0, 1.0).unwrap();
        let implicit = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius + 1.0, 1.0))),
        )
        .unwrap();
        for map in [&kernel as &dyn SamplingMapComponent, &surface, &implicit] {
            let original = map
                .forward(&[1.0, 0.3, 0.6], &mut SamplingMapContext::detached(&[]))
                .unwrap_err();
            assert!(original.downcast_ref::<SamplingEvaluationError>().is_none());
            let recovered = map
                .inverse(
                    &[1.0e20, 2.0e20, 3.0e20],
                    &mut SamplingMapContext::detached(&[]),
                )
                .unwrap_err();
            assert!(
                matches!(
                    recovered.downcast_ref::<SamplingEvaluationError>(),
                    Some(SamplingEvaluationError::Unrepresentable { .. })
                ),
                "{recovered}"
            );
        }
        let one = F::<crate::utils::QuadFloat>::default().one();
        let surface = SurfaceRadialMap::absent(3, vec![one.zero().0; 3], 1.0, 1.0).unwrap();
        let point = [1, 2, 3].map(|value| one.from_i64(value) * one.from_i64(10).powi(20));
        let inverse = surface.inverse(&point).unwrap();
        assert!(
            inverse
                .coordinates
                .iter()
                .all(|coordinate| coordinate > &one.zero() && coordinate < &one)
        );
    }

    #[test]
    fn ordinary_kernel_preserves_finite_tails_at_native_rescue_precision() {
        let kernel = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![0]),
            ParameterizationSettings {
                mode: ParameterizationMode::Cartesian,
                mapping: ParameterizationMapping::Linear,
                ..Default::default()
            },
            42.2,
            1,
        )
        .unwrap();
        let one = F::<crate::utils::ArbPrec>::default().one();
        let tiny = &one / one.from_usize(10).powi(400);
        let coordinates = [tiny, &one / one.from_usize(3), &one / one.from_usize(2)];
        let mapped = kernel.forward(&coordinates).unwrap();
        assert!(mapped.jacobian.into_f64().is_infinite());
        assert_eq!(mapped.inverse_jacobian.into_f64(), 0.0);
        let inverse = kernel.inverse(&mapped.loop_momenta).unwrap();
        assert!(
            (inverse.coordinates[0].clone() / &coordinates[0] - &one).abs()
                < &one / one.from_usize(10).powi(12)
        );
    }

    #[test]
    fn implicit_surface_solver_certifies_residual_after_iteration_or_bracket_stop() {
        let callbacks: [ImplicitSurfaceRadialEvaluator; 2] = [
            Arc::new(|_, radius| Ok((if radius < 0.7 { -1.0 } else { 1.0 }, 1.0))),
            Arc::new(|_, radius| Ok((radius - 0.9, 1.0e6))),
        ];
        for evaluator in callbacks {
            let map = ImplicitSurfaceRadialMap::new(2, vec![0.0; 2], 1.0, 1.0, evaluator).unwrap();
            let error = map.forward(&[0.3, 0.4]).unwrap_err();
            assert!(error.to_string().contains("residual"), "{error}");
        }
        for offset in [-1.0e-14, 0.0, 1.0e-14] {
            let shallow = ImplicitSurfaceRadialMap::new(
                2,
                vec![0.0; 2],
                1.0,
                1.0,
                Arc::new(move |_, radius| Ok((radius + offset, 1.0))),
            )
            .unwrap();
            let error = shallow.forward(&[0.3, 0.4]).unwrap_err();
            assert!(matches!(
                error.downcast_ref::<SamplingEvaluationError>(),
                Some(SamplingEvaluationError::UncertainGeometry { .. })
            ));
        }
        let absent = ImplicitSurfaceRadialMap::new(
            2,
            vec![0.0; 2],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius + 1.0, 1.0))),
        )
        .unwrap();
        assert!(
            absent
                .forward(&[0.3, 0.4])
                .unwrap()
                .diagnostics
                .iter()
                .any(|value| value.contains("absent_fallback"))
        );
    }

    #[test]
    fn implicit_surface_map_tracks_directional_root_and_exact_contract() {
        let evaluator: ImplicitSurfaceRadialEvaluator = Arc::new(|direction, radius| {
            let root = 2.0 + 0.35 * direction[0] - 0.2 * direction[1];
            Ok((radius - root, 1.0))
        });
        let map =
            ImplicitSurfaceRadialMap::new(3, vec![0.4, -0.3, 0.2], 1.5, 2.0, evaluator).unwrap();
        assert_eq!(map.contract().support, SamplingSupport::Full);
        assert_eq!(map.contract().jacobian, SamplingJacobian::ExactImplicit);

        let coordinates = [0.37, 0.23, 0.61];
        let forward = map.forward(&coordinates).unwrap();
        assert!(forward.jacobian.is_finite() && forward.jacobian > 0.0);
        assert!(
            forward
                .diagnostics
                .iter()
                .any(|diagnostic| diagnostic == "implicit_surface:regular_root")
        );
        let inverse = map.inverse(&forward.point).unwrap();
        assert!(inverse.residual < 1.0e-10, "{}", inverse.residual);
        for (actual, expected) in inverse.coordinates.iter().zip(coordinates) {
            assert!((actual - expected).abs() < 1.0e-10);
        }
    }

    #[test]
    fn implicit_surface_map_has_normalized_absent_fibre_fallback() {
        let evaluator: ImplicitSurfaceRadialEvaluator =
            Arc::new(|_, radius| Ok((radius + 1.0, 1.0)));
        let map = ImplicitSurfaceRadialMap::new(2, vec![0.0, 0.0], 1.0, 1.0, evaluator).unwrap();
        let forward = map.forward(&[0.29, 0.71]).unwrap();
        assert!(
            forward
                .diagnostics
                .iter()
                .any(|diagnostic| diagnostic == "implicit_surface:absent_fallback")
        );
        let inverse = map.inverse(&forward.point).unwrap();
        assert!(inverse.residual < 1.0e-10, "{}", inverse.residual);
    }

    #[test]
    fn implicit_surface_map_reports_unbracketed_root_instead_of_absent_fibre() {
        let evaluator: ImplicitSurfaceRadialEvaluator = Arc::new(|_, _| Ok((-1.0, 0.0)));
        let map = ImplicitSurfaceRadialMap::new(2, vec![0.0, 0.0], 1.0, 1.0, evaluator).unwrap();
        let error = map
            .forward(&[0.29, 0.71])
            .expect_err("negative surface value without a root must be an error");
        assert!(error.to_string().contains("bracket expansions"));
    }

    #[test]
    fn implicit_fiber_prepares_once_and_preserves_certified_fallbacks() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let preparations = Arc::new(AtomicUsize::new(0));
        let calls = preparations.clone();
        let map = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            2.0,
            1.0,
            Arc::new(|_, _| panic!("conditional evaluator required")),
        )
        .unwrap()
        .with_context_preparer(Arc::new(move |context| {
            calls.fetch_add(1, Ordering::Relaxed);
            let status = match context[0] {
                0.0 => PreparedSurfaceStatus::pinched(0.0)?,
                1.0 => PreparedSurfaceStatus::absent("certified empty fiber")?,
                _ => PreparedSurfaceStatus::existing(Some(1.0))?,
            };
            Ok((vec![context[0], 0.0, 0.0], status))
        }))
        .with_context_evaluator(Arc::new(|_, radius, center, context| {
            assert_eq!(center[0], context[0]);
            assert!(
                context[0] > 1.0,
                "certified fallback must not evaluate a root"
            );
            Ok((radius - 1.0, 1.0))
        }));
        for context in [0.0, 1.0, 2.0] {
            let count = preparations.load(Ordering::Relaxed);
            let forward = map
                .forward_with_context(&[0.2, 0.3, 0.7], &[context])
                .unwrap();
            assert_eq!(preparations.load(Ordering::Relaxed), count + 1);
            let inverse = map
                .inverse_with_context(&forward.point, &[context])
                .unwrap();
            assert_eq!(preparations.load(Ordering::Relaxed), count + 2);
            assert_eq!(
                map.inverse_density(
                    &forward.point,
                    &mut SamplingMapContext::detached(&[context]),
                )
                .unwrap(),
                Some(inverse.inverse_jacobian),
            );
            assert_eq!(preparations.load(Ordering::Relaxed), count + 3);
            assert!((forward.jacobian * inverse.inverse_jacobian - 1.0).abs() < 1.0e-12);
            if context < 2.0 {
                let ordinary = SurfaceRadialMap::new(3, vec![context, 0.0, 0.0], None, 2.0, 1.0)
                    .unwrap()
                    .forward(&[F(0.2), F(0.3), F(0.7)])
                    .unwrap();
                assert!((forward.jacobian - ordinary.jacobian.0).abs() < 1.0e-12);
            }
        }
        let before = preparations.load(Ordering::Relaxed);
        let frozen = map.freeze_context(&[2.0]).unwrap();
        assert_eq!(frozen.contract().support, SamplingSupport::Full);
        let forward = frozen.forward(&[0.2, 0.3, 0.7]).unwrap();
        frozen.inverse(&forward.point).unwrap();
        assert_eq!(preparations.load(Ordering::Relaxed), before + 1);
    }

    #[test]
    fn implicit_surface_map_accepts_conditional_complement_context() {
        let map = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            2.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
        )
        .unwrap()
        .with_context_evaluator(Arc::new(|_, radius, _, context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((radius - (1.0 + shift), 1.0))
        }));
        assert_eq!(map.contract().support, SamplingSupport::Full);
        assert!(map.contract().requires_context);
        let coordinates = [0.2, 0.27, 0.61];
        let mapped = map.forward_with_context(&coordinates, &[0.5]).unwrap();
        let radius = mapped
            .point
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        assert!((radius - 0.7).abs() < 1.0e-9);
        let inverse = map.inverse_with_context(&mapped.point, &[0.5]).unwrap();
        assert!(inverse.residual < 1.0e-9);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9)
        );
    }

    #[test]
    fn implicit_surface_map_resolves_context_dependent_center_and_root() {
        let map = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            2.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
        )
        .unwrap()
        .with_context_preparer(Arc::new(|context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((
                vec![shift, 0.0, 0.0],
                PreparedSurfaceStatus::existing(None)?,
            ))
        }))
        .with_context_evaluator(Arc::new(|_, radius, _, context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((radius - (1.0 + shift), 1.0))
        }));

        assert_eq!(map.contract().support, SamplingSupport::Full);
        assert!(map.contract().requires_context);
        let context = [0.5];
        let coordinates = [0.23, 0.31, 0.67];
        let forward = map
            .forward_with_context(&coordinates, &context)
            .expect("context-dependent directional chart");
        let displacement = forward
            .point
            .iter()
            .zip([0.5, 0.0, 0.0])
            .map(|(point, center)| (point - center).powi(2))
            .sum::<f64>()
            .sqrt();
        let split = 1.5 / (1.5 + 2.0);
        let expected_radius = 1.5 * (coordinates[0] / split);
        assert!((displacement - expected_radius).abs() < 1.0e-9);

        let inverse = map
            .inverse_with_context(&forward.point, &context)
            .expect("inverse context-dependent directional chart");
        assert!(inverse.residual < 1.0e-9, "{}", inverse.residual);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9)
        );
    }

    #[test]
    fn implicit_surface_map_rejects_invalid_context_center() {
        let wrong_dimension = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
        )
        .unwrap()
        .with_context_preparer(Arc::new(|_| {
            Ok((vec![0.0, 0.0], PreparedSurfaceStatus::existing(None)?))
        }));
        let error = wrong_dimension
            .forward_with_context(&[0.2, 0.3, 0.7], &[])
            .expect_err("wrong-dimensional context centre must fail");
        assert!(
            error
                .to_string()
                .contains("centre has dimension 2, expected 3")
        );

        let nonfinite = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
        )
        .unwrap()
        .with_context_preparer(Arc::new(|_| {
            Ok((
                vec![f64::NAN, 0.0, 0.0],
                PreparedSurfaceStatus::existing(None)?,
            ))
        }));
        let error = nonfinite
            .forward_with_context(&[0.2, 0.3, 0.7], &[])
            .expect_err("non-finite context centre must fail");
        assert!(matches!(
            error.downcast_ref::<SamplingEvaluationError>(),
            Some(SamplingEvaluationError::Unrepresentable {
                operation: "implicit surface context centre",
                ..
            })
        ));
    }

    #[test]
    fn context_center_surface_is_full_support_after_full_then_block() {
        let surface = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            2.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
        )
        .unwrap()
        .with_context_preparer(Arc::new(|context| {
            Ok((
                vec![context.first().copied().unwrap_or(0.0), 0.0, 0.0],
                PreparedSurfaceStatus::existing(None)?,
            ))
        }))
        .with_context_evaluator(Arc::new(|_, radius, _, context| {
            Ok((
                radius - (1.0 + context.first().copied().unwrap_or(0.0)),
                1.0,
            ))
        }));
        let preceding = SamplingMapAffine::new(vec![vec![1.0]], vec![0.0]).unwrap();
        let composition =
            SamplingMapComposition::then(vec![Box::new(preceding), Box::new(surface)])
                .expect("ordered conditional composition");
        assert_eq!(composition.contract().support, SamplingSupport::Full);
        assert!(!composition.contract().requires_context);
        let coordinates = [0.5, 0.23, 0.31, 0.67];
        let forward = composition
            .forward(&coordinates, &mut SamplingMapContext::detached(&[]))
            .expect("forward ordered composition");
        let inverse = composition
            .inverse(&forward.point, &mut SamplingMapContext::detached(&[]))
            .expect("inverse ordered composition")
            .expect("full-support inverse");
        assert!(inverse.residual < 1.0e-9, "{}", inverse.residual);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-9)
        );
    }

    #[test]
    fn surface_radial_map_round_trips_on_both_sides_of_regular_surface() {
        let map = SurfaceRadialMap::new(3, vec![1.0, -2.0, 0.5], Some(4.0), 2.0, 1.0)
            .expect("valid radial map");
        assert_eq!(map.contract().support, SamplingSupport::Full);
        for coordinates in [
            vec![F(0.11), F(0.23), F(0.37)],
            vec![F(0.63), F(0.41), F(0.79)],
            vec![F(0.93), F(0.71), F(0.19)],
        ] {
            let forward = map.forward(&coordinates).expect("forward map");
            assert!(forward.jacobian.0.is_finite() && forward.jacobian > F(0.0));
            let inverse = map.inverse(&forward.point).expect("inverse map");
            assert!(inverse.residual < F(1.0e-11));
            assert!(
                inverse
                    .coordinates
                    .iter()
                    .zip(&coordinates)
                    .all(|(actual, expected)| (actual - expected).abs() < F(1.0e-11))
            );
            assert!((forward.jacobian * forward.inverse_jacobian - F(1.0)).abs() < F(1.0e-11));
        }
    }

    #[test]
    fn absent_surface_fallback_has_full_radial_support() {
        let map = SurfaceRadialMap::absent(4, vec![0.0; 4], 3.0, 1.7).expect("valid fallback");
        assert_eq!(map.threshold_radius(), None);
        for coordinates in [
            vec![F(0.05), F(0.17), F(0.29), F(0.41)],
            vec![F(0.95), F(0.61), F(0.73), F(0.83)],
        ] {
            let forward = map.forward(&coordinates).expect("forward map");
            assert!(forward.radius > F(0.0));
            let inverse = map.inverse(&forward.point).expect("inverse map");
            assert!(inverse.residual < F(1.0e-10));
        }
    }

    #[test]
    fn surface_radial_map_uses_correct_sphere_measure_in_two_and_three_dimensions() {
        let map_2 = SurfaceRadialMap::absent(2, vec![0.0; 2], 2.0, 1.0).unwrap();
        let point_2 = map_2.forward(&[F(0.4), F(0.25)]).unwrap();
        // r = beta*u/(1-u), dr/du = beta/(1-u)^2 for the absent fallback.
        let radius_2 = 2.0 * 0.4 / 0.6;
        let radial_jacobian_2 = 2.0 / 0.6_f64.powi(2);
        let expected_2 = radial_jacobian_2 * 2.0 * std::f64::consts::PI * radius_2;
        assert!((point_2.jacobian.0 - expected_2).abs() < 1.0e-12);

        let map_3 = SurfaceRadialMap::absent(3, vec![0.0; 3], 2.0, 1.0).unwrap();
        let point_3 = map_3.forward(&[F(0.4), F(0.25), F(0.7)]).unwrap();
        let expected_3 = radial_jacobian_2 * 4.0 * std::f64::consts::PI * radius_2.powi(2);
        assert!((point_3.jacobian.0 - expected_3).abs() < 1.0e-12);
    }

    /// Integrate a normalized Gaussian after a complete push-forward map.
    ///
    /// This is the small, graph-independent acceptance test used when a map
    /// replaces the legacy unit-volume fixture: the cube average must include
    /// the map's exact determinant and converge to one without any process
    /// integrand machinery.
    #[test]
    fn radial_sampling_map_pushes_forward_a_normalized_gaussian() {
        let map = SurfaceRadialMap::absent(3, vec![0.0; 3], 2.0, 1.0).unwrap();
        let bins = 32usize;
        let normalisation = (2.0 * std::f64::consts::PI).powf(-1.5);
        let mut integral = 0.0;
        for i in 0..bins {
            for j in 0..bins {
                for k in 0..bins {
                    let coordinates = vec![
                        F((i as f64 + 0.5) / bins as f64),
                        F((j as f64 + 0.5) / bins as f64),
                        F((k as f64 + 0.5) / bins as f64),
                    ];
                    let point = map.forward(&coordinates).unwrap();
                    let radius_squared = point
                        .point
                        .iter()
                        .map(|component| component.0.powi(2))
                        .sum::<f64>();
                    integral += normalisation * (-0.5 * radius_squared).exp() * point.jacobian.0;
                }
            }
        }
        integral /= bins.pow(3) as f64;
        assert!(
            (integral - 1.0).abs() < 2.0e-3,
            "Gaussian integral = {integral}"
        );
    }

    #[test]
    fn ordinary_sampling_map_pushes_forward_a_normalized_gaussian() {
        let settings = ParameterizationSettings {
            mode: ParameterizationMode::Spherical,
            ..Default::default()
        };
        let map =
            SamplingMapKernel::new(SamplingMapDefinition::Lmb(vec![0]), settings, 2.0, 1).unwrap();
        let bins = 24usize;
        let normalisation = (2.0 * std::f64::consts::PI).powf(-1.5);
        let mut integral = 0.0;
        for i in 0..bins {
            for j in 0..bins {
                for k in 0..bins {
                    let coordinates = vec![
                        F((i as f64 + 0.5) / bins as f64),
                        F((j as f64 + 0.5) / bins as f64),
                        F((k as f64 + 0.5) / bins as f64),
                    ];
                    let point = map.forward(&coordinates).unwrap();
                    let radius_squared = point
                        .loop_momenta
                        .0
                        .iter()
                        .flat_map(|momentum| [momentum.px.0, momentum.py.0, momentum.pz.0])
                        .map(|component| component.powi(2))
                        .sum::<f64>();
                    integral += normalisation * (-0.5 * radius_squared).exp() * point.jacobian.0;
                }
            }
        }
        integral /= bins.pow(3) as f64;
        assert!(
            (integral - 1.0).abs() < 2.0e-2,
            "Gaussian integral = {integral}"
        );
    }

    #[test]
    fn surface_radial_map_rejects_invalid_geometry() {
        assert!(SurfaceRadialMap::new(1, vec![0.0], Some(1.0), 1.0, 1.0).is_err());
        assert!(SurfaceRadialMap::new(3, vec![0.0; 2], Some(1.0), 1.0, 1.0).is_err());
        assert!(SurfaceRadialMap::new(3, vec![0.0; 3], Some(-1.0), 1.0, 1.0).is_err());
        assert!(SurfaceRadialMap::new(3, vec![0.0; 3], Some(1.0), 0.0, 1.0).is_err());
    }

    #[test]
    fn product_composition_multiplies_block_jacobians_and_round_trips() {
        let first = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![0]),
            ParameterizationSettings::default(),
            42.2,
            1,
        )
        .unwrap();
        let second = SurfaceRadialMap::absent(2, vec![1.0, -2.0], 2.0, 1.0).unwrap();
        let first_jacobian = first.forward(&[F(0.2), F(0.4), F(0.7)]).unwrap().jacobian.0;
        let second_jacobian = second.forward(&[F(0.3), F(0.6)]).unwrap().jacobian.0;
        let composition = SamplingMapComposition::product(vec![Box::new(first), Box::new(second)])
            .expect("valid product");
        let coordinates = [0.2, 0.4, 0.7, 0.3, 0.6];
        let mapped = composition
            .forward(&coordinates, &mut SamplingMapContext::detached(&[]))
            .unwrap();
        assert_eq!(mapped.coordinates.len(), 5);
        assert_eq!(mapped.point.len(), 5);
        assert!((mapped.jacobian - first_jacobian * second_jacobian).abs() < 1.0e-10);
        assert!((mapped.jacobian * mapped.inverse_jacobian - 1.0).abs() < 1.0e-10);
        assert!(mapped.residual < 1.0e-11);
        let inverse = composition
            .inverse(&mapped.point, &mut SamplingMapContext::detached(&[]))
            .unwrap()
            .expect("full-support inverse");
        assert_eq!(
            composition
                .inverse_density(&mapped.point, &mut SamplingMapContext::detached(&[]))
                .unwrap(),
            Some(inverse.inverse_jacobian),
        );
        assert!(inverse.residual < 1.0e-10);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-10)
        );
    }

    #[test]
    fn embedded_product_permutates_master_frame_without_changing_jacobian() {
        let first = SurfaceRadialMap::absent(2, vec![0.0, 0.0], 1.3, 1.0).unwrap();
        let second = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![7]),
            ParameterizationSettings::default(),
            42.2,
            1,
        )
        .unwrap();
        let first_eval = <SurfaceRadialMap as SamplingMapComponent>::forward(
            &first,
            &[0.23, 0.71],
            &mut SamplingMapContext::detached(&[]),
        )
        .unwrap();
        let second_eval = <SamplingMapKernel as SamplingMapComponent>::forward(
            &second,
            &[0.17, 0.43, 0.89],
            &mut SamplingMapContext::detached(&[]),
        )
        .unwrap();
        let embedding = SamplingMapEmbedding::product(
            vec![Box::new(first), Box::new(second)],
            vec![2, 3, 0, 1, 4],
        )
        .expect("valid master-frame permutation");
        let coordinates = [0.23, 0.71, 0.17, 0.43, 0.89];
        let mapped = embedding
            .forward(&coordinates, &mut SamplingMapContext::detached(&[]))
            .unwrap();
        assert_eq!(mapped.point[0..2], second_eval.point[0..2]);
        assert_eq!(mapped.point[2..4], first_eval.point[0..2]);
        assert_eq!(mapped.point[4], second_eval.point[2]);
        assert!((mapped.jacobian - first_eval.jacobian * second_eval.jacobian).abs() < 1.0e-10);
        let inverse = embedding
            .inverse(&mapped.point, &mut SamplingMapContext::detached(&[]))
            .unwrap()
            .expect("full-support inverse");
        assert_eq!(
            embedding
                .inverse_density(&mapped.point, &mut SamplingMapContext::detached(&[]))
                .unwrap(),
            Some(inverse.inverse_jacobian),
        );
        assert!(inverse.residual < 1.0e-10);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-10)
        );
    }

    #[test]
    fn embedded_worker_clone_owns_nested_eager_buffers() -> Result<()> {
        use std::{sync::mpsc, time::Duration};

        crate::initialisation::test_initialise()?;
        // Retain the parent's buffer handle to test ownership through two
        // embedded compositions. Leaf maps already clone their eager buffers;
        // a shared outer composition must not prevent that clone from running.
        struct EagerProbe(Arc<Mutex<SamplingExpressionEvaluator>>);

        impl std::fmt::Debug for EagerProbe {
            fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("EagerProbe")
            }
        }

        impl Clone for EagerProbe {
            fn clone(&self) -> Self {
                Self(Arc::new(Mutex::new(self.0.lock().unwrap().clone())))
            }
        }

        impl SamplingMapComponent for EagerProbe {
            fn dimensions(&self) -> usize {
                1
            }
            fn output_dimensions(&self) -> usize {
                1
            }
            fn contract(&self) -> SamplingMapContract {
                SamplingMapContract {
                    support: SamplingSupport::Full,
                    requires_context: false,
                    requires_proposal_policy: false,
                    jacobian: SamplingJacobian::ExactForward,
                }
            }
            fn name(&self) -> &'static str {
                "eager_clone_probe"
            }
            fn forward(
                &self,
                coordinates: &[f64],
                _context: &mut SamplingMapContext<'_, f64>,
            ) -> Result<SamplingMapEvaluation> {
                let evaluated = self
                    .0
                    .lock()
                    .unwrap()
                    .evaluate_with_real_jacobian(coordinates, None)?;
                Ok(SamplingMapEvaluation {
                    coordinates: coordinates.to_vec(),
                    point: evaluated.values,
                    jacobian: evaluated.determinant,
                    inverse_jacobian: evaluated.determinant.recip(),
                    residual: 0.0,
                    support: SamplingSupport::Full,
                    diagnostics: Vec::new(),
                })
            }
            fn inverse(
                &self,
                point: &[f64],
                context: &mut SamplingMapContext<'_, f64>,
            ) -> Result<Option<SamplingMapEvaluation>> {
                self.forward(&[(point[0] - 1.0) / 3.0], context).map(Some)
            }
        }

        let program = Arc::new(Mutex::new(SamplingExpressionEvaluator::new(
            [try_parse!("3*embedded_clone::x+1").map_err(|error| eyre!(error))?],
            [try_parse!("embedded_clone::x").map_err(|error| eyre!(error))?],
            &[0],
        )?));
        let inner =
            SamplingMapEmbedding::product(vec![Box::new(EagerProbe(program.clone()))], vec![0])?;
        let parent = SamplingMapEmbedding::from_composition(
            SamplingMapComposition::then(vec![Box::new(inner)])?,
            vec![0],
        )?;
        let worker = parent.clone();
        let parent_lock = program.lock().unwrap();
        let (sender, receiver) = mpsc::channel();
        let handle = std::thread::spawn(move || -> Result<()> {
            let mapped = worker.forward(&[0.25], &mut SamplingMapContext::detached(&[]))?;
            let inverse = worker
                .inverse(&mapped.point, &mut SamplingMapContext::detached(&[]))?
                .unwrap();
            sender.send((mapped, inverse))?;
            Ok(())
        });
        let result = receiver.recv_timeout(Duration::from_secs(5));
        // Always release the parent before joining, including a regression
        // timeout, so a shared mutex produces a bounded failure rather than a hang.
        drop(parent_lock);
        handle.join().unwrap()?;
        let (mapped, inverse) = result.expect("worker used its parent's eager buffer");
        assert_eq!(mapped.point, [1.75]);
        assert_eq!(mapped.jacobian, 3.0);
        assert_eq!(inverse.coordinates, [0.25]);
        assert_eq!(mapped.jacobian * inverse.inverse_jacobian, 1.0);
        Ok(())
    }

    #[derive(Clone, Debug)]
    struct ContextShiftMap;

    impl SamplingMapComponent for ContextShiftMap {
        fn dimensions(&self) -> usize {
            1
        }

        fn output_dimensions(&self) -> usize {
            1
        }

        fn contract(&self) -> SamplingMapContract {
            SamplingMapContract {
                support: SamplingSupport::Full,
                requires_context: true,
                requires_proposal_policy: false,
                jacobian: SamplingJacobian::ExactForward,
            }
        }

        fn name(&self) -> &'static str {
            "context_shift"
        }

        fn forward(
            &self,
            coordinates: &[f64],
            context: &mut SamplingMapContext<'_>,
        ) -> Result<SamplingMapEvaluation> {
            let shift = context.previous.first().copied().unwrap_or(0.0);
            let point = coordinates[0] + shift;
            Ok(SamplingMapEvaluation {
                coordinates: coordinates.to_vec(),
                point: vec![point],
                jacobian: 1.0,
                inverse_jacobian: 1.0,
                residual: 0.0,
                support: self.contract().support,
                diagnostics: vec![format!("shift={shift}")],
            })
        }

        fn inverse(
            &self,
            point: &[f64],
            context: &mut SamplingMapContext<'_>,
        ) -> Result<Option<SamplingMapEvaluation>> {
            let shift = context.previous.first().copied().unwrap_or(0.0);
            let coordinate = point[0] - shift;
            Ok(Some(SamplingMapEvaluation {
                coordinates: vec![coordinate],
                point: point.to_vec(),
                jacobian: 1.0,
                inverse_jacobian: 1.0,
                residual: 0.0,
                support: self.contract().support,
                diagnostics: vec![format!("shift={shift}")],
            }))
        }
    }

    #[test]
    fn then_composition_passes_previous_outputs_as_context() {
        let composition = SamplingMapComposition::then(vec![
            Box::new(ContextShiftMap),
            Box::new(ContextShiftMap),
        ])
        .expect("valid conditional composition");
        assert!(composition.is_then());
        assert_eq!(composition.contract().support, SamplingSupport::Full);
        assert!(composition.contract().requires_context);
        let mapped = composition
            .forward(&[0.2, 0.4], &mut SamplingMapContext::detached(&[]))
            .unwrap();
        assert!(
            mapped
                .point
                .iter()
                .zip([0.2, 0.6])
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-15)
        );
        assert_eq!(mapped.diagnostics, vec!["shift=0", "shift=0.2"]);
        let inverse = composition
            .inverse(&mapped.point, &mut SamplingMapContext::detached(&[]))
            .unwrap()
            .expect("full-support inverse");
        assert!(
            inverse
                .coordinates
                .iter()
                .zip([0.2, 0.4])
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-15)
        );
    }

    #[test]
    fn then_composition_accepts_outer_context_for_forward_and_inverse() {
        let composition = SamplingMapComposition::then(vec![
            Box::new(ContextShiftMap),
            Box::new(ContextShiftMap),
        ])
        .expect("valid conditional composition");
        let mapped = composition
            .forward(&[0.2, 0.4], &mut SamplingMapContext::detached(&[0.1]))
            .unwrap();
        assert!(
            mapped
                .point
                .iter()
                .zip([0.3, 0.5])
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-15)
        );
        let inverse = composition
            .inverse(&mapped.point, &mut SamplingMapContext::detached(&[0.1]))
            .unwrap()
            .expect("full-support inverse");
        assert!(
            inverse
                .coordinates
                .iter()
                .zip([0.2, 0.4])
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-15)
        );
    }

    #[test]
    fn then_composition_routes_complement_context_into_implicit_surface_map() {
        let implicit = ImplicitSurfaceRadialMap::new(
            3,
            vec![0.0; 3],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0, 1.0))),
        )
        .unwrap()
        .with_context_evaluator(Arc::new(|_, radius, _, context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((radius - (1.0 + shift), 1.0))
        }));
        let composition =
            SamplingMapComposition::then(vec![Box::new(ContextShiftMap), Box::new(implicit)])
                .expect("valid conditional surface composition");
        let coordinates = [0.2, 0.3, 0.27, 0.61];
        let mapped = composition
            .forward(&coordinates, &mut SamplingMapContext::detached(&[]))
            .unwrap();
        let surface_radius = mapped.point[1..]
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        assert!((surface_radius - 0.66).abs() < 1.0e-9);
        assert!(
            mapped
                .diagnostics
                .iter()
                .any(|diagnostic| diagnostic == "implicit_surface:regular_root")
        );
        assert_eq!(mapped.support, SamplingSupport::Full);
        assert!(composition.contract().requires_context);
        let inverse = composition
            .inverse(&mapped.point, &mut SamplingMapContext::detached(&[]))
            .unwrap()
            .expect("full-support inverse");
        assert!(inverse.residual < 1.0e-10);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-10)
        );
    }

    #[test]
    fn compositions_report_empty_and_dimension_errors() {
        assert!(SamplingMapComposition::<f64>::product(Vec::new()).is_err());
        let composition = SamplingMapComposition::product(vec![Box::new(ContextShiftMap)])
            .expect("valid product");
        let error = composition
            .forward(&[0.2, 0.3], &mut SamplingMapContext::detached(&[]))
            .unwrap_err()
            .to_string();
        assert!(error.contains("received dimension 2, expected 1"));
        let error = composition
            .inverse(&[0.2, 0.3], &mut SamplingMapContext::detached(&[]))
            .unwrap_err()
            .to_string();
        assert!(error.contains("received output dimension 2, expected 1"));
    }

    #[test]
    fn generic_acceptance_harness_checks_normalized_surface_map() {
        let map = SurfaceRadialMap::absent(3, vec![0.0; 3], 2.0, 1.0).unwrap();
        let report =
            SamplingMapAcceptanceReport::normalized_gaussian(&map, 4096, 1.5, &[0.0, 0.0, 0.0])
                .unwrap();
        assert_eq!(report.finite_sample_count, report.sample_count);
        assert!((report.normalization - 1.0).abs() < 2.0e-2);
        assert!(report.normalization_stderr.is_finite());
    }

    #[test]
    fn generic_acceptance_harness_checks_embedded_partial_surface_map() {
        let surface = SurfaceRadialMap::absent(2, vec![0.0; 2], 1.3, 1.0).unwrap();
        let complement = SamplingMapKernel::new(
            SamplingMapDefinition::Lmb(vec![7]),
            ParameterizationSettings::default(),
            42.2,
            1,
        )
        .unwrap();
        let map = SamplingMapEmbedding::product(
            vec![Box::new(surface), Box::new(complement)],
            vec![2, 3, 0, 1, 4],
        )
        .unwrap();
        let report = SamplingMapAcceptanceReport::normalized_gaussian(
            &map,
            4096,
            1.5,
            &[0.0, 0.0, 0.0, 0.0, 0.0],
        )
        .unwrap();
        assert_eq!(report.finite_sample_count, report.sample_count);
        assert!((report.normalization - 1.0).abs() < 3.0e-2);
    }
}
