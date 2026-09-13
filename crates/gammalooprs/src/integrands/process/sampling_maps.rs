//! Symbolica-validated descriptions of graph-aware sampling maps.
//!
//! The graph-independent map language and the reusable ordinary-map numerical
//! kernel live here.  Graph-aware surface and cut maps are compiled by the
//! process layer. Parsing is structural: a complete Symbolica expression is
//! checked recursively, so a valid nested call cannot hide an invalid root or
//! sibling.

use std::sync::Arc;

use crate::momentum::ThreeMomentum;
use crate::momentum::sample::LoopMomenta;
use crate::settings::runtime::{
    ParameterizationMapping, ParameterizationMode, ParameterizationSettings,
};
use crate::utils::{F, FloatLike, global_inv_parameterize, global_parameterize};
use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomView},
    symbol, try_parse,
};

/// The geometric role of a resolved sampling map.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
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
            _ => Err(eyre!(
                "unknown sampling-map constructor `{name}`; expected lmb, surface, cut, soft, collinear, complement, product, intersect, then, phase_space, left or right"
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
    /// The map is valid only after an earlier conditional map has prepared data.
    Conditional,
    /// The map has explicitly enumerated regular branches.
    Branched,
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

/// Static contract carried by a future compiled sampling-map kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SamplingMapContract {
    pub support: SamplingSupport,
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

/// Result of evaluating one graph-independent sampling component.
///
/// The point is represented as a flat vector in the master raw coordinate
/// frame.  A composition therefore never has to reinterpret a child's local
/// coordinates.  `diagnostics` carries branch/support information to callers;
/// it is deliberately data rather than logging so acceptance tests can audit
/// every branch taken by a map.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingMapEvaluation {
    pub coordinates: Vec<f64>,
    pub point: Vec<f64>,
    pub jacobian: f64,
    pub inverse_jacobian: f64,
    pub residual: f64,
    pub support: SamplingSupport,
    pub diagnostics: Vec<String>,
}

/// A graph-independent map component which can participate in a product or
/// ordered conditional composition.  `context` contains previously produced
/// raw coordinates for `then` maps and is empty for independent products.
pub trait SamplingMapComponent: std::fmt::Debug + Send + Sync {
    fn dimensions(&self) -> usize;
    fn output_dimensions(&self) -> usize;
    fn contract(&self) -> SamplingMapContract;
    fn name(&self) -> &'static str;
    fn forward(&self, coordinates: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation>;
    fn inverse(&self, point: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation>;
}

/// An exact affine change of coordinates between two equal-dimensional frames.
///
/// The matrix is stored row-major and the map is `point = matrix * coordinates
/// + translation`.  This is the graph-independent numerical owner for routing
/// a selected LMB into a parent frame; graph code supplies the matrix and
/// translation after resolving the relevant edge signatures and external data.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingMapAffine {
    matrix: Vec<f64>,
    inverse_matrix: Vec<f64>,
    translation: Vec<f64>,
    dimension: usize,
    determinant: f64,
}

impl SamplingMapAffine {
    /// Construct an affine map from a square matrix and translation vector.
    ///
    /// Partial pivoting is used while constructing the inverse.  Singular and
    /// numerically rank-deficient matrices are rejected before a map can enter
    /// a sampling composition, since they have no valid push-forward density.
    pub fn new(matrix: Vec<Vec<f64>>, translation: Vec<f64>) -> Result<Self> {
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

        let mut augmented = vec![vec![0.0; 2 * dimension]; dimension];
        let mut scale = 0.0_f64;
        for (row_index, row) in matrix.iter().enumerate() {
            for (column_index, value) in row.iter().copied().enumerate() {
                augmented[row_index][column_index] = value;
                scale = scale.max(value.abs());
            }
            augmented[row_index][dimension + row_index] = 1.0;
        }
        if scale == 0.0 {
            return Err(eyre!("affine sampling map matrix is singular"));
        }
        let pivot_tolerance = scale * f64::EPSILON * dimension as f64 * 32.0;
        let mut determinant = 1.0;
        let mut row_sign = 1.0;
        for pivot_column in 0..dimension {
            let (pivot_row, pivot_abs) = (pivot_column..dimension)
                .map(|row| (row, augmented[row][pivot_column].abs()))
                .max_by(|(_, left), (_, right)| left.total_cmp(right))
                .expect("non-empty pivot range");
            if pivot_abs <= pivot_tolerance || !pivot_abs.is_finite() {
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
            let pivot = augmented[pivot_column][pivot_column];
            determinant *= pivot;
            let inverse_pivot = 1.0 / pivot;
            for value in &mut augmented[pivot_column] {
                *value *= inverse_pivot;
            }
            for row in 0..dimension {
                if row == pivot_column {
                    continue;
                }
                let factor = augmented[row][pivot_column];
                if factor == 0.0 {
                    continue;
                }
                for column in 0..2 * dimension {
                    augmented[row][column] -= factor * augmented[pivot_column][column];
                }
            }
        }
        determinant *= row_sign;
        if !determinant.is_finite() || determinant == 0.0 {
            return Err(eyre!(
                "affine sampling map determinant is not finite and non-zero: {determinant}"
            ));
        }
        let inverse_matrix = augmented
            .into_iter()
            .flat_map(|row| row.into_iter().skip(dimension))
            .collect::<Vec<_>>();
        let determinant = determinant.abs();
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(eyre!(
                "affine sampling map Jacobian is not finite and positive: {determinant}"
            ));
        }
        Ok(Self {
            matrix: matrix.into_iter().flatten().collect(),
            inverse_matrix,
            translation,
            dimension,
            determinant,
        })
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn determinant(&self) -> f64 {
        self.determinant
    }

    pub fn matrix(&self) -> &[f64] {
        &self.matrix
    }

    pub fn inverse_matrix(&self) -> &[f64] {
        &self.inverse_matrix
    }

    pub fn translation(&self) -> &[f64] {
        &self.translation
    }

    fn apply(&self, matrix: &[f64], input: &[f64], translation: &[f64]) -> Vec<f64> {
        (0..self.dimension)
            .map(|row| {
                translation[row]
                    + matrix[row * self.dimension..(row + 1) * self.dimension]
                        .iter()
                        .zip(input)
                        .map(|(coefficient, value)| coefficient * value)
                        .sum::<f64>()
            })
            .collect()
    }

    fn validate_dimension(&self, values: &[f64], role: &str) -> Result<()> {
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

impl SamplingMapComponent for SamplingMapAffine {
    fn dimensions(&self) -> usize {
        self.dimension
    }

    fn output_dimensions(&self) -> usize {
        self.dimension
    }

    fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: SamplingSupport::Full,
            jacobian: SamplingJacobian::ExactForward,
        }
    }

    fn name(&self) -> &'static str {
        "affine"
    }

    fn forward(&self, coordinates: &[f64], _context: &[f64]) -> Result<SamplingMapEvaluation> {
        self.validate_dimension(coordinates, "input")?;
        let point = self.apply(&self.matrix, coordinates, &self.translation);
        let translated = point
            .iter()
            .zip(&self.translation)
            .map(|(value, shift)| value - shift)
            .collect::<Vec<_>>();
        let recovered = self.apply(
            &self.inverse_matrix,
            &translated,
            &vec![0.0; self.dimension],
        );
        Ok(SamplingMapEvaluation {
            coordinates: coordinates.to_vec(),
            point,
            jacobian: self.determinant,
            inverse_jacobian: 1.0 / self.determinant,
            residual: max_coordinate_residual_f64(coordinates, &recovered),
            support: SamplingSupport::Full,
            diagnostics: vec![format!("affine dimension={}", self.dimension)],
        })
    }

    fn inverse(&self, point: &[f64], _context: &[f64]) -> Result<SamplingMapEvaluation> {
        self.validate_dimension(point, "point")?;
        let translated = point
            .iter()
            .zip(&self.translation)
            .map(|(value, shift)| value - shift)
            .collect::<Vec<_>>();
        let coordinates = self.apply(
            &self.inverse_matrix,
            &translated,
            &vec![0.0; self.dimension],
        );
        let recovered = self.apply(&self.matrix, &coordinates, &self.translation);
        Ok(SamplingMapEvaluation {
            coordinates,
            point: point.to_vec(),
            jacobian: self.determinant,
            inverse_jacobian: 1.0 / self.determinant,
            residual: max_coordinate_residual_f64(&recovered, point),
            support: SamplingSupport::Full,
            diagnostics: vec![format!("affine dimension={}", self.dimension)],
        })
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
            let evaluation = map.forward(&coordinates, &[])?;
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
pub struct SamplingMapComposition {
    kind: SamplingCompositionKind,
    children: Vec<Box<dyn SamplingMapComponent>>,
    dimensions: usize,
    output_dimensions: usize,
}

/// A direct-product map whose child output blocks are embedded into an
/// explicitly supplied master-frame permutation.  This is the exact generic
/// construction used for a surface subspace together with its complement:
/// each child remains in its native coordinates, while this wrapper only
/// permutes complete three-momentum blocks.  The permutation has unit
/// absolute determinant, so the product Jacobian is unchanged.
#[derive(Clone)]
pub struct SamplingMapEmbedding {
    map: Arc<SamplingMapComposition>,
    /// Product-output index -> master-frame output index.
    output_indices: Vec<usize>,
}

impl std::fmt::Debug for SamplingMapEmbedding {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("SamplingMapEmbedding")
            .field("map", &self.map)
            .field("output_indices", &self.output_indices)
            .finish()
    }
}

impl SamplingMapEmbedding {
    /// Build an embedded map from either a direct product or an ordered
    /// conditional composition.  The latter preserves the block context
    /// passed by `then(...)` while this wrapper only permutes its complete
    /// output into the master frame.
    pub fn from_composition(
        map: SamplingMapComposition,
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
            map: Arc::new(map),
            output_indices,
        })
    }

    /// Build an embedded direct product.  `output_indices` must be a complete
    /// permutation of `0..sum(child.output_dimensions())`; duplicate or
    /// missing frame coordinates are rejected before any numerical use.
    pub fn product(
        children: Vec<Box<dyn SamplingMapComponent>>,
        output_indices: Vec<usize>,
    ) -> Result<Self> {
        let map = SamplingMapComposition::product(children)?;
        Self::from_composition(map, output_indices)
    }

    pub fn output_indices(&self) -> &[usize] {
        &self.output_indices
    }

    fn embed_point(&self, product_point: &[f64]) -> Vec<f64> {
        let mut point = vec![0.0; product_point.len()];
        for (product_index, master_index) in self.output_indices.iter().enumerate() {
            point[*master_index] = product_point[product_index];
        }
        point
    }

    fn unembed_point(&self, master_point: &[f64]) -> Vec<f64> {
        self.output_indices
            .iter()
            .map(|master_index| master_point[*master_index])
            .collect()
    }
}

impl SamplingMapComponent for SamplingMapEmbedding {
    fn dimensions(&self) -> usize {
        self.map.dimensions()
    }

    fn output_dimensions(&self) -> usize {
        self.map.output_dimensions()
    }

    fn contract(&self) -> SamplingMapContract {
        self.map.contract()
    }

    fn name(&self) -> &'static str {
        "embedded_product"
    }

    fn forward(&self, coordinates: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        let mut evaluation = self.map.forward(coordinates, context)?;
        evaluation.point = self.embed_point(&evaluation.point);
        evaluation
            .diagnostics
            .push("unit-determinant master-frame embedding".to_owned());
        Ok(evaluation)
    }

    fn inverse(&self, point: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        if point.len() != self.output_dimensions() {
            return Err(eyre!(
                "sampling-map embedding inverse received output dimension {}, expected {}",
                point.len(),
                self.output_dimensions()
            ));
        }
        let mut evaluation = self.map.inverse(&self.unembed_point(point), context)?;
        // The inverse's point is the user-supplied master-frame point.
        evaluation.point = point.to_vec();
        evaluation
            .diagnostics
            .push("unit-determinant master-frame embedding".to_owned());
        Ok(evaluation)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SamplingCompositionKind {
    Product,
    Then,
}

impl std::fmt::Debug for SamplingMapComposition {
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

impl SamplingMapComposition {
    pub fn product(children: Vec<Box<dyn SamplingMapComponent>>) -> Result<Self> {
        Self::new(SamplingCompositionKind::Product, children)
    }

    pub fn then(children: Vec<Box<dyn SamplingMapComponent>>) -> Result<Self> {
        Self::new(SamplingCompositionKind::Then, children)
    }

    fn new(
        kind: SamplingCompositionKind,
        children: Vec<Box<dyn SamplingMapComponent>>,
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

    pub fn children(&self) -> &[Box<dyn SamplingMapComponent>] {
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
        coordinates: &[f64],
        initial_context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        self.validate_input(coordinates.len(), "composition")?;
        let mut offset = 0;
        let mut context = if self.is_then() {
            initial_context.to_vec()
        } else {
            Vec::new()
        };
        let mut evaluations = Vec::with_capacity(self.children.len());
        for child in &self.children {
            let end = offset + child.dimensions();
            let child_coordinates = &coordinates[offset..end];
            let evaluation = child.forward(
                child_coordinates,
                if self.is_then() { &context } else { &[] },
            )?;
            if self.is_then() {
                context.extend_from_slice(&evaluation.point);
            }
            evaluations.push(evaluation);
            offset = end;
        }
        let mut evaluation = combine_evaluations(evaluations);
        evaluation.support = self.contract().support;
        Ok(evaluation)
    }

    fn evaluate_inverse(
        &self,
        point: &[f64],
        initial_context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        self.validate_output(point.len())?;
        let mut offset = 0;
        let mut context = if self.is_then() {
            initial_context.to_vec()
        } else {
            Vec::new()
        };
        let mut evaluations = Vec::with_capacity(self.children.len());
        // In a triangular map, prior output blocks are known before later
        // blocks are inverted, so inverse traversal remains in declaration
        // order even though the determinant is block triangular.
        for child in &self.children {
            let end = offset + child.output_dimensions();
            let child_point = &point[offset..end];
            let evaluation =
                child.inverse(child_point, if self.is_then() { &context } else { &[] })?;
            if self.is_then() {
                context.extend_from_slice(child_point);
            }
            evaluations.push(evaluation);
            offset = end;
        }
        let mut evaluation = combine_evaluations(evaluations);
        evaluation.support = self.contract().support;
        Ok(evaluation)
    }
}

impl SamplingMapComponent for SamplingMapComposition {
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
                jacobian: SamplingJacobian::ExactForward,
            },
            combine_contracts,
        );
        if self.is_then() {
            // A conditional child is fully supported once a preceding full
            // block supplies its context. A conditional first block still
            // requires external preparation and remains conditional.
            contract.support = SamplingSupport::Full;
            for (index, child) in self.children.iter().enumerate() {
                match child.contract().support {
                    SamplingSupport::Branched => {
                        contract.support = SamplingSupport::Branched;
                    }
                    SamplingSupport::Conditional
                        if index == 0 && contract.support == SamplingSupport::Full =>
                    {
                        contract.support = SamplingSupport::Conditional;
                    }
                    _ => {}
                }
            }
        }
        contract
    }

    fn name(&self) -> &'static str {
        match self.kind {
            SamplingCompositionKind::Product => "product",
            SamplingCompositionKind::Then => "then",
        }
    }

    fn forward(&self, coordinates: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        self.evaluate_forward(coordinates, context)
    }

    fn inverse(&self, point: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        self.evaluate_inverse(point, context)
    }
}

pub(crate) fn combine_contracts(
    left: SamplingMapContract,
    right: SamplingMapContract,
) -> SamplingMapContract {
    let support = match (left.support, right.support) {
        (SamplingSupport::Branched, _) | (_, SamplingSupport::Branched) => {
            SamplingSupport::Branched
        }
        (SamplingSupport::Conditional, _) | (_, SamplingSupport::Conditional) => {
            SamplingSupport::Conditional
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
    SamplingMapContract { support, jacobian }
}

fn combine_evaluations(evaluations: Vec<SamplingMapEvaluation>) -> SamplingMapEvaluation {
    let mut coordinates = Vec::new();
    let mut point = Vec::new();
    let mut jacobian = 1.0;
    let mut inverse_jacobian = 1.0;
    let mut residual = 0.0_f64;
    let mut support = SamplingSupport::Full;
    let mut diagnostics = Vec::new();
    for evaluation in evaluations {
        coordinates.extend(evaluation.coordinates);
        point.extend(evaluation.point);
        jacobian *= evaluation.jacobian;
        inverse_jacobian *= evaluation.inverse_jacobian;
        residual = residual.max(evaluation.residual);
        support = match (support, evaluation.support) {
            (SamplingSupport::Branched, _) | (_, SamplingSupport::Branched) => {
                SamplingSupport::Branched
            }
            (SamplingSupport::Conditional, _) | (_, SamplingSupport::Conditional) => {
                SamplingSupport::Conditional
            }
            _ => SamplingSupport::Full,
        };
        diagnostics.extend(evaluation.diagnostics);
    }
    SamplingMapEvaluation {
        coordinates,
        point,
        jacobian,
        inverse_jacobian,
        residual,
        support,
        diagnostics,
    }
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
pub struct SurfaceRadialMap {
    dimension: usize,
    center: Vec<f64>,
    threshold_radius: Option<f64>,
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
pub type ImplicitSurfaceRadialEvaluator =
    Arc<dyn Fn(&[f64], f64) -> Result<(f64, f64)> + Send + Sync + 'static>;

/// Context-aware variant used by conditional surface charts. The context is
/// the output of an earlier `then(...)`/complement map and may contain the
/// sampled spectator loop momenta needed to solve a partial-subspace surface.
pub type ImplicitSurfaceRadialContextEvaluator =
    Arc<dyn Fn(&[f64], f64, &[f64]) -> Result<(f64, f64)> + Send + Sync + 'static>;

/// Context-dependent centre of a directional energy-surface chart.
///
/// The context is the output of an earlier ordered map (typically a sampled
/// complement block).  Keeping the centre evaluator separate from the scalar
/// surface evaluator makes the dependency explicit: derivatives with respect
/// to the active radial coordinates remain block diagonal, while centre and
/// root changes with the complement occupy only the off-diagonal Jacobian
/// block.  The returned vector must have the map's active dimension.
pub type ImplicitSurfaceCenterEvaluator =
    Arc<dyn Fn(&[f64]) -> Result<Vec<f64>> + Send + Sync + 'static>;

#[derive(Clone)]
pub struct ImplicitSurfaceRadialMap {
    dimension: usize,
    center: Vec<f64>,
    beta: f64,
    power: f64,
    evaluator: ImplicitSurfaceRadialEvaluator,
    context_evaluator: Option<ImplicitSurfaceRadialContextEvaluator>,
    center_evaluator: Option<ImplicitSurfaceCenterEvaluator>,
    root_tolerance: f64,
}

impl std::fmt::Debug for ImplicitSurfaceRadialMap {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ImplicitSurfaceRadialMap")
            .field("dimension", &self.dimension)
            .field("center", &self.center)
            .field("beta", &self.beta)
            .field("power", &self.power)
            .field("context_dependent", &self.context_evaluator.is_some())
            .field("center_context_dependent", &self.center_evaluator.is_some())
            .field("root_tolerance", &self.root_tolerance)
            .finish_non_exhaustive()
    }
}

impl ImplicitSurfaceRadialMap {
    pub fn new(
        dimension: usize,
        center: Vec<f64>,
        beta: f64,
        power: f64,
        evaluator: ImplicitSurfaceRadialEvaluator,
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
            center_evaluator: None,
            root_tolerance: 1.0e-11,
        })
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
        evaluator: ImplicitSurfaceRadialContextEvaluator,
    ) -> Self {
        self.context_evaluator = Some(evaluator);
        self
    }

    /// Attach a context-dependent centre to the chart.  The map remains
    /// exact in its active coordinates because the centre depends only on
    /// already sampled context, and therefore contributes no extra diagonal
    /// determinant factor in an ordered composition.
    pub fn with_context_center_evaluator(
        mut self,
        evaluator: ImplicitSurfaceCenterEvaluator,
    ) -> Self {
        self.center_evaluator = Some(evaluator);
        self
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn center(&self) -> &[f64] {
        &self.center
    }

    pub fn beta(&self) -> f64 {
        self.beta
    }

    pub fn power(&self) -> f64 {
        self.power
    }

    fn center_for_context(&self, context: &[f64]) -> Result<Vec<f64>> {
        let center = if let Some(evaluator) = &self.center_evaluator {
            evaluator(context)?
        } else {
            self.center.clone()
        };
        if center.len() != self.dimension {
            return Err(eyre!(
                "implicit surface radial-map centre has dimension {}, expected {}",
                center.len(),
                self.dimension
            ));
        }
        if center.iter().any(|value| !value.is_finite()) {
            return Err(eyre!(
                "implicit surface radial-map context centre must contain only finite values"
            ));
        }
        Ok(center)
    }

    pub fn contract(&self) -> SamplingMapContract {
        SamplingMapContract {
            support: if self.context_evaluator.is_some() || self.center_evaluator.is_some() {
                SamplingSupport::Conditional
            } else {
                SamplingSupport::Full
            },
            jacobian: SamplingJacobian::ExactImplicit,
        }
    }

    pub fn forward(&self, coordinates: &[f64]) -> Result<SamplingMapEvaluation> {
        self.forward_with_context(coordinates, &[])
    }

    pub fn forward_with_context(
        &self,
        coordinates: &[f64],
        context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
        self.validate_coordinates(coordinates)?;
        let (direction, angular_jacobian) = direction_from_coordinates(coordinates)?;
        let center = self.center_for_context(context)?;
        let root = self.root_for_direction(&direction, context)?;
        let (radius, radial_jacobian) = self.radius_from_coordinate(coordinates[0], root);
        if self.power != 1.0 && root.is_some_and(|(threshold, _)| radius == threshold) {
            return Err(eyre!(
                "implicit surface radial-map radius rounds to the singular threshold seam at the current precision"
            ));
        }
        let point = center
            .iter()
            .zip(&direction)
            .map(|(center, direction)| center + radius * direction)
            .collect::<Vec<_>>();
        let jacobian = radial_jacobian * angular_jacobian * radius.powi(self.dimension as i32 - 1);
        let inverse_jacobian = 1.0 / jacobian;
        if !radius.is_finite()
            || radius <= 0.0
            || point.iter().any(|component| !component.is_finite())
            || [jacobian, inverse_jacobian]
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(eyre!(
                "implicit surface radial-map point or positive Jacobian pair ({jacobian}, {inverse_jacobian}) is not representable at the current precision"
            ));
        }
        Ok(SamplingMapEvaluation {
            coordinates: coordinates.to_vec(),
            point,
            jacobian,
            inverse_jacobian,
            residual: 0.0,
            support: self.contract().support,
            diagnostics: vec![if root.is_some() {
                "implicit_surface:regular_root".to_owned()
            } else {
                "implicit_surface:absent_fallback".to_owned()
            }],
        })
    }

    pub fn inverse(&self, point: &[f64]) -> Result<SamplingMapEvaluation> {
        self.inverse_with_context(point, &[])
    }

    pub fn inverse_with_context(
        &self,
        point: &[f64],
        context: &[f64],
    ) -> Result<SamplingMapEvaluation> {
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
        let center = self.center_for_context(context)?;
        let displacement = point
            .iter()
            .zip(&center)
            .map(|(point, center)| point - center)
            .collect::<Vec<_>>();
        let radius = displacement
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        if !radius.is_finite() || radius <= 0.0 {
            return Err(eyre!(
                "implicit surface radial-map inverse is undefined at its centre"
            ));
        }
        let direction = displacement
            .iter()
            .map(|value| value / radius)
            .collect::<Vec<_>>();
        let root = self.root_for_direction(&direction, context)?;
        let coordinate = self.coordinate_from_radius(radius, root);
        let mut coordinates = coordinates_from_direction(&direction)?;
        coordinates[0] = coordinate;
        let mapped = self.forward_with_context(&coordinates, context)?;
        let residual = max_coordinate_residual_f64(point, &mapped.point);
        if !residual.is_finite() {
            return Err(eyre!(
                "implicit surface radial-map inverse residual is not representable at the current precision"
            ));
        }
        Ok(SamplingMapEvaluation {
            coordinates,
            point: point.to_vec(),
            jacobian: mapped.jacobian,
            inverse_jacobian: mapped.inverse_jacobian,
            residual,
            support: self.contract().support,
            diagnostics: mapped.diagnostics,
        })
    }

    fn validate_coordinates(&self, coordinates: &[f64]) -> Result<()> {
        if coordinates.len() != self.dimension {
            return Err(eyre!(
                "implicit surface radial-map forward received {}, expected {} coordinates",
                coordinates.len(),
                self.dimension
            ));
        }
        if coordinates
            .iter()
            .any(|coordinate| !coordinate.is_finite() || *coordinate <= 0.0 || *coordinate >= 1.0)
        {
            return Err(eyre!(
                "implicit surface radial-map coordinates must be finite and strictly inside the unit cube"
            ));
        }
        Ok(())
    }

    fn evaluate(&self, direction: &[f64], radius: f64, context: &[f64]) -> Result<(f64, f64)> {
        if let Some(evaluator) = &self.context_evaluator {
            evaluator(direction, radius, context)
        } else {
            (self.evaluator)(direction, radius)
        }
    }

    fn root_for_direction(&self, direction: &[f64], context: &[f64]) -> Result<Option<(f64, f64)>> {
        let (origin_value, origin_derivative) = self.evaluate(direction, 0.0, context)?;
        if !origin_value.is_finite() || !origin_derivative.is_finite() {
            return Err(eyre!(
                "implicit surface evaluator returned non-finite origin data"
            ));
        }
        let scale = origin_value.abs().max(1.0);
        let residual_tolerance = self.root_tolerance * scale;
        if origin_value >= 0.0 {
            return Ok(None);
        }
        if origin_value >= -residual_tolerance {
            return Err(eyre!(
                "implicit surface interior center is not certified at the current root tolerance: value={origin_value}, tolerance={residual_tolerance}"
            ));
        }
        let mut lower = 0.0;
        let mut upper = self.beta.max(1.0);
        let mut upper_value = None;
        for _ in 0..96 {
            let (value, derivative) = self.evaluate(direction, upper, context)?;
            if !value.is_finite() || !derivative.is_finite() {
                return Err(eyre!(
                    "implicit surface evaluator returned non-finite bracket data"
                ));
            }
            if value >= 0.0 {
                upper_value = Some((value, derivative));
                break;
            }
            upper *= 2.0;
            if !upper.is_finite() {
                return Err(eyre!(
                    "implicit surface radial root could not be bracketed before the radius overflowed"
                ));
            }
        }
        let Some((_, _)) = upper_value else {
            return Err(eyre!(
                "implicit surface radial root could not be bracketed after 96 radius expansions"
            ));
        };
        let mut value_at_lower = origin_value;
        let mut root = 0.5 * (lower + upper);
        for _ in 0..128 {
            let (value, derivative) = self.evaluate(direction, root, context)?;
            if !value.is_finite() || !derivative.is_finite() {
                return Err(eyre!(
                    "implicit surface evaluator returned non-finite root data"
                ));
            }
            if value.abs() <= residual_tolerance {
                break;
            }
            if value > 0.0 {
                upper = root;
            } else {
                lower = root;
                value_at_lower = value;
            }
            let newton = if derivative.abs() > f64::MIN_POSITIVE {
                root - value / derivative
            } else {
                f64::NAN
            };
            root = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * (lower + upper)
            };
        }
        let (value, derivative) = self.evaluate(direction, root, context)?;
        if !value.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
            return Err(eyre!(
                "implicit surface radial root is not regular: value={value}, derivative={derivative}, bracket_lower={value_at_lower}"
            ));
        }
        if !root.is_finite() || root <= 0.0 || value.abs() > residual_tolerance {
            return Err(eyre!(
                "implicit surface radial root failed residual certification: radius={root}, value={value}, tolerance={residual_tolerance}, bracket=[{lower}, {upper}]"
            ));
        }
        Ok(Some((root, derivative)))
    }

    fn radius_from_coordinate(&self, coordinate: f64, root: Option<(f64, f64)>) -> (f64, f64) {
        let (radius, jacobian) = SurfaceRadialMap::radius_from_coordinate(
            &F(coordinate),
            F(root.map(|(radius, _)| radius).unwrap_or(0.0)),
            F(self.beta),
            F(self.power),
        );
        (radius.0, jacobian.0)
    }

    fn coordinate_from_radius(&self, radius: f64, root: Option<(f64, f64)>) -> f64 {
        SurfaceRadialMap::coordinate_from_radius(
            &F(radius),
            F(root.map(|(radius, _)| radius).unwrap_or(0.0)),
            F(self.beta),
            F(self.power),
        )
        .0
    }
}

fn direction_from_coordinates(coordinates: &[f64]) -> Result<(Vec<f64>, f64)> {
    if coordinates.len() < 2 {
        return Err(eyre!(
            "spherical direction requires at least two coordinates"
        ));
    }
    let dimension = coordinates.len();
    let mut direction = vec![0.0; dimension];
    let mut angular_jacobian = std::f64::consts::TAU;
    let mut base = 1.0;
    for (i, coordinate) in coordinates[2..].iter().enumerate() {
        let cos_theta = -1.0 + 2.0 * coordinate;
        if cos_theta <= -1.0 || cos_theta >= 1.0 {
            return Err(eyre!(
                "polar coordinate rounds to a singular spherical chart boundary at the current precision"
            ));
        }
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        angular_jacobian *= 2.0;
        let angular_power = dimension - 3 - i;
        if angular_power > 0 {
            angular_jacobian *= sin_theta.powi(angular_power as i32);
        }
        direction[i] = base * cos_theta;
        base *= sin_theta;
    }
    let phi = std::f64::consts::TAU * coordinates[1];
    direction[dimension - 2] = base * phi.cos();
    direction[dimension - 1] = base * phi.sin();
    Ok((direction, angular_jacobian))
}

fn coordinates_from_direction(direction: &[f64]) -> Result<Vec<f64>> {
    if direction.len() < 2 {
        return Err(eyre!(
            "spherical direction requires at least two components"
        ));
    }
    let dimension = direction.len();
    let mut coordinates = vec![0.0; dimension];
    let mut base = 1.0;
    for i in 0..dimension - 2 {
        if base <= 0.0 {
            return Err(eyre!(
                "direction lies on a singular spherical chart boundary"
            ));
        }
        let cos_theta = direction[i] / base;
        if !cos_theta.is_finite() || cos_theta <= -1.0 || cos_theta >= 1.0 {
            return Err(eyre!(
                "direction lies on a singular spherical chart boundary at the current precision"
            ));
        }
        coordinates[2 + i] = (1.0 + cos_theta) / 2.0;
        base *= (1.0 - cos_theta * cos_theta).sqrt();
    }
    let mut phi = direction[dimension - 1].atan2(direction[dimension - 2]);
    if phi < 0.0 {
        phi += std::f64::consts::TAU;
    }
    coordinates[1] = phi / std::f64::consts::TAU;
    if !coordinates[1].is_finite() || coordinates[1] <= 0.0 || coordinates[1] >= 1.0 {
        return Err(eyre!(
            "direction lies on the spherical azimuth seam at the current precision"
        ));
    }
    Ok(coordinates)
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

impl SurfaceRadialMap {
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
        center: Vec<f64>,
        threshold_radius: Option<f64>,
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
        if threshold_radius.is_some_and(|radius| !radius.is_finite() || radius < 0.0) {
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
    pub fn absent(dimension: usize, center: Vec<f64>, beta: f64, power: f64) -> Result<Self> {
        Self::new(dimension, center, None, beta, power)
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn center(&self) -> &[f64] {
        &self.center
    }

    pub fn threshold_radius(&self) -> Option<f64> {
        self.threshold_radius
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
            jacobian: SamplingJacobian::ExactForward,
        }
    }

    /// Map a unit-cube point to a point around the surface centre.
    pub fn forward<T: FloatLike>(&self, coordinates: &[F<T>]) -> Result<SurfaceRadialPoint<T>> {
        self.validate_coordinates(coordinates)?;
        let (threshold, beta, power) = self.radial_parameters();
        let (radius, radial_jacobian) =
            Self::radius_from_coordinate(&coordinates[0], threshold.clone(), beta, power);
        if self.power != 1.0 && threshold > threshold.zero() && radius == threshold {
            return Err(eyre!(
                "surface radial-map radius rounds to the singular threshold seam at the current precision"
            ));
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
                return Err(eyre!(
                    "polar coordinate rounds to a singular spherical chart boundary at the current precision"
                ));
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
            point.push(component + F::from_f64(*centre));
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
            return Err(eyre!(
                "surface radial-map point or positive Jacobian pair is not representable at the current precision"
            ));
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
    pub fn inverse<T: FloatLike>(&self, point: &[F<T>]) -> Result<SurfaceRadialPoint<T>> {
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
            displacement.push(component - F::from_f64(*centre));
        }
        let mut radius_squared = displacement[0].square();
        for component in &displacement[1..] {
            radius_squared += component.square();
        }
        let radius = radius_squared.sqrt();
        let zero = radius.zero();
        if !radius.0.is_finite() || radius <= zero {
            return Err(eyre!(
                "surface radial-map inverse radius must be finite and positive; it is undefined at its centre"
            ));
        }
        let (threshold, beta, power) = self.radial_parameters();
        let radial_coordinate = Self::coordinate_from_radius(&radius, threshold, beta, power);
        let one = radius.one();
        let two = radius.from_i64(2);
        let mut coordinates = vec![zero.clone(); self.dimension];
        let mut base = radius.clone();
        for i in 0..self.dimension - 2 {
            let cos_theta = &displacement[i] / &base;
            coordinates[2 + i] = (&one + &cos_theta) / &two;
            base *= (&one - cos_theta.square()).sqrt();
        }
        let mut phi = F(displacement[self.dimension - 1]
            .0
            .atan2(&displacement[self.dimension - 2].0));
        if phi < zero {
            phi += radius.TAU();
        }
        coordinates[0] = radial_coordinate;
        coordinates[1] = phi / radius.TAU();
        let mapped = self.forward(&coordinates)?;
        let residual = max_coordinate_residual(point, &mapped.point);
        if !residual.0.is_finite() {
            return Err(eyre!(
                "surface radial-map inverse residual is not representable at the current precision"
            ));
        }
        Ok(SurfaceRadialPoint {
            coordinates,
            point: point.to_vec(),
            radius,
            jacobian: mapped.jacobian,
            inverse_jacobian: mapped.inverse_jacobian,
            residual,
        })
    }

    fn validate_coordinates<T: FloatLike>(&self, coordinates: &[F<T>]) -> Result<()> {
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

    fn radial_parameters<T: FloatLike>(&self) -> (F<T>, F<T>, F<T>) {
        (
            F::from_f64(self.threshold_radius.unwrap_or(0.0)),
            F::from_f64(self.beta),
            F::from_f64(self.power),
        )
    }

    fn radius_from_coordinate<T: FloatLike>(
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

    fn coordinate_from_radius<T: FloatLike>(
        radius: &F<T>,
        threshold: F<T>,
        beta: F<T>,
        power: F<T>,
    ) -> F<T> {
        let one = radius.one();
        let split = &threshold / (&threshold + &beta);
        if split > radius.zero() && radius < &threshold {
            &split * Self::power_complement(&(radius / &threshold), &(&one / &power))
        } else {
            let z = ((radius - &threshold) / &beta).powf(&(&one / &power));
            (&z + &split) / (&one + &z)
        }
    }

    /// Evaluate `1 - (1 - fraction)^power` without cancellation at the origin.
    /// The existing native hyperbolic operations preserve small arguments;
    /// direct subtraction is safe once its result is bounded away from zero.
    /// These native operations are not yet registered as eager primitives.
    fn power_complement<T: FloatLike>(fraction: &F<T>, power: &F<T>) -> F<T> {
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
        if let SamplingMapDefinition::Lmb(edges) = &definition {
            if edges.len() != n_loop_momenta {
                return Err(eyre!(
                    "lmb definition has {} edges, but the kernel has {} loop coordinates",
                    edges.len(),
                    n_loop_momenta
                ));
            }
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
            return Err(eyre!(
                "sampling-map forward point or positive Jacobian is not representable at the current precision"
            ));
        }
        let loop_momenta = LoopMomenta(
            raw.into_iter()
                .map(|p| ThreeMomentum::new(p[0].clone(), p[1].clone(), p[2].clone()))
                .collect(),
        );
        let (inverse_coordinates, inverse_jacobian) =
            global_inv_parameterize(&loop_momenta.0, e_cm, &self.settings);
        self.validate_coordinates(&inverse_coordinates)?;
        if !inverse_jacobian.0.is_finite() || inverse_jacobian <= inverse_jacobian.zero() {
            return Err(eyre!(
                "sampling-map inverse Jacobian must be finite and positive at the current precision"
            ));
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
        self.validate_coordinates(&coordinates)?;
        if !inverse_jacobian.0.is_finite() || inverse_jacobian <= inverse_jacobian.zero() {
            return Err(eyre!(
                "sampling-map inverse Jacobian must be finite and positive at the current precision"
            ));
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
            return Err(eyre!(
                "sampling-map inverse point, positive Jacobian or residual is not representable at the current precision"
            ));
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

impl SamplingMapComponent for SamplingMapKernel {
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

    fn forward(&self, coordinates: &[f64], _context: &[f64]) -> Result<SamplingMapEvaluation> {
        let coordinates = coordinates.iter().copied().map(F).collect::<Vec<_>>();
        let evaluation = SamplingMapKernel::forward(self, &coordinates)?;
        let point = evaluation
            .loop_momenta
            .0
            .iter()
            .flat_map(|momentum| [momentum.px.0, momentum.py.0, momentum.pz.0])
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

    fn inverse(&self, point: &[f64], _context: &[f64]) -> Result<SamplingMapEvaluation> {
        if point.len() != self.output_dimensions() {
            return Err(eyre!(
                "lmb sampling-map inverse received output dimension {}, expected {}",
                point.len(),
                self.output_dimensions()
            ));
        }
        let loop_momenta = LoopMomenta(
            point
                .chunks_exact(3)
                .map(|components| {
                    ThreeMomentum::new(F(components[0]), F(components[1]), F(components[2]))
                })
                .collect(),
        );
        let evaluation = SamplingMapKernel::inverse(self, &loop_momenta)?;
        let mapped_point = evaluation
            .loop_momenta
            .0
            .iter()
            .flat_map(|momentum| [momentum.px.0, momentum.py.0, momentum.pz.0])
            .collect();
        Ok(SamplingMapEvaluation {
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
        })
    }
}

impl SamplingMapComponent for SurfaceRadialMap {
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

    fn forward(&self, coordinates: &[f64], _context: &[f64]) -> Result<SamplingMapEvaluation> {
        let coordinates = coordinates.iter().copied().map(F).collect::<Vec<_>>();
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

    fn inverse(&self, point: &[f64], _context: &[f64]) -> Result<SamplingMapEvaluation> {
        let point_f = point.iter().copied().map(F).collect::<Vec<_>>();
        let evaluation = SurfaceRadialMap::inverse(self, &point_f)?;
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
}

impl SamplingMapComponent for ImplicitSurfaceRadialMap {
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

    fn forward(&self, coordinates: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        ImplicitSurfaceRadialMap::forward_with_context(self, coordinates, context)
    }

    fn inverse(&self, point: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
        ImplicitSurfaceRadialMap::inverse_with_context(self, point, context)
    }
}

fn max_coordinate_residual<T: FloatLike>(a: &[F<T>], b: &[F<T>]) -> F<T> {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .max_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or_else(|| F::<T>::from_f64(0.0))
}

fn max_coordinate_residual_f64(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
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
                        .radius_from_coordinate(coordinate, threshold.map(|radius| (radius, 1.0)));
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
                true,
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
    fn surface_radial_profiles_preserve_tiny_native_coordinates() {
        fn check<T: FloatLike>(decimal_exponent: i32) {
            let one = F::<T>::default().one();
            let tiny = &one / one.from_usize(10).powi(decimal_exponent);
            let tolerance = one.epsilon() * one.from_usize(1024);
            for power in [1.0, 2.0, 3.0] {
                let map = SurfaceRadialMap::new(3, vec![0.0; 3], Some(3.0), 2.0, power).unwrap();
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

        let forward = map.forward(&[0.2, -0.4], &[]).unwrap();
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
        assert!(SamplingMapAffine::new(vec![], vec![]).is_err());
        assert!(SamplingMapAffine::new(vec![vec![1.0, 0.0]], vec![0.0]).is_err());
        assert!(SamplingMapAffine::new(vec![vec![1.0, 0.0], vec![0.0, 1.0]], vec![0.0]).is_err());
        assert!(
            SamplingMapAffine::new(vec![vec![1.0, 2.0], vec![2.0, 4.0]], vec![0.0, 0.0],).is_err()
        );
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
        assert!(
            (point.jacobian.clone() * point.inverse_jacobian.clone() - F(1.0)).abs() < F(1.0e-12)
        );
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
        let shallow = ImplicitSurfaceRadialMap::new(
            2,
            vec![0.0; 2],
            1.0,
            1.0,
            Arc::new(|_, radius| Ok((radius - 1.0e-14, 1.0))),
        )
        .unwrap();
        assert!(
            shallow
                .forward(&[0.3, 0.4])
                .unwrap_err()
                .to_string()
                .contains("interior center is not certified")
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
        assert!(error.to_string().contains("could not be bracketed"));
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
        .with_context_evaluator(Arc::new(|_, radius, context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((radius - (1.0 + shift), 1.0))
        }));
        assert_eq!(map.contract().support, SamplingSupport::Conditional);
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
        .with_context_center_evaluator(Arc::new(|context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok(vec![shift, 0.0, 0.0])
        }))
        .with_context_evaluator(Arc::new(|_, radius, context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((radius - (1.0 + shift), 1.0))
        }));

        assert_eq!(map.contract().support, SamplingSupport::Conditional);
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
        .with_context_center_evaluator(Arc::new(|_| Ok(vec![0.0, 0.0])));
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
        .with_context_center_evaluator(Arc::new(|_| Ok(vec![f64::NAN, 0.0, 0.0])));
        let error = nonfinite
            .forward_with_context(&[0.2, 0.3, 0.7], &[])
            .expect_err("non-finite context centre must fail");
        assert!(
            error
                .to_string()
                .contains("must contain only finite values")
        );
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
        .with_context_center_evaluator(Arc::new(|context| {
            Ok(vec![context.first().copied().unwrap_or(0.0), 0.0, 0.0])
        }))
        .with_context_evaluator(Arc::new(|_, radius, context| {
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
        let coordinates = [0.5, 0.23, 0.31, 0.67];
        let forward = composition
            .forward(&coordinates, &[])
            .expect("forward ordered composition");
        let inverse = composition
            .inverse(&forward.point, &[])
            .expect("inverse ordered composition");
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
        let mapped = composition.forward(&coordinates, &[]).unwrap();
        assert_eq!(mapped.coordinates.len(), 5);
        assert_eq!(mapped.point.len(), 5);
        assert!((mapped.jacobian - first_jacobian * second_jacobian).abs() < 1.0e-10);
        assert!((mapped.jacobian * mapped.inverse_jacobian - 1.0).abs() < 1.0e-10);
        assert!(mapped.residual < 1.0e-11);
        let inverse = composition.inverse(&mapped.point, &[]).unwrap();
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
        let first_eval =
            <SurfaceRadialMap as SamplingMapComponent>::forward(&first, &[0.23, 0.71], &[])
                .unwrap();
        let second_eval =
            <SamplingMapKernel as SamplingMapComponent>::forward(&second, &[0.17, 0.43, 0.89], &[])
                .unwrap();
        let embedding = SamplingMapEmbedding::product(
            vec![Box::new(first), Box::new(second)],
            vec![2, 3, 0, 1, 4],
        )
        .expect("valid master-frame permutation");
        let coordinates = [0.23, 0.71, 0.17, 0.43, 0.89];
        let mapped = embedding.forward(&coordinates, &[]).unwrap();
        assert_eq!(mapped.point[0..2], second_eval.point[0..2]);
        assert_eq!(mapped.point[2..4], first_eval.point[0..2]);
        assert_eq!(mapped.point[4], second_eval.point[2]);
        assert!((mapped.jacobian - first_eval.jacobian * second_eval.jacobian).abs() < 1.0e-10);
        let inverse = embedding.inverse(&mapped.point, &[]).unwrap();
        assert!(inverse.residual < 1.0e-10);
        assert!(
            inverse
                .coordinates
                .iter()
                .zip(coordinates)
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-10)
        );
    }

    #[derive(Debug)]
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
                support: SamplingSupport::Conditional,
                jacobian: SamplingJacobian::ExactForward,
            }
        }

        fn name(&self) -> &'static str {
            "context_shift"
        }

        fn forward(&self, coordinates: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
            let shift = context.first().copied().unwrap_or(0.0);
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

        fn inverse(&self, point: &[f64], context: &[f64]) -> Result<SamplingMapEvaluation> {
            let shift = context.first().copied().unwrap_or(0.0);
            let coordinate = point[0] - shift;
            Ok(SamplingMapEvaluation {
                coordinates: vec![coordinate],
                point: point.to_vec(),
                jacobian: 1.0,
                inverse_jacobian: 1.0,
                residual: 0.0,
                support: self.contract().support,
                diagnostics: vec![format!("shift={shift}")],
            })
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
        assert_eq!(composition.contract().support, SamplingSupport::Conditional);
        let mapped = composition.forward(&[0.2, 0.4], &[]).unwrap();
        assert!(
            mapped
                .point
                .iter()
                .zip([0.2, 0.6])
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-15)
        );
        assert_eq!(mapped.diagnostics, vec!["shift=0", "shift=0.2"]);
        let inverse = composition.inverse(&mapped.point, &[]).unwrap();
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
        let mapped = composition.forward(&[0.2, 0.4], &[0.1]).unwrap();
        assert!(
            mapped
                .point
                .iter()
                .zip([0.3, 0.5])
                .all(|(actual, expected)| (actual - expected).abs() < 1.0e-15)
        );
        let inverse = composition.inverse(&mapped.point, &[0.1]).unwrap();
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
        .with_context_evaluator(Arc::new(|_, radius, context| {
            let shift = context.first().copied().unwrap_or(0.0);
            Ok((radius - (1.0 + shift), 1.0))
        }));
        let composition =
            SamplingMapComposition::then(vec![Box::new(ContextShiftMap), Box::new(implicit)])
                .expect("valid conditional surface composition");
        let coordinates = [0.2, 0.3, 0.27, 0.61];
        let mapped = composition.forward(&coordinates, &[]).unwrap();
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
        assert_eq!(mapped.support, SamplingSupport::Conditional);
        let inverse = composition.inverse(&mapped.point, &[]).unwrap();
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
        assert!(SamplingMapComposition::product(Vec::new()).is_err());
        let composition = SamplingMapComposition::product(vec![Box::new(ContextShiftMap)])
            .expect("valid product");
        let error = composition
            .forward(&[0.2, 0.3], &[])
            .unwrap_err()
            .to_string();
        assert!(error.contains("received dimension 2, expected 1"));
        let error = composition
            .inverse(&[0.2, 0.3], &[])
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
