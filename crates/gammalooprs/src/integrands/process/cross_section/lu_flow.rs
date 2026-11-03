//! Fixed points of the spatial Local Unitarity flow of one cross-section graph.
//!
//! LU dilates the spatial loop momenta about a fixed point, `y = c + t (x - c)`, and
//! localizes every cut surface `eta(y) = 0` at its radial roots with Jacobian `t^{3L}`.
//! Cuts that cross at an overlap cancel locally only if the slopes
//! `t d_t eta = grad(eta).(y - c)` share one sign on all crossed surfaces, and a root is
//! simple and well conditioned only if its slope stays away from zero. Every physical cut
//! equation is oriented as `sum(emitted) - sum(absorbed) - E`, so a positive slope on every
//! cut point is a causal flow. Convex E-surfaces are positive about any interior center;
//! thermal H-surfaces are not convex, and neither the origin nor any single center need
//! orient all of them.
//!
//! For fixed external momenta the valid centers of one surface point form a half-space, so
//! a sampled set of cut points gives a linear program in the center. The search keeps one
//! center when it covers the thermally relevant cut points. Otherwise it adds centers for
//! the remaining points and splits the integrand with a smooth partition of unity in the
//! physical point `y`, evaluated at every root of every cut. A weight vanishes wherever its
//! center is not causal on a nearby cut surface, so cuts that cancel at an overlap see one
//! map and one common smooth factor.

use clarabel::{
    algebra::CscMatrix,
    solver::{DefaultSettingsBuilder, DefaultSolver, IPSolver, NonnegativeConeT, SolverStatus},
};
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use rand::{Rng, SeedableRng, rngs::SmallRng};
use symbolica::domains::{dual::HyperDual, float::Real};
use tracing::{debug, info};

use super::CrossSectionGraphTerm;
use crate::{
    cff::{
        radial::{EnergySurfaceRay, RayGradient},
        surface::{HybridSurfaceID, HybridSurfaceRef},
    },
    graph::LoopMomentumBasis,
    momentum::{
        Rotatable, Rotation, ThreeMomentum,
        sample::{ExternalFourMomenta, ExternalThreeMomenta, LoopMomenta},
    },
    settings::{RuntimeSettings, runtime::LuCenterSearchSettings},
    utils::{
        F, FloatLike, box_muller,
        hyperdual_utils::new_constant,
        newton_solver::{RadialRootDiagnostics, RadialRootIdentity},
    },
};

/// Weight numerator of the first patch where no patch is fully causal. It keeps the
/// partition exact there; the coverage check keeps such points thermally irrelevant.
const PATCH_FLOOR: f64 = 1.0e-9;
/// Sampled cut points per surface entering one linear program.
const LP_POINTS_PER_SURFACE: usize = 160;
/// Relevant sampled cut points per surface that are constrained and checked.
const POINTS_PER_SURFACE: usize = 1000;
/// Random base points of the ray sample, beside the origin and the centers found so far.
const RANDOM_BASES: usize = 4;
/// Center coordinates stay within this many energy scales of the origin.
const CENTER_BOX: f64 = 4.0;
/// Initial and final steps of the coverage pattern search, in energy scales.
const REFINE_STEP: f64 = 0.25;
const REFINE_TOLERANCE: f64 = 1.0e-3;

/// Prepared flow centers of one graph for fixed external kinematics.
#[derive(Debug, Clone)]
pub(crate) struct LuFlowPlan {
    /// Patch centers in the generation loop-momentum basis for the warm-up external
    /// momenta; `None` is the origin.
    centers: Vec<Option<LoopMomenta<F<f64>>>>,
    gate: PatchGate,
}

/// Smooth validity gate shared by the search and the runtime weights.
#[derive(Debug, Clone, Copy)]
struct PatchGate {
    /// Energy and momentum scale of the external kinematics.
    scale: f64,
    /// Required slope `grad(eta).(y - c)` near the origin, in energy units. A cut surface
    /// whose equation is within this value of zero constrains the weights at a point.
    slope: f64,
}

impl PatchGate {
    /// Slope threshold at `|y|^2`, decaying quadratically beyond the external scale as
    /// the slopes on unbounded H-surfaces do.
    fn threshold<T: FloatLike>(&self, radius_squared: &F<T>) -> F<T> {
        let scale = exact::<T>(self.scale);
        exact::<T>(self.slope) / (radius_squared.one() + radius_squared / (&scale * &scale))
    }

    fn threshold_dual<T: FloatLike>(&self, radius_squared: &HyperDual<F<T>>) -> HyperDual<F<T>> {
        let constant = |value: f64| new_constant(radius_squared, &exact::<T>(value));
        let scale = constant(self.scale);
        constant(self.slope) / (constant(1.0) + radius_squared.clone() / (scale.clone() * &scale))
    }

    /// Gate of a center at a point on or near a cut surface with the given causal slope.
    /// The search treats every nearby surface as close to the point.
    fn validity(&self, slope: f64, radius_squared: f64) -> f64 {
        gate(F(0.0), F(slope / self.threshold(&F(radius_squared)).0)).0
    }
}

/// Infinitely differentiable step from zero at `x <= 0` to one at `x >= 1`.
fn smooth_step<T: FloatLike>(x: F<T>) -> F<T> {
    let (zero, one) = (x.zero(), x.one());
    if x <= zero {
        zero
    } else if x >= one {
        one
    } else {
        let rising = (-x.inv()).exp();
        let falling = (-(one - &x).inv()).exp();
        rising.clone() / (rising + falling)
    }
}

fn smooth_step_dual<T: FloatLike>(x: HyperDual<F<T>>) -> HyperDual<F<T>> {
    let real = x.values[0].clone();
    if real <= real.zero() {
        new_constant(&x, &real.zero())
    } else if real >= real.one() {
        new_constant(&x, &real.one())
    } else {
        let one = new_constant(&x, &real.one());
        let rising = (-(one.clone() / &x)).exp();
        let falling = (-(one.clone() / (one - &x))).exp();
        rising.clone() / (rising + falling)
    }
}

/// Validity of one center at one point with respect to one surface, from the surface
/// value and the causal slope in units of their thresholds. It is one away from the surface
/// (`|near| >= 1`) or where the slope is causal (`causal >= 1`), and zero where the surface
/// is close (`|near| <= 1/2`) and the slope is below a quarter of its threshold.
fn gate<T: FloatLike>(near: F<T>, causal: F<T>) -> F<T> {
    let (one, three, four) = (near.one(), near.from_usize(3), near.from_usize(4));
    let close = one.clone() - smooth_step((near.clone() * &near * &four - &one) / &three);
    one.clone() - close * (one.clone() - smooth_step((causal * four - one) / three))
}

fn gate_dual<T: FloatLike>(near: HyperDual<F<T>>, causal: HyperDual<F<T>>) -> HyperDual<F<T>> {
    let constant = |value: usize| new_constant(&near, &near.values[0].from_usize(value));
    let close = constant(1)
        - smooth_step_dual((near.clone() * &near * constant(4) - constant(1)) / constant(3));
    constant(1)
        - close
            * (constant(1) - smooth_step_dual((causal * constant(4) - constant(1)) / constant(3)))
}

fn exact<T: FloatLike>(value: f64) -> F<T> {
    F(T::from_f64_exact_binary(value))
}

/// One cut surface of the graph in its physical orientation.
#[derive(Clone)]
struct FlowSurface<'a> {
    surface: HybridSurfaceRef<'a>,
    /// `+1` or `-1`: the physical equation is this multiple of the represented one.
    orientation: i64,
    absorbed: &'a [EdgeIndex],
}

/// Cut surfaces of one patch with their displacements to every patch center, whose slopes
/// set the patch weights at that patch's roots.
pub(crate) struct LuPatchGeometry<'a, T: FloatLike> {
    patch: usize,
    surfaces: Vec<PatchSurface<'a, T>>,
}

/// One cut surface on a patch ray: its physical orientation and the displacements from the
/// patch center to every center, routed through its energies.
struct PatchSurface<'a, T: FloatLike> {
    orientation: F<T>,
    ray: &'a EnergySurfaceRay<T>,
    displacements: Vec<Vec<ThreeMomentum<F<T>>>>,
}

impl LuFlowPlan {
    pub(crate) fn origin() -> Self {
        Self {
            centers: vec![None],
            gate: PatchGate {
                scale: 1.0,
                slope: 1.0,
            },
        }
    }

    /// One center repeated as identical patches, whose weights split every root.
    #[cfg(test)]
    pub(crate) fn repeated(center: LoopMomenta<F<f64>>, patches: usize) -> Self {
        Self {
            centers: vec![Some(center); patches],
            gate: PatchGate {
                scale: 1.0,
                slope: 0.05,
            },
        }
    }

    pub(crate) fn patch_count(&self) -> usize {
        self.centers.len()
    }

    /// Center of one patch for the current rotation and precision. The exact binary64
    /// embedding keeps every precision on the same flow.
    pub(crate) fn center<T: FloatLike>(
        &self,
        patch: usize,
        rotation: &Rotation,
    ) -> Option<LoopMomenta<F<T>>> {
        self.centers[patch].as_ref().map(|center| {
            center
                .iter()
                .map(|momentum| momentum.map_ref(&|value| exact(value.0)).rotate(rotation))
                .collect()
        })
    }

    pub(crate) fn centers<T: FloatLike>(
        &self,
        rotation: &Rotation,
    ) -> Vec<Option<LoopMomenta<F<T>>>> {
        (0..self.patch_count())
            .map(|patch| self.center(patch, rotation))
            .collect()
    }

    /// Route the displacements from the center of `patch` to every center through the
    /// cut surfaces of that patch, given as `(orientation, ray)`.
    pub(crate) fn patch_geometry<'a, T: FloatLike>(
        &self,
        patch: usize,
        centers: &[Option<LoopMomenta<F<T>>>],
        surfaces: impl Iterator<Item = (i64, &'a EnergySurfaceRay<T>)>,
        lmb: &LoopMomentumBasis,
        zero: &F<T>,
    ) -> LuPatchGeometry<'a, T> {
        let point = |center: &Option<LoopMomenta<F<T>>>| {
            center.clone().unwrap_or_else(|| {
                LoopMomenta::from_iter(
                    lmb.loop_edges
                        .iter()
                        .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
                )
            })
        };
        let own = point(&centers[patch]);
        let displacements = centers
            .iter()
            .map(|center| &own - &point(center))
            .collect_vec();
        LuPatchGeometry {
            patch,
            surfaces: surfaces
                .map(|(orientation, ray)| PatchSurface {
                    orientation: zero.from_i64(orientation),
                    ray,
                    displacements: displacements
                        .iter()
                        .map(|displacement| ray.route_displacement(lmb, displacement))
                        .collect(),
                })
                .collect(),
        }
    }

    /// Partition weight of the geometry's patch at its root `radius` and physical point.
    pub(crate) fn weight<T: FloatLike>(
        &self,
        geometry: &LuPatchGeometry<'_, T>,
        radius: &F<T>,
        point: &LoopMomenta<F<T>>,
    ) -> F<T> {
        let one = radius.one();
        let slope = exact::<T>(self.gate.slope);
        let threshold = self.gate.threshold(&point.hyper_radius_squared(None));
        let mut numerators = vec![one.clone(); self.patch_count()];
        // One evaluation per surface; every center only adds its displacement term.
        for surface in &geometry.surfaces {
            let RayGradient {
                value,
                radial,
                weights,
            } = surface.ray.evaluate_gradient(radius);
            let near = value / &slope;
            for (numerator, displacement) in numerators.iter_mut().zip(&surface.displacements) {
                let causal = weights
                    .iter()
                    .zip(displacement)
                    .fold(radial.clone(), |sum, (weight, d)| sum + weight.clone() * d);
                *numerator *= gate(
                    near.clone(),
                    surface.orientation.clone() * causal / &threshold,
                );
            }
        }
        Self::normalize(numerators, geometry.patch, exact(PATCH_FLOOR), one)
    }

    /// Jets of [`Self::weight`] along the ray, on the integrand's dual shape.
    pub(crate) fn weight_dual<T: FloatLike>(
        &self,
        geometry: &LuPatchGeometry<'_, T>,
        radius: &HyperDual<F<T>>,
        point: &LoopMomenta<HyperDual<F<T>>>,
    ) -> HyperDual<F<T>> {
        let constant = |value: f64| new_constant(radius, &exact::<T>(value));
        let slope = constant(self.gate.slope);
        let radius_squared = point
            .iter()
            .flat_map(|momentum| [&momentum.px, &momentum.py, &momentum.pz])
            .fold(constant(0.0), |sum, value| sum + value.clone() * value);
        let threshold = self.gate.threshold_dual(&radius_squared);
        let mut numerators = vec![constant(1.0); self.patch_count()];
        for surface in &geometry.surfaces {
            let RayGradient {
                value,
                radial,
                weights,
            } = surface.ray.evaluate_gradient_dual(radius);
            let near = value / &slope;
            let orientation = new_constant(radius, &surface.orientation);
            for (numerator, displacement) in numerators.iter_mut().zip(&surface.displacements) {
                let causal =
                    weights
                        .iter()
                        .zip(displacement)
                        .fold(radial.clone(), |sum, (weight, d)| {
                            sum + weight.clone() * &d.map_ref(&|value| new_constant(radius, value))
                        });
                *numerator = numerator.clone()
                    * gate_dual(near.clone(), orientation.clone() * causal / &threshold);
            }
        }
        Self::normalize(
            numerators,
            geometry.patch,
            constant(PATCH_FLOOR),
            constant(1.0),
        )
    }

    /// `w_k = phi_k / sum_l phi_l`, with a floor on the first patch that acts only where no
    /// patch is fully causal and keeps the partition exact there.
    fn normalize<N>(numerators: Vec<N>, patch: usize, floor: N, one: N) -> N
    where
        N: Clone
            + std::ops::Add<N, Output = N>
            + std::ops::Sub<N, Output = N>
            + std::ops::Mul<N, Output = N>
            + std::ops::Div<N, Output = N>,
    {
        let floor = numerators
            .iter()
            .fold(floor, |floor, value| floor * (one.clone() - value.clone()));
        let total = numerators
            .iter()
            .cloned()
            .fold(floor.clone(), |sum, value| sum + value);
        let own = numerators[patch].clone();
        if patch == 0 {
            (own + floor) / total
        } else {
            own / total
        }
    }
}

/// One sampled point of a cut surface with its local center constraints.
struct CutPoint {
    surface: usize,
    point: Vec<f64>,
    /// `-beta` times the absorbed bath energy.
    log_weight: f64,
    /// Boltzmann factor relative to the most relevant point of the first sample over all
    /// cut surfaces, set with the constraints.
    weight: f64,
    /// `(n, n.y - threshold)` for the causal half-space `n.(y - c) >= threshold` of the
    /// point's own surface, first, and of every other cut surface near it, with `n` the
    /// oriented gradient.
    rows: Vec<(Vec<f64>, f64)>,
}

impl CutPoint {
    fn covered_by(&self, center: &[f64]) -> bool {
        self.rows.iter().all(|(normal, bound)| {
            normal.iter().zip(center).map(|(n, c)| n * c).sum::<f64>() <= *bound
        })
    }

    /// Causal slopes `n.(y - c)` of the local constraints, own surface first.
    fn slopes<'a>(&'a self, center: &'a [f64]) -> impl Iterator<Item = f64> + 'a {
        self.rows.iter().map(move |(normal, _)| {
            normal
                .iter()
                .zip(&self.point)
                .zip(center)
                .map(|((n, y), c)| n * (y - c))
                .sum::<f64>()
        })
    }

    /// Product of the smooth gates of every local constraint.
    fn validity(&self, gate: &PatchGate, center: &[f64]) -> f64 {
        let radius_squared = self.point.iter().map(|x| x * x).sum::<f64>();
        self.slopes(center)
            .map(|slope| gate.validity(slope, radius_squared))
            .product()
    }
}

/// Preparation-time search of causal flow centers for one graph.
struct CenterSearch<'a> {
    graph: &'a str,
    lmb: &'a LoopMomentumBasis,
    externals: &'a ExternalFourMomenta<F<f64>>,
    masses: &'a EdgeVec<F<f64>>,
    surfaces: Vec<FlowSurface<'a>>,
    beta: f64,
    gate: PatchGate,
    settings: &'a LuCenterSearchSettings,
    rng: SmallRng,
}

impl CrossSectionGraphTerm {
    /// Prepare the flow centers of this graph for fixed external kinematics.
    pub(crate) fn prepare_lu_flow(
        &self,
        settings: &RuntimeSettings,
        externals: &ExternalFourMomenta<F<f64>>,
        masses: &EdgeVec<F<f64>>,
    ) -> Result<LuFlowPlan> {
        let external = externals
            .iter()
            .map(|momentum| {
                momentum
                    .temporal
                    .value
                    .0
                    .abs()
                    .max(momentum.spatial.norm_squared().0.sqrt())
            })
            .fold(0.0, f64::max);
        let scale = if external > 0.0 {
            external
        } else {
            settings.kinematics.e_cm
        };
        let gate = PatchGate {
            scale,
            slope: settings.lu_flow.search.slope * scale,
        };
        if let Some(rows) = settings.lu_flow.fixed_centers.get(&self.graph.name) {
            return Ok(LuFlowPlan {
                centers: vec![Some(self.fixed_lu_center(rows, externals)?)],
                gate,
            });
        }
        if settings.lu_flow.uses_origin(&self.graph.name) {
            return Ok(LuFlowPlan::origin());
        }
        let surfaces = self
            .cut_group_data
            .cut_groups
            .iter_enumerated()
            .filter(|(id, _)| self.counterterm.cut_group_is_active(*id))
            .map(|(_, group)| FlowSurface {
                surface: match group.residue.surface_ids[0] {
                    HybridSurfaceID::Esurface(_) => {
                        HybridSurfaceRef::Esurface(&self.cut_esurface[group.cuts[0]])
                    }
                    HybridSurfaceID::Hsurface(id) => {
                        HybridSurfaceRef::Hsurface(&self.hsurfaces[id])
                    }
                    _ => unreachable!("LU residue groups contain E- or H-surfaces"),
                },
                orientation: i64::from(group.surface_sign),
                absorbed: group.absorbed_edges(&self.hsurfaces),
            })
            .collect();
        Ok(CenterSearch {
            graph: &self.graph.name,
            lmb: &self.graph.loop_momentum_basis,
            externals,
            masses,
            surfaces,
            beta: settings.general.inverse_temperature,
            gate,
            settings: &settings.lu_flow.search,
            rng: SmallRng::seed_from_u64(settings.lu_flow.search.seed),
        }
        .run())
    }

    /// Materialize an explicit center in the generation LMB from its external coefficients.
    /// Edge IDs keep the prescription independent of the order of the configured rows.
    fn fixed_lu_center(
        &self,
        rows: &[(usize, Vec<f64>)],
        externals: &ExternalFourMomenta<F<f64>>,
    ) -> Result<LoopMomenta<F<f64>>> {
        let lmb = &self.graph.loop_momentum_basis;
        let supplied = rows.iter().map(|(edge, _)| *edge).sorted().collect_vec();
        let expected = lmb
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .sorted()
            .collect_vec();
        if supplied != expected {
            return Err(eyre!(
                "LU flow center for graph '{}' must specify each generation LMB edge {:?} exactly once; received {:?}",
                self.graph.name,
                expected,
                rows,
            ));
        }
        let center = lmb
            .loop_edges
            .iter()
            .map(|edge| {
                let coefficients = &rows.iter().find(|(id, _)| *id == edge.0).unwrap().1;
                if coefficients.len() != externals.len()
                    || coefficients
                        .iter()
                        .any(|coefficient| !coefficient.is_finite())
                {
                    return Err(eyre!(
                        "LU flow center for graph '{}' loop edge {} requires {} finite external coefficients; received {:?}",
                        self.graph.name,
                        edge.0,
                        externals.len(),
                        coefficients,
                    ));
                }
                Ok(coefficients.iter().zip(externals).fold(
                    ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
                    |sum, (coefficient, momentum)| sum + momentum.spatial * F(*coefficient),
                ))
            })
            .collect::<Result<LoopMomenta<_>>>()?;
        if center
            .iter()
            .any(|p| !p.px.0.is_finite() || !p.py.0.is_finite() || !p.pz.0.is_finite())
        {
            return Err(eyre!(
                "Nonfinite LU flow center for graph '{}'",
                self.graph.name
            ));
        }
        Ok(center)
    }
}

impl CenterSearch<'_> {
    fn dimension(&self) -> usize {
        3 * self.lmb.loop_edges.len()
    }

    fn loops(&self, coordinates: &[f64]) -> LoopMomenta<F<f64>> {
        coordinates
            .chunks(3)
            .map(|v| ThreeMomentum::new(F(v[0]), F(v[1]), F(v[2])))
            .collect()
    }

    fn ray(
        &self,
        surface: &FlowSurface,
        velocity: &[f64],
        center: &[f64],
    ) -> EnergySurfaceRay<f64> {
        let (velocity, center) = (self.loops(velocity), self.loops(center));
        match surface.surface {
            HybridSurfaceRef::Esurface(esurface) => {
                esurface.routed_ray(&velocity, &center, self.externals, self.masses, self.lmb)
            }
            HybridSurfaceRef::Hsurface(hsurface) => {
                hsurface.routed_ray(&velocity, &center, self.externals, self.masses, self.lmb)
            }
            _ => unreachable!("LU cut surfaces are E- or H-surfaces"),
        }
    }

    /// Value of a cut equation at a loop-space point.
    fn value(&self, surface: &FlowSurface, point: &[f64]) -> f64 {
        self.ray(surface, &vec![0.0; point.len()], point)
            .evaluate(&F(0.0))
            .0
            .0
    }

    /// Gradient of the physically oriented cut equation at a loop-space point.
    fn oriented_gradient(&self, surface: &FlowSurface, point: &[f64]) -> Vec<f64> {
        let zero = vec![0.0; point.len()];
        let ray = self.ray(surface, &zero, point);
        let weights = ray.evaluate_gradient(&F(0.0)).weights;
        (0..point.len())
            .map(|axis| {
                let mut direction = zero.clone();
                direction[axis] = 1.0;
                let routed = ray.route_displacement(self.lmb, &self.loops(&direction));
                surface.orientation as f64
                    * weights
                        .iter()
                        .zip(&routed)
                        .map(|(weight, d)| (*weight * d).0)
                        .sum::<f64>()
            })
            .collect()
    }

    /// Absorbed bath energy of a surface at a point.
    fn absorbed_energy(&self, surface: &FlowSurface, point: &[f64]) -> f64 {
        let spatial: ExternalThreeMomenta<F<f64>> = self
            .externals
            .iter()
            .map(|momentum| momentum.spatial)
            .collect();
        let loops = self.loops(point);
        surface
            .absorbed
            .iter()
            .map(|&edge| {
                let momentum = self.lmb.edge_signatures[edge].compute_momentum(&loops, &spatial);
                (momentum.norm_squared() + self.masses[edge] * self.masses[edge])
                    .sqrt()
                    .0
            })
            .sum()
    }

    fn direction(&mut self, length: f64) -> Vec<f64> {
        let mut gaussian = vec![];
        while gaussian.len() < self.dimension() {
            let (a, b) = box_muller(
                F(1.0 - self.rng.random::<f64>()),
                F(self.rng.random::<f64>()),
            );
            gaussian.extend([a.0, b.0]);
        }
        gaussian.truncate(self.dimension());
        let norm = gaussian.iter().map(|x| x * x).sum::<f64>().sqrt();
        gaussian.iter().map(|x| x * length / norm).collect()
    }

    /// Cut points from rays through each base point in random directions.
    fn sample(&mut self, bases: &[Vec<f64>]) -> Vec<CutPoint> {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new(format!("LU center search graph '{}'", self.graph));
        let mut points = vec![];
        for base in bases {
            for index in 0..self.surfaces.len() {
                for _ in 0..self.settings.rays {
                    let velocity = self.direction(self.gate.scale);
                    let ray = self.ray(&self.surfaces[index], &velocity, base);
                    // A tangent or unresolved ray is a boundary sample; neighbouring rays
                    // cover the surface around it.
                    let Ok(roots) =
                        ray.positive_radius_roots(&F(self.gate.scale), &mut diagnostics, &identity)
                    else {
                        continue;
                    };
                    for root in roots {
                        let point = base
                            .iter()
                            .zip(&velocity)
                            .map(|(b, v)| b + root.solution.0 * v)
                            .collect_vec();
                        let log_weight =
                            -self.beta * self.absorbed_energy(&self.surfaces[index], &point);
                        points.push(CutPoint {
                            surface: index,
                            point,
                            log_weight,
                            weight: 0.0,
                            rows: vec![],
                        });
                    }
                }
            }
        }
        points
    }

    /// Keep the thermally relevant points and attach their local causal constraints.
    fn constrain(&self, mut points: Vec<CutPoint>, reference: f64) -> Vec<CutPoint> {
        let cut = self.settings.relevance.ln();
        points.retain(|point| point.log_weight - reference >= cut);
        // Thin each surface evenly across base points and rays to bound the work.
        let mut counts = vec![0_usize; self.surfaces.len()];
        for point in &points {
            counts[point.surface] += 1;
        }
        let mut seen = vec![0_usize; self.surfaces.len()];
        points.retain(|point| {
            seen[point.surface] += 1;
            (seen[point.surface] - 1)
                .is_multiple_of(counts[point.surface].div_ceil(POINTS_PER_SURFACE))
        });
        for point in &mut points {
            point.weight = (point.log_weight - reference).exp();
            let radius_squared = point.point.iter().map(|x| x * x).sum::<f64>();
            let threshold = self.gate.threshold(&F(radius_squared)).0;
            let own = point.surface;
            point.rows = std::iter::once(own)
                .chain((0..self.surfaces.len()).filter(|&index| index != own))
                .filter(|&index| {
                    index == own
                        || self.value(&self.surfaces[index], &point.point).abs() < self.gate.slope
                })
                .map(|index| {
                    let gradient = self.oriented_gradient(&self.surfaces[index], &point.point);
                    let bound = gradient
                        .iter()
                        .zip(&point.point)
                        .map(|(n, y)| n * y)
                        .sum::<f64>()
                        - threshold;
                    (gradient, bound)
                })
                .collect();
        }
        points
    }

    /// Orthonormal loop-space directions that every rotation fixing the external spatial
    /// momenta leaves invariant. For one external direction each loop center lies on its
    /// axis; at rest only the origin is invariant.
    fn invariant_basis(&self) -> Vec<Vec<f64>> {
        let mut directions: Vec<[f64; 3]> = vec![];
        for momentum in self.externals {
            let mut v = [
                momentum.spatial.px.0,
                momentum.spatial.py.0,
                momentum.spatial.pz.0,
            ];
            for d in &directions {
                let projection = v.iter().zip(d).map(|(a, b)| a * b).sum::<f64>();
                v.iter_mut().zip(d).for_each(|(a, b)| *a -= projection * b);
            }
            let norm = v.iter().map(|a| a * a).sum::<f64>().sqrt();
            if norm > 1e-12 * self.gate.scale {
                directions.push(v.map(|a| a / norm));
            }
        }
        let dimension = self.dimension();
        (0..self.lmb.loop_edges.len())
            .flat_map(|loop_index| {
                directions.iter().map(move |direction| {
                    let mut vector = vec![0.0; dimension];
                    vector[3 * loop_index..3 * loop_index + 3].copy_from_slice(direction);
                    vector
                })
            })
            .collect()
    }

    /// Center in the span of `basis` covering the largest Boltzmann-weighted share of
    /// `points`: a weighted L1 relaxation of the causal half-spaces, solved as a linear
    /// program in the basis coefficients.
    fn cover(&self, points: &[&CutPoint], basis: &[Vec<f64>]) -> Option<Vec<f64>> {
        let dimension = self.dimension();
        if basis.is_empty() {
            return Some(vec![0.0; dimension]);
        }
        let (m, n) = (basis.len(), points.len());
        let (mut rows, mut columns, mut values, mut bounds) = (vec![], vec![], vec![], vec![]);
        let mut push = |row: usize, column: usize, value: f64| {
            rows.push(row);
            columns.push(column);
            values.push(value);
        };
        let mut row = 0;
        for (slack, point) in points.iter().enumerate() {
            for (normal, bound) in &point.rows {
                for (axis, direction) in basis.iter().enumerate() {
                    push(
                        row,
                        axis,
                        normal.iter().zip(direction).map(|(a, b)| a * b).sum(),
                    );
                }
                push(row, m + slack, -1.0);
                bounds.push(*bound);
                row += 1;
            }
            push(row, m + slack, -1.0);
            bounds.push(0.0);
            row += 1;
        }
        for axis in 0..m {
            for sign in [1.0, -1.0] {
                push(row, axis, sign);
                bounds.push(CENTER_BOX * self.gate.scale);
                row += 1;
            }
        }
        let largest = points.iter().map(|point| point.weight).fold(0.0, f64::max);
        let objective = vec![0.0; m]
            .into_iter()
            .chain(points.iter().map(|point| point.weight / largest))
            .collect_vec();
        let mut solver = DefaultSolver::new(
            &CscMatrix::spalloc((m + n, m + n), 0),
            &objective,
            &CscMatrix::new_from_triplets(row, m + n, rows, columns, values),
            &bounds,
            &[NonnegativeConeT(row)],
            DefaultSettingsBuilder::default()
                .verbose(false)
                .build()
                .expect("default Clarabel settings are valid"),
        )
        .expect("the LU center linear program is well formed");
        solver.solve();
        debug!(
            graph = self.graph,
            points = n,
            status = ?solver.solution.status,
            objective = solver.solution.obj_val,
            iterations = solver.solution.iterations,
            "LU center linear program"
        );
        matches!(
            solver.solution.status,
            SolverStatus::Solved | SolverStatus::AlmostSolved
        )
        .then(|| {
            (0..dimension)
                .map(|axis| {
                    basis
                        .iter()
                        .zip(&solver.solution.x[..m])
                        .map(|(direction, coefficient)| coefficient * direction[axis])
                        .sum()
                })
                .collect()
        })
    }

    /// Pattern search on the soft coverage within the span of `basis`, starting from the
    /// linear program's center: the relaxation minimizes the total violation, not the
    /// covered Boltzmann weight.
    fn refine(&self, points: &[&CutPoint], basis: &[Vec<f64>], mut center: Vec<f64>) -> Vec<f64> {
        let coverage = |center: &[f64]| {
            points
                .iter()
                .map(|point| point.weight * point.validity(&self.gate, center))
                .sum::<f64>()
        };
        let mut best = coverage(&center);
        let mut step = REFINE_STEP * self.gate.scale;
        while step > REFINE_TOLERANCE * self.gate.scale {
            let mut improved = false;
            for direction in basis {
                for sign in [1.0, -1.0] {
                    let trial = center
                        .iter()
                        .zip(direction)
                        .map(|(c, d)| c + sign * step * d)
                        .collect_vec();
                    let value = coverage(&trial);
                    if value > best {
                        (best, center, improved) = (value, trial, true);
                    }
                }
            }
            if !improved {
                step /= 2.0;
            }
        }
        center
    }

    /// Greedy cover of the relevant points by successive centers. The valid centers of
    /// every point form a convex set, and the rotations fixing the external momenta map the
    /// complete set of cut points to itself, so the first center is chosen in their
    /// invariant subspace without loss for the relaxation. Further patches are general.
    fn greedy(&self, points: &[CutPoint]) -> Vec<Vec<f64>> {
        let invariant = self.invariant_basis();
        let full = (0..self.dimension())
            .map(|axis| {
                let mut direction = vec![0.0; self.dimension()];
                direction[axis] = 1.0;
                direction
            })
            .collect_vec();
        let total = points.iter().map(|point| point.weight).sum::<f64>();
        let mut uncovered = (0..points.len()).collect_vec();
        let mut centers: Vec<Vec<f64>> = vec![];
        while centers.len() < self.settings.max_patches && !uncovered.is_empty() {
            // Subsample each surface evenly, in a fixed order, to bound the size of the
            // linear program deterministically.
            let selected = uncovered
                .iter()
                .map(|&index| &points[index])
                .sorted_by_key(|point| point.surface)
                .chunk_by(|point| point.surface)
                .into_iter()
                .flat_map(|(_, group)| {
                    let group = group.collect_vec();
                    let stride = group.len().div_ceil(LP_POINTS_PER_SURFACE);
                    group.into_iter().step_by(stride)
                })
                .collect_vec();
            let basis = if centers.is_empty() {
                &invariant
            } else {
                &full
            };
            let Some(center) = self.cover(&selected, basis) else {
                break;
            };
            let center = self.refine(
                &uncovered.iter().map(|&index| &points[index]).collect_vec(),
                basis,
                center,
            );
            let (covered, remaining): (Vec<usize>, Vec<usize>) = uncovered
                .iter()
                .partition(|&&index| points[index].covered_by(&center));
            let gain = covered
                .iter()
                .map(|&index| points[index].weight)
                .sum::<f64>()
                / total;
            if !centers.is_empty() && gain < self.settings.min_patch_gain {
                break;
            }
            centers.push(center);
            uncovered = remaining;
        }
        centers
    }

    fn run(mut self) -> LuFlowPlan {
        let dimension = self.dimension();
        let mut bases = vec![vec![0.0; dimension]];
        for _ in 0..RANDOM_BASES {
            let length = self.gate.scale * self.rng.random::<f64>();
            bases.push(self.direction(length));
        }
        let first = self.sample(&bases);
        let reference = first
            .iter()
            .map(|point| point.log_weight)
            .fold(f64::MIN, f64::max);
        let centers = self.greedy(&self.constrain(first, reference));
        // A second round samples rays through the first centers, which expose the cut
        // points those centers see tangentially or with the wrong orientation.
        let mut second = self.sample(&bases);
        second.extend(self.sample(&centers));
        let refined = self.greedy(&self.constrain(second, reference));
        let centers = if refined.is_empty() { centers } else { refined };
        if centers.is_empty() {
            // No open cut point or no solvable program: the graph flows about the origin,
            // where its closed cuts contribute nothing.
            info!(
                "LU center search for graph '{}' found no relevant cut points; using the origin",
                self.graph
            );
            return LuFlowPlan {
                centers: vec![None],
                gate: self.gate,
            };
        }
        // Verify the soft partition on an independent sample, including rays through
        // every chosen center.
        bases.extend(centers.iter().cloned());
        let verification = self.sample(&bases);
        self.report(&centers, &self.constrain(verification, reference));
        LuFlowPlan {
            centers: centers
                .iter()
                .map(|center| {
                    center
                        .iter()
                        .any(|coordinate| *coordinate != 0.0)
                        .then(|| self.loops(center))
                })
                .collect(),
            gate: self.gate,
        }
    }

    /// Log the Boltzmann-weighted coverage of the verification sample. Points near a
    /// second cut surface are counted separately: where no patch is causal and the first
    /// center crosses the two surfaces with opposite slopes, their overlap pole is not
    /// cancelled, however small its Boltzmann weight.
    fn report(&self, centers: &[Vec<f64>], points: &[CutPoint]) {
        let total = points.iter().map(|point| point.weight).sum::<f64>();
        let (mut full, mut partial, mut negative, mut overlaps, mut opposed) =
            (0.0, 0.0, 0.0, 0, 0);
        for point in points {
            let validity = centers
                .iter()
                .map(|center| point.validity(&self.gate, center))
                .fold(0.0, f64::max);
            if validity >= 1.0 {
                full += point.weight;
            } else if validity > 0.0 {
                partial += point.weight;
            }
            let first = point.slopes(&centers[0]).collect_vec();
            if first[0] <= 0.0 {
                negative += point.weight;
            }
            if first.len() > 1 {
                overlaps += 1;
                if validity == 0.0 && first.iter().any(|slope| slope * first[0] < 0.0) {
                    opposed += 1;
                }
            }
        }
        info!(
            "LU center search for graph '{}': {} patch(es), centers {:?}; of {} relevant cut points, Boltzmann-weighted coverage {:.4} causal, {:.4} partial, {:.4} uncovered; first center misoriented on {:.4}; {} of {} overlap points keep opposite slopes",
            self.graph,
            centers.len(),
            centers
                .iter()
                .map(|center| center.iter().map(|x| format!("{x:.4}")).join(","))
                .collect_vec(),
            points.len(),
            full / total,
            partial / total,
            1.0 - (full + partial) / total,
            negative / total,
            opposed,
            overlaps,
        );
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};

    use super::*;
    use crate::{
        cff::{VertexSet, hsurface::Hsurface},
        momentum::{FourMomentum, RotationMethod, signature::LoopExtSignature},
        utils::{
            hyperdual_utils::{extract_t_derivatives, simple_n_deriv_shape},
            test_utils::dummy_hedge_graph,
        },
    };

    /// The scalar H/H overlap of the thermal vertex at spacelike `K=(0.5,0,0,1.5)`.
    /// Edges are `(K-r, r-q, q, K-q, r, K)` with masses `(H,H,L,L,H,X)`; both absorption
    /// surfaces vanish at the sample `(q, r)` below. Origin flow crosses them with opposite
    /// slopes, which leaves an uncancelled pole between the two cuts.
    struct Overlap {
        lmb: LoopMomentumBasis,
        externals: ExternalFourMomenta<F<f64>>,
        masses: EdgeVec<F<f64>>,
        hsurfaces: Vec<Hsurface>,
        absorbed: Vec<[EdgeIndex; 1]>,
        point: Vec<f64>,
    }

    impl Overlap {
        fn new() -> Self {
            let graph = dummy_hedge_graph(6);
            let signatures = [
                ([0_isize, -1], [1_isize]),
                ([-1, 1], [0]),
                ([1, 0], [0]),
                ([-1, 0], [1]),
                ([0, 1], [0]),
                ([0, 0], [1]),
            ]
            .into_iter()
            .map(|(internal, external)| {
                LoopExtSignature::from((internal.to_vec(), external.to_vec()))
            });
            let lmb = LoopMomentumBasis {
                tree: SuBitGraph::empty(0),
                loop_edges: [EdgeIndex(2), EdgeIndex(4)].into_iter().collect(),
                ext_edges: [EdgeIndex(5)].into_iter().collect(),
                edge_signatures: graph.new_edgevec_from_iter(signatures).unwrap(),
            };
            let pairs = [(2, 3), (0, 4)];
            let hsurfaces = pairs
                .map(|(positive, negative)| Hsurface {
                    positive_energies: vec![EdgeIndex(positive)],
                    negative_energies: vec![EdgeIndex(negative)],
                    external_shift: vec![(EdgeIndex(5), 1)],
                    vertex_set: VertexSet::dummy(),
                })
                .to_vec();
            Self {
                externals: [FourMomentum::from_args(F(0.5), F(0.0), F(0.0), F(1.5))]
                    .into_iter()
                    .collect(),
                masses: graph
                    .new_edgevec_from_iter([0.6, 0.6, 0.25, 0.25, 0.6, 0.0].into_iter().map(F))
                    .unwrap(),
                lmb,
                hsurfaces,
                absorbed: pairs.map(|(positive, _)| [EdgeIndex(positive)]).to_vec(),
                point: vec![0.4, 0.2, 0.4412727741192883, 0.3, -0.4, 1.1224916106437834],
            }
        }

        /// The represented equations carry `+E`; the physical ones are their negatives.
        fn surfaces(&self) -> Vec<FlowSurface<'_>> {
            self.hsurfaces
                .iter()
                .zip(&self.absorbed)
                .map(|(surface, absorbed)| FlowSurface {
                    surface: HybridSurfaceRef::Hsurface(surface),
                    orientation: -1,
                    absorbed,
                })
                .collect()
        }

        fn search<'a>(&'a self, settings: &'a LuCenterSearchSettings) -> CenterSearch<'a> {
            CenterSearch {
                graph: "overlap",
                lmb: &self.lmb,
                externals: &self.externals,
                masses: &self.masses,
                surfaces: self.surfaces(),
                beta: 1.0,
                gate: PatchGate {
                    scale: 1.5,
                    slope: settings.slope * 1.5,
                },
                settings,
                rng: SmallRng::seed_from_u64(settings.seed),
            }
        }

        /// Physical slopes `grad(eta).(y - c)` of both surfaces at the overlap.
        fn slopes(&self, search: &CenterSearch<'_>, center: &[f64]) -> Vec<f64> {
            self.surfaces()
                .iter()
                .map(|surface| {
                    let eta = search.value(surface, &self.point);
                    assert!(eta.abs() < 1e-12, "the sample lies on both surfaces");
                    search
                        .oriented_gradient(surface, &self.point)
                        .iter()
                        .zip(&self.point)
                        .zip(center)
                        .map(|((n, y), c)| n * (y - c))
                        .sum()
                })
                .collect()
        }
    }

    #[test]
    fn search_orients_the_scalar_overlap_that_origin_flow_cancels_incorrectly() {
        let overlap = Overlap::new();
        let settings = LuCenterSearchSettings {
            rays: 48,
            ..Default::default()
        };
        let search = overlap.search(&settings);
        let origin = overlap.slopes(&search, &[0.0; 6]);
        assert!(origin[0] * origin[1] < 0.0, "origin slopes {origin:?}");
        // Independent check of the routed gradients against the pilot's derivatives.
        assert!((origin[0] + 0.8109160583319578).abs() < 1e-12);
        assert!((origin[1] - 1.304509868856776).abs() < 1e-12);

        let plan = search.run();
        let center = plan.centers::<f64>(&Rotation::new(RotationMethod::Identity))[0]
            .clone()
            .unwrap();
        let center = center
            .iter()
            .flat_map(|p| [p.px.0, p.py.0, p.pz.0])
            .collect_vec();
        let search = overlap.search(&settings);
        let slopes = overlap.slopes(&search, &center);
        assert!(
            slopes.iter().all(|slope| *slope > 0.0),
            "searched center {center:?} slopes {slopes:?}"
        );
    }

    #[test]
    fn patch_weights_partition_unity_and_their_jets_follow_the_ray() {
        let overlap = Overlap::new();
        let settings = LuCenterSearchSettings::default();
        let search = overlap.search(&settings);
        let midpoint = vec![0.0, 0.0, 0.75, 0.0, 0.0, 0.75];
        let shifted = vec![0.05, 0.0, 0.75, 0.05, 0.0, 0.75];
        // Slopes of about 0.35 at the overlap sit inside the gate's smooth ramp.
        let plan = LuFlowPlan {
            centers: vec![
                None,
                Some(search.loops(&midpoint)),
                Some(search.loops(&shifted)),
            ],
            gate: PatchGate {
                scale: 1.5,
                slope: 1.1,
            },
        };
        let identity = Rotation::new(RotationMethod::Identity);
        let centers = plan.centers::<f64>(&identity);
        let point = search.loops(&overlap.point);
        let zero = LoopMomenta::from_iter(
            point
                .iter()
                .map(|_| ThreeMomentum::new(F(0.0), F(0.0), F(0.0))),
        );
        let rays = |patch: usize| {
            let center = centers[patch].clone().unwrap_or_else(|| zero.clone());
            overlap
                .hsurfaces
                .iter()
                .map(|surface| {
                    (
                        -1,
                        surface.routed_ray(
                            &(&point - &center),
                            &center,
                            &overlap.externals,
                            &overlap.masses,
                            &overlap.lmb,
                        ),
                    )
                })
                .collect_vec()
        };
        let mut total = 0.0;
        for patch in 0..3 {
            let rays = rays(patch);
            let geometry = plan.patch_geometry(
                patch,
                &centers,
                rays.iter().map(|(orientation, ray)| (*orientation, ray)),
                &overlap.lmb,
                &F(0.0),
            );
            let weight = plan.weight(&geometry, &F(1.0), &point).0;
            assert!(
                (0.0..=1.0).contains(&weight),
                "patch {patch} weight {weight}"
            );
            total += weight;
            if patch == 0 {
                assert!(weight < 1e-6, "origin is not causal at the overlap");
                continue;
            }
            assert!(weight > 0.05, "patch {patch} weight {weight}");

            // Jets along the patch ray at t=1 agree with central differences.
            let dual = HyperDual::<F<f64>>::new(simple_n_deriv_shape(2)).variable(0, F(1.0));
            let center = centers[patch].clone().unwrap();
            let dual_point = point.rescale_with_hyper_dual(&dual, None, Some(&center));
            let jets = extract_t_derivatives(plan.weight_dual(&geometry, &dual, &dual_point));
            let at = |t: f64| {
                plan.weight(&geometry, &F(t), &point.rescale(&F(t), None, Some(&center)))
                    .0
            };
            let step = 1e-4;
            assert!((jets[0].0 - weight).abs() < 1e-14);
            let first = (at(1.0 + step) - at(1.0 - step)) / (2.0 * step);
            let second = (at(1.0 + step) - 2.0 * weight + at(1.0 - step)) / (step * step);
            assert!(jets[1].0.abs() > 1e-3, "nontrivial slope {}", jets[1].0);
            assert!((jets[1].0 - first).abs() < 1e-6 * (1.0 + first.abs()));
            assert!((jets[2].0 - second).abs() < 1e-3 * (1.0 + second.abs()));
        }
        assert!((total - 1.0).abs() < 1e-14, "weights sum to {total}");
    }

    #[test]
    fn smooth_step_and_gate_limits() {
        assert_eq!(smooth_step(F(-0.5)).0, 0.0);
        assert_eq!(smooth_step(F(1.5)).0, 1.0);
        assert!((smooth_step(F(0.5)).0 - 0.5).abs() < 1e-15);
        // Far from the surface, or causal, the gate is open; close and acausal, closed.
        assert_eq!(gate(F(1.0), F(-10.0)).0, 1.0);
        assert_eq!(gate(F(0.0), F(1.0)).0, 1.0);
        assert_eq!(gate(F(0.2), F(0.25)).0, 0.0);
        assert_eq!(gate(F(0.0), F(-3.0)).0, 0.0);
    }
}
