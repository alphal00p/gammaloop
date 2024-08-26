/// Counterterm for amplitudes with constant externals.
use colored::Colorize;
use itertools::Itertools;
use ref_ops::RefNeg;
use serde::Serialize;
use spenso::complex::Complex;
use symbolica::domains::float::{NumericalFloatLike, Real};

const MAX_ITERATIONS: usize = 20;
const TOLERANCE: f64 = 10.0;

use crate::{
    cff::{
        esurface::{
            compute_esurface_cache, EsurfaceCache, EsurfaceCollection, EsurfaceID,
            ExistingEsurfaceId, ExistingEsurfaces,
        },
        expression::CFFLimit,
    },
    debug_info::DEBUG_LOGGER,
    graph::{Graph, LoopMomentumBasis},
    momentum::{FourMomentum, ThreeMomentum},
    utils::{self, into_complex_ff64, FloatLike, F},
    RotationMethod, Settings,
};

use super::overlap::OverlapStructure;

#[allow(clippy::type_complexity)]
#[derive(Debug, Clone)]
pub struct CounterTerm {
    existing_esurfaces: ExistingEsurfaces,
    maximal_overlap: OverlapStructure,
    complements_of_overlap: Vec<Vec<ExistingEsurfaceId>>,
    terms_in_counterterms: Vec<CFFLimit>,
}

impl CounterTerm {
    pub fn print_debug_data(
        &self,
        esurfaces: &EsurfaceCollection,
        external_momenta: &[FourMomentum<F<f64>>],
        lmb: &LoopMomentumBasis,
        real_mass_vector: &[F<f64>],
    ) {
        let number_of_existing_esurfaces = self.existing_esurfaces.len();
        let simplified_maximal_overlap_structure = self
            .maximal_overlap
            .overlap_groups
            .iter()
            .map(|groups| groups.existing_esurfaces.len())
            .collect_vec();

        println!("number of thresholds: {}", number_of_existing_esurfaces);
        println!(
            "overlap structure: {:?}",
            simplified_maximal_overlap_structure
        );

        for overlap_group in self.maximal_overlap.overlap_groups.iter() {
            let overlap = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            println!("overlap size: {}", overlap.len());
            println!("center: {:#?}", center);

            println!("values of esurface at the center:");

            for esurface_id in overlap.iter() {
                let esurface = &esurfaces[self.existing_esurfaces[*esurface_id]];
                println!("esurface equation: {}", esurface.string_format_in_lmb(lmb));
                println!(
                    "value: {:?}",
                    esurface.compute_from_momenta(lmb, real_mass_vector, center, external_momenta)
                );
            }
        }
    }

    pub fn construct(
        maximal_overlap: OverlapStructure,
        existing_esurfaces: &ExistingEsurfaces,
        graph: &Graph,
    ) -> Self {
        let cff = graph.get_cff();
        let complements_of_overlap = maximal_overlap
            .overlap_groups
            .iter()
            .map(|overlap_group| {
                let esurface_ids = &overlap_group.existing_esurfaces;
                existing_esurfaces
                    .iter_enumerated()
                    .map(|a| a.0)
                    .filter(|id| !esurface_ids.contains(id))
                    .collect_vec()
            })
            .collect_vec();

        let (dep_mom, dep_mom_expr) = graph.get_dep_mom_expr();

        let terms_in_counterterms = existing_esurfaces
            .iter()
            .map(|existing_esurface| {
                let limit_result =
                    cff.limit_for_esurface(*existing_esurface, dep_mom, &dep_mom_expr);
                match limit_result {
                    Ok(limit) => limit,
                    Err(error_message) => panic!(
                        "Could not generate counterterm for esurface: {:?}\n
                            esurface: {:?}\n
                            esurface_in_lmb: {}\n
                            error_message: {},
                            ",
                        existing_esurface,
                        cff.esurfaces[*existing_esurface],
                        cff.esurfaces[*existing_esurface]
                            .string_format_in_lmb(&graph.loop_momentum_basis),
                        error_message
                    ),
                }
            })
            .collect_vec();

        Self {
            existing_esurfaces: existing_esurfaces.clone(),
            maximal_overlap,
            complements_of_overlap,
            terms_in_counterterms,
        }
    }

    // the multichanneling denominator will be evaluated at r*, so we can not use the cache.
    fn evaluate_multichanneling_denominator<T: FloatLike>(
        &self,
        esurface_cache: &EsurfaceCache<F<T>>,
    ) -> F<T> {
        let const_builder = &esurface_cache[EsurfaceID::from(0usize)];

        self.complements_of_overlap
            .iter()
            .map(|esurface_id_group| {
                let group = &esurface_id_group
                    .iter()
                    .map(|esurface_id| &esurface_cache[self.existing_esurfaces[*esurface_id]])
                    .fold(const_builder.one(), |acc, e| acc * e); // Dreams of product::<T>...

                group * group
            })
            .reduce(|acc, x| &acc + &x)
            .unwrap_or_else(|| const_builder.zero())
    }

    pub fn evaluate<T: FloatLike>(
        &self,
        loop_momenta: &[ThreeMomentum<F<T>>],
        external_momenta: &[FourMomentum<F<T>>],
        graph: &Graph,
        rotation_for_overlap: RotationMethod,
        settings: &Settings,
    ) -> Complex<F<T>> {
        if settings.general.debug > 1 {
            println!("{}", "Start evaluation of threshold counterterms".green());
        }

        let real_mass_vector = graph
            .get_real_mass_vector()
            .into_iter()
            .map(F::from_ff64)
            .collect_vec();

        let e_cm = F::from_ff64(settings.kinematics.e_cm);
        let esurfaces = &graph.get_cff().esurfaces;
        let lmb = &graph.loop_momentum_basis;

        let const_builder = &loop_momenta[0].px.zero();

        let mut res = Complex::new(const_builder.zero(), const_builder.zero());

        for (overlap_group, overlap_complement) in self
            .maximal_overlap
            .overlap_groups
            .iter()
            .zip(self.complements_of_overlap.iter())
        {
            let overlap = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            if settings.general.debug > 0 {
                DEBUG_LOGGER.write(
                    "overlap_structure",
                    &(
                        overlap
                            .iter()
                            .map(|id| self.existing_esurfaces[*id])
                            .collect_vec(),
                        center,
                    ),
                );
            }

            let center_t = {
                let rotation_function = rotation_for_overlap.rotation_function();
                let rotated_center = center
                    .iter()
                    .map(rotation_function)
                    .map(|c_vec| {
                        ThreeMomentum::new(
                            F::from_ff64(c_vec.px),
                            F::from_ff64(c_vec.py),
                            F::from_ff64(c_vec.pz), // manual cast is needed for some reason?
                        )
                    })
                    .collect_vec();

                if settings.general.debug > 0 {
                    DEBUG_LOGGER.write("center", &rotated_center);
                }

                rotated_center
            };

            let shifted_momenta = loop_momenta
                .iter()
                .zip(center_t.iter())
                .map(|(momentum, center)| momentum.clone() - center.clone())
                .collect_vec();

            let (hemispherical_unit_shifted_momenta, hemispherical_radius) =
                normalize_momenta(&shifted_momenta);

            if settings.general.debug > 0 {
                DEBUG_LOGGER.write("shifted_momenta", &shifted_momenta);
                DEBUG_LOGGER.write("hemispherical_radius", &hemispherical_radius);
            }

            for existing_esurface_id in overlap.iter() {
                let esurface_id = &self.existing_esurfaces[*existing_esurface_id];
                let esurface = &esurfaces[*esurface_id];

                if settings.general.debug > 1 {
                    println!("subtracting esurface: {:?}", esurface_id);
                    println!("energies: {:?}", &esurface.energies);
                    println!("external shift: {:?}", &esurface.external_shift);
                }

                // solve the radius
                let radius_guess = esurface.get_radius_guess(
                    &hemispherical_unit_shifted_momenta,
                    external_momenta,
                    lmb,
                );

                if settings.general.debug > 1 {
                    println!("Initial condition for newton iterations: ±{}", radius_guess);
                }

                let function = |r: &_| {
                    esurface.compute_self_and_r_derivative(
                        r,
                        &hemispherical_unit_shifted_momenta,
                        &center_t,
                        external_momenta,
                        lmb,
                        &real_mass_vector,
                    )
                };

                let positive_result = newton_iteration_and_derivative(
                    &radius_guess,
                    function,
                    &F::from_f64(TOLERANCE),
                    MAX_ITERATIONS,
                    &e_cm,
                );

                let negative_result = newton_iteration_and_derivative(
                    &radius_guess.ref_neg(),
                    function,
                    &F::from_f64(TOLERANCE),
                    MAX_ITERATIONS,
                    &e_cm,
                );

                if settings.general.debug > 1 {
                    println!("Positive solution: ");
                    positive_result.debug_print(&e_cm);

                    println!("Negative solution: ");
                    negative_result.debug_print(&e_cm);
                }

                let (r_plus_eval, r_minus_eval) = (
                    self.radius_star_eval(
                        &positive_result.solution,
                        &hemispherical_unit_shifted_momenta,
                        &center_t,
                        graph,
                        external_momenta,
                        esurfaces,
                        overlap_complement,
                        *existing_esurface_id,
                    ),
                    self.radius_star_eval(
                        &negative_result.solution,
                        &hemispherical_unit_shifted_momenta,
                        &center_t,
                        graph,
                        external_momenta,
                        esurfaces,
                        overlap_complement,
                        *existing_esurface_id,
                    ),
                );

                let loop_number = loop_momenta.len();
                let (jacobian_ratio_plus, jacobian_ratio_minus) = (
                    (&positive_result.solution / &hemispherical_radius)
                        .pow(3 * loop_number as u64 - 1),
                    (&negative_result.solution / &hemispherical_radius)
                        .pow(3 * loop_number as u64 - 1),
                );

                let (uv_damper_plus, uv_damper_minus) = (
                    unnormalized_gaussian(&hemispherical_radius, &positive_result.solution, &e_cm),
                    unnormalized_gaussian(&hemispherical_radius, &negative_result.solution, &e_cm),
                );

                let (singularity_dampener_plus, singularity_damper_minus) =
                    if settings.subtraction.dampen_integrable_singularity {
                        (
                            singularity_dampener(&hemispherical_radius, &positive_result.solution),
                            singularity_dampener(&hemispherical_radius, &negative_result.solution),
                        )
                    } else {
                        (hemispherical_radius.one(), hemispherical_radius.one())
                    };

                let i = Complex::new(radius_guess.zero(), radius_guess.one());

                let radius_sign_plus = if positive_result.solution > hemispherical_radius.zero() {
                    const_builder.one()
                } else {
                    -const_builder.one()
                };

                let radius_sign_minus = if positive_result.solution > hemispherical_radius.zero() {
                    const_builder.one()
                } else {
                    -const_builder.one()
                };

                let ct_plus = Complex::new(
                    &r_plus_eval
                        * ((&hemispherical_radius - &positive_result.solution)
                            * &positive_result.derivative_at_solution)
                            .inv()
                        * &uv_damper_plus
                        * &singularity_dampener_plus
                        * &jacobian_ratio_plus,
                    radius_guess.zero(),
                );

                let minus_half = -const_builder.one() / (const_builder.one() + const_builder.one());

                let integrated_ct_plus = &i * radius_guess.PI() * &minus_half * r_plus_eval
                    / &positive_result.derivative_at_solution
                    / &positive_result.solution
                    * utils::h(
                        &positive_result.solution / &hemispherical_radius,
                        None,
                        None,
                        &settings.subtraction.integrated_ct_hfunction,
                    )
                    * &jacobian_ratio_plus
                    * &radius_sign_plus;

                let ct_minus = Complex::new(
                    &r_minus_eval
                        * ((&hemispherical_radius - &negative_result.solution)
                            * &negative_result.derivative_at_solution)
                            .inv()
                        * &uv_damper_minus
                        * &singularity_damper_minus
                        * &jacobian_ratio_minus,
                    radius_guess.zero(),
                );

                let integrated_ct_minus = &i * radius_guess.PI() * r_minus_eval * &minus_half
                    / &negative_result.derivative_at_solution
                    / &negative_result.solution
                    * utils::h(
                        &negative_result.solution / &hemispherical_radius,
                        None,
                        None,
                        &settings.subtraction.integrated_ct_hfunction,
                    )
                    * &jacobian_ratio_minus
                    * &radius_sign_minus;

                res += &ct_plus + &integrated_ct_plus + &ct_minus + &integrated_ct_minus;

                if settings.general.debug > 0 {
                    let debug_helper = DebugHelper {
                        esurface_id: *esurface_id,
                        initial_radius: radius_guess.into_ff64(),
                        plus_solution: positive_result.as_f64(),
                        minus_solution: negative_result.as_f64(),
                        jacobian_ratio_plus: jacobian_ratio_plus.into_ff64(),
                        jacobian_ratio_minus: jacobian_ratio_minus.into_ff64(),
                        uv_damper_plus: uv_damper_plus.into_ff64(),
                        uv_damper_minus: uv_damper_minus.into_ff64(),
                        singularity_dampener_plus: singularity_dampener_plus.into_ff64(),
                        singularity_dampener_minus: singularity_damper_minus.into_ff64(),
                        ct_plus: into_complex_ff64(&ct_plus),
                        ct_minus: into_complex_ff64(&ct_minus),
                        integrated_ct_plus: into_complex_ff64(&integrated_ct_plus),
                        integrated_ct_minus: into_complex_ff64(&integrated_ct_minus),
                    };

                    DEBUG_LOGGER.write("esurface_subtraction", &debug_helper);
                }
            }
        }

        // match the complex prefactor off cff
        let loop_number = graph.loop_momentum_basis.basis.len();
        let internal_vertex_number = graph.vertices.len() - graph.external_connections.len();

        let prefactor = Complex::new(const_builder.zero(), const_builder.one())
            .pow(loop_number as u64)
            * Complex::new(-const_builder.one(), const_builder.zero())
                .pow(internal_vertex_number as u64 - 1);

        res * prefactor
    }

    // evaluate radius independent part
    #[allow(clippy::too_many_arguments)]
    fn radius_star_eval<T: FloatLike>(
        &self,
        rstar: &F<T>,
        unit_loop_momenta: &[ThreeMomentum<F<T>>],
        center: &[ThreeMomentum<F<T>>],
        graph: &Graph,
        external_momenta: &[FourMomentum<F<T>>],
        esurfaces: &EsurfaceCollection,
        overlap_complement: &[ExistingEsurfaceId],
        existing_esurface_id: ExistingEsurfaceId,
    ) -> F<T> {
        let loop_momenta_at_star = unit_loop_momenta
            .iter()
            .zip(center.iter())
            .map(|(k, center)| k * rstar + center)
            .collect_vec();

        let energy_cache = graph.compute_onshell_energies(&loop_momenta_at_star, external_momenta);
        let esurface_cache = compute_esurface_cache(esurfaces, &energy_cache);
        let rstar_energy_product = graph
            .get_virtual_edges_iterator()
            .map(|(edge_id, _)| F::from_f64(2.0) * &energy_cache[edge_id])
            .fold(rstar.one(), |acc, e| acc * e);

        let multichanneling_denominator =
            self.evaluate_multichanneling_denominator(&esurface_cache);

        let multichanneling_numerator_root =
            overlap_complement.iter().fold(rstar.one(), |acc, id| {
                acc * &esurface_cache[self.existing_esurfaces[*id]]
            });

        let multichanneling_factor = &multichanneling_numerator_root
            * &multichanneling_numerator_root
            / multichanneling_denominator;

        let terms = &self.terms_in_counterterms[Into::<usize>::into(existing_esurface_id)];

        let eval_terms = terms.evaluate_from_esurface_cache(&esurface_cache, &energy_cache);

        multichanneling_factor * eval_terms / rstar_energy_product
    }
}

fn normalize_momenta<T: FloatLike>(
    momenta: &[ThreeMomentum<F<T>>],
) -> (Vec<ThreeMomentum<F<T>>>, F<T>) {
    // let the x plane of the first momenta carry all the sine functions in the definition of spherical coordinates
    // let sign_radius = momenta[0].px.signum();

    let first_momentum_x_component = &momenta[0].px;

    let sign_radius = if momenta[0].px.positive() {
        first_momentum_x_component.one()
    } else {
        first_momentum_x_component.one().neg()
    };

    let radius = sign_radius
        * momenta
            .iter()
            .map(|k| k.norm_squared())
            .reduce(|acc, x| &acc + &x)
            .unwrap()
            .sqrt();

    let normalized_momentum = momenta.iter().map(|k| k * radius.inv()).collect_vec();

    (normalized_momentum, radius)
}

/// root finding, returns the derivative at the root, so that we don't have to recompute it.
/// Also returns the value of the function whose root is being found and the number of iterations used for debug information
fn newton_iteration_and_derivative<T: FloatLike>(
    guess: &F<T>,
    f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
    tolerance: &F<T>,
    max_iterations: usize,
    e_cm: &F<T>,
) -> NewtonIterationResult<T> {
    let mut x = guess.clone();
    let (mut val_f_x, mut val_df_x) = f_x_and_df_x(&x);

    let mut iteration = 0;

    while iteration < max_iterations && val_f_x.abs() > guess.epsilon() * tolerance * e_cm {
        x -= val_f_x / val_df_x;
        (val_f_x, val_df_x) = f_x_and_df_x(&x);

        iteration += 1;
    }

    NewtonIterationResult {
        solution: x,
        derivative_at_solution: val_df_x,
        error_of_function: val_f_x,
        num_iterations_used: iteration,
    }
}

#[derive(Serialize)]
struct NewtonIterationResult<T: FloatLike> {
    solution: F<T>,
    derivative_at_solution: F<T>,
    error_of_function: F<T>,
    num_iterations_used: usize,
}

impl<T: FloatLike> NewtonIterationResult<T> {
    fn debug_print(&self, e_cm: &F<T>) {
        println!("r*: {}", &self.solution);
        println!("Derivative at r* {}", &self.derivative_at_solution);
        println!(
            "Value of esurface: {}",
            if self.error_of_function.abs()
                < self.solution.epsilon() * &F::from_f64(TOLERANCE) * e_cm
            {
                format!("{}", &self.error_of_function).green()
            } else {
                format!("{}", &self.error_of_function).red()
            }
        );
        println!(
            "Number of iterations used: {}",
            if self.num_iterations_used == MAX_ITERATIONS {
                format!("{}", &self.num_iterations_used).red()
            } else {
                format!("{}", self.num_iterations_used).green()
            }
        )
    }

    fn as_f64(&self) -> NewtonIterationResult<f64> {
        NewtonIterationResult {
            solution: self.solution.into_ff64(),
            derivative_at_solution: self.solution.into_ff64(),
            error_of_function: self.error_of_function.into_ff64(),
            num_iterations_used: self.num_iterations_used,
        }
    }
}

fn unnormalized_gaussian<T: FloatLike>(radius: &F<T>, radius_star: &F<T>, e_cm: &F<T>) -> F<T> {
    let delta_r = radius - radius_star;
    (-&delta_r * &delta_r / (e_cm * e_cm)).exp()
}

fn singularity_dampener<T: FloatLike>(radius: &F<T>, radius_star: &F<T>) -> F<T> {
    let delta_r = radius - radius_star;
    let radius_star_sq = radius_star * radius_star;
    let radius_star_4 = &radius_star_sq * &radius_star_sq;
    let delta_r_sq = &delta_r * &delta_r;

    (radius.one() - radius_star_4 / (radius_star_sq - delta_r_sq).pow(2)).exp()
}

/// Helper struct to group together all debug data related to the subtraction of a single esurface
#[derive(Serialize)]
struct DebugHelper {
    esurface_id: EsurfaceID,
    initial_radius: F<f64>,
    plus_solution: NewtonIterationResult<f64>,
    minus_solution: NewtonIterationResult<f64>,
    jacobian_ratio_plus: F<f64>,
    jacobian_ratio_minus: F<f64>,
    uv_damper_plus: F<f64>,
    uv_damper_minus: F<f64>,
    singularity_dampener_plus: F<f64>,
    singularity_dampener_minus: F<f64>,
    ct_plus: Complex<F<f64>>,
    ct_minus: Complex<F<f64>>,
    integrated_ct_plus: Complex<F<f64>>,
    integrated_ct_minus: Complex<F<f64>>,
}
