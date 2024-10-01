use bincode::{Decode, Encode};
/// Counterterm for amplitudes with constant externals.
use colored::Colorize;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use symbolica::domains::float::{NumericalFloatLike, Real};

const MAX_ITERATIONS: usize = 20;
const TOLERANCE: f64 = 10.0;

use crate::graph::BareGraph;
use crate::momentum::{Rotatable, Rotation};

use crate::numerator::{Evaluators, Numerator};
use crate::{
    cff::{
        esurface::{
            compute_esurface_cache, EsurfaceCache, EsurfaceCollection, EsurfaceID,
            ExistingEsurfaceId, ExistingEsurfaces,
        },
        expression::CFFLimit,
    },
    debug_info::DEBUG_LOGGER,
    gammaloop_integrand::DefaultSample,
    graph::{Graph, LoopMomentumBasis},
    momentum::{FourMomentum, ThreeMomentum},
    numerator::NumeratorState,
    utils::{self, into_complex_ff64, FloatLike, F},
    Settings,
};

use super::overlap::OverlapStructure;

#[allow(clippy::type_complexity)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct CounterTerm {
    #[bincode(with_serde)]
    existing_esurfaces: ExistingEsurfaces,
    #[bincode(with_serde)]
    maximal_overlap: OverlapStructure,
    #[bincode(with_serde)]
    complements_of_overlap: Vec<Vec<ExistingEsurfaceId>>,
    #[bincode(with_serde)]
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

    pub fn construct<S: NumeratorState>(
        maximal_overlap: OverlapStructure,
        existing_esurfaces: &ExistingEsurfaces,
        graph: &Graph<S>,
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

        let (dep_mom, dep_mom_expr) = graph.bare_graph.get_dep_mom_expr();

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
                            .string_format_in_lmb(&graph.bare_graph.loop_momentum_basis),
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
        sample: &DefaultSample<T>,
        graph: &BareGraph,
        esurfaces: &EsurfaceCollection,
        counterterm: &CounterTerm,
        numerator: &mut Numerator<Evaluators>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let real_mass_vector = graph
            .get_real_mass_vector()
            .into_iter()
            .map(F::from_ff64)
            .collect_vec();

        let e_cm = F::from_ff64(settings.kinematics.e_cm);
        let const_builder = sample.zero();
        let mut res = Complex::new(const_builder.zero(), const_builder.zero());

        for (overlap_group, overlap_complement) in counterterm
            .maximal_overlap
            .overlap_groups
            .iter()
            .zip(counterterm.complements_of_overlap.iter())
        {
            let overlap = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            if settings.general.debug > 0 {
                DEBUG_LOGGER.write(
                    "overlap_structure",
                    &(
                        overlap
                            .iter()
                            .map(|id| counterterm.existing_esurfaces[*id])
                            .collect_vec(),
                        center,
                    ),
                );
            }

            let (unrotated_center, rotated_center): (Vec<ThreeMomentum<F<T>>>, _) = (
                center
                    .iter()
                    .map(|c_vec| c_vec.map(&|x| F::from_ff64(x)))
                    .collect_vec(),
                center
                    .iter()
                    .map(|c_vec| c_vec.rotate(rotation_for_overlap).map(&|x| F::from_ff64(x)))
                    .collect_vec(),
            );

            if settings.general.debug > 0 {
                DEBUG_LOGGER.write("center", &rotated_center);
            }

            let shifted_momenta = sample
                .loop_moms()
                .iter()
                .zip(rotated_center.iter())
                .map(|(momentum, center)| momentum.clone() - center.clone())
                .collect_vec();

            let (hemispherical_unit_shifted_momenta, hemispherical_radius) =
                normalize_momenta(&shifted_momenta);

            if settings.general.debug > 0 {
                DEBUG_LOGGER.write("shifted_momenta", &shifted_momenta);
                DEBUG_LOGGER.write("hemispherical_radius", &hemispherical_radius);
            }

            for existing_esurface_id in overlap.iter() {
                let esurface_id = counterterm.existing_esurfaces[*existing_esurface_id];
                let esurface = &esurfaces[esurface_id];

                // solve the radius
                let (radius_guess_plus, radius_guess_negative) = esurface.get_radius_guess(
                    &hemispherical_unit_shifted_momenta,
                    sample.external_moms(),
                    &graph.loop_momentum_basis,
                    &real_mass_vector,
                );

                let function = |r: &_| {
                    esurface.compute_self_and_r_derivative(
                        r,
                        &hemispherical_unit_shifted_momenta,
                        &rotated_center,
                        sample.external_moms(),
                        &graph.loop_momentum_basis,
                        &real_mass_vector,
                    )
                };

                let positive_result = newton_iteration_and_derivative(
                    &radius_guess_plus,
                    function,
                    &F::from_f64(TOLERANCE),
                    MAX_ITERATIONS,
                    &e_cm,
                );

                let negative_result = newton_iteration_and_derivative(
                    &radius_guess_negative,
                    function,
                    &F::from_f64(TOLERANCE),
                    MAX_ITERATIONS,
                    &e_cm,
                );

                assert!(positive_result.solution > const_builder.zero());
                assert!(negative_result.solution < const_builder.zero());

                if settings.general.debug > 1 {
                    println!("Positive solution: ");
                    positive_result.debug_print(&e_cm);

                    println!("Negative solution: ");
                    negative_result.debug_print(&e_cm);
                }

                let loop_momenta_at_rstar_plus = hemispherical_unit_shifted_momenta
                    .iter()
                    .zip(&rotated_center)
                    .map(|(k, center)| k * &positive_result.solution + center)
                    .collect_vec();

                let loop_momenta_at_rstar_minus = hemispherical_unit_shifted_momenta
                    .iter()
                    .zip(&rotated_center)
                    .map(|(k, center)| k * &negative_result.solution + center)
                    .collect_vec();

                // tedious gymnastics to get the sample right, hopefully correct.
                let rplus_sample = match sample.rotated_sample {
                    None => {
                        let mut sample_to_modify = sample.clone();
                        sample_to_modify.sample.loop_moms = loop_momenta_at_rstar_plus;
                        sample_to_modify
                    }
                    Some(_) => {
                        let mut sample_to_modify = sample.clone();
                        sample_to_modify.rotated_sample.as_mut().unwrap().loop_moms =
                            loop_momenta_at_rstar_plus;

                        let new_unrotated_momenta = sample
                            .sample
                            .loop_moms
                            .iter()
                            .zip(&unrotated_center)
                            .map(|(k, center)| {
                                (k.clone() - center.clone())
                                    * hemispherical_radius.inv()
                                    * &positive_result.solution
                                    + center
                            })
                            .collect_vec();

                        sample_to_modify.sample.loop_moms = new_unrotated_momenta;
                        sample_to_modify
                    }
                };

                let rminus_sample = match sample.rotated_sample {
                    None => {
                        let mut sample_to_modify = sample.clone();
                        sample_to_modify.sample.loop_moms = loop_momenta_at_rstar_minus;
                        sample_to_modify
                    }
                    Some(_) => {
                        let mut sample_to_modify = sample.clone();
                        sample_to_modify.rotated_sample.as_mut().unwrap().loop_moms =
                            loop_momenta_at_rstar_minus;

                        let new_unrotated_momenta = sample
                            .sample
                            .loop_moms
                            .iter()
                            .zip(&unrotated_center)
                            .map(|(k, center)| {
                                (k.clone() - center.clone())
                                    * hemispherical_radius.inv()
                                    * &negative_result.solution
                                    + center
                            })
                            .collect_vec();

                        sample_to_modify.sample.loop_moms = new_unrotated_momenta;
                        sample_to_modify
                    }
                };

                let ((r_plus_eval, r_plus_energy_cache), (r_minus_eval, r_minus_energy_cache)) = (
                    CounterTerm::radius_star_eval(
                        &rplus_sample,
                        graph,
                        esurfaces,
                        counterterm,
                        numerator,
                        overlap_complement,
                        *existing_esurface_id,
                        settings,
                    ),
                    CounterTerm::radius_star_eval(
                        &rminus_sample,
                        graph,
                        esurfaces,
                        counterterm,
                        numerator,
                        overlap_complement,
                        *existing_esurface_id,
                        settings,
                    ),
                );

                let loop_number = sample.loop_moms().len();
                let (jacobian_ratio_plus, jacobian_ratio_minus) = (
                    (&positive_result.solution / &hemispherical_radius)
                        .abs()
                        .pow(3 * loop_number as u64 - 1),
                    (&negative_result.solution / &hemispherical_radius)
                        .abs()
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

                let i = Complex::new(const_builder.zero(), const_builder.one());

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

                let ct_plus = &r_plus_eval
                    * ((&hemispherical_radius - &positive_result.solution)
                        * &positive_result.derivative_at_solution)
                        .inv()
                    * &uv_damper_plus
                    * &singularity_dampener_plus
                    * &jacobian_ratio_plus;

                let minus_half = -const_builder.one() / (const_builder.one() + const_builder.one());

                let integrated_ct_plus = &i * const_builder.PI() * &minus_half * r_plus_eval
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

                let ct_minus = &r_minus_eval
                    * ((&hemispherical_radius - &negative_result.solution)
                        * &negative_result.derivative_at_solution)
                        .inv()
                    * &uv_damper_minus
                    * &singularity_damper_minus
                    * &jacobian_ratio_minus;

                let integrated_ct_minus = &i * const_builder.PI() * r_minus_eval * &minus_half
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
                        esurface_id,
                        initial_radius: radius_guess_plus.into_ff64(),
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
                        r_plus_energy_cache: r_plus_energy_cache
                            .iter()
                            .map(|x| x.into_ff64())
                            .collect_vec(),
                        r_minus_energy_cache: r_minus_energy_cache
                            .iter()
                            .map(|x| x.into_ff64())
                            .collect_vec(),
                    };

                    DEBUG_LOGGER.write("esurface_subtraction", &debug_helper);
                }
            }
        }

        // match the complex prefactor off cff
        let loop_number = graph.loop_momentum_basis.basis.len();

        let prefactor =
            Complex::new(const_builder.zero(), -const_builder.one()).pow(loop_number as u64);

        res * prefactor
    }

    // evaluate radius independent part
    #[allow(clippy::too_many_arguments)]
    fn radius_star_eval<T: FloatLike>(
        rstar_sample: &DefaultSample<T>,
        graph: &BareGraph,
        esurfaces: &EsurfaceCollection,
        counterterm: &CounterTerm,
        numerator: &mut Numerator<Evaluators>,
        overlap_complement: &[ExistingEsurfaceId],
        existing_esurface_id: ExistingEsurfaceId,
        settings: &Settings,
    ) -> (Complex<F<T>>, Vec<F<T>>) {
        let energy_cache =
            graph.compute_onshell_energies(rstar_sample.loop_moms(), rstar_sample.external_moms());

        let esurface_cache = compute_esurface_cache(esurfaces, &energy_cache);
        let rstar_energy_product = graph
            .get_virtual_edges_iterator()
            .map(|(edge_id, _)| F::from_f64(2.0) * &energy_cache[edge_id])
            .fold(energy_cache[0].one(), |acc, e| acc * e);

        let multichanneling_denominator =
            counterterm.evaluate_multichanneling_denominator(&esurface_cache);

        let multichanneling_numerator_root = overlap_complement
            .iter()
            .fold(energy_cache[0].one(), |acc, id| {
                acc * &esurface_cache[counterterm.existing_esurfaces[*id]]
            });

        let multichanneling_factor = &multichanneling_numerator_root
            * &multichanneling_numerator_root
            / multichanneling_denominator;

        let terms = &counterterm.terms_in_counterterms[Into::<usize>::into(existing_esurface_id)];

        let eval_terms = terms.evaluate_from_esurface_cache(
            graph,
            numerator,
            rstar_sample,
            &esurface_cache,
            &energy_cache,
            settings,
        );

        (
            eval_terms * &multichanneling_factor / rstar_energy_product,
            energy_cache,
        )
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
            derivative_at_solution: self.derivative_at_solution.into_ff64(),
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
    r_plus_energy_cache: Vec<F<f64>>,
    r_minus_energy_cache: Vec<F<f64>>,
}
