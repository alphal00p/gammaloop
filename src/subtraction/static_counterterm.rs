use bincode::{Decode, Encode};
/// Counterterm for amplitudes with constant externals.
use colored::Colorize;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use symbolica::domains::float::{NumericalFloatLike, Real};

const MAX_ITERATIONS: usize = 20;
const TOLERANCE: f64 = 5.0;

use crate::cff::esurface::Esurface;
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
    utils::{self, FloatLike, F},
    Settings,
};

use super::overlap::{OverlapGroup, OverlapStructure};

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
        let counterterm_builder = CounterTermBuilder::new(
            graph,
            esurfaces,
            counterterm,
            rotation_for_overlap,
            settings,
            sample,
        );

        let const_builder = sample.zero();

        let res = counterterm
            .maximal_overlap
            .overlap_groups
            .iter()
            .zip(counterterm.complements_of_overlap.iter())
            .map(|(overlap_group, overlap_complement)| {
                let overlap_builder =
                    counterterm_builder.new_overlap_builder(overlap_group, overlap_complement);

                overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|existing_esurface_id| {
                        let res = overlap_builder
                            .new_esurface_ct_builder(*existing_esurface_id)
                            .solve_rstar()
                            .new_rstar_samples()
                            .compute_residues(numerator)
                            .compute_counterterm();

                        if settings.general.debug > 0 {
                            //let ose_product = into_complex_ff64(&Complex::new(
                            //    graph.compute_energy_product(sample.loop_moms(), sample.external_moms()),
                            //    res.re.zero(),
                            //));
                            //let debug_helper = DebugHelper {
                            //    esurface_id,
                            //    edges: esurface.energies.clone(),
                            //    plus_solution: rstar_solution.solutions[0].as_f64(),
                            //    minus_solution: rstar_solution.solutions[1].as_f64(),
                            //    jacobian_ratio_plus: jacobian_ratio_plus.into_ff64(),
                            //    jacobian_ratio_minus: jacobian_ratio_minus.into_ff64(),
                            //    uv_damper_plus: uv_damper_plus.into_ff64(),
                            //    uv_damper_minus: uv_damper_minus.into_ff64(),
                            //    singularity_dampener_plus: singularity_dampener_plus.into_ff64(),
                            //    singularity_dampener_minus: singularity_damper_minus.into_ff64(),
                            //    ct_plus: into_complex_ff64(&ct_plus) * ose_product,
                            //    ct_minus: into_complex_ff64(&ct_minus) * ose_product,
                            //    integrated_ct_plus: into_complex_ff64(&integrated_ct_plus) * ose_product,
                            //    integrated_ct_minus: into_complex_ff64(&integrated_ct_minus) * ose_product,
                            //    r_plus_energy_cache: r_plus_energy_cache
                            //        .iter()
                            //        .map(|x| x.into_ff64())
                            //        .collect_vec(),
                            //    r_minus_energy_cache: r_minus_energy_cache
                            //        .iter()
                            //        .map(|x| x.into_ff64())
                            //        .collect_vec(),
                            //    r_plus_energy_product: r_plus_ose_product.into_ff64(),
                            //    r_minus_energy_product: r_minus_ose_product.into_ff64(),
                            //    r_plus_eval: into_complex_ff64(&r_plus_eval)
                            //        * r_plus_ose_product.into_ff64(),
                            //    r_minus_eval: into_complex_ff64(&r_minus_eval)
                            //        * r_minus_ose_product.into_ff64(),
                            //    h_plus: h_plus.into_ff64(),
                            //    h_minus: h_minus.into_ff64(),
                            //    r_plus_esurface_cache: r_plus_esurface_cache
                            //        .iter()
                            //        .enumerate()
                            //        .map(|(id, x)| (id, x.into_ff64()))
                            //        .collect_vec(),
                            //    r_minus_esurface_cache: r_minus_esurface_cache
                            //        .iter()
                            //        .enumerate()
                            //        .map(|(id, x)| (id, x.into_ff64()))
                            //        .collect_vec(),
                            //};
                        }
                        res
                    })
                    .fold(Complex::new_re(sample.zero()), |acc, x| acc + x)
            })
            .fold(Complex::new_re(sample.zero()), |acc, x| acc + x);

        // match the complex prefactor off cff
        let loop_number = graph.loop_momentum_basis.basis.len();

        let prefactor =
            Complex::new(const_builder.zero(), -const_builder.one()).pow(loop_number as u64);

        res * prefactor
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
    fn _debug_print(&self, e_cm: &F<T>) {
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
}

fn unnormalized_gaussian<T: FloatLike>(
    radius: &F<T>,
    radius_star: &F<T>,
    e_cm: &F<T>,
    width: &F<T>,
) -> F<T> {
    let delta_r = radius - radius_star;
    let sigma = e_cm * width;
    (-&delta_r * &delta_r / (&sigma * &sigma)).exp()
}

fn singularity_dampener<T: FloatLike>(radius: &F<T>, radius_star: &F<T>) -> F<T> {
    let delta_r = radius - radius_star;
    let radius_star_sq = radius_star * radius_star;
    let radius_star_4 = &radius_star_sq * &radius_star_sq;
    let delta_r_sq = &delta_r * &delta_r;

    (radius.one() - radius_star_4 / (radius_star_sq - delta_r_sq).pow(2)).exp()

    //(radius.one() - (radius / radius_star - radius.one()).abs())
    //    .abs()
    //    .pow(5)
}

struct CounterTermBuilder<'a, T: FloatLike> {
    real_mass_vector: Vec<F<T>>,
    e_cm: F<T>,
    graph: &'a BareGraph,
    counterterm: &'a CounterTerm,
    rotation_for_overlap: &'a Rotation,
    settings: &'a Settings,
    esurface_collection: &'a EsurfaceCollection,
    sample: &'a DefaultSample<T>,
}

impl<'a, T: FloatLike> CounterTermBuilder<'a, T> {
    fn new(
        graph: &'a BareGraph,
        esurface_collection: &'a EsurfaceCollection,
        counterterm: &'a CounterTerm,
        rotation_for_overlap: &'a Rotation,
        settings: &'a Settings,
        sample: &'a DefaultSample<T>,
    ) -> Self {
        let real_mass_vector = graph
            .get_real_mass_vector()
            .into_iter()
            .map(F::from_ff64)
            .collect_vec();

        let e_cm = F::from_ff64(settings.kinematics.e_cm);

        Self {
            real_mass_vector,
            e_cm,
            graph,
            esurface_collection,
            counterterm,
            rotation_for_overlap,
            settings,
            sample,
        }
    }

    fn new_overlap_builder(
        &'a self,
        overlap_group: &'a OverlapGroup,
        overlap_complement: &'a [ExistingEsurfaceId],
    ) -> OverlapBuilder<'a, T> {
        let overlap = &overlap_group.existing_esurfaces;
        let center = &overlap_group.center;

        if self.settings.general.debug > 0 {
            DEBUG_LOGGER.write(
                "overlap_structure",
                &(
                    overlap
                        .iter()
                        .map(|id| self.counterterm.existing_esurfaces[*id])
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
                .map(|c_vec| {
                    c_vec
                        .rotate(self.rotation_for_overlap)
                        .map(&|x| F::from_ff64(x))
                })
                .collect_vec(),
        );

        if self.settings.general.debug > 0 {
            DEBUG_LOGGER.write("center", &rotated_center);
        }

        let shifted_momenta = self
            .sample
            .loop_moms()
            .iter()
            .zip(rotated_center.iter())
            .map(|(momentum, center)| momentum.clone() - center.clone())
            .collect_vec();

        let (hemispherical_unit_shifted_momenta, hemispherical_radius) =
            normalize_momenta(&shifted_momenta);

        OverlapBuilder {
            counterterm_builder: self,
            overlap_complement,
            rotated_center,
            unrotated_center,
            hemispherical_unit_shifted_momenta,
            hemispherical_radius,
        }
    }
}

struct OverlapBuilder<'a, T: FloatLike> {
    counterterm_builder: &'a CounterTermBuilder<'a, T>,
    overlap_complement: &'a [ExistingEsurfaceId],
    rotated_center: Vec<ThreeMomentum<F<T>>>,
    unrotated_center: Vec<ThreeMomentum<F<T>>>,
    hemispherical_unit_shifted_momenta: Vec<ThreeMomentum<F<T>>>,
    hemispherical_radius: F<T>,
}

impl<'a, T: FloatLike> OverlapBuilder<'a, T> {
    fn new_esurface_ct_builder(
        &'a self,
        existing_esurface_id: ExistingEsurfaceId,
    ) -> EsurfaceCTBuilder<'_, T> {
        let esurface_id =
            self.counterterm_builder.counterterm.existing_esurfaces[existing_esurface_id];
        let esurface = &self.counterterm_builder.esurface_collection[esurface_id];

        if self.counterterm_builder.settings.general.debug > 0 {
            DEBUG_LOGGER.write("esurface_id: ", &Into::<usize>::into(esurface_id));
            DEBUG_LOGGER.write("energies: ", &esurface.energies);
        }

        EsurfaceCTBuilder {
            overlap_builder: self,
            existing_esurface_id,
            esurface,
        }
    }
}

struct EsurfaceCTBuilder<'a, T: FloatLike> {
    overlap_builder: &'a OverlapBuilder<'a, T>,
    existing_esurface_id: ExistingEsurfaceId,
    esurface: &'a Esurface,
}

impl<'a, T: FloatLike> EsurfaceCTBuilder<'a, T> {
    fn solve_rstar(self) -> RstarSolution<'a, T> {
        let settings = self.overlap_builder.counterterm_builder.settings;

        let (radius_guess_plus, radius_guess_negative) = self.esurface.get_radius_guess(
            &self.overlap_builder.hemispherical_unit_shifted_momenta,
            self.overlap_builder
                .counterterm_builder
                .sample
                .external_moms(),
            &self
                .overlap_builder
                .counterterm_builder
                .graph
                .loop_momentum_basis,
            &self.overlap_builder.counterterm_builder.real_mass_vector,
        );

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("radius_guess", &radius_guess_plus);
        }

        let function = |r: &_| {
            self.esurface.compute_self_and_r_derivative(
                r,
                &self.overlap_builder.hemispherical_unit_shifted_momenta,
                &self.overlap_builder.rotated_center,
                self.overlap_builder
                    .counterterm_builder
                    .sample
                    .external_moms(),
                &self
                    .overlap_builder
                    .counterterm_builder
                    .graph
                    .loop_momentum_basis,
                &self.overlap_builder.counterterm_builder.real_mass_vector,
            )
        };

        let solutions = [radius_guess_plus, radius_guess_negative].map(|guess| {
            newton_iteration_and_derivative(
                &guess,
                function,
                &F::from_f64(TOLERANCE),
                MAX_ITERATIONS,
                &self.overlap_builder.counterterm_builder.e_cm,
            )
        });

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("rstar_solutions", &solutions);
        }

        assert!(
            solutions[0].solution > solutions[0].solution.zero(),
            "positive result has wrong sign: ({}, {})",
            solutions[0].solution,
            solutions[1].solution
        );
        assert!(
            solutions[1].solution < solutions[1].solution.zero(),
            "negative result has wrong sign: ({}, {})",
            solutions[0].solution,
            solutions[1].solution
        );

        RstarSolution {
            solutions,
            esurface_ct_builder: self,
        }
    }
}

struct RstarSolution<'a, T: FloatLike> {
    solutions: [NewtonIterationResult<T>; 2],
    esurface_ct_builder: EsurfaceCTBuilder<'a, T>,
}

impl<'a, T: FloatLike> RstarSolution<'a, T> {
    fn new_rstar_samples(self) -> RstarSamples<'a, T> {
        let rstar_samples = self.solutions.each_ref().map(|solution| {
            let rstar_loop_momenta = self
                .esurface_ct_builder
                .overlap_builder
                .hemispherical_unit_shifted_momenta
                .iter()
                .zip(&self.esurface_ct_builder.overlap_builder.rotated_center)
                .map(|(k, center)| k * &solution.solution + center)
                .collect_vec();

            // some gymnastics to get the sample right
            // do not touch!
            let mut sample_to_modify = self
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .sample
                .clone();

            match self
                .esurface_ct_builder
                .overlap_builder
                .counterterm_builder
                .sample
                .rotated_sample
            {
                None => {
                    sample_to_modify.sample.loop_moms = rstar_loop_momenta;
                    sample_to_modify
                }
                Some(_) => {
                    sample_to_modify.rotated_sample.as_mut().unwrap().loop_moms =
                        rstar_loop_momenta;

                    let new_unrotated_momenta = self
                        .esurface_ct_builder
                        .overlap_builder
                        .counterterm_builder
                        .sample
                        .loop_moms()
                        .iter()
                        .zip(&self.esurface_ct_builder.overlap_builder.unrotated_center)
                        .map(|(k, center)| {
                            (k.clone() - center.clone())
                                * self
                                    .esurface_ct_builder
                                    .overlap_builder
                                    .hemispherical_radius
                                    .inv()
                                * &solution.solution
                                + center
                        })
                        .collect_vec();

                    sample_to_modify.sample.loop_moms = new_unrotated_momenta;
                    sample_to_modify
                }
            }
        });

        RstarSamples {
            rstar_samples,
            rstar_solutions: self,
        }
    }
}

struct RstarSamples<'a, T: FloatLike> {
    rstar_samples: [DefaultSample<T>; 2],
    rstar_solutions: RstarSolution<'a, T>,
}

impl<'a, T: FloatLike> RstarSamples<'a, T> {
    fn compute_residues(self, numerator: &mut Numerator<Evaluators>) -> ResiudeEval<'a, T> {
        let esurface_ct_builder = &self.rstar_solutions.esurface_ct_builder;
        let overlap_builder = esurface_ct_builder.overlap_builder;
        let ct_builder = overlap_builder.counterterm_builder;

        let settings = ct_builder.settings;
        let graph = ct_builder.graph;
        let esurfaces = ct_builder.esurface_collection;
        let counterterm = ct_builder.counterterm;
        let overlap_complement = overlap_builder.overlap_complement;
        let existing_esurface_id = esurface_ct_builder.existing_esurface_id;

        let residues = self.rstar_samples.each_ref().map(|rstar_sample| {
            let energy_cache = graph
                .compute_onshell_energies(rstar_sample.loop_moms(), rstar_sample.external_moms());
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

            let terms =
                &counterterm.terms_in_counterterms[Into::<usize>::into(existing_esurface_id)];

            let eval_terms = terms.evaluate_from_esurface_cache(
                graph,
                numerator,
                rstar_sample,
                &esurface_cache,
                &energy_cache,
                settings,
            );

            if settings.general.debug > 0 {
                DEBUG_LOGGER.write("rstar_res", &eval_terms);
                DEBUG_LOGGER.write("rstar_energy_cache", &energy_cache);
                DEBUG_LOGGER.write("esurface_cache", &esurface_cache);
            }

            eval_terms * &multichanneling_factor / &rstar_energy_product
        });

        ResiudeEval {
            rstar_samples: self,
            residues,
        }
    }
}

struct ResiudeEval<'a, T: FloatLike> {
    rstar_samples: RstarSamples<'a, T>,
    residues: [Complex<F<T>>; 2],
}

impl<'a, T: FloatLike> ResiudeEval<'a, T> {
    fn compute_counterterm(self) -> Complex<F<T>> {
        let r = &self
            .rstar_samples
            .rstar_solutions
            .esurface_ct_builder
            .overlap_builder
            .hemispherical_radius;

        let e_cm = &self
            .rstar_samples
            .rstar_solutions
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .e_cm;

        let settings = self
            .rstar_samples
            .rstar_solutions
            .esurface_ct_builder
            .overlap_builder
            .counterterm_builder
            .settings;

        self.residues
            .into_iter()
            .enumerate()
            .map(|(sign_index, residue)| {
                let rstar = &self.rstar_samples.rstar_solutions.solutions[sign_index].solution;
                let rstar_derivative = &self.rstar_samples.rstar_solutions.solutions[sign_index]
                    .derivative_at_solution;
                let loop_number = self.rstar_samples.rstar_samples[sign_index]
                    .loop_moms()
                    .len();

                let jacobian_ratio = (rstar / r).abs().pow(3 * loop_number as u64 - 1);
                let uv_localisation = unnormalized_gaussian(
                    r,
                    rstar,
                    e_cm,
                    &F::<T>::from_f64(settings.subtraction.ct_settings.local_ct_width),
                );

                let singularity_dampener = if settings
                    .subtraction
                    .ct_settings
                    .dampen_integrable_singularity
                {
                    singularity_dampener(r, rstar)
                } else {
                    r.one()
                };

                let local_ct = &residue
                    * rstar_derivative.inv()
                    * (r - rstar).inv()
                    * &uv_localisation
                    * &jacobian_ratio
                    * &singularity_dampener;

                let h_function = utils::h(
                    rstar / r,
                    None,
                    settings
                        .subtraction
                        .ct_settings
                        .integrated_ct_sigma
                        .map(F::<T>::from_f64),
                    &settings.subtraction.ct_settings.integrated_ct_hfunction,
                );

                let minus_half = -(r.one() + r.one()).inv();

                let i = Complex::new_im(r.one());

                let integrated_ct = &i * r.PI() * &residue / rstar_derivative / rstar
                    * &h_function
                    * &jacobian_ratio
                    * minus_half;

                if settings.general.debug > 0 {
                    DEBUG_LOGGER.write("jacobian_ratio", &jacobian_ratio);
                    DEBUG_LOGGER.write("uv_dampener", &uv_localisation);
                    DEBUG_LOGGER.write("singularity_dampener", &singularity_dampener);
                    DEBUG_LOGGER.write("h_function", &h_function);
                    DEBUG_LOGGER.write("integrated_ct", &integrated_ct);
                    DEBUG_LOGGER.write("local_ct", &local_ct);
                }

                local_ct + integrated_ct
            })
            .reduce(|sum_ct, term| sum_ct + term)
            .unwrap_or(Complex::new(r.zero(), r.zero()))
    }
}
