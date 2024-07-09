/// Counterterm for amplitudes with constant externals.
use itertools::Itertools;
use spenso::Complex;
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
    graph::{Graph, LoopMomentumBasis},
    momentum::{FourMomentum, ThreeMomentum},
    utils::{FloatLike, F},
    Settings,
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
        rotate_overlap_centers: bool,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let real_mass_vector = graph
            .get_real_mass_vector()
            .into_iter()
            .map(F::from_ff64)
            .collect_vec();

        let e_cm = F::from_ff64(settings.kinematics.e_cm);
        let sliver_width = F::from_f64(settings.subtraction.sliver_width);
        let esurfaces = &graph.get_cff().esurfaces;
        let lmb = &graph.loop_momentum_basis;

        let mut res = loop_momenta[0].px.zero();

        for (overlap_group, overlap_complement) in self
            .maximal_overlap
            .overlap_groups
            .iter()
            .zip(self.complements_of_overlap.iter())
        {
            let overlap = &overlap_group.existing_esurfaces;
            let center = &overlap_group.center;

            if settings.general.debug > 1 {
                println!("evaluating overlap structure: {:?}", overlap);
                println!("with center: {:?}", center);
            }

            let center_t = if rotate_overlap_centers {
                let rotation_function = settings.stability.rotation_axis.rotation_function();
                center
                    .iter()
                    .map(rotation_function)
                    .map(|c_vec| {
                        ThreeMomentum::new(
                            F::from_ff64(c_vec.px),
                            F::from_ff64(c_vec.py),
                            F::from_ff64(c_vec.px),
                        )
                    })
                    .collect_vec()
            } else {
                center
                    .iter()
                    .map(|c_vec| {
                        ThreeMomentum::new(
                            F::from_ff64(c_vec.px),
                            F::from_ff64(c_vec.py),
                            F::from_ff64(c_vec.px),
                        )
                    })
                    .collect_vec()
            };

            let shifted_momenta = loop_momenta
                .iter()
                .zip(center_t.iter())
                .map(|(momentum, center)| momentum.clone() - center.clone())
                .collect_vec();

            let (hemispherical_unit_shifted_momenta, hemispherical_radius) =
                normalize_momenta(&shifted_momenta);

            if settings.general.debug > 1 {
                println!("shifted_momenta: {:?}", shifted_momenta);
                println!("hemihyperradius: {:?}", hemispherical_radius);
            }

            //println!("hemispherical_radius: {}", hemispherical_radius);

            for existing_esurface_id in overlap.iter() {
                if settings.general.debug > 1 {
                    println!("subtracting esurface: {:?}", existing_esurface_id);
                }

                let esurface = &esurfaces[self.existing_esurfaces[*existing_esurface_id]];
                // solve the radius
                let radius_guess = esurface.get_radius_guess(
                    &hemispherical_unit_shifted_momenta,
                    external_momenta,
                    lmb,
                );

                if settings.general.debug > 1 {
                    println!("radius_guess: {}", radius_guess)
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

                let (radius_plus, derivative_at_r_plus) = newton_iteration_and_derivative(
                    &radius_guess,
                    function,
                    &F::from_f64(TOLERANCE),
                    MAX_ITERATIONS,
                    &e_cm,
                );

                let (radius_minus, derivative_at_r_minus) = newton_iteration_and_derivative(
                    &(-radius_guess),
                    function,
                    &F::from_f64(TOLERANCE),
                    MAX_ITERATIONS,
                    &e_cm,
                );

                if settings.general.debug > 1 {
                    println!(
                        "radius_plus: {}, radius_minus: {}",
                        radius_plus, radius_minus
                    );
                }

                let (r_plus_eval, r_minus_eval) = (
                    self.radius_star_eval(
                        &radius_plus,
                        &hemispherical_radius,
                        &hemispherical_unit_shifted_momenta,
                        &center_t,
                        graph,
                        external_momenta,
                        esurfaces,
                        overlap_complement,
                        *existing_esurface_id,
                        &e_cm,
                        &sliver_width,
                        settings.subtraction.dampen_integrable_singularity,
                    ),
                    self.radius_star_eval(
                        &radius_minus,
                        &hemispherical_radius,
                        &hemispherical_unit_shifted_momenta,
                        &center_t,
                        graph,
                        external_momenta,
                        esurfaces,
                        overlap_complement,
                        *existing_esurface_id,
                        &e_cm,
                        &sliver_width,
                        settings.subtraction.dampen_integrable_singularity,
                    ),
                );

                let ct_plus = r_plus_eval
                    * ((&hemispherical_radius - &radius_plus) * derivative_at_r_plus).inv();

                let ct_minus = r_minus_eval
                    * ((&hemispherical_radius - &radius_minus) * derivative_at_r_minus).inv();

                res += ct_plus + ct_minus;
            }
        }

        // match the complex prefactor off cff
        let loop_number = graph.loop_momentum_basis.basis.len();
        let internal_vertex_number = graph.vertices.len() - graph.external_connections.len();
        let zero = res.zero();

        let prefactor = Complex::new(zero, res.one()).pow(loop_number as u64)
            * Complex::new(-res.one(), res.zero()).pow(internal_vertex_number as u64 - 1);

        Complex::new(res.clone(), res.zero()) * prefactor
    }

    #[allow(clippy::too_many_arguments)]
    fn radius_star_eval<T: FloatLike>(
        &self,
        rstar: &F<T>,
        radius: &F<T>,
        unit_loop_momenta: &[ThreeMomentum<F<T>>],
        center: &[ThreeMomentum<F<T>>],
        graph: &Graph,
        external_momenta: &[FourMomentum<F<T>>],
        esurfaces: &EsurfaceCollection,
        overlap_complement: &[ExistingEsurfaceId],
        existing_esurface_id: ExistingEsurfaceId,
        e_cm: &F<T>,
        sliver_width: &F<T>,
        dampen_integrable_singularity: bool,
    ) -> F<T> {
        if (radius - rstar).abs() > sliver_width * e_cm {
            return radius.zero();
        }

        let loop_number = unit_loop_momenta.len();

        let jacobian_ratio = (rstar / radius).abs().powi(3 * loop_number as i32 - 1);

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
            .fold(radius.one(), |acc, e| acc * e);

        let multichanneling_denominator =
            self.evaluate_multichanneling_denominator(&esurface_cache);

        let multichanneling_numerator_root =
            overlap_complement.iter().fold(radius.one(), |acc, id| {
                acc * &esurface_cache[self.existing_esurfaces[*id]]
            });

        let multichanneling_factor = &multichanneling_numerator_root
            * &multichanneling_numerator_root
            / multichanneling_denominator;

        let terms = &self.terms_in_counterterms[Into::<usize>::into(existing_esurface_id)];

        let eval_terms = terms.evaluate_from_esurface_cache(&esurface_cache, &energy_cache);

        let r_minus_rstar = radius - rstar;
        let dampening_factor = (-&r_minus_rstar * &r_minus_rstar / (e_cm * e_cm)).exp(); // unnormalized such that the exponential is 1 at r = r*

        let singularity_dampener = if dampen_integrable_singularity {
            (radius.one()
                - (rstar * rstar * rstar * rstar
                    / (rstar * rstar - &r_minus_rstar * &r_minus_rstar).powi(2)))
            .exp()
        } else {
            radius.one()
        };

        jacobian_ratio
            * multichanneling_factor
            * eval_terms
            * dampening_factor
            * singularity_dampener
            / rstar_energy_product
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

/// root finding, also returns the derivative at the root, so that we don't have to recompute it
fn newton_iteration_and_derivative<T: FloatLike>(
    guess: &F<T>,
    f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
    tolerance: &F<T>,
    max_iterations: usize,
    e_cm: &F<T>,
) -> (F<T>, F<T>) {
    let mut x = guess.clone();
    let (mut val_f_x, mut val_df_x) = f_x_and_df_x(&x);

    let mut iteration = 0;

    while iteration < max_iterations && val_f_x.abs() > guess.epsilon() * tolerance * e_cm {
        x -= val_f_x / val_df_x;
        (val_f_x, val_df_x) = f_x_and_df_x(&x);
        // println!("x: {}, val_f_x: {}, val_df_x: {}", x, val_f_x, val_df_x);

        iteration += 1;
    }

    (x, val_df_x)
}
