use super::Settings;
use crate::graph::Graph;
use crate::integrands::HasIntegrand;
use crate::utils::{
    cast_complex, cast_lorentz_vector, format_for_compare_digits, global_parameterize, FloatLike,
};
use crate::{Precision, StabilityLevelSetting};
use itertools::Itertools;
use log::{debug, warn};
use lorentz_vector::LorentzVector;
use num::Complex;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum GammaLoopIntegrandType {
    Amplitude,
    CrossSection,
}

#[derive(Debug, Clone, Copy)]
pub struct GammaloopIntegrandTermStatistics {
    pub max_re_eval: f64,
    pub max_im_eval: f64,
    pub max_weight: f64,
    pub n_nan_evals: usize,
    pub zero_evals: usize,
    pub num_unstable_f32_points: usize,
    pub num_unstable_f64_points: usize,
    pub num_unstable_f128_points: usize,
    pub num_unstable_arb_points: usize,
}

impl GammaloopIntegrandTermStatistics {
    fn new() -> Self {
        Self {
            max_re_eval: 0.0,
            max_im_eval: 0.0,
            max_weight: 0.0,
            n_nan_evals: 0,
            zero_evals: 0,
            num_unstable_f32_points: 0,
            num_unstable_f64_points: 0,
            num_unstable_f128_points: 0,
            num_unstable_arb_points: 0,
        }
    }
}

// represents a single locally finite term
#[derive(Clone)]
pub struct GammaloopIntegrandTerm {
    pub graphs: Vec<Graph>,
    pub statistics: GammaloopIntegrandTermStatistics,
    pub n_dim: usize,
    pub n_moms: usize,
}

impl GammaloopIntegrandTerm {
    fn evaluate_in_momentum_space<T: FloatLike>(
        &self,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        use_ltd: bool,
    ) -> Complex<T> {
        if use_ltd {
            self.graphs
                .iter()
                .map(|graph| {
                    graph.evaluate_ltd_expression(loop_moms, external_moms)
                        / graph.compute_energy_product(loop_moms, external_moms)
                })
                .sum()
        } else {
            self.graphs
                .iter()
                .map(|graph| {
                    graph.evaluate_cff_expression(loop_moms, external_moms)
                        / graph.compute_energy_product(loop_moms, external_moms)
                })
                .sum()
        }
    }

    pub fn from_single_graph(graph: Graph) -> Self {
        Self {
            statistics: GammaloopIntegrandTermStatistics::new(),
            n_dim: graph.loop_momentum_basis.basis.len() * 3,
            n_moms: graph.loop_momentum_basis.basis.len(),
            graphs: vec![graph],
        }
    }
}

#[derive(Clone)]
pub struct GammaLoopIntegrand {
    pub settings: Settings,
    terms: Vec<GammaloopIntegrandTerm>,
    index_tree: IndexTree,
    integrand_type: GammaLoopIntegrandType,
}

impl Default for GammaLoopIntegrand {
    fn default() -> Self {
        let index_tree = IndexTree::Leaf(0);

        Self {
            settings: Settings::default(),
            terms: Vec::new(),
            index_tree,
            integrand_type: GammaLoopIntegrandType::Amplitude,
        }
    }
}

impl HasIntegrand for GammaLoopIntegrand {
    fn create_grid(&self) -> symbolica::numerical_integration::Grid<f64> {
        let term_dims = self.terms.iter().map(|term| term.n_dim).collect_vec();
        self.index_tree.create_grid(&self.settings, &term_dims)
    }

    #[allow(unused_variables)]
    fn evaluate_sample(
        &mut self,
        sample: &symbolica::numerical_integration::Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
    ) -> num::Complex<f64> {
        let (term, cont_sample) = self.index_tree.unpack_sample(sample);

        let x_space_point = if let Sample::Continuous(cont_weight, xs) = cont_sample.unwrap() {
            xs
        } else {
            panic!("expected continuous sample")
        };

        let (external_moms, pdf_weight) = self
            .settings
            .kinematics
            .externals
            .get_externals(x_space_point);

        debug!("external momenta: {:?}", external_moms);
        debug!("pdf weight: {}", pdf_weight);

        debug!("sampled x space point: {:?}", x_space_point);
        debug!("sampled term: {}", term);

        let n_moms = self.terms[term].n_moms;

        let (loop_moms_vec, jacobian) = match self.integrand_type {
            GammaLoopIntegrandType::Amplitude => global_parameterize(
                x_space_point,
                self.settings.kinematics.e_cm * self.settings.kinematics.e_cm,
                &self.settings,
                false,
            ),
            GammaLoopIntegrandType::CrossSection => {
                // perform t-scaling
                todo!("cross section computations not yet implemented")
            }
        };

        let loop_moms = loop_moms_vec
            .into_iter()
            .map(|p| LorentzVector::from_args(0.0, p[0], p[1], p[2]))
            .collect_vec();

        debug!("loop momenta: {:?}", loop_moms);

        debug!("jacobian: {}", jacobian);

        // setup the evaluation of the integrand in the different stability levels
        let mut results_of_stability_levels =
            Vec::with_capacity(self.settings.stability.levels.len());

        let rotation_method = self.settings.stability.rotation_axis.rotation_function();
        let wgt_ratio = wgt / self.terms[term].statistics.max_weight;

        let rotated_loop_moms = loop_moms.iter().map(&rotation_method).collect_vec();
        let rotated_external_moms = external_moms.iter().map(&rotation_method).collect_vec();

        let stab_iterator = if use_f128 {
            // overwrite the stability settings if use_f128 is enabled
            [StabilityLevelSetting {
                precision: Precision::Quad,
                required_precision_for_re: 1e-12,
                required_precision_for_im: 1e-12,
                escalate_for_large_weight_threshold: -1.,
                accepted_radius_in_x_range: [0., 1.],
            }]
            .iter()
        } else {
            self.settings.stability.levels.iter()
        };

        for stability_level in stab_iterator {
            if wgt_ratio > stability_level.escalate_for_large_weight_threshold {
                continue;
            }

            let (result, rotated_result) = self.evaluate_at_prec(
                term,
                &loop_moms,
                &external_moms,
                &rotated_loop_moms,
                &rotated_external_moms,
                stability_level.precision,
                true,
            );

            let error_real = if result.re.is_zero() && rotated_result.re.is_zero() {
                0.
            } else {
                ((result.re - rotated_result.re) / result.re).abs()
            };

            let error_im = if result.im.is_zero() && rotated_result.im.is_zero() {
                0.
            } else {
                ((result.im - rotated_result.im) / result.im).abs()
            };

            let real_stable = error_real < stability_level.required_precision_for_re;
            let im_stable = error_im < stability_level.required_precision_for_im;
            let stable = real_stable && im_stable;

            results_of_stability_levels.push((result, stable, stability_level.precision));

            if !stable {
                match stability_level.precision {
                    Precision::Single => {
                        self.terms[term].statistics.num_unstable_f32_points += 1;
                    }
                    Precision::Double => {
                        self.terms[term].statistics.num_unstable_f64_points += 1;
                    }
                    Precision::Quad => {
                        self.terms[term].statistics.num_unstable_f128_points += 1;
                    }
                    Precision::Arb(_) => {
                        self.terms[term].statistics.num_unstable_arb_points += 1;
                    }
                }

                let (real_colored, rotated_real_colored) =
                    format_for_compare_digits(result.re, rotated_result.re);

                let (im_colored, rotated_im_colored) =
                    format_for_compare_digits(result.im, rotated_result.im);

                debug!(
                    "unstable point: {:?}
                     result:         {} + {}i
                     rotated_result: {} + {}i",
                    x_space_point,
                    real_colored,
                    im_colored,
                    rotated_real_colored,
                    rotated_im_colored
                );
            } else {
                break;
            }
        }

        if wgt > self.terms[term].statistics.max_weight {
            self.terms[term].statistics.max_weight = wgt;
        }

        let (most_reliable_result, stable, precision) = results_of_stability_levels.last().unwrap();
        if !stable {
            warn!("Returning unstable point, consider adding more stability levels");
        }

        let res = most_reliable_result * jacobian * pdf_weight;

        // update the max evals
        if res.re.is_nan() {
            self.terms[term].statistics.n_nan_evals += 1;
        } else if res.re.abs() > self.terms[term].statistics.max_re_eval {
            self.terms[term].statistics.max_re_eval = res.re.abs();
        }

        if res.im.is_nan() {
            self.terms[term].statistics.n_nan_evals += 1;
        } else if res.im.abs() > self.terms[term].statistics.max_im_eval {
            self.terms[term].statistics.max_im_eval = res.im.abs();
        }

        res
    }

    fn get_event_manager_mut(&mut self) -> &mut crate::observables::EventManager {
        todo!()
    }

    // this maybe needs to change to Vec<usize>
    fn get_n_dim(&self) -> usize {
        match self.integrand_type {
            GammaLoopIntegrandType::Amplitude => {
                self.terms[0].graphs[0].loop_momentum_basis.basis.len() * 3
            }
            GammaLoopIntegrandType::CrossSection => {
                unimplemented!("cross section computations not yet implemented")
            }
        }
    }

    fn merge_results<I: HasIntegrand>(&mut self, _other: &mut I, _iter: usize) {
        todo!()
    }

    fn update_results(&mut self, _iter: usize) {
        todo!()
    }
}

#[allow(clippy::too_many_arguments)]
impl GammaLoopIntegrand {
    pub fn evaluate_at_prec(
        &self,
        term: usize,
        loop_moms: &[LorentzVector<f64>],
        external_moms: &[LorentzVector<f64>],
        rotated_loop_moms: &[LorentzVector<f64>],
        rotated_external_moms: &[LorentzVector<f64>],
        precision: Precision,
        measure_timing: bool,
    ) -> (Complex<f64>, Complex<f64>) {
        // cast the momenta to the relevant precision
        match precision {
            Precision::Single => {
                unimplemented!("From<f64> for f32 can't be implemented")
            }
            Precision::Double => {
                let before = if measure_timing {
                    Some(std::time::Instant::now())
                } else {
                    None
                };

                let result = self.terms[term].evaluate_in_momentum_space(
                    loop_moms,
                    external_moms,
                    self.settings.general.use_ltd,
                );

                let after = before.map(|t| t.elapsed());

                if measure_timing {
                    debug!(
                        "time for f64 evaluation: {:?} μs",
                        after.unwrap().as_micros()
                    );
                }

                let rotated_result = self.terms[term].evaluate_in_momentum_space(
                    rotated_loop_moms,
                    rotated_external_moms,
                    self.settings.general.use_ltd,
                );
                (result, rotated_result)
            }
            Precision::Quad => {
                let loop_moms_f128 = loop_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();
                let rotated_loop_moms_f128 = rotated_loop_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();
                let external_moms_f128 = external_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();
                let rotated_external_moms_f128 = rotated_external_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();

                let before = if measure_timing {
                    Some(std::time::Instant::now())
                } else {
                    None
                };

                let result = self.terms[term].evaluate_in_momentum_space(
                    &loop_moms_f128,
                    &external_moms_f128,
                    self.settings.general.use_ltd,
                );

                let after = before.map(|t| t.elapsed());
                if measure_timing {
                    debug!(
                        "time for f128 evaluation: {} μs",
                        after.unwrap().as_micros()
                    );
                }

                let rotated_result = self.terms[term].evaluate_in_momentum_space(
                    &rotated_loop_moms_f128,
                    &rotated_external_moms_f128,
                    self.settings.general.use_ltd,
                );

                debug!("result of term: {}, in f128: {:?}", term, result);

                // downcast back to f64
                (cast_complex(result), cast_complex(rotated_result))
            }
            Precision::Arb(_prec) => {
                unimplemented!("need better traits to use arb prec")
            }
        }
    }

    pub fn from_terms(
        terms: Vec<GammaloopIntegrandTerm>,
        integrand_type: GammaLoopIntegrandType,
        settings: Settings,
    ) -> Self {
        let index_tree = if terms.len() > 1 {
            let leafs = (0..terms.len()).map(IndexTree::Leaf).collect_vec();
            IndexTree::Node(leafs)
        } else {
            IndexTree::Leaf(0)
        };

        Self {
            terms,
            integrand_type,
            settings,
            index_tree,
        }
    }

    #[inline]
    pub fn create_sample(&self, index: &[usize], xs: Vec<f64>) -> Sample<f64> {
        self.index_tree.create_sample(index, xs)
    }
}

// struct to encode arbitrary discrete grid structure
#[derive(Clone, Debug)]
enum IndexTree {
    Leaf(usize),
    Node(Vec<IndexTree>),
}

impl IndexTree {
    // recursively build the grid according to the tree
    fn create_grid(
        &self,
        settings: &Settings,
        term_dims: &[usize],
    ) -> symbolica::numerical_integration::Grid<f64> {
        match self {
            IndexTree::Leaf(n) => Grid::Continuous(ContinuousGrid::new(
                term_dims[*n],
                settings.integrator.n_bins,
                settings.integrator.min_samples_for_update,
                settings.integrator.bin_number_evolution.clone(),
                settings.integrator.train_on_avg,
            )),
            IndexTree::Node(children) => {
                let grids = children
                    .iter()
                    .map(|child| Some(child.create_grid(settings, term_dims)))
                    .collect_vec();

                Grid::Discrete(DiscreteGrid::new(
                    grids,
                    settings.integrator.max_prob_ratio,
                    settings.integrator.train_on_avg,
                ))
            }
        }
    }

    fn unpack_sample<'a>(&self, sample: &'a Sample<f64>) -> (usize, Option<&'a Sample<f64>>) {
        match self {
            IndexTree::Leaf(n) => (*n, Some(sample)),
            IndexTree::Node(children) => {
                if let Sample::Discrete(_weight, i, Some(nested_sample)) = sample {
                    let child = &children[*i];
                    let (term, cont_sample) = child.unpack_sample(nested_sample);
                    (term, cont_sample)
                } else {
                    panic!("expected discrete sample")
                }
            }
        }
    }

    #[inline]
    #[allow(unused)]
    pub fn get_term(&self, index: &[usize]) -> usize {
        match self {
            IndexTree::Leaf(n) => *n,
            IndexTree::Node(children) => {
                let child = &children[index[0]];
                child.get_term(&index[1..])
            }
        }
    }

    #[inline]
    pub fn create_sample(&self, index: &[usize], xs: Vec<f64>) -> Sample<f64> {
        match self {
            IndexTree::Leaf(_n) => Sample::Continuous(1., xs),
            IndexTree::Node(children) => {
                let child = &children[index[0]];
                let nested_sample = child.create_sample(&index[1..], xs);
                Sample::Discrete(1., index[0], Some(Box::new(nested_sample)))
            }
        }
    }
}
