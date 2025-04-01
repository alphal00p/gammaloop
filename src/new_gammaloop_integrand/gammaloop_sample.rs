use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use typed_index_collections::TiVec;

use crate::momentum::{Rotation, ThreeMomentum};
use crate::momentum_sample::{ExternalFourMomenta, MomentumSample, PolarizationVectors};
use crate::new_graph::{FeynmanGraph, Graph};
use crate::utils::{global_parameterize, FloatLike, F};
use crate::{
    disable, DependentMomentaConstructor, DiscreteGraphSamplingSettings, Externals, Polarizations,
    SamplingSettings, Settings,
};
use symbolica::numerical_integration::Sample;

use super::IntegrandType;

// helper functions can maybe moved to utils
#[inline]
fn unwrap_cont_sample(sample: &Sample<F<f64>>) -> &[F<f64>] {
    if let Sample::Continuous(_, xs) = sample {
        xs
    } else {
        panic!("Invalid sample structure")
    }
}

#[inline]
fn unwrap_single_discrete_sample(sample: &Sample<F<f64>>) -> (usize, &[F<f64>]) {
    if let Sample::Discrete(_, index, Some(cont_sample)) = sample {
        (*index, unwrap_cont_sample(cont_sample))
    } else {
        panic!("Invalid sample structure")
    }
}

#[inline]
fn unwrap_double_discrete_sample(sample: &Sample<F<f64>>) -> (usize, (usize, &[F<f64>])) {
    if let Sample::Discrete(_, index, Some(discrete_sample)) = sample {
        (*index, unwrap_single_discrete_sample(discrete_sample))
    } else {
        panic!("Invalid sample structure")
    }
}

/// Sample whose structure depends on the sampling settings, and enforces these settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GammaLoopSample<T: FloatLike> {
    Default(MomentumSample<T>),
    MultiChanneling {
        alpha: f64,
        sample: MomentumSample<T>,
    },
    DiscreteGraph {
        graph_id: usize,
        sample: DiscreteGraphSample<T>,
    },
}

impl GammaLoopSample<f64> {
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation: &Rotation,
        rotated_externals: Externals,
        rotated_polarizations: Polarizations,
    ) -> Self {
        if rotation.is_identity() {
            return self.clone();
        }

        let rotated_externals = rotated_externals.get_indep_externals().into();
        let rotated_polarizations = match rotated_polarizations {
            Polarizations::None => TiVec::new(),
            Polarizations::Constant { polarizations } => polarizations.into(),
        };
        match self {
            GammaLoopSample::Default(sample) => {
                GammaLoopSample::Default(sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ))
            }
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.get_rotated_sample_cached(
                        rotation,
                        rotated_externals,
                        rotated_polarizations,
                    ),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.get_rotated_sample(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ),
            },
        }
    }
}

impl<T: FloatLike> GammaLoopSample<T> {
    pub fn zero(&self) -> F<T> {
        match self {
            GammaLoopSample::Default(sample) => sample.zero(),
            GammaLoopSample::MultiChanneling { sample, .. } => sample.zero(),
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.zero(),
        }
    }

    #[allow(dead_code)]
    pub fn one(&self) -> F<T> {
        match self {
            GammaLoopSample::Default(sample) => sample.one(),
            GammaLoopSample::MultiChanneling { sample, .. } => sample.one(),
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.one(),
        }
    }

    /// Cast the sample to a different precision
    #[allow(dead_code)]
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> GammaLoopSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        match self {
            GammaLoopSample::Default(sample) => GammaLoopSample::Default(sample.cast_sample()),
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.cast_sample(),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.cast_sample(),
            },
        }
    }

    fn higher_precision(&self) -> GammaLoopSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        match self {
            GammaLoopSample::Default(sample) => GammaLoopSample::Default(sample.higher_precision()),
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.higher_precision(),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.higher_precision(),
            },
        }
    }

    #[allow(dead_code)]
    fn lower_precision(&self) -> GammaLoopSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        match self {
            GammaLoopSample::Default(sample) => GammaLoopSample::Default(sample.lower_precision()),
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.lower_precision(),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.lower_precision(),
            },
        }
    }

    /// Retrieve the default sample which is contained in all types
    #[inline]
    pub fn get_default_sample(&self) -> &MomentumSample<T> {
        match self {
            GammaLoopSample::Default(sample) => sample,
            GammaLoopSample::MultiChanneling { sample, .. } => sample,
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.get_default_sample(),
        }
    }
}

/// This sample is used when importance sampling over graphs is used.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DiscreteGraphSample<T: FloatLike> {
    Default(MomentumSample<T>),
    MultiChanneling {
        alpha: f64,
        sample: MomentumSample<T>,
    },
    /// This variant is equivalent to Default, but needs to be handled differently in the evaluation.
    Tropical(MomentumSample<T>),
    DiscreteMultiChanneling {
        alpha: f64,
        channel_id: usize,
        sample: MomentumSample<T>,
    },
}

impl<T: FloatLike> DiscreteGraphSample<T> {
    pub fn zero(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default(sample) => sample.zero(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.zero(),
            DiscreteGraphSample::Tropical(sample) => sample.zero(),
            DiscreteGraphSample::DiscreteMultiChanneling { sample, .. } => sample.zero(),
        }
    }

    pub fn one(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default(sample) => sample.one(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.one(),
            DiscreteGraphSample::Tropical(sample) => sample.one(),
            DiscreteGraphSample::DiscreteMultiChanneling { sample, .. } => sample.one(),
        }
    }
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation: &Rotation,
        rotated_externals: ExternalFourMomenta<F<T>>,
        rotated_polarizations: PolarizationVectors<Complex<F<T>>>,
    ) -> Self {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ))
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.get_rotated_sample_cached(
                        rotation,
                        rotated_externals,
                        rotated_polarizations,
                    ),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ))
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ),
            },
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> DiscreteGraphSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.cast_sample())
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.cast_sample(),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.cast_sample())
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.cast_sample(),
            },
        }
    }

    fn higher_precision(&self) -> DiscreteGraphSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.higher_precision())
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.higher_precision(),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.higher_precision())
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.higher_precision(),
            },
        }
    }

    fn lower_precision(&self) -> DiscreteGraphSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.lower_precision())
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.lower_precision(),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.lower_precision())
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.lower_precision(),
            },
        }
    }

    /// Retrieve the default sample which is contained in all types
    #[inline]
    fn get_default_sample(&self) -> &MomentumSample<T> {
        match self {
            DiscreteGraphSample::Default(sample) => sample,
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample,
            DiscreteGraphSample::Tropical(sample) => sample,
            DiscreteGraphSample::DiscreteMultiChanneling { sample, .. } => sample,
        }
    }
}

#[inline]
pub fn parameterize(
    sample_point: &Sample<F<f64>>,
    polarizations: &[Polarizations],
    dependent_momenta_constructor: DependentMomentaConstructor,
    settings: &Settings,
    graphs: Option<&[Graph]>, // this is only needed for tropical sampling
) -> Result<GammaLoopSample<f64>, String> {
    match &settings.sampling {
        SamplingSettings::Default => {
            let xs = unwrap_cont_sample(sample_point);
            Ok(GammaLoopSample::Default(default_parametrize(
                xs,
                dependent_momenta_constructor,
                polarizations,
                settings,
            )))
        }
        SamplingSettings::MultiChanneling(multichanneling_settings) => {
            let xs = unwrap_cont_sample(sample_point);
            Ok(GammaLoopSample::MultiChanneling {
                alpha: multichanneling_settings.alpha,
                sample: default_parametrize(
                    xs,
                    dependent_momenta_constructor,
                    polarizations,
                    settings,
                ),
            })
        }
        SamplingSettings::DiscreteGraphs(discrete_graph_settings) => {
            match discrete_graph_settings {
                DiscreteGraphSamplingSettings::Default => {
                    let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                    Ok(GammaLoopSample::DiscreteGraph {
                        graph_id,
                        sample: DiscreteGraphSample::Default(default_parametrize(
                            xs,
                            dependent_momenta_constructor,
                            polarizations,
                            settings,
                        )),
                    })
                }
                DiscreteGraphSamplingSettings::MultiChanneling(multichanneling_settings) => {
                    let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                    Ok(GammaLoopSample::DiscreteGraph {
                        graph_id,
                        sample: DiscreteGraphSample::MultiChanneling {
                            alpha: multichanneling_settings.alpha,
                            sample: default_parametrize(
                                xs,
                                dependent_momenta_constructor,
                                polarizations,
                                settings,
                            ),
                        },
                    })
                }
                DiscreteGraphSamplingSettings::TropicalSampling(tropical_sampling_settings) => {
                    todo!("add tropical sampling support");
                    disable! {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        let externals = &self.global_data.settings.kinematics.externals;

                        let graph = match &self.graph_integrands {
                            GraphIntegrands::Amplitude(graphs) => &graphs[graph_id].graph,
                            GraphIntegrands::CrossSection(_graphs) => unimplemented!(), //,
                        };

                        let sampler = graph
                                .derived_data
                                .as_ref()
                                .unwrap()
                                .tropical_subgraph_table
                                .as_ref()
                                .expect("No tropical subgraph table present, disable tropical sampling or regenerate process with table");

                        let edge_data = graph
                            .bare_graph
                            .get_loop_edges_iterator()
                            .map(|(edge_id, edge)| {
                                let mass = edge.particle.mass.value;
                                let mass_re = mass.map(|complex_mass| complex_mass.re);

                                let shift = utils::compute_shift_part(
                                    &graph.bare_graph.loop_momentum_basis.edge_signatures[edge_id]
                                        .external,
                                    &externals.get_indep_externals(),
                                )
                                .spatial;

                                let shift_momtrop = Vector::from_array([shift.px, shift.py, shift.pz]);

                                (mass_re, shift_momtrop)
                            })
                            .collect_vec();

                        let sampling_result_result = sampler.generate_sample_from_x_space_point(
                            xs,
                            edge_data,
                            &tropical_sampling_settings.into_tropical_sampling_settings(
                                self.global_data.settings.general.debug,
                            ),
                            &DEBUG_LOGGER,
                        );

                        let sampling_result = match sampling_result_result {
                            Ok(sampling_result) => sampling_result,
                            Err(_) => {
                                return Err(String::from("tropical sampling failed"));
                            }
                        };

                        let loop_moms = sampling_result
                            .loop_momenta
                            .into_iter()
                            .map(Into::<ThreeMomentum<F<f64>>>::into)
                            .collect_vec();

                        let default_sample = MomentumSample::new(
                            loop_moms,
                            externals,
                            sampling_result.jacobian * externals.pdf(xs),
                            &self.global_data.polarizations[0],
                            &graph.bare_graph.external_in_or_out_signature(),
                        );
                        Ok(GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::Tropical(default_sample),
                        })
                    }
                }
                DiscreteGraphSamplingSettings::DiscreteMultiChanneling(
                    multichanneling_settings,
                ) => {
                    let (graph_id, (channel_id, xs)) = unwrap_double_discrete_sample(sample_point);
                    Ok(GammaLoopSample::DiscreteGraph {
                        graph_id,
                        sample: DiscreteGraphSample::DiscreteMultiChanneling {
                            alpha: multichanneling_settings.alpha,
                            channel_id,
                            sample: default_parametrize(
                                xs,
                                dependent_momenta_constructor,
                                polarizations,
                                settings,
                            ),
                        },
                    })
                }
            }
        }
    }
}

/// Default parametrize is basically everything except tropical sampling.
#[inline]
fn default_parametrize(
    xs: &[F<f64>],
    dependent_momenta_constructor: DependentMomentaConstructor,
    polarizations: &[Polarizations],
    settings: &Settings,
) -> MomentumSample<f64> {
    let externals = &settings.kinematics.externals;

    let (loop_moms_vec, param_jacobian) =
        global_parameterize(xs, settings.kinematics.e_cm.square(), settings, false);

    let loop_moms = loop_moms_vec.into_iter().map(ThreeMomentum::from).collect();

    let jacobian = param_jacobian * externals.pdf(xs);

    MomentumSample::new(
        loop_moms,
        externals,
        jacobian,
        &polarizations[0],
        dependent_momenta_constructor,
    )
}
