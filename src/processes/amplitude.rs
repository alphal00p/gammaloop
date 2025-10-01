use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::Write,
    iter,
    ops::Deref,
    path::Path,
};

use ahash::AHashSet;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use momtrop::SampleGenerator;

use idenso::color::ColorSimplifier;
use rayon::{
    iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
    ThreadPool,
};
use spenso::network::{Sequential, SmallestDegree};
use tracing::{info_span, instrument};
use tracing_indicatif::{span_ext::IndicatifSpanExt, style::ProgressStyle};
use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    cff::{
        esurface::GroupEsurfaceId,
        expression::{
            AmplitudeOrientationID, CFFExpression, OrientationData, SubgraphOrientationID,
        },
        generation::{generate_cff_expression, get_orientations_from_subgraph},
    },
    gammaloop_integrand::{
        amplitude_integrand::{AmplitudeGraphTerm, AmplitudeIntegrand, AmplitudeIntegrandData},
        LmbMultiChannelingSetup,
    },
    graph::{GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex, LoopMomentumBasis},
    model::ArcParticle,
    momentum_sample::ExternalIndex,
    numerator::symbolica_ext::AtomCoreExt,
    settings::{runtime::LockedRuntimeSettings, GlobalSettings},
    signature::SignatureLike,
    status_debug,
    subtraction::amplitude_counterterm::AmplitudeCountertermAtom,
    utils::{symbolica_ext::LOGPRINTOPTS, Length, GS, TENSORLIB, W_},
    uv::UltravioletGraph,
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::{eyre, Context};
use itertools::Itertools;
use linnet::{
    half_edge::involution::{HedgePair, Orientation},
    parser::DotGraph,
};
use log::{debug, info};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use typed_index_collections::{ti_vec, TiVec};

use crate::{
    cff::esurface::EsurfaceID,
    gammaloop_integrand::NewIntegrand,
    graph::{FeynmanGraph, Graph},
    model::Model,
    settings::global::GenerationSettings,
};

use crate::graph::parse::complete_group_parsing;

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Amplitude {
    pub name: String,
    pub integrand: Option<NewIntegrand>,
    pub graphs: Vec<AmplitudeGraph>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub external_particles: Vec<ArcParticle>,
    pub external_signature: SignatureLike<ExternalIndex>,
    pub group_derived_data: TiVec<GroupId, GroupDerivedData>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GroupDerivedData {
    pub esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
    pub esurface_atoms: TiVec<GroupEsurfaceId, Atom>,
}

impl Amplitude {
    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
          )
    )]
    pub(crate) fn warm_up(&mut self, model: &Model) -> Result<()> {
        if let Some(integrand) = &mut self.integrand {
            integrand.warm_up(model)
        } else {
            Err(eyre!(
                "Cannot warm up amplitude {} without integrand",
                self.name
            ))
        }
    }

    #[instrument(
          skip_all,
          fields(
              path = %path.as_ref().display(),
          )
    )]
    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("amp.bin"))?;
        let (mut amp, _): (Self, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        if path.as_ref().join("integrand").exists() {
            let integrand = AmplitudeIntegrand::load(path.as_ref().join("integrand"), context)?;
            amp.integrand = Some(NewIntegrand::Amplitude(integrand));
        }

        Ok(amp)
    }

    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
          )
    )]
    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        info!("Compiling amplitude {}", self.name);
        let p = path.as_ref().join(format!("amp_{}", self.name));

        let r = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                p.display()
            )
        });
        if override_existing {
            r?;
        }
        if let Some(integrand) = &mut self.integrand {
            integrand.compile(&p, override_existing, settings, thread_pool)?;
        };
        Ok(())
    }

    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
              path = %path.as_ref().display(),
          )
    )]
    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let p = path.as_ref().join(format!("amp_{}", self.name));

        let r = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                p.display()
            )
        });
        if override_existing {
            r?;
        }

        let integrand = self.integrand.take();
        if let Some(integrand) = &integrand {
            integrand.save(&p, override_existing)?;
        };

        let binary = bincode::encode_to_vec(&(*self), bincode::config::standard())?;
        if override_existing {
            fs::write(p.join("amp.bin"), binary)?;
        } else {
            let mut file = File::create_new(p.join("amp.bin"))?;
            file.write(&binary)?;
        }

        self.integrand = integrand;
        Ok(())
    }

    #[instrument(
        skip_all,
        fields(
             amplitude.name = %self.name,
        )
    )]
    pub fn preprocess(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        // preprocess each graph individually

        let preprocess_span = info_span!("Preprocessing graphs", indicatif.pb_show = true);
        preprocess_span.pb_set_style(&ProgressStyle::with_template(
            "{wide_bar} {pos}/{len} {msg}",
        )?);
        preprocess_span.pb_set_length(self.graphs.len() as u64);
        preprocess_span.pb_set_message("Preprocessing graphs");
        preprocess_span.pb_set_finish_message("all graphs preprocessed");

        let preprocess_span_enter = preprocess_span.enter();

        thread_pool.install(|| {
            self.graphs.par_iter_mut().try_for_each(|amplitude_graph| {
                let ok = amplitude_graph.preprocess(model, settings);
                preprocess_span.pb_inc(1);

                ok
            })
        })?;

        drop(preprocess_span_enter);
        drop(preprocess_span);

        self.generate_grouped_derived_data()?;

        Ok(())
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    pub fn build_integrand(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        let terms: Result<Vec<_>> = thread_pool.install(|| {
            self.graphs
                .par_iter_mut()
                .enumerate()
                .map(|(graph_id, graph)| {
                    let group_id = graph.graph.group_id.unwrap(); // should always be set
                    let esurface_map = &self.group_derived_data[group_id].esurface_map;
                    let group_pos = self.graph_group_structure[group_id]
                        .find_position(graph_id)
                        .unwrap();

                    graph.generate_term_for_graph(
                        model,
                        group_pos,
                        esurface_map.clone(),
                        global_settings,
                    )
                })
                .collect()
        });

        let amplitude_integrand = AmplitudeIntegrand {
            settings: runtime_default.into_with_modified_kinematics(
                &self.external_signature,
                &self.graphs[0].graph.get_external_masses(model),
            )?,
            data: AmplitudeIntegrandData {
                name: self.name.clone(),
                rotations: None,
                loop_cache_id: 0,
                external_cache_id: 0,
                base_external_cache_id: 0,
                graph_terms: terms?,
                external_signature: self.external_signature.clone(),
                graph_group_structure: self.graph_group_structure.clone(),
                group_derived_data: self.group_derived_data.clone(),
            },
        };
        self.integrand = Some(NewIntegrand::Amplitude(amplitude_integrand));
        Ok(())
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        for graph in &self.graphs {
            graph.write_dot(writer)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub fn generate_grouped_derived_data(&mut self) -> Result<()> {
        // for each group we must collect all inequivalent esurfaces.

        let group_derived_data = self
            .graph_group_structure
            .iter()
            .map(|group| {
                let mut group_esurface_structure =
                    BTreeMap::<Atom, TiVec<GraphGroupPosition, Option<EsurfaceID>>>::default();

                for (graph_group_position, graph_id) in group.iter_enumerated() {
                    let amplitude_graph = &self.graphs[graph_id];
                    let lmb_reps = amplitude_graph.graph.integrand_replacement(
                        &amplitude_graph.graph.full_filter(),
                        &amplitude_graph.graph.loop_momentum_basis,
                        &[W_.x___],
                    );

                    let esurfaces = &amplitude_graph
                        .derived_data
                        .cff_expression
                        .as_ref()
                        .unwrap()
                        .surfaces
                        .esurface_cache;

                    for (esurface_id, esurface) in esurfaces.iter_enumerated() {
                        let esurface_atom = esurface.lmb_atom(&amplitude_graph.graph, &lmb_reps);

                        group_esurface_structure
                            .entry(esurface_atom)
                            .or_insert(ti_vec![None; group.len()])[graph_group_position] =
                            Some(esurface_id);
                    }
                }

                let (surface_atoms, esurface_map) = group_esurface_structure.into_iter().unzip();

                GroupDerivedData {
                    esurface_map,
                    esurface_atoms: surface_atoms,
                }
            })
            .collect::<TiVec<GroupId, _>>();

        Ok(self.group_derived_data = group_derived_data)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct AmplitudeGraph {
    pub graph: Graph,
    pub derived_data: AmplitudeDerivedData,
}

impl AmplitudeGraph {
    pub(crate) fn new(graph: Graph) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {
                all_mighty_integrand: Atom::Zero,
                cff_expression: None,

                lmbs: None,
                tropical_sampler: None,
                multi_channeling_setup: None,
                threshold_counterterms: TiVec::new(),
            },
        }
    }
}

impl AmplitudeGraph {
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer)
    }

    pub(crate) fn generate_cff(&mut self) -> Result<()> {
        let shift_rewrite = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_expression = generate_cff_expression(
            &self.graph.underlying,
            &shift_rewrite,
            &self.graph.dummy_list(),
        )?;
        self.derived_data.cff_expression = Some(cff_expression);

        Ok(())
    }

    //Stage 1
    pub(crate) fn preprocess(
        &mut self,
        _model: &Model,
        settings: &GenerationSettings,
    ) -> Result<()> {
        status_debug!("Generating Cff");
        self.generate_cff()?;
        status_debug!("Building Parametric Integrand");
        self.build_parametric_integrand(settings)?;
        status_debug!("Building Tropical Sampler");
        self.build_tropical_sampler(settings)?;
        status_debug!("Building Loop Momentum Bases");
        self.build_lmbs();
        status_debug!("Building Multi-Channeling Channels");
        self.build_multi_channeling_channels();

        if settings.enable_thresholds {
            status_debug!("Building Threshold Counterterms");
            self.derived_data.threshold_counterterms =
                self.build_threshold_counterterm_parametric_integrand(settings)?;
        }

        Ok(())
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    //  pub fn ct_params(&self, model: &Model) -> ParamBuilder<f64> {
    //      let mut param_builder = self.param_builder_core(model);
    //      param_builder.uv_damp_atom(vec![
    //          Atom::var(GS.uv_damp_plus),
    //          Atom::var(GS.uv_damp_minus),
    //      ]);
    //      param_builder.derivative_at_tstar_atom(Atom::var(GS.deta));
    //      param_builder.radius_atom(Atom::var(GS.radius));
    //      param_builder.radius_star_atom(Atom::var(GS.radius_star));
    //      param_builder.h_function_atom(Atom::var(GS.hfunction));
    //      param_builder
    //  }

    //  pub(crate) fn param_builder_core(&self, model: &Model) -> ParamBuilder<f64> {
    //      // the float type does not matter here
    //      let mut param_builder = ParamBuilder::<f64>::new_empty();

    //      // this is wrong if we allow for vacuum graphs
    //      param_builder.external_energies_atom(&self.graph);
    //      param_builder.orientation_params(&self.graph);

    //      // param_builder.polarizations(&self.graph);
    //      // param_builder.polar

    //      // spatial components of external momenta
    //      param_builder.external_spatial_atom(&self.graph);
    //      param_builder.polarization_params(&self.graph);
    //      // spatial EMR
    //      param_builder.emr_spatial_atom(&self.graph);

    //      param_builder.model_parameters_atom(model);

    //      param_builder.m_uv_atom(Atom::var(GS.m_uv));

    //      param_builder.mu_r_sq_atom(Atom::var(GS.mu_r_sq));

    //      self.add_function_map(&mut param_builder);

    //      param_builder
    //  }

    // pub fn fill_in_params(
    //     &self,
    //     param_builder: &mut ParamBuilder,
    //     model: &Model,
    //     settings: &Settings,
    // ) {
    //     param_builder.external_energies_value(momentum_sample);
    // }

    // fn get_eager_const_map(&self)->HashM

    fn add_additional_factors_to_cff_atom(&self, cff_atom: &Atom) -> Atom {
        // let inverse_energy_product = self.graph.underlying.get_cff_inverse_energy_product();
        let factors_of_pi = (Atom::var(GS.pi) * 2).npow(3 * self.graph.get_loop_number() as i64);

        let result = cff_atom / factors_of_pi;
        // debug!("result: {}", result);
        result
    }

    fn new_vakint(&self) -> Vakint {
        Vakint::new(Some(VakintSettings {
            allow_unknown_integrals: false,
            evaluation_order: EvaluationOrder::alphaloop_only(),
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            run_time_decimal_precision: 32,
            number_of_terms_in_epsilon_expansion: self.graph.n_loops(&self.graph.no_dummy()) as i64
                + 1,
            // temporary_directory: Some("./form".into()),
            mu_r_sq_symbol: GS.mu_r_sq.get_name().to_string(),
            ..VakintSettings::default()
        }))
        .unwrap()
    }

    pub(crate) fn build_parametric_integrand(
        &mut self,
        settings: &GenerationSettings,
    ) -> Result<()> {
        self.derived_data.all_mighty_integrand =
            self.build_original_parametric_integrand(settings)?;
        Ok(())
    }

    fn build_threshold_counterterm_parametric_integrand(
        &self,
        settings: &GenerationSettings,
    ) -> Result<TiVec<EsurfaceID, AmplitudeCountertermAtom>> {
        let global_num = self.graph.global_network();

        let mut counterterms: TiVec<EsurfaceID, AmplitudeCountertermAtom> = TiVec::new();
        let canonize_esurface = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        for (esurface_id, esurface) in self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache
            .iter_enumerated()
        {
            if esurface.external_shift.is_empty() {
                // these will never satsify the threshold condition
                // so we can skip them
                counterterms.push(AmplitudeCountertermAtom {
                    parametric_local: Atom::new(),
                    parametric_integrated: Atom::new(),
                });
                continue;
            }

            let (circled, complement) = esurface.get_subgraph_components(&self.graph.underlying);
            let edges_in_cut = esurface.bitvec(&self.graph.underlying);

            let orientations = self
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .get_orientations_with_esurface(esurface_id);

            let first_orientation = &self
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations[orientations[0]]
                .data
                .orientation;

            let orientation_of_edges_in_esurface = esurface
                .energies
                .iter()
                .map(|e| first_orientation[*e])
                .collect_vec();

            assert!(orientations.iter().all(|o| {
                let or = &self
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap()
                    .orientations[*o]
                    .data
                    .orientation;

                let orientation_of_esurface_in_this_orientation =
                    esurface.energies.iter().map(|e| or[*e]).collect_vec();

                if orientation_of_edges_in_esurface != orientation_of_esurface_in_this_orientation {
                    println!("{:?}", orientation_of_edges_in_esurface);
                    println!("{:?}", orientation_of_esurface_in_this_orientation);
                    println!("esurface shift: {:?}", esurface.external_shift);
                    false
                } else {
                    true
                }
            }));

            let circled_wood = self.graph.wood(&circled);
            let complement_wood = self.graph.wood(&complement);

            let mut circled_forest =
                circled_wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

            let mut complement_forest =
                complement_wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

            let reverse_dangling = esurface
                .energies
                .iter()
                .zip(orientation_of_edges_in_esurface)
                .filter_map(|(e, o)| {
                    if o == Orientation::Reversed {
                        Some(*e)
                    } else {
                        None
                    }
                })
                .collect_vec();

            let circled_orientations =
                get_orientations_from_subgraph(&self.graph.underlying, &circled, &reverse_dangling)
                    .into_iter()
                    .map(|cff_graph| cff_graph.global_orientation)
                    .filter(|a| settings.orientation_pattern.filter(a))
                    .collect::<TiVec<SubgraphOrientationID, _>>();

            let complement_orientations = get_orientations_from_subgraph(
                &self.graph.underlying,
                &complement,
                &reverse_dangling,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .filter(|a| settings.orientation_pattern.filter(a))
            .collect::<TiVec<SubgraphOrientationID, _>>();

            let vakint = self.new_vakint();

            // println!("//Circled\n{}", self.graph.dot(&circled));
            // println!("//Complement\n{}", self.graph.dot(&complement));
            circled_forest.compute(
                &self.graph,
                &circled,
                &vakint,
                &circled_orientations,
                &canonize_esurface,
                &esurface.energies,
                &settings.uv,
            );

            complement_forest.compute(
                &self.graph,
                &complement,
                &vakint,
                &complement_orientations,
                &canonize_esurface,
                &esurface.energies,
                &settings.uv,
            );

            let circled_expr =
                circled_forest.orientation_parametric_expr(Some(&edges_in_cut), &self.graph);

            let complement_expr = complement_forest.orientation_parametric_expr(None, &self.graph);

            // println!("Circled Expression Network:");
            // println!("{}", circled_expr.dot_pretty());

            // println!("Complement Expression Network:");
            // println!("{}", complement_expr.dot_pretty());

            let mut product = circled_expr * complement_expr * global_num.clone();

            product
                .execute::<Sequential, SmallestDegree, _, _>(TENSORLIB.read().unwrap().deref())
                .unwrap();
            // println!("{}", product.dot_pretty());

            let scalar: Atom = product
                .result_scalar()
                .with_context(|| "in building threshold counterterm")?
                .into();

            let counterterm = scalar
                .unwrap_function(GS.color_wrap)
                .simplify_color()
                .replace(function!(GS.energy, W_.x_))
                .with(function!(GS.ose, W_.x_));

            let loop_3 = self.graph.get_loop_number() as i64 * 3;

            let grad_eta = Atom::var(GS.deta);
            let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3);
            let i = Atom::i();

            let radius = Atom::var(GS.radius);
            let radius_star = Atom::var(GS.radius_star);
            let uv_damp_plus = Atom::var(GS.uv_damp_plus);
            let uv_damp_minus = Atom::var(GS.uv_damp_minus);
            let hfunction = Atom::var(GS.hfunction);

            let delta_r_plus = &radius - &radius_star;
            let delta_r_minus = -&radius - &radius_star;

            let jacobian_ratio = (&radius_star / &radius).npow(loop_3 - 1);

            let local_prefactor = &jacobian_ratio / &factors_of_pi / &grad_eta
                * (uv_damp_plus / delta_r_plus + uv_damp_minus / delta_r_minus);

            let integrated_prefactor =
                i * Atom::var(GS.pi) * &jacobian_ratio * hfunction / factors_of_pi / grad_eta;

            let local_counterterm = local_prefactor * &counterterm;
            let integrated_counterterm = integrated_prefactor * &counterterm;

            // println!("CounterTerm{}", counterterm);
            counterterms.push(AmplitudeCountertermAtom {
                parametric_local: local_counterterm,
                parametric_integrated: integrated_counterterm,
            });
        }

        // let ct_4 = &counterterms[EsurfaceID::from(4)];
        // panic!("Counterterm 4: {}", ct_4.parametric_local);

        Ok(counterterms)
    }

    fn build_original_parametric_integrand(&self, settings: &GenerationSettings) -> Result<Atom> {
        let wood = self.graph.wood(&self.graph.no_dummy());
        debug!(
            "Wood for {}{}",
            self.graph.name,
            wood.show_graphs(&self.graph)
        );
        // debug!("{}", wood.dot(&self.graph));
        let mut forest = wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

        debug!("Forest: {}", forest.graphs(&self.graph));

        let canonize_esurface = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let orientations: TiVec<AmplitudeOrientationID, OrientationData> = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|a| a.data.clone())
            .filter(|a| settings.orientation_pattern.filter(a))
            .collect();

        let vakint = self.new_vakint();

        forest.compute(
            &self.graph,
            &self.graph.no_dummy(),
            &vakint,
            &orientations,
            &canonize_esurface,
            &[],
            &settings.uv,
        );

        let global_num = self.graph.global_network();
        let mut full = forest.orientation_parametric_expr(None, &self.graph);

        full *= global_num;

        full.execute::<Sequential, SmallestDegree, _, _>(TENSORLIB.read().unwrap().deref())
            .unwrap();

        let mut scalar: Atom = full
            .result_scalar()
            .with_context(|| format!("Failed to get scalar from network when building original paramteric integrand.")).with_note(||format!("Network: \n{}\nGraph:\n{}", full.dot_pretty(),DotGraph::from(&self.graph).debug_dot()))?
            .into();

        debug!(
            "All parametric before color atom:{}",
            scalar.printer(LOGPRINTOPTS)
        );
        scalar = scalar.unwrap_function(GS.color_wrap).simplify_color();

        scalar = self.add_additional_factors_to_cff_atom(&scalar);

        debug!(
            "All parametric integrand atom:{}",
            scalar.printer(LOGPRINTOPTS)
        );

        Ok(scalar)
    }

    fn build_lmbs(&mut self) {
        let lmbs = self
            .graph
            .generate_loop_momentum_bases(&self.graph.no_dummy());

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_tropical_sampler(&mut self, process_settings: &GenerationSettings) -> Result<()> {
        if process_settings
            .tropical_subgraph_table
            .disable_tropical_generation
        {
            debug!("Tropical subgraph table generation is disabled.");
            return Ok(());
        }
        let num_virtual_loop_edges = self.graph.iter_loop_edges().count();
        let num_loops = self.graph.loop_momentum_basis.loop_edges.len();
        let target_omega = process_settings.tropical_subgraph_table.target_omega;

        let weight = (target_omega + (3 * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64;

        debug!(
            "Building tropical subgraph table with all edge weights set to: {}",
            weight
        );

        let tropical_edges = self
            .graph
            .iter_loop_edges()
            .map(|(pair, _edge_id, edge)| {
                let is_massive = edge.data.particle.is_massive();

                let vertices = match pair {
                    HedgePair::Paired { source, sink } => (
                        self.graph.underlying.node_id(source).0 as u8,
                        self.graph.underlying.node_id(sink).0 as u8,
                    ),
                    _ => unreachable!(),
                };

                momtrop::Edge {
                    is_massive,
                    weight,
                    vertices,
                }
            })
            .collect_vec();

        let mut external_vertices_pool = AHashSet::new();

        for (pair, _, _) in self.graph.iter_non_loop_edges() {
            match pair {
                HedgePair::Paired { source, sink } => {
                    let source_id = self.graph.underlying.node_id(source).0 as u8;
                    let sink_id = self.graph.underlying.node_id(sink).0 as u8;

                    external_vertices_pool.insert(source_id);
                    external_vertices_pool.insert(sink_id);
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let id = self.graph.underlying.node_id(hedge).0 as u8;
                    external_vertices_pool.insert(id);
                }
                _ => unreachable!(),
            }
        }

        let mut external_vertices = vec![];

        for tropical_edge in &tropical_edges {
            if external_vertices_pool.contains(&tropical_edge.vertices.0) {
                external_vertices.push(tropical_edge.vertices.0);
            }

            if external_vertices_pool.contains(&tropical_edge.vertices.1) {
                external_vertices.push(tropical_edge.vertices.1);
            }
        }

        let tropical_graph = momtrop::Graph {
            edges: tropical_edges,
            externals: external_vertices,
        };

        let loop_part = self
            .graph
            .iter_loop_edges()
            .map(|(_, edge_id, _edge)| {
                self.graph.loop_momentum_basis.edge_signatures[edge_id]
                    .internal
                    .clone()
                    .to_momtrop_format()
            })
            .collect_vec();

        let sampler = tropical_graph
            .build_sampler(loop_part)
            .map_err(|e| eyre!(e))?;

        Ok(self.derived_data.tropical_sampler = Some(sampler))
    }

    // Expects cff_expression, esurface_data,
    #[instrument(
          name = "generate_term_for_graph",
          level = "info",
          skip(self, model, global_settings),
          fields(
              graph.name = %self.graph.name
          ),
          err
      )]
    fn generate_term_for_graph(
        &self,
        model: &Model,
        own_group_position: GraphGroupPosition,
        esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
        global_settings: &GlobalSettings,
    ) -> Result<AmplitudeGraphTerm> {
        AmplitudeGraphTerm::from_amplitude_graph(
            self,
            own_group_position,
            esurface_map,
            model,
            global_settings,
        )
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeDerivedData {
    pub all_mighty_integrand: Atom,
    pub threshold_counterterms: TiVec<EsurfaceID, AmplitudeCountertermAtom>,

    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub cff_expression: Option<CFFExpression<AmplitudeOrientationID>>,
}

pub trait AmplitudeState:
    Clone + std::fmt::Debug + bincode::Encode + for<'a> bincode::Decode<GammaLoopContextContainer<'a>>
{
}
impl AmplitudeState for () {}

#[derive(Clone, Encode, Decode, Debug)]
pub struct Processed {}
impl AmplitudeState for Processed {}

// #[derive(Clone, Encode, Decode, Debug)]
// pub struct ReadyForTerm {}
// impl AmplitudeState for ReadyForTerm {}

impl Amplitude {
    pub fn from_dot_string<Str: AsRef<str>>(s: Str, name: String, model: &Model) -> Result<Self> {
        let graphs = Graph::from_string(s, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub fn from_dot_file<'a, P>(p: P, name: String, model: &Model) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let graphs = Graph::from_file(p, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub fn from_graph_list(name: impl ToString, mut graphs: Vec<Graph>) -> Result<Self> {
        let mut amplitude: Amplitude = Amplitude::new(name);
        amplitude.graph_group_structure = complete_group_parsing(&mut graphs)?;

        for amplitude_graph in graphs {
            amplitude.add_graph(amplitude_graph)?;
        }
        Ok(amplitude)
    }

    pub(crate) fn new(name: impl ToString) -> Self {
        Self {
            integrand: None,
            name: name.to_string(),
            graphs: vec![],
            graph_group_structure: TiVec::new(),
            external_particles: vec![],
            external_signature: SignatureLike::from_iter(iter::empty::<i8>()),
            group_derived_data: TiVec::new(),
        }
    }

    pub(crate) fn add_graph(&mut self, graph: Graph) -> Result<()> {
        let new_external_particels = graph.get_external_partcles();
        let new_external_signature = graph.get_external_signature();

        if !self.graphs.is_empty() {
            if self.external_particles != new_external_particels {
                return Err(eyre!("amplitude graph has different number of externals")).with_note(
                    || {
                        format!(
                            "Found {} externals, expected {} for the graph {}",
                            new_external_particels.len(),
                            self.external_particles.len(),
                            DotGraph::from(&graph).debug_dot()
                        )
                    },
                );
            }

            if self.external_signature != new_external_signature {
                return Err(eyre!("wrong external signature"));
            }
        } else {
            self.external_particles = new_external_particels;
            self.external_signature = new_external_signature;
        }

        self.graphs.push(AmplitudeGraph::new(graph));

        //  TODO: validate that the graph is compatible
        Ok(())
    }
}

#[cfg(test)]
pub mod test {

    use linnet::half_edge::involution::{EdgeVec, Orientation};
    use symbolica::evaluate::OptimizationSettings;

    use crate::{
        dot,
        gammaloop_integrand::{GenericEvaluator, ParamBuilder},
        graph::parse::IntoGraph,
        initialisation::test_initialise,
        processes::AmplitudeGraph,
        utils::{test_utils::load_generic_model, GS},
    };
    #[test]
    fn amplitude_tree() {
        test_initialise().unwrap();
        let mut graph: AmplitudeGraph = dot!(digraph qqx_aaa_tree_1 {
                num="spenso::g(spenso::dind(spenso::cof(3, hedge(1))), spenso::cof(3, hedge(2)))/3"
                ext    [style=invis]
                ext -> v1:1 [particle="d" id=1];
                ext -> v3:2 [particle="d~" id=2];
                v1:3 -> ext [particle="a" id=3];
                v2:4 -> ext [particle="a" id=4];
                v3:0 -> ext [particle="a" id=0];
                v1 -> v2 [particle="d" id=5];
                v2 -> v3 [particle="d" id=6];
    })
    .unwrap();

        let model = load_generic_model("sm");

        graph.generate_cff().unwrap();
        // graph.build_parametric_integrand(&GenerationSettings::default());

        let param_builder = ParamBuilder::new(&graph.graph, &model);
        println!("{param_builder}");

        GenericEvaluator::new_from_builder(
            [GS.orientation_delta(&EdgeVec::from_iter(vec![Orientation::Default; 7]))],
            &param_builder,
            OptimizationSettings::default(),
        )
        .unwrap();
        // println!(" {}", a);
    }
}
