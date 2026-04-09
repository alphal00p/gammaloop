use std::{
    collections::{BTreeMap, HashMap},
    fmt,
    fs::{self, File},
    io::Write,
    iter,
    path::Path,
};

use ahash::AHashSet;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use momtrop::SampleGenerator;

use idenso::gamma::GammaSimplifier;
use rayon::{
    ThreadPool,
    iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::algebra::complex::Complex;
use tracing::{info_span, instrument};
use tracing_indicatif::{span_ext::IndicatifSpanExt, style::ProgressStyle};
use vakint::{EvaluationMethod, NumericalEvaluationResult, Vakint, vakint_symbol};

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::GroupEsurfaceId,
        expression::{CFFExpression, OrientationID},
    },
    graph::{
        GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex, LoopMomentumBasis,
        cuts::{CutSet, ResidueSelector},
    },
    integrands::process::{
        LmbMultiChannelingSetup,
        amplitude::{AmplitudeGraphTerm, AmplitudeIntegrand, AmplitudeIntegrandData},
        graph_to_group_id_for_group_structure,
    },
    model::ArcParticle,
    momentum::{sample::ExternalIndex, signature::SignatureLike},
    processes::{
        DotExportSettings, GraphGenerationStats, NamedGraphGenerationReport,
        StandaloneExportSettings,
    },
    settings::{GlobalSettings, RuntimeSettings, runtime::LockedRuntimeSettings},
    subtraction::amplitude_counterterm::AmplitudeCountertermAtom,
    utils::{F, GS, Length, W_},
    uv::{
        RenormalizationPart, UVgenerationSettings, UltravioletGraph,
        approx::{CutStructure, integrated::to_vakint_integrand},
        hedge_poset::Wood as NewWood,
        settings::VakintSettings,
        wood::CutWoods,
    },
};
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::{
    half_edge::{
        involution::{EdgeVec, HedgePair},
        subgraph::{ModifySubSet, SuBitGraph, SubGraphLike, SubSetOps},
    },
    parser::DotGraph,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol, Var},
    function,
};
use tracing::{debug, info};
use typed_index_collections::{TiVec, ti_vec};

use crate::{
    cff::esurface::EsurfaceID,
    graph::{FeynmanGraph, Graph},
    integrands::process::ProcessIntegrand,
    model::Model,
    settings::global::GenerationSettings,
};

use crate::graph::parse::complete_group_parsing;

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Amplitude {
    pub name: String,
    pub integrand: Option<ProcessIntegrand>,
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
    pub fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        if let Some(integrand) = &self.integrand {
            integrand.export_standalone(path, settings)?
        } else {
            return Err(eyre!(
                "Cannot warm up amplitude {} without integrand",
                self.name
            ));
        }

        Ok(())
    }

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
            amp.integrand = Some(ProcessIntegrand::Amplitude(integrand));
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
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        info!("Compiling amplitude {}", self.name);
        let p = path.as_ref().join(&self.name);

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
            let compile_times = integrand.compile(&p, override_existing, thread_pool)?;
            return Ok(compile_times
                .into_iter()
                .map(|(graph_name, duration)| NamedGraphGenerationReport {
                    integrand_name: self.name.clone(),
                    graph_name,
                    stats: GraphGenerationStats {
                        total_time: duration,
                        evaluator_compile_time: duration,
                        ..GraphGenerationStats::default()
                    },
                })
                .collect());
        };
        Ok(Vec::new())
    }

    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
              path = %path.as_ref().display(),
          )
    )]
    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let p = path.as_ref().join(&self.name);

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
            file.write_all(&binary)?;
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
        locked_runtime_settings: &LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        // preprocess each graph individually
        let integrand_name = self.name.clone();

        let preprocess_span = info_span!("Preprocessing graphs", indicatif.pb_show = true);
        preprocess_span.pb_set_style(&ProgressStyle::with_template(
            "{wide_bar} {pos}/{len} {msg}",
        )?);
        preprocess_span.pb_set_length(self.graphs.len() as u64);
        preprocess_span.pb_set_message("Preprocessing graphs");

        let preprocess_span_enter = preprocess_span.enter();

        let preprocess_reports = thread_pool.install(|| {
            let parent = preprocess_span.clone();
            self.graphs
                .par_iter_mut()
                .map(|amplitude_graph| {
                    let _guard = parent.enter();
                    let stats =
                        amplitude_graph.preprocess(model, settings, locked_runtime_settings);
                    preprocess_span.pb_inc(1);

                    stats.map(|stats| NamedGraphGenerationReport {
                        integrand_name: integrand_name.clone(),
                        graph_name: amplitude_graph.graph.name.clone(),
                        stats,
                    })
                })
                .collect::<Result<Vec<_>>>()
        })?;

        drop(preprocess_span_enter);
        drop(preprocess_span);

        self.generate_grouped_derived_data()?;

        Ok(preprocess_reports)
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name, indicatif.pb_show = true, indicatif.pb_msg = "Generating Evaluators",
          )
      )]
    pub fn build_integrand(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        let integrand_name = self.name.clone();
        let mut graph_reports = Vec::new();
        let terms: Vec<_> = thread_pool.install(|| {
            self.graphs
                .par_iter_mut()
                .enumerate()
                .map(|(graph_id, graph)| {
                    let graph_started = std::time::Instant::now();
                    let group_id = graph.graph.group_id.unwrap(); // should always be set
                    let esurface_map = &self.group_derived_data[group_id].esurface_map;
                    let group_pos = self.graph_group_structure[group_id]
                        .find_position(graph_id)
                        .unwrap();

                    let (term, mut stats) = graph.generate_term_for_graph(
                        model,
                        group_pos,
                        esurface_map.clone(),
                        global_settings,
                    )?;
                    stats.evaluator_count = term.generic_evaluator_count();
                    stats.total_time += graph_started.elapsed();
                    Ok((
                        term,
                        NamedGraphGenerationReport {
                            integrand_name: integrand_name.clone(),
                            graph_name: graph.graph.name.clone(),
                            stats,
                        },
                    ))
                })
                .collect::<Result<Vec<_>>>()
        })?;
        for (_, report) in &terms {
            graph_reports.push(report.clone());
        }
        let mut terms = terms.into_iter().map(|(term, _)| term).collect::<Vec<_>>();

        for group in self.graph_group_structure.iter() {
            let master_graph_id = group.master();
            let mc_of_master = self.graphs[master_graph_id]
                .derived_data
                .multi_channeling_setup
                .as_ref()
                .unwrap();
            let master_graph = &self.graphs[master_graph_id].graph;
            let master_external_signature = master_graph.get_external_signature();
            let master_external_pdgs = master_graph
                .get_external_partcles()
                .into_iter()
                .map(|particle| particle.pdg_code)
                .collect_vec();

            for graph_id in group.into_iter() {
                terms[graph_id].multi_channeling_setup = mc_of_master.clone();
                terms[graph_id].master_external_signature = master_external_signature.clone();
                terms[graph_id].master_external_pdgs = master_external_pdgs.clone();
            }
        }

        let mut amplitude_integrand = AmplitudeIntegrand {
            settings: runtime_default.into_with_modified_kinematics(
                &self.external_signature,
                &self.graphs[0].graph.get_external_masses(model),
            )?,
            data: AmplitudeIntegrandData {
                name: self.name.clone(),
                compilation: global_settings
                    .generation
                    .compile
                    .frozen_mode(&global_settings.generation.evaluator),
                rotations: None,
                loop_cache_id: 0,
                external_cache_id: 0,
                base_external_cache_id: 0,
                graph_terms: terms,
                external_signature: self.external_signature.clone(),
                graph_group_structure: self.graph_group_structure.clone(),
                graph_to_group_id: graph_to_group_id_for_group_structure(
                    &self.graph_group_structure,
                ),
                group_derived_data: self.group_derived_data.clone(),
            },
            event_processing_runtime: Default::default(),
            active_f64_backend: Default::default(),
        };
        let compile_times =
            amplitude_integrand.prepare_runtime_backends_after_generation_with_compile_times()?;
        for (report, compile_time) in graph_reports.iter_mut().zip(compile_times) {
            report.stats.evaluator_compile_time += compile_time;
            report.stats.total_time += compile_time;
        }
        self.integrand = Some(ProcessIntegrand::Amplitude(amplitude_integrand));
        Ok(graph_reports)
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        for graph in &self.graphs {
            graph.write_dot(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    pub fn write_dot_fmt<W: fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        for graph in &self.graphs {
            graph.write_dot_fmt(writer, settings)?;
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

        let _: () = self.group_derived_data = group_derived_data;
        Ok(())
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
    pub fn renormalization_part(
        &mut self,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        if self.derived_data.cff_expression.is_none() {
            self.generate_cff()?;
        }
        let valid_orientations: Vec<_> = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        if settings.use_legacy {
            let mut vk_settings = settings.vakint.true_settings();
            let wood = self.graph.wood_with_settings(
                &self.graph.no_dummy(),
                settings,
                &self.graph.loop_momentum_basis,
            );
            //  it needs to be the max number of loops across all divergent spinneys of that graph
            vk_settings.number_of_terms_in_epsilon_expansion = wood.max_loops as i64;

            let mut forest = wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

            let vk = (crate::utils::vakint()?, &vk_settings);
            let cuts = CutSet::empty(self.graph.n_hedges());
            forest.compute(&mut self.graph, vk, &cuts, &valid_orientations, settings)?;

            forest.pole_part_of_ends(&self.graph)
        } else {
            let cuts = CutStructure::empty(&self.graph);
            let wood = NewWood::new(cuts, &self.graph, settings);
            let mut forest = wood.unfold();
            forest.integrate(&self.graph, crate::utils::vakint()?, settings)?;

            forest.pole_part_of_ends(&self.graph)
        }
    }

    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer, settings)
    }

    pub(crate) fn write_dot_fmt<W: fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        self.graph.dot_serialize_fmt(writer, settings)
    }

    #[instrument(skip_all, fields(indicatif.pb_show = true,indicatif.pb_msg = "Generating CFF"), err)]
    pub(crate) fn generate_cff(&mut self) -> Result<()> {
        let shift_rewrite = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_expression = self.graph.generate_cff(&[], &shift_rewrite)?;
        self.derived_data.cff_expression = Some(cff_expression);

        Ok(())
    }
    #[instrument(skip_all, fields(indicatif.pb_show = true,indicatif.pb_msg = "preprocessing"), err)]
    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
    ) -> Result<GraphGenerationStats> {
        let preprocess_started = std::time::Instant::now();
        let vk = crate::utils::vakint()?;

        self.generate_cff()?;

        self.build_integrands(settings, vk)?;

        if self.graph.is_group_master {
            self.build_tropical_sampler(settings)?;
        }

        self.build_lmbs();

        if self.graph.is_group_master {
            self.build_multi_channeling_channels(settings.override_lmb_heuristics);
        }

        if settings.threshold_subtraction.enable_thresholds {
            self.derived_data.threshold_counterterms = self
                .build_threshold_counterterm_parametric_integrand(
                    settings,
                    vk,
                    locked_runtime_settings,
                    model,
                )?;
        }

        Ok(GraphGenerationStats {
            total_time: preprocess_started.elapsed(),
            ..GraphGenerationStats::default()
        })
    }

    #[instrument(skip_all, fields(indicatif.pb_show = true,indicatif.pb_msg = "Building Multi-Channeling Channels"))]
    fn build_multi_channeling_channels(&mut self, override_lmb_heuristics: bool) {
        let channels = self.graph.build_multi_channeling_channels(
            self.derived_data.lmbs.as_ref().unwrap(),
            override_lmb_heuristics,
        );

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
        let factors_of_pi = (Atom::var(GS.pi) * 2).pow(3 * self.graph.get_loop_number() as i64);

        // debug!("result: {}", result);
        cff_atom / factors_of_pi
    }

    pub fn to_numerical(
        numerical_result: AtomView,
        true_settings: &vakint::VakintSettings,
    ) -> Result<NumericalEvaluationResult> {
        Ok(NumericalEvaluationResult::from_atom(
            numerical_result,
            vakint_symbol!(&true_settings.epsilon_symbol),
            true_settings,
        )?)
    }

    pub fn analytical_evaluation<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(
        &self,
        model: &Model,
        refresh_model_values: bool,
        component: &S,
        evaluate_numerically: bool,
        vakint: &Vakint,
        true_settings: &vakint::VakintSettings,
        settings: &VakintSettings,
        run_time_settings: &RuntimeSettings,
        include_global_numerator: bool,
    ) -> Result<Atom> {
        let mut true_settings = true_settings.clone();
        true_settings.number_of_terms_in_epsilon_expansion =
            self.graph.n_loops(&self.graph.no_dummy()) as i64 + 1;
        let pysec_dec_enabled_in_vakint = true_settings.evaluation_order.0.iter().find_map(|o| {
            if let EvaluationMethod::PySecDec(opts) = o {
                Some(opts)
            } else {
                None
            }
        });

        let complex_params_vakint = if evaluate_numerically || pysec_dec_enabled_in_vakint.is_some()
        {
            let mut param_builder = self.graph.param_builder.clone(); //ParamBuilder::<f64>::new(&self.graph, model);
            if refresh_model_values {
                param_builder.update_model_values(model);
            }
            param_builder.m_uv_value(Complex::new_re(F(run_time_settings.general.m_uv)));
            param_builder.mu_r_sq_value(Complex::new_re(F(run_time_settings.general.mu_r_sq())));

            // println!("\nParamBuilder parameters:\n{}", param_builder);

            let mut complex_params: HashMap<String, symbolica::domains::float::Complex<f64>> =
                HashMap::default();
            for params in param_builder.pairs.into_iter() {
                let ps: crate::integrands::process::ParamValuePairs = params;
                for (p_name, p_value) in
                    ps.params
                        .iter()
                        .zip(ps.value_range)
                        .map(|(a, value_index)| {
                            (
                                a.to_canonical_string(),
                                param_builder.values[0][value_index],
                            )
                        })
                {
                    complex_params.insert(
                        p_name,
                        symbolica::domains::float::Complex::new(
                            p_value.re.into(),
                            p_value.im.into(),
                        ),
                    );
                }
            }
            // Make sure to remove entries already supported by vakint, as they may not match required precision
            for atom in &[
                Atom::Var(Var::new(Symbol::PI)),
                Atom::Var(Var::new(vakint_symbol!("EulerGamma"))),
                function!(Symbol::LOG, Atom::num(2)),
            ] {
                _ = complex_params.remove(&atom.to_string());
            }

            // Make sure to properly do the upcasting to required precision in vakint settings
            vakint.params_from_complex_f64(&true_settings, &complex_params)
        } else {
            HashMap::default()
        };

        if let Some(pysec_dec_opts) = pysec_dec_enabled_in_vakint {
            true_settings.evaluation_order.adjust(
                None,
                pysec_dec_opts.relative_precision,
                &HashMap::default(),
                &complex_params_vakint,
                &HashMap::default(),
            );
        }

        let mut num = self
            .graph
            .numerator(component, &self.graph.empty_subgraph());
        if include_global_numerator {
            num.state.expr *= &self.graph.global_prefactor.num;
        }

        let mut four_dimensional_integrand = num
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap()
            .simplify_gamma()
            / self
                .graph
                .denominator(component, |e| e.extra_data.vakint_edge_power.unwrap_or(1));

        // println!("Four-dimensional integrand: {}", four_dimensional_integrand);

        let mom_reps = self.graph.uv_wrapped_replacement(
            &self.graph.full_filter(),
            &self.graph.lmb_of(component),
            &[W_.x___],
        );

        // println!("Reps:");
        // for r in &mom_reps {
        //     println!("{r}");
        // }

        // rewrite the inner_t as well
        four_dimensional_integrand = four_dimensional_integrand.replace_multiple(&mom_reps);

        // println!("LMB: {}", lmb);
        // let vk_mom = vakint_symbol!("k");
        // for (i, l) in lmb.loop_edges.iter().enumerate() {
        //     four_dimensional_integrand = four_dimensional_integrand
        //         .replace(function!(GS.emr_mom, usize::from(*l) as i64))
        //         .with(function!(vk_mom, i as i64 + 1))
        //         .replace(function!(
        //             GS.emr_mom,
        //             function!(vk_mom, i as i64 + 1),
        //             W_.x___
        //         ))
        //         .with(function!(vk_mom, i as i64 + 1, W_.x___));
        // }

        let mut vakint_integrand = to_vakint_integrand(
            &four_dimensional_integrand,
            &self.graph,
            &self.graph.full_filter(),
            &self.graph.empty_subgraph::<SuBitGraph>(),
            settings,
            false,
        );

        vakint_integrand.canonicalize(&true_settings, &vakint.topologies, false)?;
        // println!("Canonized: {}", vakint_integrand);
        vakint_integrand.tensor_reduce(vakint, &true_settings)?;
        // println!("Tensor Reduced {}", vakint_integrand);
        vakint_integrand.evaluate_integral(vakint, &true_settings)?;
        // println!("Evaluated {}", vakint_integrand);
        let analytical_evaluation: Atom = vakint_integrand.into();
        // println!(
        //     "\nVakint analytical evaluation:\n{:#}",
        //     analytical_evaluation
        // );
        if !evaluate_numerically {
            Ok(analytical_evaluation)
        } else {
            let (numerical_evaluation, _error) = vakint
                .numerical_evaluation(
                    &true_settings,
                    analytical_evaluation.as_view(),
                    &HashMap::default(),
                    &complex_params_vakint,
                    None,
                )
                .unwrap();

            // println!("\nVakint numerical evaluation:\n{:#}", numerical_evaluation);

            let numerical_evaluation_atom =
                numerical_evaluation.to_atom(vakint_symbol!(true_settings.epsilon_symbol.clone()));

            Ok(numerical_evaluation_atom)
        }
    }

    #[instrument(
        skip_all,
        fields(indicatif.pb_show = true,indicatif.pb_msg = "Building Parametric Integrand"),
        err
    )]
    pub(crate) fn build_integrands(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        let valid_orientations: Vec<_> = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();
        let cutstructure = CutStructure::empty(&self.graph);
        let woods = CutWoods::new(cutstructure, &self.graph, &settings.uv);
        let mut forests = woods.unfold(&self.graph);
        forests.compute(&mut self.graph, vakint, &valid_orientations, &settings.uv)?;
        let exprs: Vec<_> = forests
            .orientation_parametric_exprs(&self.graph, false)?
            .into_iter()
            .map(|e| e.map(|a| self.add_additional_factors_to_cff_atom(&a)))
            .collect();

        self.derived_data.all_mighty_integrand = exprs[0].integrands[0].clone();
        Ok(())
    }

    #[instrument(
        skip_all,
        fields(indicatif.pb_show = true,indicatif.pb_msg = "Building Threshold Counterterms"),
        err
    )]
    fn build_threshold_counterterm_parametric_integrand(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
        locked_runtime_settings: &LockedRuntimeSettings,
        model: &Model,
    ) -> Result<TiVec<EsurfaceID, AmplitudeCountertermAtom>> {
        let mut counterterms: TiVec<EsurfaceID, AmplitudeCountertermAtom> = ti_vec![
            AmplitudeCountertermAtom::new();
            self.graph.surface_cache.esurface_cache.len()
        ];
        let valid_orientations: Vec<_> = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        let (local_prefactor, integrated_prefactor) = self.th_counterterm_prefactor_helpter();

        let global_cff = self.derived_data.cff_expression.as_ref().unwrap(); // should always be set at this point
        let esurface_raising = self
            .graph
            .determine_raised_esurfaces_from_expression(global_cff);

        let mut cuts = vec![];

        for raised_data in esurface_raising.raised_groups.into_iter() {
            assert!(
                raised_data.esurface_ids.len() == 1,
                "raise threshold subtraction not yet implemented"
            );

            assert!(
                raised_data.max_occurence == 1,
                "raise threshold subtraction not yet implemented"
            );

            let esurface_id = raised_data.esurface_ids[0];
            let esurface = &global_cff.surfaces.esurface_cache[esurface_id];

            if settings.threshold_subtraction.check_esurface_at_generation {
                let masses: EdgeVec<F<f64>> = self.graph.get_real_mass_vector(model);
                let lmb = &self.graph.loop_momentum_basis;
                if !locked_runtime_settings.existence_check(
                    esurface,
                    &masses,
                    &self.graph.get_external_signature(),
                    lmb,
                ) {
                    continue;
                }
            }

            let mut cut_union: SuBitGraph = self.graph.empty_subgraph();

            for energy in esurface.energies.iter() {
                let (_, hedge_pair) = self.graph[energy];
                match hedge_pair {
                    HedgePair::Paired { source, sink } => {
                        cut_union.add(source);
                        cut_union.add(sink);
                    }
                    _ => unreachable!(),
                }
            }

            let cutset = CutSet {
                residue_selector: ResidueSelector {
                    lu_cut: None,
                    left_th_cut: Some(raised_data.clone()),
                    right_th_cut: None,
                },
                union: cut_union,
            };

            cuts.push(cutset);
        }

        let cut_structure = CutStructure { cuts };

        let woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        let mut forests = woods.unfold(&self.graph);
        forests.compute(&mut self.graph, vakint, &valid_orientations, &settings.uv)?;

        let exprs: Vec<_> = forests.orientation_parametric_exprs(&self.graph, false)?;

        for mut expr in exprs.into_iter() {
            assert!(
                expr.integrands.len() == 1,
                "threshold counterterm generation not yet implemented for raised"
            );
            let integrand = expr.integrands.pop().unwrap();
            let local = &local_prefactor * &integrand;
            let integrated = &integrated_prefactor * integrand;

            let counterterm_atom = AmplitudeCountertermAtom {
                parametric_local: local,
                parametric_integrated: integrated,
            };
            let esurfaces = expr.cuts.residue_selector.left_th_cut.unwrap();
            let esurface_id = esurfaces.esurface_ids[0];

            counterterms[esurface_id] = counterterm_atom;
        }

        Ok(counterterms)
    }

    fn th_counterterm_prefactor_helpter(&self) -> (Atom, Atom) {
        let loop_3 = self.graph.get_loop_number() as i64 * 3;

        let grad_eta = Atom::var(GS.deta_left_th);
        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).pow(loop_3);
        let i = Atom::i();

        let radius = Atom::var(GS.radius_left);
        let radius_star = Atom::var(GS.radius_star_left);
        let uv_damp_plus = Atom::var(GS.uv_damp_plus_left);
        let uv_damp_minus = Atom::var(GS.uv_damp_minus_left);
        let hfunction = Atom::var(GS.hfunction_left_th);

        let delta_r_plus = &radius - &radius_star;
        let delta_r_minus = -&radius - &radius_star;

        let jacobian_ratio = (&radius_star / &radius).pow(loop_3 - 1);

        let local_prefactor = &jacobian_ratio / &factors_of_pi / &grad_eta
            * (uv_damp_plus / delta_r_plus + uv_damp_minus / delta_r_minus);

        let integrated_prefactor =
            -i * Atom::var(GS.pi) * &jacobian_ratio * hfunction / factors_of_pi / grad_eta;

        (local_prefactor, integrated_prefactor)
    }

    #[instrument(skip_all, fields(indicatif.pb_show = true,indicatif.pb_msg = "Building Loop Momentum Bases"))]
    fn build_lmbs(&mut self) {
        let lmbs = self
            .graph
            .generate_loop_momentum_bases_of(&self.graph.no_dummy());

        self.derived_data.lmbs = Some(lmbs)
    }

    #[instrument(skip_all, fields(indicatif.pb_show = true,indicatif.pb_msg = "Building Tropical Sampler"), err)]
    fn build_tropical_sampler(&mut self, process_settings: &GenerationSettings) -> Result<()> {
        if process_settings
            .tropical_subgraph_table
            .disable_tropical_generation
        {
            debug!("Tropical subgraph table generation is disabled.");
            return Ok(());
        }
        let num_virtual_loop_edges = self.graph.iter_loop_edges().count();

        if num_virtual_loop_edges == 0 {
            debug!("Graph has no loop edges, skipping tropical sampler generation.");
            return Ok(());
        }

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

        let _: () = self.derived_data.tropical_sampler = Some(sampler);
        Ok(())
    }

    // Expects cff_expression, esurface_data,
    #[instrument(
          name = "generate_term_for_graph",
          level = "info",
          skip(self, model, global_settings),
          fields(
              graph.name = %self.graph.name, indicatif.pb_show = true, indicatif.pb_msg = format!("Generating Evaluators for {}", self.graph.name)

          ),
          err
      )]
    fn generate_term_for_graph(
        &self,
        model: &Model,
        own_group_position: GraphGroupPosition,
        esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<EsurfaceID>>>,
        global_settings: &GlobalSettings,
    ) -> Result<(AmplitudeGraphTerm, GraphGenerationStats)> {
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
    pub cff_expression: Option<CFFExpression<OrientationID>>,
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

    pub fn from_dot_file<P>(p: P, name: String, model: &Model) -> Result<Self>
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

    fn add_graph(&mut self, graph: Graph) -> Result<()> {
        let new_external_particels = graph.get_external_partcles();
        let new_external_signature = graph.get_external_signature();

        if !self.graphs.is_empty() {
            if self.external_particles != new_external_particels {
                return Err(eyre!("amplitude graph has different externals")).with_context(|| {
                    format!(
                        "Found {} externals, expected {} for the graph {}",
                        new_external_particels.len(),
                        self.external_particles.len(),
                        DotGraph::from(&graph).debug_dot()
                    )
                });
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

    use crate::{
        dot, graph::parse::IntoGraph, initialisation::test_initialise, processes::AmplitudeGraph,
        utils::load_generic_model,
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

        let _model = load_generic_model("sm");

        graph.generate_cff().unwrap();
        // graph.build_parametric_integrand(&GenerationSettings::default());

        let param_builder = &graph.graph.param_builder;
        println!("{param_builder}");

        // GenericEvaluator::new_from_builder(
        //     [GS.orientation_delta(&EdgeVec::from_iter(vec![Orientation::Default; 7]))],
        //     param_builder,
        //     None,
        //     OptimizationSettings::default(),
        //     true,
        // )
        // .unwrap();
        // println!(" {}", a);
    }
}
