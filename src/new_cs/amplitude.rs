use std::{fs, iter, marker::PhantomData, ops::Deref, path::Path};

use ahash::AHashSet;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use momtrop::SampleGenerator;

use idenso::color::ColorSimplifier;
use spenso::{
    algebra::complex::Complex,
    network::{store::TensorScalarStoreMapping, Sequential, SmallestDegree},
    tensors::{data::StorageTensor, parametric::ParamOrConcrete},
};
use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    cff::{
        esurface::{generate_esurface_data, get_existing_esurfaces, EsurfaceDerivedData},
        expression::{
            AmplitudeOrientationID, CFFExpression, GraphOrientation, OrientationData,
            SubgraphOrientationID,
        },
        generation::{generate_cff_expression, get_orientations_from_subgraph},
    },
    model::ArcParticle,
    momentum_sample::ExternalIndex,
    new_gammaloop_integrand::{
        amplitude_integrand::{AmplitudeGraphTerm, AmplitudeIntegrand, AmplitudeIntegrandData},
        GenericEvaluator, LmbMultiChannelingSetup, ParamBuilder,
    },
    new_graph::{LMBext, LmbIndex, LoopMomentumBasis},
    numerator::symbolica_ext::AtomCoreExt,
    signature::SignatureLike,
    subtraction::{
        amplitude_counterterm::AmplitudeCountertermData,
        overlap::{find_maximal_overlap, OverlapStructure},
    },
    utils::{GS, TENSORLIB},
    uv::{approx::do_replacement_rules, UltravioletGraph},
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::{eyre, Context};
use itertools::Itertools;
use linnet::half_edge::involution::{HedgePair, Orientation};
use log::debug;
use symbolica::{atom::Atom, domains::rational::Rational, evaluate::OptimizationSettings};
use typed_index_collections::TiVec;

use crate::{
    cff::esurface::EsurfaceID,
    model::Model,
    momentum::{Rotation, RotationMethod},
    new_gammaloop_integrand::NewIntegrand,
    new_graph::{FeynmanGraph, Graph},
    DependentMomentaConstructor, GenerationSettings, RuntimeSettings,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Amplitude<S: AmplitudeState = ()> {
    pub name: String,
    pub integrand: Option<NewIntegrand>,
    pub graphs: Vec<AmplitudeGraph<S>>,
    pub external_particles: Vec<ArcParticle>,
    pub external_signature: SignatureLike<ExternalIndex>,
}

impl<S: AmplitudeState> Amplitude<S> {
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
        if let Some(integrand) = self.integrand.take() {
            integrand.save(&p, override_existing)?;
        };

        let binary = bincode::encode_to_vec(&(*self), bincode::config::standard())?;
        fs::write(p.join("amp.bin"), binary)?;
        Ok(())
    }

    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
    ) -> Result<()> {
        for amplitude_graph in self.graphs.iter_mut() {
            amplitude_graph.preprocess(model, settings)?;
        }
        Ok(())
    }

    pub(crate) fn build_integrand(
        &mut self,
        settings: RuntimeSettings,
        model: &Model,
    ) -> Result<()> {
        let terms = self
            .graphs
            .iter()
            .map(|graph| graph.generate_term_for_graph(&settings, model))
            .collect_vec();

        let rotations: Vec<Rotation> = Some(Rotation::new(RotationMethod::Identity))
            .into_iter()
            .chain(
                settings
                    .stability
                    .rotation_axis
                    .iter()
                    .map(|axis| Rotation::new(axis.rotation_method())),
            )
            .collect();

        // let orig_polarizations = self.polarizations(&settings.kinematics.externals);

        // let polarizations = rotations
        //     .iter()
        //     .map(|r| orig_polarizations.rotate(r))
        //     .collect();

        let amplitude_integrand = AmplitudeIntegrand {
            settings,
            data: AmplitudeIntegrandData {
                name: self.name.clone(),
                rotations,
                // polarizations,
                graph_terms: terms,
                external_signature: self.external_signature.clone(),
                // param_builder: self.
            },
        };
        self.integrand = Some(NewIntegrand::Amplitude(amplitude_integrand));
        Ok(())
    }

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
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct AmplitudeGraph<S: AmplitudeState = ()> {
    pub graph: Graph,
    pub derived_data: AmplitudeDerivedData<S>,
}

impl AmplitudeGraph {
    pub(crate) fn new(graph: Graph) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {
                all_mighty_integrand: Atom::Zero,
                cff_expression: None,
                state: PhantomData,
                lmbs: None,
                tropical_sampler: None,
                multi_channeling_setup: None,
                threshold_counterterm: AmplitudeCountertermData::new_empty(),
                esurface_data: None,
            },
        }
    }
}

impl<S: AmplitudeState> AmplitudeGraph<S> {
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer)
    }

    pub(crate) fn generate_cff(&mut self) -> Result<()> {
        let shift_rewrite = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_expression = generate_cff_expression(
            &self.graph.underlying,
            &shift_rewrite,
            &self.graph.underlying.dummy_list(),
        )?;
        self.derived_data.cff_expression = Some(cff_expression);

        Ok(())
    }

    //Stage 1
    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
    ) -> Result<()> {
        debug!("Generating Cff");
        self.generate_cff()?;
        debug!("Building Parametric Integrand");
        self.build_parametric_integrand();
        debug!("Building Tropical Sampler");
        self.build_tropical_sampler(settings)?;
        debug!("Building Loop Momentum Bases");
        self.build_lmbs();
        debug!("Building Multi-Channeling Channels");
        self.build_multi_channeling_channels();
        debug!("Building ESurface Derived Data");
        self.build_esurface_derived_data()?;

        if settings.enable_thresholds {
            self.build_counterterm_evaluators(model);
        }

        Ok(())
    }

    pub(crate) fn build_esurface_derived_data(&mut self) -> Result<()> {
        let lmbs = self.derived_data.lmbs.as_ref().unwrap();
        let esurfaces = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache;

        let esurface_data = generate_esurface_data(&self.graph, lmbs, esurfaces)?;
        Ok(self.derived_data.esurface_data = Some(esurface_data))
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn ct_params(&self, model: &Model) -> ParamBuilder<f64> {
        let mut param_builder = self.param_builder_core(model);
        param_builder.uv_damp_atom(vec![
            Atom::var(GS.uv_damp_plus),
            Atom::var(GS.uv_damp_minus),
        ]);
        param_builder.derivative_at_tstar_atom(Atom::var(GS.deta));
        param_builder.radius_atom(Atom::var(GS.radius));
        param_builder.radius_star_atom(Atom::var(GS.radius_star));
        param_builder
    }

    fn param_builder_core(&self, model: &Model) -> ParamBuilder<f64> {
        // the float type does not matter here
        let mut param_builder = ParamBuilder::<f64>::new();

        // this is wrong if we allow for vacuum graphs
        param_builder.external_energies_atom(&self.graph);

        // param_builder.polarizations(&self.graph);
        // param_builder.polar

        // spatial components of external momenta
        param_builder.external_spatial_atom(&self.graph);
        param_builder.polarization_params(&self.graph);
        // spatial EMR
        param_builder.emr_spatial_atom(&self.graph);

        param_builder.model_parameters_atom(model);

        param_builder.m_uv_atom(Atom::var(GS.m_uv));

        param_builder.mu_r_sq_atom(Atom::var(GS.mu_r_sq));

        self.add_function_map(&mut param_builder);

        param_builder
    }

    // pub fn fill_in_params(
    //     &self,
    //     param_builder: &mut ParamBuilder,
    //     model: &Model,
    //     settings: &Settings,
    // ) {
    //     param_builder.external_energies_value(momentum_sample);
    // }

    fn add_function_map(&self, parambuilder: &mut ParamBuilder) {
        let pi_rational = Rational::from(std::f64::consts::PI);

        for (_, e, _) in self.graph.iter_edges() {
            parambuilder
                .add_tagged_function(
                    GS.ose,
                    vec![Atom::num(e.0 as i64)],
                    format!("OSE{e}"),
                    vec![],
                    self.graph.explicit_ose_atom(e),
                )
                .unwrap()
        }
        parambuilder.add_constant(Atom::PI.into(), pi_rational.into());
        // fn_map
    }

    // fn get_eager_const_map(&self)->HashM

    fn add_additional_factors_to_cff_atom(&self, cff_atom: &Atom) -> Atom {
        // let inverse_energy_product = self.graph.underlying.get_cff_inverse_energy_product();
        let factors_of_pi =
            (Atom::var(Atom::PI) * 2).npow(3 * self.graph.underlying.get_loop_number() as i64);

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
            number_of_terms_in_epsilon_expansion: self
                .graph
                .n_loops(&self.graph.underlying.no_dummy())
                as i64
                + 1,
            // temporary_directory: Some("./form".into()),
            mu_r_sq_symbol: GS.mu_r_sq.get_name().to_string(),
            ..VakintSettings::default()
        }))
        .unwrap()
    }

    pub(crate) fn build_parametric_integrand(&mut self) {
        self.derived_data.all_mighty_integrand += self.build_original_parametric_integrand();
        for (_, a) in self
            .build_threshold_counterterm_parametric_integrand()
            .into_iter_enumerated()
        {
            // self.derived_data.all_mighty_integrand += a;
        }
    }

    fn build_threshold_counterterm_parametric_integrand(&self) -> TiVec<EsurfaceID, Atom> {
        let pols = self.graph.global_network();

        let mut counterterms: TiVec<EsurfaceID, Atom> = TiVec::new();
        let canonize_esurface = self
            .graph
            .underlying
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
                counterterms.push(Atom::new());
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
                    .collect::<TiVec<SubgraphOrientationID, _>>();

            let complement_orientations = get_orientations_from_subgraph(
                &self.graph.underlying,
                &complement,
                &reverse_dangling,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .collect::<TiVec<SubgraphOrientationID, _>>();

            let vakint = self.new_vakint();

            circled_forest.compute(
                &self.graph,
                &circled,
                &vakint,
                &circled_orientations,
                &canonize_esurface,
                &esurface.energies,
            );

            complement_forest.compute(
                &self.graph,
                &complement,
                &vakint,
                &complement_orientations,
                &canonize_esurface,
                &esurface.energies,
            );

            let circled_expr =
                circled_forest.orientation_parametric_expr(Some(&edges_in_cut), &self.graph);

            let complement_expr = complement_forest.orientation_parametric_expr(None, &self.graph);

            // println!("Circled Expression Network:");
            // println!("{}", circled_expr.dot_pretty());

            // println!("Complement Expression Network:");
            // println!("{}", complement_expr.dot_pretty());

            let mut product = circled_expr * complement_expr * pols.clone();

            product
                .execute::<Sequential, SmallestDegree, _, _>(TENSORLIB.read().unwrap().deref())
                .unwrap();
            // println!("{}", product.dot_pretty());

            let scalar: Atom = product.result_scalar().unwrap().into();

            let mut counterterm = scalar.unwrap_function(GS.color_wrap).simplify_color();

            let loop_3 = self.graph.underlying.get_loop_number() as i64 * 3;

            let grad_eta = Atom::var(GS.deta);
            let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3);

            let radius = Atom::var(GS.radius);
            let radius_star = Atom::var(GS.radius_star);
            let uv_damp_plus = Atom::var(GS.uv_damp_plus);
            let uv_damp_minus = Atom::var(GS.uv_damp_minus);

            let delta_r_plus = &radius - &radius_star;
            let delta_r_minus = -&radius - &radius_star;

            let jacobian_ratio = (&radius_star / &radius).npow(loop_3 - 1);

            let prefactor = jacobian_ratio / factors_of_pi / grad_eta
                * (uv_damp_plus / delta_r_plus + uv_damp_minus / delta_r_minus);

            counterterm *= prefactor * &counterterm;
            // println!("CounterTerm{}", counterterm);
            counterterms.push(counterterm);
        }

        counterterms
    }

    fn build_original_parametric_integrand(&self) -> Atom {
        let wood = self.graph.wood(&self.graph.underlying.no_dummy());
        let mut forest = wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

        let canonize_esurface = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let orientations: TiVec<AmplitudeOrientationID, OrientationData> = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|a| a.data.clone())
            .collect();

        let vakint = self.new_vakint();

        forest.compute(
            &self.graph,
            &self.graph.underlying.no_dummy(),
            &vakint,
            &orientations,
            &canonize_esurface,
            &[],
        );

        let mut pols = self.graph.global_network();

        let mut reps = Vec::new();
        for (p, eid, _) in self.graph.as_ref().iter_edges() {
            if p.is_paired() {
                reps.push(GS.add_parametric_sign(eid));
            }
        }

        pols = pols.replace_multiple(&reps);

        let mut full = forest.orientation_parametric_expr(None, &self.graph);

        full *= pols;

        full.execute::<Sequential, SmallestDegree, _, _>(TENSORLIB.read().unwrap().deref())
            .unwrap();

        let mut scalar: Atom = full
            .result_scalar()
            .expect(&format!(
                "Failed to get scalar from network:{}",
                full.dot_pretty()
            ))
            .into();

        scalar = scalar.unwrap_function(GS.color_wrap).simplify_color();

        scalar = self.add_additional_factors_to_cff_atom(&scalar);

        debug!("All parametric integrand atom:{:>}", scalar);

        scalar
    }

    fn build_counterterm_evaluators(&mut self, model: &Model) {
        let pols = self.graph.global_network();

        // println!("Pols:{}", pols.dot_pretty());
        let mut counterterms: TiVec<EsurfaceID, Atom> = TiVec::new();
        let canonize_esurface = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let full_orientation_list = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|o| o.data.clone())
            .collect::<TiVec<AmplitudeOrientationID, _>>();

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
                counterterms.push(Atom::new());
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
                    .collect::<TiVec<SubgraphOrientationID, _>>();

            let complement_orientations = get_orientations_from_subgraph(
                &self.graph.underlying,
                &complement,
                &reverse_dangling,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .collect::<TiVec<SubgraphOrientationID, _>>();

            let vakint = self.new_vakint();

            circled_forest.compute(
                &self.graph,
                &circled,
                &vakint,
                &circled_orientations,
                &canonize_esurface,
                &esurface.energies,
            );

            complement_forest.compute(
                &self.graph,
                &complement,
                &vakint,
                &complement_orientations,
                &canonize_esurface,
                &esurface.energies,
            );

            let mut counterterm = Atom::new();
            for orientation in orientations {
                let circled_expr = circled_forest.local_expr(
                    &full_orientation_list[orientation],
                    Some(&edges_in_cut),
                    &self.graph,
                );

                let complement_expr = complement_forest.local_expr(
                    &full_orientation_list[orientation],
                    None,
                    &self.graph,
                );

                // println!("Circled Expression Network:");
                // println!("{}", circled_expr.dot_pretty());

                // println!("Complement Expression Network:");
                // println!("{}", complement_expr.dot_pretty());

                let mut product = circled_expr
                    * complement_expr
                    * pols.map_ref(
                        |a| full_orientation_list[orientation].select(a),
                        |a| match a {
                            ParamOrConcrete::Param(a) => {
                                ParamOrConcrete::Param(a.map_data_ref_self(|a| {
                                    full_orientation_list[orientation].select(a)
                                }))
                            }
                            a => a.clone(),
                        },
                    );

                product
                    .execute::<Sequential, SmallestDegree, _, _>(TENSORLIB.read().unwrap().deref())
                    .unwrap();
                // println!("{}", product.dot_pretty());

                let mut scalar: Atom = product.result_scalar().unwrap().into();

                scalar = scalar.unwrap_function(GS.color_wrap).simplify_color();
                let orientation_result = do_replacement_rules(
                    scalar,
                    &self.graph,
                    &esurface.energies,
                    &full_orientation_list[orientation].orientation,
                );

                counterterm += orientation_result;
            }
            // println!("CounterTerm{}", counterterm);

            let loop_3 = self.graph.underlying.get_loop_number() as i64 * 3;

            let grad_eta = Atom::var(GS.deta);
            let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3);

            let radius = Atom::var(GS.radius);
            let radius_star = Atom::var(GS.radius_star);
            let uv_damp_plus = Atom::var(GS.uv_damp_plus);
            let uv_damp_minus = Atom::var(GS.uv_damp_minus);

            let delta_r_plus = &radius - &radius_star;
            let delta_r_minus = -&radius - &radius_star;

            let jacobian_ratio = (&radius_star / &radius).npow(loop_3 - 1);

            let prefactor = jacobian_ratio / factors_of_pi / grad_eta
                * (uv_damp_plus / delta_r_plus + uv_damp_minus / delta_r_minus);

            counterterm = prefactor * &counterterm;
            // println!("CounterTerm{}", counterterm);
            counterterms.push(counterterm);
        }

        let params = self.ct_params(model);
        let counterterm_evaluators = counterterms
            .into_iter()
            .map(|ct| {
                GenericEvaluator::new_from_builder([ct], &params, OptimizationSettings::default())
                    .unwrap()
            })
            .collect();

        let threshold_counterterm = AmplitudeCountertermData {
            overlap: OverlapStructure::new_empty(),
            evaluators: counterterm_evaluators,
            param_builder: params,
        };

        self.derived_data.threshold_counterterm = threshold_counterterm;
    }

    pub(crate) fn build_all_orientation_integrand_evaluator(
        &self,
        param_builder: &ParamBuilder,
    ) -> GenericEvaluator {
        debug!("Building all orientation integrand_evaluator");
        GenericEvaluator::new_from_builder(
            [self
                .orientation_atoms()
                .into_iter()
                .fold(Atom::Zero, |acc, a| acc + a)],
            &param_builder,
            OptimizationSettings::default(),
        )
        .unwrap()
    }

    pub(crate) fn build_orientations_evaluators(
        &self,
        param_builder: &ParamBuilder,
    ) -> TiVec<AmplitudeOrientationID, GenericEvaluator> {
        self.orientation_atoms()
            .into_iter()
            .map(|a| {
                GenericEvaluator::new_from_builder(
                    [a],
                    &param_builder,
                    OptimizationSettings::default(),
                )
                .unwrap()
            })
            .collect()
    }

    fn orientation_atoms(&self) -> TiVec<AmplitudeOrientationID, Atom> {
        self.derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|o| o.data.select(&self.derived_data.all_mighty_integrand))
            .collect()
    }

    fn build_lmbs(&mut self) {
        let lmbs = self
            .graph
            .generate_loop_momentum_bases(&self.graph.underlying.no_dummy());

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_tropical_sampler(&mut self, process_settings: &GenerationSettings) -> Result<()> {
        if process_settings
            .tropical_subgraph_table_settings
            .disable_tropical_generation
        {
            debug!("Tropical subgraph table generation is disabled.");
            return Ok(());
        }
        let num_virtual_loop_edges = self.graph.iter_loop_edges().count();
        let num_loops = self.graph.loop_momentum_basis.loop_edges.len();
        let target_omega = process_settings
            .tropical_subgraph_table_settings
            .target_omega;

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
    fn generate_term_for_graph(
        &self,
        settings: &RuntimeSettings,
        model: &Model,
    ) -> AmplitudeGraphTerm {
        let estimated_scale = self
            .graph
            .underlying
            .expected_scale(settings.kinematics.e_cm);

        let esurfaces = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache;

        let esurface_data = self.derived_data.esurface_data.as_ref().unwrap();
        let externals = settings.kinematics.externals.get_dependent_externals(
            DependentMomentaConstructor::Amplitude(&self.graph.get_external_signature()),
        );

        let existing_esurfaces = get_existing_esurfaces(
            esurfaces,
            esurface_data,
            &externals,
            &self.graph.loop_momentum_basis,
            settings.general.debug,
            settings.kinematics.e_cm,
        );

        let mut param_builder = self.param_builder_core(model);
        param_builder.add_external_four_mom(&externals);
        param_builder.polarizations_values(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );
        param_builder.model_parameters_value(model);
        param_builder.mu_r_sq_value(Complex::new_zero());
        param_builder.m_uv_value(Complex::new_zero());

        // (momentum_sample);

        let edge_masses = self.graph.new_edgevec(|edge, _, _| edge.mass_value());

        let overlap = find_maximal_overlap(
            &self.graph.loop_momentum_basis,
            &existing_esurfaces,
            esurfaces,
            &edge_masses,
            &externals,
            &settings,
        );

        let mut threshold_counterterm = self.derived_data.threshold_counterterm.clone();

        threshold_counterterm.overlap = overlap;

        threshold_counterterm
            .param_builder
            .add_external_four_mom(&externals);
        threshold_counterterm.param_builder.polarizations_values(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );
        threshold_counterterm
            .param_builder
            .model_parameters_value(model);
        threshold_counterterm
            .param_builder
            .mu_r_sq_value(Complex::new_zero());
        threshold_counterterm
            .param_builder
            .m_uv_value(Complex::new_zero());

        AmplitudeGraphTerm {
            integrand_evaluator_all_orientations: self
                .build_all_orientation_integrand_evaluator(&param_builder),
            integrand_evaluators: self.build_orientations_evaluators(&param_builder),
            tropical_sampler: self.derived_data.tropical_sampler.clone(),
            graph: self.graph.clone(),
            multi_channeling_setup: self
                .derived_data
                .multi_channeling_setup
                .clone()
                .expect("multi_channeling_setup should have been created"),
            lmbs: self
                .derived_data
                .lmbs
                .clone()
                .expect("lmbs should have been created"),
            threshold_counterterm,
            estimated_scale,
            esurfaces: self
                .derived_data
                .cff_expression
                .as_ref()
                .expect("cff_expression should have been created")
                .surfaces
                .esurface_cache
                .clone(),

            param_builder,
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeDerivedData<S: AmplitudeState> {
    pub all_mighty_integrand: Atom,
    pub threshold_counterterm: AmplitudeCountertermData,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub cff_expression: Option<CFFExpression<AmplitudeOrientationID>>,
    pub esurface_data: Option<EsurfaceDerivedData>,
    pub state: PhantomData<S>,
}

pub trait AmplitudeState:
    Clone + std::fmt::Debug + bincode::Encode + for<'a> bincode::Decode<GammaLoopContextContainer<'a>>
{
}
impl AmplitudeState for () {}

#[derive(Clone, Encode, Decode, Debug)]
pub struct ReadyForTerm {}
impl AmplitudeState for ReadyForTerm {}

// #[derive(Clone, Encode, Decode, Debug)]
// pub struct ReadyForTerm {}
// impl AmplitudeState for ReadyForTerm {}

impl Amplitude {
    pub(crate) fn from_dot_string<Str: AsRef<str>>(
        s: Str,
        name: String,
        model: &Model,
    ) -> Result<Self> {
        let graphs = Graph::from_string(s, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub(crate) fn from_dot_file<'a, P>(p: P, name: String, model: &Model) -> Result<Self>
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

    pub(crate) fn new(name: impl ToString) -> Self {
        Self {
            integrand: None,
            name: name.to_string(),
            graphs: vec![],
            external_particles: vec![],
            external_signature: SignatureLike::from_iter(iter::empty::<i8>()),
        }
    }

    pub(crate) fn add_graph(&mut self, graph: Graph) -> Result<()> {
        let new_external_particels = graph.get_external_partcles();
        let new_external_signature = graph.get_external_signature();

        if !self.graphs.is_empty() {
            if self.external_particles != new_external_particels {
                return Err(eyre!("amplitude graph has different number of externals"));
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
