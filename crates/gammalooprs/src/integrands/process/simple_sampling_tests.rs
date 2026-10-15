use super::{
    EvaluationSource, EvaluationTarget, GammaLoopSample, GraphTerm, MomentumSpaceEvaluationInput,
    ProcessIntegrand, ProcessIntegrandImpl, SamplingChannelId, SamplingLawProbe,
    cross_section::CrossSectionIntegrand,
    evaluate_from_source_precise, parameterize,
    sampling::{maps::SamplingEvaluationError, selection::SamplingCatalogueEntry},
};
use crate::{
    feyngen::GenerationType,
    graph::{Graph, parse::IntoGraph},
    initialisation::test_initialise,
    integrands::evaluation::EvaluationMetaData,
    model::Model,
    momentum::{
        Rotation, RotationMethod, ThreeMomentum,
        sample::{ExternalIndex, LoopMomenta},
    },
    processes::{CrossSection, ProcessDefinition},
    settings::{GlobalSettings, RuntimeSettings, runtime::Precision},
    utils::{ArbPrec, F, FloatLike, QuadFloat, SamplingPrecision},
};
use color_eyre::Result;
use spenso::algebra::complex::Complex;
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::numerical_integration::Sample;

struct SimpleSamplingFixture {
    model: Model,
    integrand: CrossSectionIntegrand,
    settings: RuntimeSettings,
}

impl SimpleSamplingFixture {
    fn new() -> Result<Self> {
        test_initialise()?;
        let model = Model::from_str(
            include_str!("../../../../../assets/models/json/scalars/scalars_2p_3p.json")
                .to_string(),
            "json",
        )?;
        let graph: Graph = r#"digraph ordinary_sampling_bubble {
            num=1; edge [pdg=1001]; node [num=1];
            ext [style=invis, is_cut=0];
            ext -> a [id=0];
            a -> b [id=1, lmb_id=0];
            a -> b [id=2];
            b -> ext [id=3];
        }"#
        .into_graph(&model)?;
        let mut settings: RuntimeSettings = toml::from_str(
            r#"
[kinematics]
e_cm = 5.0
[kinematics.externals]
type = "constant"
[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0]]
helicities = ["summed_averaged"]
"#,
        )?;
        settings.stability.rotation_axis.clear();
        settings.stability.escalate_if_exact_zero = false;
        settings.stability.loop_momenta_norm_escalation_factor = 0.0;
        settings.subtraction.disable_threshold_subtraction = true;
        let mut global = GlobalSettings::default();
        global.generation.uv.softct = false;
        global.generation.uv.subtract_uv = false;
        let definition = ProcessDefinition::from_graph_list(
            std::slice::from_ref(&graph),
            GenerationType::CrossSection,
            &model,
        )?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .stack_size(256 * 1024 * 1024)
            .build()?;
        let mut cross_section =
            CrossSection::from_graph_list("ordinary_sampling".into(), vec![graph], &model)?;
        cross_section.preprocess(
            &model,
            &definition,
            &global.generation,
            (&settings).into(),
            &pool,
        )?;
        cross_section.build_integrand(&model, &definition, &global, (&settings).into(), &pool)?;
        let ProcessIntegrand::CrossSection(integrand) = cross_section.integrand.take().unwrap()
        else {
            unreachable!("scalar cross-section fixture")
        };
        Ok(Self {
            model,
            integrand,
            settings,
        })
    }

    fn configure(&mut self, sampling: &str, indices: Vec<usize>) -> Result<Sample<F<f64>>> {
        self.integrand.settings = self.settings.clone();
        self.integrand.settings.sampling = toml::from_str(sampling)?;
        self.integrand.warm_up(&self.model)?;
        Ok(Sample::Uniform(
            F(1.0),
            indices,
            [0.23, 0.31, 0.47].map(F).to_vec(),
        ))
    }

    fn poison_higher_sampling_bridges(&mut self) {
        let setup = self.integrand.get_graph_mut(0).sampling_setup_mut();
        let failure = SamplingEvaluationError::Unrepresentable {
            operation: "ordinary-source regression",
            detail: "physical rescue must reuse its draw without remapping".into(),
        };
        setup.sampling_bridge_quad.set(Err(failure.clone()));
        setup.sampling_bridge_fixed256.set(Err(failure.clone()));
        setup.sampling_bridge_arb.set(Err(failure));
    }

    fn assert_retained_lane<T: FloatLike>(
        &mut self,
        source: &EvaluationSource<'_>,
        canonical: &GammaLoopSample<ArbPrec>,
    ) -> Result<()> {
        let actual = source.build_gamma_sample::<T, _>(
            &mut self.integrand,
            &mut EvaluationMetaData::new_empty(),
        )?;
        let expected = canonical.materialize::<T>(
            GammaLoopSample::<T>::relative_accuracy_budget(&self.integrand.settings),
        )?;
        for method in [
            RotationMethod::Identity,
            RotationMethod::EulerAngles(0.1, 0.2, 0.3),
            RotationMethod::Pi2Z,
        ] {
            let rotation = Rotation::new(method);
            let actual = actual.rotate(&rotation, 17, 23);
            let expected = expected.rotate(&rotation, 17, 23);
            assert_eq!(actual.groups.len(), expected.groups.len());
            for ((actual_group, actual_rows), (expected_group, expected_rows)) in
                actual.groups.iter().zip(&expected.groups)
            {
                assert_eq!(actual_group, expected_group);
                assert_eq!(actual_rows.len(), expected_rows.len());
                for (actual, expected) in actual_rows.iter().zip(expected_rows) {
                    assert_eq!(actual.graph_id, expected.graph_id);
                    assert_eq!(actual.channel_id, expected.channel_id);
                    assert_eq!(actual.sample.loop_moms(), expected.sample.loop_moms());
                    assert_eq!(
                        actual.sample.external_moms(),
                        expected.sample.external_moms()
                    );
                    assert_eq!(actual.sample.jacobian(), expected.sample.jacobian());
                    assert_eq!(actual.integrand_prefactor, expected.integrand_prefactor);
                    assert!(actual.prepared_lu_hosts.is_empty());
                }
            }
        }
        Ok(())
    }
}

#[test]
fn ordinary_sources_preserve_completed_rows_through_native_rescue_and_rotations() {
    std::thread::Builder::new()
        .name("ordinary-fixed-source".into())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| -> Result<()> {
            let mut fixture = SimpleSamplingFixture::new()?;
            for (graphs, channels, multichannel, weight) in [
                ("summed", "summed", false, "map_density"),
                ("monte_carlo", "summed", false, "map_density"),
                ("monte_carlo", "monte_carlo", true, "map_density"),
                ("monte_carlo", "monte_carlo", true, "ose"),
                ("monte_carlo", "summed", true, "map_density"),
                ("monte_carlo", "summed", true, "ose"),
            ] {
                let sampling = format!(
                    "graphs = '{graphs}'\nsampling_multichanneling = {multichannel}\nsampling_channels = '{channels}'\nsampling_channel_weight = '{weight}'\ndefault_channel_selection = ['auto:lmb']"
                );
                let mut indices = if graphs == "monte_carlo" {
                    vec![0]
                } else {
                    vec![]
                };
                if multichannel && channels == "monte_carlo" {
                    indices.push(0);
                }
                let input = fixture.configure(&sampling, indices)?;
                assert_eq!(
                    fixture.integrand.sampling_source_policy()?.0,
                    SamplingPrecision::Double,
                    "{sampling}"
                );
                let source = EvaluationSource::XSpace(&input);
                let direct = parameterize::<f64, _>(
                    &input,
                    &mut fixture.integrand,
                    &mut EvaluationMetaData::new_empty(),
                )?;
                let simple = source.prepare_simple_draw(
                    &mut fixture.integrand,
                    &mut EvaluationMetaData::new_empty(),
                )?;
                assert_eq!(format!("{simple:?}"), format!("{direct:?}"));
                assert!(simple.row_count() > 0);
                let canonical = simple.clone().into_canonical(
                    &fixture.integrand.settings.kinematics.externals,
                    fixture.integrand.get_dependent_momenta_constructor(),
                )?;
                let reference = fixture.integrand.get_graph(0).sampling_setup().clone();
                fixture.poison_higher_sampling_bridges();
                let prepared = EvaluationSource::PreparedSimple {
                    sample: &simple,
                    canonical: None,
                    original: &source,
                };
                fixture.assert_retained_lane::<f64>(&prepared, &canonical)?;
                fixture.assert_retained_lane::<QuadFloat>(&prepared, &canonical)?;
                fixture.assert_retained_lane::<ArbPrec>(&prepared, &canonical)?;

                if multichannel {
                    let probe = SamplingLawProbe {
                        reference: &reference,
                        channel_id: SamplingChannelId(0),
                        radial_coordinate: 0,
                        retry_lower_half: true,
                        retry_quad: true,
                        sign_coordinate: None,
                        alternating_graph_sign: false,
                        calls: std::array::from_fn(|_| AtomicUsize::new(0)),
                    };
                    let value = evaluate_from_source_precise(
                        &mut fixture.integrand,
                        EvaluationTarget::SamplingLaw(&probe),
                        source,
                        F(1.0),
                        false,
                        Complex::new_zero(),
                    )?
                    .try_into_f64()?;
                    assert_eq!(value.evaluation_metadata.final_precision(), Some(Precision::Arb));
                    for calls in &probe.calls {
                        assert_eq!(calls.load(Ordering::Relaxed), simple.row_count());
                    }
                    let expected = if channels == "summed" {
                        1.0
                    } else {
                        let row = &simple.groups[0].1[0];
                        let coordinates = row.sample.loop_moms().iter().flat_map(|momentum| {
                            [momentum.px.0, momentum.py.0, momentum.pz.0]
                        }).collect::<Vec<_>>();
                        reference.sampling_bridge::<f64>()?.partition(&coordinates)?.weight(0).unwrap()
                    };
                    assert!((value.integrand_result.re.0 - expected).abs() < 1.0e-8);
                }
            }
            Ok(())
        })
        .unwrap()
        .join()
        .unwrap()
        .unwrap();
}

#[test]
fn ordinary_source_admission_retains_strict_budget_range_and_advanced_fallbacks() {
    std::thread::Builder::new()
        .name("ordinary-source-admission".into())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| -> Result<()> {
            let mut fixture = SimpleSamplingFixture::new()?;
            for multichannel in [false, true] {
                let fallback = if multichannel {
                    SamplingPrecision::Fixed256
                } else {
                    SamplingPrecision::Arb
                };
                let sampling = format!(
                    "graphs = 'monte_carlo'\nsampling_multichanneling = {multichannel}\nsampling_channels = 'monte_carlo'\ndefault_channel_selection = ['auto:lmb']"
                );
                let indices = if multichannel { vec![0, 0] } else { vec![0] };
                fixture.configure(&sampling, indices)?;
                assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Double);
                let levels = fixture.integrand.settings.stability.levels.clone();
                fixture.integrand.settings.stability.levels[2].required_precision_for_re = 1.0e-20;
                assert_eq!(fixture.integrand.sampling_source_policy()?.0, fallback);
                fixture.integrand.settings.stability.levels = levels;

                let externals = fixture.integrand.settings.kinematics.externals.clone();
                let mut precise = externals.get_dependent_externals::<ArbPrec>(
                    fixture.integrand.get_dependent_momenta_constructor(),
                )?;
                precise[ExternalIndex(0)].spatial.px = F::<ArbPrec>::default().from_usize(2).powi(2000);
                let crate::settings::runtime::kinematic::Externals::Constant { arb_cache, .. } =
                    &mut fixture.integrand.settings.kinematics.externals;
                arb_cache.set(precise);
                assert_eq!(fixture.integrand.sampling_source_policy()?.0, fallback);
                fixture.integrand.settings.kinematics.externals = externals;
                assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Double);
            }
            let _ = fixture.configure(
                "graphs = 'monte_carlo'\nsampling_multichanneling = true\nsampling_channels = 'monte_carlo'\ndefault_channel_selection = ['auto:lmb']",
                vec![0, 0],
            )?;
            let setup = fixture.integrand.get_graph_mut(0).sampling_setup_mut();
            setup.sampling_catalogue.as_mut().unwrap().entries.push(SamplingCatalogueEntry::Surface {
                edges: vec![1, 2],
                parent_lmb: vec![1],
            });
            assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Fixed256);
            fixture.integrand.get_graph_mut(0).sampling_setup_mut().sampling_catalogue.as_mut().unwrap().entries.pop();
            assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Double);

            fixture.integrand.settings.sampling = toml::from_str(
                "graphs = 'monte_carlo'\nsampling_multichanneling = true\nsampling_channels = 'monte_carlo'\nsampling_channel_weight = 'singularity_proxy'\ndefault_channel_selection = ['auto:lmb']",
            )?;
            assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Fixed256);

            // A cut map shares the same catalogue axis with ordinary LMBs,
            // but its root-dependent proposal must keep the advanced source.
            let term = fixture.integrand.get_graph(0);
            let parent = term.graph.loop_momentum_basis.loop_edges.iter().map(|edge| edge.0).collect::<Vec<_>>();
            let (cut_id, surface) = term.cut_esurface.iter_enumerated().next().unwrap();
            let edges = surface.energies.iter().map(|edge| edge.0.to_string()).collect::<Vec<_>>().join(",");
            let sampling = format!(
                "graphs = 'monte_carlo'\nsampling_multichanneling = true\nsampling_channels = 'monte_carlo'\ndefault_channel_selection = ['auto:lmb', 'physical_cut']\n[channel_definitions.{}.physical_cut]\naround = 'phase_space(cut({edges}))'\nparent_lmb = {parent:?}\nsubspace_lmb = {parent:?}\non_cut = [{}]",
                term.graph.name, cut_id.0,
            );
            let input = fixture.configure(&sampling, vec![0, 0])?;
            assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Quad);
            assert!(!EvaluationSource::XSpace(&input).uses_binary64_source(&fixture.integrand));
            Ok(())
        })
        .unwrap()
        .join()
        .unwrap()
        .unwrap();
}

#[test]
fn ordinary_sources_round_the_original_arb_external_authority_in_every_frame() {
    std::thread::Builder::new()
        .name("ordinary-external-authority".into())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| -> Result<()> {
            let mut fixture = SimpleSamplingFixture::new()?;
            let alternate = fixture.integrand.get_graph(0).sampling_setup().all_bases
                .iter_enumerated()
                .find_map(|(id, basis)| (basis.loop_edges != fixture.integrand.get_graph(0).graph.loop_momentum_basis.loop_edges).then_some(id))
                .expect("bubble fixture must have an alternate affine LMB");
            for multichannel in [false, true] {
                let sampling = format!(
                    "graphs = 'monte_carlo'\nsampling_multichanneling = {multichannel}\nsampling_channels = 'monte_carlo'\nlmb_basis_ids = {{ ordinary_sampling_bubble = [{}] }}",
                    usize::from(alternate),
                );
                let input = fixture.configure(&sampling, if multichannel { vec![0, 0] } else { vec![0] })?;
                let stale = fixture.integrand.settings.kinematics.externals.get_dependent_externals::<f64>(
                    fixture.integrand.get_dependent_momenta_constructor(),
                )?;
                let mut original = fixture.integrand.settings.kinematics.externals.get_dependent_externals::<ArbPrec>(
                    fixture.integrand.get_dependent_momenta_constructor(),
                )?;
                let one = original[ExternalIndex(0)].spatial.px.one();
                let precise = &one / one.from_usize(8)
                    + one.from_usize(2).powi(-80)
                    + one.from_usize(2).powi(-320);
                original[ExternalIndex(0)].spatial.px = precise.clone();
                let crate::settings::runtime::kinematic::Externals::Constant { arb_cache, .. } =
                    &mut fixture.integrand.settings.kinematics.externals;
                arb_cache.set(original.clone());
                // Do not re-run external improvement: deliberately retain a
                // stale native cache and a different original Arb authority.
                fixture.integrand.warm_up_sampling()?;
                assert_eq!(fixture.integrand.sampling_source_policy()?.0, SamplingPrecision::Double);
                let rounded: crate::momentum::sample::ExternalFourMomenta<F<f64>> = original.iter().map(|momentum| momentum.map_ref(&|value| value.into_ff64())).collect();
                assert_ne!(stale, rounded);
                let source = EvaluationSource::XSpace(&input);
                let simple = source.prepare_simple_draw(&mut fixture.integrand, &mut EvaluationMetaData::new_empty())?;
                // Source rows retain map metadata. Physical lanes obtain
                // external kinematics from the original Arb authority.
                let physical = EvaluationSource::materialize_simple::<f64, _>(
                    &simple, &fixture.integrand,
                )?;
                for row in &physical.groups[0].1 {
                    assert_eq!(row.sample.external_moms(), &rounded);
                }
                let parameterization = fixture.integrand.settings.sampling.get_parameterization_settings().unwrap();
                let (coordinates, _) = crate::utils::global_parameterize(
                    &[F(0.23), F(0.31), F(0.47)],
                    F(fixture.integrand.settings.kinematics.e_cm),
                    &parameterization,
                );
                let local: LoopMomenta<F<f64>> = coordinates.into_iter().map(ThreeMomentum::from).collect();
                let term = fixture.integrand.get_graph(0);
                let basis = &term.sampling_setup().all_bases[alternate];
                let expected: LoopMomenta<F<f64>> = term.graph.loop_momentum_basis.loop_edges.iter().map(|edge| {
                    let signature = &basis.edge_signatures[*edge];
                    signature.internal.apply_typed(&local) + signature.external.apply(&rounded.raw).spatial
                }).collect();
                let stale_frame: LoopMomenta<F<f64>> = term.graph.loop_momentum_basis.loop_edges.iter().map(|edge| {
                    let signature = &basis.edge_signatures[*edge];
                    signature.internal.apply_typed(&local) + signature.external.apply(&stale.raw).spatial
                }).collect();
                assert_ne!(expected, stale_frame, "the alternate LMB must expose the external affine shift");
                for (actual, expected) in simple.get_default_sample().loop_moms().iter().zip(expected.iter()) {
                    assert!((actual.px.0 - expected.px.0).abs() < 1.0e-12);
                    assert!((actual.py.0 - expected.py.0).abs() < 1.0e-12);
                    assert!((actual.pz.0 - expected.pz.0).abs() < 1.0e-12);
                }
                let canonical = simple.clone().into_canonical(
                    &fixture.integrand.settings.kinematics.externals,
                    fixture.integrand.get_dependent_momenta_constructor(),
                )?;
                assert_eq!(canonical.get_default_sample().external_moms(), &original);
                fixture.poison_higher_sampling_bridges();
                let prepared = EvaluationSource::PreparedSimple { sample: &simple, canonical: None, original: &source };
                fixture.assert_retained_lane::<f64>(&prepared, &canonical)?;
                fixture.assert_retained_lane::<QuadFloat>(&prepared, &canonical)?;
                fixture.assert_retained_lane::<ArbPrec>(&prepared, &canonical)?;
            }

            fixture.configure("graphs = 'monte_carlo'\nsampling_multichanneling = false", vec![0])?;
            let input = MomentumSpaceEvaluationInput {
                loop_momenta: vec![ThreeMomentum::new(F(0.1), F(-0.2), F(0.3))],
                integrator_weight: F(1.0),
                graph_id: Some(0),
                group_id: None,
                orientation: None,
                channel_id: None,
            };
            let candidate = fixture.integrand.get_graph(0).sampling_setup().all_bases[alternate].clone();
            let direct = EvaluationSource::Momentum { input: &input, candidate_lmb: None };
            let profile = EvaluationSource::Momentum { input: &input, candidate_lmb: Some(&candidate) };
            assert!(direct.uses_binary64_source(&fixture.integrand));
            assert!(!profile.uses_binary64_source(&fixture.integrand));
            let canonical = profile.prepare_draw(&mut fixture.integrand, &mut EvaluationMetaData::new_empty())?.unwrap();
            let expected = super::build_direct_gamma_sample::<ArbPrec, _>(
                &mut fixture.integrand,
                &input,
                Some(&candidate),
                &mut EvaluationMetaData::new_empty(),
            )?;
            assert_eq!(canonical.get_default_sample().loop_moms(), expected.get_default_sample().loop_moms());
            Ok(())
        })
        .unwrap()
        .join()
        .unwrap()
        .unwrap();
}
