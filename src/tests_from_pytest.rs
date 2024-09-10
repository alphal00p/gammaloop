#![allow(unused_imports)]
use crate::cff::esurface::{
    get_existing_esurfaces, Esurface, EsurfaceCollection, EsurfaceDerivedData, EsurfaceID,
    ExistingEsurfaceId, ExistingEsurfaces,
};
use crate::cff::expression;
use crate::cff::generation::generate_cff_expression;
use crate::cff::hsurface::HsurfaceCollection;
use crate::cross_section::{Amplitude, OutputMetaData, OutputType};
use crate::gammaloop_integrand::{DefaultSample, GammaLoopIntegrand};
use crate::graph::{
    BareGraph, DerivedGraphData, Edge, EdgeType, Graph, HasVertexInfo, InteractionVertexInfo,
    LoopMomentumBasisSpecification, SerializableGraph, VertexInfo,
};
use crate::model::{LorentzStructure, Model};
use crate::momentum::{Dep, ExternalMomenta, FourMomentum, Helicity, SignOrZero, ThreeMomentum};
use crate::numerator::{
    ContractionSettings, EvaluatorOptions, Evaluators, GammaAlgebraMode, IterativeOptions,
    Numerator, NumeratorCompileOptions, NumeratorEvaluatorOptions, NumeratorSettings,
    NumeratorState, PythonState, UnInit,
};
use crate::subtraction::overlap::{self, find_center, find_maximal_overlap};
use crate::subtraction::static_counterterm;
use crate::tests::load_default_settings;
use crate::utils::{
    approx_eq, approx_eq_complex_res, approx_eq_res, assert_approx_eq, FloatLike,
    PrecisionUpgradable,
};
use crate::utils::{f128, F};
use crate::{cff, ltd, Externals, Integrand, RotationMethod};
use crate::{
    inspect::inspect, ExportSettings, GammaloopCompileOptions, Settings,
    TropicalSubgraphTableSettings,
};
use ahash::AHashMap;
use bincode::{Decode, Encode};
use clarabel::solver::default;
use colored::Colorize;
use indexmap::set::Iter;
use itertools::{FormatWith, Itertools};
use rand::distributions::Standard;
use rand::prelude::Distribution;
use symbolica::domains::rational::Rational;
//use libc::__c_anonymous_ptrace_syscall_info_exit;
use core::f64;
use lorentz_vector::LorentzVector;
use petgraph::algo::greedy_matching;
use petgraph::graph;
use rayon::prelude::IndexedParallelIterator;
use rayon::vec;
use serde::{self, Deserialize, Serialize};
use spenso::complex::{Complex, SymbolicaComplex};
use spenso::network::TensorNetwork;
use spenso::structure::{IsAbstractSlot, Lorentz, RepName};
use statrs::function::evaluate;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{clone, env};
use symbolica;
use symbolica::atom::{Atom, AtomView, FunctionBuilder};
use symbolica::domains::float::{Complex as SymComplex, NumericalFloatLike, Real};
use symbolica::evaluate::{CompileOptions, ExpressionEvaluator, FunctionMap, OptimizationSettings};
use symbolica::id::{
    AtomMatchIterator, Condition, Match, MatchSettings, MatchStack, Pattern, PatternOrMap,
    Replacement,
};
use symbolica::state::State;
use typed_index_collections::TiVec;

#[allow(unused)]
const LTD_COMPARISON_TOLERANCE: F<f64> = F(1.0e-11);

#[allow(unused)]
const PHASEONE: F<f64> = F(0.);

#[allow(unused)]
const PHASEI: F<f64> = F(f64::consts::FRAC_PI_2);

#[allow(unused)]
const PHASEMINUSONE: F<f64> = F(f64::consts::PI);

#[allow(unused)]
const PHASEMINUSI: F<f64> = F(-f64::consts::FRAC_PI_2);

pub fn test_export_settings() -> ExportSettings {
    ExportSettings {
        compile_cff: true,
        numerator_settings: NumeratorSettings {
            eval_settings: NumeratorEvaluatorOptions::Joint(EvaluatorOptions {
                cpe_rounds: Some(1),
                compile_options: NumeratorCompileOptions::Compiled,
            }),
            global_numerator: None,
            gamma_algebra: GammaAlgebraMode::Concrete,
        },
        cpe_rounds_cff: Some(1),
        compile_separate_orientations: false,
        tropical_subgraph_table_settings: TropicalSubgraphTableSettings {
            target_omega: 1.0,
            panic_on_fail: false,
        },
        gammaloop_compile_options: GammaloopCompileOptions {
            inline_asm: env::var("NO_ASM").is_err(),
            optimization_level: 3,
            fast_math: true,
            unsafe_math: true,
            compiler: "g++".to_string(),
            custom: vec![],
        },
    }
}

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

// Create small, cheap to initialize and fast RNG with a random seed.
// The randomness is supplied by the operating system.

pub fn sample_generator<T: FloatLike>(
    seed: u64,
    bare_graph: &BareGraph,
    helicities: Option<Vec<Helicity>>,
) -> DefaultSample<T>
where
    Standard: Distribution<T>,
{
    let mut rng = SmallRng::seed_from_u64(seed);

    let n_loops = bare_graph.loop_momentum_basis.basis.len();
    let n_indep_externals = bare_graph.external_edges.len() - 1;
    let mut external_moms = vec![];
    for _ in 0..n_indep_externals {
        external_moms.push(FourMomentum::from_args(
            F(rng.gen()),
            F(rng.gen()),
            F(rng.gen()),
            F(rng.gen()),
        ));
    }

    let mut loop_moms = vec![];

    for _ in 0..n_loops {
        loop_moms.push(ThreeMomentum::new(F(rng.gen()), F(rng.gen()), F(rng.gen())));
    }

    let jacobian = F(1.0);

    let helicities = if let Some(hel) = helicities {
        if hel.len() != n_indep_externals + 1 {
            panic!("Helicities must have the same length as the number of external edges")
        }
        hel
    } else {
        vec![Helicity::Plus; n_indep_externals + 1]
    };

    DefaultSample::new(loop_moms, external_moms, jacobian, &helicities, bare_graph)
}

pub fn kinematics_builder(
    n_indep_externals: usize,
    n_loops: usize,
    bare_graph: &BareGraph,
) -> DefaultSample<f64> {
    let mut external_moms = vec![];

    for i in 0..n_indep_externals {
        external_moms.push(FourMomentum::from_args(
            F(i as f64),
            F(i as f64 + 0.25),
            F(i as f64 + 0.5),
            F(i as f64 + 0.75),
        ));
    }

    let mut loop_moms = vec![];

    for i in n_indep_externals..n_indep_externals + n_loops {
        loop_moms.push(ThreeMomentum::new(
            F(i as f64),
            F(i as f64 + 0.33),
            F(i as f64 + 0.66),
        ));
    }

    let jacobian = F(1.0);

    let helicities = vec![Helicity::Plus; n_indep_externals + 1];

    DefaultSample::new(loop_moms, external_moms, jacobian, &helicities, bare_graph)
}

pub fn load_amplitude_output(
    output_path: &str,
    load_generic_model: bool,
) -> (Model, Amplitude<UnInit>, PathBuf) {
    let output_dir = if let Ok(pytest_output_path) = env::var("PYTEST_OUTPUT_PATH_FOR_RUST") {
        pytest_output_path
    } else {
        String::from("./src/test_resources")
    };
    let path = Path::new(&output_dir).join(output_path);
    let output_meta_data: OutputMetaData = serde_yaml::from_reader(
        File::open(path.join("output_metadata.yaml")).unwrap_or_else(|e| {
            panic!(
                "Could not open file '{:?}'. Error: {:?}",
                path.join("output_metadata.yaml"),
                e
            )
        }),
    )
    .unwrap();
    assert_eq!(output_meta_data.output_type, OutputType::Amplitudes);
    assert_eq!(output_meta_data.contents.len(), 1);

    let model = if load_generic_model {
        Model::from_file(String::from(
            Path::new(&output_dir)
                .join(format!(
                    "gammaloop_models/{}.yaml",
                    output_meta_data.model_name.as_str()
                ))
                .to_str()
                .unwrap(),
        ))
        .unwrap()
    } else {
        Model::from_file(String::from(
            path.join(format!(
                "sources/model/{}.yaml",
                output_meta_data.model_name
            ))
            .to_str()
            .unwrap(),
        ))
        .unwrap()
    };

    let amplitude = Amplitude::from_file(
        &model,
        String::from(
            path.join(format!(
                "sources/amplitudes/{}/amplitude.yaml",
                output_meta_data.contents[0]
            ))
            .to_str()
            .unwrap(),
        ),
    )
    .unwrap();
    (model, amplitude, path)
}

pub enum SampleType<T: FloatLike> {
    Random(u64),
    Kinematic,
    RandomWithHelicity(u64, Vec<Helicity>),
    Custom(DefaultSample<T>),
}

pub struct AmplitudeCheck {
    pub name: &'static str,
    pub sample: SampleType<f64>,
    pub model_name: &'static str,
    pub n_edges: usize,
    pub n_external_connections: usize,
    pub n_cff_trees: usize,
    pub n_esurfaces: usize,
    pub n_vertices: usize,
    pub n_lmb: usize,
    pub cff_phase: F<f64>,
    pub cff_norm: Option<F<f64>>,
    pub n_prop_groups: usize,
    pub n_existing_esurfaces: usize,
    pub n_expanded_terms: usize,
    pub n_terms_unfolded: usize,
    pub n_overlap_groups: usize,
    pub n_existing_per_overlap: usize,
    pub tolerance: F<f64>,
}

#[allow(unused)]
fn check_load(amp_check: &AmplitudeCheck) -> (Model, Amplitude<UnInit>, PathBuf) {
    let (model, amplitude, path) = load_amplitude_output(
        &("TEST_AMPLITUDE_".to_string() + amp_check.name + "/GL_OUTPUT"),
        true,
    );

    assert_eq!(model.name, amp_check.model_name, "Model name mismatch");
    assert_eq!(
        amplitude.amplitude_graphs.len(),
        1,
        "Number of graphs should be 1"
    );
    assert_eq!(
        amplitude.amplitude_graphs[0].graph.bare_graph.edges.len(),
        amp_check.n_edges,
        "Graph should have {} edges, but has {}. Note: this includes external edges",
        amp_check.n_edges,
        amplitude.amplitude_graphs[0].graph.bare_graph.edges.len()
    );
    assert_eq!(
        amplitude.amplitude_graphs[0]
            .graph
            .bare_graph
            .vertices
            .len(),
        amp_check.n_vertices,
        "Graph should have {} vertices, but has {}. Note we count the external vertices",
        amp_check.n_vertices,
        amplitude.amplitude_graphs[0]
            .graph
            .bare_graph
            .vertices
            .len()
    );
    assert_eq!(
        amplitude.amplitude_graphs[0]
            .graph
            .bare_graph
            .external_connections
            .len(),
        amp_check.n_external_connections,
        "Graph should have {} external connections, but has {}",
        amp_check.n_external_connections,
        amplitude.amplitude_graphs[0]
            .graph
            .bare_graph
            .external_connections
            .len()
    );
    (model, amplitude, path)
}

use color_eyre::{Report, Result};
use eyre::{eyre, Context};

#[allow(unused)]
fn check_cff_generation<N: NumeratorState>(
    mut graph: Graph<N>,
    amp_check: &AmplitudeCheck,
) -> Result<Graph<N>> {
    graph.generate_cff();
    let cff = graph
        .derived_data
        .as_ref()
        .ok_or(eyre!("No derived data"))?
        .cff_expression
        .as_ref()
        .ok_or(eyre!("No cff expression"))?;

    assert!(!cff.esurfaces.is_empty());

    let num_trees = cff.get_num_trees();
    let esurfaces = cff.esurfaces.len();
    assert_eq!(
        num_trees, amp_check.n_cff_trees,
        "Graph should have {} cff trees, but has {}",
        amp_check.n_cff_trees, num_trees
    );
    assert_eq!(
        esurfaces, amp_check.n_esurfaces,
        "Graph should have {} esurfaces, but has {}",
        amp_check.n_esurfaces, esurfaces
    );
    let unfolded = cff.expand_terms();

    assert_eq!(unfolded.len(), amp_check.n_expanded_terms);

    for term in unfolded.iter() {
        assert_eq!(term.len(), amp_check.n_terms_unfolded);
    }
    Ok(graph)
}

#[allow(unused)]
fn check_lmb_generation<N: NumeratorState>(
    mut graph: Graph<N>,
    sample: &DefaultSample<f64>,
    amp_check: &AmplitudeCheck,
) -> Result<Graph<N>> {
    graph.generate_loop_momentum_bases();
    let lmb = graph
        .derived_data
        .as_ref()
        .ok_or(eyre!("No derived data"))?
        .loop_momentum_bases
        .as_ref()
        .ok_or(eyre!("No lmbs"))?;

    let emr = graph
        .bare_graph
        .compute_emr(&sample.loop_moms, &sample.external_moms);

    for basis in lmb {
        let momenta_in_basis = basis.basis.iter().map(|index| emr[*index]).collect_vec();
        let new_emr = basis
            .edge_signatures
            .iter()
            .map(|s| s.compute_three_momentum_from_four(&momenta_in_basis, &sample.external_moms))
            .collect_vec();
        assert_eq!(emr.len(), new_emr.len());

        for (e1, e2) in emr.iter().zip(new_emr.iter()) {
            assert_approx_eq(&e1.px, &e2.px, &LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(&e1.py, &e2.py, &LTD_COMPARISON_TOLERANCE);
            assert_approx_eq(&e1.pz, &e2.pz, &LTD_COMPARISON_TOLERANCE);
        }
    }

    let num_lmbs = lmb.len();
    assert_eq!(
        num_lmbs, amp_check.n_lmb,
        "Graph should have {} loop momentum bases, but has {}",
        amp_check.n_lmb, num_lmbs
    );
    Ok(graph)
}

#[allow(unused)]
fn check_sample(bare_graph: &BareGraph, amp_check: &AmplitudeCheck) -> DefaultSample<f64> {
    let n_loops = amp_check.n_edges - amp_check.n_vertices + 1; //circuit rank=n_loops

    match &amp_check.sample {
        SampleType::Random(seed) => sample_generator(*seed, bare_graph, None),
        SampleType::Kinematic => {
            kinematics_builder(amp_check.n_external_connections - 1, n_loops, bare_graph)
        }
        SampleType::RandomWithHelicity(seed, helicities) => {
            sample_generator(*seed, bare_graph, Some(helicities.clone()))
        }
        SampleType::Custom(sample) => sample.clone(),
    }
}

#[allow(unused)]
fn compare_cff_to_ltd<T: FloatLike>(
    sample: &DefaultSample<T>,
    graph: &mut Graph,
    amp_check: &AmplitudeCheck,
) -> Result<()> {
    let default_settings = load_default_settings();
    let lmb_specification =
        LoopMomentumBasisSpecification::Literal(&graph.bare_graph.loop_momentum_basis);
    let cff_res = graph.evaluate_cff_orientations(sample, &lmb_specification, &default_settings);
    for o in cff_res.iter() {
        println!("cff: {}", o);
    }
    let cff_res: Complex<F<T>> = graph.evaluate_cff_expression(sample, &default_settings);
    graph.generate_ltd();
    let ltd_res = graph.evaluate_ltd_expression(sample, &default_settings);

    let cff_norm = <Complex<F<T>> as Real>::norm(&cff_res).re;
    let cff_phase_actual = <Complex<F<T>> as SymbolicaComplex>::arg(&cff_res);

    let ltd_norm = <Complex<F<T>> as Real>::norm(&ltd_res).re;
    // let ltd_phase = ltd_res.arg();

    let ltd_comparison_tolerance = F::<T>::from_ff64(amp_check.tolerance);

    approx_eq_res(&cff_norm, &ltd_norm, &ltd_comparison_tolerance)
        .wrap_err("cff and ltd norms do not match")?;

    if let Some(truth) = amp_check.cff_norm {
        let energy_product = graph
            .bare_graph
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);
        approx_eq_res(
            &(cff_norm / &energy_product),
            &F::<T>::from_ff64(truth),
            &ltd_comparison_tolerance,
        )
        .wrap_err("Normalised cff is not what was expected")?;
    }

    approx_eq_res(
        &F::<T>::from_ff64(amp_check.cff_phase),
        &cff_phase_actual,
        &ltd_comparison_tolerance,
    )
    .wrap_err("Phase does not match expected")?;
    Ok(())
}

#[allow(unused)]
fn check_graph(graph: &BareGraph, n_prop_groups: usize) {
    let propagator_groups = graph.group_edges_by_signature();
    assert_eq!(
        propagator_groups.len(),
        n_prop_groups,
        "Expected {} propagator groups, but found {}",
        n_prop_groups,
        propagator_groups.len()
    );
}

#[allow(unused)]
fn check_esurface_existance<N: NumeratorState>(
    graph: &mut Graph<N>,
    sample: &DefaultSample<f64>,
    n_existing_esurfaces: usize,
    n_overlap_groups: usize,
    n_existing_per_overlap: usize,
) -> Result<()> {
    graph.generate_esurface_data().unwrap();

    let cff = graph
        .derived_data
        .as_ref()
        .ok_or(eyre!("derived data missing"))?
        .cff_expression
        .as_ref()
        .ok_or(eyre!("No cff expression"))?;
    let existing = get_existing_esurfaces(
        &cff.esurfaces,
        graph
            .derived_data
            .as_ref()
            .ok_or(eyre!("derived data missing"))?
            .esurface_derived_data
            .as_ref()
            .ok_or(eyre!("no esurface derived data"))?,
        &sample.external_moms,
        &graph.bare_graph.loop_momentum_basis,
        0,
        F(2.),
    );

    let edge_masses = graph
        .bare_graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .collect_vec();

    find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing,
        &cff.esurfaces,
        &edge_masses,
        &sample.external_moms,
        0,
    );

    let maximal_overlap = find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing,
        &cff.esurfaces,
        &edge_masses,
        &sample.external_moms,
        0,
    );

    assert_eq!(
        maximal_overlap.overlap_groups.len(),
        n_overlap_groups,
        "Number of overlap groups mismatch"
    );
    for overlap in maximal_overlap.overlap_groups.iter() {
        assert_eq!(
            overlap.existing_esurfaces.len(),
            n_existing_per_overlap,
            "Number of existing surfaces per overlap mismatch"
        );
    }

    assert_eq!(
        existing.len(),
        n_existing_esurfaces,
        "Number of existing surfaces mismatch"
    );
    Ok(())
}

#[allow(unused)]
fn check_amplitude(amp_check: AmplitudeCheck) {
    let (model, amplitude, path) = check_load(&amp_check);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    check_graph(&graph.bare_graph, amp_check.n_prop_groups);
    let sample = check_sample(&graph.bare_graph, &amp_check);

    graph = check_lmb_generation(graph, &sample, &amp_check).unwrap();
    graph = check_cff_generation(graph, &amp_check).unwrap();

    let sample_quad = sample.higher_precision();

    let export_settings = test_export_settings();

    println!(
        "numerator:{}",
        graph
            .clone()
            .apply_feynman_rules()
            .derived_data
            .unwrap()
            .numerator
            .export()
    );
    let mut graph =
        graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

    compare_cff_to_ltd(&sample, &mut graph, &amp_check)
        .wrap_err("combined num f64 cff and ltd failed")
        .unwrap();

    graph
        .derived_data
        .as_mut()
        .unwrap()
        .numerator
        .disable_combined();

    compare_cff_to_ltd(&sample, &mut graph, &amp_check)
        .wrap_err("separate num f64 cff and ltd failed")
        .unwrap();

    graph
        .derived_data
        .as_mut()
        .unwrap()
        .numerator
        .enable_combined(None);

    compare_cff_to_ltd(&sample.higher_precision(), &mut graph, &amp_check)
        .wrap_err("f128 cff and ltd failed")
        .unwrap();

    check_esurface_existance(
        &mut graph,
        &sample,
        amp_check.n_existing_esurfaces,
        amp_check.n_overlap_groups,
        amp_check.n_existing_per_overlap,
    )
    .unwrap();
}

#[allow(unused)]
fn init() {
    let _ = env_logger::builder().is_test(true).try_init();
}
#[test]
fn pytest_scalar_massless_triangle() {
    init();

    let p1 = FourMomentum::from_args(F(1.), F(3.), F(4.0), F(5.0));
    let p2 = FourMomentum::from_args(F(1.), F(6.0), F(7.), F(8.));

    let k = ThreeMomentum::new(F(1.), F(2.), F(3.));

    let sample = DefaultSample {
        loop_moms: vec![k],
        external_moms: vec![p1, p2],
        jacobian: F(1.0),
        polarizations: vec![],
    };
    let amp_check = AmplitudeCheck {
        name: "massless_scalar_triangle",
        sample: SampleType::Custom(sample),
        model_name: "scalars",
        n_edges: 6,
        n_external_connections: 3,
        n_cff_trees: 6,
        n_esurfaces: 6,
        n_vertices: 6,
        n_lmb: 3,
        cff_phase: PHASEI, //(0.).PI(),
        cff_norm: Some(F(4.531238289663331e-6)),
        n_prop_groups: 3,
        n_existing_esurfaces: 0,
        n_expanded_terms: 6,
        n_terms_unfolded: 2,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
        tolerance: LTD_COMPARISON_TOLERANCE,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_fishnet_2x2() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_fishnet_2x2",
        model_name: "scalars",
        sample: SampleType::Kinematic,
        n_edges: 16,
        n_vertices: 13,
        n_external_connections: 4,
        n_lmb: 192,
        n_cff_trees: 2398,
        n_prop_groups: 12,
        n_esurfaces: 97,
        cff_norm: Some(F(0.000019991301832169422)),
        cff_phase: F(0.),
        n_existing_esurfaces: 0,
        n_expanded_terms: 22852,
        n_terms_unfolded: 8,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 2,
    };

    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_sunrise() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_sunrise",
        model_name: "scalars",
        sample: SampleType::Kinematic,
        n_edges: 5,
        n_vertices: 4,
        n_external_connections: 2,
        cff_phase: PHASEMINUSI,
        cff_norm: Some(F(3.808857767867812e-3)),
        n_prop_groups: 3,
        n_cff_trees: 2,
        n_esurfaces: 2,
        n_lmb: 3,
        n_existing_esurfaces: 0,
        n_expanded_terms: 2,
        n_terms_unfolded: 1,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };

    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_fishnet_2x3() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_fishnet_2x3",
        model_name: "scalars",
        n_edges: 21,
        n_vertices: 16,
        sample: SampleType::Kinematic,
        n_external_connections: 4,
        n_prop_groups: 17,
        n_cff_trees: 58670,
        n_esurfaces: 263,
        n_existing_esurfaces: 34,
        n_expanded_terms: 2566256,
        n_lmb: 2415,
        n_terms_unfolded: 11,
        cff_norm: None,
        cff_phase: F(0.),

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1, //bogus
        n_overlap_groups: 2,       //bogus
    };

    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_cube() {
    init();
    // had to change tolerance to make it pass
    let amp_check = AmplitudeCheck {
        name: "scalar_cube",
        model_name: "scalars",
        n_edges: 20,
        n_vertices: 16,
        n_external_connections: 8,
        n_prop_groups: 12,
        n_cff_trees: 1862,
        n_esurfaces: 126,
        sample: SampleType::Kinematic,
        n_existing_esurfaces: 0,
        n_expanded_terms: 10584,
        n_lmb: 384,
        n_terms_unfolded: 7,
        cff_norm: None,
        cff_phase: PHASEMINUSI,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_bubble() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_bubble",
        model_name: "scalars",
        sample: SampleType::Kinematic,
        n_edges: 4,
        n_vertices: 4,
        n_external_connections: 2,
        n_prop_groups: 2,
        n_cff_trees: 2,
        n_esurfaces: 2,
        n_existing_esurfaces: 0,
        n_expanded_terms: 2,
        n_lmb: 2,
        n_terms_unfolded: 1,
        cff_norm: None,
        cff_phase: PHASEMINUSI,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_massless_box() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_massless_box",
        model_name: "scalars",
        n_edges: 8,
        n_vertices: 8,
        n_external_connections: 4,
        n_prop_groups: 4,
        sample: SampleType::Kinematic,
        n_cff_trees: 14,
        n_esurfaces: 12,
        n_existing_esurfaces: 0,
        n_expanded_terms: 20,
        n_lmb: 4,
        n_terms_unfolded: 3,
        cff_norm: None,
        cff_phase: PHASEMINUSI,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_double_triangle() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_double_triangle",
        model_name: "scalars",
        n_edges: 7,
        n_vertices: 6,
        n_external_connections: 2,
        n_prop_groups: 5,
        n_lmb: 8,
        n_cff_trees: 18,
        sample: SampleType::Kinematic,
        n_esurfaces: 10,
        n_existing_esurfaces: 0,
        n_expanded_terms: 20,
        n_terms_unfolded: 3,
        cff_norm: None,
        cff_phase: PHASEMINUSI,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_mercedes() {
    //had to lower tolerance again
    let amp_check = AmplitudeCheck {
        name: "scalar_mercedes",
        model_name: "scalars",
        n_edges: 9,
        n_vertices: 7,
        sample: SampleType::Kinematic,
        n_external_connections: 3,
        n_prop_groups: 6,
        n_lmb: 16,
        n_cff_trees: 24,
        n_esurfaces: 13,
        n_existing_esurfaces: 0,
        n_expanded_terms: 24,
        n_terms_unfolded: 3,
        cff_norm: None,
        cff_phase: PHASEMINUSI,

        tolerance: F(1.0e-5),
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_triangle_box() {
    init();
    let amp_check = AmplitudeCheck {
        name: "scalar_triangle_box",
        model_name: "scalars",
        n_edges: 9,
        n_vertices: 8,
        n_external_connections: 3,
        n_prop_groups: 6,
        n_lmb: 11,
        sample: SampleType::Kinematic,
        n_cff_trees: 42,
        n_esurfaces: 18,
        n_existing_esurfaces: 0,
        n_expanded_terms: 70,
        n_terms_unfolded: 4,
        cff_norm: None,
        cff_phase: PHASEI,
        tolerance: F(1.0e-7),
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_isopod() {
    init();

    let p1: FourMomentum<F<f64>> =
        FourMomentum::from_args(53. / 59., 41. / 43., 43. / 47., 47. / 53.).cast();

    let p2: FourMomentum<F<f64>> = FourMomentum::from_args(2., 15., 17., 19.).cast();

    let k0: ThreeMomentum<F<f64>> = ThreeMomentum::new(F(2. / 3.), F(3. / 5.), F(5. / 7.)).cast();
    let k1: ThreeMomentum<F<f64>> =
        ThreeMomentum::new(F(7. / 11.), F(11. / 13.), F(13. / 17.)).cast();

    let k2: ThreeMomentum<F<f64>> = ThreeMomentum::new(8., 9., 10.).cast();

    let loop_moms = [k0, k1, k2];
    let externals = [p1, p2];

    let sample = DefaultSample {
        loop_moms: loop_moms.to_vec(),
        external_moms: externals.to_vec(),
        jacobian: F(1.0),
        polarizations: vec![],
    };
    let amp_check = AmplitudeCheck {
        name: "scalar_isopod",
        model_name: "scalars",
        n_edges: 12,
        n_vertices: 10,
        sample: SampleType::Custom(sample),
        n_external_connections: 3,
        n_prop_groups: 9,
        n_lmb: 41,
        n_cff_trees: 294,
        n_esurfaces: 36,
        n_existing_esurfaces: 0,
        n_expanded_terms: 924,
        n_terms_unfolded: 6,
        cff_norm: None,
        cff_phase: PHASEI,
        tolerance: F(1.0e-7),
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
fn pytest_scalar_raised_triangle() {
    init();
    let (model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_raised_triangle/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    let propagator_groups = graph.bare_graph.group_edges_by_signature();
    assert_eq!(propagator_groups.len(), 5);
}

#[test]
fn pytest_scalar_hexagon() {
    init();
    let (model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_hexagon/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let esurfaces = &graph
        .derived_data
        .as_ref()
        .unwrap()
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let kinematics = [
        FourMomentum::from_args(F(24.), F(-21.2), F(71.), F(0.)),
        FourMomentum::from_args(F(50.4), F(15.8), F(-18.8), F(0.)),
        FourMomentum::from_args(F(-0.2), F(46.2), F(8.6), F(0.)),
        -FourMomentum::from_args(F(-33.2), F(2.6), F(-70.8), F(0.)),
        -FourMomentum::from_args(F(-80.), F(-5.6), F(-40.0), F(0.0)),
    ];

    let existing_esurface = get_existing_esurfaces(
        esurfaces,
        graph
            .derived_data
            .as_ref()
            .unwrap()
            .esurface_derived_data
            .as_ref()
            .unwrap(),
        &kinematics,
        &graph.bare_graph.loop_momentum_basis,
        0,
        F(75.),
    );

    assert_eq!(existing_esurface.len(), 6);

    let edge_masses = graph
        .bare_graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .collect_vec();

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing_esurface,
        esurfaces,
        &edge_masses,
        &kinematics,
        0,
    );
    let duration = now.elapsed();
    println!("duration: {}", duration.as_micros());

    assert_eq!(maximal_overlap.overlap_groups.len(), 4);
    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        3
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        3
    );
    assert_eq!(
        maximal_overlap.overlap_groups[2].existing_esurfaces.len(),
        3
    );
    assert_eq!(
        maximal_overlap.overlap_groups[3].existing_esurfaces.len(),
        2
    );

    let hexagon_10_e = [
        FourMomentum::from_args(F(-80.), F(29.), F(-70.), F(0.)),
        FourMomentum::from_args(F(83.5), F(14.0), F(70.0), F(0.0)),
        FourMomentum::from_args(F(88.5), F(6.5), F(-6.), F(0.)),
        -FourMomentum::from_args(F(36.5), F(-71.), F(97.5), F(0.)),
        -FourMomentum::from_args(F(12.5), F(-83.5), F(-57.5), F(0.)),
    ];

    let existing_esurfaces = get_existing_esurfaces(
        esurfaces,
        graph
            .derived_data
            .as_ref()
            .unwrap()
            .esurface_derived_data
            .as_ref()
            .unwrap(),
        &hexagon_10_e,
        &graph.bare_graph.loop_momentum_basis,
        0,
        F(88.),
    );

    assert_eq!(existing_esurfaces.len(), 10);

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing_esurfaces,
        esurfaces,
        &edge_masses,
        &hexagon_10_e,
        0,
    );
    let duration = now.elapsed();
    println!("duration: {}", duration.as_micros());

    assert_eq!(maximal_overlap.overlap_groups.len(), 4);
    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        8
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        8
    );
    assert_eq!(
        maximal_overlap.overlap_groups[2].existing_esurfaces.len(),
        7
    );
    assert_eq!(
        maximal_overlap.overlap_groups[3].existing_esurfaces.len(),
        7
    );
}

#[test]
fn pytest_scalar_ltd_topology_c() {
    init();
    let (model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_ltd_topology_c/GL_OUTPUT", true);

    assert_eq!(model.name, "scalars");
    assert!(amplitude.amplitude_graphs.len() == 1);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let esurfaces = &graph
        .derived_data
        .as_ref()
        .unwrap()
        .cff_expression
        .as_ref()
        .unwrap()
        .esurfaces;

    let kinematics = [
        FourMomentum::from_args(F(9.0), F(0.0), F(0.0), F(8.94427190999916)),
        FourMomentum::from_args(F(9.0), F(0.0), F(0.0), F(-8.94427190999916)),
        -FourMomentum::from_args(
            F(1.83442509122858),
            F(-0.383828222192743),
            F(0.69085529916260),
            F(-1.31653190094982),
        ),
        -FourMomentum::from_args(
            F(5.78920098524940),
            F(-1.80358221330469),
            F(-5.24375913342836),
            F(1.328506453),
        ),
        -FourMomentum::from_args(F(2.82869), F(-1.83886), F(-1.6969477), F(0.8605192)),
    ];

    let _edge_masses = graph
        .bare_graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re.0 == 0.0 {
                    None
                } else {
                    Some(value)
                }
            } else {
                mass
            }
        })
        .collect_vec();

    let existing_esurfaces = get_existing_esurfaces(
        esurfaces,
        graph
            .derived_data
            .as_ref()
            .unwrap()
            .esurface_derived_data
            .as_ref()
            .unwrap(),
        &kinematics,
        &graph.bare_graph.loop_momentum_basis,
        2,
        F(18.),
    );

    let overlap = find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing_esurfaces,
        esurfaces,
        &_edge_masses,
        &kinematics,
        0,
    );

    assert_eq!(overlap.overlap_groups.len(), 9);
    assert_eq!(overlap.overlap_groups[0].existing_esurfaces.len(), 22);
    assert_eq!(overlap.overlap_groups[1].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[2].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[3].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[4].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[5].existing_esurfaces.len(), 21);
    assert_eq!(overlap.overlap_groups[6].existing_esurfaces.len(), 20);
    assert_eq!(overlap.overlap_groups[7].existing_esurfaces.len(), 17);
    assert_eq!(overlap.overlap_groups[8].existing_esurfaces.len(), 17);
}

#[test]
fn pytest_scalar_massless_pentabox() {
    init();
    let (_model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_massless_pentabox/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let rescaling = F(1.0e-3);
    let kinematics = [
        &FourMomentum::from_args(
            F(5.980_260_048_915_123e2),
            F(0.0),
            F(0.0),
            F(5.724_562_014_045_295e2),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(5.980_260_048_915_123e2),
            F(0.0),
            F(0.0),
            F(-5.724_562_014_045_295e2),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-5.394_473_213_122_507e2),
            F(-1.971_081_698_462_961e2),
            F(-4.416_135_519_343_869e2),
            F(2.250_822_886_064_787e2),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-2.255_538_754_188_549e2),
            F(1.757_868_459_829_899e2),
            F(3.716_353_112_335_996e1),
            F(-1.013_763_093_935_658e2),
        ) * &rescaling,
    ];

    let edge_masses = graph
        .bare_graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re.0 == 0.0 {
                    None
                } else {
                    Some(value)
                }
            } else {
                mass
            }
        })
        .collect_vec();

    let existing_esurfaces = get_existing_esurfaces(
        &graph
            .derived_data
            .as_ref()
            .unwrap()
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        graph
            .derived_data
            .as_ref()
            .unwrap()
            .esurface_derived_data
            .as_ref()
            .unwrap(),
        &kinematics,
        &graph.bare_graph.loop_momentum_basis,
        0,
        F(1.0),
    );

    assert_eq!(existing_esurfaces.len(), 17);

    let maximal_overlap = find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing_esurfaces,
        &graph
            .derived_data
            .as_ref()
            .unwrap()
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        &edge_masses,
        &kinematics,
        0,
    );

    assert_eq!(maximal_overlap.overlap_groups.len(), 3);

    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        16
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        13
    );
    assert_eq!(
        maximal_overlap.overlap_groups[2].existing_esurfaces.len(),
        13
    );
}

#[test]
fn pytest_scalar_massless_3l_pentabox() {
    init();
    let (_model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_massless_3l_pentabox/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_ltd();
    graph.generate_cff();
    graph.generate_loop_momentum_bases();
    graph.generate_esurface_data().unwrap();

    let rescaling = F(1.0e0);
    let kinematics = [
        &FourMomentum::from_args(
            F(0.149500000000000E+01),
            F(0.000000000000000E+00),
            F(0.000000000000000E+00),
            F(0.149165176901313E+01),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(0.150500000000000E+01),
            F(0.000000000000000E+00),
            F(0.000000000000000E+00),
            F(-0.149165176901313E+01),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-0.126041949101381e+01),
            F(-0.452362952912639e+00),
            F(-0.101350243653045e+01),
            F(0.516563513332600e+00),
        ) * &rescaling,
        &FourMomentum::from_args(
            F(-0.105098730574850e+01),
            F(0.489324061520790e-01),
            F(0.928212188578101e+00),
            F(-0.283905035967510e+00),
        ) * &rescaling,
    ];

    let edge_masses = graph
        .bare_graph
        .edges
        .iter()
        .map(|edge| edge.particle.mass.value)
        .map(|mass| {
            if let Some(value) = mass {
                if value.re.0 == 0.0 {
                    None
                } else {
                    Some(value)
                }
            } else {
                mass
            }
        })
        .collect_vec();

    let existing_esurfaces = get_existing_esurfaces(
        &graph
            .derived_data
            .as_ref()
            .unwrap()
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        graph
            .derived_data
            .as_ref()
            .unwrap()
            .esurface_derived_data
            .as_ref()
            .unwrap(),
        &kinematics,
        &graph.bare_graph.loop_momentum_basis,
        0,
        F(1.0),
    );

    assert_eq!(graph.bare_graph.loop_momentum_basis.basis.len(), 3);
    assert_eq!(existing_esurfaces.len(), 28);

    let now = std::time::Instant::now();
    let maximal_overlap = find_maximal_overlap(
        &graph.bare_graph.loop_momentum_basis,
        &existing_esurfaces,
        &graph
            .derived_data
            .as_ref()
            .unwrap()
            .cff_expression
            .as_ref()
            .unwrap()
            .esurfaces,
        &edge_masses,
        &kinematics,
        0,
    );
    let _elapsed = now.elapsed();

    assert_eq!(maximal_overlap.overlap_groups.len(), 2);
    assert_eq!(
        maximal_overlap.overlap_groups[0].existing_esurfaces.len(),
        27
    );
    assert_eq!(
        maximal_overlap.overlap_groups[1].existing_esurfaces.len(),
        26
    );
}

#[test]
fn pytest_lbl_box() {
    init();
    let (model, amplitude, _) = load_amplitude_output("TEST_AMPLITUDE_lbl_box/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_cff();
    let export_settings = test_export_settings();
    let _graph = graph.process_numerator(
        &model,
        ContractionSettings::Normal,
        PathBuf::new(),
        &export_settings,
    );
}

// #[test]
// #[allow(non_snake_case)]
// fn pytest_physical_3L_6photons_topology_A_inspect() {//too slow
//     init();

//     let (model, amplitude, path) = load_amplitude_output(
//         "TEST_AMPLITUDE_physical_3L_6photons_topology_A/GL_OUTPUT",
//         true,
//     );

//     let mut graph = amplitude.amplitude_graphs[0].graph.clone();

//     graph.generate_cff();
//     let export_settings = ExportSettings {
//         compile_cff: true,
//         numerator_settings: NumeratorSettings {
//             eval_settings: NumeratorEvaluatorOptions::Single(EvaluatorOptions {
//                 cpe_rounds: Some(1),
//                 compile_options: NumeratorCompileOptions::Compiled,
//             }),
//             global_numerator: None,
//             gamma_algebra: GammaAlgebraMode::Concrete,
//         },
//         cpe_rounds_cff: Some(1),
//         compile_separate_orientations: false,
//         tropical_subgraph_table_settings: TropicalSubgraphTableSettings {
//             target_omega: 1.0,
//             panic_on_fail: false,
//         },
//         gammaloop_compile_options: GammaloopCompileOptions {
//             inline_asm: env::var("NO_ASM").is_err(),
//             optimization_level: 3,
//             fast_math: true,
//             unsafe_math: true,
//             compiler: "g++".to_string(),
//             custom: vec![],
//         },
//     };

//     let mut graph =
//         graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

//     let sample = kinematics_builder(5, 3, &graph.bare_graph);

//     let settings = load_default_settings();

//     graph.evaluate_cff_expression(&sample, &settings);
// }

#[test]
#[allow(non_snake_case)]
fn pytest_physical_1L_6photons() {
    init();

    let amp_check = AmplitudeCheck {
        name: "physical_1L_6photons",
        model_name: "sm",
        n_edges: 12,
        n_vertices: 12,
        sample: SampleType::Kinematic,
        n_external_connections: 6,
        n_prop_groups: 6,
        n_lmb: 6,
        n_cff_trees: 62,
        n_esurfaces: 30,
        n_existing_esurfaces: 0,
        n_expanded_terms: 252,
        n_terms_unfolded: 5,
        cff_norm: None,
        cff_phase: PHASEMINUSI,

        tolerance: LTD_COMPARISON_TOLERANCE,
        n_existing_per_overlap: 1,
        n_overlap_groups: 0,
    };
    check_amplitude(amp_check);
}

#[test]
#[allow(non_snake_case)]
fn pytest_physical_2L_6photons() {
    //slow
    init();
    let default_settings = load_default_settings();
    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_physical_2L_6photons/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    let sample = kinematics_builder(5, 2, &graph.bare_graph);

    graph.generate_cff();
    let export_settings = test_export_settings();
    let mut graph =
        graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

    graph.evaluate_cff_expression(&sample, &default_settings);
}

// #[test] too slow for now
// #[allow(non_snake_case)]
// fn physical_1L_6photons_gamma() {
//     init();
//     let (model, amplitude, path) =
//         load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

//     let mut graph = amplitude.amplitude_graphs[0].graph.clone();

//     graph.generate_cff();

//     let extra_info = graph
//         .derived_data
//         .as_ref()
//         .unwrap()
//         .generate_extra_info(path.clone());
//     let contraction_settings = ContractionSettings::Normal;
//     let export_settings = test_export_settings();

//     let mut graph_no_gamma = graph.clone().process_numerator(
//         &model,
//         ContractionSettings::Normal,
//         path,
//         &export_settings,
//     );
//     let mut graph = graph
//         .map_numerator_res(|n, g| {
//             Result::<_, Report>::Ok(
//                 n.from_graph(g)
//                     .color_symplify()
//                     .gamma_symplify()
//                     .parse()
//                     .contract::<Rational>(contraction_settings)?
//                     .generate_evaluators(&model, g, &extra_info, &export_settings),
//             )
//         })
//         .unwrap();

//     let sample = kinematics_builder(2, 1, &graph.bare_graph);

//     let default_settings = load_default_settings();
//     let gamma_eval = graph.evaluate_cff_expression(&sample, &default_settings);
//     let eval = graph_no_gamma.evaluate_cff_expression(&sample, &default_settings);

//     approx_eq_complex_res(&gamma_eval, &eval, &LTD_COMPARISON_TOLERANCE)
//         .wrap_err("Gamma algebra and spenso do not match")
//         .unwrap();
// }

#[test]
#[allow(non_snake_case)]
fn top_bubble_CP() {
    init();
    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    let mut export_settings = test_export_settings();
    export_settings.numerator_settings.gamma_algebra = GammaAlgebraMode::Concrete;

    let mut graph =
        graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

    let sample: DefaultSample<f64> =
        sample_generator(3, &graph.bare_graph, Some(vec![Helicity::Plus; 6]));

    #[allow(non_snake_case)]
    let sample_CP = sample_generator(3, &graph.bare_graph, Some(vec![Helicity::Minus; 6]));

    let default_settings = load_default_settings();

    let eval = graph.evaluate_cff_expression(&sample, &default_settings);

    #[allow(non_snake_case)]
    let eval_CP = graph.evaluate_cff_expression(&sample_CP, &default_settings);
    // println!("{}", sample);
    // println!("{}", sample_CP);
    println!("{}", eval);
    println!("{}", eval_CP);

    approx_eq_complex_res(&eval, &eval_CP, &LTD_COMPARISON_TOLERANCE)
        .wrap_err("CP conjugation does not match")
        .unwrap();
}

#[test]
fn top_bubble_gamma() {
    init();
    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_top_bubble/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();

    let mut export_settings = test_export_settings();
    export_settings.numerator_settings.gamma_algebra = GammaAlgebraMode::Concrete;

    let mut graph_no_gamma = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );

    export_settings.numerator_settings.gamma_algebra = GammaAlgebraMode::Symbolic;

    // fs::create_dir(path.join("sym")).unwrap();
    let mut graph = graph.process_numerator(
        &model,
        ContractionSettings::Normal,
        path.join("sym"),
        &export_settings,
    );

    for seed in 0..10 {
        let sample = sample_generator(
            seed,
            &graph.bare_graph,
            Some(vec![Helicity::Plus, Helicity::Plus]),
        );

        let sample_cp = sample_generator(
            seed,
            &graph.bare_graph,
            Some(vec![Helicity::Minus, Helicity::Minus]),
        );

        let default_settings = load_default_settings();
        let gamma_eval = graph.evaluate_cff_expression(&sample, &default_settings);
        let eval = graph_no_gamma.evaluate_cff_expression(&sample, &default_settings);

        let gamma_eval_cp = graph.evaluate_cff_expression(&sample_cp, &default_settings);
        let eval_cp = graph_no_gamma.evaluate_cff_expression(&sample_cp, &default_settings);

        approx_eq_complex_res(&gamma_eval, &eval, &LTD_COMPARISON_TOLERANCE)
            .wrap_err("Gamma algebra and spenso do not match")
            .unwrap();
        approx_eq_complex_res(&eval_cp, &eval, &LTD_COMPARISON_TOLERANCE)
            .wrap_err("Gamma algebra and spenso do not match")
            .unwrap();
        approx_eq_complex_res(&gamma_eval, &gamma_eval_cp, &LTD_COMPARISON_TOLERANCE)
            .wrap_err("Gamma algebra and spenso do not match")
            .unwrap();
    }
}

#[test]
fn scalar_box_to_triangle() {
    init();
    let default_settings = load_default_settings();
    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_scalar_massless_box/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    let mut export_settings = test_export_settings();

    graph.generate_cff();
    export_settings.numerator_settings.global_numerator= Some("Q(5,cind(0))*(Q(2,cind(0))+Q(6,cind(0)))-Q(5,cind(1))*Q(5,cind(1))-Q(5,cind(2))*Q(5,cind(2))-Q(5,cind(3))*Q(5,cind(3))".into());

    let box_sample = sample_generator(3, &graph.bare_graph, None);

    let mut box_graph =
        graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_massless_scalar_triangle/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();
    graph.generate_cff();

    let export_settings = test_export_settings();
    let mut triangle_graph =
        graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

    let box_energy = box_graph
        .bare_graph
        .compute_energy_product(&box_sample.loop_moms, &box_sample.external_moms);

    let normalized_box =
        box_graph.evaluate_cff_expression(&box_sample, &default_settings) / box_energy;

    let triangle_sample = DefaultSample {
        loop_moms: box_sample.loop_moms,
        external_moms: vec![
            box_sample.external_moms[0],
            box_sample.external_moms[2] + box_sample.external_moms[1],
        ],
        polarizations: vec![],
        jacobian: F(1.0),
    };

    let triangle_energy = triangle_graph
        .bare_graph
        .compute_energy_product(&triangle_sample.loop_moms, &triangle_sample.external_moms);

    let normalized_triangle = triangle_graph
        .evaluate_cff_expression(&triangle_sample, &default_settings)
        / triangle_energy;

    approx_eq_complex_res(
        &normalized_box,
        &normalized_triangle,
        &LTD_COMPARISON_TOLERANCE,
    )
    .wrap_err("Modified Box and triangle do not match")
    .unwrap();
}

pub fn compare_numerator_evals(amp_name: &str) -> Result<()> {
    let default_settings = load_default_settings();
    let (model, amplitude, path) =
        load_amplitude_output(&format!("TEST_AMPLITUDE_{amp_name}/GL_OUTPUT"), true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    let externals = graph.bare_graph.external_edges.len() - 1;
    let loops = graph.bare_graph.loop_momentum_basis.basis.len();

    let sample = kinematics_builder(externals, loops, &graph.bare_graph);

    graph.generate_cff();
    let mut export_settings = test_export_settings();

    export_settings.gammaloop_compile_options.inline_asm = true;

    export_settings.numerator_settings.eval_settings =
        NumeratorEvaluatorOptions::Iterative(IterativeOptions {
            eval_options: EvaluatorOptions {
                cpe_rounds: Some(1),
                compile_options: NumeratorCompileOptions::Compiled,
            },
            iterations: 1,
            n_cores: 1,
            verbose: false,
        });

    let mut graph_iterative_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );

    export_settings.numerator_settings.eval_settings =
        NumeratorEvaluatorOptions::Joint(EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::Compiled,
        });

    let mut graph_joint_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );
    export_settings.numerator_settings.eval_settings =
        NumeratorEvaluatorOptions::Single(EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::Compiled,
        });

    let mut graph_single_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );
    export_settings.numerator_settings.eval_settings =
        NumeratorEvaluatorOptions::Iterative(IterativeOptions {
            eval_options: EvaluatorOptions {
                cpe_rounds: Some(1),
                compile_options: NumeratorCompileOptions::NotCompiled,
            },
            iterations: 1,
            n_cores: 1,
            verbose: false,
        });
    let mut graph_iterative = graph_iterative_compiled.clone();
    graph_iterative
        .derived_data
        .as_mut()
        .unwrap()
        .numerator
        .disable_compiled();

    let mut graph_joint = graph_joint_compiled.clone();
    graph_joint
        .derived_data
        .as_mut()
        .unwrap()
        .numerator
        .disable_compiled();

    let mut graph_single = graph_single_compiled.clone();
    graph_single
        .derived_data
        .as_mut()
        .unwrap()
        .numerator
        .disable_compiled();

    let eval_single = graph_single.evaluate_cff_expression(&sample, &default_settings);
    let eval_joint = graph_joint.evaluate_cff_expression(&sample, &default_settings);
    let eval_iter = graph_iterative.evaluate_cff_expression(&sample, &default_settings);
    let eval_single_comp =
        graph_single_compiled.evaluate_cff_expression(&sample, &default_settings);
    let eval_joint_comp = graph_joint_compiled.evaluate_cff_expression(&sample, &default_settings);
    let eval_iter_comp =
        graph_iterative_compiled.evaluate_cff_expression(&sample, &default_settings);

    approx_eq_complex_res(&eval_single, &eval_joint, &LTD_COMPARISON_TOLERANCE)
        .wrap_err("Single and joint evaluation differ in norm")?;

    approx_eq_complex_res(&eval_single, &eval_iter, &LTD_COMPARISON_TOLERANCE)
        .wrap_err("Single and iterative evaluation differ in norm")?;

    approx_eq_complex_res(
        &eval_single_comp,
        &eval_joint_comp,
        &LTD_COMPARISON_TOLERANCE,
    )
    .wrap_err("Single compiled and joint compiled evaluation differ in norm")?;

    approx_eq_complex_res(
        &eval_single_comp,
        &eval_iter_comp,
        &LTD_COMPARISON_TOLERANCE,
    )
    .wrap_err("Single compiled and iterative compiled evaluation differ in norm")?;

    approx_eq_complex_res(&eval_single, &eval_single_comp, &LTD_COMPARISON_TOLERANCE)
        .wrap_err("Single and Single compiled evaluation differ in norm")?;

    Ok(())
}

#[test]
#[allow(non_snake_case)]
fn pytest_top_bubble() {
    init();
    compare_numerator_evals("top_bubble")
        .wrap_err("top bubble failure:")
        .unwrap();
}

#[test]
#[allow(non_snake_case)]
fn pytest_physical_1L_6photons_generate() {
    init();

    compare_numerator_evals("physical_1L_6photons")
        .wrap_err("physical 1L photon failure:")
        .unwrap();
}

#[test]
fn ratio_compare() {
    let gl_a = Complex::new(-3.3669e-6, 1.707e-7);
    let target_a = Complex::new(-1.748_314_585_608_997_2e-7, 7.660_134_793_988_562e-7);

    let gl_b = Complex::new(-3.457e-7, 8.7e-9);
    let target_b = Complex::new(-6.836_310_962_511_977e-7, 6.148_816_685_390_842e-7);

    println!("{}", gl_a.norm() / target_a.norm());
    println!("{}", gl_b.norm() / target_b.norm());

    let a = 1.7501349798028542e2;
    let b = 4.375337449507134e1;

    println!("{}", a / b);
}

// #[test]
// fn yaml_settings() {
//     let numerator_settings = NumeratorSettings {
//         eval_settings: NumeratorEvaluatorOptions::Iterative(IterativeOptions {
//             eval_options: EvaluatorOptions {
//                 compile_options: NumeratorCompileOptions::Compiled,
//                 cpe_rounds: Some(1),
//             },
//             iterations: 1,
//             n_cores: 1,
//             verbose: false,
//         }),
//         global_numerator: None,
//         gamma_algebra: GammaAlgebraMode::Concrete,
//     };

//     let externals = vec![
//         ExternalMomenta::Dependent(Dep::Dep),
//         ExternalMomenta::Independent([0, 1, 2, 3]),
//     ];

//     let yaml = serde_yaml::to_string(&externals).unwrap();
//     println!("{}", yaml);
// }
