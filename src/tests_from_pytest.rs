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
    SerializableGraph, VertexInfo,
};
use crate::model::{LorentzStructure, Model};
use crate::momentum::{FourMomentum, Helicity, ThreeMomentum};
use crate::numerator::{
    ContractionSettings, EvaluatorOptions, Evaluators, Num, NumeratorCompileOptions,
    NumeratorEvaluatorOptions, NumeratorState, PythonState, UnInit,
};
use crate::subtraction::overlap::{self, find_center, find_maximal_overlap};
use crate::subtraction::static_counterterm;
use crate::tests::load_default_settings;
use crate::utils::{approx_eq, approx_eq_res, assert_approx_eq, FloatLike, PrecisionUpgradable};
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
use itertools::{FormatWith, Itertools};
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
use std::fs::File;
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
const LTD_COMPARISON_TOLERANCE: F<f64> = F(1.0e-8);

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
        numerator_settings: NumeratorEvaluatorOptions::Iterative {
            eval_options: EvaluatorOptions {
                cpe_rounds: Some(1),
                compile_options: NumeratorCompileOptions::Compiled,
            },
            iterations: 1,
            n_cores: 1,
            verbose: true,
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

pub struct AmplitudeCheck {
    pub name: &'static str,
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

use color_eyre::Result;
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

    kinematics_builder(amp_check.n_external_connections - 1, n_loops, bare_graph)
}

#[allow(unused)]
fn compare_cff_to_ltd<T: FloatLike>(
    sample: &DefaultSample<T>,
    graph: &mut Graph,
    amp_check: &AmplitudeCheck,
) -> Result<()> {
    let default_settings = load_default_settings();

    let cff_res: Complex<F<T>> = graph.evaluate_cff_expression(sample, &default_settings);
    graph.generate_ltd();
    let ltd_res = graph.evaluate_ltd_expression(&sample.loop_moms, &sample.external_moms);

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

    compare_cff_to_ltd(&sample_quad, &mut graph, &amp_check)
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
    let amp_check = AmplitudeCheck {
        name: "massless_scalar_triangle",
        model_name: "scalars",
        n_edges: 6,
        n_external_connections: 3,
        n_cff_trees: 6,
        n_esurfaces: 6,
        n_vertices: 6,
        n_lmb: 3,
        cff_phase: PHASEI, //(0.).PI(),
        cff_norm: Some(F(4.052725000745997e-5)),
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
    let amp_check = AmplitudeCheck {
        name: "scalar_isopod",
        model_name: "scalars",
        n_edges: 12,
        n_vertices: 10,
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

#[test]
#[allow(non_snake_case)]
fn pytest_physical_3L_6photons_topology_A_inspect() {
    init();

    let (model, amplitude, path) = load_amplitude_output(
        "TEST_AMPLITUDE_physical_3L_6photons_topology_A/GL_OUTPUT",
        true,
    );

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_cff();
    let export_settings = ExportSettings {
        compile_cff: true,
        numerator_settings: NumeratorEvaluatorOptions::Single(EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::Compiled,
        }),
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
    };

    let mut graph =
        graph.process_numerator(&model, ContractionSettings::Normal, path, &export_settings);

    let sample = kinematics_builder(5, 3, &graph.bare_graph);

    let settings = load_default_settings();

    graph.evaluate_cff_expression(&sample, &settings);
}

#[test]
#[allow(non_snake_case)]
fn pytest_physical_1L_6photons() {
    init();
    let _ = load_default_settings();
    let (model, amplitude, true_path) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

    let export_settings = test_export_settings();
    let amp = amplitude
        .export(
            true_path.to_str().as_ref().unwrap(),
            &model,
            &export_settings,
        )
        .unwrap();
    // .map(|a| a.map(|ag| ag.forget_type()));

    let mut settings: Settings = Default::default();

    match &mut settings.kinematics.externals {
        Externals::Constant {
            momenta,
            helicities,
        } => {
            *momenta = vec![
                [F(1.0), F(0.0), F(0.0), F(1.0)],
                [F(1.0), F(0.0), F(0.0), F(-1.0)],
                [F(-0.5), F(0.5), F(0.5), F(0.5)],
                [F(-0.5), F(-0.5), F(-0.5), F(-0.5)],
                [F(0.5), F(-0.5), F(-0.5), F(0.5)],
                [F(0.5), F(0.5), F(0.5), F(-0.5)],
            ];
            *helicities = vec![
                Helicity::Plus,
                Helicity::Plus,
                Helicity::Minus,
                Helicity::Minus,
                Helicity::Minus,
                Helicity::Minus,
            ];
        }
    }

    settings.stability.rotation_axis = vec![
        RotationMethod::Pi2Z,
        RotationMethod::Pi2X,
        RotationMethod::Pi2Y,
    ];

    let mut integrand = Integrand::GammaLoopIntegrand(
        GammaLoopIntegrand::amplitude_integrand_constructor(amp, settings.clone()),
    );
    let point = vec![F(0.123), F(0.3242), F(0.4233)];

    let _a = inspect(&settings, &mut integrand, point, &[], false, true, true);
    // integrand.inspect();

    // amp.generate_integrand(path_to_settings)

    // let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    // graph.generate_cff();

    // println!("{}", serde_yaml::to_string(&export_settings).unwrap());
    // let mut graph = graph.process_numerator(
    //     &model,
    //     ContractionSettings::Normal,
    //     true_path,
    //     &export_settings,
    // );

    // let sample = kinematics_builder(5, 1);
    // graph.evaluate_cff_expression(&sample, &default_settings);

    // let eval: DerivedGraphData<PythonState> = graph.derived_data.unwrap().forget_type();

    // let v = serde_json::to_string(&eval).unwrap();
    // let u: DerivedGraphData<PythonState> = serde_json::from_str(&v).unwrap();
}

#[test]
#[allow(non_snake_case)]
fn pytest_physical_2L_6photons() {
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

#[test]
#[allow(non_snake_case)]
fn pytest_top_bubble() {
    init();
    let default_settings = load_default_settings();
    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_top_bubble/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    // println!("{:?}", graph.bare_graph.loop_momentum_basis);
    let sample = kinematics_builder(5, 1, &graph.bare_graph);

    graph.generate_cff();
    let mut export_settings = test_export_settings();

    export_settings.gammaloop_compile_options.inline_asm = true;

    println!(
        "{}",
        graph
            .clone()
            .apply_feynman_rules()
            .derived_data
            .unwrap()
            .numerator
            .export()
    );

    export_settings.numerator_settings = NumeratorEvaluatorOptions::Iterative {
        eval_options: EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::Compiled,
        },
        iterations: 1,
        n_cores: 1,
        verbose: false,
    };

    let mut graph_iterative_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );

    export_settings.numerator_settings = NumeratorEvaluatorOptions::Joint(EvaluatorOptions {
        cpe_rounds: Some(1),
        compile_options: NumeratorCompileOptions::Compiled,
    });

    let mut graph_joint_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );
    export_settings.numerator_settings = NumeratorEvaluatorOptions::Single(EvaluatorOptions {
        cpe_rounds: Some(1),
        compile_options: NumeratorCompileOptions::Compiled,
    });

    let mut graph_single_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );
    export_settings.numerator_settings = NumeratorEvaluatorOptions::Iterative {
        eval_options: EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::NotCompiled,
        },
        iterations: 1,
        n_cores: 1,
        verbose: false,
    };
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

    println!(
        "Eval single{}",
        graph_single.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval joint{}",
        graph_joint.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval iterative{}",
        graph_iterative.evaluate_cff_expression(&sample, &default_settings)
    );

    println!(
        "Eval single compiled{}",
        graph_single_compiled.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval joint compiled{}",
        graph_joint_compiled.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval iterative compiled{}",
        graph_iterative_compiled.evaluate_cff_expression(&sample, &default_settings)
    );

    let a = &graph_iterative_compiled
        .derived_data
        .as_ref()
        .unwrap()
        .numerator
        .state
        .double_param_values;

    for p in a {
        println!("{:?}", p);
    }
}

#[test]
#[allow(non_snake_case)]
fn pytest_physical_1L_6photons_generate() {
    init();
    let mut default_settings = load_default_settings();
    let (model, amplitude, path) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    // println!("{:?}", graph.bare_graph.loop_momentum_basis);
    let sample = kinematics_builder(5, 1, &graph.bare_graph);

    graph.generate_cff();
    let mut export_settings = test_export_settings();

    export_settings.gammaloop_compile_options.inline_asm = true;

    println!(
        "{}",
        graph
            .clone()
            .apply_feynman_rules()
            .derived_data
            .unwrap()
            .numerator
            .export()
    );

    export_settings.numerator_settings = NumeratorEvaluatorOptions::Iterative {
        eval_options: EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::Compiled,
        },
        iterations: 1,
        n_cores: 1,
        verbose: false,
    };

    let mut graph_iterative_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );

    export_settings.numerator_settings = NumeratorEvaluatorOptions::Joint(EvaluatorOptions {
        cpe_rounds: Some(1),
        compile_options: NumeratorCompileOptions::Compiled,
    });

    let mut graph_joint_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );
    export_settings.numerator_settings = NumeratorEvaluatorOptions::Single(EvaluatorOptions {
        cpe_rounds: Some(1),
        compile_options: NumeratorCompileOptions::Compiled,
    });

    let mut graph_single_compiled = graph.clone().process_numerator(
        &model,
        ContractionSettings::Normal,
        path.clone(),
        &export_settings,
    );
    export_settings.numerator_settings = NumeratorEvaluatorOptions::Iterative {
        eval_options: EvaluatorOptions {
            cpe_rounds: Some(1),
            compile_options: NumeratorCompileOptions::NotCompiled,
        },
        iterations: 1,
        n_cores: 1,
        verbose: false,
    };
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

    default_settings.general.force_orientations = Some(vec![0]);
    println!(
        "Eval single{}",
        graph_single.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval joint{}",
        graph_joint.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval iterative{}",
        graph_iterative.evaluate_cff_expression(&sample, &default_settings)
    );

    default_settings.general.force_orientations = None;
    println!(
        "Eval single compiled{}",
        graph_single_compiled.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval joint compiled{}",
        graph_joint_compiled.evaluate_cff_expression(&sample, &default_settings)
    );
    println!(
        "Eval iterative compiled{}",
        graph_iterative_compiled.evaluate_cff_expression(&sample, &default_settings)
    );
}

#[test]
fn gamma_symplify() {
    // init();
    let (_, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    // graph.bare_graph.dot_internal_vertices()

    let num = graph.apply_feynman_rules().derived_data.unwrap().numerator;

    let mut a = num.color_symplify().state.expression;

    println!("Before: {}", a);

    a.0 = a.0.expand();

    let pats = [(
        Pattern::parse("id(aind(a_,b_))*t_(aind(d___,b_,c___))").unwrap(),
        Pattern::parse("t_(aind(d___,a_,c___))").unwrap().into(),
    )];
    let reps: Vec<Replacement> = pats
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect();
    a.replace_repeat_multiple(&reps);

    let pats = vec![
        (
            Pattern::parse("γ(aind(a_,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a_,d_,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a__,d__,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("γ(aind(a_,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a_,d__,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a__,d_,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,b_))").unwrap(),
            Pattern::parse("gamma_trace(aind(a__))").unwrap().into(),
        ),
    ];
    let reps: Vec<Replacement> = pats
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect();
    a.replace_repeat_multiple(&reps);

    let pat = Pattern::parse("gamma_trace(a__)").unwrap();

    let set = MatchSettings::default();
    let cond = Condition::default();

    let mut it = pat.pattern_match(a.0.as_view(), &cond, &set);

    let mut max_nargs = 0;
    while let Some(a) = it.next() {
        for (_, v) in a.match_stack {
            match v {
                Match::Single(s) => {
                    match s {
                        AtomView::Fun(f) => {
                            let a = f.get_nargs();
                            if a > max_nargs {
                                max_nargs = a;
                            }
                        }
                        _ => {
                            panic!("should be a function")
                        }
                    }
                    print!("{}", s)
                }
                _ => panic!("should be a single match"),
            }
            println!();
        }
    }

    let mut reps = vec![];
    for n in 1..=max_nargs {
        if n % 2 == 0 {
            let mut sum = Atom::new_num(0);

            // sum((-1)**(k+1) * d(p_[0], p_[k]) * f(*p_[1:k], *p_[k+1:l])
            for j in 1..n {
                let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));

                let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                let metric_builder = FunctionBuilder::new(State::get_symbol("g"));

                let metric_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                for k in 1..j {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }

                for k in (j + 1)..n {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }

                let metric = metric_builder
                    .add_arg(
                        &metric_builder_slots
                            .add_args(&[
                                &Atom::parse(&format!("a{}_", 0)).unwrap(),
                                &Atom::parse(&format!("a{}_", j)).unwrap(),
                            ])
                            .finish(),
                    )
                    .finish();

                let gamma = &gamma_chain_builder
                    .add_arg(&gamma_chain_builder_slots.finish())
                    .finish()
                    * &metric;

                if j % 2 == 0 {
                    sum = &sum - &gamma;
                } else {
                    sum = &sum + &gamma;
                }
            }

            let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
            let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
            for k in 0..n {
                let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
            }
            let a = gamma_chain_builder
                .add_arg(&gamma_chain_builder_slots.finish())
                .finish();

            reps.push((a.into_pattern(), sum.into_pattern().into()));
        } else {
            let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
            let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
            for k in 0..n {
                let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
            }
            let a = gamma_chain_builder
                .add_arg(&gamma_chain_builder_slots.finish())
                .finish();
            println!("{}", a);
            reps.push((a.into_pattern(), Atom::new_num(0).into_pattern().into()));
        }
    }

    reps.push((
        Pattern::parse("gamma_trace(aind())").unwrap(),
        Pattern::parse("1").unwrap().into(),
    ));

    // Dd
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(lord(a__),lord(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));
    // Du
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(lord(a__),loru(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // Uu
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(loru(a__),loru(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // Ud
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(loru(a__),lord(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));

    // dD
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(lord(b__),lord(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));
    // uD
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(loru(b__),lord(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // uU
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(loru(b__),loru(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // dU
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(lord(b__),loru(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));

    let reps = reps
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect_vec();

    a.replace_repeat_multiple(&reps);
    a.0 = a.0.expand();

    a.replace_repeat_multiple(&reps);
    a.0 = a.0.factor();
    println!("After: {}", a);
    let tn = TensorNetwork::try_from(a.0.as_view())
        .unwrap()
        .to_fully_parametric();
    println!("After parse: {}", tn.dot());
}

#[test]
fn parse() {
    init();
    let expr = Atom::parse("64/729*ee^6*MT^6*Nc*ϵ(0,aind(lord(4,7)))*ϵ(1,aind(loru(4,7)))*ϵbar(2,aind(lord(4,8)))*ϵbar(3,aind(lord(4,9)))*ϵbar(4,aind(loru(4,9)))*ϵbar(5,aind(loru(4,8)))").unwrap();

    let mut net = TensorNetwork::try_from(expr.as_view()).unwrap();
    println!("{}", net.dot());

    net.contract();

    println!("{}", net.dot());

    println!("{}", &net.graph.involution.len());

    for (ni, i) in &net.graph.involution {
        if &ni != i {
            println!("internal");
        } else {
            println!("external");
        }
    }

    // println!("{}", net.to_fully_parametric().result_tensor().unwrap());
}

#[test]
fn tsgstg() {
    let _ = load_default_settings();
    let (_, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    // graph.bare_graph.dot_internal_vertices()

    let num = graph.apply_feynman_rules().derived_data.unwrap().numerator;

    let mut a = num.color_symplify().state.expression;

    println!("Before: {}", a);

    a.0 = a.0.expand();

    let number_terms = if let Atom::Add(a) = &a.0 {
        a.get_nargs()
    } else {
        1
    };
    println!("Number of terms: {}", number_terms);

    let pats = [(
        Pattern::parse("id(aind(a_,b_))*t_(aind(d___,b_,c___))").unwrap(),
        Pattern::parse("t_(aind(d___,a_,c___))").unwrap().into(),
    )];
    let reps: Vec<Replacement> = pats
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect();
    a.replace_repeat_multiple(&reps);

    let pats = vec![
        (
            Pattern::parse("γ(aind(a_,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a_,d_,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a__,d__,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("γ(aind(a_,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a_,d__,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a__,d_,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,b_))").unwrap(),
            Pattern::parse("gamma_trace(aind(a__))").unwrap().into(),
        ),
    ];
    let reps: Vec<Replacement> = pats
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect();
    a.replace_repeat_multiple(&reps);

    println!("After preprocessing: {}", a);

    let pat = Pattern::parse("gamma_trace(a__)").unwrap();

    let set = MatchSettings::default();
    let cond = Condition::default();

    let mut it = pat.pattern_match(a.0.as_view(), &cond, &set);

    let mut max_nargs = 0;
    while let Some(a) = it.next() {
        for (_, v) in a.match_stack {
            match v {
                Match::Single(s) => {
                    match s {
                        AtomView::Fun(f) => {
                            let a = f.get_nargs();
                            if a > max_nargs {
                                max_nargs = a;
                            }
                        }
                        _ => {
                            panic!("should be a function")
                        }
                    }
                    // print!("{}", s)
                }
                _ => panic!("should be a single match"),
            }
            // println!();
        }
    }

    let mut reps: Vec<(Pattern, PatternOrMap)> = vec![];
    for n in 1..=max_nargs {
        if n % 2 == 0 {
            let mut sum = Atom::new_num(0);

            // sum((-1)**(k+1) * d(p_[0], p_[k]) * f(*p_[1:k], *p_[k+1:l])
            for j in 1..n {
                let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));

                let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                let metric_builder = FunctionBuilder::new(State::get_symbol("g"));

                let metric_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                for k in 1..j {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }

                for k in (j + 1)..n {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }

                let metric = metric_builder
                    .add_arg(
                        &metric_builder_slots
                            .add_args(&[
                                &Atom::parse(&format!("a{}_", 0)).unwrap(),
                                &Atom::parse(&format!("a{}_", j)).unwrap(),
                            ])
                            .finish(),
                    )
                    .finish();

                let gamma = &gamma_chain_builder
                    .add_arg(&gamma_chain_builder_slots.finish())
                    .finish()
                    * &metric;

                if j % 2 == 0 {
                    sum = &sum - &gamma;
                } else {
                    sum = &sum + &gamma;
                }
            }

            let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
            let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
            for k in 0..n {
                let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
            }
            let a = gamma_chain_builder
                .add_arg(&gamma_chain_builder_slots.finish())
                .finish();

            reps.push((a.into_pattern(), sum.into_pattern().into()));
        } else {
            let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
            let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
            for k in 0..n {
                let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
            }
            let a = gamma_chain_builder
                .add_arg(&gamma_chain_builder_slots.finish())
                .finish();
            // println!("{}", a);
            reps.push((a.into_pattern(), Atom::new_num(0).into_pattern().into()));
        }
    }

    reps.push((
        Pattern::parse("gamma_trace(aind())").unwrap(),
        Pattern::parse("1").unwrap().into(),
    ));

    // Dd
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(lord(a__),lord(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));
    // Du
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(lord(a__),loru(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // Uu
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(loru(a__),loru(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // Ud
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(loru(a__),lord(b__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));

    // dD
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(lord(b__),lord(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));
    // uD
    reps.push((
        Pattern::parse("f_(i___,aind(o___,loru(a__),e___))*g(aind(loru(b__),lord(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // uU
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(loru(b__),loru(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,loru(b__),e___))")
            .unwrap()
            .into(),
    ));
    // dU
    reps.push((
        Pattern::parse("f_(i___,aind(o___,lord(a__),e___))*g(aind(lord(b__),loru(a__)))").unwrap(),
        Pattern::parse("f_(i___,aind(o___,lord(b__),e___))")
            .unwrap()
            .into(),
    ));

    let reps = reps
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect_vec();

    let time = Instant::now();
    a.replace_repeat_multiple(&reps);
    a.0 = a.0.expand();

    a.replace_repeat_multiple(&reps);

    println!("Time: {}", time.elapsed().as_secs_f64());

    let number_terms = if let Atom::Add(a) = &a.0 {
        a.get_nargs()
    } else {
        1
    };
    println!("Number of terms: {}", number_terms);

    // println!("After: {}", a);
}

#[test]
fn tsgstrstg() {
    let _ = load_default_settings();
    let (_, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);

    let graph = amplitude.amplitude_graphs[0].graph.clone();

    // graph.bare_graph.dot_internal_vertices()

    let num = graph.apply_feynman_rules().derived_data.unwrap().numerator;

    let mut a = num.color_symplify().state.expression;

    println!("Before: {}", a);

    a.0 = a.0.expand();

    let number_terms = if let Atom::Add(a) = &a.0 {
        a.get_nargs()
    } else {
        1
    };
    println!("Number of terms: {}", number_terms);

    let pats = [(
        Pattern::parse("id(aind(a_,b_))*t_(aind(d___,b_,c___))").unwrap(),
        Pattern::parse("t_(aind(d___,a_,c___))").unwrap().into(),
    )];
    let reps: Vec<Replacement> = pats
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect();
    a.replace_repeat_multiple(&reps);

    let pats = vec![
        (
            Pattern::parse("γ(aind(a_,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a_,d_,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a__,d__,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("γ(aind(a_,b_,c_))*gamma_chain(aind(d__,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a_,d__,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,c_))*γ(aind(d_,c_,e_))").unwrap(),
            Pattern::parse("gamma_chain(aind(a__,d_,b_,e_))")
                .unwrap()
                .into(),
        ),
        (
            Pattern::parse("gamma_chain(aind(a__,b_,b_))").unwrap(),
            Pattern::parse("gamma_trace(aind(a__))").unwrap().into(),
        ),
    ];
    let reps: Vec<Replacement> = pats
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect();
    a.replace_repeat_multiple(&reps);

    println!("After preprocessing: {}", a);

    let pat = Pattern::parse("gamma_trace(a__)").unwrap();

    let set = MatchSettings::default();
    let cond = Condition::default();

    let mut it = pat.pattern_match(a.0.as_view(), &cond, &set);

    let mut max_nargs = 0;
    while let Some(a) = it.next() {
        for (_, v) in a.match_stack {
            match v {
                Match::Single(s) => {
                    match s {
                        AtomView::Fun(f) => {
                            let a = f.get_nargs();
                            if a > max_nargs {
                                max_nargs = a;
                            }
                        }
                        _ => {
                            panic!("should be a function")
                        }
                    }
                    // print!("{}", s)
                }
                _ => panic!("should be a single match"),
            }
            // println!();
        }
    }

    let mut reps = vec![];
    for n in 1..=max_nargs {
        if n % 2 == 0 {
            let mut sum = Atom::new_num(0);

            // sum((-1)**(k+1) * d(p_[0], p_[k]) * f(*p_[1:k], *p_[k+1:l])
            for j in 1..n {
                let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));

                let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                let metric_builder = FunctionBuilder::new(State::get_symbol("g"));

                let metric_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));

                for k in 1..j {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }

                for k in (j + 1)..n {
                    let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                    gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
                }

                let metric = metric_builder
                    .add_arg(
                        &metric_builder_slots
                            .add_args(&[
                                &Atom::parse(&format!("a{}_", 0)).unwrap(),
                                &Atom::parse(&format!("a{}_", j)).unwrap(),
                            ])
                            .finish(),
                    )
                    .finish();

                let gamma = &gamma_chain_builder
                    .add_arg(&gamma_chain_builder_slots.finish())
                    .finish()
                    * &metric;

                if j % 2 == 0 {
                    sum = &sum - &gamma;
                } else {
                    sum = &sum + &gamma;
                }
            }

            let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
            let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
            for k in 0..n {
                let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
            }
            let a = gamma_chain_builder
                .add_arg(&gamma_chain_builder_slots.finish())
                .finish();

            reps.push((a.into_pattern(), sum.into_pattern().into()));
        } else {
            let gamma_chain_builder = FunctionBuilder::new(State::get_symbol("gamma_trace"));
            let mut gamma_chain_builder_slots = FunctionBuilder::new(State::get_symbol("aind"));
            for k in 0..n {
                let mu = Atom::parse(&format!("a{}_", k)).unwrap();
                gamma_chain_builder_slots = gamma_chain_builder_slots.add_arg(&mu);
            }
            let a = gamma_chain_builder
                .add_arg(&gamma_chain_builder_slots.finish())
                .finish();
            // println!("{}", a);
            reps.push((a.into_pattern(), Atom::new_num(0).into_pattern().into()));
        }
    }

    reps.push((
        Pattern::parse("gamma_trace(aind())").unwrap(),
        Pattern::parse("1").unwrap().into(),
    ));

    // Dd
    reps.push((
        Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(lord(a__),lord(b__)))").unwrap(),
        Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
    ));
    // Du
    reps.push((
        Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(lord(a__),loru(b__)))").unwrap(),
        Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
    ));
    // Uu
    reps.push((
        Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(loru(a__),loru(b__)))").unwrap(),
        Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
    ));
    // Ud
    reps.push((
        Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(loru(a__),lord(b__)))").unwrap(),
        Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
    ));

    // dD
    reps.push((
        Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(lord(b__),lord(a__)))").unwrap(),
        Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
    ));
    // uD
    reps.push((
        Pattern::parse("f_(i_,aind(loru(a__)))*g(aind(loru(b__),lord(a__)))").unwrap(),
        Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
    ));
    // uU
    reps.push((
        Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(loru(b__),loru(a__)))").unwrap(),
        Pattern::parse("f_(i_,aind(loru(b__)))").unwrap().into(),
    ));
    // dU
    reps.push((
        Pattern::parse("f_(i_,aind(lord(a__)))*g(aind(lord(b__),loru(a__)))").unwrap(),
        Pattern::parse("f_(i_,aind(lord(b__)))").unwrap().into(),
    ));

    let reps = reps
        .iter()
        .map(|(lhs, rhs)| Replacement::new(lhs, rhs))
        .collect_vec();

    let time = Instant::now();
    a.replace_repeat_multiple(&reps);
    a.0 = a.0.expand();

    a.replace_repeat_multiple(&reps);

    println!("Time: {}", time.elapsed().as_secs_f64());

    let number_terms = if let Atom::Add(a) = &a.0 {
        a.get_nargs()
    } else {
        1
    };
    println!("Number of terms: {}", number_terms);

    // println!("After: {}", a);
}
