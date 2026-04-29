use std::collections::{BTreeMap, BTreeSet};
use std::fs;

use color_eyre::eyre::eyre;
use gammalooprs::{graph::Graph, processes::ProcessCollection};
use three_dimensional_reps::{
    EnergyEdgeIndexMap, Generate3DExpressionOptions, HybridSurfaceID, LinearSurfaceKind,
    NumeratorSamplingScaleMode, OrientationID, ParsedGraph, RepresentationMode, ThreeDExpression,
    ThreeDGraphSource, auto_numerator_expr_for_bounds,
    eval::{EvaluationInput, evaluate_expression},
    extract_signatures_and_masses_from_symbolica_expression, generate_3d_expression_from_parsed,
    generation::GenerationError,
    graph_signatures::MomentumSignature,
};

use super::*;

fn gammaloop_threedreps_dot_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("resources/graphs/threedreps")
        .join(name)
}

#[derive(Clone, Copy)]
struct ImportedGraphSpec {
    dot_name: &'static str,
    process_name: &'static str,
    internal_edges: usize,
    loop_count: usize,
    external_symbols: usize,
    external_edges: usize,
    repeated_group_sizes: &'static [usize],
}

const REPEATED_EXAMPLES: &[&str] = &[
    "box_pow3.dot",
    "sunrise_pow4.dot",
    "kite_double_nested_repeats.dot",
    "kite_nested_repeats.dot",
    "kite_sandwich_repeats.dot",
    "proper_iterated_sandwiched_bubble.dot",
    "mercedes_multi_repeats.dot",
];

const IMPORTED_GRAPH_SPECS: &[ImportedGraphSpec] = &[
    ImportedGraphSpec {
        dot_name: "box.dot",
        process_name: "threedreps_box_shape",
        internal_edges: 4,
        loop_count: 1,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[],
    },
    ImportedGraphSpec {
        dot_name: "box_pow3.dot",
        process_name: "threedreps_box_pow3_shape",
        internal_edges: 6,
        loop_count: 1,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[3],
    },
    ImportedGraphSpec {
        dot_name: "sunrise_pow4.dot",
        process_name: "threedreps_sunrise_pow4_shape",
        internal_edges: 6,
        loop_count: 2,
        external_symbols: 2,
        external_edges: 2,
        repeated_group_sizes: &[4],
    },
    ImportedGraphSpec {
        dot_name: "kite_double_nested_repeats.dot",
        process_name: "threedreps_kite_double_nested_shape",
        internal_edges: 9,
        loop_count: 2,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[2, 2, 2, 2],
    },
    ImportedGraphSpec {
        dot_name: "kite_nested_repeats.dot",
        process_name: "threedreps_kite_nested_shape",
        internal_edges: 8,
        loop_count: 2,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[2, 2, 2],
    },
    ImportedGraphSpec {
        dot_name: "kite_sandwich_repeats.dot",
        process_name: "threedreps_kite_sandwich_shape",
        internal_edges: 7,
        loop_count: 2,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[2, 2],
    },
    ImportedGraphSpec {
        dot_name: "proper_iterated_sandwiched_bubble.dot",
        process_name: "threedreps_iterated_sandwiched_shape",
        internal_edges: 8,
        loop_count: 3,
        external_symbols: 2,
        external_edges: 2,
        repeated_group_sizes: &[2, 2],
    },
    ImportedGraphSpec {
        dot_name: "mercedes_multi_repeats.dot",
        process_name: "threedreps_mercedes_multi_shape",
        internal_edges: 10,
        loop_count: 3,
        external_symbols: 0,
        external_edges: 0,
        repeated_group_sizes: &[2, 2, 3],
    },
    ImportedGraphSpec {
        dot_name: "four_loop_stress.dot",
        process_name: "threedreps_four_loop_stress_shape",
        internal_edges: 10,
        loop_count: 4,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[3],
    },
    ImportedGraphSpec {
        dot_name: "five_loop_no_repeats.dot",
        process_name: "threedreps_five_loop_no_repeats_shape",
        internal_edges: 13,
        loop_count: 5,
        external_symbols: 4,
        external_edges: 4,
        repeated_group_sizes: &[],
    },
    ImportedGraphSpec {
        dot_name: "one_loop_10_external.dot",
        process_name: "threedreps_one_loop_10_external_shape",
        internal_edges: 10,
        loop_count: 1,
        external_symbols: 10,
        external_edges: 10,
        repeated_group_sizes: &[],
    },
    ImportedGraphSpec {
        dot_name: "one_loop_15_external.dot",
        process_name: "threedreps_one_loop_15_external_shape",
        internal_edges: 15,
        loop_count: 1,
        external_symbols: 15,
        external_edges: 15,
        repeated_group_sizes: &[],
    },
];

#[derive(Clone, Copy)]
struct ProbeCase {
    name: &'static str,
    dot_name: &'static str,
    numerator: &'static str,
    bounds: &'static [(usize, usize)],
    seed: u64,
    compare: ProbeComparison,
    sampling_scale: NumeratorSamplingScaleMode,
    uniform_scale: Option<f64>,
}

#[derive(Clone, Copy)]
enum ProbeComparison {
    CffLtd { tolerance: f64 },
    CffPureLtd { tolerance: f64 },
    BuildOnly,
}

#[derive(serde::Serialize)]
struct ProbeCaseReport {
    name: String,
    dot_name: String,
    numerator: String,
    bounds: Vec<(usize, usize)>,
    comparison: String,
    sampling_scale: String,
    status: String,
    detail: String,
}

const BOUNDS_EMPTY: &[(usize, usize)] = &[];
const BOUNDS_BOX_QUADRATIC_0: &[(usize, usize)] = &[(0, 2)];
const BOUNDS_BOX_SPLIT_SPECTATORS: &[(usize, usize)] = &[(0, 1), (1, 1), (2, 1), (3, 2)];
const BOUNDS_BOX_CUBIC_0: &[(usize, usize)] = &[(0, 3)];
const BOUNDS_BOX_CUBIC_0_QUADRATIC_2: &[(usize, usize)] = &[(0, 3), (2, 2)];
const BOUNDS_BOX_HIGH_CONTACT_A: &[(usize, usize)] = &[(0, 4)];
const BOUNDS_BOX_HIGH_CONTACT_B: &[(usize, usize)] = &[(0, 3), (1, 2)];
const BOUNDS_BOX_HIGH_CONTACT_C: &[(usize, usize)] = &[(0, 2), (1, 1), (3, 3)];
const BOUNDS_BOX_HIGH_CONTACT_D: &[(usize, usize)] = &[(2, 3), (3, 3)];
const BOUNDS_BOX_ALL_AFFINE: &[(usize, usize)] = &[(0, 1), (1, 1), (2, 1), (3, 1)];
const BOUNDS_BOX_POW3_QUADRATIC: &[(usize, usize)] = &[(0, 2), (1, 1), (2, 1), (3, 2)];
const BOUNDS_BOX_POW3_CUBIC_SPECTATOR: &[(usize, usize)] = &[(0, 3), (1, 1), (2, 1), (3, 1)];
const BOUNDS_BOX_POW3_REPEATED_HIGH: &[(usize, usize)] = &[(0, 1), (1, 1), (3, 4)];
const BOUNDS_BOX_POW3_REPEATED_QUINTIC: &[(usize, usize)] = &[(3, 5)];
const BOUNDS_BOX_POW3_MIXED_QUINTIC: &[(usize, usize)] = &[(0, 2), (1, 1), (2, 1), (3, 5)];
const BOUNDS_BOX_POW3_TWO_REPEATED_HIGH: &[(usize, usize)] = &[(3, 3), (4, 4)];
const BOUNDS_BOX_POW3_AMBIGUOUS: &[(usize, usize)] = &[(0, 1), (1, 2), (2, 2), (3, 2)];
const BOUNDS_SUNRISE_QUADRATIC_0: &[(usize, usize)] = &[(0, 2)];
const BOUNDS_SUNRISE_QUINTIC_2: &[(usize, usize)] = &[(2, 5)];
const BOUNDS_SUNRISE_CUBIC_PAIR: &[(usize, usize)] = &[(2, 3), (3, 3)];
const BOUNDS_ITER_QUADRATIC_COMBO: &[(usize, usize)] = &[(0, 2), (1, 2), (2, 2)];
const BOUNDS_ITER_CUBIC_LINEAR: &[(usize, usize)] = &[(0, 3), (1, 1)];
const BOUNDS_ITER_CUBIC_PAIR: &[(usize, usize)] = &[(1, 3), (3, 3)];
const BOUNDS_FOUR_LOOP_QUADRATIC: &[(usize, usize)] = &[(0, 2), (1, 2)];
const BOUNDS_KITE_SANDWICH_HIGH: &[(usize, usize)] = &[(0, 4), (5, 3)];
const BOUNDS_KITE_NESTED_HIGH: &[(usize, usize)] = &[(0, 5)];

fn orientation_labels(expression: &ThreeDExpression<OrientationID>) -> Vec<String> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            orientation
                .data
                .label
                .clone()
                .expect("orientation label should be set")
        })
        .collect()
}

fn compact_orientation_labels(
    expression: &ThreeDExpression<OrientationID>,
    source_internal_edges: &[usize],
) -> Vec<String> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            source_internal_edges
                .iter()
                .map(|edge_id| {
                    let orientation = orientation
                        .data
                        .orientation
                        .iter()
                        .find_map(|(candidate_edge_id, orientation)| {
                            (candidate_edge_id.0 == *edge_id).then_some(orientation)
                        })
                        .expect("source edge id should be present in orientation");
                    match format!("{orientation:?}").as_str() {
                        "Default" => '+',
                        "Reversed" => '-',
                        "Undirected" => 'x',
                        other => panic!("unexpected orientation value {other}"),
                    }
                })
                .collect()
        })
        .collect()
}

fn compact_base_orientation_labels(
    expression: &ThreeDExpression<OrientationID>,
    source_internal_edges: &[usize],
) -> Vec<String> {
    compact_orientation_labels(expression, source_internal_edges)
        .into_iter()
        .map(|label| {
            label
                .split_once('|')
                .map_or(label.clone(), |(base, _)| base.to_string())
        })
        .collect()
}

fn base_orientation_labels(expression: &ThreeDExpression<OrientationID>) -> Vec<String> {
    orientation_labels(expression)
        .into_iter()
        .map(|label| {
            label
                .split_once('|')
                .map_or(label.clone(), |(base, _)| base.to_string())
        })
        .collect()
}

fn nontrivial_sign_labels(n_edges: usize) -> Vec<String> {
    (1..((1usize << n_edges) - 1))
        .map(|bitmask| {
            (0..n_edges)
                .map(|edge_id| {
                    if bitmask & (1usize << edge_id) == 0 {
                        '+'
                    } else {
                        '-'
                    }
                })
                .collect()
        })
        .collect()
}

fn assert_internal_signatures(parsed: &ParsedGraph, expected: &[(&str, [i32; 4])]) {
    let signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| {
            (
                edge.label.as_str(),
                edge.signature.loop_signature.as_slice(),
                edge.signature.external_signature.as_slice(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        signatures.len(),
        expected.len(),
        "unexpected number of internal signatures: {signatures:?}"
    );
    for ((label, loop_signature, external_signature), (expected_label, expected_external)) in
        signatures.iter().zip(expected)
    {
        assert_eq!(*label, *expected_label);
        assert_eq!(*loop_signature, &[1]);
        assert_eq!(*external_signature, expected_external);
    }
}

fn sorted_signature_pairs(signatures: impl IntoIterator<Item = MomentumSignature>) -> Vec<String> {
    let mut items = signatures
        .into_iter()
        .map(|signature| {
            format!(
                "{:?}|{:?}",
                signature.loop_signature, signature.external_signature
            )
        })
        .collect::<Vec<_>>();
    items.sort();
    items
}

fn project_external_signatures_by_name(
    parsed: &ParsedGraph,
    external_names: &[String],
) -> Vec<MomentumSignature> {
    parsed
        .internal_edges
        .iter()
        .map(|edge| {
            let mut external_signature = vec![0; external_names.len()];
            for (source_external_id, coeff) in edge.signature.external_signature.iter().enumerate()
            {
                let name = &parsed.external_names[source_external_id];
                let target_external_id = external_names
                    .iter()
                    .position(|candidate| candidate == name)
                    .unwrap_or_else(|| panic!("unexpected external name {name}"));
                external_signature[target_external_id] += coeff;
            }
            MomentumSignature {
                loop_signature: edge.signature.loop_signature.clone(),
                external_signature,
            }
        })
        .collect()
}

fn linear_denominator_surface_ids(expression: &ThreeDExpression<OrientationID>) -> Vec<Vec<usize>> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            assert_eq!(orientation.variants.len(), 1);
            orientation.variants[0]
                .denominator
                .iter_nodes()
                .filter_map(|node| match node.data {
                    HybridSurfaceID::Linear(id)
                        if expression.surfaces.linear_surface_cache[id].kind
                            == LinearSurfaceKind::Esurface =>
                    {
                        Some(id.0)
                    }
                    HybridSurfaceID::Linear(_) => None,
                    HybridSurfaceID::Unit => None,
                    other => panic!("unexpected surface id in denominator tree: {other:?}"),
                })
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect()
        })
        .collect()
}

fn variant_linear_denominator_surface_ids(
    expression: &ThreeDExpression<OrientationID>,
) -> Vec<Vec<Vec<usize>>> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            orientation
                .variants
                .iter()
                .map(|variant| {
                    variant
                        .denominator
                        .iter_nodes()
                        .filter_map(|node| match node.data {
                            HybridSurfaceID::Linear(id)
                                if expression.surfaces.linear_surface_cache[id].kind
                                    == LinearSurfaceKind::Esurface =>
                            {
                                Some(id.0)
                            }
                            HybridSurfaceID::Linear(_) => None,
                            HybridSurfaceID::Unit => None,
                            other => panic!("unexpected surface id in denominator tree: {other:?}"),
                        })
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect()
                })
                .collect()
        })
        .collect()
}

fn half_edges(expression: &ThreeDExpression<OrientationID>) -> Vec<Vec<usize>> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            assert_eq!(orientation.variants.len(), 1);
            orientation.variants[0]
                .half_edges
                .iter()
                .map(|edge| edge.0)
                .collect()
        })
        .collect()
}

fn variant_half_edges(expression: &ThreeDExpression<OrientationID>) -> Vec<Vec<Vec<usize>>> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            orientation
                .variants
                .iter()
                .map(|variant| variant.half_edges.iter().map(|edge| edge.0).collect())
                .collect()
        })
        .collect()
}

fn variant_prefactors(expression: &ThreeDExpression<OrientationID>) -> Vec<Vec<(i64, i64)>> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            orientation
                .variants
                .iter()
                .map(|variant| {
                    variant
                        .prefactor
                        .to_i64_pair()
                        .expect("test prefactors should be rational i64 pairs")
                })
                .collect()
        })
        .collect()
}

fn assert_unit_or_minus_one_prefactors(
    expression: &ThreeDExpression<OrientationID>,
    expected: (i64, i64),
) -> Result<()> {
    for orientation in &expression.orientations {
        if orientation.variants.len() != 1 {
            return Err(eyre!(
                "orientation {} has {} variants; expected 1",
                orientation.data.label.as_deref().unwrap_or("<unlabeled>"),
                orientation.variants.len()
            ));
        }
        if orientation.variants[0].prefactor.to_i64_pair() != Some(expected) {
            return Err(eyre!(
                "orientation {} has prefactor {}; expected {}/{}",
                orientation.data.label.as_deref().unwrap_or("<unlabeled>"),
                orientation.variants[0].prefactor,
                expected.0,
                expected.1
            ));
        }
    }
    Ok(())
}

fn python_reference_input(masses: Vec<f64>) -> EvaluationInput {
    let p1 = [
        0.16485399973205894,
        0.019959344163000492,
        -0.08049098445037467,
        0.051472412341362905,
    ];
    let p2 = [
        -0.46803780283570684,
        0.19462424814457885,
        -0.0697771148314748,
        0.1737676949494124,
    ];
    let p3 = [
        0.5904171138547634,
        -0.11541996721229852,
        0.2953482435514923,
        -0.17703239875345095,
    ];
    let p4 = [
        -p1[0] - p2[0] - p3[0],
        -p1[1] - p2[1] - p3[1],
        -p1[2] - p2[2] - p3[2],
        -p1[3] - p2[3] - p3[3],
    ];
    EvaluationInput {
        external_momenta: vec![p1, p2, p3, p4],
        loop_spatial_momenta: vec![[
            0.07813881659444444,
            0.20615267144468802,
            -0.1953558011893643,
        ]],
        masses,
        uniform_scale: None,
    }
}

fn repeated_box_masses_with_split(epsilon: f64) -> Vec<f64> {
    vec![
        1.138394,
        0.758865,
        0.458809,
        1.029287 - epsilon,
        1.029287,
        1.029287 + epsilon,
    ]
}

fn deterministic_input(
    parsed: &ParsedGraph,
    seed: u64,
    uniform_scale: Option<f64>,
) -> Result<EvaluationInput> {
    Ok(EvaluationInput::deterministic(
        parsed,
        seed,
        &BTreeMap::new(),
        uniform_scale,
    )?)
}

fn expression_options(
    representation: RepresentationMode,
    bounds: &[(usize, usize)],
    sampling_scale: NumeratorSamplingScaleMode,
) -> Generate3DExpressionOptions {
    Generate3DExpressionOptions {
        representation,
        energy_degree_bounds: bounds.to_vec(),
        numerator_sampling_scale: sampling_scale,
    }
}

fn compare_expression_pair(
    parsed: &ParsedGraph,
    lhs: &ThreeDExpression<OrientationID>,
    rhs: &ThreeDExpression<OrientationID>,
    numerator: &str,
    input: &EvaluationInput,
    tolerance: f64,
) -> Result<String> {
    let lhs_value = evaluate_expression(parsed, lhs, numerator, input)?.value;
    let rhs_value = evaluate_expression(parsed, rhs, numerator, input)?.value;
    let diff = (lhs_value - rhs_value).abs();
    if diff <= tolerance {
        Ok(format!(
            "ok: lhs={lhs_value:.17e}, rhs={rhs_value:.17e}, diff={diff:.3e}"
        ))
    } else {
        Err(eyre!(
            "mismatch: lhs={lhs_value:.17e}, rhs={rhs_value:.17e}, diff={diff:.3e}, tolerance={tolerance:.3e}"
        ))
    }
}

fn source_bounds_from_local(source_internal_edges: &[usize], bounds: &[(usize, usize)]) -> String {
    bounds
        .iter()
        .map(|(local_edge_id, degree)| {
            format!("{}:{degree}", source_internal_edges[*local_edge_id])
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn edge_energy_power_numerator(
    source_internal_edges: &[usize],
    bounds: &[(usize, usize)],
) -> String {
    let factors = bounds
        .iter()
        .filter(|(_, degree)| *degree > 0)
        .map(|(local_edge_id, degree)| {
            let edge_id = source_internal_edges[*local_edge_id];
            if *degree == 1 {
                format!("Q({edge_id},spenso::cind(0))")
            } else {
                format!("Q({edge_id},spenso::cind(0))^{degree}")
            }
        })
        .collect::<Vec<_>>();
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

fn loop_energy_power_numerator(loop_id: usize, degree: usize) -> String {
    if degree == 1 {
        format!("K({loop_id},spenso::cind(0))")
    } else {
        format!("K({loop_id},spenso::cind(0))^{degree}")
    }
}

fn dot_with_global_numerator(test_name: &str, dot_name: &str, numerator: &str) -> Result<PathBuf> {
    let test_root = get_tests_workspace_path().join(test_name);
    fs::create_dir_all(&test_root)?;
    let source = fs::read_to_string(gammaloop_threedreps_dot_path(dot_name))?;
    let escaped = numerator.replace('\\', "\\\\").replace('"', "\\\"");
    let updated = source.replacen(
        "particle=\"scalar_1\"",
        &format!("num=\"{escaped}\" particle=\"scalar_1\""),
        1,
    );
    let path = test_root.join(format!(
        "{}_num.dot",
        dot_name.strip_suffix(".dot").unwrap_or(dot_name)
    ));
    fs::write(&path, updated)?;
    Ok(path)
}

fn imported_graph_threedrep_workspace(test_name: &str) -> PathBuf {
    get_tests_workspace_path()
        .join(test_name)
        .join("threed_workspace")
}

struct ImportedGraphThreeDRepCase {
    expression: ThreeDExpression<OrientationID>,
    compact_expression: ThreeDExpression<OrientationID>,
}

struct ImportedGraphThreeDRepRun {
    cli: gammaloop_integration_tests::CLIState,
    manifest: JsonValue,
    parsed: ParsedGraph,
    source_internal_edges: Vec<usize>,
    cff: ImportedGraphThreeDRepCase,
    ltd: ImportedGraphThreeDRepCase,
    pure_ltd: ImportedGraphThreeDRepCase,
}

fn run_imported_graph_threedrep_test(
    test_name: &str,
    dot_name: &str,
    process_name: &str,
) -> Result<ImportedGraphThreeDRepRun> {
    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!(
        "import graphs {} -p {process_name} -i default -o",
        gammaloop_threedreps_dot_path(dot_name).display()
    ))?;
    cli.run_command(&format!(
        "3Drep test-cff-ltd -p {process_name} -i default -g 0 --workspace-path {}",
        workspace_path.display()
    ))?;

    let manifest_path = workspace_path.join("test_cff_ltd_manifest.json");
    let manifest = serde_json::from_str::<JsonValue>(&fs::read_to_string(&manifest_path)?)?;
    let graph = imported_graph_from_cli(&cli, process_name, "default", 0)?;
    let parsed = graph.to_three_d_parsed_graph()?;
    let source_internal_edges = source_internal_edges(graph);
    let cff = load_imported_graph_case(&manifest, graph, &parsed, "cff_none")?;
    let ltd = load_imported_graph_case(&manifest, graph, &parsed, "ltd_none")?;
    let pure_ltd = load_imported_graph_case(&manifest, graph, &parsed, "pureltd_none")?;

    Ok(ImportedGraphThreeDRepRun {
        cli,
        manifest,
        parsed,
        source_internal_edges,
        cff,
        ltd,
        pure_ltd,
    })
}

fn import_threedrep_graph(
    test_name: &str,
    dot_path: &Path,
    process_name: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!(
        "import graphs {} -p {process_name} -i default -o",
        dot_path.display()
    ))?;
    Ok(cli)
}

fn imported_process_name(dot_name: &str) -> Result<&'static str> {
    IMPORTED_GRAPH_SPECS
        .iter()
        .find_map(|spec| (spec.dot_name == dot_name).then_some(spec.process_name))
        .ok_or_else(|| eyre!("no imported graph spec for {dot_name}"))
}

fn old_python_probe_cases() -> Vec<ProbeCase> {
    let mut cases = vec![
        ProbeCase {
            name: "normal_box_cff_ltd_constant",
            dot_name: "box.dot",
            numerator: "1",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-12 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_quadratic_edge_energy",
            dot_name: "box.dot",
            numerator: "edges[0][0]**2",
            bounds: BOUNDS_BOX_QUADRATIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cubic_edge_energy",
            dot_name: "box.dot",
            numerator: "edges[0][0]**3",
            bounds: BOUNDS_BOX_CUBIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_4_0_0_0",
            dot_name: "box.dot",
            numerator: "edges[0][0]**4",
            bounds: BOUNDS_BOX_HIGH_CONTACT_A,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_3_2_0_0",
            dot_name: "box.dot",
            numerator: "edges[0][0]**3 * edges[1][0]**2",
            bounds: BOUNDS_BOX_HIGH_CONTACT_B,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_2_1_0_3",
            dot_name: "box.dot",
            numerator: "edges[0][0]**2 * edges[1][0] * edges[3][0]**3",
            bounds: BOUNDS_BOX_HIGH_CONTACT_C,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_0_0_3_3",
            dot_name: "box.dot",
            numerator: "edges[2][0]**3 * edges[3][0]**3",
            bounds: BOUNDS_BOX_HIGH_CONTACT_D,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_auto_sparse_caps",
            dot_name: "box.dot",
            numerator: "dot(edges[0], ext[0]) * dot(edges[0], ext[1]) * dot(edges[0], ext[2]) * dot(edges[2], ext[2]) * dot(edges[2], ext[3])",
            bounds: BOUNDS_BOX_CUBIC_0_QUADRATIC_2,
            seed: 1337,
            compare: ProbeComparison::BuildOnly,
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_cff_ltd_constant",
            dot_name: "box_pow3.dot",
            numerator: "1",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_edge_dot_numerator",
            dot_name: "box_pow3.dot",
            numerator: "dot(edges[0], ext[0]) + dot(edges[3], ext[0])",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_repeated_quadratic_bounds",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0]**2 * edges[1][0] * edges[2][0] * edges[3][0]**2",
            bounds: BOUNDS_BOX_POW3_QUADRATIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_single_cubic_nonrepeated_with_repeated_spectators",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0]**3 * edges[1][0] * edges[2][0] * edges[3][0]",
            bounds: BOUNDS_BOX_POW3_CUBIC_SPECTATOR,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_unsplit_repeated_higher_bounds",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0] * edges[1][0] * edges[3][0]**4",
            bounds: BOUNDS_BOX_POW3_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_repeated_quintic",
            dot_name: "box_pow3.dot",
            numerator: "edges[3][0]**5",
            bounds: BOUNDS_BOX_POW3_REPEATED_QUINTIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_mixed_quintic",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0]**2 * edges[1][0] * edges[2][0] * edges[3][0]**5",
            bounds: BOUNDS_BOX_POW3_MIXED_QUINTIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_two_repeated_high_bounds",
            dot_name: "box_pow3.dot",
            numerator: "dot(edges[3], ext[0]) * dot(edges[3], ext[1]) * dot(edges[3], ext[2]) * dot(edges[4], ext[0]) * dot(edges[4], ext[1]) * dot(edges[4], ext[2]) * dot(edges[4], ext[3])",
            bounds: BOUNDS_BOX_POW3_TWO_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_affine_caps_force_bounded_ltd",
            dot_name: "box_pow3.dot",
            numerator: "dot(edges[0], ext[0]) * dot(edges[1], ext[1]) * dot(edges[2], ext[2]) * dot(edges[3], ext[3])",
            bounds: BOUNDS_BOX_ALL_AFFINE,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_ambiguous_orientation_labels",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0] * edges[1][0]**2 * edges[2][0]**2 * edges[3][0]**2",
            bounds: BOUNDS_BOX_POW3_AMBIGUOUS,
            seed: 1337,
            compare: ProbeComparison::BuildOnly,
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_uniform_beyond_quadratic_m_1",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0] * edges[1][0] * edges[3][0]**4",
            bounds: BOUNDS_BOX_POW3_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
            uniform_scale: Some(1.0),
        },
        ProbeCase {
            name: "box_pow3_uniform_beyond_quadratic_m_minus_2",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0] * edges[1][0] * edges[3][0]**4",
            bounds: BOUNDS_BOX_POW3_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
            uniform_scale: Some(-2.0),
        },
        ProbeCase {
            name: "box_pow3_uniform_all_m_2_75",
            dot_name: "box_pow3.dot",
            numerator: "edges[0][0]**2 * edges[1][0] * edges[2][0] * edges[3][0]**2",
            bounds: BOUNDS_BOX_POW3_QUADRATIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::All,
            uniform_scale: Some(2.75),
        },
        ProbeCase {
            name: "sunrise_pow4_constant",
            dot_name: "sunrise_pow4.dot",
            numerator: "1",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "sunrise_pow4_edge_dot_numerator",
            dot_name: "sunrise_pow4.dot",
            numerator: "dot(edges[0], ext[0]) + dot(edges[-1], ext[0])",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "sunrise_pow4_quadratic_multiloop",
            dot_name: "sunrise_pow4.dot",
            numerator: "edges[0][0]**2",
            bounds: BOUNDS_SUNRISE_QUADRATIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "sunrise_pow4_free_lower_sector_quintic",
            dot_name: "sunrise_pow4.dot",
            numerator: "edges[2][0]**5",
            bounds: BOUNDS_SUNRISE_QUINTIC_2,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "sunrise_pow4_repeated_channel_cubic_pair",
            dot_name: "sunrise_pow4.dot",
            numerator: "edges[2][0]**3 * edges[3][0]**3",
            bounds: BOUNDS_SUNRISE_CUBIC_PAIR,
            seed: 123,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "iterated_sandwiched_linear_numerator",
            dot_name: "proper_iterated_sandwiched_bubble.dot",
            numerator: "dot(edges[1], ext[0]) + dot(edges[4], ext[0])",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "iterated_sandwiched_quadratic_combo",
            dot_name: "proper_iterated_sandwiched_bubble.dot",
            numerator: "edges[0][0]**2 * edges[1][0]**2",
            bounds: BOUNDS_ITER_QUADRATIC_COMBO,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "iterated_sandwiched_cubic_linear",
            dot_name: "proper_iterated_sandwiched_bubble.dot",
            numerator: "edges[0][0]**3 * edges[1][0]",
            bounds: BOUNDS_ITER_CUBIC_LINEAR,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "iterated_sandwiched_cubic_pair",
            dot_name: "proper_iterated_sandwiched_bubble.dot",
            numerator: "edges[1][0]**3 * edges[3][0]**3",
            bounds: BOUNDS_ITER_CUBIC_PAIR,
            seed: 123,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "kite_double_nested_linear",
            dot_name: "kite_double_nested_repeats.dot",
            numerator: "dot(edges[0], ext[0]) + dot(edges[-1], ext[0])",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "kite_nested_linear",
            dot_name: "kite_nested_repeats.dot",
            numerator: "dot(edges[0], ext[0]) + dot(edges[-1], ext[0])",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "kite_nested_high_repeated_build",
            dot_name: "kite_nested_repeats.dot",
            numerator: "edges[0][0]**5",
            bounds: BOUNDS_KITE_NESTED_HIGH,
            seed: 1337,
            compare: ProbeComparison::BuildOnly,
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "kite_sandwich_linear",
            dot_name: "kite_sandwich_repeats.dot",
            numerator: "dot(edges[0], ext[0]) + dot(edges[-1], ext[0])",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "kite_sandwich_higher_power_combo",
            dot_name: "kite_sandwich_repeats.dot",
            numerator: "edges[0][0]**4 * edges[5][0]**3",
            bounds: BOUNDS_KITE_SANDWICH_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "mercedes_multi_repeats_energy_numerator",
            dot_name: "mercedes_multi_repeats.dot",
            numerator: "edges[0][0] + edges[3][0]",
            bounds: BOUNDS_EMPTY,
            seed: 7,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "four_loop_stress_linear",
            dot_name: "four_loop_stress.dot",
            numerator: "dot(edges[0], ext[0]) + dot(edges[5], ext[1])",
            bounds: BOUNDS_EMPTY,
            seed: 5,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "four_loop_stress_quadratic",
            dot_name: "four_loop_stress.dot",
            numerator: "edges[0][0]**2 * edges[1][0]**2",
            bounds: BOUNDS_FOUR_LOOP_QUADRATIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
    ];

    for &dot_name in REPEATED_EXAMPLES {
        cases.push(ProbeCase {
            name: dot_name,
            dot_name,
            numerator: "1",
            bounds: BOUNDS_EMPTY,
            seed: 1337,
            compare: ProbeComparison::CffPureLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        });
    }

    cases
}

fn local_edge_energy_monomial(bounds: &[(usize, usize)]) -> String {
    let numerator = bounds
        .iter()
        .filter(|(_, power)| *power > 0)
        .map(|(edge_id, power)| format!("edges[{edge_id}][0]**{power}"))
        .join(" * ");
    if numerator.is_empty() {
        "1".to_string()
    } else {
        numerator
    }
}

fn run_probe_case(cli: &gammaloop_integration_tests::CLIState, case: ProbeCase) -> ProbeCaseReport {
    run_probe_case_parts(
        cli,
        case.name.to_string(),
        case.dot_name,
        case.numerator.to_string(),
        case.bounds.to_vec(),
        case.seed,
        case.compare,
        case.sampling_scale,
        case.uniform_scale,
    )
}

fn run_probe_case_parts(
    cli: &gammaloop_integration_tests::CLIState,
    name: String,
    dot_name: &str,
    numerator: String,
    bounds: Vec<(usize, usize)>,
    seed: u64,
    compare: ProbeComparison,
    sampling_scale_mode: NumeratorSamplingScaleMode,
    uniform_scale: Option<f64>,
) -> ProbeCaseReport {
    let comparison = match compare {
        ProbeComparison::CffLtd { .. } => "cff_ltd",
        ProbeComparison::CffPureLtd { .. } => "cff_pureltd",
        ProbeComparison::BuildOnly => "build_only",
    }
    .to_string();
    let sampling_scale = match sampling_scale_mode {
        NumeratorSamplingScaleMode::None => "none",
        NumeratorSamplingScaleMode::BeyondQuadratic => "beyond_quadratic",
        NumeratorSamplingScaleMode::All => "all",
    }
    .to_string();

    let run = || -> Result<String> {
        let process_name = imported_process_name(dot_name)?;
        let graph = imported_graph_from_cli(cli, process_name, "default", 0)?;
        let parsed = graph.to_three_d_parsed_graph()?;
        if !bounds.is_empty() && parsed.external_names.is_empty() {
            // Keep this diagnostic explicit for vacuum topologies: the old Python helper used
            // energy monomials, while the automatic numerator helper requires an external vector.
            let _ = auto_numerator_expr_for_bounds(
                parsed.external_names.len(),
                &bounds,
                parsed.internal_edges.len(),
            );
        }

        let cff = generate_3d_expression_from_parsed(
            &parsed,
            &expression_options(RepresentationMode::Cff, &bounds, sampling_scale_mode),
        );
        let rhs_representation = match compare {
            ProbeComparison::CffPureLtd { .. } => RepresentationMode::PureLtd,
            _ => RepresentationMode::Ltd,
        };
        let rhs = generate_3d_expression_from_parsed(
            &parsed,
            &expression_options(rhs_representation, &bounds, sampling_scale_mode),
        );

        match compare {
            ProbeComparison::BuildOnly => {
                let cff = cff?;
                let rhs = rhs?;
                Ok(format!(
                    "ok: cff_orientations={}, rhs_orientations={}, cff_terms={}, rhs_terms={}",
                    cff.orientations.len(),
                    rhs.orientations.len(),
                    cff.num_unfolded_terms(),
                    rhs.num_unfolded_terms()
                ))
            }
            ProbeComparison::CffLtd { tolerance } | ProbeComparison::CffPureLtd { tolerance } => {
                let cff = cff?;
                let rhs = rhs?;
                let input = deterministic_input(&parsed, seed, uniform_scale)?;
                compare_expression_pair(&parsed, &cff, &rhs, &numerator, &input, tolerance)
            }
        }
    };

    match run() {
        Ok(detail) => ProbeCaseReport {
            name,
            dot_name: dot_name.to_string(),
            numerator,
            bounds,
            comparison,
            sampling_scale,
            status: "ok".to_string(),
            detail,
        },
        Err(error) => ProbeCaseReport {
            name,
            dot_name: dot_name.to_string(),
            numerator,
            bounds,
            comparison,
            sampling_scale,
            status: "failed".to_string(),
            detail: classify_probe_error(error),
        },
    }
}

fn classify_probe_error(error: color_eyre::Report) -> String {
    let text = error.to_string();
    if text.contains(&GenerationError::CffHigherEnergyPowerNotImplemented.to_string()) {
        format!("unsupported cff high-power sector: {text}")
    } else if text.contains("mismatch:") {
        format!("numerical mismatch: {text}")
    } else {
        text
    }
}

fn imported_graph_from_cli<'a>(
    cli: &'a gammaloop_integration_tests::CLIState,
    process_name: &str,
    integrand_name: &str,
    graph_id: usize,
) -> Result<&'a Graph> {
    let process = cli
        .state
        .process_list
        .processes
        .iter()
        .find(|process| process.definition.folder_name == process_name)
        .ok_or_else(|| eyre!("missing imported process '{process_name}'"))?;
    match &process.collection {
        ProcessCollection::Amplitudes(amplitudes) => amplitudes
            .get(integrand_name)
            .and_then(|amplitude| amplitude.graphs.get(graph_id))
            .map(|graph| &graph.graph)
            .ok_or_else(|| {
                eyre!("missing imported amplitude graph {graph_id} in integrand '{integrand_name}'")
            }),
        ProcessCollection::CrossSections(cross_sections) => cross_sections
            .get(integrand_name)
            .and_then(|cross_section| cross_section.supergraphs.get(graph_id))
            .map(|graph| &graph.graph)
            .ok_or_else(|| {
                eyre!(
                    "missing imported cross-section graph {graph_id} in integrand '{integrand_name}'"
                )
            }),
    }
}

fn source_internal_edges(graph: &Graph) -> Vec<usize> {
    graph
        .underlying
        .iter_edges()
        .filter_map(|(pair, edge_id, _)| pair.is_paired().then_some(edge_id.0))
        .collect()
}

fn load_imported_graph_case(
    manifest: &JsonValue,
    graph: &Graph,
    parsed: &ParsedGraph,
    name: &str,
) -> Result<ImportedGraphThreeDRepCase> {
    let expression = read_manifest_case_expression(manifest, name)?;
    let edge_map = graph
        .energy_edge_index_map(parsed)
        .ok_or_else(|| eyre!("imported graph did not expose an energy-edge index map"))?;
    let inverse_edge_map = EnergyEdgeIndexMap {
        internal: edge_map
            .internal
            .iter()
            .map(|(local_edge_id, source_edge_id)| (*source_edge_id, *local_edge_id))
            .collect(),
        external: edge_map
            .external
            .iter()
            .map(|(local_edge_id, source_edge_id)| (*source_edge_id, *local_edge_id))
            .collect(),
        orientation_edge_count: parsed.internal_edges.len(),
    };
    let compact_expression = expression
        .clone()
        .remap_energy_edge_indices(&inverse_edge_map);
    Ok(ImportedGraphThreeDRepCase {
        expression,
        compact_expression,
    })
}

fn manifest_case<'a>(manifest: &'a JsonValue, name: &str) -> Result<&'a JsonValue> {
    manifest["cases"]
        .as_array()
        .ok_or_else(|| eyre!("3Drep manifest has no cases array"))?
        .iter()
        .find(|case| case["name"].as_str() == Some(name))
        .ok_or_else(|| eyre!("3Drep manifest does not contain case '{name}'"))
}

fn assert_manifest_case(
    manifest: &JsonValue,
    name: &str,
    orientation_count: usize,
    surface_count: usize,
) -> Result<()> {
    let case = manifest_case(manifest, name)?;
    assert_eq!(
        case["orientation_count"].as_u64(),
        Some(orientation_count as u64),
        "unexpected orientation count for case {name}: {case:#}"
    );
    assert_eq!(
        case["evaluator_build_status"].as_str(),
        Some("ok"),
        "unexpected evaluator build status for case {name}: {case:#}"
    );
    let expression_path = case["expression_path"]
        .as_str()
        .ok_or_else(|| eyre!("3Drep case '{name}' has no expression_path"))?;
    let expression = serde_json::from_str::<ThreeDExpression<OrientationID>>(&fs::read_to_string(
        expression_path,
    )?)?;
    assert_eq!(
        expression.surfaces.linear_surface_cache.len(),
        surface_count,
        "unexpected surface count for case {name}"
    );
    Ok(())
}

fn read_manifest_case_expression(
    manifest: &JsonValue,
    name: &str,
) -> Result<ThreeDExpression<OrientationID>> {
    let case = manifest_case(manifest, name)?;
    let expression_path = case["expression_path"]
        .as_str()
        .ok_or_else(|| eyre!("3Drep case '{name}' has no expression_path"))?;
    Ok(serde_json::from_str::<ThreeDExpression<OrientationID>>(
        &fs::read_to_string(expression_path)?,
    )?)
}

#[test]
#[serial]
fn cli_imports_old_python_equivalent_threedrep_fixtures() -> Result<()> {
    let test_name = "threedreps_cli_imports_old_python_equivalent_fixtures";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;

    for spec in IMPORTED_GRAPH_SPECS {
        cli.run_command(&format!(
            "import graphs {} -p {} -i default -o",
            gammaloop_threedreps_dot_path(spec.dot_name).display(),
            spec.process_name
        ))?;
        let graph = imported_graph_from_cli(&cli, spec.process_name, "default", 0)?;
        let parsed = graph.to_three_d_parsed_graph()?;
        let repeated_sizes = three_dimensional_reps::repeated_groups(&parsed)
            .into_iter()
            .map(|group| group.edge_ids.len())
            .sorted()
            .collect::<Vec<_>>();
        assert_eq!(
            parsed.internal_edges.len(),
            spec.internal_edges,
            "{} internal edge count",
            spec.dot_name
        );
        assert_eq!(
            parsed.loop_names.len(),
            spec.loop_count,
            "{} loop count",
            spec.dot_name
        );
        assert_eq!(
            parsed.external_names.len(),
            spec.external_symbols,
            "{} external symbol count",
            spec.dot_name
        );
        assert_eq!(
            parsed.external_edges.len(),
            spec.external_edges,
            "{} external edge count",
            spec.dot_name
        );
        assert_eq!(
            repeated_sizes, spec.repeated_group_sizes,
            "{} repeated groups",
            spec.dot_name
        );
    }

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_reconstructs_energy_power_caps_from_imported_graph_numerators() -> Result<()> {
    let cases = [
        (
            "threedreps_power_caps_box_pow3_emr",
            "box_pow3.dot",
            vec![(0usize, 1usize), (1, 1), (3, 4)],
            None,
        ),
        (
            "threedreps_power_caps_sunrise_emr",
            "sunrise_pow4.dot",
            vec![(2usize, 5usize)],
            None,
        ),
        (
            "threedreps_power_caps_box_lmb",
            "box.dot",
            vec![(0usize, 3usize)],
            Some((0usize, 3usize)),
        ),
    ];

    for (test_name, dot_name, local_bounds, lmb_bound) in cases {
        let bootstrap = import_threedrep_graph(
            &format!("{test_name}_bootstrap"),
            &gammaloop_threedreps_dot_path(dot_name),
            "bootstrap",
        )?;
        let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
        let bootstrap_source_internal_edges = source_internal_edges(bootstrap_graph);
        let numerator = if let Some((loop_id, degree)) = lmb_bound {
            loop_energy_power_numerator(loop_id, degree)
        } else {
            edge_energy_power_numerator(&bootstrap_source_internal_edges, &local_bounds)
        };
        clean_test(&bootstrap.cli_settings.state.folder);

        let dot_path = dot_with_global_numerator(test_name, dot_name, &numerator)?;
        let process_name = format!("{test_name}_process");
        let state_path = get_tests_workspace_path().join(test_name).join("state");
        let workspace_path = imported_graph_threedrep_workspace(test_name);
        let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
        cli.run_command("import model scalars-default.json")?;
        cli.run_command(&format!(
            "import graphs {} -p {process_name} -i default -o",
            dot_path.display()
        ))?;
        cli.run_command(&format!(
            "3Drep test-cff-ltd -p {process_name} -i default -g 0 --workspace-path {}",
            workspace_path.display()
        ))?;
        let manifest = serde_json::from_str::<JsonValue>(&fs::read_to_string(
            workspace_path.join("test_cff_ltd_manifest.json"),
        )?)?;
        let graph = imported_graph_from_cli(&cli, &process_name, "default", 0)?;
        let source_internal_edges = source_internal_edges(graph);
        let local_bound_map = local_bounds.iter().copied().collect::<BTreeMap<_, _>>();
        let expected = source_internal_edges
            .iter()
            .enumerate()
            .map(|(local_edge_id, source_edge_id)| {
                (
                    *source_edge_id,
                    local_bound_map.get(&local_edge_id).copied().unwrap_or(0),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            serde_json::from_value::<Vec<(usize, usize)>>(
                manifest["automatic_energy_degree_bounds"].clone()
            )?,
            expected,
            "{test_name} automatic energy bounds for numerator {numerator}"
        );

        clean_test(get_tests_workspace_path().join(test_name));
        clean_test(&cli.cli_settings.state.folder);
    }

    Ok(())
}

#[test]
#[serial]
fn cli_graph_from_signatures_round_trips_pq_mass_expression_through_import_graphs() -> Result<()> {
    let test_name = "threedreps_graph_from_signatures_pq_mass";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let signatures_path = test_root.join("pq_signatures.txt");
    let dot_path = test_root.join("pq_from_signatures.dot");
    fs::create_dir_all(&test_root)?;

    let signature_expression = "prop(k1+p1,0)*prop(k1+p1-q1,0)*prop(k1-p2+q2,0)*prop(k1,0)";
    fs::write(&signatures_path, signature_expression)?;

    let expected = extract_signatures_and_masses_from_symbolica_expression(
        signature_expression,
        "k",
        &["p".to_string(), "q".to_string()],
        "prop(q_,m_)",
    )?;

    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command(&format!(
        "3Drep graph-from-signatures --signatures-file {} --dot-output {} --external-prefixes p,q --num-vertices 4",
        signatures_path.display(),
        dot_path.display()
    ))?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!(
        "import graphs {} -p threedreps_pq_closure -i default -o",
        dot_path.display()
    ))?;

    let graph = imported_graph_from_cli(&cli, "threedreps_pq_closure", "default", 0)?;
    let parsed = graph.to_three_d_parsed_graph()?;
    assert_eq!(
        sorted_signature_pairs(project_external_signatures_by_name(
            &parsed,
            &expected.external_names
        )),
        sorted_signature_pairs(expected.signatures)
    );
    assert!(
        parsed
            .internal_edges
            .iter()
            .all(|edge| edge.mass_key.as_deref() == Some("0"))
    );

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_old_python_case_matrix_records_current_three_drep_status() -> Result<()> {
    let test_name = "threedreps_old_python_case_matrix";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;

    for spec in IMPORTED_GRAPH_SPECS {
        cli.run_command(&format!(
            "import graphs {} -p {} -i default -o",
            gammaloop_threedreps_dot_path(spec.dot_name).display(),
            spec.process_name
        ))?;
    }

    let cases = old_python_probe_cases();
    let mut reports = Vec::new();
    for case in cases {
        reports.push(run_probe_case(&cli, case));
    }

    for spec in IMPORTED_GRAPH_SPECS
        .iter()
        .filter(|spec| REPEATED_EXAMPLES.contains(&spec.dot_name))
    {
        let graph = imported_graph_from_cli(&cli, spec.process_name, "default", 0)?;
        let parsed = graph.to_three_d_parsed_graph()?;
        for edge_id in 0..parsed.internal_edges.len() {
            let numerator = if parsed.external_names.is_empty() {
                format!("edges[{edge_id}][0]")
            } else {
                format!("dot(edges[{edge_id}], ext[0])")
            };
            reports.push(run_probe_case_parts(
                &cli,
                format!(
                    "old_python_all_edge_linear_numerator::{dot_name}::edge{edge_id}",
                    dot_name = spec.dot_name
                ),
                spec.dot_name,
                numerator,
                Vec::new(),
                1337,
                ProbeComparison::CffLtd { tolerance: 1.0e-8 },
                NumeratorSamplingScaleMode::None,
                None,
            ));
        }
    }

    for e0 in 0..=6 {
        for e1 in 0..=(6 - e0) {
            for e2 in 0..=(6 - e0 - e1) {
                for e3 in 0..=(6 - e0 - e1 - e2) {
                    let bounds = [e0, e1, e2, e3]
                        .into_iter()
                        .enumerate()
                        .filter_map(|(edge_id, power)| (power > 0).then_some((edge_id, power)))
                        .collect::<Vec<_>>();
                    reports.push(run_probe_case_parts(
                        &cli,
                        format!("old_python_normal_box_all_convergent_bounds::{e0}_{e1}_{e2}_{e3}"),
                        "box.dot",
                        local_edge_energy_monomial(&bounds),
                        bounds,
                        1337,
                        ProbeComparison::CffLtd { tolerance: 1.0e-8 },
                        NumeratorSamplingScaleMode::None,
                        None,
                    ));
                }
            }
        }
    }

    let report_path = test_root.join("old_python_case_matrix_report.json");
    let report_json = serde_json::to_string_pretty(&reports)?;
    fs::create_dir_all(&test_root)?;
    fs::write(&report_path, &report_json)?;
    fs::write(
        get_tests_workspace_path().join("threedreps_old_python_case_matrix_latest.json"),
        &report_json,
    )?;
    let ok_count = reports
        .iter()
        .filter(|report| report.status == "ok")
        .count();
    let failed_count = reports.len() - ok_count;
    println!(
        "old Python 3Drep case matrix: {ok_count} ok, {failed_count} failed, {} total",
        reports.len()
    );
    for report in reports.iter().filter(|report| report.status != "ok") {
        println!(
            "  failed: {} [{}] {}",
            report.name, report.comparison, report.detail
        );
    }

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[ignore = "slow parity probe matching the old HYBRID3D_RUN_SLOW five-loop ultimate-basis test"]
#[test]
#[serial]
fn cli_old_python_slow_five_loop_ultimate_basis_case_matrix() -> Result<()> {
    let test_name = "threedreps_old_python_slow_ultimate_basis_matrix";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;

    let dot_names = [
        "five_loop_ultimate_basis1.dot",
        "five_loop_ultimate_basis2.dot",
        "five_loop_ultimate_basis3.dot",
        "five_loop_ultimate_basis4.dot",
        "five_loop_ultimate_basis5.dot",
    ];
    let mut reports = Vec::new();
    for (index, dot_name) in dot_names.iter().enumerate() {
        let process_name = format!("threedreps_ultimate_basis_{index}");
        cli.run_command(&format!(
            "import graphs {} -p {process_name} -i default -o",
            gammaloop_threedreps_dot_path(dot_name).display(),
        ))?;
        let graph = imported_graph_from_cli(&cli, &process_name, "default", 0)?;
        let parsed = graph.to_three_d_parsed_graph()?;
        let source_edges = source_internal_edges(graph);
        let numerator = (0..parsed.internal_edges.len())
            .map(|edge_id| format!("dot(edges[{edge_id}], ext[{}])", edge_id % 4))
            .join(" + ");
        let local_case = ProbeCase {
            name: dot_name,
            dot_name,
            numerator: "generated all-edge external dot numerator",
            bounds: &[],
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        };
        let cff = generate_3d_expression_from_parsed(
            &parsed,
            &expression_options(
                RepresentationMode::Cff,
                &[],
                NumeratorSamplingScaleMode::None,
            ),
        );
        let ltd = generate_3d_expression_from_parsed(
            &parsed,
            &expression_options(
                RepresentationMode::Ltd,
                &[],
                NumeratorSamplingScaleMode::None,
            ),
        );
        let status = match (cff, ltd) {
            (Ok(cff), Ok(ltd)) => {
                let input = deterministic_input(&parsed, 1337, None)?;
                match compare_expression_pair(&parsed, &cff, &ltd, &numerator, &input, 1.0e-8) {
                    Ok(detail) => ("ok".to_string(), detail),
                    Err(error) => ("failed".to_string(), error.to_string()),
                }
            }
            (Err(error), _) | (_, Err(error)) => ("failed".to_string(), error.to_string()),
        };
        reports.push(ProbeCaseReport {
            name: local_case.name.to_string(),
            dot_name: local_case.dot_name.to_string(),
            numerator,
            bounds: source_edges
                .into_iter()
                .map(|edge_id| (edge_id, 0))
                .collect(),
            comparison: "cff_ltd".to_string(),
            sampling_scale: "none".to_string(),
            status: status.0,
            detail: status.1,
        });
    }

    fs::create_dir_all(&test_root)?;
    fs::write(
        test_root.join("slow_ultimate_basis_report.json"),
        serde_json::to_string_pretty(&reports)?,
    )?;
    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_imported_box_3drep_test_uses_gammaloop_graph_path() -> Result<()> {
    let test_name = "threedreps_cli_imported_box";
    let run = run_imported_graph_threedrep_test(test_name, "box.dot", "threedreps_box")?;
    let manifest = &run.manifest;

    assert_eq!(manifest["process_id"].as_u64(), Some(0));
    assert_eq!(manifest["integrand_name"].as_str(), Some("default"));
    assert_eq!(manifest["graph_id"].as_u64(), Some(0));
    assert_eq!(manifest["cases"].as_array().map(Vec::len), Some(3));
    assert!(
        manifest["override_energy_degree_bounds"]
            .as_array()
            .is_some_and(Vec::is_empty)
    );
    assert_eq!(
        serde_json::from_value::<Vec<(usize, usize)>>(manifest["energy_degree_bounds"].clone())?,
        vec![(4, 0), (5, 0), (6, 0), (7, 0)]
    );
    assert_internal_signatures(
        &run.parsed,
        &[
            ("q0", [0, 0, 0, 0]),
            ("q1", [1, 0, 0, 0]),
            ("q2", [1, 1, 0, 0]),
            ("q3", [1, 1, 1, 0]),
        ],
    );

    assert_manifest_case(manifest, "cff_none", 14, 12)?;
    assert_manifest_case(manifest, "ltd_none", 4, 24)?;
    assert_manifest_case(manifest, "pureltd_none", 4, 24)?;

    assert_eq!(
        compact_orientation_labels(&run.cff.expression, &run.source_internal_edges),
        nontrivial_sign_labels(4)
    );
    assert_unit_or_minus_one_prefactors(&run.cff.expression, (1, 1))?;
    assert!(
        half_edges(&run.cff.expression)
            .iter()
            .all(|edges| edges == &[4, 5, 6, 7])
    );
    assert_eq!(
        linear_denominator_surface_ids(&run.cff.expression),
        vec![
            vec![0, 1, 2],
            vec![3, 4, 5],
            vec![0, 1, 4, 5],
            vec![6, 7, 8],
            vec![0, 2, 6, 8],
            vec![3, 4, 7, 8],
            vec![0, 4, 8],
            vec![9, 10, 11],
            vec![1, 2, 10, 11],
            vec![3, 5, 9, 11],
            vec![1, 5, 11],
            vec![6, 7, 9, 10],
            vec![2, 6, 10],
            vec![3, 7, 9],
        ]
    );

    assert_eq!(
        compact_orientation_labels(&run.ltd.expression, &run.source_internal_edges),
        ["+xxx", "x+xx", "xx+x", "xxx+"]
    );
    assert_unit_or_minus_one_prefactors(&run.ltd.expression, (-1, 1))?;
    assert_eq!(
        half_edges(&run.ltd.expression),
        vec![vec![4], vec![5], vec![6], vec![7]]
    );
    assert_eq!(
        linear_denominator_surface_ids(&run.ltd.expression),
        vec![
            vec![1, 3, 5],
            vec![7, 9, 11],
            vec![13, 15, 17],
            vec![19, 21, 23]
        ]
    );
    assert_eq!(
        half_edges(&run.pure_ltd.expression),
        vec![vec![4], vec![5], vec![6], vec![7]]
    );

    let input = python_reference_input(vec![1.138394, 0.758865, 0.458809, 1.029287]);
    let cff_value =
        evaluate_expression(&run.parsed, &run.cff.compact_expression, "1", &input)?.value;
    let ltd_value =
        evaluate_expression(&run.parsed, &run.ltd.compact_expression, "1", &input)?.value;
    let python_reference = 0.513_288_156_946_660_9_f64;
    assert!(
        (cff_value - python_reference).abs() < 1.0e-13,
        "cff_value={cff_value:.17e}, python_reference={python_reference:.17e}"
    );
    assert!(
        (ltd_value - python_reference).abs() < 1.0e-13,
        "ltd_value={ltd_value:.17e}, python_reference={python_reference:.17e}"
    );
    assert!(
        (cff_value - ltd_value).abs() < 1.0e-13,
        "cff_value={cff_value:.17e}, ltd_value={ltd_value:.17e}"
    );

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&run.cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_imported_box_pow3_3drep_test_uses_gammaloop_graph_path() -> Result<()> {
    let test_name = "threedreps_cli_imported_box_pow3";
    let run = run_imported_graph_threedrep_test(test_name, "box_pow3.dot", "threedreps_box_pow3")?;
    let manifest = &run.manifest;

    assert_eq!(manifest["process_id"].as_u64(), Some(0));
    assert_eq!(manifest["integrand_name"].as_str(), Some("default"));
    assert_eq!(manifest["graph_id"].as_u64(), Some(0));
    assert_eq!(manifest["cases"].as_array().map(Vec::len), Some(3));
    assert!(
        manifest["override_energy_degree_bounds"]
            .as_array()
            .is_some_and(Vec::is_empty)
    );
    assert_eq!(
        serde_json::from_value::<Vec<(usize, usize)>>(manifest["energy_degree_bounds"].clone())?,
        vec![(4, 0), (5, 0), (6, 0), (7, 0), (8, 0), (9, 0)]
    );
    assert_internal_signatures(
        &run.parsed,
        &[
            ("q0", [0, 0, 0, 0]),
            ("q1", [1, 0, 0, 0]),
            ("q2", [1, 1, 0, 0]),
            ("q3", [1, 1, 1, 0]),
            ("q4", [1, 1, 1, 0]),
            ("q5", [1, 1, 1, 0]),
        ],
    );

    assert_manifest_case(manifest, "cff_none", 62, 27)?;
    assert_manifest_case(manifest, "ltd_none", 4, 24)?;
    assert_manifest_case(manifest, "pureltd_none", 6, 57)?;

    assert_eq!(
        compact_orientation_labels(&run.cff.expression, &run.source_internal_edges),
        nontrivial_sign_labels(6)
    );
    assert_unit_or_minus_one_prefactors(&run.cff.expression, (1, 1))?;
    assert!(
        half_edges(&run.cff.expression)
            .iter()
            .all(|edges| edges == &[4, 5, 6, 7, 8, 9])
    );
    assert_eq!(
        compact_base_orientation_labels(&run.ltd.expression, &run.source_internal_edges),
        ["+xxxxx", "x+xxxx", "xx+xxx", "xxx+xx"]
    );
    assert_eq!(
        variant_prefactors(&run.ltd.expression),
        vec![
            vec![(-1, 1)],
            vec![(-1, 1)],
            vec![(-1, 1)],
            vec![(-1, 1), (-3, 1), (-6, 1)]
        ]
    );
    assert_eq!(
        variant_half_edges(&run.ltd.expression),
        vec![
            vec![vec![4]],
            vec![vec![5]],
            vec![vec![6]],
            vec![vec![7, 7, 7], vec![7, 7, 7, 7], vec![7, 7, 7, 7, 7]]
        ]
    );
    assert_eq!(
        variant_linear_denominator_surface_ids(&run.ltd.expression),
        vec![
            vec![vec![1, 3, 5]],
            vec![vec![7, 9, 11]],
            vec![vec![13, 15, 17]],
            vec![vec![19, 21, 23], vec![19, 21, 23], vec![19, 21, 23]],
        ]
    );

    assert_eq!(
        compact_base_orientation_labels(&run.pure_ltd.expression, &run.source_internal_edges),
        ["+xxxxx", "x+xxxx", "xx+xxx", "xxx+xx", "xxxx+x", "xxxxx+"]
    );
    assert_unit_or_minus_one_prefactors(&run.pure_ltd.expression, (-1, 1))?;
    assert_eq!(
        half_edges(&run.pure_ltd.expression),
        vec![vec![4], vec![5], vec![6], vec![7], vec![8], vec![9]]
    );

    let repeated_input = python_reference_input(vec![
        1.138394, 0.758865, 0.458809, 1.029287, 1.029287, 1.029287,
    ]);
    let cff_value = evaluate_expression(
        &run.parsed,
        &run.cff.compact_expression,
        "1",
        &repeated_input,
    )?
    .value;
    let ltd_value = evaluate_expression(
        &run.parsed,
        &run.ltd.compact_expression,
        "1",
        &repeated_input,
    )?
    .value;
    let python_reference = 0.338_108_873_066_020_44_f64;
    let repeated_numeric_tolerance = 1.0e-10;
    assert!(
        (cff_value - python_reference).abs() < 1.0e-13,
        "cff_value={cff_value:.17e}, python_reference={python_reference:.17e}"
    );
    assert!(
        (ltd_value - python_reference).abs() < repeated_numeric_tolerance,
        "ltd_value={ltd_value:.17e}, python_reference={python_reference:.17e}"
    );
    assert!(
        (cff_value - ltd_value).abs() < repeated_numeric_tolerance,
        "cff_value={cff_value:.17e}, ltd_value={ltd_value:.17e}"
    );

    let split_references = [
        (0.1, 0.342_116_187_424_781_67),
        (0.05, 0.339_106_080_530_878_72),
        (0.025, 0.338_357_886_377_202_07),
        (0.0125, 0.338_171_108_363_500_58),
    ];
    let mut previous_distance = f64::INFINITY;
    for (epsilon, expected) in split_references {
        let input = python_reference_input(repeated_box_masses_with_split(epsilon));
        let value =
            evaluate_expression(&run.parsed, &run.pure_ltd.compact_expression, "1", &input)?.value;
        let distance = (value - python_reference).abs();
        assert!(
            (value - expected).abs() < repeated_numeric_tolerance,
            "epsilon={epsilon}: value={value:.17e}, expected={expected:.17e}"
        );
        assert!(
            distance < previous_distance,
            "epsilon={epsilon}: pure-LTD split distance {distance:.17e} did not decrease from {previous_distance:.17e}"
        );
        previous_distance = distance;
    }

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&run.cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_graph_from_signatures_mercedes_closes_through_import_graphs() -> Result<()> {
    let test_name = "threedreps_graph_from_signatures_mercedes";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let signatures_path = test_root.join("mercedes_signatures.txt");
    let dot_path = test_root.join("mercedes_from_signatures.dot");
    fs::create_dir_all(&test_root)?;

    let signature_expression = [
        "prop(k1,0)",
        "prop(k2,0)",
        "prop(k3,0)",
        "prop(k1+k2+p1,0)",
        "prop(k2+k3+p2,0)",
        "prop(k1-k3+p3-p4,0)",
    ]
    .join("*");
    fs::write(&signatures_path, &signature_expression)?;

    let expected = extract_signatures_and_masses_from_symbolica_expression(
        &signature_expression,
        "k",
        &["p".to_string()],
        "prop(q_,m_)",
    )?;

    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command(&format!(
        "3Drep graph-from-signatures --signatures-file {} --dot-output {} --num-vertices 4",
        signatures_path.display(),
        dot_path.display()
    ))?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!(
        "import graphs {} -p threedreps_mercedes_closure -i default -o",
        dot_path.display()
    ))?;

    let graph = imported_graph_from_cli(&cli, "threedreps_mercedes_closure", "default", 0)?;
    let parsed = graph.to_three_d_parsed_graph()?;
    assert_eq!(parsed.internal_edges.len(), expected.signatures.len());
    assert_eq!(parsed.loop_names.len(), expected.loop_names.len());
    assert!(parsed.external_names.len() >= expected.external_names.len());
    assert_eq!(
        sorted_signature_pairs(project_external_signatures_by_name(
            &parsed,
            &expected.external_names
        )),
        sorted_signature_pairs(expected.signatures)
    );

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
