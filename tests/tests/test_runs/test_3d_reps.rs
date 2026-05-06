use std::collections::{BTreeMap, BTreeSet};
use std::fs;

use color_eyre::eyre::eyre;
use gammalooprs::{graph::Graph, processes::ProcessCollection};
use symbolica::atom::AtomCore;
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

fn shell_single_quote(value: &str) -> String {
    format!("'{}'", value.replace('\'', "'\\''"))
}

fn find_named_artifact(root: &Path, file_name: &str) -> Result<PathBuf> {
    let mut stack = vec![root.to_path_buf()];
    let mut matches = Vec::new();
    while let Some(path) = stack.pop() {
        for entry in fs::read_dir(&path)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.file_name().and_then(|name| name.to_str()) == Some(file_name) {
                matches.push(path);
            }
        }
    }
    matches.sort();
    match matches.as_slice() {
        [path] => Ok(path.clone()),
        [] => Err(eyre!(
            "Could not find artifact '{file_name}' under {}",
            root.display()
        )),
        many => Err(eyre!(
            "Artifact '{file_name}' under {} is ambiguous: {:?}",
            root.display(),
            many
        )),
    }
}

fn latest_oriented_expression_path(workspace: &Path) -> Result<PathBuf> {
    let pointer = workspace.join("latest_oriented_expression_path.txt");
    let path = PathBuf::from(fs::read_to_string(&pointer)?.trim());
    if path.is_absolute() {
        Ok(path)
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}

fn run_on_large_stack(name: &str, f: fn() -> Result<()>) -> Result<()> {
    std::thread::Builder::new()
        .name(name.to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(f)?
        .join()
        .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
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
        repeated_group_sizes: &[2, 3],
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
    CffPureLtdBuildOnly,
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
const DENSE_BOUNDS_NUMERATOR: &str = "auto_dense_bounds";
const DENSE_LOOP0_BOUNDS_NUMERATOR: &str = "auto_dense_loop0_bounds";
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

fn project_external_signature_by_name(
    parsed: &ParsedGraph,
    signature: &MomentumSignature,
    external_names: &[String],
) -> Result<MomentumSignature> {
    let mut external_signature = vec![0; external_names.len()];
    for (source_external_id, coeff) in signature.external_signature.iter().enumerate() {
        if *coeff == 0 {
            continue;
        }
        let name = &parsed.external_names[source_external_id];
        let target_external_id = external_names
            .iter()
            .position(|candidate| candidate == name)
            .ok_or_else(|| eyre!("unexpected external name {name}"))?;
        external_signature[target_external_id] += coeff;
    }
    Ok(MomentumSignature {
        loop_signature: signature.loop_signature.clone(),
        external_signature,
    })
}

fn signature_mass_multiset(
    signatures: impl IntoIterator<Item = MomentumSignature>,
    masses: impl IntoIterator<Item = String>,
) -> BTreeMap<String, usize> {
    let mut multiset = BTreeMap::new();
    for (signature, mass) in signatures.into_iter().zip(masses) {
        let (canonical, _) = signature.canonical_up_to_sign();
        let key = format!(
            "{:?}|{:?}|{}",
            canonical.loop_signature,
            canonical.external_signature,
            mass.trim()
        );
        *multiset.entry(key).or_default() += 1;
    }
    multiset
}

fn signature_expression_from_factors(factors: &[(&str, &str)]) -> String {
    factors
        .iter()
        .map(|(momentum, mass)| format!("prop({momentum},{mass})"))
        .join("*")
}

fn assert_graph_from_signatures_cli_closure(
    test_name: &str,
    factors: &[(&str, &str)],
    external_prefixes: &[&str],
    num_vertices: usize,
) -> Result<()> {
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let dot_path = test_root.join("from_signatures.dot");
    fs::create_dir_all(&test_root)?;

    let signature_expression = signature_expression_from_factors(factors);
    let expected = extract_signatures_and_masses_from_symbolica_expression(
        &signature_expression,
        "k",
        &external_prefixes
            .iter()
            .map(|prefix| prefix.to_string())
            .collect::<Vec<_>>(),
        "prop(q_,m_)",
    )?;
    let external_prefixes_arg = external_prefixes.join(",");
    let process_name = test_name.replace('-', "_");

    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command(&format!(
        "3Drep graph-from-signatures --signatures {} --dot-output {} --external-prefixes {} --num-vertices {}",
        shell_single_quote(&signature_expression),
        dot_path.display(),
        external_prefixes_arg,
        num_vertices,
    ))?;
    let dot = fs::read_to_string(&dot_path)?;
    assert!(
        !dot.contains("lmb_rep"),
        "graph-from-signatures output must not carry explicit lmb_rep attributes:\n{dot}"
    );

    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!(
        "import graphs --inline-dot {} -p {process_name} -i default -o",
        shell_single_quote(&dot)
    ))?;

    let graph = imported_graph_from_cli(&cli, &process_name, "default", 0)?;
    let parsed = graph.to_three_d_parsed_graph()?;
    assert_eq!(parsed.internal_edges.len(), expected.signatures.len());
    assert_eq!(parsed.loop_names.len(), expected.loop_names.len());

    let actual_signatures = parsed
        .internal_edges
        .iter()
        .map(|edge| {
            project_external_signature_by_name(&parsed, &edge.signature, &expected.external_names)
        })
        .collect::<Result<Vec<_>>>()?;
    let actual_masses = parsed
        .internal_edges
        .iter()
        .map(|edge| edge.mass_key.clone().unwrap_or_else(|| "0".to_string()))
        .collect::<Vec<_>>();
    assert_eq!(
        signature_mass_multiset(actual_signatures, actual_masses),
        signature_mass_multiset(expected.signatures, expected.masses),
        "closure mismatch for input `{signature_expression}` and generated DOT:\n{dot}"
    );

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

type SignatureClosureCase<'a> = (&'a str, &'a [(&'a str, &'a str)], &'a [&'a str], usize);

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

fn expected_prefactor_label(expected: (i64, i64)) -> String {
    if expected.1 == 1 {
        expected.0.to_string()
    } else {
        format!("{}/{}", expected.0, expected.1)
    }
}

fn variant_prefactors(expression: &ThreeDExpression<OrientationID>) -> Vec<Vec<String>> {
    expression
        .orientations
        .iter()
        .map(|orientation| {
            orientation
                .variants
                .iter()
                .map(|variant| variant.prefactor.to_canonical_string())
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
        let expected_label = expected_prefactor_label(expected);
        if orientation.variants[0].prefactor.to_canonical_string() != expected_label {
            return Err(eyre!(
                "orientation {} has prefactor {}; expected {}",
                orientation.data.label.as_deref().unwrap_or("<unlabeled>"),
                orientation.variants[0].prefactor.to_canonical_string(),
                expected_label
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
        include_cff_duplicate_signature_excess_sign: true,
        preserve_internal_edges_as_four_d_denominators: Vec::new(),
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

fn dense_edge_energy_power_numerator(
    source_internal_edges: &[usize],
    source_external_edges: &[usize],
    bounds: &[(usize, usize)],
) -> String {
    assert!(
        !source_external_edges.is_empty(),
        "dense edge-energy numerator construction needs at least one external momentum"
    );
    let mut factors = Vec::new();
    for (local_edge_id, degree) in bounds {
        let edge_id = source_internal_edges[*local_edge_id];
        for power_index in 0..*degree {
            let external_id =
                source_external_edges[(local_edge_id + power_index) % source_external_edges.len()];
            factors.push(imported_q_dot_q(edge_id, external_id));
        }
    }
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

fn dense_loop_energy_power_numerator(
    loop_id: usize,
    degree: usize,
    source_external_edges: &[usize],
) -> String {
    assert!(
        !source_external_edges.is_empty(),
        "dense loop-energy numerator construction needs at least one external momentum"
    );
    let mut factors = Vec::new();
    for power_index in 0..degree {
        let external_id =
            source_external_edges[(loop_id + power_index) % source_external_edges.len()];
        factors.push(imported_k_dot_q(loop_id, external_id));
    }
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

fn imported_q_dot_q(lhs_edge_id: usize, rhs_edge_id: usize) -> String {
    format!(
        "(Q({lhs_edge_id},spenso::cind(0))*Q({rhs_edge_id},spenso::cind(0))-Q({lhs_edge_id},spenso::cind(1))*Q({rhs_edge_id},spenso::cind(1))-Q({lhs_edge_id},spenso::cind(2))*Q({rhs_edge_id},spenso::cind(2))-Q({lhs_edge_id},spenso::cind(3))*Q({rhs_edge_id},spenso::cind(3)))"
    )
}

fn imported_k_dot_q(loop_id: usize, edge_id: usize) -> String {
    format!(
        "(K({loop_id},spenso::cind(0))*Q({edge_id},spenso::cind(0))-K({loop_id},spenso::cind(1))*Q({edge_id},spenso::cind(1))-K({loop_id},spenso::cind(2))*Q({edge_id},spenso::cind(2))-K({loop_id},spenso::cind(3))*Q({edge_id},spenso::cind(3)))"
    )
}

fn q0_component(edge_id: usize) -> String {
    format!("Q({edge_id},spenso::cind(0))")
}

fn energy_component_power(edge_id: usize, degree: usize) -> String {
    if degree == 0 {
        return "1".to_string();
    }
    std::iter::repeat_with(|| q0_component(edge_id))
        .take(degree)
        .join("*")
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

#[derive(Debug, Clone)]
struct ManifestEvaluationValue {
    id: usize,
    re: f64,
    im: f64,
    backend: String,
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
        "3Drep test-cff-ltd -p {process_name} -i default -g 0 --workspace-path {} --clean",
        workspace_path.display()
    ))?;

    let manifest_path = find_named_artifact(&workspace_path, "test_cff_ltd_manifest.json")?;
    let manifest = serde_json::from_str::<JsonValue>(&fs::read_to_string(&manifest_path)?)?;
    assert_threedrep_comparison_success(&manifest, test_name)?;
    assert_threedrep_evaluation_ids(&manifest, test_name)?;
    let cases = manifest["cases"]
        .as_array()
        .ok_or_else(|| eyre!("3Drep comparison manifest has no cases array"))?;
    for case in cases {
        if case["evaluator_build_status"]
            .as_str()
            .is_some_and(|status| status.starts_with("ok"))
        {
            let evaluations = case["evaluations"]
                .as_array()
                .ok_or_else(|| eyre!("3Drep comparison case has no evaluations array"))?;
            assert!(
                !evaluations.is_empty(),
                "3Drep comparison case should record at least one evaluation"
            );
            for evaluation in evaluations {
                assert!(
                    evaluation["value"].as_str().is_some(),
                    "3Drep comparison evaluation should record a value"
                );
                assert!(
                    !evaluation["value"]
                        .as_str()
                        .is_some_and(|value| value.contains("NaN")),
                    "3Drep comparison evaluation should avoid singular diagnostic points"
                );
                assert!(
                    evaluation["sample_evaluation_timing_seconds"]
                        .as_f64()
                        .is_some(),
                    "3Drep comparison evaluation should record wall timing"
                );
                assert!(
                    evaluation["evaluator_build_timing_seconds"]
                        .as_f64()
                        .is_some(),
                    "3Drep comparison evaluation should record evaluator build timing"
                );
            }
            if case["name"].as_str() == Some("pureltd") && dot_name.contains("pow") {
                assert!(
                    evaluations
                        .iter()
                        .any(|evaluation| evaluation["mass_shift_values"]
                            .as_array()
                            .is_some_and(|values| !values.is_empty())),
                    "pure-LTD repeated-propagator diagnostics should report the split masses"
                );
            }
        }
    }
    let graph = imported_graph_from_cli(&cli, process_name, "default", 0)?;
    let parsed = graph.to_three_d_parsed_graph()?;
    let source_internal_edges = source_internal_edges(graph);
    let cff = load_imported_graph_case(&manifest, graph, &parsed, "cff")?;
    let ltd = load_imported_graph_case(&manifest, graph, &parsed, "ltd")?;
    let pure_ltd = load_imported_graph_case(&manifest, graph, &parsed, "pureltd")?;

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

fn assert_threedrep_comparison_success(manifest: &JsonValue, context: &str) -> Result<()> {
    let status = manifest["verdict"]["status"]
        .as_str()
        .ok_or_else(|| eyre!("3Drep comparison manifest has no verdict status"))?;
    if status == "Success" {
        return Ok(());
    }

    let checks = manifest["verdict"]["checks"]
        .as_array()
        .map(|checks| {
            checks
                .iter()
                .filter(|check| check["status"].as_str() != Some("Success"))
                .map(|check| {
                    format!(
                        "{}: {} (abs {}, rel {}, tol {})",
                        check["name"].as_str().unwrap_or("<unnamed>"),
                        check["message"].as_str().unwrap_or("<no message>"),
                        check["abs_diff"].as_str().unwrap_or("-"),
                        check["rel_diff"].as_str().unwrap_or("-"),
                        check["tolerance"].as_str().unwrap_or("-"),
                    )
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();

    panic!(
        "3Drep comparison for {context} reported {status}; failing checks:\n{}",
        checks.join("\n")
    );
}

fn assert_threedrep_evaluation_ids(manifest: &JsonValue, context: &str) -> Result<()> {
    let cases = manifest["cases"]
        .as_array()
        .ok_or_else(|| eyre!("3Drep comparison manifest has no cases array"))?;
    let mut expected_id = 0u64;
    for evaluation in cases
        .iter()
        .flat_map(|case| case["evaluations"].as_array().into_iter().flatten())
    {
        assert_eq!(
            evaluation["id"].as_u64(),
            Some(expected_id),
            "3Drep comparison for {context} should assign stable sequential evaluation IDs"
        );
        expected_id += 1;
    }
    assert!(
        expected_id > 0,
        "3Drep comparison for {context} should contain evaluation IDs"
    );
    Ok(())
}

fn set_threedrep_sampling_mode(
    cli: &mut gammaloop_integration_tests::CLIState,
    mode: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set global kv global.generation.uniform_numerator_sampling_scale={mode}"
    ))
}

fn set_threedrep_compile_backend(
    cli: &mut gammaloop_integration_tests::CLIState,
    backend: Option<&str>,
) -> Result<()> {
    match backend {
        None => cli.run_command("set global kv global.generation.evaluator.compile=false"),
        Some(backend) => {
            cli.run_command(&format!(
                "set global string '\n[global.generation.evaluator]\ncompile = true\n\n[global.generation.compile]\ncompilation_mode = \"{backend}\"\noptimization_level = \"O0\"\nfast_math = false\nunsafe_math = false\ncustom = []\n'"
            ))
        }
    }
}

fn run_threedrep_test_cff_ltd_manifest(
    cli: &mut gammaloop_integration_tests::CLIState,
    workspace_path: &Path,
    process_name: &str,
    precision: &str,
    seed: u64,
    scale: f64,
) -> Result<JsonValue> {
    run_threedrep_test_cff_ltd_manifest_for_graph(
        cli,
        workspace_path,
        process_name,
        0,
        precision,
        seed,
        scale,
        None,
    )
}

fn run_threedrep_test_cff_ltd_manifest_for_graph(
    cli: &mut gammaloop_integration_tests::CLIState,
    workspace_path: &Path,
    process_name: &str,
    graph_id: usize,
    precision: &str,
    seed: u64,
    scale: f64,
    energy_degree_bounds: Option<&str>,
) -> Result<JsonValue> {
    let energy_degree_bounds = energy_degree_bounds
        .map(|bounds| format!(" --energy-degree-bounds {bounds}"))
        .unwrap_or_default();
    cli.run_command(&format!(
        "3Drep test-cff-ltd -p {process_name} -i default -g {graph_id} --workspace-path {} --precision {precision} --seed {seed} --scale {scale:.17e}{energy_degree_bounds} --clean",
        workspace_path.display()
    ))?;
    let manifest_path = find_named_artifact(workspace_path, "test_cff_ltd_manifest.json")?;
    let manifest = serde_json::from_str::<JsonValue>(&fs::read_to_string(manifest_path)?)?;
    assert_threedrep_comparison_success(&manifest, process_name)?;
    assert_threedrep_evaluation_ids(&manifest, process_name)?;
    Ok(manifest)
}

fn direct_manifest_evaluation_value(
    manifest: &JsonValue,
    case_name: &str,
) -> Result<ManifestEvaluationValue> {
    let case = manifest_case(manifest, case_name)?;
    let evaluation = case["evaluations"]
        .as_array()
        .ok_or_else(|| eyre!("3Drep case '{case_name}' has no evaluations array"))?
        .iter()
        .find(|evaluation| {
            evaluation["status"].as_str() == Some("ok")
                && evaluation["mass_shift_values"]
                    .as_array()
                    .is_some_and(Vec::is_empty)
        })
        .ok_or_else(|| eyre!("3Drep case '{case_name}' has no direct successful evaluation"))?;
    let re = evaluation["value_re"]
        .as_str()
        .ok_or_else(|| eyre!("3Drep case '{case_name}' evaluation has no real part"))?
        .parse::<f64>()?;
    let im = evaluation["value_im"]
        .as_str()
        .ok_or_else(|| eyre!("3Drep case '{case_name}' evaluation has no imaginary part"))?
        .parse::<f64>()?;
    assert!(
        re.is_finite() && im.is_finite(),
        "3Drep case {case_name} direct value should be finite: ({re}, {im})"
    );
    Ok(ManifestEvaluationValue {
        id: evaluation["id"]
            .as_u64()
            .ok_or_else(|| eyre!("3Drep case '{case_name}' evaluation has no id"))?
            as usize,
        re,
        im,
        backend: evaluation["evaluator_backend"]
            .as_str()
            .ok_or_else(|| eyre!("3Drep case '{case_name}' evaluation has no backend"))?
            .to_string(),
    })
}

fn assert_manifest_direct_values_agree(
    lhs: &ManifestEvaluationValue,
    rhs: &ManifestEvaluationValue,
    context: &str,
    tolerance: f64,
) {
    let abs_diff = (lhs.re - rhs.re).hypot(lhs.im - rhs.im);
    let rhs_abs = rhs.re.hypot(rhs.im);
    let rel_diff = abs_diff
        / if rhs_abs > 0.0 {
            rhs_abs
        } else {
            f64::MIN_POSITIVE
        };
    assert!(
        abs_diff <= 1.0e-12 || rel_diff <= tolerance,
        "{context}: values differ beyond tolerance: lhs #{}, rhs #{}, lhs=({:.17e},{:.17e}), rhs=({:.17e},{:.17e}), abs_diff={:.3e}, rel_diff={:.3e}, tolerance={:.3e}",
        lhs.id,
        rhs.id,
        lhs.re,
        lhs.im,
        rhs.re,
        rhs.im,
        abs_diff,
        rel_diff,
        tolerance
    );
}

fn assert_manifest_uses_backend(
    manifest: &JsonValue,
    expected_backend: &str,
    context: &str,
) -> Result<()> {
    let actual = manifest["settings"]["evaluator_backend"]
        .as_str()
        .ok_or_else(|| eyre!("3Drep manifest for {context} has no evaluator backend setting"))?;
    assert!(
        actual.contains(expected_backend),
        "3Drep manifest for {context} used backend '{actual}', expected it to contain '{expected_backend}'"
    );
    for case in manifest["cases"]
        .as_array()
        .ok_or_else(|| eyre!("3Drep manifest for {context} has no cases array"))?
    {
        for evaluation in case["evaluations"]
            .as_array()
            .ok_or_else(|| eyre!("3Drep manifest for {context} case has no evaluations array"))?
        {
            let backend = evaluation["evaluator_backend"].as_str().unwrap_or("-");
            assert!(
                backend.contains(expected_backend),
                "3Drep evaluation in {context} used backend '{backend}', expected it to contain '{expected_backend}'"
            );
        }
    }
    Ok(())
}

fn read_evaluate_manifest(workspace_path: &Path) -> Result<JsonValue> {
    let manifest_path = find_named_artifact(workspace_path, "evaluate_manifest.json")?;
    Ok(serde_json::from_str::<JsonValue>(&fs::read_to_string(
        manifest_path,
    )?)?)
}

fn assert_evaluate_manifest_ok(manifest: &JsonValue, context: &str) -> Result<()> {
    assert_eq!(
        manifest["evaluation"]["status"].as_str(),
        Some("ok"),
        "3Drep evaluate for {context} should succeed: {manifest:#}"
    );
    assert!(
        manifest["evaluation"]["value"].as_str().is_some(),
        "3Drep evaluate for {context} should record the evaluated value"
    );
    assert!(
        !manifest["evaluation"]["value"]
            .as_str()
            .is_some_and(|value| value.contains("NaN")),
        "3Drep evaluate for {context} should avoid singular diagnostic points"
    );
    assert!(
        manifest["evaluation"]["sample_evaluation_timing_seconds"]
            .as_f64()
            .is_some(),
        "3Drep evaluate for {context} should record sample timing"
    );
    for field in [
        "symbolica_expression_path",
        "symbolica_expression_pretty_path",
        "symbolica_expression_raw_path",
        "symbolica_expression_raw_script_path",
        "param_builder_path",
    ] {
        let path = manifest[field].as_str().ok_or_else(|| {
            eyre!("3Drep evaluate manifest for {context} has no {field}: {manifest:#}")
        })?;
        assert!(
            PathBuf::from(path).exists(),
            "3Drep evaluate manifest for {context} points {field} to a missing file: {path}"
        );
    }
    let raw_expression_path = manifest["symbolica_expression_raw_path"]
        .as_str()
        .ok_or_else(|| {
            eyre!("3Drep evaluate manifest for {context} has no symbolica_expression_raw_path")
        })?;
    assert!(
        !fs::read_to_string(raw_expression_path)?.trim().is_empty(),
        "3Drep evaluate for {context} should write a non-empty raw Symbolica evaluator input"
    );
    let raw_expression =
        serde_json::from_str::<JsonValue>(&fs::read_to_string(raw_expression_path)?)?;
    assert_eq!(
        raw_expression["schema_version"].as_u64(),
        Some(2),
        "3Drep evaluate for {context} should write a versioned raw Symbolica evaluator input"
    );
    assert!(
        raw_expression["evaluator_backend"].as_str().is_some(),
        "3Drep evaluate for {context} should record the archived GammaLoop evaluator backend"
    );
    assert!(
        raw_expression["parameters"]
            .as_array()
            .is_some_and(|parameters| !parameters.is_empty()),
        "3Drep evaluate for {context} should record raw Symbolica evaluator parameters"
    );
    assert!(
        raw_expression["calls"]
            .as_array()
            .is_some_and(|calls| !calls.is_empty()),
        "3Drep evaluate for {context} should record raw Symbolica evaluator calls"
    );
    let raw_script_path = manifest["symbolica_expression_raw_script_path"]
        .as_str()
        .ok_or_else(|| {
            eyre!(
                "3Drep evaluate manifest for {context} has no symbolica_expression_raw_script_path"
            )
        })?;
    let raw_script = fs::read_to_string(raw_script_path)?;
    assert!(
        raw_script.starts_with("#!/usr/bin/env -S rust-script --debug"),
        "3Drep evaluate for {context} should write an executable rust-script replay"
    );
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        assert!(
            fs::metadata(raw_script_path)?.permissions().mode() & 0o111 != 0,
            "3Drep evaluate for {context} should make the rust-script replay executable"
        );
    }
    Ok(())
}

fn assert_evaluate_reused_cached_evaluator(manifest: &JsonValue, context: &str) -> Result<()> {
    assert_evaluate_manifest_ok(manifest, context)?;
    assert!(
        manifest["evaluation"]["evaluator_build_timing_seconds"].is_null(),
        "3Drep evaluate for {context} should have reused the evaluator cache: {manifest:#}"
    );
    Ok(())
}

fn assert_numerator_only_q0_override(
    manifest: &JsonValue,
    edge_id: usize,
    expected_value: f64,
    context: &str,
) -> Result<()> {
    assert_eq!(
        manifest["numerator_only"].as_bool(),
        Some(true),
        "3Drep evaluate for {context} should be in numerator-only mode"
    );
    assert!(
        manifest["expression_path"].is_null(),
        "direct numerator-only evaluation for {context} should not require an oriented JSON path"
    );
    let parameter = manifest["parameters"]
        .as_array()
        .ok_or_else(|| eyre!("3Drep evaluate manifest for {context} has no parameter list"))?
        .iter()
        .find(|parameter| {
            parameter["canonical_name"]
                .as_str()
                .is_some_and(|name| name.contains(&format!("Q({edge_id},")) && name.contains("cind(0)"))
        })
        .ok_or_else(|| {
            eyre!(
                "3Drep evaluate manifest for {context} has no numerator-only Q({edge_id},0) parameter"
            )
        })?;
    assert_eq!(
        parameter["source"].as_str(),
        Some("user override"),
        "Q({edge_id},0) parameter in {context} should be user supplied"
    );
    let value = parameter["value"]
        .as_str()
        .ok_or_else(|| eyre!("Q({edge_id},0) parameter in {context} has no value"))?;
    assert!(
        value.contains(&format!("{expected_value:+.17e}")),
        "Q({edge_id},0) parameter in {context} was {value}, expected {expected_value:+.17e}"
    );
    Ok(())
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

fn import_aa_aa_threedrep_graphs(
    test_name: &str,
    process_name: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    run_commands(
        &mut cli,
        &[
            "import model sm-default.json",
            r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
numerator_interpolation_scale = 2.0

[kinematics]
e_cm = 300.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
  [500.0, 500.0, 0.0, 0.0],
  "dependent",
]
helicities = [1, 1, -1, -1]
'"#,
        ],
    )?;
    set_threedrep_sampling_mode(&mut cli, "beyond_quadratic")?;
    cli.run_command(&format!(
        "import graphs {} -p {process_name} -i default -o",
        gammaloop_threedreps_dot_path("aa_aa.dot").display()
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
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_QUADRATIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cubic_edge_energy",
            dot_name: "box.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_CUBIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cubic_lmb_carrier_energy",
            dot_name: "box.dot",
            numerator: DENSE_LOOP0_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_CUBIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_4_0_0_0",
            dot_name: "box.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_HIGH_CONTACT_A,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_3_2_0_0",
            dot_name: "box.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_HIGH_CONTACT_B,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_2_1_0_3",
            dot_name: "box.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_HIGH_CONTACT_C,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "normal_box_cff_ltd_high_contact_0_0_3_3",
            dot_name: "box.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_QUADRATIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_single_cubic_nonrepeated_with_repeated_spectators",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_CUBIC_SPECTATOR,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_unsplit_repeated_higher_bounds",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_repeated_quintic",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_REPEATED_QUINTIC,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_mixed_quintic",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_AMBIGUOUS,
            seed: 1337,
            compare: ProbeComparison::BuildOnly,
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "box_pow3_uniform_beyond_quadratic_m_1",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
            uniform_scale: Some(1.0),
        },
        ProbeCase {
            name: "box_pow3_uniform_beyond_quadratic_m_minus_2",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_BOX_POW3_REPEATED_HIGH,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::BeyondQuadratic,
            uniform_scale: Some(-2.0),
        },
        ProbeCase {
            name: "box_pow3_uniform_all_m_2_75",
            dot_name: "box_pow3.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_SUNRISE_QUADRATIC_0,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "sunrise_pow4_free_lower_sector_quintic",
            dot_name: "sunrise_pow4.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_SUNRISE_QUINTIC_2,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "sunrise_pow4_repeated_channel_cubic_pair",
            dot_name: "sunrise_pow4.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_ITER_QUADRATIC_COMBO,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "iterated_sandwiched_cubic_linear",
            dot_name: "proper_iterated_sandwiched_bubble.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
            bounds: BOUNDS_ITER_CUBIC_LINEAR,
            seed: 1337,
            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        },
        ProbeCase {
            name: "iterated_sandwiched_cubic_pair",
            dot_name: "proper_iterated_sandwiched_bubble.dot",
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            numerator: DENSE_BOUNDS_NUMERATOR,
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
            compare: ProbeComparison::CffPureLtdBuildOnly,
            sampling_scale: NumeratorSamplingScaleMode::None,
            uniform_scale: None,
        });
    }

    cases
}

fn local_edge_energy_monomial(bounds: &[(usize, usize)]) -> String {
    let numerator = bounds
        .iter()
        .flat_map(|(edge_id, power)| {
            (0..*power).map(move |power_index| {
                format!(
                    "dot(edges[{edge_id}], ext[{}])",
                    (edge_id + power_index) % 4
                )
            })
        })
        .join(" * ");
    if numerator.is_empty() {
        "1".to_string()
    } else {
        numerator
    }
}

fn local_loop_energy_monomial(loop_id: usize, power: usize, external_count: usize) -> String {
    assert!(
        external_count > 0,
        "dense loop-energy numerator construction needs at least one external momentum"
    );
    (0..power)
        .map(|power_index| {
            format!(
                "dot(loops[{loop_id}], ext[{}])",
                (loop_id + power_index) % external_count
            )
        })
        .join(" * ")
}

fn dense_probe_numerator(
    numerator: &str,
    bounds: &[(usize, usize)],
    parsed: &ParsedGraph,
) -> Result<String> {
    if bounds.is_empty()
        || (numerator != DENSE_BOUNDS_NUMERATOR
            && numerator != DENSE_LOOP0_BOUNDS_NUMERATOR
            && !numerator.contains("[0]"))
        || parsed.external_names.is_empty()
    {
        return Ok(numerator.to_string());
    }
    if numerator == DENSE_LOOP0_BOUNDS_NUMERATOR || numerator.contains("loops[0][0]") {
        let degree = bounds
            .iter()
            .find_map(|(edge_id, degree)| (*edge_id == 0).then_some(*degree))
            .unwrap_or(0);
        return Ok(local_loop_energy_monomial(
            0,
            degree,
            parsed.external_names.len(),
        ));
    }
    auto_numerator_expr_for_bounds(
        parsed.external_names.len(),
        bounds,
        parsed.internal_edges.len(),
    )
    .map_err(Into::into)
}

fn run_probe_case(cli: &gammaloop_integration_tests::CLIState, case: ProbeCase) -> ProbeCaseReport {
    run_probe_case_parts(
        cli,
        ProbeCaseRun {
            name: case.name.to_string(),
            dot_name: case.dot_name,
            numerator: case.numerator.to_string(),
            bounds: case.bounds.to_vec(),
            seed: case.seed,
            compare: case.compare,
            sampling_scale_mode: case.sampling_scale,
            uniform_scale: case.uniform_scale,
        },
    )
}

struct ProbeCaseRun<'a> {
    name: String,
    dot_name: &'a str,
    numerator: String,
    bounds: Vec<(usize, usize)>,
    seed: u64,
    compare: ProbeComparison,
    sampling_scale_mode: NumeratorSamplingScaleMode,
    uniform_scale: Option<f64>,
}

fn run_probe_case_parts(
    cli: &gammaloop_integration_tests::CLIState,
    case: ProbeCaseRun<'_>,
) -> ProbeCaseReport {
    let ProbeCaseRun {
        name,
        dot_name,
        numerator,
        bounds,
        seed,
        compare,
        sampling_scale_mode,
        uniform_scale,
    } = case;
    let comparison = match compare {
        ProbeComparison::CffLtd { .. } => "cff_ltd",
        ProbeComparison::CffPureLtd { .. } => "cff_pureltd",
        ProbeComparison::BuildOnly => "build_only",
        ProbeComparison::CffPureLtdBuildOnly => "cff_pureltd_build_only",
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
        let numerator = dense_probe_numerator(&numerator, &bounds, &parsed)?;
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
            ProbeComparison::CffPureLtd { .. } | ProbeComparison::CffPureLtdBuildOnly => {
                RepresentationMode::PureLtd
            }
            _ => RepresentationMode::Ltd,
        };
        let rhs = generate_3d_expression_from_parsed(
            &parsed,
            &expression_options(rhs_representation, &bounds, sampling_scale_mode),
        );

        match compare {
            ProbeComparison::BuildOnly | ProbeComparison::CffPureLtdBuildOnly => {
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

fn source_external_edges(graph: &Graph) -> Vec<usize> {
    graph
        .underlying
        .iter_edges()
        .filter_map(|(pair, edge_id, _)| (!pair.is_paired()).then_some(edge_id.0))
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
    assert!(
        case["evaluator_build_status"]
            .as_str()
            .is_some_and(|status| status.starts_with("ok")),
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
fn cli_validate_build_and_evaluate_use_gammaloop_graph_state() -> Result<()> {
    let test_name = "threedreps_cli_validate_build_evaluate";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let validate_path = test_root.join("validate.json");
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    let dot_string = fs::read_to_string(gammaloop_threedreps_dot_path("box.dot"))?;
    cli.run_command(&format!(
        "import graphs --inline-dot {} -p threedreps_box_build_eval -i default -o",
        shell_single_quote(&dot_string)
    ))?;

    cli.run_command(&format!(
        "3Drep validate -p threedreps_box_build_eval -i default -g 0 --json-out {}",
        validate_path.display()
    ))?;
    let validation = serde_json::from_str::<JsonValue>(&fs::read_to_string(&validate_path)?)?;
    assert_eq!(validation["graph"]["n_internal_edges"].as_u64(), Some(4));
    assert_eq!(
        validation["graph"]["loop_names"].as_array().map(Vec::len),
        Some(1)
    );
    assert_eq!(validation["validation"]["ok"].as_bool(), Some(true));

    cli.run_command(&format!(
        "3Drep build -p threedreps_box_build_eval -i default -g 0 --workspace-path {} --no-pretty",
        workspace_path.display()
    ))?;
    let expression_path = latest_oriented_expression_path(&workspace_path)?;
    let artifact = serde_json::from_str::<JsonValue>(&fs::read_to_string(&expression_path)?)?;
    assert_eq!(artifact["family"].as_str(), Some("Cff"));
    assert_eq!(artifact["graph"]["n_internal_edges"].as_u64(), Some(4));
    assert_eq!(
        artifact["expression"]["orientations"]
            .as_array()
            .map(Vec::len),
        Some(14)
    );

    cli.run_command(&format!(
        "3Drep evaluate -p threedreps_box_build_eval -i default -g 0 --representation cff --workspace-path {} --standalone-rust-only",
        workspace_path.display()
    ))?;
    let expression_dir = expression_path
        .parent()
        .ok_or_else(|| eyre!("oriented expression path has no parent"))?;
    let symbolica_expression_path = expression_dir.join("symbolica_expression.txt");
    assert!(symbolica_expression_path.exists());
    assert!(
        !fs::read_to_string(&symbolica_expression_path)?
            .trim()
            .is_empty(),
        "canonical Symbolica expression should be written"
    );
    let symbolica_expression_pretty_path = expression_dir.join("symbolica_expression_pretty.txt");
    assert!(symbolica_expression_pretty_path.exists());
    assert!(
        !fs::read_to_string(&symbolica_expression_pretty_path)?
            .trim()
            .is_empty(),
        "pretty Symbolica expression should be written"
    );
    let raw_expression_path = expression_dir.join("symbolica_expression_raw.json");
    assert!(raw_expression_path.exists());
    let raw_expression =
        serde_json::from_str::<JsonValue>(&fs::read_to_string(&raw_expression_path)?)?;
    assert_eq!(raw_expression["schema_version"].as_u64(), Some(2));
    assert!(raw_expression["evaluator_backend"].as_str().is_some());
    let raw_script_path = expression_dir.join("symbolica_expression_raw.rs");
    assert!(raw_script_path.exists());
    assert!(
        fs::read_to_string(&raw_script_path)?.starts_with("#!/usr/bin/env -S rust-script --debug")
    );
    assert!(expression_dir.join("param_builder.txt").exists());
    assert!(
        !expression_dir.join("evaluate_manifest.json").exists(),
        "standalone raw replay mode should not write the evaluate manifest"
    );

    cli.run_command(&format!(
        "3Drep evaluate -p threedreps_box_build_eval -i default -g 0 --representation cff --workspace-path {}",
        workspace_path.display()
    ))?;
    let evaluate_manifest = serde_json::from_str::<JsonValue>(&fs::read_to_string(
        expression_dir.join("evaluate_manifest.json"),
    )?)?;
    assert_eq!(
        evaluate_manifest["symbolica_expression_path"].as_str(),
        symbolica_expression_path.to_str(),
        "3Drep evaluate manifest should record the canonical Symbolica expression path"
    );
    assert_eq!(
        evaluate_manifest["symbolica_expression_pretty_path"].as_str(),
        symbolica_expression_pretty_path.to_str(),
        "3Drep evaluate manifest should record the pretty Symbolica expression path"
    );
    assert_eq!(
        evaluate_manifest["symbolica_expression_raw_path"].as_str(),
        raw_expression_path.to_str(),
        "3Drep evaluate manifest should record the raw Symbolica evaluator input path"
    );
    assert_eq!(
        evaluate_manifest["symbolica_expression_raw_script_path"].as_str(),
        raw_script_path.to_str(),
        "3Drep evaluate manifest should record the raw Symbolica evaluator replay script path"
    );
    assert!(
        evaluate_manifest["evaluation"]["value"].as_str().is_some(),
        "3Drep evaluate manifest should record the evaluated value"
    );
    assert!(
        evaluate_manifest["evaluation"]["sample_evaluation_timing_seconds"]
            .as_f64()
            .is_some(),
        "3Drep evaluate manifest should record wall timing"
    );
    assert!(
        evaluate_manifest["parameters"]
            .as_array()
            .is_some_and(|parameters| !parameters.is_empty()),
        "3Drep evaluate manifest should record input parameters"
    );

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
        clean_test(get_tests_workspace_path().join(test_name));
        let bootstrap = import_threedrep_graph(
            &format!("{test_name}_bootstrap"),
            &gammaloop_threedreps_dot_path(dot_name),
            "bootstrap",
        )?;
        let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
        let bootstrap_source_internal_edges = source_internal_edges(bootstrap_graph);
        let bootstrap_source_external_edges = source_external_edges(bootstrap_graph);
        let numerator = if let Some((loop_id, degree)) = lmb_bound {
            dense_loop_energy_power_numerator(loop_id, degree, &bootstrap_source_external_edges)
        } else {
            dense_edge_energy_power_numerator(
                &bootstrap_source_internal_edges,
                &bootstrap_source_external_edges,
                &local_bounds,
            )
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
        let manifest_path = find_named_artifact(&workspace_path, "test_cff_ltd_manifest.json")?;
        let manifest = serde_json::from_str::<JsonValue>(&fs::read_to_string(manifest_path)?)?;
        assert_threedrep_comparison_success(&manifest, test_name)?;
        assert_threedrep_evaluation_ids(&manifest, test_name)?;
        let graph = imported_graph_from_cli(&cli, &process_name, "default", 0)?;
        let expected = graph
            .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;
        assert_eq!(
            serde_json::from_value::<Vec<(usize, usize)>>(
                manifest["automatic_energy_degree_bounds"].clone()
            )?,
            expected,
            "{test_name} automatic energy bounds should match canonical graph-level analysis for numerator {numerator}"
        );

        clean_test(get_tests_workspace_path().join(test_name));
        clean_test(&cli.cli_settings.state.folder);
    }

    Ok(())
}

#[test]
#[serial]
fn cli_box_lower_contact_high_powers_match_ltd() -> Result<()> {
    let test_name = "threedreps_box_lower_contact_high_powers";
    clean_test(get_tests_workspace_path().join(test_name));
    let bootstrap = import_threedrep_graph(
        &format!("{test_name}_bootstrap"),
        &gammaloop_threedreps_dot_path("box.dot"),
        "bootstrap",
    )?;
    let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
    let source_internal_edges = source_internal_edges(bootstrap_graph);
    let source_external_edges = source_external_edges(bootstrap_graph);
    clean_test(&bootstrap.cli_settings.state.folder);
    clean_test(get_tests_workspace_path().join(format!("{test_name}_bootstrap")));

    let cases = [
        ("single_quartic", vec![(0usize, 4usize)], 1337, 1.0e-1),
        ("single_quintic", vec![(3usize, 5usize)], 1337, 1.0e-1),
        ("single_sextic", vec![(0usize, 6usize)], 1337, 1.0e-1),
        (
            "cubic_quadratic_lower_contact",
            vec![(0usize, 3usize), (1, 2)],
            1337,
            1.0e-1,
        ),
        (
            "quadratic_linear_cubic_lower_contact",
            vec![(0usize, 2usize), (1, 1), (3, 3)],
            1337,
            1.0e-1,
        ),
        (
            "quartic_linear_terminal",
            vec![(2usize, 4usize), (3, 1)],
            1337,
            1.0e-1,
        ),
        (
            "quartic_quadratic_terminal",
            vec![(0usize, 4usize), (1, 2)],
            1337,
            1.0e-1,
        ),
    ];

    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    set_threedrep_compile_backend(&mut cli, None)?;
    set_threedrep_sampling_mode(&mut cli, "none")?;

    for (case_name, local_bounds, seed, scale) in cases {
        let numerator = dense_edge_energy_power_numerator(
            &source_internal_edges,
            &source_external_edges,
            &local_bounds,
        );
        let dot_path =
            dot_with_global_numerator(&format!("{test_name}_{case_name}"), "box.dot", &numerator)?;
        let process_name = format!("{test_name}_{case_name}");
        cli.run_command(&format!(
            "import graphs {} -p {process_name} -i default -o",
            dot_path.display()
        ))?;
        if let Some(dot_dir) = dot_path.parent() {
            clean_test(dot_dir);
        }
        let case_workspace = workspace_path.join(case_name);
        let manifest = run_threedrep_test_cff_ltd_manifest(
            &mut cli,
            &case_workspace,
            &process_name,
            "Double",
            seed,
            scale,
        )?;
        let graph = imported_graph_from_cli(&cli, &process_name, "default", 0)?;
        let expected_automatic_bounds = graph
            .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;
        assert_eq!(
            serde_json::from_value::<Vec<(usize, usize)>>(
                manifest["automatic_energy_degree_bounds"].clone()
            )?,
            expected_automatic_bounds,
            "{case_name} should infer high-power energy caps from the graph numerator"
        );
        let cff_value = direct_manifest_evaluation_value(&manifest, "cff")?;
        let ltd_value = direct_manifest_evaluation_value(&manifest, "ltd")?;
        assert_manifest_direct_values_agree(
            &ltd_value,
            &cff_value,
            &format!("LTD vs CFF for lower-contact high-power box case {case_name}"),
            1.0e-7,
        );
    }

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_physical_hexagon_high_power_energy_sectors_match_ltd() -> Result<()> {
    let test_name = "threedreps_physical_hexagon_high_power_sectors";
    clean_test(get_tests_workspace_path().join(test_name));

    let bootstrap = import_threedrep_graph(
        &format!("{test_name}_bootstrap"),
        &gammaloop_threedreps_dot_path("one_loop_6_external.dot"),
        "bootstrap",
    )?;
    let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
    let source_internal_edges = source_internal_edges(bootstrap_graph);
    let source_external_edges = source_external_edges(bootstrap_graph);
    let q0 = source_internal_edges[0];
    let q1 = source_internal_edges[1];
    let p0 = source_external_edges[1];
    clean_test(&bootstrap.cli_settings.state.folder);
    clean_test(get_tests_workspace_path().join(format!("{test_name}_bootstrap")));

    let quadratic_direct = format!(
        "{}*{}",
        energy_component_power(q0, 2),
        energy_component_power(q1, 2)
    );
    let quadratic_direct_bounds = format!("{q0}:2,{q1}:2");
    let cubic_direct = format!(
        "{}*{}",
        energy_component_power(q0, 3),
        energy_component_power(q1, 3)
    );
    let cubic_direct_bounds = format!("{q0}:3,{q1}:3");
    let q1_from_q0 = format!("({}-{})", q0_component(q0), q0_component(p0));
    let quadratic_traded = format!(
        "{}*{}*{}",
        energy_component_power(q0, 2),
        q1_from_q0,
        q1_from_q0
    );
    let quadratic_traded_bounds = format!("{q0}:4");
    let cases = [
        (
            "quadratic_direct",
            quadratic_direct.as_str(),
            quadratic_direct_bounds.as_str(),
            1.0e-7,
        ),
        (
            "cubic_direct",
            cubic_direct.as_str(),
            cubic_direct_bounds.as_str(),
            1.0e-7,
        ),
        (
            "quadratic_traded",
            quadratic_traded.as_str(),
            quadratic_traded_bounds.as_str(),
            1.0e-7,
        ),
    ];

    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    set_threedrep_compile_backend(&mut cli, None)?;
    set_threedrep_sampling_mode(&mut cli, "none")?;

    for (case_name, numerator, energy_degree_bounds, tolerance) in cases {
        let dot_path = dot_with_global_numerator(
            &format!("{test_name}_{case_name}"),
            "one_loop_6_external.dot",
            numerator,
        )?;
        let process_name = format!("{test_name}_{case_name}");
        cli.run_command(&format!(
            "import graphs {} -p {process_name} -i default -o",
            dot_path.display()
        ))?;
        if let Some(dot_dir) = dot_path.parent() {
            clean_test(dot_dir);
        }
        let case_workspace = workspace_path.join(case_name);
        let manifest = run_threedrep_test_cff_ltd_manifest_for_graph(
            &mut cli,
            &case_workspace,
            &process_name,
            0,
            "Double",
            291,
            1.0e-1,
            Some(energy_degree_bounds),
        )?;
        let cff_value = direct_manifest_evaluation_value(&manifest, "cff")?;
        let ltd_value = direct_manifest_evaluation_value(&manifest, "ltd")?;
        assert_manifest_direct_values_agree(
            &ltd_value,
            &cff_value,
            &format!("LTD vs CFF for physical hexagon {case_name}"),
            tolerance,
        );
    }

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_multiloop_high_power_energy_sectors_match_ltd() -> Result<()> {
    run_on_large_stack(
        "cli_multiloop_high_power_energy_sectors_match_ltd",
        cli_multiloop_high_power_energy_sectors_match_ltd_inner,
    )
}

fn cli_multiloop_high_power_energy_sectors_match_ltd_inner() -> Result<()> {
    let test_name = "threedreps_multiloop_high_power_sectors";
    clean_test(get_tests_workspace_path().join(test_name));

    let cases: [(&str, &str, &[(usize, usize)], u64, f64); 3] = [
        (
            "sunrise_pow4_free_lower_sector_quintic",
            "sunrise_pow4.dot",
            BOUNDS_SUNRISE_QUINTIC_2,
            1337,
            1.0e-1,
        ),
        (
            "sunrise_pow4_repeated_channel_cubic_pair",
            "sunrise_pow4.dot",
            BOUNDS_SUNRISE_CUBIC_PAIR,
            123,
            1.0e-1,
        ),
        (
            "four_loop_stress_quadratic",
            "four_loop_stress.dot",
            BOUNDS_FOUR_LOOP_QUADRATIC,
            1337,
            1.0e-1,
        ),
    ];

    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    set_threedrep_compile_backend(&mut cli, None)?;
    set_threedrep_sampling_mode(&mut cli, "none")?;

    for (case_name, dot_name, local_bounds, seed, scale) in cases {
        let bootstrap = import_threedrep_graph(
            &format!("{test_name}_{case_name}_bootstrap"),
            &gammaloop_threedreps_dot_path(dot_name),
            "bootstrap",
        )?;
        let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
        let source_internal_edges = source_internal_edges(bootstrap_graph);
        let source_external_edges = source_external_edges(bootstrap_graph);
        let numerator = dense_edge_energy_power_numerator(
            &source_internal_edges,
            &source_external_edges,
            local_bounds,
        );
        let energy_degree_bounds = source_bounds_from_local(&source_internal_edges, local_bounds);
        clean_test(&bootstrap.cli_settings.state.folder);
        clean_test(get_tests_workspace_path().join(format!("{test_name}_{case_name}_bootstrap")));

        let dot_path =
            dot_with_global_numerator(&format!("{test_name}_{case_name}"), dot_name, &numerator)?;
        let process_name = format!("{test_name}_{case_name}");
        cli.run_command(&format!(
            "import graphs {} -p {process_name} -i default -o",
            dot_path.display()
        ))?;
        if let Some(dot_dir) = dot_path.parent() {
            clean_test(dot_dir);
        }

        let manifest = run_threedrep_test_cff_ltd_manifest_for_graph(
            &mut cli,
            &workspace_path.join(case_name),
            &process_name,
            0,
            "Double",
            seed,
            scale,
            Some(&energy_degree_bounds),
        )?;
        let cff_value = direct_manifest_evaluation_value(&manifest, "cff")?;
        let ltd_value = direct_manifest_evaluation_value(&manifest, "ltd")?;
        assert_manifest_direct_values_agree(
            &ltd_value,
            &cff_value,
            &format!("LTD vs CFF for multiloop high-power case {case_name}"),
            1.0e-7,
        );
    }

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_box_pow3_high_power_agrees_for_all_numerator_sample_normalizations() -> Result<()> {
    let test_name = "threedreps_box_pow3_high_power_sample_normalizations";
    clean_test(get_tests_workspace_path().join(test_name));
    let bootstrap = import_threedrep_graph(
        &format!("{test_name}_bootstrap"),
        &gammaloop_threedreps_dot_path("box_pow3.dot"),
        "bootstrap",
    )?;
    let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
    let source_internal_edges = source_internal_edges(bootstrap_graph);
    let source_external_edges = source_external_edges(bootstrap_graph);
    let numerator = dense_edge_energy_power_numerator(
        &source_internal_edges,
        &source_external_edges,
        BOUNDS_BOX_POW3_REPEATED_HIGH,
    );
    clean_test(&bootstrap.cli_settings.state.folder);

    let dot_path = dot_with_global_numerator(test_name, "box_pow3.dot", &numerator)?;
    let process_name = "threedreps_box_pow3_high_power_norms";
    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command("set default-runtime kv general.numerator_interpolation_scale=2.0")?;
    set_threedrep_compile_backend(&mut cli, None)?;
    cli.run_command(&format!(
        "import graphs {} -p {process_name} -i default -o",
        dot_path.display()
    ))?;

    let modes = [
        ("none", "NeverM"),
        ("beyond_quadratic", "MForBeyondQuadraticOnly"),
        ("all", "MForAll"),
    ];
    let mut cff_reference = None;
    for (mode, expected_manifest_mode) in modes {
        set_threedrep_sampling_mode(&mut cli, mode)?;
        let mode_workspace = workspace_path.join(format!("mode_{mode}"));
        let manifest = run_threedrep_test_cff_ltd_manifest(
            &mut cli,
            &mode_workspace,
            process_name,
            "Double",
            23,
            0.25,
        )?;
        assert_eq!(
            manifest["settings"]["numerator_samples_normalization"].as_str(),
            Some(expected_manifest_mode),
            "test-cff-ltd should use the active global numerator normalization mode"
        );
        assert_eq!(
            manifest["mass_shift_start"].as_f64(),
            Some(0.25),
            "default pure-LTD mass shift should start at 1.0 * --scale"
        );
        assert_eq!(
            serde_json::from_value::<Vec<(usize, usize)>>(
                manifest["override_energy_degree_bounds"].clone()
            )?,
            Vec::<(usize, usize)>::new(),
            "the high-power bounds should come from automatic numerator analysis"
        );
        let graph = imported_graph_from_cli(&cli, process_name, "default", 0)?;
        let expected_automatic_bounds = graph
            .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;
        assert_eq!(
            serde_json::from_value::<Vec<(usize, usize)>>(
                manifest["automatic_energy_degree_bounds"].clone()
            )?,
            expected_automatic_bounds,
            "dense four-vector numerator should use the canonical graph-level energy caps"
        );
        let cff_value = direct_manifest_evaluation_value(&manifest, "cff")?;
        let ltd_value = direct_manifest_evaluation_value(&manifest, "ltd")?;
        assert_manifest_direct_values_agree(
            &ltd_value,
            &cff_value,
            &format!("LTD vs CFF for normalization mode {mode}"),
            1.0e-7,
        );
        if let Some(reference) = &cff_reference {
            assert_manifest_direct_values_agree(
                &cff_value,
                reference,
                &format!("CFF value for normalization mode {mode} vs baseline"),
                1.0e-8,
            );
        } else {
            cff_reference = Some(cff_value);
        }
    }

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_box_pow3_high_power_compiled_double_backends_agree_with_eager() -> Result<()> {
    let test_name = "threedreps_box_pow3_high_power_compiled_backends";
    clean_test(get_tests_workspace_path().join(test_name));
    let bootstrap = import_threedrep_graph(
        &format!("{test_name}_bootstrap"),
        &gammaloop_threedreps_dot_path("box_pow3.dot"),
        "bootstrap",
    )?;
    let bootstrap_graph = imported_graph_from_cli(&bootstrap, "bootstrap", "default", 0)?;
    let source_internal_edges = source_internal_edges(bootstrap_graph);
    let source_external_edges = source_external_edges(bootstrap_graph);
    let numerator = dense_edge_energy_power_numerator(
        &source_internal_edges,
        &source_external_edges,
        BOUNDS_BOX_POW3_REPEATED_HIGH,
    );
    clean_test(&bootstrap.cli_settings.state.folder);

    let dot_path = dot_with_global_numerator(test_name, "box_pow3.dot", &numerator)?;
    let process_name = "threedreps_box_pow3_high_power_backends";
    let state_path = get_tests_workspace_path().join(test_name).join("state");
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command("set default-runtime kv general.numerator_interpolation_scale=2.0")?;
    set_threedrep_sampling_mode(&mut cli, "all")?;
    cli.run_command(&format!(
        "import graphs {} -p {process_name} -i default -o",
        dot_path.display()
    ))?;

    set_threedrep_compile_backend(&mut cli, None)?;
    let eager_manifest = run_threedrep_test_cff_ltd_manifest(
        &mut cli,
        &workspace_path.join("backend_eager"),
        process_name,
        "Double",
        23,
        0.25,
    )?;
    assert_manifest_uses_backend(&eager_manifest, "eager", "eager baseline")?;
    let eager_cff = direct_manifest_evaluation_value(&eager_manifest, "cff")?;
    let eager_ltd = direct_manifest_evaluation_value(&eager_manifest, "ltd")?;
    assert_manifest_direct_values_agree(&eager_ltd, &eager_cff, "eager LTD vs CFF", 1.0e-7);

    for backend in ["symjit", "assembly"] {
        set_threedrep_compile_backend(&mut cli, Some(backend))?;
        let manifest = run_threedrep_test_cff_ltd_manifest(
            &mut cli,
            &workspace_path.join(format!("backend_{backend}")),
            process_name,
            "Double",
            23,
            0.25,
        )?;
        assert_manifest_uses_backend(&manifest, backend, backend)?;
        let compiled_cff = direct_manifest_evaluation_value(&manifest, "cff")?;
        let compiled_ltd = direct_manifest_evaluation_value(&manifest, "ltd")?;
        assert_manifest_direct_values_agree(
            &compiled_ltd,
            &compiled_cff,
            &format!("{backend} LTD vs CFF"),
            1.0e-7,
        );
        assert_manifest_direct_values_agree(
            &compiled_cff,
            &eager_cff,
            &format!("{backend} CFF vs eager CFF"),
            1.0e-8,
        );
        assert!(
            !compiled_cff.backend.is_empty(),
            "compiled backend record should not be empty"
        );
    }

    clean_test(get_tests_workspace_path().join(test_name));
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_aa_aa_box_and_double_box_multibackend_threedrep_comparisons() -> Result<()> {
    let test_name = "threedreps_aa_aa_box_double_box_multibackend";
    let test_root = get_tests_workspace_path().join(test_name);
    clean_test(&test_root);
    let process_name = "threedreps_aa_aa_multibackend";
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = import_aa_aa_threedrep_graphs(test_name, process_name)?;
    set_threedrep_sampling_mode(&mut cli, "none")?;

    let graph_cases = [(0usize, 1usize, "box"), (1usize, 2usize, "double_box")];
    for (graph_id, loop_count, label) in graph_cases {
        let graph = imported_graph_from_cli(&cli, process_name, "default", graph_id)?;
        let parsed = graph.to_three_d_parsed_graph()?;
        assert_eq!(
            parsed.loop_names.len(),
            loop_count,
            "aa_aa {label} should have {loop_count} loop(s)"
        );
    }

    let backend_cases = [
        ("double_assembly", "Double", Some("assembly"), "assembly"),
        ("double_symjit", "Double", Some("symjit"), "symjit"),
        ("quad_eager", "Quad", None, "eager"),
    ];
    let comparison_graph_cases = [(0usize, "box")];
    for (graph_id, graph_label) in comparison_graph_cases {
        for (case_label, precision, backend, expected_backend) in backend_cases {
            set_threedrep_compile_backend(&mut cli, backend)?;
            let manifest = run_threedrep_test_cff_ltd_manifest_for_graph(
                &mut cli,
                &workspace_path.join(format!("graph_{graph_id}_{case_label}")),
                process_name,
                graph_id,
                precision,
                11,
                0.25,
                None,
            )?;
            assert_eq!(
                manifest["graph_id"].as_u64(),
                Some(graph_id as u64),
                "aa_aa {graph_label} {case_label} should target the requested graph"
            );
            assert_manifest_uses_backend(
                &manifest,
                expected_backend,
                &format!("aa_aa {graph_label} {case_label}"),
            )?;
        }
    }

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_aa_aa_double_box_energy_bounds_track_edge_energies() -> Result<()> {
    let test_name = "threedreps_aa_aa_double_box_energy_bounds";
    let test_root = get_tests_workspace_path().join(test_name);
    clean_test(&test_root);
    let process_name = "threedreps_aa_aa_double_box_bounds";
    let mut cli = import_aa_aa_threedrep_graphs(test_name, process_name)?;
    let graph = imported_graph_from_cli(&cli, process_name, "default", 1)?;
    let ltd_bounds = graph
        .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Ltd)?;
    let cff_bounds = graph
        .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;

    assert_eq!(
        ltd_bounds,
        vec![(4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (9, 1)],
        "aa_aa double-box LTD bounds should stay attached to the six fermion edge energies"
    );
    assert_eq!(
        cff_bounds, ltd_bounds,
        "aa_aa double-box CFF bounds should also stay attached to source edge energies instead of collapsing onto loop-carrier powers"
    );
    assert!(
        cff_bounds.iter().all(|(_, degree)| *degree <= 1),
        "aa_aa double-box CFF should not infer beyond-linear source-edge caps"
    );

    for (representation, expected_bounds) in [("ltd", ltd_bounds), ("cff", cff_bounds)] {
        let workspace = test_root.join(representation);
        cli.run_command(&format!(
            "3Drep build -p {process_name} -i default -g 1 --representation {representation} --numerator-samples-normalization M_for_beyond_quadratic_only --workspace-path {} --no-pretty",
            workspace.display()
        ))?;
        let expression_path = find_named_artifact(&workspace, "oriented_expression.json")?;
        let expression = serde_json::from_str::<JsonValue>(&fs::read_to_string(expression_path)?)?;
        assert_eq!(
            serde_json::from_value::<Vec<(usize, usize)>>(
                expression["automatic_energy_degree_bounds"].clone()
            )?,
            expected_bounds,
            "3Drep build should persist the graph-level {representation} energy bounds"
        );
    }

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_aa_aa_evaluate_reuses_standard_and_numerator_only_evaluator_caches() -> Result<()> {
    let test_name = "threedreps_aa_aa_evaluate_cache_reuse";
    let test_root = get_tests_workspace_path().join(test_name);
    clean_test(&test_root);
    let process_name = "threedreps_aa_aa_evaluate_cache";
    let workspace_path = imported_graph_threedrep_workspace(test_name);
    let mut cli = import_aa_aa_threedrep_graphs(test_name, process_name)?;
    set_threedrep_compile_backend(&mut cli, None)?;

    for (graph_id, graph_label) in [(0usize, "box")] {
        let standard_workspace = workspace_path.join(format!("graph_{graph_id}_standard"));
        cli.run_command(&format!(
            "3Drep build -p {process_name} -i default -g {graph_id} --representation cff --numerator-samples-normalization M_for_beyond_quadratic_only --workspace-path {} --no-pretty --clean",
            standard_workspace.display()
        ))?;
        let evaluate_command = format!(
            "3Drep evaluate -p {process_name} -i default -g {graph_id} --representation cff --numerator-samples-normalization M_for_beyond_quadratic_only --workspace-path {} --precision Double --seed 11 --scale 0.25 --eager",
            standard_workspace.display()
        );
        cli.run_command(&format!("{evaluate_command} --clean"))?;
        let first_standard = read_evaluate_manifest(&standard_workspace)?;
        assert_evaluate_manifest_ok(
            &first_standard,
            &format!("aa_aa {graph_label} standard first build"),
        )?;
        assert!(
            first_standard["evaluation"]["evaluator_build_timing_seconds"]
                .as_f64()
                .is_some(),
            "first standard aa_aa {graph_label} evaluation should build an evaluator"
        );
        let standard_value = first_standard["evaluation"]["value"]
            .as_str()
            .ok_or_else(|| eyre!("standard aa_aa {graph_label} evaluation has no value"))?
            .to_string();
        find_named_artifact(&standard_workspace, "evaluate.evaluator.bin")?;

        cli.run_command(&evaluate_command)?;
        let reused_standard = read_evaluate_manifest(&standard_workspace)?;
        assert_evaluate_reused_cached_evaluator(
            &reused_standard,
            &format!("aa_aa {graph_label} standard cached evaluation"),
        )?;
        assert_eq!(
            reused_standard["evaluation"]["value"].as_str(),
            Some(standard_value.as_str()),
            "cached standard aa_aa {graph_label} evaluation should reproduce the same value"
        );

        let numerator_only_workspace =
            workspace_path.join(format!("graph_{graph_id}_numerator_only"));
        let numerator_only_command = format!(
            "3Drep evaluate -p {process_name} -i default -g {graph_id} --representation cff --numerator-samples-normalization M_for_beyond_quadratic_only --workspace-path {} --precision Double --seed 11 --scale 0.25 --eager --numerator-only --numerator-q0 4:1.25e-1",
            numerator_only_workspace.display()
        );
        cli.run_command(&format!("{numerator_only_command} --clean"))?;
        let first_numerator_only = read_evaluate_manifest(&numerator_only_workspace)?;
        assert_evaluate_manifest_ok(
            &first_numerator_only,
            &format!("aa_aa {graph_label} numerator-only first build"),
        )?;
        assert!(
            first_numerator_only["evaluation"]["evaluator_build_timing_seconds"]
                .as_f64()
                .is_some(),
            "first numerator-only aa_aa {graph_label} evaluation should build an evaluator"
        );
        assert_numerator_only_q0_override(
            &first_numerator_only,
            4,
            0.125,
            &format!("aa_aa {graph_label} numerator-only first build"),
        )?;
        let numerator_only_value = first_numerator_only["evaluation"]["value"]
            .as_str()
            .ok_or_else(|| eyre!("numerator-only aa_aa {graph_label} evaluation has no value"))?
            .to_string();
        find_named_artifact(&numerator_only_workspace, "numerator_only.evaluator.bin")?;

        cli.run_command(&numerator_only_command)?;
        let reused_numerator_only = read_evaluate_manifest(&numerator_only_workspace)?;
        assert_evaluate_reused_cached_evaluator(
            &reused_numerator_only,
            &format!("aa_aa {graph_label} numerator-only cached evaluation"),
        )?;
        assert_numerator_only_q0_override(
            &reused_numerator_only,
            4,
            0.125,
            &format!("aa_aa {graph_label} numerator-only cached evaluation"),
        )?;
        assert_eq!(
            reused_numerator_only["evaluation"]["value"].as_str(),
            Some(numerator_only_value.as_str()),
            "cached numerator-only aa_aa {graph_label} evaluation should reproduce the same value"
        );
    }

    clean_test(test_root);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn cli_graph_from_signatures_closes_without_explicit_lmb_rep() -> Result<()> {
    let cases: &[SignatureClosureCase<'_>] = &[
        (
            "threedreps_graph_from_signatures_box_dotted",
            &[
                ("k1", "1"),
                ("k1", "11"),
                ("k1+p1", "2"),
                ("k1+p1+p2", "3"),
                ("k1+p1+p2+p3", "4"),
            ],
            &["p"],
            5,
        ),
        (
            "threedreps_graph_from_signatures_pq_mass",
            &[
                ("k1+p1", "0"),
                ("k1+p1-q1", "0"),
                ("k1-p2+q2", "0"),
                ("k1", "0"),
            ],
            &["p", "q"],
            4,
        ),
        (
            "threedreps_graph_from_signatures_two_loop_sunrise",
            &[("k1", "1"), ("k2", "2"), ("k1+k2+p1", "3")],
            &["p"],
            2,
        ),
        (
            "threedreps_graph_from_signatures_three_loop_mercedes",
            &[
                ("k1", "0"),
                ("k2", "0"),
                ("k3", "0"),
                ("k1+k2+p1", "0"),
                ("k2+k3+p2", "0"),
                ("k1-k3+p3-p4", "0"),
            ],
            &["p"],
            4,
        ),
        (
            "threedreps_graph_from_signatures_four_loop_banana",
            &[
                ("k1", "1"),
                ("k2", "2"),
                ("k3", "3"),
                ("k4", "4"),
                ("k1+k2+k3+k4+p1", "5"),
            ],
            &["p"],
            2,
        ),
        (
            "threedreps_graph_from_signatures_five_loop_banana",
            &[
                ("k1", "1"),
                ("k2", "2"),
                ("k3", "3"),
                ("k4", "4"),
                ("k5", "5"),
                ("k1+k2+k3+k4+k5+p1", "6"),
            ],
            &["p"],
            2,
        ),
    ];

    for (test_name, factors, external_prefixes, num_vertices) in cases {
        assert_graph_from_signatures_cli_closure(
            test_name,
            factors,
            external_prefixes,
            *num_vertices,
        )?;
    }

    Ok(())
}

#[test]
#[serial]
fn cli_graph_from_signatures_rejects_inconsistent_explicit_vertex_count() -> Result<()> {
    let test_name = "threedreps_graph_from_signatures_inconsistent_vertices";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let dot_path = test_root.join("from_signatures.dot");
    fs::create_dir_all(&test_root)?;

    let signature_expression = signature_expression_from_factors(&[
        ("k1", "1"),
        ("k1", "11"),
        ("k1+p1", "2"),
        ("k1+p1+p2", "3"),
        ("k1+p1+p2+p3", "4"),
    ]);
    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    let result = cli.run_command(&format!(
        "3Drep graph-from-signatures --signatures {} --dot-output {} --external-prefixes p --num-vertices 4",
        shell_single_quote(&signature_expression),
        dot_path.display(),
    ));
    assert!(
        result
            .as_ref()
            .err()
            .is_some_and(|error| error.to_string().contains("requires 5 internal vertices")),
        "expected inconsistent vertex-count error, got {result:?}",
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
                ProbeCaseRun {
                    name: format!(
                        "old_python_all_edge_linear_numerator::{dot_name}::edge{edge_id}",
                        dot_name = spec.dot_name
                    ),
                    dot_name: spec.dot_name,
                    numerator,
                    bounds: Vec::new(),
                    seed: 1337,
                    compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
                    sampling_scale_mode: NumeratorSamplingScaleMode::None,
                    uniform_scale: None,
                },
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
                        ProbeCaseRun {
                            name: format!(
                                "old_python_normal_box_all_convergent_bounds::{e0}_{e1}_{e2}_{e3}"
                            ),
                            dot_name: "box.dot",
                            numerator: local_edge_energy_monomial(&bounds),
                            bounds,
                            seed: 1337,
                            compare: ProbeComparison::CffLtd { tolerance: 1.0e-8 },
                            sampling_scale_mode: NumeratorSamplingScaleMode::None,
                            uniform_scale: None,
                        },
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
    assert_eq!(
        failed_count, 0,
        "all old Python 3Drep parity cases should pass; see {report_path:?}"
    );

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
    let expected_automatic_bounds =
        imported_graph_from_cli(&run.cli, "threedreps_box", "default", 0)?
            .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;
    assert_eq!(
        serde_json::from_value::<Vec<(usize, usize)>>(manifest["energy_degree_bounds"].clone())?,
        expected_automatic_bounds
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

    assert_manifest_case(manifest, "cff", 14, 12)?;
    assert_manifest_case(manifest, "ltd", 4, 24)?;
    assert_manifest_case(manifest, "pureltd", 4, 24)?;

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
    let expected_automatic_bounds =
        imported_graph_from_cli(&run.cli, "threedreps_box_pow3", "default", 0)?
            .automatic_numerator_energy_degree_bounds_for_3d_expression(RepresentationMode::Cff)?;
    assert_eq!(
        serde_json::from_value::<Vec<(usize, usize)>>(manifest["energy_degree_bounds"].clone())?,
        expected_automatic_bounds
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

    assert_manifest_case(manifest, "cff", 62, 27)?;
    assert_manifest_case(manifest, "ltd", 17, 24)?;
    assert_manifest_case(manifest, "pureltd", 6, 57)?;

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
        [
            "+xxxxx", "x+xxxx", "xx+xxx", "xxx+++", "xxx---", "xxx--+", "xxx--+", "xxx-+-",
            "xxx-+-", "xxx-++", "xxx-++", "xxx+--", "xxx+--", "xxx+-+", "xxx+-+", "xxx++-",
            "xxx++-",
        ]
    );
    let ltd_half_edges = variant_half_edges(&run.ltd.expression);
    assert_eq!(
        &ltd_half_edges[..3],
        [vec![vec![4]], vec![vec![5]], vec![vec![6]]]
    );
    let repeated_ltd_half_edges = vec![vec![7, 7, 7], vec![7, 7, 7, 7], vec![7, 7, 7, 7, 7]];
    assert!(
        ltd_half_edges[3..]
            .iter()
            .all(|half_edges| half_edges == &repeated_ltd_half_edges),
        "unexpected repeated LTD half-edge structure: {ltd_half_edges:?}"
    );
    let ltd_denominator_ids = variant_linear_denominator_surface_ids(&run.ltd.expression);
    assert_eq!(
        &ltd_denominator_ids[..3],
        [
            vec![vec![1, 3, 5]],
            vec![vec![7, 9, 11]],
            vec![vec![13, 15, 17]],
        ]
    );
    let repeated_ltd_denominator_ids = vec![vec![19, 21, 23]; 3];
    assert!(
        ltd_denominator_ids[3..]
            .iter()
            .all(|ids| ids == &repeated_ltd_denominator_ids),
        "unexpected repeated LTD denominator surfaces: {ltd_denominator_ids:?}"
    );
    assert!(
        variant_prefactors(&run.ltd.expression)[3..]
            .iter()
            .all(|prefactors| prefactors.len() == 3),
        "repeated LTD orientations should expose the fused three-variant derivative structure"
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
        (0.05, 0.339_106_080_530_878_7),
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
