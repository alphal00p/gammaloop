use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    time::Instant,
};

use gammalooprs::{
    initialisation::test_initialise,
    numerator::{ParsingNet, aind::Aind},
    utils::{F, FUN_LIB, TENSORLIB},
};
use spenso::{
    iterators::IteratableTensor,
    network::library::symbolic::ExplicitKey,
    network::{
        ExecutionResult, MinResultRank, Network, Sequential, SmallestDegree, Steps,
        parsing::{ParseSettings, ShorthandParsing, StructureInferenceMode},
        store::NetworkStore,
    },
    structure::{HasName, HasStructure, TensorStructure},
    tensors::{data::DataTensor, parametric::ParamOrConcrete},
};
use symbolica::{
    atom::{Atom, AtomCore, Indeterminate, Symbol},
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};

type ConcreteParsingNet = Network<
    NetworkStore<
        spenso::tensors::parametric::MixedTensor<
            F<f64>,
            spenso::network::parsing::ShadowedStructure<Aind>,
        >,
        Atom,
    >,
    ExplicitKey<Aind>,
    Symbol,
    Aind,
>;

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("gammalooprs crate should live two levels below the workspace root")
        .to_path_buf()
}

fn report(label: &str, start: Instant) {
    eprintln!("{label}: elapsed={:.3?}", start.elapsed());
}

fn strip_symbol_tags(input: &str) -> String {
    let mut stripped = String::with_capacity(input.len());
    let mut rest = input;

    while let Some(start) = rest.find("::{") {
        let (prefix, after_prefix) = rest.split_at(start);
        stripped.push_str(prefix);

        let Some(end) = after_prefix.find("}::") else {
            stripped.push_str(after_prefix);
            return stripped;
        };

        stripped.push_str("::");
        rest = &after_prefix[end + "}::".len()..];
    }

    stripped.push_str(rest);
    stripped
}

fn parse_root_input(filename: &str) -> Atom {
    let path = workspace_root().join(filename);
    let start = Instant::now();
    let input = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    report(
        &format!("read {} bytes={}", path.display(), input.len()),
        start,
    );

    let start = Instant::now();
    let input = strip_symbol_tags(&input);
    report(&format!("strip_symbol_tags bytes={}", input.len()), start);

    let start = Instant::now();
    let expr = Atom::parse(wrap_input!(&input), SymbolicaParseSettings::default())
        .expect("failed to parse Symbolica expression");
    eprintln!(
        "atom stats: terms={} bytes={}",
        expr.nterms(),
        expr.as_view().get_byte_size()
    );
    report("symbolica_parse", start);

    expr
}

fn collect_horner(expr: &Atom) -> Atom {
    let start = Instant::now();
    let collected = expr.collect_horner::<Indeterminate>(None);
    eprintln!(
        "collect_horner stats: terms={} bytes={}",
        collected.nterms(),
        collected.as_view().get_byte_size()
    );
    report("collect_horner", start);
    collected
}

fn report_actual_net_stats(label: &str, net: &ParsingNet) {
    let mut concrete_tensors = 0usize;
    let mut param_tensors = 0usize;
    let mut concrete_entries = 0usize;
    let mut param_entries = 0usize;
    let mut param_dense = 0usize;
    let mut param_sparse = 0usize;
    let mut logical_entries = 0usize;
    let mut order_counts = BTreeMap::<usize, usize>::new();
    let mut name_counts = BTreeMap::<String, usize>::new();

    for tensor in &net.store.tensors {
        let name = tensor
            .name()
            .map(|name| name.to_string())
            .unwrap_or_else(|| "<anonymous>".to_owned());
        *name_counts.entry(name).or_default() += 1;

        match tensor {
            ParamOrConcrete::Concrete(tensor) => {
                concrete_tensors += 1;
                concrete_entries += tensor.iter_flat().count();
                *order_counts.entry(tensor.structure().order()).or_default() += 1;
                logical_entries += tensor.structure().size().unwrap_or(0);
            }
            ParamOrConcrete::Param(tensor) => {
                param_tensors += 1;
                param_entries += tensor.tensor.actual_size();
                *order_counts
                    .entry(tensor.tensor.structure().order())
                    .or_default() += 1;
                logical_entries += tensor.tensor.structure().size().unwrap_or(0);
                match &tensor.tensor {
                    DataTensor::Dense(_) => param_dense += 1,
                    DataTensor::Sparse(_) => param_sparse += 1,
                }
            }
        }
    }

    let mut top_names = name_counts.into_iter().collect::<Vec<_>>();
    top_names.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    top_names.truncate(8);

    eprintln!(
        "{label} actual_net stats: graph_nodes={} graph_edges={} tensors={} concrete_tensors={} param_tensors={} scalars={} logical_entries={} concrete_entries={} param_entries={} param_dense={} param_sparse={} orders={:?} top_names={:?}",
        net.graph.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len(),
        concrete_tensors,
        param_tensors,
        net.store.scalar.len(),
        logical_entries,
        concrete_entries,
        param_entries,
        param_dense,
        param_sparse,
        order_counts,
        top_names,
    );
}

fn report_concrete_net_stats(label: &str, net: &ConcreteParsingNet) {
    let mut concrete_tensors = 0usize;
    let mut param_tensors = 0usize;
    let mut concrete_entries = 0usize;
    let mut param_entries = 0usize;
    let mut order_counts = BTreeMap::<usize, usize>::new();
    let mut name_counts = BTreeMap::<String, usize>::new();

    for tensor in &net.store.tensors {
        let name = tensor
            .name()
            .map(|name| name.to_string())
            .unwrap_or_else(|| "<anonymous>".to_owned());
        *name_counts.entry(name).or_default() += 1;

        match tensor {
            ParamOrConcrete::Concrete(tensor) => {
                concrete_tensors += 1;
                concrete_entries += tensor.iter_flat().count();
                *order_counts.entry(tensor.structure().order()).or_default() += 1;
            }
            ParamOrConcrete::Param(tensor) => {
                param_tensors += 1;
                param_entries += tensor.tensor.actual_size();
                *order_counts
                    .entry(tensor.tensor.structure().order())
                    .or_default() += 1;
            }
        }
    }

    let mut top_names = name_counts.into_iter().collect::<Vec<_>>();
    top_names.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    top_names.truncate(8);

    eprintln!(
        "{label} concrete_net stats: graph_nodes={} graph_edges={} tensors={} concrete_tensors={} param_tensors={} scalars={} concrete_entries={} param_entries={} orders={:?} top_names={:?}",
        net.graph.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len(),
        concrete_tensors,
        param_tensors,
        net.store.scalar.len(),
        concrete_entries,
        param_entries,
        order_counts,
        top_names,
    );
}

fn execute_actual_net(label: &str, mut net: ParsingNet) -> Atom {
    let lib = TENSORLIB.read().unwrap();
    let start = Instant::now();
    net.execute::<Sequential, SmallestDegree, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| panic!("{label} actual execution failed: {error}"));
    report(&format!("{label} actual_execute"), start);
    report_actual_net_stats(&format!("{label} after_execute"), &net);

    let result = net
        .result_scalar()
        .unwrap_or_else(|error| panic!("{label} actual scalar result failed: {error}"));
    let result = match result {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    };
    eprintln!(
        "{label} actual_result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    result
}

fn execute_actual_net_min_result_rank(label: &str, mut net: ParsingNet) -> Atom {
    let lib = TENSORLIB.read().unwrap();
    let start = Instant::now();
    net.execute::<Sequential, MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| panic!("{label} min-result-rank actual execution failed: {error}"));
    report(&format!("{label} min_result_rank_actual_execute"), start);
    report_actual_net_stats(&format!("{label} min_result_rank_after_execute"), &net);

    let result = net
        .result_scalar()
        .unwrap_or_else(|error| panic!("{label} min-result-rank scalar result failed: {error}"));
    let result = match result {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    };
    eprintln!(
        "{label} min_result_rank_actual_result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    result
}

fn execute_actual_net_min_result_rank_parallel(label: &str, mut net: ParsingNet) -> Atom {
    let lib = TENSORLIB.read().unwrap();
    let start = Instant::now();
    net.execute_parallel::<MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| {
            panic!("{label} parallel min-result-rank actual execution failed: {error}")
        });
    report(
        &format!("{label} parallel_min_result_rank_actual_execute"),
        start,
    );
    report_actual_net_stats(
        &format!("{label} parallel_min_result_rank_after_execute"),
        &net,
    );

    let result = net.result_scalar().unwrap_or_else(|error| {
        panic!("{label} parallel min-result-rank scalar result failed: {error}")
    });
    let result = match result {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    };
    eprintln!(
        "{label} parallel_min_result_rank_actual_result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    result
}

fn execute_actual_net_steps(label: &str, mut net: ParsingNet, steps: usize) {
    let lib = TENSORLIB.read().unwrap();
    for step in 0..steps {
        spenso::network::profile::reset();
        let start = Instant::now();
        net.execute::<Steps<1>, SmallestDegree, _, _, _>(&*lib, &*FUN_LIB)
            .unwrap_or_else(|error| panic!("{label} actual execution step {step} failed: {error}"));
        report(&format!("{label} actual_execute_step_{step}"), start);
        report_actual_net_stats(&format!("{label} after_step_{step}"), &net);
        spenso::network::profile::report(&format!("{label} after_actual_execute_step_{step}"));
    }
}

fn parse_actual_net(label: &str, expr: &Atom) -> ParsingNet {
    spenso::network::profile::reset();
    let settings = ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..Default::default()
    };
    let lib = TENSORLIB.read().unwrap();
    let start = Instant::now();
    let net = ParsingNet::try_from_view(expr.as_view(), &*lib, &settings)
        .unwrap_or_else(|error| panic!("{label} actual parse failed: {error}"));
    report(&format!("{label} actual_network_parse"), start);
    report_actual_net_stats(label, &net);
    spenso::network::profile::report(&format!("{label} after_actual_network_parse"));
    net
}

fn parse_concrete_hep_net(label: &str, expr: &Atom) -> ConcreteParsingNet {
    spenso::network::profile::reset();
    let settings = ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..Default::default()
    };
    let lib = spenso_hep_lib::hep_lib(F(1.0), F(0.0));
    let start = Instant::now();
    let net = ConcreteParsingNet::try_from_view(expr.as_view(), &lib, &settings)
        .unwrap_or_else(|error| panic!("{label} concrete HEP parse failed: {error}"));
    report(&format!("{label} concrete_hep_network_parse"), start);
    report_concrete_net_stats(label, &net);
    spenso::network::profile::report(&format!("{label} after_concrete_hep_network_parse"));
    net
}

fn execute_concrete_hep_net(label: &str, mut net: ConcreteParsingNet) {
    let lib = spenso_hep_lib::hep_lib(F(1.0), F(0.0));
    let start = Instant::now();
    net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &*FUN_LIB)
        .unwrap_or_else(|error| panic!("{label} concrete HEP execution failed: {error}"));
    report(&format!("{label} concrete_hep_execute"), start);
    report_concrete_net_stats(&format!("{label} after_concrete_hep_execute"), &net);
}

fn parse_inline_expression(input: &str) -> Atom {
    Atom::parse(wrap_input!(input), SymbolicaParseSettings::default())
        .expect("inline diagnostic expression should parse")
}

fn gamma_ladder_order_mwe_expr() -> Atom {
    parse_inline_expression(
        r#"
        spenso::g(spenso::bis(4,4),spenso::bis(4,5))
        *spenso::gamma(spenso::bis(4,4),spenso::bis(4,23),spenso::mink(4,0))
        *spenso::g(spenso::bis(4,22),spenso::bis(4,23))
        *ϵ(0,spenso::mink(4,0))
        *spenso::gamma(spenso::bis(4,5),spenso::bis(4,6),spenso::mink(4,1))
        *spenso::g(spenso::bis(4,6),spenso::bis(4,7))
        *ϵ(1,spenso::mink(4,1))
        *spenso::gamma(spenso::bis(4,7),spenso::bis(4,12),spenso::mink(4,14))
        *spenso::g(spenso::bis(4,12),spenso::bis(4,13))
        *spenso::gamma(spenso::bis(4,13),spenso::bis(4,16),spenso::mink(4,18))
        *spenso::g(spenso::bis(4,16),spenso::bis(4,17))
        *spenso::gamma(spenso::bis(4,8),spenso::bis(4,17),spenso::mink(4,2))
        *spenso::g(spenso::bis(4,8),spenso::bis(4,9))
        *ϵbar(2,spenso::mink(4,2))
        *spenso::gamma(spenso::bis(4,9),spenso::bis(4,10),spenso::mink(4,3))
        *spenso::g(spenso::bis(4,10),spenso::bis(4,11))
        *ϵbar(3,spenso::mink(4,3))
        *spenso::g(spenso::bis(4,20),spenso::bis(4,21))
        *spenso::gamma(spenso::bis(4,11),spenso::bis(4,20),spenso::mink(4,18))
        *spenso::gamma(spenso::bis(4,21),spenso::bis(4,22),spenso::mink(4,14))
        "#,
    )
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path"]
fn symbolica_expression_input_actual_network_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw", &expr);
    let _ = execute_actual_net("raw", net);
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path with MinResultRank"]
fn symbolica_expression_input_actual_network_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw_min_result_rank", &expr);
    let _ = execute_actual_net_min_result_rank("raw_min_result_rank", net);
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path with parallel MinResultRank"]
fn symbolica_expression_input_actual_network_parallel_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw_parallel_min_result_rank", &expr);
    let _ = execute_actual_net_min_result_rank_parallel("raw_parallel_min_result_rank", net);
}

#[test]
#[ignore = "diagnostic timing for the first few GammaLoop concrete/mixed tensor execution steps"]
fn symbolica_expression_input_actual_network_first_steps() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw_steps", &expr);
    execute_actual_net_steps("raw_steps", net, 8);
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path after Hornering"]
fn symbolica_expression_input_horner_actual_network_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let expr = collect_horner(&expr);
    let net = parse_actual_net("horner", &expr);
    let _ = execute_actual_net("horner", net);
}

#[test]
#[ignore = "diagnostic timing for concrete HEP tensor data execution"]
fn concrete_hep_gamma_trace_probe() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_inline_expression(
        "spenso::gamma(spenso::bis(4,1),spenso::bis(4,2),spenso::mink(4,3))*spenso::gamma(spenso::bis(4,2),spenso::bis(4,1),spenso::mink(4,3))",
    );
    let net = parse_concrete_hep_net("gamma_trace_probe", &expr);
    execute_concrete_hep_net("gamma_trace_probe", net);
}

#[test]
#[ignore = "diagnostic reproduction of the current left-deep SmallestDegree gamma ladder"]
fn gamma_ladder_mwe_smallest_degree_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = gamma_ladder_order_mwe_expr();
    let net = parse_actual_net("gamma_ladder_mwe_smallest_degree", &expr);
    let result = execute_actual_net("gamma_ladder_mwe_smallest_degree", net);
    assert!(result.nterms() > 0);
}

#[test]
#[ignore = "diagnostic comparison order that minimizes intermediate result rank"]
fn gamma_ladder_mwe_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = gamma_ladder_order_mwe_expr();
    let net = parse_actual_net("gamma_ladder_mwe_min_result_rank", &expr);
    let result = execute_actual_net_min_result_rank("gamma_ladder_mwe_min_result_rank", net);
    assert!(result.nterms() > 0);
}

#[test]
#[ignore = "diagnostic equivalence check for contraction-order output forms"]
fn gamma_ladder_mwe_order_results_match_after_expand() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = gamma_ladder_order_mwe_expr();

    let smallest = execute_actual_net(
        "gamma_ladder_mwe_match_smallest_degree",
        parse_actual_net("gamma_ladder_mwe_match_smallest_degree", &expr),
    );
    let min_result_rank = execute_actual_net_min_result_rank(
        "gamma_ladder_mwe_match_min_result_rank",
        parse_actual_net("gamma_ladder_mwe_match_min_result_rank", &expr),
    );

    let start = Instant::now();
    let diff = (smallest - min_result_rank).expand();
    report("gamma_ladder_mwe_order_diff_expand", start);
    eprintln!(
        "gamma_ladder_mwe_order_diff stats: terms={} bytes={}",
        diff.nterms(),
        diff.as_view().get_byte_size()
    );
    assert!(
        diff.is_zero(),
        "contraction orders produced different scalar values: {}",
        diff
    );
}
