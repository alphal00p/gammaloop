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
        ExecutionResult, Network, Sequential, SmallestDegree, Steps,
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

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path"]
fn symbolica_expression_input_actual_network_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw", &expr);
    let _ = execute_actual_net("raw", net);
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
