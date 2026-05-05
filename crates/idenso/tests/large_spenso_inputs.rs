use std::{
    path::{Path, PathBuf},
    time::Instant,
};

use idenso::{
    parsing_ind::Parsind,
    representations::initialize,
    tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor},
};
use spenso::{
    network::parsing::{ParseSettings, StructureFromAtom},
    structure::{
        TensorStructure,
        representation::LibraryRep,
        slot::{IsAbstractSlot, Slot},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate},
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("idenso crate should live two levels below the workspace root")
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

fn find_first_function_call<'a>(input: &'a str, name: &str) -> Option<&'a str> {
    let start = input.find(name)?;
    let open = start + name.len();
    if input.as_bytes().get(open) != Some(&b'(') {
        return None;
    }

    let mut depth = 0usize;
    for (offset, byte) in input[open..].bytes().enumerate() {
        match byte {
            b'(' => depth += 1,
            b')' => {
                depth -= 1;
                if depth == 0 {
                    return Some(&input[start..=open + offset]);
                }
            }
            _ => {}
        }
    }

    None
}

fn report_first_tensor_parse(input: &str, name: &str) {
    let Some(call) = find_first_function_call(input, name) else {
        eprintln!("sample_tensor_parse {name}: not_found");
        return;
    };

    let parsed =
        Atom::parse(wrap_input!(call), SymbolicaParseSettings::default()).expect("sample parses");
    let symbol = match parsed.as_view() {
        AtomView::Fun(fun) => fun.get_symbol().to_string(),
        other => other.to_string(),
    };
    match SymbolicTensor::<Parsind>::parse(parsed.as_view()) {
        Ok(tensor) => eprintln!(
            "sample_tensor_parse {name}: ok symbol={} order={} expression_bytes={}",
            symbol,
            tensor.structure.structure.order(),
            parsed.as_view().get_byte_size()
        ),
        Err(error) => eprintln!("sample_tensor_parse {name}: err symbol={symbol} error={error}"),
    }
}

fn report_first_slot_parse(input: &str, name: &str) {
    let Some(call) = find_first_function_call(input, name) else {
        eprintln!("sample_slot_parse {name}: not_found");
        return;
    };

    let parsed =
        Atom::parse(wrap_input!(call), SymbolicaParseSettings::default()).expect("sample parses");
    match Slot::<LibraryRep, Parsind>::try_from(parsed.as_view()) {
        Ok(slot) => eprintln!(
            "sample_slot_parse {name}: ok rep={} atom={}",
            slot.rep(),
            slot.to_atom()
        ),
        Err(error) => eprintln!(
            "sample_slot_parse {name}: err expr={} error={error}",
            parsed.as_view().to_plain_string()
        ),
    }
}

fn parse_root_input(filename: &str) -> Atom {
    initialize();

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
    report_first_slot_parse(&input, "spenso::bis");
    report_first_slot_parse(&input, "spenso::mink");
    report_first_tensor_parse(&input, "gammalooprs::Q");
    report_first_tensor_parse(&input, "gammalooprs::OSE");
    report_first_tensor_parse(&input, "spenso::gamma");
    report_first_tensor_parse(&input, "spenso::g");

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

fn parse_smallest_input() -> Atom {
    parse_root_input("symbolica_expression.txt")
}

fn time_symbolic_network_parse_only(label: &str, expr: &Atom, settings: &ParseSettings) {
    spenso::network::profile::reset();
    let start = Instant::now();
    let net = expr
        .parse_to_symbolic_net::<Parsind>(settings)
        .expect("failed to parse expression into symbolic network");
    eprintln!(
        "{label} symbolic_net stats: graph_nodes={} graph_edges={} tensors={} scalars={}",
        net.graph.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len(),
        net.store.scalar.len(),
    );
    report(&format!("{label} symbolic_network_parse"), start);
    spenso::network::profile::report(&format!("{label} after_symbolic_network_parse"));
}

fn time_symbolic_network_parse_and_execute(label: &str, expr: &Atom, settings: &ParseSettings) {
    spenso::network::profile::reset();
    let start = Instant::now();
    let net = expr
        .parse_to_symbolic_net::<Parsind>(settings)
        .expect("failed to parse expression into symbolic network");
    eprintln!(
        "{label} symbolic_net stats: graph_nodes={} graph_edges={} tensors={} scalars={}",
        net.graph.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len(),
        net.store.scalar.len(),
    );
    report(&format!("{label} symbolic_network_parse"), start);
    spenso::network::profile::report(&format!("{label} after_symbolic_network_parse"));

    let start = Instant::now();
    let result = net.simple_execute::<()>();
    eprintln!(
        "{label} symbolic_execute result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    report(&format!("{label} symbolic_execute"), start);
    spenso::network::profile::report(&format!("{label} after_symbolic_execute"));
}

fn collect_common_factors(expr: &Atom) -> Atom {
    let start = Instant::now();
    let collected = expr.collect_factors();
    eprintln!(
        "collect_factors stats: terms={} bytes={}",
        collected.nterms(),
        collected.as_view().get_byte_size()
    );
    report("collect_factors", start);
    collected
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

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn smallest_root_input_symbolic_network_parse_default() {
    let expr = parse_smallest_input();
    let settings = ParseSettings::default();
    time_symbolic_network_parse_and_execute("default", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn smallest_root_input_symbolic_network_parse_raw_full_depth_profile() {
    let expr = parse_smallest_input();
    let settings = ParseSettings::default();
    time_symbolic_network_parse_only("raw_full_depth", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn smallest_root_input_symbolic_network_parse_without_scalar_precontraction() {
    let expr = parse_smallest_input();
    let settings = ParseSettings {
        precontract_scalars: false,
        ..Default::default()
    };
    time_symbolic_network_parse_and_execute("no_precontract_scalars", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn smallest_root_input_symbolic_network_parse_top_level_tensor_factors() {
    let expr = parse_smallest_input();
    let settings = ParseSettings {
        depth_limit: Some(2),
        ..Default::default()
    };
    time_symbolic_network_parse_and_execute("top_level_tensor_factors", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn smallest_root_input_symbolic_network_parse_after_collect_factors() {
    let expr = parse_smallest_input();
    let expr = collect_common_factors(&expr);
    let settings = ParseSettings {
        depth_limit: Some(2),
        ..Default::default()
    };
    time_symbolic_network_parse_and_execute("after_collect_factors", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn smallest_root_input_symbolic_network_parse_after_collect_horner() {
    let expr = parse_smallest_input();
    let expr = collect_horner(&expr);
    let settings = ParseSettings {
        depth_limit: Some(2),
        ..Default::default()
    };
    time_symbolic_network_parse_and_execute("after_collect_horner", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn larger_root_input_symbolic_network_parse_after_collect_horner() {
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let expr = collect_horner(&expr);
    let settings = ParseSettings {
        depth_limit: Some(2),
        ..Default::default()
    };
    time_symbolic_network_parse_only("larger_after_collect_horner", &expr, &settings);
}

#[test]
#[ignore = "diagnostic timing for root-level large Spenso input"]
fn larger_root_input_symbolic_network_parse_raw_full_depth_profile() {
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let settings = ParseSettings::default();
    time_symbolic_network_parse_only("larger_raw_full_depth", &expr, &settings);
}
