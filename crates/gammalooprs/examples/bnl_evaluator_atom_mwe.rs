use std::{
    collections::BTreeSet,
    env, fs,
    io::{self, Write},
    ops::Deref,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use color_eyre::{
    Result,
    eyre::{WrapErr, bail, eyre},
};
use gammalooprs::{
    initialisation::test_initialise,
    numerator::{ParsingNet, aind::Aind, symbolica_ext::AtomCoreExt},
    processes::TensorNetworkContractionOrder,
    utils::{FUN_LIB, TENSORLIB, symbolica_ext::LogPrint},
};
use idenso::tensor::{SymbolicNet, SymbolicNetExt, SymbolicTensor};
use spenso::network::{
    DEFAULT_EXACT_JOIN_LIMIT, ExecutionResult, MinResultRank, MinResultRankWith, Network,
    PAIR_SCORE_ATOM_AWARE, PAIR_SCORE_ENTRY_AWARE, PAIR_SCORE_RESULT_RANK_ONLY, ScalarAliases,
    Sequential, SmallestDegree,
    library::{DummyLibrary, function_lib::Wrap},
    parsing::{
        ParseSettings as NetworkParseSettings, SchoonschipExpansionMode, ShorthandParsing,
        StrictTensorFilter, StructureInferenceMode,
    },
    store::NetworkStore,
};
use spenso::shadowing::TensorCollectExt;
use spenso::tensors::{data::HasTensorData, parametric::ParamOrConcrete};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    id::AliasedAtom,
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};
use tabled::{builder::Builder, settings::Style};

macro_rules! execute_concrete_network {
    ($net:expr, $order:expr) => {
        match $order {
            TensorNetworkContractionOrder::SparseAtomAware => $net.execute::<
                Sequential,
                MinResultRank,
                _,
                _,
                _,
            >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()),
            TensorNetworkContractionOrder::AtomAware => $net.execute::<
                Sequential,
                MinResultRankWith<{ PAIR_SCORE_ATOM_AWARE }, { DEFAULT_EXACT_JOIN_LIMIT }>,
                _,
                _,
                _,
            >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()),
            TensorNetworkContractionOrder::ResultRankOnly => $net.execute::<
                Sequential,
                MinResultRankWith<{ PAIR_SCORE_RESULT_RANK_ONLY }, { DEFAULT_EXACT_JOIN_LIMIT }>,
                _,
                _,
                _,
            >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()),
            TensorNetworkContractionOrder::EntryAware => $net.execute::<
                Sequential,
                MinResultRankWith<{ PAIR_SCORE_ENTRY_AWARE }, { DEFAULT_EXACT_JOIN_LIMIT }>,
                _,
                _,
                _,
            >(TENSORLIB.read().unwrap().deref(), FUN_LIB.deref()),
        }
    };
}

#[derive(Debug, Clone, Copy)]
enum TermSelection {
    All,
    First(usize),
    Term(usize),
    Range { start: usize, end: usize },
}

impl TermSelection {
    fn parse(value: &str) -> Result<Self> {
        if value == "all" {
            return Ok(Self::All);
        }
        if let Some(count) = value.strip_prefix("first:") {
            return Ok(Self::First(parse_index("first count", count)?));
        }
        if let Some(index) = value.strip_prefix("term:") {
            return Ok(Self::Term(parse_index("term index", index)?));
        }
        if let Some(range) = value.strip_prefix("range:") {
            let Some((start, end)) = range.split_once("..") else {
                bail!("range selectors use range:START..END with END exclusive");
            };
            return Ok(Self::Range {
                start: parse_index("range start", start)?,
                end: parse_index("range end", end)?,
            });
        }

        bail!("unknown selector {value:?}; use all, first:N, term:N, or range:START..END")
    }

    fn selected_label(self) -> String {
        match self {
            Self::All => "all".to_owned(),
            Self::First(count) => format!("first:{count}"),
            Self::Term(index) => format!("term:{index}"),
            Self::Range { start, end } => format!("range:{start}..{end}"),
        }
    }
}

struct Config {
    input_path: PathBuf,
    selection: TermSelection,
    execute: bool,
    dump_network_path: Option<PathBuf>,
    print_log: bool,
    unit_scalars: bool,
    kept_scalars: Option<BTreeSet<usize>>,
    scalar_details: bool,
    alias_scalar_threshold: Option<usize>,
    resolve_aliases: bool,
    dump_aliased_path: Option<PathBuf>,
    dump_scalar_aliases_md_path: Option<PathBuf>,
    symbolic_net: bool,
    symbolic_then_concrete: bool,
    collect_tensors_before_concrete: bool,
    symbolic_expand_shorthands: bool,
    symbolic_strict_filter: StrictTensorFilter,
    contraction_order: TensorNetworkContractionOrder,
}

fn main() -> Result<()> {
    test_initialise()?;

    let config = Config::from_args()?;
    let input = fs::read_to_string(&config.input_path)
        .with_context(|| format!("failed to read {}", config.input_path.display()))?;
    let mut summary = SummaryTable::default();

    print_metric_table(
        "input",
        [
            ("path", config.input_path.display().to_string()),
            ("bytes", input.len().to_string()),
        ],
    );
    summary.push(
        "input",
        "read",
        "",
        "",
        input.len().to_string(),
        format!("path={}", config.input_path.display()),
    );
    print_metric_table(
        "contraction_order",
        [(
            "order",
            contraction_order_label(config.contraction_order).to_owned(),
        )],
    );
    summary.push(
        "contraction_order",
        contraction_order_label(config.contraction_order),
        "",
        "",
        "",
        "",
    );

    let parse_started = Instant::now();
    let atom =
        Atom::parse_with_default_namespace(wrap_input!(&input), SymbolicaParseSettings::default())
            .map_err(|error| eyre!("failed to parse input atom: {error}"))?;
    let atom_stats = print_atom_stats("symbolica_parse", parse_started, &atom);
    summary.push(
        "symbolica_parse",
        "done",
        format_duration(atom_stats.elapsed),
        atom_stats.terms.to_string(),
        atom_stats.bytes.to_string(),
        "",
    );
    if config.print_log {
        println!("input_log_print_120\n{}", atom.log_print(Some(120)));
    }

    let selection_started = Instant::now();
    let selected = select_terms(&atom, config.selection)?;
    let selected_stats = AtomStats::from_atom(selection_started.elapsed(), &selected);
    if !matches!(config.selection, TermSelection::All) {
        print_metric_table(
            "selection",
            [
                ("selector", config.selection.selected_label()),
                ("elapsed", format_duration(selected_stats.elapsed)),
            ],
        );
        print_atom_shape("selected_atom", &selected);
        if config.print_log {
            println!("selected_log_print_120\n{}", selected.log_print(Some(120)));
        }
    }
    summary.push(
        "selection",
        config.selection.selected_label(),
        format_duration(selected_stats.elapsed),
        selected_stats.terms.to_string(),
        selected_stats.bytes.to_string(),
        "",
    );

    if config.symbolic_then_concrete {
        run_symbolic_then_concrete(&config, &selected, &mut summary)?;
        summary.print();
        return Ok(());
    }

    if config.symbolic_net {
        run_symbolic_net(&config, &selected, &mut summary)?;
        summary.print();
        return Ok(());
    }

    print_metric_table("network_parse", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let network_parse_started = Instant::now();
    let mut net = selected
        .parse_into_net()
        .wrap_err("failed to parse evaluator atom into tensor network")?;
    let network_parse_elapsed = network_parse_started.elapsed();
    print_metric_table(
        "network_parse",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(network_parse_elapsed)),
            ("nodes", net.graph.n_nodes().to_string()),
            ("edges", net.graph.graph.n_edges().to_string()),
            ("tensors", net.store.tensors.len().to_string()),
        ],
    );
    summary.push(
        "network_parse",
        "done",
        format_duration(network_parse_elapsed),
        "",
        "",
        format!(
            "nodes={} edges={} tensors={}",
            net.graph.n_nodes(),
            net.graph.graph.n_edges(),
            net.store.tensors.len()
        ),
    );
    let scalar_store_stats = print_scalar_store_stats(&net, config.scalar_details);
    summary.push(
        "scalar_store",
        "done",
        "",
        scalar_store_stats.total_terms.to_string(),
        scalar_store_stats.total_bytes.to_string(),
        format!(
            "count={} top_level_sums={} max_terms={} max_bytes={}",
            scalar_store_stats.count,
            scalar_store_stats.top_level_sums,
            scalar_store_stats.max_terms,
            scalar_store_stats.max_bytes
        ),
    );

    if config.unit_scalars || config.kept_scalars.is_some() {
        let replacement_stats =
            replace_scalars(&mut net, config.unit_scalars, config.kept_scalars.as_ref())?;
        summary.push(
            "scalar_replacement",
            replacement_stats.mode,
            "",
            "",
            "",
            format!(
                "replaced={} kept={}",
                replacement_stats.replaced, replacement_stats.kept
            ),
        );
    }

    let tensor_entry_stats =
        print_tensor_entry_scalar_stats("tensor_entry_scalars_before_alias", &net);
    summary.push(
        "tensor_entry_scalars",
        "before_alias",
        "",
        tensor_entry_stats.total_terms.to_string(),
        tensor_entry_stats.total_bytes.to_string(),
        tensor_entry_stats.summary_details(),
    );

    let scalar_aliases = if let Some(threshold) = config.alias_scalar_threshold {
        let aliases =
            net.alias_scalar_refs(|_, scalar| scalar.as_view().get_byte_size() >= threshold);
        let tensor_entry_stats =
            print_tensor_entry_scalar_stats("tensor_entry_scalars_after_alias", &net);
        summary.push(
            "tensor_entry_scalars",
            "after_alias",
            "",
            tensor_entry_stats.total_terms.to_string(),
            tensor_entry_stats.total_bytes.to_string(),
            tensor_entry_stats.summary_details(),
        );
        print_metric_table(
            "scalar_aliases",
            [
                ("threshold_bytes", threshold.to_string()),
                ("created", aliases.aliases_created().to_string()),
                ("terms", aliases.aliased_terms().to_string()),
                ("bytes", aliases.aliased_bytes().to_string()),
                ("max_bytes", aliases.max_aliased_bytes().to_string()),
            ],
        );
        summary.push(
            "scalar_aliases",
            "enabled",
            "",
            aliases.aliased_terms().to_string(),
            aliases.aliased_bytes().to_string(),
            format!(
                "threshold_bytes={threshold} created={} max_bytes={}",
                aliases.aliases_created(),
                aliases.max_aliased_bytes()
            ),
        );
        Some(aliases)
    } else {
        summary.push("scalar_aliases", "disabled", "", "", "", "");
        None
    };

    if let Some(path) = &config.dump_scalar_aliases_md_path {
        let Some(aliases) = &scalar_aliases else {
            bail!("--dump-scalar-aliases-md requires --alias-scalars");
        };
        let bytes = write_scalar_aliases_markdown(
            path,
            &config.input_path,
            config.alias_scalar_threshold.unwrap_or_default(),
            aliases,
            &net,
        )?;
        print_metric_table(
            "scalar_aliases_md",
            [
                ("path", path.display().to_string()),
                ("bytes", bytes.to_string()),
            ],
        );
        summary.push(
            "scalar_aliases_md",
            "written",
            "",
            "",
            bytes.to_string(),
            path.display().to_string(),
        );
    }

    if let Some(path) = config.dump_network_path {
        fs::write(&path, net.dot_pretty())
            .with_context(|| format!("failed to write {}", path.display()))?;
        print_metric_table("network_dump", [("path", path.display().to_string())]);
        summary.push(
            "network_dump",
            "written",
            "",
            "",
            "",
            path.display().to_string(),
        );
    }

    if !config.execute {
        print_metric_table("network_execute", [("skipped", "true".to_owned())]);
        summary.push("network_execute", "skipped", "", "", "", "");
        summary.print();
        return Ok(());
    }

    print_metric_table("network_execute", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let execute_started = Instant::now();
    execute_concrete_network!(net, config.contraction_order)
        .wrap_err("failed to execute evaluator atom tensor network")?;
    let execute_elapsed = execute_started.elapsed();
    print_metric_table(
        "network_execute",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(execute_elapsed)),
        ],
    );
    summary.push(
        "network_execute",
        "done",
        format_duration(execute_elapsed),
        "",
        "",
        "",
    );
    let tensor_entry_stats =
        print_tensor_entry_scalar_stats("tensor_entry_scalars_after_execute", &net);
    summary.push(
        "tensor_entry_scalars",
        "after_execute",
        "",
        tensor_entry_stats.total_terms.to_string(),
        tensor_entry_stats.total_bytes.to_string(),
        tensor_entry_stats.summary_details(),
    );

    let result_started = Instant::now();
    let result = net.result_scalar();
    let result_elapsed = result_started.elapsed();
    let result_status = if result.is_ok() { "ok" } else { "err" };
    print_metric_table(
        "result_scalar",
        [
            ("status", result_status.to_owned()),
            ("elapsed", format_duration(result_elapsed)),
        ],
    );
    summary.push(
        "result_scalar",
        result_status,
        format_duration(result_elapsed),
        "",
        "",
        "",
    );
    if let Ok(result) = result {
        let root = match result {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        };
        if let Some(aliases) = &scalar_aliases {
            let aliased = net.aliased_atom(aliases, root);
            if let Some(path) = &config.dump_aliased_path {
                let bytes = write_aliased_atom_dump(path, &aliased)?;
                print_metric_table(
                    "aliased_dump",
                    [
                        ("path", path.display().to_string()),
                        ("bytes", bytes.to_string()),
                    ],
                );
                summary.push(
                    "aliased_dump",
                    "written",
                    "",
                    "",
                    bytes.to_string(),
                    path.display().to_string(),
                );
            }
            print_metric_table(
                "result_aliased_root",
                [
                    ("aliases", aliased.get_aliases().len().to_string()),
                    ("terms", aliased.get_root().nterms().to_string()),
                    (
                        "bytes",
                        aliased.get_root().as_view().get_byte_size().to_string(),
                    ),
                ],
            );
            summary.push(
                "result_aliased_root",
                "done",
                "",
                aliased.get_root().nterms().to_string(),
                aliased.get_root().as_view().get_byte_size().to_string(),
                format!("aliases={}", aliased.get_aliases().len()),
            );
            if config.resolve_aliases {
                let resolve_started = Instant::now();
                let resolved = aliased.into_inner();
                let resolve_elapsed = resolve_started.elapsed();
                print_metric_table(
                    "result_alias_resolve",
                    [
                        ("elapsed", format_duration(resolve_elapsed)),
                        ("terms", resolved.nterms().to_string()),
                        ("bytes", resolved.as_view().get_byte_size().to_string()),
                    ],
                );
                summary.push(
                    "result_alias_resolve",
                    "done",
                    format_duration(resolve_elapsed),
                    resolved.nterms().to_string(),
                    resolved.as_view().get_byte_size().to_string(),
                    "",
                );
            } else {
                print_metric_table("result_alias_resolve", [("skipped", "true".to_owned())]);
                summary.push("result_alias_resolve", "skipped", "", "", "", "");
            }
        } else {
            if let Some(path) = &config.dump_aliased_path {
                let bytes = write_aliased_atom_dump(path, &AliasedAtom::from(root.clone()))?;
                print_metric_table(
                    "aliased_dump",
                    [
                        ("path", path.display().to_string()),
                        ("bytes", bytes.to_string()),
                    ],
                );
                summary.push(
                    "aliased_dump",
                    "written",
                    "",
                    "",
                    bytes.to_string(),
                    path.display().to_string(),
                );
            }
            print_metric_table(
                "result_atom",
                [
                    ("terms", root.nterms().to_string()),
                    ("bytes", root.as_view().get_byte_size().to_string()),
                ],
            );
            summary.push(
                "result_atom",
                "done",
                "",
                root.nterms().to_string(),
                root.as_view().get_byte_size().to_string(),
                "",
            );
        }
    }

    summary.print();
    Ok(())
}

impl Config {
    fn from_args() -> Result<Self> {
        let mut input_path =
            PathBuf::from("examples/cli/BNL/profiling/bnl_integrated_evaluator_atom_exact.sym");
        let mut selection = TermSelection::All;
        let mut execute = true;
        let mut dump_network_path = None;
        let mut print_log = false;
        let mut unit_scalars = false;
        let mut kept_scalars = None;
        let mut scalar_details = false;
        let mut alias_scalar_threshold = None;
        let mut resolve_aliases = false;
        let mut dump_aliased_path = None;
        let mut dump_scalar_aliases_md_path = None;
        let mut symbolic_net = false;
        let mut symbolic_then_concrete = false;
        let mut collect_tensors_before_concrete = false;
        let mut symbolic_expand_shorthands = false;
        let mut symbolic_strict_filter = StrictTensorFilter::ContainsReps;
        let mut contraction_order = TensorNetworkContractionOrder::default();
        let mut input_path_seen = false;
        let mut selection_seen = false;

        let mut args = env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "-h" | "--help" => {
                    print_usage();
                    std::process::exit(0);
                }
                "--skip-execute" => execute = false,
                "--print-log" => print_log = true,
                "--unit-scalars" => unit_scalars = true,
                "--scalar-details" => scalar_details = true,
                "--resolve-aliases" => resolve_aliases = true,
                "--symbolic-net" => symbolic_net = true,
                "--symbolic-then-concrete" => symbolic_then_concrete = true,
                "--collect-tensors-before-concrete" => collect_tensors_before_concrete = true,
                "--symbolic-expand-shorthands" => symbolic_expand_shorthands = true,
                "--symbolic-strict-filter" => {
                    let Some(filter) = args.next() else {
                        bail!(
                            "--symbolic-strict-filter requires tagged, tagged-checked, or contains-reps"
                        );
                    };
                    symbolic_strict_filter = parse_strict_tensor_filter(&filter)?;
                }
                "--contraction-order" => {
                    let Some(order) = args.next() else {
                        bail!(
                            "--contraction-order requires sparse-atom-aware, atom-aware, result-rank-only, or entry-aware"
                        );
                    };
                    contraction_order = parse_contraction_order(&order)?;
                }
                "--dump-aliased" => {
                    let Some(path) = args.next() else {
                        bail!("--dump-aliased requires a path");
                    };
                    dump_aliased_path = Some(PathBuf::from(path));
                }
                "--dump-scalar-aliases-md" => {
                    let Some(path) = args.next() else {
                        bail!("--dump-scalar-aliases-md requires a path");
                    };
                    dump_scalar_aliases_md_path = Some(PathBuf::from(path));
                }
                "--alias-scalars" => {
                    let Some(threshold) = args.next() else {
                        bail!("--alias-scalars requires a byte threshold");
                    };
                    alias_scalar_threshold =
                        Some(parse_index("alias scalar threshold", &threshold)?);
                }
                "--keep-scalars" => {
                    let Some(indices) = args.next() else {
                        bail!("--keep-scalars requires a comma-separated index list");
                    };
                    kept_scalars = Some(parse_index_set(&indices)?);
                }
                "--dump-network" => {
                    let Some(path) = args.next() else {
                        bail!("--dump-network requires a path");
                    };
                    dump_network_path = Some(PathBuf::from(path));
                }
                _ if looks_like_selector(&arg) => {
                    if selection_seen {
                        bail!("only one term selector can be supplied");
                    }
                    selection = TermSelection::parse(&arg)?;
                    selection_seen = true;
                }
                _ if !input_path_seen => {
                    input_path = PathBuf::from(arg);
                    input_path_seen = true;
                }
                _ => bail!("unexpected argument {arg:?}"),
            }
        }

        if unit_scalars && kept_scalars.is_some() {
            bail!("--unit-scalars and --keep-scalars are mutually exclusive");
        }
        if symbolic_net && symbolic_then_concrete {
            bail!("--symbolic-net and --symbolic-then-concrete are mutually exclusive");
        }

        Ok(Self {
            input_path,
            selection,
            execute,
            dump_network_path,
            print_log,
            unit_scalars,
            kept_scalars,
            scalar_details,
            alias_scalar_threshold,
            resolve_aliases,
            dump_aliased_path,
            dump_scalar_aliases_md_path,
            symbolic_net,
            symbolic_then_concrete,
            collect_tensors_before_concrete,
            symbolic_expand_shorthands,
            symbolic_strict_filter,
            contraction_order,
        })
    }
}

fn looks_like_selector(value: &str) -> bool {
    value == "all"
        || value.starts_with("first:")
        || value.starts_with("term:")
        || value.starts_with("range:")
}

fn parse_contraction_order(value: &str) -> Result<TensorNetworkContractionOrder> {
    match value {
        "sparse-atom-aware" => Ok(TensorNetworkContractionOrder::SparseAtomAware),
        "atom-aware" => Ok(TensorNetworkContractionOrder::AtomAware),
        "result-rank-only" => Ok(TensorNetworkContractionOrder::ResultRankOnly),
        "entry-aware" => Ok(TensorNetworkContractionOrder::EntryAware),
        _ => bail!("unknown contraction order {value:?}"),
    }
}

fn contraction_order_label(order: TensorNetworkContractionOrder) -> &'static str {
    match order {
        TensorNetworkContractionOrder::SparseAtomAware => "sparse-atom-aware",
        TensorNetworkContractionOrder::AtomAware => "atom-aware",
        TensorNetworkContractionOrder::ResultRankOnly => "result-rank-only",
        TensorNetworkContractionOrder::EntryAware => "entry-aware",
    }
}

fn print_usage() {
    println!(
        "usage: cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- [INPUT] [all|first:N|term:N|range:START..END] [--skip-execute] [--print-log] [--scalar-details] [--unit-scalars|--keep-scalars INDICES] [--alias-scalars BYTES] [--contraction-order ORDER] [--dump-network PATH] [--dump-aliased PATH] [--symbolic-net]"
    );
    println!("term indices are zero-based; range END is exclusive");
    println!(
        "--keep-scalars keeps only the listed scalar-store indices and replaces the rest with 1"
    );
    println!("--alias-scalars replaces scalar refs at or above BYTES with scalar(INDEX) aliases");
    println!("--resolve-aliases expands aliased results back into a full Symbolica atom");
    println!("--dump-aliased writes the post-execution aliased atom JSON dump");
    println!("--dump-scalar-aliases-md writes aliased scalar-store log_print markdown");
    println!(
        "--contraction-order selects sparse-atom-aware, atom-aware, result-rank-only, or entry-aware"
    );
    println!("--symbolic-net executes a SymbolicTensor network and aliases scalar-store atoms");
    println!(
        "--symbolic-then-concrete executes symbolically with aliases, keeps alias placeholders, then parses and executes a concrete tensor network"
    );
    println!(
        "--collect-tensors-before-concrete applies collect_tensors() to the aliased symbolic root before the concrete parse"
    );
    println!(
        "--symbolic-expand-shorthands uses the same chain/trace shorthand expansion mode as parse_into_net for symbolic parsing"
    );
    println!(
        "--symbolic-strict-filter selects symbolic tensor-head detection: tagged, tagged-checked, or contains-reps"
    );
}

fn parse_index(label: &str, value: &str) -> Result<usize> {
    value
        .parse()
        .with_context(|| format!("invalid {label}: {value:?}"))
}

fn parse_index_set(value: &str) -> Result<BTreeSet<usize>> {
    let mut indices = BTreeSet::new();
    for part in value.split(',') {
        if part.is_empty() {
            bail!("empty scalar index in {value:?}");
        }
        indices.insert(parse_index("scalar index", part)?);
    }
    if indices.is_empty() {
        bail!("at least one scalar index is required");
    }
    Ok(indices)
}

fn parse_strict_tensor_filter(value: &str) -> Result<StrictTensorFilter> {
    match value {
        "tagged" => Ok(StrictTensorFilter::Tagged),
        "tagged-checked" => Ok(StrictTensorFilter::TaggedChecked),
        "contains-reps" => Ok(StrictTensorFilter::ContainsReps),
        _ => bail!(
            "invalid symbolic strict tensor filter {value:?}; use tagged, tagged-checked, or contains-reps"
        ),
    }
}

fn format_duration(duration: Duration) -> String {
    format!("{duration:?}")
}

#[derive(Debug, Clone, Copy)]
struct AtomStats {
    elapsed: Duration,
    terms: usize,
    bytes: usize,
}

impl AtomStats {
    fn from_atom(elapsed: Duration, atom: &Atom) -> Self {
        Self {
            elapsed,
            terms: atom.nterms(),
            bytes: atom.as_view().get_byte_size(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct ScalarStoreStats {
    count: usize,
    top_level_sums: usize,
    total_terms: usize,
    max_terms: usize,
    total_bytes: usize,
    max_bytes: usize,
}

#[derive(Debug, Clone, Copy)]
struct TensorEntryStats {
    tensor_count: usize,
    param_tensor_count: usize,
    concrete_tensor_count: usize,
    entry_count: usize,
    total_terms: usize,
    max_terms: usize,
    total_bytes: usize,
    max_bytes: usize,
    max_tensor_index: usize,
    max_entry_index: usize,
}

impl TensorEntryStats {
    fn summary_details(self) -> String {
        format!(
            "tensors={} param_tensors={} concrete_tensors={} entries={} max_terms={} max_bytes={} max_at=tensor:{} entry:{}",
            self.tensor_count,
            self.param_tensor_count,
            self.concrete_tensor_count,
            self.entry_count,
            self.max_terms,
            self.max_bytes,
            self.max_tensor_index,
            self.max_entry_index
        )
    }
}

#[derive(Debug, Clone, Copy)]
struct ScalarReplacementStats {
    mode: &'static str,
    replaced: usize,
    kept: usize,
}

#[derive(Default)]
struct SummaryTable {
    rows: Vec<Vec<String>>,
}

impl SummaryTable {
    fn push(
        &mut self,
        stage: impl Into<String>,
        status: impl Into<String>,
        elapsed: impl Into<String>,
        terms: impl Into<String>,
        bytes: impl Into<String>,
        details: impl Into<String>,
    ) {
        self.rows.push(vec![
            stage.into(),
            status.into(),
            blank_as_dash(elapsed.into()),
            blank_as_dash(terms.into()),
            blank_as_dash(bytes.into()),
            blank_as_dash(details.into()),
        ]);
    }

    fn print(&self) {
        if self.rows.is_empty() {
            return;
        }
        println!("summary");
        print_table(
            ["stage", "status", "elapsed", "terms", "bytes", "details"],
            self.rows.clone(),
        );
    }
}

fn blank_as_dash(value: String) -> String {
    if value.is_empty() {
        "-".to_owned()
    } else {
        value
    }
}

fn print_metric_table(stage: &str, metrics: impl IntoIterator<Item = (&'static str, String)>) {
    let mut table = Builder::new();
    table.push_record(["metric", "value"]);
    for (metric, value) in metrics {
        table.push_record([metric.to_owned(), value]);
    }
    let mut table = table.build();
    table.with(Style::rounded());
    println!("{stage}");
    println!("{table}");
}

fn print_table(
    headers: impl IntoIterator<Item = impl Into<String>>,
    rows: impl IntoIterator<Item = Vec<String>>,
) {
    let mut table = Builder::new();
    table.push_record(headers.into_iter().map(Into::into));
    for row in rows {
        table.push_record(row);
    }
    let mut table = table.build();
    table.with(Style::rounded());
    println!("{table}");
}

fn print_atom_stats(label: &str, started: Instant, atom: &Atom) -> AtomStats {
    let stats = AtomStats::from_atom(started.elapsed(), atom);
    print_metric_table(
        label,
        [
            ("elapsed", format_duration(stats.elapsed)),
            ("terms", stats.terms.to_string()),
            ("bytes", stats.bytes.to_string()),
        ],
    );
    stats
}

fn print_atom_shape(label: &str, atom: &Atom) {
    print_metric_table(
        label,
        [
            ("terms", atom.nterms().to_string()),
            ("bytes", atom.as_view().get_byte_size().to_string()),
        ],
    );
}

fn parse_symbolic_network(config: &Config, selected: &Atom) -> Result<SymbolicNet<Aind>> {
    let shorthand_parsing = if config.symbolic_expand_shorthands {
        ShorthandParsing::Expand {
            schoonschip: SchoonschipExpansionMode::none(),
            trace: true,
            chain: true,
        }
    } else {
        ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        }
    };
    let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

    SymbolicNet::<Aind>::try_from_view::<SymbolicTensor<Aind>, _>(
        selected.as_view(),
        &lib,
        &NetworkParseSettings {
            shorthand_parsing,
            strict_tensor_filter: config.symbolic_strict_filter,
            ..Default::default()
        },
    )
    .wrap_err("failed to parse evaluator atom into symbolic tensor network")
}

fn run_symbolic_then_concrete(
    config: &Config,
    selected: &Atom,
    summary: &mut SummaryTable,
) -> Result<()> {
    print_metric_table("symbolic_network_parse", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let symbolic_parse_started = Instant::now();
    let mut symbolic_net = parse_symbolic_network(config, selected)?;
    let symbolic_parse_elapsed = symbolic_parse_started.elapsed();
    print_metric_table(
        "symbolic_network_parse",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(symbolic_parse_elapsed)),
            ("nodes", symbolic_net.graph.n_nodes().to_string()),
            ("edges", symbolic_net.graph.graph.n_edges().to_string()),
            ("tensors", symbolic_net.store.tensors.len().to_string()),
            (
                "strict_filter",
                format!("{:?}", config.symbolic_strict_filter),
            ),
        ],
    );
    summary.push(
        "symbolic_network_parse",
        "done",
        format_duration(symbolic_parse_elapsed),
        "",
        "",
        format!(
            "nodes={} edges={} tensors={} strict_filter={:?}",
            symbolic_net.graph.n_nodes(),
            symbolic_net.graph.graph.n_edges(),
            symbolic_net.store.tensors.len(),
            config.symbolic_strict_filter
        ),
    );

    let symbolic_scalar_store_stats =
        print_scalar_store_stats(&symbolic_net, config.scalar_details);
    summary.push(
        "symbolic_scalar_store",
        "done",
        "",
        symbolic_scalar_store_stats.total_terms.to_string(),
        symbolic_scalar_store_stats.total_bytes.to_string(),
        format!(
            "count={} top_level_sums={} max_terms={} max_bytes={}",
            symbolic_scalar_store_stats.count,
            symbolic_scalar_store_stats.top_level_sums,
            symbolic_scalar_store_stats.max_terms,
            symbolic_scalar_store_stats.max_bytes
        ),
    );

    if config.unit_scalars || config.kept_scalars.is_some() {
        let replacement_stats = replace_scalars(
            &mut symbolic_net,
            config.unit_scalars,
            config.kept_scalars.as_ref(),
        )?;
        summary.push(
            "symbolic_scalar_replacement",
            replacement_stats.mode,
            "",
            "",
            "",
            format!(
                "replaced={} kept={}",
                replacement_stats.replaced, replacement_stats.kept
            ),
        );
    }

    let symbolic_aliases = if let Some(threshold) = config.alias_scalar_threshold {
        let aliases = symbolic_net
            .alias_scalar_refs(|_, scalar| scalar.as_view().get_byte_size() >= threshold);
        print_metric_table(
            "symbolic_scalar_aliases",
            [
                ("threshold_bytes", threshold.to_string()),
                ("created", aliases.aliases_created().to_string()),
                ("terms", aliases.aliased_terms().to_string()),
                ("bytes", aliases.aliased_bytes().to_string()),
                ("max_bytes", aliases.max_aliased_bytes().to_string()),
            ],
        );
        summary.push(
            "symbolic_scalar_aliases",
            "enabled",
            "",
            aliases.aliased_terms().to_string(),
            aliases.aliased_bytes().to_string(),
            format!(
                "threshold_bytes={threshold} created={} max_bytes={}",
                aliases.aliases_created(),
                aliases.max_aliased_bytes()
            ),
        );
        Some(aliases)
    } else {
        summary.push("symbolic_scalar_aliases", "disabled", "", "", "", "");
        None
    };

    if let Some(path) = &config.dump_scalar_aliases_md_path {
        let Some(aliases) = &symbolic_aliases else {
            bail!("--dump-scalar-aliases-md requires --alias-scalars");
        };
        let bytes = write_scalar_aliases_markdown(
            path,
            &config.input_path,
            config.alias_scalar_threshold.unwrap_or_default(),
            aliases,
            &symbolic_net,
        )?;
        print_metric_table(
            "symbolic_scalar_aliases_md",
            [
                ("path", path.display().to_string()),
                ("bytes", bytes.to_string()),
            ],
        );
        summary.push(
            "symbolic_scalar_aliases_md",
            "written",
            "",
            "",
            bytes.to_string(),
            path.display().to_string(),
        );
    }

    if !config.execute {
        print_metric_table("symbolic_network_execute", [("skipped", "true".to_owned())]);
        summary.push("symbolic_network_execute", "skipped", "", "", "", "");
        return Ok(());
    }

    print_metric_table("symbolic_network_execute", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let symbolic_execute_started = Instant::now();
    let symbolic_lib = DummyLibrary::new();
    symbolic_net
        .execute::<Sequential, SmallestDegree, _, _, _>(&symbolic_lib, &Wrap {})
        .wrap_err("failed to execute symbolic tensor network")?;
    let symbolic_execute_elapsed = symbolic_execute_started.elapsed();
    print_metric_table(
        "symbolic_network_execute",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(symbolic_execute_elapsed)),
        ],
    );
    summary.push(
        "symbolic_network_execute",
        "done",
        format_duration(symbolic_execute_elapsed),
        "",
        "",
        "",
    );

    let symbolic_result_started = Instant::now();
    let symbolic_result = symbolic_net.result_tensor(&symbolic_lib);
    let symbolic_result_elapsed = symbolic_result_started.elapsed();
    let symbolic_result_status = if symbolic_result.is_ok() { "ok" } else { "err" };
    print_metric_table(
        "symbolic_result_tensor",
        [
            ("status", symbolic_result_status.to_owned()),
            ("elapsed", format_duration(symbolic_result_elapsed)),
        ],
    );
    summary.push(
        "symbolic_result_tensor",
        symbolic_result_status,
        format_duration(symbolic_result_elapsed),
        "",
        "",
        "",
    );

    let symbolic_root = match symbolic_result.wrap_err("failed to read symbolic tensor result")? {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(tensor) => tensor.expression.clone(),
    };
    let symbolic_aliased = match &symbolic_aliases {
        Some(aliases) => symbolic_net.aliased_atom(aliases, symbolic_root),
        None => AliasedAtom::from(symbolic_root),
    };
    let symbolic_operations = symbolic_aliased.count_operations();
    print_metric_table(
        "symbolic_result_aliased_root",
        [
            ("aliases", symbolic_aliased.get_aliases().len().to_string()),
            ("terms", symbolic_aliased.get_root().nterms().to_string()),
            (
                "bytes",
                symbolic_aliased
                    .get_root()
                    .as_view()
                    .get_byte_size()
                    .to_string(),
            ),
            ("adds", symbolic_operations.additions.to_string()),
            ("muls", symbolic_operations.multiplications.to_string()),
        ],
    );
    summary.push(
        "symbolic_result_aliased_root",
        "done",
        "",
        symbolic_aliased.get_root().nterms().to_string(),
        symbolic_aliased
            .get_root()
            .as_view()
            .get_byte_size()
            .to_string(),
        format!(
            "aliases={} adds={} muls={}",
            symbolic_aliased.get_aliases().len(),
            symbolic_operations.additions,
            symbolic_operations.multiplications
        ),
    );

    if let Some(path) = &config.dump_aliased_path {
        let bytes = write_aliased_atom_dump(path, &symbolic_aliased)?;
        print_metric_table(
            "symbolic_aliased_dump",
            [
                ("path", path.display().to_string()),
                ("bytes", bytes.to_string()),
            ],
        );
        summary.push(
            "symbolic_aliased_dump",
            "written",
            "",
            "",
            bytes.to_string(),
            path.display().to_string(),
        );
    }

    let root_for_concrete = symbolic_aliased.get_root().clone();
    let concrete_input = if config.collect_tensors_before_concrete {
        print_metric_table(
            "collect_tensors_before_concrete",
            [("status", "start".to_owned())],
        );
        io::stdout().flush()?;
        let collect_started = Instant::now();
        let collected = root_for_concrete.collect_tensors();
        let collect_elapsed = collect_started.elapsed();
        print_metric_table(
            "collect_tensors_before_concrete",
            [
                ("status", "done".to_owned()),
                ("elapsed", format_duration(collect_elapsed)),
                ("terms", collected.nterms().to_string()),
                ("bytes", collected.as_view().get_byte_size().to_string()),
            ],
        );
        summary.push(
            "collect_tensors_before_concrete",
            "done",
            format_duration(collect_elapsed),
            collected.nterms().to_string(),
            collected.as_view().get_byte_size().to_string(),
            "",
        );
        collected
    } else {
        summary.push("collect_tensors_before_concrete", "skipped", "", "", "", "");
        root_for_concrete
    };
    print_metric_table(
        "concrete_parse_input",
        [
            ("source", "symbolic_aliased_root".to_owned()),
            (
                "aliases_retained",
                symbolic_aliased.get_aliases().len().to_string(),
            ),
            ("terms", concrete_input.nterms().to_string()),
            (
                "bytes",
                concrete_input.as_view().get_byte_size().to_string(),
            ),
        ],
    );
    summary.push(
        "concrete_parse_input",
        "aliases_retained",
        "",
        concrete_input.nterms().to_string(),
        concrete_input.as_view().get_byte_size().to_string(),
        format!("aliases={}", symbolic_aliased.get_aliases().len()),
    );

    print_metric_table("concrete_network_parse", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let concrete_parse_started = Instant::now();
    let mut concrete_net = concrete_input
        .parse_into_net()
        .wrap_err("failed to parse aliased symbolic result into concrete tensor network")?;
    let concrete_parse_elapsed = concrete_parse_started.elapsed();
    print_metric_table(
        "concrete_network_parse",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(concrete_parse_elapsed)),
            ("nodes", concrete_net.graph.n_nodes().to_string()),
            ("edges", concrete_net.graph.graph.n_edges().to_string()),
            ("tensors", concrete_net.store.tensors.len().to_string()),
        ],
    );
    summary.push(
        "concrete_network_parse",
        "done",
        format_duration(concrete_parse_elapsed),
        "",
        "",
        format!(
            "nodes={} edges={} tensors={}",
            concrete_net.graph.n_nodes(),
            concrete_net.graph.graph.n_edges(),
            concrete_net.store.tensors.len()
        ),
    );

    let concrete_scalar_store_stats =
        print_scalar_store_stats(&concrete_net, config.scalar_details);
    summary.push(
        "concrete_scalar_store",
        "done",
        "",
        concrete_scalar_store_stats.total_terms.to_string(),
        concrete_scalar_store_stats.total_bytes.to_string(),
        format!(
            "count={} top_level_sums={} max_terms={} max_bytes={}",
            concrete_scalar_store_stats.count,
            concrete_scalar_store_stats.top_level_sums,
            concrete_scalar_store_stats.max_terms,
            concrete_scalar_store_stats.max_bytes
        ),
    );

    let concrete_tensor_entry_stats =
        print_tensor_entry_scalar_stats("concrete_tensor_entry_scalars", &concrete_net);
    summary.push(
        "concrete_tensor_entry_scalars",
        "before_execute",
        "",
        concrete_tensor_entry_stats.total_terms.to_string(),
        concrete_tensor_entry_stats.total_bytes.to_string(),
        concrete_tensor_entry_stats.summary_details(),
    );

    print_metric_table(
        "concrete_scalar_aliases",
        [
            ("skipped", "true".to_owned()),
            ("reason", "symbolic_aliases_retained".to_owned()),
        ],
    );
    summary.push(
        "concrete_scalar_aliases",
        "skipped",
        "",
        "",
        "",
        "symbolic_aliases_retained",
    );

    if let Some(path) = &config.dump_network_path {
        fs::write(path, concrete_net.dot_pretty())
            .with_context(|| format!("failed to write {}", path.display()))?;
        print_metric_table(
            "concrete_network_dump",
            [("path", path.display().to_string())],
        );
        summary.push(
            "concrete_network_dump",
            "written",
            "",
            "",
            "",
            path.display().to_string(),
        );
    }

    print_metric_table("concrete_network_execute", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let concrete_execute_started = Instant::now();
    execute_concrete_network!(concrete_net, config.contraction_order)
        .wrap_err("failed to execute concrete tensor network")?;
    let concrete_execute_elapsed = concrete_execute_started.elapsed();
    print_metric_table(
        "concrete_network_execute",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(concrete_execute_elapsed)),
        ],
    );
    summary.push(
        "concrete_network_execute",
        "done",
        format_duration(concrete_execute_elapsed),
        "",
        "",
        "",
    );
    let concrete_tensor_entry_stats = print_tensor_entry_scalar_stats(
        "concrete_tensor_entry_scalars_after_execute",
        &concrete_net,
    );
    summary.push(
        "concrete_tensor_entry_scalars",
        "after_execute",
        "",
        concrete_tensor_entry_stats.total_terms.to_string(),
        concrete_tensor_entry_stats.total_bytes.to_string(),
        concrete_tensor_entry_stats.summary_details(),
    );

    let concrete_result_started = Instant::now();
    let concrete_result = concrete_net.result_scalar();
    let concrete_result_elapsed = concrete_result_started.elapsed();
    let concrete_result_status = if concrete_result.is_ok() { "ok" } else { "err" };
    print_metric_table(
        "concrete_result_scalar",
        [
            ("status", concrete_result_status.to_owned()),
            ("elapsed", format_duration(concrete_result_elapsed)),
        ],
    );
    summary.push(
        "concrete_result_scalar",
        concrete_result_status,
        format_duration(concrete_result_elapsed),
        "",
        "",
        "",
    );

    let concrete_root = match concrete_result {
        Ok(ExecutionResult::One) => Atom::num(1),
        Ok(ExecutionResult::Zero) => Atom::Zero,
        Ok(ExecutionResult::Val(value)) => value.into_owned(),
        Err(error) => return Err(error).wrap_err("failed to read concrete scalar result"),
    };

    if symbolic_aliased.get_aliases().is_empty() {
        print_metric_table(
            "concrete_result_atom",
            [
                ("terms", concrete_root.nterms().to_string()),
                ("bytes", concrete_root.as_view().get_byte_size().to_string()),
            ],
        );
        summary.push(
            "concrete_result_atom",
            "done",
            "",
            concrete_root.nterms().to_string(),
            concrete_root.as_view().get_byte_size().to_string(),
            "",
        );
    } else {
        let mut concrete_aliased = AliasedAtom::from(concrete_root);
        for (alias, original) in symbolic_aliased.get_aliases() {
            concrete_aliased = concrete_aliased.add_alias(alias.clone(), original.clone());
        }
        print_metric_table(
            "concrete_result_aliased_root",
            [
                ("aliases", concrete_aliased.get_aliases().len().to_string()),
                ("terms", concrete_aliased.get_root().nterms().to_string()),
                (
                    "bytes",
                    concrete_aliased
                        .get_root()
                        .as_view()
                        .get_byte_size()
                        .to_string(),
                ),
            ],
        );
        summary.push(
            "concrete_result_aliased_root",
            "done",
            "",
            concrete_aliased.get_root().nterms().to_string(),
            concrete_aliased
                .get_root()
                .as_view()
                .get_byte_size()
                .to_string(),
            format!("aliases={}", concrete_aliased.get_aliases().len()),
        );
        if config.resolve_aliases {
            let resolve_started = Instant::now();
            let resolved = concrete_aliased.into_inner();
            let resolve_elapsed = resolve_started.elapsed();
            print_metric_table(
                "concrete_result_alias_resolve",
                [
                    ("elapsed", format_duration(resolve_elapsed)),
                    ("terms", resolved.nterms().to_string()),
                    ("bytes", resolved.as_view().get_byte_size().to_string()),
                ],
            );
            summary.push(
                "concrete_result_alias_resolve",
                "done",
                format_duration(resolve_elapsed),
                resolved.nterms().to_string(),
                resolved.as_view().get_byte_size().to_string(),
                "",
            );
        } else {
            print_metric_table(
                "concrete_result_alias_resolve",
                [("skipped", "true".to_owned())],
            );
            summary.push("concrete_result_alias_resolve", "skipped", "", "", "", "");
        }
    }

    Ok(())
}

fn run_symbolic_net(config: &Config, selected: &Atom, summary: &mut SummaryTable) -> Result<()> {
    print_metric_table("symbolic_network_parse", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let network_parse_started = Instant::now();
    let mut net = parse_symbolic_network(config, selected)?;
    let network_parse_elapsed = network_parse_started.elapsed();
    print_metric_table(
        "symbolic_network_parse",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(network_parse_elapsed)),
            ("nodes", net.graph.n_nodes().to_string()),
            ("edges", net.graph.graph.n_edges().to_string()),
            ("tensors", net.store.tensors.len().to_string()),
            (
                "strict_filter",
                format!("{:?}", config.symbolic_strict_filter),
            ),
        ],
    );
    summary.push(
        "symbolic_network_parse",
        "done",
        format_duration(network_parse_elapsed),
        "",
        "",
        format!(
            "nodes={} edges={} tensors={} strict_filter={:?}",
            net.graph.n_nodes(),
            net.graph.graph.n_edges(),
            net.store.tensors.len(),
            config.symbolic_strict_filter
        ),
    );

    let scalar_store_stats = print_scalar_store_stats(&net, config.scalar_details);
    summary.push(
        "scalar_store",
        "done",
        "",
        scalar_store_stats.total_terms.to_string(),
        scalar_store_stats.total_bytes.to_string(),
        format!(
            "count={} top_level_sums={} max_terms={} max_bytes={}",
            scalar_store_stats.count,
            scalar_store_stats.top_level_sums,
            scalar_store_stats.max_terms,
            scalar_store_stats.max_bytes
        ),
    );

    if config.unit_scalars || config.kept_scalars.is_some() {
        let replacement_stats =
            replace_scalars(&mut net, config.unit_scalars, config.kept_scalars.as_ref())?;
        summary.push(
            "scalar_replacement",
            replacement_stats.mode,
            "",
            "",
            "",
            format!(
                "replaced={} kept={}",
                replacement_stats.replaced, replacement_stats.kept
            ),
        );
    }

    let scalar_aliases = if let Some(threshold) = config.alias_scalar_threshold {
        let aliases =
            net.alias_scalar_refs(|_, scalar| scalar.as_view().get_byte_size() >= threshold);
        print_metric_table(
            "scalar_aliases",
            [
                ("threshold_bytes", threshold.to_string()),
                ("created", aliases.aliases_created().to_string()),
                ("terms", aliases.aliased_terms().to_string()),
                ("bytes", aliases.aliased_bytes().to_string()),
                ("max_bytes", aliases.max_aliased_bytes().to_string()),
            ],
        );
        summary.push(
            "scalar_aliases",
            "enabled",
            "",
            aliases.aliased_terms().to_string(),
            aliases.aliased_bytes().to_string(),
            format!(
                "threshold_bytes={threshold} created={} max_bytes={}",
                aliases.aliases_created(),
                aliases.max_aliased_bytes()
            ),
        );
        Some(aliases)
    } else {
        summary.push("scalar_aliases", "disabled", "", "", "", "");
        None
    };

    if let Some(path) = &config.dump_scalar_aliases_md_path {
        let Some(aliases) = &scalar_aliases else {
            bail!("--dump-scalar-aliases-md requires --alias-scalars");
        };
        let bytes = write_scalar_aliases_markdown(
            path,
            &config.input_path,
            config.alias_scalar_threshold.unwrap_or_default(),
            aliases,
            &net,
        )?;
        print_metric_table(
            "scalar_aliases_md",
            [
                ("path", path.display().to_string()),
                ("bytes", bytes.to_string()),
            ],
        );
        summary.push(
            "scalar_aliases_md",
            "written",
            "",
            "",
            bytes.to_string(),
            path.display().to_string(),
        );
    }

    if let Some(path) = &config.dump_network_path {
        fs::write(path, net.snapshot_dot())
            .with_context(|| format!("failed to write {}", path.display()))?;
        print_metric_table(
            "symbolic_network_dump",
            [("path", path.display().to_string())],
        );
        summary.push(
            "symbolic_network_dump",
            "written",
            "",
            "",
            "",
            path.display().to_string(),
        );
    }

    if !config.execute {
        print_metric_table("symbolic_network_execute", [("skipped", "true".to_owned())]);
        summary.push("symbolic_network_execute", "skipped", "", "", "", "");
        return Ok(());
    }

    print_metric_table("symbolic_network_execute", [("status", "start".to_owned())]);
    io::stdout().flush()?;
    let execute_started = Instant::now();
    let lib = DummyLibrary::new();
    net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &Wrap {})
        .wrap_err("failed to execute symbolic tensor network")?;
    let execute_elapsed = execute_started.elapsed();
    print_metric_table(
        "symbolic_network_execute",
        [
            ("status", "done".to_owned()),
            ("elapsed", format_duration(execute_elapsed)),
        ],
    );
    summary.push(
        "symbolic_network_execute",
        "done",
        format_duration(execute_elapsed),
        "",
        "",
        "",
    );

    let result_started = Instant::now();
    let result = net.result_tensor(&lib);
    let result_elapsed = result_started.elapsed();
    let result_status = if result.is_ok() { "ok" } else { "err" };
    print_metric_table(
        "symbolic_result_tensor",
        [
            ("status", result_status.to_owned()),
            ("elapsed", format_duration(result_elapsed)),
        ],
    );
    summary.push(
        "symbolic_result_tensor",
        result_status,
        format_duration(result_elapsed),
        "",
        "",
        "",
    );

    let root = match result.wrap_err("failed to read symbolic tensor result")? {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(tensor) => tensor.expression.clone(),
    };
    let aliased = match &scalar_aliases {
        Some(aliases) => net.aliased_atom(aliases, root),
        None => AliasedAtom::from(root),
    };
    let operations = aliased.count_operations();
    print_metric_table(
        "symbolic_result_aliased_root",
        [
            ("aliases", aliased.get_aliases().len().to_string()),
            ("terms", aliased.get_root().nterms().to_string()),
            (
                "bytes",
                aliased.get_root().as_view().get_byte_size().to_string(),
            ),
            ("adds", operations.additions.to_string()),
            ("muls", operations.multiplications.to_string()),
        ],
    );
    summary.push(
        "symbolic_result_aliased_root",
        "done",
        "",
        aliased.get_root().nterms().to_string(),
        aliased.get_root().as_view().get_byte_size().to_string(),
        format!(
            "aliases={} adds={} muls={}",
            aliased.get_aliases().len(),
            operations.additions,
            operations.multiplications
        ),
    );

    if let Some(path) = &config.dump_aliased_path {
        let bytes = write_aliased_atom_dump(path, &aliased)?;
        print_metric_table(
            "aliased_dump",
            [
                ("path", path.display().to_string()),
                ("bytes", bytes.to_string()),
            ],
        );
        summary.push(
            "aliased_dump",
            "written",
            "",
            "",
            bytes.to_string(),
            path.display().to_string(),
        );
    }

    println!("{}", aliased.log_print(Some(120)));
    if config.resolve_aliases {
        let resolve_started = Instant::now();
        let resolved = aliased.into_inner();
        let resolve_elapsed = resolve_started.elapsed();
        print_metric_table(
            "symbolic_result_alias_resolve",
            [
                ("elapsed", format_duration(resolve_elapsed)),
                ("terms", resolved.nterms().to_string()),
                ("bytes", resolved.as_view().get_byte_size().to_string()),
            ],
        );
        summary.push(
            "symbolic_result_alias_resolve",
            "done",
            format_duration(resolve_elapsed),
            resolved.nterms().to_string(),
            resolved.as_view().get_byte_size().to_string(),
            "",
        );
    } else {
        print_metric_table(
            "symbolic_result_alias_resolve",
            [("skipped", "true".to_owned())],
        );
        summary.push("symbolic_result_alias_resolve", "skipped", "", "", "", "");
    }

    Ok(())
}

#[derive(serde::Serialize)]
struct AliasedAtomDump {
    format: &'static str,
    version: usize,
    root: String,
    aliases: Vec<AliasedAtomDumpEntry>,
}

#[derive(serde::Serialize)]
struct AliasedAtomDumpEntry {
    alias: String,
    original: String,
}

fn write_aliased_atom_dump(path: &Path, aliased: &AliasedAtom) -> Result<usize> {
    let mut aliases = aliased
        .get_aliases()
        .iter()
        .map(|(alias, original)| AliasedAtomDumpEntry {
            alias: alias.to_plain_string(),
            original: original.to_plain_string(),
        })
        .collect::<Vec<_>>();
    aliases.sort_by(|left, right| left.alias.cmp(&right.alias));

    let dump = AliasedAtomDump {
        format: "gammaloop-aliased-atom",
        version: 1,
        root: aliased.get_root().to_plain_string(),
        aliases,
    };
    let bytes = serde_json::to_vec(&dump)?;

    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    fs::write(path, &bytes).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(bytes.len())
}

fn write_scalar_aliases_markdown<T, K, FK, A>(
    path: &Path,
    input_path: &Path,
    threshold: usize,
    aliases: &ScalarAliases,
    net: &Network<NetworkStore<T, Atom>, K, FK, A>,
) -> Result<usize> {
    let mut markdown = String::new();
    markdown.push_str("# BNL Scalar Alias Captures\n\n");
    markdown.push_str("Source atom:\n\n");
    markdown.push_str("```text\n");
    markdown.push_str(&input_path.display().to_string());
    markdown.push_str("\n```\n\n");
    markdown.push_str("This file only shows `log_print(Some(120))` fields.\n\n");
    markdown.push_str(&format!("Alias threshold: `{threshold}` bytes.\n\n"));

    for (position, index) in aliases.aliased_indices().enumerate() {
        let scalar = net
            .store
            .scalar
            .get(index)
            .ok_or_else(|| eyre!("scalar alias index {index} is out of bounds"))?;
        markdown.push_str(&format!("## Alias {}: `scalar({index})`\n\n", position + 1));
        markdown.push_str(&format!(
            "- terms: `{}`\n- bytes: `{}`\n\n",
            scalar.nterms(),
            scalar.as_view().get_byte_size()
        ));
        markdown.push_str("```text\n");
        markdown.push_str(&scalar.log_print(Some(120)).to_string());
        markdown.push_str("\n```\n\n");
    }

    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    fs::write(path, markdown.as_bytes())
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(markdown.len())
}

fn print_scalar_store_stats<T, K, FK, A>(
    net: &Network<NetworkStore<T, Atom>, K, FK, A>,
    print_details: bool,
) -> ScalarStoreStats {
    let mut top_level_sums = 0;
    let mut total_terms = 0;
    let mut max_terms = 0;
    let mut total_bytes = 0;
    let mut max_bytes = 0;
    let mut scalar_detail_rows = Vec::new();

    for (index, scalar) in net.store.scalar.iter().enumerate() {
        let view = scalar.as_view();
        let terms = scalar.nterms();
        let bytes = view.get_byte_size();
        top_level_sums += usize::from(matches!(view, AtomView::Add(_)));
        total_terms += terms;
        max_terms = max_terms.max(terms);
        total_bytes += bytes;
        max_bytes = max_bytes.max(bytes);
        if print_details {
            scalar_detail_rows.push(vec![
                index.to_string(),
                terms.to_string(),
                bytes.to_string(),
                matches!(view, AtomView::Add(_)).to_string(),
            ]);
        }
    }

    if print_details {
        print_table(
            ["index", "terms", "bytes", "top_level_sum"],
            scalar_detail_rows,
        );
    }

    print_metric_table(
        "scalar_store",
        [
            ("count", net.store.scalar.len().to_string()),
            ("top_level_sums", top_level_sums.to_string()),
            ("total_terms", total_terms.to_string()),
            ("max_terms", max_terms.to_string()),
            ("total_bytes", total_bytes.to_string()),
            ("max_bytes", max_bytes.to_string()),
        ],
    );

    ScalarStoreStats {
        count: net.store.scalar.len(),
        top_level_sums,
        total_terms,
        max_terms,
        total_bytes,
        max_bytes,
    }
}

fn print_tensor_entry_scalar_stats(stage: &str, net: &ParsingNet) -> TensorEntryStats {
    let mut stats = TensorEntryStats {
        tensor_count: net.store.tensors.len(),
        param_tensor_count: 0,
        concrete_tensor_count: 0,
        entry_count: 0,
        total_terms: 0,
        max_terms: 0,
        total_bytes: 0,
        max_bytes: 0,
        max_tensor_index: 0,
        max_entry_index: 0,
    };
    let mut top_tensors = Vec::new();

    for (tensor_index, tensor) in net.store.tensors.iter().enumerate() {
        let ParamOrConcrete::Param(tensor) = tensor else {
            stats.concrete_tensor_count += 1;
            continue;
        };
        stats.param_tensor_count += 1;

        let entries = tensor.data();
        let mut tensor_terms = 0;
        let mut tensor_bytes = 0;
        let mut tensor_max_bytes = 0;

        for (entry_index, entry) in entries.iter().enumerate() {
            let terms = entry.nterms();
            let bytes = entry.as_view().get_byte_size();
            stats.entry_count += 1;
            stats.total_terms += terms;
            stats.total_bytes += bytes;
            stats.max_terms = stats.max_terms.max(terms);

            tensor_terms += terms;
            tensor_bytes += bytes;
            tensor_max_bytes = tensor_max_bytes.max(bytes);

            if bytes > stats.max_bytes {
                stats.max_bytes = bytes;
                stats.max_tensor_index = tensor_index;
                stats.max_entry_index = entry_index;
            }
        }

        top_tensors.push((
            tensor_index,
            entries.len(),
            tensor_terms,
            tensor_bytes,
            tensor_max_bytes,
        ));
    }

    top_tensors.sort_by_key(|tensor| std::cmp::Reverse(tensor.3));

    print_metric_table(
        stage,
        [
            ("tensors", stats.tensor_count.to_string()),
            ("param_tensors", stats.param_tensor_count.to_string()),
            ("concrete_tensors", stats.concrete_tensor_count.to_string()),
            ("entries", stats.entry_count.to_string()),
            ("total_terms", stats.total_terms.to_string()),
            ("max_terms", stats.max_terms.to_string()),
            ("total_bytes", stats.total_bytes.to_string()),
            ("max_bytes", stats.max_bytes.to_string()),
            ("max_tensor", stats.max_tensor_index.to_string()),
            ("max_entry", stats.max_entry_index.to_string()),
        ],
    );

    let top_rows = top_tensors
        .into_iter()
        .take(8)
        .map(|(tensor_index, entries, terms, bytes, max_entry_bytes)| {
            vec![
                tensor_index.to_string(),
                entries.to_string(),
                terms.to_string(),
                bytes.to_string(),
                max_entry_bytes.to_string(),
            ]
        })
        .collect::<Vec<_>>();
    if !top_rows.is_empty() {
        println!("{stage}_top_tensors");
        print_table(
            ["tensor", "entries", "terms", "bytes", "max_entry_bytes"],
            top_rows,
        );
    }

    stats
}

fn replace_scalars<T, K, FK, A>(
    net: &mut Network<NetworkStore<T, Atom>, K, FK, A>,
    unit_scalars: bool,
    kept_scalars: Option<&BTreeSet<usize>>,
) -> Result<ScalarReplacementStats> {
    if let Some(indices) = kept_scalars {
        for index in indices {
            if *index >= net.store.scalar.len() {
                bail!(
                    "scalar index {index} is out of range for {} scalars",
                    net.store.scalar.len()
                );
            }
        }
    }

    let mut replaced = 0usize;
    let mut kept = 0usize;
    for (index, scalar) in net.store.scalar.iter_mut().enumerate() {
        let should_replace =
            unit_scalars || kept_scalars.is_some_and(|indices| !indices.contains(&index));
        if should_replace {
            *scalar = Atom::num(1);
            replaced += 1;
        } else {
            kept += 1;
        }
    }

    let mode = if unit_scalars {
        "unit_all"
    } else {
        "keep_selected"
    };

    print_metric_table(
        "scalar_replacement",
        [
            ("mode", mode.to_owned()),
            ("replaced", replaced.to_string()),
            ("kept", kept.to_string()),
        ],
    );
    Ok(ScalarReplacementStats {
        mode,
        replaced,
        kept,
    })
}

fn select_terms(atom: &Atom, selection: TermSelection) -> Result<Atom> {
    match selection {
        TermSelection::All => Ok(atom.clone()),
        TermSelection::First(count) => select_term_range(atom, 0, count),
        TermSelection::Term(index) => select_term_range(atom, index, index + 1),
        TermSelection::Range { start, end } => select_term_range(atom, start, end),
    }
}

fn select_term_range(atom: &Atom, start: usize, end: usize) -> Result<Atom> {
    if end <= start {
        bail!("term selection must contain at least one term");
    }

    match atom.as_view() {
        AtomView::Add(add) => {
            let total_terms = add.iter().count();
            if start >= total_terms {
                bail!("selection starts at term {start}, but atom has {total_terms} terms");
            }
            if end > total_terms {
                bail!("selection ends at term {end}, but atom has {total_terms} terms");
            }

            Ok(add
                .iter()
                .enumerate()
                .filter(|(index, _)| start <= *index && *index < end)
                .map(|(_, term)| term.to_owned())
                .fold(Atom::Zero, |sum, term| sum + term))
        }
        _ if start == 0 && end == 1 => Ok(atom.clone()),
        _ => bail!("term selection requires a top-level sum; the input atom has one term"),
    }
}
