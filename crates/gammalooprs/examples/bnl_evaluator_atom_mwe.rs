use std::{
    collections::BTreeSet,
    env, fs,
    io::{self, Write},
    ops::Deref,
    path::PathBuf,
    time::{Duration, Instant},
};

use color_eyre::{
    Result,
    eyre::{WrapErr, bail, eyre},
};
use gammalooprs::{
    initialisation::test_initialise,
    numerator::symbolica_ext::AtomCoreExt,
    utils::{FUN_LIB, TENSORLIB, symbolica_ext::LogPrint},
};
use spenso::network::{ExecutionResult, MinResultRank, Sequential};
use symbolica::{
    atom::{Atom, AtomView},
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};
use tabled::{builder::Builder, settings::Style};

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
    net.execute::<Sequential, MinResultRank, _, _, _>(
        TENSORLIB.read().unwrap().deref(),
        FUN_LIB.deref(),
    )
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
        })
    }
}

fn looks_like_selector(value: &str) -> bool {
    value == "all"
        || value.starts_with("first:")
        || value.starts_with("term:")
        || value.starts_with("range:")
}

fn print_usage() {
    println!(
        "usage: cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- [INPUT] [all|first:N|term:N|range:START..END] [--skip-execute] [--print-log] [--scalar-details] [--unit-scalars|--keep-scalars INDICES] [--alias-scalars BYTES] [--dump-network PATH]"
    );
    println!("term indices are zero-based; range END is exclusive");
    println!(
        "--keep-scalars keeps only the listed scalar-store indices and replaces the rest with 1"
    );
    println!("--alias-scalars replaces scalar refs at or above BYTES with scalar(INDEX) aliases");
    println!("--resolve-aliases expands aliased results back into a full Symbolica atom");
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

fn print_scalar_store_stats(
    net: &gammalooprs::numerator::ParsingNet,
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

fn replace_scalars(
    net: &mut gammalooprs::numerator::ParsingNet,
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
                .filter_map(|(index, term)| {
                    (start <= index && index < end).then(|| term.to_owned())
                })
                .fold(Atom::Zero, |sum, term| sum + term))
        }
        _ if start == 0 && end == 1 => Ok(atom.clone()),
        _ => bail!("term selection requires a top-level sum; the input atom has one term"),
    }
}
