use std::{
    collections::BTreeSet,
    env, fs,
    io::{self, Write},
    ops::Deref,
    path::PathBuf,
    time::Instant,
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
use spenso::network::{MinResultRank, Sequential};
use symbolica::{
    atom::{Atom, AtomView},
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};

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
}

fn main() -> Result<()> {
    test_initialise()?;

    let config = Config::from_args()?;
    let input = fs::read_to_string(&config.input_path)
        .with_context(|| format!("failed to read {}", config.input_path.display()))?;

    println!(
        "input\tpath={}\tbytes={}",
        config.input_path.display(),
        input.len()
    );

    let parse_started = Instant::now();
    let atom =
        Atom::parse_with_default_namespace(wrap_input!(&input), SymbolicaParseSettings::default())
            .map_err(|error| eyre!("failed to parse input atom: {error}"))?;
    print_atom_stats("symbolica_parse", parse_started, &atom);
    if config.print_log {
        println!("input_log_print_120\n{}", atom.log_print(Some(120)));
    }

    let selected = select_terms(&atom, config.selection)?;
    if !matches!(config.selection, TermSelection::All) {
        println!("selection\tselector={}", config.selection.selected_label());
        print_atom_stats("selected_atom", Instant::now(), &selected);
        if config.print_log {
            println!("selected_log_print_120\n{}", selected.log_print(Some(120)));
        }
    }

    println!("network_parse_start");
    io::stdout().flush()?;
    let network_parse_started = Instant::now();
    let mut net = selected
        .parse_into_net()
        .wrap_err("failed to parse evaluator atom into tensor network")?;
    println!(
        "network_parse_done\telapsed_ms={:.3}\tnodes={}\tedges={}\ttensors={}",
        network_parse_started.elapsed().as_secs_f64() * 1000.0,
        net.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len()
    );
    print_scalar_store_stats(&net, config.scalar_details);

    if config.unit_scalars || config.kept_scalars.is_some() {
        replace_scalars(&mut net, config.unit_scalars, config.kept_scalars.as_ref())?;
    }

    if let Some(path) = config.dump_network_path {
        fs::write(&path, net.dot_pretty())
            .with_context(|| format!("failed to write {}", path.display()))?;
        println!("network_dump\tpath={}", path.display());
    }

    if !config.execute {
        println!("network_execute\tskipped=true");
        return Ok(());
    }

    println!("network_execute_start");
    io::stdout().flush()?;
    let execute_started = Instant::now();
    net.execute::<Sequential, MinResultRank, _, _, _>(
        TENSORLIB.read().unwrap().deref(),
        FUN_LIB.deref(),
    )
    .wrap_err("failed to execute evaluator atom tensor network")?;
    println!(
        "network_execute_done\telapsed_ms={:.3}",
        execute_started.elapsed().as_secs_f64() * 1000.0,
    );

    let result_started = Instant::now();
    let result = net.result_scalar();
    println!(
        "result_scalar\tstatus={}\telapsed_ms={:.3}",
        if result.is_ok() { "ok" } else { "err" },
        result_started.elapsed().as_secs_f64() * 1000.0,
    );

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
        "usage: cargo run -p gammalooprs --example bnl_evaluator_atom_mwe --profile dev-optim -- [INPUT] [all|first:N|term:N|range:START..END] [--skip-execute] [--print-log] [--scalar-details] [--unit-scalars|--keep-scalars INDICES] [--dump-network PATH]"
    );
    println!("term indices are zero-based; range END is exclusive");
    println!(
        "--keep-scalars keeps only the listed scalar-store indices and replaces the rest with 1"
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

fn print_atom_stats(label: &str, started: Instant, atom: &Atom) {
    println!(
        "{label}\telapsed_ms={:.3}\tterms={}\tbytes={}",
        started.elapsed().as_secs_f64() * 1000.0,
        atom.nterms(),
        atom.as_view().get_byte_size()
    );
}

fn print_scalar_store_stats(net: &gammalooprs::numerator::ParsingNet, print_details: bool) {
    let mut top_level_sums = 0;
    let mut total_terms = 0;
    let mut max_terms = 0;
    let mut total_bytes = 0;
    let mut max_bytes = 0;

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
            println!(
                "scalar_detail\tindex={index}\tterms={terms}\tbytes={bytes}\ttop_level_sum={}",
                matches!(view, AtomView::Add(_))
            );
        }
    }

    println!(
        "scalar_store\tcount={}\ttop_level_sums={top_level_sums}\ttotal_terms={total_terms}\tmax_terms={max_terms}\ttotal_bytes={total_bytes}\tmax_bytes={max_bytes}",
        net.store.scalar.len()
    );
}

fn replace_scalars(
    net: &mut gammalooprs::numerator::ParsingNet,
    unit_scalars: bool,
    kept_scalars: Option<&BTreeSet<usize>>,
) -> Result<()> {
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

    println!(
        "scalar_replacement\tmode={}\treplaced={replaced}\tkept={kept}",
        if unit_scalars {
            "unit_all"
        } else {
            "keep_selected"
        }
    );
    Ok(())
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
