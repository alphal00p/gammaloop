use std::{env, fs, ops::Deref, path::PathBuf, time::Instant};

use color_eyre::{
    Result,
    eyre::{WrapErr, eyre},
};
use gammalooprs::{
    initialisation::test_initialise,
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{FUN_LIB, GS, TENSORLIB},
};
use idenso::{
    color::ColorSimplifier,
    dirac::GammaSimplifier,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};
use spenso::{
    network::{MinResultRank, Sequential},
    shadowing::symbolica_utils::LogPrint,
};
use symbolica::{
    atom::{Atom, AtomCore},
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};

fn main() -> Result<()> {
    test_initialise()?;

    let input_path = env::args().nth(1).map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from(
            "examples/cli/BNL/profiling/integrated_uv_start_before_gamma_simplification.sym",
        )
    });

    let input = fs::read_to_string(&input_path)
        .with_context(|| format!("failed to read {}", input_path.display()))?;
    let parse_started = Instant::now();
    let atom =
        Atom::parse_with_default_namespace(wrap_input!(&input), SymbolicaParseSettings::default())
            .map_err(|error| eyre!("failed to parse input atom: {error}"))?;
    println!(
        "symbolica_parse\telapsed_ms={:.3}\tterms={}\tbytes={}",
        parse_started.elapsed().as_secs_f64() * 1000.0,
        atom.nterms(),
        atom.as_view().get_byte_size()
    );
    println!("input_log_print_120\n{}", atom.log_print(Some(120)));

    println!("raw_network_parse\tskipped=symbolic_gammalooprs_dim_not_concrete");

    let dim4_atom = atom.replace(GS.dim).with(Atom::num(4));
    run_network("dim4_raw", &dim4_atom)?;

    let algebra_started = Instant::now();
    let simplified = atom
        .simplify_color()
        .simplify_gamma()
        .simplify_metrics()
        .to_dots();
    println!(
        "evaluator_algebra\telapsed_ms={:.3}\tterms={}\tbytes={}",
        algebra_started.elapsed().as_secs_f64() * 1000.0,
        simplified.nterms(),
        simplified.as_view().get_byte_size()
    );
    println!(
        "evaluator_algebra_log_print_120\n{}",
        simplified.log_print(Some(120))
    );
    let dim4_simplified = simplified.replace(GS.dim).with(Atom::num(4));
    run_network("dim4_evaluator_algebra", &dim4_simplified)?;

    Ok(())
}

fn run_network(label: &str, atom: &Atom) -> Result<()> {
    let parse_started = Instant::now();
    let mut net = atom.parse_into_net().wrap_err("failed to parse network")?;
    println!(
        "{label}_network_parse\telapsed_ms={:.3}\tnodes={}\tedges={}\ttensors={}",
        parse_started.elapsed().as_secs_f64() * 1000.0,
        net.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len()
    );

    let execute_started = Instant::now();
    net.execute::<Sequential, MinResultRank, _, _, _>(
        TENSORLIB.read().unwrap().deref(),
        FUN_LIB.deref(),
    )
    .wrap_err("failed to execute network")?;
    println!(
        "{label}_network_execute\telapsed_ms={:.3}",
        execute_started.elapsed().as_secs_f64() * 1000.0,
    );

    Ok(())
}
