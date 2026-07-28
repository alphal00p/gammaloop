use std::{env, fs, ops::Deref, path::PathBuf, time::Instant};

use color_eyre::{
    Result,
    eyre::{Context, eyre},
};
use gammalooprs::{
    initialisation::initialise,
    numerator::symbolica_ext::AtomCoreExt,
    utils::{FUN_LIB, TENSORLIB},
};
use spenso::network::{ExecutionResult, MinResultRank, Sequential, SequentialRef, SmallestDegree};

const DEFAULT_INPUT: &str = "/tmp/gl00_spenso_pre_network.atom";

#[derive(Debug)]
struct Options {
    input: PathBuf,
    skip_execute: bool,
    strategy: Strategy,
}

#[derive(Debug, Clone, Copy)]
enum Strategy {
    MinResultRank,
    SmallestDegree,
    RefSmallestDegree,
}

impl Options {
    fn parse() -> Result<Self> {
        let mut input = PathBuf::from(DEFAULT_INPUT);
        let mut skip_execute = false;
        let mut strategy = Strategy::MinResultRank;

        let mut args = env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--input" => {
                    input = args
                        .next()
                        .map(PathBuf::from)
                        .ok_or_else(|| eyre!("--input requires a path"))?;
                }
                "--skip-execute" => skip_execute = true,
                "--strategy" => {
                    let value = args
                        .next()
                        .ok_or_else(|| eyre!("--strategy requires a value"))?;
                    strategy = match value.as_str() {
                        "sequential-min-result-rank" => Strategy::MinResultRank,
                        "sequential-smallest-degree" => Strategy::SmallestDegree,
                        "sequential-ref-smallest-degree" => Strategy::RefSmallestDegree,
                        _ => {
                            return Err(eyre!(
                                "unknown strategy '{value}', expected sequential-min-result-rank, sequential-smallest-degree, or sequential-ref-smallest-degree"
                            ));
                        }
                    };
                }
                "-h" | "--help" => {
                    print_help();
                    std::process::exit(0);
                }
                _ if !arg.starts_with('-') => input = PathBuf::from(arg),
                _ => return Err(eyre!("unknown argument '{arg}'")),
            }
        }

        Ok(Self {
            input,
            skip_execute,
            strategy,
        })
    }
}

fn print_help() {
    println!(
        "Usage: cargo run -p gammalooprs --example gl00_spenso_execute -- [OPTIONS] [INPUT]\n\
\n\
Options:\n\
  --input <PATH>       Atom dumped by GAMMALOOP_DUMP_EVALUATOR_PRE_NETWORK_PARSE\n\
  --skip-execute       Stop after parsing the spenso network\n\
  --strategy <NAME>    sequential-min-result-rank (default), sequential-smallest-degree,\n\
                       or sequential-ref-smallest-degree\n\
  -h, --help           Show this help"
    );
}

fn main() -> Result<()> {
    initialise()?;
    drop(TENSORLIB.read().unwrap());

    let options = Options::parse()?;
    let input_started = Instant::now();
    let input = fs::read_to_string(&options.input)
        .with_context(|| format!("failed to read {}", options.input.display()))?;
    println!(
        "read input: {} bytes in {:.3}s",
        input.len(),
        input_started.elapsed().as_secs_f64()
    );

    let atom_started = Instant::now();
    let atom = symbolica::try_parse!(input.as_str()).map_err(|err| eyre!(err.to_string()))?;
    println!(
        "parsed atom: {} terms, {} bytes in {:.3}s",
        atom.nterms(),
        atom.as_view().get_byte_size(),
        atom_started.elapsed().as_secs_f64()
    );

    let network_started = Instant::now();
    let mut net = atom.parse_into_net()?;
    println!(
        "parsed spenso network in {:.3}s",
        network_started.elapsed().as_secs_f64()
    );

    if options.skip_execute {
        return Ok(());
    }

    let execute_started = Instant::now();
    match options.strategy {
        Strategy::MinResultRank => {
            net.execute::<Sequential, MinResultRank, _, _, _>(
                TENSORLIB.read().unwrap().deref(),
                FUN_LIB.deref(),
            )?;
        }
        Strategy::SmallestDegree => {
            net.execute::<Sequential, SmallestDegree, _, _, _>(
                TENSORLIB.read().unwrap().deref(),
                FUN_LIB.deref(),
            )?;
        }
        Strategy::RefSmallestDegree => {
            net.execute::<SequentialRef, SmallestDegree, _, _, _>(
                TENSORLIB.read().unwrap().deref(),
                FUN_LIB.deref(),
            )?;
        }
    }
    println!(
        "executed network with {:?} in {:.3}s",
        options.strategy,
        execute_started.elapsed().as_secs_f64()
    );

    let result = net.result_scalar()?;
    match result {
        ExecutionResult::One => println!("result: one"),
        ExecutionResult::Zero => println!("result: zero"),
        ExecutionResult::Val(value) => println!(
            "result atom: {} bytes, {} terms",
            value.as_ref().as_view().get_byte_size(),
            value.nterms()
        ),
    }

    Ok(())
}
