#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! color-eyre = { version = "0.6.5", default-features = false, features = ["capture-spantrace", "track-caller"] }
//! serde = { version = "1.0.228", features = ["derive"] }
//! serde_json = "1"
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "3638099c607d79da709989716c8dc9d5085364bd", default-features = false, features = ["serde"] }
//!
//! [patch.crates-io]
//! graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "3638099c607d79da709989716c8dc9d5085364bd" }
//! numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "3638099c607d79da709989716c8dc9d5085364bd" }
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "3638099c607d79da709989716c8dc9d5085364bd" }
//! ```

use std::{
    env, fs,
    io::{self, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use color_eyre::{
    eyre::{bail, eyre, WrapErr},
    Result,
};
use serde::Deserialize;
use symbolica::{atom::Atom, id::AliasedAtom, parser::ParseSettings, wrap_input};

#[derive(Debug)]
struct Config {
    input: PathBuf,
    output: Option<PathBuf>,
}

#[derive(Debug, Deserialize)]
struct AliasedAtomDump {
    format: String,
    version: usize,
    root: String,
    aliases: Vec<AliasedAtomDumpEntry>,
}

#[derive(Debug, Deserialize)]
struct AliasedAtomDumpEntry {
    alias: String,
    original: String,
}

fn main() -> Result<()> {
    color_eyre::install()?;

    let config = Config::from_args()?;

    let read_started = Instant::now();
    let bytes = fs::read(&config.input)
        .with_context(|| format!("failed to read {}", config.input.display()))?;
    println!(
        "read input: bytes={} elapsed={:?}",
        bytes.len(),
        read_started.elapsed()
    );

    let json_started = Instant::now();
    let dump: AliasedAtomDump = serde_json::from_slice(&bytes)
        .with_context(|| format!("failed to decode {}", config.input.display()))?;
    if dump.format != "gammaloop-aliased-atom" {
        bail!("unexpected dump format {:?}", dump.format);
    }
    if dump.version != 1 {
        bail!("unsupported aliased atom dump version {}", dump.version);
    }
    println!(
        "decoded json: aliases={} root_chars={} elapsed={:?}",
        dump.aliases.len(),
        dump.root.len(),
        json_started.elapsed()
    );

    let parse_started = Instant::now();
    let mut aliased = AliasedAtom::from(parse_atom("root", &dump.root)?);
    for (index, entry) in dump.aliases.into_iter().enumerate() {
        let alias = parse_atom(&format!("alias[{index}]"), &entry.alias)?;
        let original = parse_atom(&format!("original[{index}]"), &entry.original)?;
        aliased = aliased.add_alias(alias, original);
    }
    let root_terms = aliased.get_root().nterms();
    let root_bytes = aliased.get_root().as_view().get_byte_size();
    let (adds, muls) = aliased.count_operations();
    println!(
        "parsed aliased atom: aliases={} root_terms={} root_bytes={} adds={} muls={} elapsed={:?}",
        aliased.get_aliases().len(),
        root_terms,
        root_bytes,
        adds,
        muls,
        parse_started.elapsed()
    );

    println!("fill aliases: start");
    io::stdout().flush()?;
    let fill_started = Instant::now();
    let resolved = aliased.into_inner();
    println!(
        "fill aliases: done elapsed={:?} terms={} bytes={}",
        fill_started.elapsed(),
        resolved.nterms(),
        resolved.as_view().get_byte_size()
    );

    if let Some(output) = config.output {
        let write_started = Instant::now();
        write_atom(&output, &resolved)?;
        println!(
            "wrote output: path={} elapsed={:?}",
            output.display(),
            write_started.elapsed()
        );
    }

    Ok(())
}

impl Config {
    fn from_args() -> Result<Self> {
        let mut args = env::args().skip(1);
        let mut input = None;
        let mut output = None;

        while let Some(arg) = args.next() {
            match arg.as_str() {
                "-h" | "--help" => {
                    print_usage();
                    std::process::exit(0);
                }
                "-o" | "--output" => {
                    let Some(path) = args.next() else {
                        bail!("{arg} requires a path");
                    };
                    output = Some(PathBuf::from(path));
                }
                _ if input.is_none() => input = Some(PathBuf::from(arg)),
                _ => bail!("unexpected argument {arg:?}"),
            }
        }

        let Some(input) = input else {
            print_usage();
            bail!("missing aliased atom JSON path");
        };

        Ok(Self { input, output })
    }
}

fn parse_atom(label: &str, text: &str) -> Result<Atom> {
    Atom::parse_with_default_namespace(wrap_input!(text), ParseSettings::default())
        .map_err(|error| eyre!("failed to parse {label}: {error}"))
}

fn write_atom(path: &Path, atom: &Atom) -> Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    fs::write(path, atom.to_plain_string())
        .with_context(|| format!("failed to write {}", path.display()))
}

fn print_usage() {
    println!(
        "usage: rust-script scripts/fill_aliased_atom.rs ALIASED_ATOM.json [--output RESOLVED.sym]"
    );
}
