use std::{collections::HashMap, path::Path};

use clap::{Args, Parser, Subcommand};
use color_eyre::Result;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use gammalooprs::settings::{GlobalSettings, RuntimeSettings};
use thiserror::Error;

use super::state::State;

#[derive(Debug, Parser, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
/// Generate integrands
pub struct Generate {
    #[command(subcommand)]
    pub mode: Option<GenerateCmd>,
}

#[derive(Debug, Subcommand, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
/// Generate integrands
pub enum GenerateCmd {
    /// Cross-section (forward-scattering) generation
    Xs(SpecArgs),

    /// Amplitude generation
    Amp(SpecArgs),

    /// Reuse an already generated process
    Existing(ProcessArgs),
}

#[derive(Args, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SpecArgs {
    /// Quote-free spec:
    /// INITIAL to FINAL [veto ...] [only ...] [amp ...] [xs ...] [pert ...] [sets ...] [empty-initial] [empty-final]
    ///
    /// Example:
    ///   e+ e- to d d~ g veto u d xs QED^2==2 QCD^2<=4 pert loops=1 fsloops=2 QCD=2
    #[arg(value_name = "TOKENS", num_args = 1..)]
    pub tokens: Vec<String>,

    /// Number of threads (example common option)
    #[arg(short = 'n', long = "num-threads", default_value_t = 1)]
    pub num_threads: u32,
}

#[derive(Args, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct ProcessArgs {
    /// Stored process ID
    pub process_id: usize,

    /// Optional human name
    pub name: Option<String>,
}

impl Generate {
    pub fn run(
        &self,
        state: &mut State,
        compile_folder: impl AsRef<Path>,
        override_existing_compiled: bool,
        generation_settings: &GlobalSettings,
        runtime_settings: &RuntimeSettings,
    ) -> Result<()> {
        match &self.mode {
            Some(GenerateCmd::Xs(args)) => {
                let spec = must(parse_spec(&args.tokens));
                // hand off: mode = XS
                println!("MODE: XS\nSPEC: {spec:#?}\nthreads: {}", args.num_threads);
            }
            Some(GenerateCmd::Amp(args)) => {
                let spec = must(parse_spec(&args.tokens));
                // hand off: mode = AMP
                println!("MODE: AMP\nSPEC: {spec:#?}\nthreads: {}", args.num_threads);
            }
            Some(GenerateCmd::Existing(process)) => {
                return state.generate_integrand(
                    generation_settings,
                    runtime_settings.into(),
                    process,
                )
            }
            None => {
                state.generate_integrands(generation_settings, runtime_settings.into())?;

                if generation_settings.generation.evaluator_settings.compile {
                    state.compile_integrands(
                        compile_folder,
                        override_existing_compiled,
                        generation_settings,
                    )?
                }
                return Ok(());
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct ProcessSpec {
    pub initial: Vec<String>,
    pub final_: Vec<String>,
    pub empty_initial: bool,
    pub empty_final: bool,

    pub veto: Vec<String>, // replaces `/ ...`
    pub only: Vec<String>, // replaces `| ...`

    pub amp_couplings: Vec<String>, // QED==2, QCD<=4 (AMP)
    pub xs_couplings: Vec<String>,  // QED^2==2, QCD^2<=4 (XS)

    pub pert: Perturbative,           // {n}, {{m}}, QCD=2, QED=1
    pub final_sets: Vec<Vec<String>>, // XS multi-allowed finals
}

#[derive(Debug, Clone, Default)]
pub struct Perturbative {
    pub loops: Option<u32>,           // {n}
    pub fsloops: Option<u32>,         // {{n}}
    pub orders: HashMap<String, u32>, // QCD=2, QED=1; shorthand “QCD”≡1
}

// ----------------------- Spec tokenizer ------------------------

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("expected 'to' between initial and final states")]
    MissingTo,
    #[error("unknown clause '{0}'")]
    UnknownClause(String),
    #[error("bad value '{0}'")]
    BadValue(String),
}

pub fn parse_spec(tokens: &[String]) -> Result<ProcessSpec, ParseError> {
    use std::iter::Peekable;
    let mut it = tokens.iter().peekable();

    let mut spec = ProcessSpec {
        initial: Vec::new(),
        final_: Vec::new(),
        empty_initial: false,
        empty_final: false,
        veto: Vec::new(),
        only: Vec::new(),
        amp_couplings: Vec::new(),
        xs_couplings: Vec::new(),
        pert: Perturbative::default(),
        final_sets: Vec::new(),
    };

    // initial (or empty-initial) … to …
    if matches!(it.peek().map(|s| s.as_str()), Some("empty-initial")) {
        spec.empty_initial = true;
        it.next();
    } else {
        while let Some(p) = it.peek() {
            if p.as_str() == "to" {
                break;
            }
            spec.initial.push(it.next().unwrap().to_string());
        }
    }
    if !matches!(it.next().map(|s| s.as_str()), Some("to")) {
        return Err(ParseError::MissingTo);
    }

    // final (or empty-final)
    if matches!(it.peek().map(|s| s.as_str()), Some("empty-final")) {
        spec.empty_final = true;
        it.next();
    } else {
        while let Some(p) = it.peek() {
            let s = p.as_str();
            if is_clause(s) {
                break;
            }
            spec.final_.push(it.next().unwrap().to_string());
        }
    }

    // clauses
    while let Some(kw) = it.next() {
        match kw.as_str() {
            "veto" => collect_until_clause(&mut it, &mut spec.veto),
            "only" => collect_until_clause(&mut it, &mut spec.only),
            "amp" => collect_until_clause(&mut it, &mut spec.amp_couplings),
            "xs" => collect_until_clause(&mut it, &mut spec.xs_couplings),
            "pert" => {
                while let Some(p) = it.peek() {
                    let s = p.as_str();
                    if is_clause(s) {
                        break;
                    }
                    let term = it.next().unwrap();
                    if let Some((k, v)) = term.split_once('=') {
                        match k {
                            "loops" => {
                                spec.pert.loops = Some(
                                    v.parse().map_err(|_| ParseError::BadValue(term.clone()))?,
                                )
                            }
                            "fsloops" => {
                                spec.pert.fsloops = Some(
                                    v.parse().map_err(|_| ParseError::BadValue(term.clone()))?,
                                )
                            }
                            other => {
                                let n =
                                    v.parse().map_err(|_| ParseError::BadValue(term.clone()))?;
                                spec.pert.orders.insert(other.to_string(), n);
                            }
                        }
                    } else {
                        // shorthand: QCD ≡ QCD=1
                        spec.pert.orders.insert(term.to_string(), 1);
                    }
                }
            }
            "sets" => {
                loop {
                    let mut set = Vec::new();
                    while let Some(p) = it.peek() {
                        let s = p.as_str();
                        if s == ";" {
                            it.next();
                            break;
                        }
                        if is_clause(s) {
                            break;
                        }
                        set.push(it.next().unwrap().to_string());
                    }
                    if !set.is_empty() {
                        spec.final_sets.push(set);
                    }
                    match it.peek().map(|s| s.as_str()) {
                        Some(";") => {
                            it.next();
                            continue;
                        }
                        Some(s) if is_clause(s) => break,
                        None => break,
                        _ => {} // continue collecting another set
                    }
                }
            }
            "empty-final" => {
                spec.empty_final = true;
            }
            other => return Err(ParseError::UnknownClause(other.to_string())),
        }
    }

    Ok(spec)
}

fn collect_until_clause<'a, I: Iterator<Item = &'a String>>(
    it: &mut std::iter::Peekable<I>,
    out: &mut Vec<String>,
) {
    while let Some(p) = it.peek() {
        let s = p.as_str();
        if is_clause(s) {
            break;
        }
        out.push(it.next().unwrap().to_string());
    }
}
fn is_clause(s: &str) -> bool {
    matches!(
        s,
        "veto" | "only" | "amp" | "xs" | "pert" | "sets" | "empty-final" | "to"
    )
}

fn must<T, E: std::fmt::Display>(r: Result<T, E>) -> T {
    match r {
        Ok(v) => v,
        Err(e) => {
            eprintln!("parse error: {e}");
            std::process::exit(2);
        }
    }
}
