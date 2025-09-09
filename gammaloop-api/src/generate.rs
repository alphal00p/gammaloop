#![allow(clippy::too_many_arguments)]

use ahash::HashMap;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use clap::{Args, Parser, Subcommand, ValueEnum};
use color_eyre::Result;
use regex::Regex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use gammalooprs::feyngen::{
    diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions, GenerationType,
    NumeratorAwareGraphGroupingOption, SelfEnergyFilterOptions, SewedFilterOptions,
    SnailFilterOptions, TadpolesFilterOptions,
};
use gammalooprs::model::Model;
use gammalooprs::numerator::GlobalPrefactor;
use gammalooprs::settings::{GlobalSettings, RuntimeSettings};

use super::state::State;

// =================== CLI containers (kept close to your structure) ===================

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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq, ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum GroupingChoice {
    #[clap(name = "no_grouping")]
    NoGrouping,
    #[clap(name = "only_detect_zeroes")]
    OnlyDetectZeroes,
    #[clap(name = "group_identical_graphs_up_to_sign")]
    GroupIdenticalGraphsUpToSign,
    #[clap(name = "group_identical_graphs_up_to_scalar_rescaling")]
    GroupIdenticalGraphsUpToScalarRescaling,
}

impl GroupingChoice {
    fn as_strategy_str(self) -> &'static str {
        match self {
            GroupingChoice::NoGrouping => "no_grouping",
            GroupingChoice::OnlyDetectZeroes => "only_detect_zeroes",
            GroupingChoice::GroupIdenticalGraphsUpToSign => "group_identical_graphs_up_to_sign",
            GroupingChoice::GroupIdenticalGraphsUpToScalarRescaling => {
                "group_identical_graphs_up_to_scalar_rescaling"
            }
        }
    }
}

#[derive(Args, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct SpecArgs {
    /// Quote-free spec.
    ///
    /// Grammar (unquoted): `INITIAL to FINAL [process_options] [generation_options]`
    /// Also accepts: `INITIAL > FINAL [ ... ]`
    ///
    /// Examples:
    ///   e+ e- to d d~ g veto u d xs QED^2==2 QCD^2<=4 pert loops=1 fsloops=2 QCD=2
    ///   e+ e- > { z z, a a } [ {1} {{2}} QCD=2 QED=1 ] / u d | g ghG QED==2 QCD>=2 QCD<=4
    ///
    /// Notes:
    /// - `veto` and `only` are synonyms for `/` and `|`.
    /// - The `[ ... ]` block holds perturbative spec: `{n}`, `{{n}}`, and orders `QCD=n` with shorthand `QCD`≡1.
    #[arg(value_name = "TOKENS", num_args = 1..)]
    pub tokens: Vec<String>,

    // --------- Generation options that influence FeynGenOptions ---------
    /// Number of threads
    #[arg(short = 'n', long = "num-threads", default_value_t = 1)]
    pub num_threads: u32,

    /// Clear pre-existing processes
    #[arg(
        long = "clear-existing-processes",
        alias = "clear",
        default_value_t = false
    )]
    pub clear_existing_processes: bool,

    /// Topology filters (Option<bool> to enable smart defaults)
    #[arg(long = "filter-selfenergies")]
    pub filter_selfenergies: Option<bool>,
    #[arg(long = "filter-snails")]
    pub filter_snails: Option<bool>,
    #[arg(long = "filter-tadpoles")]
    pub filter_tadpoles: Option<bool>,

    /// Extra filters / constraints
    #[arg(long = "max-n-bridges")]
    pub max_n_bridges: Option<i32>,
    #[arg(long = "number-of-factorized-loop-subtopologies", alias = "nfactl")]
    pub number_of_factorized_loop_subtopologies: Option<i32>,
    /// Number of closed loops of anticommutating fields; -1 disables
    #[arg(long = "number-of-anticommutating-loops")]
    pub number_of_anticommutating_loops: Option<i32>,

    /// Symmetrization
    #[arg(
        long = "allow-symmetrization-of-external-fermions-in-amplitudes",
        alias = "symferm",
        default_value_t = false
    )]
    pub allow_symmetrization_of_external_fermions_in_amplitudes: bool,
    #[arg(long = "symmetrize-initial-states", default_value_t = false)]
    pub symmetrize_initial_states: bool,
    #[arg(long = "symmetrize-final-states", default_value_t = false)]
    pub symmetrize_final_states: bool,
    #[arg(long = "symmetrize-left-right-states", default_value_t = false)]
    pub symmetrize_left_right_states: bool,

    /// Numerator-aware grouping choice (restricted)
    #[arg(long = "numerator-aware-isomorphism-grouping", value_enum)]
    pub numerator_aware_isomorphism_grouping: Option<GroupingChoice>,

    /// Grouping attributes
    #[arg(long = "numerical-samples-seed")]
    pub numerical_samples_seed: Option<u16>,
    #[arg(long = "number-of-samples-for-numerator-comparisons")]
    pub number_of_samples_for_numerator_comparisons: Option<usize>,
    #[arg(
        long = "consider-internal-masses-only-in-numerator-isomorphisms",
        default_value_t = false
    )]
    pub consider_internal_masses_only_in_numerator_isomorphisms: bool,
    #[arg(
        long = "fully-numerical-substitution-when-comparing-numerators",
        default_value_t = true
    )]
    pub fully_numerical_substitution_when_comparing_numerators: bool,
    #[arg(long = "compare-canonized-numerator", default_value_t = false)]
    pub compare_canonized_numerator: bool,

    /// Graph processing toggles
    ///
    /// Format:
    ///   --loop-momentum-bases "B1=p1,p2;B2=q1,q2"
    ///   --select-graphs "g1,g2,g3"
    ///   --veto-graphs "h1,h2"
    #[arg(long = "loop-momentum-bases")]
    pub loop_momentum_bases: Option<String>,
    #[arg(long = "select-graphs")]
    pub select_graphs: Option<String>,
    #[arg(long = "veto-graphs")]
    pub veto_graphs: Option<String>,
    #[arg(long = "graph-prefix")]
    pub graph_prefix: Option<String>,

    /// Filter self-loop (explicit; default false)
    #[arg(long = "filter-self-loop")]
    pub filter_self_loop: Option<bool>,
}

#[derive(Args, Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct ProcessArgs {
    /// Stored process ID
    pub process_id: usize,

    /// Optional human name
    pub name: Option<String>,
}

// =================== Domain structs ===================

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq, Default)]
pub struct OrderRange {
    pub eq: Option<u32>,
    pub min: Option<u32>,
    pub max: Option<u32>,
}
impl OrderRange {
    fn set(&mut self, op: &str, v: u32) {
        match op {
            "==" => self.eq = Some(v),
            ">=" => self.min = Some(v),
            "<=" => self.max = Some(v),
            _ => {}
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq, PartialOrd, Ord)]
pub struct CouplingKey {
    pub name: String,
    pub power: u32, // parsed, but dropped when building FeynGenFilter::CouplingOrders
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct Perturbative {
    /// {n}
    pub loops_sum_amp_or_sum: Option<u32>,
    /// {{n}}
    pub loops_forward_graph: Option<u32>,
    /// QCD=2, QED=1; shorthand “QCD”≡1
    pub orders: BTreeMap<String, u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct ProcessSpec {
    pub initial: Vec<String>,
    pub final_: Vec<String>,
    pub empty_initial: bool,
    pub empty_final: bool,

    /// `/ ...`
    pub veto: BTreeSet<String>,
    /// `| ...`
    pub only: Option<BTreeSet<String>>,

    /// AMP amplitude coupling ranges: QED==2, QCD>=2, QCD<=4
    pub amp_couplings: BTreeMap<String, OrderRange>,
    /// XS graph-level coupling ranges: QED^2==2, QCD^2<=4  (power is dropped when building filters)
    pub xs_couplings: BTreeMap<CouplingKey, OrderRange>,

    /// [{1} {{2}} QCD=2 QED=1] or [QCD]
    pub pert: Perturbative,

    /// XS-only: final-state alternative sets `{ A, B, ... }`
    pub final_sets: Vec<Vec<String>>,

    /// Generation options captured alongside the spec
    pub feyngen: FeynGenOptions,

    /// Parsed numerator-aware grouping mode (not part of FeynGenOptions)
    pub numerator_grouping: Option<NumeratorAwareGraphGroupingOption>,
}

// =================== Runner ===================

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
                let model: &Model = &state.model;
                let spec = must(parse_spec_with_model(
                    args, /*amplitude=*/ false, model,
                ));
                // Try a generation call; report count.
                let grouping = spec
                    .numerator_grouping
                    .clone()
                    .unwrap_or(NumeratorAwareGraphGroupingOption::NoGrouping);
                let filter_self_loop = args.filter_self_loop.unwrap_or(false);
                let graph_prefix = args.graph_prefix.clone().unwrap_or_default();
                let selected_graphs = args.select_graphs.as_ref().map(|s| parse_csv_list(s));
                let vetoed_graphs = args.veto_graphs.as_ref().map(|s| parse_csv_list(s));
                let loop_momentum_bases = args
                    .loop_momentum_bases
                    .as_deref()
                    .and_then(parse_loop_momentum_bases);

                let pref = GlobalPrefactor::default();

                println!("MODE: XS\nSPEC: {:#?}", spec.clone());
                let feyngen = FeynGen::new(spec.feyngen);

                let graphs = FeynGen::generate(
                    &feyngen,
                    model,
                    &grouping,
                    filter_self_loop,
                    graph_prefix,
                    selected_graphs,
                    vetoed_graphs,
                    loop_momentum_bases,
                    pref,
                    Some(args.num_threads as usize),
                )?;
                println!("Generated {} graphs.", graphs.len());
            }
            Some(GenerateCmd::Amp(args)) => {
                let model: &Model = &state.model;
                let spec = must(parse_spec_with_model(args, /*amplitude=*/ true, model));
                let grouping = spec
                    .numerator_grouping
                    .clone()
                    .unwrap_or(NumeratorAwareGraphGroupingOption::NoGrouping);
                let filter_self_loop = args.filter_self_loop.unwrap_or(false);
                let graph_prefix = args.graph_prefix.clone().unwrap_or_default();
                let selected_graphs = args.select_graphs.as_ref().map(|s| parse_csv_list(s));
                let vetoed_graphs = args.veto_graphs.as_ref().map(|s| parse_csv_list(s));
                let loop_momentum_bases = args
                    .loop_momentum_bases
                    .as_deref()
                    .and_then(parse_loop_momentum_bases);
                let pref = GlobalPrefactor::default();

                println!("MODE: AMP\nSPEC: {:#?}", spec.clone());
                let feyngen = FeynGen::new(spec.feyngen);

                let graphs = FeynGen::generate(
                    &feyngen,
                    model,
                    &grouping,
                    filter_self_loop,
                    graph_prefix,
                    selected_graphs,
                    vetoed_graphs,
                    loop_momentum_bases,
                    pref,
                    Some(args.num_threads as usize),
                )?;
                println!("Generated {} graphs.", graphs.len());
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

// =================== Parsing ===================

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParseError {
    #[error("expected 'to' or '>' between initial and final states")]
    MissingArrow,
    #[error("unbalanced brackets/braces")]
    Unbalanced,
    #[error("invalid token '{0}'")]
    InvalidToken(String),
    #[error("unknown particle '{name}'. Valid choices: {choices}")]
    UnknownParticle { name: String, choices: String },
    #[error("unknown coupling '{name}'. Valid choices: {choices}")]
    UnknownCoupling { name: String, choices: String },
}

fn must<T, E: std::fmt::Display>(r: std::result::Result<T, E>) -> T {
    match r {
        Ok(v) => v,
        Err(e) => {
            eprintln!("parse error: {e}");
            std::process::exit(2);
        }
    }
}

pub fn parse_spec_with_model(
    args: &SpecArgs,
    amplitude: bool,
    model: &Model,
) -> std::result::Result<ProcessSpec, ParseError> {
    let raw = args.tokens.join(" ");
    let (lhs, rhs0) = split_top_level_arrow(&raw).ok_or(ParseError::MissingArrow)?;
    let (rhs_wo_bracket, pert) = extract_perturbative(rhs0.trim())?;

    // LHS initial states
    let lhs = lhs.trim();
    let (initial_names, empty_initial) = if lhs == "{}" || lhs.eq_ignore_ascii_case("empty-initial")
    {
        (Vec::new(), true)
    } else {
        (split_ws(lhs), false)
    };

    // Tokenize RHS outside the perturbative block
    let tokens = tokenize(rhs_wo_bracket)?;
    let (final_spec, empty_final, rest) = parse_final_states(tokens)?;
    let (veto_names, only, amp_couplings, xs_couplings) = parse_process_options(rest)?;

    // Validate coupling orders against model
    validate_coupling_names(model, &amp_couplings, &xs_couplings, &pert)?;

    // Convert names → PDGs using the model
    let initial_pdgs = resolve_pdgs(model, &initial_names)?;
    let primary_final_pdgs = resolve_pdgs(model, &final_spec.primary)?;
    let final_sets_pdgs = final_spec
        .sets
        .iter()
        .map(|set| resolve_pdgs(model, set))
        .collect::<Result<Vec<_>, _>>()?;

    // Resolve vetoed particles to PDGs and validate they exist
    let veto_pdgs = resolve_pdgs(model, &veto_names.iter().cloned().collect::<Vec<_>>())?;

    // Build FeynGenOptions with filters and PDGs
    let numerator_grouping = build_grouping_option(args);
    let feyngen = feyngen_from_spec_args(
        args,
        amplitude,
        &pert,
        &amp_couplings,
        &xs_couplings,
        &initial_pdgs,
        &primary_final_pdgs,
        &final_sets_pdgs,
        &veto_pdgs,
    );

    let spec = ProcessSpec {
        initial: initial_names,
        final_: final_spec.primary.clone(),
        empty_initial,
        empty_final,
        veto: veto_names,
        only,
        amp_couplings,
        xs_couplings,
        pert,
        final_sets: final_spec.sets,
        feyngen,
        numerator_grouping,
    };
    Ok(spec)
}

// ---- helpers ----

fn split_top_level_arrow(s: &str) -> Option<(&str, &str)> {
    // Accept either "to" token or '>' outside any [...] or {...}
    let mut depth_brace = 0i32;
    let mut depth_bracket = 0i32;
    let bytes = s.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        let c = bytes[i] as char;
        match c {
            '{' => depth_brace += 1,
            '}' => depth_brace -= 1,
            '[' => depth_bracket += 1,
            ']' => depth_bracket -= 1,
            '>' if depth_brace == 0 && depth_bracket == 0 => {
                return Some((&s[..i], &s[i + 1..]));
            }
            _ => {
                if depth_brace == 0
                    && depth_bracket == 0
                    && i + 2 < bytes.len()
                    && (c == 't' || c == 'T')
                    && &s[i..].to_ascii_lowercase().as_str()[..2] == "to"
                {
                    let before = if i == 0 {
                        ' '
                    } else {
                        s.as_bytes()[i - 1] as char
                    };
                    let after = if i + 2 >= s.len() {
                        ' '
                    } else {
                        s.as_bytes()[i + 2] as char
                    };
                    if before.is_whitespace() && after.is_whitespace() {
                        return Some((&s[..i], &s[i + 2..]));
                    }
                }
            }
        }
        i += 1;
    }
    None
}

fn extract_perturbative(rhs: &str) -> Result<(String, Perturbative), ParseError> {
    let mut out = rhs.to_string();
    let mut spec = Perturbative::default();
    if let Some((start, end)) = find_top_level_block(rhs, '[', ']') {
        let inside = &rhs[start + 1..end];
        spec = parse_perturbative_block(inside.trim())?;
        out.replace_range(start..=end, "");
    }
    Ok((out, spec))
}

fn parse_perturbative_block(s: &str) -> Result<Perturbative, ParseError> {
    if s.is_empty() {
        return Ok(Perturbative::default());
    }
    let re_amp = Regex::new(r"^\{(\d+)\}$").unwrap();
    let re_fwd = Regex::new(r"^\{\{(\d+)\}\}$").unwrap();
    let re_order = Regex::new(r"^([A-Za-z_][A-Za-z0-9_]*)=(\d+)$").unwrap();

    let mut spec = Perturbative::default();
    for tok in split_ws(s) {
        if let Some(c) = re_amp.captures(&tok) {
            spec.loops_sum_amp_or_sum = Some(c[1].parse().unwrap());
            continue;
        }
        if let Some(c) = re_fwd.captures(&tok) {
            spec.loops_forward_graph = Some(c[1].parse().unwrap());
            continue;
        }
        if let Some(c) = re_order.captures(&tok) {
            let n: u32 = c[2].parse().unwrap();
            spec.orders.insert(c[1].to_string(), n);
            continue;
        }
        if tok.chars().all(|ch| ch.is_ascii_alphabetic()) {
            spec.orders.insert(tok, 1);
            continue;
        }
        return Err(ParseError::InvalidToken(tok));
    }
    Ok(spec)
}

fn tokenize(s: String) -> Result<Vec<String>, ParseError> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    for ch in s.chars() {
        match ch {
            '{' => {
                depth += 1;
                cur.push(ch);
            }
            '}' => {
                depth -= 1;
                cur.push(ch);
                if depth < 0 {
                    return Err(ParseError::Unbalanced);
                }
            }
            c if c.is_whitespace() && depth == 0 => {
                if !cur.trim().is_empty() {
                    out.push(cur.trim().to_string());
                }
                cur.clear();
            }
            _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_string());
    }
    if depth != 0 {
        return Err(ParseError::Unbalanced);
    }
    Ok(out)
}

fn parse_final_states(tokens: Vec<String>) -> Result<(FinalSpec, bool, Vec<String>), ParseError> {
    if tokens.is_empty() {
        return Ok((FinalSpec::default(), false, vec![]));
    }
    let mut t = tokens;
    let mut finals = Vec::<String>::new();
    let mut empty_final = false;

    if t[0].eq_ignore_ascii_case("empty-final") || t[0] == "{}" {
        empty_final = true;
        t.remove(0);
        return Ok((FinalSpec::default(), empty_final, t));
    }

    if t[0].starts_with('{') {
        let raw = t.remove(0);
        let sets = parse_alt_sets(&raw)?;
        return Ok((
            FinalSpec {
                primary: vec![],
                sets,
            },
            empty_final,
            t,
        ));
    }

    while !t.is_empty() {
        if is_option_token(&t[0]) {
            break;
        }
        finals.push(t.remove(0));
    }
    Ok((
        FinalSpec {
            primary: finals,
            sets: vec![],
        },
        empty_final,
        t,
    ))
}

#[derive(Default)]
struct FinalSpec {
    primary: Vec<String>,
    sets: Vec<Vec<String>>,
}

fn parse_alt_sets(s: &str) -> Result<Vec<Vec<String>>, ParseError> {
    let trimmed = s.trim();
    if !(trimmed.starts_with('{') && trimmed.ends_with('}')) {
        return Err(ParseError::InvalidToken(s.to_string()));
    }
    let inner = &trimmed[1..trimmed.len() - 1];
    let mut out: Vec<Vec<String>> = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    for ch in inner.chars() {
        match ch {
            '{' => {
                depth += 1;
                cur.push(ch);
            }
            '}' => {
                depth -= 1;
                if depth < 0 {
                    return Err(ParseError::Unbalanced);
                }
                cur.push(ch);
            }
            ',' if depth == 0 => {
                let items = split_ws(cur.trim());
                if !items.is_empty() {
                    out.push(items);
                } else {
                    out.push(vec![]);
                }
                cur.clear();
            }
            _ => cur.push(ch),
        }
    }
    if !cur.trim().is_empty() {
        out.push(split_ws(cur.trim()));
    }
    Ok(out)
}

fn is_option_token(tok: &str) -> bool {
    matches!(
        tok,
        "/" | "|" | "veto" | "only" | "amp" | "xs" | "pert" | "sets" | "empty-final"
    ) || tok.starts_with('[')
        || Regex::new(r"^[A-Za-z_][A-Za-z0-9_]*(?:\^\d+)?(==|>=|<=)\d+$")
            .unwrap()
            .is_match(tok)
}

#[allow(clippy::type_complexity)]
fn parse_process_options(
    tokens: Vec<String>,
) -> Result<
    (
        BTreeSet<String>,
        Option<BTreeSet<String>>,
        BTreeMap<String, OrderRange>,
        BTreeMap<CouplingKey, OrderRange>,
    ),
    ParseError,
> {
    let mut veto = BTreeSet::<String>::new();
    let mut only: Option<BTreeSet<String>> = None;
    let mut amp = BTreeMap::<String, OrderRange>::new();
    let mut xs = BTreeMap::<CouplingKey, OrderRange>::new();

    let re = Regex::new(r"^([A-Za-z_][A-Za-z0-9_]*)(?:\^(\d+))?(==|>=|<=)(\d+)$").unwrap();

    let mut i = 0usize;
    while i < tokens.len() {
        match tokens[i].as_str() {
            "/" | "veto" => {
                i += 1;
                while i < tokens.len() && !is_option_token(&tokens[i]) {
                    veto.insert(tokens[i].clone());
                    i += 1;
                }
            }
            "|" | "only" => {
                i += 1;
                let mut s = only.take().unwrap_or_default();
                while i < tokens.len() && !is_option_token(&tokens[i]) {
                    s.insert(tokens[i].clone());
                    i += 1;
                }
                only = Some(s);
            }
            "amp" | "xs" | "pert" | "sets" | "empty-final" => {
                return Err(ParseError::InvalidToken(tokens[i].clone()));
            }
            other => {
                if let Some(c) = re.captures(other) {
                    let name = c[1].to_string();
                    let power: u32 = c.get(2).map(|m| m.as_str().parse().unwrap()).unwrap_or(1);
                    let op = c.get(3).unwrap().as_str();
                    let val: u32 = c.get(4).unwrap().as_str().parse().unwrap();
                    if power == 1 {
                        let e = amp.entry(name).or_default();
                        e.set(op, val);
                    } else {
                        let key = CouplingKey { name, power };
                        let e = xs.entry(key).or_default();
                        e.set(op, val);
                    }
                    i += 1;
                } else {
                    return Err(ParseError::InvalidToken(other.to_string()));
                }
            }
        }
    }
    Ok((veto, only, amp, xs))
}

fn find_top_level_block(s: &str, open: char, close: char) -> Option<(usize, usize)> {
    let mut depth = 0i32;
    let mut start: Option<usize> = None;
    for (i, ch) in s.char_indices() {
        if ch == open {
            if depth == 0 {
                start = Some(i);
            }
            depth += 1;
        } else if ch == close {
            depth -= 1;
            if depth == 0 {
                return start.map(|st| (st, i));
            }
        }
    }
    None
}

fn split_ws(s: &str) -> Vec<String> {
    s.split_whitespace().map(|t| t.to_string()).collect()
}

// =================== Build Feyngen options and filters ===================

fn feyngen_from_spec_args(
    a: &SpecArgs,
    amplitude: bool,
    pert: &Perturbative,
    amp_couplings: &BTreeMap<String, OrderRange>,
    xs_couplings: &BTreeMap<CouplingKey, OrderRange>,
    initial_pdgs: &[i64],
    primary_final_pdgs: &[i64],
    final_sets_pdgs: &[Vec<i64>],
    veto_pdgs: &[i64],
) -> FeynGenOptions {
    let generation_type = if amplitude {
        GenerationType::Amplitude
    } else {
        GenerationType::CrossSection
    };

    // Decide vacuum-like topology from PDGs
    let is_vacuum = if initial_pdgs.is_empty() {
        match generation_type {
            GenerationType::CrossSection => true,
            GenerationType::Amplitude => {
                primary_final_pdgs.is_empty()
                    && (final_sets_pdgs
                        .first()
                        .map(|v| v.is_empty())
                        .unwrap_or(true))
            }
        }
    } else {
        false
    };

    // Smart defaults influenced by vacuum
    let filter_tadpoles_default = !is_vacuum;
    let filter_snails_default = !is_vacuum;
    let filter_selfenergies_default = !is_vacuum;

    // Normalize filter toggles
    let filter_tadpoles = a.filter_tadpoles.unwrap_or(filter_tadpoles_default);
    let filter_snails = a.filter_snails.unwrap_or(filter_snails_default);
    let filter_selfenergies = a.filter_selfenergies.unwrap_or(filter_selfenergies_default);

    // Normalize numeric toggles with the vacuum defaults
    let number_of_factorized_loop_subtopologies = a
        .number_of_factorized_loop_subtopologies
        .or(if is_vacuum { Some(1) } else { None })
        .and_then(|n| if n < 0 { None } else { Some(n as usize) });

    let max_n_bridges = a
        .max_n_bridges
        .or(if is_vacuum { Some(0) } else { None })
        .and_then(|n| if n < 0 { None } else { Some(n as usize) });

    let number_of_fermion_loops =
        a.number_of_anticommutating_loops
            .and_then(|n| if n < 0 { None } else { Some(n as usize) });

    // Base options
    let mut fg = FeynGenOptions {
        generation_type,
        initial_pdgs: initial_pdgs.to_vec(),
        final_pdgs_lists: if final_sets_pdgs.is_empty() {
            vec![primary_final_pdgs.to_vec()]
        } else if !primary_final_pdgs.is_empty() {
            let mut sets = vec![primary_final_pdgs.to_vec()];
            sets.extend_from_slice(final_sets_pdgs);
            sets
        } else {
            final_sets_pdgs.to_vec()
        },
        loop_count_range: (1, 1), // may be overridden below
        symmetrize_initial_states: a.symmetrize_initial_states,
        symmetrize_final_states: a.symmetrize_final_states,
        symmetrize_left_right_states: a.symmetrize_left_right_states,
        allow_symmetrization_of_external_fermions_in_amplitudes: a
            .allow_symmetrization_of_external_fermions_in_amplitudes,
        max_multiplicity_for_fast_cut_filter: 6,
        amplitude_filters: FeynGenFilters(vec![]),
        cross_section_filters: FeynGenFilters(vec![]),
    };

    // Build filters
    let mut amp_filters: Vec<FeynGenFilter> = Vec::new();
    let mut xs_filters: Vec<FeynGenFilter> = Vec::new();

    // Loop counts
    if let Some(n) = pert.loops_sum_amp_or_sum {
        amp_filters.push(FeynGenFilter::LoopCountRange((n as usize, n as usize)));
        if fg.generation_type == GenerationType::Amplitude {
            fg.loop_count_range = (n as usize, n as usize);
        }
    }
    if let Some(n) = pert.loops_forward_graph {
        xs_filters.push(FeynGenFilter::LoopCountRange((n as usize, n as usize)));
        if fg.generation_type == GenerationType::CrossSection {
            fg.loop_count_range = (n as usize, n as usize);
        }
    }

    // Perturbative orders (only to matching container)
    if !pert.orders.is_empty() {
        let map: HashMap<String, usize> = pert
            .orders
            .iter()
            .map(|(k, v)| (k.clone(), *v as usize))
            .collect();
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(FeynGenFilter::PerturbativeOrders(map));
        } else {
            xs_filters.push(FeynGenFilter::PerturbativeOrders(map));
        }
    }

    // CouplingOrders
    if !amp_couplings.is_empty() {
        let cmap = coupling_orders_from_order_ranges_amp(amp_couplings);
        amp_filters.push(FeynGenFilter::CouplingOrders(cmap));
    }
    if !xs_couplings.is_empty() {
        let cmap = coupling_orders_from_order_ranges_xs(xs_couplings);
        xs_filters.push(FeynGenFilter::CouplingOrders(cmap));
    }

    // Factorized loop topologies → FactorizedLoopTopologiesCountRange((n,n))
    if let Some(n) = number_of_factorized_loop_subtopologies {
        let filt = FeynGenFilter::FactorizedLoopTopologiesCountRange((n, n));
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
        }
    }

    // Fermion loop count → FermionLoopCountRange((n,n))
    if let Some(n) = number_of_fermion_loops {
        let filt = FeynGenFilter::FermionLoopCountRange((n, n));
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
        }
    }

    // MaxNumberOfBridges
    if let Some(n) = max_n_bridges {
        let filt = FeynGenFilter::MaxNumberOfBridges(n);
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
        }
    }

    // Particle vetoes
    if !veto_pdgs.is_empty() {
        let filt = FeynGenFilter::ParticleVeto(veto_pdgs.to_vec());
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
        }
    }

    // Self-energy / snails / tadpoles:
    if filter_selfenergies {
        let filt = FeynGenFilter::SelfEnergyFilter(SelfEnergyFilterOptions::default());
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
        }
    }
    if filter_snails {
        let filt = FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions::default());
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
        }
    }
    if filter_tadpoles {
        let filt = FeynGenFilter::TadpolesFilter(TadpolesFilterOptions::default());
        if fg.generation_type == GenerationType::Amplitude {
            amp_filters.push(filt);
        } else {
            xs_filters.push(filt);
            // XS-specific sewed filter for tadpoles
            xs_filters.push(FeynGenFilter::SewedFilter(SewedFilterOptions {
                filter_tadpoles: true,
            }));
        }
    }

    fg.amplitude_filters = FeynGenFilters(amp_filters);
    fg.cross_section_filters = FeynGenFilters(xs_filters);

    fg
}

fn coupling_orders_from_order_ranges_amp(
    orders: &BTreeMap<String, OrderRange>,
) -> HashMap<String, (usize, Option<usize>)> {
    let mut out = HashMap::default();
    for (name, r) in orders {
        if let Some(eq) = r.eq {
            out.insert(name.clone(), (eq as usize, Some(eq as usize)));
        } else {
            let minv = r.min.unwrap_or(0) as usize;
            let maxv = r.max.map(|m| m as usize);
            out.insert(name.clone(), (minv, maxv));
        }
    }
    out
}

fn coupling_orders_from_order_ranges_xs(
    orders: &BTreeMap<CouplingKey, OrderRange>,
) -> HashMap<String, (usize, Option<usize>)> {
    // Drop the power and merge constraints per name.
    let mut merged: BTreeMap<String, OrderRange> = BTreeMap::new();
    for (k, r) in orders {
        let e = merged.entry(k.name.clone()).or_default();
        if let Some(eq) = r.eq {
            e.eq = Some(eq);
        }
        if let Some(min) = r.min {
            e.min = Some(e.min.map_or(min, |m| m.max(min)));
        }
        if let Some(max) = r.max {
            e.max = Some(e.max.map_or(max, |m| m.min(max)));
        }
    }
    coupling_orders_from_order_ranges_amp(&merged)
}

fn build_grouping_option(a: &SpecArgs) -> Option<NumeratorAwareGraphGroupingOption> {
    let choice = a.numerator_aware_isomorphism_grouping?;
    let strategy = choice.as_strategy_str();
    let opt = NumeratorAwareGraphGroupingOption::new_with_attributes(
        strategy,
        a.numerical_samples_seed,
        a.number_of_samples_for_numerator_comparisons,
        Some(a.consider_internal_masses_only_in_numerator_isomorphisms),
        Some(a.fully_numerical_substitution_when_comparing_numerators),
        Some(a.compare_canonized_numerator),
    );
    match opt {
        Ok(v) => Some(v),
        Err(_e) => Some(NumeratorAwareGraphGroupingOption::NoGrouping),
    }
}

// ---- Model-backed validation and resolution ----

fn all_particle_names(model: &Model) -> Vec<String> {
    model
        .particles
        .iter()
        .map(|p| p.name.clone().to_string())
        .collect::<Vec<_>>()
}

fn resolve_pdgs(model: &Model, names: &[String]) -> Result<Vec<i64>, ParseError> {
    let all_particles = all_particle_names(model);
    let mut out = Vec::with_capacity(names.len());
    'outer: for n in names {
        // Accept numeric PDG directly
        if let Ok(p) = n.parse::<i64>() {
            out.push(p);
            continue;
        }
        // Try exact name lookup via model
        if let Ok(p) = model.try_get_particle(n) {
            out.push(p.pdg_code as i64);
            continue;
        }
        // Case-insensitive fallback against both names and antinames
        for part in model.particles.iter() {
            if part.name.eq_ignore_ascii_case(n) || part.antiname.eq_ignore_ascii_case(n) {
                out.push(part.pdg_code as i64);
                continue 'outer;
            }
        }
        // Unknown
        return Err(ParseError::UnknownParticle {
            name: n.clone(),
            choices: all_particles.join(", "),
        });
    }
    Ok(out)
}

fn validate_coupling_names(
    model: &Model,
    amp: &BTreeMap<String, OrderRange>,
    xs: &BTreeMap<CouplingKey, OrderRange>,
    pert: &Perturbative,
) -> Result<(), ParseError> {
    // If the model has no coupling metadata, skip validation entirely.
    // (Some generic models omit coupling-order names.)
    let has_couplings_info = !model.couplings.is_empty();

    if !has_couplings_info {
        return Ok(());
    }

    // Collect all mentioned coupling names
    let mut names = BTreeSet::<String>::new();
    for k in amp.keys() {
        names.insert(k.clone());
    }
    for k in xs.keys() {
        names.insert(k.name.clone());
    }
    for k in pert.orders.keys() {
        names.insert(k.clone());
    }

    if names.is_empty() {
        return Ok(());
    }

    // Known couplings (case-insensitive compare)
    let known: Vec<String> = model
        .orders
        .iter()
        .map(|o| o.name.to_string())
        .collect::<Vec<_>>();
    let known_lower: BTreeSet<String> = known.iter().map(|s| s.to_ascii_lowercase()).collect();

    for n in names {
        if !known_lower.contains(&n.to_ascii_lowercase()) {
            return Err(ParseError::UnknownCoupling {
                name: n,
                choices: known.join(", "),
            });
        }
    }
    Ok(())
}

// ---- parsing of structured CLI fields ----

fn parse_csv_list(s: &str) -> Vec<String> {
    s.split(',')
        .flat_map(|x| x.split_whitespace())
        .filter(|t| !t.is_empty())
        .map(|t| t.to_string())
        .collect()
}

fn parse_loop_momentum_bases(s: &str) -> Option<HashMap<String, Vec<String>>> {
    // Format: "B1=p1,p2;B2=q1,q2"
    let mut out: HashMap<String, Vec<String>> = HashMap::default();
    for kv in s.split(';') {
        let kv = kv.trim();
        if kv.is_empty() {
            continue;
        }
        if let Some((k, v)) = kv.split_once('=') {
            let key = k.trim().to_string();
            let vals = v
                .split(',')
                .map(|x| x.trim().to_string())
                .filter(|x| !x.is_empty())
                .collect::<Vec<_>>();
            if !key.is_empty() {
                out.insert(key, vals);
            }
        } else {
            // single name without '=' is allowed but ignored
        }
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

// =================== Tests ===================

#[cfg(test)]
mod tests {
    use super::*;
    use gammalooprs::utils::test_utils::load_generic_model;
    use gammalooprs::{feyngen::GenerationType, initialisation::test_initialise};

    fn base_args(tokens: &str) -> SpecArgs {
        SpecArgs {
            tokens: tokens.split_whitespace().map(|x| x.to_string()).collect(),
            num_threads: 1,
            clear_existing_processes: false,
            filter_selfenergies: None,
            filter_snails: None,
            filter_tadpoles: None,
            max_n_bridges: None,
            number_of_factorized_loop_subtopologies: None,
            number_of_anticommutating_loops: None,
            allow_symmetrization_of_external_fermions_in_amplitudes: false,
            symmetrize_initial_states: false,
            symmetrize_final_states: false,
            symmetrize_left_right_states: false,
            numerator_aware_isomorphism_grouping: None,
            numerical_samples_seed: None,
            number_of_samples_for_numerator_comparisons: None,
            consider_internal_masses_only_in_numerator_isomorphisms: false,
            fully_numerical_substitution_when_comparing_numerators: true,
            compare_canonized_numerator: false,
            graph_prefix: None,
            loop_momentum_bases: None,
            select_graphs: None,
            veto_graphs: None,
            filter_self_loop: None,
        }
    }

    // Test helpers using a real generic model shipped with the crate
    fn parse_ok_amp(s: &str) -> ProcessSpec {
        let model = &load_generic_model("sm");
        parse_spec_with_model(&base_args(s), true, model).unwrap()
    }

    fn parse_ok_xs(s: &str) -> ProcessSpec {
        let model = &load_generic_model("sm");
        parse_spec_with_model(&base_args(s), false, model).unwrap()
    }

    #[test]
    fn basic_process_list_arrow() {
        test_initialise().unwrap();
        let ps = parse_ok_amp("e+ e- > d d~ g");
        assert_eq!(ps.initial, vec!["e+", "e-"]);
        assert_eq!(ps.final_, vec!["d", "d~", "g"]);
        assert!(!ps.empty_initial && !ps.empty_final);
    }

    #[test]
    fn basic_process_list_to() {
        test_initialise().unwrap();
        let ps = parse_ok_amp("e+ e- to d d~ g");
        assert_eq!(ps.initial, vec!["e+", "e-"]);
        assert_eq!(ps.final_, vec!["d", "d~", "g"]);
    }

    #[test]
    fn empty_sets_allowed() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("{} to {}");
        assert!(ps.empty_initial);
        assert!(ps.empty_final);
    }

    #[test]
    fn alternatives_in_final_states() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("e+ e- > { Z Z, a a }");
        assert!(ps.final_.is_empty());
        assert_eq!(ps.final_sets.len(), 2);
        assert_eq!(ps.final_sets[0], vec!["Z", "Z"]);
        assert_eq!(ps.final_sets[1], vec!["a", "a"]);
    }

    #[test]
    fn veto_and_only() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("e+ e- > mu+ mu- / u d c g ghG e- | u g ghG");
        assert!(ps.veto.contains("u"));
        assert!(ps.veto.contains("e-"));
        let sel = ps.only.as_ref().unwrap();
        assert!(sel.contains("u"));
        assert!(sel.contains("g"));
        assert!(sel.contains("ghG"));
    }

    #[test]
    fn amp_order_constraints() {
        test_initialise().unwrap();
        let ps = parse_ok_amp("e+ e- > mu+ mu- QED==2 QCD>=2 QCD<=4");
        let qcd = ps.amp_couplings.get("QCD").unwrap();
        assert_eq!(qcd.min, Some(2));
        assert_eq!(qcd.max, Some(4));
        let qed = ps.amp_couplings.get("QED").unwrap();
        assert_eq!(qed.eq, Some(2));
    }

    #[test]
    fn xs_order_constraints_powered() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("e+ e- > mu+ mu- QED^2==2 QCD^2>=2 QCD^2<=4");
        let key = CouplingKey {
            name: "QCD".into(),
            power: 2,
        };
        let qcd2 = ps.xs_couplings.get(&key).unwrap();
        assert_eq!(qcd2.min, Some(2));
        assert_eq!(qcd2.max, Some(4));
        let keyq = CouplingKey {
            name: "QED".into(),
            power: 2,
        };
        let qed2 = ps.xs_couplings.get(&keyq).unwrap();
        assert_eq!(qed2.eq, Some(2));
    }

    #[test]
    fn perturbative_block_all_forms() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("e+ e- > z [ {1} {{2}} QCD=2 QED=1 ]");
        assert_eq!(ps.pert.loops_sum_amp_or_sum, Some(1));
        assert_eq!(ps.pert.loops_forward_graph, Some(2));
        assert_eq!(ps.pert.orders.get("QCD"), Some(&2));
        assert_eq!(ps.pert.orders.get("QED"), Some(&1));
    }

    #[test]
    fn perturbative_block_shorthand() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("e+ e- > z [ QCD ]");
        assert_eq!(ps.pert.orders.get("QCD"), Some(&1));
    }

    #[test]
    fn amplitude_generation_type_is_set() {
        test_initialise().unwrap();
        let ps = parse_ok_amp("e+ e- > mu+ mu-");
        assert_eq!(ps.feyngen.generation_type, GenerationType::Amplitude);
    }

    #[test]
    fn rejects_malformed_tokens() {
        test_initialise().unwrap();
        let model = &load_generic_model("sm");

        let args = base_args("e+ > z [ {{}} ]");
        let err = parse_spec_with_model(&args, false, model).unwrap_err();
        assert_eq!(err, ParseError::InvalidToken("{{}}".into()));
    }

    #[test]
    fn missing_arrow_detected() {
        test_initialise().unwrap();
        let model = &load_generic_model("sm");

        let args = base_args("e+ e-  z z");
        let err = parse_spec_with_model(&args, false, model).unwrap_err();
        assert!(matches!(err, ParseError::MissingArrow));
    }

    #[test]
    fn vacuum_defaults_and_filters_xs() {
        test_initialise().unwrap();
        let ps = parse_ok_xs("{} to {}");
        assert_eq!(ps.feyngen.generation_type, GenerationType::CrossSection);

        let xs_filters = &ps.feyngen.cross_section_filters.0;

        // Smart defaults for vacuum-like graphs
        assert!(xs_filters
            .iter()
            .any(|f| matches!(f, FeynGenFilter::MaxNumberOfBridges(0))));
        assert!(xs_filters
            .iter()
            .any(|f| matches!(f, FeynGenFilter::FactorizedLoopTopologiesCountRange((1, 1)))));

        // No tadpoles/snails/self-energies by default for vacuum
        assert!(!xs_filters
            .iter()
            .any(|f| matches!(f, FeynGenFilter::TadpolesFilter(_))));
        assert!(!xs_filters
            .iter()
            .any(|f| matches!(f, FeynGenFilter::ZeroSnailsFilter(_))));
        assert!(!xs_filters
            .iter()
            .any(|f| matches!(f, FeynGenFilter::SelfEnergyFilter(_))));
    }

    #[test]
    fn particle_veto_resolution_cross_section() {
        use std::collections::BTreeSet;
        test_initialise().unwrap();
        // Veto both particles and anti-particles, plus a charged lepton
        let ps = parse_ok_xs("e+ e- > mu+ mu- / u d u~ e+");
        let xs_filters = &ps.feyngen.cross_section_filters.0;

        let veto_pdgs = xs_filters
            .iter()
            .find_map(|f| match f {
                FeynGenFilter::ParticleVeto(v) => Some(v.clone()),
                _ => None,
            })
            .expect("ParticleVeto filter not found");
        let set: BTreeSet<i64> = veto_pdgs.into_iter().collect();

        // u=2, d=1, u~=-2, e+=-11
        for pdg in [2_i64, 1_i64, -2_i64, -11_i64] {
            assert!(set.contains(&pdg), "missing PDG {pdg} in veto");
        }
    }

    #[test]
    fn final_sets_and_orders_and_loops_cross_section() {
        test_initialise().unwrap();
        // Cross-section with final-state alternatives and perturbative block
        let ps = parse_ok_xs("e+ e- > { Z Z, a a, H H } [ {{3}} QCD=2 QED=1 ]");

        // Final-state alternatives captured
        assert_eq!(ps.final_sets.len(), 3);
        assert_eq!(ps.feyngen.generation_type, GenerationType::CrossSection);

        // XS loop count from {{3}}
        assert_eq!(ps.feyngen.loop_count_range, (3, 3));

        // Perturbative orders end up in XS filters
        let xs_filters = &ps.feyngen.cross_section_filters.0;
        assert!(xs_filters.iter().any(|f| {
            if let FeynGenFilter::PerturbativeOrders(m) = f {
                m.get("QCD") == Some(&2usize) && m.get("QED") == Some(&1usize)
            } else {
                false
            }
        }));
    }

    #[test]
    fn amplitude_with_powered_xs_constraints_and_only() {
        test_initialise().unwrap();
        // AMP spec including "only" and XS-style powered coupling constraints
        let ps = parse_ok_amp("e+ e- > mu+ mu- | g ghG QED^2==2 QCD^2>=2 QCD^2<=4 [ {1} QCD ]");

        // AMP-side perturbative bits
        assert_eq!(ps.pert.loops_sum_amp_or_sum, Some(1));
        assert_eq!(ps.pert.orders.get("QCD"), Some(&1));

        // XS-style powered constraints recorded
        let key_qed2 = CouplingKey {
            name: "QED".into(),
            power: 2,
        };
        let key_qcd2 = CouplingKey {
            name: "QCD".into(),
            power: 2,
        };
        let qed2 = ps.xs_couplings.get(&key_qed2).unwrap();
        assert_eq!(qed2.eq, Some(2));
        let qcd2 = ps.xs_couplings.get(&key_qcd2).unwrap();
        assert_eq!(qcd2.min, Some(2));
        assert_eq!(qcd2.max, Some(4));

        // "only" selection recorded
        let only = ps.only.as_ref().expect("missing only-set");
        assert!(only.contains("g"));
        assert!(only.contains("ghG"));
    }

    #[test]
    fn amplitude_case_insensitive_and_anti_in_veto() {
        use std::collections::BTreeSet;
        test_initialise().unwrap();
        // Mixed case names and anti-particles in veto; ensure AMP-side veto filter present
        let ps = parse_ok_amp("E+ e- to Z Z / u~ d~ A");

        let amp_filters = &ps.feyngen.amplitude_filters.0;
        let veto_pdgs = amp_filters
            .iter()
            .find_map(|f| match f {
                FeynGenFilter::ParticleVeto(v) => Some(v.clone()),
                _ => None,
            })
            .expect("AMP ParticleVeto filter not found");

        let set: BTreeSet<i64> = veto_pdgs.into_iter().collect();
        // u~=-2, d~=-1, A(=a)=22
        for pdg in [-2_i64, -1_i64, 22_i64] {
            assert!(set.contains(&pdg), "missing PDG {pdg} in AMP veto");
        }
    }
}
