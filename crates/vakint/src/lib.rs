#![allow(mixed_script_confusables)]
pub mod alphaloop_numerics;
mod citations;
#[cfg(feature = "feynkit-ingress")]
pub mod diagram_integral;
pub mod fmft;
pub mod fmft_numerics;
pub mod graph;
mod hpl;
pub mod inspection;
mod lorentz;
#[cfg(test)]
mod master_precision;
pub mod matad;
pub mod matad_numerics;
pub mod precision;
mod presentation;
pub mod rustred_evaluation;
pub mod symbols;
pub mod topologies;
pub mod utils;

use ahash::RandomState;
use colored::Colorize;
use eyre::Result;
use feynkit_tensor::TensorReducer;
use graph::Graph;
pub use inspection::{
    ConstantSymbol, DataProvenance, MasterCoefficient, MasterSymbol, MasterSymbolicExpression,
    VakintConstant, VakintMaster,
};
#[allow(unused)]
use log::{debug, info, warn};
pub use precision::{
    AvailablePrecision, CoefficientPrecision, ComponentPrecision, ErrorQuality,
    NumericalCoefficient, NumericalEvaluationResult, PrecisionSource, PrecisionWarning,
};
pub use rustred::scalar_numerator::ScalarNumeratorLimits as RustRedScalarNumeratorLimits;
use rustred::solver::CandidateRoutedApplicationObserver;
pub use rustred::solver::{
    CandidateRoutedCachePolicy as RustRedCachePolicy,
    CandidateRoutedCacheStatistics as RustRedCacheWorkStatistics,
};
pub use rustred_evaluation::cache::{
    RustRedCacheFileLimits, RustRedCacheFileOptions, RustRedCacheLoadOptions,
    RustRedCacheNamespaceStatistics, RustRedCachePreparationOptions, RustRedCacheRoot,
    RustRedCacheSaveMode, RustRedCacheSaveOptions, RustRedCacheSaveReport, RustRedCacheStatistics,
    RustRedCacheWorkload, RustRedRunningCacheControl, RustRedRunningCacheOptions,
    RustRedRunningCacheReport, RustRedRunningCacheSnapshot, RustRedSavedCacheStatistics,
};
#[cfg(all(unix, not(target_arch = "wasm32")))]
pub use rustred_evaluation::cache::{RustRedCacheTransferOptions, RustRedCacheTransferReport};
pub use rustred_evaluation::five_loop::{RustRedSymbolicEvaluation, RustRedWitnessCacheStatistics};
use rustred_evaluation::monitor::MonitorContext;
pub use rustred_evaluation::monitor::{
    RustRedMonitorDisplay, RustRedMonitorOptions, RustRedSnapshotViewer,
};
#[cfg(all(unix, not(target_arch = "wasm32")))]
pub use rustred_evaluation::{
    ExternalDataLoadReport, ExternalDataLoopReport, OfficialExternalDataOptions,
};
pub use rustred_evaluation::{RustRedEvaluationError, RustRedEvaluationOptions, RustRedMasterMode};

use regex::Regex;
use rug::float::Constant;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing::{Collectable, TensorCollectFilter},
    structure::representation::{BaseRepName, Minkowski, initialize as initialize_spenso_reps},
};
use std::{
    collections::{BTreeMap, HashMap, HashSet, hash_map::Entry},
    env,
    f64::consts::LOG2_10,
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    path::PathBuf,
    process::{Command, ExitStatus, Stdio},
    sync::{
        Arc, LazyLock, Mutex,
        atomic::{AtomicU64, Ordering},
    },
    time::{SystemTime, UNIX_EPOCH},
    vec,
};
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::{
    coefficient::{Coefficient, CoefficientView},
    domains::float::{FloatLike, RealLike},
};
#[cfg(feature = "symbolica_community_module")]
pub mod symbolica_community_module;

/// Gather the Python API registered for the Vakint Symbolica community module.
#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    pyo3_stub_gen::StubInfo::from_project_root(
        "symbolica".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")),
    )
}

use crate::alphaloop_numerics::DIRECT_SUBSTITUTIONS;

#[allow(unused)]
use symbolica::{
    atom::{
        AliasedAtom, Atom, AtomCore, AtomView, FunctionBuilder, SliceType, Symbol, Var,
        representation::InlineNum,
    },
    domains::{
        atom::AtomField,
        float::{Complex, Float, Real, SingleFloat},
        integer::IntegerRing,
        rational::{Fraction, Rational},
    },
    function,
    id::{
        Condition, Match, MatchSettings, Pattern, PatternRestriction, Replacement,
        WildcardRestriction,
    },
    poly::series::Series,
    printer::{AtomPrinter, CanonicalOrderingSettings, PrintOptions},
    symbol,
    transformer::Transformer,
};
use utils::simplify_real;
use utils::vakint_macros::{vk_parse, vk_symbol};

pub use symbolica::utils::resources::ResourceLimits;
use symbolica::utils::resources::{self, ResourceBudget};
use thiserror::Error;
use topologies::{Topologies, Topology};
use version_compare::{Cmp, compare_to};

#[allow(unused)]
use symbols::SYMBOL_REGISTRY;
use symbols::{EXTERNAL_MOMENTUM_SYMBOL, LOOP_MOMENTUM_SYMBOL, METRIC_SYMBOL, S};

use phf::phf_map;

use crate::{
    symbols::DOT_SYMBOL,
    utils::{get_full_name, pysecdec_decode, pysecdec_encode, to_symbol, undress_vakint_symbols},
};

pub enum Momentum {
    Complex(ComplexMomentum),
    Real(RealMomentum),
}

pub type ComplexMomentum = (
    Complex<Float>,
    Complex<Float>,
    Complex<Float>,
    Complex<Float>,
);

pub type RealMomentum = (Float, Float, Float, Float);

pub static FORM_REPLACEMENT_INDEX_SHIFT: u64 = 13370000;

pub static NAMESPACE: &str = "vakint";
static FEYNKIT_EXTERNAL_VECTOR: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
        "vakint::feynkit_external_vector",
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    )
});
static FEYNKIT_OPAQUE_INDEX: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("vakint::feynkit_opaque_index"));
// pub static PYSECDEC_NAMESPACE_SEPARATOR: &str = "NAMESPACESEP";
// pub static PYSECDEC_UNDERSCORE: &str = "UNDERSCORE";
// pub static PYSECDEC_ATTRIBUTE_START: &str = "ATTRIBUTESTART";
// pub static PYSECDEC_ATTRIBUTE_SEPARATOR: &str = "ATTRIBUTESEP";
// pub static PYSECDEC_ATTRIBUTE_END: &str = "ATTRIBUTESEND";
pub static PYSECDEC_NAMESPACE_SEPARATOR: &str = "PSDNS";
pub static PYSECDEC_UNDERSCORE: &str = "PSDU";
pub static PYSECDEC_ATTRIBUTE_START: &str = "PSDAB";
pub static PYSECDEC_ATTRIBUTE_SEPARATOR: &str = "PSDAS";
pub static PYSECDEC_ATTRIBUTE_END: &str = "PSDAE";
static TMP_DIR_UNIQUIFIER: AtomicU64 = AtomicU64::new(0);

// pub static VAKINT_SYMBOL_UNDRESSING_RE: LazyLock<Regex> = LazyLock::new(|| {
//     Regex::new(&format!(r"{}::\{{[^{{}}]*\}}::", regex::escape(NAMESPACE)))
//         .expect("invalid regex")
// });
pub static VAKINT_SYMBOL_UNDRESSING_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(&format!(
        r"{}::(?:\{{[^{{}}]*\}}::)?",
        regex::escape(NAMESPACE)
    ))
    .expect("invalid regex")
});

fn create_unique_tmp_run_dir(base_dir: Option<&str>, prefix: &str) -> std::io::Result<PathBuf> {
    let root = base_dir.map(PathBuf::from).unwrap_or_else(env::temp_dir);
    fs::create_dir_all(&root)?;

    for _ in 0..1024 {
        let sequence = TMP_DIR_UNIQUIFIER.fetch_add(1, Ordering::Relaxed);
        let timestamp_ns = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();

        let candidate = root.join(format!(
            "{prefix}_pid{}_{}_{}",
            std::process::id(),
            timestamp_ns,
            sequence
        ));
        match fs::create_dir(&candidate) {
            Ok(()) => return Ok(candidate),
            Err(err) if err.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(err) => return Err(err),
        }
    }

    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!("Could not allocate a unique temporary directory for '{prefix}'"),
    ))
}

static MINIMAL_FORM_VERSION: &str = "4.2.1";
static MINIMAL_PYSECDEC_VERSION: &str = "1.6.4";

static FORM_SRC: phf::Map<&'static str, &'static str> = phf_map! {
    "integrateduv.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/alphaloop/integrateduv.frm"
    )),
    "tensorreduce.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/alphaloop/tensorreduce.frm"
    )),
    "pvtab10.h" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/alphaloop/pvtab10.h"
    )),
    "fmft.frm" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/fmft/fmft.frm"
    )),
    "matad-ng.hh" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/form_src/matad/matad-ng.hh"
    )),
};

static TEMPLATES: phf::Map<&'static str, &'static str> = phf_map! {
    "run_tensor_reduction.txt" =>  include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_tensor_reduction.txt"
    )),
    "run_alphaloop_integral_evaluation.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_alphaloop_integral_evaluation.txt"
    )),
    "run_pySecDec_template.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_pySecDec_template.txt"
    )),
    "run_matad.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_matad.txt"
    )),
    "run_fmft.txt" => include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/templates/run_fmft.txt"
    )),
};

#[derive(Error, Debug)]
pub enum VakintError {
    #[error("{0}")]
    ResourceLimit(#[from] resources::ResourceLimitError),
    #[error("invalid input integral found: '{0}'")]
    InvalidIntegralFormat(String),
    #[error("invalid generic expression in supported integral: {0}")]
    InvalidGenericExpression(String),
    #[error("invalid expression for the momentum of an edge: {0}")]
    InvalidMomentumExpression(String),
    #[error("invalid short expression supplied for integral: {0}")]
    InvalidShortExpression(String),
    #[error("invalid numerator expression: {0}")]
    InvalidNumerator(String),
    #[error("Could not find a method suitable for evaluating this integral up to 𝒪(ε^{1}): {0}")]
    NoEvaluationMethodFound(String, i64),
    #[error(
        "the following integral could not be identified using any of the supported topologies: {0}"
    )]
    UnreckognizedIntegral(String),
    #[error(
        "not all parts of the numerator have been identified w.r.t the canonical \
        expression (make sure to use the name '{0}' for momenta external to the topology).\n\
        This check can be disabled by setting `verify_numerator_identification` to false.\
        \nLeft-over: {1}"
    )]
    NumeratorNotReplaced(String, String),
    #[error(
        "FORM run crashed with the following error:\nstderr: {0}\nstdout: {1}\nYou can rerun the script using:\n{2}"
    )]
    FormError(String, String, String, String),
    #[error("{0}")]
    FormVersion(String),
    #[error("FORM is not installed in your system and required for vakint to work.")]
    FormUnavailable,
    #[error("PySecDec error: {0}")]
    PySecDecError(String),
    #[error("{0}")]
    PySecDecVersion(String),
    #[error(
        "PySecDec is not installed in your system and required for vakint to evaluate with PySecDec. Install it with 'pip install pysecdec'"
    )]
    PySecDecUnavailable,
    #[error(
        "Could not find FORM output file 'out.txt':\nstderr: {0}\nYou can rerun the script using:\n{1}"
    )]
    MissingFormOutput(String, String, String),
    #[error("Symbolica could not parse PySecDec output:\n{0}\nError:{1}")]
    PySecDecOutputParsingError(String, String),
    #[error("Symbolica could not parse FORM output:\n{0}\nError:{1}")]
    FormOutputParsingError(String, String),
    #[error(transparent)]
    IoError(#[from] std::io::Error), // Add this variant to convert std::io::Error to VakintError
    #[error("{0}")]
    MalformedGraph(String),
    #[error("native momentum matching exceeded {resource}: {requested} > {limit}")]
    MomentumMatchResourceLimit {
        resource: &'static str,
        requested: usize,
        limit: usize,
    },
    #[error(
        "Invalid loop normalization factor specified: {0}.\nError: {1}\nNote that the following symbols are allowed in the expression: {2}"
    )]
    InvalidLoopNormalization(String, String, String),
    #[error("Symbolica error: {0}")]
    SymbolicaError(String),
    #[error("MATAD error: {0}")]
    MATADError(String),
    #[error("FMFT error: {0}")]
    FMFTError(String),
    #[error("Numerical evaluation error: {0}")]
    EvaluationError(String),
    #[error(transparent)]
    FeynKitTensorReduction(#[from] feynkit_tensor::TensorReductionError),
    #[error("temporary FeynKit tensor-bridge syntax remained in the reduced numerator: {0}")]
    FeynKitBridgeSyntax(String),
    #[error("unsupported Spenso dot containing an integrated loop momentum: {0}")]
    FeynKitLoopDot(String),
    #[error("unknown numerator tensor reduction method '{0}'; expected 'alphaloop' or 'feynkit'")]
    InvalidTensorReductionMethod(String),
    #[error(transparent)]
    RustRedEvaluation(#[from] RustRedEvaluationError),
    #[error(transparent)]
    Inspection(#[from] inspection::InspectionError),
    #[error("numerical precision: {0}")]
    NumericalPrecision(String),
    #[error("unknown vakint error")]
    Unknown,
}

pub fn params_from_f64(
    params: &HashMap<String, f64, ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<String, Float, ahash::RandomState> {
    let binary_prec: u32 = (((decimal_prec as f64) * LOG2_10).floor() as u32).min(53);
    params
        .iter()
        .map(|(k, v)| (k.clone(), Float::with_val(binary_prec, v)))
        .collect()
}

pub fn params_from_complex_f64(
    params: &HashMap<String, Complex<f64>, ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<String, Complex<Float>, ahash::RandomState> {
    let binary_prec: u32 = (((decimal_prec as f64) * LOG2_10).floor() as u32).min(53);
    params
        .iter()
        .map(|(k, v)| {
            (
                k.clone(),
                Complex::new(
                    Float::with_val(binary_prec, v.re),
                    Float::with_val(binary_prec, v.im),
                ),
            )
        })
        .collect()
}

pub fn externals_from_f64(
    externals: &HashMap<usize, (f64, f64, f64, f64), ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<usize, Momentum, ahash::RandomState> {
    let binary_prec: u32 = (((decimal_prec as f64) * LOG2_10).floor() as u32).min(53);
    externals
        .iter()
        .map(|(&k, v)| {
            (
                k,
                Momentum::Real((
                    Float::with_val(binary_prec, v.0),
                    Float::with_val(binary_prec, v.1),
                    Float::with_val(binary_prec, v.2),
                    Float::with_val(binary_prec, v.3),
                )),
            )
        })
        .collect()
}

#[allow(clippy::type_complexity)]
pub fn externals_from_complex_f64(
    externals: &HashMap<
        usize,
        (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>),
        ahash::RandomState,
    >,
    decimal_prec: u32,
) -> HashMap<usize, Momentum, ahash::RandomState> {
    let binary_prec: u32 = (((decimal_prec as f64) * LOG2_10).floor() as u32).min(53);
    externals
        .iter()
        .map(|(&k, v)| {
            (
                k,
                Momentum::Complex((
                    Complex::new(
                        Float::with_val(binary_prec, v.0.re),
                        Float::with_val(binary_prec, v.0.im),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.1.re),
                        Float::with_val(binary_prec, v.1.im),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.2.re),
                        Float::with_val(binary_prec, v.2.im),
                    ),
                    Complex::new(
                        Float::with_val(binary_prec, v.3.re),
                        Float::with_val(binary_prec, v.3.im),
                    ),
                )),
            )
        })
        .collect()
}

fn propagators_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        let props = match m {
            Match::Single(a) => vec![*a],
            Match::Multiple(SliceType::Mul, atoms) => atoms.clone(),
            _ => return false,
        };

        let pattern = vk_parse!("prop(propID_,uedge(nl_,nr_),q_,mUVsq_,pow_)")
            .unwrap()
            .to_pattern();
        let number_node_condition = Condition::from((vk_symbol!("nl_"), ge_condition(0)))
            & Condition::from((vk_symbol!("nr_"), ge_condition(0)))
            & Condition::from((vk_symbol!("propID_"), ge_condition(0)))
            // DO NOT REQUIRE MASS TO BE A SYMBOL
            // & Condition::from((vk_symbol!("mUVsq_"), symbol_or_number()))
            & Condition::from((vk_symbol!("pow_"), symbol_or_number()));
        for p in props {
            if p.pattern_match(&pattern, Some(&number_node_condition), None)
                .next()
                .is_none()
            {
                return false;
            }
        }
        true
    }))
}

fn ge_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() >= InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn eq_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() == InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

#[allow(unused)]
fn lt_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() < InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn range_condition(min: i64, max: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() <= InlineNum::new(max, 1).as_num_view().get_coeff_view()
                && a.get_coeff_view() >= InlineNum::new(min, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
    }))
}

fn number_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        matches!(m, Match::Single(AtomView::Num(_)))
    }))
}

fn even_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        #[allow(warnings)]
        if let Match::Single(AtomView::Num(_)) = m {
            get_integer_from_atom(m.to_atom().as_view()).unwrap() % 2 == 0
        } else {
            false
        }
    }))
}

fn odd_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        #[allow(warnings)]
        if let Match::Single(AtomView::Num(_)) = m {
            get_integer_from_atom(m.to_atom().as_view()).unwrap() % 2 == 1
        } else {
            false
        }
    }))
}

fn symbol_or_number() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        matches!(
            m,
            Match::Single(AtomView::Num(_))
                | Match::Single(AtomView::Var(_))
                | Match::FunctionName(_)
        )
    }))
}

fn symbol_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        matches!(m, Match::Single(AtomView::Var(_)) | Match::FunctionName(_))
    }))
}

fn function_condition() -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| matches!(m, Match::FunctionName(_))))
}

fn apply_restriction_to_symbols(
    symbols: Vec<Symbol>,
    restriction: &WildcardRestriction,
) -> Condition<PatternRestriction> {
    symbols[1..].iter().fold(
        Condition::from((symbols[0], restriction.clone())),
        |acc, &s| acc & Condition::from((s, restriction.clone())),
    )
}

#[derive(Debug, Clone)]
pub struct ReplacementRules {
    canonical_topology: Topology,
    edge_ids_canonical_to_input_map: HashMap<usize, usize>,
    canonical_expression_substitutions: HashMap<Atom, Atom>,
    numerator_substitutions: HashMap<Atom, Atom>,
}

impl fmt::Display for ReplacementRules {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "ReplacementRules {{")?;
        writeln!(f, "  canonical_topology: {}", self.canonical_topology)?;
        let mut sorted_edge_ids_canonical_to_input_map = self
            .edge_ids_canonical_to_input_map
            .iter()
            .collect::<Vec<_>>();
        sorted_edge_ids_canonical_to_input_map.sort_by_key(|k| *k.0);
        writeln!(
            f,
            "  edge_ids_canonical_to_input_map: {{ {} }}",
            sorted_edge_ids_canonical_to_input_map
                .iter()
                .map(|(k, v)| format!("{} -> {}", k, v))
                .collect::<Vec<String>>()
                .join(", ")
        )?;
        let mut sorted_canonical_expression_substitutions = self
            .canonical_expression_substitutions
            .iter()
            .collect::<Vec<_>>();
        sorted_canonical_expression_substitutions.sort_by_key(|k| k.0);
        writeln!(
            f,
            "  canonical_expression_substitutions: {{ {} }}",
            sorted_canonical_expression_substitutions
                .iter()
                .map(|(k, v)| format!("{} -> {}", k, v))
                .collect::<Vec<String>>()
                .join(", ")
        )?;
        let mut sorted_numerator_substitutions =
            self.numerator_substitutions.iter().collect::<Vec<_>>();
        sorted_numerator_substitutions.sort_by_key(|k| k.0);
        writeln!(
            f,
            "  numerator_substitutions: {{ {} }}",
            sorted_numerator_substitutions
                .iter()
                .map(|(k, v)| format!("{} -> {}", k, v))
                .collect::<Vec<String>>()
                .join(", ")
        )?;
        write!(f, "}}")
    }
}
impl Default for ReplacementRules {
    fn default() -> Self {
        ReplacementRules {
            canonical_topology: Topology::Unknown(Integral::default()),
            edge_ids_canonical_to_input_map: HashMap::new(),
            canonical_expression_substitutions: HashMap::new(),
            numerator_substitutions: HashMap::new(),
        }
    }
}

impl ReplacementRules {
    fn apply_replacement_rules(&mut self) -> Result<(), VakintError> {
        if matches!(self.canonical_topology, Topology::Unknown(_)) {
            let integral = self.canonical_topology.get_integral_mut();
            integral.canonical_expression = Some(
                self.canonical_expression_substitutions
                    .get(&vk_parse!("integral").unwrap())
                    .unwrap()
                    .to_owned(),
            );
            let n_props: i64 = self
                .canonical_expression_substitutions
                .get(&vk_parse!("n_props").unwrap())
                .unwrap()
                .try_into()
                .unwrap();
            integral.n_props = n_props as usize;
            let n_loops: i64 = self
                .canonical_expression_substitutions
                .get(&vk_parse!("n_loops").unwrap())
                .unwrap()
                .try_into()
                .unwrap();
            integral.n_loops = n_loops as usize;
            integral.graph = Graph::new_from_atom(
                integral.canonical_expression.as_ref().unwrap().as_view(),
                integral.n_props,
            )?;
            return Ok(());
        }

        let integral = self.canonical_topology.get_integral_mut();

        for (source, target) in self.canonical_expression_substitutions.iter() {
            for expr in [
                integral.canonical_expression.as_mut(),
                integral.short_expression.as_mut(),
                integral.alphaloop_expression.as_mut(),
            ]
            .into_iter()
            .flatten()
            {
                *expr = expr.replace(source.to_pattern()).with(target.to_pattern());
            }
        }
        Ok(())
    }

    fn get_propagator_property_list(&self, property: &str) -> HashMap<usize, Atom> {
        let mut property_list = HashMap::new();
        let integral = self.canonical_topology.get_integral();
        let canonical_expression_view = integral.canonical_expression.as_ref().unwrap().as_view();
        for prop_id in 1..=integral.n_props {
            if let Some(m) = get_prop_with_id(canonical_expression_view, prop_id) {
                // println!("property {}", property);
                // println!(
                //     "key {}",
                //     &m.get(&vk_symbol!(property)).unwrap().to_atom()
                // );
                // println!("self.canonical_expression_substitution {}", self);
                let v_a = m.get(&vk_symbol!(property)).unwrap();
                property_list.insert(prop_id, v_a.to_owned());
            }
        }
        property_list
    }
}

#[derive(Debug, Clone)]
pub struct Integral {
    n_loops: usize,
    n_props: usize,
    name: String,
    generic_pattern: FullPattern,
    canonical_expression: Option<Atom>,
    // The uncontracted defining topology and k_new(i) in its loop basis.
    // Populated from the actual selected LMB during contraction, not inferred
    // later by rematching graphs or solving another momentum system.
    parent_routing: Option<(Atom, Arc<[Atom]>)>,
    short_expression: Option<Atom>,
    short_expression_pattern: Option<Pattern>,
    alphaloop_expression: Option<Atom>,
    applicable_evaluation_methods: EvaluationOrder,
    graph: Graph,
    node_pairs: HashMap<Symbol, HashSet<Symbol>>,
    unoriented_generic_pattern: FullPattern,
}

impl fmt::Display for Integral {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "name='{}', n_loops={}, n_props_top_topo={}\n   | canonical_expression='{}',\n   | canonical_pattern='{}',\n   | short_expression='{}'",
            if self.name == "UNKNOWN" {
                self.name.red()
            } else {
                self.name.green()
            },
            self.n_loops,
            self.n_props,
            self.canonical_expression
                .as_ref()
                .map(|e| format!("{}", e))
                .unwrap_or("N/A".into())
                .blue(),
            format!("{}", self.generic_pattern.pattern.to_atom().unwrap()).blue(),
            self.short_expression
                .as_ref()
                .map(|e| format!("{}", e))
                .unwrap_or("N/A".into())
                .green()
        )
    }
}

fn get_integer_from_atom(n: AtomView) -> Option<i64> {
    n.try_into().ok()
    /*
    if let Match::Single(AtomView::Num(a)) = m {
        match a.get_coeff_view() {
            CoefficientView::Natural(n, 1) => Some(n),
            _ => None,
        }
    } else {
        None
    }
    */
}

fn get_node_ids(
    match_stack: &HashMap<symbolica::atom::Symbol, Atom>,
) -> Result<(usize, usize), VakintError> {
    let id_node_left = if let Some(id) =
        get_integer_from_atom(match_stack.get(&vk_symbol!("nl_")).unwrap().as_view())
    {
        id
    } else {
        return Err(VakintError::InvalidGenericExpression(format!(
            "Left node must be an integer: {}",
            match_stack.get(&vk_symbol!("nl_")).unwrap()
        )));
    };
    let id_node_right = if let Some(id) =
        get_integer_from_atom(match_stack.get(&vk_symbol!("nr_")).unwrap().as_view())
    {
        id
    } else {
        return Err(VakintError::InvalidGenericExpression(format!(
            "Right node must be an integer: {}",
            match_stack.get(&vk_symbol!("nl_")).unwrap()
        )));
    };

    Ok((id_node_left as usize, id_node_right as usize))
}

#[allow(clippy::type_complexity)]
fn get_individual_momenta_from_atom(
    momentum: AtomView,
) -> Result<Vec<(Symbol, (Atom, i64))>, VakintError> {
    let wrapped_momentum = function!(vk_symbol!("mom"), momentum.to_owned());
    if let Some(m) = wrapped_momentum
        .pattern_match(&vk_parse!("mom(q_)").unwrap().to_pattern(), None, None)
        .next()
    {
        get_individual_momenta(m.get(&vk_symbol!("q_")).unwrap().as_view())
    } else {
        Err(VakintError::InvalidMomentumExpression(format!(
            "Edge momentum {} does not contain only 'q_'.",
            momentum
        )))
    }
}

#[allow(clippy::type_complexity)]
fn get_individual_momenta(momentum: AtomView) -> Result<Vec<(Symbol, (Atom, i64))>, VakintError> {
    let err = Err(VakintError::InvalidMomentumExpression(format!(
        "Edge momentum {} does not contain only <symbol>_(<int>).",
        momentum
    )));
    let mut res = vec![];
    let mut test = momentum.to_owned();

    while let Some(m) = test
        .pattern_match(
            &vk_parse!("q_(lmbID_)").unwrap().to_pattern(),
            Some(
                &(Condition::from((vk_symbol!("q_"), symbol_condition()))
                    & Condition::from((vk_symbol!("lmbID_"), number_condition()))),
            ),
            None,
        )
        .next_detailed()
    {
        let a_match = m.match_stack.get(vk_symbol!("lmbID_")).unwrap();
        let (mom_symbol, (atom_id, mom_id)) = (
            if let Match::FunctionName(s) = m.match_stack.get(vk_symbol!("q_")).unwrap() {
                *s
            } else {
                return err;
            },
            if let Match::Single(a) = a_match {
                if let Some(id) = get_integer_from_atom(a_match.to_atom().as_view()) {
                    (a.to_owned(), id)
                } else {
                    return err;
                }
            } else {
                return err;
            },
        );

        test = test
            .replace(function!(mom_symbol, &atom_id).to_pattern())
            .with(Atom::Zero.to_pattern());

        res.push((mom_symbol, (atom_id, mom_id)));
    }

    if !test.is_zero() {
        return err;
    }

    Ok(res)
}

fn get_prop_with_id(
    expression: AtomView,
    prop_id: usize,
) -> Option<HashMap<symbolica::atom::Symbol, Atom>> {
    let pattern: Pattern =
        vk_parse!(format!("prop({},edge(nl_,nr_),q_,mUVsq_,pow_)", prop_id).as_str())
            .unwrap()
            .to_pattern();
    let number_node_condition = Condition::from((vk_symbol!("nl_"), ge_condition(0)))
        & Condition::from((vk_symbol!("nr_"), ge_condition(0)));
    expression
        .pattern_match(&pattern, Some(&number_node_condition), None)
        .next()
        .map(|m| {
            m.iter()
                .map(|(k, v)| (*k, v.to_owned()))
                .collect::<HashMap<_, _>>()
        })
}

impl Default for Integral {
    fn default() -> Self {
        Integral::new(0, None, None, EvaluationOrder::empty()).unwrap()
    }
}

impl Integral {
    pub fn new(
        // This refers to the *total* number of propagators in the top-level topology,
        // i.e. the number of entries in the corresponding short_expression
        tot_n_props: usize,
        canonical_expression: Option<Atom>,
        short_expression: Option<Atom>,
        applicable_evaluation_methods: EvaluationOrder,
    ) -> Result<Integral, VakintError> {
        if canonical_expression.is_none() {
            // This is an unknown topology
            let all_accepting_pattern = FullPattern {
                pattern: vk_parse!("topo(props_)").unwrap().to_pattern(),
                conditions: Condition::from((vk_symbol!("props_"), propagators_condition())),
                match_settings: MatchSettings::default(),
            };
            let short_expression_pattern = short_expression.as_ref().map(|a| a.to_pattern());
            return Ok(Integral {
                name: "Unknown".into(),
                n_loops: 0,
                n_props: 0,
                generic_pattern: all_accepting_pattern.clone(),
                canonical_expression: None,
                parent_routing: None,
                short_expression,
                short_expression_pattern,
                alphaloop_expression: None,
                applicable_evaluation_methods,
                graph: Graph::default(),
                unoriented_generic_pattern: all_accepting_pattern,
                node_pairs: HashMap::new(),
            });
        }
        let e = canonical_expression.clone().unwrap();

        let graph = Graph::new_from_atom(e.as_view(), tot_n_props)?;

        let mut loop_mom_indices = HashSet::<i64>::default();
        let mut next_atom = e.clone();
        let mut old_atom = e.clone();

        let mut generic_expression = vk_parse!("1").unwrap();
        let mut generic_condition: Condition<PatternRestriction> = Condition::default();
        let mut unoriented_generic_expression = vk_parse!("1").unwrap();
        let mut unoriented_generic_condition: Condition<PatternRestriction> = Condition::default();
        let mut node_pairs = HashMap::<Symbol, HashSet<Symbol>>::default();

        for i_prop in 1..=tot_n_props {
            if let Some(m) = get_prop_with_id(e.as_view(), i_prop) {
                // Format check for the momenta
                let momenta = get_individual_momenta(m.get(&vk_symbol!("q_")).unwrap().as_view())?;

                let mass_symbol_string = if let Some(a) = m.get(&vk_symbol!("mUVsq_")) {
                    if let Some(m3) = a
                        .pattern_match(
                            &vk_parse!("msq(mid_)").unwrap().to_pattern(),
                            Some(&Condition::from((
                                vk_symbol!("mid_"),
                                range_condition(1, tot_n_props as i64),
                            ))),
                            None,
                        )
                        .next_detailed()
                    {
                        format!(
                            "msq{}_",
                            get_integer_from_atom(
                                m3.match_stack
                                    .get(vk_symbol!("mid_"))
                                    .unwrap()
                                    .to_atom()
                                    .as_view(),
                            )
                            .unwrap()
                        )
                    } else {
                        return Err(VakintError::InvalidGenericExpression(format!(
                            "Generic expression does not have masses formatted as msq(integer in [1,n_props]): {}",
                            a
                        )));
                    }
                } else {
                    return Err(VakintError::InvalidGenericExpression(
                        "Generic expression does not have masses formatted as msq(<m_integer_id>)."
                            .into(),
                    ));
                };

                // If the power is a pattern, then match it
                let power_match = m.get(&vk_symbol!("pow_")).unwrap();
                let pow_symbol_string = if let Some(pow) =
                    get_integer_from_atom(power_match.as_view())
                {
                    format!("{}", pow)
                } else {
                    let a = power_match.as_view();
                    if let Some(m3) = a
                        .pattern_match(
                            &vk_parse!("pow(pid_)").unwrap().to_pattern(),
                            Some(&Condition::from((
                                vk_symbol!("pid_"),
                                range_condition(1, tot_n_props as i64),
                            ))),
                            None,
                        )
                        .next()
                    {
                        format!(
                            "pow{}_",
                            get_integer_from_atom(m3.get(&vk_symbol!("pid_")).unwrap().as_view(),)
                                .unwrap()
                        )
                    } else {
                        return Err(VakintError::InvalidGenericExpression(format!(
                            "Generic expression does not have powers formatted as pow(integer in [1,n_props]): {}",
                            a
                        )));
                    }
                };

                for (mom_symbol, (_mom_atom_id, _mom_id)) in momenta {
                    if mom_symbol != vk_symbol!("k") {
                        return Err(VakintError::InvalidGenericExpression(
                            "Generic expression does not have momenta involving only expressions of the type k(<integer>)."
                                .to_string(),
                        ));
                    }
                    loop_mom_indices.insert(_mom_id);
                }

                let (id_node_left, id_node_right) = get_node_ids(&m)?;

                generic_expression *=
                    vk_parse!(format!(
                    "prop(id{id}_,uedge(n{id_l_node}l_,n{id_r_node}r_),q{id}__,{mass},{power})",
                    id = i_prop,
                    id_l_node = id_node_left,
                    id_r_node = id_node_right,
                    mass = mass_symbol_string,
                    power = pow_symbol_string
                )
                .as_str())
                    .unwrap();
                unoriented_generic_expression *=
                    vk_parse!(format!(
                    "prop(id{id}_,uedge(n{id_l_node}_,n{id_r_node}_),q{id}__,{mass},{power})",
                    id = i_prop,
                    id_l_node = id_node_left,
                    id_r_node = id_node_right,
                    mass = mass_symbol_string,
                    power = pow_symbol_string
                )
                .as_str())
                    .unwrap();

                let new_conditions = Condition::from((
                    vk_symbol!(format!("pow{}_", i_prop).as_str()),
                    number_condition(),
                ))
                // DO NOT REQUIRE MASS TO BE A SYMBOL
                // & Condition::from((
                //     vk_symbol!(mass_symbol_string.clone()),
                //     symbol_or_number(),
                // ))
                & Condition::from((
                    vk_symbol!(format!("id{}_", i_prop).as_str()),
                    number_condition(),
                ));
                generic_condition = generic_condition
                    & new_conditions.clone()
                    & apply_restriction_to_symbols(
                        vec![
                            vk_symbol!(format!("n{}l_", id_node_left).as_str()),
                            vk_symbol!(format!("n{}r_", id_node_right).as_str()),
                        ],
                        &symbol_or_number(),
                    );
                unoriented_generic_condition = unoriented_generic_condition
                    & new_conditions
                    & apply_restriction_to_symbols(
                        vec![
                            vk_symbol!(format!("n{}_", id_node_left).as_str()),
                            vk_symbol!(format!("n{}_", id_node_right).as_str()),
                        ],
                        &symbol_or_number(),
                    );
                for (id, side) in [(id_node_left, "l"), (id_node_right, "r")] {
                    let new_entry = vk_symbol!(format!("n{}{}_", id, side));
                    node_pairs
                        .entry(vk_symbol!(format!("n{}_", id)))
                        .and_modify(|v| {
                            v.insert(new_entry);
                        })
                        .or_insert(HashSet::from([new_entry]));
                }
                next_atom = old_atom
                    .replace(
                        vk_parse!(format!("prop({},args__)", i_prop).as_str())
                            .unwrap()
                            .to_pattern(),
                    )
                    .with(vk_parse!("1").unwrap().to_pattern());
                old_atom = next_atom.clone();
            }
        }

        if next_atom != vk_parse!("topo(1)").unwrap() {
            return Err(VakintError::InvalidGenericExpression(format!(
                "Not all propagators of the generic expression supplied have been successfully identified. Left-over: {}",
                next_atom
            )));
        }
        generic_expression = function!(vk_symbol!("topo"), &generic_expression);
        let generic_pattern = FullPattern {
            pattern: generic_expression.to_pattern(),
            conditions: generic_condition,
            match_settings: MatchSettings::default(),
        };
        let unoriented_generic_pattern = FullPattern {
            pattern: unoriented_generic_expression.to_pattern(),
            conditions: unoriented_generic_condition,
            match_settings: MatchSettings::default(),
        };

        let name = if let Some(m) = short_expression
            .as_ref()
            .unwrap()
            .pattern_match(&vk_parse!("f_(args__)").unwrap().to_pattern(), None, None)
            .next_detailed()
        {
            if let Some(Match::FunctionName(s)) = m.match_stack.get(vk_symbol!("f_")) {
                s.to_string()
            } else {
                return Err(VakintError::InvalidShortExpression(
                    short_expression.unwrap().to_string(),
                ));
            }
        } else {
            return Err(VakintError::InvalidShortExpression(
                short_expression.unwrap().to_string(),
            ));
        };

        let mut short_expression_pattern = short_expression.clone().unwrap();
        for i_prop in 1..=tot_n_props {
            short_expression_pattern = short_expression_pattern
                .replace(
                    vk_parse!(format!("pow({})", i_prop).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .allow_new_wildcards_on_rhs(true)
                .with(
                    vk_parse!(format!("pow{}_", i_prop).as_str())
                        .unwrap()
                        .to_pattern(),
                );
            short_expression_pattern = short_expression_pattern
                .replace(
                    vk_parse!(format!("msq({})", i_prop).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .allow_new_wildcards_on_rhs(true)
                .with(
                    vk_parse!(format!("msq{}_", i_prop).as_str())
                        .unwrap()
                        .to_pattern(),
                );
        }

        let mut alphaloop_expression = Atom::num(1);

        for node in graph.nodes.values() {
            let mut fb = FunctionBuilder::new(vk_symbol!("vxs"));
            for (e_id, dir) in &node.edges {
                fb = fb.add_arg(
                    &(graph.edges.get(e_id).unwrap().momentum.clone()
                        * if dir.is_incoming() { 1 } else { -1 }),
                );
            }
            alphaloop_expression *= fb.finish();
        }
        for (&e_id, edge) in graph.edges.iter() {
            alphaloop_expression *= function!(
                vk_symbol!("uvprop"),
                &edge.momentum,
                function!(vk_symbol!("pow"), Atom::num(e_id as i64))
            );
        }

        Ok(Integral {
            name,
            n_loops: loop_mom_indices.len(),
            n_props: tot_n_props,
            generic_pattern,
            canonical_expression,
            parent_routing: None,
            short_expression,
            short_expression_pattern: Some(short_expression_pattern.to_pattern()),
            alphaloop_expression: Some(alphaloop_expression),
            applicable_evaluation_methods,
            graph,
            unoriented_generic_pattern,
            node_pairs,
        })
    }

    fn match_integral_to_short_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let unwrapped_input = input
            .replace(vk_parse!("topo(integral_)").unwrap().to_pattern())
            .with(vk_parse!("integral_").unwrap().to_pattern());
        if let Some(short_expression_pattern) = self.short_expression_pattern.as_ref() {
            // Contracted or fixed powers are literals, not bound wildcards.
            // A condition on an absent wildcard remains inconclusive in the
            // native matcher and cannot authenticate a complete match.
            let mut bound_wildcards = HashSet::new();
            short_expression_pattern.visitor(&mut |pattern| {
                if let Pattern::Wildcard(symbol, _) = pattern {
                    bound_wildcards.insert(*symbol);
                }
            });
            let power_conditions = (1..=self.n_props)
                .map(|i_prop| vk_symbol!(format!("pow{}_", i_prop)))
                .filter(|symbol| bound_wildcards.contains(symbol))
                .fold(Condition::True, |conditions, symbol| {
                    conditions & Condition::from((symbol, symbol_or_number()))
                });
            if let Some(m1) = unwrapped_input
                .pattern_match(short_expression_pattern, Some(&power_conditions), None)
                .next_detailed()
            {
                let mut replacement_rules = ReplacementRules::default();
                for i_prop in 1..=self.n_props {
                    if let Some(Match::Single(a)) = m1
                        .match_stack
                        .get(vk_symbol!(format!("pow{}_", i_prop).as_str()))
                    {
                        // NO: We do not want to match propagators with zero powers in the short form,
                        // as these should be matched to the pinched version with a hardcoded zero power
                        // ZERO_POWERS: Better to keep them to allow use to match to higher-level inputs
                        // if a.is_zero() {
                        //     return Ok(None);
                        // }
                        replacement_rules.canonical_expression_substitutions.insert(
                            vk_parse!(format!("pow({})", i_prop).as_str()).unwrap(),
                            a.to_owned(),
                        );
                    }
                    // Additive masses bind as a match sequence, so retain the
                    // complete image rather than only single-atom bindings.
                    if let Some(mass) = m1.match_stack.get(vk_symbol!(format!("msq{}_", i_prop))) {
                        replacement_rules.canonical_expression_substitutions.insert(
                            vk_parse!(format!("msq({})", i_prop).as_str()).unwrap(),
                            mass.to_atom(),
                        );
                    }
                }

                // Dummy substitutions for the numerator in this case
                for i_loop in 1..=self.n_loops {
                    replacement_rules.numerator_substitutions.insert(
                        vk_parse!(format!("k({},idx___)", i_loop).as_str()).unwrap(),
                        vk_parse!(format!("k({},idx___)", i_loop).as_str()).unwrap(),
                    );
                }
                Ok(Some(replacement_rules))
            } else {
                // Make sure the user did not intend to pass a short form expression
                if let Some(m) = input
                    .pattern_match(&vk_parse!("fn_(args__)").unwrap().to_pattern(), None, None)
                    .next_detailed()
                {
                    if let Match::FunctionName(s) = m.match_stack.get(vk_symbol!("fn_")).unwrap() {
                        if *s == vk_symbol!(self.name.clone()) {
                            Err(VakintError::InvalidShortExpression(format!("{}", input)))
                        } else {
                            Ok(None)
                        }
                    } else {
                        Ok(None)
                    }
                } else {
                    Ok(None)
                }
            }
        } else {
            Ok(None)
        }
    }

    fn match_integral_to_user_input(
        &self,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        // Check if the input is a short expression
        if let Ok(Some(a)) = self.match_integral_to_short_user_input(input) {
            return Ok(Some(a));
        }

        let undirected_input = input
            .replace(vk_parse!("edge(x_,y_)").unwrap().to_pattern())
            .with(vk_parse!("uedge(x_,y_)").unwrap().to_pattern());

        // Make sure that the unoriented topology could at least match, otherwise abort immediately
        if undirected_input
            .pattern_match(
                &self.unoriented_generic_pattern.pattern,
                Some(&self.unoriented_generic_pattern.conditions),
                None,
            )
            .next()
            .is_none()
        {
            return Ok(None);
        }

        // If multiple matches are possible we will score them and try and select the one that yields the simplest LMB correspondence
        // The second Vec<(usize,usize)> in the score is used to lift degeneracies and ensure fast matching.
        #[allow(clippy::type_complexity)]
        let mut best_match_thus_far: Option<(
            ((usize, usize), (Vec<usize>, Vec<(usize, bool)>)),
            ReplacementRules,
        )> = None;

        let mut topology_matcher = undirected_input.pattern_match(
            &self.generic_pattern.pattern,
            Some(&self.generic_pattern.conditions),
            Some(&self.generic_pattern.match_settings),
        );
        'outer: while let Some(m1) = topology_matcher.next_detailed() {
            // println!("VHDEBUG match result:");
            // for (k, v) in m1.match_stack.get_matches() {
            //     println!("{} -> {}", k, v.to_atom());
            // }
            // Make sure that all matched nodes are distinct

            let mut node_matches: HashMap<Symbol, HashSet<Atom>> = HashMap::default();
            for (unoriented_node, oriented_nodes) in self.node_pairs.iter() {
                node_matches
                    .entry(*unoriented_node)
                    .or_insert(HashSet::from_iter(
                        oriented_nodes
                            .iter()
                            .map(|&on| m1.match_stack.get(on).unwrap().to_atom()),
                    ));
            }

            // This check that all oriented node point to the same value
            if node_matches.values().any(|nodes| nodes.len() != 1) {
                continue 'outer;
            }
            // This check that all unoriented nodes point to different values
            if HashSet::<Atom>::from_iter(
                node_matches
                    .values()
                    .map(|nodes| nodes.iter().next().unwrap())
                    .cloned(),
            )
            .len()
                != node_matches.len()
            {
                // We do not need a continue here because there is not enough distinct nodes w.r.t canonical expression I believe
                break 'outer;
            }

            let mut replacement_rules = ReplacementRules::default();
            let mut score_for_this_match: (usize, usize) = (0, 0);
            let mut degeneracy_lifter_for_this_match: Vec<usize> = vec![];
            let mut degeneracy_lifter_exact_for_this_match: Vec<(usize, bool)> = vec![];
            for prop_id in 1..=self.n_props {
                if let Some(canonical_prop_match) = get_prop_with_id(
                    self.canonical_expression.as_ref().unwrap().as_view(),
                    prop_id,
                ) {
                    for var_prop_id in 1..=self.n_props {
                        if let Some(mtmp) = m1
                            .match_stack
                            .get(vk_symbol!(format!("pow{}_", var_prop_id).as_str()))
                            && let Entry::Vacant(e) = replacement_rules
                                .canonical_expression_substitutions
                                .entry(vk_parse!(format!("pow({})", var_prop_id).as_str()).unwrap())
                        {
                            e.insert(mtmp.to_atom());
                        }
                        if let Some(mtmp) = m1
                            .match_stack
                            .get(vk_symbol!(format!("msq{}_", var_prop_id)))
                            && let Entry::Vacant(e) = replacement_rules
                                .canonical_expression_substitutions
                                .entry(vk_parse!(format!("msq({})", var_prop_id).as_str()).unwrap())
                        {
                            e.insert(mtmp.to_atom());
                        }
                    }

                    let input_prop_id = if let Some(i) = get_integer_from_atom(
                        m1.match_stack
                            .get(vk_symbol!(format!("id{}_", prop_id).as_str()))
                            .unwrap()
                            .to_atom()
                            .as_view(),
                    ) {
                        i as usize
                    } else {
                        panic!("Match id from input expression must be an integer.")
                    };

                    // get the node ids in user's input for this prop as well as the one in the canonical expression
                    let input_prop_match = get_prop_with_id(input, input_prop_id).unwrap();
                    let (input_id_node_left, input_id_node_right) =
                        get_node_ids(&input_prop_match).unwrap();

                    let canonical_ids = get_node_ids(&canonical_prop_match).unwrap();

                    let (canonical_id_node_left, canonical_id_node_right) = (
                        get_integer_from_atom(
                            m1.match_stack
                                .get(vk_symbol!(format!("n{}l_", canonical_ids.0).as_str()))
                                .unwrap()
                                .to_atom()
                                .as_view(),
                        )
                        .unwrap() as usize,
                        get_integer_from_atom(
                            m1.match_stack
                                .get(vk_symbol!(format!("n{}r_", canonical_ids.1).as_str()))
                                .unwrap()
                                .to_atom()
                                .as_view(),
                        )
                        .unwrap() as usize,
                    );

                    let is_edge_flipped = if (input_id_node_left, input_id_node_right)
                        == (canonical_id_node_left, canonical_id_node_right)
                    {
                        false
                    } else if (input_id_node_right, input_id_node_left)
                        == (canonical_id_node_left, canonical_id_node_right)
                    {
                        true
                    } else {
                        unreachable!(
                            "Nodes IDs should have been matches: ({},{})!=({},{})=(n{}l_,n{}r_)",
                            input_id_node_left,
                            input_id_node_right,
                            canonical_id_node_left,
                            canonical_id_node_right,
                            canonical_ids.0,
                            canonical_ids.1
                        )
                    };

                    // To lift the degeneracy between matches, we store for each propagator the difference in id between the canonical expression
                    // and then the input expression
                    if is_edge_flipped {
                        score_for_this_match.1 += 1;
                    }
                    degeneracy_lifter_for_this_match
                        .push((input_prop_id as i32 - prop_id as i32).unsigned_abs() as usize);
                    degeneracy_lifter_exact_for_this_match.push((input_prop_id, !is_edge_flipped));
                    replacement_rules
                        .edge_ids_canonical_to_input_map
                        .insert(prop_id, input_prop_id);

                    let potential_lmb_match = input_prop_match.get(&vk_symbol!("q_")).unwrap();
                    // This is an LMB that'll need replacement in the numerator
                    // if let Match::Single(AtomView::Fun(_)) = potential_lmb_match {
                    if potential_lmb_match
                        .pattern_match(&vk_parse!("q_(args__)").unwrap().to_pattern(), None, None)
                        .count()
                        == 1
                    {
                        let (input_momentum_symbol, (input_momentum_atom_id, _input_momentum_id)) =
                            get_individual_momenta(potential_lmb_match.as_view())?
                                .pop()
                                .unwrap();
                        let input_lmb_pattern = function!(
                            input_momentum_symbol,
                            &input_momentum_atom_id,
                            &vk_parse!("idx___").unwrap()
                        );

                        let mut canonical_momenta_atom_for_pattern = canonical_prop_match
                            .get(&vk_symbol!("q_"))
                            .unwrap()
                            .replace(vk_parse!("k(ilmb_)").unwrap().to_pattern())
                            .when(Condition::from((vk_symbol!("ilmb_"), number_condition())))
                            .allow_new_wildcards_on_rhs(true)
                            .with(vk_parse!("k(ilmb_,idx___)").unwrap().to_pattern());
                        let lmb_replacement_length: i32 =
                            canonical_momenta_atom_for_pattern.nterms() as i32;
                        score_for_this_match.0 += (lmb_replacement_length - 1).max(0) as usize;
                        if is_edge_flipped {
                            canonical_momenta_atom_for_pattern =
                                (canonical_momenta_atom_for_pattern * -1).expand()
                        }

                        replacement_rules
                            .numerator_substitutions
                            .insert(input_lmb_pattern, canonical_momenta_atom_for_pattern);
                    }
                } else {
                    replacement_rules.canonical_expression_substitutions.insert(
                        vk_parse!(format!("pow({})", prop_id).as_str()).unwrap(),
                        Atom::Zero,
                    );
                }
            }

            if let Some(((score, degeneracy_lifer), _repl_rules)) = best_match_thus_far.as_ref() {
                // println!(
                //     "DEBUG A {:?}, {:?}, {:?}",
                //     score, degeneracy_lifer.0, degeneracy_lifer.1
                // );
                // println!(
                //     "DEBUG B {:?}, {:?}, {:?}",
                //     score_for_this_match,
                //     degeneracy_lifter_for_this_match,
                //     degeneracy_lifter_exact_for_this_match
                // );
                // Now check the degeneracy lifter
                if *score < score_for_this_match {
                    continue 'outer;
                } else {
                    match &(
                        degeneracy_lifter_for_this_match.clone(),
                        degeneracy_lifter_exact_for_this_match.clone(),
                    )
                        .cmp(degeneracy_lifer)
                    {
                        std::cmp::Ordering::Less => {}
                        std::cmp::Ordering::Greater => continue 'outer,
                        std::cmp::Ordering::Equal => continue 'outer,
                        //Somehow there are still degeneracies, but they should be irrelevant to the user at this stage.
                        //So this should not be put in the above: panic!("The degeneracy lifter should not be equal when matching topologies."),
                    }
                }
            }
            // println!(
            //     "DEBUG C {:?}, {:?}, {:?}",
            //     score_for_this_match,
            //     degeneracy_lifter_for_this_match,
            //     degeneracy_lifter_exact_for_this_match
            // );
            best_match_thus_far = Some((
                (
                    score_for_this_match,
                    (
                        degeneracy_lifter_for_this_match.clone(),
                        degeneracy_lifter_exact_for_this_match.clone(),
                    ),
                ),
                replacement_rules.clone(),
            ));

            // Early termination if this is a perfect match
            if score_for_this_match == (0, 0)
                && degeneracy_lifter_for_this_match.iter().all(|&x| x == 0)
            {
                break 'outer;
            }
        }

        if let Some((_score, replacement_rules)) = &best_match_thus_far {
            // println!("VHDEBUG FROM integral: {}", self);
            // println!("VHDEBUG Found match: {}", replacement_rules);
            // panic!("VHDEBUG STOP");
            Ok(Some(replacement_rules.to_owned()))
        } else {
            // let msg = [
            //     format!("User topology: {}", undirected_input),
            //     format!(
            //         "Unoriented pattern: {}",
            //         self.unoriented_generic_pattern.pattern.to_atom().unwrap()
            //     ),
            //     format!(
            //         "Oriented pattern: {}",
            //         self.generic_pattern.pattern.to_atom().unwrap()
            //     ),
            // ];
            // unreachable!("There should have been a match in Vakint at this stage. This is a logic error in vakint.\n{}",msg.join("\n"));

            Ok(None)
        }
    }

    fn to_canonical(&self, replacement_rules: &ReplacementRules, short_form: bool) -> Atom {
        let mut new_expression = if short_form {
            function!(
                vk_symbol!("topo"),
                self.short_expression.as_ref().unwrap().clone()
            )
        } else {
            self.canonical_expression.as_ref().unwrap().clone()
        };
        for (source, target) in replacement_rules.canonical_expression_substitutions.iter() {
            new_expression = new_expression
                .replace(source.to_pattern())
                .with(target.to_pattern());
        }

        // NO: Remove propagators with zero powers
        // ZERO_POWERS: Better to keep them to allow use to match to higher-level inputs
        // new_expression = vk_parse!("prop(propID_,edge(nl_,nr_),q_,mUVsq_,0)")
        //     .unwrap()
        //     .replace_all(
        //         new_expression.as_view(),
        //         &vk_parse!("1").unwrap().into(),
        //         None,
        //         None,
        //     );

        new_expression
    }
}

#[derive(Debug, Clone)]
pub struct PySecDecOptions {
    pub quiet: bool,
    pub relative_precision: f64,
    pub numerical_masses: HashMap<String, f64>,
    pub numerical_external_momenta: HashMap<String, (f64, f64, f64, f64)>,
    pub min_n_evals: u64,
    pub max_n_evals: u64,
    pub reuse_existing_output: Option<String>,
}

impl Default for PySecDecOptions {
    fn default() -> Self {
        // Give some random values to the external and masses. The user is expected to change these.
        let mut numerical_masses = HashMap::default();
        numerical_masses.insert(format!("{}::muvsq", NAMESPACE), 1.0);
        let mut numerical_external_momenta = HashMap::default();
        for i in 1..=10 {
            numerical_external_momenta.insert(format!("p{}", i), (13.0, 4.0, 3.0, 12.0));
        }
        PySecDecOptions {
            quiet: true,
            relative_precision: 1.0e-7,
            numerical_masses,
            numerical_external_momenta,
            min_n_evals: 10_000,
            max_n_evals: 1_000_000_000_000,
            reuse_existing_output: None,
        }
    }
}

impl fmt::Display for PySecDecOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "quiet={}, relative_precision={:.e}, min_n_evals={}, max_n_evals={}",
            self.quiet, self.relative_precision, self.min_n_evals, self.max_n_evals
        )
    }
}

#[derive(Debug, Clone)]
pub struct FMFTOptions {
    pub expand_masters: bool,
    pub susbstitute_masters: bool,
}

impl Default for FMFTOptions {
    fn default() -> Self {
        FMFTOptions {
            expand_masters: true,
            susbstitute_masters: true,
        }
    }
}

impl fmt::Display for FMFTOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expand_masters={}, susbstitute_masters={}",
            self.expand_masters, self.susbstitute_masters
        )
    }
}

#[derive(Debug, Clone)]
pub struct AlphaLoopOptions {
    pub susbstitute_masters: bool,
}

impl fmt::Display for AlphaLoopOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "susbstitute_masters={}", self.susbstitute_masters,)
    }
}

impl Default for AlphaLoopOptions {
    fn default() -> Self {
        AlphaLoopOptions {
            susbstitute_masters: true,
        }
    }
}
#[derive(Debug, Clone)]
pub struct MATADOptions {
    pub expand_masters: bool,
    pub susbstitute_masters: bool,
    pub substitute_hpls: bool,
    pub direct_numerical_substition: bool,
}

impl Default for MATADOptions {
    fn default() -> Self {
        MATADOptions {
            expand_masters: true,
            susbstitute_masters: true,
            substitute_hpls: true,
            direct_numerical_substition: true,
        }
    }
}

impl fmt::Display for MATADOptions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expand_masters={}, susbstitute_masters={}, substitute_hpls={}, direct_numerical_substition={}",
            self.expand_masters,
            self.susbstitute_masters,
            self.substitute_hpls,
            self.direct_numerical_substition
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VakintDependency {
    FORM,
    PySecDec,
}

#[derive(Debug, Clone)]
pub enum EvaluationMethod {
    AlphaLoop(AlphaLoopOptions),
    MATAD(MATADOptions),
    FMFT(FMFTOptions),
    PySecDec(PySecDecOptions),
    /// FORM-independent scalar reduction through shipped RustRed IBP artifacts.
    RustRed(RustRedEvaluationOptions),
}

impl fmt::Display for EvaluationMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EvaluationMethod::AlphaLoop(opts) => write!(f, "AlphaLoop ({})", opts),
            EvaluationMethod::MATAD(opts) => write!(f, "MATAD ({})", opts),
            EvaluationMethod::FMFT(opts) => write!(f, "FMFT ({})", opts),
            EvaluationMethod::PySecDec(opts) => {
                write!(f, "PySecDec ({})", opts)
                //write!(f, "PySecDec")
            }
            EvaluationMethod::RustRed(opts) => write!(f, "RustRed ({})", opts),
        }
    }
}

impl EvaluationMethod {
    fn numerical_output(&self) -> bool {
        match self {
            Self::AlphaLoop(options) => options.susbstitute_masters,
            Self::MATAD(options) => {
                options.expand_masters
                    && (options.direct_numerical_substition
                        || (options.susbstitute_masters && options.substitute_hpls))
            }
            Self::FMFT(options) => options.expand_masters && options.susbstitute_masters,
            Self::PySecDec(_) => true,
            Self::RustRed(options) => options.master_mode == RustRedMasterMode::Numerical,
        }
    }

    fn stages_input_coefficients(&self) -> bool {
        matches!(self, Self::RustRed(_))
            || (self.numerical_output() && !matches!(self, Self::PySecDec(_)))
    }
    pub fn supports(&self, settings: &VakintSettings, topology: &Topology) -> bool {
        match self {
            EvaluationMethod::AlphaLoop(_) => {
                topology
                    .get_integral()
                    .applicable_evaluation_methods
                    .0
                    .iter()
                    .any(|m| matches!(m, EvaluationMethod::AlphaLoop(_)))
                    && topology.get_integral().alphaloop_expression.is_some()
                    && settings.number_of_terms_in_epsilon_expansion <= 4
                    && topology.get_integral().n_loops <= 3
            }
            EvaluationMethod::MATAD(_) => {
                topology
                    .get_integral()
                    .applicable_evaluation_methods
                    .0
                    .iter()
                    .any(|m| matches!(m, EvaluationMethod::MATAD(_)))
                    && topology.get_integral().n_loops <= 3
                    && settings.number_of_terms_in_epsilon_expansion <= 5
            }
            EvaluationMethod::FMFT(_) => {
                topology
                    .get_integral()
                    .applicable_evaluation_methods
                    .0
                    .iter()
                    .any(|m| matches!(m, EvaluationMethod::FMFT(_)))
                    && topology.get_integral().n_loops == 4
                    && settings.number_of_terms_in_epsilon_expansion <= 5
            }
            EvaluationMethod::PySecDec(_) => topology
                .get_integral()
                .applicable_evaluation_methods
                .0
                .iter()
                .any(|m| matches!(m, EvaluationMethod::PySecDec(_))),
            EvaluationMethod::RustRed(options) => {
                rustred_evaluation::supports(settings, topology, options)
            }
        }
    }

    pub fn dependencies(&self) -> Vec<VakintDependency> {
        match self {
            EvaluationMethod::AlphaLoop(_) => {
                vec![VakintDependency::FORM]
            }
            EvaluationMethod::MATAD(_) => {
                vec![VakintDependency::FORM]
            }
            EvaluationMethod::FMFT(_) => {
                vec![VakintDependency::FORM]
            }
            EvaluationMethod::PySecDec(_) => {
                vec![VakintDependency::FORM, VakintDependency::PySecDec]
            }
            EvaluationMethod::RustRed(_) => vec![],
        }
    }

    pub fn adjust(
        &mut self,
        quiet: Option<bool>,
        relative_precision: f64,
        numerical_params_real: &HashMap<String, Float, RandomState>,
        numerical_params_complex: &HashMap<String, Complex<Float>, RandomState>,
        numerical_external_momenta: &HashMap<usize, Momentum, RandomState>,
    ) {
        match self {
            EvaluationMethod::AlphaLoop(_) => {}
            EvaluationMethod::MATAD(_) => {}
            EvaluationMethod::FMFT(_) => {}
            EvaluationMethod::RustRed(_) => {}
            EvaluationMethod::PySecDec(opts) => {
                let mut f64_numerical_masses: HashMap<String, f64, std::hash::RandomState> =
                    numerical_params_complex
                        .iter()
                        .map(|(k, v)| (k.clone(), v.norm().re.to_f64()))
                        .collect();
                for (k, v) in numerical_params_real.iter() {
                    f64_numerical_masses.insert(k.clone(), v.to_f64());
                }
                let f64_numerical_external_momenta = numerical_external_momenta
                    .iter()
                    .map(|(k, v)| {
                        (
                            format!("p{}", k),
                            match v {
                                Momentum::Complex(k) => (
                                    k.0.norm().re.to_f64(),
                                    k.1.norm().re.to_f64(),
                                    k.2.norm().re.to_f64(),
                                    k.3.norm().re.to_f64(),
                                ),
                                Momentum::Real(k) => {
                                    (k.0.to_f64(), k.1.to_f64(), k.2.to_f64(), k.3.to_f64())
                                }
                            },
                        )
                    })
                    .collect();
                if let Some(is_quiet) = quiet {
                    opts.quiet = is_quiet;
                }
                opts.relative_precision = relative_precision;
                opts.numerical_masses = f64_numerical_masses;
                opts.numerical_external_momenta = f64_numerical_external_momenta;
            }
        }
    }

    pub fn evaluate_integral(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
    ) -> Result<Atom, VakintError> {
        let result = self.reduce_integral(vakint, settings, numerator, integral_specs, 0)?;
        if self.numerical_output() {
            vakint.resolve_rustred_numerical(result, settings)
        } else {
            Ok(result)
        }
    }

    fn reduce_integral(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
        source_id: usize,
    ) -> Result<Atom, VakintError> {
        let staged = if self.stages_input_coefficients() {
            Vakint::stage_input_coefficients(settings, numerator)?
        } else {
            numerator.to_owned()
        };
        let numerator = staged.as_view();
        let result = match self {
            EvaluationMethod::AlphaLoop(opts) => {
                vakint.alphaloop_evaluate(settings, numerator, integral_specs, opts)
            }
            EvaluationMethod::MATAD(opts) => {
                vakint.matad_reduce(settings, numerator, integral_specs, opts)
            }
            EvaluationMethod::FMFT(opts) => {
                vakint.fmft_reduce(settings, numerator, integral_specs, opts)
            }
            EvaluationMethod::PySecDec(opts) => {
                vakint.pysecdec_evaluate(settings, numerator, integral_specs, opts, source_id)
            }
            EvaluationMethod::RustRed(opts) => {
                vakint.rustred_evaluate(settings, numerator, integral_specs, opts)
            }
        }?;
        // Simplify logarithms and zero powers knowing that all arguments are real
        Ok(simplify_real(result.as_view()))
    }
}

#[derive(Debug, Clone)]
pub struct EvaluationOrder(
    /// Evaluation backends tried in order until one supports the requested integral.
    pub Vec<EvaluationMethod>,
);

impl Default for EvaluationOrder {
    fn default() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop(AlphaLoopOptions::default()),
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::FMFT(FMFTOptions::default()),
            EvaluationMethod::PySecDec(PySecDecOptions::default()),
        ])
    }
}

impl fmt::Display for EvaluationOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}]",
            self.0
                .iter()
                .map(|m| format!("{}", m))
                .collect::<Vec<_>>()
                .join("|")
        )
    }
}

impl EvaluationOrder {
    pub fn numerical_only() -> Self {
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions::default())])
    }
    pub fn analytic_only() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop(AlphaLoopOptions::default()),
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::FMFT(FMFTOptions::default()),
        ])
    }
    pub fn alphaloop_only() -> Self {
        EvaluationOrder(vec![EvaluationMethod::AlphaLoop(
            AlphaLoopOptions::default(),
        )])
    }
    pub fn matad_only(matad_options: Option<MATADOptions>) -> Self {
        EvaluationOrder(vec![EvaluationMethod::MATAD(
            matad_options.unwrap_or_default(),
        )])
    }
    /// Select only the FORM-independent RustRed scalar evaluator.
    pub fn rustred_only() -> Self {
        EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions::default(),
        )])
    }
    pub fn fmft_only(fmft_options: Option<FMFTOptions>) -> Self {
        EvaluationOrder(vec![EvaluationMethod::FMFT(
            fmft_options.unwrap_or_default(),
        )])
    }
    pub fn fmft_and_numerical(
        fmft_options: Option<FMFTOptions>,
        pysecdec_options: Option<PySecDecOptions>,
    ) -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::FMFT(fmft_options.unwrap_or_default()),
            EvaluationMethod::PySecDec(pysecdec_options.unwrap_or_default()),
        ])
    }
    pub fn pysecdec_only(pysecdec_options: Option<PySecDecOptions>) -> Self {
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            pysecdec_options.unwrap_or_default(),
        )])
    }
    pub fn empty() -> Self {
        EvaluationOrder(vec![])
    }
    pub fn all() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop(AlphaLoopOptions::default()),
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::FMFT(FMFTOptions::default()),
            EvaluationMethod::PySecDec(PySecDecOptions::default()),
        ])
    }
    pub fn all_but_fmft() -> Self {
        EvaluationOrder(vec![
            EvaluationMethod::AlphaLoop(AlphaLoopOptions::default()),
            EvaluationMethod::MATAD(MATADOptions::default()),
            EvaluationMethod::PySecDec(PySecDecOptions::default()),
        ])
    }

    pub fn adjust(
        &mut self,
        quiet: Option<bool>,
        relative_precision: f64,
        numerical_params_real: &HashMap<String, Float, RandomState>,
        numerical_params_complex: &HashMap<String, Complex<Float>, RandomState>,
        numerical_external_momenta: &HashMap<usize, Momentum, RandomState>,
    ) {
        for method in self.0.iter_mut() {
            method.adjust(
                quiet,
                relative_precision,
                numerical_params_real,
                numerical_params_complex,
                numerical_external_momenta,
            );
        }
    }
}

#[derive(Debug, Clone)]
pub enum InputFloatRationalizationPrecision {
    FullPrecision,
    TargetPrecision,
}

/// Strategy used to reduce tensor numerators before scalar-integral evaluation.
///
/// `AlphaLoop` invokes the historical FORM projector bundled with Vakint.
/// `FeynKit` uses FeynKit's native symmetry-aware projector and therefore does
/// not require FORM for numerator reduction.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum TensorReductionMethod {
    AlphaLoop,
    #[default]
    FeynKit,
}

impl fmt::Display for TensorReductionMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AlphaLoop => write!(f, "alphaloop"),
            Self::FeynKit => write!(f, "feynkit"),
        }
    }
}

impl std::str::FromStr for TensorReductionMethod {
    type Err = VakintError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "alphaloop" => Ok(Self::AlphaLoop),
            "feynkit" => Ok(Self::FeynKit),
            _ => Err(VakintError::InvalidTensorReductionMethod(value.to_owned())),
        }
    }
}

#[derive(Debug, Clone)]
pub struct VakintSettings {
    #[allow(unused)]
    /// Symbol used for the dimensional-regularization parameter epsilon.
    pub epsilon_symbol: String,
    /// Symbol used for the squared renormalization scale.
    pub mu_r_sq_symbol: String,
    /// Executable name or path used to invoke FORM.
    pub form_exe_path: String,
    /// Executable name or path used to invoke Python-based backends.
    pub python_exe_path: String,
    /// Whether canonicalization must account for every numerator factor.
    pub verify_numerator_identification: bool,
    /// Per-loop normalization convention applied to evaluated integrals.
    pub integral_normalization_factor: LoopNormalizationFactor,
    /// Requested arithmetic precision, not a guarantee of master-data accuracy.
    /// Known lower-precision master sources produce a warning; evaluation keeps
    /// this working precision without gaining additional accurate source digits.
    /// Missing master data and unsupported Laurent orders remain errors.
    pub run_time_decimal_precision: u32,
    /// Whether inputs may match Vakint's generic unknown-integral topology.
    pub allow_unknown_integrals: bool,
    /// Whether successful external-backend runs remove their temporary directories.
    pub clean_tmp_dir: bool,
    /// Ordered analytic and numerical backends considered for each integral.
    pub evaluation_order: EvaluationOrder,
    /// Optional externally supplied application package for symbolic five-loop RustRed.
    /// No package is downloaded or embedded; lower-loop artifacts are unchanged.
    pub rustred_artifact_path: Option<PathBuf>,
    /// Cooperative process-RSS and elapsed-time limits shared by one request.
    /// Absent ceilings impose no restriction. Opaque external calls can only
    /// be checked at their boundaries; hard deadlines require a supervisor.
    pub rustred_resources: ResourceLimits,
    /// Operational application/cache allowances; these never enlarge saved entry domains.
    pub rustred_application_options: rustred_app::SavedApplicationOptions,
    /// Per-call numerator lowering allowances for the saved five-loop lane.
    /// These do not change entry admission, legacy lower-loop programs, or
    /// the independent lower-factor product policy.
    pub rustred_scalar_numerator_limits: RustRedScalarNumeratorLimits,
    /// Separate lower-factor product application/cache allowance, not an
    /// aggregate with `rustred_application_options`.
    pub rustred_product_limits:
        rustred::reduction::terminal_normalization::VacuumProductReductionLimits,
    /// Bounded exact alternatives used to obtain an admitted five-loop routing.
    pub rustred_momentum_match_limits: topologies::MomentumMatchLimits,
    pub tensor_reduction_method: TensorReductionMethod,
    // This quantity is typically set equal to *one plus the maximum loop count* of the UV regularisation problem considered.
    // For example when considering a 2-loop problem, then:
    //   a) for the nested one-loop integrals appearing, the single pole, finite term *and* order-epsilon term will need to be considered.
    //   b) for the two-loop integrals, the double pole, single pole and finite terms will be needed, so again three terms
    /// Number of terms retained in the dimensional-regulator expansion.
    pub number_of_terms_in_epsilon_expansion: i64,
    /// Precision policy used when rationalizing floating-point input for FORM.
    pub precision_for_input_float_rationalization: InputFloatRationalizationPrecision,
    /// Whether reduced scalar products remain in dot-product notation.
    pub use_dot_product_notation: bool,
    /// Optional parent directory for uniquely named external-backend work directories.
    pub temporary_directory: Option<String>,
}

impl VakintSettings {
    pub(crate) fn with_resource_budget<T>(
        &self,
        operation: impl FnOnce() -> Result<T, VakintError>,
    ) -> Result<T, VakintError> {
        ResourceBudget::new(self.rustred_resources).run(|| {
            resources::checkpoint()?;
            let result = operation()?;
            resources::checkpoint()?;
            Ok(result)
        })
    }

    pub fn get_integral_normalization_factor_atom(&self) -> Result<Atom, VakintError> {
        self.integral_normalization_factor.to_atom(self)
    }

    /// Retain finite normalization inputs before backend algebra mixes orders.
    fn get_staged_integral_normalization_factor_atom(&self) -> Result<Atom, VakintError> {
        precision::legacy::stage_table(
            Atom::var(vk_symbol!("normalization_input")).as_view(),
            self.get_integral_normalization_factor_atom()?.as_view(),
            vk_symbol!(&self.epsilon_symbol),
            "Input",
        )
    }

    pub fn get_binary_precision(&self) -> u32 {
        ((self.run_time_decimal_precision.max(17) as f64) * LOG2_10).floor() as u32
    }

    pub fn real_to_prec(&self, re: &str) -> Float {
        let prec = self.get_binary_precision();
        Float::parse(re, Some(prec)).unwrap()
    }

    pub fn complex_to_prec(&self, re: &str, im: &str) -> Complex<Float> {
        let prec = self.get_binary_precision();
        Complex::new(
            Float::parse(re, Some(prec)).unwrap(),
            Float::parse(im, Some(prec)).unwrap(),
        )
    }
}

#[derive(Debug, Clone)]
pub enum LoopNormalizationFactor {
    #[allow(non_camel_case_types)]
    pySecDec,
    MSbar,
    FMFTandMATAD,
    Custom(String),
}

impl LoopNormalizationFactor {
    pub fn to_expression(&self) -> String {
        match self {
            LoopNormalizationFactor::pySecDec => "(𝑖*(𝜋^((4-2*eps)/2)))^(-n_loops)".into(),
            LoopNormalizationFactor::FMFTandMATAD => {
                "( 𝑖*(𝜋^((4-2*eps)/2)) * (exp(-EulerGamma))^(eps) )^(-n_loops)".into()
            }
            LoopNormalizationFactor::MSbar => {
                // We must include the 1/(2*𝜋)^D factor per loop which accompanies the text-book definition of MSbar
                // "(2*𝜋)^(-n_loops)*(2*𝜋)^(2*eps*n_loops)*(exp(log_mu_sq)/(4*𝜋*exp(-EulerGamma)))^(eps*n_loops)".into()
                // Which simplifies into the following expression
                //"(2*𝜋)^(-4*n_loops)*(exp(log_mu_sq)*𝜋*exp(EulerGamma))^(eps*n_loops)".into()
                // And where it is best to keep each term taken to an epsilon power separately so that after the expansion
                // we can easily get the expected cancellation of log(4 pi) and eulerGamma.
                "(2*𝜋)^(-4*n_loops)*exp(log_mu_sq)^(eps*n_loops)*𝜋^(eps*n_loops)*exp(EulerGamma)^(eps*n_loops)".into()
            }
            LoopNormalizationFactor::Custom(s) => s.clone(),
        }
    }

    pub fn static_allowed_symbols() -> Vec<String> {
        vec![
            "eps".into(),
            "log_mu_sq".into(),
            "EulerGamma".into(),
            "𝑖".into(),
            "I".into(),
            "𝜋".into(),
            "pi".into(),
        ]
    }

    pub fn allowed_symbols(settings: &VakintSettings) -> Vec<String> {
        let mut allowed_symbols = LoopNormalizationFactor::static_allowed_symbols();
        allowed_symbols.push(settings.epsilon_symbol.clone());
        allowed_symbols
    }

    pub fn to_atom(&self, settings: &VakintSettings) -> Result<Atom, VakintError> {
        let mut a = Atom::try_from(self)?;
        a = a
            .replace(vk_parse!("eps").unwrap().to_pattern())
            .with(vk_parse!(&settings.epsilon_symbol).unwrap().to_pattern());
        Ok(a)
    }

    pub fn validate(
        &self,
        settings: &VakintSettings,
    ) -> Result<(Atom, Series<AtomField>, NumericalEvaluationResult), VakintError> {
        let expr: Atom = self.to_atom(settings)?;
        let expanded_expr = expr
            .replace(vk_parse!("n_loops").unwrap().to_pattern())
            .with(Atom::num(1).to_pattern());

        let expanded_expr = match expanded_expr.series(
            vk_symbol!(settings.epsilon_symbol.as_str()),
            Atom::Zero.as_atom_view(),
            Rational::from(settings.number_of_terms_in_epsilon_expansion - 1),
        ) {
            Ok(a) => a,
            Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
        };

        let mut expanded_expr_atom = expanded_expr.to_atom();
        let log_mu_sq = function!(
            Symbol::LOG,
            Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );
        expanded_expr_atom = expanded_expr_atom
            .replace(vk_parse!("log_mu_sq").unwrap().to_pattern())
            .with(log_mu_sq.to_pattern());

        let mut params: HashMap<String, Float, _> = HashMap::default();
        params.insert(
            get_full_name(&vk_symbol!(settings.mu_r_sq_symbol.as_str())),
            settings.real_to_prec("1"),
        );
        let num_res = match Vakint::full_numerical_evaluation_without_error(
            settings,
            expanded_expr_atom.as_view(),
            &params,
            &HashMap::default(),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                return Err(VakintError::InvalidLoopNormalization(
                    expr.to_canonical_string(),
                    e.to_string(),
                    LoopNormalizationFactor::allowed_symbols(settings).join(","),
                ));
            }
        };
        Ok((expr, expanded_expr, num_res))
    }
}

impl fmt::Display for LoopNormalizationFactor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LoopNormalizationFactor::pySecDec => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("PySecDec").green(),
                    self.to_expression()
                )
            }
            LoopNormalizationFactor::MSbar => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("MSbar").green(),
                    self.to_expression()
                )
            }
            LoopNormalizationFactor::FMFTandMATAD => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("FMFTandMATAD").green(),
                    self.to_expression()
                )
            }
            LoopNormalizationFactor::Custom(_) => {
                write!(
                    f,
                    "{} convention [ {} ]",
                    String::from("Custom").green(),
                    self.to_expression()
                )
            }
        }
    }
}

impl TryFrom<&LoopNormalizationFactor> for Atom {
    type Error = VakintError;
    fn try_from(expr: &LoopNormalizationFactor) -> Result<Self, VakintError> {
        match vk_parse!(&expr.to_expression()) {
            Ok(a) => {
                let mut processed_a = a
                    .replace(vk_parse!("I").unwrap().to_pattern())
                    .with(Atom::var(S.cmplx_i).to_pattern());
                processed_a = processed_a
                    .replace(vk_parse!("pi").unwrap().to_pattern())
                    .with(Atom::var(Symbol::PI).to_pattern());
                Ok(processed_a)
            }
            Err(e) => Err(VakintError::InvalidLoopNormalization(
                format!("{}", expr),
                e.to_string(),
                LoopNormalizationFactor::static_allowed_symbols().join(","),
            )),
        }
    }
}

#[allow(clippy::derivable_impls)]
impl Default for VakintSettings {
    fn default() -> Self {
        VakintSettings {
            epsilon_symbol: format!("{}::ε", NAMESPACE),
            mu_r_sq_symbol: format!("{}::mursq", NAMESPACE),
            form_exe_path: env::var("FORM_PATH").unwrap_or("form".into()),
            python_exe_path: env::var("PYTHON_BIN_PATH").unwrap_or("python3".into()),
            verify_numerator_identification: false,
            run_time_decimal_precision: 32,
            integral_normalization_factor: LoopNormalizationFactor::pySecDec,
            allow_unknown_integrals: true,
            clean_tmp_dir: env::var("VAKINT_NO_CLEAN_TMP_DIR").is_err(),
            evaluation_order: EvaluationOrder::default(),
            rustred_artifact_path: None,
            rustred_resources: ResourceLimits::default(),
            rustred_application_options: rustred_app::SavedApplicationOptions::default(),
            rustred_scalar_numerator_limits: Default::default(),
            rustred_product_limits: Default::default(),
            rustred_momentum_match_limits: Default::default(),
            tensor_reduction_method: TensorReductionMethod::FeynKit,
            // Default to a three-loop UV subtraction problem, for which alphaLoop implementation can be used.
            number_of_terms_in_epsilon_expansion: 4,
            precision_for_input_float_rationalization:
                InputFloatRationalizationPrecision::FullPrecision,
            use_dot_product_notation: false,
            temporary_directory: None,
        }
    }
}

#[derive(Debug, Clone)]
struct FullPattern {
    pattern: Pattern,
    conditions: Condition<PatternRestriction>,
    match_settings: MatchSettings,
}

impl From<Pattern> for FullPattern {
    fn from(pattern: Pattern) -> FullPattern {
        FullPattern {
            pattern,
            conditions: Condition::default(),
            match_settings: MatchSettings::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Vakint {
    /// Base/custom matching registry. Successful saved-artifact loading adds a
    /// shared immutable five-loop catalogue without replacing the base entries.
    pub topologies: Topologies,
    rustred_application:
        Arc<std::sync::Mutex<Option<rustred_evaluation::five_loop::LoadedApplication>>>,
    rustred_lower_caches: Arc<std::sync::Mutex<rustred_evaluation::cache::LoadedLowerCaches>>,
    rustred_data: Option<Arc<rustred_evaluation::data::RustRedData>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VakintTerm {
    pub integral: Atom,
    pub numerator: Atom,
    pub vectors: Vec<(String, i64)>,
}

impl From<VakintTerm> for Atom {
    fn from(vakint_term: VakintTerm) -> Atom {
        vakint_term.integral * vakint_term.numerator
    }
}

impl fmt::Display for VakintTerm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({}) x {}",
            format!("{:>}", self.numerator).cyan(),
            format!("{}", self.integral).green(),
        )
    }
}

enum IntegralEvaluation {
    Complete {
        numerical: bool,
    },
    RustRed {
        loop_count: usize,
        mode: RustRedMasterMode,
    },
}

impl VakintTerm {
    fn tracks_input_coefficients(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<bool, VakintError> {
        let Some(mut matched) =
            vakint.match_topologies_for_settings(settings, self.integral.as_view())?
        else {
            return Ok(false);
        };
        // Admission checks inspect concrete powers and masses, just as they do
        // at reduction. The unmatched template still contains pow(n) slots.
        matched.apply_replacement_rules()?;
        Ok(settings
            .evaluation_order
            .0
            .iter()
            .find(|method| method.supports(settings, &matched.canonical_topology))
            .is_some_and(EvaluationMethod::stages_input_coefficients))
    }

    fn stage_numerical_inputs(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        let mut has_float = false;
        self.numerator.visitor(&mut |part| {
            if matches!(part, AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Float(_, _))) {
                has_float = true;
            }
            !has_float
        });
        if !has_float
            || !settings
                .evaluation_order
                .0
                .iter()
                .any(EvaluationMethod::stages_input_coefficients)
        {
            return Ok(());
        }
        if self.tracks_input_coefficients(vakint, settings)? {
            self.numerator = Vakint::stage_input_coefficients(settings, self.numerator.as_view())?;
        }
        Ok(())
    }

    pub fn evaluate_integral(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        settings.with_resource_budget(|| {
            let numerical = match self.reduce_integral(vakint, settings, 0, None)? {
                IntegralEvaluation::RustRed { loop_count, mode } => {
                    self.numerator = vakint.finalize_rustred(
                        self.numerator.clone(),
                        settings,
                        loop_count,
                        mode,
                        None,
                    )?;
                    mode == RustRedMasterMode::Numerical
                }
                IntegralEvaluation::Complete { numerical } => numerical,
            };
            if numerical {
                self.numerator =
                    vakint.resolve_rustred_numerical(self.numerator.clone(), settings)?;
            }
            Ok(())
        })
    }

    fn reduce_integral(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
        source_id: usize,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<IntegralEvaluation, VakintError> {
        MonitorContext::native_phase(observation.as_deref_mut(), "topology_matching", None);
        self.stage_numerical_inputs(vakint, settings)?;
        let mut integral_specs = if let Some(replacement_rules) =
            vakint.match_topologies_for_settings(settings, self.integral.as_view())?
        {
            replacement_rules
        } else {
            return Err(VakintError::UnreckognizedIntegral(
                self.integral.to_string(),
            ));
        };

        // Choose the witness before transporting the scalar numerator. A
        // finite saved union may admit only some exact family presentations;
        // ordinary/default backends keep their original matching path.
        let scalar_rustred =
            matches!(
                settings.evaluation_order.0.iter().find(|method| {
                    method.supports(settings, &integral_specs.canonical_topology)
                }),
                Some(EvaluationMethod::RustRed(_))
            );
        let scalar_five_loop =
            matches!(integral_specs.canonical_topology, Topology::FiveLoop(_)) && scalar_rustred;
        if scalar_five_loop {
            MonitorContext::native_phase(observation.as_deref_mut(), "routing_verification", None);
            integral_specs = vakint
                .match_rustred_input(settings, self.integral.as_view(), self.numerator.as_view())?
                .unwrap_or(integral_specs);
        }

        integral_specs.apply_replacement_rules()?;
        if scalar_rustred {
            self.apply_scalar_numerator_replacement_rules(&integral_specs, settings)?;
        } else {
            self.apply_numerator_replacement_rules(&integral_specs, settings)?;
        }

        for evaluation_approach in &settings.evaluation_order.0 {
            if evaluation_approach.supports(settings, &integral_specs.canonical_topology) {
                let (evaluated_integral, stage) = match evaluation_approach {
                    EvaluationMethod::RustRed(options) => (
                        vakint.rustred_reduce_observed(
                            settings,
                            self.numerator.as_view(),
                            &integral_specs,
                            observation.as_deref_mut(),
                        )?,
                        IntegralEvaluation::RustRed {
                            loop_count: integral_specs.canonical_topology.get_integral().n_loops,
                            mode: options.master_mode,
                        },
                    ),
                    _ => (
                        evaluation_approach.reduce_integral(
                            vakint,
                            settings,
                            self.numerator.as_view(),
                            &integral_specs,
                            source_id,
                        )?,
                        IntegralEvaluation::Complete {
                            numerical: evaluation_approach.numerical_output(),
                        },
                    ),
                };
                self.numerator = simplify_real(evaluated_integral.as_view());
                self.integral = Atom::num(1);
                return Ok(stage);
            }
        }
        Err(VakintError::NoEvaluationMethodFound(
            self.integral.to_string(),
            settings.number_of_terms_in_epsilon_expansion
                - integral_specs.canonical_topology.get_integral().n_loops as i64
                - 1,
        ))
    }

    pub fn apply_numerator_replacement_rules(
        &mut self,
        replacement_rules: &ReplacementRules,
        vakint_settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        let mut new_numerator = self.numerator.clone();
        let mut test = self.numerator.clone();
        // Here we must be careful that substituations can be of the form:
        // {a-> b, b->c}
        // so that we must not apply the substitutions to the outcome of the previous substitution.
        // We will use replace_all_multiple for this
        let casted_replacement_rules = replacement_rules
            .numerator_substitutions
            .iter()
            .map(|(source, target)| {
                (
                    source.clone().to_pattern(),
                    target.clone().as_view().to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let one_substitution_pattern = vk_parse!("1").unwrap().to_pattern();
        new_numerator = new_numerator.replace_multiple(
            casted_replacement_rules
                .iter()
                .map(|(source, target)| Replacement::new(source.to_owned(), target.to_owned()))
                .collect::<Vec<_>>()
                .as_slice(),
        );

        // If there is no casted_replacement_rules, this means that there was no canonization rule for the numerator, so we must manually set loop momenta to 1.
        // if casted_replacement_rules.len() == 0 {
        //     test = test
        //         .replace(
        //             vk_parse!(&format!("{}(momID_,idx___)", LOOP_MOMENTUM_SYMBOL))
        //                 .unwrap()
        //                 .to_pattern(),
        //         )
        //         .with(one_substitution_pattern.clone());
        // } else {
        //     test = test.replace_multiple(
        //         casted_replacement_rules
        //             .iter()
        //             .map(|(source, _target)| {
        //                 Replacement::new(source.to_owned(), one_substitution_pattern.clone())
        //             })
        //             .collect::<Vec<_>>()
        //             .as_slice(),
        //     );
        // }
        test = test.replace_multiple(
            casted_replacement_rules
                .iter()
                .map(|(source, _target)| {
                    Replacement::new(source.to_owned(), one_substitution_pattern.clone())
                })
                .collect::<Vec<_>>()
                .as_slice(),
        );

        // Make sure to also set all externals to zero for the test
        test = test
            .replace(
                vk_parse!(format!("{}(momID_,idx___)", EXTERNAL_MOMENTUM_SYMBOL).as_str())
                    .unwrap()
                    .to_pattern(),
            )
            .with(vk_parse!("1").unwrap().to_pattern());

        // Substitute metric as well
        test = test
            .replace(
                vk_parse!(format!("{}(idx1_,idx2_)", METRIC_SYMBOL))
                    .unwrap()
                    .to_pattern(),
            )
            .with(vk_parse!("1").unwrap().to_pattern());

        // Substitute dot product wrappers as well, with the understanding that all momenta in arguments should now have been mapped to 1 at this point.
        test = test
            .replace(
                vk_parse!(format!("{}(1,1)", DOT_SYMBOL))
                    .unwrap()
                    .to_pattern(),
            )
            .with(vk_parse!("1").unwrap().to_pattern());

        // Substitute epsilon regulator
        test = test
            .replace(
                vk_parse!(&vakint_settings.epsilon_symbol)
                    .unwrap()
                    .to_pattern(),
            )
            .with(vk_parse!("1").unwrap().to_pattern());

        if vakint_settings.verify_numerator_identification {
            // Matched propagators authenticate scalar mass parameters; retain
            // their arithmetic symbolically in this verification copy.
            let mut mass_parameters = HashSet::new();
            let masses = replacement_rules
                .canonical_topology
                .get_integral()
                .canonical_expression
                .as_ref()
                .map(|_| replacement_rules.get_propagator_property_list("mUVsq_"))
                .unwrap_or_default();
            let mass_substitutions = replacement_rules
                .canonical_expression_substitutions
                .iter()
                .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
                .collect::<Vec<_>>();
            for mass in masses.values() {
                // Canonicalization can verify against an uninstantiated
                // topology template; evaluation already applied these rules.
                let mass = mass.replace_multiple(mass_substitutions.as_slice());
                if mass.get_all_symbols(true).iter().any(|symbol| {
                    [S.k, S.p, S.dot, S.g].contains(symbol)
                        || symbol.get_name().starts_with("spenso::")
                }) {
                    return Err(VakintError::InvalidIntegralFormat(format!(
                        "matched propagator mass contains momentum or tensor notation: {mass}"
                    )));
                }
                let mut algebraic = true;
                mass.visitor(&mut |part| {
                    algebraic &= !matches!(part, AtomView::Fun(_));
                    algebraic
                });
                let parameters = mass.get_all_symbols(false);
                if algebraic {
                    mass_parameters.extend(parameters);
                }
            }
            // A vanishing mass prefactor must not hide an unmatched momentum
            // or scalar coefficient. Check structure without assigning masses
            // arbitrary values, which could also introduce spurious poles.
            let mut has_function = false;
            test.visitor(&mut |part| {
                has_function |= matches!(part, AtomView::Fun(_));
                !has_function
            });
            if has_function
                || test
                    .get_all_symbols(true)
                    .iter()
                    .any(|symbol| !mass_parameters.contains(symbol))
            {
                return Err(VakintError::NumeratorNotReplaced(
                    EXTERNAL_MOMENTUM_SYMBOL.into(),
                    test.to_string(),
                ));
            }
        }

        self.numerator = new_numerator;

        Ok(())
    }

    /// The scalar tail must restore closed contraction scopes before a linear
    /// matcher substitution can distribute indexed powers. Open tensor inputs
    /// still require the existing tensor prepass; no projection occurs here.
    pub(crate) fn apply_scalar_numerator_replacement_rules(
        &mut self,
        replacement_rules: &ReplacementRules,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        self.numerator = Self::scalar_product_form(settings, self.numerator.as_view())?;
        self.apply_numerator_replacement_rules(replacement_rules, settings)
    }

    pub fn canonicalize(
        &mut self,
        settings: &VakintSettings,
        replacement_rules: &ReplacementRules,
        short_form: bool,
    ) -> Result<(), VakintError> {
        self.integral = replacement_rules.canonical_topology.to_canonical(
            self.integral.as_view(),
            replacement_rules,
            short_form,
        );
        if matches!(replacement_rules.canonical_topology, Topology::Unknown(_),) {
            return Ok(()); // No further canonicalization possible
        }

        self.apply_numerator_replacement_rules(replacement_rules, settings)?;

        Ok(())
    }

    pub fn identify_vectors_in_numerator(
        numerator: AtomView,
    ) -> Result<Vec<(String, i64)>, VakintError> {
        let mut vectors = HashSet::new();
        // make sure the numerator is in the form of vec_(id_,idx_)
        let vector_matcher_pattern = vk_parse!("vec_(id_,idx_)").unwrap().to_pattern();
        let vector_conditions = Condition::from((vk_symbol!("vec_"), symbol_condition()))
            & Condition::from((vk_symbol!("id_"), number_condition()));
        let vector_match_settings = MatchSettings::default();
        let mut vector_matcher = numerator.pattern_match(
            &vector_matcher_pattern,
            Some(&vector_conditions),
            Some(&vector_match_settings),
        );

        while let Some(m) = vector_matcher.next_detailed() {
            if let Match::FunctionName(vec_symbol) = m.match_stack.get(vk_symbol!("vec_")).unwrap()
            {
                if *vec_symbol == vk_symbol!(LOOP_MOMENTUM_SYMBOL)
                    || *vec_symbol == vk_symbol!(EXTERNAL_MOMENTUM_SYMBOL)
                {
                    vectors.insert((
                        vec_symbol.get_stripped_name().into(),
                        get_integer_from_atom(
                            m.match_stack
                                .get(vk_symbol!("id_"))
                                .unwrap()
                                .to_atom()
                                .as_view(),
                        )
                        .unwrap(),
                    ));
                }
            } else {
                unreachable!("Vector name should be a symbol.")
            }
        }
        Ok(vectors.iter().cloned().collect::<Vec<_>>())
    }

    pub fn tensor_reduce(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        settings.with_resource_budget(|| match settings.tensor_reduction_method {
            TensorReductionMethod::AlphaLoop => self.tensor_reduce_alphaloop(vakint, settings),
            TensorReductionMethod::FeynKit => self.tensor_reduce_feynkit(settings, false),
        })
    }

    fn tensor_reduce_alphaloop(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        let mut form_numerator = self.numerator.clone();
        // Make sure to undo the dot product notation.
        // If it was not used, the command below will do nothing.
        form_numerator = Vakint::convert_from_dot_notation(form_numerator.as_view());

        // println!("VH:: A:: Reduction of numerator:\n{}", form_numerator);

        let vectors = Self::identify_vectors_in_numerator(form_numerator.as_view())?;

        let mut vector_mapping: BTreeMap<Atom, Atom> = BTreeMap::new();

        for (vec, id) in vectors.iter() {
            vector_mapping.insert(
                vk_parse!(format!("{}{}", vec, id).as_str()).unwrap(),
                vk_parse!(format!("{}({})", vec, id).as_str()).unwrap(),
            );
            // NO! Indices can move around in the reduction!, not bound to vector names!
            // let vec_pattern = vk_parse!(format!("{}({},idx_)", vec, id).as_str())
            //     .unwrap()
            //     .to_pattern();
            // let matcher = form_numerator.pattern_match(&vec_pattern, None, None);
            // for m in matcher {
            //     let idx = m.get(&vk_symbol!("idx_")).unwrap();
            //     vector_mapping.insert(
            //         vk_parse!(format!("{}{}({})", vec, id, idx.to_canonical_string()).as_str())
            //             .unwrap(),
            //         vk_parse!(format!("{}({},{})", vec, id, idx.to_canonical_string()).as_str())
            //             .unwrap(),
            //     );
            // }
            vector_mapping.insert(
                vk_parse!(format!("{}{}(idx_)", vec, id).as_str()).unwrap(),
                vk_parse!(format!("{}({},idx_)", vec, id).as_str()).unwrap(),
            );
            form_numerator = form_numerator
                .replace(
                    vk_parse!(format!("{}({},idx_)", vec, id).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .with(
                    vk_parse!(
                        format!(
                            "vec{}({}{},idx_)",
                            if *vec == EXTERNAL_MOMENTUM_SYMBOL {
                                "1"
                            } else {
                                ""
                            },
                            vec,
                            id
                        )
                        .as_str()
                    )
                    .unwrap()
                    .to_pattern(),
                );
        }

        let template =
            Template::parse_template(TEMPLATES.get("run_tensor_reduction.txt").unwrap()).unwrap();

        let mut vars: HashMap<String, String> = HashMap::new();
        let (form_header_additions, form_expression, indices) =
            vakint.prepare_expression_for_form(settings, form_numerator, true, &[])?;
        vars.insert("numerator".into(), form_expression);
        vars.insert("additional_symbols".into(), form_header_additions);

        let rendered = template
            .render(&RenderOptions {
                variables: vars,
                ..Default::default()
            })
            .unwrap();

        // println!("Rendered: {}", rendered);
        let form_result = vakint.run_form(
            settings,
            &["tensorreduce.frm".into(), "pvtab10.h".into()],
            ("run_tensor_reduction.frm".into(), rendered),
            vec![],
            settings.clean_tmp_dir,
            settings.temporary_directory.clone(),
        )?;
        // println!("Raw output: {}", form_result);

        // println!(
        //     "\n>>>> VH:: A:: Indices:\n{}",
        //     indices
        //         .clone()
        //         .iter()
        //         .map(|idx| idx.to_canonical_string())
        //         .collect::<Vec<_>>()
        //         .join(",")
        // );

        // println!(
        //     "\n>>>> VH:: A:: Vector mapping:\n{}",
        //     vector_mapping
        //         .iter()
        //         .map(|(k, v)| format!("{} -> {}", k, v.to_canonical_string()))
        //         .collect::<Vec<_>>()
        //         .join("\n")
        // );
        let mut reduced_numerator =
            vakint.process_form_output(settings, form_result, indices, vector_mapping)?;

        // println!(
        //     "\n>>>> VH:: A:: Vectors:\n{}",
        //     vectors
        //         .iter()
        //         .map(|(vec, id)| format!("{}{} -> {}({})", vec, id, vec, id))
        //         .collect::<Vec<_>>()
        //         .join("\n")
        // );

        // Map back surviving external indices
        for (vec, id) in vectors.iter() {
            reduced_numerator = reduced_numerator
                .replace(
                    vk_parse!(format!("{}{}(idx_)", vec, id).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .with(
                    vk_parse!(format!("{}({},idx_)", vec, id).as_str())
                        .unwrap()
                        .to_pattern(),
                )
        }

        if !settings.use_dot_product_notation {
            reduced_numerator = Vakint::convert_from_dot_notation(reduced_numerator.as_view());
        }

        self.numerator = reduced_numerator;
        Ok(())
    }

    /// Keep already reduced scalar contractions compact across the private
    /// five-loop RustRed handoff. The standalone tensor-reduction API continues
    /// to honor the requested indexed/dot presentation exactly as before.
    fn tensor_reduce_for_evaluation(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        self.stage_numerical_inputs(vakint, settings)?;
        let compact = (settings.rustred_artifact_path.is_some() || vakint.has_rustred_loop(5))
            && matches!(
                settings.tensor_reduction_method,
                TensorReductionMethod::FeynKit
            )
            && vakint
                .topologies
                .match_topologies_to_user_input(
                    self.integral.as_view(),
                    settings.allow_unknown_integrals,
                )?
                .is_some_and(|matched| {
                    matches!(matched.canonical_topology, Topology::FiveLoop(_))
                        && matches!(
                            settings.evaluation_order.0.iter().find(|method| {
                                method.supports(settings, &matched.canonical_topology)
                            }),
                            Some(EvaluationMethod::RustRed(_))
                        )
                });
        if compact {
            let mut scalar_settings = settings.clone();
            scalar_settings.use_dot_product_notation = true;
            self.tensor_reduce(vakint, &scalar_settings)
        } else {
            self.tensor_reduce(vakint, settings)
        }
    }

    /// Restore closed scalar products through the native index-scope boundary.
    /// This does not authorize projection of an unprocessed loop tensor.
    pub(crate) fn scalar_product_form(
        settings: &VakintSettings,
        numerator: AtomView,
    ) -> Result<Atom, VakintError> {
        // Scalar coefficients and explicit products do not require the tensor
        // registry. Only known Lorentz heads enter the checked tensor boundary.
        fn has_indexed_loop(value: AtomView<'_>) -> bool {
            if matches!(value, AtomView::Fun(function)
                if function.get_symbol() == S.k && function.get_nargs() != 1)
            {
                return true;
            }
            value.children().any(has_indexed_loop)
        }
        fn has_native_lorentz(value: AtomView<'_>) -> bool {
            if let AtomView::Fun(function) = value {
                if function.get_symbol() == SPENSO_TAG.pure_scalar {
                    return false;
                }
                if [S.k, S.p, S.g, ETS.metric].contains(&function.get_symbol()) {
                    let arguments = function.iter().collect::<Vec<_>>();
                    if arguments.iter().any(|argument| {
                        VakintTerm::spenso_minkowski_slot_parts(*argument).is_some()
                    }) || (arguments.len() >= 2
                        && VakintTerm::is_compact_spenso_minkowski_representation(
                            arguments[arguments.len() - 2],
                        )
                        && VakintTerm::is_dot_dummy_index(arguments[arguments.len() - 1]))
                    {
                        return true;
                    }
                }
            }
            value.children().any(has_native_lorentz)
        }
        let epsilon = vk_parse!(&settings.epsilon_symbol).map_err(VakintError::SymbolicaError)?;
        let dimension = Atom::num(4) - Atom::num(2) * epsilon;
        let mut incompatible_metric = None;
        numerator.visitor(&mut |value| {
            if let AtomView::Fun(function) = value {
                if function.get_symbol() == SPENSO_TAG.pure_scalar {
                    return false;
                }
                if [S.g, ETS.metric].contains(&function.get_symbol()) {
                    for argument in function.iter() {
                        if let Some((source_dimension, _)) =
                            Self::spenso_minkowski_slot_parts(argument)
                            && !(source_dimension.to_owned() - &dimension)
                                .expand()
                                .is_zero()
                        {
                            incompatible_metric = Some(source_dimension.to_owned());
                        }
                    }
                }
            }
            true
        });
        if let Some(source_dimension) = incompatible_metric {
            return Err(VakintError::InvalidNumerator(format!(
                "strict RustRed scalar admission cannot retarget native metric dimension {source_dimension} to {dimension}; run the tensor prepass explicitly"
            )));
        }
        // Close known Lorentz components before expansion; free interfaces and
        // powered dummy scopes retain the existing LorentzTensor semantics.
        let numerator = lorentz::LorentzTensor::parse(numerator, &dimension)?.expression;
        let indexed_loop = has_indexed_loop(numerator.as_view());
        if !indexed_loop && !has_native_lorentz(numerator.as_view()) {
            return Ok(numerator);
        }
        let mut scalar_settings = settings.clone();
        scalar_settings.use_dot_product_notation = true;
        let mut term = Self {
            numerator: numerator.to_owned(),
            integral: Atom::one(),
            vectors: Vec::new(),
        };
        // Native external slots belong to the existing tensor owner, including
        // their free interfaces. Integrated tensors still require the prepass.
        term.tensor_reduce_feynkit(&scalar_settings, true)?;
        Ok(term.numerator)
    }

    fn tensor_reduce_feynkit(
        &mut self,
        settings: &VakintSettings,
        scalar_only: bool,
    ) -> Result<(), VakintError> {
        initialize_spenso_reps();
        let epsilon = vk_parse!(&settings.epsilon_symbol).map_err(VakintError::SymbolicaError)?;
        let dimension = Atom::num(4) - Atom::num(2) * epsilon;
        // The native reducer operates on explicit indexed vectors.  In
        // particular, a compact contraction such as dot(p(1),k(1)) contains
        // one power of the integrated momentum and must not be treated as an
        // opaque scalar factor.  Expand Vakint's dot shorthand before adding
        // Spenso slots, just as the AlphaLoop backend does before invoking
        // FORM. Closed loop-loop products, in contrast, use the native compact
        // invariant path and must not be reopened into a spurious tensor rank.
        Self::validate_spenso_loop_dots(self.numerator.as_view())?;
        // Close independent scalar contractions before the collector distributes
        // products of tensor sums. The checked reducer below still owns loop
        // projection, scalar-only admission, and the tensor rank limits.
        let numerator =
            lorentz::LorentzTensor::parse(self.numerator.as_view(), &dimension)?.expression;
        let vakint_dots = numerator.replace_map(|term, _context, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            let arguments = function.iter().collect::<Vec<_>>();
            if function.get_symbol() == SPENSO_TAG.dot
                && arguments.len() == 2
                && Self::is_reducible_spenso_loop_dot(&arguments)
            {
                **out = FunctionBuilder::new(S.dot).add_args(arguments).finish();
            }
        });
        // Closed loop invariants are coefficients of the remaining tensor.
        // Keep them opaque through collection: distributing their scalar powers
        // would introduce a projector expansion where no tensor rank remains.
        let mut coefficient_aliases = BTreeMap::<Atom, Atom>::new();
        let coefficient_head = vk_symbol!("tensor_scalar_coefficient");
        let mut occupied_aliases = HashSet::new();
        let _ = vakint_dots.replace_map(|term, _, _| {
            if matches!(term, AtomView::Fun(function) if function.get_symbol() == coefficient_head)
            {
                occupied_aliases.insert(term.to_owned());
            }
        });
        let mut coefficient_index = 0;
        let tensor_numerator = vakint_dots.replace_map(|term, _context, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            if function.get_symbol() != S.dot || function.get_nargs() != 2 {
                return;
            }
            let arguments = function.iter().collect::<Vec<_>>();
            if !arguments.iter().all(|argument| {
                matches!(argument, AtomView::Fun(vector)
                    if vector.get_symbol() == S.k && vector.get_nargs() == 1
                        && get_integer_from_atom(vector.iter().next().unwrap()).is_some())
            }) {
                return;
            }
            // One invariant gets one alias even when it occurs in thousands
            // of scalar monomials. Existing user functions were inventoried
            // above, so alias allocation never rescans the whole numerator.
            let alias = coefficient_aliases
                .entry(term.to_owned())
                .or_insert_with(|| {
                    loop {
                        let alias = function!(coefficient_head, coefficient_index);
                        coefficient_index += 1;
                        if !occupied_aliases.contains(&alias) {
                            break alias;
                        }
                    }
                });
            **out = alias.clone();
        });
        let indexed_numerator = Vakint::convert_from_dot_notation(tensor_numerator.as_view());
        let reducer = TensorReducer::new(dimension).with_integrated_head(S.k);
        let mut next_bridge_dummy = Self::next_dot_dummy_index(indexed_numerator.as_view());
        let lorentz = TensorCollectFilter::Reps([Minkowski {}.into()]);
        // The frontend still carries Vakint's raw-index k/p/g leaves here.
        // Keep their original representation metadata together until each
        // selected sector enters the checked tensor reducer. Scalar coefficient
        // subtrees stay opaque in the existing filtered collector.
        let selected = indexed_numerator.collect_with_map(|value| match value {
            AtomView::Fun(function) => {
                [S.k, S.p, S.g].contains(&function.get_symbol())
                    || (function.get_symbol() == SPENSO_TAG.dot
                        && Self::is_reducible_spenso_loop_dot(&function.iter().collect::<Vec<_>>()))
                    || lorentz.matches(value)
            }
            AtomView::Pow(power) => {
                // A closed indexed power owns a fresh dummy scope per copy.
                // The typed reducer handles that scope; raw collection must
                // not distribute the base and merge those labels first.
                let base = power.get_base();
                [S.k, S.p, S.g]
                    .into_iter()
                    .any(|head| Self::contains_function_head(base, head))
                    || base.contains_symbol(Minkowski::selfless_symbol())
            }
            _ => false,
        });
        let mut failure = None;
        let mut reduced = selected
            .map_root(|root| {
                root.map_collects(|collected, _, out| {
                    if failure.is_some() {
                        return;
                    }
                    let AtomView::Fun(collected) = collected else {
                        unreachable!("map_collects visits only collector wrappers")
                    };
                    let term = collected.iter().next().expect("collector wrapper payload");
                    match Self::tensor_reduce_feynkit_term(
                        term,
                        &reducer,
                        &mut next_bridge_dummy,
                        scalar_only,
                    ) {
                        Ok(value) => **out = value,
                        Err(error) => failure = Some(error),
                    }
                })
                .unwrap_collect()
            })
            .into_inner();
        if let Some(error) = failure {
            return Err(error);
        }
        let mut restored = AliasedAtom::from(reduced);
        for (coefficient, alias) in coefficient_aliases {
            restored.register_alias(alias, coefficient);
        }
        reduced = restored.into_inner();
        if !settings.use_dot_product_notation {
            reduced = Vakint::convert_from_dot_notation(reduced.as_view());
            reduced = Self::repair_expanded_spenso_vector_slots(reduced.as_view());
        }
        self.numerator = reduced;
        Ok(())
    }

    fn tensor_reduce_feynkit_term(
        numerator: AtomView<'_>,
        reducer: &TensorReducer,
        next_bridge_dummy: &mut usize,
        scalar_only: bool,
    ) -> Result<Atom, VakintError> {
        let dimension = reducer.dimension();
        let opaque_minkowski_slots = Self::collect_opaque_minkowski_slots(numerator);
        let opaque_compact_vectors = Self::collect_opaque_compact_vectors(numerator, dimension);
        let numerator = Self::vakint_to_spenso_numerator(
            numerator,
            dimension,
            &opaque_minkowski_slots,
            scalar_only,
        );
        // Native admission expects atomic representation payloads. Encode the
        // original compound dimension and generated dummy indices reversibly,
        // then restore them before the existing Vakint notation bridge.
        let cooking = idenso::CookSettings::indices()
            .with_mode(idenso::CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let numerator = cooking.try_cook(numerator.as_view()).map_err(|error| {
            feynkit_tensor::TensorReductionError::TensorAlgebra(error.to_string())
        })?;
        let reduced = reducer.reduce(&reducer.admit(numerator)?)?;
        if scalar_only && reduced.max_rank() != 0 {
            return Err(VakintError::InvalidNumerator(format!(
                "RustRed scalar input still requires tensor projection of rank {}; run the tensor prepass first",
                reduced.max_rank(),
            )));
        }
        let reduced = reduced.into_expression()?.into_expression();
        let reduced = cooking.uncook(reduced.as_view());
        let reduced = Self::bridge_opaque_projectors(reduced.as_view(), next_bridge_dummy);
        let reduced =
            Self::spenso_to_vakint_numerator(reduced.as_view(), dimension, &opaque_compact_vectors);
        if Self::contains_feynkit_bridge_syntax(reduced.as_view()) {
            return Err(VakintError::FeynKitBridgeSyntax(
                reduced.to_canonical_string(),
            ));
        }
        Ok(reduced)
    }

    fn vakint_to_spenso_numerator(
        numerator: AtomView<'_>,
        dimension: &Atom,
        opaque_slots: &BTreeMap<Atom, Atom>,
        scalar_only: bool,
    ) -> Atom {
        numerator.replace_map(|term, _context, out| {
            if let AtomView::Add(sum) = term {
                // A preserved power can contain additive branches with the
                // same dummy label but different original dimensions. Attach
                // provenance within each branch before the typed reducer
                // handles that power; never distribute the surrounding graph.
                **out = Atom::add_many(
                    sum.iter()
                        .map(|summand| {
                            let slots = Self::collect_opaque_minkowski_slots(summand);
                            Self::vakint_to_spenso_numerator(
                                summand,
                                dimension,
                                &slots,
                                scalar_only,
                            )
                        })
                        .collect::<Vec<_>>(),
                );
                return;
            }
            let AtomView::Fun(function) = term else {
                return;
            };
            let arguments = function.iter().collect::<Vec<_>>();
            // Scalar admission preserves native external source interfaces.
            // The public tensor prepass keeps its existing target-dimension
            // convention; raw Vakint p(id,index) also keeps its old notation.
            let represented_external = scalar_only
                && function.get_symbol() == S.p
                && (arguments
                    .last()
                    .is_some_and(|slot| Self::spenso_minkowski_slot_parts(*slot).is_some())
                    || (arguments.len() >= 3
                        && Self::is_compact_spenso_minkowski_representation(
                            arguments[arguments.len() - 2],
                        )
                        && Self::is_dot_dummy_index(arguments[arguments.len() - 1])));
            if represented_external
                && let Some(vector) = Self::bridge_tagged_vector(
                    function.get_symbol(),
                    &arguments,
                    dimension,
                    opaque_slots,
                )
            {
                **out = vector;
                return;
            }
            if (function.get_symbol() == S.k || function.get_symbol() == S.p)
                && arguments.len() >= 3
                && Self::is_compact_spenso_minkowski_representation(arguments[arguments.len() - 2])
                && Self::is_dot_dummy_index(arguments[arguments.len() - 1])
            {
                let slot = Self::retarget_or_create_spenso_slot(
                    arguments[arguments.len() - 1],
                    dimension,
                    opaque_slots,
                );
                **out = arguments[..arguments.len() - 2]
                    .iter()
                    .fold(
                        FunctionBuilder::new(function.get_symbol()),
                        |builder, argument| builder.add_arg(*argument),
                    )
                    .add_arg(slot)
                    .finish();
            } else if (function.get_symbol() == S.k || function.get_symbol() == S.p)
                && arguments.len() == 2
            {
                if Self::is_compact_spenso_minkowski_representation(arguments[1]) {
                    **out = term.to_owned();
                    return;
                }
                let slot =
                    Self::retarget_or_create_spenso_slot(arguments[1], dimension, opaque_slots);
                **out = FunctionBuilder::new(function.get_symbol())
                    .add_arg(arguments[0])
                    .add_arg(slot)
                    .finish();
            } else if function.get_symbol() == S.g && arguments.len() == 2 {
                let left =
                    Self::retarget_or_create_spenso_slot(arguments[0], dimension, opaque_slots);
                let right =
                    Self::retarget_or_create_spenso_slot(arguments[1], dimension, opaque_slots);
                **out = ETS.metric(left, right);
            } else if function.get_symbol() == Minkowski::selfless_symbol()
                && matches!(arguments.len(), 1 | 2)
            {
                let mut representation =
                    FunctionBuilder::new(Minkowski::selfless_symbol()).add_arg(dimension);
                if let Some(index) = arguments.get(1) {
                    representation =
                        representation.add_arg(Self::protect_opaque_index(*index, opaque_slots));
                }
                **out = representation.finish();
            } else if function.get_symbol().has_tag(&SPENSO_TAG.rank1) {
                if let Some(vector) = Self::bridge_tagged_vector(
                    function.get_symbol(),
                    &arguments,
                    dimension,
                    opaque_slots,
                ) {
                    **out = vector;
                } else if arguments.last().is_some_and(|argument| {
                    Self::is_compact_spenso_minkowski_representation(*argument)
                }) {
                    // A compact tagged vector is already a scalar object from
                    // the reducer's point of view. Preserve its original
                    // representation instead of recursively retargeting it.
                    **out = term.to_owned();
                }
            } else if let Some(vector) = Self::bridge_untagged_indexed_vector(
                function.get_symbol(),
                &arguments,
                dimension,
                opaque_slots,
            ) {
                **out = vector;
            } else if arguments.len() >= 2
                && Self::is_compact_spenso_minkowski_representation(arguments[arguments.len() - 2])
                && Self::is_dot_dummy_index(arguments[arguments.len() - 1])
            {
                let original = arguments[..arguments.len() - 1]
                    .iter()
                    .fold(
                        FunctionBuilder::new(function.get_symbol()),
                        |builder, argument| builder.add_arg(*argument),
                    )
                    .finish();
                let slot = Self::retarget_or_create_spenso_slot(
                    arguments[arguments.len() - 1],
                    dimension,
                    opaque_slots,
                );
                **out = FunctionBuilder::new(*FEYNKIT_EXTERNAL_VECTOR)
                    .add_arg(
                        FunctionBuilder::new(SPENSO_TAG.pure_scalar)
                            .add_arg(original)
                            .finish(),
                    )
                    .add_arg(slot)
                    .finish();
            }
        })
    }

    fn bridge_untagged_indexed_vector(
        head: Symbol,
        arguments: &[AtomView<'_>],
        dimension: &Atom,
        opaque_slots: &BTreeMap<Atom, Atom>,
    ) -> Option<Atom> {
        let slot = *arguments.last()?;
        let (original_dimension, _index) = Self::spenso_minkowski_slot_parts(slot)?;
        if arguments[..arguments.len() - 1]
            .iter()
            .any(|argument| {
                matches!(argument, AtomView::Fun(function) if function.get_symbol().has_tag(&SPENSO_TAG.representation))
            })
        {
            return None;
        }
        let original = arguments[..arguments.len() - 1]
            .iter()
            .fold(FunctionBuilder::new(head), |builder, argument| {
                builder.add_arg(*argument)
            })
            .add_arg(
                FunctionBuilder::new(Minkowski::selfless_symbol())
                    .add_arg(original_dimension)
                    .finish(),
            )
            .finish();
        Some(
            FunctionBuilder::new(*FEYNKIT_EXTERNAL_VECTOR)
                .add_arg(
                    FunctionBuilder::new(SPENSO_TAG.pure_scalar)
                        .add_arg(original)
                        .finish(),
                )
                .add_arg(Self::retarget_or_create_spenso_slot(
                    slot,
                    dimension,
                    opaque_slots,
                ))
                .finish(),
        )
    }

    fn bridge_tagged_vector(
        head: Symbol,
        arguments: &[AtomView<'_>],
        dimension: &Atom,
        opaque_slots: &BTreeMap<Atom, Atom>,
    ) -> Option<Atom> {
        let last = *arguments.last()?;
        let (original, slot) =
            if let Some((original_dimension, _index)) = Self::spenso_minkowski_slot_parts(last) {
                let compact_representation = FunctionBuilder::new(Minkowski::selfless_symbol())
                    .add_arg(original_dimension)
                    .finish();
                let original = arguments[..arguments.len() - 1]
                    .iter()
                    .fold(FunctionBuilder::new(head), |builder, argument| {
                        builder.add_arg(*argument)
                    })
                    .add_arg(compact_representation)
                    .finish();
                let slot = Self::retarget_or_create_spenso_slot(last, dimension, opaque_slots);
                (original, slot)
            } else if arguments.len() >= 2
                && Self::is_compact_spenso_minkowski_representation(arguments[arguments.len() - 2])
                && Self::is_dot_dummy_index(last)
            {
                let original = arguments[..arguments.len() - 1]
                    .iter()
                    .fold(FunctionBuilder::new(head), |builder, argument| {
                        builder.add_arg(*argument)
                    })
                    .finish();
                let slot = Self::retarget_or_create_spenso_slot(last, dimension, opaque_slots);
                (original, slot)
            } else if Self::is_dot_dummy_index(last) {
                let original = arguments[..arguments.len() - 1]
                    .iter()
                    .fold(FunctionBuilder::new(head), |builder, argument| {
                        builder.add_arg(*argument)
                    })
                    .finish();
                let slot = Self::retarget_or_create_spenso_slot(last, dimension, opaque_slots);
                (original, slot)
            } else {
                return None;
            };

        Some(
            FunctionBuilder::new(*FEYNKIT_EXTERNAL_VECTOR)
                .add_arg(
                    FunctionBuilder::new(SPENSO_TAG.pure_scalar)
                        .add_arg(original)
                        .finish(),
                )
                .add_arg(slot)
                .finish(),
        )
    }

    fn is_spenso_minkowski_representation(argument: AtomView<'_>) -> bool {
        matches!(
            argument,
            AtomView::Fun(representation)
                if representation.get_symbol().get_namespace() == "spenso"
                    && representation.get_symbol().get_stripped_name() == "mink"
                    && matches!(representation.get_nargs(), 1 | 2)
        )
    }

    fn is_compact_spenso_minkowski_representation(argument: AtomView<'_>) -> bool {
        matches!(
            argument,
            AtomView::Fun(representation)
                if Self::is_spenso_minkowski_representation(argument)
                    && representation.get_nargs() == 1
        )
    }

    fn compact_spenso_minkowski_vector_representation(
        argument: AtomView<'_>,
    ) -> Option<AtomView<'_>> {
        let AtomView::Fun(vector) = argument else {
            return None;
        };
        let arguments = vector.iter().collect::<Vec<_>>();
        let representation = *arguments.last()?;
        if !Self::is_compact_spenso_minkowski_representation(representation)
            || arguments[..arguments.len() - 1]
                .iter()
                .any(|argument| Self::contains_spenso_representation(*argument))
        {
            return None;
        }
        Some(representation)
    }

    fn is_reducible_spenso_loop_dot(arguments: &[AtomView<'_>]) -> bool {
        let [left, right] = arguments else {
            return false;
        };
        let (Some(left_representation), Some(right_representation)) = (
            Self::compact_spenso_minkowski_vector_representation(*left),
            Self::compact_spenso_minkowski_vector_representation(*right),
        ) else {
            return false;
        };
        left_representation == right_representation
            && arguments.iter().any(
                |argument| matches!(argument, AtomView::Fun(vector) if vector.get_symbol() == S.k),
            )
    }

    fn validate_spenso_loop_dots(expression: AtomView<'_>) -> Result<(), VakintError> {
        if let AtomView::Fun(function) = expression
            && function.get_symbol() == SPENSO_TAG.dot
        {
            let arguments = function.iter().collect::<Vec<_>>();
            if arguments
                .iter()
                .any(|argument| Self::contains_function_head(*argument, S.k))
                && !Self::is_reducible_spenso_loop_dot(&arguments)
            {
                return Err(VakintError::FeynKitLoopDot(
                    expression.to_owned().to_canonical_string(),
                ));
            }
        }
        for child in expression.children() {
            Self::validate_spenso_loop_dots(child)?;
        }
        Ok(())
    }

    fn contains_spenso_representation(argument: AtomView<'_>) -> bool {
        if let AtomView::Fun(function) = argument
            && function.get_symbol().has_tag(&SPENSO_TAG.representation)
        {
            return true;
        }
        argument
            .children()
            .any(Self::contains_spenso_representation)
    }

    fn contains_function_head(argument: AtomView<'_>, head: Symbol) -> bool {
        if matches!(argument, AtomView::Fun(function) if function.get_symbol() == head) {
            return true;
        }
        argument
            .children()
            .any(|child| Self::contains_function_head(child, head))
    }

    fn is_dot_dummy_index(argument: AtomView<'_>) -> bool {
        matches!(
            argument,
            AtomView::Fun(function)
                if function.get_symbol() == S.dot_dummy_ind && function.get_nargs() == 1
        )
    }

    fn repair_expanded_spenso_vector_slots(numerator: AtomView<'_>) -> Atom {
        numerator.replace_map(|term, _context, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            let arguments = function.iter().collect::<Vec<_>>();
            if arguments.len() < 2
                || !Self::is_compact_spenso_minkowski_representation(arguments[arguments.len() - 2])
                || !Self::is_dot_dummy_index(arguments[arguments.len() - 1])
            {
                return;
            }
            let AtomView::Fun(compact_representation) = arguments[arguments.len() - 2] else {
                return;
            };
            let Some(representation_dimension) = compact_representation.iter().next() else {
                return;
            };
            let indexed_representation = FunctionBuilder::new(Minkowski::selfless_symbol())
                .add_arg(representation_dimension)
                .add_arg(arguments[arguments.len() - 1])
                .finish();
            **out = arguments[..arguments.len() - 2]
                .iter()
                .fold(
                    FunctionBuilder::new(function.get_symbol()),
                    |builder, argument| builder.add_arg(*argument),
                )
                .add_arg(indexed_representation)
                .finish();
        })
    }

    fn retarget_or_create_spenso_slot(
        argument: AtomView<'_>,
        dimension: &Atom,
        opaque_slots: &BTreeMap<Atom, Atom>,
    ) -> Atom {
        if let AtomView::Fun(representation) = argument
            && Self::is_spenso_minkowski_representation(argument)
            && representation.get_nargs() == 2
        {
            let index = representation.iter().nth(1).unwrap();
            return FunctionBuilder::new(Minkowski::selfless_symbol())
                .add_arg(dimension)
                .add_arg(Self::protect_opaque_index(index, opaque_slots))
                .finish();
        }
        FunctionBuilder::new(Minkowski::selfless_symbol())
            .add_arg(dimension)
            .add_arg(Self::protect_opaque_index(argument, opaque_slots))
            .finish()
    }

    fn protect_opaque_index(index: AtomView<'_>, opaque_slots: &BTreeMap<Atom, Atom>) -> Atom {
        opaque_slots.get(&index.to_owned()).map_or_else(
            || index.to_owned(),
            |original_slot| {
                FunctionBuilder::new(*FEYNKIT_OPAQUE_INDEX)
                    .add_arg(original_slot)
                    .finish()
            },
        )
    }

    fn original_opaque_slot(index: AtomView<'_>) -> Option<AtomView<'_>> {
        let AtomView::Fun(function) = index else {
            return None;
        };
        if function.get_symbol() != *FEYNKIT_OPAQUE_INDEX || function.get_nargs() != 1 {
            return None;
        }
        function.iter().next()
    }

    fn contains_feynkit_bridge_syntax(expression: AtomView<'_>) -> bool {
        if let AtomView::Fun(function) = expression
            && (function.get_symbol() == *FEYNKIT_EXTERNAL_VECTOR
                || function.get_symbol() == *FEYNKIT_OPAQUE_INDEX)
        {
            return true;
        }
        expression
            .children()
            .any(Self::contains_feynkit_bridge_syntax)
    }

    fn spenso_minkowski_slot_parts(slot: AtomView<'_>) -> Option<(AtomView<'_>, AtomView<'_>)> {
        let AtomView::Fun(representation) = slot else {
            return None;
        };
        if !Self::is_spenso_minkowski_representation(slot) || representation.get_nargs() != 2 {
            return None;
        }
        let arguments = representation.iter().collect::<Vec<_>>();
        Some((arguments[0], arguments[1]))
    }

    fn collect_opaque_minkowski_slots(numerator: AtomView<'_>) -> BTreeMap<Atom, Atom> {
        fn collect(expression: AtomView<'_>, slots: &mut BTreeMap<Atom, Atom>) {
            if let AtomView::Fun(function) = expression {
                if function.get_symbol() == Minkowski::selfless_symbol()
                    || function.get_symbol() == ETS.metric
                    || function.get_symbol().has_tag(&SPENSO_TAG.rank1)
                {
                    return;
                }
                for argument in function.iter() {
                    if let Some((_dimension, index)) =
                        VakintTerm::spenso_minkowski_slot_parts(argument)
                    {
                        slots
                            .entry(index.to_owned())
                            .or_insert_with(|| argument.to_owned());
                    } else {
                        collect(argument, slots);
                    }
                }
                return;
            }
            for child in expression.children() {
                collect(child, slots);
            }
        }

        let mut slots = BTreeMap::new();
        collect(numerator, &mut slots);
        slots
    }

    fn collect_opaque_compact_vectors(
        numerator: AtomView<'_>,
        dimension: &Atom,
    ) -> BTreeMap<Atom, Atom> {
        fn collect(expression: AtomView<'_>, dimension: &Atom, vectors: &mut BTreeMap<Atom, Atom>) {
            if let AtomView::Fun(function) = expression {
                if function.get_symbol() == Minkowski::selfless_symbol()
                    || function.get_symbol() == ETS.metric
                    || function.get_symbol() == S.k
                    || function.get_symbol() == S.p
                    || function.get_symbol() == S.g
                    || function.get_symbol().has_tag(&SPENSO_TAG.rank1)
                {
                    return;
                }
                let arguments = function.iter().collect::<Vec<_>>();
                let slots = arguments
                    .iter()
                    .enumerate()
                    .filter_map(|(position, argument)| {
                        VakintTerm::spenso_minkowski_slot_parts(*argument)
                            .map(|(slot_dimension, _index)| (position, slot_dimension))
                    })
                    .collect::<Vec<_>>();
                if let [(position, original_dimension)] = slots.as_slice() {
                    let normalized = arguments.iter().enumerate().fold(
                        FunctionBuilder::new(function.get_symbol()),
                        |builder, (argument_position, argument)| {
                            if argument_position == *position {
                                builder.add_arg(
                                    FunctionBuilder::new(Minkowski::selfless_symbol())
                                        .add_arg(dimension)
                                        .finish(),
                                )
                            } else {
                                builder.add_arg(*argument)
                            }
                        },
                    );
                    let original = arguments.iter().enumerate().fold(
                        FunctionBuilder::new(function.get_symbol()),
                        |builder, (argument_position, argument)| {
                            if argument_position == *position {
                                builder.add_arg(
                                    FunctionBuilder::new(Minkowski::selfless_symbol())
                                        .add_arg(*original_dimension)
                                        .finish(),
                                )
                            } else {
                                builder.add_arg(*argument)
                            }
                        },
                    );
                    vectors.insert(normalized.finish(), original.finish());
                }
                for argument in arguments {
                    if !VakintTerm::is_spenso_minkowski_representation(argument) {
                        collect(argument, dimension, vectors);
                    }
                }
                return;
            }
            for child in expression.children() {
                collect(child, dimension, vectors);
            }
        }

        let mut vectors = BTreeMap::new();
        collect(numerator, dimension, &mut vectors);
        vectors
    }

    fn next_dot_dummy_index(numerator: AtomView<'_>) -> usize {
        let dummy_pattern = function!(S.dot_dummy_ind, S.a_).to_pattern();
        numerator
            .replace(&dummy_pattern)
            .match_iter()
            .filter_map(|matched| usize::try_from(&matched[&S.a_]).ok())
            .max()
            .unwrap_or(0)
            + 1
    }

    fn bridge_opaque_projectors(numerator: AtomView<'_>, next_dummy: &mut usize) -> Atom {
        let projected = numerator.replace_map(|term, _context, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            let arguments = function.iter().collect::<Vec<_>>();

            if function.get_symbol() == ETS.metric && arguments.len() == 2 {
                let Some((_left_dimension, left_index)) =
                    Self::spenso_minkowski_slot_parts(arguments[0])
                else {
                    return;
                };
                let Some((_right_dimension, right_index)) =
                    Self::spenso_minkowski_slot_parts(arguments[1])
                else {
                    return;
                };
                let left_opaque = Self::original_opaque_slot(left_index);
                let right_opaque = Self::original_opaque_slot(right_index);
                if left_opaque.is_some() || right_opaque.is_some() {
                    **out = FunctionBuilder::new(S.g)
                        .add_arg(left_opaque.unwrap_or(arguments[0]))
                        .add_arg(right_opaque.unwrap_or(arguments[1]))
                        .finish();
                }
                return;
            }

            if function.get_symbol() != S.p && function.get_symbol() != *FEYNKIT_EXTERNAL_VECTOR {
                return;
            }
            let Some(slot) = arguments.last() else {
                return;
            };
            let Some((slot_dimension, index)) = Self::spenso_minkowski_slot_parts(*slot) else {
                return;
            };
            let Some(opaque_slot) = Self::original_opaque_slot(index) else {
                return;
            };
            let dummy = S.dot_dummy_ind(*next_dummy);
            *next_dummy += 1;
            if function.get_symbol() == *FEYNKIT_EXTERNAL_VECTOR {
                let Some(AtomView::Fun(metadata)) = arguments.first() else {
                    return;
                };
                if metadata.get_symbol() != SPENSO_TAG.pure_scalar || metadata.get_nargs() != 1 {
                    return;
                }
                let original = metadata.iter().next().unwrap();
                let original_dimension =
                    Self::compact_spenso_minkowski_vector_representation(original)
                        .and_then(|representation| match representation {
                            AtomView::Fun(representation) => representation.iter().next(),
                            _ => None,
                        })
                        .unwrap_or(slot_dimension);
                let indexed_slot = FunctionBuilder::new(Minkowski::selfless_symbol())
                    .add_arg(original_dimension)
                    .add_arg(&dummy)
                    .finish();
                let Some(vector) =
                    Self::replace_compact_minkowski_slot(original, indexed_slot.as_view())
                else {
                    return;
                };
                **out = FunctionBuilder::new(ETS.metric)
                    .add_arg(opaque_slot)
                    .add_arg(&indexed_slot)
                    .finish()
                    * vector;
                return;
            }
            // Keep Vakint's established metric/vector boundary convention.
            // GammaLoop upgrades this paired raw dummy to its topology-local
            // Spenso index before simplifying the metric. Replacing the pair
            // by a direct shared slot changes bispinor-projector semantics in
            // the existing renormalization pipeline.
            let vector = arguments[..arguments.len() - 1]
                .iter()
                .fold(
                    FunctionBuilder::new(function.get_symbol()),
                    |builder, argument| builder.add_arg(*argument),
                )
                .add_arg(&dummy)
                .finish();
            **out = FunctionBuilder::new(S.g)
                .add_arg(opaque_slot)
                .add_arg(&dummy)
                .finish()
                * vector;
        });

        projected.replace_map(|term, _context, out| {
            let Some((_slot_dimension, index)) = Self::spenso_minkowski_slot_parts(term) else {
                return;
            };
            if let Some(original_slot) = Self::original_opaque_slot(index) {
                **out = original_slot.to_owned();
            }
        })
    }

    fn replace_compact_minkowski_slot(
        vector: AtomView<'_>,
        replacement: AtomView<'_>,
    ) -> Option<Atom> {
        let AtomView::Fun(function) = vector else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        let compact_slot_position = arguments
            .last()
            .is_some_and(|argument| Self::is_compact_spenso_minkowski_representation(*argument))
            .then_some(arguments.len() - 1);
        let mut builder = FunctionBuilder::new(function.get_symbol());
        for (position, argument) in arguments.iter().enumerate() {
            if Some(position) != compact_slot_position {
                builder = builder.add_arg(*argument);
            }
        }
        Some(builder.add_arg(replacement).finish())
    }

    fn spenso_to_vakint_numerator(
        numerator: AtomView<'_>,
        dimension: &Atom,
        opaque_compact_vectors: &BTreeMap<Atom, Atom>,
    ) -> Atom {
        let dots_mapped = numerator.replace_map(|term, _context, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            if function.get_symbol() != SPENSO_TAG.dot || function.get_nargs() != 2 {
                return;
            }
            let arguments = function.iter().collect::<Vec<_>>();
            let Some(left) =
                Self::compact_spenso_vector(arguments[0], dimension, opaque_compact_vectors)
            else {
                return;
            };
            let Some(right) =
                Self::compact_spenso_vector(arguments[1], dimension, opaque_compact_vectors)
            else {
                return;
            };
            **out = FunctionBuilder::new(S.dot)
                .add_arg(left)
                .add_arg(right)
                .finish();
        });

        dots_mapped.replace_map(|term, _context, out| {
            let AtomView::Fun(function) = term else {
                return;
            };
            if let Some(vector) = opaque_compact_vectors.get(&term.to_owned()) {
                **out = vector.clone();
                return;
            }
            let arguments = function.iter().collect::<Vec<_>>();
            if function.get_symbol() == *FEYNKIT_EXTERNAL_VECTOR {
                if let Some(vector) = Self::restore_external_vector(&arguments, dimension) {
                    **out = vector;
                }
                return;
            }
            if function.get_symbol() == ETS.metric && arguments.len() == 2 {
                if arguments.iter().any(|slot| {
                    Self::spenso_minkowski_slot_parts(*slot)
                        .is_some_and(|(_dimension, index)| Self::is_dot_dummy_index(index))
                }) {
                    return;
                }
                let (Some(left), Some(right)) = (
                    Self::spenso_slot_index(arguments[0], dimension),
                    Self::spenso_slot_index(arguments[1], dimension),
                ) else {
                    return;
                };
                **out = FunctionBuilder::new(S.g)
                    .add_arg(left)
                    .add_arg(right)
                    .finish();
            } else if function.get_symbol() == S.k || function.get_symbol() == S.p {
                let Some(slot) = arguments.last() else {
                    return;
                };
                let Some(index) = Self::spenso_slot_index(*slot, dimension) else {
                    return;
                };
                **out = arguments[..arguments.len() - 1]
                    .iter()
                    .fold(
                        FunctionBuilder::new(function.get_symbol()),
                        |builder, argument| builder.add_arg(*argument),
                    )
                    .add_arg(index)
                    .finish();
            }
        })
    }

    fn compact_spenso_vector(
        vector: AtomView<'_>,
        dimension: &Atom,
        opaque_compact_vectors: &BTreeMap<Atom, Atom>,
    ) -> Option<Atom> {
        if let Some(original) = opaque_compact_vectors.get(&vector.to_owned()) {
            return Some(original.clone());
        }
        let AtomView::Fun(function) = vector else {
            return None;
        };
        let arguments = function.iter().collect::<Vec<_>>();
        if function.get_symbol() == *FEYNKIT_EXTERNAL_VECTOR {
            return Self::restore_external_vector(&arguments, dimension);
        }
        if function.get_symbol() != S.k && function.get_symbol() != S.p {
            return None;
        }
        let AtomView::Fun(representation) = *arguments.last()? else {
            return None;
        };
        if representation.get_symbol() != Minkowski::selfless_symbol()
            || representation.get_nargs() != 1
            || representation.iter().next()? != dimension.as_view()
        {
            return None;
        }
        Some(
            arguments[..arguments.len() - 1]
                .iter()
                .fold(
                    FunctionBuilder::new(function.get_symbol()),
                    |builder, argument| builder.add_arg(*argument),
                )
                .finish(),
        )
    }

    fn restore_external_vector(arguments: &[AtomView<'_>], dimension: &Atom) -> Option<Atom> {
        let [AtomView::Fun(metadata), slot] = arguments else {
            return None;
        };
        // Original vector syntax is opaque metadata, not another tensor port.
        if metadata.get_symbol() != SPENSO_TAG.pure_scalar || metadata.get_nargs() != 1 {
            return None;
        }
        let original = metadata.iter().next()?;
        let AtomView::Fun(representation) = slot else {
            return None;
        };
        if representation.get_symbol() != Minkowski::selfless_symbol()
            || !matches!(representation.get_nargs(), 1 | 2)
            || representation.iter().next()? != dimension.as_view()
        {
            return None;
        }
        if representation.get_nargs() == 1 {
            return Some(original.to_owned());
        }

        let index = representation.iter().nth(1)?;
        let AtomView::Fun(original_function) = original else {
            return None;
        };
        let original_arguments = original_function.iter().collect::<Vec<_>>();
        if original_arguments
            .last()
            .is_some_and(|argument| Self::is_compact_spenso_minkowski_representation(*argument))
        {
            let original_dimension = match *original_arguments.last().unwrap() {
                AtomView::Fun(original_representation) => original_representation.iter().next()?,
                _ => return None,
            };
            let indexed_representation = FunctionBuilder::new(Minkowski::selfless_symbol())
                .add_arg(original_dimension)
                .add_arg(index)
                .finish();
            return Some(
                original_arguments[..original_arguments.len() - 1]
                    .iter()
                    .fold(
                        FunctionBuilder::new(original_function.get_symbol()),
                        |builder, argument| builder.add_arg(*argument),
                    )
                    .add_arg(indexed_representation)
                    .finish(),
            );
        }
        Some(
            original_arguments
                .iter()
                .fold(
                    FunctionBuilder::new(original_function.get_symbol()),
                    |builder, argument| builder.add_arg(*argument),
                )
                .add_arg(index)
                .finish(),
        )
    }

    fn spenso_slot_index(slot: AtomView<'_>, dimension: &Atom) -> Option<Atom> {
        let AtomView::Fun(representation) = slot else {
            return None;
        };
        if representation.get_symbol() != Minkowski::selfless_symbol()
            || representation.get_nargs() != 2
        {
            return None;
        }
        let arguments = representation.iter().collect::<Vec<_>>();
        (arguments[0] == dimension.as_view()).then(|| arguments[1].to_owned())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VakintExpression(
    /// Expression terms, each separating a matched integral from its numerator.
    pub Vec<VakintTerm>,
);

impl fmt::Display for VakintExpression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}",
            if self.0.len() > 1 { " " } else { "" },
            self.0
                .iter()
                .map(|term| format!("{}", term))
                .collect::<Vec<_>>()
                .join("\n+")
        )
    }
}

impl VakintExpression {
    fn prepare_evaluation(
        input: AtomView,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<Self, VakintError> {
        let mut has_float = false;
        input.visitor(&mut |part| {
            if matches!(part, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("topo"))
                || precision::legacy::is_source(part)
            {
                return false;
            }
            has_float |= matches!(part, AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Float(_, _)));
            !has_float
        });
        if !has_float {
            return Self::try_from(input);
        }
        let request_id = Vakint::next_input_request_id()?;

        // Collection already combines coefficients of equal topologies. Lift
        // each finite occurrence before that algebra, then select its backend
        // using the concrete matched topology. PySecDec receives its original
        // numeric input, never these placeholders or private source leaves.
        let placeholder = (0..)
            .map(|index| vk_symbol!(format!("input_coefficient_placeholder_{index}")))
            .find(|symbol| !input.contains_symbol(*symbol))
            .unwrap();
        let mut literals = Vec::new();
        let lifted = input.replace_map(|part, _context, out| {
            if matches!(part, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("topo"))
                || precision::legacy::is_source(part)
            {
                **out = part.to_owned();
            } else if matches!(part, AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Float(_, _))) {
                let key = placeholder.call(literals.len() as i64);
                literals.push((key.clone(), part.to_owned()));
                **out = key;
            }
        });
        let mut expression = Self::try_from(lifted)?;
        for term in &mut expression.0 {
            let track = term.tracks_input_coefficients(vakint, settings)?;
            let mut replacements = Vec::new();
            for (literal_id, (key, value)) in literals.iter().enumerate() {
                if term
                    .numerator
                    .pattern_match(&key.to_pattern(), None, None)
                    .next()
                    .is_none()
                {
                    continue;
                }
                let value = if track {
                    precision::legacy::stage_table(
                        function!(
                            vk_symbol!("input_coefficient"),
                            request_id,
                            literal_id as i64
                        )
                        .as_view(),
                        value.as_view(),
                        vk_symbol!(&settings.epsilon_symbol),
                        "Input",
                    )?
                } else {
                    value.clone()
                };
                replacements.push(Replacement::new(key.clone(), value));
            }
            term.numerator = term.numerator.replace_multiple(&replacements);
        }
        expression.0.retain(|term| !term.numerator.is_zero());
        Ok(expression)
    }

    pub fn split_integrals(input: AtomView) -> Result<Vec<VakintTerm>, VakintError> {
        let mut res = vec![];

        let topo_variables = input
            .pattern_match(&vk_parse!("topo(props_)").unwrap().to_pattern(), None, None)
            .map(|m| function!(vk_symbol!("topo"), m.get(&vk_symbol!("props_")).unwrap()))
            .collect::<HashSet<_>>();

        for (integral, numerator) in
            input.coefficient_list::<i8>(&topo_variables.iter().collect::<Vec<_>>())
        {
            if integral == Atom::num(1) {
                return Err(VakintError::InvalidIntegralFormat(format!("{}", numerator)));
            }
            if utils::could_match(
                &vk_parse!("topo(props_)^n_").unwrap().to_pattern(),
                integral.as_view(),
            ) {
                return Err(VakintError::InvalidIntegralFormat(format!("{}", integral)));
            }
            res.push(VakintTerm {
                integral: integral.to_owned().clone(),
                numerator: numerator.to_owned().clone(),
                vectors: VakintTerm::identify_vectors_in_numerator(numerator.as_view())?,
            });
        }

        Ok(res)
    }

    /// Canonicalize against the supplied registry without artifact I/O. Use
    /// `Vakint::to_canonical` for settings-driven lazy artifact loading.
    pub fn canonicalize(
        &mut self,
        settings: &VakintSettings,
        topologies: &Topologies,
        short_form: bool,
    ) -> Result<(), VakintError> {
        settings.with_resource_budget(|| {
            self.canonicalize_with(settings, short_form, |input| {
                topologies.match_topologies_to_user_input(input, settings.allow_unknown_integrals)
            })
        })
    }

    fn canonicalize_with(
        &mut self,
        settings: &VakintSettings,
        short_form: bool,
        mut match_input: impl FnMut(AtomView) -> Result<Option<ReplacementRules>, VakintError>,
    ) -> Result<(), VakintError> {
        for term in self.0.iter_mut() {
            resources::checkpoint()?;
            if let Some(replacement_rules) = match_input(term.integral.as_view())? {
                //println!("replacement_rules = {}", replacement_rules,);
                term.canonicalize(settings, &replacement_rules, short_form)?;
            } else {
                return Err(VakintError::UnreckognizedIntegral(
                    term.integral.to_string(),
                ));
            }
        }
        Ok(())
    }

    pub fn tensor_reduce(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        settings.with_resource_budget(|| {
            for term in self.0.iter_mut() {
                resources::checkpoint()?;
                term.tensor_reduce(vakint, settings)?;
            }
            Ok(())
        })
    }

    pub fn evaluate_integral(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        self.evaluate_integral_observed(vakint, settings, None)
    }

    fn evaluate_integral_observed(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<(), VakintError> {
        settings.with_resource_budget(|| {
            MonitorContext::native_phase(observation.as_deref_mut(), "integral_reduction", None);
            let mut complete = Vec::new();
            let mut numerical = Atom::Zero;
            let mut numerical_vectors = std::collections::BTreeSet::new();
            let mut reductions = BTreeMap::new();
            for (source_id, term) in self.0.iter_mut().enumerate() {
                resources::checkpoint()?;
                match term.reduce_integral(
                    vakint,
                    settings,
                    source_id,
                    observation.as_deref_mut(),
                )? {
                    IntegralEvaluation::Complete { numerical: true } => {
                        numerical += &term.numerator;
                        numerical_vectors.extend(term.vectors.iter().cloned());
                    }
                    IntegralEvaluation::Complete { numerical: false } => {
                        complete.push(term.clone())
                    }
                    IntegralEvaluation::RustRed { loop_count, mode } => {
                        // Preserve the existing loop-count-specific epsilon target,
                        // but combine every matching contribution before checking
                        // missing orders, finite precision, or residual masters.
                        let (raw, vectors) = reductions
                            .entry((loop_count, mode))
                            .or_insert_with(|| (Atom::Zero, std::collections::BTreeSet::new()));
                        *raw += term.numerator.clone();
                        vectors.extend(term.vectors.iter().cloned());
                    }
                }
            }
            for ((loop_count, mode), (raw, vectors)) in reductions {
                MonitorContext::native_phase(
                    observation.as_deref_mut(),
                    match mode {
                        RustRedMasterMode::Retain => "retained_master_combination",
                        RustRedMasterMode::Symbolic => "symbolic_master_substitution",
                        RustRedMasterMode::Numerical => "numerical_master_substitution",
                    },
                    Some(raw.as_view()),
                );
                let numerator = vakint.finalize_rustred(
                    rustred_evaluation::combine_exact(raw),
                    settings,
                    loop_count,
                    mode,
                    observation.as_deref_mut(),
                )?;
                if mode == RustRedMasterMode::Numerical {
                    numerical += numerator;
                    numerical_vectors.extend(vectors);
                } else if !numerator.is_zero() {
                    complete.push(VakintTerm {
                        integral: Atom::num(1),
                        numerator,
                        vectors: vectors.into_iter().collect(),
                    });
                }
            }
            if !numerical.is_zero() {
                MonitorContext::native_phase(
                    observation,
                    "numerical_materialization",
                    Some(numerical.as_view()),
                );
                complete.push(VakintTerm {
                    integral: Atom::one(),
                    numerator: vakint.resolve_rustred_numerical(
                        rustred_evaluation::combine_exact(numerical),
                        settings,
                    )?,
                    vectors: numerical_vectors.into_iter().collect(),
                });
            }
            self.0 = complete;
            Ok(())
        })
    }

    #[allow(dead_code)]
    pub fn map<F>(&mut self, f: F) -> VakintExpression
    where
        F: Fn(&VakintTerm) -> VakintTerm,
    {
        VakintExpression(self.0.iter().map(f).collect())
    }

    #[allow(dead_code)]
    pub fn map_numerator<F>(&self, f: F) -> VakintExpression
    where
        F: Fn(AtomView) -> Atom,
    {
        VakintExpression(
            self.0
                .iter()
                .map(|term| VakintTerm {
                    integral: term.integral.clone(),
                    numerator: f(term.numerator.as_view()),
                    vectors: term.vectors.clone(),
                })
                .collect(),
        )
    }

    #[allow(dead_code)]
    pub fn map_integrals<F>(&self, f: F) -> VakintExpression
    where
        F: Fn(AtomView) -> Atom,
    {
        VakintExpression(
            self.0
                .iter()
                .map(|term| VakintTerm {
                    integral: f(term.integral.as_view()),
                    numerator: term.numerator.clone(),
                    vectors: term.vectors.clone(),
                })
                .collect(),
        )
    }
}

impl TryFrom<Atom> for VakintExpression {
    type Error = VakintError;
    fn try_from(atom: Atom) -> Result<Self, VakintError> {
        Ok(VakintExpression(VakintExpression::split_integrals(
            atom.as_view(),
        )?))
    }
}

impl TryFrom<AtomView<'_>> for VakintExpression {
    type Error = VakintError;
    fn try_from(atom_view: AtomView) -> Result<Self, VakintError> {
        Ok(VakintExpression(VakintExpression::split_integrals(
            atom_view,
        )?))
    }
}

impl From<VakintExpression> for Atom {
    fn from(vakint_expr: VakintExpression) -> Atom {
        let mut res = Atom::Zero;
        for term in vakint_expr.0.iter() {
            let t: Atom = VakintTerm::into(term.clone());
            res += t;
        }
        res
    }
}

impl Vakint {
    /// Inspect a registered RustRed constant in this engine's active data snapshot.
    pub fn get_constant(
        &self,
        symbol: &ConstantSymbol,
        prec: u32,
    ) -> Result<VakintConstant, VakintError> {
        self.rustred_data()?
            .masters
            .get_constant(symbol.clone(), prec)
    }

    /// Inspect every available Laurent coefficient and the native master graph.
    /// A bare IBP application announces identities only; its matching data pack
    /// is required for authenticated graph and coefficient inspection.
    pub fn get_master(
        &self,
        symbol: &MasterSymbol,
        prec: u32,
    ) -> Result<VakintMaster, VakintError> {
        let result = self
            .rustred_data()?
            .masters
            .get_master(symbol.clone(), prec);
        if matches!(
            &result,
            Err(VakintError::Inspection(
                inspection::InspectionError::UnknownMaster { .. }
            ))
        ) && let Some(declaration) = self.announced_bare_masters()?.get(symbol.as_atom())
        {
            return Err(inspection::InspectionError::UnavailableGraph {
                master: symbol.to_string(),
                detail: format!(
                    "declared by native application {}; load its matching external-data pack to inspect authenticated master graphs and coefficients",
                    declaration.application_path.display(),
                ),
            }
            .into());
        }
        result
    }

    /// Registered constants used by the built-in or explicitly loaded RustRed data.
    pub fn constant_symbols(&self) -> Result<Vec<ConstantSymbol>, VakintError> {
        self.rustred_data()?.masters.constant_symbols()
    }

    /// Family-qualified identities, including declared but uncomputed terminals.
    pub fn master_symbols(&self) -> Result<Vec<MasterSymbol>, VakintError> {
        let mut symbols = self.rustred_data()?.masters.master_symbols();
        for symbol in self.announced_bare_masters()?.into_keys() {
            symbols.push(MasterSymbol::try_from(symbol)?);
        }
        symbols.sort();
        symbols.dedup();
        Ok(symbols)
    }

    /// Normalize FORM's Gamma and EulerGamma spellings before the complete
    /// expression is expanded. Native Laurent series then compute every order
    /// needed by surrounding poles, without introducing an artificial Gamma tail.
    pub(crate) fn normalize_gamma_functions(expression: AtomView) -> Atom {
        expression
            .replace(vk_parse!("Gamma(x_)").unwrap().to_pattern())
            .with(function!(symbolica::transcendental::gamma(), Atom::var(S.x_)).to_pattern())
            .replace(vk_parse!("EulerGamma").unwrap().to_pattern())
            .with(Atom::var(symbolica::transcendental::euler_gamma()).to_pattern())
    }

    pub fn initialize_vakint_symbols() {
        // Force initialization of symbolica symbols with proper attributes
        LazyLock::force(&S);
    }

    pub fn new() -> Result<Self, VakintError> {
        citations::USED_CITATIONS.record(citations::CitationSource::Vakint);
        // Force initialization of symbolica symbols with proper attributes
        LazyLock::force(&S);

        let started = std::time::Instant::now();
        let topologies = Topologies::generate_topologies()?;
        log::debug!(target: "vakint::saved_application", "phase=base_catalogue_prepared elapsed_seconds={:.6}", started.elapsed().as_secs_f64());

        Ok(Vakint {
            topologies,
            rustred_application: Arc::default(),
            rustred_lower_caches: Arc::default(),
            rustred_data: None,
        })
    }

    /// Match existing lower/custom data first. Only a five-loop candidate miss
    /// with an explicitly configured artifact enters the saved loader, before
    /// the unknown fallback. Unrelated unknown lower loops do not open it.
    /// Registry-only public APIs intentionally have no implicit loader authority.
    pub(crate) fn match_topologies_for_settings(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        if let Some(matched) = self
            .topologies
            .match_topologies_to_user_input(input, false)?
        {
            return Ok(Some(matched));
        }
        let needs_five_loop = Topologies::needs_five_loop_catalogue(input)?;
        if (settings.rustred_artifact_path.is_some() || self.has_rustred_loop(5)) && needs_five_loop
        {
            self.saved_application(settings)?;
            if let Some(matched) = self
                .topologies
                .match_topologies_to_user_input(input, false)?
            {
                return Ok(Some(matched));
            }
        }
        if needs_five_loop
            && settings.rustred_artifact_path.is_none()
            && !self.has_rustred_loop(5)
            && settings
                .evaluation_order
                .0
                .iter()
                .any(|method| matches!(method, EvaluationMethod::RustRed(_)))
        {
            return Err(RustRedEvaluationError::ExternalData {
                detail: "5L RustRed data is not loaded; call load_official_external_data(5, ...) or provide a local application package".into(),
            }.into());
        }
        if settings.allow_unknown_integrals {
            self.topologies.match_topologies_to_user_input(input, true)
        } else {
            Ok(None)
        }
    }

    pub fn validate_settings(&self, settings: &VakintSettings) -> Result<(), VakintError> {
        // Verify that the chosen normalisation only contains the expected symbols
        let (_full_atom, expanded, evaluated) =
            settings.integral_normalization_factor.validate(settings)?;

        debug!(
            "Loop normalisation factor considered:
Full                          : {}
Expanded (n_loops=1)          : {}
Evaluated (n_loops=1, mu_r=1) :
{}",
            settings.integral_normalization_factor, expanded, evaluated
        );

        if settings
            .evaluation_order
            .0
            .iter()
            .any(|em| em.dependencies().contains(&VakintDependency::FORM))
        {
            let form_version = self.get_form_version(settings)?;
            match compare_to(form_version.clone(), MINIMAL_FORM_VERSION, Cmp::Ge) {
                Ok(valid) => {
                    if valid {
                        debug!(
                            "{} successfully detected with version '{}'",
                            "FORM".green(),
                            form_version.green()
                        );
                    } else {
                        return Err(VakintError::FormVersion(format!(
                            "{} version installed on your system does not meet minimal requirements: {}<{}",
                            "FORM".red(),
                            form_version.red(),
                            MINIMAL_FORM_VERSION
                        )));
                    }
                }
                Err(_) => {
                    return Err(VakintError::FormVersion(format!(
                        "Could not parse {} version '{}'.",
                        "FORM".red(),
                        form_version.red()
                    )));
                }
            };
        }

        if settings
            .evaluation_order
            .0
            .iter()
            .any(|em| em.dependencies().contains(&VakintDependency::PySecDec))
        {
            let pysecdec_version = self.get_pysecdec_version(settings)?;
            match compare_to(pysecdec_version.clone(), MINIMAL_PYSECDEC_VERSION, Cmp::Ge) {
                Ok(valid) => {
                    if valid {
                        debug!(
                            "{} successfully detected with version '{}'",
                            "PySecDec".green(),
                            pysecdec_version.green()
                        );
                    } else {
                        return Err(VakintError::FormVersion(format!(
                            "{} version installed on your system does not meet minimal requirements: {}<{}.",
                            "PySecDec".red(),
                            pysecdec_version.red(),
                            MINIMAL_PYSECDEC_VERSION
                        )));
                    }
                }
                Err(_) => {
                    return Err(VakintError::FormVersion(format!(
                        "Could not parse {} version '{}'.",
                        "PySecDec".red(),
                        pysecdec_version.red()
                    )));
                }
            };
        }

        Ok(())
    }

    pub fn params_from_f64(
        &self,
        settings: &VakintSettings,
        params: &HashMap<String, f64>,
    ) -> HashMap<String, Float, ahash::RandomState> {
        let hm = HashMap::<String, f64, ahash::RandomState>::from_iter(
            params.iter().map(|(k, v)| (k.clone(), *v)),
        );
        params_from_f64(&hm, settings.run_time_decimal_precision)
    }

    pub fn params_from_complex_f64(
        &self,
        settings: &VakintSettings,
        params: &HashMap<String, Complex<f64>>,
    ) -> HashMap<String, Complex<Float>, ahash::RandomState> {
        let hm = HashMap::<String, Complex<f64>, ahash::RandomState>::from_iter(
            params.iter().map(|(k, v)| (k.clone(), *v)),
        );
        params_from_complex_f64(&hm, settings.run_time_decimal_precision)
    }

    pub fn externals_from_f64(
        &self,
        settings: &VakintSettings,
        externals: &HashMap<usize, (f64, f64, f64, f64)>,
    ) -> HashMap<usize, Momentum, ahash::RandomState> {
        let em = HashMap::<usize, (f64, f64, f64, f64), ahash::RandomState>::from_iter(
            externals.iter().map(|(k, v)| (*k, *v)),
        );
        externals_from_f64(&em, settings.run_time_decimal_precision)
    }

    #[allow(clippy::type_complexity)]
    pub fn externals_from_complex_f64(
        &self,
        settings: &VakintSettings,
        externals: &HashMap<usize, (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>)>,
    ) -> HashMap<usize, Momentum, ahash::RandomState> {
        let em = HashMap::<
            usize,
            (Complex<f64>, Complex<f64>, Complex<f64>, Complex<f64>),
            ahash::RandomState,
        >::from_iter(externals.iter().map(|(k, v)| (*k, *v)));
        externals_from_complex_f64(&em, settings.run_time_decimal_precision)
    }

    pub fn numerical_evaluation(
        &self,
        settings: &VakintSettings,
        expression: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
        let data = self.rustred_data()?;
        let (expression, parameters) = precision::evaluation_inputs(
            expression,
            vk_symbol!(&settings.epsilon_symbol),
            params_real,
            params_complex,
            externals,
        )?;
        self.require_bare_master_coefficients(&expression)?;
        precision::evaluate_with_errors(
            expression.as_view(),
            &parameters,
            settings,
            Some(&data.masters),
        )
    }

    pub fn get_coeff_map(prec: u32) -> impl Fn(&Fraction<IntegerRing>) -> Complex<Float> {
        move |x| Complex::new(x.to_multi_prec_float(prec), Float::with_val(prec, 0.))
    }

    #[allow(clippy::complexity)]
    pub fn get_constants_map(
        settings: &VakintSettings,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<
        (
            HashMap<Atom, Float, ahash::random_state::RandomState>,
            HashMap<Atom, Complex<Float>, ahash::random_state::RandomState>,
        ),
        VakintError,
    > {
        let mut const_map_real: HashMap<Atom, Float, ahash::random_state::RandomState> =
            HashMap::default();
        let mut const_map_complex: HashMap<Atom, Complex<Float>, ahash::random_state::RandomState> =
            HashMap::default();

        let binary_prec = settings.get_binary_precision();
        const_map_real.insert(
            Atom::from(Var::new(Symbol::PI)),
            Float::with_val(binary_prec, Constant::Pi),
        );

        if let Some(ext_p) = externals.as_ref() {
            for (&i, val1) in ext_p.iter() {
                for (&j, val2) in ext_p.iter() {
                    if i > j {
                        continue;
                    }
                    let key = function!(
                        S.dot,
                        function!(S.p, Atom::num(i as i64)),
                        function!(S.p, Atom::num(j as i64))
                    );
                    match (val1, val2) {
                        (Momentum::Real(k1), Momentum::Real(k2)) => {
                            let dot_product = k1.0.to_owned() * k2.0.to_owned()
                                - k1.1.to_owned() * k2.1.to_owned()
                                - k1.2.to_owned() * k2.2.to_owned()
                                - k1.3.to_owned() * k2.3.to_owned();
                            const_map_real.insert(key, dot_product.to_owned());
                        }
                        (Momentum::Real(k1), Momentum::Complex(k2)) => {
                            let dot_product = k2.0.to_owned() * k1.0.to_owned()
                                - k2.1.to_owned() * k1.1.to_owned()
                                - k2.2.to_owned() * k1.2.to_owned()
                                - k2.3.to_owned() * k1.3.to_owned();
                            const_map_complex.insert(key, dot_product.to_owned());
                        }
                        (Momentum::Complex(k1), Momentum::Real(k2)) => {
                            let dot_product = k1.0.to_owned() * k2.0.to_owned()
                                - k1.1.to_owned() * k2.1.to_owned()
                                - k1.2.to_owned() * k2.2.to_owned()
                                - k1.3.to_owned() * k2.3.to_owned();
                            const_map_complex.insert(key, dot_product.to_owned());
                        }
                        (Momentum::Complex(k1), Momentum::Complex(k2)) => {
                            let dot_product = k1.0.to_owned() * k2.0.to_owned()
                                - k1.1.to_owned() * k2.1.to_owned()
                                - k1.2.to_owned() * k2.2.to_owned()
                                - k1.3.to_owned() * k2.3.to_owned();
                            const_map_complex.insert(key, dot_product.to_owned());
                        }
                    }
                }
            }
        }

        const_map_real.insert(
            Atom::from(Var::new(vk_symbol!("EulerGamma"))),
            Float::with_val(binary_prec, Constant::Euler),
        );

        const_map_complex.insert(
            Atom::from(Var::new(S.cmplx_i)),
            Complex::new(
                Float::with_val(binary_prec, 0),
                Float::with_val(binary_prec, 1),
            ),
        );

        const_map_real.insert(
            function!(Symbol::LOG, Atom::num(2)),
            Float::with_val(binary_prec, Constant::Log2),
        );

        for (symb, value) in params_complex.iter() {
            const_map_complex.insert(vk_parse!(symb.as_str()).unwrap(), value.clone());
        }
        for (symb, value) in params_real.iter() {
            const_map_real.insert(vk_parse!(symb.as_str()).unwrap(), value.clone());
        }

        Ok((const_map_real, const_map_complex))
    }

    pub fn partial_numerical_evaluation(
        settings: &VakintSettings,
        integral: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<Atom, VakintError> {
        let (expression, parameters) = precision::evaluation_inputs(
            integral,
            vk_symbol!(&settings.epsilon_symbol),
            params_real,
            params_complex,
            externals,
        )?;
        precision::materialize(expression.as_view(), settings, &parameters, None)
    }

    pub fn full_numerical_evaluation(
        settings: &VakintSettings,
        integral: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
        let (expression, parameters) = precision::evaluation_inputs(
            integral,
            vk_symbol!(&settings.epsilon_symbol),
            params_real,
            params_complex,
            externals,
        )?;
        precision::evaluate_with_errors(expression.as_view(), &parameters, settings, None)
    }

    pub fn full_numerical_evaluation_without_error(
        settings: &VakintSettings,
        integral: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<NumericalEvaluationResult, VakintError> {
        let (expression, parameters) = precision::evaluation_inputs(
            integral,
            vk_symbol!(&settings.epsilon_symbol),
            params_real,
            params_complex,
            externals,
        )?;
        precision::evaluate_series(
            expression.as_view(),
            vk_symbol!(&settings.epsilon_symbol),
            &parameters,
            settings,
            None,
        )
    }

    fn pysecdec_evaluate(
        &self,
        settings: &VakintSettings,
        input_numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &PySecDecOptions,
        source_id: usize,
    ) -> Result<Atom, VakintError> {
        let mut has_retained_source = false;
        input_numerator.visitor(&mut |part| {
            has_retained_source |= precision::legacy::is_source(part);
            !has_retained_source
        });
        if has_retained_source {
            return Err(VakintError::InvalidNumerator(
                "PySecDec cannot propagate the uncertainty of retained NumericalSource coefficients; use a FORM or RustRed evaluation backend for this input".into(),
            ));
        }
        let integral = integral_specs.canonical_topology.get_integral();

        let pysecdec_inputs = if options.reuse_existing_output.is_some()
            && PathBuf::from(options.reuse_existing_output.as_ref().unwrap()).exists()
        {
            vec![("run_pySecDec.py".into(), "".into())]
        } else {
            debug!(
                "Processing the following integral with {}:\n{}",
                "PySecDec".green(),
                integral
            );
            let dot_product_numerator = Vakint::convert_from_dot_notation(input_numerator);
            let vectors =
                VakintTerm::identify_vectors_in_numerator(dot_product_numerator.as_view())?;
            let numerator_atom = Vakint::convert_to_dot_notation(input_numerator);
            let numerator = numerator_atom.as_view();
            let mut processed_numerator = Vakint::convert_to_dot_notation(numerator);

            // Make sure there is no open index left
            if processed_numerator
                .pattern_match(
                    &vk_parse!("s_(id_,idx_)").unwrap().to_pattern(),
                    Some(
                        &(Condition::from((vk_symbol!("id_"), number_condition()))
                            & Condition::from((vk_symbol!("s_"), symbol_condition()))),
                    ),
                    None,
                )
                .next()
                .is_some()
            {
                return Err(VakintError::InvalidNumerator(format!(
                    "PySecDec can only handle scalar numerator. If you have open indices, make sure they are contracted with external momenta: {}",
                    processed_numerator
                )));
            }

            // Check if numerator contains additional symbols
            // First, replace functions with 1 and get all remaining symbols
            let mut numerator_additional_symbols = input_numerator
                .replace(vk_parse!("f_(args__)").unwrap().to_pattern())
                .with(vk_parse!("1").unwrap().to_pattern())
                .get_all_symbols(false);
            let eps_symbol: Symbol = vk_symbol!(settings.epsilon_symbol.clone());
            numerator_additional_symbols.retain(|&s| s != eps_symbol);

            // Convert back from dot notation
            processed_numerator = Vakint::convert_from_dot_notation(processed_numerator.as_view());

            let mut lorentz_indices = HashSet::<String>::new();
            let arc_mutex_lorentz_indices = Arc::new(Mutex::new(HashSet::new()));

            let dot_product_matcher = vk_parse!("v1_(id1_,idx_)*v2_(id2_,idx_)").unwrap();
            // Powers higher than two cannot occur as different dummy indices would have been used in
            // the call 'processed_numerator = Vakint::convert_from_dot_notation(processed_numerator.as_view(), true)'
            let square_matcher = vk_parse!("v1_(id1_,idx_)^2").unwrap();
            let mut old_processed_numerator = processed_numerator.clone();
            loop {
                let arc_mutex_lorentz_indices_sent = arc_mutex_lorentz_indices.clone();
                let dot_product_transformer = Transformer::Map(Box::new(
                    move |a_in: AtomView, _state: _, a_out: &mut Atom| {
                        if let AtomView::Fun(s) = a_in {
                            let a_in = s
                                .to_slice()
                                .iter()
                                .map(|a| {
                                    AtomPrinter::new_with_options(
                                        a,
                                        PrintOptions::file_no_namespace(),
                                    )
                                    .to_string()
                                })
                                .collect::<Vec<_>>();
                            arc_mutex_lorentz_indices_sent
                                .lock()
                                .unwrap()
                                .insert(format!("mu{}", a_in[4]));
                            *a_out = vk_parse!(
                                format!(
                                    "{}{}(mu{})*{}{}(mu{})",
                                    a_in[0], a_in[1], a_in[4], a_in[2], a_in[3], a_in[4],
                                )
                                .as_str()
                            )
                            .unwrap();
                        };

                        Ok(())
                    },
                ));

                processed_numerator = old_processed_numerator
                    .replace(dot_product_matcher.to_pattern())
                    .rhs_cache_size(0)
                    .with(Pattern::Transformer(Box::new((
                        Some(
                            vk_parse!("arg(v1_,id1_,v2_,id2_,idx_)")
                                .unwrap()
                                .to_pattern(),
                        ),
                        vec![dot_product_transformer.clone()],
                    ))));

                processed_numerator = processed_numerator
                    .replace(square_matcher.to_pattern())
                    .rhs_cache_size(0)
                    .with(Pattern::Transformer(Box::new((
                        Some(
                            vk_parse!("arg(v1_,id1_,v1_,id1_,idx_)")
                                .unwrap()
                                .to_pattern(),
                        ),
                        vec![dot_product_transformer],
                    ))));
                lorentz_indices.extend(arc_mutex_lorentz_indices.lock().unwrap().clone());
                if old_processed_numerator == processed_numerator {
                    break;
                } else {
                    old_processed_numerator = processed_numerator.clone();
                }
            }

            let mut m = Atom::new();
            let power_list_map = integral_specs.get_propagator_property_list("pow_");

            let mass_list_map = integral_specs.get_propagator_property_list("mUVsq_");
            let mut masses = HashSet::new();
            let mut power_list: Vec<i64> = vec![];

            let mut loop_momenta_ids_found = HashSet::new();
            let pysecdec_propagators = format!(
                "[{}]",
                integral
                    .graph
                    .edges
                    .iter()
                    .map(|(e_id, e)| format!(
                        "'({})**2{}'",
                        {
                            m = e.momentum.clone();
                            power_list.push(
                                power_list_map
                                    .get(e_id)
                                    .unwrap()
                                    .to_owned()
                                    .try_into()
                                    .unwrap(),
                            );
                            for i_loop in 1..=integral.n_loops {
                                let p = vk_parse!(format!("k({})", i_loop).as_str()).unwrap();
                                if utils::could_match(&p.to_pattern(), m.as_view()) {
                                    loop_momenta_ids_found.insert(i_loop);
                                }
                                m = m.replace(p.to_pattern()).with(
                                    vk_parse!(format!("k{}", i_loop).as_str())
                                        .unwrap()
                                        .to_pattern(),
                                );
                            }
                            AtomPrinter::new_with_options(
                                m.as_view(),
                                PrintOptions::file_no_namespace(),
                            )
                        },
                        {
                            //let m = mass_list_map.get(e_id).unwrap().to_owned();
                            let (m, m_sq) = match mass_list_map.get(e_id).unwrap().to_owned() {
                                Atom::Var(s) => (
                                    Atom::pow(&Atom::Var(s.clone()), Atom::num(1) / Atom::num(2)),
                                    Atom::Var(s.clone()),
                                ),
                                Atom::Pow(m) => {
                                    let base = m.to_pow_view().get_base().to_owned();
                                    let exp = m.to_pow_view().get_exp();
                                    match base {
                                        Atom::Var(s) if exp == Atom::num(2).as_view() => {
                                            (Atom::Var(s.clone()), Atom::Pow(m.clone()))
                                        }

                                        _ => {
                                            panic!(
                                                "Could not find mass symbol in integral:\n{}",
                                                integral
                                            );
                                        }
                                    }
                                }
                                _ => {
                                    panic!("Could not find mass symbol in integral:\n{}", integral);
                                }
                            };

                            if m.is_zero() {
                                "".into()
                            } else {
                                let mass_symbol_to_add = match (m.clone(), m_sq.clone()) {
                                    (Atom::Var(a), _) => a.get_symbol(),
                                    (_, Atom::Var(asq)) => asq.get_symbol(),
                                    _ => panic!(
                                        "Could not find mass symbol for integral:\n{}",
                                        integral
                                    ),
                                };
                                masses.insert(pysecdec_encode(&undress_vakint_symbols(
                                    &get_full_name(&mass_symbol_to_add),
                                )));
                                format!(
                                    "-{}",
                                    pysecdec_encode(&undress_vakint_symbols(
                                        &m_sq.to_canonical_string()
                                    ))
                                )
                            }
                        }
                    ))
                    .collect::<Vec<_>>()
                    .join(", "),
            );
            let n_loops_in_topology = loop_momenta_ids_found.len();

            // This should always hold, because pinched higher loops that reduce the loop counts must be captured by lower loop count topologies.
            // It is important that this hold otherwise MSbar normalization factor would be off.
            assert!(n_loops_in_topology == integral.n_loops);

            let mut vars: HashMap<String, String> = HashMap::new();

            vars.insert("graph_name".into(), "pySecDecRun".into());

            vars.insert("propagators".into(), pysecdec_propagators);
            let mut sorted_lorentz_indices = lorentz_indices.iter().cloned().collect::<Vec<_>>();
            sorted_lorentz_indices.sort();
            let mut lorentz_index_replacements = vec![];
            let mut sanitized_lorentz_indices = vec![];
            let mut used_sanitized_lorentz_indices = HashSet::new();
            for original_index in sorted_lorentz_indices {
                let mut sanitized_index = original_index
                    .chars()
                    .map(|c| {
                        if c.is_ascii_alphanumeric() || c == '_' {
                            c
                        } else {
                            '_'
                        }
                    })
                    .collect::<String>();
                if sanitized_index.is_empty() || sanitized_index.as_bytes()[0].is_ascii_digit() {
                    sanitized_index = format!("idx_{}", sanitized_index);
                }
                let sanitized_index_base = sanitized_index.clone();
                let mut collision_counter = 1;
                while !used_sanitized_lorentz_indices.insert(sanitized_index.clone()) {
                    collision_counter += 1;
                    sanitized_index = format!("{}_{}", sanitized_index_base, collision_counter);
                }
                lorentz_index_replacements.push((original_index, sanitized_index.clone()));
                sanitized_lorentz_indices.push(sanitized_index);
            }
            vars.insert(
                "lorentz_indices".into(),
                format!(
                    "[{}]",
                    sanitized_lorentz_indices
                        .iter()
                        .map(|li| format!("'{}'", li))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            vars.insert(
                "loop_momenta".into(),
                format!(
                    "[{}]",
                    (1..=integral.n_loops)
                        .map(|i| format!("'k{}'", i))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            let mut external_momenta = vectors
                .iter()
                .filter_map(|(v, id)| {
                    if v == EXTERNAL_MOMENTUM_SYMBOL {
                        Some(format!("{}{}", v, id))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            external_momenta.sort();
            vars.insert(
                "external_momenta".into(),
                format!(
                    "[{}]",
                    external_momenta
                        .iter()
                        .map(|em| format!("'{}'", em))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            vars.insert(
                "power_list".into(),
                format!(
                    "tuple([{}])",
                    power_list
                        .iter()
                        .map(|p| format!("{}", p))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            );
            let numerator_path = String::from("numerator.txt");
            vars.insert("numerator_path".into(), numerator_path.clone());
            // Expand the numerator around epsilon=0 to make sure it is polynomial
            processed_numerator = processed_numerator
                .series(
                    vk_symbol!(settings.epsilon_symbol.as_str()),
                    Atom::Zero.as_atom_view(),
                    Rational::from(
                        settings.number_of_terms_in_epsilon_expansion - (integral.n_loops as i64),
                    ),
                )
                .unwrap()
                .to_atom();
            // let mut numerator_string = AtomPrinter::new_with_options(
            //     processed_numerator.as_atom_view(),
            //     PrintOptions::file_no_namespace(),
            // )
            // .to_string()
            // .replace(settings.epsilon_symbol.as_str(), "eps");
            let mut numerator_string = pysecdec_encode(
                &undress_vakint_symbols(&processed_numerator.to_canonical_string()).replace(
                    &undress_vakint_symbols(settings.epsilon_symbol.as_str()),
                    "eps",
                ),
            );
            for (original_index, sanitized_index) in lorentz_index_replacements.iter() {
                if original_index != sanitized_index {
                    let encoded_original_index = pysecdec_encode(original_index);
                    let index_replacement_re = Regex::new(
                        format!(
                            r"(^|[^A-Za-z0-9_])({})([^A-Za-z0-9_]|$)",
                            regex::escape(&encoded_original_index)
                        )
                        .as_str(),
                    )
                    .unwrap();
                    numerator_string = index_replacement_re
                        .replace_all(numerator_string.as_str(), |caps: &regex::Captures<'_>| {
                            format!("{}{}{}", &caps[1], sanitized_index, &caps[3])
                        })
                        .to_string();
                }
            }

            // Powers higher than two cannot occur as different dummy indices would have been used in
            // the call 'processed_numerator = Vakint::convert_from_dot_notation(processed_numerator.as_view(), true)'
            let remove_squares_re =
                Regex::new(r"(?<vec>[\w|-|\d]+)\((?<idx>mu-?\d+)\)\^2").unwrap();
            numerator_string = remove_squares_re
                .replace_all(numerator_string.as_str(), "$vec$id($idx)*$vec$id($idx)")
                .to_string();
            // println!(
            //     "Numerator to be written to PySecDec input file:\n{}",
            //     numerator_string
            // );
            let description_path = "description.txt".into();
            let description_string = format!(
                "Integral:\n{}\nNumerator:\n{}",
                integral_specs,
                numerator.to_owned()
            );

            vars.insert("n_loops".into(), format!("{}", n_loops_in_topology));
            vars.insert("loop_additional_prefactor".into(), "1".into());
            vars.insert("additional_overall_factor".into(), "1.0".into());

            let mut masses_vec = masses.iter().cloned().collect::<Vec<_>>();
            masses_vec.sort();

            let mut real_parameters = vec![];
            let mut replacement_rules = vec![];

            // Add masses
            for m in masses_vec.iter() {
                real_parameters.push(format!("'{}'", m.clone()));
            }

            //Potential extra parameters
            let mut sorted_additional_numerator_symbols =
                numerator_additional_symbols.iter().collect::<Vec<_>>();
            sorted_additional_numerator_symbols.sort();
            for additional_param in sorted_additional_numerator_symbols.iter() {
                real_parameters.push(format!(
                    "'{}'",
                    pysecdec_encode(&undress_vakint_symbols(&get_full_name(additional_param)))
                ));
            }

            // And the external momenta
            for iv1 in 0..external_momenta.len() {
                for iv2 in iv1..external_momenta.len() {
                    let v1 = external_momenta[iv1].clone();
                    let v2 = external_momenta[iv2].clone();
                    replacement_rules.push((format!("'{}*{}'", v1, v2), format!("'{}{}'", v1, v2)));
                    real_parameters.push(format!("'{}{}'", v1, v2));
                }
            }

            vars.insert(
                "replacement_rules".into(),
                format!(
                    "[{}]",
                    replacement_rules
                        .iter()
                        .map(|(k, v)| format!("({},{})", k, v))
                        .collect::<Vec<_>>()
                        .join(",\n")
                ),
            );

            vars.insert(
                "real_parameters".into(),
                format!("[{}]", real_parameters.join(",")),
            );

            vars.insert("complex_parameters_input".into(), "[]".into());
            vars.insert("real_parameters_input".into(), "[]".into());
            let mut default_external_momenta = vec![];
            for p in external_momenta {
                if let Some(f) = options.numerical_external_momenta.get(&p) {
                    default_external_momenta.push((p, f));
                } else {
                    return Err(VakintError::EvaluationError(format!(
                        "Missing specification of numerical value for external momentum '{}'. Specify it in the PySecDecOptions of Vakint.",
                        p
                    )));
                }
            }
            vars.insert(
                "default_externals".into(),
                format!(
                    r"[{}\n]",
                    default_external_momenta
                        .iter()
                        .map(|(pi, f)| format!(
                            "({:.e},{:.e},{:.e},{:.e}) #{}",
                            f.0, f.1, f.2, f.3, pi
                        ))
                        .collect::<Vec<_>>()
                        .join(r"\n,")
                ),
            );
            let mut default_masses = vec![];
            // println!(
            //     "options.numerical_masses: {}",
            //     options
            //         .numerical_masses
            //         .iter()
            //         .map(|(k, v)| format!("{}: {}", k, v))
            //         .collect::<Vec<_>>()
            //         .join(","),
            // );
            for m in masses_vec {
                let cooked_m_symbol = undress_vakint_symbols(&pysecdec_decode(&m));
                let alternative_form = to_symbol(&cooked_m_symbol)?
                    .get_name()
                    .replace("::{}::", "::");
                // Allow for both string representation: namespace::symbol_name *and namespace::{<attributes>}::symbol_name
                // println!("Looking for mass symbol: {}", cooked_m_symbol);
                // println!(
                //     "Also trying: {}",
                //     alternative_form
                // );
                let entry = options.numerical_masses.get(&cooked_m_symbol).map_or(
                    options
                        .numerical_masses
                        .get(to_symbol(&alternative_form.clone())?.get_name()),
                    Some,
                );
                // println!("Accessing: {}", cooked_m_symbol);
                if let Some(num_m) = entry {
                    default_masses.push((m, num_m));
                } else {
                    return Err(VakintError::EvaluationError(format!(
                        "Missing specification of numerical value for mass '{}' or '{}'. Specify it in the PySecDecOptions of Vakint.",
                        cooked_m_symbol, alternative_form
                    )));
                }
            }
            for additional_param in sorted_additional_numerator_symbols.iter() {
                let param_first_form = undress_vakint_symbols(&get_full_name(additional_param));
                let alternative_form = vk_symbol!(&param_first_form.clone())
                    .get_name()
                    .replace("::{}::", "::");
                // Allow for both string representation: namespace::symbol_name *and namespace::{<attributes>}::symbol_name
                // println!(
                //     "Looking for additional numerator symbol: {}",
                //     param_first_form
                // );
                // println!("Also trying: {}",
                //     alternative_form
                // );
                let entry = options
                    .numerical_masses
                    .get(&param_first_form)
                    .map_or(options.numerical_masses.get(&alternative_form), Some);

                if let Some(num_additional_param) = entry {
                    default_masses.push((
                        pysecdec_encode(&undress_vakint_symbols(&get_full_name(additional_param))),
                        num_additional_param,
                    ));
                } else {
                    return Err(VakintError::EvaluationError(format!(
                        "Missing specification of numerical value for additional numerator symbol '{}' or '{}'. Specify it in the PySecDecOptions of Vakint.",
                        param_first_form, alternative_form
                    )));
                }
            }

            vars.insert(
                "default_masses".into(),
                format!(
                    r"[{}\n]",
                    default_masses
                        .iter()
                        .map(|(m, num_m)| format!("{:.e} #{}", num_m, m))
                        .collect::<Vec<_>>()
                        .join(r"\n,")
                ),
            );

            vars.insert(
                "max_epsilon_order".into(),
                format!(
                    "{}",
                    settings.number_of_terms_in_epsilon_expansion - (integral.n_loops as i64) - 1
                ),
            );
            vars.insert("contour_deformation".into(), "False".into());

            vars.insert(
                "drawing_input_power_list".into(),
                vars.get("power_list").unwrap().clone(),
            );
            vars.insert("drawing_input_external_lines".into(), "[]".into());
            let internal_lines_str = format!(
                "[{}]",
                integral
                    .graph
                    .edges
                    .iter()
                    .map(|(e_id, e)| {
                        format!("['e{}',[{},{}]]", e_id, e.left_node_id, e.right_node_id)
                    })
                    .collect::<Vec<_>>()
                    .join(",")
            );
            vars.insert("drawing_input_internal_lines".into(), internal_lines_str);

            vars.insert("couplings_prefactor".into(), "'1'".into());
            vars.insert("couplings_values".into(), "{}".into());

            let template =
                Template::parse_template(TEMPLATES.get("run_pySecDec_template.txt").unwrap())
                    .unwrap();

            let rendered = template
                .render(&RenderOptions {
                    variables: vars,
                    ..Default::default()
                })
                .unwrap_or_else(|e| {
                    panic!("Error while rendering template: {}", e);
                });

            vec![
                ("run_pySecDec.py".into(), rendered),
                (numerator_path, numerator_string),
                (description_path, description_string),
            ]
        };

        let mut pysecdec_options = vec!["-i".into(), "qmc".into()];
        pysecdec_options.push("--eps_rel".into());
        pysecdec_options.push(format!("{:.e}", options.relative_precision));
        pysecdec_options.push("--min_n".into());
        pysecdec_options.push(format!("{}", options.min_n_evals));
        pysecdec_options.push("--max_eval".into());
        pysecdec_options.push(format!("{}", options.max_n_evals));
        if options.quiet {
            pysecdec_options.push("-q".into());
        }
        let pysecdec_output = self.run_pysecdec(
            settings,
            &pysecdec_inputs,
            pysecdec_options,
            settings.clean_tmp_dir && options.reuse_existing_output.is_none(),
            options.reuse_existing_output.clone(),
            settings.temporary_directory.clone(),
        )?;

        if !options.quiet {
            info!(
                "PySecDec raw output:\ncentral: {}\nerror : {}",
                pysecdec_output[0], pysecdec_output[1]
            );
        }

        /*
        let pysecdec_output = [
            String::from("ep^(-1)*(2.0000000000000000e+00+(0.0000000000000000e+00)*I)+ep^(0)*(-5.4072569092295630e-01+(0.0000000000000000e+00)*I)+ep^(1)*(2.7180301350542537e+00+(0.0000000000000000e+00)*I)+ep^(2)*(-8.5638398947489147e-01+(0.0000000000000000e+00)*I)"),
            String::from("ep^(-1)*(7.9539672483176226e-17+(0.0000000000000000e+00)*I)+ep^(0)*(5.9466474629934134e-17+(0.0000000000000000e+00)*I)+ep^(1)*(1.1644579020645423e-16+(0.0000000000000000e+00)*I)+ep^(2)*(8.0813279157421749e-17+(0.0000000000000000e+00)*I)")
        ];
        */
        Self::process_pysecdec_output(
            settings,
            integral.n_loops as i64,
            source_id,
            &pysecdec_output,
        )
    }

    /// Preserve the independent finite coefficients before normalization mixes
    /// epsilon orders. The numerical integrator itself never sees private leaves.
    fn process_pysecdec_output(
        settings: &VakintSettings,
        n_loops: i64,
        source_id: usize,
        output: &[String],
    ) -> Result<Atom, VakintError> {
        if output.len() != 2 {
            return Err(VakintError::EvaluationError(format!(
                "PySecDec returned {} streams; expected central values and statistical errors",
                output.len(),
            )));
        }
        let log_mu_sq = function!(
            Symbol::LOG,
            vk_parse!(settings.mu_r_sq_symbol.as_str()).unwrap()
        );

        let pysecdec_normalization_correction = vk_parse!(
            format!(
                "(  𝑖*(𝜋^((4-2*{eps})/2))\
                )^{n_loops}",
                eps = settings.epsilon_symbol,
                n_loops = n_loops
            )
            .as_str()
        )
        .unwrap();

        let evaluated_integral = output
            .iter()
            .enumerate()
            .map(
                |(stream, out)| match vk_parse!(out.replace("I", "𝑖").as_str()) {
                    Ok(mut processed) => {
                        processed = processed
                            .replace(vk_parse!("ep").unwrap().to_pattern())
                            .with(vk_parse!(&settings.epsilon_symbol).unwrap().to_pattern());
                        processed = precision::legacy::stage_table(
                            function!(
                                vk_symbol!("pysecdec_result"),
                                source_id as i64,
                                stream as i64
                            )
                            .as_view(),
                            processed.as_view(),
                            vk_symbol!(&settings.epsilon_symbol),
                            "PySecDec",
                        )?;
                        processed = processed
                            * pysecdec_normalization_correction.to_owned()
                            * settings
                                .get_staged_integral_normalization_factor_atom()?
                                .replace(S.n_loops.to_pattern())
                                .with(Atom::num(n_loops).to_pattern());

                        let expanded_evaluation = match processed.series(
                            vk_symbol!(settings.epsilon_symbol.as_str()),
                            Atom::Zero.as_atom_view(),
                            Rational::from(
                                settings.number_of_terms_in_epsilon_expansion - n_loops - 1,
                            ),
                        ) {
                            Ok(a) => a,
                            Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
                        };

                        let mut evaluated_i = expanded_evaluation.to_atom();

                        evaluated_i = evaluated_i
                            .replace(vk_parse!("log_mu_sq").unwrap().to_pattern())
                            .with(log_mu_sq.to_pattern());

                        Ok(evaluated_i)
                    }
                    Err(err) => Err(VakintError::PySecDecOutputParsingError(out.clone(), err)),
                },
            )
            .collect::<Result<Vec<_>, VakintError>>()?;

        Ok(evaluated_integral[0].to_owned()
            + evaluated_integral[1].to_owned() * S.error_flag.to_owned())
    }

    // Identify from the short canonical expression of the integral what is the UV mass and returns
    // the atom for it and its square (one of the two will effectively be a symbol but the other one an expression).
    // Throws an error if the identification fails.
    fn identify_uv_mass_symbols(
        short_integral_expression: &Atom,
    ) -> Result<(Atom, Atom), VakintError> {
        let (muv_atom, muv_sq_atom) = if let Some(m) = short_integral_expression
            .pattern_match(
                &vk_parse!("prop(args__,m_sq_,pow_)").unwrap().to_pattern(),
                None,
                None,
            )
            .next()
        {
            match m.get(&vk_symbol!("m_sq_")).unwrap() {
                Atom::Var(s) => (
                    Atom::pow(&Atom::Var(s.clone()), Atom::num(1) / Atom::num(2)),
                    Atom::Var(s.clone()),
                ),
                Atom::Pow(m) => {
                    let base = m.to_pow_view().get_base().to_owned();
                    let exp = m.to_pow_view().get_exp();
                    match base {
                        Atom::Var(s) if exp == Atom::num(2).as_view() => {
                            (Atom::Var(s.clone()), Atom::Pow(m.clone()))
                        }

                        _ => {
                            return Err(VakintError::MalformedGraph(format!(
                                "Could not find muV in graph:\n{}",
                                short_integral_expression.to_canonical_string()
                            )));
                        }
                    }
                }
                _ => {
                    return Err(VakintError::MalformedGraph(format!(
                        "Could not find muV in graph:\n{}",
                        short_integral_expression.to_canonical_string()
                    )));
                }
            }
        } else {
            return Err(VakintError::MalformedGraph(format!(
                "Could not find muV in graph:\n{}",
                short_integral_expression.to_canonical_string()
            )));
        };

        Ok((muv_atom, muv_sq_atom))
    }

    fn alphaloop_evaluate(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
        opts: &AlphaLoopOptions,
    ) -> Result<Atom, VakintError> {
        let integral = integral_specs.canonical_topology.get_integral();

        debug!(
            "Processing the following integral with {}:\n{}",
            "AlphaLoop".green(),
            integral
        );
        // println!(
        //     "VH:: B:: Processing the following integral with {}:\n{}\n and numerator:\n{}",
        //     "AlphaLoop".green(),
        //     integral,
        //     numerator
        // );

        // DO NOT REQUIRE MASS TO BE A SYMBOL
        let (muv_atom, muv_sq_atom) =
            Vakint::identify_uv_mass_symbols(integral.canonical_expression.as_ref().unwrap())?;

        let alphaloop_expression = integral.alphaloop_expression.as_ref().unwrap().as_view();

        let mut numerator = numerator.to_owned();
        numerator = numerator.replace_multiple(&[
            Replacement::new(
                muv_sq_atom.to_pattern(),
                vk_parse!("mUV^2").unwrap().to_pattern(),
            ),
            Replacement::new(
                muv_atom.to_pattern(),
                vk_parse!("mUV").unwrap().to_pattern(),
            ),
        ]);
        // println!("Numerator : {}", numerator);
        // println!("Evaluating AlphaLoop : {}", alphaloop_expression);
        // println!("Graph:\n{}", integral.graph.to_graphviz());

        // Make sure to undo the dot product notation.
        // If it was not used, the command below will do nothing.
        let mut form_expression = numerator * alphaloop_expression.to_owned();
        form_expression = Vakint::convert_to_dot_notation(form_expression.as_view());
        // println!("Input expression with dot products : {}", form_expression);

        let mut vector_mapping: BTreeMap<Atom, Atom> = BTreeMap::new();
        let vec_pattern = vk_parse!("v_(id_)").unwrap().to_pattern();
        let vec_conditions = Condition::from((vk_symbol!("id_"), number_condition()))
            & Condition::from((vk_symbol!("v_"), symbol_condition()));
        let form_expression_clone = form_expression.clone();

        let matcher =
            form_expression_clone.pattern_match(&vec_pattern, Some(&vec_conditions), None);
        for m in matcher {
            let v = match m.get(&vk_symbol!("v_")).unwrap() {
                Atom::Var(s) => s.get_symbol(),
                _ => {
                    return Err(VakintError::MalformedGraph(format!(
                        "Could not find v in graph:\n{}",
                        integral.canonical_expression.as_ref().unwrap()
                    )));
                }
            };
            if v != S.p && v != S.k {
                continue;
            }
            let v_id = get_integer_from_atom(m.get(&vk_symbol!("id_")).unwrap().as_view()).unwrap();

            form_expression = form_expression
                .replace(
                    vk_parse!(format!("{}({})", v, v_id).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .with(
                    vk_parse!(format!("{}{}", v, v_id).as_str())
                        .unwrap()
                        .to_pattern(),
                );
            vector_mapping.insert(
                vk_parse!(format!("{}{}", v, v_id).as_str()).unwrap(),
                vk_parse!(format!("{}({})", v, v_id).as_str()).unwrap(),
            );
        }

        let vec_pattern = vk_parse!("v_(id_,idx_)").unwrap().to_pattern();
        let form_expression_clone = form_expression.clone();
        let matcher =
            form_expression_clone.pattern_match(&vec_pattern, Some(&vec_conditions), None);
        for m in matcher {
            let v = match m.get(&vk_symbol!("v_")).unwrap() {
                Atom::Var(s) => s.get_symbol(),
                _ => {
                    return Err(VakintError::MalformedGraph(format!(
                        "Could not find v in graph:\n{}",
                        integral.canonical_expression.as_ref().unwrap()
                    )));
                }
            };
            if v != S.p && v != S.k {
                continue;
            }
            let v_id = get_integer_from_atom(m.get(&vk_symbol!("id_")).unwrap().as_view()).unwrap();
            let idx = m.get(&vk_symbol!("idx_")).unwrap();
            let vector_symbol = vk_symbol!(format!("{}{}", v, v_id));
            let vector = function!(v, v_id, idx);
            form_expression = form_expression.replace(vector.to_pattern()).with(function!(
                vk_symbol!("vec1"),
                Atom::var(vector_symbol),
                idx
            ));
            vector_mapping.insert(function!(vector_symbol, idx), vector);
        }
        // println!("Input expression for FORM : {}", form_expression);

        let template = Template::parse_template(
            TEMPLATES
                .get("run_alphaloop_integral_evaluation.txt")
                .unwrap(),
        )
        .unwrap();

        let mut vars: HashMap<String, String> = HashMap::new();

        let (form_header_additions, form_expression, indices) =
            self.prepare_expression_for_form(settings, form_expression, true, &[])?;

        vars.insert("numerator".into(), form_expression);
        vars.insert("additional_symbols".into(), form_header_additions);

        let rendered = template
            .render(&RenderOptions {
                variables: vars,
                ..Default::default()
            })
            .unwrap();
        let form_result = self.run_form(
            settings,
            &["integrateduv.frm".into()],
            ("run_alphaloop_integral_evaluation.frm".into(), rendered),
            vec![
                "-D".into(),
                format!("MAXPOLE={}", integral.n_loops),
                "-D".into(),
                format!(
                    "SELECTEDEPSILONORDER={}",
                    settings.number_of_terms_in_epsilon_expansion - (integral.n_loops as i64) - 1
                ),
            ],
            settings.clean_tmp_dir,
            settings.temporary_directory.clone(),
        )?;

        // println!(
        //     "\n>>>> VH:: B:: A Indices:\n{}",
        //     indices
        //         .clone()
        //         .iter()
        //         .map(|idx| idx.to_canonical_string())
        //         .collect::<Vec<_>>()
        //         .join(",")
        // );
        // println!(
        //     "\n>>>> VH:: B:: Vector mapping:\n{}",
        //     vector_mapping
        //         .iter()
        //         .map(|(k, v)| format!("{} -> {}", k, v.to_canonical_string()))
        //         .collect::<Vec<_>>()
        //         .join("\n")
        // );

        let mut evaluated_integral =
            self.process_form_output(settings, form_result, indices, vector_mapping)?;
        debug!(
            "{}: raw result from FORM:\n{}",
            "AlphaLoop".green(),
            evaluated_integral
        );

        if !settings.use_dot_product_notation {
            evaluated_integral = Vakint::convert_from_dot_notation(evaluated_integral.as_view());
        }

        evaluated_integral = evaluated_integral
            .replace(vk_parse!("mUV").unwrap().to_pattern())
            .with(muv_atom.to_pattern());
        let log_muv_mu_sq = function!(
            Symbol::LOG,
            muv_sq_atom / Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );

        let log_mu_sq = function!(
            Symbol::LOG,
            Atom::var(vk_symbol!(settings.mu_r_sq_symbol.as_str()))
        );
        if opts.susbstitute_masters {
            let processed_constants = DIRECT_SUBSTITUTIONS
                .iter()
                .map(|(src, (trgt, condition))| {
                    let target =
                        if inspection::constants::builtin_definition(src.as_view())?.is_some() {
                            src.clone()
                        } else {
                            precision::legacy::stage_table(
                                src.as_view(),
                                trgt.as_view(),
                                vk_symbol!("ep"),
                                "AlphaLoop",
                            )?
                        };
                    Ok((src, (target, condition.clone())))
                })
                .collect::<Result<Vec<_>, VakintError>>()?;
            // for (a, (b, _c)) in processed_constants.clone() {
            //     println!("{} -> {}", a, b);
            // }

            let mut r = evaluated_integral.to_owned();
            r.repeat_map(Box::new(move |av: AtomView| {
                let mut res = av.to_owned();
                for (src, (trgt, matching_condition)) in processed_constants.iter() {
                    res = res
                        .replace(src.to_pattern())
                        .when(matching_condition)
                        .with(trgt.to_pattern());
                }
                res
            }));
            evaluated_integral = crate::hpl::normalize_constants(r.as_view())?;
        }

        //*((2*𝜋)^(2*{eps}))\
        //*((4*𝜋*exp(EulerGamma)))^{eps}
        //*((2*𝜋)^(4-2*{eps}))
        //*𝑖*(𝜋^((4-2*{eps})/2))
        // Make sure to split off the logarithmic terms with one term showing explicitely a ratio of scales and the other
        // having just a logarithm of the renormalization scale so that cancellations are symbolic when using `log_mu_sq`
        // in the normalization choice.
        // We must keep the name logmUVmu as it is reserved in the alphaloop implementation and corresponds to log(mUV^2/mu^2)
        // This is also the reason we do not simplify the expression exp(-logmUVmu+log_mu_sq)
        let alphaloop_normalization_correction = vk_parse!(
            format!(
                "(\
                    𝑖*(𝜋^((4-2*{eps})/2))\
                 * (exp(-EulerGamma))^({eps})\
                 * (exp(-logmUVmu-log_mu_sq))^({eps})\
                 )^{n_loops}",
                eps = settings.epsilon_symbol,
                n_loops = integral.n_loops
            )
            .as_str()
        )
        .unwrap();

        evaluated_integral = evaluated_integral
            * alphaloop_normalization_correction
            * settings
                .get_staged_integral_normalization_factor_atom()?
                .as_view()
                .replace(S.n_loops.to_pattern())
                .with(Atom::num(integral.n_loops as i64).to_pattern());

        let expanded_evaluation = match evaluated_integral.series(
            vk_symbol!(settings.epsilon_symbol.as_str()),
            Atom::Zero.as_atom_view(),
            Rational::from(
                settings.number_of_terms_in_epsilon_expansion - (integral.n_loops as i64) - 1,
            ),
        ) {
            Ok(a) => a,
            Err(e) => return Err(VakintError::SymbolicaError(e.to_string())),
        };
        evaluated_integral = expanded_evaluation.to_atom();

        evaluated_integral = evaluated_integral
            .replace(vk_parse!("logmUVmu").unwrap().to_pattern())
            .with(log_muv_mu_sq.to_pattern());
        evaluated_integral = evaluated_integral
            .replace(vk_parse!("log_mu_sq").unwrap().to_pattern())
            .with(log_mu_sq.to_pattern());

        Ok(evaluated_integral)
    }

    fn get_pysecdec_version(&self, settings: &VakintSettings) -> Result<String, VakintError> {
        let mut cmd = Command::new(settings.python_exe_path.as_str());
        cmd.arg("-c");
        cmd.arg("import pySecDec; print(pySecDec.__version__)");
        let output = if let Ok(o) = cmd.output() {
            o
        } else {
            return Err(VakintError::PySecDecUnavailable);
        };
        if !ExitStatus::success(&output.status) {
            return Err(VakintError::PySecDecUnavailable);
        }
        let output_str = String::from_utf8_lossy(&output.stdout).into_owned();
        let re = Regex::new(r"([\.|\d]+)").unwrap();
        let mut versions = vec![];
        for (_, [version]) in re.captures_iter(output_str.as_str()).map(|ci| ci.extract()) {
            versions.push(version);
        }
        if versions.is_empty() {
            return Err(VakintError::PySecDecVersion(format!(
                "Could not obtain PySecDec version from command:\n{:?}\nOutput was:\n{}",
                cmd, output_str
            )));
        }
        Ok(versions[0].into())
    }

    fn get_form_version(&self, settings: &VakintSettings) -> Result<String, VakintError> {
        let mut cmd = Command::new(settings.form_exe_path.as_str());
        cmd.arg("-version");
        let output = if let Ok(o) = cmd.output() {
            o
        } else {
            return Err(VakintError::FormUnavailable);
        };

        if !ExitStatus::success(&output.status) {
            return Err(VakintError::FormUnavailable);
        }
        let output_str = String::from_utf8_lossy(&output.stdout).into_owned();
        let re = Regex::new(r"FORM\s([\.|\d]+)").unwrap();
        let mut versions = vec![];
        for (_, [version]) in re.captures_iter(output_str.as_str()).map(|ci| ci.extract()) {
            versions.push(version);
        }
        if versions.is_empty() {
            return Err(VakintError::FormVersion(format!(
                "Could not obtain form version from command:\n{:?}\nOutput was:\n{}",
                cmd, output_str
            )));
        }
        Ok(versions[0].into())
    }

    pub fn convert_from_dot_notation(atom: AtomView) -> Atom {
        let dummy_pat = function!(S.dot_dummy_ind, S.a_).to_pattern();

        atom.replace(S.dot(S.a_, S.b_).pow(Atom::var(S.c_)))
            .with(S.dot_pow(S.a_, S.b_, S.c_))
            .replace(S.dot(S.a_, S.b_))
            .with(S.dot_pow(S.a_, S.b_, 1))
            .replace(S.dot_pow(S.a_, S.b_, S.c_))
            .when(S.c_.filter(|a| a < 0))
            .with(S.dot_pow(S.a_, S.b_, -Atom::var(S.c_)).pow(-1))
            .map_terms_single_core(|a| {
                let mut dummy = 1;
                for m in a.replace(&dummy_pat).match_iter() {
                    let Ok(d) = usize::try_from(&m[&S.a_]) else {
                        continue;
                    };
                    if dummy <= d {
                        dummy = d;
                    }
                }
                dummy += 1;
                let first_new_dummy = dummy;

                let mut a = a.to_owned();
                loop {
                    let dummy_wrapped = S.dot_dummy_ind(dummy);
                    let new = a
                        .replace(S.dot_pow(S.a_, S.b_, S.c_))
                        .once()
                        .with_map(move |m| {
                            let Match::Single(a) = m.get(S.a_).unwrap() else {
                                unreachable!()
                            };
                            let Match::Single(b) = m.get(S.b_).unwrap() else {
                                unreachable!()
                            };
                            let rest = S
                                .dot_pow(S.a_, S.b_, Atom::var(S.c_) - 1)
                                .to_pattern()
                                .replace_wildcards_with_matches(m);

                            // if let (Ok(a), Ok(b)) =
                            //     (FunctionBuilder::try_from(m.get(S.a_)), FunctionBuilder::try_from(b))
                            // {
                            //     a.add_arg(dummy).finish() * b.add_arg(dummy).finish() * rest
                            // }
                            match (a, b) {
                                (AtomView::Var(a), AtomView::Var(b)) => {
                                    FunctionBuilder::new(a.get_symbol())
                                        .add_arg(&dummy_wrapped)
                                        .finish()
                                        * FunctionBuilder::new(b.get_symbol())
                                            .add_arg(&dummy_wrapped)
                                            .finish()
                                        * rest
                                }
                                (AtomView::Var(a), AtomView::Fun(b)) => {
                                    FunctionBuilder::new(a.get_symbol())
                                        .add_arg(&dummy_wrapped)
                                        .finish()
                                        * FunctionBuilder::new(b.get_symbol())
                                            .add_args(b.iter().collect::<Vec<_>>())
                                            .add_arg(&dummy_wrapped)
                                            .finish()
                                        * rest
                                }
                                (AtomView::Fun(a), AtomView::Var(b)) => {
                                    FunctionBuilder::new(a.get_symbol())
                                        .add_args(a.iter().collect::<Vec<_>>())
                                        .add_arg(&dummy_wrapped)
                                        .finish()
                                        * FunctionBuilder::new(b.get_symbol())
                                            .add_arg(&dummy_wrapped)
                                            .finish()
                                        * rest
                                }
                                (AtomView::Fun(a), AtomView::Fun(b)) => {
                                    FunctionBuilder::new(a.get_symbol())
                                        .add_args(a.iter().collect::<Vec<_>>())
                                        .add_arg(&dummy_wrapped)
                                        .finish()
                                        * FunctionBuilder::new(b.get_symbol())
                                            .add_args(b.iter().collect::<Vec<_>>())
                                            .add_arg(&dummy_wrapped)
                                            .finish()
                                        * rest
                                }
                                _ => S
                                    .dot_pow(S.a_, S.b_, S.c_)
                                    .to_pattern()
                                    .replace_wildcards_with_matches(m),
                            }
                        })
                        .replace(S.dot_pow(S.a_, S.b_, 0))
                        .with(1);

                    if new == a {
                        break;
                    }
                    a = new;
                    dummy += 1;
                }
                // Each copy of a powered scalar sum owns independent dummy
                // contractions. Keep those copies factorized: expanding the
                // numerator before conversion would hide index capture.
                // FORM expands these powers; reused indices would let it contract
                // loop vectors together instead of with their external partners.
                a.replace_map_bottom_up(|view, _, out| {
                    let AtomView::Pow(power) = view else {
                        return;
                    };
                    let (base, exponent) = power.get_base_exp();
                    let Ok(exponent) = i64::try_from(exponent) else {
                        return;
                    };
                    let copies = exponent.unsigned_abs();
                    if copies < 2 || !matches!(base, AtomView::Add(_) | AtomView::Mul(_)) {
                        return;
                    }
                    let mut indices = base
                        .replace(&dummy_pat)
                        .match_iter()
                        .filter_map(|m| usize::try_from(&m[&S.a_]).ok())
                        // Existing indexed tensors keep their contractions;
                        // only dots converted in this call introduce copies.
                        .filter(|index| *index >= first_new_dummy)
                        .collect::<Vec<_>>();
                    indices.sort_unstable();
                    indices.dedup();
                    if indices.is_empty() {
                        return;
                    }
                    let mut product = Atom::one();
                    for _ in 0..copies {
                        let replacements = indices.iter().map(|index| {
                            let replacement =
                                Replacement::new(S.dot_dummy_ind(*index), S.dot_dummy_ind(dummy));
                            dummy += 1;
                            replacement
                        });
                        product *= base.replace_multiple(replacements);
                    }
                    **out = product.pow(exponent.signum());
                })
            })
    }

    pub fn convert_to_dot_notation(atom: AtomView) -> Atom {
        let mut old_expr = atom.to_owned().expand();
        // Native CAS functions are scalar even when they have two integer
        // arguments: polygamma(1, 3)^2 is not a squared indexed vector.
        let vector_head_condition = WildcardRestriction::Filter(Box::new(|matched| {
            let symbol = match matched {
                Match::Single(AtomView::Var(variable)) => variable.get_symbol(),
                Match::FunctionName(symbol) => *symbol,
                _ => return false,
            };
            !crate::hpl::is_projection_symbol(symbol)
                && (symbol.get_namespace() != "symbolica"
                    || !symbolica::state::State::is_builtin(symbol.get_stripped_name()))
        }));

        loop {
            let mut expr = old_expr
                .replace(vk_parse!("v_(id_,idx_)^n_").unwrap().to_pattern())
                .when(
                    &(Condition::from((vk_symbol!("v_"), vector_head_condition.clone()))
                        & Condition::from((vk_symbol!("id_"), number_condition()))
                        & Condition::from((vk_symbol!("n_"), even_condition()))),
                )
                .with(
                    vk_parse!("dot(v_(id_),v_(id_))^(n_/2)")
                        .unwrap()
                        .to_pattern(),
                );

            // dot products
            expr = expr
                .replace(
                    vk_parse!("v1_(id1_,idx_)*v2_(id2_,idx_)")
                        .unwrap()
                        .to_pattern(),
                )
                .when(
                    &(Condition::from((vk_symbol!("v1_"), vector_head_condition.clone()))
                        & Condition::from((vk_symbol!("v2_"), vector_head_condition.clone()))
                        & Condition::from((vk_symbol!("id1_"), number_condition()))
                        & Condition::from((vk_symbol!("id2_"), number_condition()))),
                )
                .with(vk_parse!("dot(v1_(id1_),v2_(id2_))").unwrap().to_pattern());

            // metric contraction
            expr = expr
                .replace(
                    vk_parse!("g(idx1_,idx2_)*v1_(id1_,idx1_)*v2_(id2_,idx2_)")
                        .unwrap()
                        .to_pattern(),
                )
                .when(
                    &(Condition::from((vk_symbol!("v1_"), vector_head_condition.clone()))
                        & Condition::from((vk_symbol!("v2_"), vector_head_condition.clone()))
                        & Condition::from((vk_symbol!("id1_"), number_condition()))
                        & Condition::from((vk_symbol!("id2_"), number_condition()))),
                )
                .with(vk_parse!("dot(v1_(id1_),v2_(id2_))").unwrap().to_pattern());

            if expr == old_expr {
                break;
            } else {
                old_expr = expr.to_owned();
            }
        }
        old_expr
    }

    pub fn identify_vector_indices(numerator: AtomView) -> Result<Vec<Atom>, VakintError> {
        let mut indices = HashSet::new();
        // make sure the numerator is in the form of vec_(id_,idx_)
        let vector_matcher_pattern = vk_parse!("vec_(id_,idx_)").unwrap().to_pattern();
        let vector_conditions = Condition::from((vk_symbol!("vec_"), symbol_condition()))
            & Condition::from((vk_symbol!("id_"), symbol_condition()));
        let vector_match_settings = MatchSettings::default();
        let mut vector_matcher = numerator.pattern_match(
            &vector_matcher_pattern,
            Some(&vector_conditions),
            Some(&vector_match_settings),
        );

        while let Some(m) = vector_matcher.next_detailed() {
            if let Match::FunctionName(vec_symbol) = m.match_stack.get(vk_symbol!("vec_")).unwrap()
            {
                if [
                    // LOOP_MOMENTUM_SYMBOL,
                    // EXTERNAL_MOMENTUM_SYMBOL,
                    "vec1", "vec",
                ]
                .iter()
                .any(|s| *vec_symbol == vk_symbol!(s))
                {
                    indices.insert(m.match_stack.get(vk_symbol!("idx_")).unwrap().to_atom());
                }
            } else {
                unreachable!("Vector name should be a symbol.")
            }
        }
        for metric_symbol in ["g", METRIC_SYMBOL] {
            let metric_pattern = vk_parse!(format!("{}(idx1_,idx2_)", metric_symbol))
                .unwrap()
                .to_pattern();
            let mut metric_matcher = numerator.pattern_match(&metric_pattern, None, None);
            while let Some(m) = metric_matcher.next_detailed() {
                indices.insert(m.match_stack.get(vk_symbol!("idx1_")).unwrap().to_atom());
                indices.insert(m.match_stack.get(vk_symbol!("idx2_")).unwrap().to_atom());
            }
        }
        Ok(indices.iter().cloned().collect::<Vec<_>>())
    }

    fn serialize_expression(expression: AtomView) -> String {
        // Backend input follows the exact AST, including its existing argument
        // order and antisymmetric signs. Display callbacks must never supply
        // algebra or symbol identities, and complex numbers need parentheses.
        match expression {
            AtomView::Num(number) => {
                if number.get_coeff_view().is_real() {
                    return expression.to_canonical_string();
                }
                // Serialize the two components explicitly. Display printing
                // may omit a unit imaginary numerator, while backend adapters
                // need an explicit numeric token for their imaginary unit.
                let (real, imaginary) = match number.get_coeff_view().to_owned() {
                    Coefficient::Complex(value) => (Atom::num(value.re), Atom::num(value.im)),
                    Coefficient::Float(value) => (Atom::num(value.re), Atom::num(value.im)),
                    _ => return format!("({})", expression.to_canonical_string()),
                };
                format!(
                    "({}+({})*1𝑖)",
                    real.to_canonical_string(),
                    imaginary.to_canonical_string()
                )
            }
            AtomView::Var(variable) => get_full_name(&variable.get_symbol()),
            AtomView::Fun(function) => format!(
                "{}({})",
                get_full_name(&function.get_symbol()),
                function
                    .iter()
                    .map(Self::serialize_expression)
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                format!(
                    "({})^({})",
                    Self::serialize_expression(base),
                    Self::serialize_expression(exponent)
                )
            }
            AtomView::Mul(product) => format!(
                "({})",
                product
                    .iter()
                    .map(Self::serialize_expression)
                    .collect::<Vec<_>>()
                    .join("*")
            ),
            AtomView::Add(sum) => format!(
                "({})",
                sum.iter()
                    .map(Self::serialize_expression)
                    .collect::<Vec<_>>()
                    .join("+")
            ),
        }
    }

    pub fn sanitize_user_expressions(
        &self,
        settings: &VakintSettings,
        expression: AtomView,
        substitute_indices: bool,
        additional_user_symbols: &[Symbol],
    ) -> Result<(String, String, Vec<Atom>), VakintError> {
        let mut expression = expression.to_owned();

        // If there are floats in the expression we must rationalize them to target precision
        // because the FORM version we want to support here does not necessarily support floats.
        // First wrap all coefficients within a marker function so that we can floatify them back after
        // FORM has processed them

        let binary_prec = settings.get_binary_precision();
        expression = expression.replace_map(|term, _ctx, out| {
            if let AtomView::Num(c) = term
                && let CoefficientView::Float(re, im) = c.get_coeff_view()
            {
                let mut re = re.to_float();
                let mut im = im.to_float();
                match settings.precision_for_input_float_rationalization {
                    InputFloatRationalizationPrecision::FullPrecision => {}
                    InputFloatRationalizationPrecision::TargetPrecision => {
                        re.set_prec(binary_prec.min(re.prec()));
                        im.set_prec(binary_prec.min(im.prec()));
                    }
                };
                **out = function!(
                    S.float_marker,
                    re.to_rational(),
                    im.to_rational(),
                    re.get_precision(),
                    im.get_precision()
                );
            }
        });

        // Preserve namespaces, symbol attributes and complex numeric precedence.
        let mut processed_str = Self::serialize_expression(expression.as_view());
        // println!("Original expression: {}", expression.to_canonical_string());
        // Identify user indices in p and k structures
        let mut indices = Vakint::identify_vector_indices(expression.as_view())?;
        // println!(
        //     "Indices: {}",
        //     indices
        //         .iter()
        //         .map(|i| i.to_canonical_string())
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );
        let mut form_header_indices = vec![];
        // Temporarily replace the indices to not have them appear in the symbols list
        let mut expression_no_indices = expression.to_owned();

        for (i_index, user_i) in indices.iter().enumerate() {
            if substitute_indices {
                // Reuse index atoms directly: canonical text omits symbol user data.
                let integer_index =
                    Atom::num((FORM_REPLACEMENT_INDEX_SHIFT + 1 + i_index as u64) as i64);
                for vecsymbol in [
                    vk_symbol!(LOOP_MOMENTUM_SYMBOL),
                    vk_symbol!(EXTERNAL_MOMENTUM_SYMBOL),
                    vk_symbol!("vec1"),
                    vk_symbol!("vec"),
                ] {
                    let id = Atom::var(vk_symbol!("id_"));
                    expression_no_indices = expression_no_indices
                        .replace(function!(vecsymbol, &id, user_i).to_pattern())
                        .with(function!(vecsymbol, &id, &integer_index));
                }
                for metric_symbol in ["g", METRIC_SYMBOL] {
                    let metric = vk_symbol!(metric_symbol);
                    let idx1 = Atom::var(vk_symbol!("idx1_"));
                    let idx2 = Atom::var(vk_symbol!("idx2_"));
                    expression_no_indices = expression_no_indices
                        .replace(function!(metric, &idx1, user_i).to_pattern())
                        .with(function!(metric, &idx1, &integer_index))
                        .replace(function!(metric, user_i, &idx2).to_pattern())
                        .with(function!(metric, &integer_index, &idx2));
                }
            } else {
                let litteral_form_name =
                    format!("[{}]", Vakint::serialize_expression(user_i.as_view()));
                expression_no_indices = expression_no_indices
                    .replace(user_i.to_pattern())
                    .with(Atom::num(0));
                processed_str = processed_str.replace(
                    &Vakint::serialize_expression(user_i.as_view()),
                    &litteral_form_name,
                );
                form_header_indices.push(litteral_form_name);
            }
        }
        if substitute_indices {
            processed_str = Self::serialize_expression(expression_no_indices.as_view());
        }
        processed_str = undress_vakint_symbols(&processed_str);

        // println!("processed_str= {}", processed_str);

        // // let tt = symbolica::parse!("vakint::{}::uvprop(vakint::{}::k1,1)*vakint::{}::vxs(-1*vakint::{}::k1,vakint::{}::k1)");
        // // let symbols = tt.get_all_symbols(true);
        // // println!("Symbols found: {}", symbols.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", "));
        // // let functions = tt.get_all_symbols(false);
        // // println!("Functions found: {}", functions.iter().map(|s| s.to_string()).collect::<Vec<_>>().join(", "));

        let user_symbols = expression_no_indices.get_all_symbols(true);
        let mut user_variables = expression_no_indices.get_all_symbols(false);
        let mut user_functions: HashSet<Symbol, RandomState> = user_symbols
            .iter()
            .filter(|s| !user_variables.contains(s))
            .cloned()
            .collect::<HashSet<_, _>>();

        user_variables = user_variables
            .iter()
            //.filter(|s| !SYMBOL_REGISTRY.contains(s))
            .filter(|s| S.should_symbol_be_escaped_in_form(s))
            //            .filter(|s| s.get_namespace() != NAMESPACE || !s.get_attributes().is_empty() )
            .cloned()
            .collect::<HashSet<_, _>>();
        user_functions = user_functions
            .iter()
            //.filter(|s| !SYMBOL_REGISTRY.contains(s))
            .filter(|s: &&Symbol| S.should_symbol_be_escaped_in_form(s))
            //            .filter(|s: &&Symbol| s.get_namespace() != NAMESPACE || !s.get_attributes().is_empty() )
            .cloned()
            .collect::<HashSet<_, _>>();
        // println!(
        //     "User variables: {}",
        //     user_variables
        //         .iter()
        //         .map(|s| s.to_string())
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );
        // println!(
        //     "User functions: {}",
        //     user_functions
        //         .iter()
        //         .map(|s| s.to_string())
        //         .collect::<Vec<_>>()
        //         .join(", ")
        // );

        let mut form_header_functions = vec![];
        let mut form_header_symbols = vec![];

        user_variables.insert(vk_symbol!(settings.mu_r_sq_symbol.clone()));
        for s in additional_user_symbols {
            user_variables.insert(*s);
        }
        let mut string_replacements: HashMap<String, String, ahash::RandomState> =
            HashMap::default();
        // FORM returns names to this registry; reuse symbols instead of redeclaring attributes
        // without their user data, normalization callbacks, and other registered metadata.
        for user_f in user_functions.iter() {
            let litteral_form_name = format!("[{}]", user_f.get_name());
            if user_f.get_namespace() == NAMESPACE || user_f.get_namespace() == "symbolica" {
                //            if (user_f.get_namespace() == NAMESPACE || user_f.get_namespace() == "symbolica") && user_f.get_attributes().is_empty() {
                string_replacements.insert(
                    user_f.get_stripped_name().into(),
                    litteral_form_name.clone(),
                );
                // processed_str = processed_str.replace(user_f.get_stripped_name(), &litteral_form_name);
            } else {
                string_replacements.insert(get_full_name(user_f), litteral_form_name.clone());
            }
            form_header_functions.push(litteral_form_name);
        }
        for user_v in user_variables.iter() {
            let litteral_form_name = format!("[{}]", user_v.get_name());
            if user_v.get_namespace() == NAMESPACE || user_v.get_namespace() == "symbolica" {
                //            if (user_v.get_namespace() == NAMESPACE || user_v.get_namespace() == "symbolica") && user_v.get_attributes().is_empty() {
                string_replacements.insert(
                    user_v.get_stripped_name().into(),
                    litteral_form_name.clone(),
                );
                // processed_str = processed_str.replace(user_v.get_stripped_name(), &litteral_form_name);
            } else {
                string_replacements.insert(get_full_name(user_v), litteral_form_name.clone());
            }
            form_header_symbols.push(litteral_form_name);
        }
        let mut form_header_additions = vec![];
        if !form_header_functions.is_empty() {
            form_header_additions.push(format!("CF {};", form_header_functions.join(", ")));
        }
        if !form_header_indices.is_empty() {
            form_header_additions.push(format!("I {};", form_header_indices.join(", ")));
        }
        if !form_header_symbols.is_empty() {
            form_header_additions.push(format!("S {};", form_header_symbols.join(", ")));
        }

        // If we did not substitute indices with integers, we won't have to map them back,
        // so we will return an empty list of indices here.
        if !substitute_indices {
            indices.clear();
        }
        // Numeric serialization emits an explicit coefficient before the unit.
        string_replacements.insert("𝑖".into(), "*i_".into());
        // println!("String replacements: {}", string_replacements.iter().map(|(k,v)| format!("{} -> {}", k, v)).collect::<Vec<_>>().join(", "));

        processed_str = utils::multi_string_replace(processed_str, &string_replacements);
        //println!("form_header_additions: {}", form_header_additions.join("\n"));

        Ok((form_header_additions.join("\n"), processed_str, indices))
    }

    pub fn prepare_expression_for_form(
        &self,
        settings: &VakintSettings,
        expression: Atom,
        substitute_indices: bool,
        additional_user_symbols: &[Symbol],
    ) -> Result<(String, String, Vec<Atom>), VakintError> {
        let mut processed = expression.clone();
        processed = processed
            .replace(vk_parse!(&settings.epsilon_symbol).unwrap().to_pattern())
            .with(vk_parse!("ep").unwrap().to_pattern());

        let (form_header_additions, expression_str, indices) = self.sanitize_user_expressions(
            settings,
            processed.as_view(),
            substitute_indices,
            additional_user_symbols,
        )?;

        // Index metadata bypasses FORM and also identifies vector-mapping
        // sources. Return it in the caller's epsilon convention, as at input.
        let indices = indices
            .into_iter()
            .map(|index| {
                index
                    .replace(vk_parse!("ep").unwrap().to_pattern())
                    .with(vk_parse!(&settings.epsilon_symbol).unwrap())
            })
            .collect();

        Ok((form_header_additions, expression_str, indices))
    }

    fn parse_form_output_expression(expression: &str) -> Result<Atom, String> {
        use symbolica::{
            atom::UserData,
            parser::{ParseSettings, Token},
            state::{State, Workspace},
            with_default_namespace,
        };

        // Literal indices carry canonical declarations, which omit user data.
        // Reuse exact registered declarations, retaining representation metadata
        // and the parser's rejection of conflicting attributes or tags.
        let mut symbols = State::symbol_iter()
            .filter(|(symbol, _)| !matches!(symbol.get_data(), UserData::None))
            .map(|(symbol, _)| (Atom::var(symbol).to_canonical_string().into(), symbol))
            .collect();
        let namespace = with_default_namespace!(expression, NAMESPACE);
        Workspace::get_local().with(|workspace| {
            Token::parse_with_atom_info(
                expression,
                ParseSettings::default(),
                Some((&namespace, &mut symbols, workspace)),
            )?
            .to_atom(&namespace, &mut symbols, workspace)
        })
    }

    pub fn process_form_output(
        &self,
        settings: &VakintSettings,
        form_output: String,
        indices: Vec<Atom>,
        vector_mapping: BTreeMap<Atom, Atom>,
    ) -> Result<Atom, VakintError> {
        let processed_form_str = form_output
            .replace("[", "")
            .replace("]", "")
            .replace("i_", "1𝑖")
            .replace("\\\n", "\n")
            .split("\n")
            .map(|s| s.trim())
            .collect::<Vec<_>>()
            .join("");

        // Make sure to replace the indices already in the source of the vector mappings since this replacement is done *before* substituting the indices
        let vector_mapping: BTreeMap<Atom, Atom> = vector_mapping
            .iter()
            .map(|(s, t)| {
                let mut new_s = s.to_owned();
                for (i_index, user_i) in indices.iter().enumerate() {
                    let integer_index =
                        Atom::num((FORM_REPLACEMENT_INDEX_SHIFT + 1 + (i_index as u64)) as i64);
                    new_s = new_s
                        .replace(user_i.as_view().to_pattern())
                        .with(integer_index);
                }
                (new_s, t.clone())
            })
            .collect();
        // println!(
        //     "\n>>>> VH:: C:: Processed vector mappings:\n{}",
        //     vector_mapping
        //         .iter()
        //         .map(|(s, t)| format!("{} -> {}", s, t))
        //         .collect::<Vec<_>>()
        //         .join("\n")
        // );
        // Map back the integer indices to the original expressions if substitutions took place
        match Self::parse_form_output_expression(&processed_form_str) {
            Ok(mut processed) => {
                processed = processed
                    .replace(vk_parse!("rat(x_,y_)").unwrap().to_pattern())
                    .with(vk_parse!("(x_/y_)").unwrap().to_pattern());
                processed = processed
                    .replace(vk_parse!("rat(x_)").unwrap().to_pattern())
                    .with(vk_parse!("x_").unwrap().to_pattern());
                processed = processed
                    .replace(vk_parse!("ep").unwrap().to_pattern())
                    .with(vk_parse!(&settings.epsilon_symbol).unwrap().to_pattern());
                processed = processed
                    .replace(vk_parse!("pi").unwrap().to_pattern())
                    .with(Atom::var(Symbol::PI).to_pattern());

                processed = processed
                    .replace(vk_parse!("g(idx1_,idx2_)").unwrap().to_pattern())
                    .when(
                        &(Condition::from((vk_symbol!("idx1_"), number_condition()))
                            & Condition::from((vk_symbol!("idx2_"), number_condition()))),
                    )
                    .with(
                        vk_parse!(format!("{}(idx1_,idx2_)", METRIC_SYMBOL).as_str())
                            .unwrap()
                            .to_pattern(),
                    );
                processed = processed
                    .replace(vk_parse!("g(v1_,v2_)").unwrap().to_pattern())
                    .when(
                        &(Condition::from((vk_symbol!("v1_"), symbol_condition()))
                            & Condition::from((vk_symbol!("v2_"), symbol_condition()))),
                    )
                    .with(vk_parse!("dot(v1_,v2_)").unwrap().to_pattern());
                processed = processed
                    .replace(
                        vk_parse!("g(v1_(args1_),v2_(args2_))")
                            .unwrap()
                            .to_pattern(),
                    )
                    .when(
                        &(Condition::from((vk_symbol!("v1_"), symbol_condition()))
                            & Condition::from((vk_symbol!("v2_"), symbol_condition()))),
                    )
                    .with(
                        vk_parse!("dot(v1_(args1_),v2_(args2_))")
                            .unwrap()
                            .to_pattern(),
                    );
                processed = processed
                    .replace(vk_parse!("vec1(vec_,idx_)").unwrap().to_pattern())
                    .when(Condition::from((vk_symbol!("vec_"), symbol_condition())))
                    .with(vk_parse!("vec_(idx_)").unwrap().to_pattern());

                // Undo the temporary float marker wrapping the rationalized coefficients and map them back to floats
                // TODO VHBEN
                processed = processed.replace_map(|term, _ctx, out| {
                    if let AtomView::Fun(c) = term
                        && c.get_symbol() == S.float_marker
                        && c.get_nargs() == 4
                    {
                        let mut args_iter = c.iter();
                        let re = args_iter.next().unwrap();
                        let im = args_iter.next().unwrap();
                        let re_prec = args_iter.next().unwrap();
                        let im_prec = args_iter.next().unwrap();

                        **out = Atom::num(Complex::new(
                            Rational::try_from(re)
                                .unwrap()
                                .to_multi_prec_float(u32::try_from(re_prec).unwrap()),
                            Rational::try_from(im)
                                .unwrap()
                                .to_multi_prec_float(u32::try_from(im_prec).unwrap()),
                        ));
                    }
                });

                // Convert vectors back from pi(j) notation to p(i,j) notation
                for (s, t) in vector_mapping.iter() {
                    processed = processed.replace(s.to_pattern()).with(t.to_pattern());
                }

                for (i_index, user_i) in indices.iter().enumerate() {
                    let mut replacements = vec![];
                    let integer_index =
                        Atom::num((FORM_REPLACEMENT_INDEX_SHIFT + 1 + i_index as u64) as i64);
                    let id = Atom::var(vk_symbol!("id_"));
                    for vec in [LOOP_MOMENTUM_SYMBOL, EXTERNAL_MOMENTUM_SYMBOL] {
                        let vector = vk_symbol!(vec);
                        replacements.push(Replacement::new(
                            function!(vector, &id, &integer_index).to_pattern(),
                            function!(vector, &id, user_i),
                        ));
                    }
                    let metric = vk_symbol!(METRIC_SYMBOL);
                    replacements.push(Replacement::new(
                        function!(metric, &integer_index, &id).to_pattern(),
                        function!(metric, &id, user_i),
                    ));
                    replacements.push(Replacement::new(
                        function!(metric, &id, &integer_index).to_pattern(),
                        function!(metric, user_i, &id),
                    ));
                    processed = processed.replace_multiple(&replacements);
                }

                Ok(processed)
            }
            Err(err) => Err(VakintError::FormOutputParsingError(form_output, err)),
        }
    }

    pub fn run_pysecdec(
        &self,
        settings: &VakintSettings,
        input: &[(String, String)],
        options: Vec<String>,
        clean: bool,
        reused_path: Option<String>,
        temporary_directory: Option<String>,
    ) -> Result<Vec<String>, VakintError> {
        citations::USED_CITATIONS.record(citations::CitationSource::PySecDec);
        let mut generate_pysecdec_sources = true;
        let tmp_dir = if let Some(reused_path_specified) = reused_path.as_ref() {
            let specified_dir = PathBuf::from(reused_path_specified);
            if !specified_dir.exists() {
                fs::create_dir_all(&specified_dir).map_err(|err| {
                    VakintError::PySecDecError(format!(
                        "Could not create user-specified directory to be reused for pySecDec run '{}': {}",
                        reused_path_specified.green(),
                        err
                    ))
                })?;
                warn!(
                    "User-specified directory '{}' not found, and was instead created and pySecDec sources will be regenerated.",
                    reused_path_specified
                );
            } else {
                generate_pysecdec_sources = false;
                warn!("{}",format!("User requested to re-use existing directory '{}' for the pysecdec run.\nThis is of course potentially unsafe and should be used for debugging only.\nRemove that directory to start clean.", reused_path_specified).red());
            }
            specified_dir
        } else {
            create_unique_tmp_run_dir(temporary_directory.as_deref(), "vakint_temp_pysecdec")?
        };

        if generate_pysecdec_sources {
            for input in input.iter() {
                fs::write(tmp_dir.join(&input.0), &input.1)?;
            }
        }

        let mut cmd = Command::new(settings.python_exe_path.as_str());
        cmd.arg(input[0].clone().0);
        for opt in options {
            cmd.arg(opt);
        }
        cmd.current_dir(&tmp_dir);

        if !clean {
            info!("Running {} with command: {:?}", "PySecDec".green(), cmd);
            info!(
                "You can follow the run live with 'tail -f follow_run.txt' in that temporary directory"
            );
            info!("Temporary run directory: {}", tmp_dir.display());
        } else {
            debug!("Running {} with command: {:?}", "PySecDec".green(), cmd);
            debug!(
                "You can follow the run live with 'tail -f follow_run.txt' in that temporary directory"
            );
            debug!("Temporary run directory: {}", tmp_dir.display());
        }

        let mut child = cmd.stderr(Stdio::piped()).stdout(Stdio::piped()).spawn()?;

        let stdout = child.stdout.take().unwrap();

        let reader = BufReader::new(stdout);
        let mut follow_file = File::create(tmp_dir.join("follow_run.txt"))?;

        for line in reader.lines() {
            let line_with_new_line = format!("{}\n", line?);
            follow_file.write_all(line_with_new_line.as_bytes())?;
            follow_file.flush()?;
        }

        let status = child.wait()?;

        if !ExitStatus::success(&status) {
            return Err(VakintError::FormError(
                "N/A".into(),
                "N/A".into(),
                format!("{:?}", cmd),
                tmp_dir.display().to_string(),
            ));
        }
        if !tmp_dir.join("out.txt").exists() {
            return Err(VakintError::MissingFormOutput(
                "N/A".into(),
                format!("{:?}", cmd),
                tmp_dir.display().to_string(),
            ));
        }
        let result = fs::read_to_string(tmp_dir.join("out.txt"))?;
        if clean {
            fs::remove_dir_all(&tmp_dir)?;
        }
        Ok(result.split('\n').map(|s| s.into()).collect::<Vec<_>>())
    }

    pub fn run_form(
        &self,
        settings: &VakintSettings,
        resources: &[String],
        input: (String, String),
        options: Vec<String>,
        clean: bool,
        temporary_directory: Option<String>,
    ) -> Result<String, VakintError> {
        citations::USED_CITATIONS.record(citations::CitationSource::Form);
        let tmp_dir = create_unique_tmp_run_dir(temporary_directory.as_deref(), "vakint_temp")?;

        for resource in resources.iter() {
            fs::write(tmp_dir.join(resource), FORM_SRC.get(resource).unwrap())?;
        }
        fs::write(tmp_dir.join(&input.0), &input.1)?;
        let mut cmd = Command::new(settings.form_exe_path.as_str());
        for opt in options {
            cmd.arg(opt);
        }
        cmd.arg(input.0);
        cmd.current_dir(&tmp_dir);
        if !clean {
            info!("Running {} with command: {:?}", "FORM".green(), cmd);
            info!("Temporary run directory: {}", tmp_dir.display());
        } else {
            debug!("Running {} with command: {:?}", "FORM".green(), cmd);
            debug!("Temporary run directory: {}", tmp_dir.display());
        }
        //println!("Running {} with command: {:?}", "FORM".green(), cmd);
        let output = cmd.stderr(Stdio::piped()).stdout(Stdio::piped()).output()?;
        // println!("{}", String::from_utf8_lossy(&output.stdout));
        if !ExitStatus::success(&output.status) {
            return Err(VakintError::FormError(
                String::from_utf8_lossy(&output.stderr).into(),
                String::from_utf8_lossy(&output.stdout).into(),
                format!("{:?}", cmd),
                tmp_dir.display().to_string(),
            ));
        }
        if !tmp_dir.join("out.txt").exists() {
            return Err(VakintError::MissingFormOutput(
                String::from_utf8_lossy(&output.stderr).into(),
                format!("{:?}", cmd),
                tmp_dir.display().to_string(),
            ));
        }

        let result = fs::read_to_string(tmp_dir.join("out.txt"))?;
        if clean {
            fs::remove_dir_all(&tmp_dir)?;
        }
        Ok(result)
    }

    fn next_input_request_id() -> Result<u64, VakintError> {
        static NEXT_INPUT_REQUEST: AtomicU64 = AtomicU64::new(0);
        NEXT_INPUT_REQUEST
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |next| {
                next.checked_add(1)
            })
            .map_err(|_| {
                VakintError::NumericalPrecision("input source identity space exhausted".into())
            })
    }

    fn stage_input_coefficients(
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        let mut has_float = false;
        input.visitor(&mut |part| {
            if precision::legacy::is_source(part) { return false; }
            has_float |= matches!(part, AtomView::Num(number) if matches!(number.get_coeff_view(), CoefficientView::Float(_, _)));
            !has_float
        });
        if !has_float {
            return Ok(input.to_owned());
        }
        precision::legacy::stage_table(
            function!(
                vk_symbol!("input_coefficient"),
                Self::next_input_request_id()?
            )
            .as_view(),
            input,
            vk_symbol!(&settings.epsilon_symbol),
            "Input",
        )
    }

    pub fn evaluate(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        self.evaluate_observed(settings, input, None)
    }

    /// Evaluate through the ordinary backend and finalization path with per-call diagnostics.
    pub fn evaluate_with_monitor(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        monitor: &RustRedMonitorOptions,
    ) -> Result<Atom, VakintError> {
        MonitorContext::run(monitor, |observation, _progress| {
            self.evaluate_observed(settings, input, observation)
        })
    }

    fn evaluate_observed(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            MonitorContext::native_phase(
                observation.as_deref_mut(),
                "input_preparation",
                Some(input),
            );
            self.validate_settings(settings)?;
            let mut vakint_expr = VakintExpression::prepare_evaluation(input, self, settings)?;
            MonitorContext::native_phase(observation.as_deref_mut(), "canonicalization", None);
            vakint_expr.canonicalize_with(settings, false, |input| {
                self.match_topologies_for_settings(settings, input)
            })?;
            MonitorContext::native_phase(observation.as_deref_mut(), "tensor_prepass", None);
            for term in &mut vakint_expr.0 {
                resources::checkpoint()?;
                term.tensor_reduce_for_evaluation(self, settings)?;
            }
            vakint_expr.evaluate_integral_observed(self, settings, observation.as_deref_mut())?;
            let result: Atom = vakint_expr.into();
            MonitorContext::native_phase(observation, "evaluation_output", Some(result.as_view()));
            Ok(result)
        })
    }

    pub fn to_canonical(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        short_form: bool,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            let mut vakint_expr = VakintExpression::try_from(input)?;
            vakint_expr.canonicalize_with(settings, short_form, |input| {
                self.match_topologies_for_settings(settings, input)
            })?;
            Ok(vakint_expr.into())
        })
    }

    pub fn tensor_reduce(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            let mut vakint_expr = VakintExpression::try_from(input)?;
            vakint_expr.tensor_reduce(self, settings)?;
            Ok(vakint_expr.into())
        })
    }

    pub fn evaluate_integral(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        settings.with_resource_budget(|| {
            self.validate_settings(settings)?;
            let mut vakint_expr = VakintExpression::prepare_evaluation(input, self, settings)?;
            vakint_expr.evaluate_integral(self, settings)?;
            Ok(vakint_expr.into())
        })
    }
}

#[cfg(test)]
mod tests {
    use symbolica::parse_lit;

    use super::*;

    #[test]
    fn resource_errors_expose_the_original_standard_error_source() {
        use std::error::Error;

        let original = resources::ResourceLimitError::Memory {
            rss_bytes: 0,
            additional_bytes: 2,
            limit_bytes: 1,
        };
        let error = VakintError::from(original.clone());
        assert_eq!(error.to_string(), original.to_string());
        assert_eq!(
            error
                .source()
                .unwrap()
                .downcast_ref::<resources::ResourceLimitError>(),
            Some(&original),
        );
        assert_eq!(
            RustRedEvaluationError::resource_cause(&error),
            Some(&original)
        );
    }

    #[test]
    fn resource_causes_survive_nested_io_wrappers() {
        let causes = [
            resources::ResourceLimitError::Memory {
                rss_bytes: 8,
                additional_bytes: 3,
                limit_bytes: 10,
            },
            resources::ResourceLimitError::Elapsed {
                elapsed: std::time::Duration::from_secs(2),
                limit: std::time::Duration::from_secs(1),
            },
        ];
        for original in causes {
            let wrapped = std::io::Error::other(original.clone());
            assert_eq!(
                RustRedEvaluationError::resource_cause(&wrapped),
                Some(&original)
            );
            let nested = std::io::Error::other(wrapped);
            assert_eq!(
                RustRedEvaluationError::resource_cause(&nested),
                Some(&original)
            );
        }
        assert!(
            RustRedEvaluationError::resource_cause(&std::io::Error::other("ordinary IO failure"))
                .is_none()
        );
    }

    #[test]
    fn resource_budgets_fail_gracefully_and_preserve_engine_reuse() {
        let engine = Vakint::new().unwrap();
        let input = vk_parse!("topo(I1L(mass_squared,1))").unwrap();
        let settings = VakintSettings {
            rustred_resources: ResourceLimits {
                max_elapsed: Some(std::time::Duration::ZERO),
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(matches!(
            engine.to_canonical(&settings, input.as_view(), true),
            Err(VakintError::ResourceLimit(
                resources::ResourceLimitError::Elapsed { .. }
            ))
        ));
        let settings = VakintSettings::default();
        let outer = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1),
                ..Default::default()
            },
            || Ok(2),
        );
        assert!(matches!(
            outer.run(|| engine.to_canonical(&settings, input.as_view(), true)),
            Err(VakintError::ResourceLimit(
                resources::ResourceLimitError::Memory { .. }
            ))
        ));
        assert!(
            engine
                .to_canonical(&settings, input.as_view(), true)
                .is_ok()
        );
    }

    #[test]
    fn numerator_identification_preserves_authenticated_scalar_masses() {
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            verify_numerator_identification: true,
            ..Default::default()
        };
        for (mass, numerator) in [
            ("mass_squared", "dot(k(1),k(1))-mass_squared"),
            ("mass_squared", "dot(k(1),k(1))/(mass_squared-1)"),
            ("mass^2", "dot(k(1),k(1))-mass^2"),
            ("mass_a+2*mass_b^2", "(dot(k(1),k(1))-mass_a-2*mass_b^2)^2"),
            ("mass^2", "dot(p(7),p(7))*(dot(k(1),k(1))-mass^2)"),
        ] {
            let integral = vk_parse!(format!("topo(I1L({mass},2))")).unwrap();
            let mut matched = engine
                .topologies
                .match_topologies_to_user_input(integral.as_view(), false)
                .unwrap()
                .unwrap();
            assert_eq!(
                matched
                    .canonical_expression_substitutions
                    .get(&vk_parse!("msq(1)").unwrap()),
                Some(&vk_parse!(mass).unwrap()),
                "the matched mass must be preserved before verification"
            );
            matched.apply_replacement_rules().unwrap();
            let mut term = VakintTerm {
                integral,
                numerator: vk_parse!(numerator).unwrap(),
                vectors: Vec::new(),
            };
            let input = &term.numerator * &term.integral;
            assert_eq!(
                engine
                    .to_canonical(&settings, input.as_view(), true)
                    .unwrap(),
                input
            );
            let expected = term.clone();
            term.apply_numerator_replacement_rules(&matched, &settings)
                .unwrap();
            assert_eq!(term, expected, "mass: {mass}; numerator: {numerator}");
            for parameter in vk_parse!(mass).unwrap().get_all_symbols(false) {
                assert!(term.numerator.contains_symbol(parameter));
            }
        }
    }

    #[test]
    fn numerator_identification_rejects_unmatched_momenta_and_scalar_parameters() {
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            verify_numerator_identification: true,
            ..Default::default()
        };
        let integral = vk_parse!("topo(I1L(mass_squared,2))").unwrap();
        let mut matched = engine
            .topologies
            .match_topologies_to_user_input(integral.as_view(), false)
            .unwrap()
            .unwrap();
        matched.apply_replacement_rules().unwrap();
        for numerator in [
            "dot(k(99),k(99))-mass_squared",
            "dot(q(1),q(1))-mass_squared",
            "dot(k(1),k(1))-undeclared_mass_squared",
            "(mass_squared-1)*dot(k(99),k(99))",
            "(mass_squared-1)*dot(q(1),q(1))",
            "(mass_squared-1)*undeclared_scalar",
            "(mass_squared-1)*mass_squared(1)",
        ] {
            let mut term = VakintTerm {
                integral: integral.clone(),
                numerator: vk_parse!(numerator).unwrap(),
                vectors: Vec::new(),
            };
            let original = term.clone();
            assert!(matches!(
                term.apply_numerator_replacement_rules(&matched, &settings),
                Err(VakintError::NumeratorNotReplaced(_, _))
            ));
            assert_eq!(term, original);
        }
    }

    #[test]
    fn numerator_identification_rejects_momentum_dependent_masses() {
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            verify_numerator_identification: true,
            ..Default::default()
        };
        for mass in ["dot(k(99),k(99))", "dot(p(1),p(1))"] {
            let integral = vk_parse!(format!("topo(I1L({mass},2))")).unwrap();
            let mut matched = engine
                .topologies
                .match_topologies_to_user_input(integral.as_view(), false)
                .unwrap()
                .unwrap();
            matched.apply_replacement_rules().unwrap();
            let mut term = VakintTerm {
                integral,
                numerator: vk_parse!(format!("dot(k(1),k(1))-({mass})")).unwrap(),
                vectors: Vec::new(),
            };
            let original = term.clone();
            assert!(matches!(
                term.apply_numerator_replacement_rules(&matched, &settings),
                Err(VakintError::InvalidIntegralFormat(_))
            ));
            assert_eq!(term, original);
        }
    }

    #[test]
    fn evaluation_input_tracking_precedes_topology_collection() {
        Vakint::initialize_vakint_symbols();
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]),
            run_time_decimal_precision: 40,
            ..Default::default()
        };
        let input = vk_parse!("(1.25+x)*topo(I1L(muvsq,1))+(-1.25-x)*topo(I1L(muvsq,1))").unwrap();
        let prepared =
            VakintExpression::prepare_evaluation(input.as_view(), &engine, &settings).unwrap();
        assert_eq!(prepared.0.len(), 1);
        let numerator = &prepared.0[0].numerator;
        let sources = precision::legacy::collect_sources(numerator.as_view()).unwrap();
        assert_eq!(sources.len(), 2);
        let (value, _) = precision::evaluate_coefficient(
            numerator.as_view(),
            &BTreeMap::new(),
            &std::collections::HashMap::new(),
            &settings,
        )
        .unwrap();
        assert!(value.value().re.is_zero());
        assert!(!value.precision().real().exact_zero());
        assert!(value.to_atom().is_err());
    }

    #[test]
    fn input_source_identity_survives_distribution_but_separates_occurrences_and_requests() {
        Vakint::initialize_vakint_symbols();
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Symbolic,
                },
            )]),
            ..Default::default()
        };
        let shared = vk_parse!("1.25*(topo(I1L(muvsq,1))+topo(I1L(muvsq,2)))").unwrap();
        let prepared =
            VakintExpression::prepare_evaluation(shared.as_view(), &engine, &settings).unwrap();
        assert_eq!(prepared.0.len(), 2);
        let identities = |expression: &VakintExpression| {
            expression
                .0
                .iter()
                .map(|term| {
                    let sources =
                        precision::legacy::collect_sources(term.numerator.as_view()).unwrap();
                    assert_eq!(sources.len(), 1);
                    sources.into_keys().next().unwrap()
                })
                .collect::<Vec<_>>()
        };
        let shared_ids = identities(&prepared);
        assert_eq!(shared_ids[0], shared_ids[1]);
        assert!((&prepared.0[0].numerator - &prepared.0[1].numerator).is_zero());

        let retained = Atom::from(prepared.clone());
        let cloned =
            VakintExpression::prepare_evaluation(retained.as_view(), &engine, &settings).unwrap();
        assert_eq!(identities(&cloned), shared_ids);
        let independent_request =
            VakintExpression::prepare_evaluation(shared.as_view(), &engine, &settings).unwrap();
        assert_ne!(identities(&independent_request)[0], shared_ids[0]);

        let independent = vk_parse!("1.25*topo(I1L(muvsq,1))+1.25*topo(I1L(muvsq,2))").unwrap();
        let prepared =
            VakintExpression::prepare_evaluation(independent.as_view(), &engine, &settings)
                .unwrap();
        let separate_ids = identities(&prepared);
        assert_ne!(separate_ids[0], separate_ids[1]);
    }

    #[test]
    fn lower_public_evaluation_calls_keep_independent_input_identities() {
        Vakint::initialize_vakint_symbols();
        let engine = Vakint::new().unwrap();
        let method = EvaluationMethod::RustRed(RustRedEvaluationOptions {
            master_mode: RustRedMasterMode::Retain,
        });
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![method.clone()]),
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            ..Default::default()
        };
        let topology = vk_parse!("topo(I1L(muvsq,1))").unwrap();
        let mut matched = engine
            .match_topologies_for_settings(&settings, topology.as_view())
            .unwrap()
            .unwrap();
        matched.apply_replacement_rules().unwrap();
        let source = |value: &Atom| {
            let sources = precision::legacy::collect_sources(value.as_view()).unwrap();
            assert_eq!(sources.len(), 1);
            sources.into_keys().next().unwrap()
        };
        for route in ["method", "term", "expression"] {
            let evaluate = |numerator: Atom| {
                if route == "method" {
                    return method
                        .evaluate_integral(&engine, &settings, numerator.as_view(), &matched)
                        .unwrap();
                }
                let mut term = VakintTerm {
                    integral: topology.clone(),
                    numerator,
                    vectors: Vec::new(),
                };
                if route == "term" {
                    term.evaluate_integral(&engine, &settings).unwrap();
                    term.numerator
                } else {
                    let mut expression = VakintExpression(vec![term]);
                    expression.evaluate_integral(&engine, &settings).unwrap();
                    Atom::from(expression)
                }
            };
            let first = evaluate(vk_parse!("1.25`30").unwrap());
            let second = evaluate(vk_parse!("1.25`30").unwrap());
            let first_source = source(&first);
            assert_ne!(first_source, source(&second), "{route}");
            assert!(!(first.clone() - second).is_zero(), "{route}");
            assert_eq!(
                source(&evaluate(first_source.clone())),
                first_source,
                "{route}"
            );
            assert!((first.clone() - first).is_zero(), "{route}");
        }
    }

    #[test]
    fn pysecdec_rejects_retained_sources_before_any_backend_work() {
        Vakint::initialize_vakint_symbols();
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let numerator =
            Vakint::stage_input_coefficients(&settings, vk_parse!("1.25`30").unwrap().as_view())
                .unwrap();
        let error = engine
            .pysecdec_evaluate(
                &settings,
                numerator.as_view(),
                &ReplacementRules::default(),
                &PySecDecOptions::default(),
                0,
            )
            .unwrap_err();
        assert!(matches!(error, VakintError::InvalidNumerator(_)));
        assert!(error.to_string().contains("retained NumericalSource"));
    }

    #[test]
    fn pysecdec_input_preparation_restores_floats_and_preserves_user_symbols() {
        Vakint::initialize_vakint_symbols();
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::PySecDec(
                PySecDecOptions::default(),
            )]),
            ..Default::default()
        };
        let numerator = vk_parse!("1.25+input_coefficient_placeholder_0(0)").unwrap();
        let input = &numerator * vk_parse!("topo(I1L(muvsq,1))").unwrap();
        let prepared =
            VakintExpression::prepare_evaluation(input.as_view(), &engine, &settings).unwrap();
        assert_eq!(prepared.0.len(), 1);
        assert_eq!(prepared.0[0].numerator, numerator);
        assert!(
            !prepared.0[0]
                .numerator
                .contains_symbol(vk_symbol!("NumericalSource"))
        );
        assert!(
            !prepared.0[0]
                .numerator
                .contains_symbol(vk_symbol!("input_coefficient_placeholder_1"))
        );
    }

    #[test]
    fn compound_lorentz_dimensions_and_generated_dummy_indices_roundtrip() {
        Vakint::initialize_vakint_symbols();
        for epsilon_symbol in ["ε", "eta"] {
            let expected = vk_parse!(format!(
                "dot(p(1),p(2))*dot(k(1),k(1))/(4-2*{epsilon_symbol})"
            ))
            .unwrap();
            for source in [
                "p(1,11)*k(1,11)*k(1,22)*p(2,22)",
                "dot(p(1),k(1))*dot(k(1),p(2))",
            ] {
                let input = vk_parse!(source).unwrap();
                for use_dot_product_notation in [false, true] {
                    let settings = VakintSettings {
                        epsilon_symbol: epsilon_symbol.into(),
                        use_dot_product_notation,
                        ..Default::default()
                    };
                    let mut term = VakintTerm {
                        numerator: input.clone(),
                        integral: Atom::one(),
                        vectors: Vec::new(),
                    };
                    term.tensor_reduce_feynkit(&settings, false).unwrap();
                    let actual = if use_dot_product_notation {
                        term.numerator
                    } else {
                        Vakint::convert_to_dot_notation(term.numerator.as_view())
                    };
                    assert!(
                        (actual - &expected).expand().together().is_zero(),
                        "epsilon={epsilon_symbol}, dot notation={use_dot_product_notation}, input={input}"
                    );

                    let mut scalar_only = VakintTerm {
                        numerator: input.clone(),
                        integral: Atom::one(),
                        vectors: Vec::new(),
                    };
                    assert!(scalar_only.tensor_reduce_feynkit(&settings, true).is_err());
                    assert_eq!(scalar_only.numerator, input);
                }
            }
        }
    }

    #[test]
    fn closed_scalar_coefficients_restore_collision_safe_aliases() {
        Vakint::initialize_vakint_symbols();
        let input =
            vk_parse!("tensor_scalar_coefficient(0)*dot(k(1),k(1))^3*k(2,71)*k(2,72)").unwrap();
        let expected = vk_parse!(
            "tensor_scalar_coefficient(0)*dot(k(1),k(1))^3*dot(k(2),k(2))*g(71,72)/(4-2*ε)"
        )
        .unwrap();
        for use_dot_product_notation in [false, true] {
            let settings = VakintSettings {
                use_dot_product_notation,
                ..Default::default()
            };
            let mut term = VakintTerm {
                numerator: input.clone(),
                integral: Atom::one(),
                vectors: Vec::new(),
            };
            term.tensor_reduce_feynkit(&settings, false).unwrap();
            let actual = if use_dot_product_notation {
                term.numerator
            } else {
                Vakint::convert_to_dot_notation(term.numerator.as_view())
            };
            assert!(
                (actual - &expected).expand().together().is_zero(),
                "dot notation: {use_dot_product_notation}"
            );

            let mut rejected = VakintTerm {
                numerator: input.clone(),
                integral: Atom::one(),
                vectors: Vec::new(),
            };
            assert!(rejected.tensor_reduce_feynkit(&settings, true).is_err());
            assert_eq!(rejected.numerator, input);
        }
    }

    #[test]
    fn repeated_closed_invariants_preserve_user_aliases_and_mixed_tensors() {
        Vakint::initialize_vakint_symbols();
        let scalar = (0..24)
            .map(|power| {
                vk_parse!(format!(
                    "c{power}*dot(k(1),k(1))^{power}*dot(k(1),k(2))^{}",
                    power % 3
                ))
                .unwrap()
            })
            .sum::<Atom>()
            * vk_parse!("tensor_scalar_coefficient(0)+tensor_scalar_coefficient(1)*alpha").unwrap();
        let input = &scalar * vk_parse!("k(3,71)*k(3,72)").unwrap();
        let expected = &scalar * vk_parse!("dot(k(3),k(3))*g(71,72)/(4-2*ε)").unwrap();
        for use_dot_product_notation in [false, true] {
            let settings = VakintSettings {
                use_dot_product_notation,
                ..Default::default()
            };
            let mut term = VakintTerm {
                numerator: input.clone(),
                integral: Atom::one(),
                vectors: Vec::new(),
            };
            term.tensor_reduce_feynkit(&settings, false).unwrap();
            let actual = if use_dot_product_notation {
                term.numerator
            } else {
                Vakint::convert_to_dot_notation(term.numerator.as_view())
            };
            assert!((actual - &expected).expand().together().is_zero());
            let mut rejected = VakintTerm {
                numerator: input.clone(),
                integral: Atom::one(),
                vectors: Vec::new(),
            };
            assert!(rejected.tensor_reduce_feynkit(&settings, true).is_err());
            assert_eq!(rejected.numerator, input);
        }
    }

    #[test]
    fn scalar_scope_closure_preserves_guarded_legacy_integrated_powers() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        for power in [2, 4, 6, 20] {
            let input = vk_parse!(format!("k(1,71)^{power}")).unwrap();
            let expected = vk_parse!(format!("dot(k(1),k(1))^{}", power / 2)).unwrap();
            let actual = VakintTerm::scalar_product_form(&settings, input.as_view()).unwrap();
            assert_eq!(actual, expected, "power {power}");
        }
        let over_rank = vk_parse!("k(1,71)^22").unwrap();
        assert!(matches!(
            VakintTerm::scalar_product_form(&settings, over_rank.as_view()),
            Err(VakintError::FeynKitTensorReduction(
                feynkit_tensor::TensorReductionError::UnsupportedRank {
                    rank: 22,
                    maximum: 20
                }
            ))
        ));
        for expression in [
            "k(1,71)^4*k(2,71)^2",
            "k(1,71)^4*p(1,71)",
            "k(1,71)^4*g(71,72)",
        ] {
            let input = vk_parse!(expression).unwrap();
            assert!(
                VakintTerm::scalar_product_form(&settings, input.as_view()).is_err(),
                "{expression}"
            );
        }
    }

    #[test]
    fn dot_conversion_preserves_composite_power_dummy_scopes() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let spectator = vk_parse!("(alpha+beta)^7").unwrap();
        for expression in [
            "(dot(k(1),k(1))-1)^11",
            "(dot(k(1),k(2))+dot(k(1),k(1)))^2*(dot(k(1),k(2))-1)^3",
            "(dot(k(1),k(2))+1)^-2",
        ] {
            let scalar = &spectator * vk_parse!(expression).unwrap();
            let converted = Vakint::convert_from_dot_notation(scalar.as_view());
            assert!(converted.contains(&spectator));
            assert_eq!(
                Vakint::convert_from_dot_notation(converted.as_view()),
                converted,
                "dot conversion must be idempotent: {expression}"
            );
            // Inverting a reciprocal exposes its repeated contractions without
            // asking the scalar service to integrate a non-polynomial numerator.
            let (indexed, expected) = if expression.ends_with("^-2") {
                (converted.pow(-1), scalar.pow(-1))
            } else {
                (converted, scalar)
            };
            let restored = VakintTerm::scalar_product_form(&settings, indexed.as_view()).unwrap();
            assert!(
                (restored - expected).expand().together().is_zero(),
                "{expression}"
            );
        }
        let indexed = vk_parse!("(p(5,dot_dummy_ind(9))+k(1,dot_dummy_ind(9)))^2").unwrap();
        assert_eq!(
            Vakint::convert_from_dot_notation(indexed.as_view()),
            indexed
        );
    }

    #[test]
    fn form_roundtrip_restores_epsilon_in_saved_lorentz_indices() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings {
            epsilon_symbol: "form_saved_index::epsilon".into(),
            ..Default::default()
        };
        let epsilon = vk_parse!(&settings.epsilon_symbol).unwrap();
        let mu = vk_parse!("form_saved_index::mu").unwrap();
        let index = function!(
            Minkowski::selfless_symbol(),
            Atom::num(4) - 2 * &epsilon,
            &mu
        );
        let coefficient = Atom::one() + &epsilon;
        let source = function!(vk_symbol!("p1"), &index);
        let target = function!(S.p, 1, &index);
        let metric = &coefficient * function!(vk_symbol!(METRIC_SYMBOL), &index, &mu);
        for (expression, expected, mapping) in [
            (
                &coefficient * function!(vk_symbol!("vec1"), vk_parse!("p1").unwrap(), &index),
                &coefficient * &target,
                BTreeMap::from([(source, target)]),
            ),
            (metric.clone(), metric, BTreeMap::default()),
        ] {
            for substitute_indices in [false, true] {
                let (_, serialized, indices) = vakint
                    .prepare_expression_for_form(
                        &settings,
                        expression.clone(),
                        substitute_indices,
                        &[],
                    )
                    .unwrap();
                assert_eq!(indices.is_empty(), !substitute_indices);
                assert!(
                    indices
                        .iter()
                        .all(|index| !index.contains_symbol(vk_symbol!("ep")))
                );
                let restored = vakint
                    .process_form_output(&settings, serialized, indices, mapping.clone())
                    .unwrap();
                assert_eq!(restored, expected);
                assert!(!restored.contains_symbol(vk_symbol!("ep")));
            }
        }
        let declaration = Atom::var(Minkowski::selfless_symbol()).to_canonical_string();
        assert!(declaration.contains("spenso::self_dual"));
        let conflicting =
            declaration.replacen("spenso::self_dual", "form_saved_index::conflicting_tag", 1);
        assert!(
            vakint
                .process_form_output(
                    &settings,
                    format!("{conflicting}(4,form_saved_index::mu)"),
                    Vec::new(),
                    BTreeMap::default(),
                )
                .is_err()
        );
    }

    #[test]
    fn form_serialization_handles_unit_imaginary_coefficients() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let x = vk_parse!("form_imaginary_test::x").unwrap();
        let y = vk_parse!("form_imaginary_test::y").unwrap();
        let imaginary = Atom::i();
        for expression in [
            imaginary.clone(),
            -&imaginary,
            &imaginary * &x,
            -&imaginary * &x,
            &imaginary / 3,
            &imaginary * &x * Atom::num((2, 3)),
            function!(
                vk_symbol!("form_imaginary_test"),
                &imaginary,
                1 - &imaginary
            ),
            (1 + &imaginary) * (&x + &y),
            imaginary.clone().pow(&x),
            x.pow(&imaginary),
        ] {
            for substitute_indices in [false, true] {
                let (_, serialized, indices) = vakint
                    .prepare_expression_for_form(
                        &settings,
                        expression.clone(),
                        substitute_indices,
                        &[],
                    )
                    .unwrap();
                assert!(!serialized.contains('𝑖'), "{serialized}");
                assert!(!serialized.contains("**"), "{serialized}");
                let restored = vakint
                    .process_form_output(
                        &settings,
                        serialized.clone(),
                        indices,
                        BTreeMap::default(),
                    )
                    .unwrap();
                assert_eq!(restored, expression, "FORM input: {serialized}");
            }
        }
    }

    #[test]
    fn form_serialization_preserves_complex_precedence_and_symbol_identity() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let custom = symbolica::symbol!(
            "user_space::serialization_custom",
            print = |_atom, _options, _state| Some("display_only".into())
        );
        let real = symbolica::symbol!("user_space::serialization_real"; Real);
        let antisymmetric = symbolica::symbol!("user_space::serialization_antisymmetric";
            Antisymmetric;
            print = |_atom, _options, _state| Some("display_only_antisymmetric".into()),
            tags = ["user_space::tensor_tag"]);
        let complex = Atom::one() + Atom::i() * 2;
        let custom_coefficient = function!(custom, &complex) * Atom::var(real);
        for expression in [
            complex.clone() * vk_parse!("user_space::x+user_space::y").unwrap(),
            complex.clone()
                * vk_parse!("numeric_group_0(1)+numeric_group_1_suffix+numeric_group_2").unwrap(),
            complex.clone().pow(vk_parse!("user_space::power").unwrap()),
            vk_parse!("user_space::x").unwrap().pow(&complex),
            custom_coefficient,
            function!(
                antisymmetric,
                vk_parse!("user_space::y").unwrap(),
                vk_parse!("user_space::x").unwrap()
            ) * &complex,
            complex
                * vk_parse!("tensor(user_space::chain(mink(D,mu)),mink(D,mu))*k(1,mink(D,mu))")
                    .unwrap(),
        ] {
            for substitute_indices in [false, true] {
                let (_, serialized, indices) = vakint
                    .sanitize_user_expressions(
                        &settings,
                        expression.as_view(),
                        substitute_indices,
                        &[],
                    )
                    .unwrap();
                let restored = vakint
                    .process_form_output(
                        &settings,
                        serialized.clone(),
                        indices,
                        BTreeMap::default(),
                    )
                    .unwrap();
                assert_eq!(restored, expression, "FORM input: {serialized}");
            }
        }
    }

    #[test]
    #[allow(clippy::unnecessary_operation)]
    fn undo_dots() {
        // Force lazy symbol initialization before parsing literal `dot(...)` in this test.
        S.dot;
        let expr = parse_lit!(
            -UFO::GC_10
                ^ 4 * (gammalooprs::ε - 2)
                ^ -1 * dot(p(0), p(1)) * dot(k(1), k(2))
                ^ 2 * gammalooprs::tag(0, 1) - UFO::GC_10
                ^ 4 * (gammalooprs::ε - 2)
                ^ -1 * dot(p(0), p(1)) * dot(k(1), k(2)) * dot(k(2), k(2)) * gammalooprs::tag(0, 0)
        );

        let a = Vakint::convert_from_dot_notation(expr.as_view());
        println!("a = {}", a);
    }

    #[test]
    fn rustred_scalar_evaluation_is_opt_in_and_form_independent() {
        let default_order = EvaluationOrder::default();
        assert_eq!(default_order.0.len(), 4);
        assert!(matches!(
            default_order.0.as_slice(),
            [
                EvaluationMethod::AlphaLoop(_),
                EvaluationMethod::MATAD(_),
                EvaluationMethod::FMFT(_),
                EvaluationMethod::PySecDec(_),
            ]
        ));

        let rustred_order = EvaluationOrder::rustred_only();
        let [EvaluationMethod::RustRed(options)] = rustred_order.0.as_slice() else {
            panic!("rustred_only must contain exactly the RustRed scalar evaluator");
        };
        assert!((options.master_mode != crate::RustRedMasterMode::Retain));
        assert!(rustred_order.0[0].dependencies().is_empty());
        assert_eq!(
            rustred_order.to_string(),
            "[RustRed (master_mode=numerical)]"
        );
    }

    #[test]
    fn short_form_conditions_preserve_fixed_slots_and_power_types() {
        let vakint = Vakint::new().unwrap();
        let valid = vk_parse!("topo(I3L_pinch_3_6(mass_squared,2,n,0,1,1,0))").unwrap();
        assert!(
            vakint
                .topologies
                .match_topologies_to_user_input(valid.as_view(), false)
                .unwrap()
                .is_some()
        );
        for invalid in [
            "topo(I3L_pinch_3_6(mass_squared,2,n,1,1,1,0))",
            "topo(I3L_pinch_3_6(mass_squared,2,n,0,1,1,1))",
            "topo(I3L_pinch_3_6(mass_squared,2,n+1,0,1,1,0))",
            "topo(I3L_pinch_3_6(mass_squared,2,f(n),0,1,1,0))",
        ] {
            assert!(
                vakint
                    .topologies
                    .match_topologies_to_user_input(vk_parse!(invalid).unwrap().as_view(), false)
                    .unwrap()
                    .is_none(),
                "invalid short form was admitted: {invalid}"
            );
        }
    }

    #[test]
    fn short_form_conditions_accept_all_literal_powers() {
        Vakint::initialize_vakint_symbols();
        let integral = Integral::new(
            1,
            Some(vk_parse!("topo(prop(1,edge(1,1),k(1),msq(1),pow(1)))").unwrap()),
            Some(vk_parse!("fixed_one_loop(msq(1),2)").unwrap()),
            EvaluationOrder::empty(),
        )
        .unwrap();
        for (input, expected) in [
            ("topo(fixed_one_loop(mass_squared,2))", true),
            ("topo(fixed_one_loop(mass_squared,3))", false),
        ] {
            assert_eq!(
                integral
                    .match_integral_to_short_user_input(vk_parse!(input).unwrap().as_view())
                    .unwrap()
                    .is_some(),
                expected
            );
        }
    }

    #[test]
    fn form_output_restores_vector_wrappers_and_decorated_indices() {
        let vakint = Vakint::new().unwrap();
        let index = vk_parse!("user_space::mink4(4,11)").unwrap();
        let expected = function!(S.p, Atom::num(1), &index);
        let output = vakint
            .process_form_output(
                &VakintSettings::default(),
                format!("vec1(p1,{})", FORM_REPLACEMENT_INDEX_SHIFT + 1),
                vec![index.clone()],
                BTreeMap::from([(function!(vk_symbol!("p1"), &index), expected.clone())]),
            )
            .unwrap();
        assert_eq!(output, expected);
    }

    #[test]
    fn rustred_scalar_dispatch_reduces_a_matched_one_loop_integral() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let integral = vk_parse!("topo(prop(1,edge(1,1),k(1),muvsq,1))").unwrap();
        let mut integral_specs = vakint
            .topologies
            .match_topologies_to_user_input(integral.as_view(), false)
            .unwrap()
            .unwrap();
        integral_specs.apply_replacement_rules().unwrap();
        let numerator = Atom::num(1);
        let method = EvaluationMethod::RustRed(RustRedEvaluationOptions {
            master_mode: crate::RustRedMasterMode::Retain,
        });

        let reduced = method
            .evaluate_integral(&vakint, &settings, numerator.as_view(), &integral_specs)
            .unwrap();
        assert!(reduced.to_canonical_string().contains("RustRedMaster"));
    }

    #[test]
    fn rustred_scalar_contractions_preserve_free_powers_and_foreign_functions() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let source = vk_parse!("(p(1,11)*p(2,11)+1)^2*p(3,99)*foreign_scalar(1,2)^2").unwrap();
        let expected = vk_parse!("(dot(p(1),p(2))+1)^2*p(3,99)*foreign_scalar(1,2)^2").unwrap();
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, source.as_view()).unwrap(),
            expected
        );
        for invalid in ["p(1,11)^2*p(2,11)^2", "(p(1,11)*p(2,11))^2*p(3,11)"] {
            assert!(
                VakintTerm::scalar_product_form(&settings, vk_parse!(invalid).unwrap().as_view())
                    .is_err()
            );
        }
    }

    #[test]
    fn rustred_scalar_contractions_close_native_slots_and_preserve_free_interfaces() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings {
            epsilon_symbol: "ep".into(),
            ..Default::default()
        };
        let dimension = Atom::num(4) - 2 * vk_parse!("ep").unwrap();
        let slot =
            |dimension: &Atom, index| function!(Minkowski::selfless_symbol(), dimension, index);
        let compact = |id, dimension: &Atom| {
            function!(S.p, id, function!(Minkowski::selfless_symbol(), dimension))
        };
        let mu = slot(&dimension, vk_parse!("scalar_native_mu").unwrap());
        let nu = slot(&dimension, vk_parse!("scalar_native_nu").unwrap());
        let free = function!(
            S.p,
            3,
            slot(&dimension, vk_parse!("scalar_native_free").unwrap())
        );
        let closed = function!(S.p, 1, &mu) * function!(S.g, &mu, &nu) * function!(S.p, 2, &nu);
        let dot = function!(S.dot, compact(1, &dimension), compact(2, &dimension));
        let foreign = vk_parse!("foreign_scalar(1,2)^2").unwrap();
        for power in [1, 2] {
            let numerator = (closed.clone() + 1).pow(power) * &free * &foreign;
            let expected = (dot.clone() + 1).pow(power) * &free * &foreign;
            let result = VakintTerm::scalar_product_form(&settings, numerator.as_view()).unwrap();
            assert!((result - expected).expand().is_zero());
        }
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, free.as_view()).unwrap(),
            free
        );
        let metric_trace = function!(ETS.metric, &mu, &mu);
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, metric_trace.as_view()).unwrap(),
            dimension
        );
        let dummy = S.dot_dummy_ind(91);
        let expanded = function!(
            S.p,
            1,
            function!(Minkowski::selfless_symbol(), &dimension),
            &dummy
        ) * function!(
            S.p,
            2,
            function!(Minkowski::selfless_symbol(), &dimension),
            &dummy
        );
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, expanded.as_view()).unwrap(),
            dot
        );
        // Native source dimensions need not equal the integration dimension.
        let six = Atom::num(6);
        let eight = Atom::num(8);
        let mu_six = slot(&six, vk_parse!("scalar_native_six").unwrap());
        let nu_six = slot(&six, vk_parse!("scalar_native_six_open").unwrap());
        for source in [
            function!(S.g, &mu_six, &mu_six),
            function!(S.g, &mu_six, &nu_six),
            function!(ETS.metric, &mu_six, &nu_six),
            function!(ETS.metric, &mu_six, &nu_six).pow(2),
        ] {
            let error = VakintTerm::scalar_product_form(&settings, source.as_view()).unwrap_err();
            let message = error.to_string();
            assert!(message.contains("cannot retarget native metric dimension 6"));
            assert!(message.contains(&format!("to {dimension}")));
        }
        // The native metric constructor already resolves a literal trace exactly.
        let native_trace = function!(ETS.metric, &mu_six, &mu_six);
        assert_eq!(native_trace, six);
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, native_trace.as_view()).unwrap(),
            six
        );
        let unresolved_trace = function!(ETS.metric, &mu, &nu).pow(2);
        assert!(unresolved_trace.contains_symbol(ETS.metric));
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, unresolved_trace.as_view()).unwrap(),
            dimension
        );
        let equal_dimension = vk_parse!("2*(2-ep)").unwrap();
        let equal_slot = slot(
            &equal_dimension,
            vk_parse!("scalar_native_equal_dimension").unwrap(),
        );
        let equal_trace = function!(S.g, &equal_slot, &equal_slot);
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, equal_trace.as_view()).unwrap(),
            dimension
        );
        let distinct_dimension = &dimension + vk_parse!("1/10^100").unwrap();
        let distinct_slot = slot(
            &distinct_dimension,
            vk_parse!("scalar_native_distinct_dimension").unwrap(),
        );
        let distinct_trace = function!(S.g, &distinct_slot, &distinct_slot);
        assert!(VakintTerm::scalar_product_form(&settings, distinct_trace.as_view()).is_err());
        let free_eight = function!(
            S.p,
            3,
            slot(&eight, vk_parse!("scalar_native_eight").unwrap())
        );
        let mixed = function!(S.p, 1, &mu_six) * function!(S.p, 2, &mu_six) * &free_eight;
        let expected = function!(S.dot, compact(1, &six), compact(2, &six)) * &free_eight;
        assert_eq!(
            VakintTerm::scalar_product_form(&settings, mixed.as_view()).unwrap(),
            expected
        );
    }

    #[test]
    fn rustred_scalar_contractions_native_external_dot_numerics_are_checked() {
        let engine = Vakint::new().unwrap();
        let settings = VakintSettings {
            epsilon_symbol: "ep".into(),
            run_time_decimal_precision: 16,
            ..Default::default()
        };
        let dimension = vk_parse!("4-2*ep").unwrap();
        let compact = |id| function!(S.p, id, function!(Minkowski::selfless_symbol(), &dimension));
        let native = function!(S.dot, compact(1), compact(2));
        let original = native.clone();
        let momentum = |time, space| {
            Momentum::Real((
                Float::with_val(256, time),
                Float::with_val(256, space),
                Float::with_val(256, 0),
                Float::with_val(256, 0),
            ))
        };
        let external = HashMap::from_iter([(1, momentum(3, 1)), (2, momentum(5, 4))]);
        let bare = vk_parse!("dot(p(1),p(2))").unwrap();
        let (value, statistical) = engine
            .numerical_evaluation(
                &settings,
                bare.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                Some(&external),
            )
            .unwrap();
        assert!(statistical.is_none());
        let coefficient = &value.coefficients()[&0];
        assert_eq!(coefficient.value().re.to_rational(), 11);
        assert!(coefficient.precision().real().significant_digits().unwrap() >= 16);
        assert!(coefficient.precision().imaginary().exact_zero());
        for expression in [
            native.clone(),
            function!(S.dot, compact(1), function!(S.p, 2)),
        ] {
            let (resolved, _) = engine
                .numerical_evaluation(
                    &settings,
                    expression.as_view(),
                    &HashMap::default(),
                    &HashMap::default(),
                    Some(&external),
                )
                .unwrap();
            let resolved = &resolved.coefficients()[&0];
            assert_eq!(resolved.value(), coefficient.value());
            assert_eq!(
                resolved.precision().real().significant_digits(),
                coefficient.precision().real().significant_digits()
            );
            assert_eq!(
                resolved.precision().real().absolute_error(),
                coefficient.precision().real().absolute_error()
            );
            assert!(resolved.precision().imaginary().exact_zero());
        }
        let cancelled = &native - &bare;
        let (cancelled, _) = engine
            .numerical_evaluation(
                &settings,
                cancelled.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                Some(&external),
            )
            .unwrap();
        assert!(cancelled.coefficients().is_empty());
        assert!(cancelled.to_atom(vk_symbol!("ep")).unwrap().is_zero());
        let foreign = function!(
            S.dot,
            function!(S.p, 1, function!(Minkowski::selfless_symbol(), 6)),
            function!(S.p, 2, function!(Minkowski::selfless_symbol(), 6))
        );
        let error = engine
            .numerical_evaluation(
                &settings,
                foreign.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                Some(&external),
            )
            .unwrap_err();
        assert!(matches!(&error, VakintError::EvaluationError(message)
            if message.contains(&foreign.to_canonical_string())
                && message.contains("source dimension 6")
                && message.contains(&format!("requested dimension {dimension}"))));
        let error = engine
            .numerical_evaluation(
                &settings,
                native.as_view(),
                &HashMap::default(),
                &HashMap::default(),
                None,
            )
            .unwrap_err();
        assert!(matches!(&error, VakintError::EvaluationError(message)
            if message.contains(&native.to_canonical_string())
                && message.contains("missing external momenta")));
        assert_eq!(native, original);
        // Existing expression-valued parameter keys still own their subtrees.
        for (expression, expected) in [
            (native.clone(), 17),
            (native.clone() + 1, 19),
            (bare.clone(), 23),
            (bare.clone() + 1, 29),
        ] {
            let binding = HashMap::from_iter([(
                expression.to_canonical_string(),
                Float::with_val(256, expected),
            )]);
            let (value, _) = engine
                .numerical_evaluation(
                    &settings,
                    expression.as_view(),
                    &binding,
                    &HashMap::default(),
                    Some(&external),
                )
                .unwrap();
            assert_eq!(value.coefficients()[&0].value().re.to_rational(), expected);
            assert!(
                value.coefficients()[&0]
                    .precision()
                    .real()
                    .significant_digits()
                    .unwrap()
                    >= 16
            );
        }
        let subexpression_binding =
            HashMap::from_iter([(native.to_canonical_string(), Float::with_val(256, 17))]);
        let scaled = Atom::num(2) * &native;
        let (scaled, _) = engine
            .numerical_evaluation(
                &settings,
                scaled.as_view(),
                &subexpression_binding,
                &HashMap::default(),
                Some(&external),
            )
            .unwrap();
        assert_eq!(scaled.coefficients()[&0].value().re.to_rational(), 34);
        let real = HashMap::from_iter([(native.to_canonical_string(), Float::with_val(256, 11))]);
        let complex = HashMap::from_iter([(
            native.to_canonical_string(),
            Complex::new(Float::with_val(256, 7), Float::with_val(256, 3)),
        )]);
        let (value, _) = engine
            .numerical_evaluation(
                &settings,
                native.as_view(),
                &real,
                &complex,
                Some(&external),
            )
            .unwrap();
        let coefficient = &value.coefficients()[&0];
        assert_eq!(coefficient.value().re.to_rational(), 7);
        assert_eq!(coefficient.value().im.to_rational(), 3);
        assert!(coefficient.precision().real().significant_digits().unwrap() >= 16);
        assert!(
            coefficient
                .precision()
                .imaginary()
                .significant_digits()
                .unwrap()
                >= 16
        );
    }

    #[test]
    fn rustred_scalar_contractions_match_staged_evaluation_and_requested_notation() {
        let engine = Vakint::new().unwrap();
        let topology = vk_parse!("topo(I1L(muvsq,1))").unwrap();
        let numerator = vk_parse!("(p(1,11)*p(2,11)+1)*p(3,99)").unwrap();
        let scalar = vk_parse!("(dot(p(1),p(2))+1)*p(3,99)").unwrap();
        let free = vk_parse!("p(3,99)").unwrap();
        for mode in [RustRedMasterMode::Retain, RustRedMasterMode::Symbolic] {
            for dots in [false, true] {
                let settings = VakintSettings {
                    evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                        RustRedEvaluationOptions { master_mode: mode },
                    )]),
                    tensor_reduction_method: TensorReductionMethod::FeynKit,
                    use_dot_product_notation: dots,
                    form_exe_path: "/scalar-contraction-test-must-not-use-form".into(),
                    number_of_terms_in_epsilon_expansion: 2,
                    ..Default::default()
                };
                let input = &numerator * &topology;
                let direct = engine.evaluate(&settings, input.as_view()).unwrap();
                let staged = engine.tensor_reduce(&settings, input.as_view()).unwrap();
                let staged = engine
                    .evaluate_integral(&settings, staged.as_view())
                    .unwrap();
                let expected = engine
                    .evaluate_integral(&settings, (&scalar * &topology).as_view())
                    .unwrap();
                let mut free_preserved = false;
                direct.as_view().visitor(&mut |value| {
                    free_preserved |= value == free.as_view();
                    true
                });
                assert!(free_preserved);
                assert_eq!(direct.contains_symbol(S.dot), dots);
                assert_eq!(staged.contains_symbol(S.dot), dots);
                let canonical = |value: &Atom| {
                    VakintTerm::scalar_product_form(&settings, value.as_view()).unwrap()
                };
                assert_eq!(canonical(&direct), canonical(&expected));
                assert_eq!(canonical(&staged), canonical(&expected));
            }
        }
    }

    #[test]
    fn rustred_support_is_limited_to_shipped_matcher_classes() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let rustred = EvaluationMethod::RustRed(RustRedEvaluationOptions::default());
        let alphaloop = EvaluationMethod::AlphaLoop(AlphaLoopOptions::default());
        let representative_integrals = [
            vk_parse!("topo(prop(1,edge(1,1),k(1),muvsq,1))").unwrap(),
            vk_parse!(
                "topo(\
                    prop(1,edge(1,2),k(1),muvsq,1)\
                   *prop(2,edge(1,2),k(2),muvsq,1)\
                   *prop(3,edge(2,1),k(1)+k(2),muvsq,1)\
                )"
            )
            .unwrap(),
            vk_parse!(
                "topo(\
                    prop(1,edge(1,2),k(1),muvsq,1)\
                   *prop(2,edge(2,3),k(2),muvsq,1)\
                   *prop(3,edge(3,1),k(3),muvsq,1)\
                   *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                   *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                   *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
                )"
            )
            .unwrap(),
        ];

        for integral in representative_integrals {
            let mut integral_specs = vakint
                .topologies
                .match_topologies_to_user_input(integral.as_view(), false)
                .unwrap()
                .unwrap();
            integral_specs.apply_replacement_rules().unwrap();
            assert!(rustred.supports(&settings, &integral_specs.canonical_topology));
            assert!(alphaloop.supports(&settings, &integral_specs.canonical_topology));
        }
    }

    #[test]
    fn rustred_support_checks_concrete_mass_and_power_before_coverage() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let rustred = EvaluationMethod::RustRed(RustRedEvaluationOptions::default());
        let inadmissible = [
            vk_parse!("topo(prop(1,edge(1,1),k(1),0,1))").unwrap(),
            vk_parse!("topo(prop(1,edge(1,1),k(1),muvsq,1/2))").unwrap(),
        ];

        for integral in inadmissible {
            let mut integral_specs = vakint
                .topologies
                .match_topologies_to_user_input(integral.as_view(), false)
                .unwrap()
                .unwrap();
            integral_specs.apply_replacement_rules().unwrap();
            assert!(!rustred.supports(&settings, &integral_specs.canonical_topology));
        }

        // Coverage is advisory: a valid concrete integral remains routable even
        // when its powers lie outside the campaign's certified root domain.
        let integral = vk_parse!("topo(I2L(muvsq,9223372036854775807,1,1))").unwrap();
        let mut integral_specs = vakint
            .topologies
            .match_topologies_to_user_input(integral.as_view(), false)
            .unwrap()
            .unwrap();
        integral_specs.apply_replacement_rules().unwrap();
        assert!(rustred.supports(&settings, &integral_specs.canonical_topology));
    }
}
