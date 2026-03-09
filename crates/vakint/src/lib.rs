#![allow(mixed_script_confusables)]
pub mod alphaloop_numerics;
pub mod fmft;
pub mod fmft_numerics;
pub mod graph;
pub mod matad;
pub mod matad_numerics;
pub mod symbols;
pub mod topologies;
pub mod utils;

use crate::utils::set_precision_in_polynomial_atom;
use ahash::RandomState;
use eyre::Result;
use colored::Colorize;
use graph::Graph;
#[allow(unused)]
use log::{debug, info, warn};

use regex::Regex;
use rug::float::Constant;
use std::{
    collections::{BTreeMap, HashMap, HashSet, hash_map::Entry},
    env,
    f64::consts::LOG2_10,
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, Write},
    ops::Div,
    path::PathBuf,
    process::{Command, ExitStatus, Stdio},
    sync::{Arc, LazyLock, Mutex},
    vec,
};
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::{
    coefficient::{Coefficient, CoefficientView},
    domains::float::{FloatLike, RealLike},
};
#[cfg(feature = "symbolica_community_module")]
pub mod symbolica_community_module;

use crate::alphaloop_numerics::DIRECT_SUBSTITUTIONS;

#[allow(unused)]
use symbolica::{
    atom::{
        Atom, AtomCore, AtomView, FunctionBuilder, SliceType, Symbol, Var,
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
    printer::{AtomPrinter, PrintOptions},
    transformer::Transformer,
};
use utils::simplify_real;
use utils::vakint_macros::{vk_parse, vk_symbol};

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
    #[error("unknown vakint error")]
    Unknown,
}

pub fn params_from_f64(
    params: &HashMap<String, f64, ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<String, Float, ahash::RandomState> {
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
    params
        .iter()
        .map(|(k, v)| (k.clone(), Float::with_val(binary_prec, v)))
        .collect()
}

pub fn params_from_complex_f64(
    params: &HashMap<String, Complex<f64>, ahash::RandomState>,
    decimal_prec: u32,
) -> HashMap<String, Complex<Float>, ahash::RandomState> {
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
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
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
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
    let binary_prec: u32 = ((decimal_prec as f64) * LOG2_10).floor() as u32;
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

fn gt_condition(value: i64) -> WildcardRestriction {
    WildcardRestriction::Filter(Box::new(move |m| {
        if let Match::Single(AtomView::Num(a)) = m {
            a.get_coeff_view() > InlineNum::new(value, 1).as_num_view().get_coeff_view()
        } else {
            false
        }
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

        for (_n_id, node) in graph.nodes.iter() {
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
            if let Some(m1) = unwrapped_input
                .pattern_match(
                    short_expression_pattern,
                    Some(&apply_restriction_to_symbols(
                        (1..=self.n_props)
                            .map(|i_prop| vk_symbol!(format!("pow{}_", i_prop)))
                            .collect(),
                        &symbol_or_number(),
                    )),
                    None,
                )
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
                    if let Some(Match::Single(a)) =
                        m1.match_stack.get(vk_symbol!(format!("msq{}_", i_prop)))
                    {
                        replacement_rules.canonical_expression_substitutions.insert(
                            vk_parse!(format!("msq({})", i_prop).as_str()).unwrap(),
                            a.to_owned(),
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
        }
    }
}

impl EvaluationMethod {
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
        let result = match self {
            EvaluationMethod::AlphaLoop(opts) => {
                vakint.alphaloop_evaluate(settings, numerator, integral_specs, opts)
            }
            EvaluationMethod::MATAD(opts) => {
                vakint.matad_evaluate(settings, numerator, integral_specs, opts)
            }
            EvaluationMethod::FMFT(opts) => {
                vakint.fmft_evaluate(settings, numerator, integral_specs, opts)
            }
            EvaluationMethod::PySecDec(opts) => {
                vakint.pysecdec_evaluate(settings, numerator, integral_specs, opts)
            }
        }?;
        // Simplify logarithms and zero powers knowing that all arguments are real
        Ok(simplify_real(result.as_view()))
    }
}

#[derive(Debug, Clone)]
pub struct EvaluationOrder(pub Vec<EvaluationMethod>);

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

#[derive(Debug, Clone)]
pub struct VakintSettings {
    #[allow(unused)]
    pub epsilon_symbol: String,
    pub mu_r_sq_symbol: String,
    pub form_exe_path: String,
    pub python_exe_path: String,
    pub verify_numerator_identification: bool,
    pub integral_normalization_factor: LoopNormalizationFactor,
    pub run_time_decimal_precision: u32,
    pub allow_unknown_integrals: bool,
    pub clean_tmp_dir: bool,
    pub evaluation_order: EvaluationOrder,
    // This quantity is typically set equal to *one plus the maximum loop count* of the UV regularisation problem considered.
    // For example when considering a 2-loop problem, then:
    //   a) for the nested one-loop integrals appearing, the single pole, finite term *and* order-epsilon term will need to be considered.
    //   b) for the two-loop integrals, the double pole, single pole and finite terms will be needed, so again three terms
    pub number_of_terms_in_epsilon_expansion: i64,
    pub precision_for_input_float_rationalization: InputFloatRationalizationPrecision,
    pub use_dot_product_notation: bool,
    pub temporary_directory: Option<String>,
}

impl VakintSettings {
    pub fn get_integral_normalization_factor_atom(&self) -> Result<Atom, VakintError> {
        self.integral_normalization_factor.to_atom(self)
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
                "(exp(log_mu_sq)/(4*𝜋*exp(-EulerGamma)))^(eps*n_loops)".into()
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
            true,
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
    pub topologies: Topologies,
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

impl VakintTerm {
    pub fn evaluate_integral(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        let mut integral_specs = if let Some(replacement_rules) =
            vakint.topologies.match_topologies_to_user_input(
                self.integral.as_view(),
                settings.allow_unknown_integrals,
            )? {
            replacement_rules
        } else {
            return Err(VakintError::UnreckognizedIntegral(
                self.integral.to_string(),
            ));
        };

        integral_specs.apply_replacement_rules()?;
        self.apply_numerator_replacement_rules(&integral_specs, settings)?;

        let mut could_evaluate_integral = false;
        'eval: for evaluation_approach in settings.evaluation_order.0.iter() {
            if evaluation_approach.supports(settings, &integral_specs.canonical_topology) {
                let evaluated_integral = evaluation_approach.evaluate_integral(
                    vakint,
                    settings,
                    self.numerator.as_atom_view(),
                    &integral_specs,
                )?;
                self.numerator = simplify_real(evaluated_integral.as_view());
                self.integral = Atom::num(1);
                could_evaluate_integral = true;
                break 'eval;
            }
        }
        if !could_evaluate_integral {
            return Err(VakintError::NoEvaluationMethodFound(
                self.integral.to_string(),
                settings.number_of_terms_in_epsilon_expansion
                    - (integral_specs.canonical_topology.get_integral().n_loops as i64)
                    - 1,
            ));
        }
        Ok(())
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

        if vakint_settings.verify_numerator_identification && !matches!(test, Atom::Num(_)) {
            return Err(VakintError::NumeratorNotReplaced(
                EXTERNAL_MOMENTUM_SYMBOL.into(),
                test.to_string(),
            ));
        }

        self.numerator = new_numerator;

        Ok(())
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
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VakintExpression(pub Vec<VakintTerm>);

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

    pub fn canonicalize(
        &mut self,
        settings: &VakintSettings,
        topologies: &Topologies,
        short_form: bool,
    ) -> Result<(), VakintError> {
        for term in self.0.iter_mut() {
            if let Some(replacement_rules) = topologies.match_topologies_to_user_input(
                term.integral.as_view(),
                settings.allow_unknown_integrals,
            )? {
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
        for term in self.0.iter_mut() {
            term.tensor_reduce(vakint, settings)?;
        }
        Ok(())
    }

    pub fn evaluate_integral(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(), VakintError> {
        for term in self.0.iter_mut() {
            term.evaluate_integral(vakint, settings)?;
        }
        Ok(())
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

#[derive(Debug, Clone, Default)]
pub struct NumericalEvaluationResult(pub Vec<(i64, Complex<Float>)>);

impl fmt::Display for NumericalEvaluationResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            write!(f, "{}", "Empty result".green())
        } else {
            write!(
                f,
                "{}",
                self.0
                    .iter()
                    .map(|(power, float_eval)| format!(
                        "{} : {}",
                        format!(
                            "ε^{}{}",
                            match *power {
                                power if power > 0 => "+",
                                0 => " ",
                                power if power < 0 => "",
                                _ => unreachable!(),
                            },
                            power
                        )
                        .green(),
                        format!("{}", float_eval).blue()
                    ))
                    .collect::<Vec<_>>()
                    .join("\n")
            )
        }
    }
}

impl NumericalEvaluationResult {
    pub fn to_atom(&self, epsilon_symbol: Symbol) -> Atom {
        let mut res = Atom::Zero;
        for (exp, coeff) in self.get_epsilon_coefficients() {
            // res += (Atom::num(coeff.re) + Atom::var(S.cmplx_i) * Atom::num(coeff.im))
            //     * Atom::var(epsilon_symbol).pow(Atom::num(exp))
            res += Atom::num(coeff) * Atom::var(epsilon_symbol).pow(Atom::num(exp))
        }
        res
    }

    pub fn from_atom(
        input: AtomView,
        epsilon_symbol: Symbol,
        settings: &VakintSettings,
    ) -> Result<Self, VakintError> {
        let epsilon_coeffs = input.coefficient_list::<i8>(&[Atom::var(epsilon_symbol)]);

        let epsilon_coeffs_vec = epsilon_coeffs
            .iter()
            .map(|(eps_atom, coeff)| {
                if let Some(m) = eps_atom
                    .pattern_match(
                        &vk_parse!(format!("{}^n_", epsilon_symbol).as_str())
                            .unwrap()
                            .to_pattern(),
                        Some(&Condition::from((vk_symbol!("n_"), number_condition()))),
                        None,
                    )
                    .next()
                {
                    (
                        get_integer_from_atom(m.get(&vk_symbol!("n_")).unwrap().as_view()).unwrap(),
                        coeff,
                    )
                } else if *eps_atom == Atom::var(epsilon_symbol) {
                    (1, coeff)
                } else if *eps_atom == Atom::num(1) {
                    (0, coeff)
                } else {
                    panic!("Epsilon atom should be of the form {}^n_", epsilon_symbol)
                }
            })
            .collect::<Vec<_>>();

        let binary_prec = settings.get_binary_precision();
        let mut epsilon_coeffs_vec_floats = vec![];
        for (i64, coeff) in epsilon_coeffs_vec.iter() {
            epsilon_coeffs_vec_floats.push((
                *i64,
                match coeff.evaluate::<Atom,Complex<Float>,_>(
                    |x| Complex::new(x.to_multi_prec_float(binary_prec), Float::with_val(binary_prec, 0.)),
                    &HashMap::default(),
                    &HashMap::default()
                ) {
                    Ok(x) => x,
                    Err(e) => {
                        return Err(VakintError::EvaluationError(format!(
                            "Coefficients of the Laurent series are not numbers. Expression: '{}' | Error: '{}'",
                            coeff.to_canonical_string(),
                            e
                        )));
                    }
                },
            ));
        }

        epsilon_coeffs_vec_floats.sort_by(|(i1, _), (i2, _)| i1.cmp(i2));
        Ok(NumericalEvaluationResult(epsilon_coeffs_vec_floats))
    }

    pub fn from_vec(input: Vec<(i64, (String, String))>, settings: &VakintSettings) -> Self {
        let binary_prec = settings.get_binary_precision();
        NumericalEvaluationResult(
            input
                .iter()
                .map(|(eps_pwr, (re, im))| {
                    (
                        *eps_pwr,
                        Complex::new(
                            Float::parse(re.as_str(), Some(binary_prec)).unwrap(),
                            Float::parse(im.as_str(), Some(binary_prec)).unwrap(),
                        ),
                    )
                })
                .collect::<Vec<_>>(),
        )
    }

    pub fn get_epsilon_coefficients(&self) -> Vec<(i64, Complex<Float>)> {
        self.0.clone()
    }
    pub fn get_epsilon_coefficient(&self, power: i64) -> Complex<Float> {
        match self.0.iter().find(|(i, _f)| *i == power) {
            Some((_, f)) => f.clone(),
            None => {
                if self.0.is_empty() {
                    Complex::new(Float::with_val(53, 0.), Float::with_val(53, 0.))
                } else {
                    Complex::new(self.0[0].1.re.zero(), self.0[0].1.re.zero())
                }
            }
        }
    }

    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|(_, f)| f.is_zero())
    }

    pub fn aggregate_errors(&self, other: &NumericalEvaluationResult) -> NumericalEvaluationResult {
        let mut res = NumericalEvaluationResult::default();
        let mut orders = vec![];
        for (o, _val) in self.0.iter() {
            orders.push(*o);
        }
        for (o, _val) in other.0.iter() {
            orders.push(*o);
        }
        orders.sort();
        for o in orders {
            let (a, b) = (
                &self.get_epsilon_coefficient(o),
                &other.get_epsilon_coefficient(o),
            );

            res.0.push((0, (a * a + b * b).sqrt()));
        }
        res
    }

    pub fn does_approx_match(
        &self,
        other: &NumericalEvaluationResult,
        error: Option<&NumericalEvaluationResult>,
        threshold: f64,
        max_pull: f64,
    ) -> (bool, String) {
        let mut powers = HashSet::new();
        for (power, _) in self.0.iter() {
            powers.insert(*power);
        }
        for (power, _) in other.0.iter() {
            powers.insert(*power);
        }
        let mut power_vec = powers.iter().cloned().collect::<Vec<_>>();
        power_vec.sort();

        for power in power_vec {
            let (self_val, other_val) = (
                self.get_epsilon_coefficient(power),
                other.get_epsilon_coefficient(power),
            );
            let comparisons = [
                (
                    "real",
                    self_val.re.clone(),
                    other_val.re.clone(),
                    error.map(|e| e.get_epsilon_coefficient(power).re.clone()),
                ),
                (
                    "imaginary",
                    self_val.im.clone(),
                    other_val.im.clone(),
                    error.map(|e| e.get_epsilon_coefficient(power).im.clone()),
                ),
            ];
            for (part, r, o, e) in comparisons.iter() {
                let f_max_pull = Float::with_val(r.prec(), max_pull);
                let f_threshold = Float::with_val(r.prec(), threshold);
                let delta = (r.clone() - o).norm();
                if let Some(err) = e
                    && !err.norm().is_zero()
                    && delta > f_max_pull * err.norm()
                {
                    return (
                        false,
                        format!(
                            "{} part of ε^{} coefficient does not match within max pull: {} != {} (pull = {})",
                            part,
                            power,
                            r,
                            o,
                            delta / err.norm()
                        ),
                    );
                }
                let scale = (r.norm() + o.norm()) / Float::with_val(r.prec(), 2.0);
                if scale.is_zero() {
                    if delta > f_threshold {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within abs. error required: {} != {} (abs. error = {})",
                                part, power, r, o, delta
                            ),
                        );
                    }
                } else {
                    let rel_error = delta.div(scale);
                    if rel_error > f_threshold {
                        return (
                            false,
                            format!(
                                "{} part of ε^{} coefficient does not match within rel. error required: {} != {} (rel. error = {})",
                                part, power, r, o, rel_error
                            ),
                        );
                    }
                }
            }
        }
        (true, "matches!".into())
    }
}

impl Vakint {
    pub fn initialize_vakint_symbols() {
        // Force initialization of symbolica symbols with proper attributes
        LazyLock::force(&S);
    }

    pub fn new() -> Result<Self, VakintError> {
        // Force initialization of symbolica symbols with proper attributes
        LazyLock::force(&S);

        let topologies = Topologies::generate_topologies()?;

        Ok(Vakint { topologies })
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
        Vakint::full_numerical_evaluation(
            settings,
            expression,
            params_real,
            params_complex,
            externals,
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
    ) -> Atom {
        let (const_map_real, const_map_complex) =
            Vakint::get_constants_map(settings, params_real, params_complex, externals).unwrap();

        let mut res = Vakint::convert_to_dot_notation(integral);
        for (src, trgt) in const_map_real.iter() {
            res = res.replace(src.to_pattern()).with(Atom::num(trgt.clone()));
        }
        for (src, trgt) in const_map_complex.iter() {
            // res = res.replace(src.to_pattern()).with(
            //     (Atom::num(trgt.re.clone()) + Atom::var(S.cmplx_i) * Atom::num(trgt.im.clone()))
            //         .to_pattern(),
            // );
            res = res.replace(src.to_pattern()).with(Atom::num(trgt.clone()));
        }

        res
    }

    pub fn full_numerical_evaluation(
        settings: &VakintSettings,
        integral: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<(NumericalEvaluationResult, Option<NumericalEvaluationResult>), VakintError> {
        if integral
            .pattern_match(&S.error_flag.to_pattern(), None, None)
            .next()
            .is_none()
        {
            return Ok((
                Vakint::full_numerical_evaluation_without_error(
                    settings,
                    integral,
                    params_real,
                    params_complex,
                    externals,
                )?,
                None,
            ));
        }
        let (error_atom, integral_atom) =
            utils::split_linear_atom(integral, Atom::var(S.error_flag_symbol).as_view());
        let mut central = Vakint::full_numerical_evaluation_without_error(
            settings,
            integral_atom.as_view(),
            params_real,
            params_complex,
            externals,
        )?;
        let mut error = Vakint::full_numerical_evaluation_without_error(
            settings,
            error_atom.as_view(),
            params_real,
            params_complex,
            externals,
        )?;

        for (i, eval) in central.0.iter() {
            if !error.0.iter().any(|(j, _)| i == j) {
                error
                    .0
                    .push((*i, Complex::new(eval.re.zero(), eval.re.zero())));
            }
        }
        for (i, eval) in error.0.iter() {
            if !central.0.iter().any(|(j, _)| i == j) {
                central
                    .0
                    .push((*i, Complex::new(eval.re.zero(), eval.re.zero())));
            }
        }
        Ok((central, Some(error)))
    }

    pub fn full_numerical_evaluation_without_error(
        settings: &VakintSettings,
        integral: AtomView,
        params_real: &HashMap<String, Float, ahash::RandomState>,
        params_complex: &HashMap<String, Complex<Float>, ahash::RandomState>,
        externals: Option<&HashMap<usize, Momentum, ahash::RandomState>>,
    ) -> Result<NumericalEvaluationResult, VakintError> {
        let epsilon_coeffs =
            integral.coefficient_list::<i8>(&[Atom::var(vk_symbol!(&settings.epsilon_symbol))]);
        let epsilon_coeffs_vec = epsilon_coeffs
            .iter()
            .map(|(eps_atom, coeff)| {
                if let Some(m) = eps_atom
                    .pattern_match(
                        &vk_parse!(format!("{}^n_", settings.epsilon_symbol))
                            .unwrap()
                            .to_pattern(),
                        Some(&Condition::from((vk_symbol!("n_"), number_condition()))),
                        None,
                    )
                    .next()
                {
                    (
                        get_integer_from_atom(m.get(&vk_symbol!("n_")).unwrap().as_view()).unwrap(),
                        coeff,
                    )
                } else if *eps_atom == vk_parse!(&settings.epsilon_symbol).unwrap() {
                    (1, coeff)
                } else if *eps_atom == Atom::num(1) {
                    (0, coeff)
                } else {
                    panic!("Epsilon atom should be of the form ε^n_")
                }
            })
            .collect::<Vec<_>>();

        let binary_prec = settings.get_binary_precision();

        let (map_real, map_complex) =
            Vakint::get_constants_map(settings, params_real, params_complex, externals).unwrap();
        let mut map_view: HashMap<AtomView<'_>, Complex<Float>, RandomState> = map_complex
            .iter()
            .map(|(k, v)| (k.as_view(), v.clone()))
            .collect();
        for (k, v) in map_real.iter() {
            map_view.insert(
                k.as_view(),
                Complex::new(v.clone(), Float::with_val(binary_prec, 0.)),
            );
        }
        let mut epsilon_coeffs_vec_floats = vec![];
        for (i64, coeff) in epsilon_coeffs_vec.iter() {
            let coeff_processed = Vakint::convert_to_dot_notation(coeff.as_view());
            // for (k, v) in map_view.iter() {
            //     println!("{} -> {}", k, v);
            // }
            // println!("coeff_processed={}", coeff_processed);
            epsilon_coeffs_vec_floats.push((
                *i64,
                match coeff_processed.evaluate(
                    |x| {
                        Complex::new(
                            x.to_multi_prec_float(binary_prec),
                            Float::with_val(binary_prec, 0.),
                        )
                    },
                    &map_view,
                    &HashMap::default(),
                ) {
                    Ok(x) => x,
                    Err(e) => {
                        return Err(VakintError::EvaluationError(format!(
                            "Is some tensor structure left? | Expression: '{}' | Error: '{}'",
                            coeff.to_canonical_string(),
                            e
                        )));
                    }
                },
            ));
        }

        epsilon_coeffs_vec_floats.sort_by(|(i1, _), (i2, _)| i1.cmp(i2));
        Ok(NumericalEvaluationResult(epsilon_coeffs_vec_floats))
    }

    fn pysecdec_evaluate(
        &self,
        settings: &VakintSettings,
        input_numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &PySecDecOptions,
    ) -> Result<Atom, VakintError> {
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
                                        Atom::Var(s) => {
                                            if exp == Atom::num(2).as_view() {
                                                (Atom::Var(s.clone()), Atom::Pow(m.clone()))
                                            } else {
                                                panic!(
                                                    "Could not find mass symbol in integral:\n{}",
                                                    integral
                                                );
                                            }
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
            vars.insert(
                "lorentz_indices".into(),
                format!(
                    "[{}]",
                    sorted_lorentz_indices
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
                    true,
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
        let log_mu_sq = function!(
            Symbol::LOG,
            vk_parse!(settings.mu_r_sq_symbol.as_str()).unwrap()
        );

        let pysecdec_normalization_correction = vk_parse!(
            format!(
                "(  𝑖*(𝜋^((4-2*{eps})/2))\
                )^{n_loops}",
                eps = settings.epsilon_symbol,
                n_loops = integral.n_loops
            )
            .as_str()
        )
        .unwrap();

        let evaluated_integral = pysecdec_output
            .iter()
            .map(|out| match vk_parse!(out.replace("I", "𝑖").as_str()) {
                Ok(mut processed) => {
                    processed = processed
                        .replace(vk_parse!("ep").unwrap().to_pattern())
                        .with(vk_parse!(&settings.epsilon_symbol).unwrap().to_pattern());
                    processed = processed
                        * pysecdec_normalization_correction.to_owned()
                        * settings
                            .get_integral_normalization_factor_atom()?
                            .replace(S.n_loops.to_pattern())
                            .with(Atom::num(integral.n_loops as i64).to_pattern());

                    let expanded_evaluation = match processed.series(
                        vk_symbol!(settings.epsilon_symbol.as_str()),
                        Atom::Zero.as_atom_view(),
                        Rational::from(
                            settings.number_of_terms_in_epsilon_expansion
                                - (integral.n_loops as i64)
                                - 1,
                        ),
                        true,
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
            })
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
                        Atom::Var(s) => {
                            if exp == Atom::num(2).as_view() {
                                (Atom::Var(s.clone()), Atom::Pow(m.clone()))
                            } else {
                                return Err(VakintError::MalformedGraph(format!(
                                    "Could not find muV in graph:\n{}",
                                    short_integral_expression.to_canonical_string()
                                )));
                            }
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
            form_expression = form_expression
                .replace(
                    vk_parse!(format!("{}({},{})", v, v_id, idx.to_canonical_string()).as_str())
                        .unwrap()
                        .to_pattern(),
                )
                .with(
                    vk_parse!(
                        format!("vec1({}{},{})", v, v_id, idx.to_canonical_string()).as_str()
                    )
                    .unwrap()
                    .to_pattern(),
                );
            vector_mapping.insert(
                vk_parse!(format!("{}{}({})", v, v_id, idx.to_canonical_string()).as_str())
                    .unwrap(),
                vk_parse!(format!("{}({},{})", v, v_id, idx.to_canonical_string()).as_str())
                    .unwrap(),
            );
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
                    (
                        src,
                        (
                            set_precision_in_polynomial_atom(
                                trgt.as_view(),
                                vk_symbol!("ep"),
                                settings,
                            ),
                            condition.clone(),
                        ),
                    )
                })
                .collect::<Vec<_>>();
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
            evaluated_integral = r;
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
                .get_integral_normalization_factor_atom()?
                .as_view()
                .replace(S.n_loops.to_pattern())
                .with(Atom::num(integral.n_loops as i64).to_pattern());

        let expanded_evaluation = match evaluated_integral.series(
            vk_symbol!(settings.epsilon_symbol.as_str()),
            Atom::Zero.as_atom_view(),
            Rational::from(
                settings.number_of_terms_in_epsilon_expansion - (integral.n_loops as i64) - 1,
            ),
            true,
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
            .when(S.c_.filter(|a| a.to_atom() < 0))
            .with(S.dot_pow(S.a_, S.b_, -Atom::var(S.c_)).npow(-1))
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
                                            .add_args(&b.iter().collect::<Vec<_>>())
                                            .add_arg(&dummy_wrapped)
                                            .finish()
                                        * rest
                                }
                                (AtomView::Fun(a), AtomView::Var(b)) => {
                                    FunctionBuilder::new(a.get_symbol())
                                        .add_args(&a.iter().collect::<Vec<_>>())
                                        .add_arg(&dummy_wrapped)
                                        .finish()
                                        * FunctionBuilder::new(b.get_symbol())
                                            .add_arg(&dummy_wrapped)
                                            .finish()
                                        * rest
                                }
                                (AtomView::Fun(a), AtomView::Fun(b)) => {
                                    FunctionBuilder::new(a.get_symbol())
                                        .add_args(&a.iter().collect::<Vec<_>>())
                                        .add_arg(&dummy_wrapped)
                                        .finish()
                                        * FunctionBuilder::new(b.get_symbol())
                                            .add_args(&b.iter().collect::<Vec<_>>())
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
                a
            })
    }

    pub fn convert_to_dot_notation(atom: AtomView) -> Atom {
        let mut old_expr = atom.to_owned().expand();

        loop {
            let mut expr = old_expr
                .replace(vk_parse!("v_(id_,idx_)^n_").unwrap().to_pattern())
                .when(
                    &(Condition::from((vk_symbol!("v_"), symbol_condition()))
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
                    &(Condition::from((vk_symbol!("v1_"), symbol_condition()))
                        & Condition::from((vk_symbol!("v2_"), symbol_condition()))
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
                    &(Condition::from((vk_symbol!("v1_"), symbol_condition()))
                        & Condition::from((vk_symbol!("v2_"), symbol_condition()))
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
        let binary_prec = settings.get_binary_precision();
        expression = expression.map_coefficient(|c| match c {
            CoefficientView::Float(re, im) => {
                let mut re = re.to_float();
                let mut im = im.to_float();
                match settings.precision_for_input_float_rationalization {
                    InputFloatRationalizationPrecision::FullPrecision => {}
                    InputFloatRationalizationPrecision::TargetPrecision => {
                        re.set_prec(binary_prec);
                        im.set_prec(binary_prec);
                    }
                };
                Coefficient::from(Complex::new(re.to_rational(), im.to_rational()))
            }
            _ => c.to_owned(),
        });

        // let mut processed_str =
        //     AtomPrinter::new_with_options(processed.as_view(), PrintOptions::file_no_namespace())
        //         .to_string();
        let mut processed_str = expression.to_canonical_string();
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
                for vecsymbol in [
                    vk_symbol!(LOOP_MOMENTUM_SYMBOL),
                    vk_symbol!(EXTERNAL_MOMENTUM_SYMBOL),
                    vk_symbol!("vec1"),
                    vk_symbol!("vec"),
                ] {
                    let pattern = vk_parse!(format!(
                        "{}(id_,{})",
                        vecsymbol,
                        user_i.to_canonical_string()
                    ))
                    .unwrap()
                    .to_pattern();
                    expression_no_indices = expression_no_indices.replace(pattern).with(
                        vk_parse!(format!(
                            "{}(id_,{})",
                            vecsymbol,
                            (i_index as u64) + FORM_REPLACEMENT_INDEX_SHIFT + 1
                        ))
                        .unwrap(),
                    );
                }
                for metric_symbol in ["g", METRIC_SYMBOL] {
                    let pattern = vk_parse!(format!(
                        "{}(idx1_,{})",
                        metric_symbol,
                        user_i.to_canonical_string()
                    ))
                    .unwrap()
                    .to_pattern();
                    expression_no_indices = expression_no_indices.replace(pattern).with(
                        vk_parse!(format!(
                            "{}(idx1_,{})",
                            metric_symbol,
                            (i_index as u64) + FORM_REPLACEMENT_INDEX_SHIFT + 1
                        ))
                        .unwrap(),
                    );
                    let pattern = vk_parse!(format!(
                        "{}({},idx2_)",
                        metric_symbol,
                        user_i.to_canonical_string()
                    ))
                    .unwrap()
                    .to_pattern();
                    expression_no_indices = expression_no_indices.replace(pattern).with(
                        vk_parse!(format!(
                            "{}({},idx2_)",
                            metric_symbol,
                            (i_index as u64) + FORM_REPLACEMENT_INDEX_SHIFT + 1
                        ))
                        .unwrap(),
                    );
                }
            } else {
                let litteral_form_name = format!("[{}]", user_i.to_canonical_string());
                expression_no_indices = expression_no_indices
                    .replace(user_i.to_pattern())
                    .with(Atom::num(0));
                processed_str =
                    processed_str.replace(&user_i.to_canonical_string(), &litteral_form_name);
                form_header_indices.push(litteral_form_name);
            }
        }
        if substitute_indices {
            processed_str = expression_no_indices.to_canonical_string();
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
        for user_f in user_functions.iter() {
            let litteral_form_name = format!("[{}]", get_full_name(user_f));
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
            let litteral_form_name = format!("[{}]", get_full_name(user_v));
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
        string_replacements.insert("𝑖".into(), "*i_".into());
        //processed_str = processed_str.replace("𝑖", " i_");
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

        Ok((form_header_additions, expression_str, indices))
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
        match vk_parse!(processed_form_str.as_str()) {
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
                    .when(Condition::from((vk_symbol!("v1_"), symbol_condition())))
                    .with(vk_parse!("vec_(idx_)").unwrap().to_pattern());

                // Convert vectors back from pi(j) notation to p(i,j) notation
                for (s, t) in vector_mapping.iter() {
                    processed = processed.replace(s.to_pattern()).with(t.to_pattern());
                }

                for (i_index, user_i) in indices.iter().enumerate() {
                    let mut replacements = vec![];
                    for vec in [LOOP_MOMENTUM_SYMBOL, EXTERNAL_MOMENTUM_SYMBOL] {
                        replacements.push(Replacement::new(
                            vk_parse!(format!(
                                "{}(id_,{})",
                                vec,
                                FORM_REPLACEMENT_INDEX_SHIFT + 1 + (i_index as u64)
                            ))
                            .unwrap()
                            .to_pattern(),
                            vk_parse!(format!("{}(id_,{})", vec, user_i.to_canonical_string()))
                                .unwrap(),
                        ));
                    }
                    replacements.push(Replacement::new(
                        vk_parse!(format!(
                            "{}({},id_)",
                            METRIC_SYMBOL,
                            FORM_REPLACEMENT_INDEX_SHIFT + 1 + (i_index as u64)
                        ))
                        .unwrap()
                        .to_pattern(),
                        vk_parse!(format!(
                            "{}(id_,{})",
                            METRIC_SYMBOL,
                            user_i.to_canonical_string()
                        ))
                        .unwrap(),
                    ));
                    replacements.push(Replacement::new(
                        vk_parse!(format!(
                            "{}(id_,{})",
                            METRIC_SYMBOL,
                            FORM_REPLACEMENT_INDEX_SHIFT + 1 + (i_index as u64)
                        ))
                        .unwrap()
                        .to_pattern(),
                        vk_parse!(format!(
                            "{}({},id_)",
                            METRIC_SYMBOL,
                            user_i.to_canonical_string()
                        ))
                        .unwrap(),
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
        let mut generate_pysecdec_sources = true;
        let tmp_dir = &(if let Some(reused_path_specified) = reused_path.as_ref() {
            let specified_dir = PathBuf::from(reused_path_specified);
            if !specified_dir.exists() {
                if fs::create_dir(specified_dir.clone()).is_err() {
                    return Err(VakintError::PySecDecError(format!(
                        "Could not create user-specified directory to be reused for pySecDec run '{}'.",
                        reused_path_specified.green()
                    )));
                } else {
                    warn!(
                        "User-specified directory '{}' not found, and was instead created and pySecDec sources will be regenerated.",
                        reused_path_specified
                    );
                }
            } else {
                generate_pysecdec_sources = false;
                warn!("{}",format!("User requested to re-use existing directory '{}' for the pysecdec run.\nThis is of course potentially unsafe and should be used for debugging only.\nRemove that directory to start clean.", reused_path_specified).red());
            }
            specified_dir
        } else {
            let tmp_directory = if let Some(temp_path) = temporary_directory {
                &PathBuf::from(temp_path).join("vakint_temp_pysecdec")
            } else {
                &env::temp_dir().join("vakint_temp_pysecdec")
            };
            tmp_directory.clone()
        });

        if generate_pysecdec_sources {
            if tmp_dir.exists() {
                fs::remove_dir_all(tmp_dir)?;
            }
            fs::create_dir(tmp_dir)?;
            for input in input.iter() {
                fs::write(tmp_dir.join(&input.0).to_str().unwrap(), &input.1)?;
            }
        };

        let mut cmd = Command::new(settings.python_exe_path.as_str());
        cmd.arg(input[0].clone().0);
        for opt in options {
            cmd.arg(opt);
        }
        cmd.current_dir(tmp_dir);

        if !clean {
            info!("Running {} with command: {:?}", "PySecDec".green(), cmd);
            info!(
                "You can follow the run live with 'tail -f follow_run.txt' in that temporary directory"
            );
        } else {
            debug!("Running {} with command: {:?}", "PySecDec".green(), cmd);
            debug!(
                "You can follow the run live with 'tail -f follow_run.txt' in that temporary directory"
            );
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
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        if !tmp_dir.join("out.txt").exists() {
            return Err(VakintError::MissingFormOutput(
                "N/A".into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        let result = fs::read_to_string(tmp_dir.join("out.txt"))?;
        if clean {
            fs::remove_dir_all(tmp_dir)?;
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
        let tmp_dir = if let Some(temp_path) = temporary_directory {
            &PathBuf::from(temp_path).join("vakint_temp")
        } else {
            &env::temp_dir().join("vakint_temp")
        };

        if tmp_dir.exists() {
            fs::remove_dir_all(tmp_dir)?;
        }
        fs::create_dir(tmp_dir)?;

        for resource in resources.iter() {
            fs::write(tmp_dir.join(resource), FORM_SRC.get(resource).unwrap())?;
        }
        fs::write(tmp_dir.join(&input.0).to_str().unwrap(), &input.1)?;
        let mut cmd = Command::new(settings.form_exe_path.as_str());
        for opt in options {
            cmd.arg(opt);
        }
        cmd.arg(input.0);
        cmd.current_dir(tmp_dir);
        if !clean {
            info!("Running {} with command: {:?}", "FORM".green(), cmd);
        } else {
            debug!("Running {} with command: {:?}", "FORM".green(), cmd);
        }
        //println!("Running {} with command: {:?}", "FORM".green(), cmd);
        let output = cmd.stderr(Stdio::piped()).stdout(Stdio::piped()).output()?;
        // println!("{}", String::from_utf8_lossy(&output.stdout));
        if !ExitStatus::success(&output.status) {
            return Err(VakintError::FormError(
                String::from_utf8_lossy(&output.stderr).into(),
                String::from_utf8_lossy(&output.stdout).into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }
        if !tmp_dir.join("out.txt").exists() {
            return Err(VakintError::MissingFormOutput(
                String::from_utf8_lossy(&output.stderr).into(),
                format!("{:?}", cmd),
                tmp_dir.to_str().unwrap().into(),
            ));
        }

        let result = fs::read_to_string(tmp_dir.join("out.txt"))?;
        if clean {
            fs::remove_dir_all(tmp_dir)?;
        }
        Ok(result)
    }

    pub fn evaluate(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        self.validate_settings(settings)?;
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.canonicalize(settings, &self.topologies, false)?;
        vakint_expr.tensor_reduce(self, settings)?;
        vakint_expr.evaluate_integral(self, settings)?;
        Ok(vakint_expr.into())
    }

    pub fn to_canonical(
        &self,
        settings: &VakintSettings,
        input: AtomView,
        short_form: bool,
    ) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.canonicalize(settings, &self.topologies, short_form)?;
        Ok(vakint_expr.into())
    }

    pub fn tensor_reduce(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.tensor_reduce(self, settings)?;
        Ok(vakint_expr.into())
    }

    pub fn evaluate_integral(
        &self,
        settings: &VakintSettings,
        input: AtomView,
    ) -> Result<Atom, VakintError> {
        self.validate_settings(settings)?;
        let mut vakint_expr = VakintExpression::try_from(input)?;
        vakint_expr.evaluate_integral(self, settings)?;
        Ok(vakint_expr.into())
    }
}

#[cfg(test)]
mod tests {
    use symbolica::parse_lit;

    use super::*;

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
}
