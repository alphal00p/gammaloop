//! The p master catalogue: values of p-masters beyond the vendored Forcer
//! table (the five-loop GGPPSS masters, arXiv:2104.08272), the map of the
//! RustRed residual terminals onto them, and the guards of the five-loop
//! vacuum master catalogue (`rustred_evaluation::masters`) mirrored one to one
//! (VK5L_PLAN.md section 3.5, package WP5):
//!
//! | vacuum guard (`MasterCatalog`)                   | here                                                   |
//! |--------------------------------------------------|--------------------------------------------------------|
//! | required orders `n_terms - L - 1`, per master    | the finalizer's target; each master is substituted     |
//! | depth planned from the coefficient valuations    | through `known_through`, then `ε^k·Oep(k, PMaster(m))` |
//! | bare master too shallow, spurious poles          | the surviving `Oep` after truncation:                  |
//! | after truncation (`ensure_known_orders`)         | [`crate::truncation_guard::MissingMasterOrders`]       |
//! | missing master: declared without coefficients    | [`PMasterCatalogueError::MissingCoefficients`]         |
//! | missing master: undeclared                        | [`PMasterCatalogueError::UndeclaredMasters`]          |
//! | (p only) residual terminal without a map row     | [`PMasterCatalogueError::UnmappedTerminals`]           |
//! | (p only) no catalogue at all (Retain-only today) | [`PMasterCatalogueError::NoCatalogue`]                 |
//! | finite constant without validated digits         | [`PMasterCatalogueError::Digits`] (`available: None`)  |
//! | achieved digits below the request                | [`PMasterCatalogueError::Digits`] (`available < req.`) |
//! | malformed / non-finite numerical value           | [`PMasterCatalogueError::Invalid`] at load             |
//! | Retain mode only warns about unresolved masters  | the same (`log::warn`), terminals stay `PResidual`     |
//!
//! Every error names the catalogue (`data_version`), the master or terminal,
//! the requested against the available order or digits, and the request.
//!
//! Today no five-loop catalogue ships (coordinator decision Q3: Retain-only
//! first): a `PResidual` terminal reaching `Symbolic` or `Numerical` mode is
//! [`PMasterCatalogueError::NoCatalogue`]. The values plug in later through a
//! catalogue file of schema [`P_MASTER_CATALOGUE_SCHEMA`] (WP3 fills its
//! `terminals`, `tools/p_integral_registry/ggppss_catalogue.py` its `masters`).
//!
//! Catalogue file (JSON, versioned by [`P_MASTER_CATALOGUE_VERSION`], no
//! digests):
//!
//! ```text
//! { "schema": "vakint.p-master-catalogue.v1", "version": 1,
//!   "data_version": "ggppss-5L-v1",            // named in every message
//!   "normalisation": "g-scheme",                // P_MASTER_NORMALISATION
//!   "masters": [ { "id": "ggppss_int(23)",     // the source's name
//!                  "symbol": "ggppss_int_23",  // rendered PMaster(ggppss_int_23)
//!                  "loops": 5, "graph": "<short graph text>",
//!                  "leading_order": -5,
//!                  "known_through": 3,         // optional cross-check
//!                  "unknown_from": 4,          // optional: unknown from ε^4 on
//!                  "significant_digits": 40,   // optional default of the master
//!                  "coefficients": { "-5": {"exact": "1"},
//!                                    "-4": {"numerical": "1.25", "significant_digits": 30} } } ],
//!   "terminals": [ { "node": "T162", "key": [1, 1, ...],
//!                    "status": "mapped",       // or "unmapped", "additional"
//!                    "combination": [ {"master": "ggppss_int(23)", "coefficient": "(d-4)/(2*d-7)"},
//!                                     {"master": null, "coefficient": "PG(0,0)^5"} ],
//!                    "note": "..." } ] }
//! ```
//!
//! A combination coefficient is exact in `d` (replaced by `4 - 2ε`) and may
//! carry the tags `PG(x1,x2)`, `PGamma(m)` and Forcer masters `PMaster(name)`.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use serde::Deserialize;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::domains::float::Float;
use symbolica::id::Replacement;
use symbolica::transcendental::zeta;

use super::{masters, p_gamma_symbol, p_master_symbol, p_residual_symbol, pg_symbol};
use crate::rustred_evaluation::RustRedMasterMode;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{VakintError, get_integer_from_atom};

/// The schema of a p master catalogue file.
pub const P_MASTER_CATALOGUE_SCHEMA: &str = "vakint.p-master-catalogue.v1";

/// The catalogue format version this vakint reads (simple version constant).
pub const P_MASTER_CATALOGUE_VERSION: u32 = 1;

/// The normalisation of every catalogue value: vakint's `C_G` (G-scheme at
/// `Q² = 1`, every loop divided by `εG(1,1)`, the one-loop bubble is `1/ε`),
/// which is Forcer's `Masters.prc` and the GGPPSS normalisation.
pub const P_MASTER_NORMALISATION: &str = "g-scheme";

/// At most this many entries are named in one error message.
const NAMED: usize = 8;


/// `vakint::PMasterValue(master, k)`: the `ε^k` coefficient of a catalogue
/// master that is known numerically only (or embeds a finite constant), kept
/// opaque until its digits are checked after truncation.
fn value_symbol() -> Symbol {
    vk_symbol!("PMasterValue")
}

/// The table epsilon of catalogue texts, `vakint::ep`.
fn table_epsilon() -> Symbol {
    vk_symbol!("ep")
}

/// One requested or missing digits entry of [`PMasterCatalogueError::Digits`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DigitsShortfall {
    /// The master (`id`, as in the catalogue).
    pub master: String,
    /// The Laurent order `ε^k` of the coefficient.
    pub order: i64,
    /// The validated significant digits of the coefficient, if any.
    pub available: Option<u32>,
}

/// A guard of the p master catalogue refused the request.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum PMasterCatalogueError {
    /// The catalogue file is unusable (schema, version, normalisation, a
    /// malformed or non-finite value, an inconsistent depth).
    #[error("p master catalogue {source_name}: {detail}")]
    Invalid { source_name: String, detail: String },
    /// Residual terminals reached an expanding mode while no catalogue is
    /// loaded (decision Q3: five-loop p-integrals are Retain-only until a
    /// catalogue ships).
    #[error(
        "no p master catalogue ({P_MASTER_CATALOGUE_SCHEMA}) is loaded: the reduction ended on the \
         residual terminal(s) {}, which master_mode={mode} must evaluate; master_mode=retain returns \
         them as PResidual symbols (request: {request})",
        named(terminals)
    )]
    NoCatalogue {
        mode: String,
        terminals: Vec<String>,
        request: String,
    },
    /// Residual terminals without a master map row (WP3), or additional
    /// terminals of boxed nodes that have no value yet (decision Q4).
    #[error(
        "p master data {data_version}: master_mode={mode} needs values for residual terminal(s) \
         without a master map row:\n  - {}\n(request: {request})",
        terminals.join("\n  - ")
    )]
    UnmappedTerminals {
        data_version: String,
        mode: String,
        /// `PResidual(..): <why>`.
        terminals: Vec<String>,
        request: String,
    },
    /// Masters the catalogue declares without any computed coefficient.
    #[error(
        "Unresolved p masters (data {data_version}):\n  - {}\n(request: {request})",
        masters.join("\n  - ")
    )]
    MissingCoefficients {
        data_version: String,
        /// `<id>: declared terminal has no computed coefficients; graph <graph>`.
        masters: Vec<String>,
        request: String,
    },
    /// `PMaster(name)` that neither the catalogue nor the Forcer table declares.
    #[error(
        "p master data {data_version}: undeclared p master(s) {}: no value and no declaration \
         (request: {request})",
        named(masters)
    )]
    UndeclaredMasters {
        data_version: String,
        masters: Vec<String>,
        request: String,
    },
    /// Coefficients known numerically only, in `Symbolic` mode.
    #[error(
        "p master data {data_version}: master_mode=symbolic needs exact coefficients, but {} \
         {} known numerically only; use master_mode=numerical (request: {request})",
        named(coefficients),
        if coefficients.len() == 1 { "is" } else { "are" }
    )]
    NumericalOnly {
        data_version: String,
        /// `<id> coefficient epsilon^k`.
        coefficients: Vec<String>,
        request: String,
    },
    /// Coefficients without validated digits or with fewer digits than
    /// `run_time_decimal_precision`.
    #[error(
        "p master data {data_version}: the result at {requested} digits needs master coefficients \
         that cannot support it:{} (request: {request})",
        shortfalls.iter().map(|s| match s.available {
            None => format!("\n  - master {} coefficient epsilon^{}: a finite constant without validated significant digits", s.master, s.order),
            Some(a) => format!("\n  - master {} coefficient epsilon^{}: {a} validated digits < {requested} requested", s.master, s.order),
        }).collect::<String>()
    )]
    Digits {
        data_version: String,
        requested: u32,
        shortfalls: Vec<DigitsShortfall>,
        request: String,
    },
}

fn named(items: &[String]) -> String {
    let mut text = items
        .iter()
        .take(NAMED)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");
    if items.len() > NAMED {
        text += &format!(" and {} more", items.len() - NAMED);
    }
    text
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CatalogueFile {
    schema: String,
    version: u32,
    data_version: String,
    normalisation: String,
    #[serde(default)]
    #[allow(dead_code)]
    description: Option<String>,
    masters: Vec<MasterFile>,
    #[serde(default)]
    terminals: Vec<TerminalFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MasterFile {
    id: String,
    symbol: String,
    loops: usize,
    graph: String,
    leading_order: i64,
    #[serde(default)]
    known_through: Option<i64>,
    #[serde(default)]
    unknown_from: Option<i64>,
    #[serde(default)]
    significant_digits: Option<u32>,
    #[serde(default)]
    coefficients: BTreeMap<String, CoefficientFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CoefficientFile {
    #[serde(default)]
    exact: Option<String>,
    #[serde(default)]
    numerical: Option<String>,
    #[serde(default)]
    significant_digits: Option<u32>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum TerminalStatus {
    /// A combination of masters (and tags) exact in `d`.
    Mapped,
    /// No row yet (WP3 has not reduced it onto masters).
    Unmapped,
    /// An additional terminal of a boxed node (decision Q4): a new master
    /// whose value is not known yet.
    Additional,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TerminalFile {
    node: String,
    key: Vec<i64>,
    status: TerminalStatus,
    #[serde(default)]
    combination: Vec<TermFile>,
    #[serde(default)]
    note: Option<String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TermFile {
    #[serde(default)]
    master: Option<String>,
    coefficient: String,
}

/// One Laurent coefficient of a catalogue master.
#[derive(Clone, Debug)]
pub(crate) struct PCoefficient {
    /// The exact value (constants `zeta`, `MZV`, `π`); may embed finite
    /// decimals, which then need validated digits.
    exact: Option<Atom>,
    /// A decimal literal.
    numerical: Option<String>,
    significant_digits: Option<u32>,
}

impl PCoefficient {
    /// Exact and free of finite decimals: substituted before the expansion.
    fn exact_value(&self) -> Option<&Atom> {
        self.exact
            .as_ref()
            .filter(|exact| finite_precision(exact).is_none())
    }
}

/// One master of the catalogue.
#[derive(Clone, Debug)]
pub(crate) struct PCatalogueMaster {
    pub(crate) id: String,
    pub(crate) loops: usize,
    pub(crate) graph: String,
    pub(crate) leading_order: i64,
    /// Every order `leading_order..=known_through` has a coefficient
    /// (`leading_order - 1` if none).
    pub(crate) known_through: i64,
    coefficients: BTreeMap<i64, PCoefficient>,
}

#[derive(Clone, Debug)]
enum Terminal {
    /// `Σ coefficient(ε)·master` with the table epsilon.
    Mapped(Atom),
    Unmapped(Option<String>),
    Additional(Option<String>),
}

/// A loaded p master catalogue (see the module documentation).
#[derive(Clone, Debug)]
pub struct PMasterCatalogue {
    data_version: String,
    /// By symbol name (the argument of `PMaster`).
    masters: BTreeMap<String, PCatalogueMaster>,
    /// By `(node, key)`.
    terminals: BTreeMap<(String, Vec<i64>), Terminal>,
}

/// The decimals' precision of `expression` in digits, if it embeds any.
fn finite_precision(expression: &Atom) -> Option<u32> {
    let mut bits = None;
    expression.visitor(&mut |part| {
        if let AtomView::Num(number) = part
            && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
        {
            for component in [real, imaginary] {
                let value = component.to_float();
                if !value.as_raw().is_zero() {
                    bits = Some(bits.map_or(value.prec(), |old: u32| old.min(value.prec())));
                }
            }
        }
        true
    });
    bits.map(|value| (f64::from(value) * std::f64::consts::LOG10_2).floor() as u32)
}

/// Whether `symbol` is a plain identifier (it becomes `vakint::symbol`).
fn plain_identifier(symbol: &str) -> bool {
    let mut chars = symbol.chars();
    chars.next().is_some_and(|c| c.is_ascii_alphabetic())
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// The integer power of `epsilon` in an expanded product term.
#[cfg(test)]
fn epsilon_power(term: AtomView, epsilon: Symbol) -> Option<i64> {
    match term {
        AtomView::Var(variable) if variable.get_symbol() == epsilon => Some(1),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            match base {
                AtomView::Var(variable) if variable.get_symbol() == epsilon => {
                    get_integer_from_atom(exponent)
                }
                _ if !base.contains_symbol(epsilon) => Some(0),
                _ => None,
            }
        }
        AtomView::Mul(product) => product
            .iter()
            .map(|factor| epsilon_power(factor, epsilon))
            .sum(),
        _ if !term.contains_symbol(epsilon) => Some(0),
        _ => None,
    }
}

/// `PMaster(symbol)`.
fn p_master(symbol: &str) -> Atom {
    p_master_symbol().call(&[Atom::var(vk_symbol!(symbol))][..])
}

impl PMasterCatalogue {
    /// Parse and validate a catalogue file; `source_name` names it in errors.
    pub fn from_json(source_name: &str, text: &str) -> Result<Self, PMasterCatalogueError> {
        crate::Vakint::initialize_vakint_symbols();
        let invalid = |detail: String| PMasterCatalogueError::Invalid {
            source_name: source_name.to_owned(),
            detail,
        };
        let file: CatalogueFile = serde_json::from_str(text)
            .map_err(|error| invalid(format!("not a catalogue: {error}")))?;
        if file.schema != P_MASTER_CATALOGUE_SCHEMA {
            return Err(invalid(format!(
                "schema {}, this vakint reads {P_MASTER_CATALOGUE_SCHEMA}",
                file.schema
            )));
        }
        if file.version != P_MASTER_CATALOGUE_VERSION {
            return Err(invalid(format!(
                "catalogue version {}, this vakint reads version {P_MASTER_CATALOGUE_VERSION}; \
                 regenerate it (tools/p_integral_registry/ggppss_catalogue.py)",
                file.version
            )));
        }
        if file.normalisation != P_MASTER_NORMALISATION {
            return Err(invalid(format!(
                "normalisation {:?}, vakint's p masters are {P_MASTER_NORMALISATION:?}",
                file.normalisation
            )));
        }
        if file.data_version.trim().is_empty() {
            return Err(invalid("an empty data_version".into()));
        }
        let forcer = masters::p_master_names().map_err(|error| invalid(error.to_string()))?;
        let mut catalogue = PMasterCatalogue {
            data_version: file.data_version.clone(),
            masters: BTreeMap::new(),
            terminals: BTreeMap::new(),
        };
        let mut ids = BTreeMap::new();
        for master in file.masters {
            let name = format!("master {}", master.id);
            if !plain_identifier(&master.symbol) {
                return Err(invalid(format!(
                    "{name}: symbol {:?} is not a plain identifier",
                    master.symbol
                )));
            }
            if forcer.contains(&master.symbol.as_str()) {
                return Err(invalid(format!(
                    "{name}: symbol {} is a vendored Forcer master",
                    master.symbol
                )));
            }
            if ids
                .insert(master.id.clone(), master.symbol.clone())
                .is_some()
            {
                return Err(invalid(format!("{name} is defined twice")));
            }
            let mut coefficients = BTreeMap::new();
            for (order, coefficient) in master.coefficients {
                let order: i64 = order
                    .parse()
                    .map_err(|_| invalid(format!("{name}: order {order:?} is not an integer")))?;
                let at = format!("{name} coefficient epsilon^{order}");
                if order < master.leading_order {
                    return Err(invalid(format!(
                        "{at} is below its leading order {}",
                        master.leading_order
                    )));
                }
                if master.unknown_from.is_some_and(|unknown| order >= unknown) {
                    return Err(invalid(format!(
                        "{at} lies in its unknown tail (unknown from epsilon^{})",
                        master.unknown_from.unwrap()
                    )));
                }
                let exact = coefficient
                    .exact
                    .as_deref()
                    .map(|text| {
                        parse_constant(text).map_err(|detail| invalid(format!("{at}: {detail}")))
                    })
                    .transpose()?;
                if let Some(text) = &coefficient.numerical {
                    let number = Float::parse(text, Some(64)).map_err(|error| {
                        invalid(format!("{at} has an invalid numerical value: {error}"))
                    })?;
                    if !number.as_raw().is_finite() {
                        return Err(invalid(format!("{at} has a non-finite numerical value")));
                    }
                }
                if exact.is_none() && coefficient.numerical.is_none() {
                    return Err(invalid(format!(
                        "{at} has neither an exact nor a numerical value"
                    )));
                }
                coefficients.insert(
                    order,
                    PCoefficient {
                        exact,
                        numerical: coefficient.numerical,
                        significant_digits: coefficient
                            .significant_digits
                            .or(master.significant_digits),
                    },
                );
            }
            let mut known_through = master.leading_order - 1;
            while coefficients.contains_key(&(known_through + 1)) {
                known_through += 1;
            }
            if let Some((order, _)) = coefficients.range(known_through + 1..).next() {
                return Err(invalid(format!(
                    "{name}: coefficient epsilon^{order} follows the missing epsilon^{}",
                    known_through + 1
                )));
            }
            if let Some(declared) = master.known_through
                && declared != known_through
            {
                return Err(invalid(format!(
                    "{name}: known_through {declared}, but its coefficients reach epsilon^{known_through}"
                )));
            }
            if let Some(unknown) = master.unknown_from
                && unknown != known_through + 1
                && !coefficients.is_empty()
            {
                return Err(invalid(format!(
                    "{name}: unknown from epsilon^{unknown}, but its coefficients reach epsilon^{known_through}"
                )));
            }
            if catalogue.masters.contains_key(&master.symbol) {
                return Err(invalid(format!(
                    "{name}: symbol {} is used twice",
                    master.symbol
                )));
            }
            catalogue.masters.insert(
                master.symbol.clone(),
                PCatalogueMaster {
                    id: master.id,
                    loops: master.loops,
                    graph: master.graph,
                    leading_order: master.leading_order,
                    known_through,
                    coefficients,
                },
            );
        }
        let d = Atom::var(vk_symbol!("d"));
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(table_epsilon());
        for terminal in file.terminals {
            let label = residual_label(&terminal.node, &terminal.key);
            if !plain_identifier(&terminal.node) {
                return Err(invalid(format!(
                    "terminal {label}: node is not a plain identifier"
                )));
            }
            let entry = match terminal.status {
                TerminalStatus::Mapped => {
                    if terminal.combination.is_empty() {
                        return Err(invalid(format!(
                            "terminal {label}: a mapped row without terms"
                        )));
                    }
                    let mut value = Atom::Zero;
                    for term in &terminal.combination {
                        let coefficient = vk_parse!(&term.coefficient)
                            .map_err(|error| {
                                invalid(format!(
                                    "terminal {label}: coefficient {:?}: {error}",
                                    term.coefficient
                                ))
                            })?
                            .replace(d.to_pattern())
                            .with(dimension.to_pattern());
                        check_coefficient(&coefficient)
                            .map_err(|detail| invalid(format!("terminal {label}: {detail}")))?;
                        let basis = match &term.master {
                            None => Atom::num(1),
                            Some(master) => {
                                let symbol = ids
                                    .get(master)
                                    .cloned()
                                    .or_else(|| {
                                        catalogue
                                            .masters
                                            .contains_key(master)
                                            .then(|| master.clone())
                                    })
                                    .or_else(|| {
                                        forcer.contains(&master.as_str()).then(|| master.clone())
                                    })
                                    .ok_or_else(|| {
                                        invalid(format!(
                                            "terminal {label}: undeclared master {master}"
                                        ))
                                    })?;
                                p_master(&symbol)
                            }
                        };
                        value += coefficient * basis;
                    }
                    Terminal::Mapped(value)
                }
                TerminalStatus::Unmapped => Terminal::Unmapped(terminal.note),
                TerminalStatus::Additional => Terminal::Additional(terminal.note),
            };
            if catalogue
                .terminals
                .insert((terminal.node, terminal.key), entry)
                .is_some()
            {
                return Err(invalid(format!("terminal {label} is listed twice")));
            }
        }
        Ok(catalogue)
    }

    /// Read and validate a catalogue file.
    pub fn load(path: &std::path::Path) -> Result<Self, PMasterCatalogueError> {
        let text =
            std::fs::read_to_string(path).map_err(|error| PMasterCatalogueError::Invalid {
                source_name: path.display().to_string(),
                detail: error.to_string(),
            })?;
        Self::from_json(&path.display().to_string(), &text)
    }

    /// The vendored Forcer table (`Masters.prc`) as a catalogue, under the
    /// same names: the one- to four-loop masters exercise the catalogue path
    /// with real data (the tests compare it with the Forcer path). Test only:
    /// a catalogue file may not reuse Forcer names.
    #[cfg(test)]
    pub(crate) fn from_forcer_table() -> Result<Self, VakintError> {
        let epsilon = table_epsilon();
        let mut catalogue = PMasterCatalogue {
            data_version: "forcer-masters-prc".into(),
            masters: BTreeMap::new(),
            terminals: BTreeMap::new(),
        };
        for (name, expansion) in masters::p_master_table()? {
            let mut coefficients = BTreeMap::new();
            let mut unknown_from = None;
            let expanded = expansion.expand();
            let terms: Vec<_> = match expanded.as_view() {
                AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect(),
                _ => vec![expanded.clone()],
            };
            for term in terms {
                let power = epsilon_power(term.as_view(), epsilon).ok_or_else(|| {
                    VakintError::EvaluationError(format!("Masters.prc {name}: term {term}"))
                })?;
                let coefficient = (term / Atom::var(epsilon).pow(Atom::num(power))).expand();
                if crate::truncation_guard::has_unknown_tail(coefficient.as_view()) {
                    unknown_from = Some(power);
                    continue;
                }
                *coefficients.entry(power).or_insert(Atom::Zero) += coefficient;
            }
            let unknown_from = unknown_from.ok_or_else(|| {
                VakintError::EvaluationError(format!("Masters.prc {name}: no unknown tail"))
            })?;
            let leading_order = *coefficients.keys().next().unwrap_or(&unknown_from);
            catalogue.masters.insert(
                name.clone(),
                PCatalogueMaster {
                    id: name.clone(),
                    loops: if ["no", "t105"].contains(&name.as_str()) {
                        3
                    } else {
                        4
                    },
                    graph: format!("Forcer Master({name})"),
                    leading_order,
                    known_through: unknown_from - 1,
                    coefficients: (leading_order..unknown_from)
                        .map(|power| {
                            (
                                power,
                                PCoefficient {
                                    exact: Some(
                                        coefficients.remove(&power).unwrap_or_default().expand(),
                                    ),
                                    numerical: None,
                                    significant_digits: None,
                                },
                            )
                        })
                        .collect(),
                },
            );
        }
        Ok(catalogue)
    }

    /// The catalogue's data version (named in every message).
    pub fn data_version(&self) -> &str {
        &self.data_version
    }

    /// The number of masters and of terminal rows.
    pub fn counts(&self) -> (usize, usize) {
        (self.masters.len(), self.terminals.len())
    }

    #[cfg(test)]
    pub(crate) fn master(&self, symbol: &str) -> Option<&PCatalogueMaster> {
        self.masters.get(symbol)
    }
}

/// `PResidual(node,k..)` as text.
fn residual_label(node: &str, key: &[i64]) -> String {
    format!(
        "PResidual({node},{})",
        key.iter().map(i64::to_string).collect::<Vec<_>>().join(",")
    )
}

/// A coefficient text: rational numbers and the constants `zeta(n)`,
/// `MZV(..)` and `π`; anything else (an undetermined symbol, `ep`) is refused.
fn parse_constant(text: &str) -> Result<Atom, String> {
    let parsed = vk_parse!(text).map_err(|error| error.to_string())?;
    let parsed = parsed
        .replace(Atom::var(vk_symbol!("pi")).to_pattern())
        .with(Atom::var(Symbol::PI).to_pattern())
        .replace(
            vk_symbol!("zeta")
                .call(&[Atom::var(crate::symbols::S.x_)][..])
                .to_pattern(),
        )
        .with(
            zeta()
                .call(&[Atom::var(crate::symbols::S.x_)][..])
                .to_pattern(),
        );
    let mut unknown = None;
    parsed.visitor(&mut |part| match part {
        AtomView::Var(variable) if variable.get_symbol() != Symbol::PI => {
            unknown = Some(part.to_owned());
            false
        }
        AtomView::Fun(function) if [zeta(), vk_symbol!("MZV")].contains(&function.get_symbol()) => {
            false
        }
        AtomView::Fun(_) => {
            unknown = Some(part.to_owned());
            false
        }
        _ => true,
    });
    match unknown {
        Some(part) => Err(format!(
            "unknown constant {part} (an undetermined coefficient belongs to the unknown tail)"
        )),
        None => Ok(parsed.expand()),
    }
}

/// A combination coefficient: rational in the table epsilon, times tags.
fn check_coefficient(coefficient: &Atom) -> Result<(), String> {
    let mut unknown = None;
    coefficient.visitor(&mut |part| match part {
        AtomView::Var(variable)
            if ![table_epsilon(), Symbol::PI].contains(&variable.get_symbol()) =>
        {
            unknown = Some(part.to_owned());
            false
        }
        AtomView::Fun(function)
            if [pg_symbol(), p_gamma_symbol(), p_master_symbol()]
                .contains(&function.get_symbol()) =>
        {
            false
        }
        AtomView::Fun(_) => {
            unknown = Some(part.to_owned());
            false
        }
        _ => true,
    });
    match unknown {
        Some(part) => Err(format!("coefficient symbol {part} is neither d nor a tag")),
        None => Ok(()),
    }
}

/// The p masters a finalization may substitute: the Forcer table, plus an
/// optional catalogue, and the request it serves (for messages).
///
/// [`PMasterSource::forcer`] is the 1L-4L behaviour unchanged.
#[derive(Clone, Copy, Debug)]
pub struct PMasterSource<'a> {
    catalogue: Option<&'a PMasterCatalogue>,
    request: &'a str,
}

impl<'a> PMasterSource<'a> {
    /// The vendored Forcer table only (no catalogue).
    pub fn forcer() -> Self {
        Self {
            catalogue: None,
            request: "",
        }
    }

    /// The Forcer table and `catalogue` (if any) for `request`.
    pub fn new(catalogue: Option<&'a PMasterCatalogue>, request: &'a str) -> Self {
        Self { catalogue, request }
    }

    fn request(&self) -> String {
        if self.request.is_empty() {
            "unnamed".into()
        } else {
            self.request.to_owned()
        }
    }

    /// The catalogue's master names: [`super::compare::check_symbols_with`].
    pub(crate) fn extra_names(&self) -> Vec<&'a str> {
        self.catalogue
            .map(|c| c.masters.keys().map(String::as_str).collect())
            .unwrap_or_default()
    }

    /// The backend label of a missing-order report.
    pub(crate) fn backend(&self, backend: &str) -> String {
        match self.catalogue {
            None => backend.to_owned(),
            Some(catalogue) => format!(
                "{backend}, p master data {}; request: {}",
                catalogue.data_version,
                self.request()
            ),
        }
    }

    /// A master of a missing-order report: `PMaster(name)`, then its catalogue
    /// identity.
    pub(crate) fn describe(&self, master: &Atom) -> String {
        let plain = crate::truncation_guard::plain_name(master);
        let entry = masters::p_master_name(master.as_view())
            .and_then(|name| self.catalogue.and_then(|c| c.masters.get(&name)));
        match entry {
            Some(entry) => format!(
                "{plain} [{}, {}-loop, graph {}, known ε^{}..ε^{}]",
                entry.id, entry.loops, entry.graph, entry.leading_order, entry.known_through
            ),
            None => plain,
        }
    }

    /// The expansion of `PMaster(name)` in `epsilon`: the catalogue's
    /// coefficients through `known_through` and `ε^k·Oep(k, PMaster(name))`,
    /// or the Forcer table's.
    pub(crate) fn expansion(&self, name: &str, epsilon: Symbol) -> Result<Atom, VakintError> {
        let Some(master) = self.catalogue.and_then(|c| c.masters.get(name)) else {
            return masters::p_master_expansion(name, epsilon);
        };
        let e = Atom::var(epsilon);
        let mut expansion = crate::truncation_guard::unknown_tail(
            epsilon,
            master.known_through + 1,
            &p_master(name),
        );
        for (&order, coefficient) in &master.coefficients {
            let value = match coefficient.exact_value() {
                Some(exact) => exact.clone(),
                None => value_symbol().call(&[Atom::var(vk_symbol!(name)), Atom::num(order)][..]),
            };
            expansion += value * e.pow(Atom::num(order));
        }
        Ok(expansion)
    }

    /// Replace every `PResidual` by its master combination. Without a row
    /// (or without a catalogue), `Retain` keeps it with a warning and the
    /// other modes refuse.
    pub(crate) fn resolve_terminals(
        &self,
        reduced: &Atom,
        mode: RustRedMasterMode,
        epsilon: Symbol,
    ) -> Result<Atom, VakintError> {
        let mut residuals = BTreeSet::new();
        reduced.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == p_residual_symbol()
            {
                residuals.insert(part.to_owned());
                return false;
            }
            true
        });
        if residuals.is_empty() {
            return Ok(reduced.clone());
        }
        let Some(catalogue) = self.catalogue else {
            if mode == RustRedMasterMode::Retain {
                return Ok(reduced.clone());
            }
            return Err(PMasterCatalogueError::NoCatalogue {
                mode: mode.to_string(),
                terminals: residuals
                    .iter()
                    .map(|r| crate::truncation_guard::plain_name(r))
                    .collect(),
                request: self.request(),
            }
            .into());
        };
        let mut replacements = Vec::new();
        let mut missing = Vec::new();
        for residual in residuals {
            let label = crate::truncation_guard::plain_name(&residual);
            let row = residual_key(&residual).and_then(|key| catalogue.terminals.get(&key));
            match row {
                Some(Terminal::Mapped(value)) => replacements.push(Replacement::new(
                    residual.to_pattern(),
                    value
                        .replace(Atom::var(table_epsilon()).to_pattern())
                        .with(Atom::var(epsilon).to_pattern()),
                )),
                Some(Terminal::Unmapped(note)) => missing.push(format!(
                    "{label}: no master map row yet{}",
                    note.as_ref().map(|n| format!(" ({n})")).unwrap_or_default()
                )),
                Some(Terminal::Additional(note)) => missing.push(format!(
                    "{label}: additional terminal of a boxed node, its value is not known yet{}",
                    note.as_ref().map(|n| format!(" ({n})")).unwrap_or_default()
                )),
                None => missing.push(format!("{label}: not listed in the master map")),
            }
        }
        let resolved = if replacements.is_empty() {
            reduced.clone()
        } else {
            reduced.replace_multiple(&replacements)
        };
        if missing.is_empty() {
            return Ok(resolved);
        }
        if mode == RustRedMasterMode::Retain {
            log::warn!(
                "Unresolved p terminals (data {}):\n  - {}",
                catalogue.data_version,
                missing.join("\n  - ")
            );
            return Ok(resolved);
        }
        Err(PMasterCatalogueError::UnmappedTerminals {
            data_version: catalogue.data_version.clone(),
            mode: mode.to_string(),
            terminals: missing,
            request: self.request(),
        }
        .into())
    }

    /// The missing-master guard: every `PMaster(name)` is a Forcer master or
    /// a catalogue master with coefficients. `Retain` only warns about
    /// masters declared without coefficients.
    pub(crate) fn require_coefficients(
        &self,
        reduced: &Atom,
        mode: RustRedMasterMode,
    ) -> Result<(), VakintError> {
        let Some(catalogue) = self.catalogue else {
            // The Forcer path: check_symbols refuses unknown names.
            return Ok(());
        };
        let forcer = masters::p_master_names()?;
        let mut undeclared = BTreeSet::new();
        let mut empty = BTreeSet::new();
        reduced.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == p_master_symbol()
            {
                match masters::p_master_name(part) {
                    Some(name) if forcer.contains(&name.as_str()) => {}
                    Some(name) => match catalogue.masters.get(&name) {
                        Some(master) if master.coefficients.is_empty() => {
                            empty.insert(format!(
                                "{}: declared terminal has no computed coefficients; graph {}",
                                master.id, master.graph
                            ));
                        }
                        Some(_) => {}
                        None => {
                            undeclared
                                .insert(crate::truncation_guard::plain_name(&part.to_owned()));
                        }
                    },
                    None => {
                        undeclared.insert(crate::truncation_guard::plain_name(&part.to_owned()));
                    }
                }
                return false;
            }
            true
        });
        if !undeclared.is_empty() {
            return Err(PMasterCatalogueError::UndeclaredMasters {
                data_version: catalogue.data_version.clone(),
                masters: undeclared.into_iter().collect(),
                request: self.request(),
            }
            .into());
        }
        if empty.is_empty() {
            return Ok(());
        }
        let masters: Vec<_> = empty.into_iter().collect();
        if mode == RustRedMasterMode::Retain {
            log::warn!(
                "Unresolved p masters (data {}):\n  - {}",
                catalogue.data_version,
                masters.join("\n  - ")
            );
            return Ok(());
        }
        Err(PMasterCatalogueError::MissingCoefficients {
            data_version: catalogue.data_version.clone(),
            masters,
            request: self.request(),
        }
        .into())
    }

    /// After truncation and the `Oep` guard: the coefficients that are known
    /// numerically only (or embed decimals). `Symbolic` refuses them;
    /// `Numerical` checks their validated digits against
    /// `run_time_decimal_precision` and substitutes them.
    pub(crate) fn resolve_values(
        &self,
        expression: Atom,
        mode: RustRedMasterMode,
        requested: u32,
    ) -> Result<Atom, VakintError> {
        let Some(catalogue) = self.catalogue else {
            return Ok(expression);
        };
        let mut references = BTreeSet::new();
        expression.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == value_symbol()
            {
                references.insert(part.to_owned());
                return false;
            }
            true
        });
        if references.is_empty() {
            return Ok(expression);
        }
        let mut entries = Vec::new();
        for reference in &references {
            let AtomView::Fun(function) = reference.as_view() else {
                unreachable!()
            };
            let arguments: Vec<_> = function.iter().collect();
            let found = match arguments.as_slice() {
                [AtomView::Var(name), order] => get_integer_from_atom(*order).and_then(|order| {
                    let master = catalogue
                        .masters
                        .get(name.get_symbol().get_stripped_name())?;
                    Some((master, order, master.coefficients.get(&order)?))
                }),
                _ => None,
            };
            let (master, order, coefficient) = found.ok_or_else(|| {
                VakintError::EvaluationError(format!(
                    "p master data {}: {reference} is no catalogue coefficient",
                    catalogue.data_version
                ))
            })?;
            entries.push((reference, master, order, coefficient));
        }
        if mode != RustRedMasterMode::Numerical {
            return Err(PMasterCatalogueError::NumericalOnly {
                data_version: catalogue.data_version.clone(),
                coefficients: entries
                    .iter()
                    .map(|(_, master, order, _)| {
                        format!("{} coefficient epsilon^{order}", master.id)
                    })
                    .collect(),
                request: self.request(),
            }
            .into());
        }
        let mut shortfalls = Vec::new();
        let mut replacements = Vec::new();
        for (reference, master, order, coefficient) in entries {
            let available = coefficient.significant_digits.filter(|digits| *digits > 0);
            if available.is_none_or(|digits| digits < requested) {
                shortfalls.push(DigitsShortfall {
                    master: master.id.clone(),
                    order,
                    available,
                });
                continue;
            }
            let digits = available.unwrap();
            let value = match (&coefficient.exact, &coefficient.numerical) {
                (Some(exact), _) => exact.clone(),
                (None, Some(text)) => {
                    let bits = crate::precision::decimal_bits(digits).saturating_add(1);
                    Atom::num(
                        Float::parse(text, Some(bits)).map_err(VakintError::NumericalPrecision)?,
                    )
                }
                (None, None) => unreachable!("validated at load"),
            };
            replacements.push(Replacement::new(reference.to_pattern(), value));
        }
        if !shortfalls.is_empty() {
            return Err(PMasterCatalogueError::Digits {
                data_version: catalogue.data_version.clone(),
                requested,
                shortfalls,
                request: self.request(),
            }
            .into());
        }
        Ok(expression.replace_multiple(&replacements).expand())
    }
}

/// `(node, key)` of `PResidual(node, k..)`.
fn residual_key(residual: &Atom) -> Option<(String, Vec<i64>)> {
    let AtomView::Fun(function) = residual.as_view() else {
        return None;
    };
    let mut arguments = function.iter();
    let AtomView::Var(node) = arguments.next()? else {
        return None;
    };
    let key = arguments
        .map(get_integer_from_atom)
        .collect::<Option<Vec<_>>>()?;
    Some((node.get_symbol().get_stripped_name().to_owned(), key))
}

impl fmt::Display for PMasterCatalogue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "p master catalogue {} ({} masters, {} terminal rows)",
            self.data_version,
            self.masters.len(),
            self.terminals.len()
        )
    }
}

#[cfg(test)]
#[path = "master_catalogue_tests.rs"]
mod tests;
