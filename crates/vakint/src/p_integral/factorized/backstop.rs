//! THE RESULT CHECK of a marked value (the last line behind validation): after
//! the evaluation, before an entry returns, a value that is no value is an
//! internal failure (`Evaluation`, detail prefix `result check: `):
//! - (a) vakint's INTERNAL symbols: the placeholders, parameters, pattern
//!   wildcards and input notation vakint and the factorised entries create
//!   and substitute away ([`internal`]: exact names, families of names made
//!   with format strings, every name ending in `_`, and the input-notation
//!   heads `p`, `k`, `dot`, `g`, `topo`, `prop`, .. where the entry returns the
//!   caller's notation). EVERY OTHER `vakint::` symbol is vakint's output
//!   vocabulary and is allowed (HR, GPLReal, cl2, T1ep, D6, om1, the RustRed
//!   tags, the Retain master tags, vakint's 𝑖, ..). vakint's DATA ATOMS
//!   ([`super::VAKINT_DATA_ATOMS`]: `UNKNOWN(master, k)`, `MZV(..)`) are
//!   output too, except `UNKNOWN` when the settings do not allow unknown
//!   master coefficients and no input carries it. A test (the tripwire) reads
//!   the arguments of the listed symbol constructors and parse macros in the
//!   crate's non-test code: every vakint name a literal argument creates must
//!   be internal, a data atom, in the display table or on the reviewed output
//!   list (`REVIEWED_OUTPUT`, test code), and every argument it cannot resolve
//!   to a literal must be on its reviewed site list;
//! - (b) a symbol of ANY OTHER namespace that is not among the INPUTS (every
//!   atom of the record, the prefactor argument, the settings' epsilon and
//!   scale symbols), Symbolica's own (`symbolica::`) or the caller notation's
//!   structural heads (`spenso::g`, `dot`, `mink`, `bis`, `gamma`,
//!   `gammalooprs::Q`);
//! - (c) the dimension other than as the dimension argument of a
//!   representation slot;
//! - (d) a non-finite or indeterminate number, or a power of a literal zero
//!   base whose exponent is not a positive rational number.
//!
//! What it cannot catch: a wrong value made of legitimate symbols (a wrong
//! coefficient, a missing order); an internal symbol vakint creates at run
//! time through a reviewed site (a name from data, settings, the caller or a
//! format placeholder) or through a constructor the tripwire does not list;
//! an input symbol that should have been substituted (it is an input); a zero
//! that is no number (a zero base written as an unexpanded expression).

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;

use super::{FactorizedError, FactorizedKind, FactorizedNumerator, Result, names};
use crate::VakintSettings;
use crate::utils::vakint_macros::vk_symbol;

/// vakint's internal names (without the `vakint::` namespace), exact.
const INTERNAL_EXACT: [&str; 27] = [
    // the factorised entries
    "closure_lambda",
    "slash",
    // vakint's symbol table and its input / analysis placeholders
    "partial_numerical_input",
    "rustred_input",
    "normalization_input",
    "momentum_match_parameter",
    "coordinate_limit_spectator",
    "coordinate_metadata",
    "VakintLambdaScalingAnalysis",
    "RustRedMaterializeSpectator",
    "pySecDec_integral",
    "pysecdec_result",
    "Integral",
    "LoopIntegral",
    "RelabeledLoopIntegral",
    "SeriesCoefficient",
    // precision.rs staging placeholders (substituted inside stage_table)
    "numerical_coefficient",
    "numerical_expression",
    "evaluation_input",
    // precision.rs: the error-aggregation formula's inputs (evaluated with
    // an input map in the same function)
    "aggregate_error_real_a",
    "aggregate_error_real_b",
    "aggregate_error_imaginary_a",
    "aggregate_error_imaginary_b",
    // topologies.rs / lib.rs: keys of the canonical-expression substitution
    // map, never placed in an expression
    "integral",
    "n_props",
    // the alphaLoop / FMFT / MATAD mass stand-in and normalisation log,
    // substituted back before each backend returns (lib.rs, fmft/finalize.rs,
    // matad.rs); no backend and option returns them (the result-check runs)
    "mUV",
    "logmUVmu",
];

/// Families of internal names made with format strings: a prefix and a
/// suffix around the formatted part (the code's pattern in the comment).
const INTERNAL_FAMILIES: [(&str, &str); 4] = [
    ("closure_tensor_", ""), // closure_tensor_{n} (the closure entry's slot placeholders)
    ("graph_entry_branch_", ""), // graph_entry_branch_{i} (the graph entry's colour branches)
    ("input_coefficient_placeholder_", ""), // input_coefficient_placeholder_{index} (lib.rs)
    ("coordinate_metadata_", ""), // coordinate_metadata_{i} (rational_series.rs)
];

/// vakint's input-notation heads: internal in a value returned in the
/// caller's notation (the direct entry and the vacuum entry's external
/// vectors without a leg return vakint notation: then allowed).
const NOTATION_HEADS: [&str; 12] = [
    "p",
    "k",
    "dot",
    "g",
    "topo",
    "prop",
    "edge",
    "pow",
    "uvprop",
    "EuclideanDot",
    "uedge",
    "vkdot",
];

/// The reviewed list of the OTHER names the crate's code creates in vakint's
/// namespace: vakint's output vocabulary (constants, functions, master and
/// RustRed tags, epsilon names, Retain names) and generic formatting sites;
/// allowed in a value. The tripwire requires every created name to be here or
/// internal. Test code: the check itself allows every name that is not
/// internal.
#[cfg(test)]
const REVIEWED_OUTPUT: [&str; 203] = [
    // base placeholders without evidence of internality (out of the deny table)
    "momentum_match_source_1",
    "vakint_master_transport_1",
    // names a backend can keep in a returned value (out of the deny table):
    // a Retain value's RustRed family tags, unsubstituted masters z<n>,
    // loop-momentum and rotated-momentum names, pinched topology names
    "rustred_family_1",
    "z1",
    "k1",
    "krotated1",
    "1_pinch_1",
    "D6",
    "EulerGamma",
    "GPLImag",
    "GPLReal",
    "HI",
    "HR",
    "NumericalSource",
    "PG",
    "PGamma",
    "PMaster",
    "PMasterValue",
    "PResidual",
    "PolylogImag",
    "PolylogReal",
    "RustRedCoefficient",
    "RustRedConstant",
    "RustRedFiniteLiteral",
    "RustRedMaster",
    "S2",
    "Sqrt",
    "T1ep",
    "T1ep2",
    "Zeta",
    "cl2",
    "cl4",
    "om1",
    "s6",
    "sqrt3",
    "vakintFloat",
    "𝑖",
    "ε",
    "ep",
    "d",
    "MZVvalue",
    "mzv_indices",
    "zeta",
    "pi",
    "Oep",
    "PR10D2",
    "PR10D2Eps1",
    "PR15D2",
    "PR7D2",
    "PR7D2Eps1",
    "Master",
    "GschemeConstants",
    "logMassiveClosure",
    "logQ2mu",
    "log_mu_sq",
    "p_master_basis_eps",
    "delta",
    "rat",
    "t",
    "Interval",
    "K5LogMoment",
    "K5Moment",
    // FMFT masters, their exact parts and Laurent constants (kept when the backend's expansion or substitution is off), and the arguments of their ALARM tails Oep(k, ..)
    "I_PR0",
    "I_PR1",
    "I_PR10",
    "I_PR11",
    "I_PR11d",
    "I_PR12",
    "I_PR13",
    "I_PR14",
    "I_PR15",
    "I_PR2",
    "I_PR3",
    "I_PR4",
    "I_PR4d",
    "I_PR5",
    "I_PR6",
    "I_PR7",
    "I_PR8",
    "I_PR9",
    "I_PR9d",
    "I_PR9x",
    "PR0",
    "PR0ep0",
    "PR0ep1",
    "PR1",
    "PR10",
    "PR10_exact",
    "PR10ep0",
    "PR10ep1",
    "PR11",
    "PR11_exact",
    "PR11ep0",
    "PR11ep1",
    "PR12",
    "PR12ep0",
    "PR12ep1",
    "PR13",
    "PR13_exact",
    "PR13ep0",
    "PR13ep1",
    "PR14",
    "PR14_exact",
    "PR14ep0",
    "PR14ep1",
    "PR15",
    "PR15_exact",
    "PR15ep0",
    "PR15ep1",
    "PR1_exact",
    "PR1ep1",
    "PR1ep2",
    "PR2",
    "PR2_exact",
    "PR2ep1",
    "PR2ep2",
    "PR3",
    "PR3_exact",
    "PR3ep1",
    "PR3ep2",
    "PR4",
    "PR4_exact",
    "PR4d",
    "PR4d_exact",
    "PR4dep1",
    "PR4dfin",
    "PR4ep1",
    "PR4fin",
    "PR5",
    "PR5_exact",
    "PR5ep0",
    "PR5ep1",
    "PR6",
    "PR6_exact",
    "PR6ep0",
    "PR6ep1",
    "PR7",
    "PR7_exact",
    "PR7ep1",
    "PR8",
    "PR8_exact",
    "PR8ep0",
    "PR8ep1",
    "PR9",
    "PR9_exact",
    "PR9d",
    "PR9d_exact",
    "PR9dep0",
    "PR9ep0",
    "PR9ep1",
    "PR9x",
    "T1ep4",
    "iD6",
    // MATAD master heads mi.., their coefficient heads ..bar(k) (kept when expansion or substitution is off), and the arguments of their ALARM tails Oep(k, ..)
    "BN",
    "BN1",
    "D3",
    "D4",
    "D5",
    "DM",
    "DN",
    "E3",
    "miBN",
    "miBN1",
    "miD3",
    "miD4",
    "miD5",
    "miD6",
    "miDM",
    "miDN",
    "miE3",
    "miT111",
    "BN1bar",
    "BNbar",
    "D3bar",
    "D4bar",
    "D5bar",
    "D6bar",
    "DMbar",
    "DNbar",
    "E3bar",
    "T1",
    "T111bar",
    "T111",
    // MATAD's stand-in for the UV mass (substituted back before returning, but a one-letter name a settings string can produce: not internal)
    "M",
    // short names of the known topologies (an integral left unevaluated keeps them)
    "I1L",
    "I2L",
    "I3L",
    "I4L_BMW",
    "I4L_FG",
    "I4L_H",
    "I4L_X",
    // the five-loop topology heads I5L_{sector} (the format placeholder filled by a digit)
    "I5L_1",
    // alphaLoop master constants (kept when masters are not substituted)
    "MASTER2m22EP",
    "MASTER3m001111FIN",
    "MASTER3m011101EP",
    "MASTER3m011111FIN",
    "MASTER3m111111FINPIECE",
    // the normalisation-factor notation's I and eps, and the FORM dummy index family mu{n}: names a settings string can produce, never internal
    "I",
    "eps",
    "mu1",
    // Gamma-function heads of the master expansions (turned into Symbolica's gamma only after the expansion)
    "Gam",
    "iGam",
    "Gamma",
];

/// The caller notation's structural heads.
const CALLER_HEADS: [&str; 6] = [
    names::METRIC,
    names::DOT,
    names::MINK,
    names::BIS,
    names::GAMMA,
    names::MOMENTUM,
];

/// Whether `name` (a vakint-namespace name without its namespace) is one of
/// vakint's INTERNAL symbols: exact, a family, a pattern wildcard (ending in
/// `_`), or an input-notation head where the entry returns the caller's
/// notation.
pub(crate) fn internal(name: &str, vakint_notation: bool) -> bool {
    name.ends_with('_')
        || INTERNAL_EXACT.contains(&name)
        || (!vakint_notation && NOTATION_HEADS.contains(&name))
        || INTERNAL_FAMILIES
            .iter()
            .any(|(prefix, _)| name.starts_with(prefix))
}

fn collect(atom: &Atom, out: &mut std::collections::BTreeSet<Symbol>) {
    atom.as_view().visitor(&mut |part| {
        match part {
            AtomView::Var(v) => {
                out.insert(v.get_symbol());
            }
            AtomView::Fun(f) => {
                out.insert(f.get_symbol());
            }
            _ => {}
        }
        true
    });
}

/// The symbols of every atom of `value` that enters the evaluation (edge
/// factors and masses, vertex factors, ports, the projector, the prefactor
/// fields, insertions, the closure line), of `prefactor` and of the settings.
/// The WEIGHTS are left out on purpose: provenance only, never multiplied into
/// the value, so a symbol written only in a weight must not appear in it.
fn inputs(
    value: &FactorizedNumerator,
    prefactor: &Atom,
    settings: &VakintSettings,
) -> std::collections::BTreeSet<Symbol> {
    let mut s = std::collections::BTreeSet::new();
    for edge in &value.edges {
        collect(&edge.factor, &mut s);
        collect(&edge.mass_squared, &mut s);
    }
    for vertex in &value.vertices {
        collect(&vertex.factor, &mut s);
    }
    for port in &value.ports {
        collect(port, &mut s);
    }
    if let Some(projector) = &value.projector {
        collect(&projector.tensor, &mut s);
    }
    for part in [
        &value.prefactor.overall,
        &value.prefactor.scalars,
        &value.prefactor.colour,
    ] {
        collect(part, &mut s);
    }
    for insertion in &value.insertions {
        collect(&insertion.coefficient, &mut s);
        collect(&insertion.factor, &mut s);
    }
    if let FactorizedKind::Closure(spec) = &value.kind {
        collect(&spec.mass_squared, &mut s);
        collect(&spec.auxiliary_mass_squared, &mut s);
    }
    collect(prefactor, &mut s);
    s.insert(vk_symbol!(&settings.epsilon_symbol));
    s.insert(vk_symbol!(&settings.mu_r_sq_symbol));
    s.remove(&value.dimension);
    s
}

/// Whether `name` (with its namespace) is one of vakint's data atoms
/// ([`super::VAKINT_DATA_ATOMS`]).
fn data_atom(name: &str) -> bool {
    super::VAKINT_DATA_ATOMS
        .iter()
        .any(|(atom, _)| *atom == name)
}

/// The first offending atom of `result`, as the refusal detail.
fn offence(
    view: AtomView<'_>,
    dimension: Symbol,
    known: &std::collections::BTreeSet<Symbol>,
    vakint_notation: bool,
    allow_unknown: bool,
) -> Option<String> {
    let name_ok = |symbol: Symbol| -> Option<String> {
        let name = symbol.get_name();
        if symbol == dimension {
            return Some(format!(
                "the dimension {name} outside a representation slot"
            ));
        }
        if data_atom(name) {
            // an unknown master coefficient: allowed by the settings, or
            // carried by an input
            if symbol == crate::truncation_guard::unknown_coefficient_symbol()
                && !allow_unknown
                && !known.contains(&symbol)
            {
                return Some(format!(
                    "vakint's data atom {name}, but the settings do not allow unknown master coefficients"
                ));
            }
            return None;
        }
        if let Some(local) = name.strip_prefix("vakint::") {
            if internal(local, vakint_notation) {
                return Some(format!("vakint's internal symbol {name}"));
            }
            return None;
        }
        if known.contains(&symbol)
            || name.starts_with("symbolica::")
            || CALLER_HEADS.contains(&name)
        {
            return None;
        }
        Some(format!(
            "the symbol {name}, which is neither an input nor a known output"
        ))
    };
    match view {
        AtomView::Num(_) => None,
        AtomView::Var(v) => name_ok(v.get_symbol()),
        AtomView::Fun(f) => {
            let symbol = f.get_symbol();
            if let Some(found) = name_ok(symbol) {
                return Some(found);
            }
            let rep = super::is_rep_name(symbol.get_name());
            for (i, arg) in f.iter().enumerate() {
                // the dimension argument of a representation slot
                if rep && i == 0 && matches!(arg, AtomView::Var(v) if v.get_symbol() == dimension) {
                    continue;
                }
                if let Some(found) = offence(arg, dimension, known, vakint_notation, allow_unknown)
                {
                    return Some(found);
                }
            }
            None
        }
        AtomView::Add(a) => a
            .iter()
            .find_map(|t| offence(t, dimension, known, vakint_notation, allow_unknown)),
        AtomView::Mul(m) => m
            .iter()
            .find_map(|t| offence(t, dimension, known, vakint_notation, allow_unknown)),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            offence(base, dimension, known, vakint_notation, allow_unknown)
                .or_else(|| offence(exponent, dimension, known, vakint_notation, allow_unknown))
        }
    }
}

/// A number that is no value: a non-finite or indeterminate number
/// (symbolica's complex infinity ⧞, indeterminate ¿, from a division by zero
/// or 0/0), or a power of a literal ZERO base whose exponent is not a positive
/// rational number (0^(-1), 0^(-sqrt(2)), 0^(-π), 0^y: kept unevaluated by the
/// normaliser): the first one found.
fn undefined_number(view: AtomView<'_>) -> Option<String> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        match part {
            AtomView::Num(n)
                if matches!(
                    n.get_coeff_view(),
                    CoefficientView::Infinity(..) | CoefficientView::Indeterminate
                ) =>
            {
                found = Some(format!(
                    "the non-finite number {}",
                    part.to_owned().to_canonical_string()
                ));
            }
            // a floating-point coefficient (Numerical mode) with a NaN or
            // infinite real or imaginary part: Symbolica keeps it under its
            // own Float tag, not as Infinity / Indeterminate
            AtomView::Num(n)
                if matches!(
                    n.get_coeff_view(),
                    CoefficientView::Float(re, im)
                        if !re.to_float().is_finite() || !im.to_float().is_finite()
                ) =>
            {
                found = Some(format!(
                    "the non-finite float {}",
                    part.to_owned().to_canonical_string()
                ));
            }
            // the data atom MZV with indices vakint's evaluator refuses (a
            // leading 1 diverges; non-integer or non-positive indices)
            AtomView::Fun(f)
                if f.get_symbol().get_name() == "vakint::MZV"
                    && !f
                        .iter()
                        .map(|a| i64::try_from(a).ok())
                        .collect::<Option<Vec<i64>>>()
                        .is_some_and(|n| crate::p_integral::mzv::evaluable_indices(&n)) =>
            {
                found = Some(format!(
                    "the divergent or malformed data atom {} (an MZV's first index is at least 2)",
                    part.to_owned().to_canonical_string()
                ));
            }
            AtomView::Pow(p) => {
                let (base, exponent) = p.get_base_exp();
                let positive = symbolica::domains::rational::Rational::try_from(exponent)
                    .is_ok_and(|r| !r.is_negative() && !r.is_zero());
                if base.is_zero() && !positive {
                    found = Some(format!(
                        "the undefined power {} (a zero base)",
                        part.to_owned().to_canonical_string()
                    ));
                }
            }
            _ => {}
        }
        found.is_none()
    });
    found
}

/// THE RESULT CHECK (module documentation) of the value `result` an entry is
/// about to return for the MARKED value `value` (`vakint_notation`: the
/// entry returns vakint notation). A no-op for an unmarked value.
pub(crate) fn check_result(
    value: &FactorizedNumerator,
    prefactor: &Atom,
    settings: &VakintSettings,
    result: &Atom,
    vakint_notation: bool,
) -> Result<()> {
    if !super::is_marked(value) {
        return Ok(());
    }
    if let Some(found) = undefined_number(result.as_view()) {
        return Err(FactorizedError::Evaluation {
            detail: format!("result check: {found} in the value of {}", value.label),
        });
    }
    let known = inputs(value, prefactor, settings);
    match offence(
        result.as_view(),
        value.dimension,
        &known,
        vakint_notation,
        settings.allow_unknown_master_coefficients,
    ) {
        None => Ok(()),
        Some(found) => Err(FactorizedError::Evaluation {
            detail: format!("result check: {found} in the value of {}", value.label),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(text: &str) -> Atom {
        crate::Vakint::initialize_vakint_symbols();
        crate::utils::vakint_macros::vk_parse!(text).unwrap()
    }

    fn bare(record: Option<super::super::DirectRecord>) -> FactorizedNumerator {
        FactorizedNumerator {
            label: "result check test".into(),
            kind: FactorizedKind::Direct,
            dimension: FactorizedNumerator::default_dimension(),
            loop_edges: Vec::new(),
            external_edges: Vec::new(),
            edges: Vec::new(),
            vertices: Vec::new(),
            ports: vec![
                atom("spenso::mink(gammalooprs::{}::dim,gammalooprs::hedge(1,1))"),
                atom("spenso::mink(gammalooprs::{}::dim,gammalooprs::hedge(2,1))"),
            ],
            projector: None,
            prefactor: super::super::Prefactor {
                scalars: atom("UFO::gs^2"),
                ..Default::default()
            },
            colour_branch: None,
            insertions: Vec::new(),
            record,
        }
    }

    fn marked() -> FactorizedNumerator {
        bare(Some(super::super::DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }))
    }

    fn settings() -> VakintSettings {
        VakintSettings {
            epsilon_symbol: "gammalooprs::ε".into(),
            mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
            ..Default::default()
        }
    }

    fn check(result: &str) -> Result<()> {
        let v = marked();
        check_result(&v, &v.prefactor.total(), &settings(), &atom(result), false)
    }

    fn refused(result: &str, prefix: &str) {
        match check(result) {
            Err(FactorizedError::Evaluation { detail }) => assert!(
                detail.starts_with(&format!("result check: {prefix}")),
                "{result}: {detail}"
            ),
            other => panic!("{result}: {other:?}"),
        }
    }

    /// The backstop ALONE (no validation involved): the value of the leaking
    /// Taylor parameter, a placeholder, a dimension outside a slot (inside a
    /// log, as a factor), a symbol from nowhere are refused; an input symbol,
    /// the settings' symbols, Symbolica's constants and functions, the
    /// caller notation's heads, the dimension inside a slot are accepted.
    #[test]
    fn the_result_check_alone() {
        let dim = "gammalooprs::{}::dim";
        let eps = "gammalooprs::{}::ε";
        refused(
            &format!("UFO::gs^2*log(UFO::gs*(1-vakint::closure_lambda^2))*{eps}^-1"),
            "vakint's internal symbol vakint::closure_lambda",
        );
        refused(
            "UFO::gs^2*vakint::closure_tensor_0",
            "vakint's internal symbol vakint::closure_tensor_0",
        );
        refused(
            &format!("UFO::gs^2*log(UFO::gs*({dim}-2))"),
            &format!("the dimension {} outside", "gammalooprs::dim"),
        );
        refused(&format!("UFO::gs^2*{dim}"), "the dimension");
        refused(
            "UFO::gs^2*elsewhere::nowhere",
            "the symbol elsewhere::nowhere, which is neither",
        );
        refused(
            "vakint::UNKNOWN(vakint::m1,0)*UFO::gs^2",
            "vakint's data atom vakint::UNKNOWN, but the settings do not allow",
        );
        for ok in [
            format!("UFO::gs^2*{eps}^-2*𝜋^-4*log(gammalooprs::μᵣ²)"),
            format!(
                "UFO::gs^2*spenso::dot(gammalooprs::Q(0,spenso::mink({dim})),gammalooprs::Q(0,spenso::mink({dim})))*symbolica::zeta(3)"
            ),
            format!(
                "UFO::gs^2*spenso::g(spenso::mink({dim},gammalooprs::hedge(1,1)),spenso::mink({dim},gammalooprs::hedge(2,1)))"
            ),
            format!("UFO::gs^2*vakint::PMaster(vakint::t105)*symbolica::gamma(1+{eps})*𝑒^({eps})"),
            "UFO::gs^2*vakint::HI(1,2)*vakint::Zeta(3)".to_string(),
        ] {
            check(&ok).unwrap_or_else(|e| panic!("{ok}: {e}"));
        }
        // UNKNOWN when the settings allow it
        let v = marked();
        let allow = VakintSettings {
            allow_unknown_master_coefficients: true,
            ..settings()
        };
        check_result(
            &v,
            &v.prefactor.total(),
            &allow,
            &atom("vakint::UNKNOWN(vakint::m1,0)*UFO::gs^2"),
            false,
        )
        .unwrap();
        // UNKNOWN carried by an input (here the prefactor argument) is an
        // input; MZV is output; no data atom is internal
        check_result(
            &v,
            &atom("vakint::UNKNOWN(vakint::ggppss_int_1,0)*UFO::gs^2"),
            &settings(),
            &atom("vakint::UNKNOWN(vakint::ggppss_int_1,0)*UFO::gs^2"),
            false,
        )
        .unwrap();
        check("UFO::gs^2*vakint::MZV(5,3)").unwrap();
        for (name, _) in super::super::VAKINT_DATA_ATOMS {
            let local = name.strip_prefix("vakint::").unwrap();
            assert!(!internal(local, false) && !internal(local, true), "{name}");
        }
        // an unmarked value is not checked
        let unmarked = bare(None);
        check_result(
            &unmarked,
            &Atom::num(1),
            &settings(),
            &atom("vakint::closure_lambda"),
            false,
        )
        .unwrap();
    }

    /// The backstop alone refuses a number that is no value in a marked
    /// result: the zero base under -sqrt(2) and -π (the forms an identically
    /// vanishing base under such an exponent leaves), 0^y, a division by a
    /// literal zero (complex infinity) and 0/0 (indeterminate); the twins with
    /// a base 1 (1^(-sqrt(2)) = 1) pass. Floating-point coefficients (the
    /// Numerical mode's): NaN, +inf and -inf (as real or imaginary part) are
    /// refused ("the non-finite float"), 1.5 and 1e300 pass.
    #[test]
    fn a_number_that_is_no_value_is_refused() {
        let eps = "gammalooprs::ε";
        for result in [
            format!("{eps}^-2*0^(-2^(1/2))"),
            format!("{eps}^-2*0^(-𝜋)"),
            format!("{eps}^-2*0^witness::y"),
            format!("{eps}^-2*(witness::x-witness::x)^-1"),
            format!("{eps}^-2*0/0"),
        ] {
            let error = check(&result).unwrap_err();
            match &error {
                FactorizedError::Evaluation { detail } => {
                    assert!(detail.starts_with("result check: "), "{result}: {detail}")
                }
                other => panic!("{result}: expected Evaluation, got {other:?}"),
            }
        }
        for result in [format!("{eps}^-2*1^(-2^(1/2))"), format!("{eps}^-2*2^(-𝜋)")] {
            check(&result).unwrap_or_else(|e| panic!("{result}: {e}"));
        }
        let v = marked();
        let pole = atom(&format!("{eps}^-2"));
        let with =
            |c: Atom| check_result(&v, &v.prefactor.total(), &settings(), &(&pole * c), false);
        for (what, c) in [
            ("NaN", Atom::num(f64::NAN)),
            ("+inf", Atom::num(f64::INFINITY)),
            ("-inf", Atom::num(f64::NEG_INFINITY)),
            ("i NaN", Atom::num(f64::NAN) * Atom::i()),
        ] {
            match with(c) {
                Err(FactorizedError::Evaluation { detail }) => assert!(
                    detail.starts_with("result check: the non-finite float "),
                    "{what}: {detail}"
                ),
                other => panic!("{what}: expected Evaluation, got {other:?}"),
            }
        }
        for c in [1.5, 1e300] {
            with(Atom::num(c)).unwrap_or_else(|e| panic!("{c}: {e}"));
        }
        // the data atom MZV: a leading index 1 diverges (vakint's evaluator
        // refuses it); the convergent ones pass
        for bad in ["vakint::MZV(1)", "vakint::MZV(1,2)", "vakint::MZV(0,3)"] {
            refused(
                &format!("{eps}^-2*{bad}"),
                "the divergent or malformed data atom ",
            );
        }
        for good in ["vakint::MZV(2)", "vakint::MZV(2,1)", "vakint::MZV(5,3)"] {
            check(&format!("{eps}^-2*{good}")).unwrap_or_else(|e| panic!("{good}: {e}"));
        }
    }

    /// The deny table: the series' own placeholders, the base's long
    /// placeholder names and families, the wildcards and the input-notation
    /// heads are internal; the names a backend can legitimately keep in a
    /// returned value (a Retain value's rustred_family_<hash>, an
    /// unsubstituted AlphaLoop master z<n>, short and topology names) and
    /// vakint's output vocabulary (the massive-vacuum constants, the RustRed
    /// tags) pass a marked result.
    #[test]
    fn internal_names_are_refused_and_the_output_vocabulary_passes() {
        for name in [
            "closure_tensor_12",
            "graph_entry_branch_3",
            "input_coefficient_placeholder_0",
            "pow_",
            "idx___",
            "topo",
            "partial_numerical_input",
        ] {
            assert!(internal(name, false), "{name} is internal");
        }
        assert!(
            !internal("topo", true) && !internal("dot", true),
            "vakint notation"
        );
        for name in [
            "zeta",
            "k",
            "cl2",
            "T1ep",
            "om1",
            "RustRedMaster",
            "vakintFloat",
            "t105",
            "rustred_family_1a2b",
            "z3",
            "k3",
            "krotated1",
            "I5L_32745_pinch_4_6",
            "msq",
            "x",
            "vec",
        ] {
            assert!(!internal(name, false) || name == "k", "{name} is output");
        }
        for ok in [
            "UFO::gs^2*vakint::cl2*vakint::sqrt3",
            "UFO::gs^2*(vakint::T1ep+vakint::T1ep2+vakint::D6+vakint::S2+vakint::s6+vakint::cl4)",
            "UFO::gs^2*vakint::HR(0,0,vakint::om1)",
            "UFO::gs^2*vakint::RustRedMaster(1)*vakint::RustRedCoefficient(2)*vakint::RustRedConstant(3)",
            "UFO::gs^2*vakint::RustRedFiniteLiteral(1)*vakint::vakintFloat(2)*vakint::NumericalSource(3)",
            "UFO::gs^2*vakint::𝑖",
        ] {
            check(ok).unwrap_or_else(|e| panic!("{ok}: {e}"));
        }
        refused(
            "UFO::gs^2*vakint::graph_entry_branch_0",
            "vakint's internal symbol vakint::graph_entry_branch_0",
        );
        refused(
            "UFO::gs^2*vakint::topo(1)",
            "vakint's internal symbol vakint::topo",
        );
    }

    /// The end (exclusive) of the lexical token starting at `i` when it is a
    /// string (plain, byte, raw, raw byte), a character literal or a comment.
    fn token_end(chars: &[char], i: usize) -> Option<usize> {
        let at = |k: usize| chars.get(k).copied();
        let ident = |k: usize| chars[k].is_alphanumeric() || chars[k] == '_';
        match at(i)? {
            '/' if at(i + 1) == Some('/') => Some(
                (i..chars.len())
                    .find(|&k| chars[k] == '\n')
                    .unwrap_or(chars.len()),
            ),
            '/' if at(i + 1) == Some('*') => Some(
                (i + 2..chars.len().saturating_sub(1))
                    .find(|&k| chars[k] == '*' && chars[k + 1] == '/')
                    .map_or(chars.len(), |k| k + 2),
            ),
            '"' => {
                let mut k = i + 1;
                while k < chars.len() && chars[k] != '"' {
                    k += if chars[k] == '\\' { 2 } else { 1 };
                }
                Some(k + 1)
            }
            // r".." / r#".."#, and br".." (the `b` is an identifier character)
            'r' if i == 0
                || !ident(i - 1)
                || (chars[i - 1] == 'b' && (i == 1 || !ident(i - 2))) =>
            {
                let hashes = (i + 1..chars.len())
                    .take_while(|&k| chars[k] == '#')
                    .count();
                if at(i + 1 + hashes) != Some('"') {
                    return None;
                }
                let close: Vec<char> = std::iter::once('"')
                    .chain(std::iter::repeat_n('#', hashes))
                    .collect();
                let from = i + 2 + hashes;
                Some(
                    (from..chars.len())
                        .find(|&k| chars[k..].starts_with(&close))
                        .map_or(chars.len(), |k| k + close.len()),
                )
            }
            // '\n', '\'', '\\', '\u{..}': the closing quote after the escape
            '\'' if at(i + 1) == Some('\\') => Some(
                (i + 3..chars.len())
                    .find(|&k| chars[k] == '\'')
                    .map_or(chars.len(), |k| k + 1),
            ),
            '\'' if at(i + 2) == Some('\'') => Some(i + 3),
            _ => None,
        }
    }

    /// `text` with every item under a line-leading `#[cfg(test)]` or
    /// `#[cfg(all(test, ..))]` blanked: the
    /// attribute up to the first `;` or `,` outside brackets (a `use`, a
    /// `mod x;`, a field) or through the brace block that comes first (a
    /// module, a function, an impl, a statement). Strings, raw strings,
    /// character literals and comments are skipped when counting, so a brace
    /// in them cuts nothing more; any doubt leaves text in (scanned, never
    /// hidden).
    fn without_test_items(text: &str) -> String {
        let chars: Vec<char> = text.chars().collect();
        let attribute =
            regex::Regex::new(r"(?m)^[ \t]*#\[cfg\((?:test|all\(\s*test\s*,[^\]\n]*\))\)\]")
                .unwrap();
        // char offsets of the attributes' ends
        let starts: Vec<(usize, usize)> = attribute
            .find_iter(text)
            .map(|m| {
                (
                    text[..m.start()].chars().count(),
                    text[..m.end()].chars().count(),
                )
            })
            .collect();
        let skip = |i: usize| token_end(&chars, i);
        let mut blank = vec![false; chars.len()];
        for (start, end) in starts {
            let (mut i, mut depth) = (end, 0i64);
            let stop = loop {
                if i >= chars.len() {
                    break None;
                }
                if let Some(next) = skip(i) {
                    i = next;
                    continue;
                }
                match chars[i] {
                    '(' | '[' => depth += 1,
                    ')' | ']' => depth -= 1,
                    ';' | ',' if depth == 0 => break Some(i + 1),
                    '{' if depth == 0 => {
                        let mut braces = 0i64;
                        while i < chars.len() {
                            if let Some(next) = skip(i) {
                                i = next;
                                continue;
                            }
                            match chars[i] {
                                '{' => braces += 1,
                                '}' => braces -= 1,
                                _ => {}
                            }
                            i += 1;
                            if braces == 0 {
                                break;
                            }
                        }
                        break (braces == 0).then_some(i);
                    }
                    _ => {}
                }
                if depth < 0 {
                    break None;
                }
                i += 1;
            };
            // an unbalanced item blanks only its attribute
            for b in &mut blank[start..stop.unwrap_or(end).min(chars.len())] {
                *b = true;
            }
        }
        chars
            .iter()
            .zip(blank)
            .map(|(&c, b)| if b && c != '\n' { ' ' } else { c })
            .collect()
    }

    /// The skip rule hides test items only: a test module, a test function
    /// after a `mod tests;`, a test field, a `use`; code after them and a
    /// brace inside a string, a character or a comment of a test item do not
    /// widen the cut.
    #[test]
    fn only_test_items_are_skipped() {
        let source = r##"use a::{b, c};
#[cfg(test)]
use crate::{x, y};
fn real() { let _ = symbol!("witness::kept_a"); }
#[cfg(test)]
mod tests;
#[cfg(test)]
#[test]
fn single() {
    let _ = '{';
    let _ = "}}}";
    let _ = r#"{"#; // }
    /* { */
    let _ = symbol!("witness::gone_a");
}
struct S {
    #[cfg(test)]
    inputs: Vec<(u8, u8)>,
    kept: u8,
}
#[cfg(test)]
mod inline {
    fn f() { symbol!("witness::gone_b"); }
}
fn after() { let _ = symbol!("witness::kept_b"); }
"##;
        let kept = without_test_items(source);
        for name in [
            "kept_a",
            "kept_b",
            "use a::{b, c};",
            "kept: u8,",
            "fn after()",
        ] {
            assert!(kept.contains(name), "{name} is kept:\n{kept}");
        }
        for name in [
            "gone_a",
            "gone_b",
            "crate::{x, y}",
            "mod tests;",
            "inputs",
            "fn single",
            "mod inline",
        ] {
            assert!(!kept.contains(name), "{name} is skipped:\n{kept}");
        }
        assert_eq!(kept.lines().count(), source.lines().count());
    }

    /// The symbol constructors the tripwire reads: macros (by the last segment
    /// of their path, followed by `!`) ...
    const MACRO_CONSTRUCTORS: [&str; 9] = [
        "vk_symbol",
        "vk_symbol_uncached",
        "vakint_symbol",
        "symbol",
        "vk_parse",
        "vakint_parse",
        "parse",
        "try_parse",
        "parse_lit",
    ];
    /// .. and associated functions (by their last two path segments).
    const FUNCTION_CONSTRUCTORS: [&str; 5] = [
        "NamespacedSymbol::parse",
        "NamespacedSymbol::try_parse",
        "NamespacedSymbol::try_parse_lit",
        "SymbolBuilder::new",
        "Atom::parse",
    ];
    /// The constructors whose argument is ONE symbol name (`symbol!`: every
    /// positional argument before a `;`); the others take a text that is
    /// parsed as an expression.
    const NAME_CONSTRUCTORS: [&str; 8] = [
        "vk_symbol!",
        "vk_symbol_uncached!",
        "vakint_symbol!",
        "symbol!",
        "NamespacedSymbol::parse",
        "NamespacedSymbol::try_parse",
        "NamespacedSymbol::try_parse_lit",
        "SymbolBuilder::new",
    ];

    /// The reviewed constructor sites whose argument the tripwire cannot
    /// resolve to a literal: (file under src/, constructor, argument text with
    /// whitespace collapsed, number of such sites in the file, reason). The
    /// reasons: `settings` (a symbol or a normalisation named by
    /// [`VakintSettings`]: the epsilon or the scale symbol, alone or formatted
    /// into a literal text), `data` (a name or a text from vakint's own data:
    /// tables, catalogues, registries, records, files, external tools' output,
    /// or a literal list / literal argument in vakint's code that the
    /// tripwire does not follow), `caller` (a name or a text the API caller
    /// supplies), `format` (a format string whose placeholder is a whole name
    /// or expression: a constant's name, an index, an atom's text), `macro` (a
    /// parameter `$x` of one of the
    /// listed constructor macros' own definitions: their call sites are read
    /// as that constructor). A site missing here, a count that differs and an
    /// entry no site matches all fail.
    #[rustfmt::skip]
    const REVIEWED_SITES: [(&str, &str, &str, usize, &str); 171] = [
        ("display.rs", "NamespacedSymbol::parse", "&self.name", 1, "data"),
        ("display.rs", "NamespacedSymbol::parse", "full_name", 1, "data"),
        ("display.rs", "SymbolBuilder::new", "NamespacedSymbol::parse(&self.name)", 1, "data"),
        ("fmft.rs", "vk_parse!", "&settings.epsilon_symbol", 1, "settings"),
        ("fmft.rs", "vk_parse!", "format!( \"({})^{}\", muv_sq_atom.to_canonical_string(), muv_sq_dimension ) .as_str()", 1, "format"),
        ("fmft.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("fmft.rs", "vk_symbol!", "*name", 1, "data"),
        ("fmft.rs", "vk_symbol!", "self.settings.epsilon_symbol.as_str()", 1, "settings"),
        ("fmft/finalize.rs", "vk_parse!", "format!( \"((𝑖*(𝜋^((4-2*{eps})/2)))\\ *(exp(-EulerGamma))^({eps})\\ *(exp(-logmUVmu-log_mu_sq))^({eps}))^{loop_count}\", eps = settings.epsilon_symbol, ) .as_str()", 1, "settings"),
        ("fmft/finalize.rs", "vk_symbol!", "name", 2, "data"),
        ("fmft/finalize.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 2, "settings"),
        ("fmft/finalize.rs", "vk_symbol!", "settings.mu_r_sq_symbol.as_str()", 2, "settings"),
        ("inspection.rs", "vk_parse!", "value", 1, "caller"),
        ("inspection/constants.rs", "vk_symbol!", "name", 3, "data"),
        ("inspection/definitions.rs", "vk_parse!", "&text", 1, "data"),
        ("inspection/definitions.rs", "vk_symbol!", "bar", 1, "data"),
        ("inspection/definitions.rs", "vk_symbol!", "name", 1, "data"),
        ("inspection/registry.rs", "vk_parse!", "expression", 1, "data"),
        ("inspection/registry.rs", "vk_parse!", "formula", 3, "data"),
        ("inspection/registry.rs", "vk_symbol!", "&digest", 1, "format"),
        ("lib.rs", "vk_parse!", "&expr.to_expression()", 1, "settings"),
        ("lib.rs", "vk_parse!", "&settings.epsilon_symbol", 7, "settings"),
        ("lib.rs", "vk_parse!", "&vakint_settings.epsilon_symbol", 1, "settings"),
        ("lib.rs", "vk_parse!", "format!( \"( 𝑖*(𝜋^((4-2*{eps})/2))\\ )^{n_loops}\", eps = settings.epsilon_symbol, n_loops = n_loops ) .as_str()", 1, "settings"),
        ("lib.rs", "vk_parse!", "format!( \"(\\ 𝑖*(𝜋^((4-2*{eps})/2))\\ * (exp(-EulerGamma))^({eps})\\ * (exp(-logmUVmu-log_mu_sq))^({eps})\\ )^{n_loops}\", eps = settings.epsilon_symbol, n_loops = integral.n_loops ) .as_str()", 1, "settings"),
        ("lib.rs", "vk_parse!", "format!( \"prop(id{id}_,uedge(n{id_l_node}_,n{id_r_node}_),q{id}__,{mass},{power})\", id = i_prop, id_l_node = id_node_left, id_r_node = id_node_right, mass = mass_symbol_string, power = pow_symbol_string ) .as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!( \"prop(id{id}_,uedge(n{id_l_node}l_,n{id_r_node}r_),q{id}__,{mass},{power})\", id = i_prop, id_l_node = id_node_left, id_r_node = id_node_right, mass = mass_symbol_string, power = pow_symbol_string ) .as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!( \"vec{}({}{},idx_)\", if *vec == EXTERNAL_MOMENTUM_SYMBOL { \"1\" } else { \"\" }, vec, id ) .as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!( \"{}{}(mu{})*{}{}(mu{})\", a_in[0], a_in[1], a_in[4], a_in[2], a_in[3], a_in[4], ) .as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"k({})\", i_loop).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"k({},idx___)\", i_loop).as_str()", 2, "format"),
        ("lib.rs", "vk_parse!", "format!(\"msq({})\", i_prop).as_str()", 2, "format"),
        ("lib.rs", "vk_parse!", "format!(\"msq({})\", var_prop_id).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"pow({})\", i_prop).as_str()", 2, "format"),
        ("lib.rs", "vk_parse!", "format!(\"pow({})\", prop_id).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"pow({})\", var_prop_id).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"prop({},args__)\", i_prop).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}(1,1)\", DOT_SYMBOL)", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}(idx1_,idx2_)\", METRIC_SYMBOL)", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}(idx1_,idx2_)\", METRIC_SYMBOL).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}(idx1_,idx2_)\", metric_symbol)", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}(momID_,idx___)\", EXTERNAL_MOMENTUM_SYMBOL).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}({})\", v, v_id).as_str()", 2, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}({})\", vec, id).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}({},idx_)\", vec, id).as_str()", 3, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}{}\", v, v_id).as_str()", 2, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}{}\", vec, id).as_str()", 1, "format"),
        ("lib.rs", "vk_parse!", "format!(\"{}{}(idx_)\", vec, id).as_str()", 2, "format"),
        ("lib.rs", "vk_parse!", "out.replace(\"I\", \"𝑖\").as_str()", 1, "data"),
        ("lib.rs", "vk_parse!", "settings.mu_r_sq_symbol.as_str()", 1, "settings"),
        ("lib.rs", "vk_parse!", "symb.as_str()", 2, "caller"),
        ("lib.rs", "vk_symbol!", "&param_first_form.clone()", 1, "caller"),
        ("lib.rs", "vk_symbol!", "&self.epsilon_symbol", 1, "settings"),
        ("lib.rs", "vk_symbol!", "&settings.epsilon_symbol", 9, "settings"),
        ("lib.rs", "vk_symbol!", "format!(\"{}{}\", v, v_id)", 1, "format"),
        ("lib.rs", "vk_symbol!", "metric_symbol", 1, "data"),
        ("lib.rs", "vk_symbol!", "property", 1, "data"),
        ("lib.rs", "vk_symbol!", "self.name.clone()", 1, "data"),
        ("lib.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 7, "settings"),
        ("lib.rs", "vk_symbol!", "settings.epsilon_symbol.clone()", 1, "settings"),
        ("lib.rs", "vk_symbol!", "settings.mu_r_sq_symbol.as_str()", 4, "settings"),
        ("lib.rs", "vk_symbol!", "settings.mu_r_sq_symbol.clone()", 1, "settings"),
        ("lib.rs", "vk_symbol!", "vec", 1, "data"),
        ("matad.rs", "vk_parse!", "&settings.epsilon_symbol", 1, "settings"),
        ("matad.rs", "vk_parse!", "format!( \"( (𝑖*(𝜋^((4-2*{eps})/2)))\\ * (exp(-EulerGamma))^({eps})\\ * (exp(-logmUVmu-log_mu_sq))^({eps})\\ )^{loop_count}\", eps = settings.epsilon_symbol, ) .as_str()", 1, "settings"),
        ("matad.rs", "vk_parse!", "format!(\"((-1)^{sign_power})\").as_str()", 1, "format"),
        ("matad.rs", "vk_parse!", "format!(\"({})^(pow_/2)\", muv_sq_atom.to_canonical_string()).as_str()", 1, "format"),
        ("matad.rs", "vk_parse!", "format!(\"({})^pow_\", muv_atom.to_canonical_string()).as_str()", 1, "format"),
        ("matad.rs", "vk_parse!", "format!(\"{}*dot(p{},p{})\", -orientation1 * orientation2, i_edge1, i_edge2).as_str()", 1, "format"),
        ("matad.rs", "vk_parse!", "muv_atom.to_canonical_string().as_str()", 1, "caller"),
        ("matad.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("matad.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 2, "settings"),
        ("matad.rs", "vk_symbol!", "settings.mu_r_sq_symbol.as_str()", 2, "settings"),
        ("p_integral.rs", "vk_parse!", "data.lines", 1, "data"),
        ("p_integral.rs", "vk_parse!", "family.lines", 1, "data"),
        ("p_integral.rs", "vk_symbol!", "name", 1, "data"),
        ("p_integral.rs", "vk_symbol!", "self.data.name", 1, "data"),
        ("p_integral/boundary.rs", "vk_parse!", "text", 1, "caller"),
        ("p_integral/factorized.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("p_integral/factorized.rs", "vk_symbol!", "&settings.mu_r_sq_symbol", 1, "settings"),
        ("p_integral/factorized.rs", "vk_symbol!", "epsilon", 1, "settings"),
        ("p_integral/factorized/backstop.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("p_integral/factorized/backstop.rs", "vk_symbol!", "&settings.mu_r_sq_symbol", 1, "settings"),
        ("p_integral/factorized/closure.rs", "vk_symbol!", "&settings.epsilon_symbol", 3, "settings"),
        ("p_integral/factorized/evaluate.rs", "vk_symbol!", "&settings.epsilon_symbol", 5, "settings"),
        ("p_integral/factorized/graph.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("p_integral/factorized/graph.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 1, "settings"),
        ("p_integral/factorized/text.rs", "vk_parse!", "text", 1, "caller"),
        ("p_integral/factorized/vacuum.rs", "vk_symbol!", "&settings.epsilon_symbol", 4, "settings"),
        ("p_integral/finalize.rs", "vk_parse!", "&format!( \"(𝑖*𝜋^(2-({eps}))*Gamma(1+{eps})*Gamma(1-{eps})^2/Gamma(2-2*{eps}))^{l}\\ *exp(-logQ2mu-log_mu_sq)^({eps}*{l})\", eps = settings.epsilon_symbol )", 1, "settings"),
        ("p_integral/finalize.rs", "vk_symbol!", "&settings.epsilon_symbol", 2, "settings"),
        ("p_integral/finalize.rs", "vk_symbol!", "&settings.mu_r_sq_symbol", 1, "settings"),
        ("p_integral/forcer.rs", "vk_parse!", "&cleaned", 1, "data"),
        ("p_integral/forcer.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("p_integral/forcer/batch.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("p_integral/forcer/fixtures.rs", "vk_parse!", "text", 1, "data"),
        ("p_integral/forcer/fixtures.rs", "vk_symbol!", "&settings.epsilon_symbol", 2, "settings"),
        ("p_integral/forcer/fixtures.rs", "vk_symbol!", "&settings.mu_r_sq_symbol", 1, "settings"),
        ("p_integral/forcer/fixtures.rs", "vk_symbol!", "epsilon", 1, "settings"),
        ("p_integral/input.rs", "vk_symbol!", "family.name", 1, "data"),
        ("p_integral/massive.rs", "vk_parse!", "&format!( \"(𝑖*𝜋^(2-({eps})))^{loops}\\ *(Gamma(1+{eps})*Gamma(1-{eps})^2/Gamma(2-2*{eps}))^{cut_loops}\\ *exp(-logMassiveClosure)^({eps}*{loops})\", eps = settings.epsilon_symbol, )", 1, "settings"),
        ("p_integral/massive.rs", "vk_symbol!", "&settings.epsilon_symbol", 7, "settings"),
        ("p_integral/massive.rs", "vk_symbol!", "&settings.mu_r_sq_symbol", 1, "settings"),
        ("p_integral/master_basis.rs", "vk_parse!", "&text", 1, "data"),
        ("p_integral/master_basis.rs", "vk_symbol!", "&symbol(&exchange.old)?", 1, "data"),
        ("p_integral/master_basis.rs", "vk_symbol!", "&symbol(master)?", 1, "data"),
        ("p_integral/master_catalogue.rs", "vk_parse!", "&term.coefficient", 1, "data"),
        ("p_integral/master_catalogue.rs", "vk_parse!", "text", 1, "data"),
        ("p_integral/master_catalogue.rs", "vk_symbol!", "name", 1, "data"),
        ("p_integral/master_catalogue.rs", "vk_symbol!", "symbol", 1, "data"),
        ("p_integral/masters.rs", "vk_parse!", "&body", 1, "data"),
        ("p_integral/masters.rs", "vk_symbol!", "name", 1, "data"),
        ("p_integral/parts.rs", "vk_parse!", "family.lines", 1, "data"),
        ("p_integral/parts.rs", "vk_symbol!", "name", 1, "data"),
        ("precision.rs", "vk_parse!", "name", 2, "caller"),
        ("precision.rs", "vk_symbol!", "&settings.epsilon_symbol", 2, "settings"),
        ("precision.rs", "vk_symbol!", "*name", 1, "data"),
        ("precision/legacy.rs", "vk_parse!", "&original_text", 1, "data"),
        ("precision/materialize.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("refpath.rs", "vk_symbol!", "&settings.epsilon_symbol", 2, "settings"),
        ("refpath/direct.rs", "vk_symbol!", "&settings.epsilon_symbol", 2, "settings"),
        ("rustred_evaluation/experimental/native.rs", "vakint_parse!", "label.as_str()", 2, "data"),
        ("rustred_evaluation/experimental/native.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 1, "settings"),
        ("rustred_evaluation/five_loop.rs", "vk_parse!", "format!( \"(𝑖*𝜋^(2-({eps}))*exp(-EulerGamma*({eps})))^{loop_count}\", eps = settings.epsilon_symbol )", 1, "settings"),
        ("rustred_evaluation/five_loop.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 1, "settings"),
        ("rustred_evaluation/five_loop.rs", "vk_symbol!", "settings.mu_r_sq_symbol.as_str()", 1, "settings"),
        ("rustred_evaluation/five_loop/r2_fixtures.rs", "vakint_parse!", "text", 2, "data"),
        ("rustred_evaluation/five_loop/witness_screen.rs", "vk_parse!", "&settings.epsilon_symbol", 1, "settings"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "&constant.symbol", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "exact", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "expression", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "formula", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "name", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "recorded", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_parse!", "value", 1, "data"),
        ("rustred_evaluation/masters.rs", "vk_symbol!", "&settings.epsilon_symbol", 3, "settings"),
        ("rustred_evaluation/masters.rs", "vk_symbol!", "&tag", 1, "format"),
        ("rustred_evaluation/materialize.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 1, "settings"),
        ("rustred_evaluation/mod.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("rustred_evaluation/numerator.rs", "vk_symbol!", "settings.epsilon_symbol.as_str()", 1, "settings"),
        ("rustred_evaluation/p_data.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("rustred_evaluation/p_data.rs", "vk_symbol!", "residual.node.as_str()", 1, "data"),
        ("rustred_evaluation/p_gate.rs", "try_parse!", "value", 1, "data"),
        ("rustred_evaluation/p_handover.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("rustred_evaluation/p_integral.rs", "vk_parse!", "&text", 1, "data"),
        ("rustred_evaluation/p_integral.rs", "vk_parse!", "coefficient", 1, "data"),
        ("rustred_evaluation/p_integral.rs", "vk_symbol!", "&settings.epsilon_symbol", 5, "settings"),
        ("rustred_evaluation/p_integral.rs", "vk_symbol!", "name.as_str()", 1, "data"),
        ("rustred_evaluation/p_integral.rs", "vk_symbol!", "node.trim()", 1, "data"),
        ("rustred_evaluation/p_lowering.rs", "vk_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("symbolica_community_module/inspection.rs", "vakint_parse!", "&name", 1, "caller"),
        ("symbolica_community_module/mod.rs", "vakint_symbol!", "self.settings.epsilon_symbol.clone()", 2, "settings"),
        ("symbolica_community_module/numerical.rs", "vakint_parse!", "&name", 1, "caller"),
        ("symbolica_community_module/numerical.rs", "vakint_symbol!", "&settings.epsilon_symbol", 1, "settings"),
        ("symbols.rs", "symbol!", "format!(\"{}::{}\",crate::NAMESPACE,METRIC_SYMBOL)", 1, "format"),
        ("symbols.rs", "vk_parse!", "left", 1, "data"),
        ("symbols.rs", "vk_parse!", "right", 1, "data"),
        ("symbols.rs", "vk_parse!", "value", 1, "data"),
        ("symbols.rs", "vk_symbol!", "basis", 1, "data"),
        ("symbols.rs", "vk_symbol!", "name", 3, "data"),
        ("topologies.rs", "vk_parse!", "format!( \"prop({},edge({},{}),args__)\", edge.id, edge.left_node_id, edge.right_node_id )", 1, "format"),
        ("topologies.rs", "vk_parse!", "format!(\"pow({})\", prop_id).as_str()", 1, "format"),
        ("topologies.rs", "vk_parse!", "format!(\"prop({},args__)\", prop_id).as_str()", 3, "format"),
        ("topologies.rs", "vk_parse!", "format!(\"prop({},edge(nl_,nr_),args__)\", edge.id)", 1, "format"),
        ("topologies.rs", "vk_parse!", "format!(\"prop({},edges_,q_,mUVsq_,pow_)\", prop_id).as_str()", 1, "format"),
        ("truncation_guard.rs", "vk_parse!", "master", 1, "data"),
        ("truncation_guard.rs", "vk_symbol!", "master", 1, "data"),
        ("utils.rs", "try_parse!", "$s", 7, "macro"),
        ("utils.rs", "try_parse!", "format!(\"{}\", $s)", 2, "macro"),
        ("utils.rs", "vk_parse!", "s", 2, "caller"),
        ("utils.rs", "vk_symbol_uncached!", "$s", 2, "macro"),
    ];

    /// Names in vakint's namespace that literal arguments create and that are
    /// NOT YET CLASSIFIED: whether each can legitimately appear in a marked
    /// value (then [`REVIEWED_OUTPUT`]) or is a placeholder vakint substitutes
    /// away (then [`INTERNAL_EXACT`] / [`INTERNAL_FAMILIES`]) was not decided
    /// by the tripwire's author (grouped by the files that create them). The
    /// tripwire FAILS while this list is not empty (after all its other
    /// checks); an entry no site creates fails as stale.
    #[rustfmt::skip]
    const UNDECIDED: &[&str] = &[
    ];

    /// Source files under src/ that no module declaration names (not compiled,
    /// hence not read); a new orphan and a stale entry fail.
    const REVIEWED_ORPHANS: [&str; 2] = [
        "rustred_evaluation/masters_grouping_tests.rs",
        "rustred_evaluation/masters_series_tests.rs",
    ];

    /// One non-test source file: its path under src/, its characters, and per
    /// character whether it is blanked (a comment or a test item).
    type Source = (String, Vec<char>, Vec<bool>);

    /// The tripwire's view of the crate's non-test sources: (path under src/,
    /// text). Skipped: files declared by a `#[cfg(test)]` or
    /// `#[cfg(all(test, ..))]` `mod x;` (x.rs, x/mod.rs and everything under
    /// x/, or the file a `#[path = ".."]` names and the directory of its
    /// stem), every item under such an attribute ([`without_test_items`]);
    /// comments are blanked. Returned apart: the ORPHANS, files no `mod x;`
    /// declaration (with or without `#[path]`) of any file names, which the
    /// crate does not compile (they must be on [`REVIEWED_ORPHANS`]).
    fn non_test_sources() -> (Vec<Source>, Vec<String>) {
        let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut files = Vec::new();
        let mut dirs = vec![root.clone()];
        while let Some(dir) = dirs.pop() {
            for entry in std::fs::read_dir(&dir).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    dirs.push(path);
                } else if path.extension().is_some_and(|e| e == "rs") {
                    files.push(path);
                }
            }
        }
        files.sort();
        // attributes, then `mod x;`: test-only when one attribute is
        // `cfg(test)` or `cfg(all(test, ..))`; `#[path = ".."]` names the file
        let declaration = regex::Regex::new(
            r"(?m)^((?:[ \t]*#\[[^\]]*\]\s*)*)[ \t]*(?:pub(?:\([^)]*\))?\s+)?mod\s+(\w+)\s*;",
        )
        .unwrap();
        let test_only =
            regex::Regex::new(r"#\[cfg\((?:test|all\(\s*test\s*,[^\]]*\))\)\]").unwrap();
        let path_attribute = regex::Regex::new(r#"#\[path\s*=\s*"([^"]*)"\s*\]"#).unwrap();
        let (mut test_paths, mut live_files) = (Vec::new(), vec![root.join("lib.rs")]);
        for path in &files {
            let text = std::fs::read_to_string(path).unwrap();
            let file = path.file_name().unwrap().to_string_lossy().to_string();
            let dir = if ["mod.rs", "lib.rs", "main.rs"].contains(&file.as_str()) {
                path.parent().unwrap().to_path_buf()
            } else {
                path.with_extension("")
            };
            for capture in declaration.captures_iter(&text) {
                let declared = if let Some(named) = path_attribute.captures(&capture[1]) {
                    vec![path.parent().unwrap().join(&named[1])]
                } else {
                    vec![
                        dir.join(format!("{}.rs", &capture[2])),
                        dir.join(&capture[2]).join("mod.rs"),
                    ]
                };
                if test_only.is_match(&capture[1]) {
                    for file in declared {
                        test_paths.push(file.with_extension(""));
                        if file.ends_with("mod.rs") {
                            test_paths.push(file.parent().unwrap().to_path_buf());
                        }
                        test_paths.push(file);
                    }
                } else {
                    live_files.extend(declared);
                }
            }
        }
        let relative = |path: &std::path::Path| {
            path.strip_prefix(&root)
                .unwrap()
                .to_string_lossy()
                .replace('\\', "/")
        };
        let mut orphans = Vec::new();
        let sources = files
            .into_iter()
            .filter(|path| !test_paths.iter().any(|test| path.starts_with(test)))
            .filter(|path| {
                let live = live_files.contains(path);
                if !live {
                    orphans.push(relative(path));
                }
                live
            })
            .map(|path| {
                let text = without_test_items(&std::fs::read_to_string(&path).unwrap());
                let (chars, literal) = code_view(&text);
                (relative(&path), chars, literal)
            })
            .collect();
        (sources, orphans)
    }

    /// The characters of `text` with comments blanked (newlines kept), and
    /// whether each character is inside a string or a character literal.
    fn code_view(text: &str) -> (Vec<char>, Vec<bool>) {
        let mut chars: Vec<char> = text.chars().collect();
        let mut literal = vec![false; chars.len()];
        let mut i = 0;
        while i < chars.len() {
            if let Some(end) = token_end(&chars, i) {
                let end = end.min(chars.len());
                if chars[i] == '/' {
                    for c in &mut chars[i..end] {
                        if *c != '\n' {
                            *c = ' ';
                        }
                    }
                } else {
                    for l in &mut literal[i..end] {
                        *l = true;
                    }
                }
                i = end;
            } else {
                i += 1;
            }
        }
        (chars, literal)
    }

    /// The top-level arguments of the bracket opened at `open`: (text,
    /// whether a `;` came before it), and the index after the closing bracket.
    fn arguments(chars: &[char], open: usize) -> (Vec<(String, bool)>, usize) {
        let (mut depth, mut i, mut from) = (0i64, open, open + 1);
        let mut out = Vec::new();
        let mut after_semicolon = false;
        while i < chars.len() {
            if let Some(end) = token_end(chars, i) {
                i = end;
                continue;
            }
            match chars[i] {
                '(' | '[' | '{' => depth += 1,
                ')' | ']' | '}' => {
                    depth -= 1;
                    if depth == 0 {
                        let text: String = chars[from..i].iter().collect();
                        if !text.trim().is_empty() {
                            out.push((text, after_semicolon));
                        }
                        return (out, i + 1);
                    }
                }
                ',' | ';' if depth == 1 => {
                    out.push((chars[from..i].iter().collect(), after_semicolon));
                    after_semicolon |= chars[i] == ';';
                    from = i + 1;
                }
                _ => {}
            }
            i += 1;
        }
        (out, chars.len())
    }

    /// The string literal (plain or raw) that is the whole of `text`, unescaped.
    fn string_literal(text: &str) -> Option<(String, &'static str)> {
        let text = text.trim();
        if let Some(rest) = text.strip_prefix('r') {
            let hashes = rest.chars().take_while(|&c| c == '#').count();
            let close = format!("\"{}", "#".repeat(hashes));
            let body = rest[hashes..]
                .strip_prefix('"')?
                .strip_suffix(close.as_str())?;
            return (!body.contains(close.as_str())).then(|| (body.to_string(), "raw"));
        }
        let body = text.strip_prefix('"')?.strip_suffix('"')?;
        let mut out = String::new();
        let mut chars = body.chars();
        while let Some(c) = chars.next() {
            match c {
                '"' => return None, // two literals
                '\\' => match chars.next()? {
                    'n' => out.push('\n'),
                    't' => out.push('\t'),
                    'r' => out.push('\r'),
                    '0' => out.push('\0'),
                    '\\' => out.push('\\'),
                    '"' => out.push('"'),
                    '\'' => out.push('\''),
                    '\n' => {
                        while chars.clone().next().is_some_and(char::is_whitespace) {
                            chars.next();
                        }
                    }
                    'u' => {
                        let code: String =
                            chars.by_ref().skip(1).take_while(|&c| c != '}').collect();
                        out.push(char::from_u32(u32::from_str_radix(&code, 16).ok()?)?);
                    }
                    'x' => {
                        let code: String = chars.by_ref().take(2).collect();
                        out.push(char::from(u8::from_str_radix(&code, 16).ok()?));
                    }
                    _ => return None,
                },
                _ => out.push(c),
            }
        }
        Some((out, "str"))
    }

    /// Stands for a format placeholder in a resolved format string.
    const PLACEHOLDER: char = '\u{E000}';

    /// How a constructor's argument resolves.
    enum Resolved {
        /// a literal text (placeholders as [`PLACEHOLDER`]) and its form
        Text(String, &'static str),
        /// `NamespacedSymbol::parse(<literal>)`: read at that inner site
        Inner,
        /// not resolved: the kind of argument
        Site(&'static str),
    }

    /// The `const` / `static` string items of the crate's non-test code: name
    /// -> its values (several when the name is defined more than once).
    fn string_constants(
        sources: &[(String, Vec<char>, Vec<bool>)],
    ) -> std::collections::BTreeMap<String, std::collections::BTreeSet<String>> {
        let item = regex::Regex::new(
            r#"(?:const|static)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:\s*&\s*(?:'static\s+)?str\s*=\s*((?:r#*)?"[^;]*)\s*;"#,
        )
        .unwrap();
        let mut out: std::collections::BTreeMap<String, std::collections::BTreeSet<String>> =
            Default::default();
        for (_, chars, _) in sources {
            let text: String = chars.iter().collect();
            for capture in item.captures_iter(&text) {
                let value = string_literal(&capture[2])
                    .map_or_else(|| format!("{PLACEHOLDER}unresolved"), |(value, _)| value);
                out.entry(capture[1].to_string()).or_default().insert(value);
            }
        }
        out
    }

    fn resolve(
        argument: &str,
        constants: &std::collections::BTreeMap<String, std::collections::BTreeSet<String>>,
    ) -> Resolved {
        let mut text = argument.trim();
        while let Some(rest) = text.strip_prefix('&') {
            text = rest.trim_start();
        }
        if let Some((value, form)) = string_literal(text) {
            return Resolved::Text(value, form);
        }
        // a literal or a format behind `.as_str()` / `.to_string()`
        for suffix in [".as_str()", ".to_string()"] {
            if let Some(inner) = text.strip_suffix(suffix)
                && let resolved @ Resolved::Text(..) = resolve(inner, constants)
            {
                return resolved;
            }
        }
        if text.starts_with("b\"") || text.starts_with("br\"") || text.starts_with("br#") {
            return Resolved::Site("bytes");
        }
        if text.starts_with('$') {
            return Resolved::Site("macro parameter");
        }
        for (prefix, form) in [("format!", "format"), ("NamespacedSymbol::parse", "inner")] {
            let Some(rest) = text.strip_prefix(prefix) else {
                continue;
            };
            let chars: Vec<char> = rest.trim_start().chars().collect();
            if chars.first() != Some(&'(') {
                return Resolved::Site("expression");
            }
            let (inner, end) = arguments(&chars, 0);
            if end != chars.len() || inner.is_empty() {
                return Resolved::Site("expression");
            }
            let Some((value, _)) = string_literal(&inner[0].0) else {
                return Resolved::Site(if form == "format" {
                    "format of a non-literal"
                } else {
                    "expression"
                });
            };
            if form == "inner" {
                return Resolved::Inner;
            }
            // {{ / }} are braces, every other {..} a placeholder
            let (mut out, mut chars) = (String::new(), value.chars().peekable());
            while let Some(c) = chars.next() {
                match c {
                    '{' if chars.peek() == Some(&'{') => {
                        chars.next();
                        out.push('{');
                    }
                    '}' if chars.peek() == Some(&'}') => {
                        chars.next();
                        out.push('}');
                    }
                    '{' => {
                        for c in chars.by_ref() {
                            if c == '}' {
                                break;
                            }
                        }
                        out.push(PLACEHOLDER);
                    }
                    c => out.push(c),
                }
            }
            return Resolved::Text(out, "format");
        }
        let path =
            regex::Regex::new(r"^(?:[A-Za-z_][A-Za-z0-9_]*::)*([A-Za-z_][A-Za-z0-9_]*)$").unwrap();
        if let Some(capture) = path.captures(text) {
            return match constants.get(&capture[1]) {
                Some(values) if values.len() == 1 => {
                    let value = values.first().unwrap().clone();
                    if value.contains(PLACEHOLDER) {
                        Resolved::Site("identifier")
                    } else {
                        Resolved::Text(value, "const")
                    }
                }
                _ => Resolved::Site("identifier"),
            };
        }
        if text.starts_with("concat!") {
            return Resolved::Site("concat");
        }
        Resolved::Site("expression")
    }

    /// The names a text names: (namespace, name), the namespace empty when the
    /// name is unqualified. A name constructor's text is one name; a parsed
    /// text is lexed for identifiers, namespaced or not, an attribute group
    /// `ns::{..}::name` skipped; numbers are skipped.
    fn names_in(text: &str, one_name: bool) -> Vec<(String, String)> {
        let chars: Vec<char> = text.trim().chars().collect();
        let word = |c: char| c.is_alphanumeric() || c == '_' || c == PLACEHOLDER;
        let mut out = Vec::new();
        let mut i = 0;
        while i < chars.len() {
            let c = chars[i];
            if c.is_numeric() && c != PLACEHOLDER {
                while i < chars.len() && word(chars[i]) {
                    i += 1;
                }
                continue;
            }
            if !word(c) {
                i += 1;
                continue;
            }
            let mut segments = Vec::new();
            loop {
                let start = i;
                while i < chars.len() && word(chars[i]) {
                    i += 1;
                }
                segments.push(chars[start..i].iter().collect::<String>());
                if chars.get(i) == Some(&':') && chars.get(i + 1) == Some(&':') {
                    let mut k = i + 2;
                    // an attribute group ::{..}::
                    if chars.get(k) == Some(&'{')
                        && let Some(close) = (k..chars.len()).find(|&j| chars[j] == '}')
                        && chars.get(close + 1) == Some(&':')
                        && chars.get(close + 2) == Some(&':')
                    {
                        k = close + 3;
                    }
                    if chars.get(k).is_some_and(|&c| word(c)) {
                        i = k;
                        continue;
                    }
                }
                break;
            }
            let local = segments.pop().unwrap();
            out.push((segments.join("::"), local));
        }
        if one_name && out.len() != 1 {
            // not one identifier: the whole text is the name
            let text: String = chars.iter().collect();
            return match text.rsplit_once("::") {
                Some((namespace, local)) => vec![(namespace.to_string(), local.to_string())],
                None => vec![(String::new(), text)],
            };
        }
        out
    }

    /// THE TRIPWIRE (fail closed), over the crate's non-test code
    /// ([`non_test_sources`]):
    /// - constructors read: the macros `vk_symbol!`, `vk_symbol_uncached!`,
    ///   `vakint_symbol!`, `symbol!` (every positional argument before a `;`),
    ///   `vk_parse!`, `vakint_parse!`, `parse!`, `try_parse!`, `parse_lit!`
    ///   (any path prefix, e.g. `symbolica::`, `$crate::`), and the functions
    ///   `NamespacedSymbol::parse` / `try_parse` / `try_parse_lit`,
    ///   `SymbolBuilder::new` and `Atom::parse`;
    /// - argument forms resolved: a string literal `".."`, a raw string
    ///   `r".."` / `r#".."#`, `format!("..", ..)` with a literal format string
    ///   (a placeholder before `::` is vakint's namespace, a placeholder inside
    ///   a name is filled by a digit), an identifier or path naming a
    ///   `const` / `static` `&str` item of the crate (defined once, by a
    ///   literal), each optionally behind `&`; `SymbolBuilder::new(
    ///   NamespacedSymbol::parse(<literal>))` is read at its inner site;
    /// - EVERY OTHER ARGUMENT (byte strings, `concat!`, expressions, macro
    ///   parameters, identifiers it cannot resolve, a format placeholder that
    ///   is a whole name or expression) fails unless its site is on
    ///   [`REVIEWED_SITES`] with its count and a reason; a stale entry fails;
    /// - the names: a name constructor's text is one name; a parsed text is
    ///   lexed for identifiers ([`names_in`]); an unqualified name is in the
    ///   constructor's default namespace (vakint's unless a literal other
    ///   namespace is given) except Symbolica's built-ins; every name in
    ///   vakint's namespace must be internal ([`internal`] or the display
    ///   table's `internal`), a data atom, a form or a native of the display
    ///   table (src/display_forms.json) or on [`REVIEWED_OUTPUT`]; an entry of
    ///   [`REVIEWED_OUTPUT`] no site creates fails.
    ///
    /// Not read: symbols made from strings by any other function or macro,
    /// the names a reviewed site creates at run time, and the names the
    /// crate's DEPENDENCIES create (Symbolica, spenso, RustRed: out of this
    /// test's view).
    #[test]
    fn listed_constructors_create_classified_names_and_other_arguments_are_reviewed() {
        crate::Vakint::initialize_vakint_symbols();
        let (sources, orphans) = non_test_sources();
        let constants = string_constants(&sources);
        let table: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(
                std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/display_forms.json"),
            )
            .unwrap(),
        )
        .unwrap();
        let listed = |key: &str| -> std::collections::BTreeSet<String> {
            table[key]
                .as_array()
                .unwrap()
                .iter()
                .map(|entry| entry["name"].as_str().unwrap().to_string())
                .collect()
        };
        let (display_internal, display_output) = (listed("internal"), {
            let mut s = listed("forms");
            s.extend(listed("native"));
            s
        });
        let builtin = symbolica::state::State::BUILTIN_NAMES_AND_ALIASES;
        let named = regex::Regex::new(r"^\s*[A-Za-z_]\w*\s*=[^=]").unwrap();
        let mut counts: std::collections::BTreeMap<(String, &str), usize> = Default::default();
        let mut sites: std::collections::BTreeMap<(String, String, String), (usize, &str, usize)> =
            Default::default();
        // name -> its first site
        let mut unclassified: std::collections::BTreeMap<String, String> = Default::default();
        let mut created = std::collections::BTreeSet::new();
        for (file, chars, literal) in &sources {
            let ident = |c: char| c.is_alphanumeric() || c == '_';
            let mut i = 0;
            while i < chars.len() {
                if literal[i]
                    || !ident(chars[i])
                    || chars[i].is_numeric()
                    || (i > 0 && ident(chars[i - 1]))
                {
                    i += 1;
                    continue;
                }
                let start = i;
                let mut segments = Vec::new();
                loop {
                    let from = i;
                    while i < chars.len() && ident(chars[i]) {
                        i += 1;
                    }
                    segments.push(chars[from..i].iter().collect::<String>());
                    if chars.get(i) == Some(&':')
                        && chars.get(i + 1) == Some(&':')
                        && chars.get(i + 2).is_some_and(|&c| ident(c))
                    {
                        i += 2;
                        continue;
                    }
                    break;
                }
                let after_path = i;
                let mut j = i;
                let bang = chars.get(j) == Some(&'!');
                if bang {
                    j += 1;
                }
                while chars.get(j).is_some_and(|c| c.is_whitespace()) {
                    j += 1;
                }
                let last = segments.last().unwrap().as_str();
                let pair = segments[segments.len().saturating_sub(2)..].join("::");
                let constructor = if bang
                    && MACRO_CONSTRUCTORS.contains(&last)
                    && matches!(chars.get(j), Some('(' | '[' | '{'))
                {
                    format!("{last}!")
                } else if !bang
                    && FUNCTION_CONSTRUCTORS.contains(&pair.as_str())
                    && chars.get(j) == Some(&'(')
                {
                    pair
                } else {
                    i = after_path;
                    continue;
                };
                let line = chars[..start].iter().filter(|&&c| c == '\n').count() + 1;
                let (args, _) = arguments(chars, j);
                let positional: Vec<&str> = args
                    .iter()
                    .filter(|(text, after_semicolon)| !after_semicolon && !named.is_match(text))
                    .map(|(text, _)| text.as_str())
                    .collect();
                let named_namespace = args.iter().find_map(|(text, _)| {
                    text.trim()
                        .strip_prefix("default_namespace")
                        .and_then(|rest| rest.trim_start().strip_prefix('='))
                        .map(str::trim)
                });
                let one_name = NAME_CONSTRUCTORS.contains(&constructor.as_str());
                let (texts, namespace_argument): (Vec<&str>, Option<&str>) =
                    if constructor == "symbol!" {
                        (positional.clone(), None)
                    } else if constructor == "parse_lit!" {
                        (vec![], None)
                    } else {
                        (
                            positional.first().copied().into_iter().collect(),
                            named_namespace.or(if one_name {
                                None
                            } else {
                                positional.get(1).copied()
                            }),
                        )
                    };
                // the default namespace: vakint's unless a literal or a
                // constant names another one
                let default_namespace = match namespace_argument.map(|a| resolve(a, &constants)) {
                    Some(Resolved::Text(namespace, _)) if !namespace.contains(PLACEHOLDER) => {
                        namespace
                    }
                    _ => "vakint".to_string(),
                };
                let mut resolved = Vec::new();
                if constructor == "parse_lit!" {
                    let (_, end) = arguments(chars, j);
                    resolved.push((
                        "<tokens>".to_string(),
                        Resolved::Text(chars[j + 1..end - 1].iter().collect(), "tokens"),
                    ));
                } else if texts.is_empty() {
                    resolved.push(("".to_string(), Resolved::Site("no argument")));
                } else {
                    for text in texts {
                        resolved.push((text.to_string(), resolve(text, &constants)));
                    }
                }
                for (argument, resolution) in resolved {
                    let site_kind = match resolution {
                        Resolved::Inner => {
                            *counts.entry((constructor.clone(), "inner")).or_default() += 1;
                            None
                        }
                        Resolved::Site(kind) => Some(kind),
                        Resolved::Text(text, form) => {
                            *counts.entry((constructor.clone(), form)).or_default() += 1;
                            let mut whole_placeholder = false;
                            for (namespace, local) in names_in(&text, one_name) {
                                let placeholder =
                                    !local.is_empty() && local.chars().all(|c| c == PLACEHOLDER);
                                let namespace = if namespace.is_empty() {
                                    // a parsed text names Symbolica's built-ins
                                    // unqualified; a name constructor may make
                                    // vakint's own symbol (a decorated name)
                                    if !one_name && builtin.contains(&local.as_str()) {
                                        continue;
                                    }
                                    default_namespace.clone()
                                } else if namespace.contains(PLACEHOLDER) {
                                    "vakint".to_string()
                                } else {
                                    namespace
                                };
                                if namespace != "vakint" {
                                    continue;
                                }
                                // a placeholder that is a whole name or
                                // expression: not resolved
                                if placeholder {
                                    whole_placeholder = true;
                                    continue;
                                }
                                let filled = local.replace(PLACEHOLDER, "1");
                                let full = format!("vakint::{filled}");
                                created.insert(filled.clone());
                                if !(internal(&filled, false)
                                    || display_internal.contains(&full)
                                    || data_atom(&full)
                                    || display_output.contains(&full)
                                    || REVIEWED_OUTPUT.contains(&filled.as_str()))
                                {
                                    unclassified.entry(filled.clone()).or_insert_with(|| {
                                        format!("{full} ({file}:{line}, {constructor})")
                                    });
                                }
                            }
                            whole_placeholder.then_some("format placeholder")
                        }
                    };
                    if let Some(kind) = site_kind {
                        let argument = argument.split_whitespace().collect::<Vec<_>>().join(" ");
                        let entry = sites
                            .entry((file.clone(), constructor.clone(), argument))
                            .or_insert((0, kind, line));
                        entry.0 += 1;
                        *counts
                            .entry((constructor.clone(), "REVIEWED SITE"))
                            .or_default() += 1;
                    }
                }
                i = after_path;
            }
        }
        let mut report = String::from("constructor, form: count\n");
        for ((constructor, form), count) in &counts {
            report.push_str(&format!("  {constructor}, {form}: {count}\n"));
        }
        report.push_str(&format!(
            "distinct vakint names created by resolved arguments: {}\nreviewed sites: {} distinct\n",
            created.len(),
            sites.len()
        ));
        eprintln!("{report}");
        let mut failures = Vec::new();
        for ((file, constructor, argument), (count, kind, line)) in &sites {
            match REVIEWED_SITES
                .iter()
                .find(|(f, c, a, _, _)| f == file && c == constructor && a == argument)
            {
                None => failures.push(format!(
                    "UNREVIEWED ({kind}, line {line}): (\"{file}\", \"{constructor}\", {argument:?}, {count}, \"?\"),"
                )),
                Some((_, _, _, reviewed, _)) if reviewed != count => failures.push(format!(
                    "COUNT {reviewed} -> {count}: ({file}, {constructor}, {argument:?})"
                )),
                Some((_, _, _, _, reason)) => {
                    if !["settings", "data", "caller", "format", "macro"].contains(reason) {
                        failures.push(format!("REASON {reason:?}: ({file}, {constructor}, {argument:?})"));
                    }
                }
            }
        }
        for (file, constructor, argument, _, _) in REVIEWED_SITES {
            if !sites.contains_key(&(
                file.to_string(),
                constructor.to_string(),
                argument.to_string(),
            )) {
                failures.push(format!("STALE SITE: ({file}, {constructor}, {argument:?})"));
            }
        }
        for orphan in &orphans {
            if !REVIEWED_ORPHANS.contains(&orphan.as_str()) {
                failures.push(format!("ORPHAN (not compiled, not read): {orphan:?}"));
            }
        }
        for orphan in REVIEWED_ORPHANS {
            if !orphans.iter().any(|o| o == orphan) {
                failures.push(format!("STALE ORPHAN: {orphan}"));
            }
        }
        for name in REVIEWED_OUTPUT {
            if !created.contains(name) {
                failures.push(format!("STALE REVIEWED_OUTPUT: {name}"));
            }
        }
        for name in unclassified.keys() {
            if !UNDECIDED.contains(&name.as_str()) {
                failures.push(format!("UNCLASSIFIED: {}", unclassified[name]));
            }
        }
        for name in UNDECIDED {
            if !unclassified.contains_key(*name) {
                failures.push(format!("STALE UNDECIDED: {name}"));
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
        assert!(
            UNDECIDED.is_empty(),
            "{} names in vakint's namespace are created by literal arguments and are not yet classified (UNDECIDED): {}",
            UNDECIDED.len(),
            UNDECIDED.join(", ")
        );
    }
}
