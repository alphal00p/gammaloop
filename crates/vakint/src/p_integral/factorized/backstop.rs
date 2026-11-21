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
//!   every symbol the crate's code creates and fails when a vakint name is
//!   neither in the internal table, nor a data atom, nor on the reviewed output
//!   list (`REVIEWED_OUTPUT`, test code);
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
//! coefficient, a missing order); an internal symbol vakint creates without
//! a literal or format string in its code (the tripwire reads source text);
//! an input symbol that should have been substituted (it is an input); a zero
//! that is no number (a zero base written as an unexpanded expression).

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;

use super::{FactorizedError, FactorizedKind, FactorizedNumerator, Result, names};
use crate::VakintSettings;
use crate::utils::vakint_macros::vk_symbol;

/// vakint's internal names (without the `vakint::` namespace), exact.
const INTERNAL_EXACT: [&str; 35] = [
    // the factorised entries
    "closure_lambda",
    "slash",
    // vakint's symbol table and its input / analysis placeholders
    "ERROR",
    "dot_dummy_ind",
    "dot_pow",
    "input_coefficient",
    "partial_numerical_input",
    "rustred_input",
    "normalization_input",
    "tensor_scalar_coefficient",
    "feynkit_external_vector",
    "feynkit_opaque_index",
    "momentum_match_parameter",
    "coordinate_limit_spectator",
    "coordinate_metadata",
    "VakintLambdaScalingAnalysis",
    "RustRedMaterializeSpectator",
    "x",
    "y",
    "xa",
    "ya",
    "vec",
    "vec1",
    "vxs",
    "mom",
    "n_loops",
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
];

/// Families of internal names made with format strings: a prefix and a
/// suffix around the formatted part (the code's pattern in the comment).
const INTERNAL_FAMILIES: [(&str, &str); 13] = [
    ("closure_tensor_", ""),                // closure_tensor_{n}
    ("graph_entry_branch_", ""),            // graph_entry_branch_{i}
    ("input_coefficient_placeholder_", ""), // input_coefficient_placeholder_{index}
    ("coordinate_metadata_", ""),           // coordinate_metadata_{i}
    ("momentum_match_source_", ""),         // momentum_match_source_{axis}
    ("vakint_master_transport_", ""),       // vakint_master_transport_{axis}
    ("rustred_family_", ""),                // rustred_family_{:x}
    ("I5L_", ""),                           // I5L_{}
    ("krotated", ""),                       // krotated{}
    ("z", ""),                              // z{n} (with digits only, see internal)
    ("k", ""),                              // k{} (with digits only, see internal)
    ("", "_pinch_"),                        // {}_pinch_{} (contains)
    ("msq", ""),                            // msq, msq{}_
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
const REVIEWED_OUTPUT: [&str; 58] = [
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
    "GammaSeriesTail",
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
    "{}",
    "{}{}",
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
    let digits = |rest: &str| !rest.is_empty() && rest.chars().all(|c| c.is_ascii_digit());
    name.ends_with('_')
        || INTERNAL_EXACT.contains(&name)
        || (!vakint_notation && NOTATION_HEADS.contains(&name))
        || INTERNAL_FAMILIES
            .iter()
            .any(|(prefix, contains)| match (*prefix, *contains) {
                ("z", _) | ("k", _) => name.strip_prefix(prefix).is_some_and(digits),
                ("", inner) => name.contains(inner),
                (prefix, _) => name.starts_with(prefix),
            })
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

/// The symbols of every atom of `value`, of `prefactor` and of the settings.
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
        if data_atom(&name) {
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
            "UFO::gs^2*probe::nowhere",
            "the symbol probe::nowhere, which is neither",
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

    /// The deny table: the families made with format strings and the
    /// wildcards are internal; vakint's output vocabulary, including the
    /// massive-vacuum constants and the RustRed tags, passes a marked result.
    #[test]
    fn internal_names_are_refused_and_the_output_vocabulary_passes() {
        for name in [
            "closure_tensor_12",
            "graph_entry_branch_3",
            "input_coefficient_placeholder_0",
            "z7",
            "k3",
            "pow_",
            "idx___",
            "topo",
            "dot_dummy_ind",
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

    /// `text` with every item under a line-leading `#[cfg(test)]` blanked: the
    /// attribute up to the first `;` or `,` outside brackets (a `use`, a
    /// `mod x;`, a field) or through the brace block that comes first (a
    /// module, a function, an impl, a statement). Strings, raw strings,
    /// character literals and comments are skipped when counting, so a brace
    /// in them cuts nothing more; any doubt leaves text in (scanned, never
    /// hidden).
    fn without_test_items(text: &str) -> String {
        let chars: Vec<char> = text.chars().collect();
        let attribute = regex::Regex::new(r"(?m)^[ \t]*#\[cfg\(test\)\]").unwrap();
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
        // the end (exclusive) of the lexical token starting at i, if it is a
        // string, a raw string, a character literal or a comment
        let skip = |i: usize| -> Option<usize> {
            let at = |k: usize| chars.get(k).copied();
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
                'r' if i == 0 || !(chars[i - 1].is_alphanumeric() || chars[i - 1] == '_') => {
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
                '\'' if at(i + 1) == Some('\\') => Some(
                    (i + 2..chars.len())
                        .find(|&k| chars[k] == '\'')
                        .map_or(chars.len(), |k| k + 1),
                ),
                '\'' if at(i + 2) == Some('\'') => Some(i + 3),
                _ => None,
            }
        };
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

    /// THE TRIPWIRE: every name the crate's non-test code creates in vakint's
    /// namespace with a literal or a format string (vk_symbol!, symbol!,
    /// vakint_symbol!) is either internal ([`internal`], with the format
    /// placeholders filled by a digit), a data atom or on [`REVIEWED_OUTPUT`]: a new name
    /// fails here until it is classified. Skipped: only test code (files under
    /// a tests/ directory, tests.rs, *_tests.rs, and every item under a
    /// `#[cfg(test)]`, by [`without_test_items`]); a source file with "test"
    /// elsewhere in its name is read.
    #[test]
    fn every_name_the_code_creates_is_classified() {
        let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src");
        let creating = regex::Regex::new(
            r#"(?:vk_symbol!|symbol!|vakint_symbol!)\(\s*&?(?:format!\()?\s*"([^"]*)""#,
        )
        .unwrap();
        let placeholder = regex::Regex::new(r"\{[^}]*\}").unwrap();
        let mut unclassified = Vec::new();
        let mut files = vec![root.clone()];
        let mut sources = Vec::new();
        while let Some(path) = files.pop() {
            for entry in std::fs::read_dir(&path).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    files.push(path);
                } else if path.extension().is_some_and(|e| e == "rs") {
                    sources.push(path);
                }
            }
        }
        for path in sources {
            let name = path.to_string_lossy().to_string();
            // only real test code is skipped: files under a tests/ directory,
            // tests.rs and *_tests.rs, and inline #[cfg(test)] modules
            let file = path.file_name().unwrap().to_string_lossy().to_string();
            if path.components().any(|c| c.as_os_str() == "tests")
                || file == "tests.rs"
                || file.ends_with("_tests.rs")
            {
                continue;
            }
            let text = without_test_items(&std::fs::read_to_string(&path).unwrap());
            for capture in creating.captures_iter(&text) {
                let created = &capture[1];
                let local = if let Some(rest) = created.strip_prefix("vakint::") {
                    rest
                } else if let Some(rest) = created.strip_prefix("{}::") {
                    rest
                } else if created.contains("::") {
                    continue; // another namespace: an input or a caller head
                } else {
                    created
                };
                // a format placeholder stands for a formatted part
                let filled = placeholder.replace_all(local, "1");
                if !(internal(&filled, false)
                    || data_atom(&format!("vakint::{local}"))
                    || REVIEWED_OUTPUT.contains(&local))
                {
                    unclassified.push(format!("{created} ({name})"));
                }
            }
        }
        assert!(
            unclassified.is_empty(),
            "names created in vakint's namespace that are neither internal, nor a data atom, nor reviewed output:\n{}",
            unclassified.join("\n")
        );
    }
}
