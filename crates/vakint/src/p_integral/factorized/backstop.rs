//! THE RESULT CHECK of a marked value (the last line behind validation): after
//! the evaluation, before an entry returns, the value may contain only
//! - the symbols of its INPUTS (every atom of the record, the prefactor
//!   argument, the settings' epsilon and scale symbols);
//! - the dimension only as the dimension argument of a representation slot
//!   (`mink(D, ..)`, `bis(D, ..)`, the compact `mink(D)` of a momentum);
//! - the caller notation's structural heads (`spenso::g`, `spenso::dot`,
//!   `spenso::mink`, `spenso::bis`, `spenso::gamma`, `gammalooprs::Q`), and
//!   where the entry returns vakint notation (the direct entry's value; the
//!   vacuum entry's external vectors without a leg) vakint's `p`, `k`,
//!   `dot`, `g`;
//! - Symbolica's own symbols (`symbolica::…`: 𝜋, 𝑒, γ, log, zeta, gamma, ...);
//! - vakint's output constants (`HR`, `HI`, `GPLReal`, `GPLImag`,
//!   `PolylogReal`, `PolylogImag`, `Zeta`, `Sqrt`, `EulerGamma`, `MZV`) and in
//!   Retain mode its master tags `PMaster`, `PG`, `PGamma`, `PResidual`, with
//!   their arguments; `UNKNOWN(master, k)` when the settings
//!   allow unknown master coefficients.
//!
//! Anything else (vakint's internal placeholders and parameters, e.g. the
//! Taylor parameter `vakint::closure_lambda`, the slot placeholders
//! `vakint::closure_tensor_<n>`, a dimension outside a slot, a symbol that
//! came from nowhere) is an internal failure, `Evaluation`, detail prefix
//! `result check: `, never a value.

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

use super::{FactorizedError, FactorizedKind, FactorizedNumerator, Result, names};
use crate::VakintSettings;
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;

/// vakint's own output constants (the transcendental functions and constants
/// of its values, crate::symbols: HR, HI, GPLReal, ...), with their arguments.
const VAKINT_CONSTANTS: [&str; 10] = [
    "vakint::HR",
    "vakint::HI",
    "vakint::GPLReal",
    "vakint::GPLImag",
    "vakint::PolylogReal",
    "vakint::PolylogImag",
    "vakint::Zeta",
    "vakint::Sqrt",
    "vakint::EulerGamma",
    "vakint::MZV",
];

/// Retain mode's master tags: their arguments are master names (vakint's).
const RETAIN_TAGS: [&str; 4] = [
    "vakint::PMaster",
    "vakint::PG",
    "vakint::PGamma",
    "vakint::PResidual",
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
        if known.contains(&symbol)
            || name.starts_with("symbolica::")
            || CALLER_HEADS.contains(&name)
            || VAKINT_CONSTANTS.contains(&name)
            || (vakint_notation
                && [S.p, S.k, S.dot, S.g]
                    .iter()
                    .any(|head| head.get_name() == name))
        {
            return None;
        }
        if name.starts_with("vakint::") {
            Some(format!("vakint's internal symbol {name}"))
        } else {
            Some(format!(
                "the symbol {name}, which is neither an input nor a known output"
            ))
        }
    };
    match view {
        AtomView::Num(_) => None,
        AtomView::Var(v) => name_ok(v.get_symbol()),
        AtomView::Fun(f) => {
            let symbol = f.get_symbol();
            let name = symbol.get_name();
            if RETAIN_TAGS.contains(&name)
                || VAKINT_CONSTANTS.contains(&name)
                || (allow_unknown
                    && symbol == crate::truncation_guard::unknown_coefficient_symbol())
            {
                return None;
            }
            if let Some(found) = name_ok(symbol) {
                return Some(found);
            }
            let rep = super::is_rep_name(name);
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

/// THE RESULT CHECK (module documentation) of the value `result` an entry is
/// about to return for the MARKED value `value` (`vakint_notation`: the
/// direct entry's vakint-notation value). A no-op for an unmarked value.
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
            "vakint's internal symbol vakint::UNKNOWN",
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
}
