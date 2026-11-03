//! Values of the p-masters `PMaster(name)`, parsed once from the vendored
//! Forcer `Masters.prc` (the expansions of R.N. Lee, A.V. Smirnov and
//! V.A. Smirnov, arXiv:1108.0732, converted by the Forcer authors to the
//! Forcer master basis and to the MZV datamine basis).
//!
//! Each entry is the Laurent expansion in the G-scheme `ε` known to the file,
//! followed by the unknown tail `ε^K·Oep(K, PMaster(name))` where Forcer
//! writes `ep^K*(ALARMname)`. The constants `z2..z12` become `zeta(n)`; the
//! datamine products `z5z3`, `z7z3`, `z5z3z3`, `z9z3` and `z6z4z1z1` are the
//! multiple zeta values `MZV(5,3)`, ..., `MZV(6,4,1,1)`.

use std::collections::BTreeMap;
use std::sync::LazyLock;

use regex::Regex;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::id::Replacement;
use symbolica::transcendental::zeta;

use crate::VakintError;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};

static MASTERS_PRC: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/form_src/forcer/forcer/Masters.prc"
));

/// The multiple zeta values written as products in the datamine basis.
const MZV_NAMES: [(&str, &[i64]); 5] = [
    ("z5z3", &[5, 3]),
    ("z7z3", &[7, 3]),
    ("z5z3z3", &[5, 3, 3]),
    ("z9z3", &[9, 3]),
    ("z6z4z1z1", &[6, 4, 1, 1]),
];

/// The symbol of the expansion variable inside the table, `vakint::ep`.
pub(crate) fn table_epsilon() -> Symbol {
    vk_symbol!("ep")
}

/// `MZV(n_1, ..., n_k)` (see `inspection::constants`).
pub(crate) fn mzv(indices: &[i64]) -> Atom {
    vk_symbol!("MZV").call(&indices.iter().map(|&n| Atom::num(n)).collect::<Vec<_>>())
}

fn failure(detail: impl Into<String>) -> VakintError {
    VakintError::EvaluationError(format!("Forcer Masters.prc: {}", detail.into()))
}

fn parse_table() -> Result<BTreeMap<String, Atom>, VakintError> {
    crate::Vakint::initialize_vakint_symbols();
    let masters = MASTERS_PRC
        .split(".sort:Masters-1;")
        .next()
        .ok_or_else(|| failure("no master section"))?;
    let entry = Regex::new(r"id\s+Master\((\w+)\)\s*=([^;]*);").unwrap();
    let alarm = Regex::new(r"ep\^(\d+)\*\(ALARM(\w+)\)").unwrap();
    let mut constants = (2..=12)
        .map(|n| {
            Replacement::new(
                Atom::var(vk_symbol!(format!("z{n}"))).to_pattern(),
                zeta().call(Atom::num(n)),
            )
        })
        .collect::<Vec<_>>();
    constants.extend(MZV_NAMES.iter().map(|(name, indices)| {
        Replacement::new(Atom::var(vk_symbol!(name)).to_pattern(), mzv(indices))
    }));
    let mut table = BTreeMap::new();
    for captures in entry.captures_iter(masters) {
        let name = captures[1].to_string();
        let body = captures[2]
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>();
        let mut tails = 0;
        let body = alarm.replace_all(&body, |tail: &regex::Captures| {
            tails += 1;
            if tail[2] != name {
                tails += 1;
            }
            format!("ep^{0}*Oep({0},PMaster({1}))", &tail[1], name)
        });
        if tails != 1 {
            return Err(failure(format!(
                "Master({name}) needs exactly one own ALARM tail"
            )));
        }
        let expansion = vk_parse!(&body)
            .map_err(|error| failure(format!("Master({name}): {error}")))?
            .replace_multiple(&constants);
        let epsilon = table_epsilon();
        let mut unknown = None;
        expansion.visitor(&mut |part| match part {
            // Symbolica writes ζ(2n) through π^(2n).
            AtomView::Var(variable) if ![epsilon, Symbol::PI].contains(&variable.get_symbol()) => {
                unknown = Some(part.to_owned());
                false
            }
            AtomView::Fun(function)
                if [
                    zeta(),
                    vk_symbol!("MZV"),
                    crate::truncation_guard::unknown_tail_symbol(),
                ]
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
        if let Some(part) = unknown {
            return Err(failure(format!("Master({name}): unknown constant {part}")));
        }
        if table.insert(name.clone(), expansion).is_some() {
            return Err(failure(format!("Master({name}) is defined twice")));
        }
    }
    if table.is_empty() {
        return Err(failure("no Master(...) entries"));
    }
    Ok(table)
}

/// The p-master table, parsed once.
pub(crate) fn p_master_table() -> Result<&'static BTreeMap<String, Atom>, VakintError> {
    static TABLE: LazyLock<Result<BTreeMap<String, Atom>, String>> =
        LazyLock::new(|| parse_table().map_err(|error| error.to_string()));
    TABLE
        .as_ref()
        .map_err(|error| VakintError::EvaluationError(error.clone()))
}

/// The expansion of `PMaster(name)` in `epsilon`, with its unknown tail.
pub(crate) fn p_master_expansion(name: &str, epsilon: Symbol) -> Result<Atom, VakintError> {
    let expansion = p_master_table()?.get(name).ok_or_else(|| {
        VakintError::EvaluationError(format!(
            "no value is known for the p-master {name} (vendored Forcer Masters.prc)"
        ))
    })?;
    Ok(expansion
        .replace(Atom::var(table_epsilon()).to_pattern())
        .with(Atom::var(epsilon).to_pattern()))
}
