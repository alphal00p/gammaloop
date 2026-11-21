//! THE shared helper of the tests that read the R* producer's fixture records
//! (value files with a `value` and a `valid_through` line). It decides what
//! "identical" means, so every such test calls it:
//!
//! - LOAD OR FAIL: external data comes from ONE environment variable,
//!   [`FIXTURE_ROOT_ENV`]; a test that needs it is `#[ignore]` (a plain run
//!   does not claim it) and, when run, a missing variable, directory or file
//!   or an empty record set is a FAILURE, never a silent pass;
//! - COMPARE: [`assert_equal_through`] compares the Laurent coefficients in ε
//!   exactly through `valid_through` and no further, with a SIGNED bound (a
//!   negative order is never cast to an unsigned one), and asserts that the
//!   reference is NONZERO through that order (a zero or truncated fixture
//!   cannot pass against a zero entry);
//! - COUNT: [`assert_compared`] asserts the number of records compared
//!   against the expected count.
//!
//! It depends on symbolica and std only, so that the unit tests of the crate
//! (`#[path]`) and the integration tests (`mod support`) share this one file.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};

/// The ONE environment variable of the external fixture root (the R*
/// producer's `wh-fixtures` directory: `k1-quick/`, `k5-stress/`, `k1-l5/`).
pub const FIXTURE_ROOT_ENV: &str = "VAKINT_TEST_WH_FIXTURES";

/// The external fixture root: [`FIXTURE_ROOT_ENV`], an existing directory;
/// panics otherwise (the caller is an `#[ignore]` test that was asked to run).
pub fn fixture_root() -> PathBuf {
    let Some(root) = std::env::var_os(FIXTURE_ROOT_ENV) else {
        panic!("{FIXTURE_ROOT_ENV} is not set: this test needs the external fixtures");
    };
    let root = PathBuf::from(root);
    assert!(
        root.is_dir(),
        "{FIXTURE_ROOT_ENV}={}: not a directory",
        root.display()
    );
    root
}

/// `<root>/<sub>`, an existing directory, else a failure.
pub fn fixture_dir(sub: &str) -> PathBuf {
    let dir = fixture_root().join(sub);
    assert!(dir.is_dir(), "missing fixture directory {}", dir.display());
    dir
}

/// A required environment variable (a path), else a failure.
pub fn required_env(name: &str, why: &str) -> PathBuf {
    match std::env::var_os(name) {
        Some(value) => PathBuf::from(value),
        None => panic!("{name} is not set: {why}"),
    }
}

/// Every record directory under `dir` (at depth 1 or 2) holding both
/// `input.txt` and `value.txt`, sorted; an empty set is a failure.
pub fn record_dirs(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let children = |d: &Path| -> Vec<PathBuf> {
        std::fs::read_dir(d)
            .unwrap_or_else(|e| panic!("{}: {e}", d.display()))
            .map(|e| e.unwrap().path())
            .filter(|p| p.is_dir())
            .collect()
    };
    let is_record = |p: &Path| p.join("input.txt").is_file() && p.join("value.txt").is_file();
    for child in children(dir) {
        if is_record(&child) {
            out.push(child);
        } else {
            out.extend(children(&child).into_iter().filter(|p| is_record(p)));
        }
    }
    out.sort();
    assert!(!out.is_empty(), "no record under {}", dir.display());
    out
}

/// A file's text, else a failure.
pub fn read_fixture(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The `value` line of a value file.
pub fn value_line(text: &str) -> &str {
    text.lines()
        .find_map(|l| l.strip_prefix("value\t"))
        .expect("a value line")
}

/// The order in ε through which a value file's value is valid (`valid_through
/// eps^n`, n signed).
pub fn valid_through(text: &str) -> i64 {
    let line = text
        .lines()
        .find_map(|l| l.strip_prefix("valid_through\t"))
        .expect("a valid_through line");
    line.trim_start_matches("eps^")
        .split_whitespace()
        .next()
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("valid_through {line:?}: not eps^<signed integer>"))
}

/// The Laurent coefficients of `value` in `epsilon` by SIGNED order, read off
/// the expanded sum term by term (no series call: no order bound to cast).
/// Fails if `epsilon` survives inside a coefficient (a non-polynomial
/// dependence the term-wise reading would misfile).
pub fn epsilon_coefficients(value: &Atom, epsilon: Symbol) -> BTreeMap<i64, Atom> {
    let expanded = value.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
        _ if expanded.is_zero() => Vec::new(),
        _ => vec![expanded.clone()],
    };
    let mut coefficients = BTreeMap::<i64, Atom>::new();
    for term in terms {
        let factors = match term.as_view() {
            AtomView::Mul(product) => product.iter().map(|f| f.to_owned()).collect::<Vec<_>>(),
            _ => vec![term.clone()],
        };
        let mut order = 0i64;
        let mut rest = Atom::num(1);
        for factor in factors {
            match factor.as_view() {
                AtomView::Var(v) if v.get_symbol() == epsilon => order += 1,
                AtomView::Pow(p) if matches!(p.get_base_exp().0, AtomView::Var(v) if v.get_symbol() == epsilon) =>
                {
                    let exponent = p.get_base_exp().1;
                    order += i64::try_from(exponent)
                        .unwrap_or_else(|_| panic!("a non-integer power of {epsilon:?}: {factor}"));
                }
                _ => rest *= factor,
            }
        }
        assert!(
            !rest.contains_symbol(epsilon),
            "{epsilon:?} inside a coefficient: {term}"
        );
        let entry = coefficients.entry(order).or_insert(Atom::Zero);
        *entry = &*entry + rest;
    }
    coefficients
        .into_iter()
        .map(|(n, c)| (n, c.expand()))
        .filter(|(_, c)| !c.is_zero())
        .collect()
}

/// EXACT equality of `got` and `expected` coefficient by coefficient in
/// `epsilon` from the lowest order of either through `order` (SIGNED) and NO
/// FURTHER (the fixtures' identity rule), and a NONZERO reference through
/// `order`. Fails at the first differing order (fast: no series expansion).
/// Returns the number of orders compared.
pub fn assert_equal_through(
    got: &Atom,
    expected: &Atom,
    epsilon: Symbol,
    order: i64,
    what: &str,
) -> usize {
    let got_c = epsilon_coefficients(got, epsilon);
    let expected_c = epsilon_coefficients(expected, epsilon);
    assert!(
        expected_c.keys().any(|n| *n <= order),
        "{what}: the reference is zero through eps^{order}: {expected}"
    );
    let first = got_c
        .keys()
        .chain(expected_c.keys())
        .copied()
        .min()
        .expect("a nonzero reference");
    let mut compared = 0;
    for n in first..=order {
        let g = got_c.get(&n).cloned().unwrap_or(Atom::Zero);
        let e = expected_c.get(&n).cloned().unwrap_or(Atom::Zero);
        let difference = (&g - &e).expand();
        assert!(
            difference.is_zero(),
            "{what}: the coefficient of eps^{n} differs\nentry     {g}\nreference {e}\ndifference {difference}"
        );
        compared += 1;
    }
    compared
}

/// Whether a value file DECLARES its value to be exactly zero (the `value`
/// line is the atom `0`, as the producer writes a site whose boundary
/// patterns all integrate to zero): such a fixture is compared with
/// [`assert_declared_zero`], every other one with [`assert_equal_through`]
/// (whose nonzero-reference rule then catches a value that parsed to zero).
pub fn is_declared_zero(text: &str) -> bool {
    value_line(text).trim() == "0"
}

/// `got` is zero through `order` (SIGNED): every coefficient of an order at
/// or below `order` vanishes exactly. Returns the number of orders checked.
pub fn assert_declared_zero(got: &Atom, epsilon: Symbol, order: i64, what: &str) -> usize {
    let got_c = epsilon_coefficients(got, epsilon);
    let mut checked = 0;
    for (n, c) in &got_c {
        if *n <= order {
            panic!("{what}: declared zero, but the coefficient of eps^{n} is {c}");
        }
        checked += 1;
    }
    checked
}

/// The number of records compared equals the expected count (a filter that
/// selects a subset passes `None`: then at least one).
pub fn assert_compared(what: &str, compared: usize, expected: Option<usize>) {
    match expected {
        Some(n) => assert_eq!(
            compared, n,
            "{what}: {compared} records compared, {n} expected"
        ),
        None => assert!(compared > 0, "{what}: the filter selected no record"),
    }
    eprintln!(
        "{what}: COMPARED {compared} records{}",
        if expected.is_some() {
            ""
        } else {
            " (FILTERED: count not asserted)"
        }
    );
}
