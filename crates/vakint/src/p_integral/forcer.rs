//! The Forcer backend: massless propagator integrals reduced by the vendored
//! Forcer FORM program (B. Ruijl, T. Ueda, J.A.M. Vermaseren, arXiv:1704.06650).
//!
//! Forcer receives the contracted canonical graph (positive-power lines only)
//! and returns the exact reduction in its G-scheme master basis,
//! `rat(n,d) × GschemeConstants(x1,x2) × Master(name)`. Forcer sets `Q.Q = 1`,
//! so the numerator is split into homogeneous parts before the run; every
//! part becomes one [`ReducedPropagatorValue`]. Normalization, expansion and
//! master substitution are left to the shared p-finalizer (`finalize.rs`).

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{LazyLock, Mutex};

use colored::Colorize;
use log::{debug, warn};
use regex::Regex;
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;
use symbolica::function;
use version_compare::{Cmp, compare_to};

use super::ForcerOptions;
use super::{
    FORCER_UPSTREAM_COMMIT, MINIMAL_FORCER_FORM_VERSION, PFamilyData, ReducedPropagatorValue,
    external_head,
};
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    ReplacementRules, TEMPLATES, Vakint, VakintError, VakintSettings, get_integer_from_atom,
    get_prop_with_id,
};

mod vendored {
    include!(concat!(env!("OUT_DIR"), "/vakint_forcer_sources.rs"));
}

/// Name of the directory shared by the processes of the user `owner`.
fn shared_directory_name(owner: u32) -> String {
    format!("vakint-forcer-{}-{owner}", &FORCER_UPSTREAM_COMMIT[..7])
}

/// The embedded Forcer sources, decompressed: relative path and bytes.
fn vendored_sources() -> Result<Vec<(&'static str, Vec<u8>)>, VakintError> {
    vendored::FORCER_FILES
        .iter()
        .map(|(name, compressed)| {
            let mut contents = Vec::new();
            std::io::Read::read_to_end(
                &mut flate2::read::GzDecoder::new(*compressed),
                &mut contents,
            )?;
            Ok((*name, contents))
        })
        .collect()
}

/// The owner of `path` (not following links) when it is a directory closed
/// to everyone else (no group or other permission bits).
#[cfg(unix)]
fn private_directory_owner(path: &Path) -> Option<u32> {
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    let metadata = fs::symlink_metadata(path).ok()?;
    (metadata.is_dir() && metadata.permissions().mode() & 0o077 == 0).then(|| metadata.uid())
}

#[cfg(not(unix))]
fn private_directory_owner(_path: &Path) -> Option<u32> {
    None
}

#[cfg(unix)]
fn make_private(path: &Path) -> std::io::Result<()> {
    use std::os::unix::fs::PermissionsExt;
    fs::set_permissions(path, fs::Permissions::from_mode(0o700))
}

#[cfg(not(unix))]
fn make_private(_path: &Path) -> std::io::Result<()> {
    Ok(())
}

/// Whether `directory` can be handed to FORM: a private directory of `owner`
/// holding exactly the embedded sources as regular files.
fn is_trusted(directory: &Path, owner: u32, sources: &[(&str, Vec<u8>)]) -> bool {
    private_directory_owner(directory) == Some(owner)
        && sources.iter().all(|(name, contents)| {
            let path = directory.join(name);
            fs::symlink_metadata(&path).is_ok_and(|metadata| metadata.is_file())
                && fs::read(&path).is_ok_and(|found| found == *contents)
        })
}

/// The vendored sources under `base` (`temporary_directory`, else the system
/// temporary directory), materialized once per process and base.
///
/// FORM includes `forcer.h` from this directory, so it is trusted only when
/// it is a private directory (mode 0700, not a link) of the current user
/// holding byte for byte the embedded sources. Such a directory is shared by
/// the processes of one user. Anything else at that path, e.g. a directory
/// planted by another user of a shared `/tmp` or a damaged copy, is never
/// used: the sources then go to a fresh private directory of this process.
fn materialized_sources(base: Option<&str>) -> Result<PathBuf, VakintError> {
    static DIRECTORIES: LazyLock<Mutex<HashMap<PathBuf, PathBuf>>> =
        LazyLock::new(Default::default);
    let root = base.map(PathBuf::from).unwrap_or_else(std::env::temp_dir);
    let mut directories = DIRECTORIES
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if let Some(directory) = directories.get(&root) {
        return Ok(directory.clone());
    }
    let sources = vendored_sources()?;
    let staging = crate::create_unique_tmp_run_dir(Some(&*root.to_string_lossy()), ".forcer")?;
    make_private(&staging)?;
    let directory = match private_directory_owner(&staging) {
        Some(owner) => {
            let shared = root.join(shared_directory_name(owner));
            if is_trusted(&shared, owner, &sources) {
                fs::remove_dir(&staging)?;
                shared
            } else {
                write_sources(&staging, &sources)?;
                // Atomic publication; a concurrent process may have won the race.
                let published =
                    fs::symlink_metadata(&shared).is_err() && fs::rename(&staging, &shared).is_ok();
                if published {
                    shared
                } else if is_trusted(&shared, owner, &sources) {
                    fs::remove_dir_all(&staging)?;
                    shared
                } else {
                    warn!(
                        "Ignoring the untrusted Forcer source directory '{}'; using '{}' instead.",
                        shared.display(),
                        staging.display()
                    );
                    staging
                }
            }
        }
        None => {
            write_sources(&staging, &sources)?;
            staging
        }
    };
    directories.insert(root, directory.clone());
    Ok(directory)
}

fn write_sources(directory: &Path, sources: &[(&str, Vec<u8>)]) -> Result<(), VakintError> {
    for (name, contents) in sources {
        let path = directory.join(name);
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, contents)?;
    }
    Ok(())
}

/// The FORM version at `form_exe_path` when it can run Forcer, otherwise the
/// reason why not, phrased to complete "'{form_exe_path}' ...".
pub(crate) fn forcer_form_version(form_exe_path: &str) -> Result<String, String> {
    let Some(output) = crate::probe_command_output(form_exe_path, &["-version"]) else {
        return Err("could not be run".into());
    };
    let Some(version) = Regex::new(r"FORM\s([\.\d]+)")
        .unwrap()
        .captures(&output)
        .map(|captures| captures[1].to_owned())
    else {
        return Err(format!(
            "reported no FORM version ('{}')",
            output.lines().next().unwrap_or_default().trim()
        ));
    };
    match compare_to(&version, MINIMAL_FORCER_FORM_VERSION, Cmp::Ge) {
        Ok(true) => Ok(version),
        _ => Err(format!("is FORM {version}")),
    }
}

/// One canonical line: endpoints and the integer coefficients of
/// `k(1..=L)` and of the canonical `p(1)` in its momentum.
#[derive(Clone, Debug, PartialEq, Eq)]
struct CanonicalLine {
    from: i64,
    to: i64,
    loops: Vec<i64>,
    external: i64,
}

impl CanonicalLine {
    fn momentum(&self, external: &Atom) -> Atom {
        let mut momentum = external * Atom::num(self.external);
        for (index, coefficient) in self.loops.iter().enumerate() {
            momentum += function!(S.k, Atom::num(index as i64 + 1)) * Atom::num(*coefficient);
        }
        momentum
    }
}

fn canonical_lines(family: &PFamilyData) -> Result<Vec<CanonicalLine>, VakintError> {
    let invalid = |detail: String| {
        VakintError::InvalidGenericExpression(format!("{}: {detail}", family.name))
    };
    Vakint::initialize_vakint_symbols();
    let expression = vk_parse!(family.lines).map_err(|error| invalid(error.to_string()))?;
    let heads = (1..=family.loops)
        .map(|i| function!(S.k, Atom::num(i as i64)))
        .chain([function!(S.p, Atom::num(1))])
        .collect::<Vec<_>>();
    (1..=family.line_count())
        .map(|id| {
            let prop = get_prop_with_id(expression.as_view(), id)
                .ok_or_else(|| invalid(format!("missing line {id}")))?;
            let node = |name: &str| {
                get_integer_from_atom(prop[&vk_symbol!(name)].as_view())
                    .ok_or_else(|| invalid(format!("line {id} has a non-integer vertex")))
            };
            let momentum = &prop[&vk_symbol!("q_")];
            let coefficients = heads
                .iter()
                .map(|head| {
                    let mut value = momentum.clone();
                    for other in &heads {
                        let image = if other == head {
                            Atom::num(1)
                        } else {
                            Atom::Zero
                        };
                        value = value.replace(other.to_pattern()).with(image.to_pattern());
                    }
                    get_integer_from_atom(value.expand().as_view())
                        .ok_or_else(|| invalid(format!("line {id} is not an integer combination")))
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(CanonicalLine {
                from: node("nl_")?,
                to: node("nr_")?,
                external: coefficients[family.loops],
                loops: coefficients[..family.loops].to_vec(),
            })
        })
        .collect()
}

/// Exact rational numbers for the small linear systems of a routing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Fraction(i128, i128);

impl Fraction {
    fn new(numerator: i128, denominator: i128) -> Self {
        fn gcd(a: i128, b: i128) -> i128 {
            if b == 0 { a.abs() } else { gcd(b, a % b) }
        }
        let divisor = gcd(numerator, denominator).max(1) * denominator.signum();
        Fraction(numerator / divisor, denominator / divisor)
    }
    fn is_zero(self) -> bool {
        self.0 == 0
    }
    fn sub(self, other: Self) -> Self {
        Self::new(self.0 * other.1 - other.0 * self.1, self.1 * other.1)
    }
    fn mul(self, other: Self) -> Self {
        Self::new(self.0 * other.0, self.1 * other.1)
    }
    fn div(self, other: Self) -> Self {
        Self::new(self.0 * other.1, self.1 * other.0)
    }
    /// A FORM factor such as `+p1`, `-2*Q` or `+1/2*p3`.
    fn form_term(self, vector: &str) -> String {
        let sign = if self.0 < 0 { "-" } else { "+" };
        match (self.0.abs(), self.1) {
            (1, 1) => format!("{sign}{vector}"),
            (n, 1) => format!("{sign}{n}*{vector}"),
            (n, d) => format!("{sign}{n}/{d}*{vector}"),
        }
    }
}

/// The loop momenta in terms of the remaining (positive-power) lines and the
/// external momentum: row `j` holds the coefficients of `p1..pN` and of `Q`
/// in `k(j+1)`. `None` when the remaining lines do not span the loop
/// momenta: some loop carries no propagator and the integral is scaleless.
fn loop_momenta_in_lines(lines: &[&CanonicalLine], loops: usize) -> Option<Vec<Vec<Fraction>>> {
    // Gauss-Jordan elimination on [A | I] with A the loop coefficients of
    // the remaining lines (rows), keeping the row operations to read off the
    // left inverse X with X A = 1; then k = X (p - b Q).
    let rows = lines.len();
    let mut matrix = lines
        .iter()
        .enumerate()
        .map(|(row, line)| {
            let mut entries = line
                .loops
                .iter()
                .map(|&value| Fraction::new(value as i128, 1))
                .collect::<Vec<_>>();
            entries.extend((0..rows).map(|column| Fraction::new((column == row) as i128, 1)));
            entries
        })
        .collect::<Vec<_>>();
    let mut pivot_rows = Vec::new();
    for column in 0..loops {
        let pivot = (pivot_rows.len()..rows).find(|&row| !matrix[row][column].is_zero())?;
        let target = pivot_rows.len();
        matrix.swap(target, pivot);
        let pivot_value = matrix[target][column];
        for entry in matrix[target].iter_mut() {
            *entry = entry.div(pivot_value);
        }
        for row in 0..rows {
            let factor = matrix[row][column];
            if row != target && !factor.is_zero() {
                let pivot_row = matrix[target].clone();
                for (entry, pivot_entry) in matrix[row].iter_mut().zip(pivot_row) {
                    *entry = entry.sub(factor.mul(pivot_entry));
                }
            }
        }
        pivot_rows.push(target);
    }
    Some(
        (0..loops)
            .map(|j| {
                let inverse = &matrix[j][loops..];
                let mut row = inverse.to_vec();
                let shift = inverse
                    .iter()
                    .zip(lines)
                    .fold(Fraction::new(0, 1), |sum, (x, line)| {
                        sum.sub(x.mul(Fraction::new(-(line.external as i128), 1)))
                    });
                // k = X p - (X b) Q, with `shift` = X b.
                row.push(Fraction::new(0, 1).sub(shift));
                row
            })
            .collect(),
    )
}

/// The contracted graph of a sector: union-find representatives of the
/// vertices once all pinched lines are contracted.
struct ContractedGraph {
    representative: BTreeMap<i64, i64>,
}

impl ContractedGraph {
    fn new(lines: &[CanonicalLine], pinched: &[bool]) -> Self {
        let mut representative = BTreeMap::new();
        for line in lines {
            representative.insert(line.from, line.from);
            representative.insert(line.to, line.to);
        }
        let mut graph = ContractedGraph { representative };
        for (line, &pinched) in lines.iter().zip(pinched) {
            if pinched {
                let (a, b) = (graph.find(line.from), graph.find(line.to));
                graph.representative.insert(a.max(b), a.min(b));
            }
        }
        graph
    }

    fn find(&self, mut vertex: i64) -> i64 {
        while self.representative[&vertex] != vertex {
            vertex = self.representative[&vertex];
        }
        vertex
    }

    /// The scaleless rule of `docs/p_integrals.md` (section 1, item 6) on the
    /// remaining lines: the external vertices coincide, a self-loop remains,
    /// or removing one vertex isolates a part without external vertex.
    fn is_scaleless(&self, remaining: &[&CanonicalLine]) -> bool {
        let (q_in, q_out) = (self.find(1), self.find(2));
        if q_in == q_out {
            return true;
        }
        let edges = remaining
            .iter()
            .map(|line| (self.find(line.from), self.find(line.to)))
            .collect::<Vec<_>>();
        if edges.iter().any(|(a, b)| a == b) {
            return true;
        }
        let mut vertices = edges.iter().flat_map(|&(a, b)| [a, b]).collect::<Vec<_>>();
        vertices.sort();
        vertices.dedup();
        vertices.iter().any(|&removed| {
            // Components of the graph without `removed`: each must reach an
            // external vertex (unless it is `removed`'s own external vertex).
            let mut component = BTreeMap::new();
            for &start in vertices.iter().filter(|&&v| v != removed) {
                if component.contains_key(&start) {
                    continue;
                }
                let mut stack = vec![start];
                component.insert(start, start);
                while let Some(vertex) = stack.pop() {
                    for &(a, b) in &edges {
                        for (x, y) in [(a, b), (b, a)] {
                            if x == vertex && y != removed && !component.contains_key(&y) {
                                component.insert(y, start);
                                stack.push(y);
                            }
                        }
                    }
                }
            }
            let mut roots = component.values().copied().collect::<Vec<_>>();
            roots.sort();
            roots.dedup();
            roots.iter().any(|root| {
                component
                    .iter()
                    .filter(|(_, r)| *r == root)
                    .all(|(&vertex, _)| vertex != q_in && vertex != q_out)
            })
        })
    }
}

/// One homogeneous part of the numerator: its Forcer local, the scale power
/// `w = 2L - Σa + r`, the sign `(-1)^(Σa+r)` and the factor left outside FORM.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ForcerPart {
    pub name: String,
    pub scale_power: i64,
    pub sign: i64,
    pub coefficient: Atom,
}

/// A Forcer program for one canonical integral and the parts it computes.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ForcerInput {
    pub program: String,
    pub parts: Vec<ForcerPart>,
}

fn tensor_numerator(detail: impl std::fmt::Display) -> VakintError {
    VakintError::InvalidNumerator(format!(
        "p-integral: tensor p-integrals are not supported yet: {detail}"
    ))
}

/// A real rational coefficient as an exact fraction, if it is one.
fn rational(value: AtomView) -> Option<Fraction> {
    let AtomView::Num(number) = value else {
        return None;
    };
    match number.get_coeff_view() {
        CoefficientView::Natural(numerator, denominator, 0, _) => {
            Some(Fraction::new(numerator as i128, denominator as i128))
        }
        _ => None,
    }
}

/// Split a numerator coefficient into an exact rational and the rest.
fn split_rational(coefficient: &Atom) -> (Fraction, Atom) {
    if let Some(value) = rational(coefficient.as_view()) {
        return (value, Atom::num(1));
    }
    if let AtomView::Mul(product) = coefficient.as_view() {
        let mut value = Fraction::new(1, 1);
        let mut rest = Atom::num(1);
        for factor in product.iter() {
            match rational(factor) {
                Some(number) => value = value.mul(number),
                None => rest *= factor.to_owned(),
            }
        }
        return (value, rest);
    }
    (Fraction::new(1, 1), coefficient.clone())
}

/// The Forcer program of `topo(NAME(external, powers...))` with a scalar
/// numerator in the canonical `k(j)` and the user's head; `None` when the
/// sector is scaleless (an exact zero).
pub(crate) fn forcer_input(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
) -> Result<Option<ForcerInput>, VakintError> {
    let lines = canonical_lines(family)?;
    if powers.len() != lines.len() {
        return Err(VakintError::InvalidGenericExpression(format!(
            "{}: {} powers for {} lines",
            family.name,
            powers.len(),
            lines.len()
        )));
    }
    let (head, sigma) = external_head(external)?;
    let loops = family.loops;

    // A negative power is a numerator factor in the user's momenta.
    let mut numerator = canonical_numerator.to_owned();
    for (line, &power) in lines.iter().zip(powers) {
        if power < 0 {
            let momentum = line.momentum(external);
            numerator *= function!(S.dot, &momentum, &momentum).pow(Atom::num(-power));
        }
    }
    let pinched = powers.iter().map(|&power| power <= 0).collect::<Vec<_>>();
    let graph = ContractedGraph::new(&lines, &pinched);
    let remaining = lines
        .iter()
        .zip(powers)
        .filter(|(_, power)| **power > 0)
        .collect::<Vec<_>>();
    let remaining_lines = remaining.iter().map(|(line, _)| *line).collect::<Vec<_>>();
    if graph.is_scaleless(&remaining_lines) {
        return Ok(None);
    }
    let Some(loop_momenta) = loop_momenta_in_lines(&remaining_lines, loops) else {
        return Ok(None);
    };
    let total_power: i64 = remaining.iter().map(|(_, power)| **power).sum();

    // Propagators and vertices of the contracted graph: incoming momenta
    // count positively, Q enters at canonical vertex 1 and leaves at vertex 2.
    let mut vertices = BTreeMap::<i64, Vec<String>>::new();
    let mut propagators = String::from("1");
    for (index, (line, power)) in remaining.iter().enumerate() {
        let name = format!("p{}", index + 1);
        propagators.push_str(&format!("/({name}.{name})^{power}"));
        vertices
            .entry(graph.find(line.from))
            .or_default()
            .push(format!("-{name}"));
        vertices.entry(graph.find(line.to)).or_default().push(name);
    }
    vertices.entry(graph.find(1)).or_default().push("Q".into());
    vertices.entry(graph.find(2)).or_default().push("-Q".into());
    let topology = format!(
        "{propagators}*{}",
        vertices
            .values()
            .map(|momenta| format!("vx({})", momenta.join(",")))
            .collect::<Vec<_>>()
            .join("*")
    );
    let substitutions = loop_momenta
        .iter()
        .enumerate()
        .map(|(j, row)| {
            let image = row
                .iter()
                .enumerate()
                .filter(|(_, value)| !value.is_zero())
                .map(|(column, value)| {
                    let vector = if column + 1 == row.len() {
                        "Q".to_string()
                    } else {
                        format!("p{}", column + 1)
                    };
                    value.form_term(&vector)
                })
                .collect::<String>();
            format!("id k{} = {image};", j + 1)
        })
        .collect::<Vec<_>>()
        .join("\n");

    // Homogeneous parts: loop-momentum monomials, grouped by degree and by
    // the non-rational factor of their coefficient.
    let numerator = Vakint::convert_to_dot_notation(numerator.as_view()).expand();
    let mut loop_dots = Vec::new();
    let mut error = None;
    numerator.visitor(&mut |part| {
        if error.is_some() {
            return false;
        }
        if let AtomView::Fun(function) = part
            && function.get_symbol() == S.dot
            && part.contains_symbol(S.k)
        {
            match form_dot(part, &head, sigma, loops) {
                Ok(_) => loop_dots.push(part.to_owned()),
                Err(e) => error = Some(e),
            }
            return false;
        }
        true
    });
    if let Some(error) = error {
        return Err(error);
    }
    loop_dots.sort();
    loop_dots.dedup();
    let mut groups = Vec::<((i64, Atom), Vec<String>)>::new();
    for (monomial, coefficient) in numerator.coefficient_list::<i8>(&loop_dots) {
        if coefficient.contains_symbol(S.k) {
            return Err(tensor_numerator(format!(
                "the loop momenta must appear in scalar products only, found {coefficient}"
            )));
        }
        let (degree, monomial_sign, form) = form_monomial(monomial.as_view(), &head, sigma, loops)?;
        let (value, rest) = split_rational(&coefficient);
        let term = Fraction::new(monomial_sign as i128, 1)
            .mul(value)
            .form_term(&format!("({form})"));
        let key = (degree, rest);
        match groups.iter_mut().find(|(group, _)| *group == key) {
            Some((_, terms)) => terms.push(term),
            None => groups.push((key, vec![term])),
        }
    }

    let mut parts = Vec::new();
    let mut locals = Vec::new();
    for (index, ((degree, coefficient), terms)) in groups.into_iter().enumerate() {
        let name = format!("F{}", index + 1);
        locals.push(format!("Local {name} = ({})*{topology};", terms.join("")));
        let exponent = total_power + degree;
        parts.push(ForcerPart {
            name,
            scale_power: 2 * loops as i64 - total_power + degree,
            sign: if exponent % 2 == 0 { 1 } else { -1 },
            coefficient,
        });
    }
    let template = Template::parse_template(TEMPLATES.get("run_forcer.txt").unwrap()).unwrap();
    let variables = HashMap::from([
        (
            "vectors".to_string(),
            (1..=loops)
                .map(|j| format!("k{j}"))
                .collect::<Vec<_>>()
                .join(","),
        ),
        ("locals".to_string(), locals.join("\n")),
        ("substitutions".to_string(), substitutions),
        (
            "writes".to_string(),
            parts
                .iter()
                .map(|part| format!("#write <{}.out> \"%E\", {}", part.name, part.name))
                .collect::<Vec<_>>()
                .join("\n"),
        ),
    ]);
    let program = template
        .render(&RenderOptions {
            variables,
            ..Default::default()
        })
        .map_err(|error| VakintError::InvalidGenericExpression(error.to_string()))?;
    Ok(Some(ForcerInput { program, parts }))
}

/// A loop index `j` of `k(j)` with `1 <= j <= loops`.
fn loop_index(momentum: AtomView, loops: usize) -> Option<usize> {
    let AtomView::Fun(function) = momentum else {
        return None;
    };
    if function.get_symbol() != S.k || function.get_nargs() != 1 {
        return None;
    }
    let index = get_integer_from_atom(function.iter().next().unwrap())?;
    (1..=loops as i64)
        .contains(&index)
        .then_some(index as usize)
}

/// The FORM dot product of one loop-dependent `dot(...)` and its sign: the
/// user's head `p(i)` is `σ·Q`.
fn form_dot(
    dot: AtomView,
    head: &Atom,
    sigma: i64,
    loops: usize,
) -> Result<(i64, String), VakintError> {
    let AtomView::Fun(function) = dot else {
        return Err(VakintError::InvalidNumerator(format!(
            "p-integral: unexpected numerator factor {dot}"
        )));
    };
    let arguments = function.iter().collect::<Vec<_>>();
    let [a, b] = arguments.as_slice() else {
        return Err(tensor_numerator(dot));
    };
    let (loop_argument, other) = match (loop_index(*a, loops), loop_index(*b, loops)) {
        (Some(i), Some(j)) => return Ok((1, format!("k{i}.k{j}"))),
        (Some(i), None) => (i, *b),
        (None, Some(j)) => (j, *a),
        (None, None) => {
            return Err(VakintError::InvalidNumerator(format!(
                "p-integral: unknown loop momentum in {dot}"
            )));
        }
    };
    if other == head.as_view() {
        Ok((sigma, format!("k{loop_argument}.Q")))
    } else {
        Err(tensor_numerator(format!(
            "{dot} contracts a loop momentum with a vector other than the external momentum {head}"
        )))
    }
}

/// Degree, sign and FORM text of a monomial in loop-dependent dot products.
fn form_monomial(
    monomial: AtomView,
    head: &Atom,
    sigma: i64,
    loops: usize,
) -> Result<(i64, i64, String), VakintError> {
    let factors = match monomial {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        _ => vec![monomial],
    };
    let mut degree = 0;
    let mut sign = 1;
    let mut text = Vec::new();
    for factor in factors {
        let (base, exponent) = match factor {
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                (base, get_integer_from_atom(exponent))
            }
            AtomView::Num(_) if rational(factor) == Some(Fraction::new(1, 1)) => continue,
            _ => (factor, Some(1)),
        };
        let exponent = exponent.filter(|exponent| *exponent > 0).ok_or_else(|| {
            VakintError::InvalidNumerator(format!(
                "p-integral: the numerator must be polynomial in the loop momenta, found {factor}"
            ))
        })?;
        let (factor_sign, dot) = form_dot(base, head, sigma, loops)?;
        degree += exponent;
        if exponent % 2 == 1 {
            sign *= factor_sign;
        }
        text.push(format!("({dot})^{exponent}"));
    }
    if text.is_empty() {
        text.push("1".into());
    }
    Ok((degree, sign, text.join("*")))
}

/// The exact Forcer output of one local in vakint's p-master basis: `rat(a,b)`
/// becomes `a/b`, `ep` the configured `ε`, `GschemeConstants(x1,x2)` becomes
/// `PG(x1,x2)` and `Master(n)` becomes `PMaster(n)`. Anything else is an error.
pub(crate) fn parse_forcer_output(output: &str, epsilon: Symbol) -> Result<Atom, VakintError> {
    let failure = |detail: String| VakintError::FormOutputParsingError(output.into(), detail);
    let cleaned = output
        .replace("\\\n", "")
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>();
    if cleaned.is_empty() {
        return Err(failure("empty Forcer output".into()));
    }
    let parsed = vk_parse!(&cleaned).map_err(|error| failure(error.to_string()))?;
    let (x_, y_) = (Atom::var(S.x_), Atom::var(S.y_));
    let converted = parsed
        .replace(vk_symbol!("rat").call((&x_, &y_)).to_pattern())
        .with((&x_ / &y_).to_pattern())
        .replace(vk_symbol!("GschemeConstants").call((&x_, &y_)).to_pattern())
        .with(super::pg_symbol().call((&x_, &y_)).to_pattern())
        .replace(function!(vk_symbol!("Master"), &x_).to_pattern())
        .with(function!(super::p_master_symbol(), &x_).to_pattern())
        .replace(Atom::var(vk_symbol!("ep")).to_pattern())
        .with(Atom::var(epsilon).to_pattern());
    let mut unexpected = None;
    converted.visitor(&mut |part| {
        if unexpected.is_some() {
            return false;
        }
        match part {
            AtomView::Var(variable) if variable.get_symbol() != epsilon => {
                unexpected = Some(part.to_owned());
                false
            }
            AtomView::Fun(function) if function.get_symbol() == super::pg_symbol() => {
                let integers = function.iter().all(|x| get_integer_from_atom(x).is_some());
                if !(function.get_nargs() == 2 && integers) {
                    unexpected = Some(part.to_owned());
                }
                false
            }
            AtomView::Fun(function) if function.get_symbol() == super::p_master_symbol() => {
                if !(function.get_nargs() == 1
                    && matches!(function.iter().next(), Some(AtomView::Var(_))))
                {
                    unexpected = Some(part.to_owned());
                }
                false
            }
            AtomView::Fun(_) => {
                unexpected = Some(part.to_owned());
                false
            }
            _ => true,
        }
    });
    match unexpected {
        Some(part) => Err(failure(format!(
            "unexpected symbol {part} in the Forcer result"
        ))),
        None => Ok(converted),
    }
}

/// `form_exe_path` for one core, otherwise the sibling `tform` with
/// `n_cores` workers when it runs.
fn form_command(settings: &VakintSettings) -> Command {
    let form = Path::new(&settings.form_exe_path);
    if settings.n_cores > 1 {
        let tform = match form.parent() {
            Some(directory) if !directory.as_os_str().is_empty() => {
                directory.join("tform").to_string_lossy().into_owned()
            }
            _ => "tform".to_string(),
        };
        if crate::probe_command_output(&tform, &["-version"]).is_some() {
            let mut command = Command::new(tform);
            command.arg(format!("-w{}", settings.n_cores));
            return command;
        }
    }
    Command::new(form)
}

/// Run a Forcer program and return the `%E` text of each named local.
fn run_forcer(
    settings: &VakintSettings,
    program: &str,
    locals: &[&str],
) -> Result<Vec<String>, VakintError> {
    crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Form);
    crate::citations::USED_CITATIONS.record(crate::citations::CitationSource::Forcer);
    let sources = materialized_sources(settings.temporary_directory.as_deref())?;
    let run_directory =
        crate::create_unique_tmp_run_dir(settings.temporary_directory.as_deref(), "vakint_forcer")?;
    for name in ["forcer.h", "forcer"] {
        link(&sources.join(name), &run_directory.join(name))?;
    }
    fs::write(run_directory.join("run_forcer.frm"), program)?;
    let mut command = form_command(settings);
    command.arg("run_forcer.frm").current_dir(&run_directory);
    debug!("Running {} with command: {:?}", "Forcer".green(), command);
    let output = command
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .map_err(|error| VakintError::ForcerUnavailable {
            form_exe_path: settings.form_exe_path.clone(),
            reason: format!("could not be run ({error})"),
        })?;
    let failure = |stderr: String| {
        VakintError::FormError(
            stderr,
            String::from_utf8_lossy(&output.stdout).into(),
            format!("{command:?}"),
            run_directory.display().to_string(),
        )
    };
    if !output.status.success() {
        return Err(failure(String::from_utf8_lossy(&output.stderr).into()));
    }
    let results = locals
        .iter()
        .map(|local| {
            fs::read_to_string(run_directory.join(format!("{local}.out")))
                .map_err(|error| failure(format!("missing Forcer output for {local}: {error}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if settings.clean_tmp_dir {
        fs::remove_dir_all(&run_directory)?;
    }
    Ok(results)
}

#[cfg(unix)]
fn link(source: &Path, target: &Path) -> std::io::Result<()> {
    std::os::unix::fs::symlink(source, target)
}

#[cfg(not(unix))]
fn link(source: &Path, target: &Path) -> std::io::Result<()> {
    if source.is_dir() {
        fs::create_dir_all(target)?;
        for entry in fs::read_dir(source)? {
            let entry = entry?;
            link(&entry.path(), &target.join(entry.file_name()))?;
        }
        Ok(())
    } else {
        fs::copy(source, target).map(|_| ())
    }
}

/// The backend contract of `docs/p_integrals.md` (section 5): the exact
/// reduction of `topo(NAME(external, powers...))` with a scalar numerator in
/// the canonical `k(j)` and the user's head, one value per scale power.
pub(crate) fn reduce_p(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
    settings: &VakintSettings,
) -> Result<Vec<ReducedPropagatorValue>, VakintError> {
    let Some(input) = forcer_input(family, powers, canonical_numerator, external)? else {
        return Ok(Vec::new());
    };
    if input.parts.is_empty() {
        return Ok(Vec::new());
    }
    let names = input
        .parts
        .iter()
        .map(|part| part.name.as_str())
        .collect::<Vec<_>>();
    let outputs = run_forcer(settings, &input.program, &names)?;
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    combine_parts(&input.parts, &outputs, epsilon)
}

/// Attach the raw Forcer output of each part to its sign, scale power and
/// coefficient, summing the parts of equal scale power.
pub(crate) fn combine_parts(
    parts: &[ForcerPart],
    outputs: &[String],
    epsilon: Symbol,
) -> Result<Vec<ReducedPropagatorValue>, VakintError> {
    let mut values = BTreeMap::<i64, Atom>::new();
    for (part, output) in parts.iter().zip(outputs) {
        let reduced =
            parse_forcer_output(output, epsilon)? * Atom::num(part.sign) * &part.coefficient;
        *values.entry(part.scale_power).or_default() += reduced;
    }
    Ok(values
        .into_iter()
        .map(|(scale_power, reduced)| ReducedPropagatorValue {
            scale_power,
            reduced,
        })
        .collect())
}

impl Vakint {
    /// Reduce one matched p-integral with Forcer and finalize it. The
    /// family, powers and external momentum come from the matching witness.
    pub(crate) fn forcer_evaluate(
        &self,
        settings: &VakintSettings,
        numerator: AtomView,
        integral_specs: &ReplacementRules,
        options: &ForcerOptions,
    ) -> Result<Atom, VakintError> {
        let request = integral_specs
            .propagator_request(numerator)
            .ok_or_else(|| {
                VakintError::InvalidGenericExpression(format!(
                    "Integral {} is not a p-integral supported by Forcer.",
                    integral_specs.canonical_topology.get_integral().name
                ))
            })?;
        debug!(
            "Processing {} with {}: powers {:?}, external {}",
            request.family.name,
            "Forcer".green(),
            request.powers,
            request.external
        );
        let values = reduce_p(
            request.family,
            request.powers,
            request.canonical_numerator,
            &request.external,
            settings,
        )?;
        super::finalize::finalize(
            &values,
            request.family.loops,
            &request.external,
            settings,
            options.master_mode,
            "Forcer",
        )
    }
}

#[cfg(test)]
mod tests;
