//! The Forcer backend: massless propagator integrals reduced by the vendored
//! Forcer FORM program (B. Ruijl, T. Ueda, J.A.M. Vermaseren, arXiv:1704.06650).
//!
//! Forcer receives the canonical family graph, with propagators for the
//! positive-power lines only (never contracted: Forcer silently drops a graph
//! that is none of its topologies), and returns the exact reduction in its G-scheme master basis,
//! `rat(n,d) × GschemeConstants(x1,x2) × Master(name)`. Forcer sets `Q.Q = 1`,
//! so the numerator is split into homogeneous parts before the run; every
//! part becomes one [`ReducedPropagatorValue`]. Normalization, expansion and
//! master substitution are left to the shared p-finalizer (`finalize.rs`).

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, LazyLock, Mutex};
use std::time::Instant;

use colored::Colorize;
use log::{debug, warn};
use regex::Regex;
use string_template_plus::{Render, RenderOptions, Template};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::domains::integer::{IntegerRing, Z};
use symbolica::domains::rational::Q;
use symbolica::domains::rational_polynomial::RationalPolynomial;
use symbolica::function;
use symbolica::poly::{PolyVariable, polynomial::MultivariatePolynomial};
use version_compare::{Cmp, compare_to};

use super::ForcerOptions;
use super::parts::{Fraction, PDot, PMonomial, split};
use super::{
    FORCER_UPSTREAM_COMMIT, MINIMAL_FORCER_FORM_VERSION, PFamilyData, ReducedPropagatorValue,
};
use crate::stage_timing::StageSpan;
use crate::symbols::S;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    ReplacementRules, TEMPLATES, Vakint, VakintError, VakintSettings, get_integer_from_atom,
};

mod vendored {
    include!(concat!(env!("OUT_DIR"), "/vakint_forcer_sources.rs"));
}

mod batch;
mod cache;

pub(crate) use batch::reduce_inputs;

static FORCER_CACHE: LazyLock<cache::ForcerRunCache> = LazyLock::new(Default::default);

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

/// The shared copy of the vendored sources under `base` (`temporary_directory`,
/// else the system temporary directory), materialized once per process and
/// base; `None` when that copy cannot be trusted.
///
/// FORM includes `forcer.h` from this directory, so it is trusted only when
/// it is a private directory (mode 0700, not a link) of the current user
/// holding byte for byte the embedded sources. Such a directory is shared by
/// the processes of one user. Anything else at that path, e.g. a directory
/// planted by another user of a shared `/tmp` or a damaged copy, is never
/// used: every Forcer run then writes a private copy of the sources into its
/// own run directory ([`prepare_run_directory`]). The copy goes away with the
/// run directory after a successful run (`clean_tmp_dir`) and on its own after
/// a failed one, so no private copy outlives its run unless `clean_tmp_dir`
/// is off.
fn materialized_sources(base: Option<&str>) -> Result<Option<PathBuf>, VakintError> {
    static DIRECTORIES: LazyLock<Mutex<HashMap<PathBuf, Option<PathBuf>>>> =
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
                Some(shared)
            } else {
                write_sources(&staging, &sources)?;
                // Atomic publication; a concurrent process may have won the race.
                let published =
                    fs::symlink_metadata(&shared).is_err() && fs::rename(&staging, &shared).is_ok();
                if published {
                    Some(shared)
                } else {
                    fs::remove_dir_all(&staging)?;
                    if is_trusted(&shared, owner, &sources) {
                        Some(shared)
                    } else {
                        warn!(
                            "Ignoring the untrusted Forcer source directory '{}'; every Forcer \
                             run gets a private copy of the sources instead.",
                            shared.display()
                        );
                        None
                    }
                }
            }
        }
        None => {
            fs::remove_dir(&staging)?;
            None
        }
    };
    directories.insert(root, directory.clone());
    Ok(directory)
}

/// A fresh run directory for one Forcer program and whether it holds a private
/// copy of the sources: `forcer.h` and `forcer/` link to the shared sources, or
/// are that copy (the directory made private) when there is no trusted shared
/// copy. A directory that cannot be filled is removed.
fn prepare_run_directory(settings: &VakintSettings) -> Result<(PathBuf, bool), VakintError> {
    let base = settings.temporary_directory.as_deref();
    let shared = materialized_sources(base)?;
    let run_directory = crate::create_unique_tmp_run_dir(base, "vakint_forcer")?;
    let filled = match &shared {
        Some(sources) => ["forcer.h", "forcer"].iter().try_for_each(|name| {
            link(&sources.join(name), &run_directory.join(name)).map_err(VakintError::from)
        }),
        None => make_private(&run_directory)
            .map_err(VakintError::from)
            .and_then(|()| write_sources(&run_directory, &vendored_sources()?)),
    };
    if let Err(error) = filled {
        fs::remove_dir_all(&run_directory)?;
        return Err(error);
    }
    Ok((run_directory, shared.is_none()))
}

/// Remove the private copy of the sources from a failed run directory; the
/// program and FORM's files stay for diagnosis (`VakintError::FormError` names
/// the directory).
fn remove_private_sources(run_directory: &Path) {
    let removed = fs::remove_file(run_directory.join("forcer.h"))
        .and_then(|()| fs::remove_dir_all(run_directory.join("forcer")));
    if let Err(error) = removed {
        warn!(
            "Could not remove the private Forcer sources in '{}': {error}",
            run_directory.display()
        );
    }
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
    loops: usize,
    topology: String,
    polynomials: Vec<String>,
    substitutions: String,
}

/// `+p1`, `-2*Q` or `+1/2*p3`: a FORM factor of `vector`.
fn form_term(value: Fraction, vector: &str) -> String {
    let sign = if value.0 < 0 { "-" } else { "+" };
    match (value.0.abs(), value.1) {
        (1, 1) => format!("{sign}{vector}"),
        (n, 1) => format!("{sign}{n}*{vector}"),
        (n, d) => format!("{sign}{n}/{d}*{vector}"),
    }
}

/// The FORM dot product of a canonical scalar product.
fn form_dot(dot: PDot) -> String {
    match dot {
        PDot::Loops(i, j) => format!("k{i}.k{j}"),
        PDot::External(i) => format!("k{i}.Q"),
    }
}

/// The FORM text of a monomial, `1` when it is empty.
fn form_monomial(monomial: &PMonomial) -> String {
    if monomial.0.is_empty() {
        return "1".into();
    }
    monomial
        .0
        .iter()
        .map(|(dot, exponent)| format!("({})^{exponent}", form_dot(*dot)))
        .collect::<Vec<_>>()
        .join("*")
}

/// The Forcer program of `topo(NAME(external, powers...))` with a scalar
/// numerator in the canonical `k(j)` and the user's head; `None` when the
/// sector is scaleless (an exact zero). The numerator split is the shared
/// [`split`]; Forcer only renders it.
pub(crate) fn forcer_input(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
) -> Result<Option<ForcerInput>, VakintError> {
    forcer_input_in(family, powers, canonical_numerator, external, None)
}

/// [`forcer_input`] with every rational function of `epsilon` in the part
/// coefficients reduced inside FORM, as `rat(N, D)` with Forcer's own `ep`
/// (`d = 4 - 2 ep`, the same epsilon). See [`epsilon_channels`].
pub(crate) fn forcer_input_with_epsilon(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
    epsilon: Symbol,
) -> Result<Option<ForcerInput>, VakintError> {
    forcer_input_in(family, powers, canonical_numerator, external, Some(epsilon))
}

type EpsilonRational = RationalPolynomial<IntegerRing, u32>;

/// `coefficient = sum_j S_j R_j(epsilon)`: every `S_j` is free of `epsilon`
/// and every `R_j` is a nonzero rational function of `epsilon` alone. Real
/// rational numbers join `R_j`; other numbers stay in `S_j`. `None` when a
/// factor mixes `epsilon` with other symbols.
fn epsilon_coordinates(
    coefficient: &Atom,
    epsilon: Symbol,
) -> Option<Vec<(Atom, EpsilonRational)>> {
    let expanded = coefficient.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        term => vec![term],
    };
    let mut order = Vec::<Atom>::new();
    let mut groups = HashMap::<Atom, Atom>::new();
    for term in terms {
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            factor => vec![factor],
        };
        let mut spectator = Atom::num(1);
        let mut rational = Atom::num(1);
        for factor in factors {
            if factor.contains_symbol(epsilon) {
                if factor
                    .get_all_symbols(true)
                    .iter()
                    .any(|symbol| *symbol != epsilon)
                {
                    return None;
                }
                rational *= factor.to_owned();
            } else if matches!(factor, AtomView::Num(_)) && !factor.has_complex_coefficients() {
                rational *= factor.to_owned();
            } else {
                spectator *= factor.to_owned();
            }
        }
        match groups.get_mut(&spectator) {
            Some(sum) => *sum += rational,
            None => {
                order.push(spectator.clone());
                groups.insert(spectator, rational);
            }
        }
    }
    let variables = Arc::new(vec![PolyVariable::Symbol(epsilon)]);
    let mut coordinates = Vec::with_capacity(order.len());
    for spectator in order {
        let rational: EpsilonRational = groups[&spectator]
            .try_to_rational_polynomial(&Q, &Z, variables.clone())
            .ok()?;
        if rational.numerator.variables().as_ref() != variables.as_ref()
            || rational.denominator.variables().as_ref() != variables.as_ref()
        {
            return None;
        }
        if !rational.numerator.is_zero() {
            coordinates.push((spectator, rational));
        }
    }
    Some(coordinates)
}

/// One FORCER channel of an epsilon-aware request: what stays outside FORM,
/// and the polynomial reduced inside it.
struct EpsilonChannel {
    scale_power: i64,
    sign: i64,
    coefficient: Atom,
    polynomial: String,
}

/// Few channels for the parts of one request, without elimination fill-in.
/// Within each scale power (Forcer sets `Q.Q = 1`) and Wick sign, every part
/// coefficient is a vector over its epsilon-free spectators `S_j` with
/// rational functions of epsilon as entries. Parts with proportional vectors
/// share a channel: its outside coefficient is `sum_j S_j r_j / r_first` and
/// each part enters with `rat(r_first)`, so equal and proportional
/// coefficients merge. When there are more directions than spectators, each
/// spectator gets a channel instead, and each part enters it with its own
/// entry `rat(r_j)`; transverse and longitudinal projections of arbitrary
/// epsilon dependence then share two channels. Parts whose coefficients cannot
/// be separated keep a channel of their own coefficient.
fn epsilon_channels(parts: Vec<super::parts::PPart>, epsilon: Symbol) -> Vec<EpsilonChannel> {
    let polynomial_of = |part: &super::parts::PPart| {
        part.polynomial
            .iter()
            .map(|(value, monomial)| form_term(*value, &format!("({})", form_monomial(monomial))))
            .collect::<String>()
    };
    let mut channels = Vec::<EpsilonChannel>::new();
    let mut opaque = HashMap::<(i64, i64, Atom), usize>::new();
    let mut classes = Vec::<((i64, i64), Vec<(String, Vec<(Atom, EpsilonRational)>)>)>::new();
    for part in parts {
        let polynomial = polynomial_of(&part);
        let class = (part.scale_power, part.wick_sign);
        match epsilon_coordinates(&part.coefficient, epsilon) {
            Some(coordinates) if coordinates.is_empty() => {}
            Some(coordinates) => match classes.iter_mut().find(|(key, _)| *key == class) {
                Some((_, members)) => members.push((polynomial, coordinates)),
                None => classes.push((class, vec![(polynomial, coordinates)])),
            },
            None => {
                let key = (part.scale_power, part.wick_sign, part.coefficient);
                if let Some(&index) = opaque.get(&key) {
                    channels[index].polynomial.push_str(&polynomial);
                } else {
                    opaque.insert(key.clone(), channels.len());
                    channels.push(EpsilonChannel {
                        scale_power: key.0,
                        sign: key.1,
                        coefficient: key.2,
                        polynomial,
                    });
                }
            }
        }
    }
    for ((scale_power, sign), members) in classes {
        let mut spectators = Vec::<Atom>::new();
        let mut index_of = HashMap::<Atom, usize>::new();
        for (_, coordinates) in &members {
            for (spectator, _) in coordinates {
                if !index_of.contains_key(spectator) {
                    index_of.insert(spectator.clone(), spectators.len());
                    spectators.push(spectator.clone());
                }
            }
        }
        fn render(value: &EpsilonRational) -> String {
            format!(
                "rat({},{})",
                form_epsilon_polynomial(&value.numerator),
                form_epsilon_polynomial(&value.denominator)
            )
        }
        fn add(channel: &mut EpsilonChannel, value: &EpsilonRational, polynomial: &str) {
            if value.numerator.is_one() && value.denominator.is_one() {
                channel.polynomial.push_str(polynomial);
            } else {
                channel
                    .polynomial
                    .push_str(&format!("+{}*({polynomial})", render(value)));
            }
        }
        // Normalized directions; a sorted index list makes the key canonical.
        let directions = members
            .iter()
            .map(|(_, coordinates)| {
                let mut entries = coordinates
                    .iter()
                    .map(|(spectator, value)| (index_of[spectator], value))
                    .collect::<Vec<_>>();
                entries.sort_by_key(|(index, _)| *index);
                let first = entries[0].1;
                let normalized = entries
                    .iter()
                    .map(|(index, value)| (*index, *value / first))
                    .collect::<Vec<_>>();
                let key = normalized
                    .iter()
                    .map(|(index, value)| format!("{index}:{}", render(value)))
                    .collect::<Vec<_>>()
                    .join(";");
                (key, normalized, first.clone())
            })
            .collect::<Vec<_>>();
        let distinct = directions
            .iter()
            .map(|(key, _, _)| key.as_str())
            .collect::<std::collections::HashSet<_>>()
            .len();
        if distinct <= spectators.len() {
            let mut channel_of = HashMap::<String, usize>::new();
            for ((polynomial, _), (key, normalized, first)) in members.iter().zip(directions) {
                let channel = match channel_of.get(&key) {
                    Some(&channel) => channel,
                    None => {
                        let coefficient = normalized
                            .iter()
                            .map(|(index, value)| {
                                &spectators[*index]
                                    * (value.numerator.to_expression()
                                        / value.denominator.to_expression())
                            })
                            .fold(Atom::Zero, |sum, term| sum + term);
                        channel_of.insert(key, channels.len());
                        channels.push(EpsilonChannel {
                            scale_power,
                            sign,
                            coefficient,
                            polynomial: String::new(),
                        });
                        channels.len() - 1
                    }
                };
                add(&mut channels[channel], &first, polynomial);
            }
        } else {
            let offset = channels.len();
            for spectator in &spectators {
                channels.push(EpsilonChannel {
                    scale_power,
                    sign,
                    coefficient: spectator.clone(),
                    polynomial: String::new(),
                });
            }
            for (polynomial, coordinates) in &members {
                for (spectator, value) in coordinates {
                    add(
                        &mut channels[offset + index_of[spectator]],
                        value,
                        polynomial,
                    );
                }
            }
        }
    }
    channels
}

/// Make `rat` the active PolyRatFun before the input is defined when the
/// polynomials carry epsilon-rational coefficients. Otherwise `rat` is a plain
/// function until `Forcer()` enables it, so equal integrals with different
/// coefficients stay separate terms through the substitution and the first
/// sorts: a measured four-loop request then needed 268 s and 2 GB instead of
/// 28 s and 184 MB, with identical output.
pub(crate) fn with_input_rat_coefficients(locals: String, polynomials: &[String]) -> String {
    if polynomials
        .iter()
        .any(|polynomial| polynomial.contains("rat("))
    {
        format!("PolyRatFun rat;\n{locals}")
    } else {
        locals
    }
}

/// An integer polynomial in the single variable epsilon, as Forcer's `ep`.
fn form_epsilon_polynomial(polynomial: &MultivariatePolynomial<IntegerRing, u32>) -> String {
    if polynomial.is_zero() {
        return "0".into();
    }
    let mut rendered = String::new();
    for monomial in polynomial {
        let coefficient = monomial.coefficient.to_string();
        if !rendered.is_empty() && !coefficient.starts_with('-') {
            rendered.push('+');
        }
        rendered.push_str(&coefficient);
        match monomial.exponents.first().copied().unwrap_or(0) {
            0 => {}
            1 => rendered.push_str("*ep"),
            power => rendered.push_str(&format!("*ep^{power}")),
        }
    }
    rendered
}

fn forcer_input_in(
    family: &PFamilyData,
    powers: &[i64],
    canonical_numerator: AtomView,
    external: &Atom,
    epsilon: Option<Symbol>,
) -> Result<Option<ForcerInput>, VakintError> {
    let _stage = StageSpan::new("p_integral.forcer.input_prepare");
    let Some(split) = split(family, powers, canonical_numerator, external)? else {
        return Ok(None);
    };
    let loops = family.loops;
    let remaining = &split.remaining;

    // Forcer recognizes an input only when its vertex structure is one of its
    // own topologies, whose vertices carry at most four momenta, and silently
    // drops any other input (a wrong zero). Contracting pinched lines merges
    // vertices into five or more momenta, so the input is always the
    // uncontracted family graph (cubic, a top-level Forcer topology): a pinched
    // line keeps its place in the vertices and only loses its propagator, and
    // Forcer's own reduction hands the sector over to its subtopology.
    // Incoming momenta count positively, Q enters at canonical vertex 1 and
    // leaves at vertex 2.
    let mut vertices = BTreeMap::<i64, Vec<String>>::new();
    for (index, line) in split.lines.iter().enumerate() {
        let name = format!("p{}", index + 1);
        vertices
            .entry(line.from)
            .or_default()
            .push(format!("-{name}"));
        vertices.entry(line.to).or_default().push(name);
    }
    vertices.entry(1).or_default().push("Q".into());
    vertices.entry(2).or_default().push("-Q".into());
    let mut propagators = String::from("1");
    for (index, power) in remaining {
        let name = format!("p{}", index + 1);
        propagators.push_str(&format!("/({name}.{name})^{power}"));
    }
    let topology = format!(
        "{propagators}*{}",
        vertices
            .values()
            .map(|momenta| format!("vx({})", momenta.join(",")))
            .collect::<Vec<_>>()
            .join("*")
    );
    let substitutions = split
        .loop_momenta
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
                        format!("p{}", remaining[column].0 + 1)
                    };
                    form_term(*value, &vector)
                })
                .collect::<String>();
            format!("id k{} = {image};", j + 1)
        })
        .collect::<Vec<_>>()
        .join("\n");

    // Reduce all channels in one Forcer expression: Forcer otherwise repeats
    // its setup for each local. Without an epsilon symbol, every distinct
    // part coefficient keeps its own channel; with one, the coefficients'
    // rational functions of epsilon are reduced inside FORM.
    let channels = match epsilon {
        Some(epsilon) => epsilon_channels(split.parts, epsilon),
        None => split
            .parts
            .into_iter()
            .map(|part| EpsilonChannel {
                scale_power: part.scale_power,
                sign: part.wick_sign,
                polynomial: part
                    .polynomial
                    .iter()
                    .map(|(value, monomial)| {
                        form_term(*value, &format!("({})", form_monomial(monomial)))
                    })
                    .collect(),
                coefficient: part.coefficient,
            })
            .collect(),
    };
    let mut parts = Vec::with_capacity(channels.len());
    let mut polynomials = Vec::with_capacity(channels.len());
    for (index, channel) in channels.into_iter().enumerate() {
        polynomials.push(format!("({})", channel.polynomial));
        parts.push(ForcerPart {
            name: format!("F{}", index + 1),
            scale_power: channel.scale_power,
            sign: channel.sign,
            coefficient: channel.coefficient,
        });
    }
    let channels = (0..parts.len())
        .map(|index| format!("vakintchannel{index}"))
        .collect::<Vec<_>>();
    let locals = if parts.is_empty() {
        String::new()
    } else if parts.len() == 1 {
        format!("Local F1 = {}*{topology};", polynomials[0])
    } else {
        let tagged = polynomials
            .iter()
            .zip(&channels)
            .map(|(polynomial, channel)| format!("{channel}*({polynomial})"))
            .collect::<Vec<_>>()
            .join("+");
        format!(
            "S {};\nLocal VakintChannels = ({tagged})*{topology};",
            channels.join(",")
        )
    };
    let locals = with_input_rat_coefficients(locals, &polynomials);
    let mut writes = Vec::new();
    if parts.len() > 1 {
        // Scalar tags are bookkeeping only. Extract exact coefficients before
        // attaching each part's Wick sign, scale power or spectator factor.
        writes.push(format!(
            "Bracket {};\n.sort\nDrop VakintChannels;",
            channels.join(",")
        ));
        writes.extend(
            parts
                .iter()
                .zip(&channels)
                .map(|(part, channel)| format!("Local {} = VakintChannels[{channel}];", part.name)),
        );
        writes.push(".sort".into());
    }
    writes.extend(
        parts
            .iter()
            .map(|part| format!("#write <{}.out> \"%E\", {}", part.name, part.name)),
    );
    let template = Template::parse_template(TEMPLATES.get("run_forcer.txt").unwrap()).unwrap();
    let variables = HashMap::from([
        (
            "vectors".to_string(),
            (1..=loops)
                .map(|j| format!("k{j}"))
                .collect::<Vec<_>>()
                .join(","),
        ),
        ("locals".to_string(), locals),
        ("substitutions".to_string(), substitutions.clone()),
        ("writes".to_string(), writes.join("\n")),
    ]);
    let program = template
        .render(&RenderOptions {
            variables,
            ..Default::default()
        })
        .map_err(|error| VakintError::InvalidGenericExpression(error.to_string()))?;
    Ok(Some(ForcerInput {
        program,
        parts,
        loops,
        topology,
        polynomials,
        substitutions,
    }))
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
    // Keeping backend artifacts explicitly requests a fresh run directory.
    if !settings.forcer_cache || !settings.clean_tmp_dir {
        let _stage = StageSpan::new("p_integral.forcer.cache.bypass");
        return run_forcer_uncached(settings, program, locals);
    }
    // Child spans distinguish key construction, owned work, publication and
    // actual peer waits from other cache bookkeeping.
    let _stage = StageSpan::new("p_integral.forcer.cache_lookup_wait");
    let key = {
        let _stage = StageSpan::new("p_integral.forcer.cache_key");
        cache::ForcerRunKey::new(&form_command(settings), program, locals)
    };
    FORCER_CACHE.run(key, || run_forcer_uncached(settings, program, locals))
}

fn run_forcer_uncached(
    settings: &VakintSettings,
    program: &str,
    locals: &[&str],
) -> Result<Vec<String>, VakintError> {
    let _stage = StageSpan::new("p_integral.forcer.run_directory_io");
    let (run_directory, private_sources) = prepare_run_directory(settings)?;
    let results = run_forcer_in(settings, program, locals, &run_directory);
    match &results {
        Ok(_) if settings.clean_tmp_dir => fs::remove_dir_all(&run_directory)?,
        Err(_) if private_sources => remove_private_sources(&run_directory),
        _ => {}
    }
    results
}

/// [`run_forcer`] in its prepared directory.
fn run_forcer_in(
    settings: &VakintSettings,
    program: &str,
    locals: &[&str],
    run_directory: &Path,
) -> Result<Vec<String>, VakintError> {
    fs::write(run_directory.join("run_forcer.frm"), program)?;
    let mut command = form_command(settings);
    command.arg("run_forcer.frm").current_dir(run_directory);
    debug!("Running {} with command: {:?}", "Forcer".green(), command);
    let started = Instant::now();
    let output = {
        let _stage = StageSpan::new("p_integral.forcer.subprocess");
        command
            .stderr(Stdio::piped())
            .stdout(Stdio::piped())
            .output()
            .map_err(|error| VakintError::ForcerUnavailable {
                form_exe_path: settings.form_exe_path.clone(),
                reason: format!("could not be run ({error})"),
            })?
    };
    let subprocess_elapsed = started.elapsed();
    #[cfg(feature = "experimental-rustred")]
    benchmark::record(program, subprocess_elapsed);
    tracing::info!(
        uv = true,
        rstar = true,
        timing = true,
        locals = locals.len(),
        form_subprocess_seconds = subprocess_elapsed.as_secs_f64(),
        "Executed FORCER FORM subprocess"
    );
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
    if let Some(warning) = unrecognized_topology(&String::from_utf8_lossy(&output.stdout)) {
        return Err(VakintError::ForcerUnrecognizedTopology(warning));
    }
    locals
        .iter()
        .map(|local| {
            fs::read_to_string(run_directory.join(format!("{local}.out")))
                .map_err(|error| failure(format!("missing Forcer output for {local}: {error}")))
        })
        .collect()
}

/// Forcer's warning for an input graph that matches none of its topologies
/// (`matchTopologies` in `forcer/recognition.h`), whitespace-collapsed: such
/// an input is dropped from the result.
fn unrecognized_topology(stdout: &str) -> Option<String> {
    let found = stdout.find("did not match a topology")?;
    let start = stdout[..found].rfind('\n').map_or(0, |newline| newline + 1);
    let warning = stdout[start..].split("\n\n").next().unwrap_or_default();
    Some(warning.split_whitespace().collect::<Vec<_>>().join(" "))
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
    let _stage = StageSpan::new("p_integral.forcer.reduce_dispatch");
    let started = Instant::now();
    let preparation_started = Instant::now();
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let Some(input) =
        forcer_input_with_epsilon(family, powers, canonical_numerator, external, epsilon)?
    else {
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
    let preparation_elapsed = preparation_started.elapsed();
    let run_started = Instant::now();
    let outputs = run_forcer(settings, &input.program, &names)?;
    let run_elapsed = run_started.elapsed();
    let parsing_started = Instant::now();
    let result = combine_parts(&input.parts, &outputs, epsilon)?;
    let parsing_elapsed = parsing_started.elapsed();
    tracing::info!(
        uv = true,
        rstar = true,
        timing = true,
        loops = family.loops,
        input_parts = input.parts.len(),
        output_values = result.len(),
        input_preparation_seconds = preparation_elapsed.as_secs_f64(),
        run_or_cache_seconds = run_elapsed.as_secs_f64(),
        output_parsing_seconds = parsing_elapsed.as_secs_f64(),
        seconds = started.elapsed().as_secs_f64(),
        "Reduced p-integral with FORCER"
    );
    Ok(result)
}

/// Attach the raw Forcer output of each part to its sign, scale power and
/// coefficient, summing the parts of equal scale power.
pub(crate) fn combine_parts(
    parts: &[ForcerPart],
    outputs: &[String],
    epsilon: Symbol,
) -> Result<Vec<ReducedPropagatorValue>, VakintError> {
    let _stage = StageSpan::new("p_integral.forcer.parse_combine");
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
            spectator: None,
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
        let _stage = StageSpan::new("p_integral.forcer.dispatch");
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
        let _postprocess = StageSpan::new("p_integral.forcer.postprocess");
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

/// Process-global record of the FORM runs for
/// `examples/benchmark_p_integrals.rs` (not part of the API), off unless the
/// harness enables it. The Forcer materialisation base needs no hook: it is
/// the settings' `temporary_directory`.
#[cfg(feature = "experimental-rustred")]
#[doc(hidden)]
pub mod benchmark {
    use std::sync::Mutex;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::time::Duration;

    static RECORDING: AtomicBool = AtomicBool::new(false);
    static RUNS: Mutex<(Duration, Vec<String>)> = Mutex::new((Duration::ZERO, Vec::new()));

    /// Record every FORM run from now on.
    pub fn enable() {
        RECORDING.store(true, Ordering::SeqCst);
    }

    pub(super) fn record(program: &str, elapsed: Duration) {
        if !RECORDING.load(Ordering::SeqCst) {
            return;
        }
        let mut runs = RUNS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        runs.0 += elapsed;
        runs.1.push(program.to_owned());
    }

    /// The summed wall time of the FORM subprocesses and their programs, in
    /// run order, since the previous call.
    pub fn take() -> (Duration, Vec<String>) {
        let mut runs = RUNS
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        (std::mem::take(&mut runs.0), std::mem::take(&mut runs.1))
    }
}

pub mod fixtures;

#[cfg(test)]
mod tests;
