//! Opt-in dump of the massive-closure p-integral requests of R* remainders.
//!
//! Local export patch (not upstream). `PINT6_CLOSURE_DUMP=<directory>` writes
//! one record per massive closure exactly as `integrate_single_mass` hands it
//! to `MassiveClosure::evaluate_batch`: the closure (momentum, mass, power),
//! the Laurent depth and the input `numerator * topo(prop(...))`. With
//! `PINT6_SKIP_MIN_P_LOOPS=<n>` a closure whose remaining massless p-integral
//! has at least `n` loops is written and NOT evaluated (it contributes zero).
//! Without the variables nothing changes beyond one cached environment lookup.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::{
        Mutex, OnceLock,
        atomic::{AtomicUsize, Ordering},
    },
};

use eyre::{Result, eyre};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    graph::Graph as SGraph,
};
use vakint::{
    VakintSettings,
    p_integral::{BoundaryPattern, MassiveClosure, RemainderBoundary},
    symbols::S,
};

pub(super) struct DumpSettings {
    directory: PathBuf,
    pub(super) skip_min_p_loops: usize,
}

pub(super) fn settings() -> Option<&'static DumpSettings> {
    static SETTINGS: OnceLock<Option<DumpSettings>> = OnceLock::new();
    SETTINGS
        .get_or_init(|| {
            let directory = std::env::var_os("PINT6_CLOSURE_DUMP").filter(|v| !v.is_empty())?;
            let skip_min_p_loops = std::env::var("PINT6_SKIP_MIN_P_LOOPS")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(usize::MAX);
            std::fs::create_dir_all(&directory).ok()?;
            Some(DumpSettings {
                directory: directory.into(),
                skip_min_p_loops,
            })
        })
        .as_ref()
}

/// Whether `PINT6_CLOSURE_DUMP` is set (closure records are written).
pub fn enabled() -> bool {
    settings().is_some()
}

static CONTEXT: Mutex<String> = Mutex::new(String::new());
static TERM: Mutex<String> = Mutex::new(String::new());
static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
static SETTINGS_WRITTEN: OnceLock<()> = OnceLock::new();

/// Caller label attached to every following record (diagram, branch, term).
pub fn set_context(context: impl Into<String>) {
    *CONTEXT.lock().unwrap() = context.into();
}

/// R* term label, appended to the caller context.
pub fn set_term(term: impl Into<String>) {
    *TERM.lock().unwrap() = term.into();
}

/// The root projection of the R* outer integration that requests the
/// following closures, written into each record header (`channel`,
/// `taylor_order`, `normalisation`; new keys at the coordinator's request,
/// which byte comparisons of v1 records ignore). Set around one outer
/// integration by the evaluation that runs it ([`with_channel`]).
#[derive(Clone, Debug)]
pub(super) struct ClosureChannel {
    /// `projection metric` for the closed metric root (the root closure
    /// metric contracted in the root scope, the normalisation applied after
    /// the outer Taylor step); `none` without a root projector.
    pub(super) channel: &'static str,
    /// The outer Taylor order of the plan.
    pub(super) taylor_order: i64,
    /// The root normalisation, as text.
    pub(super) normalisation: &'static str,
    /// `GAMMALOOP_UV_CHILD_TRUNCATE` when on (`None` when off, so that the
    /// off-mode records are byte-identical): `{mode}@{d}` when the remainder
    /// of this outer integration was truncated to root grade `d`, `{mode}@-`
    /// when it was not (declined, or not a support root). Written as the
    /// `uv_child_truncate` key, so that records of different modes are
    /// never mixed unknowingly.
    pub(super) uv_child_truncate: Option<String>,
}

/// The header lines of `channel` (new keys, see [`ClosureChannel`]).
fn write_channel(out: &mut impl Write, channel: &ClosureChannel) -> std::io::Result<()> {
    writeln!(out, "channel {}", channel.channel)?;
    writeln!(out, "taylor_order {}", channel.taylor_order)?;
    writeln!(out, "normalisation {}", channel.normalisation)?;
    if let Some(truncate) = &channel.uv_child_truncate {
        writeln!(out, "uv_child_truncate {truncate}")?;
    }
    Ok(())
}

thread_local! {
    static CHANNEL: std::cell::RefCell<Option<ClosureChannel>> =
        const { std::cell::RefCell::new(None) };
}

/// Run `f` (one outer integration) with `channel` recorded in the headers of
/// the closures it writes.
pub(super) fn with_channel<T>(channel: ClosureChannel, f: impl FnOnce() -> T) -> T {
    // Restored on unwinding too (a panicking diagram).
    let scope = super::ScopedSlot::enter(&CHANNEL, channel);
    let result = f();
    scope.exit();
    result
}

/// The channel of the outer integration running on this thread.
#[cfg(test)]
pub(super) fn channel_active() -> bool {
    CHANNEL.with(|slot| slot.borrow().is_some())
}

fn loop_rank(view: AtomView) -> usize {
    match view {
        AtomView::Mul(product) => product.iter().map(loop_rank).sum(),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            let exponent = i64::try_from(exponent).ok().filter(|n| *n > 0).unwrap_or(0);
            exponent as usize * loop_rank(base)
        }
        AtomView::Fun(function) if function.get_symbol() == S.dot => function
            .iter()
            .filter(|argument| argument.contains_symbol(S.k))
            .count(),
        _ => 0,
    }
}

struct Line {
    id: String,
    a: i64,
    b: i64,
    momentum: Atom,
    power: i64,
}

fn lines(integrand: &Atom) -> Result<Vec<Line>> {
    let AtomView::Fun(topology) = integrand.as_view() else {
        return Ok(Vec::new());
    };
    if topology.get_symbol() != S.topo {
        return Ok(Vec::new());
    }
    let argument = topology.iter().next().unwrap();
    let factors = match argument {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        _ => vec![argument],
    };
    factors
        .into_iter()
        .map(|factor| {
            let AtomView::Fun(prop) = factor else {
                return Err(eyre!("invalid propagator {factor}"));
            };
            let arguments = prop.iter().collect::<Vec<_>>();
            let AtomView::Fun(ends) = arguments[1] else {
                return Err(eyre!("invalid propagator ends {factor}"));
            };
            let ends = ends
                .iter()
                .map(|end| i64::try_from(end).map_err(|_| eyre!("non-integer vertex {factor}")))
                .collect::<Result<Vec<_>>>()?;
            Ok(Line {
                id: arguments[0].to_canonical_string(),
                a: ends[0],
                b: ends[1],
                momentum: arguments[2].to_owned(),
                power: i64::try_from(arguments[4]).map_err(|_| eyre!("power {factor}"))?,
            })
        })
        .collect()
}

/// Canonical keys of the positive-power line graph after merging two-point
/// chains (dots): (shape without powers, shape with merged powers, dots).
fn topology_keys(lines: &[Line], external: &Atom) -> Result<(String, String, i64, usize)> {
    // External flow at each vertex from momentum conservation of all lines.
    let mut flow: BTreeMap<i64, Atom> = BTreeMap::new();
    for line in lines {
        *flow.entry(line.a).or_insert(Atom::Zero) += &line.momentum;
        *flow.entry(line.b).or_insert(Atom::Zero) -= &line.momentum;
    }
    let mut external_vertices = BTreeSet::new();
    for (vertex, sum) in &flow {
        if !sum.expand().coefficient(external.as_view()).is_zero() {
            external_vertices.insert(*vertex);
        }
    }
    // Positive lines as an undirected multigraph with powers.
    let mut edges: Vec<(i64, i64, i64)> = lines
        .iter()
        .filter(|line| line.power > 0)
        .map(|line| (line.a.min(line.b), line.a.max(line.b), line.power))
        .collect();
    loop {
        let mut degree: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
        for (index, (a, b, _)) in edges.iter().enumerate() {
            degree.entry(*a).or_default().push(index);
            degree.entry(*b).or_default().push(index);
        }
        let Some((vertex, pair)) = degree.iter().find(|(vertex, incident)| {
            incident.len() == 2 && incident[0] != incident[1] && !external_vertices.contains(vertex)
        }) else {
            break;
        };
        let (i, j) = (pair[0], pair[1]);
        let other = |(a, b, _): (i64, i64, i64)| if a == *vertex { b } else { a };
        let (x, y) = (other(edges[i]), other(edges[j]));
        let power = edges[i].2 + edges[j].2;
        let (hi, lo) = (i.max(j), i.min(j));
        edges.remove(hi);
        edges.remove(lo);
        edges.push((x.min(y), x.max(y), power));
    }
    let vertices = edges
        .iter()
        .flat_map(|(a, b, _)| [*a, *b])
        .chain(external_vertices.iter().copied())
        .collect::<BTreeSet<_>>();
    let index = vertices
        .iter()
        .enumerate()
        .map(|(i, v)| (*v, i))
        .collect::<BTreeMap<_, _>>();
    let canonical = |with_powers: bool| -> Result<String> {
        let mut graph = SGraph::new();
        for vertex in &vertices {
            graph.add_node(i64::from(external_vertices.contains(vertex)));
        }
        for (a, b, power) in &edges {
            graph
                .add_edge(
                    index[a],
                    index[b],
                    false,
                    if with_powers { *power } else { 1 },
                )
                .map_err(|error| eyre!(error))?;
        }
        Ok(graph
            .canonize()
            .graph
            .to_string()
            .trim()
            .replace('\n', "; "))
    };
    let dots = edges.iter().map(|(_, _, power)| power - 1).sum();
    Ok((canonical(false)?, canonical(true)?, dots, edges.len()))
}

/// Largest loop number of a 2-edge-connected block of the positive lines
/// (bridges removed), as in `GraphRegion::largest_cut_component`.
fn largest_block_loops(lines: &[Line]) -> usize {
    let edges = lines
        .iter()
        .filter(|line| line.power > 0)
        .map(|line| (line.a, line.b))
        .collect::<Vec<_>>();
    let mut adjacency: BTreeMap<i64, Vec<(i64, usize)>> = BTreeMap::new();
    for (index, (a, b)) in edges.iter().enumerate() {
        adjacency.entry(*a).or_default().push((*b, index));
        adjacency.entry(*b).or_default().push((*a, index));
    }
    // Iterative-free recursion is fine: at most a few dozen vertices.
    fn visit(
        vertex: i64,
        parent_edge: Option<usize>,
        adjacency: &BTreeMap<i64, Vec<(i64, usize)>>,
        time: &mut usize,
        discovery: &mut BTreeMap<i64, usize>,
        low: &mut BTreeMap<i64, usize>,
        bridges: &mut BTreeSet<usize>,
    ) {
        *time += 1;
        discovery.insert(vertex, *time);
        low.insert(vertex, *time);
        for &(next, edge) in &adjacency[&vertex] {
            if Some(edge) == parent_edge {
                continue;
            }
            if let Some(&seen) = discovery.get(&next) {
                let value = low[&vertex].min(seen);
                low.insert(vertex, value);
            } else {
                visit(next, Some(edge), adjacency, time, discovery, low, bridges);
                let value = low[&vertex].min(low[&next]);
                low.insert(vertex, value);
                if low[&next] > discovery[&vertex] {
                    bridges.insert(edge);
                }
            }
        }
    }
    let (mut time, mut discovery, mut low, mut bridges) =
        (0, BTreeMap::new(), BTreeMap::new(), BTreeSet::new());
    for vertex in adjacency.keys() {
        if !discovery.contains_key(vertex) {
            visit(
                *vertex,
                None,
                &adjacency,
                &mut time,
                &mut discovery,
                &mut low,
                &mut bridges,
            );
        }
    }
    // Union-find over non-bridge edges; cyclomatic number per component.
    let mut parent: BTreeMap<i64, i64> = adjacency.keys().map(|v| (*v, *v)).collect();
    fn find(parent: &mut BTreeMap<i64, i64>, v: i64) -> i64 {
        let p = parent[&v];
        if p == v {
            return v;
        }
        let root = find(parent, p);
        parent.insert(v, root);
        root
    }
    for (index, (a, b)) in edges.iter().enumerate() {
        if !bridges.contains(&index) {
            let (ra, rb) = (find(&mut parent, *a), find(&mut parent, *b));
            parent.insert(ra, rb);
        }
    }
    let mut count: BTreeMap<i64, (usize, BTreeSet<i64>)> = BTreeMap::new();
    for (index, (a, b)) in edges.iter().enumerate() {
        if !bridges.contains(&index) {
            let root = find(&mut parent, *a);
            let entry = count.entry(root).or_default();
            entry.0 += 1;
            entry.1.insert(*a);
            entry.1.insert(*b);
        }
    }
    count
        .values()
        .map(|(e, v)| e + 1 - v.len())
        .max()
        .unwrap_or(0)
}

/// Write one closure record. Returns the remaining p-integral loop count and
/// the record's sequence number.
pub(super) fn write(
    settings: &DumpSettings,
    vakint_settings: &VakintSettings,
    closure: &MassiveClosure,
    numerator: &Atom,
    integrand: &Atom,
    vacuum_loops: usize,
) -> Result<(usize, usize)> {
    write_record(
        settings,
        vakint_settings,
        closure,
        numerator,
        integrand,
        vacuum_loops,
        true,
    )
}

/// [`write`], with or without the `to_vakint_counts` line: the v1 records
/// written from a boundary (DESIGN 6c) have no to_vakint call of their own,
/// and the counts of an earlier call would be stale.
fn write_record(
    settings: &DumpSettings,
    vakint_settings: &VakintSettings,
    closure: &MassiveClosure,
    numerator: &Atom,
    integrand: &Atom,
    vacuum_loops: usize,
    counts: bool,
) -> Result<(usize, usize)> {
    let p_loops = vacuum_loops.saturating_sub(1);
    let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    if SETTINGS_WRITTEN.set(()).is_ok() {
        std::fs::write(
            settings
                .directory
                .join(format!("vakint-settings-{pid}.txt")),
            format!("{vakint_settings:#?}\n"),
        )?;
    }
    let context = format!("{} {}", CONTEXT.lock().unwrap(), TERM.lock().unwrap());
    let lines = lines(integrand)?;
    let (topology, dotted, dots, positive_lines) = if lines.is_empty() {
        ("tree".to_string(), "tree".to_string(), 0, 0)
    } else {
        topology_keys(&lines, &closure.momentum)?
    };
    let mut invariants = BTreeSet::new();
    numerator.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == S.dot
            && part.contains_symbol(S.k)
        {
            invariants.insert(part.to_owned());
            return false;
        }
        true
    });
    let invariants = invariants.into_iter().collect::<Vec<_>>();
    let monomials = numerator.coefficient_list::<i16>(&invariants);
    let max_rank = monomials
        .iter()
        .map(|(monomial, _)| loop_rank(monomial.as_view()))
        .max()
        .unwrap_or(0);
    let stem = format!("closure-{pid}-{sequence:06}");
    let mut out = BufWriter::new(File::create(
        settings.directory.join(format!("{stem}.txt")),
    )?);
    writeln!(out, "schema pint6.massive-closure-request.v1")?;
    writeln!(out, "sequence {sequence}")?;
    writeln!(out, "context {context}")?;
    writeln!(out, "vacuum_loops {vacuum_loops}")?;
    writeln!(out, "p_loops {p_loops}")?;
    writeln!(out, "largest_block_loops {}", largest_block_loops(&lines))?;
    writeln!(out, "evaluated {}", p_loops < settings.skip_min_p_loops)?;
    writeln!(
        out,
        "closure_momentum {}",
        closure.momentum.to_canonical_string()
    )?;
    writeln!(
        out,
        "closure_mass_squared {}",
        closure.mass_squared.to_canonical_string()
    )?;
    writeln!(out, "closure_power {}", closure.power)?;
    writeln!(
        out,
        "epsilon_terms {}",
        vakint_settings.number_of_terms_in_epsilon_expansion
    )?;
    writeln!(out, "lines {}", lines.len())?;
    for line in &lines {
        writeln!(
            out,
            "prop {}\t{}\t{}\t{}\t{}",
            line.id,
            line.a,
            line.b,
            line.momentum.to_canonical_string(),
            line.power
        )?;
    }
    writeln!(out, "positive_lines_after_chain_merge {positive_lines}")?;
    writeln!(out, "topology {topology}")?;
    writeln!(out, "dotted_topology {dotted}")?;
    writeln!(out, "dots {dots}")?;
    writeln!(out, "numerator_monomials {}", monomials.len())?;
    writeln!(out, "max_loop_rank {max_rank}")?;
    writeln!(
        out,
        "numerator_bytes {}",
        numerator.as_view().get_byte_size()
    )?;
    // Coordinator request: the root projection of this closure (new keys).
    if let Some(channel) = CHANNEL.with(|slot| slot.borrow().clone()) {
        write_channel(&mut out, &channel)?;
    }
    // Instrumentation: the counts of the to_vakint call that produced this
    // closure (see uv::approx::counts); a new key, before the input line.
    if counts && let Some(counts) = crate::uv::approx::counts::last_counts() {
        writeln!(out, "to_vakint_counts {counts}")?;
    }
    // Exactly the Atom MassiveClosure::evaluate_batch receives with this closure.
    writeln!(
        out,
        "input {}",
        (numerator * integrand).to_canonical_string()
    )?;
    out.flush()?;
    let mut tsv = BufWriter::new(File::create(
        settings.directory.join(format!("{stem}.monomials.tsv")),
    )?);
    for (monomial, _) in &monomials {
        writeln!(
            tsv,
            "{}\t{}",
            loop_rank(monomial.as_view()),
            monomial.to_canonical_string()
        )?;
    }
    tsv.flush()?;
    Ok((p_loops, sequence))
}

// ---------------------------------------------------------------------------
// The remainder-boundary record (DESIGN 6b, 6c; Addendum B: A1-A3)
// ---------------------------------------------------------------------------

/// Boundary record numbers: a counter of their own, so the v1 record numbers
/// are the same in every mode (DESIGN 6b).
static BOUNDARY_SEQUENCE: AtomicUsize = AtomicUsize::new(0);

/// Blocks with at most this many summands and bytes of text are inlined in
/// the header (`block_atom`).
const INLINE_SUMMANDS: usize = 64;
const INLINE_BYTES: usize = 16 << 10;
/// Lines per partial sum when a block file is read back.
const READ_CHUNK: usize = 1 << 16;

/// `PINT6_BOUNDARY_V1=1`: also write today's v1 records from the boundary
/// (DESIGN 6c; a gate tool for small graphs, never production).
pub(super) fn boundary_v1() -> bool {
    static V1: OnceLock<bool> = OnceLock::new();
    *V1.get_or_init(|| std::env::var("PINT6_BOUNDARY_V1").as_deref() == Ok("1"))
}

/// `PINT6_BOUNDARY_VERIFY=1`: read every boundary record back (streamed) right
/// after writing it and compare it with the boundary that was written.
pub(super) fn boundary_verify() -> bool {
    static VERIFY: OnceLock<bool> = OnceLock::new();
    *VERIFY.get_or_init(|| std::env::var("PINT6_BOUNDARY_VERIFY").as_deref() == Ok("1"))
}

/// The canonical text of an atom, asserted to fit on one line (A2): a stray
/// newline would silently split a summand.
fn one_line(atom: AtomView<'_>, what: &dyn std::fmt::Display) -> Result<String> {
    let text = atom.to_canonical_string();
    if text.contains(['\n', '\r']) {
        return Err(eyre!(
            "pint6 boundary record: {what} does not fit on one line"
        ));
    }
    Ok(text)
}

fn write_settings_once(settings: &DumpSettings, vakint_settings: &VakintSettings) -> Result<()> {
    if SETTINGS_WRITTEN.set(()).is_ok() {
        std::fs::write(
            settings
                .directory
                .join(format!("vakint-settings-{}.txt", std::process::id())),
            format!("{vakint_settings:#?}\n"),
        )?;
    }
    Ok(())
}

/// The colour table of `PINT6_COLOUR_INVARIANTS=<file>`: TAB-separated
/// `diagram branch symbolic` lines in colour_basis v1 (the `beta::` symbols),
/// and one comment line `# colour_basis <text>` that names the basis
/// (registry and sha256). Without it every symbolic coefficient is `pending`.
struct ColourTable {
    basis: String,
    symbolic: BTreeMap<(String, String), String>,
}

fn colour_table() -> Result<&'static ColourTable> {
    static TABLE: OnceLock<std::result::Result<ColourTable, String>> = OnceLock::new();
    TABLE
        .get_or_init(|| {
            let mut table = ColourTable {
                basis: "v1".to_owned(),
                symbolic: BTreeMap::new(),
            };
            let Some(path) = std::env::var_os("PINT6_COLOUR_INVARIANTS").filter(|v| !v.is_empty())
            else {
                return Ok(table);
            };
            let text = std::fs::read_to_string(&path).map_err(|error| {
                format!(
                    "PINT6_COLOUR_INVARIANTS {}: {error}",
                    path.to_string_lossy()
                )
            })?;
            for line in text.lines() {
                if let Some(basis) = line.strip_prefix("# colour_basis ") {
                    table.basis = basis.trim().to_owned();
                    continue;
                }
                if line.starts_with('#') || line.trim().is_empty() {
                    continue;
                }
                let fields = line.split('\t').collect::<Vec<_>>();
                if fields.len() != 3 {
                    return Err(format!("PINT6_COLOUR_INVARIANTS: bad line {line:?}"));
                }
                table.symbolic.insert(
                    (fields[0].trim().to_owned(), fields[1].trim().to_owned()),
                    fields[2].trim().to_owned(),
                );
            }
            Ok(table)
        })
        .as_ref()
        .map_err(|error| eyre!("{error}"))
}

/// `key=value` of the caller context (`diagram=.. colour_branch=..
/// colour_coefficient=..`).
fn context_field<'c>(context: &'c str, key: &str) -> Option<&'c str> {
    context
        .split_whitespace()
        .find_map(|field| field.strip_prefix(key)?.strip_prefix('='))
}

/// How one block was written.
struct BlockRecord {
    summands: usize,
    text_bytes: usize,
    /// sha256 (hex) of the block file as written; `None` for an inline block.
    sha256: Option<String>,
    /// The file name, or `None` for an inline block (its text).
    file: std::result::Result<String, String>,
}

/// Write one block: inline (small), or streamed through `zstd -q -3 -T1`
/// into its own file, one summand per line; the number of lines written is
/// checked against the `summands` line (A3).
fn write_block(
    directory: &std::path::Path,
    stem: &str,
    index: usize,
    block: &Atom,
    record: &str,
) -> Result<BlockRecord> {
    let summands = match block.as_view() {
        _ if block.is_zero() => Vec::new(),
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    if summands.len() <= INLINE_SUMMANDS && block.as_view().get_byte_size() <= INLINE_BYTES {
        let text = one_line(block.as_view(), &format_args!("block {index}"))?;
        if text.len() <= INLINE_BYTES {
            return Ok(BlockRecord {
                summands: summands.len(),
                text_bytes: text.len(),
                sha256: None,
                file: Err(text),
            });
        }
    }
    let name = format!("{stem}.block-{index:03}.txt.zst");
    let path = directory.join(&name);
    let mut child = std::process::Command::new("zstd")
        .args(["-q", "-3", "-T1", "-o"])
        .arg(&path)
        .stdin(std::process::Stdio::piped())
        .spawn()
        .map_err(|error| eyre!("pint6 boundary record: cannot start zstd: {error}"))?;
    let mut text_bytes = 0;
    let mut lines = 0;
    {
        let stdin = child.stdin.take().expect("piped stdin");
        let mut out = BufWriter::with_capacity(1 << 20, stdin);
        let head = format!(
            "schema pint6.boundary-block.v1\nblock {index} record={record}\nsummands {}\n",
            summands.len()
        );
        out.write_all(head.as_bytes())?;
        text_bytes += head.len();
        for (position, summand) in summands.iter().enumerate() {
            let line = one_line(
                *summand,
                &format_args!("summand {position} of block {index}"),
            )?;
            out.write_all(line.as_bytes())?;
            out.write_all(b"\n")?;
            text_bytes += line.len() + 1;
            lines += 1;
        }
        out.flush()?;
    }
    let status = child.wait()?;
    if !status.success() {
        return Err(eyre!(
            "pint6 boundary record: zstd failed on {name}: {status}"
        ));
    }
    if lines != summands.len() {
        return Err(eyre!(
            "pint6 boundary record: block {index} has {} summands, {lines} lines written",
            summands.len()
        ));
    }
    Ok(BlockRecord {
        summands: summands.len(),
        text_bytes,
        sha256: Some(file_sha256(&path)?),
        file: Ok(name),
    })
}

/// What [`write_boundary`] wrote.
#[derive(Debug)]
pub(super) struct BoundaryWritten {
    pub(super) sequence: usize,
    pub(super) header: PathBuf,
    pub(super) block_files: usize,
    pub(super) text_bytes: usize,
}

/// Write the boundary record of one built remainder (DESIGN 6b): the block
/// files first, then the header `boundary-<pid>-<bseq>.txt`. `references`:
/// pattern order (kept and zero patterns, key order) -> sequence of the v1
/// record of the same closure (check mode, or `PINT6_BOUNDARY_V1=1`).
pub(super) fn write_boundary(
    settings: &DumpSettings,
    vakint_settings: &VakintSettings,
    built: &super::boundary::BuiltBoundary,
    counts: &serde_json::Value,
    evaluated: bool,
    references: &BTreeMap<usize, usize>,
) -> Result<BoundaryWritten> {
    write_settings_once(settings, vakint_settings)?;
    let sequence = BOUNDARY_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    let stem = format!("boundary-{pid}-{sequence:06}");
    let boundary = &built.boundary;
    let loops = built.loops;
    // Blocks first: their sizes go into the header.
    let mut blocks = Vec::with_capacity(boundary.blocks.len());
    for (index, block) in boundary.blocks.iter().enumerate() {
        blocks.push(write_block(
            &settings.directory,
            &stem,
            index,
            block,
            &format!("{}-{sequence}", built.integration_id),
        )?);
    }
    let context = CONTEXT.lock().unwrap().clone();
    let term = TERM.lock().unwrap().clone();
    let colour = colour_table()?;
    let reference = |order: usize| {
        references
            .get(&order)
            .map_or_else(|| "none".to_owned(), |sequence| sequence.to_string())
    };
    let mut classes = BTreeMap::<i64, Vec<usize>>::new();
    for (index, pattern) in boundary.patterns.iter().enumerate() {
        classes
            .entry(pattern.closure_power)
            .or_default()
            .push(index);
    }
    let class_of = classes
        .iter()
        .enumerate()
        .flat_map(|(class, (_, patterns))| patterns.iter().map(move |pattern| (*pattern, class)))
        .collect::<BTreeMap<_, _>>();
    let join = |items: &mut dyn Iterator<Item = String>| {
        let items = items.collect::<Vec<_>>();
        if items.is_empty() {
            "none".to_owned()
        } else {
            items.join(",")
        }
    };
    let path = settings.directory.join(format!("{stem}.txt"));
    let mut out = BufWriter::new(File::create(&path)?);
    writeln!(out, "schema pint6.remainder-boundary.v1")?;
    writeln!(out, "sequence {sequence}")?;
    writeln!(out, "context {context} {term}")?;
    writeln!(out, "label {}", boundary.label)?;
    writeln!(out, "integration_id {}", built.integration_id)?;
    writeln!(out, "vacuum_loops {loops}")?;
    writeln!(out, "p_loops {}", loops.saturating_sub(1))?;
    writeln!(out, "evaluated {evaluated}")?;
    writeln!(
        out,
        "epsilon_terms {}",
        vakint_settings.number_of_terms_in_epsilon_expansion
    )?;
    writeln!(out, "taylor_degree {}", built.degree)?;
    // Colour, per branch (colour_basis v1).
    let diagram = context_field(&context, "diagram").unwrap_or("none");
    let branch = context_field(&context, "colour_branch").unwrap_or("none");
    writeln!(out, "colour_basis {}", colour.basis)?;
    writeln!(out, "colour_branch {branch}")?;
    writeln!(
        out,
        "colour_coefficient {}",
        context_field(&context, "colour_coefficient").unwrap_or("none")
    )?;
    writeln!(
        out,
        "colour_coefficient_symbolic {}",
        colour
            .symbolic
            .get(&(diagram.to_owned(), branch.to_owned()))
            .map_or("pending", String::as_str)
    )?;
    let momentum = one_line(boundary.momentum.as_view(), &"closure_momentum")?;
    writeln!(out, "closure_momentum {momentum}")?;
    writeln!(
        out,
        "closure_mass_squared {}",
        one_line(boundary.mass_squared.as_view(), &"closure_mass_squared")?
    )?;
    writeln!(
        out,
        "normalisation {}",
        one_line(built.normalisation.as_view(), &"normalisation")?
    )?;
    for (j, edge) in built.basis.iter().enumerate() {
        writeln!(out, "vacuum {}\t{edge}", j + 1)?;
    }
    // The vacuum labels, then the closing momentum to itself (A1).
    writeln!(out, "routing {}", boundary.routing.len() + 1)?;
    for (label, image) in &boundary.routing {
        writeln!(
            out,
            "route {}\t{}",
            one_line(label.as_view(), &"route label")?,
            one_line(image.as_view(), &"route image")?
        )?;
    }
    writeln!(out, "route {momentum}\t{momentum}")?;
    writeln!(out, "counts {counts}")?;
    writeln!(out, "blocks {}", blocks.len())?;
    for (index, block) in blocks.iter().enumerate() {
        let grade = &built.grades[index];
        writeln!(
            out,
            "block {index}\tdegree={}\tstructure={}\tmass_power={}\tsummands={}\tatom_bytes={}\ttext_bytes={}\tfile={}\tsha256={}",
            grade.degree(),
            one_line(
                grade.structure_atom(&built.external).as_view(),
                &"block structure"
            )?,
            grade.mass_power,
            block.summands,
            boundary.blocks[index].as_view().get_byte_size(),
            block.text_bytes,
            block.file.as_ref().map_or("inline", String::as_str),
            block.sha256.as_deref().unwrap_or("inline")
        )?;
    }
    for (index, block) in blocks.iter().enumerate() {
        if let Err(text) = &block.file {
            writeln!(out, "block_atom {index}\t{text}")?;
        }
    }
    writeln!(out, "classes {}", classes.len())?;
    for (class, (power, patterns)) in classes.iter().enumerate() {
        writeln!(
            out,
            "class {class}\tclosure_power={power}\tpatterns={}\tref_closures={}",
            join(&mut patterns.iter().map(usize::to_string)),
            join(
                &mut patterns
                    .iter()
                    .filter_map(|index| references.get(&built.patterns[*index].order))
                    .map(usize::to_string)
            )
        )?;
    }
    writeln!(out, "patterns {}", boundary.patterns.len())?;
    writeln!(
        out,
        "pattern_terms {}",
        boundary
            .patterns
            .iter()
            .map(|pattern| pattern.terms.len())
            .sum::<usize>()
    )?;
    for (index, (pattern, meta)) in boundary.patterns.iter().zip(&built.patterns).enumerate() {
        let lines = lines(&pattern.topology)?;
        writeln!(
            out,
            "pattern {index}\tclass={}\tclosure_power={}\tterms={}\tlines={}\tunprojected_nonzero={}\tref_closure={}\tkey={}",
            class_of[&index],
            pattern.closure_power,
            pattern.terms.len(),
            lines.len(),
            meta.unprojected_nonzero,
            reference(meta.order),
            one_line(meta.key.as_view(), &"pattern key")?
        )?;
        writeln!(
            out,
            "pattern_integral {index}\t{}",
            one_line(meta.integral.as_view(), &"pattern integral")?
        )?;
        writeln!(
            out,
            "pattern_topology {index}\t{}",
            one_line(pattern.topology.as_view(), &"pattern topology")?
        )?;
        for line in &lines {
            writeln!(
                out,
                "pattern_prop {index}\t{}\t{}\t{}\t{}\t{}",
                line.id,
                line.a,
                line.b,
                one_line(line.momentum.as_view(), &"pattern prop momentum")?,
                line.power
            )?;
        }
        for (block, coefficient) in &pattern.terms {
            writeln!(
                out,
                "pattern_term {index}\t{block}\t{}",
                one_line(coefficient.as_view(), &"pattern term")?
            )?;
        }
    }
    // Informational: the zero patterns (no pattern of the boundary).
    writeln!(out, "zero_patterns {}", built.zero_patterns.len())?;
    for meta in &built.zero_patterns {
        writeln!(
            out,
            "zero_pattern {}\tclosure_power={}\tunprojected_nonzero={}\tref_closure={}\tkey={}",
            meta.order,
            meta.closure_power,
            meta.unprojected_nonzero,
            reference(meta.order),
            one_line(meta.key.as_view(), &"zero pattern key")?
        )?;
    }
    writeln!(out, "end")?;
    out.flush()?;
    Ok(BoundaryWritten {
        sequence,
        header: path,
        block_files: blocks.iter().filter(|block| block.file.is_ok()).count(),
        text_bytes: blocks.iter().map(|block| block.text_bytes).sum(),
    })
}

/// Today's v1 records from a boundary (DESIGN 6c): one per pattern in key
/// order, kept and zero patterns (numerator 0) alike, each numerator routed
/// and expanded as today's dump does. Returns pattern order -> sequence.
pub(super) fn write_v1_from_boundary(
    settings: &DumpSettings,
    vakint_settings: &VakintSettings,
    built: &super::boundary::BuiltBoundary,
) -> Result<BTreeMap<usize, usize>> {
    let boundary = &built.boundary;
    let mut records = built
        .patterns
        .iter()
        .enumerate()
        .map(|(index, meta)| (meta.order, Some(index), meta))
        .chain(
            built
                .zero_patterns
                .iter()
                .map(|meta| (meta.order, None, meta)),
        )
        .collect::<Vec<_>>();
    records.sort_by_key(|(order, _, _)| *order);
    let mut references = BTreeMap::new();
    for (order, index, meta) in records {
        if index.is_none() && !meta.unprojected_nonzero {
            continue;
        }
        let closure = MassiveClosure {
            momentum: boundary.momentum.clone(),
            mass_squared: boundary.mass_squared.clone(),
            power: meta.closure_power,
        };
        let numerator = match index {
            Some(index) => {
                let input = boundary.direct_input(&boundary.patterns[index]);
                let routing = input
                    .routing
                    .iter()
                    .map(|(head, image)| {
                        symbolica::id::Replacement::new(head.to_pattern(), image.clone())
                    })
                    .collect::<Vec<_>>();
                input.numerator.replace_multiple(&routing).expand()
            }
            None => Atom::Zero,
        };
        let (_, sequence) = write_record(
            settings,
            vakint_settings,
            &closure,
            &numerator,
            &meta.topology,
            built.loops,
            false,
        )?;
        references.insert(order, sequence);
    }
    Ok(references)
}

/// Parse one canonical atom of a record.
fn parse_atom(text: &str) -> Result<Atom> {
    vakint::vakint_parse!(text).map_err(|error| eyre!("cannot parse {:.80}: {error}", text))
}

/// Stream one block file back: `zstd -dcq`, one summand per line, summed in
/// chunks of [`READ_CHUNK`] lines; the line count must equal `summands`.
fn read_block_file(path: &std::path::Path, index: usize) -> Result<Atom> {
    use std::io::BufRead;
    let mut child = std::process::Command::new("zstd")
        .arg("-dcq")
        .arg(path)
        .stdout(std::process::Stdio::piped())
        .spawn()?;
    let reader =
        std::io::BufReader::with_capacity(1 << 20, child.stdout.take().expect("piped stdout"));
    let mut lines = reader.lines();
    let mut next = || -> Result<String> {
        lines
            .next()
            .ok_or_else(|| eyre!("{}: short block file", path.display()))?
            .map_err(Into::into)
    };
    if next()? != "schema pint6.boundary-block.v1" {
        return Err(eyre!(
            "{}: not a pint6.boundary-block.v1 file",
            path.display()
        ));
    }
    let block_line = next()?;
    if block_line != format!("block {index}")
        && !block_line.starts_with(&format!("block {index} record="))
    {
        return Err(eyre!("{}: wrong block number", path.display()));
    }
    let summands = next()?
        .strip_prefix("summands ")
        .and_then(|count| count.parse::<usize>().ok())
        .ok_or_else(|| eyre!("{}: no summands line", path.display()))?;
    let mut partial = Vec::new();
    let mut chunk = Vec::with_capacity(READ_CHUNK.min(summands));
    let mut read = 0;
    for line in lines.by_ref() {
        chunk.push(parse_atom(&line?)?);
        read += 1;
        if chunk.len() == READ_CHUNK {
            partial.push(Atom::add_many(chunk.iter()));
            chunk.clear();
        }
    }
    partial.push(Atom::add_many(chunk.iter()));
    let status = child.wait()?;
    if !status.success() {
        return Err(eyre!("{}: zstd -d failed: {status}", path.display()));
    }
    if read != summands {
        return Err(eyre!(
            "{}: {read} lines for {summands} summands",
            path.display()
        ));
    }
    Ok(Atom::add_many(partial.iter()))
}

/// Read a boundary record back, streaming its block files (DESIGN 6b): the
/// boundary (routing as written, the closing route included) and the
/// header's single-valued keys.
pub(super) fn read_boundary(
    header: &std::path::Path,
) -> Result<(RemainderBoundary, BTreeMap<String, String>)> {
    let directory = header.parent().unwrap_or(std::path::Path::new("."));
    let text = std::fs::read_to_string(header)?;
    let mut keys = BTreeMap::new();
    let mut routing = Vec::new();
    let mut blocks = BTreeMap::<usize, Atom>::new();
    let mut block_count = None;
    let mut patterns = BTreeMap::<usize, BoundaryPattern>::new();
    let mut ended = false;
    for line in text.lines() {
        let (key, rest) = line.split_once(' ').unwrap_or((line, ""));
        let fields = rest.split('\t').collect::<Vec<_>>();
        let field = |name: &str| {
            fields
                .iter()
                .find_map(|field| field.strip_prefix(name)?.strip_prefix('='))
                .ok_or_else(|| eyre!("{}: no {name} in {line:.120}", header.display()))
        };
        let index = || -> Result<usize> {
            fields[0]
                .parse::<usize>()
                .map_err(|_| eyre!("{}: bad index in {line:.120}", header.display()))
        };
        match key {
            "route" => routing.push((parse_atom(fields[0])?, parse_atom(fields[1])?)),
            "blocks" => block_count = Some(rest.parse::<usize>()?),
            "block" => {
                let file = field("file")?;
                if file != "inline" {
                    blocks.insert(index()?, read_block_file(&directory.join(file), index()?)?);
                }
            }
            "block_atom" => {
                blocks.insert(index()?, parse_atom(fields[1])?);
            }
            "pattern" => {
                patterns.insert(
                    index()?,
                    BoundaryPattern {
                        topology: Atom::Zero,
                        closure_power: field("closure_power")?.parse()?,
                        terms: Vec::new(),
                    },
                );
            }
            "pattern_topology" => {
                let pattern = patterns
                    .get_mut(&index()?)
                    .ok_or_else(|| eyre!("topology before pattern"))?;
                pattern.topology = parse_atom(fields[1])?;
            }
            "pattern_term" => {
                let pattern = patterns
                    .get_mut(&index()?)
                    .ok_or_else(|| eyre!("term before pattern"))?;
                pattern
                    .terms
                    .push((fields[1].parse()?, parse_atom(fields[2])?));
            }
            "end" => ended = true,
            "vacuum" | "class" | "pattern_integral" | "pattern_prop" | "zero_pattern" => {}
            _ => {
                keys.insert(key.to_owned(), rest.to_owned());
            }
        }
    }
    if !ended {
        return Err(eyre!("{}: no end line", header.display()));
    }
    let block_count = block_count.ok_or_else(|| eyre!("{}: no blocks line", header.display()))?;
    if blocks.len() != block_count || blocks.keys().copied().ne(0..block_count) {
        return Err(eyre!(
            "{}: blocks 0..{block_count} expected",
            header.display()
        ));
    }
    let route_count = keys
        .get("routing")
        .map(|count| count.parse::<usize>())
        .transpose()?;
    if route_count != Some(routing.len()) {
        return Err(eyre!(
            "{}: routing count differs from the route lines",
            header.display()
        ));
    }
    let boundary = RemainderBoundary::new(
        keys.get("label").cloned().unwrap_or_default(),
        parse_atom(
            keys.get("closure_momentum")
                .ok_or_else(|| eyre!("no closure_momentum"))?,
        )?,
        parse_atom(
            keys.get("closure_mass_squared")
                .ok_or_else(|| eyre!("no closure_mass_squared"))?,
        )?,
        routing,
        blocks.into_values().collect(),
        patterns.into_values().collect(),
    )
    .map_err(|error| eyre!("{error}"))?;
    Ok((boundary, keys))
}

/// Compare a record read back with the boundary that was written: routing
/// (plus the closing route), blocks, patterns; atoms by equality.
pub(super) fn verify_boundary(written: &RemainderBoundary, header: &std::path::Path) -> Result<()> {
    let (read, _) = read_boundary(header)?;
    let mut routing = written.routing.clone();
    routing.push((written.momentum.clone(), written.momentum.clone()));
    let differs = |what: String| {
        Err(eyre!(
            "PINT6_BOUNDARY_VERIFY {}: {what} differs",
            header.display()
        ))
    };
    if read.routing != routing {
        return differs("the routing".to_owned());
    }
    if read.momentum != written.momentum || read.mass_squared != written.mass_squared {
        return differs("the closing line".to_owned());
    }
    if read.blocks.len() != written.blocks.len() {
        return differs("the block count".to_owned());
    }
    for (index, (a, b)) in read.blocks.iter().zip(&written.blocks).enumerate() {
        if a != b {
            return differs(format!("block {index}"));
        }
    }
    if read.patterns.len() != written.patterns.len() {
        return differs("the pattern count".to_owned());
    }
    for (index, (a, b)) in read.patterns.iter().zip(&written.patterns).enumerate() {
        if a.topology != b.topology || a.closure_power != b.closure_power || a.terms != b.terms {
            return differs(format!("pattern {index}"));
        }
    }
    Ok(())
}

/// Tests: dump settings for a directory (the environment is never read).
#[cfg(test)]
pub(super) fn test_settings(directory: &std::path::Path) -> Result<DumpSettings> {
    std::fs::create_dir_all(directory)?;
    Ok(DumpSettings {
        directory: directory.to_owned(),
        skip_min_p_loops: usize::MAX,
    })
}

#[cfg(test)]
mod channel_tests {
    use super::{ClosureChannel, write_channel};

    fn lines(uv_child_truncate: Option<&str>) -> String {
        let channel = ClosureChannel {
            channel: "projection metric",
            taylor_order: 0,
            normalisation: "1/(i p^2 d_A (d-1))",
            uv_child_truncate: uv_child_truncate.map(str::to_owned),
        };
        let mut out = Vec::new();
        write_channel(&mut out, &channel).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// Off: exactly the three keys of before; on: one more, last.
    #[test]
    fn the_truncation_key_is_written_only_when_on() {
        let off = "channel projection metric\ntaylor_order 0\nnormalisation 1/(i p^2 d_A (d-1))\n";
        assert_eq!(lines(None), off);
        assert_eq!(lines(Some("1@0")), format!("{off}uv_child_truncate 1@0\n"));
        assert_eq!(
            lines(Some("graded@-")),
            format!("{off}uv_child_truncate graded@-\n")
        );
    }
}

/// The sha256 (hex) of a written file, by `sha256sum` (as the block writer
/// already uses `zstd`): the record's integrity field.
fn file_sha256(path: &std::path::Path) -> Result<String> {
    let output = std::process::Command::new("sha256sum")
        .arg(path)
        .output()
        .map_err(|error| eyre!("pint6 boundary record: cannot run sha256sum: {error}"))?;
    if !output.status.success() {
        return Err(eyre!(
            "pint6 boundary record: sha256sum failed on {}",
            path.display()
        ));
    }
    let text = String::from_utf8_lossy(&output.stdout);
    text.split_whitespace()
        .next()
        .filter(|hex| hex.len() == 64)
        .map(str::to_owned)
        .ok_or_else(|| eyre!("pint6 boundary record: no sha256 for {}", path.display()))
}
