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
use vakint::{VakintSettings, p_integral::MassiveClosure, symbols::S};

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
            incident.len() == 2
                && incident[0] != incident[1]
                && !external_vertices.contains(vertex)
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
                .add_edge(index[a], index[b], false, if with_powers { *power } else { 1 })
                .map_err(|error| eyre!(error))?;
        }
        Ok(graph.canonize().graph.to_string().trim().replace('\n', "; "))
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
            visit(*vertex, None, &adjacency, &mut time, &mut discovery, &mut low, &mut bridges);
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

/// Write one closure record. Returns the remaining p-integral loop count.
pub(super) fn write(
    settings: &DumpSettings,
    vakint_settings: &VakintSettings,
    closure: &MassiveClosure,
    numerator: &Atom,
    integrand: &Atom,
    vacuum_loops: usize,
) -> Result<usize> {
    let p_loops = vacuum_loops.saturating_sub(1);
    let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    if SETTINGS_WRITTEN.set(()).is_ok() {
        std::fs::write(
            settings.directory.join(format!("vakint-settings-{pid}.txt")),
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
    let mut out = BufWriter::new(File::create(settings.directory.join(format!("{stem}.txt")))?);
    writeln!(out, "schema pint6.massive-closure-request.v1")?;
    writeln!(out, "sequence {sequence}")?;
    writeln!(out, "context {context}")?;
    writeln!(out, "vacuum_loops {vacuum_loops}")?;
    writeln!(out, "p_loops {p_loops}")?;
    writeln!(out, "largest_block_loops {}", largest_block_loops(&lines))?;
    writeln!(out, "evaluated {}", p_loops < settings.skip_min_p_loops)?;
    writeln!(out, "closure_momentum {}", closure.momentum.to_canonical_string())?;
    writeln!(out, "closure_mass_squared {}", closure.mass_squared.to_canonical_string())?;
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
    writeln!(out, "numerator_bytes {}", numerator.as_view().get_byte_size())?;
    // Exactly the Atom MassiveClosure::evaluate_batch receives with this closure.
    writeln!(out, "input {}", (numerator * integrand).to_canonical_string())?;
    out.flush()?;
    let mut tsv =
        BufWriter::new(File::create(settings.directory.join(format!("{stem}.monomials.tsv")))?);
    for (monomial, _) in &monomials {
        writeln!(
            tsv,
            "{}\t{}",
            loop_rank(monomial.as_view()),
            monomial.to_canonical_string()
        )?;
    }
    tsv.flush()?;
    Ok(p_loops)
}
