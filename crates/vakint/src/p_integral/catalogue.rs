//! The p-integral catalogue and graph matching (`docs/p_integrals.md`, section 2).
//!
//! Every non-scaleless contraction of every family is keyed by the canonical
//! form of its multigraph, with the two external vertices coloured alike (the
//! reversal `q -> -q` is a symmetry). A normalized input is matched by one
//! canonical labelling and one lookup; its loop momenta are then solved from
//! the line correspondence (`routing::input_loop_momenta_in_canonical_basis_affine`).
//!
//! The catalogue is split into strata of one loop count and one line count,
//! each built on first use. A lookup only ever builds its own stratum, so no
//! code path relies on an eager enumeration of all contractions (64 x 2^14 at
//! five loops); a shipped index can replace a stratum without changing callers.

use std::collections::{BTreeMap, HashMap};
use std::sync::{LazyLock, OnceLock};

use symbolica::atom::{Atom, AtomCore};
use symbolica::graph::Graph as NativeGraph;

use super::families;
use super::input::{PGraph, VertexClasses, is_scaleless};
use crate::VakintError;
use crate::routing::{self, linear_momentum};
use crate::symbols::S;

/// The canonical multigraph of a contraction; vertex data 1 marks the two
/// external vertices.
type Key = NativeGraph<u8, u8>;

/// One catalogue representative: the first family in table order, then the
/// lexicographically smallest set of pinched lines.
#[derive(Clone, Debug)]
pub(crate) struct CatalogueEntry {
    pub family: usize,
    /// Pinched canonical line ids, sorted.
    pub pinched: Vec<i64>,
    /// The remaining family lines with the canonical labels of their
    /// `(from, to)` vertices.
    lines: Vec<(i64, usize, usize)>,
    /// Canonical label of the vertex where the family's `p(1)` enters.
    entry: usize,
}

/// Canonical form of a multigraph with two coloured external vertices, and the
/// canonical label of each vertex.
fn canonical(lines: &[(i64, i64)], entry: i64, exit: i64) -> (Key, BTreeMap<i64, usize>) {
    let mut vertices = lines
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .chain([entry, exit])
        .collect::<Vec<_>>();
    vertices.sort_unstable();
    vertices.dedup();
    let mut graph = Key::new();
    let nodes = vertices
        .iter()
        .map(|&vertex| {
            (
                vertex,
                graph.add_node(u8::from(vertex == entry || vertex == exit)),
            )
        })
        .collect::<BTreeMap<_, _>>();
    for (a, b) in lines {
        graph
            .add_edge(nodes[a], nodes[b], false, 0)
            .expect("vertices were added above");
    }
    let form = graph.canonize();
    let labels = nodes
        .into_iter()
        .map(|(vertex, node)| (vertex, form.vertex_map[node]))
        .collect();
    (form.graph, labels)
}

/// The catalogue strata of one loop count and one number of lines.
pub(crate) struct PCatalogue {
    strata: BTreeMap<(usize, usize), OnceLock<HashMap<Key, CatalogueEntry>>>,
}

static CATALOGUE: LazyLock<PCatalogue> = LazyLock::new(|| PCatalogue {
    strata: families()
        .iter()
        .flat_map(|family| (1..=family.lines.len()).map(|lines| (family.data.loops, lines)))
        .map(|key| (key, OnceLock::new()))
        .collect(),
});

impl PCatalogue {
    pub(crate) fn get() -> &'static Self {
        &CATALOGUE
    }

    /// The stratum of graphs with `loops` loops and `lines` lines, built on first use.
    pub(crate) fn stratum(
        &self,
        loops: usize,
        lines: usize,
    ) -> Option<&HashMap<Key, CatalogueEntry>> {
        self.strata
            .get(&(loops, lines))
            .map(|stratum| stratum.get_or_init(|| build_stratum(loops, lines)))
    }

    /// Look up the representative of a normalized graph, with the canonical
    /// label of each input vertex.
    ///
    /// The labels come from a canonical labelling of the graph with its powers
    /// as edge colours, so inputs that differ only by vertex names or routing
    /// get the same line correspondence (up to the automorphisms that preserve
    /// the powers, which are symmetries of the integral).
    pub(crate) fn lookup(&self, graph: &PGraph) -> Option<(&CatalogueEntry, BTreeMap<i64, usize>)> {
        let stratum = self.stratum(graph.loops, graph.lines.len())?;
        let mut vertices = graph
            .lines
            .iter()
            .flat_map(|line| [line.ends.0, line.ends.1])
            .chain([graph.entry, graph.exit])
            .collect::<Vec<_>>();
        vertices.sort_unstable();
        vertices.dedup();
        let mut coloured = NativeGraph::<u8, i64>::new();
        let nodes = vertices
            .iter()
            .map(|&v| {
                (
                    v,
                    coloured.add_node(u8::from(v == graph.entry || v == graph.exit)),
                )
            })
            .collect::<BTreeMap<_, _>>();
        for line in &graph.lines {
            coloured
                .add_edge(nodes[&line.ends.0], nodes[&line.ends.1], false, line.power)
                .expect("vertices were added above");
        }
        let form = coloured.canonize();
        // The power-coloured canonical graph is a function of the isomorphism
        // class: label it once more without powers to find its catalogue key.
        let ends = form
            .graph
            .edges()
            .iter()
            .map(|edge| (edge.vertices.0 as i64, edge.vertices.1 as i64))
            .collect::<Vec<_>>();
        let external = form
            .graph
            .nodes()
            .iter()
            .enumerate()
            .filter(|(_, node)| node.data == 1)
            .map(|(index, _)| index as i64)
            .collect::<Vec<_>>();
        let (key, labels) = canonical(&ends, external[0], external[external.len() - 1]);
        let entry = stratum.get(&key)?;
        let labels = nodes
            .into_iter()
            .map(|(vertex, node)| (vertex, labels[&(form.vertex_map[node] as i64)]))
            .collect();
        Some((entry, labels))
    }
}

fn build_stratum(loops: usize, lines: usize) -> HashMap<Key, CatalogueEntry> {
    let started = std::time::Instant::now();
    let mut stratum = HashMap::<Key, CatalogueEntry>::new();
    for family in families()
        .iter()
        .filter(|family| family.data.loops == loops)
    {
        let count = family.lines.len();
        if lines > count {
            continue;
        }
        for mask in 0u64..1 << count {
            if mask.count_ones() as usize != count - lines {
                continue;
            }
            let mut classes = VertexClasses::default();
            let mut pinched = Vec::new();
            let mut rest = Vec::new();
            for (position, line) in family.lines.iter().enumerate() {
                if mask & (1 << position) != 0 {
                    classes.merge(line.ends);
                    pinched.push(line.id);
                } else {
                    rest.push((line.id, line.ends));
                }
            }
            let rest = rest
                .into_iter()
                .map(|(id, (a, b))| (id, (classes.find(a), classes.find(b))))
                .collect::<Vec<_>>();
            let (entry, exit) = (classes.find(1), classes.find(2));
            let ends = rest.iter().map(|(_, ends)| *ends).collect::<Vec<_>>();
            if is_scaleless(&ends, entry, exit, loops) {
                continue;
            }
            let (key, labels) = canonical(&ends, entry, exit);
            let candidate = CatalogueEntry {
                family: family.index,
                lines: rest
                    .iter()
                    .map(|(id, (a, b))| (*id, labels[a], labels[b]))
                    .collect(),
                entry: labels[&entry],
                pinched,
            };
            match stratum.entry(key) {
                std::collections::hash_map::Entry::Vacant(slot) => {
                    slot.insert(candidate);
                }
                std::collections::hash_map::Entry::Occupied(mut slot) => {
                    let current = slot.get();
                    if current.family == candidate.family && candidate.pinched < current.pinched {
                        slot.insert(candidate);
                    }
                }
            }
        }
    }
    log::debug!(target: "vakint::p_integral", "phase=catalogue_stratum loops={loops} lines={lines} graphs={} elapsed_seconds={:.6}", stratum.len(), started.elapsed().as_secs_f64());
    stratum
}

/// A verified identification of a normalized p-integral with a family sector.
#[derive(Clone, Debug)]
pub(crate) struct PMatch {
    pub family: usize,
    /// The canonical `p(1)` is `sign·p_user`.
    pub sign: i8,
    /// The power of every canonical line, in canonical order; 0 for pinched lines.
    pub powers: Vec<i64>,
    /// `(canonical line id, input prop id)` of every remaining line.
    pub edges: Vec<(i64, i64)>,
    /// `(input k(j), affine image in the canonical k(l) and p_user)`.
    pub images: Vec<(Atom, Atom)>,
}

/// Upper bound on the assignments of parallel lines tried for one input.
/// Any assignment is an isomorphism, so the first one routes; the bound only
/// guards against pathological multigraphs.
const MAXIMAL_PARALLEL_ASSIGNMENTS: usize = 1 << 12;

/// Match a normalized p-integral graph to its catalogue representative.
/// Returns `None` when its graph is in no catalogue stratum.
pub(crate) fn match_graph(graph: &PGraph) -> Result<Option<PMatch>, VakintError> {
    let Some((entry, labels)) = PCatalogue::get().lookup(graph) else {
        return Ok(None);
    };
    let family = &families()[entry.family];
    let sign = if labels[&graph.entry] == entry.entry {
        1
    } else {
        -1
    };
    let powers = |assignment: &[(i64, i64)]| {
        family
            .lines
            .iter()
            .map(|line| {
                assignment
                    .iter()
                    .find(|(canonical, _)| *canonical == line.id)
                    .map_or(0, |(_, input)| {
                        graph.lines.iter().find(|l| l.id == *input).unwrap().power
                    })
            })
            .collect::<Vec<_>>()
    };
    // Pools of parallel lines, keyed by their unordered canonical endpoints.
    let mut pools = BTreeMap::<(usize, usize), (Vec<usize>, Vec<(i64, usize, usize)>)>::new();
    for (index, line) in graph.lines.iter().enumerate() {
        let (a, b) = (labels[&line.ends.0], labels[&line.ends.1]);
        pools.entry((a.min(b), a.max(b))).or_default().0.push(index);
    }
    for &(id, a, b) in &entry.lines {
        pools
            .entry((a.min(b), a.max(b)))
            .or_default()
            .1
            .push((id, a, b));
    }
    // Lines of one pool are assigned in order of their powers, so that the
    // correspondence only depends on the power-coloured isomorphism class.
    let pools = pools
        .into_values()
        .map(|(mut inputs, targets)| {
            inputs.sort_by_key(|&index| (graph.lines[index].power, index));
            (inputs, targets)
        })
        .collect::<Vec<_>>();
    let canonical_external = S.p.call(1);
    let mut attempts = 0;
    let mut failure = None;
    let mut assignment = vec![None; graph.lines.len()];
    let routed = assign(&pools, 0, &mut assignment, &mut |assignment| {
        if attempts >= MAXIMAL_PARALLEL_ASSIGNMENTS {
            return Some(None);
        }
        attempts += 1;
        let equations = graph
            .lines
            .iter()
            .zip(assignment)
            .map(|(line, target)| {
                let (id, a, b) = target.expect("complete assignment");
                let momentum = &family.lines.iter().find(|l| l.id == id).unwrap().momentum;
                let forward = (labels[&line.ends.0], labels[&line.ends.1]) == (a, b);
                let image = if forward { momentum.clone() } else { -momentum };
                (line.momentum.clone(), image)
            })
            .collect::<Vec<_>>();
        // An input already written in the family's routing and orientations
        // keeps its loop momenta: the solve would return the identity.
        let oriented = Atom::num(sign as i64) * &graph.external;
        if equations.iter().all(|(input, image)| {
            let image = image
                .replace(canonical_external.to_pattern())
                .with(oriented.to_pattern());
            matches!(
                (linear_momentum(input.as_view()), linear_momentum(image.as_view())),
                (Ok(input), Ok(image)) if input == image
            )
        }) {
            return Some(Some(routing::AffineRouting {
                images: (1..=graph.loops)
                    .map(|axis| (S.k.call(axis), S.k.call(axis)))
                    .collect(),
                external_sign: sign,
            }));
        }
        match routing::input_loop_momenta_in_canonical_basis_affine(
            &equations,
            Some((&graph.external, &canonical_external, sign)),
        ) {
            Ok(routing) => Some(Some(routing)),
            Err(error) => {
                failure = Some(error);
                None
            }
        }
    })
    .flatten();
    let Some(routing) = routed else {
        return Err(failure.unwrap_or_else(|| {
            VakintError::InvalidMomentumExpression(
                "p-integral: no assignment of parallel lines routes the matched family".into(),
            )
        }));
    };
    let edges = graph
        .lines
        .iter()
        .zip(&assignment)
        .map(|(line, target)| (target.expect("complete assignment").0, line.id))
        .collect::<Vec<_>>();
    Ok(Some(PMatch {
        family: entry.family,
        sign: routing.external_sign,
        powers: powers(&edges),
        edges,
        images: routing.images,
    }))
}

/// Enumerate the assignments of input lines to canonical lines within each
/// pool of parallel lines, stopping at the first `Some` of `visit`.
fn assign<T>(
    pools: &[(Vec<usize>, Vec<(i64, usize, usize)>)],
    pool: usize,
    assignment: &mut Vec<Option<(i64, usize, usize)>>,
    visit: &mut impl FnMut(&[Option<(i64, usize, usize)>]) -> Option<T>,
) -> Option<T> {
    let Some((inputs, targets)) = pools.get(pool) else {
        return visit(assignment);
    };
    debug_assert_eq!(inputs.len(), targets.len());
    let mut order = (0..targets.len()).collect::<Vec<_>>();
    loop {
        for (input, target) in inputs.iter().zip(&order) {
            assignment[*input] = Some(targets[*target]);
        }
        if let Some(found) = assign(pools, pool + 1, assignment, visit) {
            return Some(found);
        }
        if !next_permutation(&mut order) {
            return None;
        }
    }
}

/// Advance to the next lexicographic permutation; `false` after the last one.
fn next_permutation(order: &mut [usize]) -> bool {
    let Some(pivot) = (1..order.len()).rev().find(|&i| order[i - 1] < order[i]) else {
        return false;
    };
    let successor = (pivot..order.len())
        .rev()
        .find(|&i| order[i] > order[pivot - 1])
        .unwrap();
    order.swap(pivot - 1, successor);
    order[pivot..].reverse();
    true
}
