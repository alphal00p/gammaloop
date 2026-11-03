//! Catalogue, normalization and census tests of the p-integral frontend.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use symbolica::atom::{Atom, AtomCore};
use symbolica::graph::Graph as NativeGraph;

use super::catalogue::{PCatalogue, match_graph};
use super::input::{Normalized, PGraph, PLine, VertexClasses, is_scaleless};
use super::{P_FAMILIES, families};
use crate::symbols::S;

/// A routing of an undirected multigraph: line `j` of a spanning-tree
/// complement carries `k(j)`, the tree lines carry what conservation requires
/// with `p(1)` entering at `entry` and leaving at `exit`.
fn routed(ends: &[(i64, i64)], entry: i64, exit: i64) -> PGraph {
    let vertices = ends
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .chain([entry, exit])
        .collect::<BTreeSet<_>>();
    // Breadth-first spanning tree from `entry`.
    let mut parent = BTreeMap::<i64, (i64, usize)>::new();
    let mut order = vec![entry];
    let mut seen = BTreeSet::from([entry]);
    let mut position = 0;
    while position < order.len() {
        let vertex = order[position];
        position += 1;
        for (index, &(a, b)) in ends.iter().enumerate() {
            let other = if a == vertex {
                b
            } else if b == vertex {
                a
            } else {
                continue;
            };
            if seen.insert(other) {
                parent.insert(other, (vertex, index));
                order.push(other);
            }
        }
    }
    assert_eq!(seen, vertices, "the graph must be connected");
    let tree = parent
        .values()
        .map(|(_, index)| *index)
        .collect::<HashSet<_>>();
    let mut momenta = vec![Atom::Zero; ends.len()];
    // Net momentum that must leave each vertex through the tree.
    let mut excess = vertices
        .iter()
        .map(|&v| {
            let external = if v == entry {
                S.p.call(1)
            } else if v == exit {
                -S.p.call(1)
            } else {
                Atom::Zero
            };
            (v, external)
        })
        .collect::<BTreeMap<_, _>>();
    let mut loop_axis = 0;
    for (index, &(a, b)) in ends.iter().enumerate() {
        if tree.contains(&index) {
            continue;
        }
        loop_axis += 1;
        momenta[index] = S.k.call(loop_axis);
        *excess.get_mut(&a).unwrap() -= &momenta[index];
        *excess.get_mut(&b).unwrap() += &momenta[index];
    }
    // Leaves first: a tree line carries the subtree's excess towards the root.
    for &vertex in order.iter().rev() {
        if let Some(&(up, index)) = parent.get(&vertex) {
            let flow = excess[&vertex].clone();
            // The line is written ends[index]; flow goes vertex -> up.
            momenta[index] = if ends[index] == (vertex, up) {
                flow.clone()
            } else {
                -flow.clone()
            };
            *excess.get_mut(&up).unwrap() += flow;
        }
    }
    PGraph {
        lines: ends
            .iter()
            .zip(momenta)
            .enumerate()
            .map(|(index, (&ends, momentum))| PLine {
                id: index as i64 + 1,
                ends,
                momentum: momentum.expand(),
                power: 1,
            })
            .collect(),
        external: S.p.call(1),
        entry,
        exit,
        loops: loop_axis,
    }
}

/// Loopless connected multigraphs with external vertices 0 and 1 of degree 2
/// and every other vertex of degree 3: the cubic p-graphs with `2L` vertices.
fn cubic_p_graphs(loops: usize) -> Vec<Vec<(i64, i64)>> {
    let vertices = 2 * loops;
    let mut remaining = (0..vertices)
        .map(|v| if v < 2 { 2 } else { 3 })
        .collect::<Vec<_>>();
    let mut graphs = Vec::new();
    let mut seen = HashSet::new();
    let mut edges = Vec::new();
    fn extend(
        pair: (usize, usize),
        remaining: &mut [usize],
        edges: &mut Vec<(i64, i64)>,
        graphs: &mut Vec<Vec<(i64, i64)>>,
        seen: &mut HashSet<NativeGraph<u8, u8>>,
    ) {
        let n = remaining.len();
        let (i, j) = pair;
        if i == n {
            if remaining.iter().all(|&r| r == 0) {
                let mut graph = NativeGraph::<u8, u8>::new();
                for v in 0..n {
                    graph.add_node(u8::from(v < 2));
                }
                for &(a, b) in edges.iter() {
                    graph.add_edge(a as usize, b as usize, false, 0).unwrap();
                }
                let classes = edges.iter().fold(VertexClasses::default(), |mut c, &e| {
                    c.merge(e);
                    c
                });
                let mut classes = classes;
                let connected = (0..n as i64)
                    .map(|v| classes.find(v))
                    .collect::<BTreeSet<_>>()
                    .len()
                    == 1;
                if connected && seen.insert(graph.canonize().graph) {
                    graphs.push(edges.clone());
                }
            }
            return;
        }
        if j == n {
            if remaining[i] == 0 {
                extend((i + 1, i + 2), remaining, edges, graphs, seen);
            }
            return;
        }
        for multiplicity in 0..=remaining[i].min(remaining[j]) {
            remaining[i] -= multiplicity;
            remaining[j] -= multiplicity;
            for _ in 0..multiplicity {
                edges.push((i as i64, j as i64));
            }
            extend((i, j + 1), remaining, edges, graphs, seen);
            for _ in 0..multiplicity {
                edges.pop();
            }
            remaining[i] += multiplicity;
            remaining[j] += multiplicity;
        }
    }
    extend((0, 1), &mut remaining, &mut edges, &mut graphs, &mut seen);
    graphs
}

/// Section 7 census: every connected loopless cubic p-graph with `loops` loops
/// and every contraction of it is either scaleless or matched (and routed) by
/// the catalogue after normalization. Returns (graphs, scaleless, matched).
fn census(loops: usize) -> (usize, usize, usize) {
    let mut seen = HashSet::new();
    let (mut graphs, mut scaleless, mut matched) = (0, 0, 0);
    for ends in cubic_p_graphs(loops) {
        graphs += 1;
        let routed = routed(&ends, 0, 1);
        for mask in 0u32..1 << ends.len() {
            let mut graph = routed.clone();
            for (position, line) in graph.lines.iter_mut().enumerate() {
                if mask & (1 << position) != 0 {
                    line.power = 0;
                }
            }
            let normalized = graph.normalize().unwrap();
            let Normalized::Graph { graph, .. } = normalized else {
                scaleless += 1;
                continue;
            };
            let mut key = NativeGraph::<u8, u8>::new();
            let mut vertices = graph
                .lines
                .iter()
                .flat_map(|line| [line.ends.0, line.ends.1])
                .collect::<Vec<_>>();
            vertices.sort_unstable();
            vertices.dedup();
            for &v in &vertices {
                key.add_node(u8::from(v == graph.entry || v == graph.exit));
            }
            let index = |v: i64| vertices.iter().position(|&w| w == v).unwrap();
            for line in &graph.lines {
                key.add_edge(index(line.ends.0), index(line.ends.1), false, 0)
                    .unwrap();
            }
            if !seen.insert(key.canonize().graph) {
                continue;
            }
            let found = match_graph(&graph).unwrap();
            assert!(
                found.is_some(),
                "{loops}L graph {ends:?} with lines {mask:b} pinched is not covered: {}",
                graph.to_atom()
            );
            matched += 1;
        }
    }
    println!(
        "{loops}L census: {graphs} cubic p-graphs, {scaleless} scaleless contractions, {matched} distinct matched graphs"
    );
    (graphs, scaleless, matched)
}

#[test]
fn every_one_to_three_loop_p_graph_is_covered_after_normalization() {
    for loops in 1..=3 {
        let (graphs, _, matched) = census(loops);
        assert!(graphs > 0 && matched > 0, "{loops}L census is empty");
    }
}

#[test]
fn every_four_loop_p_graph_is_covered_after_normalization() {
    let (graphs, _, matched) = census(4);
    assert!(graphs > 0 && matched > 0);
}

#[test]
fn catalogue_entries_are_the_first_family_and_smallest_pinch_set() {
    for loops in 1..=3 {
        let mut expected =
            BTreeMap::<(usize, usize), Vec<(Vec<(i64, i64)>, usize, Vec<i64>)>>::new();
        for family in families()
            .iter()
            .filter(|family| family.data.loops == loops)
        {
            let count = family.lines.len();
            for mask in 0u32..1 << count {
                let mut classes = VertexClasses::default();
                let mut pinched = Vec::new();
                for (position, line) in family.lines.iter().enumerate() {
                    if mask & (1 << position) != 0 {
                        classes.merge(line.ends);
                        pinched.push(line.id);
                    }
                }
                let rest = family
                    .lines
                    .iter()
                    .filter(|line| !pinched.contains(&line.id))
                    .map(|line| (classes.find(line.ends.0), classes.find(line.ends.1)))
                    .collect::<Vec<_>>();
                let (entry, exit) = (classes.find(1), classes.find(2));
                if is_scaleless(&rest, entry, exit, loops) {
                    continue;
                }
                expected.entry((loops, rest.len())).or_default().push((
                    rest,
                    family.index,
                    pinched,
                ));
            }
        }
        for ((loops, lines), contractions) in expected {
            let stratum = PCatalogue::get().stratum(loops, lines).unwrap();
            for (rest, family, pinched) in contractions {
                let (entry, _) = PCatalogue::get().lookup(&routed(&rest, 1, 2)).unwrap();
                // The representative precedes every contraction of its class.
                assert!(
                    (entry.family, &entry.pinched) <= (family, &pinched),
                    "{} {:?} precedes the stored {} {:?}",
                    P_FAMILIES[family].name,
                    pinched,
                    P_FAMILIES[entry.family].name,
                    entry.pinched
                );
            }
            assert!(!stratum.is_empty());
        }
    }
}

#[test]
fn normalization_merges_identical_quadrics_and_contracts_tree_lines() {
    crate::Vakint::initialize_vakint_symbols();
    // A bubble with a dotted line split in two by a 2-valent vertex, and a
    // tree line carrying q into the graph.
    let integral = crate::utils::vakint_macros::vk_parse!(
        "topo(prop(1,edge(5,1),p(1),0,2)*prop(2,edge(1,3),k(1),0,1)*prop(3,edge(3,2),k(1),0,2)\
         *prop(4,edge(2,4),k(1)-p(1),0,-1)*prop(5,edge(4,1),k(1)-p(1),0,2)*prop(6,edge(2,6),p(1),0,0))"
    )
    .unwrap();
    let Normalized::Graph { graph, factor } = PGraph::parse(integral.as_view())
        .unwrap()
        .normalize()
        .unwrap()
    else {
        panic!("the bubble is not scaleless");
    };
    let lines = graph
        .lines
        .iter()
        .map(|line| (line.id, line.power))
        .collect::<Vec<_>>();
    assert_eq!(lines, [(2, 3), (4, 1)]);
    assert_eq!(factor, S.dot(S.p.call(1), S.p.call(1)).pow(Atom::num(-2)));
    let matched = match_graph(&graph).unwrap().unwrap();
    assert_eq!(P_FAMILIES[matched.family].name, "P1L");
    let mut powers = matched.powers.clone();
    powers.sort();
    assert_eq!(powers, [1, 3]);
}

#[test]
fn the_witness_supplies_the_backend_request() {
    let vakint = crate::Vakint::new().unwrap();
    let parse = |input: &str| crate::utils::vakint_macros::vk_parse!(input).unwrap();
    // A dotted bubble in the user's routing k(3), with the head p(2).
    let integral = parse("topo(prop(1,edge(1,2),k(3)+p(2),0,2)*prop(2,edge(2,1),k(3),0,1))");
    let rules = vakint
        .topologies
        .match_topologies_to_user_input(integral.as_view(), false)
        .unwrap()
        .unwrap();
    let mut term = crate::VakintTerm {
        integral,
        numerator: parse("dot(k(3),p(2))^2"),
        vectors: Vec::new(),
    };
    term.canonicalize(&crate::VakintSettings::default(), &rules, false)
        .unwrap();
    let request = rules.propagator_request(term.numerator.as_view()).unwrap();
    assert_eq!(request.family.name, "P1L");
    let mut powers = request.powers.to_vec();
    powers.sort();
    assert_eq!(powers, [1, 2]);
    let p2 = parse("p(2)");
    assert!(request.external == p2 || request.external == -&p2);
    // The canonical numerator is a polynomial in canonical scalar products.
    let numerator = term.numerator.to_string();
    assert!(
        numerator.contains("k(1)") && !numerator.contains("k(3)"),
        "{numerator}"
    );
    assert_eq!(
        super::scalar_numerator(term.numerator.as_view(), &p2).unwrap(),
        term.numerator
    );
}
