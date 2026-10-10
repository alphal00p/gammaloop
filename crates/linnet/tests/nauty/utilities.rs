//! The utility pipelines of `runalltests` whose outputs are isomorphism-class data: canonical
//! forms, iso-rejected families and their counts. Each test reads the run's answer file (or the
//! oracle-recorded inputs of its pipeline, see the fixture README) and checks that Linnet's
//! canonical forms separate exactly nauty's classes, with nauty's group sizes and orbits.
//!
//! Families that nauty iso-rejects per base graph (orientations, colourings, multiplicities) are
//! rebuilt by brute force: every labelled variant of a base is canonized, and each class must be
//! hit exactly `|Aut(base)| / |Aut(variant)|` times. That orbit-stabilizer count holds only if
//! canonical forms are complete invariants and automorphism counts are exact.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use linnet::half_edge::generate::{GraphGenerator, Stub};
use rand::{rngs::SmallRng, seq::SliceRandom, SeedableRng};

use super::{
    formats::read_graphs,
    graph::{Form, Graph},
    runs::{fixture, full_size},
};

fn rng() -> SmallRng {
    SmallRng::seed_from_u64(2026)
}

/// Canonize labelled variants of one base graph, whose automorphism group has order
/// `base_order`, and check orbit-stabilizer per class; `stabiliser` turns Linnet's group order of
/// a variant into the order of its stabiliser in the base group. Returns the number of classes.
fn classify(
    variants: impl Iterator<Item = Graph>,
    base_order: u128,
    stabiliser: impl Fn(&Graph, u128) -> u128,
) -> usize {
    let mut classes: HashMap<Form, (u128, u128)> = HashMap::new();
    for variant in variants {
        let canon = variant.canonize();
        let stabiliser = stabiliser(&variant, canon.automorphism_count().unwrap());
        assert_eq!(base_order % stabiliser, 0, "{variant:?}");
        let orbit = base_order / stabiliser;
        let (hits, expected) = classes
            .entry(canon.into_canonical_form())
            .or_insert((0, orbit));
        assert_eq!(
            *expected, orbit,
            "isomorphic variants with different groups"
        );
        *hits += 1;
    }
    for (hits, orbit) in classes.values() {
        assert_eq!(hits, orbit, "labelled variants per class");
    }
    classes.len()
}

fn same_order(_: &Graph, order: u128) -> u128 {
    order
}

/// Every orientation of the edges of `base`.
fn orientations(base: &Graph) -> impl Iterator<Item = Graph> + '_ {
    (0..1u64 << base.edges.len()).map(move |mask| {
        let arcs = base
            .edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| if mask >> i & 1 == 1 { (b, a) } else { (a, b) })
            .collect();
        Graph::directed(base.n, arcs)
    })
}

fn distinct(forms: impl IntoIterator<Item = Form>) -> usize {
    forms.into_iter().collect::<BTreeSet<_>>().len()
}

/// Runs 63-65, `geng -cd1f 10 | labelg -q[S|t]`: the 3389 connected 4-cycle-free graphs on 10
/// vertices canonically labelled by dense nauty, sparse nauty and Traces. Line `i` of the three
/// answers is one class, all classes differ, and canonical forms ignore the labelling.
#[test]
fn labelg_canonical_labellings() {
    let files = ["nautestb1.ans", "nautestb2.ans", "nautestb3.ans"].map(|f| {
        read_graphs(&fixture(f))
            .iter()
            .map(Graph::canonical_form)
            .collect::<Vec<_>>()
    });
    assert_eq!(files[0].len(), 3389);
    assert_eq!(files[0], files[1]);
    assert_eq!(files[0], files[2]);
    assert_eq!(distinct(files[0].clone()), 3389);
    let mut rng = rng();
    for (g, form) in read_graphs(&fixture("nautestb1.ans")).iter().zip(&files[0]) {
        g.check_invariance(form, 1, &mut rng);
    }
}

/// Run 67, `gentourng -q -z 8 | pickg -T3`: 14 tournaments on 8 vertices, with the group sizes
/// countg reports for them. All 2^15 labelled tournaments on 6 vertices (2^21 on 7 in
/// full-size runs, see [`full_size`]) fall into the known 56 (456) classes.
#[test]
fn gentourng_tournaments() {
    let tournaments = read_graphs(&fixture("nautestd.ans"));
    let orders: Vec<u128> = tournaments
        .iter()
        .map(|t| {
            assert!(t.directed && t.n == 8 && t.edges.len() == 28);
            let symmetry = t.symmetry();
            t.check_invariance(&symmetry.canonical_form, 2, &mut rng());
            symmetry.order
        })
        .collect();
    assert_eq!(orders, [1, 1, 1, 1, 1, 3, 1, 3, 3, 3, 1, 3, 3, 1]);
    assert_eq!(distinct(tournaments.iter().map(Graph::canonical_form)), 14);
    let sizes: &[(usize, usize, u128)] = if full_size() {
        &[(6, 56, 720), (7, 456, 5040)]
    } else {
        &[(6, 56, 720)]
    };
    for &(n, classes, factorial) in sizes {
        let complete = Graph::undirected(
            n,
            (0..n).flat_map(|b| (0..b).map(move |a| (a, b))).collect(),
        );
        assert_eq!(
            classify(orientations(&complete), factorial, same_order),
            classes
        );
    }
}

/// Free trees on `n >= 2` vertices, by the algorithm of Wright, Richmond, Odlyzko and McKay
/// (1986) that `gentreeg` implements: level sequences of trees rooted at a centre, stepped with
/// Beyer-Hedetniemi and skipping non-canonical roots.
fn free_trees(n: usize) -> Vec<Graph> {
    fn next_rooted(levels: &[usize], p: Option<usize>) -> Option<Vec<usize>> {
        let p = p.unwrap_or_else(|| levels.iter().rposition(|&l| l != 1).unwrap());
        if p == 0 {
            return None;
        }
        let q = (0..p).rev().find(|&q| levels[q] == levels[p] - 1).unwrap();
        let mut next = levels.to_vec();
        for i in p..next.len() {
            next[i] = next[i - p + q];
        }
        Some(next)
    }
    /// The subtree of the first child of the root, one level up, and the rest of the tree.
    fn split(levels: &[usize]) -> (Vec<usize>, Vec<usize>) {
        let m = (2..levels.len())
            .find(|&i| levels[i] == 1)
            .unwrap_or(levels.len());
        let left = levels[1..m].iter().map(|l| l - 1).collect();
        let rest = [0].into_iter().chain(levels[m..].iter().copied()).collect();
        (left, rest)
    }
    fn next_free(levels: Vec<usize>) -> Option<Vec<usize>> {
        let (left, rest) = split(&levels);
        let (lh, rh) = (*left.iter().max().unwrap(), *rest.iter().max().unwrap());
        let valid = rh > lh
            || (rh == lh
                && (left.len() < rest.len() || (left.len() == rest.len() && left <= rest)));
        if valid {
            return Some(levels);
        }
        let p = left.len();
        let mut next = next_rooted(&levels, Some(p))?;
        if levels[p] > 2 {
            let height = *split(&next).0.iter().max().unwrap();
            let len = next.len();
            for (k, l) in (1..=height + 1).enumerate() {
                next[len - (height + 1) + k] = l;
            }
        }
        Some(next)
    }
    let edges = |levels: &[usize]| {
        let mut stack: Vec<usize> = vec![];
        let mut edges = vec![];
        for (i, &level) in levels.iter().enumerate() {
            while stack.last().is_some_and(|&j| levels[j] >= level) {
                stack.pop();
            }
            if let Some(&j) = stack.last() {
                edges.push((j, i));
            }
            stack.push(i);
        }
        Graph::undirected(levels.len(), edges)
    };
    let mut trees = vec![];
    let mut levels: Option<Vec<usize>> = Some((0..=n / 2).chain(1..n.div_ceil(2)).collect());
    while let Some(candidate) = levels {
        let Some(tree) = next_free(candidate) else {
            break;
        };
        trees.push(edges(&tree));
        levels = next_rooted(&tree, None);
    }
    trees
}

/// Run 68, `gentreeg -q -s 18 | countg -z6 --o`: the 123867 trees on 18 vertices, those of
/// radius 6 tabulated by their number of vertex orbits. The trees are enumerated here rather than
/// generated (see below); reduced runs canonize only the 12869 trees of radius 6.
#[test]
fn gentreeg_orbit_table() {
    let answer = fixture("nauteste.ans");
    let mut expected = BTreeMap::new();
    for line in answer.lines() {
        if let Some((count, orbits)) = line.split_once(" graphs : orbits=") {
            expected.insert(
                orbits.parse::<usize>().unwrap(),
                count.trim().parse::<usize>().unwrap(),
            );
        }
    }
    assert!(answer.contains(" 12869 graphs altogether from 123867 read"));
    let trees = free_trees(18);
    assert_eq!(trees.len(), 123867);
    let full = full_size();
    let mut forms = BTreeSet::new();
    let mut table = BTreeMap::new();
    for (i, tree) in trees.iter().enumerate() {
        let radius_six = tree.radius_diameter().0 == 6;
        if !radius_six && !full {
            continue;
        }
        let canon = tree.canonize();
        if radius_six {
            let orbits = canon.node_orbits();
            let roots: BTreeSet<_> = orbits.iter().map(|(_, r)| *r).collect();
            *table.entry(roots.len()).or_insert(0) += 1;
            if i % 64 == 0 {
                let symmetry = tree.symmetry();
                assert_eq!(symmetry.node_orbits.len(), roots.len());
            }
        }
        forms.insert(canon.into_canonical_form());
    }
    let classes = if full { 123867 } else { 12869 };
    assert_eq!(forms.len(), classes);
    assert_eq!(table, expected);
}

/// Linnet's generator against the same enumeration: trees from one signature per degree, up to 9
/// vertices (11 in full-size runs). The search grows much faster than the number of trees when the
/// signatures span many degrees (up to 14 vertices with degrees 1 to 15 take minutes), so the 18
/// vertices of run 68 are out of its reach.
#[test]
fn generated_trees_match_the_free_tree_enumeration() {
    let n = if full_size() { 11 } else { 9 };
    // OEIS A000055, from 0 vertices
    let counts = [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235];
    let undirected = Stub::undirected(0u8);
    let generator = GraphGenerator::<u8, u8>::new([], (1..n).map(|d| vec![undirected.clone(); d]))
        .loops(0..=0)
        .max_vertices(n)
        .allow_zero_flow_edges(true);
    let mut generated: BTreeMap<usize, BTreeSet<Form>> = BTreeMap::new();
    for diagram in generator.iter().unwrap() {
        let tree = Graph::of_hedge_graph(diagram.graph());
        let symmetry = tree.symmetry();
        assert_eq!(diagram.automorphism_count, Some(symmetry.order));
        assert!(generated
            .entry(tree.n)
            .or_default()
            .insert(symmetry.canonical_form));
    }
    for (k, count) in counts.iter().enumerate().take(n + 1).skip(2) {
        let trees: BTreeSet<Form> = free_trees(k).iter().map(Graph::canonical_form).collect();
        assert_eq!(trees.len(), *count);
        assert_eq!(generated[&k], trees, "trees on {k} vertices");
    }
}

/// Run 69, `geng -c -qp 8 9 | directg -q -o | converseg -q | labelg -q`: the 36880 orientations
/// (no digons) of the 153 connected 5-cycle-free graphs with 8 vertices and 9 edges. Taking
/// converses permutes the orientations of each base, so the class count is unchanged.
#[test]
fn directg_orientations() {
    let bases = read_graphs(&fixture("directg-bases.g6"));
    assert_eq!(bases.len(), 153);
    let classes: usize = bases
        .iter()
        .map(|base| {
            let order = base.symmetry().order;
            let classes = classify(orientations(base), order, same_order);
            let converses = classify(orientations(base).map(|g| g.converse()), order, same_order);
            assert_eq!(classes, converses);
            classes
        })
        .sum();
    // the number of lines of nautestf.ans
    assert_eq!(classes, 36880);
}

/// Run 70, `genbg -q 4 3 10 | vcolg -m5 -T`: vertex colourings with 5 colours, up to
/// isomorphism, of the three bipartite graphs with parts of 4 and 3 vertices and 10 edges. Each
/// base takes 78125 canonizations; reduced runs check the first one only.
#[test]
fn vcolg_colourings() {
    let bases = read_graphs(&fixture("vcolg-bases.g6"));
    let checked = if full_size() { bases.len() } else { 1 };
    let mut classes = vec![];
    for base in &bases[..checked] {
        let order = base.symmetry().order;
        let colourings = (0..5u32.pow(base.n as u32)).map(|code| {
            let colours = (0..base.n).map(|v| code / 5u32.pow(v as u32) % 5).collect();
            base.clone().with_colours(colours)
        });
        classes.push((order, classify(colourings, order, same_order)));
    }
    // countg's group sizes and the colourings vcolg writes per base, 54375 lines in all
    assert_eq!(classes, [(8, 16875), (12, 13125), (4, 24375)][..checked]);
}

/// Run 71, `genposetg 8 o q`: the 16999 posets on 8 points as Hasse diagrams. They are pairwise
/// non-isomorphic, and their group sizes add up to the 431723379 labelled posets on 8 points
/// (OEIS A001035) by orbit-stabilizer: sum of 8!/|Aut(P)|.
#[test]
fn genposetg_posets() {
    let posets = read_graphs(&fixture("nautesth.ans"));
    assert_eq!(posets.len(), 16999);
    let mut rng = rng();
    let mut labelled = 0;
    let mut forms = BTreeSet::new();
    for (i, poset) in posets.iter().enumerate() {
        assert!(poset.directed && poset.n == 8);
        let canon = poset.canonize();
        labelled += 40320 / canon.automorphism_count().unwrap();
        let form = canon.into_canonical_form();
        if i % 16 == 0 {
            poset.check_invariance(&form, 1, &mut rng);
        }
        forms.insert(form);
    }
    assert_eq!(forms.len(), 16999);
    assert_eq!(labelled, 431723379);
}

/// Run 72, `genspecialg -q -c12 -b11,4 -P5,2 -T3,4,5 -X30`: a 12-cycle, K(11,4), the Petersen
/// graph, the theta graph of paths of lengths 3, 4 and 5, and a (4,7)-cage, with countg's
/// `|Aut|`, vertex orbits and edge orbits. The first four are also rebuilt here.
#[test]
fn genspecialg_graphs() {
    let graphs = read_graphs(&fixture("nautesti.ans"));
    let expected: [(u128, usize, usize); 5] = [
        (24, 1, 1),
        (958003200, 2, 1),
        (120, 1, 1),
        (2, 6, 7),
        (4, 23, 41),
    ];
    let mut rng = rng();
    for (g, &(order, orbits, edge_orbits)) in graphs.iter().zip(&expected) {
        let symmetry = g.symmetry();
        assert_eq!(
            (
                symmetry.order,
                symmetry.node_orbits.len(),
                symmetry.edge_orbits
            ),
            (order, orbits, edge_orbits)
        );
        g.check_invariance(&symmetry.canonical_form, 3, &mut rng);
    }
    let cycle = Graph::undirected(12, (0..12).map(|i| (i, (i + 1) % 12)).collect());
    let bipartite = Graph::undirected(
        15,
        (0..11)
            .flat_map(|a| (11..15).map(move |b| (a, b)))
            .collect(),
    );
    let petersen = Graph::undirected(
        10,
        (0..5)
            .flat_map(|i| [(i, (i + 1) % 5), (i, i + 5), (i + 5, (i + 2) % 5 + 5)])
            .collect(),
    );
    for (built, g) in [cycle, bipartite, petersen, theta(&[3, 4, 5], false)]
        .iter()
        .zip(&graphs)
    {
        assert_eq!(built.canonical_form(), g.canonical_form());
    }
}

/// Paths of the given lengths between two hubs, directed from hub 0 to hub 1 if `directed`.
fn theta(lengths: &[usize], directed: bool) -> Graph {
    let mut n = 2;
    let mut edges = vec![];
    for &length in lengths {
        let inner: Vec<usize> = (n..n + length - 1).collect();
        n += length - 1;
        let path: Vec<usize> = [0].into_iter().chain(inner).chain([1]).collect();
        edges.extend(path.windows(2).map(|w| (w[0], w[1])));
    }
    if directed {
        Graph::directed(n, edges)
    } else {
        Graph::undirected(n, edges)
    }
}

/// Run 73, `genspecialg -qz -c12 -b11,4 -G30,-30 -T3,4,5`: the directed versions, including a
/// 30x30 grid that is cyclic in one direction, with countg's `|Aut|` and orbits; all four are
/// rebuilt here as well.
#[test]
fn genspecialg_digraphs() {
    let graphs = read_graphs(&fixture("nautestj.ans"));
    let expected: [(u128, usize); 4] = [(12, 1), (958003200, 2), (30, 30), (1, 11)];
    for (g, &(order, orbits)) in graphs.iter().zip(&expected) {
        assert!(g.directed);
        let symmetry = g.symmetry();
        assert_eq!(
            (symmetry.order, symmetry.node_orbits.len()),
            (order, orbits)
        );
        g.check_invariance(&symmetry.canonical_form, 1, &mut rng());
    }
    let cycle = Graph::directed(12, (0..12).map(|i| (i, (i + 1) % 12)).collect());
    let bipartite = Graph::directed(
        15,
        (0..11)
            .flat_map(|a| (11..15).map(move |b| (a, b)))
            .collect(),
    );
    let grid = Graph::directed(
        900,
        (0..30)
            .flat_map(|r| (0..30).map(move |c| (r, c)))
            .flat_map(|(r, c)| {
                let open = (r < 29).then_some((30 * r + c, 30 * (r + 1) + c));
                open.into_iter()
                    .chain([(30 * r + c, 30 * r + (c + 1) % 30)])
            })
            .collect(),
    );
    for (built, g) in [cycle, bipartite, grid, theta(&[3, 4, 5], true)]
        .iter()
        .zip(&graphs)
    {
        assert_eq!(built.canonical_form(), g.canonical_form());
    }
}

/// Run 74, `gentreeg -q 3:4 | assembleg -n4:12L -cq`: the 19 disjoint unions of the trees on 3
/// and 4 vertices (P3, P4, K(1,3)) with 4 to 12 vertices. Isomorphic components are permuted:
/// `|Aut| = prod |Aut(C)|^m m!`, as countg reports.
#[test]
fn assembleg_disjoint_unions() {
    let unions = read_graphs(&fixture("nautestk.ans"));
    let orders: Vec<u128> = unions.iter().map(|g| g.symmetry().order).collect();
    assert_eq!(
        orders,
        [8, 48, 384, 16, 48, 4, 16, 24, 12, 144, 2, 8, 48, 48, 12, 144, 6, 72, 1296]
    );
    for (g, &order) in unions.iter().zip(&orders) {
        // components by shape: P3 and P4 have 2 automorphisms, K(1,3) has 6
        let mut shapes = BTreeMap::new();
        for component in g.components() {
            let degrees: Vec<usize> = component
                .iter()
                .map(|&v| g.edges.iter().filter(|&&(a, b)| a == v || b == v).count())
                .collect();
            let shape = (component.len(), *degrees.iter().max().unwrap());
            *shapes.entry(shape).or_insert(0u32) += 1;
        }
        let formula: u128 = shapes
            .iter()
            .map(|(&(_, max_degree), &m)| {
                let aut: u128 = if max_degree == 3 { 6 } else { 2 };
                aut.pow(m) * (1..=m as u128).product::<u128>()
            })
            .product();
        assert_eq!(order, formula, "{g:?}");
    }
    let (p3, p4) = (
        Graph::undirected(3, vec![(0, 1), (1, 2)]),
        Graph::undirected(4, vec![(0, 1), (1, 2), (2, 3)]),
    );
    let claw = Graph::undirected(4, vec![(0, 1), (0, 2), (0, 3)]);
    let mut built = BTreeSet::new();
    for a in 0..=4 {
        for b in 0..=3 {
            for c in 0..=3 {
                let n = 3 * a + 4 * (b + c);
                if !(4..=12).contains(&n) {
                    continue;
                }
                let mut union = Graph::undirected(0, vec![]);
                for (part, m) in [(&p3, a), (&p4, b), (&claw, c)] {
                    for _ in 0..m {
                        let offset = union.n;
                        union
                            .edges
                            .extend(part.edges.iter().map(|&(x, y)| (x + offset, y + offset)));
                        union.n += part.n;
                    }
                }
                union.colours = vec![0; union.n];
                built.insert(union.canonical_form());
            }
        }
    }
    let read: BTreeSet<Form> = unions.iter().map(Graph::canonical_form).collect();
    assert_eq!(read.len(), 19);
    assert_eq!(built, read);
}

/// Run 76, `genquarticg -q 9 | copyg -i`: the 16 connected quartic graphs on 9 vertices, read
/// from incremental sparse6, against Linnet's generator (one 4-valent signature, 10 loops,
/// simple graphs only), with countg's group sizes.
#[test]
fn genquarticg_quartic_graphs() {
    let quartic = read_graphs(&fixture("nautestm.ans"));
    assert_eq!(quartic.len(), 16);
    let mut orders: Vec<u128> = quartic.iter().map(|g| g.symmetry().order).collect();
    assert_eq!(
        orders,
        [4, 8, 12, 18, 8, 4, 32, 12, 16, 18, 2, 2, 2, 72, 8, 16]
    );
    let read: BTreeSet<Form> = quartic.iter().map(Graph::canonical_form).collect();
    assert_eq!(read.len(), 16);
    let generator = GraphGenerator::<u8, u8>::new([], [vec![Stub::undirected(0u8); 4]])
        .loops(10..=10)
        .max_vertices(9)
        .filter(|g| {
            // no parallel edges (self-loops are off by default)
            let mut edges: Vec<_> = g
                .edges()
                .map(|e| {
                    (
                        e.vertices[0].min(e.vertices[1]),
                        e.vertices[0].max(e.vertices[1]),
                    )
                })
                .collect();
            edges.sort_unstable();
            edges.windows(2).all(|w| w[0] != w[1])
        });
    let mut generated = BTreeSet::new();
    let mut generated_orders = vec![];
    for diagram in generator.iter().unwrap() {
        let g = Graph::of_hedge_graph(diagram.graph());
        assert_eq!(g.n, 9);
        let symmetry = g.symmetry();
        assert_eq!(diagram.automorphism_count, Some(symmetry.order));
        generated_orders.push(symmetry.order);
        assert!(generated.insert(symmetry.canonical_form));
    }
    assert_eq!(generated, read);
    orders.sort_unstable();
    generated_orders.sort_unstable();
    assert_eq!(generated_orders, orders);
}

/// Run 77, `gengL -d2D2 33`: the 556 2-regular graphs on 33 vertices (beyond one 32-bit word),
/// disjoint unions of cycles, one per partition of 33 into parts of at least 3, with
/// `|Aut| = prod (2k)^m m!` for `m` cycles of length `k`.
#[test]
fn gengl_two_regular_graphs() {
    let graphs = read_graphs(&fixture("nautestn.ans"));
    assert_eq!(graphs.len(), 556);
    let formula = |lengths: &[usize]| -> u128 {
        let mut counts = BTreeMap::new();
        for &k in lengths {
            *counts.entry(k).or_insert(0u32) += 1;
        }
        counts
            .iter()
            .map(|(&k, &m)| (2 * k as u128).pow(m) * (1..=m as u128).product::<u128>())
            .product()
    };
    let mut read = BTreeSet::new();
    for g in &graphs {
        let lengths: Vec<usize> = g.components().iter().map(Vec::len).collect();
        let symmetry = g.symmetry();
        assert_eq!(symmetry.order, formula(&lengths));
        read.insert(symmetry.canonical_form);
    }
    fn partitions(rest: usize, smallest: usize, parts: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if rest == 0 {
            out.push(parts.clone());
        }
        for k in smallest..=rest {
            parts.push(k);
            partitions(rest - k, k, parts, out);
            parts.pop();
        }
    }
    let mut all = vec![];
    partitions(33, 3, &mut vec![], &mut all);
    let built: BTreeSet<Form> = all
        .iter()
        .map(|lengths| {
            let mut edges = vec![];
            let mut start = 0;
            for &k in lengths {
                edges.extend((0..k).map(|i| (start + i, start + (i + 1) % k)));
                start += k;
            }
            Graph::undirected(33, edges).canonical_form()
        })
        .collect();
    assert_eq!(built.len(), 556);
    assert_eq!(built, read);
}

/// Run 80, `genktreeg -k3 -q 11 | complg -q`: complements of the 1505 3-trees on 11 vertices.
/// Complementation preserves (non-)isomorphism, so both the graphs and their complements, the
/// 3-trees with `3n - 6 = 27` edges, are pairwise distinct.
#[test]
fn genktreeg_complements() {
    let graphs = read_graphs(&fixture("nautestq.ans"));
    assert_eq!(graphs.len(), 1505);
    let mut rng = rng();
    let forms: Vec<Form> = graphs
        .iter()
        .map(|g| {
            let form = g.canonical_form();
            g.check_invariance(&form, 1, &mut rng);
            form
        })
        .collect();
    assert_eq!(distinct(forms), 1505);
    let trees: Vec<Graph> = graphs.iter().map(Graph::complement).collect();
    assert!(trees.iter().all(|t| t.edges.len() == 27));
    assert_eq!(distinct(trees.iter().map(Graph::canonical_form)), 1505);
}

/// Run 82, `geng -q -D4 9 | ranlabg -q | shortg -q`: the 15471 graphs on 9 vertices with maximum
/// degree 4 survive random relabelling and isomorph rejection. Here two random relabellings of
/// each are shuffled in, and grouping by canonical form recovers each class with its copies.
#[test]
fn shortg_random_relabellings() {
    let graphs = read_graphs(&fixture("nautests.ans"));
    assert_eq!(graphs.len(), 15471);
    let mut rng = rng();
    let mut copies: Vec<(usize, Graph)> = graphs
        .iter()
        .enumerate()
        .flat_map(|(i, g)| {
            [
                (i, g.clone()),
                (i, g.shuffled(&mut rng)),
                (i, g.shuffled(&mut rng)),
            ]
        })
        .collect();
    copies.shuffle(&mut rng);
    let mut classes: HashMap<Form, Vec<usize>> = HashMap::new();
    for (i, g) in &copies {
        classes.entry(g.canonical_form()).or_default().push(*i);
    }
    assert_eq!(classes.len(), 15471);
    assert!(classes
        .values()
        .all(|c| c.len() == 3 && c.iter().all(|&i| i == c[0])));
}

/// Run 95, `geng -Cq 8 12 | watercluster2 S Z | underlyingg -q | countg --ee,cc,zZ`: the
/// 1158026 orientations (no digons) of the 429 biconnected graphs with 8 vertices and 12 edges,
/// tabulated by the radius and diameter of the underlying graph. Reduced runs check the bases
/// of radius 3 only.
#[test]
fn watercluster_orientations() {
    let answer = fixture("nautestff.ans");
    assert!(answer.contains("Number of directed graphs: 1158026"));
    let mut expected = BTreeMap::new();
    for line in answer.lines() {
        if let Some((count, keys)) = line.split_once(" graphs : nonedges=16; components=1; radius=")
        {
            let (radius, diameter) = keys.split_once("; diameter=").unwrap();
            let key = (
                radius.parse::<usize>().unwrap(),
                diameter.parse::<usize>().unwrap(),
            );
            expected.insert(key, count.trim().parse::<usize>().unwrap());
        }
    }
    assert_eq!(expected.values().sum::<usize>(), 1158026);
    let bases = read_graphs(&fixture("watercluster2-bases.g6"));
    assert_eq!(bases.len(), 429);
    let full = full_size();
    let mut table = BTreeMap::new();
    for base in &bases {
        let key = base.radius_diameter();
        if !full && key.0 != 3 {
            continue;
        }
        let order = base.symmetry().order;
        *table.entry(key).or_insert(0) += classify(orientations(base), order, same_order);
    }
    if !full {
        expected.retain(|k, _| k.0 == 3);
    }
    assert_eq!(table, expected);
}

/// Run 96, `geng -qd3D4 8 | multig -qTm2`: multigraphs with edge multiplicities 1 or 2 on the 118
/// graphs with 8 vertices and degrees 3 to 4, iso-rejected per base (947540 in all, counted per
/// base on the oracle). Linnet also counts the swap of the two edges of each double edge, so a
/// multigraph's stabiliser in the base group is `|Aut| / 2^doubles`. Reduced runs check the cubic
/// bases only.
#[test]
fn multig_multigraphs() {
    let full = full_size();
    let mut total = 0;
    for line in fixture("multig-bases.txt").lines() {
        let (g6, count) = line.split_once(' ').unwrap();
        let count: usize = count.parse().unwrap();
        total += count;
        let base = read_graphs(g6).pop().unwrap();
        if !full && base.edges.len() > 12 {
            continue;
        }
        let order = base.symmetry().order;
        let multigraphs = (0..1u64 << base.edges.len()).map(|mask| {
            let mut edges = base.edges.clone();
            edges.extend(
                base.edges
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask >> i & 1 == 1)
                    .map(|(_, &e)| e),
            );
            Graph::undirected(base.n, edges)
        });
        let stabiliser = |g: &Graph, order: u128| order >> (g.edges.len() - base.edges.len());
        assert_eq!(classify(multigraphs, order, stabiliser), count, "{g6}");
    }
    assert_eq!(total, 947540);
}

/// Run 98, `genspecialg -q -X12 -X44 -X23 -X68 -X100 | distgraphg -qd3 | uniqg -uy`: the graphs
/// joining the vertices at distance 3 in five named graphs are pairwise non-isomorphic, with
/// countg's `|Aut|` and orbits for both the named graphs and their distance-3 graphs.
#[test]
fn distgraphg_distance_three_graphs() {
    let answer = fixture("nautestii.ans");
    assert!(answer.contains("5 graphs read from stdin, 5 unique"));
    let named = read_graphs(&fixture("distgraphg-bases.s6"));
    let expected: [(u128, usize); 5] = [(2, 30), (54, 1), (4, 16), (120, 4), (2, 10)];
    let mut forms = BTreeSet::new();
    for (g, &(order, orbits)) in named.iter().zip(&expected) {
        let distances = g.distances();
        let d3 = Graph::undirected(
            g.n,
            (0..g.n)
                .flat_map(|b| (0..b).map(move |a| (a, b)))
                .filter(|&(a, b)| distances[a][b] == 3)
                .collect(),
        );
        for h in [g, &d3] {
            let symmetry = h.symmetry();
            assert_eq!(
                (symmetry.order, symmetry.node_orbits.len()),
                (order, orbits)
            );
        }
        forms.insert(d3.canonical_form());
    }
    assert_eq!(forms.len(), 5);
}
