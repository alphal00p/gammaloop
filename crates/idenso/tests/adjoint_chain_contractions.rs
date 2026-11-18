//! Adjoint generator lines with every colour index summed except the two
//! line endpoints reduce to a group invariant times the adjoint metric.
//!
//! The local UV parity tests in GammaLoop (`factorized_growth_matches_*`)
//! compare two forms of the same integrated poles; which form a run produces
//! depends on the symbol registration order of the process. These are the
//! colour structures of the form that failed to reduce.

use idenso::color::ColorSimplifySettings;
use idenso::dirac::GammaSimplifySettings;
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use symbolica::atom::{Atom, AtomCore};

fn initialize() {
    idenso::representations::initialize();
    let _ = *idenso::epsilon::EPSILON_SYMBOL;
}

fn parse(source: &str) -> Atom {
    Atom::parse(
        source,
        "spenso",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap()
}

fn simplify(expression: &Atom) -> Atom {
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view()).unwrap())
        .unwrap()
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            color: Some(ColorSimplifySettings::default()),
            epsilon: true,
            ..Default::default()
        })
        .unwrap()
        .to_dots()
        .unwrap()
        .expanded(None, false)
        .unwrap()
        .canonize()
        .unwrap();
    cooking.uncook(tensor.expression().as_view()).expand()
}

fn generator(label: &str) -> String {
    format!("f(coad(8,{label}),in,out)")
}

fn line(start: &str, end: &str, labels: &[&str]) -> Atom {
    let generators = labels.iter().map(|l| generator(l)).collect::<Vec<_>>();
    parse(&format!(
        "chain(coad(8,{start}),coad(8,{end}),{})",
        generators.join(",")
    ))
}

/// The structure constants of `line`, written out: `F^x_{ij} = f(x, i, j)`.
fn structure_constants(start: &str, end: &str, labels: &[&str], scope: &str) -> Atom {
    let mut product = Atom::one();
    let mut previous = start.to_string();
    for (k, label) in labels.iter().enumerate() {
        let next = if k + 1 == labels.len() {
            end.to_string()
        } else {
            format!("{scope}_{k}")
        };
        product *= parse(&format!(
            "f(coad(8,{label}),coad(8,{previous}),coad(8,{next}))"
        ));
        previous = next;
    }
    product
}

fn reduced(expression: &Atom) -> bool {
    let text = expression.to_string();
    !text.contains("f(") && !text.contains("chain(") && !text.contains("trace(")
}

/// All orderings of the multiset `labels`.
fn arrangements<'a>(labels: &[&'a str]) -> Vec<Vec<&'a str>> {
    let mut out = std::collections::BTreeSet::new();
    let n = labels.len();
    let mut indices: Vec<usize> = (0..n).collect();
    // Heap's algorithm over positions.
    fn heap<'a>(
        k: usize,
        indices: &mut Vec<usize>,
        labels: &[&'a str],
        out: &mut std::collections::BTreeSet<Vec<&'a str>>,
    ) {
        if k == 1 {
            out.insert(indices.iter().map(|&i| labels[i]).collect());
            return;
        }
        heap(k - 1, indices, labels, out);
        for i in 0..k - 1 {
            if k % 2 == 0 {
                indices.swap(i, k - 1);
            } else {
                indices.swap(0, k - 1);
            }
            heap(k - 1, indices, labels, out);
        }
    }
    heap(n, &mut indices, labels, &mut out);
    out.into_iter().collect()
}

#[test]
fn four_generator_adjoint_lines_with_summed_labels_reduce_to_invariants() {
    initialize();
    let mut failures = Vec::new();
    // Label names chosen so that the summed labels sort before, between and
    // after the endpoints.
    for (a, b, x, y) in [
        ("a", "b", "x", "y"),
        ("x", "y", "a", "b"),
        ("a", "y", "b", "x"),
        ("m", "n", "c", "z"),
    ] {
        let mut cases = Vec::new();
        // Both summed labels on the line.
        for labels in arrangements(&[x, x, y, y]) {
            cases.push((a, b, labels));
        }
        // One summed label on the line and on an endpoint.
        for labels in arrangements(&[x, x, y, b]) {
            cases.push((a, y, labels.clone()));
            cases.push((y, a, labels));
        }
        for (start, end, labels) in cases {
            let chain = line(start, end, &labels);
            let explicit = structure_constants(start, end, &labels, "internal");
            let from_chain = simplify(&chain);
            let from_explicit = simplify(&explicit);
            if !reduced(&from_chain)
                || !reduced(&from_explicit)
                || !(&from_chain - &from_explicit).expand().is_zero()
            {
                failures.push(format!(
                    "{chain}\n  line -> {from_chain}\n  structure constants -> {from_explicit}"
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} unreduced or inconsistent lines:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn crossed_adjoint_line_pairs_reduce_to_the_metric() {
    initialize();
    // Two 2-generator lines joined crosswise: each line's endpoints label
    // the other's generators. Found in the integrated poles of
    // factorized_growth_matches_quartic_gluon (round 2 of the local_4d flake).
    let cases = [
        (
            vec![("a", "b", vec!["q", "c"]), ("c", "p", vec!["a", "b"])],
            true,
        ),
        (
            vec![("a", "p", vec!["q", "e"]), ("c", "e", vec!["a", "c"])],
            true,
        ),
    ];
    let mut failures = Vec::new();
    for (lines, two_ports) in cases {
        let network = Network::from_named(&lines, two_ports);
        if let Some(failure) = network.check(&NAMES[0]) {
            failures.push(failure);
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The names of a network's labels: summed labels and the two free ports.
struct Names {
    summed: [&'static str; 6],
    ports: [&'static str; 2],
    internal: &'static str,
}

/// Summed labels sorting before the ports, and after them.
const NAMES: [Names; 2] = [
    Names {
        summed: ["a", "b", "c", "d", "e", "h"],
        ports: ["q", "p"],
        internal: "z",
    },
    Names {
        summed: ["r", "s", "t", "u", "v", "w"],
        ports: ["b", "a"],
        internal: "k",
    },
];

const PORT: [usize; 2] = [usize::MAX - 1, usize::MAX];

/// Adjoint lines whose slots (start, end, generator labels) carry summed
/// label ids, each used twice, or one of the two free ports.
#[derive(Clone, Debug)]
struct Network {
    lengths: Vec<usize>,
    slots: Vec<usize>,
}

impl Network {
    fn from_named(lines: &[(&str, &str, Vec<&str>)], two_ports: bool) -> Self {
        let mut ids: Vec<String> = Vec::new();
        let mut lengths = Vec::new();
        let mut slots = Vec::new();
        let ports = ["q", "p"];
        for (start, end, labels) in lines {
            lengths.push(labels.len());
            for name in [*start, *end].into_iter().chain(labels.iter().copied()) {
                if two_ports && let Some(port) = ports.iter().position(|p| *p == name) {
                    slots.push(PORT[port]);
                    continue;
                }
                let id = ids.iter().position(|seen| seen == name).unwrap_or_else(|| {
                    ids.push(name.to_string());
                    ids.len() - 1
                });
                slots.push(id);
            }
        }
        Self { lengths, slots }
    }

    fn name(&self, id: usize, names: &Names) -> &'static str {
        match PORT.iter().position(|port| *port == id) {
            Some(port) => names.ports[port],
            None => names.summed[id],
        }
    }

    fn two_ports(&self) -> bool {
        self.slots.contains(&PORT[0])
    }

    /// Each line as (start, end, generator labels).
    fn lines(&self, names: &Names) -> Vec<(&'static str, &'static str, Vec<&'static str>)> {
        let mut offset = 0;
        self.lengths
            .iter()
            .map(|&length| {
                let slots = &self.slots[offset..offset + length + 2];
                offset += length + 2;
                (
                    self.name(slots[0], names),
                    self.name(slots[1], names),
                    slots[2..].iter().map(|&id| self.name(id, names)).collect(),
                )
            })
            .collect()
    }

    fn chains(&self, names: &Names) -> Atom {
        self.lines(names)
            .iter()
            .fold(Atom::one(), |product, (start, end, labels)| {
                product * line(start, end, labels)
            })
    }

    fn structure_constants(&self, names: &Names) -> Atom {
        self.lines(names).iter().enumerate().fold(
            Atom::one(),
            |product, (k, (start, end, labels))| {
                product * structure_constants(start, end, labels, &format!("{}{k}", names.internal))
            },
        )
    }

    /// `None` when idenso reduces the lines to group invariants (times the
    /// metric of the two ports) equal to the closed-network reduction of the
    /// written-out structure constants (closed with the ports' metric).
    fn check(&self, names: &Names) -> Option<String> {
        let chains = self.chains(names);
        let from_chains = simplify(&chains);
        let metric = parse(&format!(
            "g(coad(8,{}),coad(8,{}))",
            names.ports[0], names.ports[1]
        ));
        let mut explicit = self.structure_constants(names);
        let mut projected = from_chains.clone();
        if self.two_ports() {
            explicit = explicit * &metric;
            // c g(q,p) closed with g(q,p) is c N_A.
            let reversed = parse(&format!(
                "g(coad(8,{}),coad(8,{}))",
                names.ports[1], names.ports[0]
            ));
            projected = projected
                .replace(metric.to_pattern())
                .with(Atom::num(8).to_pattern())
                .replace(reversed.to_pattern())
                .with(Atom::num(8).to_pattern())
                .expand();
        }
        let reference = idenso::color::network::try_reduce_closed_color_network(explicit.as_view());
        let reference = match reference {
            Ok(Some(reference)) => reference,
            other => return Some(format!("{chains}\n  closed-network reducer: {other:?}")),
        };
        let leftover = |expression: &Atom| {
            let text = expression.to_string();
            names
                .ports
                .iter()
                .any(|port| text.contains(&format!("coad(8,{port})")))
        };
        if !reduced(&from_chains)
            || leftover(&projected)
            || !(&projected - &reference).expand().is_zero()
        {
            return Some(format!(
                "{chains}\n  lines -> {from_chains}\n  closed-network reducer -> {reference}"
            ));
        }
        None
    }
}

/// Every connected network of adjoint lines of the given generator counts
/// with all labels summed (no free port) or two free ports, up to
/// permutations of equally long lines and the exchange of the ports.
fn networks(lengths: &[usize]) -> Vec<Network> {
    let size: usize = lengths.iter().map(|length| length + 2).sum();
    let mut found = std::collections::BTreeMap::new();
    let mut port_choices: Vec<Option<(usize, usize)>> = vec![None];
    for i in 0..size {
        for j in i + 1..size {
            port_choices.push(Some((i, j)));
        }
    }
    let line_of: Vec<usize> = lengths
        .iter()
        .enumerate()
        .flat_map(|(k, &length)| std::iter::repeat_n(k, length + 2))
        .collect();
    let symmetries = line_permutations(lengths);
    for ports in port_choices {
        let mut slots = vec![usize::MAX - 2; size];
        if let Some((i, j)) = ports {
            slots[i] = PORT[0];
            slots[j] = PORT[1];
        }
        let mut matchings = Vec::new();
        pair_up(&mut slots, 0, &mut matchings);
        for slots in matchings {
            if !connected(&slots, &line_of, lengths.len()) {
                continue;
            }
            let network = Network {
                lengths: lengths.to_vec(),
                slots,
            };
            let key = canonical_key(&network, &symmetries);
            found.entry(key).or_insert(network);
        }
    }
    found.into_values().collect()
}

fn pair_up(slots: &mut Vec<usize>, next: usize, out: &mut Vec<Vec<usize>>) {
    let free = usize::MAX - 2;
    let Some(first) = slots.iter().position(|&slot| slot == free) else {
        out.push(slots.clone());
        return;
    };
    for partner in first + 1..slots.len() {
        if slots[partner] != free {
            continue;
        }
        slots[first] = next;
        slots[partner] = next;
        pair_up(slots, next + 1, out);
        slots[first] = free;
        slots[partner] = free;
    }
}

fn connected(slots: &[usize], line_of: &[usize], lines: usize) -> bool {
    let mut parent: Vec<usize> = (0..lines).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] == x {
            x
        } else {
            let root = find(parent, parent[x]);
            parent[x] = root;
            root
        }
    }
    for (i, &label) in slots.iter().enumerate() {
        if PORT.contains(&label) {
            continue;
        }
        for (j, &other) in slots.iter().enumerate().skip(i + 1) {
            if other == label {
                let (a, b) = (find(&mut parent, line_of[i]), find(&mut parent, line_of[j]));
                parent[a] = b;
            }
        }
    }
    (0..lines).all(|line| find(&mut parent, line) == find(&mut parent, 0))
}

/// The permutations of the lines that keep their lengths.
fn line_permutations(lengths: &[usize]) -> Vec<Vec<usize>> {
    fn permute(rest: &mut Vec<usize>, prefix: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if rest.is_empty() {
            out.push(prefix.clone());
            return;
        }
        for k in 0..rest.len() {
            let item = rest.remove(k);
            prefix.push(item);
            permute(rest, prefix, out);
            prefix.pop();
            rest.insert(k, item);
        }
    }
    let mut all = Vec::new();
    permute(&mut (0..lengths.len()).collect(), &mut Vec::new(), &mut all);
    all.retain(|order| {
        order
            .iter()
            .enumerate()
            .all(|(k, &line)| lengths[k] == lengths[line])
    });
    all
}

fn canonical_key(network: &Network, symmetries: &[Vec<usize>]) -> Vec<usize> {
    let mut offsets = vec![0];
    for length in &network.lengths {
        offsets.push(offsets.last().unwrap() + length + 2);
    }
    let mut best: Option<Vec<usize>> = None;
    for order in symmetries {
        for swap in [false, true] {
            let mut renamed = std::collections::HashMap::new();
            let key: Vec<usize> = order
                .iter()
                .flat_map(|&line| network.slots[offsets[line]..offsets[line + 1]].iter())
                .map(|&label| match PORT.iter().position(|port| *port == label) {
                    Some(port) => PORT[port ^ usize::from(swap)],
                    None => {
                        let next = renamed.len();
                        *renamed.entry(label).or_insert(next)
                    }
                })
                .collect();
            if best.as_ref().is_none_or(|best| key < *best) {
                best = Some(key);
            }
        }
    }
    best.unwrap()
}

/// Check every network of the given line lengths under both label orders;
/// writes the failures to `$ADJOINT_NETWORK_REPORT/<lengths>.txt` if set.
fn check_networks(lengths: &[usize]) {
    initialize();
    let all = networks(lengths);
    let mut failures = Vec::new();
    let (mut closed, mut open) = (0, 0);
    for network in &all {
        if network.two_ports() {
            open += 1
        } else {
            closed += 1
        }
        for names in &NAMES {
            if let Some(failure) = network.check(names) {
                failures.push(failure);
            }
        }
    }
    let tag = lengths
        .iter()
        .map(|l| l.to_string())
        .collect::<Vec<_>>()
        .join("-");
    let summary = format!(
        "lines {tag}: {} networks ({closed} closed, {open} with two ports), {} checks, {} failed",
        all.len(),
        2 * all.len(),
        failures.len()
    );
    println!("ADJOINT-NETWORKS {summary}");
    if let Some(dir) = std::env::var_os("ADJOINT_NETWORK_REPORT") {
        let path = std::path::Path::new(&dir).join(format!("lines-{tag}.txt"));
        std::fs::write(path, format!("{summary}\n{}\n", failures.join("\n"))).unwrap();
    }
    assert!(
        failures.is_empty(),
        "{summary}\n{}",
        failures
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn every_network_of_one_line_of_two_generators_reduces() {
    check_networks(&[2]);
}

#[test]
fn every_network_of_one_line_of_four_generators_reduces() {
    check_networks(&[4]);
}

#[test]
fn every_network_of_two_lines_of_two_generators_reduces() {
    check_networks(&[2, 2]);
}

#[test]
fn every_network_of_lines_of_two_and_four_generators_reduces() {
    check_networks(&[2, 4]);
}

#[test]
fn every_network_of_two_lines_of_three_generators_reduces() {
    check_networks(&[3, 3]);
}

#[test]
fn every_network_of_three_lines_of_two_generators_reduces() {
    check_networks(&[2, 2, 2]);
}
