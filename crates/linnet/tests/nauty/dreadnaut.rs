//! The `dreadtest*` runs: nauty's `nautest{1,2,3}.dre` dreadnaut scripts executed on Linnet,
//! checked against the algorithm-independent records of their `.ans` transcripts.
//!
//! Only the commands the three scripts use are interpreted; anything else panics. Each script is
//! checked against the transcript of every run that uses it (all word sizes, `MAXN` settings and
//! modes share three transcripts per script). Group sizes, orbit counts and orbits, degrees, the
//! graph header and whether two canonical-graph hashes agree must match exactly; nauty's printed
//! generators must be automorphisms of the graph and generate a group of the printed size.
//! Search statistics, canonical labellings and hash values are nauty's own and are not compared.
//!
//! `j` relabels with nauty's random number generator (KISS seeded by `ran_init(1)`), reproduced
//! here, so labelled output after `j` is comparable. Traces draws from the same generator while
//! searching, so after a `j` in Traces mode only label-free data (group size, orbit sizes) is
//! compared. Traces refuses digraphs, so the third script's Traces transcript only shows the refusal.

use std::collections::BTreeSet;

use rand::{rngs::SmallRng, SeedableRng};

use super::{
    formats::parse_graph,
    graph::{group_order, Cert, Graph, Symmetry},
    runs::{fixture, Run},
};

/// nauty's mode, from `-o As` (sparse) and `-o At` (Traces) on the `dreadtest` command line.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Mode {
    Dense,
    Sparse,
    Traces,
}

impl Mode {
    fn of_command(command: &str) -> Mode {
        if command.contains("-o At") {
            Mode::Traces
        } else if command.contains("-o As") {
            Mode::Sparse
        } else {
            Mode::Dense
        }
    }

    fn of_header(name: &str) -> Mode {
        match name {
            "dense" => Mode::Dense,
            "sparse" => Mode::Sparse,
            "Traces" => Mode::Traces,
            _ => panic!("unknown mode {name}"),
        }
    }
}

/// nauty's 64-bit KISS generator (`naurng.c`, by George Marsaglia).
struct Kiss {
    x: u64,
    y: u64,
    z: u64,
    c: u64,
}

impl Kiss {
    /// `ran_init(seed)`.
    fn new(seed: u64) -> Self {
        let mut kiss = Kiss {
            x: 1234567890987654321u64.wrapping_add(seed),
            y: 362436362436362436,
            z: 1066149217761810,
            c: 123456123456123456,
        };
        for _ in 0..1000 {
            kiss.next();
        }
        kiss
    }

    fn next(&mut self) -> u64 {
        let t = (self.x << 58).wrapping_add(self.c);
        self.c = self.x >> 6;
        self.x = self.x.wrapping_add(t);
        self.c = self.c.wrapping_add((self.x < t) as u64);
        self.y ^= self.y << 13;
        self.y ^= self.y >> 17;
        self.y ^= self.y << 43;
        self.z = self.z.wrapping_mul(6906969069).wrapping_add(1234567);
        self.x.wrapping_add(self.y).wrapping_add(self.z)
    }

    /// `ranperm(perm, n)`.
    fn permutation(&mut self, n: usize) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            perm.swap(i, j);
        }
        perm
    }
}

/// What a script makes Linnet observe, in output order.
#[derive(Debug)]
enum Observed {
    /// `?`: the graph header.
    Header {
        n: usize,
        edges: usize,
        directed: bool,
    },
    /// `v`: the (out-)degree of every vertex.
    Degrees {
        degrees: Vec<usize>,
        labelorg: usize,
    },
    /// `x` with a partition set by `f`.
    Fixing,
    /// `x`: the automorphism group of `graph`, coloured by the partition.
    Group {
        symmetry: Symmetry,
        graph: Graph,
        labelorg: usize,
        /// Whether `graph` is labelled as nauty's graph is.
        labelled: bool,
    },
    /// `o`: the orbits of the last `x`.
    Orbits {
        orbits: Vec<Vec<usize>>,
        labelorg: usize,
        labelled: bool,
    },
    /// `z`: the canonical graph of the last `x`, here its certificate.
    Hash(Cert),
}

/// The algorithm-independent records of a transcript, labels as printed.
#[derive(Debug)]
enum Expected {
    Header {
        mode: Mode,
        n: usize,
        edges: usize,
        directed: bool,
    },
    /// Vertex label ranges and their degree.
    Degrees(Vec<(usize, usize, usize)>),
    Fixing,
    Group {
        order: u128,
        orbits: usize,
        /// Printed generators in cycle notation, empty when `-a` suppresses them.
        generators: Vec<Vec<Vec<usize>>>,
    },
    Orbits(Vec<Vec<usize>>),
    Hash(String),
    /// "Traces can't handle digraphs".
    Refused,
}

/// A dreadnaut session on Linnet.
struct Session {
    mode: Mode,
    labelorg: usize,
    old_labelorg: usize,
    graph: Option<Graph>,
    /// `f=v`: the partition `[v | rest]`.
    fixed: Option<usize>,
    rng: Kiss,
    labelled: bool,
    last: Option<Symmetry>,
    observed: Vec<Observed>,
    /// For extra relabellings at every `x`, independent of nauty's generator.
    check_rng: SmallRng,
}

/// A script being read, with the character-level helpers of dreadnaut's input routines.
struct Input<'a> {
    text: &'a [u8],
    at: usize,
}

impl Input<'_> {
    fn next(&mut self) -> Option<u8> {
        let c = self.text.get(self.at).copied();
        self.at += 1;
        c
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.at).copied()
    }

    /// Skip spaces and tabs, not newlines (`GETNW`).
    fn skip_blanks(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t')) {
            self.at += 1;
        }
    }

    fn integer(&mut self) -> Option<usize> {
        let start = self.at;
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.at += 1;
        }
        std::str::from_utf8(&self.text[start..self.at])
            .unwrap()
            .parse()
            .ok()
    }

    /// `getint_sl`: an integer, optionally preceded by `=` and blanks.
    fn assigned_integer(&mut self) -> usize {
        self.skip_blanks();
        if self.peek() == Some(b'=') {
            self.at += 1;
        }
        self.skip_blanks();
        self.integer().expect("an integer argument")
    }
}

impl Session {
    fn new(mode: Mode) -> Self {
        Session {
            mode,
            labelorg: 0,
            old_labelorg: 0,
            graph: None,
            fixed: None,
            rng: Kiss::new(1),
            labelled: true,
            last: None,
            observed: vec![],
            check_rng: SmallRng::seed_from_u64(1),
        }
    }

    fn graph(&self) -> &Graph {
        self.graph.as_ref().expect("g is defined")
    }

    fn run(mut self, script: &str) -> Vec<Observed> {
        let mut input = Input {
            text: script.as_bytes(),
            at: 0,
        };
        let mut n = 0;
        while let Some(c) = input.next() {
            match c {
                b' ' | b'\t' | b'\n' | b'\r' => {}
                b'!' => {
                    while input.peek().is_some_and(|c| c != b'\n') {
                        input.at += 1;
                    }
                }
                b'n' => {
                    n = input.assigned_integer();
                    self.graph = None;
                    self.fixed = None;
                }
                b'$' => {
                    if input.peek() == Some(b'$') {
                        input.at += 1;
                        self.labelorg = self.old_labelorg;
                    } else {
                        self.old_labelorg = self.labelorg;
                        self.labelorg = input.assigned_integer();
                    }
                }
                b'g' => {
                    self.graph = Some(self.read_graph(&mut input, n));
                    self.fixed = None;
                }
                b'(' => {
                    let mut code = vec![];
                    loop {
                        match input.next().expect("a closing parenthesis") {
                            b')' => break,
                            c if c.is_ascii_whitespace() => {}
                            c => code.push(c),
                        }
                    }
                    let graph = parse_graph(std::str::from_utf8(&code).unwrap(), None);
                    n = graph.n;
                    self.graph = Some(graph);
                    self.fixed = None;
                }
                b'?' => {
                    let g = self.graph();
                    let edges = g.edges.len();
                    let (n, directed) = (g.n, g.directed);
                    self.observed.push(Observed::Header { n, edges, directed });
                }
                b'v' => {
                    let g = self.graph();
                    let mut degrees = vec![0; g.n];
                    for &(a, b) in &g.edges {
                        degrees[a] += 1;
                        if !g.directed {
                            degrees[b] += 1;
                        }
                    }
                    let labelorg = self.labelorg;
                    self.observed.push(Observed::Degrees { degrees, labelorg });
                }
                b'f' => {
                    let v = input.assigned_integer() - self.labelorg;
                    assert!(v < n, "f fixes a vertex");
                    self.fixed = Some(v);
                }
                b'x' => self.execute(),
                b'o' => {
                    let last = self.last.as_ref().expect("x before o");
                    let orbits = last.node_orbits.clone();
                    self.observed.push(Observed::Orbits {
                        orbits,
                        labelorg: self.labelorg,
                        labelled: self.labelled,
                    });
                }
                b'z' => {
                    let last = self.last.as_ref().expect("x before z");
                    self.observed.push(Observed::Hash(last.certificate.clone()));
                }
                b'j' => {
                    let perm = self.rng.permutation(n);
                    self.graph = Some(self.graph().relabelled_by(&perm));
                    self.fixed = self
                        .fixed
                        .map(|v| perm.iter().position(|&p| p == v).unwrap());
                    if self.mode == Mode::Traces {
                        self.labelled = false;
                    }
                }
                b'_' => {
                    let g = self.graph();
                    let g = if input.peek() == Some(b'_') {
                        input.at += 1;
                        g.converse()
                    } else {
                        g.complement()
                    };
                    self.graph = Some(g);
                }
                // `&` and `&&` print the partition, `i` refines it, `t` types the graph, `b` the
                // canonical labelling; `c`, `a` and `m` (possibly after `-`) only switch printing
                b'&' | b'i' | b't' | b'b' | b'c' | b'a' | b'm' | b'-' => {}
                b'q' => break,
                c => panic!(
                    "dreadnaut command {:?} is outside the translated subset",
                    c as char
                ),
            }
        }
        self.observed
    }

    /// `g`: adjacency lists, `v :` selecting a vertex, `;` advancing to the next one, ending
    /// after the last vertex or at `.` (dreadnaut's `readgraph`).
    fn read_graph(&self, input: &mut Input, n: usize) -> Graph {
        let mut edges = BTreeSet::new();
        let mut v = 0;
        loop {
            input.skip_blanks();
            match input.next().expect("the graph ends") {
                c if c.is_ascii_digit() => {
                    input.at -= 1;
                    let w = input.integer().unwrap() - self.labelorg;
                    input.skip_blanks();
                    if input.peek() == Some(b':') {
                        input.at += 1;
                        v = w;
                    } else {
                        edges.insert((v.min(w), v.max(w)));
                    }
                }
                b';' => {
                    v += 1;
                    if v == n {
                        break;
                    }
                }
                b'.' => break,
                b'\n' | b'\r' => {}
                c => panic!(
                    "graph input {:?} is outside the translated subset",
                    c as char
                ),
            }
        }
        Graph::undirected(n, edges.into_iter().collect())
    }

    /// `x`: canonize the graph coloured by the partition.
    fn execute(&mut self) {
        let mut graph = self.graph().clone();
        if let Some(v) = self.fixed {
            graph.colours = (0..graph.n).map(|w| (w != v) as u32).collect();
            if self.mode != Mode::Traces {
                self.observed.push(Observed::Fixing);
            }
        }
        let symmetry = graph.symmetry();
        graph.check_invariance(&symmetry.certificate, 2, &mut self.check_rng);
        self.last = Some(symmetry.clone());
        self.observed.push(Observed::Group {
            symmetry,
            graph,
            labelorg: self.labelorg,
            labelled: self.labelled,
        });
    }
}

/// The records of a transcript; see [`Expected`].
fn parse_answers(text: &str) -> Vec<Expected> {
    let lines: Vec<&str> = text.lines().collect();
    let mut expected = vec![];
    let mut generators = vec![];
    let mut i = 0;
    let labels = |s: &str| -> Vec<usize> {
        s.split_whitespace()
            .map(|x| x.parse().expect("a label"))
            .collect()
    };
    while i < lines.len() {
        let line = lines[i].trim();
        i += 1;
        if let Some(header) = line.strip_prefix("Mode=") {
            let field = |name: &str| {
                header
                    .split_whitespace()
                    .find_map(|f| f.strip_prefix(name))
                    .map(|v| v.parse::<usize>().unwrap())
            };
            let arcs = field("arcs=");
            expected.push(Expected::Header {
                mode: Mode::of_header(header.split_whitespace().next().unwrap()),
                n: field("n=").unwrap(),
                edges: arcs.or_else(|| field("edges=")).unwrap(),
                directed: arcs.is_some(),
            });
        } else if line == "[fixing partition]" {
            expected.push(Expected::Fixing);
        } else if line == "Traces can't handle digraphs" {
            expected.push(Expected::Refused);
        } else if line.starts_with('[') && line.ends_with(']') {
            expected.push(Expected::Hash(line.to_string()));
        } else if line.starts_with('(') || line.starts_with("Gen") {
            // a generator in cycle notation, continued on indented lines
            let mut text = line.split_once(": ").map_or(line, |(_, g)| g).to_string();
            while i < lines.len()
                && lines[i].starts_with(' ')
                && !lines[i].contains([':', ';'])
                && lines[i]
                    .trim()
                    .starts_with(|c: char| c == '(' || c.is_ascii_digit())
            {
                text.push(' ');
                text.push_str(lines[i].trim());
                i += 1;
            }
            let cycles = text
                .split(['(', ')'])
                .filter(|c| !c.trim().is_empty())
                .map(labels)
                .collect();
            generators.push(cycles);
        } else if let Some((orbits, rest)) = line
            .split_once(" orbit")
            .filter(|_| line.contains("; grpsize="))
        {
            // the summary of `x`; per-level search statistics are skipped
            let order = rest
                .split("grpsize=")
                .nth(1)
                .and_then(|r| r.split(';').next())
                .unwrap();
            expected.push(Expected::Group {
                order: order.parse().expect("an integral group size"),
                orbits: orbits.parse().unwrap(),
                generators: std::mem::take(&mut generators),
            });
            // `o` prints the orbits right after the statistics of `x`
            while i < lines.len()
                && ["cpu time", "canupdates", "invarproc"]
                    .iter()
                    .any(|p| lines[i].starts_with(p))
            {
                i += 1;
            }
            if i < lines.len() && is_orbit_list(lines[i]) {
                let mut text = lines[i].trim().to_string();
                i += 1;
                while !text.ends_with(';') {
                    text.push(' ');
                    text.push_str(lines[i].trim());
                    i += 1;
                }
                expected.push(Expected::Orbits(parse_orbits(&text)));
            }
        } else if !line.is_empty()
            && line.split_whitespace().all(|t| {
                t.split_once(':').is_some_and(|(r, d)| {
                    d.parse::<usize>().is_ok() && r.split('-').all(|v| v.parse::<usize>().is_ok())
                })
            })
        {
            for token in line.split_whitespace() {
                let (range, degree) = token.split_once(':').unwrap();
                let (from, to) = range.split_once('-').unwrap_or((range, range));
                let degree = degree.parse().unwrap();
                let range = (from.parse().unwrap(), to.parse().unwrap(), degree);
                match expected.last_mut() {
                    Some(Expected::Degrees(ranges)) => ranges.push(range),
                    _ => expected.push(Expected::Degrees(vec![range])),
                }
            }
        }
    }
    expected
}

/// Orbits look like `0:2 4 (4); 3;`, unlike adjacency lists (` 0 :  2 3;`) and the summary
/// of the next `x` (`1 orbit; grpsize=24; ...`).
fn is_orbit_list(line: &str) -> bool {
    line.trim_start().starts_with(|c: char| c.is_ascii_digit())
        && !line.contains(" : ")
        && !line.contains('=')
}

/// Cells separated by `;`, `a:b` ranges and a `(size)` after every non-trivial cell.
fn parse_orbits(text: &str) -> Vec<Vec<usize>> {
    let mut orbits: Vec<Vec<usize>> = text
        .split(';')
        .filter(|c| !c.trim().is_empty())
        .map(|cell| {
            let (members, size) = match cell.split_once('(') {
                Some((m, s)) => (m, s.trim_end_matches(')').trim().parse().unwrap()),
                None => (cell, 1),
            };
            let mut orbit = vec![];
            for token in members.split_whitespace() {
                let (a, b) = token.split_once(':').unwrap_or((token, token));
                orbit.extend(a.parse::<usize>().unwrap()..=b.parse().unwrap());
            }
            assert_eq!(orbit.len(), size, "orbit size in {cell:?}");
            orbit.sort_unstable();
            orbit
        })
        .collect();
    orbits.sort_unstable();
    orbits
}

/// nauty's labels to vertices.
fn unlabel(orbits: &[Vec<usize>], labelorg: usize) -> Vec<Vec<usize>> {
    orbits
        .iter()
        .map(|o| o.iter().map(|v| v - labelorg).collect())
        .collect()
}

fn sizes(orbits: &[Vec<usize>]) -> Vec<usize> {
    let mut sizes: Vec<usize> = orbits.iter().map(Vec::len).collect();
    sizes.sort_unstable();
    sizes
}

/// Check one run's transcript against Linnet.
fn check(run: &Run, mode: Mode, observed: &[Observed], expected: &[Expected]) {
    let context = format!("{} on {} against {}", run.command, run.input, run.answer);
    if expected.iter().any(|e| matches!(e, Expected::Refused)) {
        // Traces stops at digraphs; Linnet's view of them is checked by the other modes
        assert_eq!(mode, Mode::Traces, "{context}");
        let Some(Observed::Header { directed: true, .. }) = observed.first() else {
            panic!("only digraphs are refused: {context}");
        };
        return;
    }
    assert_eq!(observed.len(), expected.len(), "{context}: {observed:#?}");
    let mut hashes = vec![];
    for (o, e) in observed.iter().zip(expected) {
        match (o, e) {
            (
                Observed::Header { n, edges, directed },
                Expected::Header {
                    mode: m,
                    n: en,
                    edges: ee,
                    directed: ed,
                },
            ) => {
                assert_eq!((n, edges, directed), (en, ee, ed), "{context}");
                assert_eq!(*m, mode, "{context}");
            }
            (Observed::Degrees { degrees, labelorg }, Expected::Degrees(ranges)) => {
                let mut printed = vec![];
                for &(from, to, d) in ranges {
                    assert_eq!(from - labelorg, printed.len(), "{context}");
                    printed.extend((from..=to).map(|_| d));
                }
                assert_eq!(degrees, &printed, "{context}");
            }
            (Observed::Fixing, Expected::Fixing) => {}
            (
                Observed::Group {
                    symmetry,
                    graph,
                    labelorg,
                    labelled,
                },
                Expected::Group {
                    order,
                    orbits,
                    generators,
                },
            ) => {
                assert_eq!(symmetry.order, *order, "group size: {context}");
                assert_eq!(symmetry.node_orbits.len(), *orbits, "orbits: {context}");
                if *labelled && !generators.is_empty() {
                    let images: Vec<Vec<usize>> = generators
                        .iter()
                        .map(|cycles| {
                            let mut image: Vec<usize> = (0..graph.n).collect();
                            for cycle in cycles {
                                for (k, &a) in cycle.iter().enumerate() {
                                    image[a - labelorg] = cycle[(k + 1) % cycle.len()] - labelorg;
                                }
                            }
                            assert!(graph.is_automorphism(&image), "{cycles:?}: {context}");
                            image
                        })
                        .collect();
                    // nauty's generators lie in Linnet's group and generate one of its size
                    assert_eq!(group_order(&images, graph.n), *order, "{context}");
                }
            }
            (
                Observed::Orbits {
                    orbits,
                    labelorg,
                    labelled,
                },
                Expected::Orbits(printed),
            ) => {
                let printed = unlabel(printed, *labelorg);
                if *labelled {
                    assert_eq!(orbits, &printed, "{context}");
                } else {
                    assert_eq!(sizes(orbits), sizes(&printed), "{context}");
                }
            }
            (Observed::Hash(cert), Expected::Hash(hash)) => hashes.push((cert, hash)),
            (o, e) => panic!("{context}: Linnet observed {o:?} where nauty printed {e:?}"),
        }
    }
    for (a, b) in hashes.iter().zip(hashes.iter().skip(1)) {
        assert_eq!(a.0 == b.0, a.1 == b.1, "canonical graphs agree: {context}");
    }
}

/// Runs 1-61 of `runalltests`: every `dreadtest*` build in every mode, on all three scripts.
#[test]
fn dreadtest_runs_match_their_transcripts() {
    let runs: Vec<Run> = Run::all()
        .into_iter()
        .filter(|r| r.command.starts_with("./dreadtest"))
        .collect();
    assert_eq!(runs.len(), 61);
    let mut checked = BTreeSet::new();
    for run in &runs {
        let mode = Mode::of_command(&run.command);
        // the builds differ in word size and MAXN only, which Linnet has no counterpart of
        if checked.insert((run.input.clone(), run.answer.clone(), mode)) {
            let observed = Session::new(mode).run(&fixture(&run.input));
            let expected = parse_answers(&fixture(&run.answer));
            check(run, mode, &observed, &expected);
        }
    }
    // three scripts, each in dense, sparse and Traces mode
    assert_eq!(checked.len(), 9);
}

#[test]
fn kiss_matches_naurng() {
    // NEXTRAN and ranperm(perm, 12) after ran_init(1), printed by nauty 2.9.3's naurng.c
    let mut kiss = Kiss::new(1);
    let first = [kiss.next(), kiss.next(), kiss.next()];
    assert_eq!(
        first,
        [
            2409185123475645365,
            6606618937034311578,
            17769978625415337512
        ]
    );
    let perm = Kiss::new(1).permutation(12);
    assert_eq!(perm, [1, 6, 7, 11, 8, 4, 3, 0, 10, 2, 9, 5]);
}
