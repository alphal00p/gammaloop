//! nauty's graph6, sparse6 (also incremental) and digraph6 line formats, as described in
//! nauty's `formats.txt`.

use super::graph::Graph;

/// The graphs of a nauty output stream, one per line. Lines starting with `>` are the
/// utilities' messages on stderr, which `runalltests` merges into the answers.
pub fn read_graphs(text: &str) -> Vec<Graph> {
    let mut graphs: Vec<Graph> = vec![];
    for line in text.lines() {
        let line = line
            .trim_end()
            .trim_start_matches(">>graph6<<")
            .trim_start_matches(">>sparse6<<")
            .trim_start_matches(">>digraph6<<");
        if line.is_empty() || line.starts_with('>') {
            continue;
        }
        let graph = parse_graph(line, graphs.last());
        graphs.push(graph);
    }
    graphs
}

/// One graph in any of the three formats; an incremental sparse6 line (`;`) is relative to
/// `previous`.
pub fn parse_graph(line: &str, previous: Option<&Graph>) -> Graph {
    let bytes = line.as_bytes();
    match bytes[0] {
        b'&' => digraph6(&bytes[1..]),
        b':' => sparse6(&bytes[1..]),
        b';' => {
            let previous = previous.expect("an incremental sparse6 line follows a graph");
            let delta = sparse6_edges(previous.n, &bytes[1..]);
            let mut edges: Vec<_> = previous.edges.iter().map(|&(a, b)| key(a, b)).collect();
            for (a, b) in delta {
                let e = key(a, b);
                match edges.iter().position(|&f| f == e) {
                    Some(i) => {
                        edges.swap_remove(i);
                    }
                    None => edges.push(e),
                }
            }
            edges.sort_unstable();
            Graph::undirected(previous.n, edges)
        }
        _ => graph6(bytes),
    }
}

fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// `N(n)` followed by the data bytes.
fn order(bytes: &[u8]) -> (usize, &[u8]) {
    let value = |bs: &[u8]| {
        bs.iter()
            .fold(0usize, |acc, &b| (acc << 6) | (b - 63) as usize)
    };
    if bytes[0] != 126 {
        (value(&bytes[..1]), &bytes[1..])
    } else if bytes[1] != 126 {
        (value(&bytes[1..4]), &bytes[4..])
    } else {
        (value(&bytes[2..8]), &bytes[8..])
    }
}

/// The bits of `R(x)`, most significant first.
fn bits(bytes: &[u8]) -> impl Iterator<Item = bool> + '_ {
    bytes
        .iter()
        .flat_map(|&b| (0..6).rev().map(move |i| ((b - 63) >> i) & 1 == 1))
}

fn graph6(bytes: &[u8]) -> Graph {
    let (n, data) = order(bytes);
    let pairs = (1..n).flat_map(|j| (0..j).map(move |i| (i, j)));
    let edges = pairs.zip(bits(data)).filter(|&(_, b)| b).map(|(e, _)| e);
    Graph::undirected(n, edges.collect())
}

fn digraph6(bytes: &[u8]) -> Graph {
    let (n, data) = order(bytes);
    let pairs = (0..n).flat_map(|i| (0..n).map(move |j| (i, j)));
    let arcs = pairs.zip(bits(data)).filter(|&(_, b)| b).map(|(e, _)| e);
    Graph::directed(n, arcs.collect())
}

fn sparse6(bytes: &[u8]) -> Graph {
    let (n, data) = order(bytes);
    Graph::undirected(n, sparse6_edges(n, data))
}

/// The edge list of sparse6: `b[i] x[i]` pairs, `v` advancing on `b = 1` or jumping to a larger
/// `x`; an incomplete pair or one beyond `n` is padding.
fn sparse6_edges(n: usize, data: &[u8]) -> Vec<(usize, usize)> {
    let k = usize::BITS - n.saturating_sub(1).leading_zeros();
    let bits: Vec<bool> = bits(data).collect();
    let mut edges = vec![];
    let mut v = 0;
    for chunk in bits.chunks_exact(1 + k as usize) {
        let x = chunk[1..]
            .iter()
            .fold(0usize, |acc, &b| (acc << 1) | b as usize);
        if chunk[0] {
            v += 1;
        }
        if x >= n || v >= n {
            break;
        }
        if x > v {
            v = x;
        } else {
            edges.push((x, v));
        }
    }
    edges
}

#[test]
fn formats_txt_examples() {
    // the examples of nauty's formats.txt
    let mut g = parse_graph("DQc", None);
    g.edges.sort_unstable();
    assert_eq!((g.n, g.edges), (5, vec![(0, 2), (0, 4), (1, 3), (3, 4)]));
    let s = parse_graph(":Fa@x^", None);
    assert_eq!((s.n, s.edges), (7, vec![(0, 1), (0, 2), (1, 2), (5, 6)]));
    let d = parse_graph("&DI?AO?", None);
    assert!(d.directed);
    assert_eq!((d.n, d.edges), (5, vec![(0, 2), (0, 4), (3, 1), (3, 4)]));
}
