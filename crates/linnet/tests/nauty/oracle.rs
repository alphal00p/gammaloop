//! Runs that need nauty itself: pipelines too large to freeze as fixtures, executed on a nauty
//! 2.9.3 build and compared with Linnet graph by graph. The tests are ignored by default and
//! fail unless `NAUTY_DIR` names the build directory (where `make` put `geng` and the other
//! utilities):
//!
//! ```text
//! NAUTY_DIR=/path/to/nauty2_9_3 cargo nextest run --release -p linnet --test nauty \
//!     --run-ignored only oracle::
//! ```

use std::{
    collections::HashMap,
    io::Write,
    path::PathBuf,
    process::{Command, Stdio},
};

use super::{
    formats::read_graphs,
    graph::{Cert, Graph},
    runs::fixture,
};

/// The nauty build. Selecting these ignored tests without one fails rather than passing having
/// compared nothing.
fn nauty_dir() -> PathBuf {
    std::env::var_os("NAUTY_DIR")
        .map(PathBuf::from)
        .expect("NAUTY_DIR must name a nauty 2.9.3 build to compare with nauty")
}

/// Run a `runalltests`-style pipeline (`./geng ... | ./ranlabg ...`) in the nauty build, with
/// `input` on stdin; returns stdout and stderr.
fn pipeline(dir: &PathBuf, command: &str, input: &str) -> (String, String) {
    let mut child = Command::new("sh")
        .arg("-c")
        .arg(command)
        .current_dir(dir)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| panic!("{command}: {e}"));
    let mut stdin = child.stdin.take().unwrap();
    let input = input.to_string();
    let writer = std::thread::spawn(move || stdin.write_all(input.as_bytes()));
    let output = child.wait_with_output().unwrap();
    writer.join().unwrap().unwrap();
    assert!(output.status.success(), "{command} failed");
    (
        String::from_utf8(output.stdout).unwrap(),
        String::from_utf8(output.stderr).unwrap(),
    )
}

/// How often each certificate occurs.
fn census(graphs: &[Graph]) -> HashMap<Cert, usize> {
    let mut census = HashMap::new();
    for g in graphs {
        *census.entry(g.certificate()).or_insert(0) += 1;
    }
    census
}

/// The number in a `>Z ... N graphs <what>` summary line.
fn summary_count(text: &str, what: &str) -> usize {
    let line = text.lines().find(|l| l.contains(what)).unwrap();
    let before = &line[..line.find(what).unwrap()];
    before.split_whitespace().last().unwrap().parse().unwrap()
}

/// Run 62, `geng -ud1D7t 11`: 92779 triangle-free graphs on 11 vertices with degrees 1 to 7.
/// Linnet separates all of them, and a random relabelling of the stream merges back onto them.
#[test]
#[ignore = "needs a nauty 2.9.3 build in NAUTY_DIR"]
fn geng_triangle_free_graphs() {
    let dir = nauty_dir();
    let expected = summary_count(&fixture("nautesta.ans"), " graphs generated");
    assert_eq!(expected, 92779);
    let (graphs, _) = pipeline(&dir, "./geng -q -d1D7t 11", "");
    let graphs = census(&read_graphs(&graphs));
    assert_eq!(graphs.len(), expected);
    let (relabelled, _) = pipeline(&dir, "./geng -q -d1D7t 11 | ./ranlabg -q", "");
    let relabelled = census(&read_graphs(&relabelled));
    assert_eq!(relabelled, graphs);
}

/// Run 88, `geng -q -bT 12 | ranlabg -qS7654 -m5 | uniqg -uy`: five random relabellings of each
/// of the 90953 chordal bipartite graphs on 12 vertices, which Linnet groups back into 90953
/// classes of five.
#[test]
#[ignore = "needs a nauty 2.9.3 build in NAUTY_DIR"]
fn uniqg_chordal_bipartite_relabellings() {
    let dir = nauty_dir();
    let answer = fixture("nautesty.ans");
    assert!(answer.contains("454765 graphs read from stdin, 90953 unique"));
    let (copies, _) = pipeline(&dir, "./geng -q -bT 12 | ./ranlabg -qS7654 -m5", "");
    let copies = read_graphs(&copies);
    assert_eq!(copies.len(), 454765);
    let classes = census(&copies);
    assert_eq!(classes.len(), 90953);
    assert!(classes.values().all(|&c| c == 5));
}

/// `|Aut|`, vertex orbits and (for graphs) edge orbits of every graph against countg, on graph
/// streams of the suite: the random cubic graphs of run 66, the (10,12) graphs of run 75, the
/// random graphs and digraphs of runs 86 and 87, and the planar product graph of run 97.
#[test]
#[ignore = "needs a nauty 2.9.3 build in NAUTY_DIR"]
fn countg_group_sizes_and_orbits() {
    let dir = nauty_dir();
    for (stream, distinct) in [
        ("./genrang -S4321 -r3 114 100", false),
        ("./geng -q 10 12", true),
        ("./genrang -q -P3 -S1234 70 1000", false),
        ("./genrang -q -P4 -z -S1234 70 500", false),
        (
            "./geng -qtfc 9 | ./planarg -q | ./pickg -qp1:2 | ./productg -c",
            false,
        ),
    ] {
        let (text, _) = pipeline(&dir, &format!("{stream} 2>/dev/null"), "");
        let graphs = read_graphs(&text);
        assert!(!graphs.is_empty(), "{stream}");
        let directed = graphs[0].directed;
        let keys = if directed { "--a,o" } else { "--a,o,O" };
        let (table, _) = pipeline(&dir, &format!("./countg -q -V -1 {keys}"), &text);
        let rows: Vec<Vec<u128>> = table
            .lines()
            .filter_map(|l| l.strip_prefix("Graph "))
            .map(|l| {
                let (_, values) = l.split_once(" : ").unwrap();
                values
                    .split_whitespace()
                    .map(|v| v.parse().unwrap())
                    .collect()
            })
            .collect();
        assert_eq!(rows.len(), graphs.len(), "{stream}");
        for (g, row) in graphs.iter().zip(&rows) {
            let symmetry = g.symmetry();
            let mut found = vec![symmetry.order, symmetry.node_orbits.len() as u128];
            if !directed {
                found.push(symmetry.edge_orbits as u128);
            }
            assert_eq!(&found, row, "{stream}: {g:?}");
        }
        if distinct {
            assert_eq!(census(&graphs).len(), graphs.len(), "{stream}");
        }
    }
}

/// The full-size versions of runs 67 and 69: all 6880 tournaments on 8 vertices, and the
/// orientations written by directg before conversing, which Linnet classes exactly as labelg's
/// canonical forms after conversing.
#[test]
#[ignore = "needs a nauty 2.9.3 build in NAUTY_DIR"]
fn full_tournament_and_orientation_families() {
    let dir = nauty_dir();
    let (tournaments, _) = pipeline(&dir, "./gentourng -q -z 8", "");
    assert_eq!(census(&read_graphs(&tournaments)).len(), 6880);
    let (oriented, _) = pipeline(&dir, "./geng -c -qp 8 9 | ./directg -q -o", "");
    let oriented = census(&read_graphs(&oriented));
    assert_eq!(oriented.len(), 36880);
    let (conversed, _) = pipeline(
        &dir,
        "./geng -c -qp 8 9 | ./directg -q -o | ./converseg -q | ./labelg -q",
        "",
    );
    assert_eq!(census(&read_graphs(&conversed)), oriented);
}
