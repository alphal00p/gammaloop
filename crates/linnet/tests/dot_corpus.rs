//! Reads every DOT file of the repository's GammaLoop fixtures and examples.
//! Each must parse to the node, edge and half-edge counts pinned in a snapshot,
//! and both DOT writers must write each graph as DOT that reads back to a graph
//! with the same content, which is a fixed point of writing and reading.

use std::{
    fmt::Write,
    fs,
    path::{Path, PathBuf},
};

use linnet::parser::{set::DotGraphSet, DotGraph, HedgeParseError};

/// Directories, relative to the workspace root, whose `.dot` files form the
/// corpus. Nix provides them through linnet's runtime test sources.
const CORPUS: &[&str] = &[
    "crates/gammalooprs/tests/fixtures",
    "crates/linnest/examples",
    "examples/cli",
    "examples/notebooks/fixtures",
    "tests/resources/graphs",
];

/// Stale files that no DOT reader accepts, with the line and message of their
/// syntax error.
const INVALID: &[(&str, usize, &str)] = &[
    (
        "tests/resources/graphs/physical_3L_6photons_topology_A.dot",
        20,
        "`18s`, reads as two identifiers",
    ),
    (
        "tests/resources/graphs/qqx_aaa_ir_subtracted_physical.dot",
        14,
        "expected an attribute name or `]`, found `,`",
    ),
    (
        "tests/resources/graphs/qqx_aaa_subtracted.dot",
        13,
        "expected an attribute name or `]`, found `,`",
    ),
];

fn dot_files(dir: &Path, files: &mut Vec<PathBuf>) {
    let entries = fs::read_dir(dir).unwrap_or_else(|error| panic!("{}: {error}", dir.display()));
    for entry in entries {
        let path = entry.unwrap().path();
        if path.is_dir() {
            dot_files(&path, files);
        } else if path.extension().is_some_and(|extension| extension == "dot") {
            files.push(path);
        }
    }
}

/// The node and edge statements of a graph, independent of its indices.
fn content(graph: &DotGraph) -> (Vec<String>, Vec<String>) {
    let mut nodes: Vec<String> = (graph.iter_nodes())
        .map(|(_, _, node)| format!("{:?} {:?}", node.name, node.statements))
        .collect();
    let mut edges: Vec<String> = (graph.iter_edges())
        .map(|(_, _, edge)| format!("{:?}", edge.data.statements))
        .collect();
    nodes.sort();
    edges.sort();
    (nodes, edges)
}

#[test]
fn repository_dot_files_parse_and_round_trip_through_both_writers() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut files = Vec::new();
    for dir in CORPUS {
        dot_files(&root.join(dir), &mut files);
    }
    files.sort();
    assert!(files.len() >= 600, "found only {} DOT files", files.len());

    let (mut graphs, mut invalid) = (0, 0);
    // A line per file with the nodes/edges/half-edges of each graph: the round
    // trip below cannot tell a parse that is consistently wrong.
    let mut counts = String::new();
    for file in &files {
        let name = file.strip_prefix(&root).unwrap().to_str().unwrap();
        let parsed = DotGraphSet::from_file(file);
        if let Some(&(_, line, message)) = INVALID.iter().find(|(path, ..)| *path == name) {
            match parsed {
                Err(HedgeParseError::Syntax(error)) => {
                    assert_eq!(error.line, line, "{name}: {error}");
                    assert!(error.message.contains(message), "{name}: {error}");
                }
                other => panic!("{name}: expected a syntax error, got {other:?}"),
            }
            invalid += 1;
            continue;
        }
        let set = parsed.unwrap_or_else(|error| panic!("{name}: {error}"));
        write!(counts, "{name}").unwrap();
        for graph in set {
            graphs += 1;
            let sizes = (graph.n_live(), graph.n_edges(), graph.n_hedges());
            write!(counts, " {}/{}/{}", sizes.0, sizes.1, sizes.2).unwrap();
            let written = graph.debug_dot();
            let read = DotGraph::from_string(&written)
                .unwrap_or_else(|error| panic!("{name}: {error}\n{written}"));
            assert_eq!(
                (read.n_live(), read.n_edges(), read.n_hedges()),
                sizes,
                "{name}"
            );
            assert_eq!(read.global_data, graph.global_data, "{name}");
            assert_eq!(content(&read), content(&graph), "{name}");
            assert_eq!(read, read.clone().back_and_forth_dot(), "{name}");

            let mut io = Vec::new();
            graph.write_io(&mut io).unwrap();
            let read_io = DotGraph::from_string(String::from_utf8(io).unwrap())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            assert_eq!(read_io, read, "{name}: the two writers disagree");
        }
        counts.push('\n');
    }
    assert_eq!(invalid, INVALID.len());
    assert!(graphs >= 850, "read only {graphs} graphs");
    insta::assert_snapshot!("graph_counts", counts);
}
