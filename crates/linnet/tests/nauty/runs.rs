//! The runs of nauty's `runalltests` and how each one is translated.

use std::path::PathBuf;

/// The nauty-derived fixtures, with nauty's licence and attribution.
pub fn fixture_path(name: &str) -> PathBuf {
    [env!("CARGO_MANIFEST_DIR"), "tests/resources/nauty", name]
        .iter()
        .collect()
}

pub fn fixture(name: &str) -> String {
    let path = fixture_path(name);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Whether to run the full-size families, which take seconds in optimised builds and minutes at
/// opt-level 0. Debug assertions alone do not tell the two apart: CI's `ci-optim` profile keeps
/// them but optimises. So only Cargo's `dev` profile, whose binaries live under `target/debug`,
/// runs the reduced sizes.
pub fn full_size() -> bool {
    !cfg!(debug_assertions)
        || std::env::current_exe()
            .ok()
            .and_then(|exe| Some(exe.parent()?.parent()?.file_name()? != "debug"))
            .unwrap_or(true)
}

/// One `runonetest` line: a command (a pipeline), its stdin and the expected output.
#[derive(Clone, Debug)]
pub struct Run {
    pub command: String,
    /// Empty when the command reads nothing.
    pub input: String,
    pub answer: String,
}

impl Run {
    /// Every run of `runalltests.in`, including the ones only made with 128-bit support.
    pub fn all() -> Vec<Run> {
        fixture("runalltests.in")
            .lines()
            .filter(|l| l.starts_with("runonetest ") || l.starts_with("runonetestQ "))
            .map(|line| {
                let words = shell_words(line);
                let [_, command, input, answer] = &words[..] else {
                    panic!("unexpected test line {line}");
                };
                Run {
                    command: command.clone(),
                    input: input.clone(),
                    answer: answer.clone(),
                }
            })
            .collect()
    }
}

/// Words of a shell line, double quotes grouping.
fn shell_words(line: &str) -> Vec<String> {
    let mut words = vec![];
    let mut rest = line.trim();
    while !rest.is_empty() {
        let (word, tail) = match rest.strip_prefix('"') {
            Some(quoted) => quoted.split_once('"').expect("closing quote"),
            None => rest.split_once(' ').unwrap_or((rest, "")),
        };
        words.push(word.to_string());
        rest = tail.trim_start();
    }
    words
}

/// How a utility run is translated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Translation {
    /// An in-tree test, run by default (at reduced size in unoptimised builds, see [`full_size`]).
    Test(&'static str),
    /// An ignored test that runs the pipeline on a nauty build (`NAUTY_DIR`).
    Oracle(&'static str),
    /// Nothing algorithm-independent to compare with Linnet.
    NotApplicable(&'static str),
}

use Translation::{NotApplicable, Oracle, Test};

/// The utility pipelines of `runalltests`, by answer file. The `dreadtest*` runs are all
/// checked by `dreadnaut::dreadtest_runs_match_their_transcripts`.
const UTILITIES: &[(&str, Translation)] = &[
    ("nautesta.ans", Oracle("oracle::geng_triangle_free_graphs")),
    (
        "nautestb1.ans",
        Test("utilities::labelg_canonical_labellings"),
    ),
    (
        "nautestb2.ans",
        Test("utilities::labelg_canonical_labellings"),
    ),
    (
        "nautestb3.ans",
        Test("utilities::labelg_canonical_labellings"),
    ),
    (
        "nautestc.ans",
        NotApplicable("random cubic graphs tabulated by pentagon counts"),
    ),
    ("nautestd.ans", Test("utilities::gentourng_tournaments")),
    ("nauteste.ans", Test("utilities::gentreeg_orbit_table")),
    ("nautestf.ans", Test("utilities::directg_orientations")),
    ("nautestg.ans", Test("utilities::vcolg_colourings")),
    ("nautesth.ans", Test("utilities::genposetg_posets")),
    ("nautesti.ans", Test("utilities::genspecialg_graphs")),
    ("nautestj.ans", Test("utilities::genspecialg_digraphs")),
    ("nautestk.ans", Test("utilities::assembleg_disjoint_unions")),
    (
        "nautestl.ans",
        NotApplicable("countg cycle, clique and independence counts"),
    ),
    (
        "nautestm.ans",
        Test("utilities::genquarticg_quartic_graphs"),
    ),
    ("nautestn.ans", Test("utilities::gengl_two_regular_graphs")),
    (
        "nautesto.ans",
        NotApplicable("labelled output of biplabg and NRswitchg"),
    ),
    (
        "nautestp.ans",
        NotApplicable("labelled output of addedgeg and deledgeg"),
    ),
    ("nautestq.ans", Test("utilities::genktreeg_complements")),
    ("nautestr.ans", NotApplicable("countg k-tree recognition")),
    (
        "nautests.ans",
        Test("utilities::shortg_random_relabellings"),
    ),
    (
        "nautestt.ans",
        NotApplicable("countg chromatic numbers and hexagon counts"),
    ),
    (
        "nautestu.ans",
        NotApplicable("countg chromatic index and class"),
    ),
    (
        "nautestv.ans",
        NotApplicable("random subgraphs tabulated by colouring and connectivity"),
    ),
    (
        "nautestw.ans",
        NotApplicable("random graphs tabulated by connectivity"),
    ),
    (
        "nautestx.ans",
        NotApplicable("random digraphs tabulated by connectivity"),
    ),
    (
        "nautesty.ans",
        Oracle("oracle::uniqg_chordal_bipartite_relabellings"),
    ),
    (
        "nautestz.ans",
        NotApplicable("labelled 2-cores tabulated by degree"),
    ),
    (
        "nautestaa.ans",
        NotApplicable("labelled output of nbrhoodg and newedgeg"),
    ),
    (
        "nautestbb.ans",
        NotApplicable("genbg's construction ancestors as matrices"),
    ),
    (
        "nautestcc.ans",
        NotApplicable("line graphs tabulated by order and size"),
    ),
    (
        "nautestdd.ans",
        NotApplicable("countg diamond, hexagon and pentagon counts"),
    ),
    (
        "nautestee.ans",
        NotApplicable("subdivisions tabulated by degree statistics"),
    ),
    (
        "nautestff.ans",
        Test("utilities::watercluster_orientations"),
    ),
    ("nautestgg.ans", Test("utilities::multig_multigraphs")),
    (
        "nautesthh.ans",
        NotApplicable("planarity filter and Cartesian product"),
    ),
    (
        "nautestii.ans",
        Test("utilities::distgraphg_distance_three_graphs"),
    ),
];

#[test]
fn every_run_is_accounted_for() {
    let runs = Run::all();
    assert_eq!(runs.len(), 98);
    let dreadtests = runs
        .iter()
        .filter(|r| r.command.starts_with("./dreadtest"))
        .count();
    assert_eq!(dreadtests, 61);
    let utilities: Vec<&str> = runs
        .iter()
        .filter(|r| !r.command.starts_with("./dreadtest"))
        .map(|r| r.answer.as_str())
        .collect();
    let table: Vec<&str> = UTILITIES.iter().map(|(a, _)| *a).collect();
    assert_eq!(utilities, table);
    let count = |f: fn(&Translation) -> bool| UTILITIES.iter().filter(|(_, t)| f(t)).count();
    assert_eq!(count(|t| matches!(t, Test(_))), 18);
    assert_eq!(count(|t| matches!(t, Oracle(_))), 2);
    assert_eq!(count(|t| matches!(t, NotApplicable(_))), 17);
    // the tests named in the table exist
    let sources = [
        ("utilities", include_str!("utilities.rs")),
        ("oracle", include_str!("oracle.rs")),
    ];
    for (_, translation) in UTILITIES {
        let (Test(name) | Oracle(name)) = translation else {
            continue;
        };
        let (module, test) = name.split_once("::").unwrap();
        assert!(
            sources
                .iter()
                .any(|(m, s)| *m == module && s.contains(&format!("fn {test}()"))),
            "{name}"
        );
    }
}
