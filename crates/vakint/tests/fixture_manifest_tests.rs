//! The vendored fixture files of the factorized entries' tests against their
//! sha256 manifest (tests/fixtures/FIXTURE_MANIFEST.tsv): a file edited, added or
//! removed without its manifest row fails here, before any value test can be
//! affected by it. With the producers' directories given, the explicit test
//! also checks that every copied file is still byte for byte its original.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

#[path = "support/wh_fixtures.rs"]
#[allow(dead_code)]
mod wh_fixtures;

/// The directories whose files the manifest lists, exactly (relative to the crate).
const DIRECTORIES: [&str; 3] = [
    "src/p_integral/factorized/fixtures",
    "tests/fixtures/wh_ports",
    "tests/fixtures/direct_differential",
];

struct Row {
    path: String,
    sha256: String,
    origin: String,
}

fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn manifest() -> Vec<Row> {
    let path = crate_dir().join("tests/fixtures/FIXTURE_MANIFEST.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 3, "a manifest row has 3 fields: {line}");
            Row {
                path: fields[0].to_string(),
                sha256: fields[1].to_string(),
                origin: fields[2].to_string(),
            }
        })
        .collect()
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// The rows whose file (read by `read`) does not have the listed sha256.
fn mismatches(rows: &[Row], read: &dyn Fn(&str) -> Option<Vec<u8>>) -> Vec<String> {
    rows.iter()
        .filter_map(|row| match read(&row.path) {
            None => Some(format!("{}: missing", row.path)),
            Some(bytes) => {
                let got = sha256(&bytes);
                (got != row.sha256)
                    .then(|| format!("{}: sha256 {got}, listed {}", row.path, row.sha256))
            }
        })
        .collect()
}

fn read_in_crate(path: &str) -> Option<Vec<u8>> {
    std::fs::read(crate_dir().join(path)).ok()
}

/// Every listed file has its sha256; the three directories hold exactly the
/// listed files; no path is listed twice.
#[test]
fn the_vendored_fixtures_match_their_manifest() {
    let rows = manifest();
    assert_eq!(rows.len(), 55, "manifest rows");
    let listed = rows.iter().map(|r| r.path.clone()).collect::<BTreeSet<_>>();
    assert_eq!(listed.len(), rows.len(), "a path is listed twice");
    let bad = mismatches(&rows, &read_in_crate);
    assert!(
        bad.is_empty(),
        "fixtures differ from the manifest:\n{}",
        bad.join("\n")
    );
    let mut present = BTreeSet::new();
    for directory in DIRECTORIES {
        for entry in std::fs::read_dir(crate_dir().join(directory)).unwrap() {
            let name = entry.unwrap().file_name().into_string().unwrap();
            present.insert(format!("{directory}/{name}"));
        }
    }
    let unlisted = present.difference(&listed).collect::<Vec<_>>();
    let absent = listed.difference(&present).collect::<Vec<_>>();
    assert!(
        unlisted.is_empty() && absent.is_empty(),
        "files without a manifest row: {unlisted:?}; rows outside the directories: {absent:?}"
    );
}

/// The check itself: one changed byte, a missing file and a wrong listed hash
/// are each reported.
#[test]
fn the_manifest_check_reports_a_changed_file() {
    let rows = manifest();
    let first = rows[0].path.clone();
    let changed = |path: &str| {
        let mut bytes = read_in_crate(path)?;
        if path == first {
            bytes[0] ^= 1;
        }
        Some(bytes)
    };
    assert_eq!(mismatches(&rows, &changed).len(), 1, "one changed byte");
    let missing = |path: &str| {
        if path == first {
            None
        } else {
            read_in_crate(path)
        }
    };
    assert_eq!(mismatches(&rows, &missing).len(), 1, "a missing file");
    let mut wrong = manifest();
    wrong[1].sha256 = "0".repeat(64);
    assert_eq!(
        mismatches(&wrong, &read_in_crate).len(),
        1,
        "a wrong listed hash"
    );
}

/// Explicit (`--run-ignored only`): every row copied from a producer's
/// directory equals its original there byte for byte. Needs
/// VAKINT_TEST_WH_FIXTURES, VAKINT_TEST_L3_HARD_DUMP,
/// VAKINT_TEST_REVIEW_PROBE, VAKINT_TEST_HOLE_CUT_CASES and
/// VAKINT_TEST_RECORD_ROOTS (its first root).
#[test]
#[ignore = "reads the producers' directories"]
fn the_vendored_fixtures_equal_their_originals() {
    let roots = [
        ("wh-fixtures:", "VAKINT_TEST_WH_FIXTURES", 30usize),
        ("l3-hard-dump:", "VAKINT_TEST_L3_HARD_DUMP", 5),
        ("review-probe:", "VAKINT_TEST_REVIEW_PROBE", 3),
        ("hole-cut-cases:", "VAKINT_TEST_HOLE_CUT_CASES", 6),
        ("record-roots:", "VAKINT_TEST_RECORD_ROOTS", 2),
    ];
    let rows = manifest();
    for (prefix, variable, expected) in roots {
        let root = wh_fixtures::required_env(variable, "the directory of the originals");
        // VAKINT_TEST_RECORD_ROOTS lists several roots: its first one
        let root = root
            .to_string_lossy()
            .split(':')
            .next()
            .unwrap_or("")
            .to_string();
        let copied = rows
            .iter()
            .filter_map(|r| {
                r.origin.strip_prefix(prefix).map(|original| Row {
                    path: original.to_string(),
                    sha256: r.sha256.clone(),
                    origin: r.path.clone(),
                })
            })
            .collect::<Vec<_>>();
        let read = |path: &str| std::fs::read(Path::new(&root).join(path)).ok();
        let bad = mismatches(&copied, &read);
        assert!(
            bad.is_empty(),
            "{prefix} originals differ:\n{}",
            bad.join("\n")
        );
        assert_eq!(copied.len(), expected, "{prefix} rows compared");
        eprintln!("ORIGINALS {prefix} {} files identical", copied.len());
    }
}
