//! Differential oracles for diagram generation.
//!
//! Every case of [`cases::CASES`] runs with Named and Minimal output at 1, 8 and 32 threads.
//! Each policy compares full canonical keys and physical payloads across threads; the policies
//! may differ only in names and optional source identities. Ungrouped references validate each
//! group source index.
//! The runs must agree exactly (diagram
//! set, |Aut|, overall factors, frames, numerators, ids, groups, names and output order), and their
//! [`census::Summary`] must equal the case's row below. Rows include a digest of the diagram ids,
//! which depend on the diagrams alone, not on the process that computes them
//! (`ids_do_not_depend_on_threads_or_symbol_registration` checks this across fresh processes). The
//! `phi` rows are fixed by independent oracles: OEIS A005967 counts the connected cubic vacuum
//! graphs, zero-dimensional field theory gives every sum of 1/|Aut|, and FeynGraph and Graphica
//! agree with every count and |Aut| multiset (Graphica alone for the phi^3 vacua at six and seven
//! loops, where FeynGraph did not finish). The other rows were recorded with this harness; their
//! counts reproduce the profile of the native generator (e.g. 15,253 and 686,501 diagrams for
//! gg->gg at two and three loops). `failing` cases fail on the base: each must fail at every thread
//! count with an error that contains the recorded fragments, so that a change neither hides nor
//! fixes one silently; a fix records the case's row instead. The fragments leave out diagram and
//! candidate numbers, because parallel runs report whichever failure they meet first.
//!
//! Heavy cases are ignored by default; `just generation-oracle-heavy` runs them. Environment:
//! `FEYNKIT_ORACLE_THREADS=1,8,32` sets the thread counts, `FEYNKIT_ORACLE_RECORD=1` prints each
//! case's row instead of comparing it, and `FEYNKIT_ORACLE_DUMP=<dir>` writes every census as
//! text (including their complete momentum bases), for diffing two builds.

mod cases;
mod census;

use std::{collections::BTreeSet, process::Command};

use cases::{CASES, Case, REFERENCE};
use census::{Census, Scope, Summary};
use feynkit_generator::{GenerationOutput, NumeratorGrouping};
use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::domains::rational::Rational;

fn threads() -> Vec<usize> {
    std::env::var("FEYNKIT_ORACLE_THREADS").map_or(vec![1, 8, 32], |threads| {
        threads.split(',').map(|t| t.parse().unwrap()).collect()
    })
}

fn census(case: &Case, threads: usize) -> Census {
    census_with_output(case, threads, GenerationOutput::default())
}

fn census_with_output(case: &Case, threads: usize, output: GenerationOutput) -> Census {
    let model = case.model();
    let options = case.options(&model, threads).output(output);
    let result = case
        .generate_with(model, &options)
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
    assert!(result.report.completed, "{} did not complete", case.name);
    Census::of(&result)
}

/// A recorded [`Summary`], with the sums of 1/|Aut| and of the evaluated overall factors as
/// numerator and denominator in lowest terms.
struct Row {
    diagrams: usize,
    topologies: usize,
    assignments: usize,
    inverse_symmetry: (u64, u64),
    weights: (i64, i64),
    symmetry_factors: &'static [(u64, usize)],
    content: u64,
    numerators: u64,
    named: u64,
    ids: u64,
}

impl From<&Row> for Summary {
    fn from(row: &Row) -> Summary {
        Summary {
            diagrams: row.diagrams,
            topologies: row.topologies,
            assignments: row.assignments,
            inverse_symmetry: Rational::from(row.inverse_symmetry),
            weights: Rational::from(row.weights),
            symmetry_factors: row.symmetry_factors.to_vec(),
            content: row.content,
            numerators: row.numerators,
            named: row.named,
            ids: row.ids,
        }
    }
}

/// Run `name` at every thread count, require identical censuses and compare the summary with
/// `row`.
fn check(name: &str, tier: &str, row: Row) {
    let case = Case::named(name);
    let mut threads = threads().into_iter();
    let first = threads.next().expect("a thread count");
    let reference = census(case, first);
    let minimal = census_with_output(case, first, GenerationOutput::Minimal);
    reference.assert_minimal(&minimal, name);
    for count in threads {
        let other = census(case, count);
        reference.assert_same(
            &other,
            Scope::Exact,
            &format!("{name}: Named {first} vs {count} threads"),
        );
        drop(other);
        let other = census_with_output(case, count, GenerationOutput::Minimal);
        minimal.assert_same(
            &other,
            Scope::Exact,
            &format!("{name}: Minimal {first} vs {count} threads"),
        );
    }
    drop(minimal);
    // An independently ungrouped run retains even members that grouping discarded.
    // Bind Named identities to full physical sources; the policy comparison above
    // requires Minimal to preserve the same request-local indices, ratios and factors.
    let model = case.model();
    let options = case.options(&model, first);
    // The source name/index assertion is specific to these selector-free fixtures.
    let serialized = serde_json::to_value(&options).unwrap();
    assert_eq!(
        serialized["output"],
        serde_json::json!({"named": {"prefix": "D"}})
    );
    for field in ["selected_diagram_ids", "selected_diagram_names"] {
        assert!(serialized.get(field).unwrap().is_null(), "{name}: {field}");
    }
    for field in ["vetoed_diagram_ids", "vetoed_diagram_names"] {
        assert_eq!(serialized[field], serde_json::json!([]), "{name}: {field}");
    }
    let ungrouped = case
        .generate_with(
            model.clone(),
            &options
                .numerator_grouping(NumeratorGrouping::None)
                .filter_zero_color(false),
        )
        .unwrap();
    reference.assert_group_sources(&Census::of(&ungrouped), name);
    drop(ungrouped);
    if let Ok(directory) = std::env::var("FEYNKIT_ORACLE_DUMP") {
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(format!("{directory}/{name}.census"), reference.dump()).unwrap();
    }
    let summary = reference.summary();
    if std::env::var_os("FEYNKIT_ORACLE_RECORD").is_some() {
        println!("    {tier} {name}: {summary};");
    } else {
        assert_eq!(
            summary,
            Summary::from(&row),
            "{name}: summary differs from its oracle row"
        );
    }
}

// Errors of the `failing` cases, as fragments.
/// Grouping QCD numerators: spenso cannot parse the adjoint colour dimension `Nc^2 - 1`.
const COLOUR_DIMENSION: &[&str] = &[
    "failed to evaluate numerator tensors for diagram ",
    "cannot parse tensor network: ",
    "Dimension is not a var or int: -1+Nc^2",
];
/// Grouping QED numerators at two loops: a sample evaluates to an unmaterialized scaled tensor.
const SCALED_TENSOR: &[&str] = &[
    "failed to evaluate numerator tensors for diagram ",
    "result is an unmaterialized scaled tensor",
];
/// Grouping SM numerators: the square of a metric on fundamental colour indices.
const NON_SELF_DUAL_POWER: &[&str] = &[
    "failed to evaluate numerator tensors for diagram ",
    "invalid power of non-self-dual tensor",
];

/// Run `name` at every thread count and require it to fail with an error containing every
/// fragment of `error`.
fn fails(name: &str, error: &[&str]) {
    let case = Case::named(name);
    for count in threads() {
        for output in [GenerationOutput::default(), GenerationOutput::Minimal] {
            let model = case.model();
            let options = case.options(&model, count).output(output.clone());
            match case.generate_with(model, &options) {
                Ok(_) => panic!(
                    "{name} generates at {count} threads ({output:?}): record its oracle row"
                ),
                Err(found) => {
                    let found = found.to_string();
                    assert!(
                        error.iter().all(|fragment| found.contains(fragment)),
                        "{name} at {count} threads ({output:?}) fails with {found:?}, not {error:?}"
                    );
                }
            }
        }
    }
}

macro_rules! oracle_test {
    (fast $name:ident $row:expr) => {
        #[test]
        #[allow(non_snake_case)]
        fn $name() {
            check(stringify!($name), "fast", $row);
        }
    };
    (heavy $name:ident $row:expr) => {
        #[test]
        #[ignore = "heavy oracle tier: run with `just generation-oracle-heavy`"]
        #[allow(non_snake_case)]
        fn $name() {
            check(stringify!($name), "heavy", $row);
        }
    };
    (failing $name:ident $error:expr) => {
        #[test]
        #[allow(non_snake_case)]
        fn $name() {
            fails(stringify!($name), $error);
        }
    };
}

macro_rules! oracle {
    ($($tier:ident $name:ident: $row:expr;)*) => {
        $(oracle_test!($tier $name $row);)*
        /// The cases of the oracle table.
        const ORACLE: &[&str] = &[$(stringify!($name)),*];
    };
}

oracle! {
    fast phi3_vac_L2: Row { diagrams: 2, topologies: 2, assignments: 2, inverse_symmetry: (5, 24), weights: (5, 24), symmetry_factors: &[(8, 1), (12, 1)], content: 0x75492f43d8fc3c82, numerators: 0xf818ded4decf2e81, named: 0xa8bfd0dc3bc692fe, ids: 0x84de875f268883b3 };
    fast phi3_vac_L3: Row { diagrams: 5, topologies: 5, assignments: 5, inverse_symmetry: (5, 16), weights: (5, 16), symmetry_factors: &[(8, 1), (16, 2), (24, 1), (48, 1)], content: 0x54ca4e63ab0c62fd, numerators: 0xaa97877126478b1a, named: 0xfff6529fca122d30, ids: 0x03d775b39d3d1d6a };
    fast phi3_vac_L4: Row { diagrams: 17, topologies: 17, assignments: 17, inverse_symmetry: (1105, 1152), weights: (1105, 1152), symmetry_factors: &[(8, 3), (12, 1), (16, 5), (32, 4), (48, 2), (72, 1), (128, 1)], content: 0x7b4cd459d7bc1a0c, numerators: 0x7d649ade64818e41, named: 0xee49d88f08ae4dd1, ids: 0x79a4563b066f205d };
    fast phi3_vac_L5: Row { diagrams: 71, topologies: 71, assignments: 71, inverse_symmetry: (565, 128), weights: (565, 128), symmetry_factors: &[(4, 4), (8, 10), (12, 1), (16, 21), (32, 17), (48, 1), (64, 11), (96, 2), (128, 2), (256, 2)], content: 0x78f7d529862a12a8, numerators: 0x6d848bd898887500, named: 0x733e25ec335d7523, ids: 0x74084cf17998a627 };
    fast phi3_vac_L6: Row { diagrams: 388, topologies: 388, assignments: 388, inverse_symmetry: (82825, 3072), weights: (82825, 3072), symmetry_factors: &[(2, 6), (4, 25), (6, 2), (8, 60), (12, 1), (16, 88), (20, 2), (32, 89), (48, 5), (64, 67), (120, 1), (128, 27), (256, 6), (320, 2), (384, 2), (512, 4), (3072, 1)], content: 0xe012291e41427b03, numerators: 0xccd2d268f51c99e6, named: 0x23df777b7aced522, ids: 0x9cc74c2176b66c91 };
    fast phi3_vac_L7: Row { diagrams: 2592, topologies: 2592, assignments: 2592, inverse_symmetry: (19675, 96), weights: (19675, 96), symmetry_factors: &[(1, 5), (2, 70), (4, 219), (8, 408), (12, 4), (16, 529), (18, 1), (24, 3), (32, 524), (36, 1), (48, 6), (64, 425), (96, 6), (128, 246), (256, 103), (384, 3), (512, 22), (768, 2), (1024, 12), (2048, 2), (3072, 1)], content: 0x08a76811f8931455, numerators: 0x64b6e52a55ca9139, named: 0xb92ab4861772d577, ids: 0x2b2ea5c92f8a3ccd };
    fast phi4_vac_L2: Row { diagrams: 1, topologies: 1, assignments: 1, inverse_symmetry: (1, 8), weights: (1, 8), symmetry_factors: &[(8, 1)], content: 0xaaf29483b37dfc2e, numerators: 0x9ee3522ae46a98c7, named: 0x9b4091ba687c38c2, ids: 0x5a02b59c625e1892 };
    fast phi4_vac_L3: Row { diagrams: 2, topologies: 2, assignments: 2, inverse_symmetry: (1, 12), weights: (1, 12), symmetry_factors: &[(16, 1), (48, 1)], content: 0x860fea0fe2959078, numerators: 0xf1dbbf34c9c784d2, named: 0xedd9bd8540dd8fa6, ids: 0x702bd41b8ab95dab };
    fast phi4_vac_L4: Row { diagrams: 4, topologies: 4, assignments: 4, inverse_symmetry: (11, 96), weights: (11, 96), symmetry_factors: &[(24, 1), (32, 1), (48, 2)], content: 0x7c9a015c48ce66f1, numerators: 0x89321207c7e09198, named: 0xf603887e4a7d70e5, ids: 0xe52401e168c5a558 };
    fast phi4_vac_L5: Row { diagrams: 10, topologies: 10, assignments: 10, inverse_symmetry: (17, 72), weights: (17, 72), symmetry_factors: &[(16, 1), (32, 3), (48, 2), (64, 1), (128, 2), (144, 1)], content: 0xc84d3a9a80a7e93c, numerators: 0xfb17c8e5e01107ce, named: 0x165d1c5e842112fc, ids: 0x42aadd336581d6f7 };
    fast phi4_vac_L6: Row { diagrams: 28, topologies: 28, assignments: 28, inverse_symmetry: (619, 960), weights: (619, 960), symmetry_factors: &[(16, 4), (32, 6), (48, 3), (64, 4), (96, 4), (120, 1), (128, 2), (144, 1), (288, 1), (320, 2)], content: 0xaa7cbaacc8b958dc, numerators: 0x2712de1422c359eb, named: 0x07bcd38f10af35f6, ids: 0x4f84d1516ff158e7 };
    fast phi4_vac_L7: Row { diagrams: 97, topologies: 97, assignments: 97, inverse_symmetry: (709, 324), weights: (709, 324), symmetry_factors: &[(8, 4), (16, 7), (24, 2), (32, 13), (48, 6), (64, 23), (96, 13), (128, 9), (144, 1), (192, 6), (256, 3), (288, 3), (384, 1), (576, 2), (768, 3), (1296, 1)], content: 0x6876ff8f305ca35f, numerators: 0x72c37493df5cad32, named: 0xef3b02bf618e72a9, ids: 0x4c34024b28a058f6 };
    fast phi34_2pt_L1: Row { diagrams: 3, topologies: 3, assignments: 3, inverse_symmetry: (3, 2), weights: (3, 2), symmetry_factors: &[(2, 3)], content: 0x78c9ff93a7a0c1ed, numerators: 0x954202b142bde800, named: 0x4f0ba38d70b11c67, ids: 0xb42518dd8435c78d };
    fast phi34_2pt_L2: Row { diagrams: 29, topologies: 29, assignments: 29, inverse_symmetry: (25, 3), weights: (25, 3), symmetry_factors: &[(2, 6), (4, 19), (6, 2), (8, 2)], content: 0xcd7bebf282a389a7, numerators: 0x2e563a2edf198f53, named: 0xa4f3a8b11bb6c190, ids: 0x810fc00388d8ca9d };
    fast phi34_2pt_L3: Row { diagrams: 351, topologies: 351, assignments: 351, inverse_symmetry: (1741, 24), weights: (1741, 24), symmetry_factors: &[(1, 2), (2, 36), (4, 125), (6, 3), (8, 141), (12, 22), (16, 21), (48, 1)], content: 0xcad0369a6ad0e609, numerators: 0x9a49e9aaae10c880, named: 0xf73daa83990cad10, ids: 0x937a9cdf6eada069 };
    fast phi34_2pt_L4: Row { diagrams: 5076, topologies: 5076, assignments: 5076, inverse_symmetry: (80299, 96), weights: (80299, 96), symmetry_factors: &[(1, 40), (2, 370), (4, 1180), (6, 13), (8, 1719), (12, 130), (16, 1157), (24, 216), (32, 211), (36, 5), (48, 13), (64, 5), (72, 2), (96, 13), (128, 2)], content: 0xb14dd9e709143bda, numerators: 0x84fab1bbe550b67a, named: 0x6191e43a9ccb41fd, ids: 0xcc6b33da098b9b08 };
    heavy phi34_2pt_L5: Row { diagrams: 84749, topologies: 84749, assignments: 84749, inverse_symmetry: (6869123, 576), weights: (6869123, 576), symmetry_factors: &[(1, 727), (2, 5060), (4, 15053), (6, 145), (8, 23639), (12, 1092), (16, 21069), (24, 2619), (32, 10408), (36, 27), (48, 2146), (64, 2057), (72, 115), (96, 296), (128, 108), (144, 22), (192, 130), (256, 22), (288, 7), (384, 7)], content: 0xad6afc0df417f26b, numerators: 0x3687bef3749d414c, named: 0x2ba64003a44a5415, ids: 0xf00b07bbc0263ebf };
    fast phi3_2pt_L1: Row { diagrams: 1, topologies: 1, assignments: 1, inverse_symmetry: (1, 2), weights: (1, 2), symmetry_factors: &[(2, 1)], content: 0x55e64ff2da4d1e0f, numerators: 0x316859ad4f714f0a, named: 0x284ca460b3985cfe, ids: 0x957a7b87ebcffb36 };
    fast phi3_2pt_tad_L2: Row { diagrams: 10, topologies: 10, assignments: 10, inverse_symmetry: (25, 8), weights: (25, 8), symmetry_factors: &[(2, 3), (4, 6), (8, 1)], content: 0x545153d6589cc2b7, numerators: 0x305722691bab7a7f, named: 0xd9e84d3bf9cee51c, ids: 0xba30d92e49c4b2fa };
    fast phi34_2pt_identical_L3: Row { diagrams: 72, topologies: 72, assignments: 72, inverse_symmetry: (71, 3), weights: (71, 3), symmetry_factors: &[(1, 2), (2, 24), (4, 32), (6, 1), (8, 10), (12, 3)], content: 0x6454c5db351b672c, numerators: 0xe64b2d72fd7c2e4e, named: 0x4a42dbba4affe049, ids: 0xee56a211744f234e };
    fast qed_eemm_L1: Row { diagrams: 17, topologies: 10, assignments: 17, inverse_symmetry: (17, 1), weights: (1, 1), symmetry_factors: &[(1, 17)], content: 0x67015758ee0dfcf5, numerators: 0x21ff907ef41e3790, named: 0x92ff7810a94d1ff7, ids: 0x1dc5ec5f97f654ff };
    fast qed_eemm_L2: Row { diagrams: 357, topologies: 116, assignments: 357, inverse_symmetry: (357, 1), weights: (57, 1), symmetry_factors: &[(1, 357)], content: 0xb41eb8293f39b714, numerators: 0x8459c2e5cff56b07, named: 0x04ccbdf50c929329, ids: 0x566e38b1fd985ea0 };
    heavy qed_eemm_L3: Row { diagrams: 8724, topologies: 1245, assignments: 8724, inverse_symmetry: (8724, 1), weights: (-282, 1), symmetry_factors: &[(1, 8724)], content: 0xb687808b888a00cc, numerators: 0xee51106cc3fe4d1f, named: 0x661d07aebc62dae1, ids: 0xb51f888639fdfb91 };
    fast qed_ee4f_L0: Row { diagrams: 6, topologies: 7, assignments: 6, inverse_symmetry: (6, 1), weights: (-6, 1), symmetry_factors: &[(1, 6)], content: 0x66f2c9719dc08f84, numerators: 0x4efc7631fa9c125f, named: 0x8a27f32a3acc83fd, ids: 0xd5ff3aec66c9a52e };
    fast qed_ee4f_L1: Row { diagrams: 242, topologies: 304, assignments: 242, inverse_symmetry: (242, 1), weights: (10, 1), symmetry_factors: &[(1, 242)], content: 0xbac90b32e6625e27, numerators: 0xb65507056583d4ab, named: 0x76071173c740ae8c, ids: 0x8009b4ab595db4d2 };
    heavy qed_ee4f_L2: Row { diagrams: 8460, topologies: 7909, assignments: 8460, inverse_symmetry: (8460, 1), weights: (828, 1), symmetry_factors: &[(1, 8460)], content: 0x963dbcf6ba7b1c88, numerators: 0x62286f983894e5f0, named: 0x9ffdbc9c8b5864e4, ids: 0xea62a12dfd838cd6 };
    fast qed_eeaaa_L1: Row { diagrams: 384, topologies: 120, assignments: 384, inverse_symmetry: (384, 1), weights: (264, 1), symmetry_factors: &[(1, 384)], content: 0x0fb43c5f515a232f, numerators: 0x044cde5d497d3648, named: 0xd2334763314c3f01, ids: 0x045b832925723c15 };
    heavy qed_eeaaa_L2: Row { diagrams: 16344, topologies: 1581, assignments: 16344, inverse_symmetry: (16344, 1), weights: (-4356, 1), symmetry_factors: &[(1, 16344)], content: 0x98091302c044d5cf, numerators: 0x7a44d325ff935352, named: 0x5d5f16e150dd124d, ids: 0x6939a25e707794a1 };
    fast qed_eeaa_identical_L1: Row { diagrams: 60, topologies: 20, assignments: 66, inverse_symmetry: (60, 1), weights: (42, 1), symmetry_factors: &[(1, 60)], content: 0xc2f14084d3c16139, numerators: 0xd27676fa0c37edb5, named: 0x8e7e16a5e6214aa4, ids: 0xaa53fdeb847997d1 };
    fast qed_eeaa_sign_L1: Row { diagrams: 58, topologies: 20, assignments: 66, inverse_symmetry: (58, 1), weights: (36, 1), symmetry_factors: &[(1, 58)], content: 0xd5ee97232124fe0c, numerators: 0x0193aaeff8c4a870, named: 0xb0aca5b8a7de17e0, ids: 0x1aa25011dcad45d8 };
    fast qed_eeaa_scalar_L1: Row { diagrams: 57, topologies: 20, assignments: 66, inverse_symmetry: (57, 1), weights: (36, 1), symmetry_factors: &[(1, 57)], content: 0x28e8c26fd53fa6bb, numerators: 0xd1995d9a29a96d4f, named: 0xf370e1cb7c88e2c1, ids: 0x6a4ac394ef4f0517 };
    fast qed_ee4f_identical_L1: Row { diagrams: 228, topologies: 304, assignments: 242, inverse_symmetry: (228, 1), weights: (10, 1), symmetry_factors: &[(1, 228)], content: 0xac6a16e86141b488, numerators: 0xedd1b0058d7da78c, named: 0x51770c1b6d8c4003, ids: 0xa6c48dc40eff752b };
    fast qed_ee4f_sign_L1: Row { diagrams: 226, topologies: 304, assignments: 242, inverse_symmetry: (226, 1), weights: (4, 1), symmetry_factors: &[(1, 226)], content: 0x17b789011e6ec656, numerators: 0x83dc420a2d46c1cd, named: 0xae8a66f5b8d4e946, ids: 0x661bd847a204d56a };
    fast qed_ee4f_scalar_L1: Row { diagrams: 225, topologies: 304, assignments: 242, inverse_symmetry: (225, 1), weights: (4, 1), symmetry_factors: &[(1, 225)], content: 0x5ef890dc71487c0a, numerators: 0x599d7c9c1a1dbe29, named: 0x4419b71fdfc0ce0c, ids: 0x434a39f9679bff73 };
    failing qed_eemm_identical_L2: SCALED_TENSOR;
    failing qed_eemm_sign_L2: SCALED_TENSOR;
    failing qed_eemm_scalar_L2: SCALED_TENSOR;
    failing qed_eeaa_identical_L2: SCALED_TENSOR;
    failing qed_eeaa_sign_L2: SCALED_TENSOR;
    failing qed_eeaa_scalar_L2: SCALED_TENSOR;
    fast qed_ase_L3: Row { diagrams: 540, topologies: 10, assignments: 540, inverse_symmetry: (540, 1), weights: (270, 1), symmetry_factors: &[(1, 540)], content: 0xe97f90cc10f0ba33, numerators: 0x8b2ea87aafb96bb5, named: 0xba9713503bfc0b8e, ids: 0x3bb77701c9a055a5 };
    heavy qed_ase_L4: Row { diagrams: 12285, topologies: 64, assignments: 12285, inverse_symmetry: (12285, 1), weights: (189, 1), symmetry_factors: &[(1, 12285)], content: 0xdd973742675250fc, numerators: 0xae1fbd18cac0fd7b, named: 0xc54b216c87ef9aa2, ids: 0xa324ec8a029594ff };
    heavy qed_ase_L5: Row { diagrams: 322704, topologies: 519, assignments: 322704, inverse_symmetry: (322704, 1), weights: (-35478, 1), symmetry_factors: &[(1, 322704)], content: 0x8890f34a21fe1f0b, numerators: 0x49d30cb3c1b10080, named: 0xffefbae154730d8e, ids: 0x4c86ab9e11bbc114 };
    fast qed_ese_L3: Row { diagrams: 172, topologies: 10, assignments: 172, inverse_symmetry: (172, 1), weights: (-10, 1), symmetry_factors: &[(1, 172)], content: 0x200bb021aa8e3fe3, numerators: 0xa6fb39a8583e54f9, named: 0x798f5dfb32643f81, ids: 0x2391d3c49b9fa520 };
    heavy qed_ese_L4: Row { diagrams: 3485, topologies: 64, assignments: 3485, inverse_symmetry: (3485, 1), weights: (-389, 1), symmetry_factors: &[(1, 3485)], content: 0x65cd4c6b3be0e21c, numerators: 0x909048b3b2a5b609, named: 0x597790dff753eb2a, ids: 0x4b5eec491dcd879a };
    fast qcd_gggg_L1: Row { diagrams: 314, topologies: 55, assignments: 314, inverse_symmetry: (294, 1), weights: (-224, 1), symmetry_factors: &[(1, 274), (2, 40)], content: 0x76e5a6e933609cc9, numerators: 0x18e27150fbb3f7ad, named: 0xfba45d64c65d3382, ids: 0x4a8d1e8942f250a3 };
    heavy qcd_gggg_L2: Row { diagrams: 15253, topologies: 840, assignments: 15253, inverse_symmetry: (163229, 12), weights: (24377, 12), symmetry_factors: &[(1, 12086), (2, 2905), (4, 243), (6, 19)], content: 0x858578c792dbbfe3, numerators: 0xd6dbe0288601183a, named: 0xe909babfcee22c19, ids: 0x3ab1751201cb9b68 };
    heavy qcd_gggg_L3: Row { diagrams: 686501, topologies: 14340, assignments: 686501, inverse_symmetry: (4715477, 8), weights: (1169243, 24), symmetry_factors: &[(1, 505145), (2, 157099), (4, 21237), (6, 1512), (8, 1277), (12, 231)], content: 0x68c391d9d10d2457, numerators: 0xabb3ddf64ec13982, named: 0xffa312473774d4a7, ids: 0xfd917d0a7f2024c2 };
    heavy qcd_ggggg_L1: Row { diagrams: 3075, topologies: 562, assignments: 3075, inverse_symmetry: (2895, 1), weights: (-2131, 1), symmetry_factors: &[(1, 2715), (2, 360)], content: 0xcde7da42321b8460, numerators: 0x1b818b81bef6124a, named: 0x59f082315bfecb5d, ids: 0x7ec986d622fdebf8 };
    heavy qcd_ggggg_L2: Row { diagrams: 213390, topologies: 11685, assignments: 213390, inverse_symmetry: (765585, 4), weights: (131385, 4), symmetry_factors: &[(1, 170990), (2, 39280), (4, 2955), (6, 165)], content: 0x49680cb92696c541, numerators: 0x86c858892c455f99, named: 0x9a2dbe6fb9b43b21, ids: 0x371f9c146230e3a7 };
    fast qcd_uugg_L1: Row { diagrams: 99, topologies: 32, assignments: 99, inverse_symmetry: (94, 1), weights: (32, 1), symmetry_factors: &[(1, 89), (2, 10)], content: 0xcae15c82340cc7a1, numerators: 0xf3021c3c83e886f6, named: 0xa754e1c7d3d4a1c3, ids: 0xc1aed2140036cc2f };
    heavy qcd_uugg_L2: Row { diagrams: 3206, topologies: 457, assignments: 3206, inverse_symmetry: (17479, 6), weights: (1925, 6), symmetry_factors: &[(1, 2644), (2, 517), (4, 38), (6, 7)], content: 0xdeb5af9fd1a493e2, numerators: 0xb2cdc545ccb51b8c, named: 0x10c61b6b263630f4, ids: 0xe0978432262d7fb8 };
    heavy qcd_uugg_L3: Row { diagrams: 114301, topologies: 7513, assignments: 114301, inverse_symmetry: (301202, 3), weights: (-17197, 6), symmetry_factors: &[(1, 88202), (2, 22904), (4, 2686), (6, 317), (8, 152), (12, 40)], content: 0x5a4c42972dd54db7, numerators: 0x137263ecc39a8c6d, named: 0x20279d5820fdbda6, ids: 0x84af8834010c7056 };
    fast qcd_uuuu_L2: Row { diagrams: 902, topologies: 230, assignments: 902, inverse_symmetry: (5057, 6), weights: (0, 1), symmetry_factors: &[(1, 788), (2, 106), (4, 6), (6, 2)], content: 0x2ad5d9d5a3670a27, numerators: 0xf3a36ce1c109060a, named: 0xefa80f474fb7b9f4, ids: 0x0274b12c97c843ea };
    fast qcd_gse_L3: Row { diagrams: 1219, topologies: 57, assignments: 1219, inverse_symmetry: (8553, 8), weights: (405, 8), symmetry_factors: &[(1, 940), (2, 241), (4, 31), (6, 1), (8, 5), (12, 1)], content: 0x7dcb189dcb6f5aa4, numerators: 0x48f27e25d1411024, named: 0xe67a4c626946bcc4, ids: 0x2dca47e4ebb20daa };
    heavy qcd_gse_L4: Row { diagrams: 35820, topologies: 666, assignments: 35820, inverse_symmetry: (1485389, 48), weights: (81159, 16), symmetry_factors: &[(1, 26801), (2, 7660), (4, 1140), (6, 83), (8, 112), (12, 14), (16, 9), (24, 1)], content: 0x044b7b0e04b1ce58, numerators: 0xf79a9efed6e20b13, named: 0xa8aa26c9e6d32217, ids: 0xb100e7895382c522 };
    heavy qcd_gse_L5: Row { diagrams: 1267351, topologies: 9892, assignments: 1267351, inverse_symmetry: (308996023, 288), weights: (10378459, 288), symmetry_factors: &[(1, 910237), (2, 298126), (3, 14), (4, 49275), (6, 3114), (8, 5369), (12, 752), (16, 367), (24, 72), (32, 19), (36, 1), (48, 3), (72, 2)], content: 0x5670abb47f3031ca, numerators: 0xb6c8f00b6337680d, named: 0x4c952f116c1dfd32, ids: 0x63bfc0ae9193f7a5 };
    heavy qcd_gse_all_L3: Row { diagrams: 2627, topologies: 72, assignments: 2627, inverse_symmetry: (4427, 2), weights: (1997, 6), symmetry_factors: &[(1, 1871), (2, 626), (4, 102), (6, 15), (8, 10), (12, 3)], content: 0xec89d122df081617, numerators: 0x18c5f56f227468e1, named: 0xe83046e2ba632478, ids: 0xda78248eeacc9097 };
    fast qcd_tad_L2: Row { diagrams: 125, topologies: 6, assignments: 125, inverse_symmetry: (2383, 24), weights: (1207, 24), symmetry_factors: &[(1, 77), (2, 42), (4, 4), (6, 1), (8, 1)], content: 0xe6298357fd676923, numerators: 0xa53bd80db92a26d0, named: 0x55a86640b9eabf45, ids: 0x8f75a2e8cea1cdff };
    heavy qcd_tad_L3: Row { diagrams: 3490, topologies: 46, assignments: 3490, inverse_symmetry: (124199, 48), weights: (-30613, 48), symmetry_factors: &[(1, 1883), (2, 1247), (3, 7), (4, 271), (6, 29), (8, 43), (12, 5), (16, 4), (48, 1)], content: 0x6b4b08a9e6607c1b, numerators: 0x14089c877c23f095, named: 0x1a847496091e4fad, ids: 0x0ae5b9753a69c8f2 };
    fast qcd_vac_L2: Row { diagrams: 45, topologies: 3, assignments: 45, inverse_symmetry: (191, 6), weights: (107, 6), symmetry_factors: &[(1, 21), (2, 21), (8, 2), (12, 1)], content: 0x0fdd50aef970b13d, numerators: 0x6173058343f36369, named: 0x895ea4bbc0107f71, ids: 0xf3bf59ff814c96a1 };
    fast qcd_vac_L3: Row { diagrams: 649, topologies: 12, assignments: 649, inverse_symmetry: (2417, 6), weights: (-107, 1), symmetry_factors: &[(1, 231), (2, 287), (3, 7), (4, 84), (6, 14), (8, 18), (12, 1), (16, 4), (24, 1), (48, 2)], content: 0x145ae3379080724d, numerators: 0x34b5417dee2213d7, named: 0x3ebce22f2ae93d29, ids: 0xe90d219aec38cec4 };
    heavy qcd_vac_L4: Row { diagrams: 14017, topologies: 73, assignments: 14017, inverse_symmetry: (2597681, 288), weights: (284153, 288), symmetry_factors: &[(1, 5698), (2, 5383), (3, 154), (4, 1894), (6, 189), (8, 483), (12, 79), (16, 85), (24, 21), (32, 14), (48, 13), (72, 2), (128, 1), (384, 1)], content: 0xc54f17048580de58, numerators: 0xf202887af896ed26, named: 0x725a44a42dce5763, ids: 0x7213d1d60be6ae04 };
    fast qcd_vac1pi_L4: Row { diagrams: 628, topologies: 26, assignments: 628, inverse_symmetry: (11465, 36), weights: (446, 9), symmetry_factors: &[(1, 154), (2, 217), (3, 56), (4, 101), (6, 28), (8, 48), (12, 8), (16, 8), (24, 2), (32, 2), (48, 3), (72, 1)], content: 0x7966c26e2b7b8717, numerators: 0xfad36d65e1a809b0, named: 0x50b86726459b5bb6, ids: 0x338c3e1aa8433d68 };
    heavy qcd_vac1pi_L5: Row { diagrams: 9889, topologies: 166, assignments: 9889, inverse_symmetry: (7313267, 1152), weights: (227027, 1152), symmetry_factors: &[(1, 4165), (2, 3455), (4, 1494), (6, 49), (8, 472), (12, 46), (16, 159), (24, 11), (32, 22), (48, 6), (64, 4), (96, 2), (128, 3), (144, 1)], content: 0xf0295cc039fe2a58, numerators: 0x30a89b2dbb1dd465, named: 0x27a4758ab5580b0c, ids: 0xc282afcf907784ca };
    fast qcd_gggg_sym_L1: Row { diagrams: 123, topologies: 25, assignments: 123, inverse_symmetry: (147, 2), weights: (-56, 1), symmetry_factors: &[(1, 37), (2, 62), (4, 20), (8, 4)], content: 0x3ee805c0752694b4, numerators: 0xd9d2fce0e9012de3, named: 0x0c06dd4326fe31ba, ids: 0x0bd7306b15a54577 };
    failing qcd_gggg_identical_L1: COLOUR_DIMENSION;
    failing qcd_gggg_sign_L1: COLOUR_DIMENSION;
    failing qcd_gggg_scalar_L1: COLOUR_DIMENSION;
    fast qcd_gggg_ghost_filter_L1: Row { diagrams: 277, topologies: 55, assignments: 277, inverse_symmetry: (257, 1), weights: (-187, 1), symmetry_factors: &[(1, 237), (2, 40)], content: 0x90a1ce06e025a6f5, numerators: 0x66554bce91f5d5d7, named: 0x1e9b4de5869edb85, ids: 0xda3302c00009cd6d };
    fast qcd_gggg_ghost_veto_L1: Row { diagrams: 277, topologies: 55, assignments: 277, inverse_symmetry: (257, 1), weights: (-187, 1), symmetry_factors: &[(1, 237), (2, 40)], content: 0x90a1ce06e025a6f5, numerators: 0x66554bce91f5d5d7, named: 0x1e9b4de5869edb85, ids: 0xda3302c00009cd6d };
    fast ym_g6_L0: Row { diagrams: 220, topologies: 220, assignments: 220, inverse_symmetry: (220, 1), weights: (220, 1), symmetry_factors: &[(1, 220)], content: 0xe166545bc76a8cd8, numerators: 0x9f37ff3f3deffb74, named: 0x0855fe3d5967a4de, ids: 0x4a0b7e227d7ace53 };
    heavy ym_g7_L0: Row { diagrams: 2485, topologies: 2485, assignments: 2485, inverse_symmetry: (2485, 1), weights: (2485, 1), symmetry_factors: &[(1, 2485)], content: 0x057203fcdeaf182e, numerators: 0x518deeeb2941c74c, named: 0x825ef5b96a66dc8d, ids: 0x54cd7d6e49c55438 };
    heavy ym_g8_L0: Row { diagrams: 34300, topologies: 34300, assignments: 34300, inverse_symmetry: (34300, 1), weights: (34300, 1), symmetry_factors: &[(1, 34300)], content: 0x0d8cd8312f854d8b, numerators: 0xe586bec340cace76, named: 0xe921e3ce642cddae, ids: 0xa9980ccebbd89562 };
    heavy ym_gse_L4: Row { diagrams: 2571, topologies: 666, assignments: 2571, inverse_symmetry: (81653, 48), weights: (-8113, 16), symmetry_factors: &[(1, 1139), (2, 889), (4, 402), (6, 17), (8, 100), (12, 14), (16, 9), (24, 1)], content: 0xbb8866a313b845e9, numerators: 0x3ec48d93c61cce4c, named: 0xfce8c790d0d23d2d, ids: 0xabd3694977e1dd34 };
    heavy ym_gse_L5: Row { diagrams: 53201, topologies: 9892, assignments: 53201, inverse_symmetry: (9957739, 288), weights: (-2104793, 288), symmetry_factors: &[(1, 22781), (2, 18577), (3, 2), (4, 8421), (6, 315), (8, 2411), (12, 248), (16, 349), (24, 72), (32, 19), (36, 1), (48, 3), (72, 2)], content: 0x64121e08ff2ba484, numerators: 0x204eb02540fdbc12, named: 0xc2f24b12507a1efd, ids: 0x3f065e6318f32403 };
    fast ym_vac_L5: Row { diagrams: 381, topologies: 91, assignments: 381, inverse_symmetry: (95437, 576), weights: (-19595, 576), symmetry_factors: &[(1, 82), (2, 110), (4, 73), (6, 4), (8, 54), (12, 9), (16, 28), (24, 3), (32, 8), (48, 4), (64, 2), (96, 1), (128, 2), (144, 1)], content: 0x8338e2f92f7cfd9d, numerators: 0xd4db3703b3a0e21b, named: 0xbaa296b9b79eacc5, ids: 0xe0ce825b726cc352 };
    fast sm_eemm_L0: Row { diagrams: 4, topologies: 1, assignments: 4, inverse_symmetry: (4, 1), weights: (-4, 1), symmetry_factors: &[(1, 4)], content: 0x97af00ae91ef62c4, numerators: 0x7c4bccf1c63f9180, named: 0xc01c07a8158bc53a, ids: 0x91a782bcdd47cb2a };
    fast sm_eemm_L1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), weights: (-277, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0x99843c92c72133b6, numerators: 0x26f3641db6659cf2, named: 0xd3e09a8c64a70ce6, ids: 0xdd8840070bc24da1 };
    heavy sm_eemm_L2: Row { diagrams: 69697, topologies: 160, assignments: 69697, inverse_symmetry: (415391, 6), weights: (-6677, 6), symmetry_factors: &[(1, 68781), (2, 888), (4, 26), (6, 2)], content: 0xb4b43e66a248f228, numerators: 0x5a12870245dd2b8c, named: 0x3f5168c221367ed0, ids: 0xe12c04b47fd707fe };
    failing sm_eemm_identical_L1: NON_SELF_DUAL_POWER;
    failing sm_eemm_sign_L1: NON_SELF_DUAL_POWER;
    failing sm_eemm_scalar_L1: NON_SELF_DUAL_POWER;
    fast sm_ee4f_L0: Row { diagrams: 18, topologies: 7, assignments: 18, inverse_symmetry: (18, 1), weights: (-18, 1), symmetry_factors: &[(1, 18)], content: 0x7dab707954851b41, numerators: 0x71eb659079856e7a, named: 0xf3ff9500be4fb50c, ids: 0x5b5490ac1128d4fa };
    heavy sm_ee4f_L1: Row { diagrams: 5811, topologies: 322, assignments: 5811, inverse_symmetry: (5805, 1), weights: (-2685, 1), symmetry_factors: &[(1, 5799), (2, 12)], content: 0xaea64b8a548b0cf6, numerators: 0xe557f43bb1948d11, named: 0x2cd22700423a72e1, ids: 0xd66974b61312a285 };
    fast sm_eeww_L0: Row { diagrams: 4, topologies: 2, assignments: 4, inverse_symmetry: (4, 1), weights: (-4, 1), symmetry_factors: &[(1, 4)], content: 0x7477c6f9f9a39d73, numerators: 0x63db443730cf91ef, named: 0x48249b27b144364c, ids: 0xf65454f1b8b1a032 };
    fast sm_eeww_L1: Row { diagrams: 919, topologies: 23, assignments: 919, inverse_symmetry: (916, 1), weights: (-2, 1), symmetry_factors: &[(1, 913), (2, 6)], content: 0x953e54839edd5025, numerators: 0x077e524238dbfee5, named: 0x19cff4ac923cbb6c, ids: 0x60a26603fe037315 };
    fast sm_ggtt_L1: Row { diagrams: 375, topologies: 32, assignments: 375, inverse_symmetry: (370, 1), weights: (-8, 1), symmetry_factors: &[(1, 365), (2, 10)], content: 0xb57aa9d432880775, numerators: 0x6c6e9542edb615e9, named: 0x0128d204403bdd9a, ids: 0x7c92666038a148c0 };
    fast sm_ggttg_L0: Row { diagrams: 16, topologies: 16, assignments: 16, inverse_symmetry: (16, 1), weights: (16, 1), symmetry_factors: &[(1, 16)], content: 0xeb9e3621b43c55e4, numerators: 0x59b38f69f6396929, named: 0x66f155052cc8e9b4, ids: 0xfaacbb5903d5fcab };
    heavy sm_ggttg_L1: Row { diagrams: 3096, topologies: 322, assignments: 3096, inverse_symmetry: (3046, 1), weights: (-316, 1), symmetry_factors: &[(1, 2996), (2, 100)], content: 0x60920f521092c4ea, numerators: 0x96c505e3b065e987, named: 0xf724c023a54f2c22, ids: 0x8ca1cbee96df2d18 };
    fast sm_hse_L1: Row { diagrams: 17, topologies: 1, assignments: 17, inverse_symmetry: (31, 2), weights: (-5, 2), symmetry_factors: &[(1, 14), (2, 3)], content: 0x825faa38b23aed38, numerators: 0xdf018e8ddee914be, named: 0xd0ebfdcbaa07b1bf, ids: 0xa15a635f45901d97 };
    fast sm_hse_L2: Row { diagrams: 1188, topologies: 7, assignments: 1188, inverse_symmetry: (13691, 12), weights: (-4417, 12), symmetry_factors: &[(1, 1100), (2, 76), (4, 11), (6, 1)], content: 0x1c50eecb5249aa7d, numerators: 0xe7757c1488aeba4b, named: 0x446b98d9956a01bb, ids: 0xaf6c2e5f9d7199f1 };
    fast sm_ase_L2: Row { diagrams: 904, topologies: 7, assignments: 904, inverse_symmetry: (904, 1), weights: (-424, 1), symmetry_factors: &[(1, 904)], content: 0xe45d4a5e047a542c, numerators: 0x5b5f703fd30dfaa8, named: 0x39bd39bfa4d58bde, ids: 0x55b617bfcdda6709 };
    fast sm_eemm_qcd0_L1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), weights: (-277, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0x99843c92c72133b6, numerators: 0x26f3641db6659cf2, named: 0xd3e09a8c64a70ce6, ids: 0xdd8840070bc24da1 };
    heavy sm_eemm_qcd0_L2: Row { diagrams: 69193, topologies: 160, assignments: 69193, inverse_symmetry: (412367, 6), weights: (-5813, 6), symmetry_factors: &[(1, 68277), (2, 888), (4, 26), (6, 2)], content: 0xe342a202fff52526, numerators: 0x3a0eab7af9b8b1a9, named: 0xe01a8f1f431abb1e, ids: 0xc17356dba6fed0ed };
    fast sm_uudd_qed2_L1: Row { diagrams: 111, topologies: 17, assignments: 111, inverse_symmetry: (221, 2), weights: (-57, 2), symmetry_factors: &[(1, 110), (2, 1)], content: 0x422125ee75cb81f1, numerators: 0x0e06bd22fdc90a2f, named: 0xb9d3140580cff874, ids: 0x72b7c0aea620df67 };
    heavy sm_uudd_qed2_L2: Row { diagrams: 3286, topologies: 230, assignments: 3286, inverse_symmetry: (38519, 12), weights: (733, 12), symmetry_factors: &[(1, 3136), (2, 146), (4, 3), (6, 1)], content: 0xa681cec8fde8b382, numerators: 0x815c7654c050793d, named: 0xdf3098711d2e0304, ids: 0x8758ad2d4e3187f0 };
    fast ew_eemm_L1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), weights: (-277, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0xb56d549c918bbbc1, numerators: 0x2bc439d340929b34, named: 0x8dcc068955106981, ids: 0x19249f55d84eae49 };
    heavy ew_eemm_L2: Row { diagrams: 69193, topologies: 160, assignments: 69193, inverse_symmetry: (412367, 6), weights: (-5813, 6), symmetry_factors: &[(1, 68277), (2, 888), (4, 26), (6, 2)], content: 0x270481b3cce8b479, numerators: 0x6292f42e7e2fb42d, named: 0x26ee542dae6e5ef4, ids: 0x12ffd9ba3fb86474 };
    fast ct_phi3_2pt_L1: Row { diagrams: 4, topologies: 1, assignments: 4, inverse_symmetry: (2, 1), weights: (2, 1), symmetry_factors: &[(2, 4)], content: 0x792c2689f9b727df, numerators: 0x094dd49c998c2915, named: 0x3217eb0128aeea77, ids: 0x247d40d238f195f1 };
    fast ct_qcd_gggg_L1: Row { diagrams: 1607, topologies: 55, assignments: 1607, inverse_symmetry: (2883, 2), weights: (-1549, 2), symmetry_factors: &[(1, 1276), (2, 331)], content: 0x102a6a28e1c425f5, numerators: 0x66ea9c29e5d8a21d, named: 0x1ab462bc526da632, ids: 0x626e72bf1b28cdfa };
    fast ct_qcd_uugg_L1: Row { diagrams: 864, topologies: 32, assignments: 864, inverse_symmetry: (802, 1), weights: (-130, 1), symmetry_factors: &[(1, 740), (2, 124)], content: 0xbdb6f0a3f62831ef, numerators: 0xab75b6346178447f, named: 0x4c71738657302058, ids: 0xc8abf0cff4acf61b };
    fast ct_phi3_2pt_tad_L2: Row { diagrams: 148, topologies: 10, assignments: 148, inverse_symmetry: (50, 1), weights: (50, 1), symmetry_factors: &[(1, 4), (2, 44), (4, 92), (8, 8)], content: 0xf5122b0d42fad77e, numerators: 0xb1c2ba48ec6c4558, named: 0xa97896825197b385, ids: 0x6cdf274d98c848ce };
    fast ct_phi3_vac_L2: Row { diagrams: 6, topologies: 2, assignments: 6, inverse_symmetry: (5, 6), weights: (5, 6), symmetry_factors: &[(4, 1), (6, 1), (8, 2), (12, 2)], content: 0x9dad3c92417664f1, numerators: 0xf57c18fe6786fb79, named: 0x12031bbe318c7d29, ids: 0x9e59d4664215d958 };
    fast ct_phi3_vac_L3: Row { diagrams: 42, topologies: 5, assignments: 42, inverse_symmetry: (5, 1), weights: (5, 1), symmetry_factors: &[(4, 7), (6, 2), (8, 17), (16, 10), (24, 2), (48, 4)], content: 0xb5863cd09c348f6c, numerators: 0xcadc67ebf2b44f13, named: 0xc86f79cbcf475f9e, ids: 0xceb477ab49af8581 };
    fast ct_phi3_vac_L4: Row { diagrams: 505, topologies: 17, assignments: 505, inverse_symmetry: (1105, 18), weights: (1105, 18), symmetry_factors: &[(1, 1), (2, 17), (4, 76), (6, 1), (8, 161), (12, 6), (16, 153), (24, 1), (32, 70), (36, 1), (48, 6), (64, 6), (72, 2), (128, 4)], content: 0xdac00a3ffc146b39, numerators: 0x758195af180c3984, named: 0x6e909ab232f48905, ids: 0x9ddcd18b5b63a980 };
    fast ct_qcd_vac_L2: Row { diagrams: 68, topologies: 3, assignments: 68, inverse_symmetry: (1103, 24), weights: (479, 24), symmetry_factors: &[(1, 29), (2, 32), (4, 1), (6, 1), (8, 3), (12, 2)], content: 0xf69da97cd4bfb079, numerators: 0x74cb264b3044f795, named: 0x382ae441e068e697, ids: 0xbc728af284e8f8c0 };
    fast ct_qcd_vac_L3: Row { diagrams: 1816, topologies: 12, assignments: 1816, inverse_symmetry: (6341, 6), weights: (-321, 2), symmetry_factors: &[(1, 549), (2, 826), (3, 16), (4, 293), (6, 34), (8, 74), (12, 4), (16, 13), (24, 2), (48, 5)], content: 0x78eccaa3c84acd2a, numerators: 0xca96a9c3deeb724f, named: 0x69b0bd134250cfe6, ids: 0x8d11ffb4cb83725b };
    fast xs_qed_eemm_L1: Row { diagrams: 1, topologies: 17, assignments: 34, inverse_symmetry: (1, 1), weights: (-1, 1), symmetry_factors: &[(1, 1)], content: 0x831ae534b022fd0f, numerators: 0xe134f67914689408, named: 0xf40bfdc791015105, ids: 0xedd739fe2aef9d55 };
    fast xs_qed_eemm_L2: Row { diagrams: 30, topologies: 157, assignments: 714, inverse_symmetry: (30, 1), weights: (4, 1), symmetry_factors: &[(1, 30)], content: 0x3ce6c5b6af8ba125, numerators: 0xda1650e44404e0fc, named: 0x0b41c948729c87ec, ids: 0xc6537de685c0ec01 };
    fast xs_qed_eemm_L3: Row { diagrams: 824, topologies: 1536, assignments: 17448, inverse_symmetry: (824, 1), weights: (94, 1), symmetry_factors: &[(1, 824)], content: 0x362776deabca39da, numerators: 0x943ab33171ec7737, named: 0xb5205b5cda3bf8ff, ids: 0xc3905035e4cb6566 };
    fast xs_qed_eemm_identical_L2: Row { diagrams: 28, topologies: 157, assignments: 714, inverse_symmetry: (28, 1), weights: (4, 1), symmetry_factors: &[(1, 28)], content: 0x44c6a2185761acd8, numerators: 0x09c7fc2ec0007936, named: 0x3dd6bf8a6993a599, ids: 0x92036d905b3178ba };
    fast xs_qed_eemm_sign_L2: Row { diagrams: 28, topologies: 157, assignments: 714, inverse_symmetry: (28, 1), weights: (4, 1), symmetry_factors: &[(1, 28)], content: 0x44c6a2185761acd8, numerators: 0x09c7fc2ec0007936, named: 0x3dd6bf8a6993a599, ids: 0x92036d905b3178ba };
    fast xs_qed_eemm_scalar_L2: Row { diagrams: 28, topologies: 157, assignments: 714, inverse_symmetry: (28, 1), weights: (4, 1), symmetry_factors: &[(1, 28)], content: 0x44c6a2185761acd8, numerators: 0x09c7fc2ec0007936, named: 0x3dd6bf8a6993a599, ids: 0x92036d905b3178ba };
    fast xs_qed_eedd_L1: Row { diagrams: 1, topologies: 17, assignments: 34, inverse_symmetry: (1, 1), weights: (-1, 1), symmetry_factors: &[(1, 1)], content: 0xb467d6be7a12aefa, numerators: 0x52ae30868db3c38d, named: 0xdd30b2e206f30267, ids: 0x58ebcb07ba6255d3 };
    fast xs_qed_eedd_L2: Row { diagrams: 30, topologies: 157, assignments: 714, inverse_symmetry: (30, 1), weights: (4, 1), symmetry_factors: &[(1, 30)], content: 0x1a55cd5c769f869a, numerators: 0x649671e7e59d982c, named: 0x368a1cfe1b96c3be, ids: 0xfba0be1b20d86887 };
    fast xs_qed_eeaa_L2: Row { diagrams: 96, topologies: 157, assignments: 714, inverse_symmetry: (96, 1), weights: (-48, 1), symmetry_factors: &[(1, 96)], content: 0xb38bf88ee71b032b, numerators: 0x8c6fd49a5d2b9331, named: 0xe67ffcbf30548e8b, ids: 0x9d68f748335a7cfd };
    fast xs_qed_aaee_L2: Row { diagrams: 248, topologies: 75, assignments: 3456, inverse_symmetry: (248, 1), weights: (176, 1), symmetry_factors: &[(1, 248)], content: 0x35f49d2808ebf1b2, numerators: 0xa10a1ce9c9e640dd, named: 0x26940e2744fc7077, ids: 0x0ccc891930b8d526 };
    fast xs_qed_aaee_sym_L2: Row { diagrams: 124, topologies: 75, assignments: 3456, inverse_symmetry: (124, 1), weights: (176, 1), symmetry_factors: &[(1, 124)], content: 0x3ee1a214e4975eb6, numerators: 0xdd11b49adff50d91, named: 0xb0209cd25cac901b, ids: 0x91e367b69f3f8408 };
    fast xs_qed_aaee_lr_L2: Row { diagrams: 130, topologies: 75, assignments: 3456, inverse_symmetry: (130, 1), weights: (176, 1), symmetry_factors: &[(1, 130)], content: 0x326238c5cdd3bd7f, numerators: 0x7ffd06aad7022532, named: 0x6d6b618cf6a933a2, ids: 0x812de9ca53130bb4 };
    fast xs_qed_aaee_sym_lr_L2: Row { diagrams: 82, topologies: 75, assignments: 3456, inverse_symmetry: (82, 1), weights: (176, 1), symmetry_factors: &[(1, 82)], content: 0xe7473d645fa53156, numerators: 0xb7c50e7cc0dd70ab, named: 0x6839af52dd4f3be5, ids: 0x629bb7150d8c739d };
    fast xs_qed_eeee_L1: Row { diagrams: 4, topologies: 17, assignments: 34, inverse_symmetry: (4, 1), weights: (0, 1), symmetry_factors: &[(1, 4)], content: 0x29f20710d0fbfbd6, numerators: 0x099879207707ff29, named: 0xba552fc69e567351, ids: 0xfe0da547942d4483 };
    fast xs_qed_eeee_L2: Row { diagrams: 120, topologies: 157, assignments: 714, inverse_symmetry: (120, 1), weights: (0, 1), symmetry_factors: &[(1, 120)], content: 0x160b83f7a26726b5, numerators: 0x27dc41969a149d60, named: 0xfd42901c0c72c5c4, ids: 0x3cd3ebae4a2c6d84 };
    fast xs_qed_eeee_1pi_L1: Row { diagrams: 1, topologies: 3, assignments: 4, inverse_symmetry: (1, 1), weights: (-1, 1), symmetry_factors: &[(1, 1)], content: 0x8686eb17d2b0f293, numerators: 0xd82920bd552e843e, named: 0xf5ac9fd01b9eabd7, ids: 0xf777e735b39d93bf };
    fast xs_qed_eeee_1pi_L2: Row { diagrams: 31, topologies: 36, assignments: 116, inverse_symmetry: (31, 1), weights: (13, 1), symmetry_factors: &[(1, 31)], content: 0x57499dc520b5c855, numerators: 0xa51417d5ebc4fa1f, named: 0x5d7b34e127921dbe, ids: 0x35c2245f6a859c80 };
    fast xs_sm_eemm_L1: Row { diagrams: 16, topologies: 17, assignments: 856, inverse_symmetry: (16, 1), weights: (-16, 1), symmetry_factors: &[(1, 16)], content: 0x9ad7451811ae86e7, numerators: 0x01b59b2f2dfb564d, named: 0xe581437c7ccd8710, ids: 0x27354d61e6800e7c };
    fast xs_sm_eemm_L2: Row { diagrams: 3104, topologies: 230, assignments: 139394, inverse_symmetry: (3092, 1), weights: (-916, 1), symmetry_factors: &[(1, 3080), (2, 24)], content: 0xfa17e9037ca769f3, numerators: 0x3915a46d9349eac3, named: 0x18fd49c31b9c04e3, ids: 0x421c3c7d3a525aa1 };
    fast xs_sm_uudd_L1: Row { diagrams: 16, topologies: 17, assignments: 558, inverse_symmetry: (16, 1), weights: (-4, 1), symmetry_factors: &[(1, 16)], content: 0x6a0f6683f376c99d, numerators: 0x7d74ed995684f1db, named: 0x46b5bfdc898f0311, ids: 0xc4c4466a8cca1e9e };
    fast xs_qcd_gggg_L1: Row { diagrams: 10, topologies: 55, assignments: 314, inverse_symmetry: (8, 1), weights: (8, 1), symmetry_factors: &[(1, 6), (2, 4)], content: 0xf64ad313fb15b61e, numerators: 0xbde2aaa4f3b826c0, named: 0x3d5aa5fb16e7e40a, ids: 0xe2e7018b7e154d23 };
    fast xs_qcd_gggg_sym_L1: Row { diagrams: 10, topologies: 55, assignments: 314, inverse_symmetry: (8, 1), weights: (8, 1), symmetry_factors: &[(1, 6), (2, 4)], content: 0xf64ad313fb15b61e, numerators: 0xbde2aaa4f3b826c0, named: 0x3d5aa5fb16e7e40a, ids: 0xe2e7018b7e154d23 };
    fast xs_qcd_gggg_L2: Row { diagrams: 1228, topologies: 840, assignments: 15253, inverse_symmetry: (1040, 1), weights: (-808, 1), symmetry_factors: &[(1, 868), (2, 328), (4, 32)], content: 0xfcf9d3749673e682, numerators: 0x258ac63d13945e30, named: 0xa26d5fc68709e8b7, ids: 0x1ff67684a166f140 };
    fast xs_qcd_gguu_L1: Row { diagrams: 9, topologies: 55, assignments: 314, inverse_symmetry: (9, 1), weights: (-9, 1), symmetry_factors: &[(1, 9)], content: 0xf55ef4afe4fbeee8, numerators: 0x9e713702e766065d, named: 0xeadb4aa521f5d4c9, ids: 0x2ddf4705c9d52678 };
    fast xs_qcd_gguu_L2: Row { diagrams: 558, topologies: 840, assignments: 15253, inverse_symmetry: (528, 1), weights: (210, 1), symmetry_factors: &[(1, 498), (2, 60)], content: 0x455c66127358d46b, numerators: 0x0a6323254f20878e, named: 0x8349d5afc24bc7d9, ids: 0x9ae814b9b6a488d9 };
    fast xs_qcd_gguu_sym_L1: Row { diagrams: 5, topologies: 55, assignments: 314, inverse_symmetry: (5, 1), weights: (-9, 1), symmetry_factors: &[(1, 5)], content: 0xb8e35da815144186, numerators: 0xa9e3dfa7c6a49ca8, named: 0xa1706e5ab99b8939, ids: 0xde71f8568766e9f4 };
    fast xs_qcd_gguu_1pi_L1: Row { diagrams: 4, topologies: 12, assignments: 54, inverse_symmetry: (4, 1), weights: (-4, 1), symmetry_factors: &[(1, 4)], content: 0x4ee6aab30b16f1c5, numerators: 0x901412cc9fc683f1, named: 0xd40765b02a1bdf5a, ids: 0x822a3b576554f53b };
    fast sm_uudd_qcd2_qed2_L1: Row { diagrams: 93, topologies: 17, assignments: 93, inverse_symmetry: (93, 1), weights: (-25, 1), symmetry_factors: &[(1, 93)], content: 0x7b412ce40ab68cb6, numerators: 0xc13ed6023bd9fa73, named: 0x02c2096296b78f9d, ids: 0xc06602607aa85f37 };
    fast sm_uudd_qcd3_L1: Row { diagrams: 18, topologies: 17, assignments: 18, inverse_symmetry: (35, 2), weights: (-7, 2), symmetry_factors: &[(1, 17), (2, 1)], content: 0x17f0a12f131f0bca, numerators: 0xa3ec8b01207c1a45, named: 0x75d501c4462bb8d6, ids: 0x6540c51587a5e3b1 };
    fast sm_eemm_qed4_L0to1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), weights: (-277, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0x99843c92c72133b6, numerators: 0x26f3641db6659cf2, named: 0xd3e09a8c64a70ce6, ids: 0xdd8840070bc24da1 };
    fast sm_eemm_floop0_L1: Row { diagrams: 317, topologies: 10, assignments: 317, inverse_symmetry: (631, 2), weights: (-499, 2), symmetry_factors: &[(1, 314), (2, 3)], content: 0xd23895671b6f9dff, numerators: 0x517870f1cae9a209, named: 0x8bd58d385ca26fcc, ids: 0x9daa1bb29c9c73ae };
    fast qcd_gggg_floop1_L1: Row { diagrams: 222, topologies: 55, assignments: 222, inverse_symmetry: (222, 1), weights: (-222, 1), symmetry_factors: &[(1, 222)], content: 0xb72175978408ff0d, numerators: 0xc642d07c1f80bb4c, named: 0x978d5636f7e2a170, ids: 0x0c9a9ab09256e03e };
    fast qcd_gse_all_floop1_dfilter_L2: Row { diagrams: 55, topologies: 8, assignments: 55, inverse_symmetry: (95, 2), weights: (-55, 2), symmetry_factors: &[(1, 40), (2, 15)], content: 0x6a14158e96795e50, numerators: 0x4453dad44443bf23, named: 0xd35a514c8b8b41d1, ids: 0xbd4606b364c03f53 };
    fast sm_eemm_floop1_tfilter_L1: Row { diagrams: 95, topologies: 10, assignments: 95, inverse_symmetry: (95, 1), weights: (95, 1), symmetry_factors: &[(1, 95)], content: 0xdcd2cfd04f896b17, numerators: 0x5b9d199cf88e1fc2, named: 0x8df020f0055500ee, ids: 0xe2745936ddf536ec };
    fast qed_ase_floop2_L3: Row { diagrams: 405, topologies: 10, assignments: 405, inverse_symmetry: (405, 1), weights: (405, 1), symmetry_factors: &[(1, 405)], content: 0x184a252bac2dfc3a, numerators: 0xdc2e4020a7a0a27a, named: 0x109b995d1525e3de, ids: 0x6daf1e85e33ce62c };
    fast ctorder_qcd_gggg_L1: Row { diagrams: 544, topologies: 55, assignments: 544, inverse_symmetry: (492, 1), weights: (-292, 1), symmetry_factors: &[(1, 440), (2, 104)], content: 0xf8d9d872e3f3e2d5, numerators: 0x7c69f3b81e40b745, named: 0x4de45c4aa659d10e, ids: 0x6dd3cc5ef7c80e73 };
    fast ctorder_qcd_uugg_L1: Row { diagrams: 268, topologies: 32, assignments: 268, inverse_symmetry: (251, 1), weights: (13, 1), symmetry_factors: &[(1, 234), (2, 34)], content: 0x4c1ad774fc27e60a, numerators: 0x555d26de4a5d2b3e, named: 0xe620415a5190f3d2, ids: 0xfe08e4538f1fc920 };
    fast ctorder_qcd_vac_L3: Row { diagrams: 1298, topologies: 12, assignments: 1298, inverse_symmetry: (18965, 24), weights: (-3859, 24), symmetry_factors: &[(1, 431), (2, 601), (3, 14), (4, 177), (6, 30), (8, 33), (12, 3), (16, 5), (24, 1), (48, 3)], content: 0xbe06185c11dbafc3, numerators: 0x3490357e136ed1bd, named: 0xde87a43c14232163, ids: 0xe8035567cdb53120 };
    fast ctorder_phi3_vac_L4: Row { diagrams: 117, topologies: 17, assignments: 117, inverse_symmetry: (5525, 384), weights: (5525, 384), symmetry_factors: &[(2, 5), (4, 17), (8, 37), (12, 1), (16, 38), (32, 17), (64, 1), (128, 1)], content: 0xcfc6a8d4648abdc1, numerators: 0x1ceef6c80c321a4f, named: 0xbf7944ef5ddb5861, ids: 0xcfc1f9dd64584517 };
}

#[test]
fn oracle_rows_cover_every_case() {
    let cases: BTreeSet<_> = CASES.iter().map(|case| case.name).collect();
    let rows: BTreeSet<_> = ORACLE.iter().copied().collect();
    assert_eq!(cases, rows);
    assert_eq!(rows.len(), ORACLE.len(), "duplicate oracle rows");
    assert!(REFERENCE.iter().all(|name| cases.contains(name)));
}

fn is_ghost(model: &Model, diagram: &FeynmanDiagram) -> bool {
    diagram
        .edges()
        .any(|(_, _, edge)| model.particle_by_id(edge.particle).unwrap().is_ghost())
}

/// Diagram ids depend on the diagrams alone: fresh processes give the same ids at 1 thread and at
/// 8 threads after building other models first, which makes Symbolica register the symbols of the
/// cases' models in another order. The test runs its own binary as those processes.
#[test]
fn ids_do_not_depend_on_threads_or_symbol_registration() {
    const CHILD: &str = "FEYNKIT_ORACLE_ID_CHILD";
    const RUN: [&str; 3] = ["qcd_uugg_L1", "sm_eemm_L1", "xs_qed_eemm_L2"];
    if let Ok(child) = std::env::var(CHILD) {
        let (threads, preload) = match child.split_once(',') {
            Some((threads, "preload")) => (threads, true),
            _ => (child.as_str(), false),
        };
        if preload {
            for model in [
                Model::phi_3_4,
                Model::yang_mills,
                Model::electroweak,
                Model::qcd_qed,
            ] {
                model();
            }
        }
        for name in RUN {
            let case = Case::named(name);
            let model = case.model();
            let options = case.options(&model, threads.parse().unwrap());
            for diagram in case.generate_with(model, &options).unwrap().diagrams {
                println!(
                    "{CHILD} {name} {} {}",
                    diagram.name().expect("named oracle output"),
                    diagram.id().unwrap()
                );
            }
        }
        return;
    }
    let ids = |child: &str| -> Vec<String> {
        let test = "ids_do_not_depend_on_threads_or_symbol_registration";
        let output = Command::new(std::env::current_exe().unwrap())
            .args(["--exact", test, "--nocapture"])
            .env(CHILD, child)
            .output()
            .unwrap();
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(
            output.status.success(),
            "{child}: {stdout}{}",
            String::from_utf8_lossy(&output.stderr)
        );
        stdout
            .lines()
            .filter(|line| line.starts_with(CHILD))
            .map(str::to_owned)
            .collect()
    };
    let reference = ids("1");
    assert!(reference.len() > 500);
    assert_eq!(ids("8,preload"), reference);
}

/// A user filter pruning in flight, a particle veto and post-filtering the unfiltered diagrams
/// select the same ghost-free diagrams, with the same symmetry and overall factors, cuts and
/// numerators; only the names of the unfiltered run differ.
#[test]
fn in_flight_filters_select_like_vetoes_and_post_filtering() {
    let case = Case::named("qcd_gggg_L1");
    let model = case.model();
    let options = case.options(&model, 1);
    let all = case.generate_with(model.clone(), &options).unwrap();
    let ghost_free: BTreeSet<_> = all
        .diagrams
        .iter()
        .filter(|diagram| !is_ghost(&model, diagram))
        .map(|diagram| diagram.id().unwrap())
        .collect();
    assert!(ghost_free.len() < all.diagrams.len());
    let expected = Census::of(&all).restricted(|entry| ghost_free.contains(&entry.id));
    let filtered = census(Case::named("qcd_gggg_ghost_filter_L1"), 1);
    let vetoed = census(Case::named("qcd_gggg_ghost_veto_L1"), 1);
    expected.assert_same(&filtered, Scope::Content, "user filter vs post-filtering");
    filtered.assert_same(&vetoed, Scope::Content, "user filter vs particle veto");
}

/// Counterterm copies of vertex rules multiply the rule choices: colouring the `k` vertices of a
/// diagram with automorphism group `A` that use a duplicated rule in two ways gives orbits whose
/// stabilizers satisfy `sum 1/|Stab| = 2^k / |A|` (Burnside). The vacuum and tadpole cases have
/// automorphisms that move vertices with duplicated rules, so only the orbit-minimal rule choices
/// may survive.
#[test]
fn counterterm_copies_weight_every_rule_choice() {
    for (base, counterterms) in [
        ("phi3_2pt_L1", "ct_phi3_2pt_L1"),
        ("phi3_2pt_tad_L2", "ct_phi3_2pt_tad_L2"),
        ("phi3_vac_L2", "ct_phi3_vac_L2"),
        ("phi3_vac_L3", "ct_phi3_vac_L3"),
        ("phi3_vac_L4", "ct_phi3_vac_L4"),
        ("qcd_gggg_L1", "ct_qcd_gggg_L1"),
        ("qcd_uugg_L1", "ct_qcd_uugg_L1"),
        ("qcd_vac_L2", "ct_qcd_vac_L2"),
        ("qcd_vac_L3", "ct_qcd_vac_L3"),
    ] {
        let counterterms = Case::named(counterterms);
        let duplicated: BTreeSet<_> = counterterms
            .model()
            .vertex_rules()
            .iter()
            .filter_map(|rule| rule.name.strip_suffix("_CT").map(str::to_owned))
            .collect();
        let case = Case::named(base);
        let model = case.model();
        let options = case.options(&model, 1);
        let mut expected = Rational::zero();
        for diagram in case
            .generate_with(model.clone(), &options)
            .unwrap()
            .diagrams
        {
            let k = diagram
                .vertices()
                .filter(|(_, vertex)| {
                    let rule = model
                        .vertex_rule_by_id(vertex.interaction.unwrap())
                        .unwrap();
                    duplicated.contains(&rule.name)
                })
                .count();
            expected += Rational::from((1 << k, diagram.symmetry_factor()));
        }
        let actual = census(counterterms, 1).summary().inverse_symmetry;
        assert_eq!(actual, expected, "{}", counterterms.name);
    }
}

/// Counterterm copies with their own `CT` order: the rule choices with `k` copies among the `m`
/// vertices of a base diagram (automorphism group `A`) that can take one form orbits with `sum
/// 1/|Stab| = C(m, k) / |A|` (Burnside on the `k`-subsets), so a range of the order keeps exactly
/// the orbit-minimal choices with `k` in range, where the bound on partial assignments cannot tell
/// a rule from its copy.
#[test]
fn counterterm_orders_select_rule_choices_by_count() {
    for (base, counterterms, range) in [
        ("qcd_gggg_L1", "ctorder_qcd_gggg_L1", 1..=1),
        ("qcd_uugg_L1", "ctorder_qcd_uugg_L1", 1..=1),
        ("qcd_vac_L3", "ctorder_qcd_vac_L3", 0..=1),
        ("phi3_vac_L4", "ctorder_phi3_vac_L4", 2..=2),
    ] {
        let counterterms = Case::named(counterterms);
        let duplicated: BTreeSet<_> = counterterms
            .model()
            .vertex_rules()
            .iter()
            .filter_map(|rule| rule.name.strip_suffix("_CT").map(str::to_owned))
            .collect();
        let case = Case::named(base);
        let model = case.model();
        let options = case.options(&model, 1);
        let binomial = |m: u64, k: u64| (0..k).fold(1, |c, i| c * (m - i) / (i + 1));
        let mut expected = Rational::zero();
        for diagram in case
            .generate_with(model.clone(), &options)
            .unwrap()
            .diagrams
        {
            let m = diagram
                .vertices()
                .filter(|(_, vertex)| {
                    let rule = model
                        .vertex_rule_by_id(vertex.interaction.unwrap())
                        .unwrap();
                    duplicated.contains(&rule.name)
                })
                .count() as u64;
            for k in range.clone().filter(|&k| k <= m) {
                expected += Rational::from((binomial(m, k), diagram.symmetry_factor()));
            }
        }
        let actual = census(counterterms, 1).summary().inverse_symmetry;
        assert_eq!(actual, expected, "{}", counterterms.name);
    }
}
