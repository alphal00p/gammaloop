//! Differential oracles for diagram generation.
//!
//! Every successful case runs both class strategies with Named and Minimal output at
//! 1, 8 and 32 threads. The particle-only reference also certifies its historical row.
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
//! case's row instead of comparing it. `FEYNKIT_ORACLE_DUMP=<dir>` explicitly captures full
//! JSON inventories and the historical text references. Normally only mismatches write full
//! inventories; successful comparisons use exact typed payloads in anonymous files.

mod cases;
mod census;
mod classes;
mod frame_certificate;
mod frames;

use feynkit_graph::symbols;
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::process::Command;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use cases::{CASES, Case, REFERENCE};
use census::{Census, Scope, Summary};
use feynkit_generator::{
    BuiltDiagramSet, Classes, GenerationControl, GraphGroupingOptions, Grouping, Naming,
    NumeratorOptions, Run, Selection, Stage,
};
use feynkit_graph::{FeynmanDiagram, NumeratorInterface};
use feynkit_model::Model;
use idenso::{CookMode, CookSettings, CookSourceFilter};
use spenso::structure::representation::{Minkowski, RepName};
use spenso::{index_symbol, tensor_symbol};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, UserData};
use symbolica::domains::rational::Rational;
use symbolica::function;
use symbolica::state::State;

fn threads() -> Vec<usize> {
    env::var("FEYNKIT_ORACLE_THREADS").map_or(vec![1, 8, 32], |threads| {
        threads.split(',').map(|t| t.parse().unwrap()).collect()
    })
}

fn census(case: &Case, threads: usize) -> Census {
    census_with_naming(case, threads, Naming::prefix("D"), Classes::Auto)
}

fn census_with_naming(case: &Case, threads: usize, naming: Naming, classes: Classes) -> Census {
    let model = case.model();
    let request = case.request(&model).with_classes(classes);
    let run = Run::default().threads(threads);
    let result = request
        .generate(model.into(), &run)
        .and_then(|prepared| prepared.build(&case.numerators().with_naming(naming), &run))
        .and_then(|built| built.group(&case.grouping(), &run))
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
    assert!(result.report().completed, "{} did not complete", case.name);
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
    // The particle-only reference is also the historical oracle; generation and the
    // Named/Minimal builds are shared with the full class/thread comparison matrix.
    let reference = classes::check(case);
    let numerators = case.numerators().with_naming(Naming::prefix("D"));
    // The source name/index assertion is specific to these selector-free fixtures.
    let serialized = serde_json::to_value(&numerators).unwrap();
    assert_eq!(
        serialized["naming"],
        serde_json::json!({"named": {"prefix": "D"}})
    );
    let serialized = serde_json::to_value(Selection::new()).unwrap();
    for field in ["selected_diagram_ids", "selected_diagram_names"] {
        assert!(serialized.get(field).unwrap().is_null(), "{name}: {field}");
    }
    for field in ["vetoed_diagram_ids", "vetoed_diagram_names"] {
        assert_eq!(serialized[field], serde_json::json!([]), "{name}: {field}");
    }
    if let Ok(directory) = env::var("FEYNKIT_ORACLE_DUMP") {
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(format!("{directory}/{name}.census"), reference.dump()).unwrap();
    }
    let summary = reference.summary();
    if env::var_os("FEYNKIT_ORACLE_RECORD").is_some() {
        println!("    {tier} {name}: {summary};");
    } else {
        let expected = Summary::from(&row);
        if summary != expected {
            reference.capture_failure(&reference, &format!("{name}: historical row"));
        }
        assert_eq!(
            summary, expected,
            "{name}: summary differs from its oracle row"
        );
    }
}

/// Run `name` at every thread count and require it to fail with an error containing every
/// fragment of `error`.
// No case fails at present: native tensor grouping evaluates the QCD, two-loop QED and SM
// grouping cases whose samples failed before, and their rows are recorded in the heavy tier.
#[allow(dead_code)]
fn fails(name: &str, error: &[&str]) {
    let case = Case::named(name);
    for count in threads() {
        for output in [Naming::prefix("D"), Naming::None] {
            let model = case.model();
            let run = Run::default().threads(count);
            let result = case
                .request(&model)
                .generate(model.into(), &run)
                .and_then(|prepared| {
                    prepared.build(&case.numerators().with_naming(output.clone()), &run)
                })
                .and_then(|built| built.group(&case.grouping(), &run));
            match result {
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
    fast phi34_2pt_identical_L3: Row { diagrams: 72, topologies: 72, assignments: 72, inverse_symmetry: (71, 3), weights: (71, 3), symmetry_factors: &[(1, 2), (2, 24), (4, 32), (6, 1), (8, 10), (12, 3)], content: 0x6454c5db351b672c, numerators: 0xe64b2d72fd7c2e4e, named: 0x175ac2b33292f99f, ids: 0xee56a211744f234e };
    fast qed_eemm_L1: Row { diagrams: 17, topologies: 10, assignments: 17, inverse_symmetry: (17, 1), weights: (1, 1), symmetry_factors: &[(1, 17)], content: 0x67015758ee0dfcf5, numerators: 0x3866ebfda57f6fb8, named: 0x92ff7810a94d1ff7, ids: 0x1dc5ec5f97f654ff };
    fast qed_eemm_L2: Row { diagrams: 357, topologies: 116, assignments: 357, inverse_symmetry: (357, 1), weights: (57, 1), symmetry_factors: &[(1, 357)], content: 0x459953470ad8f543, numerators: 0xe5d00cfbb0473ed7, named: 0x04ccbdf50c929329, ids: 0x566e38b1fd985ea0 };
    heavy qed_eemm_L3: Row { diagrams: 8724, topologies: 1245, assignments: 8724, inverse_symmetry: (8724, 1), weights: (-282, 1), symmetry_factors: &[(1, 8724)], content: 0x03fc3b8d254ed8be, numerators: 0xa35949f16e378ad3, named: 0x661d07aebc62dae1, ids: 0xb51f888639fdfb91 };
    fast qed_ee4f_L0: Row { diagrams: 6, topologies: 7, assignments: 6, inverse_symmetry: (6, 1), weights: (-6, 1), symmetry_factors: &[(1, 6)], content: 0x66f2c9719dc08f84, numerators: 0x7f2a34ae141a6f7c, named: 0x8a27f32a3acc83fd, ids: 0xd5ff3aec66c9a52e };
    fast qed_ee4f_L1: Row { diagrams: 242, topologies: 304, assignments: 242, inverse_symmetry: (242, 1), weights: (10, 1), symmetry_factors: &[(1, 242)], content: 0x8de28c7c413ceda2, numerators: 0xaaf41914ee0bedd0, named: 0x76071173c740ae8c, ids: 0x8009b4ab595db4d2 };
    heavy qed_ee4f_L2: Row { diagrams: 8460, topologies: 7909, assignments: 8460, inverse_symmetry: (8460, 1), weights: (828, 1), symmetry_factors: &[(1, 8460)], content: 0x2b83b18b522b9b8a, numerators: 0x8c8ad7229b4ee767, named: 0x9ffdbc9c8b5864e4, ids: 0xea62a12dfd838cd6 };
    fast qed_eeaaa_L1: Row { diagrams: 384, topologies: 120, assignments: 384, inverse_symmetry: (384, 1), weights: (264, 1), symmetry_factors: &[(1, 384)], content: 0xd046854816e13604, numerators: 0xdd48fecc16d94ab7, named: 0xd2334763314c3f01, ids: 0x045b832925723c15 };
    heavy qed_eeaaa_L2: Row { diagrams: 16344, topologies: 1581, assignments: 16344, inverse_symmetry: (16344, 1), weights: (-4356, 1), symmetry_factors: &[(1, 16344)], content: 0x51d8376445158a3b, numerators: 0x537a65a9c797d70b, named: 0x5d5f16e150dd124d, ids: 0x6939a25e707794a1 };
    fast qed_eeaa_identical_L1: Row { diagrams: 66, topologies: 20, assignments: 66, inverse_symmetry: (66, 1), weights: (42, 1), symmetry_factors: &[(1, 66)], content: 0x8adeb1c02c1e722c, numerators: 0xf8be6f67ed0ba5dc, named: 0x628eec07876a65be, ids: 0x79a2cf2343cb733c };
    fast qed_eeaa_sign_L1: Row { diagrams: 48, topologies: 20, assignments: 66, inverse_symmetry: (48, 1), weights: (24, 1), symmetry_factors: &[(1, 48)], content: 0x0a00287a45ade68b, numerators: 0x284efd9f5ab9e743, named: 0x0f98ce186f74a624, ids: 0xb06d936af435481d };
    fast qed_eeaa_scalar_L1: Row { diagrams: 48, topologies: 20, assignments: 66, inverse_symmetry: (48, 1), weights: (24, 1), symmetry_factors: &[(1, 48)], content: 0x0a00287a45ade68b, numerators: 0x284efd9f5ab9e743, named: 0x0f98ce186f74a624, ids: 0xb06d936af435481d };
    fast qed_ee4f_identical_L1: Row { diagrams: 242, topologies: 304, assignments: 242, inverse_symmetry: (242, 1), weights: (10, 1), symmetry_factors: &[(1, 242)], content: 0x8de28c7c413ceda2, numerators: 0xaaf41914ee0bedd0, named: 0x092d3ac75cb11a5c, ids: 0x8009b4ab595db4d2 };
    fast qed_ee4f_sign_L1: Row { diagrams: 224, topologies: 304, assignments: 242, inverse_symmetry: (224, 1), weights: (-8, 1), symmetry_factors: &[(1, 224)], content: 0x9c93d2d3585b5388, numerators: 0xbd83a5fba62d76c8, named: 0x28f34eb5ea793feb, ids: 0x4d7b2538b8c686f5 };
    fast qed_ee4f_scalar_L1: Row { diagrams: 224, topologies: 304, assignments: 242, inverse_symmetry: (224, 1), weights: (-8, 1), symmetry_factors: &[(1, 224)], content: 0x9c93d2d3585b5388, numerators: 0xbd83a5fba62d76c8, named: 0x28f34eb5ea793feb, ids: 0x4d7b2538b8c686f5 };
    heavy qed_eemm_identical_L2: Row { diagrams: 348, topologies: 116, assignments: 357, inverse_symmetry: (348, 1), weights: (57, 1), symmetry_factors: &[(1, 348)], content: 0x4485709614bfe0d9, numerators: 0x2a0afe5fe6705dc7, named: 0xc7214f44d613ae46, ids: 0xa5953826d6302b2f };
    heavy qed_eemm_sign_L2: Row { diagrams: 312, topologies: 116, assignments: 357, inverse_symmetry: (312, 1), weights: (21, 1), symmetry_factors: &[(1, 312)], content: 0xc5360cb6f9d7e413, numerators: 0x6d7ddc98c456120b, named: 0x68cfe60bc105d7f5, ids: 0xe292899c57289190 };
    heavy qed_eemm_scalar_L2: Row { diagrams: 312, topologies: 116, assignments: 357, inverse_symmetry: (312, 1), weights: (21, 1), symmetry_factors: &[(1, 312)], content: 0xc5360cb6f9d7e413, numerators: 0x6d7ddc98c456120b, named: 0x68cfe60bc105d7f5, ids: 0xe292899c57289190 };
    heavy qed_eeaa_identical_L2: Row { diagrams: 1755, topologies: 176, assignments: 1818, inverse_symmetry: (1755, 1), weights: (-306, 1), symmetry_factors: &[(1, 1755)], content: 0x8351b3ea5e968195, numerators: 0x5bb4d863412c47df, named: 0x2bb64543e9487d90, ids: 0xdc46244e960d9a1d };
    heavy qed_eeaa_sign_L2: Row { diagrams: 999, topologies: 176, assignments: 1818, inverse_symmetry: (999, 1), weights: (-90, 1), symmetry_factors: &[(1, 999)], content: 0x3c592cbf406cde22, numerators: 0x9abd4fff2057303c, named: 0x23e0b0f899dc0308, ids: 0x7a5da2ab04c8fc81 };
    heavy qed_eeaa_scalar_L2: Row { diagrams: 999, topologies: 176, assignments: 1818, inverse_symmetry: (999, 1), weights: (-90, 1), symmetry_factors: &[(1, 999)], content: 0x3c592cbf406cde22, numerators: 0x9abd4fff2057303c, named: 0x23e0b0f899dc0308, ids: 0x7a5da2ab04c8fc81 };
    fast qed_ase_L3: Row { diagrams: 540, topologies: 10, assignments: 540, inverse_symmetry: (540, 1), weights: (270, 1), symmetry_factors: &[(1, 540)], content: 0x5f8145dd3d129246, numerators: 0x65d2b078e34de88f, named: 0xba9713503bfc0b8e, ids: 0x3bb77701c9a055a5 };
    heavy qed_ase_L4: Row { diagrams: 12285, topologies: 64, assignments: 12285, inverse_symmetry: (12285, 1), weights: (189, 1), symmetry_factors: &[(1, 12285)], content: 0x5b713868f60ec623, numerators: 0x7239a1017a75242a, named: 0xc54b216c87ef9aa2, ids: 0xa324ec8a029594ff };
    heavy qed_ase_L5: Row { diagrams: 322704, topologies: 519, assignments: 322704, inverse_symmetry: (322704, 1), weights: (-35478, 1), symmetry_factors: &[(1, 322704)], content: 0x511cabd7f8d70b65, numerators: 0x91688536d978a431, named: 0xffefbae154730d8e, ids: 0x4c86ab9e11bbc114 };
    fast qed_ese_L3: Row { diagrams: 172, topologies: 10, assignments: 172, inverse_symmetry: (172, 1), weights: (-10, 1), symmetry_factors: &[(1, 172)], content: 0x3b76d98bd234f769, numerators: 0xa7f9a5a48bbcede4, named: 0x798f5dfb32643f81, ids: 0x2391d3c49b9fa520 };
    heavy qed_ese_L4: Row { diagrams: 3485, topologies: 64, assignments: 3485, inverse_symmetry: (3485, 1), weights: (-389, 1), symmetry_factors: &[(1, 3485)], content: 0x32d9e5cb0b5d366e, numerators: 0x1327d8e11732f3df, named: 0x597790dff753eb2a, ids: 0x4b5eec491dcd879a };
    fast qcd_gggg_L1: Row { diagrams: 314, topologies: 55, assignments: 314, inverse_symmetry: (294, 1), weights: (-224, 1), symmetry_factors: &[(1, 274), (2, 40)], content: 0x5589efce8dac4b10, numerators: 0x8b23170b9f0930ac, named: 0x09c1775d49013d56, ids: 0x19ccb95e882bbda8 };
    heavy qcd_gggg_L2: Row { diagrams: 15253, topologies: 840, assignments: 15253, inverse_symmetry: (163229, 12), weights: (24377, 12), symmetry_factors: &[(1, 12086), (2, 2905), (4, 243), (6, 19)], content: 0xc3df8bb42155219d, numerators: 0xed78bbe59ad71b21, named: 0x4c7efc397c54a13f, ids: 0x2eac0a1111b9d2f9 };
    heavy qcd_gggg_L3: Row { diagrams: 686501, topologies: 14340, assignments: 686501, inverse_symmetry: (4715477, 8), weights: (1169243, 24), symmetry_factors: &[(1, 505145), (2, 157099), (4, 21237), (6, 1512), (8, 1277), (12, 231)], content: 0xbfcc378cb7bf2f09, numerators: 0x33c04d98409fa462, named: 0xe872f37eae017892, ids: 0xb37b49bc02c2fa26 };
    heavy qcd_ggggg_L1: Row { diagrams: 3075, topologies: 562, assignments: 3075, inverse_symmetry: (2895, 1), weights: (-2131, 1), symmetry_factors: &[(1, 2715), (2, 360)], content: 0x5d9d3829fa9f6a7d, numerators: 0x6af1a7a62d153fdb, named: 0x0e4596d97b1dbda4, ids: 0xa50aa7fed0d1ef87 };
    heavy qcd_ggggg_L2: Row { diagrams: 213390, topologies: 11685, assignments: 213390, inverse_symmetry: (765585, 4), weights: (131385, 4), symmetry_factors: &[(1, 170990), (2, 39280), (4, 2955), (6, 165)], content: 0x6767ccf69132de38, numerators: 0xb9a673a3d2a4af42, named: 0xa30204b0fc4fcbca, ids: 0x5c319a4b19971744 };
    fast qcd_uugg_L1: Row { diagrams: 99, topologies: 32, assignments: 99, inverse_symmetry: (94, 1), weights: (32, 1), symmetry_factors: &[(1, 89), (2, 10)], content: 0xb4b9893f4de07ca5, numerators: 0xf47bded5c185d2ee, named: 0x349113f43a3eb3e2, ids: 0x683a0083a05158be };
    heavy qcd_uugg_L2: Row { diagrams: 3206, topologies: 457, assignments: 3206, inverse_symmetry: (17479, 6), weights: (1925, 6), symmetry_factors: &[(1, 2644), (2, 517), (4, 38), (6, 7)], content: 0x9832af97da94dbd2, numerators: 0x928b57e2de93e9b2, named: 0x4eee2068848bc1db, ids: 0x1ec12bf1529e303f };
    heavy qcd_uugg_L3: Row { diagrams: 114301, topologies: 7513, assignments: 114301, inverse_symmetry: (301202, 3), weights: (-17197, 6), symmetry_factors: &[(1, 88202), (2, 22904), (4, 2686), (6, 317), (8, 152), (12, 40)], content: 0x9432fa0a003fc794, numerators: 0xa37c261b5ab81a2b, named: 0xfd5c7d2836314cb7, ids: 0xe4ce63ab3f2d728a };
    fast qcd_uuuu_L2: Row { diagrams: 902, topologies: 230, assignments: 902, inverse_symmetry: (5057, 6), weights: (0, 1), symmetry_factors: &[(1, 788), (2, 106), (4, 6), (6, 2)], content: 0xab4face31b470e62, numerators: 0x79443744b56d9ca7, named: 0x896adc3d0b635bae, ids: 0xf8751099f0ea8916 };
    fast qcd_gse_L3: Row { diagrams: 1219, topologies: 57, assignments: 1219, inverse_symmetry: (8553, 8), weights: (405, 8), symmetry_factors: &[(1, 940), (2, 241), (4, 31), (6, 1), (8, 5), (12, 1)], content: 0xdfd63c0cec5bc2f6, numerators: 0x0ab23cb9c2fb1f80, named: 0x97bd6cb32cb134cc, ids: 0x558a5800bdfc6874 };
    heavy qcd_gse_L4: Row { diagrams: 35820, topologies: 666, assignments: 35820, inverse_symmetry: (1485389, 48), weights: (81159, 16), symmetry_factors: &[(1, 26801), (2, 7660), (4, 1140), (6, 83), (8, 112), (12, 14), (16, 9), (24, 1)], content: 0x4bffa26c554c174c, numerators: 0x0d2ced572977ae4d, named: 0xbd25d6bf601f99a0, ids: 0xf35f798de1eb9af6 };
    heavy qcd_gse_L5: Row { diagrams: 1267351, topologies: 9892, assignments: 1267351, inverse_symmetry: (308996023, 288), weights: (10378459, 288), symmetry_factors: &[(1, 910237), (2, 298126), (3, 14), (4, 49275), (6, 3114), (8, 5369), (12, 752), (16, 367), (24, 72), (32, 19), (36, 1), (48, 3), (72, 2)], content: 0x4b5acdd1d4e850b7, numerators: 0xfadfbeb3dabe2d3b, named: 0xfb62c8394c2f1907, ids: 0x895057e0154d8869 };
    heavy qcd_gse_all_L3: Row { diagrams: 2627, topologies: 72, assignments: 2627, inverse_symmetry: (4427, 2), weights: (1997, 6), symmetry_factors: &[(1, 1871), (2, 626), (4, 102), (6, 15), (8, 10), (12, 3)], content: 0x6b52ae23066e5db5, numerators: 0x8f9ebf306d3c237d, named: 0x79cae111974fc394, ids: 0x39d18f9226913964 };
    fast qcd_tad_L2: Row { diagrams: 125, topologies: 6, assignments: 125, inverse_symmetry: (2383, 24), weights: (1207, 24), symmetry_factors: &[(1, 77), (2, 42), (4, 4), (6, 1), (8, 1)], content: 0x81b089b52d6a303d, numerators: 0x21fae42c6127481c, named: 0x791919f7ad5213b7, ids: 0x5d1a73141b37b5ef };
    heavy qcd_tad_L3: Row { diagrams: 3490, topologies: 46, assignments: 3490, inverse_symmetry: (124199, 48), weights: (-30613, 48), symmetry_factors: &[(1, 1883), (2, 1247), (3, 7), (4, 271), (6, 29), (8, 43), (12, 5), (16, 4), (48, 1)], content: 0x928a721d193787c7, numerators: 0x969a39ce5c58b16f, named: 0xe7f788851336165c, ids: 0xa7f746480b1ecc18 };
    fast qcd_vac_L2: Row { diagrams: 45, topologies: 3, assignments: 45, inverse_symmetry: (191, 6), weights: (107, 6), symmetry_factors: &[(1, 21), (2, 21), (8, 2), (12, 1)], content: 0x4895248b9117da28, numerators: 0x02cc567f8ebfee74, named: 0x895ea4bbc0107f71, ids: 0xf3bf59ff814c96a1 };
    fast qcd_vac_L3: Row { diagrams: 649, topologies: 12, assignments: 649, inverse_symmetry: (2417, 6), weights: (-107, 1), symmetry_factors: &[(1, 231), (2, 287), (3, 7), (4, 84), (6, 14), (8, 18), (12, 1), (16, 4), (24, 1), (48, 2)], content: 0xe96ce0588b5af67f, numerators: 0x8a03ee9b8f7deaa6, named: 0xb72efadb8b75cbff, ids: 0x4a46bf7bd28feb58 };
    heavy qcd_vac_L4: Row { diagrams: 14017, topologies: 73, assignments: 14017, inverse_symmetry: (2597681, 288), weights: (284153, 288), symmetry_factors: &[(1, 5698), (2, 5383), (3, 154), (4, 1894), (6, 189), (8, 483), (12, 79), (16, 85), (24, 21), (32, 14), (48, 13), (72, 2), (128, 1), (384, 1)], content: 0x39984039ef0665b1, numerators: 0x1d29811dc3205779, named: 0xae16b31a96c8d24a, ids: 0xfd7cf5edf27cf381 };
    fast qcd_vac1pi_L4: Row { diagrams: 628, topologies: 26, assignments: 628, inverse_symmetry: (11465, 36), weights: (446, 9), symmetry_factors: &[(1, 154), (2, 217), (3, 56), (4, 101), (6, 28), (8, 48), (12, 8), (16, 8), (24, 2), (32, 2), (48, 3), (72, 1)], content: 0xdcce6c75d4e4ddd6, numerators: 0xb06d25f1c8249a11, named: 0xf15f24a887f7696b, ids: 0xcb6b2d23d49573d1 };
    heavy qcd_vac1pi_L5: Row { diagrams: 9889, topologies: 166, assignments: 9889, inverse_symmetry: (7313267, 1152), weights: (227027, 1152), symmetry_factors: &[(1, 4165), (2, 3455), (4, 1494), (6, 49), (8, 472), (12, 46), (16, 159), (24, 11), (32, 22), (48, 6), (64, 4), (96, 2), (128, 3), (144, 1)], content: 0x25da5b1f7bcbb040, numerators: 0xf80bfb8fd9a54875, named: 0x4dfc72054c772df5, ids: 0x68f8fcac075cf35b };
    fast qcd_gggg_sym_L1: Row { diagrams: 123, topologies: 25, assignments: 123, inverse_symmetry: (147, 2), weights: (-56, 1), symmetry_factors: &[(1, 37), (2, 62), (4, 20), (8, 4)], content: 0x7ec66a5e0aa87d96, numerators: 0x6ce235818c866da5, named: 0x1dc6e0df0acac807, ids: 0xa797483cc898596b };
    heavy qcd_gggg_identical_L1: Row { diagrams: 314, topologies: 55, assignments: 314, inverse_symmetry: (294, 1), weights: (-224, 1), symmetry_factors: &[(1, 274), (2, 40)], content: 0x5589efce8dac4b10, numerators: 0x8b23170b9f0930ac, named: 0xa83f08da27fc45a4, ids: 0x19ccb95e882bbda8 };
    heavy qcd_gggg_sign_L1: Row { diagrams: 314, topologies: 55, assignments: 314, inverse_symmetry: (294, 1), weights: (-224, 1), symmetry_factors: &[(1, 274), (2, 40)], content: 0x5589efce8dac4b10, numerators: 0x8b23170b9f0930ac, named: 0xa83f08da27fc45a4, ids: 0x19ccb95e882bbda8 };
    heavy qcd_gggg_scalar_L1: Row { diagrams: 314, topologies: 55, assignments: 314, inverse_symmetry: (294, 1), weights: (-224, 1), symmetry_factors: &[(1, 274), (2, 40)], content: 0x5589efce8dac4b10, numerators: 0x8b23170b9f0930ac, named: 0xa83f08da27fc45a4, ids: 0x19ccb95e882bbda8 };
    fast qcd_gggg_ghost_filter_L1: Row { diagrams: 277, topologies: 55, assignments: 277, inverse_symmetry: (257, 1), weights: (-187, 1), symmetry_factors: &[(1, 237), (2, 40)], content: 0x90a1ce06e025a6f5, numerators: 0x441b378f711ce4fc, named: 0x1e9b4de5869edb85, ids: 0xda3302c00009cd6d };
    fast qcd_gggg_ghost_veto_L1: Row { diagrams: 277, topologies: 55, assignments: 277, inverse_symmetry: (257, 1), weights: (-187, 1), symmetry_factors: &[(1, 237), (2, 40)], content: 0x90a1ce06e025a6f5, numerators: 0x441b378f711ce4fc, named: 0x1e9b4de5869edb85, ids: 0xda3302c00009cd6d };
    fast ym_g6_L0: Row { diagrams: 220, topologies: 220, assignments: 220, inverse_symmetry: (220, 1), weights: (220, 1), symmetry_factors: &[(1, 220)], content: 0xe166545bc76a8cd8, numerators: 0xc742a6c3a2526c64, named: 0x0855fe3d5967a4de, ids: 0x4a0b7e227d7ace53 };
    heavy ym_g7_L0: Row { diagrams: 2485, topologies: 2485, assignments: 2485, inverse_symmetry: (2485, 1), weights: (2485, 1), symmetry_factors: &[(1, 2485)], content: 0x057203fcdeaf182e, numerators: 0xdb59060eb8d155e8, named: 0x825ef5b96a66dc8d, ids: 0x54cd7d6e49c55438 };
    heavy ym_g8_L0: Row { diagrams: 34300, topologies: 34300, assignments: 34300, inverse_symmetry: (34300, 1), weights: (34300, 1), symmetry_factors: &[(1, 34300)], content: 0x0d8cd8312f854d8b, numerators: 0x4be08c2dd2eb4a1e, named: 0xe921e3ce642cddae, ids: 0xa9980ccebbd89562 };
    heavy ym_gse_L4: Row { diagrams: 2571, topologies: 666, assignments: 2571, inverse_symmetry: (81653, 48), weights: (-8113, 16), symmetry_factors: &[(1, 1139), (2, 889), (4, 402), (6, 17), (8, 100), (12, 14), (16, 9), (24, 1)], content: 0xf0b98dad15d745b1, numerators: 0x6bb69664d65cdbbc, named: 0x24757abcc9b9d2ec, ids: 0x0074f4d0f4b0cda4 };
    heavy ym_gse_L5: Row { diagrams: 53201, topologies: 9892, assignments: 53201, inverse_symmetry: (9957739, 288), weights: (-2104793, 288), symmetry_factors: &[(1, 22781), (2, 18577), (3, 2), (4, 8421), (6, 315), (8, 2411), (12, 248), (16, 349), (24, 72), (32, 19), (36, 1), (48, 3), (72, 2)], content: 0x8da6be7bb6b31c5b, numerators: 0xf988e70b630fac54, named: 0x81b24addd62ad96f, ids: 0x683c20a191922cf0 };
    fast ym_ggg_zero_color_L3: Row { diagrams: 1217, topologies: 469, assignments: 1340, inverse_symmetry: (3583, 4), weights: (-1161, 4), symmetry_factors: &[(1, 649), (2, 424), (4, 132), (6, 6), (8, 6)], content: 0x19b54d6ef892c24e, numerators: 0xae5ed16f90c02201, named: 0x3569d286d4aa397d, ids: 0x1ed193f8aeb76d2f };
    fast ym_vac_L5: Row { diagrams: 381, topologies: 91, assignments: 381, inverse_symmetry: (95437, 576), weights: (-19595, 576), symmetry_factors: &[(1, 82), (2, 110), (4, 73), (6, 4), (8, 54), (12, 9), (16, 28), (24, 3), (32, 8), (48, 4), (64, 2), (96, 1), (128, 2), (144, 1)], content: 0xf7c775391f48c05c, numerators: 0xc2afd3e35908cbf1, named: 0x357bb19f8902aa0c, ids: 0x41ea11f42e0c80b0 };
    fast sm_eemm_L0: Row { diagrams: 4, topologies: 1, assignments: 4, inverse_symmetry: (4, 1), weights: (-4, 1), symmetry_factors: &[(1, 4)], content: 0x97af00ae91ef62c4, numerators: 0xc030c29958e57000, named: 0xc01c07a8158bc53a, ids: 0x91a782bcdd47cb2a };
    fast sm_eemm_L1: Row { diagrams: 464, topologies: 10, assignments: 464, inverse_symmetry: (925, 2), weights: (-205, 2), symmetry_factors: &[(1, 461), (2, 3)], content: 0x7d459f5367fbcba8, numerators: 0xdf611964612f3ec9, named: 0x7b2f048149dccbf3, ids: 0xa1e6a4f3cf7051f9 };
    heavy sm_eemm_L2: Row { diagrams: 81586, topologies: 160, assignments: 81586, inverse_symmetry: (486455, 6), weights: (26947, 6), symmetry_factors: &[(1, 80580), (2, 978), (4, 26), (6, 2)], content: 0xb858238d81f548bf, numerators: 0xd6b7c232c4c96640, named: 0xf7669cdbae2ed504, ids: 0xa7763c7a1fb23595 };
    heavy sm_eemm_identical_L1: Row { diagrams: 412, topologies: 10, assignments: 464, inverse_symmetry: (821, 2), weights: (-277, 2), symmetry_factors: &[(1, 409), (2, 3)], content: 0xb1581340d4b86f24, numerators: 0x21515b5976c706a1, named: 0xb60f369359601bab, ids: 0x406fcac875016fbb };
    heavy sm_eemm_sign_L1: Row { diagrams: 400, topologies: 10, assignments: 464, inverse_symmetry: (797, 2), weights: (-301, 2), symmetry_factors: &[(1, 397), (2, 3)], content: 0x3de986f3be431158, numerators: 0x9c6877c528487c8c, named: 0xb6f8b6c32f1dce1e, ids: 0xbc1c8afb96166bb3 };
    heavy sm_eemm_scalar_L1: Row { diagrams: 400, topologies: 10, assignments: 464, inverse_symmetry: (797, 2), weights: (-301, 2), symmetry_factors: &[(1, 397), (2, 3)], content: 0x3de986f3be431158, numerators: 0x9c6877c528487c8c, named: 0xb6f8b6c32f1dce1e, ids: 0xbc1c8afb96166bb3 };
    fast sm_ee4f_L0: Row { diagrams: 44, topologies: 7, assignments: 44, inverse_symmetry: (44, 1), weights: (-44, 1), symmetry_factors: &[(1, 44)], content: 0x3e3c8eedad8d3f82, numerators: 0xd54136b7b108d863, named: 0x30c5df45e3e970de, ids: 0x9190a3b6de2a021a };
    heavy sm_ee4f_L1: Row { diagrams: 15296, topologies: 322, assignments: 15296, inverse_symmetry: (15278, 1), weights: (-6946, 1), symmetry_factors: &[(1, 15260), (2, 36)], content: 0xe2b566a2ff0632d5, numerators: 0x04501e488fa99a50, named: 0x3a093a2587e658bd, ids: 0x26c8449dd7488b08 };
    fast sm_eeww_L0: Row { diagrams: 4, topologies: 2, assignments: 4, inverse_symmetry: (4, 1), weights: (-4, 1), symmetry_factors: &[(1, 4)], content: 0x7477c6f9f9a39d73, numerators: 0xec3019645ae8d3e5, named: 0x48249b27b144364c, ids: 0xf65454f1b8b1a032 };
    fast sm_eeww_L1: Row { diagrams: 981, topologies: 23, assignments: 981, inverse_symmetry: (978, 1), weights: (60, 1), symmetry_factors: &[(1, 975), (2, 6)], content: 0xb160c82b75af91f7, numerators: 0x726e5dc32f7fec13, named: 0xf83411efa7440cc5, ids: 0xdf39ae7ba47e78e3 };
    fast sm_ggtt_L1: Row { diagrams: 417, topologies: 32, assignments: 417, inverse_symmetry: (412, 1), weights: (-50, 1), symmetry_factors: &[(1, 407), (2, 10)], content: 0xb5a7c5fdc148c05d, numerators: 0xb5767f61570b570d, named: 0x372451345acf1071, ids: 0x0145b791db58d872 };
    fast sm_ggttg_L0: Row { diagrams: 16, topologies: 16, assignments: 16, inverse_symmetry: (16, 1), weights: (16, 1), symmetry_factors: &[(1, 16)], content: 0xeb9e3621b43c55e4, numerators: 0xc6956f2035f1c755, named: 0x66f155052cc8e9b4, ids: 0xfaacbb5903d5fcab };
    heavy sm_ggttg_L1: Row { diagrams: 3444, topologies: 322, assignments: 3444, inverse_symmetry: (3394, 1), weights: (-664, 1), symmetry_factors: &[(1, 3344), (2, 100)], content: 0x0f3d97ecdbc07803, numerators: 0x22ed8b171c1c1fd2, named: 0x13493b34f3aee8eb, ids: 0xaac450e707cb326d };
    fast sm_hse_L1: Row { diagrams: 20, topologies: 1, assignments: 20, inverse_symmetry: (37, 2), weights: (-11, 2), symmetry_factors: &[(1, 17), (2, 3)], content: 0xf75c143048b477e4, numerators: 0xe9fa747ef0e02184, named: 0x99bf9848fa996650, ids: 0x6de1b7fe0cca2268 };
    fast sm_hse_L2: Row { diagrams: 1452, topologies: 7, assignments: 1452, inverse_symmetry: (16823, 12), weights: (-7549, 12), symmetry_factors: &[(1, 1358), (2, 82), (4, 11), (6, 1)], content: 0x620a2ba4679c1969, numerators: 0x08e583c03f76c44f, named: 0x83940b678bd07c45, ids: 0xd23348e1347ccd2f };
    fast sm_ase_L2: Row { diagrams: 984, topologies: 7, assignments: 984, inverse_symmetry: (984, 1), weights: (-504, 1), symmetry_factors: &[(1, 984)], content: 0x25c64aa7157d8884, numerators: 0xcdbda25b4d328b67, named: 0x02d7af43f402ec7f, ids: 0x3d7366dbbeb2ce7c };
    fast sm_eemm_qcd0_L1: Row { diagrams: 464, topologies: 10, assignments: 464, inverse_symmetry: (925, 2), weights: (-205, 2), symmetry_factors: &[(1, 461), (2, 3)], content: 0x7d459f5367fbcba8, numerators: 0xdf611964612f3ec9, named: 0x7b2f048149dccbf3, ids: 0xa1e6a4f3cf7051f9 };
    heavy sm_eemm_qcd0_L2: Row { diagrams: 80722, topologies: 160, assignments: 80722, inverse_symmetry: (481271, 6), weights: (28675, 6), symmetry_factors: &[(1, 79716), (2, 978), (4, 26), (6, 2)], content: 0x74afb5046a18d92d, numerators: 0x7be8b1ae6036f6c0, named: 0x03cd896a8f7ad2aa, ids: 0xdf426aac3f6d7f53 };
    fast sm_uudd_qed2_L1: Row { diagrams: 186, topologies: 17, assignments: 186, inverse_symmetry: (371, 2), weights: (-71, 2), symmetry_factors: &[(1, 185), (2, 1)], content: 0x05c03352fe535f6a, numerators: 0x6997db3be8c077a4, named: 0x6f2ddb6c6a4e7b88, ids: 0x2087651fc4a2b3d4 };
    heavy sm_uudd_qed2_L2: Row { diagrams: 5461, topologies: 230, assignments: 5461, inverse_symmetry: (64169, 12), weights: (2803, 12), symmetry_factors: &[(1, 5236), (2, 221), (4, 3), (6, 1)], content: 0x9bfe72045f37142f, numerators: 0x8dd07baa70017f5c, named: 0xdada4f9094338015, ids: 0x8d40bb63c41eaa13 };
    fast ew_eemm_L1: Row { diagrams: 464, topologies: 10, assignments: 464, inverse_symmetry: (925, 2), weights: (-205, 2), symmetry_factors: &[(1, 461), (2, 3)], content: 0xa7f7ae7d2bdc72ec, numerators: 0xea66ba193c28de55, named: 0x1cbeeab0fcf33cd5, ids: 0x98f29369da185d5c };
    heavy ew_eemm_L2: Row { diagrams: 80722, topologies: 160, assignments: 80722, inverse_symmetry: (481271, 6), weights: (28675, 6), symmetry_factors: &[(1, 79716), (2, 978), (4, 26), (6, 2)], content: 0x77db78f8fc1fa24b, numerators: 0x0e37193fa8d3e4f3, named: 0x156bf7b516e462e3, ids: 0x58915abecd683c01 };
    fast ct_phi3_2pt_L1: Row { diagrams: 4, topologies: 1, assignments: 4, inverse_symmetry: (2, 1), weights: (2, 1), symmetry_factors: &[(2, 4)], content: 0x792c2689f9b727df, numerators: 0x094dd49c998c2915, named: 0x3217eb0128aeea77, ids: 0x247d40d238f195f1 };
    fast ct_qcd_gggg_L1: Row { diagrams: 1607, topologies: 55, assignments: 1607, inverse_symmetry: (2883, 2), weights: (-1549, 2), symmetry_factors: &[(1, 1276), (2, 331)], content: 0x762b4965ba414c96, numerators: 0x53082d1b0bf7fa97, named: 0xc648a3d718aff51c, ids: 0xcc8c4212ee546ab6 };
    fast ct_qcd_uugg_L1: Row { diagrams: 864, topologies: 32, assignments: 864, inverse_symmetry: (802, 1), weights: (-130, 1), symmetry_factors: &[(1, 740), (2, 124)], content: 0xc737f2677a9de6ce, numerators: 0x7e873daa2aa213e8, named: 0xdc7a0a093989054d, ids: 0xb3912abc03cca496 };
    fast ct_phi3_2pt_tad_L2: Row { diagrams: 148, topologies: 10, assignments: 148, inverse_symmetry: (50, 1), weights: (50, 1), symmetry_factors: &[(1, 4), (2, 44), (4, 92), (8, 8)], content: 0xf5122b0d42fad77e, numerators: 0xb1c2ba48ec6c4558, named: 0xa97896825197b385, ids: 0x6cdf274d98c848ce };
    fast ct_phi3_vac_L2: Row { diagrams: 6, topologies: 2, assignments: 6, inverse_symmetry: (5, 6), weights: (5, 6), symmetry_factors: &[(4, 1), (6, 1), (8, 2), (12, 2)], content: 0x9dad3c92417664f1, numerators: 0xf57c18fe6786fb79, named: 0x12031bbe318c7d29, ids: 0x9e59d4664215d958 };
    fast ct_phi3_vac_L3: Row { diagrams: 42, topologies: 5, assignments: 42, inverse_symmetry: (5, 1), weights: (5, 1), symmetry_factors: &[(4, 7), (6, 2), (8, 17), (16, 10), (24, 2), (48, 4)], content: 0x59da4aa68524faa8, numerators: 0xcadc67ebf2b44f13, named: 0xc86f79cbcf475f9e, ids: 0xceb477ab49af8581 };
    fast ct_phi3_vac_L4: Row { diagrams: 505, topologies: 17, assignments: 505, inverse_symmetry: (1105, 18), weights: (1105, 18), symmetry_factors: &[(1, 1), (2, 17), (4, 76), (6, 1), (8, 161), (12, 6), (16, 153), (24, 1), (32, 70), (36, 1), (48, 6), (64, 6), (72, 2), (128, 4)], content: 0xb78846c8c280e0f1, numerators: 0xc0711b1f8a34712b, named: 0xa9e96e243f9fd254, ids: 0xcd4c58688c00b527 };
    fast ct_qcd_vac_L2: Row { diagrams: 68, topologies: 3, assignments: 68, inverse_symmetry: (1103, 24), weights: (479, 24), symmetry_factors: &[(1, 29), (2, 32), (4, 1), (6, 1), (8, 3), (12, 2)], content: 0xb0401f2bf7eeb5c4, numerators: 0xf4dd16ce773dc439, named: 0x382ae441e068e697, ids: 0xbc728af284e8f8c0 };
    fast ct_qcd_vac_L3: Row { diagrams: 1816, topologies: 12, assignments: 1816, inverse_symmetry: (6341, 6), weights: (-321, 2), symmetry_factors: &[(1, 549), (2, 826), (3, 16), (4, 293), (6, 34), (8, 74), (12, 4), (16, 13), (24, 2), (48, 5)], content: 0x052ccef435c42ed0, numerators: 0xc901aee18233d77b, named: 0x007fe8f21acb2035, ids: 0xcbc0eae08c7c2d9e };
    fast xs_qed_eemm_L1: Row { diagrams: 1, topologies: 17, assignments: 34, inverse_symmetry: (1, 1), weights: (-1, 1), symmetry_factors: &[(1, 1)], content: 0x831ae534b022fd0f, numerators: 0x414e2a7a954beae8, named: 0xf40bfdc791015105, ids: 0xedd739fe2aef9d55 };
    fast xs_qed_eemm_L2: Row { diagrams: 30, topologies: 157, assignments: 714, inverse_symmetry: (30, 1), weights: (4, 1), symmetry_factors: &[(1, 30)], content: 0x6fcfea2ac818e3eb, numerators: 0xf889b37c1fd93143, named: 0x0b41c948729c87ec, ids: 0xc6537de685c0ec01 };
    fast xs_qed_eemm_L3: Row { diagrams: 824, topologies: 1536, assignments: 17448, inverse_symmetry: (824, 1), weights: (94, 1), symmetry_factors: &[(1, 824)], content: 0xc1a748638be6f398, numerators: 0x0f8c863152a43152, named: 0xb5205b5cda3bf8ff, ids: 0xc3905035e4cb6566 };
    fast xs_qed_eemm_identical_L2: Row { diagrams: 29, topologies: 157, assignments: 714, inverse_symmetry: (29, 1), weights: (4, 1), symmetry_factors: &[(1, 29)], content: 0xc6e7d902ad8d316f, numerators: 0x00fac78da0e0d159, named: 0x3d72c09c50112e41, ids: 0xcd1abf37940d48cf };
    fast xs_qed_eemm_sign_L2: Row { diagrams: 25, topologies: 157, assignments: 714, inverse_symmetry: (25, 1), weights: (8, 1), symmetry_factors: &[(1, 25)], content: 0x76af0acc0f7b7896, numerators: 0x5d7e712797f78751, named: 0x54e6507e4c865b62, ids: 0x52893e0efe85f881 };
    fast xs_qed_eemm_scalar_L2: Row { diagrams: 25, topologies: 157, assignments: 714, inverse_symmetry: (25, 1), weights: (8, 1), symmetry_factors: &[(1, 25)], content: 0x76af0acc0f7b7896, numerators: 0x5d7e712797f78751, named: 0x54e6507e4c865b62, ids: 0x52893e0efe85f881 };
    fast xs_qed_eedd_L1: Row { diagrams: 1, topologies: 17, assignments: 34, inverse_symmetry: (1, 1), weights: (-1, 1), symmetry_factors: &[(1, 1)], content: 0xb467d6be7a12aefa, numerators: 0xd3ad2f3a0802b6ab, named: 0xdd30b2e206f30267, ids: 0x58ebcb07ba6255d3 };
    fast xs_qed_eedd_L2: Row { diagrams: 30, topologies: 157, assignments: 714, inverse_symmetry: (30, 1), weights: (4, 1), symmetry_factors: &[(1, 30)], content: 0x1f065e0a67f8e4b1, numerators: 0x52854dced3f2ee5c, named: 0x368a1cfe1b96c3be, ids: 0xfba0be1b20d86887 };
    fast xs_qcdqed_addx_canonical_scalar_lr_L3: Row { diagrams: 82, topologies: 25, assignments: 1410, inverse_symmetry: (161, 2), weights: (-2711, 2), symmetry_factors: &[(1, 79), (2, 3)], content: 0x72adf523d45d2a92, numerators: 0x82a821cc454caa0c, named: 0x85a0b3b8988b3a99, ids: 0x84812f02ee345665 };
    fast xs_qed_eeaa_L2: Row { diagrams: 96, topologies: 157, assignments: 714, inverse_symmetry: (96, 1), weights: (-48, 1), symmetry_factors: &[(1, 96)], content: 0x8f8d4c480eadbae0, numerators: 0x8b2fe92ed5414355, named: 0xe67ffcbf30548e8b, ids: 0x9d68f748335a7cfd };
    fast xs_qed_aaee_L2: Row { diagrams: 248, topologies: 75, assignments: 3456, inverse_symmetry: (248, 1), weights: (176, 1), symmetry_factors: &[(1, 248)], content: 0xaa83568f8a1e3de4, numerators: 0xaa63603383de7b0b, named: 0x26940e2744fc7077, ids: 0x0ccc891930b8d526 };
    fast xs_qed_aaee_sym_L2: Row { diagrams: 124, topologies: 75, assignments: 3456, inverse_symmetry: (124, 1), weights: (176, 1), symmetry_factors: &[(1, 124)], content: 0x534f9729e6e7f0f1, numerators: 0x97b803df9beb77af, named: 0xb0209cd25cac901b, ids: 0x91e367b69f3f8408 };
    fast xs_qed_aaee_lr_L2: Row { diagrams: 130, topologies: 75, assignments: 3456, inverse_symmetry: (130, 1), weights: (176, 1), symmetry_factors: &[(1, 130)], content: 0xf55cebd45bdc155f, numerators: 0x9491ef44c58a58d8, named: 0x11544ab27bdc04bc, ids: 0x8e6709d74b6d01a8 };
    fast xs_qed_aaee_sym_lr_L2: Row { diagrams: 82, topologies: 75, assignments: 3456, inverse_symmetry: (82, 1), weights: (176, 1), symmetry_factors: &[(1, 82)], content: 0xf7997c216fc327b5, numerators: 0x42503dbf7a2c4825, named: 0x42f4a73ad2493141, ids: 0x3c1027b575a044ae };
    fast xs_qed_eeee_L1: Row { diagrams: 4, topologies: 17, assignments: 34, inverse_symmetry: (4, 1), weights: (0, 1), symmetry_factors: &[(1, 4)], content: 0x29f20710d0fbfbd6, numerators: 0x8eb2c3186607a1a2, named: 0xba552fc69e567351, ids: 0xfe0da547942d4483 };
    fast xs_qed_eeee_L2: Row { diagrams: 120, topologies: 157, assignments: 714, inverse_symmetry: (120, 1), weights: (0, 1), symmetry_factors: &[(1, 120)], content: 0x6442392293ae6a3a, numerators: 0x22a2c638353bfd5d, named: 0xfd42901c0c72c5c4, ids: 0x3cd3ebae4a2c6d84 };
    fast xs_qed_eeee_1pi_L1: Row { diagrams: 1, topologies: 3, assignments: 4, inverse_symmetry: (1, 1), weights: (-1, 1), symmetry_factors: &[(1, 1)], content: 0x8686eb17d2b0f293, numerators: 0x1fa84d9ad4fdf0b0, named: 0xf5ac9fd01b9eabd7, ids: 0xf777e735b39d93bf };
    fast xs_qed_eeee_1pi_L2: Row { diagrams: 31, topologies: 36, assignments: 116, inverse_symmetry: (31, 1), weights: (13, 1), symmetry_factors: &[(1, 31)], content: 0xa5f49d098f205eb2, numerators: 0xebf94f7f1a66762b, named: 0x5d7b34e127921dbe, ids: 0x35c2245f6a859c80 };
    fast xs_sm_eemm_L1: Row { diagrams: 16, topologies: 17, assignments: 928, inverse_symmetry: (16, 1), weights: (-16, 1), symmetry_factors: &[(1, 16)], content: 0x9ad7451811ae86e7, numerators: 0x4dbc4e60854caace, named: 0xe581437c7ccd8710, ids: 0x27354d61e6800e7c };
    fast xs_sm_eemm_L2: Row { diagrams: 3392, topologies: 230, assignments: 163172, inverse_symmetry: (3380, 1), weights: (-628, 1), symmetry_factors: &[(1, 3368), (2, 24)], content: 0x249c90853fb49a10, numerators: 0x7d90c9b07d93bfee, named: 0x9af22f71b53a3635, ids: 0x232d14193f0d1b4d };
    fast xs_sm_uudd_L1: Row { diagrams: 49, topologies: 17, assignments: 1652, inverse_symmetry: (49, 1), weights: (-9, 1), symmetry_factors: &[(1, 49)], content: 0x8d61c6602cd25421, numerators: 0x87213a78375152a0, named: 0x322b86a47378d1cd, ids: 0x581466731d283192 };
    fast xs_qcd_gggg_L1: Row { diagrams: 10, topologies: 55, assignments: 314, inverse_symmetry: (8, 1), weights: (8, 1), symmetry_factors: &[(1, 6), (2, 4)], content: 0x2748a62442b0c3cc, numerators: 0x2628c62dc373ac91, named: 0x3d5aa5fb16e7e40a, ids: 0xe2e7018b7e154d23 };
    fast xs_qcd_gggg_sym_L1: Row { diagrams: 10, topologies: 55, assignments: 314, inverse_symmetry: (8, 1), weights: (8, 1), symmetry_factors: &[(1, 6), (2, 4)], content: 0x2748a62442b0c3cc, numerators: 0x2628c62dc373ac91, named: 0x3d5aa5fb16e7e40a, ids: 0xe2e7018b7e154d23 };
    fast xs_qcd_gggg_L2: Row { diagrams: 1228, topologies: 840, assignments: 15253, inverse_symmetry: (1040, 1), weights: (-808, 1), symmetry_factors: &[(1, 868), (2, 328), (4, 32)], content: 0x471640cad31e9dfa, numerators: 0x6adf3bbba0e859ef, named: 0xa7634897a8bd0ea8, ids: 0x09059151749ab3a4 };
    fast xs_qcd_gguu_L1: Row { diagrams: 9, topologies: 55, assignments: 314, inverse_symmetry: (9, 1), weights: (-9, 1), symmetry_factors: &[(1, 9)], content: 0xc6a7c47bd2905b4d, numerators: 0x5edf8e5f7a11eeb3, named: 0xeadb4aa521f5d4c9, ids: 0x2ddf4705c9d52678 };
    fast xs_qcd_gguu_L2: Row { diagrams: 558, topologies: 840, assignments: 15253, inverse_symmetry: (528, 1), weights: (210, 1), symmetry_factors: &[(1, 498), (2, 60)], content: 0x58c4de67809aa876, numerators: 0xfb58d66eb84d6364, named: 0x81138e2cd8849006, ids: 0xef44682c6ce8403a };
    fast xs_qcd_gguu_sym_L1: Row { diagrams: 5, topologies: 55, assignments: 314, inverse_symmetry: (5, 1), weights: (-9, 1), symmetry_factors: &[(1, 5)], content: 0x09e371568eb58f32, numerators: 0xb9fdc23ed056fa0c, named: 0xa1706e5ab99b8939, ids: 0xde71f8568766e9f4 };
    fast xs_qcd_gguu_1pi_L1: Row { diagrams: 4, topologies: 12, assignments: 54, inverse_symmetry: (4, 1), weights: (-4, 1), symmetry_factors: &[(1, 4)], content: 0xb8787afbc1ba43bd, numerators: 0x75fbcf2f108240fb, named: 0xd40765b02a1bdf5a, ids: 0x822a3b576554f53b };
    fast sm_uudd_qcd2_qed2_L1: Row { diagrams: 168, topologies: 17, assignments: 168, inverse_symmetry: (168, 1), weights: (-32, 1), symmetry_factors: &[(1, 168)], content: 0x3ceb37034b64690e, numerators: 0x0f525334a06eeca6, named: 0xa2461904ea007f6b, ids: 0xb151dae9858525cc };
    fast sm_uudd_qcd3_L1: Row { diagrams: 18, topologies: 17, assignments: 18, inverse_symmetry: (35, 2), weights: (-7, 2), symmetry_factors: &[(1, 17), (2, 1)], content: 0xf68ee61e14441815, numerators: 0x5b5cbbd732e9fb67, named: 0x75d501c4462bb8d6, ids: 0x6540c51587a5e3b1 };
    fast sm_eemm_qed4_L0to1: Row { diagrams: 464, topologies: 10, assignments: 464, inverse_symmetry: (925, 2), weights: (-205, 2), symmetry_factors: &[(1, 461), (2, 3)], content: 0x7d459f5367fbcba8, numerators: 0xdf611964612f3ec9, named: 0x7b2f048149dccbf3, ids: 0xa1e6a4f3cf7051f9 };
    fast sm_eemm_floop0_L1: Row { diagrams: 317, topologies: 10, assignments: 317, inverse_symmetry: (631, 2), weights: (-499, 2), symmetry_factors: &[(1, 314), (2, 3)], content: 0xae30817ce19acfb8, numerators: 0x8837935e332cfed2, named: 0x8bd58d385ca26fcc, ids: 0x9daa1bb29c9c73ae };
    fast qcd_gggg_floop1_L1: Row { diagrams: 222, topologies: 55, assignments: 222, inverse_symmetry: (222, 1), weights: (-222, 1), symmetry_factors: &[(1, 222)], content: 0xb72175978408ff0d, numerators: 0x7762bcfb04e9fa1f, named: 0x978d5636f7e2a170, ids: 0x0c9a9ab09256e03e };
    fast qcd_gse_all_floop1_dfilter_L2: Row { diagrams: 55, topologies: 8, assignments: 55, inverse_symmetry: (95, 2), weights: (-55, 2), symmetry_factors: &[(1, 40), (2, 15)], content: 0x82e38ea40c358adb, numerators: 0xb8f71abe8d5dc95e, named: 0xe395dde52cbaa950, ids: 0x5f2bb4d18fa25f53 };
    fast sm_eemm_floop1_tfilter_L1: Row { diagrams: 131, topologies: 10, assignments: 131, inverse_symmetry: (131, 1), weights: (131, 1), symmetry_factors: &[(1, 131)], content: 0xe429b2e86ed27f62, numerators: 0x6ba3f3278f2e79ab, named: 0x63d9b0ffb978cddf, ids: 0x7163d47096d597d4 };
    fast qed_ase_floop2_L3: Row { diagrams: 405, topologies: 10, assignments: 405, inverse_symmetry: (405, 1), weights: (405, 1), symmetry_factors: &[(1, 405)], content: 0x5d24a56e8a9ba575, numerators: 0x43db9f7cc9f679f7, named: 0x109b995d1525e3de, ids: 0x6daf1e85e33ce62c };
    fast ctorder_qcd_gggg_L1: Row { diagrams: 544, topologies: 55, assignments: 544, inverse_symmetry: (492, 1), weights: (-292, 1), symmetry_factors: &[(1, 440), (2, 104)], content: 0x7c6df59b5530d1bb, numerators: 0x1b18ddafcee9c3d2, named: 0xffc8351857869ba7, ids: 0xeb5e38c9891720f9 };
    fast ctorder_qcd_uugg_L1: Row { diagrams: 268, topologies: 32, assignments: 268, inverse_symmetry: (251, 1), weights: (13, 1), symmetry_factors: &[(1, 234), (2, 34)], content: 0x33e52219b458ee2e, numerators: 0x3a743bf033456b6f, named: 0xc760bf1c8a759334, ids: 0x019e715c187af508 };
    fast ctorder_qcd_vac_L3: Row { diagrams: 1298, topologies: 12, assignments: 1298, inverse_symmetry: (18965, 24), weights: (-3859, 24), symmetry_factors: &[(1, 431), (2, 601), (3, 14), (4, 177), (6, 30), (8, 33), (12, 3), (16, 5), (24, 1), (48, 3)], content: 0xcc01ea12aa514317, numerators: 0x5fd0f29cbe24355f, named: 0x1912ce16ab5de4c6, ids: 0x74aeeea5e0eb56f4 };
    fast ctorder_phi3_vac_L4: Row { diagrams: 117, topologies: 17, assignments: 117, inverse_symmetry: (5525, 384), weights: (5525, 384), symmetry_factors: &[(2, 5), (4, 17), (8, 37), (12, 1), (16, 38), (32, 17), (64, 1), (128, 1)], content: 0x97fa220e063050c0, numerators: 0x1cb20ad83919649e, named: 0x9842874402915d71, ids: 0x7b8572ac486071f1 };
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

/// Fresh processes at each thread count register model and index owners in distinct orders.
/// Compare the original IDs and the complete graph, numerator and grouping payloads.
#[test]
fn ids_do_not_depend_on_threads_or_symbol_registration() {
    const CHILD: &str = "FEYNKIT_ORACLE_ID_CHILD";
    const RUN: [&str; 3] = ["qcd_uugg_L1", "sm_eemm_L1", "xs_qed_eemm_L2"];
    if let Ok(child) = env::var(CHILD) {
        let (threads, registration) = child.split_once(',').unwrap_or((&child, "ordinary"));
        // These names have no global initializer. Assert the actual registration order,
        // rather than assuming a preload changed the process-global symbol table.
        if registration == "reverse" {
            index_symbol!("feynkit_oracle_registration::b");
            index_symbol!("feynkit_oracle_registration::a");
            tensor_symbol!("feynkit_oracle_registration::B");
            tensor_symbol!("feynkit_oracle_registration::A");
        } else {
            index_symbol!("feynkit_oracle_registration::a");
            index_symbol!("feynkit_oracle_registration::b");
            tensor_symbol!("feynkit_oracle_registration::A");
            tensor_symbol!("feynkit_oracle_registration::B");
        }
        let a = index_symbol!("feynkit_oracle_registration::a");
        let b = index_symbol!("feynkit_oracle_registration::b");
        let first = tensor_symbol!("feynkit_oracle_registration::A");
        let second = tensor_symbol!("feynkit_oracle_registration::B");
        assert_eq!(a.get_id() < b.get_id(), registration != "reverse");
        assert_eq!(first.get_id() < second.get_id(), registration != "reverse");
        let vector = |head, index| {
            FunctionBuilder::new(head)
                .add_arg(Minkowski {}.new_rep(4).to_symbolic([Atom::var(index)]))
                .finish()
        };
        let tensor = NumeratorInterface::admit(vector(first, a) * vector(second, b), None).unwrap();
        let (expression, interface) = NumeratorInterface::export(&tensor);
        println!(
            "{CHILD} admitted {} {:?}",
            expression.to_canonical_string(),
            interface.map(|value| value.to_json())
        );
        if registration != "ordinary" {
            let mut owners = [
                symbols::hedge_index,
                symbols::vertex_index,
                symbols::edge_index,
            ];
            let mut models = [
                Model::phi_3_4,
                Model::yang_mills,
                Model::electroweak,
                Model::qcd_qed,
            ];
            if registration == "reverse" {
                owners.reverse();
                models.reverse();
            } else {
                assert_eq!(registration, "preload");
            }
            for owner in owners {
                owner();
            }
            for model in models {
                model();
            }
        }
        for name in RUN {
            let built = grouped(Case::named(name), threads.parse().unwrap());
            for diagram in built.iter() {
                println!(
                    "{CHILD} {name} {} {}",
                    diagram.name().expect("named oracle output"),
                    diagram.id().unwrap()
                );
            }
            let mut exact = Vec::new();
            Census::of(&built).write_exact_dump(&mut exact).unwrap();
            println!("{CHILD} {name} full {}", String::from_utf8(exact).unwrap());
        }
        return;
    }
    let test = "ids_do_not_depend_on_threads_or_symbol_registration";
    let reference = child_lines(test, CHILD, "1");
    assert!(reference.len() > 500);
    for threads in [1, 8, 32] {
        for registration in ["ordinary", "preload", "reverse"] {
            assert_eq!(
                child_lines(test, CHILD, &format!("{threads},{registration}")),
                reference,
                "threads={threads}, registration={registration}"
            );
        }
    }
}

/// Printed expressions follow Symbolica's symbol order, which must not depend on the generation
/// thread that first meets a coupling or an overall-factor head: fresh processes at 1 and at 32
/// threads print every numerator and overall factor identically. The test runs its own binary as
/// those processes.
#[test]
fn printed_numerators_do_not_depend_on_threads() {
    const CHILD: &str = "FEYNKIT_ORACLE_PRINT_CHILD";
    const RUN: [&str; 4] = [
        "qcd_uugg_L1",
        "sm_eemm_L1",
        "ym_ggg_zero_color_L3",
        "xs_qcdqed_addx_canonical_scalar_lr_L3",
    ];
    if let Ok(threads) = env::var(CHILD) {
        let threads = threads.parse().unwrap();
        for name in RUN {
            print_outputs(CHILD, name, &grouped(Case::named(name), threads));
        }
        print_outputs(
            CHILD,
            "ym_ggg_sampled_sign_subset_L3",
            &sampled_color_subset(threads),
        );
        return;
    }
    let test = "printed_numerators_do_not_depend_on_threads";
    let reference = child_lines(test, CHILD, "1");
    assert!(reference.len() > 500);
    for threads in [1, 8, 32] {
        assert_child_lines_eq(test, CHILD, &threads.to_string(), &reference);
    }
}

/// Cook the tensor-index payloads `owner(n, 1)` of the half-edge, edge and vertex owners for
/// `n < 64` in `order` (`generation` cooks none) and return the number of symbols registered
/// afterwards.
fn register_indices(order: &str) -> usize {
    let mut payloads = [
        symbols::hedge_index(),
        symbols::edge_index(),
        symbols::vertex_index(),
    ]
    .into_iter()
    .flat_map(|owner| (0..64).map(move |n| function!(owner, n, 1)))
    .collect::<Vec<_>>();
    match order {
        "generation" => payloads.clear(),
        "increasing" => {}
        "decreasing" => payloads.reverse(),
        "shuffled" => {
            // Fisher-Yates with a fixed linear congruential sequence.
            let mut state = 0x2545_f491_4f6c_dd1d_u64;
            for last in (1..payloads.len()).rev() {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                payloads.swap(last, (state >> 33) as usize % (last + 1));
            }
        }
        _ => panic!("unknown registration order {order}"),
    }
    // Native index cooking deliberately preserves Named indices after #209.
    // Explicit function encoding exercises the supported reversible aliases instead.
    let cooking = CookSettings::indices()
        .with_source_filter(CookSourceFilter::AnyFunction)
        .with_mode(CookMode::ReversibleEncoding);
    let mut registered = BTreeMap::new();
    for payload in &payloads {
        let AtomView::Fun(payload) = payload.as_view() else {
            unreachable!("index payloads are functions")
        };
        let cooked = cooking.cook_function(payload).unwrap();
        let AtomView::Var(variable) = cooked.as_view() else {
            panic!("cooked index must be a variable")
        };
        registered.insert(variable.get_symbol(), payload.as_view().to_owned());
        let head = tensor_symbol!("feynkit_oracle_registration::cooked_probe");
        let tensor = NumeratorInterface::admit(
            function!(head, Minkowski {}.new_rep(4).to_symbolic([cooked])),
            None,
        )
        .unwrap();
        assert_eq!(tensor.rank(), 1);
        assert_eq!(
            NumeratorInterface::export(&tensor).0,
            function!(
                head,
                Minkowski {}
                    .new_rep(4)
                    .to_symbolic([payload.as_view().to_owned()])
            )
        );
    }
    assert_eq!(
        registered.values().collect::<Vec<_>>(),
        payloads.iter().collect::<Vec<_>>(),
        "actual cooked index registration order"
    );
    State::symbol_iter().count()
}

/// Print the process-independent outputs of `name`, each line prefixed by
/// `prefix`: the zero-color count, every diagram's exported numerator and overall factor, and
/// every group's members with their ratios.
fn print_outputs(prefix: &str, name: &str, result: &BuiltDiagramSet) {
    println!(
        "{prefix} {name} zero {}",
        result.report().zero_numerator_count
    );
    for diagram in result.iter() {
        let (numerator, _) = NumeratorInterface::export(diagram.numerator().unwrap());
        println!(
            "{prefix} {name} {} {numerator} {}",
            diagram.name().expect("named oracle output"),
            diagram.overall_factor()
        );
    }
    let mut exact = Vec::new();
    Census::of(result).write_exact_dump(&mut exact).unwrap();
    println!("{prefix} {name} full {}", String::from_utf8(exact).unwrap());
    for group in result.groups() {
        let members = group
            .members
            .iter()
            .map(|member| format!("{}:{}", member.source_diagram, member.ratio))
            .collect::<Vec<_>>();
        println!("{prefix} {name} group {}", members.join(" "));
    }
}

/// A small witness from the full YM fixture keeps sampled-grouping registration checks cheap.
/// D0096, D0194, D0196, D0201, D0203, D0480 and D0492 vanish by color algebra; the old grouped
/// zero test misses some of them under permuted registration. D0000/D0001 exercise sampling of
/// survivors. The standalone filter is off, so only grouping can remove these seven zeros.
fn sampled_color_subset(threads: usize) -> BuiltDiagramSet {
    let case = Case::named("ym_ggg_zero_color_L3");
    let model = case.model();
    let run = Run::default().threads(threads);
    let built = case
        .request(&model)
        .generate(model.into(), &run)
        .unwrap()
        .build(
            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
            &run,
        )
        .unwrap()
        .select(
            &Selection::new().select_names([
                "D0000", "D0001", "D0096", "D0194", "D0196", "D0201", "D0203", "D0480", "D0492",
            ]),
            &run,
        )
        .unwrap();
    assert_eq!(built.len(), 9, "sampled color-zero witness selection");
    let result = built
        .group(&Grouping::UpToSign(GraphGroupingOptions::default()), &run)
        .unwrap();
    assert_eq!(result.report().zero_numerator_count.to_string(), "7");
    result
}

/// Symbolica orders terms by symbol registration, and the parallel build cooks the tensor indices
/// of half-edges, edges and vertices into symbols in the order its threads first meet them. The
/// zero-color proof (including sampled grouping with the filter off) and canonical numerator
/// grouping must not depend on that order: fresh
/// processes that register these indices first, in increasing, decreasing or shuffled order,
/// remove the same zero-color diagrams and form the same groups, with the same numerators and
/// ratios, as a process that registers them during generation. The test runs its own binary as
/// those processes.
fn check_index_registration_order(test: &str, order: &str, count: usize) {
    const CHILD: &str = "FEYNKIT_ORACLE_REGISTRATION_CHILD";
    const RUN: [&str; 2] = [
        "ym_ggg_zero_color_L3",
        "xs_qcdqed_addx_canonical_scalar_lr_L3",
    ];
    if let Ok(order) = env::var(CHILD) {
        let (order, count) = order.split_once(',').unwrap();
        let count = count.parse().unwrap();
        let registered = register_indices(order);
        for name in RUN {
            print_outputs(CHILD, name, &grouped(Case::named(name), count));
        }
        print_outputs(
            CHILD,
            "ym_ggg_sampled_sign_subset_L3",
            &sampled_color_subset(count),
        );
        // Generation must have met only pre-registered indices. Cooked index symbols keep their
        // payload; cooked dimensions encode it in their names.
        let cooked_later = State::symbol_iter()
            .skip(registered)
            .filter(|(_, name)| name.starts_with("idenso::cooked_"))
            .filter_map(|(symbol, _)| match symbol.get_data() {
                UserData::Atom(payload) => Some(payload.to_canonical_string()),
                _ => None,
            })
            .collect::<Vec<_>>();
        if order != "generation" {
            assert!(
                cooked_later.is_empty(),
                "{order}: generation cooked unregistered indices {cooked_later:?}"
            );
        }
        return;
    }
    // Each cell captures the same semantic reference in an independent fresh process.
    let reference = child_lines(test, CHILD, "generation,1");
    assert!(reference.len() > 100);
    assert_child_lines_eq(test, CHILD, &format!("{order},{count}"), &reference);
}

// One reference and one actual child per test fit the normal Nextest budget.
macro_rules! registration_cells {
    ($($name:ident: $order:literal, $count:literal;)*) => {
        $(
            #[test]
            fn $name() {
                check_index_registration_order(stringify!($name), $order, $count);
            }
        )*
    };
}

registration_cells! {
    outputs_do_not_depend_on_index_registration_order_generation_1: "generation", 1;
    outputs_do_not_depend_on_index_registration_order_increasing_1: "increasing", 1;
    outputs_do_not_depend_on_index_registration_order_decreasing_1: "decreasing", 1;
    outputs_do_not_depend_on_index_registration_order_shuffled_1: "shuffled", 1;
    outputs_do_not_depend_on_index_registration_order_generation_8: "generation", 8;
    outputs_do_not_depend_on_index_registration_order_increasing_8: "increasing", 8;
    outputs_do_not_depend_on_index_registration_order_decreasing_8: "decreasing", 8;
    outputs_do_not_depend_on_index_registration_order_shuffled_8: "shuffled", 8;
    outputs_do_not_depend_on_index_registration_order_generation_32: "generation", 32;
    outputs_do_not_depend_on_index_registration_order_increasing_32: "increasing", 32;
    outputs_do_not_depend_on_index_registration_order_decreasing_32: "decreasing", 32;
    outputs_do_not_depend_on_index_registration_order_shuffled_32: "shuffled", 32;
}

/// A child panic must fail the parent comparison, even before any payload is printed.
#[test]
#[should_panic(expected = "child")]
fn index_registration_child_failure_propagates() {
    child_lines(
        "outputs_do_not_depend_on_index_registration_order_generation_1",
        "FEYNKIT_ORACLE_REGISTRATION_CHILD",
        "generation,invalid-thread-count",
    );
}

/// The named and grouped diagrams of `case` at `threads` threads.
fn grouped(case: &Case, threads: usize) -> BuiltDiagramSet {
    let model = case.model();
    let observed = Arc::new(AtomicUsize::new(0));
    let sink = Arc::clone(&observed);
    let caller = thread::current().id();
    let run = Run::default().threads(threads).progress(move |progress| {
        if matches!(
            progress.stage,
            Stage::ZeroColor | Stage::Numerators | Stage::GroupingPreparation
        ) && progress.completed > 0usize.into()
        {
            if threads == 1 {
                assert_eq!(thread::current().id(), caller);
                assert!(rayon::current_thread_index().is_none());
            } else {
                assert_eq!(rayon::current_num_threads(), threads);
            }
            sink.fetch_add(1, Ordering::Relaxed);
        }
        GenerationControl::Continue
    });
    let result = case
        .request(&model)
        .generate(model.into(), &run)
        .unwrap()
        .build(&case.numerators().with_naming(Naming::prefix("D")), &run)
        .unwrap()
        .group(&case.grouping(), &run)
        .unwrap();
    assert!(
        observed.load(Ordering::Relaxed) > 0,
        "actual worker pool was not observed"
    );
    result
}

/// Run `test` of this binary in a fresh process with `variable=value` and return the lines it
/// prints that start with `variable`.
fn child_lines(test: &str, variable: &str, value: &str) -> Vec<String> {
    let output = Command::new(env::current_exe().unwrap())
        .args(["--exact", test, "--nocapture"])
        .env(variable, value)
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        output.status.success(),
        "{value}: child {}: stdout/stderr (first4000 characters each):\n{}\n{}",
        output.status,
        stdout.chars().take(4000).collect::<String>(),
        String::from_utf8_lossy(&output.stderr)
            .chars()
            .take(4000)
            .collect::<String>()
    );
    stdout
        .lines()
        .filter(|line| line.starts_with(variable))
        .map(str::to_owned)
        .collect()
}

/// Report the first differing line without dumping all exported numerators.
fn assert_child_lines_eq(test: &str, variable: &str, value: &str, reference: &[String]) {
    let actual = child_lines(test, variable, value);
    for (index, (actual, expected)) in actual.iter().zip(reference).enumerate() {
        assert!(
            actual == expected,
            "{value}: line {index} differs (first 240 characters):\nactual: {}\nexpected: {}",
            actual.chars().take(240).collect::<String>(),
            expected.chars().take(240).collect::<String>()
        );
    }
    assert_eq!(actual.len(), reference.len(), "{value}: output line count");
}

/// A user filter pruning in flight, a particle veto and post-filtering the unfiltered diagrams
/// select the same ghost-free diagrams, with the same symmetry and overall factors, cuts and
/// numerators; only the names of the unfiltered run differ.
#[test]
fn in_flight_filters_select_like_vetoes_and_post_filtering() {
    let case = Case::named("qcd_gggg_L1");
    let model = case.model();
    let run = Run::default().threads(1);
    let all = case
        .request(&model)
        .generate(model.clone().into(), &run)
        .unwrap()
        .build(
            &NumeratorOptions::default().with_naming(Naming::prefix("D")),
            &run,
        )
        .unwrap()
        .group(&case.grouping(), &run)
        .unwrap();
    let ghost_free: BTreeSet<_> = all
        .iter()
        .filter(|diagram| !is_ghost(&model, diagram))
        .map(|diagram| diagram.id().unwrap())
        .collect();
    assert!(ghost_free.len() < all.len());
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
        let run = Run::default().threads(1);
        let mut expected = Rational::zero();
        let generated = case
            .request(&model)
            .generate(model.clone().into(), &run)
            .unwrap()
            .build(
                &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                &run,
            )
            .unwrap()
            .group(&case.grouping(), &run)
            .unwrap();
        for diagram in generated.iter() {
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
        let run = Run::default().threads(1);
        let binomial = |m: u64, k: u64| (0..k).fold(1, |c, i| c * (m - i) / (i + 1));
        let mut expected = Rational::zero();
        let generated = case
            .request(&model)
            .generate(model.clone().into(), &run)
            .unwrap()
            .build(
                &NumeratorOptions::default().with_naming(Naming::prefix("D")),
                &run,
            )
            .unwrap()
            .group(&case.grouping(), &run)
            .unwrap();
        for diagram in generated.iter() {
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
                expected += Rational::from((u128::from(binomial(m, k)), diagram.symmetry_factor()));
            }
        }
        let actual = census(counterterms, 1).summary().inverse_symmetry;
        assert_eq!(actual, expected, "{}", counterterms.name);
    }
}
