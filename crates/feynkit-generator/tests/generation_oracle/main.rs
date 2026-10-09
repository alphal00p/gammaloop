//! Differential oracles for diagram generation.
//!
//! Every case of [`cases::CASES`] runs at 1, 8 and 32 threads. The runs must agree exactly (diagram
//! set, |Aut|, overall factors, frames, numerators, ids, groups, names and output order), and their
//! [`census::Summary`] must equal the case's row below. The `phi` rows are fixed by independent
//! oracles: OEIS A005967 counts the connected cubic vacuum graphs, zero-dimensional field theory
//! gives every sum of 1/|Aut|, and FeynGraph and Graphica agree with every count and |Aut| multiset
//! (Graphica alone for the phi^3 vacua at six and seven loops, where FeynGraph did not finish). The
//! other rows were recorded with this harness; their counts reproduce the profile of the native
//! generator (e.g. 15,253 and 686,501 diagrams for gg->gg at two and three loops). `failing` cases
//! fail on the base: each must fail at every thread count with an error that contains the recorded
//! fragments, so that a change neither hides nor fixes one silently; a fix records the case's row
//! instead. The fragments leave out diagram and candidate numbers, because parallel runs report
//! whichever failure they meet first.
//!
//! Heavy cases are ignored by default; `just generation-oracle-heavy` runs them. Environment:
//! `FEYNKIT_ORACLE_THREADS=1,8,32` sets the thread counts, `FEYNKIT_ORACLE_RECORD=1` prints each
//! case's row instead of comparing it, and `FEYNKIT_ORACLE_DUMP=<dir>` writes every census as
//! text, for diffing the diagrams of two builds.

mod cases;
mod census;

use std::collections::BTreeSet;

use cases::{CASES, Case, REFERENCE};
use census::{Census, Scope, Summary};
use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::domains::rational::Rational;

fn threads() -> Vec<usize> {
    std::env::var("FEYNKIT_ORACLE_THREADS").map_or(vec![1, 8, 32], |threads| {
        threads.split(',').map(|t| t.parse().unwrap()).collect()
    })
}

fn census(case: &Case, threads: usize) -> Census {
    let model = case.model();
    let options = case.options(&model, threads);
    let result = case
        .generate_with(model, &options)
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
    assert!(result.report.completed, "{} did not complete", case.name);
    Census::of(&result)
}

/// A recorded [`Summary`], with the sum of 1/|Aut| as numerator and denominator in lowest terms.
struct Row {
    diagrams: usize,
    topologies: usize,
    assignments: usize,
    inverse_symmetry: (u64, u64),
    symmetry_factors: &'static [(u64, usize)],
    content: u64,
    numerators: u64,
    named: u64,
}

impl From<&Row> for Summary {
    fn from(row: &Row) -> Summary {
        Summary {
            diagrams: row.diagrams,
            topologies: row.topologies,
            assignments: row.assignments,
            inverse_symmetry: Rational::from(row.inverse_symmetry),
            symmetry_factors: row.symmetry_factors.to_vec(),
            content: row.content,
            numerators: row.numerators,
            named: row.named,
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
    for count in threads {
        let other = census(case, count);
        reference.assert_same(
            &other,
            Scope::Exact,
            &format!("{name}: {first} vs {count} threads"),
        );
    }
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
/// Cross sections at two loops with gluons, and at three loops in QED: a threshold candidate's
/// cut does not match its crossing edges.
const THRESHOLD_CANDIDATE: &[&str] = &[
    "invalid topology threshold candidate ",
    ": oriented cut endpoints do not equal the crossing edges",
];
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
    for threads in threads() {
        let model = case.model();
        let options = case.options(&model, threads);
        match case.generate_with(model, &options) {
            Ok(_) => panic!("{name} generates at {threads} threads: record its oracle row"),
            Err(found) => {
                let found = found.to_string();
                assert!(
                    error.iter().all(|fragment| found.contains(fragment)),
                    "{name} at {threads} threads fails with {found:?}, not {error:?}"
                );
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
    fast phi3_vac_L2: Row { diagrams: 2, topologies: 2, assignments: 2, inverse_symmetry: (5, 24), symmetry_factors: &[(8, 1), (12, 1)], content: 0x75492f43d8fc3c82, numerators: 0xf818ded4decf2e81, named: 0x066c36fbb2e7cba6 };
    fast phi3_vac_L3: Row { diagrams: 5, topologies: 5, assignments: 5, inverse_symmetry: (5, 16), symmetry_factors: &[(8, 1), (16, 2), (24, 1), (48, 1)], content: 0x54ca4e63ab0c62fd, numerators: 0xaa97877126478b1a, named: 0xe1461100e440204c };
    fast phi3_vac_L4: Row { diagrams: 17, topologies: 17, assignments: 17, inverse_symmetry: (1105, 1152), symmetry_factors: &[(8, 3), (12, 1), (16, 5), (32, 4), (48, 2), (72, 1), (128, 1)], content: 0x7b4cd459d7bc1a0c, numerators: 0x7d649ade64818e41, named: 0x3c02c2a6679da767 };
    fast phi3_vac_L5: Row { diagrams: 71, topologies: 71, assignments: 71, inverse_symmetry: (565, 128), symmetry_factors: &[(4, 4), (8, 10), (12, 1), (16, 21), (32, 17), (48, 1), (64, 11), (96, 2), (128, 2), (256, 2)], content: 0x78f7d529862a12a8, numerators: 0x6d848bd898887500, named: 0x160d29fbeb1b8641 };
    fast phi3_vac_L6: Row { diagrams: 388, topologies: 388, assignments: 388, inverse_symmetry: (82825, 3072), symmetry_factors: &[(2, 6), (4, 25), (6, 2), (8, 60), (12, 1), (16, 88), (20, 2), (32, 89), (48, 5), (64, 67), (120, 1), (128, 27), (256, 6), (320, 2), (384, 2), (512, 4), (3072, 1)], content: 0xe012291e41427b03, numerators: 0xccd2d268f51c99e6, named: 0xcf13baf49af861a8 };
    fast phi3_vac_L7: Row { diagrams: 2592, topologies: 2592, assignments: 2592, inverse_symmetry: (19675, 96), symmetry_factors: &[(1, 5), (2, 70), (4, 219), (8, 408), (12, 4), (16, 529), (18, 1), (24, 3), (32, 524), (36, 1), (48, 6), (64, 425), (96, 6), (128, 246), (256, 103), (384, 3), (512, 22), (768, 2), (1024, 12), (2048, 2), (3072, 1)], content: 0x08a76811f8931455, numerators: 0x64b6e52a55ca9139, named: 0xd53cab502ff0b2a1 };
    fast phi4_vac_L2: Row { diagrams: 1, topologies: 1, assignments: 1, inverse_symmetry: (1, 8), symmetry_factors: &[(8, 1)], content: 0xaaf29483b37dfc2e, numerators: 0x9ee3522ae46a98c7, named: 0x9b4091ba687c38c2 };
    fast phi4_vac_L3: Row { diagrams: 2, topologies: 2, assignments: 2, inverse_symmetry: (1, 12), symmetry_factors: &[(16, 1), (48, 1)], content: 0x860fea0fe2959078, numerators: 0xf1dbbf34c9c784d2, named: 0xa22bffc660db5808 };
    fast phi4_vac_L4: Row { diagrams: 4, topologies: 4, assignments: 4, inverse_symmetry: (11, 96), symmetry_factors: &[(24, 1), (32, 1), (48, 2)], content: 0x7c9a015c48ce66f1, numerators: 0x89321207c7e09198, named: 0x2f8d97c1e73f9fc1 };
    fast phi4_vac_L5: Row { diagrams: 10, topologies: 10, assignments: 10, inverse_symmetry: (17, 72), symmetry_factors: &[(16, 1), (32, 3), (48, 2), (64, 1), (128, 2), (144, 1)], content: 0xc84d3a9a80a7e93c, numerators: 0xfb17c8e5e01107ce, named: 0xdf5b01c4c739d6fc };
    fast phi4_vac_L6: Row { diagrams: 28, topologies: 28, assignments: 28, inverse_symmetry: (619, 960), symmetry_factors: &[(16, 4), (32, 6), (48, 3), (64, 4), (96, 4), (120, 1), (128, 2), (144, 1), (288, 1), (320, 2)], content: 0xaa7cbaacc8b958dc, numerators: 0x2712de1422c359eb, named: 0x20b3e2dc1c9e938a };
    fast phi4_vac_L7: Row { diagrams: 97, topologies: 97, assignments: 97, inverse_symmetry: (709, 324), symmetry_factors: &[(8, 4), (16, 7), (24, 2), (32, 13), (48, 6), (64, 23), (96, 13), (128, 9), (144, 1), (192, 6), (256, 3), (288, 3), (384, 1), (576, 2), (768, 3), (1296, 1)], content: 0x6876ff8f305ca35f, numerators: 0x72c37493df5cad32, named: 0x66c604f0cfb9d187 };
    fast phi34_2pt_L1: Row { diagrams: 3, topologies: 3, assignments: 3, inverse_symmetry: (3, 2), symmetry_factors: &[(2, 3)], content: 0x78c9ff93a7a0c1ed, numerators: 0x954202b142bde800, named: 0xf7ef2e076d8426f3 };
    fast phi34_2pt_L2: Row { diagrams: 29, topologies: 29, assignments: 29, inverse_symmetry: (25, 3), symmetry_factors: &[(2, 6), (4, 19), (6, 2), (8, 2)], content: 0xcd7bebf282a389a7, numerators: 0x2e563a2edf198f53, named: 0x7110cf09185c128c };
    fast phi34_2pt_L3: Row { diagrams: 351, topologies: 351, assignments: 351, inverse_symmetry: (1741, 24), symmetry_factors: &[(1, 2), (2, 36), (4, 125), (6, 3), (8, 141), (12, 22), (16, 21), (48, 1)], content: 0xcad0369a6ad0e609, numerators: 0x9a49e9aaae10c880, named: 0xe32f4421dfa881a4 };
    fast phi34_2pt_L4: Row { diagrams: 5076, topologies: 5076, assignments: 5076, inverse_symmetry: (80299, 96), symmetry_factors: &[(1, 40), (2, 370), (4, 1180), (6, 13), (8, 1719), (12, 130), (16, 1157), (24, 216), (32, 211), (36, 5), (48, 13), (64, 5), (72, 2), (96, 13), (128, 2)], content: 0xb14dd9e709143bda, numerators: 0x84fab1bbe550b67a, named: 0xa08501feb637da25 };
    heavy phi34_2pt_L5: Row { diagrams: 84749, topologies: 84749, assignments: 84749, inverse_symmetry: (6869123, 576), symmetry_factors: &[(1, 727), (2, 5060), (4, 15053), (6, 145), (8, 23639), (12, 1092), (16, 21069), (24, 2619), (32, 10408), (36, 27), (48, 2146), (64, 2057), (72, 115), (96, 296), (128, 108), (144, 22), (192, 130), (256, 22), (288, 7), (384, 7)], content: 0xad6afc0df417f26b, numerators: 0x3687bef3749d414c, named: 0x6fe15fff1b3e21cf };
    fast phi3_2pt_L1: Row { diagrams: 1, topologies: 1, assignments: 1, inverse_symmetry: (1, 2), symmetry_factors: &[(2, 1)], content: 0x55e64ff2da4d1e0f, numerators: 0x316859ad4f714f0a, named: 0x284ca460b3985cfe };
    fast phi3_2pt_tad_L2: Row { diagrams: 10, topologies: 10, assignments: 10, inverse_symmetry: (25, 8), symmetry_factors: &[(2, 3), (4, 6), (8, 1)], content: 0x545153d6589cc2b7, numerators: 0x305722691bab7a7f, named: 0xc445443935c39384 };
    fast phi34_2pt_identical_L3: Row { diagrams: 72, topologies: 72, assignments: 72, inverse_symmetry: (71, 3), symmetry_factors: &[(1, 2), (2, 24), (4, 32), (6, 1), (8, 10), (12, 3)], content: 0x6454c5db351b672c, numerators: 0xe64b2d72fd7c2e4e, named: 0x4a42dbba4affe049 };
    fast qed_eemm_L1: Row { diagrams: 17, topologies: 10, assignments: 17, inverse_symmetry: (17, 1), symmetry_factors: &[(1, 17)], content: 0x67015758ee0dfcf5, numerators: 0x21ff907ef41e3790, named: 0x43028a9ad93aa749 };
    fast qed_eemm_L2: Row { diagrams: 357, topologies: 116, assignments: 357, inverse_symmetry: (357, 1), symmetry_factors: &[(1, 357)], content: 0xb41eb8293f39b714, numerators: 0x8459c2e5cff56b07, named: 0x346b21e92beaba31 };
    heavy qed_eemm_L3: Row { diagrams: 8724, topologies: 1245, assignments: 8724, inverse_symmetry: (8724, 1), symmetry_factors: &[(1, 8724)], content: 0xb687808b888a00cc, numerators: 0xee51106cc3fe4d1f, named: 0xa91e5b83ffa3ef2d };
    fast qed_ee4f_L0: Row { diagrams: 6, topologies: 7, assignments: 6, inverse_symmetry: (6, 1), symmetry_factors: &[(1, 6)], content: 0x66f2c9719dc08f84, numerators: 0x4efc7631fa9c125f, named: 0x3bd1e698ebc5dd09 };
    fast qed_ee4f_L1: Row { diagrams: 242, topologies: 304, assignments: 242, inverse_symmetry: (242, 1), symmetry_factors: &[(1, 242)], content: 0xbac90b32e6625e27, numerators: 0xb65507056583d4ab, named: 0xd83f8ba5b698b986 };
    heavy qed_ee4f_L2: Row { diagrams: 8460, topologies: 7909, assignments: 8460, inverse_symmetry: (8460, 1), symmetry_factors: &[(1, 8460)], content: 0x963dbcf6ba7b1c88, numerators: 0x62286f983894e5f0, named: 0x3bf16a0388795032 };
    fast qed_eeaaa_L1: Row { diagrams: 384, topologies: 120, assignments: 384, inverse_symmetry: (384, 1), symmetry_factors: &[(1, 384)], content: 0x0fb43c5f515a232f, numerators: 0x044cde5d497d3648, named: 0xf189f8b00263d9bd };
    heavy qed_eeaaa_L2: Row { diagrams: 16344, topologies: 1581, assignments: 16344, inverse_symmetry: (16344, 1), symmetry_factors: &[(1, 16344)], content: 0x98091302c044d5cf, numerators: 0x7a44d325ff935352, named: 0x0f719cac30ea572f };
    fast qed_eeaa_identical_L1: Row { diagrams: 60, topologies: 20, assignments: 66, inverse_symmetry: (60, 1), symmetry_factors: &[(1, 60)], content: 0xc2f14084d3c16139, numerators: 0xd27676fa0c37edb5, named: 0x8e7e16a5e6214aa4 };
    fast qed_eeaa_sign_L1: Row { diagrams: 58, topologies: 20, assignments: 66, inverse_symmetry: (58, 1), symmetry_factors: &[(1, 58)], content: 0xd5ee97232124fe0c, numerators: 0x0193aaeff8c4a870, named: 0xb0aca5b8a7de17e0 };
    fast qed_eeaa_scalar_L1: Row { diagrams: 57, topologies: 20, assignments: 66, inverse_symmetry: (57, 1), symmetry_factors: &[(1, 57)], content: 0x28e8c26fd53fa6bb, numerators: 0xd1995d9a29a96d4f, named: 0xf370e1cb7c88e2c1 };
    fast qed_ee4f_identical_L1: Row { diagrams: 228, topologies: 304, assignments: 242, inverse_symmetry: (228, 1), symmetry_factors: &[(1, 228)], content: 0xac6a16e86141b488, numerators: 0xedd1b0058d7da78c, named: 0x51770c1b6d8c4003 };
    fast qed_ee4f_sign_L1: Row { diagrams: 226, topologies: 304, assignments: 242, inverse_symmetry: (226, 1), symmetry_factors: &[(1, 226)], content: 0x17b789011e6ec656, numerators: 0x83dc420a2d46c1cd, named: 0xae8a66f5b8d4e946 };
    fast qed_ee4f_scalar_L1: Row { diagrams: 225, topologies: 304, assignments: 242, inverse_symmetry: (225, 1), symmetry_factors: &[(1, 225)], content: 0x5ef890dc71487c0a, numerators: 0x599d7c9c1a1dbe29, named: 0x4419b71fdfc0ce0c };
    failing qed_eemm_identical_L2: SCALED_TENSOR;
    failing qed_eemm_sign_L2: SCALED_TENSOR;
    failing qed_eemm_scalar_L2: SCALED_TENSOR;
    failing qed_eeaa_identical_L2: SCALED_TENSOR;
    failing qed_eeaa_sign_L2: SCALED_TENSOR;
    failing qed_eeaa_scalar_L2: SCALED_TENSOR;
    fast qed_ase_L3: Row { diagrams: 540, topologies: 10, assignments: 540, inverse_symmetry: (540, 1), symmetry_factors: &[(1, 540)], content: 0xe97f90cc10f0ba33, numerators: 0x8b2ea87aafb96bb5, named: 0xe9020f2459791f62 };
    heavy qed_ase_L4: Row { diagrams: 12285, topologies: 64, assignments: 12285, inverse_symmetry: (12285, 1), symmetry_factors: &[(1, 12285)], content: 0xdd973742675250fc, numerators: 0xae1fbd18cac0fd7b, named: 0x10e8032b64da8c6a };
    heavy qed_ase_L5: Row { diagrams: 322704, topologies: 519, assignments: 322704, inverse_symmetry: (322704, 1), symmetry_factors: &[(1, 322704)], content: 0x8890f34a21fe1f0b, numerators: 0x49d30cb3c1b10080, named: 0x5a28cadc87770696 };
    fast qed_ese_L3: Row { diagrams: 172, topologies: 10, assignments: 172, inverse_symmetry: (172, 1), symmetry_factors: &[(1, 172)], content: 0x200bb021aa8e3fe3, numerators: 0xa6fb39a8583e54f9, named: 0xefeaddfd39ba99d3 };
    heavy qed_ese_L4: Row { diagrams: 3485, topologies: 64, assignments: 3485, inverse_symmetry: (3485, 1), symmetry_factors: &[(1, 3485)], content: 0x65cd4c6b3be0e21c, numerators: 0x909048b3b2a5b609, named: 0x6b86470c3882d6c0 };
    fast qcd_gggg_L1: Row { diagrams: 314, topologies: 55, assignments: 314, inverse_symmetry: (294, 1), symmetry_factors: &[(1, 274), (2, 40)], content: 0x76e5a6e933609cc9, numerators: 0x18e27150fbb3f7ad, named: 0x10574e7bf8b60a92 };
    heavy qcd_gggg_L2: Row { diagrams: 15253, topologies: 840, assignments: 15253, inverse_symmetry: (163229, 12), symmetry_factors: &[(1, 12086), (2, 2905), (4, 243), (6, 19)], content: 0x858578c792dbbfe3, numerators: 0xd6dbe0288601183a, named: 0x9ebf94ff546a87dd };
    heavy qcd_gggg_L3: Row { diagrams: 686501, topologies: 14340, assignments: 686501, inverse_symmetry: (4715477, 8), symmetry_factors: &[(1, 505145), (2, 157099), (4, 21237), (6, 1512), (8, 1277), (12, 231)], content: 0x68c391d9d10d2457, numerators: 0xabb3ddf64ec13982, named: 0xffc3d4eba05667b3 };
    heavy qcd_ggggg_L1: Row { diagrams: 3075, topologies: 562, assignments: 3075, inverse_symmetry: (2895, 1), symmetry_factors: &[(1, 2715), (2, 360)], content: 0xcde7da42321b8460, numerators: 0x1b818b81bef6124a, named: 0x602f359ddf7e9381 };
    heavy qcd_ggggg_L2: Row { diagrams: 213390, topologies: 11685, assignments: 213390, inverse_symmetry: (765585, 4), symmetry_factors: &[(1, 170990), (2, 39280), (4, 2955), (6, 165)], content: 0x49680cb92696c541, numerators: 0x86c858892c455f99, named: 0x2f0d92e3b159dbe9 };
    fast qcd_uugg_L1: Row { diagrams: 99, topologies: 32, assignments: 99, inverse_symmetry: (94, 1), symmetry_factors: &[(1, 89), (2, 10)], content: 0xcae15c82340cc7a1, numerators: 0xf3021c3c83e886f6, named: 0x72fa4d75fb9e8be7 };
    heavy qcd_uugg_L2: Row { diagrams: 3206, topologies: 457, assignments: 3206, inverse_symmetry: (17479, 6), symmetry_factors: &[(1, 2644), (2, 517), (4, 38), (6, 7)], content: 0xdeb5af9fd1a493e2, numerators: 0xb2cdc545ccb51b8c, named: 0xe64fa32d6a7f6558 };
    heavy qcd_uugg_L3: Row { diagrams: 114301, topologies: 7513, assignments: 114301, inverse_symmetry: (301202, 3), symmetry_factors: &[(1, 88202), (2, 22904), (4, 2686), (6, 317), (8, 152), (12, 40)], content: 0x5a4c42972dd54db7, numerators: 0x137263ecc39a8c6d, named: 0xf1dd4ba547a40558 };
    fast qcd_uuuu_L2: Row { diagrams: 902, topologies: 230, assignments: 902, inverse_symmetry: (5057, 6), symmetry_factors: &[(1, 788), (2, 106), (4, 6), (6, 2)], content: 0x2ad5d9d5a3670a27, numerators: 0xf3a36ce1c109060a, named: 0x5484e8dd34275580 };
    fast qcd_gse_L3: Row { diagrams: 1219, topologies: 57, assignments: 1219, inverse_symmetry: (8553, 8), symmetry_factors: &[(1, 940), (2, 241), (4, 31), (6, 1), (8, 5), (12, 1)], content: 0x7dcb189dcb6f5aa4, numerators: 0x48f27e25d1411024, named: 0x33e84227fa4d6234 };
    heavy qcd_gse_L4: Row { diagrams: 35820, topologies: 666, assignments: 35820, inverse_symmetry: (1485389, 48), symmetry_factors: &[(1, 26801), (2, 7660), (4, 1140), (6, 83), (8, 112), (12, 14), (16, 9), (24, 1)], content: 0x044b7b0e04b1ce58, numerators: 0xf79a9efed6e20b13, named: 0x4192e2410ce67125 };
    heavy qcd_gse_L5: Row { diagrams: 1267351, topologies: 9892, assignments: 1267351, inverse_symmetry: (308996023, 288), symmetry_factors: &[(1, 910237), (2, 298126), (3, 14), (4, 49275), (6, 3114), (8, 5369), (12, 752), (16, 367), (24, 72), (32, 19), (36, 1), (48, 3), (72, 2)], content: 0x5670abb47f3031ca, numerators: 0xb6c8f00b6337680d, named: 0x25cdfa674b0e9786 };
    heavy qcd_gse_all_L3: Row { diagrams: 2627, topologies: 72, assignments: 2627, inverse_symmetry: (4427, 2), symmetry_factors: &[(1, 1871), (2, 626), (4, 102), (6, 15), (8, 10), (12, 3)], content: 0xec89d122df081617, numerators: 0x18c5f56f227468e1, named: 0x5e7eb0c29f9a6d38 };
    fast qcd_tad_L2: Row { diagrams: 125, topologies: 6, assignments: 125, inverse_symmetry: (2383, 24), symmetry_factors: &[(1, 77), (2, 42), (4, 4), (6, 1), (8, 1)], content: 0xe6298357fd676923, numerators: 0xa53bd80db92a26d0, named: 0xecc056ac57d4cda3 };
    heavy qcd_tad_L3: Row { diagrams: 3490, topologies: 46, assignments: 3490, inverse_symmetry: (124199, 48), symmetry_factors: &[(1, 1883), (2, 1247), (3, 7), (4, 271), (6, 29), (8, 43), (12, 5), (16, 4), (48, 1)], content: 0x6b4b08a9e6607c1b, numerators: 0x14089c877c23f095, named: 0x47209285313e59f9 };
    fast qcd_vac_L2: Row { diagrams: 45, topologies: 3, assignments: 45, inverse_symmetry: (191, 6), symmetry_factors: &[(1, 21), (2, 21), (8, 2), (12, 1)], content: 0x0fdd50aef970b13d, numerators: 0x6173058343f36369, named: 0x59d4c7daa0883315 };
    fast qcd_vac_L3: Row { diagrams: 649, topologies: 12, assignments: 649, inverse_symmetry: (2417, 6), symmetry_factors: &[(1, 231), (2, 287), (3, 7), (4, 84), (6, 14), (8, 18), (12, 1), (16, 4), (24, 1), (48, 2)], content: 0x145ae3379080724d, numerators: 0x34b5417dee2213d7, named: 0xe8090ea3bef23517 };
    heavy qcd_vac_L4: Row { diagrams: 14017, topologies: 73, assignments: 14017, inverse_symmetry: (2597681, 288), symmetry_factors: &[(1, 5698), (2, 5383), (3, 154), (4, 1894), (6, 189), (8, 483), (12, 79), (16, 85), (24, 21), (32, 14), (48, 13), (72, 2), (128, 1), (384, 1)], content: 0xc54f17048580de58, numerators: 0xf202887af896ed26, named: 0x46855f3bc1c010db };
    fast qcd_vac1pi_L4: Row { diagrams: 628, topologies: 26, assignments: 628, inverse_symmetry: (11465, 36), symmetry_factors: &[(1, 154), (2, 217), (3, 56), (4, 101), (6, 28), (8, 48), (12, 8), (16, 8), (24, 2), (32, 2), (48, 3), (72, 1)], content: 0x7966c26e2b7b8717, numerators: 0xfad36d65e1a809b0, named: 0xfb6251a1478c5e1e };
    heavy qcd_vac1pi_L5: Row { diagrams: 9889, topologies: 166, assignments: 9889, inverse_symmetry: (7313267, 1152), symmetry_factors: &[(1, 4165), (2, 3455), (4, 1494), (6, 49), (8, 472), (12, 46), (16, 159), (24, 11), (32, 22), (48, 6), (64, 4), (96, 2), (128, 3), (144, 1)], content: 0xf0295cc039fe2a58, numerators: 0x30a89b2dbb1dd465, named: 0x427cfb63687e0040 };
    fast qcd_gggg_sym_L1: Row { diagrams: 123, topologies: 25, assignments: 123, inverse_symmetry: (147, 2), symmetry_factors: &[(1, 37), (2, 62), (4, 20), (8, 4)], content: 0x3ee805c0752694b4, numerators: 0xd9d2fce0e9012de3, named: 0x07d863449a61cb9a };
    failing qcd_gggg_identical_L1: COLOUR_DIMENSION;
    failing qcd_gggg_sign_L1: COLOUR_DIMENSION;
    failing qcd_gggg_scalar_L1: COLOUR_DIMENSION;
    fast qcd_gggg_ghost_filter_L1: Row { diagrams: 277, topologies: 55, assignments: 277, inverse_symmetry: (257, 1), symmetry_factors: &[(1, 237), (2, 40)], content: 0x90a1ce06e025a6f5, numerators: 0x66554bce91f5d5d7, named: 0xbd140c654dd3a59f };
    fast qcd_gggg_ghost_veto_L1: Row { diagrams: 277, topologies: 55, assignments: 277, inverse_symmetry: (257, 1), symmetry_factors: &[(1, 237), (2, 40)], content: 0x90a1ce06e025a6f5, numerators: 0x66554bce91f5d5d7, named: 0xbd140c654dd3a59f };
    fast ym_g6_L0: Row { diagrams: 220, topologies: 220, assignments: 220, inverse_symmetry: (220, 1), symmetry_factors: &[(1, 220)], content: 0xe166545bc76a8cd8, numerators: 0x9f37ff3f3deffb74, named: 0xa09f9df2d844b798 };
    heavy ym_g7_L0: Row { diagrams: 2485, topologies: 2485, assignments: 2485, inverse_symmetry: (2485, 1), symmetry_factors: &[(1, 2485)], content: 0x057203fcdeaf182e, numerators: 0x518deeeb2941c74c, named: 0xf37f5c265e9865c5 };
    heavy ym_g8_L0: Row { diagrams: 34300, topologies: 34300, assignments: 34300, inverse_symmetry: (34300, 1), symmetry_factors: &[(1, 34300)], content: 0x0d8cd8312f854d8b, numerators: 0xe586bec340cace76, named: 0x526fbb59e279711c };
    heavy ym_gse_L4: Row { diagrams: 2571, topologies: 666, assignments: 2571, inverse_symmetry: (81653, 48), symmetry_factors: &[(1, 1139), (2, 889), (4, 402), (6, 17), (8, 100), (12, 14), (16, 9), (24, 1)], content: 0xbb8866a313b845e9, numerators: 0x3ec48d93c61cce4c, named: 0xc3d58015489276ab };
    heavy ym_gse_L5: Row { diagrams: 53201, topologies: 9892, assignments: 53201, inverse_symmetry: (9957739, 288), symmetry_factors: &[(1, 22781), (2, 18577), (3, 2), (4, 8421), (6, 315), (8, 2411), (12, 248), (16, 349), (24, 72), (32, 19), (36, 1), (48, 3), (72, 2)], content: 0x64121e08ff2ba484, numerators: 0x204eb02540fdbc12, named: 0xf3f645c4a0b653f1 };
    fast ym_vac_L5: Row { diagrams: 381, topologies: 91, assignments: 381, inverse_symmetry: (95437, 576), symmetry_factors: &[(1, 82), (2, 110), (4, 73), (6, 4), (8, 54), (12, 9), (16, 28), (24, 3), (32, 8), (48, 4), (64, 2), (96, 1), (128, 2), (144, 1)], content: 0x8338e2f92f7cfd9d, numerators: 0xd4db3703b3a0e21b, named: 0x787b52f109f2435f };
    fast sm_eemm_L0: Row { diagrams: 4, topologies: 1, assignments: 4, inverse_symmetry: (4, 1), symmetry_factors: &[(1, 4)], content: 0x97af00ae91ef62c4, numerators: 0x7c4bccf1c63f9180, named: 0x28df99f11af66fba };
    fast sm_eemm_L1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0x99843c92c72133b6, numerators: 0x26f3641db6659cf2, named: 0x6c75ac1959f5ae62 };
    heavy sm_eemm_L2: Row { diagrams: 69697, topologies: 160, assignments: 69697, inverse_symmetry: (415391, 6), symmetry_factors: &[(1, 68781), (2, 888), (4, 26), (6, 2)], content: 0xb4b43e66a248f228, numerators: 0x5a12870245dd2b8c, named: 0x91591ac10e9d39b2 };
    failing sm_eemm_identical_L1: NON_SELF_DUAL_POWER;
    failing sm_eemm_sign_L1: NON_SELF_DUAL_POWER;
    failing sm_eemm_scalar_L1: NON_SELF_DUAL_POWER;
    fast sm_ee4f_L0: Row { diagrams: 18, topologies: 7, assignments: 18, inverse_symmetry: (18, 1), symmetry_factors: &[(1, 18)], content: 0x7dab707954851b41, numerators: 0x71eb659079856e7a, named: 0x7f394436ced3aaa8 };
    heavy sm_ee4f_L1: Row { diagrams: 5811, topologies: 322, assignments: 5811, inverse_symmetry: (5805, 1), symmetry_factors: &[(1, 5799), (2, 12)], content: 0xaea64b8a548b0cf6, numerators: 0xe557f43bb1948d11, named: 0x16b043c30b0e21e3 };
    fast sm_eeww_L0: Row { diagrams: 4, topologies: 2, assignments: 4, inverse_symmetry: (4, 1), symmetry_factors: &[(1, 4)], content: 0x7477c6f9f9a39d73, numerators: 0x63db443730cf91ef, named: 0x48249b27b144364c };
    fast sm_eeww_L1: Row { diagrams: 919, topologies: 23, assignments: 919, inverse_symmetry: (916, 1), symmetry_factors: &[(1, 913), (2, 6)], content: 0x953e54839edd5025, numerators: 0x077e524238dbfee5, named: 0x1f717b4d5120a8b6 };
    fast sm_ggtt_L1: Row { diagrams: 375, topologies: 32, assignments: 375, inverse_symmetry: (370, 1), symmetry_factors: &[(1, 365), (2, 10)], content: 0xb57aa9d432880775, numerators: 0x6c6e9542edb615e9, named: 0xdf6790aae2a0199a };
    fast sm_ggttg_L0: Row { diagrams: 16, topologies: 16, assignments: 16, inverse_symmetry: (16, 1), symmetry_factors: &[(1, 16)], content: 0xeb9e3621b43c55e4, numerators: 0x59b38f69f6396929, named: 0x1e50e6944663b148 };
    heavy sm_ggttg_L1: Row { diagrams: 3096, topologies: 322, assignments: 3096, inverse_symmetry: (3046, 1), symmetry_factors: &[(1, 2996), (2, 100)], content: 0x60920f521092c4ea, numerators: 0x96c505e3b065e987, named: 0x17907d2e7d841a52 };
    fast sm_hse_L1: Row { diagrams: 17, topologies: 1, assignments: 17, inverse_symmetry: (31, 2), symmetry_factors: &[(1, 14), (2, 3)], content: 0x825faa38b23aed38, numerators: 0xdf018e8ddee914be, named: 0x50589e485c2a732b };
    fast sm_hse_L2: Row { diagrams: 1188, topologies: 7, assignments: 1188, inverse_symmetry: (13691, 12), symmetry_factors: &[(1, 1100), (2, 76), (4, 11), (6, 1)], content: 0x1c50eecb5249aa7d, numerators: 0xe7757c1488aeba4b, named: 0xf7cfa4073867593b };
    fast sm_ase_L2: Row { diagrams: 904, topologies: 7, assignments: 904, inverse_symmetry: (904, 1), symmetry_factors: &[(1, 904)], content: 0xe45d4a5e047a542c, numerators: 0x5b5f703fd30dfaa8, named: 0x94d4e5364dd94100 };
    fast sm_eemm_qcd0_L1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0x99843c92c72133b6, numerators: 0x26f3641db6659cf2, named: 0x6c75ac1959f5ae62 };
    heavy sm_eemm_qcd0_L2: Row { diagrams: 69193, topologies: 160, assignments: 69193, inverse_symmetry: (412367, 6), symmetry_factors: &[(1, 68277), (2, 888), (4, 26), (6, 2)], content: 0xe342a202fff52526, numerators: 0x3a0eab7af9b8b1a9, named: 0x6541eaf01bd37b2c };
    fast sm_uudd_qed2_L1: Row { diagrams: 111, topologies: 17, assignments: 111, inverse_symmetry: (221, 2), symmetry_factors: &[(1, 110), (2, 1)], content: 0x422125ee75cb81f1, numerators: 0x0e06bd22fdc90a2f, named: 0xd4111bd0134ba6da };
    heavy sm_uudd_qed2_L2: Row { diagrams: 3286, topologies: 230, assignments: 3286, inverse_symmetry: (38519, 12), symmetry_factors: &[(1, 3136), (2, 146), (4, 3), (6, 1)], content: 0xa681cec8fde8b382, numerators: 0x815c7654c050793d, named: 0x77d595df84089ea6 };
    fast ew_eemm_L1: Row { diagrams: 428, topologies: 10, assignments: 428, inverse_symmetry: (853, 2), symmetry_factors: &[(1, 425), (2, 3)], content: 0xb56d549c918bbbc1, numerators: 0x2bc439d340929b34, named: 0x927d3722af5f0c1f };
    heavy ew_eemm_L2: Row { diagrams: 69193, topologies: 160, assignments: 69193, inverse_symmetry: (412367, 6), symmetry_factors: &[(1, 68277), (2, 888), (4, 26), (6, 2)], content: 0x270481b3cce8b479, numerators: 0x6292f42e7e2fb42d, named: 0x2e9a8998ddbd7da4 };
    fast ct_phi3_2pt_L1: Row { diagrams: 4, topologies: 1, assignments: 4, inverse_symmetry: (2, 1), symmetry_factors: &[(2, 4)], content: 0x792c2689f9b727df, numerators: 0x094dd49c998c2915, named: 0x59fed4ad1262d0a7 };
    fast ct_qcd_gggg_L1: Row { diagrams: 1607, topologies: 55, assignments: 1607, inverse_symmetry: (2883, 2), symmetry_factors: &[(1, 1276), (2, 331)], content: 0x102a6a28e1c425f5, numerators: 0x66ea9c29e5d8a21d, named: 0xeac5795af3c09462 };
    fast ct_qcd_uugg_L1: Row { diagrams: 864, topologies: 32, assignments: 864, inverse_symmetry: (802, 1), symmetry_factors: &[(1, 740), (2, 124)], content: 0xbdb6f0a3f62831ef, numerators: 0xab75b6346178447f, named: 0x63ba7eed6e99be72 };
    fast ct_phi3_2pt_tad_L2: Row { diagrams: 148, topologies: 10, assignments: 148, inverse_symmetry: (50, 1), symmetry_factors: &[(1, 4), (2, 44), (4, 92), (8, 8)], content: 0xf5122b0d42fad77e, numerators: 0xb1c2ba48ec6c4558, named: 0x2d384ebf47c64577 };
    fast ct_phi3_vac_L2: Row { diagrams: 6, topologies: 2, assignments: 6, inverse_symmetry: (5, 6), symmetry_factors: &[(4, 1), (6, 1), (8, 2), (12, 2)], content: 0x9dad3c92417664f1, numerators: 0xf57c18fe6786fb79, named: 0x0814e18440d40ff7 };
    fast ct_phi3_vac_L3: Row { diagrams: 42, topologies: 5, assignments: 42, inverse_symmetry: (5, 1), symmetry_factors: &[(4, 7), (6, 2), (8, 17), (16, 10), (24, 2), (48, 4)], content: 0xb5863cd09c348f6c, numerators: 0xcadc67ebf2b44f13, named: 0x45e1c5d3b0e14d48 };
    fast ct_phi3_vac_L4: Row { diagrams: 505, topologies: 17, assignments: 505, inverse_symmetry: (1105, 18), symmetry_factors: &[(1, 1), (2, 17), (4, 76), (6, 1), (8, 161), (12, 6), (16, 153), (24, 1), (32, 70), (36, 1), (48, 6), (64, 6), (72, 2), (128, 4)], content: 0xdac00a3ffc146b39, numerators: 0x758195af180c3984, named: 0xb5fdaf4e73c26499 };
    fast ct_qcd_vac_L2: Row { diagrams: 68, topologies: 3, assignments: 68, inverse_symmetry: (1103, 24), symmetry_factors: &[(1, 29), (2, 32), (4, 1), (6, 1), (8, 3), (12, 2)], content: 0xf69da97cd4bfb079, numerators: 0x74cb264b3044f795, named: 0xaff30af9d3a7ec43 };
    fast ct_qcd_vac_L3: Row { diagrams: 1816, topologies: 12, assignments: 1816, inverse_symmetry: (6341, 6), symmetry_factors: &[(1, 549), (2, 826), (3, 16), (4, 293), (6, 34), (8, 74), (12, 4), (16, 13), (24, 2), (48, 5)], content: 0x78eccaa3c84acd2a, numerators: 0xca96a9c3deeb724f, named: 0x05091949be52c84a };
    fast xs_qed_eemm_L1: Row { diagrams: 1, topologies: 17, assignments: 34, inverse_symmetry: (1, 1), symmetry_factors: &[(1, 1)], content: 0x831ae534b022fd0f, numerators: 0xe134f67914689408, named: 0xf40bfdc791015105 };
    fast xs_qed_eemm_L2: Row { diagrams: 28, topologies: 157, assignments: 714, inverse_symmetry: (28, 1), symmetry_factors: &[(1, 28)], content: 0xf6aad76a05e6dd4f, numerators: 0x574b5a82fb6414f7, named: 0xe3328e6593163eb3 };
    failing xs_qed_eemm_L3: THRESHOLD_CANDIDATE;
    fast xs_qed_eemm_identical_L2: Row { diagrams: 26, topologies: 157, assignments: 714, inverse_symmetry: (26, 1), symmetry_factors: &[(1, 26)], content: 0xc7c61c59da4c4d71, numerators: 0x30199ed95f9efead, named: 0x7c26a4b34e3f02cd };
    fast xs_qed_eemm_sign_L2: Row { diagrams: 26, topologies: 157, assignments: 714, inverse_symmetry: (26, 1), symmetry_factors: &[(1, 26)], content: 0xc7c61c59da4c4d71, numerators: 0x30199ed95f9efead, named: 0x7c26a4b34e3f02cd };
    fast xs_qed_eemm_scalar_L2: Row { diagrams: 26, topologies: 157, assignments: 714, inverse_symmetry: (26, 1), symmetry_factors: &[(1, 26)], content: 0xc7c61c59da4c4d71, numerators: 0x30199ed95f9efead, named: 0x7c26a4b34e3f02cd };
    fast xs_qed_eedd_L1: Row { diagrams: 1, topologies: 17, assignments: 34, inverse_symmetry: (1, 1), symmetry_factors: &[(1, 1)], content: 0xb467d6be7a12aefa, numerators: 0x52ae30868db3c38d, named: 0xdd30b2e206f30267 };
    fast xs_qed_eedd_L2: Row { diagrams: 28, topologies: 157, assignments: 714, inverse_symmetry: (28, 1), symmetry_factors: &[(1, 28)], content: 0x6ad9c202edc60059, numerators: 0x9c41ac8d6a4ba8a0, named: 0x33c78f42a6595910 };
    fast xs_qed_eeaa_L2: Row { diagrams: 94, topologies: 157, assignments: 714, inverse_symmetry: (94, 1), symmetry_factors: &[(1, 94)], content: 0x632909095d92e259, numerators: 0xecb88282186d029c, named: 0xf9abb56eaadf0939 };
    fast xs_qed_eeee_L1: Row { diagrams: 4, topologies: 17, assignments: 34, inverse_symmetry: (4, 1), symmetry_factors: &[(1, 4)], content: 0x29f20710d0fbfbd6, numerators: 0x099879207707ff29, named: 0x7044ad165b20a02d };
    fast xs_qed_eeee_L2: Row { diagrams: 98, topologies: 157, assignments: 714, inverse_symmetry: (98, 1), symmetry_factors: &[(1, 98)], content: 0x84d794d7e60d51ac, numerators: 0xf23e55d60340c98f, named: 0xf986e520f619ff70 };
    fast xs_qed_eeee_1pi_L1: Row { diagrams: 1, topologies: 3, assignments: 4, inverse_symmetry: (1, 1), symmetry_factors: &[(1, 1)], content: 0x8686eb17d2b0f293, numerators: 0xd82920bd552e843e, named: 0xf5ac9fd01b9eabd7 };
    fast xs_qed_eeee_1pi_L2: Row { diagrams: 22, topologies: 36, assignments: 116, inverse_symmetry: (22, 1), symmetry_factors: &[(1, 22)], content: 0x66424d780b8aa80f, numerators: 0x92ea3d0a49dec7f4, named: 0x0ed986d818d723c1 };
    fast xs_sm_eemm_L1: Row { diagrams: 16, topologies: 17, assignments: 856, inverse_symmetry: (16, 1), symmetry_factors: &[(1, 16)], content: 0x9ad7451811ae86e7, numerators: 0x01b59b2f2dfb564d, named: 0x6aa241cf0ee01c70 };
    fast xs_sm_eemm_L2: Row { diagrams: 2912, topologies: 230, assignments: 139394, inverse_symmetry: (2900, 1), symmetry_factors: &[(1, 2888), (2, 24)], content: 0xa0b2d1d75cc290ea, numerators: 0xbd7fc115267a6235, named: 0x043c2dcd07241967 };
    fast xs_sm_uudd_L1: Row { diagrams: 16, topologies: 17, assignments: 558, inverse_symmetry: (16, 1), symmetry_factors: &[(1, 16)], content: 0x6a0f6683f376c99d, numerators: 0x7d74ed995684f1db, named: 0xb627b9d5aeb8cfe9 };
    fast xs_qcd_gggg_L1: Row { diagrams: 10, topologies: 55, assignments: 314, inverse_symmetry: (8, 1), symmetry_factors: &[(1, 6), (2, 4)], content: 0xf64ad313fb15b61e, numerators: 0xbde2aaa4f3b826c0, named: 0x17da1ea742d3ec86 };
    fast xs_qcd_gggg_sym_L1: Row { diagrams: 10, topologies: 55, assignments: 314, inverse_symmetry: (4, 1), symmetry_factors: &[(2, 6), (4, 4)], content: 0xb7b93daa191de036, numerators: 0xbde2aaa4f3b826c0, named: 0xbb1684fcd423ed44 };
    failing xs_qcd_gggg_L2: THRESHOLD_CANDIDATE;
    fast xs_qcd_gguu_L1: Row { diagrams: 9, topologies: 55, assignments: 314, inverse_symmetry: (9, 1), symmetry_factors: &[(1, 9)], content: 0xf55ef4afe4fbeee8, numerators: 0x9e713702e766065d, named: 0x967e2026a5942757 };
    failing xs_qcd_gguu_L2: THRESHOLD_CANDIDATE;
    fast xs_qcd_gguu_sym_L1: Row { diagrams: 5, topologies: 55, assignments: 314, inverse_symmetry: (9, 2), symmetry_factors: &[(1, 4), (2, 1)], content: 0xaaea5535e56a1d82, numerators: 0xa9e3dfa7c6a49ca8, named: 0x326170173160f176 };
    fast xs_qcd_gguu_1pi_L1: Row { diagrams: 4, topologies: 12, assignments: 54, inverse_symmetry: (4, 1), symmetry_factors: &[(1, 4)], content: 0x4ee6aab30b16f1c5, numerators: 0x901412cc9fc683f1, named: 0x2dec9685df5b0e1c };
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
        .map(FeynmanDiagram::id)
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
