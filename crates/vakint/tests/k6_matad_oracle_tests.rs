//! Offline MATAD oracle for RustRed's current three-loop K=6 diagnostics.
//!
//! These tests are deliberately ignored by default: they invoke FORM/MATAD and
//! are a development-time fixture/review surface, never part of Vakint's
//! production RustRed scalar path. Run them serially with `--ignored --nocapture`
//! to obtain the exact canonical records printed below.

use std::{collections::HashMap, env};

use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::SingleFloat,
};
use vakint::{
    EvaluationOrder, LoopNormalizationFactor, MATADOptions, NumericalEvaluationResult, Vakint,
    VakintSettings, vakint_parse,
};

#[derive(Clone, Copy, Debug)]
struct K6OracleCase {
    name: &'static str,
    powers: [i64; 6],
    expected_master_markers: &'static [&'static str],
    expected_leading_power: i64,
    numerical_probe: bool,
    canonical_len: usize,
    canonical_fnv1a64: u64,
}

impl K6OracleCase {
    fn input(self) -> Atom {
        let [n1, n2, n3, n4, n5, n6] = self.powers;
        vakint_parse!(format!("topo(I3L(muvsq,{n1},{n2},{n3},{n4},{n5},{n6}))").as_str()).unwrap()
    }
}

// The first three rows are the scalar corners. The remaining rows are sample
// points on the six positive-dimensional recurrence strata recorded in
// RustRed's docs/research/vakint_k6_oracle.md. They are diagnostics, not a
// claim that six sampled reductions prove all-rank closure.
const K6_ORACLE_CASES: [K6OracleCase; 9] = [
    K6OracleCase {
        name: "corner_0_1_1_1_1_0",
        powers: [0, 1, 1, 1, 1, 0],
        expected_master_markers: &["miBN"],
        expected_leading_power: -3,
        numerical_probe: true,
        canonical_len: 80,
        canonical_fnv1a64: 0x366a693f6617cb9a,
    },
    K6OracleCase {
        name: "corner_0_1_1_1_1_1",
        powers: [0, 1, 1, 1, 1, 1],
        expected_master_markers: &["Gam(1,1)", "miD5", "miT111", "miBN"],
        expected_leading_power: -3,
        numerical_probe: false,
        canonical_len: 591,
        canonical_fnv1a64: 0xeb9d4be075a9a565,
    },
    K6OracleCase {
        name: "corner_1_1_1_1_1_1",
        powers: [1, 1, 1, 1, 1, 1],
        expected_master_markers: &["miD6"],
        expected_leading_power: -1,
        numerical_probe: true,
        canonical_len: 60,
        canonical_fnv1a64: 0x2557a97baaa48ca7,
    },
    K6OracleCase {
        name: "witness_0_m1_1_2_2_1",
        powers: [0, -1, 1, 2, 2, 1],
        expected_master_markers: &["Gam(1,1)", "miT111"],
        expected_leading_power: -3,
        numerical_probe: true,
        canonical_len: 337,
        canonical_fnv1a64: 0xfe6ce8982652fc41,
    },
    K6OracleCase {
        name: "witness_0_m2_2_2_1_1",
        powers: [0, -2, 2, 2, 1, 1],
        expected_master_markers: &["Gam(1,1)", "miT111"],
        expected_leading_power: -3,
        numerical_probe: false,
        canonical_len: 447,
        canonical_fnv1a64: 0x8b7ad7fbe3afa130,
    },
    K6OracleCase {
        name: "witness_0_1_1_2_4_0",
        powers: [0, 1, 1, 2, 4, 0],
        expected_master_markers: &["Gam(1,1)", "miBN"],
        expected_leading_power: -2,
        numerical_probe: false,
        canonical_len: 646,
        canonical_fnv1a64: 0xd2ee77c7534dbb94,
    },
    K6OracleCase {
        name: "witness_0_1_1_2_5_0",
        powers: [0, 1, 1, 2, 5, 0],
        expected_master_markers: &["Gam(1,1)", "miBN"],
        expected_leading_power: -2,
        numerical_probe: true,
        canonical_len: 739,
        canonical_fnv1a64: 0xc89213fa29375533,
    },
    K6OracleCase {
        name: "witness_0_1_2_3_3_0",
        powers: [0, 1, 2, 3, 3, 0],
        expected_master_markers: &["Gam(1,1)", "miBN"],
        expected_leading_power: 0,
        numerical_probe: false,
        canonical_len: 597,
        canonical_fnv1a64: 0x8c761fe7ae99ce23,
    },
    K6OracleCase {
        name: "witness_0_1_3_2_3_0",
        powers: [0, 1, 3, 2, 3, 0],
        expected_master_markers: &["Gam(1,1)", "miBN"],
        expected_leading_power: 0,
        numerical_probe: false,
        canonical_len: 597,
        canonical_fnv1a64: 0x8c761fe7ae99ce23,
    },
];

const KNOWN_RAW_MASTER_ATOMS: [&str; 12] = [
    "Gam(1,1)",
    "iGam(1,1)",
    "miD6",
    "miD5",
    "miD4",
    "miD3",
    "miDM",
    "miDN",
    "miE3",
    "miBN",
    "miBN1",
    "miT111",
];

fn oracle_form_path() -> String {
    env::var("VAKINT_K6_ORACLE_FORM_PATH")
        .or_else(|_| env::var("FORM_PATH"))
        .expect("offline K6 oracle requires VAKINT_K6_ORACLE_FORM_PATH or FORM_PATH")
}

fn matad_settings(expand_masters: bool, precision: u32) -> VakintSettings {
    VakintSettings {
        form_exe_path: oracle_form_path(),
        evaluation_order: EvaluationOrder::matad_only(Some(MATADOptions {
            expand_masters,
            ..MATADOptions::default()
        })),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        run_time_decimal_precision: precision,
        number_of_terms_in_epsilon_expansion: 5,
        use_dot_product_notation: true,
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    }
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf29ce484222325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x100000001b3)
    })
}

fn contains_exact_atom(expression: &Atom, needle: &str) -> bool {
    let needle = vakint_parse!(needle).unwrap();
    expression
        .pattern_match(&needle.to_pattern(), None, None)
        .next()
        .is_some()
}

fn representative_laurent_fixture(
    case_name: &str,
    settings: &VakintSettings,
) -> NumericalEvaluationResult {
    let coefficients = match case_name {
        "corner_0_1_1_1_1_0" => vec![
            (-3, "2"),
            (
                -2,
                "7.666666666666666666666666666666666666666666666666666666666667",
            ),
            (
                -1,
                "22.43480220054467930941724549993807556765684970362039531320667",
            ),
            (
                0,
                "39.42929462910208211529996476007305636115461794586385093704104",
            ),
            (
                1,
                "62.92709375535910070547798992048699891696212877853418836763700",
            ),
        ],
        "corner_1_1_1_1_1_1" => vec![
            (
                -1,
                "2.404113806319188570799476323022899981529972584680997763584543",
            ),
            (
                0,
                "-10.03527847976878917191470068515890023865033349600275134054453",
            ),
            (
                1,
                "41.87670208303157617433490267097099146643159391700724902140048",
            ),
        ],
        "witness_0_m1_1_2_2_1" => vec![
            (-3, "1"),
            (-2, "1"),
            (
                -1,
                "2.686098687375853357841435320344945427463290515264912236381237",
            ),
            (
                0,
                "2.829900494804866407486639032411391520256197193198133719959162",
            ),
            (
                1,
                "2.077036576556416117978234727093020549151088290005534124488403",
            ),
        ],
        "witness_0_1_1_2_5_0" => vec![
            (
                -2,
                "-0.04166666666666666666666666666666666666666666666666666666666667",
            ),
            (
                -1,
                "-0.1041666666666666666666666666666666666666666666666666666666667",
            ),
            (
                0,
                "-0.1171667387282668606472629198333950649195973264340161704957760",
            ),
            (
                1,
                "-0.2716215592002912016676959183710054639983251635341601068823199",
            ),
        ],
        _ => panic!("no representative Laurent fixture for {case_name}"),
    };

    NumericalEvaluationResult::from_vec(
        coefficients
            .into_iter()
            .map(|(power, real)| (power, (real.to_owned(), "0".to_owned())))
            .collect(),
        settings,
    )
}

#[test]
#[ignore = "offline oracle: invokes FORM 4.3 and MATAD for all nine K=6 diagnostics"]
fn matad_k6_exact_raw_master_records() {
    let vakint = Vakint::new().unwrap();
    let settings = matad_settings(false, 80);
    let mut records = HashMap::new();

    for case in K6_ORACLE_CASES {
        let input = case.input();
        let result = vakint
            .evaluate_integral(&settings, input.as_view())
            .unwrap_or_else(|error| panic!("MATAD failed for {}: {error}", case.name));
        let canonical = result.to_canonical_string();
        let fingerprint = fnv1a64(canonical.as_bytes());

        assert!(
            !result.is_zero(),
            "{} unexpectedly reduced to zero",
            case.name
        );
        assert!(
            !canonical.contains("topo(") && !canonical.contains("I3L("),
            "{} retained an unreduced integral: {canonical}",
            case.name
        );
        for marker in KNOWN_RAW_MASTER_ATOMS {
            let expected = case.expected_master_markers.contains(&marker);
            assert_eq!(
                contains_exact_atom(&result, marker),
                expected,
                "{} exact raw-master support differs at {marker}: {canonical}",
                case.name,
            );
        }

        println!(
            "K6_MATAD_RAW_BEGIN name={} powers={:?} bytes={} fnv1a64={:#018x}\n{}\nK6_MATAD_RAW_END name={}",
            case.name,
            case.powers,
            canonical.len(),
            fingerprint,
            canonical,
            case.name
        );

        assert_eq!(canonical.len(), case.canonical_len, "{} length", case.name);
        assert_eq!(
            fingerprint, case.canonical_fnv1a64,
            "{} canonical fingerprint",
            case.name
        );
        records.insert(case.name, canonical);
    }

    assert_eq!(
        records["witness_0_1_2_3_3_0"], records["witness_0_1_3_2_3_0"],
        "the final witness pair is an exact graph-symmetry orbit at this point"
    );
}

#[test]
#[ignore = "offline oracle: invokes FORM 4.3/MATAD and expands representative K=6 Laurent series"]
fn matad_k6_representative_laurent_records() {
    const PRECISION: u32 = 80;

    let vakint = Vakint::new().unwrap();
    let settings = matad_settings(true, PRECISION);
    let parameter_values = HashMap::from([("muvsq".to_owned(), 1.0), ("mursq".to_owned(), 1.0)]);
    let parameters = vakint.params_from_f64(&settings, &parameter_values);
    let probe_all = env::var_os("VAKINT_K6_ORACLE_NUMERICAL_ALL").is_some();

    for case in K6_ORACLE_CASES
        .iter()
        .copied()
        .filter(|case| probe_all || case.numerical_probe)
    {
        let input = case.input();
        let result = vakint
            .evaluate_integral(&settings, input.as_view())
            .unwrap_or_else(|error| panic!("MATAD failed for {}: {error}", case.name));
        let (numerical, error) = Vakint::full_numerical_evaluation(
            &settings,
            result.as_view(),
            &parameters,
            &HashMap::default(),
            None,
        )
        .unwrap_or_else(|error| panic!("numerical evaluation failed for {}: {error}", case.name));
        assert!(
            error.is_none(),
            "{} unexpectedly carried an error series",
            case.name
        );

        let coefficients = numerical.get_epsilon_coefficients();
        assert!(
            !coefficients.is_empty(),
            "{} returned no Laurent terms",
            case.name
        );
        assert!(
            coefficients
                .iter()
                .all(|(power, coefficient)| (-3..=1).contains(power) && coefficient.im.is_zero()),
            "{} returned a coefficient outside epsilon^-3..epsilon^1 or with an imaginary part: {numerical}",
            case.name
        );
        assert_eq!(
            coefficients.first().map(|(power, _)| *power),
            Some(case.expected_leading_power),
            "{} leading pole depth changed",
            case.name
        );
        assert_eq!(
            coefficients.last().map(|(power, _)| *power),
            Some(1),
            "{} did not reach MATAD's three-loop epsilon^1 ceiling",
            case.name
        );

        if case.numerical_probe {
            let expected = representative_laurent_fixture(case.name, &settings);
            let (matches, detail) = numerical.does_approx_match(&expected, None, 1.0e-55, 0.0);
            assert!(matches, "{} Laurent fixture mismatch: {detail}", case.name);
        }

        println!(
            "K6_MATAD_LAURENT_BEGIN name={} powers={:?} decimal_precision={PRECISION}\n{}\nK6_MATAD_LAURENT_END name={}",
            case.name, case.powers, numerical, case.name
        );
    }
}
