//! vakint against the reference values of the original MATAD/FMFT comparisons.
//!
//! One test per physical integral. Each test writes the integral in vakint's
//! `topo(prop(...))` notation, evaluates it with vakint's default engine (RustRed
//! with its Ordered scheduler, masters substituted numerically, FORM disabled)
//! and compares every requested Laurent coefficient with frozen reference
//! values taken from the original vakint tests: analytic MATAD/alphaLoop
//! results (one to three loops), FMFT (four loops) and, where an independent
//! numerical run exists, pySecDec. A coefficient vakint or a reference does not
//! list is an exact zero (vanishing orders are written as `"0"`), so spurious
//! poles fail. The printed table shows vakint's value, each reference, the
//! relative difference and the verdict; run with `--nocapture` to see it:
//!
//! ```text
//! cargo test --release -p vakint --test reference_tables -- --nocapture
//! ```
//!
//! Four-loop cases load the official external data (`VAKINT_TEST_EXTERNAL_DATA`
//! selects a local checkout). Exact (non-numerical) checks of the same engine
//! live in `ordered_reduction_rungs.rs` and the RustRed parity tests.

#[path = "support/reference_table.rs"]
#[allow(dead_code)]
mod reference_table;
mod test_utils;

use reference_table::{LaurentTable, Reference};
use std::collections::HashMap;
use std::time::Instant;
use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

/// Decimal digits of every evaluation (the references carry 31-32 digits).
const DIGITS: u32 = 32;
/// Accepted relative difference to 32-digit analytic references (MATAD,
/// alphaLoop, FMFT): agreement to at least 28 digits.
const ANALYTIC: f64 = 1e-28;
/// Accepted relative difference to FMFT's four-loop values. FMFT's built-in
/// master constants carry about 27 significant digits, which limits its finite
/// parts (its poles agree to all 32 digits); see the independent 40-digit
/// reference of the H topology below.
const FMFT: f64 = 1e-25;
/// Accepted relative difference to pySecDec's Monte-Carlo values, the
/// threshold of the original analytic-versus-pySecDec comparisons.
const PYSECDEC: f64 = 1e-7;
/// RustRed never runs FORM; a FORM call would fail on this path.
const NO_FORM: &str = "/vakint-reference-tables-never-run-form";

/// External momenta `p(i, μ)` of the original tests, `i = 1..=n`.
#[derive(Clone, Copy)]
enum Momenta {
    None,
    /// `p_i = (0.17(i+1), 0.4(i+2), 0.3(i+3), 0.12(i+4))`.
    Small(usize),
    /// `p_i = (17(i+1), 4(i+2), 3(i+3), 12(i+4))`.
    Large(usize),
}

impl Momenta {
    fn components(self) -> Vec<(usize, (f64, f64, f64, f64))> {
        let (count, scale) = match self {
            Momenta::None => (0, [0.0; 4]),
            Momenta::Small(count) => (count, [0.17, 0.4, 0.3, 0.12]),
            Momenta::Large(count) => (count, [17.0, 4.0, 3.0, 12.0]),
        };
        (1..=count)
            .map(|i| {
                let shift = |offset: usize, factor: f64| factor * (i + offset) as f64;
                (
                    i,
                    (
                        shift(1, scale[0]),
                        shift(2, scale[1]),
                        shift(3, scale[2]),
                        shift(4, scale[3]),
                    ),
                )
            })
            .collect()
    }

    fn describe(self) -> String {
        match self {
            Momenta::None => "no external momenta".into(),
            Momenta::Small(count) => {
                format!("p_i = (0.17(i+1), 0.4(i+2), 0.3(i+3), 0.12(i+4)), i = 1..{count}")
            }
            Momenta::Large(count) => {
                format!("p_i = (17(i+1), 4(i+2), 3(i+3), 12(i+4)), i = 1..{count}")
            }
        }
    }
}

/// One vacuum integral, the point it is evaluated at and its references.
struct Integral {
    name: &'static str,
    loops: usize,
    input: &'static str,
    scheme: LoopNormalizationFactor,
    /// Number of Laurent coefficients, starting at ε^-loops.
    laurent_terms: i64,
    point: &'static [(&'static str, f64)],
    momenta: Momenta,
    /// Reference values of every requested order; an order a reference does
    /// not list must vanish.
    references: Vec<Reference>,
}

fn check(integral: Integral) {
    test_utils::run_multi_lane_acceptance(move || {
        if integral.loops >= 4 {
            test_utils::ensure_external_data_loaded();
        }
        let vakint = test_utils::get_vakint(VakintSettings {
            evaluation_order: EvaluationOrder::rustred_only(),
            integral_normalization_factor: integral.scheme.clone(),
            number_of_terms_in_epsilon_expansion: integral.laurent_terms,
            run_time_decimal_precision: DIGITS,
            use_dot_product_notation: true,
            form_exe_path: NO_FORM.into(),
            ..VakintSettings::default()
        });
        let input = vakint::vakint_parse!(integral.input).expect("integral input parses");
        let started = Instant::now();
        let evaluated = vakint
            .evaluate(input.as_view())
            .unwrap_or_else(|error| panic!("{}: vakint evaluation failed: {error}", integral.name));
        let seconds = started.elapsed().as_secs_f64();
        let point: HashMap<String, f64, ahash::RandomState> = integral
            .point
            .iter()
            .map(|(name, value)| ((*name).to_owned(), *value))
            .collect();
        let momenta: HashMap<usize, (f64, f64, f64, f64), ahash::RandomState> =
            integral.momenta.components().into_iter().collect();
        let (numerical, _) = Vakint::full_numerical_evaluation(
            &vakint.settings,
            evaluated.as_view(),
            &test_utils::params_at_reference_point(&point, DIGITS),
            &HashMap::default(),
            Some(&test_utils::externals_at_reference_point(&momenta, DIGITS)),
        )
        .unwrap_or_else(|error| panic!("{}: numerical evaluation failed: {error}", integral.name));
        let first_order = -(integral.loops as i64);
        let point = integral
            .point
            .iter()
            .map(|(name, value)| format!("{name} = {value}"))
            .collect::<Vec<_>>()
            .join(", ");
        let table = LaurentTable {
            title: integral.name.to_owned(),
            facts: vec![
                (
                    "input".into(),
                    integral.input.split_whitespace().collect::<String>(),
                ),
                (
                    "point".into(),
                    format!("{point}; {}", integral.momenta.describe()),
                ),
                (
                    "scheme".into(),
                    format!(
                        "{:?} normalization, {} Laurent terms from ε^-{}, {DIGITS} digits",
                        integral.scheme, integral.laurent_terms, integral.loops
                    ),
                ),
                (
                    "vakint".into(),
                    format!(
                        "RustRed, Ordered scheduler ({} weights, {} cores), numerical masters, \
                         no FORM; reduction {seconds:.2} s",
                        vakint.settings.ordered_weight_arithmetic, vakint.settings.n_cores
                    ),
                ),
            ],
            engine: "vakint".into(),
            values: numerical.get_epsilon_coefficients(),
            orders: Some(first_order..=first_order + integral.laurent_terms - 1),
            shown_digits: DIGITS as usize,
        };
        if let Err(failure) = table.check(&integral.references) {
            panic!("{failure}");
        }
    });
}

#[test]
fn one_loop_tadpole() {
    check(Integral {
        name: "1L tadpole",
        loops: 1,
        input: "(1)*topo(prop(1,edge(1,1),k(1),muvsq,1))",
        scheme: LoopNormalizationFactor::pySecDec,
        laurent_terms: 2,
        point: &[("muvsq", 1.0), ("mursq", 1.0)],
        momenta: Momenta::None,
        references: vec![Reference::new(
            "alphaLoop, MATAD",
            ANALYTIC,
            &[
                (-1, "1.0", "0"),
                (0, "4.227843350984671393934879099176e-1", "0"),
            ],
        )],
    });
}

#[test]
fn one_loop_doubled_tadpole_with_two_external_momentum_projections() {
    check(Integral {
        name: "1L doubled tadpole, numerator (k·p1)^2",
        loops: 1,
        input: "(k(1,11)*p(1,11)*k(1,12)*p(1,12))*topo(prop(1,edge(1,1),k(1),muvsq,2))",
        scheme: LoopNormalizationFactor::MSbar,
        laurent_terms: 5,
        point: &[("muvsq", 1.0), ("mursq", 1.0)],
        momenta: Momenta::Small(1),
        references: vec![
            Reference::new(
                "alphaLoop, MATAD",
                ANALYTIC,
                &[
                    (-1, "0", "-9.892747067878755034798287150883e-3"),
                    (0, "0", "-9.892747067878755034798287150883e-3"),
                    (1, "0", "-1.802920540121208908815751124632e-2"),
                    (2, "0", "-1.406532376649359143630292751833e-2"),
                    (3, "0", "-2.008809563968960604519605888938e-2"),
                ],
            ),
            Reference::new(
                "pySecDec (Monte Carlo)",
                PYSECDEC,
                &[
                    (-1, "0", "-9.8927470662601990625262260437e-3"),
                    (0, "0", "-9.8927470662601990625262260437e-3"),
                    (1, "0", "-1.8029205397397163324058055878e-2"),
                    (2, "0", "-1.4065323763134074397385120392e-2"),
                    (3, "0", "-2.0088095635401259642094373703e-2"),
                ],
            ),
        ],
    });
}

#[test]
fn two_loop_sunrise() {
    check(Integral {
        name: "2L sunrise, three equal masses",
        loops: 2,
        input: "(1)*topo(
              prop(1,edge(1,2),k(1),muvsq,1)
             *prop(2,edge(1,2),k(2),muvsq,1)
             *prop(3,edge(2,1),k(1)+k(2),muvsq,1))",
        scheme: LoopNormalizationFactor::MSbar,
        laurent_terms: 4,
        point: &[("muvsq", 1.0), ("mursq", 1.0)],
        momenta: Momenta::None,
        references: vec![
            Reference::new(
                "alphaLoop, MATAD",
                ANALYTIC,
                &[
                    (-2, "-6.015223977354102649894208945578e-5", "0"),
                    (-1, "-1.804567193206230794968262683673e-4", "0"),
                    (0, "-3.790208765869208760159732351246e-4", "0"),
                    (1, "-9.708141619945564042866174113010e-4", "0"),
                ],
            ),
            Reference::new(
                "pySecDec (Monte Carlo)",
                PYSECDEC,
                &[
                    (-2, "-6.0152239775845828262390568852e-5", "0"),
                    (-1, "-1.8045671933464291214477270842e-4", "0"),
                    (0, "-3.7902087660768302157521247864e-4", "0"),
                    (1, "-9.7081416197397629730403423309e-4", "0"),
                ],
            ),
        ],
    });
}

#[test]
fn three_loop_mercedes() {
    check(Integral {
        name: "3L mercedes (tetrahedron), six equal masses",
        loops: 3,
        input: "(1)*topo(
              prop(1,edge(1,2),k(1),muvsq,1)
             *prop(2,edge(2,3),k(2),muvsq,1)
             *prop(3,edge(3,1),k(3),muvsq,1)
             *prop(4,edge(1,4),k(3)-k(1),muvsq,1)
             *prop(5,edge(2,4),k(1)-k(2),muvsq,1)
             *prop(6,edge(3,4),k(2)-k(3),muvsq,1))",
        scheme: LoopNormalizationFactor::MSbar,
        laurent_terms: 5,
        point: &[("muvsq", 1.0), ("mursq", 1.0)],
        momenta: Momenta::None,
        // Overall logarithmically divergent without subdivergences: a single pole.
        references: vec![Reference::new(
            "MATAD",
            ANALYTIC,
            &[
                (-3, "0", "0"),
                (-2, "0", "0"),
                (-1, "0", "-6.105142965702799027412986390958e-7"),
                (0, "0", "2.548415539172476714996527540285e-6"),
                (1, "0", "-1.063440725963320582103097545956e-5"),
            ],
        )],
    });
}

#[test]
fn three_loop_mercedes_with_rank_four_numerator() {
    check(Integral {
        name: "3L mercedes, rank-4 numerator with two external momenta",
        loops: 3,
        input: "(   k(1,11)*k(2,11)*k(1,22)*k(2,22)
              + p(1,11)*k(3,11)*k(3,22)*p(2,22)
              + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) )
            *topo(
              prop(1,edge(1,2),k(1),muvsq,1)
             *prop(2,edge(2,3),k(2),muvsq,1)
             *prop(3,edge(3,1),k(3),muvsq,1)
             *prop(4,edge(1,4),k(3)-k(1),muvsq,1)
             *prop(5,edge(2,4),k(1)-k(2),muvsq,1)
             *prop(6,edge(3,4),k(2)-k(3),muvsq,1))",
        scheme: LoopNormalizationFactor::MSbar,
        laurent_terms: 5,
        point: &[("muvsq", 1.0), ("mursq", 2.0)],
        momenta: Momenta::Large(2),
        references: vec![Reference::new(
            "MATAD",
            ANALYTIC,
            &[
                (-3, "0", "1.314137110381316809259740852289e-3"),
                (-2, "0", "1.027363035879144648390245670512e-2"),
                (-1, "0", "4.561360375915191500552582268384e-2"),
                (0, "0", "1.409480043924165037345830826758e-1"),
                (1, "0", "5.102719154490984941123728761288e-1"),
            ],
        )],
    });
}

#[test]
fn four_loop_h_topology_with_rank_four_numerator() {
    check(Integral {
        name: "4L H topology, rank-4 numerator with two external momenta",
        loops: 4,
        input: "(   k(1,11)*k(2,11)*k(1,22)*k(2,22)
              + p(1,11)*k(3,11)*k(3,22)*p(2,22)
              + p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) )
            *topo(
              prop(1,edge(5,1),k(1),muvsq,1)
             *prop(2,edge(2,6),k(2),muvsq,1)
             *prop(3,edge(6,5),k(3),muvsq,1)
             *prop(4,edge(3,4),k(4),muvsq,1)
             *prop(5,edge(4,5),k(1)-k(3),muvsq,1)
             *prop(6,edge(6,3),k(2)-k(3),muvsq,1)
             *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)
             *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)
             *prop(9,edge(1,2),k(3)+k(4),muvsq,1))",
        scheme: LoopNormalizationFactor::MSbar,
        laurent_terms: 5,
        point: &[("muvsq", 3.0), ("mursq", 5.0)],
        momenta: Momenta::Small(2),
        references: vec![
            Reference::new(
                "FMFT",
                FMFT,
                &[
                    (-4, "1.809145974886785501452557650622e-9", "0"),
                    (-3, "1.862208677723707446525921998109e-8", "0"),
                    (-2, "8.865577059648962609058045604434e-8", "0"),
                    (-1, "4.064224364096375531168559391726e-7", "0"),
                    (0, "1.849727263972748627996238794917e-6", "0"),
                ],
            ),
            // tests/inputs/h_rank_four_external_reference.json: exact FORM
            // reduction with the released external master values, 80 digits.
            Reference::new(
                "FORM + external masters",
                ANALYTIC,
                &[
                    (-4, "1.809145974886785501452557650621587342771e-9", "0"),
                    (-3, "1.862208677723707446525921998108647961948e-8", "0"),
                    (-2, "8.865577059648962609058045604434129517102e-8", "0"),
                    (-1, "4.064224364096375531168559391726011844697e-7", "0"),
                    (0, "1.849727263972748627996238788187829600708e-6", "0"),
                ],
            ),
        ],
    });
}

#[test]
fn four_loop_master_pr9d_from_the_x_topology() {
    check(Integral {
        name: "4L X topology with two lines removed and one doubled (FMFT master PR9d)",
        loops: 4,
        input: "(1)*topo(
              prop(1,edge(5,1),k(1),muvsq,1)
             *prop(2,edge(2,6),k(2),muvsq,1)
             *prop(3,edge(6,5),k(3),muvsq,0)
             *prop(4,edge(4,3),k(4),muvsq,2)
             *prop(5,edge(3,5),k(1)-k(3),muvsq,1)
             *prop(6,edge(6,4),k(2)-k(3),muvsq,1)
             *prop(7,edge(3,2),k(3)-k(1)+k(4),muvsq,1)
             *prop(8,edge(1,4),k(3)-k(2)+k(4),muvsq,1)
             *prop(9,edge(2,1),k(3)-k(1)-k(2)+k(4),muvsq,0))",
        scheme: LoopNormalizationFactor::FMFTandMATAD,
        laurent_terms: 5,
        point: &[("muvsq", 1.0), ("mursq", 1.0)],
        momenta: Momenta::None,
        references: vec![Reference::new(
            "FMFT",
            FMFT,
            &[
                (-4, "8.333333333333333333333333333333e-2", "0"),
                (-3, "3.333333333333333333333333333333e-1", "0"),
                (-2, "-3.144646082033583725553786166618e-1", "0"),
                (-1, "5.421352941798334340259377610275", "0"),
                (0, "-28.31064373017674207211847384976", "0"),
            ],
        )],
    });
}

#[test]
fn four_loop_clover_with_numerator() {
    check(Integral {
        name: "4L clover (four one-loop tadpoles), one doubled, rank-4 numerator",
        loops: 4,
        input: "(   user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
              + user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
              + user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22) )
            *topo(
              prop(1,edge(1,1),k(1),muvsq,2)
             *prop(2,edge(1,1),k(2),muvsq,1)
             *prop(3,edge(1,1),k(3),muvsq,1)
             *prop(4,edge(1,1),k(4),muvsq,1))",
        scheme: LoopNormalizationFactor::FMFTandMATAD,
        laurent_terms: 5,
        point: &[
            ("muvsq", 0.3),
            ("mursq", 0.7),
            ("user_space::A", 3.0),
            ("user_space::B", 4.0),
            ("user_space::C", 5.0),
        ],
        momenta: Momenta::Small(2),
        references: vec![Reference::new(
            "FMFT",
            FMFT,
            &[
                (-4, "-1.897149599999999855007182247846e-1", "0"),
                (-3, "-1.495259819655131009380566817668", "0"),
                (-2, "-6.805240907875078933713325181389", "0"),
                (-1, "-22.56027900679456203938234477552", "0"),
                (0, "-60.49337040949871593265194938449", "0"),
            ],
        )],
    });
}

/// The table's own semantics, without vakint: requested orders vakint leaves
/// out are exact zeros, and an order a reference does not list must vanish.
mod laurent_table_semantics {
    use super::reference_table::{LaurentTable, Reference, complex_from_decimal};

    fn table(values: &[(i64, &str)]) -> LaurentTable {
        LaurentTable {
            title: "synthetic 3L result".into(),
            facts: Vec::new(),
            engine: "vakint".into(),
            values: values
                .iter()
                .map(|(power, imaginary)| {
                    (
                        *power,
                        complex_from_decimal("0", imaginary).expect("decimal value"),
                    )
                })
                .collect(),
            orders: Some(-3..=1),
            shown_digits: 8,
        }
    }

    const POLE_AND_FINITE: &[(i64, &str)] = &[(-1, "-6.1e-7"), (0, "2.5e-6"), (1, "-1.1e-5")];

    fn matad(coefficients: &[(i64, &str, &str)]) -> Vec<Reference> {
        vec![Reference::new("MATAD", 1e-28, coefficients)]
    }

    #[test]
    fn explicit_zero_references_match_orders_vakint_leaves_out() {
        let references = matad(&[
            (-3, "0", "0"),
            (-2, "0", "0"),
            (-1, "0", "-6.1e-7"),
            (0, "0", "2.5e-6"),
            (1, "0", "-1.1e-5"),
        ]);
        assert_eq!(table(POLE_AND_FINITE).check(&references), Ok(()));
    }

    #[test]
    fn a_spurious_pole_fails_against_an_unlisted_reference_order() {
        let mut values = POLE_AND_FINITE.to_vec();
        values.push((-2, "3e-9"));
        let references = matad(&[
            (-1, "0", "-6.1e-7"),
            (0, "0", "2.5e-6"),
            (1, "0", "-1.1e-5"),
        ]);
        let failure = table(&values)
            .check(&references)
            .expect_err("a nonzero ε^-2 without a reference value must fail");
        assert!(failure.contains("ε^-2 vs MATAD: Δabs"), "{failure}");
    }

    #[test]
    fn a_spurious_pole_fails_against_an_explicit_zero_reference() {
        let mut values = POLE_AND_FINITE.to_vec();
        values.push((-3, "1e-20"));
        let references = matad(&[
            (-3, "0", "0"),
            (-2, "0", "0"),
            (-1, "0", "-6.1e-7"),
            (0, "0", "2.5e-6"),
            (1, "0", "-1.1e-5"),
        ]);
        let failure = table(&values)
            .check(&references)
            .expect_err("a nonzero ε^-3 against a zero reference must fail");
        assert!(failure.contains("ε^-3 vs MATAD: Δabs"), "{failure}");
    }

    #[test]
    fn a_nonzero_reference_fails_when_vakint_leaves_the_order_out() {
        let references = matad(&[(-2, "0", "4e-8"), (-1, "0", "-6.1e-7"), (0, "0", "2.5e-6")]);
        let failure = table(POLE_AND_FINITE)
            .check(&references)
            .expect_err("vakint's implied zero at ε^-2 differs from the reference");
        assert!(failure.contains("ε^-2 vs MATAD: Δrel"), "{failure}");
    }

    #[test]
    fn a_reference_order_outside_the_requested_orders_is_not_computed() {
        let references = matad(&[(-1, "0", "-6.1e-7"), (2, "0", "0")]);
        let failure = table(POLE_AND_FINITE)
            .check(&references)
            .expect_err("vakint was not asked for ε^2");
        assert!(
            failure.contains("ε^2 vs MATAD: vakint did not compute this order"),
            "{failure}"
        );
    }
}
