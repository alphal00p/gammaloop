//! The literal four-loop numerical acceptance inputs, independent of reducers.
//!
//! These records mirror the existing analytic and decorated-input tests. Parent
//! labels select caller-supplied test descriptors, never engine dispatch. In
//! particular a pinched or disconnected input is not replaced by its parent.

use symbolica::atom::Atom;
use vakint::{LoopNormalizationFactor, Vakint, VakintSettings, vakint_parse};

pub type ExternalMomentum = (usize, (f64, f64, f64, f64));
pub type LaurentReference = (i64, (String, String));

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ParentDescriptor {
    H,
    X,
    Bmw,
    Fg,
}

impl ParentDescriptor {
    pub fn csv(self) -> &'static str {
        match self {
            Self::H => include_str!("../inputs/experimental_four_loop_h.csv"),
            Self::X => include_str!("../inputs/experimental_four_loop_x.csv"),
            Self::Bmw => include_str!("../inputs/experimental_four_loop_bmw.csv"),
            Self::Fg => include_str!("../inputs/experimental_four_loop_fg.csv"),
        }
    }
}

pub struct AcceptanceInput {
    pub name: &'static str,
    pub parent: ParentDescriptor,
    pub input: Atom,
    pub settings: VakintSettings,
    /// Includes spectator symbols and free-index components, not only masses.
    pub parameters: Vec<(&'static str, f64)>,
    pub external_momenta: Vec<ExternalMomentum>,
    pub expected: Vec<LaurentReference>,
}

impl AcceptanceInput {
    fn new(
        name: &'static str,
        parent: ParentDescriptor,
        input: &str,
        normalization: LoopNormalizationFactor,
        parameters: &[(&'static str, f64)],
        reference: &[(i64, &str)],
    ) -> Self {
        Self {
            name,
            parent,
            input: vakint_parse!(input).unwrap_or_else(|error| panic!("{name}: {error}")),
            settings: VakintSettings {
                integral_normalization_factor: normalization,
                number_of_terms_in_epsilon_expansion: 5,
                run_time_decimal_precision: 32,
                ..VakintSettings::default()
            },
            parameters: parameters.to_vec(),
            // Keep the upstream floating-input boundary and arithmetic exactly,
            // rather than rounding the resulting components to short decimals.
            external_momenta: (1..=2)
                .map(|i| {
                    (
                        i,
                        (
                            0.17 * (i + 1) as f64,
                            0.4 * (i + 2) as f64,
                            0.3 * (i + 3) as f64,
                            0.12 * (i + 4) as f64,
                        ),
                    )
                })
                .collect(),
            expected: reference
                .iter()
                .map(|(power, real)| (*power, ((*real).into(), "0.0".into())))
                .collect(),
        }
    }
}

const H: &str = "topo(
    prop(1,edge(5,1),k(1),muvsq,1)*prop(2,edge(2,6),k(2),muvsq,1)*
    prop(3,edge(6,5),k(3),muvsq,1)*prop(4,edge(3,4),k(4),muvsq,1)*
    prop(5,edge(4,5),k(1)-k(3),muvsq,1)*prop(6,edge(6,3),k(2)-k(3),muvsq,1)*
    prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)*
    prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)*
    prop(9,edge(1,2),k(3)+k(4),muvsq,1))";

const RANK_FOUR: &str = "k(1,11)*k(2,11)*k(1,22)*k(2,22)
    +p(1,11)*k(3,11)*k(3,22)*p(2,22)
    +p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)";
const DECORATED_RANK_FOUR: &str = "user_space::A*k(1,11)*k(2,11)*k(1,22)*k(2,22)
    +user_space::B*p(1,11)*k(3,11)*k(3,22)*p(2,22)
    +user_space::C*p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)";

const PR9D_H: &str = "topo(
    prop(1,edge(5,1),k(1),muvsq,1)*prop(2,edge(2,6),k(2),muvsq,1)*
    prop(3,edge(6,5),k(3),muvsq,0)*prop(4,edge(3,4),k(4),muvsq,2)*
    prop(5,edge(4,5),k(1)-k(3),muvsq,1)*prop(6,edge(6,3),k(2)-k(3),muvsq,1)*
    prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)*
    prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)*
    prop(9,edge(1,2),k(3)+k(4),muvsq,0))";
const PR9D_X: &str = "topo(
    prop(1,edge(5,1),k(1),muvsq,1)*prop(2,edge(2,6),k(2),muvsq,1)*
    prop(3,edge(6,5),k(3),muvsq,0)*prop(4,edge(4,3),k(4),muvsq,2)*
    prop(5,edge(3,5),k(1)-k(3),muvsq,1)*prop(6,edge(6,4),k(2)-k(3),muvsq,1)*
    prop(7,edge(3,2),k(3)-k(1)+k(4),muvsq,1)*
    prop(8,edge(1,4),k(3)-k(2)+k(4),muvsq,1)*
    prop(9,edge(2,1),k(3)-k(1)-k(2)+k(4),muvsq,0))";
const PR9D_H_PINCH: &str = "topo(
    prop(1,edge(5,2),k(1),muvsq,1)*prop(2,edge(2,5),k(2),muvsq,1)*
    prop(4,edge(3,4),k(3),muvsq,2)*prop(5,edge(4,5),k(4),muvsq,1)*
    prop(6,edge(5,3),k(4)+k(2)-k(1),muvsq,1)*
    prop(7,edge(4,2),k(3)-k(4),muvsq,1)*
    prop(8,edge(2,3),k(1)-k(2)+k(3)-k(4),muvsq,1))";
const PR9D_FG: &str = "topo(
    prop(1,edge(5,3),k(1),muvsq,1)*prop(2,edge(3,4),k(2),muvsq,2)*
    prop(3,edge(4,5),k(3),muvsq,1)*prop(4,edge(2,1),k(1)-k(3),muvsq,0)*
    prop(5,edge(5,1),k(4),muvsq,1)*prop(6,edge(4,2),k(2)-k(3),muvsq,1)*
    prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)*
    prop(8,edge(3,2),k(1)-k(2),muvsq,1))";
const PR9D_FG_PINCH: &str = "topo(
    prop(1,edge(5,3),k(1),muvsq,1)*prop(2,edge(3,4),k(2),muvsq,2)*
    prop(3,edge(4,5),k(3),muvsq,1)*prop(5,edge(5,1),k(4),muvsq,1)*
    prop(6,edge(4,1),k(2)-k(3),muvsq,1)*
    prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)*
    prop(8,edge(3,1),k(1)-k(2),muvsq,1))";
const PR11D: &str = "topo(
    prop(1,edge(1,2),k(1),muvsq,2)*prop(2,edge(2,5),k(2),muvsq,1)*
    prop(3,edge(3,4),k(3),muvsq,1)*prop(4,edge(4,5),k(4),muvsq,1)*
    prop(5,edge(2,3),k(1)-k(2),muvsq,1)*prop(6,edge(4,1),k(3)-k(4),muvsq,1)*
    prop(7,edge(5,3),k(2)+k(3)-k(1),muvsq,1)*
    prop(8,edge(1,5),k(3)-k(4)-k(1),muvsq,1))";

const CLOVER: &str = "topo(
    prop(1,edge(1,1),k(1),muvsq,1)*prop(2,edge(1,1),k(2),muvsq,1)*
    prop(3,edge(1,1),k(3),muvsq,1)*prop(4,edge(1,1),k(4),muvsq,1))";
const DOTTED_CLOVER: &str = "topo(
    prop(1,edge(1,1),k(1),muvsq,2)*prop(2,edge(1,1),k(2),muvsq,1)*
    prop(3,edge(1,1),k(3),muvsq,1)*prop(4,edge(1,1),k(4),muvsq,1))";

const FREEFORM: &str = "(
    (user_space::BBsigma(user_space::some_args)
     +user_space::{symmetric,scalar}::BBsigma2(user_space::{real}::some_args2)
     +user_space::{integer}::BBparam)
    *vakint::p(1,user_space::mink4(4,33))*vakint::p(2,user_space::mink4(4,33))
    *vakint::p(1,user_space::mink4(4,11))*vakint::p(2,user_space::mink4(4,22))
    +vakint::k(3,user_space::mink4(4,11))*vakint::k(3,user_space::mink4(4,22))
    +vakint::k(3,user_space::mink4(4,77))*vakint::p(1,user_space::mink4(4,77)))
    *vakint::topo(
    vakint::prop(9,vakint::edge(66,66),vakint::k(1),user_space::{real}::BBMUVsq,1)*
    vakint::prop(9,vakint::edge(66,66),vakint::k(2),user_space::{real}::BBMUVsq,1)*
    vakint::prop(9,vakint::edge(66,66),vakint::k(3),user_space::{real}::BBMUVsq,1)*
    vakint::prop(9,vakint::edge(66,66),vakint::k(4),user_space::{real}::BBMUVsq,1))";

/// All fourteen analytic four-loop entries and the decorated four-loop entry
/// whose historical upstream name happens to contain `1l`.
pub fn cases() -> Vec<AcceptanceInput> {
    Vakint::initialize_vakint_symbols();
    use LoopNormalizationFactor::{FMFTandMATAD, MSbar};
    use ParentDescriptor::{Bmw, Fg, H as ParentH, X};
    let unit = &[("muvsq", 1.0), ("mursq", 1.0)];
    let h_reference = &[(0, "-2.169283452273432986058475530569e-9")];
    let mut h = AcceptanceInput::new("test_integrate_4l_h", ParentH, H, MSbar, unit, h_reference);
    let mut squared = AcceptanceInput::new(
        "test_integrate_4l_h_squared_mass",
        ParentH,
        &format!(
            "(2*user_space::muv-user_space::muv^2)*{}",
            H.replace("muvsq", "user_space::muv^2")
        ),
        MSbar,
        &[("user_space::muv", 1.0), ("mursq", 1.0)],
        h_reference,
    );
    for entry in [&mut h, &mut squared] {
        entry.external_momenta = (1..=1)
            .map(|i| {
                (
                    i,
                    (
                        17.0 * (i + 1) as f64,
                        4.0 * (i + 2) as f64,
                        3.0 * (i + 3) as f64,
                        12.0 * (i + 4) as f64,
                    ),
                )
            })
            .collect();
    }
    let mut result = vec![
        h,
        squared,
        AcceptanceInput::new(
            "test_integrate_4l_h_rank_4",
            ParentH,
            &format!("({RANK_FOUR})*{H}"),
            MSbar,
            &[("muvsq", 3.0), ("mursq", 5.0)],
            &[
                (-4, "1.809145974886785501452557650622e-9"),
                (-3, "1.862208677723707446525921998109e-8"),
                (-2, "8.865577059648962609058045604434e-8"),
                (-1, "4.064224364096375531168559391726e-7"),
                (0, "7.705260630861442737917312763495e-6"),
            ],
        ),
        AcceptanceInput::new(
            "test_integrate_4l_h_rank_4_additional_symbols_numerator",
            ParentH,
            &format!("({DECORATED_RANK_FOUR})*{H}"),
            MSbar,
            &[
                ("muvsq", 3.0),
                ("mursq", 5.0),
                ("user_space::A", 5.0),
                ("user_space::B", 7.0),
                ("user_space::C", 11.0),
            ],
            &[
                (-4, "9.045729874433927507262788253108e-9"),
                (-3, "9.311043388618537232629609990544e-8"),
                (-2, "3.740608759535330337057785114301e-7"),
                (-1, "2.152059306128503126262830041472e-6"),
                (0, "6.989489067178944692297775010846e-5"),
            ],
        ),
    ];
    let pr9d = &[
        (-4, "8.333333333333333333333333333333e-2"),
        (-3, "3.333333333333333333333333333333e-1"),
        (-2, "-3.144646082033583725553786166618e-1"),
        (-1, "5.421352941798334340259377610275"),
        (0, "-28.31064373017674207211847384976"),
    ];
    for (name, parent, input) in [
        ("test_integrate_4l_PR9d_from_H", ParentH, PR9D_H),
        ("test_integrate_4l_PR9d_from_X", X, PR9D_X),
        ("test_integrate_4l_PR9d_from_H_pinch", ParentH, PR9D_H_PINCH),
        ("test_integrate_4l_PR9d_from_FG", Fg, PR9D_FG),
        ("test_integrate_4l_PR9d_from_FG_pinch", Fg, PR9D_FG_PINCH),
    ] {
        result.push(AcceptanceInput::new(
            name,
            parent,
            input,
            FMFTandMATAD,
            unit,
            pr9d,
        ));
    }
    result.extend([
        AcceptanceInput::new(
            "test_integrate_4l_PR11d",
            Bmw,
            PR11D,
            VakintSettings::default().integral_normalization_factor,
            unit,
            &[(0, "-2.906486288643112641819206002127")],
        ),
        AcceptanceInput::new(
            "test_integrate_4l_clover",
            Fg,
            CLOVER,
            FMFTandMATAD,
            unit,
            &[
                (-4, "1.000000000000000000000000000000"),
                (-3, "4.000000000000000000000000000000"),
                (-2, "13.28986813369645287294483033329"),
                (-1, "31.55672999723968577791300378449"),
                (0, "67.98165058904685502307905531744"),
            ],
        ),
        AcceptanceInput::new(
            "test_integrate_4l_clover_with_non_unit_scales",
            Fg,
            CLOVER,
            VakintSettings::default().integral_normalization_factor,
            &[("muvsq", 3.0), ("mursq", 7.0)],
            &[
                (-4, "81.00000000000000000000000000000"),
                (-3, "-218.9682569565641868485693739496"),
                (-2, "724.4490568207455247307151297051"),
                (-1, "-1446.834846767122729283863130264"),
                (0, "3106.843653546628093141699453745"),
            ],
        ),
        AcceptanceInput::new(
            "test_integrate_4l_dotted_clover",
            Fg,
            DOTTED_CLOVER,
            VakintSettings::default().integral_normalization_factor,
            &[("muvsq", 3.0), ("mursq", 1.0)],
            &[
                (-4, "27.00000000000000000000000000000"),
                (-3, "-99.98941898552139561618979131653"),
                (-2, "314.4724379257699038597615012182"),
                (-1, "-723.7613011959560846715260866565"),
                (0, "1517.892833437916940808520861336"),
            ],
        ),
        AcceptanceInput::new(
            "test_integrate_4l_clover_with_numerator",
            Fg,
            &format!("({DECORATED_RANK_FOUR})*{DOTTED_CLOVER}"),
            FMFTandMATAD,
            &[
                ("muvsq", 0.3),
                ("mursq", 0.7),
                ("user_space::A", 3.0),
                ("user_space::B", 4.0),
                ("user_space::C", 5.0),
            ],
            &[
                (-4, "-1.897149599999999855007182247846e-1"),
                (-3, "-1.495259819655131009380566817668"),
                (-2, "-6.805240907875078933713325181389"),
                (-1, "-22.56027900679456203938234477552"),
                (0, "-60.49337040949871593265194938449"),
            ],
        ),
    ]);
    let mut decorated = AcceptanceInput::new(
        "test_integrate_1l_decorated_indices_fmft",
        Fg,
        FREEFORM,
        MSbar,
        &[
            ("user_space::{real}::BBMUVsq", 1.0),
            ("some_space::{real,scalar}::BBmursq", 1.0),
            (
                "vakint::g(user_space::mink4(4,22),user_space::mink4(4,11))",
                1.0,
            ),
            ("vakint::p(1,user_space::mink4(4,11))", 1.0),
            ("vakint::p(2,user_space::mink4(4,22))", 1.0),
            ("user_space::BBsigma(user_space::some_args)", 1.0),
            (
                "user_space::{symmetric,scalar}::BBsigma2(user_space::{real}::some_args2)",
                0.0,
            ),
            ("user_space::{integer}::BBparam", 0.0),
        ],
        &[
            (-4, "-5.996072606189216818994270605252e-9"),
            (-3, "-2.378327420532500222026013157094e-8"),
            (-2, "-7.878244126888092041321958350898e-8"),
            (-1, "-1.860926787346530261093067881728e-7"),
            (0, "-3.997176154874809449646384405081e-7"),
        ],
    );
    decorated.settings.mu_r_sq_symbol = "some_space::{real,scalar}::BBmursq".into();
    result.push(decorated);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const ANALYTIC: &str = include_str!("../integral_evaluation_analytic_tests.rs");
    const FREEFORM_SOURCE: &str = include_str!("../integral_evaluation_freeform_tests.rs");

    // Narrow source-fixture check, not a Rust parser: the upstream inputs are
    // ordinary string literals with only Rust's escaped-newline continuation.
    fn original_body(name: &str) -> &str {
        let source = if name == "test_integrate_1l_decorated_indices_fmft" {
            FREEFORM_SOURCE
        } else {
            ANALYTIC
        };
        source
            .split_once(&format!("fn {name}()"))
            .unwrap_or_else(|| panic!("missing original acceptance {name}"))
            .1
            .split("\nfn ")
            .next()
            .unwrap()
    }

    fn original_input(body: &str) -> String {
        let macro_name = if body.contains("vakint_parse!(") {
            "vakint_parse!("
        } else {
            "try_parse!("
        };
        let literal = body
            .split_once(macro_name)
            .unwrap()
            .1
            .split_once('"')
            .unwrap()
            .1;
        literal.split('"').next().unwrap().replace("\\\n", "")
    }

    #[test]
    fn all_literal_inputs_and_numerical_settings_match_existing_acceptance() {
        let cases = cases();
        assert_eq!(cases.len(), 15);
        assert_eq!(ANALYTIC.matches("fn test_integrate_4l_").count(), 14);
        for case in cases {
            let body = original_body(case.name);
            let original = vakint_parse!(original_input(body)).unwrap();
            assert_eq!(case.input, original, "literal input drift: {}", case.name);
            let compact: String = body.chars().filter(|c| !c.is_whitespace()).collect();
            for (parameter, value) in &case.parameters {
                assert!(
                    compact.contains(&format!("(\"{parameter}\".into(),{value:?})")),
                    "parameter drift: {}: {parameter}={value}",
                    case.name
                );
            }
            for (power, (real, imaginary)) in &case.expected {
                assert!(
                    compact.contains(&format!(
                        "({power},(\"{real}\".into(),\"{imaginary}\".into())"
                    )),
                    "reference drift: {} epsilon^{power}",
                    case.name
                );
            }
            let expected_normalization = if body.contains("LoopNormalizationFactor::MSbar") {
                LoopNormalizationFactor::MSbar
            } else if body.contains("LoopNormalizationFactor::FMFTandMATAD") {
                LoopNormalizationFactor::FMFTandMATAD
            } else {
                VakintSettings::default().integral_normalization_factor
            };
            assert!(
                normalization_matches(
                    &case.settings.integral_normalization_factor,
                    &expected_normalization
                ),
                "normalization drift: {}",
                case.name
            );
            assert_eq!(case.settings.number_of_terms_in_epsilon_expansion, 5);
            assert_eq!(case.settings.run_time_decimal_precision, 32);
            assert!(!case.parent.csv().is_empty());
            if case.name == "test_integrate_4l_h" || case.name == "test_integrate_4l_h_squared_mass"
            {
                assert_eq!(case.external_momenta, vec![(1, (34.0, 12.0, 12.0, 60.0))]);
                assert!(compact.contains("&(1..=1).map(|i|(i,(17.0*((i+1)asf64),4.0*((i+2)asf64),3.0*((i+3)asf64),12.0*((i+4)asf64))))"));
            } else {
                assert_eq!(case.external_momenta.len(), 2);
                assert!(compact.contains("&(1..=2).map(|i|(i,(0.17*((i+1)asf64),0.4*((i+2)asf64),0.3*((i+3)asf64),0.12*((i+4)asf64))))"));
            }
        }
    }

    fn normalization_matches(
        actual: &LoopNormalizationFactor,
        expected: &LoopNormalizationFactor,
    ) -> bool {
        match (actual, expected) {
            (LoopNormalizationFactor::MSbar, LoopNormalizationFactor::MSbar)
            | (LoopNormalizationFactor::FMFTandMATAD, LoopNormalizationFactor::FMFTandMATAD)
            | (LoopNormalizationFactor::pySecDec, LoopNormalizationFactor::pySecDec) => true,
            _ => false,
        }
    }
}
