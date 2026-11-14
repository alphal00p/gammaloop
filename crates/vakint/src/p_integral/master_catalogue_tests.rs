//! The p master catalogue guards, mirroring the vacuum master catalogue tests
//! (`rustred_evaluation::masters` tests and `tests/master_truncation_guard_tests.rs`)
//! one to one on a small synthetic five-loop catalogue ([`fixture`]), and the
//! Forcer table run through the catalogue path on the frozen four-loop
//! battery (real data, the same code path).
//!
//! | vacuum test                                                            | here                                                         |
//! |------------------------------------------------------------------------|--------------------------------------------------------------|
//! | `missing_coefficients_are_not_assumed_zero`, truncation guard tests   | [`five_loop_p_master_known_through_refuses_deeper_orders`]  |
//! | `spurious_poles_require_deeper_master_coefficients`                   | [`spurious_pole_on_a_5l_p_master_needs_deeper_coefficients`] |
//! | `exact_cancellation_removes_unavailable_coefficients_and_unknowns`    | [`exact_cancellation_removes_unavailable_5l_coefficients`]   |
//! | `retained_catalogue_terminals_are_not_unknown_but_substitution_needs_coefficients` | [`missing_5l_p_master_coefficients_are_not_assumed_zero`] |
//! | `embedded_float_annotations_do_not_override_validated_precision`, finite constant without digits | [`numerical_5l_master_without_validated_digits_is_refused`] |
//! | `malformed_or_nonfinite_numerical_coefficients_report_the_master`     | [`malformed_catalogues_are_refused_with_the_master_and_order`] |
//! | Retain only warns; Numerical refuses unresolved masters               | [`unmapped_residual_is_unevaluated_in_symbolic_and_numerical`] |

use serde_json::{Value, json};
use symbolica::atom::{Atom, AtomCore};

use super::super::ReducedPropagatorValue;
use super::super::finalize::finalize_with;
use super::super::forcer::fixtures::{Finalized, ForcerFixture, I7_TERMS, finalize_settings};
use super::*;
use crate::LoopNormalizationFactor;
use crate::VakintSettings;
use crate::truncation_guard::expect_missing_orders;

/// The request named in every message of these tests.
const REQUEST: &str = "test request topo(T5(p(1),1,...))";

/// A synthetic five-loop catalogue: `syn_m1` known through `ε^0`
/// (`1/ε + 5/2`), `syn_num` numerical only (30 digits), `syn_empty` declared
/// without coefficients, `syn_float` an exact formula embedding a decimal
/// (5 validated digits), and terminal rows of every status. Shared with the
/// p path's five-loop guard tests (`rustred_evaluation::p_integral` tests).
pub(crate) fn fixture() -> Value {
    json!({
        "schema": P_MASTER_CATALOGUE_SCHEMA,
        "version": P_MASTER_CATALOGUE_VERSION,
        "data_version": "test-5L-v1",
        "normalisation": P_MASTER_NORMALISATION,
        "masters": [
            {"id": "syn-m1", "symbol": "syn_m1", "loops": 5, "graph": "syn-graph-1",
             "leading_order": -1, "unknown_from": 1,
             "coefficients": {"-1": {"exact": "1"}, "0": {"exact": "5/2"}}},
            {"id": "syn-num", "symbol": "syn_num", "loops": 5, "graph": "syn-graph-2",
             "leading_order": 0, "known_through": 0,
             "coefficients": {"0": {"numerical": "1.25", "significant_digits": 30}}},
            {"id": "syn-empty", "symbol": "syn_empty", "loops": 5, "graph": "syn-graph-3",
             "leading_order": -5},
            {"id": "syn-float", "symbol": "syn_float", "loops": 5, "graph": "syn-graph-4",
             "leading_order": 0, "significant_digits": 5,
             "coefficients": {"0": {"exact": "2.5`100"}}},
            {"id": "syn-zeta", "symbol": "syn_zeta", "loops": 5, "graph": "syn-graph-5",
             "leading_order": 0, "coefficients": {"0": {"exact": "zeta(3)+pi^2-MZV(5,3)"}}}
        ],
        "terminals": [
            {"node": "T900", "key": [1, 1, 0], "status": "mapped",
             "combination": [{"master": "syn-m1", "coefficient": "1"}]},
            {"node": "T900", "key": [1, 1, 1], "status": "mapped",
             "combination": [{"master": "syn-m1", "coefficient": "1/(d-4)"}]},
            {"node": "T900", "key": [2, 1, 0], "status": "mapped",
             "combination": [{"master": "syn_m1", "coefficient": "2"},
                             {"master": null, "coefficient": "PG(0,0)^5"}]},
            {"node": "T901", "key": [2], "status": "unmapped", "note": "WP3 row pending"},
            {"node": "T902", "key": [1], "status": "additional", "note": "spoke dot"},
            {"node": "T903", "key": [1], "status": "mapped",
             "combination": [{"master": "syn-empty", "coefficient": "1"}]},
            {"node": "T904", "key": [1], "status": "mapped",
             "combination": [{"master": "syn-num", "coefficient": "1"}]},
            {"node": "T905", "key": [1], "status": "mapped",
             "combination": [{"master": "syn-float", "coefficient": "1"}]},
            {"node": "T906", "key": [1], "status": "mapped",
             "combination": [{"master": "syn-zeta", "coefficient": "1"}]}
        ]
    })
}

fn load(value: &Value) -> Result<PMasterCatalogue, PMasterCatalogueError> {
    PMasterCatalogue::from_json("test catalogue", &value.to_string())
}

fn catalogue() -> PMasterCatalogue {
    load(&fixture()).unwrap()
}

/// `PResidual(node, key..)`.
fn residual(node: &str, key: &[i64]) -> Atom {
    let mut arguments = vec![Atom::var(vk_symbol!(node))];
    arguments.extend(key.iter().map(|&k| Atom::num(k)));
    p_residual_symbol().call(&arguments[..])
}

/// G-scheme settings: at `Q² = μ² = 1` the value is `reduced` itself.
fn settings(terms: i64, digits: u32) -> VakintSettings {
    VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        run_time_decimal_precision: digits,
        ..VakintSettings::default()
    }
}

/// Finalize the five-loop `reduced` (scale power 0) with `catalogue`.
fn run(
    catalogue: Option<&PMasterCatalogue>,
    reduced: Atom,
    mode: RustRedMasterMode,
    terms: i64,
    digits: u32,
) -> Result<Atom, VakintError> {
    crate::Vakint::initialize_vakint_symbols();
    finalize_with(
        &[ReducedPropagatorValue {
            scale_power: 0,
            reduced,
        }],
        5,
        &vk_parse!("p(1)").unwrap(),
        &settings(terms, digits),
        mode,
        "RustRed",
        PMasterSource::new(catalogue, REQUEST),
    )
}

/// The Laurent coefficients `ε^first..=ε^last` at `Q² = μ² = 1`.
fn at_unit_scale(value: &Atom, first: i64, last: i64) -> Vec<Atom> {
    super::super::finalize::tests::coefficients_at_unit_scale(
        &value
            .replace(Atom::var(vk_symbol!("mursq")).to_pattern())
            .with(Atom::num(1).to_pattern()),
        &settings(last + 6, 40),
        &vk_parse!("p(1)").unwrap(),
        first,
        last,
    )
}

/// The default epsilon (`vakint::ε`).
fn eps() -> Atom {
    Atom::var(vk_symbol!(&VakintSettings::default().epsilon_symbol))
}

fn catalogue_error(error: VakintError) -> PMasterCatalogueError {
    match error {
        VakintError::PMasterCatalogue(error) => error,
        other => panic!("expected a p master catalogue error, got: {other}"),
    }
}

fn assert_names(message: &str, parts: &[&str]) {
    for part in parts {
        assert!(message.contains(part), "{part:?} missing in: {message}");
    }
}

#[test]
fn the_synthetic_catalogue_loads_with_its_depths() {
    let catalogue = catalogue();
    assert_eq!(catalogue.counts(), (5, 9));
    assert_eq!(catalogue.data_version(), "test-5L-v1");
    let m1 = catalogue.master("syn_m1").unwrap();
    assert_eq!((m1.leading_order, m1.known_through), (-1, 0));
    assert_eq!(catalogue.master("syn_empty").unwrap().known_through, -6);
    assert!(catalogue.to_string().contains("test-5L-v1"));
}

/// ← `missing_coefficients_are_not_assumed_zero` and
/// `master_truncation_guard_tests` (`known_through`): through `ε^0` the
/// result is the known `1/ε + 5/2`; one order more names the master, its
/// first missing and last needed order, the catalogue and the request.
#[test]
fn five_loop_p_master_known_through_refuses_deeper_orders() {
    let catalogue = catalogue();
    let terminal = residual("T900", &[1, 1, 0]);
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let value = run(Some(&catalogue), terminal.clone(), mode, 6, 40).unwrap();
        assert_eq!(
            at_unit_scale(&value, -1, 0),
            [Atom::num(1), Atom::num((5, 2))],
            "{mode}"
        );
        let missing = expect_missing_orders(
            run(Some(&catalogue), terminal.clone(), mode, 7, 40).unwrap_err(),
        );
        assert_eq!(missing.requested_order, Some(1));
        let [master] = missing.masters.as_slice() else {
            panic!("one master lacks orders: {missing}");
        };
        assert!(
            master.master.starts_with("PMaster(syn_m1) [syn-m1, 5-loop"),
            "{}",
            master.master
        );
        assert_eq!(
            (
                master.first_missing_order,
                master.last_needed_order,
                master.first_affected_result_order
            ),
            (1, Some(1), 1)
        );
        assert_names(
            &missing.to_string(),
            &["test-5L-v1", REQUEST, "ε^1", "known ε^-1..ε^0"],
        );
    }
}

/// ← `spurious_poles_require_deeper_master_coefficients`: the map row's
/// `1/(d-4) = -1/(2ε)` drags the unknown `ε^1` of the master into `ε^0`.
#[test]
fn spurious_pole_on_a_5l_p_master_needs_deeper_coefficients() {
    let catalogue = catalogue();
    let missing = expect_missing_orders(
        run(
            Some(&catalogue),
            residual("T900", &[1, 1, 1]),
            RustRedMasterMode::Symbolic,
            6,
            40,
        )
        .unwrap_err(),
    );
    assert_eq!(missing.requested_order, Some(0));
    assert!(missing.backend.contains("test-5L-v1"), "{missing}");
    let [master] = missing.masters.as_slice() else {
        panic!("one master lacks orders: {missing}");
    };
    assert!(master.master.starts_with("PMaster(syn_m1)"));
    assert_eq!(master.first_missing_order, 1);
    assert_eq!(master.last_needed_order, Some(1));
    assert_eq!(master.first_affected_result_order, 0);
    // One order less is known: -1/(2ε)·(1/ε + 5/2).
    let value = run(
        Some(&catalogue),
        residual("T900", &[1, 1, 1]),
        RustRedMasterMode::Symbolic,
        5,
        40,
    )
    .unwrap();
    assert_eq!(
        at_unit_scale(&value, -2, -1),
        [Atom::num((-1, 2)), Atom::num((-5, 4))]
    );
}

/// ← `exact_cancellation_removes_unavailable_coefficients_and_unknowns`:
/// an exact cancellation needs no master order, at any depth.
#[test]
fn exact_cancellation_removes_unavailable_5l_coefficients() {
    let catalogue = catalogue();
    let master = p_master("syn_m1");
    let reduced = (residual("T900", &[1, 1, 0]) - &master) * eps().pow(Atom::num(-9));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        assert!(
            run(Some(&catalogue), reduced.clone(), mode, 12, 40)
                .unwrap()
                .is_zero()
        );
    }
    // A mapped row with a tag: 2·syn_m1 + PG(0,0)^5 = 2/ε + 5 + 1/ε^5 + ...
    let value = run(
        Some(&catalogue),
        residual("T900", &[2, 1, 0]) - Atom::num(2) * &master,
        RustRedMasterMode::Symbolic,
        12,
        40,
    )
    .unwrap();
    assert_eq!(at_unit_scale(&value, -5, -5), [Atom::num(1)]);
}

/// ← `retained_catalogue_terminals_are_not_unknown_but_substitution_needs_coefficients`
/// and `missing_coefficients_are_not_assumed_zero`: a master declared
/// without coefficients is kept by Retain and refused by the expanding modes,
/// naming its id, graph, the catalogue and the request; an undeclared master
/// is refused by name.
#[test]
fn missing_5l_p_master_coefficients_are_not_assumed_zero() {
    let catalogue = catalogue();
    let terminal = residual("T903", &[1]);
    let retained = run(
        Some(&catalogue),
        terminal.clone(),
        RustRedMasterMode::Retain,
        6,
        40,
    )
    .unwrap();
    assert!(retained.contains_symbol(p_master_symbol()));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let error =
            catalogue_error(run(Some(&catalogue), terminal.clone(), mode, 6, 40).unwrap_err());
        assert!(
            matches!(error, PMasterCatalogueError::MissingCoefficients { .. }),
            "{error}"
        );
        assert_names(
            &error.to_string(),
            &[
                "syn-empty",
                "declared terminal has no computed coefficients",
                "syn-graph-3",
                "test-5L-v1",
                REQUEST,
            ],
        );
    }
    let undeclared = p_master("syn_unknown");
    let error = catalogue_error(
        run(
            Some(&catalogue),
            undeclared,
            RustRedMasterMode::Symbolic,
            6,
            40,
        )
        .unwrap_err(),
    );
    assert!(
        matches!(error, PMasterCatalogueError::UndeclaredMasters { .. }),
        "{error}"
    );
    assert_names(
        &error.to_string(),
        &["PMaster(syn_unknown)", "test-5L-v1", REQUEST],
    );
}

/// ← `embedded_float_annotations_do_not_override_validated_precision` and
/// the "finite constant without validated significant digits" refusal: a
/// numerical-only coefficient needs validated digits (none: an error; fewer
/// than `run_time_decimal_precision`: a warning, as on the vacuum path);
/// Symbolic mode refuses it.
#[test]
fn numerical_5l_master_without_validated_digits_is_refused() {
    let catalogue = catalogue();
    let terminal = residual("T904", &[1]);
    // 30 validated digits serve a 20-digit request with the decimal value.
    let value = run(
        Some(&catalogue),
        terminal.clone(),
        RustRedMasterMode::Numerical,
        6,
        20,
    )
    .unwrap();
    let [coefficient] = at_unit_scale(&value, 0, 0).try_into().unwrap();
    let AtomView::Num(number) = coefficient.as_view() else {
        panic!("a decimal coefficient: {value}");
    };
    let CoefficientView::Float(real, _) = number.get_coeff_view() else {
        panic!("a decimal coefficient: {value}");
    };
    assert!(
        (symbolica::domains::float::RealLike::to_f64(&real.to_float()) - 1.25).abs() < 1e-15,
        "{value}"
    );
    // vk5l WP4 (as the vacuum path): a 40-digit request is served with a
    // warning, at the 30 validated digits of the coefficient.
    let value = run(
        Some(&catalogue),
        terminal.clone(),
        RustRedMasterMode::Numerical,
        6,
        40,
    )
    .unwrap();
    assert!(!value.is_zero(), "{value}");
    // Without validated digits it is refused at any precision.
    let mut file = fixture();
    file["masters"][1]["coefficients"]["0"]
        .as_object_mut()
        .unwrap()
        .remove("significant_digits");
    let undigited = load(&file).unwrap();
    let error = catalogue_error(
        run(
            Some(&undigited),
            terminal.clone(),
            RustRedMasterMode::Numerical,
            6,
            5,
        )
        .unwrap_err(),
    );
    assert_names(
        &error.to_string(),
        &["syn-num", "without validated significant digits"],
    );
    // An embedded `2.5`100` with its 5 validated digits serves a 20-digit
    // request with a warning (vacuum path), and a 5-digit one silently.
    for digits in [20, 5] {
        assert!(
            run(
                Some(&catalogue),
                residual("T905", &[1]),
                RustRedMasterMode::Numerical,
                6,
                digits
            )
            .is_ok()
        );
    }
    // Symbolic mode never substitutes a decimal.
    let error = catalogue_error(
        run(
            Some(&catalogue),
            terminal,
            RustRedMasterMode::Symbolic,
            6,
            20,
        )
        .unwrap_err(),
    );
    assert!(
        matches!(error, PMasterCatalogueError::NumericalOnly { .. }),
        "{error}"
    );
    assert_names(
        &error.to_string(),
        &["syn-num coefficient epsilon^0", "master_mode=numerical"],
    );
    // A coefficient beyond the requested order needs no digits at all.
    assert!(
        run(
            Some(&catalogue),
            residual("T904", &[1]) * eps(),
            RustRedMasterMode::Symbolic,
            6,
            40
        )
        .unwrap()
        .is_zero()
    );
}

/// Exact constants of the catalogue stay exact symbols (`ζ`, `π`, `MZV`).
#[test]
fn exact_5l_constants_stay_symbolic() {
    let catalogue = catalogue();
    let value = run(
        Some(&catalogue),
        residual("T906", &[1]),
        RustRedMasterMode::Symbolic,
        6,
        40,
    )
    .unwrap();
    let [coefficient] = at_unit_scale(&value, 0, 0).try_into().unwrap();
    let expected = parse_constant("zeta(3)+pi^2-MZV(5,3)").unwrap();
    assert!((coefficient - expected).expand().is_zero(), "{value}");
}

/// Unmapped and additional terminals (decisions Q3, Q4): Retain returns them
/// as `PResidual`; Symbolic and Numerical refuse, naming each terminal, why
/// it has no value, the catalogue and the request.
#[test]
fn unmapped_residual_is_unevaluated_in_symbolic_and_numerical() {
    let catalogue = catalogue();
    let reduced = residual("T901", &[2]) + residual("T902", &[1]) + residual("T999", &[3, 1]);
    let retained = run(
        Some(&catalogue),
        reduced.clone(),
        RustRedMasterMode::Retain,
        6,
        40,
    )
    .unwrap();
    assert!(retained.contains_symbol(p_residual_symbol()));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let error =
            catalogue_error(run(Some(&catalogue), reduced.clone(), mode, 6, 40).unwrap_err());
        let PMasterCatalogueError::UnmappedTerminals { terminals, .. } = &error else {
            panic!("expected unmapped terminals: {error}");
        };
        assert_eq!(terminals.len(), 3);
        assert_names(
            &error.to_string(),
            &[
                "PResidual(T901,2): no master map row yet (WP3 row pending)",
                "PResidual(T902,1): additional terminal of a boxed node",
                "spoke dot",
                "PResidual(T999,3,1): not listed in the master map",
                &format!("master_mode={mode}"),
                "test-5L-v1",
                REQUEST,
            ],
        );
    }
}

/// Today (no catalogue, decision Q3): five-loop terminals are Retain-only.
#[test]
fn without_a_catalogue_five_loop_terminals_are_retain_only() {
    let terminal = residual(
        "T162",
        &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
    );
    let retained = run(None, terminal.clone(), RustRedMasterMode::Retain, 6, 40).unwrap();
    assert!(retained.contains_symbol(p_residual_symbol()));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let error = catalogue_error(run(None, terminal.clone(), mode, 6, 40).unwrap_err());
        assert!(
            matches!(error, PMasterCatalogueError::NoCatalogue { .. }),
            "{error}"
        );
        assert_names(
            &error.to_string(),
            &[
                "PResidual(T162,1,1,",
                &format!("master_mode={mode}"),
                "master_mode=retain",
                REQUEST,
            ],
        );
    }
}

/// ← `malformed_or_nonfinite_numerical_coefficients_report_the_master`, plus
/// the schema, version and depth checks of the loader.
#[test]
fn malformed_catalogues_are_refused_with_the_master_and_order() {
    let refused = |edit: &dyn Fn(&mut Value), parts: &[&str]| {
        let mut file = fixture();
        edit(&mut file);
        let error = load(&file).unwrap_err();
        assert!(
            matches!(error, PMasterCatalogueError::Invalid { .. }),
            "{error}"
        );
        assert_names(&error.to_string(), parts);
    };
    for literal in ["NaN", "inf", "not-a-number"] {
        refused(
            &|file| file["masters"][1]["coefficients"]["0"]["numerical"] = json!(literal),
            &["master syn-num coefficient epsilon^0"],
        );
    }
    refused(
        &|file| file["schema"] = json!("other.v1"),
        &["schema other.v1"],
    );
    refused(
        &|file| file["version"] = json!(2),
        &["catalogue version 2", "reads version 1"],
    );
    refused(
        &|file| file["normalisation"] = json!("msbar"),
        &["normalisation"],
    );
    refused(
        &|file| file["masters"][0]["coefficients"]["1"] = json!({"exact": "1"}),
        &["master syn-m1 coefficient epsilon^1", "unknown tail"],
    );
    refused(
        &|file| {
            file["masters"][0]["coefficients"]
                .as_object_mut()
                .unwrap()
                .remove("0");
            file["masters"][0]["coefficients"]["1"] = json!({"exact": "1"});
            file["masters"][0]["unknown_from"] = json!(2);
        },
        &["master syn-m1", "epsilon^1 follows the missing epsilon^0"],
    );
    refused(
        &|file| file["masters"][1]["known_through"] = json!(3),
        &["syn-num", "known_through 3"],
    );
    refused(
        &|file| file["masters"][0]["coefficients"]["0"]["exact"] = json!("pint_undetermined(1,2)"),
        &["syn-m1 coefficient epsilon^0", "unknown constant"],
    );
    refused(
        &|file| file["masters"][0]["coefficients"]["0"]["exact"] = json!("5/2*ep"),
        &["syn-m1 coefficient epsilon^0", "unknown constant"],
    );
    refused(
        &|file| file["masters"][0]["symbol"] = json!("no1"),
        &["vendored Forcer master"],
    );
    refused(
        &|file| file["masters"][0]["symbol"] = json!("ggppss_int(1)"),
        &["plain identifier"],
    );
    refused(
        &|file| file["masters"][1]["coefficients"]["0"] = json!({}),
        &[
            "syn-num coefficient epsilon^0",
            "neither an exact nor a numerical",
        ],
    );
    refused(
        &|file| file["terminals"][0]["combination"][0]["master"] = json!("syn-missing"),
        &["PResidual(T900,1,1,0)", "undeclared master syn-missing"],
    );
    refused(
        &|file| file["terminals"][0]["combination"][0]["coefficient"] = json!("x/(d-4)"),
        &["PResidual(T900,1,1,0)", "neither d nor a tag"],
    );
    refused(
        &|file| {
            let row = file["terminals"][0].clone();
            file["terminals"].as_array_mut().unwrap().push(row);
        },
        &["PResidual(T900,1,1,0) is listed twice"],
    );
    refused(
        &|file| file["masters"][0]["extra"] = json!(1),
        &["not a catalogue", "unknown field"],
    );
}

/// The Forcer table as a catalogue (real data through the same path): every
/// case of the frozen four-loop battery (54 cases and the three-loop boundary
/// case) finalizes at 11 terms to the frozen Forcer value, or fails with the
/// frozen `MissingMasterOrders` (26 `Oep` cases), exactly as on the Forcer
/// path.
#[test]
fn the_forcer_masters_through_the_catalogue_reproduce_the_four_loop_oep_battery() {
    crate::Vakint::initialize_vakint_symbols();
    let table = PMasterCatalogue::from_forcer_table().unwrap();
    assert_eq!(table.counts().0, 16);
    let battery =
        ForcerFixture::parse(include_str!("../../tests/fixtures/forcer/p4l_battery.json")).unwrap();
    let settings = finalize_settings(I7_TERMS);
    let (mut values, mut oep) = (0, 0);
    for case in battery.cases.iter().chain(&battery.boundary) {
        let frozen = case.frozen(vk_symbol!(&settings.epsilon_symbol)).unwrap();
        let result = finalize_with(
            &frozen,
            case.loops().unwrap(),
            &vk_parse!(&case.external).unwrap(),
            &settings,
            RustRedMasterMode::Symbolic,
            "Forcer",
            PMasterSource::new(Some(&table), "battery"),
        );
        let found = match result {
            Ok(value) => {
                values += 1;
                Finalized::Symbolic(case.unit_scale_coefficients(&value, I7_TERMS).unwrap())
            }
            Err(VakintError::MissingMasterOrders(missing)) => {
                oep += 1;
                assert!(missing.backend.contains("forcer-masters-prc"), "{missing}");
                Finalized::MissingMasterOrders(
                    missing
                        .masters
                        .iter()
                        .map(|master| super::super::forcer::fixtures::MissingOrder {
                            // `PMaster(lala) [lala, 4-loop, ...]` names the master first.
                            master: master.master.split(" [").next().unwrap().to_owned(),
                            first_missing_order: master.first_missing_order,
                            first_affected_result_order: master.first_affected_result_order,
                        })
                        .collect(),
                )
            }
            Err(other) => panic!("{}: {other}", case.name),
        };
        assert_eq!(
            Some(&found),
            case.eleven_terms.as_ref(),
            "{}: the catalogue path differs from the frozen Forcer value",
            case.name
        );
    }
    assert_eq!((values, oep), (29, 26));
}

/// The GGPPSS catalogue (`tools/p_integral_registry/ggppss_catalogue.py` on
/// vakint-external-data's `p_integral/masters/5L-ggppss/masters.jsonl`), when
/// `VAKINT_TEST_P_MASTER_CATALOGUE` names it: 281 masters load, the bare
/// `ggppss_int(23)` is `ε^-5` exactly through its depth, and one order beyond
/// its first undetermined order is refused by the `Oep` guard.
#[test]
fn the_ggppss_catalogue_loads_when_named() {
    let Some(path) = std::env::var_os("VAKINT_TEST_P_MASTER_CATALOGUE") else {
        eprintln!("VAKINT_TEST_P_MASTER_CATALOGUE unset: skipped");
        return;
    };
    let catalogue = PMasterCatalogue::load(std::path::Path::new(&path)).unwrap();
    assert_eq!(catalogue.counts().0, 281);
    let int23 = catalogue.master("ggppss_int_23").unwrap();
    let depth = int23.known_through;
    let master = p_master("ggppss_int_23");
    let value = run(
        Some(&catalogue),
        master.clone(),
        RustRedMasterMode::Symbolic,
        depth + 6,
        40,
    )
    .unwrap();
    let coefficients = at_unit_scale(&value, -5, depth);
    assert_eq!(coefficients[0], Atom::num(1));
    assert!(coefficients[1..].iter().all(Atom::is_zero), "{value}");
    // A master one order beyond its known depth (below ε^14) is refused.
    let (symbol, master) = catalogue
        .masters
        .iter()
        .filter(|(_, m)| m.known_through < 14)
        .max_by_key(|(_, m)| m.known_through)
        .unwrap();
    let missing = expect_missing_orders(
        run(
            Some(&catalogue),
            p_master(symbol),
            RustRedMasterMode::Symbolic,
            master.known_through + 7,
            40,
        )
        .unwrap_err(),
    );
    assert_eq!(
        missing.masters[0].first_missing_order,
        master.known_through + 1
    );
    assert!(missing.backend.contains("ggppss-5L-v1"), "{missing}");
}
