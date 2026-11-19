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
use symbolica::function;

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
        "schema": P_MASTER_CATALOGUE_SCHEMA_V1,
        "version": 1,
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
            spectator: None,
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
    let expected = parse_constant("zeta(3)+pi^2-MZV(5,3)", false).unwrap();
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
        &|file| file["version"] = json!(3),
        &["version 3", "reads vakint.p-master-catalogue.v2 version 2"],
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

/// Lane unknown-masters: a version-2 catalogue. `u_m1` = `1/ε + 5/2 +
/// UNKNOWN(u_m1,1)·ε + UNKNOWN(u_m1,2)·ε² + ALARM(ε³)` (file default
/// `unknown_through_default: 2`); `u_deep` known beyond the default (ε^0..ε^4,
/// ALARM from ε^5); `u_own` with its own `unknown_through: 1`; `u_x` a
/// placeholder at the proven bound -10 (all orders ε^-10..ε^2 UNKNOWN), carrying the unmapped terminal
/// T911 the way helper VH 2's `X[node:key][k]` atoms do.
pub(crate) fn fixture_v2() -> Value {
    json!({
        "schema": P_MASTER_CATALOGUE_SCHEMA,
        "version": P_MASTER_CATALOGUE_VERSION,
        "data_version": "test-5L-v2",
        "normalisation": P_MASTER_NORMALISATION,
        "unknown_through_default": 2,
        "masters": [
            {"id": "u-m1", "symbol": "u_m1", "loops": 5, "graph": "u-graph-1",
             "leading_order": -1, "known_through": 0, "alarm_from": 3,
             "coefficients": {"-1": {"exact": "1"}, "0": {"exact": "5/2"}}},
            {"id": "u-deep", "symbol": "u_deep", "loops": 5, "graph": "u-graph-2",
             "leading_order": 0, "coefficients": {"0": {"exact": "1"}, "1": {"exact": "2"},
             "2": {"exact": "3"}, "3": {"exact": "4"}, "4": {"exact": "5"}}},
            {"id": "u-own", "symbol": "u_own", "loops": 5, "graph": "u-graph-3",
             "leading_order": 0, "unknown_through": 1, "coefficients": {"0": {"exact": "7"}}},
            {"id": "u-x", "symbol": "u_x", "loops": 5, "graph": "u-graph-4",
             "leading_order": -10, "placeholder": true},
            {"id": "u-dep", "symbol": "u_dep", "loops": 5, "graph": "u-graph-5",
             "leading_order": 0, "coefficients": {"0": {"exact": "1"},
             "1": {"exact": "2+3*UNKNOWN(u_m1,1)"}}},
            {"id": "u-num", "symbol": "u_num", "loops": 5, "graph": "u-graph-6",
             "leading_order": 0, "unknown_through": 1,
             "coefficients": {"0": {"numerical": "1.25", "significant_digits": 30}}},
            {"id": "u-self", "symbol": "u_self", "loops": 5, "graph": "u-graph-7",
             "leading_order": 0, "coefficients": {"0": {"exact": "1"},
             "1": {"exact": "UNKNOWN(u_self,1)"}, "2": {"exact": "4-UNKNOWN(u_self,1)"}}}
        ],
        "terminals": [
            {"node": "T910", "key": [1], "status": "mapped",
             "combination": [{"master": "u-m1", "coefficient": "1"}]},
            {"node": "T911", "key": [1], "status": "mapped",
             "combination": [{"master": "u-x", "coefficient": "1"}]}
        ]
    })
}

fn run_allowing(
    catalogue: &PMasterCatalogue,
    reduced: Atom,
    mode: RustRedMasterMode,
    allow: bool,
) -> Result<Atom, VakintError> {
    crate::Vakint::initialize_vakint_symbols();
    let mut settings = settings(6, 40);
    settings.allow_unknown_master_coefficients = allow;
    finalize_with(
        &[ReducedPropagatorValue {
            scale_power: 0,
            reduced,
            spectator: None,
        }],
        5,
        &vk_parse!("p(1)").unwrap(),
        &settings,
        mode,
        "RustRed",
        PMasterSource::new(Some(catalogue), REQUEST),
    )
}

fn unknown(master: &str, order: i64) -> Atom {
    crate::truncation_guard::unknown_coefficient(master, order)
}

fn expect_unknowns(error: VakintError) -> crate::UnknownMasterCoefficients {
    match error {
        VakintError::UnknownMasterCoefficients(unknowns) => unknowns,
        other => panic!("expected unknown master coefficients, got: {other}"),
    }
}

/// `(master, order, first affected result order)` of a refusal.
fn unknown_entries(unknowns: &crate::UnknownMasterCoefficients) -> Vec<(String, i64, i64)> {
    unknowns
        .coefficients
        .iter()
        .map(|u| (u.master.clone(), u.order, u.first_affected_result_order))
        .collect()
}

#[test]
fn version_two_catalogues_declare_expected_unknown_orders() {
    let catalogue = load(&fixture_v2()).unwrap();
    let depths = |symbol: &str| {
        let master = catalogue.master(symbol).unwrap();
        (
            master.known_through,
            master.unknown_through,
            master.placeholder,
        )
    };
    assert_eq!(depths("u_m1"), (0, 2, false));
    assert_eq!(depths("u_deep"), (4, 4, false));
    assert_eq!(depths("u_own"), (0, 1, false));
    assert_eq!(depths("u_x"), (-11, 2, true));
    // The expansion: known orders, UNKNOWN(symbol, k) up to the expected
    // depth, then the ALARM tail.
    let source = PMasterSource::new(Some(&catalogue), REQUEST);
    let e = vk_symbol!("ε");
    let ep = Atom::var(e);
    let expected = ep.pow(Atom::num(-1))
        + Atom::num(5) / Atom::num(2)
        + unknown("u_m1", 1) * &ep
        + unknown("u_m1", 2) * ep.pow(Atom::num(2))
        + crate::truncation_guard::unknown_tail(e, 3, &p_master("u_m1"));
    assert_eq!(
        (source.expansion("u_m1", e).unwrap() - expected).expand(),
        Atom::Zero
    );
    // Version 1 files keep their meaning (no UNKNOWN orders).
    let v1 = catalogue_v1_m1();
    assert_eq!(
        (v1.known_through, v1.unknown_through),
        (0, 0),
        "a version 1 master has no expected unknown orders"
    );
}

fn catalogue_v1_m1() -> PCatalogueMaster {
    catalogue().master("syn_m1").unwrap().clone()
}

/// Refused (default): a typed error naming every (master, order) that entered
/// the result; allowed: the symbols stay, with their exact coefficients, in
/// Symbolic and Numerical mode.
#[test]
fn unknown_master_coefficients_are_refused_or_kept() {
    let catalogue = load(&fixture_v2()).unwrap();
    // ε^-2·u_m1 through ε^0: needs u_m1 through ε^2 (both UNKNOWN orders).
    let reduced = p_master("u_m1") * eps().pow(Atom::num(-2));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let error = run_allowing(&catalogue, reduced.clone(), mode, false).unwrap_err();
        let message = error.to_string();
        let unknowns = expect_unknowns(error);
        assert_eq!(
            unknown_entries(&unknowns),
            [("u-m1".into(), 1, -1), ("u-m1".into(), 2, 0)],
            "{message}"
        );
        assert_names(&message, &["u-m1", "allow_unknown_master_coefficients"]);
        let value = run_allowing(&catalogue, reduced.clone(), mode, true).unwrap();
        assert!(value.contains_symbol(crate::truncation_guard::unknown_coefficient_symbol()));
        if mode == RustRedMasterMode::Symbolic {
            assert_eq!(
                at_unit_scale(&value, -3, 0),
                [
                    Atom::num(1),
                    Atom::num(5) / Atom::num(2),
                    unknown("u_m1", 1),
                    unknown("u_m1", 2)
                ]
            );
        }
    }
    // Through the known orders only: no UNKNOWN, no error, even refused.
    let known = run_allowing(
        &catalogue,
        p_master("u_m1"),
        RustRedMasterMode::Symbolic,
        false,
    );
    assert_eq!(
        at_unit_scale(&known.unwrap(), -1, 0),
        [Atom::num(1), Atom::num(5) / Atom::num(2)]
    );
}

/// ALARM (beyond the expected depth) is the MissingMasterOrders error
/// whatever the option, also for a master known beyond the default depth and
/// for one with its own depth.
#[test]
fn alarm_orders_always_fail() {
    let catalogue = load(&fixture_v2()).unwrap();
    for (reduced, master, first) in [
        (p_master("u_m1") * eps().pow(Atom::num(-3)), "u-m1", 3),
        (p_master("u_deep") * eps().pow(Atom::num(-5)), "u-deep", 5),
        (p_master("u_own") * eps().pow(Atom::num(-2)), "u-own", 2),
    ] {
        for allow in [false, true] {
            let missing = expect_missing_orders(
                run_allowing(
                    &catalogue,
                    reduced.clone(),
                    RustRedMasterMode::Symbolic,
                    allow,
                )
                .unwrap_err(),
            );
            assert_eq!(missing.masters.len(), 1, "{master}");
            assert!(missing.masters[0].master.contains(master), "{missing:?}");
            assert_eq!(missing.masters[0].first_missing_order, first);
        }
    }
    // u_deep through its known ε^4, u_own through its UNKNOWN ε^1.
    let deep = run_allowing(
        &catalogue,
        p_master("u_deep") * eps().pow(Atom::num(-4)),
        RustRedMasterMode::Symbolic,
        false,
    )
    .unwrap();
    assert_eq!(at_unit_scale(&deep, -4, 0).last().unwrap(), &Atom::num(5));
    let own = run_allowing(
        &catalogue,
        p_master("u_own") * eps().pow(Atom::num(-1)),
        RustRedMasterMode::Symbolic,
        true,
    )
    .unwrap();
    assert_eq!(
        at_unit_scale(&own, -1, 0),
        [Atom::num(7), unknown("u_own", 1)]
    );
}

/// UNKNOWN terms cancel exactly when their exact coefficients cancel (the
/// laurent-depth cancellation pitfall): `ε^-3·(u_m1 - T910)` is zero, and
/// `ε^-2·u_m1·(ε·PG(0,0) - 1 + ε³)` is `ε·u_m1` (only known orders), so both
/// pass with the option off; a sum whose UNKNOWN parts do not cancel is
/// refused with exactly the surviving pairs.
#[test]
fn unknown_terms_cancel_only_exactly() {
    let catalogue = load(&fixture_v2()).unwrap();
    let master = p_master("u_m1");
    let pg = function!(pg_symbol(), Atom::num(0), Atom::num(0));
    let zero = (residual("T910", &[1]) - &master) * eps().pow(Atom::num(-3));
    let partial = (eps().pow(Atom::num(-2))
        * &master
        * (eps() * &pg - Atom::num(1) + eps().pow(Atom::num(3))))
    .expand();
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        assert!(
            run_allowing(&catalogue, zero.clone(), mode, false)
                .unwrap()
                .is_zero()
        );
        let value = run_allowing(&catalogue, partial.clone(), mode, false).unwrap();
        assert!(!value.contains_symbol(crate::truncation_guard::unknown_coefficient_symbol()));
    }
    let value = run_allowing(&catalogue, partial, RustRedMasterMode::Symbolic, false).unwrap();
    assert_eq!(at_unit_scale(&value, -1, 0), [Atom::Zero, Atom::num(1)]);
    // Not cancelling: ε^-2·u_m1·(ε·PG(0,0) + 1) = 2·ε^-2·u_m1.
    let survivor = (eps().pow(Atom::num(-2)) * &master * (eps() * &pg + Atom::num(1))).expand();
    let unknowns = expect_unknowns(
        run_allowing(&catalogue, survivor, RustRedMasterMode::Symbolic, false).unwrap_err(),
    );
    assert_eq!(
        unknown_entries(&unknowns),
        [("u-m1".into(), 1, -1), ("u-m1".into(), 2, 0)]
    );
}

/// A placeholder master (an unmapped terminal carried as its own master):
/// every order of its leading..expected range is UNKNOWN.
#[test]
fn placeholder_masters_are_all_unknown() {
    let catalogue = load(&fixture_v2()).unwrap();
    let reduced = residual("T911", &[1]);
    let unknowns = expect_unknowns(
        run_allowing(
            &catalogue,
            reduced.clone(),
            RustRedMasterMode::Symbolic,
            false,
        )
        .unwrap_err(),
    );
    assert_eq!(
        unknown_entries(&unknowns),
        (-10..=0)
            .map(|k| ("u-x".to_string(), k, k))
            .collect::<Vec<_>>()
    );
    let value = run_allowing(&catalogue, reduced, RustRedMasterMode::Symbolic, true).unwrap();
    assert_eq!(
        at_unit_scale(&value, -10, 0),
        (-10..=0).map(|k| unknown("u_x", k)).collect::<Vec<_>>()
    );
    // The loader refuses a placeholder above the proven bound -2 x loops and a
    // mapped master whose declared leading coefficient is zero.
    let mut shallow = fixture_v2();
    shallow["masters"][3]["leading_order"] = json!(-7);
    assert_names(
        &load(&shallow).unwrap_err().to_string(),
        &["u-x", "placeholder leading_order -7", "-10"],
    );
    let mut zero_leading = fixture_v2();
    zero_leading["masters"][1]["coefficients"]["0"] = json!({"exact": "0"});
    assert_names(
        &load(&zero_leading).unwrap_err().to_string(),
        &["u-deep", "declared leading order epsilon^0 is zero"],
    );
    // The zero test is exact over a common denominator in the constants
    // (reviewer's witness: 1/a - 1/b - 1/(ab) = 0 with a = 1+zeta(3),
    // b = 2+zeta(3)), never numerical.
    let mut hidden_zero = fixture_v2();
    hidden_zero["masters"][1]["coefficients"]["0"] = json!({
        "exact": "1/(1+zeta(3))-1/(2+zeta(3))-1/((1+zeta(3))*(2+zeta(3)))"
    });
    assert_names(
        &load(&hidden_zero).unwrap_err().to_string(),
        &["u-deep", "declared leading order epsilon^0 is zero"],
    );
    // A dependence may not refer below a master's leading order.
    let mut below = fixture_v2();
    below["masters"][4]["coefficients"]["1"] = json!({"exact": "UNKNOWN(u_x,-11)"});
    assert_names(
        &load(&below).unwrap_err().to_string(),
        &["u-dep", "UNKNOWN(u_x, -11)"],
    );
}

#[test]
fn malformed_version_two_catalogues_are_refused() {
    let cases: [(&dyn Fn(&mut Value), &[&str]); 5] = [
        (
            &|file| file["masters"][0]["unknown_from"] = json!(1),
            &["u-m1", "unknown_from is a version 1 field"],
        ),
        (
            &|file| file["masters"][1]["unknown_through"] = json!(3),
            &["u-deep", "unknown_through 3 below its last known order"],
        ),
        (
            &|file| file["masters"][0]["alarm_from"] = json!(4),
            &["u-m1", "alarm_from 4"],
        ),
        (
            &|file| file["masters"][3]["coefficients"] = json!({"-2": {"exact": "1"}}),
            &["u-x", "a placeholder master has no coefficients"],
        ),
        (
            &|file| {
                file["schema"] = json!(P_MASTER_CATALOGUE_SCHEMA_V1);
                file["version"] = json!(1);
            },
            &["unknown_through_default needs schema"],
        ),
    ];
    for (edit, parts) in cases {
        let mut file = fixture_v2();
        edit(&mut file);
        let message = load(&file).unwrap_err().to_string();
        assert_names(&message, parts);
    }
}

/// Lane unknown-masters, physical example (environment driven, skipped
/// without it): a RustRed `p5l_library reduce-file` output of one five-loop
/// key (`VAKINT_TEST_UNKNOWN_KEY_FILE`: lines `Unit GammaTag([(m, e), ..]) :
/// c(d)` and `PResidual(node, [k, ..]) GammaTag(..) : c(d)`) finalized in
/// G-scheme Symbolic mode through ε^0 with a version 1 catalogue
/// (`VAKINT_TEST_UNKNOWN_CATALOGUE`) raised to version 2 with
/// `unknown_through_default` 12. Refused: prints the (master, order) pairs;
/// allowed: prints the ε^0 coefficient of every UNKNOWN symbol.
/// `VAKINT_TEST_UNKNOWN_EXPECT="ggppss_int_53:7=-7/16;ggppss_int_74:6=7/288"`
/// asserts those ε^0 coefficients and that no other UNKNOWN enters ε^0.
#[test]
fn unknown_masters_physical_key_file() {
    let (Ok(key_file), Ok(catalogue_file)) = (
        std::env::var("VAKINT_TEST_UNKNOWN_KEY_FILE"),
        std::env::var("VAKINT_TEST_UNKNOWN_CATALOGUE"),
    ) else {
        eprintln!("unknown_masters_physical_key_file: environment unset, skipped");
        return;
    };
    crate::Vakint::initialize_vakint_symbols();
    let mut file: Value =
        serde_json::from_str(&std::fs::read_to_string(&catalogue_file).unwrap()).unwrap();
    file["schema"] = json!(P_MASTER_CATALOGUE_SCHEMA);
    file["version"] = json!(P_MASTER_CATALOGUE_VERSION);
    file["unknown_through_default"] = json!(12);
    for master in file["masters"].as_array_mut().unwrap() {
        master.as_object_mut().unwrap().remove("unknown_from");
    }
    // Optional: GGPPSS's printed dependences (helper VH 2's dependences.v1.tsv:
    // `master order kind ... combination`, DEPENDENT rows `(a)*Int[i][j] + ..`)
    // and extra relations `id:order=expression;...` (expressions in
    // UNKNOWN(symbol,k)). Helper VH 2's notation is translated: zN -> zeta(N),
    // m35 / m37 / m335 (GGPPSS MZV[3,5] ..., increasing) -> MZV(5,3) / MZV(7,3)
    // / MZV(5,3,3) (the catalogue's convention), Int[i][j] ->
    // UNKNOWN(ggppss_int_i,j).
    let helper_text = |text: &str| {
        let text = regex::Regex::new(r"Int\[(\d+)\]\[(-?\d+)\]")
            .unwrap()
            .replace_all(text, "UNKNOWN(ggppss_int_$1,$2)")
            .into_owned();
        let text = regex::Regex::new(r"\bm335\b")
            .unwrap()
            .replace_all(&text, "MZV(5,3,3)")
            .into_owned();
        let text = regex::Regex::new(r"\bm35\b")
            .unwrap()
            .replace_all(&text, "MZV(5,3)")
            .into_owned();
        let text = regex::Regex::new(r"\bm37\b")
            .unwrap()
            .replace_all(&text, "MZV(7,3)")
            .into_owned();
        // Mathematica forms of the dependence file: Zeta[n], MZV[a,b,..]
        // (GGPPSS increasing indices -> the catalogue's decreasing order).
        let text = regex::Regex::new(r"\bZeta\[(\d+)\]")
            .unwrap()
            .replace_all(&text, "zeta($1)")
            .into_owned();
        let text = regex::Regex::new(r"\bMZV\[([\d,\s]+)\]")
            .unwrap()
            .replace_all(&text, |c: &regex::Captures| {
                let mut indices: Vec<&str> = c[1].split(',').map(str::trim).collect();
                indices.reverse();
                format!("MZV({})", indices.join(","))
            })
            .into_owned();
        regex::Regex::new(r"\bz(\d+)\b")
            .unwrap()
            .replace_all(&text, "zeta($1)")
            .into_owned()
    };
    let mut relations: Vec<(String, i64, String)> = Vec::new();
    if let Ok(path) = std::env::var("VAKINT_TEST_UNKNOWN_DEPENDENCES") {
        // dependences.v1.json: DEPENDENT orders = known_part + sum symbols.
        let table: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        for entry in table["dependences"].as_array().unwrap() {
            if entry["kind"] != json!("DEPENDENT") {
                continue;
            }
            let mut expression = format!("({})", entry["known_part"].as_str().unwrap());
            for (symbol, coefficient) in entry["symbols"].as_object().unwrap() {
                expression += &format!("+({})*{symbol}", coefficient.as_str().unwrap());
            }
            relations.push((
                entry["master"].as_str().unwrap().into(),
                entry["order"].as_i64().unwrap(),
                helper_text(&expression),
            ));
        }
    }
    if let Ok(extra) = std::env::var("VAKINT_TEST_UNKNOWN_EXTRA") {
        for item in extra.split(';').filter(|item| !item.is_empty()) {
            let (head, expression) = item.split_once('=').unwrap();
            let (id, order) = head.rsplit_once(':').unwrap();
            relations.push((id.into(), order.parse().unwrap(), helper_text(expression)));
        }
    }
    // The catalogue admits references to FREE coordinates only: resolve
    // relations that refer to other relations (e.g. the 53 / 74 identity
    // inside a printed dependence) by substitution, to a fixed point.
    let symbol_of = |id: &str| {
        file["masters"]
            .as_array()
            .unwrap()
            .iter()
            .find(|m| m["id"] == json!(id))
            .map(|m| m["symbol"].as_str().unwrap().to_owned())
            .unwrap_or_else(|| panic!("no master {id}"))
    };
    let rules: Vec<(Atom, Atom)> = relations
        .iter()
        .map(|(id, order, expression)| {
            (
                unknown(&symbol_of(id), *order),
                vk_parse!(expression).unwrap(),
            )
        })
        .collect();
    let relations: Vec<(String, i64, String)> = relations
        .iter()
        .zip(&rules)
        .map(|((id, order, _), (_, expression))| {
            let mut value = expression.clone();
            for _ in 0..16 {
                let mut next = value.clone();
                for (symbol, replacement) in &rules {
                    next = next
                        .replace(symbol.to_pattern())
                        .with(replacement.to_pattern());
                }
                let next = next.expand();
                if next == value {
                    break;
                }
                value = next;
            }
            (
                id.clone(),
                *order,
                // The printer spells Symbolica's ζ and π as Zeta and Pi; the
                // catalogue reads zeta and pi.
                regex::Regex::new(r"\bPi\b")
                    .unwrap()
                    .replace_all(
                        &crate::truncation_guard::plain_name(&value).replace("Zeta(", "zeta("),
                        "pi",
                    )
                    .into_owned(),
            )
        })
        .collect();
    for (id, order, expression) in &relations {
        let master = file["masters"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|m| m["id"] == json!(id))
            .unwrap_or_else(|| panic!("no master {id}"));
        let symbol = master["symbol"].as_str().unwrap().to_owned();
        master.as_object_mut().unwrap().remove("known_through");
        let coefficients = master["coefficients"].as_object_mut().unwrap();
        let last = coefficients
            .keys()
            .map(|k| k.parse::<i64>().unwrap())
            .max()
            .unwrap();
        if *order <= last {
            eprintln!("dependence {id} order {order}: already published (through {last}), skipped");
            continue;
        }
        for free in last + 1..*order {
            coefficients.insert(
                free.to_string(),
                json!({"exact": format!("UNKNOWN({symbol},{free})")}),
            );
        }
        coefficients.insert(order.to_string(), json!({"exact": expression}));
    }
    eprintln!(
        "unknown_masters_physical_key_file: {} dependent orders loaded",
        relations.len()
    );
    let catalogue = load(&file).unwrap();
    let d = vk_symbol!("d");
    let dimension = Atom::num(4) - Atom::num(2) * eps();
    let mut reduced = Atom::Zero;
    let text = std::fs::read_to_string(&key_file).unwrap();
    let tag_pattern = regex::Regex::new(r"\((-?\d+), (-?\d+)\)").unwrap();
    let mut lines = 0;
    for line in text.lines() {
        let (head, coefficient) = match line.split_once(" : ") {
            Some(parts) if line.starts_with("Unit") || line.starts_with("PResidual(") => parts,
            _ => continue,
        };
        let (terminal, tags) = head.split_once(" GammaTag(").unwrap();
        let mut term = vk_parse!(coefficient)
            .unwrap()
            .replace(Atom::var(d).to_pattern())
            .with(dimension.to_pattern());
        for tag in tag_pattern.captures_iter(tags) {
            let m: i64 = tag[1].parse().unwrap();
            let e: i64 = tag[2].parse().unwrap();
            term *= p_gamma_symbol().call(&[Atom::num(m)][..]).pow(Atom::num(e));
        }
        if let Some(rest) = terminal.strip_prefix("PResidual(") {
            let (node, key) = rest.split_once(", [").unwrap();
            let key: Vec<i64> = key
                .trim_end_matches("])")
                .split(", ")
                .map(|k| k.parse().unwrap())
                .collect();
            term *= residual(node, &key);
        }
        reduced += term;
        lines += 1;
    }
    eprintln!("unknown_masters_physical_key_file: {lines} reduced terms from {key_file}");
    // RustRed's reduce-file prints I_M(K) (the terminal rows of the catalogue
    // likewise carry their (εG(1,1))^L): the G-scheme value of the key is
    // I_G(K) = (-1)^|K| (εG(1,1))^-L I_M(K) with
    // εG(1,1) = Γ(1+ε)Γ(1-ε)²/(Γ(1-2ε)(d-3)). VAKINT_TEST_UNKNOWN_LOOPS = L and
    // VAKINT_TEST_UNKNOWN_SIGN = (-1)^|K| apply it (without them the harness
    // finalises I_M(K), whose ε-expansion carries exp(-Lγε)).
    if let Ok(loops) = std::env::var("VAKINT_TEST_UNKNOWN_LOOPS") {
        let loops: i64 = loops.parse().unwrap();
        let tag = |m: i64, e: i64| p_gamma_symbol().call(&[Atom::num(m)][..]).pow(Atom::num(e));
        reduced = reduced
            * tag(1, -loops)
            * tag(-1, -2 * loops)
            * tag(-2, loops)
            * (Atom::num(1) - Atom::num(2) * eps()).pow(Atom::num(loops));
    }
    if let Ok(sign) = std::env::var("VAKINT_TEST_UNKNOWN_SIGN") {
        reduced = reduced * Atom::num(sign.parse::<i64>().unwrap());
    }
    let reduced = reduced.expand();
    let started = std::time::Instant::now();
    let refused = run_allowing(
        &catalogue,
        reduced.clone(),
        RustRedMasterMode::Symbolic,
        false,
    );
    match &refused {
        Err(VakintError::UnknownMasterCoefficients(unknowns)) => {
            eprintln!(
                "REFUSED in {:.1} s: {unknowns}",
                started.elapsed().as_secs_f64()
            )
        }
        other => eprintln!("NOT REFUSED BY UNKNOWNS: {other:?}"),
    }
    let started = std::time::Instant::now();
    let value = run_allowing(&catalogue, reduced, RustRedMasterMode::Symbolic, true).unwrap();
    eprintln!("ALLOWED in {:.1} s", started.elapsed().as_secs_f64());
    let finite = at_unit_scale(&value, -12, 0).pop().unwrap();
    let unknowns: Vec<Atom> = {
        let mut found = std::collections::BTreeSet::new();
        finite.visitor(&mut |part| {
            if let symbolica::atom::AtomView::Fun(f) = part
                && f.get_symbol() == crate::truncation_guard::unknown_coefficient_symbol()
            {
                found.insert(part.to_owned());
                return false;
            }
            true
        });
        found.into_iter().collect()
    };
    let mut coefficients = BTreeMap::new();
    for symbol in &unknowns {
        let coefficient = finite.coefficient(symbol.as_view()).expand();
        eprintln!("eps^0 coefficient of {symbol}: {coefficient}");
        coefficients.insert(crate::truncation_guard::plain_name(symbol), coefficient);
    }
    let known_finite = finite
        .replace(
            crate::truncation_guard::unknown_coefficient_symbol()
                .call(
                    &[
                        Atom::var(crate::symbols::S.x_),
                        Atom::var(crate::symbols::S.y_),
                    ][..],
                )
                .to_pattern(),
        )
        .with(Atom::Zero.to_pattern())
        .expand();
    eprintln!("eps^0 known part: {known_finite}");
    eprintln!(
        "eps^0 contains EulerGamma: {}",
        known_finite.contains_symbol(symbolica::transcendental::euler_gamma())
    );
    if let Ok(expect) = std::env::var("VAKINT_TEST_UNKNOWN_EXPECT_FINITE") {
        let expected = parse_constant(&helper_text(&expect), false).unwrap();
        assert_eq!(
            (&known_finite - &expected).expand(),
            Atom::Zero,
            "eps^0 known part {known_finite} != {expected}"
        );
        eprintln!("eps^0 known part EQUALS the expected value");
    }
    if let Ok(expect) = std::env::var("VAKINT_TEST_UNKNOWN_EXPECT") {
        let mut expected = BTreeMap::new();
        for item in expect.split(';') {
            let (symbol, value) = item.split_once('=').unwrap();
            let (master, order) = symbol.split_once(':').unwrap();
            expected.insert(
                crate::truncation_guard::plain_name(&unknown(master, order.parse().unwrap())),
                vk_parse!(value).unwrap(),
            );
        }
        assert_eq!(
            coefficients.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>()
        );
        for (name, value) in &expected {
            assert_eq!((&coefficients[name] - value).expand(), Atom::Zero, "{name}");
        }
        let refused = refused.unwrap_err().to_string();
        for name in expected.keys() {
            let id = name.split(['(', ',']).nth(1).unwrap();
            eprintln!("refusal names {id}: {}", refused.contains(id));
        }
    }
}

/// A published or derived dependence (an expected order written in terms of
/// free UNKNOWN coordinates and known numbers, as GGPPSS prints for some deeper
/// orders, or an identity between masters): `u_dep` at ε^1 is
/// `2 + 3·UNKNOWN(u_m1,1)`; `u_self` writes its free ε^1 explicitly and its
/// ε^2 in terms of it. The result then carries only free coordinates.
#[test]
fn dependent_orders_carry_only_free_coordinates() {
    let catalogue = load(&fixture_v2()).unwrap();
    assert_eq!(catalogue.master("u_dep").unwrap().known_through, 1);
    let value = run_allowing(
        &catalogue,
        p_master("u_dep") * eps().pow(Atom::num(-1)),
        RustRedMasterMode::Symbolic,
        true,
    )
    .unwrap();
    assert_eq!(
        at_unit_scale(&value, -1, 0),
        [
            Atom::num(1),
            Atom::num(2) + Atom::num(3) * unknown("u_m1", 1)
        ]
    );
    let refused = expect_unknowns(
        run_allowing(
            &catalogue,
            p_master("u_dep") * eps().pow(Atom::num(-1)),
            RustRedMasterMode::Symbolic,
            false,
        )
        .unwrap_err(),
    );
    assert_eq!(unknown_entries(&refused), [("u-m1".into(), 1, 0)]);
    // u_self: ε^-2·u_self through ε^0 = 1/ε² + U/ε + 4 - U: the free U.
    let value = run_allowing(
        &catalogue,
        p_master("u_self") * eps().pow(Atom::num(-2)),
        RustRedMasterMode::Symbolic,
        true,
    )
    .unwrap();
    assert_eq!(
        at_unit_scale(&value, -2, 0),
        [
            Atom::num(1),
            unknown("u_self", 1),
            Atom::num(4) - unknown("u_self", 1)
        ]
    );
    // References must be free coordinates of declared masters.
    let cases: [(&dyn Fn(&mut Value), &[&str]); 4] = [
        (
            &|file| file["masters"][4]["coefficients"]["1"] = json!({"exact": "UNKNOWN(u_m1,0)"}),
            &["u-dep", "UNKNOWN(u_m1, 0)", "not a free expected order"],
        ),
        (
            &|file| file["masters"][4]["coefficients"]["1"] = json!({"exact": "UNKNOWN(u_m1,3)"}),
            &["u-dep", "UNKNOWN(u_m1, 3)"],
        ),
        (
            &|file| file["masters"][2]["coefficients"]["1"] = json!({"exact": "UNKNOWN(u_dep,1)"}),
            &["u-own", "UNKNOWN(u_dep, 1)"],
        ),
        (
            &|file| file["masters"][4]["coefficients"]["1"] = json!({"exact": "UNKNOWN(nobody,1)"}),
            &["u-dep", "UNKNOWN(nobody, 1)"],
        ),
    ];
    for (edit, parts) in cases {
        let mut file = fixture_v2();
        edit(&mut file);
        let message = load(&file).unwrap_err().to_string();
        assert_names(&message, parts);
    }
    // Version 1 files may not write UNKNOWN.
    let mut v1 = fixture();
    v1["masters"][0]["coefficients"]["0"] = json!({"exact": "UNKNOWN(syn_m1,0)"});
    assert_names(&load(&v1).unwrap_err().to_string(), &["unknown constant"]);
}

/// Products keep products of unknowns: the finite part of (M - a)²/ε² with
/// a = M's known part is UNKNOWN(m,1)²; with an independent second master it
/// is the product of their two symbols; a numerical-only coefficient times two
/// distinct UNKNOWNs keeps both symbols (Numerical mode, refusal decided on the
/// exact dependence before rounding).
#[test]
fn products_of_unknowns_and_numerical_multipliers() {
    let catalogue = load(&fixture_v2()).unwrap();
    let ep = eps();
    let shifted = p_master("u_m1") - ep.pow(Atom::num(-1)) - Atom::num(5) / Atom::num(2);
    let square = (shifted.pow(Atom::num(2)) * ep.pow(Atom::num(-2))).expand();
    let value = run_allowing(
        &catalogue,
        square.clone(),
        RustRedMasterMode::Symbolic,
        true,
    )
    .unwrap();
    assert_eq!(
        at_unit_scale(&value, 0, 0),
        [unknown("u_m1", 1).pow(Atom::num(2))]
    );
    let refused = expect_unknowns(
        run_allowing(&catalogue, square, RustRedMasterMode::Symbolic, false).unwrap_err(),
    );
    assert_eq!(unknown_entries(&refused), [("u-m1".into(), 1, 0)]);
    let pair = (&shifted * (p_master("u_own") - Atom::num(7)) * ep.pow(Atom::num(-2))).expand();
    let value = run_allowing(&catalogue, pair.clone(), RustRedMasterMode::Symbolic, true).unwrap();
    assert_eq!(
        at_unit_scale(&value, 0, 0),
        [unknown("u_m1", 1) * unknown("u_own", 1)]
    );
    let numerical = (pair * p_master("u_num")).expand();
    let refused = expect_unknowns(
        run_allowing(
            &catalogue,
            numerical.clone(),
            RustRedMasterMode::Numerical,
            false,
        )
        .unwrap_err(),
    );
    assert_eq!(
        unknown_entries(&refused),
        [("u-m1".into(), 1, 0), ("u-own".into(), 1, 0)]
    );
    let value = run_allowing(&catalogue, numerical, RustRedMasterMode::Numerical, true).unwrap();
    let product = unknown("u_m1", 1) * unknown("u_own", 1);
    let multiplier = value.expand().coefficient(product.as_view());
    let difference = (&multiplier - Atom::num(5) / Atom::num(4)).expand();
    let close = match difference.as_view() {
        symbolica::atom::AtomView::Num(_) => {
            let x: f64 = difference.to_string().parse().unwrap_or(1.0);
            x.abs() < 1e-25
        }
        _ => difference.is_zero(),
    };
    assert!(close, "{value}: multiplier {multiplier}");
}

/// ALARM wins over UNKNOWN: a result with both is MissingMasterOrders in
/// both settings.
#[test]
fn mixed_unknown_and_alarm_is_an_alarm() {
    let catalogue = load(&fixture_v2()).unwrap();
    let reduced = p_master("u_m1") * (eps().pow(Atom::num(-2)) + eps().pow(Atom::num(-3)));
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        for allow in [false, true] {
            let missing = expect_missing_orders(
                run_allowing(&catalogue, reduced.clone(), mode, allow).unwrap_err(),
            );
            assert_eq!(missing.masters[0].first_missing_order, 3);
            assert!(missing.to_string().contains("ALARM"));
        }
    }
}
