//! Public-engine acceptance of the independently certified four-loop periods.
//! This test deliberately loads the complete local pack, including its graph
//! and IBP metadata, and never contacts a producer workspace or the network.
//!
//! The data-dependent expectations below (`D2_DIGITS`, `D2_DEPENDENT_CELLS`, the
//! 520-digit cap of the `P19`/`P73` ε¹ rows and the 739 published formulas) are
//! those of `vakint-external-data` `main` at `f51bc5a2` (the closed 4L formulas
//! with 520-digit D = 2 periods). Earlier packs such as `ca49f328` cap the D = 2
//! periods at 215 digits and fail these assertions; update them together with
//! the data.

use std::{fs, path::PathBuf};

use serde_json::Value;
use symbolica::domains::{float::Float, rational::Rational};

use crate::{
    AvailablePrecision, CoefficientPrecision, ConstantSymbol, ErrorQuality, MasterSymbol,
    PrecisionWarning, Vakint, precision::decimal_bits, utils::vakint_macros::vk_parse,
};

const CERTIFIED_DIGITS: u32 = 20_000;
/// Validated digits of the independently computed D = 2 periods and their jets.
const D2_DIGITS: u32 = 520;
/// D2-dependent coefficient cells of the takeover closures (`PR7D2`, `PR10D2`,
/// `PR15D2` and their jets): the original 15 plus the 43 closed formulas' cells.
const D2_DEPENDENT_CELLS: usize = 62;

/// Parse a finite decimal as an exact rational. In particular, this comparison
/// never rounds a 20,008-digit reference or its tiny error radius to f64.
fn exact_decimal(value: &str) -> Rational {
    let (mantissa, exponent) = value
        .split_once(['e', 'E'])
        .map_or((value, 0_i64), |(mantissa, exponent)| {
            (mantissa, exponent.parse().unwrap())
        });
    let decimals = mantissa
        .split_once('.')
        .map_or(0, |(_, fraction)| fraction.len() as i64);
    let integer = mantissa.replace('.', "");
    let expression = vk_parse!(&format!("({integer})*10^({})", exponent - decimals)).unwrap();
    Rational::try_from(expression.as_view()).unwrap()
}

fn check_value(
    label: &str,
    value: &Float,
    precision: &CoefficientPrecision,
    reference: &str,
    certified_error: &str,
    direct_leaf: bool,
) {
    assert_eq!(precision.real().quality(), ErrorQuality::Bounded, "{label}");
    assert!(precision.imaginary().exact_zero(), "{label}");
    assert!(!precision.real().exact_zero(), "{label}");
    assert!(
        value.prec() <= decimal_bits(CERTIFIED_DIGITS),
        "{label}: storage precision exceeds the declared source cap"
    );
    // Available source digits and achieved output digits are distinct. The
    // conservative decimal ULP floor, binary conversion and checked final
    // rounding can reduce the latter; the returned Float must carry precisely
    // the accuracy established by its enclosure, without an invented floor.
    assert_eq!(precision.real().significant_bits(), Some(value.prec()));
    let achieved = precision.real().significant_digits().unwrap();
    assert!(
        achieved <= CERTIFIED_DIGITS,
        "{label}: optimistic precision"
    );
    // One decimal ULP can cost one relative decimal digit; the documented
    // two-bit final-rounding reserve and downward binary/decimal conversion
    // cost at most one further digit for these directly returned real leaves.
    // This only bounds gratuitous conversion loss. Scientific accuracy is
    // still checked below against the complete exact certified interval.
    if direct_leaf {
        assert!(
            achieved >= CERTIFIED_DIGITS - 2,
            "{label}: unexplained conversion loss"
        );
    }

    let reference = exact_decimal(reference);
    let certified_error = exact_decimal(certified_error);
    assert!(certified_error > 0);
    let difference = value.to_rational() - reference;
    let distance = if difference < 0 {
        -difference
    } else {
        difference
    };
    // The returned interval must enclose the complete authenticated interval,
    // including both the source error and any center shift from final rounding.
    // There is no empirical tolerance or comparison at reduced precision.
    assert!(
        distance + certified_error <= precision.real().absolute_error().to_rational(),
        "{label}: returned uncertainty fails to enclose the installed certificate"
    );
    eprintln!(
        "{label}: achieved {achieved} digits, {} stored bits",
        value.prec()
    );
}

fn check_rounding_warning(
    warnings: &[PrecisionWarning],
    precision: &CoefficientPrecision,
    requested: u32,
    version: &str,
    period: Option<&str>,
) {
    let achieved = precision.real().significant_digits().unwrap();
    if achieved < requested {
        let warning = warnings
            .iter()
            .find(|warning| {
                warning.requested_decimal_precision == requested
                    && warning.achieved_decimal_precision == Some(achieved)
                    && warning.message.contains("rounding")
                    && warning.sources.iter().any(|source| {
                        source.master_id.is_some()
                            && source.epsilon_power == Some(0)
                            && source.available_decimal_precision == Some(CERTIFIED_DIGITS)
                            && period.is_none_or(|period| source.identity.contains(period))
                    })
            })
            .expect(
                "measured rounding loss must report achieved digits and its source coefficient",
            );
        assert!(warning.data_version.split(", ").any(|item| item == version));
    }
}

fn check_cap_warning(
    warnings: &[PrecisionWarning],
    version: &str,
    period: Option<&str>,
    requested_master: Option<&str>,
) {
    let warning = warnings
        .iter()
        .find(|warning| {
            warning.requested_decimal_precision == CERTIFIED_DIGITS + 1
                && warning.achieved_decimal_precision.is_none()
                && warning.sources.iter().any(|source| {
                    source.available_decimal_precision == Some(CERTIFIED_DIGITS)
                        && source.epsilon_power == Some(0)
                        && source.master_id.is_some()
                        && period.is_none_or(|period| source.identity.contains(period))
                })
        })
        .expect("exceeding the validated cap must identify the limiting period/master coefficient");
    assert!(warning.data_version.split(", ").any(|item| item == version));
    let displayed = warning.to_string();
    for expected in [version, "20001", "20000", "ε^0"] {
        assert!(
            displayed.contains(expected),
            "warning omits {expected}: {displayed}"
        );
    }
    if let Some(period) = period {
        assert!(displayed.contains(period), "{displayed}");
    }
    if let Some(master_id) = requested_master {
        assert!(
            warnings
                .iter()
                .flat_map(|warning| &warning.sources)
                .any(|source| {
                    source.master_id.as_deref() == Some(master_id)
                        && source.epsilon_power == Some(0)
                        && source.available_decimal_precision == Some(CERTIFIED_DIGITS)
                }),
            "warnings must identify the requested master coefficient as well as any equivalent payload source"
        );
    }
}

#[test]
#[ignore = "loads the certified 20,000-digit external 4L periods; set VAKINT_TEST_EXTERNAL_DATA"]
fn public_engine_certified_periods_at_20000_digits() {
    let root = std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
        .map(PathBuf::from)
        .expect("set VAKINT_TEST_EXTERNAL_DATA to the certified external-data release checkout");
    let manifest: Value =
        serde_json::from_slice(&fs::read(root.join("4L/manifest.json")).unwrap()).unwrap();
    let version = manifest["version"].as_str().unwrap();
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(&root, 4).unwrap();

    // Primitive moment values have the direct-leaf rounding budget above.
    // Derived PR4 formulas additionally track their actual conditioning;
    // their complete certified intervals and warnings remain checked below.
    if let Some(binding) = manifest.get("constant_bindings") {
        let path = binding["path"].as_str().unwrap();
        let entries: Vec<Value> =
            serde_json::from_slice(&fs::read(root.join("4L").join(path)).unwrap()).unwrap();
        for name in [
            "K5Moment(1)",
            "K5Moment(3)",
            "K5LogMoment(1)",
            "K5LogMoment(3)",
        ] {
            let expected: ConstantSymbol = name.parse().unwrap();
            let entry = entries
                .iter()
                .find(|entry| {
                    entry["symbol"]
                        .as_str()
                        .unwrap()
                        .parse::<ConstantSymbol>()
                        .unwrap()
                        == expected
                })
                .expect("the four certified primitive moments must remain installed");
            let symbol: ConstantSymbol = entry["symbol"].as_str().unwrap().parse().unwrap();
            let constant = engine.get_constant(&symbol, CERTIFIED_DIGITS).unwrap();
            assert_eq!(
                constant.available_precision(),
                AvailablePrecision::FixedPrecision(CERTIFIED_DIGITS)
            );
            check_value(
                &symbol.to_string(),
                constant.numerical_value(),
                constant.precision(),
                entry["value"].as_str().unwrap(),
                entry["absolute_error"].as_str().unwrap(),
                true,
            );
        }
    }

    for (period, name) in [
        ("PR4fin", "P6"),
        ("PR4dfin", "P11"),
        ("PR11dep0", "P17"),
        ("PR11dep0", "P52"),
        ("PR11dep0", "P74"),
    ] {
        let metadata = manifest["masters"]
            .as_array()
            .unwrap()
            .iter()
            .find(|master| master["name"].as_str() == Some(name))
            .unwrap();
        let master_id = metadata["master_id"].as_str().unwrap();
        let coefficient = metadata["coefficients"]
            .as_array()
            .unwrap()
            .iter()
            .find(|coefficient| coefficient["epsilon_power"].as_i64() == Some(0))
            .unwrap();
        assert_eq!(
            coefficient["significant_digits"].as_u64(),
            Some(u64::from(CERTIFIED_DIGITS))
        );
        assert_eq!(coefficient["error_quality"].as_str(), Some("bounded"));
        assert!(
            coefficient["stored_significant_digits"].as_u64().unwrap()
                > u64::from(CERTIFIED_DIGITS)
        );
        let certified_error = coefficient["absolute_error"].as_str().unwrap();
        let numerical_path = metadata["files"]["numerical.json"]["path"]
            .as_str()
            .unwrap();
        let rows: Vec<(i64, String)> =
            serde_json::from_slice(&fs::read(root.join("4L").join(numerical_path)).unwrap())
                .unwrap();
        let reference = &rows.iter().find(|(power, _)| *power == 0).unwrap().1;
        let constant_symbol: ConstantSymbol = period.parse().unwrap();
        // The periods have equivalent aliases in several families. The
        // full declared ID selects the intended graph without name ambiguity.
        let master_symbol: MasterSymbol = master_id.parse().unwrap();

        for requested in [CERTIFIED_DIGITS, CERTIFIED_DIGITS + 1] {
            let constant = engine.get_constant(&constant_symbol, requested).unwrap();
            assert_eq!(constant.requested_decimal_precision(), requested);
            assert_eq!(
                constant.available_precision(),
                AvailablePrecision::FixedPrecision(CERTIFIED_DIGITS)
            );
            assert_eq!(constant.provenance().version(), version);
            assert_eq!(constant.provenance().epsilon_power(), Some(0));
            assert_eq!(constant.provenance().source_sha256().len(), 64);
            check_value(
                &format!("{period} requested {requested}"),
                constant.numerical_value(),
                constant.precision(),
                reference,
                certified_error,
                false,
            );
            assert!(constant.to_atom().is_ok());
            check_rounding_warning(
                constant.warnings(),
                constant.precision(),
                requested,
                version,
                Some(period),
            );

            let master = engine.get_master(&master_symbol, requested).unwrap();
            assert_eq!(master.master_id(), master_id);
            assert_eq!(master.requested_decimal_precision(), requested);
            assert_eq!(master.provenance().version(), version);
            master.graph().validate().unwrap();
            let finite = &master.coefficients()[&0];
            assert_eq!(
                finite.available_precision(),
                AvailablePrecision::FixedPrecision(CERTIFIED_DIGITS)
            );
            assert_eq!(finite.provenance().version(), version);
            assert_eq!(finite.provenance().epsilon_power(), Some(0));
            check_value(
                &format!("{name} epsilon^0 requested {requested}"),
                finite.numerical_value(),
                finite.precision(),
                reference,
                certified_error,
                false,
            );
            assert!(finite.to_atom().is_ok());
            check_rounding_warning(
                finite.warnings(),
                finite.precision(),
                requested,
                version,
                None,
            );
            assert_eq!(
                constant.numerical_value().to_rational(),
                finite.numerical_value().to_rational()
            );
            if requested > CERTIFIED_DIGITS {
                check_cap_warning(constant.warnings(), version, Some(period), None);
                check_cap_warning(finite.warnings(), version, None, Some(master_id));
            }
        }
    }

    // Exact definitions containing a finite period inherit its capability;
    // neither their symbolic names nor unrelated archived fallbacks make
    // them arbitrary precision. Higher orders keep their own finite limits.
    for name in ["P19", "P73"] {
        let metadata = manifest["masters"]
            .as_array()
            .unwrap()
            .iter()
            .find(|master| master["name"].as_str() == Some(name))
            .unwrap();
        let symbol: MasterSymbol = metadata["master_id"].as_str().unwrap().parse().unwrap();
        let master = engine.get_master(&symbol, 32).unwrap();
        assert_eq!(
            master.coefficients()[&0].available_precision(),
            AvailablePrecision::FixedPrecision(CERTIFIED_DIGITS),
            "{name} epsilon^0"
        );
        assert_eq!(
            master.coefficients()[&1].available_precision(),
            AvailablePrecision::FixedPrecision(D2_DIGITS),
            "{name} epsilon^1"
        );
    }
}

#[test]
#[ignore = "loads independently validated D2 periods; set VAKINT_TEST_EXTERNAL_DATA"]
fn public_engine_independent_d2_periods_and_formula_inventory() {
    let root = std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
        .map(PathBuf::from)
        .expect("set VAKINT_TEST_EXTERNAL_DATA to the external-data checkout");
    let manifest: Value =
        serde_json::from_slice(&fs::read(root.join("4L/manifest.json")).unwrap()).unwrap();
    let version = manifest["version"].as_str().unwrap();
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(&root, 4).unwrap();

    // The inspection catalogue must expose every concrete constant used by
    // the loaded formulas, including the new independently defined periods.
    let inventory = engine.constant_symbols().unwrap();
    for symbol in &inventory {
        let constant = engine.get_constant(symbol, 16).unwrap();
        assert!(constant.to_atom().is_ok(), "{symbol}");
    }
    for name in ["PR7D2", "PR10D2", "PR15D2"] {
        let symbol: ConstantSymbol = name.parse().unwrap();
        assert!(inventory.contains(&symbol));
        for requested in [16, D2_DIGITS + 1] {
            let constant = engine.get_constant(&symbol, requested).unwrap();
            assert_eq!(
                constant.available_precision(),
                AvailablePrecision::FixedPrecision(D2_DIGITS)
            );
            assert_eq!(
                constant.precision().real().quality(),
                ErrorQuality::Estimated
            );
            assert!(constant.precision().imaginary().exact_zero());
            assert_eq!(constant.provenance().version(), version);
            assert!(constant.to_atom().is_ok());
            if requested > D2_DIGITS {
                assert!(constant.warnings().iter().any(|warning| {
                    warning.requested_decimal_precision == requested
                        && warning.data_version == version
                        && warning.sources.iter().any(|source| {
                            source.identity.contains(name)
                                && source.available_decimal_precision == Some(D2_DIGITS)
                        })
                }));
                assert!(constant.numerical_value().prec() <= decimal_bits(D2_DIGITS));
            }
        }
    }

    // These archived coefficients were computed independently of the exact
    // formulas (the 43 closures against the 500-digit physical values). Compare
    // every affected cell, not an expression solved for its period.
    let mut checked = 0;
    for metadata in manifest["masters"].as_array().unwrap() {
        let symbolic_path = metadata["files"]["symbolic.json"]["path"].as_str().unwrap();
        let formulas: Vec<(i64, String)> =
            serde_json::from_slice(&fs::read(root.join("4L").join(symbolic_path)).unwrap())
                .unwrap();
        let affected = formulas
            .iter()
            .filter(|(_, formula)| {
                ["PR7D2", "PR10D2", "PR15D2"]
                    .iter()
                    .any(|name| formula.contains(name))
            })
            .collect::<Vec<_>>();
        if affected.is_empty() {
            continue;
        }
        let numerical_path = metadata["files"]["numerical.json"]["path"]
            .as_str()
            .unwrap();
        let references: Vec<(i64, String)> =
            serde_json::from_slice(&fs::read(root.join("4L").join(numerical_path)).unwrap())
                .unwrap();
        let identity: MasterSymbol = metadata["master_id"].as_str().unwrap().parse().unwrap();
        let master = engine.get_master(&identity, 80).unwrap();
        for (power, _) in affected {
            let coefficient = &master.coefficients()[power];
            // A D2 period caps the formula; a better archived value may win.
            assert!(
                matches!(coefficient.available_precision(),
                    AvailablePrecision::FixedPrecision(digits) if digits >= D2_DIGITS),
                "{identity} epsilon^{power}: {}",
                coefficient.available_precision()
            );
            let reference = exact_decimal(&references.iter().find(|(p, _)| p == power).unwrap().1);
            let difference = coefficient.numerical_value().to_rational() - &reference;
            let distance = if difference < 0 {
                -difference
            } else {
                difference
            };
            let scale = if reference < 0 { -reference } else { reference };
            let tolerance = (Rational::from(1) + scale) * exact_decimal("1e-70");
            assert!(distance < tolerance, "{identity} epsilon^{power}");
            assert!(coefficient.precision().imaginary().exact_zero());
            assert!(coefficient.to_atom().is_ok());
            checked += 1;
        }
    }
    assert_eq!(
        checked, D2_DEPENDENT_CELLS,
        "the published D2-dependent cells must all be checked"
    );
    eprintln!(
        "checked {} constants and {checked} independent D2-dependent coefficients",
        inventory.len()
    );
}

#[test]
#[ignore = "checks the definitions of every four-loop constant; set VAKINT_TEST_EXTERNAL_DATA"]
fn public_engine_four_loop_constants_have_definitions() {
    use super::definitions;
    use crate::utils::vakint_macros::vk_symbol;
    use symbolica::atom::{Atom, AtomView};
    let root = std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
        .map(PathBuf::from)
        .expect("set VAKINT_TEST_EXTERNAL_DATA to the external-data checkout");
    let manifest: Value =
        serde_json::from_slice(&fs::read(root.join("4L/manifest.json")).unwrap()).unwrap();
    let families = manifest["families"].as_array().unwrap();
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(&root, 4).unwrap();

    for symbol in engine.constant_symbols().unwrap() {
        let constant = engine.get_constant(&symbol, 16).unwrap();
        let definition = constant
            .definition()
            .unwrap_or_else(|| panic!("{symbol} has no definition"));
        definitions::validate(definition.as_view())
            .unwrap_or_else(|error| panic!("{symbol}: {error}"));
    }

    // Every external record has a definition: its own symbolic_definition, or
    // for older packs the reconstruction from its definition metadata.
    let path = manifest["constant_bindings"]["path"].as_str().unwrap();
    let entries: Vec<Value> =
        serde_json::from_slice(&fs::read(root.join("4L").join(path)).unwrap()).unwrap();
    let mut recorded = 0;
    for entry in &entries {
        let symbol: ConstantSymbol = entry["symbol"].as_str().unwrap().parse().unwrap();
        let constant = engine.get_constant(&symbol, 16).unwrap();
        let definition = constant
            .definition()
            .unwrap_or_else(|| panic!("{symbol} has no definition"));
        if let Some(text) = entry["symbolic_definition"].as_str() {
            assert_eq!(definition, &vk_parse!(text).unwrap(), "{symbol}");
            recorded += 1;
        }
        if let Some(fallback) =
            definitions::external_fallback(symbol.as_atom(), &entry["provenance"], families)
                .unwrap()
        {
            assert!(
                definitions::agree(definition, &fallback),
                "{symbol}: {definition} versus metadata {fallback}"
            );
        }
    }
    println!(
        "{} external constants, {recorded} with symbolic_definition",
        entries.len()
    );

    for name in ["PR4fin", "PR4dfin", "PR11dep0", "PR4ep1", "PR4dep1"] {
        let symbol: ConstantSymbol = name.parse().unwrap();
        let constant = engine.get_constant(&symbol, 16).unwrap();
        assert!(
            matches!(constant.definition().unwrap().as_view(), AtomView::Fun(function)
            if function.get_symbol() == vk_symbol!("SeriesCoefficient"))
        );
    }

    // The data-independent FMFT table transcribes the manifest families and powers.
    let mut transcribed = 0;
    for (name, family_h, powers) in definitions::FMFT_MASTERS {
        let digest = if family_h {
            "eb225daea1ddc8a4d61675a26b9b487f90278190fe874e314b7253aee632a5e7"
        } else {
            "693719b1bc3d2f1d4d32c9312df66ced667a542591a2ada6f6498844fc5d0de1"
        };
        let master = manifest["masters"]
            .as_array()
            .unwrap()
            .iter()
            .find(|master| {
                master["family_sha256"] == digest
                    && master["powers"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|power| power.as_i64().unwrap())
                        .eq(powers)
            })
            .unwrap_or_else(|| panic!("{name}: no manifest master with these powers"));
        assert_eq!(
            vk_parse!(master["exact_basis_expression"].as_str().unwrap()).unwrap(),
            vk_parse!(name).unwrap(),
            "{name}"
        );
        let family = families
            .iter()
            .find(|family| family["sha256"] == digest)
            .unwrap();
        assert_eq!(
            definitions::EuclideanFamily::from_json(family)
                .unwrap()
                .public_master(&powers, 0)
                .unwrap(),
            definitions::fmft_master(name, 0).unwrap(),
            "{name}"
        );
        transcribed += 1;
    }
    assert_eq!(transcribed, definitions::FMFT_MASTERS.len());
}

/// Leaves of a published formula: every symbol and every outermost function
/// call. A Symbolica polylogarithm at a complex point is not a real constant;
/// its Vakint normalization splits it into real and imaginary projections,
/// whose leaves are collected instead.
fn collect_leaves(
    expression: symbolica::atom::AtomView,
    leaves: &mut std::collections::BTreeMap<symbolica::atom::Atom, usize>,
) {
    use symbolica::atom::{AtomCore, AtomView};
    expression.visitor(&mut |part| match part {
        AtomView::Num(_) => false,
        AtomView::Var(_) => {
            *leaves.entry(part.to_owned()).or_default() += 1;
            false
        }
        AtomView::Fun(function) => {
            if function.get_symbol() == symbolica::transcendental::polylog()
                && super::constants::builtin_definition(part)
                    .unwrap()
                    .is_none()
            {
                let projected = crate::hpl::normalize_constants(part).unwrap();
                assert_ne!(
                    projected.as_view(),
                    part,
                    "unrecognized polylogarithm {part}"
                );
                collect_leaves(projected.as_view(), leaves);
            } else {
                *leaves.entry(part.to_owned()).or_default() += 1;
            }
            false
        }
        _ => true,
    });
}

/// An explicit definition names no constant: only inert heads, bound
/// variables, exact numbers, `𝜋`, `𝑒`, `γ` and native functions occur. This
/// is checked independently of `definitions::validate`.
fn assert_explicit(leaf: &symbolica::atom::Atom, definition: &symbolica::atom::Atom) {
    use crate::utils::vakint_macros::vk_symbol;
    use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
    use symbolica::transcendental::{
        bessel_i, bessel_k, euler_gamma, gamma, imaginary_part, polylog, real_part, zeta,
    };
    let heads = [
        "Integral",
        "Interval",
        "LoopIntegral",
        "SeriesCoefficient",
        "EuclideanDot",
        "x",
        "k",
    ]
    .map(|name| vk_symbol!(name));
    let variables = ["t", "ep", "delta"].map(|name| vk_symbol!(name));
    let native = [
        real_part(),
        imaginary_part(),
        Symbol::EXP,
        Symbol::LOG,
        Symbol::SIN,
        gamma(),
        bessel_k(),
        bessel_i(),
    ];
    // The standard function itself at exact arguments is its own definition.
    if let AtomView::Fun(function) = definition.as_view()
        && [zeta(), polylog()].contains(&function.get_symbol())
    {
        assert!(
            function.iter().all(|argument| argument.is_constant()),
            "{leaf}: {definition}"
        );
        return;
    }
    definition.visitor(&mut |part| {
        match part {
            AtomView::Var(variable) => {
                let symbol = variable.get_symbol();
                assert!(
                    [Symbol::PI, Symbol::E, euler_gamma()].contains(&symbol)
                        || variables.contains(&symbol),
                    "{leaf}: definition references the symbol {part}: {definition}"
                );
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                assert!(
                    heads.contains(&symbol) || native.contains(&symbol),
                    "{leaf}: definition references the function {part}: {definition}"
                );
            }
            AtomView::Num(number) => assert!(
                !matches!(
                    number.get_coeff_view(),
                    symbolica::coefficient::CoefficientView::Float(..)
                ),
                "{leaf}: inexact number in {definition}"
            ),
            _ => {}
        }
        true
    });
}

/// Gate for the published coefficients: every leaf of every 1L–4L formula of
/// the external data has an explicit definition in Vakint, built in or from
/// its `symbolic_definition`, and a normalization statement.
///
/// `VAKINT_TEST_EXTRA_FORMULAS` may name a JSON array of further formula strings,
/// for example a basis rewrite prepared for a later data release, whose leaves
/// are checked against the same loaded data. `VAKINT_TEST_LEAF_DUMP` may name a
/// JSON file that receives the definition and normalization of every leaf.
#[test]
#[ignore = "checks every leaf of every published 1L-4L formula; set VAKINT_TEST_EXTERNAL_DATA"]
fn public_engine_every_published_leaf_has_an_explicit_definition() {
    use std::collections::BTreeMap;
    use symbolica::atom::AtomCore;
    let root = std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
        .map(PathBuf::from)
        .expect("set VAKINT_TEST_EXTERNAL_DATA to the external-data checkout");
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(&root, 4).unwrap();

    let mut leaves = BTreeMap::new();
    let mut formulas = 0;
    for loops in 1..=4 {
        let directory = root.join(format!("{loops}L"));
        let manifest: Value =
            serde_json::from_slice(&fs::read(directory.join("manifest.json")).unwrap()).unwrap();
        for master in manifest["masters"].as_array().unwrap() {
            let Some(path) = master["files"]["symbolic.json"]["path"].as_str() else {
                continue;
            };
            let rows: Vec<(i64, String)> =
                serde_json::from_slice(&fs::read(directory.join(path)).unwrap()).unwrap();
            for (power, formula) in rows {
                let expression = vk_parse!(&formula)
                    .unwrap_or_else(|error| panic!("{loops}L {path} ε^{power}: {error}"));
                collect_leaves(expression.as_view(), &mut leaves);
                formulas += 1;
            }
        }
    }
    if let Some(path) = std::env::var_os("VAKINT_TEST_EXTRA_FORMULAS") {
        let extra: Vec<String> = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        for formula in &extra {
            let expression =
                vk_parse!(formula).unwrap_or_else(|error| panic!("{formula}: {error}"));
            collect_leaves(expression.as_view(), &mut leaves);
        }
        eprintln!("{} further formulas from {}", extra.len(), path.display());
        formulas += extra.len();
    }
    let mut by_class = BTreeMap::<String, usize>::new();
    let mut dump = Vec::new();
    for leaf in leaves.keys() {
        let symbol = ConstantSymbol::try_from(leaf.clone())
            .unwrap_or_else(|error| panic!("{leaf} is not a constant symbol: {error}"));
        let constant = engine
            .get_constant(&symbol, 16)
            .unwrap_or_else(|error| panic!("{leaf} is not a registered constant: {error}"));
        let definition = constant.definition().unwrap_or_else(|| {
            panic!(
                "{leaf} has no definition: {}",
                constant.definition_unavailable_reason().unwrap_or_default()
            )
        });
        super::definitions::validate(definition.as_view())
            .unwrap_or_else(|error| panic!("{leaf}: {error}"));
        assert_explicit(leaf, definition);
        let normalisation = constant.normalisation();
        assert!(
            !normalisation.is_empty() && !normalisation.starts_with("No defining integral"),
            "{leaf}"
        );
        *by_class
            .entry(
                normalisation
                    .split(": ")
                    .next()
                    .and_then(|label| label.split(". ").next())
                    .unwrap()
                    .to_owned(),
            )
            .or_default() += 1;
        dump.push(serde_json::json!({
            "leaf": leaf.to_canonical_string(),
            "definition": definition.to_canonical_string(),
            "normalisation": normalisation,
        }));
    }
    if let Some(path) = std::env::var_os("VAKINT_TEST_LEAF_DUMP") {
        fs::write(path, serde_json::to_vec_pretty(&dump).unwrap()).unwrap();
    }
    eprintln!(
        "{formulas} published 1L-4L formulas, {} distinct leaves, every one explicitly defined",
        leaves.len()
    );
    for (class, count) in &by_class {
        eprintln!("  {count:4}  {class}");
    }
    assert!(formulas >= 739, "{formulas} formulas");
    assert!(leaves.len() >= 300, "{} leaves", leaves.len());
}
