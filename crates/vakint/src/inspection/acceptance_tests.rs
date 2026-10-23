//! Public-engine acceptance of the independently certified four-loop periods.
//! This test deliberately loads the complete local pack, including its graph
//! and IBP metadata, and never contacts a producer workspace or the network.

use std::{fs, path::PathBuf};

use serde_json::Value;
use symbolica::domains::{float::Float, rational::Rational};

use crate::{
    AvailablePrecision, CoefficientPrecision, ConstantSymbol, ErrorQuality, MasterSymbol,
    PrecisionWarning, Vakint, precision::decimal_bits, utils::vakint_macros::vk_parse,
};

const CERTIFIED_DIGITS: u32 = 20_000;

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
            AvailablePrecision::FixedPrecision(200),
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
        for requested in [16, 216] {
            let constant = engine.get_constant(&symbol, requested).unwrap();
            assert_eq!(
                constant.available_precision(),
                AvailablePrecision::FixedPrecision(215)
            );
            assert_eq!(
                constant.precision().real().quality(),
                ErrorQuality::Estimated
            );
            assert!(constant.precision().imaginary().exact_zero());
            assert_eq!(constant.provenance().version(), version);
            assert!(constant.to_atom().is_ok());
            if requested == 216 {
                assert!(constant.warnings().iter().any(|warning| {
                    warning.requested_decimal_precision == requested
                        && warning.data_version == version
                        && warning.sources.iter().any(|source| {
                            source.identity.contains(name)
                                && source.available_decimal_precision == Some(215)
                        })
                }));
                assert!(constant.numerical_value().prec() <= decimal_bits(215));
            }
        }
    }

    // These archived coefficients were independently computed before the new
    // exact formulas were published. Compare all fifteen affected cells, not
    // an expression obtained by solving the new formula for its period.
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
            assert_eq!(
                coefficient.available_precision(),
                AvailablePrecision::FixedPrecision(215),
                "{identity} epsilon^{power}"
            );
            assert_eq!(
                coefficient.precision().real().quality(),
                ErrorQuality::Estimated
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
        checked, 15,
        "the published D2-dependent cells must all be checked"
    );
    eprintln!(
        "checked {} constants and {checked} independent D2-dependent coefficients",
        inventory.len()
    );
}
