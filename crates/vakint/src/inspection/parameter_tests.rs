use std::collections::HashMap;

use symbolica::{
    atom::Atom,
    domains::float::{Complex, Float},
};

use super::*;
use crate::{
    VakintSettings, precision,
    utils::vakint_macros::{vk_parse, vk_symbol},
};

#[test]
fn registered_constant_parameters_fail_before_any_normalization() {
    let (mut catalog, master_symbol) = MasterCatalog::test_catalog(None, 20);
    let master = catalog.masters.get_mut(&master_symbol).unwrap();
    master.loop_count = 4;
    master.metadata = serde_json::json!({"exact_basis_expression": "PR4"});
    let coefficient = catalog.coefficient_symbol(&master_symbol, 0);
    let settings = VakintSettings::default();
    for identity in [
        vk_parse!("EulerGamma").unwrap(),
        vk_parse!("z3").unwrap(),
        vk_parse!("HR(0,1,om1)").unwrap(),
        vk_parse!("PR4fin").unwrap(),
        coefficient,
    ] {
        let parameters = HashMap::from([(
            identity.clone(),
            Complex::new(Float::with_val(53, 2), Float::new(53)),
        )]);
        // Even an unused override must not be silently accepted.
        for error in [
            precision::evaluate_series(
                Atom::one().as_view(),
                vk_symbol!("ep"),
                &parameters,
                &settings,
                Some(&catalog),
            )
            .unwrap_err(),
            precision::materialize(
                Atom::one().as_view(),
                &settings,
                &parameters,
                Some(&catalog),
            )
            .unwrap_err(),
        ] {
            assert!(
                matches!(
                    error,
                    VakintError::Inspection(InspectionError::ReservedParameter { .. })
                ),
                "{identity}: {error}"
            );
            assert!(error.to_string().contains("registered constant"));
        }
    }
}

#[test]
fn ordinary_qualified_and_unused_low_precision_parameters_remain_allowed() {
    crate::Vakint::initialize_vakint_symbols();
    let settings = VakintSettings {
        run_time_decimal_precision: 80,
        ..VakintSettings::default()
    };
    let ordinary = vk_parse!("user_parameters::EulerGamma").unwrap();
    let parameters = HashMap::from([(
        ordinary.clone(),
        Complex::new(Float::with_val(24, 7), Float::new(24)),
    )]);
    validate_parameter(ordinary.as_view(), None).unwrap();
    let result = precision::evaluate_series(
        Atom::one().as_view(),
        vk_symbol!("ep"),
        &parameters,
        &settings,
        None,
    )
    .unwrap();
    assert!(result.warnings().is_empty());
    assert_eq!(result.coefficients()[&0].value().re, Float::with_val(80, 1));
    let used = precision::evaluate_series(
        ordinary.as_view(),
        vk_symbol!("ep"),
        &parameters,
        &settings,
        None,
    )
    .unwrap();
    assert_eq!(used.coefficients()[&0].value().re, Float::with_val(80, 7));
    assert!(!used.warnings().is_empty());
}
