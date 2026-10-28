use std::sync::Arc;

use color_eyre::Result;
use symbolica::evaluate::InliningPolicy;
use symbolica::prelude::*;

use crate::cff::CutCFFIndex;
use crate::initialisation;
use crate::integrands::process::param_builder::FnMapEntry;
use crate::uv::approx::direct_3d::DirectResidueBranches;
use crate::uv::{Integrands, ParametricIntegrandRow, ParametricIntegrandTerm};

#[test]
fn parametric_integrands_bind_coefficient_dependent_carriers_through_products() -> Result<()> {
    initialisation::test_initialise()?;
    let cut = CutCFFIndex::new_all_none();
    let coefficient = parse_lit!(parametric_binding_test::q);
    let term = |values: &[i64]| -> Result<Integrands> {
        let (_, tag) = DirectResidueBranches::numerator_scope();
        let numerator = Arc::new(FnMapEntry {
            lhs: function!(
                symbol!("gammalooprs::uv::numerator_family"),
                &tag,
                &coefficient
            ),
            rhs: coefficient.clone(),
            args: vec![Indeterminate::try_from(coefficient.clone()).unwrap()],
            tags: vec![tag],
            inlining: InliningPolicy::Never,
            is_alias: false,
        });
        let rows = values
            .iter()
            .map(|value| ParametricIntegrandRow {
                coefficients: vec![Rational::from(*value)],
                carrier: coefficient.clone(),
                source_keys: Vec::new(),
            })
            .collect();
        Integrands::from_iter([(cut, Atom::Zero)])
            .with_parametric_terms([(cut, ParametricIntegrandTerm::new(numerator, rows)?)])
    };
    let left = term(&[2, 3])?;
    let right = term(&[5, 7])?;
    // Each row contributes q*q: (4+9)*(25+49), with independent q scopes.
    assert_eq!(left.resolved()?.atom(&cut), Some(&Atom::num(13)));
    assert_eq!(right.resolved()?.atom(&cut), Some(&Atom::num(74)));
    let product = left.zip_mul(&right)?;
    assert_eq!(product.parametric_terms(&cut)[0].rows.len(), 4);
    assert_eq!(product.resolved()?.atom(&cut), Some(&Atom::num(962)));
    let triple = product.zip_mul(&term(&[11])?)?;
    assert_eq!(triple.resolved()?.atom(&cut), Some(&Atom::num(116_402)));
    Ok(())
}
