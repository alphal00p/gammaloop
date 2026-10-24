//! Project independent vacuum factors with the existing native tensor engine.

use std::collections::BTreeSet;

use feynkit_tensor::TensorReducer;
use rustred::algebra::Coefficient;
use rustred::reduction::terminal_normalization::{
    VacuumTensorProjection, VacuumTensorProjectionInput,
};
use spenso::structure::representation::{BaseRepName, Minkowski};
use symbolica::atom::{AtomCore, AtomView, FunctionBuilder};
use symbolica::id::Replacement;
use symbolica::prelude::{Q, Z};
use symbolica::utils::resources;

use super::RustRedEvaluationError;
use crate::{Vakint, VakintTerm, initialize_spenso_reps, symbols::S};

fn error(error: impl std::error::Error + 'static) -> RustRedEvaluationError {
    match RustRedEvaluationError::resource_cause(&error) {
        Some(resource) => RustRedEvaluationError::ResourceLimit(resource.clone()),
        None => RustRedEvaluationError::ScalarNumerator {
            detail: error.to_string(),
        },
    }
}

pub(super) fn project(
    input: &VacuumTensorProjectionInput,
) -> Result<VacuumTensorProjection, RustRedEvaluationError> {
    resources::checkpoint()?;
    initialize_spenso_reps();
    let labels = input
        .component_loop_momenta
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    let substitutions = labels
        .iter()
        .enumerate()
        .map(|(axis, label)| Replacement::new(label.clone(), S.k.call(axis + 1)))
        .collect::<Vec<_>>();
    let mut numerator = input
        .numerator
        .replace_multiple(&substitutions)
        .replace_map(|term, _context, out| {
            if let AtomView::Fun(function) = term
                && function.get_symbol() == input.dot_head
            {
                **out = FunctionBuilder::new(S.dot)
                    .add_args(function.iter())
                    .finish();
            }
        });
    let representation = Minkowski::selfless_symbol().call(&input.dimension);
    let mut nonzero_conditions = BTreeSet::new();
    let mut offset = 0;
    for component in &input.component_loop_momenta {
        resources::checkpoint()?;
        let mut reducer = TensorReducer::new(input.dimension.clone());
        for axis in offset..offset + component.len() {
            reducer = reducer
                .with_integrated_vector(S.k.call((axis + 1, representation.clone())))
                .map_err(error)?;
        }
        offset += component.len();
        // Other blocks are numerator spectators. Registering them as external
        // denominator momenta would request a different, transverse projector.
        let indexed = Vakint::convert_from_dot_notation(numerator.as_view());
        let mut next_dummy = VakintTerm::next_dot_dummy_index(indexed.as_view());
        let opaque_slots = VakintTerm::collect_opaque_minkowski_slots(indexed.as_view());
        let opaque_vectors =
            VakintTerm::collect_opaque_compact_vectors(indexed.as_view(), &input.dimension);
        let native = VakintTerm::vakint_to_spenso_numerator(
            indexed.as_view(),
            &input.dimension,
            &opaque_slots,
        );
        let native = reducer.admit(native).map_err(error)?;
        let reduced = reducer.reduce(&native).map_err(error)?;
        log::debug!(target: "vakint::saved_application",
            "phase=factor_tensor_projected component_loops={} residual_rank={} terms={}",
            component.len(), reduced.max_rank(), reduced.terms().len());
        // Inspect native projector terms before materializing their sum, so a
        // later exact cancellation cannot erase a projector denominator.
        for term in reduced.terms() {
            let coefficient: Coefficient = term
                .coefficient()
                .try_to_rational_polynomial(&Q, &Z, None)
                .map_err(error)?;
            nonzero_conditions.insert(coefficient.denominator.to_expression());
        }
        let reduced = reduced.into_expression().map_err(error)?.into_expression();
        let bridged = VakintTerm::bridge_opaque_projectors(reduced.as_view(), &mut next_dummy);
        let bridged = VakintTerm::spenso_to_vakint_numerator(
            bridged.as_view(),
            &input.dimension,
            &opaque_vectors,
        );
        if VakintTerm::contains_feynkit_bridge_syntax(bridged.as_view()) {
            return Err(RustRedEvaluationError::ScalarNumerator {
                detail: "factor projection left internal tensor bridge syntax".into(),
            });
        }
        numerator = Vakint::convert_to_dot_notation(bridged.as_view()).expand();
    }
    let restore = labels
        .into_iter()
        .enumerate()
        .map(|(axis, label)| Replacement::new(S.k.call(axis + 1), label))
        .collect::<Vec<_>>();
    let numerator = numerator
        .replace_multiple(&restore)
        .replace_map(|term, _context, out| {
            if let AtomView::Fun(function) = term
                && function.get_symbol() == S.dot
            {
                **out = FunctionBuilder::new(input.dot_head)
                    .add_args(function.iter())
                    .finish();
            }
        });
    resources::checkpoint()?;
    Ok(VacuumTensorProjection {
        numerator,
        nonzero_conditions: nonzero_conditions.into_iter().collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;
    use symbolica::atom::Atom;

    #[test]
    fn native_projection_keeps_other_factors_as_spectators() {
        Vakint::initialize_vakint_symbols();
        let dimension = vk_parse!("d").unwrap();
        for (numerator, expected, loops) in [
            ("dot(k(1),k(2))", "0", 2),
            ("dot(k(1),k(2))^2", "dot(k(1),k(1))*dot(k(2),k(2))/d", 2),
            (
                "dot(k(1),k(2))*dot(k(1),k(3))*dot(k(2),k(3))",
                "dot(k(1),k(1))*dot(k(2),k(2))*dot(k(3),k(3))/d^2",
                3,
            ),
        ] {
            let result = project(&VacuumTensorProjectionInput {
                numerator: vk_parse!(numerator).unwrap(),
                dot_head: S.dot,
                component_loop_momenta: (1..=loops).map(|axis| vec![S.k.call(axis)]).collect(),
                dimension: dimension.clone(),
            })
            .unwrap();
            assert!(
                (result.numerator - vk_parse!(expected).unwrap())
                    .expand()
                    .together()
                    .is_zero()
            );
            if expected != "0" {
                assert!(result.nonzero_conditions.iter().any(|condition| {
                    condition
                        .replace(dimension.to_pattern())
                        .with(Atom::Zero.to_pattern())
                        .is_zero()
                }));
            }
        }
    }
}
