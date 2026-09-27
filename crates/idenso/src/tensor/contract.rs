use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use symbolica::atom::{AliasedAtom, AtomView};

use super::{SymbolicTensor, aliases::AliasInterfaces, inference::TensorInferenceError};
use crate::shorthands::schoonschip::{Schoonschip, SlotContraction};

impl SymbolicTensor<PartialStructure> {
    /// Contract the metric/vector graph while keeping generated sums in aliases.
    /// An explicit order permutes the top-level factors after dot normalization.
    /// Materializing the result is a separate operation on the returned tensor.
    pub fn contract(
        &self,
        order: Option<&[usize]>,
    ) -> Result<SymbolicTensor<AliasInterfaces, AliasedAtom>, TensorInferenceError> {
        let source = self.expression.normalize_dots();
        if let Some(order) = order {
            let count = match source.as_view() {
                AtomView::Mul(product) => product.iter().len(),
                _ => 1,
            };
            let mut positions = order.to_vec();
            positions.sort_unstable();
            if positions != (0..count).collect::<Vec<_>>() {
                return Err(TensorInferenceError::Invalid(
                    "contraction order must list every normalized top-level factor exactly once"
                        .into(),
                ));
            }
        }
        if let Some(contracted) =
            SlotContraction::new().contract_factorized(source.as_view(), order)
        {
            let aliases = contracted
                .aliases
                .into_iter()
                .map(|(handle, body)| {
                    let scalar = PartialStructure::from_logical_slots([]);
                    Ok((
                        Self::checked_parts(handle, scalar.clone())?,
                        Self::checked_parts(body, scalar)?,
                    ))
                })
                .collect::<Result<Vec<_>, TensorInferenceError>>()?;
            // Normalizers on an emitted opaque tensor can invalidate a predicted
            // interface. Keep the established checked finishing boundary here.
            return self
                .with_rewritten_expression(contracted.root)?
                .with_aliases(aliases);
        }
        // The existing callback-sensitive contractor remains the owner of
        // non-intrinsic leaves. It preserves foreign structure and never asks
        // the bulk collector to expand an unsupported factor.
        self.with_rewritten_expression(source.schoonschip())?
            .with_aliases([])
    }
}

#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use super::*;
    use crate::shorthands::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use spenso::{g, mink, p};

    #[test]
    fn schoonschip_simplifies_metrics() {
        test_initialize();
        let expr = g!(mink!(4, mu), mink!(4, nu)) * p!(mink!(4, nu));

        assert_eq!(
            expr.schoonschip_net::<AbstractIndex>().unwrap(),
            p!(mink!(4, mu))
        );
    }

    #[test]
    fn contraction_returns_checked_aliases_and_preserves_typed_zero() {
        use spenso::structure::{
            partial::{PartialIndex, PartialStructureExt},
            representation::{LibraryRep, Minkowski, RepName},
            slot::IsAbstractSlot,
        };
        use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

        test_initialize();
        let p = spenso::vector_symbol!("contract_alias_p");
        let q = spenso::vector_symbol!("contract_alias_q");
        let r = spenso::vector_symbol!("contract_alias_r");
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let slot = rep.slot::<AbstractIndex, _>(92101).to_atom();
        let vector = |head| FunctionBuilder::new(head).add_arg(&slot).finish();
        let expression = (vector(p) + vector(q)) * vector(r);
        let source = SymbolicTensor::infer(expression).unwrap();
        let contracted = source.contract(None).unwrap();
        assert!(!contracted.expression.get_aliases().is_empty());
        assert_eq!(
            contracted.resolved().unwrap().expression.expand(),
            source
                .expression
                .schoonschip_with_settings(
                    &crate::shorthands::schoonschip::SchoonschipSettings::default()
                        .with_expanded_contracted_sums()
                )
                .expand()
        );
        let zero = SymbolicTensor::checked_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots([
                rep.slot(PartialIndex::Explicit(AbstractIndex::from(92103)))
            ]),
        )
        .unwrap();
        assert_eq!(zero.contract(None).unwrap().root(), zero);
    }

    #[test]
    fn contraction_checks_callback_induced_rank_loss() {
        use spenso::network::library::symbolic::ETS;
        use symbolica::atom::{Atom, AtomView, FunctionBuilder};

        test_initialize();
        let a = mink!(4, 92201);
        let b = mink!(4, 92203);
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "contract_alias_callback",
            norm = move |node, output| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **output = Atom::one();
                }
            }
        );
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let tensor = FunctionBuilder::new(head).add_arg(&a).finish();
        let source = SymbolicTensor::infer(metric * tensor).unwrap();
        assert!(source.contract(None).is_err());
    }
}
