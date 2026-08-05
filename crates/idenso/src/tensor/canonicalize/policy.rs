use spenso::{
    network::{
        TensorNetworkError,
        library::DummyLibrary,
        parsing::{AtomStructureExt, ParseSettings, StrictTensorFilter},
    },
    structure::slot::{AbsInd, DummyAind, ParseableAind},
};
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{CanonicalizationError, validate_tensor_symmetry};
use crate::tensor::{SymbolicNet, SymbolicTensor};

/// A symbolic network paired with the exact Atom accepted by the canonical parser policy.
///
/// Private fields and the absence of a raw-network constructor make parser provenance a
/// type invariant. General parser users may still build a [`SymbolicNet`], but only this
/// canonical-policy parse can create a canonicalization input.
pub(crate) struct CanonicalPolicyNet<Aind: AbsInd> {
    network: SymbolicNet<Aind>,
    normalized_atom: Atom,
}

impl<Aind> CanonicalPolicyNet<Aind>
where
    Aind: AbsInd + DummyAind + ParseableAind,
{
    /// Validate and parse one normalized Atom with the fixed canonical policy.
    pub(crate) fn parse(normalized_atom: Atom) -> Result<Self, CanonicalizationError> {
        validate_tensor_symmetry::<Aind>(normalized_atom.as_view())?;
        Self::validate_power_grammar(normalized_atom.as_view())?;

        let library = DummyLibrary::<SymbolicTensor<Aind>>::new();
        let network = SymbolicNet::<Aind>::try_from_view::<SymbolicTensor<Aind>, _>(
            normalized_atom.as_view(),
            &library,
            &ParseSettings::default(),
        )
        .map_err(|error| match error {
            TensorNetworkError::NegativeExponentNonScalar(context) => {
                CanonicalizationError::NegativeExponentNonScalar(context)
            }
            TensorNetworkError::NonSelfDualTensorPower(context) => {
                CanonicalizationError::NonSelfDualTensorPower(context)
            }
            error => CanonicalizationError::Network(error.to_string()),
        })?;

        Ok(Self {
            network,
            normalized_atom,
        })
    }

    fn validate_power_grammar(value: AtomView<'_>) -> Result<(), CanonicalizationError> {
        if let AtomView::Pow(power) = value {
            let (base, exponent) = power.get_base_exp();
            if i8::try_from(exponent).is_err()
                && base.contains_exposed_tensor_topology(StrictTensorFilter::Tagged)
            {
                return Err(CanonicalizationError::UnsupportedTensorPowerExponent {
                    expression: value.to_owned(),
                    base: base.to_owned(),
                    exponent: exponent.to_owned(),
                });
            }
        }

        // Explicit parser opacity boundaries own their complete payload. Once
        // one hides tensor topology, Power grammar inside it is scalar syntax
        // rather than a network Power that this policy must materialize.
        if matches!(value, AtomView::Fun(_))
            && !value.contains_exposed_tensor_topology(StrictTensorFilter::Tagged)
        {
            return Ok(());
        }

        for child in value.children() {
            Self::validate_power_grammar(child)?;
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn into_network(self) -> SymbolicNet<Aind> {
        self.network
    }

    pub(crate) fn into_atom(self) -> Atom {
        self.normalized_atom
    }
}

impl<Aind: AbsInd> CanonicalPolicyNet<Aind> {
    pub(crate) fn network(&self) -> &SymbolicNet<Aind> {
        &self.network
    }

    pub(crate) fn normalized_atom(&self) -> &Atom {
        &self.normalized_atom
    }
}

#[cfg(test)]
mod tests {
    use spenso::{
        antisym, bracket,
        network::{
            graph::{NetworkLeaf, NetworkNode, NetworkOp},
            tags::SPENSO_TAG,
        },
        slot,
        structure::{
            abstract_index::{AIND_SYMBOLS, AbstractIndex},
            slot::IsAbstractSlot,
        },
        sym, tensor_symbol,
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder, NamespacedSymbol, SymbolAttribute, SymbolBuilder},
        domains::rational::Rational,
        function, symbol,
    };

    use super::*;
    use crate::test_support::test_initialize;

    fn tensorial_base() -> Atom {
        let rep = test_initialize().mink4;
        let tensor = tensor_symbol!(canonical_policy_power_tensor);
        function!(tensor, slot!(rep, canonical_policy_power_index).to_atom())
    }

    #[test]
    fn parse_retains_exact_source_and_every_sum_term() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_policy_sum_index);
        let left = tensor_symbol!(canonical_policy_sum_left);
        let right = tensor_symbol!(canonical_policy_sum_right);
        let expression = function!(left, index.to_atom()) + function!(right, index.to_atom());

        let parsed = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
        assert_eq!(parsed.normalized_atom(), &expression);
        assert_eq!(parsed.network().store.tensors.len(), 2);
        assert!(
            parsed
                .network()
                .graph
                .graph
                .iter_nodes()
                .any(|(_, _, node)| { matches!(node, NetworkNode::Op(NetworkOp::Sum)) })
        );

        assert_eq!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone())
                .unwrap()
                .into_atom(),
            expression
        );
        assert_eq!(parsed.into_network().store.tensors.len(), 2);
    }

    #[test]
    fn explicit_minus_one_product_keeps_the_parser_owned_scalar_leaf() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_policy_minus_one_index);
        let tensor = tensor_symbol!(canonical_policy_minus_one_tensor);
        let tensor_atom = function!(tensor, index.to_atom());
        let expression = Atom::num(-1) * tensor_atom.clone();
        let assert_parser_form = |parsed: &CanonicalPolicyNet<AbstractIndex>| {
            assert_eq!(parsed.network().store.scalar, vec![Atom::num(-1)]);
            assert_eq!(parsed.network().store.tensors.len(), 1);
            assert_eq!(parsed.network().store.tensors[0].expression, tensor_atom);
            assert!(
                parsed
                    .network()
                    .graph
                    .graph
                    .iter_nodes()
                    .any(|(_, _, node)| { matches!(node, NetworkNode::Op(NetworkOp::Product)) })
            );
        };

        let parsed = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
        assert_eq!(parsed.normalized_atom(), &expression);
        assert_parser_form(&parsed);

        let canonical = parsed.canonize(AbstractIndex::Dummy).unwrap();
        assert_eq!(canonical.normalized_atom(), &expression);
        assert_parser_form(&canonical);
    }

    #[test]
    fn rejects_non_i8_exponents_only_for_tagged_tensorial_bases() {
        let base = tensorial_base();
        let exponents = [
            Atom::num(Rational::from((1, 2))),
            Atom::var(symbol!("canonical_policy_symbolic_exponent")),
            Atom::num(128),
        ];

        for exponent in exponents {
            let expression = base.clone().pow(exponent.clone());
            assert!(matches!(
                CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
                Err(CanonicalizationError::UnsupportedTensorPowerExponent {
                    expression: rejected,
                    base: rejected_base,
                    exponent: rejected_exponent,
                }) if rejected == expression
                    && rejected_base == base
                    && rejected_exponent == exponent
            ));
        }
    }

    #[test]
    fn rejects_non_i8_exponents_around_unary_wrapped_tensor_topology() {
        let inner = tensorial_base();
        let wrapper = symbol!("canonical_policy_tensor_wrapper");
        let base = function!(wrapper, inner);
        let exponents = [
            Atom::num(Rational::from((1, 2))),
            Atom::var(symbol!("canonical_policy_wrapped_exponent")),
            Atom::num(128),
        ];

        for exponent in exponents {
            let expression = base.clone().pow(exponent.clone());
            assert!(matches!(
                CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
                Err(CanonicalizationError::UnsupportedTensorPowerExponent {
                    expression: rejected,
                    base: rejected_base,
                    exponent: rejected_exponent,
                }) if rejected == expression
                    && rejected_base == base
                    && rejected_exponent == exponent
            ));
        }
    }

    #[test]
    fn accepts_opaque_scalar_powers_under_the_tagged_policy() {
        let rep = test_initialize().mink4;
        let untagged = symbol!("canonical_policy_untagged_head");
        let opaque_base = function!(
            untagged,
            slot!(rep, canonical_policy_opaque_index).to_atom()
        );
        let scalar_base = Atom::var(symbol!("canonical_policy_scalar"));
        let exponents = [
            Atom::num(Rational::from((1, 2))),
            Atom::var(symbol!("canonical_policy_opaque_exponent")),
            Atom::num(128),
        ];

        for base in [opaque_base, scalar_base] {
            for exponent in &exponents {
                let expression = base.clone().pow(exponent.clone());
                let parsed =
                    CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
                assert_eq!(parsed.normalized_atom(), &expression);
                assert!(parsed.network().store.tensors.is_empty());
            }
        }
    }

    #[test]
    fn scalar_i8_min_power_remains_one_opaque_leaf() {
        let expression =
            Atom::var(symbol!("canonical_policy_scalar_i8_min_base")).pow(Atom::num(i8::MIN));
        let parsed = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();

        assert_eq!(parsed.normalized_atom(), &expression);
        assert_eq!(parsed.network().store.scalar, vec![expression]);
        assert!(parsed.network().store.tensors.is_empty());
        assert_eq!(
            parsed
                .network()
                .graph
                .graph
                .iter_nodes()
                .filter(|(_, _, node)| matches!(node, NetworkNode::Leaf(NetworkLeaf::Scalar(_))))
                .count(),
            1
        );
        assert!(
            parsed
                .network()
                .graph
                .graph
                .iter_nodes()
                .all(|(_, _, node)| !matches!(node, NetworkNode::Op(NetworkOp::Power(_))))
        );
    }

    #[test]
    fn accepts_non_i8_powers_across_explicit_opaque_boundaries() {
        let tensor = tensorial_base();
        let wrapper = symbol!("canonical_policy_opaque_wrapper");
        let bases = [
            function!(wrapper, bracket!(tensor.clone())),
            function!(wrapper, sym!(tensor)),
        ];

        for base in bases {
            let expression = base.clone().pow(Atom::num(Rational::from((1, 2))));
            let parsed = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
            assert_eq!(parsed.normalized_atom(), &expression);
            assert!(parsed.network().store.tensors.is_empty());
        }
    }

    #[test]
    fn accepts_non_i8_powers_inside_explicit_opaque_boundaries() {
        let hidden = tensorial_base().pow(Atom::num(Rational::from((1, 2))));
        let wrapper = symbol!("canonical_policy_nested_opaque_wrapper");
        let expressions = [
            function!(wrapper, bracket!(hidden.clone())),
            function!(wrapper, sym!(hidden)),
        ];

        for expression in expressions {
            let parsed = CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()).unwrap();
            assert_eq!(parsed.normalized_atom(), &expression);
            assert!(parsed.network().store.tensors.is_empty());
        }
    }

    #[test]
    fn preserves_negative_odd_open_self_dual_power_error() {
        let base = tensorial_base();
        let expression = base.pow(-3);

        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
            Err(CanonicalizationError::NegativeExponentNonScalar(context))
                if context.contains(&expression.to_plain_string())
        ));
    }

    #[test]
    fn preserves_non_self_dual_tensor_power_error() {
        let rep = test_initialize().cof_nc;
        let tensor = tensor_symbol!(canonical_policy_non_self_dual_power_tensor);
        let base = function!(
            tensor,
            slot!(rep, canonical_policy_non_self_dual_power_index).to_atom()
        );
        let expression = base.pow(2);

        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
            Err(CanonicalizationError::NonSelfDualTensorPower(context))
                if context == expression.to_plain_string()
        ));
    }

    #[test]
    fn parse_always_applies_intrinsic_symmetry_validation() {
        let rep = test_initialize().mink4;
        let tensor = tensor_symbol!(canonical_policy_invalid_symmetric; Symmetric);
        let expression = function!(
            tensor,
            Atom::num(1),
            slot!(rep, canonical_policy_invalid_index).to_atom()
        );

        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression),
            Err(CanonicalizationError::InvalidIntrinsicArgument {
                head,
                argument: 0,
                ..
            }) if head == tensor
        ));
    }

    #[test]
    fn intrinsic_symmetry_rejects_every_non_slot_structural_form() {
        let rep = test_initialize().mink4;
        let first = slot!(rep, canonical_policy_intrinsic_grammar_first);
        let second = slot!(rep, canonical_policy_intrinsic_grammar_second);
        let tensor = tensor_symbol!(canonical_policy_intrinsic_grammar_tensor; Symmetric);
        let nested = function!(
            symbol!("canonical_policy_intrinsic_nested"),
            first.to_atom()
        );
        let bundle = function!(AIND_SYMBOLS.aind, first.to_atom(), second.to_atom());
        let partial_group = antisym!(first, second);

        for argument in [Atom::num(1), nested, bundle, partial_group] {
            let expression = function!(tensor, argument);
            assert!(matches!(
                CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
                Err(CanonicalizationError::InvalidIntrinsicArgument {
                    head,
                    argument: 0,
                    expression: rejected,
                }) if head == tensor && rejected == expression
            ));
        }
    }

    #[test]
    fn structural_partial_groups_require_direct_slot_members() {
        let rep = test_initialize().mink4;
        let slot = slot!(rep, canonical_policy_invalid_partial_group_slot);
        let tensor = tensor_symbol!(canonical_policy_invalid_partial_group_tensor);
        let expression = FunctionBuilder::new(tensor)
            .add_arg(antisym!(slot, Atom::num(1)))
            .finish();

        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
            Err(CanonicalizationError::InvalidPartialGroup {
                group,
                expression: rejected,
                ..
            }) if group == *super::super::ANTISYM && rejected == expression
        ));
    }

    #[test]
    fn dynamically_built_tagged_symbols_use_the_same_intrinsic_grammar() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_policy_dynamic_tagged_index);
        let runtime_name = String::from("idenso::canonical_policy_dynamic_symmetric");
        let tensor = SymbolBuilder::new(NamespacedSymbol::parse(&runtime_name))
            .with_attributes(vec![SymbolAttribute::Symmetric])
            .with_tags(std::slice::from_ref(&SPENSO_TAG.tensor))
            .build()
            .unwrap();
        let valid = function!(tensor, index.to_atom());

        let parsed = CanonicalPolicyNet::<AbstractIndex>::parse(valid.clone()).unwrap();
        assert_eq!(parsed.normalized_atom(), &valid);
        assert_eq!(parsed.network().store.tensors.len(), 1);

        let invalid = function!(tensor, Atom::num(1));
        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(invalid.clone()),
            Err(CanonicalizationError::InvalidIntrinsicArgument {
                head,
                argument: 0,
                expression: rejected,
            }) if head == tensor && rejected == invalid
        ));
    }

    #[test]
    fn overlapping_intrinsic_attributes_follow_symbolica_antisymmetric_precedence() {
        let tensor = SymbolBuilder::new(NamespacedSymbol::parse(
            "idenso::canonical_policy_overlapping_symmetry",
        ))
        .with_attributes(vec![
            SymbolAttribute::Symmetric,
            SymbolAttribute::Antisymmetric,
            SymbolAttribute::Cyclesymmetric,
        ])
        .with_tags(std::slice::from_ref(&SPENSO_TAG.tensor))
        .build()
        .unwrap();

        assert_eq!(
            super::super::SymmetryKind::of(tensor),
            Some(super::super::SymmetryKind::Antisymmetric)
        );
    }

    #[test]
    fn rejects_exposed_tensor_topology_hidden_in_an_ordered_argument() {
        let rep = test_initialize().mink4;
        let index = slot!(rep, canonical_policy_hidden_topology_index);
        let outer = tensor_symbol!(canonical_policy_hidden_topology_outer);
        let left = tensor_symbol!(canonical_policy_hidden_topology_left);
        let right = tensor_symbol!(canonical_policy_hidden_topology_right);
        let hidden = function!(left, index.to_atom()) * function!(right, index.to_atom());
        let parameter = Atom::var(symbol!("canonical_policy_hidden_topology_parameter"));
        let expression = function!(outer, parameter.clone(), hidden.clone());

        super::super::projection::reset_graphica_calls();
        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(expression.clone()),
            Err(CanonicalizationError::HiddenTensorTopology {
                head,
                argument: 1,
                expression: rejected,
            }) if head == outer && rejected == expression
        ));
        assert_eq!(super::super::projection::graphica_calls(), 0);

        let untagged = symbol!("canonical_policy_hidden_topology_untagged");
        let untagged_expression = function!(untagged, parameter.clone(), hidden.clone());
        assert!(matches!(
            CanonicalPolicyNet::<AbstractIndex>::parse(untagged_expression.clone()),
            Err(CanonicalizationError::HiddenTensorTopology {
                head,
                argument: 1,
                expression: rejected,
            }) if head == untagged && rejected == untagged_expression
        ));
        assert_eq!(super::super::projection::graphica_calls(), 0);

        let opaque_expression = function!(untagged, parameter.clone(), bracket!(hidden.clone()));
        let opaque = CanonicalPolicyNet::<AbstractIndex>::parse(opaque_expression.clone()).unwrap();
        assert_eq!(opaque.normalized_atom(), &opaque_expression);
        assert!(opaque.network().store.tensors.is_empty());

        let tagged_opaque_expression = function!(outer, parameter, bracket!(hidden));
        let tagged_opaque =
            CanonicalPolicyNet::<AbstractIndex>::parse(tagged_opaque_expression.clone()).unwrap();
        assert_eq!(tagged_opaque.normalized_atom(), &tagged_opaque_expression);
    }
}
