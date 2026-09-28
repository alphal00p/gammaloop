//! Root rewriting and exact literal-use registration for the shared alias owner.

use super::*;
use crate::tensor::composition::{matching_interface_position, rewrite_interface_ports};
use spenso::structure::{
    representation::LibraryRep,
    slot::{Slot, SlotMatcher},
};

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// Apply one domain pass to the original root and reachable definitions.
    /// Handles remain opaque. Newly emitted definitions wait for the next pass;
    /// disconnected definitions are preserved without invoking the operation.
    /// The operation preserves each domain's logical order. Literal observations
    /// refer to existing registry entries; generated definitions are appended
    /// after those rewrites. All literals are checked before publication.
    pub(crate) fn map_domains(
        &self,
        mut map: impl FnMut(
            SymbolicTensor<PartialStructure>,
            &mut dyn FnMut(AtomView<'_>, &Atom),
        ) -> Result<(SymbolicTensor<PartialStructure>, Vec<Definition>)>,
    ) -> Result<Self> {
        let original = self.aliases()?;
        let mut definitions = original
            .iter()
            .cloned()
            .map(|pair| (pair.0.expression.clone(), pair))
            .collect::<HashMap<_, _>>();
        let reachable = self.reachable_definitions();
        let mut literal_uses = Vec::new();
        let mut observe = |before: AtomView<'_>, after: &Atom| {
            literal_uses.push((before.to_owned(), after.clone()));
        };
        let old_root = self.root();
        let (root, mut additional) = map(old_root.clone(), &mut observe)?;
        if old_root.structure != root.structure {
            return Err(TensorInferenceError::invalid(
                "alias domain rewrite changes its established root interface",
            ));
        }
        let mut changed = root != old_root;
        for (handle, body) in original {
            if reachable.contains(&handle.expression) {
                let (next, emitted) = map(body.clone(), &mut observe)?;
                if next.structure != body.structure {
                    return Err(TensorInferenceError::invalid(
                        "alias domain rewrite changes its established definition interface",
                    ));
                }
                changed |= next != body;
                definitions.get_mut(&handle.expression).unwrap().1 = next;
                additional.extend(emitted);
            }
        }
        if !changed && literal_uses.is_empty() && additional.is_empty() {
            return Ok(self.clone());
        }
        for (source, target) in literal_uses {
            Self::register_literal_use(&source, &target, &mut definitions)?;
        }
        for pair in additional {
            if let Some(previous) = definitions.get(&pair.0.expression)
                && previous != &pair
            {
                return Err(TensorInferenceError::invalid(
                    "conflicting domain alias definitions",
                ));
            }
            definitions.insert(pair.0.expression.clone(), pair);
        }
        root.with_aliases(definitions.into_values())
    }

    /// Finish a root rewrite with its recorded literal alias uses. The original
    /// handle is mandatory: two definitions can have the same head and owner
    /// yet different bodies, so the target spelling cannot identify its source.
    /// Nested port changes are recorded by the existing positional rewriter.
    pub(crate) fn with_rewritten_root(
        &self,
        root: SymbolicTensor<PartialStructure>,
        literal_uses: &[(Atom, Atom)],
        additional: impl IntoIterator<Item = Definition>,
    ) -> Result<Self> {
        if !InterfaceInference::additive_interfaces_match(&self.structure.root, &root.structure) {
            return Err(TensorInferenceError::invalid(
                "alias root rewrite changes its established interface",
            ));
        }
        let mut additional = additional.into_iter().peekable();
        if literal_uses.is_empty()
            && additional.peek().is_none()
            && root.expression == *self.expression.get_root()
            && root.structure == self.structure.root
            && root.is_metric == self.is_metric
            && root.is_composite == self.is_composite
        {
            return Ok(self.clone());
        }
        let mut definitions = self
            .aliases()?
            .into_iter()
            .map(|pair| (pair.0.expression.clone(), pair))
            .collect::<HashMap<_, _>>();
        for (source, target) in literal_uses {
            Self::register_literal_use(source, target, &mut definitions)?;
        }
        for pair in additional {
            if let Some(previous) = definitions.get(&pair.0.expression)
                && previous != &pair
            {
                return Err(TensorInferenceError::invalid(
                    "conflicting alias definitions",
                ));
            }
            definitions.insert(pair.0.expression.clone(), pair);
        }
        // This existing boundary checks cycles, unresolved handles, callback
        // results, and encoded interfaces before a sealed result is returned.
        root.with_aliases(definitions.into_values())
    }

    pub(in crate::tensor) fn register_literal_use(
        source: &Atom,
        target: &Atom,
        definitions: &mut HashMap<Atom, Definition>,
    ) -> Result<()> {
        if source == target {
            return Ok(());
        }
        let (handle, body) = definitions.get(source).cloned().ok_or_else(|| {
            TensorInferenceError::invalid("rewritten alias use has no registered source literal")
        })?;
        let (AtomView::Fun(before), AtomView::Fun(after)) = (source.as_view(), target.as_view())
        else {
            return Err(TensorInferenceError::invalid(
                "port rewriting requires a tensor alias",
            ));
        };
        if before.get_symbol() != spenso::tensor_symbol!("idenso::tensor_alias")
            || before.get_symbol() != after.get_symbol()
            || before.get_nargs() != after.get_nargs()
        {
            return Err(TensorInferenceError::invalid(
                "alias rewriting changed its head or arity",
            ));
        }
        let handle_slots = handle.structure.logical_slots();
        let body_slots = body.structure.logical_slots();
        let mut replacements = HashMap::new();
        let mut changed_slots = HashMap::new();
        let mut changed_handle_slots = HashMap::new();
        let mut handle_claimed = vec![false; handle_slots.len()];
        let mut body_claimed = vec![false; body_slots.len()];
        for (old, new) in before.iter().zip(after.iter()) {
            // Reuse the composition owner's positional mapping. A handle's
            // written open ports need not have its declared logical order, and
            // the body may carry a different valid presentation of that order.
            let position = matching_interface_position(old, &handle_slots, &handle_claimed);
            let body_position = matching_interface_position(old, &body_slots, &body_claimed);
            if let Some(position) = position {
                handle_claimed[position] = true;
            }
            if let Some(position) = body_position {
                body_claimed[position] = true;
            }
            if old == new {
                continue;
            }
            let position = position.ok_or_else(|| {
                TensorInferenceError::invalid("alias rewrite changed non-port metadata")
            })?;
            let body_position = body_position.ok_or_else(|| {
                TensorInferenceError::invalid("alias body has no matching logical port")
            })?;
            let slot = handle_slots[position];
            if body_slots
                .get(body_position)
                .is_none_or(|candidate| candidate.rep() != slot.rep())
            {
                return Err(TensorInferenceError::invalid(
                    "alias body port representation differs",
                ));
            }
            let next = if let Ok(slot) = Slot::<LibraryRep, AbstractIndex>::try_from(new) {
                Some(slot.rep().slot(PartialIndex::Explicit(slot.aind())))
            } else {
                let interface = InterfaceInference::replacement_interface(new)?;
                if !interface.logical_slots().is_empty() {
                    let mut matcher = SlotMatcher::default();
                    let compatible = if let AtomView::Fun(vector) = new
                        && vector.get_symbol().has_tag(&SPENSO_TAG.rank1)
                        && let Some(argument) = matcher.vector_argument(vector)
                        && let Ok(representation) =
                            matcher.parse_representation::<LibraryRep>(argument)
                    {
                        representation == slot.rep()
                            && matches!(interface.logical_slots().as_slice(), [bound]
                                if bound.rep() == representation && matches!(bound.aind, PartialIndex::Open(_)))
                    } else {
                        false
                    };
                    if !compatible {
                        return Err(TensorInferenceError::invalid(
                            "alias port binding is not a compatible compact vector",
                        ));
                    }
                }
                None
            };
            replacements.insert(body_position, new.to_owned());
            changed_slots.insert(body_position, next);
            changed_handle_slots.insert(position, next);
        }
        let interface =
            PartialStructure::from_logical_slots(body_slots.into_iter().enumerate().filter_map(
                |(position, slot)| changed_slots.get(&position).copied().unwrap_or(Some(slot)),
            ));
        let handle_interface =
            PartialStructure::from_logical_slots(handle_slots.into_iter().enumerate().filter_map(
                |(position, slot)| {
                    changed_handle_slots
                        .get(&position)
                        .copied()
                        .unwrap_or(Some(slot))
                },
            ));
        let mut nested = Vec::new();
        let expression = rewrite_interface_ports(&body, &replacements, &mut |old, new| {
            if matches!(old, AtomView::Fun(function) if function.get_symbol() == spenso::tensor_symbol!("idenso::tensor_alias")) {
                nested.push((old.to_owned(), new.clone()));
            }
        }).map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
        for (old, new) in nested {
            Self::register_literal_use(&old, &new, definitions)?;
        }
        SymbolicTensor::validate_encoded_interface(&expression, &interface)?;
        let body = SymbolicTensor::checked_parts(expression, interface.clone())?;
        let handle = SymbolicTensor::checked_parts(target.clone(), handle_interface)?;
        let pair = (handle, body);
        if let Some(previous) = definitions.get(target) {
            if previous != &pair {
                return Err(TensorInferenceError::invalid(
                    "conflicting relabelled alias definitions",
                ));
            }
        } else {
            definitions.insert(target.clone(), pair);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::symbol;

    fn scalar(expression: Atom) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::checked_parts(expression, PartialStructure::from_logical_slots([])).unwrap()
    }

    #[test]
    fn domain_pass_maps_original_reachable_bodies_once_and_retains_disconnected_definitions() {
        crate::test_support::test_initialize();
        let inner = scalar(Atom::var(symbol!("domain_map_x")) + Atom::one());
        let inner_handle = inner.alias_handle().unwrap();
        let outer = scalar(&inner_handle.expression * Atom::num(2));
        let outer_handle = outer.alias_handle().unwrap();
        let disconnected = scalar(Atom::var(symbol!("domain_map_disconnected")));
        let disconnected_handle = disconnected.alias_handle().unwrap();
        let value = outer_handle
            .clone()
            .with_aliases([
                (inner_handle, inner.clone()),
                (outer_handle.clone(), outer.clone()),
                (disconnected_handle.clone(), disconnected.clone()),
            ])
            .unwrap();
        let expected = value.resolved().unwrap();
        let mut visited = Vec::new();
        let mapped = value
            .map_domains(|body, _| {
                visited.push(body.expression.clone());
                let handle = body.alias_handle()?;
                Ok((handle.clone(), vec![(handle, body)]))
            })
            .unwrap();
        assert_eq!(visited.len(), 3, "new definitions wait until the next pass");
        assert_eq!(visited[0], outer_handle.expression);
        assert!(visited.contains(&inner.expression));
        assert!(visited.contains(&outer.expression));
        assert!(!visited.contains(&disconnected.expression));
        assert_eq!(mapped.resolved().unwrap(), expected);
        assert_eq!(mapped.aliases().unwrap().len(), 6);
        assert!(
            mapped
                .aliases()
                .unwrap()
                .contains(&(disconnected_handle, disconnected))
        );
    }

    #[test]
    fn domain_pass_preserves_identity_facts_and_rejects_body_rank_loss() {
        crate::test_support::test_initialize();
        let slot = spenso::mink!(4, 91701);
        let body = SymbolicTensor::infer(
            FunctionBuilder::new(spenso::tensor_symbol!("domain_rank_tensor"))
                .add_arg(slot)
                .finish(),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let mut value = handle
            .clone()
            .with_aliases([(handle, body.clone())])
            .unwrap();
        value.proofs.contracted = true;
        let identity = value.map_domains(|body, _| Ok((body, Vec::new()))).unwrap();
        assert!(identity.proofs.contracted);
        let original_root = value.expression.get_root().clone();
        let original_aliases = value.expression.get_aliases().clone();
        let original_structure = value.structure.clone();
        assert!(
            value
                .map_domains(|domain, _| {
                    if domain.expression == body.expression {
                        Ok((scalar(Atom::one()), Vec::new()))
                    } else {
                        Ok((domain, Vec::new()))
                    }
                })
                .is_err()
        );
        assert_eq!(value.expression.get_root(), &original_root);
        assert_eq!(value.expression.get_aliases(), &original_aliases);
        assert_eq!(value.structure.root(), original_structure.root());
        assert_eq!(
            value.structure.definitions(),
            original_structure.definitions()
        );
    }

    #[test]
    fn domain_pass_does_not_reorder_an_equivalent_explicit_interface() {
        crate::test_support::test_initialize();
        let root = SymbolicTensor::infer(
            FunctionBuilder::new(spenso::tensor_symbol!("domain_order_tensor"))
                .add_arg(spenso::mink!(4, 91711))
                .add_arg(spenso::mink!(4, 91713))
                .finish(),
        )
        .unwrap();
        let value = root.with_aliases([]).unwrap();
        assert!(
            value
                .map_domains(|mut domain, _| {
                    domain.structure = PartialStructure::from_logical_slots(
                        domain.structure.logical_slots().into_iter().rev(),
                    );
                    Ok((domain, Vec::new()))
                })
                .is_err()
        );
    }

    #[test]
    fn root_rewrite_registers_compact_binding_with_distinct_logical_layouts() {
        crate::test_support::test_initialize();
        let body = SymbolicTensor::infer(
            FunctionBuilder::new(spenso::tensor_symbol!("domain_layout_tensor"))
                .add_arg(spenso::mink!(4))
                .add_arg(spenso::euc!(6))
                .finish(),
        )
        .unwrap();
        let mut handle = body.alias_handle().unwrap();
        handle.structure = PartialStructure::from_logical_slots(
            handle.structure.logical_slots().into_iter().rev(),
        );
        let AtomView::Fun(function) = handle.expression.as_view() else {
            unreachable!()
        };
        let arguments = function
            .iter()
            .map(|argument| argument.to_owned())
            .collect::<Vec<_>>();
        let compact = FunctionBuilder::new(spenso::vector_symbol!("domain_layout_p"))
            .add_arg(spenso::mink!(4))
            .finish();
        let target = FunctionBuilder::new(function.get_symbol())
            .add_arg(&arguments[0])
            .add_arg(&compact)
            .add_arg(&arguments[2])
            .finish();
        let mut definitions = HashMap::from([(handle.expression.clone(), (handle.clone(), body))]);
        SymbolicTensor::register_literal_use(&handle.expression, &target, &mut definitions)
            .unwrap();
        let (next_handle, next_body) = &definitions[&target];
        assert_eq!(
            next_body.expression,
            FunctionBuilder::new(spenso::tensor_symbol!("domain_layout_tensor"))
                .add_arg(compact)
                .add_arg(spenso::euc!(6))
                .finish()
        );
        assert_eq!(
            next_handle.structure.logical_slots(),
            next_body.structure.logical_slots()
        );
        assert_eq!(next_body.structure.logical_slots().len(), 1);
    }
}
