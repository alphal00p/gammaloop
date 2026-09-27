//! Typed interfaces for Symbolica's literal alias registry.
//!
//! Definitions live only in `AliasedAtom`. The accompanying structure records
//! logical layouts, which cannot in general be reconstructed from stored syntax.

use std::collections::{HashMap, HashSet};

use spenso::{
    network::tags::SPENSO_TAG,
    structure::{
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructure, PartialStructureExt},
        slot::IsAbstractSlot,
    },
};
use symbolica::{
    atom::{AliasedAtom, Atom, AtomCore, AtomView, FunctionBuilder},
    evaluate::EvaluatorBuilder,
};

use super::{
    SymbolicTensor,
    inference::{InterfaceInference, TensorInferenceError},
};

type Result<T> = std::result::Result<T, TensorInferenceError>;

/// Logical layouts accompanying one alias registry. Keys include port labelling.
#[derive(Clone, Debug)]
pub struct AliasInterfaces {
    root: PartialStructure,
    definitions: HashMap<Atom, (PartialStructure, PartialStructure)>,
}

impl AliasInterfaces {
    pub fn root(&self) -> &PartialStructure {
        &self.root
    }

    pub fn definitions(&self) -> &HashMap<Atom, (PartialStructure, PartialStructure)> {
        &self.definitions
    }
}

impl SymbolicTensor<PartialStructure> {
    /// Create an opaque handle with this tensor's physical ports and logical order.
    /// Each use with different labels needs its own literal alias definition.
    pub fn alias_handle(&self) -> Result<Self> {
        let owner = AbstractIndex::fresh_open_owner();
        let interface = if self.expression.is_zero() {
            self.structure.clone()
        } else {
            InterfaceInference::replacement_interface(self.expression.as_view())?
        };
        let slots = interface.logical_slots().into_iter().map(|slot| {
            let index = match slot.aind {
                PartialIndex::Explicit(index) => index,
                PartialIndex::Open(position) => AbstractIndex::Open {
                    owner,
                    axis: position.0,
                },
            };
            slot.rep().slot::<AbstractIndex, _>(index).to_atom()
        });
        let head = if self.is_scalar() {
            symbolica::symbol!("idenso::scalar_alias"; Scalar)
        } else {
            spenso::tensor_symbol!("idenso::tensor_alias")
        };
        let expression = FunctionBuilder::new(head)
            .add_arg(Atom::num(owner))
            .add_args(slots)
            .finish();
        Ok(Self::from_normalized_parts(
            expression,
            self.structure.clone(),
        ))
    }

    /// Attach checked literal definitions without resolving or expanding them.
    pub fn with_aliases(
        self,
        definitions: impl IntoIterator<Item = (Self, Self)>,
    ) -> Result<SymbolicTensor<AliasInterfaces, AliasedAtom>> {
        let mut expression = AliasedAtom::from(self.expression);
        let mut interfaces = HashMap::new();
        for (handle, body) in definitions {
            let valid_handle = match handle.expression.as_view() {
                AtomView::Var(variable) => {
                    variable.get_symbol().get_wildcard_level() == 0
                        && !variable.get_symbol().has_tag(&SPENSO_TAG.tensor)
                }
                AtomView::Fun(function) => {
                    function.get_symbol() == spenso::tensor_symbol!("idenso::tensor_alias")
                        || function.get_symbol()
                            == symbolica::symbol!("idenso::scalar_alias"; Scalar)
                        || function.get_symbol() == SPENSO_TAG.scalar
                }
                _ => false,
            };
            if !valid_handle {
                return Err(TensorInferenceError::invalid(
                    "tensor aliases require an opaque alias handle or a scalar variable",
                ));
            }
            if !InterfaceInference::additive_interfaces_match(&handle.structure, &body.structure) {
                return Err(TensorInferenceError::invalid(
                    "alias handle and body have different logical interfaces",
                ));
            }
            let encoded = InterfaceInference::replacement_interface(handle.expression.as_view())?;
            Self::validate_encoded_interface(&body.expression, &encoded)?;
            Self::validate_atom(&body.expression)?;
            if let Some(previous) = expression.get_aliases().get(&handle.expression) {
                if previous != &body.expression
                    || interfaces.get(&handle.expression)
                        != Some(&(handle.structure.clone(), body.structure.clone()))
                {
                    return Err(TensorInferenceError::invalid(
                        "conflicting alias definitions",
                    ));
                }
                continue;
            }
            interfaces.insert(
                handle.expression.clone(),
                (handle.structure, body.structure),
            );
            expression.register_alias(handle.expression, body.expression);
        }
        let value = SymbolicTensor {
            expression,
            structure: AliasInterfaces {
                root: self.structure,
                definitions: interfaces,
            },
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        };
        value.dependency_order()?;
        Ok(value)
    }
}

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    /// The unresolved root, with its established public port order.
    pub fn root(&self) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor {
            expression: self.expression.get_root().clone(),
            structure: self.structure.root.clone(),
            is_metric: self.is_metric,
            is_composite: self.is_composite,
        }
    }

    /// Typed definitions in dependency order. No body is expanded or re-inferred.
    pub fn aliases(
        &self,
    ) -> Result<
        Vec<(
            SymbolicTensor<PartialStructure>,
            SymbolicTensor<PartialStructure>,
        )>,
    > {
        self.dependency_order()?
            .into_iter()
            .map(|handle| {
                let (interface, body_interface) =
                    self.structure.definitions.get(handle).ok_or_else(|| {
                        TensorInferenceError::invalid("alias registry has no logical interface")
                    })?;
                Ok((
                    SymbolicTensor::from_normalized_parts(handle.clone(), interface.clone()),
                    SymbolicTensor::from_normalized_parts(
                        self.expression.get_aliases()[handle].clone(),
                        body_interface.clone(),
                    ),
                ))
            })
            .collect()
    }

    /// Resolve literal aliases, retaining factorization and checking callback results.
    pub fn resolved(&self) -> Result<SymbolicTensor<PartialStructure>> {
        self.dependency_order()?;
        self.root()
            .with_checked_expression(self.expression.clone().into_inner())
    }

    /// Explicit materialization. Contraction and alias construction never call this.
    pub fn expanded(&self) -> Result<SymbolicTensor<PartialStructure>> {
        self.resolved()?.expanded(None, true)
    }

    /// Apply a typed operation once to each definition; handles remain opaque.
    pub fn map_aliases(
        &self,
        mut map: impl FnMut(
            &SymbolicTensor<PartialStructure>,
            SymbolicTensor<PartialStructure>,
        ) -> Result<SymbolicTensor<PartialStructure>>,
    ) -> Result<Self> {
        let aliases = self
            .aliases()?
            .into_iter()
            .map(|(handle, body)| {
                let body = map(&handle, body)?;
                Ok((handle, body))
            })
            .collect::<Result<Vec<_>>>()?;
        self.root().with_aliases(aliases)
    }

    /// Compile the DAG directly using Symbolica's existing evaluator builder.
    pub fn evaluator<P: AtomCore>(&self, parameters: &[P]) -> Result<EvaluatorBuilder<'_>> {
        self.dependency_order()?;
        AliasedAtom::evaluator_multiple(std::slice::from_ref(&self.expression), parameters)
            .map_err(|error| TensorInferenceError::invalid(error.to_string()))
    }

    // Pinned Symbolica register_alias does not reject cycles. Check literal
    // dependencies before resolution/evaluation, including disconnected definitions.
    fn dependency_order(&self) -> Result<Vec<&Atom>> {
        let definitions = self.expression.get_aliases();
        if definitions.len() != self.structure.definitions.len()
            || definitions
                .keys()
                .any(|key| !self.structure.definitions.contains_key(key))
        {
            return Err(TensorInferenceError::invalid(
                "alias registry and interfaces disagree",
            ));
        }
        let mut keys = definitions.keys().collect::<Vec<_>>();
        keys.sort_unstable_by(|a, b| AtomView::cmp(&a.as_view(), &b.as_view()));
        let mut remaining = HashMap::new();
        let mut dependents: HashMap<&Atom, Vec<&Atom>> = HashMap::new();
        for &key in &keys {
            let mut dependencies = HashSet::new();
            definitions[key].visitor(&mut |node| {
                if let Some((handle, _)) = definitions.get_key_value(node.get_data()) {
                    dependencies.insert(handle);
                    return false;
                }
                true
            });
            remaining.insert(key, dependencies.len());
            for dependency in dependencies {
                dependents.entry(dependency).or_default().push(key);
            }
        }
        let mut order = keys
            .into_iter()
            .filter(|key| remaining[key] == 0)
            .collect::<Vec<_>>();
        let mut position = 0;
        while position < order.len() {
            if let Some(dependents) = dependents.get(order[position]) {
                for &dependent in dependents {
                    let count = remaining.get_mut(dependent).unwrap();
                    *count -= 1;
                    if *count == 0 {
                        order.push(dependent);
                    }
                }
            }
            position += 1;
        }
        if order.len() != definitions.len() {
            return Err(TensorInferenceError::invalid(
                "cyclic tensor alias definitions",
            ));
        }
        Ok(order)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use symbolica::symbol;

    fn tensor() -> SymbolicTensor<PartialStructure> {
        crate::test_support::test_initialize();
        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        SymbolicTensor::infer(
            FunctionBuilder::new(spenso::tensor_symbol!("alias_test::tensor"))
                .add_args([91001, 91003].map(|index| {
                    rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                        .to_atom()
                }))
                .finish(),
        )
        .unwrap()
    }

    #[test]
    fn literal_aliases_preserve_each_logical_layout_and_typed_zero() {
        let body = tensor();
        let mut handle = body.alias_handle().unwrap();
        handle.structure = PartialStructure::from_logical_slots(
            handle.structure.logical_slots().into_iter().rev(),
        );
        let root = handle.clone();
        let value = root
            .clone()
            .with_aliases([(handle.clone(), body.clone())])
            .unwrap();
        let pairs = value.aliases().unwrap();
        assert_eq!(pairs[0].0.structure, handle.structure);
        assert_eq!(pairs[0].1.structure, body.structure);
        assert_eq!(value.root(), root);
        let resolved = value.resolved().unwrap();
        assert_eq!(resolved.expression, body.expression);
        assert_eq!(resolved.structure, handle.structure);
        let zero = body.with_checked_expression(Atom::Zero).unwrap();
        let value = root.with_aliases([(handle, zero)]).unwrap();
        assert!(value.resolved().unwrap().expression.is_zero());
        assert_eq!(value.resolved().unwrap().structure, value.root().structure);
    }

    #[test]
    fn nested_aliases_resolve_without_expansion_and_map_once() {
        let scalar =
            |a| SymbolicTensor::checked_parts(a, PartialStructure::from_logical_slots([])).unwrap();
        let x = Atom::var(symbol!("alias_test::x"));
        let first = scalar(Atom::var(symbol!("alias_test::first")));
        let second = scalar(Atom::var(symbol!("alias_test::second")));
        let body = scalar((&x + 1).pow(3));
        let generated = body.alias_handle().unwrap();
        let scalar_alias = generated
            .clone()
            .with_aliases([(generated, body.clone())])
            .unwrap();
        assert_eq!(scalar_alias.resolved().unwrap(), body);
        let value = second
            .clone()
            .with_aliases([
                (second, scalar((&first.expression + 2).pow(2))),
                (first, body.clone()),
            ])
            .unwrap();
        assert_eq!(
            value.resolved().unwrap().expression,
            (&body.expression + 2).pow(2)
        );
        assert_eq!(
            value.expanded().unwrap().expression,
            (&body.expression + 2).pow(2).expand()
        );
        let mut calls = 0;
        let mapped = value
            .map_aliases(|_, body| {
                calls += 1;
                Ok(body)
            })
            .unwrap();
        assert_eq!(calls, 2);
        assert_eq!(mapped.resolved().unwrap(), value.resolved().unwrap());
        assert!(value.evaluator(&[x]).unwrap().build().is_ok());
    }

    #[test]
    fn alias_registry_rejects_cycles_conflicts_and_wrong_labels() {
        let scalar =
            |a| SymbolicTensor::checked_parts(a, PartialStructure::from_logical_slots([])).unwrap();
        let a = scalar(Atom::var(symbol!("alias_test::a")));
        let b = scalar(Atom::var(symbol!("alias_test::b")));
        assert!(a.clone().with_aliases([(a.clone(), a.clone())]).is_err());
        assert!(
            a.clone()
                .with_aliases([(a.clone(), b.clone()), (b.clone(), a.clone())])
                .is_err()
        );
        assert!(
            a.clone()
                .with_aliases([(a.clone(), scalar(Atom::one())), (a, scalar(Atom::num(2)))])
                .is_err()
        );
        let body = tensor();
        let handle = body.alias_handle().unwrap();
        let wrong = SymbolicTensor::infer(
            body.expression
                .replace(Atom::num(91001))
                .with(Atom::num(91005)),
        )
        .unwrap();
        assert!(
            handle
                .clone()
                .with_aliases([(handle.clone(), wrong)])
                .is_err()
        );
        let scalar_body = scalar(Atom::one());
        assert!(
            handle
                .clone()
                .with_aliases([(handle, scalar_body)])
                .is_err()
        );
    }

    #[test]
    fn alias_resolution_retains_encoded_open_ports() {
        crate::test_support::test_initialize();
        let body = SymbolicTensor::<PartialStructure>::from_signature(
            &crate::dirac::AGS.gamma_strct::<AbstractIndex>(4),
        )
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let value = handle
            .clone()
            .with_aliases([(handle, body.clone())])
            .unwrap();
        assert_eq!(value.resolved().unwrap().expression, body.expression);
        assert_eq!(value.root().structure, body.structure);

        let relabeled = body
            .reindex_interface_ports(&HashMap::from([(2, AbstractIndex::Normal(91201))]))
            .unwrap();
        let relabeled_handle = relabeled.alias_handle().unwrap();
        let uses = relabeled_handle
            .clone()
            .with_aliases([(value.root(), body), (relabeled_handle, relabeled.clone())])
            .unwrap();
        assert_eq!(uses.aliases().unwrap().len(), 2);
        assert_eq!(uses.resolved().unwrap(), relabeled);
    }

    #[test]
    fn alias_definitions_check_individual_branches_and_callback_rank_loss() {
        use crate::shorthands::schoonschip::Schoonschip;
        use spenso::network::library::symbolic::ETS;

        let body = tensor();
        let handle = body.alias_handle().unwrap();
        let wrong_branch = body
            .expression
            .replace(Atom::num(91001))
            .with(Atom::num(91005));
        let invalid = SymbolicTensor::from_normalized_parts(
            &body.expression + wrong_branch,
            body.structure.clone(),
        );
        assert!(handle.clone().with_aliases([(handle, invalid)]).is_err());

        let rep = LibraryRep::from(Minkowski {}).new_rep(4);
        let [a, b] = [91301, 91303].map(|index| {
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
                .to_atom()
        });
        let target = b.clone();
        let head = spenso::tensor_symbol!(
            "alias_test::callback_rank_loss",
            norm = move |node, out| {
                if let AtomView::Fun(function) = node
                    && function.iter().any(|argument| argument == target.as_view())
                {
                    **out = Atom::one();
                }
            }
        );
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&a)
            .add_arg(&b)
            .finish();
        let leaf = FunctionBuilder::new(head).add_arg(&a).finish();
        let body = SymbolicTensor::<PartialStructure>::infer(metric * leaf).unwrap();
        let handle = body.alias_handle().unwrap();
        let value = handle.clone().with_aliases([(handle, body)]).unwrap();
        assert!(
            value
                .map_aliases(|_, body| {
                    body.with_rewritten_expression(body.expression.schoonschip())
                })
                .is_err()
        );
    }

    #[test]
    fn borrowed_storage_keeps_the_original_payload_and_structure() {
        let body = tensor();
        let borrowed = SymbolicTensor {
            expression: body.expression.as_view(),
            structure: body.structure.clone(),
            is_metric: body.is_metric,
            is_composite: body.is_composite,
        };
        assert_eq!(
            borrowed.expression.get_data().as_ptr(),
            body.expression.as_view().get_data().as_ptr()
        );
        assert_eq!(borrowed.structure, body.structure);
    }
}
