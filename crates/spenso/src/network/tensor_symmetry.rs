//! Symmetries of disjoint argument groups on flat tensor heads.
//!
//! Group positions refer to function arguments, including any intervening scalar
//! metadata. Only the selected arguments are permuted. Symbolica normalizes each
//! group; no tensor expression is expanded into a sum of permutations.
use std::{
    collections::{BTreeSet, HashMap},
    sync::{Mutex, OnceLock},
};

use crate::{
    network::tags::{SPENSO_TAG, tensor_print},
    shadowing::FactorProjector,
    structure::{
        representation::{LibraryRep, Representation},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{
        Atom, AtomCore, AtomView, FunctionBuilder, NamespacedSymbol, Symbol, SymbolAttribute,
        SymbolBuilder, UserData, representation::FunView,
    },
    symbol,
};

const TAG_PREFIX: &str = "spenso::slot-symmetry:v1:";

/// Disjoint symmetric or antisymmetric argument groups of a tensor function.
///
/// The declaration is serialized in the symbol tags. Rust callbacks are not
/// serialized by Symbolica: restore the declaration with [`Self::register_symbol`]
/// before importing its atoms in a fresh process. Native admission rejects an
/// imported declaration whose normalization callback has not been registered.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TensorSlotSymmetry {
    groups: Vec<(FactorProjector, Vec<usize>)>,
}

fn registrations() -> &'static Mutex<HashMap<Symbol, TensorSlotSymmetry>> {
    static REGISTERED: OnceLock<Mutex<HashMap<Symbol, TensorSlotSymmetry>>> = OnceLock::new();
    REGISTERED.get_or_init(Mutex::default)
}

impl TensorSlotSymmetry {
    pub fn new(
        groups: impl IntoIterator<Item = (FactorProjector, Vec<usize>)>,
    ) -> Result<Self, String> {
        let mut groups: Vec<_> = groups.into_iter().collect();
        let mut positions = BTreeSet::new();
        for (kind, group) in &mut groups {
            if *kind == FactorProjector::Cyclic {
                return Err(
                    "tensor slot groups currently support only symmetric and antisymmetric groups"
                        .into(),
                );
            }
            if group.len() < 2 {
                return Err("a tensor slot symmetry group requires at least two arguments".into());
            }
            group.sort_unstable();
            for &position in group.iter() {
                if !positions.insert(position) {
                    return Err(format!(
                        "tensor slot symmetry groups must be disjoint; argument {position} occurs more than once"
                    ));
                }
            }
        }
        if groups.is_empty() {
            return Err("a tensor slot symmetry declaration requires at least one group".into());
        }
        groups.sort_unstable_by_key(|(_, group)| group[0]);
        Ok(Self { groups })
    }

    pub fn groups(&self) -> impl Iterator<Item = (FactorProjector, &[usize])> {
        self.groups
            .iter()
            .map(|(kind, group)| (*kind, group.as_slice()))
    }

    fn tag(&self) -> String {
        let groups: Vec<_> = self
            .groups()
            .map(|(kind, positions)| {
                (
                    match kind {
                        FactorProjector::Symmetric => "symmetric",
                        FactorProjector::Antisymmetric => "antisymmetric",
                        FactorProjector::Cyclic => unreachable!(),
                    },
                    positions,
                )
            })
            .collect();
        format!(
            "{TAG_PREFIX}{}",
            serde_json::to_string(&groups).expect("argument groups serialize")
        )
    }

    /// Read the complete persisted declaration without relying on a callback.
    pub fn from_symbol(symbol: Symbol) -> Result<Option<Self>, String> {
        Self::from_tags(symbol.get_tags())
    }

    fn from_tags(tags: &[String]) -> Result<Option<Self>, String> {
        if tags
            .iter()
            .any(|tag| tag.starts_with("spenso::slot-symmetry:") && !tag.starts_with(TAG_PREFIX))
        {
            return Err("unsupported tensor slot symmetry metadata version".into());
        }
        let mut declarations = tags.iter().filter_map(|tag| tag.strip_prefix(TAG_PREFIX));
        let Some(declaration) = declarations.next() else {
            return Ok(None);
        };
        if declarations.next().is_some() {
            return Err("tensor symbol has more than one slot symmetry declaration".into());
        }
        let groups: Vec<(String, Vec<usize>)> = serde_json::from_str(declaration)
            .map_err(|error| format!("invalid tensor slot symmetry metadata: {error}"))?;
        let groups = groups
            .into_iter()
            .map(|(kind, positions)| {
                Ok((
                    match kind.as_str() {
                        "symmetric" => FactorProjector::Symmetric,
                        "antisymmetric" => FactorProjector::Antisymmetric,
                        _ => return Err(format!("unsupported tensor slot symmetry group {kind}")),
                    },
                    positions,
                ))
            })
            .collect::<Result<Vec<_>, String>>()?;
        Self::new(groups).map(Some)
    }

    /// Restore a declaration from archived tags before importing its expressions.
    /// This reattaches the native callback from persisted scientific metadata;
    /// callers do not construct or serialize a normalization closure.
    pub fn restore_symbol(
        name: NamespacedSymbol,
        attributes: Vec<SymbolAttribute>,
        tags: Vec<String>,
        data: UserData,
    ) -> Result<Symbol, String> {
        let declaration =
            Self::from_tags(&tags)?.ok_or("archived tensor lacks a slot symmetry declaration")?;
        let tags = tags
            .into_iter()
            .filter(|tag| !tag.starts_with(TAG_PREFIX))
            .collect();
        declaration.register_symbol(name, attributes, tags, data)
    }

    /// Register a tensor declaration, retaining the caller's tags and structured
    /// data. Repeating an identical registration returns the established symbol.
    /// A previously imported head lacking its callback cannot be repaired in
    /// place; register first, then import the atom archive.
    pub fn register_symbol(
        &self,
        name: NamespacedSymbol,
        attributes: Vec<SymbolAttribute>,
        mut tags: Vec<String>,
        data: UserData,
    ) -> Result<Symbol, String> {
        if attributes.iter().any(|attribute| {
            matches!(
                attribute,
                SymbolAttribute::Symmetric
                    | SymbolAttribute::Antisymmetric
                    | SymbolAttribute::Cyclesymmetric
                    | SymbolAttribute::Linear
                    | SymbolAttribute::Flat
                    | SymbolAttribute::Scalar
            )
        }) {
            return Err("selected tensor slot symmetry cannot be combined with whole-head symmetry, Linear, Flat, or Scalar attributes".into());
        }
        if tags
            .iter()
            .any(|tag| tag.starts_with("spenso::slot-symmetry:"))
        {
            return Err("slot symmetry metadata is supplied by the native declaration".into());
        }
        if !tags.contains(&SPENSO_TAG.tensor) {
            tags.push(SPENSO_TAG.tensor.clone());
        }
        tags.push(self.tag());
        // Initialize helper symbols before taking the registration lock.
        let _ = (
            FactorProjector::Symmetric.symbol(),
            FactorProjector::Antisymmetric.symbol(),
        );
        let mut registered = registrations()
            .lock()
            .expect("tensor symmetry registration lock");
        if let Some(existing) = Symbol::get_symbol(name.clone()) {
            if registered.get(&existing) != Some(self)
                || existing.get_normalization_function().is_none()
            {
                return Err(format!(
                    "tensor symmetry callback for {} is unavailable; register its native declaration before importing atoms",
                    existing.get_name()
                ));
            }
            if existing.get_tags() != tags
                || existing.get_attributes().len() != attributes.len()
                || !attributes
                    .iter()
                    .all(|attribute| existing.get_attributes().contains(attribute))
                || existing.get_data() != &data
            {
                return Err(format!(
                    "symbol {} already has a different tensor declaration",
                    existing.get_name()
                ));
            }
            return Ok(existing);
        }
        let declaration = self.clone();
        let symbol = SymbolBuilder::new(name)
            .with_attributes(attributes)
            .with_tags(tags)
            .with_user_data(data)
            .with_print_function(tensor_print)
            .with_normalization_function(move |value, out| {
                if let AtomView::Fun(function) = value
                    && let Some(normalized) = declaration.normalize(function)
                {
                    **out = normalized;
                }
            })
            .build()
            .map_err(|error| error.to_string())?;
        registered.insert(symbol, self.clone());
        Ok(symbol)
    }

    /// Validate callback restoration and the concrete argument-group signature.
    pub fn validate_function(function: FunView<'_>) -> Result<(), String> {
        let symbol = function.get_symbol();
        let Some(declaration) = Self::from_symbol(symbol)? else {
            return Ok(());
        };
        if registrations()
            .lock()
            .expect("tensor symmetry registration lock")
            .get(&symbol)
            != Some(&declaration)
            || symbol.get_normalization_function().is_none()
        {
            return Err(format!(
                "tensor symmetry callback for {} is unavailable; register its native declaration before importing atoms",
                symbol.get_name()
            ));
        }
        let args = function.iter().collect::<Vec<_>>();
        let mut matcher = SlotMatcher::default();
        for (_, positions) in declaration.groups() {
            let mut representation = None;
            for &position in positions {
                let Some(&argument) = args.get(position) else {
                    return Err(format!(
                        "tensor symmetry argument {position} is outside a function with {} arguments",
                        args.len()
                    ));
                };
                if matches!(argument, AtomView::Var(variable)
                    if variable.get_symbol() == SPENSO_TAG.chain_in || variable.get_symbol() == SPENSO_TAG.chain_out)
                {
                    // The enclosing native matrix channel supplies this
                    // representation when it materializes its endpoints.
                    continue;
                }
                let rep = match matcher.classify(argument) {
                    SlotMatch::Explicit(slot) if !slot.is_concrete_index() => {
                        Representation::<LibraryRep>::try_from(argument)
                    }
                    _ => matcher.parse_representation::<LibraryRep>(argument),
                }
                .map_err(|_| {
                    format!("tensor symmetry argument {position} must be a direct tensor slot")
                })?;
                if representation.is_some_and(|previous| previous != rep) {
                    return Err("tensor symmetry arguments must have identical representations and dimensions".into());
                }
                representation = Some(rep);
            }
        }
        Ok(())
    }

    fn normalize(&self, function: FunView<'_>) -> Option<Atom> {
        let original = function.iter().collect::<Vec<_>>();
        let mut args = original
            .iter()
            .map(|argument| argument.to_owned())
            .collect::<Vec<_>>();
        let mut coefficient = Atom::num(1);
        let mut matcher = SlotMatcher::default();
        let mut selected_groups = Vec::new();
        for (kind, positions) in self.groups() {
            let mut representation: Option<Representation<LibraryRep>> = None;
            let mut selected = Vec::with_capacity(positions.len());
            for &position in positions {
                let argument = *original.get(position)?;
                // Occurrence-local open representations have no written index.
                // Materialization assigns them distinct identities first.
                let SlotMatch::Explicit(slot) = matcher.classify(argument) else {
                    return None;
                };
                if slot.is_concrete_index() {
                    return None;
                }
                // Index payload validation belongs to the caller's native
                // index parser, not a particular concrete label type here.
                let rep = Representation::<LibraryRep>::try_from(argument).ok()?;
                if representation.is_some_and(|previous| previous != rep) {
                    return None;
                }
                representation = Some(rep);
                selected.push(argument);
            }
            selected_groups.push((kind, positions, selected));
        }
        // A zero in one group must not hide an invalid later group.
        for (kind, positions, selected) in selected_groups {
            let normalized = kind.apply(selected);
            if normalized.is_zero() {
                return Some(Atom::Zero);
            }
            let (factor, normalized) = Self::group_parts(normalized.as_view(), kind.symbol())?;
            coefficient *= factor;
            for (&position, argument) in positions.iter().zip(normalized.iter()) {
                args[position] = argument.to_owned();
            }
        }
        if coefficient == 1
            && args
                .iter()
                .zip(&original)
                .all(|(new, old)| new.as_view() == *old)
        {
            return None;
        }
        let atom = args
            .into_iter()
            .fold(
                FunctionBuilder::new(function.get_symbol()),
                |builder, arg| builder.add_arg(arg),
            )
            .finish();
        Some(coefficient * atom)
    }

    fn group_parts(value: AtomView<'_>, head: Symbol) -> Option<(Atom, FunView<'_>)> {
        match value {
            AtomView::Fun(function) if function.get_symbol() == head => {
                Some((Atom::num(1), function))
            }
            AtomView::Mul(product) => {
                let mut coefficient = Atom::num(1);
                let mut group = None;
                for factor in product.iter() {
                    if let AtomView::Fun(function) = factor
                        && function.get_symbol() == head
                        && group.is_none()
                    {
                        group = Some(function);
                        continue;
                    }
                    if !matches!(factor, AtomView::Num(_)) {
                        return None;
                    }
                    coefficient *= factor;
                }
                group.map(|group| (coefficient, group))
            }
            _ => None,
        }
    }
}

crate::symbolica_init_lazy_static! {
    static ENVELOPE, ENVELOPE_INNER: Symbol = || symbol!("spenso::slot_symmetry_envelope"; Linear);
}
crate::symbolica_init_lazy_static! {
    static FIXED, FIXED_INNER: Symbol = || symbol!("spenso::slot_symmetry_fixed");
}

impl TensorSlotSymmetry {
    /// Expose slot subgroups to Symbolica's existing tensor canonicalizer using
    /// a compact temporary encoding. Fixed argument wrappers protect arbitrary
    /// scalar metadata from the envelope's linear coefficient extraction.
    /// This expression is an internal canonicalization view, not tensor input.
    pub fn canonical_view(function: FunView<'_>) -> Result<Atom, String> {
        let Some(declaration) = Self::from_symbol(function.get_symbol())? else {
            return Ok(function.as_view().to_owned());
        };
        Self::validate_function(function)?;
        let args = function.iter().collect::<Vec<_>>();
        let mut envelope = FunctionBuilder::new(*ENVELOPE).add_arg(
            FunctionBuilder::new(*FIXED)
                .add_arg(Atom::var(function.get_symbol()))
                .finish(),
        );
        for (position, argument) in args.iter().enumerate() {
            if let Some((kind, group)) = declaration
                .groups()
                .find(|(_, group)| group.contains(&position))
            {
                if position == group[0] {
                    envelope =
                        envelope.add_arg(kind.apply(group.iter().map(|&position| args[position])));
                }
            } else {
                envelope =
                    envelope.add_arg(FunctionBuilder::new(*FIXED).add_arg(*argument).finish());
            }
        }
        Ok(envelope.finish())
    }

    /// Remove the compact canonicalization view after Symbolica has chosen
    /// dummy labels. Rebuilding the physical head invokes its native normalizer.
    pub fn restore_canonical_view(expression: &Atom) -> Result<Atom, String> {
        if !expression.contains_symbol(*ENVELOPE) {
            return Ok(expression.clone());
        }
        let mut failure = None;
        let restored = expression.replace_map(|value, _, out| {
            let AtomView::Fun(function) = value else {
                return;
            };
            if function.get_symbol() != *ENVELOPE {
                return;
            }
            match Self::restore_envelope(function) {
                Ok(atom) => **out = atom,
                Err(error) => failure = Some(error),
            }
        });
        match failure {
            Some(error) => Err(error),
            None => Ok(restored),
        }
    }

    fn restore_envelope(function: FunView<'_>) -> Result<Atom, String> {
        let mut encoded = function.iter();
        let Some(AtomView::Fun(owner)) = encoded.next() else {
            return Err("invalid slot symmetry canonicalization owner".into());
        };
        if owner.get_symbol() != *FIXED || owner.get_nargs() != 1 {
            return Err("invalid slot symmetry canonicalization owner".into());
        }
        let Some(AtomView::Var(owner)) = owner.iter().next() else {
            return Err("invalid slot symmetry canonicalization owner".into());
        };
        let owner = owner.get_symbol();
        let declaration = Self::from_symbol(owner)?
            .ok_or("canonicalization owner lacks slot symmetry metadata")?;
        let count = function.get_nargs() - 1
            + declaration
                .groups()
                .map(|(_, group)| group.len() - 1)
                .sum::<usize>();
        if declaration
            .groups()
            .any(|(_, group)| group.iter().any(|&position| position >= count))
        {
            return Err("canonical slot group lies outside the reconstructed arguments".into());
        }
        let mut args: Vec<Option<Atom>> = vec![None; count];
        for position in 0..count {
            if let Some((kind, group)) = declaration
                .groups()
                .find(|(_, group)| group.contains(&position))
            {
                if position != group[0] {
                    continue;
                }
                let Some(AtomView::Fun(value)) = encoded.next() else {
                    return Err("invalid canonical slot group".into());
                };
                if value.get_symbol() != kind.symbol() || value.get_nargs() != group.len() {
                    return Err("invalid canonical slot group signature".into());
                }
                for (&target, value) in group.iter().zip(value.iter()) {
                    args[target] = Some(value.to_owned());
                }
            } else {
                let Some(AtomView::Fun(value)) = encoded.next() else {
                    return Err("invalid canonical fixed argument".into());
                };
                if value.get_symbol() != *FIXED || value.get_nargs() != 1 {
                    return Err("invalid canonical fixed argument signature".into());
                }
                args[position] = value.iter().next().map(|value| value.to_owned());
            }
        }
        args.into_iter()
            .try_fold(FunctionBuilder::new(owner), |builder, arg| {
                arg.map(|arg| builder.add_arg(arg))
                    .ok_or_else(|| "missing canonical tensor argument".to_string())
            })
            .map(FunctionBuilder::finish)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        network::parsing::AtomStructureExt,
        structure::{
            OrderedStructure,
            abstract_index::AbstractIndex,
            representation::{Minkowski, RepName},
            slot::IsAbstractSlot,
        },
    };
    use symbolica::{function, wrap_symbol};

    fn slot(index: usize) -> Atom {
        LibraryRep::from(Minkowski {})
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    }
    fn declaration() -> TensorSlotSymmetry {
        TensorSlotSymmetry::new([(FactorProjector::Antisymmetric, vec![1, 3])]).unwrap()
    }
    fn head(name: &str) -> Symbol {
        declaration()
            .register_symbol(
                NamespacedSymbol::parse(name),
                vec![SymbolAttribute::Real],
                vec!["test::owner".into()],
                UserData::String("kept".into()),
            )
            .unwrap()
    }

    #[test]
    fn selected_arguments_keep_metadata_and_same_representation_spectators() {
        let head = head("slot_symmetry_tests::spectators");
        let forward = function!(head, 5, slot(1), slot(7), slot(2));
        let reverse = function!(head, 5, slot(2), slot(7), slot(1));
        assert_eq!(forward, -reverse);
        assert_ne!(forward, -function!(head, 5, slot(7), slot(1), slot(2)));
        assert_eq!(function!(head, 5, slot(1), slot(7), slot(1)), Atom::Zero);
        assert_eq!(head.get_data(), &UserData::String("kept".into()));
        let function = forward.as_fun_view().unwrap();
        let encoded = TensorSlotSymmetry::canonical_view(function).unwrap();
        assert_eq!(
            TensorSlotSymmetry::restore_canonical_view(&encoded).unwrap(),
            forward
        );
        assert!(!encoded.as_view().is_zero());
    }

    #[test]
    fn declarations_reject_overlapping_and_unsupported_groups() {
        assert!(
            TensorSlotSymmetry::new([
                (FactorProjector::Antisymmetric, vec![0, 1]),
                (FactorProjector::Symmetric, vec![1, 2])
            ])
            .unwrap_err()
            .contains("disjoint")
        );
        assert!(TensorSlotSymmetry::new([(FactorProjector::Cyclic, vec![0, 1])]).is_err());
        assert!(TensorSlotSymmetry::new([(FactorProjector::Symmetric, vec![0])]).is_err());
        assert!(
            declaration()
                .register_symbol(
                    wrap_symbol!("slot_symmetry_tests::scalar_attribute"),
                    vec![SymbolAttribute::Scalar],
                    vec![],
                    UserData::None,
                )
                .is_err()
        );
    }

    #[test]
    fn variadic_derivative_spectators_keep_one_immutable_declaration() {
        let head = TensorSlotSymmetry::new([(FactorProjector::Antisymmetric, vec![1, 2])])
            .unwrap()
            .register_symbol(
                wrap_symbol!("slot_symmetry_tests::derivative"),
                vec![],
                vec![],
                UserData::None,
            )
            .unwrap();
        for count in [0, 1, 3, 17] {
            let build = |left, right| {
                FunctionBuilder::new(head)
                    .add_arg(Atom::num(count as i64))
                    .add_arg(slot(left))
                    .add_arg(slot(right))
                    .add_args((0..count).map(|index| slot(10 + index)))
                    .finish()
            };
            let value = build(1, 2);
            assert_eq!(value, -build(2, 1));
            value
                .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
                .unwrap();
            let encoded = TensorSlotSymmetry::canonical_view(value.as_fun_view().unwrap()).unwrap();
            assert_eq!(
                TensorSlotSymmetry::restore_canonical_view(&encoded).unwrap(),
                value
            );
        }
    }

    #[test]
    fn component_shadow_requires_an_explicit_realization() {
        use crate::network::parsing::ShadowedStructure;
        let value = function!(
            head("slot_symmetry_tests::component_shadow"),
            5,
            slot(1),
            slot(7),
            slot(2)
        );
        let error = value
            .infer_structure::<ShadowedStructure<AbstractIndex>>()
            .unwrap_err();
        assert!(error.to_string().contains("component-shadow realization"));
    }

    #[test]
    fn native_admission_requires_callback_and_compatible_slot_signatures() {
        let imported = SymbolBuilder::new(wrap_symbol!("slot_symmetry_tests::unrestored"))
            .with_tags([SPENSO_TAG.tensor.clone(), declaration().tag()])
            .build()
            .unwrap();
        let unrestored = function!(imported, 5, slot(1), slot(7), slot(2));
        let error = unrestored
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap_err()
            .to_string();
        assert!(error.contains("before importing atoms"), "{error}");
        let head = head("slot_symmetry_tests::bad_arguments");
        let bad = function!(
            head,
            5,
            slot(1),
            slot(7),
            LibraryRep::from(Minkowski {})
                .new_rep(3)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(2))
                .to_atom()
        );
        let error = bad
            .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
            .unwrap_err()
            .to_string();
        assert!(error.contains("identical representations"), "{error}");
        let bad = function!(head, 5, slot(1), slot(7));
        assert!(
            bad.infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
                .unwrap_err()
                .to_string()
                .contains("outside")
        );
    }

    #[test]
    fn disjoint_groups_do_not_expand_and_preserve_zero_diagnostics() {
        let declaration = TensorSlotSymmetry::new([
            (FactorProjector::Antisymmetric, vec![0, 2]),
            (FactorProjector::Symmetric, vec![1, 3]),
        ])
        .unwrap();
        let head = declaration
            .register_symbol(
                wrap_symbol!("slot_symmetry_tests::disjoint"),
                vec![],
                vec![],
                UserData::None,
            )
            .unwrap();
        let value = function!(head, slot(1), slot(3), slot(2), slot(4));
        assert_eq!(value, -function!(head, slot(2), slot(4), slot(1), slot(3)));
        assert_eq!(value, function!(head, slot(1), slot(4), slot(2), slot(3)));
        let invalid = function!(head, slot(1), slot(3), slot(1), Atom::num(6));
        assert!(!invalid.is_zero());
        assert!(
            invalid
                .infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>()
                .is_err()
        );
    }

    #[test]
    fn fresh_process_archive_child() {
        let Ok(path) = std::env::var("SPENSO_SLOT_SYMMETRY_ARCHIVE") else {
            return;
        };
        let bytes = std::fs::read(path).unwrap();
        let (tags, payload): (Vec<String>, Vec<u8>) = serde_json::from_slice(&bytes).unwrap();
        let mode = std::env::var("SPENSO_SLOT_SYMMETRY_RESTORE").unwrap();
        let restore = mode != "no";
        if mode == "rename" {
            SymbolBuilder::new(wrap_symbol!("slot_symmetry_tests::archive"))
                .with_attributes(vec![SymbolAttribute::Scalar])
                .build()
                .unwrap();
        }
        if restore {
            TensorSlotSymmetry::restore_symbol(
                NamespacedSymbol::parse(if mode == "rename" {
                    "slot_symmetry_tests::archive_renamed"
                } else {
                    "slot_symmetry_tests::archive"
                }),
                vec![SymbolAttribute::Real],
                tags,
                UserData::String("kept".into()),
            )
            .unwrap();
        }
        let imported = if mode == "rename" {
            Atom::import(
                &mut payload.as_slice(),
                Some(Box::new(|name| {
                    assert_eq!(name, "slot_symmetry_tests::archive");
                    "slot_symmetry_tests::archive_renamed".into()
                })),
            )
        } else {
            Atom::import(&mut payload.as_slice(), None)
        }
        .unwrap();
        let result = imported.infer_structure::<OrderedStructure<LibraryRep, AbstractIndex>>();
        if restore {
            result.unwrap();
            let symbol = imported.as_fun_view().unwrap().get_symbol();
            if mode == "rename" {
                assert!(symbol.get_name().ends_with("archive_renamed"));
            }
            assert_eq!(imported, -function!(symbol, 5, slot(2), slot(7), slot(1)));
            assert_eq!(function!(symbol, 5, slot(1), slot(7), slot(1)), Atom::Zero);
        } else {
            assert!(
                result
                    .unwrap_err()
                    .to_string()
                    .contains("before importing atoms")
            );
        }
    }

    #[test]
    fn fresh_process_archive_restores_metadata_before_import() {
        let head = head("slot_symmetry_tests::archive");
        let value = function!(head, 5, slot(1), slot(7), slot(2));
        let mut payload = Vec::new();
        value.export(&mut payload).unwrap();
        let path =
            std::env::temp_dir().join(format!("spenso-slot-symmetry-{}.json", std::process::id()));
        std::fs::write(
            &path,
            serde_json::to_vec(&(head.get_tags(), payload)).unwrap(),
        )
        .unwrap();
        for restore in ["no", "yes", "rename"] {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "network::tensor_symmetry::tests::fresh_process_archive_child",
                    "--nocapture",
                ])
                .env("SPENSO_SLOT_SYMMETRY_ARCHIVE", &path)
                .env("SPENSO_SLOT_SYMMETRY_RESTORE", restore)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
        }
        std::fs::remove_file(path).unwrap();
    }
}
