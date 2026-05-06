use symbolica::atom::{AtomCore, FunctionBuilder, Symbol};
use symbolica::domains::rational::Rational;

use super::*;

use library::Library;

use crate::network::library::DummyLibrary;
use crate::network::library::panicing::ErroringLibrary;
use crate::network::tags::SPENSO_TAG;

use crate::network::library::symbolic::ETS;
use crate::structure::abstract_index::{AIND_SYMBOLS, AbstractIndex};
use crate::structure::representation::{RepName, Representation};
// use crate::shadowing::Concretize;
use crate::structure::slot::{DualSlotTo, DummyAind, ParseableAind, Slot, SlotError};
use crate::structure::{
    NamedStructure, OrderedStructure, PermutedStructure, StructureContract, StructureError,
    TensorShell,
};
use crate::tensors::parametric::ParamTensor;

use std::fmt::Display;

use store::TensorScalarStore;
// use log::trace;

use symbolica::atom::{AddView, Atom, AtomView, MulView, PowView, representation::FunView};
use symbolica::id::MatchSettings;

use crate::structure::{HasStructure, TensorStructure};

use crate::{shadowing::Concretize, structure::HasName, structure::representation::LibraryRep};

pub type ShadowedStructure<Aind> = NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>;

impl<Aind: ParseableAind + AbsInd> Parse for ShadowedStructure<Aind> {
    fn parse(value: AtomView) -> Result<PermutedStructure<Self>, StructureError> {
        if let AtomView::Fun(f) = value {
            f.try_into()
        } else {
            Err(StructureError::ParsingError(value.to_plain_string()))
        }
    }
}

impl<Aind: ParseableAind + AbsInd> OrderedStructure<LibraryRep, Aind> {
    pub fn from_atomcore_unchecked<A: AtomCore>(value: A) -> Result<Self, StructureError> {
        // println!("HII");
        let view = value.as_atom_view();
        match view {
            AtomView::Add(a) => {
                // println!("Add{}", a.as_view());
                Self::from_atomcore_unchecked(a.iter().next().unwrap())
            }
            AtomView::Pow(p) => {
                // println!("Pow{}", p.as_view());
                let (base, exp) = p.get_base_exp();

                let base_strct = Self::from_atomcore_unchecked(base)?;

                if base_strct.is_scalar() {
                    Ok(base_strct)
                } else if base_strct.is_fully_self_dual()
                    && let Ok(r) = Rational::try_from(exp)
                {
                    if r.numerator() % 2 == 0 {
                        Ok(OrderedStructure::empty())
                    } else if r.denominator() == 1 {
                        Ok(base_strct)
                    } else {
                        Err(StructureError::ParsingError(format!(
                            "Invalid power of tensor {}",
                            view
                        )))
                    }
                } else {
                    Err(StructureError::ParsingError(format!(
                        "Invalid power of tensor {}",
                        view
                    )))
                }
            }
            AtomView::Mul(f) => {
                // println!("Mul{}", f.as_view());
                let mut structure = OrderedStructure::empty();

                for i in f.into_iter() {
                    match Self::from_atomcore_unchecked(i) {
                        Ok(s) => structure = structure.merge(&s)?.0,
                        // Err(StructureError::EmptyStructure(SlotError::EmptyStructure)) => {
                        //     continue;
                        // }
                        Err(_) => {
                            continue;
                        }
                    }
                }

                if structure.is_scalar() {
                    return Err(StructureError::EmptyStructure(SlotError::EmptyStructure));
                }
                Ok(structure)
            }
            AtomView::Fun(f) => match f.get_symbol() {
                s if s == AIND_SYMBOLS.aind => {
                    // println!("Fun{}", f.as_view());
                    let mut structure: Vec<Slot<LibraryRep, Aind>> = vec![];

                    for arg in f.iter() {
                        structure.push(arg.try_into()?);
                    }

                    let o = OrderedStructure::new(structure);

                    Ok(o.structure)
                }
                _ => {
                    //     println!("Fun{}", f.as_view());
                    let mut args = vec![];
                    let mut slots = vec![];
                    let mut is_structure: Option<StructureError> =
                        Some(SlotError::EmptyStructure.into());

                    for arg in f.iter() {
                        let slot: Result<Slot<LibraryRep, Aind>, _> = arg.try_into();

                        match slot {
                            Ok(slot) => {
                                is_structure = None;
                                slots.push(slot);
                            }
                            Err(e) => {
                                if let AtomView::Fun(f) = arg
                                    && f.get_symbol() == AIND_SYMBOLS.aind
                                {
                                    let internal_s = Self::from_atomcore_unchecked(f.as_view());

                                    if let Ok(s) = internal_s {
                                        slots.extend(s.structure);
                                        is_structure = None;
                                        continue;
                                    }
                                }
                                args.push(arg.to_owned());
                            }
                        }
                    }

                    if let Some(e) = is_structure {
                        Err(e)
                    } else {
                        Ok(OrderedStructure::new(slots).structure)
                    }
                }
            },
            _ => Err(StructureError::EmptyStructure(SlotError::EmptyStructure)),
        }
    }
}

impl<Aind: ParseableAind + AbsInd> TryFrom<Atom> for PermutedStructure<ShadowedStructure<Aind>> {
    type Error = StructureError;
    fn try_from(value: Atom) -> Result<Self, Self::Error> {
        ShadowedStructure::<Aind>::parse(value.as_view())
    }
}

impl<'a, Aind: ParseableAind + AbsInd> TryFrom<&'a Atom>
    for PermutedStructure<ShadowedStructure<Aind>>
{
    type Error = StructureError;
    fn try_from(value: &'a Atom) -> Result<Self, Self::Error> {
        ShadowedStructure::<Aind>::parse(value.as_view())
    }
}

impl<'a, Aind: ParseableAind + AbsInd> TryFrom<FunView<'a>>
    for PermutedStructure<ShadowedStructure<Aind>>
{
    type Error = StructureError;
    fn try_from(value: FunView<'a>) -> Result<Self, Self::Error> {
        match value.get_symbol() {
            s if s == AIND_SYMBOLS.aind => {
                let mut structure: Vec<Slot<LibraryRep, Aind>> = vec![];

                for arg in value.iter() {
                    structure.push(arg.try_into()?);
                }

                let o = OrderedStructure::new(structure);
                let p = o.map_structure(|s| s.into());

                Ok(p)
            }
            name => {
                let mut args = vec![];
                let mut slots = vec![];
                let mut is_structure: Option<StructureError> =
                    Some(SlotError::EmptyStructure.into());

                for arg in value.iter() {
                    let slot: Result<Slot<LibraryRep, Aind>, _> = arg.try_into();

                    match slot {
                        Ok(slot) => {
                            is_structure = None;
                            slots.push(slot);
                        }
                        Err(e) => {
                            if let AtomView::Fun(f) = arg
                                && f.get_symbol() == AIND_SYMBOLS.aind
                            {
                                let internal_s = Self::try_from(f);

                                if let Ok(s) = internal_s {
                                    let p = s.index_permutation;
                                    let mut v = s.structure.structure.structure;
                                    p.apply_slice_in_place_inv(&mut v); //undo sorting
                                    let p = s.rep_permutation;
                                    p.apply_slice_in_place_inv(&mut v); //undo sorting
                                    slots.extend(v);
                                    is_structure = None;
                                    continue;
                                }
                            }
                            is_structure = Some(e.into());
                            args.push(arg.to_owned());
                        }
                    }
                }

                if let Some(e) = is_structure {
                    Err(e)
                } else {
                    let mut structure: PermutedStructure<ShadowedStructure<Aind>> =
                        OrderedStructure::new(slots).map_structure(Into::into);
                    structure.structure.set_name(name);
                    if !args.is_empty() {
                        structure.structure.additional_args = Some(args);
                    }
                    Ok(structure)
                }
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct ParseSettings {
    pub precontract_scalars: bool,
    pub take_first_term_from_sum: bool,
    pub depth_limit: Option<usize>,
    pub depth_is_product_depth: bool,
    pub parse_inner_products: bool,
    pub parse_composite_scalars_as_tensors: bool,
}

impl Default for ParseSettings {
    fn default() -> Self {
        Self {
            precontract_scalars: true,
            take_first_term_from_sum: false,
            depth_limit: None,
            depth_is_product_depth: true,
            parse_inner_products: true,
            parse_composite_scalars_as_tensors: false,
        }
    }
}

#[derive(Clone, Debug, Copy)]
pub struct ParseState {
    depth: usize,
    // func_depth: usize,
}

#[allow(clippy::derivable_impls)]
impl Default for ParseState {
    fn default() -> Self {
        Self { depth: 0 }
    }
}

pub trait Parse: Sized {
    fn parse(view: AtomView) -> Result<PermutedStructure<Self>, StructureError>;

    fn parse_with_settings(
        view: AtomView,
        _settings: &ParseSettings,
    ) -> Result<PermutedStructure<Self>, StructureError> {
        Self::parse(view)
    }
}

struct ParserDummies {
    next: usize,
}

impl ParserDummies {
    fn new() -> Self {
        Self { next: 1_000_000 }
    }

    fn next(&mut self) -> AbstractIndex {
        let index = self.next;
        self.next += 1;
        AbstractIndex::new_dummy_at(index)
    }

    fn slot(&mut self, rep: &Representation<LibraryRep>) -> Slot<LibraryRep, AbstractIndex> {
        rep.slot(self.next())
    }
}

struct SchoonschipMaterialization {
    replacement: Atom,
    factors: Vec<Atom>,
    compact_self: bool,
}

fn compact_rep_pattern_match(arg: AtomView<'_>) -> Option<Representation<LibraryRep>> {
    let rep_pattern = Atom::var(SPENSO_TAG.rep_).to_pattern();
    let settings = MatchSettings {
        level_range: (0, Some(0)),
        partial: false,
        ..Default::default()
    };
    let mut matches = arg.pattern_match(&rep_pattern, None, Some(&settings));
    let matched = matches.next_detailed()?;
    let rep = rep_pattern.replace_wildcards_with_matches(matched.match_stack);
    Representation::<LibraryRep>::try_from(rep.as_view()).ok()
}

fn compact_tensor_rep_arg(value: FunView<'_>) -> Option<(usize, Representation<LibraryRep>)> {
    if value.get_symbol() == ETS.metric || value.get_symbol() == SPENSO_TAG.dot {
        return None;
    }

    if Representation::<LibraryRep>::try_from(value.as_view()).is_ok() {
        return None;
    }

    let args = value.iter().collect::<Vec<_>>();
    if args
        .iter()
        .any(|arg| Slot::<LibraryRep, AbstractIndex>::try_from(*arg).is_ok())
    {
        return None;
    }

    let rep_args = args
        .iter()
        .enumerate()
        .filter_map(|(position, arg)| compact_rep_pattern_match(*arg).map(|rep| (position, rep)))
        .collect::<Vec<_>>();

    let [(position, rep)] = rep_args.as_slice() else {
        return None;
    };

    Some((*position, *rep))
}

fn compact_vector_rep(value: AtomView<'_>) -> Option<Representation<LibraryRep>> {
    match value {
        AtomView::Fun(fun) => compact_tensor_rep_arg(fun).map(|(_, rep)| rep),
        AtomView::Add(add) => {
            let mut reps = add.iter().map(compact_vector_rep);
            let rep = reps.next()??;
            reps.all(|candidate| candidate == Some(rep)).then_some(rep)
        }
        _ => None,
    }
}

fn materialize_compact_vector_with_slot(
    value: AtomView<'_>,
    rep: &Representation<LibraryRep>,
    slot: &Atom,
) -> Option<Atom> {
    match value {
        AtomView::Fun(fun) => {
            let (position, matched_rep) = compact_tensor_rep_arg(fun)?;
            if matched_rep != *rep {
                return None;
            }

            let mut tensor = FunctionBuilder::new(fun.get_symbol());
            for (arg_position, arg) in fun.iter().enumerate() {
                if arg_position == position {
                    tensor = tensor.add_arg(slot);
                } else {
                    tensor = tensor.add_arg(arg);
                }
            }
            Some(tensor.finish())
        }
        AtomView::Add(add) => {
            let mut terms = add
                .iter()
                .map(|term| materialize_compact_vector_with_slot(term, rep, slot));
            let first = terms.next()??;
            let rest = terms.collect::<Option<Vec<_>>>()?;
            Some(rest.into_iter().fold(first, |sum, term| sum + term))
        }
        _ => None,
    }
}

fn compact_scalar_product(
    value: FunView<'_>,
    dummies: &mut ParserDummies,
) -> Option<SchoonschipMaterialization> {
    if value.get_symbol() != ETS.metric && value.get_symbol() != SPENSO_TAG.dot {
        return None;
    }

    let args = value.iter().collect::<Vec<_>>();
    let [lhs, rhs] = args.as_slice() else {
        return None;
    };

    let rep = compact_vector_rep(*lhs)?;
    if rep != compact_vector_rep(*rhs)? || !rep.rep.is_self_dual() {
        return None;
    }

    let slot = dummies.slot(&rep).to_atom();
    let lhs = materialize_compact_vector_with_slot(*lhs, &rep, &slot)?;
    let rhs = materialize_compact_vector_with_slot(*rhs, &rep, &slot)?;

    Some(SchoonschipMaterialization {
        replacement: Atom::num(1),
        factors: vec![lhs, rhs],
        compact_self: true,
    })
}

/// Turns compact Schoonschip notation such as `p(mink(4))` into an indexed
/// tensor factor and returns the fresh slot that should replace it in context.
fn compact_schoonschip_tensor(
    value: FunView<'_>,
    dummies: &mut ParserDummies,
) -> Option<SchoonschipMaterialization> {
    let (position, rep) = compact_tensor_rep_arg(value)?;
    let args = value.iter().collect::<Vec<_>>();

    let slot = dummies.slot(&rep);
    let slot_atom = slot.to_atom();
    let mut tensor = FunctionBuilder::new(value.get_symbol());
    for (arg_position, arg) in args.into_iter().enumerate() {
        if arg_position == position {
            tensor = tensor.add_arg(&slot_atom);
        } else {
            tensor = tensor.add_arg(arg);
        }
    }

    Some(SchoonschipMaterialization {
        replacement: slot_atom,
        factors: vec![tensor.finish()],
        compact_self: true,
    })
}

fn compact_dot_as_metric(value: AtomView<'_>) -> Option<Atom> {
    let AtomView::Fun(fun) = value else {
        return None;
    };

    if fun.get_symbol() != SPENSO_TAG.dot || fun.get_nargs() != 2 {
        return None;
    }

    let args = fun.iter().collect::<Vec<_>>();
    if args
        .iter()
        .all(|arg| Slot::<LibraryRep, AbstractIndex>::try_from(*arg).is_ok())
    {
        Some(
            FunctionBuilder::new(ETS.metric)
                .add_arg(args[0])
                .add_arg(args[1])
                .finish(),
        )
    } else {
        None
    }
}

fn materialize_schoonschip_arg(
    value: AtomView<'_>,
    dummies: &mut ParserDummies,
) -> Option<SchoonschipMaterialization> {
    let AtomView::Fun(fun) = value else {
        return None;
    };

    if let Some(materialized) = compact_scalar_product(fun, dummies) {
        return Some(materialized);
    }

    if let Some(materialized) = compact_schoonschip_tensor(fun, dummies) {
        return Some(materialized);
    }

    let mut changed = false;
    let mut factors = Vec::new();
    let mut rebuilt = FunctionBuilder::new(fun.get_symbol());

    for arg in fun.iter() {
        if let Some(materialized) = materialize_schoonschip_arg(arg, dummies) {
            changed = true;
            rebuilt = rebuilt.add_arg(&materialized.replacement);
            factors.extend(materialized.factors);
        } else {
            rebuilt = rebuilt.add_arg(arg);
        }
    }

    changed.then(|| {
        let replacement = rebuilt.finish();
        SchoonschipMaterialization {
            replacement: compact_dot_as_metric(replacement.as_view()).unwrap_or(replacement),
            factors,
            compact_self: false,
        }
    })
}

fn materialize_schoonschip(value: AtomView<'_>) -> Option<Atom> {
    let mut dummies = ParserDummies::new();
    let materialized = materialize_schoonschip_arg(value, &mut dummies)?;
    let mut expression = if materialized.compact_self {
        Atom::num(1)
    } else {
        materialized.replacement
    };

    for factor in materialized.factors {
        expression *= factor;
    }

    Some(expression)
}

fn replace_chain_placeholders(value: AtomView<'_>, chain_in: &Atom, chain_out: &Atom) -> Atom {
    match value {
        AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in => chain_in.clone(),
        AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out => chain_out.clone(),
        AtomView::Fun(fun) => {
            let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
            for arg in fun.iter() {
                rebuilt = rebuilt.add_arg(replace_chain_placeholders(arg, chain_in, chain_out));
            }
            rebuilt.finish()
        }
        AtomView::Add(add) => add.iter().fold(Atom::Zero, |sum, term| {
            sum + replace_chain_placeholders(term, chain_in, chain_out)
        }),
        AtomView::Mul(mul) => mul.iter().fold(Atom::num(1), |product, factor| {
            product * replace_chain_placeholders(factor, chain_in, chain_out)
        }),
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            replace_chain_placeholders(base, chain_in, chain_out).pow(exponent.to_owned())
        }
        _ => value.to_owned(),
    }
}

fn parser_product(factors: impl IntoIterator<Item = Atom>) -> Atom {
    factors
        .into_iter()
        .fold(Atom::num(1), |product, factor| product * factor)
}

impl<
    'a,
    Sc,
    T: HasStructure + TensorStructure,
    K: Clone + Display + Debug,
    // FK: Clone + Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = Sc> + Clone,
    Aind: AbsInd,
> Network<Str, K, Symbol, Aind>
where
    Sc: for<'r> TryFrom<AtomView<'r>> + Clone,
    TensorNetworkError<K, Symbol>: for<'r> From<<Sc as TryFrom<AtomView<'r>>>::Error>,
{
    #[allow(clippy::result_large_err)]
    pub fn try_from_view<S, Lib: Library<S, Key = K>>(
        value: AtomView<'a>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let state = ParseState::default();
        Self::try_from_view_impl(value, state, library, settings)
    }

    #[allow(clippy::result_large_err)]
    fn try_from_view_impl<S, Lib: Library<S, Key = K>>(
        value: AtomView<'a>,
        state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        match value {
            AtomView::Mul(m) => Self::try_from_mul(m, state, library, settings),
            AtomView::Fun(f) => Self::try_from_fun(f, state, library, settings),
            AtomView::Add(a) => Self::try_from_add(a, state, library, settings),
            AtomView::Pow(p) => Self::try_from_pow(p, state, library, settings),
            a => Ok(Network::from_scalar(a.try_into()?)),
        }
    }

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    fn as_leaf<S>(
        value: AtomView<'a>,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let s: Result<PermutedStructure<S>, _> = S::parse_with_settings(value, settings);

        // println!("Looking at :{}", value);
        return if let Ok(s) = s {
            Ok(Self::from_tensor(
                s.structure
                    .to_shell()
                    .concretize(Some(s.index_permutation.inverse())),
            ))
        } else {
            Ok(Self::from_scalar(value.try_into()?))
        };
    }

    #[allow(clippy::result_large_err)]
    fn try_from_mul<S, Lib: Library<S, Key = K>>(
        value: MulView<'a>,
        mut state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        // println!("Mul");
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            // println!("Mul leaf");
            return Self::as_leaf::<S>(value.as_view(), settings);
        }

        state.depth += 1;
        // println!("{} for mul {}", state.depth, value.as_view());
        let mut iter = value.iter();
        let first_atom = iter.next().unwrap();
        let first = Self::try_from_view_impl(first_atom, state, library, settings)?;

        // state

        if settings.precontract_scalars {
            let mut scalars = Atom::num(1);

            let rest: Result<Vec<_>, _> = iter
                .filter_map(
                    |a| match Self::try_from_view_impl(a, state, library, settings) {
                        Ok(n) => {
                            if let NetworkState::PureScalar = n.state {
                                scalars *= a;
                                None
                            } else {
                                Some(Ok(n))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    },
                )
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = first.state {
                scalars *= first_atom;
            } else {
                res.push(first);
            }

            if res.is_empty() {
                Ok(Self::from_scalar(value.as_view().try_into()?))
            } else {
                let s = if scalars != Atom::num(1) {
                    Self::from_scalar(scalars.as_view().try_into()?)
                } else {
                    res.pop().unwrap()
                };

                Ok(s.n_mul(res))
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(|a| Self::try_from_view_impl(a, state, library, settings))
                .collect();

            Ok(first.n_mul(rest?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_fun<S, Lib: Library<S, Key = K>>(
        value: FunView<'a>,
        state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        // <PermutedStructure<S>>::Error: Debug,
    {
        let symbol = value.get_symbol();

        if symbol == SPENSO_TAG.bracket {
            let mut n_muls = value
                .iter()
                .map(|a| Self::try_from_view_impl(a, state, library, settings))
                .collect::<Result<Vec<_>, _>>()?;

            Ok(n_muls.pop().unwrap().n_mul(n_muls))
        } else if symbol == SPENSO_TAG.chain {
            Self::parse_chain(value, state, library, settings)
        } else if symbol == SPENSO_TAG.trace {
            Self::parse_trace(value, state, library, settings)
        } else if symbol == SPENSO_TAG.pure_scalar {
            if value.get_nargs() != 1 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            Ok(Self::from_scalar(value.iter().next().unwrap().try_into()?))
        } else if let Some(materialized) = materialize_schoonschip(value.as_view()) {
            Self::try_from_view_impl(materialized.as_view(), state, library, settings)
        } else if symbol.has_tag(&SPENSO_TAG.tag) {
            if value.get_nargs() != 1 {
                return Err(TensorNetworkError::TooManyArgsFunction(
                    value.as_view().to_plain_string(),
                ));
            }

            let inner = value.iter().next().unwrap();
            let inner_tensor = Self::try_from_view_impl(inner, state, library, settings)?;

            Ok(inner_tensor.fun(symbol))
        } else {
            let s: Result<PermutedStructure<S>, _> =
                S::parse_with_settings(value.as_view(), settings);

            if let Ok(s) = s {
                // println!("Perm:{}", s.index_permutation);
                // let s = s;
                match library.key_for_structure(&s) {
                    Ok(key) => {
                        // println!("Adding lib");
                        // let t = library.get(&key).unwrap();
                        Ok(Self::library_tensor(
                            &s.structure,
                            PermutedStructure {
                                structure: key,
                                rep_permutation: s.rep_permutation,
                                index_permutation: s.index_permutation,
                            },
                        ))
                    }
                    Err(_) => Ok(Self::from_tensor(
                        s.structure
                            .to_shell()
                            .concretize(Some(s.index_permutation.inverse())),
                    )),
                }
            } else {
                Ok(Self::from_scalar(value.as_view().try_into()?))
            }
        }
    }

    #[allow(clippy::result_large_err)]
    fn parse_chain<S, Lib: Library<S, Key = K>>(
        value: FunView<'a>,
        state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let args = value.iter().collect::<Vec<_>>();
        if args.len() < 2 {
            return Err(TensorNetworkError::TooManyArgsFunction(
                value.as_view().to_plain_string(),
            ));
        }

        let start = args[0].to_owned();
        let end = args[1].to_owned();
        let start_slot = Slot::<LibraryRep, AbstractIndex>::try_from(args[0])
            .map_err(|err| eyre!("invalid chain start `{}`: {err}", args[0]))?;

        let factors = &args[2..];
        if factors.is_empty() {
            let metric = FunctionBuilder::new(ETS.metric)
                .add_arg(&start)
                .add_arg(&end)
                .finish();
            return Self::try_from_view_impl(metric.as_view(), state, library, settings);
        }

        // The chain head is inert at expression level; parsing materializes it
        // into ordinary indexed tensor factors with fresh contracted links.
        let mut dummies = ParserDummies::new();
        let mut left = start;
        let mut materialized_factors = Vec::with_capacity(factors.len());
        for (position, factor) in factors.iter().enumerate() {
            let (right, next_left) = if position + 1 == factors.len() {
                (end.clone(), None)
            } else {
                let link = start_slot.reindex(dummies.next());
                (link.dual().to_atom(), Some(link.to_atom()))
            };

            let factor = replace_chain_placeholders(*factor, &left, &right);
            materialized_factors.push(materialize_schoonschip(factor.as_view()).unwrap_or(factor));
            if let Some(next_left) = next_left {
                left = next_left;
            }
        }

        let expression = parser_product(materialized_factors);
        Self::try_from_view_impl(expression.as_view(), state, library, settings)
    }

    #[allow(clippy::result_large_err)]
    fn parse_trace<S, Lib: Library<S, Key = K>>(
        value: FunView<'a>,
        state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let args = value.iter().collect::<Vec<_>>();
        let Some(rep_view) = args.first() else {
            return Err(TensorNetworkError::TooManyArgsFunction(
                value.as_view().to_plain_string(),
            ));
        };

        let rep = Representation::<LibraryRep>::try_from(*rep_view)
            .map_err(|err| eyre!("invalid trace representation `{rep_view}`: {err}"))?;
        let factors = &args[1..];

        if factors.is_empty() {
            return Self::try_from_view_impl(
                rep.dim.to_symbolic().as_view(),
                state,
                library,
                settings,
            );
        }

        // A trace is parsed as a closed chain over fresh slots of the declared
        // representation, so the factors themselves can stay endpoint-free.
        let mut dummies = ParserDummies::new();
        if let [factor] = factors {
            let left = dummies.slot(&rep);
            let right = dummies.slot(&rep);
            let bridge = dummies.slot(&rep);
            let materialized_factor =
                replace_chain_placeholders(*factor, &left.to_atom(), &right.dual().to_atom());
            let materialized_factor = materialize_schoonschip(materialized_factor.as_view())
                .unwrap_or(materialized_factor);
            let left_closure = FunctionBuilder::new(ETS.metric)
                .add_arg(bridge.to_atom())
                .add_arg(left.dual().to_atom())
                .finish();
            let right_closure = FunctionBuilder::new(ETS.metric)
                .add_arg(right.to_atom())
                .add_arg(bridge.dual().to_atom())
                .finish();
            let factor_net =
                Self::try_from_view_impl(materialized_factor.as_view(), state, library, settings)?;
            let left_net =
                Self::try_from_view_impl(left_closure.as_view(), state, library, settings)?;
            let right_net =
                Self::try_from_view_impl(right_closure.as_view(), state, library, settings)?;
            return Ok(factor_net.n_mul([right_net, left_net]));
        }

        let links = (0..factors.len())
            .map(|_| dummies.slot(&rep))
            .collect::<Vec<_>>();

        let materialized_factors = factors
            .iter()
            .enumerate()
            .map(|(position, factor)| {
                let left = links[position].to_atom();
                let right = links[(position + 1) % links.len()].dual().to_atom();
                let factor = replace_chain_placeholders(*factor, &left, &right);
                materialize_schoonschip(factor.as_view()).unwrap_or(factor)
            })
            .collect::<Vec<_>>();

        let expression = parser_product(materialized_factors);
        Self::try_from_view_impl(expression.as_view(), state, library, settings)
    }

    #[allow(clippy::result_large_err)]
    fn try_from_pow<S, Lib: Library<S, Key = K>>(
        value: PowView<'a>,
        mut state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> std::result::Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            return Self::as_leaf::<S>(value.as_view(), settings);
        }

        if !settings.depth_is_product_depth {
            state.depth += 1;
        }
        let (base, exp) = value.get_base_exp();

        if let Ok(n) = i8::try_from(exp) {
            // println!("base:{base}");
            let base = Self::try_from_view_impl(base, state, library, settings)?;

            // println!("base state {:?}", base.state);
            if settings.precontract_scalars
                && let NetworkState::PureScalar = base.state
            {
                // println!("Pure");
                return Ok(Self::from_scalar(value.as_view().try_into()?));
            }

            if let NetworkState::Tensor = base.state {
                Err(TensorNetworkError::NonSelfDualTensorPower(
                    value.as_view().to_plain_string(),
                ))
            } else if n < 0 {
                // An even power of a self_dual tensor, or scalar is a scalar
                if n % 2 == 0 || base.state.is_scalar() {
                    let out = base.pow(n);
                    // println!("{:?}", out.state);
                    Ok(out)
                } else {
                    Err(TensorNetworkError::NegativeExponentNonScalar(format!(
                        "Atom:{},graph of base: {}, dangling indices: {:?}",
                        value.as_view().to_plain_string(),
                        base.dot(),
                        base.graph.dangling_indices()
                    )))
                }
            } else {
                let out = base.pow(n);
                // println!("{:?}", out.state);
                Ok(out)
            }
        } else {
            Ok(Self::from_scalar(value.as_view().try_into()?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_add<S, Lib: Library<S, Key = K>>(
        value: AddView<'a>,
        mut state: ParseState,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        if let Some(a) = settings.depth_limit
            && a <= state.depth
        {
            return Self::as_leaf::<S>(value.as_view(), settings);
        }

        if !settings.depth_is_product_depth {
            state.depth += 1;
        }
        let mut iter = value.iter();

        let first_atom = iter.next().unwrap();

        let first = Self::try_from_view_impl(first_atom, state, library, settings)?;
        if settings.take_first_term_from_sum {
            Ok(first)
        } else if settings.precontract_scalars {
            let mut scalars = Atom::Zero;

            let rest: Result<Vec<_>, _> = iter
                .filter_map(
                    |a| match Self::try_from_view_impl(a, state, library, settings) {
                        Ok(n) => {
                            if n.state.is_compatible(&first.state) {
                                if let NetworkState::PureScalar = n.state {
                                    scalars += a;
                                    None
                                } else {
                                    Some(Ok(n))
                                }
                            } else {
                                Some(Err(TensorNetworkError::IncompatibleSummand(format!(
                                    "{} is {:?} vs {} is {:?}",
                                    a, n.state, first_atom, first.state
                                ))))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    },
                )
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = first.state {
                scalars += first_atom;
            } else {
                res.push(first);
            }

            if res.is_empty() {
                Ok(Self::from_scalar(value.as_view().try_into()?))
            } else {
                let s = if scalars != Atom::Zero {
                    Self::from_scalar(scalars.as_view().try_into()?)
                } else {
                    res.pop().unwrap()
                };
                Ok(s.n_add(res))
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(
                    |a| match Self::try_from_view_impl(a, state, library, settings) {
                        Ok(n) => {
                            if n.state.is_compatible(&first.state) {
                                Ok(n)
                            } else {
                                Err(TensorNetworkError::IncompatibleSummand(format!(
                                    "{} is {:?} vs {} is {:?}",
                                    a, n.state, first_atom, first.state
                                )))
                            }
                        }
                        Err(e) => Err(e),
                    },
                )
                .collect();

            Ok(first.n_add(rest?))
        }
    }
}
pub type ParamNet<Aind> =
    Network<NetworkStore<ParamTensor<ShadowedStructure<Aind>>, Atom>, DummyKey, Symbol, Aind>;

impl<Aind: AbsInd + DummyAind + ParseableAind + 'static> ParamNet<Aind> {
    pub fn simple_execute(&mut self) {
        let lib = DummyLibrary::<_>::new();

        self.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &ErroringLibrary::new())
            .unwrap();
    }
}
pub trait NetworkParse {
    #[allow(clippy::result_large_err)]
    fn parse_to_atom_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;
}

impl NetworkParse for Atom {
    fn parse_to_atom_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_atom_net::<Aind>(settings)
    }
}

impl NetworkParse for AtomView<'_> {
    fn parse_to_atom_net<Aind: AbsInd + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<ParamTensor<ShadowedStructure<Aind>>>::new();

        ParamNet::<Aind>::try_from_view::<ShadowedStructure<Aind>, _>(*self, &lib, settings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::NetworkState;
    use crate::structure::representation::{Lorentz, Minkowski, RepName};
    use crate::{chain, slot, trace};
    use symbolica::{function, symbol};

    fn mink4() -> Representation<Minkowski> {
        Minkowski {}.new_rep(4)
    }

    fn chain_factor(name: Symbol) -> Atom {
        FunctionBuilder::new(name)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .finish()
    }

    #[test]
    fn parse_chain_materializes_placeholder_links() {
        let rep = mink4();
        let expr = chain!(
            slot!(rep, i),
            slot!(rep, j),
            chain_factor(symbol!("f")),
            chain_factor(symbol!("g")),
        );

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(parsed.state, NetworkState::SelfDualTensor);
        assert_eq!(parsed.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn parse_trace_materializes_closed_links() {
        let rep = mink4();
        let expr = trace!(&rep, chain_factor(symbol!("f")), chain_factor(symbol!("g")));

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_single_factor_trace_closes_dualizable_links() {
        let rep = Lorentz {}.new_rep(4);
        let expr = trace!(&rep, chain_factor(symbol!("f")));

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_empty_trace_unwraps_representation_dimension() {
        let rep = mink4();
        let expr = trace!(&rep);

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(parsed.state, NetworkState::PureScalar);
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_empty_chain_as_endpoint_metric() {
        let rep = mink4();
        let expr = chain!(slot!(rep, i), slot!(rep, j));

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(parsed.state, NetworkState::SelfDualTensor);
        assert_eq!(parsed.graph.dangling_indices().len(), 2);
    }

    #[test]
    fn compact_schoonschip_tensor_uses_direct_rep_pattern_match() {
        let rep = mink4();
        let expr = function!(symbol!("p"), rep.to_symbolic([]));
        let AtomView::Fun(fun) = expr.as_view() else {
            panic!("expected function");
        };

        let materialized = compact_schoonschip_tensor(fun, &mut ParserDummies::new()).unwrap();

        assert!(
            Slot::<LibraryRep, AbstractIndex>::try_from(materialized.replacement.as_view()).is_ok()
        );
        assert_eq!(materialized.factors.len(), 1);
    }

    #[test]
    fn compact_schoonschip_tensor_ignores_nested_rep_pattern_match() {
        let rep = mink4();
        let expr = function!(
            symbol!("p"),
            function!(symbol!("label"), rep.to_symbolic([]))
        );
        let AtomView::Fun(fun) = expr.as_view() else {
            panic!("expected function");
        };

        assert!(compact_schoonschip_tensor(fun, &mut ParserDummies::new()).is_none());
    }

    #[test]
    fn materialize_schoonschip_vectors_parse_as_dot_product() {
        let rep = mink4();
        let expr = function!(symbol!("p"), rep.to_symbolic([]))
            * function!(symbol!("q"), rep.to_symbolic([]));

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn three_argument_metric_inner_product_is_not_parser_syntax() {
        let rep = mink4();
        let expr = function!(
            ETS.metric,
            rep.to_symbolic([]),
            Atom::var(symbol!("p")),
            Atom::var(symbol!("q"))
        );

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(parsed.state, NetworkState::PureScalar);
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_schoonschipped_metric_product() {
        let rep = mink4();
        let expr = function!(
            ETS.metric,
            function!(symbol!("p"), rep.to_symbolic([])),
            function!(symbol!("q"), rep.to_symbolic([]))
        );

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_schoonschipped_metric_sum_product() {
        let rep = mink4();
        let k1 = function!(symbol!("k"), Atom::num(1), rep.to_symbolic([]));
        let k2 = function!(symbol!("k"), Atom::num(2), rep.to_symbolic([]));
        let p = function!(symbol!("p"), Atom::num(3), rep.to_symbolic([]));
        let expr = function!(ETS.metric, k1 + k2, p);

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_schoonschipped_dot_product() {
        let rep = mink4();
        let expr = function!(
            SPENSO_TAG.dot,
            function!(symbol!("p"), rep.to_symbolic([])),
            function!(symbol!("q"), rep.to_symbolic([]))
        );

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_linear_schoonschipped_dot_product() {
        let rep = mink4();
        let p = function!(symbol!("p"), rep.to_symbolic([]));
        let q = function!(symbol!("q"), rep.to_symbolic([]));
        let r = function!(symbol!("r"), rep.to_symbolic([]));
        let expr = function!(SPENSO_TAG.dot, p + q, r);

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert!(parsed.state.is_scalar());
        assert!(parsed.graph.dangling_indices().is_empty());
    }

    #[test]
    fn parse_chain_materializes_schoonschip_factor_arguments() {
        let rep = mink4();
        let compact_vector = function!(symbol!("p"), rep.to_symbolic([]));
        let factor = FunctionBuilder::new(symbol!("f"))
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .add_arg(&compact_vector)
            .finish();
        let expr = chain!(slot!(rep, i), slot!(rep, j), factor);

        let parsed = expr
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();

        assert_eq!(parsed.state, NetworkState::SelfDualTensor);
        assert_eq!(parsed.graph.dangling_indices().len(), 2);
    }
}
