use symbolica::atom::{AtomCore, FunctionBuilder, Symbol};

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
    NamedStructure, OrderedStructure, PermutedStructure, StructureError, TensorShell,
};
use crate::tensors::parametric::ParamTensor;

use std::{cell::Cell, fmt::Display, marker::PhantomData, rc::Rc};

use store::TensorScalarStore;
// use log::trace;

use symbolica::atom::{AddView, Atom, AtomView, MulView, PowView, representation::FunView};
use symbolica::id::MatchSettings;

use crate::structure::{HasStructure, TensorStructure};

use crate::{shadowing::Concretize, structure::HasName, structure::representation::LibraryRep};

pub type ShadowedStructure<Aind> = NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>;

mod structure_inference;
pub use structure_inference::{AtomStructureExt, StructureFromAtom, StructureInferenceMode};

impl<Aind: ParseableAind + AbsInd> Parse for ShadowedStructure<Aind> {
    fn parse(value: AtomView) -> Result<PermutedStructure<Self>, StructureError> {
        if let AtomView::Fun(f) = value {
            f.try_into()
        } else {
            Err(StructureError::ParsingError(value.to_plain_string()))
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
                            if slots.is_empty() {
                                is_structure = Some(e.into());
                            }
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ShorthandParsing {
    /// Expand shorthand notation into an explicit network.
    Expand,
    /// Keep shorthand notation as a leaf and infer its exposed structure.
    Opaque { inference: StructureInferenceMode },
}

impl ShorthandParsing {
    pub fn expands(self) -> bool {
        matches!(self, Self::Expand)
    }

    fn opaque_inference_for(self, symbol: Symbol) -> Option<StructureInferenceMode> {
        match self {
            Self::Expand => None,
            Self::Opaque { inference } if Self::is_shorthand_symbol(symbol) => Some(inference),
            Self::Opaque { .. } => None,
        }
    }

    fn is_shorthand_symbol(symbol: Symbol) -> bool {
        symbol == SPENSO_TAG.chain || symbol == SPENSO_TAG.trace || symbol == SPENSO_TAG.dot
    }
}

#[derive(Clone, Debug)]
pub struct ParseSettings {
    /// Fold factors that parse as pure scalars into one scalar factor while
    /// parsing a product.
    ///
    /// This keeps scalar-only subexpressions out of the tensor graph when they
    /// cannot affect contraction topology. Disable it when the caller needs
    /// every product factor to remain represented separately in the network.
    pub precontract_scalars: bool,

    /// Parse only the first summand of an addition.
    ///
    /// This is meant for structure discovery paths where all summands are
    /// expected to have the same external structure, for example listing
    /// dangling indices without building a full sum network.
    pub take_first_term_from_sum: bool,

    /// Stop recursive parsing once the current parse depth reaches this value.
    ///
    /// At the limit, the current expression is handed to `Parse::parse_with_settings`
    /// as a leaf. `None` means there is no depth limit.
    pub depth_limit: Option<usize>,

    /// Selects what contributes to `depth_limit`.
    ///
    /// When true, only product nesting increments parse depth. When false,
    /// additions and powers also increment depth before their children are
    /// parsed.
    pub depth_is_product_depth: bool,

    /// Selects how shorthand notation is represented in the parsed network.
    ///
    /// `Expand` turns shorthand such as compact dot products, chains, traces,
    /// and Schoonschip vector notation into explicit graph structure. Fresh
    /// dummies created by this lowering are local to the expansion.
    ///
    /// `Opaque` keeps shorthand as a leaf tensor or scalar. Its inference mode
    /// controls whether the exposed structure comes from a fast syntactic walk
    /// or from expanding the shorthand and reading the resulting dangling slots.
    pub shorthand_parsing: ShorthandParsing,

    /// Allow parser implementations to treat composite scalar expressions as
    /// scalar-structured tensors.
    ///
    /// This flag is passed through to `Parse::parse_with_settings`. The
    /// symbolic tensor parser uses it for additive or multiplicative scalar
    /// composites so recursive contraction passes can inspect their internals
    /// instead of storing them as opaque pure scalars.
    pub parse_composite_scalars_as_tensors: bool,
}

impl Default for ParseSettings {
    fn default() -> Self {
        Self {
            precontract_scalars: true,
            take_first_term_from_sum: false,
            depth_limit: None,
            depth_is_product_depth: true,
            shorthand_parsing: ShorthandParsing::Expand,
            parse_composite_scalars_as_tensors: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ParseState<Aind = AbstractIndex> {
    depth: usize,
    next_dummy: Rc<Cell<usize>>,
    _aind: PhantomData<fn() -> Aind>,
}

#[allow(clippy::derivable_impls)]
impl<Aind> Default for ParseState<Aind> {
    fn default() -> Self {
        Self {
            depth: 0,
            next_dummy: Rc::new(Cell::new(1_000_000)),
            _aind: PhantomData,
        }
    }
}

impl<Aind: DummyAind> ParseState<Aind> {
    fn next(&self) -> Aind {
        let index = self.next_dummy.get();
        self.next_dummy.set(index + 1);
        Aind::new_dummy_at(index)
    }

    fn slot(&self, rep: &Representation<LibraryRep>) -> Slot<LibraryRep, Aind> {
        rep.slot(self.next())
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

struct SchoonschipMaterialization {
    replacement: Atom,
    factors: Vec<Atom>,
    compact_self: bool,
}

struct ShorthandMaterializer<'a, Aind> {
    state: &'a ParseState<Aind>,
    settings: &'a ParseSettings,
}

impl<'a, Aind: AbsInd + DummyAind + ParseableAind> ShorthandMaterializer<'a, Aind> {
    fn new(state: &'a ParseState<Aind>, settings: &'a ParseSettings) -> Self {
        Self { state, settings }
    }

    fn materialize_compact_vector_product(&self, value: MulView<'_>) -> Option<Atom> {
        let args = value.iter().collect::<Vec<_>>();
        let compact_positions = args
            .iter()
            .enumerate()
            .filter_map(|(position, arg)| Self::compact_vector_rep(*arg).map(|rep| (position, rep)))
            .collect::<Vec<_>>();
        let [(left_position, rep), (right_position, right_rep)] = compact_positions.as_slice()
        else {
            return None;
        };

        if rep != right_rep || !rep.rep.is_self_dual() {
            return None;
        }

        let slot = self.state.slot(rep).to_atom();
        let mut materialized = Atom::num(1);
        for (position, arg) in args.into_iter().enumerate() {
            let factor = if position == *left_position || position == *right_position {
                Self::materialize_compact_vector_with_slot(arg, rep, &slot)?
            } else {
                arg.to_owned()
            };
            materialized *= factor;
        }

        Some(materialized)
    }

    fn materialize_schoonschip(&self, value: AtomView<'_>) -> Option<Atom> {
        let materialized = self.materialize_schoonschip_arg(value)?;
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

    fn materialize_schoonschip_arg(
        &self,
        value: AtomView<'_>,
    ) -> Option<SchoonschipMaterialization> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        if self.settings.shorthand_parsing.expands()
            && let Some(materialized) = self.compact_scalar_product(fun)
        {
            return Some(materialized);
        }

        let mut changed = false;
        let mut factors = Vec::new();
        let mut rebuilt = FunctionBuilder::new(fun.get_symbol());

        for arg in fun.iter() {
            if let Some(materialized) = self.materialize_schoonschip_arg(arg) {
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
                replacement: Self::compact_dot_as_metric(replacement.as_view())
                    .unwrap_or(replacement),
                factors,
                compact_self: false,
            }
        })
    }

    fn compact_scalar_product(&self, value: FunView<'_>) -> Option<SchoonschipMaterialization> {
        if value.get_symbol() != ETS.metric && value.get_symbol() != SPENSO_TAG.dot {
            return None;
        }

        let args = value.iter().collect::<Vec<_>>();
        let [lhs, rhs] = args.as_slice() else {
            return None;
        };

        let rep = Self::compact_vector_rep(*lhs)?;
        if rep != Self::compact_vector_rep(*rhs)? || !rep.rep.is_self_dual() {
            return None;
        }

        let slot = self.state.slot(&rep).to_atom();
        let lhs = Self::materialize_compact_vector_with_slot(*lhs, &rep, &slot)?;
        let rhs = Self::materialize_compact_vector_with_slot(*rhs, &rep, &slot)?;

        Some(SchoonschipMaterialization {
            replacement: Atom::num(1),
            factors: vec![lhs, rhs],
            compact_self: true,
        })
    }

    fn compact_vector_rep(value: AtomView<'_>) -> Option<Representation<LibraryRep>> {
        match value {
            AtomView::Fun(fun) => Self::compact_tensor_rep_arg(fun).map(|(_, rep)| rep),
            AtomView::Add(add) => {
                let mut reps = add.iter().map(Self::compact_vector_rep);
                let rep = reps.next()??;
                reps.all(|candidate| candidate == Some(rep)).then_some(rep)
            }
            _ => None,
        }
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
            .any(|arg| Slot::<LibraryRep, Aind>::try_from(*arg).is_ok())
        {
            return None;
        }

        let rep_args = args
            .iter()
            .enumerate()
            .filter_map(|(position, arg)| {
                Self::compact_rep_pattern_match(*arg).map(|rep| (position, rep))
            })
            .collect::<Vec<_>>();

        let [(position, rep)] = rep_args.as_slice() else {
            return None;
        };

        Some((*position, *rep))
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

    fn materialize_compact_vector_with_slot(
        value: AtomView<'_>,
        rep: &Representation<LibraryRep>,
        slot: &Atom,
    ) -> Option<Atom> {
        match value {
            AtomView::Fun(fun) => {
                let (position, matched_rep) = Self::compact_tensor_rep_arg(fun)?;
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
                    .map(|term| Self::materialize_compact_vector_with_slot(term, rep, slot));
                let first = terms.next()??;
                let rest = terms.collect::<Option<Vec<_>>>()?;
                Some(rest.into_iter().fold(first, |sum, term| sum + term))
            }
            _ => None,
        }
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
            .all(|arg| Slot::<LibraryRep, Aind>::try_from(*arg).is_ok())
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
}

struct ChainExpansion;

impl ChainExpansion {
    fn replace_placeholders(value: AtomView<'_>, chain_in: &Atom, chain_out: &Atom) -> Atom {
        match value {
            AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in => chain_in.clone(),
            AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_out => chain_out.clone(),
            AtomView::Fun(fun) => {
                let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
                for arg in fun.iter() {
                    rebuilt = rebuilt.add_arg(Self::replace_placeholders(arg, chain_in, chain_out));
                }
                rebuilt.finish()
            }
            AtomView::Add(add) => add.iter().fold(Atom::Zero, |sum, term| {
                sum + Self::replace_placeholders(term, chain_in, chain_out)
            }),
            AtomView::Mul(mul) => mul.iter().fold(Atom::num(1), |product, factor| {
                product * Self::replace_placeholders(factor, chain_in, chain_out)
            }),
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                Self::replace_placeholders(base, chain_in, chain_out).pow(exponent.to_owned())
            }
            _ => value.to_owned(),
        }
    }

    fn scale_first_factor(value: AtomView<'_>, scalar: &Atom) -> Option<Atom> {
        let AtomView::Fun(fun) = value else {
            return None;
        };

        let offset = if fun.get_symbol() == SPENSO_TAG.chain {
            2
        } else if fun.get_symbol() == SPENSO_TAG.trace {
            1
        } else {
            return None;
        };

        if fun.get_nargs() <= offset {
            return None;
        }

        let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
        for (position, arg) in fun.iter().enumerate() {
            if position == offset {
                rebuilt = rebuilt.add_arg(scalar.clone() * arg.to_owned());
            } else {
                rebuilt = rebuilt.add_arg(arg);
            }
        }
        Some(rebuilt.finish())
    }
}

impl<
    'a,
    Sc,
    T: HasStructure + TensorStructure,
    K: Clone + Display + Debug,
    // FK: Clone + Display + Debug,
    Str: TensorScalarStore<Tensor = T, Scalar = Sc> + Clone,
    Aind: AbsInd + DummyAind + ParseableAind,
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
        S: TensorStructure + Clone + Parse + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let state = ParseState::<Aind>::default();
        Self::try_from_view_impl(value, state, library, settings)
    }

    #[allow(clippy::result_large_err)]
    fn try_from_view_impl<S, Lib: Library<S, Key = K>>(
        value: AtomView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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

    #[allow(clippy::type_complexity, clippy::result_large_err)]
    fn as_inferred_leaf<S>(
        value: AtomView<'a>,
        mode: StructureInferenceMode,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let s = S::structure_from_atom(value, mode);

        if let Ok(s) = s {
            Ok(Self::from_tensor(
                s.structure
                    .to_shell()
                    .concretize(Some(s.index_permutation.inverse())),
            ))
        } else {
            Ok(Self::from_scalar(value.try_into()?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_mul<S, Lib: Library<S, Key = K>>(
        value: MulView<'a>,
        mut state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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

        if settings.shorthand_parsing.expands()
            && let Some(materialized) = ShorthandMaterializer::<Aind>::new(&state, settings)
                .materialize_compact_vector_product(value)
        {
            return Self::try_from_view_impl(materialized.as_view(), state, library, settings);
        }

        state.depth += 1;
        // println!("{} for mul {}", state.depth, value.as_view());
        let mut iter = value.iter();
        let first_atom = iter.next().unwrap();
        let first = Self::try_from_view_impl(first_atom, state.clone(), library, settings)?;

        // state

        if settings.precontract_scalars {
            let mut scalars = Atom::num(1);

            let rest: Result<Vec<_>, _> = iter
                .filter_map(|a| {
                    match Self::try_from_view_impl(a, state.clone(), library, settings) {
                        Ok(n) => {
                            if let NetworkState::PureScalar = n.state {
                                scalars *= a;
                                None
                            } else {
                                Some(Ok((a.to_owned(), n)))
                            }
                        }
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect();

            let mut res = rest?;

            if let NetworkState::PureScalar = first.state {
                scalars *= first_atom;
            } else {
                res.push((first_atom.to_owned(), first));
            }

            if res.is_empty() {
                Ok(Self::from_scalar(value.as_view().try_into()?))
            } else if scalars != Atom::num(1)
                && res.len() == 1
                && let Some(scaled) =
                    ChainExpansion::scale_first_factor(res[0].0.as_view(), &scalars)
            {
                Self::try_from_view_impl(scaled.as_view(), state.clone(), library, settings)
            } else {
                let s = if scalars != Atom::num(1) {
                    Self::from_scalar(scalars.as_view().try_into()?)
                } else {
                    res.pop().unwrap().1
                };

                Ok(s.n_mul(res.into_iter().map(|(_, net)| net)))
            }
        } else {
            let rest: Result<Vec<_>, _> = iter
                .map(|a| Self::try_from_view_impl(a, state.clone(), library, settings))
                .collect();

            Ok(first.n_mul(rest?))
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_fun<S, Lib: Library<S, Key = K>>(
        value: FunView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
        // <PermutedStructure<S>>::Error: Debug,
    {
        let symbol = value.get_symbol();

        if let Some(inference) = settings.shorthand_parsing.opaque_inference_for(symbol) {
            return Self::as_inferred_leaf::<S>(value.as_view(), inference);
        }

        if symbol == SPENSO_TAG.bracket {
            let mut n_muls = value
                .iter()
                .map(|a| Self::try_from_view_impl(a, state.clone(), library, settings))
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
            } else if settings.shorthand_parsing.expands()
                && let Some(materialized) = ShorthandMaterializer::<Aind>::new(&state, settings)
                    .materialize_schoonschip(value.as_view())
            {
                Self::try_from_view_impl(materialized.as_view(), state, library, settings)
            } else {
                Ok(Self::from_scalar(value.as_view().try_into()?))
            }
        }
    }

    #[allow(clippy::result_large_err)]
    fn parse_chain<S, Lib: Library<S, Key = K>>(
        value: FunView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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

        let start = Slot::<LibraryRep, Aind>::try_from(args[0])
            .map_err(|err| eyre!("invalid chain start `{}`: {err}", args[0]))?;
        let end = Slot::<LibraryRep, Aind>::try_from(args[1])
            .map_err(|err| eyre!("invalid chain end `{}`: {err}", args[1]))?;
        let factors = &args[2..];

        if factors.is_empty() {
            let metric = FunctionBuilder::new(ETS.metric)
                .add_arg(start.to_atom())
                .add_arg(end.to_atom())
                .finish();
            return Self::try_from_view_impl(metric.as_view(), state, library, settings);
        }

        let mut factor_networks = Vec::new();
        let mut left = start;
        for (position, factor) in factors.iter().enumerate() {
            let right = if position + 1 == factors.len() {
                end.dual()
            } else {
                state.slot(&left.rep)
            };
            let factor = ChainExpansion::replace_placeholders(
                *factor,
                &left.to_atom(),
                &right.dual().to_atom(),
            );
            let factor = ShorthandMaterializer::<Aind>::new(&state, settings)
                .materialize_schoonschip(factor.as_view())
                .unwrap_or(factor);
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib>(
                factor,
                state.clone(),
                library,
                settings,
            )?);
            left = right;
        }

        Ok(Self::chain_like_network_product(factor_networks))
    }

    #[allow(clippy::result_large_err)]
    fn parse_trace<S, Lib: Library<S, Key = K>>(
        value: FunView<'a>,
        state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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

        if let [factor] = factors {
            let left = state.slot(&rep);
            let right = state.slot(&rep);
            let bridge = state.slot(&rep);
            let materialized_factor = ChainExpansion::replace_placeholders(
                *factor,
                &left.to_atom(),
                &right.dual().to_atom(),
            );
            let materialized_factor = ShorthandMaterializer::<Aind>::new(&state, settings)
                .materialize_schoonschip(materialized_factor.as_view())
                .unwrap_or(materialized_factor);
            let left_closure = FunctionBuilder::new(ETS.metric)
                .add_arg(bridge.to_atom())
                .add_arg(left.dual().to_atom())
                .finish();
            let right_closure = FunctionBuilder::new(ETS.metric)
                .add_arg(right.to_atom())
                .add_arg(bridge.dual().to_atom())
                .finish();
            let factor_net = Self::try_from_view_impl(
                materialized_factor.as_view(),
                state.clone(),
                library,
                settings,
            )?;
            let left_net =
                Self::try_from_view_impl(left_closure.as_view(), state.clone(), library, settings)?;
            let right_net = Self::try_from_view_impl(
                right_closure.as_view(),
                state.clone(),
                library,
                settings,
            )?;
            return Ok(factor_net.n_mul([right_net, left_net]));
        }

        let links = (0..factors.len())
            .map(|_| state.slot(&rep))
            .collect::<Vec<_>>();

        let materialized_factors = factors
            .iter()
            .enumerate()
            .map(|(position, factor)| {
                let left = links[position].to_atom();
                let right = links[(position + 1) % factors.len()].dual().to_atom();
                let factor = ChainExpansion::replace_placeholders(*factor, &left, &right);
                ShorthandMaterializer::<Aind>::new(&state, settings)
                    .materialize_schoonschip(factor.as_view())
                    .unwrap_or(factor)
            })
            .collect::<Vec<_>>();

        let mut factor_networks = Vec::new();
        for factor in materialized_factors {
            factor_networks.extend(Self::parse_chain_like_factor_networks::<S, Lib>(
                factor,
                state.clone(),
                library,
                settings,
            )?);
        }
        Ok(Self::chain_like_network_product(factor_networks))
    }

    #[allow(clippy::result_large_err)]
    fn parse_chain_like_factor_networks<S, Lib: Library<S, Key = K>>(
        factor: Atom,
        state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Vec<Self>, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
        TensorShell<S>: Concretize<T>,
        S::Slot: IsAbstractSlot<Aind = Aind>,
        T::Slot: IsAbstractSlot<Aind = Aind>,
    {
        let AtomView::Mul(product) = factor.as_view() else {
            return Ok(vec![Self::try_from_view_impl(
                factor.as_view(),
                state.clone(),
                library,
                settings,
            )?]);
        };

        let mut scalars = Vec::new();
        let mut tensors = Vec::new();
        for arg in product.iter() {
            let network = Self::try_from_view_impl(arg, state.clone(), library, settings)?;
            if network.state == NetworkState::PureScalar {
                scalars.push(network);
            } else {
                tensors.push(network);
            }
        }

        if scalars.is_empty() || tensors.len() != 1 {
            return Ok(vec![Self::try_from_view_impl(
                factor.as_view(),
                state.clone(),
                library,
                settings,
            )?]);
        }

        scalars.extend(tensors);
        Ok(scalars)
    }

    fn chain_like_network_product(networks: Vec<Self>) -> Self {
        let mut networks = networks.into_iter();
        let first = networks.next().unwrap();
        let rest = networks.collect::<Vec<_>>();
        if rest.is_empty() {
            first
        } else {
            first.n_mul(rest)
        }
    }

    #[allow(clippy::result_large_err)]
    fn try_from_pow<S, Lib: Library<S, Key = K>>(
        value: PowView<'a>,
        mut state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> std::result::Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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
        mut state: ParseState<Aind>,
        library: &Lib,
        settings: &ParseSettings,
    ) -> Result<Self, TensorNetworkError<K, Symbol>>
    where
        S: TensorStructure + Clone + Parse + StructureFromAtom,
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

        let first = Self::try_from_view_impl(first_atom, state.clone(), library, settings)?;
        if settings.take_first_term_from_sum {
            Ok(first)
        } else if settings.precontract_scalars {
            let mut scalars = Atom::Zero;

            let rest: Result<Vec<_>, _> = iter
                .filter_map(|a| {
                    match Self::try_from_view_impl(a, state.clone(), library, settings) {
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
                    }
                })
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
                    |a| match Self::try_from_view_impl(a, state.clone(), library, settings) {
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
    fn parse_to_atom_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>>;
}

impl NetworkParse for Atom {
    fn parse_to_atom_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        self.as_view().parse_to_atom_net::<Aind>(settings)
    }
}

impl NetworkParse for AtomView<'_> {
    fn parse_to_atom_net<Aind: AbsInd + DummyAind + ParseableAind>(
        &self,
        settings: &ParseSettings,
    ) -> Result<ParamNet<Aind>, TensorNetworkError<DummyKey, Symbol>> {
        let lib = DummyLibrary::<ParamTensor<ShadowedStructure<Aind>>>::new();

        ParamNet::<Aind>::try_from_view::<ShadowedStructure<Aind>, _>(*self, &lib, settings)
    }
}

#[cfg(test)]
mod test;
